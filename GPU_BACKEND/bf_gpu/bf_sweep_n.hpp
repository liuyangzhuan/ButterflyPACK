#pragma once

#ifdef H2_HAVE_GPU

// The butterfly matvec sweep (GPU_BACKEND/bf_gpu).  The index arithmetic here
// is taken from BF_block_MVP_dat_nonbatch and is the part that has to match
// exactly: a wrong mapping gives plausible numbers.
//
// The sweep has two halves that meet at level_half, and they are NOT the same
// loop.  Below the midpoint it runs over OUTPUT blocks and gathers two inputs
// into each; above it it runs over INPUT blocks and scatters each into two
// outputs.  At one process the two collapse into one formula, which is why a
// single-process version can be written without the distinction.  They do not
// collapse once the blocks are distributed: the halves order the blocks
// differently (row-wise below the midpoint, column-wise above, and the
// reverse when transposed), each rank owns a strided subset of a level, and
// which blocks a rank must exchange follows from that ordering.  The grids
// therefore come from the Fortran (bf_ranges.hpp), never from 2^l arithmetic
// here.
//
// The four loop forms, as the four branches of the CPU sweep write them.  All
// indices are global unless named _loc; `w` is 0 or 1 and selects one of the
// pair.  Note that the CPU takes the second of a pair as the next LOCAL block
// (index_jj_loc + 1), not as the local index of the next global block: the
// pair is always a consecutive pair, which is what the range doubling below
// the midpoint arranges.  Mirrored here, because with inc > 1 the two differ.
//
//   'N' below   out(i, j)            = sum_w K(i, 2j-1+w)   * in((i+1)/2, 2j-1+w)
//   'N' above   out(2i-1+w, (j+1)/2) +=      K(2i-1+w, j)   * in(i, j)
//   'T' below   out(i, j)            = sum_w K(2i-1+w, j)^T * in(2i-1+w, (j+1)/2)
//   'T' above   out((i+1)/2, 2j-1+w) +=      K(i, 2j-1+w)^T * in(i, j)
//
// Three communications sit between the levels, all mirrored on the device
// (bf_exchange.hpp):
//
//   * below the midpoint, after each level, a pair-wise BROADCAST that hands
//     a rank the second of the two input blocks its gather needs;
//   * at the midpoint, a redistribution between the row-wise and column-wise
//     orderings (BF_all2all_vec_n_ker);
//   * above the midpoint, after each level, a pair-wise REDUCTION that sums
//     the partial scatters of the ranks sharing an output block.
//
// The CPU allocates a vector block on first write and so uses beta = 0 there
// and beta = 1 afterwards.  Every block of a level exists here, so a level is
// zeroed once and every product accumulates with beta = 1.  That is also what
// the exchanges need: they add into the target, because a rank may already
// hold its own contribution to it.
//
// Within one batched GEMM no two items may write the same block.  Below the
// midpoint the pair writes one output, so the two halves of the pair go in
// separate passes; above it the pair writes two outputs but two different
// loop positions share one output, so the passes split on the parity of the
// index that halves.

#include "bf_exchange.hpp"
#include "bf_ranges.hpp"
#include "bf_sweep.hpp"

namespace bf {
namespace gpu {

namespace detail {

// The local index of a global block index, as the CPU computes it (one based
// there, zero based here).  The CPU divides without checking, because in a
// correct sweep the index always lands on a local block; one that does not is
// a disagreement between these grids and the Fortran's, so it is raised
// rather than silently folded onto the wrong block.
inline int local_of(int index, int idx, int inc, int n, const char* what) {
    const int d = index - idx;
    if (inc > 0 && d >= 0 && d % inc == 0) {
        const int k = d / inc;
        if (k < n) return k;
    }
    throw std::runtime_error(std::string("bf::gpu: ") + what + ": global index " + std::to_string(index) +
                             " is outside the local range (idx " + std::to_string(idx) + ", inc " +
                             std::to_string(inc) + ", n " + std::to_string(n) + ")");
}

// One product of the sweep: a kernel block, an input block and the output
// block it accumulates into, all as local indices.
struct Product {
    int ki = 0, kj = 0;
    int in_i = 0, in_j = 0;
    int out_i = 0, out_j = 0;
};

// The pair runs along the kernel's columns; `on_input` says whether the
// second of the pair is a second input ('N' below the midpoint: a gather) or
// a second output ('T' above it: a scatter).
inline Product pair_along_columns(const LevelDesc& k, const VecLevel& in, const VecLevel& out, int index_i,
                                  int index_j, int w, bool on_input) {
    Product p;
    p.kj = local_of(2 * index_j - 1, k.idx_c, k.inc_c, k.nc, "kernel column") + w;
    p.ki = local_of(index_i, k.idx_r, k.inc_r, k.nr, "kernel row");
    if (on_input) {  // 'N' below: out(i, j) gathers in((i+1)/2, 2j-1+w)
        p.in_i = local_of((index_i + 1) / 2, in.idx_r, in.inc_r, in.nr, "input row");
        p.in_j = local_of(2 * index_j - 1, in.idx_c, in.inc_c, in.nc, "input column") + w;
        p.out_i = local_of(index_i, out.idx_r, out.inc_r, out.nr, "output row");
        p.out_j = local_of(index_j, out.idx_c, out.inc_c, out.nc, "output column");
    } else {         // 'T' above: out((i+1)/2, 2j-1+w) takes in(i, j)
        p.in_i = local_of(index_i, in.idx_r, in.inc_r, in.nr, "input row");
        p.in_j = local_of(index_j, in.idx_c, in.inc_c, in.nc, "input column");
        p.out_i = local_of((index_i + 1) / 2, out.idx_r, out.inc_r, out.nr, "output row");
        p.out_j = local_of(2 * index_j - 1, out.idx_c, out.inc_c, out.nc, "output column") + w;
    }
    return p;
}

// The pair runs along the kernel's rows: the transpose of the above.
inline Product pair_along_rows(const LevelDesc& k, const VecLevel& in, const VecLevel& out, int index_i,
                               int index_j, int w, bool on_input) {
    Product p;
    p.ki = local_of(2 * index_i - 1, k.idx_r, k.inc_r, k.nr, "kernel row") + w;
    p.kj = local_of(index_j, k.idx_c, k.inc_c, k.nc, "kernel column");
    if (on_input) {  // 'T' below: out(i, j) gathers in(2i-1+w, (j+1)/2)
        p.in_i = local_of(2 * index_i - 1, in.idx_r, in.inc_r, in.nr, "input row") + w;
        p.in_j = local_of((index_j + 1) / 2, in.idx_c, in.inc_c, in.nc, "input column");
        p.out_i = local_of(index_i, out.idx_r, out.inc_r, out.nr, "output row");
        p.out_j = local_of(index_j, out.idx_c, out.inc_c, out.nc, "output column");
    } else {         // 'N' above: out(2i-1+w, (j+1)/2) takes in(i, j)
        p.in_i = local_of(index_i, in.idx_r, in.inc_r, in.nr, "input row");
        p.in_j = local_of(index_j, in.idx_c, in.inc_c, in.nc, "input column");
        p.out_i = local_of(2 * index_i - 1, out.idx_r, out.inc_r, out.nr, "output row") + w;
        p.out_j = local_of((index_j + 1) / 2, out.idx_c, out.inc_c, out.nc, "output column");
    }
    return p;
}

// Which loop form a sweep level takes, and which grid it loops over.
struct LevelForm {
    bool below = false;       // below the midpoint: a gather
    bool along_columns = false;  // the pair runs along the kernel's columns
    bool transposed = false;  // the kernel applies transposed
    int in_level = 0, out_level = 0;
};

inline LevelForm form_of(const BfDevice& bf, bool forward, int level) {
    const int L = bf.level_butterfly;
    LevelForm f;
    f.transposed = !forward;
    f.below = forward ? (level <= bf.level_half) : (level >= bf.level_half + 1);
    // 'N' below and 'T' above pair along the columns; the other two along the
    // rows.  Equivalently: the forward sweep pairs along columns exactly when
    // it is below the midpoint.
    f.along_columns = forward ? f.below : !f.below;
    f.in_level = forward ? level : L - level + 1;
    f.out_level = forward ? level + 1 : L - level + 2;
    return f;
}

// The grid a level's loop runs over: the kernel ranges in the ordering that
// half of the sweep uses.
inline const SweepRanges& loop_grid(const BfDevice& bf, bool forward, bool below) {
    const bool row_wise = forward ? below : !below;
    return row_wise ? bf.kernel_r : bf.kernel_c;
}

inline void set_grid(VecLevel& v, const SweepRanges& r, int l) {
    const size_t u = static_cast<size_t>(l);
    v.idx_r = r.idx_r[u]; v.inc_r = r.inc_r[u]; v.nr = r.nr[u];
    v.idx_c = r.idx_c[u]; v.inc_c = r.inc_c[u]; v.nc = r.nc[u];
    const size_t blocks = static_cast<size_t>(std::max(v.nr, 0)) * static_cast<size_t>(std::max(v.nc, 0));
    v.rows.assign(blocks, 0);
}

}  // namespace detail

// A sweep's plan for one orientation: the vector levels and the exchange
// plans of every level.  Built once per butterfly block and reused, because
// building an exchange plan costs two small collectives and a matvec is
// called tens of thousands of times in a STRUMPACK factorization.  Only the
// element offsets depend on nrhs, and those are a host-side prefix sum.
struct SweepPlan {
    char chara = 0;
    int nrhs = 0;
    std::vector<VecLevel> vec;
    // The midpoint level has two grids: the ordering it is written in by the
    // first half, and the one the redistribution leaves it in, which is what
    // the second half reads.  `vec[mid_level]` holds the second; `mid_src`
    // the first.
    VecLevel mid_src;
    int mid_level = -1;
    std::vector<ExchangePlan> pair;   // the exchange after sweep level l
    ExchangePlan mid;                 // the midpoint redistribution
    bool built = false;

    // The traffic the plan carries: records and rows sent and received, over
    // the pair-wise exchanges and the midpoint.  A sweep that matches the CPU
    // while moving nothing has only checked its index maps, so a validation
    // run reports this alongside the comparison.
    struct Traffic {
        int64_t pair_send = 0, pair_recv = 0, mid_send = 0, mid_recv = 0;
        int64_t pair_rows = 0, mid_rows = 0;
        int levels_exchanging = 0;
    };
    Traffic traffic() const {
        Traffic t;
        for (const ExchangePlan& p : pair) {
            if (p.send.empty() && p.recv.empty()) continue;
            t.levels_exchanging++;
            t.pair_send += static_cast<int64_t>(p.send.size());
            t.pair_recv += static_cast<int64_t>(p.recv.size());
            for (const ExchangeRecord& r : p.recv) t.pair_rows += r.rows;
        }
        t.mid_send = static_cast<int64_t>(mid.send.size());
        t.mid_recv = static_cast<int64_t>(mid.recv.size());
        for (const ExchangeRecord& r : mid.recv) t.mid_rows += r.rows;
        return t;
    }

    void release() {
        for (ExchangePlan& p : pair) release_plan(p);
        release_plan(mid);
        pair.clear();
        vec.clear();
        built = false;
        chara = 0;
        nrhs = 0;
        mid_level = -1;
    }
};

// Build a sweep's plan: the vector grids, every block's height, and the
// exchange plans.  The steps are ordered because they depend on each other -
// the midpoint's target heights are only known from the arrivals, and the
// second half's grids depend on the ordering the redistribution leaves.
//
// A block's height is the row count (or, transposed, the column count) of the
// kernel that writes it, and which kernel that is depends on the grids.
// Walking the sweep's own loops is what keeps the two in step; at one process
// every kernel in a row shares its row count, so a shortcut works there, but
// not with strided ownership.
//
// Collective over the block's communicator.  Called once per butterfly block
// and orientation; a later matvec with a different nrhs only re-packs.
inline void build_sweep_plan(const BfBridge& bridge, const BfDevice& bf, char chara, int nrhs,
                             SweepPlan& plan) {
    if (!bridge.usable()) {
        throw std::invalid_argument("bf::gpu::build_sweep_plan: the bridge is missing a rule");
    }
    const int L = bf.level_butterfly;
    const int lh = bf.level_half;
    const bool forward = (chara == 'N' || chara == 'n');
    plan.release();
    plan.chara = forward ? 'N' : 'T';
    plan.nrhs = nrhs;

    const SweepRanges& vr = forward ? bf.vec_n : bf.vec_t;
    plan.vec.assign(static_cast<size_t>(L) + 3, VecLevel());
    for (int l = 0; l <= L + 2; ++l) detail::set_grid(plan.vec[static_cast<size_t>(l)], vr, l);
    plan.mid_level = forward ? lh + 1 : L - lh + 1;
    plan.mid_src = plan.vec[static_cast<size_t>(plan.mid_level)];
    plan.pair.assign(static_cast<size_t>(L) + 2, ExchangePlan());

    // the last level of the first half writes the midpoint level in the
    // ordering it is produced in, not the one the second half reads
    const int boundary = forward ? lh : lh + 1;
    auto out_level_of = [&](const detail::LevelForm& f, int level) -> VecLevel& {
        if (f.below && level == boundary) return plan.mid_src;
        return plan.vec[static_cast<size_t>(f.out_level)];
    };

    // ---- the leaf product: V^T forward, U^T transposed, into vector level 1
    {
        const int leaf_level = forward ? 0 : L + 1;
        const detail::LevelForm f = detail::form_of(bf, forward, leaf_level);
        VecLevel& out = out_level_of(f, leaf_level);
        const LevelDesc& leaf = bf.level(leaf_level);
        for (int b = 0; b < leaf.nr * leaf.nc; ++b) {
            const BlockDesc& d = bf.blocks[static_cast<size_t>(leaf.first + b)];
            if (d.rows == 0) continue;
            const int index_i = forward ? 1 : b * leaf.inc_r + leaf.idx_r;
            const int index_j = forward ? b * leaf.inc_c + leaf.idx_c : 1;
            const int il = detail::local_of(index_i, out.idx_r, out.inc_r, out.nr, "leaf output row");
            const int jl = detail::local_of(index_j, out.idx_c, out.inc_c, out.nc, "leaf output column");
            out.rows[static_cast<size_t>(il) * out.nc + jl] = d.cols;
        }
    }

    // ---- the kernel levels of the first half
    auto size_level = [&](int level) {
        const detail::LevelForm f = detail::form_of(bf, forward, level);
        const LevelDesc& kerl = bf.level(level);
        const SweepRanges& g = detail::loop_grid(bf, forward, f.below);
        const size_t lu = static_cast<size_t>(level);
        const VecLevel& in = plan.vec[static_cast<size_t>(f.in_level)];
        VecLevel& out = out_level_of(f, level);
        for (int il = 0; il < g.nr[lu]; ++il) {
            for (int jl = 0; jl < g.nc[lu]; ++jl) {
                const int index_i = il * g.inc_r[lu] + g.idx_r[lu];
                const int index_j = jl * g.inc_c[lu] + g.idx_c[lu];
                for (int w = 0; w < 2; ++w) {
                    const detail::Product p =
                        f.along_columns
                            ? detail::pair_along_columns(kerl, in, out, index_i, index_j, w, f.below)
                            : detail::pair_along_rows(kerl, in, out, index_i, index_j, w, f.below);
                    const BlockDesc& k = bf.block(level, p.ki, p.kj);
                    out.rows[static_cast<size_t>(p.out_i) * out.nc + p.out_j] = f.transposed ? k.cols : k.rows;
                }
            }
        }
    };
    for (int step = 1; step <= L; ++step) {
        const int level = forward ? step : L - step + 1;
        if (!detail::form_of(bf, forward, level).below) continue;
        size_level(level);
    }
    for (VecLevel& v : plan.vec) v.pack(nrhs);
    plan.mid_src.pack(nrhs);

    // ---- the midpoint redistribution, which also sizes its target
    plan.mid = build_all2all_plan(bridge, boundary, forward ? 'C' : 'R', plan.mid_src,
                                  plan.vec[static_cast<size_t>(plan.mid_level)]);

    // ---- the kernel levels of the second half, now that the midpoint grid
    // is known
    for (int step = 1; step <= L; ++step) {
        const int level = forward ? step : L - step + 1;
        if (detail::form_of(bf, forward, level).below) continue;
        size_level(level);
    }
    for (int l = 0; l <= L + 2; ++l) {
        if (l != plan.mid_level) plan.vec[static_cast<size_t>(l)].pack(nrhs);
    }

    // ---- the closing product's sub-group.  When it spans more than one
    // rank the factor has a single local block and the vector block it
    // multiplies lives only on the head, which the CPU broadcasts before
    // applying it.  That broadcast is not mirrored yet, and without it the
    // non-head ranks would contribute zeros rather than fail, so the case is
    // refused instead of answered wrongly.
    if (bf.nproc > 1) {
        const int close_level = forward ? L + 1 : 0;
        const SweepRanges& g = forward ? bf.kernel_c : bf.kernel_r;
        const size_t cu = static_cast<size_t>(close_level);
        if (g.nr[cu] > 0 && g.nc[cu] > 0) {
            const BfBridge::Group grp = bridge.block_group(close_level, g.idx_r[cu], g.idx_c[cu],
                                                           forward ? 'C' : 'R');
            if (grp.nproc > 1) {
                throw std::runtime_error(
                    "bf::gpu::build_sweep_plan: the closing factor's sub-group spans " +
                    std::to_string(grp.nproc) +
                    " ranks, which needs the broadcast the CPU sweep does there and this does not yet");
            }
        }
    }

    // ---- the pair-wise exchanges, with the mode and collect the CPU passes
    if (bf.nproc > 1) {
        const char mode = forward ? 'R' : 'C';
        for (int step = 0; step <= L + 1; ++step) {
            const int level = forward ? step : L + 1 - step;
            const detail::LevelForm f = detail::form_of(bf, forward, level);
            // the CPU skips the exchange at the level that ends each half
            const bool ends_half = f.below ? (level == boundary)
                                           : (forward ? level == L + 1 : level == 0);
            if (ends_half) continue;
            VecLevel& v = out_level_of(f, level);
            plan.pair[static_cast<size_t>(level)] =
                build_plan(bridge, level, mode, f.below ? 'B' : 'R', v, v);
        }
    }
    plan.built = true;
}

// Re-pack a built plan for a different nrhs.  Only the element offsets and
// the message sizes change with it; the records, the peers and the grids do
// not, which is why a plan outlives a matvec.
inline void repack(SweepPlan& plan, int nrhs) {
    if (plan.nrhs == nrhs) return;
    for (VecLevel& v : plan.vec) v.pack(nrhs);
    plan.mid_src.pack(nrhs);
    plan.nrhs = nrhs;
}

// The sweep.  Vout = alpha * op(bf) * Vin + beta * Vout for 'N', and
// Vout = alpha * op(bf)^T * Vin + beta * Vout for 'T' and 'C'; every pointer
// device resident.  `plan` must have been built for this block and this
// orientation.
//
// The caller's alpha and beta belong to the closing product alone: the
// interior runs with alpha 1 and beta 1 over zeroed levels, and Vout is the
// only buffer the caller already holds a value in.
inline double sweep(const BfDevice& bf, SweepPlan& plan, int nrhs, const DataType* vin, int ldi,
                    DataType* vout, int ldo, DataType alpha, DataType beta, SweepScratch& s_scratch,
                    cudaStream_t stream = nullptr) {
    if (!plan.built || plan.vec.empty()) {
        throw std::invalid_argument("bf::gpu::sweep: the plan is not built");
    }
    const int L = bf.level_butterfly;
    const bool forward = (plan.chara == 'N');
    repack(plan, nrhs);

    // the input and output row partitions: the accumulated heights of this
    // rank's V and U blocks (arr_acc_n and arr_acc_m on the host)
    const LevelDesc& dV = bf.level(0);
    const LevelDesc& dU = bf.level(L + 1);
    std::vector<int64_t> acc_n(static_cast<size_t>(dV.nc) + 1, 0);
    for (int j = 0; j < dV.nc; ++j)
        acc_n[static_cast<size_t>(j) + 1] = acc_n[static_cast<size_t>(j)] + bf.block(0, 0, j).rows;
    std::vector<int64_t> acc_m(static_cast<size_t>(dU.nr) + 1, 0);
    for (int i = 0; i < dU.nr; ++i)
        acc_m[static_cast<size_t>(i) + 1] = acc_m[static_cast<size_t>(i)] + bf.block(L + 1, i, 0).rows;

    // one arena block per vector level; every level is zeroed, because the
    // products accumulate with beta 1 and so do the exchanges
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    std::vector<DataType*> owned;
    auto attach = [&](VecLevel& v) {
        if (v.elements == 0) { v.data = nullptr; return; }
        v.data = heap.alloc<DataType>(static_cast<size_t>(v.elements) * sizeof(DataType),
                                      fmm::gpu::ResidentKind::Matvec);
        owned.push_back(v.data);
        fmm::gpu::check_cuda(cudaMemsetAsync(v.data, 0, static_cast<size_t>(v.elements) * sizeof(DataType),
                                             stream), "bf sweep zero");
    };
    for (VecLevel& v : plan.vec) attach(v);
    attach(plan.mid_src);

    double flops = 0.0;
    std::vector<GemmItem> items;
    const int boundary = forward ? bf.level_half : bf.level_half + 1;
    auto out_level_of = [&](const detail::LevelForm& f, int level) -> VecLevel& {
        if (f.below && level == boundary) return plan.mid_src;
        return plan.vec[static_cast<size_t>(f.out_level)];
    };

    // ---- the leaf product: the factor applied transposed to the input
    auto leaf_product = [&](int leaf_level) {
        const detail::LevelForm f = detail::form_of(bf, forward, leaf_level);
        VecLevel& out = out_level_of(f, leaf_level);
        const LevelDesc& leaf = bf.level(leaf_level);
        const std::vector<int64_t>& acc = forward ? acc_n : acc_m;
        items.clear();
        for (int b = 0; b < leaf.nr * leaf.nc; ++b) {
            const BlockDesc& d = bf.blocks[static_cast<size_t>(leaf.first + b)];
            if (d.rows == 0 || d.cols == 0) continue;
            const int index_i = forward ? 1 : b * leaf.inc_r + leaf.idx_r;
            const int index_j = forward ? b * leaf.inc_c + leaf.idx_c : 1;
            const int il = detail::local_of(index_i, out.idx_r, out.inc_r, out.nr, "leaf output row");
            const int jl = detail::local_of(index_j, out.idx_c, out.inc_c, out.nc, "leaf output column");
            const size_t ob = static_cast<size_t>(il) * out.nc + jl;
            GemmItem it;
            it.a = bf.values + d.offset;  it.lda = d.lda;
            it.b = vin + acc[static_cast<size_t>(b)];  it.ldb = ldi;
            it.c = out.data + out.offset[ob];  it.ldc = out.rows[ob];
            it.m = d.cols; it.n = nrhs; it.k = d.rows;   // A^T: (cols x rows)
            items.push_back(it);
        }
        flops += issue(items, DataType(1), DataType(1), /*trans_a=*/true, s_scratch);
    };

    // ---- one kernel level
    auto run_level = [&](int level) {
        const detail::LevelForm f = detail::form_of(bf, forward, level);
        const LevelDesc& kerl = bf.level(level);
        const SweepRanges& g = detail::loop_grid(bf, forward, f.below);
        const size_t lu = static_cast<size_t>(level);
        const VecLevel& in = plan.vec[static_cast<size_t>(f.in_level)];
        VecLevel& out = out_level_of(f, level);
        if (out.elements == 0) return;
        // Below the midpoint the pair writes one output, so its two halves go
        // in separate batches; above it the pair writes two outputs but two
        // loop positions share one, so the batches split on the parity of the
        // index that halves.
        for (int pass = 0; pass < 2; ++pass) {
            items.clear();
            for (int il = 0; il < g.nr[lu]; ++il) {
                for (int jl = 0; jl < g.nc[lu]; ++jl) {
                    const int index_i = il * g.inc_r[lu] + g.idx_r[lu];
                    const int index_j = jl * g.inc_c[lu] + g.idx_c[lu];
                    if (!f.below) {
                        const int split = f.along_columns ? index_i : index_j;
                        if (split % 2 != 1 - pass) continue;
                    }
                    for (int w = 0; w < 2; ++w) {
                        if (f.below && w != pass) continue;
                        const detail::Product p =
                            f.along_columns
                                ? detail::pair_along_columns(kerl, in, out, index_i, index_j, w, f.below)
                                : detail::pair_along_rows(kerl, in, out, index_i, index_j, w, f.below);
                        const BlockDesc& k = bf.block(level, p.ki, p.kj);
                        if (k.rows == 0 || k.cols == 0) continue;
                        const size_t ib = static_cast<size_t>(p.in_i) * in.nc + p.in_j;
                        const size_t ob = static_cast<size_t>(p.out_i) * out.nc + p.out_j;
                        if (in.rows[ib] == 0 || out.rows[ob] == 0) continue;
                        GemmItem it;
                        it.a = bf.values + k.offset;  it.lda = k.lda;
                        it.b = in.data + in.offset[ib];  it.ldb = in.rows[ib];
                        it.c = out.data + out.offset[ob];  it.ldc = out.rows[ob];
                        it.m = f.transposed ? k.cols : k.rows;
                        it.n = nrhs;
                        it.k = f.transposed ? k.rows : k.cols;
                        items.push_back(it);
                    }
                }
            }
            flops += issue(items, DataType(1), DataType(1), f.transposed, s_scratch);
        }
    };

    // ---- the closing product, with the caller's scaling
    auto closing_product = [&](int close_level) {
        const LevelDesc& d = bf.level(close_level);
        const VecLevel& in = plan.vec[static_cast<size_t>(L) + 1];
        const std::vector<int64_t>& acc = forward ? acc_m : acc_n;
        items.clear();
        for (int b = 0; b < d.nr * d.nc; ++b) {
            const BlockDesc& k = bf.blocks[static_cast<size_t>(d.first + b)];
            if (k.rows == 0 || k.cols == 0) continue;
            const int index_i = forward ? b * d.inc_r + d.idx_r : 1;
            const int index_j = forward ? 1 : b * d.inc_c + d.idx_c;
            const int il = detail::local_of(index_i, in.idx_r, in.inc_r, in.nr, "closing input row");
            const int jl = detail::local_of(index_j, in.idx_c, in.inc_c, in.nc, "closing input column");
            const size_t ib = static_cast<size_t>(il) * in.nc + jl;
            if (in.rows[ib] == 0) continue;
            GemmItem it;
            it.a = bf.values + k.offset;  it.lda = k.lda;
            it.b = in.data + in.offset[ib];  it.ldb = in.rows[ib];
            it.c = vout + acc[static_cast<size_t>(b)];  it.ldc = ldo;
            it.m = k.rows; it.n = nrhs; it.k = k.cols;
            items.push_back(it);
        }
        flops += issue(items, alpha, beta, false, s_scratch);
    };

    auto exchange_after = [&](int level, int tag) {
        ExchangePlan& ex = plan.pair[static_cast<size_t>(level)];
        if (ex.empty()) return;
        const detail::LevelForm f = detail::form_of(bf, forward, level);
        VecLevel& v = out_level_of(f, level);
        exchange_level(ex, v, v, nrhs, tag, stream);
    };

    // ---- the first half: the leaf product, then its kernel levels, each
    // followed by the pair-wise broadcast.  The level that ends the half has
    // no exchange - the redistribution follows it instead.
    const int step = forward ? 1 : -1;
    const int first = forward ? 0 : L + 1;
    for (int level = first;; level += step) {
        if (level == first) leaf_product(level);
        else run_level(level);
        if (level == boundary) break;
        exchange_after(level, 100 + level);
    }

    // ---- the midpoint: between the two orderings
    if (!plan.mid.empty()) {
        exchange_level(plan.mid, plan.mid_src, plan.vec[static_cast<size_t>(plan.mid_level)], nrhs, 200, stream);
    } else if (plan.mid_src.data != nullptr && plan.mid_level >= 0) {
        // one rank: the orderings coincide, so the level is already in place
        VecLevel& dst = plan.vec[static_cast<size_t>(plan.mid_level)];
        if (dst.data != nullptr && dst.elements == plan.mid_src.elements) {
            fmm::gpu::check_cuda(cudaMemcpyAsync(dst.data, plan.mid_src.data,
                                                 static_cast<size_t>(dst.elements) * sizeof(DataType),
                                                 cudaMemcpyDeviceToDevice, stream), "bf sweep midpoint");
        }
    }

    // ---- the second half: its kernel levels, each followed by the pair-wise
    // reduction, then the closing product
    const int last = forward ? L + 1 : 0;
    for (int level = boundary + step;; level += step) {
        if (level == last) closing_product(level);
        else run_level(level);
        if (level == last) break;
        exchange_after(level, 300 + level);
    }

    for (DataType* p : owned) heap.free(p);
    for (VecLevel& v : plan.vec) v.data = nullptr;
    plan.mid_src.data = nullptr;
    return flops;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
