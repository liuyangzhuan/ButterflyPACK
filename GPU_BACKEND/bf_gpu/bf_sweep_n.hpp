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
// and beta = 1 afterwards.  Here, on one process, the first batch writing a
// level writes all of it, so it takes beta = 0 and the next beta = 1; with
// several, a level is zeroed before its first write and every product
// accumulates with beta = 1, because the exchanges add into the target - a
// rank may already hold its own contribution to it - and some blocks only
// ever receive.
//
// Within one batched GEMM no two items may write the same block.  Below the
// midpoint the pair writes one output, so the two halves of the pair go in
// separate passes; above it the pair writes two outputs but two different
// loop positions share one output, so the passes split on the parity of the
// index that halves.

#include <cstring>
#include <map>
#include <memory>

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
inline Product pair_along_columns(const LevelDesc& k, const VecGrid& in, const VecGrid& out, int index_i,
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
inline Product pair_along_rows(const LevelDesc& k, const VecGrid& in, const VecGrid& out, int index_i,
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

template<typename T>
inline LevelForm form_of(const BfDevice<T>& bf, bool forward, int level) {
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
template<typename T>
inline const SweepRanges& loop_grid(const BfDevice<T>& bf, bool forward, bool below) {
    const bool row_wise = forward ? below : !below;
    return row_wise ? bf.kernel_r : bf.kernel_c;
}

inline void set_grid(VecGrid& v, const SweepRanges& r, int l) {
    const size_t u = static_cast<size_t>(l);
    v.idx_r = r.idx_r[u]; v.inc_r = r.inc_r[u]; v.nr = r.nr[u];
    v.idx_c = r.idx_c[u]; v.inc_c = r.inc_c[u]; v.nc = r.nc[u];
    const size_t blocks = static_cast<size_t>(std::max(v.nr, 0)) * static_cast<size_t>(std::max(v.nc, 0));
    v.rows.assign(blocks, 0);
}

}  // namespace detail

// ---------------------------------------------------------------------------
// A sweep recorded once per plan and nrhs, replayed by every call.
//
// Measured on one A100 (262144^2, L = 12, 2026-10-25): rebuilding every
// level's GEMM list on the host and uploading its nine arrays per call cost
// 0.6 ms of copies and 1.6 ms of host time at one column, and zeroing every
// vector level 3.5 ms at 256 columns.  Recorded, a call runs one kernel that
// turns the vector offsets into pointers, then the GEMMs and the exchanges;
// the vector levels take turns in two buffers instead of one each.
// ---------------------------------------------------------------------------

// One GEMM: its factor block by offset into the factors, its input and output
// by buffer and element offset.  Buffers: 0 the caller's x, 1 and 2 the two
// the vector levels take turns in, 3 the caller's y.
struct SweepItem {
    const void* a = nullptr;  // the factor block itself (a T*), when it lives outside the factors' block
    int64_t a_off = 0, b_off = 0, c_off = 0;
    int b_buf = 0, c_buf = 0;
    int m = 0, n = 0, k = 0, lda = 0, ldb = 0, ldc = 0;
};

// One batched GEMM: items [first, first + count).  beta is 0 where the batch
// is the first to write a level and writes all of it (one process), else 1
// over a level zeroed first; the closing product takes the caller's alpha
// and beta.
struct SweepBatch {
    int first = 0, count = 0;
    bool trans_a = false;
    bool closing = false;
    double beta = 1.0;  // 0 or 1
    int max_m = 0, max_n = 0, max_k = 0;
    int tc_first = 0, tc32 = 0, tc48 = 0, tc64 = 0;  // its tensor-core tiles in the exec's list (32s, 48s, 64s)
    double flops = 0.0;
};

// The recorded sweep: the steps in order, the batches, which buffer each
// vector level lives in, and the GEMMs' arrays on the device.
// A step of a recorded sweep.
enum class SweepOp { Launch, Zero, Exchange, Midpoint, Copy, Group };
struct SweepStep {
    SweepOp op = SweepOp::Launch;
    int arg = 0;     // Launch: the batch; Zero, Copy: the level written (index into levels);
                     // Exchange: the sweep level whose pair exchange this is; Group: the
                     // plan's group
    int level = -1;  // Exchange, Group: the level exchanged; Copy: the level copied from
    int tag = 0;
};

template<typename T>
struct SweepExec {
    using DataType = T;
    int nrhs = 0;
    std::vector<SweepStep> program;
    std::vector<SweepBatch> batches;
    std::vector<VecLevel<T>*> levels;  // the vector levels, in the order the sweep writes them
    std::vector<int> buffer;         // 1 or 2, per entry of levels
    int64_t buffer_elements = 0;
    int count = 0;
    double flops = 0.0;
    char* block = nullptr;           // the device arrays, one block of the arena
    const DataType** d_a = nullptr;
    const DataType** d_b = nullptr;
    DataType** d_c = nullptr;
    magma_int_t *d_m = nullptr, *d_n = nullptr, *d_k = nullptr;
    magma_int_t *d_lda = nullptr, *d_ldb = nullptr, *d_ldc = nullptr;
    int64_t *d_boff = nullptr, *d_coff = nullptr;
    signed char *d_bbuf = nullptr, *d_cbuf = nullptr;
    int2* d_tc = nullptr;

    SweepExec() = default;
    SweepExec(const SweepExec&) = delete;
    SweepExec& operator=(const SweepExec&) = delete;
    ~SweepExec() { release(); }
    void release() {
        if (block != nullptr) fmm::gpu::DeviceHeap::instance().free(block);
        block = nullptr;
    }
};

// Batches of few columns (bf_sweep_exec.cu): every element of C summed in
// MAGMA's order, so the result is MAGMA's to the bit.  sweep_narrow_kind:
// 0 none fits, else the kernel launch_sweep_narrow will use.
int sweep_narrow_kind(bool trans_a, int max_m, int max_n);
void launch_sweep_narrow(bool trans_a, int count, const int* m, const int* n, const int* k, double alpha,
                         const double* const* a, const int* lda, const double* const* b, const int* ldb,
                         double beta, double* const* c, const int* ldc, int max_m, int max_n, cudaStream_t stream);
// The narrow kernel of a batch of scalar T: they are real kernels, so a
// complex batch has none.
template<typename T>
inline int narrow_kind(bool trans_a, int max_m, int max_n) {
    return std::is_same<T, double>::value ? sweep_narrow_kind(trans_a, max_m, max_n) : 0;
}

// The kernel of a replay (bf_sweep_exec.cu): every GEMM's input and output
// pointer from its offset and buffer, and the leading dimensions of x and y.
// dst = conj(src), rows x cols, column major, in place when dst is src: the
// conjugations around a 'T' sweep that make it 'C' on complex data
// (op^H x = conj(op^T conj(x))).  (These two defined for double and
// fmm::gpu::dcomplex.)
template<typename T>
void launch_conj(T* dst, int64_t ldd, const T* src, int64_t lds, int64_t rows, int64_t cols, cudaStream_t stream);
template<typename T>
void launch_sweep_patch(int count, const int64_t* b_off, const int64_t* c_off, const signed char* b_buf,
                        const signed char* c_buf, const T* vin, T* buf0, T* buf1, T* vout, int ldi, int ldo,
                        const T** b, T** c, int* ldb, int* ldc, cudaStream_t stream);

// A sweep's plan for one orientation: the vector levels and the exchange
// plans of every level.  Built once per butterfly block and reused, because
// building an exchange plan costs two small collectives and a matvec is
// called tens of thousands of times in a STRUMPACK factorization.  Only the
// element offsets depend on nrhs, and those are a host-side prefix sum.
template<typename T>
struct SweepPlan {
    using DataType = T;
    char chara = 0;
    int nrhs = 0;
    std::vector<VecLevel<T>> vec;
    // The midpoint level has two grids: the ordering it is written in by the
    // first half, and the one the redistribution leaves it in, which is what
    // the second half reads.  `vec[mid_level]` holds the second; `mid_src`
    // the first.
    VecLevel<T> mid_src;
    int mid_level = -1;
    std::vector<ExchangePlan> pair;   // the exchange after sweep level l
    ExchangePlan mid;                 // the midpoint redistribution
    // a leaf several ranks share: the sum after the opening factor and the
    // broadcast before the closing one (bf_exchange.hpp), -1 for none
    std::vector<GroupPlan> groups;
    int open_group = -1, close_group = -1;
    // a partial plan (build_sweep_plan's stop): its last level, and the
    // vector level it leaves its result in
    int stop = -1, result = -1;
    bool built = false;
    // the sweep recorded for each nrhs a caller has used (sweep_exec)
    std::map<int, std::unique_ptr<SweepExec<T>>> execs;

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
        execs.clear();
        for (ExchangePlan& p : pair) release_plan(p);
        release_plan(mid);
        pair.clear();
        vec.clear();
        groups.clear();
        open_group = close_group = -1;
        stop = result = -1;
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
// `stop`, a level of the first half, asks for BF_block_MVP_partial's sweep
// instead: the leaf product and the first half's levels up to `stop`, each
// with its pair-wise broadcast, or with the midpoint redistribution when
// `stop` ends the half.  The construction from matvecs sweeps its half-built
// butterfly so, and the levels past `stop` need not exist yet.  The result is
// plan.vec[plan.result].
//
// Collective over the block's communicator.  Called once per butterfly block
// and orientation; a later matvec with a different nrhs only re-packs.
template<typename T>
inline void build_sweep_plan(const BfBridge& bridge, const BfDevice<T>& bf, char chara, int nrhs,
                             SweepPlan<T>& plan, int stop = -1) {
    using VecLevel = bf::gpu::VecLevel<T>;
    if (!bridge.usable()) {
        throw std::invalid_argument("bf::gpu::build_sweep_plan: the bridge is missing a rule");
    }
    const int L = bf.level_butterfly;
    const int lh = bf.level_half;
    const bool forward = (chara == 'N' || chara == 'n');
    plan.release();
    plan.chara = forward ? 'N' : 'T';
    plan.nrhs = nrhs;
    plan.stop = stop;

    const SweepRanges& vr = forward ? bf.vec_n : bf.vec_t;
    plan.vec.assign(static_cast<size_t>(L) + 3, VecLevel());
    for (int l = 0; l <= L + 2; ++l) detail::set_grid(plan.vec[static_cast<size_t>(l)], vr, l);
    plan.mid_level = forward ? lh + 1 : L - lh + 1;
    plan.mid_src = plan.vec[static_cast<size_t>(plan.mid_level)];
    plan.pair.assign(static_cast<size_t>(L) + 2, ExchangePlan());

    // the last level of the first half writes the midpoint level in the
    // ordering it is produced in, not the one the second half reads
    const int boundary = forward ? lh : lh + 1;
    if (stop >= 0 && (forward ? stop > boundary : (stop < boundary || stop > L + 1))) {
        throw std::invalid_argument("bf::gpu::build_sweep_plan: stop " + std::to_string(stop) +
                                    " is not in the first half");
    }
    const bool full = stop < 0;
    // the first half's levels the plan covers
    auto covered = [&](int level) { return full || (forward ? level <= stop : level >= stop); };
    auto out_level_of = [&](const detail::LevelForm& f, int level) -> VecLevel& {
        if (f.below && level == boundary) return plan.mid_src;
        return plan.vec[static_cast<size_t>(f.out_level)];
    };
    auto vec_index = [&](const VecLevel& v) {
        return &v == &plan.mid_src ? -1 : static_cast<int>(&v - plan.vec.data());
    };

    // ---- the leaf product: V^T forward, U^T transposed, into vector level 1.
    // A block of no local rows still has its height: on a leaf several ranks
    // share, the rank's zeros go into the group's sum.
    const int leaf_level = forward ? 0 : L + 1;
    {
        const detail::LevelForm f = detail::form_of(bf, forward, leaf_level);
        VecLevel& out = out_level_of(f, leaf_level);
        const LevelDesc& leaf = bf.level(leaf_level);
        for (int b = 0; b < leaf.nr * leaf.nc; ++b) {
            const BlockDesc& d = bf.blocks[static_cast<size_t>(leaf.first + b)];
            if (d.cols == 0) continue;
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
        if (!detail::form_of(bf, forward, level).below || !covered(level)) continue;
        size_level(level);
    }
    for (VecLevel& v : plan.vec) v.pack(nrhs);
    plan.mid_src.pack(nrhs);

    // ---- a leaf several ranks share: the sum of the opening product over
    // the leaf's group.  Such a rank holds the one leaf.
    if (bf.nproc > 1) {
        const LevelDesc& leaf = bf.level(leaf_level);
        if (leaf.nr * leaf.nc == 1) {
            const int index_i = forward ? 1 : leaf.idx_r, index_j = forward ? leaf.idx_c : 1;
            const BfBridge::Group grp = bridge.block_group(leaf_level, index_i, index_j, forward ? 'R' : 'C');
            if (grp.nproc > 1) {
                const VecLevel& out = out_level_of(detail::form_of(bf, forward, leaf_level), leaf_level);
                const size_t b = static_cast<size_t>(out.row_local(index_i)) * out.nc + out.col_local(index_j);
                plan.open_group = static_cast<int>(plan.groups.size());
                plan.groups.push_back(build_group_plan(grp.comm, true, vec_index(out), b));
            }
        }
    }

    // ---- the midpoint redistribution, which also sizes its target
    if (full || stop == boundary) {
        plan.mid = build_all2all_plan(bridge, boundary, forward ? 'C' : 'R', plan.mid_src,
                                      plan.vec[static_cast<size_t>(plan.mid_level)]);
    }

    if (full) {
        // ---- the kernel levels of the second half, now that the midpoint
        // grid is known
        for (int step = 1; step <= L; ++step) {
            const int level = forward ? step : L - step + 1;
            if (detail::form_of(bf, forward, level).below) continue;
            size_level(level);
        }
        for (int l = 0; l <= L + 2; ++l) {
            if (l != plan.mid_level) plan.vec[static_cast<size_t>(l)].pack(nrhs);
        }

        // ---- the closing factor's group.  When it spans more than one rank
        // the factor has a single local block, and the vector block it
        // multiplies is on the group's head alone, which broadcasts it.  The
        // other ranks size the block by the factor, since no product of
        // theirs writes it.
        if (bf.nproc > 1) {
            const int close_level = forward ? L + 1 : 0;
            const LevelDesc& d = bf.level(close_level);
            if (d.nr * d.nc == 1) {
                const int index_i = forward ? d.idx_r : 1, index_j = forward ? 1 : d.idx_c;
                const BfBridge::Group grp = bridge.block_group(close_level, index_i, index_j, forward ? 'C' : 'R');
                if (grp.nproc > 1) {
                    VecLevel& in = plan.vec[static_cast<size_t>(L) + 1];
                    const size_t b = static_cast<size_t>(in.row_local(index_i)) * in.nc + in.col_local(index_j);
                    const int h = bf.blocks[static_cast<size_t>(d.first)].cols;
                    if (in.rows[b] == 0) {
                        in.rows[b] = h;
                        in.pack(nrhs);
                        if (L + 1 == plan.mid_level) {
                            // the redistribution delivered nothing into this
                            // block, but may have into others of its level
                            for (ExchangeRecord& r : plan.mid.recv) {
                                r.row_offset = in.row_offset[static_cast<size_t>(in.row_local(r.index_i)) * in.nc +
                                                             in.col_local(r.index_j)];
                            }
                            if (!plan.mid.recv.empty()) {
                                fmm::gpu::check_cuda(cudaMemcpy(plan.mid.d_recv, plan.mid.recv.data(),
                                                                plan.mid.recv.size() * sizeof(ExchangeRecord),
                                                                cudaMemcpyHostToDevice),
                                                     "bf sweep plan: the midpoint's records");
                            }
                        }
                    } else if (in.rows[b] != h) {
                        throw std::runtime_error("bf::gpu::build_sweep_plan: the closing block has " +
                                                 std::to_string(in.rows[b]) + " rows against the factor's " +
                                                 std::to_string(h));
                    }
                    plan.close_group = static_cast<int>(plan.groups.size());
                    plan.groups.push_back(build_group_plan(grp.comm, false, L + 1, b));
                }
            }
        }
    }

    // ---- the pair-wise exchanges, with the mode and collect the CPU passes
    if (bf.nproc > 1) {
        const char mode = forward ? 'R' : 'C';
        for (int step = 0; step <= L + 1; ++step) {
            const int level = forward ? step : L + 1 - step;
            const detail::LevelForm f = detail::form_of(bf, forward, level);
            if (!(f.below ? covered(level) : full)) continue;
            // the CPU skips the exchange at the level that ends each half
            const bool ends_half = f.below ? (level == boundary)
                                           : (forward ? level == L + 1 : level == 0);
            if (ends_half) continue;
            VecLevel& v = out_level_of(f, level);
            plan.pair[static_cast<size_t>(level)] =
                build_plan(bridge, level, mode, f.below ? 'B' : 'R', v, v);
        }
    }
    if (!full) {
        plan.result = stop == boundary ? plan.mid_level
                                        : detail::form_of(bf, forward, stop).out_level;
    }
    plan.built = true;
}

// Re-pack a built plan for a different nrhs.  Only the element offsets and
// the message sizes change with it; the records, the peers and the grids do
// not, which is why a plan outlives a matvec.
template<typename T>
inline void repack(SweepPlan<T>& plan, int nrhs) {
    if (plan.nrhs == nrhs) return;
    for (VecLevel<T>& v : plan.vec) v.pack(nrhs);
    plan.mid_src.pack(nrhs);
    plan.nrhs = nrhs;
}

// The leaf and closing products' row partitions: where each of this rank's
// V and U blocks starts in x and y (arr_acc_n and arr_acc_m on the host).
template<typename T>
inline void leaf_offsets(const BfDevice<T>& bf, std::vector<int64_t>& acc_n, std::vector<int64_t>& acc_m) {
    const int L = bf.level_butterfly;
    const LevelDesc& dV = bf.level(0);
    const LevelDesc& dU = bf.level(L + 1);
    acc_n.assign(static_cast<size_t>(dV.nc) + 1, 0);
    for (int j = 0; j < dV.nc; ++j)
        acc_n[static_cast<size_t>(j) + 1] = acc_n[static_cast<size_t>(j)] + bf.block(0, 0, j).rows;
    acc_m.assign(static_cast<size_t>(dU.nr) + 1, 0);
    for (int i = 0; i < dU.nr; ++i)
        acc_m[static_cast<size_t>(i) + 1] = acc_m[static_cast<size_t>(i)] + bf.block(L + 1, i, 0).rows;
}

// Whether a batch runs on the FP64 tensor cores (fmm::gpu::launch_dgemm_vbatched_tc,
// which rounds as MAGMA does) when no narrow kernel takes it: transposed
// always; untransposed from 25 columns, where MAGMA's GEMMs slow by half
// (8192 GEMMs of 50 x 50 on one A100: 0.27 ms at 24 columns, 0.47 at 28;
// the tensor cores 0.30 and 0.32), and for blocks of more than 64 rows and
// 32 inner columns at any width (16 columns of 128 x 64: 0.83 against 0.47;
// of 128 x 24, MAGMA's 0.21 against 0.34).  Measured on the same A100
// (262144^2, L = 12, ms a matvec, with 48 and 64 tiles): untransposed 9.2
// against 6.6 at 64 columns and 33 against 24 at 256.
template<typename T>
inline bool sweep_tensor_cores(bool trans_a, int max_m, int max_n, int max_k) {
    if (!std::is_same<T, double>::value) return false;
    return trans_a || max_n >= 25 || (max_m > 64 && max_k > 32);
}

// Up to this many columns a batch's tensor-core tiles are 32 x 32, which
// leave the least of a tile empty there (8192 GEMMs of 50 x 50 at 32
// columns: 0.34 ms against 0.42 with one 64 tile each); wider, each GEMM
// takes 48s or 64s by what it leaves empty.
constexpr int kTcNarrowCols = 32;

// Record the sweep of `plan` for nrhs: the same walk as the CPU sweep's
// levels (the index arithmetic above), but every GEMM kept as offsets, and
// the order of the GEMMs, zeroings and exchanges kept as a program.
template<typename T>
inline SweepExec<T>& sweep_exec(const BfDevice<T>& bf, SweepPlan<T>& plan, int nrhs) {
    using DataType = T;
    using VecLevel = bf::gpu::VecLevel<T>;
    using SweepExec = bf::gpu::SweepExec<T>;
    auto found = plan.execs.find(nrhs);
    if (found != plan.execs.end()) return *found->second;
    repack(plan, nrhs);
    auto owned = std::make_unique<SweepExec>();
    SweepExec& ex = *owned;
    ex.nrhs = nrhs;
    const int L = bf.level_butterfly;
    const bool forward = (plan.chara == 'N');
    // with several ranks an exchange adds into blocks no local product writes,
    // so every level is zeroed before its first write, as the CPU's are
    const bool zeroing = bf.nproc > 1;
    std::vector<int64_t> acc_n, acc_m;
    leaf_offsets(bf, acc_n, acc_m);

    std::vector<SweepItem> items;
    std::vector<char> written;
    auto level_index = [&](const VecLevel* v) {
        for (size_t i = 0; i < ex.levels.size(); ++i)
            if (ex.levels[i] == v) return static_cast<int>(i);
        return -1;
    };
    // a level takes the buffer the level written before it does not use
    auto take_level = [&](VecLevel* v, int buf) {
        ex.levels.push_back(v);
        ex.buffer.push_back(buf);
        written.push_back(0);
        ex.buffer_elements = std::max(ex.buffer_elements, v->elements);
        return static_cast<int>(ex.levels.size()) - 1;
    };
    auto level_of = [&](VecLevel* v) {
        const int i = level_index(v);
        if (i >= 0) return i;
        return take_level(v, ex.buffer.empty() ? 1 : 3 - ex.buffer.back());
    };
    auto zero_first = [&](int li) {
        if (written[static_cast<size_t>(li)]) return;
        ex.program.push_back(SweepStep{SweepOp::Zero, li, -1, 0});
        written[static_cast<size_t>(li)] = 1;
    };
    // the items since `first` as one batch writing level li (-1: y)
    auto close_batch = [&](size_t first, bool trans_a, int li) {
        if (items.size() == first) return;
        SweepBatch b;
        b.first = static_cast<int>(first);
        b.count = static_cast<int>(items.size() - first);
        b.trans_a = trans_a;
        b.closing = li < 0;
        for (size_t i = first; i < items.size(); ++i) {
            const SweepItem& it = items[i];
            b.max_m = std::max(b.max_m, it.m);
            b.max_n = std::max(b.max_n, it.n);
            b.max_k = std::max(b.max_k, it.k);
            b.flops += 2.0 * it.m * it.n * it.k;
        }
        if (!b.closing) {
            const size_t lu = static_cast<size_t>(li);
            if (written[lu]) {
                b.beta = 1.0;
            } else {
                // the batch writes every block of the level with a nonzero
                // height, once each: no zeroing, beta 0
                const VecLevel& v = *ex.levels[lu];
                int64_t blocks = 0;
                for (int r : v.rows) blocks += r > 0 ? 1 : 0;
                if (!zeroing && b.count == blocks) {
                    b.beta = 0.0;
                    written[lu] = 1;
                } else {
                    zero_first(li);
                    b.beta = 1.0;
                }
            }
        }
        ex.flops += b.flops;
        ex.program.push_back(SweepStep{SweepOp::Launch, static_cast<int>(ex.batches.size()), -1, 0});
        ex.batches.push_back(b);
    };

    const int boundary = forward ? bf.level_half : bf.level_half + 1;
    auto out_level_of = [&](const detail::LevelForm& f, int level) -> VecLevel& {
        if (f.below && level == boundary) return plan.mid_src;
        return plan.vec[static_cast<size_t>(f.out_level)];
    };

    // ---- the leaf product: the factor applied transposed to the input
    auto leaf_product = [&](int leaf_level) {
        const detail::LevelForm f = detail::form_of(bf, forward, leaf_level);
        VecLevel& out = out_level_of(f, leaf_level);
        if (out.elements == 0) return;
        const int li = level_of(&out);
        const LevelDesc& leaf = bf.level(leaf_level);
        const std::vector<int64_t>& acc = forward ? acc_n : acc_m;
        const size_t first = items.size();
        for (int b = 0; b < leaf.nr * leaf.nc; ++b) {
            const BlockDesc& d = bf.blocks[static_cast<size_t>(leaf.first + b)];
            if (d.rows == 0 || d.cols == 0) continue;
            const int index_i = forward ? 1 : b * leaf.inc_r + leaf.idx_r;
            const int index_j = forward ? b * leaf.inc_c + leaf.idx_c : 1;
            const int il = detail::local_of(index_i, out.idx_r, out.inc_r, out.nr, "leaf output row");
            const int jl = detail::local_of(index_j, out.idx_c, out.inc_c, out.nc, "leaf output column");
            const size_t ob = static_cast<size_t>(il) * out.nc + jl;
            SweepItem it;
            it.a = d.ptr;  it.a_off = d.offset;  it.lda = d.lda;
            it.b_buf = 0;  it.b_off = acc[static_cast<size_t>(b)];  it.ldb = 0;  // x and its ld, per call
            it.c_buf = ex.buffer[static_cast<size_t>(li)];  it.c_off = out.offset[ob];  it.ldc = out.rows[ob];
            it.m = d.cols; it.n = nrhs; it.k = d.rows;   // A^T: (cols x rows)
            items.push_back(it);
        }
        close_batch(first, /*trans_a=*/true, li);
    };

    // ---- one kernel level: below the midpoint the pair writes one output, so
    // its two halves go in separate batches; above it the pair writes two
    // outputs but two loop positions share one, so the batches split on the
    // parity of the index that halves.
    auto run_level = [&](int level) {
        const detail::LevelForm f = detail::form_of(bf, forward, level);
        const LevelDesc& kerl = bf.level(level);
        const SweepRanges& g = detail::loop_grid(bf, forward, f.below);
        const size_t lu = static_cast<size_t>(level);
        const VecLevel& in = plan.vec[static_cast<size_t>(f.in_level)];
        VecLevel& out = out_level_of(f, level);
        if (out.elements == 0) return;
        const int in_li = level_index(&in);
        const int li = level_of(&out);
        for (int pass = 0; pass < 2; ++pass) {
            const size_t first = items.size();
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
                        if (in_li < 0) {
                            throw std::runtime_error("bf::gpu::sweep_exec: a level is read before it is written");
                        }
                        SweepItem it;
                        it.a = k.ptr;  it.a_off = k.offset;  it.lda = k.lda;
                        it.b_buf = ex.buffer[static_cast<size_t>(in_li)];  it.b_off = in.offset[ib];
                        it.ldb = in.rows[ib];
                        it.c_buf = ex.buffer[static_cast<size_t>(li)];  it.c_off = out.offset[ob];
                        it.ldc = out.rows[ob];
                        it.m = f.transposed ? k.cols : k.rows;
                        it.n = nrhs;
                        it.k = f.transposed ? k.rows : k.cols;
                        items.push_back(it);
                    }
                }
            }
            close_batch(first, f.transposed, li);
        }
    };

    // ---- the closing product, with the caller's scaling
    auto closing_product = [&](int close_level) {
        const LevelDesc& d = bf.level(close_level);
        const VecLevel& in = plan.vec[static_cast<size_t>(L) + 1];
        const int in_li = level_index(&in);
        const std::vector<int64_t>& acc = forward ? acc_m : acc_n;
        const size_t first = items.size();
        for (int b = 0; b < d.nr * d.nc; ++b) {
            const BlockDesc& k = bf.blocks[static_cast<size_t>(d.first + b)];
            if (k.rows == 0 || k.cols == 0) continue;
            const int index_i = forward ? b * d.inc_r + d.idx_r : 1;
            const int index_j = forward ? 1 : b * d.inc_c + d.idx_c;
            const int il = detail::local_of(index_i, in.idx_r, in.inc_r, in.nr, "closing input row");
            const int jl = detail::local_of(index_j, in.idx_c, in.inc_c, in.nc, "closing input column");
            const size_t ib = static_cast<size_t>(il) * in.nc + jl;
            if (in.rows[ib] == 0) continue;
            if (in_li < 0) throw std::runtime_error("bf::gpu::sweep_exec: the closing level is never written");
            SweepItem it;
            it.a = k.ptr;  it.a_off = k.offset;  it.lda = k.lda;
            it.b_buf = ex.buffer[static_cast<size_t>(in_li)];  it.b_off = in.offset[ib];  it.ldb = in.rows[ib];
            it.c_buf = 3;  it.c_off = acc[static_cast<size_t>(b)];  it.ldc = 0;  // y and its ld, per call
            it.m = k.rows; it.n = nrhs; it.k = k.cols;
            items.push_back(it);
        }
        close_batch(first, false, -1);
    };

    auto exchange_after = [&](int level, int tag) {
        ExchangePlan& xp = plan.pair[static_cast<size_t>(level)];
        if (xp.empty()) return;
        const detail::LevelForm f = detail::form_of(bf, forward, level);
        const int li = level_of(&out_level_of(f, level));
        zero_first(li);  // a rank with no product of its own at this level
        ex.program.push_back(SweepStep{SweepOp::Exchange, level, li, tag});
    };

    // ---- the first half: the leaf product, then its kernel levels, each
    // followed by the pair-wise broadcast.  The level that ends the half has
    // no exchange - the redistribution follows it instead.
    const int step = forward ? 1 : -1;
    const int first_level = forward ? 0 : L + 1;
    const bool full = plan.stop < 0;
    for (int level = first_level;; level += step) {
        if (level == first_level) {
            leaf_product(level);
            if (plan.open_group >= 0) {
                // the group's sum, over the rank's part of the leaf
                const int li = level_of(&out_level_of(detail::form_of(bf, forward, level), level));
                zero_first(li);
                ex.program.push_back(SweepStep{SweepOp::Group, plan.open_group, li, 0});
            }
        } else {
            run_level(level);
        }
        if (level == boundary) break;
        exchange_after(level, 100 + level);
        if (level == plan.stop) break;
    }

    // ---- the midpoint: between the two orderings
    VecLevel& mid_dst = plan.vec[static_cast<size_t>(plan.mid_level)];
    if (!full && plan.stop != boundary) {
        // a partial plan ending before the midpoint
    } else if (!plan.mid.empty()) {
        const int src = level_of(&plan.mid_src);
        zero_first(src);
        const int dst = level_of(&mid_dst);
        zero_first(dst);  // the arrivals add into it
        ex.program.push_back(SweepStep{SweepOp::Midpoint, 0, -1, 200});
    } else if (plan.mid_src.elements > 0 && plan.mid_level >= 0) {
        // one rank: the orderings coincide, so the level is already in place
        const int src = level_index(&plan.mid_src);
        if (src >= 0 && mid_dst.elements == plan.mid_src.elements && mid_dst.offset == plan.mid_src.offset) {
            const int dst = take_level(&mid_dst, ex.buffer[static_cast<size_t>(src)]);
            written[static_cast<size_t>(dst)] = 1;
        } else if (src >= 0 && mid_dst.elements == plan.mid_src.elements) {
            const int dst = level_of(&mid_dst);
            written[static_cast<size_t>(dst)] = 1;
            ex.program.push_back(SweepStep{SweepOp::Copy, dst, src, 0});
        }
    }

    // ---- the second half: its kernel levels, each followed by the pair-wise
    // reduction, then the closing product
    const int last = forward ? L + 1 : 0;
    for (int level = boundary + step; full; level += step) {
        if (level == last) {
            if (plan.close_group >= 0) {
                // the head's block to the group's other ranks, which no
                // product of theirs writes
                const int li = level_of(&plan.vec[static_cast<size_t>(L) + 1]);
                zero_first(li);
                ex.program.push_back(SweepStep{SweepOp::Group, plan.close_group, li, 0});
            }
            closing_product(level);
            break;
        }
        run_level(level);
        exchange_after(level, 300 + level);
    }

    // ---- the tensor-core tiles of the batches that run on them: those no
    // narrow kernel takes (the closing product's narrow kernel needs the
    // caller's alpha 1 and beta 0 or 1, so it keeps its tiles)
    std::vector<int2> tiles;
    const bool tc = std::is_same<DataType, double>::value;
    for (SweepBatch& b : ex.batches) {
        b.tc_first = static_cast<int>(tiles.size());
        b.tc32 = b.tc48 = b.tc64 = 0;
        if (!tc || !sweep_tensor_cores<T>(b.trans_a, b.max_m, b.max_n, b.max_k)) continue;
        if (!b.closing && narrow_kind<T>(b.trans_a, b.max_m, b.max_n) != 0) continue;
        std::vector<int2> t32, t48, t64;
        for (int i = 0; i < b.count; ++i) {
            const SweepItem& it = items[static_cast<size_t>(b.first + i)];
            if (it.m <= 0 || it.n <= 0) continue;
            typename fmm::gpu::VBatch<double>::Entry e{};
            e.m = it.m;
            e.n = it.n;
            const int tile = b.max_n <= kTcNarrowCols ? 32 : fmm::gpu::VBatch<double>::tc_tile_of(e);
            std::vector<int2>& out = tile == 32 ? t32 : (tile == 48 ? t48 : t64);
            const int rows = (it.m + tile - 1) / tile, cols = (it.n + tile - 1) / tile;
            for (int c = 0; c < cols; ++c)
                for (int r = 0; r < rows; ++r) out.push_back(make_int2(i, (r << 16) | c));
        }
        b.tc32 = static_cast<int>(t32.size());
        b.tc48 = static_cast<int>(t48.size());
        b.tc64 = static_cast<int>(t64.size());
        tiles.insert(tiles.end(), t32.begin(), t32.end());
        tiles.insert(tiles.end(), t48.begin(), t48.end());
        tiles.insert(tiles.end(), t64.begin(), t64.end());
    }

    // ---- the arrays, in one block of the arena and one copy
    static_assert(sizeof(magma_int_t) == sizeof(int), "the tensor-core GEMM takes 32-bit sizes");
    const size_t n = items.size();
    ex.count = static_cast<int>(n);
    size_t at = 0;
    auto place = [&](size_t bytes) {
        const size_t here = at;
        at = fmm::gpu::align_up(at + std::max<size_t>(bytes, 1));
        return here;
    };
    const size_t o_a = place(n * sizeof(DataType*)), o_b = place(n * sizeof(DataType*));
    const size_t o_c = place(n * sizeof(DataType*));
    const size_t o_sizes = place(6 * (n + 1) * sizeof(magma_int_t));
    const size_t o_boff = place(n * sizeof(int64_t)), o_coff = place(n * sizeof(int64_t));
    const size_t o_bbuf = place(n), o_cbuf = place(n);
    const size_t o_tc = place(tiles.size() * sizeof(int2));
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    ex.block = heap.try_alloc(at, fmm::gpu::ResidentKind::Matvec);
    if (ex.block == nullptr) throw std::runtime_error("bf::gpu::sweep_exec: no room on the device for the sweep's arrays");
    std::vector<char> image(at, 0);
    const DataType** ha = reinterpret_cast<const DataType**>(image.data() + o_a);
    magma_int_t* hs = reinterpret_cast<magma_int_t*>(image.data() + o_sizes);
    int64_t* hboff = reinterpret_cast<int64_t*>(image.data() + o_boff);
    int64_t* hcoff = reinterpret_cast<int64_t*>(image.data() + o_coff);
    signed char* hbbuf = reinterpret_cast<signed char*>(image.data() + o_bbuf);
    signed char* hcbuf = reinterpret_cast<signed char*>(image.data() + o_cbuf);
    for (size_t i = 0; i < n; ++i) {
        const SweepItem& it = items[i];
        ha[i] = it.a != nullptr ? static_cast<const DataType*>(it.a) : bf.values + it.a_off;
        hs[0 * (n + 1) + i] = it.m;
        hs[1 * (n + 1) + i] = it.n;
        hs[2 * (n + 1) + i] = it.k;
        hs[3 * (n + 1) + i] = it.lda;
        hs[4 * (n + 1) + i] = it.ldb;
        hs[5 * (n + 1) + i] = it.ldc;
        hboff[i] = it.b_off;
        hcoff[i] = it.c_off;
        hbbuf[i] = static_cast<signed char>(it.b_buf);
        hcbuf[i] = static_cast<signed char>(it.c_buf);
    }
    if (!tiles.empty()) std::memcpy(image.data() + o_tc, tiles.data(), tiles.size() * sizeof(int2));
    fmm::gpu::check_cuda(cudaMemcpy(ex.block, image.data(), at, cudaMemcpyHostToDevice), "bf sweep record");
    ex.d_a = reinterpret_cast<const DataType**>(ex.block + o_a);
    ex.d_b = reinterpret_cast<const DataType**>(ex.block + o_b);
    ex.d_c = reinterpret_cast<DataType**>(ex.block + o_c);
    magma_int_t* ds = reinterpret_cast<magma_int_t*>(ex.block + o_sizes);
    ex.d_m = ds;
    ex.d_n = ds + (n + 1);
    ex.d_k = ds + 2 * (n + 1);
    ex.d_lda = ds + 3 * (n + 1);
    ex.d_ldb = ds + 4 * (n + 1);
    ex.d_ldc = ds + 5 * (n + 1);
    ex.d_boff = reinterpret_cast<int64_t*>(ex.block + o_boff);
    ex.d_coff = reinterpret_cast<int64_t*>(ex.block + o_coff);
    ex.d_bbuf = reinterpret_cast<signed char*>(ex.block + o_bbuf);
    ex.d_cbuf = reinterpret_cast<signed char*>(ex.block + o_cbuf);
    ex.d_tc = reinterpret_cast<int2*>(ex.block + o_tc);

    // a caller alternating among a few nrhs keeps them all; more than that
    // and the oldest go
    constexpr size_t kKeep = 4;
    while (plan.execs.size() >= kKeep) plan.execs.erase(plan.execs.begin());
    SweepExec& ref = *owned;
    plan.execs[nrhs] = std::move(owned);
    return ref;
}

// One batch of a replay: a narrow kernel for few columns, the tensor cores
// for many, else MAGMA - all three sum every element in the same order.
template<typename T>
inline void run_batch(const SweepExec<T>& ex, const SweepBatch& b, T alpha, T beta, cudaStream_t stream) {
    using DataType = T;
    const DataType al = b.closing ? alpha : DataType(1);
    const DataType be = b.closing ? beta : DataType(b.beta);
    const size_t f = static_cast<size_t>(b.first);
    if constexpr (std::is_same<DataType, double>::value) {
        if (al == DataType(1) && (be == DataType(0) || be == DataType(1)) &&
            sweep_narrow_kind(b.trans_a, b.max_m, b.max_n) != 0) {
            launch_sweep_narrow(b.trans_a, b.count, ex.d_m + f, ex.d_n + f, ex.d_k + f, al, ex.d_a + f,
                                ex.d_lda + f, ex.d_b + f, ex.d_ldb + f, be, ex.d_c + f, ex.d_ldc + f, b.max_m,
                                b.max_n, stream);
            return;
        }
    }
    if (b.tc32 + b.tc48 + b.tc64 > 0) {
        if constexpr (std::is_same<DataType, double>::value) {
            const int parts[3] = {b.tc32, b.tc48, b.tc64};
            int at = b.tc_first;
            for (int part = 0; part < 3; ++part) {
                if (parts[part] > 0) {
                    fmm::gpu::launch_dgemm_vbatched_tc(b.trans_a, false, ex.d_m + f, ex.d_n + f, ex.d_k + f, al,
                                                       ex.d_a + f, ex.d_lda + f, ex.d_b + f, ex.d_ldb + f, be,
                                                       ex.d_c + f, ex.d_ldc + f, ex.d_tc + at, parts[part],
                                                       32 + 16 * part, stream);
                }
                at += parts[part];
            }
            return;
        }
    }
    fmm::gpu::gemm_vbatched<DataType>(b.trans_a ? MagmaTrans : MagmaNoTrans, MagmaNoTrans, ex.d_m + f, ex.d_n + f,
                                      ex.d_k + f, al, ex.d_a + f, ex.d_lda + f, ex.d_b + f, ex.d_ldb + f, be,
                                      ex.d_c + f, ex.d_ldc + f, static_cast<magma_int_t>(b.count),
                                      fmm::gpu::Context::instance().queue(), b.max_m, b.max_n, b.max_k);
}

namespace detail {

// A replay of ex: its two buffers from the arena, the kernel that turns the
// GEMMs' offsets into pointers, then the steps.
template<typename T>
inline void replay(SweepPlan<T>& plan, SweepExec<T>& ex, int nrhs, const T* vin, int ldi, T* vout, int ldo, T alpha,
                   T beta, T* buf[3], cudaStream_t stream) {
    using DataType = T;
    using VecLevel = bf::gpu::VecLevel<T>;
    using SweepExec = bf::gpu::SweepExec<T>;
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    buf[0] = buf[1] = buf[2] = nullptr;
    if (ex.buffer_elements > 0) {
        const size_t bytes = static_cast<size_t>(ex.buffer_elements) * sizeof(DataType);
        buf[1] = heap.alloc<DataType>(bytes, fmm::gpu::ResidentKind::Matvec);
        buf[2] = heap.alloc<DataType>(bytes, fmm::gpu::ResidentKind::Matvec);
    }
    for (size_t i = 0; i < ex.levels.size(); ++i) ex.levels[i]->data = buf[ex.buffer[i]];
    launch_sweep_patch(ex.count, ex.d_boff, ex.d_coff, ex.d_bbuf, ex.d_cbuf, vin, buf[1], buf[2], vout, ldi, ldo,
                       ex.d_b, ex.d_c, ex.d_ldb, ex.d_ldc, stream);

    for (const SweepStep& st : ex.program) {
        switch (st.op) {
            case SweepOp::Launch:
                run_batch(ex, ex.batches[static_cast<size_t>(st.arg)], alpha, beta, stream);
                break;
            case SweepOp::Zero: {
                const VecLevel& v = *ex.levels[static_cast<size_t>(st.arg)];
                if (v.elements > 0) {
                    fmm::gpu::check_cuda(cudaMemsetAsync(v.data, 0, static_cast<size_t>(v.elements) * sizeof(DataType),
                                                         stream), "bf sweep zero");
                }
                break;
            }
            case SweepOp::Exchange: {
                VecLevel& v = *ex.levels[static_cast<size_t>(st.level)];
                exchange_level(plan.pair[static_cast<size_t>(st.arg)], v, v, nrhs, st.tag, stream);
                break;
            }
            case SweepOp::Midpoint:
                exchange_level(plan.mid, plan.mid_src, plan.vec[static_cast<size_t>(plan.mid_level)], nrhs, st.tag,
                               stream);
                break;
            case SweepOp::Copy: {
                const VecLevel& dst = *ex.levels[static_cast<size_t>(st.arg)];
                const VecLevel& src = *ex.levels[static_cast<size_t>(st.level)];
                fmm::gpu::check_cuda(cudaMemcpyAsync(dst.data, src.data, static_cast<size_t>(dst.elements) * sizeof(DataType),
                                                     cudaMemcpyDeviceToDevice, stream), "bf sweep midpoint");
                break;
            }
            case SweepOp::Group: {
                const GroupPlan& g = plan.groups[static_cast<size_t>(st.arg)];
                VecLevel& v = *ex.levels[static_cast<size_t>(st.level)];
                group_level(g, v.data + v.offset[g.block], static_cast<int64_t>(v.rows[g.block]) * nrhs, stream);
                break;
            }
        }
    }
}

}  // namespace detail

// The sweep.  Vout = alpha * op(bf) * Vin + beta * Vout for 'N', and
// Vout = alpha * op(bf)^T * Vin + beta * Vout for 'T' and 'C'; every pointer
// device resident.  `plan` must have been built for this block and this
// orientation; the first call for an nrhs records the sweep (sweep_exec).
//
// The caller's alpha and beta belong to the closing product alone.  The
// stream must be the backend's (fmm::gpu::Context), on which MAGMA runs.
template<typename T>
inline double sweep(const BfDevice<T>& bf, SweepPlan<T>& plan, int nrhs, const T* vin, int ldi, T* vout, int ldo,
                    T alpha, T beta, SweepScratch<T>& /*unused*/, cudaStream_t stream = nullptr) {
    using DataType = T;
    using VecLevel = bf::gpu::VecLevel<T>;
    if (!plan.built || plan.vec.empty() || plan.stop >= 0) {
        throw std::invalid_argument("bf::gpu::sweep: the plan is not built, or is a partial one");
    }
    if (stream == nullptr) stream = fmm::gpu::Context::instance().stream();
    SweepExec<T>& ex = sweep_exec(bf, plan, nrhs);
    repack(plan, nrhs);  // the levels' offsets, which the exchanges read
    DataType* buf[3];
    detail::replay<T>(plan, ex, nrhs, vin, ldi, vout, ldo, alpha, beta, buf, stream);

    // the buffers go back with the work still queued: whatever reuses them
    // runs after it, on the same stream
    auto& heap = fmm::gpu::DeviceHeap::instance();
    for (int b = 1; b <= 2; ++b)
        if (buf[b] != nullptr) heap.free(buf[b]);
    for (VecLevel& v : plan.vec) v.data = nullptr;
    plan.mid_src.data = nullptr;
    return ex.flops;
}

// The sweep of a partial plan (build_sweep_plan's stop) over x (N_loc x nrhs
// forward, M_loc x nrhs transposed, ld ldi): the level it ends in,
// plan.vec[plan.result], whose data is a buffer of the arena the caller gives
// back (DeviceHeap::free) once the work that reads it is queued.  flops: what
// the sweep did.
template<typename T>
inline VecLevel<T>& sweep_partial(const BfDevice<T>& bf, SweepPlan<T>& plan, int nrhs, const T* vin, int ldi,
                                  double* flops = nullptr, cudaStream_t stream = nullptr) {
    using DataType = T;
    using VecLevel = bf::gpu::VecLevel<T>;
    if (!plan.built || plan.vec.empty() || plan.stop < 0) {
        throw std::invalid_argument("bf::gpu::sweep_partial: the plan is not built, or is a full one");
    }
    if (stream == nullptr) stream = fmm::gpu::Context::instance().stream();
    SweepExec<T>& ex = sweep_exec(bf, plan, nrhs);
    repack(plan, nrhs);
    DataType* buf[3];
    detail::replay<T>(plan, ex, nrhs, vin, ldi, nullptr, 0, DataType(1), DataType(0), buf, stream);
    VecLevel& out = plan.vec[static_cast<size_t>(plan.result)];
    DataType* kept = out.elements > 0 ? out.data : nullptr;
    if (out.elements > 0 && kept == nullptr) {
        throw std::runtime_error("bf::gpu::sweep_partial: the plan's last level is never written");
    }
    auto& heap = fmm::gpu::DeviceHeap::instance();
    for (int b = 1; b <= 2; ++b)
        if (buf[b] != nullptr && buf[b] != kept) heap.free(buf[b]);
    for (VecLevel& v : plan.vec) v.data = nullptr;
    plan.mid_src.data = nullptr;
    out.data = kept;
    if (flops != nullptr) *flops = ex.flops;
    return out;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
