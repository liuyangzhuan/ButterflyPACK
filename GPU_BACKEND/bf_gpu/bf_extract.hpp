#pragma once

#ifdef H2_HAVE_GPU

// Entry extraction from a device-resident butterfly (GPU_BACKEND/bf_gpu): the
// blocks BF(I_k, J_k) of a list of intersections, as BF_block_extraction
// computes them on the host, with the bookkeeping in arrays and each level's
// products one batched GEMM.
//
// Why not the host's structure: measured on the requests a butterfly's own
// construction makes (thousands of small intersections a call, BF2 from BF1
// in EXAMPLE/Butterfly_Driver.cpp), three quarters of BF_block_extraction's
// time is its linked-list bookkeeping, not its GEMMs.  So moving the GEMMs
// alone would have bought little.
//
// Every entry follows one path through the factors.  With row i in leaf a and
// column j in leaf b, L = level_butterfly, and Kerl(t) the 2^t x 2^(L-t+1)
// grid of level t (0-based block indices),
//
//   BF(i, j) = U_a(i,:) Kerl(L)(a, b>>(L-1)) ... Kerl(t)(a>>(L-t), b>>(t-1)) ... Kerl(1)(a>>(L-1), b) V_b(j,:)^T
//
// so the product meets in the middle, at h = level_half, as the host's does:
// the right half depends on j and on a's group at level h only, the left half
// on i and on b's group at level h only.  A level's partial products are
// therefore kept per (intersection, row group, column group), each a matrix
// over the intersection's columns (right half) or rows (left half) that lie
// in its groups.  Within an intersection the rows and the columns are put in
// leaf order first (a counting sort), so every group is a contiguous run and
// the two children of a group split it in two: a product writes a sub-block of
// its parent's matrix and no two products of a level write the same one, so a
// level is one batch with beta 0.  The entries come out in leaf order and a
// final scatter puts them in the caller's.
//
// The ranks follow the vector levels of the 'N' sweep: level t's block
// (F, C), F < 2^(t-1), C < 2^(L-t+1), has the rank of V_C's columns at t = 1
// and of Kerl(t-1)(F, 2C)'s rows above (U is Kerl(L+1)).
//
// One process only: a butterfly on several ranks declines (the caller runs
// the host path), as does a request any of whose blocks another rank owns.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "bf_device.hpp"
#include "../h2_gpu/gpu_runtime.hpp"

namespace bf {
namespace gpu {

// One gather of factor rows into a partial product: the right half's first
// level takes V_b's rows transposed (r x count), the left half's U_a's rows
// as they are (count x r).
struct ExtractGather {
    int64_t src = 0;     // the factor block, offset in BfDevice::values
    int64_t idx = 0;     // its local rows to take, offset in the index array
    int64_t dst = 0;     // offset in the partial buffer
    int ld_src = 0;      // rows of the factor block
    int count = 0;       // rows taken
    int r = 0;           // columns of the factor block
    int transposed = 0;  // 1: dst is r x count; 0: count x r
};

// One intersection's block, from leaf order to the caller's.
struct ExtractScatter {
    int64_t src = 0, dst = 0;     // both m x n, column major
    int64_t pos_r = 0, pos_c = 0;  // offsets in the index array: each caller row's (column's) leaf-order position
    int m = 0, n = 0;
};

void launch_extract_gather(const ExtractGather* items, int count, const int* idx, const DataType* values,
                           DataType* dst, cudaStream_t stream);
void launch_extract_scatter(const ExtractScatter* items, int count, const int* idx, const DataType* src,
                            DataType* dst, cudaStream_t stream);

// The leaves of a resident butterfly on one process: which leaf each row
// and column is in and where each leaf starts.  Made at the first extraction
// and kept with the factors.
struct ExtractLeaves {
    std::vector<int> row_leaf, col_leaf;
    std::vector<int> row_head, col_head;  // 2^L + 1 each
    bool ready = false;
};

// Seconds and counts of one extraction, for BPACK_TRACE=phase.
struct ExtractProfile {
    double plan = 0, run = 0, copy = 0;
    int64_t items = 0, entries = 0;
    double flops = 0;
};

namespace detail {

// A run of a sorted list: the elements [start, start + count) are in group g.
struct Run {
    int g = 0, start = 0, count = 0;
};

// The runs of lf >> shift over a list sorted by leaf.
inline void runs_of(const int* lf, int n, int shift, std::vector<Run>& out) {
    out.clear();
    for (int p = 0; p < n;) {
        const int g = lf[p] >> shift;
        int q = p + 1;
        while (q < n && (lf[q] >> shift) == g) ++q;
        out.push_back(Run{g, p, q - p});
        p = q;
    }
}

// The index of group g among runs (sorted by group).
inline int run_index(const std::vector<Run>& runs, int g) {
    int lo = 0, hi = static_cast<int>(runs.size()) - 1;
    while (lo <= hi) {
        const int mid = (lo + hi) / 2;
        if (runs[mid].g == g) return mid;
        if (runs[mid].g < g) lo = mid + 1;
        else hi = mid - 1;
    }
    throw std::logic_error("bf::gpu::extract: a parent group is missing");
}

}  // namespace detail

inline void prepare_leaves(const BfDevice& bf, ExtractLeaves& lv) {
    if (lv.ready) return;
    const int L = bf.level_butterfly;
    const int nleaf = 1 << L;
    const LevelDesc& dv = bf.level(0);
    const LevelDesc& du = bf.level(L + 1);
    if (dv.nr != 1 || dv.nc != nleaf || du.nr != nleaf || du.nc != 1) {
        throw std::invalid_argument("bf::gpu::extract: the factors are not all on this process");
    }
    lv.col_head.assign(static_cast<size_t>(nleaf) + 1, 0);
    lv.row_head.assign(static_cast<size_t>(nleaf) + 1, 0);
    for (int b = 0; b < nleaf; ++b) lv.col_head[b + 1] = lv.col_head[b] + bf.block(0, 0, b).rows;
    for (int a = 0; a < nleaf; ++a) lv.row_head[a + 1] = lv.row_head[a] + bf.block(L + 1, a, 0).rows;
    if (lv.col_head[nleaf] != bf.N || lv.row_head[nleaf] != bf.M) {
        throw std::invalid_argument("bf::gpu::extract: the leaves do not cover the block");
    }
    lv.col_leaf.resize(static_cast<size_t>(bf.N));
    lv.row_leaf.resize(static_cast<size_t>(bf.M));
    for (int b = 0; b < nleaf; ++b)
        for (int j = lv.col_head[b]; j < lv.col_head[b + 1]; ++j) lv.col_leaf[j] = b;
    for (int a = 0; a < nleaf; ++a)
        for (int i = lv.row_head[a]; i < lv.row_head[a + 1]; ++i) lv.row_leaf[i] = a;
    lv.ready = true;
}

// BF(I_k, J_k) for k < count: rowidx[k] rows from `rows` and colidx[k]
// columns from `cols` (0-based, the butterfly's own order, stacked), into
// `out` (host), the blocks stacked in order, each column major.  Returns
// false when the device has no room for the partial products (nothing
// written; the caller runs the host path).
inline bool extract(const BfDevice& bf, ExtractLeaves& leaves, int count, const int* rowidx, const int* colidx,
                    const int* rows, const int* cols, DataType* out, ExtractProfile& prof,
                    fmm::gpu::MetaBuilder& meta, fmm::gpu::DeviceBuffer& meta_device) {
    using clock = std::chrono::steady_clock;
    using detail::Run;
    const auto t0 = clock::now();
    prepare_leaves(bf, leaves);
    const int L = bf.level_butterfly;
    const int h = std::min(std::max(bf.level_half, 0), L);

    // the rank of vector level t's block (F, C)
    auto rank_of = [&](int t, int F, int C) -> int {
        if (t == 1) return bf.block(0, 0, C).cols;
        return bf.block(t - 1, F, 2 * C).rows;
    };
    auto factor = [&](int l, int i, int j) -> const BlockDesc& { return bf.block(l, i, j); };

    // ---- per intersection: the rows and columns in leaf order, and the runs
    // of their groups at every level
    std::vector<int64_t> row_at(static_cast<size_t>(count) + 1, 0), col_at(static_cast<size_t>(count) + 1, 0);
    for (int k = 0; k < count; ++k) {
        row_at[k + 1] = row_at[k] + rowidx[k];
        col_at[k + 1] = col_at[k] + colidx[k];
    }
    // leaf order: srow[p'] is the caller's row at leaf-order position p'
    std::vector<int> srow(static_cast<size_t>(row_at[count])), scol(static_cast<size_t>(col_at[count]));
    std::vector<int> rleaf(srow.size()), cleaf(scol.size());
    std::vector<int> pos_r(srow.size()), pos_c(scol.size());
    for (int k = 0; k < count; ++k) {
        auto order = [&](const int* idx, int64_t at, int n, const std::vector<int>& leaf_of, std::vector<int>& sorted,
                         std::vector<int>& lf, std::vector<int>& pos) {
            std::vector<int> p(static_cast<size_t>(n));
            std::iota(p.begin(), p.end(), 0);
            std::stable_sort(p.begin(), p.end(), [&](int x, int y) {
                return leaf_of[idx[at + x]] < leaf_of[idx[at + y]];
            });
            for (int q = 0; q < n; ++q) {
                sorted[at + q] = idx[at + p[q]];
                lf[at + q] = leaf_of[sorted[at + q]];
                pos[at + p[q]] = q;
            }
        };
        order(rows, row_at[k], rowidx[k], leaves.row_leaf, srow, rleaf, pos_r);
        order(cols, col_at[k], colidx[k], leaves.col_leaf, scol, cleaf, pos_c);
    }
    // runs[k][t]: the row groups (shift L - t) and the column groups (shift t)
    // of level t, t = 0 .. L
    std::vector<std::vector<std::vector<Run>>> rrun(static_cast<size_t>(count)), crun(static_cast<size_t>(count));
    for (int k = 0; k < count; ++k) {
        rrun[k].resize(static_cast<size_t>(L) + 1);
        crun[k].resize(static_cast<size_t>(L) + 1);
        for (int t = 0; t <= L; ++t) {
            detail::runs_of(rleaf.data() + row_at[k], rowidx[k], L - t, rrun[k][t]);
            detail::runs_of(cleaf.data() + col_at[k], colidx[k], t, crun[k][t]);
        }
    }

    // ---- the partial products: offsets per (level, intersection, row run,
    // column run).  Right half: vector levels 1 .. h+1, matrices
    // rank x (columns in the run); left half: L+1 down to h+1, (rows in the
    // run) x rank.  Level t's groups are those of tree level t-1.
    struct Slots {
        std::vector<std::vector<int64_t>> at;  // [k][fi * nC + ci]
        int64_t total = 0;
    };
    auto layout = [&](int t, bool right, Slots& s) {
        s.at.assign(static_cast<size_t>(count), {});
        s.total = 0;
        for (int k = 0; k < count; ++k) {
            const auto& fr = rrun[k][t - 1];
            const auto& cr = crun[k][t - 1];
            auto& at = s.at[k];
            at.resize(fr.size() * cr.size());
            for (size_t fi = 0; fi < fr.size(); ++fi) {
                for (size_t ci = 0; ci < cr.size(); ++ci) {
                    const int r = rank_of(t, fr[fi].g, cr[ci].g);
                    at[fi * cr.size() + ci] = s.total;
                    s.total += static_cast<int64_t>(r) * (right ? cr[ci].count : fr[fi].count);
                }
            }
        }
    };
    std::vector<Slots> R(static_cast<size_t>(h) + 2), W(static_cast<size_t>(L) + 2);
    for (int t = 1; t <= h + 1; ++t) layout(t, true, R[t]);
    for (int t = L + 1; t >= h + 1; --t) layout(t, false, W[t]);

    // the device buffers: one per level and half, the leaf-order entries, the
    // caller-order entries
    int64_t entries = 0;
    std::vector<int64_t> out_at(static_cast<size_t>(count) + 1, 0);
    for (int k = 0; k < count; ++k) out_at[k + 1] = out_at[k] + static_cast<int64_t>(rowidx[k]) * colidx[k];
    entries = out_at[count];
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    std::vector<char*> owned;
    bool room = true;
    auto take = [&](int64_t elements) -> DataType* {
        if (!room) return nullptr;
        char* p = heap.try_alloc(static_cast<size_t>(std::max<int64_t>(elements, 1)) * sizeof(DataType),
                                 fmm::gpu::ResidentKind::Matvec);
        if (p == nullptr) { room = false; return nullptr; }
        owned.push_back(p);
        return reinterpret_cast<DataType*>(p);
    };
    std::vector<DataType*> dR(static_cast<size_t>(h) + 2, nullptr), dW(static_cast<size_t>(L) + 2, nullptr);
    for (int t = 1; t <= h + 1; ++t) dR[t] = take(R[t].total);
    for (int t = L + 1; t >= h + 1; --t) dW[t] = take(W[t].total);
    DataType* sorted = take(entries);
    DataType* d_out = take(entries);
    if (!room) {
        for (char* p : owned) heap.free(p);
        return false;
    }

    // ---- the work, all staged in one metadata image
    meta.clear();
    std::vector<ExtractGather> gathers;
    std::vector<int> index;
    // V: level 1 of the right half, one row group (all rows) per intersection
    for (int k = 0; k < count; ++k) {
        const auto& cr = crun[k][0];
        for (size_t ci = 0; ci < cr.size(); ++ci) {
            const int b = cr[ci].g;
            const BlockDesc& v = factor(0, 0, b);
            ExtractGather g;
            g.src = v.offset;
            g.ld_src = v.lda;
            g.r = v.cols;
            g.count = cr[ci].count;
            g.transposed = 1;
            g.dst = R[1].at[k][ci];
            g.idx = static_cast<int64_t>(index.size());
            for (int q = 0; q < cr[ci].count; ++q)
                index.push_back(scol[col_at[k] + cr[ci].start + q] - leaves.col_head[b]);
            gathers.push_back(g);
        }
    }
    const int nv = static_cast<int>(gathers.size());
    // U: level L+1 of the left half, one column group (all columns)
    for (int k = 0; k < count; ++k) {
        const auto& fr = rrun[k][L];
        for (size_t fi = 0; fi < fr.size(); ++fi) {
            const int a = fr[fi].g;
            const BlockDesc& u = factor(L + 1, a, 0);
            ExtractGather g;
            g.src = u.offset;
            g.ld_src = u.lda;
            g.r = u.cols;
            g.count = fr[fi].count;
            g.transposed = 0;
            g.dst = W[L + 1].at[k][fi];
            g.idx = static_cast<int64_t>(index.size());
            for (int q = 0; q < fr[fi].count; ++q)
                index.push_back(srow[row_at[k] + fr[fi].start + q] - leaves.row_head[a]);
            gathers.push_back(g);
        }
    }
    const int nu = static_cast<int>(gathers.size()) - nv;

    // GEMM batches, at most kMaxBatch entries each
    constexpr size_t kMaxBatch = 65535;
    std::vector<fmm::gpu::VBatch<DataType>> batches;
    std::vector<int> batch_level;  // which stage a batch belongs to, in launch order
    auto add = [&](int stage, const DataType* a, int lda, const DataType* b, int ldb, DataType* c, int ldc, int m,
                   int n, int kk) {
        if (m <= 0 || n <= 0) return;
        if (batches.empty() || batch_level.back() != stage || batches.back().entries.size() >= kMaxBatch) {
            batches.emplace_back();
            batch_level.push_back(stage);
        }
        typename fmm::gpu::VBatch<DataType>::Entry e;
        e.a = a; e.b = b; e.c = c;
        e.m = m; e.n = n; e.k = kk;
        e.lda = std::max(lda, 1); e.ldb = std::max(ldb, 1); e.ldc = std::max(ldc, 1);
        batches.back().entries.push_back(e);
        prof.flops += 2.0 * m * n * kk;
    };
    // right half: vector level t -> t+1 through Kerl(t), t = 1 .. h
    for (int t = 1; t <= h; ++t) {
        for (int k = 0; k < count; ++k) {
            const auto& f_in = rrun[k][t - 1];
            const auto& f_out = rrun[k][t];
            const auto& c_in = crun[k][t - 1];
            const auto& c_out = crun[k][t];
            for (size_t fo = 0; fo < f_out.size(); ++fo) {
                const int F = f_out[fo].g;
                const int fi = detail::run_index(f_in, F >> 1);
                for (size_t ci = 0; ci < c_in.size(); ++ci) {
                    const int C = c_in[ci].g;
                    const int co = detail::run_index(c_out, C >> 1);
                    const BlockDesc& K = factor(t, F, C);
                    const int r_out = rank_of(t + 1, F, C >> 1);
                    const int r_in = rank_of(t, F >> 1, C);
                    const DataType* B = dR[t] + R[t].at[k][fi * c_in.size() + ci];
                    DataType* Cp = dR[t + 1] + R[t + 1].at[k][fo * c_out.size() + co] +
                                   static_cast<int64_t>(c_in[ci].start - c_out[co].start) * r_out;
                    add(t, bf.values + K.offset, K.lda, B, r_in, Cp, r_out, r_out, c_in[ci].count, r_in);
                }
            }
        }
    }
    // left half: vector level t -> t-1 through Kerl(t-1), t = L+1 .. h+2
    for (int t = L + 1; t >= h + 2; --t) {
        for (int k = 0; k < count; ++k) {
            const auto& f_in = rrun[k][t - 1];
            const auto& f_out = rrun[k][t - 2];
            const auto& c_in = crun[k][t - 1];
            const auto& c_out = crun[k][t - 2];
            for (size_t fi = 0; fi < f_in.size(); ++fi) {
                const int F = f_in[fi].g;
                const int fo = detail::run_index(f_out, F >> 1);
                for (size_t co = 0; co < c_out.size(); ++co) {
                    const int J = c_out[co].g;
                    const int ci = detail::run_index(c_in, J >> 1);
                    const BlockDesc& K = factor(t - 1, F, J);
                    const int r_in = rank_of(t, F, J >> 1);
                    const int r_out = rank_of(t - 1, F >> 1, J);
                    const DataType* A = dW[t] + W[t].at[k][fi * c_in.size() + ci];
                    DataType* Cp = dW[t - 1] + W[t - 1].at[k][fo * c_out.size() + co] +
                                   (f_in[fi].start - f_out[fo].start);
                    add(1000 + t, A, f_in[fi].count, bf.values + K.offset, K.lda, Cp, f_out[fo].count,
                        f_in[fi].count, r_out, r_in);
                }
            }
        }
    }
    // the middle: entries in leaf order
    for (int k = 0; k < count; ++k) {
        const auto& fr = rrun[k][h];
        const auto& cr = crun[k][h];
        const int m = rowidx[k];
        for (size_t fi = 0; fi < fr.size(); ++fi) {
            for (size_t ci = 0; ci < cr.size(); ++ci) {
                const int r = rank_of(h + 1, fr[fi].g, cr[ci].g);
                const DataType* A = dW[h + 1] + W[h + 1].at[k][fi * cr.size() + ci];
                const DataType* B = dR[h + 1] + R[h + 1].at[k][fi * cr.size() + ci];
                DataType* Cp = sorted + out_at[k] + fr[fi].start + static_cast<int64_t>(cr[ci].start) * m;
                add(2000, A, fr[fi].count, B, r, Cp, m, fr[fi].count, cr[ci].count, r);
            }
        }
    }
    // the scatter to the caller's order
    std::vector<ExtractScatter> scatters(static_cast<size_t>(count));
    for (int k = 0; k < count; ++k) {
        ExtractScatter& s = scatters[k];
        s.src = s.dst = out_at[k];
        s.m = rowidx[k];
        s.n = colidx[k];
        s.pos_r = static_cast<int64_t>(index.size());
        index.insert(index.end(), pos_r.begin() + row_at[k], pos_r.begin() + row_at[k + 1]);
        s.pos_c = static_cast<int64_t>(index.size());
        index.insert(index.end(), pos_c.begin() + col_at[k], pos_c.begin() + col_at[k + 1]);
    }

    const size_t at_gather = meta.append(gathers);
    const size_t at_scatter = meta.append(scatters);
    const size_t at_index = meta.append(index);
    for (auto& b : batches) b.stage(meta);
    for (const auto& b : batches) prof.items += static_cast<int64_t>(b.count());
    const auto t1 = clock::now();

    // ---- launch, in dependence order on the backend's stream
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    magma_queue_t queue = fmm::gpu::Context::instance().queue();
    char* dm = meta.upload(meta_device, stream);
    const auto* d_gather = reinterpret_cast<const ExtractGather*>(dm + at_gather);
    const auto* d_scatter = reinterpret_cast<const ExtractScatter*>(dm + at_scatter);
    const auto* d_index = reinterpret_cast<const int*>(dm + at_index);
    // a rank of zero leaves its products unwritten: the entries start at zero
    fmm::gpu::check_cuda(cudaMemsetAsync(sorted, 0, static_cast<size_t>(std::max<int64_t>(entries, 1)) * sizeof(DataType),
                                         stream), "bf extract zero");
    if (nv > 0) launch_extract_gather(d_gather, nv, d_index, bf.values, dR[1], stream);
    if (nu > 0) launch_extract_gather(d_gather + nv, nu, d_index, bf.values, dW[L + 1], stream);
    for (const auto& b : batches) b.gemm(dm, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(0), queue);
    if (count > 0) launch_extract_scatter(d_scatter, count, d_index, sorted, d_out, stream);
    if (entries > 0) {
        fmm::gpu::check_cuda(cudaMemcpyAsync(out, d_out, static_cast<size_t>(entries) * sizeof(DataType),
                                             cudaMemcpyDeviceToHost, stream), "bf extract download");
    }
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "bf extract");
    const auto t2 = clock::now();
    for (char* p : owned) heap.free(p);
    prof.entries = entries;
    prof.plan = std::chrono::duration<double>(t1 - t0).count();
    prof.run = std::chrono::duration<double>(t2 - t1).count();
    return true;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
