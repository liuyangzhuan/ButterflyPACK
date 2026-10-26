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
// On one process extract() returns the whole blocks.  On several ranks
// extract_dist() returns this rank's rows of each block, as the host's
// BF_block_extraction does inside BF_ExtractElement, whose redistribution
// then sends them where they are wanted; see there.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "bf_device.hpp"
#include "bf_exchange.hpp"
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

// n elements from src + src_off to dst + dst_off
struct ExtractCopy {
    int64_t src = 0, dst = 0, n = 0;
};

void launch_extract_gather(const ExtractGather* items, int count, const int* idx, const DataType* values,
                           DataType* dst, cudaStream_t stream);
void launch_extract_copies(const ExtractCopy* items, int count, const DataType* src, DataType* dst,
                           cudaStream_t stream);
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
    // several ranks: the midpoint's elements sent and received, and its seconds
    int64_t sent = 0, received = 0;
    double mid = 0;
};

namespace detail {

// A run of a sorted list: the elements [start, start + count) are in group g.
struct Run {
    int g = 0, start = 0, count = 0;
};

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

// One stage's GEMMs, in the arrays MAGMA's variable-size batch reads.
struct ExtractStage {
    std::vector<magma_int_t> m, n, k, lda, ldb, ldc;
    std::vector<const DataType*> a, b;
    std::vector<DataType*> c;
    void resize(size_t count) {
        m.resize(count); n.resize(count); k.resize(count);
        lda.resize(count); ldb.resize(count); ldc.resize(count);
        a.resize(count); b.resize(count); c.resize(count);
    }
    size_t size() const { return m.size(); }
};

// BF(I_k, J_k) for k < count: rowidx[k] rows from `rows` and colidx[k]
// columns from `cols` (0-based, the butterfly's own order, stacked), into
// `out` (host), the blocks stacked in order, each column major.  Returns
// false when the device has no room for the partial products (nothing
// written; the caller runs the host path).
//
// The planning is per intersection and runs in parallel over them: every
// list below is flat, laid out by a counting pass and filled by a second.
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

    std::vector<int64_t> row_at(static_cast<size_t>(count) + 1, 0), col_at(row_at.size(), 0),
        out_at(row_at.size(), 0);
    for (int k = 0; k < count; ++k) {
        row_at[k + 1] = row_at[k] + rowidx[k];
        col_at[k + 1] = col_at[k] + colidx[k];
        out_at[k + 1] = out_at[k] + static_cast<int64_t>(rowidx[k]) * colidx[k];
    }
    const int64_t nrow = row_at[count], ncol = col_at[count], entries = out_at[count];

    // ---- each intersection's rows and columns in leaf order: srow[p'] is the
    // caller's row at leaf-order position p', pos_r its inverse
    std::vector<int> srow(static_cast<size_t>(nrow)), rleaf(srow.size()), pos_r(srow.size());
    std::vector<int> scol(static_cast<size_t>(ncol)), cleaf(scol.size()), pos_c(scol.size());
#pragma omp parallel
    {
        std::vector<int> p;
        auto order = [&](const int* idx, int64_t at, int n, const std::vector<int>& leaf_of, int* sorted, int* lf,
                         int* pos) {
            p.resize(static_cast<size_t>(n));
            std::iota(p.begin(), p.end(), 0);
            std::stable_sort(p.begin(), p.end(),
                             [&](int x, int y) { return leaf_of[idx[at + x]] < leaf_of[idx[at + y]]; });
            for (int q = 0; q < n; ++q) {
                sorted[at + q] = idx[at + p[q]];
                lf[at + q] = leaf_of[sorted[at + q]];
                pos[at + p[q]] = q;
            }
        };
#pragma omp for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            order(rows, row_at[k], rowidx[k], leaves.row_leaf, srow.data(), rleaf.data(), pos_r.data());
            order(cols, col_at[k], colidx[k], leaves.col_leaf, scol.data(), cleaf.data(), pos_c.data());
        }
    }

    // ---- the runs of every tree level t = 0 .. L: intersection k's row
    // groups (leaf >> (L - t)) are rr[t][ro[t][k] .. ro[t][k+1]), its column
    // groups (leaf >> t) cr[t][co[t][k] .. co[t][k+1])
    auto count_runs = [](const int* lf, int n, int shift) {
        int runs = 0;
        for (int q = 0; q < n; ++q)
            if (q == 0 || (lf[q] >> shift) != (lf[q - 1] >> shift)) ++runs;
        return runs;
    };
    std::vector<std::vector<Run>> rr(static_cast<size_t>(L) + 1), cr(static_cast<size_t>(L) + 1);
    std::vector<std::vector<int64_t>> ro(static_cast<size_t>(L) + 1, std::vector<int64_t>(row_at.size(), 0));
    std::vector<std::vector<int64_t>> co(ro);
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int t = 0; t <= L; ++t) {
            ro[t][k + 1] = count_runs(rleaf.data() + row_at[k], rowidx[k], L - t);
            co[t][k + 1] = count_runs(cleaf.data() + col_at[k], colidx[k], t);
        }
    }
    for (int t = 0; t <= L; ++t) {
        for (int k = 0; k < count; ++k) {
            ro[t][k + 1] += ro[t][k];
            co[t][k + 1] += co[t][k];
        }
        rr[t].resize(static_cast<size_t>(ro[t][count]));
        cr[t].resize(static_cast<size_t>(co[t][count]));
    }
    auto fill_runs = [](const int* lf, int n, int shift, Run* out) {
        int at = 0;
        for (int p = 0; p < n;) {
            const int g = lf[p] >> shift;
            int q = p + 1;
            while (q < n && (lf[q] >> shift) == g) ++q;
            out[at++] = Run{g, p, q - p};
            p = q;
        }
    };
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int t = 0; t <= L; ++t) {
            fill_runs(rleaf.data() + row_at[k], rowidx[k], L - t, rr[t].data() + ro[t][k]);
            fill_runs(cleaf.data() + col_at[k], colidx[k], t, cr[t].data() + co[t][k]);
        }
    }
    // the index of group g among an intersection's runs (sorted by group)
    auto find = [](const Run* r, int64_t n, int g) -> int {
        int64_t lo = 0, hi = n - 1;
        while (lo <= hi) {
            const int64_t mid = (lo + hi) / 2;
            if (r[mid].g == g) return static_cast<int>(mid);
            if (r[mid].g < g) lo = mid + 1;
            else hi = mid - 1;
        }
        throw std::logic_error("bf::gpu::extract: a parent group is missing");
    };

    // ---- the partial products of vector level t (groups of tree level t-1):
    // slot (k, fi, ci) is at[base[k] + fi * nC + ci], an offset into the
    // level's buffer.  Right half: rank x (columns in the run); left half:
    // (rows in the run) x rank.
    struct Level {
        std::vector<int64_t> base, at;
        int64_t total = 0;
    };
    auto layout = [&](int t, bool right, Level& s) {
        const auto& RR = rr[t - 1];
        const auto& CR = cr[t - 1];
        const auto& RO = ro[t - 1];
        const auto& CO = co[t - 1];
        s.base.assign(row_at.size(), 0);
        std::vector<int64_t> elems(row_at.size(), 0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            const int64_t nF = RO[k + 1] - RO[k], nC = CO[k + 1] - CO[k];
            int64_t e = 0;
            for (int64_t fi = 0; fi < nF; ++fi)
                for (int64_t ci = 0; ci < nC; ++ci)
                    e += static_cast<int64_t>(rank_of(t, RR[RO[k] + fi].g, CR[CO[k] + ci].g)) *
                         (right ? CR[CO[k] + ci].count : RR[RO[k] + fi].count);
            s.base[k + 1] = nF * nC;
            elems[k + 1] = e;
        }
        for (int k = 0; k < count; ++k) {
            s.base[k + 1] += s.base[k];
            elems[k + 1] += elems[k];
        }
        s.at.resize(static_cast<size_t>(s.base[count]));
        s.total = elems[count];
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            const int64_t nF = RO[k + 1] - RO[k], nC = CO[k + 1] - CO[k];
            int64_t off = elems[k], idx = s.base[k];
            for (int64_t fi = 0; fi < nF; ++fi) {
                for (int64_t ci = 0; ci < nC; ++ci) {
                    s.at[idx++] = off;
                    off += static_cast<int64_t>(rank_of(t, RR[RO[k] + fi].g, CR[CO[k] + ci].g)) *
                           (right ? CR[CO[k] + ci].count : RR[RO[k] + fi].count);
                }
            }
        }
    };
    std::vector<Level> R(static_cast<size_t>(h) + 2), W(static_cast<size_t>(L) + 2);
    for (int t = 1; t <= h + 1; ++t) layout(t, true, R[t]);
    for (int t = L + 1; t >= h + 1; --t) layout(t, false, W[t]);

    // ---- the device buffers: a level and half each, the leaf-order entries,
    // the caller-order entries
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    std::vector<std::pair<char*, int64_t>> owned;
    bool room = true;
    auto take = [&](int64_t elements) -> DataType* {
        if (!room) return nullptr;
        const int64_t n = std::max<int64_t>(elements, 1);
        char* p = heap.try_alloc(static_cast<size_t>(n) * sizeof(DataType), fmm::gpu::ResidentKind::Matvec);
        if (p == nullptr) {
            room = false;
            return nullptr;
        }
        owned.push_back({p, n});
        return reinterpret_cast<DataType*>(p);
    };
    std::vector<DataType*> dR(static_cast<size_t>(h) + 2, nullptr), dW(static_cast<size_t>(L) + 2, nullptr);
    for (int t = 1; t <= h + 1; ++t) dR[t] = take(R[t].total);
    for (int t = L + 1; t >= h + 1; --t) dW[t] = take(W[t].total);
    DataType* sorted = take(entries);
    DataType* d_out = take(entries);
    if (!room) {
        for (auto& o : owned) heap.free(o.first);
        return false;
    }

    // ---- the GEMM stages, in launch order: the right half's levels, the
    // left half's, the middle.  Each counts its products per intersection,
    // then fills them in parallel.
    std::vector<ExtractStage> stages;
    auto build = [&](auto items_of, auto fill_one) {
        std::vector<int64_t> at(row_at.size(), 0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) at[k + 1] = items_of(k);
        for (int k = 0; k < count; ++k) at[k + 1] += at[k];
        stages.emplace_back();
        ExtractStage& s = stages.back();
        s.resize(static_cast<size_t>(at[count]));
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) fill_one(k, s, at[k]);
    };
    auto put = [](ExtractStage& s, int64_t i, const DataType* a, int lda, const DataType* b, int ldb, DataType* c,
                  int ldc, int m, int n, int kk) {
        s.a[i] = a; s.b[i] = b; s.c[i] = c;
        s.m[i] = m; s.n[i] = n; s.k[i] = kk;
        s.lda[i] = std::max(lda, 1); s.ldb[i] = std::max(ldb, 1); s.ldc[i] = std::max(ldc, 1);
    };
    // right half: vector level t -> t+1 through Kerl(t), t = 1 .. h
    for (int t = 1; t <= h; ++t) {
        build([&](int k) { return (ro[t][k + 1] - ro[t][k]) * (co[t - 1][k + 1] - co[t - 1][k]); },
              [&](int k, ExtractStage& s, int64_t i) {
                  const Run* f_in = rr[t - 1].data() + ro[t - 1][k];
                  const int64_t nf_in = ro[t - 1][k + 1] - ro[t - 1][k];
                  const Run* f_out = rr[t].data() + ro[t][k];
                  const int64_t nf_out = ro[t][k + 1] - ro[t][k];
                  const Run* c_in = cr[t - 1].data() + co[t - 1][k];
                  const int64_t nc_in = co[t - 1][k + 1] - co[t - 1][k];
                  const Run* c_out = cr[t].data() + co[t][k];
                  const int64_t nc_out = co[t][k + 1] - co[t][k];
                  for (int64_t fo = 0; fo < nf_out; ++fo) {
                      const int F = f_out[fo].g;
                      const int fi = find(f_in, nf_in, F >> 1);
                      for (int64_t ci = 0; ci < nc_in; ++ci) {
                          const int C = c_in[ci].g;
                          const int cp = find(c_out, nc_out, C >> 1);
                          const BlockDesc& K = bf.block(t, F, C);
                          const int r_out = rank_of(t + 1, F, C >> 1);
                          const int r_in = rank_of(t, F >> 1, C);
                          const DataType* B = dR[t] + R[t].at[R[t].base[k] + fi * nc_in + ci];
                          DataType* Cp = dR[t + 1] + R[t + 1].at[R[t + 1].base[k] + fo * nc_out + cp] +
                                         static_cast<int64_t>(c_in[ci].start - c_out[cp].start) * r_out;
                          put(s, i++, bf.values + K.offset, K.lda, B, r_in, Cp, r_out, r_out, c_in[ci].count, r_in);
                      }
                  }
              });
    }
    // left half: vector level t -> t-1 through Kerl(t-1), t = L+1 .. h+2
    for (int t = L + 1; t >= h + 2; --t) {
        build([&](int k) { return (ro[t - 1][k + 1] - ro[t - 1][k]) * (co[t - 2][k + 1] - co[t - 2][k]); },
              [&](int k, ExtractStage& s, int64_t i) {
                  const Run* f_in = rr[t - 1].data() + ro[t - 1][k];
                  const int64_t nf_in = ro[t - 1][k + 1] - ro[t - 1][k];
                  const Run* f_out = rr[t - 2].data() + ro[t - 2][k];
                  const int64_t nf_out = ro[t - 2][k + 1] - ro[t - 2][k];
                  const Run* c_in = cr[t - 1].data() + co[t - 1][k];
                  const int64_t nc_in = co[t - 1][k + 1] - co[t - 1][k];
                  const Run* c_out = cr[t - 2].data() + co[t - 2][k];
                  const int64_t nc_out = co[t - 2][k + 1] - co[t - 2][k];
                  for (int64_t fi = 0; fi < nf_in; ++fi) {
                      const int F = f_in[fi].g;
                      const int fp = find(f_out, nf_out, F >> 1);
                      for (int64_t cj = 0; cj < nc_out; ++cj) {
                          const int J = c_out[cj].g;
                          const int ci = find(c_in, nc_in, J >> 1);
                          const BlockDesc& K = bf.block(t - 1, F, J);
                          const int r_in = rank_of(t, F, J >> 1);
                          const int r_out = rank_of(t - 1, F >> 1, J);
                          const DataType* A = dW[t] + W[t].at[W[t].base[k] + fi * nc_in + ci];
                          DataType* Cp = dW[t - 1] + W[t - 1].at[W[t - 1].base[k] + fp * nc_out + cj] +
                                         (f_in[fi].start - f_out[fp].start);
                          put(s, i++, A, f_in[fi].count, bf.values + K.offset, K.lda, Cp, f_out[fp].count,
                              f_in[fi].count, r_out, r_in);
                      }
                  }
              });
    }
    // the middle: entries in leaf order
    build([&](int k) { return (ro[h][k + 1] - ro[h][k]) * (co[h][k + 1] - co[h][k]); },
          [&](int k, ExtractStage& s, int64_t i) {
              const Run* fr = rr[h].data() + ro[h][k];
              const int64_t nf = ro[h][k + 1] - ro[h][k];
              const Run* cc = cr[h].data() + co[h][k];
              const int64_t nc = co[h][k + 1] - co[h][k];
              const int m = rowidx[k];
              for (int64_t fi = 0; fi < nf; ++fi) {
                  for (int64_t ci = 0; ci < nc; ++ci) {
                      const int r = rank_of(h + 1, fr[fi].g, cc[ci].g);
                      const DataType* A = dW[h + 1] + W[h + 1].at[W[h + 1].base[k] + fi * nc + ci];
                      const DataType* B = dR[h + 1] + R[h + 1].at[R[h + 1].base[k] + fi * nc + ci];
                      DataType* Cp = sorted + out_at[k] + fr[fi].start + static_cast<int64_t>(cc[ci].start) * m;
                      put(s, i++, A, fr[fi].count, B, r, Cp, m, fr[fi].count, cc[ci].count, r);
                  }
              }
          });

    // ---- gathers (V: one per column leaf of an intersection; U: one per row
    // leaf), their local indices, and the scatter's positions
    std::vector<int64_t> v_at(row_at.size(), 0), u_at(row_at.size(), 0);
    for (int k = 0; k < count; ++k) {
        v_at[k + 1] = v_at[k] + (co[0][k + 1] - co[0][k]);
        u_at[k + 1] = u_at[k] + (ro[L][k + 1] - ro[L][k]);
    }
    const int nv = static_cast<int>(v_at[count]), nu = static_cast<int>(u_at[count]);
    // the index array: V's local columns, U's local rows, pos_r, pos_c
    const int64_t idx_v = 0, idx_u = ncol, idx_pr = ncol + nrow, idx_pc = ncol + 2 * nrow;
    const int64_t nindex = 2 * (nrow + ncol);

    // ---- the metadata image: laid out first (it may move while it grows),
    // then filled
    constexpr int64_t kMaxBatch = 65535;
    struct Chunk {
        size_t stage = 0;
        int64_t first = 0, count = 0;
        size_t sizes = 0, pointers = 0;
        magma_int_t max_m = 0, max_n = 0, max_k = 0;
    };
    std::vector<Chunk> chunks;
    meta.clear();
    const size_t at_gather = meta.alloc(static_cast<size_t>(nv + nu) * sizeof(ExtractGather));
    const size_t at_scatter = meta.alloc(static_cast<size_t>(count) * sizeof(ExtractScatter));
    const size_t at_index = meta.alloc(static_cast<size_t>(nindex) * sizeof(int));
    for (size_t s = 0; s < stages.size(); ++s) {
        const int64_t n = static_cast<int64_t>(stages[s].size());
        for (int64_t first = 0; first < n; first += kMaxBatch) {
            Chunk c;
            c.stage = s;
            c.first = first;
            c.count = std::min(kMaxBatch, n - first);
            c.sizes = meta.alloc(static_cast<size_t>(6 * (c.count + 1)) * sizeof(magma_int_t));
            c.pointers = meta.alloc(static_cast<size_t>(3 * c.count) * sizeof(void*));
            chunks.push_back(c);
        }
    }
    auto* gathers = reinterpret_cast<ExtractGather*>(meta.host(at_gather));
    auto* scatters = reinterpret_cast<ExtractScatter*>(meta.host(at_scatter));
    int* index = reinterpret_cast<int*>(meta.host(at_index));
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int64_t q = 0; q < colidx[k]; ++q) {
            const int64_t p = col_at[k] + q;
            index[idx_v + p] = scol[p] - leaves.col_head[cleaf[p]];
            index[idx_pc + p] = pos_c[p];
        }
        for (int64_t q = 0; q < rowidx[k]; ++q) {
            const int64_t p = row_at[k] + q;
            index[idx_u + p] = srow[p] - leaves.row_head[rleaf[p]];
            index[idx_pr + p] = pos_r[p];
        }
        int64_t g = v_at[k];
        for (int64_t ci = co[0][k]; ci < co[0][k + 1]; ++ci) {
            const Run& run = cr[0][ci];
            const BlockDesc& v = bf.block(0, 0, run.g);
            ExtractGather e;
            e.src = v.offset; e.ld_src = v.lda; e.r = v.cols;
            e.count = run.count; e.transposed = 1;
            e.dst = R[1].at[R[1].base[k] + (ci - co[0][k])];
            e.idx = idx_v + col_at[k] + run.start;
            gathers[g++] = e;
        }
        g = nv + u_at[k];
        for (int64_t fi = ro[L][k]; fi < ro[L][k + 1]; ++fi) {
            const Run& run = rr[L][fi];
            const BlockDesc& u = bf.block(L + 1, run.g, 0);
            ExtractGather e;
            e.src = u.offset; e.ld_src = u.lda; e.r = u.cols;
            e.count = run.count; e.transposed = 0;
            e.dst = W[L + 1].at[W[L + 1].base[k] + (fi - ro[L][k])];
            e.idx = idx_u + row_at[k] + run.start;
            gathers[g++] = e;
        }
        ExtractScatter sc;
        sc.src = sc.dst = out_at[k];
        sc.m = rowidx[k];
        sc.n = colidx[k];
        sc.pos_r = idx_pr + row_at[k];
        sc.pos_c = idx_pc + col_at[k];
        scatters[k] = sc;
    }
    double flops = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : flops)
    for (size_t ch = 0; ch < chunks.size(); ++ch) {
        Chunk& c = chunks[ch];
        const ExtractStage& s = stages[c.stage];
        auto* sz = reinterpret_cast<magma_int_t*>(meta.host(c.sizes));
        auto* pt = reinterpret_cast<const void**>(meta.host(c.pointers));
        const int64_t n = c.count;
        for (int64_t q = 0; q < n; ++q) {
            const int64_t i = c.first + q;
            sz[0 * (n + 1) + q] = s.m[i];
            sz[1 * (n + 1) + q] = s.n[i];
            sz[2 * (n + 1) + q] = s.k[i];
            sz[3 * (n + 1) + q] = s.lda[i];
            sz[4 * (n + 1) + q] = s.ldb[i];
            sz[5 * (n + 1) + q] = s.ldc[i];
            pt[q] = s.a[i];
            pt[n + q] = s.b[i];
            pt[2 * n + q] = s.c[i];
            c.max_m = std::max(c.max_m, s.m[i]);
            c.max_n = std::max(c.max_n, s.n[i]);
            c.max_k = std::max(c.max_k, s.k[i]);
            flops += 2.0 * s.m[i] * s.n[i] * s.k[i];
        }
        for (int w = 0; w < 6; ++w) sz[w * (n + 1) + n] = 0;
    }
    prof.flops += flops;
    for (const ExtractStage& s : stages) prof.items += static_cast<int64_t>(s.size());
    const auto t1 = clock::now();

    // ---- launch, in dependence order on the backend's stream
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    magma_queue_t queue = fmm::gpu::Context::instance().queue();
    char* dm = meta.upload(meta_device, stream);
    // a product of rank zero writes nothing: everything starts at zero
    for (const auto& o : owned) {
        fmm::gpu::check_cuda(cudaMemsetAsync(o.first, 0, static_cast<size_t>(o.second) * sizeof(DataType), stream),
                             "bf extract zero");
    }
    if (nv > 0) launch_extract_gather(reinterpret_cast<const ExtractGather*>(dm + at_gather), nv,
                                      reinterpret_cast<const int*>(dm + at_index), bf.values, dR[1], stream);
    if (nu > 0) launch_extract_gather(reinterpret_cast<const ExtractGather*>(dm + at_gather) + nv, nu,
                                      reinterpret_cast<const int*>(dm + at_index), bf.values, dW[L + 1], stream);
    for (const Chunk& c : chunks) {
        if (c.max_m <= 0 || c.max_n <= 0) continue;
        auto* sz = reinterpret_cast<magma_int_t*>(dm + c.sizes);
        auto* pt = reinterpret_cast<void**>(dm + c.pointers);
        const int64_t n = c.count;
        fmm::gpu::gemm_vbatched<DataType>(
            MagmaNoTrans, MagmaNoTrans, sz, sz + (n + 1), sz + 2 * (n + 1), DataType(1),
            reinterpret_cast<const DataType* const*>(pt), sz + 3 * (n + 1),
            reinterpret_cast<const DataType* const*>(pt + n), sz + 4 * (n + 1), DataType(0),
            reinterpret_cast<DataType**>(pt + 2 * n), sz + 5 * (n + 1), static_cast<magma_int_t>(n), queue,
            c.max_m, c.max_n, std::max<magma_int_t>(c.max_k, 1));
    }
    if (count > 0) launch_extract_scatter(reinterpret_cast<const ExtractScatter*>(dm + at_scatter), count,
                                          reinterpret_cast<const int*>(dm + at_index), sorted, d_out, stream);
    if (entries > 0) {
        fmm::gpu::check_cuda(cudaMemcpyAsync(out, d_out, static_cast<size_t>(entries) * sizeof(DataType),
                                             cudaMemcpyDeviceToHost, stream), "bf extract download");
    }
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "bf extract");
    const auto t2 = clock::now();
    for (auto& o : owned) heap.free(o.first);
    prof.entries = entries;
    prof.plan = std::chrono::duration<double>(t1 - t0).count();
    prof.run = std::chrono::duration<double>(t2 - t1).count();
    return true;
}

// ---------------------------------------------------------------------------
// Several ranks.  The host's BF_block_extraction sweeps the requested entries
// through the butterfly as its matvec sweeps a vector: the right half (V,
// Kerl(1 .. h)) in the row-wise ordering, with a pair-wise broadcast after
// each level, the midpoint redistribution, the left half in the column-wise
// ordering with a pair-wise reduction after each level, and U; each rank ends
// with its rows of every intersection.  When every column group of tree level
// h has its leaves on one rank and so does every row group, the pair-wise
// exchanges carry nothing: a rank holds the whole right half of its column
// groups and the whole left half of its row groups.  That is the case here
// (2^h groups each way a rank at most, so up to 2^min(h, L-h) ranks): the
// right half's partial products are computed where their columns are, the
// left half's (from U, as extract() does) where their rows are, and the one
// exchange is the midpoint's - the right half's products of tree level h
// going to the rank of their row group, which multiplies them with its left
// half's.  Any other distribution declines, on every rank together, and the
// host runs.
//
// The request is the host's: every rank passes every intersection (checked),
// each row and column by its leaf and its position in the leaf, which of an
// intersection's rows this rank holds (in the order its rows of the block are
// wanted), and the rank holding each leaf.
struct ExtractRequest {
    int count = 0;
    const int* nr = nullptr;         // rows of each intersection
    const int* nc = nullptr;         // columns of each
    const int* row_leaf = nullptr;   // each row's leaf (zero based), the intersections' stacked
    const int* row_pos = nullptr;    // its position in the leaf
    const int* col_leaf = nullptr;   // each column's
    const int* col_pos = nullptr;
    const int* nloc = nullptr;       // this rank's rows of each intersection
    const int* local = nullptr;      // their positions in its rows, stacked
    const int* row_owner = nullptr;  // the rank (in the block's communicator) of each row leaf, 2^L
    const int* col_owner = nullptr;  // of each column leaf
};

// This rank's rows of every intersection into out (host), the blocks
// nloc[k] x nc[k], column major, stacked in order.  0: done; 1: declined on
// every rank (another distribution, no room on some rank).  Collective over
// the block's communicator.
inline int extract_dist(const BfDevice& bf, MPI_Comm comm, const ExtractRequest& q, DataType* out,
                        ExtractProfile& prof, fmm::gpu::MetaBuilder& meta, fmm::gpu::DeviceBuffer& meta_device) {
    using clock = std::chrono::steady_clock;
    using detail::Run;
    const auto t0 = clock::now();
    const int L = bf.level_butterfly;
    const int h = std::min(std::max(bf.level_half, 0), L);
    const int count = q.count;
    const int nleaf = 1 << L;
    int me = 0, nproc = 1;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &nproc);
    int ok = 1;

    std::vector<int64_t> row_at(static_cast<size_t>(count) + 1, 0), col_at(row_at.size(), 0),
        loc_at(row_at.size(), 0), out_at(row_at.size(), 0);
    for (int k = 0; k < count; ++k) {
        row_at[k + 1] = row_at[k] + q.nr[k];
        col_at[k + 1] = col_at[k] + q.nc[k];
        loc_at[k + 1] = loc_at[k] + q.nloc[k];
        out_at[k + 1] = out_at[k] + static_cast<int64_t>(q.nloc[k]) * q.nc[k];
    }
    const int64_t nrow = row_at[count], ncol = col_at[count], nlocal = loc_at[count], entries = out_at[count];

    // ---- the same request on every rank: a hash of it, its least and most
    {
        uint64_t hv = 1469598103934665603ull;
        auto mix = [&](const int* a, int64_t n) {
            for (int64_t i = 0; i < n; ++i) hv = (hv ^ static_cast<uint64_t>(static_cast<uint32_t>(a[i]))) * 1099511628211ull;
        };
        mix(&count, 1);
        mix(q.nr, count);
        mix(q.nc, count);
        mix(q.row_leaf, nrow);
        mix(q.row_pos, nrow);
        mix(q.col_leaf, ncol);
        mix(q.col_pos, ncol);
        mix(q.row_owner, nleaf);
        mix(q.col_owner, nleaf);
        unsigned long long lo = hv, hi = hv;
        MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, comm);
        MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm);
        if (lo != hi) ok = 0;
    }
    // ---- the rank of each row group and column group of tree level h, -1
    // where a group's leaves are on several ranks (or a leaf on none)
    std::vector<int> rgo(static_cast<size_t>(1) << h, -2), cgo(static_cast<size_t>(1) << (L - h), -2);
    for (int a = 0; a < nleaf; ++a) {
        int& g = rgo[static_cast<size_t>(a >> (L - h))];
        const int o = q.row_owner[a];
        g = (o < 0 || o >= nproc || (g != -2 && g != o)) ? -1 : o;
    }
    for (int b = 0; b < nleaf; ++b) {
        int& g = cgo[static_cast<size_t>(b >> h)];
        const int o = q.col_owner[b];
        g = (o < 0 || o >= nproc || (g != -2 && g != o)) ? -1 : o;
    }
    for (int g : rgo) if (g < 0) ok = 0;
    for (int g : cgo) if (g < 0) ok = 0;
    // a row group g of tree level t >= h, and a column group of level t <= h
    auto row_mine = [&](int t, int g) { return rgo[static_cast<size_t>(g >> (t - h))] == me; };
    auto col_mine = [&](int t, int g) { return cgo[static_cast<size_t>(g >> (h - t))] == me; };
    // a factor block this rank holds, by its zero-based global indices, or null
    auto block_at = [&](int l, int gi, int gj) -> const BlockDesc* {
        const LevelDesc& d = bf.level(l);
        const int di = gi + 1 - d.idx_r, dj = gj + 1 - d.idx_c;
        if (d.inc_r <= 0 || d.inc_c <= 0 || di < 0 || dj < 0 || di % d.inc_r != 0 || dj % d.inc_c != 0) return nullptr;
        const int i = di / d.inc_r, j = dj / d.inc_c;
        if (i >= d.nr || j >= d.nc) return nullptr;
        const BlockDesc& b = bf.blocks[static_cast<size_t>(d.first + static_cast<int64_t>(i) * d.nc + j)];
        return &b;
    };
    // the rank of vector level t's block (F, C), from the half that holds it:
    // the right half from the factor below it, the left half from the one
    // above; -1 when that factor is not here
    auto rank_r = [&](int t, int F, int C) -> int {
        const BlockDesc* b = t == 1 ? block_at(0, 0, C) : block_at(t - 1, F, 2 * C);
        return b == nullptr ? -1 : (t == 1 ? b->cols : b->rows);
    };
    auto rank_w = [&](int t, int F, int C) -> int {
        const BlockDesc* b = t == L + 1 ? block_at(L + 1, F, 0) : block_at(t, 2 * F, C);
        return b == nullptr ? -1 : b->cols;
    };

    // ---- each intersection's rows and columns in leaf order (stable), and
    // their runs at every tree level, as extract() does
    std::vector<int> srow(static_cast<size_t>(nrow)), rleaf(srow.size()), pos_r(srow.size());
    std::vector<int> scol(static_cast<size_t>(ncol)), cleaf(scol.size()), pos_c(scol.size());
#pragma omp parallel
    {
        std::vector<int> p;
        auto order = [&](const int* leaf, int64_t at, int n, int* sorted, int* lf, int* pos) {
            p.resize(static_cast<size_t>(n));
            std::iota(p.begin(), p.end(), 0);
            std::stable_sort(p.begin(), p.end(), [&](int x, int y) { return leaf[at + x] < leaf[at + y]; });
            for (int r = 0; r < n; ++r) {
                sorted[at + r] = p[r];  // the position in the caller's order
                lf[at + r] = leaf[at + p[r]];
                pos[at + p[r]] = r;
            }
        };
#pragma omp for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            order(q.row_leaf, row_at[k], q.nr[k], srow.data(), rleaf.data(), pos_r.data());
            order(q.col_leaf, col_at[k], q.nc[k], scol.data(), cleaf.data(), pos_c.data());
        }
    }
    auto count_runs = [](const int* lf, int n, int shift) {
        int runs = 0;
        for (int r = 0; r < n; ++r)
            if (r == 0 || (lf[r] >> shift) != (lf[r - 1] >> shift)) ++runs;
        return runs;
    };
    std::vector<std::vector<Run>> rr(static_cast<size_t>(L) + 1), cr(static_cast<size_t>(L) + 1);
    std::vector<std::vector<int64_t>> ro(static_cast<size_t>(L) + 1, std::vector<int64_t>(row_at.size(), 0));
    std::vector<std::vector<int64_t>> co(ro);
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int t = 0; t <= L; ++t) {
            ro[t][k + 1] = count_runs(rleaf.data() + row_at[k], q.nr[k], L - t);
            co[t][k + 1] = count_runs(cleaf.data() + col_at[k], q.nc[k], t);
        }
    }
    for (int t = 0; t <= L; ++t) {
        for (int k = 0; k < count; ++k) {
            ro[t][k + 1] += ro[t][k];
            co[t][k + 1] += co[t][k];
        }
        rr[t].resize(static_cast<size_t>(ro[t][count]));
        cr[t].resize(static_cast<size_t>(co[t][count]));
    }
    auto fill_runs = [](const int* lf, int n, int shift, Run* r_out) {
        int at = 0;
        for (int p = 0; p < n;) {
            const int g = lf[p] >> shift;
            int r = p + 1;
            while (r < n && (lf[r] >> shift) == g) ++r;
            r_out[at++] = Run{g, p, r - p};
            p = r;
        }
    };
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int t = 0; t <= L; ++t) {
            fill_runs(rleaf.data() + row_at[k], q.nr[k], L - t, rr[t].data() + ro[t][k]);
            fill_runs(cleaf.data() + col_at[k], q.nc[k], t, cr[t].data() + co[t][k]);
        }
    }
    auto find = [](const Run* r, int64_t n, int g) -> int {
        int64_t lo = 0, hi = n - 1;
        while (lo <= hi) {
            const int64_t mid = (lo + hi) / 2;
            if (r[mid].g == g) return static_cast<int>(mid);
            if (r[mid].g < g) lo = mid + 1;
            else hi = mid - 1;
        }
        throw std::logic_error("bf::gpu::extract_dist: a parent group is missing");
    };
    // this rank's rows of each intersection: a contiguous run of its leaf
    // order (its leaves are), starting at first[k]; pos_l maps each wanted row
    // to its place in that run
    std::vector<int64_t> first(row_at.size(), 0);
    std::vector<int> pos_l(static_cast<size_t>(std::max<int64_t>(nlocal, 1)));
    int bad_local = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : bad_local)
    for (int k = 0; k < count; ++k) {
        int64_t lo = q.nr[k];
        for (int r = 0; r < q.nr[k]; ++r)
            if (q.row_owner[rleaf[row_at[k] + r]] == me) { lo = r; break; }
        first[k] = lo;
        int seen = 0;
        for (int r = 0; r < q.nr[k]; ++r) seen += q.row_owner[rleaf[row_at[k] + r]] == me ? 1 : 0;
        if (seen != q.nloc[k]) ++bad_local;
        for (int r = 0; r < q.nloc[k]; ++r) {
            const int p = q.local[loc_at[k] + r];
            const int s = p >= 0 && p < q.nr[k] ? pos_r[row_at[k] + p] : -1;
            if (s < lo || s >= lo + q.nloc[k] || q.row_owner[q.row_leaf[row_at[k] + p]] != me) ++bad_local;
            pos_l[loc_at[k] + r] = static_cast<int>(s - lo);
        }
    }
    if (bad_local > 0) ok = 0;

    // ---- the partial products: slot (k, fi, ci) of a level at
    // at[base[k] + fi * nC + ci], sized by `size` (0: not this rank's)
    struct Level {
        std::vector<int64_t> base, at;
        int64_t total = 0;
    };
    int bad_rank = 0;
    auto layout = [&](int tl, Level& s, auto size) {  // tl: the tree level of the runs
        const auto& RR = rr[tl];
        const auto& CR = cr[tl];
        const auto& RO = ro[tl];
        const auto& CO = co[tl];
        s.base.assign(row_at.size(), 0);
        std::vector<int64_t> elems(row_at.size(), 0);
        for (int k = 0; k < count; ++k) {
            s.base[k + 1] = s.base[k] + (RO[k + 1] - RO[k]) * (CO[k + 1] - CO[k]);
        }
        s.at.assign(static_cast<size_t>(s.base[count]), -1);
        int bad = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : bad)
        for (int k = 0; k < count; ++k) {
            const int64_t nF = RO[k + 1] - RO[k], nC = CO[k + 1] - CO[k];
            int64_t e = 0;
            for (int64_t fi = 0; fi < nF; ++fi)
                for (int64_t ci = 0; ci < nC; ++ci) {
                    const int64_t z = size(RR[RO[k] + fi], CR[CO[k] + ci]);
                    if (z < 0) ++bad;
                    else e += z;
                }
            elems[k + 1] = e;
        }
        bad_rank += bad;
        for (int k = 0; k < count; ++k) elems[k + 1] += elems[k];
        s.total = elems[count];
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            const int64_t nF = RO[k + 1] - RO[k], nC = CO[k + 1] - CO[k];
            int64_t off = elems[k], idx = s.base[k];
            for (int64_t fi = 0; fi < nF; ++fi)
                for (int64_t ci = 0; ci < nC; ++ci, ++idx) {
                    const int64_t z = size(RR[RO[k] + fi], CR[CO[k] + ci]);
                    if (z > 0) {
                        s.at[static_cast<size_t>(idx)] = off;
                        off += z;
                    }
                }
        }
    };
    // right half, vector levels 1 .. h+1 (runs of tree level t-1): this rank's column groups
    std::vector<Level> R(static_cast<size_t>(h) + 2), W(static_cast<size_t>(L) + 2);
    Level Rmid;
    for (int t = 1; t <= h + 1; ++t) {
        layout(t - 1, R[t], [&](const Run& f, const Run& c) -> int64_t {
            if (!col_mine(t - 1, c.g)) return 0;
            const int r = rank_r(t, f.g, c.g);
            return r < 0 ? -1 : static_cast<int64_t>(r) * c.count;
        });
    }
    // left half, vector levels L+1 .. h+1: this rank's row groups
    for (int t = L + 1; t >= h + 1; --t) {
        layout(t - 1, W[t], [&](const Run& f, const Run& c) -> int64_t {
            if (!row_mine(t - 1, f.g)) return 0;
            const int r = rank_w(t, f.g, c.g);
            return r < 0 ? -1 : static_cast<int64_t>(r) * f.count;
        });
    }
    // the right half's products of vector level h+1 that meet this rank's rows
    layout(h, Rmid, [&](const Run& f, const Run& c) -> int64_t {
        if (!row_mine(h, f.g)) return 0;
        const int r = rank_w(h + 1, f.g, c.g);
        return r < 0 ? -1 : static_cast<int64_t>(r) * c.count;
    });
    if (bad_rank > 0) ok = 0;

    // ---- the GEMM stages, as extract()'s, over this rank's halves
    std::vector<ExtractStage> stages;
    int bad_block = 0;
    auto build = [&](auto items_of, auto fill_one) {
        std::vector<int64_t> at(row_at.size(), 0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) at[k + 1] = items_of(k);
        for (int k = 0; k < count; ++k) at[k + 1] += at[k];
        stages.emplace_back();
        ExtractStage& s = stages.back();
        s.resize(static_cast<size_t>(at[count]));
        int bad = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : bad)
        for (int k = 0; k < count; ++k) bad += fill_one(k, s, at[k]);
        bad_block += bad;
    };
    std::vector<DataType*> dR(static_cast<size_t>(h) + 2, nullptr), dW(static_cast<size_t>(L) + 2, nullptr);
    DataType* dRmid = nullptr;
    DataType* sorted = nullptr;
    DataType* d_out = nullptr;
    auto put = [](ExtractStage& s, int64_t i, const DataType* a, int lda, const DataType* b, int ldb, DataType* c,
                  int ldc, int m, int n, int kk) {
        s.a[i] = a; s.b[i] = b; s.c[i] = c;
        s.m[i] = m; s.n[i] = n; s.k[i] = kk;
        s.lda[i] = std::max(lda, 1); s.ldb[i] = std::max(ldb, 1); s.ldc[i] = std::max(ldc, 1);
    };
    // The buffers come first, since the items hold pointers into them; on no
    // room the plan goes on without them and the group declines below.
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    std::vector<std::pair<char*, int64_t>> owned;
    bool room = true;
    auto take = [&](int64_t elements) -> DataType* {
        const int64_t n = std::max<int64_t>(elements, 1);
        char* p = room ? heap.try_alloc(static_cast<size_t>(n) * sizeof(DataType), fmm::gpu::ResidentKind::Matvec)
                       : nullptr;
        if (p == nullptr) {
            room = false;
            return nullptr;
        }
        owned.push_back({p, n});
        return reinterpret_cast<DataType*>(p);
    };
    for (int t = 1; t <= h + 1; ++t) dR[t] = take(R[t].total);
    for (int t = L + 1; t >= h + 1; --t) dW[t] = take(W[t].total);
    dRmid = take(Rmid.total);
    sorted = take(entries);
    d_out = take(entries);
    if (!room) ok = 0;
    const bool go = ok != 0;  // pointers into missing buffers are never followed: nothing launches unless all is well

    if (go) {
        // right half: vector level t -> t+1 through Kerl(t), over this rank's column groups
        for (int t = 1; t <= h; ++t) {
            build([&](int k) {
                      int64_t n = 0;
                      for (int64_t ci = co[t - 1][k]; ci < co[t - 1][k + 1]; ++ci)
                          if (col_mine(t - 1, cr[t - 1][ci].g)) n += ro[t][k + 1] - ro[t][k];
                      return n;
                  },
                  [&](int k, ExtractStage& s, int64_t i) {
                      int bad = 0;
                      const Run* f_in = rr[t - 1].data() + ro[t - 1][k];
                      const int64_t nf_in = ro[t - 1][k + 1] - ro[t - 1][k];
                      const Run* f_out = rr[t].data() + ro[t][k];
                      const int64_t nf_out = ro[t][k + 1] - ro[t][k];
                      const Run* c_in = cr[t - 1].data() + co[t - 1][k];
                      const int64_t nc_in = co[t - 1][k + 1] - co[t - 1][k];
                      const Run* c_out = cr[t].data() + co[t][k];
                      const int64_t nc_out = co[t][k + 1] - co[t][k];
                      for (int64_t ci = 0; ci < nc_in; ++ci) {
                          const int C = c_in[ci].g;
                          if (!col_mine(t - 1, C)) continue;
                          const int cp = find(c_out, nc_out, C >> 1);
                          for (int64_t fo = 0; fo < nf_out; ++fo) {
                              const int F = f_out[fo].g;
                              const int fi = find(f_in, nf_in, F >> 1);
                              const BlockDesc* K = block_at(t, F, C);
                              const int r_out = rank_r(t + 1, F, C >> 1);
                              const int r_in = rank_r(t, F >> 1, C);
                              if (K == nullptr || r_out < 0 || r_in < 0) {
                                  ++bad;
                                  put(s, i++, nullptr, 1, nullptr, 1, nullptr, 1, 0, 0, 0);
                                  continue;
                              }
                              const DataType* B = dR[t] + R[t].at[R[t].base[k] + fi * nc_in + ci];
                              DataType* Cp = dR[t + 1] + R[t + 1].at[R[t + 1].base[k] + fo * nc_out + cp] +
                                             static_cast<int64_t>(c_in[ci].start - c_out[cp].start) * r_out;
                              put(s, i++, bf.values + K->offset, K->lda, B, r_in, Cp, r_out, r_out, c_in[ci].count, r_in);
                          }
                      }
                      return bad;
                  });
        }
        // left half: vector level t -> t-1 through Kerl(t-1), over this rank's row groups
        for (int t = L + 1; t >= h + 2; --t) {
            build([&](int k) {
                      int64_t n = 0;
                      for (int64_t fi = ro[t - 1][k]; fi < ro[t - 1][k + 1]; ++fi)
                          if (row_mine(t - 1, rr[t - 1][fi].g)) n += co[t - 2][k + 1] - co[t - 2][k];
                      return n;
                  },
                  [&](int k, ExtractStage& s, int64_t i) {
                      int bad = 0;
                      const Run* f_in = rr[t - 1].data() + ro[t - 1][k];
                      const int64_t nf_in = ro[t - 1][k + 1] - ro[t - 1][k];
                      const Run* f_out = rr[t - 2].data() + ro[t - 2][k];
                      const int64_t nf_out = ro[t - 2][k + 1] - ro[t - 2][k];
                      const Run* c_in = cr[t - 1].data() + co[t - 1][k];
                      const int64_t nc_in = co[t - 1][k + 1] - co[t - 1][k];
                      const Run* c_out = cr[t - 2].data() + co[t - 2][k];
                      const int64_t nc_out = co[t - 2][k + 1] - co[t - 2][k];
                      for (int64_t fi = 0; fi < nf_in; ++fi) {
                          const int F = f_in[fi].g;
                          if (!row_mine(t - 1, F)) continue;
                          const int fp = find(f_out, nf_out, F >> 1);
                          for (int64_t cj = 0; cj < nc_out; ++cj) {
                              const int J = c_out[cj].g;
                              const int ci = find(c_in, nc_in, J >> 1);
                              const BlockDesc* K = block_at(t - 1, F, J);
                              const int r_in = rank_w(t, F, J >> 1);
                              const int r_out = rank_w(t - 1, F >> 1, J);
                              if (K == nullptr || r_in < 0 || r_out < 0) {
                                  ++bad;
                                  put(s, i++, nullptr, 1, nullptr, 1, nullptr, 1, 0, 0, 0);
                                  continue;
                              }
                              const DataType* A = dW[t] + W[t].at[W[t].base[k] + fi * nc_in + ci];
                              DataType* Cp = dW[t - 1] + W[t - 1].at[W[t - 1].base[k] + fp * nc_out + cj] +
                                             (f_in[fi].start - f_out[fp].start);
                              put(s, i++, A, f_in[fi].count, bf.values + K->offset, K->lda, Cp, f_out[fp].count,
                                  f_in[fi].count, r_out, r_in);
                          }
                      }
                      return bad;
                  });
        }
        // the middle: this rank's rows, every column, in leaf order
        build([&](int k) {
                  int64_t n = 0;
                  for (int64_t fi = ro[h][k]; fi < ro[h][k + 1]; ++fi)
                      if (row_mine(h, rr[h][fi].g)) n += co[h][k + 1] - co[h][k];
                  return n;
              },
              [&](int k, ExtractStage& s, int64_t i) {
                  const Run* fr = rr[h].data() + ro[h][k];
                  const int64_t nf = ro[h][k + 1] - ro[h][k];
                  const Run* cc = cr[h].data() + co[h][k];
                  const int64_t nc = co[h][k + 1] - co[h][k];
                  const int m = q.nloc[k];
                  for (int64_t fi = 0; fi < nf; ++fi) {
                      if (!row_mine(h, fr[fi].g)) continue;
                      for (int64_t ci = 0; ci < nc; ++ci) {
                          const int r = rank_w(h + 1, fr[fi].g, cc[ci].g);
                          const DataType* A = dW[h + 1] + W[h + 1].at[W[h + 1].base[k] + fi * nc + ci];
                          const DataType* B = dRmid + Rmid.at[Rmid.base[k] + fi * nc + ci];
                          DataType* Cp = sorted + out_at[k] + (fr[fi].start - first[k]) +
                                         static_cast<int64_t>(cc[ci].start) * m;
                          put(s, i++, A, fr[fi].count, B, r, Cp, m, fr[fi].count, cc[ci].count, r);
                      }
                  }
                  return 0;
              });
    }
    if (bad_block > 0) ok = 0;

    // ---- the midpoint: the right half's products of vector level h+1 to
    // the rank of their row group, in one order on both ends (intersection,
    // row run, column run); each end sizes a product from its own factor
    std::vector<ExtractCopy> pack, unpack;
    std::vector<int64_t> scount(static_cast<size_t>(nproc), 0), rcount(scount.size(), 0);
    if (go) {
        std::vector<std::vector<ExtractCopy>> sends(static_cast<size_t>(nproc)), recvs(static_cast<size_t>(nproc));
        for (int k = 0; k < count; ++k) {
            const Run* fr = rr[h].data() + ro[h][k];
            const int64_t nf = ro[h][k + 1] - ro[h][k];
            const Run* cc = cr[h].data() + co[h][k];
            const int64_t nc = co[h][k + 1] - co[h][k];
            for (int64_t fi = 0; fi < nf; ++fi) {
                const int dst = rgo[static_cast<size_t>(fr[fi].g)];
                for (int64_t ci = 0; ci < nc; ++ci) {
                    const int src = cgo[static_cast<size_t>(cc[ci].g)];
                    const int64_t slot = fi * nc + ci;
                    if (src == me) {
                        const int64_t n = static_cast<int64_t>(rank_r(h + 1, fr[fi].g, cc[ci].g)) * cc[ci].count;
                        sends[static_cast<size_t>(dst)].push_back(ExtractCopy{R[h + 1].at[R[h + 1].base[k] + slot], 0, n});
                    }
                    if (dst == me) {
                        const int64_t n = static_cast<int64_t>(rank_w(h + 1, fr[fi].g, cc[ci].g)) * cc[ci].count;
                        recvs[static_cast<size_t>(src)].push_back(ExtractCopy{0, Rmid.at[Rmid.base[k] + slot], n});
                    }
                }
            }
        }
        int64_t sat = 0, rat = 0;
        for (int p = 0; p < nproc; ++p) {
            for (ExtractCopy& c : sends[static_cast<size_t>(p)]) {
                c.dst = sat;
                sat += c.n;
                scount[static_cast<size_t>(p)] += c.n;
                if (c.n > 0) pack.push_back(c);
            }
            for (ExtractCopy& c : recvs[static_cast<size_t>(p)]) {
                c.src = rat;
                rat += c.n;
                rcount[static_cast<size_t>(p)] += c.n;
                if (c.n > 0) unpack.push_back(c);
            }
        }
    }
    // the midpoint's buffers, with the others: a rank without room declines
    // with the group
    int64_t stotal = 0, rtotal = 0;
    for (int p = 0; p < nproc; ++p) {
        stotal += scount[static_cast<size_t>(p)];
        rtotal += rcount[static_cast<size_t>(p)];
    }
    DataType* dsend = take(stotal);
    DataType* drecv = take(rtotal);
    if (!room) ok = 0;
    // both ends must agree on what moves (they size it from different factors)
    {
        std::vector<int64_t> told(scount.size(), 0);
        MPI_Alltoall(scount.data(), 1, MPI_INT64_T, told.data(), 1, MPI_INT64_T, comm);
        for (int p = 0; p < nproc; ++p) {
            if (told[static_cast<size_t>(p)] != rcount[static_cast<size_t>(p)] && go) ok = 0;
            if (scount[static_cast<size_t>(p)] > INT32_MAX || rcount[static_cast<size_t>(p)] > INT32_MAX) ok = 0;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, comm);
    if (!ok) {
        for (auto& o : owned) heap.free(o.first);
        return 1;
    }

    // ---- gathers (V: this rank's column leaves; U: its row leaves), the
    // index array, the scatter of this rank's rows
    std::vector<int64_t> v_at(row_at.size(), 0), u_at(row_at.size(), 0);
    for (int k = 0; k < count; ++k) {
        int64_t nv = 0, nu = 0;
        for (int64_t ci = co[0][k]; ci < co[0][k + 1]; ++ci) nv += col_mine(0, cr[0][ci].g) ? 1 : 0;
        for (int64_t fi = ro[L][k]; fi < ro[L][k + 1]; ++fi) nu += row_mine(L, rr[L][fi].g) ? 1 : 0;
        v_at[k + 1] = v_at[k] + nv;
        u_at[k + 1] = u_at[k] + nu;
    }
    const int nv = static_cast<int>(v_at[count]), nu = static_cast<int>(u_at[count]);
    // the index array: V's rows to take (by column), U's (by row), pos_l, pos_c
    const int64_t idx_v = 0, idx_u = ncol, idx_pr = ncol + nrow, idx_pc = ncol + nrow + nlocal;
    const int64_t nindex = ncol + nrow + nlocal + ncol;
    constexpr int64_t kMaxBatch = 65535;
    struct Chunk {
        size_t stage = 0;
        int64_t first = 0, count = 0;
        size_t sizes = 0, pointers = 0;
        magma_int_t max_m = 0, max_n = 0, max_k = 0;
    };
    std::vector<Chunk> chunks;
    meta.clear();
    const size_t at_gather = meta.alloc(static_cast<size_t>(nv + nu) * sizeof(ExtractGather));
    const size_t at_scatter = meta.alloc(static_cast<size_t>(count) * sizeof(ExtractScatter));
    const size_t at_index = meta.alloc(static_cast<size_t>(std::max<int64_t>(nindex, 1)) * sizeof(int));
    const size_t at_pack = meta.alloc(pack.size() * sizeof(ExtractCopy));
    const size_t at_unpack = meta.alloc(unpack.size() * sizeof(ExtractCopy));
    for (size_t s = 0; s < stages.size(); ++s) {
        const int64_t n = static_cast<int64_t>(stages[s].size());
        for (int64_t f = 0; f < n; f += kMaxBatch) {
            Chunk c;
            c.stage = s;
            c.first = f;
            c.count = std::min(kMaxBatch, n - f);
            c.sizes = meta.alloc(static_cast<size_t>(6 * (c.count + 1)) * sizeof(magma_int_t));
            c.pointers = meta.alloc(static_cast<size_t>(3 * c.count) * sizeof(void*));
            chunks.push_back(c);
        }
    }
    auto* gathers = reinterpret_cast<ExtractGather*>(meta.host(at_gather));
    auto* scatters = reinterpret_cast<ExtractScatter*>(meta.host(at_scatter));
    int* index = reinterpret_cast<int*>(meta.host(at_index));
    if (!pack.empty()) std::memcpy(meta.host(at_pack), pack.data(), pack.size() * sizeof(ExtractCopy));
    if (!unpack.empty()) std::memcpy(meta.host(at_unpack), unpack.data(), unpack.size() * sizeof(ExtractCopy));
#pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < count; ++k) {
        for (int64_t c = 0; c < q.nc[k]; ++c) {
            const int64_t p = col_at[k] + c;
            index[idx_v + p] = q.col_pos[col_at[k] + scol[p]];
            index[idx_pc + p] = pos_c[p];
        }
        for (int64_t r = 0; r < q.nr[k]; ++r) {
            const int64_t p = row_at[k] + r;
            index[idx_u + p] = q.row_pos[row_at[k] + srow[p]];
        }
        for (int64_t r = 0; r < q.nloc[k]; ++r) index[idx_pr + loc_at[k] + r] = pos_l[loc_at[k] + r];
        int64_t g = v_at[k];
        for (int64_t ci = co[0][k]; ci < co[0][k + 1]; ++ci) {
            const Run& run = cr[0][ci];
            if (!col_mine(0, run.g)) continue;
            const BlockDesc* v = block_at(0, 0, run.g);
            ExtractGather e;
            e.src = v->offset; e.ld_src = v->lda; e.r = v->cols;
            e.count = run.count; e.transposed = 1;
            e.dst = R[1].at[R[1].base[k] + (ci - co[0][k])];
            e.idx = idx_v + col_at[k] + run.start;
            gathers[g++] = e;
        }
        g = nv + u_at[k];
        for (int64_t fi = ro[L][k]; fi < ro[L][k + 1]; ++fi) {
            const Run& run = rr[L][fi];
            if (!row_mine(L, run.g)) continue;
            const BlockDesc* u = block_at(L + 1, run.g, 0);
            ExtractGather e;
            e.src = u->offset; e.ld_src = u->lda; e.r = u->cols;
            e.count = run.count; e.transposed = 0;
            e.dst = W[L + 1].at[W[L + 1].base[k] + (fi - ro[L][k])];
            e.idx = idx_u + row_at[k] + run.start;
            gathers[g++] = e;
        }
        ExtractScatter sc;
        sc.src = sc.dst = out_at[k];
        sc.m = q.nloc[k];
        sc.n = q.nc[k];
        sc.pos_r = idx_pr + loc_at[k];
        sc.pos_c = idx_pc + col_at[k];
        scatters[k] = sc;
    }
    double flops = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : flops)
    for (size_t ch = 0; ch < chunks.size(); ++ch) {
        Chunk& c = chunks[ch];
        const ExtractStage& s = stages[c.stage];
        auto* sz = reinterpret_cast<magma_int_t*>(meta.host(c.sizes));
        auto* pt = reinterpret_cast<const void**>(meta.host(c.pointers));
        const int64_t n = c.count;
        for (int64_t r = 0; r < n; ++r) {
            const int64_t i = c.first + r;
            sz[0 * (n + 1) + r] = s.m[i];
            sz[1 * (n + 1) + r] = s.n[i];
            sz[2 * (n + 1) + r] = s.k[i];
            sz[3 * (n + 1) + r] = s.lda[i];
            sz[4 * (n + 1) + r] = s.ldb[i];
            sz[5 * (n + 1) + r] = s.ldc[i];
            pt[r] = s.a[i];
            pt[n + r] = s.b[i];
            pt[2 * n + r] = s.c[i];
            c.max_m = std::max(c.max_m, s.m[i]);
            c.max_n = std::max(c.max_n, s.n[i]);
            c.max_k = std::max(c.max_k, s.k[i]);
            flops += 2.0 * s.m[i] * s.n[i] * s.k[i];
        }
        for (int w = 0; w < 6; ++w) sz[w * (n + 1) + n] = 0;
    }
    prof.flops += flops;
    for (const ExtractStage& s : stages) prof.items += static_cast<int64_t>(s.size());
    const auto t1 = clock::now();

    // ---- launch: the halves, the midpoint, the middle, the scatter
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    magma_queue_t queue = fmm::gpu::Context::instance().queue();
    char* dm = meta.upload(meta_device, stream);
    for (const auto& o : owned) {
        fmm::gpu::check_cuda(cudaMemsetAsync(o.first, 0, static_cast<size_t>(o.second) * sizeof(DataType), stream),
                             "bf extract zero");
    }
    if (nv > 0) launch_extract_gather(reinterpret_cast<const ExtractGather*>(dm + at_gather), nv,
                                      reinterpret_cast<const int*>(dm + at_index), bf.values, dR[1], stream);
    if (nu > 0) launch_extract_gather(reinterpret_cast<const ExtractGather*>(dm + at_gather) + nv, nu,
                                      reinterpret_cast<const int*>(dm + at_index), bf.values, dW[L + 1], stream);
    auto run_stage = [&](size_t s) {
        for (const Chunk& c : chunks) {
            if (c.stage != s || c.max_m <= 0 || c.max_n <= 0) continue;
            auto* sz = reinterpret_cast<magma_int_t*>(dm + c.sizes);
            auto* pt = reinterpret_cast<void**>(dm + c.pointers);
            const int64_t n = c.count;
            fmm::gpu::gemm_vbatched<DataType>(
                MagmaNoTrans, MagmaNoTrans, sz, sz + (n + 1), sz + 2 * (n + 1), DataType(1),
                reinterpret_cast<const DataType* const*>(pt), sz + 3 * (n + 1),
                reinterpret_cast<const DataType* const*>(pt + n), sz + 4 * (n + 1), DataType(0),
                reinterpret_cast<DataType**>(pt + 2 * n), sz + 5 * (n + 1), static_cast<magma_int_t>(n), queue,
                c.max_m, c.max_n, std::max<magma_int_t>(c.max_k, 1));
        }
    };
    // the right half's and the left half's stages, then (last) the middle's
    for (size_t s = 0; s + 1 < stages.size(); ++s) run_stage(s);
    // the midpoint: packed, exchanged between the devices (alltoallv_device:
    // NCCL on one node), unpacked
    const auto tm = clock::now();
    launch_extract_copies(reinterpret_cast<const ExtractCopy*>(dm + at_pack), static_cast<int>(pack.size()), dR[h + 1],
                          dsend, stream);
    alltoallv_device(comm, dsend, scount.data(), drecv, rcount.data(), stream);
    launch_extract_copies(reinterpret_cast<const ExtractCopy*>(dm + at_unpack), static_cast<int>(unpack.size()), drecv,
                          dRmid, stream);
    prof.mid = std::chrono::duration<double>(clock::now() - tm).count();
    prof.sent = stotal;
    prof.received = rtotal;
    if (!stages.empty()) run_stage(stages.size() - 1);
    if (count > 0) launch_extract_scatter(reinterpret_cast<const ExtractScatter*>(dm + at_scatter), count,
                                          reinterpret_cast<const int*>(dm + at_index), sorted, d_out, stream);
    if (entries > 0) {
        fmm::gpu::check_cuda(cudaMemcpyAsync(out, d_out, static_cast<size_t>(entries) * sizeof(DataType),
                                             cudaMemcpyDeviceToHost, stream), "bf extract download");
    }
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "bf extract");
    const auto t2 = clock::now();
    for (auto& o : owned) heap.free(o.first);
    prof.entries = entries;
    prof.plan = std::chrono::duration<double>(t1 - t0).count();
    prof.run = std::chrono::duration<double>(t2 - t1).count();
    return 0;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
