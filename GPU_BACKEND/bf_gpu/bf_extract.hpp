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
#include <cstring>
#include <memory>
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

// n elements from src + src_off to dst + dst_off: contiguous, or with
// rows > 0 a rows x (n / rows) block into a matrix of leading dimension ld
struct ExtractCopy {
    int64_t src = 0, dst = 0, n = 0;
    int64_t ld = 0;
    int rows = 0;
};

// (each defined for double and fmm::gpu::dcomplex in bf_extract.cu)
template<typename T>
void launch_extract_gather(const ExtractGather* items, int count, const int* idx, const T* values, T* dst,
                           cudaStream_t stream);
template<typename T>
void launch_extract_copies(const ExtractCopy* items, int count, const T* src, T* dst, cudaStream_t stream);
template<typename T>
void launch_extract_scatter(const ExtractScatter* items, int count, const int* idx, const T* src, T* dst,
                            cudaStream_t stream);

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
    // several ranks: the exchanges (the midpoint's and any other), their
    // elements sent and received, and the midpoint's seconds
    int exchanges = 0;
    int64_t sent = 0, received = 0;
    double mid = 0;
};

namespace detail {

// A run of a sorted list: the elements [start, start + count) are in group g.
struct Run {
    int g = 0, start = 0, count = 0;
};

}  // namespace detail

template<typename T>
inline void prepare_leaves(const BfDevice<T>& bf, ExtractLeaves& lv) {
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
template<typename T>
struct ExtractStage {
    using DataType = T;
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
template<typename T>
inline bool extract(const BfDevice<T>& bf, ExtractLeaves& leaves, int count, const int* rowidx, const int* colidx,
                    const int* rows, const int* cols, T* out, ExtractProfile& prof,
                    fmm::gpu::MetaBuilder& meta, fmm::gpu::DeviceBuffer& meta_device) {
    using DataType = T;
    using ExtractStage = bf::gpu::ExtractStage<T>;
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
// with its rows of every intersection.  Here the halves meet in the middle as
// in extract(): each of the right half's partial products is made on the
// rank holding the factor block that makes it, from V up, and each of the
// left half's likewise from U down; the right half's products of vector
// level h+1 go to the rank holding the left half's product each meets, which
// multiplies the two.
//
// Who holds a block follows from the leaves' owners, since GetBlockPID maps
// a block to a leaf and the leaf to the head of its process group: in the
// row-wise ordering, block (F, G) of level t - row F, column group G (the
// kernel's columns 2G and 2G+1) - to column leaf (G << t) | rev_t(F), rev_t
// reversing F's t bits; in the column-wise ordering, block (G, J) of level
// t - row group G (the kernel's rows 2G and 2G+1), column J - to row leaf
// (G << (L-t+1)) | rev_{L-t+1}(J).  The right half's product (F, C) of
// vector level t is used by level t's blocks (2F, C>>1) and (2F+1, C>>1),
// whose leaves are its maker's and the one differing from that in bit t-1;
// the left half's product (F, C) of vector level t by level t-1's (F>>1, 2C)
// and (F>>1, 2C+1), its maker's leaf and the one differing in bit L-t+1.  So
// a level's products all stay with their makers unless two leaves one bit
// apart there have different owners: with 2^p ranks of 2^(L-p) consecutive
// leaves each, the right half's levels t > L-p and the left half's
// t <= p+1, none of them when p <= min(h, L-h), the common case.  Where they
// do, each product goes to its other user once made - the host's pair-wise
// exchanges.  And when a row group of tree level h has several owners
// (p > L-h), the rank meeting its halves sends each owner its rows.
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

namespace detail {

// x's low `bits` bits in reverse order
inline int reverse_bits(int x, int bits) {
    if (bits <= 0) return 0;
    uint32_t v = static_cast<uint32_t>(x);
    v = ((v >> 1) & 0x55555555u) | ((v & 0x55555555u) << 1);
    v = ((v >> 2) & 0x33333333u) | ((v & 0x33333333u) << 2);
    v = ((v >> 4) & 0x0F0F0F0Fu) | ((v & 0x0F0F0F0Fu) << 4);
    v = ((v >> 8) & 0x00FF00FFu) | ((v & 0x00FF00FFu) << 8);
    v = (v >> 16) | (v << 16);
    return static_cast<int>(v >> (32 - bits));
}

// The lowest bit at which two leaves differing only there can have different
// owners (L when one rank owns all 2^L leaves); -1 when an owner is not a
// rank of the group or a rank's leaves are not consecutive.  With consecutive
// leaves it is the least number of trailing zeros of a leaf that starts a
// rank's leaves: such a leaf and the one below it in that bit have different
// owners, and a split with more trailing zeros cuts no pair at a lower bit.
inline int owner_split_bit(const int* owner, int L, int nproc) {
    int low = L;
    std::vector<char> seen(static_cast<size_t>(nproc), 0);
    for (int a = 0; a < (1 << L); ++a) {
        const int o = owner[a];
        if (o < 0 || o >= nproc) return -1;
        if (a > 0 && o != owner[a - 1]) {
            if (seen[static_cast<size_t>(o)]) return -1;
            low = std::min(low, __builtin_ctz(static_cast<unsigned>(a)));
        }
        seen[static_cast<size_t>(o)] = 1;
    }
    return low;
}

// The runs [lo, hi) of a list sorted by group whose groups are parent's
// children (2 parent and 2 parent + 1).
inline void child_runs(const Run* r, int64_t n, int parent, int64_t& lo, int64_t& hi) {
    lo = std::lower_bound(r, r + n, 2 * parent, [](const Run& x, int g) { return x.g < g; }) - r;
    hi = lo;
    while (hi < n && (r[hi].g >> 1) == parent) ++hi;
}

}  // namespace detail

// This rank's rows of every intersection into out (host), the blocks
// nloc[k] x nc[k], column major, stacked in order.  0: done; 1: declined on
// every rank (leaves not split among the ranks in consecutive runs, no room
// on some rank).  Collective over the block's communicator.
template<typename T>
inline int extract_dist(const BfDevice<T>& bf, MPI_Comm comm, const ExtractRequest& q, T* out, ExtractProfile& prof,
                        fmm::gpu::MetaBuilder& meta, fmm::gpu::DeviceBuffer& meta_device) {
    using DataType = T;
    using ExtractStage = bf::gpu::ExtractStage<T>;
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
        if (lo != hi) return 1;  // every rank sees it
    }
    // ---- the leaves' owners, the same on every rank now, so every rank
    // declines here together or none does
    const int low_row = detail::owner_split_bit(q.row_owner, L, nproc);
    const int low_col = detail::owner_split_bit(q.col_owner, L, nproc);
    if (low_row < 0 || low_col < 0) return 1;
    // where a product has a second user (see above): the right half's of
    // vector level t <= h, the left half's of level t >= h+2, and the rows
    // of a row group of tree level h
    auto shared_r = [&](int t) { return t - 1 >= low_col; };
    auto shared_w = [&](int t) { return L - t + 1 >= low_row; };
    const bool shared_out = low_row < L - h;
    // the maker of the right half's product (F, C) of vector level t (block
    // (F, C) of row-wise level t-1) and of the left half's (block (F, C) of
    // column-wise level t); where no split falls among the leaves F (C)
    // picks from, any of them names the rank
    auto own_r = [&](int t, int F, int C) {
        const int b = t - 1;
        return q.col_owner[(C << b) | (b > low_col ? detail::reverse_bits(F, b) : 0)];
    };
    auto own_w = [&](int t, int F, int C) {
        const int b = L - t + 1;
        return q.row_owner[(F << b) | (b > low_row ? detail::reverse_bits(C, b) : 0)];
    };
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
    // the rank of vector level t's block (F, C), from the factor that makes
    // it: the right half's from the one below, the left half's from the one
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
    // the owners' parts of rows [p0, p1) of intersection k's leaf order:
    // part(r0, r1, owner) for each, in order
    auto by_owner = [&](int k, int64_t p0, int64_t p1, auto part) {
        for (int64_t r = p0; r < p1;) {
            const int o = q.row_owner[rleaf[row_at[k] + r]];
            int64_t e = r + 1;
            while (e < p1 && q.row_owner[rleaf[row_at[k] + e]] == o) ++e;
            part(r, e, o);
            r = e;
        }
    };

    // ---- the partial products: slot (k, fi, ci) of a level at
    // at[base[k] + fi * nC + ci] (-1: not on this rank), size[...] elements.
    // A level of a large request has millions of slots, most of them other
    // ranks', so the two arrays are written only by the layout's parallel
    // passes, never filled first (a serial fill of one cost 3-6 ms a call).
    struct Level {
        std::vector<int64_t> base;
        std::unique_ptr<int64_t[]> at, size;
        int64_t total = 0;
    };
    int bad_rank = 0;
    auto layout = [&](int tl, Level& s, auto size_of) {  // tl: the tree level of the runs
        const auto& RR = rr[tl];
        const auto& CR = cr[tl];
        const auto& RO = ro[tl];
        const auto& CO = co[tl];
        s.base.assign(row_at.size(), 0);
        for (int k = 0; k < count; ++k) s.base[k + 1] = s.base[k] + (RO[k + 1] - RO[k]) * (CO[k + 1] - CO[k]);
        const size_t slots = static_cast<size_t>(std::max<int64_t>(s.base[count], 1));
        s.at.reset(new int64_t[slots]);
        s.size.reset(new int64_t[slots]);
        std::vector<int64_t> elems(row_at.size(), 0);
        int bad = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : bad)
        for (int k = 0; k < count; ++k) {
            const int64_t nF = RO[k + 1] - RO[k], nC = CO[k + 1] - CO[k];
            int64_t e = 0, idx = s.base[k];
            for (int64_t fi = 0; fi < nF; ++fi)
                for (int64_t ci = 0; ci < nC; ++ci, ++idx) {
                    int64_t z = size_of(k, RR[RO[k] + fi], CR[CO[k] + ci]);
                    if (z < 0) {
                        ++bad;
                        z = 0;
                    }
                    s.size[idx] = z;
                    e += z;
                }
            elems[k + 1] = e;
        }
        bad_rank += bad;
        for (int k = 0; k < count; ++k) elems[k + 1] += elems[k];
        s.total = elems[count];
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            int64_t off = elems[k];
            for (int64_t idx = s.base[k]; idx < s.base[k + 1]; ++idx) {
                const int64_t z = s.size[idx];
                s.at[idx] = z > 0 ? off : -1;
                off += z;
            }
        }
    };
    // right half, vector levels 1 .. h+1 (runs of tree level t-1): the
    // products this rank makes, and at a shared level those it uses, sized
    // from the factor that uses them
    std::vector<Level> R(static_cast<size_t>(h) + 2), W(static_cast<size_t>(L) + 2);
    Level Rmid;
    for (int t = 1; t <= h + 1; ++t) {
        const bool used_here = t <= h && shared_r(t);
        layout(t - 1, R[t], [&](int k, const Run& f, const Run& c) -> int64_t {
            if (own_r(t, f.g, c.g) == me) {
                const int r = rank_r(t, f.g, c.g);
                return r < 0 ? -1 : static_cast<int64_t>(r) * c.count;
            }
            if (!used_here) return 0;
            const Run* fr = rr[t].data() + ro[t][k];
            int64_t lo = 0, hi = 0;
            detail::child_runs(fr, ro[t][k + 1] - ro[t][k], f.g, lo, hi);
            for (int64_t i = lo; i < hi; ++i) {
                if (own_r(t + 1, fr[i].g, c.g >> 1) != me) continue;
                const BlockDesc* b = block_at(t, fr[i].g, c.g);
                return b == nullptr ? -1 : static_cast<int64_t>(b->cols) * c.count;
            }
            return 0;
        });
    }
    // left half, vector levels L+1 .. h+1 (runs of tree level t-1), likewise
    for (int t = L + 1; t >= h + 1; --t) {
        const bool used_here = t >= h + 2 && shared_w(t);
        layout(t - 1, W[t], [&](int k, const Run& f, const Run& c) -> int64_t {
            if (own_w(t, f.g, c.g) == me) {
                const int r = rank_w(t, f.g, c.g);
                return r < 0 ? -1 : static_cast<int64_t>(r) * f.count;
            }
            if (!used_here) return 0;
            const Run* cc = cr[t - 2].data() + co[t - 2][k];
            int64_t lo = 0, hi = 0;
            detail::child_runs(cc, co[t - 2][k + 1] - co[t - 2][k], c.g, lo, hi);
            for (int64_t i = lo; i < hi; ++i) {
                if (own_w(t - 1, f.g >> 1, cc[i].g) != me) continue;
                const BlockDesc* b = block_at(t - 1, f.g, cc[i].g);
                return b == nullptr ? -1 : static_cast<int64_t>(b->rows) * f.count;
            }
            return 0;
        });
    }
    // the right half's products of vector level h+1 that meet this rank's
    // left half's
    layout(h, Rmid, [&](int, const Run& f, const Run& c) -> int64_t {
        if (own_w(h + 1, f.g, c.g) != me) return 0;
        const int r = rank_w(h + 1, f.g, c.g);
        return r < 0 ? -1 : static_cast<int64_t>(r) * c.count;
    });
    if (bad_rank > 0) ok = 0;
    // the meet's rows that other ranks own, in order (intersection, row run,
    // column run, owner), staged before they go: where each intersection's
    // start
    std::vector<int64_t> stage_at(row_at.size(), 0);
    if (shared_out) {
#pragma omp parallel for schedule(dynamic, 16)
        for (int k = 0; k < count; ++k) {
            int64_t n = 0;
            for (int64_t fi = ro[h][k]; fi < ro[h][k + 1]; ++fi) {
                const Run& f = rr[h][fi];
                for (int64_t ci = co[h][k]; ci < co[h][k + 1]; ++ci) {
                    const Run& c = cr[h][ci];
                    if (own_w(h + 1, f.g, c.g) != me) continue;
                    by_owner(k, f.start, f.start + f.count, [&](int64_t r0, int64_t r1, int o) {
                        if (o != me) n += (r1 - r0) * c.count;
                    });
                }
            }
            stage_at[k + 1] = n;
        }
        for (int k = 0; k < count; ++k) stage_at[k + 1] += stage_at[k];
    }

    // ---- the GEMM stages, as extract()'s, over the products this rank makes
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
    DataType* dstage = nullptr;
    auto put = [](ExtractStage& s, int64_t i, const DataType* a, int lda, const DataType* b, int ldb, DataType* c,
                  int ldc, int m, int n, int kk) {
        s.a[i] = a; s.b[i] = b; s.c[i] = c;
        s.m[i] = m; s.n[i] = n; s.k[i] = kk;
        s.lda[i] = std::max(lda, 1); s.ldb[i] = std::max(ldb, 1); s.ldc[i] = std::max(ldc, 1);
    };
    auto put_bad = [&](ExtractStage& s, int64_t i) { put(s, i, nullptr, 1, nullptr, 1, nullptr, 1, 0, 0, 0); };
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
    if (shared_out) dstage = take(stage_at[count]);
    if (!room) ok = 0;
    const bool go = ok != 0;  // pointers into missing buffers are never followed: nothing launches unless all is well

    if (go) {
        // right half: vector level t -> t+1 through Kerl(t).  Where no split
        // falls among a column group's leaves (t <= low_col) one rank makes
        // all its products, whatever their rows.
        for (int t = 1; t <= h; ++t) {
            const bool whole = t <= low_col;
            build([&](int k) {
                      int64_t n = 0;
                      for (int64_t ci = co[t - 1][k]; ci < co[t - 1][k + 1]; ++ci) {
                          const int C = cr[t - 1][ci].g;
                          if (whole) {
                              if (own_r(t + 1, 0, C >> 1) == me) n += ro[t][k + 1] - ro[t][k];
                              continue;
                          }
                          for (int64_t fo = ro[t][k]; fo < ro[t][k + 1]; ++fo)
                              if (own_r(t + 1, rr[t][fo].g, C >> 1) == me) ++n;
                      }
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
                          if (whole && own_r(t + 1, 0, C >> 1) != me) continue;
                          const int cp = find(c_out, nc_out, C >> 1);
                          for (int64_t fo = 0; fo < nf_out; ++fo) {
                              const int F = f_out[fo].g;
                              if (!whole && own_r(t + 1, F, C >> 1) != me) continue;
                              const int fi = find(f_in, nf_in, F >> 1);
                              const BlockDesc* K = block_at(t, F, C);  // r_out x r_in
                              const int64_t a_in = R[t].at[R[t].base[k] + fi * nc_in + ci];
                              const int64_t a_out = R[t + 1].at[R[t + 1].base[k] + fo * nc_out + cp];
                              if (K == nullptr || a_in < 0 || a_out < 0) {
                                  ++bad;
                                  put_bad(s, i++);
                                  continue;
                              }
                              DataType* Cp = dR[t + 1] + a_out + static_cast<int64_t>(c_in[ci].start - c_out[cp].start) * K->rows;
                              put(s, i++, bf.values + K->offset, K->lda, dR[t] + a_in, K->cols, Cp, K->rows, K->rows,
                                  c_in[ci].count, K->cols);
                          }
                      }
                      return bad;
                  });
        }
        // left half: vector level t -> t-1 through Kerl(t-1).  Where no split
        // falls among a row group's leaves (L-t+2 <= low_row) one rank makes
        // all its products, whatever their columns.
        for (int t = L + 1; t >= h + 2; --t) {
            const bool whole = L - t + 2 <= low_row;
            build([&](int k) {
                      int64_t n = 0;
                      for (int64_t fi = ro[t - 1][k]; fi < ro[t - 1][k + 1]; ++fi) {
                          const int F = rr[t - 1][fi].g;
                          if (whole) {
                              if (own_w(t - 1, F >> 1, 0) == me) n += co[t - 2][k + 1] - co[t - 2][k];
                              continue;
                          }
                          for (int64_t cj = co[t - 2][k]; cj < co[t - 2][k + 1]; ++cj)
                              if (own_w(t - 1, F >> 1, cr[t - 2][cj].g) == me) ++n;
                      }
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
                          if (whole && own_w(t - 1, F >> 1, 0) != me) continue;
                          const int fp = find(f_out, nf_out, F >> 1);
                          for (int64_t cj = 0; cj < nc_out; ++cj) {
                              const int J = c_out[cj].g;
                              if (!whole && own_w(t - 1, F >> 1, J) != me) continue;
                              const int ci = find(c_in, nc_in, J >> 1);
                              const BlockDesc* K = block_at(t - 1, F, J);  // r_in x r_out
                              const int64_t a_in = W[t].at[W[t].base[k] + fi * nc_in + ci];
                              const int64_t a_out = W[t - 1].at[W[t - 1].base[k] + fp * nc_out + cj];
                              if (K == nullptr || a_in < 0 || a_out < 0) {
                                  ++bad;
                                  put_bad(s, i++);
                                  continue;
                              }
                              DataType* Cp = dW[t - 1] + a_out + (f_in[fi].start - f_out[fp].start);
                              put(s, i++, dW[t] + a_in, f_in[fi].count, bf.values + K->offset, K->lda, Cp,
                                  f_out[fp].count, f_in[fi].count, K->cols, K->rows);
                          }
                      }
                      return bad;
                  });
        }
        // the middle: the halves this rank meets, into its rows of the
        // result, and the rows other ranks own into the staging buffer
        build([&](int k) {
                  int64_t n = 0;
                  const int64_t nc = co[h][k + 1] - co[h][k];
                  for (int64_t fi = ro[h][k]; fi < ro[h][k + 1]; ++fi) {
                      const Run& f = rr[h][fi];
                      if (!shared_out) {
                          if (own_w(h + 1, f.g, 0) == me) n += nc;
                          continue;
                      }
                      for (int64_t ci = co[h][k]; ci < co[h][k + 1]; ++ci) {
                          if (own_w(h + 1, f.g, cr[h][ci].g) != me) continue;
                          by_owner(k, f.start, f.start + f.count, [&](int64_t, int64_t, int) { ++n; });
                      }
                  }
                  return n;
              },
              [&](int k, ExtractStage& s, int64_t i) {
                  int bad = 0;
                  const Run* fr = rr[h].data() + ro[h][k];
                  const int64_t nf = ro[h][k + 1] - ro[h][k];
                  const Run* cc = cr[h].data() + co[h][k];
                  const int64_t nc = co[h][k + 1] - co[h][k];
                  const int m = q.nloc[k];
                  int64_t staged = stage_at[k];
                  for (int64_t fi = 0; fi < nf; ++fi) {
                      if (!shared_out && own_w(h + 1, fr[fi].g, 0) != me) continue;
                      for (int64_t ci = 0; ci < nc; ++ci) {
                          if (shared_out && own_w(h + 1, fr[fi].g, cc[ci].g) != me) continue;
                          const int r = rank_w(h + 1, fr[fi].g, cc[ci].g);
                          const int64_t aw = W[h + 1].at[W[h + 1].base[k] + fi * nc + ci];
                          const int64_t am = Rmid.at[Rmid.base[k] + fi * nc + ci];
                          auto part = [&](int64_t r0, int64_t r1, int o) {
                              if (r < 0 || aw < 0 || am < 0) {
                                  ++bad;
                                  put_bad(s, i++);
                                  return;
                              }
                              const DataType* A = dW[h + 1] + aw + (r0 - fr[fi].start);
                              const int rows = static_cast<int>(r1 - r0);
                              if (o == me) {
                                  DataType* Cp = sorted + out_at[k] + (r0 - first[k]) + static_cast<int64_t>(cc[ci].start) * m;
                                  put(s, i++, A, fr[fi].count, dRmid + am, r, Cp, m, rows, cc[ci].count, r);
                              } else {
                                  put(s, i++, A, fr[fi].count, dRmid + am, r, dstage + staged, rows, rows, cc[ci].count, r);
                                  staged += static_cast<int64_t>(rows) * cc[ci].count;
                              }
                          };
                          if (shared_out) by_owner(k, fr[fi].start, fr[fi].start + fr[fi].count, part);
                          else part(fr[fi].start, fr[fi].start + fr[fi].count, me);
                      }
                  }
                  return bad;
              });
    }
    if (bad_block > 0) ok = 0;

    // ---- the exchanges, in the order they run, the same on every rank: the
    // right half's shared levels up, the left half's down, the midpoint, the
    // meet's rows other ranks own.  Each packs slots of one buffer, sends
    // them, and unpacks them into slots of another, in one order on both
    // ends; each end sizes a slot from its own factor.
    struct Exchange {
        std::vector<ExtractCopy> pack, unpack;
        std::vector<int64_t> scount, rcount;
        int64_t stotal = 0, rtotal = 0;
        size_t at_pack = 0, at_unpack = 0;
    };
    std::vector<Exchange> xs;
    std::vector<int> xr(static_cast<size_t>(h) + 2, -1), xw(static_cast<size_t>(L) + 2, -1);
    for (int t = 1; t <= h; ++t) {
        if (!shared_r(t)) continue;
        xr[static_cast<size_t>(t)] = static_cast<int>(xs.size());
        xs.emplace_back();
    }
    for (int t = L + 1; t >= h + 2; --t) {
        if (!shared_w(t)) continue;
        xw[static_cast<size_t>(t)] = static_cast<int>(xs.size());
        xs.emplace_back();
    }
    const int xmid = static_cast<int>(xs.size());
    xs.emplace_back();
    const int xout = shared_out ? static_cast<int>(xs.size()) : -1;
    if (shared_out) xs.emplace_back();
    for (Exchange& x : xs) {
        x.scount.assign(static_cast<size_t>(nproc), 0);
        x.rcount.assign(static_cast<size_t>(nproc), 0);
    }
    using Lists = std::vector<std::vector<ExtractCopy>>;
    auto finish = [&](Exchange& x, Lists& sends, Lists& recvs) {
        for (int p = 0; p < nproc; ++p) {
            for (ExtractCopy c : sends[static_cast<size_t>(p)]) {
                c.dst = x.stotal;
                x.stotal += c.n;
                x.scount[static_cast<size_t>(p)] += c.n;
                if (c.n > 0) x.pack.push_back(c);
            }
            for (ExtractCopy c : recvs[static_cast<size_t>(p)]) {
                c.src = x.rtotal;
                x.rtotal += c.n;
                x.rcount[static_cast<size_t>(p)] += c.n;
                if (c.n > 0) x.unpack.push_back(c);
            }
        }
    };
    // a shared level's products to their other users: `users(k, f, c, maker,
    // u)` writes those of slot (f, c) other than its maker into u
    auto level_exchange = [&](Exchange& x, const Level& s, int tl, auto maker_of, auto users) {
        Lists sends(static_cast<size_t>(nproc)), recvs(static_cast<size_t>(nproc));
        for (int k = 0; k < count; ++k) {
            const int64_t nF = ro[tl][k + 1] - ro[tl][k], nC = co[tl][k + 1] - co[tl][k];
            for (int64_t fi = 0; fi < nF; ++fi) {
                const Run& f = rr[tl][ro[tl][k] + fi];
                for (int64_t ci = 0; ci < nC; ++ci) {
                    const Run& c = cr[tl][co[tl][k] + ci];
                    const int maker = maker_of(f.g, c.g);
                    int u[2];
                    const int nu = users(k, f.g, c.g, maker, u);
                    if (nu == 0 || (maker != me && u[0] != me && (nu < 2 || u[1] != me))) continue;
                    const size_t idx = static_cast<size_t>(s.base[k] + fi * nC + ci);
                    for (int j = 0; j < nu; ++j) {
                        if (maker == me) sends[static_cast<size_t>(u[j])].push_back(ExtractCopy{s.at[idx], 0, s.size[idx]});
                        else if (u[j] == me) recvs[static_cast<size_t>(maker)].push_back(ExtractCopy{0, s.at[idx], s.size[idx]});
                    }
                }
            }
        }
        finish(x, sends, recvs);
    };
    if (go) {
        for (int t = 1; t <= h; ++t) {
            if (xr[static_cast<size_t>(t)] < 0) continue;
            level_exchange(
                xs[static_cast<size_t>(xr[static_cast<size_t>(t)])], R[t], t - 1,
                [&](int F, int C) { return own_r(t, F, C); },
                [&](int k, int F, int C, int maker, int* u) {
                    const Run* fr = rr[t].data() + ro[t][k];
                    int64_t lo = 0, hi = 0;
                    detail::child_runs(fr, ro[t][k + 1] - ro[t][k], F, lo, hi);
                    int n = 0;
                    for (int64_t i = lo; i < hi; ++i) {
                        const int x = own_r(t + 1, fr[i].g, C >> 1);
                        if (x != maker && (n == 0 || u[0] != x)) u[n++] = x;
                    }
                    return n;
                });
        }
        for (int t = L + 1; t >= h + 2; --t) {
            if (xw[static_cast<size_t>(t)] < 0) continue;
            level_exchange(
                xs[static_cast<size_t>(xw[static_cast<size_t>(t)])], W[t], t - 1,
                [&](int F, int C) { return own_w(t, F, C); },
                [&](int k, int F, int C, int maker, int* u) {
                    const Run* cc = cr[t - 2].data() + co[t - 2][k];
                    int64_t lo = 0, hi = 0;
                    detail::child_runs(cc, co[t - 2][k + 1] - co[t - 2][k], C, lo, hi);
                    int n = 0;
                    for (int64_t i = lo; i < hi; ++i) {
                        const int x = own_w(t - 1, F >> 1, cc[i].g);
                        if (x != maker && (n == 0 || u[0] != x)) u[n++] = x;
                    }
                    return n;
                });
        }
        // the midpoint: the right half's products of vector level h+1 to the
        // rank meeting them
        {
            Lists sends(static_cast<size_t>(nproc)), recvs(static_cast<size_t>(nproc));
            for (int k = 0; k < count; ++k) {
                const int64_t nf = ro[h][k + 1] - ro[h][k], nc = co[h][k + 1] - co[h][k];
                for (int64_t fi = 0; fi < nf; ++fi) {
                    const Run& f = rr[h][ro[h][k] + fi];
                    for (int64_t ci = 0; ci < nc; ++ci) {
                        const Run& c = cr[h][co[h][k] + ci];
                        const int src = own_r(h + 1, f.g, c.g), dst = own_w(h + 1, f.g, c.g);
                        const size_t idx = static_cast<size_t>(R[h + 1].base[k] + fi * nc + ci);
                        if (src == me) {
                            sends[static_cast<size_t>(dst)].push_back(ExtractCopy{R[h + 1].at[idx], 0, R[h + 1].size[idx]});
                        }
                        if (dst == me) {
                            recvs[static_cast<size_t>(src)].push_back(ExtractCopy{0, Rmid.at[idx], Rmid.size[idx]});
                        }
                    }
                }
            }
            finish(xs[static_cast<size_t>(xmid)], sends, recvs);
        }
        // the meet's rows other ranks own, from the staging buffer to their
        // places in the owners' results
        if (xout >= 0) {
            Lists sends(static_cast<size_t>(nproc)), recvs(static_cast<size_t>(nproc));
            for (int k = 0; k < count; ++k) {
                const int m = q.nloc[k];
                const int64_t mine0 = first[k], mine1 = first[k] + q.nloc[k];
                int64_t staged = stage_at[k];
                for (int64_t fi = ro[h][k]; fi < ro[h][k + 1]; ++fi) {
                    const Run& f = rr[h][fi];
                    for (int64_t ci = co[h][k]; ci < co[h][k + 1]; ++ci) {
                        const Run& c = cr[h][ci];
                        const int w = own_w(h + 1, f.g, c.g);
                        if (w == me) {
                            by_owner(k, f.start, f.start + f.count, [&](int64_t r0, int64_t r1, int o) {
                                if (o == me) return;
                                const int64_t n = (r1 - r0) * c.count;
                                sends[static_cast<size_t>(o)].push_back(ExtractCopy{staged, 0, n});
                                staged += n;
                            });
                            continue;
                        }
                        const int64_t r0 = std::max<int64_t>(f.start, mine0), r1 = std::min<int64_t>(f.start + f.count, mine1);
                        if (r0 >= r1) continue;
                        recvs[static_cast<size_t>(w)].push_back(ExtractCopy{
                            0, out_at[k] + (r0 - mine0) + static_cast<int64_t>(c.start) * m, (r1 - r0) * c.count, m,
                            static_cast<int>(r1 - r0)});
                    }
                }
            }
            finish(xs[static_cast<size_t>(xout)], sends, recvs);
        }
    }
    // the exchanges' buffers, one pair for all (they run one after another
    // on the stream), with the others: a rank without room declines with the
    // group
    int64_t max_s = 0, max_r = 0, stotal = 0, rtotal = 0;
    for (const Exchange& x : xs) {
        max_s = std::max(max_s, x.stotal);
        max_r = std::max(max_r, x.rtotal);
        stotal += x.stotal;
        rtotal += x.rtotal;
    }
    DataType* dsend = take(max_s);
    DataType* drecv = take(max_r);
    if (!room) ok = 0;
    // both ends of every exchange must agree on what moves (they size it from
    // different factors)
    {
        const int nx = static_cast<int>(xs.size());
        std::vector<int64_t> mine(static_cast<size_t>(nproc) * nx, 0), told(mine.size(), 0);
        for (int p = 0; p < nproc; ++p)
            for (int e = 0; e < nx; ++e) mine[static_cast<size_t>(p) * nx + e] = xs[static_cast<size_t>(e)].scount[static_cast<size_t>(p)];
        MPI_Alltoall(mine.data(), nx, MPI_INT64_T, told.data(), nx, MPI_INT64_T, comm);
        for (int p = 0; p < nproc; ++p) {
            for (int e = 0; e < nx; ++e) {
                const Exchange& x = xs[static_cast<size_t>(e)];
                if (told[static_cast<size_t>(p) * nx + e] != x.rcount[static_cast<size_t>(p)] && go) ok = 0;
                if (x.scount[static_cast<size_t>(p)] > INT32_MAX || x.rcount[static_cast<size_t>(p)] > INT32_MAX) ok = 0;
            }
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
        for (int64_t ci = co[0][k]; ci < co[0][k + 1]; ++ci) nv += q.col_owner[cr[0][ci].g] == me ? 1 : 0;
        for (int64_t fi = ro[L][k]; fi < ro[L][k + 1]; ++fi) nu += q.row_owner[rr[L][fi].g] == me ? 1 : 0;
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
    for (Exchange& x : xs) {
        x.at_pack = meta.alloc(x.pack.size() * sizeof(ExtractCopy));
        x.at_unpack = meta.alloc(x.unpack.size() * sizeof(ExtractCopy));
    }
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
    for (const Exchange& x : xs) {
        if (!x.pack.empty()) std::memcpy(meta.host(x.at_pack), x.pack.data(), x.pack.size() * sizeof(ExtractCopy));
        if (!x.unpack.empty()) std::memcpy(meta.host(x.at_unpack), x.unpack.data(), x.unpack.size() * sizeof(ExtractCopy));
    }
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
            if (q.col_owner[run.g] != me) continue;
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
            if (q.row_owner[run.g] != me) continue;
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

    // ---- launch: the halves level by level, each shared level's exchange
    // before the level that uses it, the midpoint, the middle, the rows other
    // ranks own, the scatter
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
    // pack, exchange between the devices (alltoallv_device: NCCL where the
    // group has it), unpack
    auto exchange = [&](const Exchange& x, const DataType* from, DataType* to) {
        launch_extract_copies(reinterpret_cast<const ExtractCopy*>(dm + x.at_pack), static_cast<int>(x.pack.size()),
                              from, dsend, stream);
        alltoallv_device(comm, dsend, x.scount.data(), drecv, x.rcount.data(), stream);
        launch_extract_copies(reinterpret_cast<const ExtractCopy*>(dm + x.at_unpack), static_cast<int>(x.unpack.size()),
                              drecv, to, stream);
    };
    for (int t = 1; t <= h; ++t) {
        if (xr[static_cast<size_t>(t)] >= 0) exchange(xs[static_cast<size_t>(xr[static_cast<size_t>(t)])], dR[t], dR[t]);
        run_stage(static_cast<size_t>(t - 1));
    }
    for (int t = L + 1; t >= h + 2; --t) {
        if (xw[static_cast<size_t>(t)] >= 0) exchange(xs[static_cast<size_t>(xw[static_cast<size_t>(t)])], dW[t], dW[t]);
        run_stage(static_cast<size_t>(h + L + 1 - t));
    }
    const auto tm = clock::now();
    exchange(xs[static_cast<size_t>(xmid)], dR[h + 1], dRmid);
    prof.mid = std::chrono::duration<double>(clock::now() - tm).count();
    run_stage(static_cast<size_t>(L));
    if (xout >= 0) exchange(xs[static_cast<size_t>(xout)], dstage, sorted);
    prof.exchanges = static_cast<int>(xs.size());
    prof.sent = stotal;
    prof.received = rtotal;
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
