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

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
