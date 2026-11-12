#pragma once
// ---------------------------------------------------------------------------
// LU with partial pivoting of one large dense block on the calling process's
// OpenMP team, over the sequential BLAS.  The root's Schur complement was one
// getrf call, so one of the node's cores: 304 s of a 789 s factorization at
// 160^3 (PPW-10 Helmholtz, 64 CPU nodes, 2026-11-01).
//
// Recursive as LAPACK's getrf2: the left half of the columns factored, its
// row swaps and the triangular solve applied to the right half, the trailing
// block updated by one GEMM of inner dimension n/2, the trailing block
// factored, its swaps applied back to the left half.  Panels of at most
// `leaf` columns go to the sequential getrf; the swaps, solves and GEMMs are
// split across the team in column chunks and tiles.  (A right-looking split
// in 128-column panels made every update a GEMM of inner dimension 128, which
// the 128 cores of a node run at memory bandwidth: 2.8 of 4.4 s at n = 9559.)
//
// The result is getrf's: L (unit diagonal) below and U on and above the
// diagonal in place, 1-based pivots over the whole block, info the first
// exactly zero pivot (1-based) or 0, so getrs and the logdet read it as
// before.  The pivots are partial pivoting's; the rounding is not getrf's bit
// for bit.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <complex>
#include <cstdint>
#include <type_traits>

#include <omp.h>

#include "blas_declare.hpp"

namespace fmm {

template <typename DataType>
inline void getrf_sequential(int m, int n, DataType* a, int lda, int* ipiv, int* info) {
    if constexpr (std::is_same_v<DataType, double>) {
        dgetrf_(&m, &n, a, &lda, ipiv, info);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgetrf_(&m, &n, a, &lda, ipiv, info);
    } else {
        static_assert(sizeof(DataType) == 0, "getrf_sequential: double or complex<double>");
    }
}

template <typename DataType>
inline void laswp_sequential(int n, DataType* a, int lda, int k1, int k2, const int* ipiv) {
    const int one = 1;
    if constexpr (std::is_same_v<DataType, double>) {
        dlaswp_(&n, a, &lda, &k1, &k2, ipiv, &one);
    } else {
        zlaswp_(&n, a, &lda, &k1, &k2, ipiv, &one);
    }
}

namespace parallel_getrf_detail {

constexpr int leaf = 64;         // widest panel the sequential getrf takes (64: 3.6 s at n = 9559, 128: 4.0, 256: 4.1)
constexpr int swap_chunk = 256;  // columns a task swaps
constexpr int solve_chunk = 64;  // least columns a task solves

inline int64_t ceil_div(int64_t a, int64_t b) { return (a + b - 1) / b; }

// rows k1..k2 (1-based) of ipiv swapped in the n columns at a
template <typename DataType>
void swap_rows(int n, DataType* a, int lda, int k1, int k2, const int* ipiv, int threads) {
    const int chunks = static_cast<int>(ceil_div(n, swap_chunk));
    #pragma omp parallel for num_threads(std::min(threads, chunks)) schedule(dynamic)
    for (int c = 0; c < chunks; ++c) {
        const int c0 = c * swap_chunk;
        laswp_sequential(std::min(swap_chunk, n - c0), a + static_cast<int64_t>(c0) * lda, lda, k1, k2, ipiv);
    }
}

// the n1 rows at a12 = L11^{-1} a12 (L11 unit lower at a11), the n2 columns split
template <typename DataType>
void solve_unit_lower(int n1, int n2, const DataType* a11, DataType* a12, int lda, int threads) {
    const DataType one(1.0);
    const int width = static_cast<int>(std::max<int64_t>(solve_chunk, ceil_div(n2, threads)));
    const int chunks = static_cast<int>(ceil_div(n2, width));
    #pragma omp parallel for num_threads(std::min(threads, chunks)) schedule(dynamic)
    for (int c = 0; c < chunks; ++c) {
        const int c0 = c * width;
        const int w = std::min(width, n2 - c0);
        trsm_("L", "L", "N", "U", &n1, &w, &one, a11, &lda, a12 + static_cast<int64_t>(c0) * lda, &lda);
    }
}

// c -= a b (rows x cols, inner k), in tiles of at most 256 columns and 64..256 rows, about two a thread
template <typename DataType>
void update_trailing(int rows, int cols, int k, const DataType* a, const DataType* b, DataType* c,
                     int lda, int threads) {
    const DataType one(1.0);
    const DataType minus_one(-1.0);
    const int tile_cols = cols <= 512 ? cols : 256;
    const int64_t col_tiles = ceil_div(cols, tile_cols);
    int tile_rows = 256;
    while (tile_rows > 64 && ceil_div(rows, tile_rows) * col_tiles < 2 * threads) tile_rows /= 2;
    const int64_t row_tiles = ceil_div(rows, tile_rows);
    const int64_t tiles = row_tiles * col_tiles;
    #pragma omp parallel for num_threads(static_cast<int>(std::min<int64_t>(threads, tiles))) schedule(dynamic)
    for (int64_t t = 0; t < tiles; ++t) {
        const int r0 = static_cast<int>((t % row_tiles) * tile_rows);
        const int c0 = static_cast<int>((t / row_tiles) * tile_cols);
        const int mr = std::min(tile_rows, rows - r0);
        const int nc = std::min(tile_cols, cols - c0);
        gemm_("N", "N", &mr, &nc, &k, &minus_one, a + r0, &lda,
              b + static_cast<int64_t>(c0) * lda, &lda, &one,
              c + r0 + static_cast<int64_t>(c0) * lda, &lda);
    }
}

// getrf2 of the m x n panel at a (m >= n)
template <typename DataType>
void factor(int m, int n, DataType* a, int lda, int* ipiv, int* info, int threads) {
    if (n <= leaf) {
        getrf_sequential(m, n, a, lda, ipiv, info);
        return;
    }
    auto at = [&](int64_t row, int64_t col) { return a + row + col * static_cast<int64_t>(lda); };
    const int n1 = std::min(m, n) / 2;
    const int n2 = n - n1;
    int left_info = 0;
    factor(m, n1, a, lda, ipiv, &left_info, threads);
    if (left_info != 0 && *info == 0) *info = left_info;
    swap_rows(n2, at(0, n1), lda, 1, n1, ipiv, threads);
    solve_unit_lower(n1, n2, at(0, 0), at(0, n1), lda, threads);
    update_trailing(m - n1, n2, n1, at(n1, 0), at(0, n1), at(n1, n1), lda, threads);
    int right_info = 0;
    factor(m - n1, n2, at(n1, n1), lda, ipiv + n1, &right_info, threads);
    if (right_info != 0 && *info == 0) *info = right_info + n1;
    const int pivots = std::min(m, n);
    for (int i = n1; i < pivots; ++i) ipiv[i] += n1;
    swap_rows(n1, at(0, 0), lda, n1 + 1, pivots, ipiv, threads);
}

} // namespace parallel_getrf_detail

// The n x n block `a` (leading dimension lda) factored in place; ipiv holds n
// entries.  Outside a parallel region with more than one thread and a block
// wider than two leaves, the team does it; otherwise one getrf call.
template <typename DataType>
void parallel_getrf(int n, DataType* a, int lda, int* ipiv, int* info) {
    *info = 0;
    const int threads = omp_in_parallel() ? 1 : omp_get_max_threads();
    if (threads <= 1 || n <= 2 * parallel_getrf_detail::leaf) {
        getrf_sequential(n, n, a, lda, ipiv, info);
        return;
    }
    parallel_getrf_detail::factor(n, n, a, lda, ipiv, info, threads);
}

} // namespace fmm
