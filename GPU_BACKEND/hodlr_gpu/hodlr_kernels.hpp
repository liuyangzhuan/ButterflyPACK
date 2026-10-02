#pragma once
// Small batched kernels of the HODLR GPU backend (hodlr_kernels.cu), for
// double and fmm::gpu::dcomplex.  Matrices are column major; each item of a
// launch describes one matrix (or pair), and a launch runs all its items.

#ifdef H2_HAVE_GPU

#include "gpu_common/device_scalar.hpp"

#include <cuda_runtime.h>

#include <cstdint>

namespace bpack {
namespace gpu {

// Row interchanges of LAPACK getrf (1-based pivots) applied to b (m x n):
// for i = 0 .. m-1, swap rows i and ipiv[i] - 1 (LAPACK laswp, forward).
template<typename T>
struct RowSwapItem {
    T* b;
    int ldb;
    int m;
    int n;
    const int* ipiv;
};
template<typename T>
void launch_row_swaps(const RowSwapItem<T>* items, int count, int max_n, cudaStream_t stream);

// The column permutation of LAPACK's column-pivoted QR (?geqp3 on small
// matrices: ?laqp2, with its norm downdating) of the r x L matrix A, first p
// entries (0-based; pivots for the first min(r, L) of them), where
// A(i, l) = a[i + l * lda], or a[l + i * lda] when transposed.  The columns
// listed at mask_ptr (device), or at meta + mask_offset when mask_ptr is
// null (nmask int32), count as zero.  work: r x L scratch of the item,
// norms: 3 L doubles, perm: L ints; out: p ints.
template<typename T>
struct PivotItem {
    const T* a;
    int lda;
    int transposed;
    int r;
    int L;
    int p;
    int nmask;
    int64_t mask_offset;
    const int* mask_ptr;
    T* work;
    double* norms;
    int* perm;
    int* out;
};
// one thread block per item (for many items of moderate L)
template<typename T>
void launch_qrcp_pivots(const PivotItem<T>* items, int count, const char* meta, cudaStream_t stream);
// the same pivots with the columns of each item spread over thread blocks,
// two launches per step (for items of large L, L >= 1024); max_r, max_L,
// max_steps (min(p, r, L)) and max_p over the items
template<typename T>
void launch_qrcp_pivots_wide(const PivotItem<T>* items, int count, int max_r, int max_L, int max_steps, int max_p,
                             const char* meta, cudaStream_t stream);

// The first n columns of a (Householder vectors of ?geqrf below the
// diagonal) as a unit lower trapezoidal matrix: a(i, j) = 0 for i < j and
// a(j, j) = 1, j < n.
template<typename T>
struct UnitLowerItem {
    T* a;
    int lda;
    int n;
};
template<typename T>
void launch_unit_lower(const UnitLowerItem<T>* items, int count, int max_n, cudaStream_t stream);

// a (n x n): zero below the diagonal (the R of a QR from its output)
template<typename T>
struct TriuItem {
    T* a;
    int lda;
    int n;
};
template<typename T>
void launch_zero_lower(const TriuItem<T>* items, int count, int max_n, cudaStream_t stream);

// c(:, j) = a(:, j) s[j] (m x n; s may be null: 1), conjugated if conj
template<typename T>
struct ColScaleItem {
    T* c;
    int ldc;
    const T* a;
    int lda;
    int m;
    int n;
    const double* s;
    int conj;
};
template<typename T>
void launch_col_scale(const ColScaleItem<T>* items, int count, int max_m, int max_n, cudaStream_t stream);

// t (k x k) = T^-1 of the compact WY form of ?geqrf's Q = I - Y T Y^H:
// striu(gram) + diag(1 / tau), zero below the diagonal (gram = Y^H Y)
template<typename T>
struct TinvItem {
    T* t;
    const T* gram;
    int ldg;
    const T* tau;
    int k;
};
template<typename T>
void launch_tinv(const TinvItem<T>* items, int count, int max_k, cudaStream_t stream);

// c = alpha a + beta b (m x n); b may be null (then c = alpha a), and c may be
// a or b.
template<typename T>
struct AxpbyItem {
    T* c;
    int ldc;
    const T* a;
    int lda;
    const T* b;
    int ldb;
    int m;
    int n;
    T alpha;
    T beta;
};
template<typename T>
void launch_axpby(const AxpbyItem<T>* items, int count, int max_m, int max_n, cudaStream_t stream);

// a(i, i) += value for i < n
template<typename T>
struct DiagItem {
    T* a;
    int lda;
    int n;
    T value;
};
template<typename T>
void launch_add_diagonal(const DiagItem<T>* items, int count, int max_n, cudaStream_t stream);

// a = (a + a^T) / 2 (n x n, plain transpose)
template<typename T>
struct SquareItem {
    T* a;
    int lda;
    int n;
};
template<typename T>
void launch_symmetrize(const SquareItem<T>* items, int count, int max_n, cudaStream_t stream);

// out[0 .. n-1] = diagonal of a (for the log-determinant on the host)
template<typename T>
struct DiagCopyItem {
    const T* a;
    int lda;
    int n;
    T* out;
};
template<typename T>
void launch_copy_diagonal(const DiagCopyItem<T>* items, int count, int max_n, cudaStream_t stream);

// *norm = Frobenius norm of a (m x n), summed in a fixed order
template<typename T>
struct NormItem {
    const T* a;
    int lda;
    int m;
    int n;
    double* norm;
};
template<typename T>
void launch_fnorm(const NormItem<T>* items, int count, cudaStream_t stream);

// out (n x k, leading dimension n) = the rows rows[0..n) (from 0) of a
// (k columns, leading dimension lda)
template<typename T>
void launch_gather_rows(const T* a, int lda, int k, const int* rows, int n, T* out, cudaStream_t stream);

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
