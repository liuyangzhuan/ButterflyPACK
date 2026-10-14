// Kernels of the device solve.  See solve_kernels.hpp.

#include "solve_kernels.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace fmm {
namespace gpu {

namespace {

constexpr int kThreads = 256;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

// The block's dynamic shared memory as an array of T.
template<typename T>
__device__ __forceinline__ T* shared_array() {
    extern __shared__ __align__(16) unsigned char solve_shared[];
    return reinterpret_cast<T*>(solve_shared);
}

// The right-hand sides go through a box in passes of up to kMaxCols columns:
// each pass reads the box's factors once and applies every element to all of
// its columns (column q's x_S, x_R and, where it fits, x_N in shared memory at
// q times the column's stride).  A column's arithmetic is the single-column
// order - the same partial sums, summed the same way - so its result does not
// depend on nrhs or on how the columns are grouped.
template<typename T>
constexpr int kMaxCols = std::is_same_v<T, double> ? 8 : 4;

// Whether unrolled column q is in the pass (column 0 always is: the compiler
// drops the test, and a one-column launch is the single-column code).
__device__ __forceinline__ bool live(int q, int nc) { return q == 0 || q < nc; }

// s_q += sum_{j0 <= j < j1} a[j * lda] x_q[j] for the columns x_q = x + q * xs
// (q < nc), in the order of j; the loads of four j are issued before their
// products, so they are in flight together.
template<typename T, int CB>
__device__ __forceinline__ void accumulate_row(const T* a, int64_t lda, const T* x, int xs, int nc, int j0, int j1,
                                               T s[CB]) {
    int j = j0;
    for (; j + 3 < j1; j += 4) {
        const T a0 = a[j * lda], a1 = a[(j + 1) * lda], a2 = a[(j + 2) * lda], a3 = a[(j + 3) * lda];
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) {
                const T* xq = x + q * xs;
                s[q] += a0 * xq[j];
                s[q] += a1 * xq[j + 1];
                s[q] += a2 * xq[j + 2];
                s[q] += a3 * xq[j + 3];
            }
        }
    }
    for (; j < j1; ++j) {
        const T a0 = a[j * lda];
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) s[q] += a0 * x[q * xs + j];
        }
    }
}

// row_dot for the columns x_q = x + q * xs (q < nc): out[q] = sum_j a[j * lda]
// x_q[j], four independent partial sums a column (loads in flight)
template<typename T, int CB>
__device__ __forceinline__ void row_dot_cols(const T* a, int lda, const T* x, int xs, int nc, int n, T out[CB]) {
    T s0[CB], s1[CB], s2[CB], s3[CB];
#pragma unroll
    for (int q = 0; q < CB; ++q) {
        s0[q] = T(0.0);
        s1[q] = T(0.0);
        s2[q] = T(0.0);
        s3[q] = T(0.0);
    }
    const int64_t ld = lda;
    int j = 0;
    for (; j + 3 < n; j += 4) {
        const T a0 = a[j * ld], a1 = a[(j + 1) * ld], a2 = a[(j + 2) * ld], a3 = a[(j + 3) * ld];
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) {
                const T* xq = x + q * xs;
                s0[q] += a0 * xq[j];
                s1[q] += a1 * xq[j + 1];
                s2[q] += a2 * xq[j + 2];
                s3[q] += a3 * xq[j + 3];
            }
        }
    }
    for (; j < n; ++j) {
        const T a0 = a[j * ld];
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) s0[q] += a0 * x[q * xs + j];
        }
    }
#pragma unroll
    for (int q = 0; q < CB; ++q) out[q] = (s0[q] + s1[q]) + (s2[q] + s3[q]);
}

__device__ __forceinline__ double warp_sum(double v) {
    for (int offset = 16; offset > 0; offset >>= 1) v += __shfl_down_sync(0xffffffffu, v, offset);
    return v;
}
__device__ __forceinline__ dcomplex warp_sum(dcomplex v) { return dcomplex(warp_sum(v.re), warp_sum(v.im)); }

// ---- one column at a time: the kernels as they were before the passes (each
// column in turn, the factors read again for each).  A launch whose passes
// would take one column - nrhs 1, or boxes too large for two columns of
// shared memory - runs these, so its code is what it was; the passes'
// kernels further down take two columns or more (and one column when x_N is
// in the work buffer: they stream it through shared memory in chunks).

// sum_j a[j * lda] x[j], four independent partial sums (loads in flight)
template<typename T>
__device__ __forceinline__ T row_dot(const T* a, int lda, const T* x, int n) {
    T s0 = T(0.0), s1 = T(0.0), s2 = T(0.0), s3 = T(0.0);
    const int64_t ld = lda;
    int j = 0;
    for (; j + 3 < n; j += 4) {
        s0 += a[j * ld] * x[j];
        s1 += a[(j + 1) * ld] * x[j + 1];
        s2 += a[(j + 2) * ld] * x[j + 2];
        s3 += a[(j + 3) * ld] * x[j + 3];
    }
    for (; j < n; ++j) s0 += a[j * ld] * x[j];
    return (s0 + s1) + (s2 + s3);
}

template<typename T>
__global__ void __launch_bounds__(kThreads) solve_forward_one_kernel(const SolveBox* boxes, const int* wave, T* vec,
                                                                 T* work, int nrhs) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    T* xs = sh;
    T* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (bT != nullptr) {  // x_R -= T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                const T* col = bT + static_cast<int64_t>(j) * b.k;
                T s = T(0.0);
                for (int i = lane; i < b.k; i += 32) s += col[i] * xs[i];
                s = warp_sum(s);
                if (lane == 0) xr[j] -= s;
            }
            __syncthreads();
        }
        for (int j = tid; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
        if (bxsr != nullptr) {  // x_S += X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s = T(0.0);
                for (int j = 0; j < b.r; ++j) s += bxsr[i + static_cast<int64_t>(j) * b.k] * xr[j];
                v[b.skel[i]] = xs[i] + s;
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates X_NR x_R
            T* u = work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                T s = T(0.0);
                for (int j = 0; j < b.r; ++j) s += bxnr[i + static_cast<int64_t>(j) * b.ntot] * xr[j];
                u[i] = s;
            }
        }
        __syncthreads();
    }
}

// Columns j .. j+3 of A^T x (A: rows x cols, ld): lane partial sums.
template<typename T>
__device__ __forceinline__ void cols4_dot(const T* a, int64_t ld, const T* x, int rows, int j, int cols, int lane,
                                          T out[4]) {
    T s0 = T(0.0), s1 = T(0.0), s2 = T(0.0), s3 = T(0.0);
    const T* c0 = a + j * ld;
    const bool has1 = j + 1 < cols, has2 = j + 2 < cols, has3 = j + 3 < cols;
    for (int i = lane; i < rows; i += 32) {
        const T xi = x[i];
        s0 += c0[i] * xi;
        if (has1) s1 += c0[ld + i] * xi;
        if (has2) s2 += c0[2 * ld + i] * xi;
        if (has3) s3 += c0[3 * ld + i] * xi;
    }
    out[0] = warp_sum(s0);
    out[1] = warp_sum(s1);
    out[2] = warp_sum(s2);
    out[3] = warp_sum(s3);
}

// x_N in shared memory (after x_S, x_R) when kSharedXn, else in the work buffer
template<typename T, bool kSharedXn>
__global__ void __launch_bounds__(kThreads) solve_backward_one_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                                  const int* wave, T* vec, const T* ghost, T* work,
                                                                  int nrhs) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    T* xs = sh;
    T* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        T* xn = kSharedXn ? sh + b.k + b.r : work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        // x_N: the neighbors' current rows, in X_NR's row order
        for (int a = 0; a < b.nslots; ++a) {
            const SolveSlot sl = slots[b.slot0 + a];
            const T* src = (sl.ghost ? ghost : vec) + sl.vec * nrhs + static_cast<int64_t>(c) * sl.n;
            for (int i = tid; i < sl.count; i += blockDim.x) xn[sl.row0 + i] = sl.full ? src[i] : src[sl.skel[i]];
        }
        __syncthreads();
        // x_R += X_SR^T x_S, then x_R += X_NR^T x_N (two products, as on the host);
        // a warp takes four columns at a time
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            T s1[4] = {T(0.0), T(0.0), T(0.0), T(0.0)}, s2[4];
            if (bxsr != nullptr) cols4_dot(bxsr, b.k, xs, b.k, j, b.r, lane, s1);
            cols4_dot(bxnr, b.ntot, xn, b.ntot, j, b.r, lane, s2);
            if (lane == 0) {
                for (int q = 0; q < 4 && j + q < b.r; ++q) xr[j + q] = (xr[j + q] + s1[q]) + s2[q];
            }
        }
        __syncthreads();
        for (int j = tid; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
        if (bT != nullptr) {  // x_S -= T x_R
            for (int i = tid; i < b.k; i += blockDim.x) v[b.skel[i]] = xs[i] - row_dot(bT + i, b.k, xr, b.r);
        }
        __syncthreads();
    }
}

// LAPACK getrs ('N') with the 1-based interchanges of getrf.
template<typename T>
__global__ void __launch_bounds__(kThreads) solve_diagonal_one_kernel(const SolveBox* boxes, const int* list, T* vec,
                                                                  int nrhs) {
    T* xr = shared_array<T>();
    const SolveBox b = boxes[list[blockIdx.x]];
    const T* lu = static_cast<const T*>(b.lu);
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int j = tid; j < r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (tid == 0) {
            for (int i = 0; i < r; ++i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const T t = xr[i];
                    xr[i] = xr[p];
                    xr[p] = t;
                }
            }
        }
        for (int i = 0; i < r; ++i) {  // unit lower L
            __syncthreads();
            const T xi = xr[i];
            const T* col = lu + static_cast<int64_t>(i) * r;
            for (int l = i + 1 + tid; l < r; l += blockDim.x) xr[l] -= col[l] * xi;
        }
        for (int i = r - 1; i >= 0; --i) {  // upper U
            __syncthreads();
            if (tid == 0) xr[i] = xr[i] / lu[i + static_cast<int64_t>(i) * r];
            __syncthreads();
            const T xi = xr[i];
            const T* col = lu + static_cast<int64_t>(i) * r;
            for (int l = tid; l < i; l += blockDim.x) xr[l] -= col[l] * xi;
        }
        __syncthreads();
        for (int j = tid; j < r; j += blockDim.x) v[b.red[j]] = xr[j];
        __syncthreads();
    }
}

// ---- the multiply

template<typename T, bool kSharedXn>
__global__ void __launch_bounds__(kThreads) mul_forward_one_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                               const int* wave, T* vec, const T* ghost, T* work,
                                                               int nrhs) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    T* xs = sh;
    T* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        T* xn = kSharedXn ? sh + b.k + b.r : work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (bT != nullptr) {  // x_S += T x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s = T(0.0);
                for (int j = 0; j < b.r; ++j) s += bT[i + static_cast<int64_t>(j) * b.k] * xr[j];
                xs[i] += s;
                v[b.skel[i]] = xs[i];
            }
        }
        for (int a = 0; a < b.nslots; ++a) {  // x_N, in X_NR's row order
            const SolveSlot sl = slots[b.slot0 + a];
            const T* src = (sl.ghost ? ghost : vec) + sl.vec * nrhs + static_cast<int64_t>(c) * sl.n;
            for (int i = tid; i < sl.count; i += blockDim.x) xn[sl.row0 + i] = sl.full ? src[i] : src[sl.skel[i]];
        }
        __syncthreads();
        // x_R -= X_SR^T x_S, then x_R -= X_NR^T x_N (two products, as on the host)
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            T s1[4] = {T(0.0), T(0.0), T(0.0), T(0.0)}, s2[4];
            if (bxsr != nullptr) cols4_dot(bxsr, b.k, xs, b.k, j, b.r, lane, s1);
            cols4_dot(bxnr, b.ntot, xn, b.ntot, j, b.r, lane, s2);
            if (lane == 0) {
                for (int q = 0; q < 4 && j + q < b.r; ++q) xr[j + q] = (xr[j + q] - s1[q]) - s2[q];
            }
        }
        __syncthreads();
        for (int j = tid; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
        __syncthreads();
    }
}

// x_R = P L U x_R (LAPACK trmm U, trmm unit L, laswp with the interchanges reversed)
template<typename T>
__global__ void __launch_bounds__(kThreads) mul_diagonal_one_kernel(const SolveBox* boxes, const int* list, T* vec,
                                                                int nrhs) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[list[blockIdx.x]];
    const T* lu = static_cast<const T*>(b.lu);
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    T* x = sh;
    T* y = sh + r;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int j = tid; j < r; j += blockDim.x) x[j] = v[b.red[j]];
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // y = U x
            T s = T(0.0);
            for (int j = i; j < r; ++j) s += lu[i + static_cast<int64_t>(j) * r] * x[j];
            y[i] = s;
        }
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // x = L y (unit diagonal)
            T s = y[i];
            for (int j = 0; j < i; ++j) s += lu[i + static_cast<int64_t>(j) * r] * y[j];
            x[i] = s;
        }
        __syncthreads();
        if (tid == 0) {
            for (int i = r - 1; i >= 0; --i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const T t = x[i];
                    x[i] = x[p];
                    x[p] = t;
                }
            }
        }
        __syncthreads();
        for (int j = tid; j < r; j += blockDim.x) v[b.red[j]] = x[j];
        __syncthreads();
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) mul_backward_one_kernel(const SolveBox* boxes, const int* wave, T* vec,
                                                                T* work, int nrhs) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    T* xs = sh;
    T* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (bxsr != nullptr) {  // x_S -= X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s = T(0.0);
                for (int j = 0; j < b.r; ++j) s += bxsr[i + static_cast<int64_t>(j) * b.k] * xr[j];
                xs[i] -= s;
                v[b.skel[i]] = xs[i];
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates -X_NR x_R
            T* u = work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                T s = T(0.0);
                for (int j = 0; j < b.r; ++j) s += bxnr[i + static_cast<int64_t>(j) * b.ntot] * xr[j];
                u[i] = -s;
            }
        }
        __syncthreads();
        if (bT != nullptr) {  // x_R += T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                const T* col = bT + static_cast<int64_t>(j) * b.k;
                T s = T(0.0);
                for (int i = lane; i < b.k; i += 32) s += col[i] * xs[i];
                s = warp_sum(s);
                if (lane == 0) v[b.red[j]] = xr[j] + s;
            }
        }
        __syncthreads();
    }
}

// ---- passes of several columns

// sum_i a[i] x_q[i] for the columns x_q = x + q * xs (q < nc), lanes striding
// i; lane 0 gets the sums
template<typename T, int CB>
__device__ __forceinline__ void col_dot_cols(const T* a, const T* x, int xs, int nc, int rows, int lane, T s[CB]) {
#pragma unroll
    for (int q = 0; q < CB; ++q) s[q] = T(0.0);
    for (int i = lane; i < rows; i += 32) {
        const T t = a[i];
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) s[q] += t * x[q * xs + i];
        }
    }
#pragma unroll
    for (int q = 0; q < CB; ++q) {
        if (live(q, nc)) s[q] = warp_sum(s[q]);
    }
}

// Columns j .. j+3 of A^T x_q for the columns x_q = x + q * xs (q < nc; A:
// rows x cols, ld): lane partial sums, summed over the warp (lane 0 gets them)
template<typename T, int CB>
__device__ __forceinline__ void cols4_dot_cols(const T* a, int64_t ld, const T* x, int64_t xs, int nc, int rows, int j,
                                               int cols, int lane, T s[4][CB]) {
#pragma unroll
    for (int p = 0; p < 4; ++p) {
#pragma unroll
        for (int q = 0; q < CB; ++q) s[p][q] = T(0.0);
    }
    const T* c0 = a + j * ld;
    const bool has1 = j + 1 < cols, has2 = j + 2 < cols, has3 = j + 3 < cols;
    for (int i = lane; i < rows; i += 32) {
        const T a0 = c0[i];
        const T a1 = has1 ? c0[ld + i] : T(0.0);
        const T a2 = has2 ? c0[2 * ld + i] : T(0.0);
        const T a3 = has3 ? c0[3 * ld + i] : T(0.0);
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) {
                const T xi = x[q * xs + i];
                s[0][q] += a0 * xi;
                if (has1) s[1][q] += a1 * xi;
                if (has2) s[2][q] += a2 * xi;
                if (has3) s[3][q] += a3 * xi;
            }
        }
    }
#pragma unroll
    for (int p = 0; p < 4; ++p) {
#pragma unroll
        for (int q = 0; q < CB; ++q) {
            if (live(q, nc)) s[p][q] = warp_sum(s[p][q]);
        }
    }
}

// x_R +/-= X_NR^T x_N for the pass's columns when x_N is in the work buffer
// (column q at xn + q * ntot; x_R of column q at xr + q * stride).  The rows
// of x_N go through shared memory in chunks of kChunkRows, and a warp keeps
// its four columns' lane sums over all the chunks: lane l takes the rows
// i = l (mod 32) in increasing order, the single-column order.  x_N is read
// once a round of the warps' columns instead of once for every four.
constexpr int kChunkRows = 256;

template<typename T, int CB, bool kSubtract>
__device__ __forceinline__ void xnr_t_chunked(const SolveBox& b, const T* bxnr, const T* xn, T* chunk, T* xr,
                                              int stride, int nc) {
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarps = blockDim.x >> 5;
    const int64_t ld = b.ntot;
    for (int jr = 0; jr < b.r; jr += 4 * nwarps) {
        const int j = jr + 4 * warp;
        const bool has0 = j < b.r, has1 = j + 1 < b.r, has2 = j + 2 < b.r, has3 = j + 3 < b.r;
        const T* c0 = bxnr + static_cast<int64_t>(j) * ld;
        T s[4][CB];
#pragma unroll
        for (int p = 0; p < 4; ++p) {
#pragma unroll
            for (int q = 0; q < CB; ++q) s[p][q] = T(0.0);
        }
        for (int i0 = 0; i0 < b.ntot; i0 += kChunkRows) {
            const int rows = min(kChunkRows, b.ntot - i0);
            __syncthreads();  // (the last chunk is read)
            for (int q = 0; q < nc; ++q) {
                const T* src = xn + q * ld + i0;
                T* dst = chunk + q * kChunkRows;
                for (int i = threadIdx.x; i < rows; i += blockDim.x) dst[i] = src[i];
            }
            __syncthreads();
            if (has0) {
                for (int i = lane; i < rows; i += 32) {
                    const int64_t gi = i0 + i;
                    const T a0 = c0[gi];
                    const T a1 = has1 ? c0[ld + gi] : T(0.0);
                    const T a2 = has2 ? c0[2 * ld + gi] : T(0.0);
                    const T a3 = has3 ? c0[3 * ld + gi] : T(0.0);
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc)) {
                            const T xi = chunk[q * kChunkRows + i];
                            s[0][q] += a0 * xi;
                            if (has1) s[1][q] += a1 * xi;
                            if (has2) s[2][q] += a2 * xi;
                            if (has3) s[3][q] += a3 * xi;
                        }
                    }
                }
            }
        }
        if (has0) {
#pragma unroll
            for (int p = 0; p < 4; ++p) {
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) s[p][q] = warp_sum(s[p][q]);
                }
            }
            if (lane == 0) {
#pragma unroll
                for (int p = 0; p < 4; ++p) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc) && j + p < b.r) {
                            T& x = xr[q * stride + j + p];
                            x = kSubtract ? x - s[p][q] : x + s[p][q];
                        }
                    }
                }
            }
        }
    }
}

// x_S and x_R of the pass's columns into shared memory (column q at sh + q * stride)
template<typename T>
__device__ __forceinline__ void load_sr(const SolveBox& b, const T* v0, T* sh, int stride, int nc) {
    for (int q = 0; q < nc; ++q) {
        const T* v = v0 + static_cast<int64_t>(q) * b.n;
        T* xs = sh + q * stride;
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = threadIdx.x; j < b.r; j += blockDim.x) xs[b.k + j] = v[b.red[j]];
    }
}
// x_R of the pass's columns back into the vector
template<typename T>
__device__ __forceinline__ void store_r(const SolveBox& b, T* v0, const T* sh, int stride, int nc) {
    for (int q = 0; q < nc; ++q) {
        T* v = v0 + static_cast<int64_t>(q) * b.n;
        const T* xr = sh + q * stride + b.k;
        for (int j = threadIdx.x; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
    }
}
// x_N of the pass's columns (column q at xn + q * xs), the neighbors' current
// rows in X_NR's row order
template<typename T>
__device__ __forceinline__ void gather_n(const SolveBox& b, const SolveSlot* slots, const T* vec, const T* ghost,
                                         int nrhs, int c0, T* xn, int64_t xs, int nc) {
    for (int q = 0; q < nc; ++q) {
        T* x = xn + q * xs;
        for (int a = 0; a < b.nslots; ++a) {
            const SolveSlot sl = slots[b.slot0 + a];
            const T* src = (sl.ghost ? ghost : vec) + sl.vec * nrhs + static_cast<int64_t>(c0 + q) * sl.n;
            for (int i = threadIdx.x; i < sl.count; i += blockDim.x) {
                x[sl.row0 + i] = sl.full ? src[i] : src[sl.skel[i]];
            }
        }
    }
}

template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) solve_forward_kernel(const SolveBox* boxes, const int* wave, T* vec,
                                                                 T* work, int nrhs, int cols) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    const int stride = b.k + b.r;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        load_sr(b, v0, sh, stride, nc);
        __syncthreads();
        if (bT != nullptr) {  // x_R -= T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                T s[CB];
                col_dot_cols<T, CB>(bT + static_cast<int64_t>(j) * b.k, sh, stride, nc, b.k, lane, s);
                if (lane == 0) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc)) sh[q * stride + b.k + j] -= s[q];
                    }
                }
            }
            __syncthreads();
        }
        store_r(b, v0, sh, stride, nc);
        if (bxsr != nullptr) {  // x_S += X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                accumulate_row<T, CB>(bxsr + i, b.k, sh + b.k, stride, nc, 0, b.r, s);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) v0[static_cast<int64_t>(q) * b.n + b.skel[i]] = sh[q * stride + i] + s[q];
                }
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates X_NR x_R
            T* u0 = work + b.work * nrhs + static_cast<int64_t>(c0) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                T s[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                accumulate_row<T, CB>(bxnr + i, b.ntot, sh + b.k, stride, nc, 0, b.r, s);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) u0[static_cast<int64_t>(q) * b.ntot + i] = s[q];
                }
            }
        }
        __syncthreads();
    }
}

// x_N in shared memory (after each column's x_S, x_R) when kSharedXn, else
// in the work buffer
template<typename T, int CB, bool kSharedXn>
__global__ void __launch_bounds__(kThreads) solve_backward_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                                  const int* wave, T* vec, const T* ghost, T* work,
                                                                  int nrhs, int cols) {
    // with x_N in the work buffer, its chunks first, then the columns
    T* const chunk = shared_array<T>();
    T* sh = kSharedXn ? chunk : chunk + kChunkRows * CB;
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    const int kr = b.k + b.r;
    const int stride = kSharedXn ? kr + b.ntot : kr;
    const int64_t xn_stride = kSharedXn ? stride : b.ntot;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        T* xn = kSharedXn ? sh + kr : work + b.work * nrhs + static_cast<int64_t>(c0) * b.ntot;
        load_sr(b, v0, sh, stride, nc);
        gather_n(b, slots, static_cast<const T*>(vec), ghost, nrhs, c0, xn, xn_stride, nc);
        __syncthreads();
        // x_R += X_SR^T x_S, then x_R += X_NR^T x_N (two products, as on the host);
        // a warp takes four columns at a time
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            T s[4][CB];
            if (bxsr != nullptr) {
                cols4_dot_cols<T, CB>(bxsr, b.k, sh, stride, nc, b.k, j, b.r, lane, s);
            } else {
#pragma unroll
                for (int p = 0; p < 4; ++p) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) s[p][q] = T(0.0);
                }
            }
            if (lane == 0) {
#pragma unroll
                for (int p = 0; p < 4; ++p) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc) && j + p < b.r) {
                            T& x = sh[q * stride + b.k + j + p];
                            x = x + s[p][q];
                        }
                    }
                }
            }
            if constexpr (kSharedXn) {
                cols4_dot_cols<T, CB>(bxnr, b.ntot, xn, xn_stride, nc, b.ntot, j, b.r, lane, s);
                if (lane == 0) {
#pragma unroll
                    for (int p = 0; p < 4; ++p) {
#pragma unroll
                        for (int q = 0; q < CB; ++q) {
                            if (live(q, nc) && j + p < b.r) {
                            T& x = sh[q * stride + b.k + j + p];
                            x = x + s[p][q];
                        }
                        }
                    }
                }
            }
        }
        if constexpr (!kSharedXn) xnr_t_chunked<T, CB, false>(b, bxnr, xn, chunk, sh + b.k, stride, nc);
        __syncthreads();
        store_r(b, v0, sh, stride, nc);
        if (bT != nullptr) {  // x_S -= T x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T d[CB];
                row_dot_cols<T, CB>(bT + i, b.k, sh + b.k, stride, nc, b.r, d);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) v0[static_cast<int64_t>(q) * b.n + b.skel[i]] = sh[q * stride + i] - d[q];
                }
            }
        }
        __syncthreads();
    }
}

// LAPACK getrs ('N') with the 1-based interchanges of getrf.
template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) solve_diagonal_kernel(const SolveBox* boxes, const int* list, T* vec,
                                                                  int nrhs, int cols) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[list[blockIdx.x]];
    const T* lu = static_cast<const T*>(b.lu);
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        for (int q = 0; q < nc; ++q) {
            const T* v = v0 + static_cast<int64_t>(q) * b.n;
            for (int j = tid; j < r; j += blockDim.x) sh[q * r + j] = v[b.red[j]];
        }
        __syncthreads();
        if (tid < nc) {  // a thread a column
            T* xr = sh + tid * r;
            for (int i = 0; i < r; ++i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const T t = xr[i];
                    xr[i] = xr[p];
                    xr[p] = t;
                }
            }
        }
        for (int i = 0; i < r; ++i) {  // unit lower L
            __syncthreads();
            T xi[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) xi[q] = live(q, nc) ? sh[q * r + i] : T(0.0);
            const T* col = lu + static_cast<int64_t>(i) * r;
            for (int l = i + 1 + tid; l < r; l += blockDim.x) {
                const T a = col[l];
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) sh[q * r + l] -= a * xi[q];
                }
            }
        }
        for (int i = r - 1; i >= 0; --i) {  // upper U
            __syncthreads();
            if (tid < nc) sh[tid * r + i] = sh[tid * r + i] / lu[i + static_cast<int64_t>(i) * r];
            __syncthreads();
            T xi[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) xi[q] = live(q, nc) ? sh[q * r + i] : T(0.0);
            const T* col = lu + static_cast<int64_t>(i) * r;
            for (int l = tid; l < i; l += blockDim.x) {
                const T a = col[l];
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) sh[q * r + l] -= a * xi[q];
                }
            }
        }
        __syncthreads();
        for (int q = 0; q < nc; ++q) {
            T* v = v0 + static_cast<int64_t>(q) * b.n;
            for (int j = tid; j < r; j += blockDim.x) v[b.red[j]] = sh[q * r + j];
        }
        __syncthreads();
    }
}

// ---- the multiply

template<typename T, int CB, bool kSharedXn>
__global__ void __launch_bounds__(kThreads) mul_forward_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                               const int* wave, T* vec, const T* ghost, T* work,
                                                               int nrhs, int cols) {
    // with x_N in the work buffer, its chunks first, then the columns
    T* const chunk = shared_array<T>();
    T* sh = kSharedXn ? chunk : chunk + kChunkRows * CB;
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    const int kr = b.k + b.r;
    const int stride = kSharedXn ? kr + b.ntot : kr;
    const int64_t xn_stride = kSharedXn ? stride : b.ntot;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        T* xn = kSharedXn ? sh + kr : work + b.work * nrhs + static_cast<int64_t>(c0) * b.ntot;
        load_sr(b, v0, sh, stride, nc);
        __syncthreads();
        if (bT != nullptr) {  // x_S += T x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                accumulate_row<T, CB>(bT + i, b.k, sh + b.k, stride, nc, 0, b.r, s);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) {
                        sh[q * stride + i] += s[q];
                        v0[static_cast<int64_t>(q) * b.n + b.skel[i]] = sh[q * stride + i];
                    }
                }
            }
        }
        gather_n(b, slots, static_cast<const T*>(vec), ghost, nrhs, c0, xn, xn_stride, nc);
        __syncthreads();
        // x_R -= X_SR^T x_S, then x_R -= X_NR^T x_N (two products, as on the host)
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            T s[4][CB];
            if (bxsr != nullptr) {
                cols4_dot_cols<T, CB>(bxsr, b.k, sh, stride, nc, b.k, j, b.r, lane, s);
            } else {
#pragma unroll
                for (int p = 0; p < 4; ++p) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) s[p][q] = T(0.0);
                }
            }
            if (lane == 0) {
#pragma unroll
                for (int p = 0; p < 4; ++p) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc) && j + p < b.r) {
                            T& x = sh[q * stride + b.k + j + p];
                            x = x - s[p][q];
                        }
                    }
                }
            }
            if constexpr (kSharedXn) {
                cols4_dot_cols<T, CB>(bxnr, b.ntot, xn, xn_stride, nc, b.ntot, j, b.r, lane, s);
                if (lane == 0) {
#pragma unroll
                    for (int p = 0; p < 4; ++p) {
#pragma unroll
                        for (int q = 0; q < CB; ++q) {
                            if (live(q, nc) && j + p < b.r) {
                            T& x = sh[q * stride + b.k + j + p];
                            x = x - s[p][q];
                        }
                        }
                    }
                }
            }
        }
        if constexpr (!kSharedXn) xnr_t_chunked<T, CB, true>(b, bxnr, xn, chunk, sh + b.k, stride, nc);
        __syncthreads();
        store_r(b, v0, sh, stride, nc);
        __syncthreads();
    }
}

// x_R = P L U x_R (LAPACK trmm U, trmm unit L, laswp with the interchanges reversed)
template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) mul_diagonal_kernel(const SolveBox* boxes, const int* list, T* vec,
                                                                int nrhs, int cols) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[list[blockIdx.x]];
    const T* lu = static_cast<const T*>(b.lu);
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    const int stride = 2 * r;  // a column's x, then its y
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        for (int q = 0; q < nc; ++q) {
            const T* v = v0 + static_cast<int64_t>(q) * b.n;
            for (int j = tid; j < r; j += blockDim.x) sh[q * stride + j] = v[b.red[j]];
        }
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // y = U x
            T s[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) s[q] = T(0.0);
            accumulate_row<T, CB>(lu + i, r, sh, stride, nc, i, r, s);
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                if (live(q, nc)) sh[q * stride + r + i] = s[q];
            }
        }
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // x = L y (unit diagonal)
            T s[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) s[q] = live(q, nc) ? sh[q * stride + r + i] : T(0.0);
            accumulate_row<T, CB>(lu + i, r, sh + r, stride, nc, 0, i, s);
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                if (live(q, nc)) sh[q * stride + i] = s[q];
            }
        }
        __syncthreads();
        if (tid < nc) {  // a thread a column
            T* x = sh + tid * stride;
            for (int i = r - 1; i >= 0; --i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const T t = x[i];
                    x[i] = x[p];
                    x[p] = t;
                }
            }
        }
        __syncthreads();
        for (int q = 0; q < nc; ++q) {
            T* v = v0 + static_cast<int64_t>(q) * b.n;
            for (int j = tid; j < r; j += blockDim.x) v[b.red[j]] = sh[q * stride + j];
        }
        __syncthreads();
    }
}

template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) mul_backward_kernel(const SolveBox* boxes, const int* wave, T* vec,
                                                                T* work, int nrhs, int cols) {
    T* sh = shared_array<T>();
    const SolveBox b = boxes[wave[blockIdx.x]];
    const T* bT = static_cast<const T*>(b.T);
    const T* bxsr = static_cast<const T*>(b.xsr);
    const T* bxnr = static_cast<const T*>(b.xnr);
    const int stride = b.k + b.r;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + b.vec * nrhs + static_cast<int64_t>(c0) * b.n;
        load_sr(b, v0, sh, stride, nc);
        __syncthreads();
        if (bxsr != nullptr) {  // x_S -= X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                T s[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                accumulate_row<T, CB>(bxsr + i, b.k, sh + b.k, stride, nc, 0, b.r, s);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) {
                        sh[q * stride + i] -= s[q];
                        v0[static_cast<int64_t>(q) * b.n + b.skel[i]] = sh[q * stride + i];
                    }
                }
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates -X_NR x_R
            T* u0 = work + b.work * nrhs + static_cast<int64_t>(c0) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                T s[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                accumulate_row<T, CB>(bxnr + i, b.ntot, sh + b.k, stride, nc, 0, b.r, s);
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) u0[static_cast<int64_t>(q) * b.ntot + i] = -s[q];
                }
            }
        }
        __syncthreads();
        if (bT != nullptr) {  // x_R += T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                T s[CB];
                col_dot_cols<T, CB>(bT + static_cast<int64_t>(j) * b.k, sh, stride, nc, b.k, lane, s);
                if (lane == 0) {
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc)) v0[static_cast<int64_t>(q) * b.n + b.red[j]] = sh[q * stride + b.k + j] + s[q];
                    }
                }
            }
        }
        __syncthreads();
    }
}

template<typename T>
__global__ void solve_accum_kernel(const SolveAccum* items, const SolvePart* parts, int nrhs, T* vec, const T* work,
                                   T* outbox, const T* inbox) {
    const SolveAccum it = items[blockIdx.x];
    const SolvePart* p = parts + it.part0;
    const int64_t total = static_cast<int64_t>(it.count) * nrhs;
    T* dst = (it.outbox ? outbox : vec) + it.dst * nrhs;
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        const int i = static_cast<int>(e % it.count);
        const int64_t c = e / it.count;
        T s = T(0.0);
        for (int q = 0; q < it.nparts; ++q) {
            const SolvePart pq = p[q];
            const T* src = (pq.inbox ? inbox : work) + pq.base * nrhs + pq.row0;
            const T x = src[i + c * pq.ld];
            s = q == 0 ? x : s + x;
        }
        const int64_t row = it.rows != nullptr ? it.rows[i] : i;
        T* d = dst + row + c * it.ld;
        if (it.outbox) {
            *d = s;
        } else {
            *d += s;
        }
    }
}

template<typename T>
__global__ void solve_copy_kernel(const SolveCopy* items, int nrhs, const T* src, T* dst) {
    const SolveCopy it = items[blockIdx.x];
    const int64_t count = static_cast<int64_t>(it.n) * nrhs;
    const T* s = src + it.src * nrhs;
    T* d = dst + it.dst * nrhs;
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < count;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        d[e] = s[e];
    }
}

__global__ void copy_spans_kernel(const CopySpan* items) {
    const CopySpan it = items[blockIdx.x];
    const uint32_t* s = static_cast<const uint32_t*>(it.src);
    uint32_t* d = static_cast<uint32_t*>(it.dst);
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < it.words;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        d[e] = s[e];
    }
}

// the per-box kernels may use up to kSolveSharedLimit bytes
template<typename Kernel>
void allow_shared(Kernel kernel, size_t bytes) {
    if (bytes > kSolveSharedLimit) throw std::runtime_error("device solve: box too large for shared memory");
    if (bytes > 48 * 1024) {
        cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(kSolveSharedLimit));
    }
}

unsigned grid_y(int64_t max_count) {
    const int64_t blocks = (max_count + kThreads - 1) / kThreads;
    return static_cast<unsigned>(std::max<int64_t>(1, std::min<int64_t>(blocks, 64)));
}

// The columns a pass of a per-box kernel takes when each needs `per_column`
// bytes of shared memory: one, or up to kMaxCols<T> within half the limit
// (so two blocks still fit on an SM).
template<typename T>
int pass_columns(int nrhs, size_t per_column) {
    if (nrhs <= 1) return 1;
    const size_t fit = (kSolveSharedLimit / 2) / std::max<size_t>(per_column, 1);
    const size_t cols = std::min<size_t>({static_cast<size_t>(nrhs), static_cast<size_t>(kMaxCols<T>), fit});
    return static_cast<int>(std::max<size_t>(cols, 1));
}
// f(std::integral_constant<int, CB>) with the smallest instantiated CB (1, 2,
// 4 or 8) that holds `cols` columns
template<typename T, typename F>
void with_cols(int cols, F&& f) {
    if (cols <= 1) {
        f(std::integral_constant<int, 1>{});
    } else if (cols <= 2) {
        f(std::integral_constant<int, 2>{});
    } else if (cols <= 4 || kMaxCols<T> < 8) {
        f(std::integral_constant<int, 4>{});
    } else {
        f(std::integral_constant<int, (kMaxCols<T> < 8 ? 4 : 8)>{});
    }
}

}  // namespace

template<typename T>
void launch_solve_forward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                          cudaStream_t stream) {
    if (count <= 0) return;
    const size_t per_column = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(solve_forward_one_kernel<T>, shared);
        solve_forward_one_kernel<T><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(solve_forward_kernel<T, CB>, shared);
            solve_forward_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs, cols);
        });
    }
    check_launch("solve_forward_kernel");
}

template<typename T>
void launch_solve_backward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                           const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    // x_N in shared memory when that takes as many columns a pass as the work
    // buffer would (one column: whenever it fits, as before)
    const size_t per_column = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const size_t per_column_xn = per_column + static_cast<size_t>(max_ntot) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const int cols_xn = per_column_xn <= kSolveSharedLimit ? pass_columns<T>(nrhs, per_column_xn) : 0;
    if (cols_xn >= cols && cols_xn == 1) {
        allow_shared(solve_backward_one_kernel<T, true>, per_column_xn);
        solve_backward_one_kernel<T, true><<<count, kThreads, per_column_xn, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                       work, nrhs);
    } else {
        with_cols<T>(cols_xn >= cols ? cols_xn : cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            if (cols_xn >= cols) {
                const size_t shared = static_cast<size_t>(cols_xn) * per_column_xn;
                allow_shared(solve_backward_kernel<T, CB, true>, shared);
                solve_backward_kernel<T, CB, true><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                        work, nrhs, cols_xn);
            } else {  // (x_N in the work buffer, through shared memory in chunks)
                const size_t shared =
                    static_cast<size_t>(kChunkRows) * CB * sizeof(T) + static_cast<size_t>(cols) * per_column;
                allow_shared(solve_backward_kernel<T, CB, false>, shared);
                solve_backward_kernel<T, CB, false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                         work, nrhs, cols);
            }
        });
    }
    check_launch("solve_backward_kernel");
}

template<typename T>
void launch_solve_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                           cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t per_column = static_cast<size_t>(max_r) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(solve_diagonal_one_kernel<T>, shared);
        solve_diagonal_one_kernel<T><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(solve_diagonal_kernel<T, CB>, shared);
            solve_diagonal_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs, cols);
        });
    }
    check_launch("solve_diagonal_kernel");
}

template<typename T>
void launch_solve_accum(const SolveAccum* items, const SolvePart* parts, int count, int max_count, int nrhs, T* vec,
                        const T* work, T* outbox, const T* inbox, cudaStream_t stream) {
    if (count <= 0 || max_count <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(static_cast<int64_t>(max_count) * nrhs));
    solve_accum_kernel<T><<<grid, kThreads, 0, stream>>>(items, parts, nrhs, vec, work, outbox, inbox);
    check_launch("solve_accum_kernel");
}

void launch_copy_spans(const CopySpan* items, int count, int64_t max_words, cudaStream_t stream) {
    if (count <= 0 || max_words <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(max_words));
    copy_spans_kernel<<<grid, kThreads, 0, stream>>>(items);
    check_launch("copy_spans_kernel");
}

__global__ void patch_xnr_kernel(SolveBox* boxes, const int* idx, void* const* dsts, void* const* t_dsts,
                                 void* const* xsr_dsts, int count) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    SolveBox& b = boxes[idx[i]];
    if (dsts[i] != nullptr) b.xnr = dsts[i];
    if (t_dsts[i] != nullptr) b.T = t_dsts[i];
    if (xsr_dsts[i] != nullptr) b.xsr = xsr_dsts[i];
}

void launch_patch_xnr(SolveBox* boxes, const int* idx, void* const* dsts, void* const* t_dsts,
                      void* const* xsr_dsts, int count, cudaStream_t stream) {
    if (count <= 0) return;
    const int threads = 128;
    patch_xnr_kernel<<<(count + threads - 1) / threads, threads, 0, stream>>>(boxes, idx, dsts, t_dsts, xsr_dsts,
                                                                              count);
    check_launch("patch_xnr_kernel");
}

// One block per item: the threads stride the diagonal, then a fixed-order
// reduction (the sums do not depend on the schedule).
template<typename T>
__global__ void logdet_lu_kernel(const LogdetItem* items, double* out, int* bad) {
    constexpr double kPi = 3.14159265358979323846;
    __shared__ double s_logabs[kThreads], s_arg[kThreads];
    __shared__ int s_bad[kThreads];
    const LogdetItem item = items[blockIdx.x];
    const T* lu = static_cast<const T*>(item.lu);
    const int r = item.r;
    double logabs = 0.0, arg = 0.0;
    int first_bad = 0;
    for (int k = threadIdx.x; k < r; k += blockDim.x) {
        const T d = lu[k + static_cast<int64_t>(k) * r];
        double magnitude, phase;
        if constexpr (std::is_same_v<T, double>) {
            magnitude = fabs(d);
            phase = d < 0.0 ? kPi : 0.0;  // std::arg of a real
        } else {
            magnitude = hypot(d.re, d.im);
            phase = atan2(d.im, d.re);
        }
        const int pivot = item.ipiv[k];
        if (!(magnitude > 1e-14) || pivot < 1 || pivot > r) {
            if (first_bad == 0) first_bad = k + 1;
        }
        logabs += log(magnitude);
        arg += phase;
        if (pivot != k + 1) arg += kPi;
    }
    s_logabs[threadIdx.x] = logabs;
    s_arg[threadIdx.x] = arg;
    s_bad[threadIdx.x] = first_bad;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            s_logabs[threadIdx.x] += s_logabs[threadIdx.x + stride];
            s_arg[threadIdx.x] += s_arg[threadIdx.x + stride];
            const int a = s_bad[threadIdx.x], b = s_bad[threadIdx.x + stride];
            s_bad[threadIdx.x] = a == 0 ? b : (b == 0 ? a : min(a, b));
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        out[2 * blockIdx.x] = s_logabs[0];
        out[2 * blockIdx.x + 1] = s_arg[0];
        bad[blockIdx.x] = s_bad[0];
    }
}

template<typename T>
void launch_logdet_lu(const LogdetItem* items, int count, double* out, int* bad, cudaStream_t stream) {
    if (count <= 0) return;
    logdet_lu_kernel<T><<<count, kThreads, 0, stream>>>(items, out, bad);
    check_launch("logdet_lu_kernel");
}

template<typename T>
void launch_mul_forward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                        const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    // x_N placed as in launch_solve_backward
    const size_t per_column = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const size_t per_column_xn = per_column + static_cast<size_t>(max_ntot) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const int cols_xn = per_column_xn <= kSolveSharedLimit ? pass_columns<T>(nrhs, per_column_xn) : 0;
    if (cols_xn >= cols && cols_xn == 1) {
        allow_shared(mul_forward_one_kernel<T, true>, per_column_xn);
        mul_forward_one_kernel<T, true><<<count, kThreads, per_column_xn, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                    work, nrhs);
    } else {
        with_cols<T>(cols_xn >= cols ? cols_xn : cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            if (cols_xn >= cols) {
                const size_t shared = static_cast<size_t>(cols_xn) * per_column_xn;
                allow_shared(mul_forward_kernel<T, CB, true>, shared);
                mul_forward_kernel<T, CB, true><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                     work, nrhs, cols_xn);
            } else {  // (x_N in the work buffer, through shared memory in chunks)
                const size_t shared =
                    static_cast<size_t>(kChunkRows) * CB * sizeof(T) + static_cast<size_t>(cols) * per_column;
                allow_shared(mul_forward_kernel<T, CB, false>, shared);
                mul_forward_kernel<T, CB, false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost,
                                                                                      work, nrhs, cols);
            }
        });
    }
    check_launch("mul_forward_kernel");
}

template<typename T>
void launch_mul_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                         cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t per_column = 2 * static_cast<size_t>(max_r) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(mul_diagonal_one_kernel<T>, shared);
        mul_diagonal_one_kernel<T><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(mul_diagonal_kernel<T, CB>, shared);
            mul_diagonal_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs, cols);
        });
    }
    check_launch("mul_diagonal_kernel");
}

template<typename T>
void launch_mul_backward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                         cudaStream_t stream) {
    if (count <= 0) return;
    const size_t per_column = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(mul_backward_one_kernel<T>, shared);
        mul_backward_one_kernel<T><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(mul_backward_kernel<T, CB>, shared);
            mul_backward_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs, cols);
        });
    }
    check_launch("mul_backward_kernel");
}

template<typename T>
void launch_solve_copy(const SolveCopy* items, int count, int max_n, int nrhs, const T* src, T* dst,
                       cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(static_cast<int64_t>(max_n) * nrhs));
    solve_copy_kernel<T><<<grid, kThreads, 0, stream>>>(items, nrhs, src, dst);
    check_launch("solve_copy_kernel");
}

#define H2_SOLVE_LAUNCHES(T)                                                                                          \
    template void launch_solve_forward<T>(const SolveBox*, const int*, int, T*, T*, int, int, cudaStream_t);          \
    template void launch_solve_backward<T>(const SolveBox*, const SolveSlot*, const int*, int, T*, const T*, T*, int, \
                                           int, int, cudaStream_t);                                                   \
    template void launch_solve_diagonal<T>(const SolveBox*, const int*, int, T*, int, int, cudaStream_t);             \
    template void launch_solve_accum<T>(const SolveAccum*, const SolvePart*, int, int, int, T*, const T*, T*,         \
                                        const T*, cudaStream_t);                                                      \
    template void launch_solve_copy<T>(const SolveCopy*, int, int, int, const T*, T*, cudaStream_t);                  \
    template void launch_mul_forward<T>(const SolveBox*, const SolveSlot*, const int*, int, T*, const T*, T*, int,    \
                                        int, int, cudaStream_t);                                                      \
    template void launch_mul_diagonal<T>(const SolveBox*, const int*, int, T*, int, int, cudaStream_t);               \
    template void launch_mul_backward<T>(const SolveBox*, const int*, int, T*, T*, int, int, cudaStream_t);        \
    template void launch_logdet_lu<T>(const LogdetItem*, int, double*, int*, cudaStream_t);
H2_SOLVE_LAUNCHES(double)
H2_SOLVE_LAUNCHES(dcomplex)
#undef H2_SOLVE_LAUNCHES

}  // namespace gpu
}  // namespace fmm
