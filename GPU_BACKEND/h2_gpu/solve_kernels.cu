// Kernels of the device solve.  See solve_kernels.hpp.

#include "solve_kernels.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

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

__device__ __forceinline__ double warp_sum(double v) {
    for (int offset = 16; offset > 0; offset >>= 1) v += __shfl_down_sync(0xffffffffu, v, offset);
    return v;
}
__device__ __forceinline__ dcomplex warp_sum(dcomplex v) { return dcomplex(warp_sum(v.re), warp_sum(v.im)); }

template<typename T>
__global__ void __launch_bounds__(kThreads) solve_forward_kernel(const SolveBox* boxes, const int* wave, T* vec,
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
__global__ void __launch_bounds__(kThreads) solve_backward_kernel(const SolveBox* boxes, const SolveSlot* slots,
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
__global__ void __launch_bounds__(kThreads) solve_diagonal_kernel(const SolveBox* boxes, const int* list, T* vec,
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
__global__ void __launch_bounds__(kThreads) mul_forward_kernel(const SolveBox* boxes, const SolveSlot* slots,
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
__global__ void __launch_bounds__(kThreads) mul_diagonal_kernel(const SolveBox* boxes, const int* list, T* vec,
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
__global__ void __launch_bounds__(kThreads) mul_backward_kernel(const SolveBox* boxes, const int* wave, T* vec,
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

}  // namespace

template<typename T>
void launch_solve_forward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                          cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    allow_shared(solve_forward_kernel<T>, shared);
    solve_forward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
    check_launch("solve_forward_kernel");
}

template<typename T>
void launch_solve_backward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                           const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const size_t with_xn = shared + static_cast<size_t>(max_ntot) * sizeof(T);
    if (with_xn <= kSolveSharedLimit) {
        allow_shared(solve_backward_kernel<T, true>, with_xn);
        solve_backward_kernel<T, true><<<count, kThreads, with_xn, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    } else {
        allow_shared(solve_backward_kernel<T, false>, shared);
        solve_backward_kernel<T, false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    }
    check_launch("solve_backward_kernel");
}

template<typename T>
void launch_solve_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                           cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t shared = static_cast<size_t>(max_r) * sizeof(T);
    allow_shared(solve_diagonal_kernel<T>, shared);
    solve_diagonal_kernel<T><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
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

template<typename T>
void launch_mul_forward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                        const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    const size_t with_xn = shared + static_cast<size_t>(max_ntot) * sizeof(T);
    if (with_xn <= kSolveSharedLimit) {
        allow_shared(mul_forward_kernel<T, true>, with_xn);
        mul_forward_kernel<T, true><<<count, kThreads, with_xn, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    } else {
        allow_shared(mul_forward_kernel<T, false>, shared);
        mul_forward_kernel<T, false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    }
    check_launch("mul_forward_kernel");
}

template<typename T>
void launch_mul_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                         cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t shared = 2 * static_cast<size_t>(max_r) * sizeof(T);
    allow_shared(mul_diagonal_kernel<T>, shared);
    mul_diagonal_kernel<T><<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
    check_launch("mul_diagonal_kernel");
}

template<typename T>
void launch_mul_backward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                         cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(T);
    allow_shared(mul_backward_kernel<T>, shared);
    mul_backward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
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
    template void launch_mul_backward<T>(const SolveBox*, const int*, int, T*, T*, int, int, cudaStream_t);
H2_SOLVE_LAUNCHES(double)
H2_SOLVE_LAUNCHES(dcomplex)
#undef H2_SOLVE_LAUNCHES

}  // namespace gpu
}  // namespace fmm
