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

// sum_j a[j * lda] x[j], four independent partial sums (loads in flight)
__device__ __forceinline__ double row_dot(const double* a, int lda, const double* x, int n) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
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

__global__ void __launch_bounds__(kThreads) solve_forward_kernel(const SolveBox* boxes, const int* wave,
                                                                 double* vec, double* work, int nrhs) {
    extern __shared__ double sh[];
    const SolveBox b = boxes[wave[blockIdx.x]];
    double* xs = sh;
    double* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (b.T != nullptr) {  // x_R -= T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                const double* col = b.T + static_cast<int64_t>(j) * b.k;
                double s = 0.0;
                for (int i = lane; i < b.k; i += 32) s += col[i] * xs[i];
                s = warp_sum(s);
                if (lane == 0) xr[j] -= s;
            }
            __syncthreads();
        }
        for (int j = tid; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
        if (b.xsr != nullptr) {  // x_S += X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                double s = 0.0;
                for (int j = 0; j < b.r; ++j) s += b.xsr[i + static_cast<int64_t>(j) * b.k] * xr[j];
                v[b.skel[i]] = xs[i] + s;
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates X_NR x_R
            double* u = work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                double s = 0.0;
                for (int j = 0; j < b.r; ++j) s += b.xnr[i + static_cast<int64_t>(j) * b.ntot] * xr[j];
                u[i] = s;
            }
        }
        __syncthreads();
    }
}

// Columns j .. j+3 of A^T x (A: rows x cols, ld): lane partial sums.
__device__ __forceinline__ void cols4_dot(const double* a, int64_t ld, const double* x, int rows, int j, int cols,
                                          int lane, double out[4]) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    const double* c0 = a + j * ld;
    const bool has1 = j + 1 < cols, has2 = j + 2 < cols, has3 = j + 3 < cols;
    for (int i = lane; i < rows; i += 32) {
        const double xi = x[i];
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
template<bool kSharedXn>
__global__ void __launch_bounds__(kThreads) solve_backward_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                                  const int* wave, double* vec, const double* ghost,
                                                                  double* work, int nrhs) {
    extern __shared__ double sh[];
    const SolveBox b = boxes[wave[blockIdx.x]];
    double* xs = sh;
    double* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        double* xn = kSharedXn ? sh + b.k + b.r : work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        // x_N: the neighbors' current rows, in X_NR's row order
        for (int a = 0; a < b.nslots; ++a) {
            const SolveSlot sl = slots[b.slot0 + a];
            const double* src = (sl.ghost ? ghost : vec) + sl.vec * nrhs + static_cast<int64_t>(c) * sl.n;
            for (int i = tid; i < sl.count; i += blockDim.x) xn[sl.row0 + i] = sl.full ? src[i] : src[sl.skel[i]];
        }
        __syncthreads();
        // x_R += X_SR^T x_S, then x_R += X_NR^T x_N (two products, as on the host);
        // a warp takes four columns at a time
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            double s1[4] = {0.0, 0.0, 0.0, 0.0}, s2[4];
            if (b.xsr != nullptr) cols4_dot(b.xsr, b.k, xs, b.k, j, b.r, lane, s1);
            cols4_dot(b.xnr, b.ntot, xn, b.ntot, j, b.r, lane, s2);
            if (lane == 0) {
                for (int q = 0; q < 4 && j + q < b.r; ++q) xr[j + q] = (xr[j + q] + s1[q]) + s2[q];
            }
        }
        __syncthreads();
        for (int j = tid; j < b.r; j += blockDim.x) v[b.red[j]] = xr[j];
        if (b.T != nullptr) {  // x_S -= T x_R
            for (int i = tid; i < b.k; i += blockDim.x) v[b.skel[i]] = xs[i] - row_dot(b.T + i, b.k, xr, b.r);
        }
        __syncthreads();
    }
}

// LAPACK getrs ('N') with the 1-based interchanges of getrf.
__global__ void __launch_bounds__(kThreads) solve_diagonal_kernel(const SolveBox* boxes, const int* list, double* vec,
                                                                  int nrhs) {
    extern __shared__ double xr[];
    const SolveBox b = boxes[list[blockIdx.x]];
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int j = tid; j < r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (tid == 0) {
            for (int i = 0; i < r; ++i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const double t = xr[i];
                    xr[i] = xr[p];
                    xr[p] = t;
                }
            }
        }
        for (int i = 0; i < r; ++i) {  // unit lower L
            __syncthreads();
            const double xi = xr[i];
            const double* col = b.lu + static_cast<int64_t>(i) * r;
            for (int l = i + 1 + tid; l < r; l += blockDim.x) xr[l] -= col[l] * xi;
        }
        for (int i = r - 1; i >= 0; --i) {  // upper U
            __syncthreads();
            if (tid == 0) xr[i] /= b.lu[i + static_cast<int64_t>(i) * r];
            __syncthreads();
            const double xi = xr[i];
            const double* col = b.lu + static_cast<int64_t>(i) * r;
            for (int l = tid; l < i; l += blockDim.x) xr[l] -= col[l] * xi;
        }
        __syncthreads();
        for (int j = tid; j < r; j += blockDim.x) v[b.red[j]] = xr[j];
        __syncthreads();
    }
}

// ---- the multiply

template<bool kSharedXn>
__global__ void __launch_bounds__(kThreads) mul_forward_kernel(const SolveBox* boxes, const SolveSlot* slots,
                                                               const int* wave, double* vec, const double* ghost,
                                                               double* work, int nrhs) {
    extern __shared__ double sh[];
    const SolveBox b = boxes[wave[blockIdx.x]];
    double* xs = sh;
    double* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        double* xn = kSharedXn ? sh + b.k + b.r : work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (b.T != nullptr) {  // x_S += T x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                double s = 0.0;
                for (int j = 0; j < b.r; ++j) s += b.T[i + static_cast<int64_t>(j) * b.k] * xr[j];
                xs[i] += s;
                v[b.skel[i]] = xs[i];
            }
        }
        for (int a = 0; a < b.nslots; ++a) {  // x_N, in X_NR's row order
            const SolveSlot sl = slots[b.slot0 + a];
            const double* src = (sl.ghost ? ghost : vec) + sl.vec * nrhs + static_cast<int64_t>(c) * sl.n;
            for (int i = tid; i < sl.count; i += blockDim.x) xn[sl.row0 + i] = sl.full ? src[i] : src[sl.skel[i]];
        }
        __syncthreads();
        // x_R -= X_SR^T x_S, then x_R -= X_NR^T x_N (two products, as on the host)
        for (int j = 4 * warp; j < b.r; j += 4 * nwarps) {
            double s1[4] = {0.0, 0.0, 0.0, 0.0}, s2[4];
            if (b.xsr != nullptr) cols4_dot(b.xsr, b.k, xs, b.k, j, b.r, lane, s1);
            cols4_dot(b.xnr, b.ntot, xn, b.ntot, j, b.r, lane, s2);
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
__global__ void __launch_bounds__(kThreads) mul_diagonal_kernel(const SolveBox* boxes, const int* list, double* vec,
                                                                int nrhs) {
    extern __shared__ double sh[];
    const SolveBox b = boxes[list[blockIdx.x]];
    const int r = b.r, tid = threadIdx.x;
    if (r == 0) return;
    double* x = sh;
    double* y = sh + r;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int j = tid; j < r; j += blockDim.x) x[j] = v[b.red[j]];
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // y = U x
            double s = 0.0;
            for (int j = i; j < r; ++j) s += b.lu[i + static_cast<int64_t>(j) * r] * x[j];
            y[i] = s;
        }
        __syncthreads();
        for (int i = tid; i < r; i += blockDim.x) {  // x = L y (unit diagonal)
            double s = y[i];
            for (int j = 0; j < i; ++j) s += b.lu[i + static_cast<int64_t>(j) * r] * y[j];
            x[i] = s;
        }
        __syncthreads();
        if (tid == 0) {
            for (int i = r - 1; i >= 0; --i) {
                const int p = b.ipiv[i] - 1;
                if (p != i) {
                    const double t = x[i];
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

__global__ void __launch_bounds__(kThreads) mul_backward_kernel(const SolveBox* boxes, const int* wave, double* vec,
                                                                double* work, int nrhs) {
    extern __shared__ double sh[];
    const SolveBox b = boxes[wave[blockIdx.x]];
    double* xs = sh;
    double* xr = sh + b.k;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        double* v = vec + b.vec * nrhs + static_cast<int64_t>(c) * b.n;
        for (int i = tid; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        for (int j = tid; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        if (b.xsr != nullptr) {  // x_S -= X_SR x_R
            for (int i = tid; i < b.k; i += blockDim.x) {
                double s = 0.0;
                for (int j = 0; j < b.r; ++j) s += b.xsr[i + static_cast<int64_t>(j) * b.k] * xr[j];
                xs[i] -= s;
                v[b.skel[i]] = xs[i];
            }
        }
        if (b.ntot > 0) {  // the neighbors' updates -X_NR x_R
            double* u = work + b.work * nrhs + static_cast<int64_t>(c) * b.ntot;
            for (int i = tid; i < b.ntot; i += blockDim.x) {
                double s = 0.0;
                for (int j = 0; j < b.r; ++j) s += b.xnr[i + static_cast<int64_t>(j) * b.ntot] * xr[j];
                u[i] = -s;
            }
        }
        __syncthreads();
        if (b.T != nullptr) {  // x_R += T^T x_S
            for (int j = warp; j < b.r; j += nwarps) {
                const double* col = b.T + static_cast<int64_t>(j) * b.k;
                double s = 0.0;
                for (int i = lane; i < b.k; i += 32) s += col[i] * xs[i];
                s = warp_sum(s);
                if (lane == 0) v[b.red[j]] = xr[j] + s;
            }
        }
        __syncthreads();
    }
}

__global__ void solve_accum_kernel(const SolveAccum* items, const SolvePart* parts, int nrhs, double* vec,
                                   const double* work, double* outbox, const double* inbox) {
    const SolveAccum it = items[blockIdx.x];
    const SolvePart* p = parts + it.part0;
    const int64_t total = static_cast<int64_t>(it.count) * nrhs;
    double* dst = (it.outbox ? outbox : vec) + it.dst * nrhs;
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        const int i = static_cast<int>(e % it.count);
        const int64_t c = e / it.count;
        double s = 0.0;
        for (int q = 0; q < it.nparts; ++q) {
            const SolvePart pq = p[q];
            const double* src = (pq.inbox ? inbox : work) + pq.base * nrhs + pq.row0;
            const double x = src[i + c * pq.ld];
            s = q == 0 ? x : s + x;
        }
        const int64_t row = it.rows != nullptr ? it.rows[i] : i;
        double* d = dst + row + c * it.ld;
        if (it.outbox) {
            *d = s;
        } else {
            *d += s;
        }
    }
}

__global__ void solve_copy_kernel(const SolveCopy* items, int nrhs, const double* src, double* dst) {
    const SolveCopy it = items[blockIdx.x];
    const int64_t count = static_cast<int64_t>(it.n) * nrhs;
    const double* s = src + it.src * nrhs;
    double* d = dst + it.dst * nrhs;
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

void launch_solve_forward(const SolveBox* boxes, const int* wave, int count, double* vec, double* work, int nrhs,
                          int max_n, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(double);
    allow_shared(solve_forward_kernel, shared);
    solve_forward_kernel<<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
    check_launch("solve_forward_kernel");
}

void launch_solve_backward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, double* vec,
                           const double* ghost, double* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(double);
    const size_t with_xn = shared + static_cast<size_t>(max_ntot) * sizeof(double);
    if (with_xn <= kSolveSharedLimit) {
        allow_shared(solve_backward_kernel<true>, with_xn);
        solve_backward_kernel<true><<<count, kThreads, with_xn, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    } else {
        allow_shared(solve_backward_kernel<false>, shared);
        solve_backward_kernel<false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    }
    check_launch("solve_backward_kernel");
}

void launch_solve_diagonal(const SolveBox* boxes, const int* list, int count, double* vec, int nrhs, int max_r,
                           cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t shared = static_cast<size_t>(max_r) * sizeof(double);
    allow_shared(solve_diagonal_kernel, shared);
    solve_diagonal_kernel<<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
    check_launch("solve_diagonal_kernel");
}

void launch_solve_accum(const SolveAccum* items, const SolvePart* parts, int count, int max_count, int nrhs,
                        double* vec, const double* work, double* outbox, const double* inbox, cudaStream_t stream) {
    if (count <= 0 || max_count <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(static_cast<int64_t>(max_count) * nrhs));
    solve_accum_kernel<<<grid, kThreads, 0, stream>>>(items, parts, nrhs, vec, work, outbox, inbox);
    check_launch("solve_accum_kernel");
}

void launch_copy_spans(const CopySpan* items, int count, int64_t max_words, cudaStream_t stream) {
    if (count <= 0 || max_words <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(max_words));
    copy_spans_kernel<<<grid, kThreads, 0, stream>>>(items);
    check_launch("copy_spans_kernel");
}

void launch_mul_forward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, double* vec,
                        const double* ghost, double* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(double);
    const size_t with_xn = shared + static_cast<size_t>(max_ntot) * sizeof(double);
    if (with_xn <= kSolveSharedLimit) {
        allow_shared(mul_forward_kernel<true>, with_xn);
        mul_forward_kernel<true><<<count, kThreads, with_xn, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    } else {
        allow_shared(mul_forward_kernel<false>, shared);
        mul_forward_kernel<false><<<count, kThreads, shared, stream>>>(boxes, slots, wave, vec, ghost, work, nrhs);
    }
    check_launch("mul_forward_kernel");
}

void launch_mul_diagonal(const SolveBox* boxes, const int* list, int count, double* vec, int nrhs, int max_r,
                         cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t shared = 2 * static_cast<size_t>(max_r) * sizeof(double);
    allow_shared(mul_diagonal_kernel, shared);
    mul_diagonal_kernel<<<count, kThreads, shared, stream>>>(boxes, list, vec, nrhs);
    check_launch("mul_diagonal_kernel");
}

void launch_mul_backward(const SolveBox* boxes, const int* wave, int count, double* vec, double* work, int nrhs,
                         int max_n, cudaStream_t stream) {
    if (count <= 0) return;
    const size_t shared = static_cast<size_t>(std::max(max_n, 1)) * sizeof(double);
    allow_shared(mul_backward_kernel, shared);
    mul_backward_kernel<<<count, kThreads, shared, stream>>>(boxes, wave, vec, work, nrhs);
    check_launch("mul_backward_kernel");
}

void launch_solve_copy(const SolveCopy* items, int count, int max_n, int nrhs, const double* src, double* dst,
                       cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), grid_y(static_cast<int64_t>(max_n) * nrhs));
    solve_copy_kernel<<<grid, kThreads, 0, stream>>>(items, nrhs, src, dst);
    check_launch("solve_copy_kernel");
}

}  // namespace gpu
}  // namespace fmm
