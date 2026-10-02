// Kernels of the device matvec of the compression-only H2.  See
// h2_matvec_kernels.hpp.

#include "h2_matvec_kernels.hpp"

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

template<typename T>
__device__ __forceinline__ T* shared_array() {
    extern __shared__ __align__(16) unsigned char h2_shared[];
    return reinterpret_cast<T*>(h2_shared);
}

__device__ __forceinline__ double warp_sum(double v) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xffffffffu, v, o);
    return v;
}
__device__ __forceinline__ dcomplex warp_sum(dcomplex v) { return dcomplex(warp_sum(v.re), warp_sum(v.im)); }

// x_S += T x_R, x_R of a column in shared memory
template<typename T>
__global__ void __launch_bounds__(kThreads) h2_upward_kernel(const H2Box* boxes, T* vec, int64_t ld, int nrhs) {
    const H2Box b = boxes[blockIdx.x];
    if (b.r == 0 || b.T == nullptr) return;
    T* xr = shared_array<T>();
    const T* Tm = static_cast<const T*>(b.T);
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + static_cast<int64_t>(c) * ld + b.vec;
        for (int j = threadIdx.x; j < b.r; j += blockDim.x) xr[j] = v[b.red[j]];
        __syncthreads();
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) {
            T s = T(0.0);
            for (int j = 0; j < b.r; ++j) s += Tm[i + static_cast<int64_t>(j) * b.k] * xr[j];
            v[b.skel[i]] += s;
        }
        __syncthreads();
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) h2_gather_q_kernel(const H2Box* boxes, const T* vec, int64_t ld, T* q,
                                                               int64_t ldq, int nrhs) {
    const H2Box b = boxes[blockIdx.x];
    for (int c = 0; c < nrhs; ++c) {
        const T* v = vec + static_cast<int64_t>(c) * ld + b.vec;
        T* o = q + static_cast<int64_t>(c) * ldq + b.q;
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) o[i] = v[b.skel[i]];
    }
}

// y[rows] += sum_j K_j x_j for up to 32 * kRows rows: a warp per block
// (lane l takes rows l, l + 32, ...), warp w taking blocks w, w + 8, ...;
// the warps' sums are added in warp order, so the result does not depend on
// scheduling.  No barrier inside the loop over blocks.
template<typename T, int kRows>
__global__ void __launch_bounds__(kH2AccumThreads) h2_accumulate_warp_kernel(const H2Box* boxes,
                                                                             const H2BlockRef* blocks, T* vec,
                                                                             int64_t ld, const T* local, const T* ghost,
                                                                             const T* partial, int nrhs, int near) {
    constexpr int kWarps = kH2AccumThreads / 32;
    const H2Box b = boxes[blockIdx.x];
    const int m = near ? b.n : b.k;
    const int first = near ? b.near0 : b.block0;
    const int count = near ? b.nnear : b.nblocks;
    if (m == 0 || m > 32 * kRows || count == 0) return;  // (larger boxes: h2_accumulate_kernel)
    __shared__ T part[kWarps][32 * kRows];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T s[kRows];
        #pragma unroll
        for (int q = 0; q < kRows; ++q) s[q] = T(0.0);
        for (int j = warp; j < count; j += kWarps) {
            const H2BlockRef br = blocks[first + j];
            const T* K = static_cast<const T*>(br.K);
            const T* x = (br.ghost == 0 ? local : br.ghost == 1 ? ghost : partial) + br.base * nrhs +
                         static_cast<int64_t>(c) * br.stride + br.off;
            T t[kRows];
            #pragma unroll
            for (int q = 0; q < kRows; ++q) t[q] = T(0.0);
            if (K == nullptr) {  // a precomputed vector
                #pragma unroll
                for (int q = 0; q < kRows; ++q) {
                    const int i = lane + 32 * q;
                    if (i < m) t[q] = x[i];
                }
            }
            for (int cc = 0; K != nullptr && cc < br.cols; ++cc) {
                const T xc = x[cc];
                const T* col = K + static_cast<int64_t>(cc) * m;
                #pragma unroll
                for (int q = 0; q < kRows; ++q) {
                    const int i = lane + 32 * q;
                    if (i < m) t[q] += col[i] * xc;
                }
            }
            #pragma unroll
            for (int q = 0; q < kRows; ++q) s[q] += t[q];
        }
        #pragma unroll
        for (int q = 0; q < kRows; ++q) part[warp][lane + 32 * q] = s[q];
        __syncthreads();
        T* y = vec + static_cast<int64_t>(c) * ld + b.vec;
        for (int i = threadIdx.x; i < m; i += blockDim.x) {
            T add = part[0][i];
            for (int w = 1; w < kWarps; ++w) add += part[w][i];
            T& yi = near ? y[i] : y[b.skel[i]];
            yi += add;
        }
        __syncthreads();
    }
}

// y[rows] += sum_j K_j x_j.  The rows of a pass are split over column
// groups whose partial sums are added in group order, and the blocks in
// list order: the sums do not depend on scheduling.
template<typename T>
__global__ void __launch_bounds__(kH2AccumThreads) h2_accumulate_kernel(const H2Box* boxes, const H2BlockRef* blocks,
                                                                        T* vec, int64_t ld, const T* local,
                                                                        const T* ghost, const T* partial, int nrhs,
                                                                        int near, int skip_small) {
    const H2Box b = boxes[blockIdx.x];
    const int m = near ? b.n : b.k;
    const int first = near ? b.near0 : b.block0;
    const int count = near ? b.nnear : b.nblocks;
    if (m == 0 || count == 0 || m <= skip_small) return;  // (smaller boxes: h2_accumulate_warp_kernel)
    T* acc = shared_array<T>();          // m
    T* part = acc + m;                   // kH2AccumThreads
    T* xs = part + kH2AccumThreads;      // the current block's source vector
    const int tid = threadIdx.x;
    const int R = m < kH2AccumThreads ? m : kH2AccumThreads;  // rows per pass
    const int G = kH2AccumThreads / R;                         // column groups
    const int ri = tid % R, g = tid / R;
    for (int c = 0; c < nrhs; ++c) {
        T* y = vec + static_cast<int64_t>(c) * ld + b.vec;
        for (int i = tid; i < m; i += blockDim.x) acc[i] = near ? y[i] : y[b.skel[i]];
        for (int j = 0; j < count; ++j) {
            const H2BlockRef br = blocks[first + j];
            const T* K = static_cast<const T*>(br.K);
            const T* x = (br.ghost == 0 ? local : br.ghost == 1 ? ghost : partial) + br.base * nrhs +
                         static_cast<int64_t>(c) * br.stride + br.off;
            if (K == nullptr) {  // a precomputed vector
                __syncthreads();
                for (int i = tid; i < m; i += blockDim.x) acc[i] += x[i];
                continue;
            }
            __syncthreads();  // the previous block is done with xs
            for (int e = tid; e < br.cols; e += blockDim.x) xs[e] = x[e];
            __syncthreads();
            for (int r0 = 0; r0 < m; r0 += R) {
                const int i = r0 + ri;
                if (g < G && i < m) {
                    T s = T(0.0);
                    for (int cc = g; cc < br.cols; cc += G) s += K[i + static_cast<int64_t>(cc) * m] * xs[cc];
                    part[g * R + ri] = s;
                }
                __syncthreads();
                if (g == 0 && i < m) {
                    T s = part[ri];
                    for (int h = 1; h < G; ++h) s += part[h * R + ri];
                    acc[i] += s;
                }
                __syncthreads();
            }
        }
        __syncthreads();
        for (int i = tid; i < m; i += blockDim.x) {
            if (near) {
                y[i] = acc[i];
            } else {
                y[b.skel[i]] = acc[i];
            }
        }
        __syncthreads();
    }
}

// x_R += T^T x_S: a warp per redundant entry
template<typename T>
__global__ void __launch_bounds__(kThreads) h2_downward_kernel(const H2Box* boxes, T* vec, int64_t ld, int nrhs) {
    const H2Box b = boxes[blockIdx.x];
    if (b.r == 0 || b.T == nullptr) return;
    T* xs = shared_array<T>();
    const T* Tm = static_cast<const T*>(b.T);
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarps = blockDim.x >> 5;
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + static_cast<int64_t>(c) * ld + b.vec;
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) xs[i] = v[b.skel[i]];
        __syncthreads();
        for (int j = warp; j < b.r; j += nwarps) {
            const T* col = Tm + static_cast<int64_t>(j) * b.k;
            T s = T(0.0);
            for (int i = lane; i < b.k; i += 32) s += col[i] * xs[i];
            s = warp_sum(s);
            if (lane == 0) v[b.red[j]] += s;
        }
        __syncthreads();
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) h2_handoff_kernel(const H2Handoff* items, T* child, int64_t ldc, T* parent,
                                                              int64_t ldp, int nrhs, int gather) {
    const H2Handoff it = items[blockIdx.x];
    for (int c = 0; c < nrhs; ++c) {
        T* x = child + static_cast<int64_t>(c) * ldc + it.child;
        T* p = parent + static_cast<int64_t>(c) * ldp + it.parent;
        for (int i = threadIdx.x; i < it.k; i += blockDim.x) {
            if (gather) {
                p[i] = x[it.skel[i]];
            } else {
                x[it.skel[i]] = p[i];
            }
        }
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) h2_pack_kernel(const H2Span* spans, const T* src, T* dst, int nrhs) {
    const H2Span sp = spans[blockIdx.x];
    for (int c = 0; c < nrhs; ++c) {
        const T* s = src + static_cast<int64_t>(c) * sp.src_stride + sp.src;
        T* d = dst + sp.dst_base * nrhs + static_cast<int64_t>(c) * sp.dst_stride + sp.dst_off;
        for (int i = threadIdx.x; i < sp.len; i += blockDim.x) d[i] = s[i];
    }
}

__device__ __forceinline__ double warp_reduce(double v) { return warp_sum(v); }
__device__ __forceinline__ dcomplex warp_reduce(dcomplex v) { return warp_sum(v); }

// One pair per thread block, K read once: warp w takes the row strips w, w
// + 8, ... (32 rows); a = K x_s row by row, b = K^T x_t from the strips'
// column sums, added in warp order.
template<typename T>
__global__ void __launch_bounds__(kThreads) h2_pairs_kernel(const H2Pair* pairs, const T* x, int64_t ldx, T* partial,
                                                            int64_t ldp, int nrhs) {
    constexpr int kWarps = kThreads / 32;
    const H2Pair pr = pairs[blockIdx.x];
    const T* K = static_cast<const T*>(pr.K);
    T* bsum = shared_array<T>();  // kWarps x cols
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    const int m = pr.rows, n = pr.cols;
    for (int c = 0; c < nrhs; ++c) {
        const T* xs = x + static_cast<int64_t>(c) * ldx + pr.xs;
        const T* xt = x + static_cast<int64_t>(c) * ldx + pr.xt;
        T* a = partial + static_cast<int64_t>(c) * ldp + pr.pa;
        T* b = partial + static_cast<int64_t>(c) * ldp + pr.pb;
        for (int j = threadIdx.x; j < kWarps * n; j += blockDim.x) bsum[j] = T(0.0);
        __syncthreads();
        for (int r0 = 32 * warp; r0 < m; r0 += 32 * kWarps) {
            const int i = r0 + lane;
            const T xi = i < m ? xt[i] : T(0.0);
            T ai = T(0.0);
            for (int cc = 0; cc < n; ++cc) {
                const T v = i < m ? K[i + static_cast<int64_t>(cc) * m] : T(0.0);
                ai += v * xs[cc];
                const T col = warp_reduce(v * xi);
                if (lane == 0) bsum[warp * n + cc] += col;
            }
            if (i < m) a[i] = ai;
        }
        __syncthreads();
        for (int j = threadIdx.x; j < n; j += blockDim.x) {
            T sum = bsum[j];
            for (int w = 1; w < kWarps; ++w) sum += bsum[w * n + j];
            b[j] = sum;
        }
        __syncthreads();
    }
}

template<typename Kernel>
void allow_shared(Kernel kernel, size_t bytes) {
    if (bytes > kH2MatvecSharedLimit) throw std::runtime_error("device H2 matvec: box too large for shared memory");
    if (bytes > 48 * 1024) {
        cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize,
                             static_cast<int>(kH2MatvecSharedLimit));
    }
}

}  // namespace

template<typename T>
void launch_h2_upward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_r, cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t shared = static_cast<size_t>(max_r) * sizeof(T);
    allow_shared(h2_upward_kernel<T>, shared);
    h2_upward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs);
    check_launch("h2_upward_kernel");
}

template<typename T>
void launch_h2_gather_q(const H2Box* boxes, int count, const T* vec, int64_t ld, T* q, int64_t ldq, int nrhs,
                        cudaStream_t stream) {
    if (count <= 0) return;
    h2_gather_q_kernel<T><<<count, kThreads, 0, stream>>>(boxes, vec, ld, q, ldq, nrhs);
    check_launch("h2_gather_q_kernel");
}

template<typename T>
void launch_h2_accumulate(const H2Box* boxes, int count, const H2BlockRef* blocks, T* vec, int64_t ld, const T* local,
                          const T* ghost, const T* partial, int nrhs, bool near, int max_rows, int max_cols,
                          cudaStream_t stream) {
    if (count <= 0 || max_rows <= 0) return;
    // boxes of up to 256 rows: a warp per block; larger ones: a pass per block
    const int nearf = near ? 1 : 0;
    int covered = 0;
    if (max_rows <= 32) {
        h2_accumulate_warp_kernel<T, 1><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf);
        covered = 32;
    } else if (max_rows <= 64) {
        h2_accumulate_warp_kernel<T, 2><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf);
        covered = 64;
    } else if (max_rows <= 128) {
        h2_accumulate_warp_kernel<T, 4><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf);
        covered = 128;
    } else {
        h2_accumulate_warp_kernel<T, 8><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf);
        covered = 256;
    }
    check_launch("h2_accumulate_warp_kernel");
    if (max_rows <= covered) return;
    const size_t shared = h2_accumulate_shared<T>(max_rows, max_cols);
    allow_shared(h2_accumulate_kernel<T>, shared);
    h2_accumulate_kernel<T><<<count, kH2AccumThreads, shared, stream>>>(boxes, blocks, vec, ld, local, ghost, partial,
                                                                         nrhs, nearf, covered);
    check_launch("h2_accumulate_kernel");
}

template<typename T>
void launch_h2_downward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_k, cudaStream_t stream) {
    if (count <= 0 || max_k <= 0) return;
    const size_t shared = static_cast<size_t>(max_k) * sizeof(T);
    allow_shared(h2_downward_kernel<T>, shared);
    h2_downward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs);
    check_launch("h2_downward_kernel");
}

template<typename T>
void launch_h2_pairs(const H2Pair* pairs, int count, const T* x, int64_t ldx, T* partial, int64_t ldp, int nrhs,
                     int max_cols, cudaStream_t stream) {
    if (count <= 0 || max_cols <= 0) return;
    if (max_cols > kH2PairMaxCols) throw std::runtime_error("device H2 matvec: near block too wide for the pair kernel");
    const size_t shared = static_cast<size_t>(kThreads / 32) * max_cols * sizeof(T);
    allow_shared(h2_pairs_kernel<T>, shared);
    h2_pairs_kernel<T><<<count, kThreads, shared, stream>>>(pairs, x, ldx, partial, ldp, nrhs);
    check_launch("h2_pairs_kernel");
}

template<typename T>
void launch_h2_handoff(const H2Handoff* items, int count, T* child, int64_t ldc, T* parent, int64_t ldp, int nrhs,
                       bool gather, cudaStream_t stream) {
    if (count <= 0) return;
    h2_handoff_kernel<T><<<count, kThreads, 0, stream>>>(items, child, ldc, parent, ldp, nrhs, gather ? 1 : 0);
    check_launch("h2_handoff_kernel");
}

template<typename T>
void launch_h2_pack(const H2Span* spans, int count, const T* src, T* dst, int nrhs, int max_len, cudaStream_t stream) {
    if (count <= 0 || max_len <= 0) return;
    h2_pack_kernel<T><<<count, kThreads, 0, stream>>>(spans, src, dst, nrhs);
    check_launch("h2_pack_kernel");
}

#define H2_MATVEC_LAUNCHES(T)                                                                                     \
    template void launch_h2_upward<T>(const H2Box*, int, T*, int64_t, int, int, cudaStream_t);                   \
    template void launch_h2_gather_q<T>(const H2Box*, int, const T*, int64_t, T*, int64_t, int, cudaStream_t);   \
    template void launch_h2_accumulate<T>(const H2Box*, int, const H2BlockRef*, T*, int64_t, const T*, const T*,  \
                                          const T*, int, bool, int, int, cudaStream_t);                          \
    template void launch_h2_pairs<T>(const H2Pair*, int, const T*, int64_t, T*, int64_t, int, int, cudaStream_t); \
    template void launch_h2_downward<T>(const H2Box*, int, T*, int64_t, int, int, cudaStream_t);                 \
    template void launch_h2_handoff<T>(const H2Handoff*, int, T*, int64_t, T*, int64_t, int, bool, cudaStream_t); \
    template void launch_h2_pack<T>(const H2Span*, int, const T*, T*, int, int, cudaStream_t);
H2_MATVEC_LAUNCHES(double)
H2_MATVEC_LAUNCHES(dcomplex)
#undef H2_MATVEC_LAUNCHES

}  // namespace gpu
}  // namespace fmm
