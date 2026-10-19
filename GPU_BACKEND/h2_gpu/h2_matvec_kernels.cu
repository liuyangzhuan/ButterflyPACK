// Kernels of the device matvec of the compression-only H2.  See
// h2_matvec_kernels.hpp.

#include "h2_matvec_kernels.hpp"

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

// A target of the transposed products a thread block: warp w takes its
// entries w, w + 8, ..., lane l the rows l, l + 32, ... of every block (in
// order), the lanes' sums added by warp_sum: the result does not depend on
// scheduling.
template<typename T>
__global__ void __launch_bounds__(kThreads) h2_accumulate_t_kernel(const H2TTarget* targets, const H2TBlock* blocks,
                                                                   const T* x, int64_t ldx, T* local, T* send,
                                                                   int nrhs) {
    const H2TTarget t = targets[blockIdx.x];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarps = blockDim.x >> 5;
    T* out0 = (t.dest == 0 ? local : send) + t.base * nrhs + t.off;
    for (int c = 0; c < nrhs; ++c) {
        T* out = out0 + static_cast<int64_t>(c) * t.stride;
        for (int i = warp; i < t.len; i += nwarps) {
            T s = T(0.0);
            for (int j = t.block0; j < t.block0 + t.nblocks; ++j) {
                const H2TBlock b = blocks[j];
                const T* col = static_cast<const T*>(b.K) + static_cast<int64_t>(i) * b.m;
                const T* xs = x + static_cast<int64_t>(c) * ldx + b.x;
                for (int r = lane; r < b.m; r += 32) s += col[r] * xs[r];
            }
            s = warp_sum(s);
            if (lane == 0) out[i] += s;
        }
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) h2_unpack_add_kernel(const H2Span* spans, const T* dst, T* src, int nrhs) {
    const H2Span sp = spans[blockIdx.x];
    for (int c = 0; c < nrhs; ++c) {
        T* s = src + static_cast<int64_t>(c) * sp.src_stride + sp.src;
        const T* d = dst + sp.dst_base * nrhs + static_cast<int64_t>(c) * sp.dst_stride + sp.dst_off;
        for (int i = threadIdx.x; i < sp.len; i += blockDim.x) s[i] += d[i];
    }
}

template<typename T>
__global__ void __launch_bounds__(kThreads) h2_scatter_add_q_kernel(const H2Box* boxes, const T* q, int64_t ldq, T* vec,
                                                                    int64_t ld, int nrhs) {
    const H2Box b = boxes[blockIdx.x];
    for (int c = 0; c < nrhs; ++c) {
        T* v = vec + static_cast<int64_t>(c) * ld + b.vec;
        const T* p = q + static_cast<int64_t>(c) * ldq + b.q;
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) v[b.skel[i]] += p[i];
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

// ---- passes of several columns
//
// The kernels above take the columns one at a time and read a box's blocks
// again for each; these take up to kMaxCols at once (column q's vectors in
// shared memory at q times a stride) and apply each block element to all of
// them.  A column's arithmetic is the single-column order - the same sums
// in the same order - so its result does not depend on how the columns are
// grouped.  A launch whose passes would take one column runs the kernels
// above.

template<typename T>
constexpr int kMaxCols = std::is_same_v<T, double> ? 8 : 4;

// Whether unrolled column q is in the pass (column 0 always is).
__device__ __forceinline__ bool live(int q, int nc) { return q == 0 || q < nc; }

// x_S += T x_R
template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) h2_upward_cols_kernel(const H2Box* boxes, T* vec, int64_t ld, int nrhs,
                                                                  int cols) {
    const H2Box b = boxes[blockIdx.x];
    if (b.r == 0 || b.T == nullptr) return;
    T* xr = shared_array<T>();  // column q at xr + q * r
    const T* Tm = static_cast<const T*>(b.T);
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + static_cast<int64_t>(c0) * ld + b.vec;
        for (int q = 0; q < nc; ++q) {
            const T* v = v0 + static_cast<int64_t>(q) * ld;
            for (int j = threadIdx.x; j < b.r; j += blockDim.x) xr[q * b.r + j] = v[b.red[j]];
        }
        __syncthreads();
        for (int i = threadIdx.x; i < b.k; i += blockDim.x) {
            T s[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) s[q] = T(0.0);
            const T* a = Tm + i;
            const int64_t lda = b.k;
            int j = 0;
            for (; j + 3 < b.r; j += 4) {  // (the loads of four j in flight together)
                const T a0 = a[j * lda], a1 = a[(j + 1) * lda], a2 = a[(j + 2) * lda], a3 = a[(j + 3) * lda];
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) {
                        const T* x = xr + q * b.r;
                        s[q] += a0 * x[j];
                        s[q] += a1 * x[j + 1];
                        s[q] += a2 * x[j + 2];
                        s[q] += a3 * x[j + 3];
                    }
                }
            }
            for (; j < b.r; ++j) {
                const T a0 = a[j * lda];
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) s[q] += a0 * xr[q * b.r + j];
                }
            }
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                if (live(q, nc)) v0[static_cast<int64_t>(q) * ld + b.skel[i]] += s[q];
            }
        }
        __syncthreads();
    }
}

// y[rows] += sum_j K_j x_j for boxes of up to `covered` rows, as
// h2_accumulate_warp_kernel a column: a warp per block (lane l takes row l of
// a 32-row strip), warp w taking blocks w, w + 8, ..., the warps' sums added
// in warp order.  The strips go one after another (all warps alike), each
// through all of the box's blocks, so a lane holds one strip's sums.
template<typename T, int CB>
__global__ void __launch_bounds__(kH2AccumThreads) h2_accumulate_warp_cols_kernel(
    const H2Box* boxes, const H2BlockRef* blocks, T* vec, int64_t ld, const T* local, const T* ghost,
    const T* partial, int nrhs, int near, int covered, int cols) {
    constexpr int kWarps = kH2AccumThreads / 32;
    const H2Box b = boxes[blockIdx.x];
    const int m = near ? b.n : b.k;
    const int first = near ? b.near0 : b.block0;
    const int count = near ? b.nnear : b.nblocks;
    if (m == 0 || m > covered || count == 0) return;  // (larger boxes: h2_accumulate_cols_kernel)
    __shared__ T part[kWarps][CB][32];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        for (int r0 = 0; r0 < m; r0 += 32) {
            const int i = r0 + lane;
            T s[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) s[q] = T(0.0);
            for (int j = warp; j < count; j += kWarps) {
                const H2BlockRef br = blocks[first + j];
                const T* K = static_cast<const T*>(br.K);
                const T* x = (br.ghost == 0 ? local : br.ghost == 1 ? ghost : partial) + br.base * nrhs +
                             static_cast<int64_t>(c0) * br.stride + br.off;  // column q at x + q * stride
                T t[CB];
#pragma unroll
                for (int q = 0; q < CB; ++q) t[q] = T(0.0);
                if (K == nullptr) {  // a precomputed vector
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc) && i < m) t[q] = x[q * br.stride + i];
                    }
                }
                for (int cc = 0; K != nullptr && cc < br.cols; ++cc) {
                    const T kv = i < m ? K[i + static_cast<int64_t>(cc) * m] : T(0.0);
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc)) {
                            const T xc = x[q * br.stride + cc];
                            if (i < m) t[q] += kv * xc;
                        }
                    }
                }
#pragma unroll
                for (int q = 0; q < CB; ++q) s[q] += t[q];
            }
#pragma unroll
            for (int q = 0; q < CB; ++q) part[warp][q][lane] = s[q];
            __syncthreads();
            for (int e = threadIdx.x; e < 32 * nc; e += blockDim.x) {
                const int l = e & 31, q = e >> 5, row = r0 + l;
                if (row >= m) continue;
                T add = part[0][q][l];
                for (int w = 1; w < kWarps; ++w) add += part[w][q][l];
                T* y = vec + static_cast<int64_t>(c0 + q) * ld + b.vec;
                T& yi = near ? y[row] : y[b.skel[row]];
                yi += add;
            }
            __syncthreads();
        }
    }
}

// y[rows] += sum_j K_j x_j for the larger boxes, as h2_accumulate_kernel a
// column (column q's acc at q * mr, its part at q * kH2AccumThreads, its
// source vector at q * mc; mr, mc: the launch's largest rows and columns).
template<typename T, int CB>
__global__ void __launch_bounds__(kH2AccumThreads) h2_accumulate_cols_kernel(
    const H2Box* boxes, const H2BlockRef* blocks, T* vec, int64_t ld, const T* local, const T* ghost,
    const T* partial, int nrhs, int near, int skip_small, int cols, int mr, int mc) {
    const H2Box b = boxes[blockIdx.x];
    const int m = near ? b.n : b.k;
    const int first = near ? b.near0 : b.block0;
    const int count = near ? b.nnear : b.nblocks;
    if (m == 0 || count == 0 || m <= skip_small) return;  // (smaller boxes: h2_accumulate_warp_cols_kernel)
    T* acc = shared_array<T>();
    T* part = acc + static_cast<int64_t>(cols) * mr;
    T* xs = part + static_cast<int64_t>(cols) * kH2AccumThreads;
    const int tid = threadIdx.x;
    const int R = m < kH2AccumThreads ? m : kH2AccumThreads;  // rows per pass
    const int G = kH2AccumThreads / R;                         // column groups
    const int ri = tid % R, g = tid / R;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* y0 = vec + static_cast<int64_t>(c0) * ld + b.vec;
        for (int q = 0; q < nc; ++q) {
            const T* y = y0 + static_cast<int64_t>(q) * ld;
            for (int i = tid; i < m; i += blockDim.x) acc[q * mr + i] = near ? y[i] : y[b.skel[i]];
        }
        for (int j = 0; j < count; ++j) {
            const H2BlockRef br = blocks[first + j];
            const T* K = static_cast<const T*>(br.K);
            const T* x = (br.ghost == 0 ? local : br.ghost == 1 ? ghost : partial) + br.base * nrhs +
                         static_cast<int64_t>(c0) * br.stride + br.off;  // column q at x + q * stride
            if (K == nullptr) {  // a precomputed vector
                __syncthreads();
                for (int q = 0; q < nc; ++q) {
                    for (int i = tid; i < m; i += blockDim.x) acc[q * mr + i] += x[q * br.stride + i];
                }
                continue;
            }
            __syncthreads();  // the previous block is done with xs
            for (int q = 0; q < nc; ++q) {
                for (int e = tid; e < br.cols; e += blockDim.x) xs[q * mc + e] = x[q * br.stride + e];
            }
            __syncthreads();
            for (int r0 = 0; r0 < m; r0 += R) {
                const int i = r0 + ri;
                if (g < G && i < m) {
                    T s[CB];
#pragma unroll
                    for (int q = 0; q < CB; ++q) s[q] = T(0.0);
                    for (int cc = g; cc < br.cols; cc += G) {
                        const T kv = K[i + static_cast<int64_t>(cc) * m];
#pragma unroll
                        for (int q = 0; q < CB; ++q) {
                            if (live(q, nc)) s[q] += kv * xs[q * mc + cc];
                        }
                    }
#pragma unroll
                    for (int q = 0; q < CB; ++q) {
                        if (live(q, nc)) part[q * kH2AccumThreads + g * R + ri] = s[q];
                    }
                }
                __syncthreads();
                if (g == 0 && i < m) {
                    for (int q = 0; q < nc; ++q) {
                        T s = part[q * kH2AccumThreads + ri];
                        for (int h = 1; h < G; ++h) s += part[q * kH2AccumThreads + h * R + ri];
                        acc[q * mr + i] += s;
                    }
                }
                __syncthreads();
            }
        }
        __syncthreads();
        for (int q = 0; q < nc; ++q) {
            T* y = y0 + static_cast<int64_t>(q) * ld;
            for (int i = tid; i < m; i += blockDim.x) {
                if (near) {
                    y[i] = acc[q * mr + i];
                } else {
                    y[b.skel[i]] = acc[q * mr + i];
                }
            }
        }
        __syncthreads();
    }
}

// x_R += T^T x_S: a warp per redundant entry
template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) h2_downward_cols_kernel(const H2Box* boxes, T* vec, int64_t ld, int nrhs,
                                                                    int cols) {
    const H2Box b = boxes[blockIdx.x];
    if (b.r == 0 || b.T == nullptr) return;
    T* xs = shared_array<T>();  // column q at xs + q * k
    const T* Tm = static_cast<const T*>(b.T);
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarps = blockDim.x >> 5;
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        T* v0 = vec + static_cast<int64_t>(c0) * ld + b.vec;
        for (int q = 0; q < nc; ++q) {
            const T* v = v0 + static_cast<int64_t>(q) * ld;
            for (int i = threadIdx.x; i < b.k; i += blockDim.x) xs[q * b.k + i] = v[b.skel[i]];
        }
        __syncthreads();
        for (int j = warp; j < b.r; j += nwarps) {
            const T* col = Tm + static_cast<int64_t>(j) * b.k;
            T s[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) s[q] = T(0.0);
            for (int i = lane; i < b.k; i += 32) {
                const T t = col[i];
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) s[q] += t * xs[q * b.k + i];
                }
            }
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                if (live(q, nc)) s[q] = warp_sum(s[q]);
            }
            if (lane == 0) {
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (live(q, nc)) v0[static_cast<int64_t>(q) * ld + b.red[j]] += s[q];
                }
            }
        }
        __syncthreads();
    }
}

__device__ __forceinline__ double shfl_xor(double v, int o) { return __shfl_xor_sync(0xffffffffu, v, o); }
__device__ __forceinline__ dcomplex shfl_xor(dcomplex v, int o) {
    return dcomplex(__shfl_xor_sync(0xffffffffu, v.re, o), __shfl_xor_sync(0xffffffffu, v.im, o));
}

// The sums over the warp of p[0..15] (lane l's values), each bitwise what
// warp_sum gives it - the xor butterfly, every lane adding its partner's
// value to its own, offsets 16, 8, 4, 2, 1; every lane ends with the same
// sum - in 31 shuffles instead of 80: offset 16 adds all sixteen, and each
// later offset halves what a lane keeps (the half its bit of the offset
// selects), so lane l ends holding the sum of index l & 15.
template<int O, typename T>
__device__ __forceinline__ void halve_sum(T p[16], int lane) {
    const bool upper = (lane & O) != 0;
#pragma unroll
    for (int u = 0; u < O; ++u) {
        const T send = upper ? p[u] : p[u + O];
        const T keep = upper ? p[u + O] : p[u];
        p[u] = keep + shfl_xor(send, O);
    }
}
template<typename T>
__device__ __forceinline__ T transposed_sum16(T p[16], int lane) {
#pragma unroll
    for (int u = 0; u < 16; ++u) p[u] = p[u] + shfl_xor(p[u], 16);
    halve_sum<8>(p, lane);
    halve_sum<4>(p, lane);
    halve_sum<2>(p, lane);
    halve_sum<1>(p, lane);
    return p[0];
}

// The pairs, as h2_pairs_kernel a column (column q's column sums at bsum + q
// * kWarps * mc; mc: the launch's widest block).  K is read once a pass, in
// groups of 16 columns whose b sums (K^T x_t over the strip) go through
// transposed_sum16: what dominated the pairs was a warp reduction for every
// block column and right-hand side.
template<typename T, int CB>
__global__ void __launch_bounds__(kThreads) h2_pairs_cols_kernel(const H2Pair* pairs, const T* x, int64_t ldx,
                                                                 T* partial, int64_t ldp, int nrhs, int cols, int mc) {
    constexpr int kWarps = kThreads / 32, kG = 16;
    const H2Pair pr = pairs[blockIdx.x];
    const T* K = static_cast<const T*>(pr.K);
    T* bsum = shared_array<T>();
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    const int m = pr.rows, n = pr.cols;
    const int64_t bs = static_cast<int64_t>(kWarps) * mc;  // a column's column sums
    for (int c0 = 0; c0 < nrhs; c0 += cols) {
        const int nc = min(cols, nrhs - c0);
        for (int q = 0; q < nc; ++q) {
            for (int j = threadIdx.x; j < kWarps * n; j += blockDim.x) bsum[q * bs + j] = T(0.0);
        }
        __syncthreads();
        for (int r0 = 32 * warp; r0 < m; r0 += 32 * kWarps) {
            const int i = r0 + lane;
            T xi[CB], ai[CB];
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                xi[q] = live(q, nc) && i < m ? x[static_cast<int64_t>(c0 + q) * ldx + pr.xt + i] : T(0.0);
                ai[q] = T(0.0);
            }
            for (int cc0 = 0; cc0 < n; cc0 += kG) {
                T kv[kG];
#pragma unroll
                for (int u = 0; u < kG; ++u) {
                    kv[u] = i < m && cc0 + u < n ? K[i + static_cast<int64_t>(cc0 + u) * m] : T(0.0);
                }
#pragma unroll
                for (int q = 0; q < CB; ++q) {
                    if (!live(q, nc)) continue;
                    const T* xs = x + static_cast<int64_t>(c0 + q) * ldx + pr.xs + cc0;
#pragma unroll
                    for (int u = 0; u < kG; ++u) {
                        if (cc0 + u < n) ai[q] += kv[u] * xs[u];
                    }
                    T pp[kG];
#pragma unroll
                    for (int u = 0; u < kG; ++u) pp[u] = kv[u] * xi[q];
                    const T col = transposed_sum16(pp, lane);
                    const int cc = cc0 + (lane & 15);
                    if (lane < 16 && cc < n) bsum[q * bs + warp * n + cc] += col;
                }
            }
#pragma unroll
            for (int q = 0; q < CB; ++q) {
                if (live(q, nc) && i < m) partial[static_cast<int64_t>(c0 + q) * ldp + pr.pa + i] = ai[q];
            }
        }
        __syncthreads();
        for (int q = 0; q < nc; ++q) {
            for (int j = threadIdx.x; j < n; j += blockDim.x) {
                T sum = bsum[q * bs + j];
                for (int w = 1; w < kWarps; ++w) sum += bsum[q * bs + w * n + j];
                partial[static_cast<int64_t>(c0 + q) * ldp + pr.pb + j] = sum;
            }
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

// The columns a pass takes when each needs `per_column` bytes of shared
// memory: one, or up to kMaxCols<T> within half the limit (two blocks still
// fit on an SM).
template<typename T>
int pass_columns(int nrhs, size_t per_column) {
    if (nrhs <= 1) return 1;
    const size_t fit = (kH2MatvecSharedLimit / 2) / std::max<size_t>(per_column, 1);
    const size_t cols = std::min<size_t>({static_cast<size_t>(nrhs), static_cast<size_t>(kMaxCols<T>), fit});
    return static_cast<int>(std::max<size_t>(cols, 1));
}
// f(std::integral_constant<int, CB>) with the smallest instantiated CB (2, 4
// or 8) that holds `cols` (>= 2) columns
template<typename T, typename F>
void with_cols(int cols, F&& f) {
    if (cols <= 2) {
        f(std::integral_constant<int, 2>{});
    } else if (cols <= 4 || kMaxCols<T> < 8) {
        f(std::integral_constant<int, 4>{});
    } else {
        f(std::integral_constant<int, (kMaxCols<T> < 8 ? 4 : 8)>{});
    }
}

}  // namespace

template<typename T>
void launch_h2_upward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_r, cudaStream_t stream) {
    if (count <= 0 || max_r <= 0) return;
    const size_t per_column = static_cast<size_t>(max_r) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(h2_upward_kernel<T>, shared);
        h2_upward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(h2_upward_cols_kernel<T, CB>, shared);
            h2_upward_cols_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs, cols);
        });
    }
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
    const int covered = max_rows <= 32 ? 32 : max_rows <= 64 ? 64 : max_rows <= 128 ? 128 : 256;
    if (nrhs == 1) {
        if (max_rows <= 32) {
            h2_accumulate_warp_kernel<T, 1><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local,
                                                                                   ghost, partial, nrhs, nearf);
        } else if (max_rows <= 64) {
            h2_accumulate_warp_kernel<T, 2><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local,
                                                                                   ghost, partial, nrhs, nearf);
        } else if (max_rows <= 128) {
            h2_accumulate_warp_kernel<T, 4><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local,
                                                                                   ghost, partial, nrhs, nearf);
        } else {
            h2_accumulate_warp_kernel<T, 8><<<count, kH2AccumThreads, 0, stream>>>(boxes, blocks, vec, ld, local,
                                                                                   ghost, partial, nrhs, nearf);
        }
    } else {  // (its shared memory is static: a pass takes kMaxCols columns)
        const int cols = std::min(nrhs, kMaxCols<T>);
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            h2_accumulate_warp_cols_kernel<T, CB><<<count, kH2AccumThreads, 0, stream>>>(
                boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf, covered, cols);
        });
    }
    check_launch("h2_accumulate_warp_kernel");
    if (max_rows <= covered) return;
    const size_t per_column = h2_accumulate_shared<T>(max_rows, max_cols);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(h2_accumulate_kernel<T>, shared);
        h2_accumulate_kernel<T><<<count, kH2AccumThreads, shared, stream>>>(boxes, blocks, vec, ld, local, ghost,
                                                                             partial, nrhs, nearf, covered);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(h2_accumulate_cols_kernel<T, CB>, shared);
            h2_accumulate_cols_kernel<T, CB><<<count, kH2AccumThreads, shared, stream>>>(
                boxes, blocks, vec, ld, local, ghost, partial, nrhs, nearf, covered, cols, max_rows, max_cols);
        });
    }
    check_launch("h2_accumulate_kernel");
}

template<typename T>
void launch_h2_downward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_k, cudaStream_t stream) {
    if (count <= 0 || max_k <= 0) return;
    const size_t per_column = static_cast<size_t>(max_k) * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(h2_downward_kernel<T>, shared);
        h2_downward_kernel<T><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(h2_downward_cols_kernel<T, CB>, shared);
            h2_downward_cols_kernel<T, CB><<<count, kThreads, shared, stream>>>(boxes, vec, ld, nrhs, cols);
        });
    }
    check_launch("h2_downward_kernel");
}

template<typename T>
void launch_h2_pairs(const H2Pair* pairs, int count, const T* x, int64_t ldx, T* partial, int64_t ldp, int nrhs,
                     int max_cols, cudaStream_t stream) {
    if (count <= 0 || max_cols <= 0) return;
    if (max_cols > kH2PairMaxCols) throw std::runtime_error("device H2 matvec: near block too wide for the pair kernel");
    const size_t per_column = static_cast<size_t>(kThreads / 32) * max_cols * sizeof(T);
    const int cols = pass_columns<T>(nrhs, per_column);
    const size_t shared = static_cast<size_t>(cols) * per_column;
    if (cols == 1) {
        allow_shared(h2_pairs_kernel<T>, shared);
        h2_pairs_kernel<T><<<count, kThreads, shared, stream>>>(pairs, x, ldx, partial, ldp, nrhs);
    } else {
        with_cols<T>(cols, [&](auto cb) {
            constexpr int CB = decltype(cb)::value;
            allow_shared(h2_pairs_cols_kernel<T, CB>, shared);
            h2_pairs_cols_kernel<T, CB><<<count, kThreads, shared, stream>>>(pairs, x, ldx, partial, ldp, nrhs, cols,
                                                                             max_cols);
        });
    }
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

template<typename T>
void launch_h2_accumulate_t(const H2TTarget* targets, int count, const H2TBlock* blocks, const T* x, int64_t ldx,
                            T* local, T* send, int nrhs, cudaStream_t stream) {
    if (count <= 0) return;
    h2_accumulate_t_kernel<T><<<count, kThreads, 0, stream>>>(targets, blocks, x, ldx, local, send, nrhs);
    check_launch("h2_accumulate_t_kernel");
}

template<typename T>
void launch_h2_unpack_add(const H2Span* spans, int count, const T* dst, T* src, int nrhs, int max_len,
                          cudaStream_t stream) {
    if (count <= 0 || max_len <= 0) return;
    h2_unpack_add_kernel<T><<<count, kThreads, 0, stream>>>(spans, dst, src, nrhs);
    check_launch("h2_unpack_add_kernel");
}

template<typename T>
void launch_h2_scatter_add_q(const H2Box* boxes, int count, const T* q, int64_t ldq, T* vec, int64_t ld, int nrhs,
                             cudaStream_t stream) {
    if (count <= 0) return;
    h2_scatter_add_q_kernel<T><<<count, kThreads, 0, stream>>>(boxes, q, ldq, vec, ld, nrhs);
    check_launch("h2_scatter_add_q_kernel");
}

#define H2_MATVEC_LAUNCHES(T)                                                                                     \
    template void launch_h2_upward<T>(const H2Box*, int, T*, int64_t, int, int, cudaStream_t);                   \
    template void launch_h2_gather_q<T>(const H2Box*, int, const T*, int64_t, T*, int64_t, int, cudaStream_t);   \
    template void launch_h2_accumulate<T>(const H2Box*, int, const H2BlockRef*, T*, int64_t, const T*, const T*,  \
                                          const T*, int, bool, int, int, cudaStream_t);                          \
    template void launch_h2_pairs<T>(const H2Pair*, int, const T*, int64_t, T*, int64_t, int, int, cudaStream_t); \
    template void launch_h2_downward<T>(const H2Box*, int, T*, int64_t, int, int, cudaStream_t);                 \
    template void launch_h2_handoff<T>(const H2Handoff*, int, T*, int64_t, T*, int64_t, int, bool, cudaStream_t); \
    template void launch_h2_pack<T>(const H2Span*, int, const T*, T*, int, int, cudaStream_t);                   \
    template void launch_h2_accumulate_t<T>(const H2TTarget*, int, const H2TBlock*, const T*, int64_t, T*, T*,  \
                                            int, cudaStream_t);                                                  \
    template void launch_h2_unpack_add<T>(const H2Span*, int, const T*, T*, int, int, cudaStream_t);            \
    template void launch_h2_scatter_add_q<T>(const H2Box*, int, const T*, int64_t, T*, int64_t, int, cudaStream_t);
H2_MATVEC_LAUNCHES(double)
H2_MATVEC_LAUNCHES(dcomplex)
#undef H2_MATVEC_LAUNCHES

}  // namespace gpu
}  // namespace fmm
