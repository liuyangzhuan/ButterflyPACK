// The one kernel a recorded sweep runs per call (bf_sweep_n.hpp, SweepExec):
// its GEMMs' vector pointers, from offsets into the caller's x and y and the
// two ping-pong buffers, and the leading dimensions of x and y.

#include "bf_device.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace bf {
namespace gpu {

namespace {

constexpr int kThreads = 256;
// the real narrow kernels' scalar (the complex ones below take double2, the
// layout of fmm::gpu::dcomplex)
using DataType = double;

// Buffers: 0 the caller's x, 1 and 2 the ping-pong buffers, 3 the caller's y.
template<typename T>
__global__ void patch_kernel(int count, const int64_t* b_off, const int64_t* c_off, const signed char* b_buf,
                             const signed char* c_buf, const T* vin, T* buf0, T* buf1, T* vout, int ldi, int ldo,
                             const T** b, T** c, int* ldb, int* ldc) {
    for (int i = blockIdx.x * kThreads + threadIdx.x; i < count; i += gridDim.x * kThreads) {
        const int bb = b_buf[i], cb = c_buf[i];
        b[i] = (bb == 0 ? vin : (bb == 1 ? buf0 : buf1)) + b_off[i];
        c[i] = (cb == 3 ? vout : (cb == 1 ? buf0 : buf1)) + c_off[i];
        if (bb == 0) ldb[i] = ldi;
        if (cb == 3) ldc[i] = ldo;
    }
}

// C = alpha A^T B + beta C for at most kStagedMaxCols columns and
// kStagedMaxRows rows: one thread block per entry, its threads splitting C's
// columns in groups of kGroupCols, a thread per row of C within a group.  k
// goes in chunks of 32 staged in shared memory: a warp per column of A, then
// of B, its lanes along k, with kLoads columns' loads in flight before it
// stores any.  (Staging one column at a time made a chunk a chain of memory
// latencies: a 262144^2 matvec at 16 columns took 5.0 ms against the
// untransposed sweep's 3.4.)  A^T is stored with an odd pitch, so neither the
// stores nor the reads collide on a bank; B k-major, so a thread reads its
// columns as 16-byte broadcasts.  Each element of C accumulates its k
// products in order from zero, one rounding each, then takes alpha and beta
// as MAGMA's kernels do - so the result is MAGMA's, and the tensor-core
// GEMM's, to the bit.  kAccum reads C before the k loop rather than after it,
// and gets its own register bound: at 64 registers its C values spill.
constexpr int kKC = 32;
constexpr int kGroupCols = 8;
constexpr int kStagedMaxRows = 64;
constexpr int kStagedMaxCols = 32;
constexpr int kStagedMaxThreads = 256;  // 64 rows x 4 groups
constexpr int kLoads = 6;

template<bool kAccum>
__global__ void __launch_bounds__(kStagedMaxThreads, kAccum ? 3 : 4)
staged_kernel(const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a, const int* lda,
              const DataType* const* b, const int* ldb, DataType beta, DataType* const* c, const int* ldc,
              int pitch, int pitch_b, int rows_pad) {
    static_assert(std::is_same<DataType, double>::value, "B is read as double2");
    extern __shared__ DataType smem[];
    DataType* sA = smem;                                       // A^T(r, kk) at kk * pitch + r
    DataType* sB = smem + static_cast<size_t>(kKC) * pitch;    // B(kk, j) at kk * pitch_b + j
    const int e = blockIdx.x;
    const int M = m[e], N = n[e], K = k[e];
    if (M <= 0 || N <= 0) return;
    const DataType* __restrict__ A = a[e];
    const DataType* __restrict__ B = b[e];
    const int la = lda[e], lb = ldb[e];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, warps = blockDim.x >> 5;
    const int r = tid % rows_pad, j0 = (tid / rows_pad) * kGroupCols;
    const bool mine = r < M && j0 < N;
    DataType* C = c[e] + r;
    const int lc = ldc[e];
    DataType cv[kAccum ? kGroupCols : 1];
    if (kAccum && mine) {
#pragma unroll
        for (int j = 0; j < kGroupCols; ++j) cv[j] = j0 + j < N ? C[static_cast<int64_t>(j0 + j) * lc] : DataType(0);
    }
    DataType acc[kGroupCols];
#pragma unroll
    for (int j = 0; j < kGroupCols; ++j) acc[j] = DataType(0);
    const int ncols = M + N;  // a chunk's columns: A's M, then B's N
    for (int k0 = 0; k0 < K; k0 += kKC) {
        const int kc = min(kKC, K - k0);
        __syncthreads();
        if (lane < kc) {
            for (int c0 = warp; c0 < ncols; c0 += warps * kLoads) {
                DataType v[kLoads];
#pragma unroll
                for (int u = 0; u < kLoads; ++u) {
                    const int col = c0 + u * warps;
                    const DataType* src = col < M ? A + static_cast<int64_t>(col) * la
                                                  : B + static_cast<int64_t>(col - M) * lb;
                    v[u] = col < ncols ? __ldg(src + k0 + lane) : DataType(0);
                }
#pragma unroll
                for (int u = 0; u < kLoads; ++u) {
                    const int col = c0 + u * warps;
                    if (col < M) sA[lane * pitch + col] = v[u];
                    else if (col < ncols) sB[lane * pitch_b + col - M] = v[u];
                }
            }
        }
        __syncthreads();
        if (mine) {
            const DataType* ap = sA + r;
            const double2* bp = reinterpret_cast<const double2*>(sB + j0);
            const int pb2 = pitch_b / 2;
#pragma unroll 4
            for (int kk = 0; kk < kc; ++kk) {
                const DataType av = ap[kk * pitch];
#pragma unroll
                for (int q = 0; q < kGroupCols / 2; ++q) {
                    const double2 bv = bp[kk * pb2 + q];
                    acc[2 * q] = fma(av, bv.x, acc[2 * q]);
                    acc[2 * q + 1] = fma(av, bv.y, acc[2 * q + 1]);
                }
            }
        }
    }
    if (!mine) return;
#pragma unroll
    for (int j = 0; j < kGroupCols; ++j) {
        if (j0 + j >= N) continue;
        C[static_cast<int64_t>(j0 + j) * lc] = kAccum ? alpha * acc[j] + beta * cv[j] : alpha * acc[j];
    }
}

// The same products without staging, for C = alpha A B + beta C: a thread per element of
// C, reading its row of A and its column of B straight from memory (a warp's
// threads read consecutive rows of A, so every load is one line), in the
// same order.  Faster than staging for one column, where a level is a pass
// over its factors and nothing else (measured: 1.6 against 2.3 ms for a
// 262144^2 butterfly).
__global__ void __launch_bounds__(256)
direct_kernel(const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a, const int* lda,
              const DataType* const* b, const int* ldb, DataType beta, DataType* const* c, const int* ldc) {
    const int e = blockIdx.y;
    const int M = m[e], N = n[e], K = k[e];
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= M * N) return;
    const int r = t % M, j = t / M;
    const int la = lda[e];
    const DataType* bp = b[e] + static_cast<int64_t>(j) * ldb[e];
    const DataType* ap = a[e] + r;
    DataType acc = DataType(0);
#pragma unroll 8
    for (int kk = 0; kk < K; ++kk) acc = fma(__ldg(ap + static_cast<int64_t>(kk) * la), __ldg(bp + kk), acc);
    DataType* cp = c[e] + r + static_cast<int64_t>(j) * ldc[e];
    *cp = beta == DataType(0) ? alpha * acc : alpha * acc + beta * *cp;
}

// At least two warps a block: one warp alone stages too little at a time (a
// 24-row leaf factor at one column ran 10% slower).  The widest block's
// shared memory, 32 x (65 + 34) values, is under the 48 KB a launch gets.
void staged_with(int count, const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a,
                 const int* lda, const DataType* const* b, const int* ldb, DataType beta, DataType* const* c,
                 const int* ldc, int max_m, int max_n, cudaStream_t stream) {
    const int pitch = max_m | 1;  // odd
    const int rows_pad = (max_m + 31) / 32 * 32;
    const int groups = (max_n + kGroupCols - 1) / kGroupCols;
    const int pitch_b = groups * kGroupCols + 2;  // even, for the 16-byte reads
    const int threads = std::max(rows_pad * groups, 64);
    const size_t shared = static_cast<size_t>(kKC) * (pitch + pitch_b) * sizeof(DataType);
    if (beta == DataType(0)) {
        staged_kernel<false><<<count, threads, shared, stream>>>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, pitch,
                                                                 pitch_b, rows_pad);
    } else {
        staged_kernel<true><<<count, threads, shared, stream>>>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, pitch,
                                                                pitch_b, rows_pad);
    }
}

// ---- double complex: an untransposed form and the staged transposed one,
// alpha 1 and beta 0 or 1.  Each element sums its k products from zero in
// cuCfma's order (re after a.re b.re, then -a.im b.im; im after a.re b.im,
// then a.im b.re), MAGMA's zgemm's, so the result is MAGMA's to the bit
// (checked on 4096 GEMMs of 25 x 25, 26 x 64, 64 x 26 and 50 x 50 at 1-32
// columns, beta 0 and 1).  Four columns a thread (kZGroupCols): a complex
// accumulator is two doubles.
constexpr int kZGroupCols = 4;

__device__ __forceinline__ void zfma(double2& c, double2 a, double2 b) {
    c.x = __fma_rn(a.x, b.x, c.x);
    c.x = __fma_rn(-a.y, b.y, c.x);
    c.y = __fma_rn(a.x, b.y, c.y);
    c.y = __fma_rn(a.y, b.x, c.y);
}

// C = A B (+ C) for at most kZRowsMaxCols columns: a thread per row of C and
// group of columns, reading A straight from memory (a warp's threads read
// consecutive rows) kZRowsLoads values at a time, B staged in chunks of kKC.
// A batch of a few hundred GEMMs leaves most of the GPU idle, and then the
// time is the chain of A's loads from memory: with four in flight, a
// 65536^2 matvec on four A100s (256 GEMMs of 25 x 25 a batch) took 13.8 us
// a batch against MAGMA's 11.4 (nsys); with eight it is the faster.  (A
// benchmark that reuses its blocks from L2 does not show this.)  At 16
// columns MAGMA is as fast.
constexpr int kZRowsMaxCols = 8;
constexpr int kZRowsLoads = 8;

__global__ void __launch_bounds__(kStagedMaxThreads)
zrows_kernel(const int* m, const int* n, const int* k, const double2* const* a, const int* lda,
             const double2* const* b, const int* ldb, bool accum, double2* const* c, const int* ldc, int rows_pad,
             int pitch_b) {
    extern __shared__ double2 zshared[];  // B(kk, j) at kk * pitch_b + j
    const int e = blockIdx.x;
    const int M = m[e], N = n[e], K = k[e];
    if (M <= 0 || N <= 0) return;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, warps = blockDim.x >> 5;
    const int r = tid % rows_pad, j0 = (tid / rows_pad) * kZGroupCols;
    const bool mine = r < M && j0 < N;
    const double2* __restrict__ A = a[e] + r;
    const double2* __restrict__ B = b[e];
    const int la = lda[e], lb = ldb[e];
    double2 acc[kZGroupCols];
#pragma unroll
    for (int j = 0; j < kZGroupCols; ++j) acc[j] = make_double2(0.0, 0.0);
    for (int k0 = 0; k0 < K; k0 += kKC) {
        const int kc = min(kKC, K - k0);
        __syncthreads();
        if (lane < kc) {
            for (int col = warp; col < N; col += warps) {
                zshared[lane * pitch_b + col] = __ldg(B + static_cast<int64_t>(col) * lb + k0 + lane);
            }
        }
        __syncthreads();
        if (mine) {
            // kZRowsLoads of A's values in flight before any is used
            int kk = 0;
            for (; kk + kZRowsLoads <= kc; kk += kZRowsLoads) {
                double2 av[kZRowsLoads];
#pragma unroll
                for (int u = 0; u < kZRowsLoads; ++u) av[u] = __ldg(A + static_cast<int64_t>(k0 + kk + u) * la);
#pragma unroll
                for (int u = 0; u < kZRowsLoads; ++u)
#pragma unroll
                    for (int q = 0; q < kZGroupCols; ++q) zfma(acc[q], av[u], zshared[(kk + u) * pitch_b + j0 + q]);
            }
            for (; kk < kc; ++kk) {
                const double2 av = __ldg(A + static_cast<int64_t>(k0 + kk) * la);
#pragma unroll
                for (int q = 0; q < kZGroupCols; ++q) zfma(acc[q], av, zshared[kk * pitch_b + j0 + q]);
            }
        }
    }
    if (!mine) return;
    double2* C = c[e] + r;
    const int lc = ldc[e];
#pragma unroll
    for (int j = 0; j < kZGroupCols; ++j) {
        if (j0 + j >= N) continue;
        double2 v = acc[j];
        if (accum) {
            const double2 o = C[static_cast<int64_t>(j0 + j) * lc];
            v.x += o.x;
            v.y += o.y;
        }
        C[static_cast<int64_t>(j0 + j) * lc] = v;
    }
}

// C = A^T B (+ C), staged_kernel's form for complex values, at most
// kZStagedMaxCols columns.  An odd pitch keeps a quarter warp's 16-byte
// stores on distinct banks, as it does a warp's 8-byte ones.  4096 GEMMs of
// 25 x 25 on one A100, ms MAGMA -> this at beta 0: 0.119 -> 0.048 at 1
// column, 0.122 -> 0.066 at 8, 0.128 -> 0.088 at 16 (the tensor cores 0.093).
constexpr int kZStagedMaxCols = 16;

template<bool kAccum>
__global__ void __launch_bounds__(kStagedMaxThreads)
zstaged_kernel(const int* m, const int* n, const int* k, const double2* const* a, const int* lda,
               const double2* const* b, const int* ldb, double2* const* c, const int* ldc, int pitch, int pitch_b,
               int rows_pad) {
    extern __shared__ double2 zshared[];
    double2* sA = zshared;                                     // A^T(r, kk) at kk * pitch + r
    double2* sB = zshared + static_cast<size_t>(kKC) * pitch;  // B(kk, j) at kk * pitch_b + j
    const int e = blockIdx.x;
    const int M = m[e], N = n[e], K = k[e];
    if (M <= 0 || N <= 0) return;
    const double2* __restrict__ A = a[e];
    const double2* __restrict__ B = b[e];
    const int la = lda[e], lb = ldb[e];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5, warps = blockDim.x >> 5;
    const int r = tid % rows_pad, j0 = (tid / rows_pad) * kZGroupCols;
    const bool mine = r < M && j0 < N;
    double2* C = c[e] + r;
    const int lc = ldc[e];
    double2 cv[kAccum ? kZGroupCols : 1];
    if (kAccum && mine) {
#pragma unroll
        for (int j = 0; j < kZGroupCols; ++j) {
            cv[j] = j0 + j < N ? C[static_cast<int64_t>(j0 + j) * lc] : make_double2(0.0, 0.0);
        }
    }
    double2 acc[kZGroupCols];
#pragma unroll
    for (int j = 0; j < kZGroupCols; ++j) acc[j] = make_double2(0.0, 0.0);
    const int ncols = M + N;
    for (int k0 = 0; k0 < K; k0 += kKC) {
        const int kc = min(kKC, K - k0);
        __syncthreads();
        if (lane < kc) {
            for (int c0 = warp; c0 < ncols; c0 += warps * kLoads) {
                double2 v[kLoads];
#pragma unroll
                for (int u = 0; u < kLoads; ++u) {
                    const int col = c0 + u * warps;
                    const double2* src = col < M ? A + static_cast<int64_t>(col) * la
                                                 : B + static_cast<int64_t>(col - M) * lb;
                    v[u] = col < ncols ? __ldg(src + k0 + lane) : make_double2(0.0, 0.0);
                }
#pragma unroll
                for (int u = 0; u < kLoads; ++u) {
                    const int col = c0 + u * warps;
                    if (col < M) sA[lane * pitch + col] = v[u];
                    else if (col < ncols) sB[lane * pitch_b + col - M] = v[u];
                }
            }
        }
        __syncthreads();
        if (mine) {
            const double2* ap = sA + r;
            const double2* bp = sB + j0;
#pragma unroll 4
            for (int kk = 0; kk < kc; ++kk) {
                const double2 av = ap[kk * pitch];
#pragma unroll
                for (int q = 0; q < kZGroupCols; ++q) zfma(acc[q], av, bp[kk * pitch_b + q]);
            }
        }
    }
    if (!mine) return;
#pragma unroll
    for (int j = 0; j < kZGroupCols; ++j) {
        if (j0 + j >= N) continue;
        double2 v = acc[j];
        if (kAccum) {
            v.x += cv[j].x;
            v.y += cv[j].y;
        }
        C[static_cast<int64_t>(j0 + j) * lc] = v;
    }
}

}  // namespace

// Which narrow kernel takes a batch, if one does (bf_sweep_n.hpp): the
// direct one for one untransposed column, the staged one for transposed
// batches of up to kStagedMaxCols columns and kStagedMaxRows rows; MAGMA and
// the tensor cores take the rest.  Measured on one A100, 8192 GEMMs of
// 50 x 50 (ms, beta 0 / 1): at 16 columns the staged kernel 0.23 / 0.28, the
// tensor cores 0.28 / 0.37 and untransposed MAGMA 0.25 / 0.32; at 32
// columns 0.35 / 0.43, 0.34 / 0.54 and 0.47 / 0.57.  From 100 rows the
// tensor cores are the faster (160 x 160 at 16 columns: 0.42 against 0.47).
// Untransposed, MAGMA is the faster from 2 columns on.
int sweep_narrow_kind(bool trans_a, int max_m, int max_n) {
    if (!trans_a) return max_n == 1 ? 1 : 0;
    return (max_n <= kStagedMaxCols && max_m <= kStagedMaxRows) ? 2 : 0;
}

void launch_sweep_narrow(bool trans_a, int count, const int* m, const int* n, const int* k, DataType alpha,
                         const DataType* const* a, const int* lda, const DataType* const* b, const int* ldb,
                         DataType beta, DataType* const* c, const int* ldc, int max_m, int max_n,
                         cudaStream_t stream) {
    if (count <= 0) return;
    const int kind = sweep_narrow_kind(trans_a, max_m, max_n);
    if (kind == 1) {
        const int outputs = max_m * max_n;
        const int bs = std::min(256, (outputs + 31) / 32 * 32);
        for (int first = 0; first < count; first += 65535) {
            const int part = std::min(count - first, 65535);
            const dim3 grid((outputs + bs - 1) / bs, part);
            direct_kernel<<<grid, bs, 0, stream>>>(m + first, n + first, k + first, alpha, a + first, lda + first,
                                                   b + first, ldb + first, beta, c + first, ldc + first);
        }
    } else if (kind == 2) {
        staged_with(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, max_m, max_n, stream);
    } else {
        throw std::invalid_argument("launch_sweep_narrow: no narrow kernel takes this batch");
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in launch_sweep_narrow: ") +
                                 cudaGetErrorString(status));
    }
}

// The complex batches' narrow kernels: the rows-direct one for up to
// kZRowsMaxCols untransposed columns, the staged one for up to
// kZStagedMaxCols transposed; at most kStagedMaxRows rows either way.
// Untransposed from 9 columns MAGMA and the tensor cores are the faster.
int sweep_narrow_kind_z(bool trans_a, int max_m, int max_n) {
    if (max_m > kStagedMaxRows) return 0;
    if (!trans_a) return max_n <= kZRowsMaxCols ? 3 : 0;
    return max_n <= kZStagedMaxCols ? 4 : 0;
}

void launch_sweep_narrow(bool trans_a, int count, const int* m, const int* n, const int* k, fmm::gpu::dcomplex alpha,
                         const fmm::gpu::dcomplex* const* a, const int* lda, const fmm::gpu::dcomplex* const* b,
                         const int* ldb, fmm::gpu::dcomplex beta, fmm::gpu::dcomplex* const* c, const int* ldc,
                         int max_m, int max_n, cudaStream_t stream) {
    if (count <= 0) return;
    const bool accum = beta == fmm::gpu::dcomplex(1.0, 0.0);
    if (alpha != fmm::gpu::dcomplex(1.0, 0.0) || (!accum && beta != fmm::gpu::dcomplex(0.0, 0.0))) {
        throw std::invalid_argument("launch_sweep_narrow: complex batches take alpha 1 and beta 0 or 1");
    }
    const auto* za = reinterpret_cast<const double2* const*>(a);
    const auto* zb = reinterpret_cast<const double2* const*>(b);
    auto* zc = reinterpret_cast<double2* const*>(c);
    const int kind = sweep_narrow_kind_z(trans_a, max_m, max_n);
    const int rows_pad = (max_m + 31) / 32 * 32;
    const int groups = (max_n + kZGroupCols - 1) / kZGroupCols;
    const int threads = std::max(rows_pad * groups, 64);
    const int pitch_b = groups * kZGroupCols | 1;  // odd
    if (kind == 3) {
        const size_t shared = static_cast<size_t>(kKC) * pitch_b * sizeof(double2);
        zrows_kernel<<<count, threads, shared, stream>>>(m, n, k, za, lda, zb, ldb, accum, zc, ldc, rows_pad,
                                                         pitch_b);
    } else if (kind == 4) {
        // the widest block's shared memory, 32 x (65 + 17) complex values, is under 48 KB
        const int pitch = max_m | 1;
        const size_t shared = static_cast<size_t>(kKC) * (pitch + pitch_b) * sizeof(double2);
        if (accum) {
            zstaged_kernel<true><<<count, threads, shared, stream>>>(m, n, k, za, lda, zb, ldb, zc, ldc, pitch,
                                                                     pitch_b, rows_pad);
        } else {
            zstaged_kernel<false><<<count, threads, shared, stream>>>(m, n, k, za, lda, zb, ldb, zc, ldc, pitch,
                                                                      pitch_b, rows_pad);
        }
    } else {
        throw std::invalid_argument("launch_sweep_narrow: no narrow kernel takes this complex batch");
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in launch_sweep_narrow: ") +
                                 cudaGetErrorString(status));
    }
}

template<typename T>
void launch_sweep_patch(int count, const int64_t* b_off, const int64_t* c_off, const signed char* b_buf,
                        const signed char* c_buf, const T* vin, T* buf0, T* buf1, T* vout, int ldi, int ldo,
                        const T** b, T** c, int* ldb, int* ldc, cudaStream_t stream) {
    if (count <= 0) return;
    const int blocks = std::min((count + kThreads - 1) / kThreads, 1024);
    patch_kernel<T><<<blocks, kThreads, 0, stream>>>(count, b_off, c_off, b_buf, c_buf, vin, buf0, buf1, vout, ldi,
                                                  ldo, b, c, ldb, ldc);
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in launch_sweep_patch: ") +
                                 cudaGetErrorString(status));
    }
}

namespace {

// dst = conj(src), rows x cols, src and dst column major (in place when they
// are the same)
template<typename T>
__global__ void conj_kernel(T* dst, int64_t ldd, const T* src, int64_t lds, int64_t rows, int64_t cols) {
    using fmm::gpu::conj;
    const int64_t total = rows * cols;
    for (int64_t t = blockIdx.x * static_cast<int64_t>(blockDim.x) + threadIdx.x; t < total;
         t += static_cast<int64_t>(gridDim.x) * blockDim.x) {
        const int64_t r = t % rows, c = t / rows;
        dst[r + c * ldd] = conj(src[r + c * lds]);
    }
}

}  // namespace

template<typename T>
void launch_conj(T* dst, int64_t ldd, const T* src, int64_t lds, int64_t rows, int64_t cols, cudaStream_t stream) {
    if (rows <= 0 || cols <= 0) return;
    const int64_t blocks = std::min<int64_t>((rows * cols + kThreads - 1) / kThreads, 4096);
    conj_kernel<T><<<static_cast<int>(blocks), kThreads, 0, stream>>>(dst, ldd, src, lds, rows, cols);
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in launch_conj: ") + cudaGetErrorString(status));
    }
}

// the double and the double complex libraries
#define BF_PATCH_INSTANTIATE(T)                                                                                       \
    template void launch_sweep_patch<T>(int, const int64_t*, const int64_t*, const signed char*, const signed char*,  \
                                        const T*, T*, T*, T*, int, int, const T**, T**, int*, int*, cudaStream_t); \
    template void launch_conj<T>(T*, int64_t, const T*, int64_t, int64_t, int64_t, cudaStream_t);
BF_PATCH_INSTANTIATE(double)
BF_PATCH_INSTANTIATE(fmm::gpu::dcomplex)
#undef BF_PATCH_INSTANTIATE

}  // namespace gpu
}  // namespace bf
