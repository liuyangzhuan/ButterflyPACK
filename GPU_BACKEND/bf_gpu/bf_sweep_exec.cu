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

// Buffers: 0 the caller's x, 1 and 2 the ping-pong buffers, 3 the caller's y.
__global__ void patch_kernel(int count, const int64_t* b_off, const int64_t* c_off, const signed char* b_buf,
                             const signed char* c_buf, const DataType* vin, DataType* buf0, DataType* buf1,
                             DataType* vout, int ldi, int ldo, const DataType** b, DataType** c, int* ldb,
                             int* ldc) {
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

void launch_sweep_patch(int count, const int64_t* b_off, const int64_t* c_off, const signed char* b_buf,
                        const signed char* c_buf, const DataType* vin, DataType* buf0, DataType* buf1,
                        DataType* vout, int ldi, int ldo, const DataType** b, DataType** c, int* ldb, int* ldc,
                        cudaStream_t stream) {
    if (count <= 0) return;
    const int blocks = std::min((count + kThreads - 1) / kThreads, 1024);
    patch_kernel<<<blocks, kThreads, 0, stream>>>(count, b_off, c_off, b_buf, c_buf, vin, buf0, buf1, vout, ldi,
                                                  ldo, b, c, ldb, ldc);
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in launch_sweep_patch: ") +
                                 cudaGetErrorString(status));
    }
}

}  // namespace gpu
}  // namespace bf
