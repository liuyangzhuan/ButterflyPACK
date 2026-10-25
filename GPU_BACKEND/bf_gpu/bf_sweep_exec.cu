// The one kernel a recorded sweep runs per call (bf_sweep_n.hpp, SweepExec):
// its GEMMs' vector pointers, from offsets into the caller's x and y and the
// two ping-pong buffers, and the leading dimensions of x and y.

#include "bf_device.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>

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

// C = alpha A^T B + beta C for at most kCols columns: one thread block per
// entry and a thread per row of C, k in chunks of 32 staged in shared memory
// (a warp per column of A, its lanes along k, stored with an odd pitch so
// neither the stores nor the reads collide on a bank).  Each element of C
// accumulates its k products in order from zero, one rounding each, then
// takes alpha and beta as MAGMA's kernels do - so the result is MAGMA's, and
// the tensor-core GEMM's, to the bit.  A batch of few columns streams the
// factors once instead of through GEMM tiles that a column or two leave
// mostly empty.
constexpr int kKC = 32;
constexpr int kSkinnyMaxRows = 160;  // the pitch keeps a block's shared memory under 48 KB

template<int kCols, int kThreads>
__global__ void __launch_bounds__(kThreads)
skinny_kernel(const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a, const int* lda,
              const DataType* const* b, const int* ldb, DataType beta, DataType* const* c, const int* ldc,
              int pitch) {
    constexpr int kWarps = kThreads / 32;
    extern __shared__ DataType smem[];
    DataType* sA = smem;                                       // A^T(r, kk) at kk * pitch + r
    DataType* sB = smem + static_cast<size_t>(kKC) * pitch;    // B(kk, j) at j * kKC + kk
    const int e = blockIdx.x;
    const int M = m[e], N = n[e], K = k[e];
    if (M <= 0 || N <= 0) return;
    const DataType* A = a[e];
    const DataType* B = b[e];
    const int la = lda[e], lb = ldb[e];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    DataType acc[kCols];
#pragma unroll
    for (int j = 0; j < kCols; ++j) acc[j] = DataType(0);
    for (int k0 = 0; k0 < K; k0 += kKC) {
        const int kc = min(kKC, K - k0);
        __syncthreads();
        // A^T(r, kk) = A(kk, r) at A[kk + r la]: a warp per r, lanes along kk
        for (int r = warp; r < M; r += kWarps) {
            if (lane < kc) sA[lane * pitch + r] = A[(k0 + lane) + static_cast<int64_t>(r) * la];
        }
        for (int j = warp; j < N; j += kWarps) {
            if (lane < kc) sB[j * kKC + lane] = B[(k0 + lane) + static_cast<int64_t>(j) * lb];
        }
        __syncthreads();
        if (tid < M) {
            for (int kk = 0; kk < kc; ++kk) {
                const DataType av = sA[kk * pitch + tid];
#pragma unroll
                for (int j = 0; j < kCols; ++j) {
                    if (j < N) acc[j] = fma(av, sB[j * kKC + kk], acc[j]);
                }
            }
        }
    }
    if (tid >= M) return;
    DataType* C = c[e] + tid;
    const int lc = ldc[e];
#pragma unroll
    for (int j = 0; j < kCols; ++j) {
        if (j >= N) continue;
        DataType* cp = C + static_cast<int64_t>(j) * lc;
        *cp = beta == DataType(0) ? alpha * acc[j] : alpha * acc[j] + beta * *cp;
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

template<int kCols, int kThreads>
void skinny_threads(int count, const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a,
                    const int* lda, const DataType* const* b, const int* ldb, DataType beta, DataType* const* c,
                    const int* ldc, int pitch, cudaStream_t stream) {
    const size_t shared = (static_cast<size_t>(kKC) * pitch + static_cast<size_t>(kCols) * kKC) * sizeof(DataType);
    skinny_kernel<kCols, kThreads><<<count, kThreads, shared, stream>>>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc,
                                                                       pitch);
}

template<int kCols>
void skinny_with(int count, const int* m, const int* n, const int* k, DataType alpha, const DataType* const* a,
                 const int* lda, const DataType* const* b, const int* ldb, DataType beta, DataType* const* c,
                 const int* ldc, int max_m, cudaStream_t stream) {
    const int pitch = max_m | 1;  // odd
    if (max_m <= 64) {
        skinny_threads<kCols, 64>(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, pitch, stream);
    } else if (max_m <= 128) {
        skinny_threads<kCols, 128>(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, pitch, stream);
    } else {
        skinny_threads<kCols, 192>(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, pitch, stream);
    }
}

}  // namespace

// Which narrow kernel takes a batch, if one does (bf_sweep_n.hpp): the
// direct one for one untransposed column, the staged one for up to 8
// transposed columns; MAGMA and the tensor cores take the rest.  Measured on
// one A100 (262144^2, L = 12): MAGMA's transposed batches run at 4.6 ms a
// matvec whatever the column count up to 16, the staged kernel 2.3-2.9 ms up
// to 8 columns; untransposed, MAGMA is the faster from 2 columns on.
int sweep_narrow_kind(bool trans_a, int max_m, int max_n) {
    if (!trans_a) return max_n == 1 ? 1 : 0;
    return (max_n <= 8 && max_m <= kSkinnyMaxRows) ? 2 : 0;
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
        if (max_n <= 4) skinny_with<4>(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, max_m, stream);
        else skinny_with<8>(count, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, max_m, stream);
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
