// Batched DGEMM of variable sizes on the FP64 tensor cores (sm_80 DMMA,
// mma.sync m8n8k4) for H2_use_gpu=2.  See launch_dgemm_vbatched_tc in
// device_kernels.hpp.
//
// A block computes a BM x BM tile of one C_i with WM x WN warps, each warp a
// (BM / WM) x (BM / WN) part as 8 x 8 DMMA tiles.  K advances 16 at a time
// through kStages shared buffers filled by cp.async, so several slices are in
// flight while one is multiplied.  Both operands are kept in shared memory as
// [row][k] (A as op(A)[i][k], B as op(B)^T[j][k]) with a padded pitch, so
// every fragment load of a half-warp hits distinct banks.  The host lists the
// (matrix, tile) of every block, the tiles of a matrix together (so its
// panels are reread from L2) and no block without work.
//
// Each element of C accumulates its k products in order with one rounding
// each (DMMA is IEEE fp64), as MAGMA's kernels do: the results are the same.

#include "device_kernels.hpp"

#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

constexpr int kBK = 16;
constexpr int kStages = 3;
constexpr int kPitch = kBK + 4;  // doubles per shared row

// d += a * b on one 8 x 8 x 4 product (fragment layouts of mma.m8n8k4.f64:
// a = A[lane / 4][lane % 4], b = B[lane % 4][lane / 4],
// d = D[lane / 4][2 (lane % 4) + {0, 1}]).
__device__ __forceinline__ void dmma(double (&d)[2], double a, double b) {
    asm("mma.sync.aligned.m8n8k4.row.col.f64.f64.f64.f64 {%0, %1}, {%2}, {%3}, {%0, %1};\n"
        : "+d"(d[0]), "+d"(d[1])
        : "d"(a), "d"(b));
}

// 8 bytes global -> shared, zero-filled when !valid (src is then not read).
__device__ __forceinline__ void cp_async8(unsigned dst, const double* src, bool valid) {
    asm volatile("cp.async.ca.shared.global [%0], [%1], 8, %2;\n" ::"r"(dst), "l"(src), "r"(valid ? 8 : 0));
}
__device__ __forceinline__ void cp_async_commit() { asm volatile("cp.async.commit_group;\n" ::); }
template<int kPending>
__device__ __forceinline__ void cp_async_wait() {
    asm volatile("cp.async.wait_group %0;\n" ::"n"(kPending));
}

template<int BM, int WM, int WN>
struct TileConfig {
    static constexpr int kTile = BM;
    static constexpr int kThreads = 32 * WM * WN;
    static constexpr int kWarpsN = WN;
    static constexpr int kWarpRows = BM / WM, kWarpCols = BM / WN;
    static constexpr int kTI = kWarpRows / 8, kTJ = kWarpCols / 8;  // DMMA tiles per warp
    static constexpr int kLoads = BM * kBK / kThreads;               // per operand, thread and slice
    static constexpr int kSliceDoubles = BM * kPitch;
    static constexpr size_t kSharedBytes = 2 * kStages * kSliceDoubles * sizeof(double);
    static_assert(kWarpRows % 8 == 0 && kWarpCols % 8 == 0, "warp tiles are made of 8 x 8 DMMA tiles");
    static_assert(kThreads % BM == 0 && kThreads % kBK == 0 && (BM * kBK) % kThreads == 0,
                  "the slice loaders need whole rows of threads");
};

// One operand's part of every BM x 16 slice, in slice order.  Element q of a
// thread, along the operand's contiguous direction so each warp reads whole
// lines (T threads):
//   rows contiguous (op(A) = A, op(B)^T = B^T): row tid % BM, k tid / BM + (T / BM) q,
//                                                 at X[row + k ld]
//   k contiguous    (op(A) = A^T, op(B)^T = B): k tid % 16, row tid / 16 + (T / 16) q,
//                                                 at X[k + row ld]
// Elements outside rows x k are zero.
template<class C, bool kRowsContiguous>
struct SliceLoader {
    static constexpr int kQRows = C::kThreads / kBK;    // row step of q (k contiguous)
    static constexpr int kQK = C::kThreads / C::kTile;  // k step of q (rows contiguous)
    const double* base;  // the matrix (a valid address for masked copies)
    const double* next;  // element q = 0 of the next slice
    int64_t q_step;      // between elements q and q + 1
    int64_t slice_step;  // between slices
    unsigned dst0;       // shared byte offset of element q = 0 in a slice
    int k_first;         // k of element q = 0 within a slice
    unsigned row_ok;     // bit q: element q's row is inside the matrix

    __device__ SliceLoader(const double* x, int ld, int rows, int row0, int tid) : base(x) {
        int row, kk;
        if (kRowsContiguous) {
            row = tid % C::kTile;
            kk = tid / C::kTile;
            q_step = kQK * static_cast<int64_t>(ld);
            slice_step = kBK * static_cast<int64_t>(ld);
            next = x + (row0 + row) + static_cast<int64_t>(kk) * ld;
            row_ok = row0 + row < rows ? 0xffffu : 0u;
        } else {
            kk = tid % kBK;
            row = tid / kBK;
            q_step = kQRows * static_cast<int64_t>(ld);
            slice_step = kBK;
            next = x + kk + static_cast<int64_t>(row0 + row) * ld;
            row_ok = 0;
            for (int q = 0; q < C::kLoads; ++q) {
                if (row0 + row + kQRows * q < rows) row_ok |= 1u << q;
            }
        }
        dst0 = static_cast<unsigned>((row * kPitch + kk) * sizeof(double));
        k_first = kk;
    }

    // Copy the next slice (starting at k0) into the slice buffer at `smem`.
    __device__ __forceinline__ void copy(unsigned smem, int k0, int k) {
#pragma unroll
        for (int q = 0; q < C::kLoads; ++q) {
            const int gk = k0 + k_first + (kRowsContiguous ? kQK * q : 0);
            const bool valid = ((row_ok >> q) & 1u) && gk < k;
            const unsigned dst = smem + dst0 +
                                 static_cast<unsigned>((kRowsContiguous ? kQK * q : kQRows * q * kPitch) *
                                                       sizeof(double));
            cp_async8(dst, valid ? next + q * q_step : base, valid);
        }
        next += slice_step;
    }
};

// C_i = alpha op(A_i) op(B_i) + beta C_i (C_i not read when beta = 0), block
// x on the tile blocks[x] = {i, tile row << 16 | tile column}.
template<class C, bool kTransA, bool kTransB>
__global__ void __launch_bounds__(C::kThreads)
dgemm_vbatched_tc_kernel(const int* __restrict__ ms, const int* __restrict__ ns, const int* __restrict__ ks,
                         double alpha, const double* const* __restrict__ as, const int* __restrict__ ldas,
                         const double* const* __restrict__ bs, const int* __restrict__ ldbs, double beta,
                         double* const* __restrict__ cs, const int* __restrict__ ldcs,
                         const int2* __restrict__ blocks) {
    constexpr int BM = C::kTile, TI = C::kTI, TJ = C::kTJ;
    extern __shared__ double shared[];  // stage s: A slice at 2 s, B slice at 2 s + 1
    const int2 block = blocks[blockIdx.x];
    const int item = block.x;
    const int m = ms[item], n = ns[item], k = ks[item];
    const int i0 = (block.y >> 16) * BM, j0 = (block.y & 0xffff) * BM;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const int wm = (warp / C::kWarpsN) * C::kWarpRows, wn = (warp % C::kWarpsN) * C::kWarpCols;
    // op(A) has contiguous rows when A is not transposed; op(B)^T has
    // contiguous rows (j) when B is transposed.
    SliceLoader<C, !kTransA> load_a(as[item], ldas[item], m, i0, tid);
    SliceLoader<C, kTransB> load_b(bs[item], ldbs[item], n, j0, tid);
    const unsigned smem = static_cast<unsigned>(__cvta_generic_to_shared(shared));
    constexpr unsigned kSliceBytes = C::kSliceDoubles * sizeof(double);

    const int slices = (k + kBK - 1) / kBK;
    int issued = 0;
    auto issue = [&] {  // the next slice, or an empty group past the end
        if (issued < slices) {
            const unsigned stage = smem + static_cast<unsigned>(issued % kStages) * 2 * kSliceBytes;
            load_a.copy(stage, issued * kBK, k);
            load_b.copy(stage + kSliceBytes, issued * kBK, k);
        }
        ++issued;
        cp_async_commit();
    };
#pragma unroll
    for (int s = 0; s < kStages - 1; ++s) issue();

    double acc[TI][TJ][2];
#pragma unroll
    for (int ti = 0; ti < TI; ++ti)
#pragma unroll
        for (int tj = 0; tj < TJ; ++tj) acc[ti][tj][0] = acc[ti][tj][1] = 0.0;

    const int frag_row = lane >> 2, frag_k = lane & 3;
    for (int t = 0; t < slices; ++t) {
        cp_async_wait<kStages - 2>();  // slice t has landed (this thread's copies)
        __syncthreads();                // everyone's copies; slice t - 1 is no longer read
        issue();                        // into the buffer of slice t - 1
        const double* sa = shared + (t % kStages) * 2 * C::kSliceDoubles;
        const double* pa = sa + (wm + frag_row) * kPitch + frag_k;
        const double* pb = sa + C::kSliceDoubles + (wn + frag_row) * kPitch + frag_k;
#pragma unroll
        for (int kk = 0; kk < kBK; kk += 4) {
            double fa[TI], fb[TJ];
#pragma unroll
            for (int u = 0; u < TI; ++u) fa[u] = pa[u * 8 * kPitch + kk];
#pragma unroll
            for (int u = 0; u < TJ; ++u) fb[u] = pb[u * 8 * kPitch + kk];
#pragma unroll
            for (int ti = 0; ti < TI; ++ti)
#pragma unroll
                for (int tj = 0; tj < TJ; ++tj) dmma(acc[ti][tj], fa[ti], fb[tj]);
        }
    }
    cp_async_wait<0>();

    // C by rows of DMMA tiles: the old values of a thread's row are all loaded
    // before any is stored, so the loads overlap
    double* __restrict__ c = cs[item];
    const int ldc = ldcs[item];
#pragma unroll
    for (int ti = 0; ti < TI; ++ti) {
        const int row = i0 + wm + ti * 8 + frag_row;
        if (row >= m) continue;
        double* __restrict__ crow = c + row;
        double old[TJ][2];
#pragma unroll
        for (int tj = 0; tj < TJ; ++tj)
#pragma unroll
            for (int e = 0; e < 2; ++e) {
                const int col = j0 + wn + tj * 8 + 2 * frag_k + e;
                old[tj][e] = beta != 0.0 && col < n ? crow[static_cast<int64_t>(col) * ldc] : 0.0;
            }
#pragma unroll
        for (int tj = 0; tj < TJ; ++tj)
#pragma unroll
            for (int e = 0; e < 2; ++e) {
                const int col = j0 + wn + tj * 8 + 2 * frag_k + e;
                if (col >= n) continue;
                crow[static_cast<int64_t>(col) * ldc] =
                    beta == 0.0 ? alpha * acc[ti][tj][e] : fma(alpha, acc[ti][tj][e], beta * old[tj][e]);
            }
    }
}

template<class C, bool kTransA, bool kTransB>
void launch_tc(const int* m, const int* n, const int* k, double alpha, const double* const* a, const int* lda,
               const double* const* b, const int* ldb, double beta, double* const* c, const int* ldc,
               const int2* blocks, int num_blocks, cudaStream_t stream) {
    auto kernel = dgemm_vbatched_tc_kernel<C, kTransA, kTransB>;
    static const cudaError_t attr = cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                                         static_cast<int>(C::kSharedBytes));
    if (attr != cudaSuccess) {
        throw std::runtime_error(std::string("launch_dgemm_vbatched_tc: ") + cudaGetErrorString(attr));
    }
    kernel<<<static_cast<unsigned>(num_blocks), C::kThreads, C::kSharedBytes, stream>>>(m, n, k, alpha, a, lda, b,
                                                                                       ldb, beta, c, ldc, blocks);
}

template<class C>
void launch_tc_config(bool trans_a, bool trans_b, const int* m, const int* n, const int* k, double alpha,
                      const double* const* a, const int* lda, const double* const* b, const int* ldb, double beta,
                      double* const* c, const int* ldc, const int2* blocks, int num_blocks, cudaStream_t stream) {
    if (!trans_a && !trans_b) {
        launch_tc<C, false, false>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, blocks, num_blocks, stream);
    } else if (!trans_a && trans_b) {
        launch_tc<C, false, true>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, blocks, num_blocks, stream);
    } else if (trans_a && !trans_b) {
        launch_tc<C, true, false>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, blocks, num_blocks, stream);
    } else {
        launch_tc<C, true, true>(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, blocks, num_blocks, stream);
    }
}

}  // namespace

void launch_dgemm_vbatched_tc(bool trans_a, bool trans_b, const int* m, const int* n, const int* k, double alpha,
                              const double* const* a, const int* lda, const double* const* b, const int* ldb,
                              double beta, double* const* c, const int* ldc, const int2* blocks, int num_blocks,
                              int tile, cudaStream_t stream) {
    if (num_blocks <= 0) return;
    switch (tile) {
        case 32:
            launch_tc_config<TileConfig<32, 2, 2>>(trans_a, trans_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc,
                                                   blocks, num_blocks, stream);
            break;
        case 48:
            launch_tc_config<TileConfig<48, 3, 2>>(trans_a, trans_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc,
                                                   blocks, num_blocks, stream);
            break;
        case 64:
            launch_tc_config<TileConfig<64, 2, 2>>(trans_a, trans_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc,
                                                   blocks, num_blocks, stream);
            break;
        default:
            throw std::invalid_argument("launch_dgemm_vbatched_tc: tile must be 32, 48 or 64");
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in dgemm_vbatched_tc_kernel: ") +
                                 cudaGetErrorString(status));
    }
}

}  // namespace gpu
}  // namespace fmm
