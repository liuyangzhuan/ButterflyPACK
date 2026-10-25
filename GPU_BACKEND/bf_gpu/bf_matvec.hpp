#pragma once

#ifdef H2_HAVE_GPU

// Butterfly matvec on the device (GPU_BACKEND/bf_gpu).
//
// 'N': Vout = a * (U * K_L * .. * K_1 * V^T) * Vin + b * Vout, swept from V to
// U; 'T'/'C' sweeps the other way with the factors transposed.  A level's
// blocks are independent, so each level is one variable-size batched GEMM
// (fmm::gpu::gemm_vbatched) rather than a call per block: the deep index tree
// leaves many small blocks per level, which is the shape a batch wants and
// the shape that starved the CPU's taskloop.
//
// This header covers the single-process case: the per-level exchange and the
// all-to-all at level_half are where half the CPU matvec's time goes, and
// they are added next, on the device, so the vectors never return to the host.

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "bf_device.hpp"
#include "../h2_gpu/gpu_runtime.hpp"

namespace bf {
namespace gpu {

// Scratch for one sweep: the pointer and size arrays a batched GEMM reads,
// plus the two vector buffers it ping-pongs between.  Kept across calls
// because sample_Schur applies the same operators thousands of times (47432
// inverse calls at Poisson 400^3), and re-allocating per call would cost more
// than the arithmetic.
struct SweepScratch {
    DataType* vec[2] = {nullptr, nullptr};   // transient, from the arena
    int64_t vec_elements = 0;
    const DataType** d_a = nullptr;          // batch arrays, device side
    const DataType** d_b = nullptr;
    DataType** d_c = nullptr;
    magma_int_t *d_m = nullptr, *d_n = nullptr, *d_k = nullptr;
    magma_int_t *d_lda = nullptr, *d_ldb = nullptr, *d_ldc = nullptr;
    int64_t batch_capacity = 0;
    // host staging for those arrays (one copy per level, not per block)
    std::vector<const DataType*> h_a, h_b;
    std::vector<DataType*> h_c;
    std::vector<magma_int_t> h_m, h_n, h_k, h_lda, h_ldb, h_ldc;

    void reserve(int64_t batch, int64_t elements);
    void release();
    ~SweepScratch() { release(); }
};

// Vout(M_loc, nrhs) = a * op(bf) * Vin(N_loc, nrhs) + b * Vout, all device
// resident.  `op` is 'N', 'T' or 'C'.  Returns the flops performed, so a
// caller can account for them the way the CPU path does.
double matvec(const BfDevice& bf, char op, int nrhs,
              const DataType* vin, int ldi, DataType* vout, int ldo,
              DataType alpha, DataType beta, SweepScratch& scratch);

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
