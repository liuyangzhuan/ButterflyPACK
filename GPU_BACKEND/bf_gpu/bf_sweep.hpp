#pragma once

#ifdef H2_HAVE_GPU

// The 'N' sweep of a butterfly on the device, one process
// (GPU_BACKEND/bf_gpu).  Mirrors BF_block_MVP_dat_nonbatch in
// SRC/Bplus_utilities.f90; the multi-rank exchange and the all-to-all at
// level_half are not here yet.
//
// The sweep has two halves that meet at level_half, because the butterfly's
// row-wise and column-wise orderings meet there:
//
//   levels 0 .. level_half            GATHER: one output block is the sum of
//                                     two input blocks, through two kernel
//                                     blocks in the same row - two GEMMs into
//                                     the same C, beta 0 then beta 1.
//   levels level_half+1 .. L+1        SCATTER: one input block feeds two
//                                     output blocks, through two kernel
//                                     blocks in the same column - two GEMMs
//                                     with different C.
//
// Level 0 applies V^T to slices of the input (partitioned by the accumulated
// local column sizes, arr_acc_n on the host), and level L+1 applies U into
// slices of the output with the caller's alpha and beta.  Everything between
// is kernel blocks only.
//
// Each level is issued as one or two variable-size batched GEMMs, since a
// level's blocks are independent.  That is the point of the port: the deep
// index tree leaves many small blocks per level, which a batch likes and
// which starved the CPU's taskloop.

#include <algorithm>
#include <vector>

#include "bf_matvec_impl.hpp"
#include "bf_vec.hpp"

namespace bf {
namespace gpu {

// One entry of a level's batch.
struct GemmItem {
    const DataType* a;
    const DataType* b;
    DataType* c;
    int m, n, k, lda, ldb, ldc;
};

// Issue a batch: C = alpha * op(A) * B + beta * C for every entry.  MAGMA
// takes one alpha and one beta per call, so a level that accumulates two
// products into the same C issues twice, beta 0 then beta 1, and only the
// last level of a sweep carries the caller's alpha.
inline double issue(std::vector<GemmItem>& items, DataType alpha, DataType beta,
                    bool trans_a, SweepScratch& s) {
    if (items.empty()) return 0.0;
    const int64_t batch = static_cast<int64_t>(items.size());
    s.reserve(batch, 0);
    magma_int_t max_m = 0, max_n = 0, max_k = 0;
    double flops = 0.0;
    for (int64_t i = 0; i < batch; ++i) {
        const GemmItem& it = items[static_cast<size_t>(i)];
        s.h_a[static_cast<size_t>(i)] = it.a;
        s.h_b[static_cast<size_t>(i)] = it.b;
        s.h_c[static_cast<size_t>(i)] = it.c;
        s.h_m[static_cast<size_t>(i)] = it.m;
        s.h_n[static_cast<size_t>(i)] = it.n;
        s.h_k[static_cast<size_t>(i)] = it.k;
        s.h_lda[static_cast<size_t>(i)] = it.lda;
        s.h_ldb[static_cast<size_t>(i)] = it.ldb;
        s.h_ldc[static_cast<size_t>(i)] = it.ldc;
        max_m = std::max(max_m, static_cast<magma_int_t>(it.m));
        max_n = std::max(max_n, static_cast<magma_int_t>(it.n));
        max_k = std::max(max_k, static_cast<magma_int_t>(it.k));
        flops += 2.0 * it.m * it.n * it.k;
    }
    auto stream = fmm::gpu::Context::instance().stream();
    auto put = [&](const void* src, void* dst, size_t bytes) {
        fmm::gpu::check_cuda(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream),
                             "bf::gpu::issue");
    };
    put(s.h_a.data(), s.d_a, static_cast<size_t>(batch) * sizeof(DataType*));
    put(s.h_b.data(), s.d_b, static_cast<size_t>(batch) * sizeof(DataType*));
    put(s.h_c.data(), s.d_c, static_cast<size_t>(batch) * sizeof(DataType*));
    put(s.h_m.data(), s.d_m, static_cast<size_t>(batch) * sizeof(magma_int_t));
    put(s.h_n.data(), s.d_n, static_cast<size_t>(batch) * sizeof(magma_int_t));
    put(s.h_k.data(), s.d_k, static_cast<size_t>(batch) * sizeof(magma_int_t));
    put(s.h_lda.data(), s.d_lda, static_cast<size_t>(batch) * sizeof(magma_int_t));
    put(s.h_ldb.data(), s.d_ldb, static_cast<size_t>(batch) * sizeof(magma_int_t));
    put(s.h_ldc.data(), s.d_ldc, static_cast<size_t>(batch) * sizeof(magma_int_t));
    fmm::gpu::gemm_vbatched<DataType>(
        trans_a ? MagmaTrans : MagmaNoTrans, MagmaNoTrans,
        s.d_m, s.d_n, s.d_k, alpha, s.d_a, s.d_lda, s.d_b, s.d_ldb,
        beta, s.d_c, s.d_ldc, static_cast<magma_int_t>(batch),
        fmm::gpu::Context::instance().queue(), max_m, max_n, max_k);
    items.clear();
    return flops;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
