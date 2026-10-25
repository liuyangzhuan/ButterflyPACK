#pragma once

#ifdef H2_HAVE_GPU

// The device butterfly matvec (declared in bf_matvec.hpp).

#include "bf_matvec.hpp"

namespace bf {
namespace gpu {

inline void SweepScratch::release() {
    auto& heap = fmm::gpu::DeviceHeap::instance();
    for (auto*& p : vec) { if (p) { heap.free(p); p = nullptr; } }
    for (void* p : {static_cast<void*>(d_a), static_cast<void*>(d_b),
                    static_cast<void*>(d_c), static_cast<void*>(d_m),
                    static_cast<void*>(d_n), static_cast<void*>(d_k),
                    static_cast<void*>(d_lda), static_cast<void*>(d_ldb),
                    static_cast<void*>(d_ldc)}) {
        if (p) heap.free(p);
    }
    d_a = d_b = nullptr; d_c = nullptr;
    d_m = d_n = d_k = d_lda = d_ldb = d_ldc = nullptr;
    vec_elements = 0; batch_capacity = 0;
}

inline void SweepScratch::reserve(int64_t batch, int64_t elements) {
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    // MAGMA's vbatched reads batch + 1 entries of each size array
    const int64_t nsz = batch + 1;
    if (batch > batch_capacity) {
        release();
        d_a = heap.alloc<const DataType*>(static_cast<size_t>(batch) * sizeof(DataType*),
                                          fmm::gpu::ResidentKind::Matvec);
        d_b = heap.alloc<const DataType*>(static_cast<size_t>(batch) * sizeof(DataType*),
                                          fmm::gpu::ResidentKind::Matvec);
        d_c = heap.alloc<DataType*>(static_cast<size_t>(batch) * sizeof(DataType*),
                                    fmm::gpu::ResidentKind::Matvec);
        auto ints = [&]() {
            return heap.alloc<magma_int_t>(static_cast<size_t>(nsz) * sizeof(magma_int_t),
                                           fmm::gpu::ResidentKind::Matvec);
        };
        d_m = ints(); d_n = ints(); d_k = ints();
        d_lda = ints(); d_ldb = ints(); d_ldc = ints();
        batch_capacity = batch;
    }
    if (elements > vec_elements) {
        for (auto*& p : vec) { if (p) { heap.free(p); p = nullptr; } }
        for (auto*& p : vec) {
            p = heap.alloc<DataType>(static_cast<size_t>(elements) * sizeof(DataType),
                                     fmm::gpu::ResidentKind::Matvec);
        }
        vec_elements = elements;
    }
    h_a.resize(static_cast<size_t>(batch)); h_b.resize(static_cast<size_t>(batch));
    h_c.resize(static_cast<size_t>(batch));
    h_m.resize(static_cast<size_t>(nsz)); h_n.resize(static_cast<size_t>(nsz));
    h_k.resize(static_cast<size_t>(nsz)); h_lda.resize(static_cast<size_t>(nsz));
    h_ldb.resize(static_cast<size_t>(nsz)); h_ldc.resize(static_cast<size_t>(nsz));
}

// One level's batch: C_i = op(A_i) * B_i, A_i the level's blocks, B_i and C_i
// slices of the ping-pong buffers.  Returns the flops.
inline double apply_level(const BfDevice& bf, int l, bool transposed, int nrhs,
                          const DataType* in, DataType* out,
                          const int* in_off, const int* out_off,
                          SweepScratch& s) {
    const LevelDesc& d = bf.level(l);
    const int64_t nblk = static_cast<int64_t>(d.nr) * d.nc;
    s.reserve(nblk, 0);
    magma_int_t max_m = 0, max_n = 0, max_k = 0;
    int64_t batch = 0;
    double flops = 0.0;
    for (int64_t b = 0; b < nblk; ++b) {
        const BlockDesc& bd = bf.blocks[static_cast<size_t>(d.first + b)];
        if (bd.rows == 0 || bd.cols == 0) continue;
        // 'N' applies A^T for V and A for the rest; a transposed sweep flips it
        const int m = transposed ? bd.cols : bd.rows;
        const int k = transposed ? bd.rows : bd.cols;
        s.h_a[static_cast<size_t>(batch)] = bf.values + bd.offset;
        s.h_lda[static_cast<size_t>(batch)] = bd.lda;
        s.h_b[static_cast<size_t>(batch)] = in + in_off[b];
        s.h_ldb[static_cast<size_t>(batch)] = k;
        s.h_c[static_cast<size_t>(batch)] = out + out_off[b];
        s.h_ldc[static_cast<size_t>(batch)] = m;
        s.h_m[static_cast<size_t>(batch)] = m;
        s.h_n[static_cast<size_t>(batch)] = nrhs;
        s.h_k[static_cast<size_t>(batch)] = k;
        max_m = std::max(max_m, static_cast<magma_int_t>(m));
        max_n = std::max(max_n, static_cast<magma_int_t>(nrhs));
        max_k = std::max(max_k, static_cast<magma_int_t>(k));
        flops += 2.0 * m * nrhs * k;
        ++batch;
    }
    if (batch == 0) return 0.0;
    auto stream = fmm::gpu::Context::instance().stream();
    auto put = [&](const void* src, void* dst, size_t bytes) {
        fmm::gpu::check_cuda(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, stream),
                             "bf::gpu::apply_level (batch arrays)");
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
        transposed ? MagmaTrans : MagmaNoTrans, MagmaNoTrans,
        s.d_m, s.d_n, s.d_k, DataType(1), s.d_a, s.d_lda, s.d_b, s.d_ldb,
        DataType(0), s.d_c, s.d_ldc, static_cast<magma_int_t>(batch),
        fmm::gpu::Context::instance().queue(), max_m, max_n, max_k);
    return flops;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
