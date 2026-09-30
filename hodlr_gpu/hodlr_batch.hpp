#pragma once
// Host descriptions of the batched operations of the HODLR GPU backend,
// staged in one pinned metadata image (fmm::gpu::MetaBuilder) and uploaded
// with a single copy per step: MAGMA GEMM batches (fmm::gpu::VBatch),
// triangular-solve and LU batches, and the item lists of hodlr_kernels.cu.

#ifdef H2_HAVE_GPU

#include "color_gpu/gpu_runtime.hpp"
#include "hodlr_kernels.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace bpack {
namespace gpu {

using fmm::gpu::MetaBuilder;

// Items of one launch of a hodlr_kernels.cu kernel.
template<typename Item>
struct ItemList {
    std::vector<Item> items;
    size_t offset = 0;
    int max_m = 0, max_n = 0;

    void add(const Item& item, int m, int n) {
        items.push_back(item);
        max_m = std::max(max_m, m);
        max_n = std::max(max_n, n);
    }
    int count() const { return static_cast<int>(items.size()); }
    void stage(MetaBuilder& meta) { offset = meta.append(items); }
    const Item* device(char* meta_device) const { return reinterpret_cast<const Item*>(meta_device + offset); }
};

// B_i = op(A_i)^{-1} B_i from the left, A_i m x m triangular, B_i m x n.
template<typename T>
struct TrsmBatch {
    struct Entry {
        T* a;
        T* b;
        magma_int_t m, n, lda, ldb;
    };
    std::vector<Entry> entries;
    size_t sizes_offset = 0, pointers_offset = 0;
    magma_int_t max_m = 0, max_n = 0;

    int count() const { return static_cast<int>(entries.size()); }
    void add(T* a, int lda, T* b, int ldb, int m, int n) {
        if (m <= 0 || n <= 0) return;
        entries.push_back({a, b, m, n, lda, ldb});
    }
    void stage(MetaBuilder& meta) {
        const size_t c = entries.size();
        std::vector<magma_int_t> sizes(4 * (c + 1), 0);
        std::vector<T*> pointers(2 * std::max<size_t>(c, 1), nullptr);
        max_m = max_n = 0;
        for (size_t i = 0; i < c; ++i) {
            const Entry& e = entries[i];
            max_m = std::max(max_m, e.m);
            max_n = std::max(max_n, e.n);
            sizes[0 * (c + 1) + i] = e.m;
            sizes[1 * (c + 1) + i] = e.n;
            sizes[2 * (c + 1) + i] = e.lda;
            sizes[3 * (c + 1) + i] = e.ldb;
            pointers[i] = e.a;
            pointers[c + i] = e.b;
        }
        sizes_offset = meta.append(sizes);
        pointers_offset = meta.append(pointers);
    }
    void solve(char* meta_device, magma_uplo_t uplo, magma_diag_t diag, magma_queue_t queue) const {
        const size_t c = entries.size();
        if (c == 0) return;
        magma_int_t* s = reinterpret_cast<magma_int_t*>(meta_device + sizes_offset);
        T** p = reinterpret_cast<T**>(meta_device + pointers_offset);
        fmm::gpu::trsm_vbatched<T>(MagmaLeft, uplo, MagmaNoTrans, diag, max_m, max_n, s, s + (c + 1), T(1.0), p,
                                   s + 2 * (c + 1), p + c, s + 3 * (c + 1), static_cast<magma_int_t>(c), queue);
    }
};

// LU with partial pivoting of square matrices; info and 1-based pivots on
// the device (info[i] > 0: exact zero pivot).
template<typename T>
struct GetrfBatch {
    struct Entry {
        T* a;
        magma_int_t n, lda;
        magma_int_t* ipiv;
    };
    std::vector<Entry> entries;
    size_t sizes_offset = 0, pointers_offset = 0;
    magma_int_t max_n = 0;

    int count() const { return static_cast<int>(entries.size()); }
    void add(T* a, int lda, int n, int* ipiv) { entries.push_back({a, n, lda, ipiv}); }
    void stage(MetaBuilder& meta) {
        const size_t c = entries.size();
        std::vector<magma_int_t> sizes(3 * (c + 1), 0);
        std::vector<void*> pointers(2 * std::max<size_t>(c, 1), nullptr);
        max_n = 0;
        for (size_t i = 0; i < c; ++i) {
            const Entry& e = entries[i];
            max_n = std::max(max_n, e.n);
            sizes[0 * (c + 1) + i] = e.n;
            sizes[1 * (c + 1) + i] = e.n;
            sizes[2 * (c + 1) + i] = e.lda;
            pointers[i] = e.a;
            pointers[c + i] = e.ipiv;
        }
        sizes_offset = meta.append(sizes);
        pointers_offset = meta.append(pointers);
    }
    // info: device array of count() entries
    void factor(char* meta_device, magma_int_t* info, fmm::gpu::DeviceBuffer& work, magma_queue_t queue) const {
        const size_t c = entries.size();
        if (c == 0) return;
        magma_int_t* s = reinterpret_cast<magma_int_t*>(meta_device + sizes_offset);
        void** p = reinterpret_cast<void**>(meta_device + pointers_offset);
        fmm::gpu::getrf_vbatched<T>(max_n, s, s + (c + 1), reinterpret_cast<T**>(p), s + 2 * (c + 1),
                                    reinterpret_cast<magma_int_t**>(p + c), info, static_cast<magma_int_t>(c), work,
                                    queue);
    }
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
