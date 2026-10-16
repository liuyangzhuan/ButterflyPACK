#pragma once
// The device evaluator of entries through a collective callback
// (elem_extract=1 with H2_use_gpu; h2_parallel/core/collective_entries.hpp):
// a block evaluator that reads the blocks the rank has fetched; those it has
// not are fetched in one round on the spot (CollectiveEntries::fetch_now),
// which the other ranks join at their own next round (a fetch of theirs, a
// quiesce, a phase of the host code, or serving).  Each block goes to the
// device once and is then dropped from the fetched ones.  It replaces the
// application's evaluator when the H2 is built from such a callback
// (C_BPACK_wrapper.cpp) and is never warmed up (that would be a fetch).

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "evaluator.hpp"
#include "gpu_runtime.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace fmm {
namespace gpu {

template<typename DataType, typename Entries>
void collective_block_evaluator(const bpack_gpu_block_batch* batch, void* user) {
    Entries* entries = static_cast<Entries*>(user);
    cudaStream_t stream = static_cast<cudaStream_t>(batch->stream);
    const int count = batch->count;
    std::vector<size_t> at(static_cast<size_t>(count) + 1, 0);
    for (int b = 0; b < count; ++b) {
        const bpack_gpu_block& blk = batch->blocks[b];
        at[static_cast<size_t>(b) + 1] = at[static_cast<size_t>(b)] + static_cast<size_t>(blk.m) * blk.n;
    }
    if (at.back() == 0) return;
    // the blocks packed column by column, those not fetched yet fetched now
    std::vector<DataType> image(at.back());
    std::vector<int> missing;
    for (int b = 0; b < count; ++b) {
        const bpack_gpu_block& blk = batch->blocks[b];
        if (blk.m <= 0 || blk.n <= 0) continue;
        if (!entries->lookup(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n, image.data() + at[static_cast<size_t>(b)],
                             blk.m)) {
            entries->record(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n);
            missing.push_back(b);
        }
    }
    if (!missing.empty()) {
        entries->fetch_now();
        for (int b : missing) {
            const bpack_gpu_block& blk = batch->blocks[b];
            if (!entries->lookup(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n,
                                 image.data() + at[static_cast<size_t>(b)], blk.m)) {
                throw std::logic_error("collective GPU evaluator: a block its round did not fetch");
            }
        }
    }
    for (int b = 0; b < count; ++b) {
        const bpack_gpu_block& blk = batch->blocks[b];
        if (blk.m > 0 && blk.n > 0) entries->release(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n);
    }
    // one upload (pageable: staged before the call returns), then into the
    // blocks on the batch's stream
    const size_t bytes = image.size() * sizeof(DataType);
    char* d_image = static_cast<char*>(batch->allocate(batch->allocator, bytes));
    if (d_image == nullptr) throw std::runtime_error("collective GPU evaluator: no device memory for a batch");
    check_cuda(cudaMemcpyAsync(d_image, image.data(), bytes, cudaMemcpyHostToDevice, stream),
               "collective GPU evaluator upload");
    for (int b = 0; b < count; ++b) {
        const bpack_gpu_block& blk = batch->blocks[b];
        if (blk.m <= 0 || blk.n <= 0) continue;
        const size_t col = static_cast<size_t>(blk.m) * sizeof(DataType);
        check_cuda(cudaMemcpy2DAsync(blk.out, static_cast<size_t>(blk.ld) * sizeof(DataType),
                                     d_image + at[static_cast<size_t>(b)] * sizeof(DataType), col, col,
                                     static_cast<size_t>(blk.n), cudaMemcpyDeviceToDevice, stream),
                   "collective GPU evaluator copy");
    }
    batch->release(batch->allocator, d_image);  // (later launches are ordered after the copies)
}

// The evaluator over `entries` (a CollectiveEntries<DataType, ...>, which
// outlives it: both belong to the H2's kernel).
template<typename DataType, typename Entries>
std::shared_ptr<Evaluator> make_collective_evaluator(Entries* entries) {
    return Evaluator::from_collective(&collective_block_evaluator<DataType, Entries>, entries,
                                      static_cast<int>(sizeof(DataType)));
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
