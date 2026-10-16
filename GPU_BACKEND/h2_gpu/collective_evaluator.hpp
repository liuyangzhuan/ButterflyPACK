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
// With a sketch callback too, the ID's sketches of kernel rows are asked of
// the source sketched (collective_sketches), d x n for each box.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "evaluator.hpp"
#include "gpu_runtime.hpp"

#include <cstdint>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <vector>

namespace fmm {
namespace gpu {

template<typename DataType, typename Entries>
void collective_block_evaluator(const bpack_gpu_block_batch* batch, void* user) {
    Entries* entries = static_cast<Entries*>(user);
    cudaStream_t stream = static_cast<cudaStream_t>(batch->stream);
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    double fetching = 0.0;  // (the round: counted with the rounds)
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
        const auto t_fetch = clock::now();
        entries->fetch_now();
        fetching = std::chrono::duration<double>(clock::now() - t_fetch).count();
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
    entries->count_delivery(std::chrono::duration<double>(clock::now() - t0).count() - fetching);
}

// The sketches of a batch of kernel rows (Evaluator::sketch_rows), from the
// fetched ones, those not fetched yet fetched in one round on the spot, then
// into their slots on the device; each dropped once read.
template<typename DataType, typename Entries>
void collective_sketches(Entries* entries, std::vector<Evaluator::SketchRequest>& requests, cudaStream_t stream) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    double fetching = 0.0;
    std::vector<size_t> at(requests.size() + 1, 0);
    std::vector<typename Entries::SketchKey> keys;
    for (size_t q = 0; q < requests.size(); ++q) {
        const Evaluator::SketchRequest& r = requests[q];
        at[q + 1] = at[q] + static_cast<size_t>(r.d) * r.cols.size();
        keys.push_back(typename Entries::SketchKey{r.rows.data(), static_cast<int64_t>(r.rows.size()), r.cols.data(),
                                                   static_cast<int64_t>(r.cols.size()), r.d,
                                                   static_cast<int64_t>(r.orow.size()), r.orow.data(), r.odest.data(),
                                                   r.oval.data()});
    }
    if (at.back() == 0) return;
    std::vector<DataType> image(at.back());
    std::vector<size_t> missing;
    for (size_t q = 0; q < requests.size(); ++q) {
        if (!entries->lookup_sketch(keys[q], image.data() + at[q], requests[q].d)) {
            entries->record_sketch(keys[q]);
            missing.push_back(q);
        }
    }
    if (!missing.empty()) {
        const auto t_fetch = clock::now();
        entries->fetch_now();
        fetching = std::chrono::duration<double>(clock::now() - t_fetch).count();
        for (size_t q : missing) {
            if (!entries->lookup_sketch(keys[q], image.data() + at[q], requests[q].d)) {
                throw std::logic_error("collective GPU evaluator: a sketch its round did not fetch");
            }
        }
    }
    for (const auto& k : keys) entries->release_sketch(k);
    DeviceHeap& heap = DeviceHeap::instance();
    char* d_image = heap.alloc(image.size() * sizeof(DataType));
    check_cuda(cudaMemcpyAsync(d_image, image.data(), image.size() * sizeof(DataType), cudaMemcpyHostToDevice, stream),
               "collective GPU sketches upload");
    for (size_t q = 0; q < requests.size(); ++q) {
        const Evaluator::SketchRequest& r = requests[q];
        if (r.d <= 0 || r.cols.empty()) continue;
        const size_t col = static_cast<size_t>(r.d) * sizeof(DataType);
        check_cuda(cudaMemcpy2DAsync(r.out, static_cast<size_t>(r.ldo) * sizeof(DataType),
                                     d_image + at[q] * sizeof(DataType), col, col, r.cols.size(),
                                     cudaMemcpyDeviceToDevice, stream),
                   "collective GPU sketches copy");
    }
    heap.free(d_image);  // (later launches are ordered after the copies)
    entries->count_delivery(std::chrono::duration<double>(clock::now() - t0).count() - fetching);
}

// The evaluator over `entries` (a CollectiveEntries<DataType, ...>, which
// outlives it: both belong to the H2's kernel).
template<typename DataType, typename Entries>
std::shared_ptr<Evaluator> make_collective_evaluator(Entries* entries) {
    return Evaluator::from_collective(
        &collective_block_evaluator<DataType, Entries>, entries, static_cast<int>(sizeof(DataType)),
        [entries] { return entries->has_sketch(); },
        [entries](std::vector<Evaluator::SketchRequest>& requests, cudaStream_t stream) {
            collective_sketches<DataType>(entries, requests, stream);
        });
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
