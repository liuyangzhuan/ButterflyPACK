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
// With a sketch callback too (and its transposed one, for a nonsymmetric
// matrix), the ID's sketches of kernel rows are asked of
// the source sketched (collective_sketches), d x n for each box.

#ifdef H2_HAVE_GPU

#include "core/tree.hpp"
#include "device_heap.hpp"
#include "evaluator.hpp"
#include "gpu_runtime.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <vector>

namespace fmm {
namespace gpu {

// The host image of a batch's values, uploaded in one copy: one of the pool's
// two pinned image buffers in turn, taken once the upload from it two images
// ago is done.  The callback writes the values into pinned memory and the
// upload is one DMA; from pageable memory the driver staged it on the calling
// thread (~1 s of a critical rank's 23 s of sequential fronts, Poisson 400^3,
// 64 ranks).  An image taken while another is held (by a nested batch or
// another thread; none is expected) is heap memory, staged as before.
class CollectiveImage {
public:
    explicit CollectiveImage(size_t bytes) {
        State& s = state();
        bool idle = false;
        if (!s.busy.compare_exchange_strong(idle, true)) {
            heap_.resize(bytes);
            data_ = heap_.data();
            return;
        }
        turn_ = s.turn;
        s.turn ^= 1;
        check_cuda(cudaEventSynchronize(s.done[turn_]), "collective image");
        data_ = static_cast<char*>(pinned_pool().image[turn_].reserve(std::max<size_t>(bytes, 1)));
    }
    ~CollectiveImage() {
        if (turn_ >= 0) state().busy = false;
    }
    CollectiveImage(const CollectiveImage&) = delete;
    CollectiveImage& operator=(const CollectiveImage&) = delete;
    template<typename T>
    T* as() const { return reinterpret_cast<T*>(data_); }
    // bytes of it to dst on the stream
    void upload(void* dst, size_t bytes, cudaStream_t stream) {
        check_cuda(cudaMemcpyAsync(dst, data_, bytes, cudaMemcpyHostToDevice, stream), "collective image upload");
        // (from heap memory the copy is staged before the call returns)
        if (turn_ >= 0) check_cuda(cudaEventRecord(state().done[turn_], stream), "cudaEventRecord");
    }

private:
    struct State {
        cudaEvent_t done[2] = {nullptr, nullptr};
        int turn = 0;
        std::atomic<bool> busy{false};
    };
    // (never freed, as the pool: the CUDA runtime's own shutdown would race it)
    static State& state() {
        static State* s = [] {
            State* p = new State;
            for (cudaEvent_t& e : p->done) check_cuda(cudaEventCreateWithFlags(&e, cudaEventDisableTiming), "cudaEventCreate");
            return p;
        }();
        return *s;
    }
    char* data_ = nullptr;
    int turn_ = -1;
    HostBuffer<char> heap_;
};

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
    // (one rank: all of them, straight into the image, zeroed for the
    // callback)
    CollectiveImage img(at.back() * sizeof(DataType));
    DataType* const image = img.as<DataType>();
    if (entries->alone()) {
        parallel_zero(image, at.back());
        std::vector<typename Entries::BlockRef> refs;
        refs.reserve(static_cast<size_t>(count));
        for (int b = 0; b < count; ++b) {
            const bpack_gpu_block& blk = batch->blocks[b];
            if (blk.m > 0 && blk.n > 0) refs.push_back({blk.row_ids_host, blk.m, blk.col_ids_host, blk.n});
        }
        const auto t_fetch = clock::now();
        entries->fetch_blocks_direct(refs, image);
        fetching = std::chrono::duration<double>(clock::now() - t_fetch).count();
    } else {
        std::vector<int> missing;
        for (int b = 0; b < count; ++b) {
            const bpack_gpu_block& blk = batch->blocks[b];
            if (blk.m <= 0 || blk.n <= 0) continue;
            if (!entries->lookup(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n,
                                 image + at[static_cast<size_t>(b)], blk.m)) {
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
                                     image + at[static_cast<size_t>(b)], blk.m)) {
                    throw std::logic_error("collective GPU evaluator: a block its round did not fetch");
                }
            }
        }
        for (int b = 0; b < count; ++b) {
            const bpack_gpu_block& blk = batch->blocks[b];
            if (blk.m > 0 && blk.n > 0) entries->release(blk.row_ids_host, blk.m, blk.col_ids_host, blk.n);
        }
    }
    // one upload, then into the blocks on the batch's stream
    const size_t bytes = at.back() * sizeof(DataType);
    char* d_image = static_cast<char*>(batch->allocate(batch->allocator, bytes));
    if (d_image == nullptr) throw std::runtime_error("collective GPU evaluator: no device memory for a batch");
    img.upload(d_image, bytes, stream);
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
// into their slots on the device; each dropped once read.  A request's rows
// are asked for as Omega K(rows, cols) for a symmetric matrix, else as the
// transposed sketch Omega K(cols, rows)^T, and its swap rows (krows) as Omega
// K(krows, cols), the two parts summed.
template<typename DataType, typename Entries>
void collective_sketches(Entries* entries, std::vector<Evaluator::SketchRequest>& requests, cudaStream_t stream) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    double fetching = 0.0;
    const bool transposed_rows = !entries->symmetric();
    std::vector<size_t> at(requests.size() + 1, 0);
    // per request its parts' keys (rows, then krows; d = 0: none) and images
    std::vector<typename Entries::SketchKey> keys;
    std::vector<std::pair<size_t, int>> key_of;  // (request, part)
    for (size_t q = 0; q < requests.size(); ++q) {
        const Evaluator::SketchRequest& r = requests[q];
        at[q + 1] = at[q] + static_cast<size_t>(r.d) * r.cols.size();
        if (!r.rows.empty()) {
            keys.push_back(typename Entries::SketchKey{r.rows.data(), static_cast<int64_t>(r.rows.size()), r.cols.data(),
                                                       static_cast<int64_t>(r.cols.size()), r.d,
                                                       static_cast<int64_t>(r.orow.size()), r.orow.data(), r.odest.data(),
                                                       r.oval.data(), transposed_rows});
            key_of.emplace_back(q, 0);
        }
        if (!r.krows.empty()) {
            keys.push_back(typename Entries::SketchKey{r.krows.data(), static_cast<int64_t>(r.krows.size()),
                                                       r.cols.data(), static_cast<int64_t>(r.cols.size()), r.d,
                                                       static_cast<int64_t>(r.korow.size()), r.korow.data(),
                                                       r.kodest.data(), r.koval.data(), false});
            key_of.emplace_back(q, 1);
        }
    }
    if (at.back() == 0) return;
    std::vector<int> parts(requests.size(), 0);
    for (const auto& kq : key_of) ++parts[kq.first];
    if (entries->alone() && entries->has_sketch_device()) {
        // One rank with a device sketch callback: each request's sketch
        // added straight into its output, zeroed here, an orientation's keys
        // in one call (a request's two parts in one orientation would share
        // their output: the host's then).  Should an orientation decline, the
        // host path below writes every output in full.
        bool clash = false;
        for (int t = 0; t < 2 && !clash; ++t) {
            std::vector<int> seen(requests.size(), 0);
            for (size_t k = 0; k < keys.size(); ++k) {
                if (keys[k].transposed == (t == 1) && seen[key_of[k].first]++ > 0) clash = true;
            }
        }
        if (!clash) {
            for (const Evaluator::SketchRequest& r : requests) {
                if (r.d <= 0 || r.cols.empty()) continue;
                check_cuda(cudaMemset2DAsync(r.out, static_cast<size_t>(r.ldo) * sizeof(DataType), 0,
                                             static_cast<size_t>(r.d) * sizeof(DataType), r.cols.size(), stream),
                           "collective device sketches zero");
            }
            check_cuda(cudaStreamSynchronize(stream), "collective device sketches zero");  // (the callee's streams write next)
            bool done = true;
            for (int t = 0; t < 2 && done; ++t) {
                std::vector<typename Entries::SketchKey> ks;
                std::vector<void*> outs;
                std::vector<int64_t> lds;
                for (size_t k = 0; k < keys.size(); ++k) {
                    if (keys[k].transposed != (t == 1)) continue;
                    const Evaluator::SketchRequest& r = requests[key_of[k].first];
                    ks.push_back(keys[k]);
                    outs.push_back(r.out);
                    lds.push_back(r.ldo);
                }
                const auto t_fetch = clock::now();
                done = entries->fetch_sketches_device(ks, t == 1, outs.data(), lds.data(), stream);
                fetching += std::chrono::duration<double>(clock::now() - t_fetch).count();
            }
            if (done) {
                entries->count_delivery(std::chrono::duration<double>(clock::now() - t0).count() - fetching);
                return;
            }
        }
    }
    CollectiveImage img(at.back() * sizeof(DataType));
    DataType* const image = img.as<DataType>();
    parallel_zero(image, at.back());
    if (entries->alone()) {
        // One rank: each orientation's parts in one callback call, straight
        // into the image when they lie there back to back, one a request
        // (a symmetric matrix's), else added in from a packed copy.
        for (int t = 0; t < 2; ++t) {
            std::vector<typename Entries::SketchKey> ks;
            std::vector<size_t> qs;
            for (size_t k = 0; k < keys.size(); ++k) {
                if (keys[k].transposed != (t == 1)) continue;
                ks.push_back(keys[k]);
                qs.push_back(key_of[k].first);
            }
            if (ks.empty()) continue;
            const size_t base = at[qs[0]];
            size_t total = 0;
            bool in_place = true;
            for (size_t q : qs) {
                in_place = in_place && parts[q] == 1 && at[q] == base + total;
                total += static_cast<size_t>(requests[q].d) * requests[q].cols.size();
            }
            struct PackedTag;  // (kept over the calls: a fresh one faulted its pages in again)
            KeptBuffer<DataType, PackedTag> packed(in_place ? 0 : total);
            parallel_zero(packed.data(), in_place ? 0 : total);
            const auto t_fetch = clock::now();
            entries->fetch_sketches_direct(ks, t == 1, in_place ? image + base : packed.data());
            fetching += std::chrono::duration<double>(clock::now() - t_fetch).count();
            if (in_place) continue;
            // each request's parts added into its slot, a piece of a slot a thread (a request's parts, next
            // to each other in qs, in one piece and in their order: the sums as one thread made them)
            auto len = [&](size_t q) { return static_cast<size_t>(requests[q].d) * requests[q].cols.size(); };
            std::vector<size_t> run_q, run_at{0}, part_off;  // (runs of one request: its parts' offsets in packed)
            for (size_t i = 0, off = 0; i < qs.size(); off += len(qs[i]), ++i) {
                if (run_q.empty() || run_q.back() != qs[i]) {
                    if (!run_q.empty()) run_at.push_back(part_off.size());
                    run_q.push_back(qs[i]);
                }
                part_off.push_back(off);
            }
            run_at.push_back(part_off.size());
            std::vector<size_t> sorted_q(run_q);
            std::sort(sorted_q.begin(), sorted_q.end());
            const bool distinct = std::adjacent_find(sorted_q.begin(), sorted_q.end()) == sorted_q.end();
            constexpr size_t kPiece = size_t{1} << 15;
            std::vector<std::pair<size_t, size_t>> pieces;  // (run, first entry)
            for (size_t r = 0; r < run_q.size(); ++r)
                for (size_t a = 0; a < len(run_q[r]); a += kPiece) pieces.push_back({r, a});
#pragma omp parallel for schedule(dynamic, 4) if (distinct && pieces.size() > 1 && !omp_in_parallel())
            for (size_t p = 0; p < pieces.size(); ++p) {
                const size_t r = pieces[p].first, a = pieces[p].second, b = std::min(len(run_q[r]), a + kPiece);
                DataType* y = image + at[run_q[r]];
                for (size_t j = run_at[r]; j < run_at[r + 1]; ++j) {
                    const DataType* x = packed.data() + part_off[j];
                    for (size_t e = a; e < b; ++e) y[e] += x[e];
                }
            }
        }
    } else {
        std::vector<HostBuffer<DataType>> part(keys.size());  // (each written in full by its lookup)
        std::vector<size_t> missing;
        for (size_t k = 0; k < keys.size(); ++k) {
            const size_t q = key_of[k].first;
            part[k].resize(static_cast<size_t>(requests[q].d) * requests[q].cols.size());
            if (!entries->lookup_sketch(keys[k], part[k].data(), requests[q].d)) {
                entries->record_sketch(keys[k]);
                missing.push_back(k);
            }
        }
        if (!missing.empty()) {
            const auto t_fetch = clock::now();
            entries->fetch_now();
            fetching = std::chrono::duration<double>(clock::now() - t_fetch).count();
            for (size_t k : missing) {
                if (!entries->lookup_sketch(keys[k], part[k].data(), requests[key_of[k].first].d)) {
                    throw std::logic_error("collective GPU evaluator: a sketch its round did not fetch");
                }
            }
        }
        for (const auto& k : keys) entries->release_sketch(k);
        // (a request of one part: its values as they came)
        for (size_t k = 0; k < keys.size(); ++k) {
            DataType* y = image + at[key_of[k].first];
            if (parts[key_of[k].first] == 1) {
                std::copy(part[k].begin(), part[k].end(), y);
            } else {
                for (size_t i = 0; i < part[k].size(); ++i) y[i] += part[k][i];
            }
        }
    }
    DeviceHeap& heap = DeviceHeap::instance();
    char* d_image = heap.alloc(at.back() * sizeof(DataType));
    img.upload(d_image, at.back() * sizeof(DataType), stream);
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
