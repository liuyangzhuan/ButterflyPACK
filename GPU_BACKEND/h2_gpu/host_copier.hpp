#pragma once
// Background copies of device results into host vectors (H2 GPU backend).
//
// A worker thread owns a copy stream and two fixed pinned slots.  Each job is
// a device block, the event after which it is complete, and the destination
// segments inside it; queued jobs stream through the slots as one chunk
// sequence, and the worker appends every chunk to its destination vectors
// (with a small OpenMP team) while the next chunk is in flight.  The job's
// finalize step then runs on the worker (dimensions, checks).
//
// The device heap is not thread-safe, so the owner frees the blocks of
// finished jobs (take_finished) on its own thread.  A job's block is released
// once its last chunk is copied, before the host store.  Errors raised on the
// worker are rethrown by take_finished / wait_released / wait_all.

#ifdef H2_HAVE_GPU

#include "gpu_runtime.hpp"

#include <omp.h>

#include <algorithm>
#include <complex>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace fmm {
namespace gpu {

class HostCopier {
public:
    struct Segment {
        std::vector<double>* doubles = nullptr;  // exactly one destination
        std::vector<int>* ints = nullptr;
        std::vector<std::complex<double>>* complexes = nullptr;
        size_t offset = 0;                       // bytes into the device block
        size_t bytes = 0;
        Segment() = default;
        Segment(std::vector<double>* d, std::vector<int>* i, size_t off, size_t b)
            : doubles(d), ints(i), offset(off), bytes(b) {}
        Segment(std::vector<std::complex<double>>* c, std::vector<int>* i, size_t off, size_t b)
            : ints(i), complexes(c), offset(off), bytes(b) {}
        Segment(std::nullptr_t, std::vector<int>* i, size_t off, size_t b) : ints(i), offset(off), bytes(b) {}
    };
    struct Job {
        char* device = nullptr;
        size_t bytes = 0;
        cudaEvent_t ready = nullptr;             // recorded by the owner; destroyed here
        std::vector<Segment> segments;
        std::function<void()> finalize;
        bool released = false;                   // device block handed back
    };

    // `ring` selects the pinned slot pair (one per concurrent copier).
    // `threads`: at most this many copy threads, and at most half the
    // rank's OpenMP threads, so the copies leave its cores to the host work
    // of the waves (4 ranks per node: 16 OpenMP threads, 8 copy threads; 16
    // spun the waves' OpenMP loops for up to 13 ms each and made the
    // elimination 8-12% slower).
    explicit HostCopier(int device, int ring = 0, int threads = 16) : device_(device), ring_(ring) {
        threads_ = std::max(1, std::min(threads, omp_get_max_threads() / 2));
        worker_ = std::thread([this] { run(); });
    }
    HostCopier(const HostCopier&) = delete;
    HostCopier& operator=(const HostCopier&) = delete;
    ~HostCopier() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        cv_.notify_all();
        worker_.join();
        for (Job& job : queue_) {
            if (job.ready) cudaEventDestroy(job.ready);
        }
    }

    void submit(Job&& job) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            queue_.push_back(std::move(job));
            ++pending_;
            ++unreleased_;
        }
        cv_.notify_all();
    }

    // Device blocks of finished jobs, for the owner to free.
    std::vector<char*> take_finished() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (error_) std::rethrow_exception(error_);
        std::vector<char*> done;
        done.swap(finished_);
        return done;
    }

    // Block until some device block is released (or none is outstanding).
    std::vector<char*> wait_released() {
        std::unique_lock<std::mutex> lock(mutex_);
        idle_cv_.wait(lock, [this] { return !finished_.empty() || unreleased_ == 0 || error_; });
        if (error_) std::rethrow_exception(error_);
        std::vector<char*> done;
        done.swap(finished_);
        return done;
    }

    // Block until every submitted job is finished.
    std::vector<char*> wait_all() {
        std::unique_lock<std::mutex> lock(mutex_);
        idle_cv_.wait(lock, [this] { return pending_ == 0 || error_; });
        if (error_) std::rethrow_exception(error_);
        std::vector<char*> done;
        done.swap(finished_);
        return done;
    }

    double busy_seconds() const { return busy_; }
    double wait_seconds() const { return wait_; }
    double bytes() const { return bytes_; }

private:
    void run() {
        cudaSetDevice(device_);
        cudaStream_t stream = nullptr;
        cudaEvent_t copied[2] = {nullptr, nullptr};
        char* slots[2] = {nullptr, nullptr};
        try {
            check_cuda(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking), "cudaStreamCreate");
            for (auto& ev : copied) check_cuda(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming), "cudaEventCreate");
            for (int s = 0; s < 2; ++s) {
                slots[s] = static_cast<char*>(pinned_pool().ring[2 * ring_ + s].reserve(PinnedPool::kRingSlot));
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(mutex_);
            error_ = std::current_exception();
            idle_cv_.notify_all();
            return;
        }
        while (true) {
            std::deque<Job> jobs;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                cv_.wait(lock, [this] { return stop_ || !queue_.empty(); });
                if (queue_.empty() && stop_) break;
                jobs.swap(queue_);
            }
            try {
                process(jobs, stream, copied, slots);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!error_) error_ = std::current_exception();
            }
            {
                std::lock_guard<std::mutex> lock(mutex_);
                for (Job& job : jobs) {
                    if (!job.released) {  // an error stopped the job early
                        finished_.push_back(job.device);
                        --unreleased_;
                    }
                }
                pending_ -= static_cast<int>(jobs.size());
            }
            idle_cv_.notify_all();
        }
        for (auto& ev : copied) if (ev) cudaEventDestroy(ev);
        if (stream) cudaStreamDestroy(stream);
    }

    void process(std::deque<Job>& jobs, cudaStream_t stream, cudaEvent_t copied[2], char* slots[2]) {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        constexpr size_t slot_bytes = PinnedPool::kRingSlot;
        for (Job& job : jobs) {
            std::sort(job.segments.begin(), job.segments.end(),
                      [](const Segment& a, const Segment& b) { return a.offset < b.offset; });
            #pragma omp parallel for num_threads(threads_) schedule(dynamic, 16)
            for (int64_t i = 0; i < static_cast<int64_t>(job.segments.size()); ++i) {
                Segment& seg = job.segments[static_cast<size_t>(i)];
                if (seg.doubles) {
                    seg.doubles->clear();
                    seg.doubles->reserve(seg.bytes / sizeof(double));
                } else if (seg.complexes) {
                    seg.complexes->clear();
                    seg.complexes->reserve(seg.bytes / sizeof(std::complex<double>));
                } else {
                    seg.ints->clear();
                    seg.ints->reserve(seg.bytes / sizeof(int));
                }
            }
        }
        struct Chunk { size_t job; size_t begin; size_t end; bool first; bool last; };
        std::vector<Chunk> chunks;
        for (size_t j = 0; j < jobs.size(); ++j) {
            const size_t bytes = jobs[j].bytes;
            if (bytes == 0) {
                chunks.push_back({j, 0, 0, true, true});
                continue;
            }
            for (size_t begin = 0; begin < bytes; begin += slot_bytes) {
                const size_t end = std::min(bytes, begin + slot_bytes);
                chunks.push_back({j, begin, end, begin == 0, end == bytes});
            }
        }
        // A chunk's bytes that some segment reads, at their offsets in the
        // slot (ranges less than kGap apart copied together): the parts of a
        // block no segment wants (e.g. the X_NR the eliminator keeps on the
        // device only) stay there.
        constexpr size_t kGap = size_t{128} << 10;
        auto issue = [&](size_t c) {
            const Chunk& ch = chunks[c];
            const Job& job = jobs[ch.job];
            if (ch.first && job.ready) check_cuda(cudaStreamWaitEvent(stream, job.ready, 0), "cudaStreamWaitEvent");
            auto copy = [&](size_t lo, size_t hi) {
                check_cuda(cudaMemcpyAsync(slots[c % 2] + (lo - ch.begin), job.device + lo, hi - lo,
                                           cudaMemcpyDeviceToHost, stream), "host copy");
                bytes_ += static_cast<double>(hi - lo);
            };
            if (ch.end > ch.begin) {
                auto first = std::upper_bound(job.segments.begin(), job.segments.end(), ch.begin,
                                              [](size_t v, const Segment& seg) { return v < seg.offset + seg.bytes; });
                auto last = std::lower_bound(job.segments.begin(), job.segments.end(), ch.end,
                                             [](const Segment& seg, size_t v) { return seg.offset < v; });
                size_t lo = 0, hi = 0;
                bool open = false;
                for (auto it = first; it != last; ++it) {
                    const size_t a = std::max(ch.begin, it->offset), b = std::min(ch.end, it->offset + it->bytes);
                    if (a >= b) continue;
                    if (open && a <= hi + kGap) {
                        hi = std::max(hi, b);
                        continue;
                    }
                    if (open) copy(lo, hi);
                    lo = a;
                    hi = b;
                    open = true;
                }
                if (open) copy(lo, hi);
            }
            check_cuda(cudaEventRecord(copied[c % 2], stream), "cudaEventRecord");
        };
        if (!chunks.empty()) issue(0);
        for (size_t c = 0; c < chunks.size(); ++c) {
            const auto tw = clock::now();
            check_cuda(cudaEventSynchronize(copied[c % 2]), "host copy");
            wait_ += std::chrono::duration<double>(clock::now() - tw).count();
            if (c + 1 < chunks.size()) issue(c + 1);  // the other slot: chunk c - 1 is stored
            const Chunk& ch = chunks[c];
            Job& job = jobs[ch.job];
            if (ch.last) {  // all of the job's device data is in host memory
                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    finished_.push_back(job.device);
                    --unreleased_;
                }
                job.released = true;
                idle_cv_.notify_all();
            }
            const char* slot = slots[c % 2];
            auto first = std::upper_bound(job.segments.begin(), job.segments.end(), ch.begin,
                                          [](size_t v, const Segment& seg) { return v < seg.offset + seg.bytes; });
            auto last = std::lower_bound(job.segments.begin(), job.segments.end(), ch.end,
                                         [](const Segment& seg, size_t v) { return seg.offset < v; });
            const int64_t overlap = last - first;
            #pragma omp parallel for num_threads(threads_) schedule(dynamic, 4)
            for (int64_t i = 0; i < overlap; ++i) {
                const Segment& seg = *(first + i);
                const size_t lo = std::max(ch.begin, seg.offset), hi = std::min(ch.end, seg.offset + seg.bytes);
                const char* src = slot + (lo - ch.begin);
                if (seg.doubles) {
                    const double* p = reinterpret_cast<const double*>(src);
                    seg.doubles->insert(seg.doubles->end(), p, p + (hi - lo) / sizeof(double));
                } else if (seg.complexes) {
                    // (chunk boundaries and segment offsets are multiples of 16 bytes)
                    const std::complex<double>* p = reinterpret_cast<const std::complex<double>*>(src);
                    seg.complexes->insert(seg.complexes->end(), p, p + (hi - lo) / sizeof(std::complex<double>));
                } else {
                    const int* p = reinterpret_cast<const int*>(src);
                    seg.ints->insert(seg.ints->end(), p, p + (hi - lo) / sizeof(int));
                }
            }
            if (ch.last) {
                if (job.ready) {
                    check_cuda(cudaEventDestroy(job.ready), "cudaEventDestroy");
                    job.ready = nullptr;
                }
                if (job.finalize) job.finalize();
            }
        }
        busy_ += std::chrono::duration<double>(clock::now() - t0).count();
    }

    int device_;
    int ring_ = 0;
    int threads_ = 16;
    std::thread worker_;
    std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable idle_cv_;
    std::deque<Job> queue_;
    std::vector<char*> finished_;
    int pending_ = 0;     // jobs not yet stored and finalized
    int unreleased_ = 0;  // jobs whose device block is not yet handed back
    bool stop_ = false;
    std::exception_ptr error_;
    double busy_ = 0.0, wait_ = 0.0, bytes_ = 0.0;
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
