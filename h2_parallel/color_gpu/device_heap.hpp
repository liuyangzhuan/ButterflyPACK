#pragma once
// Device memory of the H2 Color GPU backend: one large allocation per rank,
// carved by a host-side allocator.  Near-field and Schur blocks change shape
// as boxes are eliminated, so a level allocates and frees tens of thousands
// of blocks per wave; host bookkeeping keeps that off the CUDA allocator.
// All device work runs on one stream, so memory freed on the host may be
// reused by any later launch without a synchronization.
//
// Blocks are placed by lifetime: resident ones (the level's Schur and
// near-field blocks, fill sources, point tables) take the best-fitting free
// range (the lowest among equals), per-wave buffers the highest range that
// fits.  Holes left by freed resident blocks are then refilled by later
// resident blocks, which grow from the bottom, and the wave buffers find
// contiguous space at the top instead of the gaps between long-lived blocks.
// Free ranges are indexed by address (coalescing) and by size (best fit).

#ifdef H2_HAVE_GPU

#include "gpu_runtime.hpp"

#include <cstdlib>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

class DeviceHeap {
public:
    static DeviceHeap& instance() {
        static DeviceHeap heap;
        return heap;
    }
    // A second, small arena for the buffers MPI reads and writes (device
    // exchange): CUDA-aware MPI maps or registers the whole allocation that
    // holds a buffer when it first meets it, at a cost that grows with the
    // allocation.  Reserved before the main arena when possible.
    static DeviceHeap& exchange_arena() {
        static DeviceHeap heap;
        return heap;
    }

    DeviceHeap(const DeviceHeap&) = delete;
    DeviceHeap& operator=(const DeviceHeap&) = delete;

    // Reserve the arena on first use: H2_GPU_HEAP_GB GiB if set, else a
    // fraction (H2_GPU_HEAP_FRACTION, default 0.85) of the device memory free
    // at that point, or, when several ranks share the device
    // (set_device_share), an equal part of that fraction of the memory they
    // found free together.  MAGMA's workspaces and the pinned staging
    // buffers stay outside the arena.
    void ensure_initialized() {
        if (base_ != nullptr) return;
        Context::instance().activate();
        size_t free_bytes = 0, total_bytes = 0;
        check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
        if (const char* env = std::getenv("H2_GPU_HEAP_GB")) {
            const double bytes = std::atof(env) * static_cast<double>(size_t{1} << 30);
            if (!(bytes > 0.0) || bytes > static_cast<double>(free_bytes)) {
                throw std::invalid_argument("H2_GPU_HEAP_GB must be positive and at most the free device memory (" +
                                            std::to_string(static_cast<double>(free_bytes) / (size_t{1} << 30)) + " GiB)");
            }
            capacity_ = static_cast<size_t>(bytes) / kAlign * kAlign;
        } else {
            double fraction = 0.85;
            if (const char* frac = std::getenv("H2_GPU_HEAP_FRACTION")) fraction = std::atof(frac);
            if (!(fraction > 0.0 && fraction < 1.0)) {
                throw std::invalid_argument("H2_GPU_HEAP_FRACTION must lie in (0, 1)");
            }
            if (share_ranks_ > 1) {
                free_bytes = std::min(free_bytes, share_free_);
                capacity_ = static_cast<size_t>(fraction * static_cast<double>(share_free_) / share_ranks_) / kAlign * kAlign;
                capacity_ = std::min(capacity_, static_cast<size_t>(fraction * static_cast<double>(free_bytes)) / kAlign * kAlign);
            } else {
                capacity_ = static_cast<size_t>(fraction * static_cast<double>(free_bytes)) / kAlign * kAlign;
            }
        }
        void* ptr = nullptr;
        check_cuda(cudaMalloc(&ptr, capacity_), "cudaMalloc (H2 GPU heap)");
        base_ = static_cast<char*>(ptr);
        reset();
    }

    // A fixed-size arena (no-op once reserved); false if the device has no
    // room for it.
    bool initialize_fixed(size_t bytes) {
        if (base_ != nullptr) return true;
        Context::instance().activate();
        void* ptr = nullptr;
        const size_t capacity = (bytes + kAlign - 1) / kAlign * kAlign;
        if (capacity == 0 || cudaMalloc(&ptr, capacity) != cudaSuccess) {
            cudaGetLastError();  // clear the failed allocation
            return false;
        }
        base_ = static_cast<char*>(ptr);
        capacity_ = capacity;
        reset();
        return true;
    }
    bool initialized() const { return base_ != nullptr; }
    // `ranks` ranks (this one included) use this device and found
    // free_bytes free once all their contexts existed: the arena reserved
    // later takes an equal part (before ensure_initialized; no effect after)
    void set_device_share(int ranks, size_t free_bytes) {
        share_ranks_ = std::max(ranks, 1);
        share_free_ = free_bytes;
    }
    bool owns(const void* ptr) const {
        const char* p = static_cast<const char*>(ptr);
        return base_ != nullptr && p >= base_ && p < base_ + capacity_;
    }
    // A transient block if it fits now (no reclaim), else nullptr.
    char* try_alloc(size_t bytes) {
        if (base_ == nullptr) return nullptr;
        const size_t size = std::max<size_t>(kAlign, (bytes + kAlign - 1) / kAlign * kAlign);
        size_t offset = 0;
        if (!take(size, false, offset)) return nullptr;
        allocated_.emplace(offset, size);
        used_ += size;
        peak_ = std::max(peak_, used_);
        return base_ + offset;
    }

    // A per-wave (transient) block, from the top of the arena.
    template<typename T = char>
    T* alloc(size_t bytes) {
        return reinterpret_cast<T*>(alloc_bytes(bytes, false));
    }
    // A block that lives for a level or longer, from the bottom.
    template<typename T = char>
    T* alloc_resident(size_t bytes) {
        return reinterpret_cast<T*>(alloc_bytes(bytes, true));
    }

    // Called when an allocation does not fit: frees what it can (blocks
    // whose background copies are done) and returns whether it freed any.
    void set_reclaimer(std::function<bool()> reclaimer) { reclaimer_ = std::move(reclaimer); }

    char* alloc_bytes(size_t bytes, bool resident) {
        ensure_initialized();
        const size_t size = std::max<size_t>(kAlign, (bytes + kAlign - 1) / kAlign * kAlign);
        size_t offset = 0;
        while (!take(size, resident, offset)) {
            if (!(reclaimer_ && reclaimer_())) {
                std::ostringstream oss;
                oss << "H2 GPU heap exhausted: request " << size << " bytes, in use " << used_ << " of "
                    << capacity_ << " bytes (largest free block " << largest_free()
                    << "); use more ranks or a smaller problem per GPU";
                throw std::runtime_error(oss.str());
            }
            ++reclaims_;
        }
        allocated_.emplace(offset, size);
        used_ += size;
        peak_ = std::max(peak_, used_);
        return base_ + offset;
    }

    void free(const void* ptr) {
        if (ptr == nullptr) return;
        const size_t offset = static_cast<size_t>(static_cast<const char*>(ptr) - base_);
        auto it = allocated_.find(offset);
        if (it == allocated_.end()) throw std::runtime_error("H2 GPU heap: free of an unknown pointer");
        size_t start = offset, length = it->second;
        allocated_.erase(it);
        used_ -= length;
        auto next = free_.lower_bound(start);
        if (next != free_.begin()) {
            auto prev = std::prev(next);
            if (prev->first + prev->second == start) {
                start = prev->first;
                length += prev->second;
                remove_free(prev);
            }
        }
        if (next != free_.end() && start + length == next->first) {
            length += next->second;
            remove_free(next);
        }
        add_free(start, length);
    }

    // Free every block (start of a factorization).
    void reset() {
        free_.clear();
        by_size_.clear();
        allocated_.clear();
        used_ = 0;
        if (capacity_ > 0) add_free(0, capacity_);
    }

    size_t used() const { return used_; }
    size_t largest_free() const { return by_size_.empty() ? 0 : by_size_.rbegin()->first; }
    size_t peak() const { return peak_; }
    size_t capacity() const { return capacity_; }
    void reset_peak() { peak_ = used_; }
    int64_t reclaims() const { return reclaims_; }

private:
    static constexpr size_t kAlign = 256;

    DeviceHeap() = default;
    ~DeviceHeap() = default;  // released at process exit (see Context)

    void add_free(size_t offset, size_t length) {
        free_.emplace(offset, length);
        by_size_.emplace(length, offset);
    }
    void remove_free(std::map<size_t, size_t>::iterator it) {
        by_size_.erase({it->second, it->first});
        free_.erase(it);
    }

    // Best fit (lowest address among equal sizes), cut from the start of the
    // range (resident), or the highest range that fits, cut from its end
    // (per-wave).
    bool take(size_t size, bool resident, size_t& offset) {
        if (resident) {
            auto best = by_size_.lower_bound({size, 0});
            if (best == by_size_.end()) return false;
            const size_t start = best->second, length = best->first;
            remove_free(free_.find(start));
            offset = start;
            if (length > size) add_free(start + size, length - size);
            return true;
        }
        for (auto it = free_.rbegin(); it != free_.rend(); ++it) {
            if (it->second < size) continue;
            const size_t start = it->first, length = it->second;
            remove_free(free_.find(start));
            offset = start + length - size;
            if (length > size) add_free(start, length - size);
            return true;
        }
        return false;
    }

    int share_ranks_ = 1;
    size_t share_free_ = 0;
    char* base_ = nullptr;
    size_t capacity_ = 0;
    size_t used_ = 0;
    size_t peak_ = 0;
    std::map<size_t, size_t> free_;                  // offset -> length, address order
    std::set<std::pair<size_t, size_t>> by_size_;    // (length, offset) of the same ranges
    std::unordered_map<size_t, size_t> allocated_;   // offset -> length
    std::function<bool()> reclaimer_;
    int64_t reclaims_ = 0;
};

// ---------------------------------------------------------------------------
// Several H2 operators in one process.  With H2_GPU_KEEP_OPERATORS=1 (the
// default), the device data a factorization or compression leaves for the
// solve and the matvec (device_solve_store(), device_matvec_store()) is kept
// per operator (its tree), so a process that alternates between operators (a
// Gaussian process with its covariance matrix and the derivative matrices of
// its gradient) finds their device data again instead of falling back to the
// host.  With a single operator nothing changes: a new factorization or
// compression of an operator replaces its own data.  0: one set of device
// data, replaced by every factorization or compression.
//
// The operators of a process are assumed to share its ranks, which call
// their H2 operations in the same order: the use order below is then the
// same on every rank, and so are the evictions.  The least recently used
// operators' data is released collectively at the start of a build while the
// other operators hold more than half of some rank's heap (an eviction on
// one rank only would split the ranks between the device and host paths of
// that operator); an allocation that still does not fit fails as with 0.
inline bool keep_operators() {
    static const bool keep = [] {
        const char* v = std::getenv("H2_GPU_KEEP_OPERATORS");
        return v == nullptr || std::atoi(v) != 0;
    }();
    return keep;
}

struct OperatorContext {
    const void* active = nullptr;                         // the operator of the current H2 call
    uint64_t clock = 0;
    std::unordered_map<const void*, uint64_t> last_use;   // keep_operators(): operator -> last use
    // per kind of store (solve, matvec): release an operator's data; bytes it holds
    std::vector<std::function<void(const void*)>> releasers;
    std::vector<std::function<size_t(const void*)>> holders;
};

inline OperatorContext& operator_context() {
    static OperatorContext context;
    return context;
}

// An H2 call on operator `tree` starts: the stores now serve it.
inline void activate_operator(const void* tree) {
    OperatorContext& c = operator_context();
    c.active = tree;
    if (keep_operators()) c.last_use[tree] = ++c.clock;
}

// Device bytes the stores hold for operator `tree`.
inline size_t operator_bytes(const void* tree) {
    size_t bytes = 0;
    for (const auto& held : operator_context().holders) bytes += held(tree);
    return bytes;
}

// Whether operators other than the active one hold device data.
inline bool other_operators_hold_data() {
    if (!keep_operators()) return false;
    const OperatorContext& c = operator_context();
    for (const auto& entry : c.last_use) {
        if (entry.first != c.active && operator_bytes(entry.first) > 0) return true;
    }
    return false;
}

// The device data of operator `tree` is released (its destruction, on every rank).
inline void release_operator(const void* tree) {
    OperatorContext& c = operator_context();
    if (keep_operators()) {
        for (const auto& release : c.releasers) release(tree);
        c.last_use.erase(tree);
    }
    if (c.active == tree) c.active = nullptr;
}

// Start of a factorization or compression of operator `tree` (collective
// over comm): activates it and makes room from the other operators' data.
inline void begin_operator_build(const void* tree, MPI_Comm comm) {
    activate_operator(tree);
    if (!keep_operators()) return;
    OperatorContext& c = operator_context();
    const DeviceHeap& heap = DeviceHeap::instance();
    for (;;) {
        size_t others = 0;
        const void* oldest = nullptr;
        uint64_t oldest_use = 0;
        for (const auto& [t, use] : c.last_use) {
            if (t == tree) continue;
            others += operator_bytes(t);
            if (oldest == nullptr || use < oldest_use) {
                oldest = t;
                oldest_use = use;
            }
        }
        int over = others > heap.capacity() / 2 ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &over, 1, MPI_INT, MPI_MAX, comm);
        if (!over || oldest == nullptr) return;
        for (const auto& release : c.releasers) release(oldest);
        c.last_use.erase(oldest);
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
