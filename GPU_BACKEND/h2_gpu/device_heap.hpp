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
// Free ranges are indexed by address (coalescing), by size (best fit) and by
// size class (the highest range that fits, without a walk: see take).

#ifdef H2_HAVE_GPU

#include "gpu_runtime.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <array>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <dlfcn.h>
#include <sstream>
#include <iomanip>
#include <string>

#include "../gpu_common/bpack_env.hpp"
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// What a resident block is for, so a level's resident peak can be broken
// down (BPACK_TRACE=phase).  Resident blocks are the ones that live for a
// level or longer, which is what a level's memory ceiling is made of; the
// transient ones come and go within a wave.  A block with no tag counts as
// Other.
enum class ResidentKind : uint8_t {
    Other = 0,
    FillXrr,       // keep_fill_source: X_RR_full, per box (never compressed)
    FillRaw,       // ... its temp2^T, while it is still stored raw
    FillPacked,    // ... its temp2^T, compressed a slot at a time
    LevelSchur,    // this level's own Schur and edge blocks, accumulated into by its waves
    Transition,    // the parents' Schur and edge blocks, built by the transition
    Near,          // step 5's dense neighbor blocks
    Points,        // the level's coordinates and point ids
    Lists,         // skeleton lists, training lists, index tables
    KeptFactors,   // the solve factors a level keeps (H2_compress_factors)
    Restore,       // a transition chunk's restored fill sources
    CaHalo,        // a CA level's received halo pieces
    Reduction,     // the process reduction's pack and receive buffers
    Matvec,        // the matvec's device blocks
    // Three things used to share Transition and so printed as "parent
    // blocks", which made the leaf's 9.47 GB look like parent blocks at a
    // level that adopts none:
    Uploaded,      // begin(): the host's blocks uploaded for a level
    Adopted,       // begin(): the parent blocks of the level below, restored
    OwnerNear,     // the owner pass's new near-field blocks (its own edges_)
    OwnerSrc,      // X_NR^T slots, the owner pass's other operand, built and freed per wave
    // Transient kinds (the top of the arena), so the transient half of the
    // peak can be attributed the way the resident half already is:
    WaveWork,      // a box's tmp1, tmp2, X_SR, A_NS, X_NS
    WaveOutput,    // its X_NR and temp2, until their stores exist
    WaveSource,    // the source block the host copier takes (X_RR_full, X_NR orig)
    WaveResult,    // the result block (LU, pivots, temp1, X_RS, T)
    Codec,         // compressor staging and decompression destinations
    SrcWindow,     // the owner pass's window of source slots
    Sketch,        // the wave's sketch buffers, split four ways because
    SketchFill,    // ... the decompressed fill slots a group reads
    SketchVec,    // ... the random vectors it applies the matrix to
    SketchLists,   // ... per box, the source lists and index tables
    SketchPW,      // ... per box, [P; W]
    TrRestore,     // the transition chunk's fill sources, restored from the host
    TrBuffers,     // ... and its G and P buffers
    Exchange,      // halo and generator buffers that did not fit the exchange arena
    Gather,        // the near store's chunk homes, gathered contiguous to be read
    Adaptive,      // the adaptive ID's frames and status
    GenUpload,     // remote generators, uploaded: temp2, X_NR, X_RS, temp2^T, X_RR_full
    GenWork,       // ... and the buffers of the products they drive
    Count
};
inline const char* resident_kind_name(ResidentKind k) {
    switch (k) {
        case ResidentKind::FillXrr: return "fill X_RR_full";
        case ResidentKind::FillRaw: return "fill temp2^T raw";
        case ResidentKind::FillPacked: return "fill temp2^T packed";
        case ResidentKind::LevelSchur: return "level Schur/edges";
        case ResidentKind::Transition: return "parent blocks";
        case ResidentKind::Near: return "near blocks";
        case ResidentKind::Points: return "points";
        case ResidentKind::Lists: return "lists";
        case ResidentKind::KeptFactors: return "kept factors";
        case ResidentKind::Restore: return "restores";
        case ResidentKind::CaHalo: return "CA halo";
        case ResidentKind::OwnerSrc: return "owner source slots";
        case ResidentKind::WaveWork: return "wave work";
        case ResidentKind::WaveOutput: return "wave outputs";
        case ResidentKind::WaveSource: return "wave sources";
        case ResidentKind::WaveResult: return "wave result";
        case ResidentKind::Codec: return "codec staging";
        case ResidentKind::SrcWindow: return "source window";
        case ResidentKind::Sketch: return "sketch";
        case ResidentKind::SketchFill: return "sketch fill slots";
        case ResidentKind::SketchVec: return "sketch vectors";
        case ResidentKind::SketchLists: return "sketch lists";
        case ResidentKind::SketchPW: return "sketch P;W";
        case ResidentKind::TrRestore: return "transition restores";
        case ResidentKind::TrBuffers: return "transition G;P";
        case ResidentKind::Exchange: return "exchange overflow";
        case ResidentKind::Gather: return "store gathers";
        case ResidentKind::Adaptive: return "adaptive ID";
        case ResidentKind::GenUpload: return "generator uploads";
        case ResidentKind::GenWork: return "generator products";
        case ResidentKind::Reduction: return "reduction";
        case ResidentKind::Matvec: return "matvec";
        case ResidentKind::Uploaded: return "uploaded blocks";
        case ResidentKind::Adopted: return "adopted blocks";
        case ResidentKind::OwnerNear: return "owner near blocks";
        default: return "other";
    }
}

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

    // Reserve the arena on first use: a fraction (BPACK_GPU_HEAP_FRACTION,
    // default 0.85) of the device memory free at that point, or, when
    // several ranks share the device (share_device), an equal part of that
    // fraction of the memory they found free together.  MAGMA's workspaces
    // and the pinned staging buffers stay outside the arena.
    void ensure_initialized() {
        if (base_ != nullptr) return;
        Context::instance().activate();
        size_t free_bytes = 0, total_bytes = 0;
        check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
        double fraction = 0.85;
        if (const char* env = std::getenv("BPACK_GPU_HEAP_FRACTION")) fraction = std::atof(env);
        if (!(fraction > 0.0 && fraction < 1.0)) {
            throw std::invalid_argument("BPACK_GPU_HEAP_FRACTION must lie in (0, 1)");
        }
        if (share_ranks_ > 1) {
            free_bytes = std::min(free_bytes, share_free_);
            capacity_ = static_cast<size_t>(fraction * static_cast<double>(share_free_) / share_ranks_) / kAlign * kAlign;
            capacity_ = std::min(capacity_, static_cast<size_t>(fraction * static_cast<double>(free_bytes)) / kAlign * kAlign);
        } else {
            capacity_ = static_cast<size_t>(fraction * static_cast<double>(free_bytes)) / kAlign * kAlign;
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
    // Whether this pointer is the start of a block the heap currently holds.
    // owns() only says the address is inside the arena, which does not tell a
    // foreign pointer from one that has already been freed.
    bool is_allocated(const void* ptr) const {
        if (!owns(ptr)) return false;
        const size_t offset = static_cast<size_t>(static_cast<const char*>(ptr) - base_);
        return allocated_.find(offset) != allocated_.end();
    }
    bool owns(const void* ptr) const {
        const char* p = static_cast<const char*>(ptr);
        return base_ != nullptr && p >= base_ && p < base_ + capacity_;
    }
    // The length of a block this heap handed out (0: not one of its blocks).
    size_t bytes_of(const void* ptr) const {
        if (!owns(ptr)) return 0;
        const auto it = allocated_.find(static_cast<size_t>(static_cast<const char*>(ptr) - base_));
        return it == allocated_.end() ? 0 : it->second.first;
    }
    // A transient block if it fits now (no reclaim), else nullptr.
    char* try_alloc(size_t bytes, ResidentKind kind = ResidentKind::Other) {
        if (base_ == nullptr) return nullptr;
        const size_t size = std::max<size_t>(kAlign, (bytes + kAlign - 1) / kAlign * kAlign);
        size_t offset = 0;
        if (!take(size, false, offset)) return nullptr;
        allocated_.emplace(offset, std::make_pair(size, false));
        used_ += size;
        transient_used_ += size;
        note_transient(offset, size, kind);
        note_peak();
        return base_ + offset;
    }

    // A per-wave (transient) block, from the top of the arena.
    template<typename T = char>
    T* alloc(size_t bytes, ResidentKind kind = ResidentKind::Other) {
        return reinterpret_cast<T*>(alloc_bytes(bytes, false, kind));
    }
    // A resident block if it fits now (no reclaim), else nullptr: for a
    // caller that has a correct fallback, such as storing a fill source raw
    // when its compressed copy will not fit.
    char* try_alloc_resident(size_t bytes, ResidentKind kind = ResidentKind::Other) {
        if (base_ == nullptr) return nullptr;
        const size_t size = std::max<size_t>(kAlign, (bytes + kAlign - 1) / kAlign * kAlign);
        size_t offset = 0;
        if (!take(size, true, offset)) return nullptr;
        allocated_.emplace(offset, std::make_pair(size, true));
        used_ += size;
        resident_used_ += size;
        resident_kind_[offset] = kind;
        by_kind_[static_cast<size_t>(kind)] += size;
        ++alloc_n_[static_cast<size_t>(kind)];
        alloc_b_[static_cast<size_t>(kind)] += size;
        note_peak();
        return base_ + offset;
    }

    // A block that lives for a level or longer, from the bottom.  `kind` is
    // only for the resident breakdown.
    template<typename T = char>
    T* alloc_resident(size_t bytes, ResidentKind kind = ResidentKind::Other) {
        return reinterpret_cast<T*>(alloc_bytes(bytes, true, kind));
    }

    // Called when an allocation does not fit: frees what it can (blocks
    // whose background copies are done) and returns whether it freed any.
    void set_reclaimer(std::function<bool()> reclaimer) { reclaimer_ = std::move(reclaimer); }

    char* alloc_bytes(size_t bytes, bool resident, ResidentKind kind = ResidentKind::Other) {
        ensure_initialized();
        const size_t size = std::max<size_t>(kAlign, (bytes + kAlign - 1) / kAlign * kAlign);
        size_t offset = 0;
        while (!take(size, resident, offset)) {
            if (!(reclaimer_ && reclaimer_())) {
                std::ostringstream oss;
                oss << "H2 GPU heap exhausted: request " << size << " bytes (" << resident_kind_name(kind)
                    << (resident ? ", resident" : ", transient") << "), in use " << used_ << " of "
                    << capacity_ << " bytes (largest free block " << largest_free()
                    << "); use more ranks or a smaller problem per GPU";
                // what holds it, the largest kinds first (a failure is the time to know)
                std::vector<std::pair<size_t, size_t>> kinds;
                for (size_t k = 0; k < by_kind_.size(); ++k) if (by_kind_[k] + by_tkind_[k] > (size_t{64} << 20)) kinds.push_back({by_kind_[k] + by_tkind_[k], k});
                std::sort(kinds.rbegin(), kinds.rend());
                oss << "; held: ";
                for (size_t i = 0; i < kinds.size() && i < 8; ++i) {
                    const size_t k = kinds[i].second;
                    oss << (i ? ", " : "") << resident_kind_name(static_cast<ResidentKind>(k)) << " " << std::fixed
                        << std::setprecision(2) << static_cast<double>(by_kind_[k]) / 1e9 << "+"
                        << static_cast<double>(by_tkind_[k]) / 1e9 << " GB";
                }
                oss << " (resident+transient); " << free_.size() << " free ranges";
                throw std::runtime_error(oss.str());
            }
            ++reclaims_;
        }
        allocated_.emplace(offset, std::make_pair(size, resident));
        if (heap_trace_) {
            alloc_site_[offset] = __builtin_return_address(0);
            free_site_.erase(offset);
        }
        used_ += size;
        if (resident) {
            resident_used_ += size;
            resident_kind_[offset] = kind;
            by_kind_[static_cast<size_t>(kind)] += size;
            ++alloc_n_[static_cast<size_t>(kind)];
            alloc_b_[static_cast<size_t>(kind)] += size;
        } else {
            transient_used_ += size;
            note_transient(offset, size, kind);
        }
        note_peak();
        return base_ + offset;
    }

    void free(const void* ptr) {
        if (ptr == nullptr) return;
        const size_t offset = static_cast<size_t>(static_cast<const char*>(ptr) - base_);
        auto it = allocated_.find(offset);
        if (it == allocated_.end()) {
            if (heap_trace_) {
                std::ostringstream t;
                t << "H2 GPU heap: double free at offset " << offset << " (" << ptr << ")";
                const auto a = alloc_site_.find(offset);
                const auto f = free_site_.find(offset);
                auto name_of = [](const void* addr) {
                    Dl_info info{};
                    if (addr != nullptr && dladdr(addr, &info) != 0 && info.dli_sname != nullptr) {
                        return std::string(info.dli_sname) + "+" +
                               std::to_string(static_cast<const char*>(addr) -
                                              static_cast<const char*>(info.dli_saddr));
                    }
                    return std::string("?");
                };
                const auto k = freed_kind_.find(offset);
                if (k != freed_kind_.end()) {
                    t << ", " << k->second.second << " bytes of resident kind "
                      << static_cast<int>(k->second.first);
                }
                t << "\n  allocated by  " << name_of(a != alloc_site_.end() ? a->second : nullptr)
                  << "\n  last freed by " << name_of(f != free_site_.end() ? f->second : nullptr)
                  << "\n  freed again by " << name_of(__builtin_return_address(0));
                throw std::runtime_error(t.str());
            }
            throw std::runtime_error("H2 GPU heap: free of an unknown pointer");
        }
        if (heap_trace_) {
            free_site_[offset] = __builtin_return_address(0);
            ResidentKind fk = ResidentKind::Other;
            if (auto kit = resident_kind_.find(offset); kit != resident_kind_.end()) fk = kit->second;
            freed_kind_[offset] = {fk, it->second.first};
        }
        size_t start = offset, length = it->second.first;
        const bool was_resident = it->second.second;
        allocated_.erase(it);
        used_ -= length;
        if (was_resident) {
            resident_used_ -= length;
            if (auto kit = resident_kind_.find(offset); kit != resident_kind_.end()) {
                by_kind_[static_cast<size_t>(kit->second)] -= length;
                ++free_n_[static_cast<size_t>(kit->second)];
                free_b_[static_cast<size_t>(kit->second)] += length;
                resident_kind_.erase(kit);
            }
        } else {
            transient_used_ -= length;
            if (auto tit = transient_kind_.find(offset); tit != transient_kind_.end()) {
                by_tkind_[static_cast<size_t>(tit->second)] -= length;
                transient_kind_.erase(tit);
            }
        }
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
        for (auto& b : top_bins_) b.clear();
        allocated_.clear();
        used_ = 0;
        resident_used_ = 0;
        transient_used_ = 0;
        resident_kind_.clear();
        by_kind_.fill(0);
        by_kind_at_peak_.fill(0);
        transient_kind_.clear();
        by_tkind_.fill(0);
        by_tkind_at_peak_.fill(0);
        if (capacity_ > 0) add_free(0, capacity_);
    }

    size_t used() const { return used_; }
    size_t largest_free() const { return by_size_.empty() ? 0 : by_size_.rbegin()->first; }
    size_t peak() const { return peak_; }
    // the allocator's work: free ranges now, transient takes and the ranges they scanned from the top
    size_t free_ranges() const { return free_.size(); }
    size_t take_calls() const { return take_calls_; }
    size_t take_scans() const { return take_scans_; }
    size_t resident_used() const { return resident_used_; }
    size_t transient_used() const { return transient_used_; }
    // The split of peak(): what each class held when it was reached.
    size_t resident_at_peak() const { return resident_at_peak_; }
    size_t transient_at_peak() const { return transient_at_peak_; }
    size_t capacity() const { return capacity_; }
    void reset_peak() {
        peak_ = used_;
        resident_at_peak_ = resident_used_;
        transient_at_peak_ = transient_used_;
        by_kind_at_peak_ = by_kind_;
        by_tkind_at_peak_ = by_tkind_;
    }
    // What each kind of resident block held when the peak was reached.
    size_t resident_at_peak(ResidentKind k) const { return by_kind_at_peak_[static_cast<size_t>(k)]; }
    // The same for the transient half, which had no such accounting: a wave's
    // buffers were reported as what it plans to allocate, and reading that as
    // the composition of the peak cost a design.
    size_t transient_at_peak(ResidentKind k) const { return by_tkind_at_peak_[static_cast<size_t>(k)]; }
    size_t resident_now(ResidentKind k) const { return by_kind_[static_cast<size_t>(k)]; }
    // Cumulative, so a kind's current bytes can be checked against what was
    // ever allocated to it.  A kind reporting bytes at the peak that it was
    // never allocated means the accounting is wrong, not the memory.
    // Compaction counts as a free and an allocation, so alloc - free = now.
    size_t kind_alloc_bytes(ResidentKind k) const { return alloc_b_[static_cast<size_t>(k)]; }
    size_t kind_free_bytes(ResidentKind k) const { return free_b_[static_cast<size_t>(k)]; }
    size_t kind_alloc_count(ResidentKind k) const { return alloc_n_[static_cast<size_t>(k)]; }
    size_t kind_free_count(ResidentKind k) const { return free_n_[static_cast<size_t>(k)]; }
    int64_t reclaims() const { return reclaims_; }

    // Compaction of long-lived blocks.  Blocks that outlive a level (the
    // solve factors a level keeps, one block per wave, and the parents'
    // blocks of a device transition) are placed among the level's own
    // blocks, and once those are freed they peg the free memory into pieces
    // (640^3 on 512 ranks: 10 GB free in pieces under 1.3 GB).  The listed
    // blocks slide towards the start of the arena, in address order: each to
    // the lowest free range that holds it, or, when the range right below it
    // is free, down to its start (an overlapping move).  `move(dst, src,
    // bytes)` copies a block (dst < src; the ranges overlap when dst + bytes
    // > src).  Unlisted blocks stay where they are, and the free space
    // between them merges.  `blocks` holds the new addresses on return; the
    // moves made, by old address, are returned for pointers into the blocks.
    struct Move {
        char* from;
        char* to;
        size_t bytes;
    };
    std::vector<Move> compact(std::vector<char*>& blocks,
                              const std::function<void(char*, const char*, size_t)>& move) {
        std::vector<Move> moves;
        std::vector<size_t> order(blocks.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return blocks[a] < blocks[b]; });
        for (size_t i : order) {
            if (blocks[i] == nullptr) continue;
            const size_t offset = static_cast<size_t>(blocks[i] - base_);
            auto alloc = allocated_.find(offset);
            if (alloc == allocated_.end()) throw std::runtime_error("H2 GPU heap: compaction of an unknown block");
            const size_t size = alloc->second.first;
            const bool was_resident = alloc->second.second;
            // the move frees and re-registers the block, so carry its kind
            ResidentKind kind = ResidentKind::Other;
            if (auto kit = resident_kind_.find(offset); kit != resident_kind_.end()) kind = kit->second;
            size_t dest = offset;
            for (auto it = free_.begin(); it != free_.end() && it->first < offset; ++it) {
                if (it->second >= size) {
                    dest = it->first;
                    break;
                }
            }
            if (dest == offset) {
                auto next = free_.lower_bound(offset);  // (the block's own range is not free)
                if (next != free_.begin()) {
                    auto below = std::prev(next);
                    if (below->first + below->second == offset) dest = below->first;
                }
            }
            if (dest == offset) continue;
            move(base_ + dest, base_ + offset, size);
            free(blocks[i]);  // (merges with the destination range when adjacent)
            auto range = free_.find(dest);
            if (range == free_.end() || range->second < size) {
                throw std::runtime_error("H2 GPU heap: compaction lost its destination range");
            }
            const size_t length = range->second;
            remove_free(range);
            if (length > size) add_free(dest + size, length - size);
            allocated_.emplace(dest, std::make_pair(size, was_resident));
            used_ += size;
            if (was_resident) {
                resident_used_ += size;
                resident_kind_[dest] = kind;
                by_kind_[static_cast<size_t>(kind)] += size;
                ++alloc_n_[static_cast<size_t>(kind)];
                alloc_b_[static_cast<size_t>(kind)] += size;
            } else {
                transient_used_ += size;
                note_transient(dest, size, kind);
            }
            moves.push_back({blocks[i], base_ + dest, size});
            blocks[i] = base_ + dest;
        }
        std::sort(moves.begin(), moves.end(), [](const Move& a, const Move& b) { return a.from < b.from; });
        return moves;
    }
    // Pointer `p` into a moved block, after the move (unchanged otherwise);
    // `moves` as compact returns them.
    template<typename T>
    static T* relocated(T* p, const std::vector<Move>& moves) {
        const char* c = reinterpret_cast<const char*>(p);
        auto it = std::upper_bound(moves.begin(), moves.end(), c,
                                   [](const char* x, const Move& m) { return x < m.from; });
        if (it == moves.begin()) return p;
        --it;
        if (c < it->from || c >= it->from + it->bytes) return p;
        return reinterpret_cast<T*>(const_cast<char*>(it->to + (c - it->from)));
    }
    // Free bytes in all.
    size_t available() const { return capacity_ - used_; }

private:
    void note_peak() {
        if (used_ <= peak_) return;
        peak_ = used_;
        resident_at_peak_ = resident_used_;
        transient_at_peak_ = transient_used_;
        by_kind_at_peak_ = by_kind_;
        by_tkind_at_peak_ = by_tkind_;
    }
    void note_transient(size_t offset, size_t size, ResidentKind kind) {
        transient_kind_[offset] = kind;
        by_tkind_[static_cast<size_t>(kind)] += size;
    }

    static constexpr size_t kAlign = 256;

    DeviceHeap() = default;
    ~DeviceHeap() = default;  // released at process exit (see Context)

    void add_free(size_t offset, size_t length) {
        free_.emplace(offset, length);
        by_size_.emplace(length, offset);
        top_bins_[bin_of(length)].insert(offset);
    }
    void remove_free(std::map<size_t, size_t>::iterator it) {
        by_size_.erase({it->second, it->first});
        top_bins_[bin_of(it->second)].erase(it->first);
        free_.erase(it);
    }
    // The size class of a length: floor(log2), so a class spans a factor of two.
    static int bin_of(size_t length) { return 63 - __builtin_clzll(length | 1); }

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
        ++take_calls_;
        // The highest range that fits, found through the size classes: in the
        // request's own class the highest range of at least the size (a class
        // spans a factor of two, so some of it may be too short), in the
        // classes above every range fits and the highest is the class's last.
        // The work is bounded by one class's population, not the free range
        // count: a walk from the top visited 2.5e9 ranges over the 142k takes
        // of the 768^3 leaf on 80 GB cards, where the level packs less and the
        // raw blocks leave 8224 free ranges (39.7 s of host time in the wave
        // buffers' allocation alone, the leaf 2.2x slower than on 40 GB).
        const int k = bin_of(size);
        bool found = false;
        size_t best = 0;
        for (auto it = top_bins_[static_cast<size_t>(k)].rbegin(); it != top_bins_[static_cast<size_t>(k)].rend(); ++it) {
            ++take_scans_;
            if (free_.find(*it)->second >= size) {
                best = *it;
                found = true;
                break;
            }
        }
        for (size_t b = static_cast<size_t>(k) + 1; b < top_bins_.size(); ++b) {
            if (top_bins_[b].empty()) continue;
            ++take_scans_;
            const size_t off = *top_bins_[b].rbegin();
            if (!found || off > best) {
                best = off;
                found = true;
            }
        }
        if (!found) return false;
        auto it = free_.find(best);
        const size_t start = it->first, length = it->second;
        remove_free(it);
        offset = start + length - size;
        if (length > size) add_free(start, length - size);
        return true;
    }

    int share_ranks_ = 1;
    size_t share_free_ = 0;
    char* base_ = nullptr;
    size_t capacity_ = 0;
    size_t used_ = 0;
    size_t peak_ = 0;
    size_t take_calls_ = 0, take_scans_ = 0;
    // Resident blocks (a level or longer: chiefly the retained solve factors
    // and the fill sources) versus per-wave transient blocks, and what each
    // held when the overall peak was reached.
    size_t resident_used_ = 0;
    std::unordered_map<size_t, ResidentKind> resident_kind_;  // offset -> what it is for
    std::array<size_t, static_cast<size_t>(ResidentKind::Count)> by_kind_{};
    std::array<size_t, static_cast<size_t>(ResidentKind::Count)> alloc_n_{}, alloc_b_{}, free_n_{}, free_b_{};
    std::array<size_t, static_cast<size_t>(ResidentKind::Count)> by_kind_at_peak_{};
    std::unordered_map<size_t, ResidentKind> transient_kind_;  // offset -> what it is for
    std::array<size_t, static_cast<size_t>(ResidentKind::Count)> by_tkind_{};
    std::array<size_t, static_cast<size_t>(ResidentKind::Count)> by_tkind_at_peak_{};
    size_t transient_used_ = 0;
    size_t resident_at_peak_ = 0;
    size_t transient_at_peak_ = 0;
    std::map<size_t, size_t> free_;                  // offset -> length, address order
    std::set<std::pair<size_t, size_t>> by_size_;    // (length, offset) of the same ranges
    std::array<std::set<size_t>, 64> top_bins_;      // offsets of the same ranges by size class
    // offset -> (length, resident)
    // BPACK_CHECK=heap: who allocated each live block, and who last freed each
    // offset.  A double free then names both sites instead of only the second,
    // which is what hypothesis-bisecting a 25 minute run cannot do.
    bool heap_trace_ = env::check("heap");
    std::unordered_map<size_t, const void*> alloc_site_;
    std::unordered_map<size_t, const void*> free_site_;
    std::unordered_map<size_t, std::pair<ResidentKind, size_t>> freed_kind_;
    std::unordered_map<size_t, std::pair<size_t, bool>> allocated_;   // offset -> length
    std::function<bool()> reclaimer_;
    int64_t reclaims_ = 0;
};

// The sub-buffers of a pass laid out back to back in one logical range
// (their offsets 256-aligned, in increasing order), placed in pieces of at
// most kPieceBytes, each its own transient heap block (a larger sub-buffer
// gets a piece of its own), instead of one block for the whole range: after
// a level the heap's free memory is in pieces (between the kept solve
// factors and the parent level's blocks), and one block of the pass's total
// may not fit although the total does.
class PiecedBuffer {
public:
    static constexpr size_t kPieceBytes = size_t{256} << 20;

    PiecedBuffer() = default;
    PiecedBuffer(const PiecedBuffer&) = delete;
    PiecedBuffer& operator=(const PiecedBuffer&) = delete;
    ~PiecedBuffer() { release(); }

    // `starts`: the sub-buffers' offsets; `end`: the end of the last one.
    void allocate(DeviceHeap& heap, const std::vector<size_t>& starts, size_t end,
                  ResidentKind kind = ResidentKind::Other) {
        release();
        heap_ = &heap;
        for (size_t j = 0; j < starts.size(); ++j) {
            const size_t stop = j + 1 < starts.size() ? starts[j + 1] : end;
            if (start_.empty() || (stop - start_.back() > kPieceBytes && starts[j] > start_.back())) {
                start_.push_back(starts[j]);
            }
        }
        for (size_t p = 0; p < start_.size(); ++p) {
            const size_t stop = p + 1 < start_.size() ? start_[p + 1] : end;
            base_.push_back(stop > start_[p] ? heap.alloc(stop - start_[p], kind) : nullptr);
        }
    }
    // Device address of a logical offset (inside a sub-buffer).
    char* at(size_t offset) const {
        const size_t p = static_cast<size_t>(std::upper_bound(start_.begin(), start_.end(), offset) - start_.begin()) - 1;
        return base_[p] + (offset - start_[p]);
    }
    // (later launches are ordered after the pass's reads: one stream)
    void release() {
        for (char* b : base_) {
            if (b != nullptr) heap_->free(b);
        }
        base_.clear();
        start_.clear();
    }

private:
    DeviceHeap* heap_ = nullptr;
    std::vector<size_t> start_;  // logical offset of each piece
    std::vector<char*> base_;    // its block
};

// ---------------------------------------------------------------------------
// Ranks sharing a GPU.  Collective over comm, at every call (a process may
// build matrices on several communicators, a sub-communicator first, and
// every rank of each must take part): the ranks of comm on a node that use
// the same GPU (by device UUID).  The first call sets the split of the GPU's
// memory, the device pool of each rank taking an equal part
// (DeviceHeap::set_device_share, before the pool is reserved).  Returns how
// many ranks of comm use this rank's GPU.
inline int share_device(MPI_Comm comm) {
    static bool applied = false;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm node = MPI_COMM_NULL;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
    int nrank = 0, nsize = 1;
    MPI_Comm_rank(node, &nrank);
    MPI_Comm_size(node, &nsize);
    Context::instance().activate();
    int dev = 0;
    check_cuda(cudaGetDevice(&dev), "cudaGetDevice");
    cudaDeviceProp prop;
    check_cuda(cudaGetDeviceProperties(&prop, dev), "cudaGetDeviceProperties");
    std::vector<char> uuids(16 * static_cast<size_t>(nsize));
    MPI_Allgather(prop.uuid.bytes, 16, MPI_CHAR, uuids.data(), 16, MPI_CHAR, node);
    int color = -1, count = 0;
    for (int r = 0; r < nsize; ++r) {
        if (std::memcmp(uuids.data() + 16 * r, prop.uuid.bytes, 16) != 0) continue;
        if (color < 0) color = r;
        ++count;
    }
    MPI_Comm dev_comm = MPI_COMM_NULL;
    MPI_Comm_split(node, color, nrank, &dev_comm);
    MPI_Barrier(dev_comm);  // (every context on the device exists)
    size_t free_bytes = 0, total_bytes = 0;
    check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
    unsigned long long f = free_bytes, fmin = free_bytes;
    MPI_Allreduce(&f, &fmin, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, dev_comm);
    if (count > 1 && !applied) DeviceHeap::instance().set_device_share(count, static_cast<size_t>(fmin));
    applied = true;
    MPI_Comm_free(&dev_comm);
    MPI_Comm_free(&node);
    return count;
}

// The exchange arena (DeviceHeap::exchange_arena()) of a run of nproc ranks
// with GPU-aware MPI, before the device pool takes its fraction of the free
// memory: BPACK_GPU_EXCHANGE_MB MB, default min(2560 MB, a sixteenth of the
// GPU), or `fraction` of the GPU when that is more (a component-owner CA
// level holds a colour's arrivals until its end: up to 2.4 GB per rank
// each way at the leaf of 320^3 on 64 ranks), divided among the `share`
// ranks that use the GPU (share_device).  MPI sends from and receives into
// it (registered once for MPI; buffers of the main pool made the transfers
// ~30x slower); a message that does not fit goes through the host.  No-op
// once reserved.
inline void reserve_exchange_arena(int nproc, int share, double fraction = 0.0) {
    DeviceHeap& arena = DeviceHeap::exchange_arena();
    if (nproc <= 1 || !device_exchange_enabled() || arena.initialized()) return;
    Context::instance().activate();
    size_t free_bytes = 0, total_bytes = 0;
    check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes), "cudaMemGetInfo");
    size_t mb = std::min<size_t>(2560, (total_bytes >> 20) / 16);
    if (fraction > 0.0) mb = std::max(mb, static_cast<size_t>(fraction * static_cast<double>(total_bytes >> 20)));
    mb /= static_cast<size_t>(std::max(share, 1));
    if (const char* env = std::getenv("BPACK_GPU_EXCHANGE_MB")) mb = static_cast<size_t>(std::atoll(env));
    arena.initialize_fixed(mb << 20);
}

// ---------------------------------------------------------------------------
// Several H2 operators in one process.  The device data a factorization or
// compression leaves for the solve and the matvec (device_solve_store(),
// device_matvec_store()) is kept per operator (its tree), so a process that
// alternates between operators (a Gaussian process with its covariance
// matrix and the derivative matrices of its gradient) finds their device
// data again instead of falling back to the host.  With a single operator
// nothing changes: a new factorization or compression of an operator
// replaces its own data.
//
// The operators of a process are assumed to share its ranks, which call
// their H2 operations in the same order: the use order below is then the
// same on every rank, and so are the evictions.  The least recently used
// operators' data is released collectively at the start of a build while the
// other operators hold more than half of some rank's heap (an eviction on
// one rank only would split the ranks between the device and host paths of
// that operator); an allocation that still does not fit fails.
struct OperatorContext {
    const void* active = nullptr;                         // the operator of the current H2 call
    uint64_t clock = 0;
    std::unordered_map<const void*, uint64_t> last_use;   // operator -> last use
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
    c.last_use[tree] = ++c.clock;
}

// Device bytes the stores hold for operator `tree`.
inline size_t operator_bytes(const void* tree) {
    size_t bytes = 0;
    for (const auto& held : operator_context().holders) bytes += held(tree);
    return bytes;
}

// Whether operators other than the active one hold device data.
inline bool other_operators_hold_data() {
    const OperatorContext& c = operator_context();
    for (const auto& entry : c.last_use) {
        if (entry.first != c.active && operator_bytes(entry.first) > 0) return true;
    }
    return false;
}

// The device data of operator `tree` is released (its destruction, on every rank).
inline void release_operator(const void* tree) {
    OperatorContext& c = operator_context();
    for (const auto& release : c.releasers) release(tree);
    c.last_use.erase(tree);
    if (c.active == tree) c.active = nullptr;
}

// Start of a factorization or compression of operator `tree` (collective
// over comm): the GPU sharing and the exchange arena (the first time, before
// the device pool; `exchange_fraction`: the arena's share of the GPU when
// more than the default, see reserve_exchange_arena), then activates the
// operator and makes room from the other operators' data.
inline void begin_operator_build(const void* tree, MPI_Comm comm, double exchange_fraction = 0.0) {
    int nproc = 1;
    MPI_Comm_size(comm, &nproc);
    reserve_exchange_arena(nproc, share_device(comm), exchange_fraction);
    activate_operator(tree);
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
