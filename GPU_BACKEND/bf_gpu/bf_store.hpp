#pragma once

#ifdef H2_HAVE_GPU

// What a butterfly keeps on the device between matvecs (GPU_BACKEND/bf_gpu).
//
// A device matvec is useful to its caller only if the factors and the sweep
// plans survive the call: STRUMPACK's Schur sampling applies F12, then F11's
// inverse, then F21, tens of thousands of times per front, so uploading the
// factors or rebuilding a plan per call would cost far more than the
// arithmetic.  The plans in particular cost two small collectives each to
// build (bf_exchange.hpp), which is why they are built once per block and
// orientation.
//
// Keyed on the Fortran matrixblock's address, which is stable for the life of
// the block.  An address can be reused after a block is freed, so every
// lookup re-reads the block's layout and discards the entry when it does not
// match: the check is a handful of integer reads against a cached fingerprint,
// which is nothing beside a matvec, and it turns a stale-store answer - wrong
// numbers, no error - into a re-upload.

#include <cstdint>
#include <map>
#include <memory>

#include "bf_device.hpp"
#include "bf_extract.hpp"
#include "bf_matvec.hpp"
#include "bf_ranges.hpp"
#include "bf_sweep_n.hpp"
#include "bf_upload.hpp"

namespace bf {
namespace gpu {

// Enough of a block's shape to notice that the address now holds a different
// block.
struct BfFingerprint {
    int level_butterfly = -1, level_half = -1;
    int M = 0, N = 0, M_loc = 0, N_loc = 0;
    int nlevels = 0, nblocks = 0;

    bool operator==(const BfFingerprint& o) const {
        return level_butterfly == o.level_butterfly && level_half == o.level_half && M == o.M &&
               N == o.N && M_loc == o.M_loc && N_loc == o.N_loc && nlevels == o.nlevels &&
               nblocks == o.nblocks;
    }
    bool operator!=(const BfFingerprint& o) const { return !(*this == o); }
};

struct BfResident {
    BfFingerprint print;
    BfDevice dev;
    BfBridge bridge;
    SweepPlan plan_n, plan_t;
    SweepScratch scratch;
    // the entry extraction's (bf_extract.hpp): the leaves, and the metadata
    // image of a call with its device copy, both reused between calls
    ExtractLeaves leaves;
    fmm::gpu::MetaBuilder extract_meta;
    fmm::gpu::DeviceBuffer extract_meta_device;

    ~BfResident() { release(); }
    void release() {
        plan_n.release();
        plan_t.release();
        scratch.release();
        extract_meta_device.release();
        bf::gpu::release(dev);
    }
    // the plan of one orientation, built on first use
    SweepPlan& plan(char op) { return (op == 'N' || op == 'n') ? plan_n : plan_t; }
};

// The process's resident butterflies.  Not thread safe: the callers are the
// matvec entry points, which STRUMPACK reaches from one thread per front.
inline std::map<const void*, std::unique_ptr<BfResident>>& bf_store() {
    static std::map<const void*, std::unique_ptr<BfResident>> store;
    return store;
}

inline BfResident* bf_store_find(const void* block, const BfFingerprint& print) {
    auto it = bf_store().find(block);
    if (it == bf_store().end()) return nullptr;
    if (it->second->print != print) {    // the address now holds another block
        bf_store().erase(it);
        return nullptr;
    }
    return it->second.get();
}

inline BfResident& bf_store_add(const void* block, const BfFingerprint& print) {
    auto entry = std::make_unique<BfResident>();
    entry->print = print;
    BfResident& ref = *entry;
    bf_store()[block] = std::move(entry);
    return ref;
}

inline void bf_store_release(const void* block) { bf_store().erase(block); }
inline void bf_store_release_all() { bf_store().clear(); }

// Device bytes the store holds, for the memory trace.
inline size_t bf_store_bytes() {
    size_t bytes = 0;
    for (const auto& e : bf_store()) {
        bytes += static_cast<size_t>(e.second->dev.elements) * sizeof(DataType);
    }
    return bytes;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
