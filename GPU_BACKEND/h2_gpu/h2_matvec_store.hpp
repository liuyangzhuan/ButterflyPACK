#pragma once
// Device data of the compression-only H2 matvec (h2_matvec.hpp): per level,
// the box and block tables and the device copies of the coupling and near
// blocks the compression kept.  The matvec runs on the device only when
// every rank kept all of its levels (decided collectively at the end of the
// compression); otherwise the host blocks serve it.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "device_kernels.hpp"   // IndexList, PointTable
#include "h2_matvec_kernels.hpp"
#ifdef H2_HAVE_NVCOMP
#include <nvcomp/native/bitcomp.h>
#endif
#include <memory>
#include <vector>

#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// The compression keeps the blocks while the heap stays under this
// fraction of its capacity.
constexpr double kMatvecKeepFraction = 0.9;
// BPACK_CHECK=matvec: also run the host matvec and print the difference.
inline bool device_matvec_check() { return env::check("matvec"); }

// Messages of one step with the neighbor ranks, fixed by the compression:
// per peer, a send and a receive segment (points; a message is
// size * nrhs values at base * nrhs), and the spans that pack the sends.
// The transposed products (A^T x) run it backwards: what a rank receives
// in it, it sends (the sums of its contributions to those boxes), and what
// it sends, it receives and adds, a peer's spans at a time (span_begin).
struct H2Exchange {
    std::vector<int> peers;
    std::vector<int64_t> send_base, send_size, recv_base, recv_size;
    int64_t send_points = 0, recv_points = 0;
    const H2Span* d_spans = nullptr;
    int nspans = 0, max_len = 0;
    std::vector<int> span_begin;  // peer p's spans: span_begin[p] .. span_begin[p + 1]
};

// How a compression-only H2 keeps its blocks for the device matvec.  raw:
// on the device as evaluated (and on the host).  With H2_compress_factors
// (core/h2_block_codec.hpp): evaluated, none stored - the matvec evaluates
// each chunk of blocks into the store's scratch before the launches that
// read it (a GPU entry evaluator); compressed, as Bitcomp chunks on the
// device (and on the host) that the matvec decompresses into the scratch.
// A 900K space-time GP's derivative matrices held 20-25 GB of blocks a rank
// each on the host, raw.
enum class H2BlockStorage { raw = 0, evaluated = 1, compressed = 2 };

// A block the matvec puts into the scratch (evaluated or compressed): its
// place there, its shape, and its source.
struct H2MatBlock {
    size_t off = 0;                 // bytes from the scratch's start
    int m = 0, n = 0;
    IndexList rows{}, cols{};       // evaluated: the evaluator's lists (in the level's lists image)
    const char* packed = nullptr;   // compressed: its chunks end to end, on the device
    std::vector<uint32_t> strides;
};

// The blocks one group of launches reads, put into the scratch first.
struct H2MatChunk {
    int phase = 0;                  // 0: coupling, 1: leaf pairs (pairs i0..i1), 2: leaf near
    int i0 = 0, i1 = 0;
    // phases 0 and 2: a slice of every box's blocks (box, its first block in its list, count), so that
    // each launch has all the level's boxes (whole boxes a chunk left 4 of 256 level-5 boxes a launch)
    std::vector<int> slice_box, slice_first, slice_count;
    const H2Box* d_boxes = nullptr;  // the slices as boxes (the commit's table)
    int nboxes = 0;
    size_t bytes = 0;               // of the scratch, from its start
    std::vector<size_t> blocks;     // in H2MatvecLevel::mat_blocks
    // evaluated: the evaluator's items (the host copy for their count and sizes)
    std::shared_ptr<void> h_items;  // std::vector<EvalItemT<S>>
    const void* d_items = nullptr;
    int max_m = 0, max_n = 0;
    // compressed: one batch decompression
    const void** d_srcs = nullptr;
    void** d_dsts = nullptr;
    size_t nchunks = 0;
#ifdef H2_HAVE_NVCOMP
    bitcompHandle_t plan{};
    bool have_plan = false;
#endif
};

struct H2MatvecLevel {
    bool active = false;
    std::vector<H2Box> boxes;             // local boxes, in order
    int64_t points = 0;                   // level vector length (per right-hand side)
    int64_t q_points = 0;                 // compact skeleton vectors of the local boxes
    int max_n = 0, max_k = 0, max_r = 0, max_q_cols = 0, max_near_cols = 0;
    // the blocks of the boxes (coupling blocks of every box in order, then
    // near blocks): data, columns, source box, and for a near block of two
    // local boxes its partial vector (offset, or -1)
    std::vector<const void*> block_K;
    std::vector<int> block_cols;
    std::vector<int64_t> block_source;
    std::vector<int64_t> block_partial;
    int ncoupling = 0;
    // pairs of local boxes, one kept near block each (symmetric kernels)
    std::vector<H2Pair> pairs;
    int64_t partial_points = 0;
    int max_pair_cols = 0;
    const H2Pair* d_pairs = nullptr;
    // the plan (commit_device_matvec)
    H2Exchange coupling, near;            // remote skeleton vectors; remote leaf vectors
    bool handoff_device = false;          // to the parent level, on the device
    const H2Handoff* d_handoff = nullptr;
    int nhandoff = 0;
    std::vector<char*> allocations;       // blocks, factors, tables (heap)
    const H2Box* d_boxes = nullptr;
    const H2BlockRef* d_blocks = nullptr;
    // the transposed products (a nonsymmetric kernel's A^T x): the coupling's
    // targets (local boxes' compact skeleton vectors, remote boxes' send
    // entries) and the near blocks' (local boxes' vectors, remote boxes')
    const H2TTarget* d_tq_targets = nullptr;
    const H2TBlock* d_tq_blocks = nullptr;
    int ntq = 0;
    const H2TTarget* d_tx_targets = nullptr;
    const H2TBlock* d_tx_blocks = nullptr;
    int ntx = 0;
    double bytes = 0.0;
    // H2_compress_factors (storage other than raw): block_K and the pairs'
    // K hold scratch offsets + 1 until the commit points them into the
    // scratch; the chunks put the blocks there in each matvec
    std::vector<H2MatBlock> mat_blocks;
    std::vector<H2MatChunk> mat_chunks;
    const char* d_lists = nullptr;  // evaluated: the skeleton lists the items index (the meta base)
    PointTable point_table{};       // evaluated: the level's points (local and remote boxes)
};

struct DeviceMatvecStore {
    const void* tree = nullptr;
    bool building = false;   // the compression keeps blocks
    bool failed = false;     // this rank could not keep everything
    bool usable = false;     // all ranks kept everything
    bool transposable = false;  // the transposed products are planned (a nonsymmetric kernel)
    std::vector<H2MatvecLevel> levels;
    double bytes = 0.0;
    // H2_compress_factors (H2BlockStorage): the chunks' scratch (the largest
    // chunk of any level), and the evaluator an evaluated store calls
    H2BlockStorage storage = H2BlockStorage::raw;
    char* scratch = nullptr;
    size_t scratch_bytes = 0;
    std::shared_ptr<void> evaluator;

    bool empty() const { return levels.empty(); }
    void release() {
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& lv : levels) {
            for (char* p : lv.allocations) heap.free(p);
#ifdef H2_HAVE_NVCOMP
            for (auto& c : lv.mat_chunks) {
                if (c.have_plan) bitcompDestroyPlan(c.plan);
                c.have_plan = false;
            }
#endif
        }
        if (scratch != nullptr) heap.free(scratch);
        scratch = nullptr;
        scratch_bytes = 0;
        storage = H2BlockStorage::raw;
        evaluator.reset();
        levels.clear();
        tree = nullptr;
        building = false;
        failed = false;
        usable = false;
        transposable = false;
        bytes = 0.0;
    }
};

// The store of the active operator (one per operator, device_heap.hpp).
inline DeviceMatvecStore& device_matvec_store() {
    static std::unordered_map<const void*, DeviceMatvecStore>* const stores = [] {
        auto* s = new std::unordered_map<const void*, DeviceMatvecStore>;  // (never freed: see DeviceHeap)
        OperatorContext& c = operator_context();
        c.releasers.push_back([s](const void* tree) {
            auto it = s->find(tree);
            if (it == s->end()) return;
            // an evaluated store holds no blocks (points, lists, T, the
            // chunks' tables) and has no host copy to fall back to: an
            // eviction leaves it (the operator's destruction still frees it)
            if (operator_context().evicting && it->second.storage == H2BlockStorage::evaluated) return;
            it->second.release();
            s->erase(it);
        });
        c.holders.push_back([s](const void* tree) -> size_t {
            auto it = s->find(tree);
            return it == s->end() ? 0 : static_cast<size_t>(it->second.bytes);
        });
        return s;
    }();
    return (*stores)[operator_context().active];
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
