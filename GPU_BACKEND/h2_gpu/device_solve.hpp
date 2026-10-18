#pragma once
// Device solve of the H2 Color factorization (GPU backend).
//
// The factors the host solve reads (T, X_SR, X_NR, the LU of X_RR and its
// pivots, index lists) are uploaded once, at the first solve after a
// factorization, into the device heap, where they stay for later solves.
// Whether the device solve runs is decided collectively: every rank must be
// able to hold its factors (and every level must be a Color level with an
// LU-factored X_RR), else all ranks keep the host solve.
//
// A solve runs each level's Color waves on the device:
//   forward   x_R -= T^T x_S;  x_S += X_SR x_R;  neighbors += X_NR x_R
//             (neighbor updates of a wave summed in box order and applied
//             after it; those of other ranks' boxes sent after each wave)
//   diagonal  x_R = X_RR^{-1} x_R, right after the level's forward waves
//             (it touches only redundant entries, final by then)
//   backward  x_R += X_SR^T x_S + X_NR^T x_N;  x_S -= T x_R, waves reversed,
//             other ranks' neighbor vectors refreshed before each wave (a
//             nonsymmetric matrix: x_R += X_RS x_S + X_RN x_N, its X_RS and
//             X_RN stored transposed as X_SR and X_NR are: SolveBox::xrs_t)
// The hand-off between levels (gather_skeleton_to_parent,
// scatter_solution_to_children, with the process reductions) and the root
// solve stay on the host: a level's vectors go to the device and back once
// per sweep.  The arithmetic matches the host solve up to rounding (its
// sums of neighbor updates follow the host's thread chunks).
//
// A replicated CA level (factored by the device's CA mode, GPU_CA_PLAN.md
// M4b) runs the host CA solve's scheme instead: its level vector holds the
// local boxes, then the ghosts (eliminated again here, with their factors),
// and a read-only area holds the assisting boxes; the waves are the CA
// groups x Morton parities, updates go to local or ghost boxes only, and the
// ghosts' and assisting boxes' vectors come from their owners through the
// host gather (gather_CA_boxes_solve): once before the forward sweep, after
// each backward group but the last.  No per-wave messages.

#ifdef H2_HAVE_GPU

#include "core/factorization.hpp"  // the tree and box types
#include "device_heap.hpp"
#include "gpu_runtime.hpp"
#include "h2_matvec_store.hpp"
#include "solve_kernels.hpp"

#include <algorithm>
#include <limits>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include "factor_compress.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// A level keeps its factors on the device during the factorization when an
// estimate of their bytes fits under this fraction of the heap on every
// rank, decided before its waves (factorization_driver.hpp
// decide_solve_keep); the rest of the heap is the level's working memory.
// The others are uploaded at the first solve.
constexpr double kSolveKeepFraction = 0.8;

// H2_compress_factors: the buffers a level's X_NR decompression rotates
// through (DeviceSolveLevel::xnr_scratch).  Two is enough to hide the
// decompression of the next wave behind the current wave's kernels and
// messages; each buffer holds the widest wave, so the scratch costs twice
// that per level.
constexpr size_t kXnrBuffers = 2;

// The second buffer is taken only while the heap keeps this much of itself
// free afterwards: the overlap is worth about a third of the decompression,
// never a level that then does not fit.
constexpr double kXnrSecondBufferHeadroom = 0.1;

// BPACK_CHECK=solve: also run the host solve (and multiply) and print the
// difference.
inline bool device_solve_check() { return env::check("solve"); }
// The host solve of a check runs with the device solve suspended.
inline bool& device_solve_suspended() {
    static bool suspended = false;
    return suspended;
}

// One level of the device solve (on this rank): factors, and the launch
// and message tables of its sweeps (built once, at the first solve).
struct DeviceSolveLevel {
    bool active = false;
    int level = -1;
    bool streamed = false;   // (H2 8x) its kept factors come by wave: the diagonal solves run with each wave
    std::vector<char*> blocks;          // factors and tables (heap)
    // H2_compress_factors: the waves' X_NR, decompressed just before their
    // sweep kernels into kXnrBuffers buffers in turn (a wave uses the buffer
    // of its parity), so a wave decompresses while the one before it is still
    // being read.
    char* xnr_scratch = nullptr;
    size_t xnr_scratch_bytes = 0;   // xnr_buffers * xnr_buffer_bytes
    size_t xnr_buffer_bytes = 0;    // one buffer: the widest wave
    size_t xnr_buffers = 1;         // kXnrBuffers, or 1 when the heap is tight
    std::vector<SolveXnrWave> xnr_waves;   // parallel to `waves`
    bool xnr_compressed = false;           // any wave has compressed X_NR
    // staged by the box walk, consumed by the table builders (which know the
    // wave order): per box, each compressed factor's bytes, uncompressed size
    // and the padded lengths of the chunks the compressed bytes are made of
    // (src null: that factor is read in place)
    struct PackedPart {
        const void* src = nullptr;
        size_t raw_bytes = 0;
        const uint32_t* strides = nullptr;
        uint32_t nchunks = 0;
    };
    std::vector<PackedPart> pending_xnr, pending_t, pending_xsr;
    std::vector<SolveBox> boxes;        // local boxes, device pointers
    std::vector<SolveSlot> slots;       // their one-hop slots
    std::vector<int64_t> slot_morton;   // the neighbor of each slot
    int64_t points = 0;                 // level vector length (per right-hand side)
    int64_t work_points = 0;            // largest wave's X_NR rows
    int max_n = 0, max_r = 0;
    // other ranks' boxes one hop from local ones (multi-rank levels)
    std::vector<int> peers;
    struct Ghost { int64_t morton; int owner; int n; int64_t offset; };
    std::vector<Ghost> ghosts;
    std::unordered_map<int64_t, int> ghost_of;
    int64_t ghost_points = 0;
    std::vector<std::vector<int>> readers;  // per local box: peers (indices) that read its vector
    // A Color wave.  Message sizes are per peer, in points.
    struct Wave {
        int order0 = 0, count = 0;           // its boxes (forward order; reversed for the backward sweep)
        int max_ntot = 0;                    // largest X_NR row count of its boxes
        int diag0 = 0, ndiag = 0;            // its local boxes' diagonal solves (L.d_diag), run with the wave
        // H2 8x (2026-10-10): the wave's boxes in segments whose raw X_NR
        // fits one buffer of the solve's scratch (first, count in the wave's
        // forward order; one segment when the wave fits), seg0 indexing
        // L.xnr_waves - at the 384^3 leaf a wave's X_NR is 9 GB raw
        struct Seg { int first = 0, count = 0; };
        std::vector<Seg> segs;
        size_t seg0 = 0;
        int item0 = 0, nitems = 0, max_count = 0;        // forward accumulation (local, outbox)
        int in_item0 = 0, in_nitems = 0, in_max_count = 0;  // updates from other ranks
        std::vector<int64_t> out, in;        // forward messages
        // other ranks' vectors refreshed before this wave: ref[0] for waves in
        // reverse (the solve's backward sweep), ref[1] in order (the multiply's
        // forward sweep)
        struct Refresh {
            int send0 = 0, nsend = 0, send_max_n = 0;
            int recv0 = 0, nrecv = 0, recv_max_n = 0;
            std::vector<int64_t> send, recv;
        };
        Refresh ref[2];
    };
    std::vector<Wave> waves;
    // a replicated CA level: the level vector holds the local boxes (first
    // local_points), then the ghosts (ghost_box: their Morton indices, in
    // box order after the local ones); `ghosts` above are the assisting boxes
    // (read only); group_end[g]: the wave after CA group g
    bool ca = false;
    bool host = false;  // a CA level whose factors did not fit: its steps run on the host (host_level)
    int64_t local_points = 0;
    std::vector<int64_t> ghost_box;
    std::vector<int> group_end;
    // its halo exchanges in device memory (build_ca_halo), with `peers`;
    // kinds: 0 the ghosts, 1 the ghosts and the assisting boxes, 2 the
    // assisting boxes; the copies in d_copies
    struct CaHalo {
        std::vector<int64_t> out, in;  // points per peer
        int64_t out_total = 0, in_total = 0;
        int send0 = 0, nsend = 0, send_max_n = 0;  // level vector -> message
        int vec0 = 0, nvec = 0, vec_max_n = 0;     // message -> level vector (ghosts)
        int area0 = 0, narea = 0, area_max_n = 0;  // message -> read-only area (assisting boxes)
    };
    bool device_halo = false;
    CaHalo halo[3];
    int ndiag = 0;
    int64_t out_max = 0, in_max = 0, send_max = 0, recv_max = 0;  // largest messages (points, all peers)
    // device tables
    const SolveBox* d_boxes = nullptr;
    const SolveSlot* d_slots = nullptr;
    const int* d_order = nullptr;       // forward order, wave after wave
    const int* d_rorder = nullptr;      // backward order (each wave reversed)
    const int* d_diag = nullptr;
    const SolveAccum* d_items = nullptr;
    const SolvePart* d_parts = nullptr;
    const SolveCopy* d_copies = nullptr;
};

// Factors of a level kept on the device during the factorization (all ranks
// of the level agreed), consumed by the first solve's preparation.
struct KeptSolveBox {
    const void* T = nullptr;  // elements of the factorization's type
    const void* xsr = nullptr;
    const void* xnr = nullptr;
    // a nonsymmetric matrix: X_RS^T (k x r) and X_RN^T (ntot x r), the
    // backward sweep's factors laid out as X_SR and X_NR (SolveBox::xrs_t)
    const void* xrs_t = nullptr;
    const void* xrn_t = nullptr;
    const void* lu = nullptr;
    const int* ipiv = nullptr;
    int k = 0, r = 0, ntot = 0;
    // H2_compress_factors: xnr holds this many bytes of lossy-compressed
    // data instead of ntot * r elements (0: uncompressed, as before), and
    // was quantized with this delta.  The bytes are the box's chunks laid end
    // to end (kXnrChunkBytes of elements each, the last one short); their
    // padded lengths are xnr_nchunks entries from xnr_chunk0 of the level's
    // chunk_strides, which is how a chunk is found.
    size_t xnr_bytes = 0;
    double xnr_delta = 0.0;
    uint32_t xnr_chunk0 = 0;
    uint32_t xnr_nchunks = 0;
    // H2_compress_factors: T and X_SR likewise, k * r elements each, their
    // chunks found the same way (bytes 0: uncompressed, read in place).
    size_t t_bytes = 0, xsr_bytes = 0;
    uint32_t t_chunk0 = 0, t_nchunks = 0, xsr_chunk0 = 0, xsr_nchunks = 0;
    // H2_compress_factors: set between the wave that kept this box and the
    // transition that transposes its X_NR out of the fill source
    // (transpose_kept_xnr).  It must never still be set when the solve runs,
    // and a null xnr would otherwise read as "not eliminated".
    bool xnr_pending = false;
};
struct KeptSolveLevel {
    std::vector<char*> blocks;
    std::unordered_map<int64_t, KeptSolveBox> boxes;  // Morton index -> factors
    std::vector<uint32_t> chunk_strides;  // every kept box's chunks, in the order they were kept
    double bytes = 0.0;
    // Host staging.  Once a level's transition is done nothing of the
    // factorization reads its store again, so its blocks go to pinned host
    // memory (stage_kept_level, at the start of the next level) and come
    // back at the first solve (upload_kept_level).  The boxes keep the device
    // addresses the blocks had and are rebased on the upload; a read in
    // between goes through kept_fetch.  The kept-LU logdet is summed before
    // the LU leaves (device_logdet_kept_level keeps it here).
    struct StagedBlock { uintptr_t was = 0; size_t off = 0, bytes = 0; char* host = nullptr; char* home = nullptr; };
    std::vector<StagedBlock> staged;  // by `was`; empty: the blocks are on the device
    // the staged bytes: chunks of the pinned pool, as a rule the ones the
    // level's stores just gave back (H2 host memory, 2026-10-10; pageable
    // memory before, as pinning anew ran at 1.9 GB/s)
    std::unique_ptr<PinnedArena> host;
    size_t staged_bytes = 0;
    // H2 8x (2026-10-11): a staged level the solve streams instead of
    // keeping resident.  Its blocks are numbered by the first solve wave
    // that reads them and each takes a slot of the store's ring (number mod
    // slots); a box's pointers are rebased to its blocks' slots once, at
    // the set-up (the staged form kept in boxes_staged), and a wave's
    // blocks are copied into their slots from the host chunks ahead of its
    // step (kept_wave_resident), the next wave's under this one's kernels.
    bool streamed = false;   // (a wave's blocks come into the ring ahead of its step; the diagonal solves run per wave)
    std::unordered_map<int64_t, KeptSolveBox> boxes_staged;   // the boxes with their staged-form pointers
    std::vector<std::vector<int>> wave_blocks;   // per solve wave: the staged blocks it reads (ids into `staged`)
    std::vector<size_t> block_off;               // per staged block: its place in the ring (bytes)
    size_t ring_bytes = 0;                       // the ring this level needs: no two blocks of a wave or of consecutive waves overlap
    std::vector<double> logdet_per_box;  // (log|det|, arg) per kept LU, in local box order
    int64_t logdet_boxes = -1;        // -1: not summed yet
    bool on_host() const { return host != nullptr; }
};

inline void rebase_kept_level(KeptSolveLevel& kl, const std::vector<char*>& targets);  // (below: a streamed level's pointers)

// The X_NR scratch of the solves (H2_compress_factors): one per process,
// shared by the set-ups of its operators.  They solve one at a time on the
// context's stream, and a solve decompresses every segment it reads into the
// scratch anew (DeviceSolveRun's buffer state lives for one solve), so one
// scratch serves them all; a scratch per operator held kXnrBuffers * 2 GB
// each (a quarter of the card's room), 16.7 GB for 8 diagonal blocks of a
// block-diagonal preconditioner on one GPU.  Allocated by the first set-up
// that decompresses, released with the last set-up that refers to it.  Its
// buffers are sized by that set-up's room; an operator with a box they cannot
// hold (a later, larger one, the first set up on a fuller card) gets a
// scratch of its own, sized by its own room as before the scratch was
// shared.  The shared one is not regrown: every set-up's decompression plans
// hold addresses in it.
struct SharedXnrScratch {
    char* ptr = nullptr;
    size_t bytes = 0, segment = 0;
    int users = 0;
};
inline SharedXnrScratch& shared_xnr_scratch() {
    static SharedXnrScratch scratch;
    return scratch;
}
inline void release_shared_xnr_scratch() {
    SharedXnrScratch& s = shared_xnr_scratch();
    if (s.users > 0 && --s.users == 0 && s.ptr != nullptr) {
        DeviceHeap::instance().free(s.ptr);
        s = SharedXnrScratch{};
    }
}

struct DeviceSolveStore {
    std::map<int, KeptSolveLevel> kept;  // level -> factors kept during the factorization
    std::set<int> ca_levels;             // replicated CA levels factored on the device (ghost factors kept or on the host)
    // their assisting boxes (points, skeleton positions), recorded after the
    // post-elimination assisting gather (the level clears its own copy)
    struct AssistingBox { int n = 0; std::vector<int> skeleton; };
    std::map<int, std::unordered_map<int64_t, AssistingBox>> ca_assisting;
    const void* tree = nullptr;
    bool decided = false;
    bool usable = false;
    std::vector<DeviceSolveLevel> levels;
    double bytes = 0.0;
    // H2 8x: the device home the streamed levels share (the largest of
    // them), and what the last sweeps streamed
    char* stream_block = nullptr;
    size_t stream_bytes = 0;
    int stream_levels = 0;
    // H2 8x: one X_NR scratch for every level's sweeps, kXnrBuffers buffers
    // of a segment each (the levels sweep one at a time)
    char* xnr_scratch = nullptr;
    size_t xnr_scratch_bytes = 0, xnr_segment_bytes = 0;
    bool xnr_scratch_own = false;  // this operator's own, not the process's (SharedXnrScratch)
    double stream_seconds = 0.0, stream_gb = 0.0;
    struct ResidentBlock { int level; int block; size_t pos, bytes; };
    std::vector<ResidentBlock> resident;             // what the ring holds now (a copy evicts what it overlaps)
    PinnedBuffer host_vec, host_out, host_in;  // staging of the solves (kept: pinning is slow)
    // the kept factors into the operator's host boxes (materialize_host_factors,
    // set by the factorization that kept them): before an eviction
    std::function<void(DeviceSolveStore&)> spill;
    bool empty() const { return levels.empty() && kept.empty(); }

    // A preparation's tables and uploaded factors (the kept factors stay).
    void release_prepared() {
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& lv : levels) {
            for (char* b : lv.blocks) heap.free(b);
        }
        levels.clear();
        for (auto& [level, kl] : kept) {
            if (kl.streamed) rebase_kept_level(kl, {});  // (its pointers back to the staged form)
        }
        if (stream_block != nullptr) heap.free(stream_block);
        stream_block = nullptr;
        stream_bytes = 0;
        stream_levels = 0;
        resident.clear();
        if (xnr_scratch != nullptr) {
            if (xnr_scratch_own) heap.free(xnr_scratch);
            else release_shared_xnr_scratch();
        }
        xnr_scratch = nullptr;
        xnr_scratch_bytes = xnr_segment_bytes = 0;
        xnr_scratch_own = false;
        tree = nullptr;
        decided = false;
        usable = false;
        bytes = 0.0;
    }
    void release() {
        release_prepared();
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& [level, kl] : kept) {
            for (char* b : kl.blocks) heap.free(b);
            kl.host.reset();
        }
        kept.clear();
        ca_levels.clear();
        ca_assisting.clear();
        tree = nullptr;
        decided = false;
        usable = false;
        bytes = 0.0;
    }
};

// Staged blocks keep a page-aligned start, on the host and once uploaded.
constexpr size_t kStagedAlign = size_t{64} << 10;

// The staged block holding `p` (an address a kept block had on the device),
// or null.
inline const KeptSolveLevel::StagedBlock* staged_block_of(const KeptSolveLevel& kl, const void* p) {
    const uintptr_t a = reinterpret_cast<uintptr_t>(p);
    auto it = std::upper_bound(kl.staged.begin(), kl.staged.end(), a,
                               [](uintptr_t v, const KeptSolveLevel::StagedBlock& b) { return v < b.was; });
    if (it == kl.staged.begin()) return nullptr;
    --it;
    return a < it->was + it->bytes ? &*it : nullptr;
}

// `bytes` of a kept factor into host memory: from the device, or from the
// staged copy when the level is on the host.
inline void kept_fetch(const KeptSolveLevel& kl, const void* src, size_t bytes, void* dst, const char* what) {
    if (bytes == 0) return;
    if (!kl.on_host()) {
        check_cuda(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost), what);
        return;
    }
    // (a streamed level: its callers read the boxes' staged-form pointers, boxes_staged)
    const KeptSolveLevel::StagedBlock* b = staged_block_of(kl, src);
    if (b == nullptr || reinterpret_cast<uintptr_t>(src) + bytes > b->was + b->bytes) {
        throw std::runtime_error(std::string("H2 GPU: ") + what + ": not in the level's staged kept store");
    }
    std::memcpy(dst, b->host + (reinterpret_cast<uintptr_t>(src) - b->was), bytes);
}

// A level's kept blocks to host memory, and off the device: chunks of the
// pinned pool, which as a rule are the ones the level's stores gave back
// (pinned already; pinning anew took 1.2 s for 2.3 GB where the copy took
// 0.09, so pageable memory was used before the pool).  The allocation and
// the copy are timed apart.
// Every box's pointers are checked to fall in a block, and the first bytes
// of each block are compared after the copy, so a wrong length or offset
// fails here and not in the solve.
inline void stage_kept_level(KeptSolveLevel& kl, cudaStream_t stream, double& alloc_seconds, double& copy_seconds) {
    using clock = std::chrono::steady_clock;
    if (kl.on_host() || kl.blocks.empty()) return;
    DeviceHeap& heap = DeviceHeap::instance();
    std::vector<KeptSolveLevel::StagedBlock> staged;
    staged.reserve(kl.blocks.size());
    size_t total = 0;
    for (char* b : kl.blocks) {
        const size_t n = heap.bytes_of(b);
        if (n == 0) throw std::runtime_error("H2 GPU: a kept block is not on the heap");
        staged.push_back(KeptSolveLevel::StagedBlock{reinterpret_cast<uintptr_t>(b), total, n});
        total = align_up(total + n, kStagedAlign);  // (each block keeps a page-aligned start)
    }
    std::vector<KeptSolveLevel::StagedBlock> by_was = staged;
    std::sort(by_was.begin(), by_was.end(),
              [](const KeptSolveLevel::StagedBlock& a, const KeptSolveLevel::StagedBlock& b) { return a.was < b.was; });
    for (size_t i = 1; i < by_was.size(); ++i) {  // the blocks must not overlap or repeat (the lookup needs it)
        if (by_was[i].was < by_was[i - 1].was + by_was[i - 1].bytes) {
            throw std::runtime_error("H2 GPU: kept blocks overlap: " + std::to_string(by_was[i - 1].was) + "+" +
                                     std::to_string(by_was[i - 1].bytes) + " and " + std::to_string(by_was[i].was) +
                                     "+" + std::to_string(by_was[i].bytes) + " (" + std::to_string(by_was.size()) +
                                     " blocks)");
        }
    }
    {
        KeptSolveLevel probe;  // (the lookup over the sorted list, before the copy)
        probe.staged = by_was;
        auto inside = [&](const void* p, size_t bytes, int64_t morton, const char* what) {
            if (p == nullptr) return;
            const KeptSolveLevel::StagedBlock* b = staged_block_of(probe, p);
            if (b == nullptr || reinterpret_cast<uintptr_t>(p) + bytes > b->was + b->bytes) {
                throw std::runtime_error(std::string("H2 GPU: the kept ") + what + " of box " + std::to_string(morton) +
                                         " is not inside the level's kept blocks");
            }
        };
        for (const auto& [morton, kb] : kl.boxes) {
            inside(kb.T, 0, morton, "T");
            inside(kb.xsr, 0, morton, "X_SR");
            inside(kb.xnr, kb.xnr_bytes, morton, "X_NR");
            inside(kb.xrs_t, 0, morton, "X_RS");
            inside(kb.xrn_t, 0, morton, "X_RN");
            inside(kb.lu, 0, morton, "LU");
            inside(kb.ipiv, static_cast<size_t>(kb.r) * sizeof(int), morton, "pivots");
        }
    }
    auto t0 = clock::now();
    auto host = std::make_unique<PinnedArena>();
    for (KeptSolveLevel::StagedBlock& b : staged) b.host = host->alloc(std::max<size_t>(b.bytes, 1));
    by_was = staged;
    std::sort(by_was.begin(), by_was.end(),
              [](const KeptSolveLevel::StagedBlock& a, const KeptSolveLevel::StagedBlock& b) { return a.was < b.was; });
    alloc_seconds += std::chrono::duration<double>(clock::now() - t0).count();
    t0 = clock::now();
    check_cuda(cudaStreamSynchronize(stream), "kept factors: to the host");  // (the last writes of the blocks)
    for (size_t i = 0; i < kl.blocks.size(); ++i) {
        check_cuda(cudaMemcpy(staged[i].host, kl.blocks[i], staged[i].bytes, cudaMemcpyDeviceToHost),
                   "kept factors: to the host");
    }
    copy_seconds += std::chrono::duration<double>(clock::now() - t0).count();
    {   // the copy against the device at each block's start, and the bytes against the store's accounting
        double summed = 0.0;
        for (size_t i = 0; i < kl.blocks.size(); ++i) {
            summed += static_cast<double>(staged[i].bytes);
            char head[256];
            const size_t n = std::min<size_t>(256, staged[i].bytes);
            check_cuda(cudaMemcpy(head, kl.blocks[i], n, cudaMemcpyDeviceToHost), "kept factors: check");
            if (std::memcmp(head, staged[i].host, n) != 0) {
                throw std::runtime_error("H2 GPU: staged kept block " + std::to_string(i) + " differs from the device's");
            }
        }
        if (summed + 1.0 < kl.bytes) {
            throw std::runtime_error("H2 GPU: the staged kept blocks hold fewer bytes than the store counts");
        }
    }
    for (char* b : kl.blocks) heap.free(b);
    kl.blocks.clear();
    kl.staged = std::move(by_was);
    kl.host = std::move(host);
    kl.staged_bytes = total;
}

// H2 8x (2026-10-11): the staged block a kept pointer falls in, by the
// block's current base (its slot while streamed, else the device address
// it had).  Two blocks of a streamed level may share a slot, so this is
// only asked about staged-form pointers (boxes_staged) or at the rebase.
inline const KeptSolveLevel::StagedBlock* kept_block_of(const KeptSolveLevel& kl, const void* p, size_t* which) {
    const uintptr_t a = reinterpret_cast<uintptr_t>(p);
    for (size_t i = 0; i < kl.staged.size(); ++i) {
        const KeptSolveLevel::StagedBlock& b = kl.staged[i];
        if (a >= b.was && a < b.was + b.bytes) {
            if (which != nullptr) *which = i;
            return &b;
        }
    }
    return nullptr;
}

// A staged level's box pointers moved from the staged form (the device
// addresses the blocks had, `was`) to the slots `targets[b]` of the ring
// (empty: back to the staged form, from boxes_staged).  The bytes are not
// moved: kept_wave_resident copies a wave's blocks in ahead of its step.
inline void rebase_kept_level(KeptSolveLevel& kl, const std::vector<char*>& targets) {
    if (targets.empty()) {
        if (kl.streamed) kl.boxes = kl.boxes_staged;
        kl.boxes_staged.clear();
        for (auto& b : kl.staged) b.home = nullptr;
        kl.streamed = false;
        return;
    }
    kl.boxes_staged = kl.boxes;
    auto move = [&](const void* p) -> const char* {
        if (p == nullptr) return nullptr;
        size_t i = 0;
        const KeptSolveLevel::StagedBlock* b = kept_block_of(kl, p, &i);
        if (b == nullptr) throw std::runtime_error("H2 GPU: a kept factor is outside its level's staged blocks (stream)");
        return targets[i] + (reinterpret_cast<uintptr_t>(p) - b->was);
    };
    for (auto& [morton, kb] : kl.boxes) {
        kb.T = move(kb.T);
        kb.xsr = move(kb.xsr);
        kb.xnr = move(kb.xnr);
        kb.xrs_t = move(kb.xrs_t);
        kb.xrn_t = move(kb.xrn_t);
        kb.lu = move(kb.lu);
        kb.ipiv = reinterpret_cast<const int*>(move(kb.ipiv));
    }
    for (size_t i = 0; i < kl.staged.size(); ++i) kl.staged[i].home = targets[i];
    kl.streamed = true;
}

// A wave's blocks into their slots on `stream`, those not there already
// (the ring remembers what each slot holds).  `sync`: wait for the copies.
inline void kept_wave_resident(DeviceSolveStore& store, int level, KeptSolveLevel& kl, size_t w, cudaStream_t stream,
                               bool sync) {
    if (!kl.streamed || w >= kl.wave_blocks.size()) return;
    for (int b : kl.wave_blocks[w]) {
        bool there = false;
        for (const auto& r : store.resident) {
            if (r.level == level && r.block == b) { there = true; break; }
        }
        if (there) continue;
        const KeptSolveLevel::StagedBlock& sb = kl.staged[static_cast<size_t>(b)];
        const size_t pos = kl.block_off[static_cast<size_t>(b)];
        // (the copy evicts whatever overlaps its range)
        store.resident.erase(std::remove_if(store.resident.begin(), store.resident.end(), [&](const auto& r) {
                                 return r.pos < pos + sb.bytes && pos < r.pos + r.bytes;
                             }), store.resident.end());
        check_cuda(cudaMemcpyAsync(sb.home, sb.host, sb.bytes, cudaMemcpyHostToDevice, stream), "streamed kept factors");
        store.resident.push_back({level, b, pos, sb.bytes});
        store.stream_gb += static_cast<double>(sb.bytes) / 1e9;
    }
    if (sync) check_cuda(cudaStreamSynchronize(stream), "streamed kept factors");
}

namespace solve_detail {
template<typename CoordType, typename DataType>
std::vector<std::vector<int>> color_waves(const TreeLevel<CoordType, DataType>& lvl, int dimension);
}
// The solve's waves of a level as lists of Morton indices, in the order
// the sweeps take them: the Color binning (color_waves) or the CA one (the
// groups x parities of build_ca_solve_tables).
template<typename CoordType, typename DataType>
std::vector<std::vector<int64_t>> solve_wave_mortons(ParallelTree<CoordType, DataType>* tree,
                                                     const TreeLevel<CoordType, DataType>& lvl, bool replicated) {
    std::vector<std::vector<int64_t>> waves;
    if (!replicated) {
        for (const auto& bin : solve_detail::color_waves(lvl, tree->dimension)) {
            std::vector<int64_t> m;
            m.reserve(bin.size());
            for (int b : bin) m.push_back(lvl.local_morton_start + b);
            waves.push_back(std::move(m));
        }
        return waves;
    }
    const int num_waves = 1 << tree->dimension;
    std::vector<std::vector<int64_t>> groups;
    {
        std::vector<int64_t> interior;
        for (int64_t local_idx : lvl.interior_id) interior.push_back(local_idx + lvl.local_morton_start);
        if (lvl.num_active_processes == 1) {
            std::vector<int64_t> boundary;
            for (int64_t local_idx : lvl.boundary_id) boundary.push_back(local_idx + lvl.local_morton_start);
            groups.push_back(std::move(boundary));
        } else {
            std::vector<int64_t> blue_orange = lvl.blue;
            blue_orange.insert(blue_orange.end(), lvl.orange.begin(), lvl.orange.end());
            groups.push_back(std::move(blue_orange));
            if (tree->dimension == 3) groups.push_back(lvl.purple);
            groups.push_back(lvl.green);
        }
        groups.push_back(std::move(interior));
    }
    for (const auto& group : groups) {
        std::vector<std::vector<int64_t>> bins(static_cast<size_t>(num_waves));
        for (int64_t m : group) bins[static_cast<size_t>(m & (num_waves - 1))].push_back(m);
        for (auto& bin : bins) waves.push_back(std::move(bin));
    }
    return waves;
}

// The streaming plan of a staged level over the solve's waves: the blocks
// each wave reads, the blocks numbered by first use, a slot each (number
// mod slots) with the ring grown until no two blocks of a wave or of two
// consecutive waves share a slot (the prefetch is one wave deep).
// Consecutive among the waves that read blocks: an owner level's home order
// leaves a rank waves without boxes, and the prefetch passes over them.
inline void plan_kept_streaming(KeptSolveLevel& kl, const std::vector<std::vector<int64_t>>& waves) {
    const size_t nb = kl.staged.size();
    std::vector<int> first_use(nb, -1);
    std::vector<std::vector<int>> wave_blocks(waves.size());
    for (size_t w = 0; w < waves.size(); ++w) {
        std::vector<int> used;
        for (int64_t m : waves[w]) {
            auto it = kl.boxes.find(m);
            if (it == kl.boxes.end()) continue;
            const KeptSolveBox& kb = it->second;
            for (const void* ptr : {kb.T, kb.xsr, kb.xnr, kb.xrs_t, kb.xrn_t, kb.lu, static_cast<const void*>(kb.ipiv)}) {
                if (ptr == nullptr) continue;
                size_t i = 0;
                if (kept_block_of(kl, ptr, &i) == nullptr) {
                    throw std::runtime_error("H2 GPU solve: a kept factor is outside its level's staged blocks (plan)");
                }
                used.push_back(static_cast<int>(i));
            }
        }
        std::sort(used.begin(), used.end());
        used.erase(std::unique(used.begin(), used.end()), used.end());
        for (int b : used) {
            if (first_use[static_cast<size_t>(b)] < 0) first_use[static_cast<size_t>(b)] = static_cast<int>(w);
        }
        wave_blocks[w] = std::move(used);
    }
    std::vector<int> ordered(nb);
    for (size_t i = 0; i < nb; ++i) ordered[i] = static_cast<int>(i);
    const int never = std::numeric_limits<int>::max();
    std::stable_sort(ordered.begin(), ordered.end(), [&](int a, int b) {
        const int fa = first_use[static_cast<size_t>(a)] < 0 ? never : first_use[static_cast<size_t>(a)];
        const int fb = first_use[static_cast<size_t>(b)] < 0 ? never : first_use[static_cast<size_t>(b)];
        return fa < fb;
    });
    // placed in the ring in that order, each block whole (a block that would
    // cross the end starts over at 0); the ring grows until no two distinct
    // blocks of a wave or of two consecutive waves overlap (the prefetch is
    // one wave deep).  The first try is twice the largest pair.
    std::vector<size_t> reading;   // the waves that read blocks, in order
    for (size_t w = 0; w < wave_blocks.size(); ++w) {
        if (!wave_blocks[w].empty()) reading.push_back(w);
    }
    size_t largest_pair = 0, largest_block = 0;
    for (size_t i = 0; i < reading.size(); ++i) {
        std::vector<int> blocks;
        for (size_t j = i; j < std::min(i + 2, reading.size()); ++j) {
            for (int b : wave_blocks[reading[j]]) blocks.push_back(b);
        }
        std::sort(blocks.begin(), blocks.end());
        blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
        size_t bytes = 0;
        for (int b : blocks) bytes += align_up(kl.staged[static_cast<size_t>(b)].bytes, kStagedAlign);
        largest_pair = std::max(largest_pair, bytes);
    }
    for (const auto& b : kl.staged) largest_block = std::max(largest_block, align_up(b.bytes, kStagedAlign));
    std::vector<size_t> off(nb, 0);
    auto place = [&](size_t ring) {
        size_t at = 0;
        for (size_t r = 0; r < nb; ++r) {
            const size_t b = static_cast<size_t>(ordered[r]);
            const size_t bytes = align_up(kl.staged[b].bytes, kStagedAlign);
            if (at + bytes > ring) at = 0;
            off[b] = at;
            at += bytes;
        }
    };
    auto overlaps = [&]() {
        for (size_t i = 0; i < reading.size(); ++i) {
            std::vector<int> blocks;
            for (size_t j = i; j < std::min(i + 2, reading.size()); ++j) {
                for (int b : wave_blocks[reading[j]]) blocks.push_back(b);
            }
            std::sort(blocks.begin(), blocks.end());
            blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
            for (size_t x = 0; x < blocks.size(); ++x) {
                const size_t bx = static_cast<size_t>(blocks[x]);
                for (size_t y = x + 1; y < blocks.size(); ++y) {
                    const size_t by = static_cast<size_t>(blocks[y]);
                    if (off[bx] < off[by] + kl.staged[by].bytes && off[by] < off[bx] + kl.staged[bx].bytes) return true;
                }
            }
        }
        return false;
    };
    size_t ring = std::max(2 * largest_pair, largest_block);
    size_t total = 0;
    for (const auto& b : kl.staged) total += align_up(b.bytes, kStagedAlign);
    for (;;) {
        place(ring);
        if (!overlaps() || ring >= total) break;
        ring = std::min(total, ring + std::max(ring / 8, largest_block));
    }
    kl.ring_bytes = ring;
    kl.block_off = std::move(off);
    kl.wave_blocks = std::move(wave_blocks);
}

// A staged level back on the device, in one block; every box's pointers
// rebased to it.
inline void upload_kept_level(KeptSolveLevel& kl, cudaStream_t stream) {
    if (!kl.on_host()) return;
    DeviceHeap& heap = DeviceHeap::instance();
    char* raw = heap.alloc_resident(std::max<size_t>(kl.staged_bytes + kStagedAlign, 1), ResidentKind::KeptFactors);
    char* dev = reinterpret_cast<char*>(align_up(reinterpret_cast<uintptr_t>(raw), kStagedAlign));
    check_cuda(cudaStreamSynchronize(stream), "kept factors: upload");
    for (const KeptSolveLevel::StagedBlock& b : kl.staged) {
        check_cuda(cudaMemcpy(dev + b.off, b.host, b.bytes, cudaMemcpyHostToDevice), "kept factors: upload");
    }
    auto rebase = [&](const void* p) -> const char* {
        if (p == nullptr) return nullptr;
        const KeptSolveLevel::StagedBlock* b = staged_block_of(kl, p);
        if (b == nullptr) throw std::runtime_error("H2 GPU: a kept factor is outside its level's staged blocks");
        return dev + b->off + (reinterpret_cast<uintptr_t>(p) - b->was);
    };
    for (auto& [morton, kb] : kl.boxes) {
        kb.T = rebase(kb.T);
        kb.xsr = rebase(kb.xsr);
        kb.xnr = rebase(kb.xnr);
        kb.xrs_t = rebase(kb.xrs_t);
        kb.xrn_t = rebase(kb.xrn_t);
        kb.lu = rebase(kb.lu);
        kb.ipiv = reinterpret_cast<const int*>(rebase(kb.ipiv));
    }
    kl.host.reset();   // its chunks back to the pinned pool
    pinned_chunk_pool().trim();
    kl.staged.clear();
    kl.staged_bytes = 0;
    kl.blocks.push_back(raw);
}

// The store of the active operator (one per operator, device_heap.hpp).
inline DeviceSolveStore& device_solve_store() {
    static std::unordered_map<const void*, DeviceSolveStore>* const stores = [] {
        auto* s = new std::unordered_map<const void*, DeviceSolveStore>;  // (never freed: see DeviceHeap)
        OperatorContext& c = operator_context();
        c.releasers.push_back([s](const void* tree) {
            auto it = s->find(tree);
            if (it == s->end()) return;
            it->second.release();
            s->erase(it);
        });
        c.holders.push_back([s](const void* tree) -> size_t {
            auto it = s->find(tree);
            if (it == s->end()) return 0;
            double bytes = it->second.bytes;
            for (const auto& [level, kl] : it->second.kept) bytes += kl.bytes;
            return static_cast<size_t>(bytes);
        });
        c.spillers.push_back([s](const void* tree) {
            auto it = s->find(tree);
            if (it != s->end() && it->second.spill && !it->second.kept.empty()) it->second.spill(it->second);
        });
        return s;
    }();
    return (*stores)[operator_context().active];
}

// A factorization or compression replaces the factors: drop the device
// copies (the solve's factors, and the blocks of a compression-only H2).
inline void invalidate_device_solve() {
    device_solve_store().release();
    device_matvec_store().release();
}

// ---------------------------------------------------------------------------
// Messages with the one-hop neighbor ranks: a header of int64 and a payload
// of doubles per peer, in two rounds (sizes, then all messages at once).
struct PeerMessage {
    std::vector<int64_t> head;
    std::vector<double> data;
};

inline double exchange_peer_messages(MPI_Comm comm, const std::vector<int>& peers, const std::vector<PeerMessage>& out,
                                     std::vector<PeerMessage>& in, int tag) {
    const auto t0 = std::chrono::steady_clock::now();
    const size_t np = peers.size();
    in.assign(np, PeerMessage{});
    std::vector<uint64_t> sizes_out(2 * np), sizes_in(2 * np, 0);
    for (size_t i = 0; i < np; ++i) {
        sizes_out[2 * i] = out[i].head.size();
        sizes_out[2 * i + 1] = out[i].data.size();
    }
    std::vector<MPI_Request> reqs;
    reqs.reserve(4 * np);
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        MPI_Irecv(&sizes_in[2 * i], 2, MPI_UINT64_T, peers[i], tag, comm, &reqs.back());
    }
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        MPI_Isend(&sizes_out[2 * i], 2, MPI_UINT64_T, peers[i], tag, comm, &reqs.back());
    }
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.clear();
    for (size_t i = 0; i < np; ++i) {
        in[i].head.resize(sizes_in[2 * i]);
        in[i].data.resize(sizes_in[2 * i + 1]);
        if (!in[i].head.empty()) {
            reqs.emplace_back();
            MPI_Irecv(in[i].head.data(), static_cast<int>(in[i].head.size()), MPI_INT64_T, peers[i], tag + 1, comm,
                      &reqs.back());
        }
        if (!in[i].data.empty()) {
            reqs.emplace_back();
            MPI_Irecv(in[i].data.data(), static_cast<int>(in[i].data.size()), MPI_DOUBLE, peers[i], tag + 2, comm,
                      &reqs.back());
        }
    }
    for (size_t i = 0; i < np; ++i) {
        if (!out[i].head.empty()) {
            reqs.emplace_back();
            MPI_Isend(out[i].head.data(), static_cast<int>(out[i].head.size()), MPI_INT64_T, peers[i], tag + 1, comm,
                      &reqs.back());
        }
        if (!out[i].data.empty()) {
            reqs.emplace_back();
            MPI_Isend(out[i].data.data(), static_cast<int>(out[i].data.size()), MPI_DOUBLE, peers[i], tag + 2, comm,
                      &reqs.back());
        }
    }
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// Preparation (collective over the tree's ranks)
// ---------------------------------------------------------------------------
namespace solve_detail {

template<typename CoordType, typename DataType>
int owner_of(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl, int level, int64_t morton) {
    (void)tree;
    (void)level;
    return lvl.owner_rank_of(morton);  // (the level's region ranges)
}

// Color waves as the host solve bins them: boundary boxes by color, then
// interior boxes by color, Morton order within a wave.
template<typename CoordType, typename DataType>
std::vector<std::vector<int>> color_waves(const TreeLevel<CoordType, DataType>& lvl, int dimension) {
    const int num_colors = 1 << dimension;
    std::vector<std::vector<int>> waves(static_cast<size_t>(2 * num_colors));
    if (lvl.ca_home_factors) {
        // a component-owner level with its factors at home: the level's
        // elimination order (the schedule's sequence numbers), not parities
        waves.clear();
        for (size_t b = 0; b < lvl.local_boxes.size(); ++b) {
            const int64_t morton = lvl.local_morton_start + static_cast<int64_t>(b);
            const auto it = lvl.ca_home_order.find(morton);
            if (it == lvl.ca_home_order.end()) {
                throw std::runtime_error("device solve: box " + std::to_string(morton) + " of an owner level has no wave");
            }
            const size_t w = static_cast<size_t>(it->second);
            if (waves.size() <= w) waves.resize(w + 1);
            waves[w].push_back(static_cast<int>(b));
        }
        return waves;
    }
    for (size_t b = 0; b < lvl.local_boxes.size(); ++b) {
        const int64_t morton = lvl.local_morton_start + static_cast<int64_t>(b);
        const int color = static_cast<int>(morton & (num_colors - 1));
        waves[static_cast<size_t>(color + (lvl.local_boxes[b].on_boundary ? 0 : num_colors))].push_back(
            static_cast<int>(b));
    }
    return waves;
}

// Per-box byte layout of the factors in a device block.
struct BoxLayout {
    size_t T = 0, xsr = 0, xnr = 0, lu = 0, ipiv = 0, skel = 0, red = 0, bytes = 0;
    size_t xrs_t = 0, xrn_t = 0;   // (a nonsymmetric matrix)
    bool has_T = false, has_xsr = false, has_xrs = false;
};

}  // namespace solve_detail

namespace solve_detail {

// BPACK_TRACE=wave: what the per-wave X_NR decompression has to work with at
// this level - how many waves, how much each one holds, and the scratch the
// largest one sizes.  A wave far below the batch size that fills the GPU runs
// the decompressor at a fraction of its rate, so this is what decides whether
// several waves need batching together (GPU_CA_PLAN.md C1).
inline void report_xnr_wave_sizes(const DeviceSolveLevel& L, int level) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) return;
    std::vector<size_t> bytes;
    std::vector<size_t> blocks;
    size_t total = 0, boxes_total = 0;
    for (const auto& xw : L.xnr_waves) {
        if (xw.empty()) continue;
        size_t b = 0;
        boxes_total += xw.boxes.size();
        for (size_t i = 0; i < xw.boxes.size(); ++i) {
            const size_t bx = static_cast<size_t>(xw.boxes[i]);
            b += L.pending_xnr[bx].raw_bytes + L.pending_t[bx].raw_bytes + L.pending_xsr[bx].raw_bytes;
        }
        bytes.push_back(b);
        blocks.push_back(xw.h_srcs.size());  // the batch's blocks: chunks, not boxes
        total += b;
    }
    if (bytes.empty()) return;
    std::vector<size_t> sb = bytes, sk = blocks;
    std::sort(sb.begin(), sb.end());
    std::sort(sk.begin(), sk.end());
    const auto mib = [](size_t b) { return static_cast<double>(b) / (1024.0 * 1024.0); };
    std::cout << "H2_XNR_WAVES level " << level << ": " << bytes.size() << " waves, "
              << std::fixed << std::setprecision(2) << mib(total) << " MiB of X_NR, T and X_SR; per wave "
              << mib(sb.front()) << " / " << mib(sb[sb.size() / 2]) << " / " << mib(sb.back())
              << " MiB (min/med/max), " << sk.front() << " / " << sk[sk.size() / 2] << " / " << sk.back()
              << " blocks of " << boxes_total << " boxes; scratch " << mib(L.xnr_scratch_bytes) << " MiB in "
              << L.xnr_buffers << std::endl;
}

// H2_compress_factors: the per-wave batches that decompress a wave's X_NR
// into L.xnr_scratch just before its sweep kernels.  The scratch is sized for
// the largest wave, so only one wave's X_NR is ever uncompressed; both sweeps
// read X_NR, so a wave is decompressed twice per solve.  The batch plans are
// built once here and reused by every sweep of every solve.
inline void build_solve_xnr_waves(DeviceSolveLevel& L, const std::vector<int>& order, cudaStream_t stream,
                                  int level) {
#ifdef H2_HAVE_NVCOMP
    for (size_t w = 0; w < L.waves.size(); ++w) {   // (one segment a wave unless cut below)
        L.waves[w].segs.assign(1, {0, L.waves[w].count});
        L.waves[w].seg0 = w;
    }
    if (!L.xnr_compressed) {
        if (env::trace("wave")) {
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            size_t boxes_with_xnr = 0, on_heap = 0;
            double elems = 0.0;
            const DeviceHeap& heap = DeviceHeap::instance();
            for (const SolveBox& sb : L.boxes) {
                if (sb.xnr == nullptr) continue;
                ++boxes_with_xnr;
                if (heap.owns(sb.xnr)) ++on_heap;
                elems += static_cast<double>(sb.ntot) * static_cast<double>(sb.r);
            }
            if (rank == 0) {
                std::cout << "H2_XNR_WAVES level " << level << ": X_NR read in place (not compressed), "
                          << boxes_with_xnr << " of " << L.boxes.size() << " boxes have X_NR, " << on_heap
                          << " on the device heap; " << std::fixed << std::setprecision(2)
                          << elems * 8.0 / 1e9 << " GB at 8 bytes an element" << std::endl;
            }
        }
        return;
    }
    // A box's factors were compressed in chunks, so the batch has one entry
    // per chunk: its compressed bytes sit `stride` apart in the kept block,
    // and it decompresses kXnrChunkBytes further into the scratch than the one
    // before.  A box's compressed factors land one after another, T, X_SR,
    // then X_NR, each at its own offset (kNone: that factor is read in place
    // and its pointer is left alone); the box's entry for the pointer patch
    // is each factor's first chunk.
    // H2 8x: a wave is cut into segments whose raw bytes fit one buffer of
    // the store's shared scratch (store.xnr_segment_bytes); a segment is the
    // unit the sweeps decompress and launch.
    DeviceSolveStore& store = device_solve_store();
    if (store.xnr_scratch == nullptr || store.xnr_segment_bytes == 0) {
        throw std::runtime_error("H2 GPU solve: the X_NR scratch was not set up before the level's waves");
    }
    const size_t seg_cap = store.xnr_segment_bytes;
    L.xnr_waves.clear();
    const size_t chunk_bytes_raw = kXnrChunkBytes;
    constexpr size_t kNone = static_cast<size_t>(-1);
    struct BoxOffsets { size_t xnr, t, xsr; };
    std::vector<std::vector<BoxOffsets>> offsets;   // per segment
    size_t widest = 0;
    auto raw_of = [&](size_t b) {
        size_t r = 0;
        for (const DeviceSolveLevel::PackedPart* part : {&L.pending_t[b], &L.pending_xsr[b], &L.pending_xnr[b]}) {
            if (part->src != nullptr) r += align_up(part->raw_bytes);
        }
        return r;
    };
    for (size_t w = 0; w < L.waves.size(); ++w) {
        auto& W = L.waves[w];
        W.segs.clear();
        W.seg0 = L.xnr_waves.size();
        int first = 0;
        while (first < W.count) {
            int last = first;
            size_t seg_bytes = 0;
            while (last < W.count) {
                const size_t rb = raw_of(static_cast<size_t>(order[static_cast<size_t>(W.order0 + last)]));
                if (last > first && seg_bytes + rb > seg_cap) break;
                seg_bytes += rb;
                ++last;
            }
            SolveXnrWave xw;
            std::vector<BoxOffsets> segoff;
            size_t bytes = 0;
            auto place = [&](const DeviceSolveLevel::PackedPart& part, const char* what) {
                if (part.src == nullptr) return kNone;
                if (part.strides == nullptr || part.nchunks == 0) {
                    throw std::runtime_error(std::string("device solve: compressed ") + what +
                                             " without its chunk lengths");
                }
                const size_t at = bytes;
                const char* src = static_cast<const char*>(part.src);
                for (size_t c = 0; c < part.nchunks; ++c) {
                    xw.h_srcs.push_back(src);
                    xw.h_dsts.push_back(reinterpret_cast<void*>(at + c * chunk_bytes_raw));  // offset for now
                    src += part.strides[c];
                }
                bytes = align_up(bytes + part.raw_bytes);
                return at;
            };
            for (int i = first; i < last; ++i) {
                const size_t b = static_cast<size_t>(order[static_cast<size_t>(W.order0 + i)]);
                if (L.pending_xnr[b].src == nullptr && L.pending_t[b].src == nullptr &&
                    L.pending_xsr[b].src == nullptr) {
                    continue;
                }
                xw.boxes.push_back(static_cast<int>(b));
                BoxOffsets o{};
                o.t = place(L.pending_t[b], "T");
                o.xsr = place(L.pending_xsr[b], "X_SR");
                o.xnr = place(L.pending_xnr[b], "X_NR");
                segoff.push_back(o);
            }
            xw.raw_bytes = bytes;
            widest = std::max(widest, bytes);
            W.segs.push_back({first, last - first});
            L.xnr_waves.push_back(std::move(xw));
            offsets.push_back(std::move(segoff));
            first = last;
        }
    }
    if (widest == 0) { L.xnr_compressed = false; return; }
    if (widest > seg_cap) {
        throw std::runtime_error("H2 GPU solve: a box's compressed factors (" + std::to_string(widest >> 20) +
                                 " MB raw) exceed the X_NR scratch's buffer (" + std::to_string(seg_cap >> 20) + " MB)");
    }
    // The shared scratch: kXnrBuffers buffers of a segment each; a segment is
    // decompressed into the buffer of its parity while the previous one's
    // kernels and messages are still using the other (prepare_wave_xnr).
    L.xnr_buffers = kXnrBuffers;
    L.xnr_buffer_bytes = seg_cap;
    L.xnr_scratch = store.xnr_scratch;
    L.xnr_scratch_bytes = store.xnr_scratch_bytes;
    if (env::trace("wave")) report_xnr_wave_sizes(L, level);
    KeptSolveLevel* klp = nullptr;
    if (auto kit = store.kept.find(level); kit != store.kept.end() && kit->second.streamed) klp = &kit->second;
    for (size_t w = 0; w < L.waves.size(); ++w) {
      // (a streamed level: the wave's blocks in their slots while its plans
      // are made - the plan reads the compressed headers)
      if (klp != nullptr) kept_wave_resident(store, level, *klp, w, stream, true);
      for (size_t s = 0; s < L.waves[w].segs.size(); ++s) {
        const size_t g = L.waves[w].seg0 + s;
        auto& xw = L.xnr_waves[g];
        if (xw.empty()) continue;
        char* base = L.xnr_scratch + (g % L.xnr_buffers) * L.xnr_buffer_bytes;
        for (auto& d : xw.h_dsts) {
            d = base + reinterpret_cast<size_t>(d);  // offset -> address
        }
        auto address = [&](size_t off) { return off == kNone ? nullptr : static_cast<void*>(base + off); };
        xw.h_box_dsts.resize(xw.boxes.size());
        xw.h_box_t_dsts.resize(xw.boxes.size());
        xw.h_box_xsr_dsts.resize(xw.boxes.size());
        for (size_t i = 0; i < xw.boxes.size(); ++i) {
            xw.h_box_dsts[i] = address(offsets[g][i].xnr);
            xw.h_box_t_dsts[i] = address(offsets[g][i].t);
            xw.h_box_xsr_dsts[i] = address(offsets[g][i].xsr);
        }
        int code = 0;
        if (!build_solve_xnr_wave(xw, DeviceHeap::instance(), stream, &code)) {
            throw std::runtime_error("H2 GPU solve: could not plan segment " + std::to_string(g) + "'s X_NR decompression (" +
                                     std::to_string(xw.h_srcs.size()) + " chunks of " + std::to_string(xw.boxes.size()) +
                                     " boxes, bitcomp " + std::to_string(code) + ")");
        }
      }
    }
#else
    (void)L; (void)order; (void)stream;
#endif
}

// A level's tables on the device (one block): the other ranks' skeleton
// lists, slots, boxes, wave orders, diagonal list, accumulations, copies.
// Returns their bytes.
inline size_t upload_solve_tables(DeviceSolveLevel& L, const std::vector<int>& remote_lists,
                                  const std::unordered_map<int64_t, size_t>& remote_list_at,
                                  const std::vector<int>& order, const std::vector<int>& rorder,
                                  const std::vector<int>& diag, const std::vector<SolveAccum>& items,
                                  const std::vector<SolvePart>& parts, const std::vector<SolveCopy>& copies) {
    const size_t I = sizeof(int);
    const size_t t_lists = 0;
    const size_t t_slots = align_up(t_lists + remote_lists.size() * I);
    const size_t t_boxes = align_up(t_slots + L.slots.size() * sizeof(SolveSlot));
    const size_t t_order = align_up(t_boxes + L.boxes.size() * sizeof(SolveBox));
    const size_t t_rorder = align_up(t_order + order.size() * I);
    const size_t t_diag = align_up(t_rorder + rorder.size() * I);
    const size_t t_items = align_up(t_diag + diag.size() * I);
    const size_t t_parts = align_up(t_items + items.size() * sizeof(SolveAccum));
    const size_t t_copies = align_up(t_parts + parts.size() * sizeof(SolvePart));
    const size_t t_end = align_up(t_copies + copies.size() * sizeof(SolveCopy));
    char* tables = DeviceHeap::instance().alloc_resident(std::max<size_t>(t_end, 1));
    L.blocks.push_back(tables);
    for (size_t s = 0; s < L.slots.size(); ++s) {
        SolveSlot& sl = L.slots[s];
        if (sl.ghost && !sl.full) {
            sl.skel = reinterpret_cast<const int*>(tables + t_lists) + remote_list_at.at(L.slot_morton[s]);
        }
    }
    std::vector<char> image(t_end, 0);
    auto put = [&](size_t at, const void* data, size_t bytes) {
        if (bytes > 0) std::memcpy(image.data() + at, data, bytes);
    };
    put(t_lists, remote_lists.data(), remote_lists.size() * I);
    put(t_slots, L.slots.data(), L.slots.size() * sizeof(SolveSlot));
    put(t_boxes, L.boxes.data(), L.boxes.size() * sizeof(SolveBox));
    put(t_order, order.data(), order.size() * I);
    put(t_rorder, rorder.data(), rorder.size() * I);
    put(t_diag, diag.data(), diag.size() * I);
    put(t_items, items.data(), items.size() * sizeof(SolveAccum));
    put(t_parts, parts.data(), parts.size() * sizeof(SolvePart));
    put(t_copies, copies.data(), copies.size() * sizeof(SolveCopy));
    check_cuda(cudaMemcpy(tables, image.data(), t_end, cudaMemcpyHostToDevice), "solve tables");
    L.d_slots = reinterpret_cast<const SolveSlot*>(tables + t_slots);
    L.d_boxes = reinterpret_cast<const SolveBox*>(tables + t_boxes);
    L.d_order = reinterpret_cast<const int*>(tables + t_order);
    L.d_rorder = reinterpret_cast<const int*>(tables + t_rorder);
    L.d_diag = reinterpret_cast<const int*>(tables + t_diag);
    L.d_items = reinterpret_cast<const SolveAccum*>(tables + t_items);
    L.d_parts = reinterpret_cast<const SolvePart*>(tables + t_parts);
    L.d_copies = reinterpret_cast<const SolveCopy*>(tables + t_copies);
    return t_end;
}

}  // namespace solve_detail

// The launch and message tables of a level's sweeps.  Each wave's
// accumulations (local neighbors, and per-peer blocks of the outbox) and
// vector refreshes are fixed by the factorization, so they are built once;
// one exchange with the neighbor ranks gives each rank the layout of what
// it will receive.  Returns the bytes of the tables on the device.
template<typename CoordType, typename DataType>
size_t build_solve_tables(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl, int level,
                          DeviceSolveLevel& L, const std::vector<int>& remote_lists,
                          const std::unordered_map<int64_t, size_t>& remote_list_at, double& comm_seconds) {
    const size_t I = sizeof(int);
    const size_t nb = L.boxes.size();
    const size_t np = L.peers.size();
    auto is_local = [&](int64_t m) {
        return m >= lvl.local_morton_start && m < lvl.local_morton_start + static_cast<int64_t>(nb);
    };
    std::unordered_map<int, size_t> peer_index;
    for (size_t i = 0; i < np; ++i) peer_index.emplace(L.peers[i], i);
    auto local_box = [&](int64_t m) -> const SolveBox& { return L.boxes[static_cast<size_t>(m - lvl.local_morton_start)]; };

    // waves (eliminated boxes only), work offsets
    const auto bins = solve_detail::color_waves(lvl, tree->dimension);
    std::vector<int> order, rorder, diag;
    L.waves.assign(bins.size(), DeviceSolveLevel::Wave{});
    for (size_t w = 0; w < bins.size(); ++w) {
        auto& W = L.waves[w];
        W.order0 = static_cast<int>(order.size());
        W.diag0 = static_cast<int>(diag.size());
        int64_t work = 0;
        for (int b : bins[w]) {
            SolveBox& sb = L.boxes[static_cast<size_t>(b)];
            if (sb.k == 0 || sb.r == 0) continue;
            sb.work = work;
            work += sb.ntot;
            W.max_ntot = std::max(W.max_ntot, sb.ntot);
            order.push_back(b);
            diag.push_back(b);   // (its diagonal solve runs with its wave: H2 8x, the LU streamed by wave)
        }
        W.count = static_cast<int>(order.size()) - W.order0;
        W.ndiag = static_cast<int>(diag.size()) - W.diag0;
        rorder.insert(rorder.end(), order.rbegin(), order.rbegin() + W.count);
        L.work_points = std::max(L.work_points, work);
    }
    L.ndiag = static_cast<int>(diag.size());

    // forward: each wave's targets, in first-contribution order
    std::vector<SolveAccum> items;
    std::vector<SolvePart> parts;
    std::vector<std::vector<int64_t>> head(np);  // to each peer: (wave, morton, full, count)... then (wave, morton)...
    for (size_t w = 0; w < L.waves.size(); ++w) {
        auto& W = L.waves[w];
        struct Target { int64_t morton; int count; int full; std::vector<SolvePart> parts; };
        std::vector<Target> targets;
        std::unordered_map<int64_t, size_t> target_of;
        for (int i = W.order0; i < W.order0 + W.count; ++i) {
            const SolveBox& sb = L.boxes[static_cast<size_t>(order[static_cast<size_t>(i)])];
            for (int a = 0; a < sb.nslots; ++a) {
                const size_t s = static_cast<size_t>(sb.slot0 + a);
                const SolveSlot& sl = L.slots[s];
                const int64_t m = L.slot_morton[s];
                auto it = target_of.find(m);
                if (it == target_of.end()) {
                    it = target_of.emplace(m, targets.size()).first;
                    targets.push_back(Target{m, sl.count, sl.full, {}});
                }
                Target& t = targets[it->second];
                if (t.count != sl.count || t.full != sl.full) throw std::runtime_error("device solve: inconsistent updates");
                t.parts.push_back(SolvePart{sb.work, sl.row0, sb.ntot, 0});
            }
        }
        W.out.assign(np, 0);
        std::vector<std::pair<size_t, int64_t>> remote(targets.size(), {0, -1});  // peer, offset in its block
        for (size_t t = 0; t < targets.size(); ++t) {
            if (is_local(targets[t].morton)) continue;
            auto pit = peer_index.find(solve_detail::owner_of(tree, lvl, level, targets[t].morton));
            if (pit == peer_index.end()) throw std::runtime_error("device solve: update for a non-neighbor rank");
            remote[t] = {pit->second, W.out[pit->second]};
            W.out[pit->second] += targets[t].count;
            auto& h = head[pit->second];
            h.insert(h.end(), {static_cast<int64_t>(w), targets[t].morton, targets[t].full, targets[t].count});
        }
        std::vector<int64_t> prefix(np + 1, 0);
        for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + W.out[q];
        L.out_max = std::max(L.out_max, prefix[np]);
        W.item0 = static_cast<int>(items.size());
        for (size_t t = 0; t < targets.size(); ++t) {
            const Target& tg = targets[t];
            SolveAccum it{};
            it.count = tg.count;
            it.nparts = static_cast<int>(tg.parts.size());
            it.part0 = static_cast<int64_t>(parts.size());
            parts.insert(parts.end(), tg.parts.begin(), tg.parts.end());
            if (remote[t].second < 0) {
                const SolveBox& target = local_box(tg.morton);
                it.dst = target.vec;
                it.ld = target.n;
                it.rows = tg.full ? nullptr : target.skel;
                it.outbox = 0;
            } else {
                it.dst = prefix[remote[t].first] + remote[t].second;
                it.ld = tg.count;
                it.rows = nullptr;
                it.outbox = 1;
            }
            W.max_count = std::max(W.max_count, tg.count);
            items.push_back(it);
        }
        W.nitems = static_cast<int>(items.size()) - W.item0;
    }

    // Before each wave of a sweep that reads neighbors, the vectors other
    // ranks read that changed since the last refresh (all of them before the
    // sweep's first wave): waves in reverse (kind 0, the solve's backward
    // sweep) or in order (kind 1, the multiply's forward sweep).
    std::vector<SolveCopy> copies;
    const int nw = static_cast<int>(L.waves.size());
    std::vector<std::vector<int64_t>> head_ref[2] = {std::vector<std::vector<int64_t>>(np),
                                                     std::vector<std::vector<int64_t>>(np)};
    for (int kind = 0; kind < 2; ++kind) {
        for (int step = 0; step < nw; ++step) {
            const int w = kind == 0 ? nw - 1 - step : step;
            auto& R = L.waves[static_cast<size_t>(w)].ref[kind];
            R.send.assign(np, 0);
            std::vector<int> changed;
            if (step == 0) {
                for (size_t b = 0; b < nb; ++b) changed.push_back(static_cast<int>(b));
            } else {
                const auto& P = L.waves[static_cast<size_t>(kind == 0 ? w + 1 : w - 1)];
                changed.assign(order.begin() + P.order0, order.begin() + P.order0 + P.count);
                std::sort(changed.begin(), changed.end());
            }
            std::vector<std::vector<int>> per_peer(np);
            for (int b : changed) {
                for (int q : L.readers[static_cast<size_t>(b)]) per_peer[static_cast<size_t>(q)].push_back(b);
            }
            std::vector<int64_t> prefix(np + 1, 0);
            for (size_t q = 0; q < np; ++q) {
                for (int b : per_peer[q]) R.send[q] += L.boxes[static_cast<size_t>(b)].n;
                prefix[q + 1] = prefix[q] + R.send[q];
            }
            L.send_max = std::max(L.send_max, prefix[np]);
            R.send0 = static_cast<int>(copies.size());
            for (size_t q = 0; q < np; ++q) {
                int64_t at = prefix[q];
                for (int b : per_peer[q]) {
                    const SolveBox& sb = L.boxes[static_cast<size_t>(b)];
                    copies.push_back(SolveCopy{sb.vec, at, sb.n});
                    at += sb.n;
                    R.send_max_n = std::max(R.send_max_n, sb.n);
                    head_ref[kind][q].insert(head_ref[kind][q].end(), {static_cast<int64_t>(w), lvl.local_morton_start + b});
                }
            }
            R.nsend = static_cast<int>(copies.size()) - R.send0;
        }
    }

    // what the neighbor ranks send us, and where it goes
    if (np > 0) {
        std::vector<PeerMessage> out(np), in;
        for (size_t q = 0; q < np; ++q) {
            auto& h = out[q].head;
            h.push_back(static_cast<int64_t>(head[q].size() / 4));
            h.insert(h.end(), head[q].begin(), head[q].end());
            for (int kind = 0; kind < 2; ++kind) {
                h.push_back(static_cast<int64_t>(head_ref[kind][q].size() / 2));
                h.insert(h.end(), head_ref[kind][q].begin(), head_ref[kind][q].end());
            }
        }
        comm_seconds += exchange_peer_messages(tree->comm, L.peers, out, in, 746);
        struct Fwd { int64_t morton; int full, count; };
        std::vector<std::vector<std::vector<Fwd>>> fwd(L.waves.size(), std::vector<std::vector<Fwd>>(np));
        std::vector<std::vector<std::vector<int64_t>>> bwd[2];
        for (auto& v : bwd) v.assign(L.waves.size(), std::vector<std::vector<int64_t>>(np));
        for (size_t q = 0; q < np; ++q) {
            const auto& h = in[q].head;
            size_t at = 0;
            const int64_t nf = h.at(at++);
            for (int64_t e = 0; e < nf; ++e, at += 4) {
                fwd.at(static_cast<size_t>(h.at(at)))[q].push_back(
                    Fwd{h.at(at + 1), static_cast<int>(h.at(at + 2)), static_cast<int>(h.at(at + 3))});
            }
            for (int kind = 0; kind < 2; ++kind) {
                const int64_t nbk = h.at(at++);
                for (int64_t e = 0; e < nbk; ++e, at += 2) {
                    bwd[kind].at(static_cast<size_t>(h.at(at)))[q].push_back(h.at(at + 1));
                }
            }
            if (at != h.size()) throw std::runtime_error("device solve: malformed layout message");
        }
        for (size_t w = 0; w < L.waves.size(); ++w) {
            auto& W = L.waves[w];
            // forward: updates of local boxes, parts in peer order
            W.in.assign(np, 0);
            for (size_t q = 0; q < np; ++q) {
                for (const Fwd& f : fwd[w][q]) W.in[q] += f.count;
            }
            std::vector<int64_t> prefix(np + 1, 0);
            for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + W.in[q];
            L.in_max = std::max(L.in_max, prefix[np]);
            struct Target { int64_t morton; int count; int full; std::vector<SolvePart> parts; };
            std::vector<Target> targets;
            std::unordered_map<int64_t, size_t> target_of;
            for (size_t q = 0; q < np; ++q) {
                int64_t at = prefix[q];
                for (const Fwd& f : fwd[w][q]) {
                    if (!is_local(f.morton)) throw std::runtime_error("device solve: update for a box of another rank");
                    auto it = target_of.find(f.morton);
                    if (it == target_of.end()) {
                        it = target_of.emplace(f.morton, targets.size()).first;
                        targets.push_back(Target{f.morton, f.count, f.full, {}});
                    }
                    targets[it->second].parts.push_back(SolvePart{at, 0, f.count, 1});
                    at += f.count;
                }
            }
            W.in_item0 = static_cast<int>(items.size());
            for (const Target& t : targets) {
                const SolveBox& target = local_box(t.morton);
                if (t.full ? t.count != target.n : t.count != target.k) {
                    throw std::runtime_error("device solve: update of box " + std::to_string(t.morton) + " has another size");
                }
                SolveAccum it{};
                it.dst = target.vec;
                it.ld = target.n;
                it.count = t.count;
                it.rows = t.full ? nullptr : target.skel;
                it.outbox = 0;
                it.nparts = static_cast<int>(t.parts.size());
                it.part0 = static_cast<int64_t>(parts.size());
                parts.insert(parts.end(), t.parts.begin(), t.parts.end());
                W.in_max_count = std::max(W.in_max_count, t.count);
                items.push_back(it);
            }
            W.in_nitems = static_cast<int>(items.size()) - W.in_item0;
            // refreshed vectors into the ghost area, for either order
            for (int kind = 0; kind < 2; ++kind) {
                auto& R = W.ref[kind];
                R.recv.assign(np, 0);
                R.recv0 = static_cast<int>(copies.size());
                for (size_t q = 0; q < np; ++q) {
                    for (int64_t m : bwd[kind][w][q]) R.recv[q] += L.ghosts[static_cast<size_t>(L.ghost_of.at(m))].n;
                }
                prefix.assign(np + 1, 0);
                for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + R.recv[q];
                L.recv_max = std::max(L.recv_max, prefix[np]);
                for (size_t q = 0; q < np; ++q) {
                    int64_t at = prefix[q];
                    for (int64_t m : bwd[kind][w][q]) {
                        const auto& g = L.ghosts[static_cast<size_t>(L.ghost_of.at(m))];
                        copies.push_back(SolveCopy{at, g.offset, g.n});
                        at += g.n;
                        R.recv_max_n = std::max(R.recv_max_n, g.n);
                    }
                }
                R.nrecv = static_cast<int>(copies.size()) - R.recv0;
            }
        }
    }

    solve_detail::build_solve_xnr_waves(L, order, Context::instance().stream(), lvl.level);
    return solve_detail::upload_solve_tables(L, remote_lists, remote_list_at, order, rorder, diag, items, parts, copies);
}

// The tables of a replicated CA level's sweeps: the CA groups (the host's
// make_CA_box_groups) x Morton parities over the local boxes and the ghosts
// (L.boxes: local first, then ghosts; ghost_at: Morton -> box index), with
// group_end marking where each group's waves end; each wave's updates go to
// its local or ghost neighbors (assisting boxes are read only), summed in
// the order of first contribution.  The diagonal solves are the local
// boxes'.  No messages.
template<typename CoordType, typename DataType>
size_t build_ca_solve_tables(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl,
                             DeviceSolveLevel& L, size_t nb, const std::unordered_map<int64_t, size_t>& ghost_at,
                             const std::vector<int>& remote_lists,
                             const std::unordered_map<int64_t, size_t>& remote_list_at) {
    const int num_waves = 1 << tree->dimension;
    auto writable = [&](int64_t m) -> int64_t {
        const int64_t local = m - lvl.local_morton_start;
        if (local >= 0 && local < static_cast<int64_t>(nb)) return local;
        const auto it = ghost_at.find(m);
        return it == ghost_at.end() ? -1 : static_cast<int64_t>(it->second);
    };
    std::vector<int> order, rorder, diag;
    L.waves.clear();
    L.group_end.clear();
    // the CA groups (as butterfly::make_CA_box_groups and the factorization's
    // factorize_CA_level_gpu): blue/orange, purple (3D), green, interior; on
    // one rank boundary, interior
    std::vector<std::vector<int64_t>> groups;
    {
        std::vector<int64_t> interior;
        for (int64_t local_idx : lvl.interior_id) interior.push_back(local_idx + lvl.local_morton_start);
        if (lvl.num_active_processes == 1) {
            std::vector<int64_t> boundary;
            for (int64_t local_idx : lvl.boundary_id) boundary.push_back(local_idx + lvl.local_morton_start);
            groups.push_back(std::move(boundary));
        } else {
            std::vector<int64_t> blue_orange = lvl.blue;
            blue_orange.insert(blue_orange.end(), lvl.orange.begin(), lvl.orange.end());
            groups.push_back(std::move(blue_orange));
            if (tree->dimension == 3) groups.push_back(lvl.purple);
            groups.push_back(lvl.green);
        }
        groups.push_back(std::move(interior));
    }
    for (const auto& group : groups) {
        std::vector<std::vector<int>> bins(static_cast<size_t>(num_waves));
        for (int64_t m : group) {
            const int64_t b = writable(m);
            if (b < 0) throw std::runtime_error("device solve: CA box " + std::to_string(m) + " neither local nor a ghost");
            bins[static_cast<size_t>(m & (num_waves - 1))].push_back(static_cast<int>(b));
        }
        for (const auto& bin : bins) {
            DeviceSolveLevel::Wave W;
            W.order0 = static_cast<int>(order.size());
            W.diag0 = static_cast<int>(diag.size());
            int64_t work = 0;
            for (int b : bin) {
                SolveBox& sb = L.boxes[static_cast<size_t>(b)];
                if (sb.k == 0 || sb.r == 0) continue;
                sb.work = work;
                work += sb.ntot;
                W.max_ntot = std::max(W.max_ntot, sb.ntot);
                order.push_back(b);
                if (b < static_cast<int>(nb)) diag.push_back(b);   // (the local boxes' diagonal solves, with their wave)
            }
            W.count = static_cast<int>(order.size()) - W.order0;
            W.ndiag = static_cast<int>(diag.size()) - W.diag0;
            rorder.insert(rorder.end(), order.rbegin(), order.rbegin() + W.count);
            L.work_points = std::max(L.work_points, work);
            L.waves.push_back(std::move(W));
        }
        L.group_end.push_back(static_cast<int>(L.waves.size()));
    }
    L.ndiag = static_cast<int>(diag.size());

    std::vector<SolveAccum> items;
    std::vector<SolvePart> parts;
    for (auto& W : L.waves) {
        struct Target { int64_t morton; int count; int full; std::vector<SolvePart> parts; };
        std::vector<Target> targets;
        std::unordered_map<int64_t, size_t> target_of;
        for (int i = W.order0; i < W.order0 + W.count; ++i) {
            const SolveBox& sb = L.boxes[static_cast<size_t>(order[static_cast<size_t>(i)])];
            for (int a = 0; a < sb.nslots; ++a) {
                const size_t sidx = static_cast<size_t>(sb.slot0 + a);
                const SolveSlot& sl = L.slots[sidx];
                if (sl.ghost) continue;  // an assisting box: its owner updates it
                const int64_t m = L.slot_morton[sidx];
                auto it = target_of.find(m);
                if (it == target_of.end()) {
                    it = target_of.emplace(m, targets.size()).first;
                    targets.push_back(Target{m, sl.count, sl.full, {}});
                }
                Target& t = targets[it->second];
                if (t.count != sl.count || t.full != sl.full) throw std::runtime_error("device solve: inconsistent updates");
                t.parts.push_back(SolvePart{sb.work, sl.row0, sb.ntot, 0});
            }
        }
        W.item0 = static_cast<int>(items.size());
        for (const Target& tg : targets) {
            const SolveBox& target = L.boxes[static_cast<size_t>(writable(tg.morton))];
            SolveAccum it{};
            it.dst = target.vec;
            it.ld = target.n;
            it.count = tg.count;
            it.rows = tg.full ? nullptr : target.skel;
            it.outbox = 0;
            it.nparts = static_cast<int>(tg.parts.size());
            it.part0 = static_cast<int64_t>(parts.size());
            parts.insert(parts.end(), tg.parts.begin(), tg.parts.end());
            W.max_count = std::max(W.max_count, tg.count);
            items.push_back(it);
        }
        W.nitems = static_cast<int>(items.size()) - W.item0;
    }
    solve_detail::build_solve_xnr_waves(L, order, Context::instance().stream(), lvl.level);
    return solve_detail::upload_solve_tables(L, remote_lists, remote_list_at, order, rorder, diag, items, parts, {});
}

// An upper bound of the bytes a level's kept solve factors take (T, X_SR,
// X_NR, the LU of X_RR, pivots; kept_solve_factors): per box of n points with
// r redundant, k = n - r skeleton ones and neighbors of ntot points,
// 2 k r + r^2 = r (2 n - r) <= n^2 and ntot r <= ntot n.  The boxes the level
// eliminates: the local ones, and a replicated CA level's ghosts.  A
// neighbor of unknown size (another rank's) counts as the box's own size.
// The factors usually take about half of it (r about n / 2, and neighbors
// eliminated before the box have shrunk to their skeletons): the keep
// decision scales it (factorization_driver.hpp decide_solve_keep).
template<typename CoordType, typename DataType>
double solve_keep_bound(const TreeLevel<CoordType, DataType>& lvl, bool ca, double xnr_divisor = 1.0,
                        bool nonsym = false) {
    const double D = static_cast<double>(sizeof(DataType)), I = sizeof(int);
    std::unordered_map<int64_t, int64_t> size;
    for (const auto& box : lvl.local_boxes) size.emplace(box.morton_index, box.num_points);
    if (ca) {
        for (const auto& box : lvl.ghost_boxes) size.emplace(box.morton_index, box.num_points);
    }
    double bytes = 0.0;
    auto add = [&](const BoxData<CoordType, DataType>& box) {
        const double n = static_cast<double>(box.num_points);
        if (n <= 0.0) return;
        double ntot = 0.0;
        for (int64_t m : box.one_hop) {
            const auto it = size.find(m);
            ntot += it != size.end() ? static_cast<double>(it->second) : n;
        }
        // (a nonsymmetric matrix: X_RS and X_RN as well, as large as X_SR and X_NR)
        bytes += (nonsym ? 2.0 : 1.0) * (n * n + ntot * n / xnr_divisor) * D + n * I + 7.0 * 256.0;
    };
    for (const auto& box : lvl.local_boxes) add(box);
    if (ca) {
        for (const auto& box : lvl.ghost_boxes) add(box);
    }
    return bytes;
}

// Kept boxes whose factors have no host copy (the eliminator skips the
// downloads of X_NR, X_SR and the LU of X_RR with its pivots for the boxes it
// keeps): copied from the kept factors, before a host solve or multiply
// reads them, or the kept factors are released.
template<typename CoordType, typename DataType>
void materialize_host_factors_level(ParallelTree<CoordType, DataType>* tree, int level, const KeptSolveLevel& kl) {
    const auto& kl_boxes = kl.streamed ? kl.boxes_staged : kl.boxes;   // (a streamed level: its staged-form pointers)
    if (level < 0 || level >= tree->num_levels) return;
    Context::instance().activate();
    check_cuda(cudaStreamSynchronize(Context::instance().stream()), "host factors");
    auto& lvl = tree->levels[static_cast<size_t>(level)];
    auto fill = [&](BoxData<CoordType, DataType>& box) {
        const auto it = kl_boxes.find(box.morton_index);
        if (it == kl_boxes.end()) return;
        const KeptSolveBox& kb = it->second;
        const size_t k = static_cast<size_t>(kb.k), r = static_cast<size_t>(kb.r);
        auto copy = [&](const void* src, size_t count) {
            std::vector<DataType> v(count);
            kept_fetch(kl, src, count * sizeof(DataType), v.data(), "host factors");
            return v;
        };
        if (kb.xnr_pending) {
            throw std::runtime_error("H2 GPU: the kept X_NR of box " + std::to_string(box.morton_index) +
                                     " was never transposed out of its fill source");
        }
        if (kb.ntot > 0 && !box.X_NR.is_allocated()) {
            // H2_compress_factors: kb.xnr is xnr_bytes of compressed data, not
            // ntot * r elements.  Bring the compressed bytes over (far fewer)
            // and decompress on the host, which uses the same format.
            if (kb.xnr_bytes > 0) {
#ifdef H2_HAVE_NVCOMP
                std::vector<char> packed(kb.xnr_bytes);
                kept_fetch(kl, kb.xnr, kb.xnr_bytes, packed.data(), "host factors (compressed X_NR)");
                std::vector<DataType> v(static_cast<size_t>(kb.ntot) * r);
                decompress_xnr_host<typename DeviceScalar<DataType>::type>(
                    packed.data(), reinterpret_cast<typename DeviceScalar<DataType>::type*>(v.data()),
                    kl.chunk_strides.data() + kb.xnr_chunk0, kb.xnr_nchunks);
                box.X_NR.set_owned(kb.ntot, kb.r, std::move(v), MatrixStorage<DataType>::FULL);
#else
                throw std::runtime_error(
                    "H2 GPU: compressed X_NR on the host needs a build with nvCOMP");
#endif
            } else
            box.X_NR.set_owned(kb.ntot, kb.r, copy(kb.xnr, static_cast<size_t>(kb.ntot) * r),
                               MatrixStorage<DataType>::FULL);
        }
        if (!box.X_SR.is_allocated()) {
            if (kb.xsr_bytes > 0) {  // H2_compress_factors: compressed like X_NR
#ifdef H2_HAVE_NVCOMP
                std::vector<char> packed(kb.xsr_bytes);
                kept_fetch(kl, kb.xsr, kb.xsr_bytes, packed.data(), "host factors (compressed X_SR)");
                std::vector<DataType> v(k * r);
                decompress_xnr_host<typename DeviceScalar<DataType>::type>(
                    packed.data(), reinterpret_cast<typename DeviceScalar<DataType>::type*>(v.data()),
                    kl.chunk_strides.data() + kb.xsr_chunk0, kb.xsr_nchunks);
                box.X_SR.set_owned(kb.k, kb.r, std::move(v), MatrixStorage<DataType>::FULL);
#else
                throw std::runtime_error(
                    "H2 GPU: compressed X_SR on the host needs a build with nvCOMP");
#endif
            } else
            box.X_SR.set_owned(kb.k, kb.r, copy(kb.xsr, k * r), MatrixStorage<DataType>::FULL);
        }
        // (a nonsymmetric matrix: X_RS (r x k) and X_RN (r x ntot), from their
        // transposed kept copies)
        auto transposed = [&](const void* src, size_t rows, size_t cols) {  // src rows x cols -> cols x rows
            const std::vector<DataType> t = copy(src, rows * cols);
            std::vector<DataType> v(rows * cols);
            for (size_t j = 0; j < cols; ++j)
                for (size_t i = 0; i < rows; ++i) v[j + i * cols] = t[i + j * rows];
            return v;
        };
        if (kb.xrs_t != nullptr && !box.X_RS.is_allocated()) {
            box.X_RS.set_owned(kb.r, kb.k, transposed(kb.xrs_t, k, r), MatrixStorage<DataType>::FULL);
        }
        if (kb.xrn_t != nullptr && kb.ntot > 0 && !box.X_RN.is_allocated()) {
            box.X_RN.set_owned(kb.r, kb.ntot, transposed(kb.xrn_t, static_cast<size_t>(kb.ntot), r),
                               MatrixStorage<DataType>::FULL);
        }
        if (!box.X_RR.is_allocated()) {
            box.X_RR.set_owned(kb.r, kb.r, copy(kb.lu, r * r), MatrixStorage<DataType>::LU_FACTORED);
            box.X_RR_pivots.resize(r);
            kept_fetch(kl, kb.ipiv, r * sizeof(int), box.X_RR_pivots.data(), "host factors");
        }
    };
    for (auto& box : lvl.local_boxes) fill(box);
    for (auto& box : lvl.ghost_boxes) fill(box);
}

template<typename CoordType, typename DataType>
void materialize_host_factors(ParallelTree<CoordType, DataType>* tree, DeviceSolveStore& store) {
    for (auto& [level, kl] : store.kept) materialize_host_factors_level(tree, level, kl);
}
template<typename CoordType, typename DataType>
void materialize_host_factors(ParallelTree<CoordType, DataType>* tree) {
    materialize_host_factors(tree, device_solve_store());
}

// The log-determinant contributions (the sum of log|u_kk| and of the
// arguments of the diagonal of the LU of X_RR, LAPACK pivots counted, as
// core/butterfly_logdet.hpp accumulate_logdet_lu) of the local boxes of
// `level` whose LU stayed on the device without a host copy, from the kept
// factors.  Returns how many boxes it covered (0 when the level's factors
// are not kept); the host sums the others.
template<typename CoordType, typename DataType>
int64_t device_logdet_kept_level(ParallelTree<CoordType, DataType>* tree, int level, double& logabs, double& arg) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)tree; (void)level; (void)logabs; (void)arg;
    return 0;
  } else {
    using S = typename DeviceScalar<DataType>::type;
    DeviceSolveStore& store = device_solve_store();
    const auto kit = store.kept.find(level);
    if (kit == store.kept.end() || level < 0 || level >= tree->num_levels) return 0;
    KeptSolveLevel& kl = kit->second;
    const auto& lvl = tree->levels[static_cast<size_t>(level)];
    std::vector<LogdetItem> items;
    std::vector<int64_t> mortons;
    for (const auto& box : lvl.local_boxes) {
        const size_t r = box.redundant_indices.size();
        if (r == 0 || box.X_RR.is_allocated()) continue;
        const auto it = kl.boxes.find(box.morton_index);
        if (it == kl.boxes.end() || it->second.r != static_cast<int>(r)) {
            throw std::runtime_error("logdet: box " + std::to_string(box.morton_index) + " of level " +
                                     std::to_string(level) + " has neither a host LU nor a kept one");
        }
        items.push_back(LogdetItem{it->second.lu, it->second.ipiv, static_cast<int>(r)});
        mortons.push_back(box.morton_index);
    }
    if (items.empty()) return 0;
    if (kl.on_host()) {
        // the sum was taken before the LU left the device (stage_kept_level
        // runs this first); the boxes it covered are those without a host
        // LU, which is what `items` holds now as well
        if (kl.logdet_boxes != static_cast<int64_t>(items.size())) {
            throw std::runtime_error("logdet: level " + std::to_string(level) + " was staged with " +
                                     std::to_string(kl.logdet_boxes) + " kept LUs, " +
                                     std::to_string(items.size()) + " are asked for");
        }
        for (int64_t i = 0; i < kl.logdet_boxes; ++i) {  // box by box, as the device path adds them
            logabs += kl.logdet_per_box[static_cast<size_t>(2 * i)];
            arg += kl.logdet_per_box[static_cast<size_t>(2 * i + 1)];
        }
        return kl.logdet_boxes;
    }
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    const size_t n = items.size();
    const size_t item_bytes = align_up(n * sizeof(LogdetItem)), out_bytes = align_up(2 * n * sizeof(double));
    DeviceBuffer buffer;
    char* d = static_cast<char*>(buffer.reserve(item_bytes + out_bytes + align_up(n * sizeof(int))));
    check_cuda(cudaMemcpyAsync(d, items.data(), n * sizeof(LogdetItem), cudaMemcpyHostToDevice, stream), "logdet");
    launch_logdet_lu<S>(reinterpret_cast<const LogdetItem*>(d), static_cast<int>(n),
                        reinterpret_cast<double*>(d + item_bytes), reinterpret_cast<int*>(d + item_bytes + out_bytes),
                        stream);
    std::vector<double> out(2 * n);
    std::vector<int> bad(n);
    check_cuda(cudaMemcpyAsync(out.data(), d + item_bytes, 2 * n * sizeof(double), cudaMemcpyDeviceToHost, stream),
               "logdet");
    check_cuda(cudaMemcpyAsync(bad.data(), d + item_bytes + out_bytes, n * sizeof(int), cudaMemcpyDeviceToHost, stream),
               "logdet");
    check_cuda(cudaStreamSynchronize(stream), "logdet");
    for (size_t i = 0; i < n; ++i) {
        if (bad[i] != 0) {
            throw std::runtime_error("logdet: singular/near-singular or invalid pivot at k = " +
                                     std::to_string(bad[i] - 1) + " of box " + std::to_string(mortons[i]) +
                                     " (level " + std::to_string(level) + ", on the device)");
        }
        logabs += out[2 * i];
        arg += out[2 * i + 1];
    }
    kl.logdet_per_box = out;  // (kept for after the level is staged)
    kl.logdet_boxes = static_cast<int64_t>(n);
    return static_cast<int64_t>(n);
  }
}

// A replicated CA level's halo exchanges in device memory, in place of the
// host gather_CA_boxes_solve: the owners' vectors of the ghosts (into the
// level vector) and of the assisting boxes (into the read-only area), packed
// from and unpacked into device memory by copies.  Each rank sends the
// owners its requests (Morton, ghost or assisting, points) once; a message
// holds the requested boxes in the requester's order.  Collective over the
// tree's ranks (also those inactive on the level).  Returns the bytes of the
// copy table on the device.
template<typename CoordType, typename DataType>
size_t build_ca_halo(ParallelTree<CoordType, DataType>* tree, int level, DeviceSolveLevel& L, double& comm_seconds) {
    const auto t0 = std::chrono::steady_clock::now();
    auto& lvl = tree->levels[static_cast<size_t>(level)];
    int P = 0, me = 0;
    MPI_Comm_size(tree->comm, &P);
    MPI_Comm_rank(tree->comm, &me);
    const bool on = L.active && !L.host;
    const size_t nb = L.boxes.size() - L.ghost_box.size();
    struct Want { int64_t m; int n; bool ghost; int64_t dst; };
    std::vector<Want> want;
    if (on) {
        for (size_t j = 0; j < L.ghost_box.size(); ++j) {
            want.push_back({L.ghost_box[j], L.boxes[nb + j].n, true, L.boxes[nb + j].vec});
        }
        for (const auto& g : L.ghosts) want.push_back({g.morton, g.n, false, g.offset});
    }
    std::sort(want.begin(), want.end(), [](const Want& a, const Want& b) { return a.m < b.m; });
    std::vector<std::vector<size_t>> by_owner(static_cast<size_t>(P));
    if (!want.empty()) {
        std::vector<uint64_t> ms;
        ms.reserve(want.size());
        for (const Want& w : want) ms.push_back(static_cast<uint64_t>(w.m));
        const std::vector<uint32_t> region =
            morton::assign_to_processes_nd(tree->dimension, ms, lvl.num_active_processes, 1u << level);
        for (size_t i = 0; i < want.size(); ++i) {
            const int owner = lvl.morton_to_rank.at(static_cast<int>(region[i]));
            if (owner == me || owner < 0 || owner >= P) throw std::runtime_error("device solve: CA halo box owner");
            by_owner[static_cast<size_t>(owner)].push_back(i);
        }
    }
    // requests: (2 m + ghost, n) per box
    std::vector<int> scount(static_cast<size_t>(P)), rcount(static_cast<size_t>(P)), sdispl(static_cast<size_t>(P)),
        rdispl(static_cast<size_t>(P));
    for (int q = 0; q < P; ++q) scount[static_cast<size_t>(q)] = 2 * static_cast<int>(by_owner[static_cast<size_t>(q)].size());
    MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, tree->comm);
    for (int q = 1; q < P; ++q) {
        sdispl[static_cast<size_t>(q)] = sdispl[static_cast<size_t>(q - 1)] + scount[static_cast<size_t>(q - 1)];
        rdispl[static_cast<size_t>(q)] = rdispl[static_cast<size_t>(q - 1)] + rcount[static_cast<size_t>(q - 1)];
    }
    std::vector<int64_t> sbuf(static_cast<size_t>(sdispl.back() + scount.back()));
    std::vector<int64_t> rbuf(static_cast<size_t>(std::max(rdispl.back() + rcount.back(), 1)));
    for (int q = 0; q < P; ++q) {
        size_t at = static_cast<size_t>(sdispl[static_cast<size_t>(q)]);
        for (size_t i : by_owner[static_cast<size_t>(q)]) {
            sbuf[at++] = 2 * want[i].m + (want[i].ghost ? 1 : 0);
            sbuf[at++] = want[i].n;
        }
    }
    MPI_Alltoallv(sbuf.data(), scount.data(), sdispl.data(), MPI_INT64_T, rbuf.data(), rcount.data(), rdispl.data(),
                  MPI_INT64_T, tree->comm);
    comm_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!on) return 0;

    L.peers.clear();
    for (int q = 0; q < P; ++q) {
        if (scount[static_cast<size_t>(q)] > 0 || rcount[static_cast<size_t>(q)] > 0) L.peers.push_back(q);
    }
    const size_t np = L.peers.size();
    std::vector<SolveCopy> copies;
    L.send_max = L.recv_max = 0;
    for (int kind = 0; kind < 3; ++kind) {
        auto wanted = [&](bool ghost) { return kind == 1 || ghost == (kind == 0); };
        DeviceSolveLevel::CaHalo& H = L.halo[kind];
        H.out.assign(np, 0);
        H.in.assign(np, 0);
        // sends: the requests of each peer, in its order
        H.send0 = static_cast<int>(copies.size());
        int64_t at = 0;
        for (size_t p = 0; p < np; ++p) {
            const int q = L.peers[p];
            for (int e = rdispl[static_cast<size_t>(q)]; e < rdispl[static_cast<size_t>(q)] + rcount[static_cast<size_t>(q)];
                 e += 2) {
                const int64_t code = rbuf[static_cast<size_t>(e)], n = rbuf[static_cast<size_t>(e + 1)];
                if (!wanted((code & 1) != 0)) continue;
                const int64_t b = (code >> 1) - lvl.local_morton_start;
                if (b < 0 || b >= static_cast<int64_t>(nb) || L.boxes[static_cast<size_t>(b)].n != n) {
                    throw std::runtime_error("device solve: CA halo request for box " + std::to_string(code >> 1));
                }
                copies.push_back(SolveCopy{L.boxes[static_cast<size_t>(b)].vec, at, static_cast<int>(n)});
                H.send_max_n = std::max(H.send_max_n, static_cast<int>(n));
                H.out[p] += n;
                at += n;
            }
        }
        H.nsend = static_cast<int>(copies.size()) - H.send0;
        H.out_total = at;
        // receives: the ghosts into the level vector, the assisting boxes into the area
        std::vector<SolveCopy> to_vec, to_area;
        at = 0;
        for (size_t p = 0; p < np; ++p) {
            for (size_t i : by_owner[static_cast<size_t>(L.peers[p])]) {
                const Want& w = want[i];
                if (!wanted(w.ghost)) continue;
                (w.ghost ? to_vec : to_area).push_back(SolveCopy{at, w.dst, w.n});
                (w.ghost ? H.vec_max_n : H.area_max_n) = std::max(w.ghost ? H.vec_max_n : H.area_max_n, w.n);
                H.in[p] += w.n;
                at += w.n;
            }
        }
        H.in_total = at;
        H.vec0 = static_cast<int>(copies.size());
        H.nvec = static_cast<int>(to_vec.size());
        copies.insert(copies.end(), to_vec.begin(), to_vec.end());
        H.area0 = static_cast<int>(copies.size());
        H.narea = static_cast<int>(to_area.size());
        copies.insert(copies.end(), to_area.begin(), to_area.end());
        L.send_max = std::max(L.send_max, H.out_total);
        L.recv_max = std::max(L.recv_max, H.in_total);
    }
    const size_t bytes = align_up(std::max<size_t>(copies.size() * sizeof(SolveCopy), 1));
    char* table = DeviceHeap::instance().alloc_resident(bytes);
    L.blocks.push_back(table);
    if (!copies.empty()) {
        check_cuda(cudaMemcpy(table, copies.data(), copies.size() * sizeof(SolveCopy), cudaMemcpyHostToDevice),
                   "CA halo table");
    }
    L.d_copies = reinterpret_cast<const SolveCopy*>(table);
    L.device_halo = true;
    return bytes;
}

template<typename CoordType, typename DataType>
bool prepare_device_solve(ParallelTree<CoordType, DataType>* tree, int verbosity, bool host_levels_allowed) {
    DeviceSolveStore& store = device_solve_store();
    if (store.decided && store.tree == tree) return store.usable;
    store.release_prepared();
    store.tree = tree;
    store.decided = true;
    {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        const size_t D = sizeof(DataType), I = sizeof(int);
        const int leaf = tree->num_levels - 1;
        // a nonsymmetric matrix: X_RS and X_RN too (the backward sweep's)
        const bool nonsym = !tree->factor_symmetric;
        int rank = 0;
        MPI_Comm_rank(tree->comm, &rank);

        // (a component-owner level with its factors at home solves as a Color level)
        auto replicated = [&](int level) {
            return tree->level_uses_CA(level) && !tree->levels[static_cast<size_t>(level)].ca_home_factors;
        };
        // ---- eligibility and size, on every rank
        std::string reason;
        std::vector<double> level_need(static_cast<size_t>(leaf + 1), 0.0);
        std::vector<double> level_staged(static_cast<size_t>(leaf + 1), 0.0);  // ... of which staged kept factors (H2 8x)
        std::vector<int> host_only(static_cast<size_t>(leaf + 1), 0);  // CA levels factored on the host
        for (int level = 2; level <= leaf && reason.empty(); ++level) {
            auto& lvl = tree->levels[static_cast<size_t>(level)];
            if (!lvl.is_process_active) continue;
            double& need = level_need[static_cast<size_t>(level)];
            const bool ca = tree->level_uses_CA(level) && !lvl.ca_home_factors;
            if (ca && store.ca_levels.count(level) == 0) {
                if (host_levels_allowed) {
                    host_only[static_cast<size_t>(level)] = 1;
                    continue;
                }
                reason = "CA levels not factored by the device's replicated CA";
                break;
            }
            const bool kept_level = store.kept.count(level) != 0;
            if (kept_level && store.kept.at(level).on_host()) {  // (uploaded again at the setup, or streamed)
                need += static_cast<double>(store.kept.at(level).staged_bytes);
                level_staged[static_cast<size_t>(level)] = static_cast<double>(store.kept.at(level).staged_bytes);
            }
            auto check = [&](const BoxData<CoordType, DataType>& box, bool kept_box) {
                const int64_t n = box.num_points, k = static_cast<int64_t>(box.skeleton_indices.size()),
                              r = static_cast<int64_t>(box.redundant_indices.size());
                if (static_cast<size_t>(n) * D > kSolveSharedLimit) {
                    reason = "boxes too large";
                    return;
                }
                // (a kept box's LU may have no host copy; the kept one is an LU)
                if (r > 0 && k > 0 && !(kept_box && !box.X_RR.is_allocated()) &&
                    box.X_RR.format != MatrixStorage<DataType>::LU_FACTORED) {
                    reason = "X_RR not LU-factored (use H2_XRR_factor=1)";
                    return;
                }
                const int64_t ntot = box.X_NR.is_allocated() ? box.X_NR.rows : 0;
                const int64_t factors = 2 * k * r + ntot * r + r * r + (nonsym ? k * r + ntot * r : 0);
                need += static_cast<double>((kept_box ? 0 : factors * D) + (r + k + r) * I + 1024);
            };
            for (const auto& box : lvl.local_boxes) {
                if (reason.empty()) check(box, kept_level);
            }
            if (ca) {
                for (const auto& box : lvl.ghost_boxes) {  // (eliminated here again: their factors too)
                    if (reason.empty()) check(box, kept_level);
                }
            }
        }
        DeviceHeap& heap = DeviceHeap::instance();
        heap.ensure_initialized();
        // Residency, the same on every rank (each level's largest need, the
        // smallest room): every Color level on the device, then the CA levels
        // while they fit, leaf first (the most work); a CA level that does
        // not fit runs its sweeps on the host (run_device_solve's host_level).
        double room = static_cast<double>(heap.capacity() - heap.used());
        int bad = reason.empty() ? 0 : 1;
        MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, tree->comm);
        std::vector<double> need_max = level_need;
        MPI_Allreduce(MPI_IN_PLACE, need_max.data(), leaf + 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, host_only.data(), leaf + 1, MPI_INT, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &room, 1, MPI_DOUBLE, MPI_MIN, tree->comm);
        std::vector<double> staged_max = level_staged;
        MPI_Allreduce(MPI_IN_PLACE, staged_max.data(), leaf + 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        const double budget = room - static_cast<double>(size_t{256} << 20);
        // H2 8x: the X_NR scratch every level's sweeps share - kXnrBuffers
        // buffers of one segment, a segment at most 2 GB (an eighth of the
        // room when that is less), so a level's widest wave no longer sets it.
        // Only for compressed kept factors (H2_compress_factors; the others
        // are read in place), and the process's scratch when it has one
        // (shared_xnr_scratch) and its buffers hold this operator's widest
        // box - else one of its own: the plan counts a new one only, the
        // largest over the ranks.  A box's width is what build_solve_xnr_waves
        // puts in a buffer for it, its compressed factors uncompressed (here
        // an upper bound: T and X_SR are counted where a level reads neither).
        int compressed = 0;
        size_t widest_box = 0;
        for (const auto& [level, kl] : store.kept) {
            for (const auto& [morton, kb] : kl.boxes) {
                if (kb.xnr_bytes == 0 && kb.t_bytes == 0 && kb.xsr_bytes == 0) continue;
                compressed = 1;
                const size_t kr = static_cast<size_t>(kb.k) * static_cast<size_t>(kb.r) * D;
                size_t raw = kb.xnr_bytes > 0
                    ? align_up(static_cast<size_t>(kb.ntot) * static_cast<size_t>(kb.r) * D) : 0;
                if (kb.t_bytes > 0) raw += align_up(kr);
                if (kb.xsr_bytes > 0) raw += align_up(kr);
                widest_box = std::max(widest_box, raw);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &compressed, 1, MPI_INT, MPI_MAX, tree->comm);
        const size_t segment = std::min<size_t>(size_t{2} << 30,
                                                std::max<size_t>(size_t{256} << 20, static_cast<size_t>(budget / 8)));
        const SharedXnrScratch& shared_now = shared_xnr_scratch();
        const bool own_scratch = compressed && shared_now.ptr != nullptr && widest_box > shared_now.segment;
        double scratch = compressed && (shared_now.ptr == nullptr || own_scratch)
                             ? static_cast<double>(kXnrBuffers * segment) : 0.0;
        MPI_Allreduce(MPI_IN_PLACE, &scratch, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        double planned = scratch;
        for (int level = 2; level <= leaf; ++level) {
            if (!replicated(level)) planned += 1.05 * need_max[static_cast<size_t>(level)];
        }
        // H2 8x (2026-10-11): the levels whose staged kept factors do not fit
        // beside the rest are streamed - a wave's blocks come into a ring of
        // slots ahead of its step instead of the level staying resident -
        // largest first until the rest fits.  The same on every rank
        // (allreduced sizes).  At 384^3 on 8 ranks the kept factors are 27 GB
        // a rank against a 24 GB budget; the ring is a few GB.
        std::vector<double> home_need(static_cast<size_t>(leaf + 1), 0.0);
        for (int level = 2; level <= leaf; ++level) {
            auto kit = store.kept.find(level);
            if (kit == store.kept.end() || !kit->second.on_host() || kit->second.staged_bytes == 0) continue;
            if (host_only[static_cast<size_t>(level)]) continue;
            auto& lvl = tree->levels[static_cast<size_t>(level)];
            if (!lvl.is_process_active) continue;
            plan_kept_streaming(kit->second, solve_wave_mortons(tree, lvl, replicated(level)));
            home_need[static_cast<size_t>(level)] = static_cast<double>(kit->second.ring_bytes + kStagedAlign);
        }
        MPI_Allreduce(MPI_IN_PLACE, home_need.data(), leaf + 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        std::set<int> streamed;
        double stream_home = 0.0;
        if (!bad && planned > budget) {
            std::vector<std::pair<double, int>> cand;
            for (int level = 2; level <= leaf; ++level) {
                if (!replicated(level) && staged_max[static_cast<size_t>(level)] > 0.0) {
                    cand.emplace_back(staged_max[static_cast<size_t>(level)], level);
                }
            }
            std::sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
            for (const auto& [bytes, level] : cand) {
                if (planned + stream_home <= budget) break;
                if (home_need[static_cast<size_t>(level)] <= 0.0) continue;
                streamed.insert(level);
                planned -= 1.05 * bytes;
                stream_home = std::max(stream_home, home_need[static_cast<size_t>(level)]);
            }
        }
        std::set<int> host_levels;
        if (!bad && planned + stream_home > budget) {
            bad = 1;
            if (reason.empty()) reason = "not enough device memory";
        }
        for (int level = leaf; level >= 2 && !bad; --level) {
            if (!replicated(level)) continue;
            const double n = 1.05 * need_max[static_cast<size_t>(level)];
            const double st = 1.05 * staged_max[static_cast<size_t>(level)];
            const double home = std::max(stream_home, home_need[static_cast<size_t>(level)]);
            if (host_only[static_cast<size_t>(level)]) {
                host_levels.insert(level);
            } else if (planned + n + stream_home <= budget) {
                planned += n;
            } else if (st > 0.0 && home_need[static_cast<size_t>(level)] > 0.0 && planned + (n - st) + home <= budget) {
                // H2 8x: a replicated level streams its staged factors like a
                // Color one before it goes to the host (whose fallback
                // uncompresses the X_NR there)
                streamed.insert(level);
                planned += n - st;
                stream_home = home;
            } else if (host_levels_allowed) {
                host_levels.insert(level);
            } else {
                bad = 1;
                if (reason.empty()) reason = "not enough device memory";
            }
        }
        const int all_ok = bad ? 0 : 1;
        if (rank == 0) {   // the residency decision, with its numbers (H2 8x, 2026-10-10)
            std::string needs;
            for (int level = 2; level <= leaf; ++level) {
                char b[64];
                std::snprintf(b, sizeof b, "%s%d: %.2f", needs.empty() ? "" : ", ", level, need_max[static_cast<size_t>(level)] / 1e9);
                needs += b;
            }
            std::string streams;
            for (int level : streamed) streams += (streams.empty() ? "" : ",") + std::to_string(level);
            {
                const DeviceHeap& h = DeviceHeap::instance();
                std::string held;
                for (size_t k = 0; k < static_cast<size_t>(ResidentKind::Count); ++k) {
                    const size_t bytes = h.resident_now(static_cast<ResidentKind>(k));
                    if (bytes < (size_t{64} << 20)) continue;
                    char b[96];
                    std::snprintf(b, sizeof b, "%s%s %.2f", held.empty() ? "" : ", ", resident_kind_name(static_cast<ResidentKind>(k)),
                                  static_cast<double>(bytes) / 1e9);
                    held += b;
                }
                std::printf("GPU solve: the card at the set-up (rank 0): %.2f GB in use of %.2f; by kind (GB): %s\n",
                            static_cast<double>(h.used()) / 1e9, static_cast<double>(h.capacity()) / 1e9, held.c_str());
            }
            std::printf("GPU solve: residency: planned %.2f GB (X_NR scratch %.2f) against a budget of %.2f (room %.2f GB on "
                        "the tightest rank); needs by level (GB): %s%s%s%s%s%s\n", planned / 1e9, scratch / 1e9, budget / 1e9,
                        room / 1e9, needs.c_str(),
                        streams.empty() ? "" : "; streamed through a ", streams.empty() ? "" : (std::to_string(stream_home / 1e9).substr(0, 5) + " GB home: levels ").c_str(),
                        streams.c_str(), all_ok ? "" : "; off: ", all_ok ? "" : (reason.empty() ? "another rank" : reason.c_str()));
            std::fflush(stdout);
        }
        if (!all_ok) {
            materialize_host_factors(tree);  // (the host solve reads them)
            store.release();  // the kept factors too
            store.tree = tree;
            store.decided = true;
            if (rank == 0) {  // once per factorization, also for a silent first solve
                std::printf("GPU solve: off (%s on some rank); host solve\n", reason.empty() ? "another rank" : reason.c_str());
                std::fflush(stdout);
            }
            return false;
        }

        Context& ctx = Context::instance();
        ctx.activate();
        cudaStream_t stream = ctx.stream();
        store.levels.resize(static_cast<size_t>(leaf + 1));
        if (compressed && own_scratch) {   // (the process's scratch cannot hold its widest box)
            store.xnr_segment_bytes = segment;
            store.xnr_scratch_bytes = kXnrBuffers * segment;
            store.xnr_scratch = DeviceHeap::instance().alloc_resident(store.xnr_scratch_bytes);
            store.xnr_scratch_own = true;
        } else if (compressed) {
            SharedXnrScratch& shared = shared_xnr_scratch();
            if (shared.ptr == nullptr) {
                shared.segment = segment;
                shared.bytes = kXnrBuffers * segment;
                shared.ptr = DeviceHeap::instance().alloc_resident(shared.bytes);
            }
            ++shared.users;
            store.xnr_segment_bytes = shared.segment;
            store.xnr_scratch_bytes = shared.bytes;
            store.xnr_scratch = shared.ptr;
        }
        if (!streamed.empty()) {   // the ring the streamed levels share (H2 8x)
            size_t ring = 0;
            for (int l : streamed) {
                if (auto k2 = store.kept.find(l); k2 != store.kept.end()) ring = std::max(ring, k2->second.ring_bytes);
            }
            store.stream_bytes = std::max<size_t>(ring, 1) + kStagedAlign;
            store.stream_block = DeviceHeap::instance().alloc_resident(store.stream_bytes, ResidentKind::KeptFactors);
            store.resident.clear();
        }
        double upload_seconds = 0.0, register_seconds = 0.0, pack_seconds = 0.0, wait_seconds = 0.0;
        double kept_bytes = 0.0, staged_bytes = 0.0;
        for (const auto& [level, kl] : store.kept) kept_bytes += kl.bytes;
        cudaEvent_t staged[2] = {nullptr, nullptr};
        int turn = 0;

        auto card_line = [&](const char* when, int level) {   // (H2 8x: every rank's card, per level of the set-up; BPACK_TRACE=wave)
            if (!env::trace("wave")) return;
            const DeviceHeap& h = DeviceHeap::instance();
            std::string held;
            for (size_t k = 0; k < static_cast<size_t>(ResidentKind::Count); ++k) {
                const size_t bytes = h.resident_now(static_cast<ResidentKind>(k));
                if (bytes < (size_t{256} << 20)) continue;
                char b[96];
                std::snprintf(b, sizeof b, "%s%s %.2f", held.empty() ? "" : ", ", resident_kind_name(static_cast<ResidentKind>(k)),
                              static_cast<double>(bytes) / 1e9);
                held += b;
            }
            std::printf("GPU solve: rank %d %s level %d: %.2f GB in use (largest free %.2f); %s\n", rank, when, level,
                        static_cast<double>(h.used()) / 1e9, static_cast<double>(h.largest_free()) / 1e9, held.c_str());
            std::fflush(stdout);
        };
        card_line("before", 0);
        for (int level = 2; level <= leaf; ++level) {
            auto& lvl = tree->levels[static_cast<size_t>(level)];
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (!lvl.is_process_active) continue;
            L.active = true;
            L.level = level;
            card_line("at", level);
            if (host_levels.count(level)) {  // (no factors, no tables)
                L.ca = true;
                L.host = true;
                // (a kept level solved on the host reads its factors there)
                if (auto kit = store.kept.find(level); kit != store.kept.end()) {
                    materialize_host_factors_level(tree, level, kit->second);
                }
                continue;
            }
            if (auto kit = store.kept.find(level); kit != store.kept.end() && kit->second.on_host()) {
                if (streamed.count(level) != 0) {
                    // streamed: its blocks' pointers to their slots of the ring
                    // now, the bytes at each sweep step (H2 8x)
                    KeptSolveLevel& kl = kit->second;
                    char* ring = reinterpret_cast<char*>(align_up(reinterpret_cast<uintptr_t>(store.stream_block), kStagedAlign));
                    std::vector<char*> targets(kl.staged.size());
                    for (size_t b = 0; b < kl.staged.size(); ++b) targets[b] = ring + kl.block_off[b];
                    rebase_kept_level(kl, targets);
                    L.streamed = true;
                    ++store.stream_levels;
                } else {
                    // staged to the host after its level: back on the device for the sweeps
                    const auto tu = clock::now();
                    staged_bytes += static_cast<double>(kit->second.staged_bytes);
                    upload_kept_level(kit->second, Context::instance().stream());
                    upload_seconds += std::chrono::duration<double>(clock::now() - tu).count();
                }
            }
            const size_t nb = lvl.local_boxes.size();
            auto is_local = [&](int64_t m) {
                return m >= lvl.local_morton_start && m < lvl.local_morton_start + static_cast<int64_t>(nb);
            };
            // the boxes of the level vector: the local ones, then (a replicated
            // CA level) the ghosts, in Morton order
            using Box = BoxData<CoordType, DataType>;
            L.ca = replicated(level);
            std::vector<const Box*> bx;
            bx.reserve(nb + (L.ca ? lvl.ghost_boxes.size() : 0));
            for (const Box& box : lvl.local_boxes) bx.push_back(&box);
            std::unordered_map<int64_t, size_t> ghost_at;  // Morton -> index in bx
            if (L.ca) {
                std::vector<std::pair<int64_t, size_t>> ghosts;
                for (size_t g = 0; g < lvl.ghost_boxes.size(); ++g) ghosts.emplace_back(lvl.ghost_boxes[g].morton_index, g);
                std::sort(ghosts.begin(), ghosts.end());
                for (const auto& [m, g] : ghosts) {
                    ghost_at.emplace(m, bx.size());
                    L.ghost_box.push_back(m);
                    bx.push_back(&lvl.ghost_boxes[g]);
                }
            }
            const size_t nw = bx.size();  // boxes the level eliminates (local, ghosts)
            auto writable = [&](int64_t m) -> int64_t {
                if (is_local(m)) return m - lvl.local_morton_start;
                const auto it = ghost_at.find(m);
                return it == ghost_at.end() ? -1 : static_cast<int64_t>(it->second);
            };

            // ---- other ranks' neighbors: sizes and skeletons, and who reads ours
            // (a CA level: the assisting boxes, read only, from the
            // factorization's post-elimination assisting gather)
            const auto tr = clock::now();
            L.readers.assign(nb, {});
            std::map<int64_t, std::vector<int64_t>> remote_skel;  // morton -> skeleton
            if (L.ca) {
                std::set<int64_t> wanted;
                for (const Box* box : bx) {
                    if (box->num_points == 0) continue;
                    for (int64_t m : box->one_hop) {
                        if (writable(m) < 0) wanted.insert(m);
                    }
                }
                const auto& assisting = store.ca_assisting[level];
                for (int64_t m : wanted) {
                    const auto it = assisting.find(m);
                    if (it == assisting.end()) {
                        throw std::runtime_error("device solve: assisting box " + std::to_string(m) + " unknown");
                    }
                    const int n = it->second.n;
                    L.ghost_of.emplace(m, static_cast<int>(L.ghosts.size()));
                    L.ghosts.push_back({m, -1, n, L.ghost_points});
                    L.ghost_points += n;
                    remote_skel[m].assign(it->second.skeleton.begin(), it->second.skeleton.end());
                }
            } else if (lvl.num_active_processes > 1) {
                L.peers = compute_one_hop_neighbor_ranks(tree, lvl, level);
                std::unordered_map<int, size_t> peer_index;
                for (size_t i = 0; i < L.peers.size(); ++i) peer_index.emplace(L.peers[i], i);
                std::vector<PeerMessage> out(L.peers.size()), in;
                std::map<int64_t, int> wanted;
                for (const auto& box : lvl.local_boxes) {
                    // (an empty box of an unstructured tree keeps unfiltered neighbor lists)
                    if (box.num_points == 0) continue;
                    for (int64_t m : box.one_hop) {
                        if (!is_local(m)) wanted.emplace(m, solve_detail::owner_of(tree, lvl, level, m));
                    }
                }
                for (const auto& [m, owner] : wanted) {
                    auto it = peer_index.find(owner);
                    if (it == peer_index.end()) throw std::runtime_error("device solve: neighbor on a non-neighbor rank");
                    out[it->second].head.push_back(m);
                }
                register_seconds += exchange_peer_messages(tree->comm, L.peers, out, in, 740);
                // replies: n, k, skeleton per requested box
                std::vector<PeerMessage> reply(L.peers.size()), got;
                for (size_t i = 0; i < L.peers.size(); ++i) {
                    for (int64_t m : in[i].head) {
                        if (!is_local(m)) throw std::runtime_error("device solve: request for a box of another rank");
                        const auto& box = lvl.local_boxes[static_cast<size_t>(m - lvl.local_morton_start)];
                        L.readers[static_cast<size_t>(m - lvl.local_morton_start)].push_back(static_cast<int>(i));
                        reply[i].head.push_back(m);
                        reply[i].head.push_back(box.num_points);
                        reply[i].head.push_back(static_cast<int64_t>(box.skeleton_indices.size()));
                        reply[i].head.insert(reply[i].head.end(), box.skeleton_indices.begin(), box.skeleton_indices.end());
                    }
                }
                register_seconds += exchange_peer_messages(tree->comm, L.peers, reply, got, 743);
                for (size_t i = 0; i < L.peers.size(); ++i) {
                    const auto& h = got[i].head;
                    for (size_t at = 0; at < h.size();) {
                        const int64_t m = h[at], n = h[at + 1], k = h[at + 2];
                        L.ghost_of.emplace(m, static_cast<int>(L.ghosts.size()));
                        L.ghosts.push_back({m, L.peers[i], static_cast<int>(n), L.ghost_points});
                        L.ghost_points += n;
                        remote_skel[m].assign(h.begin() + static_cast<int64_t>(at) + 3,
                                              h.begin() + static_cast<int64_t>(at) + 3 + k);
                        at += 3 + static_cast<size_t>(k);
                    }
                }
            }
            register_seconds += std::chrono::duration<double>(clock::now() - tr).count();

            // ---- layout: per box [T | X_SR | X_NR | LU] doubles, [ipiv | skel | red] ints
            // (lists only for a level whose factors stayed on the device)
            const KeptSolveLevel* kept = nullptr;
            if (auto kit = store.kept.find(level); kit != store.kept.end()) kept = &kit->second;
            std::vector<const KeptSolveBox*> kept_box(nw, nullptr);
            L.pending_xnr.assign(nw, DeviceSolveLevel::PackedPart{});
            L.pending_t.assign(nw, DeviceSolveLevel::PackedPart{});
            L.pending_xsr.assign(nw, DeviceSolveLevel::PackedPart{});
            std::vector<int64_t> base(nw);
            std::vector<solve_detail::BoxLayout> lay(nw);
            std::vector<int> ntot(nw, 0);
            std::vector<size_t> chunk_of(nw), offset_in(nw);
            std::vector<size_t> chunk_bytes;
            constexpr size_t kChunk = size_t{256} << 20;
            for (size_t b = 0; b < nw; ++b) {
                const Box& box = *bx[b];
                if (b == nb) L.local_points = L.points;
                base[b] = L.points;
                L.points += box.num_points;
                const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                const bool eliminated = k > 0 && r > 0;
                solve_detail::BoxLayout& l = lay[b];
                // (a kept box's X_NR, X_SR and LU may have no host copy: its
                // shapes from the kept one)
                const KeptSolveBox* kb = nullptr;
                if (kept != nullptr && eliminated) {
                    auto it = kept->boxes.find(box.morton_index);
                    if (it != kept->boxes.end()) kb = &it->second;
                }
                ntot[b] = eliminated && box.X_NR.is_allocated() && !box.one_hop.empty() ? static_cast<int>(box.X_NR.rows) : 0;
                if (kb != nullptr && !box.X_NR.is_allocated() && !box.one_hop.empty()) ntot[b] = kb->ntot;
                l.has_T = eliminated && (kb != nullptr ? kb->T != nullptr : box.interpolation_matrix.is_allocated());
                l.has_xsr = eliminated && (kb != nullptr ? kb->xsr != nullptr : box.X_SR.is_allocated());
                // a nonsymmetric matrix: X_RS (r x k) and X_RN (r x ntot) beside them
                l.has_xrs = nonsym && eliminated && (kb != nullptr ? kb->xrs_t != nullptr : box.X_RS.is_allocated());
                if (nonsym && eliminated && kb == nullptr &&
                    ((l.has_xrs && (box.X_RS.rows != static_cast<int64_t>(r) || box.X_RS.cols != static_cast<int64_t>(k))) ||
                     (ntot[b] > 0 && (!box.X_RN.is_allocated() || box.X_RN.rows != static_cast<int64_t>(r) ||
                                      box.X_RN.cols != ntot[b])))) {
                    throw std::runtime_error("device solve: X_RS / X_RN of box " + std::to_string(box.morton_index) +
                                             " do not match");
                }
                if (kept != nullptr && eliminated) {
                    if (kb == nullptr || kb->k != static_cast<int>(k) || kb->r != static_cast<int>(r) ||
                        kb->ntot != ntot[b]) {
                        throw std::runtime_error("device solve: kept factors of box " + std::to_string(box.morton_index) +
                                                 " do not match");
                    }
                    kept_box[b] = kb;
                }
                size_t at = 0;
                if (eliminated && kept_box[b] == nullptr) {
                    l.T = at;    at = align_up(at + (l.has_T ? k * r * D : 0));
                    l.xsr = at;  at = align_up(at + (l.has_xsr ? k * r * D : 0));
                    l.xnr = at;  at = align_up(at + static_cast<size_t>(ntot[b]) * r * D);
                    if (nonsym) {
                        l.xrs_t = at;  at = align_up(at + (l.has_xrs ? k * r * D : 0));
                        l.xrn_t = at;  at = align_up(at + static_cast<size_t>(ntot[b]) * r * D);
                    }
                    l.lu = at;   at = align_up(at + r * r * D);
                    l.ipiv = at; at = align_up(at + r * I);
                }
                l.skel = at; at = align_up(at + k * I);
                l.red = at;  at = align_up(at + r * I);
                l.bytes = at;
                if (chunk_bytes.empty() || (chunk_bytes.back() > 0 && chunk_bytes.back() + at > kChunk)) chunk_bytes.push_back(0);
                chunk_of[b] = chunk_bytes.size() - 1;
                offset_in[b] = chunk_bytes.back();
                chunk_bytes.back() += at;
                L.max_n = std::max(L.max_n, static_cast<int>(box.num_points));
                L.max_r = std::max(L.max_r, static_cast<int>(r));
            }
            if (nw == nb) L.local_points = L.points;
            std::vector<size_t> chunk_first(chunk_bytes.size() + 1, nw);  // boxes of chunk c: [first[c], first[c + 1])
            for (size_t b = nw; b-- > 0;) chunk_first[chunk_of[b]] = b;
            std::vector<char*> chunks(chunk_bytes.size());
            for (size_t c = 0; c < chunks.size(); ++c) {
                chunks[c] = heap.alloc_resident(std::max<size_t>(chunk_bytes[c], 1));
                L.blocks.push_back(chunks[c]);
            }
            // ---- upload, chunk by chunk through two pinned staging buffers
            const auto tu = clock::now();
            for (size_t c = 0; c < chunks.size(); ++c) {
                if (staged[turn] == nullptr) {
                    check_cuda(cudaEventCreateWithFlags(&staged[turn], cudaEventDisableTiming), "cudaEventCreate");
                } else {
                    const auto tw = clock::now();
                    check_cuda(cudaEventSynchronize(staged[turn]), "solve upload");
                    wait_seconds += std::chrono::duration<double>(clock::now() - tw).count();
                }
                char* h = static_cast<char*>(pinned_pool().staging[turn].reserve(std::max<size_t>(chunk_bytes[c], 1)));
                const auto tp = clock::now();
                // pieces of about 1 MB (a matrix's column range, or a box's lists), so all threads copy
                struct Piece { size_t box; int what; size_t col0, col1; };  // what: 0 T, 1 X_SR, 2 X_NR, 3 LU, 4 lists,
                                                                            // 5 X_RS^T, 6 X_RN^T
                std::vector<Piece> pieces;
                for (size_t b = chunk_first[c]; b < chunk_first[c + 1]; ++b) {
                    const Box& box = *bx[b];
                    const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                    pieces.push_back(Piece{b, 4, 0, 0});
                    if (k == 0 || r == 0 || kept_box[b] != nullptr) continue;
                    auto split = [&](int what, size_t rows) {
                        const size_t step = std::max<size_t>(1, (size_t{1} << 20) / std::max<size_t>(rows * D, 1));
                        for (size_t j = 0; j < r; j += step) pieces.push_back(Piece{b, what, j, std::min(r, j + step)});
                    };
                    if (lay[b].has_T) split(0, k);
                    if (lay[b].has_xsr) split(1, k);
                    if (ntot[b] > 0) split(2, static_cast<size_t>(ntot[b]));
                    if (lay[b].has_xrs) split(5, k);
                    if (nonsym && ntot[b] > 0) split(6, static_cast<size_t>(ntot[b]));
                    split(3, r);
                }
                #pragma omp parallel for schedule(dynamic, 1)
                for (int64_t pi = 0; pi < static_cast<int64_t>(pieces.size()); ++pi) {
                    const Piece& pc = pieces[static_cast<size_t>(pi)];
                    const Box& box = *bx[pc.box];
                    const solve_detail::BoxLayout& l = lay[pc.box];
                    char* p = h + offset_in[pc.box];
                    const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                    auto put = [&](size_t at, const MatrixStorage<DataType>& m, size_t rows) {
                        DataType* d = reinterpret_cast<DataType*>(p + at);
                        for (size_t j = pc.col0; j < pc.col1; ++j)
                            std::memcpy(d + j * rows, m.data.data() + j * static_cast<size_t>(m.lda), rows * D);
                    };
                    // m^T (m: cols x rows, ld m.lda) as rows x r, columns col0 .. col1
                    auto put_t = [&](size_t at, const MatrixStorage<DataType>& m, size_t rows) {
                        DataType* d = reinterpret_cast<DataType*>(p + at);
                        const size_t ld = static_cast<size_t>(m.lda);
                        for (size_t j = pc.col0; j < pc.col1; ++j)
                            for (size_t i = 0; i < rows; ++i) d[i + j * rows] = m.data[j + i * ld];
                    };
                    switch (pc.what) {
                        case 0: put(l.T, box.interpolation_matrix, k); break;
                        case 1: put(l.xsr, box.X_SR, k); break;
                        case 2: put(l.xnr, box.X_NR, static_cast<size_t>(ntot[pc.box])); break;
                        case 3: put(l.lu, box.X_RR, r); break;
                        case 5: put_t(l.xrs_t, box.X_RS, k); break;
                        case 6: put_t(l.xrn_t, box.X_RN, static_cast<size_t>(ntot[pc.box])); break;
                        default: {
                            if (k > 0 && r > 0 && kept_box[pc.box] == nullptr) {
                                std::memcpy(p + l.ipiv, box.X_RR_pivots.data(), r * I);
                            }
                            int* sk = reinterpret_cast<int*>(p + l.skel);
                            for (size_t i = 0; i < k; ++i) sk[i] = static_cast<int>(box.skeleton_indices[i]);
                            int* rd = reinterpret_cast<int*>(p + l.red);
                            for (size_t i = 0; i < r; ++i) rd[i] = static_cast<int>(box.redundant_indices[i]);
                        }
                    }
                }
                pack_seconds += std::chrono::duration<double>(clock::now() - tp).count();
                check_cuda(cudaMemcpyAsync(chunks[c], h, chunk_bytes[c], cudaMemcpyHostToDevice, stream), "solve upload");
                check_cuda(cudaEventRecord(staged[turn], stream), "cudaEventRecord");
                turn ^= 1;
                store.bytes += static_cast<double>(chunk_bytes[c]);
            }
            upload_seconds += std::chrono::duration<double>(clock::now() - tu).count();

            // ---- boxes and slots
            L.boxes.resize(nw);
            std::vector<int> remote_lists;  // skeletons of other ranks' neighbors
            std::unordered_map<int64_t, size_t> remote_list_at;
            for (const auto& [m, sk] : remote_skel) {
                remote_list_at.emplace(m, remote_lists.size());
                for (int64_t x : sk) remote_lists.push_back(static_cast<int>(x));
            }
            for (size_t b = 0; b < nw; ++b) {
                const Box& box = *bx[b];
                const solve_detail::BoxLayout& l = lay[b];
                const char* p = chunks[chunk_of[b]] + offset_in[b];
                SolveBox& sb = L.boxes[b];
                sb.vec = base[b];
                sb.n = static_cast<int>(box.num_points);
                sb.k = static_cast<int>(box.skeleton_indices.size());
                sb.r = static_cast<int>(box.redundant_indices.size());
                sb.ntot = ntot[b];
                sb.skel = reinterpret_cast<const int*>(p + l.skel);
                sb.red = reinterpret_cast<const int*>(p + l.red);
                const bool eliminated = sb.k > 0 && sb.r > 0;
                if (const KeptSolveBox* kb = kept_box[b]) {
                    if (kb->xnr_pending) {
                        throw std::runtime_error("H2 GPU solve: the kept X_NR of box " +
                                                 std::to_string(box.morton_index) +
                                                 " was never transposed out of its fill source");
                    }
                    if (device_solve_check() && nonsym && kb->xnr_bytes == 0) {
                        // BPACK_CHECK=solve: a nonsymmetric box's kept factors against
                        // its host copies (kept beside them), the largest difference
                        // relative to the largest entry, per factor
                        auto diff = [&](const void* dev, const MatrixStorage<DataType>& m, bool transposed,
                                        size_t rows, size_t cols) -> double {   // dev: rows x cols
                            if (dev == nullptr || !m.is_allocated() || rows * cols == 0) return -1.0;
                            std::vector<DataType> d(rows * cols);
                            check_cuda(cudaMemcpy(d.data(), dev, d.size() * D, cudaMemcpyDeviceToHost), "kept check");
                            double e = 0.0, a = 0.0;
                            for (size_t j = 0; j < cols; ++j) {
                                for (size_t i = 0; i < rows; ++i) {
                                    const DataType h = transposed ? m.data[j + i * static_cast<size_t>(m.lda)]
                                                                  : m.data[i + j * static_cast<size_t>(m.lda)];
                                    e = std::max(e, static_cast<double>(std::abs(d[i + j * rows] - h)));
                                    a = std::max(a, static_cast<double>(std::abs(h)));
                                }
                            }
                            return a > 0.0 ? e / a : e;
                        };
                        const size_t kk = static_cast<size_t>(kb->k), rr = static_cast<size_t>(kb->r),
                                     nn = static_cast<size_t>(kb->ntot);
                        const double e_t = diff(kb->T, box.interpolation_matrix, false, kk, rr);
                        const double e_sr = diff(kb->xsr, box.X_SR, false, kk, rr);
                        const double e_nr = diff(kb->xnr, box.X_NR, false, nn, rr);
                        const double e_rs = diff(kb->xrs_t, box.X_RS, true, kk, rr);
                        const double e_rn = diff(kb->xrn_t, box.X_RN, true, nn, rr);
                        if (std::max({e_t, e_sr, e_nr, e_rs, e_rn}) > 1e-12) {
                            std::printf("  [gpu] level %d box %lld kept vs host: T %.2e, X_SR %.2e, X_NR %.2e, "
                                        "X_RS %.2e, X_RN %.2e\n", level, static_cast<long long>(box.morton_index), e_t,
                                        e_sr, e_nr, e_rs, e_rn);
                            std::fflush(stdout);
                        }
                    }
                    sb.T = l.has_T ? kb->T : nullptr;
                    sb.xsr = l.has_xsr ? kb->xsr : nullptr;
                    sb.xnr = kb->xnr;
                    sb.xrs_t = l.has_xrs ? kb->xrs_t : nullptr;
                    sb.xrn_t = nonsym ? kb->xrn_t : nullptr;
                    if (nonsym && sb.ntot > 0 && sb.xrn_t == nullptr) {
                        throw std::runtime_error("device solve: kept box " + std::to_string(box.morton_index) +
                                                 " has no X_RN");
                    }
                    sb.lu = kb->lu;
                    sb.ipiv = kb->ipiv;
                    // A compressed factor cannot be read in place: note it,
                    // and the sweeps decompress its wave into L.xnr_scratch.
                    auto pending = [&](DeviceSolveLevel::PackedPart& part, const void* src, size_t bytes,
                                       size_t raw_bytes, uint32_t chunk0, uint32_t nchunks) {
                        if (bytes == 0 || raw_bytes == 0) return;
                        L.xnr_compressed = true;
                        part.src = src;
                        part.raw_bytes = raw_bytes;
                        part.strides = kept->chunk_strides.data() + chunk0;
                        part.nchunks = nchunks;
                    };
                    const size_t rr = static_cast<size_t>(sb.r), kk = static_cast<size_t>(sb.k);
                    pending(L.pending_xnr[b], kb->xnr, kb->xnr_bytes, static_cast<size_t>(ntot[b]) * rr * D,
                            kb->xnr_chunk0, kb->xnr_nchunks);
                    if (l.has_T) {
                        pending(L.pending_t[b], kb->T, kb->t_bytes, kk * rr * D, kb->t_chunk0, kb->t_nchunks);
                    }
                    if (l.has_xsr) {
                        pending(L.pending_xsr[b], kb->xsr, kb->xsr_bytes, kk * rr * D, kb->xsr_chunk0,
                                kb->xsr_nchunks);
                    }
                } else {
                    sb.T = eliminated && l.has_T ? p + l.T : nullptr;
                    sb.xsr = eliminated && l.has_xsr ? p + l.xsr : nullptr;
                    sb.xnr = eliminated ? p + l.xnr : nullptr;
                    sb.xrs_t = eliminated && l.has_xrs ? p + l.xrs_t : nullptr;
                    sb.xrn_t = eliminated && nonsym ? p + l.xrn_t : nullptr;
                    sb.lu = eliminated ? p + l.lu : nullptr;
                    sb.ipiv = eliminated ? reinterpret_cast<const int*>(p + l.ipiv) : nullptr;
                }
                sb.slot0 = static_cast<int>(L.slots.size());
                sb.nslots = 0;
                sb.work = 0;
                if (ntot[b] == 0) continue;  // no neighbor terms (or no elimination: the vector only)
                if (box.use_full_set.size() != box.one_hop.size()) {
                    throw std::runtime_error("device solve: use_full_set missing for box " + std::to_string(box.morton_index));
                }
                // neighbor counts: checked against those cached at factorization
                static const std::vector<int64_t> none;
                const auto& cached = b < nb && b < lvl.solve_neighbor_size.size() ? lvl.solve_neighbor_size[b] : none;
                int row = 0;
                for (size_t a = 0; a < box.one_hop.size(); ++a) {
                    const int64_t m = box.one_hop[a];
                    const bool full = box.use_full_set[a] == 1;
                    SolveSlot sl{};
                    sl.full = full ? 1 : 0;
                    sl.row0 = row;
                    if (const int64_t wi = writable(m); wi >= 0) {  // (a CA level: local or ghost)
                        const size_t nbx = static_cast<size_t>(wi);
                        const Box& nbox = *bx[nbx];
                        sl.ghost = 0;
                        sl.vec = base[nbx];
                        sl.n = static_cast<int>(nbox.num_points);
                        sl.count = full ? sl.n : static_cast<int>(nbox.skeleton_indices.size());
                        sl.skel = reinterpret_cast<const int*>(chunks[chunk_of[nbx]] + offset_in[nbx] + lay[nbx].skel);
                    } else {
                        auto it = L.ghost_of.find(m);
                        if (it == L.ghost_of.end()) throw std::runtime_error("device solve: neighbor " + std::to_string(m) + " unknown");
                        const auto& g = L.ghosts[static_cast<size_t>(it->second)];
                        sl.ghost = 1;
                        sl.vec = g.offset;
                        sl.n = g.n;
                        sl.count = full ? g.n : static_cast<int>(remote_skel.at(m).size());
                        sl.skel = nullptr;  // set once the remote lists are on the device
                    }
                    if (cached.size() == box.one_hop.size() && cached[a] != sl.count) {
                        throw std::runtime_error("device solve: neighbor size mismatch for box " +
                                                 std::to_string(box.morton_index));
                    }
                    row += sl.count;
                    L.slots.push_back(sl);
                    L.slot_morton.push_back(m);
                    ++sb.nslots;
                }
                if (row != ntot[b]) throw std::runtime_error("device solve: X_NR rows do not match the neighbors");
            }
            // ---- the sweeps' tables
            store.bytes += static_cast<double>(
                L.ca ? build_ca_solve_tables(tree, lvl, L, nb, ghost_at, remote_lists, remote_list_at)
                     : build_solve_tables(tree, lvl, level, L, remote_lists, remote_list_at, register_seconds));
        }
        check_cuda(cudaStreamSynchronize(stream), "solve upload");
        for (cudaEvent_t e : staged) {
            if (e != nullptr) cudaEventDestroy(e);
        }
        // the replicated CA levels' halos in device memory (all ranks)
        for (int level = 2; level <= leaf; ++level) {
            if (!replicated(level) || host_levels.count(level) != 0) continue;
            store.bytes += static_cast<double>(
                build_ca_halo(tree, level, store.levels[static_cast<size_t>(level)], register_seconds));
        }
        store.usable = true;
        double seconds = std::chrono::duration<double>(clock::now() - t0).count();
        double bytes = store.bytes;
        MPI_Allreduce(MPI_IN_PLACE, &seconds, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &kept_bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &staged_bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        std::string kept_levels, on_host;
        for (const auto& [level, kl] : store.kept) kept_levels += (kept_levels.empty() ? "" : ",") + std::to_string(level);
        for (int level : host_levels) on_host += (on_host.empty() ? "" : ",") + std::to_string(level);
        (void)verbosity;
        if (rank == 0) {  // once per factorization, also for a silent first solve
            std::printf("GPU solve: set up in %.2f s; factors kept during the factorization: %s "
                        "(up to %.2f GB per rank, of which %.2f GB staged to the host at the levels' ends and brought "
                        "back here), uploaded now: up to %.2f GB per rank (upload %.2f s [host packing "
                        "%.2f, waiting for the copies %.2f], neighbor registration %.2f s on rank 0)\n",
                        seconds, kept_levels.empty() ? "none" : ("levels " + kept_levels).c_str(), kept_bytes / 1e9,
                        staged_bytes / 1e9, bytes / 1e9, upload_seconds, pack_seconds, wait_seconds, register_seconds);
            if (!host_levels.empty()) {
                std::string sizes;
                for (int level : host_levels) {
                    char buf[64];
                    std::snprintf(buf, sizeof(buf), "%s%d: %.2f GB", sizes.empty() ? "" : ", ", level,
                                  host_only[static_cast<size_t>(level)] ? 0.0 : need_max[static_cast<size_t>(level)] / 1e9);
                    sizes += buf;
                }
                std::printf("GPU solve: CA levels %s on the host (%s; device room %.2f GB per rank, "
                            "resident %.2f GB)\n",
                            on_host.c_str(), sizes.c_str(), room / 1e9, planned / 1e9);
            }
            std::fflush(stdout);
        }
        return true;
    }
}

// ---------------------------------------------------------------------------
// The sweeps
// ---------------------------------------------------------------------------
struct DeviceSolveTimes {
    double forward = 0.0, diagonal = 0.0, backward = 0.0, transfer = 0.0, comm = 0.0;
    double halo = 0.0;  // CA halo packing and unpacking (device copies)
    int64_t direct_messages = 0, staged_messages = 0;  // exchanges moving data: in device memory, through the host
    // H2_compress_factors: the waves' X_NR decompression, timed on the stream
    // (CUDA events, so the sweeps are not serialized to measure it)
    double xnr = 0.0;
    int64_t xnr_calls = 0;
    double xnr_bytes = 0.0;  // uncompressed bytes produced
};

template<typename CoordType, typename DataType>
class DeviceSolveRun {
    using S = typename DeviceScalar<DataType>::type;  // vector elements on the device
    static constexpr int kWords = static_cast<int>(sizeof(S) / sizeof(double));  // MPI_DOUBLE per element

public:
    DeviceSolveRun(ParallelTree<CoordType, DataType>* tree, int nrhs) : tree_(tree), nrhs_(nrhs) {}
    ~DeviceSolveRun() {
        for (XnrSlot& s : xnr_slots_) {
            if (!s.made) continue;
            cudaEventDestroy(s.start);
            cudaEventDestroy(s.stop);
        }
        if (xnr_stream_ != nullptr) {
            cudaStreamSynchronize(xnr_stream_);
            for (size_t b = 0; b < kXnrBuffers; ++b) {
                cudaEventDestroy(xnr_done_[b]);
                cudaEventDestroy(xnr_used_[b]);
            }
            for (size_t b = 0; b < 2; ++b) {
                if (kept_ready_[b] != nullptr) cudaEventDestroy(kept_ready_[b]);
            }
            if (kept_done_ != nullptr) {
                cudaEventDestroy(kept_done_);
            }
            cudaStreamDestroy(xnr_stream_);
        }
    }
    DeviceSolveRun(const DeviceSolveRun&) = delete;
    DeviceSolveRun& operator=(const DeviceSolveRun&) = delete;

    // H2_compress_factors: a wave's X_NR, decompressed into the level's

    // scratch and its SolveBox pointers moved there, before the wave's kernels.

    // Everything is on one stream, so wave w's kernels see w's pointers even

    // though wave w+1 overwrites them afterwards.

    // A ring of event pairs around the decompression calls.  A slot is read
    // back only when it comes round again, hundreds of calls later, so the
    // measurement costs no synchronization of its own.
    static constexpr size_t kXnrSlots = 256;
    struct XnrSlot { cudaEvent_t start{}, stop{}; bool made = false, pending = false; double bytes = 0.0; };
    std::vector<XnrSlot> xnr_slots_{kXnrSlots};
    size_t xnr_at_ = 0;

    void xnr_slot_collect(XnrSlot& s) {
        if (!s.pending) return;
        float ms = 0.0f;
        check_cuda(cudaEventSynchronize(s.stop), "xnr timer");
        check_cuda(cudaEventElapsedTime(&ms, s.start, s.stop), "xnr timer");
        times.xnr += static_cast<double>(ms) * 1e-3;
        times.xnr_bytes += s.bytes;
        ++times.xnr_calls;
        s.pending = false;
    }
    int xnr_timer_begin(cudaStream_t on) {
        XnrSlot& s = xnr_slots_[xnr_at_];
        xnr_slot_collect(s);
        if (!s.made) {
            check_cuda(cudaEventCreate(&s.start), "xnr timer");
            check_cuda(cudaEventCreate(&s.stop), "xnr timer");
            s.made = true;
        }
        check_cuda(cudaEventRecord(s.start, on), "xnr timer");
        return static_cast<int>(xnr_at_);
    }
    void xnr_timer_end(int slot, double bytes, cudaStream_t on) {
        XnrSlot& s = xnr_slots_[static_cast<size_t>(slot)];
        check_cuda(cudaEventRecord(s.stop, on), "xnr timer");
        s.pending = true;
        s.bytes = bytes;
        xnr_at_ = (static_cast<size_t>(slot) + 1) % kXnrSlots;
    }
    // The buffer rotation (prepare_wave_xnr): which wave each buffer holds,
    // whether its last reader has been marked, and the wave before this one.
    static constexpr size_t kNoWave = static_cast<size_t>(-1);
    cudaStream_t xnr_stream_ = nullptr;
    std::array<cudaEvent_t, kXnrBuffers> xnr_done_{};   // its decompression finished
    std::array<cudaEvent_t, kXnrBuffers> xnr_used_{};   // its reader finished
    std::array<bool, kXnrBuffers> xnr_used_valid_{};
    std::array<size_t, kXnrBuffers> xnr_ready_{kNoWave, kNoWave};
    size_t xnr_prev_ = kNoWave;
    const DeviceSolveLevel* xnr_level_ = nullptr;
    // Every pending slot, at the end of the sweeps.
    void xnr_timers_drain() {
        for (XnrSlot& s : xnr_slots_) xnr_slot_collect(s);
    }

    // Wave `w` of a sweep that goes on to wave w + step (step +1 forward, -1
    // backward): its X_NR is decompressed in the buffer of its parity, and
    // the next wave's decompression is started here, on a stream of its own,
    // so it runs under this wave's kernels and messages.  The buffer the next
    // wave will overwrite was read by the wave before this one, so the
    // decompression stream waits for that read first.
    void prepare_wave_xnr(DeviceSolveLevel& L, size_t w, int step) {
#ifdef H2_HAVE_NVCOMP

        if (!L.xnr_compressed || w >= L.xnr_waves.size()) return;

        SolveXnrWave& xw = L.xnr_waves[w];

        if (xw.empty()) return;

        if (xnr_level_ != &L) {  // another level: its buffers hold nothing of ours

            xnr_sync_buffers();

            xnr_level_ = &L;

            xnr_ready_.fill(kNoWave);

            xnr_used_valid_.fill(false);

            xnr_prev_ = kNoWave;

        }

        const size_t b = w % L.xnr_buffers;

        // The previous wave's kernels are on the stream by now, so an event
        // recorded here is after its reads of its buffer.

        if (xnr_prev_ != kNoWave) {

            const size_t pb = xnr_prev_ % L.xnr_buffers;

            check_cuda(cudaEventRecord(xnr_used_[pb], stream_), "xnr buffer use");

            xnr_used_valid_[pb] = true;

        }

        if (xnr_ready_[b] != w) xnr_decompress(L, w);  // not prefetched (the sweep's first wave)

        check_cuda(cudaStreamWaitEvent(stream_, xnr_done_[b], 0), "xnr buffer ready");

        // d_boxes is const for the kernels; the table block is writable.  The

        // wave's box indices and destinations are already on the device, so

        // this is one small kernel, not a copy of the whole table.

        launch_patch_xnr(const_cast<SolveBox*>(L.d_boxes), xw.d_box_idx, xw.d_box_dsts, xw.d_box_t_dsts,
                         xw.d_box_xsr_dsts, static_cast<int>(xw.boxes.size()), stream_);

        // With one buffer there is nowhere to put the next wave: it would
        // overwrite this one before its kernels have run.

        const size_t next = w + static_cast<size_t>(step);

        if (L.xnr_buffers > 1 && step != 0 && next < L.xnr_waves.size() && !L.xnr_waves[next].empty()) {

            xnr_decompress(L, next);

        }

        xnr_prev_ = w;
#else

        (void)L; (void)w; (void)step;
#endif

    }

#ifdef H2_HAVE_NVCOMP
    // Wave `w` into the buffer of its parity, on the decompression stream,
    // once whoever read that buffer last is done with it.
    void xnr_decompress(DeviceSolveLevel& L, size_t w) {
        const size_t b = w % L.xnr_buffers;
        xnr_make_resources();
        if (xnr_used_valid_[b]) {
            check_cuda(cudaStreamWaitEvent(xnr_stream_, xnr_used_[b], 0), "xnr buffer free");
            xnr_used_valid_[b] = false;
        }
        const int slot = xnr_timer_begin(xnr_stream_);
        decompress_solve_xnr_wave(L.xnr_waves[w], xnr_stream_);
        xnr_timer_end(slot, L.xnr_waves[w].raw_bytes, xnr_stream_);
        check_cuda(cudaEventRecord(xnr_done_[b], xnr_stream_), "xnr buffer ready");
        xnr_ready_[b] = w;
    }
    void xnr_make_resources() {
        if (xnr_stream_ != nullptr) return;
        check_cuda(cudaStreamCreateWithFlags(&xnr_stream_, cudaStreamNonBlocking), "xnr stream");
        for (size_t b = 0; b < kXnrBuffers; ++b) {
            check_cuda(cudaEventCreateWithFlags(&xnr_done_[b], cudaEventDisableTiming), "xnr events");
            check_cuda(cudaEventCreateWithFlags(&xnr_used_[b], cudaEventDisableTiming), "xnr events");
        }
        for (size_t b = 0; b < 2; ++b) check_cuda(cudaEventCreateWithFlags(&kept_ready_[b], cudaEventDisableTiming), "kept events");
        check_cuda(cudaEventCreateWithFlags(&kept_done_, cudaEventDisableTiming), "kept events");
    }
    // H2 8x: a streamed level's blocks for wave `w` in their slots before
    // its kernels, and the next wave's (w + step) copied in under them on
    // the decompression stream, which also orders them before the next
    // wave's X_NR decompression.  The copies for the next wave wait for the
    // kernels enqueued so far (the previous wave's), which is what the ring
    // guarantees they may overwrite.
    cudaEvent_t kept_ready_[2] = {nullptr, nullptr};
    cudaEvent_t kept_done_ = nullptr;
    int kept_level_ = -1;
    size_t kept_next_ = kNoWave;
    void prepare_wave_kept(DeviceSolveLevel& L, size_t w, int step) {
        DeviceSolveStore& store = device_solve_store();
        auto it = store.kept.find(L.level);
        if (it == store.kept.end() || !it->second.streamed) return;
        KeptSolveLevel& kl = it->second;
        xnr_make_resources();
        const auto t0 = std::chrono::steady_clock::now();
        if (kept_level_ != L.level || kept_next_ != w) {   // not prefetched: this wave now
            kept_wave_resident(store, L.level, kl, w, xnr_stream_, false);
            check_cuda(cudaEventRecord(kept_ready_[w % 2], xnr_stream_), "kept ready");
        }
        check_cuda(cudaStreamWaitEvent(stream_, kept_ready_[w % 2], 0), "kept ready");
        check_cuda(cudaEventRecord(kept_done_, stream_), "kept done");
        // the next wave that reads blocks: an owner level's home order leaves
        // a rank waves without boxes, whose neighbours the X_NR prefetch
        // reads across (its segments skip them), and the ring keeps the
        // reading waves apart (plan_kept_streaming)
        size_t next = kNoWave;
        if (step > 0) {
            for (size_t v = w + 1; v < kl.wave_blocks.size() && next == kNoWave; ++v) {
                if (!kl.wave_blocks[v].empty()) next = v;
            }
        } else {
            for (size_t v = std::min(w, kl.wave_blocks.size()); v-- > 0 && next == kNoWave;) {
                if (!kl.wave_blocks[v].empty()) next = v;
            }
        }
        if (next != kNoWave) {
            check_cuda(cudaStreamWaitEvent(xnr_stream_, kept_done_, 0), "kept done");
            kept_wave_resident(store, L.level, kl, next, xnr_stream_, false);
            check_cuda(cudaEventRecord(kept_ready_[next % 2], xnr_stream_), "kept ready");
            kept_next_ = next;
        } else {
            kept_next_ = kNoWave;
        }
        kept_level_ = L.level;
        store.stream_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    // Before another level takes the buffers over.
    void xnr_sync_buffers() {
        if (xnr_stream_ == nullptr) return;
        check_cuda(cudaStreamSynchronize(xnr_stream_), "xnr buffers");
    }
#else
    void xnr_sync_buffers() {}
#endif


    // a level's vectors: host solve data -> device
    S* upload(DeviceSolveLevel& L, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t len = static_cast<size_t>(L.points) * nrhs_;
        S* d = heap_.alloc<S>(std::max<size_t>(len, 1) * sizeof(S));
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        for (size_t b = 0; b < L.boxes.size(); ++b) {
            const auto& v = data[b].left_side;
            if (v.size() != static_cast<size_t>(L.boxes[b].n) * nrhs_) throw std::runtime_error("device solve: vector size");
            std::memcpy(h + L.boxes[b].vec * nrhs_, v.data(), v.size() * sizeof(DataType));
        }
        check_cuda(cudaMemcpyAsync(d, h, len * sizeof(S), cudaMemcpyHostToDevice, stream_), "solve upload");
        check_cuda(cudaStreamSynchronize(stream_), "solve upload");
        times.transfer += seconds_since(t0);
        return d;
    }
    void download(DeviceSolveLevel& L, S* d, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t len = static_cast<size_t>(L.points) * nrhs_;
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        check_cuda(cudaMemcpyAsync(h, d, len * sizeof(S), cudaMemcpyDeviceToHost, stream_), "solve download");
        check_cuda(cudaStreamSynchronize(stream_), "solve download");
        for (size_t b = 0; b < L.boxes.size(); ++b) {
            auto& v = data[b].left_side;
            std::memcpy(v.data(), h + L.boxes[b].vec * nrhs_, v.size() * sizeof(DataType));
        }
        heap_.free(d);
        times.transfer += seconds_since(t0);
    }

    // The forward waves of a level, then its diagonal solves.
    void forward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* outbox = alloc_message(L.out_max);
        S* inbox = alloc_message(L.in_max);
        for (size_t wi = 0; wi < L.waves.size(); ++wi) {
            const auto& W = L.waves[wi];
            prepare_wave_kept(L, wi, +1);
            for (size_t s = 0; s < W.segs.size(); ++s) {  // (H2 8x: the wave's segments, each through the scratch)
                prepare_wave_xnr(L, W.seg0 + s, +1);
                launch_solve_forward(L.d_boxes, L.d_order + W.order0 + W.segs[s].first, W.segs[s].count, vec, work, nrhs_, L.max_n, stream_);
            }
            if (L.streamed) launch_solve_diagonal(L.d_boxes, L.d_diag + W.diag0, W.ndiag, vec, nrhs_, L.max_r, stream_);  // (with its wave: H2 8x)
            launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, outbox, inbox,
                               stream_);
            if (L.peers.empty()) continue;
            exchange(L, W.out, W.in, outbox, inbox, 750);
            launch_solve_accum(L.d_items + W.in_item0, L.d_parts, W.in_nitems, W.in_max_count, nrhs_, vec, work, outbox,
                               inbox, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "forward sweep");
        free_points(work);
        free_message(outbox);
        free_message(inbox);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        if (!L.streamed) launch_solve_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal solves");
        times.diagonal += seconds_since(t1);
    }

    // The backward waves of a level, in reverse, other ranks' neighbor
    // vectors refreshed before each.
    void backward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* ghost = alloc_points(L.ghost_points);
        S* sendbox = alloc_message(L.send_max);
        S* inbox = alloc_message(L.recv_max);
        for (int w = static_cast<int>(L.waves.size()) - 1; w >= 0; --w) {
            const auto& W = L.waves[static_cast<size_t>(w)];
            if (!L.peers.empty()) refresh(L, W.ref[0], vec, ghost, sendbox, inbox);
            prepare_wave_kept(L, static_cast<size_t>(w), -1);
            for (size_t s = W.segs.size(); s-- > 0;) {  // (H2 8x: the wave's segments, in reverse)
                prepare_wave_xnr(L, W.seg0 + s, -1);
                launch_solve_backward(L.d_boxes, L.d_slots, L.d_rorder + W.order0 + (W.count - W.segs[s].first - W.segs[s].count), W.segs[s].count, vec, ghost, work, nrhs_, L.max_n,
                                  W.max_ntot, stream_);
            }
        }
        check_cuda(cudaStreamSynchronize(stream_), "backward sweep");
        free_points(work);
        free_points(ghost);
        free_message(sendbox);
        free_message(inbox);
        times.backward += seconds_since(t0);
    }

    // The multiply F x (core/apply_mul.hpp): a level's forward W waves
    // in order (reading neighbors, other ranks' vectors refreshed before
    // each), then its diagonal multiplies (x_R final by then).
    void mul_forward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* ghost = alloc_points(L.ghost_points);
        S* sendbox = alloc_message(L.send_max);
        S* inbox = alloc_message(L.recv_max);
        for (size_t wi = 0; wi < L.waves.size(); ++wi) {
            const auto& W = L.waves[wi];
            if (!L.peers.empty()) refresh(L, W.ref[1], vec, ghost, sendbox, inbox);
            prepare_wave_kept(L, wi, +1);
            for (size_t s = 0; s < W.segs.size(); ++s) {  // (H2 8x: the wave's segments, each through the scratch)
                prepare_wave_xnr(L, W.seg0 + s, +1);
                launch_mul_forward(L.d_boxes, L.d_slots, L.d_order + W.order0 + W.segs[s].first, W.segs[s].count, vec, ghost, work, nrhs_, L.max_n,
                               W.max_ntot, stream_);
            }
            if (L.streamed) launch_mul_diagonal(L.d_boxes, L.d_diag + W.diag0, W.ndiag, vec, nrhs_, L.max_r, stream_);  // (with its wave: H2 8x)
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply forward");
        free_points(work);
        free_points(ghost);
        free_message(sendbox);
        free_message(inbox);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        if (!L.streamed) launch_mul_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal multiplies");
        times.diagonal += seconds_since(t1);
    }
    // A level's backward V waves in reverse (reverse Morton order within a
    // wave), their neighbor updates applied after each as in the solve's
    // forward sweep.
    void mul_backward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* outbox = alloc_message(L.out_max);
        S* inbox = alloc_message(L.in_max);
        for (int w = static_cast<int>(L.waves.size()) - 1; w >= 0; --w) {
            const auto& W = L.waves[static_cast<size_t>(w)];
            prepare_wave_kept(L, static_cast<size_t>(w), -1);
            for (size_t s = W.segs.size(); s-- > 0;) {  // (H2 8x: the wave's segments, in reverse)
                prepare_wave_xnr(L, W.seg0 + s, -1);
                launch_mul_backward(L.d_boxes, L.d_rorder + W.order0 + (W.count - W.segs[s].first - W.segs[s].count), W.segs[s].count, vec, work, nrhs_, L.max_n, stream_);
            }
            launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, outbox, inbox,
                               stream_);
            if (L.peers.empty()) continue;
            exchange(L, W.out, W.in, outbox, inbox, 770);
            launch_solve_accum(L.d_items + W.in_item0, L.d_parts, W.in_nitems, W.in_max_count, nrhs_, vec, work, outbox,
                               inbox, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply backward");
        free_points(work);
        free_message(outbox);
        free_message(inbox);
        times.backward += seconds_since(t0);
    }

    // ---- a replicated CA level (the host CA solve's scheme; the gathers
    //      are the host's gather_CA_boxes_solve, collective over `comm`)

    // The level vector with the local boxes' vectors (the ghosts' part is
    // filled by put_ca_halo).
    S* upload_ca_local(DeviceSolveLevel& L, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t nb = L.boxes.size() - L.ghost_box.size();
        S* d = heap_.alloc<S>(std::max<size_t>(static_cast<size_t>(L.points) * nrhs_, 1) * sizeof(S));
        const size_t len = static_cast<size_t>(L.local_points) * nrhs_;
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        for (size_t b = 0; b < nb; ++b) {
            const auto& v = data[b].left_side;
            if (v.size() != static_cast<size_t>(L.boxes[b].n) * nrhs_) throw std::runtime_error("device solve: vector size");
            std::memcpy(h + L.boxes[b].vec * nrhs_, v.data(), v.size() * sizeof(DataType));
        }
        if (len > 0) check_cuda(cudaMemcpyAsync(d, h, len * sizeof(S), cudaMemcpyHostToDevice, stream_), "solve upload");
        check_cuda(cudaStreamSynchronize(stream_), "solve upload");
        times.transfer += seconds_since(t0);
        return d;
    }
    // The local boxes' vectors back into the solve data (the level vector stays).
    void download_ca_local(DeviceSolveLevel& L, const S* d, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t nb = L.boxes.size() - L.ghost_box.size();
        const size_t len = static_cast<size_t>(L.local_points) * nrhs_;
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        if (len > 0) check_cuda(cudaMemcpyAsync(h, d, len * sizeof(S), cudaMemcpyDeviceToHost, stream_), "solve download");
        check_cuda(cudaStreamSynchronize(stream_), "solve download");
        for (size_t b = 0; b < nb; ++b) {
            auto& v = data[b].left_side;
            std::memcpy(v.data(), h + L.boxes[b].vec * nrhs_, v.size() * sizeof(DataType));
        }
        times.transfer += seconds_since(t0);
    }
    // The host gather's copies of the ghosts (into the level vector) and of
    // the assisting boxes (into the read-only area).
    void put_ca_halo(DeviceSolveLevel& L, TreeLevel<CoordType, DataType>& lvl, S* vec, S* area, bool ghosts,
                     bool assisting) {
        const auto t0 = std::chrono::steady_clock::now();
        auto vector_of = [&](int64_t m, int n) -> const std::vector<DataType>& {
            const auto it = lvl.ghost_and_assisting_box_points_for_solve_map.find(m);
            if (it == lvl.ghost_and_assisting_box_points_for_solve_map.end()) {
                throw std::runtime_error("device solve: no gathered vector for box " + std::to_string(m));
            }
            const auto& v = lvl.ghost_and_assisting_boxes_for_solve[static_cast<size_t>(it->second)].left_side;
            if (v.size() != static_cast<size_t>(n) * nrhs_) throw std::runtime_error("device solve: gathered vector size");
            return v;
        };
        const size_t nb = L.boxes.size() - L.ghost_box.size();
        const size_t glen = static_cast<size_t>(L.points - L.local_points) * nrhs_;
        const size_t alen = static_cast<size_t>(L.ghost_points) * nrhs_;
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(glen + alen, 1) * sizeof(S)));
        if (ghosts && glen > 0) {
            for (size_t j = 0; j < L.ghost_box.size(); ++j) {
                const SolveBox& sb = L.boxes[nb + j];
                const auto& v = vector_of(L.ghost_box[j], sb.n);
                std::memcpy(h + (sb.vec - L.local_points) * nrhs_, v.data(), v.size() * sizeof(DataType));
            }
            check_cuda(cudaMemcpyAsync(vec + static_cast<size_t>(L.local_points) * nrhs_, h, glen * sizeof(S),
                                       cudaMemcpyHostToDevice, stream_),
                       "solve halo upload");
        }
        if (assisting && alen > 0) {
            for (const auto& g : L.ghosts) {
                const auto& v = vector_of(g.morton, g.n);
                std::memcpy(h + glen + g.offset * nrhs_, v.data(), v.size() * sizeof(DataType));
            }
            check_cuda(cudaMemcpyAsync(area, h + glen, alen * sizeof(S), cudaMemcpyHostToDevice, stream_),
                       "solve halo upload");
        }
        check_cuda(cudaStreamSynchronize(stream_), "solve halo upload");
        times.transfer += seconds_since(t0);
    }
    // The halo in device memory (build_ca_halo; kind 0 the ghosts, 1 the
    // ghosts and the assisting boxes, 2 the assisting boxes): the owners'
    // vectors into the level vector's ghost part and the read-only area.
    void ca_halo(DeviceSolveLevel& L, int kind, S* vec, S* area) {
        if (L.peers.empty()) return;
        const auto t0 = std::chrono::steady_clock::now();
        const double comm0 = times.comm;
        const auto& H = L.halo[kind];
        S* sendbox = alloc_message(H.out_total);
        S* inbox = alloc_message(H.in_total);
        launch_solve_copy(L.d_copies + H.send0, H.nsend, H.send_max_n, nrhs_, static_cast<const S*>(vec), sendbox, stream_);
        exchange(L, H.out, H.in, sendbox, inbox, 770 + kind);
        launch_solve_copy(L.d_copies + H.vec0, H.nvec, H.vec_max_n, nrhs_, static_cast<const S*>(inbox), vec, stream_);
        if (area != nullptr) {
            launch_solve_copy(L.d_copies + H.area0, H.narea, H.area_max_n, nrhs_, static_cast<const S*>(inbox), area,
                              stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "CA halo");
        free_message(sendbox);
        free_message(inbox);
        times.halo += seconds_since(t0) - (times.comm - comm0);
    }
    template<typename Gather>
    void timed_gather(Gather&& gather) {
        const auto t0 = std::chrono::steady_clock::now();
        gather();
        times.comm += seconds_since(t0);
        ++times.staged_messages;
    }

    // Solve, forward: the ghosts' and assisting boxes' vectors gathered once,
    // the CA waves over the local boxes and the ghosts, then the local
    // boxes' diagonal solves.
    void ca_forward(DeviceSolveLevel& L, TreeLevel<CoordType, DataType>& lvl, int level,
                    std::vector<SolveDataRequest<CoordType, DataType>>& data, MPI_Comm comm) {
        if (!L.device_halo) timed_gather([&] { gather_CA_boxes_solve(tree_, level, data, comm); });
        S* vec = upload_ca_local(L, data);
        if (L.device_halo) {
            ca_halo(L, 0, vec, nullptr);
        } else {
            put_ca_halo(L, lvl, vec, nullptr, true, false);
        }
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        for (size_t wi = 0; wi < L.waves.size(); ++wi) {
            const auto& W = L.waves[wi];
            if (W.count == 0) continue;
            prepare_wave_kept(L, wi, +1);
            for (size_t s = 0; s < W.segs.size(); ++s) {  // (H2 8x: the wave's segments, each through the scratch)
                prepare_wave_xnr(L, W.seg0 + s, +1);
                launch_solve_forward(L.d_boxes, L.d_order + W.order0 + W.segs[s].first, W.segs[s].count, vec, work, nrhs_, L.max_n, stream_);
            }
            if (L.streamed) launch_solve_diagonal(L.d_boxes, L.d_diag + W.diag0, W.ndiag, vec, nrhs_, L.max_r, stream_);  // (with its wave: H2 8x)
            launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, static_cast<S*>(nullptr),
                               static_cast<const S*>(nullptr), stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "forward sweep");
        free_points(work);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        if (!L.streamed) launch_solve_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal solves");
        times.diagonal += seconds_since(t1);
        download_ca_local(L, vec, data);
        heap_.free(vec);
    }
    // Solve, backward: the groups in reverse; after each but the last (on
    // several ranks), the owners' vectors of the ghosts (the first time) and
    // of the assisting boxes, for the groups still to come.
    void ca_backward(DeviceSolveLevel& L, TreeLevel<CoordType, DataType>& lvl, int level,
                     std::vector<SolveDataRequest<CoordType, DataType>>& data, MPI_Comm comm) {
        S* vec = upload_ca_local(L, data);
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* area = alloc_points(L.ghost_points);
        const int groups = static_cast<int>(L.group_end.size());
        for (int g = groups - 1, done = 0; g >= 0; --g) {
            const int w0 = g == 0 ? 0 : L.group_end[static_cast<size_t>(g - 1)];
            for (int w = L.group_end[static_cast<size_t>(g)] - 1; w >= w0; --w) {
                const auto& W = L.waves[static_cast<size_t>(w)];
                if (W.count == 0) continue;
                prepare_wave_kept(L, static_cast<size_t>(w), -1);
                for (size_t s = W.segs.size(); s-- > 0;) {  // (H2 8x: the wave's segments, in reverse)
                    prepare_wave_xnr(L, W.seg0 + s, -1);
                    launch_solve_backward(L.d_boxes, L.d_slots, L.d_rorder + W.order0 + (W.count - W.segs[s].first - W.segs[s].count), W.segs[s].count, vec, area, work, nrhs_,
                                      L.max_n, W.max_ntot, stream_);
                }
            }
            if (lvl.num_active_processes > 1 && ++done < groups) {
                if (L.device_halo) {
                    ca_halo(L, done == 1 ? 1 : 2, vec, area);
                } else {
                    download_ca_local(L, vec, data);
                    timed_gather([&] { gather_CA_boxes_solve(tree_, level, data, comm, /*assist_only=*/done > 1); });
                    put_ca_halo(L, lvl, vec, area, /*ghosts=*/done == 1, /*assisting=*/true);
                }
            }
        }
        check_cuda(cudaStreamSynchronize(stream_), "backward sweep");
        free_points(work);
        free_points(area);
        times.backward += seconds_since(t0);
        download_ca_local(L, vec, data);
        heap_.free(vec);
    }
    // Multiply, forward W: everything gathered once, the CA waves in order
    // (reading neighbors), then the local boxes' diagonal multiplies.
    void ca_mul_forward(DeviceSolveLevel& L, TreeLevel<CoordType, DataType>& lvl, int level,
                        std::vector<SolveDataRequest<CoordType, DataType>>& data, MPI_Comm comm) {
        if (!L.device_halo) timed_gather([&] { gather_CA_boxes_solve(tree_, level, data, comm); });
        S* vec = upload_ca_local(L, data);
        S* area = alloc_points(L.ghost_points);
        if (L.device_halo) {
            ca_halo(L, 1, vec, area);
        } else {
            put_ca_halo(L, lvl, vec, area, true, true);
        }
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        for (size_t wi = 0; wi < L.waves.size(); ++wi) {
            const auto& W = L.waves[wi];
            if (W.count == 0) continue;
            prepare_wave_kept(L, wi, +1);
            for (size_t s = 0; s < W.segs.size(); ++s) {  // (H2 8x: the wave's segments, each through the scratch)
                prepare_wave_xnr(L, W.seg0 + s, +1);
                launch_mul_forward(L.d_boxes, L.d_slots, L.d_order + W.order0 + W.segs[s].first, W.segs[s].count, vec, area, work, nrhs_, L.max_n,
                               W.max_ntot, stream_);
            }
            if (L.streamed) launch_mul_diagonal(L.d_boxes, L.d_diag + W.diag0, W.ndiag, vec, nrhs_, L.max_r, stream_);  // (with its wave: H2 8x)
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply forward");
        free_points(work);
        free_points(area);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        if (!L.streamed) launch_mul_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal multiplies");
        times.diagonal += seconds_since(t1);
        download_ca_local(L, vec, data);
        heap_.free(vec);
    }
    // Multiply, backward V: everything gathered, the groups in reverse with
    // updates to local or ghost boxes, everything gathered again after each
    // group but the last (on several ranks).
    void ca_mul_backward(DeviceSolveLevel& L, TreeLevel<CoordType, DataType>& lvl, int level,
                         std::vector<SolveDataRequest<CoordType, DataType>>& data, MPI_Comm comm) {
        if (!L.device_halo) timed_gather([&] { gather_CA_boxes_solve(tree_, level, data, comm); });
        S* vec = upload_ca_local(L, data);
        S* area = alloc_points(L.ghost_points);
        if (L.device_halo) {
            ca_halo(L, 1, vec, area);
        } else {
            put_ca_halo(L, lvl, vec, area, true, true);
        }
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        const int groups = static_cast<int>(L.group_end.size());
        for (int g = groups - 1, done = 0; g >= 0; --g) {
            const int w0 = g == 0 ? 0 : L.group_end[static_cast<size_t>(g - 1)];
            for (int w = L.group_end[static_cast<size_t>(g)] - 1; w >= w0; --w) {
                const auto& W = L.waves[static_cast<size_t>(w)];
                if (W.count == 0) continue;
                prepare_wave_kept(L, static_cast<size_t>(w), -1);
                for (size_t s = W.segs.size(); s-- > 0;) {  // (H2 8x: the wave's segments, in reverse)
                    prepare_wave_xnr(L, W.seg0 + s, -1);
                    launch_mul_backward(L.d_boxes, L.d_rorder + W.order0 + (W.count - W.segs[s].first - W.segs[s].count), W.segs[s].count, vec, work, nrhs_, L.max_n, stream_);
                }
                launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, static_cast<S*>(nullptr),
                                   static_cast<const S*>(nullptr), stream_);
            }
            if (lvl.num_active_processes > 1 && ++done < groups) {
                if (L.device_halo) {
                    ca_halo(L, 1, vec, area);
                } else {
                    download_ca_local(L, vec, data);
                    timed_gather([&] { gather_CA_boxes_solve(tree_, level, data, comm); });
                    put_ca_halo(L, lvl, vec, area, true, true);
                }
            }
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply backward");
        free_points(work);
        free_points(area);
        times.backward += seconds_since(t0);
        download_ca_local(L, vec, data);
        heap_.free(vec);
    }

    DeviceSolveTimes times;

private:
    static double seconds_since(std::chrono::steady_clock::time_point t) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    }
    S* alloc_points(int64_t points) {
        return points > 0 ? heap_.alloc<S>(static_cast<size_t>(points) * nrhs_ * sizeof(S)) : nullptr;
    }
    void free_points(S* p) {
        if (p != nullptr) heap_.free(p);
    }
    // Message buffers: in the MPI exchange arena when MPI is CUDA-aware
    // (sent and received in device memory), else staged through the host.
    S* alloc_message(int64_t points) {
        if (points <= 0) return nullptr;
        const size_t bytes = static_cast<size_t>(points) * nrhs_ * sizeof(S);
        if (device_exchange_enabled()) {
            if (char* p = DeviceHeap::exchange_arena().try_alloc(bytes)) return reinterpret_cast<S*>(p);
        }
        return heap_.alloc<S>(bytes);
    }
    void free_message(S* p) {
        if (p == nullptr) return;
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (arena.owns(p)) {
            arena.free(p);
        } else {
            heap_.free(p);
        }
    }
    // Other ranks' vectors that changed, into the ghost area.
    void refresh(DeviceSolveLevel& L, const DeviceSolveLevel::Wave::Refresh& R, S* vec, S* ghost,
                 S* sendbox, S* inbox) {
        launch_solve_copy(L.d_copies + R.send0, R.nsend, R.send_max_n, nrhs_, vec, sendbox, stream_);
        exchange(L, R.send, R.recv, sendbox, inbox, 760);
        launch_solve_copy(L.d_copies + R.recv0, R.nrecv, R.recv_max_n, nrhs_, inbox, ghost, stream_);
    }

    // Messages of known size with the neighbor ranks: out[q] and in[q] points
    // to and from peer q, in device memory when both buffers are in the MPI
    // exchange arena, else staged through the host.
    void exchange(DeviceSolveLevel& L, const std::vector<int64_t>& out, const std::vector<int64_t>& in,
                  const S* d_out, S* d_in, int tag) {
        const size_t np = L.peers.size();
        int64_t total_out = 0, total_in = 0;
        for (size_t q = 0; q < np; ++q) {
            total_out += out[q];
            total_in += in[q];
        }
        check_cuda(cudaStreamSynchronize(stream_), "solve messages");
        const auto t0 = std::chrono::steady_clock::now();
        const DeviceHeap& arena = DeviceHeap::exchange_arena();
        if ((total_out == 0 || arena.owns(d_out)) && (total_in == 0 || arena.owns(d_in))) {
            // CUDA-aware MPI: straight from and into device memory
            std::vector<MPI_Request> reqs;
            int64_t at = 0;
            for (size_t q = 0; q < np; ++q) {
                if (in[q] > 0) {
                    reqs.emplace_back();
                    MPI_Irecv(d_in + at * nrhs_, static_cast<int>(in[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                              &reqs.back());
                }
                at += in[q];
            }
            at = 0;
            for (size_t q = 0; q < np; ++q) {
                if (out[q] > 0) {
                    reqs.emplace_back();
                    MPI_Isend(d_out + at * nrhs_, static_cast<int>(out[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag,
                              tree_->comm, &reqs.back());
                }
                at += out[q];
            }
            MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
            if (total_out + total_in > 0) ++times.direct_messages;
            times.comm += seconds_since(t0);
            return;
        }
        if (total_out + total_in > 0) ++times.staged_messages;
        S* h_out = static_cast<S*>(host_out_.reserve(std::max<int64_t>(total_out, 1) * nrhs_ * sizeof(S)));
        S* h_in = static_cast<S*>(host_in_.reserve(std::max<int64_t>(total_in, 1) * nrhs_ * sizeof(S)));
        if (total_out > 0) {
            check_cuda(cudaMemcpy(h_out, d_out, static_cast<size_t>(total_out) * nrhs_ * sizeof(S),
                                  cudaMemcpyDeviceToHost), "solve messages");
        }
        std::vector<MPI_Request> reqs;
        int64_t at = 0;
        for (size_t q = 0; q < np; ++q) {
            if (in[q] > 0) {
                reqs.emplace_back();
                MPI_Irecv(h_in + at * nrhs_, static_cast<int>(in[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                          &reqs.back());
            }
            at += in[q];
        }
        at = 0;
        for (size_t q = 0; q < np; ++q) {
            if (out[q] > 0) {
                reqs.emplace_back();
                MPI_Isend(h_out + at * nrhs_, static_cast<int>(out[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                          &reqs.back());
            }
            at += out[q];
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        if (total_in > 0) {
            check_cuda(cudaMemcpyAsync(d_in, h_in, static_cast<size_t>(total_in) * nrhs_ * sizeof(S),
                                       cudaMemcpyHostToDevice, stream_), "solve messages");
        }
        times.comm += seconds_since(t0);
    }

    ParallelTree<CoordType, DataType>* tree_;
    int nrhs_;
    DeviceHeap& heap_ = DeviceHeap::instance();
    cudaStream_t stream_ = Context::instance().stream();
    PinnedBuffer& host_vec_ = device_solve_store().host_vec;
    PinnedBuffer& host_out_ = device_solve_store().host_out;
    PinnedBuffer& host_in_ = device_solve_store().host_in;
};

// The steps of a CA level run on the host (DeviceSolveLevel::host): the
// host CA solve's code for that level on its solve data, as the callback
// `host_level(level, step)` of run_device_solve.
enum class HostStep { SolveForward, SolveBackward, MulForward, MulBackward };
using HostLevelStep = std::function<void(int, HostStep)>;

// A replicated CA level's communicator for the host gathers.
inline MPI_Comm ca_comm(const std::vector<MPI_Comm>* level_comms, int level) {
    if (level_comms == nullptr || static_cast<size_t>(level) >= level_comms->size()) {
        throw std::runtime_error("device solve: a CA level without its communicator");
    }
    return (*level_comms)[static_cast<size_t>(level)];
}

// The solve sweeps on the device (solve data initialized as by the host
// solve).  Collective over the tree's ranks.
template<typename CoordType, typename DataType>
void device_solve_sweeps(ParallelTree<CoordType, DataType>* tree,
                         std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data, int nrhs,
                         int verbosity, const std::vector<MPI_Comm>* level_comms, const HostLevelStep* host_level) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    DeviceSolveStore& store = device_solve_store();
    const int leaf = tree->num_levels - 1;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    DeviceSolveRun<CoordType, DataType> run(tree, nrhs);
    double host_seconds = 0.0;
    store.stream_seconds = store.stream_gb = 0.0;
    // forward (with each level's diagonal solves)
    for (int level = leaf; level >= 1; --level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (L.host) {
                const auto tl = clock::now();
                (*host_level)(level, HostStep::SolveForward);
                host_seconds += std::chrono::duration<double>(clock::now() - tl).count();
            } else if (L.ca) {
                run.ca_forward(L, lvl, level, solve_data[static_cast<size_t>(level)], ca_comm(level_comms, level));
            } else {
                auto* vec = run.upload(L, solve_data[static_cast<size_t>(level)]);
                run.forward(L, vec);
                run.download(L, vec, solve_data[static_cast<size_t>(level)]);
            }
        }
        const auto th = clock::now();
        gather_skeleton_to_parent(lvl, tree->levels[static_cast<size_t>(level - 1)], solve_data[static_cast<size_t>(level)],
                                  solve_data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    }
    const auto th = clock::now();
    auto& root = tree->levels[0];
    if (root.is_process_active && !root.local_boxes.empty()) apply_diagonal_solve(root, solve_data[0][0], false);
    host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    // backward
    for (int level = 1; level <= leaf; ++level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        auto& parent = tree->levels[static_cast<size_t>(level - 1)];
        const auto ts = clock::now();
        if (lvl.is_process_active || parent.is_process_active) {
            scatter_solution_to_children(lvl, parent, solve_data[static_cast<size_t>(level)],
                                         solve_data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        }
        host_seconds += std::chrono::duration<double>(clock::now() - ts).count();
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (L.host) {
                const auto tl = clock::now();
                (*host_level)(level, HostStep::SolveBackward);
                host_seconds += std::chrono::duration<double>(clock::now() - tl).count();
            } else if (L.ca) {
                run.ca_backward(L, lvl, level, solve_data[static_cast<size_t>(level)], ca_comm(level_comms, level));
            } else {
                auto* vec = run.upload(L, solve_data[static_cast<size_t>(level)]);
                run.backward(L, vec);
                run.download(L, vec, solve_data[static_cast<size_t>(level)]);
            }
        }
    }
    const double total = std::chrono::duration<double>(clock::now() - t0).count();
    run.xnr_timers_drain();
    if (verbosity >= 0 && rank == 0) {
        const auto& t = run.times;
        std::printf("GPU solve: %.3f s (rank 0: forward %.3f, diagonal %.3f, backward %.3f, of which MPI %.3f%s; "
                    "vector transfers %.3f; CA halo copies %.3f; host level hand-off and root %.3f)\n",
                    total, t.forward, t.diagonal, t.backward, t.comm,
                    t.direct_messages > 0 ? (t.staged_messages > 0 ? " partly in device memory" : " in device memory")
                                          : (t.staged_messages > 0 ? " through the host" : ""),
                    t.transfer, t.halo, host_seconds);
        if (t.xnr_calls > 0) {
            std::printf("  X_NR decompression: %.3f s in %lld calls, %.2f GB uncompressed, %.1f GB/s\n", t.xnr,
                        static_cast<long long>(t.xnr_calls), t.xnr_bytes / 1e9,
                        t.xnr > 0.0 ? t.xnr_bytes / t.xnr / 1e9 : 0.0);
        }
        if (store.stream_levels > 0) {
            std::printf("  streamed kept factors: %d levels, %.2f GB through a %.2f GB ring (host-side %.3f s)\n",
                        store.stream_levels, store.stream_gb, static_cast<double>(store.stream_bytes) / 1e9,
                        store.stream_seconds);
        }
        std::fflush(stdout);
    }
}

// The multiply F x on the device (data initialized as by the host multiply;
// hierarchical_mul_parallel).  Collective over the tree's ranks.
template<typename CoordType, typename DataType>
void device_mul_sweeps(ParallelTree<CoordType, DataType>* tree,
                       std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& data, int nrhs, bool verbose,
                       const std::vector<MPI_Comm>* level_comms, const HostLevelStep* host_level) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    DeviceSolveStore& store = device_solve_store();
    const int leaf = tree->num_levels - 1;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    DeviceSolveRun<CoordType, DataType> run(tree, nrhs);
    double host_seconds = 0.0;
    store.stream_seconds = store.stream_gb = 0.0;
    // forward W (with each level's diagonal multiplies)
    for (int level = leaf; level >= 1; --level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (L.host) {
                const auto tl = clock::now();
                (*host_level)(level, HostStep::MulForward);
                host_seconds += std::chrono::duration<double>(clock::now() - tl).count();
            } else if (L.ca) {
                run.ca_mul_forward(L, lvl, level, data[static_cast<size_t>(level)], ca_comm(level_comms, level));
            } else {
                auto* vec = run.upload(L, data[static_cast<size_t>(level)]);
                run.mul_forward(L, vec);
                run.download(L, vec, data[static_cast<size_t>(level)]);
            }
        }
        const auto th = clock::now();
        gather_skeleton_to_parent(lvl, tree->levels[static_cast<size_t>(level - 1)], data[static_cast<size_t>(level)],
                                  data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    }
    const auto th = clock::now();
    auto& root = tree->levels[0];
    if (root.is_process_active && !root.local_boxes.empty()) apply_diagonal_multiply(root, data[0][0], false);
    host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    // backward V
    for (int level = 1; level <= leaf; ++level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        auto& parent = tree->levels[static_cast<size_t>(level - 1)];
        const auto ts = clock::now();
        if (lvl.is_process_active || parent.is_process_active) {
            scatter_solution_to_children(lvl, parent, data[static_cast<size_t>(level)],
                                         data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        }
        host_seconds += std::chrono::duration<double>(clock::now() - ts).count();
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (L.host) {
                const auto tl = clock::now();
                (*host_level)(level, HostStep::MulBackward);
                host_seconds += std::chrono::duration<double>(clock::now() - tl).count();
            } else if (L.ca) {
                run.ca_mul_backward(L, lvl, level, data[static_cast<size_t>(level)], ca_comm(level_comms, level));
            } else {
                auto* vec = run.upload(L, data[static_cast<size_t>(level)]);
                run.mul_backward(L, vec);
                run.download(L, vec, data[static_cast<size_t>(level)]);
            }
        }
    }
    const double total = std::chrono::duration<double>(clock::now() - t0).count();
    run.xnr_timers_drain();
    if (verbose && rank == 0) {
        const auto& t = run.times;
        std::printf("GPU multiply: %.3f s (rank 0: forward %.3f, diagonal %.3f, backward %.3f, of which MPI %.3f%s; "
                    "vector transfers %.3f; CA halo copies %.3f; host level hand-off and root %.3f)\n",
                    total, t.forward, t.diagonal, t.backward, t.comm,
                    t.direct_messages > 0 ? (t.staged_messages > 0 ? " partly in device memory" : " in device memory")
                                          : (t.staged_messages > 0 ? " through the host" : ""),
                    t.transfer, t.halo, host_seconds);
        if (t.xnr_calls > 0) {
            std::printf("  X_NR decompression: %.3f s in %lld calls, %.2f GB uncompressed, %.1f GB/s\n", t.xnr,
                        static_cast<long long>(t.xnr_calls), t.xnr_bytes / 1e9,
                        t.xnr > 0.0 ? t.xnr_bytes / t.xnr / 1e9 : 0.0);
        }
        if (store.stream_levels > 0) {
            std::printf("  streamed kept factors: %d levels, %.2f GB through a %.2f GB ring (host-side %.3f s)\n",
                        store.stream_levels, store.stream_gb, static_cast<double>(store.stream_bytes) / 1e9,
                        store.stream_seconds);
        }
        std::fflush(stdout);
    }
}

// The solve (multiply = false) or the multiply F x of the host solvers
// (full-grid and unstructured) on the device, when the device solve runs
// (decided collectively): solve_data initialized as the host path does, the
// result left in its leaf left_side.  With BPACK_CHECK=solve the host
// path runs too (host(host_data), with the device path suspended), and rank
// 0 prints their difference.  Returns false when the host path must run.
// `level_comms`: the solve's per-level communicators (the host CA solve's
// gathers on replicated CA levels; null where there are none).
// `host_level`: a CA level's step on the host, on solve_data (for the CA
// levels whose factors do not fit in device memory; null: none can run so).
template<typename CoordType, typename DataType, typename HostRun>
bool run_device_solve(ParallelTree<CoordType, DataType>* tree,
                      std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data, int nrhs,
                      bool multiply, int verbosity, HostRun&& host, const std::vector<MPI_Comm>* level_comms = nullptr,
                      const HostLevelStep* host_level = nullptr) {
    activate_operator(tree);  // (before the store is read: one per operator)
    // (a nonsymmetric factor: the host solve, whose backward sweep reads X_RS
    // and X_RN; the device kernels read X_SR^T and X_NR^T)
    if (device_solve_suspended() || !color_gpu_enabled()) {
        materialize_host_factors(tree);  // (the host path reads them)
        return false;
    }
    if (!prepare_device_solve(tree, verbosity, host_level != nullptr)) return false;
    if (multiply) {
        device_mul_sweeps(tree, solve_data, nrhs, verbosity >= 1, level_comms, host_level);
    } else {
        device_solve_sweeps(tree, solve_data, nrhs, verbosity, level_comms, host_level);
    }
    if (device_solve_check()) {
        const int leaf = tree->num_levels - 1;
        std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> host_data(
            static_cast<size_t>(tree->num_levels));
        device_solve_suspended() = true;
        host(host_data);
        device_solve_suspended() = false;
        double sums[2] = {0.0, 0.0};
        if (tree->levels[static_cast<size_t>(leaf)].is_process_active) {
            for (size_t b = 0; b < solve_data[static_cast<size_t>(leaf)].size(); ++b) {
                const auto& x = solve_data[static_cast<size_t>(leaf)][b].left_side;
                const auto& y = host_data[static_cast<size_t>(leaf)][b].left_side;
                for (size_t i = 0; i < x.size(); ++i) {
                    sums[0] += std::norm(x[i] - y[i]);
                    sums[1] += std::norm(y[i]);
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, tree->comm);
        int rank = 0;
        MPI_Comm_rank(tree->comm, &rank);
        if (rank == 0) {
            std::printf(multiply ? "GPU multiply check: |y_gpu - y_host| / |y_host| = %.3e\n"
                                 : "GPU solve check: |x_gpu - x_host| / |x_host| = %.3e\n",
                        sums[1] > 0.0 ? std::sqrt(sums[0] / sums[1]) : 0.0);
            std::fflush(stdout);
        }
    }
    return true;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
