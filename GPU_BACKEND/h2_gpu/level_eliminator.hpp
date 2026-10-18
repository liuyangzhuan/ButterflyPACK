#pragma once
// Device-resident elimination of the Color waves of one level (H2 GPU backend).
//
// For the whole level the GPU holds the matrix blocks that later waves modify:
// the near-field block of every one-hop pair (one physical matrix per edge,
// rows = the higher Morton box, columns = the lower one) and the self (Schur)
// block of every box.  Blocks always have their current shape: a box's rows
// or columns are its full point set until it is eliminated and its skeleton
// afterwards, so a block never needs slicing when it is read.
//
// Per wave, the host computes the sketch and the ID of every box (the CPU
// code), plans all device work, uploads one metadata image, and the device
// runs the elimination batched over the wave's boxes: X_RR, X_SR, X_NR and
// the LU of X_RR, temp1 / temp2, the Schur block, the near-block update of
// step 5, and the deferred X_NN owner pass.  The host then receives the
// factors each box keeps for the solve and for the sketches of later waves
// (LU and pivots, temp1, temp2, X_RR_full, X_RS_entry).  At level end the
// blocks come back to the host so the CPU level transition runs unchanged.
//
// Scope (checked by level_eliminator_supported): real or complex symmetric
// kernels with a registered device form, streamed sketches with lazy far
// fill, owner-deferred X_NN updates, LU of X_RR.  Complex data is symmetric,
// not Hermitian: every transpose is plain.  The arithmetic matches the CPU
// path up to rounding (GEMM summation order differs).
//
// Replicated CA levels (`ca`, h2_gpu/GPU_CA_PLAN.md M1): the ghost boxes
// the halo brought in are eliminated here as well (writable boxes: point
// slots, Schur and near blocks on the device; their blocks from the host
// halo gather, or, after another device level, adopted from the device halo
// of h2_gpu/ca_halo.hpp, M2), in the CA groups' waves; nothing is
// transported during the level (no generators);
// a pair of writable boxes is updated once, into its single block.  Every
// copy of a box, on its owner and on the ranks holding it as a ghost, must
// come out bitwise identical: the device ID is (every launch variant gives
// the same bits), and the X_RR LU and its solves run per size class
// (tests/batch_determinism.cpp).

#ifdef H2_HAVE_GPU

#include "adaptive_rows.hpp"
#include "device_heap.hpp"
#include "factor_compress.hpp"
#include "device_kernels.hpp"
#include "device_solve.hpp"
#include "solve_kernels.hpp"
#include "gpu_runtime.hpp"
#include "host_copier.hpp"
#include "evaluator.hpp"
#include "owner_plan.hpp"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <deque>
#include <map>
#include <exception>
#include <cmath>
#include <iomanip>
#include <limits>
#include <list>
#include <random>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>
#include <optional>

namespace fmm {
namespace gpu {

// Per-level totals of the device box path (seconds and bytes, this rank).
// The fill store's host copies (H2 host staging): one arena for the process,
// reset at each level's end.
inline PinnedArena& fill_host_arena() {
    static PinnedArena* const arena = new PinnedArena;  // (lives to the process's end, as the heap does)
    return *arena;
}

// The near store's host copies (H2 host staging (c)): two arenas by level
// parity, because a level's transition packs the parent blocks for the
// level below while its own store is still read.  A level's region holds
// its store and the runs it adopted from the level above, and resets when
// its eliminator goes (release_near_store).
inline PinnedArena& near_host_arena(int level) {
    static PinnedArena* const arenas[2] = {new PinnedArena, new PinnedArena};  // (live to the process's end)
    return *arenas[level & 1];
}
// A level's near store compacted at its end (compact_near_store): its live
// blocks, out of the store's region, until the level's eliminator goes.
inline PinnedArena& near_compact_arena() {
    static PinnedArena* const arena = new PinnedArena;  // (lives to the process's end)
    return *arena;
}
// The eliminator's hold limit (memory pressure), published for the driver.
inline size_t g_hold_limit_bytes = 0;
inline bool device_under_pressure() {
    const DeviceHeap& h = DeviceHeap::instance();
    return g_hold_limit_bytes > 0 && h.used() > g_hold_limit_bytes;
}

struct EliminatorStats {
    double begin = 0.0;       // point table and block upload
    double sketch = 0.0;      // device sketch (plan, launches, Y download)
    double sketch_plan = 0.0; //   host index lists
    double sketch_gpu = 0.0;  //   metadata build and upload, kernels, Y download
    double sk_wait = 0.0;     // host blocked on the sketch (and the previous wave's device work)
    double sk_meta = 0.0, sk_upload = 0.0, sk_rows = 0.0, sk_stored = 0.0, sk_p = 0.0, sk_fill = 0.0,
           sk_id = 0.0, sk_download = 0.0;  //   ... its parts (device times from events)
    double el[8] = {0, 0, 0, 0, 0, 0, 0, 0};  // device: fills, X_RR/X_SR, LU, X_NR, solves, Schur+step 5, owner targets, owner GEMMs
    double finish_sources = 0.0;  // level end: wait for the background copies
    double finish_store = 0.0;    // background copier: busy time (overlapped)
    double copier_wait = 0.0;     //   ... of which waiting for device chunks
    double sk_adaptive = 0.0;     // adaptive ID rows (H2_ID_proxy 2): the rounds after the first (host time)
    int64_t adaptive_rounds = 0;  //   ... rounds of all waves; boxes that took more than one node
    int64_t adaptive_multi = 0;
    double id = 0.0;          // host ID (and host sketch without the device sketch)
    double id_loop = 0.0;     //   ... device ID: its per-box store loop
    double plan = 0.0;        // host planning and metadata
    double plan_boxes = 0.0, plan_owner = 0.0;  //   ... box region, owner pass
    // their parts: buffers (layouts and allocations), launch items (of which
    // the near blocks' allocations); owner candidates, pairs (parallel, with
    // a CA level's order by source), the join (new targets), batches
    double pb_buffers = 0.0, pb_items = 0.0;
    // BPACK_TRACE=wave: the heap's own host time in the box region's buffers,
    // and the box region's planning while the stream had nothing queued (the
    // GPU waits for the host there)
    double pb_heap = 0.0, pb_idle_plan = 0.0;
    int64_t pb_heap_ops = 0, pb_idle_chunks = 0, pb_chunks = 0;
    // the owner groups' raw blocks in one heap block a group (group_arena_):
    // groups, blocks, groups with no range for it, blocks moved out raw
    int64_t og_arena_groups = 0, og_arena_blocks = 0, og_arena_declined = 0, og_arena_moved = 0;
    // the owner pass's host time by part of the group loop
    double og_mat = 0.0, og_win = 0.0, og_stage = 0.0, og_launch = 0.0, og_pack = 0.0;
    double pb_near = 0.0, pb_near_alloc = 0.0;  // step 5's near blocks (after the box region's launch), their allocations
    double po_candidates = 0.0, po_pairs = 0.0, po_order = 0.0, po_batches = 0.0;
    // The fill sources kept for later waves' sketches (keep_fill_source: per
    // box X_RR_full and temp2^T, the second one ntot x r like X_NR): what
    // they hold now and the most they held, which is the level's resident
    // peak less the kept factors and the tables.
    double persist_bytes = 0.0;
    double persist_peak = 0.0;
    double persist_xrr_bytes = 0.0;   // ... of which X_RR_full
    double persist_gen_bytes = 0.0;   // ... of which the generator boxes' (other ranks' assisting boxes)
    int64_t xr_dropped = 0;           // boxes whose X_RR_full is rebuilt from their LU instead of stored
    int64_t xr_rebuilt = 0;           // ... and the rebuilds that cost
    double ghost_lu_bytes = 0.0;      // H2_CA_solve home: the ghosts' LU, held to the level's end
    // H2 host memory (2026-10-10): the fill sources freed without a host
    // copy made, the X_RR_full the transition took from the host for them,
    // and the X_RR_full rebuilt from the LU when a keep failed
    int64_t rc_sources = 0, tr_xrr_host = 0, xr_restored = 0;
    // apply_remote_generators: the partial products of other ranks' generators,
    // all staged before they are summed into their targets
    double gw_bytes = 0.0;            // the staging that costs
    double gw_target_bytes = 0.0;     // the distinct targets they sum into
    int64_t gw_products = 0;
    int64_t gw_targets = 0;
    int64_t gw_max_parts = 0;         // the most partials one target takes
    int64_t persist_blocks = 0;
    double reclaim_wait = 0.0;    // full heap: waits for background copies to free blocks
    double exchange_wait = 0.0;   // full exchange arena: waits for receive buffers' host copies
    double launch = 0.0;                        // host time issuing a wave's device work
    double device = 0.0;      // device elimination and owner pass (device time)
    double download = 0.0;    // factor download
    double store = 0.0;       // host BoxData update
    double finish = 0.0;      // level end (wait for the background copies)
    double transition = 0.0;  // device transition to the parent level
    double tr_plan = 0.0, tr_blocks = 0.0, tr_p = 0.0, tr_fill = 0.0;  // its parts
    double compaction = 0.0, compaction_bytes = 0.0;  // level end: kept blocks slid down the heap
    double reduction = 0.0;   // device process reduction (the parents' blocks to their owner)
    double tr_restore = 0.0, tr_restore_bytes = 0.0;  // restored fill sources (host gather, in plan)
    int64_t tr_chunks = 0;
    // Per transition chunk, the blocks its parents' pairs write into: the
    // largest chunk's bytes, every chunk's summed, and how many times a block
    // is written by a chunk after the first (an edge between parents in
    // different chunks).  If the sum is close to the total and the revisits
    // few, a chunk's blocks could be decompressed, updated and recompressed
    // for the chunk only, and the rest of the level's would stay compressed.
    // Per sketch wave, the fill sources its boxes read: the largest wave's
    // distinct temp2^T bytes (the scratch a per-wave decompression would
    // need), the sum over waves (what it would decompress in all) and the
    // (wave, source) pairs against the distinct sources (how many waves read
    // a source, so how many times it would be decompressed).
    // ... and the slices the readers actually take: a box reads a source's
    // temp2^T only over its own columns, which is a contiguous byte range, so
    // a decompression could cover the chunks of that range alone.
    // H2_compress_factors: the fill sources' temp2^T, compressed one block per
    // slot (a slot is one one-hop neighbor's rows, which is what a reader
    // takes), measured as the compression runs.  BPACK_CHECK=fillsrc also
    // decompresses each wave's and reports the worst error against the
    // original.
    double fs_raw_bytes = 0.0;
    double fs_stored_bytes = 0.0;
    int64_t fs_slots = 0;
    int64_t fs_boxes = 0;
    // H2_compress_factors: the kept X_NR made from the fill source's temp2^T
    // rather than from a second compression of temp2 (transpose_kept_xnr)
    int64_t kx_boxes = 0;          // handed over
    int64_t kx_groups = 0;
    int64_t kx_host = 0;           // boxes that had to go round through the host
    double kx_raw_bytes = 0.0;     // uncompressed bytes transposed
    double kx_stored_bytes = 0.0;  // what their kept X_NR came to
    double kx_seconds = 0.0;
    // BPACK_CHECK=fillsrc: what one owner group distinctly needs of the wave's
    // two n_tot x r outputs, which is what a bounded window would hold
    double os_union_max = 0.0;
    double os_union_sum = 0.0;
    double os_raw_max = 0.0;      // the same group's slot bytes with repeats
    int64_t os_slots = 0;
    int64_t os_groups = 0;
    // the X_NR^T store: the owner pass's other operand, per wave
    double xt_raw_bytes = 0.0;
    double xt_stored_bytes = 0.0;
    int64_t xt_slots = 0;
    int64_t xt_boxes = 0;
    int64_t xt_short = 0;      // groups the staging did not fit
    // the owner pass's source window
    double ow_bytes = 0.0;     // decompressed into it
    double ow_max = 0.0;       // the largest group's window
    int64_t ow_groups = 0;
    int64_t ow_slots = 0;
    int64_t fs_waves_plain = 0;   // a wave whose staging did not fit
    int64_t fs_packed_short = 0;  // a box whose packed block did not fit (it stays raw)
    int64_t np_shared_blocks = 0;  // replicated CA levels: near blocks on more than one rank, held packed always
    double np_shared_raw = 0.0;    // ... their raw bytes
    // the scratch a wave's sketch needs for the slots it reads, and what it
    // decompressed in all
    double fs_host_bytes = 0.0;    // the packed temp2^T on the host, at most (H2 host staging)
    double ns_host_bytes = 0.0;    // the near store on the host, at most
    double ns_compact_live = 0.0, ns_compact_held = 0.0, ns_compact_s = 0.0;  // compacted at the level end
    int64_t fs_read_dev[3] = {0, 0, 0}, fs_read_host[3] = {0, 0, 0};  // slot chunks read by the window, the sketch, the transition: from the device / the host
    double fs_wave_scratch = 0.0;
    double fs_wave_bytes = 0.0;
    int64_t fs_wave_slots = 0;
    int64_t fs_wave_groups = 0;
    // the transition putting compressed fill sources back whole
    // What a wave reads of the blocks the previous transition left: a box's own
    // Schur block, read when it is eliminated, and the edge blocks of its
    // one-hop neighbors.  The union per wave is the scratch a per-wave
    // decompression would need; the pairs against the distinct blocks say how
    // many waves read a block, so how often it would be decompressed.
    // H2_compress_factors: what the parent blocks this transition built would
    // come to compressed.  Measured over all of them at the end of
    // build_parent, where they are complete; they are still stored raw.
    // the adopted blocks the waves materialized, and what it cost
    double mat_bytes = 0.0;
    double mat_seconds = 0.0;
    int64_t mat_blocks = 0, mat_batches = 0, mat_left = 0;
    int64_t mat_prefetched = 0;   // group gathers run under the previous group's GEMMs
    double mat_site_bytes[6] = {0, 0, 0, 0, 0, 0};   // which caller decompressed what
    int64_t mat_site_calls[6] = {0, 0, 0, 0, 0, 0};
    int64_t owner_groups = 0;   // the owner pass's bounded target groups
    // the near blocks held packed between groups
    double np_raw = 0.0, np_stored = 0.0, np_seconds = 0.0;
    int64_t np_blocks = 0, np_groups = 0, np_groups_plain = 0, np_reclaimed = 0;
    int64_t np_tails = 0, np_no_room = 0;   // chunks that grew; blocks left uncompressed
    // ... of which the adopted blocks the remote-generator passes put back,
    // re-packed right after each pass rather than left raw for a later reader
    int64_t np_gen_blocks = 0;
    double np_gen_raw = 0.0;
    // ... and the targets those passes made from the kernel, packed with them
    int64_t np_gen_fresh_blocks = 0;
    double np_gen_fresh_raw = 0.0;
    // the owner groups that packed under pressure, and those left raw
    int64_t np_pressure_groups = 0, np_raw_groups = 0;
    int64_t np_waves = 0, np_waves_chunked = 0;      // waves whose production ran in chunks (for room)
    int64_t tr_pack_chunks = 0, tr_pack_chunks_packed = 0;   // transition chunks, and those whose parents were packed (under pressure)
    int64_t fs_staged_late = 0;                      // fill sources copied to the host on demand
    int64_t fw_no_window = 0, fw_fallbacks = 0;      // waves without a gather window; blocks put back one by one instead
    // ... of which step 5's fresh blocks, packed as soon as add_store has
    // written them rather than waiting for an owner group to come by
    double np5_raw = 0.0, np5_stored = 0.0, np5_seconds = 0.0;
    int64_t np5_blocks = 0, np5_groups = 0;
    // the wave's gather window: what it put back, in how many passes
    double fw_bytes = 0.0, fw_seconds = 0.0;
    int64_t fw_batches = 0, fw_passes = 0, fw_peak = 0;
    // the compressed length trajectory of the near-field blocks
    double nl_raw = 0.0, nl_first = 0.0, nl_max = 0.0, nl_last = 0.0, nl_worst = 0.0;
    double nl_bin[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int64_t nl_blocks = 0;
    // What a wave writes among the near-field blocks: the working set a
    // decompress-update-recompress scheme would hold, and the traffic it
    // would pay.  A window is one eliminate_wave and the remote generators
    // applied before it.
    // The wave's transient, by array, taken at the wave whose total is the
    // largest: which part of it tiling nt could actually remove (ans and xns
    // are scratch; xnr and temp2 are outputs that must be produced in full).
    // How well the level's own near-field blocks compress: a bounded sample,
    // at the tolerance fraction in use and at a four times tighter one.
    // compression throughput, in situ (decompression is already known: 262 GB/s
    // at the transition's materialization, 365 GB/s in the solve)
    double nb_comp_seconds = 0.0, nb_comp_bytes = 0.0;
    double nb_tight_seconds = 0.0, nb_tight_bytes = 0.0;
    double nb_raw = 0.0, nb_loose = 0.0, nb_tight = 0.0;   // the Schur blocks
    double ne_raw = 0.0, ne_loose = 0.0, ne_tight = 0.0;   // the near-field edges
    int64_t nb_blocks = 0, ne_blocks = 0;
    double tw_max = 0.0, tw_ans = 0.0, tw_xns = 0.0, tw_xnr = 0.0, tw_temp2 = 0.0,
           tw_xnr_orig = 0.0, tw_result = 0.0, tw_t = 0.0, tw_small = 0.0;
    double wr_max_window = 0.0, wr_sum = 0.0, wr_distinct = 0.0;
    int64_t wr_windows = 0, wr_blocks = 0, wr_updates = 0;
    double wr_final_q[4] = {0.0, 0.0, 0.0, 0.0};   // bytes by the quarter of the level their last write fell in
    double pb_raw_bytes = 0.0;
    double pb_stored_bytes = 0.0;
    int64_t pb_blocks = 0;
    int64_t pb_packed_short = 0;  // a block whose packed copy did not fit
    int64_t pb_host_short = 0;    // a group of parent blocks left raw for want of host room
    double pb_check_err = 0.0;
    double pb_check_tol = 0.0;
    double lsch_union_max_wave = 0.0;
    double lsch_union_sum = 0.0;
    double ledge_union_max_wave = 0.0;
    double ledge_union_sum = 0.0;
    int64_t lsch_pairs = 0, lsch_distinct = 0;
    int64_t ledge_pairs = 0, ledge_distinct = 0;
    double fs_tr_bytes = 0.0;
    int64_t fs_tr_sources = 0;
    double fs_check_err = 0.0;    // BPACK_CHECK=fillsrc: worst |decompressed - original|
    double fs_check_tol = 0.0;    //   ... and the tolerance it was held to
    double sk_union_bytes_max_wave = 0.0;
    double sk_union_bytes_sum = 0.0;
    double sk_run_bytes_max_wave = 0.0;
    double sk_run_bytes_sum = 0.0;
    int64_t sk_runs = 0;
    double sk_slice_bytes_max_wave = 0.0;
    double sk_slice_bytes_sum = 0.0;
    double sk_src_bytes_max_wave = 0.0;
    double sk_src_bytes_sum = 0.0;
    int64_t sk_src_pairs = 0;
    int64_t sk_src_distinct = 0;
    int64_t sk_src_waves = 0;
    double tr_block_bytes_total = 0.0;
    double tr_block_bytes_max_chunk = 0.0;
    double tr_block_bytes_sum = 0.0;
    int64_t tr_block_revisits = 0;
    int64_t tr_blocks_max_chunk = 0;
    int64_t transition_fill_gemms = 0;
    double owner_flops = 0.0, transition_flops = 0.0, solve_flops = 0.0;
    // BPACK_TRACE=wave: the floating-point work by device phase, for rates per
    // wave - the sketch (rows and stored blocks, its GEMMs, the ID's QR with
    // pivoting at its full size, an upper bound) and the box region by the
    // marks el[] is timed with (kernel values, gathers and adds count none)
    double fl_sketch = 0.0, fl_sketch_id = 0.0;
    double fl_box[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // replicated CA levels: owner-pass updates read nowhere on this rank, not
    // made - pairs of two eliminated ghosts (since before 2026-10-12), an
    // eliminated ghost's own Schur block and its pairs with boxes beyond the
    // halo (2026-10-12) - their floating-point work and counts
    double dead_pair_flops = 0.0, dead_schur_flops = 0.0, dead_beyond_flops = 0.0;
    int64_t dead_pair_tasks = 0, dead_schur_tasks = 0, dead_beyond_tasks = 0;
    double bytes_up = 0.0;
    double bytes_down = 0.0;
    int64_t boxes = 0;
    int64_t owner_gemms = 0;
    int64_t owner_batches = 0;
    int64_t new_targets = 0;
    int64_t heap_reclaims = 0;  // allocations that waited for the background copies
    double remote = 0.0;          // host time of the remote-generator updates (multi-rank)
    double emit = 0.0, exchange = 0.0;  // generator emission, transport between ranks
    bool device_exchange = false;       // generators sent from device memory
    int64_t exchange_fallbacks = 0;     // MPI buffers the exchange arena could not hold
    int64_t exchange_host = 0;          // ... nor the main heap: in pinned host memory
    int64_t keep_splits = 0;            // keep groups split for want of one range
    int64_t owner_early_installs = 0;   // owner payloads installed at a wave boundary, before their colour's end
    double exchange_host_bytes = 0.0;
    int64_t remote_generators = 0;
    // H2_compress_factors: waves whose X_NR was compressed, waves that fell
    // back to storing it uncompressed, and why.
    // What the kept blocks are made of: X_NR uncompressed and as stored,
    // against everything else (T, X_SR, the LU, pivots); of that, T and X_SR
    // compress with X_NR, the LU and pivots never do and so bound the
    // whole-block ratio.
    double xnr_raw_bytes = 0.0;
    double xnr_stored_bytes = 0.0;
    double other_raw_bytes = 0.0;
    double tx_raw_bytes = 0.0;     // of other_raw_bytes: T and X_SR, and what they came to
    double tx_stored_bytes = 0.0;
    int64_t xnr_waves_compressed = 0;
    int64_t xnr_waves_plain = 0;
    int64_t xnr_not_requested = 0;
    size_t heap_peak = 0;
    // The split of heap_peak: retained (a level or longer, chiefly the
    // solve factors and the fill sources) versus per-wave transient.
    size_t heap_peak_resident = 0;
    size_t heap_peak_transient = 0;
};

inline EliminatorStats& eliminator_stats() {
    static EliminatorStats stats;
    return stats;
}

template<typename T>
struct DeviceMatrixT {
    T* ptr = nullptr;
    int rows = 0;
    int cols = 0;
};

// Device blocks of one level, in the eliminator's layout: the Schur block of
// every box (by Morton index) and the near-field block of every one-hop pair
// (by (lo, hi) key: rows = hi's points, columns = lo's).  Built by the device
// transition and adopted by the next level's eliminator, or copied to the
// host (download_level_blocks) when that level runs on the host.
// H2_compress_factors: a block stored compressed.  Its bytes are the chunks
// end to end, each padded to the stride beside it; `strides` of the level
// holds them from chunk0.  rows and cols are the uncompressed shape, which is
// what the adopting level checks and what a decompression needs.
template<typename T>
struct PackedBlockT {
    char* packed = nullptr;
    size_t bytes = 0;
    uint32_t chunk0 = 0;
    uint32_t nchunks = 0;
    int rows = 0;
    int cols = 0;
    // Where the bytes are.  `packed != nullptr` is one contiguous run, which
    // is what build_parent writes once per level.  `packed == nullptr` means
    // the bytes are in the near-block store, as `nchunks` chunk homes from
    // `home0` - each its own run in B_compact plus an optional tail - because
    // those blocks are re-compressed every time their subset runs and a run
    // that may grow cannot stay contiguous.
    char* owner = nullptr;
    uint32_t home0 = 0;
    bool live = false;      // the bytes are current; false once decompressed
    uint32_t written_seq = 0;  // the near store's scatter that last wrote its homes (0: none)
};

template<typename T>
struct DeviceLevelBlocksT {
    std::unordered_map<int64_t, DeviceMatrixT<T>> schur;
    std::unordered_map<uint64_t, DeviceMatrixT<T>> edges;
    // a nonsymmetric matrix: the near-field blocks are directed, K(row, col)
    // under (col << 32) | row (rows = row's points), both orientations
    bool directed = false;
    // ... the same blocks compressed, when H2_compress_factors asked for it.
    // A key in both is stored raw and compressed (while the raw form is still
    // what the next level reads); a key only here is stored compressed.
    std::unordered_map<int64_t, PackedBlockT<T>> schur_packed;
    std::unordered_map<uint64_t, PackedBlockT<T>> edges_packed;
    std::vector<uint32_t> strides;
    // the allocations those point into, one per group rather than one per
    // block: 7099 separate blocks fragmented the heap badly enough that the
    // transition could not get its next buffer
    std::vector<char*> packed_owned;
};
// The device blocks of a factorization's DataType (double or dcomplex elements).
template<typename DataType>
using DeviceLevelBlocks = DeviceLevelBlocksT<typename DeviceScalar<DataType>::type>;

// Handle of a level's eliminator, independent of the kernel type.
template<typename CoordType, typename DataType>
class LevelEliminatorBase {
public:
    virtual ~LevelEliminatorBase() = default;
    // Eliminates one wave; returns the number of boundary boxes in it.
    virtual int eliminate_wave(const std::vector<int64_t>& wave, int wave_index) = 0;
    // Waits for the factor copies; the level's blocks stay on the device.
    virtual void finish() = 0;
    // Near-field and Schur blocks into the host BoxData (host transition).
    virtual void download_blocks() = 0;
    // The parent level's blocks, assembled on the device from this level's
    // (the device transition); this level's device data is released.
    virtual bool can_build_parent() const = 0;
    // A level that is not eliminated (level 1): its boxes keep all points as
    // skeleton, so the adopted blocks go straight into the transition.
    virtual void adopt_without_elimination() = 0;
    // Multi-rank levels, around each transport of the Color loop: the
    // generators of the last wave's boxes with remote neighbors (their host
    // copies are complete on return), then the device updates from the
    // installed remote generators.
    virtual void emit_generators(PendingFactorUpdates<DataType>& pending) = 0;
    virtual void receive_remote(const std::vector<int64_t>& installed) = 0;
    // Multi-rank levels with CUDA-aware MPI: the whole transport instead,
    // generators sent from and received into device memory.  Returns the
    // MPI time.
    virtual bool device_exchange() const = 0;
    virtual std::chrono::high_resolution_clock::duration exchange() = 0;
    virtual void refresh_remote_skeletons() = 0;
    // A component-owner CA level (owner_plan.hpp): the plan before the
    // waves; after a colour's owned waves, its transport (the colour's
    // payloads out, the colour's arrivals installed on the device).
    virtual void set_owner_plan(const OwnerPlanView* plan) {
        (void)plan;
        throw std::runtime_error("LevelEliminator: owner levels are not supported by this eliminator");
    }
    virtual void owner_exchange(int colour) {
        (void)colour;
        throw std::runtime_error("LevelEliminator: owner levels are not supported by this eliminator");
    }
    // The solve's factors copied to a device block as each wave produces
    // them (device_solve.hpp), while the heap has room: whether this rank
    // still holds all of the level's, then the verdict all ranks agreed on.
    virtual bool solve_factors_kept() const = 0;
    virtual void commit_solve_factors(bool keep) = 0;
    // Before the waves, the level's keep verdict (all ranks alike, from an
    // estimate of the factors' bytes; factorization_driver.hpp
    // decide_solve_keep) and its budget: a level that keeps them skips the
    // host copies of their X_NR, and its waves keep factors while they stay
    // within `budget` bytes (else the level's keep fails).
    // Replicated CA levels: whether this rank expects its near set not to fit
    // raw (asked after set_solve_keep), and the level's decision, taken on
    // every rank alike, to hold the blocks on more than one rank packed
    virtual bool shared_packing_expected() { return false; }
    virtual void set_pack_shared(bool /*pack*/, bool /*fallback*/) {}
    // ... and the end of its ghosts: from the interior group on no rank
    // eliminates a ghost, so no copy of a block on another rank is read
    // again and every block packs by this rank's pressure alone.  Called
    // at the same wave on every rank.
    virtual void replicas_done() {}
    // H2_CA_solve home: the level is solved from each box's home factors, so a
    // ghost's factors are not kept here (as an owner level's owned ghost's)
    virtual void set_home_solve(bool /*home*/) {}
    virtual double shared_packing_estimate() const { return 0.0; }
    virtual double shared_packing_room() const { return 0.0; }
    virtual void set_solve_keep(bool keep, double budget, double expected_kept = 0.0) {
        (void)keep;
        (void)budget;
        (void)expected_kept;
    }
    // The bytes this level's solve factors take (kept or not: the next
    // level's keep estimate), and why its keep failed during the waves (0:
    // it did not; 1: the factors outgrew the budget; 2: no free range of the
    // heap held a wave's factors).
    virtual double level_factor_bytes() const { return 0.0; }
    // What the kept blocks actually hold (compressed, with
    // H2_compress_factors), against level_factor_bytes()'s uncompressed size.
    virtual double level_kept_bytes() const { return 0.0; }
    virtual int solve_keep_failure() const { return 0; }
    // The host copies of the X_NR whose download was skipped, before a host
    // transition reads them.
    virtual void restore_host_copies() {}
    // Whether build_parent may store the parent blocks compressed
    // (H2_compress_factors).  Only the next level's eliminator reads them
    // through their compressed form; the process reduction and the host
    // download both read the raw blocks, so the transition decides.
    virtual void set_pack_parent_blocks(bool pack) { (void)pack; }
    virtual std::unique_ptr<DeviceLevelBlocks<DataType>> build_parent(
        std::vector<BoxData<CoordType, DataType>>& parents) = 0;
};

// The device box path covers this level (the caller has already selected a
// streamed, lazy, owner-deferred Color level).
template<typename CoordType, typename DataType, typename KernelType>
// `lazy_schur` >= 0: the level's H2_lazy_schur mode, for a level ahead (the
// runtime is still configured for the current one); -1: the runtime's.
bool level_eliminator_supported(const TreeLevel<CoordType, DataType>& level, const KernelType* kernel,
                                int dimension, FactorizationMethod method, std::string* reason, bool ca = false,
                                int lazy_schur = -1) {
    auto fail = [&](const char* why) {
        if (reason) *reason = why;
        return false;
    };
    static_assert(std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>,
                  "the GPU box path supports double and complex<double>");
    (void)dimension;
    const char* why = nullptr;
    const bool symmetric = kernel_symmetric(kernel, 0);
    const Evaluator* evaluator = registered_evaluator<DataType>(kernel->gpu_evaluator, &why, symmetric);
    if (evaluator == nullptr) return fail(why);
    if (method != FactorizationMethod::LU) return fail("X_RR is factored by LU on the GPU (use H2_XRR_factor=1)");
    if (!symmetric) {
        // (the device sketch evaluates the K(F, B) rows as K(row, col): inside
        // the sketch with an entry evaluator, else as blocks of their own)
        if (evaluator->collective()) return fail("entries through a collective callback need a symmetric matrix");
        if (ca) return fail("CA levels need a symmetric matrix");
        // (several ranks: the generators carry temp4 and X_RN, sent device to device only)
        if (level.num_active_processes != 1 && !device_exchange_enabled()) {
            return fail("a nonsymmetric matrix's levels on several ranks need the device exchange (CUDA-aware MPI)");
        }
    }
    const bool lazy_far = lazy_schur < 0 ? lazy_far_field_mode() == LazyFarFieldMode::LAZY : lazy_schur > 0;
    const bool generator_near = lazy_schur < 0 ? generator_near_enabled() : lazy_schur > 1;
    if (ca) {
        // replicated CA: lazy far fill, explicit near updates (lazy mode 1)
        if (!lazy_far) return fail("CA levels need H2_lazy_schur >= 1");
        return true;
    }
    if (level.num_active_processes != 1) {
        // other ranks exchange lazy generators only (generated near blocks)
        if (!level.ghost_boxes.empty()) return fail("ghost boxes outside a CA level");
        if (!(generator_near && lazy_far)) {
            return fail("multi-rank levels need H2_lazy_schur=2");
        }
    }
    return true;
}

template<typename CoordType, typename DataType, typename KernelType>
class LevelEliminator : public LevelEliminatorBase<CoordType, DataType> {
    // Device element type (double, or dcomplex for std::complex<double>)
    // and the launch items and blocks over it.
    using S = typename DeviceScalar<DataType>::type;
    using DeviceMatrix = DeviceMatrixT<S>;
    using PackedBlock = PackedBlockT<S>;
    using Blocks = DeviceLevelBlocksT<S>;
    using EvalItem = EvalItemT<S>;
    using GatherItem = GatherItemT<S>;
    using AddStoreItem = AddStoreItemT<S>;
    using IdentityItem = IdentityItemT<S>;
    using SymAddItem = SymAddItemT<S>;
    using ColumnSwapItem = ColumnSwapItemT<S>;
    using OrderedSketchItem = OrderedSketchItemT<S>;
    using TransposeItem = TransposeItemT<S>;
    using QrcpItem = QrcpItemT<S>;
    using SumAddItem = SumAddItemT<S>;
    // flops per multiply-add of S, in units of real ones (statistics)
    static constexpr double kFlopScale = is_complex_scalar<S> ? 4.0 : 1.0;

public:
    using Box = BoxData<CoordType, DataType>;
    using Level = TreeLevel<CoordType, DataType>;
    using Tree = ParallelTree<CoordType, DataType>;

    // occupancy: a tree of occupied boxes (color_unstructured), whose local
    // slabs also hold empty boxes; they are never in a wave nor one hop from
    // a box, and stay eliminated with an empty skeleton.  Its transports
    // follow the unstructured host transport's protocol (exchange()).
    // ca: a replicated CA level (see the notes at the top); owner: a
    // component-owner CA level (owner_plan.hpp), whose plan set_owner_plan
    // supplies before the waves.
    LevelEliminator(Tree* tree, int level_index, KernelType* kernel, double tolerance,
                    std::unique_ptr<Blocks> adopt = nullptr, bool occupancy = false, bool device_sketch = true,
                    bool ca = false, bool owner = false, int compress_factors = 0)
        : tree_(tree), level_(tree->levels[static_cast<size_t>(level_index)]), kernel_(kernel),
          tolerance_(tolerance), heap_(DeviceHeap::instance()), adopt_(std::move(adopt)), occupancy_(occupancy),
          ca_(ca || owner), owner_mode_(owner), compress_factors_(compress_factors)
#ifdef H2_HAVE_NVCOMP
          , xnr_compressor_(DeviceHeap::instance())
          , fill_compressor_(DeviceHeap::instance())
#endif
    {
        // A nonsymmetric matrix: directed near blocks (dkey), X_RS and X_RN
        // beside X_SR and X_NR, temp4 a second fill-source store.  Without
        // the compressed stores and the kept factors yet.
        nonsym_ = !kernel_symmetric(kernel, 0);
        if (nonsym_) compress_factors_ = 0;
        evaluator_ = evaluator_of(kernel->gpu_evaluator);
        gpu_sketch_ = device_sketch && device_sketch_supported(tree, *evaluator_);
        if (nonsym_ && !gpu_sketch_) {
            throw std::runtime_error("LevelEliminator: a nonsymmetric matrix needs the device sketch");
        }
        if (ca_ && !gpu_sketch_) throw std::runtime_error("LevelEliminator: CA levels need the device sketch");
        early_free_sources_ = level_.num_active_processes > 1;
        level_index_ = level_index;
        // (the host sketch would read remote generators' X_NR during the level)
        device_exchange_ = (!ca_ || owner_mode_) && level_.num_active_processes > 1 && device_exchange_enabled() &&
                           gpu_sketch_;
        if (owner_mode_ && !device_exchange_) {
            throw std::runtime_error("LevelEliminator: owner levels need the device exchange (CUDA-aware MPI)");
        }
        if (nonsym_ && level_.num_active_processes > 1 && !device_exchange_) {
            throw std::runtime_error("LevelEliminator: a nonsymmetric matrix's multi-rank levels need the device exchange");
        }
#ifdef H2_HAVE_NVCOMP
        // H2 host memory (2026-10-10): with the fill sources stored compressed
        // and the device sketch, the transition puts a freed source back from
        // its packed temp2^T (collect_fill_source) and rebuilds X_RR_full from
        // the kept LU, so the host copies that served that restore are not
        // made: the remote generators' temp2 (4.94 GB a rank at the 192^3
        // leaf, x4-8 at 384^3) and the kept boxes' X_RR_full (0.99 GB).  No
        // run had read them (0 reclaims in every pressure run).  Where
        // neither device form is there, the restore reads the X_RR_full a
        // remote generator still keeps on the host, or rebuilds it from the
        // LU a failed keep rescued (host_xrr_full_from_lu).
        packed_restore_ = gpu_sketch_ && compress_factors_ >= 1;
        replicated_ = ca_ && !owner_mode_ && compress_factors_ >= 1;
#endif
    }
    void set_owner_plan(const OwnerPlanView* plan) override {
        if (!owner_mode_) throw std::runtime_error("LevelEliminator: a plan for a level that is not an owner level");
        if (plan == nullptr) throw std::runtime_error("LevelEliminator: an owner level without a plan");
        owner_ = plan;
        begin();  // (the level's start reads the plan: which ghosts are eliminated here)
    }
    void owner_exchange(int colour) override;

    bool device_sketch() const { return gpu_sketch_; }

    LevelEliminator(const LevelEliminator&) = delete;
    LevelEliminator& operator=(const LevelEliminator&) = delete;

    ~LevelEliminator() override {
        if (heap_initialized_) heap_.set_reclaimer(nullptr);
        copier_.reset();  // joins the worker
        exchange_copier_.reset();
        if (fill_copy_stream_ != nullptr) {
            evict_fill_copies(true);  // (release_level_data did, unless a throw got here first)
            cudaStreamDestroy(fill_copy_stream_);
            fill_copy_stream_ = nullptr;
        }
        if (heap_initialized_) {
            // a level that never took the device transition hands its kept
            // X_NR over here instead; a throw must not leave the destructor
            try {
                transpose_kept_xnr();
            } catch (const std::exception& e) {
                std::printf("[gpu] level %d: the kept X_NR could not be transposed: %s\n", level_index_, e.what());
                std::fflush(stdout);
            }
        }
        drop_kept(false, "teardown");
        release_level_data();
        if (near_copy_stream_ != nullptr) {   // (release_level_data waited for the last scatter)
            fill_compressor_.hold_staging(nullptr);
            if (near_scatter_done_ != nullptr) cudaEventDestroy(near_scatter_done_);
            if (near_prefetch_done_ != nullptr) cudaEventDestroy(near_prefetch_done_);
            if (near_reuse_ != nullptr) cudaEventDestroy(near_reuse_);
            cudaStreamDestroy(near_copy_stream_);
            near_copy_stream_ = nullptr;
        }
    }

    void begin();
    int eliminate_wave(const std::vector<int64_t>& wave, int wave_index) override;
    void finish() override;
    void download_blocks() override { download_blocks_to_host(); }
    // needs the device fill sources (freed ones are restored from their
    // host copies), unless nothing was eliminated
    bool can_build_parent() const override { return gpu_sketch_ || !eliminated_any_; }
    void adopt_without_elimination() override {
        for (BoxState& st : states_) {
            st.eliminated = true;
            st.skeleton.resize(static_cast<size_t>(st.n));
            for (int i = 0; i < st.n; ++i) st.skeleton[static_cast<size_t>(i)] = i;
            st.skeleton_identity = true;
        }
    }
    void set_pack_parent_blocks(bool pack) override { pack_parents_ = pack; }
    std::unique_ptr<Blocks> build_parent(std::vector<Box>& parents) override;
    void emit_generators(PendingFactorUpdates<DataType>& pending) override;
    void receive_remote(const std::vector<int64_t>& installed) override;
    bool device_exchange() const override { return device_exchange_; }
    // CA levels, after the host's post-elimination assisting gather: the
    // skeletons of the assisting boxes, for the transition.
    void refresh_remote_skeletons() override { sync_remote_boxes(); }
    std::chrono::high_resolution_clock::duration exchange() override;
    bool solve_factors_kept() const override { return keep_solve_ && !keep_failed_; }
    void commit_solve_factors(bool keep) override {
        if (keep && keep_solve_ && !keep_failed_) {
            DeviceSolveStore& store = device_solve_store();
            store.kept[level_index_] = std::move(kept_);
            kept_ = KeptSolveLevel{};
            // (an eviction of this operator's device data copies them to its host boxes first)
            store.spill = [tree = tree_](DeviceSolveStore& s) { materialize_host_factors(tree, s); };
            // the pledged boxes' X_NR is handed to this copy at the transition
            committed_kept_ = &device_solve_store().kept[level_index_];
            skipped_factors_.clear();  // (read from the store from now on)
        } else {
            rescue_skipped_xnr(true);  // (another rank did not keep its factors)
            drop_kept(true, "commit: another rank did not keep");
        }
    }
    void set_solve_keep(bool keep, double budget, double expected_kept = 0.0) override {
        if (!keep) {
            drop_kept(true, "keep declined at the level's start");
            keep_solve_ = false;
        }
        keep_budgeted_ = keep_solve_;
        keep_budget_ = budget;
        kept_expected_ = keep ? expected_kept : 0.0;   // (the hold limit reserves it)
        // (the host copies stay where the host reads X_NR during the level:
        // generators sent through the host, the replica check)
        // (a nonsymmetric matrix: kept as well, for now, with X_RS and X_RN)
        skip_xnr_ = keep_solve_ && !ca_replica_check_enabled() && !nonsym_ &&
                    (device_exchange_ || ca_ || level_.num_active_processes == 1);
    }
    void restore_host_copies() override { rescue_skipped_xnr(); }
    double level_factor_bytes() const override { return factor_bytes_; }
    double level_kept_bytes() const override { return kept_.bytes; }
    int solve_keep_failure() const override { return keep_failure_; }

private:
    struct BoxState {
        int slot = 0;                 // first point slot of the box
        int n = 0;
        bool eliminated = false;
        bool skeleton_identity = true;
        std::vector<int> skeleton;    // local positions, when eliminated
        DeviceMatrix schur;           // current x current, or none
        int64_t skeleton_stamp = -1;  // wave whose metadata holds the skeleton list
        int64_t skeleton_offset = -1;
        const int* d_skeleton = nullptr;  // device skeleton positions (device ID), unless identity
        // lazy fill source (eliminated, device sketch): temp2 stored row-major
        // (r x ntot column-major, i.e. temp2^T) and X_RR_full (r x r)
        const S* temp2t = nullptr;
        const S* xrr_full = nullptr;
        // X_RR_full is the symmetrized X_RR as it stood before getrf factored
        // it in place, and a kept box's LU is resident for the solve, so the
        // stored copy is redundant: launch_lu_expand rebuilds P L U where a
        // reader wants it.  Null for a box whose factors are not kept, and for
        // a remote box, whose LU is on its owner's rank.
        const S* xrr_lu = nullptr;
        const int* xrr_piv = nullptr;
        int ntot = 0;
        int r = 0;
        // boxes of other ranks (assisting boxes): points appended to the
        // table; eliminated once their skeleton has arrived
        bool remote = false;
        // CA levels: a ghost box (another rank's, eliminated here too)
        bool ghost = false;
        bool on_boundary = false;
        size_t persist_bytes = 0;  // of the fill-source block, for the stats
        // multi-rank levels: the fill-source block (temp2^T, X_RR_full) is
        // freed once every local neighbor, its only reader, is eliminated
        char* persist = nullptr;      // X_RR_full
        char* persist_t2t = nullptr;  // temp2^T, in a block of its own so it can go once compressed
        size_t persist_t2t_bytes = 0;
        // a nonsymmetric matrix: temp4 = -X_RR^{-1} X_RN (r x ntot), the
        // column side of the fill, laid out as temp2^T; on a multi-rank level
        // in a block of its own (persist_t4), released with temp2^T
        const S* temp4 = nullptr;
        char* persist_t4 = nullptr;
        size_t persist_t4_bytes = 0;
        // H2_compress_factors: temp2^T compressed a slot at a time, the slots
        // laid end to end in t2t_packed.  A reader names a slot by its first
        // row (that is what lazy_far_locate_endpoint_rows gives it), so
        // t2t_rows maps a first row to its slot; the slot says where its bytes
        // start and which chunks they are, because a slot over kXnrChunkBytes
        // is more than one chunk (about 500 KB a slot at level 4, 2.2 MB at
        // level 2, against a 512 KB chunk).
        struct T2tSlot {
            uint64_t off = 0;      // of the slot's bytes in t2t_packed
            uint32_t chunk0 = 0;   // its first chunk in t2t_strides
            uint32_t nchunks = 0;
            int32_t rows = 0;      // uncompressed rows, for the destination
        };
        char* t2t_packed = nullptr;        // on the device: until this wave's owner pass has read it
        size_t t2t_packed_bytes = 0;
        char* t2t_host = nullptr;          // the same bytes in pinned host memory (fill_host_arena),
                                           // which the later readers decompress from over PCIe
        std::vector<T2tSlot> t2t_slots;
        std::vector<uint32_t> t2t_strides;
        // X_NR^T, the owner pass's other operand, cut into the same slots.
        // Only the owner pass of the wave that made it reads it, so it is
        // built, used and freed inside that wave; t2t_rows indexes both,
        // because the row partition is the same.
        char* xnrt_packed = nullptr;
        size_t xnrt_packed_bytes = 0;
        std::vector<T2tSlot> xnrt_slots;
        std::vector<uint32_t> xnrt_strides;
        std::unordered_map<int, int> t2t_rows;  // first row -> slot
        // H2_compress_factors: this box is kept and its compressed temp2^T is
        // the only copy of the solve factor, so the bytes have to reach the
        // transition, where transpose_kept_xnr() makes the kept X_NR from
        // them.  `reclaimed` means the readers are done and the box behaves as
        // if the fill source were gone - only the bytes are still there.
        bool t2t_pledged = false;
        bool t2t_reclaimed = false;
        int readers = 0;
        // static ID training rows of the box: slots after training_slot0_
        const int* d_training = nullptr;
        int ntraining = 0;
    };

    struct WaveBox {
        Box* box = nullptr;
        BoxState* state = nullptr;
        int n = 0, k = 0, r = 0, ntot = 0;
        std::vector<int64_t> counts;   // rows per one-hop neighbor (current sizes)
        std::vector<int> row0;         // first X_NR row of each neighbor
        // the slot a neighbor position holds in the temp2^T and X_NR^T stores,
        // which skip the neighbors with no rows; -1 where there is none
        std::vector<int32_t> slot_of_pos;
        int64_t skeleton_offset = -1;  // metadata offsets of the local index lists
        int64_t redundant_offset = -1;
        // device buffers
        S* T = nullptr;         // interpolation matrix (k x r), leading dimension ldt
        int ldt = 0;
        S* xrr = nullptr;       // A_RR -> X_RR -> LU(X_RR)
        S* xrr_full = nullptr;  // X_RR before factorization
        S* xsr = nullptr;       // A_SR -> X_SR
        S* s = nullptr;         // A_SS -> new Schur block
        S* tmp1 = nullptr;      // -A_SR^T T
        S* tmp2 = nullptr;      // A_SS T
        S* ans = nullptr;       // A_NS of all neighbors (ntot x k)
        S* xnr = nullptr;       // A_NR -> X_NR (ntot x r)
        S* xns = nullptr;       // temp2 X_SR^T (ntot x k)
        S* temp1 = nullptr;
        S* temp2 = nullptr;
        int* piv = nullptr;
        S* persist_xrr = nullptr;     // X_RR_full, kept for later sketches
        S* persist_temp2t = nullptr;  // temp2^T (r x ntot), kept for later sketches
        // offsets inside the result, source and persist blocks
        size_t off_xrr = 0, off_piv = 0, off_temp1 = 0, off_xrs = 0, off_temp2 = 0, off_xrr_full = 0;
        size_t off_t = 0;            // compact T for the host (device ID)
        int source_group = 0;
        bool has_remote = false;     // a neighbor on another rank: its generator is sent
        size_t off_xnr_orig = 0;     // original X_NR in the source group (generator boxes)
        // a nonsymmetric matrix: A_RS -> X_RS (r x k), A_SN (k x ntot), A_RN
        // -> X_RN (r x ntot, read by the owner pass), the X_SN update X_SR
        // temp4 (k x ntot), the row-swapped copies of X_RS and X_RN, temp3
        // (result block, off_temp3) and temp4 (temp2's group, off_temp4,
        // then the fill source persist_temp4)
        S* xrs = nullptr;
        S* asn = nullptr;
        S* xrn = nullptr;
        S* xsn = nullptr;
        S* xrs_p = nullptr;
        S* xrn_p = nullptr;
        S* temp3 = nullptr;
        S* temp4 = nullptr;
        S* persist_temp4 = nullptr;
        size_t off_temp3 = 0, off_temp4 = 0;
    };

    // (lo, hi) Morton pair of a near-field edge.
    static uint64_t edge_key(int64_t a, int64_t b) {
        const int64_t lo = std::min(a, b), hi = std::max(a, b);
        return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    }
    // A nonsymmetric matrix's near block K(row, col) (rows = row's points):
    // the key of K(hi, lo) is edge_key(lo, hi), that of K(lo, hi) the swapped
    // one.  A symmetric level stores K(hi, lo) only.
    static uint64_t dkey(int64_t row, int64_t col) {
        return (static_cast<uint64_t>(col) << 32) | static_cast<uint64_t>(row);
    }

    bool is_local(int64_t morton) const {
        const int64_t idx = morton - level_.local_morton_start;
        return idx >= 0 && idx < static_cast<int64_t>(states_.size());
    }
    // A box eliminated on this rank: local (an owner level: unless another
    // rank eliminates it), or a ghost of a CA level (an owner level: one of
    // its owned components).
    bool is_writable(int64_t morton) {
        if (is_local(morton)) return owner_ == nullptr || !owner_->foreign(morton);
        if (!ca_) return false;
        auto it = remote_index_.find(morton);
        if (it == remote_index_.end() || !remote_states_[it->second].ghost) return false;
        return owner_ == nullptr || owner_->owned_ghost(morton);
    }
    // A box whose pairs this rank keeps up to date (the owner pass's targets):
    // local, or a ghost eliminated here.  (A replicated CA level: every ghost.)
    bool is_interest(int64_t morton) {
        if (is_local(morton)) return true;
        if (!ca_) return false;
        auto it = remote_index_.find(morton);
        if (it == remote_index_.end() || !remote_states_[it->second].ghost) return false;
        return owner_ == nullptr || owner_->owned_ghost(morton);
    }
    // The host box of a local box or a ghost (a CA level), else null.
    Box* held_box(int64_t morton) {
        if (Box* box = level_.find_local_box(morton)) return box;
        return ca_ ? level_.find_ghost_box(morton) : nullptr;
    }
    Box* writable_box(int64_t morton) { return is_writable(morton) ? held_box(morton) : nullptr; }
    bool has_state(int64_t morton) const { return is_local(morton) || remote_index_.count(morton) != 0; }
    BoxState& state_of(int64_t morton) {
        if (is_local(morton)) return states_[static_cast<size_t>(morton - level_.local_morton_start)];
        auto it = remote_index_.find(morton);
        if (it == remote_index_.end()) {
            throw std::runtime_error("LevelEliminator: box " + std::to_string(morton) +
                                     " is neither local nor an assisting box");
        }
        return remote_states_[it->second];
    }
    // Rows of a ring box of the sketch: its full points or its skeleton.
    struct RingInfo {
        int64_t full = 0;
        const std::vector<int64_t>* skeleton = nullptr;
        bool on_boundary = false;
        bool remote = false;
    };
    RingInfo ring_info(int64_t morton) {
        RingInfo info;
        if (Box* nb = held_box(morton)) {
            info.full = nb->num_points;
            info.skeleton = &nb->skeleton_indices;
            info.on_boundary = nb->on_boundary;
            return info;
        }
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(morton);
        if (it == level_.assisting_box_points_for_kernel_evaluation.end()) {
            throw std::runtime_error("LevelEliminator: ring box " + std::to_string(morton) + " has no point data");
        }
        const auto& assist = level_.assisting_boxes[static_cast<size_t>(it->second)];
        info.full = static_cast<int64_t>(assist.indices.size());
        info.skeleton = &assist.skel_indices;
        info.on_boundary = assist.on_boundary;
        info.remote = true;
        return info;
    }
    // A fill-source block of box `st` whose local neighbors among `hop` are
    // its readers (multi-rank levels), or a level allocation.
    // Whatever a box still holds of its temp2^T: the raw block, the compressed
    // one, or both while the readers have not moved over.
    // Whether a box still has its temp2^T, in either form.
    static bool has_fill_source(const BoxState& st) {
        if (st.t2t_reclaimed) return false;  // held for the kept X_NR, not for reading
        return st.temp2t != nullptr || packed_of(st) != nullptr;
    }
    // The packed temp2^T a reader decompresses from: the device copy while it
    // is there (the wave that made it), else the pinned host copy, which the
    // device reads over PCIe (UVA).
    static const char* packed_of(const BoxState& st) { return st.t2t_packed != nullptr ? st.t2t_packed : st.t2t_host; }
    // The device copies whose host copies are done go back to the heap (at a
    // wave's end, after its owner pass); `all`: wait for them (the level's end).
    struct PendingEvict { cudaEvent_t done = nullptr; std::vector<BoxState*> boxes; };
    std::deque<PendingEvict> pending_evictions_;
    // The copies to the host run on a stream of their own, under the wave's
    // GEMMs: on the main stream they sat in the owner pass's pipeline (+0.5 s
    // a leaf).  The gather they read is synchronized before they start.
    cudaStream_t fill_copy_stream_ = nullptr;
    cudaStream_t fill_copy_stream() {
        if (fill_copy_stream_ == nullptr) {
            check_cuda(cudaStreamCreateWithFlags(&fill_copy_stream_, cudaStreamNonBlocking), "fill copy stream");
        }
        return fill_copy_stream_;
    }
    // The packed temp2^T still without a host copy - the producer makes one
    // only under pressure - copied now that the room is needed, in one batch
    // with one event; evict_fill_copies then frees the device copies as the
    // copies land.
    void stage_fill_copies() {
        std::vector<BoxState*> todo;
        auto collect = [&](BoxState& st) {
            if (st.t2t_packed != nullptr && st.t2t_host == nullptr && st.t2t_packed_bytes > 0) todo.push_back(&st);
        };
        for (BoxState& st : states_) collect(st);
        for (BoxState& st : remote_states_) collect(st);
        if (todo.empty()) return;
        auto& stats = eliminator_stats();
        PinnedArena& arena = fill_host_arena();
        cudaStream_t cs = fill_copy_stream();
        cudaEvent_t ready = make_event();   // the device copies' last writes, on the compute stream
        check_cuda(cudaEventRecord(ready, Context::instance().stream()), "fill store to the host");
        check_cuda(cudaStreamWaitEvent(cs, ready, 0), "fill store to the host");
        check_cuda(cudaEventDestroy(ready), "fill store event");
        for (BoxState* st : todo) {
            st->t2t_host = arena.alloc(st->t2t_packed_bytes);
            check_cuda(cudaMemcpyAsync(st->t2t_host, st->t2t_packed, st->t2t_packed_bytes, cudaMemcpyDeviceToHost, cs),
                       "fill store to the host");
        }
        stats.fs_host_bytes = std::max(stats.fs_host_bytes, static_cast<double>(arena.used()));
        stats.fs_staged_late += static_cast<int64_t>(todo.size());
        PendingEvict pe;
        check_cuda(cudaEventCreateWithFlags(&pe.done, cudaEventDisableTiming), "fill store event");
        check_cuda(cudaEventRecord(pe.done, cs), "fill store event");
        pe.boxes = std::move(todo);
        pending_evictions_.push_back(std::move(pe));
    }
    void evict_fill_copies(bool all) {
        auto& stats = eliminator_stats();
        while (!pending_evictions_.empty()) {
            PendingEvict& pe = pending_evictions_.front();
            if (all) {
                check_cuda(cudaEventSynchronize(pe.done), "fill store event");
            } else {
                const cudaError_t q = cudaEventQuery(pe.done);
                if (q == cudaErrorNotReady) break;  // (later ones were recorded after it)
                check_cuda(q, "fill store event");
            }
            for (BoxState* st : pe.boxes) {
                if (st->t2t_packed == nullptr || st->t2t_host == nullptr) continue;  // (released since)
                stats.persist_bytes -= static_cast<double>(st->t2t_packed_bytes);
                heap_.free(st->t2t_packed);
                st->t2t_packed = nullptr;
            }
            cudaEventDestroy(pe.done);
            pending_evictions_.pop_front();
        }
    }
    // `final` at the level's end, where nothing is kept back.  Mid-level a
    // pledged box gives up its raw temp2^T and its readers but keeps the
    // compressed bytes, which are its kept X_NR until the transition.
    void release_fill_temp2t(BoxState& st, bool final = true) {
        if (st.persist_t4 != nullptr) {  // (a nonsymmetric matrix: the fill's column side goes with it)
            eliminator_stats().persist_bytes -= static_cast<double>(st.persist_t4_bytes);
            heap_.free(st.persist_t4);
            st.persist_t4 = nullptr;
            st.persist_t4_bytes = 0;
        }
        st.temp4 = nullptr;
        if (st.persist_t2t != nullptr) {
            eliminator_stats().persist_bytes -= static_cast<double>(st.persist_t2t_bytes);
            heap_.free(st.persist_t2t);
            st.persist_t2t = nullptr;
            st.persist_t2t_bytes = 0;
        }
        if (!final && st.t2t_pledged && packed_of(st) != nullptr) {
            st.t2t_reclaimed = true;
            st.temp2t = nullptr;
            return;
        }
        if (st.t2t_packed != nullptr) {
            eliminator_stats().persist_bytes -= static_cast<double>(st.t2t_packed_bytes);
            heap_.free(st.t2t_packed);
            st.t2t_packed = nullptr;
        }
        st.t2t_host = nullptr;  // (the arena is reset at the level's end)
        st.t2t_packed_bytes = 0;
        st.t2t_slots.clear();
        st.t2t_strides.clear();
        st.t2t_rows.clear();
        st.temp2t = nullptr;
    }
    void keep_fill_source(BoxState& st, char* block, size_t bytes, char* t2t, size_t t2t_bytes,
                          const std::vector<int64_t>& hop) {
        if (!early_free_sources_) {
            level_allocs_.push_back(block);
            if (t2t != nullptr) level_allocs_.push_back(t2t);
            return;
        }
        st.persist = block;
        st.persist_bytes = bytes;
        st.persist_t2t = t2t;
        st.persist_t2t_bytes = t2t_bytes;
        st.readers = 0;
        for (int64_t m : hop) {
            if (is_writable(m) && !state_of(m).eliminated) ++st.readers;
        }
    }
    void set_remote_skeletons(const std::vector<std::pair<BoxState*, const std::vector<int64_t>*>>& skeletons);
    void sync_remote_boxes(bool skeletons = true);
    // a generator received into device memory (device exchange)
    struct Incoming { int r = 0, k = 0, ntot = 0, nx = 0; const char* bulk = nullptr; };
    void apply_remote_generators(const std::vector<int64_t>& installed,
                                 const std::unordered_map<int64_t, Incoming>* incoming = nullptr);
    // Bulk of a generator sent between devices: temp2 (ntot x r), the rows
    // of the original X_NR that belong to the destination's boxes (nx x r,
    // its one-hop slots in order), X_RS (r x k) and X_RR_full (r x r).  A
    // nonsymmetric matrix sends no X_NR (nx = 0) but, after these, the
    // original X_RN whole (r x ntot: the pair updates K(L, N) read the
    // columns of every N one hop from a destination box), temp4 (r x ntot,
    // the fill's column side) and the solved X_SR, temp1 (k x r).
    struct GenLayout {
        size_t xnr, xrs, xrr, xrn = 0, temp4 = 0, xsr_solved = 0, bytes;
        GenLayout(size_t ntot, size_t r, size_t k, size_t nx, bool nonsym = false) {
            xnr = align_up(ntot * r * sizeof(S));
            xrs = xnr + align_up(nx * r * sizeof(S));
            xrr = xrs + align_up(r * k * sizeof(S));
            bytes = xrr + align_up(r * r * sizeof(S));
            if (nonsym) {
                xrn = bytes;
                temp4 = xrn + align_up(r * ntot * sizeof(S));
                xsr_solved = temp4 + align_up(r * ntot * sizeof(S));
                bytes = xsr_solved + align_up(k * r * sizeof(S));
            }
        }
    };
    // One wave's generator pieces: planned once, because the headers and the
    // outbox offsets are the wave's, and copied a production chunk at a time,
    // because that is where its temp2 is still raw.
    struct GenRun { int row, len; };
    struct GenPiece { int peer; size_t off; const WaveBox* wb; int nx; std::vector<GenRun> runs; };
    std::vector<GenPiece> gen_pieces_;
    void plan_generators(const std::vector<WaveBox>& boxes, const std::vector<int>& active);
    void copy_generators(const std::vector<WaveBox>& boxes, const std::vector<int>& active,
                         const std::vector<char*>& d_sources, const std::vector<char*>& d_temp2,
                         const char* d_xnr, char* d_result, cudaStream_t stream);
    void pack_generators(const std::vector<WaveBox>& boxes, const std::vector<int>& active,
                         const std::vector<char*>& d_sources,
                         const std::vector<char*>& d_temp2, const char* d_xnr,
                         char* d_result, cudaStream_t stream);
    // An owner level (owner_plan.hpp): a box's piece for a recipient is the
    // generator (GenLayout) and, for the box's home rank, the solve's
    // factors after it.
    struct OwnerLayout : GenLayout {
        size_t t, xsr, lu, piv, schur, total;
        OwnerLayout(size_t ntot, size_t r, size_t k, size_t nx, bool full) : GenLayout(ntot, r, k, nx) {
            t = this->bytes;
            xsr = t + (full ? align_up(k * r * sizeof(S)) : 0);
            lu = xsr + (full ? align_up(k * r * sizeof(S)) : 0);
            piv = lu + (full ? align_up(r * r * sizeof(S)) : 0);
            schur = piv + (full ? align_up(r * sizeof(int)) : 0);  // the box's Schur block after its elimination (k x k)
            total = schur + (full ? align_up(k * k * sizeof(S)) : 0);
        }
    };
    // An owner level's wave message: (box, recipient) pieces, planned once a
    // wave from the boxes' shapes, copied a production chunk at a time, sent
    // once the last chunk's copies are on the stream.
    struct OwnerRun { int row, len; };
    struct OwnerPiece { int peer; size_t off; const WaveBox* wb; int nx; bool full; std::vector<OwnerRun> runs; };
    void plan_owner_generators(const std::vector<WaveBox>& boxes, const std::vector<int>& active);
    void copy_owner_generators(const std::vector<WaveBox>& boxes, const std::vector<int>& active,
                               const std::vector<char*>& d_sources, const std::vector<char*>& d_temp2,
                               const char* d_xnr, char* d_result, cudaStream_t stream);
    void finish_owner_generators(cudaStream_t stream);
    struct OwnerEntry {
        int64_t morton = -1;
        int32_t seq = 0;
        OwnerPayload kind = OwnerPayload::COMPACT;
        int r = 0, k = 0, ntot = 0, nx = 0;
        const char* bulk = nullptr;
        std::vector<int64_t> one_hop, counts, skeleton, use_full_set;
    };
    void land_full_factors(const std::vector<const OwnerEntry*>& entries);  // into the kept store, else the host boxes
    // MPI buffers: from the exchange arena, or the main one when it is full
    char* exchange_alloc(size_t bytes) {
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (char* p = arena.try_alloc(bytes)) return p;
        // receive buffers the copier still holds, as their host copies finish
        const auto tw = std::chrono::steady_clock::now();
        auto waited = [&] {
            eliminator_stats().exchange_wait +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - tw).count();
        };
        while (exchange_copier_ && arena.initialized() && arena.capacity() >= bytes) {
            std::vector<char*> done = exchange_copier_->wait_released();
            if (done.empty()) break;
            free_copied_blocks(std::move(done));
            if (char* p = arena.try_alloc(bytes)) {
                waited();
                return p;
            }
        }
        waited();
        ++eliminator_stats().exchange_fallbacks;
        if (char* p = heap_.try_alloc(bytes, ResidentKind::Exchange)) return p;
        // No range of the main heap holds it either: pinned host memory, which
        // the device reads and writes at the same address (2026-10-11: CA3's
        // owner-wave payloads at 8x are 2-7 GB, one block a wave, and a heap
        // with 8-12 GB free in pieces had no range for one).  Slower to read
        // than device memory, and only where nothing else fits.
        PinnedChunkPool::Chunk c = pinned_chunk_pool().take(bytes);
        exchange_host_[c.ptr] = c;
        ++eliminator_stats().exchange_host;
        eliminator_stats().exchange_host_bytes += static_cast<double>(bytes);
        return c.ptr;
    }
    std::unordered_map<const void*, PinnedChunkPool::Chunk> exchange_host_;
    void free_block(const void* ptr) {
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (arena.owns(ptr)) {
            arena.free(ptr);
        } else if (auto h = exchange_host_.find(ptr); h != exchange_host_.end()) {
            // (the pool is not stream-ordered: what still reads it finishes first)
            check_cuda(cudaStreamSynchronize(Context::instance().stream()), "host exchange buffer");
            pinned_chunk_pool().give(h->second);
            exchange_host_.erase(h);
        } else {
            heap_.free(ptr);
        }
    }
    void keep_solve_factors(const std::vector<WaveBox>& boxes, const std::vector<int>& active, char* d_result,
                            const std::vector<char*>& d_temp2, cudaStream_t stream);
    bool keep_solve_factor_group(const std::vector<const WaveBox*>& kept_boxes, char* d_result,
                                 const std::vector<char*>& d_temp2, cudaStream_t stream);
    // H2_CA_solve home: a ghost's LU and pivots, held to the level's end (its
    // owner keeps them for the solve).  Its stored X_RR_full then goes and is
    // rebuilt from them wherever it is read, as its owner's is: the rebuild is
    // not bitwise the stored X_RR_full, and a ghost reading the latter fed
    // its later neighbours a different fill (192^3, 8 ranks: 8374 of 9488
    // ghost copies differed from their owners').
    void keep_ghost_lu(const std::vector<const WaveBox*>& ghosts, char* d_result, cudaStream_t stream) {
        if (ghosts.empty()) return;
        const size_t D = sizeof(S);
        size_t bytes = 0;
        for (const WaveBox* wb : ghosts) {
            const size_t r = static_cast<size_t>(wb->r);
            bytes += align_up(r * r * D) + align_up(r * sizeof(int));
        }
        char* block = heap_.alloc_resident(bytes, ResidentKind::FillXrr);
        level_allocs_.push_back(block);
        std::vector<CopySpan> items;
        int64_t max_words = 0;
        size_t at = 0;
        for (const WaveBox* wb : ghosts) {
            const size_t r = static_cast<size_t>(wb->r);
            char* lu = block + at;
            at += align_up(r * r * D);
            char* piv = block + at;
            at += align_up(r * sizeof(int));
            items.push_back(CopySpan{d_result + wb->off_xrr, lu, static_cast<int64_t>(r * r * D / 4)});
            items.push_back(CopySpan{d_result + wb->off_piv, piv, static_cast<int64_t>(r * sizeof(int) / 4)});
            max_words = std::max(max_words, static_cast<int64_t>(r * r * D / 4));
            wb->state->xrr_lu = reinterpret_cast<const S*>(lu);
            wb->state->xrr_piv = reinterpret_cast<const int*>(piv);
        }
        const size_t list_bytes = items.size() * sizeof(CopySpan);
        if (list_bytes > solve_meta_device_.capacity()) {
            check_cuda(cudaStreamSynchronize(stream), "ghost LU");  // the old list may still be read
            solve_meta_device_.reserve(list_bytes);
        }
        auto* d_items = static_cast<CopySpan*>(solve_meta_device_.reserve(list_bytes));
        check_cuda(cudaMemcpyAsync(d_items, items.data(), list_bytes, cudaMemcpyHostToDevice, stream), "ghost LU");
        launch_copy_spans(d_items, static_cast<int>(items.size()), max_words, stream);
        eliminator_stats().ghost_lu_bytes += static_cast<double>(bytes);
    }
    bool keep_retry_smaller_ = false;
    bool keep_may_split_ = false;   // (only a wave's own keep groups split; landed factors do not)
    // H2_compress_factors: a box whose kept X_NR is still the fill source's
    // compressed temp2^T.  The two are one matrix - temp2, and temp2
    // transposed - so only one is stored, and the transition transposes it
    // into the form the solve reads.
    // A box's factors to keep, wherever they are: a wave's in its result and
    // source blocks, a landed box's in its FULL payload.  `st` is the box's
    // state when the level's X_RR_full rebuild may read the kept LU (null for
    // a landed box, whose X_RR_full arrives with its payload); `pledge` keeps
    // no X_NR here, the transition makes it from the fill source.
    // (a nonsymmetric matrix: xrs and xrn, X_RS (r x k) and X_RN (r x ntot),
    // kept transposed)
    struct KeptSource { int64_t morton = -1; BoxState* st = nullptr; const char* t = nullptr; const char* xsr = nullptr;
                        const char* xnr = nullptr; const char* lu = nullptr; const char* piv = nullptr;
                        const char* xrs = nullptr; const char* xrn = nullptr;
                        int k = 0, r = 0, ntot = 0; bool pledge = false; };
    bool keep_factor_sources(const std::vector<KeptSource>& srcs, cudaStream_t stream);
    bool keep_sources_or_halves(const std::vector<KeptSource>& srcs, cudaStream_t stream);
    struct PledgedXnr { int64_t morton = 0; BoxState* st = nullptr; int ntot = 0; int r = 0;
                        size_t bytes = 0; };  // what its fill source stores now
    bool can_pledge_fill_source(const WaveBox& wb) const {
        if (compress_factors_ < 1 || !gpu_sketch_) return false;
        const BoxState* st = wb.state;
        if (st == nullptr || packed_of(*st) == nullptr || st->t2t_pledged) return false;
        if (st->t2t_slots.empty() || wb.r <= 0 || wb.ntot <= 0) return false;
        int64_t rows = 0;  // the slots have to be all of X_NR's rows, in order
        for (const auto& sl : st->t2t_slots) rows += sl.rows;
        return rows == wb.ntot;
    }
    // Rebuilt rather than stored: see BoxState::xrr_lu.
    static bool xrr_rebuilt(const BoxState& st) {
        return st.xrr_full == nullptr && st.xrr_lu != nullptr && st.xrr_piv != nullptr && st.r > 0;
    }
    // H2 host memory: whether the host can supply a source's X_RR_full to the
    // transition - its copy, or the LU and pivots a failed keep rescued.
    static bool host_xrr_available(const Box& b, int r) {
        if (b.X_RR_full.is_allocated() && b.X_RR_full.rows == r && b.X_RR_full.cols == r) return true;
        return b.X_RR.is_allocated() && b.X_RR.format == MatrixStorage<DataType>::LU_FACTORED && b.X_RR.rows == r &&
               b.X_RR.cols == r && b.X_RR_pivots.size() >= static_cast<size_t>(r);
    }
    // X_RR_full = P L U from getrf's one block (unit lower L, upper U, 1-based
    // pivots undone on the rows last to first, as launch_row_swaps does):
    // r x r, ld r.  r^3/3 flops a box, on the rare paths only.
    static std::vector<DataType> xrr_full_from_lu(const DataType* lu, int64_t ld, const int* piv, int r) {
        using Hs = std::conditional_t<std::is_same_v<DataType, double>, double, std::complex<double>>;
        static_assert(sizeof(Hs) == sizeof(DataType), "host scalar shares the layout");
        std::vector<DataType> out(static_cast<size_t>(r) * static_cast<size_t>(r));
        Hs* o = reinterpret_cast<Hs*>(out.data());
        const Hs* a = reinterpret_cast<const Hs*>(lu);
        for (int j = 0; j < r; ++j) {
            for (int i = 0; i < r; ++i) {
                Hs sum{};
                const int kmax = std::min(i, j);
                for (int k = 0; k <= kmax; ++k) {
                    const Hs l = k < i ? a[i + static_cast<size_t>(k) * ld] : Hs{1.0};
                    sum += l * a[k + static_cast<size_t>(j) * ld];
                }
                o[i + static_cast<size_t>(j) * r] = sum;
            }
        }
        for (int i = r - 1; i >= 0; --i) {
            const int p = piv[i] - 1;
            if (p == i) continue;
            for (int j = 0; j < r; ++j) std::swap(o[i + static_cast<size_t>(j) * r], o[p + static_cast<size_t>(j) * r]);
        }
        return out;
    }
    static std::vector<DataType> host_xrr_full_from_lu(const Box& b) {
        return xrr_full_from_lu(b.X_RR.data.data(), b.X_RR.lda, b.X_RR_pivots.data(), static_cast<int>(b.X_RR.rows));
    }
    // The host does not need this box's X_NR: its factors are kept on the
    // device, and rescue_skipped_xnr reads them back if the host ever asks.
    // The download and the early free of temp2 have to agree about this, so
    // both ask here rather than each spelling it out.
    bool host_skips_xnr(const WaveBox& wb) const {
        return skip_xnr_ && !keep_failed_ && wb.ntot > 0 &&
               kept_.boxes.count(wb.box->morton_index) != 0 &&
               !ghost_factors_elsewhere(wb.box->morton_index);
    }
    // A ghost eliminated here whose factors live with its owner: an owner
    // level's owned ghost (they go home with the FULL payload), or a
    // replicated level's ghost under H2_CA_solve home (the level is solved
    // from its home factors).  Not kept here for the solve.
    bool home_solve_ = false;
    bool ghost_factors_elsewhere(int64_t morton) const {
        if (owner_ != nullptr) return owner_->owned_ghost(morton);
        if (!home_solve_ || !ca_) return false;
        const auto it = remote_index_.find(morton);
        return it != remote_index_.end() && remote_states_[it->second].ghost;
    }
    // ... nor its host copies (the solve factors, its X_NR where the packed
    // store has its temp2^T), except under the replica check, which hashes
    // the ghost copies against their owners'
    bool ghost_host_copies_spent(int64_t morton) const {
        if (owner_ != nullptr) return owner_->owned_ghost(morton);
        return ghost_factors_elsewhere(morton) && !ca_replica_check_enabled();
    }
    // An owner level's owned ghost whose temp2^T is in the packed fill store:
    // the transition puts that back (packed_restore_, as for a remote
    // generator; its payload goes home device to device), so no host copy of
    // its X_NR (temp2) is made.  Its X_RR_full still goes to the host (no
    // kept LU rebuilds it).  Not a solve factor here, so not on the rescue
    // list either.  (CA3's leaf at 768^3 on 64 ranks held ~14 GB a rank more
    // outside the pinned chunks than CA0's, and the host ran out.)
    bool host_drops_ghost_xnr(const WaveBox& wb) const {
        return packed_restore_ && (owner_ == nullptr || device_exchange_) && wb.ntot > 0 && wb.state != nullptr &&
               ghost_host_copies_spent(wb.box->morton_index) && packed_of(*wb.state) != nullptr;
    }
    // H2 host memory: nor its X_RR_full - the transition rebuilds it from the
    // kept LU (an H2_CA_solve home ghost's from the LU it holds to the level's
    // end).  Not on a level whose generators go through the host
    // (emit_generators reads the copy).
    bool host_skips_xrr(const WaveBox& wb) const {
        const int64_t m = wb.box->morton_index;
        const bool lu_here = ghost_factors_elsewhere(m) ? owner_ == nullptr && wb.state != nullptr &&
                                                              wb.state->xrr_lu != nullptr
                                                        : kept_.boxes.count(m) != 0;
        return packed_restore_ && skip_xnr_ && !keep_failed_ && lu_here && !(wb.has_remote && !device_exchange_);
    }
    // The slot a neighbor position holds in the two stores, or -1.
    static int32_t slot_at(const WaveBox& wb, size_t pos) {
        return pos < wb.slot_of_pos.size() ? wb.slot_of_pos[pos] : -1;
    }
    // One owner group's window of source slots: the union it reads of temp2^T
    // and X_NR^T, 0.82 GB at the 192^3 leaf against the 5.14 GB the wave
    // holds raw (measured, see GPU_CA_PLAN.md step 1).
    static constexpr size_t kOwnerWindowBytes = size_t{2} << 30;   // (x2, 2026-10-10: headroom for fewer groups)
    // How many chunks the production half of a wave runs in.  A box's dense
    // output is 2 n_tot r, and the batch used to be the whole wave, so this
    // divides the largest item in the transient.  The consumers still see the
    // whole wave.
    static constexpr size_t kWaveChunkBytes = size_t{2} << 30;   // a chunk's outputs and work at most (the count was a fixed 4 a wave until 2026-10-10)
    void transpose_kept_xnr();
    bool transpose_kept_xnr_group(const std::vector<PledgedXnr>& todo, size_t first, size_t last, size_t raw,
                                  KeptSolveLevel& kl, cudaStream_t stream);
    bool kept_xnr_through_host(const PledgedXnr& p, KeptSolveLevel& kl);
    // The raw bytes one handover group decompresses and transposes.  Twice
    // this is the staging it needs, at the point in the level where the heap
    // has most room.
    static constexpr size_t kKeptTransposeBudget = size_t{512} << 20;   // (x2, 2026-10-10)
    // One fill source to compress: its temp2^T and the slot row counts that
    // partition it (a slot per one-hop neighbor).
    struct FillSrcRef {
        const S* t2t = nullptr;
        int r = 0;
        const std::vector<int64_t>* counts = nullptr;
        BoxState* st = nullptr;   // where the packed form is recorded
    };
    void compress_fill_sources(std::vector<WaveBox>& boxes, const std::vector<int>& active, cudaStream_t stream);
    void compress_fill_source_set(const std::vector<FillSrcRef>& items, cudaStream_t stream,
                                  bool into_xnrt = false);
    // Every box's X_NR^T store, at the end of the wave that read it.
    void release_owner_sources();
    // The chunks of a compressed fill source's temp2^T, appended to a batch
    // that puts it back whole into `dst` (ntot x r elements, as it was laid
    // out).  Returns false if it is not stored compressed.
    bool collect_fill_source(const BoxState& st, S* dst, std::vector<const void*>& h_src,
                             std::vector<void*>& h_dst);
    // One batch for everything collected.
#ifdef H2_HAVE_NVCOMP
    // Whether this block is adopted and still packed (its ptr is null and the
    // key is in the packed store).  A null ptr otherwise means "no block
    // here", which the readers answer from the kernel instead, so the two
    // must never be confused.
    bool is_packed_schur(int64_t morton) const {
        const auto it = packed_schur_.find(morton);
        return it != packed_schur_.end() && it->second.live;
    }
    bool is_packed_edge(uint64_t key) const {
        const auto it = packed_edges_.find(key);
        return it != packed_edges_.end() && it->second.live;
    }
    struct AdoptedBlock { bool schur; uint64_t key; DeviceMatrix* m; };
    enum MatSite { MatBeginSchur = 0, MatOwnerJoin, MatGroup, MatDownload, MatGenerator,
                   MatBuildParent, MatSiteCount };
    // `take`: where a block's raw memory comes from instead of a heap block of
    // its own (null from it: the heap)
    using RawTake = std::function<char*(const AdoptedBlock&, size_t)>;
    void materialize_adopted(const std::vector<AdoptedBlock>& want, cudaStream_t stream,
                             MatSite site = MatSiteCount, const RawTake* take = nullptr);
    // H2 host staging (c), asynchronous: the gather half of a materialization
    // kept apart, so a group's gather can run on the near copy stream under
    // the previous group's GEMMs and be decompressed when its turn comes.
    struct GatherPlan {
        std::vector<AdoptedBlock> want;
        char* gathered = nullptr;          // the chunk-home pieces, contiguous
        std::vector<const void*> src;      // per chunk: its piece, or a contiguous run's chunk (host)
        std::vector<uint32_t> dst_block;   // ... which block of `want`
        std::vector<uint32_t> dst_chunk;   // ... and which of its chunks
        cudaEvent_t done = nullptr;        // recorded after a gather on another stream
    };
    void plan_gather(GatherPlan& g, cudaStream_t stream, cudaEvent_t done);
    void materialize_from_gather(GatherPlan& g, cudaStream_t stream, MatSite site, const RawTake* take = nullptr);
    // An owner group's raw blocks that its pack takes back, in one heap block
    // cut by a bump pointer (2026-10-12): the leaf's owner pass made a heap
    // block for each of them and gave it back after the pack - 2.85M blocks
    // of ~105 KB a leaf, 5.3 s of the heap's host time with the frees (384^3
    // on 8 ranks, 40 GB cards).  Sizes round as the heap's do, so what the
    // heap counts at every decision is unchanged.
    struct GroupArena { char* base = nullptr; size_t bytes = 0, used = 0; };
    GroupArena group_arena_;
    bool in_group_arena(const void* p) const {
        const char* c = static_cast<const char*>(p);
        return group_arena_.base != nullptr && c >= group_arena_.base && c < group_arena_.base + group_arena_.bytes;
    }
    static size_t arena_round(size_t bytes) {
        const size_t a = DeviceHeap::block_alignment();
        return std::max<size_t>(a, (bytes + a - 1) / a * a);
    }
    char* group_arena_take(size_t bytes) {
        const size_t size = arena_round(bytes);
        if (group_arena_.base == nullptr || group_arena_.used + size > group_arena_.bytes) return nullptr;
        char* p = group_arena_.base + group_arena_.used;
        group_arena_.used += size;
        return p;
    }
    // Everywhere outside the wave materialization a null ptr means "no block
    // here", and the reader answers from the kernel or makes a fresh block.
    // It must not also mean "adopted and still packed": that reads as the
    // original A rather than the updated block, with no error.  Any site that
    // acts on a null ptr says so here.
    void assert_not_packed(bool schur, uint64_t key, const char* site) const {
        if (schur ? is_packed_schur(static_cast<int64_t>(key)) : is_packed_edge(key)) {
            throw std::runtime_error(std::string("LevelEliminator: ") + site +
                                     " reached an adopted block that is still packed (key " +
                                     std::to_string(key) + ")");
        }
    }
#else
    void assert_not_packed(bool, uint64_t, const char*) const {}
#endif
    // A near-field block this window writes, keyed by its map slot (stable
    // across the replacement of the block it holds).
    void note_write(const void* slot, size_t bytes) {
        if (slot == nullptr || bytes == 0) return;
        if (wave_writes_.insert(slot).second) wave_write_bytes_ += bytes;
        size_t& seen = level_writes_[slot];
        seen = std::max(seen, bytes);
        auto& lw = level_last_wave_[slot];
        lw = std::max(lw, cur_wave_);
        ++eliminator_stats().wr_updates;
    }
    void close_write_window() {
        auto& stats = eliminator_stats();
        if (!wave_writes_.empty()) {
            stats.wr_max_window = std::max(stats.wr_max_window, static_cast<double>(wave_write_bytes_));
            stats.wr_sum += static_cast<double>(wave_write_bytes_);
            ++stats.wr_windows;
        }
        wave_writes_.clear();
        wave_write_bytes_ = 0;
    }
    void close_write_level() {
        close_write_window();
        auto& stats = eliminator_stats();
        int last = 0;
        for (const auto& [slot, w] : level_last_wave_) last = std::max(last, w);
        for (const auto& [slot, bytes] : level_writes_) {
            stats.wr_distinct += static_cast<double>(bytes);
            // which quarter of the level's waves saw this block's last write
            const auto it = level_last_wave_.find(slot);
            const int w = it == level_last_wave_.end() ? 0 : it->second;
            const int q = last > 0 ? std::min(3, (4 * w) / (last + 1)) : 0;
            stats.wr_final_q[static_cast<size_t>(q)] += static_cast<double>(bytes);
        }
        stats.wr_blocks += static_cast<int64_t>(level_writes_.size());
        level_writes_.clear();
        level_last_wave_.clear();
    }
    std::unordered_set<const void*> wave_writes_;
    size_t wave_write_bytes_ = 0;
    std::unordered_map<const void*, size_t> level_writes_;
    std::unordered_map<const void*, int> level_last_wave_;
    int cur_wave_ = 0;
#ifdef H2_HAVE_NVCOMP
    // Holding the near-field blocks packed turns on how well they compress,
    // which is not the parent blocks' 7.2x by assumption: they are a
    // different matrix at a different point in the elimination.  A bounded
    // sample, because the compressor stages at least the input size.  Both
    // tolerances, because feeding the fill GEMMs lossy inputs argues for a
    // tighter one than the factors get.
    // half a group each, the group size step 3b would use
    static constexpr size_t kNearSampleBudget = size_t{256} << 20;
    void measure_near_blocks(cudaStream_t stream) {
        // measured and recorded; behind a check so it cannot inflate a timing
        if (near_sampled_ || compress_factors_ < 1 || !env::check("fillsrc")) return;
        near_sampled_ = true;
        auto& stats = eliminator_stats();
        const double base = tolerance_ * kXnrToleranceFraction;
        // each kind on its own half of the budget, and reported apart
        for (int kind = 0; kind < 2; ++kind) {
            std::vector<const S*> srcs;
            std::vector<size_t> elems;
            size_t raw = 0;
            auto take = [&](const DeviceMatrix& m) {
                if (m.ptr == nullptr || m.rows <= 0 || m.cols <= 0) return;
                const size_t bytes = static_cast<size_t>(m.rows) * m.cols * sizeof(S);
                if (raw + bytes > kNearSampleBudget) return;
                srcs.push_back(m.ptr);
                elems.push_back(static_cast<size_t>(m.rows) * m.cols);
                raw += bytes;
            };
            if (kind == 0) {
                for (const BoxState& st : states_) take(st.schur);
            } else {
                for (const auto& kv : edges_) take(kv.second);
            }
            if (srcs.empty()) continue;
            double loose = 0.0, tight = 0.0;
            {   // warm-up, discarded: the first call allocates the staging and
                // builds the plan, which showed up as 1.5 GB/s against 114
                std::vector<CompressedXnr> warm;
                fill_compressor_.compress(srcs, elems, base, warm, stream);
                check_cuda(cudaStreamSynchronize(stream), "near block warm-up");
            }
            for (int t = 0; t < 2; ++t) {
                std::vector<CompressedXnr> out;
                // in situ, around the call as it runs: the one rate in this
                // work never measured, and the whole cost of holding the near
                // blocks packed turns on it
                const auto c0 = std::chrono::steady_clock::now();
                const bool ok = fill_compressor_.compress(srcs, elems, t == 0 ? base : base * 0.25, out, stream);
                check_cuda(cudaStreamSynchronize(stream), "near block compression timing");
                const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
                if (!ok) continue;
                stats.nb_comp_seconds += secs;
                stats.nb_comp_bytes += static_cast<double>(raw);
                if (t == 1) {
                    stats.nb_tight_seconds += secs;
                    stats.nb_tight_bytes += static_cast<double>(raw);
                }
                double stored = 0.0;
                for (const auto& c : out) stored += static_cast<double>(c.bytes);
                if (t == 0) loose = stored; else tight = stored;
            }
            if (kind == 0) {
                stats.nb_raw = static_cast<double>(raw);
                stats.nb_loose = loose; stats.nb_tight = tight;
                stats.nb_blocks = static_cast<int64_t>(srcs.size());
            } else {
                stats.ne_raw = static_cast<double>(raw);
                stats.ne_loose = loose; stats.ne_tight = tight;
                stats.ne_blocks = static_cast<int64_t>(srcs.size());
            }
        }
        fill_compressor_.release();
    }
    bool near_sampled_ = false;
    uint64_t owner_group_stamp_ = 1ull << 40;   // skeleton lists of a group's image

    // BPACK_CHECK=heap: every device pointer a launch item carries must lie in
    // the heap's arena.  assert_not_packed() can only catch a site that
    // *branches* on a null pointer; a gather, an eval, a sum-add or a GEMM
    // simply takes the pointer and uses it, so a block nobody materialized
    // reaches the GPU as null + offset and faults in the kernel instead of
    // being reported at its source.  This checks them on the host, naming the
    // list and the index, so one run gives every such site at once.
    int bad_ptrs_ = 0;
    void want_device(const void* ptr, const char* list, const char* field, size_t i) {
        if (ptr != nullptr && heap_.owns(ptr)) return;
        if (bad_ptrs_++ < 16) {
            std::printf("[heap] level %d: %s[%zu].%s = %p is not heap device memory\n",
                        level_index_, list, i, field, ptr);
            std::fflush(stdout);
        }
    }
    void check_items(const std::vector<GatherItem>& v, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i].m <= 0 || v[i].n <= 0) continue;
            want_device(v[i].out, name, "out", i);
            want_device(v[i].src, name, "src", i);
        }
    }
    void check_items(const std::vector<EvalItem>& v, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i].m <= 0 || v[i].n <= 0) continue;
            want_device(v[i].out, name, "out", i);
        }
    }
    void check_items(const std::vector<SumAddItem>& v, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i].rows <= 0 || v[i].cols <= 0) continue;
            want_device(v[i].target, name, "target", i);
        }
    }
    void check_items(const std::vector<AddStoreItem>& v, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i].m <= 0 || v[i].n <= 0) continue;
            want_device(v[i].out, name, "out", i);
            want_device(v[i].a, name, "a", i);
            want_device(v[i].b, name, "b", i);
        }
    }
    // Only the target: a generator's a and b operands are the incoming
    // payload's bulk or an upload of the host box's matrices, both outside
    // the DeviceHeap, so owns() says nothing useful about them.
    void check_targets(const VBatch<S>& b, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < b.entries.size(); ++i) {
            const auto& e = b.entries[i];
            if (e.m <= 0 || e.n <= 0 || e.k <= 0) continue;
            want_device(e.c, name, "c", i);
        }
    }
    void check_batch(const VBatch<S>& b, const char* name) {
        if (!env::check("heap")) return;
        for (size_t i = 0; i < b.entries.size(); ++i) {
            const auto& e = b.entries[i];
            if (e.m <= 0 || e.n <= 0 || e.k <= 0) continue;
            want_device(e.a, name, "a", i);
            want_device(e.b, name, "b", i);
            want_device(e.c, name, "c", i);
        }
    }
    std::unordered_map<char*, size_t> packed_live_;   // group allocation -> live blocks in it

    // ---- the near-block compressed store (GPU_CA_PLAN.md step 3c)
    //
    // B_compact holds each chunk's compressed bytes at a home it keeps for
    // the level, sized at its first compression with no margin.  When a
    // re-compressed chunk is longer, only the excess goes to a B_newtail
    // slot; the home is never abandoned and nothing accumulates.  Slots are
    // recorded by index, so B_newtail can be enlarged by copying slot k to
    // slot k without rewriting any record.
    static constexpr uint32_t kNoSlot = 0xffffffffu;
    // Extents double, from small, so the store costs about what it holds
    // whatever the problem size.  A flat 512 MB minimum made 96^3 worse than
    // no compression at all: the level's compressed near blocks come to
    // 0.13 GB there and the one extent charged 0.54.
    static constexpr size_t kCompactFirst = size_t{32} << 20;
    static constexpr size_t kCompactMax = size_t{512} << 20;
    size_t compact_total_ = 0;
    static constexpr size_t kTailSlot = size_t{64} << 10;   // > the worst growth measured (17 KB)
    struct ChunkHome {
        char* at = nullptr;     // its run in B_compact
        uint32_t here = 0;      // bytes stored there: its first compressed length
        uint32_t total = 0;     // its current compressed length
        uint32_t tail = kNoSlot;
    };
    std::vector<ChunkHome> chunk_homes_;
    std::vector<char*> compact_extents_;
    char* compact_base_ = nullptr;
    size_t compact_at_ = 0, compact_end_ = 0;
    // B_newtail grows by adding an extent, never by reallocating: a copy
    // kernel built earlier in the same group already holds tail_at() pointers,
    // and moving the slots under them is an illegal access (and frees a base
    // the heap then sees freed twice).  Slot k keeps its address for the
    // level, which is also why a slot index is a stable record.
    static constexpr size_t kTailSlotsPerExtent = 256;   // 16 MB at 64 KB a slot
    std::vector<char*> tail_extents_;
    std::vector<uint32_t> tail_free_;

    char* compact_reserve(size_t bytes) {
        if (compact_base_ == nullptr || compact_at_ + bytes > compact_end_) {
            const size_t grown = std::min(kCompactMax, std::max(kCompactFirst, compact_total_));
            const size_t want = std::max(grown, align_up(bytes));
            // H2 host staging (c): the extents are host memory, device-visible;
            // the scatter writes them and the gather reads them over PCIe
            char* e = near_host_arena(level_index_).alloc(want);
            auto& ns = eliminator_stats();
            ns.ns_host_bytes = std::max(ns.ns_host_bytes, static_cast<double>(near_host_arena(level_index_).used()));
            compact_extents_.push_back(e);
            compact_base_ = e;
            compact_at_ = 0;
            compact_end_ = want;
            compact_total_ += want;
        }
        char* p = compact_base_ + compact_at_;
        compact_at_ += bytes;
        return p;
    }
    bool tail_grow() {
        char* e = near_host_arena(level_index_).alloc(kTailSlotsPerExtent * kTailSlot);
        auto& ns = eliminator_stats();
        ns.ns_host_bytes = std::max(ns.ns_host_bytes, static_cast<double>(near_host_arena(level_index_).used()));
        const size_t base = tail_extents_.size() * kTailSlotsPerExtent;
        tail_extents_.push_back(e);
        for (size_t i = 0; i < kTailSlotsPerExtent; ++i) {
            tail_free_.push_back(static_cast<uint32_t>(base + i));
        }
        return true;
    }
    void tail_release(ChunkHome& h) {
        if (h.tail == kNoSlot) return;
        tail_free_.push_back(h.tail);
        h.tail = kNoSlot;
    }
    char* tail_at(uint32_t slot) const {
        return tail_extents_[slot / kTailSlotsPerExtent] + (slot % kTailSlotsPerExtent) * kTailSlot;
    }
    // The near store at a level's end (2026-10-11): a block's chunk homes keep
    // the length of its first compressed write, when it was largest, and the
    // blocks shrink to their skeletons as their boxes are eliminated, so the
    // store's region held 28 GB at the 384^3 CA0 leaf with 2.2 GB live (25 / 1.7
    // for Color) - pinned host memory the transition then sat on while it
    // packed the parents beside it, which killed CA0's interior ranks at
    // 768^3.  The live blocks go into a region of their own as contiguous
    // runs (the form build_parent writes, which every reader takes) and the
    // store's region back to the pinned pool, for the parents and the staged
    // factors to reuse.  A stale entry (its block raw again) restarts fresh.
    void compact_near_store() {
        if (chunk_homes_.empty()) return;
        const auto t0 = std::chrono::steady_clock::now();
        struct Move { PackedBlock* pb; bool store; size_t bytes; };
        std::vector<Move> moves;
        double live = 0.0;
        auto consider = [&](PackedBlock& pb) {
            if (pb.packed == nullptr && pb.nchunks > 0 && !pb.live) {   // stale: its homes go with the region
                pb = PackedBlock{};
                return;
            }
            if (!pb.live) return;
            size_t bytes = 0;
            if (pb.packed == nullptr) {
                for (uint32_t c = 0; c < pb.nchunks; ++c) bytes += chunk_homes_[pb.home0 + c].total;
                moves.push_back(Move{&pb, true, bytes});
            } else if (!heap_.owns(pb.packed)) {   // a run adopted into this level's host region
                for (uint32_t c = 0; c < pb.nchunks; ++c) bytes += packed_strides_[pb.chunk0 + c];
                moves.push_back(Move{&pb, false, bytes});
            }
            live += static_cast<double>(bytes);
        };
        const double held = static_cast<double>(near_host_arena(level_index_).used());
        // (a pass over the entries first: a stale one must not be read below)
        for (auto& kv : packed_schur_) consider(kv.second);
        for (auto& kv : packed_edges_) consider(kv.second);
        if (held < 2.0 * live + 512e6) {   // not worth the copy: leave the store as it is
            return;
        }
        wait_near_scatter_host();
        if (near_copy_stream_ != nullptr) check_cuda(cudaStreamSynchronize(near_copy_stream_), "near store compaction");
        check_cuda(cudaStreamSynchronize(Context::instance().stream()), "near store compaction");
        // Through pageable memory, so the region's own chunks take the live
        // bytes back: copying straight into a new region pinned 2 GB more on a
        // full node, and a rank that pinned slowly held the level end 21 s
        // (768^3 CA0; pinning slows 10x on a fragmented node).
        std::vector<size_t> at(moves.size() + 1, 0);
        for (size_t i = 0; i < moves.size(); ++i) at[i + 1] = at[i] + moves[i].bytes;
        std::vector<char> staged(at.back());
        #pragma omp parallel for schedule(dynamic, 64)
        for (int64_t i = 0; i < static_cast<int64_t>(moves.size()); ++i) {
            const Move& mv = moves[static_cast<size_t>(i)];
            char* o = staged.data() + at[static_cast<size_t>(i)];
            if (!mv.store) {
                std::memcpy(o, mv.pb->packed, mv.bytes);
                continue;
            }
            for (uint32_t c = 0; c < mv.pb->nchunks; ++c) {
                const ChunkHome& h = chunk_homes_[mv.pb->home0 + c];
                const uint32_t in_home = std::min(h.here, h.total);
                std::memcpy(o, h.at, in_home);
                if (h.total > in_home) std::memcpy(o + in_home, tail_at(h.tail), h.total - in_home);
                o += h.total;
            }
        }
        std::vector<uint32_t> chunk0(moves.size(), 0);
        for (size_t i = 0; i < moves.size(); ++i) {
            if (!moves[i].store) continue;
            const PackedBlock& pb = *moves[i].pb;
            chunk0[i] = static_cast<uint32_t>(packed_strides_.size());
            for (uint32_t c = 0; c < pb.nchunks; ++c) packed_strides_.push_back(chunk_homes_[pb.home0 + c].total);
        }
        compact_extents_.clear();
        compact_base_ = nullptr;
        compact_at_ = compact_end_ = 0;
        compact_total_ = 0;
        tail_extents_.clear();
        tail_free_.clear();
        chunk_homes_.clear();
        near_host_arena(level_index_).reset();   // its chunks idle in the pool, for the region below
        PinnedArena& dst = near_compact_arena();
        std::vector<char*> to(moves.size());
        for (size_t i = 0; i < moves.size(); ++i) to[i] = dst.alloc(moves[i].bytes);
        #pragma omp parallel for schedule(dynamic, 64)
        for (int64_t i = 0; i < static_cast<int64_t>(moves.size()); ++i) {
            std::memcpy(to[static_cast<size_t>(i)], staged.data() + at[static_cast<size_t>(i)], moves[static_cast<size_t>(i)].bytes);
        }
        for (size_t i = 0; i < moves.size(); ++i) {
            PackedBlock& pb = *moves[i].pb;
            if (moves[i].store) {
                pb.chunk0 = chunk0[i];
                pb.home0 = 0;
            }
            pb.packed = to[i];
            pb.bytes = moves[i].bytes;
            pb.owner = nullptr;   // (its region resets with the level)
            pb.written_seq = 0;
        }
        // the host owners of the adopted runs went with them
        for (auto it = packed_live_.begin(); it != packed_live_.end();) {
            if (!heap_.owns(it->first)) it = packed_live_.erase(it); else ++it;
        }
        packed_owned_.erase(std::remove_if(packed_owned_.begin(), packed_owned_.end(),
                                           [&](char* p) { return !heap_.owns(p); }),
                            packed_owned_.end());
        auto& ns = eliminator_stats();
        ns.ns_compact_live = live;
        ns.ns_compact_held = held;
        ns.ns_compact_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    void release_near_store() {
        wait_near_scatter_host();   // the last group's scatter, before its region resets
        // (the extents are the level's host region, reset below, as are the
        // runs adopted from the level above)
        compact_extents_.clear();
        compact_base_ = nullptr;
        compact_at_ = compact_end_ = 0;
        compact_total_ = 0;
        tail_extents_.clear();
        tail_free_.clear();
        chunk_homes_.clear();
        near_host_arena(level_index_).reset();
        near_compact_arena().reset();
    }
    // l_first, l_max and l_last of each near-field block's compressed length,
    // over the waves that write it.  A slot allocator would give each block
    // sigma = m * l_first bytes once and overwrite in place, so the margin m
    // it needs is exactly max(l_max / l_first) - and whether that is 1.0 or
    // 2.0 decides whether the scheme bounds memory or not.
    struct NearLen { size_t first = 0, last = 0, max = 0, raw = 0; };
    std::unordered_map<const void*, NearLen> near_len_;
    template<typename TargetT>
    void measure_near_sizes(std::vector<TargetT>& targets, size_t f, size_t l, cudaStream_t stream) {
        std::vector<const S*> srcs;
        std::vector<size_t> elems;
        std::vector<size_t> which;
        for (size_t i = f; i < l; ++i) {
            const TargetT& t = targets[i];
            if (t.block == nullptr || t.block->ptr == nullptr || t.rows <= 0 || t.cols <= 0) continue;
            srcs.push_back(t.block->ptr);
            elems.push_back(static_cast<size_t>(t.rows) * t.cols);
            which.push_back(i);
        }
        if (srcs.empty()) return;
        std::vector<CompressedXnr> comp;
        const double tol = tolerance_ * kXnrToleranceFraction * 0.25;
        if (!fill_compressor_.compress(srcs, elems, tol, comp, stream)) return;
        for (size_t j = 0; j < which.size(); ++j) {
            if (comp[j].bytes == 0) continue;
            NearLen& rec = near_len_[targets[which[j]].block];
            if (rec.first == 0) rec.first = comp[j].bytes;
            rec.last = comp[j].bytes;
            rec.max = std::max(rec.max, comp[j].bytes);
            rec.raw = std::max(rec.raw, elems[j] * sizeof(S));
        }
        fill_compressor_.release();   // the compressed bytes are not kept: measuring only
    }
    void report_near_sizes() {
        if (near_len_.empty()) return;
        auto& st = eliminator_stats();
        for (const auto& [slot, rec] : near_len_) {
            if (rec.first == 0) continue;
            st.nl_raw += static_cast<double>(rec.raw);
            st.nl_first += static_cast<double>(rec.first);
            st.nl_max += static_cast<double>(rec.max);
            st.nl_last += static_cast<double>(rec.last);
            const double ratio = static_cast<double>(rec.max) / static_cast<double>(rec.first);
            st.nl_worst = std::max(st.nl_worst, ratio);
            const size_t bin = ratio <= 1.0001 ? 0 : (ratio <= 1.1 ? 1 : (ratio <= 1.25 ? 2 : (ratio <= 1.5 ? 3 : 4)));
            st.nl_bin[bin] += static_cast<double>(rec.first);
            ++st.nl_blocks;
        }
        near_len_.clear();
    }
    // One fewer block in this allocation; the allocation goes back at zero.
    void release_packed(char* owner) {
        if (owner == nullptr) return;
        auto it = packed_live_.find(owner);
        if (it == packed_live_.end()) return;
        if (--it->second > 0) return;
        packed_live_.erase(it);
        if (heap_.owns(owner)) heap_.free(owner);  // (a host run goes with its level's region)
        packed_owned_.erase(std::remove(packed_owned_.begin(), packed_owned_.end(), owner),
                            packed_owned_.end());
        ++eliminator_stats().np_reclaimed;
    }
    // The packed entry of a block that is about to hold its own bytes again.
    void drop_packed(bool schur, uint64_t key) {
        if (schur) {
            auto it = packed_schur_.find(static_cast<int64_t>(key));
            if (it == packed_schur_.end()) return;
            char* owner = it->second.owner;
            packed_schur_.erase(it);
            release_packed(owner);
        } else {
            auto it = packed_edges_.find(key);
            if (it == packed_edges_.end()) return;
            char* owner = it->second.owner;
            packed_edges_.erase(it);
            release_packed(owner);
        }
    }
    // Which of a group's blocks to pack (replicated CA levels; see block_shared)
    enum class PackOnly { All, Shared, Private };
    // Replicated CA levels: the blocks on more than one rank by the level's
    // decision (pack_shared_), this rank's own under its pressure.  Should
    // the rank still be short with its own packed on a level that decided
    // otherwise, it packs the shared ones too and says so: its ghost copies
    // may then differ from their owners' (a skeleton among them would stop
    // the level above), which beats running out of memory here.
    template<typename TargetT>
    void pack_near_by_rule(std::vector<TargetT>& targets, size_t first, size_t last, cudaStream_t stream,
                           bool pressure) {
        if (replicas_done_) {   // (the interior group: no ghost reads a copy any more)
            if (pressure) pack_near_group(targets, first, last, stream);
            return;
        }
        if (pack_shared_) pack_near_group(targets, first, last, stream, 1.0, PackOnly::Shared);
        if (!pressure) return;
        pack_near_group(targets, first, last, stream, 1.0, PackOnly::Private);
        if (!pack_shared_ && shared_fallback_ && under_pressure()) {
            warn_replica_fallback("short of room with its own blocks packed, the shared ones packed too");
            pack_near_group(targets, first, last, stream, 1.0, PackOnly::Shared);
        }
    }
    // A finished group's near blocks, compressed into the packed store, their
    // raw memory returned.  The tolerance is the tighter one: these feed the
    // fill GEMMs, where the factors only feed the solve, and four times
    // tighter costs about a sixth of the ratio (6.6x -> 5.4x at the leaf).
    // Keeping them raw when the group cannot be packed is a valid degraded
    // form - less saving, same answer - so a failure here is not fatal.
    template<typename TargetT>
    void pack_near_group(std::vector<TargetT>& targets, size_t first, size_t last, cudaStream_t stream,
                         double tol_scale = 1.0, PackOnly only = PackOnly::All) {
        std::vector<const S*> srcs;
        std::vector<size_t> elems;
        std::vector<size_t> which;
        for (size_t i = first; i < last; ++i) {
            const TargetT& t = targets[i];
            if (t.block == nullptr || t.block->ptr == nullptr || t.rows <= 0 || t.cols <= 0) continue;
            if (only != PackOnly::All && block_shared(static_cast<uint64_t>(t.key)) != (only == PackOnly::Shared)) continue;
            if (only == PackOnly::Shared) {
                ++eliminator_stats().np_shared_blocks;
                eliminator_stats().np_shared_raw += static_cast<double>(static_cast<size_t>(t.rows) * t.cols * sizeof(S));
            }
            srcs.push_back(t.block->ptr);
            elems.push_back(static_cast<size_t>(t.rows) * t.cols);
            which.push_back(i);
        }
        if (srcs.empty()) return;
        auto& stats = eliminator_stats();
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<CompressedXnr> comp;
        // the tighter tolerance: these feed the fill GEMMs, not just the solve
        // (tol_scale < 1: a pack that would otherwise add a quantization of
        // its own to a block the owner pass packs again, finer so it does not)
        const double tol = tolerance_ * kXnrToleranceFraction * 0.25 * tol_scale;
        if (!fill_compressor_.compress(srcs, elems, tol, comp, stream)) {
            ++stats.np_groups_plain;
            if (only == PackOnly::Shared) warn_replica_fallback("a group of near blocks stayed raw (no staging)");
            return;
        }
        const auto& cks = fill_compressor_.chunks();
        const char* staging = fill_compressor_.staging();
        std::vector<CopySpan> spans;
        int64_t max_words = 0;
        std::vector<bool> ok(which.size(), false);

        for (size_t j = 0; j < which.size(); ++j) {
            if (comp[j].bytes == 0) continue;
            const TargetT& t = targets[which[j]];
            const bool schur = static_cast<int64_t>(t.key) < 0;
            const uint64_t pkey = schur ? static_cast<uint64_t>(static_cast<int64_t>(~t.key)) : t.key;
            PackedBlock* pb = nullptr;
            if (schur) {
                pb = &packed_schur_[static_cast<int64_t>(pkey)];
            } else {
                pb = &packed_edges_[pkey];
            }
            if (pb->packed != nullptr) {
                // it still holds build_parent's contiguous run; release that
                // and start chunk homes, so home0 is never read as if a
                // contiguous block had one
                char* owner = pb->owner;
                *pb = PackedBlock{};
                release_packed(owner);
            }
            if (pb->nchunks == 0) {
                // first sight: a run of chunk homes.  The chunk count is
                // ceil(elems / chunk_elems) and elems only shrinks as the
                // block's box is eliminated, so this run is big enough for
                // the level.
                pb->home0 = static_cast<uint32_t>(chunk_homes_.size());
                chunk_homes_.resize(chunk_homes_.size() + comp[j].nchunks);
            }
            const size_t spans_before = spans.size();
            bool stored = comp[j].nchunks <= pb_home_room(pb, comp[j].nchunks);
            for (size_t c = 0; stored && c < comp[j].nchunks; ++c) {
                ChunkHome& h = chunk_homes_[pb->home0 + c];
                const auto& ck = cks[comp[j].chunk0 + c];
                const uint32_t len = static_cast<uint32_t>(ck.stride);
                const char* src = staging + ck.offset;
                tail_release(h);              // its old excess, before the new
                if (h.at == nullptr) {
                    h.at = compact_reserve(len);
                    if (h.at == nullptr) { stored = false; break; }
                    h.here = len;
                }
                const uint32_t fits = std::min(len, h.here);
                spans.push_back(CopySpan{src, h.at, static_cast<int64_t>(fits / 4)});
                max_words = std::max(max_words, static_cast<int64_t>(fits / 4));
                if (len > h.here) {
                    const uint32_t excess = len - h.here;
                    if (excess > kTailSlot) { stored = false; break; }
                    if (tail_free_.empty() && !tail_grow()) { stored = false; break; }
                    h.tail = tail_free_.back();
                    tail_free_.pop_back();
                    spans.push_back(CopySpan{src + h.here, tail_at(h.tail),
                                             static_cast<int64_t>(excess / 4)});
                    max_words = std::max(max_words, static_cast<int64_t>(excess / 4));
                    ++stats.np_tails;
                }
                h.total = len;
            }
            if (!stored) {
                // no room: this block stays uncompressed, which is what every
                // block is without this feature.  Its stale entry must go, or
                // a reader would decompress bytes that no longer match it.
                spans.resize(spans_before);
                if (schur) packed_schur_.erase(static_cast<int64_t>(pkey)); else packed_edges_.erase(pkey);
                ++stats.np_no_room;
                continue;
            }
            pb->packed = nullptr;
            pb->live = true;
            pb->bytes = comp[j].bytes;
            pb->nchunks = static_cast<uint32_t>(comp[j].nchunks);
            pb->rows = t.rows;
            pb->cols = t.cols;
            pb->written_seq = near_scatter_seq_ + 1;   // the scatter about to be issued
            ok[j] = true;
            stats.np_raw += static_cast<double>(elems[j] * sizeof(S));
            stats.np_stored += static_cast<double>(comp[j].bytes);
            ++stats.np_blocks;
        }

        // The scatter into the host homes runs on the near copy stream under
        // the next group's work (the compression synchronized the stream, so
        // the staging is complete); the raw blocks were read by the
        // compression and go back now, and the staging is held until the
        // scatter is done (hold_staging).
        scatter_near_async(spans, max_words);
        for (size_t j = 0; j < which.size(); ++j) {
            if (!ok[j]) continue;
            TargetT& t = targets[which[j]];
            if (t.block->ptr == nullptr) continue;
            if (!in_group_arena(t.block->ptr)) heap_.free(t.block->ptr);   // (the group's block: at its end)
            t.block->ptr = nullptr;
        }
        ++stats.np_groups;
        stats.np_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
    // How many chunk homes this block's run can hold.  It was sized at the
    // block's first compression and the count only shrinks, so this is a
    // check against a bug rather than an expected path.
    size_t pb_home_room(const PackedBlock* pb, size_t want) const {
        return pb->home0 + want <= chunk_homes_.size() ? want : 0;
    }
#endif
#ifdef H2_HAVE_NVCOMP
    // A block the wave's gather will read, placed at a known offset of the
    // gather window.  The offset is fixed before the items are staged, which
    // is what lets the gather keep one up-front staging: a GatherItem's src
    // is baked into the uploaded image and cannot be patched afterwards, so
    // the decompressed copy has to land at an address chosen in advance.
    struct WinBlock { bool schur; uint64_t key; size_t offset; int rows; int cols; };
    static constexpr size_t kFillWindowBytes = size_t{1} << 30;   // (x2, 2026-10-10)
    // Holding the near blocks packed, and with it the bounded gather window
    // and the bounded update groups.  On, the leaf peak falls 20.72 -> 17.12 GB
    // on Color and 22.51 -> 19.57 GB on CA0 (192^3, 8 ranks), for 14% and 19%
    // of Factor time and no change in accuracy.  Off restores the single-group
    // owner pass that predates all of this; see GPU_CA_PLAN.md step 3c.
    static constexpr bool kHoldNearBlocksPacked = true;
    // The bytes of near blocks one group writes - the owner pass's groups, and
    // now step 5's packing.  Halving it halves both the store's packing
    // staging and the owner window, the transient's two largest items, for
    // twice the groups: more launches, and a block written in more groups is
    // re-packed more often.  512 -> 256 MB measured -0.48 GB for +0.1 s.
    static constexpr size_t kNearGroupBudget = kHoldNearBlocksPacked ? (size_t{512} << 20)   // (x2, 2026-10-10)
                                                                     : std::numeric_limits<size_t>::max();
    // Memory pressure (2026-10-10).  The near store's packing - an owner
    // group's targets, step 5's fresh blocks, a generator chunk's blocks -
    // and the fill store's eviction to the host each cost time and are only
    // needed when the card is short of room.  A block may be raw or packed
    // at any moment, every reader handles both, so each pack site asks
    // whether the heap is above the hold limit now: the capacity less what
    // the transient budgets may still take.  At 192^3 on 8 ranks nothing
    // packs and the level runs as it did with the store off; at 8x the
    // boxes a rank everything packs.  BPACK_GPU_HEAP_FRACTION shrinks the
    // heap and so forces the packed regime for a test.
    static constexpr size_t kTransientReserve = kOwnerWindowBytes + kFillWindowBytes + 2 * kNearGroupBudget +
                                                kFillWaveScratch + kXnrStagingBudget + 2 * kKeptTransposeBudget;
    // ... less the solve factors this level will accumulate on the card
    // before the next level starts (kept_expected_, the keep rule's own
    // estimate): the raw near set would otherwise take that room first and
    // the factors, growing through the level and at its handover, overflow
    // it.  At 192^3 the term is 2.5 GB (23.1 -> 20.6 GB, nothing packs); at
    // 8x the boxes a rank it is ~20 GB and the level packs from its first
    // wave (2026-10-10).
    size_t hold_limit() const {
        const size_t cap = heap_.capacity();
        const size_t reserve = kTransientReserve + cap / 10 + static_cast<size_t>(std::max(0.0, kept_expected_));
        return cap > reserve ? cap - reserve : 0;
    }
    bool under_pressure() const {
        g_hold_limit_bytes = hold_limit();
        return heap_.used() > g_hold_limit_bytes;
    }
    // The owner pass's and the generator pass's group budgets.  Not scaled
    // with room: a wave's groups are formed while the heap is still under
    // the hold limit, and a 4x window taken then left a 13.8 GB heap without
    // room for the next chunk's gather window (CA3 at a 35% heap, 2026-10-10).
    // The transient reserve behind the hold limit assumes these sizes.
    size_t near_group_budget() const { return kNearGroupBudget; }
    size_t owner_window_bytes() const { return kOwnerWindowBytes; }
    void materialize_into_window(const std::vector<WinBlock>& list, char* window, cudaStream_t stream) {
        if (list.empty() || window == nullptr) return;
        auto& stats = eliminator_stats();
        const auto t0 = std::chrono::steady_clock::now();
        const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
        auto packed_of = [&](const WinBlock& w) -> const PackedBlock* {
            if (w.schur) {
                auto it = packed_schur_.find(static_cast<int64_t>(w.key));
                return it == packed_schur_.end() ? nullptr : &it->second;
            }
            auto it = packed_edges_.find(w.key);
            return it == packed_edges_.end() ? nullptr : &it->second;
        };
        // the chunk-home blocks are in pieces, so gather them contiguous first
        size_t gather_bytes = 0;
        uint32_t need = 0;   // the last scatter that wrote any of them
        for (const WinBlock& w : list) {
            const PackedBlock* pb = packed_of(w);
            if (pb == nullptr || pb->packed != nullptr) continue;
            for (uint32_t c = 0; c < pb->nchunks; ++c) gather_bytes += chunk_homes_[pb->home0 + c].total;
            need = std::max(need, pb->written_seq);
        }
        wait_near_scatter(stream, need);
        char* gathered = gather_bytes > 0 ? heap_.alloc(align_up(gather_bytes), ResidentKind::Gather) : nullptr;
        std::vector<CopySpan> gspans;
        std::vector<const void*> h_src;
        std::vector<void*> h_dst;
        int64_t gmax = 0;
        size_t gat = 0;
        for (const WinBlock& w : list) {
            const PackedBlock* pb = packed_of(w);
            if (pb == nullptr) continue;
            S* dst0 = reinterpret_cast<S*>(window + w.offset);
            if (pb->packed != nullptr) {          // build_parent's contiguous run
                const char* at = pb->packed;
                for (uint32_t c = 0; c < pb->nchunks; ++c) {
                    h_src.push_back(at);
                    h_dst.push_back(dst0 + static_cast<size_t>(c) * chunk_elems);
                    at += packed_strides_[pb->chunk0 + c];
                }
                continue;
            }
            if (gathered == nullptr) continue;
            for (uint32_t c = 0; c < pb->nchunks; ++c) {
                const ChunkHome& h = chunk_homes_[pb->home0 + c];
                char* piece = gathered + gat;
                const uint32_t in_home = std::min(h.here, h.total);
                gspans.push_back(CopySpan{h.at, piece, static_cast<int64_t>(in_home / 4)});
                gmax = std::max(gmax, static_cast<int64_t>(in_home / 4));
                if (h.total > in_home) {
                    const uint32_t rest = h.total - in_home;
                    gspans.push_back(CopySpan{tail_at(h.tail), piece + in_home,
                                              static_cast<int64_t>(rest / 4)});
                    gmax = std::max(gmax, static_cast<int64_t>(rest / 4));
                }
                h_src.push_back(piece);
                h_dst.push_back(dst0 + static_cast<size_t>(c) * chunk_elems);
                gat += h.total;
            }
        }
        if (!gspans.empty()) {
            const size_t list_bytes = gspans.size() * sizeof(CopySpan);
            if (list_bytes > fill_meta_device_.capacity()) {
                check_cuda(cudaStreamSynchronize(stream), "gather window");
                fill_meta_device_.reserve(list_bytes);
            }
            auto* d_items = static_cast<CopySpan*>(fill_meta_device_.reserve(list_bytes));
            check_cuda(cudaMemcpyAsync(d_items, gspans.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                       "gather window");
            launch_copy_spans(d_items, static_cast<int>(gspans.size()), gmax, stream);
            check_cuda(cudaStreamSynchronize(stream), "gather window");
        }
        if (!h_src.empty()) {
            run_fill_batch(h_src, h_dst, stream);
            ++stats.fw_batches;
            for (const WinBlock& w : list) {
                stats.fw_bytes += static_cast<double>(static_cast<size_t>(w.rows) * w.cols * sizeof(S));
            }
        }
        heap_.free(gathered);
        stats.fw_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
#endif
    void run_fill_batch(const std::vector<const void*>& h_src, const std::vector<void*>& h_dst,
                        cudaStream_t stream);
    // The X_NR of boxes kept without a host copy, copied back from the kept
    // factors (the level's own, or committed to the store): before a fill
    // source leaves the device (the transition restores it from the host
    // copy), before the kept factors are dropped, and for a host transition.
    // The level then downloads X_NR again.  `all`: their X_SR and LU too
    // (before the kept factors are dropped; otherwise those stay on the
    // device, materialize_host_factors).
    void rescue_skipped_xnr(bool all = false) {
        skip_xnr_ = false;
        if (skipped_xnr_.empty() && (!all || skipped_factors_.empty())) return;
        check_cuda(cudaStreamSynchronize(Context::instance().stream()), "X_NR host copies");
        const KeptSolveLevel* committed = nullptr;
        auto& store_kept = device_solve_store().kept;
        if (auto it = store_kept.find(level_index_); it != store_kept.end()) committed = &it->second;
        struct Kept { const KeptSolveBox* box; const KeptSolveLevel* level; };
        auto kept_of = [&](const Box* box, const char* what) {
            Kept k{nullptr, nullptr};
            if (auto it = kept_.boxes.find(box->morton_index); it != kept_.boxes.end()) {
                k = Kept{&it->second, &kept_};
            }
            if (k.box == nullptr && committed != nullptr) {
                if (auto it = committed->boxes.find(box->morton_index); it != committed->boxes.end()) {
                    k = Kept{&it->second, committed};
                }
            }
            if (k.box == nullptr) {
                throw std::runtime_error(std::string("LevelEliminator: ") + what + " of box " +
                                         std::to_string(box->morton_index) + " neither on the host nor kept");
            }
            return k;
        };
        auto copy = [&](const void* src, size_t count, const char* what) {
            std::vector<DataType> v(count);
            check_cuda(cudaMemcpy(v.data(), src, count * sizeof(DataType), cudaMemcpyDeviceToHost), what);
            return v;
        };
        int n_plain = 0, n_packed = 0, n_had = 0;
        for (Box* box : skipped_xnr_) {
            if (box->X_NR.is_allocated()) { ++n_had; continue; }
#ifdef H2_HAVE_NVCOMP
            // No kept record: a keep that failed after this box's chunk had
            // already given its raw temp2 back.  Its temp2^T store is still
            // there, and that is X_NR transposed.
            if (kept_.boxes.count(box->morton_index) == 0 &&
                (committed == nullptr || committed->boxes.count(box->morton_index) == 0)) {
                const BoxState& st = state_of(box->morton_index);
                if (packed_of(st) != nullptr && !st.t2t_slots.empty() && st.ntot > 0 && st.r > 0) {
                    const size_t nt = static_cast<size_t>(st.ntot), r = static_cast<size_t>(st.r);
                    std::vector<char> packed(st.t2t_packed_bytes);
                    check_cuda(cudaMemcpy(packed.data(), packed_of(st), packed.size(), cudaMemcpyDefault),
                               "X_NR host copies (store)");
                    std::vector<S> t(nt * r);
                    size_t o = 0;
                    for (const auto& sl : st.t2t_slots) {
                        decompress_xnr_host<S>(packed.data() + sl.off, t.data() + o,
                                               st.t2t_strides.data() + sl.chunk0, sl.nchunks);
                        o += static_cast<size_t>(sl.rows) * r;
                    }
                    std::vector<DataType> v(nt * r);
                    S* x = reinterpret_cast<S*>(v.data());
                    for (size_t j = 0; j < r; ++j) {
                        for (size_t i = 0; i < nt; ++i) x[i + j * nt] = t[j + i * r];
                    }
                    box->X_NR.set_owned(st.ntot, st.r, std::move(v), MatrixStorage<DataType>::FULL);
                    ++n_packed;
                    continue;
                }
            }
#endif
            const auto kept = kept_of(box, "X_NR");
            const KeptSolveBox* kb = kept.box;
#ifdef H2_HAVE_NVCOMP
            if (kb->xnr_pending) {
                // Mid-level: this box's kept X_NR is still the fill source's
                // temp2^T, which the transition has yet to transpose.  Read it
                // from there, a slot at a time, and transpose on the way in.
                const BoxState& st = state_of(box->morton_index);
                if (packed_of(st) == nullptr || st.t2t_slots.empty()) {
                    throw std::runtime_error("LevelEliminator: no X_NR to rescue for box " +
                                             std::to_string(box->morton_index));
                }
                const size_t nt = static_cast<size_t>(kb->ntot), r = static_cast<size_t>(kb->r);
                std::vector<char> packed(st.t2t_packed_bytes);
                check_cuda(cudaMemcpy(packed.data(), packed_of(st), packed.size(), cudaMemcpyDefault),
                           "X_NR host copies (pledged)");
                std::vector<S> t(nt * r);
                size_t o = 0;
                for (const auto& sl : st.t2t_slots) {
                    decompress_xnr_host<S>(packed.data() + sl.off, t.data() + o,
                                           st.t2t_strides.data() + sl.chunk0, sl.nchunks);
                    o += static_cast<size_t>(sl.rows) * r;
                }
                std::vector<DataType> v(nt * r);
                S* x = reinterpret_cast<S*>(v.data());
                for (size_t j = 0; j < r; ++j) {
                    for (size_t i = 0; i < nt; ++i) x[i + j * nt] = t[j + i * r];
                }
                box->X_NR.set_owned(kb->ntot, kb->r, std::move(v), MatrixStorage<DataType>::FULL);
                ++n_packed;
                continue;
            }
#endif
            if (kb->xnr_bytes > 0) {
                // H2_compress_factors: kb->xnr is xnr_bytes of compressed
                // chunks, so reading ntot * r elements from it would run off
                // the end and hand the host nonsense.  Bring the compressed
                // bytes over (far fewer) and decompress, as the solve's
                // materialize_host_factors does.
#ifdef H2_HAVE_NVCOMP
                std::vector<char> packed(kb->xnr_bytes);
                check_cuda(cudaMemcpy(packed.data(), kb->xnr, kb->xnr_bytes, cudaMemcpyDeviceToHost),
                           "X_NR host copies (compressed)");
                std::vector<DataType> v(static_cast<size_t>(kb->ntot) * static_cast<size_t>(kb->r));
                decompress_xnr_host<S>(packed.data(), reinterpret_cast<S*>(v.data()),
                                       kept.level->chunk_strides.data() + kb->xnr_chunk0, kb->xnr_nchunks);
                box->X_NR.set_owned(kb->ntot, kb->r, std::move(v), MatrixStorage<DataType>::FULL);
                ++n_packed;
#else
                throw std::runtime_error("H2 GPU: compressed X_NR on the host needs a build with nvCOMP");
#endif
                continue;
            }
            box->X_NR.set_owned(kb->ntot, kb->r, copy(kb->xnr, static_cast<size_t>(kb->ntot) * kb->r, "X_NR host copies"),
                                MatrixStorage<DataType>::FULL);
            ++n_plain;
        }
        if (env::trace("wave")) {
            int rk = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rk);
            std::printf("  [gpu] rescue on rank %d level %d: X_NR %d packed, %d plain, %d already on the host; "
                        "%zu factor boxes\n", rk, level_index_, n_packed, n_plain, n_had,
                        all ? skipped_factors_.size() : size_t{0});
            std::fflush(stdout);
        }
        skipped_xnr_.clear();
        if (!all) return;
        for (Box* box : skipped_factors_) {
            const auto kept = kept_of(box, "X_SR");
            const KeptSolveBox* kb = kept.box;
            const size_t k = static_cast<size_t>(kb->k), r = static_cast<size_t>(kb->r);
            if (!box->X_SR.is_allocated()) {
                if (kb->xsr_bytes > 0) {
                    // H2_compress_factors: compressed like X_NR, so the bytes
                    // come over and are decompressed here
#ifdef H2_HAVE_NVCOMP
                    std::vector<char> packed(kb->xsr_bytes);
                    check_cuda(cudaMemcpy(packed.data(), kb->xsr, kb->xsr_bytes, cudaMemcpyDeviceToHost),
                               "X_SR host copies (compressed)");
                    std::vector<DataType> v(k * r);
                    decompress_xnr_host<S>(packed.data(), reinterpret_cast<S*>(v.data()),
                                           kept.level->chunk_strides.data() + kb->xsr_chunk0, kb->xsr_nchunks);
                    box->X_SR.set_owned(kb->k, kb->r, std::move(v), MatrixStorage<DataType>::FULL);
#else
                    throw std::runtime_error("H2 GPU: compressed X_SR on the host needs a build with nvCOMP");
#endif
                } else {
                    box->X_SR.set_owned(kb->k, kb->r, copy(kb->xsr, k * r, "X_SR host copies"),
                                        MatrixStorage<DataType>::FULL);
                }
            }
            if (!box->X_RR.is_allocated()) {
                box->X_RR.set_owned(kb->r, kb->r, copy(kb->lu, r * r, "LU host copies"),
                                    MatrixStorage<DataType>::LU_FACTORED);
                box->X_RR_pivots.resize(r);
                check_cuda(cudaMemcpy(box->X_RR_pivots.data(), kb->ipiv, r * sizeof(int), cudaMemcpyDeviceToHost),
                           "pivot host copies");
            }
        }
        skipped_factors_.clear();
    }
    double kept_bytes_total() const {  // (what is on the device: a staged level's store is not)
        double bytes = kept_.bytes;
        for (const auto& [lvl, kl] : device_solve_store().kept) {
            if (!kl.on_host()) bytes += kl.bytes;
        }
        return bytes;
    }
    // The boxes that dropped X_RR_full for their kept LU get it back as a
    // block of their own before the LU goes (drop_kept): rebuilt on the host
    // from the LU, r x r a box, on this rare path only.  Without it the
    // transition (and the sketches still to come) would read the freed LU.
    void restore_dropped_xrr(const char* who) {
        std::vector<BoxState*> todo;
        for (BoxState& st : states_) if (xrr_rebuilt(st)) todo.push_back(&st);
        for (BoxState& st : remote_states_) if (xrr_rebuilt(st)) todo.push_back(&st);
        if (todo.empty()) return;
        {   // (a rare path by design: said whenever it runs)
            int rk = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rk);
            std::printf("  [gpu] level %d rank %d: X_RR_full of %zu boxes rebuilt from the kept LU before it goes (%s)\n",
                        level_index_, rk, todo.size(), who);
            std::fflush(stdout);
        }
        check_cuda(cudaStreamSynchronize(Context::instance().stream()), "X_RR_full from the kept LU");
        auto& stats = eliminator_stats();
        const size_t D = sizeof(DataType);
        for (BoxState* st : todo) {
            const size_t r = static_cast<size_t>(st->r), rr = r * r;
            std::vector<DataType> lu(rr);
            std::vector<int> piv(r);
            check_cuda(cudaMemcpy(lu.data(), st->xrr_lu, rr * D, cudaMemcpyDeviceToHost), "X_RR_full from the kept LU");
            check_cuda(cudaMemcpy(piv.data(), st->xrr_piv, r * sizeof(int), cudaMemcpyDeviceToHost), "X_RR_full from the kept LU");
            const std::vector<DataType> full = xrr_full_from_lu(lu.data(), static_cast<int64_t>(r), piv.data(), st->r);
            char* blk = heap_.alloc_resident(rr * D, ResidentKind::FillXrr);
            check_cuda(cudaMemcpy(blk, full.data(), rr * D, cudaMemcpyHostToDevice), "X_RR_full from the kept LU");
            st->xrr_full = reinterpret_cast<const S*>(blk);
            if (early_free_sources_ && st->persist == nullptr) {
                st->persist = blk;
                st->persist_bytes = rr * D;
            } else {
                level_allocs_.push_back(blk);
            }
            stats.persist_bytes += static_cast<double>(rr * D);
            stats.persist_xrr_bytes += static_cast<double>(rr * D);
            ++stats.xr_restored;
        }
    }
    // `restore_xrr`: false at the teardown, where nothing reads X_RR_full again.
    void drop_kept(bool restore_xrr, const char* who) {
        if (restore_xrr) restore_dropped_xrr(who);
        for (char* b : kept_.blocks) heap_.free(b);
        kept_ = KeptSolveLevel{};
        for (BoxState& st : states_) { st.xrr_lu = nullptr; st.xrr_piv = nullptr; }   // (they pointed into the blocks)
        for (BoxState& st : remote_states_) { st.xrr_lu = nullptr; st.xrr_piv = nullptr; }
        // nothing is kept, so a pledged box owns its compressed temp2^T again
        // and release_level_data() frees it with the rest
        for (PledgedXnr& p : pledged_xnr_) {
            if (p.st != nullptr) p.st->t2t_pledged = false;
        }
        pledged_xnr_.clear();
        committed_kept_ = nullptr;
    }
    int owner_of(int64_t morton) {
        auto it = owner_cache_.find(morton);
        if (it != owner_cache_.end()) return it->second;
        const int rank = level_.owner_rank_of(morton);  // (the level's region ranges)
        owner_cache_.emplace(morton, rank);
        return rank;
    }
    static int current_size(const BoxState& s) {
        return s.eliminated ? static_cast<int>(s.skeleton.size()) : s.n;
    }
    // Point slots of a box's current rows or columns.
    IndexList current_slots(BoxState& s) { return current_slots_in(s, meta_, wave_stamp_); }
    // Skeleton lists without a device copy go into the image `meta`, once
    // per `stamp` (which names the image being built).
    IndexList current_slots_in(BoxState& s, MetaBuilder& meta, uint64_t stamp) {
        IndexList list;
        list.base = s.slot;
        if (s.eliminated && !s.skeleton_identity && s.d_skeleton != nullptr) {
            list.ptr = s.d_skeleton;
        } else if (s.eliminated && !s.skeleton_identity) {
            if (s.skeleton_stamp != stamp) {
                s.skeleton_offset = static_cast<int64_t>(meta.append(s.skeleton));
                s.skeleton_stamp = stamp;
            }
            list.offset = s.skeleton_offset;
        }
        return list;
    }
    // Device times of the last wave's elimination, once its work is done.
    void collect_elimination_marks();

    void run_ids(const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count);
    void sketch_wave(const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count);
    void free_copied_blocks(std::vector<char*>&& blocks);
    void download_blocks_to_host();
    void release_level_data();
    void check_device_kernel();
    // Kernel entries of `items` (in the metadata image `meta`, on the device
    // at md + offset): the application's evaluator (evaluator.hpp), whose
    // scratch comes from the heap: no block freed ahead of a pending read may
    // be live.
    void eval_blocks(const std::vector<EvalItem>& items, MetaBuilder& meta, const char* md, size_t offset, int max_m,
                     int max_n, cudaStream_t stream);

    Tree* tree_;
    Level& level_;
    KernelType* kernel_;
    double tolerance_;
    DeviceHeap& heap_;
    const Evaluator* evaluator_ = nullptr;  // the application's (H2Kernel::gpu_evaluator)
    bool heap_initialized_ = false;
    std::unique_ptr<Blocks> adopt_;  // blocks from the device transition, if any
    std::unordered_set<const void*> lsch_seen_, ledge_seen_;  // blocks any wave read (the stats)
    std::unordered_set<int64_t> sk_src_all_;         // fill sources any wave read (the stats)
    std::unordered_set<const void*> tr_block_seen_;  // transition blocks a chunk has written (the stats)
    std::vector<char*> level_allocs_;            // level-lifetime allocations besides the blocks
    bool gpu_sketch_ = false;
    bool occupancy_ = false;           // empty boxes in the local slab (see the constructor)
    bool ca_ = false;                  // a replicated CA level (ghosts eliminated here too)
    bool early_free_sources_ = false;  // multi-rank level: fill sources released once read
    bool nonsym_ = false;              // a nonsymmetric matrix (see the constructor)
#ifdef H2_HAVE_NVCOMP
    // H2_compress_factors: the blocks adopted from the level below, still in
    // their compressed form.  packed_owned_ holds the group allocations they
    // point into and is freed with the level.
    std::unordered_map<int64_t, PackedBlock> packed_schur_;
    std::unordered_map<uint64_t, PackedBlock> packed_edges_;
    std::vector<uint32_t> packed_strides_;
    std::vector<char*> packed_owned_;
#endif
    int compress_factors_ = 0;         // H2_compress_factors
    bool pack_parents_ = false;        // ... and the parents' blocks are adopted from the device
#ifdef H2_HAVE_NVCOMP
    WaveXnrCompressor<S> xnr_compressor_;
    WaveXnrCompressor<S> fill_compressor_;  // the fill sources' temp2^T, per slot
#endif
    // Fill sources whose local readers have all sketched, oldest first.  They
    // stay on the device for the transition (which otherwise restores them
    // from their host copies) until the heap needs their memory.
    std::deque<std::pair<BoxState*, int64_t>> released_sources_;  // (state, Morton)
    bool reclaim_released_source() {
        while (!released_sources_.empty()) {
            auto [es, morton] = released_sources_.front();
            released_sources_.pop_front();
            // (a box whose X_RR_full was dropped for its LU has no persist
            // block, but its temp2^T is still there to retire)
            if ((es->persist == nullptr && !xrr_rebuilt(*es)) || es->readers > 0) continue;
#ifdef H2_HAVE_NVCOMP
            if (packed_restore_) {
                // H2 host memory: no host copy is made.  The transition puts
                // the temp2^T back from its packed form, and X_RR_full from
                // the kept LU or the copy a remote generator keeps on the
                // host; what it could not put back stays on the device.
                auto& stats = eliminator_stats();
                const Box* hb = es->remote ? level_.find_generator_box(morton) : held_box(morton);
                const bool xrr_ok = xrr_rebuilt(*es) || (hb != nullptr && host_xrr_available(*hb, es->r));
                bool freed = false;
                if (es->persist != nullptr && xrr_ok) {
                    stats.persist_bytes -= static_cast<double>(es->persist_bytes);
                    heap_.free(es->persist);  // its readers' launches are all ordered before later ones
                    es->persist = nullptr;
                    es->xrr_full = nullptr;
                    freed = true;
                }
                if (es->persist_t2t != nullptr && packed_of(*es) != nullptr) {
                    stats.persist_bytes -= static_cast<double>(es->persist_t2t_bytes);
                    heap_.free(es->persist_t2t);
                    es->persist_t2t = nullptr;
                    es->persist_t2t_bytes = 0;
                    es->temp2t = nullptr;  // the readers take the packed form now
                    freed = true;
                }
                if (!freed) continue;
                ++stats.rc_sources;
                return true;
            }
#endif
            rescue_skipped_xnr();  // (a restore reads the host copies)
            eliminator_stats().persist_bytes -= static_cast<double>(es->persist_bytes);
            heap_.free(es->persist);  // its readers' launches are all ordered before later ones
            es->persist = nullptr;
            release_fill_temp2t(*es, false);
            es->temp2t = nullptr;
            es->xrr_full = nullptr;
            return true;
        }
        return false;
    }
    bool eliminated_any_ = false;  // a wave ran (the transition may need fill sources)
    char* wave_sketch_ = nullptr;  // sketches of the current wave; T in place after the device ID
    std::vector<char*> wave_sketch_extra_;  // ... blocks of boxes whose ID target grew (adaptive rows)
    std::unique_ptr<HostCopier> copier_;  // background copies of the factors into BoxData
    // device exchange: host copies of the received generators, apart from the
    // wave downloads so the receive buffers return to the exchange arena early
    std::unique_ptr<HostCopier> exchange_copier_;

    std::vector<BoxState> states_;
    std::deque<BoxState> remote_states_;               // stable addresses (plans keep pointers)
    std::unordered_map<int64_t, size_t> remote_index_;  // Morton -> remote_states_
    std::unordered_map<uint64_t, DeviceMatrix> edges_;
    PointTable points_;
    int64_t num_slots_ = 0;                            // point slots: local, training, then remote
    int training_slot0_ = 0;                           // first slot of the static ID training points
    // slots for the far points an adaptive round samples (H2_ID_proxy 2),
    // rewritten every round; a larger region is taken at the table's end
    int scratch_slot0_ = 0, scratch_slots_ = 0;
    std::vector<double> scratch_xyz_;
    int place_scratch_points(const std::vector<int64_t>& global_ids, cudaStream_t stream);
    double* d_xyz_ = nullptr;
    int64_t* d_ids_ = nullptr;
    std::vector<int64_t> last_wave_;                   // boxes of the last wave (their generators)
    int last_wave_index_ = -1;

    // device exchange (multi-rank levels, CUDA-aware MPI)
    bool device_exchange_ = false;
    // H2 host memory: a freed fill source is put back from its packed
    // temp2^T and the kept LU, so no temp2 of a remote generator and no
    // X_RR_full of a kept box is copied to the host (see the constructor)
    bool packed_restore_ = false;
    int level_index_ = 0;
    double fill_frac_ = 0.2;   // the fill store's quantization step as a fraction of tol_comp (decide_fill_precision)
    // Replicated CA levels with compression (2026-10-11): a near block that
    // exists on more than one rank - both its boxes are ghosts here or held
    // by another rank (its ghosts or assisting boxes) - is held packed after
    // every write, whatever the rank's room, and one private to this rank
    // only under this rank's pressure.  Taken by each rank from its own heap
    // for every block, the decision packed a ghost's blocks on one rank and
    // left the owner's raw, and the quantization made the copies differ
    // (192^3 CA0 at a 35% heap: 5219 of 9488 ghost copies, skeletons among
    // them) - the crash of CA0 above its leaf at 8x.  A level-wide decision
    // taken a wave late could not follow the green group's first wave, which
    // grows the near set by ~20 GB at 8x.  A private block cannot reach a
    // replicated box's elimination: the other copy's rank does not have it.
    bool replicated_ = false;
    // Whether this level holds its shared blocks packed: decided at the level's
    // start on every rank alike (the driver's MAX of shared_packing_expected),
    // because packing them after every write costs ~1.6 s of a 4 s leaf at
    // 192^3 on 8 ranks, where nothing needs packing.
    bool pack_shared_ = false;
    bool replicas_done_ = false;       // the interior group: see replicas_done()
    double shared_grow_ = 0.0;         // the estimate behind the decision (for the level's report)
    bool shared_packing_expected() override {
        if (!replicated_) return false;
        // Above the leaf the level's blocks are on the card already (adopted from
        // the transition or uploaded here) and the owner pass makes few new
        // ones, so the room it needs is a wave's working set.  At the leaf they
        // are made during the level: the near set raw at its peak, estimated as
        // the full-size blocks of the boxes and their one-hop pairs - of the
        // boxes the boundary groups eliminate (the rank's boundary boxes and its
        // ghosts), since the rule ends with them (replicas_done) and the
        // interior group packs by pressure.  Not half of it: a rank that runs
        // short with the shared blocks raw packs them alone and its ghost
        // copies differ (768^3 on 64 ranks, 80 GB cards: half the count, 37.0
        // GB against 40.9 GB of room, let 38 of 64 ranks fall back in the green
        // group, whose raw peak was 41.5 GB at the wave ends, and the level
        // above stopped on a near block shape).
        const bool leaf = level_index_ == static_cast<int>(tree_->levels.size()) - 1;
        double grow = 0.1 * static_cast<double>(heap_.capacity());
        if (leaf) {
            double full = 0.0;
            auto add = [&](const Box& b) {
                const double n = static_cast<double>(b.num_points);
                full += n * n * (1.0 + 0.5 * static_cast<double>(b.one_hop.size()));
            };
            for (int64_t li : level_.boundary_id) add(level_.local_boxes[static_cast<size_t>(li)]);
            for (const Box& g : level_.ghost_boxes) add(g);
            grow = full * static_cast<double>(sizeof(S));
        }
        shared_grow_ = grow;
        return static_cast<double>(heap_.used()) + grow > static_cast<double>(hold_limit());
    }
    double shared_packing_estimate() const override { return shared_grow_; }
    double shared_packing_room() const override {
        const double room = static_cast<double>(hold_limit()) - static_cast<double>(heap_.used());
        return room;
    }
    void set_pack_shared(bool pack, bool fallback) override {
        pack_shared_ = replicated_ && pack;
        shared_fallback_ = fallback;
    }
    // (H2_CA_pack_shared never: a rank short of room does not pack them alone)
    bool shared_fallback_ = true;
    void replicas_done() override {
        replicas_done_ = true;
        pack_shared_ = false;
    }
    void set_home_solve(bool home) override { home_solve_ = home && ca_ && owner_ == nullptr; }
    std::vector<char> local_shared_;   // by local index: on another rank too
    bool local_shared_ready_ = false, local_shared_all_ = false;
    bool box_shared(int64_t m) {
        if (!is_local(m)) return true;
        if (!local_shared_ready_) {
            local_shared_ready_ = true;
            local_shared_.assign(states_.size(), 0);
            local_shared_all_ = !level_.ca_halo_sends_valid;   // (unknown: every block shared)
            for (const auto& send : level_.ca_halo_sends) {
                for (const auto* list : {&send.ghosts, &send.assisting}) {
                    for (int64_t b : *list) {
                        if (is_local(b)) local_shared_[static_cast<size_t>(b - level_.local_morton_start)] = 1;
                    }
                }
            }
        }
        return local_shared_all_ || local_shared_[static_cast<size_t>(m - level_.local_morton_start)] != 0;
    }
    bool block_shared(uint64_t key) {
        const int64_t sk = static_cast<int64_t>(key);
        if (sk < 0) return box_shared(~sk);   // a Schur block
        return box_shared(static_cast<int64_t>(key >> 32)) && box_shared(static_cast<int64_t>(key & 0xffffffffu));
    }
    bool replica_warned_ = false;
    static double tr_rss_gb() { double rss = 0.0, hwm = 0.0; host_memory_gb(rss, hwm); return rss; }
    static int tr_rank() { int r = 0; MPI_Comm_rank(MPI_COMM_WORLD, &r); return r; }
    void warn_replica_fallback(const char* what) {
        if (!replicated_ || replica_warned_) return;
        replica_warned_ = true;
        int rk = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        std::printf("[gpu] level %d rank %d: %s on a replicated CA level; its ghost copies may differ from "
                    "their owners'\n", level_index_, rk, what);
        std::fflush(stdout);
    }
    bool fill_frac_decided_ = false;
    // The fill store's step (2026-10-11): kFillToleranceFraction where the
    // level can afford it, else the kept X_NR's.  Decided at the level's first
    // pack, when the keep rule's estimate is in the hold limit and the level's
    // working set is on the card: the store's bound at the fine step (27 n^2
    // an n-point box, over kFillFineRatio) must fit in half the room under the
    // hold limit.  The leaf never: its own skeletons do not depend on its fill
    // store (0.06% measured) and its store is the largest.  What the rule
    // protects is the host: under pressure the store goes to the host, and
    // at 768^3 on 64 ranks the level above the leaf packed fine put 6 GB a rank
    // more there and the node ran out (the level's bound 110 GB, the room 10).
    void decide_fill_precision() {
        if (fill_frac_decided_) return;
        fill_frac_decided_ = true;
        fill_frac_ = kXnrToleranceFraction;
        const bool leaf = level_index_ == static_cast<int>(tree_->levels.size()) - 1;
        double bound = 0.0;
        for (const Box& b : level_.local_boxes) {
            // a box and its neighbors (27 on a full grid; with occupancy the
            // box's own one-hop list, which a sparse tree keeps much shorter)
            const double nb = occupancy_ ? static_cast<double>(b.one_hop.size()) + 1.0 : 27.0;
            bound += nb * static_cast<double>(b.num_points) * static_cast<double>(b.num_points) * sizeof(S);
        }
        const size_t limit = hold_limit();
        const size_t used = heap_.used();
        const double room = limit > used ? static_cast<double>(limit - used) : 0.0;
        const bool fine = !leaf && bound / kFillFineRatio < room / 2.0;
        if (fine) fill_frac_ = kFillToleranceFraction;
        int rk = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        if (rk == 0) {
            std::printf("  [gpu] level %d fill store: packed at %.2g of tol_comp (%s; bound %.1f GB at the fine step, room %.1f GB under the hold limit)\n",
                        level_index_, fill_frac_, leaf ? "the leaf" : fine ? "fits" : "does not fit",
                        bound / kFillFineRatio / 1e9, room / 1e9);
            std::fflush(stdout);
        }
    }
    struct PeerOut { std::vector<int64_t> header; size_t offset = 0, bytes = 0; };
    // An owner level's transport (owner_progress): an owned wave's payloads
    // are one block of the exchange arena with a segment per recipient,
    // sent (a header and the bulk per recipient) once the wave's device
    // copies are done, and released once sent.
    struct OwnerWaveSend {
        struct Message { int peer; std::vector<int64_t> header; const char* ptr; size_t bytes; };
        int colour = 0;
        char* block = nullptr;
        cudaEvent_t ready = nullptr;
        bool issued = false;
        std::vector<Message> messages;
        std::vector<MPI_Request> reqs;
    };
    // The wave message being built (plan_owner_generators): its pieces and
    // their segment offsets, the send with its block and headers, finished
    // once the wave's last production chunk has copied into it.
    struct OwnerWavePlan { std::vector<OwnerPiece> pieces; std::map<int, size_t> seg_off; OwnerWaveSend send;
                           bool open = false; };
    OwnerWavePlan owner_plan_;
    // A wave payload of another owner, received as its header arrives, its
    // bulk straight into the exchange arena; installed at its colour's end.
    struct OwnerWaveRecv {
        int source = -1, colour = 0;
        int32_t seq = 0;
        size_t bytes = 0;
        char* buf = nullptr;
        std::vector<int64_t> header;
        std::vector<MPI_Request> reqs;
    };
    // receive buffers read by launches still in flight: freed once they ran
    struct OwnerPendingFree { cudaEvent_t done = nullptr; std::vector<char*> blocks; };
    void owner_progress();
    const OwnerPlanView* owner_ = nullptr;
    bool owner_mode_ = false;
    std::list<OwnerWaveSend> owner_sends_;
    std::list<OwnerWaveRecv> owner_recvs_;
    void install_owner_messages(int colour, const std::vector<OwnerWaveRecv*>& msgs, double t_wait, size_t sent_to,
                                size_t from_ranks);
    void owner_install_arrived(int colour);
    std::list<OwnerPendingFree> owner_pending_free_;
    std::map<std::pair<int, int>, int> owner_sent_count_;  // (peer, colour) -> wave messages sent
    std::map<std::pair<int, int>, int> owner_recv_count_;  // (source, colour) -> wave messages received
    std::map<std::pair<int, int>, int> owner_marks_;       // (source, colour) -> wave messages it sent (its colour-end marker)
    std::map<int, uint64_t> owner_sent_bytes_;             // colour -> payload bytes sent
    int64_t owner_seq_ = 0;                                 // colour transports so far
    char* outbox_ = nullptr;                           // generator bulk of the last wave, a segment per rank
    std::map<int, PeerOut> peer_out_;                  // rank -> its generators (header, segment)
    std::vector<int64_t> since_transport_;             // local boxes eliminated since the last transport
    std::unordered_map<int64_t, std::vector<int>> requesters_;  // local box -> ranks holding it (assisting)
    bool registered_ = false;                          // requests exchanged (first transport of the level)
    int64_t exchange_seq_ = 0;                         // transports of the level so far (in every header)
    struct RemoteHost { std::vector<DataType> temp2, xrr, temp4; int ntot = 0, r = 0; };  // (temp4: nonsymmetric)
    std::unordered_map<int64_t, std::unique_ptr<RemoteHost>> remote_host_;  // background host copies
    std::unordered_map<int64_t, int> owner_cache_;

    // the solve's factors kept on the device (device_solve.hpp)
    bool keep_solve_ = false;
    bool keep_failed_ = false;
    bool keep_budgeted_ = false;  // decided before the waves (set_solve_keep)
    double keep_budget_ = 0.0;    // then: the bytes its kept factors may take
    double kept_expected_ = 0.0;  // what its kept store is expected to hold on the device by its end
    double factor_bytes_ = 0.0;   // the eliminated boxes' factors, as keep_solve_factors lays them out
    int keep_failure_ = 0;        // solve_keep_failure
    bool skip_xnr_ = false;       // kept boxes' X_NR, X_SR and LU not downloaded
    std::vector<Box*> skipped_xnr_;      // ... the boxes whose X_NR stayed on the device (rescue_skipped_xnr)
    std::vector<Box*> skipped_factors_;  // ... whose X_SR and LU did
    KeptSolveLevel kept_;
    std::vector<PledgedXnr> pledged_xnr_;       // in the order they were kept
    KeptSolveLevel* committed_kept_ = nullptr;  // the store's copy, after the commit
    DeviceBuffer solve_meta_device_;
    DeviceBuffer fill_meta_device_;  // span lists of the fill-source gathers
    // H2 host staging (c), asynchronous.  The near store's scatters (device
    // staging -> host homes) run on a stream of their own under the next
    // group's work, and a group's gather (host homes -> device staging) is
    // prefetched under the previous group's GEMMs.  A reader waits only for
    // the scatter that last wrote one of its blocks (PackedBlock::written_seq
    // against the scatters the host knows complete); the compressor holds its
    // staging until the scatter reading it is done.
    cudaStream_t near_copy_stream_ = nullptr;
    cudaStream_t near_copy_stream() {
        if (near_copy_stream_ == nullptr) {
            check_cuda(cudaStreamCreateWithFlags(&near_copy_stream_, cudaStreamNonBlocking), "near copy stream");
        }
        return near_copy_stream_;
    }
    cudaEvent_t near_scatter_done_ = nullptr;    // recorded after each scatter
    cudaEvent_t near_prefetch_done_ = nullptr;   // ... after a prefetched gather
    cudaEvent_t near_reuse_ = nullptr;           // the compute stream past a decompression: its scratch may be reused
    uint32_t near_scatter_seq_ = 0;              // scatters issued
    uint32_t near_scatter_synced_ = 0;           // ... the host has waited for
    DeviceBuffer near_scatter_meta_;             // the scatter's span list (the gathers' is reused under it)
    DeviceBuffer near_prefetch_meta_;            // the prefetched gather's
    static cudaEvent_t make_event() {
        cudaEvent_t e = nullptr;
        check_cuda(cudaEventCreateWithFlags(&e, cudaEventDisableTiming), "near copy event");
        return e;
    }
    void scatter_near_async(const std::vector<CopySpan>& spans, int64_t max_words) {
        if (spans.empty()) return;
        wait_near_scatter_host();   // its span list and the compressor's staging are reused
        const size_t list_bytes = spans.size() * sizeof(CopySpan);
        auto* d_items = static_cast<CopySpan*>(near_scatter_meta_.reserve(list_bytes));
        cudaStream_t cs = near_copy_stream();
        if (near_scatter_done_ == nullptr) near_scatter_done_ = make_event();
        check_cuda(cudaMemcpyAsync(d_items, spans.data(), list_bytes, cudaMemcpyHostToDevice, cs),
                   "near blocks stored");  // (pageable source: staged before the call returns)
        launch_copy_spans(d_items, static_cast<int>(spans.size()), max_words, cs);
        check_cuda(cudaEventRecord(near_scatter_done_, cs), "near scatter event");
        ++near_scatter_seq_;
        fill_compressor_.hold_staging(near_scatter_done_);
    }
    void wait_near_scatter_host() {
        if (near_scatter_synced_ == near_scatter_seq_) return;
        check_cuda(cudaEventSynchronize(near_scatter_done_), "near scatter");
        near_scatter_synced_ = near_scatter_seq_;
    }
    // `stream` is about to read blocks last written by scatter `need`
    void wait_near_scatter(cudaStream_t stream, uint32_t need) {
        if (need == 0 || need <= near_scatter_synced_ || near_scatter_done_ == nullptr) return;
        check_cuda(cudaStreamWaitEvent(stream, near_scatter_done_, 0), "near scatter");
    }
    DeviceBuffer t2t_meta_device_;   // transpose lists of the fill-source sub-groups

    int64_t wave_stamp_ = 0;
    MetaBuilder& meta_ = pinned_pool().meta;
    DeviceBuffer meta_device_;
    MetaBuilder& owner_meta_ = pinned_pool().owner_meta;  // a wave's owner pass
    DeviceBuffer owner_meta_device_;
    // Of the last wave, read at the next synchronization.  One set per
    // production chunk (phases 0-5) and one for the tail (6-7), each with the
    // phase its first mark belongs to.
    std::vector<std::pair<std::unique_ptr<StreamMarks>, int>> elim_marks_;
    DeviceBuffer getrf_work_;
};

// ---------------------------------------------------------------------------
// Level start: point table, and the Schur and near-field blocks the level
// transition assembled on the host.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::begin() {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context& ctx = Context::instance();
    ctx.activate();
    // (the exchange arena, for MPI from device memory, came first:
    // begin_operator_build)
    heap_.ensure_initialized();
    // a fresh factorization (unless earlier levels' solve factors are kept, or
    // other operators keep their data, device_heap.hpp); otherwise the parent
    // blocks live on
    if (!adopt_ && device_solve_store().empty() && !other_operators_hold_data()) heap_.reset();
    // (a replicated CA level: its ghosts' factors too, which the device CA
    // solve eliminates again; keyed by Morton index like the local ones)
    keep_solve_ = gpu_sketch_;
    if (keep_solve_) solve_meta_device_.reserve(size_t{1} << 20);  // span lists of the waves
    heap_.reset_peak();
    heap_initialized_ = true;
    copier_ = std::make_unique<HostCopier>(ctx.device());
    if (device_exchange_) exchange_copier_ = std::make_unique<HostCopier>(ctx.device(), 1, 8);
    // a full heap waits for the next background copy, whose block then frees
    heap_.set_reclaimer([this] {
        const auto tw = std::chrono::steady_clock::now();
        bool freed = false;
        for (HostCopier* c : {copier_.get(), exchange_copier_.get()}) {
            if (c == nullptr) continue;
            std::vector<char*> done = c->take_finished();
            freed = freed || !done.empty();
            free_copied_blocks(std::move(done));
        }
        if (!freed) freed = reclaim_released_source();
        if (!freed && copier_) {
            std::vector<char*> done = copier_->wait_released();
            freed = !done.empty();
            free_copied_blocks(std::move(done));
        }
        if (!freed && nonsym_ && keep_solve_ && !keep_failed_ && !kept_.blocks.empty()) {
            // A nonsymmetric level keeps its host copies (skip_xnr_ is off), so
            // its kept factors can go when the level is short: the first solve
            // uploads them from the host instead.  (Its X_RR_full stays with the
            // fill sources: nothing has to be rebuilt from the LU here.)
            drop_kept(false, "the heap is short (nonsymmetric: the host copies serve)");
            keep_failed_ = true;
            keep_failure_ = 2;
            freed = true;
        }
        eliminator_stats().reclaim_wait += std::chrono::duration<double>(std::chrono::steady_clock::now() - tw).count();
        return freed;
    });
    cudaStream_t stream = ctx.stream();
    auto& stats = eliminator_stats();

    // ---- point table: slots are the local boxes' points in box order
    const size_t num_boxes = level_.local_boxes.size();
    states_.assign(num_boxes, BoxState{});
    int64_t num_slots = 0;
    for (size_t b = 0; b < num_boxes; ++b) {
        const Box& box = level_.local_boxes[b];
        if (box.morton_index >= (int64_t{1} << 32)) throw std::runtime_error("LevelEliminator: Morton index exceeds 32 bits");
        if (box.num_points <= 0) {
            if (!occupancy_) throw std::runtime_error("LevelEliminator: empty box " + std::to_string(box.morton_index));
            states_[b].slot = static_cast<int>(num_slots);
            states_[b].eliminated = true;  // empty skeleton
            continue;
        }
        if (!box.far_field_modified_interactions.empty()) {
            throw std::runtime_error("LevelEliminator: stored far-field blocks are not supported (lazy mode expected)");
        }
        states_[b].slot = static_cast<int>(num_slots);
        states_[b].n = static_cast<int>(box.num_points);
        num_slots += box.num_points;
    }
    // CA levels: the ghost boxes (other ranks', eliminated here too), in
    // Morton order, as remote states marked ghost with slots after the local
    // ones.  The boxes that sketch here are the local ones, then the ghosts.
    std::vector<Box*> ghosts;
    if (ca_) {
        std::vector<std::pair<int64_t, size_t>> order;
        for (size_t g = 0; g < level_.ghost_boxes.size(); ++g) order.emplace_back(level_.ghost_boxes[g].morton_index, g);
        std::sort(order.begin(), order.end());
        for (const auto& [morton, g] : order) {
            Box& box = level_.ghost_boxes[g];
            if (box.num_points <= 0) throw std::runtime_error("LevelEliminator: empty ghost box " + std::to_string(morton));
            if (!box.far_field_modified_interactions.empty()) {
                throw std::runtime_error("LevelEliminator: stored far-field blocks are not supported (lazy mode expected)");
            }
            BoxState st;
            st.ghost = true;
            st.slot = static_cast<int>(num_slots);
            st.n = static_cast<int>(box.num_points);
            st.on_boundary = box.on_boundary;
            remote_index_[morton] = remote_states_.size();
            remote_states_.push_back(std::move(st));
            ghosts.push_back(&box);
            num_slots += box.num_points;
        }
    }
    const size_t num_sketching = num_boxes + ghosts.size();
    auto sketching_box = [&](size_t i) -> Box& { return i < num_boxes ? level_.local_boxes[i] : *ghosts[i - num_boxes]; };
    auto sketching_state = [&](size_t i) -> BoxState& {
        return i < num_boxes ? states_[i] : state_of(ghosts[i - num_boxes]->morton_index);
    };
    // static ID training rows of the device sketch (H2_ID_radius > 2,
    // H2_ID_proxy 1), as the host selects them: their points after the local
    // ones, each once, and a list of their slots per box
    std::vector<std::vector<int64_t>> training(num_sketching);
    std::vector<int64_t> training_points;  // global indices, in slot order
    std::vector<int> training_lists;       // per box, slots relative to training_slot0_
    std::vector<size_t> training_at(num_sketching, 0);
    if (gpu_sketch_ && (tree_->id_neighborhood_radius > 2 || tree_->id_proxy_mode == 1)) {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        #pragma omp parallel for schedule(dynamic)
        for (int64_t b = 0; b < static_cast<int64_t>(num_sketching); ++b) {
            if (sketching_state(static_cast<size_t>(b)).n == 0) continue;
            try {
                training[static_cast<size_t>(b)] =
                    select_static_id_training_indices(tree_, &sketching_box(static_cast<size_t>(b)));
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        std::unordered_map<int64_t, int> slot_of;
        for (size_t b = 0; b < num_sketching; ++b) {
            training_at[b] = training_lists.size();
            for (int64_t index : training[b]) {
                auto it = slot_of.emplace(index, static_cast<int>(training_points.size())).first;
                if (it->second == static_cast<int>(training_points.size())) training_points.push_back(index);
                training_lists.push_back(it->second);
            }
        }
    }
    training_slot0_ = static_cast<int>(num_slots);
    num_slots += static_cast<int64_t>(training_points.size());
    if (num_slots >= std::numeric_limits<int>::max()) throw std::runtime_error("LevelEliminator: too many points per rank");
    const size_t dim = static_cast<size_t>(tree_->dimension);
    std::vector<double> xyz(dim * static_cast<size_t>(num_slots));
    std::vector<int64_t> ids(static_cast<size_t>(num_slots));
    for (size_t b = 0; b < num_sketching; ++b) {
        const Box& box = sketching_box(b);
        const size_t slot = static_cast<size_t>(sketching_state(b).slot);
        for (int64_t i = 0; i < box.num_points; ++i) {
            for (size_t d = 0; d < dim; ++d) {
                xyz[dim * (slot + static_cast<size_t>(i)) + d] = static_cast<double>(box.point_coords[dim * i + d]);
            }
            ids[slot + i] = box.point_indices[static_cast<size_t>(i)];
        }
    }
    {
        // (an evaluator that reads only the ids needs no coordinates; the
        // global ones exist for H2_ID_proxy 1)
        const auto& coords = tree_->id_source_point_coords;
        const bool have_coords = coords.size() == static_cast<size_t>(tree_->num_points) * dim;
        for (size_t t = 0; t < training_points.size(); ++t) {
            const size_t slot = static_cast<size_t>(training_slot0_) + t;
            const int64_t index = training_points[t];
            ids[slot] = index;
            for (size_t d = 0; d < dim; ++d) {
                xyz[dim * slot + d] = have_coords ? static_cast<double>(coords[static_cast<size_t>(index) * dim + d]) : 0.0;
            }
        }
    }
    double* d_xyz = heap_.alloc_resident<double>(xyz.size() * sizeof(double), ResidentKind::Points);
    int64_t* d_ids = heap_.alloc_resident<int64_t>(ids.size() * sizeof(int64_t), ResidentKind::Points);
    check_cuda(cudaMemcpyAsync(d_xyz, xyz.data(), xyz.size() * sizeof(double), cudaMemcpyHostToDevice, stream), "upload points");
    check_cuda(cudaMemcpyAsync(d_ids, ids.data(), ids.size() * sizeof(int64_t), cudaMemcpyHostToDevice, stream), "upload ids");
    points_.xyz = d_xyz;
    points_.ids = d_ids;
    points_.dim = tree_->dimension;
    d_xyz_ = d_xyz;  // the table grows when remote boxes arrive (sync_remote_boxes)
    d_ids_ = d_ids;
    num_slots_ = num_slots;
    stats.bytes_up += static_cast<double>(xyz.size() * sizeof(double) + ids.size() * sizeof(int64_t));
    if (!training_lists.empty()) {
        int* d_lists = heap_.alloc_resident<int>(training_lists.size() * sizeof(int), ResidentKind::Lists);
        level_allocs_.push_back(reinterpret_cast<char*>(d_lists));
        check_cuda(cudaMemcpyAsync(d_lists, training_lists.data(), training_lists.size() * sizeof(int),
                                   cudaMemcpyHostToDevice, stream), "upload training lists");
        for (size_t b = 0; b < num_sketching; ++b) {
            BoxState& st = sketching_state(b);
            st.d_training = d_lists + training_at[b];
            st.ntraining = static_cast<int>(training[b].size());
        }
        stats.bytes_up += static_cast<double>(training_lists.size() * sizeof(int));
    }

    // CA levels: the assisting boxes the halo gather listed (points only;
    // their skeletons arrive after the level, refresh_remote_skeletons)
    if (ca_) sync_remote_boxes(false);

    // BPACK_CHECK=kernel: the GPU evaluator against the host's entries
    // (entries through a collective callback: no host entries of its own to compare)
    if (env::check("kernel") && !evaluator_->collective()) check_device_kernel();

    // ---- Schur and near-field blocks: adopted from the device transition,
    //      or packed on the host (in parallel), uploaded once, and scattered
    //      into their own heap allocations with one batched copy
#ifdef H2_HAVE_NVCOMP
    // H2_compress_factors: the transition stored its blocks compressed and
    // handed the raw ones back, so put them back before they are adopted.
    //
    // In three passes, and that matters: alloc_resident can run the reclaimer,
    // which frees other blocks, so no block address may be held across an
    // allocation.  Allocating and collecting in one loop gave an illegal
    // access in the next level's transposes.
    if (adopt_ && !(adopt_->schur_packed.empty() && adopt_->edges_packed.empty())) {
        // The level keeps them packed and materializes each block in the wave
        // that reads it (materialize_adopted), instead of putting the whole
        // set back here: at 192^3 the adopted set is 6.29 GB of a 7.74 GB
        // level, against a per-wave read union of 0.61 GB.  Every block is
        // read once and replaced by the wave that reads it, so its entry goes
        // when it is materialized and the key is free for the fresh block.
        packed_schur_ = std::move(adopt_->schur_packed);
        packed_edges_ = std::move(adopt_->edges_packed);
        packed_strides_ = std::move(adopt_->strides);
        packed_owned_ = std::move(adopt_->packed_owned);
        // the adopted entries' owners, so a materialization can reclaim them
        packed_live_.clear();
        for (auto& [morton, pb] : packed_schur_) {
            if (pb.owner == nullptr && packed_owned_.size() == 1) pb.owner = packed_owned_.front();
            if (pb.owner != nullptr) ++packed_live_[pb.owner];
        }
        for (auto& [key, pb] : packed_edges_) {
            if (pb.owner == nullptr && packed_owned_.size() == 1) pb.owner = packed_owned_.front();
            if (pb.owner != nullptr) ++packed_live_[pb.owner];
        }
        adopt_->schur_packed.clear();
        adopt_->edges_packed.clear();
        adopt_->strides.clear();
        adopt_->packed_owned.clear();
        // The adopted Schur blocks stay packed too, like the near blocks, and
        // go back in the wave that reads them.  They used to be materialized
        // here, all of them: a box's Schur block is reinstalled by several
        // writers - ensure_schur builds it from the kernel when it finds no
        // block, the wave replaces it with the skeleton block after
        // elimination, the exchange installs one - and when each of them
        // treated a null pointer as "there is no block", a packed one was
        // dropped and its Schur-complement contents replaced by kernel values:
        // a 1.2e-2 error in logdet with no crash (50 of 863 blocks at level 3
        // of 96^3).  Since then every such site either names the block to the
        // materialization that precedes it (the wave's gather reads it through
        // its window, the remote-generator pass and build_parent list it) or
        // asserts it is not packed (assert_not_packed), so a missed site is a
        // throw, not a wrong number.  The slot keeps the block's shape with a
        // null pointer; after the box's elimination the skeleton block takes
        // the slot and the stale packed entry is never consulted again, as
        // with a replaced near block (the packer resets it when it next packs
        // the slot).  At 192^3 the eager set was 0.46-0.48 GB of the 1.8 GB
        // of adopted blocks resident at level 4's peak.
    }
#endif
    if (adopt_) {
        // (a CA level: the ghosts' blocks too, received from their owners'
        // devices, h2_gpu/ca_halo.hpp)
        for (size_t b = 0; b < num_sketching; ++b) {
            const Box& box = sketching_box(b);
            BoxState& st = sketching_state(b);
            if (st.n == 0) continue;  // an empty box (occupancy)
            auto it = adopt_->schur.find(box.morton_index);
            if (it == adopt_->schur.end()) {
                // (an owner level: only the boxes eliminated here need their
                // Schur block; the others' are never read)
                if (owner_ != nullptr && !is_writable(box.morton_index)) continue;
                throw std::runtime_error("LevelEliminator::begin: adopted Schur block of box " +
                                         std::to_string(box.morton_index) + " is missing");
            }
            if (it->second.rows != st.n || it->second.cols != st.n) {
                throw std::runtime_error("LevelEliminator::begin: adopted Schur block of box " +
                                         std::to_string(box.morton_index) + " is mis-shaped");
            }
            st.schur = it->second;
            it->second.ptr = nullptr;  // st.schur owns it now, and only it
        }
        if (adopt_->directed != nonsym_) {
            throw std::runtime_error("LevelEliminator::begin: adopted blocks of the other symmetry");
        }
        for (const auto& kv : adopt_->edges) {
            // (rows = the key's low half: the higher box, or a directed block's row box)
            const int64_t lo = static_cast<int64_t>(kv.first >> 32), hi = static_cast<int64_t>(kv.first & 0xffffffffu);
            // a box of another rank is known only from its block here
            if ((is_interest(hi) && kv.second.rows != state_of(hi).n) ||
                (is_interest(lo) && kv.second.cols != state_of(lo).n) || (!is_interest(lo) && !is_interest(hi))) {
                throw std::runtime_error("LevelEliminator::begin: adopted near block is mis-shaped");
            }
        }
        edges_ = std::move(adopt_->edges);
        adopt_.reset();
        check_cuda(cudaStreamSynchronize(stream), "level upload");
        stats.begin += std::chrono::duration<double>(clock::now() - t0).count();
        return;
    }
    struct Placement {
        const MatrixStorage<DataType>* schur = nullptr;      // or
        const ModifiedBlock<DataType>* block = nullptr;      // its view, transposed when `transpose`
        bool transpose = false;
        bool a_sn = false;                                    // (nonsymmetric: its A_SN view)
        DeviceMatrix* dst = nullptr;
        size_t offset = 0;                                    // in elements
    };
    std::vector<Placement> placements;
    size_t image_elems = 0;
    auto place = [&](Placement p, int rows, int cols) {
        p.dst->rows = rows;
        p.dst->cols = cols;
        p.dst->ptr = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * sizeof(S), ResidentKind::Uploaded);
        p.offset = image_elems;
        image_elems += align_up(static_cast<size_t>(rows) * cols, 32);
        placements.push_back(p);
    };
    for (size_t b = 0; b < num_sketching; ++b) {
        Box& box = sketching_box(b);
        BoxState& st = sketching_state(b);
        if (st.n == 0) continue;  // an empty box (occupancy): no blocks
        if (box.schur_complement.is_allocated()) {
            if (box.schur_complement.rows != st.n || box.schur_complement.cols != st.n ||
                box.schur_complement.lda != st.n) {
                throw std::runtime_error("LevelEliminator::begin: Schur block of box " +
                                         std::to_string(box.morton_index) + " is not n x n");
            }
            Placement p;
            p.schur = &box.schur_complement;
            p.dst = &st.schur;
            place(p, st.n, st.n);
        }
        if (nonsym_) {
            // both directed views, untransposed: A_NS[other] = K(other, box),
            // A_SN[other] = K(box, other); a pair of local boxes has each
            // block twice, the first one taken
            auto place_view = [&](int64_t other, const ModifiedBlock<DataType>& block, bool a_sn) {
                const int64_t view_rows = a_sn ? block.A_SN.rows : block.a_ns_rows();
                const int64_t view_cols = a_sn ? block.A_SN.cols : block.a_ns_cols();
                const int64_t other_dim = a_sn ? view_cols : view_rows;
                const int64_t own_dim = a_sn ? view_rows : view_cols;
                const int other_n = is_writable(other) ? state_of(other).n : static_cast<int>(other_dim);
                if (other_dim != other_n || own_dim != st.n) {
                    std::ostringstream oss;
                    oss << "LevelEliminator::begin: near " << (a_sn ? "A_SN" : "A_NS") << " view (" << box.morton_index
                        << ", " << other << ") is " << view_rows << " x " << view_cols;
                    throw std::runtime_error(oss.str());
                }
                const uint64_t key = a_sn ? dkey(box.morton_index, other) : dkey(other, box.morton_index);
                if (edges_.count(key)) return;
                Placement p;
                p.block = &block;
                p.a_sn = a_sn;
                p.dst = &edges_[key];
                place(p, static_cast<int>(view_rows), static_cast<int>(view_cols));
            };
            for (const auto& entry : box.near_field_interaction_map) {
                const auto& block = box.near_field_modified_interactions[static_cast<size_t>(entry.second)];
                if (block.a_ns_is_allocated()) place_view(entry.first, block, false);
            }
            for (const auto& entry : box.near_field_interaction_map_nonsymmetry) {
                const auto& block = box.near_field_modified_interactions[static_cast<size_t>(entry.second)];
                if (block.A_SN.is_allocated()) place_view(entry.first, block, true);
            }
            continue;
        }
        for (const auto& entry : box.near_field_interaction_map) {
            const int64_t other = entry.first;
            const auto& block = box.near_field_modified_interactions[static_cast<size_t>(entry.second)];
            if (!block.a_ns_is_allocated()) continue;
            const uint64_t key = edge_key(box.morton_index, other);
            if (edges_.count(key)) continue;
            // view A_NS of this box for `other`: rows = other's points, cols =
            // this box's (a remote box's size is known only from its block)
            const int other_n = is_writable(other) ? state_of(other).n : static_cast<int>(block.a_ns_rows());
            if (block.a_ns_rows() != other_n || block.a_ns_cols() != st.n) {
                std::ostringstream oss;
                oss << "LevelEliminator::begin: near block (" << box.morton_index << ", " << other
                    << ") is " << block.a_ns_rows() << " x " << block.a_ns_cols() << ", expected "
                    << other_n << " x " << st.n;
                throw std::runtime_error(oss.str());
            }
            Placement p;
            p.block = &block;
            p.transpose = box.morton_index > other;  // physical rows = the higher Morton box
            p.dst = &edges_[key];
            if (p.transpose) place(p, st.n, other_n); else place(p, other_n, st.n);
        }
    }
    if (!placements.empty()) {
        // Pageable image: pinning a multi-GB buffer costs more than the
        // slower pageable copy of a once-per-level upload.
        std::unique_ptr<DataType[]> image_owner(new DataType[image_elems]);
        DataType* h_image = image_owner.get();
        #pragma omp parallel for schedule(dynamic)
        for (int64_t i = 0; i < static_cast<int64_t>(placements.size()); ++i) {
            const Placement& p = placements[static_cast<size_t>(i)];
            DataType* out = h_image + p.offset;
            const int rows = p.dst->rows, cols = p.dst->cols;
            if (p.schur != nullptr) {
                std::memcpy(out, p.schur->data.data(), static_cast<size_t>(rows) * cols * sizeof(DataType));
            } else if (p.a_sn) {
                for (int j = 0; j < cols; ++j)
                    for (int r = 0; r < rows; ++r) out[r + static_cast<size_t>(j) * rows] = p.block->A_SN(r, j);
            } else if (!p.transpose) {
                for (int j = 0; j < cols; ++j)
                    for (int r = 0; r < rows; ++r) out[r + static_cast<size_t>(j) * rows] = p.block->a_ns(r, j);
            } else {
                for (int j = 0; j < cols; ++j)
                    for (int r = 0; r < rows; ++r) out[r + static_cast<size_t>(j) * rows] = p.block->a_ns(j, r);
            }
        }
        S* d_image = heap_.alloc<S>(image_elems * sizeof(S));
        check_cuda(cudaMemcpyAsync(d_image, h_image, image_elems * sizeof(S), cudaMemcpyHostToDevice, stream),
                   "upload blocks");
        std::vector<GatherItem> copies;
        int max_m = 0, max_n = 0;
        for (const Placement& p : placements) {
            copies.push_back(GatherItem{p.dst->ptr, p.dst->rows, p.dst->rows, p.dst->cols, d_image + p.offset, 1,
                                        p.dst->rows, IndexList{}, IndexList{}});
            max_m = std::max(max_m, p.dst->rows);
            max_n = std::max(max_n, p.dst->cols);
        }
        meta_.clear();
        const size_t off = meta_.append(copies);
        char* md = meta_.upload(meta_device_, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off), static_cast<int>(copies.size()), max_m, max_n,
                      md, stream);
        check_cuda(cudaStreamSynchronize(stream), "block upload");
        heap_.free(d_image);
        stats.bytes_up += static_cast<double>(image_elems * sizeof(S));
    }
    check_cuda(cudaStreamSynchronize(stream), "level upload");

    // The device copies are now authoritative; drop the host ones so a stale
    // read fails loudly instead of using outdated values.
    for (size_t b = 0; b < num_sketching; ++b) {
        Box& box = sketching_box(b);
        box.schur_complement = MatrixStorage<DataType>{};
        std::vector<ModifiedBlock<DataType>>().swap(box.near_field_modified_interactions);
        box.near_field_interaction_map.clear();
        box.near_field_interaction_map_nonsymmetry.clear();
    }
    stats.begin += std::chrono::duration<double>(clock::now() - t0).count();
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::eval_blocks(const std::vector<EvalItem>& items, MetaBuilder&,
                                                                    const char* md, size_t offset, int max_m, int max_n,
                                                                    cudaStream_t stream) {
    evaluator_->eval<S>(items, reinterpret_cast<const EvalItem*>(md + offset), max_m, max_n, md, points_, stream);
}

// Entries of the application's device evaluator against its host callback's,
// up to 32 x 32 of each of three blocks: a box with itself (the singular
// self terms), with a local neighbor, and with the last local box.  Prints
// the largest difference relative to the block's largest entry.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::check_device_kernel() {
    cudaStream_t stream = Context::instance().stream();
    int b0 = -1, b2 = -1;
    for (size_t b = 0; b < states_.size(); ++b) {
        if (states_[b].n == 0) continue;
        if (b0 < 0) b0 = static_cast<int>(b);
        b2 = static_cast<int>(b);
    }
    if (b0 < 0) return;
    int b1 = b0;
    for (int64_t m : level_.local_boxes[static_cast<size_t>(b0)].one_hop) {
        if (m != level_.local_boxes[static_cast<size_t>(b0)].morton_index && is_local(m) && state_of(m).n > 0) {
            b1 = static_cast<int>(m - level_.local_morton_start);
            break;
        }
    }
    const int pairs[3][2] = {{b0, b0}, {b0, b1}, {b0, b2}};
    const char* names[3] = {"self", "near", "far"};
    std::vector<EvalItem> evals;
    std::vector<size_t> offset(4, 0);
    int max_m = 0, max_n = 0;
    for (int p = 0; p < 3; ++p) {
        const BoxState& a = states_[static_cast<size_t>(pairs[p][0])];
        const BoxState& c = states_[static_cast<size_t>(pairs[p][1])];
        const int m = std::min(a.n, 32), n = std::min(c.n, 32);
        offset[static_cast<size_t>(p) + 1] = offset[static_cast<size_t>(p)] + static_cast<size_t>(m) * n;
        max_m = std::max(max_m, m);
        max_n = std::max(max_n, n);
    }
    S* d_out = heap_.alloc<S>(offset[3] * sizeof(S));
    for (int p = 0; p < 3; ++p) {
        const BoxState& a = states_[static_cast<size_t>(pairs[p][0])];
        const BoxState& c = states_[static_cast<size_t>(pairs[p][1])];
        evals.push_back(EvalItem{d_out + offset[static_cast<size_t>(p)], std::min(a.n, 32), std::min(a.n, 32),
                                 std::min(c.n, 32), IndexList{-1, a.slot, nullptr}, IndexList{-1, c.slot, nullptr}});
    }
    meta_.clear();
    const size_t off = meta_.append(evals);
    char* md = meta_.upload(meta_device_, stream);
    evaluator_->eval<S>(evals, reinterpret_cast<const EvalItem*>(md + off), max_m, max_n, md, points_, stream);
    std::vector<DataType> dev(offset[3]);
    check_cuda(cudaMemcpyAsync(dev.data(), d_out, offset[3] * sizeof(S), cudaMemcpyDeviceToHost, stream), "kernel check");
    check_cuda(cudaStreamSynchronize(stream), "kernel check");
    heap_.free(d_out);
    std::ostringstream oss;
    oss << "  [gpu] level " << level_index_ << " kernel check (rank " << tree_->mpi_rank << ", "
        << evaluator_->origin() << "): max |K_gpu - K_host| / max |K_host|";
    for (int p = 0; p < 3; ++p) {
        const Box& a = level_.local_boxes[static_cast<size_t>(pairs[p][0])];
        const Box& c = level_.local_boxes[static_cast<size_t>(pairs[p][1])];
        const int m = std::min(static_cast<int>(a.num_points), 32), n = std::min(static_cast<int>(c.num_points), 32);
        std::vector<DataType> host(static_cast<size_t>(m) * n);
        kernel_->evaluate_block_by_index(a.point_indices.data(), m, c.point_indices.data(), n, host.data(), m);
        double diff = 0.0, scale = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                const DataType h = host[static_cast<size_t>(i + j * m)];
                const size_t at = offset[static_cast<size_t>(p)] + static_cast<size_t>(i + j * m);
                diff = std::max(diff, static_cast<double>(std::abs(dev[at] - h)));
                scale = std::max(scale, static_cast<double>(std::abs(h)));
            }
        }
        oss << (p == 0 ? " " : ", ") << names[p] << " " << std::scientific << std::setprecision(2)
            << (scale > 0.0 ? diff / scale : diff);
        oss << " (" << m << "x" << n << ")";
    }
    std::printf("%s\n", oss.str().c_str());
    std::fflush(stdout);
}

// ---------------------------------------------------------------------------
// Host sketch and ID of every box of the wave (the CPU code).  With entries
// through a collective callback the boxes whose sketch asked for blocks not
// fetched yet sketch again after a round (run_entry_phase), as on the host.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::run_ids(
    const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count) {
    std::exception_ptr failure;
    std::mutex failure_mutex;
    int boundary = 0;
    const int team = std::max(1, omp_get_max_threads());
    const int split = split_threads_for(static_cast<int64_t>(wave.size()), team);
    std::vector<int64_t> wave_items(wave.size());
    std::iota(wave_items.begin(), wave_items.end(), int64_t{0});
    run_entry_phase(kernel_entries(kernel_, 0), std::move(wave_items), [&](const std::vector<int64_t>& items) {
    std::vector<int64_t> abandoned;
    #pragma omp parallel default(shared) reduction(+ : boundary)
    {
        FactorizationThreadScratch<CoordType, DataType> scratch;
        scratch.split_threads = split;
        #pragma omp for schedule(dynamic)
        for (int64_t item = 0; item < static_cast<int64_t>(items.size()); ++item) {
            const int64_t bi = items[static_cast<size_t>(item)];
            try {
                Box* box = writable_box(wave[static_cast<size_t>(bi)]);
                if (box == nullptr) throw std::runtime_error("LevelEliminator: wave box is not local");
                entry_item_begin();
                gather_id_target_streamed(tree_, box, level_, kernel_, scratch, box->on_boundary, tolerance_);
                entry_item_end();
                if (!scratch.streamed_sketch_valid) {
                    throw std::runtime_error("LevelEliminator: streamed sketch unavailable for box " +
                                             std::to_string(box->morton_index));
                }
                boundary += box->on_boundary ? 1 : 0;
                IDResult<DataType> id = fmm::compute_id_complex(
                    scratch.sketch_storage.data(), scratch.streamed_sketch_rows, scratch.workspace_cols,
                    scratch.streamed_sketch_rows, tolerance_, 0);
                ensure_nonempty_box_id_rank(id, scratch.workspace_cols);
                box->skeleton_indices = id.skeleton_indices;
                box->redundant_indices = id.redundant_indices;
                box->interpolation_matrix = std::move(id.interpolation);
                if (h2_id_trace_enabled()) {
                    auto wave_it = level_.elimination_wave.find(box->morton_index);
                    h2_id_trace_write(
                        "L" + std::to_string(level_.level) + " m=" + std::to_string(box->morton_index) +
                        " local=" + std::to_string(is_local(box->morton_index) ? 1 : 0) +
                        " ob=" + std::to_string(box->on_boundary) +
                        " n=" + std::to_string(box->num_points) +
                        " k=" + std::to_string(box->skeleton_indices.size()) +
                        " wave=" + std::to_string(wave_it == level_.elimination_wave.end() ? -1 : wave_it->second) +
                        " |" + scratch.id_trace + h2_id_adaptive_trace(scratch.id_adaptive));
                }
                WaveBox& wb = boxes[static_cast<size_t>(bi)];
                wb.box = box;
                wb.n = static_cast<int>(box->num_points);
                wb.k = static_cast<int>(box->skeleton_indices.size());
                wb.r = static_cast<int>(box->redundant_indices.size());
            } catch (const H2EntryMiss&) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                abandoned.push_back(bi);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (entry_item_abandoned()) {
                    abandoned.push_back(bi);
                } else if (!failure) {
                    failure = std::current_exception();
                }
            }
        }
    }
    if (failure) std::rethrow_exception(failure);
    return abandoned;
    });
    boundary_count = boundary;
}

// ---------------------------------------------------------------------------
// Device sketch of the wave's ID targets, then the host ID.
//
// Mirrors gather_id_target_streamed: the rows are the distance-2 ring (full
// or skeleton rows by the same policy), sketched with the same sparse sign
// draws (4 per row, stream seeded with morton + 1).  Every destination row
// sums its contributions in the CPU order (row, then draw), so the kernel
// part of Y matches the host sketch up to kernel rounding.  The lazy fill is
// applied in sketch space as on the host: Y -= W_E P_E per fill source E, in
// the host's source order, with P_E = X_RR_full_E temp2_E[B rows]^T and
// W_E the sketch of E's temp2 rows of the ring blocks that list E.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::sketch_wave(
    const std::vector<int64_t>& wave, std::vector<WaveBox>& boxes, int& boundary_count) {
    using clock = std::chrono::steady_clock;
    EntryOriginScope origin(EntryOrigin::sketch);  // (collective entries' profile)
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    auto t0 = clock::now();

    struct Source {
        int64_t morton = -1;
        const BoxState* state = nullptr;
        int slot_b = 0;             // first temp2 row of this box
        // The temp2^T the readers index: the source's resident block, or, when
        // it is stored compressed, its region of the wave's scratch holding
        // just the slots this wave reads.  slot_b and the runs' src_offset are
        // relative to it either way.
        const S* t2t_base = nullptr;
        int group = 0;              // which group's scratch holds its slots
        std::vector<RunDesc> runs;  // its stored rows, by ring block
        int64_t entries = 0;        // sketched stored rows * sk
        int offset = 0;             // first column of W / row of P in the box's fill buffers
        int* ptr = nullptr;         // device lists (built on the device)
        int* list = nullptr;
        int* row_index = nullptr;   // source row of each sketched row
        int rows = 0;
        // a nonsymmetric matrix: this term fills the K(B, F)^T rows (W from
        // temp4's columns, P = X_RR^T temp2^T[B]), else the K(F, B) ones (W
        // from temp2's rows, P = X_RR temp4[:, B]); p_base is P's right operand
        bool temp4_rows = false;
        const S* p_base = nullptr;
        bool p_trans = false;
    };
    struct Plan {
        int n = 0, d = 0, sk = 0;
        double scale = 0.0;
        int64_t total_rows = 0;
        // a nonsymmetric matrix: the K(F, B) and K(T, B) rows, of the stacked
        // [K(F, B); K(B, F)^T; K(T, B); K(B, T)^T], evaluated as K(row, col)
        int swap_lo0 = 0, swap_hi0 = 0, swap_lo1 = 0, swap_hi1 = 0;
        std::vector<RowBlockDesc> blocks;  // ring blocks: rows and their point slots
        std::vector<Source> sources;
        // device lists (built on the device): draws, row slots, destination lists
        int* draws = nullptr;
        int* rows = nullptr;
        int* ptr = nullptr;
        int* list = nullptr;
        char* lists_block = nullptr;
        std::string trace;
        S* Y = nullptr;        // d x n
        int ldy = 0;           // its leading dimension: d, or the room for the adaptive rows (H2_ID_proxy 2)
        size_t off_y = 0;
        int fill_rank = 0;     // sum of the sources' r, over every group
        // P_E (r_E x n) and W_E (d x r_E) live in the group's shared block,
        // not in a per-box stack: see the fill batches below.
    };
    std::vector<Plan> plans(wave.size());
    const bool trace = h2_id_trace_enabled();
    // the fill sources this wave reads, for the stats (one entry per source,
    // however many of the wave's boxes read it)
    std::mutex src_seen_mutex;
    std::unordered_map<int64_t, size_t> src_seen;
    size_t src_slice_bytes = 0;  // the readers' P slices of them, summed
    // every reader takes one whole slot, so the scratch a wave needs is the
    // union of the (source, slot) pairs it reads, not the sum of its requests
    std::set<std::pair<int64_t, int>> src_union;
    size_t src_union_bytes = 0;
    size_t src_run_bytes = 0;    // ... and their W runs
    size_t src_runs_count = 0;

    // (the stats of the wave's fill sources, once its plans are built)
    struct SrcStats {
        const std::unordered_map<int64_t, size_t>& seen;
        const size_t& slices;
        const size_t& runs;
        const size_t& nruns;
        const size_t& uni;
        std::unordered_set<int64_t>& all;
        ~SrcStats() {
            auto& st = eliminator_stats();
            if (seen.empty()) return;
            size_t bytes = 0;
            for (const auto& [m, b] : seen) {
                bytes += b;
                all.insert(m);
            }
            st.sk_union_bytes_sum += static_cast<double>(uni);
            st.sk_union_bytes_max_wave = std::max(st.sk_union_bytes_max_wave, static_cast<double>(uni));
            st.sk_run_bytes_sum += static_cast<double>(runs);
            st.sk_run_bytes_max_wave = std::max(st.sk_run_bytes_max_wave, static_cast<double>(runs));
            st.sk_runs += static_cast<int64_t>(nruns);
            st.sk_slice_bytes_sum += static_cast<double>(slices);
            st.sk_slice_bytes_max_wave = std::max(st.sk_slice_bytes_max_wave, static_cast<double>(slices));
            st.sk_src_bytes_sum += static_cast<double>(bytes);
            st.sk_src_bytes_max_wave = std::max(st.sk_src_bytes_max_wave, static_cast<double>(bytes));
            st.sk_src_pairs += static_cast<int64_t>(seen.size());
            st.sk_src_distinct = static_cast<int64_t>(all.size());
            ++st.sk_src_waves;
        }
    } src_stats{src_seen, src_slice_bytes, src_run_bytes, src_runs_count, src_union_bytes, sk_src_all_};

    // ---- host plans (index lists only)
    {
        std::exception_ptr failure;
        std::mutex failure_mutex;
        int boundary = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+ : boundary)
        for (int64_t bi = 0; bi < static_cast<int64_t>(wave.size()); ++bi) {
            try {
                Box* box = writable_box(wave[static_cast<size_t>(bi)]);
                if (box == nullptr) throw std::runtime_error("LevelEliminator: wave box is not local");
                // (an owner level: the sketch counts a box as eliminated
                // only in this box's component or a lower colour)
                if (owner_ != nullptr && owner_->set_sketch_context) {
                    owner_->set_sketch_context(owner_->comp_of(box->morton_index));
                }
                boundary += box->on_boundary ? 1 : 0;
                Plan& p = plans[static_cast<size_t>(bi)];
                p.n = static_cast<int>(box->num_points);
                for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
                    if (sketch_box_eliminated(level_, box->one_hop[idx])) box->use_full_set[idx] = 0;
                }

                // ring blocks: row policy and pair fill sources (a box of
                // another rank is an assisting box: skeleton rows once
                // eliminated and its skeleton has arrived, else all points)
                struct Block { int64_t morton; RingInfo info; int row_base; int a; bool skeleton_rows; std::vector<int64_t> sources; };
                std::vector<Block> blocks(box->two_hop.size());
                std::vector<LazyFarSource> pair_sources;
                bool any_sources = false;
                int row_base = 0;
                for (size_t j = 0; j < box->two_hop.size(); ++j) {
                    const int64_t rm = box->two_hop[j];
                    const RingInfo ri = ring_info(rm);
                    collect_lazy_far_sources(level_, box, rm, level_.dimension, pair_sources);
                    Block& blk = blocks[j];
                    blk.morton = rm;
                    blk.info = ri;
                    for (const auto& src : pair_sources) blk.sources.push_back(src.morton);
                    any_sources = any_sources || !pair_sources.empty();
                    const bool both_on_boundary = box->on_boundary && ri.on_boundary;
                    if (ri.remote) {
                        blk.skeleton_rows = sketch_box_eliminated(level_, rm) && !ri.skeleton->empty();
                    } else {
                        blk.skeleton_rows = !(both_on_boundary && pair_sources.empty()) &&
                                            sketch_box_eliminated(level_, rm);
                    }
                    blk.a = blk.skeleton_rows ? static_cast<int>(ri.skeleton->size()) : static_cast<int>(ri.full);
                    if (blk.a == 0) throw std::runtime_error("LevelEliminator: empty ring block");
                    blk.row_base = row_base;
                    const BoxState& rs = state_of(rm);
                    if (blk.skeleton_rows && !rs.skeleton_identity && rs.d_skeleton == nullptr) {
                        throw std::runtime_error("LevelEliminator: ring skeleton not on the device");
                    }
                    p.blocks.push_back(RowBlockDesc{row_base, blk.a, rs.slot,
                                                    blk.skeleton_rows && !rs.skeleton_identity ? rs.d_skeleton : nullptr});
                    row_base += blk.a;
                    if (trace) {
                        char policy = ri.remote ? (blk.skeleton_rows ? 's' : 'a')
                                                : (blk.skeleton_rows ? 'S' : ((both_on_boundary && pair_sources.empty()) ? 'B' : 'F'));
                        p.trace += ' ' + std::to_string(rm) + ':' + std::to_string(blk.a) + policy;
                        for (int64_t m : blk.sources) p.trace += '+' + std::to_string(m);
                    }
                }
                // static training rows after the ring: kernel rows without fill
                // (their draws continue the box's stream, as on the host)
                const int ring_rows = row_base;
                if (nonsym_) {
                    // the ring again, its K(B, F)^T rows (the kernel's own order)
                    const size_t ring_blocks = p.blocks.size();
                    for (size_t j = 0; j < ring_blocks; ++j) {
                        RowBlockDesc rb = p.blocks[j];
                        rb.row_base += ring_rows;
                        p.blocks.push_back(rb);
                    }
                    row_base += ring_rows;
                    p.swap_lo0 = 0;
                    p.swap_hi0 = ring_rows;
                }
                if (const BoxState& bs = state_of(box->morton_index); bs.ntraining > 0) {
                    p.blocks.push_back(RowBlockDesc{row_base, bs.ntraining, training_slot0_, bs.d_training});
                    if (nonsym_) {
                        p.swap_lo1 = row_base;
                        p.swap_hi1 = row_base + bs.ntraining;
                        p.blocks.push_back(RowBlockDesc{row_base + bs.ntraining, bs.ntraining, training_slot0_,
                                                        bs.d_training});
                        row_base += bs.ntraining;
                    }
                    row_base += bs.ntraining;
                }
                p.total_rows = row_base;

                // sketch parameters (compute_id_sparse_sketch's); the draws
                // themselves are generated on the device
                int64_t d = static_cast<int64_t>(std::ceil(1.0 * p.n));
                d = std::min<int64_t>(d, p.total_rows);
                d = std::max<int64_t>(d, p.n);
                p.d = static_cast<int>(d);
                p.sk = std::min(4, p.d);
                p.scale = 1.0 / std::sqrt(static_cast<double>(p.sk));

                // fill sources of the box, in the host order (one_hop order),
                // with the runs of ring rows each one contributes to
                if (any_sources) {
                    LazyFarEndpoint<CoordType> col_end;
                    col_end.morton = box->morton_index;
                    col_end.full_size = p.n;
                    col_end.skeleton = &box->skeleton_indices;
                    col_end.wanted = nullptr;
                    col_end.wanted_count = p.n;
                    std::vector<int64_t> positions;
                    for (int64_t cm : box->one_hop) {
                        if (!sketch_box_eliminated(level_, cm)) continue;
                        if (level_.elimination_wave.find(cm) == level_.elimination_wave.end()) continue;
                        const BoxState& ss = state_of(cm);
                        Box* src_box = ss.remote ? level_.find_generator_box(cm) : held_box(cm);
                        if (src_box == nullptr) {
                            throw std::runtime_error("LevelEliminator: fill source " + std::to_string(cm) + " not found");
                        }
                        if (ss.r == 0 || !has_fill_source(ss)) continue;  // no Schur contribution
                        {
                            std::lock_guard<std::mutex> lock(src_seen_mutex);
                            src_seen.emplace(cm, static_cast<size_t>(ss.ntot) * static_cast<size_t>(ss.r) * sizeof(S));
                        }
                        int64_t slot_offset = 0, slot_count = 0;
                        if (!lazy_far_locate_endpoint_rows(src_box, col_end, slot_offset, slot_count, positions)) continue;
                        if (slot_count != p.n) throw std::runtime_error("LevelEliminator: fill source holds a sliced box");
                        {
                            std::lock_guard<std::mutex> lock(src_seen_mutex);
                            src_slice_bytes += static_cast<size_t>(slot_count) * static_cast<size_t>(ss.r) * sizeof(S);
                            if (src_union.emplace(cm, static_cast<int>(slot_offset)).second) {
                                src_union_bytes += static_cast<size_t>(slot_count) * static_cast<size_t>(ss.r) * sizeof(S);
                            }
                        }
                        Source src;
                        src.morton = cm;
                        src.state = &ss;
                        src.slot_b = static_cast<int>(slot_offset);
                        for (const Block& blk : blocks) {
                            if (std::find(blk.sources.begin(), blk.sources.end(), cm) == blk.sources.end()) continue;
                            LazyFarEndpoint<CoordType> row_end;
                            row_end.morton = blk.morton;
                            row_end.full_size = blk.info.full;
                            row_end.skeleton = blk.info.skeleton;
                            row_end.wanted = blk.skeleton_rows ? blk.info.skeleton : nullptr;
                            row_end.wanted_count = blk.a;
                            int64_t row_offset = 0, row_count = 0;
                            if (!lazy_far_locate_endpoint_rows(src_box, row_end, row_offset, row_count, positions)) continue;
                            // rows of a full-size slot are the wanted positions
                            // (the skeleton for skeleton rows); a skeleton-size
                            // slot is already the skeleton, in order
                            const BoxState& rs = state_of(blk.morton);
                            const bool by_skeleton = row_count == blk.info.full && blk.skeleton_rows &&
                                                     !rs.skeleton_identity;
                            src.runs.push_back(RunDesc{blk.row_base, blk.a, static_cast<int>(row_offset),
                                                       by_skeleton ? rs.d_skeleton : nullptr});
                            src.entries += static_cast<int64_t>(blk.a) * p.sk;
                        }
                        if (src.runs.empty()) continue;
                        {
                            // the W side reads a contiguous run of E's temp2^T
                            // per ring block, on top of the P side's own slot
                            size_t run_rows = 0;
                            for (const RunDesc& rd : src.runs) run_rows += static_cast<size_t>(rd.count);
                            std::lock_guard<std::mutex> lock(src_seen_mutex);
                            src_run_bytes += run_rows * static_cast<size_t>(ss.r) * sizeof(S);
                            src_runs_count += src.runs.size();
                            for (const RunDesc& rd : src.runs) {
                                if (src_union.emplace(cm, rd.src_offset).second) {
                                    src_union_bytes += static_cast<size_t>(rd.count) * static_cast<size_t>(ss.r) * sizeof(S);
                                }
                            }
                        }
                        if (src.runs.size() > static_cast<size_t>(kMaxSourceRuns)) {
                            throw std::runtime_error("LevelEliminator: too many ring blocks per fill source");
                        }
                        if (nonsym_) {
                            // its second term: the same rows of the K(B, F)^T half
                            Source second = src;
                            second.temp4_rows = true;
                            for (RunDesc& rd : second.runs) rd.row_base += ring_rows;
                            p.sources.push_back(std::move(src));
                            p.sources.push_back(std::move(second));
                            continue;
                        }
                        p.sources.push_back(std::move(src));
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        boundary_count = boundary;
        if (owner_ != nullptr && owner_->set_sketch_context) owner_->set_sketch_context(-1);
    }

    // ---- H2_compress_factors: the slots this wave reads out of fill sources
    // that are stored compressed, decompressed into a scratch.
    //
    // Both readers take whole slots - the P side the box's own, the W side one
    // run per ring block - so the scratch holds, per source, just the slots
    // this wave names, laid end to end.  slot_b and each run's src_offset are
    // rewritten to that region, so the row_index the device builds from the
    // runs comes out relative to it and the kernels need no change.  A source
    // still held raw keeps its own pointer and its original offsets.
    char* fill_scratch = nullptr;
    // per group, the chunks to decompress into the scratch before its batches;
    // the pointer arrays of every group go to the device in one upload and
    // each group's plan is built once, so a group costs one batch call
    struct FillGroup {
        std::vector<const void*> src;
        std::vector<void*> dst;
        size_t at = 0;  // its range in the uploaded arrays
#ifdef H2_HAVE_NVCOMP
        bitcompHandle_t plan{};
        bool have_plan = false;
#endif
    };
    std::vector<FillGroup> fill_groups;
    char* fill_ptrs = nullptr;
    const void** fill_d_src = nullptr;
    void** fill_d_dst = nullptr;
#ifdef H2_HAVE_NVCOMP
    {
        // the slots this wave names in each compressed source
        struct Region { const BoxState* st = nullptr; int group = 0; size_t base_rows = 0; std::map<int, size_t> at; };
        std::map<int64_t, Region> regions;
        for (Plan& p : plans) {
            for (Source& src : p.sources) {
                const BoxState& ss = *src.state;
                if (packed_of(ss) == nullptr) continue;  // still raw
                Region& rg = regions[src.morton];
                rg.st = &ss;
                rg.at.emplace(src.slot_b, 0);
                for (const RunDesc& rd : src.runs) rg.at.emplace(rd.src_offset, 0);
            }
        }
        // sources in groups whose slots fit one scratch, which every group
        // then reuses; the regions are relative to the group
        const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
        size_t group_elems = 0, widest = 0;
        int group = 0;
        for (auto& [m, rg] : regions) {
            const BoxState& ss = *rg.st;
            size_t rows = 0;
            for (const auto& [first_row, off] : rg.at) {
                const auto it = ss.t2t_rows.find(first_row);
                if (it == ss.t2t_rows.end()) {
                    throw std::runtime_error("LevelEliminator: fill source slot " + std::to_string(first_row) +
                                             " of box " + std::to_string(m) + " was not compressed");
                }
                rows += static_cast<size_t>(ss.t2t_slots[static_cast<size_t>(it->second)].rows);
            }
            const size_t elems = rows * static_cast<size_t>(ss.r);
            if (group_elems > 0 && (group_elems + elems) * sizeof(S) > kFillWaveScratch) {
                widest = std::max(widest, group_elems);
                ++group;
                group_elems = 0;
            }
            rg.group = group;
            rg.base_rows = group_elems;  // in elements, from the group's start
            size_t at_rows = 0;
            for (auto& [first_row, off] : rg.at) {
                const auto& sl = ss.t2t_slots[static_cast<size_t>(ss.t2t_rows.at(first_row))];
                off = at_rows;
                at_rows += static_cast<size_t>(sl.rows);
            }
            group_elems += elems;
            stats.fs_wave_bytes += static_cast<double>(elems * sizeof(S));
        }
        widest = std::max(widest, group_elems);
        if (widest > 0) {
            fill_scratch = heap_.alloc(widest * sizeof(S), ResidentKind::SketchFill);
            fill_groups.resize(static_cast<size_t>(group) + 1);
            for (auto& [m, rg] : regions) {
                const BoxState& ss = *rg.st;
                S* region = reinterpret_cast<S*>(fill_scratch) + rg.base_rows;
                FillGroup& fg = fill_groups[static_cast<size_t>(rg.group)];
                for (const auto& [first_row, off] : rg.at) {
                    const auto& sl = ss.t2t_slots[static_cast<size_t>(ss.t2t_rows.at(first_row))];
                    const char* csrc = packed_of(ss) + sl.off;
                    if (ss.t2t_packed != nullptr) ++stats.fs_read_dev[1]; else ++stats.fs_read_host[1];
                    S* dst = region + off * static_cast<size_t>(ss.r);
                    for (uint32_t c = 0; c < sl.nchunks; ++c) {
                        fg.src.push_back(csrc);
                        fg.dst.push_back(dst + static_cast<size_t>(c) * chunk_elems);
                        csrc += ss.t2t_strides[sl.chunk0 + c];
                    }
                }
            }
            size_t chunks = 0;
            for (FillGroup& fg : fill_groups) {
                fg.at = chunks;
                chunks += fg.src.size();
            }
            stats.fs_wave_scratch = std::max(stats.fs_wave_scratch, static_cast<double>(widest * sizeof(S)));
            stats.fs_wave_slots += static_cast<int64_t>(chunks);
            stats.fs_wave_groups += static_cast<int64_t>(fill_groups.size());
            if (chunks > 0) {
                std::vector<const void*> h_src;
                std::vector<void*> h_dst;
                h_src.reserve(chunks);
                h_dst.reserve(chunks);
                for (const FillGroup& fg : fill_groups) {
                    h_src.insert(h_src.end(), fg.src.begin(), fg.src.end());
                    h_dst.insert(h_dst.end(), fg.dst.begin(), fg.dst.end());
                }
                const size_t span = align_up(chunks * sizeof(void*));
                fill_ptrs = heap_.alloc(span * 2);
                fill_d_src = reinterpret_cast<const void**>(fill_ptrs);
                fill_d_dst = reinterpret_cast<void**>(fill_ptrs + span);
                check_cuda(cudaMemcpyAsync(fill_d_src, h_src.data(), chunks * sizeof(void*),
                                           cudaMemcpyHostToDevice, stream), "fill slots");
                check_cuda(cudaMemcpyAsync(fill_d_dst, h_dst.data(), chunks * sizeof(void*),
                                           cudaMemcpyHostToDevice, stream), "fill slots");
                check_cuda(cudaStreamSynchronize(stream), "fill slots");  // (once, for every group)
                for (FillGroup& fg : fill_groups) {
                    if (fg.src.empty()) continue;
                    if (bitcompCreateBatchPlanFromCompressedData(&fg.plan, fill_d_src + fg.at, fg.src.size(),
                                                                 nullptr, stream) != BITCOMP_SUCCESS) {
                        throw std::runtime_error("LevelEliminator: could not plan the fill slots' decompression");
                    }
                    fg.have_plan = true;
                }
            }
        }
        // the readers' bases and offsets, relative to their group's scratch
        for (Plan& p : plans) {
            for (Source& src : p.sources) {
                const BoxState& ss = *src.state;
                if (packed_of(ss) == nullptr) {
                    src.t2t_base = ss.temp2t;  // raw, as before
                    // the group of the packed source just below it in Morton
                    // order, so a box's terms run in Morton order whatever
                    // the groups (see the fill passes below)
                    auto below = regions.lower_bound(src.morton);
                    src.group = below == regions.begin() ? 0 : std::prev(below)->second.group;
                    continue;
                }
                const Region& rg = regions.at(src.morton);
                src.t2t_base = reinterpret_cast<const S*>(fill_scratch) + rg.base_rows;
                src.group = rg.group;
                src.slot_b = static_cast<int>(rg.at.at(src.slot_b));
                for (RunDesc& rd : src.runs) rd.src_offset = static_cast<int>(rg.at.at(rd.src_offset));
            }
        }
    }
#else
    for (Plan& p : plans) {
        for (Source& src : p.sources) src.t2t_base = src.state->temp2t;
    }
#endif
    // P's right operand; a nonsymmetric matrix's terms: temp2^T's rows with
    // P = X_RR temp4[:, B], temp4's with P = X_RR^T temp2^T[B]
    for (Plan& p : plans) {
        for (Source& src : p.sources) {
            src.p_base = src.t2t_base;
            if (!nonsym_) continue;
            if (src.state->temp4 == nullptr || src.state->temp2t == nullptr) {
                throw std::runtime_error("LevelEliminator: fill source " + std::to_string(src.morton) +
                                         " without its temp2^T and temp4");
            }
            if (src.temp4_rows) {
                src.t2t_base = src.state->temp4;
                src.p_base = src.state->temp2t;
                src.p_trans = true;
            } else {
                src.p_base = src.state->temp4;
            }
        }
    }

    stats.sketch_plan += std::chrono::duration<double>(clock::now() - t0).count();
    const auto t_gpu = clock::now();

    // ---- device buffers and launch items
    ++wave_stamp_;
    meta_.clear();
    const size_t D = sizeof(S);
    // adaptive ID rows (H2_ID_proxy 2): each sketch heads a taller block that
    // the selected far rows join, and the ID runs in a copy (d_f)
    using Selector = AdaptiveRowSelector<CoordType, DataType>;
    const bool adaptive = tree_->id_proxy_mode == 2;
    size_t y_bytes = 0;
    for (Plan& p : plans) {
        p.off_y = y_bytes;
        p.ldy = adaptive ? Selector::first_capacity(tree_, p.n) : p.d;
        y_bytes = align_up(y_bytes + static_cast<size_t>(p.ldy) * p.n * D);
        for (Source& src : p.sources) {
            src.offset = p.fill_rank;
            p.fill_rank += src.state->r;
        }
    }
    char* d_y = heap_.alloc(std::max<size_t>(y_bytes, 1), ResidentKind::SketchVec);
    char* d_f = adaptive ? heap_.alloc(std::max<size_t>(y_bytes, 1), ResidentKind::SketchVec) : nullptr;
    std::vector<char*> fill_blocks;  // per box [P; W] and lists, freed after the sketch
    std::vector<SketchListsItem> list_items;
    std::vector<SourceListsItem> source_list_items;
    int max_row_blocks = 0;
    for (size_t bi = 0; bi < plans.size(); ++bi) {
        Plan& p = plans[bi];
        const size_t I = sizeof(int);
        if (p.total_rows > kEntryRowMask || p.d > kEntryDestMask) {
            throw std::runtime_error("LevelEliminator: sketch too large for the packed list entries");
        }
        const size_t draws_n = static_cast<size_t>(p.total_rows) * p.sk;
        size_t bytes = 0;
        const size_t o_draws = bytes; bytes = align_up(bytes + draws_n * I);
        const size_t o_rows = bytes;  bytes = align_up(bytes + static_cast<size_t>(p.total_rows) * I);
        const size_t o_ptr = bytes;   bytes = align_up(bytes + (kSketchOwners + 1) * I);
        const size_t o_list = bytes;  bytes = align_up(bytes + draws_n * I);
        std::vector<size_t> o_src(p.sources.size());
        for (size_t j = 0; j < p.sources.size(); ++j) {
            o_src[j] = bytes;
            bytes = align_up(bytes + (kSketchOwners + 1) * I);
            bytes = align_up(bytes + static_cast<size_t>(p.sources[j].entries) * I);
            bytes = align_up(bytes + static_cast<size_t>(p.sources[j].entries / p.sk) * I);  // row_index
        }
        char* block = heap_.alloc(bytes, ResidentKind::SketchLists);
        fill_blocks.push_back(block);
        p.draws = reinterpret_cast<int*>(block + o_draws);
        p.rows = reinterpret_cast<int*>(block + o_rows);
        p.ptr = reinterpret_cast<int*>(block + o_ptr);
        p.list = reinterpret_cast<int*>(block + o_list);
        SketchListsItem li;
        li.seed = static_cast<uint64_t>(wave[bi] + 1);
        li.d = p.d;
        li.sk = p.sk;
        li.rows = static_cast<int>(p.total_rows);
        li.num_blocks = static_cast<int>(p.blocks.size());
        li.blocks_offset = static_cast<int64_t>(meta_.append(p.blocks));
        li.draws = p.draws;
        li.rows_out = p.rows;
        li.ptr = p.ptr;
        li.entries = p.list;
        list_items.push_back(li);
        max_row_blocks = std::max(max_row_blocks, li.num_blocks);
        for (size_t j = 0; j < p.sources.size(); ++j) {
            Source& src = p.sources[j];
            char* at = block + o_src[j];
            src.ptr = reinterpret_cast<int*>(at);
            at += align_up((kSketchOwners + 1) * I);
            src.list = reinterpret_cast<int*>(at);
            at += align_up(static_cast<size_t>(src.entries) * I);
            src.row_index = reinterpret_cast<int*>(at);
            src.rows = static_cast<int>(src.entries / p.sk);
            SourceListsItem si;
            si.box = static_cast<int>(bi);
            si.num_runs = static_cast<int>(src.runs.size());
            si.runs_offset = static_cast<int64_t>(meta_.append(src.runs));
            si.ptr = src.ptr;
            si.entries = src.list;
            si.row_index = src.row_index;
            source_list_items.push_back(si);
        }
    }
    std::vector<OrderedSketchItem> row_items;     // kernel rows -> Y
    // The two batches that read a fill source's temp2^T are split by group: a
    // group's slots are only in the scratch while its batches run.  The rest
    // read W and P, which persist across the groups.
    const size_t num_groups = std::max<size_t>(fill_groups.size(), 1);
    std::vector<std::vector<OrderedSketchItem>> stored_by_group(num_groups);
    std::vector<VBatch<S>> p_by_group(num_groups);
    std::vector<VBatch<S>> p_t_by_group(num_groups);  // (a nonsymmetric matrix's P = X_RR^T temp2^T[B])
    // per group, its fill passes (see below)
    std::vector<std::vector<VBatch<S>>> fill_by_group(num_groups);
    // With the fill sources compressed, one term a source in Morton order
    // (2026-10-11).  The groups are cut by a scratch budget over the whole
    // wave, so a box's sources fell into groups by what else its wave held,
    // and one product a (group, box) summed them in a different order on a
    // ghost's replica than on its owner: the replicated CA levels' ghost
    // copies stopped being their owners' (96^3 CA0, BPACK_CHECK=replica: 1948
    // of 3536 differed, 0 uncompressed), and the rank decisions that flipped
    // with it killed CA0 at the level above its leaf.  A term a source, the
    // terms of a box in Morton order (the groups are Morton ranges, raw
    // sources placed by Morton too), is the same sequence on every rank.  A
    // box's terms in one group are passes: pass j holds every box's j-th.
    const bool per_source_fill = compress_factors_ >= 1;
    int max_d = 0, max_n = 0, max_r = 0;
    // Y -= sum_E W_E P_E is a sum over the box's fill sources, so a group's
    // terms can be accumulated into Y and the scratch reused.  Stacking every
    // source of every box first - one product a box - was 5.41 GB at the
    // 192^3 leaf, 46x the Y it produces, because the whole wave's stacks were
    // live at once.  The groups already exist: a group's temp2^T slots are
    // only in the fill scratch while its batches run, and the same bound now
    // holds for P and W.  The sum changes order, so the result moves by
    // rounding, which is what the host has always done (it applies the
    // sources one by one).
    auto group_of = [&](const Source& src) {
        return static_cast<size_t>(src.group) < num_groups ? static_cast<size_t>(src.group) : size_t{0};
    };
    std::vector<std::vector<int>> g_rank(num_groups, std::vector<int>(plans.size(), 0));
    std::vector<std::vector<size_t>> g_poff(num_groups, std::vector<size_t>(plans.size(), 0));
    std::vector<std::vector<size_t>> g_woff(num_groups, std::vector<size_t>(plans.size(), 0));
    std::vector<std::vector<int>> g_at(num_groups, std::vector<int>(plans.size(), 0));
    size_t pw_elems = 0;  // the widest group's P and W, which every group reuses
    for (size_t g = 0; g < num_groups; ++g) {
        size_t at = 0;
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            const Plan& p = plans[bi];
            int rank = 0;
            for (const Source& src : p.sources) {
                if (group_of(src) == g) rank += src.state->r;
            }
            g_rank[g][bi] = rank;
            if (rank == 0) continue;
            g_poff[g][bi] = at;
            at += static_cast<size_t>(rank) * static_cast<size_t>(p.n);
        }
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            if (g_rank[g][bi] == 0) continue;
            g_woff[g][bi] = at;
            at += static_cast<size_t>(plans[bi].d) * static_cast<size_t>(g_rank[g][bi]);
        }
        pw_elems = std::max(pw_elems, at);
    }
    // X_RR_full of the sources that no longer store it, rebuilt from their LU
    // into a window of this group's own: L and U split out, multiplied by a
    // batched GEMM, the pivots undone on the rows.  The product is magma's and
    // not hand-rolled - the first version was a dot product a thread with
    // stride-ldlu loads and cost 2.1 s a run.  One rebuild a (group, source),
    // not one per reader.
    struct XrrRebuild { size_t gi; const S* lu; const int* piv; int r; size_t off; };
    std::vector<XrrRebuild> xrr_jobs;
    std::map<std::pair<size_t, int64_t>, size_t> xrr_off;  // (group, source) -> its result
    size_t xrr_elems = 0;
    {
        std::vector<size_t> at(num_groups, 0);
        for (const Plan& p : plans) {
            if (p.fill_rank == 0) continue;
            for (const Source& src : p.sources) {
                if (!xrr_rebuilt(*src.state)) continue;
                const size_t gi = group_of(src);
                const auto key = std::make_pair(gi, src.morton);
                if (xrr_off.count(key) != 0) continue;
                const BoxState& ss = *src.state;
                xrr_off.emplace(key, at[gi]);
                xrr_jobs.push_back(XrrRebuild{gi, ss.xrr_lu, ss.xrr_piv, ss.r, at[gi]});
                at[gi] += 3 * static_cast<size_t>(ss.r) * static_cast<size_t>(ss.r);  // result, L, U
            }
        }
        for (size_t g = 0; g < num_groups; ++g) xrr_elems = std::max(xrr_elems, at[g]);
    }
    std::vector<std::vector<LuSplitItemT<S>>> xrr_split(num_groups);
    std::vector<std::vector<RowSwapItemT<S>>> xrr_swap(num_groups);
    std::vector<VBatch<S>> xrr_gemm(num_groups);
    int max_xrr_r = 0;
    S* xrr_win = nullptr;
    if (xrr_elems > 0) {
        char* block = heap_.alloc(xrr_elems * D, ResidentKind::SketchPW);
        fill_blocks.push_back(block);
        xrr_win = reinterpret_cast<S*>(block);
        for (const XrrRebuild& jb : xrr_jobs) {
            const size_t rr = static_cast<size_t>(jb.r) * static_cast<size_t>(jb.r);
            S* dst = xrr_win + jb.off;
            S* l = dst + rr;
            S* u = l + rr;
            xrr_split[jb.gi].push_back(LuSplitItemT<S>{jb.lu, jb.r, l, u, jb.r});
            xrr_gemm[jb.gi].entries.push_back({l, u, dst, jb.r, jb.r, jb.r, jb.r, jb.r, jb.r});
            xrr_swap[jb.gi].push_back(RowSwapItemT<S>{dst, jb.r, jb.r, jb.piv, jb.r});
            max_xrr_r = std::max(max_xrr_r, jb.r);
        }
        eliminator_stats().xr_rebuilt += static_cast<int64_t>(xrr_jobs.size());
    }
    S* pw = nullptr;
    if (pw_elems > 0) {
        char* block = heap_.alloc(pw_elems * D, ResidentKind::SketchPW);
        fill_blocks.push_back(block);
        pw = reinterpret_cast<S*>(block);
    }
    {
        for (Plan& p : plans) {
            p.Y = reinterpret_cast<S*>(d_y + p.off_y);
            OrderedSketchItem item{};
            item.out = p.Y;
            item.ldo = p.ldy;
            item.d = p.d;
            item.ncols = p.n;
            item.rows = static_cast<int>(p.total_rows);
            item.ptr = p.ptr;
            item.entries = p.list;
            item.scale = p.scale;
            item.row_slots = p.rows;
            item.col_base = state_of(wave[static_cast<size_t>(&p - plans.data())]).slot;
            item.swap_lo0 = p.swap_lo0;
            item.swap_hi0 = p.swap_hi0;
            item.swap_lo1 = p.swap_lo1;
            item.swap_hi1 = p.swap_hi1;
            row_items.push_back(item);
            max_d = std::max(max_d, p.d);
            max_n = std::max(max_n, p.n);
            if (p.fill_rank == 0) continue;
            const size_t bi = static_cast<size_t>(&p - plans.data());
            struct Term { size_t g; int64_t morton; int off; int r; };
            std::vector<Term> terms;
            for (Source& src : p.sources) {
                const BoxState& ss = *src.state;
                const int r = ss.r;
                // P_E = X_RR_full temp2[B rows]^T; temp2 is resident row-major (NN)
                const size_t gi = group_of(src);
                const S* xrr = ss.xrr_full;
                if (xrr == nullptr) {  // rebuilt into this group's window
                    const auto it = xrr_off.find(std::make_pair(gi, src.morton));
                    if (it == xrr_off.end()) {
                        throw std::runtime_error("LevelEliminator: fill source " + std::to_string(src.morton) +
                                                 " has no X_RR_full and no LU to rebuild it from");
                    }
                    xrr = xrr_win + it->second;
                }
                const int rank_g = g_rank[gi][bi];
                S* Pg = pw + g_poff[gi][bi];
                S* Wg = pw + g_woff[gi][bi];
                const int off_g = g_at[gi][bi];  // this source's rows of the group's P
                g_at[gi][bi] += r;
                (src.p_trans ? p_t_by_group : p_by_group)[gi].entries.push_back(
                    {xrr, src.p_base + static_cast<size_t>(src.slot_b) * r, Pg + off_g, r, p.n, r, r, r, rank_g});
                OrderedSketchItem w{};
                w.out = Wg + static_cast<size_t>(off_g) * p.d;
                w.ldo = p.d;
                w.d = p.d;
                w.ncols = r;
                w.rows = src.rows;
                w.ptr = src.ptr;
                w.entries = src.list;
                w.scale = p.scale;
                w.src = src.t2t_base;
                w.row_stride = r;
                w.row_index = src.row_index;
                stored_by_group[gi].push_back(w);
                max_r = std::max(max_r, r);
                terms.push_back(Term{gi, src.morton, off_g, r});
            }
            if (per_source_fill) {
                // Y -= W_E P_E a source at a time, in Morton order
                std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
                    return a.g != b.g ? a.g < b.g : a.morton < b.morton;
                });
                size_t pass = 0, cur = num_groups;
                for (const Term& tm : terms) {
                    if (tm.g != cur) {
                        cur = tm.g;
                        pass = 0;
                    }
                    std::vector<VBatch<S>>& passes = fill_by_group[tm.g];
                    if (passes.size() <= pass) passes.resize(pass + 1);
                    passes[pass].entries.push_back({pw + g_woff[tm.g][bi] + static_cast<size_t>(tm.off) * p.d,
                                                    pw + g_poff[tm.g][bi] + tm.off, p.Y, p.d, p.n, tm.r, p.d,
                                                    g_rank[tm.g][bi], p.ldy});
                    ++pass;
                }
            } else {
                // Y -= [W_E1 ...] [P_E1; ...]: uncompressed there is one group,
                // so this is one product a box on every rank
                for (size_t g = 0; g < num_groups; ++g) {
                    const int rank_g = g_rank[g][bi];
                    if (rank_g == 0) continue;
                    if (fill_by_group[g].empty()) fill_by_group[g].resize(1);
                    fill_by_group[g][0].entries.push_back({pw + g_woff[g][bi], pw + g_poff[g][bi], p.Y,
                                                           p.d, p.n, rank_g, p.d, rank_g, p.ldy});
                }
            }
        }
    }
    // device ID: rank, flag, norm and pivots of every box, in one small block
    size_t id_bytes = 0;
    const size_t id_rank = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(int));
    const size_t id_flag = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(int));
    const size_t id_norm = id_bytes; id_bytes = align_up(id_bytes + plans.size() * sizeof(double));
    std::vector<size_t> id_jpvt(plans.size());
    int max_id_n = 0;
    for (size_t bi = 0; bi < plans.size(); ++bi) {
        id_jpvt[bi] = id_bytes;
        id_bytes = align_up(id_bytes + static_cast<size_t>(plans[bi].n) * sizeof(int));
        max_id_n = std::max(max_id_n, plans[bi].n);
    }
    char* d_id = heap_.alloc_resident(std::max<size_t>(id_bytes, 1), ResidentKind::Lists);  // pivots: the level's skeleton lists
    level_allocs_.push_back(d_id);  // its pivots are the device skeleton lists
    std::vector<QrcpItem> id_items(adaptive ? 0 : plans.size());
    std::unique_ptr<Selector> selector;
    if (adaptive) {
        typename Selector::Hooks hooks;
        hooks.place_points = [this](const std::vector<int64_t>& ids, cudaStream_t s) {
            return place_scratch_points(ids, s);
        };
        hooks.eval = [this](const std::vector<EvalItem>& items, MetaBuilder& meta, const char* image, size_t offset,
                            int max_m, int max_n, cudaStream_t s) {
            eval_blocks(items, meta, image, offset, max_m, max_n, s);
        };
        selector = std::make_unique<Selector>(tree_, tolerance_, heap_, std::move(hooks));
        std::vector<typename Selector::Target> targets(plans.size());
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            const Plan& p = plans[bi];
            typename Selector::Target& target = targets[bi];
            target.box = writable_box(wave[bi]);
            target.n = p.n;
            target.cols.base = state_of(wave[bi]).slot;
            target.w = p.Y;
            target.f = reinterpret_cast<S*>(d_f + p.off_y);
            target.cap = p.ldy;
            target.jpvt = reinterpret_cast<int*>(d_id + id_jpvt[bi]);
            target.rank = reinterpret_cast<int*>(d_id + id_rank) + bi;
            target.norm = reinterpret_cast<double*>(d_id + id_norm) + bi;
            target.flag = reinterpret_cast<int*>(d_id + id_flag) + bi;
        }
        selector->start(std::move(targets));
    }
    for (size_t bi = 0; bi < id_items.size(); ++bi) {
        const Plan& p = plans[bi];
        id_items[bi] = QrcpItem{p.Y, p.d, p.n, p.d, reinterpret_cast<int*>(d_id + id_jpvt[bi]),
                                reinterpret_cast<int*>(d_id + id_rank) + bi, reinterpret_cast<double*>(d_id + id_norm) + bi,
                                reinterpret_cast<int*>(d_id + id_flag) + bi};
    }
    const size_t off_list_items = meta_.append(list_items);
    const size_t off_source_list_items = meta_.append(source_list_items);
    const size_t off_rows = meta_.append(row_items);
    std::vector<size_t> off_stored(num_groups, 0);
    for (size_t g = 0; g < num_groups; ++g) off_stored[g] = meta_.append(stored_by_group[g]);
    std::vector<size_t> off_xrr_split(num_groups, 0), off_xrr_swap(num_groups, 0);
    for (size_t g = 0; g < num_groups; ++g) {
        off_xrr_split[g] = meta_.append(xrr_split[g]);
        off_xrr_swap[g] = meta_.append(xrr_swap[g]);
    }
    for (VBatch<S>& b : xrr_gemm) b.stage(meta_);
    const size_t off_id = meta_.append(id_items);
    for (VBatch<S>& b : p_by_group) b.stage(meta_);
    for (VBatch<S>& b : p_t_by_group) b.stage(meta_);
    for (std::vector<VBatch<S>>& passes : fill_by_group) {
        for (VBatch<S>& b : passes) b.stage(meta_);
    }
    if (adaptive) selector->plan(meta_);  // its first round: the IDs of the sketches, the root nodes

    {   // (the sketch's floating-point work, for BPACK_TRACE=wave)
        auto vb = [](const VBatch<S>& b) {
            double f = 0.0;
            for (const auto& e : b.entries) f += 2.0 * e.m * e.n * e.k;
            return f;
        };
        double f = 0.0;
        for (const OrderedSketchItem& it : row_items) f += 2.0 * it.d * static_cast<double>(it.rows) * it.ncols;
        for (const auto& grp : stored_by_group) {
            for (const OrderedSketchItem& it : grp) f += 2.0 * it.d * static_cast<double>(it.rows) * it.ncols;
        }
        for (const VBatch<S>& b : xrr_gemm) f += vb(b);
        for (const VBatch<S>& b : p_by_group) f += vb(b);
        for (const VBatch<S>& b : p_t_by_group) f += vb(b);
        for (const auto& passes : fill_by_group) {
            for (const VBatch<S>& b : passes) f += vb(b);
        }
        double fid = 0.0;
        for (const QrcpItem& it : id_items) {
            const double m = it.m, n = it.n, k = std::min(m, n);
            fid += 4.0 * m * n * k - 4.0 / 3.0 * k * k * k;
        }
        stats.fl_sketch += kFlopScale * f;
        stats.fl_sketch_id += kFlopScale * fid;
    }
    // ---- launches, then Y to the host
    stats.sk_meta += std::chrono::duration<double>(clock::now() - t_gpu).count();
    StreamMarks marks;
    marks.mark(stream);
    char* md = meta_.upload(meta_device_, stream);
    stats.bytes_up += static_cast<double>(meta_.size());
    launch_sketch_lists(reinterpret_cast<const SketchListsItem*>(md + off_list_items),
                        static_cast<int>(list_items.size()), max_row_blocks,
                        reinterpret_cast<const SourceListsItem*>(md + off_source_list_items),
                        static_cast<int>(source_list_items.size()), md, stream);
    marks.mark(stream);
    // (the row slots are the sketch lists' output: read on the device)
    evaluator_->sketch_rows<S>(row_items, reinterpret_cast<const OrderedSketchItem*>(md + off_rows), max_d, max_n,
                               points_, stream);
    marks.mark(stream);
    // H2_compress_factors: a group's slots are decompressed into the scratch,
    // then the batches that read them run, then the next group overwrites it.
    // Everything is on one stream, so that order holds without any waiting.
    for (size_t g = 0; g < num_groups; ++g) {
#ifdef H2_HAVE_NVCOMP
        if (g < fill_groups.size() && fill_groups[g].have_plan) {
            FillGroup& fg = fill_groups[g];
            bitcompSetStream(fg.plan, stream);
            if (bitcompBatchUncompress(fg.plan, fill_d_src + fg.at, fill_d_dst + fg.at) != BITCOMP_SUCCESS) {
                throw std::runtime_error("LevelEliminator: decompressing the fill slots failed");
            }
        }
#endif
        if (!stored_by_group[g].empty()) {
            launch_ordered_sketch_stored(reinterpret_cast<const OrderedSketchItem*>(md + off_stored[g]),
                                         static_cast<int>(stored_by_group[g].size()), max_d, max_r, stream);
        }
        if (!xrr_split[g].empty()) {  // X_RR_full = P L U, before its readers
            launch_lu_split(reinterpret_cast<const LuSplitItemT<S>*>(md + off_xrr_split[g]),
                            static_cast<int>(xrr_split[g].size()), max_xrr_r, stream);
            xrr_gemm[g].gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
            launch_row_swaps(reinterpret_cast<const RowSwapItemT<S>*>(md + off_xrr_swap[g]),
                             static_cast<int>(xrr_swap[g].size()), stream);
        }
        if (!p_by_group[g].entries.empty()) {
            p_by_group[g].gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
        }
        if (!p_t_by_group[g].entries.empty()) {
            p_t_by_group[g].gemm(md, MagmaTrans, MagmaNoTrans, 1.0, 0.0, queue);
        }
        // and collapsed into Y before the next group reuses the scratch (both
        // are on this stream, so the ordering is the stream's)
        for (VBatch<S>& b : fill_by_group[g]) {
            if (!b.entries.empty()) b.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
        }
    }
    marks.mark(stream);
    marks.mark(stream);
    marks.mark(stream);
    if (adaptive) {
        selector->launch(md, stream);
    } else {
        char* d_qrcp_work = nullptr;
        if (const size_t w = qrcp_work_bytes<S>(static_cast<int>(id_items.size()), max_id_n)) {
            d_qrcp_work = heap_.alloc(w, ResidentKind::Adaptive);
        }
        launch_qrcp(reinterpret_cast<const QrcpItem*>(md + off_id), static_cast<int>(id_items.size()), max_id_n,
                    tolerance_, d_qrcp_work, stream);  // (every launch variant gives the same bits)
        heap_.free(d_qrcp_work);  // later launches are ordered after the ID
    }
    marks.mark(stream);
    // written by the SMs into mapped pinned memory: a copy-engine transfer
    // would queue behind the background downloads of the previous wave
    const size_t status_bytes = adaptive ? selector->status_bytes() : 0;
    char* h_status = static_cast<char*>(pinned_pool().result.reserve(std::max<size_t>(status_bytes + id_bytes, 256)));
    char* h_id = h_status + status_bytes;
    if (adaptive) launch_copy_bytes(h_status, selector->status(), status_bytes, stream);
    launch_copy_bytes(h_id, d_id, id_bytes, stream);
    marks.mark(stream);
    {
        // also waits for the previous wave's elimination, queued before
        const auto tw = clock::now();
        check_cuda(cudaStreamSynchronize(stream), "device sketch");
        stats.sk_wait += std::chrono::duration<double>(clock::now() - tw).count();
    }
    // BPACK_TRACE=wave at the top levels: the ID's R diagonal per box, which
    // tells a rank set by the matrix (the diagonal falls through the tolerance)
    // from one set by noise (it sits within a few multiples of it; see the
    // fill store's quantization, 2026-10-11)
    if (!adaptive && level_index_ <= 4 && env::trace("wave")) {
        const int* ranks0 = reinterpret_cast<const int*>(h_id + id_rank);
        int rk = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        auto mag = [](const S& v) -> double {
            if constexpr (std::is_same_v<S, double>) return std::fabs(v);
            else return std::sqrt(v.re * v.re + v.im * v.im);
        };
        std::vector<S> diag;
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            const Plan& p = plans[bi];
            const int K = ranks0[bi];
            if (K <= 0) continue;
            diag.resize(static_cast<size_t>(K));
            check_cuda(cudaMemcpy2D(diag.data(), sizeof(S), p.Y, (static_cast<size_t>(p.d) + 1) * sizeof(S), sizeof(S),
                                    static_cast<size_t>(K), cudaMemcpyDeviceToHost), "ID diagonal");
            const double r00 = mag(diag[0]);
            int below10 = 0, below3 = 0, first_1e1 = -1, first_1e2 = -1;
            double last = 0.0;
            for (int i = 0; i < K; ++i) {
                const double q = r00 > 0.0 ? mag(diag[static_cast<size_t>(i)]) / r00 : 0.0;
                if (q < 10.0 * tolerance_) ++below10;
                if (q < 3.0 * tolerance_) ++below3;
                if (first_1e1 < 0 && q < 0.1) first_1e1 = i;
                if (first_1e2 < 0 && q < 0.01) first_1e2 = i;
                if (i == K - 1) last = q;
            }
            std::printf("  [gpu] level %d rank %d id box %lld: n %d d %d rank %d (%.0f%%) | R_ii/R_00 below 0.1 from %d, below 0.01 from %d, last %.2e"
                        " | kept within 10 tol %d, within 3 tol %d | sources %zu\n", level_index_, rk, static_cast<long long>(wave[bi]),
                        p.n, p.d, K, 100.0 * K / p.n, first_1e1, first_1e2, last, below10, below3, p.sources.size());
        }
        std::fflush(stdout);
    }
    stats.sk_upload += marks.seconds(0);
    stats.sk_rows += marks.seconds(1);
    stats.sk_stored += marks.seconds(2);
    stats.sk_p += marks.seconds(3);
    stats.sk_fill += marks.seconds(4);
    stats.sk_id += marks.seconds(5);
    stats.sk_download += marks.seconds(6);
    stats.bytes_down += static_cast<double>(id_bytes);
    for (char* fill : fill_blocks) heap_.free(fill);
#ifdef H2_HAVE_NVCOMP
    for (FillGroup& fg : fill_groups) {
        if (fg.have_plan) bitcompDestroyPlan(fg.plan);
    }
#endif
    if (fill_ptrs != nullptr) heap_.free(fill_ptrs);
    if (fill_scratch != nullptr) heap_.free(fill_scratch);  // the wave's decompressed slots
    if (adaptive) {
        // the boxes whose root node did not settle them: further rounds, each
        // an image of its own (the sketch's is no longer read)
        const auto ta = clock::now();
        selector->finish(h_status);
        while (selector->active()) {
            ++wave_stamp_;
            meta_.clear();
            selector->plan(meta_);
            md = meta_.upload(meta_device_, stream);
            stats.bytes_up += static_cast<double>(meta_.size());
            selector->launch(md, stream);
            launch_copy_bytes(h_status, selector->status(), status_bytes, stream);
            launch_copy_bytes(h_id, d_id, id_bytes, stream);
            check_cuda(cudaStreamSynchronize(stream), "adaptive ID rows");
            stats.bytes_down += static_cast<double>(status_bytes + id_bytes);
            selector->finish(h_status);
        }
        stats.sk_adaptive += std::chrono::duration<double>(clock::now() - ta).count();
        stats.adaptive_rounds += selector->rounds();
        for (size_t bi = 0; bi < plans.size(); ++bi) {
            if (selector->stats(bi).frames > 1) ++stats.adaptive_multi;
        }
    }
    // d_id stays: its pivot arrays are the device skeleton lists of the boxes
    if (adaptive) {
        heap_.free(d_y);        // the targets (their launches have run)
        wave_sketch_ = d_f;     // T lives here until the elimination has run,
        wave_sketch_extra_ = selector->take_factor_blocks();  // or in the block of a target that grew
    } else {
        wave_sketch_ = d_y;  // T lives here until the elimination has run
    }
    stats.sketch_gpu += std::chrono::duration<double>(clock::now() - t_gpu).count();
    stats.sketch += std::chrono::duration<double>(clock::now() - t0).count();

    // ---- skeletons from the device ID (compute_id_complex's conventions:
    // full rank keeps the natural order; rank 0 keeps column 0 with a zero T)
    t0 = clock::now();
    std::exception_ptr failure;
    std::mutex failure_mutex;
    const int* ranks = reinterpret_cast<const int*>(h_id + id_rank);
    const int* flags = reinterpret_cast<const int*>(h_id + id_flag);
    const double* norms = reinterpret_cast<const double*>(h_id + id_norm);
    #pragma omp parallel for schedule(dynamic, 8)
    for (int64_t bi = 0; bi < static_cast<int64_t>(wave.size()); ++bi) {
        try {
            const Plan& p = plans[static_cast<size_t>(bi)];
            Box* box = writable_box(wave[static_cast<size_t>(bi)]);
            if (flags[bi] == 1) throw std::runtime_error("ID input contains NaN or Inf (box " + std::to_string(box->morton_index) + ")");
            if (flags[bi] == 2) throw std::runtime_error("Triangular solve produced a non-finite interpolation matrix (box " + std::to_string(box->morton_index) + ")");
            const int n = p.n;
            const int K = ranks[bi];
            const int* jpvt = reinterpret_cast<const int*>(h_id + id_jpvt[static_cast<size_t>(bi)]);
            WaveBox& wb = boxes[static_cast<size_t>(bi)];
            box->skeleton_indices.clear();
            box->redundant_indices.clear();
            box->interpolation_matrix = MatrixStorage<DataType>{};  // arrives with the factors
            if (K == n) {
                for (int i = 0; i < n; ++i) box->skeleton_indices.push_back(i);
                box->interpolation_matrix.allocate(n, 0);
            } else if (K == 0) {
                box->skeleton_indices.push_back(0);
                for (int i = 1; i < n; ++i) box->redundant_indices.push_back(i);
            } else {
                box->skeleton_indices.assign(jpvt, jpvt + K);
                box->redundant_indices.assign(jpvt + K, jpvt + n);
            }
            wb.box = box;
            wb.n = n;
            wb.k = static_cast<int>(box->skeleton_indices.size());
            wb.r = static_cast<int>(box->redundant_indices.size());
            // R12 columns, solved in place (adaptive rows: in the ID's block)
            S* factor = adaptive ? const_cast<S*>(selector->id_factor(static_cast<size_t>(bi))) : p.Y;
            wb.ldt = adaptive ? selector->id_ld(static_cast<size_t>(bi)) : p.d;
            wb.T = factor + static_cast<size_t>(wb.k) * wb.ldt;
            wb.state = &state_of(box->morton_index);
            wb.state->d_skeleton = K == n ? nullptr : reinterpret_cast<const int*>(d_id + id_jpvt[static_cast<size_t>(bi)]);
            if (trace) {
                std::ostringstream tail;
                tail << std::setprecision(10) << " | rows=" << p.total_rows << " d=" << p.d
                     << " sketch_norm=" << norms[bi] << " box_sources=";
                for (const Source& src : p.sources) tail << src.morton << ',';
                if (adaptive) tail << h2_id_adaptive_trace(selector->stats(static_cast<size_t>(bi)));
                h2_id_trace_write(
                    "L" + std::to_string(level_.level) + " m=" + std::to_string(box->morton_index) +
                    " local=" + std::to_string(is_local(box->morton_index) ? 1 : 0) +
                    " ob=" + std::to_string(box->on_boundary) + " n=" + std::to_string(box->num_points) +
                    " k=" + std::to_string(wb.k) + " wave=-1 |" + p.trace + tail.str());
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(failure_mutex);
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
    stats.id_loop += std::chrono::duration<double>(clock::now() - t0).count();
    stats.id += std::chrono::duration<double>(clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// One wave.  Returns the number of boundary boxes in it.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
int LevelEliminator<CoordType, DataType, KernelType>::eliminate_wave(
    const std::vector<int64_t>& wave, int wave_index) {
    using clock = std::chrono::steady_clock;
    auto& stats = eliminator_stats();
    // (BPACK_TRACE=wave: the wave's own host time by phase, printed at its end)
    const auto t_wave0 = clock::now();
    const double w_sketch0 = stats.sketch, w_plan0 = stats.plan, w_launch0 = stats.launch, w_store0 = stats.store,
                 w_wait0 = stats.sk_wait;
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    eliminated_any_ = true;
    // a transport must have sent the generators of the previous wave's
    // boxes with remote neighbors (the transport schedule of the loop;
    // CA levels transport nothing)
    close_write_window();
    cur_wave_ = wave_index;
#ifdef H2_HAVE_NVCOMP
    if (wave_index >= 1) measure_near_blocks(stream);
#endif
    if (ca_) last_wave_.clear();
    for (int64_t m : last_wave_) {
        const Box* box = level_.find_local_box(m);
        if (box == nullptr) continue;
        for (int64_t nb : box->one_hop) {
            if (!is_local(nb)) {
                throw std::runtime_error("LevelEliminator: wave " + std::to_string(wave_index) +
                                         " starts before the generators of wave " +
                                         std::to_string(last_wave_index_) + " were sent");
            }
        }
    }
    last_wave_ = wave;  // their generators go out at the next transport
    last_wave_index_ = wave_index;
    since_transport_.insert(since_transport_.end(), wave.begin(), wave.end());
    if (owner_ != nullptr) {
        owner_progress();  // (an owner level: other owners' payloads, this rank's completed sends)
        owner_install_arrived(wave_index / std::max(owner_->num_waves, 1));
    }

    // ---- 1. sketch and ID (device, or host where the device sketch is not supported)
    auto t0 = clock::now();
    std::vector<WaveBox> boxes(wave.size());
    int boundary_count = 0;
    if (gpu_sketch_) {
        sketch_wave(wave, boxes, boundary_count);
    } else {
        run_ids(wave, boxes, boundary_count);
        stats.id += std::chrono::duration<double>(clock::now() - t0).count();
    }
    collect_elimination_marks();
    free_copied_blocks(copier_->take_finished());  // copies that finished during the sketch
    if (exchange_copier_) free_copied_blocks(exchange_copier_->take_finished());
    stats.boxes += static_cast<int64_t>(wave.size());

    // ---- 2. plan: shapes, buffers, launch items
    t0 = clock::now();
    ++wave_stamp_;
    meta_.clear();
    std::unordered_map<int64_t, size_t> wave_index_of;  // morton -> index in boxes
    for (size_t i = 0; i < boxes.size(); ++i) {
        WaveBox& wb = boxes[i];
        wb.state = &state_of(wb.box->morton_index);
        wave_index_of.emplace(wb.box->morton_index, i);
        if (wb.state->eliminated) throw std::runtime_error("LevelEliminator: box eliminated twice");
    }

    // Current neighbor sizes, X_NR row layout and per-box buffer sizes.
    // Buffers of the wave: result (LU, pivots, temp1, X_SR, T; downloaded
    // now), source (temp2 and X_RR_full for the host, in groups of at most
    // kSourceGroup bytes; downloaded while the next wave runs), persist
    // (X_RR_full and temp2^T per box, fill sources of later device sketches;
    // kept to level end), work (per box, freed at wave end).  Only the small
    // result block is one piece: multi-GB blocks would fragment the heap.
    constexpr size_t kSourceGroup = size_t{256} << 20;
    size_t result_bytes = 0, t_bytes = 0;
    // The production half runs in chunks of the wave's boxes (see the loop
    // below), and a source or temp2 group may not straddle two of them, so the
    // chunk a box belongs to is known here, before the groups are formed.
    size_t n_active_boxes = 0;
    for (const WaveBox& wb : boxes) {
        if (wb.r != 0) ++n_active_boxes;
    }
    // ... and only when the wave's outputs would not fit beside what the
    // heap holds: a chunk costs a compression pass (fill sources, kept
    // factors, X_NR^T) with its syncs, 0.9 s a leaf level at 192^3 for four
    size_t wave_bytes = 0;
    for (const WaveBox& wb : boxes) {
        if (wb.r == 0) continue;
        size_t ntot = 0;
        for (int64_t nm : wb.box->one_hop) ntot += static_cast<size_t>(current_size(state_of(nm)));
        wave_bytes += 2 * ntot * static_cast<size_t>(wb.state->n) * sizeof(S);
    }
    const bool chunk_for_room = heap_.used() + wave_bytes > hold_limit();
    const bool chunk_production = compress_factors_ >= 1 && gpu_sketch_ && skip_xnr_ &&
                                  keep_solve_ && !keep_failed_ && chunk_for_room;
    ++stats.np_waves;
    if (chunk_production) ++stats.np_waves_chunked;
    // as many chunks as the wave's bytes need at kWaveChunkBytes each, so a
    // chunk's outputs stay bounded whatever the wave's size (a 56 GB CA wave
    // at 8x the boxes a rank gets 28 chunks, not 4)
    const size_t wave_chunks = chunk_production
        ? std::min(std::max<size_t>(n_active_boxes, 1),
                   std::max<size_t>(1, (wave_bytes + kWaveChunkBytes - 1) / kWaveChunkBytes))
        : 1;
    const size_t chunk_boxes = (std::max<size_t>(n_active_boxes, 1) + wave_chunks - 1) / wave_chunks;
    std::vector<int> group_chunk_of;   // per source group, the chunk it belongs to
    std::vector<size_t> source_group_bytes;
    // temp2 in a block of its own, grouped exactly as the source blocks are:
    // it is the owner pass's operand and the solve factor, and once both are
    // read out of their stores it can go while the source block is still
    // being copied to the host (the copier owns that one).
    std::vector<size_t> temp2_group_bytes;
    // original X_NR of generator boxes, device exchange (packed, not
    // downloaded): a block per chunk, read by the chunk's generator copies
    // only (one block for the wave was 4.27 GB, CA3 at 768^3 on 64 ranks)
    std::vector<size_t> xnr_chunk_bytes;
    const size_t D = sizeof(S);
    std::vector<int> active;  // boxes with a redundant part
    for (size_t i = 0; i < boxes.size(); ++i) {
        WaveBox& wb = boxes[i];
        Box* box = wb.box;
        if (wb.r == 0) {
            // Full rank: no elimination; owner payloads still carry one
            // (empty) row count per one-hop slot.
            box->deferred_xnn_neighbor_point_counts.assign(box->one_hop.size(), 0);
            continue;
        }
        active.push_back(static_cast<int>(i));
        wb.counts.resize(box->one_hop.size());
        wb.row0.resize(box->one_hop.size());
        // (an owner level: the box's state goes to its recipients)
        wb.has_remote = owner_ != nullptr && owner_->recipients.count(box->morton_index) != 0;
        int rows = 0;
        for (size_t a = 0; a < box->one_hop.size(); ++a) {
            const BoxState& ns = state_of(box->one_hop[a]);
            wb.row0[a] = rows;
            wb.counts[a] = current_size(ns);
            rows += current_size(ns);
            wb.has_remote = wb.has_remote || (ns.remote && !ca_);  // (CA levels send no generators)
        }
        wb.ntot = rows;
        const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r);
        wb.off_xrr = result_bytes;      result_bytes = align_up(result_bytes + r * r * D);
        wb.off_piv = result_bytes;      result_bytes = align_up(result_bytes + r * sizeof(int));
        wb.off_temp1 = result_bytes;    result_bytes = align_up(result_bytes + k * r * D);
        wb.off_xrs = result_bytes;      result_bytes = align_up(result_bytes + r * k * D);
        if (nonsym_) {
            wb.off_temp3 = result_bytes; result_bytes = align_up(result_bytes + r * k * D);
        }
        if (gpu_sketch_) {
            wb.off_t = result_bytes;    result_bytes = align_up(result_bytes + k * r * D);
        } else {
            t_bytes = align_up(t_bytes + k * r * D);
        }
    }
    // CA levels: the LUs and the solves for X_RR^{-1} by size class (every
    // copy of a box must factor identically): the active boxes grouped by
    // class, in wave order within
    if (ca_) {
        std::stable_sort(active.begin(), active.end(), [&](int a, int b) {
            return lu_size_class(boxes[static_cast<size_t>(a)].r) < lu_size_class(boxes[static_cast<size_t>(b)].r);
        });
    }
    // The source and temp2 groups, in the order the chunks take the boxes
    // (after the sort: a box's chunk is its place in `active`, and its
    // groups' blocks are allocated by that chunk)
    for (size_t p = 0; p < active.size(); ++p) {
        WaveBox& wb = boxes[static_cast<size_t>(active[p])];
        const size_t r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
        const size_t bytes = align_up(nt * r * D) + align_up(r * r * D) +
                             (wb.has_remote && !device_exchange_ ? align_up(nt * r * D) : 0);
        const int my_chunk = static_cast<int>(p / chunk_boxes);
        // (both blocks bounded: with the device exchange a box adds only its
        // X_RR_full to the source block, so a bound on that alone let one temp2
        // block take an unchunked wave's whole temp2 - 7.3 GB, CA3 at 8x)
        const size_t t2_add = align_up(nt * r * D) * (nonsym_ ? 2 : 1);  // (temp4 beside temp2)
        if (source_group_bytes.empty() || group_chunk_of.back() != my_chunk ||
            (source_group_bytes.back() > 0 && source_group_bytes.back() + bytes > kSourceGroup) ||
            (temp2_group_bytes.back() > 0 && temp2_group_bytes.back() + t2_add > kSourceGroup)) {
            source_group_bytes.push_back(0);
            temp2_group_bytes.push_back(0);
            group_chunk_of.push_back(my_chunk);
        }
        wb.source_group = static_cast<int>(source_group_bytes.size()) - 1;
        size_t& group = source_group_bytes.back();
        size_t& t2 = temp2_group_bytes.back();
        wb.off_temp2 = t2;       t2 = align_up(t2 + nt * r * D);
        if (nonsym_) {
            wb.off_temp4 = t2;   t2 = align_up(t2 + nt * r * D);
        }
        wb.off_xrr_full = group; group = align_up(group + r * r * D);
        if (wb.has_remote && !device_exchange_) {
            wb.off_xnr_orig = group;
            group = align_up(group + nt * r * D);
        } else if (wb.has_remote && !nonsym_) {  // (a nonsymmetric generator sends X_RN, not X_NR)
            if (xnr_chunk_bytes.size() <= static_cast<size_t>(my_chunk)) xnr_chunk_bytes.resize(static_cast<size_t>(my_chunk) + 1, 0);
            size_t& xb = xnr_chunk_bytes[static_cast<size_t>(my_chunk)];
            wb.off_xnr_orig = xb;
            xb = align_up(xb + nt * r * D);
        }
    }
    // (the blocks this wave reads of the previous transition's, for the stats)
    std::set<const void*> wave_schur, wave_edge;
    size_t wave_schur_bytes = 0, wave_edge_bytes = 0;
    const size_t info_offset = result_bytes;
    result_bytes = align_up(result_bytes + std::max<size_t>(active.size(), 1) * sizeof(int));

    char* d_result = heap_.alloc(std::max<size_t>(result_bytes, 1), ResidentKind::WaveResult);
    std::vector<char*> d_sources(source_group_bytes.size());
    for (size_t g = 0; g < d_sources.size(); ++g) d_sources[g] = heap_.alloc(std::max<size_t>(source_group_bytes[g], 1), ResidentKind::WaveSource);
    std::vector<char*> d_temp2(temp2_group_bytes.size(), nullptr);  // allocated by their chunk
    std::vector<char> temp2_to_copier(d_temp2.size(), 0);  // the copier frees these
    char* d_xnr = nullptr;  // (the chunk's, below)
    std::vector<char*> work_blocks;
    std::vector<char*> xnr_blocks;
    char* d_t = heap_.alloc(std::max<size_t>(t_bytes, 1), ResidentKind::WaveWork);
    int* d_info = reinterpret_cast<int*>(d_result + info_offset);

    // ---- The production half, in chunks of the wave's boxes.
    //
    // Everything a box's elimination produces - A_NS, A_NR, X_NR, temp2 and
    // its factors - is dense FP64 that the batched calls write into, and the
    // batch was the whole wave, so all of it was live at once: 4.52 GB of the
    // 192^3 leaf's 6.73 GB transient.  Compressing cannot help, because the
    // values must exist before they can be compressed; so the producer is
    // chunked.  The consumers - the near pass and the owner pass, which read
    // the compressed stores - still run once over the whole wave, which keeps
    // a block's contributions inside one owner pass and so leaves the near
    // store's re-packing where it is.
    bool from_store = compress_factors_ >= 1;  // and every chunk's stores, below
    elim_marks_.clear();
    // Only where nothing outside a chunk reads its raw output: the owner pass
    // has to read the stores (so compression must be on), and the host must
    // not want X_NR, which it does not while the factors are kept.
    const bool can_chunk = chunk_production;
    if (device_exchange_) {  // (the generators' headers and offsets are the wave's; the copies are the chunks')
        if (owner_ != nullptr) {
            plan_owner_generators(boxes, active);
        } else {
            // A full-rank box sends its generator too, without fill (its wave
            // and one-hop list), as the host transports do.  Other ranks mark
            // it eliminated without it - the structured loop by its colour
            // (mark_assisting_eliminated), the unstructured one when its
            // skeleton arrives - and the lazy far fill then wants the wave of
            // every eliminated neighbor (collect_lazy_far_sources), which
            // threw for one that never sent it: the 128K EFIE sphere on 8
            // ranks with H2_compress_factors 1 (unstructured), and the RBF
            // kernel at tol 1e-10 (structured, every leaf box full rank and
            // some at level 3).
            std::vector<int> all(boxes.size());
            for (size_t i = 0; i < boxes.size(); ++i) all[i] = static_cast<int>(i);
            plan_generators(boxes, all);
        }
    }
    const int nchunks = can_chunk ? static_cast<int>(std::min(wave_chunks, std::max<size_t>(active.size(), 1))) : 1;
    for (int chunk = 0; chunk < nchunks; ++chunk) {
        const size_t c0 = std::min(active.size(), static_cast<size_t>(chunk) * chunk_boxes);
        const size_t c1 = nchunks == 1 ? active.size() : std::min(active.size(), c0 + chunk_boxes);
        if (c0 >= c1) continue;
        // (the chunk's planning starts here: the buffers timer measured from
        // the previous chunk's launch, so its compressions and keeps with
        // their syncs were counted as buffers)
        if (chunk > 0) t0 = clock::now();
        const auto t_chunk = clock::now();
        const bool gpu_idle_at_plan = cudaStreamQuery(stream) == cudaSuccess;
        const double heap_s0 = heap_.op_seconds();
        const size_t heap_n0 = heap_.alloc_calls() + heap_.free_calls();
        const std::vector<int> chunk_active(active.begin() + static_cast<ptrdiff_t>(c0),
                                            active.begin() + static_cast<ptrdiff_t>(c1));
        elim_marks_.emplace_back(std::make_unique<StreamMarks>(), 0);
        StreamMarks& marks = *elim_marks_.back().first;
        // this chunk's temp2 blocks; its X_NR blocks come with its boxes
        for (size_t g = 0; g < d_temp2.size(); ++g) {
            if (group_chunk_of[g] != chunk) continue;
            d_temp2[g] = heap_.alloc(std::max<size_t>(temp2_group_bytes[g], 1), ResidentKind::WaveOutput);
        }
        if (static_cast<size_t>(chunk) < xnr_chunk_bytes.size() && xnr_chunk_bytes[static_cast<size_t>(chunk)] > 0) {
            d_xnr = heap_.alloc(xnr_chunk_bytes[static_cast<size_t>(chunk)], ResidentKind::WaveSource);
        }
        const size_t xnr_first = xnr_blocks.size();
        const size_t work_first = work_blocks.size();
        std::vector<AddStoreItem> stores;   // step 5's, for this chunk
        int max_store_m = 0, max_store_n = 0;
        // ... and the blocks they write, to pack once written.  All the
        // store's packer reads of a target is these four fields.
        struct FreshNear { DeviceMatrix* block; uint64_t key; int rows, cols; };
        std::vector<FreshNear> fresh_near;
        // Host image of all interpolation matrices, laid out as d_t - only
        // on the host-sketch path.  Sized for the whole wave and allocated
        // once per production chunk, it was a 60 MB zero-filled vector 128
        // times a leaf level with the device sketch never reading it: the
        // 1.3 s of "plan boxes [buffers]" at 192^3 (2026-10-10).
        std::vector<char> t_image(gpu_sketch_ ? 0 : t_bytes);
        {
            size_t t_off = 0;
            for (int i : chunk_active) {
                WaveBox& wb = boxes[static_cast<size_t>(i)];
                const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
                if (!gpu_sketch_) {
                    const auto& T = wb.box->interpolation_matrix;
                    if (T.rows != wb.k || T.cols != wb.r || T.lda != wb.k) throw std::runtime_error("LevelEliminator: unexpected T shape");
                    std::memcpy(t_image.data() + t_off, T.data.data(), k * r * D);
                    wb.T = reinterpret_cast<S*>(d_t + t_off);
                    wb.ldt = wb.k;
                    t_off = align_up(t_off + k * r * D);
                }
                wb.xrr = reinterpret_cast<S*>(d_result + wb.off_xrr);
                wb.piv = reinterpret_cast<int*>(d_result + wb.off_piv);
                wb.temp1 = reinterpret_cast<S*>(d_result + wb.off_temp1);

                char* source = d_sources[static_cast<size_t>(wb.source_group)];
                wb.temp2 = reinterpret_cast<S*>(d_temp2[static_cast<size_t>(wb.source_group)] + wb.off_temp2);
                wb.xrr_full = reinterpret_cast<S*>(source + wb.off_xrr_full);
                if (gpu_sketch_) {
                    // X_RR_full gets its block now.  temp2^T only when its
                    // sub-group is transposed and compressed near the end of the
                    // wave, so a wave never holds more than one group of it raw
                    // (2.26 GB a wave at the Color leaf against about 0.1 GB a
                    // group) - but that only works because each group is freed
                    // once packed.  Without compression nothing frees it, so it
                    // has to be allocated up front as it always was.
                    const size_t xrr_bytes = align_up(r * r * D);
                    char* persist = heap_.alloc_resident(xrr_bytes, ResidentKind::FillXrr);
                    char* persist_t2t = nullptr;
                    size_t t2t_bytes = 0;
                    if (compress_factors_ < 1) {
                        t2t_bytes = std::max<size_t>(nt * r * D, 1);
                        persist_t2t = heap_.alloc_resident(t2t_bytes, ResidentKind::FillRaw);
                        wb.persist_temp2t = reinterpret_cast<S*>(persist_t2t);
                    }
                    {
                        auto& es = eliminator_stats();
                        es.persist_bytes += static_cast<double>(xrr_bytes + t2t_bytes);
                        es.persist_xrr_bytes += static_cast<double>(xrr_bytes);
                        es.persist_peak = std::max(es.persist_peak, es.persist_bytes);
                        ++es.persist_blocks;
                    }
                    keep_fill_source(*wb.state, persist, xrr_bytes, persist_t2t, t2t_bytes, wb.box->one_hop);
                    wb.persist_xrr = reinterpret_cast<S*>(persist);
                    if (nonsym_) {
                        // temp4, the fill's column side, kept as temp2^T is: to the
                        // level's end, or (several ranks) until its readers are done
                        const size_t t4_bytes = std::max<size_t>(nt * r * D, 1);
                        char* persist_t4 = heap_.alloc_resident(t4_bytes, ResidentKind::FillRaw);
                        if (early_free_sources_) {
                            wb.state->persist_t4 = persist_t4;
                            wb.state->persist_t4_bytes = t4_bytes;
                        } else {
                            level_allocs_.push_back(persist_t4);
                        }
                        wb.persist_temp4 = reinterpret_cast<S*>(persist_t4);
                        eliminator_stats().persist_bytes += static_cast<double>(t4_bytes);
                    }
                }
                size_t w_off = 0;
                const size_t o_tmp1 = w_off; w_off = align_up(w_off + r * r * D);
                const size_t o_tmp2 = w_off; w_off = align_up(w_off + k * r * D);
                const size_t o_xsr = w_off;  w_off = align_up(w_off + k * r * D);
                const size_t o_ans = w_off;  w_off = align_up(w_off + nt * k * D);
                const size_t o_xns = w_off;  w_off = align_up(w_off + nt * k * D);
                size_t o_xrs = 0, o_asn = 0, o_xsn = 0, o_xrsp = 0, o_xrnp = 0;
                if (nonsym_) {
                    o_xrs = w_off;  w_off = align_up(w_off + r * k * D);
                    o_asn = w_off;  w_off = align_up(w_off + k * nt * D);
                    o_xsn = w_off;  w_off = align_up(w_off + k * nt * D);
                    o_xrsp = w_off; w_off = align_up(w_off + r * k * D);
                    o_xrnp = w_off; w_off = align_up(w_off + r * nt * D);
                }
                char* work = heap_.alloc(w_off, ResidentKind::WaveWork);
                work_blocks.push_back(work);
                if (nonsym_) {
                    wb.xrs = reinterpret_cast<S*>(work + o_xrs);
                    wb.asn = reinterpret_cast<S*>(work + o_asn);
                    wb.xsn = reinterpret_cast<S*>(work + o_xsn);
                    wb.xrs_p = reinterpret_cast<S*>(work + o_xrsp);
                    wb.xrn_p = reinterpret_cast<S*>(work + o_xrnp);
                    wb.temp3 = reinterpret_cast<S*>(d_result + wb.off_temp3);
                    wb.temp4 = reinterpret_cast<S*>(d_temp2[static_cast<size_t>(wb.source_group)] + wb.off_temp4);
                    // X_RN, the owner pass's operand: its own block, as X_NR
                    char* xrn_block = heap_.alloc(std::max<size_t>(nt * r * D, 1), ResidentKind::WaveOutput);
                    xnr_blocks.push_back(xrn_block);
                    wb.xrn = reinterpret_cast<S*>(xrn_block);
                }
                // X_NR likewise: the near pass reads A_NS and X_NS out of the work
                // block after X_NR is in its store, so X_NR cannot share with them
                char* xnr_block = heap_.alloc(std::max<size_t>(nt * r * D, 1), ResidentKind::WaveOutput);
                xnr_blocks.push_back(xnr_block);
                wb.tmp1 = reinterpret_cast<S*>(work + o_tmp1);
                wb.tmp2 = reinterpret_cast<S*>(work + o_tmp2);
                wb.xsr = reinterpret_cast<S*>(work + o_xsr);
                wb.ans = reinterpret_cast<S*>(work + o_ans);
                wb.xnr = reinterpret_cast<S*>(xnr_block);
                wb.xns = reinterpret_cast<S*>(work + o_xns);
                wb.s = heap_.alloc_resident<S>(std::max<size_t>(k * k * D, 1), ResidentKind::Near);  // becomes the Schur block
            }
        }

        {   // the wave's transient by array (unaligned sums, so a little under the
            // allocations themselves)
            double w_ans = 0, w_xns = 0, w_xnr = 0, w_temp2 = 0, w_orig = 0, w_res = 0, w_t = 0, w_small = 0;
            for (int i : chunk_active) {
                const WaveBox& wb = boxes[static_cast<size_t>(i)];
                const double k = wb.k, r = wb.r, nt = wb.ntot, d = static_cast<double>(D);
                w_ans += nt * k * d;                    // A_NS, gathered neighbour rows
                w_xns += nt * k * d;                    // the X_NS update written back
                w_xnr += nt * r * d;                    // X_NR = A_NR - A_NS T
                w_temp2 += nt * r * d;                  // temp2 = -X_NR X_RR^-1
                if (wb.has_remote) w_orig += nt * r * d;  // the original X_NR for the exchange
                w_small += r * r * d + 2 * k * r * d;   // tmp1, tmp2, xsr
                w_small += r * r * d;                   // X_RR_full
                w_res += r * r * d + r * 4 + 2 * k * r * d;  // X_RR, pivots, temp1, X_RS
                if (gpu_sketch_) w_res += k * r * d; else w_t += k * r * d;
            }
            const double total = w_ans + w_xns + w_xnr + w_temp2 + w_orig + w_res + w_t + w_small;
            if (total > stats.tw_max) {
                stats.tw_max = total;
                stats.tw_ans = w_ans; stats.tw_xns = w_xns; stats.tw_xnr = w_xnr;
                stats.tw_temp2 = w_temp2; stats.tw_xnr_orig = w_orig;
                stats.tw_result = w_res; stats.tw_t = w_t; stats.tw_small = w_small;
            }
        }
        stats.pb_buffers += std::chrono::duration<double>(clock::now() - t0).count();
        stats.pb_heap += heap_.op_seconds() - heap_s0;
        stats.pb_heap_ops += static_cast<int64_t>(heap_.alloc_calls() + heap_.free_calls() - heap_n0);
        const auto t_items = clock::now();
        // Launch items of the box region.
        std::vector<EvalItem> evals;
        std::vector<GatherItem> fills;
        std::vector<GatherItem> copy_full;   // X_RR_full <- X_RR
        std::vector<IdentityItem> identities;  // X_RR^{-1} = I U^{-1} L^{-1}, before the row swaps
        std::vector<IdentityItem> symmetrizes;  // X_RR := its symmetric part, before the LU (as the host)
        std::vector<GatherItem> copy_xnr;      // original X_NR of generator boxes, for the host
        int max_xnr_m = 0, max_xnr_n = 0;
        std::vector<SymAddItem> sym_adds;
        std::vector<ColumnSwapItem> swaps;
        // (stores: declared before the chunks, filled by the near pass)
        std::vector<TransposeItem> transposes;  // temp2^T into the persist block
        int max_tr_m = 0, max_tr_n = 0;
        VBatch<S> g1, g2, g3, g4, g5, g6, g7, trsm, getrf, gsolve;
        // a nonsymmetric matrix: X_RR -= A_RS T (NN), X_RR -= T^T A_SR (TN),
        // X_RS and X_RN -= T^T (A_SS, A_SN) (TN), S += temp1 X_RS (NN), the X_NS
        // and X_SN updates (NN); P^T on the copies of X_RS and X_RN; temp4 into
        // its fill source
        VBatch<S> gn_rs_t, gn_t_sr, gn_t_rows, g6n, g7n;
        std::vector<RowSwapItemT<S>> row_swaps_fwd;
        std::vector<GatherItem> copy_t4;
        int max_t4_m = 0, max_t4_n = 0;
        int max_eval_m = 0, max_eval_n = 0, max_fill_m = 0, max_fill_n = 0, max_copy_m = 0, max_copy_n = 0;
        int max_r = 0;
        std::vector<DeviceMatrix> retired;  // blocks replaced by this wave

        auto add_eval = [&](S* out, int ld, int m, int n, IndexList rows, IndexList cols) {
            if (m <= 0 || n <= 0) return;
            evals.push_back(EvalItem{out, ld, m, n, rows, cols});
            max_eval_m = std::max(max_eval_m, m);
            max_eval_n = std::max(max_eval_n, n);
        };
        auto add_gather = [&](std::vector<GatherItem>& list, int& mm, int& mn, S* out, int ld, int m, int n,
                              const S* src, int64_t rs, int64_t cs, IndexList rows, IndexList cols) {
            if (m <= 0 || n <= 0) return;
            list.push_back(GatherItem{out, ld, m, n, src, rs, cs, rows, cols});
            mm = std::max(mm, m);
            mn = std::max(mn, n);
        };
        auto push = [](VBatch<S>& batch, const S* a, int lda, const S* b, int ldb, S* c, int ldc,
                       int m, int n, int k) {
            if (m <= 0 || n <= 0 || k <= 0) return;
            batch.entries.push_back({a, b, c, m, n, k, std::max(lda, 1), std::max(ldb, 1), std::max(ldc, 1)});
        };
        const IndexList identity{};

    #ifdef H2_HAVE_NVCOMP
        // H2_compress_factors: the blocks this wave's gather reads come back into
        // a bounded window, a pass at a time, rather than all at once - the
        // unbounded version left 6.16 GB of decompressed copies resident at the
        // 192^3 leaf, which is why holding the near blocks packed saved nothing.
        // A block's offset in the window is chosen here, before the items are
        // staged, because a GatherItem's src is baked into the uploaded image and
        // cannot be patched later.
        char* fill_window = nullptr;
        std::vector<std::vector<WinBlock>> win_passes;
        std::unordered_map<const void*, std::pair<size_t, size_t>> win_of;   // block -> (pass, offset)
        size_t win_at = 0;
        size_t fill_window_bytes = kFillWindowBytes;
        if (kHoldNearBlocksPacked && (!packed_schur_.empty() || !packed_edges_.empty())) {
            // halved until it fits: a tight heap gets more passes, not a throw
            while ((fill_window = heap_.try_alloc(fill_window_bytes, ResidentKind::SrcWindow)) == nullptr &&
                   fill_window_bytes > (size_t{64} << 20)) {
                fill_window_bytes /= 2;
            }
            if (fill_window != nullptr) win_passes.emplace_back();
            else ++stats.fw_no_window;
        }
        std::vector<size_t> fill_pass;   // the pass each fills item belongs to
        auto gather_src = [&](DeviceMatrix& m, bool schur, uint64_t key) -> const S* {
            if (m.ptr != nullptr) return m.ptr;
            if (!(schur ? is_packed_schur(static_cast<int64_t>(key)) : is_packed_edge(key))) return nullptr;
            if (fill_window == nullptr) {
                // no window could be allocated (a tight heap): the block put
                // back on its own, resident, as the owner pass does
                std::vector<AdoptedBlock> one{AdoptedBlock{schur, key, &m}};
                materialize_adopted(one, stream, MatOwnerJoin);
                ++stats.fw_fallbacks;
                return m.ptr;
            }
            const auto it = win_of.find(&m);
            if (it != win_of.end()) return reinterpret_cast<const S*>(fill_window + it->second.second);
            const size_t bytes = align_up(static_cast<size_t>(m.rows) * m.cols * sizeof(S));
            if (bytes > fill_window_bytes) return nullptr;      // cannot be windowed; left to the guards
            if (win_at + bytes > fill_window_bytes) {
                win_passes.emplace_back();
                win_at = 0;
            }
            const size_t off = win_at;
            win_at += bytes;
            win_passes.back().push_back(WinBlock{schur, key, off, m.rows, m.cols});
            win_of.emplace(&m, std::make_pair(win_passes.size() - 1, off));
            return reinterpret_cast<const S*>(fill_window + off);
        };
        auto pass_of = [&](const DeviceMatrix& m) {
            const auto it = win_of.find(&m);
            return it == win_of.end() ? size_t{0} : it->second.first;
        };
        auto tag_fills = [&](size_t pass) {
            while (fill_pass.size() < fills.size()) fill_pass.push_back(pass);
        };
    #endif

        for (int i : chunk_active) {
            WaveBox& wb = boxes[static_cast<size_t>(i)];
            Box* box = wb.box;
            BoxState& st = *wb.state;
            const int n = wb.n, k = wb.k, r = wb.r, nt = wb.ntot;
            std::vector<int> skel(box->skeleton_indices.begin(), box->skeleton_indices.end());
            std::vector<int> red(box->redundant_indices.begin(), box->redundant_indices.end());
            wb.skeleton_offset = static_cast<int64_t>(meta_.append(skel));
            wb.redundant_offset = static_cast<int64_t>(meta_.append(red));
            const IndexList skel_local{wb.skeleton_offset, 0}, red_local{wb.redundant_offset, 0};
            const IndexList skel_slots{wb.skeleton_offset, st.slot}, red_slots{wb.redundant_offset, st.slot};

            // A_RR, A_SR, A_SS from the Schur block or the kernel
            const S* schur_src = gather_src(st.schur, true, static_cast<uint64_t>(box->morton_index));
            if (schur_src != nullptr) {
                if (wave_schur.insert(schur_src).second) {
                    wave_schur_bytes += static_cast<size_t>(st.schur.rows) * st.schur.cols * D;
                }
                if (st.schur.rows != n || st.schur.cols != n) throw std::runtime_error("LevelEliminator: Schur block shape");
                add_gather(fills, max_fill_m, max_fill_n, wb.xrr, r, r, r, schur_src, 1, n, red_local, red_local);
                add_gather(fills, max_fill_m, max_fill_n, wb.xsr, k, k, r, schur_src, 1, n, skel_local, red_local);
                add_gather(fills, max_fill_m, max_fill_n, wb.s, k, k, k, schur_src, 1, n, skel_local, skel_local);
                if (nonsym_) add_gather(fills, max_fill_m, max_fill_n, wb.xrs, r, r, k, schur_src, 1, n, red_local, skel_local);
                tag_fills(pass_of(st.schur));
                retired.push_back(st.schur);
            } else {
                assert_not_packed(true, static_cast<uint64_t>(box->morton_index), "the wave's Schur gather");
                add_eval(wb.xrr, r, r, r, red_slots, red_slots);
                add_eval(wb.xsr, k, k, r, skel_slots, red_slots);
                add_eval(wb.s, k, k, k, skel_slots, skel_slots);
                if (nonsym_) add_eval(wb.xrs, r, r, k, red_slots, skel_slots);
            }

            // A_NS / A_NR of every neighbor: stored block or kernel
            for (size_t a = 0; a < box->one_hop.size(); ++a) {
                const int64_t nm = box->one_hop[a];
                BoxState& ns = state_of(nm);
                const int c = static_cast<int>(wb.counts[a]);
                if (c == 0) continue;
                if (nonsym_) {
                    // A_NS, A_NR from K(N, B) (c x n); A_SN, A_RN from K(B, N) (n x c)
                    const IndexList nrows = current_slots(ns);
                    auto it_ns = edges_.find(dkey(nm, box->morton_index));
                    if (it_ns != edges_.end() && it_ns->second.ptr != nullptr) {
                        const DeviceMatrix& e = it_ns->second;
                        if (e.rows != c || e.cols != n) throw std::runtime_error("LevelEliminator: K(N, B) block shape");
                        add_gather(fills, max_fill_m, max_fill_n, wb.ans + wb.row0[a], nt, c, k, e.ptr, 1, e.rows,
                                   identity, skel_local);
                        add_gather(fills, max_fill_m, max_fill_n, wb.xnr + wb.row0[a], nt, c, r, e.ptr, 1, e.rows,
                                   identity, red_local);
                        retired.push_back(e);
                    } else {
                        add_eval(wb.ans + wb.row0[a], nt, c, k, nrows, skel_slots);
                        add_eval(wb.xnr + wb.row0[a], nt, c, r, nrows, red_slots);
                    }
                    auto it_sn = edges_.find(dkey(box->morton_index, nm));
                    S* asn_cols = wb.asn + static_cast<size_t>(wb.row0[a]) * k;
                    S* arn_cols = wb.xrn + static_cast<size_t>(wb.row0[a]) * r;
                    if (it_sn != edges_.end() && it_sn->second.ptr != nullptr) {
                        const DeviceMatrix& e = it_sn->second;
                        if (e.rows != n || e.cols != c) throw std::runtime_error("LevelEliminator: K(B, N) block shape");
                        add_gather(fills, max_fill_m, max_fill_n, asn_cols, k, k, c, e.ptr, 1, e.rows, skel_local,
                                   identity);
                        add_gather(fills, max_fill_m, max_fill_n, arn_cols, r, r, c, e.ptr, 1, e.rows, red_local,
                                   identity);
                        retired.push_back(e);
                    } else {
                        add_eval(asn_cols, k, k, c, skel_slots, nrows);
                        add_eval(arn_cols, r, r, c, red_slots, nrows);
                    }
                    continue;
                }
                S* ans_rows = wb.ans + wb.row0[a];
                S* anr_rows = wb.xnr + wb.row0[a];
                auto it = edges_.find(edge_key(box->morton_index, nm));
                const S* edge_src = it == edges_.end()
                                        ? nullptr
                                        : gather_src(it->second, false, edge_key(box->morton_index, nm));
                if (it != edges_.end() && edge_src != nullptr) {
                    const DeviceMatrix& e = it->second;
                    if (wave_edge.insert(edge_src).second) {
                        wave_edge_bytes += static_cast<size_t>(e.rows) * e.cols * D;
                    }
                    int64_t rs = 1, cs = e.rows;  // view rows = neighbor, cols = this box
                    if (box->morton_index < nm) {
                        if (e.rows != c || e.cols != n) {
                            throw std::runtime_error("LevelEliminator: near block shape (lo): block " + std::to_string(e.rows) +
                                                     " x " + std::to_string(e.cols) + ", box " + std::to_string(box->morton_index) +
                                                     " n " + std::to_string(n) + ", neighbor " + std::to_string(nm) + " c " +
                                                     std::to_string(c) + (e.ptr == nullptr ? " (windowed)" : " (raw)") +
                                                     (is_packed_edge(edge_key(box->morton_index, nm)) ? ", packed entry live" : ""));
                        }
                    } else {
                        if (e.rows != n || e.cols != c) {
                            throw std::runtime_error("LevelEliminator: near block shape (hi): block " + std::to_string(e.rows) +
                                                     " x " + std::to_string(e.cols) + ", box " + std::to_string(box->morton_index) +
                                                     " n " + std::to_string(n) + ", neighbor " + std::to_string(nm) + " c " +
                                                     std::to_string(c) + (e.ptr == nullptr ? " (windowed)" : " (raw)") +
                                                     (is_packed_edge(edge_key(box->morton_index, nm)) ? ", packed entry live" : ""));
                        }
                        rs = e.rows;
                        cs = 1;
                    }
                    add_gather(fills, max_fill_m, max_fill_n, ans_rows, nt, c, k, edge_src, rs, cs, identity, skel_local);
                    add_gather(fills, max_fill_m, max_fill_n, anr_rows, nt, c, r, edge_src, rs, cs, identity, red_local);
                    tag_fills(pass_of(it->second));
                    retired.push_back(e);
                } else {
                    assert_not_packed(false, edge_key(box->morton_index, nm), "the wave's near gather");
                    const IndexList rows = current_slots(ns);
                    add_eval(ans_rows, nt, c, k, rows, skel_slots);
                    add_eval(anr_rows, nt, c, r, rows, red_slots);
                }
            }

            if (nonsym_) {
                // X_RR = A_RR - A_RS T - T^T A_SR + T^T A_SS T,  X_SR = A_SR - A_SS T,
                // X_RS = A_RS - T^T A_SS,  X_RN = A_RN - T^T A_SN (no symmetric part)
                push(g2, wb.s, k, wb.T, wb.ldt, wb.tmp2, k, k, r, k);        // tmp2 = A_SS T       (NN)
                push(gn_rs_t, wb.xrs, r, wb.T, wb.ldt, wb.xrr, r, r, r, k);  // X_RR -= A_RS T      (NN)
                push(g3, wb.T, wb.ldt, wb.tmp2, k, wb.xrr, r, r, r, k);      // X_RR += T^T tmp2    (TN)
                push(gn_t_sr, wb.T, wb.ldt, wb.xsr, k, wb.xrr, r, r, r, k);  // X_RR -= T^T A_SR    (TN)
                push(g4, wb.s, k, wb.T, wb.ldt, wb.xsr, k, k, r, k);         // X_SR -= A_SS T      (NN)
                push(gn_t_rows, wb.T, wb.ldt, wb.s, k, wb.xrs, r, r, k, k);  // X_RS -= T^T A_SS    (TN)
                push(gn_t_rows, wb.T, wb.ldt, wb.asn, k, wb.xrn, r, r, nt, k);  // X_RN -= T^T A_SN
                add_gather(copy_full, max_copy_m, max_copy_n, wb.xrr_full, r, r, r, wb.xrr, 1, r, identity, identity);
                // the original X_RS (the host's X_RS_entry), and the copies of
                // X_RS and X_RN that P^T is applied to
                add_gather(copy_full, max_copy_m, max_copy_n, reinterpret_cast<S*>(d_result + wb.off_xrs), r, r, k,
                           wb.xrs, 1, r, identity, identity);
                add_gather(copy_full, max_copy_m, max_copy_n, wb.xrs_p, r, r, k, wb.xrs, 1, r, identity, identity);
                add_gather(copy_full, max_copy_m, max_copy_n, wb.xrn_p, r, r, nt, wb.xrn, 1, r, identity, identity);
            } else {
            // X_RR = A_RR - A_SR^T T - T^T A_SR + T^T A_SS T,  X_SR = A_SR - A_SS T
            push(g1, wb.xsr, k, wb.T, wb.ldt, wb.tmp1, r, r, r, k);      // tmp1 = -A_SR^T T   (TN)
            sym_adds.push_back(SymAddItem{wb.xrr, r, wb.tmp1, r, r});
            push(g2, wb.s, k, wb.T, wb.ldt, wb.tmp2, k, k, r, k);        // tmp2 = A_SS T       (NN)
            push(g3, wb.T, wb.ldt, wb.tmp2, k, wb.xrr, r, r, r, k);      // X_RR += T^T tmp2    (TN)
            push(g4, wb.s, k, wb.T, wb.ldt, wb.xsr, k, k, r, k);         // X_SR -= A_SS T      (NN)
            // the LU reads both triangles; the symmetric path assumes X_RR^T = X_RR
            // (compute_and_modify factors the symmetric part too)
            symmetrizes.push_back(IdentityItem{wb.xrr, r, r});
            add_gather(copy_full, max_copy_m, max_copy_n, wb.xrr_full, r, r, r, wb.xrr, 1, r, identity, identity);
            transposes.push_back(TransposeItem{wb.xsr, k, reinterpret_cast<S*>(d_result + wb.off_xrs), r, k, r});
            max_tr_m = std::max(max_tr_m, k);
            max_tr_n = std::max(max_tr_n, r);
            }
            if (gpu_sketch_) {
                add_gather(copy_full, max_copy_m, max_copy_n, reinterpret_cast<S*>(d_result + wb.off_t), k, k, r,
                           wb.T, 1, wb.ldt, identity, identity);
                add_gather(copy_full, max_copy_m, max_copy_n, wb.persist_xrr, r, r, r, wb.xrr, 1, r, identity, identity);
                if (compress_factors_ < 1 && nt > 0) {
                    // (with compression it runs with its sub-group, below)
                    transposes.push_back(TransposeItem{wb.temp2, nt, wb.persist_temp2t, r, nt, r});
                    max_tr_m = std::max(max_tr_m, nt);
                    max_tr_n = std::max(max_tr_n, r);
                }
            }
            push(getrf, wb.xrr, r, nullptr, 1, nullptr, 1, r, r, 1);
            // X_NR = A_NR - A_NS T
            push(g5, wb.ans, nt, wb.T, wb.ldt, wb.xnr, nt, nt, r, k);
            if (wb.has_remote && nt > 0 && !nonsym_) {
                char* source = d_xnr != nullptr ? d_xnr : d_sources[static_cast<size_t>(wb.source_group)];
                add_gather(copy_xnr, max_xnr_m, max_xnr_n, reinterpret_cast<S*>(source + wb.off_xnr_orig), nt, nt, r,
                           wb.xnr, 1, nt, identity, identity);
            }
            // temp1 = -X_SR X_RR^{-1},  temp2 = -X_NR X_RR^{-1}: with X_RR = P L U,
            // W = U^{-1} L^{-1} from two small solves (tmp1 is free after the
            // X_RR update), then one product each, then the row swaps of P
            identities.push_back(IdentityItem{wb.tmp1, r, r});
            push(trsm, wb.xrr, r, nullptr, 1, wb.tmp1, r, r, r, 1);
            push(gsolve, wb.xsr, k, wb.tmp1, r, wb.temp1, k, k, r, r);
            push(gsolve, wb.xnr, nt, wb.tmp1, r, wb.temp2, nt, nt, r, r);
            swaps.push_back(ColumnSwapItem{wb.temp1, k, k, wb.piv, r});
            if (nt > 0) swaps.push_back(ColumnSwapItem{wb.temp2, nt, nt, wb.piv, r});
            max_r = std::max(max_r, r);
            if (nonsym_) {
                // temp3 = -X_RR^{-1} X_RS = -W (P^T X_RS),  temp4 = -W (P^T X_RN)
                row_swaps_fwd.push_back(RowSwapItemT<S>{wb.xrs_p, r, k, wb.piv, r, 1});
                if (nt > 0) row_swaps_fwd.push_back(RowSwapItemT<S>{wb.xrn_p, r, nt, wb.piv, r, 1});
                push(gsolve, wb.tmp1, r, wb.xrs_p, r, wb.temp3, r, r, k, r);
                push(gsolve, wb.tmp1, r, wb.xrn_p, r, wb.temp4, r, r, nt, r);
                if (nt > 0) {
                    add_gather(copy_t4, max_t4_m, max_t4_n, wb.persist_temp4, r, r, nt, wb.temp4, 1, r, identity,
                               identity);
                }
                // S = A_SS + temp1 X_RS,  X_NS update = temp2 X_RS,  X_SN update = X_SR temp4
                push(g6n, wb.temp1, k, wb.xrs, r, wb.s, k, k, k, r);
                push(g7n, wb.temp2, nt, wb.xrs, r, wb.xns, nt, nt, k, r);
                push(g7n, wb.xsr, k, wb.temp4, r, wb.xsn, k, k, nt, r);
            } else {
            // S = A_SS + temp1 X_SR^T,  X_NS update = temp2 X_SR^T
            push(g6, wb.temp1, k, wb.xsr, k, wb.s, k, k, k, r);
            push(g7, wb.temp2, nt, wb.xsr, k, wb.xns, nt, nt, k, r);
            }

            st.schur = DeviceMatrix{wb.s, k, k};
            note_write(&st.schur, static_cast<size_t>(k) * k * D);
        }

        stats.pb_items += std::chrono::duration<double>(clock::now() - t_items).count();
        stats.plan_boxes += std::chrono::duration<double>(clock::now() - t0).count();
        // ---- metadata image: item arrays and MAGMA arrays
        const size_t off_evals = meta_.append(evals);
        // the window's passes run one after another, so the items must be
        // contiguous by pass: each pass loads the window, then gathers from it
        std::vector<size_t> pass_count;
        if (fill_window != nullptr) {
            while (fill_pass.size() < fills.size()) fill_pass.push_back(0);
            std::vector<size_t> order(fills.size());
            std::iota(order.begin(), order.end(), size_t{0});
            std::stable_sort(order.begin(), order.end(),
                             [&](size_t a, size_t b) { return fill_pass[a] < fill_pass[b]; });
            std::vector<GatherItem> sorted;
            sorted.reserve(fills.size());
            pass_count.assign(win_passes.size(), 0);
            for (size_t i : order) {
                sorted.push_back(fills[i]);
                if (fill_pass[i] < pass_count.size()) ++pass_count[fill_pass[i]];
            }
            fills.swap(sorted);
        }
        const size_t off_fills = meta_.append(fills);
        const size_t off_copy_full = meta_.append(copy_full);
        const size_t off_identities = meta_.append(identities);
        const size_t off_symmetrizes = meta_.append(symmetrizes);
        const size_t off_copy_xnr = meta_.append(copy_xnr);
        const size_t off_sym = meta_.append(sym_adds);
        const size_t off_swaps = meta_.append(swaps);
        const size_t off_transposes = meta_.append(transposes);
        const size_t off_row_swaps_fwd = meta_.append(row_swaps_fwd);
        const size_t off_copy_t4 = meta_.append(copy_t4);
        for (VBatch<S>* b : {&g1, &g2, &g3, &g4, &g5, &g6, &g7, &trsm, &getrf, &gsolve}) b->stage(meta_);
        for (VBatch<S>* b : {&gn_rs_t, &gn_t_sr, &gn_t_rows, &g6n, &g7n}) b->stage(meta_);
        // CA levels: one LU launch, and one launch of each solve for
        // X_RR^{-1}, per size class (getrf's and trsm's entries are in the chunk's
        // order, grouped by class above)
        std::vector<VBatch<S>> lu_classes, trsm_classes;
        std::vector<size_t> lu_first;  // first getrf entry of each class
        if (ca_) {
            if (getrf.entries.size() != chunk_active.size()) throw std::runtime_error("LevelEliminator: one LU per chunk_active box expected");
            auto split = [](const VBatch<S>& batch, std::vector<VBatch<S>>& classes, std::vector<size_t>* first) {
                for (size_t e = 0; e < batch.entries.size(); ++e) {
                    const int c = lu_size_class(batch.entries[e].n);
                    if (classes.empty() || lu_size_class(classes.back().entries.front().n) != c) {
                        classes.emplace_back();
                        if (first) first->push_back(e);
                    }
                    classes.back().entries.push_back(batch.entries[e]);
                }
            };
            split(getrf, lu_classes, &lu_first);
            split(trsm, trsm_classes, nullptr);
            for (VBatch<S>& b : lu_classes) b.stage(meta_);
            for (VBatch<S>& b : trsm_classes) b.stage(meta_);
        }
        // getrf pivot pointers: one per matrix, in getrf order
        std::vector<int*> piv_ptrs;
        for (int i : chunk_active) piv_ptrs.push_back(boxes[static_cast<size_t>(i)].piv);
        const size_t off_piv_ptrs = meta_.append(piv_ptrs);
        stats.plan += std::chrono::duration<double>(clock::now() - t0).count();
        ++stats.pb_chunks;
        if (gpu_idle_at_plan) {
            ++stats.pb_idle_chunks;
            stats.pb_idle_plan += std::chrono::duration<double>(clock::now() - t_chunk).count();
        }

        // ---- 3. device work, in stream order: the box region, then the owner
        // pass, planned on the host while the box region runs.  Nothing waits
        // for the device here: the next wave's sketch synchronizes.
        t0 = clock::now();
        char* md = meta_.upload(meta_device_, stream);
        stats.bytes_up += static_cast<double>(meta_.size() + t_bytes);
        if (t_bytes > 0) {
            check_cuda(cudaMemcpyAsync(d_t, t_image.data(), t_bytes, cudaMemcpyHostToDevice, stream), "upload T");
        }
        auto items = [&](size_t offset) { return md + offset; };
        // (the marks are made before the chunks and span them)

        marks.mark(stream);
        check_items(evals, "wave evals");
        check_items(fills, "wave fills");
        {
            EntryOriginScope origin(EntryOrigin::box);  // (collective entries' profile)
            eval_blocks(evals, meta_, md, off_evals, max_eval_m, max_eval_n, stream);
        }
        if (fill_window != nullptr) {
            auto& st_fw = eliminator_stats();
            st_fw.fw_passes += static_cast<int64_t>(win_passes.size());
            size_t at = 0;
            for (size_t pass = 0; pass < win_passes.size(); ++pass) {
                materialize_into_window(win_passes[pass], fill_window, stream);
                const size_t n = pass < pass_count.size() ? pass_count[pass] : 0;
                if (n > 0) {
                    launch_gather(reinterpret_cast<const GatherItem*>(items(off_fills)) + at,
                                  static_cast<int>(n), max_fill_m, max_fill_n, md, stream);
                }
                // the next pass overwrites the window, so this one's reads finish first
                check_cuda(cudaStreamSynchronize(stream), "gather window pass");
                at += n;
            }
            heap_.free(fill_window);
            fill_window = nullptr;
        } else {
            launch_gather(reinterpret_cast<const GatherItem*>(items(off_fills)), static_cast<int>(fills.size()),
                          max_fill_m, max_fill_n, md, stream);
        }
        // The fills were the last reads of the blocks this wave replaces: later
        // allocations may reuse them, every later launch being ordered after the
        // reads on the stream (not before: eval_blocks allocates and writes).
        for (const DeviceMatrix& m : retired) heap_.free(m.ptr);
        marks.mark(stream);  // 0: fills
        if (!g1.entries.empty()) g1.gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 0.0, queue);
        if (!sym_adds.empty()) {
            launch_sym_add(reinterpret_cast<const SymAddItem*>(items(off_sym)), static_cast<int>(sym_adds.size()),
                           max_r, stream);
        }
        g2.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
        if (!gn_rs_t.entries.empty()) gn_rs_t.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
        g3.gemm(md, MagmaTrans, MagmaNoTrans, 1.0, 1.0, queue);
        if (!gn_t_sr.entries.empty()) gn_t_sr.gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
        g4.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
        if (!gn_t_rows.entries.empty()) gn_t_rows.gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
        launch_symmetrize(reinterpret_cast<const IdentityItem*>(items(off_symmetrizes)), static_cast<int>(symmetrizes.size()),
                          max_r, stream);
        marks.mark(stream);  // 1: X_RR, X_SR
        launch_gather(reinterpret_cast<const GatherItem*>(items(off_copy_full)), static_cast<int>(copy_full.size()),
                      max_copy_m, max_copy_n, md, stream);
        if (ca_) {
            for (size_t c = 0; c < lu_classes.size(); ++c) {
                const VBatch<S>& b = lu_classes[c];
                getrf_vbatched<S>(lu_size_class(b.max_n), b.size_array(md, 0), b.size_array(md, 1),
                                  b.template pointer_array<S*>(md, 0), b.size_array(md, 3),
                                  reinterpret_cast<magma_int_t**>(md + off_piv_ptrs) + lu_first[c],
                                  d_info + c0 + lu_first[c],
                                  static_cast<magma_int_t>(b.count()), getrf_work_, queue);
            }
        } else {
            getrf_vbatched<S>(getrf.max_n, getrf.size_array(md, 0), getrf.size_array(md, 1),
                              getrf.template pointer_array<S*>(md, 0), getrf.size_array(md, 3),
                              reinterpret_cast<magma_int_t**>(md + off_piv_ptrs), d_info + c0,
                              static_cast<magma_int_t>(getrf.count()), getrf_work_, queue);
        }
        marks.mark(stream);  // 2: copies, LU
        g5.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 1.0, queue);
        launch_gather(reinterpret_cast<const GatherItem*>(items(off_copy_xnr)), static_cast<int>(copy_xnr.size()),
                      max_xnr_m, max_xnr_n, md, stream);
        marks.mark(stream);  // 3: X_NR
        launch_identity(reinterpret_cast<const IdentityItem*>(items(off_identities)), static_cast<int>(identities.size()),
                        max_r, stream);
        // W = U^{-1} L^{-1}; on CA levels per size class, with max_m = max_n =
        // the class (the lower solve's blocking follows the batch's largest box)
        auto solve_w = [&](const VBatch<S>& b, magma_int_t max_m, magma_int_t max_n) {
            trsm_vbatched<S>(MagmaRight, MagmaUpper, MagmaNoTrans, MagmaNonUnit, max_m, max_n,
                             b.size_array(md, 0), b.size_array(md, 1), 1.0,
                             b.template pointer_array<S*>(md, 0), b.size_array(md, 3),
                             b.template pointer_array<S*>(md, 2), b.size_array(md, 5),
                             static_cast<magma_int_t>(b.count()), queue);
            trsm_vbatched<S>(MagmaRight, MagmaLower, MagmaNoTrans, MagmaUnit, max_m, max_n,
                             b.size_array(md, 0), b.size_array(md, 1), 1.0,
                             b.template pointer_array<S*>(md, 0), b.size_array(md, 3),
                             b.template pointer_array<S*>(md, 2), b.size_array(md, 5),
                             static_cast<magma_int_t>(b.count()), queue);
        };
        if (ca_) {
            for (const VBatch<S>& b : trsm_classes) solve_w(b, lu_size_class(b.max_n), lu_size_class(b.max_n));
        } else {
            solve_w(trsm, trsm.max_m, trsm.max_n);
        }
        if (!row_swaps_fwd.empty()) {  // (P^T on the copies of X_RS and X_RN, before their solves)
            launch_row_swaps(reinterpret_cast<const RowSwapItemT<S>*>(items(off_row_swaps_fwd)),
                             static_cast<int>(row_swaps_fwd.size()), stream);
        }
        gsolve.gemm(md, MagmaNoTrans, MagmaNoTrans, -1.0, 0.0, queue);
        launch_column_swaps(reinterpret_cast<const ColumnSwapItem*>(items(off_swaps)), static_cast<int>(swaps.size()), stream);
        if (!copy_t4.empty()) {
            launch_gather(reinterpret_cast<const GatherItem*>(items(off_copy_t4)), static_cast<int>(copy_t4.size()),
                          max_t4_m, max_t4_n, md, stream);
        }
        marks.mark(stream);  // 4: temp1, temp2
        for (const auto& e : trsm.entries) stats.solve_flops += kFlopScale * 2.0 * e.m * e.n * e.n;  // two triangular solves
        for (const auto& e : gsolve.entries) stats.solve_flops += kFlopScale * 2.0 * e.m * e.n * e.k;
        {   // (the box region's work by the marks it is timed with, BPACK_TRACE=wave)
            auto vb = [](const VBatch<S>& b) {
                double f = 0.0;
                for (const auto& e : b.entries) f += 2.0 * e.m * e.n * e.k;
                return f;
            };
            stats.fl_box[1] += kFlopScale * (vb(g1) + vb(g2) + vb(g3) + vb(g4));
            for (const auto& e : getrf.entries) stats.fl_box[2] += kFlopScale * 2.0 / 3.0 * e.m * e.m * e.m;
            stats.fl_box[3] += kFlopScale * vb(g5);
            for (const auto& e : trsm.entries) stats.fl_box[4] += kFlopScale * 2.0 * e.m * e.n * e.n;
            stats.fl_box[4] += kFlopScale * vb(gsolve);
            stats.fl_box[5] += kFlopScale * (vb(g6) + vb(g7));
        }
        launch_transpose(reinterpret_cast<const TransposeItem*>(items(off_transposes)), static_cast<int>(transposes.size()),
                         max_tr_m, max_tr_n, stream);
        if (!g6.entries.empty()) g6.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 1.0, queue);
        if (!g7.entries.empty()) g7.gemm(md, MagmaNoTrans, MagmaTrans, 1.0, 0.0, queue);
        if (!g6n.entries.empty()) g6n.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 1.0, queue);
        if (!g7n.entries.empty()) g7n.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
        marks.mark(stream);  // 5: transposes, Schur
        stats.launch += std::chrono::duration<double>(clock::now() - t0).count();

        // ---- the wave's fill sources, compressed a slot at a time (the readers
        // take one slot's rows, so a slot is the block that decompresses alone).
        // The transposes above wrote temp2^T and run on this stream, so the
        // compression is ordered after them.  Nothing reads the compressed form
        // yet: this measures the ratio and, with BPACK_CHECK=fillsrc, checks the
        // round trip against the original.
        if (compress_factors_ >= 1 && gpu_sketch_) compress_fill_sources(boxes, chunk_active, stream);

        // ---- The generators, this wave's kept factors, and then the raw outputs.
        //
        // These three used to run after the owner pass, which was fine while the
        // owner pass read temp2 and X_NR raw.  Now that it reads them out of their
        // stores, they are the last readers of the raw form, and moving them here
        // is what lets the 5.14 GB of a leaf wave's 6.29 GB transient go before
        // the pass that used to keep it alive.  None of the three depends on the
        // owner pass: they read the solves' output.
        if (device_exchange_) {
            if (owner_ != nullptr) {
                copy_owner_generators(boxes, chunk_active, d_sources, d_temp2, d_xnr, d_result, stream);
            } else {
                copy_generators(boxes, chunk_active, d_sources, d_temp2, d_xnr, d_result, stream);
            }
            if (d_xnr != nullptr) {
                heap_.free(d_xnr);  // read by the chunk's packs only, and they have run
                d_xnr = nullptr;
            }
        }
        for (int i : chunk_active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            if (wb.k == 0 || wb.r == 0) continue;
            if (ghost_factors_elsewhere(wb.box->morton_index)) continue;
            const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
            factor_bytes_ += static_cast<double>(2 * align_up(k * r * sizeof(S)) + align_up(nt * r * sizeof(S)) +
                                                 align_up(r * r * sizeof(S)) + align_up(r * sizeof(int)));
        }
        if (keep_solve_ && !keep_failed_) keep_solve_factors(boxes, chunk_active, d_result, d_temp2, stream);
        // Whether the host still wants a group's X_NR: with the factors kept on
        // the device it does not (rescue_skipped_xnr reads them back if it ever
        // needs them), and then temp2 has no reader left outside its store.
        bool chunk_stored = compress_factors_ >= 1;
        for (int i : chunk_active) {
            if (!chunk_stored) break;
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            if (wb.r == 0 || wb.ntot == 0) continue;
            const BoxState* sst = wb.state;
            if (sst == nullptr || packed_of(*sst) == nullptr || sst->xnrt_packed == nullptr ||
                sst->t2t_slots.size() != sst->xnrt_slots.size()) {
                chunk_stored = false;
            }
        }
        // the owner pass reads the whole wave, so it may only take the stores
        // if every chunk stored its own
        from_store = from_store && chunk_stored;
        std::vector<char> temp2_needed(d_temp2.size(), 0);
        for (int i : chunk_active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            if (wb.ntot == 0) continue;
            if (!host_skips_xnr(wb) && !host_drops_ghost_xnr(wb)) temp2_needed[static_cast<size_t>(wb.source_group)] = 1;
        }
        // This chunk's raw output, now that its stores are built and the last
        // readers of the raw form - the generators and the keep - have run.
        if (chunk_stored) {
            for (size_t b = xnr_first; b < xnr_blocks.size(); ++b) heap_.free(xnr_blocks[b]);
            xnr_blocks.resize(xnr_first);
            for (size_t g = 0; g < d_temp2.size(); ++g) {
                if (group_chunk_of[g] != chunk || temp2_needed[g] || d_temp2[g] == nullptr) continue;
                heap_.free(d_temp2[g]);
                d_temp2[g] = nullptr;
            }
            for (int i : chunk_active) {  // nothing may read the raw form again
                WaveBox& wb = boxes[static_cast<size_t>(i)];
                wb.xnr = nullptr;
                wb.temp2 = nullptr;
            }
        }
        // ---- Step 5, planned while the device runs the box region: the near
        // block of every neighbor becomes A_NS + update (current neighbor rows x
        // this box's skeleton), stored in the canonical orientation.  Its tens of
        // thousands of block allocations per wave kept the device idle when they
        // came before the launches above; they may now take the blocks those
        // launches retired.  The items go with the owner pass's image.
        {
            const auto t_near = clock::now();
            for (int i : chunk_active) {
                WaveBox& wb = boxes[static_cast<size_t>(i)];
                Box* box = wb.box;
                const int k = wb.k, nt = wb.ntot;
                for (size_t a = 0; a < box->one_hop.size(); ++a) {
                    const int64_t nm = box->one_hop[a];
                    const int c = static_cast<int>(wb.counts[a]);
                    if (c == 0) continue;
                    if (nonsym_) {
                        // K(N, B) = A_NS + temp2 X_RS (c x k), K(B, N) = A_SN + X_SR temp4 (k x c)
                        DeviceMatrix ns_block{heap_.alloc_resident<S>(static_cast<size_t>(c) * k * D, ResidentKind::LevelSchur),
                                              c, k};
                        DeviceMatrix sn_block{heap_.alloc_resident<S>(static_cast<size_t>(c) * k * D, ResidentKind::LevelSchur),
                                              k, c};
                        stores.push_back(AddStoreItem{ns_block.ptr, 1, c, c, k, wb.ans + wb.row0[a], nt,
                                                      wb.xns + wb.row0[a], nt});
                        const size_t col0 = static_cast<size_t>(wb.row0[a]) * k;
                        stores.push_back(AddStoreItem{sn_block.ptr, 1, k, k, c, wb.asn + col0, k, wb.xsn + col0, k});
                        max_store_m = std::max({max_store_m, c, k});
                        max_store_n = std::max({max_store_n, c, k});
                        for (const auto& [key, block] : {std::make_pair(dkey(nm, box->morton_index), ns_block),
                                                         std::make_pair(dkey(box->morton_index, nm), sn_block)}) {
                            edges_[key] = block;
                            note_write(&edges_[key], static_cast<size_t>(block.rows) * block.cols * D);
                        }
                        continue;
                    }
                    DeviceMatrix fresh;
                    const auto t_alloc = clock::now();
                    fresh.ptr = heap_.alloc_resident<S>(static_cast<size_t>(c) * k * D, ResidentKind::LevelSchur);
                    stats.pb_near_alloc += std::chrono::duration<double>(clock::now() - t_alloc).count();
                    int64_t ors = 1, ocs = c;
                    if (box->morton_index < nm) {
                        fresh.rows = c;  // rows = neighbor (hi)
                        fresh.cols = k;
                    } else {
                        fresh.rows = k;  // rows = this box (hi)
                        fresh.cols = c;
                        ors = k;
                        ocs = 1;
                    }
                    stores.push_back(AddStoreItem{fresh.ptr, ors, ocs, c, k, wb.ans + wb.row0[a], nt, wb.xns + wb.row0[a], nt});
                    max_store_m = std::max(max_store_m, c);
                    max_store_n = std::max(max_store_n, k);
                    const uint64_t fkey = edge_key(box->morton_index, nm);
                    edges_[fkey] = fresh;
                    note_write(&edges_[fkey], static_cast<size_t>(fresh.rows) * fresh.cols * D);
                    fresh_near.push_back(FreshNear{&edges_[fkey], fkey, fresh.rows, fresh.cols});
                }
            }
            const double seconds = std::chrono::duration<double>(clock::now() - t_near).count();
            stats.pb_near += seconds;
            stats.plan += seconds;
        }

        // The near pass's stores first, in an image of their own: they write the
        // fresh blocks of step 5 and carry no metadata offsets.
        {
            owner_meta_.clear();
            const size_t off_stores = owner_meta_.append(stores);
            char* smd = owner_meta_.upload(owner_meta_device_, stream);
            stats.bytes_up += static_cast<double>(owner_meta_.size());
            check_items(stores, "near stores");
            launch_add_store(reinterpret_cast<const AddStoreItem*>(smd + off_stores), static_cast<int>(stores.size()),
                             max_store_m, max_store_n, stream);
            check_cuda(cudaStreamSynchronize(stream), "owner stores");  // the image is reused below
        }
#ifdef H2_HAVE_NVCOMP
        // Step 5's blocks into the near store now, in groups bounded as the
        // owner pass's are.  They used to sit raw until the first owner group
        // to update them came by, which at the 192^3 leaf was 0.88 GB of the
        // resident set at the peak.  Nothing in this wave writes them again:
        // the owner pass writes pairs of a source's neighbours, and these are
        // the source's own.  A later wave's owner group puts one back exactly
        // as it puts back any other packed target.
        const bool np5_pressure = under_pressure();
        if (kHoldNearBlocksPacked && compress_factors_ >= 1 && !fresh_near.empty() && (np5_pressure || pack_shared_)) {
            const auto t_np5 = std::chrono::steady_clock::now();
            const double raw0 = stats.np_raw, stored0 = stats.np_stored;
            const int64_t blocks0 = stats.np_blocks;
            size_t at = 0;
            while (at < fresh_near.size()) {
                size_t bytes = 0, end = at;
                while (end < fresh_near.size() && (end == at || bytes < kNearGroupBudget)) {
                    bytes += static_cast<size_t>(fresh_near[end].rows) * fresh_near[end].cols * D;
                    ++end;
                }
                if (replicated_) pack_near_by_rule(fresh_near, at, end, stream, np5_pressure);
                else pack_near_group(fresh_near, at, end, stream);
                ++stats.np5_groups;
                at = end;
            }
            stats.np5_raw += stats.np_raw - raw0;
            stats.np5_stored += stats.np_stored - stored0;
            stats.np5_blocks += stats.np_blocks - blocks0;
            stats.np5_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_np5).count();
        }
#endif
        marks.mark(stream);  // 6: step 5 (near blocks)
        // A_NS and X_NS were the last of this chunk's raw output with a
        // reader, and the stores above have been launched and waited for.
        for (size_t b = work_first; b < work_blocks.size(); ++b) heap_.free(work_blocks[b]);
        work_blocks.resize(work_first);
        // (an owner level: the transport between the chunks too.  A payload
        // lived about three waves when only the wave boundaries moved it -
        // issued at the sender's next boundary, received at the receiver's,
        // installed at the one after - 14.8 GB a rank, CA3 at 768^3 on 64
        // ranks.  The near store is consistent here (the chunk's packs have
        // run), and an install touches no block of this colour's waves here.)
        if (owner_ != nullptr) {
            owner_progress();
            owner_install_arrived(wave_index / std::max(owner_->num_waves, 1));
        }
    }
    // (after every chunk: a last chunk with no boxes is skipped above)
    if (device_exchange_ && owner_ != nullptr) finish_owner_generators(stream);
    if (d_xnr != nullptr) {
        heap_.free(d_xnr);
        d_xnr = nullptr;
    }
    elim_marks_.emplace_back(std::make_unique<StreamMarks>(), 7);  // the owner GEMMs
    StreamMarks& marks = *elim_marks_.back().first;
    marks.mark(stream);  // where the production ends


    const auto t_owner = clock::now();
    // The slot each neighbor position holds in the two stores: both skip the
    // neighbors with no rows, so a position is not a slot index.
    for (WaveBox& wb : boxes) {
        wb.slot_of_pos.assign(wb.counts.size(), -1);
        int32_t next = 0;
        for (size_t a = 0; a < wb.counts.size(); ++a) {
            if (wb.counts[a] > 0) wb.slot_of_pos[a] = next++;
        }
    }
    auto win_key = [](int32_t src, int32_t slot, int which) {
        return (static_cast<uint64_t>(src) << 33) | (static_cast<uint64_t>(slot) << 1) |
               static_cast<uint64_t>(which);
    };
    // the slot's uncompressed bytes, which is where it lands in the window
    auto slot_bytes = [&](int32_t src, int32_t slot, int which) -> size_t {
        const WaveBox& sb = boxes[static_cast<size_t>(src)];
        const auto& slots = which == 0 ? sb.state->t2t_slots : sb.state->xnrt_slots;
        if (static_cast<size_t>(slot) >= slots.size()) return 0;
        return align_up(static_cast<size_t>(slots[static_cast<size_t>(slot)].rows) *
                        static_cast<size_t>(sb.r) * D);
    };
    owner_meta_.clear();
    const uint64_t owner_stamp = ~wave_stamp_;  // skeleton lists of this image
    // ---- deferred X_NN owner pass: for every source E of the wave and every
    // pair of its neighbors, target(lo, hi) += temp2_E[hi rows] X_NR_E[lo rows]^T,
    // and Schur(C) += temp2_E[C rows] X_NR_E[C rows]^T.  Contributions to a
    // target arrive in the CPU order: sources in the owner's one_hop order.
    // CA levels: the candidates are the writable boxes (local and ghosts); a
    // pair of two eliminated ghosts is never read here and is skipped.  A
    // pair's contributions arrive by source Morton there: its copies on other
    // ranks may be reached from the other endpoint.
    std::vector<int64_t> candidates;
    {
        std::unordered_set<int64_t> seen;
        for (int i : active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            if (wb.ntot == 0) continue;
            for (int64_t c : wb.box->one_hop) {
                if (is_interest(c) && wave_index_of.count(c) == 0 && seen.insert(c).second) candidates.push_back(c);
            }
        }
        std::sort(candidates.begin(), candidates.end());
    }
    auto t_phase = clock::now();
    auto owner_lap = [&](double& phase) {
        const auto now = clock::now();
        phase += std::chrono::duration<double>(now - t_phase).count();
        t_phase = now;
    };
    stats.po_candidates += std::chrono::duration<double>(t_phase - t_owner).count();
    std::vector<EvalItem> owner_evals;
    int max_owner_m = 0, max_owner_n = 0;
    // Every target belongs to one candidate (its Schur block; a pair with a
    // higher box, or with another rank's box): the candidates' targets and
    // tasks are found in parallel, each in the order of the sequential loop,
    // then created (blocks, their kernel values, the edge map) and joined in
    // candidate order, which is the sequential loop's target order.
    // Which slot of which source box a task's two operands are, so its group
    // can point them at the window it decompressed into.
    struct SrcSlot { int32_t src = -1; int32_t sa = -1; int32_t sb = -1; };
    struct CandidateTarget {
        uint64_t key = 0;
        DeviceMatrix* block = nullptr;  // the Schur block or an existing near block; else created at the join
        BoxState* row_box = nullptr;
        BoxState* col_box = nullptr;
        int rows = 0, cols = 0;
        std::vector<typename VBatch<S>::Entry> tasks;  // C set at the join
        std::vector<int64_t> sources;                 // of the tasks
        std::vector<SrcSlot> src_slots;               // ... and their slots
    };
    std::vector<std::vector<CandidateTarget>> candidate_targets(candidates.size());
    // (per candidate, summed after the loop: the dead updates not made, by kind
    // - two eliminated ghosts, an eliminated ghost's Schur block, with a box
    // beyond the halo - flops and tasks)
    std::vector<std::array<double, 6>> dead_by_candidate(candidates.size(), std::array<double, 6>{});
    std::exception_ptr owner_failure;
    std::mutex owner_failure_mutex;
    #pragma omp parallel for schedule(dynamic, 4)
    for (int64_t ci = 0; ci < static_cast<int64_t>(candidates.size()); ++ci) {
        try {
            const int64_t cm = candidates[static_cast<size_t>(ci)];
            auto& mine = candidate_targets[static_cast<size_t>(ci)];
            Box* cbox = held_box(cm);
            BoxState& cs = state_of(cm);
            // (a pair of two eliminated ghosts is never read; an owner level:
            // of two eliminated boxes of other ranks)
            const bool eliminated_ghost = cs.eliminated && (owner_ != nullptr ? !is_local(cm) : cs.ghost);
            auto target_for = [&](uint64_t key, DeviceMatrix* block, int rows, int cols, BoxState& row_box,
                                  BoxState& col_box) -> CandidateTarget& {
                for (CandidateTarget& t : mine) {
                    if (t.key == key) return t;
                }
                if (block != nullptr && block->ptr != nullptr && (block->rows != rows || block->cols != cols)) {
                    std::ostringstream oss;
                    oss << "LevelEliminator: owner target is " << block->rows << " x " << block->cols
                        << ", expected " << rows << " x " << cols;
                    throw std::runtime_error(oss.str());
                }
                mine.push_back(CandidateTarget{key, block, &row_box, &col_box, rows, cols, {}, {}});
                return mine.back();
            };
            auto near_block = [&](uint64_t key) -> DeviceMatrix* {
                auto it = edges_.find(key);
                return it != edges_.end() ? &it->second : nullptr;
            };
            for (int64_t em : cbox->one_hop) {
                auto wit = wave_index_of.find(em);
                if (wit == wave_index_of.end()) continue;
                const int32_t src_idx = static_cast<int32_t>(wit->second);
                const WaveBox& src = boxes[wit->second];
                if (src.r == 0 || src.ntot == 0) continue;
                // the raw outputs are gone when the operands come from the
                // stores, so the entries are built with null bases and the
                // group loop points them at its window
                auto op = [](const S* base, int row) -> const S* {
                    return base == nullptr ? nullptr : base + row;
                };
                const auto& src_hop = src.box->one_hop;
                const size_t c_pos = static_cast<size_t>(std::find(src_hop.begin(), src_hop.end(), cm) - src_hop.begin());
                if (c_pos == src_hop.size()) throw std::runtime_error("LevelEliminator: asymmetric one_hop lists");
                const int n_c = static_cast<int>(src.counts[c_pos]);
                if (n_c == 0) continue;
                const S* x_nr_c = op(src.xnr, src.row0[c_pos]);
                for (size_t a = 0; a < src_hop.size(); ++a) {
                    const int64_t am = src_hop[a];
                    const int n_a = static_cast<int>(src.counts[a]);
                    if (n_a == 0) continue;
                    auto& dead = dead_by_candidate[static_cast<size_t>(ci)];
                    const double task_flops = 2.0 * n_c * n_a * src.r;
                    if (eliminated_ghost && am != cm) {
                        const BoxState& as = state_of(am);
                        if (as.eliminated && (owner_ != nullptr ? !is_local(am) : as.ghost)) {
                            if (cm < am || !is_interest(am)) { dead[0] += task_flops; dead[1] += 1; }  // (once a pair)
                            continue;
                        }
                    }
                    // A replicated level: nor an eliminated ghost's own Schur block,
                    // nor its pairs with boxes beyond this rank's halo.  Its
                    // parent is another rank's, and the transition builds only
                    // this rank's parents (pairs with a box of this rank); the
                    // box beyond the halo is never eliminated here.  The owner
                    // makes both for real.
                    if (eliminated_ghost && owner_ == nullptr && ca_ && (am == cm || !is_interest(am))) {
                        const int k = am == cm ? 2 : 4;
                        dead[k] += task_flops;
                        dead[k + 1] += 1;
                        continue;
                    }
                    CandidateTarget* target = nullptr;
                    // (an owner level: a block is updated here only if it is read
                    // here, see apply_remote_generators; an owned ghost's blocks
                    // are spent after its elimination here)
                    if (owner_ != nullptr && !is_local(cm) && cs.eliminated) {
                        if (am == cm || !is_local(am)) {
                            if (am == cm || !(is_writable(am) && !state_of(am).eliminated)) continue;
                        }
                    }
                    if (nonsym_) {
                        // every ordered pair: K(x, y) += temp2_E[x rows] X_RN_E[:, y cols] (NN)
                        auto task = [&](CandidateTarget& t, size_t x_pos, int nx, size_t y_pos, int ny) {
                            t.tasks.push_back({src.temp2 + src.row0[x_pos],
                                               src.xrn + static_cast<size_t>(src.row0[y_pos]) * src.r, nullptr,
                                               nx, ny, src.r, src.ntot, src.r, nx});
                            t.sources.push_back(em);
                            t.src_slots.push_back(SrcSlot{src_idx, -1, -1});
                        };
                        if (am == cm) {
                            if (!is_writable(cm) && !(is_local(cm) && cs.eliminated)) continue;
                            task(target_for(~static_cast<uint64_t>(cm), &cs.schur, n_c, n_c, cs, cs), c_pos, n_c,
                                 c_pos, n_c);
                        } else if (deferred_xnn_boxes_are_one_hop(level_.dimension, cm, am) &&
                                   (cm < am || !is_interest(am))) {
                            // both of the pair's blocks, with the lower candidate (or this
                            // rank's copies of a pair with another rank's box)
                            BoxState& as = state_of(am);
                            task(target_for(dkey(am, cm), near_block(dkey(am, cm)), n_a, n_c, as, cs), a, n_a, c_pos,
                                 n_c);
                            task(target_for(dkey(cm, am), near_block(dkey(cm, am)), n_c, n_a, cs, as), c_pos, n_c, a,
                                 n_a);
                        }
                        continue;
                    }
                    if (am == cm) {
                        // (a box eliminated elsewhere: its Schur block is read once its
                        // owner's has landed here, i.e. once it is eliminated)
                        if (!is_writable(cm) && !(is_local(cm) && cs.eliminated)) continue;
                        target = &target_for(~static_cast<uint64_t>(cm), &cs.schur, n_c, n_c, cs, cs);
                    } else if (cm < am && deferred_xnn_boxes_are_one_hop(level_.dimension, cm, am)) {
                        BoxState& as = state_of(am);
                        const uint64_t key = edge_key(cm, am);
                        target = &target_for(key, near_block(key), n_a, n_c, as, cs);
                    } else if (!is_interest(am) && deferred_xnn_boxes_are_one_hop(level_.dimension, cm, am)) {
                        // A pair with a lower box of another rank: this rank's
                        // copy (the other rank's comes from the generator), rows =
                        // the candidate (hi): temp2_E[c rows] X_NR_E[a rows]^T, as
                        // every pair (the transpose of the lo-rows product would
                        // round differently: its terms differ)
                        BoxState& as = state_of(am);
                        const uint64_t key = edge_key(cm, am);
                        target = &target_for(key, near_block(key), n_c, n_a, cs, as);
                        target->tasks.push_back({op(src.temp2, src.row0[c_pos]), op(src.xnr, src.row0[a]), nullptr,
                                                 n_c, n_a, src.r, src.ntot, src.ntot, n_c});
                        target->sources.push_back(em);
                        target->src_slots.push_back(SrcSlot{src_idx, slot_at(src, c_pos), slot_at(src, a)});
                        continue;
                    } else {
                        continue;  // owned by the other endpoint, or a lazy far pair
                    }
                    target->tasks.push_back({op(src.temp2, src.row0[a]), x_nr_c, nullptr,
                                             n_a, n_c, src.r, src.ntot, src.ntot, n_a});
                    target->sources.push_back(em);
                    target->src_slots.push_back(SrcSlot{src_idx, slot_at(src, a), slot_at(src, c_pos)});
                }
            }
            if (ca_) {  // (a pair's contributions by source Morton)
                for (CandidateTarget& t : mine) {
                    if (t.tasks.size() < 2) continue;
                    std::vector<size_t> order(t.tasks.size());
                    std::iota(order.begin(), order.end(), size_t{0});
                    std::stable_sort(order.begin(), order.end(),
                                     [&](size_t x, size_t y) { return t.sources[x] < t.sources[y]; });
                    std::vector<typename VBatch<S>::Entry> sorted;
                    sorted.reserve(order.size());
                    for (size_t x : order) sorted.push_back(t.tasks[x]);
                    t.tasks.swap(sorted);
                    // The sources and their slots travel with their tasks.  A
                    // task's operands are named by its slot pair, so
                    // permuting the tasks alone hands the GEMM another box's
                    // rows - which is what it did: NaN in the next level's ID
                    // input, on the CA path only, because only it sorts.
                    if (t.src_slots.size() == t.tasks.size()) {
                        std::vector<SrcSlot> ss;
                        ss.reserve(order.size());
                        for (size_t x : order) ss.push_back(t.src_slots[x]);
                        t.src_slots.swap(ss);
                    }
                    std::vector<int64_t> src;
                    src.reserve(order.size());
                    for (size_t x : order) src.push_back(t.sources[x]);
                    t.sources.swap(src);
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(owner_failure_mutex);
            if (!owner_failure) owner_failure = std::current_exception();
        }
    }
    if (owner_failure) std::rethrow_exception(owner_failure);
    for (const auto& d : dead_by_candidate) {
        stats.dead_pair_flops += kFlopScale * d[0];
        stats.dead_pair_tasks += static_cast<int64_t>(d[1]);
        stats.dead_schur_flops += kFlopScale * d[2];
        stats.dead_schur_tasks += static_cast<int64_t>(d[3]);
        stats.dead_beyond_flops += kFlopScale * d[4];
        stats.dead_beyond_tasks += static_cast<int64_t>(d[5]);
    }
    owner_lap(stats.po_pairs);
    // the join: blocks of new targets (their kernel values first), in the
    // sequential loop's order
    // The join records a target; its block is allocated or materialized, its
    // kernel values evaluated and its C pointers filled in when its group
    // runs, so that only one group's blocks need to be raw at a time.
    struct Target {
        DeviceMatrix* block = nullptr;
        std::vector<typename VBatch<S>::Entry> tasks;
        uint64_t key = 0;
        int rows = 0, cols = 0;
        BoxState* row_box = nullptr;
        BoxState* col_box = nullptr;
        bool fresh = false;   // allocated by its group: its kernel values first
        std::vector<SrcSlot> src_slots;
    };
    std::vector<Target> targets;
    // (piece 4 put every adopted block back here, before the join, because the
    // join allocated them.  Step 3b-i moved allocation into the group loop
    // and the join only records shapes, so that pre-pass became redundant -
    // and harmful: it decompressed the whole set up front, 6.16 GB live at
    // the 192^3 leaf, which is exactly the bounding the groups exist to do.)
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        for (CandidateTarget& t : candidate_targets[ci]) {
            DeviceMatrix* block = t.block != nullptr ? t.block : &edges_[t.key];
            if (block->ptr == nullptr) {
                // allocated, or put back, when its group runs
            } else if (block->rows != t.rows || block->cols != t.cols) {
                std::ostringstream oss;
                oss << "LevelEliminator: owner target is " << block->rows << " x " << block->cols
                    << ", expected " << t.rows << " x " << t.cols;
                throw std::runtime_error(oss.str());
            }
            targets.push_back(Target{block, std::move(t.tasks), t.key, t.rows, t.cols, t.row_box, t.col_box,
                                     false, std::move(t.src_slots)});
        }
    }
    owner_lap(stats.po_order);
    // Sub-batch j holds the j-th contribution of every target, so no launch
    // writes a block twice and each block sees its contributions in order.
    //
    // In groups of a bounded set of blocks, so that holding the near-field
    // blocks compressed needs one group decompressed at a time rather than
    // the 10.61 GB a wave touches.  Each target lies wholly within one group
    // and keeps its contributions' order inside that group's sub-batches, so
    // every block still sees the same contributions in the same order and the
    // bits do not move; only which GEMMs share a launch changes, which the
    // kernels' batch-independence covers (tests/batch_determinism.cpp).  The
    // group boundaries follow the targets' canonical order and their shapes,
    // both rank-independent, so replicas group alike.
    // Holding the near blocks packed is **off**: measured at 192^3 it raises
    // the leaf peak (20.72 -> 21.01 GB Color, 22.51 -> 23.48 CA0) and costs
    // 0.76 / 1.70 s, because the per-wave materialization is still unbounded
    // (6.16 GB of decompressed copies at the peak) while the store itself
    // works (near blocks 9.47 -> 3.24 GB).  GPU_CA_PLAN.md step 3c lists what
    // remains: bounded windows for the three readers, and a group budget
    // chosen for the GEMM launch count rather than guessed.
    //
    // With it off the targets form a single group, which is the launch
    // structure the owner pass had before any of this: the budget alone takes
    // the leaf from 120 batched GEMM calls to 730 and costs 0.64 s, whether
    // anything is compressed or not.

    std::vector<std::pair<size_t, size_t>> owner_groups;  // [first, last) of targets
    {
        // Two budgets: the targets the group writes, and - when the operands
        // come from their stores - the source slots it has to decompress, so
        // the window is bounded by construction rather than by what the data
        // happens to need.
        std::unordered_set<uint64_t> seen;
        const size_t group_budget = near_group_budget(), window_budget = owner_window_bytes();
        size_t at = 0;
        while (at < targets.size()) {
            size_t bytes = 0, win = 0, end = at;
            seen.clear();
            while (end < targets.size()) {
                if (end > at && (bytes >= group_budget || win >= window_budget)) break;
                bytes += static_cast<size_t>(targets[end].rows) * targets[end].cols * D;
                if (from_store) {
                    for (const SrcSlot& ss : targets[end].src_slots) {
                        for (int which = 0; which < 2; ++which) {
                            const int32_t slot = which == 0 ? ss.sa : ss.sb;
                            if (ss.src < 0 || slot < 0) continue;
                            if (!seen.insert(win_key(ss.src, slot, which)).second) continue;
                            win += slot_bytes(ss.src, slot, which);
                        }
                    }
                }
                ++end;
            }
            owner_groups.push_back({at, end});
            stats.owner_groups += 1;
            at = end;
        }
    }

    owner_lap(stats.po_batches);
    stats.plan_owner += std::chrono::duration<double>(clock::now() - t_owner).count();
    stats.plan += std::chrono::duration<double>(clock::now() - t_owner).count();
    t0 = clock::now();
    // Then one group at a time: its blocks raw, its kernel values, its GEMMs.
    // Each group stages its own image, because a slot list is an offset into
    // the image it was appended to (current_slots_in) and so cannot be shared
    // across uploads; the stamp is bumped per group to force the re-append.
    std::optional<GatherPlan> next_plan;   // the next group's gather, run under this group's GEMMs
    for (size_t gi = 0; gi < owner_groups.size(); ++gi) {
        const auto [first, last] = owner_groups[gi];
        auto tg = clock::now();   // the group's host time, by part
        auto lap = [&](double& acc) {
            const auto now = clock::now();
            acc += std::chrono::duration<double>(now - tg).count();
            tg = now;
        };
        if (env::check("fillsrc")) {
            // What this group distinctly needs of temp2 and X_NR.  A slot is
            // named by its address, which is unique per (source, neighbor), so
            // the map is the union a window would have to hold; summed over
            // groups it is what decompressing per group would cost.  Compare
            // with the wave's outputs, reported as the largest wave's
            // transient.
            std::unordered_map<const S*, size_t> slots;
            double with_repeats = 0.0;
            for (size_t i = first; i < last; ++i) {
                for (const auto& e : targets[i].tasks) {
                    const size_t a_bytes = static_cast<size_t>(e.m) * static_cast<size_t>(e.k) * D;
                    const size_t b_bytes = static_cast<size_t>(e.n) * static_cast<size_t>(e.k) * D;
                    slots.emplace(e.a, a_bytes);
                    slots.emplace(e.b, b_bytes);
                    with_repeats += static_cast<double>(a_bytes + b_bytes);
                }
            }
            size_t bytes = 0;
            for (const auto& kv : slots) bytes += kv.second;
            auto& st = eliminator_stats();
            st.os_union_max = std::max(st.os_union_max, static_cast<double>(bytes));
            st.os_union_sum += static_cast<double>(bytes);
            st.os_raw_max = std::max(st.os_raw_max, with_repeats);
            st.os_slots += static_cast<int64_t>(slots.size());
            ++st.os_groups;
        }
        // The blocks this group's pack will take back - by the rule the pack
        // applies after the GEMMs, with the pressure foreseen from what the
        // group adds to the heap - in one heap block (group_arena_).  A block
        // the pack leaves raw after all moves to a block of its own before
        // that one goes back, so the decision itself stays where it was.
        std::unordered_set<const DeviceMatrix*> arena_set;
        if (kHoldNearBlocksPacked && compress_factors_ >= 1) {
            size_t group_raw = 0;
            for (size_t i = first; i < last; ++i) {
                const Target& t = targets[i];
                if (t.block->ptr == nullptr) group_raw += arena_round(static_cast<size_t>(t.rows) * t.cols * D);
            }
            const bool pressure = heap_.used() + group_raw > hold_limit();
            size_t arena_bytes = 0;
            for (size_t i = first; i < last; ++i) {
                const Target& t = targets[i];
                if (t.block->ptr != nullptr || t.rows <= 0 || t.cols <= 0) continue;
                bool packs = pressure;
                if (replicated_ && !replicas_done_) {
                    packs = block_shared(static_cast<uint64_t>(t.key)) ? pack_shared_ : pressure;
                }
                if (!packs) continue;
                arena_set.insert(t.block);
                arena_bytes += arena_round(static_cast<size_t>(t.rows) * t.cols * D);
            }
            if (arena_bytes > 0) {
                group_arena_.base = heap_.try_alloc_resident(arena_bytes, ResidentKind::Adopted);
                group_arena_.bytes = group_arena_.base != nullptr ? arena_bytes : 0;
                group_arena_.used = 0;
                if (group_arena_.base != nullptr) {
                    ++stats.og_arena_groups;
                    stats.og_arena_blocks += static_cast<int64_t>(arena_set.size());
                } else {
                    ++stats.og_arena_declined;   // (no range for it: a heap block each, as before)
                    arena_set.clear();
                }
            }
        }
        const RawTake arena_take = [&](const AdoptedBlock& w, size_t bytes) -> char* {
            return arena_set.count(w.m) != 0 ? group_arena_take(bytes) : nullptr;
        };
        {   // every block of the group raw: the ones held packed put back in one
            // batch, the ones that do not exist yet allocated and marked for
            // their kernel values
            std::vector<AdoptedBlock> want;
            for (size_t i = first; i < last; ++i) {
                Target& t = targets[i];
                if (t.block->ptr != nullptr) continue;
                const bool schur = static_cast<int64_t>(t.key) < 0;
                const int64_t morton = schur ? static_cast<int64_t>(~t.key) : 0;
                if (schur ? is_packed_schur(morton) : is_packed_edge(t.key)) {
                    t.block->rows = t.rows;
                    t.block->cols = t.cols;
                    want.push_back(AdoptedBlock{schur, schur ? static_cast<uint64_t>(morton) : t.key, t.block});
                    continue;
                }
                const size_t raw = static_cast<size_t>(t.rows) * t.cols * D;
                char* taken = arena_set.count(t.block) != 0 ? group_arena_take(raw) : nullptr;
                t.block->ptr = taken != nullptr ? reinterpret_cast<S*>(taken)
                                                : heap_.alloc_resident<S>(raw, ResidentKind::OwnerNear);
                t.block->rows = t.rows;
                t.block->cols = t.cols;
                t.fresh = true;
                ++stats.new_targets;
            }
            if (next_plan) {   // gathered already, under the previous group's GEMMs
                bool same = next_plan->want.size() == want.size();
                for (size_t i = 0; same && i < want.size(); ++i) same = next_plan->want[i].m == want[i].m;
                if (!same) throw std::runtime_error("LevelEliminator: the prefetched gather does not match its group");
                materialize_from_gather(*next_plan, stream, MatGroup, &arena_take);
                next_plan.reset();
            } else {
                materialize_adopted(want, stream, MatGroup, &arena_take);
            }
        }
        lap(stats.og_mat);
        // The group's source slots, decompressed into a window of their own.
        // This is what replaces holding the wave's two n_tot x r outputs raw:
        // a slot is one neighbor's rows, which is exactly what a task reads,
        // and the group's union of them is bounded by kOwnerWindowBytes.
        char* window = nullptr;
        std::unordered_map<uint64_t, size_t> win_at;
        if (from_store) {
            struct Need { int32_t src; int32_t slot; int which; size_t off; };
            std::vector<Need> needs;   // in traversal order, so every rank agrees
            size_t need_bytes = 0;
            for (size_t i = first; i < last; ++i) {
                for (const SrcSlot& ss : targets[i].src_slots) {
                    for (int which = 0; which < 2; ++which) {
                        const int32_t slot = which == 0 ? ss.sa : ss.sb;
                        if (ss.src < 0 || slot < 0) continue;
                        const auto [it, fresh_slot] = win_at.emplace(win_key(ss.src, slot, which), need_bytes);
                        if (!fresh_slot) continue;
                        needs.push_back(Need{ss.src, slot, which, need_bytes});
                        need_bytes += slot_bytes(ss.src, slot, which);
                    }
                }
            }
            window = heap_.try_alloc(std::max<size_t>(need_bytes, 1), ResidentKind::SrcWindow);
            if (window == nullptr) {
                throw std::runtime_error("H2 GPU: no room for the owner pass's source window (" +
                                         std::to_string(need_bytes / (1u << 20)) + " MB)");
            }
            SolveXnrWave w;
            for (const Need& nd : needs) {
                const WaveBox& sb = boxes[static_cast<size_t>(nd.src)];
                const BoxState* sst = sb.state;
                const char* packed = nd.which == 0 ? packed_of(*sst) : sst->xnrt_packed;
                if (nd.which == 0) { if (sst->t2t_packed != nullptr) ++stats.fs_read_dev[0]; else ++stats.fs_read_host[0]; }
                const auto& slots = nd.which == 0 ? sst->t2t_slots : sst->xnrt_slots;
                const auto& strides = nd.which == 0 ? sst->t2t_strides : sst->xnrt_strides;
                const auto& sl = slots[static_cast<size_t>(nd.slot)];
                const char* src = packed + sl.off;
                w.boxes.push_back(static_cast<int>(w.boxes.size()));
                w.h_box_dsts.push_back(window + nd.off);
                for (uint32_t c = 0; c < sl.nchunks; ++c) {
                    w.h_srcs.push_back(src);
                    w.h_dsts.push_back(window + nd.off + static_cast<size_t>(c) * kXnrChunkBytes);
                    src += strides[sl.chunk0 + c];
                }
            }
            if (!w.h_srcs.empty()) {
                if (!build_solve_xnr_wave(w, heap_, stream)) {
                    throw std::runtime_error("H2 GPU: could not plan the owner window's decompression");
                }
                decompress_solve_xnr_wave(w, stream);
                check_cuda(cudaStreamSynchronize(stream), "owner window");
                release_solve_xnr_wave(w);
                heap_.free(reinterpret_cast<char*>(const_cast<void**>(w.d_srcs)));
            }
            stats.ow_bytes += static_cast<double>(need_bytes);
            stats.ow_max = std::max(stats.ow_max, static_cast<double>(need_bytes));
            stats.ow_slots += static_cast<int64_t>(needs.size());
            ++stats.ow_groups;
        }
        lap(stats.og_win);
        owner_meta_.clear();
        ++owner_group_stamp_;
        std::vector<EvalItem> evals;
        int ev_m = 0, ev_n = 0;
        size_t depth = 0;
        for (size_t i = first; i < last; ++i) {
            Target& t = targets[i];
            note_write(t.block, static_cast<size_t>(t.rows) * t.cols * D);
            depth = std::max(depth, t.tasks.size());
            if (!t.fresh) continue;
            evals.push_back(EvalItem{t.block->ptr, t.rows, t.rows, t.cols,
                                     current_slots_in(*t.row_box, owner_meta_, owner_group_stamp_),
                                     current_slots_in(*t.col_box, owner_meta_, owner_group_stamp_)});
            ev_m = std::max(ev_m, t.rows);
            ev_n = std::max(ev_n, t.cols);
        }
        const size_t off_ev = owner_meta_.append(evals);
        std::vector<VBatch<S>> gb(depth);
        for (size_t i = first; i < last; ++i) {
            Target& t = targets[i];
            for (auto& e : t.tasks) e.c = t.block->ptr;
            if (from_store) {
                // the slots are stored transposed, r x rows, so the operands
                // have leading dimension r and the batch runs (Trans, NoTrans)
                if (t.src_slots.size() != t.tasks.size()) {
                    throw std::runtime_error("LevelEliminator: owner target has " +
                                             std::to_string(t.tasks.size()) + " tasks but " +
                                             std::to_string(t.src_slots.size()) + " source slots");
                }
                for (size_t j = 0; j < t.tasks.size(); ++j) {
                    const SrcSlot& ss = t.src_slots[j];
                    const magma_int_t r = static_cast<magma_int_t>(boxes[static_cast<size_t>(ss.src)].r);
                    t.tasks[j].a = reinterpret_cast<const S*>(window + win_at.at(win_key(ss.src, ss.sa, 0)));
                    t.tasks[j].lda = r;
                    t.tasks[j].b = reinterpret_cast<const S*>(window + win_at.at(win_key(ss.src, ss.sb, 1)));
                    t.tasks[j].ldb = r;
                }
            }
            for (size_t j = 0; j < t.tasks.size(); ++j) gb[j].entries.push_back(t.tasks[j]);
        }
        for (auto& b : gb) b.stage(owner_meta_);
        char* gmd = owner_meta_.upload(owner_meta_device_, stream);
        stats.bytes_up += static_cast<double>(owner_meta_.size());
        lap(stats.og_stage);
        check_items(evals, "group evals");
        for (const auto& b : gb) check_batch(b, "group gemms");
        {
            EntryOriginScope origin(EntryOrigin::fill);  // (collective entries' profile)
            eval_blocks(evals, owner_meta_, gmd, off_ev, ev_m, ev_n, stream);
        }
        for (const auto& b : gb) {
            // (a nonsymmetric matrix: temp2 rows times X_RN columns)
            b.gemm(gmd, from_store ? MagmaTrans : MagmaNoTrans,
                   from_store || nonsym_ ? MagmaNoTrans : MagmaTrans, 1.0, 1.0, queue);
            stats.owner_gemms += static_cast<int64_t>(b.count());
            for (const auto& e : b.entries) {
                stats.owner_flops += kFlopScale * 2.0 * e.m * e.n * e.k;
                stats.fl_box[7] += kFlopScale * 2.0 * e.m * e.n * e.k;
            }
            ++stats.owner_batches;
        }
        // The next group's packed blocks gathered now, on the near copy
        // stream, while this group's GEMMs run.  A wave's groups partition
        // its targets, so nothing this group writes is among them.
        if (gi + 1 < owner_groups.size() && kHoldNearBlocksPacked && compress_factors_ >= 1) {
            const auto [nf, nl] = owner_groups[gi + 1];
            GatherPlan np;
            for (size_t i = nf; i < nl; ++i) {
                Target& t = targets[i];
                if (t.block->ptr != nullptr) continue;
                const bool schur = static_cast<int64_t>(t.key) < 0;
                const int64_t morton = schur ? static_cast<int64_t>(~t.key) : 0;
                if (!(schur ? is_packed_schur(morton) : is_packed_edge(t.key))) continue;
                t.block->rows = t.rows;
                t.block->cols = t.cols;
                np.want.push_back(AdoptedBlock{schur, schur ? static_cast<uint64_t>(morton) : t.key, t.block});
            }
            if (!np.want.empty()) {
                if (near_prefetch_done_ == nullptr) near_prefetch_done_ = make_event();
                check_cuda(cudaEventSynchronize(near_prefetch_done_), "prefetched gather");  // (its list is reused)
                plan_gather(np, near_copy_stream(), near_prefetch_done_);
                next_plan = std::move(np);
                ++stats.mat_prefetched;
            }
        }
        check_cuda(cudaStreamSynchronize(stream), "owner group");  // the image is reused
        heap_.free(window);
        lap(stats.og_launch);
        // The policy the store exists for: every near block this group just
        // wrote goes back compressed, so the resident near-field set is the
        // packed one (9.47 -> 3.25 GB at the 192^3 Color leaf) and raw copies
        // live only inside the wave that reads them.
        //
        // This measured negative twice before.  Both times the loss was the
        // piece-4 owner-join pre-pass, which decompressed the whole adopted
        // set up front and so held the raw bytes anyway, defeating the bound
        // the groups exist to impose; with it removed the policy gains 3.60 GB
        // on Color and 2.94 GB on CA0.  What remains is re-packing a block
        // once per writing wave, 8.9 times on average: 70-83% of the bytes
        // have their last write in the final quarter of the waves, so
        // compressing once per block is not available.
        const bool og_pressure = under_pressure();
        if (kHoldNearBlocksPacked && compress_factors_ >= 1 && replicated_) {
            pack_near_by_rule(targets, first, last, stream, og_pressure);
            ++(og_pressure ? stats.np_pressure_groups : stats.np_raw_groups);
        } else if (kHoldNearBlocksPacked && compress_factors_ >= 1 && og_pressure) {
            pack_near_group(targets, first, last, stream);
            ++stats.np_pressure_groups;
        } else if (kHoldNearBlocksPacked && compress_factors_ >= 1) {
            ++stats.np_raw_groups;   // room enough: the group's blocks stay raw
        }
        // measuring only: compress to learn the lengths, keep the full
        // precision blocks, discard the compressed bytes
        if (compress_factors_ >= 1 && env::check("fillsrc")) {
            measure_near_sizes(targets, first, last, stream);
        }
        if (group_arena_.base != nullptr) {
            // what the pack left raw (the pressure foreseen wrongly, or no room
            // in the store) to a block of its own, then the group's block back
            bool moved = false;
            for (size_t i = first; i < last; ++i) {
                Target& t = targets[i];
                if (t.block->ptr == nullptr || !in_group_arena(t.block->ptr)) continue;
                const size_t raw = static_cast<size_t>(t.rows) * t.cols * D;
                S* own = heap_.alloc_resident<S>(raw, t.fresh ? ResidentKind::OwnerNear : ResidentKind::Adopted);
                check_cuda(cudaMemcpyAsync(own, t.block->ptr, raw, cudaMemcpyDeviceToDevice, stream), "owner group block");
                t.block->ptr = own;
                ++stats.og_arena_moved;
                moved = true;
            }
            // (the pack's compression synchronized the stream; the moves have
            // not, and the heap may hand the block to a writer on another
            // stream - the near prefetch - before they run)
            if (moved) check_cuda(cudaStreamSynchronize(stream), "owner group block");
            heap_.free(group_arena_.base);
            group_arena_ = GroupArena{};
        }
        lap(stats.og_pack);
    }
    marks.mark(stream);  // 7: owner GEMMs
    if (env::check("heap")) {
        // A near block whose slot still names memory the heap has given back:
        // the wave read and retired it, nothing replaced the entry, and the
        // level's end will free it a second time.  Report the first few keys
        // so the pairs can be identified rather than guessed at.
        int dangling = 0;
        for (const auto& kv : edges_) {
            if (kv.second.ptr == nullptr || heap_.is_allocated(kv.second.ptr)) continue;
            if (dangling++ < 6) {
                std::printf("[heap] level %d wave %d: edges_(%lld,%lld) names freed memory, %dx%d\n",
                            level_index_, wave_index, static_cast<long long>(kv.first >> 32),
                            static_cast<long long>(kv.first & 0xffffffffu), kv.second.rows, kv.second.cols);
            }
        }
        for (size_t b = 0; b < states_.size(); ++b) {
            const DeviceMatrix& m = states_[b].schur;
            if (m.ptr == nullptr || heap_.is_allocated(m.ptr)) continue;
            if (dangling++ < 10) {
                std::printf("[heap] level %d wave %d: schur of local box %zu names freed memory\n",
                            level_index_, wave_index, b);
            }
        }
        if (dangling > 0) {
            std::printf("[heap] level %d wave %d: %d dangling block slots\n", level_index_, wave_index, dangling);
            std::fflush(stdout);
        }
    }
    // ---- 4. factors to the host, in the background: the result block (LU,
    // pivots, temp1, X_RS, T; LU status checked on arrival) and the source
    // groups (temp2 as X_NR, X_RR_full).
    {
        auto ready = [&] {
            cudaEvent_t ev = nullptr;
            check_cuda(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming), "cudaEventCreate");
            check_cuda(cudaEventRecord(ev, stream), "cudaEventRecord");
            return ev;
        };
        using Seg = HostCopier::Segment;
        const bool with_t = gpu_sketch_;
        HostCopier::Job result;
        result.device = d_result;
        result.bytes = result_bytes;
        result.ready = ready();
        auto info_host = std::make_shared<std::vector<int>>();
        result.segments.push_back(Seg{nullptr, info_host.get(), info_offset, active.size() * sizeof(int)});
        struct Dims { Box* box; int k, r; bool host; };
        auto dims = std::make_shared<std::vector<Dims>>();
        for (int i : active) {
            const WaveBox& wb = boxes[static_cast<size_t>(i)];
            Box* box = wb.box;
            const size_t k = static_cast<size_t>(wb.k), r = static_cast<size_t>(wb.r);
            // (ghosts too: the CA solve applies them again)
            // A kept box's X_SR and LU stay on the device only, as its X_NR
            // (the device solve and logdet read them there;
            // materialize_host_factors brings them back for a host use;
            // rescue_skipped_xnr when the keep fails); its X_RS goes unread
            // (the host transport, which reads it, is not on).  T stays on
            // the host for the matvec.
            // An owner level's owned ghost: its factors go to its home rank
            // (the FULL payload), nothing stays here; and a replicated level's
            // ghost under H2_CA_solve home (solved from home factors).
            const bool owned_ghost = ghost_host_copies_spent(box->morton_index);
            const bool host = !owned_ghost && !(skip_xnr_ && !keep_failed_ && kept_.boxes.count(box->morton_index) != 0);
            dims->push_back({box, wb.k, wb.r, host});
            box->X_RR = MatrixStorage<DataType>{};
            box->X_SR = MatrixStorage<DataType>{};
            box->X_RS_entry = MatrixStorage<DataType>{};
            if (nonsym_) box->X_RS = MatrixStorage<DataType>{};
            if (with_t) box->interpolation_matrix = MatrixStorage<DataType>{};
            if (owned_ghost) continue;
            if (host) {
                result.segments.push_back(Seg{&box->X_RR.data, nullptr, wb.off_xrr, r * r * D});
                result.segments.push_back(Seg{nullptr, &box->X_RR_pivots, wb.off_piv, r * sizeof(int)});
                result.segments.push_back(Seg{&box->X_SR.data, nullptr, wb.off_temp1, k * r * D});
                result.segments.push_back(Seg{&box->X_RS_entry.data, nullptr, wb.off_xrs, r * k * D});
                // (a nonsymmetric matrix: temp3, the solve's X_RS)
                if (nonsym_) result.segments.push_back(Seg{&box->X_RS.data, nullptr, wb.off_temp3, r * k * D});
            } else {
                box->X_RR_pivots.clear();
                skipped_factors_.push_back(box);
            }
            if (with_t) result.segments.push_back(Seg{&box->interpolation_matrix.data, nullptr, wb.off_t, k * r * D});
        }
        result.finalize = [info_host, dims, with_t, nonsym = nonsym_] {
            for (size_t j = 0; j < dims->size(); ++j) {
                const Dims& d = (*dims)[j];
                if ((*info_host)[j] != 0) {
                    throw std::runtime_error("LevelEliminator: LU of X_RR failed with INFO = " +
                                             std::to_string((*info_host)[j]) + " for box " +
                                             std::to_string(d.box->morton_index));
                }
                auto shape = [](MatrixStorage<DataType>& m, int rows, int cols, typename MatrixStorage<DataType>::Format f) {
                    m.rows = rows;
                    m.cols = cols;
                    m.lda = rows;
                    m.format = f;
                };
                if (d.host) {
                    shape(d.box->X_RR, d.r, d.r, MatrixStorage<DataType>::LU_FACTORED);
                    shape(d.box->X_SR, d.k, d.r, MatrixStorage<DataType>::FULL);
                    shape(d.box->X_RS_entry, d.r, d.k, MatrixStorage<DataType>::FULL);
                    if (nonsym) shape(d.box->X_RS, d.r, d.k, MatrixStorage<DataType>::FULL);
                }
                if (with_t) shape(d.box->interpolation_matrix, d.k, d.r, MatrixStorage<DataType>::FULL);
            }
        };
        copier_->submit(std::move(result));
        for (size_t g = 0; g < d_sources.size(); ++g) {
            HostCopier::Job job;
            job.device = d_sources[g];
            job.bytes = source_group_bytes[g];
            job.ready = ready();
            // temp2 is a block of its own now, so the X_NR the host still
            // wants is a second job; when nothing wants it, the block was
            // given back before the owner pass and this job stays empty.
            HostCopier::Job t2job;
            t2job.device = d_temp2[g];
            t2job.bytes = temp2_group_bytes[g];
            struct Shape { Box* box; int ntot, r; bool xrr; };
            auto shapes = std::make_shared<std::vector<Shape>>();
            for (int i : active) {
                const WaveBox& wb = boxes[static_cast<size_t>(i)];
                if (wb.source_group != static_cast<int>(g)) continue;
                const size_t nt = static_cast<size_t>(wb.ntot), r = static_cast<size_t>(wb.r);
                wb.box->X_NR = MatrixStorage<DataType>{};
                wb.box->X_RR_full = MatrixStorage<DataType>{};
                if (nonsym_) {
                    // temp4, the solve's X_RN (beside temp2 in its group)
                    wb.box->X_RN = MatrixStorage<DataType>{};
                    if (nt > 0 && d_temp2[g] != nullptr) {
                        t2job.segments.push_back(Seg{&wb.box->X_RN.data, nullptr, wb.off_temp4, r * nt * D});
                    }
                }
                // (a kept box's X_NR stays on the device only: rescue_skipped_xnr)
                // (an owned ghost keeps its host copies: a fill source the
                // transition may restore, not a solve factor here - its X_NR
                // only where its temp2^T is not in the packed store)
                const bool skip = host_skips_xnr(wb);
                const bool ghost_drop = host_drops_ghost_xnr(wb);
                if (skip) skipped_xnr_.push_back(wb.box);
                const bool skip_xrr = host_skips_xrr(wb);  // (H2 host memory: rebuilt from the kept LU instead)
                if (nt > 0 && !skip && !ghost_drop) {
                    if (d_temp2[g] == nullptr) {
                        // its chunk gave the raw form back before the keep
                        // failed; rescue_skipped_xnr reads the store instead
                        skipped_xnr_.push_back(wb.box);
                    } else {
                        t2job.segments.push_back(Seg{&wb.box->X_NR.data, nullptr, wb.off_temp2, nt * r * D});
                    }
                }
                if (!skip_xrr) job.segments.push_back(Seg{&wb.box->X_RR_full.data, nullptr, wb.off_xrr_full, r * r * D});
                if (wb.has_remote && nt > 0 && !device_exchange_) {
                    job.segments.push_back(Seg{&wb.box->lazy_original_x_nr, nullptr, wb.off_xnr_orig, nt * r * D});
                }
                shapes->push_back({wb.box, skip || ghost_drop ? 0 : wb.ntot, wb.r, !skip_xrr});
            }
            job.finalize = [shapes, nonsym = nonsym_] {
                for (const Shape& sh : *shapes) {
                    if (sh.ntot > 0) {
                        sh.box->X_NR.rows = sh.ntot;
                        sh.box->X_NR.cols = sh.r;
                        sh.box->X_NR.lda = sh.ntot;
                        sh.box->X_NR.format = MatrixStorage<DataType>::FULL;
                        if (nonsym) {
                            sh.box->X_RN.rows = sh.r;
                            sh.box->X_RN.cols = sh.ntot;
                            sh.box->X_RN.lda = sh.r;
                            sh.box->X_RN.format = MatrixStorage<DataType>::FULL;
                        }
                    }
                    if (!sh.xrr) continue;
                    sh.box->X_RR_full.rows = sh.r;
                    sh.box->X_RR_full.cols = sh.r;
                    sh.box->X_RR_full.lda = sh.r;
                    sh.box->X_RR_full.format = MatrixStorage<DataType>::FULL;
                }
            };
            copier_->submit(std::move(job));
            if (!t2job.segments.empty()) {
                t2job.ready = ready();
                t2job.finalize = [] {};  // (the shapes are set by the source job's)
                temp2_to_copier[g] = 1;
                copier_->submit(std::move(t2job));
            }
        }
    }
    for (int i : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(i)];
        Box* box = wb.box;
        box->deferred_xnn_neighbor_point_counts = wb.counts;
        std::vector<DataType>().swap(box->deferred_xnn_temp2);
        if (is_local(box->morton_index)) {  // (a ghost's solve data lives with its owner)
            level_.solve_neighbor_size[static_cast<size_t>(box->morton_index - level_.local_morton_start)] = wb.counts;
        }
        if (gpu_sketch_) {
            BoxState& st = *wb.state;
            if (packed_of(st) == nullptr) st.temp2t = wb.persist_temp2t;  // else stored compressed
            if (nonsym_) st.temp4 = wb.persist_temp4;
            st.xrr_full = wb.persist_xrr;
            st.ntot = wb.ntot;
            st.r = wb.r;
            // Its LU is kept, so the stored X_RR_full is redundant and goes
            // back here - r x r a box, 0.90 GB of the 192^3 Color leaf's
            // resident set.  Only when the box owns the block: without
            // early_free_sources_ it belongs to level_allocs_.
            // (not for a nonsymmetric matrix, whose kept factors may go mid-level: see the reclaimer)
            if (!nonsym_ && st.xrr_lu != nullptr && st.persist != nullptr &&
                st.persist == reinterpret_cast<char*>(wb.persist_xrr)) {
                auto& es = eliminator_stats();
                es.persist_bytes -= static_cast<double>(st.persist_bytes);
                heap_.free(st.persist);
                st.persist = nullptr;
                st.persist_bytes = 0;
                st.xrr_full = nullptr;
                ++es.xr_dropped;
            }
        }
    }
    stats.launch += std::chrono::duration<double>(clock::now() - t0).count();
    t0 = clock::now();
    if (!gpu_sketch_) {
        free_copied_blocks(copier_->wait_all());  // the host sketch of the next wave reads them
    } else {
        free_copied_blocks(copier_->take_finished());
    }
    for (WaveBox& wb : boxes) {
        BoxState& st = *wb.state;
        st.eliminated = true;
        st.skeleton.assign(wb.box->skeleton_indices.begin(), wb.box->skeleton_indices.end());
        st.skeleton_identity = true;
        for (size_t i = 0; i < st.skeleton.size(); ++i) {
            if (st.skeleton[i] != static_cast<int>(i)) {
                st.skeleton_identity = false;
                break;
            }
        }
        if (st.skeleton_identity && static_cast<int>(st.skeleton.size()) != st.n) st.skeleton_identity = false;
    }
    {
        auto& es = eliminator_stats();
        es.lsch_union_sum += static_cast<double>(wave_schur_bytes);
        es.lsch_union_max_wave = std::max(es.lsch_union_max_wave, static_cast<double>(wave_schur_bytes));
        es.ledge_union_sum += static_cast<double>(wave_edge_bytes);
        es.ledge_union_max_wave = std::max(es.ledge_union_max_wave, static_cast<double>(wave_edge_bytes));
        es.lsch_pairs += static_cast<int64_t>(wave_schur.size());
        es.ledge_pairs += static_cast<int64_t>(wave_edge.size());
        for (const void* q : wave_schur) lsch_seen_.insert(q);
        for (const void* q : wave_edge) ledge_seen_.insert(q);
        es.lsch_distinct = static_cast<int64_t>(lsch_seen_.size());
        es.ledge_distinct = static_cast<int64_t>(ledge_seen_.size());
    }
    for (char* work : work_blocks) heap_.free(work);
    for (char* b : xnr_blocks) heap_.free(b);  // (empty when they went before the owner pass)
    for (size_t g = 0; g < d_temp2.size(); ++g) {
        if (d_temp2[g] != nullptr && !temp2_to_copier[g]) heap_.free(d_temp2[g]);
    }
    heap_.free(d_t);
    if (early_free_sources_) {
        // the wave's sketches were the last reads of sources whose local
        // neighbors are now all eliminated: the heap may take them back
        for (const WaveBox& wb : boxes) {
            for (int64_t m : wb.box->one_hop) {
                BoxState& es = state_of(m);
                if (es.persist != nullptr && --es.readers == 0) released_sources_.push_back({&es, m});
            }
        }
    }
    if (wave_sketch_ != nullptr) {
        heap_.free(wave_sketch_);
        wave_sketch_ = nullptr;
    }
    for (char* p : wave_sketch_extra_) heap_.free(p);
    wave_sketch_extra_.clear();
    release_owner_sources();  // read by this wave's owner pass only
    if (under_pressure()) {   // the packed temp2^T to the host and off the device, when the room is needed
        stage_fill_copies();
        evict_fill_copies(false);
    }
    if (owner_ != nullptr) {  // (this wave's payload goes out now that its copies are done, not at the next wave)
        owner_progress();
        owner_install_arrived(wave_index / std::max(owner_->num_waves, 1));
    }
    stats.store += std::chrono::duration<double>(clock::now() - t0).count();
    if (env::trace("wave")) {   // the host side, wave by wave (BPACK_TRACE=wave)
        double rss = 0.0, hwm = 0.0;
        host_memory_gb(rss, hwm);
        double sends = 0.0, recvs = 0.0;   // (an owner level: the payloads still held)
        for (const OwnerWaveSend& w : owner_sends_) {
            for (const auto& m : w.messages) sends += static_cast<double>(m.bytes);
        }
        for (const OwnerWaveRecv& r : owner_recvs_) recvs += static_cast<double>(r.bytes);
        std::printf("  [gpu] level %d wave %d: host RSS %.2f GB (high-water %.2f), heap used %.2f GB, pinned chunks %.2f GB "
                    "(%.2f idle)%s\n", level_index_, last_wave_index_, rss, hwm, static_cast<double>(heap_.used()) / 1e9,
                    static_cast<double>(pinned_chunk_pool().pinned()) / 1e9, static_cast<double>(pinned_chunk_pool().idle()) / 1e9,
                    owner_ != nullptr ? (", owner payloads held: sends " + std::to_string(sends / 1e9).substr(0, 5) +
                                         " GB in " + std::to_string(owner_sends_.size()) + " waves, receives " +
                                         std::to_string(recvs / 1e9).substr(0, 5) + " GB in " +
                                         std::to_string(owner_recvs_.size())).c_str() : "");
        std::printf("  [gpu] level %d wave %d: %zu boxes in %.3f s of host time: sketch %.3f (of which waiting %.3f), plan %.3f, "
                    "launch %.3f, store %.3f\n", level_index_, wave_index, wave.size(),
                    std::chrono::duration<double>(clock::now() - t_wave0).count(), stats.sketch - w_sketch0,
                    stats.sk_wait - w_wait0, stats.plan - w_plan0, stats.launch - w_launch0, stats.store - w_store0);
        std::fflush(stdout);
    }
    if (heap_.peak() > stats.heap_peak) {
        stats.heap_peak = heap_.peak();
        stats.heap_peak_resident = heap_.resident_at_peak();
        stats.heap_peak_transient = heap_.transient_at_peak();
    }
    return boundary_count;
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::free_copied_blocks(std::vector<char*>&& blocks) {
    for (char* block : blocks) free_block(block);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::collect_elimination_marks() {
    if (elim_marks_.empty()) return;
    auto& stats = eliminator_stats();
    for (const auto& [m, base] : elim_marks_) {
        if (!m || m->count() < 2) continue;
        m->synchronize();
        for (size_t i = 0; i + 1 < m->count(); ++i) {
            const size_t phase = static_cast<size_t>(base) + i;
            if (phase >= 8) break;
            const double t = m->seconds(i);
            stats.el[phase] += t;
            stats.device += t;
        }
    }
    elim_marks_.clear();
}

// ---------------------------------------------------------------------------
// Level end: Schur and near-field blocks back into the host BoxData.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::download_blocks_to_host() {
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    auto& stats = eliminator_stats();

#ifdef H2_HAVE_NVCOMP
    // The download skips a null ptr, so anything still packed would be
    // silently lost.  Near blocks are held packed between groups now, so this
    // puts the remainder back rather than refusing.
    if (!packed_schur_.empty() || !packed_edges_.empty()) {
        std::vector<AdoptedBlock> want;
        for (const auto& [morton, pb] : packed_schur_) {
            if (!is_interest(morton)) continue;
            BoxState& st = state_of(morton);
            if (st.schur.ptr == nullptr) {
                want.push_back(AdoptedBlock{true, static_cast<uint64_t>(morton), &st.schur});
            }
        }
        for (const auto& [key, pb] : packed_edges_) {
            auto it = edges_.find(key);
            if (it != edges_.end() && it->second.ptr == nullptr) {
                want.push_back(AdoptedBlock{false, key, &it->second});
            }
        }
        materialize_adopted(want, Context::instance().stream(), MatDownload);
    }
#endif
    struct Piece { const DeviceMatrix* m; size_t offset; };
    std::vector<Piece> pieces;
    size_t total = 0;
    for (const BoxState& st : states_) {
        if (st.schur.ptr) { pieces.push_back({&st.schur, total}); total += static_cast<size_t>(st.schur.rows) * st.schur.cols; }
    }
    // (CA levels: pairs of two ghosts stay behind, no local box keeps them)
    std::vector<std::pair<uint64_t, const DeviceMatrix*>> edge_list;
    edge_list.reserve(edges_.size());
    for (const auto& kv : edges_) {
        const int64_t lo = static_cast<int64_t>(kv.first >> 32), hi = static_cast<int64_t>(kv.first & 0xffffffffu);
        if (ca_ && !is_local(lo) && !is_local(hi)) continue;
        edge_list.push_back({kv.first, &kv.second});
    }
    std::sort(edge_list.begin(), edge_list.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& kv : edge_list) {
        if (kv.second->ptr == nullptr) throw std::runtime_error("LevelEliminator: unallocated edge at level end");
        pieces.push_back({kv.second, total});
        total += static_cast<size_t>(kv.second->rows) * kv.second->cols;
    }
    if (total == 0) return;

    // pack on the device with one batched copy, then one download
    S* d_pack = heap_.alloc<S>(total * sizeof(S));
    {
        std::vector<GatherItem> copies;
        int max_m = 0, max_n = 0;
        for (const Piece& p : pieces) {
            copies.push_back(GatherItem{d_pack + p.offset, p.m->rows, p.m->rows, p.m->cols, p.m->ptr, 1, p.m->rows,
                                        IndexList{}, IndexList{}});
            max_m = std::max(max_m, p.m->rows);
            max_n = std::max(max_n, p.m->cols);
        }
        meta_.clear();
        const size_t off = meta_.append(copies);
        char* md = meta_.upload(meta_device_, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off), static_cast<int>(copies.size()), max_m, max_n,
                      md, stream);
    }
    DataType* h_pack = static_cast<DataType*>(pinned_pool().result.reserve(total * sizeof(S)));
    check_cuda(cudaMemcpyAsync(h_pack, d_pack, total * sizeof(S), cudaMemcpyDeviceToHost, stream), "download blocks");
    check_cuda(cudaStreamSynchronize(stream), "download blocks");
    stats.bytes_down += static_cast<double>(total * sizeof(S));

    // host copies in parallel, then the (serial) map updates
    std::vector<std::vector<DataType>> data(pieces.size());
    #pragma omp parallel for schedule(dynamic, 64)
    for (int64_t i = 0; i < static_cast<int64_t>(pieces.size()); ++i) {
        const Piece& p = pieces[static_cast<size_t>(i)];
        data[static_cast<size_t>(i)].assign(h_pack + p.offset,
                                            h_pack + p.offset + static_cast<size_t>(p.m->rows) * p.m->cols);
    }
    size_t piece = 0;
    for (size_t b = 0; b < states_.size(); ++b) {
        const BoxState& st = states_[b];
        if (!st.schur.ptr) continue;
        level_.local_boxes[b].schur_complement.set_owned(st.schur.rows, st.schur.cols, std::move(data[piece++]),
                                                         MatrixStorage<DataType>::FULL);
    }
    for (const auto& kv : edge_list) {
        std::vector<DataType>& block_data = data[piece++];
        if (nonsym_) {
            // K(row, col): the row box's A_SN[col] and the column box's A_NS[row]
            const int64_t col = static_cast<int64_t>(kv.first >> 32);
            const int64_t row = static_cast<int64_t>(kv.first & 0xffffffffu);
            const DeviceMatrix& m = *kv.second;
            Box* row_box = level_.find_local_box(row);
            Box* col_box = level_.find_local_box(col);
            if (row_box == nullptr && col_box == nullptr) {
                throw std::runtime_error("LevelEliminator: near block without a local box at level end");
            }
            if (col_box != nullptr) {
                ModifiedBlock<DataType> view;
                view.neighbor_morton = row;
                std::vector<DataType> copy = row_box != nullptr ? block_data : std::vector<DataType>{};
                view.set_a_ns_owned(m.rows, m.cols, row_box != nullptr ? std::move(copy) : std::move(block_data),
                                    MatrixStorage<DataType>::FULL);
                col_box->near_field_interaction_map[row] =
                    static_cast<int64_t>(col_box->near_field_modified_interactions.size());
                col_box->near_field_modified_interactions.push_back(std::move(view));
            }
            if (row_box != nullptr) {
                ModifiedBlock<DataType> view;
                view.neighbor_morton = col;
                view.A_SN.set_owned(m.rows, m.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
                row_box->near_field_interaction_map_nonsymmetry[col] =
                    static_cast<int64_t>(row_box->near_field_modified_interactions.size());
                row_box->near_field_modified_interactions.push_back(std::move(view));
            }
            continue;
        }
        const int64_t lo = static_cast<int64_t>(kv.first >> 32);
        const int64_t hi = static_cast<int64_t>(kv.first & 0xffffffffu);
        const DeviceMatrix& m = *kv.second;
        Box* lo_box = level_.find_local_box(lo);
        Box* hi_box = level_.find_local_box(hi);
        if (lo_box == nullptr || hi_box == nullptr) {
            // a pair with a box of another rank: this rank's copy, as the
            // local box's view (rows = the other box's points)
            if (lo_box == nullptr && hi_box == nullptr) {
                throw std::runtime_error("LevelEliminator: edge without a local box at level end");
            }
            ModifiedBlock<DataType> view;
            if (lo_box != nullptr) {
                view.neighbor_morton = hi;
                view.set_a_ns_owned(m.rows, m.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
                lo_box->near_field_interaction_map[hi] = static_cast<int64_t>(lo_box->near_field_modified_interactions.size());
                lo_box->near_field_modified_interactions.push_back(std::move(view));
            } else {
                std::vector<DataType> transposed(block_data.size());
                for (int j = 0; j < m.cols; ++j)
                    for (int i = 0; i < m.rows; ++i)
                        transposed[j + static_cast<size_t>(i) * m.cols] = block_data[i + static_cast<size_t>(j) * m.rows];
                view.neighbor_morton = lo;
                view.set_a_ns_owned(m.cols, m.rows, std::move(transposed), MatrixStorage<DataType>::FULL);
                hi_box->near_field_interaction_map[lo] = static_cast<int64_t>(hi_box->near_field_modified_interactions.size());
                hi_box->near_field_modified_interactions.push_back(std::move(view));
            }
            continue;
        }
        // lo's view (rows = hi, cols = lo) is the physical block; hi's view
        // shares it transposed (share_symmetric_level_edges below).
        ModifiedBlock<DataType> lo_block;
        lo_block.neighbor_morton = hi;
        lo_block.set_a_ns_owned(m.rows, m.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
        lo_box->near_field_interaction_map[hi] = static_cast<int64_t>(lo_box->near_field_modified_interactions.size());
        lo_box->near_field_modified_interactions.push_back(std::move(lo_block));
        ModifiedBlock<DataType> hi_block;
        hi_block.neighbor_morton = lo;
        hi_box->near_field_interaction_map[lo] = static_cast<int64_t>(hi_box->near_field_modified_interactions.size());
        hi_box->near_field_modified_interactions.push_back(std::move(hi_block));
    }
    if (!nonsym_) share_symmetric_level_edges(level_);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::finish() {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    for (const BoxState& st : states_) {
        if (!st.eliminated) throw std::runtime_error("LevelEliminator::finish: a box was not eliminated");
    }
    auto& stats = eliminator_stats();
    if (env::trace("wave")) {   // (8x diagnostics) the near store's live bytes against what its host region holds
        double live = 0.0, runs = 0.0;
        size_t blocks = 0;
        auto count = [&](const PackedBlock& pb) {
            if (!pb.live) return;
            ++blocks;
            if (pb.packed != nullptr) { runs += static_cast<double>(pb.bytes); return; }
            for (uint32_t c = 0; c < pb.nchunks; ++c) live += static_cast<double>(chunk_homes_[pb.home0 + c].total);
        };
        for (const auto& kv : packed_schur_) count(kv.second);
        for (const auto& kv : packed_edges_) count(kv.second);
        double held = 0.0;
        for (const ChunkHome& h : chunk_homes_) held += static_cast<double>(h.here);
        int rk = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rk);
        std::printf("  [gpu] level %d rank %d near store at the level end: %zu live blocks, %.2f GB live in the store "
                    "(+%.2f GB in runs) against %.2f GB of homes, %.2f GB compact extents, %zu tail extents; host region %.2f GB\n",
                    level_index_, rk, blocks, live / 1e9, runs / 1e9, held / 1e9, static_cast<double>(compact_total_) / 1e9,
                    tail_extents_.size(), static_cast<double>(near_host_arena(level_index_).used()) / 1e9);
        std::fflush(stdout);
    }
    if (owner_ != nullptr) {  // the last payloads out, the receive buffers back
        check_cuda(cudaStreamSynchronize(Context::instance().stream()), "owner level: last sends");
        owner_progress();
        for (OwnerWaveSend& w : owner_sends_) {
            if (!w.issued) throw std::runtime_error("LevelEliminator::finish: an owned wave's payload was never sent");
            MPI_Waitall(static_cast<int>(w.reqs.size()), w.reqs.data(), MPI_STATUSES_IGNORE);
            if (w.block != nullptr) free_block(w.block);
            cudaEventDestroy(w.ready);
        }
        owner_sends_.clear();
        if (!owner_recvs_.empty()) throw std::runtime_error("LevelEliminator::finish: payloads of an owner level not installed");
        owner_progress();  // (the pending frees: the stream is idle)
        if (!owner_pending_free_.empty()) throw std::runtime_error("LevelEliminator::finish: receive buffers still pending");
    }
    {
        const auto t1 = clock::now();
        free_copied_blocks(copier_->wait_all());
        if (exchange_copier_) free_copied_blocks(exchange_copier_->wait_all());
        collect_elimination_marks();
        // host copies of the remote generators (restored fill sources, host transition)
        for (auto& [m, rh] : remote_host_) {
            Box* gb = level_.find_generator_box(m);
            if (gb == nullptr || (!packed_restore_ && rh->temp2.size() != static_cast<size_t>(rh->ntot) * rh->r) ||
                rh->xrr.size() != static_cast<size_t>(rh->r) * rh->r) {
                throw std::runtime_error("LevelEliminator::finish: host copy of generator " + std::to_string(m) +
                                         " incomplete");
            }
            if (!packed_restore_) gb->X_NR.set_owned(rh->ntot, rh->r, std::move(rh->temp2), MatrixStorage<DataType>::FULL);
            gb->X_RR_full.set_owned(rh->r, rh->r, std::move(rh->xrr), MatrixStorage<DataType>::FULL);
            if (nonsym_) {
                // X_RN holds temp4, as a finalized source does (install_remote_generators)
                if (rh->temp4.size() != static_cast<size_t>(rh->ntot) * rh->r) {
                    throw std::runtime_error("LevelEliminator::finish: host copy of generator " + std::to_string(m) +
                                             " has no temp4");
                }
                gb->X_RN.set_owned(rh->r, rh->ntot, std::move(rh->temp4), MatrixStorage<DataType>::FULL);
            }
        }
        remote_host_.clear();
        stats.finish_sources += std::chrono::duration<double>(clock::now() - t1).count();
        stats.finish_store = copier_->busy_seconds();
        stats.copier_wait = copier_->wait_seconds();
        stats.bytes_down += copier_->bytes();
        if (exchange_copier_) stats.bytes_down += exchange_copier_->bytes();
        heap_.set_reclaimer(nullptr);
        stats.heap_reclaims = heap_.reclaims();
        copier_.reset();
        exchange_copier_.reset();
    }
    if (heap_.peak() > stats.heap_peak) {
        stats.heap_peak = heap_.peak();
        stats.heap_peak_resident = heap_.resident_at_peak();
        stats.heap_peak_transient = heap_.transient_at_peak();
    }
    compact_near_store();   // before the transition packs the parents beside it (into the freed chunks)
    stats.finish += std::chrono::duration<double>(clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// Multi-rank levels.  Boxes of other ranks within two hops (the level's
// assisting boxes) get states and point slots; their skeleton, and so their
// current rows, arrive with the assisting data of a transport.  Remote
// sources arrive as generators (lazy Schur mode 2): their updates of this
// rank's blocks are formed here on the device (the host's
// form_near_updates_from_generators + apply_updates_with_kernel_symmetric),
// and they stay resident as fill sources of later sketches.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::set_remote_skeletons(
    const std::vector<std::pair<BoxState*, const std::vector<int64_t>*>>& skeletons) {
    std::vector<int> all;
    std::vector<std::pair<BoxState*, size_t>> listed;
    for (const auto& [st, skeleton] : skeletons) {
        if (st->eliminated) continue;
        st->eliminated = true;
        st->skeleton.assign(skeleton->begin(), skeleton->end());
        st->skeleton_identity = static_cast<int>(st->skeleton.size()) == st->n;
        for (size_t i = 0; st->skeleton_identity && i < st->skeleton.size(); ++i) {
            st->skeleton_identity = st->skeleton[i] == static_cast<int>(i);
        }
        if (st->skeleton_identity) continue;
        listed.emplace_back(st, all.size());
        all.insert(all.end(), st->skeleton.begin(), st->skeleton.end());
    }
    if (all.empty()) return;
    int* d_lists = heap_.alloc_resident<int>(all.size() * sizeof(int), ResidentKind::Lists);
    level_allocs_.push_back(reinterpret_cast<char*>(d_lists));
    check_cuda(cudaMemcpyAsync(d_lists, all.data(), all.size() * sizeof(int), cudaMemcpyHostToDevice,
                               Context::instance().stream()),
               "upload remote skeletons");
    for (const auto& [st, start] : listed) {
        st->d_skeleton = d_lists + start;
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::sync_remote_boxes(bool skeletons) {
    cudaStream_t stream = Context::instance().stream();
    auto& stats = eliminator_stats();
    std::vector<std::pair<int64_t, int64_t>> fresh;  // Morton, assisting index
    for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) {
        if (remote_index_.count(kv.first) == 0) fresh.emplace_back(kv.first, kv.second);
    }
    std::sort(fresh.begin(), fresh.end());
    if (!fresh.empty()) {
        int64_t extra = 0;
        for (const auto& [morton, index] : fresh) {
            const auto& assist = level_.assisting_boxes[static_cast<size_t>(index)];
            if (assist.indices.empty() ||
                assist.coords.size() != static_cast<size_t>(tree_->dimension) * assist.indices.size()) {
                throw std::runtime_error("LevelEliminator: assisting box " + std::to_string(morton) +
                                         " has no point data");
            }
            extra += static_cast<int64_t>(assist.indices.size());
        }
        if (num_slots_ + extra >= std::numeric_limits<int>::max()) {
            throw std::runtime_error("LevelEliminator: too many points per rank");
        }
        const int64_t dim = tree_->dimension;
        std::vector<double> xyz(static_cast<size_t>(dim * extra));
        std::vector<int64_t> ids(static_cast<size_t>(extra));
        int64_t at = 0;
        for (const auto& [morton, index] : fresh) {
            const auto& assist = level_.assisting_boxes[static_cast<size_t>(index)];
            BoxState st;
            st.remote = true;
            st.slot = static_cast<int>(num_slots_ + at);
            st.n = static_cast<int>(assist.indices.size());
            st.on_boundary = assist.on_boundary;
            for (size_t i = 0; i < assist.indices.size(); ++i) {
                for (int64_t d = 0; d < dim; ++d) {
                    xyz[static_cast<size_t>(dim * (at + static_cast<int64_t>(i)) + d)] =
                        static_cast<double>(assist.coords[static_cast<size_t>(dim) * i + static_cast<size_t>(d)]);
                }
                ids[static_cast<size_t>(at) + i] = assist.indices[i];
            }
            remote_index_[morton] = remote_states_.size();
            remote_states_.push_back(std::move(st));
            at += static_cast<int64_t>(assist.indices.size());
        }
        // a larger table: the local slots and earlier remote ones keep their places
        const int64_t total = num_slots_ + extra;
        double* xyz_new = heap_.alloc_resident<double>(static_cast<size_t>(dim * total) * sizeof(double), ResidentKind::Points);
        int64_t* ids_new = heap_.alloc_resident<int64_t>(static_cast<size_t>(total) * sizeof(int64_t), ResidentKind::Points);
        check_cuda(cudaMemcpyAsync(xyz_new, d_xyz_, static_cast<size_t>(dim * num_slots_) * sizeof(double),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        check_cuda(cudaMemcpyAsync(ids_new, d_ids_, static_cast<size_t>(num_slots_) * sizeof(int64_t),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        check_cuda(cudaMemcpyAsync(xyz_new + dim * num_slots_, xyz.data(), xyz.size() * sizeof(double),
                                   cudaMemcpyHostToDevice, stream), "remote points");
        check_cuda(cudaMemcpyAsync(ids_new + num_slots_, ids.data(), ids.size() * sizeof(int64_t),
                                   cudaMemcpyHostToDevice, stream), "remote ids");
        heap_.free(d_xyz_);  // later launches are ordered after the copies
        heap_.free(d_ids_);
        d_xyz_ = xyz_new;
        d_ids_ = ids_new;
        points_.xyz = xyz_new;
        points_.ids = ids_new;
        num_slots_ = total;
        stats.bytes_up += static_cast<double>(xyz.size() * sizeof(double) + ids.size() * sizeof(int64_t));
    }
    // boxes eliminated since the last transport (their skeleton came with it)
    if (!skeletons) return;
    std::vector<std::pair<BoxState*, const std::vector<int64_t>*>> fresh_skeletons;
    for (const auto& [morton, index] : remote_index_) {
        BoxState& st = remote_states_[index];
        if (st.eliminated || st.ghost) continue;
        const auto& assist = level_.assisting_boxes[static_cast<size_t>(
            level_.assisting_box_points_for_kernel_evaluation.at(morton))];
        if (!assist.skel_indices.empty()) fresh_skeletons.emplace_back(&st, &assist.skel_indices);
    }
    set_remote_skeletons(fresh_skeletons);
}

// Far points sampled by a round of the adaptive ID rows (H2_ID_proxy 2), by
// global index: their coordinates (the tree's global ones, dim per point; an
// evaluator that reads only the ids needs none) and ids go to the scratch
// slots, which a later round overwrites.  Returns the first slot.
template<typename CoordType, typename DataType, typename KernelType>
int LevelEliminator<CoordType, DataType, KernelType>::place_scratch_points(const std::vector<int64_t>& global_ids,
                                                                           cudaStream_t stream) {
    const int64_t count = static_cast<int64_t>(global_ids.size());
    const int64_t dim = tree_->dimension;
    if (count > scratch_slots_) {
        // a larger table with a new scratch region at its end: the slots in
        // use keep their places (the old region stays unused)
        const int64_t slots = count + count / 2;
        const int64_t total = num_slots_ + slots;
        if (total >= std::numeric_limits<int>::max()) throw std::runtime_error("LevelEliminator: too many points per rank");
        double* xyz_new = heap_.alloc_resident<double>(static_cast<size_t>(dim * total) * sizeof(double), ResidentKind::Points);
        int64_t* ids_new = heap_.alloc_resident<int64_t>(static_cast<size_t>(total) * sizeof(int64_t), ResidentKind::Points);
        check_cuda(cudaMemcpyAsync(xyz_new, d_xyz_, static_cast<size_t>(dim * num_slots_) * sizeof(double),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        check_cuda(cudaMemcpyAsync(ids_new, d_ids_, static_cast<size_t>(num_slots_) * sizeof(int64_t),
                                   cudaMemcpyDeviceToDevice, stream), "point table");
        heap_.free(d_xyz_);  // later launches are ordered after the copies
        heap_.free(d_ids_);
        d_xyz_ = xyz_new;
        d_ids_ = ids_new;
        points_.xyz = xyz_new;
        points_.ids = ids_new;
        scratch_slot0_ = static_cast<int>(num_slots_);
        scratch_slots_ = static_cast<int>(slots);
        num_slots_ = total;
    }
    const auto& coords = tree_->id_source_point_coords;
    const bool have_coords = coords.size() == static_cast<size_t>(tree_->num_points) * static_cast<size_t>(dim);
    scratch_xyz_.assign(static_cast<size_t>(dim * count), 0.0);
    if (have_coords) {
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < count; ++i) {
            const size_t at = static_cast<size_t>(global_ids[static_cast<size_t>(i)]) * static_cast<size_t>(dim);
            for (int64_t d = 0; d < dim; ++d) {
                scratch_xyz_[static_cast<size_t>(dim * i + d)] = static_cast<double>(coords[at + static_cast<size_t>(d)]);
            }
        }
    }
    check_cuda(cudaMemcpyAsync(d_xyz_ + dim * static_cast<int64_t>(scratch_slot0_), scratch_xyz_.data(),
                               scratch_xyz_.size() * sizeof(double), cudaMemcpyHostToDevice, stream), "sampled points");
    check_cuda(cudaMemcpyAsync(d_ids_ + scratch_slot0_, global_ids.data(), global_ids.size() * sizeof(int64_t),
                               cudaMemcpyHostToDevice, stream), "sampled ids");
    eliminator_stats().bytes_up += static_cast<double>(scratch_xyz_.size() * sizeof(double) +
                                                      global_ids.size() * sizeof(int64_t));
    return scratch_slot0_;
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::receive_remote(const std::vector<int64_t>& installed) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    sync_remote_boxes();
    apply_remote_generators(installed);
    eliminator_stats().remote += std::chrono::duration<double>(clock::now() - t0).count();
}

// For each installed source E (wave, then Morton order) and each local
// neighbor L of E (rows of E's temp2 at the counts of E's elimination):
//   Schur(L)  += temp2_L X_NR_L^T
//   edge(L,E) := A(E skeleton, L) + X_RS^T temp2_L^T        (E's rows shrink)
//   edge(L,N) += temp2_N X_NR_L^T  for N one hop from L (a local N once, from
//                                  its lower endpoint)
// with X_NR the original coupling.  Deltas of a block from several sources
// are summed in source order before they are added, as on the host.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::apply_remote_generators(
    const std::vector<int64_t>& installed, const std::unordered_map<int64_t, Incoming>* incoming) {
    if (installed.empty()) return;
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    const size_t D = sizeof(S);

    std::vector<int64_t> order(installed);
    std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
        const int32_t wa = level_.elimination_wave.at(a), wb = level_.elimination_wave.at(b);
        return wa != wb ? wa < wb : a < b;
    });

    struct Gen {
        int64_t morton = -1;
        Box* box = nullptr;
        BoxState* st = nullptr;
        int r = 0, k = 0, ntot = 0;
        std::vector<int> offset;  // first temp2 row of each one-hop neighbor
        std::vector<int> xnr_offset;  // first X_NR row of each local one-hop neighbor
        int xnr_ld = 0;
        const S* temp2 = nullptr;
        const S* xnr = nullptr;
        const S* xrs = nullptr;
        // a nonsymmetric matrix: the original X_RN (r x ntot), temp4 (r x
        // ntot) and the solved X_SR, temp1 (k x r)
        const S* xrn = nullptr;
        const S* temp4 = nullptr;
        const S* xsr_solved = nullptr;
        const char* bulk = nullptr;  // received into device memory (GenLayout)
    };
    if (nonsym_ && incoming == nullptr) {
        throw std::runtime_error("LevelEliminator: a nonsymmetric matrix's generators come through the device exchange");
    }
    std::vector<Gen> gens;
    std::vector<std::pair<BoxState*, const std::vector<int64_t>*>> skeletons;
    // (an owner level: the installed boxes are this rank's ghosts, or its
    // own boxes another rank eliminated; the generator's rows and updates
    // cover the boxes of interest here, not only the local ones)
    auto interest = [&](int64_t m) { return owner_ != nullptr ? is_interest(m) : is_local(m); };
    for (int64_t m : order) {
        Box* gb = owner_ != nullptr ? held_box(m) : level_.find_generator_box(m);
        if (gb == nullptr) throw std::runtime_error("LevelEliminator: generator " + std::to_string(m) + " missing");
        BoxState& st = state_of(m);
        if (owner_ == nullptr && !st.remote) throw std::runtime_error("LevelEliminator: generator of local box " + std::to_string(m));
        if (!gb->skeleton_indices.empty()) skeletons.emplace_back(&st, &gb->skeleton_indices);
        Gen g;
        g.morton = m;
        g.box = gb;
        g.st = &st;
        ++stats.remote_generators;
        const auto& counts = gb->deferred_xnn_neighbor_point_counts;
        if (incoming != nullptr) {  // dimensions from the header, data in the receive buffer
            auto it = incoming->find(m);
            if (it == incoming->end()) continue;  // full rank: nothing to add, no fill
            g.r = it->second.r;
            g.k = it->second.k;
            g.ntot = it->second.ntot;
            g.xnr_ld = it->second.nx;
            g.bulk = it->second.bulk;
            if (counts.size() != gb->one_hop.size() || g.k != static_cast<int>(gb->skeleton_indices.size())) {
                throw std::runtime_error("LevelEliminator: incomplete generator " + std::to_string(m));
            }
        } else {
            g.r = gb->X_NR.is_allocated() ? static_cast<int>(gb->X_NR.cols) : 0;
            if (g.r == 0) continue;  // full rank: nothing to add, no fill
            g.ntot = static_cast<int>(gb->X_NR.rows);
            g.k = gb->X_RS.is_allocated() ? static_cast<int>(gb->X_RS.cols) : 0;
            if (counts.size() != gb->one_hop.size() || gb->X_NR.lda != gb->X_NR.rows ||
                gb->lazy_original_x_nr.size() != static_cast<size_t>(g.ntot) * g.r ||
                !gb->X_RR_full.is_allocated() || gb->X_RR_full.rows != g.r || gb->X_RR_full.cols != g.r ||
                (g.k > 0 && (gb->X_RS.rows != g.r || gb->X_RS.lda != g.r)) ||
                g.k != static_cast<int>(gb->skeleton_indices.size())) {
                throw std::runtime_error("LevelEliminator: incomplete generator " + std::to_string(m));
            }
        }
        g.offset.resize(counts.size());
        int rows = 0;
        for (size_t i = 0; i < counts.size(); ++i) {
            g.offset[i] = rows;
            rows += static_cast<int>(counts[i]);
        }
        if (rows != g.ntot) throw std::runtime_error("LevelEliminator: generator row counts " + std::to_string(m));
        if (owner_ != nullptr) {
            // (the generator's view of its neighbors of interest here must be
            // this rank's: the boxes of no interest it saw eliminated may
            // have no state here)
            for (size_t a = 0; a < counts.size(); ++a) {
                const int64_t nm = gb->one_hop[a];
                if (counts[a] == 0 || !is_interest(nm)) continue;
                const BoxState& ns = state_of(nm);
                if (current_size(ns) != static_cast<int>(counts[a])) {
                    const auto nw = level_.elimination_wave.find(nm);
                    throw std::runtime_error("LevelEliminator: installed generator " + std::to_string(m) + " (wave " +
                                             std::to_string(level_.elimination_wave.at(m)) + ") saw its neighbor " +
                                             std::to_string(nm) + " (" + (is_local(nm) ? "local" : ns.ghost ? "ghost" : "assisting") +
                                             (is_writable(nm) ? ", eliminated here" : "") + ") at " + std::to_string(counts[a]) +
                                             " rows; here it has " + std::to_string(current_size(ns)) + " (" +
                                             (ns.eliminated ? "eliminated in wave " + std::to_string(nw != level_.elimination_wave.end() ? nw->second : -1) : "not eliminated") +
                                             ")");
                }
            }
        }
        if (g.bulk != nullptr && nonsym_) {  // (no X_NR: X_RN whole)
            g.xnr_offset = g.offset;
            g.xnr_ld = g.ntot;
        } else if (g.bulk != nullptr) {  // X_NR: rows of the local neighbors only
            g.xnr_offset.assign(counts.size(), -1);
            int local_rows = 0;
            for (size_t a = 0; a < counts.size(); ++a) {
                if (counts[a] == 0 || !interest(gb->one_hop[a])) continue;
                g.xnr_offset[a] = local_rows;
                local_rows += static_cast<int>(counts[a]);
            }
            if (local_rows != g.xnr_ld) {
                throw std::runtime_error("LevelEliminator: generator " + std::to_string(m) + " X_NR rows mismatch");
            }
        } else {
            g.xnr_offset = g.offset;
            g.xnr_ld = g.ntot;
        }
        gens.push_back(std::move(g));
    }
    set_remote_skeletons(skeletons);
    if (gens.empty()) return;

    // uploads (unless received into device memory): [temp2 | original X_NR
    // | X_RS] for this pass, and [temp2^T | X_RR_full], kept as a fill
    // source of later sketches
    // (one sub-buffer per generator, in pieces: PiecedBuffer)
    size_t tmp_bytes = 0;
    std::vector<size_t> tmp_off(gens.size()), tmp_starts;
    for (size_t i = 0; i < gens.size(); ++i) {
        if (gens[i].bulk != nullptr) continue;
        const size_t nr = static_cast<size_t>(gens[i].ntot) * gens[i].r * D;
        tmp_off[i] = tmp_bytes;
        tmp_starts.push_back(tmp_bytes);
        tmp_bytes = align_up(align_up(tmp_bytes + nr) + nr) + align_up(static_cast<size_t>(gens[i].r) * gens[i].k * D);
    }
    PiecedBuffer tmp;
    tmp.allocate(heap_, tmp_starts, tmp_bytes, ResidentKind::GenUpload);
    std::vector<TransposeItem> transposes;
    // With compression on, the generators' temp2^T is transposed into one
    // staging a group at a time and packed group by group, instead of a raw
    // block each until a single pass at the end.  A group's transposes are
    // contiguous in the list, one a generator in gens order.
    struct GenGroup { std::vector<size_t> gens; std::vector<size_t> off; size_t bytes = 0; size_t first_tr = 0; };
    std::vector<GenGroup> gen_groups;
    size_t gen_staging_max = 0;
    int max_tr_m = 0, max_tr_n = 0;
    for (size_t i = 0; i < gens.size(); ++i) {
        Gen& g = gens[i];
        Box* gb = g.box;
        const size_t nr = static_cast<size_t>(g.ntot) * g.r * D;
        const char* xrr_device = nullptr;
        if (g.bulk != nullptr) {
            const GenLayout lay(static_cast<size_t>(g.ntot), static_cast<size_t>(g.r), static_cast<size_t>(g.k),
                                static_cast<size_t>(nonsym_ ? 0 : g.xnr_ld), nonsym_);
            g.temp2 = reinterpret_cast<const S*>(g.bulk);
            g.xnr = reinterpret_cast<const S*>(g.bulk + lay.xnr);
            g.xrs = reinterpret_cast<const S*>(g.bulk + lay.xrs);
            xrr_device = g.bulk + lay.xrr;
            if (nonsym_) {
                g.xrn = reinterpret_cast<const S*>(g.bulk + lay.xrn);
                g.temp4 = reinterpret_cast<const S*>(g.bulk + lay.temp4);
                g.xsr_solved = reinterpret_cast<const S*>(g.bulk + lay.xsr_solved);
            }
        } else {
            char* base = tmp.at(tmp_off[i]);
            g.temp2 = reinterpret_cast<const S*>(base);
            g.xnr = reinterpret_cast<const S*>(base + align_up(nr));
            g.xrs = reinterpret_cast<const S*>(base + align_up(align_up(nr) + nr));
            // pageable copies: the host data is staged before each call returns
            check_cuda(cudaMemcpyAsync(base, gb->X_NR.data.data(), nr, cudaMemcpyHostToDevice, stream),
                       "generator temp2");
            check_cuda(cudaMemcpyAsync(base + align_up(nr), gb->lazy_original_x_nr.data(), nr, cudaMemcpyHostToDevice,
                                       stream), "generator X_NR");
            if (g.k > 0) {
                check_cuda(cudaMemcpyAsync(base + align_up(align_up(nr) + nr), gb->X_RS.data.data(),
                                           static_cast<size_t>(g.r) * g.k * D, cudaMemcpyHostToDevice, stream),
                           "generator X_RS");
            }
            std::vector<DataType>().swap(gb->lazy_original_x_nr);  // needed only here (as on the host)
            stats.bytes_up += static_cast<double>(2 * nr + static_cast<size_t>(g.r) * (g.r + g.k) * D);
        }
        const size_t t_bytes = align_up(nr);
        const size_t xrr_bytes = static_cast<size_t>(g.r) * g.r * D;
        char* persist = heap_.alloc_resident(std::max<size_t>(xrr_bytes, 1), ResidentKind::FillXrr);
        // temp2^T: raw only until its slots are packed, so with compression on
        // it goes into the group staging below rather than a block of its own.
        // One block a generator, all waiting on a single pass at the end, was
        // 1.18 GB of Color's level 4 at the peak.
        char* persist_t2t = nullptr;
        if (compress_factors_ < 1) {
            persist_t2t = heap_.alloc_resident(std::max<size_t>(t_bytes, 1), ResidentKind::FillRaw);
        } else {
            if (gen_groups.empty() || (gen_groups.back().bytes > 0 &&
                                       gen_groups.back().bytes + t_bytes > kFillStagingBudget)) {
                gen_groups.push_back(GenGroup{});
                gen_groups.back().first_tr = transposes.size();
            }
            GenGroup& gg = gen_groups.back();
            gg.gens.push_back(i);
            gg.off.push_back(gg.bytes);
            gg.bytes += t_bytes;
            gen_staging_max = std::max(gen_staging_max, gg.bytes);
        }
        {
            auto& es = eliminator_stats();
            const size_t raw = persist_t2t != nullptr ? t_bytes : 0;
            es.persist_bytes += static_cast<double>(xrr_bytes + raw);
            es.persist_gen_bytes += static_cast<double>(xrr_bytes + raw);
            es.persist_peak = std::max(es.persist_peak, es.persist_bytes);
            ++es.persist_blocks;
        }
        keep_fill_source(*g.st, persist, xrr_bytes, persist_t2t, persist_t2t != nullptr ? t_bytes : 0, gb->one_hop);
        if (owner_ != nullptr) g.st->readers = std::numeric_limits<int>::max() / 2;  // (no host copy to restore it from)
        if (nonsym_) {
            // temp4, the fill's column side, as the source's own rank keeps it
            const size_t t4_bytes = std::max<size_t>(nr, 1);
            char* persist_t4 = heap_.alloc_resident(t4_bytes, ResidentKind::FillRaw);
            if (early_free_sources_) {
                g.st->persist_t4 = persist_t4;
                g.st->persist_t4_bytes = t4_bytes;
            } else {
                level_allocs_.push_back(persist_t4);
            }
            if (nr > 0) {
                check_cuda(cudaMemcpyAsync(persist_t4, g.temp4, nr, cudaMemcpyDeviceToDevice, stream), "generator temp4");
            }
            g.st->temp4 = reinterpret_cast<const S*>(persist_t4);
            auto& es = eliminator_stats();
            es.persist_bytes += static_cast<double>(t4_bytes);
            es.persist_gen_bytes += static_cast<double>(t4_bytes);
            es.persist_peak = std::max(es.persist_peak, es.persist_bytes);
        }
        g.st->temp2t = reinterpret_cast<const S*>(persist_t2t);
        g.st->xrr_full = reinterpret_cast<const S*>(persist);
        g.st->ntot = g.ntot;
        g.st->r = g.r;
        if (xrr_device != nullptr) {
            check_cuda(cudaMemcpyAsync(persist, xrr_device, static_cast<size_t>(g.r) * g.r * D,
                                       cudaMemcpyDeviceToDevice, stream), "generator X_RR_full");
        } else {
            check_cuda(cudaMemcpyAsync(persist, gb->X_RR_full.data.data(),
                                       static_cast<size_t>(g.r) * g.r * D, cudaMemcpyHostToDevice, stream),
                       "generator X_RR_full");
        }
        transposes.push_back(TransposeItem{g.temp2, g.ntot, reinterpret_cast<S*>(persist_t2t), g.r, g.ntot, g.r});
        max_tr_m = std::max(max_tr_m, g.ntot);
        max_tr_n = std::max(max_tr_n, g.r);
    }

    ++wave_stamp_;
    meta_.clear();
    std::vector<EvalItem> evals;
    std::vector<GatherItem> gathers;
    int max_ev_m = 0, max_ev_n = 0, max_ga_m = 0, max_ga_n = 0;
    auto add_eval = [&](S* out, int ld, int m, int n, IndexList rows, IndexList cols) {
        evals.push_back(EvalItem{out, ld, m, n, rows, cols});
        max_ev_m = std::max(max_ev_m, m);
        max_ev_n = std::max(max_ev_n, n);
    };
    VBatch<S> rep_tt, rep_nn;  // replaced blocks: + X_RS^T temp2_L^T (or its transpose)
    struct Product { size_t off; const S* a; int lda; const S* b; int ldb; int m, n, k; };
    std::vector<Product> products;  // a b^T into the work block
    struct Accum { S* target; int ld, rows, cols; std::vector<size_t> parts; };
    std::vector<Accum> accums;
    std::unordered_map<uint64_t, size_t> accum_of;  // edge key, or ~Morton for a Schur block
    std::unordered_set<uint64_t> replaced;
    std::vector<char*> retired;
    size_t work_bytes = 0;

    // The targets a chunk makes from the kernel (no stored block yet): packed
    // at the chunk's end with the adopted ones, so they live in the store
    // from the start.  Left raw, they were 2.95 GB of the Color leaf's peak
    // at 192^3 (34987 targets of 7.75 GB over the level) and 0.59 of CA3's.
    std::vector<AdoptedBlock> fresh_made;
    auto ensure_schur = [&](int64_t morton, BoxState& st) -> DeviceMatrix& {
        if (st.schur.ptr == nullptr) {
            const int c = current_size(st);
            st.schur = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(c) * c * D, ResidentKind::LevelSchur), c, c};
            add_eval(st.schur.ptr, c, c, c, current_slots(st), current_slots(st));
            fresh_made.push_back(AdoptedBlock{true, static_cast<uint64_t>(morton), &st.schur});
        }
        return st.schur;
    };
    auto ensure_edge = [&](int64_t lo, int64_t hi) -> DeviceMatrix& {
        const uint64_t key = edge_key(lo, hi);
        auto it = edges_.find(key);
        if (it != edges_.end()) {
            if (it->second.ptr == nullptr) assert_not_packed(false, key, "ensure_edge");
            return it->second;
        }
        BoxState& ls = state_of(lo);
        BoxState& hs = state_of(hi);
        const int rows = current_size(hs), cols = current_size(ls);
        DeviceMatrix& e = edges_[key];
        e = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D, ResidentKind::LevelSchur), rows, cols};
        add_eval(e.ptr, rows, rows, cols, current_slots(hs), current_slots(ls));
        fresh_made.push_back(AdoptedBlock{false, key, &e});
        return e;
    };
    // (a nonsymmetric matrix: the directed block K(row, col))
    auto ensure_dblock = [&](int64_t row, int64_t col) -> DeviceMatrix& {
        const uint64_t key = dkey(row, col);
        auto it = edges_.find(key);
        if (it != edges_.end()) {
            if (it->second.ptr == nullptr) assert_not_packed(false, key, "ensure_dblock");
            return it->second;
        }
        BoxState& rs = state_of(row);
        BoxState& cs = state_of(col);
        const int rows = current_size(rs), cols = current_size(cs);
        DeviceMatrix& e = edges_[key];
        e = DeviceMatrix{heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D, ResidentKind::LevelSchur), rows, cols};
        add_eval(e.ptr, rows, rows, cols, current_slots(rs), current_slots(cs));
        fresh_made.push_back(AdoptedBlock{false, key, &e});
        return e;
    };
    auto accumulate = [&](uint64_t key, const DeviceMatrix& target, int rows, int cols, const S* a, int lda,
                          const S* b, int ldb, int k) {
        note_write(&target, static_cast<size_t>(target.rows) * target.cols * D);
        if (target.rows != rows || target.cols != cols) {
            const int64_t lo = static_cast<int64_t>(key >> 32), hi = static_cast<int64_t>(key & 0xffffffffu);
            std::ostringstream oss;
            oss << "LevelEliminator: remote update of block (" << lo << ", " << hi << ") is " << rows << " x " << cols
                << " but the block is " << target.rows << " x " << target.cols;
            if (static_cast<int64_t>(key) < 0) {
                oss << " (a Schur block)";
            } else {
                oss << " (" << lo << ": " << (is_local(lo) ? "local" : "ghost") << (is_writable(lo) ? " eliminated here" : "")
                    << (state_of(lo).eliminated ? " eliminated" : "") << " size " << current_size(state_of(lo)) << "; " << hi
                    << ": " << (is_local(hi) ? "local" : "ghost") << (is_writable(hi) ? " eliminated here" : "")
                    << (state_of(hi).eliminated ? " eliminated" : "") << " size " << current_size(state_of(hi)) << ")";
            }
            throw std::runtime_error(oss.str());
        }
        // (work_bytes is what staging these would have cost, for the report)
        work_bytes = align_up(work_bytes + static_cast<size_t>(rows) * cols * D);
        products.push_back(Product{0, a, lda, b, ldb, rows, cols, k});
        auto it = accum_of.find(key);
        if (it == accum_of.end()) {
            it = accum_of.emplace(key, accums.size()).first;
            accums.push_back(Accum{target.ptr, target.rows, rows, cols, {}});
        }
        accums[it->second].parts.push_back(products.size() - 1);
    };

    // An owner level: a block is kept up to date here only if it is read
    // here: a Schur block by its box's elimination here or, a box homed here,
    // by the transition; a pair by an endpoint's elimination here, or by the
    // transition when an endpoint is homed here.  (An owned ghost's other
    // blocks are spent after its elimination here; a box of no interest may
    // have no state.)  One definition, so the pass below and the
    // materialization above it cannot drift apart.
    auto block_spent = [&](int64_t lm, const BoxState& ls) {
        return owner_ != nullptr && !is_local(lm) && ls.eliminated;
    };
    auto pair_read = [&](int64_t lm, const BoxState& ls, int64_t nm) {
        if (owner_ == nullptr) return true;
        if (is_local(lm) || is_local(nm) || !ls.eliminated) return true;
        return is_writable(nm) && !state_of(nm).eliminated;
    };
    auto schur_written_here = [&](int64_t lm, const BoxState& ls) {
        return (is_writable(lm) && !block_spent(lm, ls)) ||
               (owner_ != nullptr && is_local(lm) && ls.eliminated);
    };

#ifdef H2_HAVE_NVCOMP
    // H2_compress_factors: the adopted blocks one generator's updates read -
    // ensure_schur's Schur blocks as well as the near blocks replace_edge()
    // and ensure_edge() gather from - as candidates, less what `seen` holds.
    // The interest gate is the same one the loop applies, so this names what
    // the generator reads and not more.  A chunk of generators is
    // materialized in one batch (run_fill_batch() plans and synchronizes
    // twice per call, so a call per pair is not an option).
    // `fresh` counts what the generator would make from the kernel (a target
    // with no stored block, and every replace_edge): those are raw from the
    // chunk's start to its pack as well, so the chunk's budget covers them
    // too.  A block two generators of a chunk would both make is counted
    // twice, which only makes the chunk smaller.
    auto want_of_gen = [&](const Gen& g, const std::unordered_set<const void*>& seen,
                           std::vector<AdoptedBlock>& cand, size_t& fresh) {
        const auto& hop = g.box->one_hop;
        const auto& counts = g.box->deferred_xnn_neighbor_point_counts;
        std::unordered_set<const void*> mine;   // (a block twice in one generator)
        auto want_edge = [&](int64_t a, int64_t b) {
            const uint64_t key = edge_key(a, b);
            auto it = edges_.find(key);
            if (it == edges_.end() || it->second.ptr != nullptr) return;
            if (is_packed_edge(key) && seen.count(&it->second) == 0 && mine.insert(&it->second).second) {
                cand.push_back(AdoptedBlock{false, key, &it->second});
            }
        };
        auto replace_bytes = [&](size_t li) {   // replace_edge's fresh block for pair (hop[li], g)
            const bool e_hi = g.morton > hop[li];
            const size_t c_l = static_cast<size_t>(counts[li]), k = static_cast<size_t>(g.k);
            return (e_hi ? k : c_l) * (e_hi ? c_l : k) * D;
        };
        for (size_t li = 0; li < hop.size(); ++li) {
            if (counts[li] == 0) continue;
            const int64_t lm = hop[li];
            if (owner_ == nullptr && !has_state(lm)) continue;  // (see process_gen)
            BoxState& ls = state_of(lm);
            if (!interest(lm)) {
                // replace_edge(false) for a box of no interest, on an owner level
                if (owner_ != nullptr && g.k > 0 && is_interest(g.morton)) {
                    want_edge(lm, g.morton);
                    fresh += replace_bytes(li);
                }
                continue;
            }
            // ensure_schur(ls)
            if (schur_written_here(lm, ls) && ls.schur.ptr == nullptr) {
                if (is_packed_schur(lm)) {
                    if (seen.count(&ls.schur) == 0 && mine.insert(&ls.schur).second) {
                        cand.push_back(AdoptedBlock{true, static_cast<uint64_t>(lm), &ls.schur});
                    }
                } else {
                    const size_t c = static_cast<size_t>(current_size(ls));
                    fresh += c * c * D;
                }
            }
            // replace_edge(false)
            if (g.k > 0 && pair_read(lm, ls, g.morton)) {
                want_edge(lm, g.morton);
                fresh += replace_bytes(li);
            }
            // ensure_edge(min, max) of the inner loop
            for (size_t ni = 0; ni < hop.size(); ++ni) {
                if (ni == li || counts[ni] == 0) continue;
                const int64_t nm = hop[ni];
                if (!deferred_xnn_boxes_are_one_hop(level_.dimension, lm, nm)) continue;
                if (interest(nm) && nm < lm) continue;
                if (!pair_read(lm, ls, nm)) continue;
                want_edge(lm, nm);
                if (edges_.find(edge_key(lm, nm)) == edges_.end()) {
                    const BoxState& hs = state_of(std::max(lm, nm));
                    const BoxState& los = state_of(std::min(lm, nm));
                    fresh += static_cast<size_t>(current_size(hs)) * static_cast<size_t>(current_size(los)) * D;
                }
            }
        }
    };
#endif
    // One generator's updates, into the chunk's batches.
    auto process_gen = [&](const Gen& g) {
        const auto& hop = g.box->one_hop;
        const auto& counts = g.box->deferred_xnn_neighbor_point_counts;
        for (size_t li = 0; li < hop.size(); ++li) {
            const int64_t lm = hop[li];
            const int c_l = static_cast<int>(counts[li]);
            if (c_l == 0) continue;
            // A Color level updates nothing of a neighbor of no interest here,
            // and it may have no state here at all: the unstructured exchange
            // sends a generator to every rank holding it as an assisting box
            // (the lazy requesters), not only to its neighbors' ranks.  On a
            // full grid every neighbor of a received generator is known.
            if (owner_ == nullptr && !has_state(lm)) continue;
            BoxState& ls = state_of(lm);
            if (current_size(ls) != c_l) {
                const auto gw = level_.elimination_wave.find(g.morton);
                const auto lw = level_.elimination_wave.find(lm);
                throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) + " (wave " +
                                         std::to_string(gw != level_.elimination_wave.end() ? gw->second : -1) +
                                         ") sees box " + std::to_string(lm) + " at another size: " +
                                         std::to_string(c_l) + " against " + std::to_string(current_size(ls)) +
                                         " here (eliminated in wave " +
                                         std::to_string(lw != level_.elimination_wave.end() ? lw->second : -1) +
                                         "; the last wave is " + std::to_string(last_wave_index_) + ")");
            }
            const S* t2_l = g.temp2 + g.offset[li];
            // edge (L, E) := A(E skeleton, L) + X_RS^T temp2_L^T: E's rows shrink
            // to its skeleton (a stored block gathers its rows; none: the kernel)
            auto replace_edge = [&](bool require_stored) {
                const uint64_t key = edge_key(lm, g.morton);
                const bool e_hi = g.morton > lm;
                const int rows = e_hi ? g.k : c_l, cols = e_hi ? c_l : g.k;
                auto it = edges_.find(key);
                if (it != edges_.end() && it->second.ptr == nullptr) {
                    assert_not_packed(false, key, "the stored-block pass");
                }
                if (it == edges_.end() && require_stored) return;
                S* fresh = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D, ResidentKind::LevelSchur);
                IndexList skeleton = current_slots(*g.st);
                skeleton.base = 0;  // positions within E
                if (it != edges_.end()) {
                    const DeviceMatrix& old = it->second;
                    if (old.rows != (e_hi ? g.st->n : c_l) || old.cols != (e_hi ? c_l : g.st->n)) {
                        throw std::runtime_error("LevelEliminator: stored block of remote source " + std::to_string(g.morton) +
                                                 " (wave " + std::to_string(level_.elimination_wave.at(g.morton)) + ", n " +
                                                 std::to_string(g.st->n) + ", k " + std::to_string(g.k) + ") with box " +
                                                 std::to_string(lm) + " (" + (is_local(lm) ? "local" : "ghost") +
                                                 (is_writable(lm) ? ", eliminated here" : "") + ", count " + std::to_string(c_l) +
                                                 ", size " + std::to_string(current_size(ls)) + ") is " + std::to_string(old.rows) +
                                                 " x " + std::to_string(old.cols) + " (rows = the higher box)");
                    }
                    gathers.push_back(e_hi ? GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, skeleton, IndexList{}}
                                           : GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, IndexList{}, skeleton});
                    max_ga_m = std::max(max_ga_m, rows);
                    max_ga_n = std::max(max_ga_n, cols);
                    retired.push_back(reinterpret_cast<char*>(old.ptr));
                } else if (e_hi) {
                    add_eval(fresh, rows, rows, cols, current_slots(*g.st), current_slots(ls));
                } else {
                    add_eval(fresh, rows, rows, cols, current_slots(ls), current_slots(*g.st));
                }
                if (e_hi) {
                    rep_tt.entries.push_back({g.xrs, t2_l, fresh, g.k, c_l, g.r, g.r, g.ntot, rows});
                } else {
                    rep_nn.entries.push_back({t2_l, g.xrs, fresh, c_l, g.k, g.r, g.ntot, g.r, rows});
                }
                edges_[key] = DeviceMatrix{fresh, rows, cols};
                fresh_made.push_back(AdoptedBlock{false, key, &edges_[key]});
                note_write(&edges_[key], static_cast<size_t>(rows) * cols * D);
                replaced.insert(key);
            };
            if (!interest(lm)) {
                // (an owner level: this rank's own box eliminated elsewhere also
                // gets its pairs with the boxes of no interest here, as its
                // owner's elimination made them; the transition reads them)
                if (owner_ != nullptr && g.k > 0 && is_interest(g.morton)) replace_edge(false);
                continue;
            }
            const S* x_l = g.xnr + g.xnr_offset[li];
            // An owner level: a block is kept up to date here only if it is
            // read here: a Schur block by its box's elimination here or, a
            // box homed here, by the transition; a pair by an endpoint's
            // elimination here, or by the transition when an endpoint is
            // homed here.  (An owned ghost's other blocks are spent after
            // its elimination here; a box of no interest may have no state.)
            const bool spent = block_spent(lm, ls);
            auto pair_read_here = [&](int64_t nm) { return pair_read(lm, ls, nm); };
            if (schur_written_here(lm, ls)) {
                // ensure_schur() makes a fresh block from the kernel on a null
                // ptr, so an adopted block still packed here is a wrong answer
                // rather than an error
                if (ls.schur.ptr == nullptr) assert_not_packed(true, static_cast<uint64_t>(lm), "ensure_schur");
                accumulate(~static_cast<uint64_t>(lm), ensure_schur(lm, ls), c_l, c_l, t2_l, g.ntot, x_l, g.xnr_ld, g.r);
            }
            if (g.k > 0 && pair_read_here(g.morton)) replace_edge(false);

            for (size_t ni = 0; ni < hop.size(); ++ni) {
                if (ni == li) continue;
                const int64_t nm = hop[ni];
                const int c_n = static_cast<int>(counts[ni]);
                if (c_n == 0 || !deferred_xnn_boxes_are_one_hop(level_.dimension, lm, nm)) continue;
                if (interest(nm) && nm < lm) continue;  // from the lower endpoint of interest
                if (!pair_read_here(nm)) continue;
                if (current_size(state_of(nm)) != c_n) {
                    const auto gw = level_.elimination_wave.find(g.morton);
                    const auto nw = level_.elimination_wave.find(nm);
                    throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) + " (wave " +
                                             std::to_string(gw != level_.elimination_wave.end() ? gw->second : -1) +
                                             ") sees box " + std::to_string(nm) + " at another size: " +
                                             std::to_string(c_n) + " against " +
                                             std::to_string(current_size(state_of(nm))) + " here (" +
                                             (is_local(nm) ? "local" : "remote") + ", eliminated in wave " +
                                             std::to_string(nw != level_.elimination_wave.end() ? nw->second : -1) +
                                             ", next to local box " + std::to_string(lm) + "; the last wave is " +
                                             std::to_string(last_wave_index_) + ")");
                }
                const uint64_t key = edge_key(lm, nm);
                if (replaced.count(key)) throw std::runtime_error("LevelEliminator: a replaced block is also updated");
                const DeviceMatrix& target = ensure_edge(std::min(lm, nm), std::max(lm, nm));
                const S* t2_n = g.temp2 + g.offset[ni];
                if (lm < nm) {
                    accumulate(key, target, c_n, c_l, t2_n, g.ntot, x_l, g.xnr_ld, g.r);  // rows: N
                } else {
                    accumulate(key, target, c_l, c_n, x_l, g.xnr_ld, t2_n, g.ntot, g.r);  // rows: L
                }
            }
        }
    };

    // A nonsymmetric matrix's generator (the host's
    // form_near_updates_from_generators_nonsymmetric), all products NN, with
    // X_RN the source's original coupling (r x ntot) and temp1 its solved
    // X_SR; for each neighbor L of interest here:
    //   Schur(L) += temp2_L X_RN[:, L]
    //   K(B, L)  := K(S_B, L) + temp1 X_RN[:, L]     (k x c_l; B's rows shrink)
    //   K(L, B)  := K(L, S_B) + temp2_L X_RS         (c_l x k; its columns)
    //   K(L, N)  += temp2_L X_RN[:, N],  K(N, L) += temp2_N X_RN[:, L]
    //               for N one hop from L (a pair of interest once, from its lower box)
    auto process_gen_nonsym = [&](const Gen& g) {
        const auto& hop = g.box->one_hop;
        const auto& counts = g.box->deferred_xnn_neighbor_point_counts;
        const size_t r = static_cast<size_t>(g.r);
        for (size_t li = 0; li < hop.size(); ++li) {
            const int64_t lm = hop[li];
            const int c_l = static_cast<int>(counts[li]);
            if (c_l == 0 || !has_state(lm) || !interest(lm)) continue;  // (see process_gen)
            BoxState& ls = state_of(lm);
            if (current_size(ls) != c_l) {
                throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) + " sees box " +
                                         std::to_string(lm) + " at another size: " + std::to_string(c_l) +
                                         " against " + std::to_string(current_size(ls)) + " here");
            }
            const S* t2_l = g.temp2 + g.offset[li];
            const S* xrn_l = g.xrn + static_cast<size_t>(g.offset[li]) * r;
            if (ls.schur.ptr == nullptr) assert_not_packed(true, static_cast<uint64_t>(lm), "ensure_schur");
            accumulate(~static_cast<uint64_t>(lm), ensure_schur(lm, ls), c_l, c_l, t2_l, g.ntot, xrn_l, g.r, g.r);
            if (g.k > 0) {
                IndexList skeleton = current_slots(*g.st);
                skeleton.base = 0;  // positions within B
                // b_rows: K(B, L), B's rows at its skeleton; else K(L, B), its columns
                auto replace_view = [&](bool b_rows) {
                    const uint64_t key = b_rows ? dkey(g.morton, lm) : dkey(lm, g.morton);
                    const int rows = b_rows ? g.k : c_l, cols = b_rows ? c_l : g.k;
                    auto it = edges_.find(key);
                    if (it != edges_.end() && it->second.ptr == nullptr) assert_not_packed(false, key, "the stored-block pass");
                    S* fresh = heap_.alloc_resident<S>(static_cast<size_t>(rows) * cols * D, ResidentKind::LevelSchur);
                    if (it != edges_.end()) {
                        const DeviceMatrix& old = it->second;
                        if (old.rows != (b_rows ? g.st->n : c_l) || old.cols != (b_rows ? c_l : g.st->n)) {
                            throw std::runtime_error("LevelEliminator: stored block K(" +
                                                     std::to_string(b_rows ? g.morton : lm) + ", " +
                                                     std::to_string(b_rows ? lm : g.morton) + ") of remote source " +
                                                     std::to_string(g.morton) + " is " + std::to_string(old.rows) +
                                                     " x " + std::to_string(old.cols));
                        }
                        gathers.push_back(b_rows ? GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, skeleton, IndexList{}}
                                                 : GatherItem{fresh, rows, rows, cols, old.ptr, 1, old.rows, IndexList{}, skeleton});
                        max_ga_m = std::max(max_ga_m, rows);
                        max_ga_n = std::max(max_ga_n, cols);
                        retired.push_back(reinterpret_cast<char*>(old.ptr));
                    } else if (b_rows) {
                        add_eval(fresh, rows, rows, cols, current_slots(*g.st), current_slots(ls));
                    } else {
                        add_eval(fresh, rows, rows, cols, current_slots(ls), current_slots(*g.st));
                    }
                    if (b_rows) {
                        rep_nn.entries.push_back({g.xsr_solved, xrn_l, fresh, g.k, c_l, g.r, g.k, g.r, rows});
                    } else {
                        rep_nn.entries.push_back({t2_l, g.xrs, fresh, c_l, g.k, g.r, g.ntot, g.r, rows});
                    }
                    edges_[key] = DeviceMatrix{fresh, rows, cols};
                    fresh_made.push_back(AdoptedBlock{false, key, &edges_[key]});
                    note_write(&edges_[key], static_cast<size_t>(rows) * cols * D);
                    replaced.insert(key);
                };
                replace_view(true);
                replace_view(false);
            }
            for (size_t ni = 0; ni < hop.size(); ++ni) {
                if (ni == li) continue;
                const int64_t nm = hop[ni];
                const int c_n = static_cast<int>(counts[ni]);
                if (c_n == 0 || !deferred_xnn_boxes_are_one_hop(level_.dimension, lm, nm)) continue;
                if (interest(nm) && nm < lm) continue;  // from the lower endpoint of interest
                if (current_size(state_of(nm)) != c_n) {
                    throw std::runtime_error("LevelEliminator: generator " + std::to_string(g.morton) + " sees box " +
                                             std::to_string(nm) + " at another size: " + std::to_string(c_n) +
                                             " against " + std::to_string(current_size(state_of(nm))) +
                                             " here (next to box " + std::to_string(lm) + ")");
                }
                const uint64_t k_ln = dkey(lm, nm), k_nl = dkey(nm, lm);
                if (replaced.count(k_ln) || replaced.count(k_nl)) {
                    throw std::runtime_error("LevelEliminator: a replaced block is also updated");
                }
                const S* t2_n = g.temp2 + g.offset[ni];
                const S* xrn_n = g.xrn + static_cast<size_t>(g.offset[ni]) * r;
                accumulate(k_ln, ensure_dblock(lm, nm), c_l, c_n, t2_l, g.ntot, xrn_n, g.r, g.r);
                accumulate(k_nl, ensure_dblock(nm, lm), c_n, c_l, t2_n, g.ntot, xrn_l, g.r, g.r);
            }
        }
    };

    // ---- the generators' temp2^T first, in an image of their own: their
    // transposes into the staging and the slots packed group by group (the
    // second keep_fill_source site).  None of it touches the adopted blocks,
    // so it runs once for the pass, before the chunks below.
    // The staging, and the grouped generators' temp2^T addresses into it,
    // before the transposes are baked into the image.  No room for it means
    // the old way: a block each, packed in one pass at the end.
    char* gen_staging = nullptr;
    if (!gen_groups.empty()) {
        gen_staging = heap_.try_alloc(std::max<size_t>(gen_staging_max, 1), ResidentKind::Codec);
        for (GenGroup& gg : gen_groups) {
            for (size_t j = 0; j < gg.gens.size(); ++j) {
                Gen& g = gens[gg.gens[j]];
                char* dst = gen_staging != nullptr ? gen_staging + gg.off[j] : nullptr;
                if (dst == nullptr) {
                    const size_t t_bytes = align_up(static_cast<size_t>(g.ntot) * g.r * D);
                    dst = heap_.alloc_resident(std::max<size_t>(t_bytes, 1), ResidentKind::FillRaw);
                    auto& es = eliminator_stats();
                    es.persist_bytes += static_cast<double>(t_bytes);
                    es.persist_gen_bytes += static_cast<double>(t_bytes);
                    if (early_free_sources_) {
                        g.st->persist_t2t = dst;
                        g.st->persist_t2t_bytes = t_bytes;
                    } else {
                        level_allocs_.push_back(dst);
                    }
                }
                g.st->temp2t = reinterpret_cast<const S*>(dst);
                transposes[gg.first_tr + j].dst = reinterpret_cast<S*>(dst);
            }
        }
        if (gen_staging == nullptr) gen_groups.clear();  // one pass at the end, as before
    }
    const size_t off_transposes = meta_.append(transposes);
    {
        char* md = meta_.upload(meta_device_, stream);
        stats.bytes_up += static_cast<double>(meta_.size());
        // A generator's fill-source slots, for the packer.
        auto fill_item = [&](size_t gi, std::vector<int64_t>& counts, std::vector<FillSrcRef>& items) {
            const Gen& g = gens[gi];
            if (g.st == nullptr || g.st->temp2t == nullptr || g.r == 0 || g.ntot == 0) return;
            // (an owner level's generator may be one of this rank's own boxes,
            // landed FULL: the box found when the generator was set up, not the
            // generator table, which has only the remote ones)
            const Box* gb = g.box;
            if (gb == nullptr) return;
            // its rows per neighbor as the generator laid them out (its
            // header's counts): a rank that holds it as an assisting box
            // (unstructured: the lazy requesters) need not know every
            // neighbor of it, and on a full grid the two views agree
            counts.assign(gb->deferred_xnn_neighbor_point_counts.begin(),
                          gb->deferred_xnn_neighbor_point_counts.end());
            items.push_back(FillSrcRef{g.st->temp2t, g.r, &counts, g.st});
        };
        if (gen_groups.empty()) {
            launch_transpose(reinterpret_cast<const TransposeItem*>(md + off_transposes),
                             static_cast<int>(transposes.size()), max_tr_m, max_tr_n, stream);
            // ---- the generator fill sources, compressed a slot at a time now
            // that their transposes have written temp2^T (the second
            // keep_fill_source site)
            if (compress_factors_ >= 1) {
                std::vector<std::vector<int64_t>> counts(gens.size());
                std::vector<FillSrcRef> items;
                items.reserve(gens.size());
                for (size_t gi = 0; gi < gens.size(); ++gi) fill_item(gi, counts[gi], items);
                compress_fill_source_set(items, stream);
            }
        } else {
            // ---- group by group: the group's transposes into the staging, its
            // slots packed (which waits for the staging to be read), the staging
            // reused by the next.  A generator whose pack found no room keeps a
            // raw block of its own, copied out before the staging moves on.
            for (const GenGroup& gg : gen_groups) {
                launch_transpose(reinterpret_cast<const TransposeItem*>(md + off_transposes) + gg.first_tr,
                                 static_cast<int>(gg.gens.size()), max_tr_m, max_tr_n, stream);
                std::vector<std::vector<int64_t>> counts(gg.gens.size());
                std::vector<FillSrcRef> items;
                items.reserve(gg.gens.size());
                for (size_t j = 0; j < gg.gens.size(); ++j) fill_item(gg.gens[j], counts[j], items);
                compress_fill_source_set(items, stream);
                for (size_t j = 0; j < gg.gens.size(); ++j) {
                    Gen& g = gens[gg.gens[j]];
                    if (g.st == nullptr) continue;
                    if (packed_of(*g.st) != nullptr) {
                        g.st->temp2t = nullptr;  // the packed form is the fill source now
                        continue;
                    }
                    const size_t t_bytes = align_up(static_cast<size_t>(g.ntot) * g.r * D);
                    char* own = heap_.alloc_resident(std::max<size_t>(t_bytes, 1), ResidentKind::FillRaw);
                    check_cuda(cudaMemcpyAsync(own, gen_staging + gg.off[j], t_bytes, cudaMemcpyDeviceToDevice, stream),
                               "generator temp2^T kept raw");
                    auto& es = eliminator_stats();
                    es.persist_bytes += static_cast<double>(t_bytes);
                    es.persist_gen_bytes += static_cast<double>(t_bytes);
                    if (early_free_sources_) {
                        g.st->persist_t2t = own;
                        g.st->persist_t2t_bytes = t_bytes;
                    } else {
                        level_allocs_.push_back(own);
                    }
                    g.st->temp2t = reinterpret_cast<const S*>(own);
                }
                check_cuda(cudaStreamSynchronize(stream), "generator staging");  // before the next group overwrites it
            }
            heap_.free(gen_staging);
        }
        check_cuda(cudaStreamSynchronize(stream), "generator transposes");  // (the image is reused below)
    }
    // ---- then the updates, in chunks of generators whose adopted read set
    // fits kNearGroupBudget raw.  A chunk's blocks are put back, written and
    // packed again before the next chunk's, so a pass's union is never
    // resident at once: it was 2.2 GB of the 2.8 GB level-4 peak at 192^3,
    // where a wave's own gather needs 0.6.  A block two chunks read is
    // packed twice, which the finer step of these packs allows.
    size_t g0 = 0;
    while (g0 < gens.size()) {
        std::vector<AdoptedBlock> adopted_read;   // what this chunk puts back; packed again at its end
        size_t g1 = g0;
#ifdef H2_HAVE_NVCOMP
        if (!(packed_schur_.empty() && packed_edges_.empty())) {
            std::unordered_set<const void*> seen;
            size_t bytes = 0;
            while (g1 < gens.size()) {
                std::vector<AdoptedBlock> cand;
                size_t added = 0;   // what the generator puts back, and what it makes
                want_of_gen(gens[g1], seen, cand, added);
                for (const AdoptedBlock& w : cand) added += static_cast<size_t>(w.m->rows) * w.m->cols * D;
                if (g1 > g0 && bytes + added > near_group_budget()) break;
                for (const AdoptedBlock& w : cand) {
                    seen.insert(w.m);
                    adopted_read.push_back(w);
                }
                bytes += added;
                ++g1;
            }
            materialize_adopted(adopted_read, stream, MatGenerator);
        } else {
            g1 = gens.size();
        }
#else
        g1 = gens.size();
#endif
        ++wave_stamp_;
        meta_.clear();
        evals.clear();
        gathers.clear();
        products.clear();
        accums.clear();
        accum_of.clear();
        rep_tt.entries.clear();
        rep_nn.entries.clear();
        fresh_made.clear();
        max_ev_m = max_ev_n = max_ga_m = max_ga_n = 0;
        work_bytes = 0;
        for (size_t gi = g0; gi < g1; ++gi) {
            if (nonsym_) process_gen_nonsym(gens[gi]); else process_gen(gens[gi]);
        }

        // Every partial product added straight into its target, sub-batch j
        // holding the j-th partial of every target, so no launch writes a block
        // twice.  This is what the local owner pass does; staging them all first
        // and summing afterwards cost 4.62 GB at the 192^3 Color leaf, 1.4x the
        // targets it summed into, because most targets take a single partial.
        size_t depth = 0;
        for (const Accum& a : accums) depth = std::max(depth, a.parts.size());
        {
            auto& st = eliminator_stats();
            st.gw_bytes += static_cast<double>(work_bytes);  // what staging would have cost
            st.gw_products += static_cast<int64_t>(products.size());
            st.gw_targets += static_cast<int64_t>(accums.size());
            st.gw_max_parts = std::max(st.gw_max_parts, static_cast<int64_t>(depth));
            for (const Accum& a : accums) {
                st.gw_target_bytes += static_cast<double>(static_cast<size_t>(a.rows) * a.cols * D);
            }
        }
        std::vector<VBatch<S>> prods(depth);
        for (const Accum& acc : accums) {
            for (size_t j = 0; j < acc.parts.size(); ++j) {
                const Product& p = products[acc.parts[j]];
                prods[j].entries.push_back({p.a, p.b, acc.target, p.m, p.n, p.k, p.lda, p.ldb, acc.ld});
            }
        }
        const size_t off_evals = meta_.append(evals);
        const size_t off_gathers = meta_.append(gathers);
        for (VBatch<S>& b : prods) b.stage(meta_);
        rep_tt.stage(meta_);
        rep_nn.stage(meta_);
        char* md = meta_.upload(meta_device_, stream);
        stats.bytes_up += static_cast<double>(meta_.size());
        {
            EntryOriginScope origin(EntryOrigin::remote);  // (collective entries' profile)
            eval_blocks(evals, meta_, md, off_evals, max_ev_m, max_ev_n, stream);
        }
        launch_gather(reinterpret_cast<const GatherItem*>(md + off_gathers), static_cast<int>(gathers.size()), max_ga_m,
                      max_ga_n, md, stream);
        for (const VBatch<S>& b : prods) {
            check_targets(b, "generator products");
            b.gemm(md, MagmaNoTrans, nonsym_ ? MagmaNoTrans : MagmaTrans, 1.0, 1.0, queue);
        }
        rep_tt.gemm(md, MagmaTrans, MagmaTrans, 1.0, 1.0, queue);
        rep_nn.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 1.0, queue);
        check_targets(rep_tt, "generator rep_tt");
        check_targets(rep_nn, "generator rep_nn");
        check_cuda(cudaStreamSynchronize(stream), "generator chunk");  // the image is reused by the next chunk
    #ifdef H2_HAVE_NVCOMP
        // The adopted blocks this pass put back go packed again right away, the
        // way the owner pass's groups do.  Left raw, they stayed resident until
        // some later reader happened to retire them: 1.76-1.87 GB of "adopted
        // blocks" at the 192^3 level-4 peak, against a 0.6 GB per-wave read
        // union (GPU_CA_PLAN.md, piece 4).  The map slots are what the pass
        // wrote into, whether it accumulated into the adopted block or replaced
        // it with a fresh one, so the slot's current block is what goes back.
        if (kHoldNearBlocksPacked && compress_factors_ >= 1 && !(adopted_read.empty() && fresh_made.empty()) &&
            under_pressure()) {
            struct PackItem { DeviceMatrix* block; uint64_t key; int rows, cols; };
            std::vector<PackItem> items;
            items.reserve(adopted_read.size() + fresh_made.size());
            std::unordered_set<const DeviceMatrix*> slots;   // (a replaced block is in both lists)
            for (const AdoptedBlock& w : adopted_read) {
                if (w.m->ptr == nullptr || w.m->rows <= 0 || w.m->cols <= 0 || !slots.insert(w.m).second) continue;
                items.push_back(PackItem{w.m, w.schur ? ~w.key : w.key, w.m->rows, w.m->cols});
                ++stats.np_gen_blocks;
                stats.np_gen_raw += static_cast<double>(static_cast<size_t>(w.m->rows) * w.m->cols * D);
            }
            for (const AdoptedBlock& w : fresh_made) {
                if (w.m->ptr == nullptr || w.m->rows <= 0 || w.m->cols <= 0 || !slots.insert(w.m).second) continue;
                items.push_back(PackItem{w.m, w.schur ? ~w.key : w.key, w.m->rows, w.m->cols});
                ++stats.np_gen_fresh_blocks;
                stats.np_gen_fresh_raw += static_cast<double>(static_cast<size_t>(w.m->rows) * w.m->cols * D);
            }
            size_t at = 0;
            while (at < items.size()) {
                size_t bytes = 0, end = at;
                while (end < items.size() && (end == at || bytes < kNearGroupBudget)) {
                    bytes += static_cast<size_t>(items[end].rows) * items[end].cols * D;
                    ++end;
                }
                pack_near_group(items, at, end, stream, 0.25);   // (finer: these are packed again by the owner pass)
                at = end;
            }
        }
    #endif
        // the replaced blocks this chunk's gathers read (later launches are
        // ordered after them)
        for (char* p : retired) heap_.free(p);
        retired.clear();
        g0 = g1;
    }
    tmp.release();

    if (heap_.peak() > stats.heap_peak) {
        stats.heap_peak = heap_.peak();
        stats.heap_peak_resident = heap_.resident_at_peak();
        stats.heap_peak_transient = heap_.transient_at_peak();
    }
}

// Generators of the last wave's boxes with a neighbor on another rank, as
// the host's emit_lazy_generators packages them (X_NR here is the original
// coupling, downloaded for these boxes; the host X_NR holds temp2).
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::emit_generators(PendingFactorUpdates<DataType>& pending) {
    if (last_wave_.empty()) return;
    if (copier_) free_copied_blocks(copier_->wait_all());  // the host copies of the wave
    for (int64_t m : last_wave_) {
        Box* box = level_.find_local_box(m);
        if (box == nullptr) continue;
        bool has_remote = false;
        for (int64_t nb : box->one_hop) has_remote = has_remote || !is_local(nb);
        if (!has_remote) continue;
        auto generator = std::make_shared<GeneratorPayload<DataType>>();
        generator->wave = static_cast<int32_t>(last_wave_index_);
        generator->one_hop = box->one_hop;
        const bool has_fill = box->X_NR.is_allocated() && box->X_NR.cols > 0;
        if (has_fill) {
            const int64_t r = box->X_NR.cols, rows = box->X_NR.rows;
            const int64_t k = static_cast<int64_t>(box->skeleton_indices.size());
            if (box->deferred_xnn_neighbor_point_counts.size() != box->one_hop.size() ||
                static_cast<int64_t>(box->lazy_original_x_nr.size()) != rows * r ||
                box->X_RR_full.rows != r || box->X_RS_entry.rows != r || box->X_RS_entry.cols != k) {
                throw std::runtime_error("LevelEliminator::emit_generators: incomplete state for source " +
                                         std::to_string(m));
            }
            generator->r = r;
            generator->total_rows = rows;
            generator->neighbor_point_counts = box->deferred_xnn_neighbor_point_counts;
            generator->temp2.assign(box->X_NR.data.begin(), box->X_NR.data.begin() + rows * r);
            generator->original_x_nr = std::move(box->lazy_original_x_nr);
            box->lazy_original_x_nr.clear();
            generator->x_rr_full.assign(box->X_RR_full.data.begin(), box->X_RR_full.data.begin() + r * r);
            generator->k = k;
            generator->skeleton_indices = box->skeleton_indices;
            generator->x_rs.assign(box->X_RS_entry.data.begin(), box->X_RS_entry.data.begin() + r * k);
        }
        pending.generators[m] = std::move(generator);
    }
    last_wave_.clear();
}

// ---------------------------------------------------------------------------
// The solve's factors of a wave (T, X_SR, X_NR, the LU of X_RR and its
// pivots, as the host keeps them) copied into a device block of their own,
// while the heap stays under kSolveKeepFraction of its capacity;
// past that, this rank drops the level's copies.  Launched before the
// wave's result and source blocks can be released.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// A compressed fill source put back whole, slot by slot in their original
// order, which is how temp2^T was laid out.  The transition reads a source's
// rows from all over it, so it takes the whole thing rather than slots; its
// block lives for the chunk, like the sources restored from the host.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::collect_fill_source(
    const BoxState& st, S* dst, std::vector<const void*>& h_src, std::vector<void*>& h_dst) {
#ifdef H2_HAVE_NVCOMP
    if (packed_of(st) == nullptr || st.t2t_slots.empty()) return false;
    const size_t r = static_cast<size_t>(st.r);
    const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
    size_t row = 0;
    for (const auto& sl : st.t2t_slots) {
        const char* csrc = packed_of(st) + sl.off;
        if (st.t2t_packed != nullptr) ++eliminator_stats().fs_read_dev[2]; else ++eliminator_stats().fs_read_host[2];
        S* out = dst + row * r;
        for (uint32_t c = 0; c < sl.nchunks; ++c) {
            h_src.push_back(csrc);
            h_dst.push_back(out + static_cast<size_t>(c) * chunk_elems);
            csrc += st.t2t_strides[sl.chunk0 + c];
        }
        row += static_cast<size_t>(sl.rows);
    }
    eliminator_stats().fs_tr_bytes += static_cast<double>(static_cast<size_t>(st.ntot) * r * sizeof(S));
    ++eliminator_stats().fs_tr_sources;
    return true;
#else
    (void)st; (void)dst; (void)h_src; (void)h_dst;
    return false;
#endif
}

// ---------------------------------------------------------------------------
// One batch over every chunk collected for a transition chunk's fill sources.
// Per source this used to be its own plan, pointer allocation and two stream
// synchronizations, which the transition paid for every source of every chunk
// - 38% of the factorization.  The same mistake the sketch's group loop made.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::run_fill_batch(
    const std::vector<const void*>& h_src, const std::vector<void*>& h_dst, cudaStream_t stream) {
#ifdef H2_HAVE_NVCOMP
    const size_t n = h_src.size();
    if (n == 0) return;
    const size_t span = align_up(n * sizeof(void*));
    char* ptrs = heap_.alloc(span * 2);
    auto** d_src = reinterpret_cast<const void**>(ptrs);
    auto** d_dst = reinterpret_cast<void**>(ptrs + span);
    check_cuda(cudaMemcpyAsync(d_src, h_src.data(), n * sizeof(void*), cudaMemcpyHostToDevice, stream), "fill sources");
    check_cuda(cudaMemcpyAsync(d_dst, h_dst.data(), n * sizeof(void*), cudaMemcpyHostToDevice, stream), "fill sources");
    check_cuda(cudaStreamSynchronize(stream), "fill sources");
    bitcompHandle_t plan{};
    if (bitcompCreateBatchPlanFromCompressedData(&plan, d_src, n, nullptr, stream) != BITCOMP_SUCCESS) {
        heap_.free(ptrs);
        throw std::runtime_error("LevelEliminator: could not plan the fill sources' decompression");
    }
    bitcompSetStream(plan, stream);
    const bitcompResult_t rc = bitcompBatchUncompress(plan, d_src, d_dst);
    bitcompDestroyPlan(plan);
    if (rc != BITCOMP_SUCCESS) {
        heap_.free(ptrs);
        throw std::runtime_error("LevelEliminator: decompressing the fill sources failed");
    }
    check_cuda(cudaStreamSynchronize(stream), "fill sources");
    heap_.free(ptrs);
#else
    (void)h_src; (void)h_dst; (void)stream;
#endif
}

#ifdef H2_HAVE_NVCOMP
// ---------------------------------------------------------------------------
// The adopted blocks a wave is about to read, put back in one batch.  Every
// block is read once and replaced by that wave, so its packed entry goes with
// it and the key is free for the fresh block.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::materialize_adopted(
    const std::vector<AdoptedBlock>& want, cudaStream_t stream, MatSite site, const RawTake* take) {
    if (want.empty()) return;
    GatherPlan g;
    g.want = want;
    plan_gather(g, stream, nullptr);
    materialize_from_gather(g, stream, site, take);
}

// The gather half.  A block written by build_parent is one contiguous run
// (in the parent level's host region), which bitcomp reads in place.  A
// near block is a run of chunk homes - each a piece in B_compact plus an
// optional tail in B_newtail, both host memory - so its bytes are gathered
// into one device scratch first, by the span copy on `stream`.  On the near
// copy stream (a prefetch) the copy waits for the compute stream to pass
// its last decompression, whose scratch the heap may have handed out again,
// and `done` is recorded after it.  Nothing is allocated for the raw blocks
// here, so a prefetched plan costs only the packed bytes.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::plan_gather(GatherPlan& g, cudaStream_t stream,
                                                                     cudaEvent_t done) {
    auto& stats = eliminator_stats();
    const auto t0 = std::chrono::steady_clock::now();
    auto packed_of = [&](const AdoptedBlock& w) -> const PackedBlock* {
        if (w.schur) {
            auto it = packed_schur_.find(static_cast<int64_t>(w.key));
            return it == packed_schur_.end() ? nullptr : &it->second;
        }
        auto it = packed_edges_.find(w.key);
        return it == packed_edges_.end() ? nullptr : &it->second;
    };
    const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
    size_t gather_bytes = 0;
    uint32_t need = 0;   // the last scatter that wrote any of them
    for (const AdoptedBlock& w : g.want) {
        if (w.m->ptr != nullptr) continue;
        const PackedBlock* pb = packed_of(w);
        if (pb == nullptr || pb->packed != nullptr || pb->nchunks == 0) continue;
        for (uint32_t c = 0; c < pb->nchunks; ++c) gather_bytes += chunk_homes_[pb->home0 + c].total;
        need = std::max(need, pb->written_seq);
    }
    g.gathered = gather_bytes > 0 ? heap_.alloc(align_up(gather_bytes), ResidentKind::Gather) : nullptr;
    std::vector<CopySpan> gspans;
    int64_t gmax = 0;
    size_t gat = 0;
    for (uint32_t bi = 0; bi < g.want.size(); ++bi) {
        const AdoptedBlock& w = g.want[bi];
        if (w.m->ptr != nullptr) continue;
        const PackedBlock* pb = packed_of(w);
        if (pb == nullptr || (pb->packed == nullptr && pb->nchunks == 0)) continue;
        if (pb->packed != nullptr) {          // build_parent's contiguous run: read in place
            const char* at = pb->packed;
            for (uint32_t c = 0; c < pb->nchunks; ++c) {
                g.src.push_back(at);
                g.dst_block.push_back(bi);
                g.dst_chunk.push_back(c);
                at += packed_strides_[pb->chunk0 + c];
            }
            continue;
        }
        if (g.gathered == nullptr) continue;
        for (uint32_t c = 0; c < pb->nchunks; ++c) {
            const ChunkHome& h = chunk_homes_[pb->home0 + c];
            char* dst = g.gathered + gat;
            const uint32_t in_home = std::min(h.here, h.total);
            gspans.push_back(CopySpan{h.at, dst, static_cast<int64_t>(in_home / 4)});
            gmax = std::max(gmax, static_cast<int64_t>(in_home / 4));
            if (h.total > in_home) {
                const uint32_t rest = h.total - in_home;
                gspans.push_back(CopySpan{tail_at(h.tail), dst + in_home, static_cast<int64_t>(rest / 4)});
                gmax = std::max(gmax, static_cast<int64_t>(rest / 4));
            }
            g.src.push_back(dst);
            g.dst_block.push_back(bi);
            g.dst_chunk.push_back(c);
            gat += h.total;
        }
    }
    (void)chunk_elems;
    if (!gspans.empty()) {
        const bool prefetch = done != nullptr;
        DeviceBuffer& meta = prefetch ? near_prefetch_meta_ : fill_meta_device_;
        const size_t list_bytes = gspans.size() * sizeof(CopySpan);
        if (!prefetch && list_bytes > fill_meta_device_.capacity()) {
            check_cuda(cudaStreamSynchronize(stream), "near blocks gathered");
        }
        auto* d_items = static_cast<CopySpan*>(meta.reserve(list_bytes));
        if (prefetch && near_reuse_ != nullptr) {
            check_cuda(cudaStreamWaitEvent(stream, near_reuse_, 0), "near blocks gathered");
        }
        wait_near_scatter(stream, need);
        check_cuda(cudaMemcpyAsync(d_items, gspans.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                   "near blocks gathered");
        launch_copy_spans(d_items, static_cast<int>(gspans.size()), gmax, stream);
        if (prefetch) {
            check_cuda(cudaEventRecord(done, stream), "near blocks gathered");
            g.done = done;
        } else {
            check_cuda(cudaStreamSynchronize(stream), "near blocks gathered");
        }
    }
    stats.mat_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// The decompression half: the raw blocks allocated, the batch run from the
// gathered pieces and the in-place runs, the scratch returned.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::materialize_from_gather(GatherPlan& g, cudaStream_t stream,
                                                                                 MatSite site, const RawTake* take) {
    auto& stats = eliminator_stats();
    const auto t0 = std::chrono::steady_clock::now();
    auto packed_of = [&](const AdoptedBlock& w) -> const PackedBlock* {
        if (w.schur) {
            auto it = packed_schur_.find(static_cast<int64_t>(w.key));
            return it == packed_schur_.end() ? nullptr : &it->second;
        }
        auto it = packed_edges_.find(w.key);
        return it == packed_edges_.end() ? nullptr : &it->second;
    };
    const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
    // 1. every block allocated first: alloc_resident can run the reclaimer,
    //    so no block address may be held across an allocation
    for (const AdoptedBlock& w : g.want) {
        if (w.m->ptr != nullptr) continue;
        const PackedBlock* pb = packed_of(w);
        // `packed == nullptr` used to mean "no bytes".  In the near-block
        // store it means the opposite - the bytes are in chunk homes - so
        // testing it alone skipped every near block, left its ptr null, and
        // handed null + offset to whatever read it next (compute-sanitizer:
        // an invalid 8-byte read at 0x21c0 in sum_add_kernel, from
        // apply_remote_generators).  A block has bytes if either form is set.
        if (pb == nullptr || (pb->packed == nullptr && pb->nchunks == 0)) continue;
        if (pb->rows != w.m->rows || pb->cols != w.m->cols) {
            throw std::runtime_error("LevelEliminator::materialize_adopted: shape of the packed block changed");
        }
        const size_t raw = static_cast<size_t>(pb->rows) * pb->cols * sizeof(S);
        char* taken = take != nullptr ? (*take)(w, raw) : nullptr;
        w.m->ptr = taken != nullptr ? reinterpret_cast<S*>(taken) : heap_.alloc_resident<S>(raw, ResidentKind::Adopted);
        const double bytes = static_cast<double>(raw);
        stats.mat_bytes += bytes;
        if (site < MatSiteCount) stats.mat_site_bytes[site] += bytes;
        ++stats.mat_blocks;
    }
    // 2. their chunks' destinations
    std::vector<void*> h_dst(g.src.size());
    for (size_t i = 0; i < g.src.size(); ++i) {
        const AdoptedBlock& w = g.want[g.dst_block[i]];
        if (w.m->ptr == nullptr) throw std::runtime_error("LevelEliminator: a gathered block has no raw memory");
        h_dst[i] = reinterpret_cast<S*>(w.m->ptr) + static_cast<size_t>(g.dst_chunk[i]) * chunk_elems;
    }
    // 3. one batch for all of them, after the gather if it ran elsewhere
    if (g.done != nullptr) check_cuda(cudaStreamWaitEvent(stream, g.done, 0), "prefetched gather");
    if (!g.src.empty()) {
        run_fill_batch(g.src, h_dst, stream);
        ++stats.mat_batches;
        if (site < MatSiteCount) ++stats.mat_site_calls[site];
        if (near_reuse_ == nullptr) near_reuse_ = make_event();
        check_cuda(cudaEventRecord(near_reuse_, stream), "near blocks gathered");
    }
    heap_.free(g.gathered);
    g.gathered = nullptr;
    for (const AdoptedBlock& w : g.want) {
        PackedBlock* pb = w.schur ? (packed_schur_.count(static_cast<int64_t>(w.key))
                                         ? &packed_schur_[static_cast<int64_t>(w.key)] : nullptr)
                                  : (packed_edges_.count(w.key) ? &packed_edges_[w.key] : nullptr);
        if (pb == nullptr) continue;
        if (pb->packed != nullptr) {
            drop_packed(w.schur, w.key);   // build_parent's run: erase and reclaim
        } else {
            pb->live = false;              // near block: keep its homes for the next compression
        }
    }
    stats.mat_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
#endif

// ---------------------------------------------------------------------------
// The wave's local fill sources, a sub-group at a time: a group's temp2^T is
// given its block, transposed into it, compressed and the block handed back
// before the next group starts.  So a wave holds one group of raw temp2^T, not
// all of it - 2.26 GB a wave at the Color leaf against about 0.1 GB a group.
//
// The group is the same size the compressor's staging allows
// (kFillStagingBudget), since the two go together: a group is transposed,
// staged, packed and freed as one.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::compress_fill_sources(
    std::vector<WaveBox>& boxes, const std::vector<int>& active, cudaStream_t stream) {
    const size_t D = sizeof(S);
    std::vector<int> group;          // boxes of the current group
    std::vector<TransposeItem> tr;   // their temp2^T transposes
    int max_m = 0, max_n = 0;
    size_t group_raw = 0;

    auto run = [&]() {
        if (group.empty()) return;
        // the group's blocks, and the transposes that fill them
        tr.clear();
        max_m = 0;
        max_n = 0;
        std::vector<FillSrcRef> items;
        items.reserve(group.size());
        for (int i : group) {
            WaveBox& wb = boxes[static_cast<size_t>(i)];
            const size_t r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
            const size_t bytes = std::max<size_t>(nt * r * D, 1);
            // A box without its temp2^T has no fill source and quietly changes
            // the answer, so this is not a place to degrade: alloc_resident
            // throws if the group does not fit, and the group is bounded.
            char* blk = heap_.alloc_resident(bytes, ResidentKind::FillRaw);
            auto& es = eliminator_stats();
            es.persist_bytes += static_cast<double>(bytes);
            es.persist_peak = std::max(es.persist_peak, es.persist_bytes);
            wb.state->persist_t2t = blk;
            wb.state->persist_t2t_bytes = bytes;
            wb.persist_temp2t = reinterpret_cast<S*>(blk);
            tr.push_back(TransposeItem{wb.temp2, static_cast<int>(nt), wb.persist_temp2t,
                                       static_cast<int>(r), static_cast<int>(nt), static_cast<int>(r)});
            max_m = std::max(max_m, static_cast<int>(nt));
            max_n = std::max(max_n, static_cast<int>(r));
            items.push_back(FillSrcRef{wb.persist_temp2t, wb.r, &wb.counts, wb.state});
        }
        if (!tr.empty()) {
            const size_t list_bytes = tr.size() * sizeof(TransposeItem);
            if (list_bytes > t2t_meta_device_.capacity()) {
                check_cuda(cudaStreamSynchronize(stream), "fill transposes");  // the old list may still be read
                t2t_meta_device_.reserve(list_bytes);
            }
            auto* d_tr = static_cast<TransposeItem*>(t2t_meta_device_.reserve(list_bytes));
            check_cuda(cudaMemcpyAsync(d_tr, tr.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                       "fill transposes");
            launch_transpose(d_tr, static_cast<int>(tr.size()), max_m, max_n, stream);
        }
        // compressed and the raw blocks handed back (compress_fill_source_set
        // frees them once its gather is on the stream); without compression
        // the raw blocks are what the readers use, so they stay
        if (compress_factors_ >= 1) compress_fill_source_set(items, stream);
        // X_NR^T, the owner pass's other operand, cut into the same slots.
        // Only this wave's owner pass reads it, so the transposed staging is
        // transient and goes back here, and the store goes at the end of the
        // wave (release_owner_sources).
        if (compress_factors_ >= 1) {
            std::vector<FillSrcRef> xitems;
            std::vector<char*> staging;
            std::vector<TransposeItemT<S>> xtr;
            int xm = 0, xn = 0;
            bool room = true;
            for (int i : group) {
                WaveBox& wb = boxes[static_cast<size_t>(i)];
                const size_t r = static_cast<size_t>(wb.r), nt = static_cast<size_t>(wb.ntot);
                char* blk = heap_.try_alloc(std::max<size_t>(nt * r * D, 1), ResidentKind::Codec);
                if (blk == nullptr) { room = false; break; }
                staging.push_back(blk);
                xtr.push_back(TransposeItemT<S>{wb.xnr, static_cast<int>(nt), reinterpret_cast<S*>(blk),
                                                static_cast<int>(r), static_cast<int>(nt), static_cast<int>(r)});
                xm = std::max(xm, static_cast<int>(nt));
                xn = std::max(xn, static_cast<int>(r));
                xitems.push_back(FillSrcRef{reinterpret_cast<S*>(blk), wb.r, &wb.counts, wb.state});
            }
            if (room && !xtr.empty()) {
                const size_t list_bytes = xtr.size() * sizeof(TransposeItemT<S>);
                if (list_bytes > t2t_meta_device_.capacity()) {
                    check_cuda(cudaStreamSynchronize(stream), "owner source transposes");
                    t2t_meta_device_.reserve(list_bytes);
                }
                auto* d_tr = static_cast<TransposeItemT<S>*>(t2t_meta_device_.reserve(list_bytes));
                check_cuda(cudaMemcpyAsync(d_tr, xtr.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                           "owner source transposes");
                launch_transpose(d_tr, static_cast<int>(xtr.size()), xm, xn, stream);
                compress_fill_source_set(xitems, stream, true);  // (it syncs before returning)
            } else {
                ++eliminator_stats().xt_short;
            }
            for (char* b : staging) heap_.free(b);
        }
        group.clear();
        group_raw = 0;
    };

    for (int i : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(i)];
        if (wb.r == 0 || wb.ntot == 0 || wb.counts.empty()) continue;
        const size_t bytes = static_cast<size_t>(wb.ntot) * static_cast<size_t>(wb.r) * D;
        if (!group.empty() && group_raw + bytes > kFillStagingBudget) run();
        group.push_back(i);
        group_raw += bytes;
    }
    run();
}

// ---------------------------------------------------------------------------
// The fill sources of `items`, compressed a slot at a time, in groups whose
// raw bytes fit the staging the heap can spare.  The slot counts partition
// temp2^T by one-hop neighbor, which is the granularity the readers take.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::compress_fill_source_set(
    const std::vector<FillSrcRef>& items, cudaStream_t stream, bool into_xnrt) {
#ifdef H2_HAVE_NVCOMP
    const size_t D = sizeof(S);
    const bool check = env::check("fillsrc");
    auto& stats = eliminator_stats();

    decide_fill_precision();
    std::vector<const S*> srcs;
    std::vector<size_t> elems;
    std::vector<int> slot_rows;  // per block, its uncompressed rows
    std::vector<int> slot_box;   // ... and its box, whose maximum sets its step
    // temp2^T or X_NR^T: the same packing into different fields of the box.
    struct Dest {
        char*& packed;
        size_t& bytes;
        std::vector<typename BoxState::T2tSlot>& slots;
        std::vector<uint32_t>& strides;
    };
    auto dest_of = [into_xnrt](BoxState& st) -> Dest {
        if (into_xnrt) {
            return Dest{st.xnrt_packed, st.xnrt_packed_bytes, st.xnrt_slots, st.xnrt_strides};
        }
        return Dest{st.t2t_packed, st.t2t_packed_bytes, st.t2t_slots, st.t2t_strides};
    };
    const ResidentKind dest_kind = into_xnrt ? ResidentKind::OwnerSrc : ResidentKind::FillPacked;
    struct GroupItem { BoxState* st; size_t first; size_t last; };
    std::vector<GroupItem> group_of;
    size_t group_raw = 0;

    // Gather one group's compressed slots out of the compressor's staging into
    // a packed block per box, recording where each slot's bytes start and
    // which chunks they are.  A box whose packed block does not fit keeps its
    // raw temp2^T, which is correct and only costs memory.
    auto pack = [&](const std::vector<CompressedXnr>& comp) {
        const auto& cks = fill_compressor_.chunks();
        const char* staging = fill_compressor_.staging();
        std::vector<CopySpan> spans;
        std::vector<BoxState*> raw_to_free;
        std::vector<BoxState*> packed_now;  // (temp2^T only: X_NR^T stays on the device for its wave)
        int64_t max_words = 0;
        for (const GroupItem& gi : group_of) {
            BoxState& st = *gi.st;
            Dest dst = dest_of(st);
            size_t total = 0;
            for (size_t i = gi.first; i < gi.last; ++i) total += comp[i].bytes;
            if (total == 0) continue;
            char* packed = heap_.try_alloc_resident(total, dest_kind);
            if (packed == nullptr) {  // no room: this box stays raw
                ++stats.fs_packed_short;
                warn_replica_fallback("a fill source stayed raw (no room for its packed block)");
                continue;
            }
            if (dst.packed != nullptr) {  // (a box compressed again drops the old one)
                stats.persist_bytes -= static_cast<double>(dst.bytes);
                heap_.free(dst.packed);
                dst.slots.clear();
                dst.strides.clear();
                if (!into_xnrt) st.t2t_rows.clear();
            }
            dst.packed = packed;
            dst.bytes = total;
            if (!into_xnrt) packed_now.push_back(&st);
            stats.persist_bytes += static_cast<double>(total);
            stats.persist_peak = std::max(stats.persist_peak, stats.persist_bytes);
            uint64_t off = 0;
            int row = 0;
            for (size_t i = gi.first; i < gi.last; ++i) {
                typename BoxState::T2tSlot sl;
                sl.off = off;
                sl.chunk0 = static_cast<uint32_t>(dst.strides.size());
                sl.nchunks = static_cast<uint32_t>(comp[i].nchunks);
                sl.rows = slot_rows[i];
                for (size_t c = 0; c < comp[i].nchunks; ++c) {
                    const auto& ck = cks[comp[i].chunk0 + c];
                    spans.push_back(CopySpan{staging + ck.offset, packed + off,
                                             static_cast<int64_t>(ck.stride / 4)});
                    max_words = std::max(max_words, static_cast<int64_t>(ck.stride / 4));
                    dst.strides.push_back(static_cast<uint32_t>(ck.stride));
                    off += ck.stride;
                }
                // the row map is the partition, which both stores share
                if (!into_xnrt) st.t2t_rows.emplace(row, static_cast<int>(dst.slots.size()));
                dst.slots.push_back(sl);
                row += sl.rows;
            }
            // The packed block now holds what the raw one did - the gather
            // below is synchronized before the staging goes back - so the raw
            // block can go.  st.temp2t going null is what sends both readers
            // to the packed form.
            if (!into_xnrt) raw_to_free.push_back(&st);  // X_NR^T's staging is the caller's
        }
        if (spans.empty()) return;
        const size_t list_bytes = spans.size() * sizeof(CopySpan);
        if (list_bytes > fill_meta_device_.capacity()) {
            check_cuda(cudaStreamSynchronize(stream), "fill sources");  // the old list may still be read
            fill_meta_device_.reserve(list_bytes);
        }
        auto* d_items = static_cast<CopySpan*>(fill_meta_device_.reserve(list_bytes));
        check_cuda(cudaMemcpyAsync(d_items, spans.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                   "fill sources");
        launch_copy_spans(d_items, static_cast<int>(spans.size()), max_words, stream);
        // the gather reads the staging, which goes back to the heap next
        check_cuda(cudaStreamSynchronize(stream), "fill sources");
        // ---- the host copy (H2 host staging of the fill store): the same
        // bytes in pinned memory, which the later readers decompress from
        // over PCIe; the device copy goes once this wave's owner pass has
        // read it (evict_fill_copies, when the copies are done)
        if (!packed_now.empty() && under_pressure()) {   // (else stage_fill_copies, when the room is needed)
            PinnedArena& arena = fill_host_arena();
            cudaStream_t cs = fill_copy_stream();
            for (BoxState* st : packed_now) {
                st->t2t_host = arena.alloc(st->t2t_packed_bytes);
                check_cuda(cudaMemcpyAsync(st->t2t_host, st->t2t_packed, st->t2t_packed_bytes, cudaMemcpyDeviceToHost,
                                           cs), "fill store to the host");
            }
            stats.fs_host_bytes = std::max(stats.fs_host_bytes, static_cast<double>(arena.used()));
            PendingEvict pe;
            check_cuda(cudaEventCreateWithFlags(&pe.done, cudaEventDisableTiming), "fill store event");
            check_cuda(cudaEventRecord(pe.done, cs), "fill store event");
            pe.boxes = packed_now;
            pending_evictions_.push_back(std::move(pe));
        }
        // (the raw blocks go back below, after the check has compared against
        // them - BPACK_CHECK=fillsrc read freed memory otherwise)

        auto free_raw = [&]() {
            for (BoxState* st : raw_to_free) {
                if (st->persist_t2t == nullptr) continue;
                stats.persist_bytes -= static_cast<double>(st->persist_t2t_bytes);
                heap_.free(st->persist_t2t);
                st->persist_t2t = nullptr;
                st->persist_t2t_bytes = 0;
                st->temp2t = nullptr;  // the readers take the packed form now
            }
            raw_to_free.clear();
        };

        // BPACK_CHECK=fillsrc: read every slot back out of the packed block,
        // through the offsets and strides just recorded, and compare it with
        // what it was made from.  This is what tests the gather and the
        // offset arithmetic - a wrong stride would otherwise surface as a
        // plausible but wrong answer much later.
        if (!check) {
            free_raw();
            return;
        }
        size_t widest = 0;
        for (size_t e : elems) widest = std::max(widest, e);
        char* back = heap_.try_alloc(std::max<size_t>(widest * sizeof(S), 1));
        if (back == nullptr) {
            free_raw();
            return;
        }
        std::vector<S> got, want;
        const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
        for (const GroupItem& gi : group_of) {
            BoxState& st = *gi.st;
            Dest dst = dest_of(st);
            if (dst.packed == nullptr) continue;
            const size_t slot0 = dst.slots.size() - (gi.last - gi.first);
            for (size_t i = gi.first; i < gi.last; ++i) {
                const size_t si = slot0 + (i - gi.first);
                if (si >= dst.slots.size()) break;
                const auto& sl = dst.slots[si];
                decompress_xnr<S>(dst.packed + sl.off, reinterpret_cast<S*>(back),
                                  dst.strides.data() + sl.chunk0, sl.nchunks, stream);
                check_cuda(cudaStreamSynchronize(stream), "fillsrc check");
                got.assign(elems[i], S{});
                want.assign(elems[i], S{});
                check_cuda(cudaMemcpy(got.data(), back, elems[i] * sizeof(S), cudaMemcpyDeviceToHost),
                           "fillsrc check");
                check_cuda(cudaMemcpy(want.data(), srcs[i], elems[i] * sizeof(S), cudaMemcpyDeviceToHost),
                           "fillsrc check");
                const double* g = reinterpret_cast<const double*>(got.data());
                const double* w = reinterpret_cast<const double*>(want.data());
                const size_t nd = elems[i] * (sizeof(S) / sizeof(double));
                for (size_t j = 0; j < nd; ++j) {
                    stats.fs_check_err = std::max(stats.fs_check_err, std::fabs(g[j] - w[j]));
                }
                stats.fs_check_tol = std::max(stats.fs_check_tol, comp[i].delta * 0.5);
            }
        }
        heap_.free(back);
        free_raw();
    };

    auto flush = [&]() {
        if (!srcs.empty()) {
            std::vector<CompressedXnr> comp;
            if (fill_compressor_.compress(srcs, elems, tolerance_ * fill_frac_, comp, stream,
                                          &slot_box)) {
                size_t raw = 0, stored = 0;
                for (size_t i = 0; i < srcs.size(); ++i) {
                    raw += elems[i] * D;
                    stored += comp[i].bytes;
                }
                if (into_xnrt) {
                    stats.xt_raw_bytes += static_cast<double>(raw);
                    stats.xt_stored_bytes += static_cast<double>(stored);
                    stats.xt_slots += static_cast<int64_t>(srcs.size());
                    stats.xt_boxes += static_cast<int64_t>(group_of.size());
                } else {
                    stats.fs_raw_bytes += static_cast<double>(raw);
                    stats.fs_stored_bytes += static_cast<double>(stored);
                    stats.fs_slots += static_cast<int64_t>(srcs.size());
                    stats.fs_boxes += static_cast<int64_t>(group_of.size());
                }
                pack(comp);
            } else {
                ++stats.fs_waves_plain;
                warn_replica_fallback("a group of fill sources stayed raw (no staging)");
            }
            fill_compressor_.release();
        }
        srcs.clear();
        elems.clear();
        slot_rows.clear();
        slot_box.clear();
        group_of.clear();
        group_raw = 0;
    };

    for (const FillSrcRef& it : items) {
        if (it.t2t == nullptr || it.r == 0 || it.st == nullptr || it.counts == nullptr || it.counts->empty()) continue;
        const size_t r = static_cast<size_t>(it.r);
        size_t bytes = 0;
        for (int64_t c : *it.counts) bytes += static_cast<size_t>(std::max<int64_t>(c, 0)) * r * D;
        if (bytes == 0) continue;
        // The staging never exceeds what the heap can spare: this runs while
        // the wave's own working memory is still allocated.
        const size_t spare = heap_.largest_free() / 2;
        const size_t budget = std::max<size_t>(std::min(kFillStagingBudget, spare), 1);
        if (!srcs.empty() && group_raw + bytes > budget) flush();
        const size_t first = srcs.size();
        size_t row = 0;
        for (int64_t c : *it.counts) {
            const size_t rows = static_cast<size_t>(std::max<int64_t>(c, 0));
            if (rows > 0) {
                srcs.push_back(it.t2t + row * r);
                elems.push_back(rows * r);
                slot_rows.push_back(static_cast<int>(rows));
                slot_box.push_back(static_cast<int>(group_of.size()));  // the box being packed
            }
            row += rows;
        }
        if (srcs.size() > first) group_of.push_back(GroupItem{it.st, first, srcs.size()});
        group_raw += bytes;
    }
    flush();
#else
    (void)items; (void)stream;
#endif
}

// ---------------------------------------------------------------------------
// One group of a wave's boxes: their X_NR compressed, then the group's
// factors gathered into a device block of its own.  A wave is kept in groups
// so the compressor's staging stays bounded - it needs room for the whole
// group at once (bitcompMaxBuflen is at least the input size however small
// the chunks are), and a whole leaf wave does not fit beside the level's own
// working set.  Returns false when the keep failed and the level dropped it.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::keep_solve_factor_group(
    const std::vector<const WaveBox*>& kept_boxes, char* d_result, const std::vector<char*>& d_temp2,
    cudaStream_t stream) {
    std::vector<KeptSource> srcs;
    srcs.reserve(kept_boxes.size());
    for (const WaveBox* wb : kept_boxes) {
        const char* source = d_temp2[static_cast<size_t>(wb->source_group)];
        KeptSource s;
        s.morton = wb->box->morton_index;
        s.st = wb->state;
        s.t = d_result + wb->off_t;
        s.xsr = d_result + wb->off_temp1;
        s.xnr = source + wb->off_temp2;
        s.lu = d_result + wb->off_xrr;
        s.piv = d_result + wb->off_piv;
        if (nonsym_) {  // temp3 and temp4, the solve's X_RS and X_RN
            s.xrs = d_result + wb->off_temp3;
            s.xrn = source + wb->off_temp4;
        }
        s.k = wb->k;
        s.r = wb->r;
        s.ntot = wb->ntot;
        // This box's fill source is the same matrix transposed, and the wave
        // compressed it a moment ago (compress_fill_sources runs first).
        // Pledge that copy rather than make a second one.
        s.pledge = can_pledge_fill_source(*wb);
        srcs.push_back(s);
    }
    return keep_sources_or_halves(srcs, stream);
}

// One kept block for `srcs_in`: per box T and X_SR, X_NR (unless pledged),
// the LU of X_RR and its pivots; with H2_compress_factors on, T, X_SR and
// X_NR go through one compressor batch (the LU and pivots stay exact).  A
// wave's boxes come through keep_solve_factor_group, an owner level's landed
// boxes through land_full_factors.  False when the keep fails, after the
// rescue: nothing was kept, and keep_failed_ is set.
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::keep_factor_sources(
    const std::vector<KeptSource>& srcs_in, cudaStream_t stream) {
    const size_t D = sizeof(S), I = sizeof(int);
    struct Span { size_t off; const char* src; size_t bytes; };
    std::vector<Span> spans;
    struct Made { int64_t morton; size_t T, xsr, xnr, lu, piv; int k, r, ntot;
                  size_t xnr_bytes; double xnr_delta;
                  uint32_t xnr_chunk0; uint32_t xnr_nchunks; bool xnr_pledged;
                  size_t t_bytes, xsr_bytes; uint32_t t_chunk0, t_nchunks, xsr_chunk0, xsr_nchunks;
                  size_t xrs_t = SIZE_MAX, xrn_t = SIZE_MAX; };
    std::vector<Made> made;
    // (a nonsymmetric matrix: X_RS and X_RN transposed into the block, after its spans)
    struct Flip { const char* src; size_t off; int m, n; };
    std::vector<Flip> flips;

    std::vector<CompressedXnr> comp;
    std::vector<int> comp_at(srcs_in.size(), -1);  // box -> its block of comp
    std::vector<int> t_at(srcs_in.size(), -1), xsr_at(srcs_in.size(), -1);  // T and X_SR likewise
    bool compressed = false;
#ifdef H2_HAVE_NVCOMP
    if (compress_factors_ >= 1 && !srcs_in.empty()) {
        std::vector<const S*> srcs;
        std::vector<size_t> elems;
        std::vector<double> tols;  // per block: T and X_SR take the finer step
        srcs.reserve(3 * srcs_in.size());
        elems.reserve(3 * srcs_in.size());
        tols.reserve(3 * srcs_in.size());
        for (size_t bi = 0; bi < srcs_in.size(); ++bi) {
            const KeptSource& ks = srcs_in[bi];
            const size_t nt = static_cast<size_t>(ks.ntot), r = static_cast<size_t>(ks.r);
            const size_t k = static_cast<size_t>(ks.k);
            // T and X_SR (k x r each) go in the same batch: small next to
            // X_NR, but once that is compressed they are a quarter of what
            // the solve keeps
            t_at[bi] = static_cast<int>(srcs.size());
            srcs.push_back(reinterpret_cast<const S*>(ks.t));
            elems.push_back(k * r);
            tols.push_back(tolerance_ * kSmallFactorToleranceFraction);
            xsr_at[bi] = static_cast<int>(srcs.size());
            srcs.push_back(reinterpret_cast<const S*>(ks.xsr));
            elems.push_back(k * r);
            tols.push_back(tolerance_ * kSmallFactorToleranceFraction);
            if (ks.pledge) continue;
            comp_at[bi] = static_cast<int>(srcs.size());
            srcs.push_back(reinterpret_cast<const S*>(ks.xnr));
            elems.push_back(nt * r);
            tols.push_back(tolerance_ * kXnrToleranceFraction);
        }
        if (!srcs.empty()) {
            compressed = xnr_compressor_.compress(
                srcs, elems, tolerance_ * kXnrToleranceFraction, comp, stream, nullptr, &tols);
            if (compressed) ++eliminator_stats().xnr_waves_compressed;
            else ++eliminator_stats().xnr_waves_plain;
        }
    }
#endif

    size_t total = 0;
    // One factor of the block: compressed, one span per chunk, end to end,
    // each at the 256-byte start bitcomp needs to read it, and those padded
    // lengths go into the level's pool, which is how the solve walks to a
    // chunk; otherwise the raw bytes.  Returns the factor's offset.
    auto place = [&](int at_comp, const char* raw, size_t raw_bytes, size_t& bytes, uint32_t& chunk0,
                     uint32_t& nchunks) {
        const size_t off = total;
        if (compressed && at_comp >= 0 && raw_bytes > 0 && comp[static_cast<size_t>(at_comp)].bytes > 0) {
            const CompressedXnr& ci = comp[static_cast<size_t>(at_comp)];
            const auto& cks = xnr_compressor_.chunks();
            chunk0 = static_cast<uint32_t>(kept_.chunk_strides.size());
            nchunks = static_cast<uint32_t>(ci.nchunks);
            size_t at = total;
            for (size_t c = 0; c < ci.nchunks; ++c) {
                const auto& ck = cks[ci.chunk0 + c];
                spans.push_back({at, xnr_compressor_.staging() + ck.offset, ck.stride});
                kept_.chunk_strides.push_back(static_cast<uint32_t>(ck.stride));
                at += ck.stride;
            }
            bytes = ci.bytes;  // the chunks' padded lengths, summed
            total = align_up(at);
        } else {
            if (raw_bytes > 0) spans.push_back({total, raw, raw_bytes});
            total = align_up(total + raw_bytes);
        }
        return off;
    };
    for (size_t bi = 0; bi < srcs_in.size(); ++bi) {
        const KeptSource& ks = srcs_in[bi];
        const size_t k = static_cast<size_t>(ks.k), r = static_cast<size_t>(ks.r), nt = static_cast<size_t>(ks.ntot);
        Made m{ks.morton, 0, 0, 0, 0, 0, ks.k, ks.r, ks.ntot, 0, 0.0, 0, 0, false, 0, 0, 0, 0, 0, 0};
        m.T = place(t_at[bi], ks.t, k * r * D, m.t_bytes, m.t_chunk0, m.t_nchunks);
        m.xsr = place(xsr_at[bi], ks.xsr, k * r * D, m.xsr_bytes, m.xsr_chunk0, m.xsr_nchunks);
        if (ks.pledge) {
            // No room taken here: the transition transposes the fill source's
            // bytes and puts the kept X_NR in a block of its own.
            m.xnr = total;
            m.xnr_pledged = true;
        } else {
            m.xnr = place(comp_at[bi], ks.xnr, nt * r * D, m.xnr_bytes, m.xnr_chunk0, m.xnr_nchunks);
            if (m.xnr_bytes > 0) m.xnr_delta = comp[static_cast<size_t>(comp_at[bi])].delta;
        }
        m.lu = total;  spans.push_back({total, ks.lu, r * r * D});   total = align_up(total + r * r * D);
        m.piv = total; spans.push_back({total, ks.piv, r * I});      total = align_up(total + r * I);
        if (ks.xrs != nullptr && k > 0) {
            m.xrs_t = total;
            flips.push_back(Flip{ks.xrs, total, ks.r, ks.k});
            total = align_up(total + k * r * D);
        }
        if (ks.xrn != nullptr && nt > 0) {
            m.xrn_t = total;
            flips.push_back(Flip{ks.xrn, total, ks.r, ks.ntot});
            total = align_up(total + nt * r * D);
        }
        {
            auto& es = eliminator_stats();
            es.xnr_raw_bytes += static_cast<double>(align_up(nt * r * D));
            // a pledged box's X_NR is its fill source until the transition, so
            // that is what it stores; the handover corrects the difference
            es.xnr_stored_bytes += static_cast<double>(
                m.xnr_pledged ? ks.st->t2t_packed_bytes
                              : (m.xnr_bytes > 0 ? align_up(m.xnr_bytes) : align_up(nt * r * D)));
            es.other_raw_bytes += static_cast<double>(
                2 * align_up(k * r * D) + align_up(r * r * D) + align_up(r * I));
            es.tx_raw_bytes += static_cast<double>(2 * align_up(k * r * D));
            es.tx_stored_bytes += static_cast<double>(
                (m.t_bytes > 0 ? align_up(m.t_bytes) : align_up(k * r * D)) +
                (m.xsr_bytes > 0 ? align_up(m.xsr_bytes) : align_up(k * r * D)));
        }
        made.push_back(m);
    }
    if (total == 0) return true;
    // (a level kept by the up-front decision: within its budget; otherwise
    // under the keep fraction of the heap)
    const bool over = keep_budgeted_ ? kept_.bytes + static_cast<double>(total) > keep_budget_
                                     : static_cast<double>(heap_.used() + total) >
                                           kSolveKeepFraction * static_cast<double>(heap_.capacity());
    if (keep_may_split_ && !over && heap_.largest_free() < total && srcs_in.size() > 1) {
        // within the budget, short only of one range for the whole group:
        // the caller tries its halves, each in a block of its own
#ifdef H2_HAVE_NVCOMP
        if (compressed) xnr_compressor_.release();
#endif
        keep_retry_smaller_ = true;
        ++eliminator_stats().keep_splits;
        return false;
    }
    if (over || heap_.largest_free() < total) {
        if (env::trace("wave")) {  // what the wave asked for against what the heap had
            int r = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &r);
            std::printf("  [gpu] keep failed on rank %d: wave wants %.1f MB, largest free %.1f MB, heap %.2f of "
                        "%.2f GB, kept so far %.2f GB, budget %.2f GB (%s)\n", r,
                        static_cast<double>(total) / 1e6, static_cast<double>(heap_.largest_free()) / 1e6,
                        static_cast<double>(heap_.used()) / 1e9, static_cast<double>(heap_.capacity()) / 1e9,
                        kept_.bytes / 1e9, keep_budget_ / 1e9, over ? "over budget" : "no contiguous range");
            std::fflush(stdout);
        }
#ifdef H2_HAVE_NVCOMP
        if (compressed) xnr_compressor_.release();
#endif
        rescue_skipped_xnr(true);
        drop_kept(true, "keep failed mid-level");
        keep_failed_ = true;
        keep_failure_ = over ? 1 : 2;
        return false;
    }
    char* block = heap_.alloc_resident(total, ResidentKind::KeptFactors);
    if (keep_failed_) {  // (the reclaimer dropped the kept store to find this room: a nonsymmetric level)
#ifdef H2_HAVE_NVCOMP
        if (compressed) xnr_compressor_.release();
#endif
        heap_.free(block);
        return false;
    }
    kept_.blocks.push_back(block);
    kept_.bytes += static_cast<double>(total);
    for (size_t bi = 0; bi < made.size(); ++bi) {  // the pledged X_NR is kept too
        if (made[bi].xnr_pledged) kept_.bytes += static_cast<double>(srcs_in[bi].st->t2t_packed_bytes);
    }
    // the span list is a few KB: a pageable copy (staged before the call
    // returns) into a buffer reserved at level start, so no allocation lands
    // in the wave's pipeline
    std::vector<CopySpan> items;
    int64_t max_words = 0;
    for (const Span& sp : spans) {
        items.push_back(CopySpan{sp.src, block + sp.off, static_cast<int64_t>(sp.bytes / 4)});
        max_words = std::max(max_words, static_cast<int64_t>(sp.bytes / 4));
    }
    // (a nonsymmetric matrix: X_RS and X_RN's transposes, their list after the spans')
    std::vector<TransposeItem> tr;
    int max_tm = 0, max_tn = 0;
    for (const Flip& f : flips) {
        tr.push_back(TransposeItem{reinterpret_cast<const S*>(f.src), f.m, reinterpret_cast<S*>(block + f.off), f.n,
                                   f.m, f.n});
        max_tm = std::max(max_tm, f.m);
        max_tn = std::max(max_tn, f.n);
    }
    const size_t span_bytes = align_up(items.size() * sizeof(CopySpan));
    const size_t list_bytes = span_bytes + tr.size() * sizeof(TransposeItem);
    if (list_bytes > solve_meta_device_.capacity()) {
        check_cuda(cudaStreamSynchronize(stream), "kept factors");  // the old list may still be read
        solve_meta_device_.reserve(list_bytes);
    }
    char* d_list = static_cast<char*>(solve_meta_device_.reserve(list_bytes));
    std::vector<char> h_list(list_bytes);
    std::memcpy(h_list.data(), items.data(), items.size() * sizeof(CopySpan));
    if (!tr.empty()) std::memcpy(h_list.data() + span_bytes, tr.data(), tr.size() * sizeof(TransposeItem));
    check_cuda(cudaMemcpyAsync(d_list, h_list.data(), list_bytes, cudaMemcpyHostToDevice, stream), "kept factors");
    launch_copy_spans(reinterpret_cast<const CopySpan*>(d_list), static_cast<int>(items.size()), max_words, stream);
    if (!tr.empty()) {
        launch_transpose(reinterpret_cast<const TransposeItem*>(d_list + span_bytes), static_cast<int>(tr.size()), max_tm,
                         max_tn, stream);
    }
#ifdef H2_HAVE_NVCOMP
    if (compressed) {
        // the gather reads the compressor's staging: it may only go back to
        // the heap once the copy is done
        check_cuda(cudaStreamSynchronize(stream), "kept factors (compressed)");
        xnr_compressor_.release();
    }
#endif
    for (size_t bi = 0; bi < made.size(); ++bi) {
        const Made& m = made[bi];
        KeptSolveBox kb;
        kb.T = block + m.T;
        kb.xsr = block + m.xsr;
        kb.xnr = m.xnr_pledged ? nullptr : block + m.xnr;
        kb.xnr_pending = m.xnr_pledged;
        kb.lu = block + m.lu;
        kb.ipiv = reinterpret_cast<const int*>(block + m.piv);
        kb.xrs_t = m.xrs_t != SIZE_MAX ? block + m.xrs_t : nullptr;
        kb.xrn_t = m.xrn_t != SIZE_MAX ? block + m.xrn_t : nullptr;
        kb.k = m.k;
        kb.r = m.r;
        kb.ntot = m.ntot;
        kb.xnr_bytes = m.xnr_bytes;
        kb.xnr_delta = m.xnr_delta;
        kb.xnr_chunk0 = m.xnr_chunk0;
        kb.xnr_nchunks = m.xnr_nchunks;
        kb.t_bytes = m.t_bytes;
        kb.t_chunk0 = m.t_chunk0;
        kb.t_nchunks = m.t_nchunks;
        kb.xsr_bytes = m.xsr_bytes;
        kb.xsr_chunk0 = m.xsr_chunk0;
        kb.xsr_nchunks = m.xsr_nchunks;
        kept_.boxes[m.morton] = kb;
        if (BoxState* st = srcs_in[bi].st; st != nullptr) {
            st->xrr_lu = reinterpret_cast<const S*>(kb.lu);   // X_RR_full = P L U of this
            st->xrr_piv = kb.ipiv;
        }
        if (m.xnr_pledged) {
            BoxState* st = srcs_in[bi].st;
            st->t2t_pledged = true;
            pledged_xnr_.push_back(PledgedXnr{m.morton, st, m.ntot, m.r, st->t2t_packed_bytes});
        }
    }
    return true;
}

// Every box's X_NR^T store.  The owner pass of the wave that built it is its
// only reader, so this runs at the end of each wave.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::release_owner_sources() {
    if (!heap_initialized_) return;
    auto& stats = eliminator_stats();
    auto drop = [&](BoxState& st) {
        if (st.xnrt_packed == nullptr) return;
        stats.persist_bytes -= static_cast<double>(st.xnrt_packed_bytes);
        heap_.free(st.xnrt_packed);
        st.xnrt_packed = nullptr;
        st.xnrt_packed_bytes = 0;
        st.xnrt_slots.clear();
        st.xnrt_strides.clear();
    };
    for (BoxState& st : states_) drop(st);
    for (BoxState& st : remote_states_) drop(st);
}

// ---------------------------------------------------------------------------
// The kept X_NR of the pledged boxes, made from their fill sources.
//
// temp2 is stored once instead of twice.  Through the level it is the fill
// source: temp2^T, r x ntot, compressed a slot - one one-hop neighbor - at a
// time, which is the form and the granularity the sketches read.  Here, the
// level done and the heap at its emptiest, each pledged box's slots are
// decompressed end to end (they are consecutive column ranges, so that is the
// whole r x ntot), transposed into the ntot x r the solve's kernels read, and
// compressed again into a block of the kept store.
//
// The second quantization is at the same tolerance as the first, on blocks
// cut differently (whole matrices, not slots), so it can add up to another
// delta/2 to a kept factor.  Groups are bounded by kKeptTransposeBudget; a
// group that does not fit is retried a box at a time, and a box that still
// does not fit goes round through the host.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::transpose_kept_xnr() {
#ifdef H2_HAVE_NVCOMP
    if (pledged_xnr_.empty()) return;
    std::vector<PledgedXnr> todo;
    todo.swap(pledged_xnr_);  // (a throw must not leave them to be tried twice)
    KeptSolveLevel& kl = committed_kept_ != nullptr ? *committed_kept_ : kept_;
    const size_t D = sizeof(S);
    cudaStream_t stream = Context::instance().stream();
    auto& stats = eliminator_stats();
    const auto t0 = std::chrono::steady_clock::now();
    size_t at = 0;
    while (at < todo.size()) {
        size_t raw = 0, end = at;
        while (end < todo.size()) {
            const size_t bytes = align_up(static_cast<size_t>(todo[end].ntot) *
                                          static_cast<size_t>(todo[end].r) * D);
            if (end > at && raw + bytes > kKeptTransposeBudget) break;
            raw += bytes;
            ++end;
        }
        if (!transpose_kept_xnr_group(todo, at, end, raw, kl, stream)) {
            for (size_t i = at; i < end; ++i) {
                const size_t one = align_up(static_cast<size_t>(todo[i].ntot) *
                                            static_cast<size_t>(todo[i].r) * D);
                if (transpose_kept_xnr_group(todo, i, i + 1, one, kl, stream)) continue;
                if (kept_xnr_through_host(todo[i], kl)) continue;
                throw std::runtime_error("H2 GPU: no room for the kept X_NR of box " +
                                         std::to_string(todo[i].morton));
            }
        }
        at = end;
    }
    stats.kx_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
#endif
}

// One group of pledged boxes.  False means it did not fit, with nothing
// changed: the box keeps its fill source and the caller tries a smaller group
// or the host.
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::transpose_kept_xnr_group(
    const std::vector<PledgedXnr>& todo, size_t first, size_t last, size_t raw, KeptSolveLevel& kl,
    cudaStream_t stream) {
#ifdef H2_HAVE_NVCOMP
    const size_t D = sizeof(S);
    auto& stats = eliminator_stats();
    const size_t n = last - first;
    if (n == 0) return true;
    char* a = heap_.try_alloc(std::max<size_t>(raw, 1), ResidentKind::Codec);     // temp2^T
    char* b = a != nullptr ? heap_.try_alloc(std::max<size_t>(raw, 1), ResidentKind::Codec) : nullptr;  // temp2
    if (a == nullptr || b == nullptr) {
        heap_.free(b);
        heap_.free(a);
        return false;
    }
    std::vector<size_t> off(n, 0);
    {
        size_t o = 0;
        for (size_t i = first; i < last; ++i) {
            off[i - first] = o;
            o = align_up(o + static_cast<size_t>(todo[i].ntot) * static_cast<size_t>(todo[i].r) * D);
        }
    }
    // 1. every slot of every box of the group, decompressed in one batch
    SolveXnrWave w;
    for (size_t i = first; i < last; ++i) {
        const PledgedXnr& p = todo[i];
        const BoxState& st = *p.st;
        size_t so = off[i - first];
        w.boxes.push_back(static_cast<int>(i - first));
        w.h_box_dsts.push_back(a + so);
        for (const auto& sl : st.t2t_slots) {
            const char* src = packed_of(st) + sl.off;
            for (uint32_t c = 0; c < sl.nchunks; ++c) {
                w.h_srcs.push_back(src);
                w.h_dsts.push_back(a + so + static_cast<size_t>(c) * kXnrChunkBytes);
                src += st.t2t_strides[sl.chunk0 + c];
            }
            so += static_cast<size_t>(sl.rows) * static_cast<size_t>(p.r) * D;
        }
    }
    if (w.h_srcs.empty() || !build_solve_xnr_wave(w, heap_, stream)) {
        release_solve_xnr_wave(w);
        if (w.d_srcs != nullptr) heap_.free(reinterpret_cast<char*>(const_cast<void**>(w.d_srcs)));
        heap_.free(b);
        heap_.free(a);
        return false;
    }
    decompress_solve_xnr_wave(w, stream);
    check_cuda(cudaStreamSynchronize(stream), "kept X_NR: decompress");
    release_solve_xnr_wave(w);
    heap_.free(reinterpret_cast<char*>(const_cast<void**>(w.d_srcs)));
    // 2. transposed into the shape the solve reads
    {
        std::vector<TransposeItemT<S>> tr;
        tr.reserve(n);
        int max_m = 0, max_n = 0;
        for (size_t i = first; i < last; ++i) {
            const PledgedXnr& p = todo[i];
            const size_t o = off[i - first];
            tr.push_back(TransposeItemT<S>{reinterpret_cast<const S*>(a + o), p.r,
                                           reinterpret_cast<S*>(b + o), p.ntot, p.r, p.ntot});
            max_m = std::max(max_m, p.r);
            max_n = std::max(max_n, p.ntot);
        }
        const size_t list_bytes = tr.size() * sizeof(TransposeItemT<S>);
        if (list_bytes > t2t_meta_device_.capacity()) {
            check_cuda(cudaStreamSynchronize(stream), "kept X_NR: transposes");
            t2t_meta_device_.reserve(list_bytes);
        }
        auto* d_tr = static_cast<TransposeItemT<S>*>(t2t_meta_device_.reserve(list_bytes));
        check_cuda(cudaMemcpyAsync(d_tr, tr.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                   "kept X_NR: transposes");
        launch_transpose(d_tr, static_cast<int>(tr.size()), max_m, max_n, stream);
    }
    // 3. compressed again, into the kept store's layout
    std::vector<const S*> srcs;
    std::vector<size_t> elems;
    srcs.reserve(n);
    elems.reserve(n);
    for (size_t i = first; i < last; ++i) {
        srcs.push_back(reinterpret_cast<const S*>(b + off[i - first]));
        elems.push_back(static_cast<size_t>(todo[i].ntot) * static_cast<size_t>(todo[i].r));
    }
    std::vector<CompressedXnr> comp;
    if (!xnr_compressor_.compress(srcs, elems, tolerance_ * kXnrToleranceFraction, comp, stream)) {
        xnr_compressor_.release();
        heap_.free(b);
        heap_.free(a);
        return false;
    }
    size_t total = 0;
    std::vector<size_t> dst(n, 0);
    for (size_t i = 0; i < n; ++i) {
        dst[i] = total;
        total = align_up(total + comp[i].bytes);
    }
    char* block = total > 0 ? heap_.try_alloc_resident(total, ResidentKind::KeptFactors) : nullptr;
    if (block == nullptr) {
        xnr_compressor_.release();
        heap_.free(b);
        heap_.free(a);
        return false;
    }
    const auto& cks = xnr_compressor_.chunks();
    std::vector<CopySpan> items;
    int64_t max_words = 0;
    for (size_t i = 0; i < n; ++i) {
        const PledgedXnr& p = todo[first + i];
        auto it = kl.boxes.find(p.morton);
        if (it == kl.boxes.end()) {
            throw std::runtime_error("H2 GPU: pledged box " + std::to_string(p.morton) +
                                     " is not in the kept store");
        }
        KeptSolveBox& kb = it->second;
        kb.xnr = block + dst[i];
        kb.xnr_bytes = comp[i].bytes;
        kb.xnr_delta = comp[i].delta;
        kb.xnr_chunk0 = static_cast<uint32_t>(kl.chunk_strides.size());
        kb.xnr_nchunks = static_cast<uint32_t>(comp[i].nchunks);
        kb.xnr_pending = false;
        size_t at = dst[i];
        for (size_t c = 0; c < comp[i].nchunks; ++c) {
            const auto& ck = cks[comp[i].chunk0 + c];
            items.push_back(CopySpan{xnr_compressor_.staging() + ck.offset, block + at,
                                     static_cast<int64_t>(ck.stride / 4)});
            max_words = std::max(max_words, static_cast<int64_t>(ck.stride / 4));
            kl.chunk_strides.push_back(static_cast<uint32_t>(ck.stride));
            at += ck.stride;
        }
    }
    kl.blocks.push_back(block);
    {
        double was = 0.0;
        for (size_t i = first; i < last; ++i) was += static_cast<double>(todo[i].bytes);
        kl.bytes += static_cast<double>(total) - was;
        stats.xnr_stored_bytes += static_cast<double>(total) - was;
    }
    {
        const size_t list_bytes = items.size() * sizeof(CopySpan);
        if (list_bytes > solve_meta_device_.capacity()) {
            check_cuda(cudaStreamSynchronize(stream), "kept X_NR: gather");
            solve_meta_device_.reserve(list_bytes);
        }
        auto* d_items = static_cast<CopySpan*>(solve_meta_device_.reserve(list_bytes));
        check_cuda(cudaMemcpyAsync(d_items, items.data(), list_bytes, cudaMemcpyHostToDevice, stream),
                   "kept X_NR: gather");
        launch_copy_spans(d_items, static_cast<int>(items.size()), max_words, stream);
    }
    check_cuda(cudaStreamSynchronize(stream), "kept X_NR");  // the staging goes back now
    xnr_compressor_.release();
    heap_.free(b);
    heap_.free(a);
    // 4. the fill source's copy is the one that goes
    for (size_t i = first; i < last; ++i) {
        BoxState& st = *todo[i].st;
        if (st.t2t_packed != nullptr) {
            stats.persist_bytes -= static_cast<double>(st.t2t_packed_bytes);
            heap_.free(st.t2t_packed);
        }
        st.t2t_packed = nullptr;
        st.t2t_host = nullptr;
        st.t2t_packed_bytes = 0;
        st.t2t_slots.clear();
        st.t2t_strides.clear();
        st.t2t_rows.clear();
        st.t2t_pledged = false;
        st.t2t_reclaimed = false;
        st.temp2t = nullptr;
    }
    stats.kx_boxes += static_cast<int64_t>(n);
    ++stats.kx_groups;
    stats.kx_raw_bytes += static_cast<double>(raw);
    stats.kx_stored_bytes += static_cast<double>(total);
    return true;
#else
    (void)todo; (void)first; (void)last; (void)raw; (void)kl; (void)stream;
    return true;
#endif
}

// A box whose transpose did not fit on the device: its slots are decompressed
// and transposed on the host and come back uncompressed, which is a form the
// solve already reads (xnr_bytes 0).  Bigger, and only for the corner where
// the device had no room at all.
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::kept_xnr_through_host(const PledgedXnr& p,
                                                                             KeptSolveLevel& kl) {
#ifdef H2_HAVE_NVCOMP
    BoxState& st = *p.st;
    if (packed_of(st) == nullptr || st.t2t_slots.empty()) return false;
    const size_t D = sizeof(S);
    const size_t nt = static_cast<size_t>(p.ntot), r = static_cast<size_t>(p.r);
    const size_t bytes = nt * r * D;
    std::vector<char> packed(st.t2t_packed_bytes);
    check_cuda(cudaMemcpy(packed.data(), packed_of(st), packed.size(), cudaMemcpyDefault),
               "kept X_NR through the host");
    std::vector<S> t(nt * r), x(nt * r);
    size_t o = 0;
    for (const auto& sl : st.t2t_slots) {
        decompress_xnr_host<S>(packed.data() + sl.off, t.data() + o,
                               st.t2t_strides.data() + sl.chunk0, sl.nchunks);
        o += static_cast<size_t>(sl.rows) * r;
    }
    for (size_t j = 0; j < r; ++j) {
        for (size_t i = 0; i < nt; ++i) x[i + j * nt] = t[j + i * r];
    }
    char* blk = heap_.try_alloc_resident(bytes, ResidentKind::KeptFactors);
    if (blk == nullptr) return false;
    check_cuda(cudaMemcpy(blk, x.data(), bytes, cudaMemcpyHostToDevice), "kept X_NR through the host");
    auto it = kl.boxes.find(p.morton);
    if (it == kl.boxes.end()) {
        heap_.free(blk);
        return false;
    }
    KeptSolveBox& kb = it->second;
    kb.xnr = blk;  // uncompressed, which is what xnr_bytes 0 means
    kb.xnr_bytes = 0;
    kb.xnr_delta = 0.0;
    kb.xnr_chunk0 = 0;
    kb.xnr_nchunks = 0;
    kb.xnr_pending = false;
    kl.blocks.push_back(blk);
    kl.bytes += static_cast<double>(bytes) - static_cast<double>(p.bytes);
    auto& stats = eliminator_stats();
    stats.xnr_stored_bytes += static_cast<double>(bytes) - static_cast<double>(p.bytes);
    if (st.t2t_packed != nullptr) {
        stats.persist_bytes -= static_cast<double>(st.t2t_packed_bytes);
        heap_.free(st.t2t_packed);
    }
    st.t2t_packed = nullptr;
    st.t2t_host = nullptr;
    st.t2t_packed_bytes = 0;
    st.t2t_slots.clear();
    st.t2t_strides.clear();
    st.t2t_rows.clear();
    st.t2t_pledged = false;
    st.t2t_reclaimed = false;
    st.temp2t = nullptr;
    ++stats.kx_host;
    ++stats.kx_boxes;
    return true;
#else
    (void)p; (void)kl;
    return false;
#endif
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::keep_solve_factors(
    const std::vector<WaveBox>& boxes, const std::vector<int>& active, char* d_result,
    const std::vector<char*>& d_temp2, cudaStream_t stream) {
    const size_t D = sizeof(S);

    // ---- this wave's boxes whose factors are worth keeping
    std::vector<const WaveBox*> kept_boxes, home_ghosts;
    for (int i : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(i)];
        if (wb.k == 0 || wb.r == 0) continue;
        if (ghost_factors_elsewhere(wb.box->morton_index)) {
            if (owner_ == nullptr && wb.state != nullptr) home_ghosts.push_back(&wb);
            continue;
        }
        kept_boxes.push_back(&wb);
    }
    keep_ghost_lu(home_ghosts, d_result, stream);
    if (kept_boxes.empty()) return;
    if (compress_factors_ < 1) ++eliminator_stats().xnr_not_requested;

    // ---- in groups whose X_NR fits the staging budget, so neither the
    // compressor's scratch nor the block the group is gathered into is as
    // large as the whole wave.  A group takes at least one box, whatever its
    // size.  Without compression a group is the whole wave, as before.
    size_t at = 0;
    while (at < kept_boxes.size()) {
        std::vector<const WaveBox*> group;
        size_t raw = 0;
        while (at < kept_boxes.size()) {
            const WaveBox* wb = kept_boxes[at];
            const size_t bytes = static_cast<size_t>(wb->ntot) * static_cast<size_t>(wb->r) * D;
            if (compress_factors_ >= 1 && !group.empty() && raw + bytes > kXnrStagingBudget) break;
            group.push_back(wb);
            raw += bytes;
            ++at;
        }
        if (!keep_solve_factor_group(group, d_result, d_temp2, stream)) return;
    }
}

// Kept factors, else their halves (2026-10-11): a group's kept factors are
// one block, and on a fragmented heap there was no range for 1.1 GB with 6
// GB free (CA3 at 8x, the landed owners' factors).  The keep then failed for
// the level, and with it the chunked production, so the next wave asked for
// its outputs whole.
template<typename CoordType, typename DataType, typename KernelType>
bool LevelEliminator<CoordType, DataType, KernelType>::keep_sources_or_halves(
    const std::vector<KeptSource>& srcs, cudaStream_t stream) {
    keep_retry_smaller_ = false;
    keep_may_split_ = true;
    const bool kept = keep_factor_sources(srcs, stream);
    keep_may_split_ = false;
    if (kept) return true;
    if (!keep_retry_smaller_) return false;   // the keep failed (and said so)
    keep_retry_smaller_ = false;
    const size_t half = srcs.size() / 2;
    const std::vector<KeptSource> a(srcs.begin(), srcs.begin() + static_cast<std::ptrdiff_t>(half));
    const std::vector<KeptSource> b(srcs.begin() + static_cast<std::ptrdiff_t>(half), srcs.end());
    return keep_sources_or_halves(a, stream) && keep_sources_or_halves(b, stream);
}


// ---------------------------------------------------------------------------
// Device exchange (multi-rank levels, CUDA-aware MPI).
//
// Sender: the generators of a wave's boxes with a neighbor on another rank
// go into the outbox, one segment per destination rank (GenLayout each),
// copied from the wave's source groups and result block before the copier
// releases them; host headers describe them.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::plan_generators(
    const std::vector<WaveBox>& boxes, const std::vector<int>& active) {
    if (outbox_ != nullptr || !peer_out_.empty()) {
        throw std::runtime_error("LevelEliminator: generators of two waves wait for a transport");
    }
    using Run = GenRun;        // rows of the original X_NR for the destination
    using Piece = GenPiece;
    std::vector<Piece>& pieces = gen_pieces_;
    pieces.clear();
    std::vector<int> peers;
    for (int bi_ : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(bi_)];
        const Box* box = wb.box;
        peers.clear();
        for (int64_t nb : box->one_hop) {
            if (is_local(nb)) continue;
            const int p = owner_of(nb);
            if (std::find(peers.begin(), peers.end(), p) == peers.end()) peers.push_back(p);
        }
        // occupancy: also the ranks holding the box as an assisting box (the
        // unstructured host transport's lazy_generator_requesters)
        auto rq = occupancy_ && !peers.empty() ? requesters_.find(box->morton_index) : requesters_.end();
        if (rq != requesters_.end()) {
            for (int p : rq->second) {
                if (std::find(peers.begin(), peers.end(), p) == peers.end()) peers.push_back(p);
            }
        }
        if (peers.empty()) continue;
        const bool fill = wb.r > 0;
        if (fill && !wb.has_remote) throw std::runtime_error("LevelEliminator: generator without its X_NR copy");
        const auto& skeleton = box->skeleton_indices;
        for (int p : peers) {
            Piece piece{p, 0, &wb, 0, {}};
            if (fill && !nonsym_) {  // (a nonsymmetric generator sends X_RN whole instead)
                for (size_t a = 0; a < box->one_hop.size(); ++a) {
                    const int len = wb.counts[a];
                    if (len == 0 || is_local(box->one_hop[a]) || owner_of(box->one_hop[a]) != p) continue;
                    if (!piece.runs.empty() && piece.runs.back().row + piece.runs.back().len == wb.row0[a]) {
                        piece.runs.back().len += len;
                    } else {
                        piece.runs.push_back(Run{wb.row0[a], len});
                    }
                    piece.nx += len;
                }
            }
            const GenLayout lay(static_cast<size_t>(wb.ntot), static_cast<size_t>(wb.r), static_cast<size_t>(wb.k),
                                static_cast<size_t>(piece.nx), nonsym_);
            PeerOut& po = peer_out_[p];
            std::vector<int64_t>& h = po.header;
            h.push_back(box->morton_index);
            h.push_back(last_wave_index_);
            h.push_back(fill ? wb.r : 0);
            h.push_back(fill ? wb.k : 0);
            h.push_back(fill ? wb.ntot : 0);
            h.push_back(piece.nx);
            h.push_back(fill ? static_cast<int64_t>(po.bytes) : -1);
            h.push_back(static_cast<int64_t>(box->one_hop.size()));
            h.insert(h.end(), box->one_hop.begin(), box->one_hop.end());
            if (fill) {
                h.insert(h.end(), wb.counts.begin(), wb.counts.end());
                h.push_back(static_cast<int64_t>(skeleton.size()));
                h.insert(h.end(), skeleton.begin(), skeleton.end());
                piece.off = po.bytes;
                pieces.push_back(std::move(piece));
                po.bytes += lay.bytes;
            }
        }
    }
    size_t total = 0;
    for (auto& kv : peer_out_) {
        kv.second.offset = total;
        total = align_up(total + kv.second.bytes);
    }
    if (total == 0) return;
    outbox_ = exchange_alloc(total);
}

// The pieces of the boxes in `active`, copied into the outbox the plan
// allocated.  Their temp2 is raw only inside their own production chunk.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::copy_generators(
    const std::vector<WaveBox>& boxes, const std::vector<int>& active, const std::vector<char*>& d_sources,
    const std::vector<char*>& d_temp2, const char* d_xnr, char* d_result, cudaStream_t stream) {
    if (outbox_ == nullptr || gen_pieces_.empty()) return;
    const size_t D = sizeof(S);
    std::unordered_set<const WaveBox*> mine;
    for (int i : active) mine.insert(&boxes[static_cast<size_t>(i)]);
    auto copy = [&](char* to, const char* from, size_t bytes, const char* what) {
        if (bytes > 0) check_cuda(cudaMemcpyAsync(to, from, bytes, cudaMemcpyDefault, stream), what);
    };
    for (const GenPiece& pc : gen_pieces_) {
        if (mine.count(pc.wb) == 0) continue;
        const WaveBox& wb = *pc.wb;
        const size_t r = static_cast<size_t>(wb.r), k = static_cast<size_t>(wb.k), nt = static_cast<size_t>(wb.ntot);
        const size_t nx = static_cast<size_t>(pc.nx);
        const GenLayout lay(nt, r, k, nx, nonsym_);
        char* dst = outbox_ + peer_out_.at(pc.peer).offset + pc.off;
        const char* source = d_sources[static_cast<size_t>(wb.source_group)];
        const char* xnr = d_xnr != nullptr ? d_xnr : source;
        copy(dst, d_temp2[static_cast<size_t>(wb.source_group)] + wb.off_temp2, nt * r * D, "pack temp2");
        size_t row = 0;
        for (const GenRun& run : pc.runs) {  // the destination's rows of every column
            check_cuda(cudaMemcpy2DAsync(dst + lay.xnr + row * D, nx * D, xnr + wb.off_xnr_orig + run.row * D, nt * D,
                                         static_cast<size_t>(run.len) * D, r, cudaMemcpyDefault, stream),
                       "pack X_NR");
            row += static_cast<size_t>(run.len);
        }
        copy(dst + lay.xrs, d_result + wb.off_xrs, r * k * D, "pack X_RS");
        copy(dst + lay.xrr, source + wb.off_xrr_full, r * r * D, "pack X_RR_full");
        if (nonsym_) {
            copy(dst + lay.xrn, reinterpret_cast<const char*>(wb.xrn), r * nt * D, "pack X_RN");
            copy(dst + lay.temp4, d_temp2[static_cast<size_t>(wb.source_group)] + wb.off_temp4, r * nt * D,
                 "pack temp4");
            copy(dst + lay.xsr_solved, d_result + wb.off_temp1, k * r * D, "pack X_SR");
        }
    }
}

// The transport of the Color loop.  Each rank sends each neighbor rank, all
// at once, its sizes, a host header -- the transport number; on the first
// transport the assisting boxes it holds of that rank (whose owner then
// sends their skeletons as they are eliminated); the generators of the last
// wave; and the skeletons of the rank's assisting boxes eliminated since
// the last transport -- and the generators' bulk, device to device.  The
// receiver installs the generator boxes (metadata), refreshes the assisting
// skeletons, forms the updates of its blocks on the device straight from the
// receive buffer, and copies the generators' temp2 and X_RR_full to the host
// in the background (attached at finish()).  Point data of the assisting
// boxes does not change within a level: it comes once, with the first
// transport.
template<typename CoordType, typename DataType, typename KernelType>
std::chrono::high_resolution_clock::duration LevelEliminator<CoordType, DataType, KernelType>::exchange() {
    using hclock = std::chrono::high_resolution_clock;
    const auto t_start = hclock::now();
    auto& stats = eliminator_stats();
    auto& timers = transport_timers();
    auto seconds_since = [](hclock::time_point t) { return std::chrono::duration<double>(hclock::now() - t).count(); };
    const size_t D = sizeof(S);
    hclock::duration comm_time{};
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    // the outbox is packed, and no queued launch still uses heap memory the
    // receive buffer may take
    check_cuda(cudaStreamSynchronize(stream), "exchange: last wave");

    MPI_Comm comm = tree_->comm;
    const std::vector<int> peers = compute_one_hop_neighbor_ranks(tree_, level_, level_index_);
    const size_t np = peers.size();
    std::unordered_map<int, size_t> peer_index;
    for (size_t i = 0; i < np; ++i) peer_index.emplace(peers[i], i);
    auto index_of = [&](int rank, const char* what) {
        auto it = peer_index.find(rank);
        if (it == peer_index.end()) {
            throw std::runtime_error(std::string("LevelEliminator::exchange: ") + what + " on rank " +
                                     std::to_string(rank) + ", not a neighbor rank");
        }
        return it->second;
    };

    // headers: [requests] [generators] [skeletons], each after its length
    std::vector<std::vector<int64_t>> requests(np), skeletons(np);
    const bool first = !registered_;
    if (first) {
        if (!since_transport_.empty()) {
            throw std::runtime_error("LevelEliminator::exchange: boxes eliminated before the first transport");
        }
        for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) {
            requests[index_of(owner_of(kv.first), "assisting box")].push_back(kv.first);
        }
        for (auto& v : requests) std::sort(v.begin(), v.end());
    }
    for (int64_t m : since_transport_) {
        auto it = requesters_.find(m);
        if (it == requesters_.end()) continue;
        const Box* box = level_.find_local_box(m);
        for (int p : it->second) {
            auto& sk = skeletons[index_of(p, "requester")];
            sk.push_back(m);
            sk.push_back(static_cast<int64_t>(box->skeleton_indices.size()));
            sk.insert(sk.end(), box->skeleton_indices.begin(), box->skeleton_indices.end());
        }
    }
    for (const auto& kv : peer_out_) index_of(kv.first, "generator destination");
    // one header per peer, always: [transport number] [requests] [generators]
    // [skeletons], each list after its length
    std::vector<std::vector<int64_t>> head_out(np), head_in(np);
    std::vector<size_t> bulk_out(np, 0), bulk_off(np, 0);
    const std::vector<int64_t> none;
    for (size_t i = 0; i < np; ++i) {
        auto& h = head_out[i];
        auto append = [&](const std::vector<int64_t>& part) {
            h.push_back(static_cast<int64_t>(part.size()));
            h.insert(h.end(), part.begin(), part.end());
        };
        auto po = peer_out_.find(peers[i]);
        if (po != peer_out_.end()) {
            bulk_out[i] = po->second.bytes;
            bulk_off[i] = po->second.offset;
        }
        h.push_back(exchange_seq_);  // (the receiver checks it is at the same transport)
        append(requests[i]);
        append(po != peer_out_.end() ? po->second.header : none);
        append(skeletons[i]);
    }

    auto check_mpi = [](int err, const char* what) {
        if (err != MPI_SUCCESS) throw std::runtime_error(std::string("LevelEliminator::exchange: ") + what + " failed");
    };
    constexpr int kTagSizes = 720, kTagHead = 721, kTagBulk = 722;
    std::vector<MPI_Request> reqs;
    static const bool trace = env::trace("exchange");
    const int64_t fallbacks0 = stats.exchange_fallbacks;
    const bool outbox_in_arena = outbox_ == nullptr || DeviceHeap::exchange_arena().owns(outbox_);
    // Everything goes out at once: a peer's sizes (its header's length, its
    // bulk's bytes), the header (host) and the generator bulk (device to
    // device); no rank waits for its peers before sending.  A peer's header
    // and bulk are received as soon as its sizes are in, so the transfers
    // from early peers overlap the wait for the late ones.  (Explicit
    // receives throughout: probing the headers for their lengths matched
    // them out of order on this MPI.)
    std::vector<uint64_t> sizes_out(2 * np, 0), sizes_in(2 * np, 0);
    std::vector<MPI_Request> size_reqs(np, MPI_REQUEST_NULL);
    auto t_comm = hclock::now();
    timers.presend += std::chrono::duration<double>(t_comm - t_start).count();  // (the last wave's device tail, the headers)
    for (size_t i = 0; i < np; ++i) {
        check_mpi(MPI_Irecv(&sizes_in[2 * i], 2, MPI_UINT64_T, peers[i], kTagSizes, comm, &size_reqs[i]), "size Irecv");
    }
    for (size_t i = 0; i < np; ++i) {
        sizes_out[2 * i] = head_out[i].size();
        sizes_out[2 * i + 1] = bulk_out[i];
        reqs.emplace_back();
        check_mpi(MPI_Isend(&sizes_out[2 * i], 2, MPI_UINT64_T, peers[i], kTagSizes, comm, &reqs.back()), "size Isend");
        reqs.emplace_back();
        check_mpi(MPI_Isend(head_out[i].data(), static_cast<int>(head_out[i].size()), MPI_INT64_T, peers[i], kTagHead,
                            comm, &reqs.back()), "header Isend");
        if (bulk_out[i] > 0) {
            check_mpi(MPI_Isend_large(outbox_ + bulk_off[i], bulk_out[i], MPI_BYTE, peers[i], kTagBulk, comm, reqs),
                      "bulk Isend");
        }
    }
    std::vector<char*> d_recv(np, nullptr);
    std::vector<size_t> bulk_in(np, 0);
    for (size_t done = 0; done < np; ++done) {
        int idx = MPI_UNDEFINED;
        check_mpi(MPI_Waitany(static_cast<int>(np), size_reqs.data(), &idx, MPI_STATUS_IGNORE), "size Waitany");
        if (idx == MPI_UNDEFINED) throw std::runtime_error("LevelEliminator::exchange: a peer's sizes never came");
        const size_t i = static_cast<size_t>(idx);
        head_in[i].resize(static_cast<size_t>(sizes_in[2 * i]));
        if (head_in[i].empty()) throw std::runtime_error("LevelEliminator::exchange: empty header");
        reqs.emplace_back();
        check_mpi(MPI_Irecv(head_in[i].data(), static_cast<int>(head_in[i].size()), MPI_INT64_T, peers[i], kTagHead,
                            comm, &reqs.back()), "header Irecv");
        bulk_in[i] = static_cast<size_t>(sizes_in[2 * i + 1]);
        if (bulk_in[i] > 0) {
            d_recv[i] = exchange_alloc(bulk_in[i]);
            check_mpi(MPI_Irecv_large(d_recv[i], bulk_in[i], MPI_BYTE, peers[i], kTagBulk, comm, reqs), "bulk Irecv");
        }
    }
    const double t_sizes = seconds_since(t_comm);  // the sizes' arrival: the slowest peer's wave
    timers.sizes += t_sizes;
    const auto t_bulk = hclock::now();
    check_mpi(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE), "exchange Waitall");
    reqs.clear();
    comm_time += hclock::now() - t_comm;
    const double t_payload = seconds_since(t_bulk);
    timers.payload += t_payload;
    for (size_t i = 0; i < np; ++i) {
        if (head_in[i][0] != exchange_seq_) {
            throw std::runtime_error("LevelEliminator::exchange: rank " + std::to_string(tree_->mpi_rank) +
                                     " at transport " + std::to_string(exchange_seq_) + " of level " +
                                     std::to_string(level_index_) + " received the header of transport " +
                                     std::to_string(head_in[i][0]) + " from rank " + std::to_string(peers[i]));
        }
    }
    for (size_t i = 0; i < np; ++i) {
        timers.bytes_sent += static_cast<double>(head_out[i].size() * sizeof(int64_t) + bulk_out[i]);
        timers.bytes_received += static_cast<double>(head_in[i].size() * sizeof(int64_t) + bulk_in[i]);
    }
    if (trace) {
        uint64_t out = 0, in = 0;
        for (size_t i = 0; i < np; ++i) {
            out += bulk_out[i];
            in += bulk_in[i];
        }
        const DeviceHeap& arena = DeviceHeap::exchange_arena();
        std::printf("[exchange] rank %d level %d wave %d: sizes %.3f s, rest %.3f s, out %.3f GB%s, in %.3f GB%s, "
                    "arena used %.2f of %.2f GB\n",
                    tree_->mpi_rank, level_index_, last_wave_index_, t_sizes, t_payload, out / 1e9,
                    outbox_in_arena ? "" : " (main heap)", in / 1e9,
                    stats.exchange_fallbacks > fallbacks0 ? " (main heap)" : "", arena.used() / 1e9,
                    arena.capacity() / 1e9);
        std::fflush(stdout);
    }
    free_block(outbox_);
    outbox_ = nullptr;
    peer_out_.clear();

    // headers: requesters, generator boxes, assisting skeletons
    const auto t_parse = hclock::now();
    std::vector<int64_t> installed;
    std::unordered_map<int64_t, Incoming> incoming;
    std::vector<std::pair<int64_t, std::vector<int64_t>>> fresh;  // assisting skeletons
    for (size_t i = 0; i < np; ++i) {
        const std::vector<int64_t>& h = head_in[i];
        size_t at = 1;  // (after the transport number)
        auto next = [&]() -> int64_t {
            if (at >= h.size()) throw std::runtime_error("LevelEliminator::exchange: truncated header");
            return h[at++];
        };
        const int64_t num_requests = next();
        for (int64_t q = 0; q < num_requests; ++q) {
            const int64_t m = next();
            if (!is_local(m)) throw std::runtime_error("LevelEliminator::exchange: request for a box of another rank");
            requesters_[m].push_back(peers[i]);
        }
        const size_t generators_end = static_cast<size_t>(next()) + at;
        while (at < generators_end) {
            const int64_t m = next();
            const int32_t wave = static_cast<int32_t>(next());
            const int r = static_cast<int>(next()), k = static_cast<int>(next()), ntot = static_cast<int>(next());
            const int nx = static_cast<int>(next());
            const int64_t off = next();
            Box box;
            box.morton_index = m;
            box.one_hop.resize(static_cast<size_t>(next()));
            for (auto& x : box.one_hop) x = next();
            std::vector<int64_t> skeleton;
            if (r > 0) {
                box.deferred_xnn_neighbor_point_counts.resize(box.one_hop.size());
                for (auto& c : box.deferred_xnn_neighbor_point_counts) c = next();
                skeleton.resize(static_cast<size_t>(next()));
                for (auto& x : skeleton) x = next();
                box.skeleton_indices.assign(skeleton.begin(), skeleton.end());
            }
            if (level_.generator_id_to_index.count(m) != 0) continue;
            level_.generator_id_to_index[m] = static_cast<int64_t>(level_.generator_boxes.size());
            level_.generator_boxes.push_back(std::move(box));
            level_.elimination_wave[m] = wave;
            installed.push_back(m);
            if (r > 0) {
                const GenLayout lay(static_cast<size_t>(ntot), static_cast<size_t>(r), static_cast<size_t>(k),
                                    static_cast<size_t>(nx), nonsym_);
                if (off < 0 || static_cast<uint64_t>(off) + lay.bytes > bulk_in[i]) {
                    throw std::runtime_error("LevelEliminator::exchange: generator outside its bulk message");
                }
                incoming.emplace(m, Incoming{r, k, ntot, nx, d_recv[i] + static_cast<size_t>(off)});
                fresh.emplace_back(m, std::move(skeleton));
            }
        }
        if (at != generators_end) throw std::runtime_error("LevelEliminator::exchange: malformed generator header");
        const size_t skeletons_end = static_cast<size_t>(next()) + at;
        while (at < skeletons_end) {
            const int64_t m = next();
            std::vector<int64_t> skeleton(static_cast<size_t>(next()));
            for (auto& x : skeleton) x = next();
            fresh.emplace_back(m, std::move(skeleton));
        }
        if (at != h.size()) throw std::runtime_error("LevelEliminator::exchange: malformed header");
    }
    for (auto& [m, skeleton] : fresh) {
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(m);
        if (it == level_.assisting_box_points_for_kernel_evaluation.end()) continue;
        auto& assist = level_.assisting_boxes[static_cast<size_t>(it->second)];
        assist.skel_indices.assign(skeleton.begin(), skeleton.end());
    }
    timers.deserialize += seconds_since(t_parse);
    if (first) {
        // point data of the assisting boxes, once per level
        const auto t_assisting = hclock::now();
        std::vector<int64_t> need_assist;
        for (const auto& kv : level_.assisting_box_points_for_kernel_evaluation) need_assist.push_back(kv.first);
        if (occupancy_) {
            // every rank enters (an empty process region has no requests but
            // serves its neighbors'), as the unstructured host transport does
            comm_time += exchange_assisting_for_mortons_onehop(tree_, level_, level_index_, peers, need_assist, {},
                                                               true);
        } else if (!need_assist.empty()) {
            comm_time += exchange_assisting_for_mortons_onehop(tree_, level_, level_index_, peers, need_assist);
        }
        kernel_->register_level_coordinates(level_);
        timers.assisting += seconds_since(t_assisting);
    }
    registered_ = true;

    // device updates, straight from the receive buffer
    const auto t_apply = hclock::now();
    sync_remote_boxes();
    apply_remote_generators(installed, &incoming);
    for (size_t i = 0; i < np; ++i) {
        if (d_recv[i] == nullptr) continue;
        // host copies of X_RR_full, and of temp2 unless the transition puts a
        // freed generator back from the packed store (packed_restore_, H2
        // host memory); the copier releases the buffer
        HostCopier::Job job;
        job.device = d_recv[i];
        job.bytes = bulk_in[i];
        for (int64_t m : installed) {
            auto it = incoming.find(m);
            if (it == incoming.end()) continue;
            const Incoming& in = it->second;
            if (in.bulk < d_recv[i] || in.bulk >= d_recv[i] + bulk_in[i]) continue;  // (another peer's buffer)
            const GenLayout lay(static_cast<size_t>(in.ntot), static_cast<size_t>(in.r), static_cast<size_t>(in.k),
                                static_cast<size_t>(in.nx), nonsym_);
            auto rh = std::make_unique<RemoteHost>();
            rh->ntot = in.ntot;
            rh->r = in.r;
            const size_t base = static_cast<size_t>(in.bulk - d_recv[i]);
            const size_t nr = static_cast<size_t>(in.ntot) * in.r * D;
            if (nr > 0 && !packed_restore_) job.segments.push_back(HostCopier::Segment{&rh->temp2, nullptr, base, nr});
            job.segments.push_back(HostCopier::Segment{&rh->xrr, nullptr, base + lay.xrr, static_cast<size_t>(in.r) * in.r * D});
            // (a nonsymmetric matrix: temp4, the host transition's column side)
            if (nonsym_ && nr > 0) job.segments.push_back(HostCopier::Segment{&rh->temp4, nullptr, base + lay.temp4, nr});
            remote_host_[m] = std::move(rh);
        }
        if (job.segments.empty()) {
            free_block(d_recv[i]);  // later launches are ordered after its reads
        } else {
            check_cuda(cudaEventCreateWithFlags(&job.ready, cudaEventDisableTiming), "cudaEventCreate");
            check_cuda(cudaEventRecord(job.ready, stream), "cudaEventRecord");
            exchange_copier_->submit(std::move(job));
        }
    }
    since_transport_.clear();
    last_wave_.clear();
    ++exchange_seq_;
    stats.remote += seconds_since(t_apply);
    stats.exchange += std::chrono::duration<double>(t_apply - t_start).count();
    stats.device_exchange = true;
    return comm_time;
}

// ---------------------------------------------------------------------------
// Component-owner CA levels (owner_plan.hpp; h2_gpu/GPU_CA_PLAN.md, M5).
//
// Sender: after each owned wave, the state of the wave's boxes for their
// recipients goes into a device block of the exchange arena, one piece per
// (box, recipient): the generator (GenLayout: temp2, the recipient's rows of
// the original X_NR, X_RS, X_RR_full) and, for the box's home rank (FULL),
// the solve's factors after it (T, X_SR, the LU of X_RR, its pivots); a
// SKELETON recipient gets a header entry only.  A colour's pieces are
// concatenated per recipient at its transport (owner_exchange), which sends
// every recipient its header and bulk, receives the colour's payloads from
// the owners of the components installed here, installs them wave by wave
// (the skeletons, then the device updates from the generators, as the
// Color transport does for one wave) and lands the FULL factors of this
// rank's boxes eliminated elsewhere in the kept store or the host boxes.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::plan_owner_generators(
    const std::vector<WaveBox>& boxes, const std::vector<int>& active) {
    if (owner_plan_.open) throw std::runtime_error("LevelEliminator: an owner wave's message was never finished");
    owner_plan_ = OwnerWavePlan{};
    std::vector<OwnerPiece>& pieces = owner_plan_.pieces;
    const int colour = last_wave_index_ / std::max(owner_->num_waves, 1);
    // per recipient: the wave message's header [colour][wave][bulk bytes]
    // [entries...] and the bytes of its segment of the wave's block
    std::map<int, std::vector<int64_t>> headers;
    std::map<int, size_t> peer_bytes;
    for (int bi_ : active) {
        const WaveBox& wb = boxes[static_cast<size_t>(bi_)];
        const Box* box = wb.box;
        const auto rit = owner_->recipients.find(box->morton_index);
        if (rit == owner_->recipients.end()) continue;
        const bool fill = wb.r > 0;
        if (fill && !wb.has_remote) throw std::runtime_error("LevelEliminator: owned box without its X_NR copy");
        const auto& skeleton = box->skeleton_indices;
        const auto& redundant = box->redundant_indices;
        for (const auto& [peer, kind] : rit->second) {
            auto hit = headers.find(peer);
            if (hit == headers.end()) {
                hit = headers.emplace(peer, std::vector<int64_t>{colour, last_wave_index_, 0}).first;
                peer_bytes[peer] = 0;
            }
            std::vector<int64_t>& h = hit->second;
            size_t& bytes = peer_bytes[peer];
            h.push_back(box->morton_index);
            h.push_back(last_wave_index_);
            h.push_back(static_cast<int64_t>(kind));
            if (kind == OwnerPayload::SKELETON) {
                h.push_back(static_cast<int64_t>(skeleton.size()));
                h.insert(h.end(), skeleton.begin(), skeleton.end());
                continue;
            }
            OwnerPiece piece{peer, 0, &wb, 0, kind == OwnerPayload::FULL, {}};
            if (fill) {
                for (size_t a = 0; a < box->one_hop.size(); ++a) {
                    const int len = static_cast<int>(wb.counts[a]);
                    if (len == 0 || !owner_->interest_of(peer, box->one_hop[a])) continue;
                    if (!piece.runs.empty() && piece.runs.back().row + piece.runs.back().len == wb.row0[a]) {
                        piece.runs.back().len += len;
                    } else {
                        piece.runs.push_back(OwnerRun{wb.row0[a], len});
                    }
                    piece.nx += len;
                }
            }
            const OwnerLayout lay(static_cast<size_t>(wb.ntot), static_cast<size_t>(wb.r), static_cast<size_t>(wb.k),
                                  static_cast<size_t>(piece.nx), piece.full);
            h.push_back(fill ? wb.r : 0);
            h.push_back(fill ? wb.k : 0);
            h.push_back(fill ? wb.ntot : 0);
            h.push_back(piece.nx);
            h.push_back(fill ? static_cast<int64_t>(bytes) : -1);  // (the offset in the segment)
            h.push_back(static_cast<int64_t>(box->one_hop.size()));
            h.insert(h.end(), box->one_hop.begin(), box->one_hop.end());
            h.push_back(static_cast<int64_t>(skeleton.size()));
            h.insert(h.end(), skeleton.begin(), skeleton.end());
            h.push_back(static_cast<int64_t>(redundant.size()));
            h.insert(h.end(), redundant.begin(), redundant.end());
            if (fill) {
                h.insert(h.end(), wb.counts.begin(), wb.counts.end());
                if (piece.full) {  // (the home's solve reads which neighbors' rows are the full sets)
                    for (size_t a = 0; a < box->one_hop.size(); ++a) {
                        h.push_back(a < box->use_full_set.size() ? static_cast<int64_t>(box->use_full_set[a]) : 1);
                    }
                }
                piece.off = bytes;
                pieces.push_back(std::move(piece));
                bytes += lay.total;
            }
        }
    }
    if (headers.empty()) return;  // (nothing to copy or send: the plan stays closed)
    // one block of the exchange arena for the wave's pieces, a segment per
    // recipient (the sends come straight from it)
    size_t total = 0;
    std::map<int, size_t>& seg_off = owner_plan_.seg_off;
    for (const auto& [peer, bytes] : peer_bytes) {
        seg_off[peer] = total;
        total = align_up(total + bytes);
    }
    OwnerWaveSend& send = owner_plan_.send;
    send.colour = colour;
    if (total > 0) send.block = exchange_alloc(total);
    auto& timers = transport_timers();
    for (auto& [peer, h] : headers) {
        const size_t bytes = peer_bytes.at(peer);
        h[2] = static_cast<int64_t>(bytes);
        timers.bytes_sent += static_cast<double>(h.size() * sizeof(int64_t) + bytes);
        owner_sent_bytes_[colour] += bytes;
        send.messages.push_back({peer, std::move(h), bytes > 0 ? send.block + seg_off.at(peer) : nullptr, bytes});
        ++owner_sent_count_[{peer, colour}];
    }
    owner_plan_.open = true;
}

// The pieces of `active` (a production chunk) into the wave's block: temp2,
// the recipient's rows of the original X_NR, X_RS, X_RR_full, and for the
// home rank the solve's factors and the Schur block after them.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::copy_owner_generators(
    const std::vector<WaveBox>& boxes, const std::vector<int>& active, const std::vector<char*>& d_sources,
    const std::vector<char*>& d_temp2, const char* d_xnr, char* d_result, cudaStream_t stream) {
    if (!owner_plan_.open || owner_plan_.send.block == nullptr) return;
    const size_t D = sizeof(S), I = sizeof(int);
    std::unordered_set<const WaveBox*> mine;
    for (int i : active) mine.insert(&boxes[static_cast<size_t>(i)]);
    auto copy = [&](char* to, const char* from, size_t bytes, const char* what) {
        if (bytes > 0) check_cuda(cudaMemcpyAsync(to, from, bytes, cudaMemcpyDefault, stream), what);
    };
    for (const OwnerPiece& pc : owner_plan_.pieces) {
        if (mine.count(pc.wb) == 0) continue;
        const WaveBox& wb = *pc.wb;
        const size_t r = static_cast<size_t>(wb.r), k = static_cast<size_t>(wb.k), nt = static_cast<size_t>(wb.ntot);
        const size_t nx = static_cast<size_t>(pc.nx);
        const OwnerLayout lay(nt, r, k, nx, pc.full);
        char* dst = owner_plan_.send.block + owner_plan_.seg_off.at(pc.peer) + pc.off;
        const char* source = d_sources[static_cast<size_t>(wb.source_group)];
        const char* xnr = d_xnr != nullptr ? d_xnr : source;
        copy(dst, d_temp2[static_cast<size_t>(wb.source_group)] + wb.off_temp2, nt * r * D, "pack temp2");
        size_t row = 0;
        for (const OwnerRun& run : pc.runs) {  // the recipient's rows of every column
            check_cuda(cudaMemcpy2DAsync(dst + lay.xnr + row * D, nx * D, xnr + wb.off_xnr_orig + run.row * D, nt * D,
                                         static_cast<size_t>(run.len) * D, r, cudaMemcpyDefault, stream),
                       "pack X_NR");
            row += static_cast<size_t>(run.len);
        }
        copy(dst + lay.xrs, d_result + wb.off_xrs, r * k * D, "pack X_RS");
        copy(dst + lay.xrr, source + wb.off_xrr_full, r * r * D, "pack X_RR_full");
        if (pc.full) {
            copy(dst + lay.t, d_result + wb.off_t, k * r * D, "pack T");
            copy(dst + lay.xsr, d_result + wb.off_temp1, k * r * D, "pack X_SR");
            copy(dst + lay.lu, d_result + wb.off_xrr, r * r * D, "pack LU");
            copy(dst + lay.piv, d_result + wb.off_piv, r * I, "pack pivots");
            copy(dst + lay.schur, reinterpret_cast<const char*>(wb.s), k * k * D, "pack Schur");
        }
    }
}

// The wave's message goes to the transport once its last copies are on the
// stream (owner_progress sends it when the event has passed).
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::finish_owner_generators(cudaStream_t stream) {
    if (!owner_plan_.open) return;
    OwnerWaveSend send = std::move(owner_plan_.send);
    owner_plan_ = OwnerWavePlan{};
    check_cuda(cudaEventCreateWithFlags(&send.ready, cudaEventDisableTiming), "cudaEventCreate");
    check_cuda(cudaEventRecord(send.ready, stream), "cudaEventRecord");
    owner_sends_.push_back(std::move(send));
}

// The owner level's transport between the waves (and in a colour's wait):
// the waves' payloads whose device copies are done go out (a header, then
// the bulk in pieces, per recipient), sent ones release their block; other
// owners' wave headers are received as they arrive, their bulk straight
// into the exchange arena, and their colour-end markers (how many wave
// messages a peer sent this colour) noted; receive buffers whose installs
// have run return to the arena.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::owner_progress() {
    MPI_Comm comm = tree_->comm;
    constexpr int kTagHead = 733, kTagBulk = 734, kTagMark = 735;
    constexpr size_t kChunk = size_t{256} << 20;  // (one message per piece, in order)
    auto check_mpi = [](int err, const char* what) {
        if (err != MPI_SUCCESS) throw std::runtime_error(std::string("LevelEliminator::owner_progress: ") + what + " failed");
    };
    for (OwnerWaveSend& w : owner_sends_) {
        if (w.issued) continue;
        const cudaError_t q = cudaEventQuery(w.ready);
        if (q == cudaErrorNotReady) break;  // (later waves' copies are not done either)
        check_cuda(q, "owner send event");
        for (auto& m : w.messages) {
            w.reqs.emplace_back();
            check_mpi(MPI_Isend(m.header.data(), static_cast<int>(m.header.size()), MPI_INT64_T, m.peer, kTagHead, comm,
                                &w.reqs.back()), "header Isend");
            for (size_t off = 0; off < m.bytes; off += kChunk) {
                const size_t piece = std::min(kChunk, m.bytes - off);
                w.reqs.emplace_back();
                check_mpi(MPI_Isend(m.ptr + off, static_cast<int>(piece), MPI_BYTE, m.peer, kTagBulk, comm, &w.reqs.back()),
                          "bulk Isend");
            }
        }
        w.issued = true;
    }
    for (auto it = owner_sends_.begin(); it != owner_sends_.end();) {
        if (!it->issued) break;
        int done = 0;
        check_mpi(MPI_Testall(static_cast<int>(it->reqs.size()), it->reqs.data(), &done, MPI_STATUSES_IGNORE), "send Testall");
        if (!done) {
            ++it;
            continue;
        }
        if (it->block != nullptr) free_block(it->block);
        cudaEventDestroy(it->ready);
        it = owner_sends_.erase(it);
    }
    while (true) {
        int flag = 0;
        MPI_Status status;
        check_mpi(MPI_Iprobe(MPI_ANY_SOURCE, kTagHead, comm, &flag, &status), "header Iprobe");
        if (!flag) break;
        int count = 0;
        MPI_Get_count(&status, MPI_INT64_T, &count);
        OwnerWaveRecv r;
        r.source = status.MPI_SOURCE;
        r.header.resize(static_cast<size_t>(std::max(count, 0)));
        check_mpi(MPI_Recv(r.header.data(), count, MPI_INT64_T, r.source, kTagHead, comm, MPI_STATUS_IGNORE), "header Recv");
        if (count < 3) throw std::runtime_error("LevelEliminator::owner_progress: short header");
        r.colour = static_cast<int>(r.header[0]);
        r.seq = static_cast<int32_t>(r.header[1]);
        r.bytes = static_cast<size_t>(r.header[2]);
        if (r.bytes > 0) {
            r.buf = exchange_alloc(r.bytes);
            if (!DeviceHeap::exchange_arena().owns(r.buf)) {
                // A block of the main heap: its memory is reused in the
                // stream's order (a freed block may still be read by launches
                // in flight), but the receive writes it as the data arrives.
                check_cuda(cudaStreamSynchronize(Context::instance().stream()), "owner receive buffer");
            }
            for (size_t off = 0; off < r.bytes; off += kChunk) {
                const size_t piece = std::min(kChunk, r.bytes - off);
                r.reqs.emplace_back();
                check_mpi(MPI_Irecv(r.buf + off, static_cast<int>(piece), MPI_BYTE, r.source, kTagBulk, comm, &r.reqs.back()),
                          "bulk Irecv");
            }
        }
        ++owner_recv_count_[{r.source, r.colour}];
        owner_recvs_.push_back(std::move(r));
    }
    while (true) {
        int flag = 0;
        MPI_Status status;
        check_mpi(MPI_Iprobe(MPI_ANY_SOURCE, kTagMark, comm, &flag, &status), "marker Iprobe");
        if (!flag) break;
        int64_t mark[2] = {0, 0};
        check_mpi(MPI_Recv(mark, 2, MPI_INT64_T, status.MPI_SOURCE, kTagMark, comm, MPI_STATUS_IGNORE), "marker Recv");
        owner_marks_[{status.MPI_SOURCE, static_cast<int>(mark[0])}] = static_cast<int>(mark[1]);
    }
    for (auto it = owner_pending_free_.begin(); it != owner_pending_free_.end();) {
        const cudaError_t q = cudaEventQuery(it->done);
        if (q == cudaErrorNotReady) break;  // (later ones ran after it)
        check_cuda(q, "owner pending free");
        for (char* b : it->blocks) free_block(b);
        cudaEventDestroy(it->done);
        it = owner_pending_free_.erase(it);
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::land_full_factors(
    const std::vector<const OwnerEntry*>& entries) {
    if (entries.empty()) return;
    const size_t D = sizeof(S), I = sizeof(int);
    cudaStream_t stream = Context::instance().stream();
    auto host = [&](const char* src, size_t count, const char* what) {
        std::vector<DataType> v(count);
        if (count > 0) check_cuda(cudaMemcpy(v.data(), src, count * D, cudaMemcpyDefault), what);
        return v;
    };
    struct Landed { const OwnerEntry* e; Box* box; OwnerLayout lay; size_t k, r, nt; };
    std::vector<Landed> landed;
    landed.reserve(entries.size());
    for (const OwnerEntry* ep : entries) {
        const OwnerEntry& e = *ep;
        Box* box = held_box(e.morton);
        if (box == nullptr || !is_local(e.morton)) {
            throw std::runtime_error("LevelEliminator: FULL state of box " + std::to_string(e.morton) + ", not local");
        }
        const size_t k = static_cast<size_t>(e.k), r = static_cast<size_t>(e.r), nt = static_cast<size_t>(e.ntot);
        const OwnerLayout lay(nt, r, k, static_cast<size_t>(e.nx), true);
        // T stays on the host (the matvec reads it there)
        box->interpolation_matrix.set_owned(e.k, e.r, host(e.bulk + lay.t, k * r, "landed T"), MatrixStorage<DataType>::FULL);
        box->X_RR = MatrixStorage<DataType>{};
        box->X_SR = MatrixStorage<DataType>{};
        box->X_RS_entry = MatrixStorage<DataType>{};
        box->X_RR_pivots.clear();
        factor_bytes_ += static_cast<double>(2 * align_up(k * r * D) + align_up(nt * r * D) + align_up(r * r * D) +
                                             align_up(r * I));
        landed.push_back(Landed{ep, box, lay, k, r, nt});
    }
    // ---- into the kept store, in groups whose X_NR fits the compressor's
    // staging (as keep_solve_factors groups a wave), straight from the
    // payloads; a group the keep refuses goes to the host with the rest
    std::vector<char> kept(landed.size(), 0);
    size_t at = 0;
    while (at < landed.size() && keep_solve_ && !keep_failed_) {
        std::vector<KeptSource> srcs;
        size_t raw = 0;
        const size_t first = at;
        while (at < landed.size()) {
            const Landed& L = landed[at];
            const size_t bytes = L.nt * L.r * D;
            if (compress_factors_ >= 1 && !srcs.empty() && raw + bytes > kXnrStagingBudget) break;
            KeptSource s;
            s.morton = L.e->morton;
            s.st = nullptr;  // (its X_RR_full came with the payload: no rebuild from this LU)
            s.t = L.e->bulk + L.lay.t;
            s.xsr = L.e->bulk + L.lay.xsr;
            s.xnr = L.e->bulk;
            s.lu = L.e->bulk + L.lay.lu;
            s.piv = L.e->bulk + L.lay.piv;
            s.k = L.e->k;
            s.r = L.e->r;
            s.ntot = L.e->ntot;
            s.pledge = false;
            srcs.push_back(s);
            raw += bytes;
            ++at;
        }
        if (!keep_sources_or_halves(srcs, stream)) break;  // (keep_failed_ is set: the rest go to the host)
        for (size_t i = first; i < at; ++i) kept[i] = 1;
    }
    for (size_t i = 0; i < landed.size(); ++i) {
        const Landed& L = landed[i];
        const OwnerEntry& e = *L.e;
        Box* box = L.box;
        const size_t k = L.k, r = L.r, nt = L.nt;
        const OwnerLayout& lay = L.lay;
        if (kept[i]) {
            skipped_factors_.push_back(box);
            if (nt > 0 && skip_xnr_) skipped_xnr_.push_back(box);
            box->X_NR = MatrixStorage<DataType>{};
        } else {
            box->X_RR.set_owned(e.r, e.r, host(e.bulk + lay.lu, r * r, "landed LU"), MatrixStorage<DataType>::LU_FACTORED);
            box->X_RR_pivots.resize(r);
            if (r > 0) check_cuda(cudaMemcpy(box->X_RR_pivots.data(), e.bulk + lay.piv, r * I, cudaMemcpyDefault), "landed pivots");
            box->X_SR.set_owned(e.k, e.r, host(e.bulk + lay.xsr, k * r, "landed X_SR"), MatrixStorage<DataType>::FULL);
            box->X_NR = MatrixStorage<DataType>{};
            if (nt > 0) box->X_NR.set_owned(e.ntot, e.r, host(e.bulk, nt * r, "landed X_NR"), MatrixStorage<DataType>::FULL);
        }
        box->X_RR_full.set_owned(e.r, e.r, host(e.bulk + lay.xrr, r * r, "landed X_RR_full"), MatrixStorage<DataType>::FULL);
        // the box's Schur block after its elimination: this rank's copy takes
        // the later sources' updates, the transition reads it
        {
            BoxState& st = state_of(e.morton);
            S* blk = heap_.alloc_resident<S>(std::max<size_t>(k * k * D, 1), ResidentKind::Near);
            if (k > 0) check_cuda(cudaMemcpyAsync(blk, e.bulk + lay.schur, k * k * D, cudaMemcpyDefault, stream), "landed Schur");
            if (st.schur.ptr != nullptr) heap_.free(st.schur.ptr);  // (later launches are ordered after its reads)
            st.schur = DeviceMatrix{blk, e.k, e.k};
        }
    }
}

// Received owner payloads installed (2026-10-11: split out of owner_exchange):
// their boxes' headers, the FULL factors of this rank's boxes, the skeleton
// marks and the generators' updates, wave by wave in the global order, then
// their receive buffers released.  owner_exchange installs a colour's
// remaining payloads at its end; owner_install_arrived those that have fully
// arrived at every wave boundary of the colour, so a colour's payloads are not
// all held to its end (~18 GB a rank for the green colour of CA3 at 8x).
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::install_owner_messages(
    int colour, const std::vector<OwnerWaveRecv*>& msgs, double t_wait, size_t sent_to, size_t from_ranks) {
    using hclock = std::chrono::high_resolution_clock;
    auto& stats = eliminator_stats();
    auto& timers = transport_timers();
    auto seconds_since = [](hclock::time_point t) { return std::chrono::duration<double>(hclock::now() - t).count(); };
    cudaStream_t stream = Context::instance().stream();
    if (msgs.empty()) return;
    // the colour's messages in arrival order (the bulk pointers are theirs)
    uint64_t bytes_in_total = 0;
    for (OwnerWaveRecv* r : msgs) {
        bytes_in_total += r->bytes;
        timers.bytes_received += static_cast<double>(r->header.size() * sizeof(int64_t) + r->bytes);
    }
    if (env::trace("exchange")) {
        const auto ob = owner_sent_bytes_.find(colour);
        std::printf("[exchange] rank %d level %d colour %d: wait %.3f s, out %.3f GB to %zu ranks, in %.3f GB in %zu wave messages from %zu ranks\n",
                    tree_->mpi_rank, level_index_, colour, t_wait, (ob == owner_sent_bytes_.end() ? 0 : ob->second) / 1e9, sent_to,
                    bytes_in_total / 1e9, msgs.size(), from_ranks);
        std::fflush(stdout);
    }
    // ---- the headers: the installed boxes' metadata
    const auto t_parse = hclock::now();
    std::vector<OwnerEntry> entries;
    for (size_t i = 0; i < msgs.size(); ++i) {
        const OwnerWaveRecv& msg = *msgs[i];
        const std::vector<int64_t>& h = msg.header;
        size_t at = 3;  // (after the colour, the wave and the bytes)
        auto next = [&]() -> int64_t {
            if (at >= h.size()) throw std::runtime_error("LevelEliminator::owner_exchange: truncated header");
            return h[at++];
        };
        auto list = [&](std::vector<int64_t>& v) {
            v.resize(static_cast<size_t>(next()));
            for (auto& x : v) x = next();
        };
        while (at < h.size()) {
            OwnerEntry e;
            e.morton = next();
            e.seq = static_cast<int32_t>(next());
            e.kind = static_cast<OwnerPayload>(next());
            if (e.kind == OwnerPayload::SKELETON) {
                list(e.skeleton);
                entries.push_back(std::move(e));
                continue;
            }
            e.r = static_cast<int>(next());
            e.k = static_cast<int>(next());
            e.ntot = static_cast<int>(next());
            e.nx = static_cast<int>(next());
            const int64_t off = next();
            list(e.one_hop);
            list(e.skeleton);
            std::vector<int64_t> redundant;
            list(redundant);
            if (e.r > 0) {
                e.counts.resize(e.one_hop.size());
                for (auto& c : e.counts) c = next();
                if (e.kind == OwnerPayload::FULL) {
                    e.use_full_set.resize(e.one_hop.size());
                    for (auto& u : e.use_full_set) u = next();
                }
                const OwnerLayout lay(static_cast<size_t>(e.ntot), static_cast<size_t>(e.r), static_cast<size_t>(e.k),
                                      static_cast<size_t>(e.nx), e.kind == OwnerPayload::FULL);
                if (off < 0 || static_cast<uint64_t>(off) + lay.total > msg.bytes) {
                    throw std::runtime_error("LevelEliminator::owner_exchange: a box's state outside its wave message");
                }
                e.bulk = msg.buf + static_cast<size_t>(off);
            }
            Box* b = held_box(e.morton);
            if (b == nullptr) {
                throw std::runtime_error("LevelEliminator::owner_exchange: rank " + std::to_string(tree_->mpi_rank) +
                                         " does not hold box " + std::to_string(e.morton));
            }
            if (is_writable(e.morton)) {
                throw std::runtime_error("LevelEliminator::owner_exchange: the state of box " + std::to_string(e.morton) +
                                         ", eliminated here, arrived from rank " + std::to_string(msg.source));
            }
            if (e.seq != msg.seq) throw std::runtime_error("LevelEliminator::owner_exchange: a box of another wave in a wave message");
            if (b->one_hop != e.one_hop) throw std::runtime_error("LevelEliminator::owner_exchange: one-hop lists differ");
            b->skeleton_indices.assign(e.skeleton.begin(), e.skeleton.end());
            b->redundant_indices.assign(redundant.begin(), redundant.end());
            b->deferred_xnn_neighbor_point_counts.assign(e.one_hop.size(), 0);
            if (e.r > 0) b->deferred_xnn_neighbor_point_counts = e.counts;
            if (e.kind == OwnerPayload::FULL && e.r > 0) b->use_full_set = e.use_full_set;
            if (is_local(e.morton)) {
                level_.solve_neighbor_size[static_cast<size_t>(e.morton - level_.local_morton_start)] =
                    b->deferred_xnn_neighbor_point_counts;
            }
            entries.push_back(std::move(e));
        }
    }
    // skeleton-only boxes: their host record (a ghost box, or an assisting box)
    for (OwnerEntry& e : entries) {
        if (e.kind != OwnerPayload::SKELETON) continue;
        if (Box* b = held_box(e.morton)) {
            b->skeleton_indices.assign(e.skeleton.begin(), e.skeleton.end());
        } else {
            auto it = level_.assisting_box_points_for_kernel_evaluation.find(e.morton);
            if (it == level_.assisting_box_points_for_kernel_evaluation.end()) {
                throw std::runtime_error("LevelEliminator::owner_exchange: skeleton of box " + std::to_string(e.morton) +
                                         ", unknown here");
            }
            level_.assisting_boxes[static_cast<size_t>(it->second)].skel_indices.assign(e.skeleton.begin(), e.skeleton.end());
        }
    }
    for (const OwnerEntry& e : entries) {
        level_.elimination_wave[e.morton] = e.seq;
        level_.eliminated_boxes.insert(e.morton);
    }
    if (env::trace("exchange")) {
        size_t full = 0, compact = 0, skel = 0, expected = 0, missing = 0;
        std::unordered_set<int64_t> got;
        for (const OwnerEntry& e : entries) {
            got.insert(e.morton);
            if (e.kind == OwnerPayload::FULL) ++full;
            else if (e.kind == OwnerPayload::COMPACT) ++compact;
            else ++skel;
        }
        for (const auto& [m, seq] : owner_->seq_of_box) {
            if (seq / std::max(owner_->num_waves, 1) != colour || is_writable(m) || held_box(m) == nullptr) continue;
            ++expected;
            if (got.count(m) == 0) ++missing;
        }
        std::printf("[exchange] rank %d level %d colour %d: entries %zu (full %zu, compact %zu, skeleton %zu); held boxes of the colour %zu, without a payload %zu; arena fallbacks %lld\n",
                    tree_->mpi_rank, level_index_, colour, entries.size(), full, compact, skel, expected, missing,
                    static_cast<long long>(stats.exchange_fallbacks));
        std::fflush(stdout);
    }
    timers.deserialize += seconds_since(t_parse);
    // ---- the device updates, wave by wave in the global order: a wave's
    //      skeletons (its full-rank and skeleton-only boxes), then the
    //      updates from its generators (which mark their own skeletons)
    const auto t_apply = hclock::now();
    std::unordered_map<int64_t, Incoming> incoming;
    std::map<int32_t, std::vector<size_t>> by_seq;
    for (size_t i = 0; i < entries.size(); ++i) {
        const OwnerEntry& e = entries[i];
        by_seq[e.seq].push_back(i);
        if (e.kind != OwnerPayload::SKELETON && e.r > 0) incoming.emplace(e.morton, Incoming{e.r, e.k, e.ntot, e.nx, e.bulk});
    }
    auto skeleton_of = [&](const OwnerEntry& e) -> const std::vector<int64_t>* {
        if (Box* b = held_box(e.morton)) return &b->skeleton_indices;
        return &level_.assisting_boxes[static_cast<size_t>(level_.assisting_box_points_for_kernel_evaluation.at(e.morton))]
                    .skel_indices;
    };
    for (const auto& [seq, idx] : by_seq) {
        // the FULL state of this rank's boxes first: their Schur blocks take
        // the wave's and later waves' updates
        {
            std::vector<const OwnerEntry*> full;
            for (size_t i : idx) {
                const OwnerEntry& e = entries[i];
                if (e.kind == OwnerPayload::FULL && e.r > 0) full.push_back(&e);
            }
            land_full_factors(full);  // (one batch: their factors are compressed together)
        }
        std::vector<std::pair<BoxState*, const std::vector<int64_t>*>> marks;
        std::vector<int64_t> gens;
        for (size_t i : idx) {
            const OwnerEntry& e = entries[i];
            if (e.kind != OwnerPayload::SKELETON && e.r > 0) gens.push_back(e.morton);
            else marks.emplace_back(&state_of(e.morton), skeleton_of(e));
        }
        set_remote_skeletons(marks);
        if (!gens.empty()) apply_remote_generators(gens, &incoming);
    }
    // (the installed generators stay resident for the level: no host copies
    // to restore them from; the receive buffers return to the arena once
    // the installs ran, as the receives of later waves may land in them)
    for (const OwnerEntry& e : entries) {
        if (e.kind == OwnerPayload::SKELETON || e.r == 0) continue;
        Box* b = held_box(e.morton);
        std::vector<DataType>().swap(b->deferred_xnn_temp2);
    }
    {
        OwnerPendingFree pf;
        const std::unordered_set<const OwnerWaveRecv*> here(msgs.begin(), msgs.end());
        for (auto it = owner_recvs_.begin(); it != owner_recvs_.end();) {
            if (here.count(&*it) == 0) { ++it; continue; }
            if (it->buf != nullptr) pf.blocks.push_back(it->buf);
            it = owner_recvs_.erase(it);
        }
        if (!pf.blocks.empty()) {
            check_cuda(cudaEventCreateWithFlags(&pf.done, cudaEventDisableTiming), "cudaEventCreate");
            check_cuda(cudaEventRecord(pf.done, stream), "cudaEventRecord");
            owner_pending_free_.push_back(std::move(pf));
        }
    }
    stats.remote += seconds_since(t_apply);
}

// At a wave boundary of an owner level: this colour's payloads whose bulk
// has fully arrived, installed now (see install_owner_messages).
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::owner_install_arrived(int colour) {
    std::vector<OwnerWaveRecv*> msgs;
    for (OwnerWaveRecv& r : owner_recvs_) {
        if (r.colour != colour) continue;
        if (!r.reqs.empty()) {
            int ok = 0;
            if (MPI_Testall(static_cast<int>(r.reqs.size()), r.reqs.data(), &ok, MPI_STATUSES_IGNORE) != MPI_SUCCESS) {
                throw std::runtime_error("LevelEliminator::owner_install_arrived: bulk Testall failed");
            }
            if (!ok) continue;
        }
        msgs.push_back(&r);
    }
    if (msgs.empty()) return;
    ++eliminator_stats().owner_early_installs;
    install_owner_messages(colour, msgs, 0.0, 0, 0);
}

template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::owner_exchange(int colour) {
    using hclock = std::chrono::high_resolution_clock;
    const auto t_start = hclock::now();
    auto& stats = eliminator_stats();
    auto& timers = transport_timers();
    auto seconds_since = [](hclock::time_point t) { return std::chrono::duration<double>(hclock::now() - t).count(); };
    if (owner_ == nullptr) throw std::runtime_error("LevelEliminator::owner_exchange: no plan");
    const OwnerPlanView::Colour& col = owner_->colours.at(static_cast<size_t>(colour));
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    MPI_Comm comm = tree_->comm;
    auto check_mpi = [](int err, const char* what) {
        if (err != MPI_SUCCESS) throw std::runtime_error(std::string("LevelEliminator::owner_exchange: ") + what + " failed");
    };
    constexpr int kTagMark = 735;
    const std::vector<int>& send_peers = col.send_peers;
    const std::vector<int>& recv_peers = col.recv_peers;
    const size_t ns = send_peers.size();
    // ---- this colour's remaining payloads out (their copies are done once
    //      the stream is), then every recipient's colour-end marker: how
    //      many wave messages it got from this rank this colour
    check_cuda(cudaStreamSynchronize(stream), "owner exchange: payloads");
    owner_progress();
    for (const auto& [key, count] : owner_sent_count_) {
        if (key.second == colour && !std::binary_search(send_peers.begin(), send_peers.end(), key.first)) {
            throw std::runtime_error("LevelEliminator::owner_exchange: a payload for rank " + std::to_string(key.first) +
                                     ", not a recipient of colour " + std::to_string(colour));
        }
    }
    std::vector<int64_t> marks_out(2 * ns, 0);
    std::vector<MPI_Request> mark_reqs(ns, MPI_REQUEST_NULL);
    for (size_t i = 0; i < ns; ++i) {
        marks_out[2 * i] = colour;
        const auto it = owner_sent_count_.find({send_peers[i], colour});
        marks_out[2 * i + 1] = it == owner_sent_count_.end() ? 0 : it->second;
        check_mpi(MPI_Isend(&marks_out[2 * i], 2, MPI_INT64_T, send_peers[i], kTagMark, comm, &mark_reqs[i]), "marker Isend");
    }
    // ---- the colour's arrivals: every receive peer's marker, as many wave
    //      messages, their bulk landed
    const auto t_comm = hclock::now();
    while (true) {
        owner_progress();
        bool done = true;
        for (int peer : recv_peers) {
            const auto mk = owner_marks_.find({peer, colour});
            if (mk == owner_marks_.end()) { done = false; break; }
            const auto rc = owner_recv_count_.find({peer, colour});
            const int got = rc == owner_recv_count_.end() ? 0 : rc->second;
            if (got < mk->second) { done = false; break; }
            if (got > mk->second) {
                throw std::runtime_error("LevelEliminator::owner_exchange: more wave messages from rank " + std::to_string(peer) +
                                         " than its marker says");
            }
        }
        for (auto it = owner_recvs_.begin(); done && it != owner_recvs_.end(); ++it) {
            if (it->colour != colour || it->reqs.empty()) continue;
            int ok = 0;
            check_mpi(MPI_Testall(static_cast<int>(it->reqs.size()), it->reqs.data(), &ok, MPI_STATUSES_IGNORE), "bulk Testall");
            if (!ok) done = false;
        }
        if (done) break;
    }
    check_mpi(MPI_Waitall(static_cast<int>(mark_reqs.size()), mark_reqs.data(), MPI_STATUSES_IGNORE), "marker Waitall");
    const double t_wait = seconds_since(t_comm);
    timers.payload += t_wait;
    for (const auto& [key, count] : owner_recv_count_) {
        if (key.second == colour && !std::binary_search(recv_peers.begin(), recv_peers.end(), key.first)) {
            throw std::runtime_error("LevelEliminator::owner_exchange: a payload from rank " + std::to_string(key.first) +
                                     ", not a sender of colour " + std::to_string(colour));
        }
    }
    // the colour's messages not yet installed, in arrival order
    std::vector<OwnerWaveRecv*> msgs;
    for (OwnerWaveRecv& r : owner_recvs_) {
        if (r.colour == colour) msgs.push_back(&r);
    }
    const auto t_apply = hclock::now();
    install_owner_messages(colour, msgs, t_wait, ns, recv_peers.size());
    for (auto it = owner_marks_.begin(); it != owner_marks_.end();) {
        it = it->first.second == colour ? owner_marks_.erase(it) : std::next(it);
    }
    for (auto it = owner_recv_count_.begin(); it != owner_recv_count_.end();) {
        it = it->first.second == colour ? owner_recv_count_.erase(it) : std::next(it);
    }
    ++owner_seq_;
    stats.exchange += std::chrono::duration<double>(t_apply - t_start).count();
    stats.device_exchange = true;
    if (heap_.peak() > stats.heap_peak) {
        stats.heap_peak = heap_.peak();
        stats.heap_peak_resident = heap_.resident_at_peak();
        stats.heap_peak_transient = heap_.transient_at_peak();
    }
}

// ---------------------------------------------------------------------------
// Device transition: the parent level's Schur and near-field blocks, from
// this level's blocks, fill sources and kernel (build_parent_level_
// interactions of the host).  A parent block is the 8 x 8 array of child
// blocks C(i, j) (skeleton x skeleton):
//   i == j            the child's Schur block (kernel if it has none)
//   one-hop           the stored near-field block
//   two-hop, sources  K - sum_E temp2_E[i] X_RR_full_E temp2_E[j]^T, in the
//                     host's source order (lazy far fill)
//   otherwise         K
// The parent Schur block holds C(i, j) of the parent's own children; the
// near-field block of parents P < Q holds C(i, j)^T for i in P, j in Q (rows
// = Q's points).  The lazy fill of C^T is P_{E,j}^T G_{E,i} with
// G_{E,x} = temp2_E[x rows]^T (columns of the resident temp2^T) and
// P_{E,j} = X_RR_full_E G_{E,j}.  Parents run in chunks that bound the
// memory of G and P.  A tree of occupied boxes (occupancy) builds only the
// occupied parents (num_children > 0) and their occupied neighbors
// (one_hop), and skips the pairs with an empty child (no skeleton), as the
// host's build_parent_level_interactions_unstructured does.
// ---------------------------------------------------------------------------
template<typename CoordType, typename DataType, typename KernelType>
std::unique_ptr<DeviceLevelBlocks<DataType>> LevelEliminator<CoordType, DataType, KernelType>::build_parent(
    std::vector<Box>& parents) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    auto& stats = eliminator_stats();
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    magma_queue_t queue = ctx.queue();
    if (!can_build_parent()) throw std::runtime_error("LevelEliminator::build_parent: needs the resident fill sources");
    // Released fill sources still on the device are read in place, unless
    // keeping them would leave less than half the heap for the parents'
    // blocks and the chunk buffers (the parents' blocks are at most the size
    // of this level's).  None is reclaimed during the transition: its plan
    // holds their pointers.
    // (the solve factors kept on the device are not counted)
    while (static_cast<double>(heap_.used()) - kept_bytes_total() > static_cast<double>(heap_.capacity() / 2) &&
           reclaim_released_source()) {
    }
    released_sources_.clear();
#ifdef H2_HAVE_NVCOMP
    // H2_compress_factors: whatever the waves did not read, the transition
    // does - a child's Schur block for a diagonal pair, its near blocks for
    // the rest - and it reads them without testing the pointer, falling back
    // to the kernel instead (that was a 1.2e-2 error in logdet, with no crash
    // and no guard).  This is the remainder after the waves (50 of 863 blocks
    // at level 3 of 96^3), so one batch here costs nothing and settles which
    // transition path reads what.
    if (!packed_schur_.empty() || !packed_edges_.empty()) {
        std::vector<AdoptedBlock> want;
        for (const auto& [morton, pb] : packed_schur_) {
            if (!is_interest(morton)) continue;
            BoxState& st = state_of(morton);
            if (st.schur.ptr == nullptr) {
                want.push_back(AdoptedBlock{true, static_cast<uint64_t>(morton), &st.schur});
            }
        }
        for (const auto& [key, pb] : packed_edges_) {
            auto it = edges_.find(key);
            if (it != edges_.end() && it->second.ptr == nullptr) {
                want.push_back(AdoptedBlock{false, key, &it->second});
            }
        }
        materialize_adopted(want, Context::instance().stream(), MatBuildParent);
    }
#endif
    const int dim = level_.dimension;
    const int nc = morton::children_per_box(dim);
    const size_t np = parents.size();
    if (np * static_cast<size_t>(nc) != level_.local_boxes.size()) {
        throw std::runtime_error("LevelEliminator::build_parent: parent/child count mismatch");
    }
    const size_t D = sizeof(S);
    const int64_t parent_start = parents.front().morton_index;
    const uint32_t grid = 1u << (level_.level - 1);
    auto is_local_parent = [&](int64_t pm) { return pm >= parent_start && pm < parent_start + static_cast<int64_t>(np); };
    auto built = [&](size_t p) { return !occupancy_ || parents[p].num_children > 0; };
    // with occupancy, an empty child: a local box without points, or a box
    // of another rank without point data (its record is the zero-sized slot
    // of add_empty_parent_transition_assisting_slots)
    auto empty_child = [&](int64_t cm) {
        if (!occupancy_) return false;
        if (is_local(cm)) return state_of(cm).n == 0;
        if (remote_index_.count(cm) != 0) return false;
        auto it = level_.assisting_box_points_for_kernel_evaluation.find(cm);
        if (it != level_.assisting_box_points_for_kernel_evaluation.end() &&
            !level_.assisting_boxes[static_cast<size_t>(it->second)].indices.empty()) {
            throw std::runtime_error("LevelEliminator::build_parent: child " + std::to_string(cm) +
                                     " has point data but no device state");
        }
        return true;
    };

    // Parents: the local ones and their neighbors on other ranks (whose
    // children are assisting boxes, all eliminated by now), with the parent
    // point offsets of their children (skeletons, in child order).
    struct Parent {
        int64_t morton = -1;
        std::vector<int64_t> child;
        std::vector<int> offset;  // nc + 1
        int n = 0;
    };
    std::unordered_map<int64_t, Parent> parent_of;
    auto parent = [&](int64_t pm) -> Parent& {
        auto it = parent_of.find(pm);
        if (it != parent_of.end()) return it->second;
        Parent& pa = parent_of[pm];
        pa.morton = pm;
        pa.offset.assign(static_cast<size_t>(nc) + 1, 0);
        for (int c = 0; c < nc; ++c) {
            const int64_t cm = pm * nc + c;
            pa.child.push_back(cm);
            if (empty_child(cm)) {
                pa.offset[static_cast<size_t>(c) + 1] = pa.offset[static_cast<size_t>(c)];
                continue;
            }
            const BoxState& cs = state_of(cm);
            if (!cs.eliminated) throw std::runtime_error("LevelEliminator::build_parent: child " + std::to_string(cm) + " not eliminated");
            pa.offset[static_cast<size_t>(c) + 1] = pa.offset[static_cast<size_t>(c)] + static_cast<int>(cs.skeleton.size());
        }
        pa.n = pa.offset[static_cast<size_t>(nc)];
        return pa;
    };
    for (size_t p = 0; p < np; ++p) {
        if (!built(p)) continue;
        if (parent(parents[p].morton_index).n != parents[p].num_points) {
            throw std::runtime_error("LevelEliminator::build_parent: parent point count mismatch");
        }
    }
    // host skeleton (positions) of a child, for the fill-row lookups
    auto host_skeleton = [&](int64_t cm) -> const std::vector<int64_t>* { return ring_info(cm).skeleton; };
    auto host_full = [&](int64_t cm) -> int64_t { return ring_info(cm).full; };
    // The lazy far fill's sources of a child pair (i, j): the eliminated
    // one-hop neighbors of i next to j, in the host's order (wave, Morton;
    // collect_lazy_far_sources).  Listed once per local child here, with
    // their coordinates, so a pair filters a short list instead of decoding
    // and looking up each of its 26 neighbors (the leaf plans hundreds of
    // thousands of two-hop pairs).  A neighbor without a recorded wave is
    // left to collect_lazy_far_sources, which reports it.  Each carries its
    // device state and a dense index of it (the local boxes, then the remote
    // ones): the chunks look up a source's G and P blocks by that index and
    // the direction from the source to the child (one of its 3^dim
    // neighbors) instead of hashing Morton pairs, a million times at the leaf.
    struct SourceCandidate {
        int32_t wave;
        int64_t morton;
        int64_t coords[3];
        bool known;
        const BoxState* state;  // (null: no device state; state_of reports it when the source is used)
        int id;
    };
    const size_t num_source_ids = states_.size() + remote_states_.size();
    constexpr int kDirections = 27;
    std::vector<std::vector<SourceCandidate>> candidates_of(level_.local_boxes.size());
    for (size_t b = 0; b < level_.local_boxes.size(); ++b) {
        auto& list = candidates_of[b];
        for (int64_t candidate : level_.local_boxes[b].one_hop) {
            if (!sketch_box_eliminated(level_, candidate)) continue;
            const auto wave_it = level_.elimination_wave.find(candidate);
            SourceCandidate c{wave_it != level_.elimination_wave.end() ? wave_it->second : 0, candidate, {0, 0, 0},
                              wave_it != level_.elimination_wave.end(), nullptr, -1};
            lazy_far_decode_coords(dim, candidate, c.coords);
            if (is_local(candidate)) {
                c.id = static_cast<int>(candidate - level_.local_morton_start);
                c.state = &states_[static_cast<size_t>(c.id)];
            } else if (auto rit = remote_index_.find(candidate); rit != remote_index_.end()) {
                c.id = static_cast<int>(states_.size() + rit->second);
                c.state = &remote_states_[rit->second];
            }
            list.push_back(c);
        }
        std::sort(list.begin(), list.end(), [](const SourceCandidate& a, const SourceCandidate& b) {
            return a.wave != b.wave ? a.wave < b.wave : a.morton < b.morton;
        });
    }
    auto sources_of = [&](int64_t im, int64_t jm, const int64_t cj3[3], std::vector<const SourceCandidate*>& out) {
        out.clear();
        for (const SourceCandidate& c : candidates_of[static_cast<size_t>(im - level_.local_morton_start)]) {
            int64_t chebyshev = 0;
            for (int d = 0; d < dim; ++d) chebyshev = std::max(chebyshev, std::abs(c.coords[d] - cj3[d]));
            if (chebyshev != 1) continue;
            if (!c.known) {
                std::vector<LazyFarSource> unused;
                collect_lazy_far_sources(level_, level_.find_local_box(im), jm, dim, unused);  // (throws)
                throw std::runtime_error("LevelEliminator::build_parent: fill source " + std::to_string(c.morton) +
                                         " has no elimination wave");
            }
            if (c.state == nullptr) state_of(c.morton);  // (throws)
            out.push_back(&c);
        }
    };
    // the direction from a source to one of its neighbors: 0 .. 3^dim - 1
    auto direction = [&](const int64_t from[3], const int64_t to[3]) {
        int slot = 0, stride = 1;
        for (int d = 0; d < dim; ++d) {
            const int64_t delta = to[d] - from[d];
            if (delta < -1 || delta > 1) throw std::runtime_error("LevelEliminator::build_parent: a fill source is not next to its endpoint");
            slot += static_cast<int>(delta + 1) * stride;
            stride *= 3;
        }
        return slot;
    };
    std::vector<int> g_index, p_index;  // per chunk: (source id, direction) -> its G, its P (-1: none yet)
    std::vector<int> g4_index;          // ... and (a nonsymmetric matrix) its temp4 columns

    // Parent blocks: the Schur block of each local parent, and a near-field
    // block for each neighbor parent (a local pair once, from its lower
    // parent; a pair with a parent of another rank always: this rank's copy).
    auto out = std::make_unique<Blocks>();
    out->directed = nonsym_;
    // (a nonsymmetric matrix: K(P, Q) and K(Q, P) of a pair, rows_q telling which)
    struct Pair { int64_t pm, qm; DeviceMatrix* block; bool rows_q = false; };  // qm == pm: Schur block
    std::vector<std::vector<Pair>> pairs_of(np);
    for (size_t p = 0; p < np; ++p) {
        if (!built(p) || (occupancy_ && parents[p].num_points == 0)) continue;
        const int64_t pm = parents[p].morton_index;
        const int n = parents[p].num_points;
        // the shape now, the block when its chunk is planned: every block is
        // written by exactly one chunk, so allocating here would hold the
        // whole set (9.47 GB at 192^3, 74% of the level's peak) for the sake
        // of one chunk's worth of it
        DeviceMatrix& sm = out->schur[pm];
        sm = DeviceMatrix{nullptr, n, n};
        pairs_of[p].push_back(Pair{pm, pm, &sm});
    }
    for (size_t p = 0; p < np; ++p) {
        if (!built(p) || (occupancy_ && parents[p].num_points == 0)) continue;
        const int64_t pm = parents[p].morton_index;
        std::vector<int64_t> neighbors;
        if (occupancy_) {
            neighbors = parents[p].one_hop;  // the occupied ones
        } else {
            for (uint64_t qu : morton::neighbors_nd(dim, static_cast<uint64_t>(pm), grid)) {
                neighbors.push_back(static_cast<int64_t>(qu));
            }
        }
        for (const int64_t qm : neighbors) {
            if (is_local_parent(qm) && qm <= pm) continue;
            const Parent& q = parent(qm);
            if (occupancy_ && q.n == 0) continue;
            if (nonsym_) {
                const int pn = static_cast<int>(parents[p].num_points);
                DeviceMatrix& pq = out->edges[dkey(pm, qm)];
                pq = DeviceMatrix{nullptr, pn, q.n};
                pairs_of[p].push_back(Pair{pm, qm, &pq, false});
                DeviceMatrix& qp = out->edges[dkey(qm, pm)];
                qp = DeviceMatrix{nullptr, q.n, pn};
                pairs_of[p].push_back(Pair{pm, qm, &qp, true});
                continue;
            }
            const int64_t lo = std::min(pm, qm), hi = std::max(pm, qm);
            const int rows = hi == qm ? q.n : parents[p].num_points, cols = hi == qm ? parents[p].num_points : q.n;
            DeviceMatrix& em = out->edges[edge_key(lo, hi)];
            em = DeviceMatrix{nullptr, rows, cols};
            pairs_of[p].push_back(Pair{pm, qm, &em});
        }
    }

    // Fill sources freed during the waves (multi-rank levels) are restored
    // once, in the chunk of their first reader, as [G_all | X_RR_full]
    // (ld r): G_all holds the columns G_{E,x} = temp2_E(rows of x's
    // skeleton, :)^T of every neighbor x.  The host gathers the rows (and
    // X_RR_full^T) column by column into pinned slices, uploaded as they
    // fill, and the device transposes them into the source's block.  Each
    // block is freed after the chunk of its last reader (the parent of the
    // source's last local neighbor).
    struct Restored {
        const S* gall = nullptr;
        const S* xrr = nullptr;
        const S* g4all = nullptr;  // (nonsymmetric: temp4's columns at the same points as gall)
        std::unordered_map<int64_t, int> col_of;  // neighbor -> first column of its G
        char* block = nullptr;
        char* xrr_block = nullptr;   // (H2 host memory) its X_RR_full alone, restored from the host
        size_t last_parent = 0;
    };
    std::unordered_map<int64_t, Restored> restored;
    // X_RR_full of the sources that dropped it for their LU, rebuilt once a
    // chunk into blocks of its own (3 r x r: the product and its L and U),
    // freed with the chunk.  Freed at the transition's end before (H2 8x,
    // 2026-10-10): over the 384^3 leaf's 100+ chunks they came to 14 GB
    // and starved the chunk budget to 50 MB.
    std::vector<char*> xrr_exp_blocks;
    cudaEvent_t staged[2] = {nullptr, nullptr};
    int stage_turn = 0;
    auto last_reader = [&](const Box* sb) {
        size_t last = 0;
        for (int64_t m : sb->one_hop) {
            if (is_local(m)) last = std::max(last, static_cast<size_t>(m - level_.local_morton_start) / static_cast<size_t>(nc));
        }
        return last;
    };

    size_t p0 = 0;
    std::vector<const SourceCandidate*> sources;
    std::vector<int64_t> positions;
    int64_t fill_gemms = 0;
    // Chunks alternate between two metadata images and are not waited for:
    // the host plans chunk c + 1 while the device runs chunk c.
    std::vector<std::unique_ptr<StreamMarks>> chunk_marks;
    while (p0 < np) {
        const auto t_plan = clock::now();
        ++wave_stamp_;
        const bool second = chunk_marks.size() % 2 == 1;
        MetaBuilder& meta = second ? owner_meta_ : meta_;
        DeviceBuffer& meta_device = second ? owner_meta_device_ : meta_device_;
        meta.clear();
        std::vector<EvalItem> evals;
        std::vector<GatherItem> gathers, g_gathers;
        int max_eval_m = 0, max_eval_n = 0, max_gather_m = 0, max_gather_n = 0, max_g_m = 0, max_g_n = 0;
        // chunk budget for G, P and restored fill sources (the free space is
        // fragmented by the level's blocks)
        const size_t budget = std::min({size_t{2} << 30, (heap_.capacity() - heap_.used()) / 4, heap_.largest_free() / 2});
        size_t buffer_bytes = 0;
        // A fill source: its resident temp2^T and X_RR_full, or its restored
        // G columns (restored in an earlier chunk, or in this one: offsets in
        // this chunk's restore block).
        struct Source {
            const S* t2t = nullptr;   // resident temp2^T
            const S* gall = nullptr;  // restored G columns
            const S* xrr = nullptr;
            const S* t4 = nullptr;    // (nonsymmetric) resident temp4, and
            const S* g4all = nullptr; // its restored columns
            size_t gall_off = SIZE_MAX, xrr_off = SIZE_MAX, g4_off = SIZE_MAX;  // in this chunk's restore block
            const std::unordered_map<int64_t, int>* col_of = nullptr;
            int r = 0;
            Box* box = nullptr;
        };
        std::unordered_map<int64_t, Source> source_of;
        // the chunk's compressed fill sources, decompressed in one batch below
        std::vector<const void*> packed_src;
        std::vector<void*> packed_dst;
        std::vector<LuSplitItemT<S>> xrr_split;
        std::vector<RowSwapItemT<S>> xrr_swap;
        VBatch<S> xrr_gemm;
        int max_xrr_r = 0;
        // restored in this chunk, at increasing offsets of the staging buffer
        struct Pending { int64_t morton; Box* box; int r; size_t gall_off, xrr_off, end; size_t last_parent; std::vector<int64_t> rows; bool xrr_only = false; size_t g4_off = SIZE_MAX; };
        std::vector<Pending> pending;
        std::unordered_map<int64_t, std::unordered_map<int64_t, int>> pending_cols;
        size_t restore_bytes = 0;
        // the packed sources put back whole in this chunk (H2 8x, 2026-10-10):
        // counted against the budget like the host restores - at the 384^3
        // leaf every source is packed and they are the level's raw temp2^T,
        // 146 GB a rank, so a chunk sized by its parents alone filled the card
        size_t packed_restore_bytes = 0;
        size_t xrr_rebuild_bytes = 0;    // ... and the X_RR_full rebuilt from the LU for this chunk
        std::vector<int64_t> rows_of;
        auto source = [&](int64_t em) -> Source& {
            auto it = source_of.find(em);
            if (it != source_of.end()) return it->second;
            const BoxState& es = state_of(em);
            Source& s = source_of[em];
            s.r = es.r;
            s.box = es.remote ? level_.find_generator_box(em) : held_box(em);
            // X_RR_full: stored, or rebuilt from the LU this box kept
            const S* xrr_of_es = es.xrr_full;
            if (xrr_of_es == nullptr && xrr_rebuilt(es)) {
                const size_t rr = static_cast<size_t>(es.r) * static_cast<size_t>(es.r);
                char* blk = heap_.try_alloc_resident(3 * rr * sizeof(S), ResidentKind::Restore);
                if (blk == nullptr) {
                    throw std::runtime_error("LevelEliminator::build_parent: no room to rebuild X_RR_full of " +
                                             std::to_string(em));
                }
                xrr_exp_blocks.push_back(blk);
                xrr_rebuild_bytes += 3 * rr * sizeof(S);
                S* dst = reinterpret_cast<S*>(blk);
                S* l = dst + rr;
                S* u = l + rr;
                xrr_split.push_back(LuSplitItemT<S>{es.xrr_lu, es.r, l, u, es.r});
                xrr_gemm.entries.push_back({l, u, dst, es.r, es.r, es.r, es.r, es.r, es.r});
                xrr_swap.push_back(RowSwapItemT<S>{dst, es.r, es.r, es.xrr_piv, es.r});
                max_xrr_r = std::max(max_xrr_r, es.r);
                xrr_of_es = dst;
                ++eliminator_stats().xr_rebuilt;
            }
            if (s.box == nullptr) throw std::runtime_error("LevelEliminator::build_parent: fill source " + std::to_string(em) + " not found");
            if (es.temp2t != nullptr) {
                s.t2t = es.temp2t;
                s.xrr = xrr_of_es;
                if (nonsym_) {
                    if (es.temp4 == nullptr) {
                        throw std::runtime_error("LevelEliminator::build_parent: fill source " + std::to_string(em) +
                                                 " has its temp2^T but not its temp4");
                    }
                    s.t4 = es.temp4;
                }
                return s;
            }
            auto rit = restored.find(em);
            if (rit != restored.end()) {
                s.gall = rit->second.gall;
                s.g4all = rit->second.g4all;
                s.xrr = rit->second.xrr;
                s.col_of = &rit->second.col_of;
                if (rit->second.gall == nullptr) s.t2t = reinterpret_cast<const S*>(rit->second.block);
                return s;
            }
#ifdef H2_HAVE_NVCOMP
            // H2_compress_factors: stored compressed, so put it back whole
            // into a block of this chunk and read it as a resident one.  That
            // is both simpler than gathering its rows out of the slots and
            // faster than the host restore below (0.53 s for 2.40 GB at the
            // leaf, against decompression at 365 GB/s).
            if (packed_of(es) != nullptr && es.ntot > 0 && es.r > 0) {
                const size_t elems = static_cast<size_t>(es.ntot) * static_cast<size_t>(es.r);
                char* blk = heap_.try_alloc_resident(elems * sizeof(S), ResidentKind::Restore);
                if (blk != nullptr) {
                    // collected now, decompressed in one batch before the
                    // chunk's launches (the pointer is only read there)
                    if (collect_fill_source(es, reinterpret_cast<S*>(blk), packed_src, packed_dst)) {
                        packed_restore_bytes += elems * sizeof(S);
                        Restored& rs = restored[em];
                        rs.block = blk;
                        rs.gall = nullptr;
                        rs.xrr = nullptr;
                        rs.last_parent = p0;  // freed at the end of this chunk (p0 < p1)
                        s.t2t = reinterpret_cast<const S*>(blk);
                        s.xrr = xrr_of_es;
                        if (s.xrr == nullptr) {
                            // H2 host memory: its X_RR_full left the device and
                            // no LU rebuilds it (a remote generator; a box whose
                            // keep failed): the host copy, or the rescued host
                            // LU, through this chunk's restore block as a
                            // source with no temp2 rows.
                            if (!host_xrr_available(*s.box, s.r)) {
                                throw std::runtime_error("LevelEliminator::build_parent: X_RR_full of fill source " +
                                                         std::to_string(em) + " is neither on the device nor on the host");
                            }
                            Pending pd;
                            pd.morton = em;
                            pd.box = s.box;
                            pd.r = s.r;
                            pd.xrr_only = true;
                            pd.gall_off = restore_bytes;
                            pd.xrr_off = restore_bytes;
                            restore_bytes = align_up(pd.xrr_off + static_cast<size_t>(s.r) * s.r * D);
                            pd.end = restore_bytes;
                            pd.last_parent = p0;
                            s.xrr_off = pd.xrr_off;
                            pending.push_back(std::move(pd));
                            ++eliminator_stats().tr_xrr_host;
                        }
                        return s;
                    }
                    heap_.free(blk);
                }
            }
#endif

            // restore in this chunk: the skeleton rows of each neighbor with
            // point data here (the only ones a pair of this rank can read)
            Box* sb = s.box;
            if (sb->X_NR.rows != es.ntot || sb->X_NR.cols != s.r || !host_xrr_available(*sb, s.r) ||
                (nonsym_ && (sb->X_RN.rows != s.r || sb->X_RN.cols != es.ntot || sb->X_RN.lda != s.r))) {
#ifdef H2_HAVE_NVCOMP
                if (packed_of(es) != nullptr && es.ntot > 0 && es.r > 0 && !sb->X_NR.is_allocated()) {
                    // (no host copy under packed_restore_: the packed put-back
                    // above found no room, which is the chunking's problem)
                    throw std::runtime_error("LevelEliminator::build_parent: no room to put back fill source " +
                                             std::to_string(em) + " from the packed store (" +
                                             std::to_string(static_cast<size_t>(es.ntot) * es.r * sizeof(S)) +
                                             " bytes wanted; heap used " + std::to_string(heap_.used() >> 20) + " of " +
                                             std::to_string(heap_.capacity() >> 20) + " MB, largest free " +
                                             std::to_string(heap_.largest_free() >> 20) + " MB, chunk budget " +
                                             std::to_string(budget >> 20) + " MB)");
                }
#endif
                throw std::runtime_error("LevelEliminator::build_parent: host copy of fill source " +
                                         std::to_string(em) + " missing");
            }
            Pending pd;
            pd.morton = em;
            pd.box = sb;
            pd.r = s.r;
            auto& cols = pending_cols[em];
            for (int64_t xm : sb->one_hop) {
                if (!is_local(xm) && remote_index_.count(xm) == 0) continue;
                LazyFarEndpoint<CoordType> endpoint;
                endpoint.morton = xm;
                endpoint.full_size = host_full(xm);
                endpoint.skeleton = host_skeleton(xm);
                endpoint.wanted = endpoint.skeleton;
                endpoint.wanted_count = static_cast<int64_t>(endpoint.skeleton->size());
                int64_t slot_offset = 0, slot_count = 0;
                if (endpoint.wanted_count == 0 ||
                    !lazy_far_locate_endpoint_rows(sb, endpoint, slot_offset, slot_count, rows_of)) {
                    continue;
                }
                cols[xm] = static_cast<int>(pd.rows.size());
                for (int64_t q : rows_of) pd.rows.push_back(slot_offset + q);
            }
            pd.gall_off = restore_bytes;
            pd.xrr_off = restore_bytes + static_cast<size_t>(s.r) * pd.rows.size() * D;
            restore_bytes = align_up(pd.xrr_off + static_cast<size_t>(s.r) * s.r * D);
            if (nonsym_) {  // [G_all | X_RR_full | temp4 at the same points]
                pd.g4_off = pd.xrr_off + static_cast<size_t>(s.r) * s.r * D;
                restore_bytes = align_up(pd.g4_off + static_cast<size_t>(s.r) * pd.rows.size() * D);
            }
            pd.end = restore_bytes;
            pd.last_parent = last_reader(sb);
            s.gall_off = pd.gall_off;
            s.xrr_off = pd.xrr_off;
            s.g4_off = pd.g4_off;
            s.col_of = &cols;
            pending.push_back(std::move(pd));
            return s;
        };
        // G_{E,x} and P_{E,j}, deduplicated per chunk; a pointer is absolute,
        // or an offset in the chunk's buffers (G gathers, P) or in its
        // restore block
        struct GInfo { const S* ptr; int ld; size_t buffer = SIZE_MAX; size_t restore = SIZE_MAX; };
        std::vector<GInfo> g_infos;       // the chunk's G blocks, and
        std::vector<size_t> p_offsets;    // its P blocks (buffer offsets), by g_index / p_index
        g_index.assign(num_source_ids * kDirections, -1);
        p_index.assign(num_source_ids * kDirections, -1);
        if (nonsym_) g4_index.assign(num_source_ids * kDirections, -1);
        struct GPending { size_t offset; const S* src; int r; int k; IndexList cols; };
        std::vector<GPending> g_pending;
        struct PPending { size_t offset; int g; const S* xrr; size_t xrr_restore; int r; int k; };
        std::vector<PPending> p_pending;
        // out -= P_{E,j}^T G_{E,i} (rows j) or G_{E,i}^T P_{E,j} (rows i)
        struct FillTask { int g_i; int p_j; S* out; int ld; int m; int n; int r; bool rows_i; };
        std::vector<std::vector<FillTask>> fills;  // per target child pair, in source order

        // G_{E,x}: temp2_E rows of child x (its skeleton), as columns of temp2_E^T.
        // col_side (a nonsymmetric matrix): temp4_E's columns of x instead,
        // laid out as temp2^T is.
        auto g_for = [&](const SourceCandidate& e, int64_t xm, const int64_t x3[3], bool col_side = false) -> int {
            const size_t slot = static_cast<size_t>(e.id) * kDirections + static_cast<size_t>(direction(e.coords, x3));
            std::vector<int>& index = col_side ? g4_index : g_index;
            if (index[slot] >= 0) return index[slot];
            const int64_t em = e.morton;
            auto record = [&](const GInfo& info) {
                g_infos.push_back(info);
                index[slot] = static_cast<int>(g_infos.size()) - 1;
                return index[slot];
            };
            Source& s = source(em);
            BoxState& xs = state_of(xm);
            if (s.t2t == nullptr) {  // restored: contiguous columns
                auto cit = s.col_of->find(xm);
                if (cit == s.col_of->end()) {
                    throw std::runtime_error("LevelEliminator::build_parent: restored fill source without rows for its endpoint");
                }
                GInfo info;
                info.ld = s.r;
                const S* all = col_side ? s.g4all : s.gall;
                const size_t all_off = col_side ? s.g4_off : s.gall_off;
                if (s.gall != nullptr) {
                    info.ptr = all + static_cast<size_t>(cit->second) * s.r;
                } else {
                    info.ptr = nullptr;
                    info.restore = all_off + static_cast<size_t>(cit->second) * s.r * D;
                }
                return record(info);
            }
            const S* base = col_side ? s.t4 : s.t2t;
            LazyFarEndpoint<CoordType> endpoint;
            endpoint.morton = xm;
            endpoint.full_size = host_full(xm);
            endpoint.skeleton = host_skeleton(xm);
            endpoint.wanted = endpoint.skeleton;
            endpoint.wanted_count = static_cast<int64_t>(endpoint.skeleton->size());
            int64_t slot_offset = 0, slot_count = 0;
            if (!lazy_far_locate_endpoint_rows(s.box, endpoint, slot_offset, slot_count, positions)) {
                throw std::runtime_error("LevelEliminator::build_parent: fill source without rows for its endpoint");
            }
            GInfo info;
            info.ld = s.r;
            const int k = static_cast<int>(endpoint.skeleton->size());
            const bool by_skeleton = slot_count == endpoint.full_size && !xs.skeleton_identity;
            if (!by_skeleton) {  // contiguous columns
                info.ptr = base + static_cast<size_t>(slot_offset) * s.r;
            } else {
                IndexList cols = current_slots_in(xs, meta, wave_stamp_);
                cols.base = static_cast<int>(slot_offset);
                info.ptr = nullptr;
                info.buffer = buffer_bytes;
                g_pending.push_back(GPending{buffer_bytes, base, s.r, k, cols});
                buffer_bytes = align_up(buffer_bytes + static_cast<size_t>(s.r) * k * D);
            }
            return record(info);
        };
        auto p_for = [&](const SourceCandidate& e, int64_t jm, const int64_t j3[3]) -> int {
            const size_t slot = static_cast<size_t>(e.id) * kDirections + static_cast<size_t>(direction(e.coords, j3));
            if (p_index[slot] >= 0) return p_index[slot];
            // (a nonsymmetric matrix: P_{E,j} = X_RR_full_E temp4_E[:, j], the fill's column side)
            const int g = g_for(e, jm, j3, nonsym_);
            const Source& s = source(e.morton);
            const int k = static_cast<int>(state_of(jm).skeleton.size());
            p_offsets.push_back(buffer_bytes);
            p_index[slot] = static_cast<int>(p_offsets.size()) - 1;
            p_pending.push_back(PPending{buffer_bytes, g, s.xrr, s.xrr == nullptr ? s.xrr_off : SIZE_MAX, s.r, k});
            buffer_bytes = align_up(buffer_bytes + static_cast<size_t>(s.r) * k * D);
            return p_index[slot];
        };

        // plan the chunk's parents
        // Keyed by the map slot, not by the block's address: the blocks are
        // allocated and freed per chunk now, so an address is reused and
        // would both merge two blocks here and count a false revisit.
        std::unordered_map<const void*, size_t> chunk_blocks;  // target slot -> its bytes
        size_t chunk_block_bytes = 0;
        // ... and which key it is under, so the compressed form can be stored
        struct ChunkRef { bool schur; int64_t morton; uint64_t key; const DeviceMatrix* m; };
        std::vector<ChunkRef> chunk_refs;
        size_t p1 = p0;
        while (p1 < np && (p1 == p0 || buffer_bytes + 2 * restore_bytes + packed_restore_bytes + xrr_rebuild_bytes +
                                              chunk_block_bytes < budget)) {
            const size_t p = p1++;
            for (const Pair& pr : pairs_of[p]) {
                DeviceMatrix& M = *pr.block;
                const size_t m_bytes = static_cast<size_t>(M.rows) * M.cols * D;
                if (m_bytes > 0 && chunk_blocks.emplace(pr.block, m_bytes).second) {
                    if (M.ptr == nullptr) {
                        M.ptr = heap_.alloc_resident<S>(m_bytes, ResidentKind::Transition);
                    }
                    chunk_block_bytes += m_bytes;
                    const bool is_schur = pr.qm == pr.pm;
                    chunk_refs.push_back(ChunkRef{is_schur, pr.pm,
                                                  is_schur ? 0
                                                  : nonsym_ ? (pr.rows_q ? dkey(pr.qm, pr.pm) : dkey(pr.pm, pr.qm))
                                                            : edge_key(pr.pm, pr.qm),
                                                  &M});
                }
                const bool diag = pr.qm == pr.pm;
                const Parent& pa = parent(pr.pm);
                const Parent& qa = parent(pr.qm);
                // rows of M: the higher parent (Q for a local pair); a
                // nonsymmetric matrix's: the pair's orientation
                const bool rows_q = nonsym_ ? pr.rows_q : !diag && pr.qm > pr.pm;
                for (int ci = 0; ci < nc; ++ci) {
                    const int64_t im = pa.child[static_cast<size_t>(ci)];
                    for (int cj = 0; cj < nc; ++cj) {
                        const int64_t jm = qa.child[static_cast<size_t>(cj)];
                        // target view: rows R, columns S, at (row0, col0) of M
                        const int64_t Rm = rows_q ? jm : im, Sm = rows_q ? im : jm;
                        if (empty_child(Rm) || empty_child(Sm)) continue;
                        const int row0 = rows_q ? qa.offset[static_cast<size_t>(cj)] : pa.offset[static_cast<size_t>(ci)];
                        const int col0 = rows_q ? pa.offset[static_cast<size_t>(ci)] : qa.offset[static_cast<size_t>(cj)];
                        BoxState& rs = state_of(Rm);
                        BoxState& ss = state_of(Sm);
                        const int kr = static_cast<int>(rs.skeleton.size()), ks = static_cast<int>(ss.skeleton.size());
                        if (kr == 0 || ks == 0) continue;  // (an empty block)
                        S* target = M.ptr + row0 + static_cast<size_t>(col0) * M.rows;
                        auto eval = [&] {
                            evals.push_back(EvalItem{target, M.rows, kr, ks, current_slots_in(rs, meta, wave_stamp_),
                                                     current_slots_in(ss, meta, wave_stamp_)});
                            max_eval_m = std::max(max_eval_m, kr);
                            max_eval_n = std::max(max_eval_n, ks);
                        };
                        auto gather = [&](const S* src, int64_t rsd, int64_t csd) {
                            gathers.push_back(GatherItem{target, M.rows, kr, ks, src, rsd, csd, IndexList{}, IndexList{}});
                            max_gather_m = std::max(max_gather_m, kr);
                            max_gather_n = std::max(max_gather_n, ks);
                        };
                        if (im == jm) {
                            if (rs.schur.ptr != nullptr) {
                                if (rs.schur.rows != kr || rs.schur.cols != kr) {
                                    throw std::runtime_error("LevelEliminator::build_parent: child Schur block shape");
                                }
                                gather(rs.schur.ptr, 1, kr);
                            } else {
                                assert_not_packed(true, static_cast<uint64_t>(im), "build_parent's child Schur");
                                eval();
                            }
                            continue;
                        }
                        int64_t chebyshev = 0;
                        bool reaches_two = false;
                        int64_t ci3[3], cj3[3];
                        {
                            lazy_far_decode_coords(dim, im, ci3);
                            lazy_far_decode_coords(dim, jm, cj3);
                            for (int d = 0; d < dim; ++d) {
                                const int64_t dist = std::abs(ci3[d] - cj3[d]);
                                chebyshev = std::max(chebyshev, dist);
                                reaches_two = reaches_two || dist == 2;
                            }
                        }
                        if (chebyshev <= 1 && nonsym_) {
                            // the stored K(R, S) (rows = R's points), else the kernel
                            auto it = edges_.find(dkey(Rm, Sm));
                            if (it != edges_.end() && it->second.ptr == nullptr) {
                                assert_not_packed(false, dkey(Rm, Sm), "build_parent");
                            }
                            if (it == edges_.end()) {
                                eval();
                                continue;
                            }
                            const DeviceMatrix& E = it->second;
                            if (E.rows != kr || E.cols != ks) throw std::runtime_error("LevelEliminator::build_parent: child near block shape");
                            gather(E.ptr, 1, E.rows);
                            continue;
                        }
                        if (chebyshev <= 1) {
                            auto it = edges_.find(edge_key(im, jm));
                            if (it != edges_.end() && it->second.ptr == nullptr) {
                                assert_not_packed(false, edge_key(im, jm), "build_parent");
                            }
                            if (it == edges_.end()) {
                                // a pair no update reached (e.g. of full-rank
                                // boxes, or with a box of another rank) keeps the
                                // kernel, as on the host (extract_child_interaction)
                                eval();
                                continue;
                            }
                            const DeviceMatrix& E = it->second;  // rows = hi, columns = lo
                            if (Rm > Sm) {
                                if (E.rows != kr || E.cols != ks) throw std::runtime_error("LevelEliminator::build_parent: child near block shape");
                                gather(E.ptr, 1, E.rows);
                            } else {
                                if (E.rows != ks || E.cols != kr) throw std::runtime_error("LevelEliminator::build_parent: child near block shape");
                                gather(E.ptr, E.rows, 1);
                            }
                            continue;
                        }
                        eval();
                        if (chebyshev != 2 || !reaches_two) continue;  // beyond two hops: kernel only
                        // two-hop: lazy fill of C(i, j) from i's eliminated neighbors next to j.
                        // CA levels: the pair's copy on the other rank must come out
                        // bitwise alike, so i is the lower child (not the local one)
                        // and the sources come by Morton (the common eliminated
                        // neighbors, writable on both ranks).
                        sources_of(im, jm, cj3, sources);
                        if (nonsym_) {
                            // C(R, S) -= temp2_E[R] X_RR_full_E temp4_E[:, S] = G_{E,R}^T P_{E,S}
                            const int64_t* r3 = rows_q ? cj3 : ci3;
                            const int64_t* s3 = rows_q ? ci3 : cj3;
                            size_t t = 0;
                            for (const SourceCandidate* src : sources) {
                                const BoxState& es = *src->state;
                                if (es.r == 0) continue;
                                if (fills.size() <= t) fills.resize(t + 1);
                                fills[t].push_back(FillTask{g_for(*src, Rm, r3), p_for(*src, Sm, s3), target, M.rows, kr,
                                                            ks, es.r, true});
                                ++t;
                            }
                            continue;
                        }
                        const bool i_first = !ca_ || im < jm;
                        const int64_t gi = i_first ? im : jm, pj = i_first ? jm : im;
                        const int64_t* gi3 = i_first ? ci3 : cj3;
                        const int64_t* pj3 = i_first ? cj3 : ci3;
                        if (ca_) {
                            std::sort(sources.begin(), sources.end(),
                                      [](const SourceCandidate* x, const SourceCandidate* y) { return x->morton < y->morton; });
                        }
                        size_t t = 0;
                        for (const SourceCandidate* src : sources) {
                            const BoxState& es = *src->state;
                            if (es.r == 0) continue;  // no Schur contribution
                            if (fills.size() <= t) fills.resize(t + 1);
                            fills[t].push_back(FillTask{g_for(*src, gi, gi3), p_for(*src, pj, pj3), target, M.rows, kr, ks,
                                                        es.r, Rm == gi});
                            ++t;
                        }
                    }
                }
            }
        }

        // the chunk's restored sources: [temp2 rows; X_RR_full^T] of each
        // (L = columns + r rows, r columns) gathered into pinned slices,
        // uploaded into a raw buffer and transposed into the source's block
        std::vector<char*> d_restore(pending.size());
        auto restore_ptr = [&](size_t off) -> const S* {  // staging offset -> device
            size_t lo = 0, hi = pending.size();
            while (hi - lo > 1) {
                const size_t mid = (lo + hi) / 2;
                (pending[mid].gall_off <= off ? lo : hi) = mid;
            }
            return reinterpret_cast<const S*>(d_restore[lo] + (off - pending[lo].gall_off));
        };
        char* d_raw = nullptr;
        std::vector<TransposeItem> restores;
        int max_rs_m = 0, max_rs_n = 0;
        if (!pending.empty()) {
            const auto t_restore = clock::now();
            d_raw = heap_.alloc(restore_bytes, ResidentKind::TrRestore);
            for (size_t q = 0; q < pending.size(); ++q) {
                const Pending& pd = pending[q];
                d_restore[q] = heap_.alloc_resident(pd.end - pd.gall_off, ResidentKind::Restore);
                const int L = static_cast<int>(pd.rows.size()) * (nonsym_ && !pd.xrr_only ? 2 : 1) + pd.r;
                restores.push_back(TransposeItem{reinterpret_cast<const S*>(d_raw + pd.gall_off), L,
                                                 reinterpret_cast<S*>(d_restore[q]), pd.r, L, pd.r});
                max_rs_m = std::max(max_rs_m, L);
                max_rs_n = std::max(max_rs_n, pd.r);
                Restored& rs = restored[pd.morton];
                if (pd.xrr_only) {  // (the packed restore's record: only X_RR_full comes from here)
                    rs.xrr = restore_ptr(pd.xrr_off);
                    rs.xrr_block = d_restore[q];
                    continue;
                }
                rs.gall = restore_ptr(pd.gall_off);
                rs.xrr = restore_ptr(pd.xrr_off);
                if (nonsym_) rs.g4all = restore_ptr(pd.g4_off);
                rs.col_of = std::move(pending_cols[pd.morton]);
                rs.block = d_restore[q];
                rs.last_parent = pd.last_parent;
            }
            constexpr size_t kSlice = size_t{128} << 20;
            for (size_t q0 = 0; q0 < pending.size();) {
                const size_t base = pending[q0].gall_off;
                size_t q1 = q0 + 1;
                while (q1 < pending.size() && pending[q1].end - base <= kSlice) ++q1;
                const size_t bytes = pending[q1 - 1].end - base;
                if (staged[stage_turn] == nullptr) {
                    check_cuda(cudaEventCreateWithFlags(&staged[stage_turn], cudaEventDisableTiming), "cudaEventCreate");
                } else {
                    check_cuda(cudaEventSynchronize(staged[stage_turn]), "restore staging");  // its last upload is done
                }
                char* h_stage = static_cast<char*>(pinned_pool().staging[stage_turn].reserve(bytes));
                #pragma omp parallel for schedule(dynamic)
                for (int64_t q = static_cast<int64_t>(q0); q < static_cast<int64_t>(q1); ++q) {
                    const Pending& pd = pending[static_cast<size_t>(q)];
                    const DataType* temp2 = pd.box->X_NR.data.data();
                    const int64_t ld = pd.box->X_NR.lda;
                    std::vector<DataType> rebuilt;  // (H2 host memory: from the host LU when the copy is not there)
                    if (!pd.box->X_RR_full.is_allocated()) rebuilt = host_xrr_full_from_lu(*pd.box);
                    const DataType* xrr = rebuilt.empty() ? pd.box->X_RR_full.data.data() : rebuilt.data();
                    const int64_t ldx = rebuilt.empty() ? pd.box->X_RR_full.lda : static_cast<int64_t>(pd.r);
                    const size_t nrows = pd.rows.size();
                    const bool with_t4 = nonsym_ && !pd.xrr_only;  // [temp2 rows; X_RR_full^T; temp4 columns^T]
                    const size_t L = nrows * (with_t4 ? 2 : 1) + static_cast<size_t>(pd.r);
                    const DataType* temp4 = with_t4 ? pd.box->X_RN.data.data() : nullptr;
                    const size_t ld4 = with_t4 ? static_cast<size_t>(pd.box->X_RN.lda) : 0;
                    DataType* h = reinterpret_cast<DataType*>(h_stage + (pd.gall_off - base));
                    for (int t = 0; t < pd.r; ++t) {
                        DataType* col = h + static_cast<size_t>(t) * L;
                        const DataType* src = temp2 + static_cast<size_t>(t) * ld;
                        for (size_t c = 0; c < nrows; ++c) col[c] = src[pd.rows[c]];
                        for (int u = 0; u < pd.r; ++u) col[nrows + u] = xrr[t + static_cast<size_t>(u) * ldx];
                        if (with_t4) {
                            DataType* c4 = col + nrows + static_cast<size_t>(pd.r);
                            for (size_t c = 0; c < nrows; ++c) c4[c] = temp4[static_cast<size_t>(t) + static_cast<size_t>(pd.rows[c]) * ld4];
                        }
                    }
                }
                check_cuda(cudaMemcpyAsync(d_raw + base, h_stage, bytes, cudaMemcpyHostToDevice, stream),
                           "restore fill sources");
                check_cuda(cudaEventRecord(staged[stage_turn], stream), "cudaEventRecord");
                stage_turn ^= 1;
                q0 = q1;
            }
            stats.bytes_up += static_cast<double>(restore_bytes);
            stats.tr_restore_bytes += static_cast<double>(restore_bytes);
            stats.tr_restore += std::chrono::duration<double>(clock::now() - t_restore).count();
        }

        // buffers of G and P for the chunk
        char* d_buffers = heap_.alloc(std::max<size_t>(buffer_bytes, 1), ResidentKind::TrBuffers);
        auto at = [&](size_t off) { return reinterpret_cast<S*>(d_buffers + off); };
        auto g_ptr = [&](int index) -> const S* {
            const GInfo& g = g_infos[static_cast<size_t>(index)];
            if (g.restore != SIZE_MAX) return restore_ptr(g.restore);
            return g.buffer == SIZE_MAX ? g.ptr : at(g.buffer);
        };
        for (const GPending& g : g_pending) {
            g_gathers.push_back(GatherItem{at(g.offset), g.r, g.r, g.k, g.src, 1, g.r, IndexList{}, g.cols});
            max_g_m = std::max(max_g_m, g.r);
            max_g_n = std::max(max_g_n, g.k);
        }
        VBatch<S> p_batch;
        for (const PPending& pp : p_pending) {
            const S* xrr = pp.xrr_restore == SIZE_MAX ? pp.xrr : restore_ptr(pp.xrr_restore);
            p_batch.entries.push_back({xrr, g_ptr(pp.g), at(pp.offset), pp.r, pp.k, pp.r, pp.r,
                                       g_infos[static_cast<size_t>(pp.g)].ld, pp.r});
        }
        // two batches per source rank: fills with rows j and fills with rows i
        std::vector<VBatch<S>> fill_j(fills.size()), fill_i(fills.size());
        for (size_t t = 0; t < fills.size(); ++t) {
            for (const FillTask& f : fills[t]) {
                const S* P = reinterpret_cast<const S*>(d_buffers + p_offsets[static_cast<size_t>(f.p_j)]);
                const int g_ld = g_infos[static_cast<size_t>(f.g_i)].ld;
                if (!f.rows_i) {
                    // out (m x n) -= P_{E,j}^T (m x r) G_{E,i} (r x n)
                    fill_j[t].entries.push_back({P, g_ptr(f.g_i), f.out, f.m, f.n, f.r, f.r, g_ld, f.ld});
                } else {
                    // out (m x n) -= G_{E,i}^T (m x r) P_{E,j} (r x n): the same terms, transposed
                    fill_i[t].entries.push_back({g_ptr(f.g_i), P, f.out, f.m, f.n, f.r, g_ld, f.r, f.ld});
                }
                ++fill_gemms;
                stats.transition_flops += kFlopScale * 2.0 * f.m * f.n * f.r;
            }
        }
        const size_t off_restores = meta.append(restores);
        const size_t off_xrr_split = meta.append(xrr_split);
        const size_t off_xrr_swap = meta.append(xrr_swap);
        xrr_gemm.stage(meta);
        const size_t off_evals = meta.append(evals);
        const size_t off_gathers = meta.append(gathers);
        const size_t off_g = meta.append(g_gathers);
        p_batch.stage(meta);
        for (auto& b : fill_j) b.stage(meta);
        for (auto& b : fill_i) b.stage(meta);
        stats.tr_plan += std::chrono::duration<double>(clock::now() - t_plan).count();
        {
            size_t bytes = 0;
            for (const auto& [ptr, b] : chunk_blocks) {
                bytes += b;
                if (!tr_block_seen_.emplace(ptr).second) ++stats.tr_block_revisits;
            }
            stats.tr_block_bytes_sum += static_cast<double>(bytes);
            stats.tr_block_bytes_max_chunk = std::max(stats.tr_block_bytes_max_chunk, static_cast<double>(bytes));
            stats.tr_blocks_max_chunk = std::max(stats.tr_blocks_max_chunk,
                                                 static_cast<int64_t>(chunk_blocks.size()));
        }
        ++stats.tr_chunks;
        char* md = meta.upload(meta_device, stream);
        stats.bytes_up += static_cast<double>(meta.size());
        chunk_marks.push_back(std::make_unique<StreamMarks>());
        StreamMarks& tmarks = *chunk_marks.back();
        tmarks.mark(stream);
        run_fill_batch(packed_src, packed_dst, stream);  // the chunk's compressed sources, put back
        if (!xrr_split.empty()) {  // and the X_RR_full its P products read
            launch_lu_split(reinterpret_cast<const LuSplitItemT<S>*>(md + off_xrr_split),
                            static_cast<int>(xrr_split.size()), max_xrr_r, stream);
            xrr_gemm.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
            launch_row_swaps(reinterpret_cast<const RowSwapItemT<S>*>(md + off_xrr_swap),
                             static_cast<int>(xrr_swap.size()), stream);
        }
        launch_transpose(reinterpret_cast<const TransposeItem*>(md + off_restores), static_cast<int>(restores.size()),
                         max_rs_m, max_rs_n, stream);
        {
            EntryOriginScope origin(EntryOrigin::transition);  // (collective entries' profile)
            eval_blocks(evals, meta, md, off_evals, max_eval_m, max_eval_n, stream);
        }
        launch_gather(reinterpret_cast<const GatherItem*>(md + off_gathers), static_cast<int>(gathers.size()),
                      max_gather_m, max_gather_n, md, stream);
        launch_gather(reinterpret_cast<const GatherItem*>(md + off_g), static_cast<int>(g_gathers.size()), max_g_m,
                      max_g_n, md, stream);
        tmarks.mark(stream);
        p_batch.gemm(md, MagmaNoTrans, MagmaNoTrans, 1.0, 0.0, queue);
        tmarks.mark(stream);
        for (size_t t = 0; t < fills.size(); ++t) {
            fill_j[t].gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
            fill_i[t].gemm(md, MagmaTrans, MagmaNoTrans, -1.0, 1.0, queue);
        }
        tmarks.mark(stream);
#ifdef H2_HAVE_NVCOMP
        // H2_compress_factors: this chunk's parent blocks, compressed.  Every
        // block is written by exactly one chunk (the parents are partitioned
        // into contiguous ranges, 0 revisits measured), so a chunk's blocks
        // are final here.  The raw blocks stay for now - the next level still
        // reads them - so this must not change any result.
        ++stats.tr_pack_chunks;
        if (compress_factors_ >= 1 && pack_parents_ && !chunk_refs.empty() && under_pressure()) {
            ++stats.tr_pack_chunks_packed;
            check_cuda(cudaStreamSynchronize(stream), "parent blocks");  // written by the launches above
            std::vector<const S*> srcs;
            std::vector<size_t> elems;
            std::vector<size_t> ref_of;
            size_t group = 0;
            auto flush = [&]() {
                if (srcs.empty()) return;
                std::vector<CompressedXnr> comp;
                if (fill_compressor_.compress(srcs, elems, tolerance_ * kXnrToleranceFraction, comp, stream)) {
                    const auto& cks = fill_compressor_.chunks();
                    const char* staging = fill_compressor_.staging();
                    std::vector<CopySpan> spans;
                    int64_t max_words = 0;
                    std::vector<std::pair<bool, uint64_t>> raw_to_free;  // (schur?, morton or edge key)
                    size_t group_packed = 0;
                    for (size_t j = 0; j < srcs.size(); ++j) group_packed += comp[j].bytes;
                    // (the level below adopts these: its host region)
                    // the host's room first: these copies only save device memory
                    // at the next level, and here the level's own near and fill
                    // stores are still on the host (768^3 CA0 on 64 ranks: 7 GB
                    // more a rank over the transition, the interior ranks
                    // killed at 62 GB of a 60 GB share)
                    const bool host_room = group_packed == 0 ||
                                           host_room_bytes() + static_cast<double>(pinned_chunk_pool().idle()) >
                                               kHostReserveBytes + static_cast<double>(group_packed);
                    if (!host_room) ++stats.pb_host_short;
                    char* owned = group_packed > 0 && host_room ? near_host_arena(level_index_ - 1).alloc(group_packed) : nullptr;
                    if (owned != nullptr) out->packed_owned.push_back(owned);
                    size_t owned_at = 0;
                    for (size_t j = 0; j < srcs.size(); ++j) {
                        const ChunkRef& ref = chunk_refs[ref_of[j]];
                        if (comp[j].bytes == 0) continue;
                        if (owned == nullptr) {  // no room: these stay raw only
                            ++stats.pb_packed_short;
                            continue;
                        }
                        char* packed = owned + owned_at;
                        owned_at += comp[j].bytes;
                        PackedBlock pb;
                        pb.packed = packed;
                        pb.live = true;   // its bytes are current: without this
                        // every adopted block reads as "not packed" to
                        // is_packed_edge(), the wave never materializes it,
                        // and the fills gather from a null pointer
                        pb.bytes = comp[j].bytes;
                        pb.chunk0 = static_cast<uint32_t>(out->strides.size());
                        pb.nchunks = static_cast<uint32_t>(comp[j].nchunks);
                        pb.rows = ref.m->rows;
                        pb.cols = ref.m->cols;
                        size_t off = 0;
                        for (size_t c = 0; c < comp[j].nchunks; ++c) {
                            const auto& ck = cks[comp[j].chunk0 + c];
                            spans.push_back(CopySpan{staging + ck.offset, packed + off,
                                                     static_cast<int64_t>(ck.stride / 4)});
                            max_words = std::max(max_words, static_cast<int64_t>(ck.stride / 4));
                            out->strides.push_back(static_cast<uint32_t>(ck.stride));
                            off += ck.stride;
                        }
                        if (ref.schur) {
                            out->schur_packed[ref.morton] = pb;
                        } else {
                            out->edges_packed[ref.key] = pb;
                        }
                        raw_to_free.push_back({ref.schur, ref.schur ? static_cast<uint64_t>(ref.morton) : ref.key});
                        stats.pb_raw_bytes += static_cast<double>(elems[j] * sizeof(S));
                        stats.pb_stored_bytes += static_cast<double>(comp[j].bytes);
                        ++stats.pb_blocks;
                    }
                    if (!spans.empty()) {
                        const size_t list_bytes = spans.size() * sizeof(CopySpan);
                        if (list_bytes > fill_meta_device_.capacity()) {
                            check_cuda(cudaStreamSynchronize(stream), "parent blocks");
                            fill_meta_device_.reserve(list_bytes);
                        }
                        auto* d_items = static_cast<CopySpan*>(fill_meta_device_.reserve(list_bytes));
                        check_cuda(cudaMemcpyAsync(d_items, spans.data(), list_bytes, cudaMemcpyHostToDevice,
                                                   stream), "parent blocks");
                        launch_copy_spans(d_items, static_cast<int>(spans.size()), max_words, stream);
                        // the gather reads the staging, which goes back next
                        check_cuda(cudaStreamSynchronize(stream), "parent blocks");
                        // and the raw blocks: by key, and only what the heap
                        // still holds, so a repeat is a no-op rather than an
                        // error (measured clean - 7099 freed, 0 stale)
                        for (const auto& [is_schur, key] : raw_to_free) {
                            DeviceMatrix* m = nullptr;
                            if (is_schur) {
                                auto it = out->schur.find(static_cast<int64_t>(key));
                                if (it != out->schur.end()) m = &it->second;
                            } else {
                                auto it = out->edges.find(key);
                                if (it != out->edges.end()) m = &it->second;
                            }
                            if (m == nullptr || m->ptr == nullptr) continue;
                            if (!heap_.is_allocated(m->ptr)) { m->ptr = nullptr; continue; }
                            heap_.free(m->ptr);
                            m->ptr = nullptr;
                        }
                    }
                }
                fill_compressor_.release();
                srcs.clear();
                elems.clear();
                ref_of.clear();
                group = 0;
            };
            for (size_t i = 0; i < chunk_refs.size(); ++i) {
                const DeviceMatrix& m = *chunk_refs[i].m;
                if (m.ptr == nullptr || m.rows <= 0 || m.cols <= 0) continue;
                const size_t bytes = static_cast<size_t>(m.rows) * m.cols * D;
                if (!srcs.empty() && group + bytes > kFillStagingBudget) flush();
                srcs.push_back(m.ptr);
                elems.push_back(static_cast<size_t>(m.rows) * m.cols);
                ref_of.push_back(i);
                group += bytes;
            }
            flush();
        }
#endif
        heap_.free(d_buffers);  // later launches are ordered after this chunk
        heap_.free(d_raw);
        for (char* b : xrr_exp_blocks) heap_.free(b);  // (the chunk's P products read them; the next chunk rebuilds its own)
        xrr_exp_blocks.clear();
        // restored sources whose last reader was in this chunk or before
        for (auto it = restored.begin(); it != restored.end();) {
            if (it->second.last_parent < p1) {
                heap_.free(it->second.block);
                if (it->second.xrr_block != nullptr) heap_.free(it->second.xrr_block);
                it = restored.erase(it);
            } else {
                ++it;
            }
        }
        if (heap_.peak() > stats.heap_peak) {
        stats.heap_peak = heap_.peak();
        stats.heap_peak_resident = heap_.resident_at_peak();
        stats.heap_peak_transient = heap_.transient_at_peak();
    }
        if (env::trace("wave")) {   // the chunk against its budget (BPACK_TRACE=wave)
            std::printf("  [gpu] level %d transition chunk %zu: parents %zu..%zu of %zu, budget %.2f GB = buffers %.2f + "
                        "host restore %.2f + packed restore %.2f + X_RR_full rebuilt %.2f + parent blocks %.2f; heap used "
                        "%.2f GB (restored %.2f, parents %.2f, kept %.2f), largest free %.2f; host RSS %.2f GB, pinned %.2f, rank %d\n", level_index_,
                        chunk_marks.size(), p0, p1, np,
                        static_cast<double>(budget) / 1e9, static_cast<double>(buffer_bytes) / 1e9,
                        static_cast<double>(restore_bytes) / 1e9, static_cast<double>(packed_restore_bytes) / 1e9,
                        static_cast<double>(xrr_rebuild_bytes) / 1e9,
                        static_cast<double>(chunk_block_bytes) / 1e9, static_cast<double>(heap_.used()) / 1e9,
                        static_cast<double>(heap_.resident_now(ResidentKind::Restore)) / 1e9,
                        static_cast<double>(heap_.resident_now(ResidentKind::Transition)) / 1e9,
                        static_cast<double>(heap_.resident_now(ResidentKind::KeptFactors)) / 1e9,
                        static_cast<double>(heap_.largest_free()) / 1e9,
                        tr_rss_gb(), static_cast<double>(pinned_chunk_pool().pinned()) / 1e9, tr_rank());
            std::fflush(stdout);
        }
        p0 = p1;
    }
    for (auto& kv : restored) {
        heap_.free(kv.second.block);
        if (kv.second.xrr_block != nullptr) heap_.free(kv.second.xrr_block);
    }
    for (char* b : xrr_exp_blocks) heap_.free(b);
    for (cudaEvent_t e : staged) {
        if (e != nullptr) cudaEventDestroy(e);
    }
    stats.transition_fill_gemms += fill_gemms;
    check_cuda(cudaStreamSynchronize(stream), "device transition");
    for (const auto& m : chunk_marks) {
        stats.tr_blocks += m->seconds(0);
        stats.tr_p += m->seconds(1);
        stats.tr_fill += m->seconds(2);
    }

    // the child level is done: release its device data (not the new blocks)
    release_level_data();
    stats.transition += std::chrono::duration<double>(clock::now() - t0).count();
    return out;
}

// Free this level's device data (blocks, fill sources, point table).  Blocks
// handed to the parent level are not the level's any more.
template<typename CoordType, typename DataType, typename KernelType>
void LevelEliminator<CoordType, DataType, KernelType>::release_level_data() {
    released_sources_.clear();  // their blocks are freed with the states below
    if (!heap_initialized_) return;
    // The kept X_NR hand-over (transpose_kept_xnr) reads only the fill
    // sources and needs room for the boxes' compressed X_NR, so the level's
    // blocks go first (2026-10-11: at 768^3 CA0 an interior rank came to it
    // with 8.4 GB free and 12 GB to place, the near and Schur blocks, the
    // sketch buffers and the adopted blocks still on the card).
    evict_fill_copies(true);
    release_owner_sources();
    if (outbox_ != nullptr) free_block(outbox_);
    outbox_ = nullptr;
    for (OwnerWaveSend& w : owner_sends_) {
        if (w.issued && !w.reqs.empty()) MPI_Waitall(static_cast<int>(w.reqs.size()), w.reqs.data(), MPI_STATUSES_IGNORE);
        if (w.block != nullptr) free_block(w.block);
        if (w.ready != nullptr) cudaEventDestroy(w.ready);
    }
    owner_sends_.clear();
    for (OwnerWaveRecv& r : owner_recvs_) {
        if (!r.reqs.empty()) MPI_Waitall(static_cast<int>(r.reqs.size()), r.reqs.data(), MPI_STATUSES_IGNORE);
        if (r.buf != nullptr) free_block(r.buf);
    }
    owner_recvs_.clear();
    for (OwnerPendingFree& pf : owner_pending_free_) {
        if (pf.done != nullptr) {
            cudaEventSynchronize(pf.done);
            cudaEventDestroy(pf.done);
        }
        for (char* b : pf.blocks) free_block(b);
    }
    owner_pending_free_.clear();
    owner_sent_count_.clear();
    owner_recv_count_.clear();
    owner_marks_.clear();
    owner_sent_bytes_.clear();
#ifdef H2_HAVE_NVCOMP
    close_write_level();
    report_near_sizes();
    eliminator_stats().mat_left += static_cast<int64_t>(packed_schur_.size() + packed_edges_.size());
    for (char* b : packed_owned_) {
        if (heap_.owns(b)) heap_.free(b);  // (a host run goes with the region, reset by release_near_store)
    }
    packed_owned_.clear();
    packed_live_.clear();
    release_near_store();
    packed_schur_.clear();
    packed_edges_.clear();
    packed_strides_.clear();
#endif
    for (auto& kv : edges_) heap_.free(kv.second.ptr);
    edges_.clear();
    for (BoxState& st : states_) {
        heap_.free(st.schur.ptr);
        st.schur = DeviceMatrix{};
    }
    for (BoxState& st : remote_states_) {
        heap_.free(st.schur.ptr);
        st.schur = DeviceMatrix{};
    }
    if (wave_sketch_ != nullptr) {
        heap_.free(wave_sketch_);
        wave_sketch_ = nullptr;
    }
    for (char* p : wave_sketch_extra_) heap_.free(p);
    wave_sketch_extra_.clear();
    if (adopt_) {
        for (auto& kv : adopt_->schur) heap_.free(kv.second.ptr);
        for (auto& kv : adopt_->edges) heap_.free(kv.second.ptr);
        adopt_.reset();
    }
    transpose_kept_xnr();  // the kept X_NR, while the fill sources are still here
    for (BoxState& st : states_) {
        heap_.free(st.persist);
        st.persist = nullptr;
        release_fill_temp2t(st);  // its temp2^T block, raw or compressed
    }
    for (BoxState& st : remote_states_) {
        heap_.free(st.persist);
        st.persist = nullptr;
        release_fill_temp2t(st);
    }
    fill_host_arena().reset();  // the level's packed temp2^T on the host, all of it read by now
    for (char* p : level_allocs_) heap_.free(p);
    level_allocs_.clear();
    heap_.free(d_xyz_);
    heap_.free(d_ids_);
    d_xyz_ = nullptr;
    d_ids_ = nullptr;
    heap_initialized_ = false;
}

// Device blocks of a level into its host boxes (for a level that runs on
// the host, e.g. level 1 before the root, or a CA level's halo gather):
// Schur blocks, and each near-field block as the lower box's view (rows =
// the higher box) with an empty reciprocal for share_symmetric_level_edges.
// The device copies are freed.
//
// The blocks stream through two pinned ring slots of the (then idle) level
// copier: one slot's copies overlap the host unpacking of the other's, in
// parallel over blocks.  A block over a slot is copied on its own.
template<typename CoordType, typename DataType>
void download_level_blocks(DeviceLevelBlocks<DataType>& blocks, std::vector<BoxData<CoordType, DataType>>& boxes) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)blocks; (void)boxes;
    throw std::runtime_error("download_level_blocks: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    static_assert(sizeof(S) == sizeof(DataType), "device and host elements share their layout");
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    if (boxes.empty()) return;
    const int64_t start = boxes.front().morton_index;
    auto box_of = [&](int64_t morton) -> BoxData<CoordType, DataType>& {
        const int64_t idx = morton - start;
        if (idx < 0 || idx >= static_cast<int64_t>(boxes.size())) {
            throw std::runtime_error("download_level_blocks: box is not local");
        }
        return boxes[static_cast<size_t>(idx)];
    };
    check_cuda(cudaStreamSynchronize(stream), "download level blocks");

    // every block's data on the host: the Schur blocks, then the edges
    std::vector<const DeviceMatrixT<S>*> pieces;
    pieces.reserve(blocks.schur.size() + blocks.edges.size());
    for (auto& kv : blocks.schur) pieces.push_back(&kv.second);
    for (auto& kv : blocks.edges) pieces.push_back(&kv.second);
    std::vector<std::vector<DataType>> data(pieces.size());
    {
        auto bytes_of = [&](size_t i) { return static_cast<size_t>(pieces[i]->rows) * pieces[i]->cols * sizeof(S); };
        constexpr size_t slot_bytes = PinnedPool::kRingSlot;
        char* slot[2] = {static_cast<char*>(pinned_pool().ring[0].reserve(slot_bytes)),
                         static_cast<char*>(pinned_pool().ring[1].reserve(slot_bytes))};
        cudaEvent_t copied[2] = {nullptr, nullptr};
        for (auto& ev : copied) check_cuda(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming), "cudaEventCreate");
        // the pieces [first, last) in slot s, at consecutive offsets
        auto unpack = [&](int s, size_t first, size_t last) {
            check_cuda(cudaEventSynchronize(copied[s]), "download level blocks");
            std::vector<size_t> offset(last - first + 1, 0);
            for (size_t i = first; i < last; ++i) offset[i - first + 1] = offset[i - first] + bytes_of(i);
            #pragma omp parallel for schedule(dynamic, 8)
            for (int64_t i = static_cast<int64_t>(first); i < static_cast<int64_t>(last); ++i) {
                const DataType* src = reinterpret_cast<const DataType*>(slot[s] + offset[static_cast<size_t>(i) - first]);
                data[static_cast<size_t>(i)].assign(src, src + bytes_of(static_cast<size_t>(i)) / sizeof(S));
            }
        };
        int s = 0;
        bool pending = false;
        size_t pending_first = 0, pending_last = 0;
        size_t i = 0;
        while (i < pieces.size()) {
            if (bytes_of(i) > slot_bytes) {
                data[i].resize(bytes_of(i) / sizeof(S));
                check_cuda(cudaMemcpy(data[i].data(), pieces[i]->ptr, bytes_of(i), cudaMemcpyDeviceToHost),
                           "download level blocks");
                ++i;
                continue;
            }
            const size_t first = i;
            size_t used = 0;
            while (i < pieces.size() && bytes_of(i) <= slot_bytes && used + bytes_of(i) <= slot_bytes) {
                check_cuda(cudaMemcpyAsync(slot[s] + used, pieces[i]->ptr, bytes_of(i), cudaMemcpyDeviceToHost, stream),
                           "download level blocks");
                used += bytes_of(i);
                ++i;
            }
            check_cuda(cudaEventRecord(copied[s], stream), "cudaEventRecord");
            if (pending) unpack(1 - s, pending_first, pending_last);
            pending = true;
            pending_first = first;
            pending_last = i;
            s = 1 - s;
        }
        if (pending) unpack(1 - s, pending_first, pending_last);
        for (auto& ev : copied) cudaEventDestroy(ev);
    }

    size_t piece = 0;
    for (auto& kv : blocks.schur) {
        box_of(kv.first).schur_complement.set_owned(kv.second.rows, kv.second.cols, std::move(data[piece++]),
                                                    MatrixStorage<DataType>::FULL);
        heap.free(kv.second.ptr);
    }
    auto local_box = [&](int64_t morton) -> BoxData<CoordType, DataType>* {
        const int64_t idx = morton - start;
        return idx >= 0 && idx < static_cast<int64_t>(boxes.size()) ? &boxes[static_cast<size_t>(idx)] : nullptr;
    };
    for (auto& kv : blocks.edges) {
        const int64_t lo = static_cast<int64_t>(kv.first >> 32), hi = static_cast<int64_t>(kv.first & 0xffffffffu);
        std::vector<DataType> block_data = std::move(data[piece++]);
        if (blocks.directed) {
            // K(row, col) (row = the low half): the column box's A_NS[row], the
            // row box's A_SN[col]
            const int64_t col = lo, row = hi;
            auto* col_box = local_box(col);
            auto* row_box = local_box(row);
            if (col_box == nullptr && row_box == nullptr) {
                throw std::runtime_error("download_level_blocks: near block without a local box");
            }
            if (col_box != nullptr) {
                ModifiedBlock<DataType> view;
                view.neighbor_morton = row;
                std::vector<DataType> copy = row_box != nullptr ? block_data : std::vector<DataType>{};
                view.set_a_ns_owned(kv.second.rows, kv.second.cols,
                                    row_box != nullptr ? std::move(copy) : std::move(block_data),
                                    MatrixStorage<DataType>::FULL);
                col_box->near_field_interaction_map[row] =
                    static_cast<int64_t>(col_box->near_field_modified_interactions.size());
                col_box->near_field_modified_interactions.push_back(std::move(view));
            }
            if (row_box != nullptr) {
                ModifiedBlock<DataType> view;
                view.neighbor_morton = col;
                view.A_SN.set_owned(kv.second.rows, kv.second.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
                row_box->near_field_interaction_map_nonsymmetry[col] =
                    static_cast<int64_t>(row_box->near_field_modified_interactions.size());
                row_box->near_field_modified_interactions.push_back(std::move(view));
            }
            heap.free(kv.second.ptr);
            continue;
        }
        if (local_box(lo) == nullptr || local_box(hi) == nullptr) {
            // a pair with a box of another rank: the local box's view
            // (rows = the other box's points)
            ModifiedBlock<DataType> view;
            if (auto* lb = local_box(lo)) {
                view.neighbor_morton = hi;
                view.set_a_ns_owned(kv.second.rows, kv.second.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
                lb->near_field_interaction_map[hi] = static_cast<int64_t>(lb->near_field_modified_interactions.size());
                lb->near_field_modified_interactions.push_back(std::move(view));
            } else if (auto* hb = local_box(hi)) {
                std::vector<DataType> transposed(block_data.size());
                for (int j = 0; j < kv.second.cols; ++j)
                    for (int r = 0; r < kv.second.rows; ++r)
                        transposed[j + static_cast<size_t>(r) * kv.second.cols] = block_data[r + static_cast<size_t>(j) * kv.second.rows];
                view.neighbor_morton = lo;
                view.set_a_ns_owned(kv.second.cols, kv.second.rows, std::move(transposed), MatrixStorage<DataType>::FULL);
                hb->near_field_interaction_map[lo] = static_cast<int64_t>(hb->near_field_modified_interactions.size());
                hb->near_field_modified_interactions.push_back(std::move(view));
            } else {
                throw std::runtime_error("download_level_blocks: edge without a local box");
            }
            heap.free(kv.second.ptr);
            continue;
        }
        auto& lo_box = box_of(lo);
        auto& hi_box = box_of(hi);
        ModifiedBlock<DataType> lo_block;
        lo_block.neighbor_morton = hi;
        lo_block.set_a_ns_owned(kv.second.rows, kv.second.cols, std::move(block_data), MatrixStorage<DataType>::FULL);
        lo_box.near_field_interaction_map[hi] = static_cast<int64_t>(lo_box.near_field_modified_interactions.size());
        lo_box.near_field_modified_interactions.push_back(std::move(lo_block));
        ModifiedBlock<DataType> hi_block;
        hi_block.neighbor_morton = lo;
        hi_box.near_field_interaction_map[lo] = static_cast<int64_t>(hi_box.near_field_modified_interactions.size());
        hi_box.near_field_modified_interactions.push_back(std::move(hi_block));
        heap.free(kv.second.ptr);
    }
    blocks.schur.clear();
    blocks.edges.clear();
  }
}

// LU of the root block, left on the device by the device transition of
// level 1, with the factors copied to the root box (as the host does:
// X_RR = LU, 1-based pivots; of the block's symmetric part when `symmetric`).
// The device block is freed.
template<typename CoordType, typename DataType>
void factor_root_on_device(DeviceLevelBlocks<DataType>& blocks, BoxData<CoordType, DataType>& root, bool symmetric) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)blocks; (void)root; (void)symmetric;
    throw std::runtime_error("factor_root_on_device: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    auto it = blocks.schur.find(root.morton_index);
    if (it == blocks.schur.end() || blocks.schur.size() != 1 || !blocks.edges.empty()) {
        throw std::runtime_error("factor_root_on_device: expected the root block alone");
    }
    const DeviceMatrixT<S> m = it->second;
    if (m.rows != root.num_points || m.cols != root.num_points) {
        throw std::runtime_error("factor_root_on_device: root block shape");
    }
    Context& ctx = Context::instance();
    ctx.activate();
    const magma_int_t n = m.rows;
    if (symmetric) launch_symmetrize(m.ptr, static_cast<int>(n), static_cast<int>(n), ctx.stream());
    check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    std::vector<magma_int_t> piv(static_cast<size_t>(n));
    magma_int_t info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        magma_dgetrf_gpu(n, n, m.ptr, n, piv.data(), &info);
    } else {
        magma_zgetrf_gpu(n, n, reinterpret_cast<magmaDoubleComplex*>(m.ptr), n, piv.data(), &info);
    }
    if (info != 0) {
        throw std::runtime_error("factor_root_on_device: LU factorization of root failed with INFO = " +
                                 std::to_string(info));
    }
    std::vector<DataType> lu(static_cast<size_t>(n) * static_cast<size_t>(n));
    check_cuda(cudaMemcpy(lu.data(), m.ptr, lu.size() * sizeof(S), cudaMemcpyDeviceToHost), "root factors");
    DeviceHeap::instance().free(m.ptr);
    blocks.schur.clear();
    root.X_RR.set_owned(n, n, std::move(lu), MatrixStorage<DataType>::LU_FACTORED);
    root.X_RR_pivots.assign(piv.begin(), piv.end());
  }
}

// LU of a root block assembled on the host (multi-rank runs, whose level-1
// blocks pass through the host for the process reduction) on the device,
// when the heap has room for it: the root box gets the factors as the host
// LU leaves them (X_RR = LU, 1-based pivots; of the block's symmetric part
// when `symmetric`).  Returns false, with nothing changed, otherwise.
template<typename CoordType, typename DataType>
bool factor_host_root_on_device(BoxData<CoordType, DataType>& root, bool symmetric) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)root; (void)symmetric;
    return false;
  } else {
    using S = typename DeviceScalar<DataType>::type;
    const MatrixStorage<DataType>& A = root.schur_complement;
    if (!A.is_allocated() || A.rows <= 0 || A.cols != A.rows || A.lda != A.rows) return false;
    DeviceHeap& heap = DeviceHeap::instance();
    if (!heap.initialized()) return false;
    const magma_int_t n = static_cast<magma_int_t>(A.rows);
    const size_t bytes = static_cast<size_t>(n) * static_cast<size_t>(n) * sizeof(S);
    Context& ctx = Context::instance();
    ctx.activate();
    check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    char* d = heap.try_alloc(bytes);
    if (d == nullptr) return false;
    check_cuda(cudaMemcpy(d, A.data.data(), bytes, cudaMemcpyHostToDevice), "root block");
    if (symmetric) {
        launch_symmetrize(reinterpret_cast<S*>(d), static_cast<int>(n), static_cast<int>(n), ctx.stream());
        check_cuda(cudaStreamSynchronize(ctx.stream()), "root block");
    }
    std::vector<magma_int_t> piv(static_cast<size_t>(n));
    magma_int_t info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        magma_dgetrf_gpu(n, n, reinterpret_cast<double*>(d), n, piv.data(), &info);
    } else {
        magma_zgetrf_gpu(n, n, reinterpret_cast<magmaDoubleComplex*>(d), n, piv.data(), &info);
    }
    if (info != 0) {
        heap.free(d);
        throw std::runtime_error("factor_host_root_on_device: LU factorization of root failed with INFO = " +
                                 std::to_string(info));
    }
    std::vector<DataType> lu(static_cast<size_t>(n) * static_cast<size_t>(n));
    check_cuda(cudaMemcpy(lu.data(), d, bytes, cudaMemcpyDeviceToHost), "root factors");
    heap.free(d);
    root.X_RR.set_owned(n, n, std::move(lu), MatrixStorage<DataType>::LU_FACTORED);
    root.X_RR_pivots.assign(piv.begin(), piv.end());
    return true;
  }
}

// Whether the device box path would take level `level_index` (same test as
// make_level_eliminator), for deciding where the parent blocks go.
template<typename CoordType, typename DataType, typename KernelType>
bool level_eliminator_would_run(const ParallelTree<CoordType, DataType>* tree, int level_index, const KernelType* kernel,
                                FactorizationMethod method, int lazy_schur = -1) {
    if constexpr (gpu_data_type<DataType>) {
        return level_eliminator_supported(tree->levels[static_cast<size_t>(level_index)], kernel, tree->dimension,
                                          method, nullptr, false, lazy_schur);
    } else {
        (void)tree; (void)level_index; (void)kernel; (void)method; (void)lazy_schur;
        return false;
    }
}

// The device box path for this level, uploaded and ready, or nullptr (with
// the reason) when the level stays on the host.
template<typename CoordType, typename DataType, typename KernelType>
std::unique_ptr<LevelEliminatorBase<CoordType, DataType>> make_level_eliminator(
    ParallelTree<CoordType, DataType>* tree, int level_index, KernelType* kernel, double tolerance,
    FactorizationMethod method, std::string* reason, std::unique_ptr<DeviceLevelBlocks<DataType>> adopt = nullptr,
    bool occupancy = false, bool device_sketch = true, bool ca = false, bool owner = false,
    int compress_factors = 0) {
    if constexpr (gpu_data_type<DataType>) {
        const auto& level = tree->levels[static_cast<size_t>(level_index)];
        if (!level_eliminator_supported(level, kernel, tree->dimension, method, reason, ca || owner)) {
            if (adopt) throw std::runtime_error("make_level_eliminator: device blocks for a level the device cannot run");
            return nullptr;
        }
        if (!kernel_symmetric(kernel, 0) && !device_sketch) {
            if (reason) *reason = "a nonsymmetric matrix's boxes over the device sketch's size run on the host";
            if (adopt) throw std::runtime_error("make_level_eliminator: device blocks for a level the device cannot run");
            return nullptr;
        }
        auto eliminator = std::make_unique<LevelEliminator<CoordType, DataType, KernelType>>(
            tree, level_index, kernel, tolerance, std::move(adopt), occupancy, device_sketch, ca, owner,
            compress_factors);
        if (!owner) eliminator->begin();  // (an owner level begins after its plan)
        return eliminator;
    } else {
        (void)tree; (void)level_index; (void)kernel; (void)tolerance; (void)method; (void)adopt; (void)occupancy;
        (void)device_sketch; (void)ca; (void)owner; (void)compress_factors;
        if (reason) *reason = "single precision is not run on the GPU";
        return nullptr;
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
