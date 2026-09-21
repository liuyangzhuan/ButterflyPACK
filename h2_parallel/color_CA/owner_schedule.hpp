#pragma once
// ---------------------------------------------------------------------------
// owner_schedule.hpp — component-owner elimination, asynchronous per-component
// schedule (H2_CA_owner_component=3). Design: concise-algorithm-ref/
// component-owner-dataflow-design.md §6–§10.
//
// Every boundary component (blue cluster incl. orange, purple run, green face)
// is eliminated once by one owner.  Execution is counter driven
// (dataflow::Runtime): an owned component starts when its predecessors'
// generators are available here; a non-owned shared component is "installed"
// (its shipped pre-pass snapshots written into the local/ghost boxes and its
// passes replayed) when its payload has arrived and its predecessors are
// available.  No per-wave barriers; colour barriers only on coarse levels
// (synchronous colour schedule, owner_schedule_run_boundary).
//
// Batches: everything that is runnable at once and of one colour — owned
// components to eliminate and arrived components to install — is processed
// as ONE batch with merged waves (same-colour components never share a live
// pair, so this is exact) and one R2 replay.  This keeps the thread team busy
// where per-component waves would starve it.  On coarse levels (fewer than
// 32 owned boxes per thread) a colour is instead processed as exactly two
// batches — all owned components, then all shared copies — after waiting for
// completeness, so the waves are as full as they can be.
//
// Waves: inside a colour, boxes are eliminated in waves of pairwise
// non-adjacent boxes.  The waves are a per-component greedy partition
// (owner_wave_partition: a face sheet needs ~4 classes, an edge 2) instead of
// the 8 Morton-parity classes; the partition is a pure function of the
// shared component graph, hence identical on all ranks, and it defines the
// elimination sequence (seq_of) that every deferred replay orders by.
//
// Determinism (bitwise equal to the replicated CA order for every pair of
// interest) comes from the pair-liveness rules — see `pair_watermark`:
//   R1  while a batch B is processed its sources update, in wave order,
//       exactly the pairs with an endpoint in B;
//   R2  before B starts, every already-available source S not in B is applied
//       in canonical (colour, wave, morton) order to the pairs with an
//       endpoint in B — but only those S above the pair's watermark, i.e. not
//       already applied when the pair's other endpoint's component started;
//   R3  at the interior start all pairs of interest receive, in canonical
//       order, every source above their watermark; boundary sources are then
//       finalized (temp2 -> X_NR) and the interior runs as today.
// A pair is "of interest" when an endpoint is local or in a component this
// rank owns; other pairs are never touched here.
//
// Routing: a box's state goes to its home rank (FULL: everything the
// transition and the solve read), to every rank holding a pair of interest
// with it — a rank with a local or owned box adjacent to it — as a COMPACT
// state (no geometry, near blocks only toward that rank's interest boxes),
// and to every rank with an interest box within TWO hops as a SKELETON
// section (index sets only): the lazy far-field regeneration of a pair
// (B, F) reads the far endpoint F's skeleton indices, and F is two hops from
// B.  Recipients install only what they were sent and mark only that as
// eliminated.
//
// Include order: after factorization.hpp (uses its pass functions).
//
// INVARIANTS — re-read before touching the event loop, the payload rules or
// the filters (the unit test covers 1, 2, 4; the serial oracle covers 3):
//  1. Routing consistency: every recipient of a box's state expects it — a
//     sharer of the component (FULL/COMPACT/SKELETON) or a rank listing the
//     component in auxiliary_components() (SKELETON).  dataflow::
//     payload_recipients is the single source of truth for both sides.
//  2. Readiness: an owned component starts only when its 1-hop predecessors
//     (pairs) and 2-hop lower-colour predecessors (far skeletons) are
//     available here; a shared component is installed only when its 1-hop
//     predecessors are.  Availability = installed/eliminated here (or, for
//     auxiliary components, skeleton arrived).
//  2b. Entry state at install: an installed box reads no local pair copy —
//     its X_NR and its partner near blocks are shipped in full and its kernel
//     near blocks are rebuilt from the kernel (see the payload note below).
//     Hence all copies of a batch install in one round, in any order.  (Until
//     ledger 30 the row diet rebuilt them from the pair copies, which had to
//     be at the source's entry state: R2 before the batch, R1 wave by wave.)
//  3. Determinism: numbers never depend on timing.  Thread splits: every
//     K-way split works in fixed chunks (split_columns_chunked, the chunked
//     right-inverse) so numbers are bitwise independent of K — K depends on
//     the wave size, hence on timing, and OpenBLAS rounds differently per
//     call shape (2026-08-27: fast != serial until this held).  Pairs: R1/R2/R3 filters
//     (pair_watermark).  Sketch: sketch_eliminated_filter — during a batch
//     of colour c only own-component or lower-colour boxes count as
//     eliminated.  Batch membership is free to vary.
//  4. Payload framing: [morton][kind][size][data]; COMPACT carries all near
//     blocks and is shared by every 1-hop reader; receivers install and mark
//     only what they were sent.
//  5. Reaction latency (design §9): arrivals are serviced at every wave
//     boundary (owner_schedule_run_boundary polls before each wave); the
//     runnable set is drained one colour at a time, lowest first; a
//     runnable component of lower colour than the batch in progress
//     preempts it after the current wave (batches are resumable).
// ---------------------------------------------------------------------------

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <climits>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <mpi.h>
#include <omp.h>

#include "owner_exchange.hpp"

namespace fmm {

/// Component of the box currently being sketched on this thread (-1: none);
/// read by the level's sketch_eliminated_filter during a batch.
inline thread_local int32_t owner_current_sketch_comp = -1;

// Payload policy (measured 2026-08-27, design §13 items 13–14): receivers
// always recompute temp2 (split over K threads per box) and always rebuild
// the REPLACEd near blocks; shipping either was slower. No knobs.

inline bool owner_serial_enabled() {
    return ca_owner_serial_enabled();
}

// ===========================================================================
// Compact eliminated state: what a receiver needs from a source box.
// ===========================================================================
// Fields (no geometry — the receiver's copy already has it):
//   skeleton/redundant indices, use_full_set, X_RR pivots, T, X_RR, X_RR_full,
//   X_RS, X_SR, Schur, X_RN, the neighbour row counts and block flags, the
//   X_NR row-diet section, then the near blocks passing `keep_neighbor` and
//   all far blocks (none under lazy).  temp2 is never shipped (recomputed).
//
// X_NR and the near blocks are shipped in full (2026-09-08, ledger 30).  The
// row diet of design §16 (regenerate X_NR rows and rebuild near blocks from
// the receiver's own pair copies) cut the bytes by half, but it made every
// install depend on the receiver's pair copies being at the source's entry
// state (invariant 2b), which forces copies to be installed wave by wave —
// 8 thin rounds per colour that, at 64 ranks, cost more than the bytes
// (ledger 29).  With everything shipped an installed box reads nothing local:
//   * all X_NR rows are shipped, as the raw column-major matrix, straight
//     out of the owner's box storage (no packing copy: the send reads the
//     matrix in place, owner_schedule_push);
//   * near blocks of BLOCK_PARTNER slots the receiver needs are shipped
//     (BLOCK_SHIPPED); BLOCK_KERNEL slots are still rebuilt from the kernel
//     + update, which reads no pair copy and is bitwise the owner's;
//   * temp2 is recomputed from the shipped X_NR and X_RR.
// So all copies of a batch are installed in ONE parallel round; their passes
// still run wave by wave (owner_schedule_batch_wave), and a pass never
// touches a copy installed after its source (the copy already carries that
// source's effect — factorization.hpp, newer-candidate rules).  The row-diet
// reader (owner_regenerate_x_nr) and the per-slot payload format are kept
// for reference; the snapshot never selects them.

using dataflow::PayloadKind;

namespace owner_detail {

template<typename T>
inline size_t vec_size(const std::vector<T>& v) { return sizeof(size_t) + v.size() * sizeof(T); }
template<typename T>
inline char* vec_put(const std::vector<T>& v, char* p) {
    const size_t n = v.size();
    std::memcpy(p, &n, sizeof(size_t));
    p += sizeof(size_t);
    if (n) { std::memcpy(p, v.data(), n * sizeof(T)); p += n * sizeof(T); }
    return p;
}
template<typename T>
inline const char* vec_get(std::vector<T>& v, const char* p) {
    size_t n = 0;
    std::memcpy(&n, p, sizeof(size_t));
    p += sizeof(size_t);
    v.resize(n);
    if (n) { std::memcpy(v.data(), p, n * sizeof(T)); p += n * sizeof(T); }
    return p;
}

}  // namespace owner_detail

/// Per-neighbour block provenance shipped with a box (index = one_hop slot):
///   0  no stored block on the owner before its elimination (receiver
///      rebuilds from the kernel + update), 1 stored (receiver rebuilds from
///      the partner's copy + update), 2 shipped in this payload.
enum : int8_t { BLOCK_KERNEL = 0, BLOCK_PARTNER = 1, BLOCK_SHIPPED = 2 };

/// Row range [off, off + n) of X_NR taken by 1-hop slot `slot` (rows are
/// stored slot by slot, in one_hop order, as compute_and_modify assembles them).
template<typename CoordType, typename DataType>
inline int64_t x_nr_slot_row_offset(const BoxData<CoordType, DataType>& box, size_t slot) {
    int64_t off = 0;
    for (size_t i = 0; i < slot; ++i) off += box.deferred_xnn_neighbor_point_counts[i];
    return off;
}

/// X_NR layout byte in the payload.
enum : int8_t { XNR_ROW_BLOCKS = 0, XNR_FULL_RAW = 1 };

/// Every row block shipped: the matrix goes as one raw column-major block.
inline bool xnr_all_shipped(const std::vector<int8_t>& rows_shipped) {
    for (int8_t v : rows_shipped) if (!v) return false;
    return true;
}

/// Bytes of the X_NR data part of a payload section.
template<typename CoordType, typename DataType>
size_t compact_state_xnr_bytes(const BoxData<CoordType, DataType>& box, const std::vector<int8_t>& rows_shipped) {
    if (!box.X_NR.is_allocated()) return 0;
    if (xnr_all_shipped(rows_shipped)) return static_cast<size_t>(box.X_NR.rows * box.X_NR.cols) * sizeof(DataType);
    size_t s = 0;
    for (size_t i = 0; i < rows_shipped.size(); ++i)
        if (rows_shipped[i]) s += static_cast<size_t>(box.deferred_xnn_neighbor_point_counts[i] * box.X_NR.cols) * sizeof(DataType);
    return s;
}

/// Bytes of a payload section before its X_NR data part (the zero-copy send
/// splits the section there: prefix and suffix from the packed buffer, the
/// matrix from the box).
template<typename CoordType, typename DataType>
size_t compact_state_xnr_prefix(const BoxData<CoordType, DataType>& box, const std::vector<int8_t>& block_flags,
                                const std::vector<int8_t>& rows_shipped) {
    using namespace owner_detail;
    size_t s = 0;
    s += vec_size(box.skeleton_indices) + vec_size(box.redundant_indices) + vec_size(box.use_full_set) +
         vec_size(box.X_RR_pivots);
    s += get_serialized_size(box.interpolation_matrix) + get_serialized_size(box.X_RR) +
         get_serialized_size(box.X_RR_full) + get_serialized_size(box.X_RS) + get_serialized_size(box.X_RS_entry) +
         get_serialized_size(box.X_SR) + get_serialized_size(box.schur_complement) + get_serialized_size(box.X_RN);
    s += vec_size(box.deferred_xnn_neighbor_point_counts) + vec_size(block_flags) + vec_size(rows_shipped);
    s += 2 * sizeof(int64_t) + sizeof(int8_t);
    return s;
}

/// The X_NR matrix as one contiguous block in the box's storage (lda == rows),
/// which the send can read in place.
template<typename CoordType, typename DataType>
inline bool xnr_contiguous(const BoxData<CoordType, DataType>& box) {
    return box.X_NR.is_allocated() && box.X_NR.lda == box.X_NR.rows &&
           box.X_NR.data.size() >= static_cast<size_t>(box.X_NR.rows * box.X_NR.cols);
}

template<typename CoordType, typename DataType>
size_t compact_state_size(const BoxData<CoordType, DataType>& box, const std::function<bool(int64_t)>& keep_neighbor,
                          const std::vector<int8_t>& block_flags, const std::vector<int8_t>& rows_shipped) {
    using namespace owner_detail;
    size_t s = 0;
    s += vec_size(box.skeleton_indices) + vec_size(box.redundant_indices) + vec_size(box.use_full_set) +
         vec_size(box.X_RR_pivots);
    s += get_serialized_size(box.interpolation_matrix) + get_serialized_size(box.X_RR) +
         get_serialized_size(box.X_RR_full) + get_serialized_size(box.X_RS) + get_serialized_size(box.X_RS_entry) +
         get_serialized_size(box.X_SR) + get_serialized_size(box.schur_complement) + get_serialized_size(box.X_RN);
    s += vec_size(box.deferred_xnn_neighbor_point_counts) + vec_size(block_flags) + vec_size(rows_shipped);
    // X_NR section: dims, layout byte, then the matrix (full raw) or the
    // shipped row blocks (row diet)
    s += 2 * sizeof(int64_t) + sizeof(int8_t);
    s += compact_state_xnr_bytes(box, rows_shipped);
    s += sizeof(size_t);
    for (const auto& b : box.near_field_modified_interactions)
        if (keep_neighbor(b.neighbor_morton)) s += get_serialized_size(b);
    s += sizeof(size_t);
    for (const auto& b : box.far_field_modified_interactions) s += get_serialized_size(b);
    return s;
}

/// Serialize a box's section.  When every X_NR row block is shipped and
/// `xnr_in_place` is set, the matrix bytes are NOT written: the section is
/// packed without them (the caller sends the matrix straight from the box
/// between the packed prefix and suffix) and `p` is not advanced over them.
template<typename CoordType, typename DataType>
char* serialize_compact_state(const BoxData<CoordType, DataType>& box, const std::function<bool(int64_t)>& keep_neighbor,
                              char* p, const std::vector<int8_t>& block_flags, const std::vector<int8_t>& rows_shipped,
                              bool xnr_in_place = false) {
    using namespace owner_detail;
    if (block_flags.size() != box.one_hop.size() || rows_shipped.size() != box.one_hop.size() ||
        box.deferred_xnn_neighbor_point_counts.size() != box.one_hop.size())
        throw std::runtime_error("serialize_compact_state: per-slot vectors do not match one_hop for box " +
                                 std::to_string(box.morton_index));
    p = vec_put(box.skeleton_indices, p);
    p = vec_put(box.redundant_indices, p);
    p = vec_put(box.use_full_set, p);
    p = vec_put(box.X_RR_pivots, p);
    p = serialize(box.interpolation_matrix, p);
    p = serialize(box.X_RR, p);
    p = serialize(box.X_RR_full, p);
    p = serialize(box.X_RS, p);
    p = serialize(box.X_RS_entry, p);
    p = serialize(box.X_SR, p);
    p = serialize(box.schur_complement, p);
    p = serialize(box.X_RN, p);
    p = vec_put(box.deferred_xnn_neighbor_point_counts, p);
    p = vec_put(block_flags, p);
    p = vec_put(rows_shipped, p);
    {
        const int64_t rows = box.X_NR.is_allocated() ? box.X_NR.rows : 0;
        const int64_t cols = box.X_NR.is_allocated() ? box.X_NR.cols : 0;
        std::memcpy(p, &rows, sizeof(int64_t));
        p += sizeof(int64_t);
        std::memcpy(p, &cols, sizeof(int64_t));
        p += sizeof(int64_t);
        const int8_t layout = xnr_all_shipped(rows_shipped) ? XNR_FULL_RAW : XNR_ROW_BLOCKS;
        std::memcpy(p, &layout, sizeof(int8_t));
        p += sizeof(int8_t);
        if (layout == XNR_FULL_RAW) {
            if (xnr_in_place) {
                // bytes come from the box itself (owner_schedule_push); nothing written here
            } else if (rows > 0 && cols > 0) {
                for (int64_t j = 0; j < cols; ++j) {   // column by column: contiguous whatever lda is
                    std::memcpy(p, box.X_NR.data.data() + j * box.X_NR.lda, static_cast<size_t>(rows) * sizeof(DataType));
                    p += static_cast<size_t>(rows) * sizeof(DataType);
                }
            }
        } else {
            int64_t off = 0;
            for (size_t i = 0; i < rows_shipped.size(); ++i) {
                const int64_t n_i = box.deferred_xnn_neighbor_point_counts[i];
                if (rows_shipped[i] && n_i > 0) {
                    if (off + n_i > rows) throw std::runtime_error("serialize_compact_state: X_NR rows short of the slot counts");
                    for (int64_t j = 0; j < cols; ++j) {   // column-major: n_i contiguous values per column
                        std::memcpy(p, box.X_NR.data.data() + j * box.X_NR.lda + off, static_cast<size_t>(n_i) * sizeof(DataType));
                        p += static_cast<size_t>(n_i) * sizeof(DataType);
                    }
                }
                off += n_i;
            }
        }
    }
    size_t n = 0;
    for (const auto& b : box.near_field_modified_interactions) if (keep_neighbor(b.neighbor_morton)) ++n;
    std::memcpy(p, &n, sizeof(size_t));
    p += sizeof(size_t);
    for (const auto& b : box.near_field_modified_interactions)
        if (keep_neighbor(b.neighbor_morton)) p = serialize(b, p);
    n = box.far_field_modified_interactions.size();
    std::memcpy(p, &n, sizeof(size_t));
    p += sizeof(size_t);
    for (const auto& b : box.far_field_modified_interactions) p = serialize(b, p);
    return p;
}

/// Inverse of serialize_compact_state.  X_NR comes back allocated (rows x
/// cols) with the shipped row blocks in place and zeros elsewhere; the
/// caller regenerates the other rows (owner_regenerate_x_nr) using the
/// returned block flags and shipped-row mask.
template<typename CoordType, typename DataType>
const char* deserialize_compact_state(BoxData<CoordType, DataType>& box, const char* p,
                                      std::vector<int8_t>& block_flags_out, std::vector<int8_t>& rows_shipped_out) {
    using namespace owner_detail;
    p = vec_get(box.skeleton_indices, p);
    p = vec_get(box.redundant_indices, p);
    p = vec_get(box.use_full_set, p);
    p = vec_get(box.X_RR_pivots, p);
    p = deserialize(box.interpolation_matrix, p);
    p = deserialize(box.X_RR, p);
    p = deserialize(box.X_RR_full, p);
    p = deserialize(box.X_RS, p);
    p = deserialize(box.X_RS_entry, p);
    p = deserialize(box.X_SR, p);
    p = deserialize(box.schur_complement, p);
    p = deserialize(box.X_RN, p);
    p = vec_get(box.deferred_xnn_neighbor_point_counts, p);
    p = vec_get(block_flags_out, p);
    p = vec_get(rows_shipped_out, p);
    box.deferred_xnn_temp2.clear();   // never shipped; rebuilt by the receiver
    {
        int64_t rows = 0, cols = 0;
        std::memcpy(&rows, p, sizeof(int64_t));
        p += sizeof(int64_t);
        std::memcpy(&cols, p, sizeof(int64_t));
        p += sizeof(int64_t);
        int8_t layout = XNR_ROW_BLOCKS;
        std::memcpy(&layout, p, sizeof(int8_t));
        p += sizeof(int8_t);
        box.X_NR.allocate(rows, cols, MatrixStorage<DataType>::FULL);   // lda == rows
        if (layout == XNR_FULL_RAW) {
            const size_t bytes = static_cast<size_t>(rows * cols) * sizeof(DataType);
            if (bytes) std::memcpy(box.X_NR.data.data(), p, bytes);
            p += bytes;
        } else {
            // allocate() keeps old contents when not shrinking; the regeneration
            // relies on unshipped rows (and their A_NS rows) starting from zero
            std::fill(box.X_NR.data.begin(), box.X_NR.data.end(), DataType{0});
            int64_t off = 0;
            for (size_t i = 0; i < rows_shipped_out.size(); ++i) {
                const int64_t n_i = box.deferred_xnn_neighbor_point_counts[i];
                if (rows_shipped_out[i] && n_i > 0) {
                    if (off + n_i > rows) throw std::runtime_error("deserialize_compact_state: X_NR rows short of the slot counts");
                    for (int64_t j = 0; j < cols; ++j) {
                        std::memcpy(box.X_NR.data.data() + j * box.X_NR.lda + off, p, static_cast<size_t>(n_i) * sizeof(DataType));
                        p += static_cast<size_t>(n_i) * sizeof(DataType);
                    }
                }
                off += n_i;
            }
        }
    }
    auto install_blocks = [&](std::vector<ModifiedBlock<DataType>>& blocks, std::unordered_map<int64_t, int64_t>& map) {
        size_t n = 0;
        std::memcpy(&n, p, sizeof(size_t));
        p += sizeof(size_t);
        for (size_t i = 0; i < n; ++i) {
            ModifiedBlock<DataType> tmp;
            p = deserialize(tmp, p);
            auto it = map.find(tmp.neighbor_morton);
            if (it == map.end()) {
                const int64_t idx = static_cast<int64_t>(blocks.size());
                blocks.push_back(std::move(tmp));
                map[blocks.back().neighbor_morton] = idx;
            } else {
                blocks[static_cast<size_t>(it->second)] = std::move(tmp);
            }
        }
    };
    install_blocks(box.near_field_modified_interactions, box.near_field_interaction_map);
    box.num_near_field_interactions = static_cast<int64_t>(box.near_field_modified_interactions.size());
    install_blocks(box.far_field_modified_interactions, box.far_field_interaction_map);
    box.num_far_field_interactions = static_cast<int64_t>(box.far_field_modified_interactions.size());
    return p;
}

template<typename CoordType, typename DataType>
size_t skeleton_section_size(const BoxData<CoordType, DataType>& box) {
    using namespace owner_detail;
    return vec_size(box.skeleton_indices) + vec_size(box.redundant_indices) + vec_size(box.use_full_set);
}
template<typename CoordType, typename DataType>
char* serialize_skeleton_section(const BoxData<CoordType, DataType>& box, char* p) {
    using namespace owner_detail;
    p = vec_put(box.skeleton_indices, p);
    p = vec_put(box.redundant_indices, p);
    p = vec_put(box.use_full_set, p);
    return p;
}
template<typename CoordType, typename DataType>
const char* deserialize_skeleton_section(BoxData<CoordType, DataType>& box, const char* p) {
    using namespace owner_detail;
    p = vec_get(box.skeleton_indices, p);
    p = vec_get(box.redundant_indices, p);
    p = vec_get(box.use_full_set, p);
    return p;
}

// ===========================================================================
// Receiver-side regeneration of a shipped source's X_NR (row diet, design §16)
// ===========================================================================
/// Fill the rows of `box.X_NR` that were not shipped, exactly as
/// compute_and_modify assembled them on the owner: per 1-hop slot the block
/// A(G active rows, X points) — from the kernel for BLOCK_KERNEL slots
/// (bitwise identical to the owner's), from this rank's copy of the pair
/// (G,X), transposed, for slots whose neighbour is of interest here (equal to
/// the owner's entry block up to rounding, as for the block rebuild) — split
/// into its S columns (A_NS) and R columns (A_NR), then one GEMM
///   X_NR = A_NR - A_NS * T
/// over all rows with the owner's call shape.  Shipped row blocks are already
/// in X_NR and have zero A_NS rows, so they pass through unchanged.
template<typename CoordType, typename DataType, typename KernelType>
void owner_regenerate_x_nr(TreeLevel<CoordType, DataType>& level, KernelType* kernel,
                           BoxData<CoordType, DataType>& box, const std::vector<int8_t>& flags,
                           const std::vector<int8_t>& rows_shipped, int dimension, int split, double* t_ms = nullptr) {
    using rg_clock = std::chrono::high_resolution_clock;
    const auto t0 = rg_clock::now();
    const int64_t k = static_cast<int64_t>(box.skeleton_indices.size());
    const int64_t r = static_cast<int64_t>(box.redundant_indices.size());
    const int64_t n = box.num_points;
    const int64_t total = box.X_NR.is_allocated() ? box.X_NR.rows : 0;
    if (k == 0 || r == 0 || total == 0) return;
    const std::string who = "owner_regenerate_x_nr: box " + std::to_string(box.morton_index) + ": ";
    if (box.X_NR.cols != r) throw std::runtime_error(who + "X_NR cols != r");
    const auto& counts = box.deferred_xnn_neighbor_point_counts;
    if (counts.size() != box.one_hop.size() || flags.size() != box.one_hop.size() || rows_shipped.size() != box.one_hop.size())
        throw std::runtime_error(who + "per-slot vectors do not match one_hop");
    const auto& T = box.interpolation_matrix;
    if (!T.is_allocated() || T.rows != k || T.cols != r) throw std::runtime_error(who + "interpolation matrix missing");
    bool any = false;
    for (size_t i = 0; i < counts.size(); ++i) if (!rows_shipped[i] && counts[i] > 0) any = true;
    if (!any) return;   // every row block was shipped

    // Per-thread scratch with retained capacity: no allocation churn per box.
    static thread_local std::vector<DataType> A_NS_all, A_NB;
    static thread_local std::vector<int64_t> selected_indices;
    A_NS_all.resize(static_cast<size_t>(total * k));
    // Rows of shipped (or empty) slots take no part in the GEMM below: their
    // A_NS rows must be zero (X_NR already holds their shipped values).
    {
        int64_t o = 0;
        for (size_t i = 0; i < box.one_hop.size(); ++i) {
            const int64_t n_g = counts[i];
            if (n_g > 0 && rows_shipped[i])
                for (int64_t j = 0; j < k; ++j)
                    std::fill_n(A_NS_all.begin() + (o + j * total), n_g, DataType{0});
            o += n_g;
        }
    }
    int64_t off = 0;
    for (size_t i = 0; i < box.one_hop.size(); ++i) {
        const int64_t g = box.one_hop[i];
        const int64_t n_g = counts[i];
        if (n_g == 0 || rows_shipped[i]) { off += n_g; continue; }
        BoxData<CoordType, DataType>* gbox = level.find_local_box(g);
        if (gbox == nullptr) gbox = level.find_ghost_box(g);
        if (flags[i] == BLOCK_KERNEL) {
            // pure kernel rows: G's active points (full or skeleton, told by the
            // slot count) against all of X's points, as the owner evaluated them
            const int64_t* gindices = nullptr;
            if (gbox != nullptr) {
                if (n_g == gbox->num_points) {
                    gindices = gbox->point_indices.data();
                } else if (n_g == static_cast<int64_t>(gbox->skeleton_indices.size())) {
                    selected_indices.resize(static_cast<size_t>(n_g));
                    for (int64_t q = 0; q < n_g; ++q) {
                        const int64_t src = gbox->skeleton_indices[static_cast<size_t>(q)];
                        selected_indices[static_cast<size_t>(q)] =
                            gbox->point_indices[static_cast<size_t>(src)];
                    }
                    gindices = selected_indices.data();
                } else {
                    throw std::runtime_error(who + "slot count of neighbour " + std::to_string(g) + " matches neither full nor skeleton");
                }
            } else {
                auto it = level.assisting_box_points_for_kernel_evaluation.find(g);
                if (it == level.assisting_box_points_for_kernel_evaluation.end())
                    throw std::runtime_error(who + "neighbour " + std::to_string(g) + " not held (no shell)");
                const auto& ab = level.assisting_boxes[static_cast<size_t>(it->second)];
                if (static_cast<int64_t>(ab.indices.size()) != n_g)
                    throw std::runtime_error(who + "assisting neighbour " + std::to_string(g) + " point count != slot count");
                gindices = ab.indices.data();
            }
            A_NB.resize(static_cast<size_t>(n_g * n));
            evaluate_block_by_index_split<DataType>(
                kernel, gindices, n_g,
                box.point_indices.data(), n,
                A_NB.data(), n_g, split);
            for (int64_t j = 0; j < k; ++j) {
                const int64_t src = box.skeleton_indices[static_cast<size_t>(j)];
                for (int64_t q = 0; q < n_g; ++q)
                    A_NS_all[static_cast<size_t>((off + q) + j * total)] = A_NB[static_cast<size_t>(q + src * n_g)];
            }
            for (int64_t j = 0; j < r; ++j) {
                const int64_t src = box.redundant_indices[static_cast<size_t>(j)];
                for (int64_t q = 0; q < n_g; ++q)
                    box.X_NR.data[static_cast<size_t>((off + q) + j * box.X_NR.lda)] = A_NB[static_cast<size_t>(q + src * n_g)];
            }
        } else {
            // modified rows: this rank's copy of the pair on G's side, stored as
            // (X points x G points), is the owner's entry block transposed
            // (R2/R1 brought it to that state).  Read it once, straight into
            // the S columns (A_NS) and R columns (X_NR) of X's row block — no
            // transposed temporary, no second copy (ledger 20e: four memory
            // passes per block before).  Rows are taken in groups of four so
            // the strided writes fill whole cache lines.
            if (gbox == nullptr) throw std::runtime_error(who + "neighbour " + std::to_string(g) + " (partner copy holder) not held");
            auto it = gbox->near_field_interaction_map.find(box.morton_index);
            if (it == gbox->near_field_interaction_map.end())
                throw std::runtime_error(who + "partner copy on neighbour " + std::to_string(g) + " missing");
            const auto& gblk =
                gbox->near_field_modified_interactions[
                    static_cast<size_t>(it->second)];
            if (!gblk.a_ns_is_allocated()) throw std::runtime_error(who + "partner copy not allocated");
            if (gblk.a_ns_rows() != n)
                throw std::runtime_error(who + "partner copy of neighbour " + std::to_string(g) + " has " +
                                         std::to_string(gblk.a_ns_rows()) + " X rows, expected " + std::to_string(n));
            // which columns of G's copy (G's points) form this slot's rows
            const int64_t* gsel = nullptr;
            if (n_g == gblk.a_ns_cols()) {
                gsel = nullptr;   // identity
            } else if (n_g == static_cast<int64_t>(gbox->skeleton_indices.size()) && gblk.a_ns_cols() == gbox->num_points) {
                gsel = gbox->skeleton_indices.data();   // G eliminated after the copy was made: skeleton rows
            } else {
                throw std::runtime_error(who + "partner copy of neighbour " + std::to_string(g) + " has " +
                                         std::to_string(gblk.a_ns_cols()) + " G columns; slot count " + std::to_string(n_g) +
                                         ", skeleton " + std::to_string(gbox->skeleton_indices.size()));
            }
            const int64_t* skel = box.skeleton_indices.data();
            const int64_t* red = box.redundant_indices.data();
            DataType* xnr = box.X_NR.data.data();
            const int64_t ldx = box.X_NR.lda;
            for (int64_t q0 = 0; q0 < n_g; q0 += 4) {
                const int64_t qn = std::min<int64_t>(4, n_g - q0);
                for (int64_t j = 0; j < k; ++j) {
                    DataType* dst = A_NS_all.data() + (off + q0) + j * total;
                    const int64_t sj = skel[j];
                    for (int64_t t = 0; t < qn; ++t) {
                        const int64_t gc = gsel ? gsel[q0 + t] : (q0 + t);
                        dst[t] = gblk.a_ns(sj, gc);
                    }
                }
                for (int64_t j = 0; j < r; ++j) {
                    DataType* dst = xnr + (off + q0) + j * ldx;
                    const int64_t rj = red[j];
                    for (int64_t t = 0; t < qn; ++t) {
                        const int64_t gc = gsel ? gsel[q0 + t] : (q0 + t);
                        dst[t] = gblk.a_ns(rj, gc);
                    }
                }
            }
        }
        off += n_g;
    }
    if (off != total) throw std::runtime_error(who + "slot counts do not sum to X_NR rows");
    {   // X_NR = A_NR - A_NS * T (the owner's GEMM: same shape, same rounding)
        int M = static_cast<int>(total), N = static_cast<int>(r), K = static_cast<int>(k);
        int ldt = static_cast<int>(T.lda), ldx = static_cast<int>(box.X_NR.lda);
        DataType alpha = -1.0, beta = 1.0;
        gemm_("N", "N", &M, &N, &K, &alpha, A_NS_all.data(), &M, T.data.data(), &ldt, &beta, box.X_NR.data.data(), &ldx);
    }
    if (t_ms) *t_ms += std::chrono::duration<double, std::milli>(rg_clock::now() - t0).count();
}

// ===========================================================================
// Receiver-side rebuild of a shipped source (payload diet)
// ===========================================================================
// The payload carries X_NR (pre-finalize; only the row blocks the receiver
// cannot regenerate, see owner_regenerate_x_nr) but not temp2, and no near
// blocks except those the home rank cannot rebuild.  The receiver redoes
// step two's temp2 solve and step 5 with the owner's routine and operation
// order:
//   temp2 = -(X_NR * X_RR^{-1})            (apply_right_inverse_in_place)
//   update = temp2 * X_RS_entry             (one GEMM, all neighbours)
//   block(X,G) = base + update_G, base = slice(transpose(G's copy)) for a
//   pair that had a stored block, kernel(G active pts, X pts)[:, S_X] otherwise.
// Only neighbours of interest here are rebuilt (others are never read).

/// `rebuild_all`: the home rank needs the box's blocks toward every
/// neighbour (the transition reads them), not only toward its interest boxes.
template<typename CoordType, typename DataType, typename KernelType>
void owner_rebuild_shipped_source(const std::function<bool(int64_t)>& interest, TreeLevel<CoordType, DataType>& level,
                                  KernelType* kernel, FactorizationMethod method, BoxData<CoordType, DataType>& box,
                                  const std::vector<int8_t>& flags, int dimension, bool rebuild_all,
                                  double* t_temp2_ms = nullptr, double* t_blocks_ms = nullptr, int split = 1) {
    using rb_clock = std::chrono::high_resolution_clock;
    const auto rb_t0 = rb_clock::now();
    const int64_t k = static_cast<int64_t>(box.skeleton_indices.size());
    const int64_t r = static_cast<int64_t>(box.redundant_indices.size());
    if (k == 0 || r == 0 || !box.X_NR.is_allocated() || box.deferred_xnn_neighbor_point_counts.empty()) return;
    const int64_t total = box.X_NR.rows;
    if (box.X_NR.cols != r) throw std::runtime_error("owner_rebuild_shipped_source: X_NR cols != r");
    if (flags.size() != box.one_hop.size())
        throw std::runtime_error("owner_rebuild_shipped_source: block flags missing for box " + std::to_string(box.morton_index));

    if (box.deferred_xnn_temp2.empty() && total > 0) {
        std::vector<DataType> temp2(box.X_NR.data.begin(), box.X_NR.data.end());
        apply_right_inverse_in_place(box.X_RR, box.X_RR_pivots, method, temp2, total, r,
                                     "owner_rebuild_shipped_source temp2", split);
        for (auto& v : temp2) v = -v;
        box.deferred_xnn_temp2 = std::move(temp2);
    }
    const auto rb_t1 = rb_clock::now();
    if (t_temp2_ms) *t_temp2_ms += std::chrono::duration<double, std::milli>(rb_t1 - rb_t0).count();
    if (!box.X_RS_entry.is_allocated()) throw std::runtime_error("owner_rebuild_shipped_source: X_RS_entry missing");

    // update = temp2 * X_RS_entry, in fixed row chunks over K threads (the
    // chunk, not K, fixes the per-element rounding — see apply_right_inverse)
    std::vector<DataType> update(static_cast<size_t>(total * k));
    {
        constexpr int64_t CHUNK = 512;
        const int64_t nchunks = (total + CHUNK - 1) / CHUNK;
        auto rows_gemm = [&](int64_t r0, int64_t r1) {
            int m = static_cast<int>(r1 - r0), n = static_cast<int>(k), kk = static_cast<int>(r), ld = static_cast<int>(total);
            DataType alpha = 1.0, beta = 0.0;
            gemm_("N", "N", &m, &n, &kk, &alpha, box.deferred_xnn_temp2.data() + r0, &ld, box.X_RS_entry.data.data(), &kk,
                  &beta, update.data() + r0, &ld);
        };
        if (split > 1 && nchunks > 1) {
            const int K = static_cast<int>(std::min<int64_t>(split, nchunks));
            #pragma omp taskloop num_tasks(K) default(shared)
            for (int64_t c = 0; c < nchunks; ++c) rows_gemm(c * CHUNK, std::min<int64_t>(total, (c + 1) * CHUNK));
        } else {
            for (int64_t c = 0; c < nchunks; ++c) rows_gemm(c * CHUNK, std::min<int64_t>(total, (c + 1) * CHUNK));
        }
    }

    int64_t current_row = 0;
    for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
        const int64_t g = box.one_hop[idx];
        const int64_t n_g = box.deferred_xnn_neighbor_point_counts[idx];
        const int8_t flag = flags[idx];
        if (n_g == 0 || flag == BLOCK_SHIPPED || (!rebuild_all && !interest(g))) { current_row += n_g; continue; }
        BoxData<CoordType, DataType>* gbox = level.find_local_box(g);
        if (gbox == nullptr) gbox = level.find_ghost_box(g);
        if (gbox == nullptr)
            throw std::runtime_error("owner_rebuild_shipped_source: neighbour " + std::to_string(g) + " of box " +
                                     std::to_string(box.morton_index) + " not held");
        std::vector<DataType> base;
        if (flag == BLOCK_PARTNER) {
            auto it = gbox->near_field_interaction_map.find(box.morton_index);
            if (it == gbox->near_field_interaction_map.end())
                throw std::runtime_error("owner_rebuild_shipped_source: partner copy missing: box " +
                                         std::to_string(box.morton_index) + " neighbour " + std::to_string(g));
            const auto& gblk = gbox->near_field_modified_interactions[static_cast<size_t>(it->second)];
            if (!gblk.a_ns_is_allocated())
                throw std::runtime_error("owner_rebuild_shipped_source: partner block not allocated");
            ModifiedBlock<DataType> tmp;   // X's own copy is the transpose of G's copy
            tmp.neighbor_morton = g;
            const int64_t transposed_rows = gblk.a_ns_cols();
            const int64_t transposed_cols = gblk.a_ns_rows();
            std::vector<DataType> transposed(
                static_cast<size_t>(transposed_rows * transposed_cols));
            for (int64_t col = 0; col < gblk.a_ns_cols(); ++col) {
                for (int64_t row = 0; row < gblk.a_ns_rows(); ++row) {
                    transposed[static_cast<size_t>(col + row * transposed_rows)] =
                        gblk.a_ns(row, col);
                }
            }
            tmp.set_a_ns_owned(
                transposed_rows, transposed_cols,
                std::move(transposed), MatrixStorage<DataType>::FULL);
            base = slice_modified_block_both_directions<CoordType, DataType>(tmp, level, g, box.skeleton_indices, true, n_g);
        } else {
            const bool use_skeleton = (n_g != gbox->num_points);
            if (use_skeleton && n_g != static_cast<int64_t>(gbox->skeleton_indices.size()))
                throw std::runtime_error("owner_rebuild_shipped_source: slot count of neighbour " + std::to_string(g) +
                                         " matches neither full nor skeleton");
            std::vector<int64_t> indices(static_cast<size_t>(n_g));
            for (int64_t i = 0; i < n_g; ++i) {
                const int64_t src = use_skeleton ? gbox->skeleton_indices[static_cast<size_t>(i)] : i;
                indices[static_cast<size_t>(i)] =
                    gbox->point_indices[static_cast<size_t>(src)];
            }
            std::vector<DataType> A_NB(static_cast<size_t>(n_g * box.num_points));
            evaluate_block_by_index_split<DataType>(
                kernel, indices.data(), n_g,
                box.point_indices.data(), box.num_points,
                A_NB.data(), n_g, split);
            base.resize(static_cast<size_t>(n_g * k));
            for (int64_t j = 0; j < k; ++j) {
                const int64_t src_col = box.skeleton_indices[static_cast<size_t>(j)];
                for (int64_t i = 0; i < n_g; ++i)
                    base[static_cast<size_t>(i + j * n_g)] = A_NB[static_cast<size_t>(i + src_col * n_g)];
            }
        }
        if (static_cast<int64_t>(base.size()) != n_g * k)
            throw std::runtime_error("owner_rebuild_shipped_source: base size mismatch");
        for (int64_t j = 0; j < k; ++j)
            for (int64_t i = 0; i < n_g; ++i)
                base[static_cast<size_t>(i + j * n_g)] += update[static_cast<size_t>((current_row + i) + j * total)];
        auto it = box.near_field_interaction_map.find(g);
        if (it == box.near_field_interaction_map.end()) {
            ModifiedBlock<DataType> nb;
            nb.neighbor_morton = g;
            nb.set_a_ns_owned(
                n_g, k, std::move(base), MatrixStorage<DataType>::FULL);
            const int64_t nidx = static_cast<int64_t>(box.near_field_modified_interactions.size());
            box.near_field_modified_interactions.push_back(std::move(nb));
            box.near_field_interaction_map[g] = nidx;
        } else {
            box.near_field_modified_interactions[
                static_cast<size_t>(it->second)].set_a_ns_owned(
                    n_g, k, std::move(base),
                    MatrixStorage<DataType>::FULL);
        }
        current_row += n_g;
    }
    box.num_near_field_interactions = static_cast<int64_t>(box.near_field_modified_interactions.size());
    if (t_blocks_ms) *t_blocks_ms += std::chrono::duration<double, std::milli>(rb_clock::now() - rb_t1).count();
}

// ===========================================================================
// Buffer pool
// ===========================================================================
/// Process-wide pool of large byte buffers (payload send/receive buffers, the
/// reduction's parent-box buffers).  The engine allocates and frees GB-scale
/// buffers every batch; freeing them makes tcmalloc return whole spans to the
/// OS, after which every large allocation in the process — the next batch's
/// buffers, the parents received at the reduction, the next level's boxes —
/// page-faults fresh 4 KB pages under all threads at once and loses its huge
/// pages (design ledger items 19/20: +21 s in one reduction and +8 s at L3
/// for Helmholtz 192³).  Buffers taken from here keep their capacity for the
/// life of the process, so the allocator never sees the churn.  Numbers are
/// unaffected: the same bytes are written, into memory that is already mapped.
/// Set by allocator_policy_apply() when the process allocator keeps freed
/// pages resident (no release to the OS); then the pool is drained at the end
/// of each level so its memory is reusable by everything else.
inline bool owner_allocator_keeps_freed_pages = false;

class OwnerBufferPool {
public:
    static OwnerBufferPool& instance() {
        static OwnerBufferPool p;
        return p;
    }
    /// A buffer of size() == bytes.  Best fit among the free buffers whose
    /// capacity covers max(bytes, reserve_hint); otherwise the largest free
    /// one, grown; otherwise a new one.  reserve_hint lets a buffer that will
    /// grow (a batch's destination payload) start at the largest size seen.
    std::vector<char> acquire(size_t bytes, size_t reserve_hint = 0) {
        std::lock_guard<std::mutex> lock(mu_);
        const size_t want = std::max(bytes, reserve_hint);
        size_t best = free_.size(), best_cap = 0;
        for (size_t i = 0; i < free_.size(); ++i) {
            const size_t c = free_[i].capacity();
            if (c >= want && (best == free_.size() || c < best_cap)) { best = i; best_cap = c; }
        }
        if (best == free_.size())
            for (size_t i = 0; i < free_.size(); ++i)
                if (free_[i].capacity() > best_cap) { best = i; best_cap = free_[i].capacity(); }
        std::vector<char> v;
        if (best < free_.size()) {
            v = std::move(free_[best]);
            free_[best] = std::move(free_.back());
            free_.pop_back();
        }
        if (v.capacity() < want) v.reserve(want);
        v.resize(bytes);
        return v;
    }
    // (An MADV_HUGEPAGE hint on these buffers was tried on 2026-09-08, ledger
    // 33/34: the compute nodes run THP 'always' with defrag 'madvise', the
    // hint left the huge-page coverage unchanged and can only add compaction
    // stalls there, so it was removed.  The THP mode is printed in the banner.)
    void release(std::vector<char>&& v) {
        if (v.capacity() == 0) return;
        std::lock_guard<std::mutex> lock(mu_);
        max_seen_ = std::max(max_seen_, v.capacity());
        v.clear();
        free_.push_back(std::move(v));
    }
    size_t max_seen() const { std::lock_guard<std::mutex> lock(mu_); return max_seen_; }
    /// Free every pooled buffer to the allocator.  Used at the end of a
    /// level's boundary phase when the allocator keeps freed pages resident
    /// (tcmalloc with release rate 0): the buffers then serve the reduction's
    /// parent boxes instead of sitting idle (ledger 20e: rank 0 entered the
    /// reduction with 17 GB pinned here and faulted 14 GB for the parents).
    void drain() {
        std::lock_guard<std::mutex> lock(mu_);
        free_.clear();
        free_.shrink_to_fit();
    }
    size_t held_bytes() const {
        std::lock_guard<std::mutex> lock(mu_);
        size_t b = 0;
        for (const auto& v : free_) b += v.capacity();
        return b;
    }
private:
    mutable std::mutex mu_;
    std::vector<std::vector<char>> free_;
    size_t max_seen_ = 0;
};

// ===========================================================================
// State
// ===========================================================================

template<typename CoordType, typename DataType>
struct OwnerScheduleState {
    bool active = false;
    int dimension = 3;
    int num_waves = 8;
    int num_boundary_colors = 3;
    int32_t level_index = 0;
    dataflow::ProcessGrid grid;
    dataflow::ComponentGraph graph;
    dataflow::Runtime rt;
    uint32_t my_pid = 0;
    int my_rank = -1;
    MPI_Comm comm = MPI_COMM_NULL;
    ParallelTree<CoordType, DataType>* tree = nullptr;
    std::unordered_map<int, int> pid_to_rank;
    std::unordered_set<int64_t> interest_boxes;   ///< local boxes + boxes of owned components
    std::vector<int32_t> comp_seq_color;          ///< colour index of each component
    /// Wave of every boundary box (component-aware partition, see wave_of);
    /// identical on every rank because it is derived from the shared graph.
    std::unordered_map<int64_t, int32_t> box_wave;
    bool serial = false;
    bool sync_colors = false;                     ///< this level ran the synchronous colour batches (coarse level, fast schedule)
    bool coarse_level = false;                    ///< fewer than 8 home boundary boxes per thread: per-colour balancing

    struct OutMsg {
        std::vector<char> buf;                                  ///< packed part of the payload
        std::vector<int64_t> header;
        std::vector<MPI_Request> reqs;                          ///< header + one per wire part
    };
    std::list<OutMsg> outbox;
    /// Arrived payloads awaiting their install event (component id -> buffer).
    std::unordered_map<int32_t, std::vector<char>> inbox;
    /// Payloads whose parts are landing in the background (owner_schedule_receive_one).
    struct PendingRecv {
        int32_t id = -1;
        uint32_t from = 0;
        size_t total = 0;
        std::vector<char> buf;
        std::vector<MPI_Request> reqs;
    };
    std::list<PendingRecv> pending_recv;
    /// Auxiliary (skeleton-only, non-shared) payloads and which are installed.
    std::unordered_map<int32_t, std::vector<char>> aux_inbox;
    std::unordered_set<int32_t> aux_installed;
    /// Sketch filter installed on the level during a batch's box region.
    std::function<bool(int64_t)> sketch_filter;
    /// Interest predicate installed on the level for the whole level (pair
    /// ownership: the interest side of a mixed pair accumulates).
    std::function<bool(int64_t)> interest_filter;
    /// Ghost boxes whose entry BLOCKS this rank needs: the remote boxes of the
    /// components it eliminates (owned-only entry-state fetch, design §7).
    std::unordered_set<int64_t> block_ghosts;

    // stats
    int64_t comps_eliminated = 0, comps_installed = 0, boxes_eliminated = 0, batches = 0, preemptions = 0;
    size_t bytes_sent = 0, bytes_recv = 0;
    double t_wait_ms = 0, t_pack_ms = 0, t_install_ms = 0, t_replay_ms = 0, t_elim_ms = 0, t_final_ms = 0, t_passes_ms = 0;
    double t_recv_ms = 0;                        ///< landing payload parts (MPI_Recv after the header), main thread
    int64_t parts_sent = 0, parts_recv = 0;      ///< wire parts (segments split at PART_BYTES)
    /// X_NR buffers of this rank's eliminated boxes whose storage the final
    /// event hands over to temp2: kept alive here until the sends that read
    /// them in place have completed (released after the outbox drain).
    std::vector<std::vector<DataType>> xnr_hold;
    double t_rebuild_temp2_ms = 0, t_rebuild_blocks_ms = 0;   ///< inside t_install_ms (summed over threads)
    double t_rebuild_xnr_ms = 0;                              ///< X_NR row regeneration, inside t_install_ms (thread-sum)
    /// X_NR bytes by row class, summed over recipients: shipped, regenerated
    /// from the kernel, regenerated from the receiver's partner copy.
    size_t xnr_bytes_shipped = 0, xnr_bytes_regen_kernel = 0, xnr_bytes_regen_partner = 0;

    static constexpr int TAG_HDR = 3301;
    static constexpr int TAG_PART = 3302;
    static constexpr size_t PART_BYTES = size_t(1) << 30;

    int32_t comp_of(int64_t m) const { return graph.comp_id_of_box(m); }
    /// Wave (elimination class inside its colour) of a boundary box: the
    /// component-aware partition computed at setup (owner_wave_partition), so
    /// that a component's boxes fall into as few classes as its own adjacency
    /// needs instead of the 8 Morton-parity classes of the replicated code.
    int32_t wave_of(int64_t m) const {
        auto it = box_wave.find(m);
        if (it == box_wave.end()) throw std::runtime_error("owner_schedule: box " + std::to_string(m) + " has no wave");
        return it->second;
    }
    int32_t seq_of(int64_t m) const {
        const int32_t c = comp_of(m);
        const int32_t ci = c >= 0 ? comp_seq_color[static_cast<size_t>(c)] : num_boundary_colors;
        return ci * num_waves + wave_of(m);
    }
    bool interest(int64_t x, int64_t y) const { return interest_boxes.count(x) != 0 || interest_boxes.count(y) != 0; }
    /// Is box n of interest to process pid (its local box, or in a component it owns)?
    bool interest_of(uint32_t pid, int64_t n) const {
        if (grid.proc_of_box(n) == pid) return true;
        const int32_t c = comp_of(n);
        return c >= 0 && graph.comps[static_cast<size_t>(c)].owner == static_cast<int32_t>(pid);
    }
    /// Recipients of box m's state (see dataflow::payload_recipients).
    std::vector<std::pair<uint32_t, PayloadKind>> recipients_of_box(int64_t m) const {
        return dataflow::payload_recipients(graph, my_pid, m);
    }
};

// ===========================================================================
// Component-aware wave partition
// ===========================================================================
//
// Inside one colour the boxes are eliminated in waves; a wave must contain
// only pairwise non-adjacent boxes (paper §3.3: Schur updates of one wave
// never touch each other's rows), and the order of the waves is what every
// rank's deferred replay sorts by, so the partition must be identical
// everywhere.  The replicated code uses the 8 Morton-parity classes of a
// full 3D grid; a boundary component needs far fewer: a 2D face sheet is
// 4-colourable, an edge line 2-colourable, so 8 classes leave most waves
// nearly empty and stretch the level over many small box regions.
//
// Here every component is partitioned on its own 26-adjacency graph, taking
// the best of three valid candidates (its Morton parity classes renumbered
// densely, and two greedy colourings in Morton order); classes are numbered
// by descending size so that all components' largest classes share wave 0.
// Components never share a wave's box region with each other (they are
// non-adjacent), so classes are local to a component and the global wave
// count is just the largest class count over all components (<= 8).
inline std::unordered_map<int64_t, int32_t> owner_wave_partition(const dataflow::ComponentGraph& graph,
                                                                 const dataflow::ProcessGrid& grid, int& num_waves) {
    std::unordered_map<int64_t, int32_t> wave;
    size_t total = 0;
    for (const auto& k : graph.comps) total += k.boxes.size();
    wave.reserve(total);
    num_waves = 1;
    const int parity_classes = grid.dimension == 2 ? 4 : 8;
    // Three candidate partitions per component; the one with the fewest
    // classes wins (ties: the flattest).  Morton parity is always valid and
    // optimal for dense clumps (blue+orange: 8) and for sheets (4) and edges
    // (2); the greedy variants can do better on irregular strips (purple).
    // A plain greedy alone can need 13 classes on a blue clump, so the parity
    // candidate is the safety net that keeps every colour at <= 8 waves.
    std::vector<int64_t> boxes;
    std::vector<int32_t> cls;   // class of boxes[i] in the candidate being built
    std::unordered_map<int64_t, size_t> index_of;
    for (const auto& k : graph.comps) {
        boxes.assign(k.boxes.begin(), k.boxes.end());
        std::sort(boxes.begin(), boxes.end());
        index_of.clear();
        for (size_t i = 0; i < boxes.size(); ++i) index_of.emplace(boxes[i], i);
        // adjacency inside the component, by index
        std::vector<std::vector<size_t>> adj(boxes.size());
        for (size_t i = 0; i < boxes.size(); ++i)
            for (uint64_t n : grid.neighbors(boxes[i])) {
                auto it = index_of.find(static_cast<int64_t>(n));
                if (it != index_of.end()) adj[i].push_back(it->second);
            }
        // candidate 0: Morton parity, classes renumbered densely by size
        std::vector<int32_t> best;
        int best_classes = 0;
        int best_max = 0;
        auto consider = [&](std::vector<int32_t>& cand) {
            int nc = 0;
            for (int32_t c : cand) nc = std::max(nc, c + 1);
            std::vector<int> size(static_cast<size_t>(nc), 0);
            for (int32_t c : cand) ++size[static_cast<size_t>(c)];
            // renumber classes by descending size (stable) so that every
            // component's largest class is wave 0: merged batches then have
            // their fullest wave first and their thinnest last
            std::vector<int> order(static_cast<size_t>(nc));
            for (int c = 0; c < nc; ++c) order[static_cast<size_t>(c)] = c;
            std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                return size[static_cast<size_t>(a)] > size[static_cast<size_t>(b)];
            });
            std::vector<int32_t> rank(static_cast<size_t>(nc));
            for (int r = 0; r < nc; ++r) rank[static_cast<size_t>(order[static_cast<size_t>(r)])] = r;
            for (int32_t& c : cand) c = rank[static_cast<size_t>(c)];
            const int mx = size.empty() ? 0 : *std::max_element(size.begin(), size.end());
            if (best.empty() || nc < best_classes || (nc == best_classes && mx < best_max)) {
                best = cand;
                best_classes = nc;
                best_max = mx;
            }
        };
        {
            cls.resize(boxes.size());
            std::vector<int32_t> seen(static_cast<size_t>(parity_classes), -1);
            int32_t next = 0;
            for (size_t i = 0; i < boxes.size(); ++i) {
                const int32_t par = static_cast<int32_t>(boxes[i] % parity_classes);
                if (seen[static_cast<size_t>(par)] < 0) seen[static_cast<size_t>(par)] = next++;
                cls[i] = seen[static_cast<size_t>(par)];
            }
            consider(cls);
        }
        // candidates 1 and 2: greedy in Morton order — first free class, and
        // least-filled free class
        for (int variant = 0; variant < 2; ++variant) {
            cls.assign(boxes.size(), -1);
            std::vector<int32_t> members;
            std::vector<char> used;
            for (size_t i = 0; i < boxes.size(); ++i) {
                used.assign(members.size(), 0);
                for (size_t j : adj[i]) if (cls[j] >= 0) used[static_cast<size_t>(cls[j])] = 1;
                int32_t pick = -1;
                for (size_t c = 0; c < members.size(); ++c) {
                    if (used[c]) continue;
                    if (pick < 0 || (variant == 1 && members[c] < members[static_cast<size_t>(pick)])) pick = static_cast<int32_t>(c);
                    if (variant == 0) break;
                }
                if (pick < 0) {
                    pick = static_cast<int32_t>(members.size());
                    members.push_back(0);
                }
                ++members[static_cast<size_t>(pick)];
                cls[i] = pick;
            }
            consider(cls);
        }
        for (size_t i = 0; i < boxes.size(); ++i) wave.emplace(boxes[i], best[i]);
        num_waves = std::max(num_waves, best_classes);
    }
    return wave;
}

// ===========================================================================
// Setup
// ===========================================================================

template<typename CoordType, typename DataType>
void owner_schedule_setup(
                          ParallelTree<CoordType, DataType>* tree,
                          int32_t level_index,
                          MPI_Comm level_comm,
                          OwnerScheduleState<CoordType, DataType>& st) {
    using namespace dataflow;
    st = OwnerScheduleState<CoordType, DataType>{};
    if (tree == nullptr) return;
    auto& level = tree->levels[static_cast<size_t>(level_index)];
    if (!level.is_process_active || level.num_active_processes <= 1 || level.my_morton_id < 0) return;
    if (lazy_far_field_mode() != LazyFarFieldMode::LAZY)
        throw std::runtime_error(
            "H2_CA_owner_component=3 requires H2_lazy_schur=1 or 2");

    st.tree = tree;
    st.dimension = tree->dimension;
    const int dimension = st.dimension;
    st.num_waves = dimension == 2 ? 4 : 8;   // replaced below by the component-aware partition
    st.num_boundary_colors = num_boundary_colors(dimension);
    st.level_index = level_index;
    st.serial = owner_serial_enabled();
    st.grid = ProcessGrid(dimension, 1u << level_index, static_cast<uint32_t>(level.num_active_processes));
    ColorOf color_of = make_ca_classifier(
        st.grid, level_index,
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_blue_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, int32_t dim, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_orange_box(m, lvl, dim, gs, lgs, off);
        },
        [](int64_t m, int32_t lvl, uint32_t gs, uint32_t lgs, const uint32_t* off) {
            return is_purple_box(m, lvl, gs, lgs, off);
        });
    st.graph = build_component_graph(st.grid, color_of);
    // Schedule and balancing objective are decided together, before the
    // assignment (ledger 33).  A coarse level — fewer than 8 local boundary
    // boxes per thread — runs the synchronous colour batches, whose level
    // time is the sum over colours of the slowest rank's work, so each colour
    // is balanced on its own.  A fine level runs arrival-driven, where a rank
    // that finishes its blue early moves straight into its green, so the
    // totals over the colours are balanced instead (a rank heavy in green is
    // made light in blue and purple).  The count of boundary boxes whose home
    // is this rank is known before the assignment and equals the owned count
    // on average.  (The level's colour lists also hold the halo's boundary
    // boxes, so they are not the measure — ledger 35.)
    {
        int64_t local_boundary = 0;
        for (const auto& k : st.graph.comps)
            for (int64_t m : k.boxes)
                if (st.grid.proc_of_box(m) == static_cast<uint32_t>(level.my_morton_id)) ++local_boundary;
        // The ownership must not depend on the schedule oracle: the serial
        // oracle runs one component per batch instead of the synchronous
        // batches, but it must own exactly what the fast schedule owns, or
        // the accumulating side of a pair (interest rule) differs and the
        // checksums with it (ledger 46).  `coarse` decides both.
        const bool coarse = local_boundary < 8 * static_cast<int64_t>(std::max(1, omp_get_max_threads()));
        st.sync_colors = !st.serial && coarse;
        st.coarse_level = coarse;
    }
    assign_owners(st.graph, st.coarse_level ? AssignPolicy::LPT_PER_COLOR : AssignPolicy::LPT_CUMULATIVE);
    {
        const std::string e = validate_graph(st.graph);
        if (!e.empty()) throw std::runtime_error("owner_schedule_setup: invalid graph: " + e);
        const std::string l = validate_against_lists(st.graph, level.blue, level.orange, level.purple, level.green);
        if (!l.empty()) throw std::runtime_error("owner_schedule_setup: graph disagrees with level lists: " + l);
    }
    st.my_pid = static_cast<uint32_t>(level.my_morton_id);
    st.comm = level_comm;
    MPI_Comm_rank(level_comm, &st.my_rank);
    int comm_size = 0;
    MPI_Comm_size(level_comm, &comm_size);
    std::vector<int> pids(static_cast<size_t>(comm_size), -1);
    const int my_pid = level.my_morton_id;
    MPI_Allgather(
        &my_pid, 1, MPI_INT,
        pids.data(), 1, MPI_INT,
        level_comm);
    for (int comm_rank = 0; comm_rank < comm_size; ++comm_rank) {
        st.pid_to_rank.emplace(pids[static_cast<size_t>(comm_rank)], comm_rank);
    }

    st.comp_seq_color.resize(st.graph.comps.size());
    const std::vector<Color> order = boundary_colors(dimension);
    for (const auto& k : st.graph.comps) {
        auto it = std::find(order.begin(), order.end(), k.color);
        st.comp_seq_color[static_cast<size_t>(k.id)] = static_cast<int32_t>(it - order.begin());
    }
    // waves: component-aware partition (replaces morton % 8); num_waves is
    // the largest class count of any component and is the same on every rank
    st.box_wave = owner_wave_partition(st.graph, st.grid, st.num_waves);

    RuntimeOptions opts;
    opts.entry_state_prefetched = true;
    opts.factors_with_generators = true;
    opts.generators_for_shared_only = true;
    opts.two_hop_dependencies = true;
    opts.aux_comps = auxiliary_components(st.graph, st.my_pid);
    st.rt.reset(&st.graph, st.my_pid, opts);
    st.rt.restamp_ready(omp_get_wtime());

    for (const auto& box : level.local_boxes) st.interest_boxes.insert(box.morton_index);
    for (int32_t id : st.rt.owned()) {
        const Component& k = st.graph.comps[static_cast<size_t>(id)];
        for (int64_t m : k.boxes) {
            st.interest_boxes.insert(m);
            if (level.find_local_box(m) == nullptr) st.block_ghosts.insert(m);
        }
    }
    st.interest_filter = [&st](int64_t m) -> bool { return st.interest_boxes.count(m) != 0; };
    level.pair_interest_filter = &st.interest_filter;
    {
        std::vector<int32_t> mine = st.rt.owned();
        mine.insert(mine.end(), st.rt.shared_not_owned().begin(), st.rt.shared_not_owned().end());
        std::set<int32_t> mine_set(mine.begin(), mine.end());
        for (int32_t id : mine)
            for (int32_t p : st.graph.comps[static_cast<size_t>(id)].preds)
                if (!mine_set.count(p))
                    throw std::runtime_error("owner_schedule_setup: predecessor " + std::to_string(p) + " of component " +
                                             std::to_string(id) + " is not shared by rank " + std::to_string(tree->mpi_rank));
        // 2-hop predecessors of OWNED components must be processed here or
        // arrive as auxiliary skeleton payloads.
        for (int32_t id : st.rt.owned())
            for (int32_t p : st.graph.comps[static_cast<size_t>(id)].preds2)
                if (!mine_set.count(p) && !st.rt.is_aux(p))
                    throw std::runtime_error("owner_schedule_setup: 2-hop predecessor " + std::to_string(p) +
                                             " of owned component " + std::to_string(id) +
                                             " is neither shared nor auxiliary on rank " + std::to_string(tree->mpi_rank));
    }
    // sketch filter: during a batch of colour c, a marked box counts as
    // eliminated for the sketch iff it is in the box's own component or of a
    // lower colour (its skeleton is then known here by the dependency rules)
    st.sketch_filter = [&st](int64_t m) -> bool {
        const int32_t c = st.comp_of(m);
        if (c < 0) return true;
        if (c == owner_current_sketch_comp) return true;
        const int32_t cur = owner_current_sketch_comp;
        const int32_t batch_color = cur >= 0 ? st.comp_seq_color[static_cast<size_t>(cur)] : st.num_boundary_colors;
        return st.comp_seq_color[static_cast<size_t>(c)] < batch_color;
    };
    st.active = true;
}

/// After the halo (stage 0 at least): every owned or shared box must be held
/// here (local or ghost shell).  Owned ghosts also need their entry blocks,
/// which the owned-only fetch requested.
template<typename CoordType, typename DataType>
void owner_schedule_check_presence(TreeLevel<CoordType, DataType>& level, OwnerScheduleState<CoordType, DataType>& st) {
    if (!st.active) return;
    for (int32_t id : st.rt.owned())
        for (int64_t m : st.graph.comps[static_cast<size_t>(id)].boxes)
            if (level.find_local_box(m) == nullptr && level.find_ghost_box(m) == nullptr)
                throw std::runtime_error("owner_schedule: owned box " + std::to_string(m) + " of component " +
                                         std::to_string(id) + " not held by rank " + std::to_string(st.my_rank));
    for (int32_t id : st.rt.shared_not_owned())
        for (int64_t m : st.graph.comps[static_cast<size_t>(id)].boxes)
            if (level.find_local_box(m) == nullptr && level.find_ghost_box(m) == nullptr)
                throw std::runtime_error("owner_schedule: shared box " + std::to_string(m) + " of component " +
                                         std::to_string(id) + " not held by rank " + std::to_string(st.my_rank));
}

template<typename CoordType, typename DataType>
std::string owner_schedule_describe(const OwnerScheduleState<CoordType, DataType>& st) {
    if (!st.active) return "owner-schedule: inactive";
    size_t naux = 0;
    for (const auto& e : st.rt.expected()) if (e.kind == dataflow::Arrival::SKELETON) ++naux;
    return std::string("owner-schedule") + (st.serial ? " (SERIAL oracle)" : "") + ": " + dataflow::describe(st.graph) +
           " | this rank owns " + std::to_string(st.rt.owned().size()) + " components, installs " +
           std::to_string(st.rt.shared_not_owned().size()) + ", auxiliary skeleton payloads " + std::to_string(naux) +
           " | entry blocks requested for " + std::to_string(st.block_ghosts.size()) + " ghosts";
}

// ===========================================================================
// Pair liveness
// ===========================================================================

/// A batch = the components of one colour processed together (eliminated or
/// installed).  `in_event[c]` marks membership; `event_color` is their colour
/// index, or num_boundary_colors for the final event.
struct OwnerEvent {
    std::vector<char> in_event;
    int32_t event_color = 0;
    bool final = false;
    bool contains(int32_t c) const { return c >= 0 && in_event[static_cast<size_t>(c)] != 0; }
};

template<typename CoordType, typename DataType>
OwnerEvent make_event(const OwnerScheduleState<CoordType, DataType>& st, const std::vector<int32_t>& comps) {
    OwnerEvent e;
    e.in_event.assign(st.graph.comps.size(), 0);
    for (int32_t c : comps) e.in_event[static_cast<size_t>(c)] = 1;
    e.event_color = comps.empty() ? st.num_boundary_colors : st.comp_seq_color[static_cast<size_t>(comps.front())];
    e.final = comps.empty();
    for (int32_t c : comps)
        if (st.comp_seq_color[static_cast<size_t>(c)] != e.event_color)
            throw std::runtime_error("owner_schedule: mixed-colour batch");
    return e;
}

/// Highest colour index among the pair's boundary endpoints whose components
/// were already consumed (not in the event, colour below the event's; every
/// boundary endpoint at the final event).  Sources with colour <= watermark
/// were applied at that earlier event; the endpoints' own components' sources
/// came in via R1.
template<typename CoordType, typename DataType>
int32_t pair_watermark(const OwnerScheduleState<CoordType, DataType>& st, const OwnerEvent& e, int64_t x, int64_t y) {
    int32_t wm = -1;
    for (int64_t z : {x, y}) {
        const int32_t c = st.comp_of(z);
        if (c < 0 || e.contains(c)) continue;
        const int32_t cc = st.comp_seq_color[static_cast<size_t>(c)];
        if (cc < e.event_color) wm = std::max(wm, cc);
    }
    return wm;
}

/// R2 / R3 filter.
template<typename CoordType, typename DataType>
DeferredPairFilter make_consumption_filter(const OwnerScheduleState<CoordType, DataType>& st, const OwnerEvent& e) {
    return [&st, &e](int64_t s, int64_t x, int64_t y) -> bool {
        if (!st.interest(x, y)) return false;
        const int32_t cx = st.comp_of(x), cy = st.comp_of(y), cs = st.comp_of(s);
        if (!e.final && !e.contains(cx) && !e.contains(cy)) return false;
        if (cs == cx || cs == cy) return false;
        return st.comp_seq_color[static_cast<size_t>(cs)] > pair_watermark(st, e, x, y);
    };
}

/// R1 filter: the batch's sources update exactly the pairs of interest with
/// an endpoint in the batch.
template<typename CoordType, typename DataType>
DeferredPairFilter make_batch_filter(const OwnerScheduleState<CoordType, DataType>& st, const OwnerEvent& e) {
    return [&st, &e](int64_t /*s*/, int64_t x, int64_t y) -> bool {
        if (!st.interest(x, y)) return false;
        return e.contains(st.comp_of(x)) || e.contains(st.comp_of(y));
    };
}

// ===========================================================================
// Passes (candidates -> owner -> mirror) for one source group
// ===========================================================================

template<typename CoordType, typename DataType, typename KernelType>
void owner_schedule_run_passes(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                               KernelType* kernel, const std::vector<int64_t>& sources,
                               const DeferredPairFilter* filter) {
    using clock = std::chrono::high_resolution_clock;
    if (sources.empty()) return;
    const auto t0 = clock::now();
    std::unordered_set<int64_t> wave_box_set(sources.begin(), sources.end());
    const int max_threads = std::max(1, omp_get_max_threads());

    std::vector<std::vector<int64_t>> thread_candidates(static_cast<size_t>(max_threads));
    {
        std::exception_ptr ex;
        std::mutex ex_mutex;
        std::atomic<bool> failed{false};
        #pragma omp parallel default(shared) if (sources.size() > 1)
        {
            auto& local = thread_candidates[static_cast<size_t>(omp_get_thread_num())];
            #pragma omp for schedule(static)
            for (int64_t i = 0; i < static_cast<int64_t>(sources.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    BoxData<CoordType, DataType>* box = level.find_local_box(sources[static_cast<size_t>(i)]);
                    if (box == nullptr) box = level.find_ghost_box(sources[static_cast<size_t>(i)]);
                    if (box == nullptr) continue;
                    collect_owner_deferred_xnn_candidates_for_source_box(
                        box, level, wave_box_set, local,
                        /*include_ghosts=*/true);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
        }
        if (ex) std::rethrow_exception(ex);
    }
    std::vector<int64_t> candidates;
    for (auto& v : thread_candidates) candidates.insert(candidates.end(), v.begin(), v.end());
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) return;
    FMM_PHASE_LAP_BEGIN(pass_lap);

    std::vector<std::vector<DeferredXnnTargetKey>> mirror_targets(candidates.size());
    {
        std::exception_ptr ex;
        std::mutex ex_mutex;
        std::atomic<bool> failed{false};
        const int split = split_threads_for(static_cast<int64_t>(candidates.size()), omp_get_max_threads());
        #pragma omp parallel default(shared)
        {
            DeferredXnnOwnerScratch<DataType> scratch;
            scratch.split_threads = split;
            #pragma omp for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(candidates.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    BoxData<CoordType, DataType>* box = level.find_local_box(candidates[static_cast<size_t>(i)]);
                    if (box == nullptr) box = level.find_ghost_box(candidates[static_cast<size_t>(i)]);
                    if (box == nullptr) continue;
                    apply_owner_deferred_xnn_updates_for_candidate_box(
                        box->morton_index, level, kernel, wave_box_set,
                        scratch, mirror_targets[static_cast<size_t>(i)],
                        static_cast<PendingFactorUpdates<DataType>*>(nullptr),
                        /*include_ghosts=*/true, filter);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
        }
        if (ex) std::rethrow_exception(ex);
    }
    FMM_PHASE_LAP(pass_lap, WALL_OWNER);

    {
        std::exception_ptr ex;
        std::mutex ex_mutex;
        std::atomic<bool> failed{false};
        #pragma omp parallel default(shared) if (candidates.size() > 1)
        {
            #pragma omp for schedule(static)
            for (int64_t i = 0; i < static_cast<int64_t>(candidates.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    BoxData<CoordType, DataType>* box = level.find_local_box(candidates[static_cast<size_t>(i)]);
                    if (box == nullptr) box = level.find_ghost_box(candidates[static_cast<size_t>(i)]);
                    if (box == nullptr) continue;
                    apply_symmetric_owner_deferred_xnn_updates_for_candidate_box(
                        box, mirror_targets[static_cast<size_t>(i)], level,
                        /*include_ghosts=*/true);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
        }
        if (ex) std::rethrow_exception(ex);
    }
    FMM_PHASE_LAP(pass_lap, WALL_MIRROR);
    st.t_passes_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

/// Replay `sources` under `filter` in ONE pass (design ledger item 5).  The
/// sources span several elimination waves; the canonical (wave, neighbour)
/// order is enforced per candidate inside
/// apply_owner_deferred_xnn_updates_for_candidate_box, so the result is
/// bitwise the one the earlier wave-by-wave replay produced (up to 24 small
/// passes per batch at L4, each with three parallel regions over a handful
/// of candidates — 7.5 s of R2 replay at L4 Helmholtz).  The mirror pass runs
/// once after all owner-side accumulation: the owner pass never reads the
/// mirror copies, so the order of mirror writes does not matter.
template<typename CoordType, typename DataType, typename KernelType>
void owner_schedule_replay(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                           KernelType* kernel, std::vector<int64_t> sources, const DeferredPairFilter& filter) {
    if (sources.empty()) return;
    std::sort(sources.begin(), sources.end(), [&](int64_t a, int64_t b) {
        const int32_t sa = st.seq_of(a), sb = st.seq_of(b);
        return sa != sb ? sa < sb : a < b;
    });
    owner_schedule_run_passes(st, level, kernel, sources, &filter);
}

/// Eliminated sources (present here) adjacent to any box of the batch's
/// components, outside the batch.
template<typename CoordType, typename DataType>
std::vector<int64_t> owner_schedule_pending_sources(const OwnerScheduleState<CoordType, DataType>& st,
                                                    TreeLevel<CoordType, DataType>& level, const OwnerEvent& e,
                                                    const std::vector<int32_t>& comps) {
    std::set<int64_t> out;
    for (int32_t id : comps) {
        const dataflow::Component& k = st.graph.comps[static_cast<size_t>(id)];
        for (int64_t m : k.boxes) {
            BoxData<CoordType, DataType>* box = level.find_local_box(m);
            if (box == nullptr) box = level.find_ghost_box(m);
            if (box == nullptr) continue;
            for (int64_t n : box->one_hop) {
                if (level.eliminated_boxes.count(n) == 0) continue;
                const int32_t cn = st.comp_of(n);
                if (cn < 0 || e.contains(cn)) continue;
                out.insert(n);
            }
        }
    }
    return std::vector<int64_t>(out.begin(), out.end());
}

// ===========================================================================
// Transport
// ===========================================================================

template<typename CoordType, typename DataType>
void owner_schedule_progress_outbox(OwnerScheduleState<CoordType, DataType>& st, bool wait_all) {
    for (auto it = st.outbox.begin(); it != st.outbox.end();) {
        int done = 0;
        if (wait_all) {
            owner_mpi_check(MPI_Waitall(static_cast<int>(it->reqs.size()), it->reqs.data(), MPI_STATUSES_IGNORE),
                            "MPI_Waitall(outbox)");
            done = 1;
        } else {
            owner_mpi_check(MPI_Testall(static_cast<int>(it->reqs.size()), it->reqs.data(), &done, MPI_STATUSES_IGNORE),
                            "MPI_Testall(outbox)");
        }
        if (done) {
            OwnerBufferPool::instance().release(std::move(it->buf));   // back to the pool, capacity kept
            it = st.outbox.erase(it);
        } else {
            ++it;
        }
    }
}

/// One piece of a payload on the wire: either a range of the packed buffer
/// (`ext == nullptr`, at `off`) or bytes read in place from a box (`ext`).
struct OwnerWireSeg {
    size_t off = 0;
    const char* ext = nullptr;
    size_t len = 0;
};

/// Per-component, per-destination payload under construction.  The receiver
/// sees the concatenation of the segments: [morton][kind][size][section]* with
/// every X_NR matrix in place inside its section.
struct OwnerPayloads {
    std::map<int32_t, std::map<uint32_t, std::vector<char>>> bufs;      ///< packed bytes
    std::map<int32_t, std::map<uint32_t, std::vector<OwnerWireSeg>>> segs;   ///< wire order
    std::map<int32_t, std::map<uint32_t, int64_t>> counts;               ///< boxes per payload
};

template<typename CoordType, typename DataType>
int owner_schedule_post_headers(OwnerScheduleState<CoordType, DataType>& st);

/// Snapshot the just-eliminated `boxes` (pre-pass) into their recipients'
/// payloads: per (box, recipient) a section tailored to that recipient (X_NR
/// row diet, home-rank blocks), serialized directly into the destination
/// buffer at a precomputed offset.
template<typename CoordType, typename DataType>
void owner_schedule_snapshot(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                             const std::vector<int64_t>& boxes, OwnerPayloads& out,
                             const std::unordered_map<int64_t, std::vector<int8_t>>& pre_flags) {
    using clock = std::chrono::high_resolution_clock;
    const auto t0 = clock::now();
    // One plan per (box, recipient): what to write and where.  Sizes are
    // computed first (parallel over boxes), each destination buffer is grown
    // once (pooled: its capacity survives the batch), then every plan is
    // serialized straight into its buffer at its offset (parallel over plans)
    // — no per-box blobs and no second copy.
    struct Plan {
        int64_t m = 0;
        int32_t comp = -1;
        uint32_t dest = 0;
        PayloadKind kind = PayloadKind::COMPACT;
        size_t size = 0;                      ///< section bytes on the wire, without the frame header
        size_t xnr_prefix = 0;                ///< packed bytes of the section before the X_NR matrix
        size_t xnr_ext = 0;                   ///< X_NR bytes sent in place from the box (0: packed)
        const char* xnr_ptr = nullptr;
        std::vector<int8_t> flags, rows_shipped;
        std::vector<int64_t> ship_blocks;     ///< near blocks shipped (sorted)
        std::vector<char>* dst = nullptr;
        size_t off = 0;                       ///< offset of the frame header in the packed buffer
        size_t packed() const { return size - xnr_ext; }
    };
    constexpr size_t HEADER = sizeof(int64_t) + sizeof(int8_t) + sizeof(size_t);
    std::vector<std::vector<Plan>> plans(boxes.size());
    std::atomic<bool> failed{false};
    std::atomic<size_t> xnr_ship{0}, xnr_kern{0}, xnr_part{0};   // X_NR bytes by row class
    #pragma omp parallel for schedule(dynamic)
    for (int64_t j = 0; j < static_cast<int64_t>(boxes.size()); ++j) {
        if (failed.load(std::memory_order_relaxed)) continue;
        const int64_t m = boxes[static_cast<size_t>(j)];
        BoxData<CoordType, DataType>* box = level.find_local_box(m);
        if (box == nullptr) box = level.find_ghost_box(m);
        if (box == nullptr) { failed.store(true); continue; }
        auto pf = pre_flags.find(m);
        if (pf == pre_flags.end()) { failed.store(true); continue; }
        const std::vector<int8_t>& flags = pf->second;   // BLOCK_KERNEL / BLOCK_PARTNER per one_hop slot
        const auto& counts = box->deferred_xnn_neighbor_point_counts;
        if (flags.size() != box->one_hop.size() || counts.size() != box->one_hop.size()) { failed.store(true); continue; }
        const int64_t r_cols = box->X_NR.is_allocated() ? box->X_NR.cols : 0;
        const int32_t comp = st.comp_of(m);
        for (const auto& [dest, kind] : st.recipients_of_box(m)) {
            Plan p;
            p.m = m;
            p.comp = comp;
            p.dest = dest;
            p.kind = kind;
            if (kind == PayloadKind::SKELETON) {
                p.size = skeleton_section_size(*box);
                plans[static_cast<size_t>(j)].push_back(std::move(p));
                continue;
            }
            // Full shipping (see the payload note at the top): every X_NR row
            // block, and the near block of every BLOCK_PARTNER slot the
            // recipient will hold (all slots for the home rank, slots of
            // interest otherwise).  BLOCK_KERNEL slots are rebuilt there from
            // the kernel.  Nothing at install reads a local pair copy.
            p.flags = flags;
            p.rows_shipped.assign(box->one_hop.size(), 1);
            for (size_t i = 0; i < box->one_hop.size(); ++i) {
                const int64_t g = box->one_hop[i];
                if (flags[i] == BLOCK_PARTNER && (kind == PayloadKind::FULL || st.interest_of(dest, g))) {
                    p.ship_blocks.push_back(g);
                    p.flags[i] = BLOCK_SHIPPED;
                }
                xnr_ship += static_cast<size_t>(counts[i] * r_cols) * sizeof(DataType);
            }
            std::sort(p.ship_blocks.begin(), p.ship_blocks.end());
            const auto& sb = p.ship_blocks;
            const std::function<bool(int64_t)> keep = [&sb](int64_t nb) { return std::binary_search(sb.begin(), sb.end(), nb); };
            p.size = compact_state_size(*box, keep, p.flags, p.rows_shipped);
            p.xnr_prefix = compact_state_xnr_prefix(*box, p.flags, p.rows_shipped);
            if (xnr_contiguous(*box)) {   // send the matrix in place, no copy
                p.xnr_ext = compact_state_xnr_bytes(*box, p.rows_shipped);
                p.xnr_ptr = reinterpret_cast<const char*>(box->X_NR.data.data());
            }
            plans[static_cast<size_t>(j)].push_back(std::move(p));
        }
    }
    if (failed.load()) throw std::runtime_error("owner_schedule_snapshot: box or its pre-elimination flags missing");
    st.xnr_bytes_shipped += xnr_ship.load();
    st.xnr_bytes_regen_kernel += xnr_kern.load();
    st.xnr_bytes_regen_partner += xnr_part.load();

    // Place: one growth per destination buffer.  A new destination buffer
    // comes from the pool with the largest capacity seen so far reserved, so
    // that later waves of the batch append without reallocating.  The wire
    // order is recorded as segments: packed bytes up to the X_NR matrix, the
    // matrix from the box, packed bytes after it (adjacent packed segments
    // merged).  Offsets, not pointers: the buffer may still grow.
    std::vector<Plan*> all;
    std::map<std::vector<char>*, size_t> grow;
    auto add_packed = [](std::vector<OwnerWireSeg>& segs, size_t off, size_t len) {
        if (len == 0) return;
        if (!segs.empty() && segs.back().ext == nullptr && segs.back().off + segs.back().len == off) segs.back().len += len;
        else segs.push_back({off, nullptr, len});
    };
    for (auto& per_box : plans)
        for (Plan& p : per_box) {
            std::vector<char>& dst = out.bufs[p.comp][p.dest];
            if (dst.capacity() == 0) dst = OwnerBufferPool::instance().acquire(0, OwnerBufferPool::instance().max_seen());
            p.dst = &dst;
            size_t& g = grow[&dst];
            p.off = dst.size() + g;
            g += HEADER + p.packed();
            out.counts[p.comp][p.dest] += 1;
            std::vector<OwnerWireSeg>& segs = out.segs[p.comp][p.dest];
            if (p.xnr_ext) {
                add_packed(segs, p.off, HEADER + p.xnr_prefix);
                segs.push_back({0, p.xnr_ptr, p.xnr_ext});
                add_packed(segs, p.off + HEADER + p.xnr_prefix, p.packed() - p.xnr_prefix);
            } else {
                add_packed(segs, p.off, HEADER + p.packed());
            }
            all.push_back(&p);
        }
    for (auto& kv : grow) kv.first->resize(kv.first->size() + kv.second);

    // Serialize every plan in place: [morton][kind][size][section]
    #pragma omp parallel for schedule(dynamic)
    for (int64_t i = 0; i < static_cast<int64_t>(all.size()); ++i) {
        if (failed.load(std::memory_order_relaxed)) continue;
        const Plan& p = *all[static_cast<size_t>(i)];
        BoxData<CoordType, DataType>* box = level.find_local_box(p.m);
        if (box == nullptr) box = level.find_ghost_box(p.m);
        if (box == nullptr) { failed.store(true); continue; }
        char* q = p.dst->data() + p.off;
        std::memcpy(q, &p.m, sizeof(int64_t));
        q += sizeof(int64_t);
        const int8_t k8 = static_cast<int8_t>(p.kind);
        std::memcpy(q, &k8, sizeof(int8_t));
        q += sizeof(int8_t);
        std::memcpy(q, &p.size, sizeof(size_t));
        q += sizeof(size_t);
        char* end = nullptr;
        if (p.kind == PayloadKind::SKELETON) {
            end = serialize_skeleton_section(*box, q);
        } else {
            const auto& sb = p.ship_blocks;
            const std::function<bool(int64_t)> keep = [&sb](int64_t nb) { return std::binary_search(sb.begin(), sb.end(), nb); };
            end = serialize_compact_state(*box, keep, q, p.flags, p.rows_shipped, p.xnr_ext != 0);
        }
        if (end != q + p.packed()) failed.store(true);
    }
    if (failed.load()) throw std::runtime_error("owner_schedule_snapshot: pack size mismatch");
    st.t_pack_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

/// Ship every (component, destination) payload of the batch.
template<typename CoordType, typename DataType>
void owner_schedule_push(OwnerScheduleState<CoordType, DataType>& st, OwnerPayloads& payloads) {
    for (auto& per_comp : payloads.bufs) {
        const int32_t id = per_comp.first;
        for (auto& per_dest : per_comp.second) {
            // A batch can contain thousands of large wire segments.  Keep
            // older sends moving and post peers' matching receives before
            // adding another destination to the outbox.
            owner_schedule_progress_outbox(st, false);
            owner_schedule_post_headers(st);
            const uint32_t dest_pid = per_dest.first;
            auto it = st.pid_to_rank.find(static_cast<int>(dest_pid));
            if (it == st.pid_to_rank.end()) throw std::runtime_error("owner_schedule_push: no rank for pid " + std::to_string(dest_pid));
            const int dest = it->second;
            std::vector<char>& buf = per_dest.second;
            const std::vector<OwnerWireSeg>& segs = payloads.segs[id][dest_pid];
            // wire parts: every segment, split at PART_BYTES; X_NR segments
            // are sent straight from the box's storage (kept unchanged and
            // alive until the outbox is drained at the end of the level)
            size_t total = 0;
            int64_t nparts = 0;
            for (const auto& sg : segs) {
                total += sg.len;
                nparts += static_cast<int64_t>((sg.len + st.PART_BYTES - 1) / st.PART_BYTES);
            }
            st.outbox.emplace_back();
            auto& msg = st.outbox.back();
            msg.buf = std::move(buf);
            // header: [id][boxes][total bytes][nparts][len of every part] — the
            // lengths let the receiver post all its receives at once
            msg.header = {static_cast<int64_t>(id), payloads.counts[id][dest_pid], static_cast<int64_t>(total), nparts};
            for (const auto& sg : segs)
                for (size_t off = 0; off < sg.len; off += st.PART_BYTES)
                    msg.header.push_back(static_cast<int64_t>(std::min(st.PART_BYTES, sg.len - off)));
            MPI_Request r;
            owner_schedule_post_headers(st);
            owner_mpi_check(MPI_Isend(msg.header.data(), static_cast<int>(msg.header.size()), MPI_INT64_T, dest, st.TAG_HDR,
                                      st.comm, &r),
                            "MPI_Isend(header)");
            msg.reqs.push_back(r);
            for (const auto& sg : segs) {
                const char* base = sg.ext ? sg.ext : msg.buf.data() + sg.off;
                for (size_t off = 0; off < sg.len; off += st.PART_BYTES) {
                    const size_t len = std::min(st.PART_BYTES, sg.len - off);
                    // MPI_Isend may internally wait for transport resources.
                    // Polling headers between postings prevents a set of
                    // owners from filling those resources while none of them
                    // reaches the receive-posting path.
                    owner_schedule_post_headers(st);
                    owner_mpi_check(MPI_Isend(base + off, static_cast<int>(len), MPI_CHAR, dest, st.TAG_PART,
                                              st.comm, &r),
                                    "MPI_Isend(part)");
                    msg.reqs.push_back(r);
                }
            }
            st.bytes_sent += total;
            st.parts_sent += nparts;
        }
    }
    owner_schedule_progress_outbox(st, false);
    owner_schedule_post_headers(st);
    payloads.bufs.clear();
    payloads.segs.clear();
    payloads.counts.clear();
}

/// A payload whose header has arrived and whose parts are landing in the
/// background: the receives are posted at once (the header lists every
/// part's length), the MPI progress thread moves the bytes while this rank
/// computes, and owner_schedule_complete_receives hands the buffer to the
/// inbox once every part is in.
template<typename CoordType, typename DataType>
void owner_schedule_receive_one(OwnerScheduleState<CoordType, DataType>& st, const MPI_Status& hdr_status) {
    const int source = hdr_status.MPI_SOURCE;
    int hdr_len = 0;
    owner_mpi_check(MPI_Get_count(&hdr_status, MPI_INT64_T, &hdr_len), "MPI_Get_count(header)");
    if (hdr_len < 4) throw std::runtime_error("owner_schedule_receive_one: short header");
    std::vector<int64_t> header(static_cast<size_t>(hdr_len));
    owner_mpi_check(MPI_Recv(header.data(), hdr_len, MPI_INT64_T, source, st.TAG_HDR, st.comm, MPI_STATUS_IGNORE),
                    "MPI_Recv(header)");
    const int32_t id = static_cast<int32_t>(header[0]);
    const size_t total = static_cast<size_t>(header[2]);
    const int64_t nparts = header[3];
    if (hdr_len != 4 + nparts) throw std::runtime_error("owner_schedule_receive_one: header does not list every part");
    typename OwnerScheduleState<CoordType, DataType>::PendingRecv pr;
    pr.id = id;
    pr.total = total;
    pr.buf = OwnerBufferPool::instance().acquire(total);   // pooled: no fresh pages
    bool found = false;
    for (const auto& kv : st.pid_to_rank) if (kv.second == source) { pr.from = static_cast<uint32_t>(kv.first); found = true; break; }
    if (!found) throw std::runtime_error("owner_schedule_receive_one: unknown sender rank");
    // parts match the sends in posting order (same source, same tag), so the
    // receives are posted in header order at their offsets
    size_t off = 0;
    pr.reqs.resize(static_cast<size_t>(nparts));
    for (int64_t p = 0; p < nparts; ++p) {
        const int64_t len = header[static_cast<size_t>(4 + p)];
        if (len < 0 || off + static_cast<size_t>(len) > total)
            throw std::runtime_error("owner_schedule_receive_one: payload part overruns the announced size");
        owner_mpi_check(MPI_Irecv(pr.buf.data() + off, static_cast<int>(len), MPI_CHAR, source, st.TAG_PART, st.comm,
                                  &pr.reqs[static_cast<size_t>(p)]),
                        "MPI_Irecv(part)");
        off += static_cast<size_t>(len);
    }
    if (off != total) throw std::runtime_error("owner_schedule_receive_one: payload parts short of the announced size");
    st.parts_recv += nparts;
    st.pending_recv.push_back(std::move(pr));
}

/// Post matching receives for every owner payload header currently available.
/// This is separate from owner_schedule_poll so the send-posting loop can make
/// receive-side progress without touching the outbox entry it is still building.
template<typename CoordType, typename DataType>
int owner_schedule_post_headers(OwnerScheduleState<CoordType, DataType>& st) {
    int posted = 0;
    while (true) {
        int flag = 0;
        MPI_Status status;
        owner_mpi_check(
            MPI_Iprobe(MPI_ANY_SOURCE, st.TAG_HDR, st.comm, &flag, &status),
            "MPI_Iprobe(header)");
        if (!flag) break;
        owner_schedule_receive_one(st, status);
        ++posted;
    }
    return posted;
}

/// Hand every fully landed payload to the inbox and the runtime.  Returns
/// how many completed.
template<typename CoordType, typename DataType>
int owner_schedule_complete_receives(OwnerScheduleState<CoordType, DataType>& st) {
    int done_count = 0;
    for (auto it = st.pending_recv.begin(); it != st.pending_recv.end();) {
        int done = 0;
        owner_mpi_check(MPI_Testall(static_cast<int>(it->reqs.size()), it->reqs.data(), &done, MPI_STATUSES_IGNORE),
                        "MPI_Testall(parts)");
        if (!done) { ++it; continue; }
        st.bytes_recv += it->total;
        const int32_t id = it->id;
        if (st.rt.is_aux(id)) {
            if (st.aux_inbox.count(id)) throw std::runtime_error("owner_schedule: duplicate auxiliary payload for component " + std::to_string(id));
            st.aux_inbox[id] = std::move(it->buf);
            st.rt.on_aux_arrived(id, it->from, omp_get_wtime());
        } else {
            if (st.inbox.count(id)) throw std::runtime_error("owner_schedule: duplicate payload for component " + std::to_string(id));
            st.inbox[id] = std::move(it->buf);
            st.rt.on_payload_arrived(id, it->from);
        }
        it = st.pending_recv.erase(it);
        ++done_count;
    }
    return done_count;
}

/// Service the wire: complete sends, post receives for every announced
/// payload, hand landed payloads over.  `block`: do not return before at
/// least one payload has landed (or a header has been posted when nothing
/// was pending) — the callers re-check their condition and call again.
template<typename CoordType, typename DataType>
void owner_schedule_poll(OwnerScheduleState<CoordType, DataType>& st, bool block) {
    using clock = std::chrono::high_resolution_clock;
    owner_schedule_progress_outbox(st, false);
    owner_schedule_post_headers(st);
    if (owner_schedule_complete_receives(st) > 0 || !block) return;
    const auto t0 = clock::now();
    if (st.pending_recv.empty()) {
        // nothing in flight: sleep in the probe until a header shows up
        MPI_Status status;
        owner_mpi_check(MPI_Probe(MPI_ANY_SOURCE, st.TAG_HDR, st.comm, &status), "MPI_Probe(header)");
        st.t_wait_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        owner_schedule_receive_one(st, status);
        owner_schedule_post_headers(st);
        owner_schedule_complete_receives(st);
        return;
    }
    // parts in flight: spin on their completion (the main thread has nothing
    // else to do), taking new headers as they come
    while (true) {
        owner_schedule_post_headers(st);
        if (owner_schedule_complete_receives(st) > 0) break;
        if (st.pending_recv.empty()) break;
        int idx = 0, flag = 0;
        owner_mpi_check(MPI_Testany(static_cast<int>(st.pending_recv.front().reqs.size()), st.pending_recv.front().reqs.data(),
                                    &idx, &flag, MPI_STATUS_IGNORE),
                        "MPI_Testany(parts)");   // drives progress; completion is read by Testall above
    }
    st.t_recv_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

// ===========================================================================
// Batch event
// ===========================================================================

template<typename CoordType, typename DataType>
void owner_schedule_mark_eliminated(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                                    const std::vector<int64_t>& boxes) {
    for (int64_t m : boxes) {
        level.eliminated_boxes.insert(m);
        level.elimination_wave[m] = st.seq_of(m);
    }
}

/// Install every arrived auxiliary payload of colour < `below_color`
/// (all of them for below_color >= num_boundary_colors): skeleton sections
/// only, boxes marked eliminated.
template<typename CoordType, typename DataType>
void owner_schedule_install_aux(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                                int32_t below_color) {
    for (auto it = st.aux_inbox.begin(); it != st.aux_inbox.end();) {
        const int32_t id = it->first;
        if (st.comp_seq_color[static_cast<size_t>(id)] >= below_color) { ++it; continue; }
        const char* p = it->second.data();
        const char* end = p + it->second.size();
        std::vector<int64_t> boxes;
        while (p < end) {
            int64_t m = 0;
            int8_t k8 = 0;
            size_t sz = 0;
            std::memcpy(&m, p, sizeof(int64_t));
            p += sizeof(int64_t);
            std::memcpy(&k8, p, sizeof(int8_t));
            p += sizeof(int8_t);
            std::memcpy(&sz, p, sizeof(size_t));
            p += sizeof(size_t);
            if (static_cast<PayloadKind>(k8) != PayloadKind::SKELETON)
                throw std::runtime_error("owner_schedule_install_aux: non-skeleton entry in auxiliary payload");
            BoxData<CoordType, DataType>* box = level.find_local_box(m);
            if (box == nullptr) box = level.find_ghost_box(m);
            if (box == nullptr)
                throw std::runtime_error("owner_schedule_install_aux: rank " + std::to_string(st.my_rank) +
                                         " does not hold box " + std::to_string(m));
            if (deserialize_skeleton_section(*box, p) != p + sz)
                throw std::runtime_error("owner_schedule_install_aux: unpack size mismatch");
            p += sz;
            boxes.push_back(m);
        }
        owner_schedule_mark_eliminated(st, level, boxes);
        st.aux_installed.insert(id);
        OwnerBufferPool::instance().release(std::move(it->second));
        it = st.aux_inbox.erase(it);
    }
}

/// A batch in progress: `owned` components are eliminated here, `installed`
/// components' payloads are installed here; all of one colour.  Resumable
/// wave by wave so that the scheduler can poll for arrivals between waves
/// and let a lower-colour component preempt the remaining waves (design §9).
template<typename CoordType, typename DataType>
struct OwnerBatch {
    std::vector<int32_t> owned, installed;
    OwnerEvent ev;
    int32_t color = 0;
    struct Entry { const char* p; size_t size; PayloadKind kind; };
    std::unordered_map<int64_t, Entry> entries;   ///< installed components' payload entries
    DeferredPairFilter f1;                        ///< R1 filter (captures ev by reference)
    OwnerPayloads payloads;
    int next_wave = 0;
    int waves_run = 0;
    bool installs_done = false;                          ///< all copies are installed in the batch's first round
    std::vector<std::vector<int64_t>> installed_by_wave; ///< installed sources per wave, for the per-wave passes
};

template<typename CoordType, typename DataType, typename KernelType>
std::unique_ptr<OwnerBatch<CoordType, DataType>> owner_schedule_batch_begin(
    OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level, KernelType* kernel,
    const std::vector<int32_t>& owned, const std::vector<int32_t>& installed) {
    using clock = std::chrono::high_resolution_clock;
    auto b = std::make_unique<OwnerBatch<CoordType, DataType>>();
    b->owned = owned;
    b->installed = installed;
    std::vector<int32_t> comps = owned;
    comps.insert(comps.end(), installed.begin(), installed.end());
    if (comps.empty()) throw std::runtime_error("owner_schedule_batch_begin: empty batch");
    ++st.batches;
    b->ev = make_event(st, comps);
    b->color = b->ev.event_color;
    owner_schedule_install_aux(st, level, b->ev.event_color);   // lower-colour far endpoints

    for (int32_t id : installed) {
        auto it = st.inbox.find(id);
        if (it == st.inbox.end()) throw std::runtime_error("owner_schedule_batch_begin: payload missing for component " + std::to_string(id));
        const char* p = it->second.data();
        const char* end = p + it->second.size();
        while (p < end) {
            int64_t m = 0;
            int8_t k8 = 0;
            size_t sz = 0;
            std::memcpy(&m, p, sizeof(int64_t));
            p += sizeof(int64_t);
            std::memcpy(&k8, p, sizeof(int8_t));
            p += sizeof(int8_t);
            std::memcpy(&sz, p, sizeof(size_t));
            p += sizeof(size_t);
            b->entries[m] = {p, sz, static_cast<PayloadKind>(k8)};
            p += sz;
        }
    }

    // R2: pending lower sources onto this rank's copies of pairs with an endpoint in the batch
    {
        const auto t0 = clock::now();
        owner_schedule_replay(st, level, kernel, owner_schedule_pending_sources(st, level, b->ev, comps),
                              make_consumption_filter(st, b->ev));
        st.t_replay_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    }
    b->f1 = make_batch_filter(st, b->ev);
    return b;
}

/// Run the batch's next wave.  Returns false when no waves remain.
template<typename CoordType, typename DataType, typename KernelType>
bool owner_schedule_batch_wave(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                               KernelType* kernel, OwnerBatch<CoordType, DataType>& b,
                               const std::vector<CoordType>& unit_proxy_points, int num_proxy, CoordType proxy_radius,
                               double tolerance, bool is_symmetric, bool is_hermitian,
                               FactorizationMethod factorization_method, int use_sketch) {
    using clock = std::chrono::high_resolution_clock;
    using Entry = typename OwnerBatch<CoordType, DataType>::Entry;
    while (b.next_wave < st.num_waves) {
        const int w = b.next_wave++;
        std::vector<int64_t> wave_owned, wave_installed, wave_skeleton;
        for (int32_t id : b.owned)
            for (int64_t m : st.graph.comps[static_cast<size_t>(id)].boxes)
                if (st.wave_of(m) == w) wave_owned.push_back(m);
        // Copies: ALL of them are installed (payload written, temp2 and kernel
        // blocks rebuilt, marked eliminated) in the batch's first round — an
        // installed box reads nothing local (payload note at the top), so this
        // is one parallel loop over every copy instead of 8 thin ones.  Their
        // PASSES still run wave by wave (`installed_by_wave`): a copy's own
        // "replace" of its pairs and a later source's accumulation onto the
        // same pair must not meet in one pass (the deferred mirror write at
        // the end of a pass would clobber the accumulation).  Older sources
        // never touch a copy installed after them: see the newer-candidate
        // rules in apply_owner_deferred_xnn_updates_for_candidate_box.
        if (!b.installs_done) {
            b.installs_done = true;
            b.installed_by_wave.assign(static_cast<size_t>(st.num_waves), {});
            for (int32_t id : b.installed)
                for (int64_t m : st.graph.comps[static_cast<size_t>(id)].boxes) {
                    auto e = b.entries.find(m);
                    if (e == b.entries.end()) continue;
                    if (e->second.kind == PayloadKind::SKELETON) {
                        wave_skeleton.push_back(m);
                    } else {
                        wave_installed.push_back(m);
                        b.installed_by_wave[static_cast<size_t>(st.wave_of(m))].push_back(m);
                    }
                }
            std::sort(wave_installed.begin(), wave_installed.end());
            std::sort(wave_skeleton.begin(), wave_skeleton.end());
            for (auto& v : b.installed_by_wave) std::sort(v.begin(), v.end());
        }
        const std::vector<int64_t> pass_installed =
            b.installed_by_wave.empty() ? std::vector<int64_t>{} : b.installed_by_wave[static_cast<size_t>(w)];
        if (wave_owned.empty() && wave_installed.empty() && wave_skeleton.empty() && pass_installed.empty()) continue;
        ++b.waves_run;

        // owned: per-box region
        if (!wave_owned.empty()) {
            const auto t0 = clock::now();
            FMM_PHASE_NOTE_WAVE(wave_owned.size());
            FMM_PHASE_LAP_BEGIN(wave_lap);
            std::unordered_set<int64_t> wave_box_set(wave_owned.begin(), wave_owned.end());
            std::exception_ptr ex;
            std::mutex ex_mutex;
            std::atomic<bool> failed{false};
            const int wave_split = split_threads_for(static_cast<int64_t>(wave_owned.size()), omp_get_max_threads());
            std::vector<std::vector<int8_t>> pre(wave_owned.size());   // stored-block provenance before elimination
            level.sketch_eliminated_filter = &st.sketch_filter;
            #pragma omp parallel default(shared)
            {
                FactorizationThreadScratch<CoordType, DataType> scratch;
                scratch.split_threads = wave_split;
                #pragma omp for schedule(dynamic)
                for (int64_t i = 0; i < static_cast<int64_t>(wave_owned.size()); ++i) {
                    if (failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const int64_t m = wave_owned[static_cast<size_t>(i)];
                        owner_current_sketch_comp = st.comp_of(m);
                        FMM_PHASE_BOX_SCOPE();
                        FMM_PHASE_LAP_BEGIN(box_lap);
                        BoxData<CoordType, DataType>* box_ptr = level.find_local_box(m);
                        if (box_ptr == nullptr) box_ptr = level.find_ghost_box(m);
                        if (box_ptr == nullptr) throw std::runtime_error("owner_schedule: owned box missing at elimination");
                        auto& box = *box_ptr;
                        {
                            auto& f = pre[static_cast<size_t>(i)];
                            f.assign(box.one_hop.size(), BLOCK_KERNEL);
                            for (size_t q = 0; q < box.one_hop.size(); ++q)
                                if (box.near_field_interaction_map.count(box.one_hop[q])) f[q] = BLOCK_PARTNER;
                        }
                        if (use_sketch == 2) {
                            gather_id_target_streamed(
                                st.tree, &box, level, kernel,
                                scratch, box.on_boundary);
                        } else {
                            scratch.streamed_sketch_valid = false;
                            gather_id_workspace(
                                st.tree, &box, level, kernel, tolerance,
                                unit_proxy_points.data(), num_proxy,
                                proxy_radius, is_symmetric,
                                scratch.workspace, scratch.workspace_rows,
                                scratch.workspace_cols, 0,
                                box.on_boundary,
                                /*use_CA_boundary_semantics=*/true);
                        }
                        FMM_PHASE_LAP(box_lap, BOX_SKETCH);
                        compute_and_modify(
                            st.dimension, &box, level, kernel, scratch,
                            tolerance, use_sketch,
                            is_symmetric, is_hermitian,
                            static_cast<PendingFactorUpdates<DataType>*>(nullptr),
                            factorization_method, true);
                    } catch (...) {
                        if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                    }
                }
            }
            level.sketch_eliminated_filter = nullptr;
            owner_current_sketch_comp = -1;
            FMM_PHASE_LAP(wave_lap, WALL_BOX);
            FMM_PHASE_WAVE_BOX_DONE();
            if (ex) std::rethrow_exception(ex);
            st.boxes_eliminated += static_cast<int64_t>(wave_owned.size());
            st.t_elim_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
            std::unordered_map<int64_t, std::vector<int8_t>> pre_flags;
            for (size_t i = 0; i < wave_owned.size(); ++i) pre_flags[wave_owned[i]] = std::move(pre[i]);
            owner_schedule_snapshot(st, level, wave_owned, b.payloads, pre_flags);   // pre-pass
        }

        // installed: write the wave's snapshots (or skeleton sections) into this
        // rank's copies, then rebuild temp2 and the needed near blocks (diet)
        if (!wave_installed.empty() || !wave_skeleton.empty()) {
            const auto t0 = clock::now();
            for (int64_t m : wave_skeleton) {
                BoxData<CoordType, DataType>* box = level.find_local_box(m);
                if (box == nullptr) box = level.find_ghost_box(m);
                if (box == nullptr)
                    throw std::runtime_error("owner_schedule_batch_wave: rank " + std::to_string(st.my_rank) +
                                             " does not hold shipped box " + std::to_string(m));
                const Entry e = b.entries.at(m);
                if (deserialize_skeleton_section(*box, e.p) != e.p + e.size)
                    throw std::runtime_error("owner_schedule_batch_wave: unpack size mismatch (skeleton)");
            }
            std::exception_ptr ex;
            std::mutex ex_mutex;
            std::atomic<bool> failed{false};
            std::vector<double> th_temp2(static_cast<size_t>(std::max(1, omp_get_max_threads())), 0.0), th_blocks = th_temp2,
                th_xnr = th_temp2;
            // Split each box's temp2 solve and kernel blocks over K threads
            // (taskloops on this team) when the round has fewer boxes than
            // threads, as the box region does.
            const int install_split = split_threads_for(static_cast<int64_t>(wave_installed.size()), omp_get_max_threads());
            #pragma omp parallel for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(wave_installed.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    const int64_t m = wave_installed[static_cast<size_t>(i)];
                    BoxData<CoordType, DataType>* box = level.find_local_box(m);
                    if (box == nullptr) box = level.find_ghost_box(m);
                    if (box == nullptr)
                        throw std::runtime_error("owner_schedule_batch_wave: rank " + std::to_string(st.my_rank) +
                                                 " does not hold shipped box " + std::to_string(m));
                    const Entry e = b.entries.at(m);
                    std::vector<int8_t> flags, rows_shipped;
                    if (deserialize_compact_state(*box, e.p, flags, rows_shipped) != e.p + e.size)
                        throw std::runtime_error("owner_schedule_batch_wave: unpack size mismatch");
                    // X_NR arrives complete (regeneration returns at once); then
                    // temp2 and the kernel-slot near blocks toward this rank's boxes
                    owner_regenerate_x_nr(level, kernel, *box, flags, rows_shipped, st.dimension, install_split,
                                          &th_xnr[static_cast<size_t>(omp_get_thread_num())]);
                    owner_rebuild_shipped_source(st.interest_filter, level, kernel, factorization_method, *box, flags,
                                                 st.dimension, e.kind == PayloadKind::FULL,
                                                 &th_temp2[static_cast<size_t>(omp_get_thread_num())],
                                                 &th_blocks[static_cast<size_t>(omp_get_thread_num())], install_split);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
            if (ex) std::rethrow_exception(ex);
            for (size_t t = 0; t < th_temp2.size(); ++t) {
                st.t_rebuild_temp2_ms += th_temp2[t];
                st.t_rebuild_blocks_ms += th_blocks[t];
                st.t_rebuild_xnr_ms += th_xnr[t];
            }
            st.t_install_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        }

        // sources for this wave's passes: everything with temp2 here; skeleton-
        // only boxes are marked eliminated (their skeleton is what 2-hop readers
        // need) but are never sources on this rank
        // installed copies are marked as soon as they are installed (all in
        // the first round), owned boxes at their wave; this wave's passes take
        // the owned sources of the wave and the installed sources of the wave
        owner_schedule_mark_eliminated(st, level, wave_installed);
        owner_schedule_mark_eliminated(st, level, wave_skeleton);
        owner_schedule_mark_eliminated(st, level, wave_owned);
        std::vector<int64_t> wave_all = wave_owned;
        wave_all.insert(wave_all.end(), pass_installed.begin(), pass_installed.end());
        owner_schedule_run_passes(st, level, kernel, wave_all, &b.f1);
        return b.next_wave < st.num_waves;   // one wave per call: the caller polls in between
    }
    return false;
}

/// Ship the batch's payloads and tell the runtime its components are done.
template<typename CoordType, typename DataType>
void owner_schedule_batch_finish(OwnerScheduleState<CoordType, DataType>& st, OwnerBatch<CoordType, DataType>& b) {
    owner_schedule_push(st, b.payloads);
    for (int32_t id : b.installed) {
        auto it = st.inbox.find(id);
        if (it != st.inbox.end()) {
            OwnerBufferPool::instance().release(std::move(it->second));
            st.inbox.erase(it);
        }
    }
    const double now = omp_get_wtime();
    for (int32_t id : b.owned) { st.rt.on_owned_finished(id, now); ++st.comps_eliminated; }
    for (int32_t id : b.installed) { st.rt.on_installed(id, now); ++st.comps_installed; }
}

/// Final event: every pair of interest receives its pending sources in
/// canonical order; boundary sources are finalized (temp2 -> X_NR).
template<typename CoordType, typename DataType, typename KernelType>
void owner_schedule_final_event(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                                KernelType* kernel) {
    using clock = std::chrono::high_resolution_clock;
    const auto t0 = clock::now();
    const OwnerEvent ev = make_event(st, {});
    owner_schedule_install_aux(st, level, st.num_boundary_colors);   // everything still pending
    std::vector<int64_t> sources;
    for (int64_t m : level.eliminated_boxes) if (st.comp_of(m) >= 0) sources.push_back(m);
    owner_schedule_replay(st, level, kernel, sources, make_consumption_filter(st, ev));
    #pragma omp parallel for schedule(dynamic)
    for (int64_t i = 0; i < static_cast<int64_t>(sources.size()); ++i) {
        BoxData<CoordType, DataType>* box = level.find_local_box(sources[static_cast<size_t>(i)]);
        if (box == nullptr) box = level.find_ghost_box(sources[static_cast<size_t>(i)]);
        if (box == nullptr) continue;
        // finalize_deferred_xnn_source_box, but the X_NR buffer is handed to
        // st.xnr_hold instead of being freed: sends still read it in place
        // (released after the outbox drain in owner_schedule_run_boundary)
        if (!box->deferred_xnn_temp2.empty()) {
            if (!box->X_NR.is_allocated()) throw std::runtime_error("owner_schedule_final_event: source X_NR missing");
            std::vector<DataType> old;
            old.swap(box->X_NR.data);
            box->X_NR.data = std::move(box->deferred_xnn_temp2);
            #pragma omp critical(owner_xnr_hold)
            st.xnr_hold.push_back(std::move(old));
        }
        std::vector<DataType>().swap(box->deferred_xnn_temp2);
    }
    st.t_final_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
}

// ===========================================================================
// The boundary phase
// ===========================================================================

template<typename CoordType, typename DataType, typename KernelType>
void owner_schedule_run_boundary(OwnerScheduleState<CoordType, DataType>& st, TreeLevel<CoordType, DataType>& level,
                                 KernelType* kernel, const std::vector<CoordType>& unit_proxy_points, int num_proxy,
                                 CoordType proxy_radius, double tolerance, bool is_symmetric, bool is_hermitian,
                                 FactorizationMethod factorization_method, int use_sketch) {
    if (!st.active) return;
    auto& rt = st.rt;
    using Batch = OwnerBatch<CoordType, DataType>;
    auto wave = [&](Batch& b) {
        return owner_schedule_batch_wave(st, level, kernel, b, unit_proxy_points, num_proxy, proxy_radius, tolerance,
                                         is_symmetric, is_hermitian, factorization_method, use_sketch);
    };
    auto all_done = [&]() {
        if (!rt.all_owned_done()) return false;
        for (int32_t id : rt.shared_not_owned()) if (!rt.is_done(id)) return false;
        return rt.boundary_complete();   // every expected payload (incl. auxiliary) consumed
    };
    auto color_of = [&](int32_t id) { return st.comp_seq_color[static_cast<size_t>(id)]; };

    // Coarse levels run the SYNCHRONOUS colour schedule below: one batch per
    // colour with every owned component of the colour, then one batch with
    // every shared copy of the colour, waiting for completeness before each.
    // Arrival-driven batching would split the same work into many small box
    // regions (ledger 27: 40-46 rounds per level against colour's 16), and at
    // a coarse level a full wave is what keeps the threads busy (ledger 28b:
    // at 4.5 boxes/thread the rounds went 47 -> 32 for 0.4 s of waiting).
    // Fine levels run arrival-driven with the cumulative balancing, so that a
    // rank light in blue starts its (heavier) green early.  The choice is made
    // in owner_schedule_setup (st.sync_colors), before the owner assignment,
    // because the balancing objective depends on it.

    if (st.serial) {
        // Debug oracle: components in id order, one per batch, blocking waits,
        // no preemption.  Same numerics as the fast path.
        std::vector<int32_t> mine = rt.owned();
        mine.insert(mine.end(), rt.shared_not_owned().begin(), rt.shared_not_owned().end());
        std::sort(mine.begin(), mine.end());
        for (int32_t id : mine) {
            const bool owned = st.graph.comps[static_cast<size_t>(id)].owner == static_cast<int32_t>(st.my_pid);
            while (true) {
                owner_schedule_poll(st, false);
                bool got = false;
                std::vector<int32_t> stash;
                if (owned) {
                    while (rt.has_ready()) {
                        const int32_t r = rt.pop_ready(omp_get_wtime());
                        if (r == id) { got = true; break; }
                        stash.push_back(r);
                    }
                    for (int32_t r : stash) rt.requeue(r);
                } else {
                    while (rt.has_installable()) {
                        const int32_t r = rt.pop_installable();
                        if (r == id) { got = true; break; }
                        stash.push_back(r);
                    }
                    for (int32_t r : stash) rt.requeue_installable(r);
                }
                if (got) break;
                owner_schedule_poll(st, true);
            }
            auto b = owned ? owner_schedule_batch_begin(st, level, kernel, std::vector<int32_t>{id}, {})
                           : owner_schedule_batch_begin(st, level, kernel, {}, std::vector<int32_t>{id});
            while (wave(*b)) {}
            owner_schedule_batch_finish(st, *b);
        }
        while (!rt.boundary_complete()) owner_schedule_poll(st, true);
    } else if (st.sync_colors) {
        // Synchronous colour schedule (coarse levels, see above).  Per colour:
        // (1) wait until every owned component of the colour is runnable and
        //     eliminate them all in one batch (copies of the colour that have
        //     already arrived join it);
        // (2) wait until every remaining shared copy of the colour has arrived
        //     and install them all in one batch.
        // Both waits only depend on lower colours finishing elsewhere, so no
        // rank can wait on a rank that is waiting on it.  Same numerics as the
        // arrival-driven path (invariant 3: numerics do not depend on batching).
        auto take_color = [&](int32_t c, std::vector<int32_t>& owned, std::vector<int32_t>& installed) {
            std::vector<int32_t> stash;
            while (rt.has_ready()) { const int32_t r = rt.pop_ready(omp_get_wtime()); (color_of(r) == c ? owned : stash).push_back(r); }
            for (int32_t r : stash) rt.requeue(r);
            stash.clear();
            while (rt.has_installable()) { const int32_t r = rt.pop_installable(); (color_of(r) == c ? installed : stash).push_back(r); }
            for (int32_t r : stash) rt.requeue_installable(r);
            std::sort(owned.begin(), owned.end());
            std::sort(installed.begin(), installed.end());
        };
        auto give_back = [&](std::vector<int32_t>& owned, std::vector<int32_t>& installed) {
            for (int32_t r : owned) rt.requeue(r);
            for (int32_t r : installed) rt.requeue_installable(r);
            owned.clear();
            installed.clear();
        };
        for (int32_t c = 0; c < st.num_boundary_colors; ++c) {
            size_t want_owned = 0;
            for (int32_t id : rt.owned()) if (color_of(id) == c) ++want_owned;
            std::vector<int32_t> owned, installed;
            if (want_owned > 0) {
                while (true) {
                    owner_schedule_poll(st, false);
                    take_color(c, owned, installed);
                    if (owned.size() == want_owned) break;
                    give_back(owned, installed);
                    owner_schedule_poll(st, true);
                }
                auto b = owner_schedule_batch_begin(st, level, kernel, owned, installed);
                while (wave(*b)) {}
                owner_schedule_batch_finish(st, *b);
                owned.clear();
                installed.clear();
            }
            while (true) {
                owner_schedule_poll(st, false);
                size_t pending = 0;   // shared copies of colour c still to install
                for (int32_t id : rt.shared_not_owned()) if (color_of(id) == c && !rt.is_done(id)) ++pending;
                if (pending == 0) break;
                take_color(c, owned, installed);
                if (!owned.empty())
                    throw std::runtime_error("owner_schedule: owned component runnable after its colour's batch");
                if (installed.size() == pending) break;
                give_back(owned, installed);
                owner_schedule_poll(st, true);
            }
            if (!installed.empty()) {
                auto b = owner_schedule_batch_begin(st, level, kernel, {}, installed);
                while (wave(*b)) {}
                owner_schedule_batch_finish(st, *b);
            }
        }
        while (!rt.boundary_complete()) owner_schedule_poll(st, true);
        if (!all_done()) throw std::runtime_error("owner_schedule: synchronous schedule ended with work pending");
    } else {
        // Fast path (design §9): arrivals are serviced at every wave boundary;
        // the runnable set is drained one colour at a time, lowest colour
        // first (its results unblock the most); a runnable component of lower
        // colour than the batch in progress preempts it — the batch is
        // suspended after its current wave and resumed afterwards.  Safe:
        // a preempting component is never within two hops of the suspended
        // batch (it would be its predecessor), so neither pairs nor sketches
        // of the two interact.  Block only when nothing is runnable.
        std::vector<std::unique_ptr<Batch>> stack;
        auto best_runnable_color = [&]() -> int32_t {
            int32_t best = INT32_MAX;
            std::vector<int32_t> stash;
            while (rt.has_ready()) { const int32_t r = rt.pop_ready(omp_get_wtime()); best = std::min(best, color_of(r)); stash.push_back(r); }
            for (int32_t r : stash) rt.requeue(r);
            stash.clear();
            while (rt.has_installable()) { const int32_t r = rt.pop_installable(); best = std::min(best, color_of(r)); stash.push_back(r); }
            for (int32_t r : stash) rt.requeue_installable(r);
            return best;
        };
        auto drain_color = [&](int32_t c, std::vector<int32_t>& owned, std::vector<int32_t>& installed) {
            std::vector<int32_t> stash;
            while (rt.has_ready()) { const int32_t r = rt.pop_ready(omp_get_wtime()); (color_of(r) == c ? owned : stash).push_back(r); }
            for (int32_t r : stash) rt.requeue(r);
            stash.clear();
            while (rt.has_installable()) { const int32_t r = rt.pop_installable(); (color_of(r) == c ? installed : stash).push_back(r); }
            for (int32_t r : stash) rt.requeue_installable(r);
            std::sort(owned.begin(), owned.end());
            std::sort(installed.begin(), installed.end());
        };
        while (!(all_done() && stack.empty())) {
            owner_schedule_poll(st, false);
            const int32_t best = best_runnable_color();
            const bool preempt = !stack.empty() && best < stack.back()->color;
            if (best != INT32_MAX && (stack.empty() || preempt)) {
                std::vector<int32_t> owned, installed;
                drain_color(best, owned, installed);
                if (preempt) ++st.preemptions;
                stack.push_back(owner_schedule_batch_begin(st, level, kernel, owned, installed));
                continue;
            }
            if (!stack.empty()) {
                Batch& b = *stack.back();
                if (!wave(b)) {
                    owner_schedule_batch_finish(st, b);
                    stack.pop_back();
                }
                continue;
            }
            owner_schedule_poll(st, true);   // nothing runnable, nothing in progress
        }
    }
    if (!rt.interior_ready()) {
        rt.dump_pending(stdout);
        throw std::runtime_error("owner_schedule: boundary complete but interior counter not zero");
    }
    owner_schedule_final_event(st, level, kernel);   // keeps the X_NR buffers alive in st.xnr_hold
    // Drain the sends (they read X_NR in place).  Every peer keeps receiving
    // until it has consumed all payloads it expects, so the wait cannot deadlock.
    owner_schedule_progress_outbox(st, true);
    if (!st.pending_recv.empty()) throw std::runtime_error("owner_schedule: payload receives still pending at the end of the level");
    st.xnr_hold.clear();
    if (owner_allocator_keeps_freed_pages) OwnerBufferPool::instance().drain();
}

template<typename CoordType, typename DataType>
void owner_schedule_print_timing(const OwnerScheduleState<CoordType, DataType>& st, int level_index) {
    if (!st.active) return;
    const auto sm = st.rt.summary();
    std::printf("  [owner-schedule] level %d (%s, %d waves/colour): %lld batches (%lld preemptions), eliminated %lld comps (%lld boxes), installed %lld comps, "
                "sent %.1f MB, recv %.1f MB | ms: box-region %.1f, passes(all) %.1f, R2 replay %.1f, final %.1f, "
                "pack %.1f, install %.1f [thread-sum: X_NR regen %.1f, temp2 solve %.1f, block rebuild %.1f], idle-wait %.1f, landing-wait %.1f (%lld/%lld parts sent/recv) | "
                "ready->start wait max %.0f ms over %d comps | X_NR rows MB: shipped %.1f, regen kernel %.1f, regen partner %.1f | buffer pool %.1f GB\n",
                level_index, st.serial ? "serial oracle" : st.sync_colors ? "synchronous colours" : "arrival-driven",
                st.num_waves, static_cast<long long>(st.batches), static_cast<long long>(st.preemptions), static_cast<long long>(st.comps_eliminated),
                static_cast<long long>(st.boxes_eliminated), static_cast<long long>(st.comps_installed),
                st.bytes_sent / 1048576.0, st.bytes_recv / 1048576.0, st.t_elim_ms, st.t_passes_ms, st.t_replay_ms,
                st.t_final_ms, st.t_pack_ms, st.t_install_ms, st.t_rebuild_xnr_ms, st.t_rebuild_temp2_ms, st.t_rebuild_blocks_ms,
                st.t_wait_ms, st.t_recv_ms, static_cast<long long>(st.parts_sent), static_cast<long long>(st.parts_recv),
                1000.0 * sm.max_wait_after_ready, sm.n,
                st.xnr_bytes_shipped / 1048576.0, st.xnr_bytes_regen_kernel / 1048576.0, st.xnr_bytes_regen_partner / 1048576.0,
                OwnerBufferPool::instance().held_bytes() / (1024.0 * 1024.0 * 1024.0));
    std::fflush(stdout);
}

/// Collective: the slowest rank's view of the same fields (MPI_MAX per field,
/// so the numbers may come from different ranks).  Every rank calls it;
/// inactive ranks contribute zeros.
template<typename CoordType, typename DataType>
void owner_schedule_print_timing_max(const OwnerScheduleState<CoordType, DataType>& st, int level_index, MPI_Comm comm,
                                     int print_rank, int rank) {
    double v[12] = {0};
    if (st.active) {
        v[0] = st.t_elim_ms; v[1] = st.t_passes_ms; v[2] = st.t_replay_ms; v[3] = st.t_final_ms; v[4] = st.t_pack_ms;
        v[5] = st.t_install_ms; v[6] = st.t_wait_ms; v[7] = st.t_recv_ms; v[8] = st.bytes_sent / 1048576.0;
        v[9] = st.bytes_recv / 1048576.0; v[10] = static_cast<double>(st.boxes_eliminated); v[11] = static_cast<double>(st.comps_installed);
    }
    double m[12] = {0};
    MPI_Reduce(v, m, 12, MPI_DOUBLE, MPI_MAX, print_rank, comm);
    int any = st.active ? 1 : 0, any_all = 0;
    MPI_Reduce(&any, &any_all, 1, MPI_INT, MPI_MAX, print_rank, comm);
    if (rank == print_rank && any_all)
        std::printf("  [owner-schedule] level %d max over ranks: box-region %.1f, passes(all) %.1f, R2 replay %.1f, final %.1f, "
                    "pack %.1f, install %.1f, idle-wait %.1f, landing-wait %.1f ms | sent %.1f MB, recv %.1f MB | eliminated %.0f boxes, installed %.0f comps\n",
                    level_index, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
}

}  // namespace fmm
