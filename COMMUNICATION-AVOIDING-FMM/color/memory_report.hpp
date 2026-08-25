#ifndef MEMORY_REPORT_HPP
#define MEMORY_REPORT_HPP

/**
 * @file memory_report.hpp
 * @brief Per-level memory instrumentation for the H2 factorization.
 *
 * Purely additive: nothing here is required by the algorithm, and the whole
 * header can be deleted without touching factorization or solve code.
 *
 * Usage (inside the per-level loop of hierarchical_factorization_parallel):
 *
 *     report_level_memory(level, current_level, "post-elim",
 *                         tree->comm, level_print_rank, verbose);
 *
 * IMPORTANT: report_level_memory performs MPI collectives over `comm`.
 * EVERY rank in the communicator must call it, including ranks that are
 * inactive at this level. Never place a call inside an
 * `if (level.is_process_active)` guard.
 */

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "tree.hpp"

namespace fmm {

// ---------------------------------------------------------------------------
// Per-box breakdown
// ---------------------------------------------------------------------------

/**
 * @brief Memory footprint of one BoxData, split by category.
 *
 * All members are size_t byte counts and must stay that way: raw() reinterprets
 * the object as a flat array so the whole struct can go through a single
 * MPI_Reduce without a derived datatype.
 */
struct BoxMemoryBreakdown {
    size_t fixed       = 0;  ///< sizeof(BoxData): inline fields + container headers
    size_t points      = 0;  ///< point_indices, point_coords
    size_t indices     = 0;  ///< redundant/skeleton/one_hop/two_hop/use_full_set/pivots
    size_t interp      = 0;  ///< interpolation_matrix (T)
    size_t x_rr        = 0;  ///< X_RR
    size_t x_rs_sr     = 0;  ///< X_RS + X_SR
    size_t x_rn_nr     = 0;  ///< X_RN + X_NR + deferred_xnn_* scratch
    size_t schur       = 0;  ///< schur_complement
    size_t near_blocks = 0;  ///< near_field_modified_interactions (+ A_SN/A_NS heap)
    size_t far_blocks  = 0;  ///< far_field_modified_interactions  (+ A_SN/A_NS heap)
    size_t maps        = 0;  ///< all four interaction index maps

    static constexpr int N = 11;

    size_t total() const {
        return fixed + points + indices + interp + x_rr + x_rs_sr
             + x_rn_nr + schur + near_blocks + far_blocks + maps;
    }

    BoxMemoryBreakdown& operator+=(const BoxMemoryBreakdown& o) {
        fixed       += o.fixed;       points     += o.points;
        indices     += o.indices;     interp     += o.interp;
        x_rr        += o.x_rr;        x_rs_sr    += o.x_rs_sr;
        x_rn_nr     += o.x_rn_nr;     schur      += o.schur;
        near_blocks += o.near_blocks; far_blocks += o.far_blocks;
        maps        += o.maps;
        return *this;
    }

    size_t*       raw()       { return reinterpret_cast<size_t*>(this); }
    const size_t* raw() const { return reinterpret_cast<const size_t*>(this); }
};

static_assert(sizeof(BoxMemoryBreakdown) == BoxMemoryBreakdown::N * sizeof(size_t),
              "BoxMemoryBreakdown must remain a flat array of size_t for MPI_Reduce");
static_assert(sizeof(size_t) == sizeof(uint64_t),
              "size_t must match MPI_UINT64_T on this platform");

// ---------------------------------------------------------------------------
// Sizing helpers (mirrors calculate_box_data_size in tree.hpp)
// ---------------------------------------------------------------------------

namespace detail {

/// Heap bytes owned by a MatrixStorage (the vector allocation only).
template<typename DataType>
inline size_t mr_matrix_heap(const MatrixStorage<DataType>& m) {
    return m.data.capacity() * sizeof(DataType);
}

/// Heap bytes owned by the element storage of a vector (header is inline).
template<typename Vec>
inline size_t mr_vec_heap(const Vec& v) {
    return v.capacity() * sizeof(typename Vec::value_type);
}

/// Rough heap estimate for an unordered_map / unordered_set:
/// bucket array plus a per-element node (value + next pointer + cached hash).
template<typename Map>
inline size_t mr_map_heap(const Map& m) {
    using value_type = typename Map::value_type;
    return m.bucket_count() * sizeof(void*)
         + m.size() * (sizeof(value_type) + sizeof(void*) + sizeof(size_t));
}

/// Heap inside one ModifiedBlock; its sizeof is covered by the vector capacity.
template<typename DataType>
inline size_t mr_block_heap(const ModifiedBlock<DataType>& mb) {
    return mr_matrix_heap(mb.A_SN) + mr_matrix_heap(mb.A_NS);
}

template<typename DataType>
inline size_t mr_blocks_total(const std::vector<ModifiedBlock<DataType>>& v) {
    size_t total = mr_vec_heap(v);
    for (const auto& mb : v) total += mr_block_heap(mb);
    return total;
}

} // namespace detail

/**
 * @brief Categorized memory footprint of a single box.
 *
 * Superset of calculate_box_data_size(): additionally counts
 * deferred_xnn_neighbor_point_counts and deferred_xnn_temp2, which are live
 * during elimination but omitted by the original function.
 */
template<typename CoordType, typename DataType>
BoxMemoryBreakdown calculate_box_memory_breakdown(const BoxData<CoordType, DataType>& box) {
    using namespace detail;
    BoxMemoryBreakdown b;

    b.fixed = sizeof(BoxData<CoordType, DataType>);

    b.points = mr_vec_heap(box.point_indices)
             + mr_vec_heap(box.point_coords);

    b.indices = mr_vec_heap(box.redundant_indices)
              + mr_vec_heap(box.skeleton_indices)
              + mr_vec_heap(box.one_hop)
              + mr_vec_heap(box.two_hop)
              + mr_vec_heap(box.use_full_set)
              + mr_vec_heap(box.X_RR_pivots);

    b.interp = mr_matrix_heap(box.interpolation_matrix);
    b.x_rr   = mr_matrix_heap(box.X_RR);

    b.x_rs_sr = mr_matrix_heap(box.X_RS)
              + mr_matrix_heap(box.X_SR);

    b.x_rn_nr = mr_matrix_heap(box.X_RN)
              + mr_matrix_heap(box.X_NR)
              + mr_vec_heap(box.deferred_xnn_neighbor_point_counts)
              + mr_vec_heap(box.deferred_xnn_temp2);

    b.schur = mr_matrix_heap(box.schur_complement);

    b.near_blocks = mr_blocks_total(box.near_field_modified_interactions);
    b.far_blocks  = mr_blocks_total(box.far_field_modified_interactions);

    b.maps = mr_map_heap(box.near_field_interaction_map)
           + mr_map_heap(box.far_field_interaction_map)
           + mr_map_heap(box.near_field_interaction_map_nonsymmetry)
           + mr_map_heap(box.far_field_interaction_map_nonsymmetry);

    return b;
}

/**
 * @brief Bytes held by the containers owned by TreeLevel itself, not by any box.
 *
 * Per-box accounting cannot see these; without them a level total is an unknown
 * fraction of the real footprint. box_locks additionally owns one
 * heap-allocated omp_lock_t per entry.
 */
template<typename CoordType, typename DataType>
size_t level_overhead_bytes(const TreeLevel<CoordType, DataType>& level) {
    using namespace detail;
    size_t total = 0;

    // Color ordering
    total += mr_vec_heap(level.blue)   + mr_vec_heap(level.orange)
           + mr_vec_heap(level.purple) + mr_vec_heap(level.green);

    // Box id partitions
    total += mr_vec_heap(level.ghost_id)
           + mr_vec_heap(level.interior_id)
           + mr_vec_heap(level.boundary_id);

    // Lookup maps
    total += mr_map_heap(level.assisting_box_points_for_kernel_evaluation)
           + mr_map_heap(level.ghost_id_to_index)
           + mr_map_heap(level.ghost_and_assisting_box_points_for_solve_map)
           + mr_map_heap(level.morton_to_rank)
           + mr_map_heap(level.rank_to_morton)
           + mr_map_heap(level.eliminated_boxes);

    // Per-box OpenMP locks: map node plus the heap-allocated lock itself
    total += mr_map_heap(level.box_locks)
           + level.box_locks.size() * sizeof(omp_lock_t);

    // Solve-phase bookkeeping
    total += mr_vec_heap(level.is_ghost_solve)
           + mr_vec_heap(level.active_process_ranks)
           + mr_vec_heap(level.solve_neighbor_size);
    for (const auto& v : level.solve_neighbor_size) total += mr_vec_heap(v);

    return total;
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

inline std::string mr_fmt_bytes(size_t bytes) {
    char buf[64];
    const double x = static_cast<double>(bytes);
    if      (x >= 1073741824.0) std::snprintf(buf, sizeof buf, "%.2f GB", x / 1073741824.0);
    else if (x >= 1048576.0)    std::snprintf(buf, sizeof buf, "%.2f MB", x / 1048576.0);
    else if (x >= 1024.0)       std::snprintf(buf, sizeof buf, "%.2f KB", x / 1024.0);
    else                        std::snprintf(buf, sizeof buf, "%zu B", bytes);
    return std::string(buf);
}

/**
 * @brief Reduce and print per-level memory usage.
 *
 * COLLECTIVE over `comm`: all ranks must call it, active or not. Inactive ranks
 * contribute zeros naturally because local_boxes is empty.
 *
 * @param tag short label for the sample point, e.g. "post-elim" / "post-clear"
 */
template<typename CoordType, typename DataType>
void report_level_memory(const TreeLevel<CoordType, DataType>& level,
                         int level_idx,
                         const char* tag,
                         MPI_Comm comm,
                         int print_rank,
                         bool verbose)
{
    constexpr size_t kNoMin = std::numeric_limits<size_t>::max();

    BoxMemoryBreakdown local_sum;
    size_t local_max = 0;
    size_t local_min = kNoMin;

    for (const auto& box : level.local_boxes) {
        const BoxMemoryBreakdown bd = calculate_box_memory_breakdown(box);
        local_sum += bd;
        const size_t t = bd.total();
        if (t > local_max) local_max = t;
        if (t < local_min) local_min = t;
    }

    constexpr int NPACK = BoxMemoryBreakdown::N + 2;  // + level overhead + box count
    size_t sendbuf[NPACK];
    std::copy(local_sum.raw(), local_sum.raw() + BoxMemoryBreakdown::N, sendbuf);
    sendbuf[BoxMemoryBreakdown::N + 0] = level_overhead_bytes(level);
    sendbuf[BoxMemoryBreakdown::N + 1] = level.local_boxes.size();

    size_t recvbuf[NPACK];
    std::fill(recvbuf, recvbuf + NPACK, size_t{0});
    size_t global_max = 0;
    size_t global_min = 0;

    // ---- collectives: every rank in comm must reach these ----
    MPI_Reduce(sendbuf, recvbuf, NPACK, MPI_UINT64_T, MPI_SUM, print_rank, comm);
    MPI_Reduce(&local_max, &global_max, 1, MPI_UINT64_T, MPI_MAX, print_rank, comm);
    MPI_Reduce(&local_min, &global_min, 1, MPI_UINT64_T, MPI_MIN, print_rank, comm);

    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank != print_rank || !verbose) return;

    const size_t overhead  = recvbuf[BoxMemoryBreakdown::N + 0];
    const size_t box_count = recvbuf[BoxMemoryBreakdown::N + 1];

    size_t box_total = 0;
    for (int i = 0; i < BoxMemoryBreakdown::N; ++i) box_total += recvbuf[i];
    const size_t grand_total = box_total + overhead;

    static const char* const kNames[BoxMemoryBreakdown::N] = {
        "fixed", "points", "indices", "interp", "X_RR", "X_RS/X_SR",
        "X_RN/X_NR", "schur", "near_blocks", "far_blocks", "maps"
    };

    std::vector<std::pair<size_t, const char*>> rows;
    rows.reserve(BoxMemoryBreakdown::N);
    for (int i = 0; i < BoxMemoryBreakdown::N; ++i) rows.emplace_back(recvbuf[i], kNames[i]);
    std::sort(rows.begin(), rows.end(),
              [](const std::pair<size_t, const char*>& a,
                 const std::pair<size_t, const char*>& b) { return a.first > b.first; });

    std::printf("\n=== Level %d  [%s]  %zu boxes (all ranks) ===\n",
                level_idx, tag, box_count);

    if (box_count > 0) {
        const size_t mean = box_total / box_count;
        const std::string min_str =
            (global_min == kNoMin) ? std::string("n/a") : mr_fmt_bytes(global_min);
        std::printf("  per-box:  min %s   max %s   mean %s\n",
                    min_str.c_str(),
                    mr_fmt_bytes(global_max).c_str(),
                    mr_fmt_bytes(mean).c_str());
    } else {
        std::printf("  per-box:  n/a (no boxes)\n");
    }

    for (int i = 0; i < BoxMemoryBreakdown::N; ++i) {
        if (rows[i].first == 0) continue;
        const double pct = (grand_total > 0)
            ? 100.0 * static_cast<double>(rows[i].first) / static_cast<double>(grand_total)
            : 0.0;
        std::printf("  %-12s %12s  %5.1f%%\n",
                    rows[i].second, mr_fmt_bytes(rows[i].first).c_str(), pct);
    }
    if (overhead > 0) {
        const double pct = (grand_total > 0)
            ? 100.0 * static_cast<double>(overhead) / static_cast<double>(grand_total)
            : 0.0;
        std::printf("  %-12s %12s  %5.1f%%\n",
                    "level_overhd", mr_fmt_bytes(overhead).c_str(), pct);
    }
    std::printf("  %-12s %12s\n", "TOTAL", mr_fmt_bytes(grand_total).c_str());
    std::fflush(stdout);
}

} // namespace fmm

#endif // MEMORY_REPORT_HPP
