#pragma once

#include "gpu_common/bpack_env.hpp"
#include "factorization.hpp"  // the deferred X_NN types
#include "tree.hpp"

#include <algorithm>
#include <cstdio>
#include <functional>
#include <cstdlib>
#include <string>
#include <sys/resource.h>
#include <unordered_set>
#include <unistd.h>
#include <vector>

namespace butterfly {

struct H2MemoryCategories {
    size_t local_metadata = 0;
    size_t halo_metadata = 0;
    size_t local_factors = 0;
    size_t halo_factors = 0;
    size_t local_near = 0;
    size_t halo_near = 0;
    size_t local_far = 0;
    size_t halo_far = 0;
    size_t local_deferred = 0;
    size_t halo_deferred = 0;
    size_t h2_blocks = 0;
    size_t assisting = 0;
    size_t solve = 0;
    size_t level_metadata = 0;
    size_t pending = 0;
    size_t scratch = 0;
    size_t staged_parent = 0;
    size_t communication = 0;
    size_t gpu_host_staging = 0;

    // A breakdown of local_factors / halo_factors by the role each block plays
    // in the solve, reported on its own H2_MEM_FACTOR line.  These are already
    // counted in local_factors / halo_factors, so total() must NOT add them.
    size_t local_factor_pivot = 0;
    size_t halo_factor_pivot = 0;
    size_t local_factor_basis = 0;
    size_t halo_factor_basis = 0;
    size_t local_factor_skeleton = 0;
    size_t halo_factor_skeleton = 0;
    size_t local_factor_near_coupling = 0;
    size_t halo_factor_near_coupling = 0;
    size_t local_factor_schur = 0;
    size_t halo_factor_schur = 0;
    // Which array the near coupling sits in, and whether the box is a ghost
    // (a replica of another rank's box) or a lazy generator (regeneration
    // state, which may be for the rank's own boxes).  Already counted above.
    size_t near_x_nr = 0;
    size_t near_x_rn = 0;
    size_t near_lazy_original = 0;
    size_t ghost_factors_only = 0;
    size_t generator_factors_only = 0;
    size_t generator_near_coupling = 0;

    size_t total() const {
        return local_metadata + halo_metadata + local_factors + halo_factors +
               local_near + halo_near + local_far + halo_far +
               local_deferred + halo_deferred + h2_blocks + assisting + solve +
               level_metadata + pending + scratch + staged_parent +
               communication + gpu_host_staging;
    }
};

struct H2MemorySnapshot {
    int level = -1;
    std::string algorithm;
    std::string phase;
    H2MemoryCategories categories;
    size_t current_rss = 0;
    size_t high_water_rss = 0;
};

// Pinned host bytes of the GPU backend, which live outside every tree
// structure.  The GPU path installs this once (core must not include the
// backend's headers); without it the bytes land in `unattributed`.
inline std::function<size_t()>& h2_diag_host_staging_probe() {
    static std::function<size_t()> probe;
    return probe;
}

template<typename VectorType>
inline size_t h2_diag_vector_bytes(const VectorType& values) {
    return values.capacity() * sizeof(typename VectorType::value_type);
}

template<typename MapType>
inline size_t h2_diag_map_bytes(const MapType& values) {
    using value_type = typename MapType::value_type;
    return values.bucket_count() * sizeof(void*) +
           values.size() * (sizeof(value_type) + sizeof(void*) + sizeof(size_t));
}

template<typename DataType>
inline size_t h2_diag_matrix_bytes(const fmm::MatrixStorage<DataType>& matrix) {
    return matrix.data.capacity() * sizeof(DataType);
}

template<typename DataType>
inline size_t h2_diag_modified_block_bytes(
    const std::vector<fmm::ModifiedBlock<DataType>>& blocks) {
    size_t bytes = h2_diag_vector_bytes(blocks);
    for (const auto& block : blocks) {
        bytes += block.a_ns_heap_bytes();
        bytes += h2_diag_matrix_bytes(block.A_SN);
    }
    return bytes;
}

template<typename DataType>
inline size_t h2_diag_h2_block_bytes(
    const std::vector<fmm::H2Block<DataType>>& blocks) {
    size_t bytes = h2_diag_vector_bytes(blocks);
    for (const auto& block : blocks) {
        bytes += h2_diag_matrix_bytes(block.matrix);
    }
    return bytes;
}

template<typename CoordType, typename DataType>
inline void h2_diag_accumulate_box(
    const fmm::BoxData<CoordType, DataType>& box,
    bool halo,
    H2MemoryCategories& categories,
    bool is_generator = false) {
    size_t metadata = sizeof(box);
    metadata += h2_diag_vector_bytes(box.point_indices);
    metadata += h2_diag_vector_bytes(box.point_coords);
    metadata += h2_diag_vector_bytes(box.redundant_indices);
    metadata += h2_diag_vector_bytes(box.skeleton_indices);
    metadata += h2_diag_vector_bytes(box.one_hop);
    metadata += h2_diag_vector_bytes(box.use_full_set);
    metadata += h2_diag_vector_bytes(box.two_hop);

    // The factor bucket split by solve role: the pivot blocks drive the
    // triangular solves, the ID basis and the coupling blocks are the bulk.
    // The five sub-totals add up to `factors` by construction.
    const size_t factor_pivot =
        h2_diag_vector_bytes(box.X_RR_pivots) +
        h2_diag_matrix_bytes(box.X_RR) +
        h2_diag_matrix_bytes(box.X_RR_full);
    const size_t factor_basis = h2_diag_matrix_bytes(box.interpolation_matrix);
    const size_t factor_skeleton =
        h2_diag_matrix_bytes(box.X_RS) +
        h2_diag_matrix_bytes(box.X_SR) +
        h2_diag_matrix_bytes(box.X_RS_entry);
    const size_t near_x_rn = h2_diag_matrix_bytes(box.X_RN);
    const size_t near_x_nr = h2_diag_matrix_bytes(box.X_NR);
    // The lazy generator's copy of the pre-finalize X_NR.
    const size_t near_lazy_original =
        h2_diag_vector_bytes(box.lazy_original_x_nr);
    const size_t factor_near_coupling =
        near_x_rn + near_x_nr + near_lazy_original;
    const size_t factor_schur = h2_diag_matrix_bytes(box.schur_complement);
    const size_t factors = factor_pivot + factor_basis + factor_skeleton +
                           factor_near_coupling + factor_schur;

    categories.near_x_rn += near_x_rn;
    categories.near_x_nr += near_x_nr;
    categories.near_lazy_original += near_lazy_original;
    if (is_generator) {
        categories.generator_factors_only += factors;
        categories.generator_near_coupling += factor_near_coupling;
    } else if (halo) {
        categories.ghost_factors_only += factors;
    }

    size_t near = h2_diag_modified_block_bytes(
        box.near_field_modified_interactions);
    near += h2_diag_map_bytes(box.near_field_interaction_map);
    near += h2_diag_map_bytes(box.near_field_interaction_map_nonsymmetry);

    size_t far = h2_diag_modified_block_bytes(
        box.far_field_modified_interactions);
    far += h2_diag_map_bytes(box.far_field_interaction_map);
    far += h2_diag_map_bytes(box.far_field_interaction_map_nonsymmetry);

    const size_t deferred =
        h2_diag_vector_bytes(box.deferred_xnn_neighbor_point_counts) +
        h2_diag_vector_bytes(box.deferred_xnn_temp2);

    categories.h2_blocks += h2_diag_h2_block_bytes(box.h2_interaction_blocks);
    categories.h2_blocks += h2_diag_h2_block_bytes(box.h2_near_blocks);

    if (halo) {
        categories.halo_metadata += metadata;
        categories.halo_factors += factors;
        categories.halo_factor_pivot += factor_pivot;
        categories.halo_factor_basis += factor_basis;
        categories.halo_factor_skeleton += factor_skeleton;
        categories.halo_factor_near_coupling += factor_near_coupling;
        categories.halo_factor_schur += factor_schur;
        categories.halo_near += near;
        categories.halo_far += far;
        categories.halo_deferred += deferred;
    } else {
        categories.local_metadata += metadata;
        categories.local_factors += factors;
        categories.local_factor_pivot += factor_pivot;
        categories.local_factor_basis += factor_basis;
        categories.local_factor_skeleton += factor_skeleton;
        categories.local_factor_near_coupling += factor_near_coupling;
        categories.local_factor_schur += factor_schur;
        categories.local_near += near;
        categories.local_far += far;
        categories.local_deferred += deferred;
    }
}

template<typename CoordType>
inline size_t h2_diag_assisting_request_bytes(
    const fmm::PointDataRequest<CoordType>& request) {
    return sizeof(request) + h2_diag_vector_bytes(request.coords) +
           h2_diag_vector_bytes(request.indices) +
           h2_diag_vector_bytes(request.skel_indices);
}

template<typename CoordType, typename DataType>
inline size_t h2_diag_solve_request_bytes(
    const fmm::SolveDataRequest<CoordType, DataType>& request) {
    size_t bytes = sizeof(request);
    bytes += h2_diag_vector_bytes(request.right_side);
    bytes += h2_diag_vector_bytes(request.left_side);
    bytes += h2_diag_vector_bytes(request.redundant_indices);
    bytes += h2_diag_vector_bytes(request.skeleton_indices);
    bytes += h2_diag_vector_bytes(request.one_hop);
    bytes += h2_diag_vector_bytes(request.use_full_set);
    bytes += h2_diag_vector_bytes(request.X_RR_pivots);
    bytes += h2_diag_matrix_bytes(request.interpolation_matrix);
    bytes += h2_diag_matrix_bytes(request.X_RR);
    bytes += h2_diag_matrix_bytes(request.X_RS);
    bytes += h2_diag_matrix_bytes(request.X_SR);
    bytes += h2_diag_matrix_bytes(request.schur_complement);
    bytes += h2_diag_matrix_bytes(request.X_RN);
    bytes += h2_diag_matrix_bytes(request.X_NR);
    return bytes;
}

template<typename CoordType, typename DataType>
inline void h2_diag_accumulate_level_metadata(
    const fmm::TreeLevel<CoordType, DataType>& level,
    H2MemoryCategories& categories) {
    size_t bytes = sizeof(level);
    bytes += h2_diag_vector_bytes(level.blue);
    bytes += h2_diag_vector_bytes(level.orange);
    bytes += h2_diag_vector_bytes(level.purple);
    bytes += h2_diag_vector_bytes(level.green);
    bytes += h2_diag_vector_bytes(level.ghost_id);
    bytes += h2_diag_vector_bytes(level.interior_id);
    bytes += h2_diag_vector_bytes(level.boundary_id);
    bytes += h2_diag_vector_bytes(level.active_process_ranks);
    bytes += h2_diag_vector_bytes(level.children_senders);
    bytes += h2_diag_vector_bytes(level.is_ghost_solve);
    bytes += h2_diag_map_bytes(level.assisting_box_points_for_kernel_evaluation);
    bytes += h2_diag_map_bytes(level.ghost_id_to_index);
    bytes += h2_diag_map_bytes(level.box_locks);
    bytes += h2_diag_map_bytes(level.eliminated_boxes);
    bytes += h2_diag_map_bytes(level.ghost_and_assisting_box_points_for_solve_map);
    bytes += h2_diag_map_bytes(level.morton_to_rank);
    bytes += h2_diag_map_bytes(level.rank_to_morton);
    bytes += h2_diag_map_bytes(level.elimination_wave);
    bytes += h2_diag_map_bytes(level.generator_id_to_index);
    bytes += h2_diag_vector_bytes(level.solve_neighbor_size);
    for (const auto& neighbor_sizes : level.solve_neighbor_size) {
        bytes += h2_diag_vector_bytes(neighbor_sizes);
    }
    categories.level_metadata += bytes;

    categories.assisting += h2_diag_vector_bytes(level.assisting_boxes);
    for (const auto& request : level.assisting_boxes) {
        categories.assisting += h2_diag_assisting_request_bytes(request) - sizeof(request);
    }

    categories.solve += h2_diag_vector_bytes(
        level.ghost_and_assisting_boxes_for_solve);
    for (const auto& request : level.ghost_and_assisting_boxes_for_solve) {
        categories.solve += h2_diag_solve_request_bytes(request) - sizeof(request);
    }
}

template<typename DataType>
inline size_t h2_diag_pending_bytes(
    const fmm::PendingFactorUpdates<DataType>& pending) {
    size_t bytes = h2_diag_map_bytes(pending.replace_blocks);
    for (const auto& entry : pending.replace_blocks) {
        bytes += h2_diag_vector_bytes(entry.second.data);
    }
    bytes += h2_diag_map_bytes(pending.accumulated_deltas);
    for (const auto& entry : pending.accumulated_deltas) {
        bytes += h2_diag_vector_bytes(entry.second.data);
    }
    bytes += h2_diag_map_bytes(pending.generators);
    std::unordered_set<const fmm::GeneratorPayload<DataType>*> seen;
    for (const auto& entry : pending.generators) {
        const auto* generator = entry.second.get();
        if (generator == nullptr || !seen.insert(generator).second) continue;
        bytes += sizeof(*generator);
        bytes += h2_diag_vector_bytes(generator->one_hop);
        bytes += h2_diag_vector_bytes(generator->neighbor_point_counts);
        bytes += h2_diag_vector_bytes(generator->temp2);
        bytes += h2_diag_vector_bytes(generator->x_rr_full);
        bytes += h2_diag_vector_bytes(generator->skeleton_indices);
        bytes += h2_diag_vector_bytes(generator->x_rs);
    }
    return bytes;
}

template<typename CoordType, typename DataType>
inline size_t h2_diag_factorization_scratch_bytes(
    const fmm::FactorizationThreadScratch<CoordType, DataType>& scratch) {
    return h2_diag_vector_bytes(scratch.workspace) +
           h2_diag_vector_bytes(scratch.x_nn_full) +
           h2_diag_vector_bytes(scratch.sketch_storage) +
           h2_diag_vector_bytes(scratch.x_bb) +
           h2_diag_vector_bytes(scratch.neighbor_point_counts) +
           h2_diag_vector_bytes(scratch.a_ns_all) +
           h2_diag_vector_bytes(scratch.a_sn_all) +
           h2_diag_vector_bytes(scratch.temp1) +
           h2_diag_vector_bytes(scratch.temp2) +
           h2_diag_vector_bytes(scratch.temp3) +
           h2_diag_vector_bytes(scratch.temp4) +
           h2_diag_vector_bytes(scratch.x_rs_original) +
           h2_diag_vector_bytes(scratch.x_ns_update) +
           h2_diag_vector_bytes(scratch.x_sn_update) +
           h2_diag_vector_bytes(scratch.update_buffer) +
           h2_diag_vector_bytes(scratch.eval_buffer) +
           h2_diag_vector_bytes(scratch.coord_buffer) +
           h2_diag_vector_bytes(scratch.index_buffer) +
           h2_diag_vector_bytes(scratch.stream_block) +
           h2_diag_vector_bytes(scratch.stream_sketch_acc) +
           h2_diag_vector_bytes(scratch.stream_stored_tmp) +
           h2_diag_vector_bytes(scratch.stream_counts) +
           h2_diag_vector_bytes(scratch.stream_row_indices) +
           h2_diag_vector_bytes(scratch.stream_pair_source_mortons) +
           h2_diag_vector_bytes(scratch.stream_pair_source_offsets) +
           h2_diag_vector_bytes(scratch.stream_box_source_mortons) +
           h2_diag_vector_bytes(scratch.stream_box_source_boxes) +
           h2_diag_vector_bytes(scratch.stream_box_source_r) +
           h2_diag_vector_bytes(scratch.stream_P_all) +
           h2_diag_vector_bytes(scratch.stream_Trow) +
           h2_diag_vector_bytes(scratch.stream_TrowT) +
           h2_diag_vector_bytes(scratch.stream_W) +
           h2_diag_vector_bytes(scratch.stream_row_idx) +
           h2_diag_vector_bytes(scratch.stream_row_sign) +
           h2_diag_vector_bytes(scratch.stream_blk_draw_idx) +
           h2_diag_vector_bytes(scratch.stream_blk_draw_sign) +
           h2_diag_vector_bytes(scratch.stream_blk_row_base) +
           h2_diag_vector_bytes(scratch.stream_blk_full_size) +
           h2_diag_vector_bytes(scratch.stream_blk_skeleton) +
           h2_diag_vector_bytes(scratch.stream_blk_wanted) +
           h2_diag_vector_bytes(scratch.stream_positions);
}

template<typename DataType>
inline size_t h2_diag_deferred_owner_scratch_bytes(
    const fmm::DeferredXnnOwnerScratch<DataType>& scratch) {
    size_t bytes = h2_diag_vector_bytes(scratch.packed_temp2_rows);
    bytes += h2_diag_vector_bytes(scratch.packed_updates);
    bytes += h2_diag_vector_bytes(scratch.owned_row_blocks);
    bytes += h2_diag_vector_bytes(scratch.accumulated_targets);
    for (const auto& target : scratch.accumulated_targets) {
        bytes += h2_diag_vector_bytes(target.data);
    }
    bytes += h2_diag_map_bytes(scratch.accumulated_target_indices);
    bytes += h2_diag_map_bytes(scratch.preallocated_mirror_targets);
    return bytes;
}

template<typename CoordType, typename DataType>
inline size_t h2_diag_box_vector_bytes(
    const std::vector<fmm::BoxData<CoordType, DataType>>& boxes) {
    H2MemoryCategories categories;
    for (const auto& box : boxes) {
        h2_diag_accumulate_box(box, false, categories);
    }
    return categories.total();
}

inline size_t h2_diag_current_rss_bytes() {
    FILE* statm = std::fopen("/proc/self/statm", "r");
    if (statm == nullptr) return 0;
    unsigned long total_pages = 0;
    unsigned long resident_pages = 0;
    const int parsed = std::fscanf(statm, "%lu %lu", &total_pages, &resident_pages);
    std::fclose(statm);
    if (parsed != 2) return 0;
    const long page_size = sysconf(_SC_PAGESIZE);
    return page_size > 0 ? static_cast<size_t>(resident_pages) *
                               static_cast<size_t>(page_size)
                         : 0;
}

inline size_t h2_diag_high_water_rss_bytes() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
    return static_cast<size_t>(usage.ru_maxrss);
#else
    return static_cast<size_t>(usage.ru_maxrss) * 1024;
#endif
}

class H2FactorizationMemoryDiagnostics {
public:
    // BPACK_TRACE=memory
    H2FactorizationMemoryDiagnostics() : enabled_(fmm::env::trace("memory")) {}

    bool enabled() const { return enabled_; }

    template<typename CoordType, typename DataType>
    void record(
        const fmm::ParallelTree<CoordType, DataType>* tree,
        int level,
        const std::string& algorithm,
        const std::string& phase,
        size_t pending = 0,
        size_t scratch = 0,
        size_t staged_parent = 0,
        size_t communication = 0) {
        if (!enabled_) return;

        H2MemorySnapshot snapshot;
        snapshot.level = level;
        snapshot.algorithm = algorithm;
        snapshot.phase = phase;
        snapshot.categories.pending = pending;
        snapshot.categories.scratch = scratch;
        snapshot.categories.staged_parent = staged_parent;
        snapshot.categories.communication = communication;

        for (const auto& tree_level : tree->levels) {
            h2_diag_accumulate_level_metadata(tree_level, snapshot.categories);
            for (const auto& box : tree_level.local_boxes) {
                h2_diag_accumulate_box(box, false, snapshot.categories);
            }
            for (const auto& box : tree_level.ghost_boxes) {
                h2_diag_accumulate_box(box, true, snapshot.categories);
            }
            for (const auto& box : tree_level.generator_boxes) {
                h2_diag_accumulate_box(box, true, snapshot.categories, true);
            }
        }

        const auto& host_staging_probe = h2_diag_host_staging_probe();
        snapshot.categories.gpu_host_staging =
            host_staging_probe ? host_staging_probe() : 0;

        snapshot.current_rss = h2_diag_current_rss_bytes();
        snapshot.high_water_rss = h2_diag_high_water_rss_bytes();
        snapshots_.push_back(std::move(snapshot));
    }

    void print(MPI_Comm comm, int rank, int size) const {
        if (!enabled_) return;

        for (int output_rank = 0; output_rank < size; ++output_rank) {
            MPI_Barrier(comm);
            if (rank != output_rank) continue;

            for (size_t sequence = 0; sequence < snapshots_.size(); ++sequence) {
                const auto& snapshot = snapshots_[sequence];
                const auto& c = snapshot.categories;
                const size_t accounted = c.total();
                const size_t unattributed = snapshot.current_rss > accounted
                    ? snapshot.current_rss - accounted
                    : 0;
                constexpr double gib = 1024.0 * 1024.0 * 1024.0;
                std::printf(
                    "H2_MEM_DIAG rank=%d seq=%zu level=%d algorithm=%s phase=%s "
                    "rss=%.6f hwm=%.6f accounted=%.6f unattributed=%.6f "
                    "local_meta=%.6f halo_meta=%.6f local_factor=%.6f halo_factor=%.6f "
                    "local_near=%.6f halo_near=%.6f local_far=%.6f halo_far=%.6f "
                    "local_deferred=%.6f halo_deferred=%.6f h2_blocks=%.6f "
                    "assisting=%.6f solve=%.6f level_meta=%.6f pending=%.6f "
                    "scratch=%.6f staged_parent=%.6f communication=%.6f "
                    "gpu_host_staging=%.6f\n",
                    rank, sequence, snapshot.level,
                    snapshot.algorithm.c_str(), snapshot.phase.c_str(),
                    snapshot.current_rss / gib, snapshot.high_water_rss / gib,
                    accounted / gib, unattributed / gib,
                    c.local_metadata / gib, c.halo_metadata / gib,
                    c.local_factors / gib, c.halo_factors / gib,
                    c.local_near / gib, c.halo_near / gib,
                    c.local_far / gib, c.halo_far / gib,
                    c.local_deferred / gib, c.halo_deferred / gib,
                    c.h2_blocks / gib, c.assisting / gib, c.solve / gib,
                    c.level_metadata / gib, c.pending / gib, c.scratch / gib,
                    c.staged_parent / gib, c.communication / gib,
                    c.gpu_host_staging / gib);
                // How local_factors / halo_factors above split by solve role.
                std::printf(
                    "H2_MEM_FACTOR rank=%d seq=%zu level=%d algorithm=%s phase=%s "
                    "local_pivot=%.6f local_basis=%.6f local_skeleton=%.6f "
                    "local_near_coupling=%.6f local_schur=%.6f "
                    "halo_pivot=%.6f halo_basis=%.6f halo_skeleton=%.6f "
                    "halo_near_coupling=%.6f halo_schur=%.6f "
                    "near_x_nr=%.6f near_x_rn=%.6f near_lazy_original=%.6f "
                    "ghost_factors_only=%.6f generator_factors_only=%.6f "
                    "generator_near_coupling=%.6f\n",
                    rank, sequence, snapshot.level,
                    snapshot.algorithm.c_str(), snapshot.phase.c_str(),
                    c.local_factor_pivot / gib, c.local_factor_basis / gib,
                    c.local_factor_skeleton / gib,
                    c.local_factor_near_coupling / gib, c.local_factor_schur / gib,
                    c.halo_factor_pivot / gib, c.halo_factor_basis / gib,
                    c.halo_factor_skeleton / gib,
                    c.halo_factor_near_coupling / gib, c.halo_factor_schur / gib,
                    c.near_x_nr / gib, c.near_x_rn / gib,
                    c.near_lazy_original / gib, c.ghost_factors_only / gib,
                    c.generator_factors_only / gib,
                    c.generator_near_coupling / gib);
            }
            std::fflush(stdout);
        }
        MPI_Barrier(comm);
    }

private:
    bool enabled_ = false;
    std::vector<H2MemorySnapshot> snapshots_;
};

} // namespace butterfly
