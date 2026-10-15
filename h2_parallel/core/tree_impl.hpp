#ifndef TREE_IMPL_HPP
#define TREE_IMPL_HPP

#include "tree.hpp"
#include "halo_geometry.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <ostream>
#include <iostream>
#include <set>
#include <malloc.h>

namespace fmm {



// ===== BoxData Implementation =====

template<typename CoordType, typename DataType>
BoxData<CoordType, DataType>::BoxData() 
    : morton_index(-1), level(-1), size(0.0),
      parent_morton(-1), num_children(0),
      num_points(0),
      num_near_field_interactions(0),
      num_far_field_interactions(0) {
    
    grid_coords[0] = grid_coords[1] = grid_coords[2] = -1;
    center[0] = center[1] = center[2] = 0.0;
    for (int i = 0; i < 6; ++i) bounds[i] = 0.0;
    for (int i = 0; i < 8; ++i) children_morton[i] = -1;
}

// Destructor now handled by default - vectors clean up automatically

// ===== TreeLevel Implementation =====

template<typename CoordType, typename DataType>
TreeLevel<CoordType, DataType>::TreeLevel()
    : level(-1), num_boxes_global(0), num_boxes_local(0),
      local_morton_start(-1), local_morton_end(-1),
      num_active_processes(0),
      is_process_active(false), my_morton_id(-1), parent_level_owner(-1) {
}

// Destructor now handled by default - vectors clean up automatically

template<typename CoordType, typename DataType>
BoxData<CoordType, DataType>* TreeLevel<CoordType, DataType>::find_local_box(int64_t morton_index) {
    // Simple O(1) bracket indexing for contiguous Morton range
    if (morton_index >= local_morton_start && morton_index <= local_morton_end) {
        int64_t offset = morton_index - local_morton_start;
        return &local_boxes[offset];
    }
    return nullptr;
}

template<typename CoordType, typename DataType>
BoxData<CoordType, DataType>* TreeLevel<CoordType, DataType>::find_ghost_box(int64_t morton_index) {

    auto it = ghost_id_to_index.find(morton_index);

    if (it != ghost_id_to_index.end()) {
        return &ghost_boxes[it->second];
    } else {
        return nullptr;
    }
    
}




// ===== ParallelTree Implementation =====

template<typename CoordType, typename DataType>
ParallelTree<CoordType, DataType>::ParallelTree()
    : dimension(0), num_levels(0), num_points(0),
      mpi_rank(0), mpi_size(1), comm(MPI_COMM_WORLD),
      reduction_pattern(ReductionPattern::UNIFORM), reduction_threshold(0),
      CA_level(10000), id_neighborhood_radius(2), id_proxy_mode(0),
      id_proxy_points(8), id_adaptive_batch(64) {
    for (int i = 0; i < 6; ++i) global_bounds[i] = 0.0;
}


/**
 * @brief Clears memory-intensive interaction matrices after factorization.
 * 
 * This function iterates through all local and ghost boxes in a given tree level
 * and deallocates the data associated with near-field and far-field modified 
 * interaction blocks (Ã matrices). 
 * 
 * This is typically called after the factorization phase is complete and before 
 * the solve phase begins, as these matrices are not required for the solve itself,
 * and clearing them significantly reduces the memory footprint.
 * 
 * @tparam CoordType The coordinate data type (e.g., double, float).
 * @tparam DataType The matrix element data type (e.g., double, std::complex<double>).
 * @param level The tree level whose interaction matrices will be cleared.
 */
template<typename CoordType, typename DataType>
void clear_modified_interaction_matrices(
    TreeLevel<CoordType, DataType>& level,
    bool preserve_ghost_boxes = false) {

    // Helper lambda to clear interaction matrices for a single box.
    // This avoids code duplication for local and ghost boxes.
    auto clear_box_interactions = [](BoxData<CoordType, DataType>& box) {
        
        // --- Clear Near-Field Interactions ---
        // The ModifiedBlock destructor will call the MatrixStorage destructor,
        // which in turn calls the std::vector destructor for its `data` member.
        // Using the swap idiom is a guaranteed way to release the vector's capacity.
        std::vector<ModifiedBlock<DataType>>().swap(box.near_field_modified_interactions);
        box.near_field_interaction_map.clear();

        // Also reset the associated counter to maintain a consistent state.
        box.num_near_field_interactions = 0;

        // --- Clear Far-Field Interactions ---
        std::vector<ModifiedBlock<DataType>>().swap(box.far_field_modified_interactions);
        box.far_field_interaction_map.clear();

        // Reset the counter.
        box.num_far_field_interactions = 0;
        box.X_RR_full = MatrixStorage<DataType>{};
        box.X_RS_entry = MatrixStorage<DataType>{};
        std::vector<int64_t>().swap(box.deferred_xnn_neighbor_point_counts);
        std::vector<DataType>().swap(box.lazy_original_x_nr);
    };

    // Iterate through all local boxes and clear their interaction data.
    for (auto& local_box : level.local_boxes) {
        clear_box_interactions(local_box);
    }

    // Iterate through all ghost boxes and clear their interaction data as well.
    for (auto& ghost_box : level.ghost_boxes) {
        clear_box_interactions(ghost_box);
    }

    if (!preserve_ghost_boxes) {
        std::vector<BoxData<CoordType, DataType>>().swap(level.ghost_boxes);
        level.ghost_id_to_index.clear();
    }
    std::vector<PointDataRequest<CoordType>>().swap(level.assisting_boxes);
    level.assisting_box_points_for_kernel_evaluation.clear();
    level.assisting_plan = {};
    std::vector<BoxData<CoordType, DataType>>().swap(level.generator_boxes);
    level.generator_id_to_index.clear();
    level.elimination_wave.clear();
    level.lazy_generator_requesters.clear();
    // malloc_trim(0);
}

template<typename CoordType, typename DataType>
void clear_ghosts(TreeLevel<CoordType, DataType>& level) {

    // Helper lambda to clear interaction matrices for a single box.
    // This avoids code duplication for local and ghost boxes.
    auto clear_box_interactions = [](BoxData<CoordType, DataType>& box) {
        
        // --- Clear Near-Field Interactions ---
        // The ModifiedBlock destructor will call the MatrixStorage destructor,
        // which in turn calls the std::vector destructor for its `data` member.
        // Using the swap idiom is a guaranteed way to release the vector's capacity.
        std::vector<ModifiedBlock<DataType>>().swap(box.near_field_modified_interactions);
        box.near_field_interaction_map.clear();

        // Also reset the associated counter to maintain a consistent state.
        box.num_near_field_interactions = 0;

        // --- Clear Far-Field Interactions ---
        std::vector<ModifiedBlock<DataType>>().swap(box.far_field_modified_interactions);
        box.far_field_interaction_map.clear();

        // Reset the counter.
        box.num_far_field_interactions = 0;
    };

    // Iterate through all ghost boxes and clear their interaction data as well.
    for (auto& ghost_box : level.ghost_boxes) {
        clear_box_interactions(ghost_box);
    }
    // malloc_trim(0);
}

/**
 * @brief Clear a PendingFactorUpdates container and aggressively release its memory.
 *
 * This clears both maps, swaps with empty maps to release bucket memory, and also frees
 * the internal DenseBlock data vectors by moving out of the maps first.
 *
 * Use this after each wave if you want to avoid carrying peak memory to the next wave.
 */
template <typename DataType>
void clear_pending_factor_updates_memory(PendingFactorUpdates<DataType>& p)
{
    std::unordered_map<int64_t, std::shared_ptr<GeneratorPayload<DataType>>>().swap(
        p.generators);
    // First free DenseBlock vectors (values) explicitly, then drop the hash tables.
    for (auto& kv : p.replace_blocks) {
        auto& b = kv.second;
        std::vector<DataType>().swap(b.data);
        b.rows = 0;
        b.cols = 0;
    }
    for (auto& kv : p.accumulated_deltas) {
        auto& b = kv.second;
        std::vector<DataType>().swap(b.data);
        b.rows = 0;
        b.cols = 0;
    }

    // Now release unordered_map bucket memory.
    std::unordered_map<ReplaceKey, DenseBlock<DataType>, ReplaceKeyHash>().swap(p.replace_blocks);
    std::unordered_map<EdgeKey, DenseBlock<DataType>, EdgeKeyHash>().swap(p.accumulated_deltas);
}

/**
 * @brief Clear a PendingSolveUpdates container and aggressively release its memory.
 *
 * Frees all stored update vectors and releases unordered_map bucket memory
 * by swapping with empty maps.
 */
template <typename DataType>
void clear_pending_solve_updates_memory(PendingSolveUpdates<DataType>& p)
{
    // Free vectors first (values).
    for (auto& kv : p.full_updates) {
        std::vector<DataType>().swap(kv.second);
    }
    for (auto& kv : p.skel_updates) {
        std::vector<DataType>().swap(kv.second);
    }

    // Release unordered_map bucket memory.
    std::unordered_map<int64_t, std::vector<DataType>>().swap(p.full_updates);
    std::unordered_map<int64_t, std::vector<DataType>>().swap(p.skel_updates);
}

// Destructor now handled by default - vectors clean up automatically

// ===== Helper Functions =====

/**
 * @brief Compute box center, size, and bounds from Morton index (FMM3D.pdf Section 4.4.3)
 * 
 * Ensures uniform proxy surface geometry by using MAXIMUM dimension range.
 * This creates square (2D) or cubic (3D) boxes, with some potentially extending
 * beyond the physical domain for rectangular domains.
 */
template<typename CoordType>
void compute_box_geometry(
    int64_t morton_index,
    int32_t level,
    const CoordType global_bounds[6],
    int32_t dimension,
    CoordType center[3],
    CoordType& size,
    CoordType bounds[6]) {
    
    // Decode Morton index to grid coordinates.
    uint32_t x = 0, y = 0, z = 0;
    morton::decode_nd(dimension, morton_index, x, y, z);
    
    // Calculate MAXIMUM dimension range (Section 4.4.3)
    // This ensures boxes are square/cubic for uniform proxy surfaces
    CoordType max_range = global_bounds[1] - global_bounds[0];
    for (int d = 1; d < dimension; ++d) {
        max_range = std::max(
            max_range, global_bounds[2 * d + 1] - global_bounds[2 * d]);
    }
    
    // Number of boxes per dimension at this level
    uint32_t boxes_per_dim = 1U << level;  // 2^level
    
    // Box size based on MAXIMUM dimension (ensures square/cubic boxes)
    CoordType box_size = max_range / boxes_per_dim;
    
    // Compute box bounds (positioned relative to domain minimum)
    const uint32_t grid_coordinates[3] = {x, y, z};
    for (int d = 0; d < 3; ++d) {
        if (d < dimension) {
            bounds[2 * d] = global_bounds[2 * d] +
                            grid_coordinates[d] * box_size;
            bounds[2 * d + 1] = bounds[2 * d] + box_size;
            center[d] = bounds[2 * d] + CoordType(0.5) * box_size;
        } else {
            bounds[2 * d] = CoordType(0);
            bounds[2 * d + 1] = CoordType(0);
            center[d] = CoordType(0);
        }
    }
    
    // Store size
    size = box_size;
}

/**
 * @brief Map point to grid coordinates and Morton index (FMM3D.pdf Section 4.4.4)
 */
template<typename CoordType>
int64_t point_to_morton(
    const CoordType* point,
    int32_t dimension,
    const CoordType global_bounds[6],
    int32_t level,
    int32_t grid_coords[3]) {
    
    // Calculate MAXIMUM dimension range
    CoordType max_range = global_bounds[1] - global_bounds[0];
    for (int d = 1; d < dimension; ++d) {
        max_range = std::max(
            max_range, global_bounds[2 * d + 1] - global_bounds[2 * d]);
    }
    
    // Number of boxes per dimension
    uint32_t boxes_per_dim = 1U << level;
    CoordType box_size = max_range / boxes_per_dim;
    
    uint32_t coordinates[3] = {0, 0, 0};
    for (int d = 0; d < dimension; ++d) {
        const auto raw = static_cast<int64_t>(
            std::floor((point[d] - global_bounds[2 * d]) / box_size));
        const auto clamped = std::max<int64_t>(
            0, std::min<int64_t>(raw, boxes_per_dim - 1));
        coordinates[d] = static_cast<uint32_t>(clamped);
    }

    for (int d = 0; d < 3; ++d) {
        grid_coords[d] = static_cast<int32_t>(coordinates[d]);
    }

    return morton::encode_nd(
        dimension, coordinates[0], coordinates[1], coordinates[2]);
}

/**
 * @brief The points in each box of the region weight level (region_partition.hpp),
 * for splitting the levels' ranges by points: every point binned at the leaf, as
 * the trees bin them, and counted in its box of that level.  `distributed`: each
 * rank holds its own points (the counts are summed over `comm`); else every rank
 * holds all of them and rank 0's count is broadcast.  Size: at most
 * kRegionWeightBoxes boxes a rank.
 */
template<typename CoordType>
RegionWeights region_weights_of_points(
    const CoordType* coords,
    int64_t num_points,
    int32_t dimension,
    int32_t num_levels,
    const CoordType global_bounds[6],
    MPI_Comm comm,
    bool distributed) {
    int rank = 0, size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const int level = region_weight_level(dimension, num_levels, size);
    const int leaf = num_levels - 1;
    const int shift = dimension * (leaf - level);
    std::vector<int64_t> count(size_t{1} << (dimension * level), 0);
    if (distributed || rank == 0) {
        if (num_points > 0 && coords == nullptr) {
            throw std::invalid_argument("region_weights_of_points: null coordinates");
        }
        for (int64_t i = 0; i < num_points; ++i) {
            CoordType point[3] = {CoordType(0), CoordType(0), CoordType(0)};
            for (int d = 0; d < dimension; ++d) point[d] = coords[i * dimension + d];
            int32_t grid[3];
            const int64_t morton = point_to_morton(point, dimension, global_bounds, leaf, grid);
            ++count[static_cast<size_t>(morton >> shift)];
        }
    }
    if (count.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("region_weights_of_points: too many boxes");
    }
    const int n = static_cast<int>(count.size());
    if (distributed) {
        MPI_Allreduce(MPI_IN_PLACE, count.data(), n, MPI_INT64_T, MPI_SUM, comm);
    } else {
        MPI_Bcast(count.data(), n, MPI_INT64_T, 0, comm);
    }
    return RegionWeights::of_counts(level, count);
}



/**
 * @brief Compute which processes are active at each level (FMM3D.pdf Section 5.3)
 * 
 * Also builds Morton ID mapping and communication patterns for process reduction.
 */
template<typename CoordType, typename DataType>
void compute_active_processes(ParallelTree<CoordType, DataType>* tree,
                              const RegionWeights& region_weights = RegionWeights{}) {
     /**
     * Process Reduction Logic (FMM3D.pdf Section 6.4):
     * 
     * When boxes_per_process < threshold, reduce active processes by factor of:
     * - 4 for 2D (2×2 grid of processes merge)
     * - 8 for 3D (2×2×2 cube of processes merge)
     * 
     * Active-process counts still follow the same reduction rule. When a
     * reduction happens, parent owners are selected in a node-aware way:
     *   1. Prefer currently inactive ranks, cycling across nodes.
     *   2. If every rank is still active (the first reduction), rotate within
     *      each child group instead of always taking the first child.
     * This lets coarse levels move across nodes without assigning a parent box
     * to a child rank from the wrong Morton group.
     */
    
    // Any rank count for the unstructured Color path (any_process_count):
    // the regions are then contiguous Morton ranges of unequal length
    // (region_partition.hpp); the structured and CA levels keep equal cubes.
    if (!tree->any_process_count && !morton::is_valid_process_count(tree->dimension, tree->mpi_size)) {
        throw std::invalid_argument(
            std::to_string(tree->dimension) +
            "D tree requires process count to be " +
            std::to_string(1 << tree->dimension) + "^k (any count with H2_unstructured=1)");
    }
    // the counts and ranges of every level (region_partition.hpp), the
    // ranges split by the points in the boxes where the caller gave them (the
    // unstructured path); possibly with the first reduction a level finer
    // than reduction_threshold asks, so that the leaf regions are even
    const RegionPlan region_plan = plan_regions(tree->dimension, tree->num_levels, tree->mpi_size,
                                                tree->reduction_threshold, region_weights);
    const std::vector<int32_t>& region_count = region_plan.count;
    const std::vector<std::vector<int64_t>>& region_start_of = region_plan.start;
    if (tree->mpi_rank == 0 && (!region_weights.empty() || !region_plan.bricks)) {
        std::cout << "H2 tree: " << tree->mpi_size << " ranks, "
                  << (region_plan.bricks ? "equal cubes" : "Morton ranges") << ": the largest leaf region holds "
                  << region_plan.leaf_imbalance << " of the mean " << (region_weights.empty() ? "boxes" : "points");
        if (region_plan.threshold != tree->reduction_threshold) {
            std::cout << "; reduction_threshold " << tree->reduction_threshold << " -> " << region_plan.threshold
                      << " (the first reduction a level finer, for whole parent boxes)";
        }
        std::cout << "; active ranks by level:";
        for (int32_t c : region_count) std::cout << " " << c;
        std::cout << std::endl;
        if (region_plan.empty_region) {
            std::cout << "H2 tree: warning: a leaf region has no points (no layout of these "
                      << tree->mpi_size << " ranks avoids it)" << std::endl;
        }
    }
    tree->reduction_threshold = region_plan.threshold;

    std::vector<int> gathered_node_roots(tree->mpi_size, -1);
    std::vector<int> gathered_local_ranks(tree->mpi_size, -1);
    MPI_Comm shared_comm = MPI_COMM_NULL;
    MPI_Comm_split_type(
        tree->comm, MPI_COMM_TYPE_SHARED, tree->mpi_rank, MPI_INFO_NULL, &shared_comm);

    int local_rank_on_node = 0;
    MPI_Comm_rank(shared_comm, &local_rank_on_node);

    int node_root_rank = tree->mpi_rank;
    MPI_Allreduce(
        &tree->mpi_rank, &node_root_rank, 1, MPI_INT, MPI_MIN, shared_comm);

    MPI_Allgather(
        &node_root_rank, 1, MPI_INT,
        gathered_node_roots.data(), 1, MPI_INT, tree->comm);
    MPI_Allgather(
        &local_rank_on_node, 1, MPI_INT,
        gathered_local_ranks.data(), 1, MPI_INT, tree->comm);

    if (shared_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&shared_comm);
    }

    std::vector<int> unique_node_roots = gathered_node_roots;
    std::sort(unique_node_roots.begin(), unique_node_roots.end());
    unique_node_roots.erase(
        std::unique(unique_node_roots.begin(), unique_node_roots.end()),
        unique_node_roots.end());

    std::vector<std::vector<int>> ranks_by_node(unique_node_roots.size());
    int max_ranks_per_node = 0;
    for (size_t node_idx = 0; node_idx < unique_node_roots.size(); ++node_idx) {
        for (int rank = 0; rank < tree->mpi_size; ++rank) {
            if (gathered_node_roots[rank] == unique_node_roots[node_idx]) {
                ranks_by_node[node_idx].push_back(rank);
            }
        }
        std::sort(
            ranks_by_node[node_idx].begin(),
            ranks_by_node[node_idx].end(),
            [&](int lhs, int rhs) {
                if (gathered_local_ranks[lhs] != gathered_local_ranks[rhs]) {
                    return gathered_local_ranks[lhs] < gathered_local_ranks[rhs];
                }
                return lhs < rhs;
            });
        max_ranks_per_node = std::max(
            max_ranks_per_node,
            static_cast<int>(ranks_by_node[node_idx].size()));
    }

    std::vector<int> active_rank_cycle;
    active_rank_cycle.reserve(static_cast<size_t>(tree->mpi_size));
    for (int local_slot = 0; local_slot < max_ranks_per_node; ++local_slot) {
        for (const auto& node_ranks : ranks_by_node) {
            if (local_slot < static_cast<int>(node_ranks.size())) {
                active_rank_cycle.push_back(node_ranks[static_cast<size_t>(local_slot)]);
            }
        }
    }

    if (active_rank_cycle.size() != static_cast<size_t>(tree->mpi_size)) {
        throw std::runtime_error("compute_active_processes: invalid active-rank cycle size");
    }

    std::vector<int> active_rank_cycle_position(tree->mpi_size, -1);
    for (size_t idx = 0; idx < active_rank_cycle.size(); ++idx) {
        active_rank_cycle_position[active_rank_cycle[idx]] = static_cast<int>(idx);
    }

    size_t inactive_rank_cycle_offset = 0;
    int reduction_event = 0;

    int32_t leaf_level = tree->num_levels - 1;
    tree->levels[leaf_level].num_active_processes = tree->mpi_size;
    tree->levels[leaf_level].active_process_ranks.resize(tree->mpi_size);
    for (int32_t i = 0; i < tree->mpi_size; ++i) {
        tree->levels[leaf_level].active_process_ranks[i] = i;
    }
    tree->levels[leaf_level].is_process_active = (tree->mpi_rank < tree->mpi_size);
    
    for (int32_t i = 0; i < tree->mpi_size; ++i) {
        tree->levels[leaf_level].morton_to_rank[i] = i;
        tree->levels[leaf_level].rank_to_morton[i] = i;
    }
    tree->levels[leaf_level].my_morton_id = tree->mpi_rank;
    tree->levels[leaf_level].parent_level_owner = -1;
    tree->levels[leaf_level].region_start = region_start_of[static_cast<size_t>(leaf_level)];
    tree->levels[leaf_level].regions_cubic =
        region_starts_are_bricks(tree->dimension, leaf_level, region_start_of[static_cast<size_t>(leaf_level)]);
    
    for (int32_t level = leaf_level - 1; level >= 0; --level) {
        auto& lvl = tree->levels[level];
        auto& child_lvl = tree->levels[level + 1];
        lvl.active_process_ranks.clear();
        lvl.morton_to_rank.clear();
        lvl.rank_to_morton.clear();
        lvl.children_senders.clear();
        
        int32_t prev_active = child_lvl.num_active_processes;

        switch (tree->reduction_pattern) {
            case ReductionPattern::UNIFORM:
                // fewer than reduction_threshold boxes per process of the level
                // below: prev_active / 2^d, rounded up (region_counts; the same
                // as before for (2^d)^k)
                lvl.num_active_processes = region_count[static_cast<size_t>(level)];
                break;
            case ReductionPattern::ADAPTIVE:
            case ReductionPattern::WORK_PRESERVING:
            case ReductionPattern::AGGRESSIVE:
                throw std::runtime_error("Reduction pattern not yet implemented");
        }

        lvl.active_process_ranks.resize(lvl.num_active_processes);
        
        if (lvl.num_active_processes < prev_active) {
            if (lvl.num_active_processes > tree->mpi_size) {
                throw std::runtime_error("compute_active_processes: active process count exceeds MPI size");
            }

            std::vector<int> inactive_ranks;
            inactive_ranks.reserve(static_cast<size_t>(tree->mpi_size - prev_active));
            for (int rank : active_rank_cycle) {
                if (child_lvl.rank_to_morton.find(rank) == child_lvl.rank_to_morton.end()) {
                    inactive_ranks.push_back(rank);
                }
            }

            const bool use_inactive_owners =
                inactive_ranks.size() >= static_cast<size_t>(lvl.num_active_processes);
            for (int32_t i = 0; i < lvl.num_active_processes; ++i) {
                int parent_rank = -1;
                if (use_inactive_owners) {
                    const size_t cycle_idx =
                        (inactive_rank_cycle_offset + static_cast<size_t>(i)) % inactive_ranks.size();
                    parent_rank = inactive_ranks[cycle_idx];
                } else {
                    // the child regions under region i (region_first_child:
                    // groups of 2^d for (2^d)^k, else as even as can be)
                    std::vector<int> group_ranks;
                    const int child_morton_start = region_first_child(i, lvl.num_active_processes, prev_active);
                    const int child_morton_end = region_first_child(i + 1, lvl.num_active_processes, prev_active);
                    group_ranks.reserve(static_cast<size_t>(child_morton_end - child_morton_start));
                    for (int child_morton = child_morton_start;
                         child_morton < child_morton_end;
                         ++child_morton) {
                        auto child_it = child_lvl.morton_to_rank.find(child_morton);
                        if (child_it == child_lvl.morton_to_rank.end()) {
                            throw std::runtime_error(
                                "compute_active_processes: missing child Morton owner for reduction group");
                        }
                        group_ranks.push_back(child_it->second);
                    }
                    std::sort(
                        group_ranks.begin(),
                        group_ranks.end(),
                        [&](int lhs, int rhs) {
                            return active_rank_cycle_position[lhs] < active_rank_cycle_position[rhs];
                        });

                    const size_t group_offset =
                        static_cast<size_t>(reduction_event + i) % group_ranks.size();
                    parent_rank = group_ranks[group_offset];
                }
                lvl.active_process_ranks[i] = parent_rank;
                lvl.morton_to_rank[i] = parent_rank;
                lvl.rank_to_morton[parent_rank] = i;
            }
            if (use_inactive_owners) {
                inactive_rank_cycle_offset =
                    (inactive_rank_cycle_offset + static_cast<size_t>(lvl.num_active_processes)) %
                    inactive_ranks.size();
            }
            ++reduction_event;
        } else {
            for (int32_t i = 0; i < lvl.num_active_processes; ++i) {
                lvl.active_process_ranks[i] = child_lvl.active_process_ranks[i];
                lvl.morton_to_rank[i] = child_lvl.morton_to_rank[i];
                int rank = child_lvl.morton_to_rank[i];
                lvl.rank_to_morton[rank] = i;
            }
        }
        
        lvl.region_start = region_start_of[static_cast<size_t>(level)];
        lvl.regions_cubic = region_starts_are_bricks(tree->dimension, level, region_start_of[static_cast<size_t>(level)]);

        auto rank_it = lvl.rank_to_morton.find(tree->mpi_rank);
        if (rank_it != lvl.rank_to_morton.end()) {
            lvl.is_process_active = true;
            lvl.my_morton_id = rank_it->second;
        } else {
            lvl.is_process_active = false;
            lvl.my_morton_id = -1;
        }
        
        if (child_lvl.is_process_active) {
            if (lvl.num_active_processes < child_lvl.num_active_processes) {
                int my_child_morton = child_lvl.my_morton_id;
                int parent_morton =
                    region_parent(my_child_morton, lvl.num_active_processes, child_lvl.num_active_processes);
                
                auto parent_it = lvl.morton_to_rank.find(parent_morton);
                if (parent_it != lvl.morton_to_rank.end()) {
                    child_lvl.parent_level_owner = parent_it->second;
                } else {
                    throw std::runtime_error(
                        "compute_active_processes: Parent morton " + std::to_string(parent_morton) +
                        " not found in level " + std::to_string(level));
                }
            } else {
                child_lvl.parent_level_owner = tree->mpi_rank;
            }
        }
        
        if (lvl.is_process_active) {
            lvl.children_senders.clear();
            
            if (lvl.num_active_processes < child_lvl.num_active_processes) {
                // printf("mpi rank: %d, Level %d: Active Processes = %d, This Rank Active = %s\n", tree->mpi_rank,
                //     level, lvl.num_active_processes,
                //     lvl.is_process_active ? "Yes" : "No");
                int parent_morton = lvl.my_morton_id;
                int child_morton_start =
                    region_first_child(parent_morton, lvl.num_active_processes, child_lvl.num_active_processes);
                int child_morton_end =
                    region_first_child(parent_morton + 1, lvl.num_active_processes, child_lvl.num_active_processes) - 1;
                // printf("child morton start: %d, child morton end: %d\n", child_morton_start, child_morton_end);
                
                for (int child_morton = child_morton_start; child_morton <= child_morton_end; ++child_morton) {
                    int child_rank = child_lvl.morton_to_rank[child_morton];
                    lvl.children_senders.push_back(child_rank);
                }
            }
            // printf("mpi rank: %d, Level %d: children size: %zu\n", tree->mpi_rank, level, lvl.children_senders.size());
        }
    }
    fflush(stdout);
}

/**
 * @brief Initialize local boxes with geometry and hierarchy information
 */
template<typename CoordType, typename DataType>
void initialize_local_boxes(
    std::vector<BoxData<CoordType, DataType>>& local_boxes,
    int64_t local_morton_start,
    int32_t level,
    int dimension,
    const CoordType global_bounds[6]) {
    
    int64_t num_boxes_local = local_boxes.size();
    
    for (int64_t i = 0; i < num_boxes_local; ++i) {
        int64_t morton_idx = local_morton_start + i;
        auto& box = local_boxes[i];
        
        box.morton_index = morton_idx;
        box.level = level;
        
        uint32_t x = 0, y = 0, z = 0;
        morton::decode_nd(dimension, morton_idx, x, y, z);
        box.grid_coords[0] = static_cast<int32_t>(x);
        box.grid_coords[1] = static_cast<int32_t>(y);
        box.grid_coords[2] = static_cast<int32_t>(z);
        
        // Compute geometry using helper function
        compute_box_geometry(morton_idx, level, global_bounds,
                           dimension, box.center, box.size, box.bounds);
        
        // Set parent (if not root)
        if (level > 0) {
            box.parent_morton =
                morton_idx / morton::children_per_box(dimension);
        } else {
            box.parent_morton = -1;
        }
        
        // Initialize children to -1
        box.num_children = 0;
        for (int j = 0; j < 8; ++j) {
            box.children_morton[j] = -1;
        }
    }
}

/**
 * @brief Distribute boxes to processes using Morton ordering (FMM3D.pdf Section 5.1.3)
 */
template<typename CoordType, typename DataType>
void distribute_boxes(
    int32_t level,
    ParallelTree<CoordType, DataType>* tree) {
    
    auto& lvl = tree->levels[level];
    lvl.level = level;
    int64_t total_boxes = 1LL << (tree->dimension * level);
    lvl.num_boxes_global = total_boxes;
    
    // Check if this process is active at this level using rank_to_morton
    auto rank_it = lvl.rank_to_morton.find(tree->mpi_rank);
    if (rank_it == lvl.rank_to_morton.end()) {
        // Process is inactive at this level
        lvl.num_boxes_local = 0;
        lvl.local_morton_start = -1;
        lvl.local_morton_end = -1;
        lvl.local_boxes.clear();
        return;
    }
    
    // Get this process's Morton region ID (O(1) lookup)
    int32_t active_index = rank_it->second;  // This is my_morton_id
    
    // The region's Morton range (region_partition.hpp: the equal ranges of
    // Section 5.1.3 for (2^d)^k regions, else unequal ones)
    if (lvl.region_start.size() != static_cast<size_t>(lvl.num_active_processes) + 1 ||
        lvl.region_start.back() != total_boxes) {
        throw std::runtime_error("distribute_boxes: level " + std::to_string(level) + " has no region starts");
    }
    lvl.local_morton_start = lvl.region_start[static_cast<size_t>(active_index)];
    lvl.local_morton_end = lvl.region_start[static_cast<size_t>(active_index) + 1] - 1;
    lvl.num_boxes_local = lvl.local_morton_end - lvl.local_morton_start + 1;

    // Create local boxes and initialize them
    lvl.local_boxes.resize(lvl.num_boxes_local);
    initialize_local_boxes(lvl.local_boxes, lvl.local_morton_start, 
                          level, tree->dimension, tree->global_bounds);

    lvl.solve_neighbor_size.resize(lvl.num_boxes_local);
    for(auto& inner_list: lvl.solve_neighbor_size){
        inner_list.clear();
    }

}


/**
 * @brief Broadcast points from root to all processes and assign to boxes
 * 
 * @param point_coords Point coordinates in column-major format (only valid on root)
 * @param num_points Number of points
 * @param tree Parallel tree structure
 */
template<typename CoordType, typename DataType>
void assign_points_to_boxes(
    const CoordType* point_coords,
    int64_t num_points,
    std::vector<int>& idx_map,
    ParallelTree<CoordType, DataType>* tree) {
    
    if (num_points == 0) {
        return;
    }
    
    int rank;
    MPI_Comm_rank(tree->comm, &rank);
    
    // Broadcast points to all processes
    std::vector<CoordType> all_points;
    if (rank == 0) {
        if (point_coords == nullptr) {
            throw std::invalid_argument("Point coordinates null on root process");
        }
        all_points.assign(point_coords, point_coords + tree->dimension * num_points);
    } else {
        all_points.resize(tree->dimension * num_points);
    }
    
    MPI_Datatype mpi_type = std::is_same_v<CoordType, float> ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Bcast(all_points.data(), tree->dimension * num_points, mpi_type, 0, tree->comm);
    
    int32_t leaf_level = tree->num_levels - 1;
    auto& leaf_lvl = tree->levels[leaf_level];
    
    if (!leaf_lvl.is_process_active || leaf_lvl.num_boxes_local == 0) {
        return;
    }
    
    // Build map of Morton index to local box index
    std::unordered_map<int64_t, int64_t> morton_to_local;
    for (int64_t i = 0; i < leaf_lvl.num_boxes_local; ++i) {
        morton_to_local[leaf_lvl.local_boxes[i].morton_index] = i;
    }
    
    // Temporary vectors to collect points per box
    std::vector<std::vector<int64_t>> box_point_indices(leaf_lvl.num_boxes_local);
    std::vector<std::vector<CoordType>> box_point_coords(leaf_lvl.num_boxes_local);
    
    // Assign each point to its box
    for (int64_t pt_idx = 0; pt_idx < num_points; ++pt_idx) {
        CoordType point[3] = {CoordType(0), CoordType(0), CoordType(0)};
        for (int d = 0; d < tree->dimension; ++d) {
            point[d] = all_points[pt_idx * tree->dimension + d];
        }
        
        // Compute Morton index for this point
        int32_t grid_coords[3];
        int64_t morton_idx = point_to_morton(point, tree->dimension,
                                            tree->global_bounds, leaf_level, grid_coords);
        
        // Check if this box is local
        auto it = morton_to_local.find(morton_idx);
        if (it != morton_to_local.end()) {
            int64_t local_idx = it->second;
            box_point_indices[local_idx].push_back(pt_idx);
            
            // Store coordinates in column-major format
            for (int d = 0; d < tree->dimension; ++d) {
                box_point_coords[local_idx].push_back(point[d]);
            }
        }
    }
    
    // Copy data into boxes
    for (int64_t i = 0; i < leaf_lvl.num_boxes_local; ++i) {
        auto& box = leaf_lvl.local_boxes[i];
        int64_t n = box_point_indices[i].size();
        
        box.num_points = n;
        if (n > 0) {
            box.point_indices = box_point_indices[i];
            box.point_coords = box_point_coords[i];
            for (int global_pt_idx : box.point_indices) {
                idx_map.push_back(global_pt_idx);
            }
        }
    }

}

/**
 * @brief Assign uniformly distributed points to boxes (no MPI needed)
 * 
 * Each process independently generates the uniform grid points that fall
 * within its local boxes. No communication required.
 * 
 * @param N Total number of points (will create n_per_dim^dimension grid)
 *          For 2D: sqrt(N) × sqrt(N) grid
 *          For 3D: cbrt(N) × cbrt(N) × cbrt(N) grid
 * @param tree Parallel tree structure
 */
template<typename CoordType, typename DataType>
void assign_uniform_points(
    int64_t N,
    ParallelTree<CoordType, DataType>* tree) {
    
    int32_t leaf_level = tree->num_levels - 1;
    auto& leaf_lvl = tree->levels[leaf_level];
    
    if (!leaf_lvl.is_process_active) {
        return;
    }
    
    // Compute grid dimension.
    int n_per_dim;
    if (tree->dimension == 1) {
        n_per_dim = static_cast<int>(N);
    } else if (tree->dimension == 2) {
        n_per_dim = static_cast<int>(std::sqrt(static_cast<double>(N)));
        // Ensure N = n_per_dim^2
        if (n_per_dim * n_per_dim != N) {
            n_per_dim = static_cast<int>(std::sqrt(static_cast<double>(N))) + 1;
        }
    } else {
        n_per_dim = static_cast<int>(std::cbrt(static_cast<double>(N)));
        // Ensure N = n_per_dim^3
        if (n_per_dim * n_per_dim * n_per_dim != N) {
            n_per_dim = static_cast<int>(std::cbrt(static_cast<double>(N))) + 1;
        }
    }
    
    // Grid spacing
    CoordType hx = (tree->global_bounds[1] - tree->global_bounds[0]) / n_per_dim;
    CoordType hy = (tree->dimension >= 2)
        ? (tree->global_bounds[3] - tree->global_bounds[2]) / n_per_dim
        : CoordType(0);
    CoordType hz = (tree->dimension == 3) ? 
                   (tree->global_bounds[5] - tree->global_bounds[4]) / n_per_dim : 0;
    
    // For each local box, generate points that fall inside
    for (int64_t box_idx = 0; box_idx < leaf_lvl.num_boxes_local; ++box_idx) {
        auto& box = leaf_lvl.local_boxes[box_idx];
        
        // Box bounds
        CoordType x_min = box.bounds[0];
        CoordType x_max = box.bounds[1];
        // Find grid indices that overlap with this box
        int i_min = static_cast<int>(std::floor((x_min - tree->global_bounds[0]) / hx));
        int i_max = static_cast<int>(std::floor((x_max - tree->global_bounds[0]) / hx));
        
        // Clamp to valid range
        i_min = std::max(0, std::min(n_per_dim - 1, i_min));
        i_max = std::max(0, std::min(n_per_dim - 1, i_max));

        if (tree->dimension == 1) {
            for (int i = i_min; i <= i_max; ++i) {
                CoordType x = tree->global_bounds[0] + (i + 0.5) * hx;
                if (x >= x_min && x < x_max) {
                    box.point_indices.push_back(i);
                    box.point_coords.push_back(x);
                    box.num_points++;
                }
            }
            continue;
        }

        CoordType y_min = box.bounds[2];
        CoordType y_max = box.bounds[3];
        int j_min = static_cast<int>(std::floor((y_min - tree->global_bounds[2]) / hy));
        int j_max = static_cast<int>(std::floor((y_max - tree->global_bounds[2]) / hy));
        j_min = std::max(0, std::min(n_per_dim - 1, j_min));
        j_max = std::max(0, std::min(n_per_dim - 1, j_max));
        
        if (tree->dimension == 2) {
            // 2D case
            for (int i = i_min; i <= i_max; ++i) {
                for (int j = j_min; j <= j_max; ++j) {
                    // Point at grid cell center
                    CoordType x = tree->global_bounds[0] + (i + 0.5) * hx;
                    CoordType y = tree->global_bounds[2] + (j + 0.5) * hy;
                    
                    // Check if point is in box (use < for upper bound to avoid duplicates)
                    if (x >= x_min && x < x_max && y >= y_min && y < y_max) {
                        int64_t pt_idx = i * n_per_dim + j;
                        box.point_indices.push_back(pt_idx);
                        box.point_coords.push_back(x);
                        box.point_coords.push_back(y);
                        box.num_points++;
                    }
                }
            }
        } else {
            // 3D case
            CoordType z_min = box.bounds[4];
            CoordType z_max = box.bounds[5];
            
            int k_min = static_cast<int>(std::floor((z_min - tree->global_bounds[4]) / hz));
            int k_max = static_cast<int>(std::floor((z_max - tree->global_bounds[4]) / hz));
            k_min = std::max(0, std::min(n_per_dim - 1, k_min));
            k_max = std::max(0, std::min(n_per_dim - 1, k_max));
            
            for (int i = i_min; i <= i_max; ++i) {
                for (int j = j_min; j <= j_max; ++j) {
                    for (int k = k_min; k <= k_max; ++k) {
                        CoordType x = tree->global_bounds[0] + (i + 0.5) * hx;
                        CoordType y = tree->global_bounds[2] + (j + 0.5) * hy;
                        CoordType z = tree->global_bounds[4] + (k + 0.5) * hz;
                        
                        if (x >= x_min && x < x_max && 
                            y >= y_min && y < y_max && 
                            z >= z_min && z < z_max) {
                            int64_t pt_idx = (i * n_per_dim + j) * n_per_dim + k;
                            box.point_indices.push_back(pt_idx);
                            box.point_coords.push_back(x);
                            box.point_coords.push_back(y);
                            box.point_coords.push_back(z);
                            box.num_points++;
                        }
                    }
                }
            }
        }
    }
}




/**
 * @brief Check if a box is on the boundary
 * 
 * A box is on the boundary if:
 *   - It has at least one 1-hop neighbor on a different process, OR
 *   - It is on the global domain boundary (coordinates at 0 or grid_size-1)
 * 
 * @return true if box is on boundary, false if interior
 */
bool check_boundary_condition(
    uint64_t morton_index,
    const auto& grid_coords,  // Box grid coordinates
    uint32_t grid_size,
    int dimension,
    int64_t local_morton_start,
    int64_t local_morton_end) {
    
    // Check criterion 1: Has neighbors on different processes
    const std::vector<uint64_t> neighbors =
        morton::neighbors_nd(dimension, morton_index, grid_size);
    
    for (uint64_t neighbor_morton_u : neighbors) {
        int64_t neighbor_morton = static_cast<int64_t>(neighbor_morton_u);
        
        // If neighbor is outside our local range, it's on a different process
        if (neighbor_morton < local_morton_start || neighbor_morton > local_morton_end) {
            return true;
        }
    }
    
    // Check criterion 2: On global domain boundary
    for (int d = 0; d < dimension; ++d) {
        if (grid_coords[d] == 0 || 
            grid_coords[d] == static_cast<int32_t>(grid_size) - 1) {
            return true;
        }
    }
    
    return false;
}

/**
 * @brief Populate interior and boundary box classification
 * 
 * Boundary box: 
 *   - Has at least one 1-hop neighbor on a different process, OR
 *   - Is on the global domain boundary (coordinates at 0 or grid_size-1)
 * Interior box: All other boxes
 * 
 * This ensures boxes near global corners/edges are always classified as boundary,
 * which is necessary for proper color classification in the FMM algorithm.
 */
template<typename CoordType, typename DataType>
void populate_interior_boundary(
    TreeLevel<CoordType, DataType>& lvl,
    ParallelTree<CoordType, DataType>* tree) {
    
    lvl.interior_id.clear();
    lvl.boundary_id.clear();
    lvl.eliminated_boxes.clear();
    
    uint32_t grid_size = 1 << lvl.level;
    
    // Loop through all local boxes
    for (int64_t local_idx = 0; local_idx < lvl.num_boxes_local; ++local_idx) {
        auto& box = lvl.local_boxes[local_idx];
        
        // Check if box is on boundary
        bool is_boundary = check_boundary_condition(
            box.morton_index,
            box.grid_coords,
            grid_size,
            tree->dimension,
            lvl.local_morton_start,
            lvl.local_morton_end
        );
        
        // Classify box
        if (is_boundary) {
            lvl.boundary_id.push_back(local_idx);
            lvl.local_boxes[local_idx].on_boundary = true;
        } else {
            lvl.interior_id.push_back(local_idx);
            lvl.local_boxes[local_idx].on_boundary = false;
        }
    }
}

/**
 * @brief Compute and store 1-hop and 2-hop neighbor lists
 * 
 * 1-hop: All boxes in 3×3 (2D) or 3×3×3 (3D) region, excluding self
 * 2-hop: All boxes in 5×5 (2D) or 5×5×5 (3D) region, excluding self and 1-hop
 */
template<typename CoordType, typename DataType>
void compute_neighbor_lists(
    std::vector<BoxData<CoordType, DataType>> &local_boxes,
    int dimension, int level) {
    
    uint32_t grid_size = 1 << level;
    
    for (int64_t i = 0; i < local_boxes.size(); ++i) {
        auto& box = local_boxes[i];
        
        // ===== Compute 1-hop neighbors =====
        const std::vector<uint64_t> neighbors_1hop =
            morton::neighbors_nd(dimension, box.morton_index, grid_size);
        
        // Convert and store
        box.one_hop.clear();
        box.one_hop.reserve(neighbors_1hop.size());
        for (uint64_t n : neighbors_1hop) {
            box.one_hop.push_back(static_cast<int64_t>(n));
        }
        box.use_full_set.clear();
        box.use_full_set.resize(box.one_hop.size(), 1);
        
        // ===== Compute 2-hop neighbors =====
        const std::vector<uint64_t> neighbors_2hop =
            morton::neighbors_2hop_nd(dimension, box.morton_index, grid_size);
        
        // Convert and store
        box.two_hop.clear();
        box.two_hop.reserve(neighbors_2hop.size());
        for (uint64_t n : neighbors_2hop) {
            box.two_hop.push_back(static_cast<int64_t>(n));
        }
    }
}




/**
 * @brief Compute assisting boxes for kernel evaluation and solve phase
 * 
 * Assisting boxes are remote boxes needed for accurate kernel evaluation.
 * For each ghost box, we compute its 1-hop and 2-hop neighbors that
 * are not on this process and not already in the ghost set.
 * 
 * The map stores: morton_index -> index in future points vector
 * This index will be used when receiving points via MPI.
 * 
 * @param level Tree level to process
 * @param tree Parallel tree structure
 */
template<typename CoordType, typename DataType>
void compute_assisting_boxes(
    int32_t level,
    ParallelTree<CoordType, DataType>* tree) {
    
    auto& lvl = tree->levels[level];
    
    if (!lvl.is_process_active || lvl.num_boxes_local == 0) {
        return;
    }
    
    const bool use_CA = tree->level_uses_CA(level);
    // the seeds: the ghosts of a CA level, else the boundary boxes
    std::vector<int64_t> seeds;
    if (use_CA) {
        seeds = lvl.ghost_id;
    } else {
        for (int64_t local_idx : lvl.boundary_id) seeds.push_back(lvl.local_boxes[local_idx].morton_index);
        std::sort(seeds.begin(), seeds.end());
    }
    RegionAssisting ra;
    if (!lvl.regions_cubic) {
        // Regions of unequal length (any rank count; the unstructured Color
        // path, never CA): the boxes of other ranks within one and two hops
        // of the boundary boxes, by the level's ranges - all the boxes of
        // other ranks within two hops of any box here, since a path out of
        // this rank's range leaves it at a boundary box.
        if (use_CA) {
            throw std::runtime_error("compute_assisting_boxes: level " + std::to_string(level) + " is a CA level but its " +
                                     std::to_string(lvl.num_active_processes) + " regions are not equal cubes");
        }
        const uint32_t grid = 1u << level;
        std::set<int64_t> all, one;
        auto other = [&](int64_t n) { return n < lvl.local_morton_start || n > lvl.local_morton_end; };
        for (int64_t seed : seeds) {
            for (uint64_t nu : morton::neighbors_nd(tree->dimension, static_cast<uint64_t>(seed), grid)) {
                const int64_t n = static_cast<int64_t>(nu);
                if (other(n)) {
                    all.insert(n);
                    one.insert(n);
                }
            }
            for (uint64_t nu : morton::neighbors_2hop_nd(tree->dimension, static_cast<uint64_t>(seed), grid)) {
                const int64_t n = static_cast<int64_t>(nu);
                if (other(n)) all.insert(n);
            }
        }
        ra.all.assign(all.begin(), all.end());
        ra.one_hop.assign(one.begin(), one.end());
        // The ranks this one exchanges with (compute_one_hop_neighbor_ranks):
        // the owners of the boxes of other ranks within three hops of the
        // boundary boxes, i.e. of any box here (the same argument) - two hops
        // for the factorization's assisting boxes, three for the interaction
        // lists of the H2 construction and the matvec (the children of the
        // parent's neighbors).  Symmetric, as the pairwise exchanges need.
        constexpr int kPeerHops = 3;
        std::set<int> peers;
        const int dim = tree->dimension;
        for (int64_t seed : seeds) {
            uint32_t x = 0, y = 0, z = 0;
            morton::decode_nd(dim, static_cast<uint64_t>(seed), x, y, z);
            const int zr = dim >= 3 ? kPeerHops : 0, yr = dim >= 2 ? kPeerHops : 0;
            for (int dz = -zr; dz <= zr; ++dz) {
                for (int dy = -yr; dy <= yr; ++dy) {
                    for (int dx = -kPeerHops; dx <= kPeerHops; ++dx) {
                        const int64_t nx = static_cast<int64_t>(x) + dx, ny = static_cast<int64_t>(y) + dy,
                                      nz = static_cast<int64_t>(z) + dz;
                        if (nx < 0 || ny < 0 || nz < 0 || nx >= grid || ny >= grid || nz >= grid) continue;
                        const int64_t n = static_cast<int64_t>(morton::encode_nd(
                            dim, static_cast<uint32_t>(nx), static_cast<uint32_t>(ny), static_cast<uint32_t>(nz)));
                        if (other(n)) peers.insert(lvl.owner_rank_of(n));
                    }
                }
            }
        }
        peers.erase(tree->mpi_rank);
        lvl.region_peers.assign(peers.begin(), peers.end());
        lvl.ca_halo_sends.clear();
        lvl.ca_halo_sends_valid = false;
    } else {
        const RegionLayout L = RegionLayout::of(tree->dimension, level, lvl.num_active_processes);
        const uint32_t my_process_id = static_cast<uint32_t>(lvl.my_morton_id);
        if (!L.uniform || lvl.local_morton_start != L.first(my_process_id) ||
            lvl.num_boxes_local != L.boxes_per_region) {
            throw std::runtime_error("compute_assisting_boxes: level " + std::to_string(level) +
                                     " has no uniform region layout (" + std::to_string(lvl.num_active_processes) +
                                     " active processes, " + std::to_string(lvl.num_boxes_local) + " local boxes)");
        }
        const std::vector<int64_t> no_ghosts;
        ra = assisting_of_region_reference(L, my_process_id, seeds, use_CA ? lvl.ghost_id : no_ghosts);
        if (use_CA) {
            // the set by region of the boxes, as the level-entry halo gather's
            // senders compute it for this region (halo_geometry.hpp), must agree
            std::vector<int64_t> by_region;
            for_each_region_within(L, my_process_id, halo_region_radius(L), [&](uint32_t r) {
                const std::vector<int64_t> part = assisting_of_region_in(L, my_process_id, r);
                by_region.insert(by_region.end(), part.begin(), part.end());
            });
            std::sort(by_region.begin(), by_region.end());
            if (by_region != ra.all) {
                throw std::runtime_error("compute_assisting_boxes: the assisting sets by region and by enumeration "
                                         "differ (level " + std::to_string(level) + ")");
            }
            // what the regions around this one need from it (the level-entry
            // halo gather's sends): here, not in the gather, where a rank with
            // 26 neighboring regions would keep those with 7 waiting
            std::vector<uint32_t> regions;
            for_each_region_within(L, my_process_id, halo_region_radius(L), [&](uint32_t q) { regions.push_back(q); });
            using Send = typename TreeLevel<CoordType, DataType>::CAHaloSend;
            std::vector<Send> sends(regions.size());
            #pragma omp parallel for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(regions.size()); ++i) {
                const uint32_t q = regions[static_cast<size_t>(i)];
                const RegionExtra e = region_extra(L, q);
                Send& send = sends[static_cast<size_t>(i)];
                send.region = q;
                send.ghosts = ghosts_of_region_in(L, q, my_process_id, e);
                send.assisting = assisting_of_region_in(L, q, my_process_id, e);
            }
            lvl.ca_halo_sends.clear();
            for (Send& send : sends) {
                if (!send.ghosts.empty() || !send.assisting.empty()) lvl.ca_halo_sends.push_back(std::move(send));
            }
            lvl.ca_halo_sends_valid = true;
        } else {
            lvl.ca_halo_sends.clear();
            lvl.ca_halo_sends_valid = false;
        }
    }
    
    lvl.assisting_box_points_for_kernel_evaluation.clear();
    for (size_t i = 0; i < ra.all.size(); ++i) {
        lvl.assisting_box_points_for_kernel_evaluation[ra.all[i]] = i;
    }

    // the set used during the solve phase: the one-hop part, and the ghosts of a CA level
    std::set<int64_t> assisting_and_ghost_set_for_solve(ra.one_hop.begin(), ra.one_hop.end());
    std::set<int64_t> solve_assist_set;
    if (use_CA) {
        solve_assist_set.insert(ra.one_hop.begin(), ra.one_hop.end());
        assisting_and_ghost_set_for_solve.insert(lvl.ghost_id.begin(), lvl.ghost_id.end());
    }
    const std::vector<int64_t> assisting_and_ghost_vec_for_solve(assisting_and_ghost_set_for_solve.begin(),
                                                                 assisting_and_ghost_set_for_solve.end());
    lvl.ghost_and_assisting_box_points_for_solve_map.clear();
    lvl.is_ghost_solve.clear();
    for (size_t i = 0; i < assisting_and_ghost_vec_for_solve.size(); ++i) {
        lvl.ghost_and_assisting_box_points_for_solve_map[assisting_and_ghost_vec_for_solve[i]] = i;
        lvl.is_ghost_solve.push_back(
            use_CA && solve_assist_set.count(assisting_and_ghost_vec_for_solve[i]) == 0);
    }
    lvl.assisting_boxes.resize(lvl.assisting_box_points_for_kernel_evaluation.size());

}

template<typename CoordType, typename DataType>
void compute_ghost_boxes_and_colors(
    int32_t level,
    ParallelTree<CoordType, DataType>* tree);

template<typename CoordType, typename DataType>
ParallelTree<CoordType, DataType>* create_uniform_tree(
    const CoordType* point_coords,
    int64_t num_points,
    int32_t num_levels,
    const CoordType global_bounds[6],
    int32_t dimension = 3,
    MPI_Comm comm = MPI_COMM_WORLD,
    int64_t reduction_threshold = 0,
    int CA_level = 10000,
    std::vector<int>& idx_map = std::vector<int>(), 
    ReductionPattern pattern = ReductionPattern::UNIFORM,  // NEW parameter
    bool any_process_count = false,  // any rank count (unstructured Color); else (2^dimension)^k
    const RegionWeights& region_weights = RegionWeights{}) {  // the boxes' points (region_weights_of_points), or none
    
    // Validate inputs
    if (dimension < 1 || dimension > 3) {
        throw std::invalid_argument("Dimension must be 1, 2, or 3");
    }
    if (num_levels <= 0) {
        throw std::invalid_argument("Number of levels must be positive");
    }
    
    // Only UNIFORM pattern supported for now
    if (pattern != ReductionPattern::UNIFORM) {
        throw std::invalid_argument("Only UNIFORM reduction pattern currently supported");
    }
    
    // Create tree structure
    auto tree = new ParallelTree<CoordType, DataType>();
    tree->dimension = dimension;
    tree->num_levels = num_levels;
    tree->num_points = num_points;
    tree->reduction_threshold = (reduction_threshold > 0)
        ? reduction_threshold
        : morton::children_per_box(dimension);
    tree->CA_level = (dimension == 1 || CA_level < 0) ? num_levels : CA_level;
    tree->reduction_pattern = pattern;  // NEW
    tree->comm = comm;                  // NEW
    tree->any_process_count = any_process_count;
    
    // Copy global bounds
    std::memcpy(tree->global_bounds, global_bounds, 6 * sizeof(CoordType));
    
    // Get MPI info
    MPI_Comm_rank(comm, &tree->mpi_rank);
    MPI_Comm_size(comm, &tree->mpi_size);
    
    // Allocate levels - use resize instead of new[]
    tree->levels.resize(num_levels);
    
    // Step 1: Determine active processes at each level (Section 5.1.3)
    compute_active_processes(tree, region_weights);

    // Root assembly requires all level-1 sibling boxes on one process. Process
    // reduction from level 2 to level 1 is supported, but level 1 to level 0 is not.
    const int reduction_factor = morton::children_per_box(tree->dimension);
    if (tree->num_levels > 1 && tree->levels[1].num_active_processes != 1) {
        long long cap = 1; 
        for (int i = 0; i < tree->num_levels - 2; ++i) cap *= reduction_factor;
        throw std::invalid_argument(
            "H2 tree: level 1 has " + std::to_string(tree->levels[1].num_active_processes) +
            " active ranks (must be 1). Either reduce ranks to <= " + std::to_string(cap) +
            " (= reduction_factor^(num_levels-2), num_levels=" + std::to_string(tree->num_levels) +
            "), or add tree levels by decreasing Nmin_leaf/increasing the problem size. "
            "Raising reduction_threshold helps only when a level is not already reducing.");
    }
    
    // Step 2-3: Build tree level by level
    for (int32_t level = 0; level < num_levels; ++level) {
        // Distribute boxes to processes (Section 5.1.3)
        distribute_boxes(level, tree);
        
    }

    // Step 4: Assign points to leaf boxes
    if (point_coords != nullptr && num_points > 0) {
        assign_points_to_boxes(point_coords, num_points, idx_map, tree);
    }
    else if (num_points > 0) {
        assign_uniform_points(num_points, tree);
        // std::cout << "num points in box 0: " << tree->levels[num_levels - 1].local_boxes[0].point_indices.size() << std::endl;
    }
    
 

    // Classify interior and boundary boxes for all active levels
    for (int32_t level = 0; level < tree->num_levels; ++level) {
        if (tree->levels[level].is_process_active) {
            populate_interior_boundary(tree->levels[level], tree);
            compute_neighbor_lists(tree->levels[level].local_boxes, tree->dimension, level);
            if (tree->level_uses_CA(level)) {
                compute_ghost_boxes_and_colors(level, tree);
            }
            compute_assisting_boxes(level, tree);
            tree->levels[level].dimension = tree->dimension;  // Store dimension in level for later use
        }
    }
    
    return tree;
}

template<typename CoordType, typename DataType>
void compute_ghost_boxes_and_colors(
    int32_t level,
    ParallelTree<CoordType, DataType>* tree) {
    
    auto& lvl = tree->levels[level];
    
    if (!lvl.is_process_active || lvl.num_boxes_local == 0) {
        return;
    }
    
    // The ghost set from the level's region layout (halo_geometry.hpp): the
    // functions the level-entry halo gather's senders use for this region,
    // so what they send is what this rank expects.
    const RegionLayout L = RegionLayout::of(tree->dimension, level, lvl.num_active_processes);
    const uint32_t my_process_id = static_cast<uint32_t>(lvl.my_morton_id);
    if (!L.uniform || lvl.local_morton_start != L.first(my_process_id) ||
        lvl.num_boxes_local != L.boxes_per_region) {
        throw std::runtime_error("compute_ghost_boxes_and_colors: level " + std::to_string(level) +
                                 " has no uniform region layout (" + std::to_string(lvl.num_active_processes) +
                                 " active processes, " + std::to_string(lvl.num_boxes_local) + " local boxes)");
    }
    // (the region's faces are the boundary boxes populate_interior_boundary found)
    {
        std::vector<int64_t> boundary;
        for (int64_t local_idx : lvl.boundary_id) boundary.push_back(lvl.local_boxes[local_idx].morton_index);
        std::sort(boundary.begin(), boundary.end());
        if (boundary != L.boundary(my_process_id)) {
            throw std::runtime_error("compute_ghost_boxes_and_colors: the boundary boxes are not the region's faces");
        }
    }
    const uint32_t grid_size = L.grid;
    const uint32_t local_grid_size = L.width;
    const std::array<uint32_t, 3> local_offset_array = L.offset(my_process_id);
    const uint32_t* local_offset = local_offset_array.data();
    
    // Helper lambda to get offset for any box's owning process
    auto get_box_owner_offset = [&](int64_t morton_idx) -> std::array<uint32_t, 3> {
        return L.offset_of(morton_idx);
    };
    
    // ========== Steps 1-3: the ghost set, its blue ghosts and orange completions ==========
    // (the one-hop shell of the region, in 3D the two-hop boxes on faces of
    // their regions, and the orange neighbors of the blue ghosts), by region
    // of the ghosts, as the gather's senders compute it for this region;
    // checked against the region's own form and against the original
    // enumeration from the boundary boxes
    std::set<int64_t> temp_ghost_set;
    std::set<int64_t> blue_boxes;
    std::set<int64_t> orange_boxes;
    for_each_region_within(L, my_process_id, halo_region_radius(L), [&](uint32_t r) {
        for (int64_t g : ghosts_of_region_in(L, my_process_id, r)) temp_ghost_set.insert(g);
    });
    {
        const RegionGhosts rg = ghosts_of_region(L, my_process_id);
        const RegionGhosts ref = ghosts_of_region_reference(L, my_process_id);
        if (std::vector<int64_t>(temp_ghost_set.begin(), temp_ghost_set.end()) != rg.ghosts ||
            rg.ghosts != ref.ghosts || rg.blue != ref.blue || rg.orange != ref.orange) {
            throw std::runtime_error("compute_ghost_boxes_and_colors: the ghost sets by region, of the region and "
                                     "by enumeration differ (level " + std::to_string(level) + ")");
        }
        blue_boxes.insert(rg.blue.begin(), rg.blue.end());
        orange_boxes.insert(rg.orange.begin(), rg.orange.end());
    }
    
    // ========== Step 4: Classify remaining ghost boxes (green and purple) ==========
    
    std::set<int64_t> green_boxes;
    std::set<int64_t> purple_boxes;
    
    for (int64_t ghost_morton : temp_ghost_set) {
        // Skip if already classified as blue or orange
        if (blue_boxes.count(ghost_morton) > 0 || orange_boxes.count(ghost_morton) > 0) {
            continue;
        }
        
        // Classify from owning process's perspective
        auto ghost_offset = get_box_owner_offset(ghost_morton);
        
        if (tree->dimension == 3 && is_purple_box(ghost_morton, level, grid_size, local_grid_size, ghost_offset.data())) {
            purple_boxes.insert(ghost_morton);
        } else if (is_orange_box(
                       ghost_morton, level, tree->dimension, grid_size,
                       local_grid_size, ghost_offset.data())) {
            orange_boxes.insert(ghost_morton);
        } else {
            green_boxes.insert(ghost_morton);
        }
    }
    
    // ========== Step 5: Classify local boundary boxes ==========
    
    for (int64_t local_idx : lvl.boundary_id) {
        int64_t box_morton = lvl.local_boxes[local_idx].morton_index;
        
        if (is_blue_box(box_morton, level, tree->dimension, grid_size, local_grid_size, local_offset)) {
            blue_boxes.insert(box_morton);
        } else if (is_orange_box(box_morton, level, tree->dimension, grid_size, local_grid_size, local_offset)) {
            orange_boxes.insert(box_morton);
        } else if (tree->dimension == 3 && is_purple_box(box_morton, level, grid_size, local_grid_size, local_offset)) {
            purple_boxes.insert(box_morton);
        } else {
            green_boxes.insert(box_morton);
        }
    }
    
    // ========== Step 6: Populate ordering vectors ==========
    
    lvl.blue.assign(blue_boxes.begin(), blue_boxes.end());
    lvl.orange.assign(orange_boxes.begin(), orange_boxes.end());
    lvl.green.assign(green_boxes.begin(), green_boxes.end());
    if (tree->dimension == 3) {
        lvl.purple.assign(purple_boxes.begin(), purple_boxes.end());
    }
    std::sort(lvl.blue.begin(), lvl.blue.end());
    std::sort(lvl.orange.begin(), lvl.orange.end());
    std::sort(lvl.green.begin(), lvl.green.end());
    std::sort(lvl.purple.begin(), lvl.purple.end());
    
    // ========== Populate ghost box data structures ==========
    
    // Convert final ghost set to vector
    std::vector<int64_t> temp_ghost(temp_ghost_set.begin(), temp_ghost_set.end());
    std::sort(temp_ghost.begin(), temp_ghost.end());
    
    lvl.ghost_id.clear();
    lvl.ghost_id_to_index.clear();
    
    for (size_t i = 0; i < temp_ghost.size(); ++i) {
        lvl.ghost_id.push_back(temp_ghost[i]);
        lvl.ghost_id_to_index[temp_ghost[i]] = i;
    }
    
    // Allocate ghost box storage
    lvl.ghost_boxes.resize(temp_ghost.size());
    
    // Initialize basic ghost box data
    for (size_t i = 0; i < temp_ghost.size(); ++i) {
        auto& ghost_box = lvl.ghost_boxes[i];
        int64_t morton_idx = temp_ghost[i];
        
        ghost_box.morton_index = morton_idx;
        ghost_box.level = level;
        
        // Decode grid coordinates
        if (tree->dimension == 2) {
            uint32_t x, y;
            morton::decode_2d(morton_idx, x, y);
            ghost_box.grid_coords[0] = static_cast<int32_t>(x);
            ghost_box.grid_coords[1] = static_cast<int32_t>(y);
            ghost_box.grid_coords[2] = 0;
        } else {
            uint32_t x, y, z;
            morton::decode_3d(morton_idx, x, y, z);
            ghost_box.grid_coords[0] = static_cast<int32_t>(x);
            ghost_box.grid_coords[1] = static_cast<int32_t>(y);
            ghost_box.grid_coords[2] = static_cast<int32_t>(z);
        }
        
        // Compute geometry
        compute_box_geometry(morton_idx, level, tree->global_bounds,
                           tree->dimension, ghost_box.center, 
                           ghost_box.size, ghost_box.bounds);
    }
}

} // namespace fmm

#endif // TREE_IMPL_HPP
