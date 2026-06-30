#ifndef TREE_IMPL_HPP
#define TREE_IMPL_HPP

#include "tree.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>
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
      reduction_pattern(ReductionPattern::UNIFORM), reduction_threshold(64) {
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
void clear_modified_interaction_matrices(TreeLevel<CoordType, DataType>& level) {

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

    // Iterate through all local boxes and clear their interaction data.
    for (auto& local_box : level.local_boxes) {
        clear_box_interactions(local_box);
    }

    // Iterate through all ghost boxes and clear their interaction data as well.
    for (auto& ghost_box : level.ghost_boxes) {
        clear_box_interactions(ghost_box);
    }

    std::vector<BoxData<CoordType, DataType>>().swap(level.ghost_boxes);
    level.ghost_id_to_index.clear();
    std::vector<PointDataRequest<CoordType>>().swap(level.assisting_boxes);
    level.assisting_box_points_for_kernel_evaluation.clear();
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
    
    // Decode Morton index to grid coordinates
    uint32_t x, y, z;
    if (dimension == 2) {
        morton::decode_2d(morton_index, x, y);
        z = 0;
    } else {
        morton::decode_3d(morton_index, x, y, z);
    }
    
    // Calculate MAXIMUM dimension range (Section 4.4.3)
    // This ensures boxes are square/cubic for uniform proxy surfaces
    CoordType max_range = global_bounds[1] - global_bounds[0];  // X range
    max_range = std::max(max_range, global_bounds[3] - global_bounds[2]);  // Y range
    if (dimension == 3) {
        max_range = std::max(max_range, global_bounds[5] - global_bounds[4]);  // Z range
    }
    
    // Number of boxes per dimension at this level
    uint32_t boxes_per_dim = 1U << level;  // 2^level
    
    // Box size based on MAXIMUM dimension (ensures square/cubic boxes)
    CoordType box_size = max_range / boxes_per_dim;
    
    // Compute box bounds (positioned relative to domain minimum)
    bounds[0] = global_bounds[0] + x * box_size;  // xmin
    bounds[1] = bounds[0] + box_size;              // xmax
    bounds[2] = global_bounds[2] + y * box_size;  // ymin
    bounds[3] = bounds[2] + box_size;              // ymax
    
    if (dimension == 3) {
        bounds[4] = global_bounds[4] + z * box_size;  // zmin
        bounds[5] = bounds[4] + box_size;              // zmax
    } else {
        bounds[4] = 0.0;
        bounds[5] = 0.0;
    }
    
    // Compute center
    center[0] = bounds[0] + 0.5 * box_size;
    center[1] = bounds[2] + 0.5 * box_size;
    center[2] = (dimension == 3) ? (bounds[4] + 0.5 * box_size) : 0.0;
    
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
    max_range = std::max(max_range, global_bounds[3] - global_bounds[2]);
    if (dimension == 3) {
        max_range = std::max(max_range, global_bounds[5] - global_bounds[4]);
    }
    
    // Number of boxes per dimension
    uint32_t boxes_per_dim = 1U << level;
    CoordType box_size = max_range / boxes_per_dim;
    
    // Compute grid coordinates
    uint32_t x = static_cast<uint32_t>((point[0] - global_bounds[0]) / box_size);
    uint32_t y = static_cast<uint32_t>((point[1] - global_bounds[2]) / box_size);
    uint32_t z = 0;
    
    if (dimension == 3) {
        z = static_cast<uint32_t>((point[2] - global_bounds[4]) / box_size);
    }
    
    // Clamp to valid range
    x = std::min(x, boxes_per_dim - 1);
    y = std::min(y, boxes_per_dim - 1);
    if (dimension == 3) {
        z = std::min(z, boxes_per_dim - 1);
    }
    
    // Store grid coordinates
    grid_coords[0] = static_cast<int32_t>(x);
    grid_coords[1] = static_cast<int32_t>(y);
    grid_coords[2] = static_cast<int32_t>(z);
    
    // Encode to Morton index
    if (dimension == 2) {
        return morton::encode_2d(x, y);
    } else {
        return morton::encode_3d(x, y, z);
    }
}



/**
 * @brief Compute which processes are active at each level (FMM3D.pdf Section 5.3)
 * 
 * Also builds Morton ID mapping and communication patterns for process reduction.
 */
template<typename CoordType, typename DataType>
void compute_active_processes(ParallelTree<CoordType, DataType>* tree) {
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
    
    if (tree->dimension == 2) {
        if (!morton::is_power_of_4(tree->mpi_size)) {
            throw std::invalid_argument("2D tree requires process count to be 4^k (1, 4, 16, 64, ...)");
        }
    } else {
        if (!morton::is_power_of_8(tree->mpi_size)) {
            throw std::invalid_argument("3D tree requires process count to be 8^k (1, 8, 64, 512, ...)");
        }
    }

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
    
    for (int32_t level = leaf_level - 1; level >= 0; --level) {
        auto& lvl = tree->levels[level];
        auto& child_lvl = tree->levels[level + 1];
        lvl.active_process_ranks.clear();
        lvl.morton_to_rank.clear();
        lvl.rank_to_morton.clear();
        lvl.children_senders.clear();
        
        int64_t boxes_at_level = 1LL << (tree->dimension * level);
        int32_t prev_active = child_lvl.num_active_processes;
        int64_t boxes_per_process = boxes_at_level / prev_active;
        int reduction_factor = (tree->dimension == 2) ? 4 : 8;

        switch (tree->reduction_pattern) {
            case ReductionPattern::UNIFORM:
                if (boxes_per_process < tree->reduction_threshold) {
                    lvl.num_active_processes = std::max(1, prev_active / reduction_factor);
                    if (tree->dimension == 3) {
                        if(lvl.num_active_processes == 1 && boxes_per_process >= 512){
                            lvl.num_active_processes = prev_active;
                        }
                    }
                } else {
                    lvl.num_active_processes = prev_active;
                }
                // printf("mpi rank: %d, Level %d: boxes at level: %lld, boxes per process: %lld, active processes: %d\n", tree->mpi_rank,
                //     level, boxes_at_level, boxes_at_level / lvl.num_active_processes, lvl.num_active_processes);
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
            int stride = prev_active / lvl.num_active_processes;

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
                    std::vector<int> group_ranks;
                    group_ranks.reserve(static_cast<size_t>(stride));
                    const int child_morton_start = i * stride;
                    const int child_morton_end = child_morton_start + stride;
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
                int parent_morton = my_child_morton / reduction_factor;
                
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
                int child_morton_start = parent_morton * reduction_factor;
                int child_morton_end = std::min(child_morton_start + reduction_factor - 1,
                                            child_lvl.num_active_processes - 1);
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
        
        // Decode grid coordinates
        if (dimension == 2) {
            uint32_t x, y;
            morton::decode_2d(morton_idx, x, y);
            box.grid_coords[0] = static_cast<int32_t>(x);
            box.grid_coords[1] = static_cast<int32_t>(y);
            box.grid_coords[2] = 0;
        } else {
            uint32_t x, y, z;
            morton::decode_3d(morton_idx, x, y, z);
            box.grid_coords[0] = static_cast<int32_t>(x);
            box.grid_coords[1] = static_cast<int32_t>(y);
            box.grid_coords[2] = static_cast<int32_t>(z);
        }
        
        // Compute geometry using helper function
        compute_box_geometry(morton_idx, level, global_bounds,
                           dimension, box.center, box.size, box.bounds);
        
        // Set parent (if not root)
        if (level > 0) {
            if (dimension == 2) {
                uint32_t x = box.grid_coords[0];
                uint32_t y = box.grid_coords[1];
                box.parent_morton = morton::encode_2d(x / 2, y / 2);
            } else {
                uint32_t x = box.grid_coords[0];
                uint32_t y = box.grid_coords[1];
                uint32_t z = box.grid_coords[2];
                box.parent_morton = morton::encode_3d(x / 2, y / 2, z / 2);
            }
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
    
    // Compute Morton range (Equations 17-21 from Section 5.1.3)
    int64_t boxes_per_process = total_boxes / lvl.num_active_processes;
    int64_t remainder = total_boxes % lvl.num_active_processes;
    
    if (active_index < remainder) {
        //std::cout << "remainder: " << remainder << ", active_index: " << active_index << ", boxes_per_process: " << boxes_per_process << ", total_boxes: " << total_boxes << ", num_active_processes: " << lvl.num_active_processes << std::endl;
        //throw std::runtime_error("Error in box distribution logic, this shouldn't happen");
        std::cerr << "Warning: potentially having fewer boxes than active processes." << std::endl;
        // First 'remainder' processes get one extra box
        lvl.local_morton_start = active_index * (boxes_per_process + 1);
        lvl.local_morton_end = lvl.local_morton_start + boxes_per_process;  // +1 already included
        lvl.num_boxes_local = boxes_per_process + 1;
    } else {
        // Remaining processes get standard allocation
        lvl.local_morton_start = remainder * (boxes_per_process + 1) + 
                                 (active_index - remainder) * boxes_per_process;
        lvl.local_morton_end = lvl.local_morton_start + boxes_per_process - 1;
        lvl.num_boxes_local = boxes_per_process;
    }
    
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
    ParallelTree<CoordType, DataType>* tree) {
    
    if (num_points == 0) {
        return;
    }
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
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
    MPI_Bcast(all_points.data(), tree->dimension * num_points, mpi_type, 0, MPI_COMM_WORLD);
    
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
        CoordType point[3];
        // Column-major (interleaved) format
        point[0] = all_points[pt_idx * tree->dimension];
        point[1] = all_points[pt_idx * tree->dimension + 1];
        point[2] = (tree->dimension == 3) ? all_points[pt_idx * tree->dimension + 2] : 0;
        
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
            box_point_coords[local_idx].push_back(point[0]);
            box_point_coords[local_idx].push_back(point[1]);
            if (tree->dimension == 3) {
                box_point_coords[local_idx].push_back(point[2]);
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
    
    // Compute grid dimension
    int n_per_dim;
    if (tree->dimension == 2) {
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
    CoordType hy = (tree->global_bounds[3] - tree->global_bounds[2]) / n_per_dim;
    CoordType hz = (tree->dimension == 3) ? 
                   (tree->global_bounds[5] - tree->global_bounds[4]) / n_per_dim : 0;
    
    // For each local box, generate points that fall inside
    for (int64_t box_idx = 0; box_idx < leaf_lvl.num_boxes_local; ++box_idx) {
        auto& box = leaf_lvl.local_boxes[box_idx];
        
        // Box bounds
        CoordType x_min = box.bounds[0];
        CoordType x_max = box.bounds[1];
        CoordType y_min = box.bounds[2];
        CoordType y_max = box.bounds[3];
        
        // Find grid indices that overlap with this box
        int i_min = static_cast<int>(std::floor((x_min - tree->global_bounds[0]) / hx));
        int i_max = static_cast<int>(std::floor((x_max - tree->global_bounds[0]) / hx));
        int j_min = static_cast<int>(std::floor((y_min - tree->global_bounds[2]) / hy));
        int j_max = static_cast<int>(std::floor((y_max - tree->global_bounds[2]) / hy));
        
        // Clamp to valid range
        i_min = std::max(0, std::min(n_per_dim - 1, i_min));
        i_max = std::max(0, std::min(n_per_dim - 1, i_max));
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
    std::vector<uint64_t> neighbors;
    if (dimension == 2) {
        neighbors = morton::neighbors_2d(morton_index, grid_size);
    } else {
        neighbors = morton::neighbors_3d(morton_index, grid_size);
    }
    
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
        std::vector<uint64_t> neighbors_1hop;
        if (dimension == 2) {
            neighbors_1hop = morton::neighbors_2d(box.morton_index, grid_size);
        } else {
            neighbors_1hop = morton::neighbors_3d(box.morton_index, grid_size);
        }
        
        // Convert and store
        box.one_hop.clear();
        box.one_hop.reserve(neighbors_1hop.size());
        for (uint64_t n : neighbors_1hop) {
            box.one_hop.push_back(static_cast<int64_t>(n));
        }
        box.use_full_set.clear();
        box.use_full_set.resize(box.one_hop.size(), 1);
        
        // ===== Compute 2-hop neighbors =====
        std::vector<uint64_t> neighbors_2hop;
        if (dimension == 2) {
            neighbors_2hop = morton::neighbors_2hop_2d(box.morton_index, grid_size);
        } else {
            neighbors_2hop = morton::neighbors_2hop_3d(box.morton_index, grid_size);
        }
        
        // Convert and store
        box.two_hop.clear();
        box.two_hop.reserve(neighbors_2hop.size());
        for (uint64_t n : neighbors_2hop) {
            box.two_hop.push_back(static_cast<int64_t>(n));
        }
    }
}

/**
 * @brief Check if a box is a blue box (corner neighborhood) using LOCAL grid coordinates
 * 
 * 2D: Box is one of the 4 corners of the LOCAL process grid
 * 3D: Box is in a 2×2×2 neighborhood around any LOCAL cube corner, excluding diagonal interior
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param dimension 2 or 3
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is blue
 */
inline bool is_blue_box(int64_t morton_index, int32_t level, int32_t dimension, uint32_t grid_size,
                        uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global = 0;
    
    if (dimension == 2) {
        morton::decode_2d(morton_index, x_global, y_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // In 2D, blue boxes are the 4 corners of the LOCAL grid
        return (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) && 
               (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1);
        
    } else {
        morton::decode_3d(morton_index, x_global, y_global, z_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        int32_t z = z_global - local_offset[2];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
            z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Local corners: {0, local_grid_size-1} for each dimension
        uint32_t corners[8][3] = {
            {0, 0, 0}, {local_grid_size-1, 0, 0},
            {0, local_grid_size-1, 0}, {local_grid_size-1, local_grid_size-1, 0},
            {0, 0, local_grid_size-1}, {local_grid_size-1, 0, local_grid_size-1},
            {0, local_grid_size-1, local_grid_size-1}, {local_grid_size-1, local_grid_size-1, local_grid_size-1}
        };
        
        for (int c = 0; c < 8; ++c) {
            uint32_t cx = corners[c][0];
            uint32_t cy = corners[c][1];
            uint32_t cz = corners[c][2];
            
            // Check if box is in 2×2×2 neighborhood (using local coordinates)
            int32_t dx = std::abs(x - static_cast<int32_t>(cx));
            int32_t dy = std::abs(y - static_cast<int32_t>(cy));
            int32_t dz = std::abs(z - static_cast<int32_t>(cz));
            
            if (dx <= 1 && dy <= 1 && dz <= 1) {
                // Exclude diagonal interior box (all distances = 1)
                if (dx == 1 && dy == 1 && dz == 1) {
                    continue;
                }
                return true;
            }
        }
        
        return false;
    }
}


/**
 * @brief Get Manhattan distance-1 neighbors (face-adjacent only)
 * Returns up to 6 neighbors in 3D, up to 4 neighbors in 2D
 */
inline std::vector<uint64_t> get_manhattan_neighbors_3d(uint64_t morton_index, uint32_t grid_size) {
    uint32_t x, y, z;
    morton::decode_3d(morton_index, x, y, z);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(6);
    
    // Check all 6 face directions
    if (x > 0) neighbors.push_back(morton::encode_3d(x-1, y, z));
    if (x < grid_size-1) neighbors.push_back(morton::encode_3d(x+1, y, z));
    if (y > 0) neighbors.push_back(morton::encode_3d(x, y-1, z));
    if (y < grid_size-1) neighbors.push_back(morton::encode_3d(x, y+1, z));
    if (z > 0) neighbors.push_back(morton::encode_3d(x, y, z-1));
    if (z < grid_size-1) neighbors.push_back(morton::encode_3d(x, y, z+1));
    
    return neighbors;
}

inline std::vector<uint64_t> get_manhattan_neighbors_2d(uint64_t morton_index, uint32_t grid_size) {
    uint32_t x, y;
    morton::decode_2d(morton_index, x, y);
    
    std::vector<uint64_t> neighbors;
    neighbors.reserve(4);
    
    if (x > 0) neighbors.push_back(morton::encode_2d(x-1, y));
    if (x < grid_size-1) neighbors.push_back(morton::encode_2d(x+1, y));
    if (y > 0) neighbors.push_back(morton::encode_2d(x, y-1));
    if (y < grid_size-1) neighbors.push_back(morton::encode_2d(x, y+1));
    
    return neighbors;
}

/**
 * @brief Check if a box is an orange box (boundary box adjacent to blue) using LOCAL grid
 * 
 * 2D: LOCAL boundary box at Manhattan distance 1 from any LOCAL corner
 * 3D: ANY boundary box (edge or face) adjacent to at least one blue box
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param dimension 2 or 3
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is orange
 */
inline bool is_orange_box(int64_t morton_index, int32_t level, int32_t dimension, uint32_t grid_size,
                          uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global = 0;
    
    if (dimension == 2) {
        morton::decode_2d(morton_index, x_global, y_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Count LOCAL boundary coordinates
        int boundary_count = 0;
        if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        
        // Must be on boundary
        if (boundary_count == 0) return false;
        
        // Must not be a corner (blue)
        if (boundary_count == 2) return false;
        
        // Check Manhattan distance to each LOCAL corner
        uint32_t corners[4][2] = {
            {0, 0}, {local_grid_size-1, 0}, 
            {0, local_grid_size-1}, {local_grid_size-1, local_grid_size-1}
        };
        
        for (int c = 0; c < 4; ++c) {
            int32_t dist = std::abs(x - static_cast<int32_t>(corners[c][0])) +
                          std::abs(y - static_cast<int32_t>(corners[c][1]));
            if (dist == 1) return true;
        }
        
        return false;
        
    } else {
        morton::decode_3d(morton_index, x_global, y_global, z_global);
        
        // Convert to local coordinates
        int32_t x = x_global - local_offset[0];
        int32_t y = y_global - local_offset[1];
        int32_t z = z_global - local_offset[2];
        
        // Check if within this process's local region
        if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
            y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
            z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
            return false;
        }
        
        // Count LOCAL boundary coordinates
        int boundary_count = 0;
        if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        if (z == 0 || z == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
        
        // Must be on boundary (edge or face)
        if (boundary_count == 0) return false;
        
        // Must be adjacent (Manhattan distance 1) to at least one blue box
        // Check only face-adjacent neighbors (Manhattan distance 1)
        std::vector<uint64_t> neighbors = get_manhattan_neighbors_3d(morton_index, grid_size);
        
        for (uint64_t n_u : neighbors) {
            int64_t n = static_cast<int64_t>(n_u);
            if (is_blue_box(n, level, dimension, grid_size, local_grid_size, local_offset)) {
                return true;
            }
        }
        
        return false;
    }
}

/**
 * @brief Check if a box is a purple box (3D only: edges + face 2nd layer, not orange) using LOCAL grid
 * 
 * 3D: Edge boxes (b=2) OR face 2nd layer boxes (b=1, min_face_dist=1), 
 *     but NOT adjacent to any blue box
 * 
 * @param morton_index Morton index of the box (global)
 * @param level Tree level
 * @param grid_size Global grid size
 * @param local_grid_size Size of local grid per process
 * @param local_offset Starting global coordinates of this process's region
 * @return true if box is purple
 */
inline bool is_purple_box(int64_t morton_index, int32_t level, uint32_t grid_size,
                          uint32_t local_grid_size, const uint32_t local_offset[3]) {
    uint32_t x_global, y_global, z_global;
    morton::decode_3d(morton_index, x_global, y_global, z_global);
    
    // Convert to local coordinates
    int32_t x = x_global - local_offset[0];
    int32_t y = y_global - local_offset[1];
    int32_t z = z_global - local_offset[2];
    
    // Check if within this process's local region
    if (x < 0 || x >= static_cast<int32_t>(local_grid_size) || 
        y < 0 || y >= static_cast<int32_t>(local_grid_size) ||
        z < 0 || z >= static_cast<int32_t>(local_grid_size)) {
        return false;
    }
    
    // Count LOCAL boundary coordinates
    int boundary_count = 0;
    if (x == 0 || x == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    if (y == 0 || y == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    if (z == 0 || z == static_cast<int32_t>(local_grid_size) - 1) boundary_count++;
    
    bool is_candidate = false;
    
    // Case 1: On an edge (boundary_count == 2)
    if (boundary_count == 2) {
        is_candidate = true;
    }
    // Case 2: On a face in the second-to-last layer (boundary_count == 1)
    else if (boundary_count == 1) {
        // Calculate minimum distance to edge within the face (LOCAL coordinates)
        auto dist_to_boundary = [local_grid_size](int32_t coord) -> uint32_t {
            return std::min(coord, static_cast<int32_t>(local_grid_size) - 1 - coord);
        };
        
        uint32_t min_face_dist = local_grid_size;
        if (x != 0 && x != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(x));
        }
        if (y != 0 && y != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(y));
        }
        if (z != 0 && z != static_cast<int32_t>(local_grid_size) - 1) {
            min_face_dist = std::min(min_face_dist, dist_to_boundary(z));
        }
        
        // Must be in second-to-last layer
        if (min_face_dist == 1) {
            is_candidate = true;
        }
    }
    
    if (!is_candidate) return false;
    
    // Must NOT be adjacent to any blue box (otherwise it would be orange)
    std::vector<uint64_t> neighbors = get_manhattan_neighbors_3d(morton_index, grid_size);
    
    for (uint64_t n_u : neighbors) {
        int64_t n = static_cast<int64_t>(n_u);
        if (is_blue_box(n, level, 3, grid_size, local_grid_size, local_offset)) {
            return false;  // Adjacent to blue means it's orange, not purple
        }
    }
    
    return true;
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
    
    uint32_t grid_size = 1 << level;
    uint32_t my_process_id = static_cast<uint32_t>(lvl.my_morton_id);
    
    // Use set to avoid duplicates and for fast lookup
    std::unordered_set<int64_t> assisting_set_for_factorize;
    std::unordered_set<int64_t> assisting_and_ghost_set_for_solve;
    
    // For each ghost box
    for (int64_t ghost_morton : lvl.boundary_id) {
        ghost_morton += lvl.local_morton_start;
        // Get 1-hop neighbors
        std::vector<uint64_t> neighbors_1hop;
        if (tree->dimension == 2) {
            neighbors_1hop = morton::neighbors_2d(ghost_morton, grid_size);
        } else {
            neighbors_1hop = morton::neighbors_3d(ghost_morton, grid_size);
        }
        
        // Batch process assignment for 1-hop
        std::vector<uint32_t> processes_1hop;
        if (tree->dimension == 2) {
            processes_1hop = morton::assign_to_processes_2d(
                neighbors_1hop, lvl.num_active_processes, grid_size);
        } else {
            processes_1hop = morton::assign_to_processes_3d(
                neighbors_1hop, lvl.num_active_processes, grid_size);
        }
        
        // Add 1-hop neighbors that qualify
        for (size_t i = 0; i < neighbors_1hop.size(); ++i) {
            int64_t n = static_cast<int64_t>(neighbors_1hop[i]);
            
            if (processes_1hop[i] != my_process_id) {
                assisting_set_for_factorize.insert(n);
                assisting_and_ghost_set_for_solve.insert(n);
            }
        }
        
        // Get 2-hop neighbors
        std::vector<uint64_t> neighbors_2hop;
        if (tree->dimension == 2) {
            neighbors_2hop = morton::neighbors_2hop_2d(ghost_morton, grid_size);
        } else {
            neighbors_2hop = morton::neighbors_2hop_3d(ghost_morton, grid_size);
        }
        
        // Batch process assignment for 2-hop
        std::vector<uint32_t> processes_2hop;
        if (tree->dimension == 2) {
            processes_2hop = morton::assign_to_processes_2d(
                neighbors_2hop, lvl.num_active_processes, grid_size);
        } else {
            processes_2hop = morton::assign_to_processes_3d(
                neighbors_2hop, lvl.num_active_processes, grid_size);
        }
        
        // Add 2-hop neighbors that qualify for factorization assistance only.
        for (size_t i = 0; i < neighbors_2hop.size(); ++i) {
            int64_t n = static_cast<int64_t>(neighbors_2hop[i]);
            
            if (processes_2hop[i] != my_process_id) {
                assisting_set_for_factorize.insert(n);
            }
        }
    }
    
    // Convert to sorted vector and create morton_index -> vector_index mapping
    std::vector<int64_t> assisting_vec_for_factorize(assisting_set_for_factorize.begin(), assisting_set_for_factorize.end());
    std::sort(assisting_vec_for_factorize.begin(), assisting_vec_for_factorize.end());
    
    lvl.assisting_box_points_for_kernel_evaluation.clear();
    for (size_t i = 0; i < assisting_vec_for_factorize.size(); ++i) {
        lvl.assisting_box_points_for_kernel_evaluation[assisting_vec_for_factorize[i]] = i;
    }

    // construct the set used during solve phase
    std::vector<int64_t> assisting_and_ghost_vec_for_solve(assisting_and_ghost_set_for_solve.begin(), assisting_and_ghost_set_for_solve.end());
    std::sort(assisting_and_ghost_vec_for_solve.begin(), assisting_and_ghost_vec_for_solve.end());
    lvl.ghost_and_assisting_box_points_for_solve_map.clear();
    lvl.is_ghost_solve.clear();
    for (size_t i = 0; i < assisting_and_ghost_vec_for_solve.size(); ++i) {
        lvl.ghost_and_assisting_box_points_for_solve_map[assisting_and_ghost_vec_for_solve[i]] = i;
        lvl.is_ghost_solve.push_back(false);
    }
    lvl.assisting_boxes.resize(lvl.assisting_box_points_for_kernel_evaluation.size());

}

template<typename CoordType, typename DataType>
ParallelTree<CoordType, DataType>* create_uniform_tree(
    const CoordType* point_coords,
    int64_t num_points,
    int32_t num_levels,
    const CoordType global_bounds[6],
    int32_t dimension = 3,
    MPI_Comm comm = MPI_COMM_WORLD,
    int64_t reduction_threshold = 64,
    ReductionPattern pattern = ReductionPattern::UNIFORM) {  // NEW parameter
    
    // Validate inputs
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("Dimension must be 2 or 3");
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
    tree->reduction_threshold = reduction_threshold;
    tree->reduction_pattern = pattern;  // NEW
    tree->comm = comm;                  // NEW
    
    // Copy global bounds
    std::memcpy(tree->global_bounds, global_bounds, 6 * sizeof(CoordType));
    
    // Get MPI info
    MPI_Comm_rank(comm, &tree->mpi_rank);
    MPI_Comm_size(comm, &tree->mpi_size);
    
    // Allocate levels - use resize instead of new[]
    tree->levels.resize(num_levels);
    
    // Step 1: Determine active processes at each level (Section 5.1.3)
    compute_active_processes(tree);
    
    // Step 2-3: Build tree level by level
    for (int32_t level = 0; level < num_levels; ++level) {
        // Distribute boxes to processes (Section 5.1.3)
        distribute_boxes(level, tree);
        
    }
    
    // Step 4: Assign points to leaf boxes
    if (point_coords != nullptr && num_points > 0) {
        assign_points_to_boxes(point_coords, num_points, tree);
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
            compute_assisting_boxes(level, tree);
            tree->levels[level].dimension = tree->dimension;  // Store dimension in level for later use
        }
    }
    
    return tree;
}

} // namespace fmm

#endif // TREE_IMPL_HPP





/*
template<typename CoordType, typename DataType>
void compute_ghost_boxes_and_colors(
    int32_t level,
    ParallelTree<CoordType, DataType>* tree) {
    
    auto& lvl = tree->levels[level];
    
    if (!lvl.is_process_active || lvl.num_boxes_local == 0) {
        return;
    }
    
    uint32_t grid_size = 1 << level;
    std::set<int64_t> temp_ghost_set;  // Use set to avoid duplicates
    
    // Get this process's Morton ID for comparison
    uint32_t my_process_id = static_cast<uint32_t>(lvl.my_morton_id);
    
    // ========== Diagnostic counters ==========
    std::map<uint32_t, int> ghosts_from_step1;  // process_id -> count
    std::map<uint32_t, int> ghosts_from_step2;  // process_id -> count
    std::map<uint32_t, int> ghosts_from_step3;  // process_id -> count
    
    // ========== Compute local grid parameters ==========
    
    // Calculate number of processes per dimension
    uint32_t num_procs_per_dim;
    if (tree->dimension == 2) {
        num_procs_per_dim = 1 << (__builtin_ctz(lvl.num_active_processes) / 2);
    } else {
        num_procs_per_dim = 1 << (__builtin_ctz(lvl.num_active_processes) / 3);
    }
    
    // Local grid size per process
    uint32_t local_grid_size = grid_size / num_procs_per_dim;
    
    // Compute this process's offset in the global grid
    uint32_t local_offset[3] = {0, 0, 0};
    if (lvl.num_boxes_local > 0) {
        // Find minimum coordinates among all local boxes
        local_offset[0] = lvl.local_boxes[0].grid_coords[0];
        local_offset[1] = lvl.local_boxes[0].grid_coords[1];
        local_offset[2] = lvl.local_boxes[0].grid_coords[2];
        
        for (int64_t i = 1; i < lvl.num_boxes_local; ++i) {
            local_offset[0] = std::min(local_offset[0], static_cast<uint32_t>(lvl.local_boxes[i].grid_coords[0]));
            local_offset[1] = std::min(local_offset[1], static_cast<uint32_t>(lvl.local_boxes[i].grid_coords[1]));
            local_offset[2] = std::min(local_offset[2], static_cast<uint32_t>(lvl.local_boxes[i].grid_coords[2]));
        }
    }
    
    // Helper lambda to get offset for any box's owning process
    auto get_box_owner_offset = [&](int64_t morton_idx) -> std::array<uint32_t, 3> {
        uint32_t gx, gy, gz = 0;
        if (tree->dimension == 2) {
            morton::decode_2d(morton_idx, gx, gy);
        } else {
            morton::decode_3d(morton_idx, gx, gy, gz);
        }
        
        return {
            (gx / local_grid_size) * local_grid_size,
            (gy / local_grid_size) * local_grid_size,
            (gz / local_grid_size) * local_grid_size
        };
    };
    
    // ========== Step 1: Get 1-hop neighbors ==========
    
    for (int64_t local_idx : lvl.boundary_id) {
        auto& box = lvl.local_boxes[local_idx];
        
        // Get 1-hop neighbors
        std::vector<uint64_t> neighbors_1hop;
        if (tree->dimension == 2) {
            neighbors_1hop = morton::neighbors_2d(box.morton_index, grid_size);
        } else {
            neighbors_1hop = morton::neighbors_3d(box.morton_index, grid_size);
        }
        
        // Batch process assignment for 1-hop neighbors
        std::vector<uint32_t> neighbor_processes;
        if (tree->dimension == 2) {
            neighbor_processes = morton::assign_to_processes_2d(
                neighbors_1hop, lvl.num_active_processes, grid_size);
        } else {
            neighbor_processes = morton::assign_to_processes_3d(
                neighbors_1hop, lvl.num_active_processes, grid_size);
        }
        
        // Add 1-hop neighbors that are not on this process
        for (size_t i = 0; i < neighbors_1hop.size(); ++i) {
            if (neighbor_processes[i] != my_process_id) {
                int64_t n = static_cast<int64_t>(neighbors_1hop[i]);
                if (temp_ghost_set.insert(n).second) {
                    // Was newly inserted
                    ghosts_from_step1[neighbor_processes[i]]++;
                }
            }
        }
    }
    
    size_t after_step1 = temp_ghost_set.size();
    
    // ========== Step 2: For 3D, add 2-hop neighbors on boundary of their process ==========
    
    if (tree->dimension == 3) {
        for (int64_t local_idx : lvl.boundary_id) {
            auto& box = lvl.local_boxes[local_idx];
            
            std::vector<uint64_t> neighbors_2hop = morton::neighbors_2hop_3d(box.morton_index, grid_size);
            
            // Batch process assignment for 2-hop neighbors
            std::vector<uint32_t> neighbor_2hop_processes = morton::assign_to_processes_3d(
                neighbors_2hop, lvl.num_active_processes, grid_size);
            
            // For each 2-hop neighbor, check if it's on boundary of its process
            for (size_t i = 0; i < neighbors_2hop.size(); ++i) {
                uint64_t n = neighbors_2hop[i];
                uint32_t n_process = neighbor_2hop_processes[i];
                
                // Skip if on this process
                if (n_process == my_process_id) continue;
                
                // Decode global coordinates
                uint32_t nx, ny, nz;
                morton::decode_3d(n, nx, ny, nz);
                
                // Compute local coordinates within the owning process
                uint32_t n_local_x = nx % local_grid_size;
                uint32_t n_local_y = ny % local_grid_size;
                uint32_t n_local_z = nz % local_grid_size;
                
                // Check if on one of the 6 faces of the local cube (boundary of its process)
                bool is_on_face = (n_local_x == 0 || n_local_x == local_grid_size - 1) ||
                                  (n_local_y == 0 || n_local_y == local_grid_size - 1) ||
                                  (n_local_z == 0 || n_local_z == local_grid_size - 1);
                
                if (is_on_face) {
                    if (temp_ghost_set.insert(static_cast<int64_t>(n)).second) {
                        // Was newly inserted
                        ghosts_from_step2[n_process]++;
                    }
                }
            }
        }
    }
    
    size_t after_step2 = temp_ghost_set.size();
    
    // ========== Step 3: Find missing orange boxes ==========
    
    // Convert to vector temporarily to iterate
    std::vector<int64_t> temp_ghost_after_step2(temp_ghost_set.begin(), temp_ghost_set.end());
    
    // Track blue and orange boxes found in step 3
    std::set<int64_t> blue_boxes;
    std::set<int64_t> orange_boxes;
    
    // Find blue boxes in temp_ghost (from their owning process's perspective)
    for (int64_t ghost_morton : temp_ghost_after_step2) {
        auto ghost_offset = get_box_owner_offset(ghost_morton);
        
        if (is_blue_box(ghost_morton, level, tree->dimension, grid_size, local_grid_size, ghost_offset.data())) {
            blue_boxes.insert(ghost_morton);
        }
    }
    
    // For each blue box, find 1-hop neighbors that are orange
    for (int64_t blue_morton : blue_boxes) {
        std::vector<uint64_t> neighbors;
        if (tree->dimension == 2) {
            neighbors = morton::neighbors_2d(blue_morton, grid_size);
        } else {
            neighbors = morton::neighbors_3d(blue_morton, grid_size);
        }
        
        // Batch process assignment for neighbors
        std::vector<uint32_t> neighbor_processes;
        if (tree->dimension == 2) {
            neighbor_processes = morton::assign_to_processes_2d(
                neighbors, lvl.num_active_processes, grid_size);
        } else {
            neighbor_processes = morton::assign_to_processes_3d(
                neighbors, lvl.num_active_processes, grid_size);
        }
        
        for (size_t i = 0; i < neighbors.size(); ++i) {
            int64_t n = static_cast<int64_t>(neighbors[i]);
            
            // Skip if already in temp_ghost
            if (temp_ghost_set.count(n) > 0) continue;
            
            // Skip if on this process
            if (neighbor_processes[i] == my_process_id) continue;
            
            // Check if this neighbor is an orange box (from its owning process's perspective)
            auto n_offset = get_box_owner_offset(n);
            if (is_orange_box(n, level, tree->dimension, grid_size, local_grid_size, n_offset.data())) {
                if (temp_ghost_set.insert(n).second) {
                    // Was newly inserted
                    orange_boxes.insert(n);
                    ghosts_from_step3[neighbor_processes[i]]++;
                }
            }
        }
    }
    
    size_t after_step3 = temp_ghost_set.size();
    
    // ========== Print diagnostics ==========
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        std::cout << "\n========== Ghost Box Diagnostics (Rank 0) ==========" << std::endl;
        std::cout << "After Step 1 (1-hop neighbors): " << after_step1 << " ghosts" << std::endl;
        std::cout << "  Breakdown by source process:" << std::endl;
        for (auto& p : ghosts_from_step1) {
            std::cout << "    Process " << p.first << ": " << p.second << " boxes" << std::endl;
        }
        
        if (tree->dimension == 3) {
            std::cout << "\nAfter Step 2 (2-hop boundary neighbors): " << after_step2 << " ghosts (+";
            std::cout << (after_step2 - after_step1) << ")" << std::endl;
            std::cout << "  Breakdown by source process:" << std::endl;
            for (auto& p : ghosts_from_step2) {
                std::cout << "    Process " << p.first << ": " << p.second << " boxes" << std::endl;
            }
        }
        
        std::cout << "\nAfter Step 3 (orange neighbors of blue): " << after_step3 << " ghosts (+";
        std::cout << (after_step3 - after_step2) << ")" << std::endl;
        std::cout << "  Blue boxes found in ghost region: " << blue_boxes.size() << std::endl;
        std::cout << "  Breakdown by source process:" << std::endl;
        for (auto& p : ghosts_from_step3) {
            std::cout << "    Process " << p.first << ": " << p.second << " boxes" << std::endl;
        }
        
        std::cout << "\nTotal ghost boxes: " << temp_ghost_set.size() << std::endl;
        std::cout << "====================================================\n" << std::endl;
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
    
    // ========== Populate ghost box data structures ==========
    
    // Convert final ghost set to vector
    std::vector<int64_t> temp_ghost(temp_ghost_set.begin(), temp_ghost_set.end());
    
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
*/
