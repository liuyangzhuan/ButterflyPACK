#pragma once

#include "occupied_topology.hpp"

namespace butterfly {
namespace color_unstructured {
using namespace fmm;

template<typename CoordType, typename DataType, typename KernelType>
std::vector<BoxData<CoordType, DataType>> build_parent_level_interactions_unstructured(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    int dimension,
    bool is_symmetric,
    bool is_hermitian,
    KernelType* kernel,
    const CoordType global_bounds[6],
    const OccupiedTopology& occupied_topology) {

    const int num_children = morton::children_per_box(dimension);

    // Calculate number of parent boxes this process owns
    int64_t num_parent_boxes = child_level.local_boxes.size() / num_children;

    if (child_level.local_boxes.size() % num_children != 0) {
        throw std::runtime_error(
            "build_parent_level: Child boxes not evenly divisible by " +
            std::to_string(num_children));
    }

    std::vector<BoxData<CoordType, DataType>> parent_boxes;
    parent_boxes.resize(num_parent_boxes);

    // Calculate parent level number and grid size
    int32_t parent_level_num = child_level.level - 1;
    uint32_t grid_size = 1 << parent_level_num;

    // Calculate parent Morton range for this process
    int64_t local_morton_start = child_level.local_boxes[0].morton_index / num_children;
    int64_t local_morton_end = child_level.local_boxes.back().morton_index / num_children;

    // ===== Step 1: Initialize parent boxes =====

    initialize_local_boxes(
        parent_boxes,
        local_morton_start,
        parent_level_num,
        dimension,
        global_bounds);

    compute_occupied_neighbor_lists(
        parent_boxes, parent_level_num, dimension, occupied_topology);
    const auto occupied_parent_indices = occupied_box_indices_from_topology(
        parent_boxes, parent_level_num, occupied_topology);

    // ===== Step 2: Check boundary conditions =====
    // if no reduction happen, then no multiplier needed, otherwise multiply by 4 (2D) or 8 (3D)
    int transition_multiplier = (parent_level.num_active_processes == child_level.num_active_processes) ? 1 : num_children;
    int64_t start_accounted_for_reduction = parent_level.rank_to_morton[child_level.parent_level_owner] * num_parent_boxes * transition_multiplier;
    int64_t end_accounted_for_reduction = (parent_level.rank_to_morton[child_level.parent_level_owner] + 1) * (num_parent_boxes) * transition_multiplier - 1;
    for (int64_t parent_idx : occupied_parent_indices) {
        auto& parent_box = parent_boxes[static_cast<size_t>(parent_idx)];
        parent_box.on_boundary = check_boundary_condition(
            parent_box.morton_index,
            parent_box.grid_coords,
            grid_size,
            dimension,
            start_accounted_for_reduction,
            end_accounted_for_reduction);
    }

    // for (auto& parent_box : parent_boxes) {
    //     parent_box.on_boundary = check_boundary_condition(
    //         parent_box.morton_index,
    //         parent_box.grid_coords,
    //         grid_size,
    //         dimension,
    //         local_morton_start,
    //         local_morton_end);
    // }

    // ===== Step 3: Accumulate skeleton points from children =====
    for (int64_t parent_idx : occupied_parent_indices) {
        auto& parent_box = parent_boxes[static_cast<size_t>(parent_idx)];
        int64_t first_child_idx = parent_idx * num_children;

        int64_t total_skeleton_points = 0;
        for (int c = 0; c < num_children; ++c) {
            auto& child = child_level.local_boxes[first_child_idx + c];
            total_skeleton_points += child.skeleton_indices.size();
        }

        parent_box.point_indices.reserve(total_skeleton_points);
        parent_box.point_coords.reserve(total_skeleton_points * dimension);

        for (int c = 0; c < num_children; ++c) {
            auto& child = child_level.local_boxes[first_child_idx + c];


            for (int64_t skel_idx : child.skeleton_indices) {
                parent_box.point_indices.push_back(child.point_indices[skel_idx]);
                for (int d = 0; d < dimension; ++d) {
                    parent_box.point_coords.push_back(
                        child.point_coords[skel_idx * dimension + d]);
                }
            }
        }

        parent_box.num_points = total_skeleton_points;
        parent_box.parent_morton = parent_box.morton_index / num_children;

        parent_box.num_children = num_children;
        for (int c = 0; c < num_children; ++c) {
            parent_box.children_morton[c] = child_level.local_boxes[first_child_idx + c].morton_index;
        }
        for (int c = num_children; c < 8; ++c) {
            parent_box.children_morton[c] = -1;
        }
    }

    // ===== Step 4: Build parent-level modified interactions =====

    // Helper struct to store child info (either from ghost or assisting boxes)
    struct ChildInfo {
        BoxData<CoordType, DataType>* box_ptr;  // nullptr if assisting box
        std::vector<int64_t> indices_ptr;
        int64_t num_points;
        bool is_ghost;  // true if from ghost_boxes, false if from assisting_boxes
    };

    struct DeferredParentBlock {
        size_t target_index;
        ModifiedBlock<DataType> block;
    };
    std::vector<std::vector<DeferredParentBlock>> deferred_parent_blocks(
        occupied_parent_indices.size());
    std::exception_ptr build_exception;
    std::mutex build_exception_mutex;
    std::atomic<bool> build_failed{false};
    const bool parallel_build = is_symmetric || is_hermitian;

    #pragma omp parallel for schedule(dynamic) if (parallel_build && occupied_parent_indices.size() > 1)
    for (int64_t occupied_slot = 0;
         occupied_slot < static_cast<int64_t>(occupied_parent_indices.size());
         ++occupied_slot) {
        const size_t b1_idx = static_cast<size_t>(
            occupied_parent_indices[static_cast<size_t>(occupied_slot)]);
        if (build_failed.load(std::memory_order_relaxed)) continue;
        try {
        auto& B1 = parent_boxes[b1_idx];
        if (!occupied_topology.contains(parent_level_num, B1.morton_index)) {
            continue;
        }
        // printf("Assembling interactions for parent box %lu (morton %lu)\n",
        //        b1_idx, B1.morton_index);
        // fflush(stdout);
        // Get relevant neighbors: [self, 1-hop]
        std::vector<uint64_t> relevant_neighbors;
        relevant_neighbors.push_back(B1.morton_index);

        std::vector<uint64_t> one_hop =
            morton::neighbors_nd(dimension, B1.morton_index, grid_size);
        relevant_neighbors.insert(relevant_neighbors.end(), one_hop.begin(), one_hop.end());
        relevant_neighbors.erase(
            std::remove_if(
                relevant_neighbors.begin(), relevant_neighbors.end(),
                [&](uint64_t neighbor_morton) {
                    return !occupied_topology.contains(
                        parent_level_num,
                        static_cast<int64_t>(neighbor_morton));
                }),
            relevant_neighbors.end());

        for (uint64_t neighbor_morton : relevant_neighbors) {
            LazyFarTransitionCache<DataType> pair_lazy_cache;



            // Check if neighbor is on-process or off-process
            bool is_on_process = (neighbor_morton >= local_morton_start &&
                                 neighbor_morton <= local_morton_end);

            // ===== DIAGONAL CASE =====
            if (neighbor_morton == B1.morton_index) {
                int64_t total_rows = B1.num_points;
                int64_t total_cols = B1.num_points;
                std::vector<DataType> I_B1_B1(total_rows * total_cols);

                int64_t row_offset = 0;
                int64_t first_child_b1 = b1_idx * num_children;

                for (int ci = 0; ci < num_children; ++ci) {
                    auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                    int64_t n_i = child_i.skeleton_indices.size();

                    int64_t col_offset = 0;
                    for (int cj = 0; cj < num_children; ++cj) {
                        auto& child_j = child_level.local_boxes[first_child_b1 + cj];
                        int64_t n_j = child_j.skeleton_indices.size();

                        std::vector<DataType> C_block = extract_child_interaction(
                            &child_i, &child_j, child_level, dimension, kernel,
                            false, &pair_lazy_cache);

                        for (int64_t col = 0; col < n_j; ++col) {
                            for (int64_t row = 0; row < n_i; ++row) {
                                I_B1_B1[(row_offset + row) + (col_offset + col) * total_rows] =
                                    C_block[row + col * n_i];
                            }
                        }
                        col_offset += n_j;
                    }
                    row_offset += n_i;
                }

                B1.schur_complement.set_owned(
                    total_rows, total_cols, std::move(I_B1_B1), MatrixStorage<DataType>::FULL);
                continue;
            }

            // ===== OFF-DIAGONAL CASES =====

            if (is_symmetric || is_hermitian) {
                // ========== SYMMETRIC CASE ==========

                if (is_on_process) {
                    // Both B1 and B2 on process - lower triangular optimization
                    if (neighbor_morton < B1.morton_index) {
                        continue; // Skip, already handled
                    }

                    size_t b2_idx = neighbor_morton - local_morton_start;
                    auto& B2 = parent_boxes[b2_idx];

                    int64_t total_rows = B1.num_points;
                    int64_t total_cols = B2.num_points;
                    std::vector<DataType> I_B1_B2(total_rows * total_cols);

                    // Assemble I(B1, B2)
                    int64_t row_offset = 0;
                    int64_t first_child_b1 = b1_idx * num_children;

                    for (int ci = 0; ci < num_children; ++ci) {
                        auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                        int64_t n_i = child_i.skeleton_indices.size();

                        int64_t col_offset = 0;
                        int64_t first_child_b2 = b2_idx * num_children;

                        for (int cj = 0; cj < num_children; ++cj) {
                            auto& child_j = child_level.local_boxes[first_child_b2 + cj];
                            int64_t n_j = child_j.skeleton_indices.size();

                            std::vector<DataType> C_block = extract_child_interaction(
                                &child_i, &child_j, child_level, dimension, kernel,
                                false, &pair_lazy_cache);

                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B1_B2[(row_offset + row) + (col_offset + col) * total_rows] =
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }

                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = B2.morton_index;
                    ModifiedBlock<DataType> block_b2;
                    block_b2.neighbor_morton = B1.morton_index;

                    if (is_symmetric && !is_hermitian) {
                        // B2 owns the already assembled payload; B1 is its
                        // transpose view of the same allocation.
                        block_b2.set_a_ns_owned(
                            total_rows, total_cols, std::move(I_B1_B2),
                            MatrixStorage<DataType>::FULL);
                        if (!share_symmetric_a_ns_pair(block_b2, block_b1)) {
                            throw std::runtime_error(
                                "build_parent_level: cannot share symmetric parent edge");
                        }
                    } else {
                        // Preserve the directed Hermitian path until conjugate
                        // views are represented explicitly.
                        std::vector<DataType> I_transposed(
                            static_cast<size_t>(total_rows * total_cols));
                        for (int64_t i = 0; i < total_rows; ++i) {
                            for (int64_t j = 0; j < total_cols; ++j) {
                                I_transposed[static_cast<size_t>(
                                    j + i * total_cols)] =
                                    I_B1_B2[static_cast<size_t>(
                                        i + j * total_rows)];
                            }
                        }
                        block_b1.set_a_ns_owned(
                            total_cols, total_rows, std::move(I_transposed),
                            MatrixStorage<DataType>::FULL);
                        block_b2.set_a_ns_owned(
                            total_rows, total_cols, std::move(I_B1_B2),
                            MatrixStorage<DataType>::FULL);
                    }

                    int64_t block_idx_b1 =
                        B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(
                        std::move(block_b1));
                    B1.near_field_interaction_map[B2.morton_index] =
                        block_idx_b1;

                    deferred_parent_blocks[static_cast<size_t>(occupied_slot)].push_back(
                        DeferredParentBlock{b2_idx, std::move(block_b2)});

                } else {

                    // B2 is off-process (ghost or assisting) - store only B1's view

                    // Find B2's children (check ghost first, then assisting)
                    std::vector<ChildInfo> b2_children;
                    for (int c = 0; c < num_children; ++c) {
                        int64_t child_morton = neighbor_morton * num_children + c;

                        // Check ghost boxes first
                        auto ghost_it = child_level.ghost_id_to_index.find(child_morton);
                        if (ghost_it != child_level.ghost_id_to_index.end()) {
                            auto& ghost_box = child_level.ghost_boxes[ghost_it->second];
                            ChildInfo info;
                            info.box_ptr = &ghost_box;
                            info.indices_ptr = ghost_box.point_indices;
                            info.num_points = ghost_box.skeleton_indices.size();
                            info.is_ghost = true;
                            b2_children.push_back(info);
                            continue;
                        }

                        // Check assisting boxes
                        // printf("map size: %lu\n", child_level.assisting_box_points_for_kernel_evaluation.size());
                        auto assist_it = child_level.assisting_box_points_for_kernel_evaluation.find(child_morton);
                        if (assist_it != child_level.assisting_box_points_for_kernel_evaluation.end()) {

                            auto& assist_box = child_level.assisting_boxes[assist_it->second];
                            ChildInfo info;
                            info.box_ptr = nullptr;
                            info.indices_ptr = extract_skeleton_indices(assist_box, dimension);
                            info.num_points = assist_box.skel_indices.size();
                            info.is_ghost = false;
                            b2_children.push_back(info);
                            continue;
                        }



                        // Not found anywhere
                        throw std::runtime_error(
                            "build_parent_level: Child " + std::to_string(child_morton) +
                            " not found in ghost_boxes or assisting_boxes");
                    }

                    // Calculate B2's total points
                    int64_t total_cols = 0;
                    for (const auto& info : b2_children) {
                        total_cols += info.num_points;
                    }

                    int64_t total_rows = B1.num_points;
                    std::vector<DataType> I_B1_B2(total_rows * total_cols);

                    // Assemble I(B1, B2) using only B1's children
                    int64_t row_offset = 0;
                    int64_t first_child_b1 = b1_idx * num_children;

                    std::vector<DataType> C_block;
                    for (int ci = 0; ci < num_children; ++ci) {
                        auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                        int64_t n_i = child_i.skeleton_indices.size();

                        int64_t col_offset = 0;
                        for (int cj = 0; cj < num_children; ++cj) {
                            const auto& child_j_info = b2_children[cj];
                            int64_t n_j = child_j_info.num_points;



                            if (child_j_info.is_ghost) {
                                // Use extract_child_interaction for ghost boxes
                                C_block = extract_child_interaction(
                                    &child_i, child_j_info.box_ptr, child_level,
                                    dimension, kernel, false, &pair_lazy_cache);

                            } else {
                                // // Direct kernel evaluation for assisting boxes
                                // C_block.resize(n_i * n_j);

                                // // Extract skeleton coords from child_i
                                // std::vector<CoordType> child_i_skel_coords(n_i * dimension);
                                // for (int64_t idx = 0; idx < n_i; ++idx) {
                                //     int64_t skel_idx = child_i.skeleton_indices[idx];
                                //     for (int d = 0; d < dimension; ++d) {
                                //         child_i_skel_coords[idx * dimension + d] =
                                //             child_i.point_coords[skel_idx * dimension + d];
                                //     }
                                // }

                                // kernel->evaluate_block(
                                //     child_i_skel_coords.data(), n_i,
                                //     child_j_info.coords_ptr.data(), n_j,
                                //     C_block.data(), n_i);
                                // Assisting box - try to find in child_i's modified interactions first
                                int64_t child_j_morton = neighbor_morton * num_children + cj;
                                C_block = extract_or_evaluate_child_interaction_for_assisting(
                                    &child_i,           // Looking in child_i's storage
                                    child_j_morton,     // For neighbor child_j
                                    child_j_info.indices_ptr,
                                    n_i,
                                    n_j,
                                    child_level,
                                    dimension,
                                    kernel,
                                    true, // transpose source-owned A_NS
                                    &pair_lazy_cache
                                );
                            }

                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B1_B2[(row_offset + row) + (col_offset + col) * total_rows] =
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }

                    // Transpose for B1's A_NS view (B2 × B1)
                    std::vector<DataType> I_transposed(total_rows * total_cols);
                    for (int64_t i = 0; i < total_rows; ++i) {
                        for (int64_t j = 0; j < total_cols; ++j) {
                            I_transposed[j + i * total_cols] = I_B1_B2[i + j * total_rows];
                        }
                    }

                    // Store only B1's view
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = neighbor_morton;
                    block_b1.set_a_ns_owned(
                        total_cols, total_rows, std::move(I_transposed), MatrixStorage<DataType>::FULL);

                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map[neighbor_morton] = block_idx_b1;
                }

            } else {
                // ========== NONSYMMETRIC CASE ==========

                if (is_on_process) {
                    // Both on process - always process (no triangular optimization)

                    size_t b2_idx = neighbor_morton - local_morton_start;
                    auto& B2 = parent_boxes[b2_idx];

                    int64_t total_rows = B1.num_points;
                    int64_t total_cols = B2.num_points;
                    std::vector<DataType> I_B1_B2(total_rows * total_cols);

                    // Assemble I(B1, B2)
                    int64_t row_offset = 0;
                    int64_t first_child_b1 = b1_idx * num_children;

                    for (int ci = 0; ci < num_children; ++ci) {
                        auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                        int64_t n_i = child_i.skeleton_indices.size();

                        int64_t col_offset = 0;
                        int64_t first_child_b2 = b2_idx * num_children;

                        for (int cj = 0; cj < num_children; ++cj) {
                            auto& child_j = child_level.local_boxes[first_child_b2 + cj];
                            int64_t n_j = child_j.skeleton_indices.size();

                            std::vector<DataType> C_block = extract_child_interaction(
                                &child_i, &child_j, child_level, dimension, kernel,
                                false, &pair_lazy_cache);

                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B1_B2[(row_offset + row) + (col_offset + col) * total_rows] =
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }

                    // Transpose for A_NS (B2 × B1)
                    std::vector<DataType> I_transposed(total_rows * total_cols);
                    for (int64_t i = 0; i < total_rows; ++i) {
                        for (int64_t j = 0; j < total_cols; ++j) {
                            I_transposed[j + i * total_cols] = I_B1_B2[i + j * total_rows];
                        }
                    }

                    // Store A_NS for B1
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = B2.morton_index;
                    block_b1.A_NS.allocate(total_cols, total_rows, MatrixStorage<DataType>::FULL);
                    block_b1.A_NS.data = I_transposed;

                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map_nonsymmetry[B2.morton_index] = block_idx_b1;

                    // Store A_SN for B2
                    ModifiedBlock<DataType> block_b2;
                    block_b2.neighbor_morton = B1.morton_index;
                    block_b2.A_SN.set_owned(
                        total_cols, total_rows, std::move(I_transposed), MatrixStorage<DataType>::FULL);

                    int64_t block_idx_b2 = B2.near_field_modified_interactions.size();
                    B2.near_field_modified_interactions.push_back(std::move(block_b2));
                    B2.near_field_interaction_map_nonsymmetry[B1.morton_index] = block_idx_b2;

                } else {
                    // B2 is off-process (ghost or assisting) - store both A_NS and A_SN for B1

                    // Find B2's children (check ghost first, then assisting)
                    std::vector<ChildInfo> b2_children;
                    for (int c = 0; c < num_children; ++c) {
                        int64_t child_morton = neighbor_morton * num_children + c;

                        // Check ghost boxes first
                        auto ghost_it = child_level.ghost_id_to_index.find(child_morton);
                        if (ghost_it != child_level.ghost_id_to_index.end()) {
                            auto& ghost_box = child_level.ghost_boxes[ghost_it->second];
                            ChildInfo info;
                            info.box_ptr = &ghost_box;
                            info.indices_ptr = ghost_box.point_indices;
                            info.num_points = ghost_box.skeleton_indices.size();
                            info.is_ghost = true;
                            b2_children.push_back(info);
                            continue;
                        }

                        // Check assisting boxes
                        auto assist_it = child_level.assisting_box_points_for_kernel_evaluation.find(child_morton);
                        if (assist_it != child_level.assisting_box_points_for_kernel_evaluation.end()) {
                            auto& assist_box = child_level.assisting_boxes[assist_it->second];
                            ChildInfo info;
                            info.box_ptr = nullptr;
                            info.indices_ptr = extract_skeleton_indices(assist_box, dimension);
                            info.num_points = assist_box.skel_indices.size();
                            info.is_ghost = false;
                            b2_children.push_back(info);
                            continue;
                        }

                        // Not found anywhere
                        throw std::runtime_error(
                            "build_parent_level: Child " + std::to_string(child_morton) +
                            " not found in ghost_boxes or assisting_boxes");
                    }

                    // Calculate B2's total points
                    int64_t total_cols = 0;
                    for (const auto& info : b2_children) {
                        total_cols += info.num_points;
                    }

                    int64_t total_rows = B1.num_points;
                    int64_t first_child_b1 = b1_idx * num_children;

                    // ===== Assemble I(B1, B2) for A_SN =====
                    std::vector<DataType> I_B1_B2(total_rows * total_cols);

                    int64_t row_offset = 0;
                    for (int ci = 0; ci < num_children; ++ci) {
                        auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                        int64_t n_i = child_i.skeleton_indices.size();

                        // Extract skeleton coords from child_i
                        std::vector<CoordType> child_i_skel_coords(n_i * dimension);
                        for (int64_t idx = 0; idx < n_i; ++idx) {
                            int64_t skel_idx = child_i.skeleton_indices[idx];
                            for (int d = 0; d < dimension; ++d) {
                                child_i_skel_coords[idx * dimension + d] =
                                    child_i.point_coords[skel_idx * dimension + d];
                            }
                        }

                        int64_t col_offset = 0;
                        for (int cj = 0; cj < num_children; ++cj) {
                            const auto& child_j_info = b2_children[cj];
                            int64_t n_j = child_j_info.num_points;

                            std::vector<DataType> C_block;

                            if (child_j_info.is_ghost) {
                                C_block = extract_child_interaction(
                                    &child_i, child_j_info.box_ptr, child_level,
                                    dimension, kernel, false, &pair_lazy_cache);
                            } else {
                                // // Direct kernel evaluation
                                // C_block.resize(n_i * n_j);
                                // kernel->evaluate_block(
                                //     child_i_skel_coords.data(), n_i,
                                //     child_j_info.coords_ptr.data(), n_j,
                                //     C_block.data(), n_i);
                                // Assisting box - try to find in child_i's modified interactions first
                                int64_t child_j_morton = neighbor_morton * num_children + cj;
                                C_block = extract_or_evaluate_child_interaction_for_assisting(
                                    &child_i,           // Looking in child_i's storage
                                    child_j_morton,     // For neighbor child_j
                                    child_j_info.indices_ptr,
                                    n_i,
                                    n_j,
                                    child_level,
                                    dimension,
                                    kernel,
                                    true,
                                    &pair_lazy_cache
                                );
                            }

                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B1_B2[(row_offset + row) + (col_offset + col) * total_rows] =
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }

                    // ===== Assemble I(B2, B1) for A_NS =====
                    std::vector<DataType> I_B2_B1(total_cols * total_rows);

                    row_offset = 0;
                    for (int ci = 0; ci < num_children; ++ci) {
                        const auto& child_i_info = b2_children[ci];
                        int64_t n_i = child_i_info.num_points;

                        int64_t col_offset = 0;
                        for (int cj = 0; cj < num_children; ++cj) {
                            auto& child_j = child_level.local_boxes[first_child_b1 + cj];
                            int64_t n_j = child_j.skeleton_indices.size();

                            // Extract skeleton coords from child_j
                            std::vector<CoordType> child_j_skel_coords(n_j * dimension);
                            for (int64_t idx = 0; idx < n_j; ++idx) {
                                int64_t skel_idx = child_j.skeleton_indices[idx];
                                for (int d = 0; d < dimension; ++d) {
                                    child_j_skel_coords[idx * dimension + d] =
                                        child_j.point_coords[skel_idx * dimension + d];
                                }
                            }

                            std::vector<DataType> C_block;

                            if (child_i_info.is_ghost) {
                                C_block = extract_child_interaction(
                                    child_i_info.box_ptr, &child_j, child_level,
                                    dimension, kernel, false, &pair_lazy_cache);
                            } else {
                                // // Direct kernel evaluation
                                // C_block.resize(n_i * n_j);
                                // kernel->evaluate_block(
                                //     child_i_info.coords_ptr.data(), n_i,
                                //     child_j_skel_coords.data(), n_j,
                                //     C_block.data(), n_i);
                                // Assisting box - try to find in child_j's modified interactions first
                                int64_t child_i_morton = neighbor_morton * num_children + ci;

                                // Looking in child_j for child_i
                                // child_j.A_NS[child_i] is (n_i × n_j)
                                // Function call: extract from child_j (source), looking for child_i (target)
                                // Returns: (n_source × n_target) = (n_j × n_i) after transpose
                                C_block = extract_or_evaluate_child_interaction_for_assisting(
                                    &child_j,           // Looking in child_j's storage
                                    child_i_morton,     // For neighbor child_i
                                    child_i_info.indices_ptr,
                                    n_j,                // n_source = child_j skeleton size
                                    n_i,                // n_target = child_i skeleton size
                                    child_level,
                                    dimension,
                                    kernel,
                                    true,
                                    &pair_lazy_cache
                                );

                                // C_block is now (n_j × n_i), but we need (n_i × n_j) for I(B2, B1)
                                // Transpose again
                                std::vector<DataType> C_block_transposed(n_i * n_j);
                                for (int64_t i = 0; i < n_j; ++i) {
                                    for (int64_t j = 0; j < n_i; ++j) {
                                        C_block_transposed[j + i * n_i] = C_block[i + j * n_j];
                                    }
                                }
                                C_block = std::move(C_block_transposed);
                            }

                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B2_B1[(row_offset + row) + (col_offset + col) * total_cols] =
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }

                    // Transpose I_B2_B1 for A_NS (B1 × B2)
                    std::vector<DataType> I_B2_B1_transposed(total_cols * total_rows);
                    for (int64_t i = 0; i < total_cols; ++i) {
                        for (int64_t j = 0; j < total_rows; ++j) {
                            I_B2_B1_transposed[j + i * total_rows] = I_B2_B1[i + j * total_cols];
                        }
                    }

                    // Store both blocks for B1
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = neighbor_morton;

                    // A_NS (B1 × B2)
                    block_b1.A_NS.set_owned(
                        total_rows, total_cols, std::move(I_B2_B1_transposed), MatrixStorage<DataType>::FULL);

                    // A_SN (B1 × B2)
                    block_b1.A_SN.set_owned(
                        total_rows, total_cols, std::move(I_B1_B2), MatrixStorage<DataType>::FULL);

                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map_nonsymmetry[neighbor_morton] = block_idx_b1;
                }
            }
        }
        } catch (...) {
            if (!build_failed.exchange(true, std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(build_exception_mutex);
                build_exception = std::current_exception();
            }
        }
    }

    if (build_exception) std::rethrow_exception(build_exception);

    if (parallel_build) {
        for (size_t occupied_slot = 0;
             occupied_slot < deferred_parent_blocks.size(); ++occupied_slot) {
            const size_t b1_idx = static_cast<size_t>(
                occupied_parent_indices[occupied_slot]);
            const int64_t source_morton = parent_boxes[b1_idx].morton_index;
            for (auto& deferred : deferred_parent_blocks[occupied_slot]) {
                auto& target = parent_boxes[deferred.target_index];
                const int64_t block_index = static_cast<int64_t>(
                    target.near_field_modified_interactions.size());
                target.near_field_modified_interactions.push_back(
                    std::move(deferred.block));
                target.near_field_interaction_map[source_morton] = block_index;
            }
        }
    }

    return parent_boxes;
}

}  // namespace color_unstructured
}  // namespace butterfly
