#pragma once

#include "occupied_topology.hpp"

#include <numeric>

namespace butterfly {
namespace color_unstructured {
using namespace fmm;

// Steps 1-3 of the transition: the parent boxes of this rank's children,
// their occupied neighbor lists, boundary flags and skeleton points, without
// their blocks.
template<typename CoordType, typename DataType>
std::vector<BoxData<CoordType, DataType>> build_parent_level_structure_unstructured(
    const TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    int dimension,
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
    // the range of the parent region these parents go to (this rank's, or with
    // a reduction its group's): the parent level's region ranges
    // (region_partition.hpp; the equal ones for (2^d)^k ranks)
    const int parent_region = parent_level.rank_to_morton.at(child_level.parent_level_owner);
    const int64_t start_accounted_for_reduction = parent_level.region_start.at(static_cast<size_t>(parent_region));
    const int64_t end_accounted_for_reduction = parent_level.region_start.at(static_cast<size_t>(parent_region) + 1) - 1;
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
            const auto& child = child_level.local_boxes[first_child_idx + c];
            total_skeleton_points += child.skeleton_indices.size();
        }

        parent_box.point_indices.reserve(total_skeleton_points);
        parent_box.point_coords.reserve(total_skeleton_points * dimension);

        for (int c = 0; c < num_children; ++c) {
            const auto& child = child_level.local_boxes[first_child_idx + c];


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
    return parent_boxes;
}

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

    // ===== Steps 1-3: parent boxes, neighbors, boundary flags, points =====
    std::vector<BoxData<CoordType, DataType>> parent_boxes =
        build_parent_level_structure_unstructured(
            child_level, parent_level, dimension, global_bounds,
            occupied_topology);
    int32_t parent_level_num = child_level.level - 1;
    uint32_t grid_size = 1 << parent_level_num;
    int64_t local_morton_start = child_level.local_boxes[0].morton_index / num_children;
    int64_t local_morton_end = child_level.local_boxes.back().morton_index / num_children;
    const auto occupied_parent_indices = occupied_box_indices_from_topology(
        parent_boxes, parent_level_num, occupied_topology);

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
        bool a_sn = false;  // (a nonsymmetric matrix's A_SN view)
    };
    std::vector<std::vector<DeferredParentBlock>> deferred_parent_blocks(
        occupied_parent_indices.size());
    std::exception_ptr build_exception;
    std::mutex build_exception_mutex;
    std::atomic<bool> build_failed{false};
    // (a nonsymmetric matrix: a parent builds its views, and those of an
    // on-process neighbor above it)
    const bool parallel_build = true;
    // (its kernel blocks split over the threads a few parents leave idle)
    const int transition_split = split_threads_for(
        static_cast<int64_t>(occupied_parent_indices.size()), std::max(1, omp_get_max_threads()));

    // With entries through a collective callback, the parents that asked
    // for a block not fetched yet run again after a round, from a fresh box:
    // an iteration writes only its parent's Schur complement and near blocks
    // (and its own deferred blocks for the others, a nonsymmetric matrix's
    // views of an on-process neighbor among them), which start empty.
    auto* entries = kernel_entries(kernel, 0);
    if (entries != nullptr) {
        for (int64_t parent_index : occupied_parent_indices) {
            const auto& B = parent_boxes[static_cast<size_t>(parent_index)];
            if (B.schur_complement.is_allocated() || !B.near_field_modified_interactions.empty() ||
                !B.near_field_interaction_map.empty() || !B.near_field_interaction_map_nonsymmetry.empty()) {
                throw std::logic_error(
                    "build_parent_level_interactions_unstructured: a parent box not fresh");
            }
        }
    }
    std::vector<int64_t> parent_items(occupied_parent_indices.size());
    std::iota(parent_items.begin(), parent_items.end(), int64_t{0});
    run_entry_phase(entries, std::move(parent_items), [&](const std::vector<int64_t>& items) {
    std::vector<int64_t> abandoned;
    std::mutex abandoned_mutex;
    #pragma omp parallel for schedule(dynamic) if (parallel_build && items.size() > 1)
    for (int64_t item = 0; item < static_cast<int64_t>(items.size()); ++item) {
        const int64_t occupied_slot = items[static_cast<size_t>(item)];
        const size_t b1_idx = static_cast<size_t>(
            occupied_parent_indices[static_cast<size_t>(occupied_slot)]);
        if (build_failed.load(std::memory_order_relaxed)) continue;
        try {
        auto& B1 = parent_boxes[b1_idx];
        if (!occupied_topology.contains(parent_level_num, B1.morton_index)) {
            continue;
        }
        if (entries != nullptr) {
            B1.schur_complement = MatrixStorage<DataType>{};
            B1.near_field_modified_interactions.clear();
            B1.near_field_interaction_map.clear();
            B1.near_field_interaction_map_nonsymmetry.clear();
            deferred_parent_blocks[static_cast<size_t>(occupied_slot)].clear();
        }
        entry_item_begin();
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

        // A nonsymmetric matrix: the children of B1 and their skeletons, and
        // the block of a child pair from the children's directed views
        struct ChildRef {
            int64_t morton = -1;
            BoxData<CoordType, DataType>* box = nullptr;  // local or ghost
            std::vector<int64_t> skeleton;                  // global indices
            int64_t full_size = 0;                          // (its points)
            const std::vector<int64_t>* skeleton_positions = nullptr;
        };
        auto child_ref = [&](BoxData<CoordType, DataType>* box) {
            ChildRef ref;
            ref.morton = box->morton_index;
            ref.box = box;
            ref.full_size = box->num_points;
            ref.skeleton_positions = &box->skeleton_indices;
            ref.skeleton.resize(box->skeleton_indices.size());
            for (size_t s_idx = 0; s_idx < box->skeleton_indices.size(); ++s_idx) {
                ref.skeleton[s_idx] = box->point_indices[static_cast<size_t>(box->skeleton_indices[s_idx])];
            }
            return ref;
        };
        std::vector<ChildRef> b1_children;
        const int64_t first_child_b1_ns = static_cast<int64_t>(b1_idx) * num_children;
        for (int c = 0; c < num_children && !(is_symmetric || is_hermitian); ++c) {
            b1_children.push_back(child_ref(&child_level.local_boxes[static_cast<size_t>(first_child_b1_ns + c)]));
        }
        // K(a, b) at the skeletons, n_a x n_b
        std::vector<DataType> C_block;
        auto child_pair_block = [&](const ChildRef& a, const ChildRef& b) {
            const int64_t n_a = static_cast<int64_t>(a.skeleton.size());
            const int64_t n_b = static_cast<int64_t>(b.skeleton.size());
            if (n_a == 0 || n_b == 0) {
                C_block.clear();
                return;
            }
            if (a.morton == b.morton) {
                // the child's Schur block, untransposed
                C_block.resize(static_cast<size_t>(n_a * n_a));
                const auto* schur = a.box != nullptr ? &a.box->schur_complement : nullptr;
                if (schur != nullptr && schur->is_allocated()) {
                    if (schur->rows != n_a || schur->cols != n_a) {
                        throw std::runtime_error("build_parent_level: child " + std::to_string(a.morton) +
                                                 " Schur block not at its skeleton");
                    }
                    for (int64_t col = 0; col < n_a; ++col) {
                        for (int64_t row = 0; row < n_a; ++row) {
                            C_block[static_cast<size_t>(row + col * n_a)] = (*schur)(row, col);
                        }
                    }
                } else {
                    kernel->evaluate_block_by_index(
                        a.skeleton.data(), n_a, a.skeleton.data(), n_a, C_block.data(), n_a);
                }
                return;
            }
            uint32_t ac[3] = {0, 0, 0};
            uint32_t bc[3] = {0, 0, 0};
            morton::decode_nd(dimension, a.morton, ac[0], ac[1], ac[2]);
            morton::decode_nd(dimension, b.morton, bc[0], bc[1], bc[2]);
            int64_t hops = 0;
            for (int d = 0; d < dimension; ++d) {
                hops = std::max<int64_t>(
                    hops, std::abs(static_cast<int64_t>(ac[d]) - static_cast<int64_t>(bc[d])));
            }
            if (hops == 2 && lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                // lazy fill: no far view is stored; the kernel less the
                // fill of the pair's eliminated common neighbors
                std::vector<LazyFarSource> sources;
                if (a.box != nullptr) {
                    collect_lazy_far_sources(child_level, a.box, b.morton, dimension, sources);
                } else if (b.box != nullptr) {
                    collect_lazy_far_sources(child_level, b.box, a.morton, dimension, sources);
                }
                auto endpoint = [](const ChildRef& c) {
                    LazyFarEndpoint<CoordType> e;
                    e.morton = c.morton;
                    e.full_size = c.full_size;
                    e.skeleton = c.skeleton_positions;
                    e.wanted = c.skeleton_positions;
                    e.wanted_count = static_cast<int64_t>(c.skeleton.size());
                    e.indices = c.skeleton.data();
                    return e;
                };
                C_block.resize(static_cast<size_t>(n_a * n_b));
                regenerate_far_block_into(
                    child_level, kernel, endpoint(a), endpoint(b), sources, C_block.data(), n_a);
                return;
            }
            if (hops <= 2) {
                const EdgeKind kind = hops <= 1 ? EdgeKind::Near : EdgeKind::Far;
                if (a.box != nullptr &&
                    read_directed_view(child_level, *a.box, b.morton, kind, true, n_a, n_b, C_block)) {
                    return;
                }
                if (b.box != nullptr &&
                    read_directed_view(child_level, *b.box, a.morton, kind, false, n_b, n_a, C_block)) {
                    return;
                }
            }
            C_block.resize(static_cast<size_t>(n_a * n_b));
            evaluate_block_by_index_split<DataType>(
                kernel, a.skeleton.data(), n_a, b.skeleton.data(), n_b, C_block.data(), n_a,
                transition_split);
        };
        auto assemble = [&](const std::vector<ChildRef>& row_children,
                            const std::vector<ChildRef>& col_children) {
            int64_t total_rows = 0;
            int64_t total_cols = 0;
            for (const auto& c : row_children) total_rows += static_cast<int64_t>(c.skeleton.size());
            for (const auto& c : col_children) total_cols += static_cast<int64_t>(c.skeleton.size());
            std::vector<DataType> I(static_cast<size_t>(total_rows * total_cols));
            int64_t row_offset = 0;
            for (const auto& a : row_children) {
                const int64_t n_a = static_cast<int64_t>(a.skeleton.size());
                int64_t col_offset = 0;
                for (const auto& b : col_children) {
                    const int64_t n_b = static_cast<int64_t>(b.skeleton.size());
                    child_pair_block(a, b);
                    for (int64_t col = 0; col < n_b; ++col) {
                        for (int64_t row = 0; row < n_a; ++row) {
                            I[static_cast<size_t>((row_offset + row) + (col_offset + col) * total_rows)] =
                                C_block[static_cast<size_t>(row + col * n_a)];
                        }
                    }
                    col_offset += n_b;
                }
                row_offset += n_a;
            }
            MatrixStorage<DataType> m;
            m.set_owned(total_rows, total_cols, std::move(I), MatrixStorage<DataType>::FULL);
            return m;
        };

        for (uint64_t neighbor_morton : relevant_neighbors) {
            LazyFarTransitionCache<DataType> pair_lazy_cache;



            // Check if neighbor is on-process or off-process
            bool is_on_process = (neighbor_morton >= local_morton_start &&
                                 neighbor_morton <= local_morton_end);

            // ===== DIAGONAL CASE =====
            if (neighbor_morton == B1.morton_index && !(is_symmetric || is_hermitian)) {
                B1.schur_complement = assemble(b1_children, b1_children);
                if (B1.schur_complement.rows != B1.num_points) {
                    throw std::runtime_error("build_parent_level: parent " + std::to_string(B1.morton_index) +
                                             " has " + std::to_string(B1.num_points) + " points, its children " +
                                             std::to_string(B1.schur_complement.rows) + " skeleton points");
                }
                continue;
            }
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
                // B1's two views of B2, A_SN[B2] = K(B1, B2) and A_NS[B2] =
                // K(B2, B1), at the children's skeletons; with B2 on process,
                // the lower box builds both boxes' (B2's are B1's swapped).
                // A child pair's block comes from the views of a child on
                // process (or a ghost), else the kernel.
                if (is_on_process && static_cast<int64_t>(neighbor_morton) < B1.morton_index) {
                    continue;
                }
                std::vector<ChildRef> b2_children;
                for (int c = 0; c < num_children; ++c) {
                    const int64_t child_morton = static_cast<int64_t>(neighbor_morton) * num_children + c;
                    if (is_on_process) {
                        const int64_t b2_idx = static_cast<int64_t>(neighbor_morton) - local_morton_start;
                        b2_children.push_back(
                            child_ref(&child_level.local_boxes[static_cast<size_t>(b2_idx * num_children + c)]));
                        continue;
                    }
                    auto ghost_it = child_level.ghost_id_to_index.find(child_morton);
                    if (ghost_it != child_level.ghost_id_to_index.end()) {
                        b2_children.push_back(child_ref(&child_level.ghost_boxes[ghost_it->second]));
                        continue;
                    }
                    auto assist_it = child_level.assisting_box_points_for_kernel_evaluation.find(child_morton);
                    if (assist_it == child_level.assisting_box_points_for_kernel_evaluation.end()) {
                        throw std::runtime_error(
                            "build_parent_level: Child " + std::to_string(child_morton) +
                            " not found in ghost_boxes or assisting_boxes");
                    }
                    const auto& assist = child_level.assisting_boxes[assist_it->second];
                    ChildRef ref;
                    ref.morton = child_morton;
                    ref.skeleton = extract_skeleton_indices(assist, dimension);
                    ref.full_size = static_cast<int64_t>(assist.indices.size());
                    ref.skeleton_positions = &assist.skel_indices;
                    b2_children.push_back(std::move(ref));
                }

                MatrixStorage<DataType> I_B1_B2 = assemble(b1_children, b2_children);
                MatrixStorage<DataType> I_B2_B1 = assemble(b2_children, b1_children);
                if (I_B1_B2.rows != B1.num_points) {
                    throw std::runtime_error("build_parent_level: parent " + std::to_string(B1.morton_index) +
                                             " has " + std::to_string(B1.num_points) + " points, its children " +
                                             std::to_string(I_B1_B2.rows) + " skeleton points");
                }
                if (is_on_process) {
                    const size_t b2_idx = static_cast<size_t>(static_cast<int64_t>(neighbor_morton) - local_morton_start);
                    DeferredParentBlock sn{b2_idx, ModifiedBlock<DataType>{}, true};
                    sn.block.neighbor_morton = B1.morton_index;
                    sn.block.A_SN.set_owned(I_B2_B1.rows, I_B2_B1.cols, std::vector<DataType>(I_B2_B1.data),
                                            MatrixStorage<DataType>::FULL);
                    DeferredParentBlock ns{b2_idx, ModifiedBlock<DataType>{}, false};
                    ns.block.neighbor_morton = B1.morton_index;
                    ns.block.A_NS.set_owned(I_B1_B2.rows, I_B1_B2.cols, std::vector<DataType>(I_B1_B2.data),
                                            MatrixStorage<DataType>::FULL);
                    deferred_parent_blocks[static_cast<size_t>(occupied_slot)].push_back(std::move(sn));
                    deferred_parent_blocks[static_cast<size_t>(occupied_slot)].push_back(std::move(ns));
                }
                upsert_directed_block(B1, static_cast<int64_t>(neighbor_morton), EdgeKind::Near, true).A_SN =
                    std::move(I_B1_B2);
                upsert_directed_block(B1, static_cast<int64_t>(neighbor_morton), EdgeKind::Near, false).A_NS =
                    std::move(I_B2_B1);
            }
        }
        entry_item_end();
        } catch (const H2EntryMiss&) {
            std::lock_guard<std::mutex> lock(abandoned_mutex);
            abandoned.push_back(occupied_slot);
        } catch (...) {
            if (entry_item_abandoned()) {
                std::lock_guard<std::mutex> lock(abandoned_mutex);
                abandoned.push_back(occupied_slot);
                continue;
            }
            if (!build_failed.exchange(true, std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(build_exception_mutex);
                build_exception = std::current_exception();
            }
        }
    }
    if (build_exception) std::rethrow_exception(build_exception);
    return abandoned;
    });

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
                (deferred.a_sn ? target.near_field_interaction_map_nonsymmetry
                               : target.near_field_interaction_map)[source_morton] = block_index;
            }
        }
    }

    return parent_boxes;
}

}  // namespace color_unstructured
}  // namespace butterfly
