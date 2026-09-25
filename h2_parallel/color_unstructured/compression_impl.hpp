#pragma once

#include "occupied_topology.hpp"

namespace butterfly {
namespace color_unstructured {
using namespace fmm;

template<typename CoordType, typename DataType>
std::vector<int64_t> h2_interaction_list_unstructured(
    const ParallelTree<CoordType, DataType>* tree,
    const BoxData<CoordType, DataType>& box,
    const OccupiedTopology& occupied_topology) {
    if (box.level < 2 ||
        !occupied_topology.contains(box.level, box.morton_index)) {
        return {};
    }

    const int children_per_parent =
        morton::children_per_box(tree->dimension);
    const uint32_t parent_grid_size = 1u << (box.level - 1);
    const int64_t parent_morton = box.morton_index / children_per_parent;
    std::vector<uint64_t> parent_neighbors = morton::neighbors_nd(
        tree->dimension, parent_morton, parent_grid_size);
    parent_neighbors.push_back(static_cast<uint64_t>(parent_morton));

    std::unordered_set<int64_t> near;
    near.reserve(box.one_hop.size() + 1);
    near.insert(box.morton_index);
    near.insert(box.one_hop.begin(), box.one_hop.end());

    std::vector<int64_t> interaction;
    interaction.reserve(parent_neighbors.size() * children_per_parent);
    const int64_t boxes_at_level = 1LL << (tree->dimension * box.level);
    for (uint64_t parent : parent_neighbors) {
        const int64_t first_child =
            static_cast<int64_t>(parent) * children_per_parent;
        for (int child = 0; child < children_per_parent; ++child) {
            const int64_t candidate = first_child + child;
            if (candidate >= 0 && candidate < boxes_at_level &&
                occupied_topology.contains(box.level, candidate) &&
                !near.count(candidate)) {
                interaction.push_back(candidate);
            }
        }
    }

    std::sort(interaction.begin(), interaction.end());
    interaction.erase(
        std::unique(interaction.begin(), interaction.end()), interaction.end());
    return interaction;
}

template<typename CoordType, typename DataType>
void exchange_h2_point_metadata_unstructured(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    bool need_id_points,
    bool need_interaction_points,
    bool need_near_points,
    const OccupiedTopology& occupied_topology) {
    auto& level = tree->levels[level_number];
    if (!level.is_process_active) return;

    std::vector<int64_t> needed;
    const auto occupied_indices = occupied_local_indices(level);
    for (int64_t local_index : occupied_indices) {
        const auto& box = level.local_boxes[static_cast<size_t>(local_index)];
        if (need_id_points) {
            needed.insert(needed.end(), box.two_hop.begin(), box.two_hop.end());
        }
        if (need_interaction_points) {
            auto interaction = h2_interaction_list_unstructured(
                tree, box, occupied_topology);
            needed.insert(needed.end(), interaction.begin(), interaction.end());
        }
        if (need_near_points) {
            needed.insert(needed.end(), box.one_hop.begin(), box.one_hop.end());
        }
    }

    needed.erase(
        std::remove_if(
            needed.begin(), needed.end(),
            [&](int64_t morton_index) {
                return !occupied_topology.contains(level_number, morton_index) ||
                    level.find_local_box(morton_index) != nullptr;
            }),
        needed.end());
    std::sort(needed.begin(), needed.end());
    needed.erase(std::unique(needed.begin(), needed.end()), needed.end());

    ensure_h2_assisting_slots(level, needed);
    const std::vector<int> neighbor_ranks =
        compute_one_hop_neighbor_ranks(tree, level, level_number);
    exchange_assisting_for_mortons_onehop(
        tree, level, level_number, neighbor_ranks, needed);
}

template<typename CoordType, typename DataType, typename KernelType>
void build_h2_blocks_for_level_unstructured(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    KernelType* kernel,
    bool build_interactions,
    bool build_leaf_near,
    const OccupiedTopology& occupied_topology) {
    auto& level = tree->levels[level_number];
    if (!level.is_process_active) return;

    std::exception_ptr build_exception;
    std::mutex build_exception_mutex;
    std::atomic<bool> build_failed{false};

    const auto occupied_indices = occupied_local_indices(level);
    #pragma omp parallel for schedule(dynamic) if (occupied_indices.size() > 1)
    for (int64_t occupied_slot = 0;
         occupied_slot < static_cast<int64_t>(occupied_indices.size());
         ++occupied_slot) {
        if (build_failed.load(std::memory_order_relaxed)) continue;
        try {
            const int64_t box_index =
                occupied_indices[static_cast<size_t>(occupied_slot)];
            auto& target = level.local_boxes[static_cast<size_t>(box_index)];
            target.h2_interaction_blocks.clear();
            target.h2_near_blocks.clear();

            if (build_interactions) {
                const auto target_indices =
                    h2_box_global_indices(level, target.morton_index, true);
                const auto interaction = h2_interaction_list_unstructured(
                    tree, target, occupied_topology);
                target.h2_interaction_blocks.reserve(interaction.size());
                for (int64_t source_morton : interaction) {
                    const auto source_indices =
                        h2_box_global_indices(level, source_morton, true);
                    if (target_indices.empty() || source_indices.empty()) continue;

                    H2Block<DataType> block;
                    block.source_morton = source_morton;
                    block.matrix.allocate(
                        static_cast<int64_t>(target_indices.size()),
                        static_cast<int64_t>(source_indices.size()),
                        MatrixStorage<DataType>::FULL);
                    kernel->evaluate_block_by_index(
                        target_indices.data(),
                        static_cast<int64_t>(target_indices.size()),
                        source_indices.data(),
                        static_cast<int64_t>(source_indices.size()),
                        block.matrix.data.data(), block.matrix.lda);
                    target.h2_interaction_blocks.push_back(std::move(block));
                }
            }

            if (build_leaf_near) {
                std::vector<int64_t> near = target.one_hop;
                near.push_back(target.morton_index);
                std::sort(near.begin(), near.end());
                near.erase(std::unique(near.begin(), near.end()), near.end());
                target.h2_near_blocks.reserve(near.size());
                for (int64_t source_morton : near) {
                    const auto source_indices =
                        h2_box_global_indices(level, source_morton, false);
                    if (target.point_indices.empty() || source_indices.empty()) {
                        continue;
                    }

                    H2Block<DataType> block;
                    block.source_morton = source_morton;
                    block.matrix.allocate(
                        target.num_points,
                        static_cast<int64_t>(source_indices.size()),
                        MatrixStorage<DataType>::FULL);
                    kernel->evaluate_block_by_index(
                        target.point_indices.data(), target.num_points,
                        source_indices.data(),
                        static_cast<int64_t>(source_indices.size()),
                        block.matrix.data.data(), block.matrix.lda);
                    target.h2_near_blocks.push_back(std::move(block));
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
}

template<typename CoordType, typename DataType>
std::vector<BoxData<CoordType, DataType>> build_h2_parent_boxes_unstructured(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    int dimension,
    const CoordType global_bounds[6],
    const OccupiedTopology& occupied_topology) {
    if (!child_level.is_process_active || child_level.local_boxes.empty()) {
        return {};
    }

    const int children_per_parent = morton::children_per_box(dimension);
    if (child_level.local_boxes.size() % children_per_parent != 0 ||
        child_level.local_boxes.front().morton_index % children_per_parent != 0) {
        throw std::runtime_error(
            "build_h2_parent_boxes_unstructured: child slab is not parent aligned");
    }

    const int parent_level_index = child_level.level - 1;
    const int64_t parent_count =
        static_cast<int64_t>(child_level.local_boxes.size()) /
        children_per_parent;
    const int64_t parent_start =
        child_level.local_boxes.front().morton_index / children_per_parent;
    std::vector<BoxData<CoordType, DataType>> parents(
        static_cast<size_t>(parent_count));
    initialize_local_boxes(
        parents, parent_start, parent_level_index, dimension, global_bounds);
    compute_occupied_neighbor_lists(
        parents, parent_level_index, dimension, occupied_topology);

    const auto occupied_parent_indices = occupied_box_indices_from_topology(
        parents, parent_level_index, occupied_topology);
    for (int64_t parent_index : occupied_parent_indices) {
        auto& parent = parents[static_cast<size_t>(parent_index)];
        const int64_t first_child = parent_index * children_per_parent;

        int64_t point_count = 0;
        for (int child = 0; child < children_per_parent; ++child) {
            point_count += static_cast<int64_t>(
                child_level.local_boxes[
                    static_cast<size_t>(first_child + child)]
                    .skeleton_indices.size());
        }
        parent.point_indices.reserve(static_cast<size_t>(point_count));
        parent.point_coords.reserve(
            static_cast<size_t>(point_count * dimension));

        for (int child = 0; child < children_per_parent; ++child) {
            auto& child_box = child_level.local_boxes[
                static_cast<size_t>(first_child + child)];
            parent.children_morton[child] = child_box.morton_index;
            for (int64_t skeleton_index : child_box.skeleton_indices) {
                parent.point_indices.push_back(
                    child_box.point_indices.at(
                        static_cast<size_t>(skeleton_index)));
                for (int d = 0; d < dimension; ++d) {
                    parent.point_coords.push_back(
                        child_box.point_coords.at(static_cast<size_t>(
                            skeleton_index * dimension + d)));
                }
            }
        }
        parent.num_children = children_per_parent;
        parent.num_points = point_count;

        if (auto* existing =
                parent_level.find_local_box(parent.morton_index)) {
            parent.on_boundary = existing->on_boundary;
        }
    }
    return parents;
}

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_compression_unstructured(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    bool use_sketch,
    bool verbose,
    const OccupiedTopology& occupied_topology) {

    const int rank = tree->mpi_rank;
    const int leaf_level = tree->num_levels - 1;
    int64_t local_rankmax = 0;

    if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
        std::cout << "\n========================================\n"
                  << "H2 Compression Only (Parallel MPI)\n"
                  << "========================================" << std::endl;
    }

    if (leaf_level < 2) {
        exchange_h2_point_metadata_unstructured(
            tree, leaf_level, false, false, true, occupied_topology);
        kernel->register_level_coordinates(tree->levels[leaf_level]);
        build_h2_blocks_for_level_unstructured(
            tree, leaf_level, kernel, false, true, occupied_topology);
    } else {
        for (int level_number = leaf_level; level_number >= 2; --level_number) {
            auto& level = tree->levels[level_number];
            filter_level_topology(tree, level_number, occupied_topology);
            const int print_rank = smallest_active_rank(level);
            const double level_start = MPI_Wtime();

            if (level.is_process_active && !level.eliminated_boxes.empty()) {
                throw std::runtime_error(
                    "hierarchical_compression_parallel: tree already contains elimination state");
            }

            exchange_h2_point_metadata_unstructured(
                tree, level_number, true, true, level_number == leaf_level,
                occupied_topology);
            kernel->register_level_coordinates(level);

            if (level.is_process_active) {
                std::exception_ptr id_exception;
                std::mutex id_exception_mutex;
                std::atomic<bool> id_failed{false};

                const auto occupied_indices = occupied_local_indices(level);
                #pragma omp parallel for schedule(dynamic) if (occupied_indices.size() > 1)
                for (int64_t occupied_slot = 0;
                     occupied_slot < static_cast<int64_t>(occupied_indices.size());
                     ++occupied_slot) {
                    if (id_failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const int64_t box_index =
                            occupied_indices[static_cast<size_t>(occupied_slot)];
                        auto& box = level.local_boxes[static_cast<size_t>(box_index)];
                        h2_skeletonize_box(
                            tree,
                            &box,
                            level, kernel, tolerance, use_sketch);
                    } catch (...) {
                        if (!id_failed.exchange(true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(id_exception_mutex);
                            id_exception = std::current_exception();
                        }
                    }
                }
                if (id_exception) std::rethrow_exception(id_exception);

                for (int64_t local_index : occupied_indices) {
                    const auto& box =
                        level.local_boxes[static_cast<size_t>(local_index)];
                    local_rankmax = std::max<int64_t>(
                        local_rankmax, static_cast<int64_t>(
                            box.skeleton_indices.size()));
                }
            }

            // Refresh remote records after all owners have selected skeletons.
            exchange_h2_point_metadata_unstructured(
                tree, level_number, true, true, level_number == leaf_level,
                occupied_topology);
            kernel->register_level_coordinates(level);
            build_h2_blocks_for_level_unstructured(
                tree, level_number, kernel, true,
                level_number == leaf_level, occupied_topology);

            if (level_number > 2) {
                auto parents = build_h2_parent_boxes_unstructured(
                    level, tree->levels[level_number - 1],
                    tree->dimension, tree->global_bounds,
                    occupied_topology);
                install_h2_parent_boxes(tree, level_number, std::move(parents));
            }

            if (verbose) {
                int64_t local_skeletons = 0;
                int64_t local_points = 0;
                for (int64_t local_index : occupied_local_indices(level)) {
                    const auto& box =
                        level.local_boxes[static_cast<size_t>(local_index)];
                    local_skeletons += static_cast<int64_t>(
                        box.skeleton_indices.size());
                    local_points += box.num_points;
                }

                int64_t global_skeletons = 0;
                int64_t global_points = 0;
                double level_elapsed = MPI_Wtime() - level_start;
                double max_level_elapsed = 0.0;
                MPI_Reduce(
                    &local_skeletons, &global_skeletons, 1, MPI_INT64_T,
                    MPI_SUM, print_rank, tree->comm);
                MPI_Reduce(
                    &local_points, &global_points, 1, MPI_INT64_T,
                    MPI_SUM, print_rank, tree->comm);
                MPI_Reduce(
                    &level_elapsed, &max_level_elapsed, 1, MPI_DOUBLE,
                    MPI_MAX, print_rank, tree->comm);

                if (rank == print_rank) {
                    const double ratio = global_points > 0
                        ? static_cast<double>(global_skeletons) /
                            static_cast<double>(global_points)
                        : 0.0;
                    std::cout << "  Level " << level_number
                              << ": compression ratio=" << ratio
                              << ", time=" << max_level_elapsed << " s"
                              << std::endl;
                }
            }
        }
    }

    size_t local_memory = 0;
    for (int level_number = 0; level_number <= leaf_level; ++level_number) {
        const auto& level = tree->levels[level_number];
        for (const auto& box : level.local_boxes) {
            local_memory += calculate_box_data_size(box);
        }
    }
    if (out_rankmax) *out_rankmax = local_rankmax;
    if (memory_per_rank) *memory_per_rank = local_memory;

}


}  // namespace color_unstructured
}  // namespace butterfly
