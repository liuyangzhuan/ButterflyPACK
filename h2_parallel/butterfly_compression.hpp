#pragma once

#include "butterfly_types.hpp"
#include "butterfly_solve.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace butterfly {
using namespace fmm;

template<typename CoordType, typename DataType>
std::vector<int64_t> h2_interaction_list(
    const ParallelTree<CoordType, DataType>* tree,
    const BoxData<CoordType, DataType>& box) {

    if (box.level < 2) return {};

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
        const int64_t first_child = static_cast<int64_t>(parent) * children_per_parent;
        for (int child = 0; child < children_per_parent; ++child) {
            const int64_t candidate = first_child + child;
            if (candidate >= 0 && candidate < boxes_at_level && !near.count(candidate)) {
                interaction.push_back(candidate);
            }
        }
    }

    std::sort(interaction.begin(), interaction.end());
    interaction.erase(std::unique(interaction.begin(), interaction.end()), interaction.end());
    return interaction;
}

template<typename CoordType, typename DataType>
void ensure_h2_assisting_slots(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<int64_t>& mortons) {

    for (int64_t morton_index : mortons) {
        if (level.find_local_box(morton_index) != nullptr) continue;
        if (level.assisting_box_points_for_kernel_evaluation.count(morton_index)) continue;

        const int64_t slot = static_cast<int64_t>(level.assisting_boxes.size());
        PointDataRequest<CoordType> request;
        request.morton_index = morton_index;
        level.assisting_boxes.push_back(std::move(request));
        level.assisting_box_points_for_kernel_evaluation[morton_index] = slot;
    }
}

template<typename CoordType, typename DataType>
std::vector<int64_t> h2_box_global_indices(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton_index,
    bool skeleton_only) {

    if (auto* box = level.find_local_box(morton_index)) {
        if (!skeleton_only) return box->point_indices;

        std::vector<int64_t> result;
        result.reserve(box->skeleton_indices.size());
        for (int64_t local_index : box->skeleton_indices) {
            if (local_index < 0 || local_index >= static_cast<int64_t>(box->point_indices.size())) {
                throw std::runtime_error("h2_box_global_indices: invalid local skeleton index");
            }
            result.push_back(box->point_indices[static_cast<size_t>(local_index)]);
        }
        return result;
    }

    auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(morton_index);
    if (assist_it == level.assisting_box_points_for_kernel_evaluation.end()) {
        throw std::runtime_error(
            "h2_box_global_indices: missing remote box " + std::to_string(morton_index));
    }
    const auto& request = level.assisting_boxes.at(static_cast<size_t>(assist_it->second));
    if (!skeleton_only) {
        if (request.indices.empty()) {
            throw std::runtime_error(
                "h2_box_global_indices: remote box has no point indices " +
                std::to_string(morton_index));
        }
        return request.indices;
    }

    std::vector<int64_t> result;
    result.reserve(request.skel_indices.size());
    for (int64_t local_index : request.skel_indices) {
        if (local_index < 0 || local_index >= static_cast<int64_t>(request.indices.size())) {
            throw std::runtime_error("h2_box_global_indices: invalid remote skeleton index");
        }
        result.push_back(request.indices[static_cast<size_t>(local_index)]);
    }
    return result;
}

template<typename CoordType, typename DataType>
void exchange_h2_point_metadata(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    bool need_id_points,
    bool need_interaction_points,
    bool need_near_points) {

    auto& level = tree->levels[level_number];
    if (!level.is_process_active) return;

    std::vector<int64_t> needed;
    for (const auto& box : level.local_boxes) {
        if (need_id_points) {
            needed.insert(needed.end(), box.two_hop.begin(), box.two_hop.end());
        }
        if (need_interaction_points) {
            auto interaction = h2_interaction_list(tree, box);
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
                return level.find_local_box(morton_index) != nullptr;
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
void h2_skeletonize_box(
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    double tolerance,
    bool use_sketch) {

    FactorizationThreadScratch<CoordType, DataType> scratch;
    gather_id_workspace(
        box, level, kernel,
        static_cast<const CoordType*>(nullptr), 0, CoordType{0}, true,
        scratch.workspace, scratch.workspace_rows, scratch.workspace_cols,
        0, box->on_boundary);

    if (scratch.workspace_cols == 0) {
        box->skeleton_indices.clear();
        box->redundant_indices.clear();
        box->interpolation_matrix.allocate(0, 0, MatrixStorage<DataType>::NONE);
        return;
    }
    if (scratch.workspace_rows == 0) {
        box->skeleton_indices.resize(static_cast<size_t>(box->num_points));
        std::iota(box->skeleton_indices.begin(), box->skeleton_indices.end(), int64_t{0});
        box->redundant_indices.clear();
        box->interpolation_matrix.allocate(box->num_points, 0);
        return;
    }

    constexpr double sketch_factor = 1.0;
    constexpr int sketch_nonzeros = 4;
    IDResult<DataType> id;
    if (use_sketch) {
        id = compute_id_sparse_sketch(
            scratch.workspace.data(), scratch.sketch_storage,
            scratch.workspace_rows, scratch.workspace_cols, scratch.workspace_rows,
            tolerance, sketch_factor, sketch_nonzeros,
            static_cast<uint64_t>(box->morton_index + 1));
    } else {
        id = compute_id_complex(
            scratch.workspace.data(),
            scratch.workspace_rows, scratch.workspace_cols, scratch.workspace_rows,
            tolerance, 0);
    }

    box->skeleton_indices = std::move(id.skeleton_indices);
    box->redundant_indices = std::move(id.redundant_indices);
    box->interpolation_matrix = std::move(id.interpolation);
}

template<typename CoordType, typename DataType, typename KernelType>
void build_h2_blocks_for_level(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    KernelType* kernel,
    bool build_interactions,
    bool build_leaf_near) {

    auto& level = tree->levels[level_number];
    if (!level.is_process_active) return;

    std::exception_ptr build_exception;
    std::mutex build_exception_mutex;
    std::atomic<bool> build_failed{false};

    #pragma omp parallel for schedule(dynamic) if (level.local_boxes.size() > 1)
    for (int64_t box_index = 0;
         box_index < static_cast<int64_t>(level.local_boxes.size());
         ++box_index) {
        if (build_failed.load(std::memory_order_relaxed)) continue;

        try {
            auto& target = level.local_boxes[static_cast<size_t>(box_index)];
            target.h2_interaction_blocks.clear();
            target.h2_near_blocks.clear();

            if (build_interactions) {
                const auto target_indices =
                    h2_box_global_indices(level, target.morton_index, true);
                const auto interaction = h2_interaction_list(tree, target);
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
                        target_indices.data(), static_cast<int64_t>(target_indices.size()),
                        source_indices.data(), static_cast<int64_t>(source_indices.size()),
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
                    if (target.point_indices.empty() || source_indices.empty()) continue;

                    H2Block<DataType> block;
                    block.source_morton = source_morton;
                    block.matrix.allocate(
                        target.num_points,
                        static_cast<int64_t>(source_indices.size()),
                        MatrixStorage<DataType>::FULL);
                    kernel->evaluate_block_by_index(
                        target.point_indices.data(), target.num_points,
                        source_indices.data(), static_cast<int64_t>(source_indices.size()),
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
std::vector<BoxData<CoordType, DataType>> build_h2_parent_boxes(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    int dimension,
    const CoordType global_bounds[6]) {

    if (!child_level.is_process_active || child_level.local_boxes.empty()) return {};

    const int children_per_parent = morton::children_per_box(dimension);
    if (child_level.local_boxes.size() % children_per_parent != 0 ||
        child_level.local_boxes.front().morton_index % children_per_parent != 0) {
        throw std::runtime_error("build_h2_parent_boxes: child slab is not parent aligned");
    }

    const int64_t parent_count =
        static_cast<int64_t>(child_level.local_boxes.size()) / children_per_parent;
    const int64_t parent_start =
        child_level.local_boxes.front().morton_index / children_per_parent;
    std::vector<BoxData<CoordType, DataType>> parents(static_cast<size_t>(parent_count));
    initialize_local_boxes(
        parents, parent_start, child_level.level - 1, dimension, global_bounds);
    compute_neighbor_lists(parents, dimension, child_level.level - 1);

    for (int64_t parent_index = 0; parent_index < parent_count; ++parent_index) {
        auto& parent = parents[static_cast<size_t>(parent_index)];
        const int64_t first_child = parent_index * children_per_parent;

        int64_t point_count = 0;
        for (int child = 0; child < children_per_parent; ++child) {
            point_count += static_cast<int64_t>(
                child_level.local_boxes[static_cast<size_t>(first_child + child)]
                    .skeleton_indices.size());
        }
        parent.point_indices.reserve(static_cast<size_t>(point_count));
        parent.point_coords.reserve(static_cast<size_t>(point_count * dimension));

        for (int child = 0; child < children_per_parent; ++child) {
            auto& child_box =
                child_level.local_boxes[static_cast<size_t>(first_child + child)];
            parent.children_morton[child] = child_box.morton_index;
            for (int64_t skeleton_index : child_box.skeleton_indices) {
                parent.point_indices.push_back(
                    child_box.point_indices.at(static_cast<size_t>(skeleton_index)));
                for (int d = 0; d < dimension; ++d) {
                    parent.point_coords.push_back(
                        child_box.point_coords.at(
                            static_cast<size_t>(skeleton_index * dimension + d)));
                }
            }
        }
        parent.num_children = children_per_parent;
        parent.num_points = point_count;

        if (auto* existing = parent_level.find_local_box(parent.morton_index)) {
            parent.on_boundary = existing->on_boundary;
        }
    }
    return parents;
}

template<typename CoordType, typename DataType>
void install_h2_parent_boxes(
    ParallelTree<CoordType, DataType>* tree,
    int child_level_number,
    std::vector<BoxData<CoordType, DataType>> local_parents) {

    auto& child_level = tree->levels[child_level_number];
    auto& parent_level = tree->levels[child_level_number - 1];
    const int rank = tree->mpi_rank;
    const bool reduction =
        parent_level.num_active_processes != child_level.num_active_processes;

    if (!reduction) {
        if (parent_level.is_process_active) {
            parent_level.local_boxes = std::move(local_parents);
        }
        return;
    }

    const bool keep_local =
        child_level.is_process_active && child_level.parent_level_owner == rank;
    const bool send_parents =
        child_level.is_process_active && child_level.parent_level_owner != rank;
    const int size_tag = 7200 + 2 * child_level_number;
    const int data_tag = size_tag + 1;

    if (parent_level.is_process_active) {
        std::vector<BoxData<CoordType, DataType>> gathered;
        for (int child_rank : parent_level.children_senders) {
            if (keep_local && child_rank == rank) {
                gathered.insert(
                    gathered.end(),
                    std::make_move_iterator(local_parents.begin()),
                    std::make_move_iterator(local_parents.end()));
                continue;
            }

            int64_t buffer_size = 0;
            MPI_Status status;
            MPI_Recv(
                &buffer_size, 1, MPI_INT64_T, child_rank, size_tag, tree->comm, &status);
            std::vector<char> buffer(static_cast<size_t>(buffer_size));
            MPI_Recv_large(
                buffer.data(), static_cast<size_t>(buffer_size), MPI_CHAR,
                child_rank, data_tag, tree->comm, &status);
            auto received = deserialize_boxes<CoordType, DataType>(buffer);
            gathered.insert(
                gathered.end(),
                std::make_move_iterator(received.begin()),
                std::make_move_iterator(received.end()));
        }
        std::sort(
            gathered.begin(), gathered.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.morton_index < rhs.morton_index;
            });
        parent_level.local_boxes = std::move(gathered);
    }

    if (send_parents) {
        std::vector<char> buffer = serialize_boxes(local_parents);
        const int64_t buffer_size = static_cast<int64_t>(buffer.size());
        MPI_Send(
            &buffer_size, 1, MPI_INT64_T,
            child_level.parent_level_owner, size_tag, tree->comm);
        MPI_Send_large(
            buffer.data(), buffer.size(), MPI_CHAR,
            child_level.parent_level_owner, data_tag, tree->comm);
    }

    if (parent_level.is_process_active &&
        parent_level.local_boxes.size() != static_cast<size_t>(parent_level.num_boxes_local)) {
        throw std::runtime_error(
            "install_h2_parent_boxes: received the wrong number of parent boxes");
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_compression_parallel(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    bool use_sketch,
    bool verbose = true) {

    const int rank = tree->mpi_rank;
    const int leaf_level = tree->num_levels - 1;
    int64_t local_rankmax = 0;

    if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
        std::cout << "\n========================================\n"
                  << "H2 Compression Only (Parallel MPI)\n"
                  << "========================================" << std::endl;
    }

    if (leaf_level < 2) {
        exchange_h2_point_metadata(tree, leaf_level, false, false, true);
        build_h2_blocks_for_level(tree, leaf_level, kernel, false, true);
    } else {
        for (int level_number = leaf_level; level_number >= 2; --level_number) {
            auto& level = tree->levels[level_number];
            const int print_rank = smallest_active_rank(level);
            const double level_start = MPI_Wtime();

            if (level.is_process_active && !level.eliminated_boxes.empty()) {
                throw std::runtime_error(
                    "hierarchical_compression_parallel: tree already contains elimination state");
            }

            exchange_h2_point_metadata(
                tree, level_number, true, true, level_number == leaf_level);

            if (level.is_process_active) {
                std::exception_ptr id_exception;
                std::mutex id_exception_mutex;
                std::atomic<bool> id_failed{false};

                #pragma omp parallel for schedule(dynamic) if (level.local_boxes.size() > 1)
                for (int64_t box_index = 0;
                     box_index < static_cast<int64_t>(level.local_boxes.size());
                     ++box_index) {
                    if (id_failed.load(std::memory_order_relaxed)) continue;
                    try {
                        h2_skeletonize_box(
                            &level.local_boxes[static_cast<size_t>(box_index)],
                            level, kernel, tolerance, use_sketch);
                    } catch (...) {
                        if (!id_failed.exchange(true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(id_exception_mutex);
                            id_exception = std::current_exception();
                        }
                    }
                }
                if (id_exception) std::rethrow_exception(id_exception);

                for (const auto& box : level.local_boxes) {
                    local_rankmax = std::max<int64_t>(
                        local_rankmax, static_cast<int64_t>(box.skeleton_indices.size()));
                }
            }

            // Refresh remote records after all owners have selected skeletons.
            exchange_h2_point_metadata(
                tree, level_number, true, true, level_number == leaf_level);
            build_h2_blocks_for_level(
                tree, level_number, kernel, true, level_number == leaf_level);

            if (level_number > 2) {
                auto parents = build_h2_parent_boxes(
                    level, tree->levels[level_number - 1],
                    tree->dimension, tree->global_bounds);
                install_h2_parent_boxes(tree, level_number, std::move(parents));
            }

            if (verbose) {
                int64_t local_skeletons = 0;
                int64_t local_points = 0;
                for (const auto& box : level.local_boxes) {
                    local_skeletons += static_cast<int64_t>(box.skeleton_indices.size());
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

    if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
        std::cout << "H2 compression-only construction complete" << std::endl;
    }
}

template<typename CoordType, typename DataType>
std::vector<DataType> h2_local_vector_for_morton(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<SolveDataRequest<CoordType, DataType>>& level_data,
    int64_t morton_index,
    bool skeleton_only) {

    auto* box = level.find_local_box(morton_index);
    if (box == nullptr) {
        throw std::runtime_error(
            "h2_local_vector_for_morton: requested box is not local");
    }
    const int64_t local_index = morton_index - level.local_morton_start;
    if (local_index < 0 || local_index >= static_cast<int64_t>(level_data.size())) {
        throw std::runtime_error(
            "h2_local_vector_for_morton: missing local vector data");
    }
    const auto& data = level_data[static_cast<size_t>(local_index)];

    if (!skeleton_only) return data.right_side;

    const int64_t skeleton_count =
        static_cast<int64_t>(box->skeleton_indices.size());
    std::vector<DataType> result(
        static_cast<size_t>(skeleton_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            result[static_cast<size_t>(i + column * skeleton_count)] =
                data.left_side.at(static_cast<size_t>(
                    box->skeleton_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    return result;
}

template<typename CoordType, typename DataType>
std::unordered_map<int64_t, std::vector<DataType>> exchange_h2_vectors_onehop(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    const std::vector<SolveDataRequest<CoordType, DataType>>& level_data,
    std::vector<int64_t> needed_remote_mortons,
    bool skeleton_only,
    int tag_base) {

    auto& level = tree->levels[level_number];
    std::unordered_map<int64_t, std::vector<DataType>> received_vectors;
    if (!level.is_process_active) return received_vectors;

    const int rank = tree->mpi_rank;
    const auto neighbor_ranks =
        compute_one_hop_neighbor_ranks(tree, level, level_number);
    const std::unordered_set<int> neighbor_set(
        neighbor_ranks.begin(), neighbor_ranks.end());

    const uint32_t grid_size = 1u << level_number;
    auto owner_of_morton = [&](int64_t morton_index) {
        std::vector<uint64_t> one{static_cast<uint64_t>(morton_index)};
        const auto regions = morton::assign_to_processes_nd(
            tree->dimension, one, level.num_active_processes, grid_size);
        return level.morton_to_rank.at(static_cast<int>(regions.front()));
    };

    std::sort(needed_remote_mortons.begin(), needed_remote_mortons.end());
    needed_remote_mortons.erase(
        std::unique(needed_remote_mortons.begin(), needed_remote_mortons.end()),
        needed_remote_mortons.end());

    std::unordered_map<int, std::vector<int64_t>> requests_to_send;
    for (int64_t morton_index : needed_remote_mortons) {
        const int owner = owner_of_morton(morton_index);
        if (owner == rank) continue;
        if (!neighbor_set.count(owner)) {
            throw std::runtime_error(
                "exchange_h2_vectors_onehop: source owner is not a neighboring rank");
        }
        requests_to_send[owner].push_back(morton_index);
    }

    std::vector<int> send_counts(neighbor_ranks.size(), 0);
    std::vector<int> recv_counts(neighbor_ranks.size(), 0);
    std::vector<MPI_Request> mpi_requests;
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        auto it = requests_to_send.find(neighbor_ranks[i]);
        if (it != requests_to_send.end()) {
            send_counts[i] = static_cast<int>(it->second.size());
        }
        MPI_Request request;
        MPI_Irecv(
            &recv_counts[i], 1, MPI_INT, neighbor_ranks[i],
            tag_base, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Send(
            &send_counts[i], 1, MPI_INT, neighbor_ranks[i],
            tag_base, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<int64_t>> requests_received(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_counts[i] <= 0) continue;
        requests_received[i].resize(static_cast<size_t>(recv_counts[i]));
        MPI_Request request;
        MPI_Irecv(
            requests_received[i].data(), recv_counts[i], MPI_INT64_T,
            neighbor_ranks[i], tag_base + 1, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_counts[i] <= 0) continue;
        MPI_Send(
            requests_to_send[neighbor_ranks[i]].data(), send_counts[i], MPI_INT64_T,
            neighbor_ranks[i], tag_base + 1, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<char>> send_buffers(neighbor_ranks.size());
    std::vector<uint64_t> send_sizes(neighbor_ranks.size(), 0);
    std::vector<uint64_t> recv_sizes(neighbor_ranks.size(), 0);
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        size_t total_size = 0;
        std::vector<std::vector<DataType>> values;
        values.reserve(requests_received[i].size());
        for (int64_t morton_index : requests_received[i]) {
            values.push_back(h2_local_vector_for_morton(
                level, level_data, morton_index, skeleton_only));
            total_size += 2 * sizeof(int64_t) +
                values.back().size() * sizeof(DataType);
        }

        send_buffers[i].resize(total_size);
        char* pointer = send_buffers[i].data();
        for (size_t record = 0; record < requests_received[i].size(); ++record) {
            const int64_t morton_index = requests_received[i][record];
            const int64_t value_count =
                static_cast<int64_t>(values[record].size());
            std::memcpy(pointer, &morton_index, sizeof(int64_t));
            pointer += sizeof(int64_t);
            std::memcpy(pointer, &value_count, sizeof(int64_t));
            pointer += sizeof(int64_t);
            if (value_count > 0) {
                std::memcpy(
                    pointer, values[record].data(),
                    static_cast<size_t>(value_count) * sizeof(DataType));
                pointer += static_cast<size_t>(value_count) * sizeof(DataType);
            }
        }
        send_sizes[i] = static_cast<uint64_t>(total_size);
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request request;
        MPI_Irecv(
            &recv_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i],
            tag_base + 2, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Send(
            &send_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i],
            tag_base + 2, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<char>> recv_buffers(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        recv_buffers[i].resize(static_cast<size_t>(recv_sizes[i]));
        if (recv_sizes[i] == 0) continue;
        MPI_Irecv_large(
            recv_buffers[i].data(), recv_buffers[i].size(), MPI_CHAR,
            neighbor_ranks[i], tag_base + 3, tree->comm, mpi_requests);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_sizes[i] == 0) continue;
        MPI_Send_large(
            send_buffers[i].data(), send_buffers[i].size(), MPI_CHAR,
            neighbor_ranks[i], tag_base + 3, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
    }

    for (const auto& buffer : recv_buffers) {
        const char* pointer = buffer.data();
        const char* end = pointer + buffer.size();
        while (pointer < end) {
            if (end - pointer < static_cast<ptrdiff_t>(2 * sizeof(int64_t))) {
                throw std::runtime_error(
                    "exchange_h2_vectors_onehop: truncated record header");
            }
            int64_t morton_index = -1;
            int64_t value_count = -1;
            std::memcpy(&morton_index, pointer, sizeof(int64_t));
            pointer += sizeof(int64_t);
            std::memcpy(&value_count, pointer, sizeof(int64_t));
            pointer += sizeof(int64_t);
            if (value_count < 0 ||
                end - pointer < static_cast<ptrdiff_t>(
                    static_cast<size_t>(value_count) * sizeof(DataType))) {
                throw std::runtime_error(
                    "exchange_h2_vectors_onehop: invalid record size");
            }
            auto& values = received_vectors[morton_index];
            values.resize(static_cast<size_t>(value_count));
            if (value_count > 0) {
                std::memcpy(
                    values.data(), pointer,
                    static_cast<size_t>(value_count) * sizeof(DataType));
                pointer += static_cast<size_t>(value_count) * sizeof(DataType);
            }
        }
    }
    return received_vectors;
}

template<typename DataType>
void h2_matrix_vector_product(
    const MatrixStorage<DataType>& matrix,
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    int64_t nrhs,
    char transpose = 'N') {

    const bool transposed = transpose != 'N' && transpose != 'n';
    const int64_t input_size = transposed ? matrix.rows : matrix.cols;
    const int64_t output_size = transposed ? matrix.cols : matrix.rows;
    if (nrhs <= 0 || static_cast<int64_t>(input.size()) != input_size * nrhs) {
        throw std::runtime_error("h2_matrix_vector_product: input dimension mismatch");
    }
    output.assign(static_cast<size_t>(output_size * nrhs), DataType{0});
    if (input_size == 0 || output_size == 0) return;

    const int m = static_cast<int>(output_size);
    const int n = static_cast<int>(nrhs);
    const int k = static_cast<int>(input_size);
    const int lda = static_cast<int>(matrix.lda);
    const int ldb = static_cast<int>(input_size);
    const int ldc = static_cast<int>(output_size);
    const DataType alpha = DataType{1};
    const DataType beta = DataType{0};
    const char trans_b = 'N';
    gemm_(&transpose, &trans_b, &m, &n, &k, &alpha,
          matrix.data.data(), &lda, input.data(), &ldb,
          &beta, output.data(), &ldc);
}

template<typename CoordType, typename DataType>
void apply_h2_upward_projection(
    const BoxData<CoordType, DataType>& box,
    SolveDataRequest<CoordType, DataType>& data) {

    if (box.redundant_indices.empty()) return;
    const int64_t redundant_count =
        static_cast<int64_t>(box.redundant_indices.size());
    const int64_t skeleton_count =
        static_cast<int64_t>(box.skeleton_indices.size());
    std::vector<DataType> redundant(
        static_cast<size_t>(redundant_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < redundant_count; ++i) {
            redundant[static_cast<size_t>(i + column * redundant_count)] =
                data.left_side.at(static_cast<size_t>(
                    box.redundant_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    std::vector<DataType> projected;
    h2_matrix_vector_product(
        box.interpolation_matrix, redundant, projected, data.nrhs, 'N');
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            data.left_side.at(static_cast<size_t>(
                box.skeleton_indices[static_cast<size_t>(i)] +
                column * data.num_points)) +=
                projected[static_cast<size_t>(i + column * skeleton_count)];
        }
    }
}

template<typename CoordType, typename DataType>
void apply_h2_downward_interpolation(
    const BoxData<CoordType, DataType>& box,
    SolveDataRequest<CoordType, DataType>& data) {

    if (box.redundant_indices.empty()) return;
    const int64_t skeleton_count =
        static_cast<int64_t>(box.skeleton_indices.size());
    const int64_t redundant_count =
        static_cast<int64_t>(box.redundant_indices.size());
    std::vector<DataType> skeleton(
        static_cast<size_t>(skeleton_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            skeleton[static_cast<size_t>(i + column * skeleton_count)] =
                data.left_side.at(static_cast<size_t>(
                    box.skeleton_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    std::vector<DataType> interpolated;
    h2_matrix_vector_product(
        box.interpolation_matrix, skeleton, interpolated, data.nrhs, 'T');
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < redundant_count; ++i) {
            data.left_side.at(static_cast<size_t>(
                box.redundant_indices[static_cast<size_t>(i)] +
                column * data.num_points)) +=
                interpolated[static_cast<size_t>(i + column * redundant_count)];
        }
    }
}

template<typename CoordType, typename DataType>
void apply_h2_interactions(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<SolveDataRequest<CoordType, DataType>>& source_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& target_data,
    const std::unordered_map<int64_t, std::vector<DataType>>& remote_sources) {

    std::exception_ptr apply_exception;
    std::mutex apply_exception_mutex;
    std::atomic<bool> apply_failed{false};

    #pragma omp parallel for schedule(dynamic) if (level.local_boxes.size() > 1)
    for (int64_t box_index = 0;
         box_index < static_cast<int64_t>(level.local_boxes.size());
         ++box_index) {
        if (apply_failed.load(std::memory_order_relaxed)) continue;
        try {
            const auto& target = level.local_boxes[static_cast<size_t>(box_index)];
            auto& target_vector = target_data[static_cast<size_t>(box_index)].left_side;
            for (const auto& block : target.h2_interaction_blocks) {
                std::vector<DataType> source;
                if (level.find_local_box(block.source_morton) != nullptr) {
                    source = h2_local_vector_for_morton(
                        level, source_data, block.source_morton, true);
                } else {
                    auto it = remote_sources.find(block.source_morton);
                    if (it == remote_sources.end()) {
                        throw std::runtime_error(
                            "apply_h2_interactions: missing remote multipole vector");
                    }
                    source = it->second;
                }

                std::vector<DataType> contribution;
                const int64_t nrhs =
                    target_data[static_cast<size_t>(box_index)].nrhs;
                h2_matrix_vector_product(
                    block.matrix, source, contribution, nrhs, 'N');
                const int64_t skeleton_count =
                    static_cast<int64_t>(target.skeleton_indices.size());
                if (contribution.size() !=
                    static_cast<size_t>(skeleton_count * nrhs)) {
                    throw std::runtime_error(
                        "apply_h2_interactions: target rank mismatch");
                }
                const int64_t target_rows =
                    target_data[static_cast<size_t>(box_index)].num_points;
                for (int64_t column = 0; column < nrhs; ++column) {
                    for (int64_t i = 0; i < skeleton_count; ++i) {
                        target_vector.at(static_cast<size_t>(
                            target.skeleton_indices[static_cast<size_t>(i)] +
                            column * target_rows)) += contribution[static_cast<size_t>(
                                i + column * skeleton_count)];
                    }
                }
            }
        } catch (...) {
            if (!apply_failed.exchange(true, std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(apply_exception_mutex);
                apply_exception = std::current_exception();
            }
        }
    }
    if (apply_exception) std::rethrow_exception(apply_exception);
}

template<typename CoordType, typename DataType>
void apply_h2_leaf_near(
    TreeLevel<CoordType, DataType>& leaf,
    const std::vector<SolveDataRequest<CoordType, DataType>>& source_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& target_data,
    const std::unordered_map<int64_t, std::vector<DataType>>& remote_sources) {

    std::exception_ptr apply_exception;
    std::mutex apply_exception_mutex;
    std::atomic<bool> apply_failed{false};

    #pragma omp parallel for schedule(dynamic) if (leaf.local_boxes.size() > 1)
    for (int64_t box_index = 0;
         box_index < static_cast<int64_t>(leaf.local_boxes.size());
         ++box_index) {
        if (apply_failed.load(std::memory_order_relaxed)) continue;
        try {
            const auto& target = leaf.local_boxes[static_cast<size_t>(box_index)];
            auto& output = target_data[static_cast<size_t>(box_index)].left_side;
            for (const auto& block : target.h2_near_blocks) {
                std::vector<DataType> source;
                if (leaf.find_local_box(block.source_morton) != nullptr) {
                    source = h2_local_vector_for_morton(
                        leaf, source_data, block.source_morton, false);
                } else {
                    auto it = remote_sources.find(block.source_morton);
                    if (it == remote_sources.end()) {
                        throw std::runtime_error(
                            "apply_h2_leaf_near: missing remote source vector");
                    }
                    source = it->second;
                }
                std::vector<DataType> contribution;
                const int64_t nrhs =
                    target_data[static_cast<size_t>(box_index)].nrhs;
                h2_matrix_vector_product(
                    block.matrix, source, contribution, nrhs, 'N');
                if (contribution.size() != output.size()) {
                    throw std::runtime_error(
                        "apply_h2_leaf_near: target size mismatch");
                }
                for (size_t i = 0; i < contribution.size(); ++i) {
                    output[i] += contribution[i];
                }
            }
        } catch (...) {
            if (!apply_failed.exchange(true, std::memory_order_relaxed)) {
                std::lock_guard<std::mutex> lock(apply_exception_mutex);
                apply_exception = std::current_exception();
            }
        }
    }
    if (apply_exception) std::rethrow_exception(apply_exception);
}

template<typename CoordType, typename DataType>
void hierarchical_h2_mul_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    int nrhs,
    bool verbose) {

    if (nrhs <= 0 || input.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::invalid_argument(
            "hierarchical_h2_mul_parallel: invalid batched input dimensions");
    }
    const int64_t local_points =
        static_cast<int64_t>(input.size() / static_cast<size_t>(nrhs));

    const int rank = tree->mpi_rank;
    const int leaf_level = tree->num_levels - 1;
    const int first_h2_level = std::min(2, leaf_level);
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> source_data(
        static_cast<size_t>(tree->num_levels));
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> target_data(
        static_cast<size_t>(tree->num_levels));

    int64_t input_offset = 0;
    for (int level_number = first_h2_level;
         level_number <= leaf_level;
         ++level_number) {
        auto& level = tree->levels[level_number];
        if (!level.is_process_active) continue;
        source_data[level_number].resize(level.local_boxes.size());
        target_data[level_number].resize(level.local_boxes.size());
        for (size_t box_index = 0; box_index < level.local_boxes.size(); ++box_index) {
            const auto& box = level.local_boxes[box_index];
            auto& source = source_data[level_number][box_index];
            auto& target = target_data[level_number][box_index];
            source.initialize(box.morton_index, rank, box.num_points, nrhs);
            target.initialize(box.morton_index, rank, box.num_points, nrhs);
            source.skeleton_indices = box.skeleton_indices;
            source.redundant_indices = box.redundant_indices;
            target.skeleton_indices = box.skeleton_indices;
            target.redundant_indices = box.redundant_indices;

            if (level_number == leaf_level) {
                for (int column = 0; column < nrhs; ++column) {
                    for (int64_t i = 0; i < box.num_points; ++i) {
                        source.left_side[static_cast<size_t>(
                            i + static_cast<int64_t>(column) * box.num_points)] =
                            input[static_cast<size_t>(
                                input_offset + i +
                                static_cast<int64_t>(column) * local_points)];
                    }
                }
                input_offset += box.num_points;
                source.right_side = source.left_side;
            }
        }
    }
    if (input_offset != local_points) {
        throw std::runtime_error(
            "hierarchical_h2_mul_parallel: local input length does not match leaf DOFs");
    }

    if (leaf_level >= 2) {
        // Upward pass: q_B = x_B[S] + T_B x_B[R].
        for (int level_number = leaf_level; level_number >= 2; --level_number) {
            auto& level = tree->levels[level_number];
            if (level.is_process_active) {
                #pragma omp parallel for schedule(static) if (level.local_boxes.size() > 1)
                for (int64_t box_index = 0;
                     box_index < static_cast<int64_t>(level.local_boxes.size());
                     ++box_index) {
                    apply_h2_upward_projection(
                        level.local_boxes[static_cast<size_t>(box_index)],
                        source_data[level_number][static_cast<size_t>(box_index)]);
                }
            }
            if (level_number > 2) {
                gather_skeleton_to_parent(
                    level, tree->levels[level_number - 1],
                    source_data[level_number], source_data[level_number - 1],
                    tree->dimension, tree->comm);
            }
        }

        // Interaction and downward passes. Child skeleton entries are still
        // zero when parent data is scattered, so the existing replace scatter
        // gives the required inherited value before local interactions add in.
        for (int level_number = 2; level_number <= leaf_level; ++level_number) {
            auto& level = tree->levels[level_number];
            std::vector<int64_t> needed;
            if (level.is_process_active) {
                for (const auto& box : level.local_boxes) {
                    for (const auto& block : box.h2_interaction_blocks) {
                        if (level.find_local_box(block.source_morton) == nullptr) {
                            needed.push_back(block.source_morton);
                        }
                    }
                }
            }
            const auto remote_sources = exchange_h2_vectors_onehop(
                tree, level_number, source_data[level_number],
                std::move(needed), true, 800 + 8 * level_number);

            if (level.is_process_active) {
                apply_h2_interactions(
                    level, source_data[level_number], target_data[level_number],
                    remote_sources);
                #pragma omp parallel for schedule(static) if (level.local_boxes.size() > 1)
                for (int64_t box_index = 0;
                     box_index < static_cast<int64_t>(level.local_boxes.size());
                     ++box_index) {
                    apply_h2_downward_interpolation(
                        level.local_boxes[static_cast<size_t>(box_index)],
                        target_data[level_number][static_cast<size_t>(box_index)]);
                }
            }

            if (level_number < leaf_level) {
                scatter_solution_to_children(
                    tree->levels[level_number + 1], level,
                    target_data[level_number + 1], target_data[level_number],
                    tree->dimension, tree->comm);
            }
        }
    }

    auto& leaf = tree->levels[leaf_level];
    std::vector<int64_t> near_needed;
    if (leaf.is_process_active) {
        for (const auto& box : leaf.local_boxes) {
            for (const auto& block : box.h2_near_blocks) {
                if (leaf.find_local_box(block.source_morton) == nullptr) {
                    near_needed.push_back(block.source_morton);
                }
            }
        }
    }
    const auto remote_near_sources = exchange_h2_vectors_onehop(
        tree, leaf_level, source_data[leaf_level],
        std::move(near_needed), false, 1800 + 8 * leaf_level);
    if (leaf.is_process_active) {
        apply_h2_leaf_near(
            leaf, source_data[leaf_level], target_data[leaf_level],
            remote_near_sources);
    }

    output.assign(input.size(), DataType{0});
    if (leaf.is_process_active) {
        int64_t local_row = 0;
        for (const auto& box_data : target_data[leaf_level]) {
            for (int column = 0; column < nrhs; ++column) {
                for (int64_t i = 0; i < box_data.num_points; ++i) {
                    output[static_cast<size_t>(
                        local_row + i +
                        static_cast<int64_t>(column) * local_points)] =
                        box_data.left_side[static_cast<size_t>(
                            i + static_cast<int64_t>(column) *
                                box_data.num_points)];
                }
            }
            local_row += box_data.num_points;
        }
        if (local_row != local_points) {
            throw std::runtime_error(
                "hierarchical_h2_mul_parallel: local output row mismatch");
        }
    }
    if (output.size() != input.size()) {
        throw std::runtime_error(
            "hierarchical_h2_mul_parallel: local output length mismatch");
    }

    if (verbose && rank == smallest_active_rank(leaf)) {
        std::cout << "H2 compression-only multiply complete" << std::endl;
    }
}

template<typename DataType>
std::vector<DataType> h2_global_column_dots(
    const std::vector<DataType>& lhs,
    const std::vector<DataType>& rhs,
    int nrhs,
    MPI_Comm comm) {

    if (nrhs <= 0 || lhs.size() != rhs.size() ||
        lhs.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::runtime_error(
            "h2_global_column_dots: invalid batched vector dimensions");
    }
    const size_t local_rows = lhs.size() / static_cast<size_t>(nrhs);
    std::vector<DataType> local(static_cast<size_t>(nrhs), DataType{0});
    for (int column = 0; column < nrhs; ++column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        for (size_t row = 0; row < local_rows; ++row) {
            if constexpr (is_complex_v<DataType>) {
                local[static_cast<size_t>(column)] +=
                    std::conj(lhs[offset + row]) * rhs[offset + row];
            } else {
                local[static_cast<size_t>(column)] +=
                    lhs[offset + row] * rhs[offset + row];
            }
        }
    }
    std::vector<DataType> global(static_cast<size_t>(nrhs), DataType{0});
    MPI_Allreduce(
        local.data(), global.data(), nrhs, mpi_datatype_for<DataType>(),
        MPI_SUM, comm);
    return global;
}

template<typename DataType>
std::vector<double> h2_global_column_norms(
    const std::vector<DataType>& values,
    int nrhs,
    MPI_Comm comm) {

    if (nrhs <= 0 || values.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::runtime_error(
            "h2_global_column_norms: invalid batched vector dimensions");
    }
    const size_t local_rows = values.size() / static_cast<size_t>(nrhs);
    std::vector<double> local(static_cast<size_t>(nrhs), 0.0);
    for (int column = 0; column < nrhs; ++column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        for (size_t row = 0; row < local_rows; ++row) {
            if constexpr (is_complex_v<DataType>) {
                local[static_cast<size_t>(column)] +=
                    std::norm(values[offset + row]);
            } else {
                const double value = static_cast<double>(values[offset + row]);
                local[static_cast<size_t>(column)] += value * value;
            }
        }
    }
    std::vector<double> global(static_cast<size_t>(nrhs), 0.0);
    MPI_Allreduce(
        local.data(), global.data(), nrhs, MPI_DOUBLE, MPI_SUM, comm);
    for (double& value : global) value = std::sqrt(value);
    return global;
}

template<typename CoordType, typename DataType>
void hierarchical_h2_bicgstab_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<DataType>& solution,
    int nrhs,
    double tolerance,
    int max_iterations,
    int* completed_iterations = nullptr,
    double* final_relative_residual = nullptr,
    bool verbose = true) {

    if (tolerance <= 0.0) {
        throw std::invalid_argument(
            "hierarchical_h2_bicgstab_parallel: tolerance must be positive");
    }
    if (max_iterations <= 0) {
        throw std::invalid_argument(
            "hierarchical_h2_bicgstab_parallel: max_iterations must be positive");
    }
    if (nrhs <= 0 || rhs.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::invalid_argument(
            "hierarchical_h2_bicgstab_parallel: invalid batched RHS dimensions");
    }

    const size_t value_count = rhs.size();
    const size_t local_rows = value_count / static_cast<size_t>(nrhs);
    solution.assign(value_count, DataType{0});
    std::vector<DataType> residual = rhs;
    std::vector<DataType> shadow = residual;
    std::vector<DataType> direction(value_count, DataType{0});
    std::vector<DataType> image(value_count, DataType{0});
    std::vector<DataType> intermediate(value_count, DataType{0});
    std::vector<DataType> intermediate_image(value_count, DataType{0});

    const std::vector<double> rhs_norm =
        h2_global_column_norms(rhs, nrhs, tree->comm);
    std::vector<double> relative_residual(static_cast<size_t>(nrhs), 0.0);
    std::vector<bool> active(static_cast<size_t>(nrhs), false);
    std::vector<int> column_iterations(static_cast<size_t>(nrhs), 0);
    for (int column = 0; column < nrhs; ++column) {
        if (rhs_norm[static_cast<size_t>(column)] > 0.0) {
            active[static_cast<size_t>(column)] = true;
            relative_residual[static_cast<size_t>(column)] = 1.0;
        }
    }

    std::vector<DataType> rho_previous(
        static_cast<size_t>(nrhs), DataType{1});
    std::vector<DataType> alpha(static_cast<size_t>(nrhs), DataType{1});
    std::vector<DataType> omega(static_cast<size_t>(nrhs), DataType{1});
    auto unusable_scalar = [](const DataType& value) {
        const double magnitude = std::abs(value);
        return magnitude == 0.0 || !std::isfinite(magnitude);
    };
    auto any_active = [&]() {
        return std::any_of(active.begin(), active.end(), [](bool value) {
            return value;
        });
    };
    auto clear_column = [&](std::vector<DataType>& values, int column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        std::fill(
            values.begin() + static_cast<ptrdiff_t>(offset),
            values.begin() + static_cast<ptrdiff_t>(offset + local_rows),
            DataType{0});
    };

    for (int iteration = 1;
         iteration <= max_iterations && any_active();
         ++iteration) {
        const std::vector<DataType> rho =
            h2_global_column_dots(shadow, residual, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) {
                clear_column(direction, column);
                continue;
            }
            if (unusable_scalar(rho[static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "H2 BiCGSTAB breakdown: rho is zero for RHS " +
                    std::to_string(column));
            }
            const size_t offset = static_cast<size_t>(column) * local_rows;
            if (iteration == 1) {
                std::copy_n(
                    residual.begin() + static_cast<ptrdiff_t>(offset),
                    local_rows,
                    direction.begin() + static_cast<ptrdiff_t>(offset));
            } else {
                if (unusable_scalar(omega[static_cast<size_t>(column)])) {
                    throw std::runtime_error(
                        "H2 BiCGSTAB breakdown: omega is zero for RHS " +
                        std::to_string(column));
                }
                const DataType beta =
                    (rho[static_cast<size_t>(column)] /
                     rho_previous[static_cast<size_t>(column)]) *
                    (alpha[static_cast<size_t>(column)] /
                     omega[static_cast<size_t>(column)]);
                for (size_t row = 0; row < local_rows; ++row) {
                    const size_t index = offset + row;
                    direction[index] = residual[index] + beta *
                        (direction[index] -
                         omega[static_cast<size_t>(column)] * image[index]);
                }
            }
        }

        hierarchical_h2_mul_parallel(
            tree, direction, image, nrhs, false);
        const std::vector<DataType> denominator =
            h2_global_column_dots(shadow, image, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) {
                clear_column(intermediate, column);
                continue;
            }
            if (unusable_scalar(denominator[static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "H2 BiCGSTAB breakdown: alpha denominator is zero for RHS " +
                    std::to_string(column));
            }
            alpha[static_cast<size_t>(column)] =
                rho[static_cast<size_t>(column)] /
                denominator[static_cast<size_t>(column)];
            const size_t offset = static_cast<size_t>(column) * local_rows;
            for (size_t row = 0; row < local_rows; ++row) {
                const size_t index = offset + row;
                intermediate[index] = residual[index] -
                    alpha[static_cast<size_t>(column)] * image[index];
            }
        }

        const std::vector<double> intermediate_norm =
            h2_global_column_norms(intermediate, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            relative_residual[static_cast<size_t>(column)] =
                intermediate_norm[static_cast<size_t>(column)] /
                rhs_norm[static_cast<size_t>(column)];
            if (relative_residual[static_cast<size_t>(column)] <= tolerance) {
                const size_t offset = static_cast<size_t>(column) * local_rows;
                for (size_t row = 0; row < local_rows; ++row) {
                    const size_t index = offset + row;
                    solution[index] += alpha[static_cast<size_t>(column)] *
                        direction[index];
                }
                active[static_cast<size_t>(column)] = false;
                column_iterations[static_cast<size_t>(column)] = iteration;
                clear_column(intermediate, column);
            }
        }
        if (!any_active()) break;

        hierarchical_h2_mul_parallel(
            tree, intermediate, intermediate_image, nrhs, false);
        const std::vector<DataType> omega_numerator =
            h2_global_column_dots(
                intermediate_image, intermediate, nrhs, tree->comm);
        const std::vector<DataType> omega_denominator =
            h2_global_column_dots(
                intermediate_image, intermediate_image, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            if (unusable_scalar(
                    omega_denominator[static_cast<size_t>(column)])) {
                throw std::runtime_error(
                    "H2 BiCGSTAB breakdown: omega denominator is zero for RHS " +
                    std::to_string(column));
            }
            omega[static_cast<size_t>(column)] =
                omega_numerator[static_cast<size_t>(column)] /
                omega_denominator[static_cast<size_t>(column)];
            const size_t offset = static_cast<size_t>(column) * local_rows;
            for (size_t row = 0; row < local_rows; ++row) {
                const size_t index = offset + row;
                solution[index] += alpha[static_cast<size_t>(column)] *
                    direction[index] + omega[static_cast<size_t>(column)] *
                    intermediate[index];
                residual[index] = intermediate[index] -
                    omega[static_cast<size_t>(column)] *
                    intermediate_image[index];
            }
        }

        const std::vector<double> residual_norm =
            h2_global_column_norms(residual, nrhs, tree->comm);
        for (int column = 0; column < nrhs; ++column) {
            if (!active[static_cast<size_t>(column)]) continue;
            relative_residual[static_cast<size_t>(column)] =
                residual_norm[static_cast<size_t>(column)] /
                rhs_norm[static_cast<size_t>(column)];
            if (relative_residual[static_cast<size_t>(column)] <= tolerance) {
                active[static_cast<size_t>(column)] = false;
                column_iterations[static_cast<size_t>(column)] = iteration;
            } else {
                rho_previous[static_cast<size_t>(column)] =
                    rho[static_cast<size_t>(column)];
            }
        }
    }

    const int iterations_done = *std::max_element(
        column_iterations.begin(), column_iterations.end());
    const double max_relative_residual = *std::max_element(
        relative_residual.begin(), relative_residual.end());
    if (completed_iterations) *completed_iterations = iterations_done;
    if (final_relative_residual) {
        *final_relative_residual = max_relative_residual;
    }
    if (any_active()) {
        throw std::runtime_error(
            "H2 BiCGSTAB did not converge in " +
            std::to_string(max_iterations) +
            " iterations; maximum relative residual=" +
            std::to_string(max_relative_residual));
    }

    if (verbose && tree->mpi_rank == smallest_active_rank(
            tree->levels[tree->num_levels - 1])) {
        std::cout << "H2 BiCGSTAB converged in " << iterations_done
                  << " iterations, maximum relative residual="
                  << max_relative_residual
                  << std::endl;
    }
}

template<typename CoordType, typename DataType>
void butterfly_compression_parallel(
    H2<CoordType, DataType>* solver,
    double* compression_time,
    double* entryeval_time) {

    if constexpr (!(std::is_same_v<DataType, double> ||
                    std::is_same_v<DataType, std::complex<double>>)) {
        throw std::runtime_error(
            "H2 compression-only construction supports double precision only");
    } else {
        if (solver->build_state == H2BuildState::H2_COMPRESSED) {
            if (compression_time) *compression_time = 0.0;
            if (entryeval_time) *entryeval_time = 0.0;
            return;
        }
        if (solver->build_state == H2BuildState::RS_FACTORIZED) {
            throw std::runtime_error(
                "butterfly_compression_parallel: solver already contains an RS-S factorization");
        }

        solver->kernel.entryeval_time_per_thread.assign(omp_get_max_threads(), 0.0);
        const double start = MPI_Wtime();
        hierarchical_compression_parallel(
            solver->tree.get(), &solver->kernel, solver->options.tolerance,
            &solver->last_factor_rankmax, &solver->factorization_memory,
            solver->options.use_sketch,
            solver->options.verbosity >= 1);
        double elapsed = MPI_Wtime() - start;
        MPI_Allreduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
        if (compression_time) *compression_time = elapsed;

        double entry_time = 0.0;
        if (!solver->kernel.entryeval_time_per_thread.empty()) {
            entry_time = *std::max_element(
                solver->kernel.entryeval_time_per_thread.begin(),
                solver->kernel.entryeval_time_per_thread.end());
        }
        MPI_Allreduce(MPI_IN_PLACE, &entry_time, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
        if (entryeval_time) *entryeval_time = entry_time;

        MPI_Allreduce(
            MPI_IN_PLACE, &solver->last_factor_rankmax,
            1, MPI_INT64_T, MPI_MAX, solver->comm);
        solver->build_state = H2BuildState::H2_COMPRESSED;
        solver->factorized = false;
    }
}

} // namespace butterfly
