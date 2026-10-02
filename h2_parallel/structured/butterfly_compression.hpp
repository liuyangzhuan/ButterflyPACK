#pragma once

#include "gpu_common/bpack_env.hpp"
#include "core/butterfly_types.hpp"
#include "butterfly_solve.hpp"
#include "core/butterfly_matvec.hpp"
#ifdef H2_HAVE_GPU
#include "h2_gpu/compression_gpu.hpp"
#include "h2_gpu/h2_matvec.hpp"
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
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
        tree, level, level_number, neighbor_ranks, needed, {}, false, true);
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

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_compression_parallel(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    int use_sketch,  // H2_use_sketch
    bool verbose = true) {

    const int rank = tree->mpi_rank;
    const int leaf_level = tree->num_levels - 1;
    int64_t local_rankmax = 0;

    if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
        std::cout << "\n========================================\n"
                  << "H2 Compression Only (Parallel MPI)\n"
                  << "========================================" << std::endl;
    }

    // H2_use_gpu: the IDs and blocks of every level on the device
    // (h2_gpu/compression_gpu.hpp); the conditions are the same on every rank.
    bool gpu_blocks = false, gpu_ids = false;
#ifdef H2_HAVE_GPU
    if (color_gpu_enabled()) {
        std::string reason;
        gpu_blocks = gpu::compression_supported(tree, kernel, &reason);
        // (adaptive ID rows, H2_ID_proxy 2: on the device when they are
        // selected against the sketch, H2_use_sketch 2)
        gpu_ids = gpu_blocks && use_sketch && (tree->id_proxy_mode != 2 || use_sketch == 2) &&
                  gpu::device_sketch_supported(tree, *gpu::evaluator_of(kernel->gpu_evaluator));
        if (gpu_blocks) gpu::begin_device_compression(tree->num_levels);
        // the application's evaluator before its first use (not in the
        // compression time, butterfly_compression_parallel)
        double evaluator_warm_up = 0.0;
        if constexpr (gpu::gpu_data_type<DataType>) {
            if (gpu_blocks) {
                evaluator_warm_up = gpu::warm_up_registered_evaluator<DataType>(tree, kernel->gpu_evaluator);
            }
        }
        h2_warmup_seconds() = evaluator_warm_up;
        if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
            std::cout << "  GPU compression: "
                      << (gpu_blocks ? (gpu_ids ? "IDs and blocks on the device"
                                                : "blocks on the device, IDs on the host")
                                     : "off (" + reason + ")")
                      << std::endl;
            if (evaluator_warm_up > 0.0) {
                std::printf("  [gpu] warm-up: evaluator %.2f s, before the levels (not in the compression time)\n",
                            evaluator_warm_up);
                std::fflush(stdout);
            }
        }
    }
#endif
    (void)gpu_ids;

    if (leaf_level < 2) {
        exchange_h2_point_metadata(tree, leaf_level, false, false, true);
        kernel->register_level_coordinates(tree->levels[leaf_level]);
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
            kernel->register_level_coordinates(level);

            if (level.is_process_active) {
                std::exception_ptr id_exception;
                std::mutex id_exception_mutex;
                std::atomic<bool> id_failed{false};

                // boxes whose ID runs here: all of them, or those the device left
                std::vector<int64_t> host_boxes;
#ifdef H2_HAVE_GPU
                if (gpu_ids) {
                    gpu::compression_stats() = gpu::CompressionStats{};
                    host_boxes = gpu::compress_level_ids(tree, level_number, kernel, tolerance);
                } else
#endif
                {
                    host_boxes.resize(level.local_boxes.size());
                    std::iota(host_boxes.begin(), host_boxes.end(), int64_t{0});
                }

                #pragma omp parallel for schedule(dynamic) if (host_boxes.size() > 1)
                for (int64_t host_index = 0;
                     host_index < static_cast<int64_t>(host_boxes.size());
                     ++host_index) {
                    if (id_failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const int64_t box_index = host_boxes[static_cast<size_t>(host_index)];
                        h2_skeletonize_box(
                            tree,
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
            // The exchange publishes a box's skeleton only once its owner has
            // marked the box eliminated (the factorization's rule); here every
            // local skeleton is final, so the boxes are marked for the
            // exchange (without it, coupling blocks with remote sources would
            // be dropped as empty).
            if (level.is_process_active) {
                for (const auto& box : level.local_boxes) {
                    level.eliminated_boxes.insert(box.morton_index);
                }
            }
            exchange_h2_point_metadata(
                tree, level_number, true, true, level_number == leaf_level);
            level.eliminated_boxes.clear();
            kernel->register_level_coordinates(level);
#ifdef H2_HAVE_GPU
            if (gpu_blocks) {
                if (!gpu_ids) gpu::compression_stats() = gpu::CompressionStats{};
                std::vector<std::vector<int64_t>> sources(level.local_boxes.size());
                if (level.is_process_active) {
                    for (size_t b = 0; b < level.local_boxes.size(); ++b) {
                        sources[b] = h2_interaction_list(tree, level.local_boxes[b]);
                    }
                }
                gpu::build_level_blocks(tree, level_number, kernel, sources, level_number == leaf_level);
            } else
#endif
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
#ifdef H2_HAVE_GPU
                    if (gpu_blocks) {
                        const auto& g = gpu::compression_stats();
                        std::printf("  [gpu] level %d compression (rank %d): IDs %.2f s (plan %.2f, device %.2f, store %.2f; "
                                    "%lld boxes, %lld on the host), blocks %.2f s (plan %.2f, device and copies %.2f, "
                                    "of which waiting for copies %.2f; %lld coupling, %lld near) | up %.2f GB, down %.2f GB, "
                                    "heap peak %.2f GB\n",
                                    level_number, rank, g.ids, g.id_plan, g.id_device, g.id_store,
                                    static_cast<long long>(g.id_boxes), static_cast<long long>(g.host_id_boxes), g.blocks,
                                    g.blocks_plan, g.blocks_device, g.blocks_wait,
                                    static_cast<long long>(g.interaction_blocks), static_cast<long long>(g.near_blocks),
                                    g.bytes_up / 1e9, g.bytes_down / 1e9, g.heap_peak / 1e9);
                        std::fflush(stdout);
                    }
#endif
                }
            }
        }
    }

#ifdef H2_HAVE_GPU
    // the device matvec needs every rank's blocks on its device
    if (gpu_blocks) gpu::commit_device_matvec(tree, verbose);
#endif

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

#ifdef H2_HAVE_GPU
    // the whole matvec on the device (h2_gpu/h2_matvec.hpp) when the
    // compression kept the blocks there on every rank
    if (gpu::run_device_h2_mul(tree, input, output, nrhs)) {
        if (gpu::device_matvec_check()) {  // BPACK_CHECK=matvec: the host matvec as the reference
            std::vector<DataType> host_output;
            gpu::device_matvec_suspended() = true;
            hierarchical_h2_mul_parallel(tree, input, host_output, nrhs, false);
            gpu::device_matvec_suspended() = false;
            double sums[2] = {0.0, 0.0};
            for (size_t i = 0; i < output.size(); ++i) {
                sums[0] += std::norm(output[i] - host_output[i]);
                sums[1] += std::norm(host_output[i]);
            }
            MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, tree->comm);
            if (rank == 0) {
                std::printf("GPU matvec check: |y_gpu - y_host| / |y_host| = %.3e\n",
                            sums[1] > 0.0 ? std::sqrt(sums[0] / sums[1]) : 0.0);
                std::fflush(stdout);
            }
        }
        if (verbose && rank == smallest_active_rank(tree->levels[leaf_level])) {
            std::cout << "H2 compression-only multiply complete" << std::endl;
        }
        return;
    }
#endif
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
#ifdef H2_HAVE_GPU
        if (gpu::device_matvec_usable(tree)) {
            const auto& m = gpu::matvec_stats();
            std::printf("  GPU matvec (rank %d, %lld calls since the compression): %.2f s (upward %.2f, interactions and "
                        "downward %.2f, near %.2f; of these, messages %.2f and host hand-offs %.2f; input/output "
                        "transfers %.2f)\n",
                        tree->mpi_rank, static_cast<long long>(m.calls), m.total, m.upward, m.coupling, m.near, m.mpi,
                        m.host_handoff, m.transfer);
            std::fflush(stdout);
        }
#endif
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
        configure_color_gpu(solver->options.use_gpu != 0);
        fmm::env::report_environment(solver->tree->comm);
#ifdef H2_HAVE_GPU
        fmm::gpu::tensor_core_gemm() = solver->options.use_gpu == 2;
        if (solver->options.use_gpu != 0) {
            fmm::gpu::begin_operator_build(solver->tree.get(), solver->tree->comm);
            fmm::gpu::invalidate_device_solve();
        }
#endif
        h2_warmup_seconds() = 0.0;
        const double start = MPI_Wtime();
        hierarchical_compression_parallel(
            solver->tree.get(), &solver->kernel, solver->options.tolerance,
            &solver->last_factor_rankmax, &solver->factorization_memory,
            solver->options.use_sketch,
            solver->options.verbosity >= 1);
        double elapsed = MPI_Wtime() - start - h2_warmup_seconds();  // (the GPU warm-up is not compression)
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

        if (solver->options.verbosity >= 0 && solver->tree->mpi_rank == 0) {
            std::cout << "\n========================================\n"
                      << "H2 Compression Only Complete\n"
                      << "  total time: "
                      << std::llround(elapsed * 1000.0) << " ms\n"
                      << "========================================\n"
                      << std::endl;
        }
    }
}

} // namespace butterfly
