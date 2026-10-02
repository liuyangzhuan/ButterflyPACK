#pragma once

// H2 (format-7) integration: CA (communication-avoiding) level machinery
// shared by both backends - CA box grouping and scheduling, and the
// per-level CA factorization driver.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "factorization.hpp"
#include "memory_diagnostics.hpp"
#include "morton.hpp"
#include "owner_schedule.hpp"
#include "staged_halo.hpp"
#include "tree.hpp"

namespace butterfly {
using namespace fmm;

template<typename CoordType, typename DataType>
std::vector<std::vector<int64_t>> make_CA_box_groups(
    const TreeLevel<CoordType, DataType>& level,
    int dimension,
    bool reverse_order,
    bool include_boundary = true) {
    std::vector<std::vector<int64_t>> groups;
    std::vector<int64_t> interior;
    interior.reserve(level.interior_id.size());
    for (int64_t local_idx : level.interior_id) {
        interior.push_back(local_idx + level.local_morton_start);
    }

    if (!include_boundary) {
        groups.push_back(std::move(interior));
    } else if (level.num_active_processes == 1) {
        std::vector<int64_t> boundary;
        boundary.reserve(level.boundary_id.size());
        for (int64_t local_idx : level.boundary_id) {
            boundary.push_back(local_idx + level.local_morton_start);
        }
        groups.push_back(std::move(boundary));
        groups.push_back(std::move(interior));
    } else {
        std::vector<int64_t> blue_orange = level.blue;
        blue_orange.insert(
            blue_orange.end(), level.orange.begin(), level.orange.end());
        groups.push_back(std::move(blue_orange));
        if (dimension == 3) {
            groups.push_back(level.purple);
        }
        groups.push_back(level.green);
        groups.push_back(std::move(interior));
    }

    if (reverse_order) {
        std::reverse(groups.begin(), groups.end());
    }
    return groups;
}

template<typename CoordType, typename DataType>
SolveDataRequest<CoordType, DataType>* resolve_CA_solve_box(
    TreeLevel<CoordType, DataType>& level,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int64_t morton,
    bool& is_ghost) {
    if (level.is_box_on_process(morton)) {
        is_ghost = false;
        return &level_solve_data[static_cast<size_t>(
            morton - level.local_morton_start)];
    }

    auto it = level.ghost_and_assisting_box_points_for_solve_map.find(morton);
    if (it == level.ghost_and_assisting_box_points_for_solve_map.end() ||
        !level.is_ghost_solve[it->second]) {
        return nullptr;
    }
    is_ghost = true;
    return &level.ghost_and_assisting_boxes_for_solve[it->second];
}

template<typename CoordType, typename DataType, typename BoxOperation,
         typename GroupComplete>
void apply_CA_level_schedule(
    TreeLevel<CoordType, DataType>& level,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int dimension,
    bool reverse_order,
    bool writes_neighbors,
    BoxOperation&& operation,
    GroupComplete&& group_complete,
    bool include_boundary = true) {
    if (!level.is_process_active) return;

    auto groups = make_CA_box_groups(
        level, dimension, reverse_order, include_boundary);
    const int num_waves = 1 << dimension;
    const int max_threads = std::max(1, omp_get_max_threads());

    for (size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {
        std::vector<std::vector<int64_t>> waves(
            static_cast<size_t>(num_waves));
        for (int64_t morton : groups[group_idx]) {
            waves[static_cast<size_t>(morton & (num_waves - 1))].push_back(morton);
        }

        for (int wave_step = 0; wave_step < num_waves; ++wave_step) {
            const int wave_idx = reverse_order
                ? num_waves - 1 - wave_step
                : wave_step;
            auto& wave = waves[static_cast<size_t>(wave_idx)];
            if (wave.empty()) continue;

            std::vector<PendingSolveUpdates<DataType>> thread_pending(
                static_cast<size_t>(max_threads));
            std::exception_ptr wave_exception;
            std::mutex wave_exception_mutex;
            std::atomic<bool> wave_failed{false};

            #pragma omp parallel default(shared) if (wave.size() > 1)
            {
                const int tid = omp_get_thread_num();
                auto& pending = thread_pending[static_cast<size_t>(tid)];

                #pragma omp for schedule(static)
                for (int64_t idx = 0;
                     idx < static_cast<int64_t>(wave.size());
                     ++idx) {
                    if (wave_failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const size_t ordered_idx = reverse_order
                            ? wave.size() - 1 - static_cast<size_t>(idx)
                            : static_cast<size_t>(idx);
                        const int64_t morton = wave[ordered_idx];
                        bool is_ghost = false;
                        auto* solve_box = resolve_CA_solve_box(
                            level, level_solve_data, morton, is_ghost);
                        if (solve_box == nullptr) {
                            throw std::runtime_error(
                                "CA solve schedule cannot resolve Morton box " +
                                std::to_string(morton));
                        }
                        operation(*solve_box, is_ghost, pending);
                    } catch (...) {
                        if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(
                                wave_exception_mutex);
                            wave_exception = std::current_exception();
                        }
                    }
                }
            }

            if (wave_exception) std::rethrow_exception(wave_exception);

            if (writes_neighbors) {
                PendingSolveUpdates<DataType> merged;
                for (int tid = 0; tid < max_threads; ++tid) {
                    merge_pending_solve(
                        merged, thread_pending[static_cast<size_t>(tid)]);
                    clear_pending_solve_updates_memory(
                        thread_pending[static_cast<size_t>(tid)]);
                }
                apply_pending_solve_updates_local_or_ghost(
                    level, level_solve_data, merged);
                clear_pending_solve_updates_memory(merged);
            }
        }

        group_complete(group_idx, groups.size());
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void factorize_CA_level(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int current_level,
    KernelType* kernel,
    double tolerance,
    int use_sketch,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t& total_skeleton,
    int64_t& total_redundant,
    int64_t& local_max_skel,
    bool print_detail,
    int level_print_rank,
    H2FactorizationMemoryDiagnostics* memory_diagnostics,
    OwnerScheduleState<CoordType, DataType>* owner_schedule,
    StagedHaloState<CoordType, DataType>* staged_state,
    bool staged_overlap_scheduling) {
    auto& level = tree->levels[current_level];
    const int dimension = tree->dimension;
    const int num_waves = 1 << dimension;

    std::vector<std::string> group_names;
    std::vector<const std::vector<int64_t>*> group_lists;
    std::vector<int64_t> boundary;
    std::vector<int64_t> interior;
    std::vector<int64_t> blue_orange;

    interior.reserve(level.interior_id.size());
    for (int64_t local_idx : level.interior_id) {
        interior.push_back(local_idx + level.local_morton_start);
    }

    const bool use_owner_schedule =
        owner_schedule != nullptr && owner_schedule->active;

    if (use_owner_schedule) {
        group_names = {"interior"};
        group_lists = {&interior};
    } else if (level.num_active_processes == 1) {
        boundary.reserve(level.boundary_id.size());
        for (int64_t local_idx : level.boundary_id) {
            boundary.push_back(local_idx + level.local_morton_start);
        }
        group_names = {"boundary", "interior"};
        group_lists = {&boundary, &interior};
    } else {
        blue_orange = level.blue;
        blue_orange.insert(
            blue_orange.end(), level.orange.begin(), level.orange.end());
        group_names.push_back("blue/orange");
        group_lists.push_back(&blue_orange);
        if (dimension == 3) {
            group_names.push_back("purple");
            group_lists.push_back(&level.purple);
        }
        group_names.push_back("green");
        group_lists.push_back(&level.green);
        group_names.push_back("interior");
        group_lists.push_back(&interior);
    }

    auto complete_staged_stage = [&](int stage) {
        if (!staged_overlap_scheduling || staged_state == nullptr ||
            staged_state->stage_arrived[stage]) {
            return;
        }
        staged_halo_wait_stage(level, *staged_state, stage);
        kernel->register_level_coordinates(level);
        const auto merge_start = std::chrono::high_resolution_clock::now();
        staged_halo_merge_stage(
            level, *staged_state, stage, kernel);
        staged_state->t_merge_ms[stage] +=
            std::chrono::duration<double, std::milli>(
                std::chrono::high_resolution_clock::now() - merge_start)
                .count();
    };

    if (use_owner_schedule) {
        for (int stage = 1; stage < STAGED_HALO_STAGES; ++stage) {
            complete_staged_stage(stage);
        }
        owner_schedule_check_presence(level, *owner_schedule);
        owner_schedule_run_boundary(
            *owner_schedule, level, kernel,
            unit_proxy_points, num_proxy, proxy_radius,
            tolerance, is_symmetric, is_hermitian,
            factorization_method, use_sketch);
    }

    for (size_t group_idx = 0; group_idx < group_names.size(); ++group_idx) {
        const std::string& group_name = group_names[group_idx];
        const auto& group = *group_lists[group_idx];

        if (staged_overlap_scheduling) {
            if (group_name == "blue/orange") {
                complete_staged_stage(1);
            } else if (group_name == "boundary") {
                complete_staged_stage(1);
                complete_staged_stage(2);
                complete_staged_stage(3);
            } else if (group_name == "purple") {
                complete_staged_stage(2);
            } else if (group_name == "green") {
                complete_staged_stage(2);
                complete_staged_stage(3);
            } else if (group_name == "interior") {
                complete_staged_stage(2);
                complete_staged_stage(3);
                if (!level.staged_unarrived.empty()) {
                    throw std::runtime_error(
                        "staged halo overlap: unarrived ghosts at interior");
                }
            }
        }
        if (print_detail && tree->mpi_rank == level_print_rank) {
            std::cout << "  Processing CA " << group_name << " group ("
                      << group.size() << " boxes)..." << std::endl;
        }

        std::vector<std::vector<int64_t>> waves(
            static_cast<size_t>(num_waves));
        for (int64_t morton : group) {
            waves[static_cast<size_t>(morton & (num_waves - 1))].push_back(morton);
        }

        if (group_name == "interior") {
            clear_ghosts(level);
        }

        for (size_t wave_idx = 0; wave_idx < waves.size(); ++wave_idx) {
            const auto& wave = waves[wave_idx];
            if (wave.empty()) continue;
            std::unordered_set<int64_t> wave_boxes(wave.begin(), wave.end());
            for (int64_t morton : wave) {
                auto* box = level.find_local_box(morton);
                if (box == nullptr) box = level.find_ghost_box(morton);
                if (box == nullptr) continue;
                for (int64_t neighbor : box->one_hop) {
                    if (neighbor != morton && wave_boxes.count(neighbor) != 0) {
                        throw std::runtime_error(
                            "CA wave contains one-hop neighbors at level " +
                            std::to_string(current_level));
                    }
                }
            }

            const bool use_owner_deferred_xnn =
                is_symmetric && !is_hermitian;
            const int max_wave_threads = std::max(1, omp_get_max_threads());
            std::vector<std::vector<int64_t>> thread_xnn_candidate_boxes;
            std::vector<size_t> candidate_box_offsets;
            std::vector<int64_t> wave_xnn_candidate_boxes;
            std::vector<std::vector<DeferredXnnTargetKey>>
                wave_xnn_mirror_targets;
            if (use_owner_deferred_xnn) {
                thread_xnn_candidate_boxes.resize(
                    static_cast<size_t>(max_wave_threads));
                candidate_box_offsets.resize(
                    thread_xnn_candidate_boxes.size() + 1, 0);
            }

            std::exception_ptr wave_exception;
            std::mutex wave_exception_mutex;
            std::atomic<bool> wave_failed{false};
            size_t wave_scratch_bytes = 0;

            const int wave_split = split_threads_for(
                static_cast<int64_t>(wave.size()), max_wave_threads);

            #pragma omp parallel default(shared)
            {
                FactorizationThreadScratch<CoordType, DataType> scratch;
                scratch.split_threads = wave_split;
                #pragma omp for schedule(dynamic)
                for (int64_t idx = 0; idx < static_cast<int64_t>(wave.size()); ++idx) {
                    if (wave_failed.load(std::memory_order_relaxed)) continue;
                    try {
                        const int64_t morton = wave[static_cast<size_t>(idx)];
                        auto* box = level.find_local_box(morton);
                        if (box == nullptr) box = level.find_ghost_box(morton);
                        if (box == nullptr) {
                            throw std::runtime_error(
                                "CA factorization cannot resolve Morton box " +
                                std::to_string(morton));
                        }

                        if (stream_sketch_enabled()) {
                            gather_id_target_streamed(
                                tree, box, level, kernel,
                                scratch, box->on_boundary, tolerance);
                        } else {
                            scratch.streamed_sketch_valid = false;
                            gather_id_workspace(
                                tree,
                                box, level, kernel, tolerance,
                                unit_proxy_points.data(), num_proxy,
                                proxy_radius, is_symmetric,
                                scratch.workspace, scratch.workspace_rows,
                                scratch.workspace_cols, 0, box->on_boundary,
                                /*use_CA_boundary_semantics=*/true,
                                &scratch.id_adaptive);
                        }
                        compute_and_modify(
                            dimension, box, level, kernel, scratch,
                            tolerance, use_sketch, is_symmetric, is_hermitian,
                            static_cast<PendingFactorUpdates<DataType>*>(nullptr),
                            factorization_method, use_owner_deferred_xnn);
                    } catch (...) {
                        if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(wave_exception_mutex);
                            wave_exception = std::current_exception();
                        }
                    }
                }

                if (memory_diagnostics != nullptr && memory_diagnostics->enabled()) {
                    const size_t thread_scratch_bytes =
                        h2_diag_factorization_scratch_bytes(scratch);
                    #pragma omp atomic update
                    wave_scratch_bytes += thread_scratch_bytes;
                    #pragma omp barrier
                    #pragma omp single
                    {
                        memory_diagnostics->record(
                            tree, current_level, "CA",
                            "group" + std::to_string(group_idx) + "_wave" +
                                std::to_string(wave_idx) + "_compute",
                            0, wave_scratch_bytes);
                    }
                }
            }

            if (wave_exception) std::rethrow_exception(wave_exception);

            const int32_t wave_sequence = static_cast<int32_t>(
                use_owner_schedule
                    ? owner_schedule->num_boundary_colors *
                          owner_schedule->num_waves +
                          static_cast<int>(wave_idx)
                    : static_cast<int>(group_idx) * num_waves +
                          static_cast<int>(wave_idx));
            for (int64_t morton : wave) {
                level.eliminated_boxes.insert(morton);
                level.elimination_wave[morton] = wave_sequence;
            }
            slice_far_field_blocks(level, is_symmetric, is_hermitian);

            if (use_owner_deferred_xnn) {
                std::exception_ptr candidate_exception;
                std::mutex candidate_exception_mutex;
                std::atomic<bool> candidate_failed{false};

                #pragma omp parallel default(shared) if (wave.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    auto& local_candidates =
                        thread_xnn_candidate_boxes[static_cast<size_t>(tid)];

                    #pragma omp for schedule(static)
                    for (int64_t idx = 0;
                         idx < static_cast<int64_t>(wave.size());
                         ++idx) {
                        if (candidate_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }
                        try {
                            const int64_t morton =
                                wave[static_cast<size_t>(idx)];
                            auto* source_box = level.find_local_box(morton);
                            if (source_box == nullptr) {
                                source_box = level.find_ghost_box(morton);
                            }
                            if (source_box == nullptr) continue;

                            collect_owner_deferred_xnn_candidates_for_source_box(
                                source_box, level, wave_boxes,
                                local_candidates,
                                /*include_ghosts=*/true);
                        } catch (...) {
                            if (!candidate_failed.exchange(
                                    true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(
                                    candidate_exception_mutex);
                                candidate_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (candidate_exception) {
                    std::rethrow_exception(candidate_exception);
                }

                for (size_t tid = 0;
                     tid < thread_xnn_candidate_boxes.size();
                     ++tid) {
                    candidate_box_offsets[tid + 1] =
                        candidate_box_offsets[tid] +
                        thread_xnn_candidate_boxes[tid].size();
                }
                wave_xnn_candidate_boxes.resize(
                    candidate_box_offsets.back());
                for (size_t tid = 0;
                     tid < thread_xnn_candidate_boxes.size();
                     ++tid) {
                    size_t out_idx = candidate_box_offsets[tid];
                    for (int64_t candidate :
                         thread_xnn_candidate_boxes[tid]) {
                        wave_xnn_candidate_boxes[out_idx++] = candidate;
                    }
                }
                std::sort(
                    wave_xnn_candidate_boxes.begin(),
                    wave_xnn_candidate_boxes.end());
                wave_xnn_candidate_boxes.erase(
                    std::unique(
                        wave_xnn_candidate_boxes.begin(),
                        wave_xnn_candidate_boxes.end()),
                    wave_xnn_candidate_boxes.end());
                wave_xnn_mirror_targets.resize(
                    wave_xnn_candidate_boxes.size());

                std::exception_ptr owner_exception;
                std::mutex owner_exception_mutex;
                std::atomic<bool> owner_failed{false};
                size_t owner_scratch_bytes = 0;

                const int owner_split = split_threads_for(
                    static_cast<int64_t>(wave_xnn_candidate_boxes.size()),
                    max_wave_threads);

                #pragma omp parallel default(shared)
                {
                    DeferredXnnOwnerScratch<DataType> owner_scratch;
                    owner_scratch.split_threads = owner_split;

                    #pragma omp for schedule(dynamic)
                    for (int64_t idx = 0;
                         idx < static_cast<int64_t>(
                                   wave_xnn_candidate_boxes.size());
                         ++idx) {
                        if (owner_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }
                        try {
                            const int64_t candidate =
                                wave_xnn_candidate_boxes[
                                    static_cast<size_t>(idx)];
                            apply_owner_deferred_xnn_updates_for_candidate_box(
                                candidate, level, kernel, wave_boxes,
                                owner_scratch,
                                wave_xnn_mirror_targets[
                                    static_cast<size_t>(idx)],
                                static_cast<PendingFactorUpdates<DataType>*>(
                                    nullptr),
                                /*include_ghosts=*/true);
                        } catch (...) {
                            if (!owner_failed.exchange(
                                    true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(
                                    owner_exception_mutex);
                                owner_exception = std::current_exception();
                            }
                        }
                    }

                    if (memory_diagnostics != nullptr && memory_diagnostics->enabled()) {
                        const size_t thread_scratch_bytes =
                            h2_diag_deferred_owner_scratch_bytes(owner_scratch);
                        #pragma omp atomic update
                        owner_scratch_bytes += thread_scratch_bytes;
                        #pragma omp barrier
                        #pragma omp single
                        {
                            memory_diagnostics->record(
                                tree, current_level, "CA",
                                "group" + std::to_string(group_idx) + "_wave" +
                                    std::to_string(wave_idx) + "_owner",
                                0, owner_scratch_bytes);
                        }
                    }
                }

                if (owner_exception) {
                    std::rethrow_exception(owner_exception);
                }

                std::exception_ptr mirror_exception;
                std::mutex mirror_exception_mutex;
                std::atomic<bool> mirror_failed{false};

                #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0;
                         idx < static_cast<int64_t>(
                                   wave_xnn_candidate_boxes.size());
                         ++idx) {
                        if (mirror_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }
                        try {
                            const int64_t candidate =
                                wave_xnn_candidate_boxes[
                                    static_cast<size_t>(idx)];
                            auto* candidate_box =
                                level.find_local_box(candidate);
                            if (candidate_box == nullptr) {
                                candidate_box =
                                    level.find_ghost_box(candidate);
                            }
                            if (candidate_box == nullptr) continue;

                            apply_symmetric_owner_deferred_xnn_updates_for_candidate_box(
                                candidate_box,
                                wave_xnn_mirror_targets[
                                    static_cast<size_t>(idx)],
                                level,
                                /*include_ghosts=*/true);
                        } catch (...) {
                            if (!mirror_failed.exchange(
                                    true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(
                                    mirror_exception_mutex);
                                mirror_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (mirror_exception) {
                    std::rethrow_exception(mirror_exception);
                }

                if (!use_owner_schedule) {
                    share_symmetric_level_edges(level);
                }

                std::exception_ptr finalize_exception;
                std::mutex finalize_exception_mutex;
                std::atomic<bool> finalize_failed{false};

                #pragma omp parallel default(shared) if (wave.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0;
                         idx < static_cast<int64_t>(wave.size());
                         ++idx) {
                        if (finalize_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }
                        try {
                            const int64_t morton =
                                wave[static_cast<size_t>(idx)];
                            auto* source_box = level.find_local_box(morton);
                            if (source_box == nullptr) {
                                source_box = level.find_ghost_box(morton);
                            }
                            if (source_box == nullptr) continue;
                            finalize_deferred_xnn_source_box(source_box);
                        } catch (...) {
                            if (!finalize_failed.exchange(
                                    true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(
                                    finalize_exception_mutex);
                                finalize_exception =
                                    std::current_exception();
                            }
                        }
                    }
                }

                if (finalize_exception) {
                    std::rethrow_exception(finalize_exception);
                }
            }

            if (memory_diagnostics != nullptr && memory_diagnostics->enabled()) {
                memory_diagnostics->record(
                    tree, current_level, "CA",
                    "group" + std::to_string(group_idx) + "_wave" +
                        std::to_string(wave_idx) + "_complete");
            }
        }
    }

    for (const auto& box : level.local_boxes) {
        total_skeleton += static_cast<int64_t>(box.skeleton_indices.size());
        total_redundant += static_cast<int64_t>(box.redundant_indices.size());
        local_max_skel = std::max<int64_t>(
            local_max_skel, box.skeleton_indices.size());
    }
}

} // namespace butterfly
