#pragma once

#include "occupied_topology.hpp"
#include "parent_transition.hpp"
#include "owner_deferred.hpp"
#include "verification.hpp"

namespace butterfly {
namespace color_unstructured {
using namespace fmm;

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_unstructured(
    fmm::ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    int use_sketch,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    int lazy_schur,
    int gemm_split,
    int ca_staged_halo,
    int ca_owner_component,
    int ca_owner_serial,
    int verbosity,
    const OccupiedTopology& occupied_topology) {

    // To Do: NEED TO FIX KERNEL!!!!!
    using clock = std::chrono::high_resolution_clock;

    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    const bool print_summary = verbosity >= 0;
    const bool print_detail = verbosity >= 1;
    const bool print_trace = verbosity >= 2;
    H2FactorizationMemoryDiagnostics memory_diagnostics;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    FactorizationCommunicatorSet factorization_comms =
        make_factorization_communicators(tree);

    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int num_children = morton::children_per_box(dimension);
    const int factorization_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);

    if (use_sketch < 0 || use_sketch > 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: use_sketch must be 0, 1, or 2");
    }
    if (lazy_schur < 0 || lazy_schur > 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: lazy_schur must be 0, 1, or 2");
    }
    if (lazy_schur != 0 && use_sketch != 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: lazy_schur requires use_sketch=2");
    }
    if (gemm_split < 0) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: gemm_split must be nonnegative");
    }
    if (ca_staged_halo != 0 && ca_staged_halo != 2) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_staged_halo must be 0 or 2");
    }
    if (ca_owner_component != 0 && ca_owner_component != 3) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_component must be 0 or 3");
    }
    if (ca_owner_serial != 0 && ca_owner_serial != 1) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_serial must be 0 or 1");
    }
    if (ca_owner_component == 3 && lazy_schur == 0) {
        throw std::invalid_argument(
            "hierarchical_factorization_parallel: ca_owner_component=3 requires lazy_schur=1 or 2");
    }

    struct FactorizationRuntimeReset {
        ~FactorizationRuntimeReset() {
            configure_color_factorization_runtime(false, 0, 0);
            configure_ca_factorization_runtime(0, 0, false, false);
        }
    } factorization_runtime_reset;
    owner_solve_records().clear();

    if (dynamic_threading.enabled &&
        !tree->levels[leaf_level].is_process_active) {
        park_inactive_rank_on_service_cpu(dynamic_threading);
    }

    if (print_summary && rank == factorization_header_rank) {
        const char* factorization_name = "Unknown";
        switch (factorization_method) {
            case FactorizationMethod::CHOLESKY:
                factorization_name = "Cholesky";
                break;
            case FactorizationMethod::LU:
                factorization_name = "LU";
                break;
            case FactorizationMethod::BUNCH_KAUFMAN:
                factorization_name = "Bunch-Kaufman";
                break;
            case FactorizationMethod::NONE:
                factorization_name = "None";
                break;
        }
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Factorization (Parallel)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "Levels: " << num_levels << std::endl;
        std::cout << "Leaf level: " << leaf_level << std::endl;
        std::cout << "Dimension: " << dimension << "D" << std::endl;
        std::cout << "Tolerance: " << tolerance << std::endl;
        std::cout << "Matrix property: " << (is_symmetric ? "Symmetric" : "Nonsymmetric") << std::endl;
        std::cout << "Factorization: " << factorization_name << std::endl;
        printf("max thread: %d\n", omp_get_max_threads());
        std::cout << "========================================\n" << std::endl;
    }

    // ===== Main factorization loop: leaf level down to level 1 =====
    auto total_time = clock::now();
    auto segment_start = clock::now();
    clock::duration total_data_exchange_time{};
    clock::duration total_reduction_time{};
    int64_t local_max_skel = 0;

    enum DeepFactorizationPhase : size_t {
        DEEP_TRANSPORT_WALL = 0,
        DEEP_TRANSPORT_MPI,
        DEEP_POST_TRANSPORT,
        DEEP_WAVE_SETUP,
        DEEP_PRIMARY,
        DEEP_OWNER_PREP,
        DEEP_OWNER_REPLAY,
        DEEP_MIRROR,
        DEEP_SHARE,
        DEEP_GENERATOR,
        DEEP_FINALIZE,
        DEEP_BOOKKEEPING,
        DEEP_OTHER,
        DEEP_PHASE_COUNT
    };
    const std::array<const char*, DEEP_PHASE_COUNT> deep_phase_names = {{
        "transport wall",
        "  MPI inside transport",
        "post-transport refresh/slicing",
        "wave setup",
        "primary ID/elimination",
        "owner candidate preparation",
        "owner replay",
        "symmetric mirror",
        "symmetric edge sharing",
        "lazy generator packaging",
        "source finalization",
        "pending merge/bookkeeping",
        "other/unaccounted"
    }};

    memory_diagnostics.record(tree, leaf_level, "setup", "factor_start");

    for (int current_level = leaf_level; current_level >= 1; current_level--) {
        auto level_start = std::chrono::high_resolution_clock::now();
        clock::duration level_data_exchange{};
        clock::duration level_reduction{};
        std::array<double, DEEP_PHASE_COUNT> deep_phase_ms{};
        double deep_elimination_local_ms = 0.0;
        const bool profile_deep_level =
            print_detail &&
            current_level >= std::max(2, leaf_level - 2);
        auto record_deep_phase = [&](DeepFactorizationPhase phase,
                                     const clock::time_point& start) {
            if (!profile_deep_level) return;
            deep_phase_ms[static_cast<size_t>(phase)] +=
                std::chrono::duration<double, std::milli>(
                    clock::now() - start).count();
        };

        auto& level = tree->levels[current_level];
        auto& parent_level = tree->levels[current_level - 1];
        filter_level_topology(tree, current_level, occupied_topology);
        const bool participates_in_transition =
            level.is_process_active || parent_level.is_process_active;
        if (!participates_in_transition) {
            continue;
        }

        MPI_Comm level_comm =
            factorization_comms.level[static_cast<size_t>(current_level)];
        MPI_Comm transition_comm =
            factorization_comms.transition[static_cast<size_t>(current_level)];
        const int level_print_rank = smallest_active_rank(level);
        LevelThreadPlan level_thread_plan;
        if (level.is_process_active) {
            level_thread_plan = configure_level_thread_plan(
                dynamic_threading, level_comm, level);
        }
        const bool CA_requested = tree->level_requests_CA(current_level);
        const bool use_CA_level = tree->level_uses_CA(current_level);
        const bool use_streamed_level =
            current_level > 1 && use_sketch == 2 &&
            is_symmetric && !is_hermitian && tree->id_proxy_mode != 2;
        const int level_lazy_schur = use_streamed_level
            ? (use_CA_level ? std::min(lazy_schur, 1) : lazy_schur)
            : 0;
        configure_color_factorization_runtime(
            use_streamed_level, level_lazy_schur, gemm_split);
        configure_ca_factorization_runtime(
            use_CA_level ? ca_staged_halo : 0,
            use_CA_level ? ca_owner_component : 0,
            use_CA_level && ca_owner_serial != 0,
            use_CA_level && is_symmetric);

        OwnerScheduleState<CoordType, DataType> owner_schedule;
        StagedHaloState<CoordType, DataType> staged_state;
        bool staged_overlap_scheduling = false;
        if (level.is_process_active && use_CA_level &&
            current_level > 1 && level.num_active_processes > 1 &&
            ca_owner_component == 3) {
            if (!(is_symmetric && !is_hermitian)) {
                throw std::runtime_error(
                    "H2_CA_owner_component=3 requires a symmetric, non-Hermitian factorization");
            }
            if (!use_streamed_level ||
                lazy_far_field_mode() != LazyFarFieldMode::LAZY) {
                throw std::runtime_error(
                    "H2_CA_owner_component=3 requires streamed sketching and lazy Schur updates");
            }
            owner_schedule_setup(
                tree, current_level, level_comm, owner_schedule);
        }

        if (level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            FactorizationMemoryDiagnosticCallback gather_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                gather_memory_diagnostic =
                    [&](const char* phase, size_t pending, size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "CA",
                            std::string("initial_gather_") + phase,
                            pending, 0, 0, communication);
                    };
            }
            if (ca_staged_halo == 2 && is_symmetric && !is_hermitian) {
                initiate_staged_halo_gather(
                    tree, current_level, staged_state, level_comm,
                    owner_schedule.active
                        ? &owner_schedule.block_ghosts
                        : nullptr);
                staged_halo_wait_stage(level, staged_state, 0);

                const bool overlap_ok =
                    current_level > 1 &&
                    lazy_far_field_mode() == LazyFarFieldMode::LAZY;
                if (overlap_ok) {
                    build_staged_stage_map(
                        level, level.staged_stage_map);
                    level.staged_unarrived.clear();
                    for (int stage = 1;
                         stage < STAGED_HALO_STAGES; ++stage) {
                        for (const auto& peer_mortons :
                             staged_state.recv_mortons[stage]) {
                            level.staged_unarrived.insert(
                                peer_mortons.begin(),
                                peer_mortons.end());
                        }
                    }
                    level.staged_pending = std::make_unique<
                        typename TreeLevel<CoordType, DataType>::
                            StagedOverlapPending>();
                    level.staged_overlap_on = true;
                    staged_overlap_scheduling = true;
                } else {
                    for (int stage = 1;
                         stage < STAGED_HALO_STAGES; ++stage) {
                        staged_halo_wait_stage(
                            level, staged_state, stage);
                    }
                    staged_halo_finish(staged_state);
                }
            } else {
                gather_CA_factorization_data(
                    tree, current_level, level_comm,
                    is_symmetric && !is_hermitian,
                    gather_memory_diagnostic);
            }
            const auto gather_duration = clock::now() - segment_start;
            level_data_exchange += gather_duration;
            if (print_detail && rank == level_print_rank) {
                std::cout << "  CA initial ghost/assisting gather time: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 gather_duration).count()
                          << " ms" << std::endl;
            }
        }

        if (level.is_process_active) {
            kernel->register_level_coordinates(level);
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            level.is_process_active && use_CA_level ? "post_CA_gather" :
                                                       "level_start");

        // declare locks
        std::unordered_map<int64_t, omp_lock_t*> box_locks;

         // initialize locks
        if (level.is_process_active) {
            auto add_occupied_lock = [&](int64_t local_index) {
                auto& box = level.local_boxes[static_cast<size_t>(local_index)];
                auto* lock = new omp_lock_t;
                omp_init_lock(lock);
                box_locks[box.morton_index] = lock;
            };
            for (int64_t local_index : level.boundary_id) {
                add_occupied_lock(local_index);
            }
            for (int64_t local_index : level.interior_id) {
                add_occupied_lock(local_index);
            }
            if (use_CA_level) {
                for (auto& box : level.ghost_boxes) {
                    auto* lock = new omp_lock_t;
                    omp_init_lock(lock);
                    box_locks[box.morton_index] = lock;
                }
            }

            {
                auto* global_lock1 = new omp_lock_t;
                omp_init_lock(global_lock1);
                box_locks[-1] = global_lock1;
                auto* global_lock2 = new omp_lock_t;
                omp_init_lock(global_lock2);
                box_locks[-2] = global_lock2;
            }
            level.box_locks = std::move(box_locks);
        }

        if (print_summary && rank == level_print_rank) {
            std::cout << "\n===== Level " << current_level << " =====" << std::endl;
            std::cout << "Active processes: " << level.num_active_processes << " from rank: " << rank << std::endl;
            std::cout << "Boxes per process: " << level.num_boxes_local << std::endl;
        }
        print_level_thread_plan(
            dynamic_threading,
            current_level,
            "factorization",
            level_print_rank,
            rank,
            level_thread_plan,
            print_detail);
        if (print_detail && rank == level_print_rank && current_level > 1) {
            std::cout << "  Level algorithm: "
                      << (use_CA_level ? "CA" : "color");
            if (CA_requested && !use_CA_level) {
                std::cout << " (automatic CA fallback: "
                          << tree->boxes_per_active_process(current_level)
                          << " boxes per active process; distributed CA requires at least "
                          << tree->minimum_CA_boxes_per_active_process() << ")";
            }
            std::cout << std::endl;
        }

        // ===== Step 1: Gather ghost and assisting boxes =====



        // ===== Step 2: Eliminate boxes in colored order =====
        std::chrono::milliseconds elim_duration{0};
        int64_t total_skeleton = 0;
        int64_t total_redundant = 0;

        if (current_level > 1 && level.is_process_active) {
            // Regular levels: do ID and elimination in colored order

            auto elim_start = std::chrono::high_resolution_clock::now();

            if (use_CA_level) {
                factorize_CA_level(
                    tree, current_level, kernel, tolerance,
                    use_sketch,
                    is_symmetric, is_hermitian, factorization_method,
                    unit_proxy_points, num_proxy, proxy_radius,
                    total_skeleton, total_redundant, local_max_skel,
                    print_detail, level_print_rank, &memory_diagnostics,
                    &owner_schedule, &staged_state,
                    staged_overlap_scheduling);
            } else {
                const int num_colors = 1 << dimension;

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Both use morton & (num_colors-1) coloring, guaranteeing that
            // boxes within the same wave are non-adjacent.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx : level.boundary_id) {
                const int64_t morton_idx =
                    level.local_boxes[static_cast<size_t>(local_idx)].morton_index;
                const int color_id = static_cast<int>(
                    morton_idx & (num_colors - 1));
                color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
            }
            for (int64_t local_idx : level.interior_id) {
                const int64_t morton_idx =
                    level.local_boxes[static_cast<size_t>(local_idx)].morton_index;
                const int color_id = static_cast<int>(
                    morton_idx & (num_colors - 1));
                interior_sub_bins[static_cast<size_t>(color_id)].push_back(
                    morton_idx);
            }

            // Append interior sub-waves after boundary waves
            const int interior_start_loc = static_cast<int>(color_bins.size());  // == num_colors
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout:
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves

            PendingFactorUpdates<DataType> pending_updates;
            auto pending_updates_empty = [&]() {
                return pending_updates.replace_blocks.empty() &&
                    pending_updates.accumulated_deltas.empty() &&
                    pending_updates.generators.empty();
            };
            std::vector<int> lazy_interior_flush_after_wave(
                static_cast<size_t>(num_colors), 0);
            bool lazy_interior_schedule_ready =
                lazy_far_field_mode() != LazyFarFieldMode::LAZY;
            auto build_lazy_interior_flush_schedule = [&]() {
                std::vector<int> local_flush_after_wave(
                    static_cast<size_t>(num_colors), 0);

                // An interior wave can publish state only when one of its
                // occupied boxes crosses a process boundary. Generators use
                // this exact test in emit_lazy_generators; replacement and
                // additive updates are likewise routed through one-hop
                // owners. Requesters discovered by assisting exchanges can
                // add generator destinations, but cannot make a purely local
                // source emit a generator.
                for (int color = 0; color < num_colors; ++color) {
                    const auto& wave = color_bins[static_cast<size_t>(
                        interior_start_loc + color)];
                    for (int64_t morton : wave) {
                        const auto* box = level.find_local_box(morton);
                        if (box == nullptr) continue;
                        for (int64_t neighbor : box->one_hop) {
                            if (level.find_local_box(neighbor) == nullptr) {
                                local_flush_after_wave[
                                    static_cast<size_t>(color)] = 1;
                                break;
                            }
                        }
                        if (local_flush_after_wave[
                                static_cast<size_t>(color)] != 0) {
                            break;
                        }
                    }
                }

                const auto decision_start = clock::now();
                MPI_Allreduce(
                    local_flush_after_wave.data(),
                    lazy_interior_flush_after_wave.data(),
                    num_colors, MPI_INT, MPI_MAX, level_comm);
                const auto decision_duration =
                    clock::now() - decision_start;
                level_reduction += decision_duration;
                if (profile_deep_level) {
                    const double decision_ms =
                        std::chrono::duration<double, std::milli>(
                            decision_duration).count();
                    deep_phase_ms[DEEP_TRANSPORT_WALL] += decision_ms;
                    deep_phase_ms[DEEP_TRANSPORT_MPI] += decision_ms;
                }
                lazy_interior_schedule_ready = true;
            };
            auto transport_color_updates =
                [&](const FactorizationMemoryDiagnosticCallback& diagnostic) {
                    return transport_and_apply_factor_updates_symmetric_onehop(
                        tree, current_level, kernel, pending_updates,
                        false, diagnostic, true);
                };
            auto refresh_installed_lazy_remote_state = [&]() {
                if (lazy_far_field_mode() != LazyFarFieldMode::LAZY) return;

                // A generator proves that the remote source has completed
                // elimination. Mirror its skeleton into an overlapping ghost
                // copy as well: ghost geometry is gathered once at setup,
                // whereas lazy generators arrive wave by wave.
                for (const auto& entry : level.generator_id_to_index) {
                    level.eliminated_boxes.insert(entry.first);
                    auto ghost_it = level.ghost_id_to_index.find(entry.first);
                    if (ghost_it != level.ghost_id_to_index.end()) {
                        auto& ghost = level.ghost_boxes[
                            static_cast<size_t>(ghost_it->second)];
                        const auto& generator = level.generator_boxes[
                            static_cast<size_t>(entry.second)];
                        ghost.skeleton_indices = generator.skeleton_indices;
                    }
                }

                // The preceding assisting exchange refreshes every remote
                // two-hop endpoint from its owner. A nonempty skeleton proves
                // that endpoint has been eliminated, even when this rank does
                // not need (and is not sent) that endpoint's generator.
                // Prefer current owner metadata over an overlapping setup-time
                // ghost so lazy temp2 slots and row selection agree.
                for (const auto& entry :
                     level.assisting_box_points_for_kernel_evaluation) {
                    const auto& assist = level.assisting_boxes[
                        static_cast<size_t>(entry.second)];
                    if (assist.skel_indices.empty()) continue;

                    level.eliminated_boxes.insert(entry.first);
                    auto ghost_it = level.ghost_id_to_index.find(entry.first);
                    if (ghost_it != level.ghost_id_to_index.end()) {
                        auto& ghost = level.ghost_boxes[
                            static_cast<size_t>(ghost_it->second)];
                        ghost.skeleton_indices = assist.skel_indices;
                        ghost.on_boundary = assist.on_boundary;
                    }
                }
            };
            int boundary_count = 0;
            bool to_store = true;
            const bool store_interior_wave = true;

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {

                const int  color_id_mod    = counter % num_colors;
                const bool is_interior     = (counter >= interior_start_loc);

                // ----------------------------------------------------------------
                // Communication / transport step (single-threaded)
                // ----------------------------------------------------------------
                // Eager Color can batch the interior sub-waves because an
                // interior box has no remote numerical state to publish. In
                // an unstructured tree, however, owner and assisting copies
                // can classify the same occupied box differently at a sparse
                // process boundary. Preserve the corrected lazy ordering for
                // every interior wave that can publish remote state. The
                // schedule is agreed once per level after the mandatory first
                // interior transport has completed its assisting exchange.
                bool transport_required =
                    !is_interior || counter == interior_start_loc;
                if (!transport_required &&
                    lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                    if (!lazy_interior_schedule_ready) {
                        throw std::runtime_error(
                            "lazy interior transport schedule is not ready");
                    }
                    const int previous_interior_color =
                        counter - interior_start_loc - 1;
                    transport_required =
                        lazy_interior_flush_after_wave[static_cast<size_t>(
                            previous_interior_color)] != 0;
                    if (!transport_required && !pending_updates_empty()) {
                        throw std::runtime_error(
                            "lazy interior transport schedule missed pending updates");
                    }
                }
                if (transport_required) {
                    if (memory_diagnostics.enabled()) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            "wave" + std::to_string(counter) + "_pre_transport",
                            h2_diag_pending_bytes(pending_updates));
                    }
                    FactorizationMemoryDiagnosticCallback transport_memory_diagnostic;
                    if (memory_diagnostics.enabled()) {
                        transport_memory_diagnostic =
                            [&](const char* phase, size_t additional_pending,
                                size_t communication) {
                                memory_diagnostics.record(
                                    tree, current_level, "color",
                                    "wave" + std::to_string(counter) +
                                        "_transport_" + phase,
                                    h2_diag_pending_bytes(pending_updates) +
                                        additional_pending,
                                    0, 0, communication);
                            };
                    }
                    const auto transport_wall_start = clock::now();
                    const auto comm_duration_raw =
                        transport_color_updates(transport_memory_diagnostic);
                    record_deep_phase(
                        DEEP_TRANSPORT_WALL, transport_wall_start);
                    if (profile_deep_level) {
                        deep_phase_ms[DEEP_TRANSPORT_MPI] +=
                            std::chrono::duration<double, std::milli>(
                                comm_duration_raw).count();
                    }
                    const auto post_transport_start = clock::now();
                    refresh_installed_lazy_remote_state();
                    level_data_exchange += comm_duration_raw;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                    record_deep_phase(
                        DEEP_POST_TRANSPORT, post_transport_start);
                    if (memory_diagnostics.enabled()) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            "wave" + std::to_string(counter) + "_post_transport",
                            h2_diag_pending_bytes(pending_updates));
                    }
                    auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                    if (print_detail && rank == level_print_rank) {
                        std::cout << "  Comm time: "
                                << comm_duration.count()
                                << " ms" << std::endl;
                    }
                }

                if (lazy_far_field_mode() == LazyFarFieldMode::LAZY &&
                    counter == interior_start_loc) {
                    build_lazy_interior_flush_schedule();
                }

                const auto& color_list = color_bins[static_cast<size_t>(counter)];
                auto mark_assisting_boxes_eliminated = [&]() {
                    for (const auto& kv :
                         level.assisting_box_points_for_kernel_evaluation) {
                        const int assisting_color = static_cast<int>(
                            kv.first & (num_colors - 1));
                        const bool assisting_is_boundary =
                            level.assisting_boxes[static_cast<size_t>(kv.second)]
                                .on_boundary;
                        const bool mark =
                            (!is_interior && assisting_is_boundary &&
                             assisting_color == color_id_mod) ||
                            (is_interior && !assisting_is_boundary &&
                             assisting_color == color_id_mod);
                        // A lazy remote source is usable only once its
                        // generator and elimination wave have arrived. Do not
                        // advance it from locally reconstructed boundary
                        // metadata, which can differ across sparse regions.
                        if (mark &&
                            lazy_far_field_mode() != LazyFarFieldMode::LAZY) {
                            level.eliminated_boxes.insert(kv.first);
                        }
                    }
                };

                if (print_detail && rank == level_print_rank) {
                    std::cout << "  Processing "
                            << (is_interior ? "interior sub-wave " : "boundary color ")
                            << color_id_mod
                            << " (" << color_list.size() << " boxes)..." << std::endl;
                }

                if (color_list.empty()) {
                    // All ranks retain the same wave schedule and have already
                    // participated in this wave's transport.  Avoid numerical
                    // work while still advancing remote assisting-box state.
                    mark_assisting_boxes_eliminated();
                    continue;
                }

                const auto wave_setup_start = clock::now();
                const bool enable_deferred_xnn =
                    to_store && (store_interior_wave || !is_interior);
                const bool use_owner_deferred_xnn =
                    enable_deferred_xnn && is_symmetric && !is_hermitian;

                std::unordered_set<int64_t> wave_box_set;
                wave_box_set.insert(color_list.begin(), color_list.end());

                for (int64_t morton_idx : color_list) {
                    BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                    if (box_ptr == nullptr) {
                        throw std::runtime_error(
                            "Morton index " + std::to_string(morton_idx) +
                            " not found at level " + std::to_string(current_level) +
                            " (counter: " + std::to_string(counter) + ")");
                    }

                    for (int64_t neighbor_morton : box_ptr->one_hop) {
                        if (neighbor_morton != morton_idx &&
                            wave_box_set.count(neighbor_morton) != 0) {
                            std::ostringstream oss;
                            oss << "Wave construction error: one-hop conflict"
                                << " level=" << current_level
                                << " color_id_mod=" << color_id_mod
                                << " counter=" << counter
                                << " box=" << morton_idx
                                << " neighbor=" << neighbor_morton;
                            throw std::runtime_error(oss.str());
                        }
                    }
                }

                const int max_threads = std::max(1, omp_get_max_threads());
                std::vector<PendingFactorUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_threads));
                std::vector<std::vector<int64_t>> thread_xnn_candidate_boxes;
                std::vector<size_t> candidate_box_offsets;
                std::vector<int64_t> wave_xnn_candidate_boxes;
                std::vector<std::vector<DeferredXnnTargetKey>> wave_xnn_mirror_targets;
                std::vector<int> thread_boundary_counts(static_cast<size_t>(max_threads), 0);

                if (use_owner_deferred_xnn) {
                    thread_xnn_candidate_boxes.resize(static_cast<size_t>(max_threads));
                    candidate_box_offsets.resize(thread_xnn_candidate_boxes.size() + 1, 0);
                }

                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};
                size_t wave_scratch_bytes = 0;

                const int wave_team = std::max(1, omp_get_max_threads());
                const int wave_split = split_threads_for(
                    static_cast<int64_t>(color_list.size()), wave_team);
                record_deep_phase(DEEP_WAVE_SETUP, wave_setup_start);
                const auto primary_start = clock::now();

                #pragma omp parallel default(shared)
                {
                    const int tid = omp_get_thread_num();
                    FactorizationThreadScratch<CoordType, DataType> scratch;
                    scratch.split_threads = wave_split;

                    #pragma omp for schedule(dynamic)
                    for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                        if (wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx = color_list[static_cast<size_t>(bi)];
                            BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                            if (box_ptr == nullptr) {
                                throw std::runtime_error(
                                    "Morton index " + std::to_string(morton_idx) +
                                    " not found at level " + std::to_string(current_level) +
                                    " (counter: " + std::to_string(counter) + ")");
                            }

                            auto& box = *box_ptr;

                            // Keep the geometric owner slot in the Color wave
                            // so local and remote elimination state advances in
                            // lockstep, but do no numerical work for an empty
                            // box. Its filtered neighbor lists are also empty.
                            if (!occupied_topology.contains(
                                    current_level, morton_idx)) {
                                box.skeleton_indices.clear();
                                box.redundant_indices.clear();
                                continue;
                            }

                            if (use_streamed_level) {
                                gather_id_target_streamed(
                                    tree, &box, level, kernel, scratch,
                                    box.on_boundary);
                            } else {
                                scratch.streamed_sketch_valid = false;
                                gather_id_workspace(
                                    tree,
                                    &box, level, kernel, tolerance,
                                    unit_proxy_points.data(), num_proxy,
                                    proxy_radius, is_symmetric,
                                    scratch.workspace, scratch.workspace_rows,
                                    scratch.workspace_cols, 0,
                                    box.on_boundary);
                            }

                            thread_boundary_counts[static_cast<size_t>(tid)] += box.on_boundary;

                            compute_and_modify(dimension,
                                &box, level, kernel,
                                scratch,
                                tolerance, use_sketch, is_symmetric, is_hermitian,
                                &thread_pending[tid],
                                factorization_method, use_owner_deferred_xnn
                            );

                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }

                    if (memory_diagnostics.enabled()) {
                        const size_t thread_scratch_bytes =
                            h2_diag_factorization_scratch_bytes(scratch);
                        #pragma omp atomic update
                        wave_scratch_bytes += thread_scratch_bytes;
                        #pragma omp barrier
                        #pragma omp single
                        {
                            size_t all_pending_bytes =
                                h2_diag_pending_bytes(pending_updates);
                            for (const auto& thread_updates : thread_pending) {
                                all_pending_bytes +=
                                    h2_diag_pending_bytes(thread_updates);
                            }
                            memory_diagnostics.record(
                                tree, current_level, "color",
                                "wave" + std::to_string(counter) + "_compute",
                                all_pending_bytes, wave_scratch_bytes);
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }

                for (int local_boundary_count : thread_boundary_counts) {
                    boundary_count += local_boundary_count;
                }
                record_deep_phase(DEEP_PRIMARY, primary_start);

                if (use_owner_deferred_xnn) {
                    const auto owner_prep_start = clock::now();
                    for (int64_t morton_idx : color_list) {
                        level.eliminated_boxes.insert(morton_idx);
                        level.elimination_wave[morton_idx] =
                            static_cast<int32_t>(counter);
                    }

                    slice_far_field_blocks(level, is_symmetric, is_hermitian);

                    std::exception_ptr candidate_exception;
                    std::mutex candidate_exception_mutex;
                    std::atomic<bool> candidate_failed{false};

                    #pragma omp parallel default(shared) if (color_list.size() > 1)
                    {
                        const int tid = omp_get_thread_num();
                        auto& local_candidates =
                            thread_xnn_candidate_boxes[static_cast<size_t>(tid)];

                        #pragma omp for schedule(static)
                        for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                            if (candidate_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t morton_idx = color_list[static_cast<size_t>(bi)];
                                if (!occupied_topology.contains(
                                        current_level, morton_idx)) {
                                    continue;
                                }
                                BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                                if (box_ptr == nullptr) {
                                    continue;
                                }

                                collect_owner_deferred_xnn_candidates_for_source_box(
                                    box_ptr,
                                    level,
                                    wave_box_set,
                                    local_candidates);
                            } catch (...) {
                                if (!candidate_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(candidate_exception_mutex);
                                    candidate_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (candidate_exception) {
                        std::rethrow_exception(candidate_exception);
                    }

                    for (size_t tid = 0; tid < thread_xnn_candidate_boxes.size(); ++tid) {
                        candidate_box_offsets[tid + 1] =
                            candidate_box_offsets[tid] + thread_xnn_candidate_boxes[tid].size();
                    }

                    wave_xnn_candidate_boxes.resize(candidate_box_offsets.back());
                    for (int tid_idx = 0; tid_idx < static_cast<int>(thread_xnn_candidate_boxes.size()); ++tid_idx) {
                        size_t out_idx = candidate_box_offsets[static_cast<size_t>(tid_idx)];
                        for (int64_t candidate_morton : thread_xnn_candidate_boxes[static_cast<size_t>(tid_idx)]) {
                            wave_xnn_candidate_boxes[out_idx++] = candidate_morton;
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
                    wave_xnn_mirror_targets.resize(wave_xnn_candidate_boxes.size());

                    std::exception_ptr owner_exception;
                    std::mutex owner_exception_mutex;
                    std::atomic<bool> owner_failed{false};
                    size_t owner_scratch_bytes = 0;

                    const int owner_team = std::max(1, omp_get_max_threads());
                    // A candidate replay can touch every pair in the
                    // candidate's two-hop neighborhood. Use a modulo-three
                    // spatial coloring so candidates executed together have
                    // disjoint writable neighborhoods.
                    const int owner_color_count =
                        (dimension == 3) ? 27 : 9;
                    std::vector<std::vector<int64_t>> owner_color_bins(
                        static_cast<size_t>(owner_color_count));
                    for (int64_t idx = 0;
                         idx < static_cast<int64_t>(
                             wave_xnn_candidate_boxes.size());
                         ++idx) {
                        uint32_t x = 0;
                        uint32_t y = 0;
                        uint32_t z = 0;
                        morton::decode_nd(
                            dimension,
                            wave_xnn_candidate_boxes[
                                static_cast<size_t>(idx)],
                            x, y, z);
                        const int owner_color =
                            static_cast<int>(x % 3) +
                            3 * static_cast<int>(y % 3) +
                            (dimension == 3
                                ? 9 * static_cast<int>(z % 3)
                                : 0);
                        owner_color_bins[
                            static_cast<size_t>(owner_color)].push_back(idx);
                    }

                    record_deep_phase(DEEP_OWNER_PREP, owner_prep_start);
                    const auto owner_replay_start = clock::now();
                    for (const auto& owner_color_bin : owner_color_bins) {
                        if (owner_color_bin.empty()) continue;
                        const int owner_split = split_threads_for(
                            static_cast<int64_t>(owner_color_bin.size()),
                            owner_team);
                        #pragma omp parallel default(shared)
                        {
                            const int tid = omp_get_thread_num();
                            DeferredXnnOwnerScratch<DataType> scratch;
                            scratch.split_threads = owner_split;

                            #pragma omp for schedule(dynamic)
                            for (int64_t bin_slot = 0;
                                 bin_slot < static_cast<int64_t>(
                                     owner_color_bin.size());
                                 ++bin_slot) {
                                if (owner_failed.load(
                                        std::memory_order_relaxed)) {
                                    continue;
                                }

                                try {
                                    const int64_t idx = owner_color_bin[
                                        static_cast<size_t>(bin_slot)];
                                    const int64_t candidate_morton =
                                        wave_xnn_candidate_boxes[
                                            static_cast<size_t>(idx)];

                                    // Deferred store=true replay:
                                    //   local-local   -> local owner replay + local mirror
                                    //   local-remote  -> local replay + remote transport
                                    //   remote-remote -> remote transport only
                                    apply_unstructured_owner_deferred_xnn_updates_for_candidate_box(
                                        candidate_morton,
                                        level,
                                        kernel,
                                        wave_box_set,
                                        scratch,
                                        wave_xnn_mirror_targets[
                                            static_cast<size_t>(idx)],
                                        &thread_pending[
                                            static_cast<size_t>(tid)]);
                                } catch (...) {
                                    if (!owner_failed.exchange(
                                            true,
                                            std::memory_order_relaxed)) {
                                        std::lock_guard<std::mutex> lock(
                                            owner_exception_mutex);
                                        owner_exception =
                                            std::current_exception();
                                    }
                                }
                            }

                            if (memory_diagnostics.enabled()) {
                                const size_t thread_scratch_bytes =
                                    h2_diag_deferred_owner_scratch_bytes(
                                        scratch);
                                #pragma omp atomic update
                                owner_scratch_bytes += thread_scratch_bytes;
                                #pragma omp barrier
                                #pragma omp single
                                {
                                    size_t all_pending_bytes =
                                        h2_diag_pending_bytes(pending_updates);
                                    for (const auto& thread_updates :
                                         thread_pending) {
                                        all_pending_bytes +=
                                            h2_diag_pending_bytes(
                                                thread_updates);
                                    }
                                    memory_diagnostics.record(
                                        tree, current_level, "color",
                                        "wave" + std::to_string(counter) +
                                            "_owner",
                                        all_pending_bytes,
                                        owner_scratch_bytes);
                                }
                            }
                        }
                    }

                    record_deep_phase(
                        DEEP_OWNER_REPLAY, owner_replay_start);
                    if (owner_exception) {
                        std::rethrow_exception(owner_exception);
                    }

                    const auto mirror_start = clock::now();
                    std::exception_ptr mirror_exception;
                    std::mutex mirror_exception_mutex;
                    std::atomic<bool> mirror_failed{false};

                    #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                    {
                        #pragma omp for schedule(static)
                        for (int64_t idx = 0; idx < static_cast<int64_t>(wave_xnn_candidate_boxes.size()); ++idx) {
                            if (mirror_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t candidate_morton =
                                    wave_xnn_candidate_boxes[static_cast<size_t>(idx)];
                                BoxData<CoordType, DataType>* candidate_box =
                                    level.find_local_box(candidate_morton);
                                if (candidate_box == nullptr) {
                                    candidate_box = level.find_ghost_box(candidate_morton);
                                }
                                if (candidate_box == nullptr) {
                                    continue;
                                }

                                apply_symmetric_owner_deferred_xnn_updates_for_candidate_box(
                                    candidate_box,
                                    wave_xnn_mirror_targets[static_cast<size_t>(idx)],
                                    level);
                            } catch (...) {
                                if (!mirror_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(mirror_exception_mutex);
                                    mirror_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (mirror_exception) {
                        std::rethrow_exception(mirror_exception);
                    }
                    record_deep_phase(DEEP_MIRROR, mirror_start);

                    const auto share_start = clock::now();
                    share_symmetric_level_edges(level);
                    record_deep_phase(DEEP_SHARE, share_start);

                    const auto generator_start = clock::now();
                    // Package the generator while X_NR still contains the
                    // original one-hop coupling. Finalization replaces X_NR
                    // with temp2, which remains the representation retained
                    // by local and installed generator boxes.
                    if (lazy_far_field_mode() == LazyFarFieldMode::LAZY) {
                        std::vector<int64_t> occupied_color_list;
                        occupied_color_list.reserve(color_list.size());
                        for (int64_t morton_idx : color_list) {
                            if (occupied_topology.contains(
                                    current_level, morton_idx)) {
                                occupied_color_list.push_back(morton_idx);
                            }
                        }
                        emit_lazy_generators(
                            level, occupied_color_list, counter,
                            pending_updates);
                    }
                    record_deep_phase(DEEP_GENERATOR, generator_start);

                    const auto finalize_start = clock::now();
                    std::exception_ptr finalize_exception;
                    std::mutex finalize_exception_mutex;
                    std::atomic<bool> finalize_failed{false};

                    #pragma omp parallel default(shared) if (color_list.size() > 1)
                    {
                        #pragma omp for schedule(static)
                        for (int64_t bi = 0; bi < static_cast<int64_t>(color_list.size()); ++bi) {
                            if (finalize_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t source_morton = color_list[static_cast<size_t>(bi)];
                                if (!occupied_topology.contains(
                                        current_level, source_morton)) {
                                    continue;
                                }
                                BoxData<CoordType, DataType>* source_box = level.find_local_box(source_morton);
                                if (source_box == nullptr) {
                                    continue;
                                }

                                finalize_deferred_xnn_source_box(source_box);
                            } catch (...) {
                                if (!finalize_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(finalize_exception_mutex);
                                    finalize_exception = std::current_exception();
                                }
                            }
                        }
                    }

                    if (finalize_exception) {
                        std::rethrow_exception(finalize_exception);
                    }
                    record_deep_phase(DEEP_FINALIZE, finalize_start);

                } else {
                    slice_far_field_blocks(level, is_symmetric, is_hermitian);
                }
                const auto bookkeeping_start = clock::now();
                for (int t = 0; t < max_threads; ++t) {
                    merge_pending(pending_updates, thread_pending[t]);
                    clear_pending_factor_updates_memory(thread_pending[t]);
                }

                if (!use_owner_deferred_xnn) {
                    for (size_t bi = 0; bi < color_list.size(); ++bi) {
                        const int64_t morton_idx = color_list[bi];
                        level.eliminated_boxes.insert(morton_idx);
                        level.elimination_wave[morton_idx] =
                            static_cast<int32_t>(counter);
                    }
                }

                // Advance remote assisting boxes after local numerical work.
                mark_assisting_boxes_eliminated();
                record_deep_phase(DEEP_BOOKKEEPING, bookkeeping_start);

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, "color",
                        "wave" + std::to_string(counter) + "_complete",
                        h2_diag_pending_bytes(pending_updates));
                }

            }

            // Keep the final flush outside the wave body so empty trailing bins
            // cannot strand updates produced by the last nonempty wave.
            if (memory_diagnostics.enabled()) {
                memory_diagnostics.record(
                    tree, current_level, "color", "final_pre_transport",
                    h2_diag_pending_bytes(pending_updates));
            }
            FactorizationMemoryDiagnosticCallback final_transport_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                final_transport_memory_diagnostic =
                    [&](const char* phase, size_t additional_pending,
                        size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "color",
                            std::string("final_transport_") + phase,
                            h2_diag_pending_bytes(pending_updates) +
                                additional_pending,
                            0, 0, communication);
                    };
            }
            const auto final_transport_wall_start = clock::now();
            const auto final_comm_duration =
                transport_color_updates(final_transport_memory_diagnostic);
            record_deep_phase(
                DEEP_TRANSPORT_WALL, final_transport_wall_start);
            if (profile_deep_level) {
                deep_phase_ms[DEEP_TRANSPORT_MPI] +=
                    std::chrono::duration<double, std::milli>(
                        final_comm_duration).count();
            }
            const auto final_post_transport_start = clock::now();
            refresh_installed_lazy_remote_state();
            level_data_exchange += final_comm_duration;
            update_neighbor_slicing_for_level(level, is_symmetric);
            record_deep_phase(
                DEEP_POST_TRANSPORT, final_post_transport_start);
            if (memory_diagnostics.enabled()) {
                memory_diagnostics.record(
                    tree, current_level, "color", "final_post_transport",
                    h2_diag_pending_bytes(pending_updates));
            }

            const auto occupied_indices = occupied_local_indices(level);
            for (int64_t local_index : occupied_indices) {
                const auto& box =
                    level.local_boxes[static_cast<size_t>(local_index)];
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
                local_max_skel = std::max<int64_t>(
                    local_max_skel, box.skeleton_indices.size());
            }
            }

            const auto elim_end = std::chrono::high_resolution_clock::now();
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                elim_end - elim_start);
            if (profile_deep_level) {
                deep_elimination_local_ms =
                    std::chrono::duration<double, std::milli>(
                        elim_end - elim_start).count();
            }
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_elimination");

        if (current_level > 1 && level.is_process_active) {
            if (profile_deep_level && !use_CA_level) {
                double exclusive_ms = 0.0;
                for (size_t phase = 0; phase < DEEP_PHASE_COUNT; ++phase) {
                    if (phase != DEEP_TRANSPORT_MPI &&
                        phase != DEEP_OTHER) {
                        exclusive_ms += deep_phase_ms[phase];
                    }
                }
                deep_phase_ms[DEEP_OTHER] = std::max(
                    0.0, deep_elimination_local_ms - exclusive_ms);

                std::array<double, DEEP_PHASE_COUNT> min_phase_ms{};
                std::array<double, DEEP_PHASE_COUNT> max_phase_ms{};
                std::array<double, DEEP_PHASE_COUNT> sum_phase_ms{};
                MPI_Reduce(
                    deep_phase_ms.data(), min_phase_ms.data(),
                    DEEP_PHASE_COUNT, MPI_DOUBLE, MPI_MIN, 0, level_comm);
                MPI_Reduce(
                    deep_phase_ms.data(), max_phase_ms.data(),
                    DEEP_PHASE_COUNT, MPI_DOUBLE, MPI_MAX, 0, level_comm);
                MPI_Reduce(
                    deep_phase_ms.data(), sum_phase_ms.data(),
                    DEEP_PHASE_COUNT, MPI_DOUBLE, MPI_SUM, 0, level_comm);

                if (rank == level_print_rank) {
                    std::cout
                        << "  Temporary color_unstructured phase profile"
                        << " (min / average / max across "
                        << level.num_active_processes << " active ranks):"
                        << std::endl;
                    for (size_t phase = 0;
                         phase < DEEP_PHASE_COUNT; ++phase) {
                        const double average =
                            sum_phase_ms[phase] /
                            static_cast<double>(
                                level.num_active_processes);
                        std::printf(
                            "    %-34s %10.3f / %10.3f / %10.3f ms\n",
                            deep_phase_names[phase],
                            min_phase_ms[phase], average,
                            max_phase_ms[phase]);
                    }
                    std::cout
                        << "    (MPI inside transport is a subset of"
                        << " transport wall.)" << std::endl;
                }
            }

            double min_elim_ms = 0.0;
            double max_elim_ms = 0.0;
            reduce_active_duration_bounds_ms(
                level_comm,
                0,
                true,
                elim_duration,
                min_elim_ms,
                max_elim_ms);

            const int64_t local_occupied_count =
                static_cast<int64_t>(level.boundary_id.size()) +
                static_cast<int64_t>(level.interior_id.size());
            const int64_t local_compression_counts[3] = {
                total_skeleton,
                total_redundant,
                local_occupied_count
            };
            int64_t global_compression_counts[3] = {0, 0, 0};
            MPI_Reduce(
                local_compression_counts,
                global_compression_counts,
                3,
                MPI_INT64_T,
                MPI_SUM,
                0,
                level_comm);

            if (print_detail && rank == level_print_rank) {
                std::cout << "  Elimination time: shortest=" << std::llround(min_elim_ms)
                          << " ms, longest=" << std::llround(max_elim_ms) << " ms" << std::endl;

                const int64_t global_skeleton = global_compression_counts[0];
                const int64_t global_redundant = global_compression_counts[1];
                const int64_t global_occupied_count =
                    global_compression_counts[2];
                const int64_t global_original =
                    global_skeleton + global_redundant;
                const double compression = global_original > 0
                    ? static_cast<double>(global_skeleton) /
                        static_cast<double>(global_original)
                    : 0.0;
                std::cout << "  Compression ratio: " << compression
                          << " (" << (compression * 100) << "%)" << std::endl;
                std::cout << "  Average skeleton size: "
                          << (global_occupied_count > 0
                              ? static_cast<double>(global_skeleton) /
                                  static_cast<double>(global_occupied_count)
                              : 0.0)
                          << std::endl;
            }
        }

        if (current_level > 1) {
            const bool owner_engine_level =
                use_CA_level && ca_owner_component == 3 &&
                level.num_active_processes > 1;
            owner_solve_record_level(
                owner_schedule, current_level, owner_engine_level);
            if (level.is_process_active) {
                if (print_detail && rank == level_print_rank &&
                    owner_schedule.active) {
                    owner_schedule_print_timing(
                        owner_schedule, current_level);
                }
                level.pair_interest_filter = nullptr;
                level.sketch_eliminated_filter = nullptr;
            }
        }

        if (level.is_process_active && staged_overlap_scheduling) {
            if (level.staged_pending != nullptr &&
                (!level.staged_pending->deltas.empty() ||
                 !level.staged_pending->mirrors.empty())) {
                throw std::runtime_error(
                    "staged halo overlap: pending updates not drained at level end");
            }
            level.staged_overlap_on = false;
            level.staged_unarrived.clear();
            level.staged_stage_map.clear();
            level.staged_pending.reset();
        }
        if (level.is_process_active && staged_state.active) {
            if (print_detail && rank == level_print_rank) {
                staged_halo_print_timing(
                    staged_state, current_level);
            }
            staged_halo_finish(staged_state);
        }

        // ===== Step 3: Gather assisting boxes post-elimination =====
        if (current_level > 1 && level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            FactorizationMemoryDiagnosticCallback assisting_memory_diagnostic;
            if (memory_diagnostics.enabled()) {
                assisting_memory_diagnostic =
                    [&](const char* phase, size_t pending, size_t communication) {
                        memory_diagnostics.record(
                            tree, current_level, "CA",
                            std::string("post_elimination_assisting_") + phase,
                            pending, 0, 0, communication);
                    };
            }
            gather_CA_assisting_boxes_factorization(
                tree, current_level, level_comm, 412,
                assisting_memory_diagnostic);
            const auto gather_duration = clock::now() - segment_start;
            level_data_exchange += gather_duration;
            update_neighbor_slicing_for_level(
                level, is_symmetric,
                /*use_received_assisting_skeletons=*/true);
            if (print_detail && rank == level_print_rank) {
                std::cout << "  CA post-elimination assisting gather time: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 gather_duration).count()
                          << " ms" << std::endl;
            }
        } else {
            clear_ghosts(level);
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_assisting_gather");

        if (level.is_process_active && !use_CA_level) {
            // The shared dense parent assembler indexes all geometric
            // children. Supply zero-sized records only for its empty remote
            // children; occupied remote children must already carry the
            // skeleton metadata exchanged by the Color waves.
            add_empty_parent_transition_assisting_slots(
                tree, current_level, occupied_topology);
        }




        // ===== Step 4: Build parent level interactions =====

        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1 && level.is_process_active) {
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Setting all DOFs as skeleton (no redundant DOFs at level 1)" << std::endl;
            }

            for (int64_t local_index : occupied_local_indices(level)) {
                auto& box =
                    level.local_boxes[static_cast<size_t>(local_index)];
                // All points are skeleton (no elimination at level 1)
                box.skeleton_indices.resize(box.num_points);
                for (int64_t i = 0; i < box.num_points; ++i) {
                    box.skeleton_indices[i] = i;
                }

                // No redundant DOFs
                box.redundant_indices.clear();
            }
        }

        auto transition_start = std::chrono::high_resolution_clock::now();
        // if(rank == 1)
        // {
        //     build_parent_level_interactions<CoordType, DataType, KernelType>(
        //         level,          // child_level
        //         dimension,
        //         is_symmetric,
        //         is_hermitian,
        //         kernel,
        //         tree->global_bounds
        //     );

        // }
        // MPI_Barrier(tree->comm);
        // exit(0);
        std::vector<BoxData<CoordType, DataType>> parent_boxes;
        if (level.is_process_active) {
            parent_boxes =
                build_parent_level_interactions_unstructured<
                    CoordType, DataType, KernelType>(
                    level,
                    tree->levels[current_level - 1],
                    dimension,
                    is_symmetric,
                    is_hermitian,
                    kernel,
                    tree->global_bounds,
                    occupied_topology);
        }

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_build", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes));
        }




        auto transition_end = std::chrono::high_resolution_clock::now();
        auto transition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(transition_end - transition_start);

        if (print_detail && rank == level_print_rank) {
            std::cout << "  Level transition time: " << transition_duration.count() << " ms" << std::endl;
            std::cout << "  Parent boxes created: " << parent_boxes.size() << std::endl;
        }

        // ===== Step 5: Handle process reduction =====

        const bool reduction_occurred =
            (parent_level.num_active_processes != level.num_active_processes);
        const bool keep_local_parent_boxes =
            level.is_process_active && (level.parent_level_owner == rank);
        const bool sends_parent_boxes =
            reduction_occurred && level.is_process_active && level.parent_level_owner != rank;

        std::vector<char> send_buffer;
        int64_t send_buffer_size = 0;
        if (sends_parent_boxes) {
            // Pre-pack parent boxes before the synchronization point so the
            // communication timer excludes sender-side serialization work.
            send_buffer = serialize_boxes(parent_boxes);
            send_buffer_size = static_cast<int64_t>(send_buffer.size());
        }

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_pack", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes),
                h2_diag_vector_bytes(send_buffer));
        }

        if (reduction_occurred) {
            // Synchronize after parent-box construction/packing so the reduction
            // communication timer does not count time spent waiting for slower
            // ranks that are still preparing their payloads.
            MPI_Barrier(transition_comm);
        }

        if (!reduction_occurred) {
            if (parent_level.is_process_active) {
                parent_level.local_boxes = std::move(parent_boxes);
            }
        } else if (parent_level.is_process_active) {
            std::vector<BoxData<CoordType, DataType>> all_parent_boxes;

            for (int child_rank : parent_level.children_senders) {
                if (keep_local_parent_boxes && child_rank == rank) {
                    all_parent_boxes.insert(
                        all_parent_boxes.end(),
                        std::make_move_iterator(parent_boxes.begin()),
                        std::make_move_iterator(parent_boxes.end())
                    );
                    continue;
                }

                std::vector<char> recv_buffer;
                int64_t buffer_size = 0;
                MPI_Status status;
                segment_start = clock::now();
                MPI_Recv(&buffer_size, 1, MPI_INT64_T, child_rank, 0, tree->comm, &status);

                recv_buffer.resize(buffer_size);
                MPI_Recv_large(recv_buffer.data(), buffer_size, MPI_CHAR, child_rank, 1, tree->comm, &status);
                level_reduction += (clock::now() - segment_start);

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, use_CA_level ? "CA" : "color",
                        "parent_receive_buffer", 0, 0,
                        h2_diag_box_vector_bytes(parent_boxes) +
                            h2_diag_box_vector_bytes(all_parent_boxes),
                        h2_diag_vector_bytes(recv_buffer));
                }

                std::vector<BoxData<CoordType, DataType>> child_parent_boxes =
                    deserialize_boxes<CoordType, DataType>(recv_buffer);

                if (memory_diagnostics.enabled()) {
                    memory_diagnostics.record(
                        tree, current_level, use_CA_level ? "CA" : "color",
                        "parent_deserialized", 0, 0,
                        h2_diag_box_vector_bytes(parent_boxes) +
                            h2_diag_box_vector_bytes(all_parent_boxes) +
                            h2_diag_box_vector_bytes(child_parent_boxes),
                        h2_diag_vector_bytes(recv_buffer));
                }

                all_parent_boxes.insert(
                    all_parent_boxes.end(),
                    std::make_move_iterator(child_parent_boxes.begin()),
                    std::make_move_iterator(child_parent_boxes.end())
                );
            }

            parent_level.local_boxes = std::move(all_parent_boxes);
        }

        if (parent_level.is_process_active &&
            is_symmetric && !is_hermitian) {
            share_symmetric_level_edges(parent_level);
        }

        if (sends_parent_boxes) {
            // This process does NOT own the parent - send to parent_level_owner
            segment_start = clock::now();
            MPI_Send(&send_buffer_size, 1, MPI_INT64_T, level.parent_level_owner, 0, tree->comm);

            // Send data
            MPI_Send_large(send_buffer.data(), send_buffer_size, MPI_CHAR, level.parent_level_owner, 1, tree->comm);
            level_reduction += (clock::now() - segment_start);

            // This process will become inactive at parent level
            parent_level.local_boxes.clear();

            if (print_trace) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size()
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }

        if (memory_diagnostics.enabled()) {
            memory_diagnostics.record(
                tree, current_level, use_CA_level ? "CA" : "color",
                "post_parent_redistribution", 0, 0,
                h2_diag_box_vector_bytes(parent_boxes),
                h2_diag_vector_bytes(send_buffer));
        }

        // Clear modified interaction matrices to free memory
        if (level.is_process_active) {
            clear_modified_interaction_matrices(level, use_CA_level);
        }

        memory_diagnostics.record(
            tree, current_level, use_CA_level ? "CA" : "color",
            "post_level_clear");

        // Teardown
        if (level.is_process_active) {
            for (auto& [morton, lock] : level.box_locks) {
                omp_destroy_lock(lock);
                delete lock;
            }
        }

        {
            double local_deltas[2] = {
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(level_data_exchange).count(),
                std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(level_reduction).count()
            };
            double max_deltas[2] = {0.0, 0.0};
            int transition_rank = -1;
            MPI_Comm_rank(transition_comm, &transition_rank);
            MPI_Reduce(
                local_deltas,
                max_deltas,
                2,
                MPI_DOUBLE,
                MPI_MAX,
                0,
                transition_comm);
            if (transition_rank == 0) {
                total_data_exchange_time += std::chrono::duration_cast<clock::duration>(
                    std::chrono::duration<double, std::milli>(max_deltas[0]));
                total_reduction_time += std::chrono::duration_cast<clock::duration>(
                    std::chrono::duration<double, std::milli>(max_deltas[1]));
            }
        }

        auto level_end = std::chrono::high_resolution_clock::now();
        auto level_duration = std::chrono::duration_cast<std::chrono::milliseconds>(level_end - level_start);

        if (print_detail && rank == level_print_rank) {
            std::cout << "  Total level time: " << level_duration.count() << " ms" << std::endl;
        }

        if (dynamic_threading.enabled &&
            level.is_process_active && !parent_level.is_process_active) {
            park_inactive_rank_on_service_cpu(dynamic_threading);
        }

    }


    if (out_rankmax) *out_rankmax = local_max_skel;


    // ===== Special handling for level 0 (root) =====

    const int root_print_rank = smallest_active_rank(tree->levels[0]);
    if (print_detail && rank == root_print_rank) {
        std::cout << "\n===== Level 0 (Root) =====" << std::endl;
    }

    auto& root_level = tree->levels[0];

    if (!root_level.is_process_active || root_level.local_boxes.empty()) {
        if (print_trace) {
            std::cout << "  Process " << rank << " has no boxes at root level" << std::endl;
        }
    } else {
        if (root_level.local_boxes.size() != 1) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Expected exactly 1 box at root level, got " +
                std::to_string(root_level.local_boxes.size()));
        }

        auto& root_box = root_level.local_boxes[0];

        if (print_detail && rank == root_print_rank) {
            std::cout << "  Root box points: " << root_box.num_points << std::endl;
        }

        // At level 0, the assembled matrix is just the schur complement
        if (!root_box.schur_complement.is_allocated()) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Root box schur complement not allocated");
        }

        int64_t n = root_box.schur_complement.rows;

        if (print_detail && rank == root_print_rank) {
            std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
        }

        // Factorize the root schur complement for diagonal solve
        root_box.X_RR.allocate(n, n, MatrixStorage<DataType>::FULL);

        if (factorization_method == FactorizationMethod::CHOLESKY) {
            root_box.X_RR_pivots.clear();
            // Copy schur complement to X_RR
            root_box.X_RR.data = root_box.schur_complement.data;

            // Perform Cholesky factorization in-place
            char uplo = 'L';
            int nn = n;
            int info = 0;

            if constexpr (std::is_same_v<DataType, double>) {
                dpotrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn, &info);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zsychol_(&uplo, &nn, root_box.X_RR.data.data(), &nn, &info);
            }

            if (info != 0) {
                throw std::runtime_error(
                    "hierarchical_factorization_parallel: Cholesky factorization of root failed at pivot " +
                    std::to_string(info));
            }

            root_box.X_RR.format = MatrixStorage<DataType>::CHOLESKY_L;

            if (print_detail && rank == root_print_rank) {
                std::cout << "  ✓ Root Cholesky factorization complete" << std::endl;
            }

        } else if (factorization_method == FactorizationMethod::LU) {
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR_pivots.resize(static_cast<size_t>(n));

            int nn = n;
            int info = 0;

            if constexpr (std::is_same_v<DataType, double>) {
                dgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), &info);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), &info);
            }

            if (info != 0) {
                throw std::runtime_error(
                    "hierarchical_factorization_parallel: LU factorization of root failed at pivot " +
                    std::to_string(info));
            }

            root_box.X_RR.format = MatrixStorage<DataType>::LU_FACTORED;

            if (print_detail && rank == root_print_rank) {
                std::cout << "  ✓ Root LU factorization complete" << std::endl;
            }

        } else if (factorization_method == FactorizationMethod::BUNCH_KAUFMAN) {
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR_pivots.resize(static_cast<size_t>(n));

            char uplo = 'L';
            int nn = n;
            int lwork = -1;
            int info = 0;
            std::vector<DataType> work(1);
            sytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                   root_box.X_RR_pivots.data(), work.data(), &lwork, &info);
            if (info != 0) {
                throw std::runtime_error(
                    "hierarchical_factorization_parallel: Bunch-Kaufman root workspace query failed with INFO = " +
                    std::to_string(info));
            }
            lwork = std::max(1, static_cast<int>(std::real(work[0])));
            work.resize(static_cast<size_t>(lwork));
            sytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                   root_box.X_RR_pivots.data(), work.data(), &lwork, &info);

            if (info != 0) {
                throw std::runtime_error(
                    "hierarchical_factorization_parallel: Bunch-Kaufman factorization of root failed at pivot " +
                    std::to_string(info));
            }

            root_box.X_RR.format = MatrixStorage<DataType>::BUNCH_KAUFMAN;

            if (print_detail && rank == root_print_rank) {
                std::cout << "  ✓ Root Bunch-Kaufman factorization complete" << std::endl;
            }

        } else {
            // No factorization: just copy schur complement to X_RR
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR.format = MatrixStorage<DataType>::FULL;
            root_box.X_RR_pivots.clear();

            if (print_detail && rank == root_print_rank) {
                std::cout << "  ✓ Root matrix copied (no factorization)" << std::endl;
            }
        }

        // Mark root as skeleton only (no redundant DOFs at this level)
        root_box.skeleton_indices.resize(n);
        for (int64_t i = 0; i < n; ++i) {
            root_box.skeleton_indices[i] = i;
        }
        root_box.redundant_indices.clear();
    }

    memory_diagnostics.record(tree, 0, "root", "post_root_factorization");

    passive_mpi_barrier(tree->comm);

    const double local_timing_ms[2] = {
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_data_exchange_time).count(),
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_reduction_time).count()
    };
    double global_timing_ms[2] = {0.0, 0.0};
    MPI_Reduce(
        local_timing_ms,
        global_timing_ms,
        2,
        MPI_DOUBLE,
        MPI_SUM,
        root_print_rank,
        tree->comm);

    if (print_summary && rank == root_print_rank) {
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - total_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ Hierarchical Factorization Complete" << std::endl;
        std::cout << "  total time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  data exchange communication time: " << std::llround(global_timing_ms[0]) << " ms" << std::endl;
        std::cout << "  process reduction communication time: " << std::llround(global_timing_ms[1]) << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    destroy_factorization_communicators(factorization_comms);
    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);


    size_t local_memory_usage = 0;
    size_t CA_halo_memory_usage = 0;
    for (int current_level = leaf_level; current_level >= 1; current_level--) {

        auto& level = tree->levels[current_level];
        for (auto& box : level.local_boxes) {
            local_memory_usage += calculate_box_data_size(box);
        }
        if (tree->level_uses_CA(current_level)) {
            for (const auto& box : level.ghost_boxes) {
                CA_halo_memory_usage += calculate_box_data_size(box);
            }
        }
    }
    if (print_detail) {
        printf("factorization memory usage on rank %d: %.10f GB local + %.10f GB CA halo\n",
               rank,
               local_memory_usage / (1024.0 * 1024.0 * 1024.0),
               CA_halo_memory_usage / (1024.0 * 1024.0 * 1024.0));
        fflush(stdout);
    }

    *memory_per_rank = local_memory_usage + CA_halo_memory_usage;

    memory_diagnostics.record(tree, 0, "retained", "factorization_complete");

    if (print_summary) {
        double logabsdet;
        DataType phase;
        hierarchical_logdet_parallel(tree, &logabsdet, &phase);
        if (rank == 0) {
            std::cout.precision(17);
            std::cout << "logdet: " << phase << " " << logabsdet << std::endl;
        }

        (void)h2_quick_verification_unstructured(tree, kernel);
    }

    memory_diagnostics.record(tree, 0, "verification", "post_quick_verification");
    memory_diagnostics.print(tree->comm, rank, size);
}


}  // namespace color_unstructured
}  // namespace butterfly
