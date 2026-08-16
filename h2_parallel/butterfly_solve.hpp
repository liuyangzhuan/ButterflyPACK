#pragma once

#include "butterfly_types.hpp"

namespace butterfly {
using namespace fmm;

// Gather local solution for both solve and mult
// solution is stored on the left_side regardless of whether it is a solve or mult operation 
template<typename CoordType, typename DataType> 
void gather_local_solution(ParallelTree<CoordType, DataType>* tree, 
  std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data, 
  DataType* x, int* Nloc) {
  const int leaf_level = tree->num_levels - 1;
  auto& leaf = tree->levels[leaf_level];

  int64_t g = 0;
  for (int64_t box_idx = 0; box_idx < leaf.num_boxes_local; ++box_idx) {
    const auto& solve_box = solve_data[leaf_level][box_idx];
    for (size_t i = 0; i < solve_box.left_side.size(); ++i) {   // == box.num_points
      x[g++] = solve_box.left_side[i];
    }
  }
  // sanity: g should equal *Nloc (== myseg on this rank)
  if (g != *Nloc) {
    throw std::runtime_error(
      "c_bpack_solve (format 7): solution DOF count mismatch — wrote " +
      std::to_string(g) + " values but Nloc = " + std::to_string(*Nloc) +
      " (leaf-box ordering inconsistent with construct_init's myseg)");
  }
}

template<typename CoordType, typename DataType>
void hierarchical_solve_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    int verbosity = 1) {
    
    int rank;
    MPI_Comm_rank(tree->comm, &rank);
    const bool print_summary = verbosity >= 0;
    const bool print_detail = verbosity >= 1;
    using clock = std::chrono::high_resolution_clock;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    auto get_data_start = clock::now();
    clock::duration communication_total_forward{};
    clock::duration communication_total_backward{};
    
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int solve_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    const int solve_root_rank =
        smallest_active_rank(tree->levels[0]);
    const LevelThreadPlan solve_thread_plan =
        configure_static_process_thread_plan(dynamic_threading);
    
    if (print_summary && rank == solve_header_rank) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Solve (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
        if (print_detail && dynamic_threading.enabled) {
            std::cout << "Solve thread plan: threads=" << solve_thread_plan.threads;
            if (solve_thread_plan.cpu_begin >= 0 && solve_thread_plan.cpu_end >= 0) {
                std::cout << ", cpus=[" << solve_thread_plan.cpu_begin
                          << ", " << solve_thread_plan.cpu_end << "]";
            }
            std::cout << std::endl;
        }
    }
    
    // ===== Initialize solve data structures =====
    
    
    for (int level = 0; level < num_levels; ++level) {
        auto& tree_level = tree->levels[level];
        
        if (!tree_level.is_process_active) {
            continue;
        }
        
        solve_data[level].resize(tree_level.num_boxes_local);
        int64_t global_idx = 0;
        for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
            auto& box = tree_level.local_boxes[box_idx];
            auto& solve_box = solve_data[level][box_idx];
            
            // solve_box.morton_index = box.morton_index;
            // solve_box.source_rank = rank;
            // solve_box.right_side.resize(box.num_points, DataType{0.0});
            // solve_box.left_side.resize(box.num_points, DataType{0.0});
            solve_box.initialize(box.morton_index, rank, box.num_points);
            
            solve_box.skeleton_indices = box.skeleton_indices;
            solve_box.redundant_indices = box.redundant_indices;

            
            
            // To Do: add righthand side with matrix
            // Set RHS from global vector (only at leaf level)
            if (level == leaf_level) {
                for (int64_t i = 0; i < box.num_points; ++i) {
                    solve_box.right_side[i] = rhs[global_idx];
                    global_idx++;
                }
                solve_box.left_side = solve_box.right_side;
            }
        }
    }

    
    if (print_detail && rank == solve_header_rank) {
        std::cout << "Initialized solve data structures" << std::endl;
    }
    
    // ===== Phase 1: Forward Sweep (V^{-1}) with level transitions =====
    
    if (print_detail && rank == solve_header_rank) {
        std::cout << "\n===== Phase 1: Forward Sweep (V^{-1}) =====" << std::endl;
    }
    
    auto forward_start = std::chrono::high_resolution_clock::now();
    PendingSolveUpdates<DataType> pending_solve;

    for (int level = leaf_level; level >= 1; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        // Error check: level 2→1 should not trigger reduction
        if (level == 2) {
            auto& parent_level = tree->levels[1];
            if (tree_level.num_active_processes != parent_level.num_active_processes) {
                throw std::runtime_error("Reduction between level 2 and level 1 is not allowed");
            }
        }
        if (print_detail && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << " on rank: " << rank << std::endl;
        }
       
        
        // // 1. Apply V^{-1} by color (skip level 1 - no far-field elimination)
        // if (level >= 2) {
        //     // 2^d color-wave forward sweep (local boxes only)
        //     const int num_colors = 1 << tree->dimension;

        //     // Bin ALL local boxes by color
        //     // Boundary-only colored bins
        //     std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

        //     // Interior boxes (local, in increasing Morton order by construction)
        //     std::vector<int64_t> interior_boxes;
        //     interior_boxes.reserve(tree_level.local_boxes.size());

        //     for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
        //         const int64_t morton_idx = tree_level.local_morton_start + local_idx;
        //         const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];

        //         if (box.on_boundary) {
        //             // 4-color in 2D / 8-color in 3D, based on morton % num_colors
        //             const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
        //             color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
        //         } else {
        //             // Interior: keep a simple list (already in morton order)
        //             interior_boxes.push_back(morton_idx);
        //         }
        //     }
        //     color_bins.push_back(std::move(interior_boxes));

        //     // for(auto &each_color : color_bins){
        //     //     each_color.clear();
        //     // }

        //     // for (size_t c = 0; c < num_colors; ++c) {
        //     //     int64_t start = (static_cast<int64_t>(c)     * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);
        //     //     int64_t end   = (static_cast<int64_t>(c + 1) * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);

        //     //     color_bins[c].reserve(static_cast<size_t>(end - start));
        //     //     for (int64_t x = start; x < end; ++x) {
        //     //         color_bins[c].push_back(x);
        //     //     }
                
        //     // }

        //     // pending updates generated by forward elimination to nonlocal neighbors
            

        //     // Process waves color-by-color
        //     for (int color_id = 0; color_id < color_bins.size(); ++color_id) {

        //         // Apply elimination for all local boxes of this color
        //         for (int64_t morton_idx : color_bins[(size_t)color_id]) {
        //             const int64_t local_idx = morton_idx - tree_level.local_morton_start;

        //             apply_forward_elimination(
        //                 tree_level,
        //                 solve_data[level][(size_t)local_idx],
        //                 solve_data[level],
        //                 fmm::MatrixProperty::SYMMETRIC,
        //                 pending_solve,
        //                 /*is_ghost=*/false
        //             );
        //         }

        //         // After each color wave, ship accumulated nonlocal neighbor updates and apply received ones.
        //         // (This matches the "wave" semantics you used for factor updates.)
        //         transport_and_apply_solve_updates_onehop(
        //             tree,
        //             tree_level,
        //             level,
        //             pending_solve,
        //             solve_data
        //         );
                
        //     }
            
           
        // }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;
            const int max_forward_threads = std::max(1, omp_get_max_threads());

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Both use morton & (num_colors-1), matching factorization order.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto&   box        = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id   = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }

            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout:
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {
                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (color_list.empty()) {
                    // Still need to transport even if no local boxes in this wave
                    communication_total_forward += transport_and_apply_solve_updates_onehop(
                        tree, tree_level, level, pending_solve, solve_data);
                    continue;
                }

                std::vector<PendingSolveUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_forward_threads));
                std::exception_ptr forward_wave_exception;
                std::mutex forward_wave_exception_mutex;
                std::atomic<bool> forward_wave_failed{false};

                #pragma omp parallel default(shared) if (color_list.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    auto& local_pending = thread_pending[static_cast<size_t>(tid)];

                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(color_list.size()); ++idx) {
                        if (forward_wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx = color_list[static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            apply_forward_elimination(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                local_pending,
                                /*is_ghost=*/false,
                                /*defer_local_updates=*/true
                            );
                        } catch (...) {
                            if (!forward_wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(forward_wave_exception_mutex);
                                forward_wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (forward_wave_exception) {
                    std::rethrow_exception(forward_wave_exception);
                }

                for (int t = 0; t < max_forward_threads; ++t) {
                    merge_pending_solve(pending_solve, thread_pending[static_cast<size_t>(t)]);
                    clear_pending_solve_updates_memory(thread_pending[static_cast<size_t>(t)]);
                }

                communication_total_forward += transport_and_apply_solve_updates_onehop(
                    tree, tree_level, level, pending_solve, solve_data);
            }
        }
        
        
        // Gather skeleton to parent level
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        gather_skeleton_to_parent(
            tree_level,
            parent_level,
            solve_data[level],
            solve_data[level - 1],
            tree->dimension
        );
        communication_total_forward += (clock::now() - get_data_start);

       
   
    }
    

    auto forward_end = std::chrono::high_resolution_clock::now();
    auto forward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(forward_end - forward_start - communication_total_forward);
    
    if (print_detail && rank == solve_root_rank) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }
    
    // ===== Phase 2: Diagonal Solve (D^{-1}) =====
    
    if (print_detail && rank == solve_root_rank) {
        std::cout << "\n===== Phase 2: Diagonal Solve (D^{-1}) =====" << std::endl;
    }
    
    auto diagonal_start = std::chrono::high_resolution_clock::now();
    
    int64_t num_diagonal_solves = 0;
    
    // Solve all X_RR blocks (levels N-1 down to 2)
    for (int level = leaf_level; level >= 2; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        // std::cout << "  Diagonal solves before from rank: " << rank << std::endl;
        if (!tree_level.is_process_active) {
            continue;
        }

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            int64_t local_diagonal_solves = 0;

            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) {
                    continue;
                }

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];
                    auto& solve_box = solve_data[level][static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) {
                        continue;
                    }

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    std::vector<DataType> b_R(static_cast<size_t>(r));
                    for (int64_t i = 0; i < r; ++i) {
                        b_R[static_cast<size_t>(i)] =
                            solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]];
                    }

                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        char uplo = 'L';
                        int n = static_cast<int>(r), nrhs = 1;
                        int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;

                        if constexpr (std::is_same_v<DataType, double>) {
                            dpotrs_(&uplo, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    b_R.data(), &ldb, &info);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            zsychol_solve_(&uplo, &n, &nrhs,
                                           box.X_RR.data.data(), &lda,
                                           b_R.data(), &ldb, &info);
                        }

                        if (info != 0) {
                            throw std::runtime_error("Diagonal solve failed for X_RR");
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal solve missing LU pivots for X_RR");
                        }

                        char trans = 'N';
                        int n = static_cast<int>(r), nrhs = 1;
                        int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;

                        if constexpr (std::is_same_v<DataType, double>) {
                            dgetrs_(&trans, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            zgetrs_(&trans, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                        }

                        if (info != 0) {
                            throw std::runtime_error("LU diagonal solve failed for X_RR");
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal solve missing Bunch-Kaufman pivots for X_RR");
                        }

                        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char uplo = 'L';
                            int n = static_cast<int>(r), nrhs = 1;
                            int lda = static_cast<int>(r), ldb = static_cast<int>(r), info = 0;
                            zsytrs_(&uplo, &n, &nrhs,
                                    box.X_RR.data.data(), &lda,
                                    box.X_RR_pivots.data(),
                                    b_R.data(), &ldb, &info);
                            if (info != 0) {
                                throw std::runtime_error("Bunch-Kaufman diagonal solve failed for X_RR");
                            }
                        } else {
                            throw std::runtime_error("BUNCH_KAUFMAN format only supported for complex<double>");
                        }
                    } else {
                        throw std::runtime_error("Unsupported X_RR format in diagonal solve");
                    }

                    for (int64_t i = 0; i < r; ++i) {
                        solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]] =
                            b_R[static_cast<size_t>(i)];
                    }

                    local_diagonal_solves++;
                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            num_diagonal_solves += local_diagonal_solves;
        }

        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }
    
    // Solve root X_RR (level 0)
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        apply_diagonal_solve(root_level, solve_data[0][0], false);
        num_diagonal_solves++;
    }
    
    auto diagonal_end = std::chrono::high_resolution_clock::now();
    auto diagonal_duration = std::chrono::duration_cast<std::chrono::milliseconds>(diagonal_end - diagonal_start);
    
    if (print_detail && rank == solve_root_rank) {
        std::cout << "  Diagonal solves: " << num_diagonal_solves  << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal solve time: " << diagonal_duration.count() << " ms" << std::endl;
    }
    
    
    // ===== Phase 3: Backward Sweep (W^{-1}) with level transitions =====
    
    if (print_detail && rank == solve_root_rank) {
        std::cout << "\n===== Phase 3: Backward Sweep (W^{-1}) =====" << std::endl;
    }
    
    auto backward_start = std::chrono::high_resolution_clock::now();
    
    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        
        if (print_detail && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << " from rank: " << rank << std::endl;
        }

        // for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
        //     auto& box = tree_level.local_boxes[box_idx];
        //     auto& solve_box = solve_data[level][box_idx];
            
        //     for(int gg = 0; gg < box.skeleton_indices.size(); gg++)
        //     {
        //         assert(solve_box.skeleton_indices[gg] == box.skeleton_indices[gg]);
        //     }
        // }
        
        // 1. Scatter solution from parent
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        if (tree_level.is_process_active || parent_level.is_process_active) {
            scatter_solution_to_children(
                tree_level,
                parent_level,
                solve_data[level],
                solve_data[level - 1],
                tree->dimension
            );
        }
        communication_total_backward += (clock::now() - get_data_start);

        // for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
        //     auto& box = tree_level.local_boxes[box_idx];
        //     auto& solve_box = solve_data[level][box_idx];
            
        //     for(int gg = 0; gg < box.skeleton_indices.size(); gg++)
        //     {
        //         assert(solve_box.skeleton_indices[gg] == box.skeleton_indices[gg]);
        //     }
        // }
        
        if (print_detail && rank == level_print_rank) {
            std::cout << "    ← Scattered from level " << (level - 1) << std::endl;
        }
        
        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            // ----------------------------------------------------------------
            // Build boundary color bins [0, num_colors-1]
            // Build interior sub-wave bins [num_colors, 2*num_colors-1]
            // Matches factorization binning exactly - reversed on iteration.
            // ----------------------------------------------------------------
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto&   box        = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id   = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }

            // Append interior sub-waves after boundary waves, matching forward order
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }
            // color_bins layout (same as forward):
            //   [0,           num_colors-1]: boundary waves
            //   [num_colors, 2*num_colors-1]: interior sub-waves
            //
            // Backward sweep: iterate in reverse, interior sub-waves first,
            // then boundary waves. Within each wave, reverse Morton order.

            for (int counter = static_cast<int>(color_bins.size()) - 1; counter >= 0; --counter) {
                // Gather new updates before each wave (same position as forward transport)
                get_data_start = clock::now();
                gather_boxes_solve(tree, level, solve_data[level], true);
                communication_total_backward += (clock::now() - get_data_start);

                auto& bins = color_bins[static_cast<size_t>(counter)];

                if (bins.empty()) continue;

                std::exception_ptr backward_wave_exception;
                std::mutex backward_wave_exception_mutex;
                std::atomic<bool> backward_wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (backward_wave_failed.load(std::memory_order_relaxed)) {
                            continue;
                        }

                        try {
                            const int64_t morton_idx =
                                bins[bins.size() - 1 - static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            apply_backward_substitution(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                /*is_ghost=*/false
                            );
                        } catch (...) {
                            if (!backward_wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(backward_wave_exception_mutex);
                                backward_wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (backward_wave_exception) {
                    std::rethrow_exception(backward_wave_exception);
                }
            }
        }
        // if (level >= 2) {
            
        //    // 2^d color-wave backward sweep (reverse order vs forward):
        //     // colors: (2^d - 1) ... 0
        //     // within each color: reverse Morton order
            
        //     const int num_colors = 1 << tree->dimension;

        //     // Bin ALL local boxes by color
        //     // Boundary-only colored bins
        //     std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

        //     // Interior boxes (local, in increasing Morton order by construction)
        //     std::vector<int64_t> interior_boxes;
        //     interior_boxes.reserve(tree_level.local_boxes.size());

        //     for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
        //         const int64_t morton_idx = tree_level.local_morton_start + local_idx;
        //         const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];

        //         if (box.on_boundary) {
        //             // 4-color in 2D / 8-color in 3D, based on morton % num_colors
        //             const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
        //             color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
        //         } else {
        //             // Interior: keep a simple list (already in morton order)
        //             interior_boxes.push_back(morton_idx);
        //         }
        //     }
        //     color_bins.push_back(std::move(interior_boxes));

        //     // for(auto &each_color : color_bins){
        //     //     each_color.clear();
        //     // }

        //     // for (size_t c = 0; c < num_colors; ++c) {
        //     //     int64_t start = (static_cast<int64_t>(c)     * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);
        //     //     int64_t end   = (static_cast<int64_t>(c + 1) * tree_level.local_boxes.size()) / static_cast<int64_t>(num_colors);

        //     //     color_bins[c].reserve(static_cast<size_t>(end - start));
        //     //     for (int64_t x = start; x < end; ++x) {
        //     //         color_bins[c].push_back(x);
        //     //     }
        //     // }

        //     // Backward: reverse color order, and reverse within-color order
        //     for (int color_id = color_bins.size() - 1; color_id >= 0; --color_id) {
        //         // gather new updates after every wave
        //         gather_boxes_solve(tree, level, solve_data[level], true);
                
        //         auto& bins = color_bins[(size_t)color_id];
        //         for (auto it = bins.rbegin(); it != bins.rend(); ++it) {
        //             const int64_t morton_idx = *it;
        //             const int64_t local_idx  = morton_idx - tree_level.local_morton_start;

        //             apply_backward_substitution(
        //                 tree_level,
        //                 solve_data[level][(size_t)local_idx],
        //                 solve_data[level],
        //                 fmm::MatrixProperty::SYMMETRIC,
        //                 /*is_ghost=*/false
        //             );
        //         }

                
        //     }
            
           
        // }
    }
    
    auto backward_end = std::chrono::high_resolution_clock::now();
    auto backward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(backward_end - backward_start - communication_total_backward);
    
    if (print_detail && rank == solve_root_rank) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }
    
    
    if (print_summary && rank == solve_root_rank) {
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - forward_start);
        auto communication_duration_forward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_forward);
        auto communication_duration_backward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_backward);
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total solve time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  Forward:  " << forward_duration.count() << " ms" << std::endl;
        std::cout << "  Diagonal: " << diagonal_duration.count() << " ms" << std::endl;
        std::cout << "  Backward: " << backward_duration.count() << " ms" << std::endl;
        std::cout << "  Communication (forward):  " << communication_duration_forward.count() << " ms" << std::endl;
        std::cout << "  Communication (backward):  " << communication_duration_backward.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);
}

/**
 * @brief Hierarchical multiply: y = F * x, with F the approximate factorization.
 *
 * F = (prod_i V_i) * P_l * D * P_l^* * (prod_i W_i)
 *
 * Phase 1 (forward, leaf -> root): apply W operators. W reads from neighbors
 *   but does not write to them -> use the same gather-before-wave pattern as
 *   solve's backward sweep (W^{-1}).
 *
 * Phase 2 (diagonal): multiply X_RR blocks.
 *
 * Phase 3 (backward, root -> leaf): apply V operators. V writes to neighbors
 *   -> use the same pending-update / transport pattern as solve's forward
 *   sweep (V^{-1}). Waves are traversed in reverse order.
 */
template<typename CoordType, typename DataType>
void hierarchical_mul_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& input_vec,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    bool verbose = true) {

    int rank;
    MPI_Comm_rank(tree->comm, &rank);
    using clock = std::chrono::high_resolution_clock;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    auto get_data_start = clock::now();
    clock::duration communication_total_forward{};
    clock::duration communication_total_backward{};

    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    const int mul_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    const int mul_root_rank =
        smallest_active_rank(tree->levels[0]);
    const LevelThreadPlan mul_thread_plan =
        configure_static_process_thread_plan(dynamic_threading);

    if (verbose && rank == mul_header_rank) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Multiply (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
        if (dynamic_threading.enabled) {
            std::cout << "Multiply thread plan: threads=" << mul_thread_plan.threads;
            if (mul_thread_plan.cpu_begin >= 0 && mul_thread_plan.cpu_end >= 0) {
                std::cout << ", cpus=[" << mul_thread_plan.cpu_begin
                          << ", " << mul_thread_plan.cpu_end << "]";
            }
            std::cout << std::endl;
        }
    }

    // ===== Initialize solve data structures =====
    for (int level = 0; level < num_levels; ++level) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        solve_data[level].resize(tree_level.num_boxes_local);
        int64_t global_idx = 0;
        for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
            auto& box = tree_level.local_boxes[box_idx];
            auto& solve_box = solve_data[level][box_idx];

            solve_box.initialize(box.morton_index, rank, box.num_points);
            solve_box.skeleton_indices = box.skeleton_indices;
            solve_box.redundant_indices = box.redundant_indices;

            if (level == leaf_level) {
                for (int64_t i = 0; i < box.num_points; ++i) {
                    solve_box.right_side[i] = input_vec[global_idx];
                    global_idx++;
                }
                solve_box.left_side = solve_box.right_side;
            }
        }
    }

    if (verbose && rank == mul_header_rank) {
        std::cout << "Initialized multiply data structures" << std::endl;
    }

    // ===== Phase 1: Forward Sweep (W) =====
    if (verbose && rank == mul_header_rank) {
        std::cout << "\n===== Phase 1: Forward Sweep (W) =====" << std::endl;
    }

    auto forward_start = std::chrono::high_resolution_clock::now();

    for (int level = leaf_level; level >= 1; level--) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);

        if (level == 2) {
            auto& parent_level = tree->levels[1];
            if (tree_level.num_active_processes != parent_level.num_active_processes) {
                throw std::runtime_error("Reduction between level 2 and level 1 is not allowed");
            }
        }
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local
                      << " boxes on rank: " << rank << std::endl;
        }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            // Build boundary-then-interior sub-wave bins (same as solve forward).
            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int color_id = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }

            // Apply W in forward color order (reversed of solve backward which
            // applies W^{-1}). Within each wave, gather fresh neighbor data
            // BEFORE applying (W is a gather, not scatter).
            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {
                // Fresh neighbor data for this wave.
                get_data_start = clock::now();
                gather_boxes_solve(tree, level, solve_data[level], true);
                communication_total_forward += (clock::now() - get_data_start);

                auto& bins = color_bins[static_cast<size_t>(counter)];
                if (bins.empty()) continue;

                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (wave_failed.load(std::memory_order_relaxed)) continue;

                        try {
                            const int64_t morton_idx = bins[static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            fmm::apply_mul_forward_W(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                /*is_ghost=*/false);
                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }
            }
        }

        // Gather skeleton values up to parent level (same as solve forward).
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        gather_skeleton_to_parent(
            tree_level,
            parent_level,
            solve_data[level],
            solve_data[level - 1],
            tree->dimension);
        communication_total_forward += (clock::now() - get_data_start);
    }

    auto forward_end = std::chrono::high_resolution_clock::now();
    auto forward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(forward_end - forward_start - communication_total_forward);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }

    // ===== Phase 2: Diagonal Multiply (D) =====
    if (verbose && rank == mul_root_rank) {
        std::cout << "\n===== Phase 2: Diagonal Multiply (D) =====" << std::endl;
    }

    auto diagonal_start = std::chrono::high_resolution_clock::now();

    int64_t num_diagonal_muls = 0;

    for (int level = leaf_level; level >= 2; level--) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            int64_t local_diagonal_muls = 0;

            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) continue;

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];
                    auto& solve_box = solve_data[level][static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) continue;

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    std::vector<DataType> x_R(static_cast<size_t>(r));
                    for (int64_t i = 0; i < r; ++i) {
                        x_R[static_cast<size_t>(i)] =
                            solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]];
                    }

                    int n = static_cast<int>(r);
                    int incx = 1;

                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        char uplo = 'L', diag_N = 'N';
                        int lda = static_cast<int>(r);

                        if constexpr (std::is_same_v<DataType, double>) {
                            char trans_T = 'T';
                            dtrmv_(&uplo, &trans_T, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char trans_N = 'N';
                            dtrmv_(&uplo, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char trans_T = 'T';
                            ztrmv_(&uplo, &trans_T, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char trans_N = 'N';
                            ztrmv_(&uplo, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal multiply missing LU pivots");
                        }
                        int lda = static_cast<int>(r);

                        if constexpr (std::is_same_v<DataType, double>) {
                            char uplo_U = 'U', trans_N = 'N', diag_N = 'N';
                            dtrmv_(&uplo_U, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char uplo_L = 'L', diag_U = 'U';
                            dtrmv_(&uplo_L, &trans_N, &diag_U, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            // P * z (not P^T * z) — use INCX=-1.
                            int nrhs = 1, k1 = 1, k2 = n, inc_rev = -1;
                            dlaswp_(&nrhs, x_R.data(), &n,
                                    &k1, &k2, box.X_RR_pivots.data(), &inc_rev);
                        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            char uplo_U = 'U', trans_N = 'N', diag_N = 'N';
                            ztrmv_(&uplo_U, &trans_N, &diag_N, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            char uplo_L = 'L', diag_U = 'U';
                            ztrmv_(&uplo_L, &trans_N, &diag_U, &n,
                                   box.X_RR.data.data(), &lda,
                                   x_R.data(), &incx);
                            int nrhs = 1, k1 = 1, k2 = n, inc_rev = -1;
                            zlaswp_(&nrhs, x_R.data(), &n,
                                    &k1, &k2, box.X_RR_pivots.data(), &inc_rev);
                        }
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() < static_cast<size_t>(r)) {
                            throw std::runtime_error("Diagonal multiply missing Bunch-Kaufman pivots");
                        }
                        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                            fmm::bunch_kaufman_multiply(n, box.X_RR.data.data(), n,
                                                       box.X_RR_pivots.data(), x_R.data());
                        } else {
                            throw std::runtime_error("BUNCH_KAUFMAN format only supported for complex<double>");
                        }
                    } else {
                        throw std::runtime_error("Diagonal multiply: unsupported X_RR format");
                    }

                    for (int64_t i = 0; i < r; ++i) {
                        solve_box.left_side[box.redundant_indices[static_cast<size_t>(i)]] =
                            x_R[static_cast<size_t>(i)];
                    }

                    local_diagonal_muls++;
                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            num_diagonal_muls += local_diagonal_muls;
        }

        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }

    // Root X_RR multiply
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        fmm::apply_diagonal_multiply(root_level, solve_data[0][0], false);
        num_diagonal_muls++;
    }

    auto diagonal_end = std::chrono::high_resolution_clock::now();
    auto diagonal_duration = std::chrono::duration_cast<std::chrono::milliseconds>(diagonal_end - diagonal_start);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Diagonal multiplies: " << num_diagonal_muls << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal multiply time: " << diagonal_duration.count() << " ms" << std::endl;
    }

    // ===== Phase 3: Backward Sweep (V) =====
    if (verbose && rank == mul_root_rank) {
        std::cout << "\n===== Phase 3: Backward Sweep (V) =====" << std::endl;
    }

    auto backward_start = std::chrono::high_resolution_clock::now();

    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);

        if (verbose && rank == level_print_rank) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local
                      << " boxes from rank: " << rank << std::endl;
        }

        // Scatter solution from parent (same as solve backward).
        auto& parent_level = tree->levels[level - 1];
        get_data_start = clock::now();
        if (tree_level.is_process_active || parent_level.is_process_active) {
            scatter_solution_to_children(
                tree_level,
                parent_level,
                solve_data[level],
                solve_data[level - 1],
                tree->dimension);
        }
        communication_total_backward += (clock::now() - get_data_start);

        if (verbose && rank == level_print_rank) {
            std::cout << "    <- Scattered from level " << (level - 1) << std::endl;
        }

        if (level >= 2 && tree_level.is_process_active) {
            const int num_colors = 1 << tree->dimension;

            std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));
            std::vector<std::vector<int64_t>> interior_sub_bins(static_cast<size_t>(num_colors));

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(tree_level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx = tree_level.local_morton_start + local_idx;
                const auto& box = tree_level.local_boxes[static_cast<size_t>(local_idx)];
                const int color_id = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
            }
            for (int c = 0; c < num_colors; ++c) {
                color_bins.push_back(std::move(interior_sub_bins[c]));
            }

            // Apply V in REVERSE wave order (mirror of solve forward's V^{-1}).
            // V is a scatter: generate pending updates, then transport.
            const int max_backward_threads = std::max(1, omp_get_max_threads());
            PendingSolveUpdates<DataType> pending_mul;

            for (int counter = static_cast<int>(color_bins.size()) - 1; counter >= 0; --counter) {
                auto& bins = color_bins[static_cast<size_t>(counter)];

                if (bins.empty()) {
                    // Still participate in transport so neighbors don't deadlock.
                    communication_total_backward += transport_and_apply_solve_updates_onehop(
                        tree, tree_level, level, pending_mul, solve_data);
                    continue;
                }

                std::vector<PendingSolveUpdates<DataType>> thread_pending(
                    static_cast<size_t>(max_backward_threads));
                std::exception_ptr wave_exception;
                std::mutex wave_exception_mutex;
                std::atomic<bool> wave_failed{false};

                #pragma omp parallel default(shared) if (bins.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    auto& local_pending = thread_pending[static_cast<size_t>(tid)];

                    #pragma omp for schedule(static)
                    for (int64_t idx = 0; idx < static_cast<int64_t>(bins.size()); ++idx) {
                        if (wave_failed.load(std::memory_order_relaxed)) continue;

                        try {
                            // Mirror solve backward's reverse-within-wave order.
                            const int64_t morton_idx =
                                bins[bins.size() - 1 - static_cast<size_t>(idx)];
                            const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                            fmm::apply_mul_backward_V_with_pending(
                                tree_level,
                                solve_data[level][static_cast<size_t>(local_idx)],
                                solve_data[level],
                                fmm::MatrixProperty::SYMMETRIC,
                                local_pending,
                                /*is_ghost=*/false);
                        } catch (...) {
                            if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                                std::lock_guard<std::mutex> lock(wave_exception_mutex);
                                wave_exception = std::current_exception();
                            }
                        }
                    }
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }

                for (int t = 0; t < max_backward_threads; ++t) {
                    merge_pending_solve(pending_mul, thread_pending[static_cast<size_t>(t)]);
                    clear_pending_solve_updates_memory(thread_pending[static_cast<size_t>(t)]);
                }

                communication_total_backward += transport_and_apply_solve_updates_onehop(
                    tree, tree_level, level, pending_mul, solve_data);
            }
        }
    }

    auto backward_end = std::chrono::high_resolution_clock::now();
    auto backward_duration = std::chrono::duration_cast<std::chrono::milliseconds>(backward_end - backward_start - communication_total_backward);

    if (verbose && rank == mul_root_rank) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }

    if (verbose && rank == mul_root_rank) {
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - forward_start);
        auto communication_duration_forward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_forward);
        auto communication_duration_backward = std::chrono::duration_cast<std::chrono::milliseconds>(communication_total_backward);
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total multiply time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  Forward:  " << forward_duration.count() << " ms" << std::endl;
        std::cout << "  Diagonal: " << diagonal_duration.count() << " ms" << std::endl;
        std::cout << "  Backward: " << backward_duration.count() << " ms" << std::endl;
        std::cout << "  Communication (forward):  " << communication_duration_forward.count() << " ms" << std::endl;
        std::cout << "  Communication (backward):  " << communication_duration_backward.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);
}


/**
 * @brief Verify solution using direct BLAS matrix-vector product
 * 
 * Builds full dense matrix A and computes A*x using BLAS.
 * WARNING: O(N²) memory and O(N²) time - only for small problems!
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @tparam KernelType Kernel evaluator type
 * @param kernel Kernel evaluator
 * @param rhs Original right-hand side vector
 * @param solution Computed solution vector
 * @param N Number of points
 * @param verbose Print detailed output
 * @return Relative residual norm
 */

} // namespace butterfly
