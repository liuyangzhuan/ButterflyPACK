#pragma once

#include "butterfly_types.hpp"
#include "butterfly_solve.hpp"
#include "butterfly_verification.hpp"

namespace butterfly {
using namespace fmm;

struct FactorizationCommunicatorSet {
    std::vector<MPI_Comm> level;
    std::vector<MPI_Comm> transition;
};

template<typename CoordType, typename DataType>
FactorizationCommunicatorSet make_factorization_communicators(
    fmm::ParallelTree<CoordType, DataType>* tree) {
    FactorizationCommunicatorSet communicators;
    communicators.level.assign(
        static_cast<size_t>(tree->num_levels), MPI_COMM_NULL);
    communicators.transition.assign(
        static_cast<size_t>(tree->num_levels), MPI_COMM_NULL);

    for (int level = 0; level < tree->num_levels; ++level) {
        const int color = tree->levels[level].is_process_active
            ? 0
            : MPI_UNDEFINED;
        MPI_Comm_split(
            tree->comm,
            color,
            tree->mpi_rank,
            &communicators.level[static_cast<size_t>(level)]);
    }

    for (int child_level = 1;
         child_level < tree->num_levels;
         ++child_level) {
        const bool participates =
            tree->levels[child_level].is_process_active ||
            tree->levels[child_level - 1].is_process_active;
        MPI_Comm_split(
            tree->comm,
            participates ? 0 : MPI_UNDEFINED,
            tree->mpi_rank,
            &communicators.transition[static_cast<size_t>(child_level)]);
    }

    return communicators;
}

inline void destroy_factorization_communicators(
    FactorizationCommunicatorSet& communicators) {
    for (MPI_Comm& comm : communicators.transition) {
        if (comm != MPI_COMM_NULL) {
            MPI_Comm_free(&comm);
        }
    }
    for (MPI_Comm& comm : communicators.level) {
        if (comm != MPI_COMM_NULL) {
            MPI_Comm_free(&comm);
        }
    }
}

// To Do: need to find better place for this function
// thresh: block-determinant magnitudes at or below this are treated as singular
// (throws instead of feeding 0 into std::log and producing -inf). Default 0.0
// guards only the exactly-zero case; pass a small positive value to also catch
// near-singular pivots.
template<typename DataType>
inline void accumulate_logdet_bunch_kaufman(const DataType* A, int64_t r, int64_t lda,
                                            const std::vector<int>& pivots,
                                            double& logabs, double& arg,
                                            double thresh = 1e-14) {
    int64_t k = 0;
    while (k < r) {
        if (pivots[static_cast<size_t>(k)] > 0) {          // 1×1 pivot
            DataType d = A[k + k * lda];
            double ad = std::abs(d);
            if (ad <= thresh) {
                throw std::runtime_error(
                    "accumulate_logdet_bunch_kaufman: singular/near-singular 1x1 pivot (|d| = " +
                    std::to_string(ad) + ")");
            }
            logabs += std::log(ad);
            arg    += std::arg(d);
            k += 1;
        } else {                                            // 2×2 pivot (pivots[k]==pivots[k+1]<0)
            if (k + 1 >= r) {
                throw std::runtime_error(
                    "accumulate_logdet_bunch_kaufman: 2x2 pivot at last index k=" +
                    std::to_string(k) + " (r=" + std::to_string(r) + "), malformed pivots");
            }
            DataType a = A[k       + k       * lda];
            DataType c = A[(k + 1) + (k + 1) * lda];
            DataType b = A[(k + 1) + k       * lda];        // sub-diagonal
            DataType det2 = a * c - b * b;
            double ad2 = std::abs(det2);
            if (ad2 <= thresh) {
                throw std::runtime_error(
                    "accumulate_logdet_bunch_kaufman: singular/near-singular 2x2 pivot (|det2| = " +
                    std::to_string(ad2) + ")");
            }
            logabs += std::log(ad2);
            arg    += std::arg(det2);
            k += 2;
        }
    }
}

template<typename DataType>
inline void accumulate_logdet_lu(const DataType* A, int64_t r, int64_t lda,
                                 const std::vector<int>& pivots,
                                 double& logabs, double& arg,
                                 double thresh = 1e-14) {
    if (pivots.size() != static_cast<size_t>(r)) {
        throw std::runtime_error(
            "accumulate_logdet_lu: pivot count (" + std::to_string(pivots.size()) +
            ") does not match matrix size (" + std::to_string(r) + ")");
    }

    for (int64_t k = 0; k < r; ++k) {
        const DataType diagonal = A[k + k * lda];
        const double magnitude = std::abs(diagonal);
        if (magnitude <= thresh) {
            throw std::runtime_error(
                "accumulate_logdet_lu: singular/near-singular pivot (|u_kk| = " +
                std::to_string(magnitude) + ")");
        }

        const int pivot = pivots[static_cast<size_t>(k)];
        if (pivot < 1 || pivot > r) {
            throw std::runtime_error(
                "accumulate_logdet_lu: invalid LAPACK pivot " +
                std::to_string(pivot) + " at index " + std::to_string(k));
        }

        logabs += std::log(magnitude);
        arg += std::arg(diagonal);
        if (pivot != k + 1) {
            arg += std::acos(-1.0);
        }
    }
}


template<typename CoordType, typename DataType>
void hierarchical_logdet_parallel(fmm::ParallelTree<CoordType, DataType>* tree,
                                  double* logabsdet, DataType* phase) {

    int leaf_level = tree->num_levels - 1;

    double logabs_local = 0.0;   // Σ log|det|
    double arg_local    = 0.0;   // Σ arg(det)

    // iterate through all the levels
    for (int level = leaf_level; level >= 2; level--) {

        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            double t_logabs = 0.0, t_arg = 0.0;   // thread-local

            // iterate through all the boxes (in parallel)
            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) {
                    continue;
                }

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) {
                        continue;
                    }

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    int64_t lda = box.X_RR.rows;
                    if (r != lda) {
                        throw std::runtime_error(
                            "logdet: redundant_indices.size() (" + std::to_string(r) +
                            ") != X_RR.rows (" + std::to_string(lda) + ")");
                    }
                    const DataType* A = box.X_RR.data.data();

                    // compute log|det| and phase for each box
                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        throw std::runtime_error(
                            "logdet: CHOLESKY_L is not implemented yet");
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        accumulate_logdet_lu(
                            A, r, lda, box.X_RR_pivots, t_logabs, t_arg);
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() != static_cast<size_t>(r)) {
                            throw std::runtime_error(
                                "logdet: pivots is incorrect (X_RR_pivots.size() = " +
                                std::to_string(box.X_RR_pivots.size()) +
                                ", expected " + std::to_string(r) + ")");
                        }
                        accumulate_logdet_bunch_kaufman(A, r, lda, box.X_RR_pivots, t_logabs, t_arg);
                    } else {
                        throw std::runtime_error("Diagonal multiply: unsupported X_RR format");
                    }

                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            logabs_local += t_logabs;
            #pragma omp atomic
            arg_local    += t_arg;
        }
        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }

    // handle root level
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        auto& rb = root_level.local_boxes[0];
        if (rb.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
            accumulate_logdet_bunch_kaufman(rb.X_RR.data.data(), rb.X_RR.rows, rb.X_RR.rows,
                                            rb.X_RR_pivots, logabs_local, arg_local);
        } else if (rb.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
            accumulate_logdet_lu(rb.X_RR.data.data(), rb.X_RR.rows, rb.X_RR.rows,
                                 rb.X_RR_pivots, logabs_local, arg_local);
        } else {
            throw std::runtime_error(
                "logdet (root): unsupported X_RR factorization format");
        }
    }

    // accumulate logabsdet and phase
    double buf[2] = { logabs_local, arg_local };
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, tree->comm);
    *logabsdet = buf[0];
    if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        *phase = DataType(std::cos(buf[1]), std::sin(buf[1]));   // e^{iθ}, any angle
    } else {
        *phase = std::cos(buf[1]);   // ±1: for real dets, θ is a multiple of π, so cos θ = ±1, sin θ ≈ 0
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void factorize_CA_level(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int current_level,
    KernelType* kernel,
    double tolerance,
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
    int level_print_rank) {
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

    if (level.num_active_processes == 1) {
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

    for (size_t group_idx = 0; group_idx < group_names.size(); ++group_idx) {
        const std::string& group_name = group_names[group_idx];
        const auto& group = *group_lists[group_idx];
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

        for (const auto& wave : waves) {
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

            #pragma omp parallel default(shared) if (wave.size() > 1)
            {
                FactorizationThreadScratch<CoordType, DataType> scratch;
                #pragma omp for schedule(static)
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

                        gather_id_workspace(
                            box, level, kernel,
                            unit_proxy_points.data(), num_proxy,
                            proxy_radius, is_symmetric,
                            scratch.workspace, scratch.workspace_rows,
                            scratch.workspace_cols, 0, box->on_boundary,
                            /*use_CA_boundary_semantics=*/true);
                        compute_and_modify(
                            dimension, box, level, kernel, scratch,
                            tolerance, is_symmetric, is_hermitian,
                            static_cast<PendingFactorUpdates<DataType>*>(nullptr),
                            factorization_method, use_owner_deferred_xnn);
                    } catch (...) {
                        if (!wave_failed.exchange(true, std::memory_order_relaxed)) {
                            std::lock_guard<std::mutex> lock(wave_exception_mutex);
                            wave_exception = std::current_exception();
                        }
                    }
                }
            }

            if (wave_exception) std::rethrow_exception(wave_exception);

            for (int64_t morton : wave) {
                level.eliminated_boxes.insert(morton);
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

                #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                {
                    DeferredXnnOwnerScratch<DataType> owner_scratch;

                    #pragma omp for schedule(static)
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
        }
    }

    for (const auto& box : level.local_boxes) {
        total_skeleton += static_cast<int64_t>(box.skeleton_indices.size());
        total_redundant += static_cast<int64_t>(box.redundant_indices.size());
        local_max_skel = std::max<int64_t>(
            local_max_skel, box.skeleton_indices.size());
    }
}

/**
 * @brief Hierarchical factorization routine (parallel)
 * 
 * Performs bottom-up hierarchical factorization from leaf level to root
 * with MPI communication and process reduction.
 * 
 * Algorithm:
 * 1. For each level from leaf (N-1) down to 1:
 *    a) Gather ghost and assisting boxes (MPI communication)
 *    b) Eliminate boxes in colored order: blue → orange → purple (3D) → green → interior
 *    c) Gather assisting boxes post-elimination (MPI communication)
 *    d) Build parent level interactions
 *    e) Handle process reduction (send/receive parent boxes if needed)
 * 2. At level 1 → 0 transition:
 *    - No elimination at level 1 (only 4/8 boxes, no far-field)
 *    - Aggregate to level 0 (single box or few boxes)
 *    - Factorize level 0 schur complement for diagonal solve
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @tparam KernelType Kernel evaluator type
 * @param tree Hierarchical tree structure
 * @param kernel Kernel evaluator
 * @param tolerance ID tolerance
 * @param is_symmetric Matrix symmetry flag
 * @param is_hermitian Matrix Hermitian flag
 * @param factorization_method Dense diagonal-block factorization method
 * @param unit_proxy_points Unit sphere proxy points
 * @param num_proxy Number of proxy points
 * @param proxy_radius Proxy sphere radius multiplier
 * @param verbosity -1: quiet, 0: summary, 1+: detailed progress
 */
template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel(
    fmm::ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t* out_rankmax,
    size_t* memory_per_rank,
    int verbosity = 1) {

    // To Do: NEED TO FIX KERNEL!!!!!
    using clock = std::chrono::high_resolution_clock;
    
    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    const bool print_summary = verbosity >= 0;
    const bool print_detail = verbosity >= 1;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    FactorizationCommunicatorSet factorization_comms =
        make_factorization_communicators(tree);
    
    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    int num_children = (dimension == 2) ? 4 : 8;
    const int factorization_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    
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

    for (int current_level = leaf_level; current_level >= 1; current_level--) {
        auto level_start = std::chrono::high_resolution_clock::now();
        clock::duration level_data_exchange{};
        clock::duration level_reduction{};
        
        auto& level = tree->levels[current_level];
        auto& parent_level = tree->levels[current_level - 1];
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

        if (level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            gather_CA_factorization_data(
                tree, current_level, level_comm,
                is_symmetric && !is_hermitian);
            const auto gather_duration = clock::now() - segment_start;
            level_data_exchange += gather_duration;
            if (print_detail && rank == level_print_rank) {
                std::cout << "  CA ghost/assisting gather time: "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 gather_duration).count()
                          << " ms" << std::endl;
            }
        }

        // declare locks
        std::unordered_map<int64_t, omp_lock_t*> box_locks;

         // initialize locks
        if (level.is_process_active) {
            for (auto& box : level.local_boxes) {
                auto* lock = new omp_lock_t;
                omp_init_lock(lock);
                box_locks[box.morton_index] = lock;
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
                    is_symmetric, is_hermitian, factorization_method,
                    unit_proxy_points, num_proxy, proxy_radius,
                    total_skeleton, total_redundant, local_max_skel,
                    print_detail, level_print_rank);
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

            for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(level.local_boxes.size()); ++local_idx) {
                const int64_t morton_idx  = level.local_morton_start + local_idx;
                const auto&   box         = level.local_boxes[static_cast<size_t>(local_idx)];
                const int     color_id    = static_cast<int>(morton_idx & (num_colors - 1));

                if (box.on_boundary) {
                    color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                } else {
                    interior_sub_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
                }
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
            int boundary_count = 0;
            bool to_store = true;
            const bool store_interior_wave = true;

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {

                const int  color_id_mod    = counter % num_colors;
                const bool is_interior     = (counter >= interior_start_loc);

                // ----------------------------------------------------------------
                // Communication / transport step (single-threaded)
                // ----------------------------------------------------------------
                if (!is_interior || counter == interior_start_loc) {
                    const auto comm_duration_raw = transport_and_apply_factor_updates_symmetric_onehop(
                        tree, current_level, kernel, pending_updates, false);
                    level_data_exchange += comm_duration_raw;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                    auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                    if (print_detail && rank == level_print_rank) {
                        std::cout << "  Comm time: "
                                << comm_duration.count()
                                << " ms" << std::endl;
                    }
                }
                

              

                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (print_detail && rank == level_print_rank) {
                    std::cout << "  Processing "
                            << (is_interior ? "interior sub-wave " : "boundary color ")
                            << color_id_mod
                            << " (" << color_list.size() << " boxes)..." << std::endl;
                }

                if (color_list.empty()) continue;

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

                #pragma omp parallel default(shared) if (color_list.size() > 1)
                {
                    const int tid = omp_get_thread_num();
                    FactorizationThreadScratch<CoordType, DataType> scratch;

                    #pragma omp for schedule(static)
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

                            gather_id_workspace(
                                &box, level, kernel,
                                unit_proxy_points.data(), num_proxy,
                                proxy_radius, is_symmetric,
                                scratch.workspace, scratch.workspace_rows, scratch.workspace_cols,
                                0,
                                box.on_boundary
                            );

                            thread_boundary_counts[static_cast<size_t>(tid)] += box.on_boundary;

                            compute_and_modify(dimension,
                                &box, level, kernel,
                                scratch,
                                tolerance, is_symmetric, is_hermitian,
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
                }

                if (wave_exception) {
                    std::rethrow_exception(wave_exception);
                }

                for (int local_boundary_count : thread_boundary_counts) {
                    boundary_count += local_boundary_count;
                }

                if (use_owner_deferred_xnn) {
                    for (int64_t morton_idx : color_list) {
                        level.eliminated_boxes.insert(morton_idx);
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

                    #pragma omp parallel default(shared) if (wave_xnn_candidate_boxes.size() > 1)
                    {
                        const int tid = omp_get_thread_num();
                        DeferredXnnOwnerScratch<DataType> scratch;

                        #pragma omp for schedule(static)
                        for (int64_t idx = 0; idx < static_cast<int64_t>(wave_xnn_candidate_boxes.size()); ++idx) {
                            if (owner_failed.load(std::memory_order_relaxed)) {
                                continue;
                            }

                            try {
                                const int64_t candidate_morton =
                                    wave_xnn_candidate_boxes[static_cast<size_t>(idx)];

                                // Deferred store=true replay:
                                //   local-local   -> local owner replay + local mirror
                                //   local-remote  -> local replay + remote transport
                                //   remote-remote -> remote transport only
                                apply_owner_deferred_xnn_updates_for_candidate_box(
                                    candidate_morton,
                                    level,
                                    kernel,
                                    wave_box_set,
                                    scratch,
                                    wave_xnn_mirror_targets[static_cast<size_t>(idx)],
                                    &thread_pending[static_cast<size_t>(tid)]);
                            } catch (...) {
                                if (!owner_failed.exchange(true, std::memory_order_relaxed)) {
                                    std::lock_guard<std::mutex> lock(owner_exception_mutex);
                                    owner_exception = std::current_exception();
                                }
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
                } else {
                    slice_far_field_blocks(level, is_symmetric, is_hermitian);
                }

                for (int t = 0; t < max_threads; ++t) {
                    merge_pending(pending_updates, thread_pending[t]);
                    clear_pending_factor_updates_memory(thread_pending[t]);
                }

                if (!use_owner_deferred_xnn) {
                    for (size_t bi = 0; bi < color_list.size(); ++bi) {
                        const int64_t morton_idx = color_list[bi];
                        level.eliminated_boxes.insert(morton_idx);
                    }
                }

                // ----------------------------------------------------------------
                // Mark assisting boxes as eliminated after this wave.
                // Boundary assisting boxes are marked by their boundary color wave.
                // Interior assisting boxes are marked by their interior sub-wave.
                // ----------------------------------------------------------------
                for (const auto& kv : level.assisting_box_points_for_kernel_evaluation) {
                    const int  assisting_color    = static_cast<int>(kv.first & (num_colors - 1));
                    const bool is_boundary_assist = level.assisting_boxes[kv.second].on_boundary;

                    const bool mark =
                        (!is_interior && is_boundary_assist  && assisting_color == color_id_mod) ||
                        ( is_interior && !is_boundary_assist && assisting_color == color_id_mod);

                    if (mark) {
                        level.eliminated_boxes.insert(kv.first);
                    }
                }

            }

            // Keep the final flush outside the wave body so empty trailing bins
            // cannot strand updates produced by the last nonempty wave.
            const auto final_comm_duration = transport_and_apply_factor_updates_symmetric_onehop(
                tree, current_level, kernel, pending_updates, false);
            level_data_exchange += final_comm_duration;
            update_neighbor_slicing_for_level(level, is_symmetric);
            
            for (const auto& box : level.local_boxes) {
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
                local_max_skel = std::max<int64_t>(local_max_skel, box.skeleton_indices.size());
            }
            }

            const auto elim_end = std::chrono::high_resolution_clock::now();
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                elim_end - elim_start);
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }

        if (current_level > 1 && level.is_process_active) {
            double min_elim_ms = 0.0;
            double max_elim_ms = 0.0;
            reduce_active_duration_bounds_ms(
                level_comm,
                0,
                true,
                elim_duration,
                min_elim_ms,
                max_elim_ms);
            
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Elimination time: shortest=" << std::llround(min_elim_ms)
                          << " ms, longest=" << std::llround(max_elim_ms) << " ms" << std::endl;
                
                const double compression =
                    static_cast<double>(total_skeleton) / (total_skeleton + total_redundant);
                std::cout << "  Compression ratio: " << compression
                          << " (" << (compression * 100) << "%)" << std::endl;
                std::cout << "  Average skeleton size: "
                          << static_cast<double>(total_skeleton) / level.local_boxes.size()
                          << std::endl;
            }
        }
        
        // ===== Step 3: Gather assisting boxes post-elimination =====
        if (current_level > 1 && level.is_process_active && use_CA_level) {
            segment_start = clock::now();
            gather_CA_assisting_boxes_factorization(
                tree, current_level, level_comm, 412);
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
        
        

        
        // ===== Step 4: Build parent level interactions =====
        
        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1 && level.is_process_active) {
            if (print_detail && rank == level_print_rank) {
                std::cout << "  Setting all DOFs as skeleton (no redundant DOFs at level 1)" << std::endl;
            }
            
            for (auto& box : level.local_boxes) {
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
            parent_boxes = build_parent_level_interactions<CoordType, DataType, KernelType>(
                level,
                tree->levels[current_level - 1],
                dimension,
                is_symmetric,
                is_hermitian,
                kernel,
                tree->global_bounds
            );
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
                
                std::vector<BoxData<CoordType, DataType>> child_parent_boxes =
                    deserialize_boxes<CoordType, DataType>(recv_buffer);
                
                all_parent_boxes.insert(
                    all_parent_boxes.end(),
                    std::make_move_iterator(child_parent_boxes.begin()),
                    std::make_move_iterator(child_parent_boxes.end())
                );
            }

            parent_level.local_boxes = std::move(all_parent_boxes);
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

            if (print_detail) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size() 
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }
        
        // Clear modified interaction matrices to free memory
        if (level.is_process_active) {
            clear_modified_interaction_matrices(level, use_CA_level);
        }

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
        
    }


    if (out_rankmax) *out_rankmax = local_max_skel;

    
    // ===== Special handling for level 0 (root) =====
    
    const int root_print_rank = smallest_active_rank(tree->levels[0]);
    if (print_detail && rank == root_print_rank) {
        std::cout << "\n===== Level 0 (Root) =====" << std::endl;
    }
    
    auto& root_level = tree->levels[0];
    
    if (!root_level.is_process_active || root_level.local_boxes.empty()) {
        if (print_detail) {
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

    if (print_summary) {
        double logabsdet;
        DataType phase;
        hierarchical_logdet_parallel(tree, &logabsdet, &phase);
        if (rank == 0) {
            std::cout.precision(17);
            std::cout << "logdet: " << phase << " " << logabsdet << std::endl;
        }

        (void)butterfly::h2_quick_verification(tree, kernel);
    }
}

template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel_if_supported(
    fmm::ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    int64_t* out_rankmax = nullptr,
    size_t* memory_per_rank = nullptr,
    int verbosity = 1) {
    if constexpr (std::is_same_v<DataType, double> ||
                  std::is_same_v<DataType, std::complex<double>>) {
        hierarchical_factorization_parallel<CoordType, DataType, KernelType>(
            tree, 
            kernel, 
            tolerance, 
            is_symmetric, 
            is_hermitian, 
            factorization_method, 
            unit_proxy_points, 
            num_proxy, 
            proxy_radius, 
            out_rankmax,
            memory_per_rank,
            verbosity);   // instantiated ONLY for double types
    } else {
        throw std::runtime_error("H2/FMM only supports double / std::complex<double>");
    }
}

template<typename CoordType, typename DataType>
void butterfly_factorization_parallel(H2<CoordType,DataType>* solver, double* factorization_time, double* entryeval_time) {

  if (solver->build_state == H2BuildState::RS_FACTORIZED) {
    if (factorization_time) *factorization_time = 0.0;
    if (entryeval_time) *entryeval_time = 0.0;
    return;
  }
  if (solver->build_state == H2BuildState::H2_COMPRESSED) {
    throw std::runtime_error(
      "butterfly_factorization_parallel: solver already contains a compression-only H2 representation");
  }

  int rank = 0;
  MPI_Comm_rank(solver->comm, &rank);

  const auto factorization_method =
    fmm::FactorizationMethod::BUNCH_KAUFMAN;
  fmm::HierarchicalFactorization<CoordType, DataType, butterfly::H2Kernel<CoordType, DataType>> factorizer(
    solver->options.N,
    fmm::MatrixProperty::SYMMETRIC,
    &solver->kernel,
    solver->options.dimension,
    factorization_method,
    solver->options.num_proxy);

  const auto& unit_proxy = factorizer.get_unit_proxy_points();
  const int num_proxy = factorizer.get_num_proxy_points();

  const bool is_symmetric = true;
  const bool is_hermitian = false;
  
  solver->kernel.entryeval_time_per_thread.assign(omp_get_max_threads(), 0.0);
  auto& ev = solver->kernel.entryeval_time_per_thread;

//   auto total_start = std::chrono::high_resolution_clock::now();
  double t0 = MPI_Wtime();
  butterfly::hierarchical_factorization_parallel_if_supported(
    solver->tree.get(),
    &solver->kernel,
    solver->options.tolerance,
    is_symmetric,
    is_hermitian,
    factorization_method,
    unit_proxy,
    num_proxy,
    0.0, // proxy_radius = 0
    &solver->last_factor_rankmax,
    &solver->factorization_memory,
    solver->options.verbosity);
  double tf = MPI_Wtime() - t0;
  MPI_Allreduce(MPI_IN_PLACE, &tf, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
  *factorization_time = tf;
  
  double t_entry = 0.0;
  if (!ev.empty()) {
    t_entry = *std::max_element(ev.begin(), ev.end());
  }
  MPI_Allreduce(MPI_IN_PLACE, &t_entry, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
  *entryeval_time = t_entry;
  solver->factorized = true;
  solver->build_state = H2BuildState::RS_FACTORIZED;

  MPI_Allreduce(MPI_IN_PLACE, &solver->last_factor_rankmax, 1, MPI_INT64_T, MPI_MAX, solver->comm);

//   auto total_end = std::chrono::high_resolution_clock::now();
//   auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//     total_end - total_start);


//   if (rank == 0) {
//     std::cout << "\n========================================" << std::endl;
//     std::cout << "Total factorization time: " << total_duration.count() << " ms" << std::endl;
//     std::cout << "========================================\n" << std::endl;
//   }

}

/**
 * @brief Hierarchical solve routine (parallel MPI version)
 */

} // namespace butterfly
