#pragma once

#include "butterfly_types.hpp"
#include "butterfly_solve.hpp"

namespace butterfly {
using namespace fmm;

template<typename DataType, typename KernelType>
double verify_solution_direct(
    MPI_Comm comm,
    KernelType* kernel,
    const std::vector<DataType>& rhs,
    const std::vector<DataType>& solution,
    int64_t N,
    int dimension = 2,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    if (rank != 0) {
        return 0.0;  // Only verify on rank 0
    }
    
    if (verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Solution Verification (Direct BLAS)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Total DOFs: " << N << std::endl;
        
        double matrix_memory_mb = (static_cast<double>(N) * N * sizeof(DataType)) / (1024.0 * 1024.0);
        std::cout << "Matrix memory: " << matrix_memory_mb << " MB" << std::endl;
        
        if (N > 10000) {
            std::cout << "⚠ WARNING: N = " << N << " is large for direct method!" << std::endl;
        }
    }
    
    // ===== Step 1: Build full matrix A =====
    
    auto build_start = std::chrono::high_resolution_clock::now();
    
    std::vector<DataType> A(N * N);
    
    // Build A from the SAME entry evaluation the H2 solver compressed/factored (Zmn via
    // evaluate_by_index), stored column-major (A[i + j*N]) to match the BLAS matvec below.
    // No /N scaling and no special diagonal case: evaluate_by_index returns the exact matrix
    // entry for every (i, j), including the self term when i == j.
    #pragma omp parallel for collapse(2) if(N > 1000)
    for (int64_t i = 0; i < N; ++i) {
        for (int64_t j = 0; j < N; ++j) {
            kernel->evaluate_by_index(i, j, A.data(), N);   // writes A[i + j*N] = K(i, j)
        }
    }
    
    auto build_end = std::chrono::high_resolution_clock::now();
    auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start);
    
    if (verbose) {
        std::cout << "Matrix assembly time: " << build_duration.count() << " ms" << std::endl;
    }
    
    // ===== Step 2: Compute A*x using BLAS =====
    
    auto matvec_start = std::chrono::high_resolution_clock::now();
    
    std::vector<DataType> Ax(N, DataType{0.0});
    
    if constexpr (std::is_same_v<DataType, double>) {
        // DGEMV: y = alpha*A*x + beta*y
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        double alpha = 1.0;
        double beta = 0.0;
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        dgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, float>) {
        // SGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        float alpha = 1.0f;
        float beta = 0.0f;
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        sgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        // ZGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        zgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        // CGEMV
        char trans = 'N';
        int m = static_cast<int>(N);
        int n = static_cast<int>(N);
        std::complex<float> alpha(1.0f, 0.0f);
        std::complex<float> beta(0.0f, 0.0f);
        int lda = static_cast<int>(N);
        int incx = 1;
        int incy = 1;

        cgemv_(&trans, &m, &n, &alpha, A.data(), &lda,
            solution.data(), &incx, &beta, Ax.data(), &incy);

    } else {
        throw std::runtime_error("Unsupported DataType for BLAS");
    }
    
    auto matvec_end = std::chrono::high_resolution_clock::now();
    auto matvec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(matvec_end - matvec_start);
    
    if (verbose) {
        std::cout << "BLAS matvec time: " << matvec_duration.count() << " ms" << std::endl;
    }
    
    // ===== Step 3: Compute residual =====
    
    // Deduce the underlying real type
    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;

    std::vector<DataType> residual(N);
    for (int64_t i = 0; i < N; ++i) {
        residual[i] = Ax[i] - rhs[i];
    }

    // Compute norms using BLAS
    RealType residual_norm = 0.0;
    RealType rhs_norm = 0.0;

    if constexpr (std::is_same_v<DataType, double>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = dnrm2_(&n, residual.data(), &inc);
        rhs_norm = dnrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, float>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = snrm2_(&n, residual.data(), &inc);
        rhs_norm = snrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = dznrm2_(&n, residual.data(), &inc);
        rhs_norm = dznrm2_(&n, rhs.data(), &inc);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        int n = static_cast<int>(N);
        int inc = 1;
        residual_norm = scnrm2_(&n, residual.data(), &inc);
        rhs_norm = scnrm2_(&n, rhs.data(), &inc);
    }

    RealType relative_error = residual_norm / rhs_norm;

    if (verbose) {
        std::cout << "\nBackward Error Analysis:" << std::endl;
        std::cout << "  ||Ax - b||₂ = " << std::scientific << std::setprecision(6)
                  << residual_norm << std::endl;
        std::cout << "  ||b||₂      = " << rhs_norm << std::endl;
        std::cout << "  Relative residual = " << relative_error << std::endl;

        // Additional statistics
        RealType max_residual = 0.0;
        int64_t max_idx = 0;
        for (int64_t i = 0; i < N; ++i) {
            RealType abs_res = std::abs(residual[i]);
            if (abs_res > max_residual) {
                max_residual = abs_res;
                max_idx = i;
            }
        }

        std::cout << "  Max residual    = " << max_residual
                  << " (at index " << max_idx << ")" << std::endl;

        // Solution quality assessment
        if (relative_error < static_cast<RealType>(1e-10)) {
            std::cout << "\n  ✓ EXCELLENT: Relative error < 1e-10" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-6)) {
            std::cout << "\n  ✓ VERY GOOD: Relative error < 1e-6" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-3)) {
            std::cout << "\n  ✓ GOOD: Relative error < 1e-3" << std::endl;
        } else if (relative_error < static_cast<RealType>(1e-1)) {
            std::cout << "\n  ⚠ WARNING: Relative error > 1e-3" << std::endl;
        } else {
            std::cout << "\n  ✗ ERROR: Relative error > 1e-1 (solution likely incorrect)" << std::endl;
        }

        std::cout << "========================================\n" << std::endl;
    }
    
    return relative_error;
}

/**
 * @brief Gather solution and RHS from all processes to process 0
 * 
 * @tparam CoordType Coordinate type (float or double)
 * @tparam DataType Data type for matrix entries
 * 
 * @param tree Parallel tree structure
 * @param solve_data Solve data at all levels (only leaf level used)
 * @param solution Output solution vector (only populated on rank 0)
 * @param aggregated_rhs Output RHS vector (only populated on rank 0)
 * 
 * Gathers point_indices, left_side, and right_side from all leaf boxes 
 * across all processes and assembles into global vectors on rank 0.
 */
template<typename CoordType, typename DataType>
void gather_solution_to_root(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data,
    std::vector<DataType>& solution,
    std::vector<DataType>& aggregated_rhs) {
    
    int rank, size;
    MPI_Comm_rank(tree->comm, &rank);
    MPI_Comm_size(tree->comm, &size);
    
    int leaf_level = tree->num_levels - 1;
    auto& leaf_level_ref = tree->levels[leaf_level];

    const int64_t local_solve_levels = static_cast<int64_t>(solve_data.size());
    if (leaf_level < 0 || leaf_level >= static_cast<int>(solve_data.size())) {
        std::fprintf(stderr,
                     "gather_solution_to_root: rank=%d invalid leaf level %d for solve_data size %lld\n",
                     rank,
                     leaf_level,
                     static_cast<long long>(local_solve_levels));
        std::fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int64_t expected_leaf_boxes = leaf_level_ref.num_boxes_local;
    const int64_t actual_leaf_boxes = static_cast<int64_t>(solve_data[leaf_level].size());
    if (actual_leaf_boxes != expected_leaf_boxes) {
        std::fprintf(stderr,
                     "gather_solution_to_root: rank=%d solve_data leaf mismatch: expected %lld boxes, got %lld entries\n",
                     rank,
                     static_cast<long long>(expected_leaf_boxes),
                     static_cast<long long>(actual_leaf_boxes));
        std::fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // ===== Step 1: Collect local indices, solution values, and RHS values =====
    
    std::vector<int64_t> local_indices;
    std::vector<DataType> local_solution;
    std::vector<DataType> local_rhs;
    
    for (int64_t box_idx = 0; box_idx < leaf_level_ref.num_boxes_local; ++box_idx) {
        auto& box = leaf_level_ref.local_boxes[box_idx];
        auto& solve_box = solve_data[leaf_level][box_idx];
        
        for (int64_t i = 0; i < box.num_points; ++i) {
            local_indices.push_back(box.point_indices[i]);
            local_solution.push_back(solve_box.left_side[i]);
            local_rhs.push_back(solve_box.right_side[i]);
        }
    }
    
    // ===== Step 2: Gather counts on root =====

    const size_t local_count = local_indices.size();
    if (local_count > static_cast<size_t>(INT_MAX)) {
        throw std::runtime_error(
            "gather_solution_to_root: local verification gather count exceeds INT_MAX on rank " +
            std::to_string(rank) + ": " + std::to_string(local_count));
    }

    const int local_count_int = static_cast<int>(local_count);
    std::vector<int> recv_counts(rank == 0 ? size : 0, 0);
    int err = MPI_Gather(&local_count_int,
                         1,
                         MPI_INT,
                         rank == 0 ? recv_counts.data() : nullptr,
                         1,
                         MPI_INT,
                         0,
                         tree->comm);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gather(counts) failed with error " +
            std::string(errbuf, errlen));
    }
    
    // ===== Step 3: Calculate displacements on rank 0 =====
    
    int total_count = 0;
    std::vector<int> recv_displs(rank == 0 ? size : 0, 0);
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            if (i > 0) {
                recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
            }
            if (recv_counts[i] < 0 || total_count > INT_MAX - recv_counts[i]) {
                throw std::runtime_error(
                    "gather_solution_to_root: total verification gather count exceeds INT_MAX");
            }
            total_count += recv_counts[i];
        }
    }
    
    // ===== Step 4: Gather indices, solution, and RHS on root =====
    
    const MPI_Datatype mpi_datatype = mpi_datatype_for<DataType>();

    std::vector<int64_t> all_indices(rank == 0 ? static_cast<size_t>(total_count) : 0);
    std::vector<DataType> all_solution(rank == 0 ? static_cast<size_t>(total_count) : 0);
    std::vector<DataType> all_rhs(rank == 0 ? static_cast<size_t>(total_count) : 0);

    err = MPI_Gatherv(local_indices.empty() ? nullptr : local_indices.data(),
                      local_count_int,
                      MPI_INT64_T,
                      rank == 0 ? all_indices.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      MPI_INT64_T,
                      0,
                      tree->comm);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(indices) failed with error " +
            std::string(errbuf, errlen));
    }

    err = MPI_Gatherv(local_solution.empty() ? nullptr : local_solution.data(),
                      local_count_int,
                      mpi_datatype,
                      rank == 0 ? all_solution.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      mpi_datatype,
                      0,
                      tree->comm);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(solution) failed with error " +
            std::string(errbuf, errlen));
    }

    err = MPI_Gatherv(local_rhs.empty() ? nullptr : local_rhs.data(),
                      local_count_int,
                      mpi_datatype,
                      rank == 0 ? all_rhs.data() : nullptr,
                      rank == 0 ? recv_counts.data() : nullptr,
                      rank == 0 ? recv_displs.data() : nullptr,
                      mpi_datatype,
                      0,
                      tree->comm);
    if (err != MPI_SUCCESS) {
        char errbuf[MPI_MAX_ERROR_STRING];
        int errlen = 0;
        MPI_Error_string(err, errbuf, &errlen);
        throw std::runtime_error(
            "gather_solution_to_root: MPI_Gatherv(rhs) failed with error " +
            std::string(errbuf, errlen));
    }
    
    // ===== Step 5: Assemble solution and RHS on rank 0 =====
    
    if (rank == 0) {
        solution.resize(tree->num_points, DataType{0.0});
        aggregated_rhs.resize(tree->num_points, DataType{0.0});
        
        for (int i = 0; i < total_count; ++i) {
            int64_t global_idx = all_indices[i];
            
            if (global_idx < 0 || global_idx >= tree->num_points) {
                throw std::runtime_error(
                    "gather_solution_to_root: Invalid global index " + 
                    std::to_string(global_idx));
            }
            
            solution[global_idx] = all_solution[i];
            aggregated_rhs[global_idx] = all_rhs[i];
        }
    }
}

//provide booleans for whether to do certain checks
template<typename CoordType, typename DataType, typename KernelType>
inline int h2_direct_verification(ParallelTree<CoordType, DataType>* tree, 
    const std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data,
    KernelType* kernel, 
    const ProgramOptions& options) {

    int rank; 
    MPI_Comm_rank(tree->comm, &rank);

    // Invariant: options.N (the Npo handed to construct_init) must equal the built tree's
    // point count. Both are derived from the same Npo, so a mismatch means the options struct
    // and the tree are inconsistent -- fail loudly rather than size the verification wrong.
    // (Global, rank-consistent scalars, so every rank takes the same branch: abort is safe.)
    if (options.N != tree->num_points) {
      if (rank == 0) {
        std::cerr << "h2_direct_verification: options.N (" << options.N
                  << ") != tree->num_points (" << tree->num_points
                  << "); inconsistent solver state." << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int64_t N = tree->num_points;

    // verify solution x from Ax=b is reasonable, this is the ground truth comparison for smaller matrices
    if (N <= 8000) {

      // Verification is purely index-based (evaluate_by_index): it needs no coordinates and
      // imposes no grid-shape requirement, so any N <= 4096 is verifiable.
      std::vector<DataType> solution;
      std::vector<DataType> aggregated_rhs;

      butterfly::gather_solution_to_root(tree, solve_data, solution, aggregated_rhs);
      
      const double direct_error = butterfly::verify_solution_direct(
        tree->comm,
        kernel,
        aggregated_rhs,
        solution,
        N,
        options.dimension,
        true);
      if (rank == 0) {
        std::cout << "verifying with direct matrix vector multiply since N <= 4096" << std::endl;
        std::cout << "Direct error: " << direct_error << std::endl;
      }
    } else {
        if (rank ==0) {
            std::cout << "no direct verification. N > 4096." <<std::endl;
        }
        
    }

    return 0;
}

inline uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

template<typename RealType>
RealType centered_uniform_from_hash(uint64_t key) {
    constexpr long double scale =
        1.0L / static_cast<long double>(std::numeric_limits<uint64_t>::max());
    const long double u =
        static_cast<long double>(splitmix64(key)) * scale;
    return static_cast<RealType>(u - 0.5L);
}

template<typename DataType>
DataType make_verification_entry(int64_t global_index, uint64_t seed) {
    const uint64_t base =
        splitmix64(seed ^ static_cast<uint64_t>(global_index));
    if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        return DataType(
            centered_uniform_from_hash<double>(base),
            centered_uniform_from_hash<double>(base ^ 0x6a09e667f3bcc909ULL));
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        return DataType(
            centered_uniform_from_hash<float>(base),
            centered_uniform_from_hash<float>(base ^ 0x6a09e667f3bcc909ULL));
    } else {
        return centered_uniform_from_hash<DataType>(base);
    }
}

// Step 1: pick `num_src` distinct random columns in [0, N) and give each a random weight.
// Fully deterministic in (N, num_src, seed): every rank calls this and gets the identical
// result, so no MPI_Bcast is needed to agree on x.
template<typename DataType>
SparseTestVector<DataType> make_sparse_test_vector(int64_t N, int64_t Npt_src, uint64_t seed) {
    SparseTestVector<DataType> stv;

    stv.idx.reserve(static_cast<size_t>(Npt_src));
    stv.weight.reserve(static_cast<size_t>(Npt_src));

    std::unordered_set<int64_t> seen;
    seen.reserve(static_cast<size_t>(Npt_src) * 2);

    // draw distinct candidate columns via the hash; dedup until we have Npt_src of them
    for (uint64_t t = 0; static_cast<int64_t>(stv.idx.size()) < Npt_src; ++t) {
        const int64_t cand = static_cast<int64_t>(
            splitmix64(seed ^ (0xD1B54A32D192ED03ULL * (t + 1))) % static_cast<uint64_t>(N));
        if (seen.insert(cand).second) {
            stv.idx.push_back(cand);
            stv.weight.push_back(make_verification_entry<DataType>(cand, seed));
        }
    }
    return stv;
}


template<typename CoordType, typename DataType, typename KernelType>
double h2_verification_smvp(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,               // H2Kernel with evaluate_block_by_index
    int num_src = 20,                 // sparse support size (Fortran uses 20)
    unsigned seed = 12345,            // fixed → deterministic regression check
    bool verbose = true) {

    
    int rank, size;
    MPI_Comm_rank(tree->comm, &rank);
    MPI_Comm_size(tree->comm, &size);

    int64_t N = tree->num_points;
    const int leaf_level = tree->num_levels - 1;
    auto& leaf = tree->levels[leaf_level];


    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<std::is_same_v<DataType, std::complex<float>>, float, DataType>>;

    auto mag2 = [](const DataType& z) -> RealType {
        RealType a = std::abs(z);   // works for real and complex
        return a * a;
    };

    // step 0: setting up sparse vector and distribute into local segments

    // generate sparse test vector
    const int64_t Npt_src = std::min<int64_t>(num_src, N);
    const SparseTestVector<DataType> stv = make_sparse_test_vector<DataType>(N, Npt_src, seed);

    // define input_vec, local segments for each MPI rank 
    std::unordered_set<int64_t> src_set(stv.idx.begin(), stv.idx.end());

    std::vector<DataType> input_vec;
    std::vector<int64_t>  local_rows;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& box = leaf.local_boxes[b];
        for (int64_t i = 0; i < box.num_points; ++i) {
            const int64_t g = box.point_indices[i];
            local_rows.push_back(g);
            // nonzero only on the support; value is the SAME deterministic weight used above
            input_vec.push_back(src_set.count(g) ? make_verification_entry<DataType>(g, seed)
                                                 : DataType{0.0});
        }
    }
    const int64_t Nloc = static_cast<int64_t>(input_vec.size());

    // step 1: calculate A*x for each x --> b
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>
        mul_data(tree->num_levels);

    const double t0 = MPI_Wtime();
    butterfly::hierarchical_mul_parallel(tree, input_vec, mul_data, false);
    double t_mul = MPI_Wtime() - t0;

    // step 2: calculate A_true*x for each x
    std::vector<DataType> block(static_cast<size_t>(Nloc) * static_cast<size_t>(Npt_src));
    if (Nloc > 0) {
        kernel->evaluate_block_by_index(local_rows.data(), Nloc,
                                        stv.idx.data(),    Npt_src,
                                        block.data(),      Nloc);
    }

    RealType approx_Ax = 0, exact_Ax = 0, diff = 0;
    int64_t pos = 0;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& solve_box = mul_data[leaf_level][b];
        for (size_t i = 0; i < solve_box.left_side.size(); ++i, ++pos) {
            DataType y_ref{0.0};
            for (int64_t k = 0; k < Npt_src; ++k)
                y_ref += stv.weight[k] * block[pos + k * Nloc];   // exact weighted column sum

            const DataType y_h2 = solve_box.left_side[i];     // compressed A x
            approx_Ax += mag2(y_h2);
            exact_Ax += mag2(y_ref);
            diff += mag2(y_h2 - y_ref);
        }
    }

    // step 3: find error \|A_true*x - b\|

    RealType localv[3]  = { approx_Ax, exact_Ax, diff };
    RealType globalv[3] = { 0, 0, 0 };
    const MPI_Datatype real_mpi =
        std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(localv, globalv, 3, real_mpi, MPI_SUM, tree->comm);

    MPI_Allreduce(MPI_IN_PLACE, &t_mul, 1, MPI_DOUBLE, MPI_MAX, tree->comm);

    const RealType nrm_h2  = std::sqrt(globalv[0]);
    const RealType nrm_ref = std::sqrt(globalv[1]);
    const RealType acc     = (globalv[0] > RealType(0))
                             ? std::sqrt(globalv[2] / globalv[1]) : RealType(0);

    if (rank == 0 && verbose) {
        std::cout << "H2_CheckError(mvp): fnorm: "
                  << std::scientific << std::setprecision(7)
                  << nrm_h2 << "  " << nrm_ref
                  << "  acc: " << acc
                  << "  time: " << t_mul << std::endl;
    }
    return static_cast<double>(acc);
}

} // namespace butterfly
