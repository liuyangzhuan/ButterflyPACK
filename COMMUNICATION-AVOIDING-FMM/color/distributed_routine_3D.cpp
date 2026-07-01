#include "factorization.hpp"
#include "solver.hpp"
#include "tree_impl.hpp"
#include "kernel.hpp"
#include "id_decomposition.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <random>
#include <complex>
#include <sched.h>

namespace fmm {


/**
* @brief Print vector values for debugging
*/
template<typename DataType>
void print_vector(
    const std::vector<DataType>& vec,
    const std::string& name,
    int64_t max_print = 20) {
    
    std::cout << "\n" << name << " (size " << vec.size() << "):" << std::endl;
    
    // Print first few entries
    int64_t n_print = std::min(max_print, static_cast<int64_t>(vec.size()));
    
    std::cout << "  First " << n_print << " entries:" << std::endl;
    for (int64_t i = 0; i < n_print; ++i) {
        std::cout << "    [" << std::setw(4) << i << "] = " 
                  << std::scientific << std::setprecision(6) << vec[i] << std::endl;
    }
    
    // // Print last few entries
    // if (vec.size() > max_print * 2) {
    //     std::cout << "  ..." << std::endl;
    //     std::cout << "  Last " << n_print << " entries:" << std::endl;
    //     for (int64_t i = vec.size() - n_print; i < vec.size(); ++i) {
    //         std::cout << "    [" << std::setw(4) << i << "] = " 
    //                   << std::scientific << std::setprecision(6) << vec[i] << std::endl;
    //     }
    // }
    
    // Compute statistics
    DataType min_val = *std::min_element(vec.begin(), vec.end());
    DataType max_val = *std::max_element(vec.begin(), vec.end());
    DataType sum = std::accumulate(vec.begin(), vec.end(), DataType{0.0});
    DataType mean = sum / vec.size();
    
    DataType variance = 0.0;
    for (const auto& val : vec) {
        variance += (val - mean) * (val - mean);
    }
    variance /= vec.size();
    DataType std_dev = std::sqrt(variance);
    
    std::cout << "\n  Statistics:" << std::endl;
    std::cout << "    Min:    " << std::scientific << min_val << std::endl;
    std::cout << "    Max:    " << max_val << std::endl;
    std::cout << "    Mean:   " << mean << std::endl;
    std::cout << "    Std:    " << std_dev << std::endl;
    std::cout << "    L2norm: " << std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0)) << std::endl;
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
 * @param grid_points Grid point coordinates
 * @param N Number of points
 * @param verbose Print detailed output
 * @return Relative residual norm
 */
template<typename CoordType, typename DataType, typename KernelType>
DataType verify_solution_direct(
    KernelType* kernel,
    const std::vector<DataType>& rhs,
    const std::vector<DataType>& solution,
    const std::vector<CoordType>& grid_points,
    int64_t N,
    int dimension = 2,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
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
    
    #pragma omp parallel for collapse(2) if(N > 1000)
    for (int64_t i = 0; i < N; ++i) {
        for (int64_t j = 0; j < N; ++j) {
            if (i == j) {
                // Diagonal entry
                A[i * N + j] = kernel->evaluate_diagonal(N);
            } else {
                // Off-diagonal entry
                CoordType xi[3], xj[3];
                for (int d = 0; d < dimension; ++d) {
                    xi[d] = grid_points[i * dimension + d];
                    xj[d] = grid_points[j * dimension + d];
                }
                
                A[i * N + j] = kernel->evaluate(xi, xj) / static_cast<DataType>(N);
            }
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
        std::cout << "\nResidual Analysis:" << std::endl;
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
 * @brief Verify solution using FFT-based matrix-vector product
 * 
 * Computes relative residual ||Ax - b|| / ||b|| where A*x is computed via FFT.
 * 
 * @tparam CoordType Coordinate data type
 * @tparam DataType Matrix data type
 * @param tree Hierarchical tree with solution
 * @param rhs Original right-hand side vector
 * @param solution Computed solution vector
 * @param grid_points Grid point coordinates (for FFT kernel)
 * @param grid_size Grid dimension n (n² = N)
 * @param verbose Print detailed output
 * @return Relative residual norm
 */
template<typename CoordType, typename DataType>
DataType verify_solution_fft(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    const std::vector<DataType>& solution,
    const std::vector<CoordType>& grid_points,
    int64_t grid_size,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank != 0) {
        return 0.0;  // Only verify on rank 0
    }
    
    int64_t N = solution.size();
    
    if (verbose) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Solution Verification (FFT)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Grid size: " << grid_size << " × " << grid_size << std::endl;
        std::cout << "Total DOFs: " << N << std::endl;
    }
    
    // Create FFT kernel
    kernel::LaplaceKernel3D_FFT<DataType> fft_kernel(
        grid_points.data(),
        grid_size,
        N
    );
    
    if (verbose) {
        std::cout << "FFT kernel memory: " 
                  << (fft_kernel.memory_usage() / 1024.0 / 1024.0) 
                  << " MB" << std::endl;
    }
    
    // Compute A*x using FFT
    std::vector<DataType> Ax(N);
    
    auto matvec_start = std::chrono::high_resolution_clock::now();
    fft_kernel.matvec(solution.data(), Ax.data());
    auto matvec_end = std::chrono::high_resolution_clock::now();
    auto matvec_duration = std::chrono::duration_cast<std::chrono::milliseconds>(matvec_end - matvec_start);
    
    if (verbose) {
        std::cout << "FFT matvec time: " << matvec_duration.count() << " ms" << std::endl;
    }
    
    // Compute residual: r = Ax - b
    std::vector<DataType> residual(N);
    for (int64_t i = 0; i < N; ++i) {
        residual[i] = Ax[i] - rhs[i];
    }
    
    // Compute norms
    DataType residual_norm = 0.0;
    DataType rhs_norm = 0.0;
    
    for (int64_t i = 0; i < N; ++i) {
        residual_norm += residual[i] * residual[i];
        rhs_norm += rhs[i] * rhs[i];
    }
    
    residual_norm = std::sqrt(residual_norm);
    rhs_norm = std::sqrt(rhs_norm);
    
    DataType relative_error = residual_norm / rhs_norm;
    
    if (verbose) {
        std::cout << "\nResidual Analysis:" << std::endl;
        std::cout << "  ||Ax - b||₂ = " << std::scientific << std::setprecision(6) 
                  << residual_norm << std::endl;
        std::cout << "  ||b||₂      = " << rhs_norm << std::endl;
        std::cout << "  Relative residual = " << relative_error << std::endl;
        
        // Additional statistics
        DataType max_residual = 0.0;
        int64_t max_idx = 0;
        for (int64_t i = 0; i < N; ++i) {
            DataType abs_res = std::abs(residual[i]);
            if (abs_res > max_residual) {
                max_residual = abs_res;
                max_idx = i;
            }
        }
        
        std::cout << "  Max residual    = " << max_residual 
                  << " (at index " << max_idx << ")" << std::endl;
        
        // Check if solution is reasonable
        if (relative_error < 1e-6) {
            std::cout << "\n  ✓ EXCELLENT: Relative error < 1e-6" << std::endl;
        } else if (relative_error < 1e-3) {
            std::cout << "\n  ✓ GOOD: Relative error < 1e-3" << std::endl;
        } else if (relative_error < 1e-1) {
            std::cout << "\n  ⚠ WARNING: Relative error > 1e-3" << std::endl;
        } else {
            std::cout << "\n  ✗ ERROR: Relative error > 1e-1 (solution likely incorrect)" << std::endl;
        }
        
        std::cout << "========================================\n" << std::endl;
    }
    
    return relative_error;
}



/**
 * @brief Build grid points for FFT verification
 * 
 * Extracts grid coordinates from leaf-level boxes across all processes.
 * Gathers to rank 0. Assumes uniform grid structure.
 */
template<typename CoordType, typename DataType>
std::vector<CoordType> build_grid_points_from_tree(
    ParallelTree<CoordType, DataType>* tree,
    int64_t grid_size) {
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int leaf_level = tree->num_levels - 1;
    auto& leaf_level_ref = tree->levels[leaf_level];
    int dimension = tree->dimension;
    
    // ===== Step 1: Collect local indices and coordinates =====
    
    std::vector<int64_t> local_indices;
    std::vector<CoordType> local_coords;
    
    for (int64_t box_idx = 0; box_idx < leaf_level_ref.num_boxes_local; ++box_idx) {
        auto& box = leaf_level_ref.local_boxes[box_idx];
        
        for (int64_t i = 0; i < box.num_points; ++i) {
            local_indices.push_back(box.point_indices[i]);
            
            // Extract coordinates (dimension values per point)
            for (int d = 0; d < dimension; ++d) {
                local_coords.push_back(box.point_coords[i * dimension + d]);
            }
        }
    }
    
    // ===== Step 2: Gather counts to root =====
    
    int local_count = static_cast<int>(local_indices.size());
    std::vector<int> recv_counts(rank == 0 ? size : 0);
    std::vector<int> recv_displs(rank == 0 ? size : 0);
    
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // ===== Step 3: Calculate displacements on rank 0 =====
    
    int total_count = 0;
    if (rank == 0) {
        recv_displs[0] = 0;
        for (int i = 0; i < size; ++i) {
            if (i > 0) {
                recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            }
            total_count += recv_counts[i];
        }
    }
    
    // ===== Step 4: Gather indices =====
    
    std::vector<int64_t> all_indices(rank == 0 ? total_count : 0);
    MPI_Gatherv(
        local_indices.data(), local_count, MPI_INT64_T,
        all_indices.data(), recv_counts.data(), recv_displs.data(), MPI_INT64_T,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 5: Gather coordinates =====
    
    MPI_Datatype mpi_coord_type;
    if constexpr (std::is_same_v<CoordType, double>) {
        mpi_coord_type = MPI_DOUBLE;
    } else if constexpr (std::is_same_v<CoordType, float>) {
        mpi_coord_type = MPI_FLOAT;
    } else {
        throw std::runtime_error("Unsupported CoordType for MPI");
    }
    
    // Adjust counts/displs for coordinate data (dimension values per point)
    std::vector<int> coord_recv_counts(rank == 0 ? size : 0);
    std::vector<int> coord_recv_displs(rank == 0 ? size : 0);
    
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            coord_recv_counts[i] = recv_counts[i] * dimension;
            coord_recv_displs[i] = recv_displs[i] * dimension;
        }
    }
    
    int local_coord_count = local_count * dimension;
    std::vector<CoordType> all_coords(rank == 0 ? total_count * dimension : 0);
    
    MPI_Gatherv(
        local_coords.data(), local_coord_count, mpi_coord_type,
        all_coords.data(), coord_recv_counts.data(), coord_recv_displs.data(), mpi_coord_type,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 6: Assemble grid points on rank 0 =====
    
    std::vector<CoordType> grid_points;
    
    if (rank == 0) {
        int64_t N = (dimension == 2) ? (grid_size * grid_size) : (grid_size * grid_size * grid_size);
        grid_points.resize(N * dimension, CoordType{0.0});
        
        for (int i = 0; i < total_count; ++i) {
            int64_t global_idx = all_indices[i];
            
            if (global_idx < 0 || global_idx >= N) {
                throw std::runtime_error(
                    "build_grid_points_from_tree: Invalid global index " + 
                    std::to_string(global_idx));
            }
            
            // Copy coordinates (interleaved: [x0, y0, x1, y1, ...] or [x0, y0, z0, x1, y1, z1, ...])
            for (int d = 0; d < dimension; ++d) {
                grid_points[global_idx * dimension + d] = all_coords[i * dimension + d];
            }
        }
        // print_vector(grid_points, "Grid Points", std::min<int64_t>(10, grid_points.size()));
    }
    
    
    return grid_points;
}


// /**
//  * @brief Serialize vector of BoxData into buffer for MPI communication
//  * @param boxes Vector of boxes to serialize
//  * @return Buffer containing serialized data
//  */
// template<typename CoordType, typename DataType>
// std::vector<char> serialize_boxes(const std::vector<BoxData<CoordType, DataType>>& boxes) {
//     // First pass: calculate total size needed
//     int64_t num_boxes = boxes.size();
//     int64_t total_size = sizeof(int64_t);  // Space for num_boxes
    
//     // Calculate size for each box (we need to do this by actually serializing, 
//     // or having a size calculation function - for now, use a temporary buffer approach)
//     std::vector<int64_t> box_sizes(num_boxes);
    
//     for (int64_t i = 0; i < num_boxes; ++i) {
//         // Estimate size conservatively (can be refined with actual size calculation)
//         // For now, serialize to temporary buffer to get exact size
//         std::vector<char> temp_buffer(1024 * 1024);  // 1MB temp buffer per box
//         char* end_ptr = serialize(boxes[i], temp_buffer.data());
//         box_sizes[i] = end_ptr - temp_buffer.data();
//         total_size += box_sizes[i];
//     }
    
//     // Allocate buffer
//     std::vector<char> buffer(total_size);
//     char* ptr = buffer.data();
    
//     // Write number of boxes
//     std::memcpy(ptr, &num_boxes, sizeof(int64_t));
//     ptr += sizeof(int64_t);
    
//     // Serialize each box
//     for (int64_t i = 0; i < num_boxes; ++i) {
//         ptr = serialize(boxes[i], ptr);
//     }
    
//     // Resize to actual used size
//     buffer.resize(ptr - buffer.data());
    
//     return buffer;
// }

/**
 * @brief Serialize vector of BoxData into buffer for MPI communication
 * @param boxes Vector of boxes to serialize
 * @return Buffer containing serialized data
 */
template<typename CoordType, typename DataType>
std::vector<char> serialize_boxes(const std::vector<BoxData<CoordType, DataType>>& boxes) {
    // First pass: calculate total size needed using get_serialized_size
    int64_t num_boxes = boxes.size();
    size_t total_size = sizeof(int64_t);  // Space for num_boxes
    
    std::vector<size_t> box_sizes(num_boxes);
    
    for (int64_t i = 0; i < num_boxes; ++i) {
        box_sizes[i] = get_serialized_size(boxes[i]);
        total_size += box_sizes[i];
    }
    
    // Allocate exact buffer size
    std::vector<char> buffer(total_size);
    char* ptr = buffer.data();
    
    // Write number of boxes
    std::memcpy(ptr, &num_boxes, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // Serialize each box
    for (int64_t i = 0; i < num_boxes; ++i) {
        char* ptr_before = ptr;
        ptr = serialize(boxes[i], ptr);
        
        // Verify we wrote the expected size
        size_t bytes_written = ptr - ptr_before;
        if (bytes_written != box_sizes[i]) {
            throw std::runtime_error("Box " + std::to_string(i) + 
                                   " size mismatch: expected " + 
                                   std::to_string(box_sizes[i]) + 
                                   " but wrote " + std::to_string(bytes_written));
        }
    }
    
    // Verify total size
    size_t actual_size = ptr - buffer.data();
    if (actual_size != total_size) {
        throw std::runtime_error("Total size mismatch: expected " + 
                               std::to_string(total_size) + 
                               " but wrote " + std::to_string(actual_size));
    }
    
    return buffer;
}

/**
 * @brief Deserialize vector of BoxData from buffer
 * @param buffer Buffer containing serialized data
 * @return Vector of deserialized boxes
 */
template<typename CoordType, typename DataType>
std::vector<BoxData<CoordType, DataType>> deserialize_boxes(const std::vector<char>& buffer) {
    const char* ptr = buffer.data();
    
    // Read number of boxes
    int64_t num_boxes = 0;
    std::memcpy(&num_boxes, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);
    
    // Deserialize each box
    std::vector<BoxData<CoordType, DataType>> boxes;
    boxes.reserve(num_boxes);
    
    for (int64_t i = 0; i < num_boxes; ++i) {
        BoxData<CoordType, DataType> box;
        ptr = deserialize(box, ptr);
        boxes.push_back(std::move(box));
    }
    
    return boxes;
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
 * @param factorization_method Cholesky or None
 * @param unit_proxy_points Unit sphere proxy points
 * @param num_proxy Number of proxy points
 * @param proxy_radius Proxy sphere radius multiplier
 * @param verbose Print progress messages
 */
template<typename CoordType, typename DataType, typename KernelType>
void hierarchical_factorization_parallel(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
    bool verbose = true) {
    using clock = std::chrono::high_resolution_clock;
    
    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    
    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    int num_children = (dimension == 2) ? 4 : 8;
    
    if (verbose && rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Factorization (Parallel)" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "MPI processes: " << size << std::endl;
        std::cout << "Levels: " << num_levels << std::endl;
        std::cout << "Leaf level: " << leaf_level << std::endl;
        std::cout << "Dimension: " << dimension << "D" << std::endl;
        std::cout << "Tolerance: " << tolerance << std::endl;
        std::cout << "Matrix property: " << (is_symmetric ? "Symmetric" : "Nonsymmetric") << std::endl;
        std::cout << "Factorization: " << (factorization_method == FactorizationMethod::CHOLESKY ? "Cholesky" : "None") << std::endl;
        printf("max thread: %d\n", omp_get_max_threads());
        std::cout << "========================================\n" << std::endl;
    }
    
    // ===== Main factorization loop: leaf level down to level 1 =====
    auto total_time = clock::now();
    auto segment_start = clock::now();
    clock::duration total_communication_time{};
    
    for (int current_level = leaf_level; current_level >= 1; current_level--) {
        auto level_start = std::chrono::high_resolution_clock::now();
        
        auto& level = tree->levels[current_level];
        
        if (!level.is_process_active) {
            if (verbose && rank == 0) {
                std::cout << "\n===== Level " << current_level << " (Process " << rank << " inactive) =====" << std::endl;
            }
            continue;
        }

        // declare locks
        std::unordered_map<int64_t, omp_lock_t*> box_locks;

         // initialize locks
        for (auto& box : level.local_boxes) {
            auto* lock = new omp_lock_t;
            omp_init_lock(lock);
            box_locks[box.morton_index] = lock;
        }
        for (auto& box : level.ghost_boxes) {
            auto* lock = new omp_lock_t;
            omp_init_lock(lock);
            box_locks[box.morton_index] = lock;
        }
        for (auto& box : level.ghost_and_assisting_boxes_for_solve) {
            auto* lock = new omp_lock_t;
            omp_init_lock(lock);
            box_locks[box.morton_index] = lock;
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
        
        if (verbose && rank == 0) {
            std::cout << "\n===== Level " << current_level << " =====" << std::endl;
            std::cout << "Active processes: " << level.num_active_processes << std::endl;
            std::cout << "Boxes per process: " << level.num_boxes_local << std::endl;
        }
        
        // ===== Step 1: Gather ghost and assisting boxes =====
        

        
        // ===== Step 2: Eliminate boxes in colored order =====
        
        if (current_level > 1) {
            // Regular levels: do ID and elimination in colored order

            auto elim_start = std::chrono::high_resolution_clock::now();

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

            PendingFactorUpdates<double> pending_updates;
            int boundary_count = 0;

            for (int counter = 0; counter < static_cast<int>(color_bins.size()); ++counter) {

                const int  color_id_mod    = counter % num_colors;
                const bool is_interior     = (counter >= interior_start_loc);

                // ----------------------------------------------------------------
                // Communication / transport step (single-threaded)
                // ----------------------------------------------------------------
                const auto comm_duration_raw = transport_and_apply_factor_updates_symmetric_onehop(
                    tree, current_level, kernel, pending_updates, false);
                total_communication_time += comm_duration_raw;
                update_neighbor_slicing_for_level(level, is_symmetric);
                auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                if (verbose && rank == 0) {
                    std::cout << "  Comm time: "
                            << comm_duration.count()
                            << " ms" << std::endl;
                }

                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (verbose && rank == 0) {
                    std::cout << "  Processing "
                            << (is_interior ? "interior sub-wave " : "boundary color ")
                            << color_id_mod
                            << " (" << color_list.size() << " boxes)..." << std::endl;
                }

                if (color_list.empty()) continue;

                // ----------------------------------------------------------------
                // Per-thread pending_updates: avoids contention on the shared
                // accumulator. Merged back after the parallel section.
                // ----------------------------------------------------------------
                const int max_threads = omp_get_max_threads();
                std::vector<PendingFactorUpdates<double>> thread_pending(max_threads);
                
                int local_boundary_count = 0;

                #pragma omp parallel reduction(+:local_boundary_count)
                {
                    const int tid = omp_get_thread_num();
                    // cpu_set_t cpuset;
                    // CPU_ZERO(&cpuset);
                    // CPU_SET(tid, &cpuset);  // pin thread tid to CPU tid

                    // // int rc = sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);

                    // // // verify
                    // int cpu = sched_getcpu();
                    // printf("thread %d pinned to CPU %d\n", tid, cpu);

                    FactorizationThreadScratch<CoordType, DataType> scratch;

                    #pragma omp for schedule(static)
                    for (size_t bi = 0; bi < color_list.size(); ++bi) {


                        try {
                            const int64_t morton_idx = color_list[bi];

                            BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
                            if (box_ptr == nullptr) box_ptr = level.find_ghost_box(morton_idx);
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

                            local_boundary_count += box.on_boundary;

                            compute_and_modify(dimension,
                                &box, level, kernel,
                                scratch,
                                tolerance, is_symmetric, is_hermitian,
                                &thread_pending[tid],
                                factorization_method
                            );

                        } catch (const std::exception& e) {
                            throw std::runtime_error(
                                "Error in OpenMP parallel section at level " + std::to_string(current_level) +
                                ", color_id_mod " + std::to_string(color_id_mod) +
                                ", morton_idx " + std::to_string(color_list[bi]) +
                                ": " + e.what()
                            );
                        }
                    }
                } // end parallel

                

                boundary_count += local_boundary_count;

                slice_far_field_blocks(level, is_symmetric, is_hermitian);
                // Merge per-thread pending_updates into global accumulator
                for (int t = 0; t < max_threads; ++t) {
                    merge_pending(pending_updates, thread_pending[t]);
                    clear_pending_factor_updates_memory(thread_pending[t]);
                }

                // ----------------------------------------------------------------
                // Mark assisting boxes as eliminated after this wave.
                // Boundary assisting boxes are marked by their boundary color wave.
                // Interior assisting boxes are marked by their interior sub-wave.
                // ----------------------------------------------------------------
                for (size_t bi = 0; bi < color_list.size(); ++bi) {
                    const int64_t morton_idx = color_list[bi];
                    level.eliminated_boxes.insert(morton_idx);
                }
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

                // Final transport after the very last wave
                if (counter == static_cast<int>(color_bins.size()) - 1) {
                    const auto final_comm_duration = transport_and_apply_factor_updates_symmetric_onehop(
                        tree, current_level, kernel, pending_updates, false);
                    total_communication_time += final_comm_duration;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                }
            }
            // // Regular levels: do ID and elimination in colored order
            
            // auto elim_start = std::chrono::high_resolution_clock::now();
            
            // std::vector<double> workspace;
            // std::vector<double> X_NN_full;
            // std::vector<double> sketch_storage;
            // PendingFactorUpdates<double> pending_updates;
            // int64_t ws_rows, ws_cols;
            
            // // Define color order: blue → orange → purple (3D) → green → interior
            // // New coloring scheme: 2^d colors (4 in 2D, 8 in 3D)
            // // color_id = morton_idx mod 2^d  (use bitmask since 2^d is power of 2)
            // const int num_colors = 1 << dimension;

            // // Bin ALL local boxes by color
            // // Boundary-only colored bins
            // std::vector<std::vector<int64_t>> color_bins(static_cast<size_t>(num_colors));

            // // Interior boxes (local, in increasing Morton order by construction)
            // std::vector<int64_t> interior_boxes;
            // interior_boxes.reserve(level.local_boxes.size());

            // for (int64_t local_idx = 0; local_idx < static_cast<int64_t>(level.local_boxes.size()); ++local_idx) {
            //     const int64_t morton_idx = level.local_morton_start + local_idx;
            //     const auto& box = level.local_boxes[static_cast<size_t>(local_idx)];

            //     if (box.on_boundary) {
            //         // 4-color in 2D / 8-color in 3D, based on morton % num_colors
            //         const int color_id = static_cast<int>(morton_idx & (num_colors - 1));
            //         color_bins[static_cast<size_t>(color_id)].push_back(morton_idx);
            //     } else {
            //         // Interior: keep a simple list (already in morton order)
            //         interior_boxes.push_back(morton_idx);
            //     }
            // }
            // color_bins.push_back(std::move(interior_boxes));
            // int64_t interior_loc = color_bins.size() - 1;

            // // for(auto &each_color : color_bins){
            // //     each_color.clear();
            // // }

            // // for (size_t c = 0; c < num_colors; ++c) {
            // //     int64_t start = (static_cast<int64_t>(c)     * level.local_boxes.size()) / static_cast<int64_t>(num_colors);
            // //     int64_t end   = (static_cast<int64_t>(c + 1) * level.local_boxes.size()) / static_cast<int64_t>(num_colors);

            // //     color_bins[c].reserve(static_cast<size_t>(end - start));
            // //     for (int64_t x = start; x < end; ++x) {
            // //         color_bins[c].push_back(x);
            // //     }
            // // }

            // // Process each color in order
            // int boundary_count = 0;
            // for (int counter = 0; counter < color_bins.size(); ++counter) {
            //     int color_id = counter;
                
                
            //     // If your algorithm requires communication/synchronization between colors,
            //     // this is the natural place to do it (per-color).
            //     // e.g. gather_ghost_and_assisting_boxes_factorization(...);
            //     // e.g. transport_and_apply_factor_updates_symmetric_onehop(...);
            //     auto comm_start = std::chrono::high_resolution_clock::now();
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel, pending_updates, false);
            //     update_neighbor_slicing_for_level(level, is_symmetric);
            //     auto comm_end = std::chrono::high_resolution_clock::now();
            //     auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_end - comm_start);
            //     // if(counter == 1 && rank == 0){
            //     //     fflush(stdout);
            //     //     assert(false);
            //     // }
            //     if (verbose && rank == 0) {
            //         std::cout << "  Ghost/assisting gather time: " << comm_duration.count() << " ms" << std::endl;
            //     }
            //     const auto& color_list = color_bins[static_cast<size_t>(color_id)];

            //     if (verbose && rank == 0) {
            //         std::cout << "  Processing color " << color_id
            //                 << " (" << color_list.size() << " boxes)..." << std::endl;
            //     }

            //     for (int64_t morton_idx : color_list) {
            //         // For this loop morton_idx should always be local, but keep the ghost fallback for safety
            //         BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);
            //         if (box_ptr == nullptr) {
            //             box_ptr = level.find_ghost_box(morton_idx);
            //         }
            //         if (box_ptr == nullptr) {
            //             throw std::runtime_error(
            //                 "hierarchical_factorization_parallel: Morton index " +
            //                 std::to_string(morton_idx) + " not found in local_boxes or ghost_boxes at level " +
            //                 std::to_string(current_level) + " (color: " + std::to_string(color_id) + ")");
            //         }

            //         auto& box = *box_ptr;

            //         gather_id_workspace(
            //             &box, level, kernel,
            //             unit_proxy_points.data(), num_proxy,
            //             proxy_radius, is_symmetric,
            //             workspace, ws_rows, ws_cols,
            //             0,  // DEBUG
            //             box.on_boundary
            //         );

            //         boundary_count += box.on_boundary;

            //         compute_and_modify(dimension,
            //             &box, level, kernel,
            //             workspace, ws_rows, ws_cols,
            //             tolerance, is_symmetric, is_hermitian, X_NN_full, sketch_storage, &pending_updates,
            //             factorization_method
            //         );
                    
                    
            //         // if(box.morton_index == 65){
            //         //     fflush(stdout);
            //         //     assert(false);
            //         // }
                    
            //     }

            //     // update eliminated_boxes to include the incoming assisting box ids, otherwise it's going to think that they aren't eliminated
            //     for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         if(counter == interior_loc){
            //             level.eliminated_boxes.insert(kv.first);
            //         }else{
            //             if ((kv.first % (1ULL << tree->dimension)) == color_id && level.assisting_boxes[kv.second].on_boundary){
            //                 level.eliminated_boxes.insert(kv.first);
            //             }
            //         }
                    
            //     }


            //     auto assist_start = std::chrono::high_resolution_clock::now();
                
                
            //     // verify_modified_interaction_symmetry(tree->levels[current_level], 1e-16, true);
            //     if(counter == color_bins.size() - 1) {
            //         transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel, pending_updates, false);
            //         update_neighbor_slicing_for_level(level, is_symmetric);
            //     }
                
            //     auto assist_end = std::chrono::high_resolution_clock::now();
            //     auto assist_duration = std::chrono::duration_cast<std::chrono::milliseconds>(assist_end - assist_start);
            //     if (verbose && rank == 0) {
            //         std::cout << "  Assisting boxes gather time: " << assist_duration.count() << " ms" << std::endl;
            //     }
                
                
            // }
            
            // int world_rank=0, world_size=0;
            // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            // for (int counter = 0; counter < num_colors; ++counter) {
            //     int color_id = perm[counter];
            //     if (current_level == 4) color_id = counter;

            //     const auto& color_list = color_bins[(size_t)color_id];

            //     // (Optional) flush anything already pending before starting this color
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel,
            //                                                         pending_updates, /*is_hermitian=*/false);
                

            //     for (int turn = 0; turn < level.num_active_processes; ++turn) {
                    
            //         // ---------------- compute phase: ONLY rank==turn ----------------
            //         if (world_rank == turn) {
            //             if (verbose) {
            //                 std::cout << "[rank " << world_rank << "] color " << color_id
            //                         << " turn " << turn << " computing...\n";
            //             }
                        
            //             for (int64_t morton_idx : color_list) {
            //                 BoxData<CoordType, DataType>* box_ptr = level.find_local_box(morton_idx);

            //                 // IMPORTANT: do NOT fall back to ghost here
            //                 if (!box_ptr) continue; // not owned by this rank

            //                 auto& box = *box_ptr;
            //                 printf("rank %d processing box id: %d, color: %d, num_points: %d\n", world_rank, box.morton_index, color_id, static_cast<int>(box.num_points));
            //                 fflush(stdout);
                           
            //                 gather_id_workspace(&box, level, kernel,
            //                                     unit_proxy_points.data(), num_proxy,
            //                                     proxy_radius, is_symmetric,
            //                                     workspace, ws_rows, ws_cols,
            //                                     0, box.on_boundary);

            //                 compute_and_modify(dimension, &box, level, kernel,
            //                                 workspace, ws_rows, ws_cols,
            //                                 tolerance, is_symmetric, is_hermitian,
            //                                 X_NN_full, sketch_storage,
            //                                 &pending_updates,
            //                                 factorization_method);
            //                 // if(current_level == 4 && color_id == 2 && box.morton_index == 22){
            //                 //     fflush(stdout);
            //                 //     assert(false);
            //                 // }
            //                 if(box.morton_index == 65){
            //                     print_workspace_stats(level.local_boxes[1].far_field_modified_interactions[level.local_boxes[1].far_field_interaction_map[21]].A_NS.data, "A_NSaaaaaaaaaa (after gather_id_workspace)");
            //                     fflush(stdout);
            //                     // assert(false);
            //                 }
                            
            //             }
                        

            //             // Optional: inspect what I'm about to send
            //             // print_pending_factor_updates(pending_updates, std::cout,
            //             //     ("rank " + std::to_string(world_rank) + " outgoing: ").c_str());
            //         }
            //         if(world_rank == 1 && turn != 0){
            //                     print_workspace_stats(level.local_boxes[1].far_field_modified_interactions[level.local_boxes[1].far_field_interaction_map[21]].A_NS.data, "A_NSaaaaaaaaaa (after gather_id_workspace)");
            //                     fflush(stdout);
                                
            //             }
                    
            //         // update eliminated_boxes to include the incoming assisting box ids, otherwise it's going to think that they aren't eliminated
            //         // for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         //     if ((kv.first % (1ULL << tree->dimension)) == color_id && (kv.first < ((turn + 1) * level.num_boxes_local))){
            //         //         level.eliminated_boxes.insert(kv.first);
            //         //     }
            //         // }
            //         if(current_level == 4){
            //             MPI_Barrier(tree->comm); // synchronize after each turn (optional, depends on algorithm requirements)
            //         }
                    

            //     }
                
            //     for(const auto& kv : level.assisting_box_points_for_kernel_evaluation){
            //         if ((kv.first % (1ULL << tree->dimension)) == color_id){
            //             level.eliminated_boxes.insert(kv.first);
            //         }
            //     }
 
                
            //     // ---------------- comm/apply phase: ALL ranks participate ----------------
            //     transport_and_apply_factor_updates_symmetric_onehop(tree, current_level, kernel,
            //                                                         pending_updates, /*is_hermitian=*/false);
                
            //     // If your slicing update depends on newly-arrived data, do it here (per turn)
            //     update_neighbor_slicing_for_level(level, is_symmetric);

                

            // }
            // printf("level: %d, on boundary sum: %d from rank: %d\n", current_level, boundary_count, rank);
            
            auto elim_end = std::chrono::high_resolution_clock::now();
            auto elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(elim_end - elim_start);
            
            if (verbose && rank == 0) {
                std::cout << "  Elimination time: " << elim_duration.count() << " ms" << std::endl;
                
                // Print compression statistics
                int64_t total_skeleton = 0;
                int64_t total_redundant = 0;
                for (const auto& box : level.local_boxes) {
                    total_skeleton += box.skeleton_indices.size();
                    total_redundant += box.redundant_indices.size();
                    // printf("  Box %d: skeleton size = %d, redundant size = %d, total = %d\n", 
                    //        box.morton_index, 
                    //        static_cast<int>(box.skeleton_indices.size()), 
                    //        static_cast<int>(box.redundant_indices.size()), 
                    //        static_cast<int>(box.skeleton_indices.size() + box.redundant_indices.size()));
                }
                
                double compression = static_cast<double>(total_skeleton) / (total_skeleton + total_redundant);
                std::cout << "  Compression ratio: " << compression << " (" << (compression * 100) << "%)" << std::endl;
                std::cout << "  Average skeleton size: " << static_cast<double>(total_skeleton) / level.local_boxes.size() << std::endl;
            }
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (verbose && rank == 0) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }
        
        // ===== Step 3: Gather assisting boxes post-elimination =====
        clear_ghosts(level);
        
        

        
        // ===== Step 4: Build parent level interactions =====
        
        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1) {
            if (verbose && rank == 0) {
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
        std::vector<BoxData<CoordType, DataType>> parent_boxes = build_parent_level_interactions<CoordType, DataType, KernelType>(
            level,          // child_level
            tree->levels[current_level - 1],   // parent_level
            dimension,
            is_symmetric,
            is_hermitian,
            kernel,
            tree->global_bounds
        );
        
        
        
        
        auto transition_end = std::chrono::high_resolution_clock::now();
        auto transition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(transition_end - transition_start);
        
        if (verbose && rank == 0) {
            std::cout << "  Level transition time: " << transition_duration.count() << " ms" << std::endl;
            std::cout << "  Parent boxes created: " << parent_boxes.size() << std::endl;
        }
        
        // ===== Step 5: Handle process reduction =====
        
        auto& parent_level = tree->levels[current_level - 1];
        
        // Error check: Level 2 → 1 should NOT trigger reduction
        if (current_level == 2 && level.parent_level_owner != rank) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Level 2 → 1 triggered process reduction! "
                "This is not allowed. parent_level_owner = " + std::to_string(level.parent_level_owner) +
                ", current rank = " + std::to_string(rank));
        }


        if (level.parent_level_owner == rank) {
            // This process owns the parent - no send/receive needed
            if (parent_level.children_senders.size() > 1) {
                // This process is receiving from children
                
                // if (verbose && rank == 0) {
                //     std::cout << "  Process " << rank << " receiving from " << parent_level.children_senders.size() 
                //               << " children processes" << std::endl;
                // }
                
                // Start with local parent boxes
                std::vector<BoxData<CoordType, DataType>> all_parent_boxes = std::move(parent_boxes);
                
                // Receive from each child in children_senders order
                for (int child_rank : parent_level.children_senders) {
                    if (child_rank == rank) {
                        // Already have local boxes
                        continue;
                    }
                    
                    // Receive serialized parent boxes from child_rank
                    std::vector<char> recv_buffer;
                    
                    // First receive size
                    int64_t buffer_size = 0;
                    MPI_Status status;
                    segment_start = clock::now();
                    MPI_Recv(&buffer_size, 1, MPI_INT64_T, child_rank, 0, tree->comm, &status);
                    
                    // Receive data
                    recv_buffer.resize(buffer_size);
                    MPI_Recv_large(recv_buffer.data(), buffer_size, MPI_CHAR, child_rank, 1, tree->comm, &status);
                    total_communication_time += (clock::now() - segment_start);
                    
                    // Deserialize parent boxes
                    std::vector<BoxData<CoordType, DataType>> child_parent_boxes = 
                        deserialize_boxes<CoordType, DataType>(recv_buffer);
                    
                    // Append to all_parent_boxes
                    all_parent_boxes.insert(
                        all_parent_boxes.end(),
                        std::make_move_iterator(child_parent_boxes.begin()),
                        std::make_move_iterator(child_parent_boxes.end())
                    );
                }
                
                // Assign concatenated result to parent_level.local_boxes
                parent_level.local_boxes = std::move(all_parent_boxes);
                
                
            } else {
                // No children to receive from - just use local parent boxes
                parent_level.local_boxes = std::move(parent_boxes);
            }
            
        } else {
            // This process does NOT own the parent - send to parent_level_owner
            

            
            // Serialize parent boxes
            std::vector<char> send_buffer = serialize_boxes(parent_boxes);
            
            // Send size first
            int64_t buffer_size = send_buffer.size();
            MPI_Send(&buffer_size, 1, MPI_INT64_T, level.parent_level_owner, 0, tree->comm);
            
            // Send data
            MPI_Send_large(send_buffer.data(), buffer_size, MPI_CHAR, level.parent_level_owner, 1, tree->comm);
            
            // This process will become inactive at parent level
            parent_level.local_boxes.clear();

            if (verbose) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size() 
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }
        
        // Clear modified interaction matrices to free memory
        clear_modified_interaction_matrices(level);

        // Teardown
        for (auto& [morton, lock] : level.box_locks) {
            omp_destroy_lock(lock);
            delete lock;
        }
        
        auto level_end = std::chrono::high_resolution_clock::now();
        auto level_duration = std::chrono::duration_cast<std::chrono::milliseconds>(level_end - level_start);
        
        if (verbose && rank == 0) {
            std::cout << "  Total level time: " << level_duration.count() << " ms" << std::endl;
        }
        
    }
    
    // ===== Special handling for level 0 (root) =====
    
    if (verbose && rank == 0) {
        std::cout << "\n===== Level 0 (Root) =====" << std::endl;
    }
    
    auto& root_level = tree->levels[0];
    
    if (!root_level.is_process_active || root_level.local_boxes.empty()) {
        if (verbose) {
            std::cout << "  Process " << rank << " has no boxes at root level" << std::endl;
        }
    } else {
        if (root_level.local_boxes.size() != 1) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Expected exactly 1 box at root level, got " +
                std::to_string(root_level.local_boxes.size()));
        }
        
        auto& root_box = root_level.local_boxes[0];
        
        if (verbose && rank == 0) {
            std::cout << "  Root box points: " << root_box.num_points << std::endl;
        }
        
        // At level 0, the assembled matrix is just the schur complement
        if (!root_box.schur_complement.is_allocated()) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Root box schur complement not allocated");
        }
        
        int64_t n = root_box.schur_complement.rows;
        
        if (verbose && rank == 0) {
            std::cout << "  Schur complement size: " << n << " × " << n << std::endl;
        }
        
        // Factorize the root schur complement for diagonal solve
        root_box.X_RR.allocate(n, n, MatrixStorage<DataType>::FULL);
        
        if (factorization_method == FactorizationMethod::CHOLESKY) {
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
            
            if (verbose && rank == 0) {
                std::cout << "  ✓ Root Cholesky factorization complete" << std::endl;
            }
            
        } else {
            // No factorization: just copy schur complement to X_RR
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR.format = MatrixStorage<DataType>::FULL;
            
            if (verbose && rank == 0) {
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
    
    MPI_Barrier(tree->comm);
    
    if (verbose && rank == 0) {
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - total_time);
        auto total_comm_time_milli = std::chrono::duration_cast<std::chrono::milliseconds>(total_communication_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ Hierarchical Factorization Complete" << std::endl;
        std::cout << "  total time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  communication time: " << total_comm_time_milli.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    
}


/**
 * @brief Hierarchical solve routine (parallel MPI version)
 */
template<typename CoordType, typename DataType>
void hierarchical_solve_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    using clock = std::chrono::high_resolution_clock;
    auto get_data_start = clock::now();
    clock::duration communication_total_forward{};
    clock::duration communication_total_backward{};
    
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    
    if (verbose && rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Solve (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
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

            
            
            
            // Set RHS from global vector (only at leaf level)
            if (level == leaf_level) {
                for (int64_t i = 0; i < box.num_points; ++i) {
                    // int64_t global_idx = box.point_indices[i];
                    // solve_box.right_side[i] = rhs[global_idx];
                    solve_box.right_side[i] = rhs[global_idx];
                    global_idx++;
                }
                solve_box.left_side = solve_box.right_side;
            }
        }
    }

    
    if (verbose && rank == 0) {
        std::cout << "Initialized solve data structures" << std::endl;
    }
    
    // ===== Phase 1: Forward Sweep (V^{-1}) with level transitions =====
    
    if (verbose && rank == 0) {
        std::cout << "\n===== Phase 1: Forward Sweep (V^{-1}) =====" << std::endl;
    }
    
    auto forward_start = std::chrono::high_resolution_clock::now();
    PendingSolveUpdates<DataType> pending_solve;

    for (int level = leaf_level; level >= 1; level--) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) {
            continue;
        }
        // Error check: level 2→1 should not trigger reduction
        if (level == 2) {
            auto& parent_level = tree->levels[1];
            if (tree_level.num_active_processes != parent_level.num_active_processes) {
                throw std::runtime_error("Reduction between level 2 and level 1 is not allowed");
            }
        }
        if (verbose && rank == 0) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << std::endl;
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

        if (level >= 2) {
            const int num_colors = 1 << tree->dimension;

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

            const int interior_start_loc = static_cast<int>(color_bins.size()); // == num_colors
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

                for (int64_t morton_idx : color_list) {
                    const int64_t local_idx = morton_idx - tree_level.local_morton_start;

                    apply_forward_elimination(
                        tree_level,
                        solve_data[level][static_cast<size_t>(local_idx)],
                        solve_data[level],
                        fmm::MatrixProperty::SYMMETRIC,
                        pending_solve,
                        /*is_ghost=*/false
                    );
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
    
    if (verbose && rank == 0) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }
    
    // ===== Phase 2: Diagonal Solve (D^{-1}) =====
    
    if (verbose && rank == 0) {
        std::cout << "\n===== Phase 2: Diagonal Solve (D^{-1}) =====" << std::endl;
    }
    
    auto diagonal_start = std::chrono::high_resolution_clock::now();
    
    int64_t num_diagonal_solves = 0;
    
    // Solve all X_RR blocks (levels N-1 down to 2)
    for (int level = leaf_level; level >= 2; level--) {
        auto& tree_level = tree->levels[level];
        // std::cout << "  Diagonal solves before from rank: " << rank << std::endl;
        if (!tree_level.is_process_active) {
            continue;
        }
        
        for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
            auto& box = tree_level.local_boxes[box_idx];
            auto& solve_box = solve_data[level][box_idx];
            
            if (box.redundant_indices.empty()) {
                continue;
            }
            // std::cout << "  Diagonal solves after from rank: " << rank << "with size: " << tree_level.num_boxes_local << std::endl;
            int64_t r = box.redundant_indices.size();
            
            std::vector<DataType> b_R(r);
            for (int64_t i = 0; i < r; ++i) {
                b_R[i] = solve_box.left_side[box.redundant_indices[i]];
            }
            
            if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                char uplo = 'L';
                int n = r, nrhs = 1;
                int lda = r, ldb = r, info = 0;

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
            }
            
            for (int64_t i = 0; i < r; ++i) {
                solve_box.left_side[box.redundant_indices[i]] = b_R[i];
            }
            
            num_diagonal_solves++;
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
    
    if (verbose) {
        std::cout << "  Diagonal solves: " << num_diagonal_solves  << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal solve time: " << diagonal_duration.count() << " ms" << std::endl;
    }
    
    
    // ===== Phase 3: Backward Sweep (W^{-1}) with level transitions =====
    
    if (verbose && rank == 0) {
        std::cout << "\n===== Phase 3: Backward Sweep (W^{-1}) =====" << std::endl;
    }
    
    auto backward_start = std::chrono::high_resolution_clock::now();
    
    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) {
            continue;
        }
        
        if (verbose && rank == 0) {
            std::cout << "  Level " << level << ": " << tree_level.num_boxes_local << " boxes" << std::endl;
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
        scatter_solution_to_children(
            tree_level,
            parent_level,
            solve_data[level],
            solve_data[level - 1],
            tree->dimension
        );
        communication_total_backward += (clock::now() - get_data_start);

        // for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
        //     auto& box = tree_level.local_boxes[box_idx];
        //     auto& solve_box = solve_data[level][box_idx];
            
        //     for(int gg = 0; gg < box.skeleton_indices.size(); gg++)
        //     {
        //         assert(solve_box.skeleton_indices[gg] == box.skeleton_indices[gg]);
        //     }
        // }
        
        if (verbose && rank == 0) {
            std::cout << "    ← Scattered from level " << (level - 1) << std::endl;
        }
        
        if (level >= 2) {
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

                // Reverse Morton order within each wave
                for (auto it = bins.rbegin(); it != bins.rend(); ++it) {
                    const int64_t morton_idx = *it;
                    const int64_t local_idx  = morton_idx - tree_level.local_morton_start;

                    apply_backward_substitution(
                        tree_level,
                        solve_data[level][static_cast<size_t>(local_idx)],
                        solve_data[level],
                        fmm::MatrixProperty::SYMMETRIC,
                        /*is_ghost=*/false
                    );
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
    
    if (verbose && rank == 0) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }
    
    
    if (verbose && rank == 0) {
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int leaf_level = tree->num_levels - 1;
    auto& leaf_level_ref = tree->levels[leaf_level];
    
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
    
    // ===== Step 2: Gather counts to root =====
    
    int local_count = static_cast<int>(local_indices.size());
    std::vector<int> recv_counts(rank == 0 ? size : 0);
    std::vector<int> recv_displs(rank == 0 ? size : 0);
    
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // ===== Step 3: Calculate displacements on rank 0 =====
    
    int total_count = 0;
    if (rank == 0) {
        recv_displs[0] = 0;
        for (int i = 0; i < size; ++i) {
            if (i > 0) {
                recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            }
            total_count += recv_counts[i];
        }
    }
    
    // ===== Step 4: Gather indices =====
    
    std::vector<int64_t> all_indices(rank == 0 ? total_count : 0);
    MPI_Gatherv(
        local_indices.data(), local_count, MPI_INT64_T,
        all_indices.data(), recv_counts.data(), recv_displs.data(), MPI_INT64_T,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 5: Gather solution values =====
    
    MPI_Datatype mpi_datatype;
    if constexpr (std::is_same_v<DataType, double>) {
        mpi_datatype = MPI_DOUBLE;
    } else if constexpr (std::is_same_v<DataType, float>) {
        mpi_datatype = MPI_FLOAT;
    } else {
        throw std::runtime_error("Unsupported DataType for MPI");
    }
    
    std::vector<DataType> all_solution(rank == 0 ? total_count : 0);
    MPI_Gatherv(
        local_solution.data(), local_count, mpi_datatype,
        all_solution.data(), recv_counts.data(), recv_displs.data(), mpi_datatype,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 6: Gather RHS values =====
    
    std::vector<DataType> all_rhs(rank == 0 ? total_count : 0);
    MPI_Gatherv(
        local_rhs.data(), local_count, mpi_datatype,
        all_rhs.data(), recv_counts.data(), recv_displs.data(), mpi_datatype,
        0, MPI_COMM_WORLD
    );
    
    // ===== Step 7: Assemble solution and RHS on rank 0 =====
    
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

} // namespace fmm


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] 
                  << " <num_levels> <N> <grid_size> <tolerance>\n";
        return 1;
    }

    int num_levels = std::atoi(argv[1]);
    int64_t N = std::atoll(argv[2]);
    int64_t grid_size = std::atoll(argv[3]);
    double tolerance = std::atof(argv[4]);
    int dimension = 3;
    
    try {
        if (rank == 0) {
            std::cout << "=== Hierarchical Factorization Test (3D Laplace Kernel) ===" << std::endl;
        }
        
        // Create a 3D tree

        double bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
        int64_t reduction_threshold = (dimension == 2) ? 64 : 4096;
     
        std::vector<int> idx_map;
        auto tree = fmm::create_uniform_tree<double, double>(
            nullptr, N, num_levels, bounds, dimension, MPI_COMM_WORLD, idx_map, reduction_threshold);
        
        int leaf_level = num_levels - 1;
        
        if (rank == 0) {
            std::cout << "\nTree Setup:" << std::endl;
            std::cout << "Levels: " << num_levels << std::endl;
            std::cout << "Leaf level: " << leaf_level << std::endl;
            std::cout << "Total points: " << N << std::endl;
            std::cout << "MPI processes: " << size << std::endl;
        }
        

        // ===== NEW: Create kernel and factorizer for proxy points =====
        
        kernel::LaplaceKernel3D<double> kernel(N);
        
        fmm::HierarchicalFactorization<double, double, kernel::LaplaceKernel3D<double>> 
            factorizer(N, fmm::MatrixProperty::SYMMETRIC, &kernel, dimension, 
                      fmm::FactorizationMethod::CHOLESKY, 256);
        
        const auto& unit_proxy = factorizer.get_unit_proxy_points();
        int num_proxy = factorizer.get_num_proxy_points();
        
        // ===== NEW: Run hierarchical factorization =====
        

        bool is_symmetric = true;
        bool is_hermitian = false;
        
        MPI_Barrier(MPI_COMM_WORLD);


        
        // Set to 1 thread for reproducibility
        // openblas_set_num_threads(4);
        // omp_set_num_threads(4);
        
        auto total_start = std::chrono::high_resolution_clock::now();
        
        fmm::hierarchical_factorization_parallel(
            tree,
            &kernel,
            tolerance,
            is_symmetric,
            is_hermitian,
            fmm::FactorizationMethod::CHOLESKY,
            unit_proxy,
            num_proxy,
            2.5,     // proxy_radius
            true     // verbose
        );
        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start);
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total factorization time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    
        
      
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        if (rank == 0) {
            std::cout << "\n=== Factorization Complete ===" << std::endl;
        }

        // ===== After factorization =====
        
        
        int64_t total_point_coords = 0;
        int64_t min_ind = N;
        int64_t max_ind = 0;
        for (const auto& box : tree->levels[num_levels - 1].local_boxes) {
            total_point_coords += box.point_coords.size();
            for(int i = 0; i < box.point_indices.size(); i++)
            {
                min_ind = std::min(min_ind, box.point_indices[i]);
                max_ind = std::max(max_ind, box.point_indices[i]);
            }
        }
        // // Generate random RHS
        std::vector<double> rhs(total_point_coords);

        std::mt19937_64 gen(12345);  // Mersenne Twister
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (int64_t i = 0; i < total_point_coords; ++i) {
            rhs[i] = dist(gen);
        }
        // std::vector<double> rhs(total_point_coords);

        // std::mt19937_64 gen(12345);  // Mersenne Twister
        // std::uniform_real_distribution<double> dist(0.0, 1.0);

        // // Calculate quarter size
        // int64_t quarter_size = total_point_coords / 4;

        // // Fill first quarter with random values
        // for (int64_t i = 0; i < quarter_size; ++i) {
        //     rhs[i] = dist(gen);
        // }

        // // Copy first quarter to the remaining 3 quarters
        // std::copy(rhs.begin(), rhs.begin() + quarter_size, rhs.begin() + quarter_size);
        // std::copy(rhs.begin(), rhs.begin() + quarter_size, rhs.begin() + 2 * quarter_size);
        // std::copy(rhs.begin(), rhs.begin() + quarter_size, rhs.begin() + 3 * quarter_size);

        // // Handle any remaining elements if total_point_coords is not divisible by 4
        // int64_t remaining = total_point_coords - 4 * quarter_size;
        // for (int64_t i = 0; i < remaining; ++i) {
        //     rhs[4 * quarter_size + i] = rhs[i];
        // }

        // ===== NEW: Print RHS =====
        // fmm::print_vector(rhs, "Right-hand side (b)", 10);
        
        // Solve
        std::vector<std::vector<fmm::SolveDataRequest<double, double>>> solve_data(num_levels);
        fmm::hierarchical_solve_parallel(tree, rhs, solve_data, true);
        

        // ===== NEW: Print solution =====
        // fmm::print_vector(solution, "Solution (x)", 10);

        // ===== NEW: Verify solution using FFT =====
        std::vector<double> solution;
        std::vector<double> aggregated_rhs;
        fmm::gather_solution_to_root(tree, solve_data, solution, aggregated_rhs);
        // if (rank == 0)
        // {
        //     fmm::print_vector(solution, "Gathered Solution (x)", 16384);
        // }
        

        
        // Build grid points from tree
        auto grid_points = fmm::build_grid_points_from_tree(tree, grid_size);
        
        if (!grid_points.empty()) {
            // Verify solution
            double relative_error = fmm::verify_solution_fft(
                tree,
                aggregated_rhs,
                solution,
                grid_points,
                grid_size,
                true  // verbose
            );
                // Verify with direct BLAS (slow but exact, only for small N)
        
            if(N <= 4096)
            {
                printf("verifying with direct matrix vector multiply since N <= 4096");
                double direct_error = fmm::verify_solution_direct(
                    &kernel, aggregated_rhs, solution, grid_points, 
                    tree->num_points, dimension, true
                );
                std::cout << "Direct error: " << direct_error << std::endl;
            }
           

            
            std::cout << "\nFinal relative error: " << std::scientific << std::setprecision(10)
                        << relative_error << std::endl;
            
        }
        
        
        
    
        
        delete tree;
        
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

  
    
    MPI_Finalize();
    return 0;
}
