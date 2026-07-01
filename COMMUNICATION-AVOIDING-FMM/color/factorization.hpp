#ifndef FACTORIZATION_HPP
#define FACTORIZATION_HPP

#include "tree_impl.hpp"
#include "kernel.hpp"
#include "morton.hpp"
#include "serialization.hpp"
#include "id_decomposition.hpp"
#include <mpi.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <complex>
#include <sstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include "blas_declare.hpp"

/* used blas functions */
// dpotrf_, dgetrf_, dgetri_, dgemm_/gemm_, trsm_

namespace fmm {

template <typename T>
bool is_nan(const T& val) {
    return std::isnan(val);
}

template <typename T>
bool is_nan(const std::complex<T>& val) {
    return std::isnan(val.real()) || std::isnan(val.imag());
}

template <typename T>
bool is_finite_value(const T& val) {
    return std::isfinite(val);
}

template <typename T>
bool is_finite_value(const std::complex<T>& val) {
    return std::isfinite(val.real()) && std::isfinite(val.imag());
}

template <typename T>
double value_abs(const T& val) {
    return std::abs(val);
}

template <typename T>
double value_sq_norm(const T& val) {
    const double mag = std::abs(val);
    return mag * mag;
}

template <typename T>
double value_sq_norm(const std::complex<T>& val) {
    return std::norm(val);
}

struct MatrixDiagnosticSummary {
    int64_t rows = 0;
    int64_t cols = 0;
    size_t expected_size = 0;
    size_t actual_size = 0;
    size_t finite_count = 0;
    size_t nan_count = 0;
    size_t inf_count = 0;
    size_t zero_columns = 0;
    double min_abs = 0.0;
    double max_abs = 0.0;
    double col_norm_min = 0.0;
    double col_norm_avg = 0.0;
    double col_norm_max = 0.0;
    int64_t col_norm_argmin = -1;
    int64_t col_norm_argmax = -1;
};

template <typename DataType>
MatrixDiagnosticSummary summarize_colmajor_matrix(
    const std::vector<DataType>& data,
    int64_t rows,
    int64_t cols) {
    MatrixDiagnosticSummary summary;
    summary.rows = rows;
    summary.cols = cols;
    summary.actual_size = data.size();

    if (rows <= 0 || cols <= 0) {
        return summary;
    }

    summary.expected_size = static_cast<size_t>(rows) * static_cast<size_t>(cols);

    const size_t usable_size = std::min(summary.expected_size, summary.actual_size);
    std::vector<long double> column_sq_norms(static_cast<size_t>(cols), 0.0L);

    bool saw_finite = false;
    double min_abs = std::numeric_limits<double>::infinity();
    double max_abs = 0.0;

    for (int64_t j = 0; j < cols; ++j) {
        const size_t col_base = static_cast<size_t>(j) * static_cast<size_t>(rows);
        if (col_base >= usable_size) {
            continue;
        }

        const int64_t available_rows = static_cast<int64_t>(
            std::min(static_cast<size_t>(rows), usable_size - col_base));
        for (int64_t i = 0; i < available_rows; ++i) {
            const auto& val = data[col_base + static_cast<size_t>(i)];
            if (is_nan(val)) {
                ++summary.nan_count;
                continue;
            }
            if (!is_finite_value(val)) {
                ++summary.inf_count;
                continue;
            }

            ++summary.finite_count;
            const double mag = value_abs(val);
            min_abs = std::min(min_abs, mag);
            max_abs = std::max(max_abs, mag);
            column_sq_norms[static_cast<size_t>(j)] += static_cast<long double>(value_sq_norm(val));
            saw_finite = true;
        }
    }

    if (saw_finite) {
        summary.min_abs = min_abs;
        summary.max_abs = max_abs;
    }

    if (!column_sq_norms.empty()) {
        long double norm_sum = 0.0L;
        bool saw_col_norm = false;
        double col_norm_min = std::numeric_limits<double>::infinity();
        double col_norm_max = 0.0;

        for (size_t j = 0; j < column_sq_norms.size(); ++j) {
            const double norm = std::sqrt(static_cast<double>(column_sq_norms[j]));
            norm_sum += static_cast<long double>(norm);
            if (norm == 0.0) {
                ++summary.zero_columns;
            }
            if (!saw_col_norm || norm < col_norm_min) {
                col_norm_min = norm;
                summary.col_norm_argmin = static_cast<int64_t>(j);
            }
            if (!saw_col_norm || norm > col_norm_max) {
                col_norm_max = norm;
                summary.col_norm_argmax = static_cast<int64_t>(j);
            }
            saw_col_norm = true;
        }

        if (saw_col_norm) {
            summary.col_norm_min = col_norm_min;
            summary.col_norm_max = col_norm_max;
            summary.col_norm_avg = static_cast<double>(norm_sum / static_cast<long double>(column_sq_norms.size()));
        }
    }

    return summary;
}

inline std::string format_matrix_diagnostic_summary(
    const std::string& label,
    const MatrixDiagnosticSummary& summary) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << label
        << ": rows=" << summary.rows
        << " cols=" << summary.cols
        << " expected_size=" << summary.expected_size
        << " actual_size=" << summary.actual_size
        << " finite=" << summary.finite_count
        << " nan=" << summary.nan_count
        << " inf=" << summary.inf_count
        << " min_abs=" << summary.min_abs
        << " max_abs=" << summary.max_abs
        << " col_norm[min/avg/max]="
        << summary.col_norm_min << "/"
        << summary.col_norm_avg << "/"
        << summary.col_norm_max
        << " col_argmin=" << summary.col_norm_argmin
        << " col_argmax=" << summary.col_norm_argmax
        << " zero_cols=" << summary.zero_columns;
    return oss.str();
}

/**
 * @brief Matrix symmetry property enumeration
 */
enum class MatrixProperty {
    SYMMETRIC,      ///< symmetric: A = A^T, may be complex
    HERMITIAN,      ///< Complex Hermitian: A = A^H
    NONSYMMETRIC    ///< No symmetry assumption
};

/**
 * @brief Factorization method for applying matrix inverses
 * 
 * Determines how to factorize and apply inverse of matrices like X_RR:
 * - CHOLESKY: Use Cholesky factorization (A = L*L^T or U^T*U)
 *             Requires symmetric/Hermitian positive definite
 *             Most efficient for SPD matrices
 * - LU: Use LU factorization with partial pivoting
 *       Works for general nonsymmetric matrices
 * - NONE: Store explicit inverse (A^{-1})
 *         Fastest solve but most memory
 */
enum class FactorizationMethod {
    CHOLESKY,    ///< Cholesky factorization (for SPD matrices)
    LU,          ///< LU factorization with pivoting
    COMPLEX_SYM, ///< Bunch-Kaufman for complex symmetric matrices (A = P*L*D*L^T*P^T)
    NONE         ///< Explicit inverse (no factorization)
};

template <typename DataType>
void solve_lu_factored_system_in_place(
    const MatrixStorage<DataType>& matrix,
    const std::vector<int>& pivots,
    char trans,
    DataType* rhs,
    int n,
    int nrhs,
    int ldb,
    const char* context) {
    if (n == 0 || nrhs == 0) {
        return;
    }
    if (matrix.rows != n || matrix.cols != n) {
        throw std::runtime_error(std::string(context) + ": LU solve dimension mismatch");
    }
    if (pivots.size() < static_cast<size_t>(n)) {
        throw std::runtime_error(std::string(context) + ": missing LU pivots");
    }

    int lda = static_cast<int>(matrix.lda);
    int info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        dgetrs_(&trans, &n, &nrhs,
                const_cast<double*>(matrix.data.data()), &lda,
                const_cast<int*>(pivots.data()),
                rhs, &ldb, &info);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgetrs_(&trans, &n, &nrhs,
                const_cast<std::complex<double>*>(matrix.data.data()), &lda,
                const_cast<int*>(pivots.data()),
                rhs, &ldb, &info);
    } else {
        static_assert(sizeof(DataType) == 0, "Unsupported DataType for LU solve");
    }

    if (info != 0) {
        throw std::runtime_error(
            std::string(context) + ": LU solve failed with INFO = " + std::to_string(info));
    }
}

template <typename DataType>
std::vector<DataType> transpose_colmajor_copy(
    const std::vector<DataType>& matrix,
    int64_t rows,
    int64_t cols) {
    std::vector<DataType> transposed(static_cast<size_t>(rows * cols));
    for (int64_t j = 0; j < cols; ++j) {
        for (int64_t i = 0; i < rows; ++i) {
            transposed[static_cast<size_t>(j + i * cols)] =
                matrix[static_cast<size_t>(i + j * rows)];
        }
    }
    return transposed;
}

template <typename DataType>
void transpose_colmajor_back(
    const std::vector<DataType>& transposed,
    std::vector<DataType>& matrix,
    int64_t rows,
    int64_t cols) {
    matrix.resize(static_cast<size_t>(rows * cols));
    for (int64_t j = 0; j < cols; ++j) {
        for (int64_t i = 0; i < rows; ++i) {
            matrix[static_cast<size_t>(i + j * rows)] =
                transposed[static_cast<size_t>(j + i * cols)];
        }
    }
}

template <typename DataType>
void apply_right_inverse_in_place(
    const MatrixStorage<DataType>& inverse_or_factor,
    const std::vector<int>& pivots,
    FactorizationMethod factorization_method,
    std::vector<DataType>& matrix,
    int64_t rows,
    int64_t cols,
    const char* context) {
    if (rows == 0 || cols == 0) {
        return;
    }
    if (inverse_or_factor.rows != cols || inverse_or_factor.cols != cols) {
        throw std::runtime_error(std::string(context) + ": right-inverse dimension mismatch");
    }

    if (factorization_method == FactorizationMethod::CHOLESKY) {
        char side = 'R';
        char uplo = 'L';
        char transa = 'T';
        char diag = 'N';
        int m = static_cast<int>(rows);
        int n = static_cast<int>(cols);
        int lda = static_cast<int>(inverse_or_factor.lda);
        DataType alpha = 1.0;

        trsm_(&side, &uplo, &transa, &diag,
              &m, &n, &alpha,
              inverse_or_factor.data.data(), &lda,
              matrix.data(), &m);

        transa = 'N';
        trsm_(&side, &uplo, &transa, &diag,
              &m, &n, &alpha,
              inverse_or_factor.data.data(), &lda,
              matrix.data(), &m);
    } else if (factorization_method == FactorizationMethod::LU) {
        std::vector<DataType> transposed = transpose_colmajor_copy(matrix, rows, cols);
        const int n = static_cast<int>(cols);
        const int nrhs = static_cast<int>(rows);
        const int ldb = static_cast<int>(cols);
        char trans = 'T';
        solve_lu_factored_system_in_place(
            inverse_or_factor, pivots, trans, transposed.data(), n, nrhs, ldb, context);
        transpose_colmajor_back(transposed, matrix, rows, cols);
    } else if (factorization_method == FactorizationMethod::COMPLEX_SYM) {
        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            std::vector<DataType> transposed = transpose_colmajor_copy(matrix, rows, cols);
            const int n = static_cast<int>(cols);
            const int nrhs = static_cast<int>(rows);
            const int ldb = n;
            char uplo = 'L';
            int info = 0;
            zsytrs_(&uplo, &n, &nrhs,
                    inverse_or_factor.data.data(), &n,
                    const_cast<int*>(pivots.data()),
                    transposed.data(), &ldb, &info);
            if (info != 0)
                throw std::runtime_error(std::string(context) +
                    ": zsytrs right-inverse failed with INFO = " + std::to_string(info));
            transpose_colmajor_back(transposed, matrix, rows, cols);
        } else {
            throw std::runtime_error("COMPLEX_SYM factorization only supported for complex<double>");
        }
    } else if (factorization_method == FactorizationMethod::NONE) {
        std::vector<DataType> temp = matrix;
        int m = static_cast<int>(rows);
        int n = static_cast<int>(cols);
        int kk = static_cast<int>(cols);
        int lda = static_cast<int>(inverse_or_factor.lda);
        DataType alpha = 1.0;
        DataType beta = 0.0;

        gemm_("N", "N", &m, &n, &kk,
              &alpha, temp.data(), &m,
              inverse_or_factor.data.data(), &lda,
              &beta, matrix.data(), &m);
    }
}

template <typename DataType>
void apply_left_inverse_in_place(
    const MatrixStorage<DataType>& inverse_or_factor,
    const std::vector<int>& pivots,
    FactorizationMethod factorization_method,
    std::vector<DataType>& matrix,
    int64_t rows,
    int64_t cols,
    const char* context) {
    if (rows == 0 || cols == 0) {
        return;
    }
    if (inverse_or_factor.rows != rows || inverse_or_factor.cols != rows) {
        throw std::runtime_error(std::string(context) + ": left-inverse dimension mismatch");
    }

    if (factorization_method == FactorizationMethod::CHOLESKY) {
        char side = 'L';
        char uplo = 'L';
        char transa = 'N';
        char diag = 'N';
        int m = static_cast<int>(rows);
        int n = static_cast<int>(cols);
        int lda = static_cast<int>(inverse_or_factor.lda);
        DataType alpha = 1.0;

        trsm_(&side, &uplo, &transa, &diag,
              &m, &n, &alpha,
              inverse_or_factor.data.data(), &lda,
              matrix.data(), &m);

        transa = 'T';
        trsm_(&side, &uplo, &transa, &diag,
              &m, &n, &alpha,
              inverse_or_factor.data.data(), &lda,
              matrix.data(), &m);
    } else if (factorization_method == FactorizationMethod::LU) {
        const int n = static_cast<int>(rows);
        const int nrhs = static_cast<int>(cols);
        const int ldb = static_cast<int>(rows);
        char trans = 'N';
        solve_lu_factored_system_in_place(
            inverse_or_factor, pivots, trans, matrix.data(), n, nrhs, ldb, context);
    } else if (factorization_method == FactorizationMethod::COMPLEX_SYM) {
        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            const int n = static_cast<int>(rows);
            const int nrhs = static_cast<int>(cols);
            const int ldb = n;
            char uplo = 'L';
            int info = 0;
            zsytrs_(&uplo, &n, &nrhs,
                    inverse_or_factor.data.data(), &n,
                    const_cast<int*>(pivots.data()),
                    matrix.data(), &ldb, &info);
            if (info != 0)
                throw std::runtime_error(std::string(context) +
                    ": zsytrs left-inverse failed with INFO = " + std::to_string(info));
        } else {
            throw std::runtime_error("COMPLEX_SYM factorization only supported for complex<double>");
        }
    } else if (factorization_method == FactorizationMethod::NONE) {
        std::vector<DataType> temp = matrix;
        int m = static_cast<int>(rows);
        int n = static_cast<int>(cols);
        int kk = static_cast<int>(rows);
        int lda = static_cast<int>(inverse_or_factor.lda);
        DataType alpha = 1.0;
        DataType beta = 0.0;

        gemm_("N", "N", &m, &n, &kk,
              &alpha, inverse_or_factor.data.data(), &lda,
              temp.data(), &kk,
              &beta, matrix.data(), &m);
    }
}



/**
 * @brief Hierarchical matrix factorization class
 * 
 * Performs RS factorization (Section 2 of FMM3D.pdf) with proxy-based
 * interpolative decomposition.
 * 
 * @tparam CoordType Coordinate type (float or double)
 * @tparam DataType Matrix data type (can be complex)
 */
template<typename CoordType, typename DataType>
class HierarchicalFactorization {
private:
    int64_t N_points;              ///< Total number of points
    MatrixProperty property;        ///< Matrix symmetry property
    FactorizationMethod factorization_method;  ///< Method for matrix factorization
    int dimension;                 ///< Spatial dimension (2 or 3)
    void (*kernel)(int*, int*, DataType*, void*); /// Kernel evaluation function
    // Proxy point configuration
    int num_proxy_points;          ///< Number of proxy points per box
    std::vector<CoordType> unit_proxy_points;  ///< Unit circle/sphere proxy points (dim × num_proxy)
    
    // Compression parameters
    double tolerance;            ///< ID compression tolerance
    CoordType proxy_radius_factor; ///< Proxy surface radius = factor × box_size
    
public:
    // In the public section of HierarchicalFactorization class:

    
    /**
     * @brief Constructor for hierarchical factorization
     * 
     * @param N Total number of points in the problem
     * @param prop Matrix property (symmetric, hermitian, or nonsymmetric)
     * @param kernel_func Kernel evaluator
     * @param dim Spatial dimension (2 or 3)
     * @param factorization_type Method for factorizing/inverting matrices (default: CHOLESKY)
     * @param num_proxy Number of proxy points (-1 uses default 32 for 2D, 256 for 3D)
     * @param tol Compression tolerance (default: 1e-6)
     * @param proxy_factor Proxy surface radius factor (default: 2.5)
     */
    HierarchicalFactorization(
        int64_t N,
        MatrixProperty prop,
        void (*kernel_func)(int*, int*, DataType*, void*),
        int dim,
        FactorizationMethod factorization_type = FactorizationMethod::CHOLESKY,
        int num_proxy = -1,
        double tol = 1e-6,
        CoordType proxy_factor = 2.5)
        : N_points(N), 
          property(prop), 
          factorization_method(factorization_type),
          kernel(kernel_func), 
          dimension(dim), 
          tolerance(tol),
          proxy_radius_factor(proxy_factor) {
        
        // if (kernel == nullptr) {
        //     throw std::invalid_argument("Kernel cannot be null");
        // }
        
        if (dim != 2 && dim != 3) {
            throw std::invalid_argument("Dimension must be 2 or 3");
        }

        if (num_proxy < -1) {
            throw std::invalid_argument(
                "num_proxy must be -1 (default), 0, or a positive integer");
        }
        if (dim == 3 && num_proxy == 1) {
            throw std::invalid_argument(
                "num_proxy=1 is invalid in 3D. Use 0 to disable proxies or >= 2.");
        }
        
        // Validate factorization method matches matrix property
        if (factorization_type == FactorizationMethod::CHOLESKY) {
            if (prop == MatrixProperty::NONSYMMETRIC) {
                throw std::invalid_argument(
                    "Cholesky factorization requires symmetric or Hermitian matrix property");
            }
        }
        
        // Set default proxy points based on dimension if not specified
        if (num_proxy == -1) {
            num_proxy_points = (dim == 2) ? 32 : 256;
        } else {
            num_proxy_points = num_proxy;
        }
        
        // Initialize unit proxy points
        initialize_proxy_points();
    }
    
    /**
     * @brief Initialize unit circle/sphere proxy points
     * 
     * For 2D: Uniform points on unit circle
     * For 3D: Fibonacci sphere or icosahedral subdivision
     */
    void initialize_proxy_points() {
        if (num_proxy_points == 0) {
            unit_proxy_points.clear();
            return;
        }

        unit_proxy_points.resize(dimension * num_proxy_points);
        
        if (dimension == 2) {
            // Uniform points on unit circle
            for (int i = 0; i < num_proxy_points; ++i) {
                CoordType theta = 2.0 * M_PI * i / num_proxy_points;
                unit_proxy_points[i * 2] = std::cos(theta);
                unit_proxy_points[i * 2 + 1] = std::sin(theta);
            }
        } else {
            // Fibonacci sphere for 3D
            CoordType phi = (1.0 + std::sqrt(5.0)) / 2.0;  // Golden ratio
            
            for (int i = 0; i < num_proxy_points; ++i) {
                CoordType y = 1.0 - 2.0 * i / (num_proxy_points - 1.0);
                CoordType radius = std::sqrt(1.0 - y * y);
                CoordType theta = 2.0 * M_PI * i / phi;
                
                unit_proxy_points[i * 3] = radius * std::cos(theta);
                unit_proxy_points[i * 3 + 1] = y;
                unit_proxy_points[i * 3 + 2] = radius * std::sin(theta);
            }
        }
    }
    
    /**
     * @brief Get scaled and centered proxy points for a box
     * 
     * @param box_center Center of the box [x, y, z]
     * @param box_size Half-width of the box
     * @return Proxy point coordinates (column-major: dim × num_proxy_points)
     */
    std::vector<CoordType> get_box_proxy_points(
        const CoordType* box_center,
        CoordType box_size) const {
        
        std::vector<CoordType> proxy_points(dimension * num_proxy_points);
        CoordType radius = proxy_radius_factor * box_size;
        
        for (int i = 0; i < num_proxy_points; ++i) {
            for (int d = 0; d < dimension; ++d) {
                proxy_points[i * dimension + d] = 
                    box_center[d] + radius * unit_proxy_points[i * dimension + d];
            }
        }
        
        return proxy_points;
    }
    
    /**
     * @brief Perform factorization at a given level
     * 
     * Implements Algorithm 1 from Section 2 (upward pass):
     * 1. Evaluate proxy interactions
     * 2. Compute interpolative decomposition (skeleton/redundant partition)
     * 3. Compress near-field interactions
     * 4. Form Schur complement
     * 
     * @param tree Parallel tree structure
     * @param level Level to factorize (0 = root, L-1 = leaf)
     */
    void factorize_level(
        ParallelTree<CoordType, DataType>* tree,
        int level);
    
    /**
     * @brief Perform complete hierarchical factorization
     * 
     * Factorizes from leaf level to root (upward pass).
     * Must be called after gather_ghost_and_assisting_boxes.
     * 
     * @param tree Parallel tree structure
     */
    void factorize(ParallelTree<CoordType, DataType>* tree) {
        // Factorize from leaf to root (upward pass)
        for (int level = tree->num_levels - 1; level >= 0; --level) {
            factorize_level(tree, level);
            
            // Synchronize between levels
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    
    // Getters
    MatrixProperty get_property() const { return property; }
    FactorizationMethod get_factorization_method() const { return factorization_method; }
    int64_t get_num_points() const { return N_points; }
    int get_dimension() const { return dimension; }
    double get_tolerance() const { return tolerance; }
    int get_num_proxy_points() const { return num_proxy_points; }
    
    /**
    * @brief Get unit proxy points (on unit circle/sphere)
    */
    const std::vector<CoordType>& get_unit_proxy_points() const {
        return unit_proxy_points;
    }

    /**
     * @brief Check if matrix is symmetric (real or complex)
     */
    bool is_symmetric() const {
        return property == MatrixProperty::SYMMETRIC || 
               property == MatrixProperty::HERMITIAN;
    }

    
};

/**
* @brief Extract submatrix from column-major matrix
*/
template<typename T>
void extract_submatrix(
    const T* A, int64_t lda,
    const std::vector<int64_t>& row_indices,
    const std::vector<int64_t>& col_indices,
    T* out, int64_t ldout) {
    
    int64_t nrows = row_indices.size();
    int64_t ncols = col_indices.size();
    
    for (int64_t j = 0; j < ncols; ++j) {
        int64_t src_col = col_indices[j];
        for (int64_t i = 0; i < nrows; ++i) {
            int64_t src_row = row_indices[i];
            out[i + j * ldout] = A[src_row + src_col * lda];
        }
    }
}


// /**
//  * @brief Helper: Get neighbor interaction block with skeleton slicing if needed
//  * 
//  * Handles the case where a neighbor has been eliminated and only skeleton DOFs remain.
//  * 
//  * @param neighbor_morton Morton index of neighbor
//  * @param stored_block The stored ModifiedBlock (either A_NS or A_SN)
//  * @param level Tree level for looking up neighbor boxes and elimination status
//  * @param is_A_NS true for A_NS (n_neighbor × m_box), false for A_SN (m_box × n_neighbor)
//  * @param box_size Number of points in current box (for column/row count)
//  * @return Properly sliced matrix data (column-major)
//  */
// template<typename CoordType, typename DataType>
// std::vector<DataType> get_sliced_neighbor_block(
//     int64_t neighbor_morton,
//     const ModifiedBlock<DataType>& stored_block,
//     TreeLevel<CoordType, DataType>& level,
//     bool is_A_NS,
//     int64_t box_size,
//     const std::vector<int64_t>* current_box_skeleton = nullptr) {
    
//     const auto& source_matrix = is_A_NS ? stored_block.A_NS : stored_block.A_SN;
    
//     if (!source_matrix.is_allocated()) {
//         throw std::runtime_error(
//             "get_sliced_neighbor_block: Matrix not allocated for neighbor " + 
//             std::to_string(neighbor_morton));
//     }
    
//     // Check if NEIGHBOR has been eliminated (for row/col slicing of neighbor dimension)
//     bool neighbor_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
//                                 level.eliminated_boxes.end());
    
//     BoxData<CoordType, DataType>* neighbor_box = nullptr;
//     if (neighbor_eliminated) {
//         neighbor_box = level.find_local_box(neighbor_morton);
//         if (neighbor_box == nullptr) {
//             auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
//             if (ghost_it != level.ghost_id_to_index.end()) {
//                 neighbor_box = &level.ghost_boxes[ghost_it->second];
//             }
//         }
        
//         if (neighbor_box == nullptr) {
//             throw std::runtime_error(
//                 "get_sliced_neighbor_block: Cannot find eliminated neighbor box " + 
//                 std::to_string(neighbor_morton));
//         }
//     }
    
//     if (is_A_NS) {
//         // A_NS is (n_neighbor × box_size)
//         int64_t stored_rows = source_matrix.rows;
//         int64_t stored_cols = source_matrix.cols;
        
//         // Determine final dimensions
//         // Row slicing: based on whether NEIGHBOR is eliminated
//         int64_t final_rows = neighbor_eliminated ? neighbor_box->skeleton_indices.size() : stored_rows;
//         // Column slicing: based on whether current_box_skeleton is PROVIDED (not eliminated_boxes!)
//         int64_t final_cols = current_box_skeleton ? current_box_skeleton->size() : stored_cols;
        
//         // If already correct size, return as-is
//         if (stored_rows == final_rows && stored_cols == final_cols) {
//             return source_matrix.data;
//         }
        
//         // Need to slice
//         std::vector<DataType> sliced(final_rows * final_cols);
        
//         // Row indices: slice if neighbor eliminated
//         std::vector<int64_t> row_indices = neighbor_eliminated ? 
//             neighbor_box->skeleton_indices : 
//             std::vector<int64_t>(stored_rows);
//         if (!neighbor_eliminated) {
//             std::iota(row_indices.begin(), row_indices.end(), 0);
//         }
        
//         // Column indices: slice if skeleton provided
//         std::vector<int64_t> col_indices = current_box_skeleton ? 
//             *current_box_skeleton : 
//             std::vector<int64_t>(stored_cols);
//         if (!current_box_skeleton) {
//             std::iota(col_indices.begin(), col_indices.end(), 0);
//         }
        
//         extract_submatrix(
//             source_matrix.data.data(), stored_rows,
//             row_indices,
//             col_indices,
//             sliced.data(), final_rows
//         );
        
//         return sliced;
        
//     } else {
//         // A_SN is (box_size × n_neighbor)
//         int64_t stored_rows = source_matrix.rows;
//         int64_t stored_cols = source_matrix.cols;
        
//         // Row slicing: based on whether current_box_skeleton is PROVIDED
//         int64_t final_rows = current_box_skeleton ? current_box_skeleton->size() : stored_rows;
//         // Column slicing: based on whether NEIGHBOR is eliminated
//         int64_t final_cols = neighbor_eliminated ? neighbor_box->skeleton_indices.size() : stored_cols;
        
//         if (stored_rows == final_rows && stored_cols == final_cols) {
//             return source_matrix.data;
//         }
        
//         std::vector<DataType> sliced(final_rows * final_cols);
        
//         // Row indices: slice if skeleton provided
//         std::vector<int64_t> row_indices = current_box_skeleton ? 
//             *current_box_skeleton : 
//             std::vector<int64_t>(stored_rows);
//         if (!current_box_skeleton) {
//             std::iota(row_indices.begin(), row_indices.end(), 0);
//         }
        
//         // Column indices: slice if neighbor eliminated
//         std::vector<int64_t> col_indices = neighbor_eliminated ? 
//             neighbor_box->skeleton_indices : 
//             std::vector<int64_t>(stored_cols);
//         if (!neighbor_eliminated) {
//             std::iota(col_indices.begin(), col_indices.end(), 0);
//         }
        
//         extract_submatrix(
//             source_matrix.data.data(), stored_rows,
//             row_indices,
//             col_indices,
//             sliced.data(), final_rows
//         );
        
//         return sliced;
//     }
// }

/**
 * @brief Helper: Get neighbor interaction block with skeleton slicing if needed
 * 
 * Handles the case where a neighbor has been eliminated and only skeleton DOFs remain.
 * 
 * @param neighbor_morton Morton index of neighbor
 * @param stored_block The stored ModifiedBlock (either A_NS or A_SN)
 * @param level Tree level for looking up neighbor boxes and elimination status
 * @param is_A_NS true for A_NS (n_neighbor × m_box), false for A_SN (m_box × n_neighbor)
 * @param box_size Number of points in current box (for column/row count)
 * @return Properly sliced matrix data (column-major)
 */
template<typename CoordType, typename DataType>
void get_sliced_neighbor_block_into(
    int64_t neighbor_morton,
    const ModifiedBlock<DataType>& stored_block,
    TreeLevel<CoordType, DataType>& level,
    bool is_A_NS,
    int64_t box_size,
    std::vector<DataType>& out)
{
    const auto& M = is_A_NS ? stored_block.A_NS : stored_block.A_SN;

    if (!M.is_allocated()) {
        throw std::runtime_error(
            "get_sliced_neighbor_block: Matrix not allocated for neighbor " +
            std::to_string(neighbor_morton));
    }

    auto pack_no_padding_into = [&](const MatrixStorage<DataType>& A) {
        out.resize(static_cast<size_t>(A.rows * A.cols));
        for (int64_t j = 0; j < A.cols; ++j)
            for (int64_t i = 0; i < A.rows; ++i)
                out[static_cast<size_t>(i + j * A.rows)] = A(i, j);
    };

    const bool neighbor_eliminated =
        (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

    // Validate the "box_size" side dimension
    if (is_A_NS) {
        if (M.cols != box_size) {
            throw std::runtime_error(
                "get_sliced_neighbor_block(A_NS): Column mismatch - stored=" +
                std::to_string(M.cols) + " expected=" + std::to_string(box_size) +
                " neighbor=" + std::to_string(neighbor_morton));
        }
    } else {
        if (M.rows != box_size) {
            throw std::runtime_error(
                "get_sliced_neighbor_block(A_SN): Row mismatch - stored=" +
                std::to_string(M.rows) + " expected=" + std::to_string(box_size) +
                " neighbor=" + std::to_string(neighbor_morton));
        }
    }

    if (!neighbor_eliminated) {
        pack_no_padding_into(M);
        return;
    }

    // ---- fetch skeleton indices: local -> ghost -> assisting ----
    const std::vector<int64_t>* skel = nullptr;

    if (auto* nb = level.find_local_box(neighbor_morton)) {
        skel = &nb->skeleton_indices;
    } else {
        auto git = level.ghost_id_to_index.find(neighbor_morton);
        if (git != level.ghost_id_to_index.end()) {
            skel = &level.ghost_boxes[static_cast<size_t>(git->second)].skeleton_indices;
        } else {
            auto ait = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
            if (ait != level.assisting_box_points_for_kernel_evaluation.end()) {
                const auto idx = static_cast<size_t>(ait->second);
                if (idx >= level.assisting_boxes.size())
                    throw std::runtime_error("get_sliced_neighbor_block: assisting index OOB");
                skel = &level.assisting_boxes[idx].skel_indices;
            }
        }
    }

    if (!skel) {
        throw std::runtime_error(
            "get_sliced_neighbor_block: Cannot find eliminated neighbor box (local/ghost/assist) " +
            std::to_string(neighbor_morton));
    }

    const int64_t skeleton_size = static_cast<int64_t>(skel->size());
    if (skeleton_size <= 0) {
        throw std::runtime_error(
            "get_sliced_neighbor_block: eliminated neighbor has empty skeleton_indices: " +
            std::to_string(neighbor_morton));
    }

    // The dimension corresponding to the neighbor DOFs in the stored block:
    // - A_NS: rows correspond to neighbor
    // - A_SN: cols correspond to neighbor
    const int64_t neighbor_dim = is_A_NS ? M.rows : M.cols;

    // If the stored block is already at skeleton resolution, do NOT interpret skel_indices
    // as indices into this matrix. Just return it.
    if (neighbor_dim == skeleton_size) {
        pack_no_padding_into(M);
        return;
    }

    // If stored neighbor dimension is smaller than skeleton size, state is inconsistent.
    if (neighbor_dim < skeleton_size) {
        throw std::runtime_error(
            "get_sliced_neighbor_block: Stored neighbor_dim < skeleton_size. neighbor=" +
            std::to_string(neighbor_morton) + " neighbor_dim=" + std::to_string(neighbor_dim) +
            " skeleton_size=" + std::to_string(skeleton_size));
    }

    // Only when we actually slice (neighbor_dim > skeleton_size) do we need index-value checks.
    {
        int64_t max_idx = -1;
        int64_t min_idx = std::numeric_limits<int64_t>::max();
        for (int64_t idx : *skel) { max_idx = std::max(max_idx, idx); min_idx = std::min(min_idx, idx); }

        if (min_idx < 0 || max_idx >= neighbor_dim) {
            throw std::runtime_error(
                "get_sliced_neighbor_block: skeleton index out of range for stored block. neighbor=" +
                std::to_string(neighbor_morton) + " min_skel=" + std::to_string(min_idx) +
                " max_skel=" + std::to_string(max_idx) +
                " stored_neighbor_dim=" + std::to_string(neighbor_dim) +
                " M.rows=" + std::to_string(M.rows) + " M.cols=" + std::to_string(M.cols));
        }
    }

    // ---- slice ----
    if (is_A_NS) {
        // Slice rows: (neighbor_dim x box_size) -> (skeleton_size x box_size)
        out.resize(static_cast<size_t>(skeleton_size * box_size));
        for (int64_t j = 0; j < box_size; ++j) {
            for (int64_t i = 0; i < skeleton_size; ++i) {
                out[static_cast<size_t>(i + j * skeleton_size)] =
                    M((*skel)[static_cast<size_t>(i)], j);
            }
        }
        return;
    } else {
        // Slice cols: (box_size x neighbor_dim) -> (box_size x skeleton_size)
        out.resize(static_cast<size_t>(box_size * skeleton_size));
        for (int64_t j = 0; j < skeleton_size; ++j) {
            const int64_t src_col = (*skel)[static_cast<size_t>(j)];
            for (int64_t i = 0; i < box_size; ++i) {
                out[static_cast<size_t>(i + j * box_size)] = M(i, src_col);
            }
        }
        return;
    }
}

template<typename CoordType, typename DataType>
std::vector<DataType> get_sliced_neighbor_block(
    int64_t neighbor_morton,
    const ModifiedBlock<DataType>& stored_block,
    TreeLevel<CoordType, DataType>& level,
    bool is_A_NS,
    int64_t box_size)
{
    std::vector<DataType> out;
    get_sliced_neighbor_block_into(
        neighbor_morton, stored_block, level, is_A_NS, box_size, out);
    return out;
}

/**
 * @brief Helper: Slice modified block in both dimensions
 * 
 * Handles double slicing:
 * 1. Slice columns/rows to current box's skeleton (just computed)
 * 2. Slice rows/columns to neighbor's skeleton (if eliminated)
 * 
 * @param stored_block The stored ModifiedBlock
 * @param level Tree level for neighbor lookup
 * @param neighbor_morton Neighbor box Morton ID
 * @param current_box_skeleton Current box's skeleton indices
 * @param is_A_NS true for A_NS (n_neighbor × m_box), false for A_SN (m_box × n_neighbor)
 * @param expected_neighbor_size Expected neighbor size after elimination check
 * @return Properly sliced matrix (column-major)
 */
/**
 * @brief Helper: Slice modified block in both dimensions (current box + neighbor box),
 *        with support for assisting boxes.
 *
 * Handles double slicing:
 *  1) Slice the current-box dimension to current_box_skeleton (just computed)
 *  2) Slice the neighbor-box dimension to neighbor skeleton IF the neighbor was eliminated
 *     AND the stored block is still at full neighbor resolution.
 *
 * This version supports obtaining neighbor skeleton indices from:
 *  - local boxes
 *  - ghost boxes
 *  - assisting boxes (level.assisting_box_points_for_kernel_evaluation + level.assisting_boxes[idx].skel_indices)
 *
 * It also avoids assuming MatrixStorage::lda == rows by packing matrices to contiguous
 * (rows*cols) when returning or when using them as intermediate buffers.
 */

template<typename CoordType, typename DataType>
std::vector<DataType> slice_modified_block_both_directions(
    const ModifiedBlock<DataType>& stored_block,
    TreeLevel<CoordType, DataType>& level,
    int64_t neighbor_morton,
    const std::vector<int64_t>& current_box_skeleton,
    bool is_A_NS,
    int64_t expected_neighbor_size)
{
    const auto& source = is_A_NS ? stored_block.A_NS : stored_block.A_SN;
    const int64_t k = static_cast<int64_t>(current_box_skeleton.size());

    if (!source.is_allocated()) {
        throw std::runtime_error("slice_modified_block_both_directions: Source not allocated");
    }
    if (k <= 0) {
        throw std::runtime_error("slice_modified_block_both_directions: empty current_box_skeleton");
    }
    if (expected_neighbor_size <= 0) {
        throw std::runtime_error("slice_modified_block_both_directions: expected_neighbor_size <= 0");
    }

    auto pack_no_padding = [&](const MatrixStorage<DataType>& A) -> std::vector<DataType> {
        std::vector<DataType> out(static_cast<size_t>(A.rows * A.cols));
        for (int64_t j = 0; j < A.cols; ++j) {
            for (int64_t i = 0; i < A.rows; ++i) {
                out[static_cast<size_t>(i + j * A.rows)] = A(i, j);
            }
        }
        return out;
    };

    const bool neighbor_eliminated =
        (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

    auto get_neighbor_skeleton = [&]() -> const std::vector<int64_t>* {
        if (auto* nb = level.find_local_box(neighbor_morton)) {
            return &nb->skeleton_indices;
        }
        auto git = level.ghost_id_to_index.find(neighbor_morton);
        if (git != level.ghost_id_to_index.end()) {
            return &level.ghost_boxes[static_cast<size_t>(git->second)].skeleton_indices;
        }
        auto ait = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
        if (ait != level.assisting_box_points_for_kernel_evaluation.end()) {
            const size_t idx = static_cast<size_t>(ait->second);
            if (idx >= level.assisting_boxes.size()) {
                throw std::runtime_error("slice_modified_block_both_directions: assisting index OOB");
            }
            return &level.assisting_boxes[idx].skel_indices;
        }
        return nullptr;
    };

    if (is_A_NS) {
        // A_NS: (n_neighbor × num_points_box), want: (n_neighbor_skel × k)
        const int64_t stored_rows = source.rows;
        const int64_t stored_cols = source.cols;
        const bool current_side_resolved = (stored_cols == k);
        const bool neighbor_side_resolved = (stored_rows == expected_neighbor_size);

        if (current_side_resolved && neighbor_side_resolved) {
            return pack_no_padding(source);
        }
        if (stored_rows < expected_neighbor_size) {
            throw std::runtime_error(
                "slice_modified_block_both_directions(A_NS): stored_rows < expected_neighbor_size. "
                "stored_rows=" + std::to_string(stored_rows) +
                " expected=" + std::to_string(expected_neighbor_size));
        }

        const std::vector<int64_t>* neighbor_skeleton = nullptr;
        if (!neighbor_side_resolved) {
            if (!neighbor_eliminated) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_NS): row size mismatch but neighbor not eliminated. "
                    "stored_rows=" + std::to_string(stored_rows) +
                    " expected=" + std::to_string(expected_neighbor_size) +
                    " neighbor=" + std::to_string(neighbor_morton));
            }

            neighbor_skeleton = get_neighbor_skeleton();
            if (!neighbor_skeleton) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_NS): cannot find neighbor skeleton (local/ghost/assist). "
                    "neighbor=" + std::to_string(neighbor_morton));
            }
            if (static_cast<int64_t>(neighbor_skeleton->size()) < expected_neighbor_size) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_NS): neighbor skeleton size < expected_neighbor_size. "
                    "skel_size=" + std::to_string(neighbor_skeleton->size()) +
                    " expected=" + std::to_string(expected_neighbor_size) +
                    " neighbor=" + std::to_string(neighbor_morton));
            }

            int64_t min_idx = std::numeric_limits<int64_t>::max();
            int64_t max_idx = -1;
            for (int64_t i = 0; i < expected_neighbor_size; ++i) {
                const int64_t idx = (*neighbor_skeleton)[static_cast<size_t>(i)];
                min_idx = std::min(min_idx, idx);
                max_idx = std::max(max_idx, idx);
            }
            if (min_idx < 0 || max_idx >= stored_rows) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_NS): neighbor skeleton index out of bounds. "
                    "min=" + std::to_string(min_idx) + " max=" + std::to_string(max_idx) +
                    " stored_rows=" + std::to_string(stored_rows) +
                    " neighbor=" + std::to_string(neighbor_morton));
            }
        }

        std::vector<DataType> result(static_cast<size_t>(expected_neighbor_size * k));
        for (int64_t j = 0; j < k; ++j) {
            const int64_t src_col = current_side_resolved ? j : current_box_skeleton[static_cast<size_t>(j)];
            if (src_col < 0 || src_col >= stored_cols) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_NS): current skeleton col out of bounds. "
                    "src_col=" + std::to_string(src_col) + " stored_cols=" + std::to_string(stored_cols));
            }
            for (int64_t i = 0; i < expected_neighbor_size; ++i) {
                const int64_t src_row = neighbor_side_resolved ? i : (*neighbor_skeleton)[static_cast<size_t>(i)];
                result[static_cast<size_t>(i + j * expected_neighbor_size)] = source(src_row, src_col);
            }
        }
        return result;
    }

    // A_SN: (num_points_box × n_neighbor), want: (k × n_neighbor_skel)
    const int64_t stored_rows = source.rows;
    const int64_t stored_cols = source.cols;
    const bool current_side_resolved = (stored_rows == k);
    const bool neighbor_side_resolved = (stored_cols == expected_neighbor_size);

    if (current_side_resolved && neighbor_side_resolved) {
        return pack_no_padding(source);
    }
    if (stored_cols < expected_neighbor_size) {
        throw std::runtime_error(
            "slice_modified_block_both_directions(A_SN): stored_cols < expected_neighbor_size. "
            "stored_cols=" + std::to_string(stored_cols) +
            " expected=" + std::to_string(expected_neighbor_size));
    }

    const std::vector<int64_t>* neighbor_skeleton = nullptr;
    if (!neighbor_side_resolved) {
        if (!neighbor_eliminated) {
            throw std::runtime_error(
                "slice_modified_block_both_directions(A_SN): column size mismatch but neighbor not eliminated. "
                "stored_cols=" + std::to_string(stored_cols) +
                " expected=" + std::to_string(expected_neighbor_size) +
                " neighbor=" + std::to_string(neighbor_morton));
        }

        neighbor_skeleton = get_neighbor_skeleton();
        if (!neighbor_skeleton) {
            throw std::runtime_error(
                "slice_modified_block_both_directions(A_SN): cannot find neighbor skeleton (local/ghost/assist). "
                "neighbor=" + std::to_string(neighbor_morton));
        }
        if (static_cast<int64_t>(neighbor_skeleton->size()) < expected_neighbor_size) {
            throw std::runtime_error(
                "slice_modified_block_both_directions(A_SN): neighbor skeleton size < expected_neighbor_size. "
                "skel_size=" + std::to_string(neighbor_skeleton->size()) +
                " expected=" + std::to_string(expected_neighbor_size) +
                " neighbor=" + std::to_string(neighbor_morton));
        }

        int64_t min_idx = std::numeric_limits<int64_t>::max();
        int64_t max_idx = -1;
        for (int64_t j = 0; j < expected_neighbor_size; ++j) {
            const int64_t idx = (*neighbor_skeleton)[static_cast<size_t>(j)];
            min_idx = std::min(min_idx, idx);
            max_idx = std::max(max_idx, idx);
        }
        if (min_idx < 0 || max_idx >= stored_cols) {
            throw std::runtime_error(
                "slice_modified_block_both_directions(A_SN): neighbor skeleton index out of bounds. "
                "min=" + std::to_string(min_idx) + " max=" + std::to_string(max_idx) +
                " stored_cols=" + std::to_string(stored_cols) +
                " neighbor=" + std::to_string(neighbor_morton));
        }
    }

    std::vector<DataType> result(static_cast<size_t>(k * expected_neighbor_size));
    for (int64_t j = 0; j < expected_neighbor_size; ++j) {
        const int64_t src_col = neighbor_side_resolved ? j : (*neighbor_skeleton)[static_cast<size_t>(j)];
        for (int64_t i = 0; i < k; ++i) {
            const int64_t src_row = current_side_resolved ? i : current_box_skeleton[static_cast<size_t>(i)];
            if (src_row < 0 || src_row >= stored_rows) {
                throw std::runtime_error(
                    "slice_modified_block_both_directions(A_SN): current skeleton row out of bounds. "
                    "src_row=" + std::to_string(src_row) + " stored_rows=" + std::to_string(stored_rows));
            }
            result[static_cast<size_t>(i + j * k)] = source(src_row, src_col);
        }
    }
    return result;
}



struct SymmetryErrorReport {
    double max_abs = 0.0;     // max |A_ij - A_ji|
    double rel_frob = 0.0;    // ||A-A^T||_F / ||A||_F
    int64_t imax = -1, jmax = -1;
    double a_ij = 0.0, a_ji = 0.0;

    double min_diag =  std::numeric_limits<double>::infinity();
    double max_diag = -std::numeric_limits<double>::infinity();
    bool has_nan_or_inf = false;
};

// A is column-major with leading dimension lda (stride between rows in a column).
static SymmetryErrorReport symmetry_error_colmajor_real(
    const double* A, int64_t rows, int64_t cols, int64_t lda)
{
    if (!A) throw std::runtime_error("symmetry_error_colmajor_real: A==nullptr");
    if (rows <= 0 || cols <= 0) throw std::runtime_error("symmetry_error_colmajor_real: bad dims");
    if (lda < rows) throw std::runtime_error("symmetry_error_colmajor_real: lda < rows");

    SymmetryErrorReport rep;

    const int64_t n = std::min(rows, cols);
    long double ss_diff = 0.0L;
    long double ss_A    = 0.0L;

    // diag stats + ||A||_F
    for (int64_t j = 0; j < cols; ++j) {
        for (int64_t i = 0; i < rows; ++i) {
            const double a = A[i + j * lda];
            if (!std::isfinite(a)) rep.has_nan_or_inf = true;
            ss_A += (long double)a * (long double)a;
        }
    }
    for (int64_t i = 0; i < n; ++i) {
        const double d = A[i + i * lda];
        rep.min_diag = std::min(rep.min_diag, d);
        rep.max_diag = std::max(rep.max_diag, d);
    }

    // symmetry error (only over the square overlap)
    for (int64_t j = 0; j < n; ++j) {
        for (int64_t i = j + 1; i < n; ++i) {
            const double aij = A[i + j * lda];
            const double aji = A[j + i * lda];
            const double diff = aij - aji;
            const double ad = std::abs(diff);

            ss_diff += (long double)diff * (long double)diff * 2.0L; // count (i,j) and (j,i)

            if (ad > rep.max_abs) {
                rep.max_abs = ad;
                rep.imax = i; rep.jmax = j;
                rep.a_ij = aij; rep.a_ji = aji;
            }
        }
    }

    const long double nA = std::sqrt(ss_A);
    const long double nd = std::sqrt(ss_diff);
    rep.rel_frob = (nA > 0.0L) ? (double)(nd / nA) : (nd == 0.0L ? 0.0 : std::numeric_limits<double>::infinity());
    return rep;
}

static void print_symmetry_report(
    const std::string& name,
    const SymmetryErrorReport& r,
    std::ostream& os)
{
    os << name << " symmetry check:\n";
    os << "  max |A_ij - A_ji| = " << r.max_abs
       << " at (i=" << r.imax << ", j=" << r.jmax << ")"
       << " with A_ij=" << r.a_ij << " A_ji=" << r.a_ji << "\n";
    os << "  ||A-A^T||_F / ||A||_F = " << r.rel_frob << "\n";
    os << "  diag min=" << r.min_diag << " max=" << r.max_diag << "\n";
    os << "  has NaN/Inf: " << (r.has_nan_or_inf ? "YES" : "NO") << "\n";
}

/*
struct BlockStats {
    double min =  std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    double max_abs = 0.0;

    // location of extrema in (row, col) within the *full* workspace matrix
    int64_t min_row = -1, min_col = -1;
    int64_t max_row = -1, max_col = -1;
    int64_t maxabs_row = -1, maxabs_col = -1;

    long double sum = 0.0L;
    long double sum_sq = 0.0L;

    bool has_nan = false;
    bool has_inf = false;

    int64_t nrows = 0;
    int64_t ncols = 0;
};

template<typename DataType>
static BlockStats stats_block_colmajor(const std::vector<double>& W,
                                       int64_t ld,          // workspace_rows
                                       int64_t row0,        // starting row of the block
                                       int64_t nrows,
                                       int64_t ncols)
{
    if (ld <= 0 || nrows < 0 || ncols < 0) {
        throw std::runtime_error("stats_block_colmajor: bad dims");
    }
    if (nrows == 0 || ncols == 0) {
        BlockStats s; s.nrows = nrows; s.ncols = ncols; return s;
    }
    const int64_t rows_total = ld;
    const int64_t cols_total = (int64_t)(W.size() / (size_t)ld); // best-effort; caller should pass correct ncols

    // Bounds check against vector size (assumes full matrix is ld x ncols)
    const size_t need = (size_t)(ld * ncols);
    if (W.size() < need) {
        throw std::runtime_error("stats_block_colmajor: W.size too small for ld*ncols");
    }
    if (row0 < 0 || row0 + nrows > rows_total) {
        throw std::runtime_error("stats_block_colmajor: row range out of bounds");
    }
    (void)cols_total;

    BlockStats s;
    s.nrows = nrows;
    s.ncols = ncols;

    for (int64_t col = 0; col < ncols; ++col) {
        const int64_t base = col * ld + row0;
        for (int64_t r = 0; r < nrows; ++r) {
            const double x = W[(size_t)(base + r)];

            if (std::isnan(x)) { s.has_nan = true; continue; }
            if (!std::isfinite(x)) { s.has_inf = true; continue; }

            if (x < s.min) { s.min = x; s.min_row = row0 + r; s.min_col = col; }
            if (x > s.max) { s.max = x; s.max_row = row0 + r; s.max_col = col; }

            const double ax = std::abs(x);
            if (ax > s.max_abs) { s.max_abs = ax; s.maxabs_row = row0 + r; s.maxabs_col = col; }

            s.sum += (long double)x;
            s.sum_sq += (long double)x * (long double)x;
        }
    }
    return s;
}

static void print_block_stats(const BlockStats& s, const std::string& label, std::ostream& os)
{
    os << label << ":\n";
    os << "  shape=" << s.nrows << " x " << s.ncols
       << "  has_nan=" << (s.has_nan ? "YES" : "NO")
       << "  has_inf=" << (s.has_inf ? "YES" : "NO") << "\n";

    const int64_t count = s.nrows * s.ncols;
    if (count <= 0) return;

    if (std::isfinite(s.min) && std::isfinite(s.max)) {
        const long double mean = s.sum / (long double)count;
        const long double var = std::max((long double)0.0,
                                         s.sum_sq / (long double)count - mean * mean);
        const long double rms = std::sqrt(s.sum_sq / (long double)count);

        os << "  min=" << s.min << " @ (row=" << s.min_row << ", col=" << s.min_col << ")\n";
        os << "  max=" << s.max << " @ (row=" << s.max_row << ", col=" << s.max_col << ")\n";
        os << "  max_abs=" << s.max_abs << " @ (row=" << s.maxabs_row << ", col=" << s.maxabs_col << ")\n";
        os << "  mean=" << (double)mean << "  rms=" << (double)rms
           << "  std=" << (double)std::sqrt(var) << "\n";
    } else {
        os << "  min/max not finite (block may be all NaN/Inf)\n";
    }
}

// Summarize each neighbor row-segment: rows [offset, offset+n_i)
// If you also have neighbor morton IDs, pass them for labeling.
template<typename DataType>
static void print_workspace_segment_stats(const std::vector<DataType>& workspace,
                                         int64_t workspace_rows,
                                         int64_t workspace_cols,
                                         const std::vector<int64_t>& neighbor_point_counts,
                                         const std::vector<int64_t>* neighbor_mortons = nullptr,
                                         int64_t start_row_offset = 0,
                                         std::ostream& os = std::cerr)
{
    if (workspace_rows <= 0 || workspace_cols <= 0)
        throw std::runtime_error("print_workspace_segment_stats: bad workspace dims");

    const size_t need = (size_t)(workspace_rows * workspace_cols);
    if (workspace.size() < need)
        throw std::runtime_error("print_workspace_segment_stats: workspace.size < rows*cols");

    // Optional sanity check on labeling vector
    if (neighbor_mortons && neighbor_mortons->size() != neighbor_point_counts.size())
        throw std::runtime_error("print_workspace_segment_stats: neighbor_mortons size mismatch");

    // Global stats for the used region
    os << "workspace global used-region stats:\n";
    auto global = stats_block_colmajor(workspace, workspace_rows, 0, workspace_rows, workspace_cols);
    print_block_stats(global, "  ALL rows", os);

    int64_t offset = start_row_offset;
    int64_t sum_rows = 0;

    os << "workspace per-neighbor row segments:\n";
    for (size_t i = 0; i < neighbor_point_counts.size(); ++i) {
        const int64_t nrows = neighbor_point_counts[i];
        if (nrows < 0) throw std::runtime_error("neighbor_point_counts has negative entry");

        const int64_t row0 = offset;
        const int64_t row1 = offset + nrows; // exclusive

        std::string label = "  seg[" + std::to_string(i) + "]";
        if (neighbor_mortons) label += " morton=" + std::to_string((*neighbor_mortons)[i]);
        label += " rows=[" + std::to_string(row0) + "," + std::to_string(row1) + ")";

        if (nrows == 0) {
            os << label << ": (empty)\n";
        } else {
            auto s = stats_block_colmajor(workspace, workspace_rows, row0, nrows, workspace_cols);
            print_block_stats(s, label, os);
        }

        offset += nrows;
        sum_rows += nrows;
    }

    os << "segment row accounting:\n"
       << "  start_row_offset=" << start_row_offset
       << " sum(neighbor_point_counts)=" << sum_rows
       << " end_offset=" << (start_row_offset + sum_rows)
       << " workspace_rows=" << workspace_rows << "\n";

    if (start_row_offset + sum_rows != workspace_rows) {
        os << "  WARNING: segments do not exactly cover workspace_rows\n";
    }
}
*/

enum class DeferredXnnTargetKind : uint8_t {
    SCHUR = 0,
    NEAR_A_NS = 1,
    FAR_A_NS = 2
};

struct DeferredXnnTargetKey {
    int64_t box_morton = -1;
    int64_t neighbor_morton = -1;
    DeferredXnnTargetKind kind = DeferredXnnTargetKind::SCHUR;

    bool operator==(const DeferredXnnTargetKey& other) const noexcept {
        return box_morton == other.box_morton &&
               neighbor_morton == other.neighbor_morton &&
               kind == other.kind;
    }
};

struct DeferredXnnTargetKeyHash {
    size_t operator()(const DeferredXnnTargetKey& key) const noexcept {
        size_t h = std::hash<int64_t>{}(key.box_morton);
        h ^= std::hash<int64_t>{}(key.neighbor_morton) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(static_cast<int>(key.kind)) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

template<typename CoordType, typename DataType>
struct DeferredXnnEndpoint {
    BoxData<CoordType, DataType>* box = nullptr;
    PointDataRequest<CoordType>* assisting = nullptr;
    bool is_local = false;
    bool is_ghost = false;
    bool is_assisting = false;
};

template<typename CoordType, typename DataType>
struct DeferredXnnPairTask {
    uint64_t pair_order = 0;
    int64_t LG_morton = -1;
    int64_t RG_morton = -1;
    int64_t n_LG = 0;
    int64_t n_RG = 0;
    int64_t row_offset = 0;
    int64_t col_offset = 0;
    bool is_diagonal = false;
    bool is_one_hop = false;
    bool LG_is_assisting = false;
    bool RG_is_assisting = false;
};

template<typename DataType>
struct DeferredXnnPendingUpdate {
    uint64_t source_morton = 0;
    uint64_t local_order = 0;
    DeferredXnnTargetKey target;
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<DataType> delta;
};

struct DeferredXnnBaseTask {
    DeferredXnnTargetKey target;
    int64_t rows = 0;
    int64_t cols = 0;
};

template<typename DataType>
struct DeferredXnnBaseMatrix {
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<DataType> data;
};

template<typename DataType>
struct DeferredXnnBoxBatch {
    int64_t box_morton = -1;
    std::vector<DeferredXnnPendingUpdate<DataType>> updates;
    std::unordered_map<DeferredXnnTargetKey,
                       DeferredXnnBaseMatrix<DataType>,
                       DeferredXnnTargetKeyHash> base_matrices;
};

template<typename DataType>
struct DeferredXnnBatchWorkspace {
    std::vector<DeferredXnnBaseTask> unique_missing_targets;
    std::vector<DeferredXnnBaseMatrix<DataType>> computed_bases;
    std::vector<DeferredXnnBoxBatch<DataType>> box_batches;
};

struct DeferredXnnOwnedRowBlock {
    uint64_t local_order = 0;
    int64_t neighbor_morton = -1;
    DeferredXnnTargetKind kind = DeferredXnnTargetKind::SCHUR;
    int64_t source_row_offset = 0;
    int64_t packed_row_offset = 0;  // Row offset inside the GEMM output buffer.
    int64_t rows = 0;
    // A deferred row block can contribute to the local candidate storage, to a
    // transported canonical ADD delta, or to both at once.
    bool accumulate_locally = false;
    bool emit_remote_add = false;
};

template<typename DataType>
struct DeferredXnnAccumulatedTarget {
    DeferredXnnTargetKey target;
    int64_t rows = 0;
    int64_t cols = 0;
    std::vector<DataType> data;
};

template<typename DataType>
struct DeferredXnnOwnerScratch {
    std::vector<DataType> packed_temp2_rows;
    std::vector<DataType> packed_updates;
    std::vector<DeferredXnnOwnedRowBlock> owned_row_blocks;
    // Accumulate all source contributions for the same destination block here so
    // the box state is materialized and written back only once per candidate.
    std::vector<DeferredXnnAccumulatedTarget<DataType>> accumulated_targets;
    std::unordered_map<DeferredXnnTargetKey, size_t, DeferredXnnTargetKeyHash>
        accumulated_target_indices;
    // Track which mirrored half-pairs have already been preallocated for the
    // current candidate box so repeated source boxes do not recheck them.
    std::unordered_set<DeferredXnnTargetKey, DeferredXnnTargetKeyHash>
        preallocated_mirror_targets;
};

template<typename CoordType, typename DataType>
DeferredXnnEndpoint<CoordType, DataType> resolve_deferred_xnn_endpoint(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton) {
    DeferredXnnEndpoint<CoordType, DataType> endpoint;
    endpoint.box = level.find_local_box(morton);
    if (endpoint.box == nullptr) {
        endpoint.box = level.find_ghost_box(morton);
        endpoint.is_ghost = (endpoint.box != nullptr);
    } else {
        endpoint.is_local = true;
    }
    if (endpoint.box == nullptr) {
        auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(morton);
        if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
            endpoint.assisting = &level.assisting_boxes[assist_it->second];
            endpoint.is_assisting = true;
        }
    }
    return endpoint;
}

template<typename CoordType, typename DataType>
BoxData<CoordType, DataType>* resolve_deferred_xnn_box(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton) {
    BoxData<CoordType, DataType>* box = level.find_local_box(morton);
    if (box == nullptr) {
        box = level.find_ghost_box(morton);
    }
    return box;
}

template<typename CoordType, typename DataType>
bool deferred_xnn_box_has_nonlocal_endpoint(
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level) {
    if (box == nullptr) {
        return false;
    }

    for (int64_t neighbor_morton : box->one_hop) {
        if (level.find_local_box(neighbor_morton) == nullptr) {
            return true;
        }
    }

    return false;
}

template<typename CoordType, typename DataType>
bool deferred_xnn_box_is_eliminated(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton) {
    return level.eliminated_boxes.find(morton) != level.eliminated_boxes.end();
}

template<typename CoordType, typename DataType>
std::vector<CoordType> gather_deferred_xnn_coords(
    const BoxData<CoordType, DataType>* box,
    bool eliminated,
    int dimension) {
    if (!eliminated) {
        return box->point_coords;
    }

    std::vector<CoordType> coords(box->skeleton_indices.size() * dimension);
    for (size_t i = 0; i < box->skeleton_indices.size(); ++i) {
        int64_t src_idx = box->skeleton_indices[i];
        for (int d = 0; d < dimension; ++d) {
            coords[i * dimension + d] = box->point_coords[src_idx * dimension + d];
        }
    }
    return coords;
}

template<typename CoordType, typename DataType>
bool deferred_xnn_target_exists(
    const DeferredXnnTargetKey& target,
    BoxData<CoordType, DataType>* target_box) {
    if (target.kind == DeferredXnnTargetKind::SCHUR) {
        return target_box->schur_complement.is_allocated();
    }

    const auto& interaction_map =
        (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_interaction_map :
            target_box->far_field_interaction_map;
    return interaction_map.find(target.neighbor_morton) != interaction_map.end();
}

template<typename CoordType, typename DataType, typename KernelType>
std::vector<DataType> compute_deferred_xnn_base_matrix(
    const DeferredXnnBaseTask& task,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    int dimension) {
    if (task.rows == 0 || task.cols == 0) {
        return {};
    }

    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, task.target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error("compute_deferred_xnn_base_matrix: target box missing");
    }

    const bool target_eliminated = deferred_xnn_box_is_eliminated(level, task.target.box_morton);
    std::vector<CoordType> target_coords =
        gather_deferred_xnn_coords(target_box, target_eliminated, dimension);

    if (static_cast<int64_t>(target_coords.size()) != task.cols * dimension) {
        throw std::runtime_error("compute_deferred_xnn_base_matrix: target coordinate size mismatch");
    }

    std::vector<DataType> base(static_cast<size_t>(task.rows * task.cols));

    if (task.target.kind == DeferredXnnTargetKind::SCHUR) {
        kernel->evaluate_block(
            target_coords.data(), task.rows,
            target_coords.data(), task.cols,
            base.data(), task.rows);
        return base;
    }

    DeferredXnnEndpoint<CoordType, DataType> neighbor =
        resolve_deferred_xnn_endpoint(level, task.target.neighbor_morton);
    if (neighbor.box == nullptr && !neighbor.is_assisting) {
        throw std::runtime_error("compute_deferred_xnn_base_matrix: neighbor missing");
    }

    std::vector<CoordType> neighbor_coords;
    const CoordType* neighbor_coords_ptr =
        coords_ptr_maybe_sliced(
            level,
            task.target.neighbor_morton,
            task.rows,
            dimension,
            /*box_or_null=*/neighbor.is_assisting ? nullptr : neighbor.box,
            /*assist_or_null=*/neighbor.is_assisting ? neighbor.assisting : nullptr,
            neighbor_coords);

    kernel->evaluate_block(
        neighbor_coords_ptr, task.rows,
        target_coords.data(), task.cols,
        base.data(), task.rows);
    return base;
}

template<typename DataType>
void accumulate_deferred_xnn_matrix_in_place(
    std::vector<DataType>& target,
    const std::vector<DataType>& delta) {
    if (target.size() != delta.size()) {
        throw std::runtime_error("accumulate_deferred_xnn_matrix_in_place: size mismatch");
    }

    // Isolated accumulation hook so the reduction type can be upgraded later.
    for (size_t i = 0; i < target.size(); ++i) {
        target[i] += delta[i];
    }
}

template<typename DataType>
void accumulate_deferred_xnn_matrix_slice_in_place(
    std::vector<DataType>& target,
    const DataType* source,
    int64_t source_ld,
    int64_t source_row_offset,
    int64_t rows,
    int64_t cols) {
    if (static_cast<int64_t>(target.size()) != rows * cols) {
        throw std::runtime_error("accumulate_deferred_xnn_matrix_slice_in_place: size mismatch");
    }

    for (int64_t col = 0; col < cols; ++col) {
        const DataType* src_col = source + source_row_offset + col * source_ld;
        DataType* dst_col = target.data() + col * rows;
        for (int64_t row = 0; row < rows; ++row) {
            dst_col[row] += src_col[row];
        }
    }
}

// Standard GEMM can read a contiguous row span directly from column-major storage.
// This helper lets the owner-side path skip packing when the selected rows already
// form one span in the deferred temp2 matrix.
inline bool deferred_xnn_owned_rows_form_single_span(
    const std::vector<DeferredXnnOwnedRowBlock>& blocks,
    int64_t& span_start,
    int64_t& span_rows) {
    if (blocks.empty()) {
        span_start = 0;
        span_rows = 0;
        return false;
    }

    span_start = blocks.front().source_row_offset;
    int64_t next_offset = span_start;
    for (const auto& block : blocks) {
        if (block.source_row_offset != next_offset) {
            span_rows = 0;
            return false;
        }
        next_offset += block.rows;
    }

    span_rows = next_offset - span_start;
    return true;
}

// Load the current destination block into a scratch matrix before grouped
// accumulation starts. If the block does not exist yet, synthesize the kernel
// base so later source contributions can be folded into the same buffer.
template<typename CoordType, typename DataType, typename KernelType>
std::vector<DataType> materialize_deferred_xnn_target_matrix_for_accumulation(
    const DeferredXnnTargetKey& target,
    int64_t rows,
    int64_t cols,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    int dimension) {
    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error(
            "materialize_deferred_xnn_target_matrix_for_accumulation: target box missing");
    }

    if (target.kind == DeferredXnnTargetKind::SCHUR) {
        if (target_box->schur_complement.is_allocated()) {
            if (target_box->schur_complement.rows != rows ||
                target_box->schur_complement.cols != cols) {
                throw std::runtime_error(
                    "materialize_deferred_xnn_target_matrix_for_accumulation: Schur dimension mismatch");
            }
            return target_box->schur_complement.data;
        }

        DeferredXnnBaseTask task;
        task.target = target;
        task.rows = rows;
        task.cols = cols;
        // A missing block starts from the kernel interaction before deferred
        // X_NN updates are added on top.
        return compute_deferred_xnn_base_matrix(task, level, kernel, dimension);
    }

    auto& interaction_map =
        (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_interaction_map :
            target_box->far_field_interaction_map;
    auto& modified_interactions =
        (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_modified_interactions :
            target_box->far_field_modified_interactions;

    auto it = interaction_map.find(target.neighbor_morton);
    if (it == interaction_map.end()) {
        DeferredXnnBaseTask task;
        task.target = target;
        task.rows = rows;
        task.cols = cols;
        // No modified block exists yet, so materialize the unsliced kernel block.
        return compute_deferred_xnn_base_matrix(task, level, kernel, dimension);
    }

    auto& block = modified_interactions[it->second];
    if (!block.A_NS.is_allocated()) {
        throw std::runtime_error(
            "materialize_deferred_xnn_target_matrix_for_accumulation: A_NS missing");
    }

    if (block.A_NS.rows != rows || block.A_NS.cols != cols) {
        std::ostringstream oss;
        oss << "materialize_deferred_xnn_target_matrix_for_accumulation: existing block dimension mismatch"
            << " target_box=" << target.box_morton
            << " neighbor_box=" << target.neighbor_morton
            << " target_kind="
            << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                    "NEAR_A_NS" : "FAR_A_NS")
            << " stored_rows=" << block.A_NS.rows
            << " stored_cols=" << block.A_NS.cols
            << " expected_rows=" << rows
            << " expected_cols=" << cols;
        throw std::runtime_error(oss.str());
    }

    const size_t expected_size = static_cast<size_t>(rows * cols);
    if (block.A_NS.data.size() != expected_size) {
        std::ostringstream oss;
        oss << "materialize_deferred_xnn_target_matrix_for_accumulation: existing block data size mismatch"
            << " target_box=" << target.box_morton
            << " neighbor_box=" << target.neighbor_morton
            << " target_kind="
            << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                    "NEAR_A_NS" : "FAR_A_NS")
            << " data_size=" << block.A_NS.data.size()
            << " expected_size=" << expected_size;
        throw std::runtime_error(oss.str());
    }

    return block.A_NS.data;
}

// Write one fully accumulated destination block back to box storage. At this
// point every contributing source for the current candidate box has already been
// folded into target_state.data in deterministic order. This owner-side path
// replaces the current stored view with the grouped accumulated view; the
// symmetric mirror overwrite below uses a stricter helper that validates and
// reuses existing mirrored storage when possible.
template<typename CoordType, typename DataType>
void flush_deferred_xnn_target_matrix_from_accumulation(
    DeferredXnnAccumulatedTarget<DataType>& target_state,
    TreeLevel<CoordType, DataType>& level) {
    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, target_state.target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error(
            "flush_deferred_xnn_target_matrix_from_accumulation: target box missing");
    }

    if (target_state.target.kind == DeferredXnnTargetKind::SCHUR) {
        target_box->schur_complement.set_owned(
            target_state.rows, target_state.cols, std::move(target_state.data),
            MatrixStorage<DataType>::FULL);
        return;
    }

    auto& interaction_map =
        (target_state.target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_interaction_map :
            target_box->far_field_interaction_map;
    auto& modified_interactions =
        (target_state.target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_modified_interactions :
            target_box->far_field_modified_interactions;

    auto it = interaction_map.find(target_state.target.neighbor_morton);
    if (it == interaction_map.end()) {
        ModifiedBlock<DataType> new_block;
        new_block.neighbor_morton = target_state.target.neighbor_morton;
        new_block.A_NS.set_owned(
            target_state.rows, target_state.cols, std::move(target_state.data),
            MatrixStorage<DataType>::FULL);

        int64_t new_idx = static_cast<int64_t>(modified_interactions.size());
        modified_interactions.push_back(std::move(new_block));
        interaction_map[target_state.target.neighbor_morton] = new_idx;
        return;
    }

    auto& block = modified_interactions[it->second];
    block.A_NS.set_owned(
        target_state.rows, target_state.cols, std::move(target_state.data),
        MatrixStorage<DataType>::FULL);
}

// Ensure that the non-owner half of a symmetric pair already has concrete
// storage before the mirror-copy phase starts. This keeps the later parallel
// mirror pass free of map insertions and vector growth.
template<typename CoordType, typename DataType>
void ensure_symmetric_owner_deferred_xnn_target_matrix_storage(
    const DeferredXnnTargetKey& target,
    int64_t rows,
    int64_t cols,
    TreeLevel<CoordType, DataType>& level) {
    if (target.kind == DeferredXnnTargetKind::SCHUR) {
        throw std::runtime_error(
            "ensure_symmetric_owner_deferred_xnn_target_matrix_storage: unexpected Schur target");
    }

    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error(
            "ensure_symmetric_owner_deferred_xnn_target_matrix_storage: target box missing");
    }

    auto& interaction_map =
        (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_interaction_map :
            target_box->far_field_interaction_map;
    auto& modified_interactions =
        (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_modified_interactions :
            target_box->far_field_modified_interactions;

    auto it = interaction_map.find(target.neighbor_morton);
    if (it == interaction_map.end()) {
        ModifiedBlock<DataType> new_block;
        new_block.neighbor_morton = target.neighbor_morton;
        new_block.A_NS.allocate(rows, cols, MatrixStorage<DataType>::FULL);

        int64_t new_idx = static_cast<int64_t>(modified_interactions.size());
        modified_interactions.push_back(std::move(new_block));
        interaction_map[target.neighbor_morton] = new_idx;
        return;
    }

    auto& block = modified_interactions[it->second];
    if (!block.A_NS.is_allocated()) {
        block.A_NS.allocate(rows, cols, MatrixStorage<DataType>::FULL);
        return;
    }

    if (block.A_NS.rows != rows || block.A_NS.cols != cols) {
        std::ostringstream oss;
        oss << "ensure_symmetric_owner_deferred_xnn_target_matrix_storage: preallocated mirror block dimension mismatch"
            << " target_box=" << target.box_morton
            << " neighbor_box=" << target.neighbor_morton
            << " target_kind="
            << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                    "NEAR_A_NS" : "FAR_A_NS")
            << " existing_rows=" << block.A_NS.rows
            << " existing_cols=" << block.A_NS.cols
            << " requested_rows=" << rows
            << " requested_cols=" << cols;
        throw std::runtime_error(oss.str());
    }

    const size_t expected_size = static_cast<size_t>(rows * cols);
    if (block.A_NS.data.size() != expected_size) {
        std::ostringstream oss;
        oss << "ensure_symmetric_owner_deferred_xnn_target_matrix_storage: preallocated mirror block data size mismatch"
            << " target_box=" << target.box_morton
            << " neighbor_box=" << target.neighbor_morton
            << " target_kind="
            << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                    "NEAR_A_NS" : "FAR_A_NS")
            << " existing_data_size=" << block.A_NS.data.size()
            << " expected_size=" << expected_size;
        throw std::runtime_error(oss.str());
    }
}

template<typename CoordType, typename DataType>
void apply_deferred_xnn_update(
    const DeferredXnnPendingUpdate<DataType>& update,
    TreeLevel<CoordType, DataType>& level,
    std::unordered_map<DeferredXnnTargetKey,
                       DeferredXnnBaseMatrix<DataType>,
                       DeferredXnnTargetKeyHash>& base_matrices) {
    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, update.target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error("apply_deferred_xnn_update: target box missing");
    }

    if (update.target.kind == DeferredXnnTargetKind::SCHUR) {
        if (!target_box->schur_complement.is_allocated()) {
            auto base_it = base_matrices.find(update.target);
            if (base_it == base_matrices.end()) {
                throw std::runtime_error("apply_deferred_xnn_update: missing base matrix for Schur target");
            }

            target_box->schur_complement.allocate(
                update.rows, update.cols, MatrixStorage<DataType>::FULL);
            target_box->schur_complement.data = std::move(base_it->second.data);
            base_matrices.erase(base_it);
        }

        if (target_box->schur_complement.rows != update.rows ||
            target_box->schur_complement.cols != update.cols) {
            throw std::runtime_error("apply_deferred_xnn_update: Schur dimension mismatch");
        }

        accumulate_deferred_xnn_matrix_in_place(
            target_box->schur_complement.data, update.delta);
        return;
    }

    auto& interaction_map =
        (update.target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_interaction_map :
            target_box->far_field_interaction_map;
    auto& modified_interactions =
        (update.target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
            target_box->near_field_modified_interactions :
            target_box->far_field_modified_interactions;

    auto it = interaction_map.find(update.target.neighbor_morton);
    if (it == interaction_map.end()) {
        auto base_it = base_matrices.find(update.target);
        if (base_it == base_matrices.end()) {
            throw std::runtime_error("apply_deferred_xnn_update: missing base matrix for interaction target");
        }

        ModifiedBlock<DataType> new_block;
        new_block.neighbor_morton = update.target.neighbor_morton;
        new_block.A_NS.set_owned(
            update.rows, update.cols, std::move(base_it->second.data),
            MatrixStorage<DataType>::FULL);
        accumulate_deferred_xnn_matrix_in_place(new_block.A_NS.data, update.delta);

        int64_t new_idx = static_cast<int64_t>(modified_interactions.size());
        modified_interactions.push_back(std::move(new_block));
        interaction_map[update.target.neighbor_morton] = new_idx;
        base_matrices.erase(base_it);
        return;
    }

    const bool target_eliminated =
        deferred_xnn_box_is_eliminated(level, update.target.box_morton);
    auto& block = modified_interactions[it->second];

    std::vector<DataType> current;
    try {
        if (target_eliminated) {
            current = slice_modified_block_both_directions<CoordType, DataType>(
                block, level, update.target.neighbor_morton,
                target_box->skeleton_indices, true, update.rows);
        } else {
            current = get_sliced_neighbor_block<CoordType, DataType>(
                update.target.neighbor_morton, block, level, true, target_box->num_points);
        }
    } catch (const std::runtime_error& e) {
        std::ostringstream oss;
        oss << "apply_deferred_xnn_update: block slicing failed"
            << " target_box=" << update.target.box_morton
            << " neighbor_box=" << update.target.neighbor_morton
            << " source_box=" << update.source_morton
            << " update_rows=" << update.rows
            << " update_cols=" << update.cols
            << " target_eliminated=" << target_eliminated
            << " target_num_points=" << target_box->num_points
            << " target_skeleton_size=" << target_box->skeleton_indices.size()
            << " stored_A_NS_rows=" << block.A_NS.rows
            << " stored_A_NS_cols=" << block.A_NS.cols
            << " cause={" << e.what() << "}";
        throw std::runtime_error(oss.str());
    }

    if (static_cast<int64_t>(current.size()) != update.rows * update.cols) {
        throw std::runtime_error("apply_deferred_xnn_update: sliced A_NS dimension mismatch");
    }

    accumulate_deferred_xnn_matrix_in_place(current, update.delta);
    block.A_NS.set_owned(
        update.rows, update.cols, std::move(current), MatrixStorage<DataType>::FULL);
}

inline bool deferred_xnn_boxes_are_one_hop(
    int dimension,
    int64_t lhs_morton,
    int64_t rhs_morton) {
    if (dimension == 2) {
        uint32_t lhs_x, lhs_y, rhs_x, rhs_y;
        morton::decode_2d(lhs_morton, lhs_x, lhs_y);
        morton::decode_2d(rhs_morton, rhs_x, rhs_y);
        return (std::abs(static_cast<int64_t>(lhs_x) - static_cast<int64_t>(rhs_x)) <= 1 &&
                std::abs(static_cast<int64_t>(lhs_y) - static_cast<int64_t>(rhs_y)) <= 1);
    }

    uint32_t lhs_x, lhs_y, lhs_z, rhs_x, rhs_y, rhs_z;
    morton::decode_3d(lhs_morton, lhs_x, lhs_y, lhs_z);
    morton::decode_3d(rhs_morton, rhs_x, rhs_y, rhs_z);
    return (std::abs(static_cast<int64_t>(lhs_x) - static_cast<int64_t>(rhs_x)) <= 1 &&
            std::abs(static_cast<int64_t>(lhs_y) - static_cast<int64_t>(rhs_y)) <= 1 &&
            std::abs(static_cast<int64_t>(lhs_z) - static_cast<int64_t>(rhs_z)) <= 1);
}

// Owner-side accumulation may still need to run on eliminated ghost boxes: the
// global ownership rule can assign them as the canonical source block for a
// local mirror copy on this rank. The owner phase therefore skips only missing
// or assisting boxes.
template<typename CoordType, typename DataType>
bool deferred_xnn_should_skip_owner_candidate_box(
    BoxData<CoordType, DataType>* candidate_box,
    TreeLevel<CoordType, DataType>& level) {
    if (candidate_box == nullptr) {
        return true;
    }

    DeferredXnnEndpoint<CoordType, DataType> endpoint =
        resolve_deferred_xnn_endpoint(level, candidate_box->morton_index);
    if (endpoint.box == nullptr || !endpoint.is_local) {
        return true;
    }

    return false;
}

// The mirror-writeback phase should follow the same box-availability rules as
// the owner phase. Even eliminated ghosts may be read again on this rank, so
// they still need their mirrored interaction view populated.
template<typename CoordType, typename DataType>
bool deferred_xnn_should_skip_mirror_candidate_box(
    BoxData<CoordType, DataType>* candidate_box,
    TreeLevel<CoordType, DataType>& level) {
    return deferred_xnn_should_skip_owner_candidate_box(candidate_box, level);
}

struct DeferredXnnOwnershipKey {
    int color_order = -1;
    int64_t morton = -1;
};

inline bool operator<(
    const DeferredXnnOwnershipKey& lhs,
    const DeferredXnnOwnershipKey& rhs) {
    if (lhs.color_order != rhs.color_order) {
        return lhs.color_order < rhs.color_order;
    }
    return lhs.morton < rhs.morton;
}

template<typename CoordType, typename DataType>
int deferred_xnn_color_order(
    TreeLevel<CoordType, DataType>& level,
    BoxData<CoordType, DataType>* box,
    int64_t morton) {
    if (std::binary_search(level.blue.begin(), level.blue.end(), morton)) {
        return 0;
    }
    if (std::binary_search(level.orange.begin(), level.orange.end(), morton)) {
        return 1;
    }
    if (level.dimension == 3 &&
        std::binary_search(level.purple.begin(), level.purple.end(), morton)) {
        return 2;
    }
    if (std::binary_search(level.green.begin(), level.green.end(), morton)) {
        return (level.dimension == 3) ? 3 : 2;
    }

    const int interior_color = (level.dimension == 3) ? 4 : 3;
    if (box != nullptr && !box->on_boundary) {
        return interior_color;
    }

    std::ostringstream oss;
    oss << "deferred_xnn_color_order: box missing color classification"
        << " box=" << morton
        << " level=" << level.level
        << " on_boundary=" << (box != nullptr ? box->on_boundary : -1);
    throw std::runtime_error(oss.str());
}

template<typename CoordType, typename DataType>
DeferredXnnOwnershipKey deferred_xnn_global_ownership_key_for_morton(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton) {
    if (level.num_active_processes <= 0) {
        throw std::runtime_error(
            "deferred_xnn_global_ownership_key_for_morton: invalid active process count");
    }

    const uint32_t grid_size = 1u << level.level;
    const uint32_t active_processes = static_cast<uint32_t>(level.num_active_processes);
    const uint32_t procs_per_dim =
        (level.dimension == 2) ?
            (1u << (__builtin_ctz(active_processes) / 2)) :
            (1u << (__builtin_ctz(active_processes) / 3));
    const uint32_t local_grid_size = grid_size / procs_per_dim;

    uint32_t proc_x = 0;
    uint32_t proc_y = 0;
    uint32_t proc_z = 0;
    int64_t owner_process_morton = -1;

    if (level.dimension == 2) {
        uint32_t x = 0;
        uint32_t y = 0;
        morton::decode_2d(static_cast<uint64_t>(morton), x, y);
        proc_x = x / local_grid_size;
        proc_y = y / local_grid_size;
        owner_process_morton = static_cast<int64_t>(morton::encode_2d(proc_x, proc_y));
    } else {
        uint32_t x = 0;
        uint32_t y = 0;
        uint32_t z = 0;
        morton::decode_3d(static_cast<uint64_t>(morton), x, y, z);
        proc_x = x / local_grid_size;
        proc_y = y / local_grid_size;
        proc_z = z / local_grid_size;
        owner_process_morton = static_cast<int64_t>(morton::encode_3d(proc_x, proc_y, proc_z));
    }

    std::array<uint32_t, 3> owner_offset = {
        proc_x * local_grid_size,
        proc_y * local_grid_size,
        proc_z * local_grid_size
    };

    DeferredXnnOwnershipKey key;
    if (is_blue_box(
            morton, level.level, level.dimension,
            grid_size, local_grid_size, owner_offset.data())) {
        key.color_order = 0;
    } else if (is_orange_box(
                   morton, level.level, level.dimension,
                   grid_size, local_grid_size, owner_offset.data())) {
        key.color_order = 1;
    } else if (level.dimension == 3 &&
               is_purple_box(
                   morton, level.level, grid_size,
                   local_grid_size, owner_offset.data())) {
        key.color_order = 2;
    } else {
        std::array<int32_t, 3> grid_coords = {0, 0, 0};
        if (level.dimension == 2) {
            uint32_t x = 0;
            uint32_t y = 0;
            morton::decode_2d(static_cast<uint64_t>(morton), x, y);
            grid_coords[0] = static_cast<int32_t>(x);
            grid_coords[1] = static_cast<int32_t>(y);
        } else {
            uint32_t x = 0;
            uint32_t y = 0;
            uint32_t z = 0;
            morton::decode_3d(static_cast<uint64_t>(morton), x, y, z);
            grid_coords[0] = static_cast<int32_t>(x);
            grid_coords[1] = static_cast<int32_t>(y);
            grid_coords[2] = static_cast<int32_t>(z);
        }

        const int64_t boxes_per_process = level.num_boxes_global / level.num_active_processes;
        const int64_t remainder = level.num_boxes_global % level.num_active_processes;

        int64_t owner_start = -1;
        int64_t owner_end = -1;
        if (owner_process_morton < remainder) {
            owner_start = owner_process_morton * (boxes_per_process + 1);
            owner_end = owner_start + boxes_per_process;
        } else {
            owner_start = remainder * (boxes_per_process + 1) +
                (owner_process_morton - remainder) * boxes_per_process;
            owner_end = owner_start + boxes_per_process - 1;
        }

        const bool on_boundary = check_boundary_condition(
            static_cast<uint64_t>(morton),
            grid_coords,
            grid_size,
            level.dimension,
            owner_start,
            owner_end);
        key.color_order = on_boundary ?
            ((level.dimension == 3) ? 3 : 2) :
            ((level.dimension == 3) ? 4 : 3);
    }

    key.morton = morton;
    return key;
}

// Retain the global ordering key helpers for debugging; the active deferred
// owner rule compares endpoints by color first and Morton second.
template<typename CoordType, typename DataType>
DeferredXnnOwnershipKey deferred_xnn_ownership_key(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton) {
    BoxData<CoordType, DataType>* box = resolve_deferred_xnn_box(level, morton);
    if (box == nullptr) {
        throw std::runtime_error(
            "deferred_xnn_ownership_key: box missing for ownership rule");
    }

    DeferredXnnOwnershipKey key;
    key.color_order = deferred_xnn_color_order(level, box, morton);
    key.morton = morton;
    return key;
}

template<typename CoordType, typename DataType>
bool deferred_xnn_first_box_owns_pair(
    TreeLevel<CoordType, DataType>& level,
    int64_t first_morton,
    int64_t second_morton,
    bool second_is_assisting) {
    (void)level;
    if (first_morton == second_morton) {
        return true;
    }
    // The color implementation has no ghost boxes. If the second endpoint is
    // assisting, then the first endpoint is the only local writable box here.
    if (second_is_assisting) {
        return true;
    }
    // Otherwise both endpoints are local on this rank, so break ties by Morton.
    return first_morton < second_morton;
}

template<typename DataType>
void deferred_xnn_accumulate_dense(
    DenseBlock<DataType>& dst,
    int64_t rows,
    int64_t cols,
    const std::vector<DataType>& delta) {
    if (delta.size() != static_cast<size_t>(rows * cols)) {
        throw std::runtime_error("deferred_xnn_accumulate_dense: delta size mismatch");
    }
    if (dst.data.empty()) {
        dst.rows = rows;
        dst.cols = cols;
        dst.data = delta;
        return;
    }
    if (dst.rows != rows || dst.cols != cols) {
        throw std::runtime_error("deferred_xnn_accumulate_dense: dimension mismatch");
    }
    for (size_t i = 0; i < dst.data.size(); ++i) {
        dst.data[i] += delta[i];
    }
}

template<typename DataType>
void deferred_xnn_accumulate_canonical_edge_delta_from_slice(
    PendingFactorUpdates<DataType>& pending,
    int64_t first_morton,
    int64_t second_morton,
    EdgeKind kind,
    const DataType* source,
    int64_t source_ld,
    int64_t source_row_offset,
    int64_t rows_second,
    int64_t cols_first) {
    if (rows_second == 0 || cols_first == 0) {
        return;
    }

    const int64_t lo = std::min(first_morton, second_morton);
    const int64_t hi = std::max(first_morton, second_morton);

    std::vector<DataType> delta_lo_hi;
    int64_t dst_rows = 0;
    int64_t dst_cols = 0;

    if (first_morton == lo) {
        dst_rows = rows_second;
        dst_cols = cols_first;
        delta_lo_hi.resize(static_cast<size_t>(dst_rows * dst_cols));
        for (int64_t col = 0; col < cols_first; ++col) {
            const DataType* src_col = source + source_row_offset + col * source_ld;
            DataType* dst_col = delta_lo_hi.data() + col * rows_second;
            std::copy(src_col, src_col + rows_second, dst_col);
        }
    } else {
        dst_rows = cols_first;
        dst_cols = rows_second;
        delta_lo_hi.resize(static_cast<size_t>(dst_rows * dst_cols));
        for (int64_t col = 0; col < rows_second; ++col) {
            for (int64_t row = 0; row < cols_first; ++row) {
                delta_lo_hi[static_cast<size_t>(row + col * dst_rows)] =
                    source[static_cast<size_t>(source_row_offset + col + row * source_ld)];
            }
        }
    }

    auto& dst = pending.accumulated_deltas[EdgeKey{lo, hi, kind}];
    deferred_xnn_accumulate_dense(dst, dst_rows, dst_cols, delta_lo_hi);
}

template<typename CoordType, typename DataType>
void deferred_xnn_cache_remote_source_pair_replace_from_source_box(
    BoxData<CoordType, DataType>* source_box,
    int64_t candidate_morton,
    int64_t n_candidate,
    PendingFactorUpdates<DataType>& pending) {
    if (source_box == nullptr) {
        throw std::runtime_error(
            "deferred_xnn_cache_remote_source_pair_replace_from_source_box: source box missing");
    }

    auto it = source_box->near_field_interaction_map.find(candidate_morton);
    if (it == source_box->near_field_interaction_map.end()) {
        throw std::runtime_error(
            "deferred_xnn_cache_remote_source_pair_replace_from_source_box: source pair missing");
    }

    const auto& block = source_box->near_field_modified_interactions[it->second];
    if (!block.A_NS.is_allocated()) {
        throw std::runtime_error(
            "deferred_xnn_cache_remote_source_pair_replace_from_source_box: source A_NS missing");
    }

    const int64_t k_source = static_cast<int64_t>(source_box->skeleton_indices.size());
    if (block.A_NS.rows != n_candidate || block.A_NS.cols != k_source) {
        std::ostringstream oss;
        oss << "deferred_xnn_cache_remote_source_pair_replace_from_source_box: unexpected source block dimensions"
            << " source_box=" << source_box->morton_index
            << " candidate_box=" << candidate_morton
            << " stored_rows=" << block.A_NS.rows
            << " stored_cols=" << block.A_NS.cols
            << " expected_rows=" << n_candidate
            << " expected_cols=" << k_source;
        throw std::runtime_error(oss.str());
    }

    std::vector<DataType> transpose(static_cast<size_t>(k_source * n_candidate));
    for (int64_t col = 0; col < k_source; ++col) {
        for (int64_t row = 0; row < n_candidate; ++row) {
            transpose[static_cast<size_t>(col + row * k_source)] =
                block.A_NS.data[static_cast<size_t>(row + col * n_candidate)];
        }
    }

    pending.replace_blocks[ReplaceKey{candidate_morton, source_box->morton_index}] =
        DenseBlock<DataType>{k_source, n_candidate, std::move(transpose)};
}

// Previous-wave source boxes contribute only to non-wave candidate endpoints.
// In the no-halo color path this includes both:
//   1. local boxes, which can receive lock-free deferred writes directly, and
//   2. assisting boxes, which only need deferred transport payloads.
template<typename CoordType, typename DataType>
void collect_owner_deferred_xnn_candidates_for_source_box(
    BoxData<CoordType, DataType>* source_box,
    TreeLevel<CoordType, DataType>& level,
    const std::unordered_set<int64_t>& wave_box_set,
    std::vector<int64_t>& candidate_boxes) {
    if (source_box == nullptr || source_box->deferred_xnn_temp2.empty()) {
        return;
    }

    for (int64_t candidate_morton : source_box->one_hop) {
        if (wave_box_set.find(candidate_morton) != wave_box_set.end()) {
            continue;
        }

        DeferredXnnEndpoint<CoordType, DataType> endpoint =
            resolve_deferred_xnn_endpoint(level, candidate_morton);
        if (endpoint.is_assisting) {
            candidate_boxes.push_back(candidate_morton);
            continue;
        }
        if (endpoint.box == nullptr || !endpoint.is_local) {
            continue;
        }

        candidate_boxes.push_back(candidate_morton);
    }
}

template<typename CoordType, typename DataType>
void ensure_symmetric_owner_deferred_source_pair_target_matrix_storage(
    const DeferredXnnTargetKey& target,
    int64_t rows,
    int64_t cols,
    TreeLevel<CoordType, DataType>& level) {
    if (target.kind != DeferredXnnTargetKind::NEAR_A_NS) {
        throw std::runtime_error(
            "ensure_symmetric_owner_deferred_source_pair_target_matrix_storage: unexpected target kind");
    }

    BoxData<CoordType, DataType>* target_box =
        resolve_deferred_xnn_box(level, target.box_morton);
    if (target_box == nullptr) {
        throw std::runtime_error(
            "ensure_symmetric_owner_deferred_source_pair_target_matrix_storage: target box missing");
    }

    auto& interaction_map = target_box->near_field_interaction_map;
    auto& modified_interactions = target_box->near_field_modified_interactions;

    auto it = interaction_map.find(target.neighbor_morton);
    if (it == interaction_map.end()) {
        ModifiedBlock<DataType> new_block;
        new_block.neighbor_morton = target.neighbor_morton;
        new_block.A_NS.allocate(rows, cols, MatrixStorage<DataType>::FULL);

        int64_t new_idx = static_cast<int64_t>(modified_interactions.size());
        modified_interactions.push_back(std::move(new_block));
        interaction_map[target.neighbor_morton] = new_idx;
        return;
    }

    auto& block = modified_interactions[it->second];
    block.A_NS.allocate(rows, cols, MatrixStorage<DataType>::FULL);
}

// For one candidate box, replay each previous-wave source in deterministic order.
// Each source still uses one GEMM, but the resulting slices are merged into
// per-target accumulators so repeated slice/writeback of the same block is
// avoided. The same pass also preallocates the non-owner mirror blocks on self,
// which lets the later symmetric copy run without mutating the interaction
// containers.
template<typename CoordType, typename DataType, typename KernelType>
void apply_owner_deferred_xnn_updates_for_candidate_box(
    int64_t candidate_morton,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    const std::unordered_set<int64_t>& wave_box_set,
    DeferredXnnOwnerScratch<DataType>& scratch,
    std::vector<DeferredXnnTargetKey>& mirror_targets,
    PendingFactorUpdates<DataType>* pending) {
    mirror_targets.clear();

    DeferredXnnEndpoint<CoordType, DataType> candidate_endpoint =
        resolve_deferred_xnn_endpoint(level, candidate_morton);
    BoxData<CoordType, DataType>* candidate_box = candidate_endpoint.box;
    const bool candidate_is_local = candidate_endpoint.is_local;
    const bool candidate_is_assisting = candidate_endpoint.is_assisting;

    if (!candidate_is_local && !candidate_is_assisting) {
        return;
    }
    if (candidate_is_local &&
        (candidate_box == nullptr || candidate_box->one_hop.empty() ||
         deferred_xnn_should_skip_owner_candidate_box(candidate_box, level))) {
        return;
    }

    const int dimension = candidate_is_local
        ? ((candidate_box->bounds[4] == candidate_box->bounds[5]) ? 2 : 3)
        : level.dimension;
    scratch.accumulated_targets.clear();
    scratch.accumulated_target_indices.clear();
    scratch.preallocated_mirror_targets.clear();

    std::vector<int64_t> candidate_sources;
    if (candidate_is_local) {
        for (int64_t source_morton : candidate_box->one_hop) {
            if (wave_box_set.find(source_morton) != wave_box_set.end()) {
                candidate_sources.push_back(source_morton);
            }
        }
    } else {
        for (int64_t source_morton : wave_box_set) {
            BoxData<CoordType, DataType>* source_box =
                resolve_deferred_xnn_box(level, source_morton);
            if (source_box == nullptr || source_box->deferred_xnn_temp2.empty()) {
                continue;
            }
            if (std::find(
                    source_box->one_hop.begin(),
                    source_box->one_hop.end(),
                    candidate_morton) != source_box->one_hop.end()) {
                candidate_sources.push_back(source_morton);
            }
        }
    }

    for (int64_t source_morton : candidate_sources) {
        // The source keeps both deferred temp2 and the original X_NR until this
        // post-wave pass consumes them.
        BoxData<CoordType, DataType>* source_box =
            resolve_deferred_xnn_box(level, source_morton);
        if (source_box == nullptr || source_box->deferred_xnn_temp2.empty()) {
            continue;
        }
        if (!source_box->X_NR.is_allocated()) {
            throw std::runtime_error(
                "apply_owner_deferred_xnn_updates_for_candidate_box: source X_NR missing");
        }

        const auto& source_neighbor_counts = source_box->deferred_xnn_neighbor_point_counts;
        if (source_neighbor_counts.size() != source_box->one_hop.size()) {
            throw std::runtime_error(
                "apply_owner_deferred_xnn_updates_for_candidate_box: neighbor count size mismatch");
        }

        size_t source_candidate_idx = source_box->one_hop.size();
        int64_t candidate_col_offset = 0;
        for (size_t idx = 0; idx < source_box->one_hop.size(); ++idx) {
            if (source_box->one_hop[idx] == candidate_morton) {
                source_candidate_idx = idx;
                break;
            }
            // Source X_NR stores neighbor blocks contiguously in one_hop order.
            candidate_col_offset += source_neighbor_counts[idx];
        }

        if (source_candidate_idx == source_box->one_hop.size()) {
            throw std::runtime_error(
                "apply_owner_deferred_xnn_updates_for_candidate_box: candidate not found in source one_hop");
        }

        const int64_t n_candidate = source_neighbor_counts[source_candidate_idx];
        if (n_candidate == 0) {
            continue;
        }

        // Step 7 is fully deferred in store=true mode. If the source pair target
        // is local, mirror the already-updated Step-5 block later without a GEMM.
        // If the source pair target is remote, serialize that transpose now as a
        // REPLACE payload for the target owner.
        if (candidate_is_local) {
            DeferredXnnTargetKey source_mirror_target;
            source_mirror_target.box_morton = candidate_morton;
            source_mirror_target.neighbor_morton = source_morton;
            source_mirror_target.kind = DeferredXnnTargetKind::NEAR_A_NS;
            auto source_mirror_insert =
                scratch.preallocated_mirror_targets.insert(source_mirror_target);
            if (source_mirror_insert.second) {
                ensure_symmetric_owner_deferred_source_pair_target_matrix_storage(
                    source_mirror_target,
                    static_cast<int64_t>(source_box->skeleton_indices.size()),
                    n_candidate,
                    level);
                mirror_targets.push_back(source_mirror_target);
            }
        } else {
            if (pending == nullptr) {
                throw std::runtime_error(
                    "apply_owner_deferred_xnn_updates_for_candidate_box: remote source pair requires pending updates");
            }
            deferred_xnn_cache_remote_source_pair_replace_from_source_box(
                source_box, candidate_morton, n_candidate, *pending);
        }

        const int64_t total_neighbor_points = source_box->X_NR.rows;
        const int64_t r = source_box->X_NR.cols;
        if (static_cast<int64_t>(source_box->deferred_xnn_temp2.size()) != total_neighbor_points * r) {
            throw std::runtime_error(
                "apply_owner_deferred_xnn_updates_for_candidate_box: deferred temp2 size mismatch");
        }

        scratch.owned_row_blocks.clear();
        int64_t packed_total_rows = 0;
        int64_t row_offset = 0;
        for (size_t row_idx = 0; row_idx < source_box->one_hop.size(); ++row_idx) {
            const int64_t neighbor_morton = source_box->one_hop[row_idx];
            const int64_t n_neighbor = source_neighbor_counts[row_idx];
            const bool is_diagonal = (neighbor_morton == candidate_morton);

            DeferredXnnTargetKind target_kind = DeferredXnnTargetKind::SCHUR;
            bool accumulate_locally = false;
            bool emit_remote_add = false;

            if (is_diagonal) {
                // The candidate Schur update is local only when the candidate is
                // local. A remote-only diagonal stays transport-only.
                accumulate_locally = candidate_is_local;
                emit_remote_add = !candidate_is_local;
            } else if (n_neighbor > 0) {
                DeferredXnnEndpoint<CoordType, DataType> neighbor_endpoint =
                    resolve_deferred_xnn_endpoint(level, neighbor_morton);
                if (neighbor_endpoint.box == nullptr && !neighbor_endpoint.is_assisting) {
                    throw std::runtime_error(
                        "apply_owner_deferred_xnn_updates_for_candidate_box: neighbor box missing");
                }

                target_kind = deferred_xnn_boxes_are_one_hop(
                    dimension, candidate_morton, neighbor_morton) ?
                        DeferredXnnTargetKind::NEAR_A_NS :
                        DeferredXnnTargetKind::FAR_A_NS;

                if (candidate_is_local) {
                    if (neighbor_endpoint.is_local) {
                        const bool candidate_owns_pair = deferred_xnn_first_box_owns_pair(
                            level,
                            candidate_morton,
                            neighbor_morton,
                            false);

                        if (!candidate_owns_pair) {
                            DeferredXnnTargetKey mirror_target;
                            mirror_target.box_morton = candidate_morton;
                            mirror_target.neighbor_morton = neighbor_morton;
                            mirror_target.kind = target_kind;

                            auto insert_result =
                                scratch.preallocated_mirror_targets.insert(mirror_target);
                            if (insert_result.second) {
                                ensure_symmetric_owner_deferred_xnn_target_matrix_storage(
                                    mirror_target, n_neighbor, n_candidate, level);
                                mirror_targets.push_back(mirror_target);
                            }
                        }

                        // Case 1: local-local pair. The canonical local owner keeps
                        // the block here and the mirror pass reconstructs the other
                        // local endpoint later.
                        accumulate_locally = candidate_owns_pair;
                    } else {
                        // Case 2: local-remote pair. The local candidate owns the
                        // local storage update regardless of Morton order, and the
                        // same slice is also exported as a canonical remote ADD.
                        accumulate_locally = true;
                        emit_remote_add = true;
                    }
                } else {
                    if (!neighbor_endpoint.is_local) {
                        // Case 3: remote-remote pair. There is no local storage to
                        // update, so only one assisting endpoint emits the canonical
                        // ADD. Smaller Morton is the deterministic emitter.
                        emit_remote_add = (candidate_morton < neighbor_morton);
                    }
                    // candidate remote + neighbor local: skip here because the local
                    // candidate will rebuild its own storage and emit the remote ADD.
                }
            }

            if ((accumulate_locally || emit_remote_add) && n_neighbor > 0) {
                DeferredXnnOwnedRowBlock block;
                block.local_order = static_cast<uint64_t>(row_idx);
                block.neighbor_morton = neighbor_morton;
                block.kind = is_diagonal ?
                    DeferredXnnTargetKind::SCHUR : target_kind;
                block.source_row_offset = row_offset;
                block.packed_row_offset = packed_total_rows;
                block.rows = n_neighbor;
                block.accumulate_locally = accumulate_locally;
                block.emit_remote_add = emit_remote_add;
                scratch.owned_row_blocks.push_back(block);
                packed_total_rows += n_neighbor;
            }

            row_offset += n_neighbor;
        }

        if (row_offset != total_neighbor_points) {
            throw std::runtime_error(
                "apply_owner_deferred_xnn_updates_for_candidate_box: row offset mismatch");
        }
        if (packed_total_rows == 0) {
            continue;
        }

        for (const auto& block : scratch.owned_row_blocks) {
            std::vector<DataType> block_update(
                static_cast<size_t>(block.rows * n_candidate));
            int m = static_cast<int>(block.rows);
            int n = static_cast<int>(n_candidate);
            int k = static_cast<int>(r);
            int lda = static_cast<int>(total_neighbor_points);
            int ldb = static_cast<int>(total_neighbor_points);
            int ldc = static_cast<int>(block.rows);
            DataType alpha = 1.0;
            DataType beta = 0.0;
            const DataType* source_rows_ptr =
                source_box->deferred_xnn_temp2.data() + block.source_row_offset;

            gemm_("N", "T", &m, &n, &k,
                &alpha,
                source_rows_ptr, &lda,
                source_box->X_NR.data.data() + candidate_col_offset, &ldb,
                &beta,
                block_update.data(), &ldc);

            if (block.accumulate_locally) {
                DeferredXnnTargetKey target;
                target.box_morton = candidate_morton;
                target.neighbor_morton = block.neighbor_morton;
                target.kind = block.kind;

                auto state_it = scratch.accumulated_target_indices.find(target);
                if (state_it == scratch.accumulated_target_indices.end()) {
                    size_t new_idx = scratch.accumulated_targets.size();
                    scratch.accumulated_target_indices.emplace(target, new_idx);

                    DeferredXnnAccumulatedTarget<DataType> target_state;
                    target_state.target = target;
                    target_state.rows = block.rows;
                    target_state.cols = n_candidate;
                    target_state.data =
                        materialize_deferred_xnn_target_matrix_for_accumulation(
                            target, block.rows, n_candidate, level, kernel, dimension);
                    scratch.accumulated_targets.push_back(std::move(target_state));
                    state_it = scratch.accumulated_target_indices.find(target);
                } else {
                    auto& existing = scratch.accumulated_targets[state_it->second];
                    if (existing.rows != block.rows || existing.cols != n_candidate) {
                        throw std::runtime_error(
                            "apply_owner_deferred_xnn_updates_for_candidate_box: inconsistent target dimensions");
                    }
                }

                auto& target_state = scratch.accumulated_targets[state_it->second];
                accumulate_deferred_xnn_matrix_in_place(
                    target_state.data, block_update);
            }

            if (block.emit_remote_add) {
                if (pending == nullptr) {
                    throw std::runtime_error(
                        "apply_owner_deferred_xnn_updates_for_candidate_box: remote ADD requires pending updates");
                }

                const EdgeKind edge_kind =
                    (block.kind == DeferredXnnTargetKind::SCHUR) ? EdgeKind::Diag :
                    (block.kind == DeferredXnnTargetKind::NEAR_A_NS) ? EdgeKind::Near :
                                                                       EdgeKind::Far;
                deferred_xnn_accumulate_canonical_edge_delta_from_slice(
                    *pending,
                    candidate_morton,
                    block.neighbor_morton,
                    edge_kind,
                    block_update.data(),
                    block.rows,
                    0,
                    block.rows,
                    n_candidate);
            }
        }
    }

    for (auto& target_state : scratch.accumulated_targets) {
        flush_deferred_xnn_target_matrix_from_accumulation(target_state, level);
    }
}

// Copy the owner half-pair directly into a preallocated non-owner block. By the
// time this runs, the owner-side accumulation pass is finished, and every
// mirrored destination block has already been created with the correct shape.
template<typename CoordType, typename DataType>
void apply_symmetric_owner_deferred_xnn_updates_for_candidate_box(
    BoxData<CoordType, DataType>* candidate_box,
    const std::vector<DeferredXnnTargetKey>& mirror_targets,
    TreeLevel<CoordType, DataType>& level) {
    if (candidate_box == nullptr ||
        deferred_xnn_should_skip_mirror_candidate_box(candidate_box, level)) {
        return;
    }

    for (const auto& target : mirror_targets) {
        if (target.kind == DeferredXnnTargetKind::SCHUR) {
            throw std::runtime_error(
                "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: unexpected Schur target");
        }
        if (target.box_morton != candidate_box->morton_index) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: mirror target assigned to wrong candidate"
                << " candidate_box=" << candidate_box->morton_index
                << " target_box=" << target.box_morton
                << " neighbor_box=" << target.neighbor_morton;
            throw std::runtime_error(oss.str());
        }

        DeferredXnnEndpoint<CoordType, DataType> owner_endpoint =
            resolve_deferred_xnn_endpoint(level, target.neighbor_morton);
        if (owner_endpoint.box == nullptr || owner_endpoint.is_assisting) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: owner box missing"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton;
            throw std::runtime_error(oss.str());
        }

        auto& owner_interaction_map =
            (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
                owner_endpoint.box->near_field_interaction_map :
                owner_endpoint.box->far_field_interaction_map;
        auto& owner_interactions =
            (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
                owner_endpoint.box->near_field_modified_interactions :
                owner_endpoint.box->far_field_modified_interactions;

        auto owner_it = owner_interaction_map.find(candidate_box->morton_index);
        if (owner_it == owner_interaction_map.end()) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: owner block missing"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS");
            throw std::runtime_error(oss.str());
        }

        const auto& owner_block = owner_interactions[owner_it->second];
        if (!owner_block.A_NS.is_allocated()) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: owner block storage missing"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS");
            throw std::runtime_error(oss.str());
        }

        auto& target_interaction_map =
            (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
                candidate_box->near_field_interaction_map :
                candidate_box->far_field_interaction_map;
        auto& target_interactions =
            (target.kind == DeferredXnnTargetKind::NEAR_A_NS) ?
                candidate_box->near_field_modified_interactions :
                candidate_box->far_field_modified_interactions;

        auto target_it = target_interaction_map.find(target.neighbor_morton);
        if (target_it == target_interaction_map.end()) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: preallocated mirror block missing"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS");
            throw std::runtime_error(oss.str());
        }

        auto& target_block = target_interactions[target_it->second];
        if (!target_block.A_NS.is_allocated()) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: preallocated mirror storage missing"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS");
            throw std::runtime_error(oss.str());
        }

        const int64_t expected_rows = owner_block.A_NS.cols;
        const int64_t expected_cols = owner_block.A_NS.rows;
        if (target_block.A_NS.rows != expected_rows ||
            target_block.A_NS.cols != expected_cols) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: preallocated mirror block dimension mismatch"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS")
                << " existing_rows=" << target_block.A_NS.rows
                << " existing_cols=" << target_block.A_NS.cols
                << " expected_rows=" << expected_rows
                << " expected_cols=" << expected_cols;
            throw std::runtime_error(oss.str());
        }

        const size_t expected_size = static_cast<size_t>(expected_rows * expected_cols);
        if (target_block.A_NS.data.size() != expected_size ||
            owner_block.A_NS.data.size() != expected_size) {
            std::ostringstream oss;
            oss << "apply_symmetric_owner_deferred_xnn_updates_for_candidate_box: mirror data size mismatch"
                << " candidate_box=" << candidate_box->morton_index
                << " owner_box=" << target.neighbor_morton
                << " target_kind="
                << (target.kind == DeferredXnnTargetKind::NEAR_A_NS ?
                        "NEAR_A_NS" : "FAR_A_NS")
                << " target_data_size=" << target_block.A_NS.data.size()
                << " owner_data_size=" << owner_block.A_NS.data.size()
                << " expected_size=" << expected_size;
            throw std::runtime_error(oss.str());
        }

        // Copy the transpose directly into the already allocated destination
        // buffer. No map insertions or vector growth happen in this phase.
        for (int64_t col = 0; col < owner_block.A_NS.cols; ++col) {
            for (int64_t row = 0; row < owner_block.A_NS.rows; ++row) {
                target_block.A_NS.data[static_cast<size_t>(col + row * owner_block.A_NS.cols)] =
                    owner_block.A_NS.data[static_cast<size_t>(row + col * owner_block.A_NS.rows)];
            }
        }
    }
}

// Once all deferred owner-side updates have been applied and mirrored, restore
// the source box to the usual post-step-two state by replacing X_NR with temp2.
template<typename CoordType, typename DataType>
void finalize_deferred_xnn_source_box(
    BoxData<CoordType, DataType>* source_box) {
    if (source_box == nullptr) {
        return;
    }

    if (!source_box->deferred_xnn_temp2.empty()) {
        if (!source_box->X_NR.is_allocated()) {
            throw std::runtime_error(
                "finalize_deferred_xnn_source_box: source X_NR missing");
        }
        source_box->X_NR.data = std::move(source_box->deferred_xnn_temp2);
    }

    std::vector<DataType>().swap(source_box->deferred_xnn_temp2);
    std::vector<int64_t>().swap(source_box->deferred_xnn_neighbor_point_counts);
}


template<typename CoordType, typename DataType>
inline bool is_local_box(TreeLevel<CoordType, DataType>& level, BoxData<CoordType, DataType>* box){
    return box->morton_index >= level.local_morton_start && box->morton_index <= level.local_morton_end;
}

// template<typename CoordType, typename DataType, typename KernelType>
// void gather_id_workspace(
//     BoxData<CoordType, DataType>* box,
//     TreeLevel<CoordType, DataType>& level,
//     KernelType* kernel,
//     const CoordType* unit_proxy_points,
//     int64_t num_proxy,
//     CoordType proxy_radius_factor,
//     bool is_symmetric,
//     std::vector<DataType>& workspace,
//     int64_t& workspace_rows,
//     int64_t& workspace_cols, 
//     int DEBUG = 0,
//     bool on_boundary = false) {
    
//     // determine neighbor box eliminated status
//     for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
//         int64_t neighbor_morton = box->one_hop[idx];
//         if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
//             box->use_full_set[idx] = 0;
//         }
//     }

    
//     if (box == nullptr || box->num_points == 0) {
//         workspace_rows = 0;
//         workspace_cols = 0;
//         workspace.clear();
//         return;
//     }
    
//     int dimension = (box->bounds[4] == box->bounds[5]) ? 2 : 3;
//     workspace_cols = box->num_points;
    
//     // Transform proxy points
//     CoordType proxy_radius = proxy_radius_factor * box->size;
//     std::vector<CoordType> transformed_proxy(dimension * num_proxy);
    
//     for (int64_t i = 0; i < num_proxy; ++i) {
//         for (int d = 0; d < dimension; ++d) {
//             transformed_proxy[i * dimension + d] = 
//                 box->center[d] + proxy_radius * unit_proxy_points[i * dimension + d];
//         }
//     }
//     // printf("box center, x: %f, y: %f, z: %f\n", box->center[0], box->center[1], box->center[2]);

    
//     // Step 1: Count total rows (accounting for eliminated boxes)
//     int64_t total_rows = 0;
//     std::vector<int64_t> neighbor_point_counts;

    
//     for (int64_t neighbor_morton : box->two_hop) {
//         BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
        
//         if (neighbor_box == nullptr) {
//             auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
//             if (ghost_it != level.ghost_id_to_index.end()) {
//                 neighbor_box = &level.ghost_boxes[ghost_it->second];
//             }
//         }
        
//         int64_t n_neighbor = 0;
        
//         if (neighbor_box != nullptr) {
            
//             // Check if both boxes are on boundary
//             bool both_on_boundary = on_boundary && neighbor_box->on_boundary;
            
//             // Check if neighbor is in far_field_interaction_map
//             bool in_far_field_map = (box->far_field_interaction_map.find(neighbor_morton) != 
//                                      box->far_field_interaction_map.end());
            
//             // If both on boundary AND not in far_field map, use full point set
//             if (both_on_boundary && !in_far_field_map) {
//                 n_neighbor = neighbor_box->num_points;
//             } else if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
//                 // Otherwise check if neighbor is eliminated
//                 n_neighbor = neighbor_box->skeleton_indices.size();
//             } else {
//                 n_neighbor = neighbor_box->num_points;
//             }
//             if(box->morton_index == 5){
//                 // assert(false);
//                 printf("Box %ld: nei num_points=%ld, on_boundary=%d, two hop size: %d， far_field_interaction_map size: %d, near_field_interaction_map size: %d, neighbor morton: %ld\n", box->morton_index, n_neighbor, box->on_boundary, box->two_hop.size(), box->far_field_interaction_map.size(), box->near_field_interaction_map.size(), neighbor_box->morton_index);
//             }
//         } else {
//             // Check assisting boxes
//             auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
//             if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
//                 int64_t assist_idx = assist_it->second;
//                 const auto& assist_box = level.assisting_boxes[assist_idx];
//                 n_neighbor = assist_box.coords.size() / dimension;
//                 // override size if it's an eliminated box
//                 if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
//                     n_neighbor = assist_box.skel_indices.size();
//                 }
                
//                 // For assisting boxes, check on_boundary flag
//                 bool both_on_boundary = on_boundary && assist_box.on_boundary;
//                 bool in_far_field_map = (box->far_field_interaction_map.find(neighbor_morton) != 
//                                          box->far_field_interaction_map.end());
                
//                 // If both on boundary AND not in far_field map, use full size (already set)
//                 // Otherwise, we can't slice assisting boxes anyway, so keep full size
//                 if(box->morton_index == 104){
//                     printf("assist box num: %d, num points: %d, skel points: %d\n", assist_box.morton_index, assist_box.coords.size() / dimension, assist_box.skel_indices.size());
//                 }
//             }
//         }
        
//         if (n_neighbor == 0) {
//             throw std::runtime_error(
//                 "gather_id_workspace: 2-hop neighbor box " + 
//                 std::to_string(neighbor_morton) + " not found or has no points");
//         }
        
//         neighbor_point_counts.push_back(n_neighbor);
//         total_rows += n_neighbor;
//     }
    
//     // Add neighbor transposes (if nonsymmetric)
//     if (!is_symmetric) {
//         for (int64_t n_neighbor : neighbor_point_counts) {
//             total_rows += n_neighbor;
//         }
//     }
    
//     // Add proxy points at the end
//     total_rows += num_proxy;
//     if (!is_symmetric) {
//         total_rows += num_proxy;
//     }
    
//     workspace_rows = total_rows;
//     workspace.resize(workspace_rows * workspace_cols, DataType{0.0});
    
//     // Step 2: Fill A_NS blocks (use far_field_interaction_map)
//     int64_t current_row_offset = 0;
    
    
//     for (size_t neighbor_idx = 0; neighbor_idx < box->two_hop.size(); ++neighbor_idx) {
//         int64_t neighbor_morton = box->two_hop[neighbor_idx];
//         int64_t n_neighbor = neighbor_point_counts[neighbor_idx];
        
//         auto it = box->far_field_interaction_map.find(neighbor_morton);
        
        
//         if (it != box->far_field_interaction_map.end()) {
//             // Found in modified blocks - get with slicing
//             int64_t block_idx = it->second;
//             const auto& modified_block = box->far_field_modified_interactions[block_idx];
//             if(box->morton_index == 5){
                
//                 printf("Box %ld: neighbor %ld in far_field_interaction_map, evaluating kernel directly. rows=%d, cols=%d\n", box->morton_index, neighbor_morton, modified_block.A_NS.rows, modified_block.A_NS.cols);
             
//             }
//             std::vector<DataType> sliced_block = get_sliced_neighbor_block(
//                 neighbor_morton, modified_block, level, 
//                 true,  // is_A_NS
//                 workspace_cols
//             );

            
//             // Copy into workspace
//             int64_t actual_rows = sliced_block.size() / workspace_cols;
//             for (int64_t col = 0; col < workspace_cols; ++col) {
//                 for (int64_t row = 0; row < actual_rows; ++row) {
//                     workspace[(current_row_offset + row) + col * workspace_rows] = 
//                         sliced_block[row + col * actual_rows];
//                 }
//             }
//         } else {
            
//             // Not in modified blocks - evaluate kernel
//             const CoordType* neighbor_coords = nullptr;
//             int64_t actual_n_neighbor = n_neighbor;
            
//             BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
            
//             if (neighbor_box == nullptr) {
//                 auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
//                 if (ghost_it != level.ghost_id_to_index.end()) {
//                     neighbor_box = &level.ghost_boxes[ghost_it->second];
//                 }
//             }
            
//             if (neighbor_box != nullptr) {
                
//                 // Check if both boxes are on boundary
//                 bool both_on_boundary = on_boundary && neighbor_box->on_boundary;
                
//                 // Check if eliminated
//                 bool is_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
//                                      level.eliminated_boxes.end());
                
//                 // Use skeleton ONLY if eliminated AND NOT both on boundary
//                 if (is_eliminated && !both_on_boundary) {
//                     if(box->morton_index == 5){
//                         printf("Box %ld: neighbor %ld not in far_field_interaction_map, evaluating kernel directly. on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                    
//                     }
//                     // Use skeleton points only
//                     actual_n_neighbor = neighbor_box->skeleton_indices.size();
//                     std::vector<CoordType> skeleton_coords(actual_n_neighbor * dimension);
//                     for (int64_t i = 0; i < actual_n_neighbor; ++i) {
//                         int64_t src_idx = neighbor_box->skeleton_indices[i];
//                         for (int d = 0; d < dimension; ++d) {
//                             skeleton_coords[i * dimension + d] = 
//                                 neighbor_box->point_coords[src_idx * dimension + d];
//                         }
//                     }
                    
//                     std::vector<DataType> A_NB(actual_n_neighbor * workspace_cols);
//                     kernel->evaluate_block(
//                         skeleton_coords.data(), actual_n_neighbor,
//                         box->point_coords.data(), workspace_cols,
//                         A_NB.data(), actual_n_neighbor
//                     );
                    
//                     for (int64_t col = 0; col < workspace_cols; ++col) {
//                         for (int64_t row = 0; row < actual_n_neighbor; ++row) {
//                             workspace[(current_row_offset + row) + col * workspace_rows] = 
//                                 A_NB[row + col * actual_n_neighbor];
//                         }
//                     }
//                     if(box->morton_index == 64){
//                         printf("skel Box %ld: neighbor %ld on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                        
//                     }
//                 } else {
//                     // Use all points (either not eliminated OR both on boundary)
//                     neighbor_coords = neighbor_box->point_coords.data();
                    
//                     std::vector<DataType> A_NB(n_neighbor * workspace_cols);
//                     kernel->evaluate_block(
//                         neighbor_coords, n_neighbor,
//                         box->point_coords.data(), workspace_cols,
//                         A_NB.data(), n_neighbor
//                     );
                    
//                     for (int64_t col = 0; col < workspace_cols; ++col) {
//                         for (int64_t row = 0; row < n_neighbor; ++row) {
//                             workspace[(current_row_offset + row) + col * workspace_rows] = 
//                                 A_NB[row + col * n_neighbor];
//                         }
//                     }
//                     if(box->morton_index == 64){
//                         printf("Box %ld: neighbor %ld on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                        
//                     }
//                 }
//             } else {
                
//                 // assisting neighbor: use full coords if not eliminated; use skeleton coords if eliminated
//                 auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
//                 if (assist_it == level.assisting_box_points_for_kernel_evaluation.end()) {
//                     throw std::runtime_error("gather_id_workspace: missing assisting data for neighbor " +
//                                             std::to_string(neighbor_morton));
//                 }

//                 const int64_t assist_idx = assist_it->second;
//                 if (assist_idx < 0 || (size_t)assist_idx >= level.assisting_boxes.size()) {
//                     throw std::runtime_error("gather_id_workspace: assisting index out of range for neighbor " +
//                                             std::to_string(neighbor_morton));
//                 }

//                 auto& req = level.assisting_boxes[(size_t)assist_idx];


//                 // Is the neighbor eliminated? If yes, use its skeleton points only.
//                 const bool neighbor_eliminated =
//                     (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

//                 const CoordType* neighbor_coords_ptr = nullptr;
//                 std::vector<CoordType> neighbor_skel_coords; // storage if we need to gather skeleton coords

//                 if (neighbor_eliminated) {
//                     if (req.skel_indices.empty()) {
//                         throw std::runtime_error("gather_id_workspace: neighbor " + std::to_string(neighbor_morton) +
//                                                 " is eliminated but assisting skel_indices is empty");
//                     }
//                     if ((int64_t)req.skel_indices.size() < n_neighbor) {
//                         throw std::runtime_error("gather_id_workspace: assisting skel_indices too small for neighbor " +
//                                                 std::to_string(neighbor_morton) + " (need n_neighbor=" +
//                                                 std::to_string(n_neighbor) + ", have " +
//                                                 std::to_string(req.skel_indices.size()) + ")");
//                     }
//                     if ((int64_t)req.coords.size() < dimension) {
//                         throw std::runtime_error("gather_id_workspace: assisting coords empty for neighbor " +
//                                                 std::to_string(neighbor_morton));
//                     }

//                     neighbor_skel_coords.resize((size_t)n_neighbor * (size_t)dimension);
//                     for (int64_t i = 0; i < n_neighbor; ++i) {
//                         const int64_t src = req.skel_indices[(size_t)i];
//                         if (src < 0 || (size_t)(src * dimension + (dimension - 1)) >= req.coords.size()) {
//                             throw std::runtime_error("gather_id_workspace: bad skel index " + std::to_string(src) +
//                                                     " for neighbor " + std::to_string(neighbor_morton));
//                         }
//                         for (int d = 0; d < dimension; ++d) {
//                             neighbor_skel_coords[(size_t)i * dimension + (size_t)d] =
//                                 req.coords[(size_t)src * dimension + (size_t)d];
//                         }
//                     }
//                     neighbor_coords_ptr = neighbor_skel_coords.data();
//                     if(box->morton_index == 64){
//                         printf("eliminate assisting Box %ld: neighbor %ld on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                        
//                     }
//                 } else {
//                     if(box->morton_index == 64){
//                         printf("non-eliminate assisting Box %ld: neighbor %ld on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                        
//                     }
//                     // full (uneliminated) neighbor coordinates
//                     // (Optional sanity check: req.coords.size() == n_neighbor*dim)
//                     neighbor_coords_ptr = req.coords.data();
//                 }

//                 // Evaluate block A(neighbor, box_points)
//                 std::vector<DataType> A_NB((size_t)n_neighbor * (size_t)workspace_cols);
//                 kernel->evaluate_block(
//                     neighbor_coords_ptr, n_neighbor,
//                     box->point_coords.data(), workspace_cols,
//                     A_NB.data(), n_neighbor
//                 );

//                 // Scatter into workspace (column-major)
//                 for (int64_t col = 0; col < workspace_cols; ++col) {
//                     for (int64_t row = 0; row < n_neighbor; ++row) {
//                         workspace[(current_row_offset + row) + col * workspace_rows] =
//                             A_NB[(size_t)row + (size_t)col * (size_t)n_neighbor];
//                     }
//                 }
//             }
//         }
        
//         current_row_offset += n_neighbor;
//     }
    
//     // Step 3: Fill A_SN transpose blocks (nonsymmetric only, use far_field_interaction_map_nonsymmetry)
//     if (!is_symmetric) {
//         for (size_t neighbor_idx = 0; neighbor_idx < box->two_hop.size(); ++neighbor_idx) {
//             int64_t neighbor_morton = box->two_hop[neighbor_idx];
//             int64_t n_neighbor = neighbor_point_counts[neighbor_idx];
            
//             auto it = box->far_field_interaction_map_nonsymmetry.find(neighbor_morton);
            
//             if (it != box->far_field_interaction_map_nonsymmetry.end()) {
//                 // Found in modified blocks - get with slicing
//                 int64_t block_idx = it->second;
//                 const auto& modified_block = box->far_field_modified_interactions[block_idx];
                
//                 if (!modified_block.A_SN.is_allocated()) {
//                     throw std::runtime_error(
//                         "gather_id_workspace: Non-symmetric problem but A_SN not allocated for neighbor " + 
//                         std::to_string(neighbor_morton));
//                 }
                
//                 std::vector<DataType> sliced_block = get_sliced_neighbor_block(
//                     neighbor_morton, modified_block, level,
//                     false,  // is_A_NS = false (this is A_SN)
//                     workspace_cols
//                 );
                
//                 // A_SN is (workspace_cols × n_neighbor), need to transpose for stacking
//                 int64_t actual_cols = sliced_block.size() / workspace_cols;
//                 for (int64_t col = 0; col < workspace_cols; ++col) {
//                     for (int64_t row = 0; row < actual_cols; ++row) {
//                         workspace[(current_row_offset + row) + col * workspace_rows] = 
//                             sliced_block[col + row * workspace_cols];  // Transpose
//                     }
//                 }
//             } else {
//                 // Not in modified blocks - evaluate kernel
//                 const CoordType* neighbor_coords = nullptr;
//                 int64_t actual_n_neighbor = n_neighbor;
                
//                 BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
                
//                 if (neighbor_box == nullptr) {
//                     auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
//                     if (ghost_it != level.ghost_id_to_index.end()) {
//                         neighbor_box = &level.ghost_boxes[ghost_it->second];
//                     }
//                 }
                
//                 if (neighbor_box != nullptr) {
//                     // Check if both boxes are on boundary
//                     bool both_on_boundary = on_boundary && neighbor_box->on_boundary;
                    
//                     // Check if eliminated
//                     bool is_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
//                                          level.eliminated_boxes.end());
                    
//                     // Use skeleton ONLY if eliminated AND NOT both on boundary
//                     if (is_eliminated && !both_on_boundary) {
//                         // Use skeleton points only
//                         actual_n_neighbor = neighbor_box->skeleton_indices.size();
//                         std::vector<CoordType> skeleton_coords(actual_n_neighbor * dimension);
//                         for (int64_t i = 0; i < actual_n_neighbor; ++i) {
//                             int64_t src_idx = neighbor_box->skeleton_indices[i];
//                             for (int d = 0; d < dimension; ++d) {
//                                 skeleton_coords[i * dimension + d] = 
//                                     neighbor_box->point_coords[src_idx * dimension + d];
//                             }
//                         }
                        
//                         std::vector<DataType> A_BN(workspace_cols * actual_n_neighbor);
//                         kernel->evaluate_block(
//                             box->point_coords.data(), workspace_cols,
//                             skeleton_coords.data(), actual_n_neighbor,
//                             A_BN.data(), workspace_cols
//                         );
                        
//                         for (int64_t col = 0; col < workspace_cols; ++col) {
//                             for (int64_t row = 0; row < actual_n_neighbor; ++row) {
//                                 workspace[(current_row_offset + row) + col * workspace_rows] = 
//                                     A_BN[col + row * workspace_cols];
//                             }
//                         }
//                     } else {
//                         // Use all points (either not eliminated OR both on boundary)
//                         neighbor_coords = neighbor_box->point_coords.data();
                        
//                         std::vector<DataType> A_BN(workspace_cols * n_neighbor);
//                         kernel->evaluate_block(
//                             box->point_coords.data(), workspace_cols,
//                             neighbor_coords, n_neighbor,
//                             A_BN.data(), workspace_cols
//                         );
                        
//                         for (int64_t col = 0; col < workspace_cols; ++col) {
//                             for (int64_t row = 0; row < n_neighbor; ++row) {
//                                 workspace[(current_row_offset + row) + col * workspace_rows] = 
//                                     A_BN[col + row * workspace_cols];
//                             }
//                         }
//                     }
//                 } else {
//                     // assisting neighbor: use full coords if not eliminated; use skeleton coords if eliminated
//                     auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
//                     if (assist_it == level.assisting_box_points_for_kernel_evaluation.end()) {
//                         throw std::runtime_error("gather_id_workspace: missing assisting data for neighbor " +
//                                                 std::to_string(neighbor_morton));
//                     }

//                     const int64_t assist_idx = assist_it->second;
//                     if (assist_idx < 0 || (size_t)assist_idx >= level.assisting_boxes.size()) {
//                         throw std::runtime_error("gather_id_workspace: assisting index out of range for neighbor " +
//                                                 std::to_string(neighbor_morton));
//                     }

//                     auto& req = level.assisting_boxes[(size_t)assist_idx];


//                     // Is the neighbor eliminated? If yes, use its skeleton points only.
//                     const bool neighbor_eliminated =
//                         (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

//                     const CoordType* neighbor_coords_ptr = nullptr;
//                     std::vector<CoordType> neighbor_skel_coords; // storage if we need to gather skeleton coords

//                     if (neighbor_eliminated) {
//                         if (req.skel_indices.empty()) {
//                             throw std::runtime_error("gather_id_workspace: neighbor " + std::to_string(neighbor_morton) +
//                                                     " is eliminated but assisting skel_indices is empty");
//                         }
//                         if ((int64_t)req.skel_indices.size() < n_neighbor) {
//                             throw std::runtime_error("gather_id_workspace: assisting skel_indices too small for neighbor " +
//                                                     std::to_string(neighbor_morton) + " (need n_neighbor=" +
//                                                     std::to_string(n_neighbor) + ", have " +
//                                                     std::to_string(req.skel_indices.size()) + ")");
//                         }
//                         if ((int64_t)req.coords.size() < dimension) {
//                             throw std::runtime_error("gather_id_workspace: assisting coords empty for neighbor " +
//                                                     std::to_string(neighbor_morton));
//                         }

//                         neighbor_skel_coords.resize((size_t)n_neighbor * (size_t)dimension);
//                         for (int64_t i = 0; i < n_neighbor; ++i) {
//                             const int64_t src = req.skel_indices[(size_t)i];
//                             if (src < 0 || (size_t)(src * dimension + (dimension - 1)) >= req.coords.size()) {
//                                 throw std::runtime_error("gather_id_workspace: bad skel index " + std::to_string(src) +
//                                                         " for neighbor " + std::to_string(neighbor_morton));
//                             }
//                             for (int d = 0; d < dimension; ++d) {
//                                 neighbor_skel_coords[(size_t)i * dimension + (size_t)d] =
//                                     req.coords[(size_t)src * dimension + (size_t)d];
//                             }
//                         }
//                         neighbor_coords_ptr = neighbor_skel_coords.data();
//                     } else {
//                         // full (uneliminated) neighbor coordinates
//                         // (Optional sanity check: req.coords.size() == n_neighbor*dim)
//                         neighbor_coords_ptr = req.coords.data();
//                     }

//                     // Evaluate block A(neighbor, box_points)
//                     std::vector<DataType> A_NB((size_t)n_neighbor * (size_t)workspace_cols);
//                     kernel->evaluate_block(
//                         neighbor_coords_ptr, n_neighbor,
//                         box->point_coords.data(), workspace_cols,
//                         A_NB.data(), n_neighbor
//                     );

//                     // Scatter into workspace (column-major)
//                     for (int64_t col = 0; col < workspace_cols; ++col) {
//                         for (int64_t row = 0; row < n_neighbor; ++row) {
//                             workspace[(current_row_offset + row) + col * workspace_rows] =
//                                 A_NB[(size_t)row + (size_t)col * (size_t)n_neighbor];
//                         }
//                     }
//                 }
//             }
            
//             current_row_offset += n_neighbor;
//         }
//     }
    
//     // Step 4: Add proxy blocks at the end (no changes needed - not 2-hop neighbor interaction)
//     std::vector<DataType> A_proxy_B(num_proxy * workspace_cols);
//     kernel->evaluate_block(
//         transformed_proxy.data(), num_proxy,
//         box->point_coords.data(), workspace_cols,
//         A_proxy_B.data(), num_proxy
//     );
    
//     for (int64_t col = 0; col < workspace_cols; ++col) {
//         for (int64_t row = 0; row < num_proxy; ++row) {
//             workspace[(current_row_offset + row) + col * workspace_rows] = 
//                 A_proxy_B[row + col * num_proxy];
//         }
//     }
    
//     current_row_offset += num_proxy;
    
//     if (!is_symmetric) {
//         std::vector<DataType> A_B_proxy(workspace_cols * num_proxy);
//         kernel->evaluate_block(
//             box->point_coords.data(), workspace_cols,
//             transformed_proxy.data(), num_proxy,
//             A_B_proxy.data(), workspace_cols
//         );
        
//         for (int64_t col = 0; col < workspace_cols; ++col) {
//             for (int64_t row = 0; row < num_proxy; ++row) {
//                 workspace[(current_row_offset + row) + col * workspace_rows] = 
//                     A_B_proxy[col + row * workspace_cols];
//             }
//         }
        
//         current_row_offset += num_proxy;
//     }
    
//     if (current_row_offset != workspace_rows) {
//         throw std::runtime_error(
//             "gather_id_workspace: Row count mismatch. Expected " + 
//             std::to_string(workspace_rows) + ", got " + 
//             std::to_string(current_row_offset));
//     }

//     if (DEBUG){
//         print_workspace_segment_stats(workspace, workspace_rows, workspace_cols,
//                 neighbor_point_counts,
//                 nullptr, 0, std::cerr);
//         for (int64_t neighbor_morton : box->two_hop) {
//             printf("segment neighbor morton: %ld\n", neighbor_morton);
//         }
//     }
   
    
//     // // DEBUG: Show if any near-field blocks were created
//     // if (DEBUG) {
//     //     std::cout << "DEBUG gather_id_workspace END: Box " << box->morton_index << std::endl;
//     //     std::cout << "  near_field_modified_interactions: " << box->near_field_modified_interactions.size() << std::endl;
//     //     for (const auto& [morton, idx] : box->near_field_interaction_map) {
//     //         std::cout << "    Neighbor " << morton << std::endl;
//     //     }
//     // }
// }

struct VecStats {
    double min =  std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    double max_abs = 0.0;
    long double sum = 0.0L;
    long double sum_sq = 0.0L;
    int64_t imin = -1, imax = -1, imaxabs = -1;
    bool has_nan = false;
    bool has_inf = false;
    int64_t count = 0;
};

static VecStats stats_workspace(const std::vector<double>& v) {
    VecStats s;
    s.count = (int64_t)v.size();
    for (int64_t i = 0; i < (int64_t)v.size(); ++i) {
        const double x = v[(size_t)i];
        if (std::isnan(x)) { s.has_nan = true; continue; }
        if (!std::isfinite(x)) { s.has_inf = true; continue; }

        if (x < s.min) { s.min = x; s.imin = i; }
        if (x > s.max) { s.max = x; s.imax = i; }
        const double ax = std::abs(x);
        if (ax > s.max_abs) { s.max_abs = ax; s.imaxabs = i; }

        s.sum += (long double)x;
        s.sum_sq += (long double)x * (long double)x;
    }
    return s;
}

static void print_workspace_stats(const std::vector<double>& v,
                                 const std::string& name = "workspace",
                                 std::ostream& os = std::cerr) {
    auto s = stats_workspace(v);

    os << name << " stats:\n";
    os << "  size=" << s.count
       << " has_nan=" << (s.has_nan ? "YES" : "NO")
       << " has_inf=" << (s.has_inf ? "YES" : "NO") << "\n";

    if (s.count == 0) return;

    if (std::isfinite(s.min) && std::isfinite(s.max)) {
        const long double mean = s.sum / (long double)s.count;
        const long double var = std::max((long double)0.0, s.sum_sq / (long double)s.count - mean * mean);
        const long double rms = std::sqrt(s.sum_sq / (long double)s.count);

        os << "  min=" << s.min << " @ " << s.imin << "\n";
        os << "  max=" << s.max << " @ " << s.imax << "\n";
        os << "  max_abs=" << s.max_abs << " @ " << s.imaxabs << "\n";
        os << "  mean=" << (double)mean << "\n";
        os << "  rms=" << (double)rms << "\n";
        os << "  std=" << (double)std::sqrt(var) << "\n";
    } else {
        os << "  min/max not finite (likely all NaN/Inf)\n";
    }
}


template<typename CoordType, typename DataType, typename KernelType>
void gather_id_workspace(
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    const CoordType* unit_proxy_points,
    int64_t num_proxy,
    CoordType proxy_radius_factor,
    bool is_symmetric,
    std::vector<DataType>& workspace,
    int64_t& workspace_rows,
    int64_t& workspace_cols, 
    int DEBUG = 0,
    bool on_boundary = false) {

    if (box == nullptr || box->num_points == 0) {
        workspace_rows = 0;
        workspace_cols = 0;
        workspace.clear();
        return;
    }
    
    // determine neighbor box eliminated status
    for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
        int64_t neighbor_morton = box->one_hop[idx];
        if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
            box->use_full_set[idx] = 0;
        }
    }

    
   
    
    int dimension = level.dimension;
    workspace_cols = box->num_points;
    
    // Transform proxy points
    std::vector<CoordType> transformed_proxy;
    if (num_proxy > 0) {
        CoordType proxy_radius = proxy_radius_factor * box->size;
        transformed_proxy.resize(dimension * num_proxy);

        for (int64_t i = 0; i < num_proxy; ++i) {
            for (int d = 0; d < dimension; ++d) {
                transformed_proxy[i * dimension + d] =
                    box->center[d] + proxy_radius * unit_proxy_points[i * dimension + d];
            }
        }
    }
    // printf("box center, x: %f, y: %f, z: %f\n", box->center[0], box->center[1], box->center[2]);

    
    // Step 1: Count total rows (accounting for eliminated boxes)
    int64_t total_rows = 0;
    std::vector<int64_t> neighbor_point_counts;

    
    for (int64_t neighbor_morton : box->two_hop) {
        BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
        
        if (neighbor_box == nullptr) {
            auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
            if (ghost_it != level.ghost_id_to_index.end()) {
                neighbor_box = &level.ghost_boxes[ghost_it->second];
            }
        }
        
        int64_t n_neighbor = 0;
        
        if (neighbor_box != nullptr) {
            

            // // Check if neighbor is in far_field_interaction_map
            // bool in_far_field_map = (box->far_field_interaction_map.find(neighbor_morton) != 
            //                          box->far_field_interaction_map.end());
            
            // If either on boundary AND not in far_field map, use full point set
            if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                // Otherwise check if neighbor is eliminated
                n_neighbor = neighbor_box->skeleton_indices.size();
            } else {
                n_neighbor = neighbor_box->num_points;
            }
            // gggggggggggggggggggggg
            // if(!in_far_field_map){
            //     n_neighbor = neighbor_box->num_points;
            // }

        } else {
            // Check assisting boxes
            auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
            if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
                int64_t assist_idx = assist_it->second;
                const auto& assist_box = level.assisting_boxes[assist_idx];
                n_neighbor = assist_box.coords.size() / dimension;
                // override size if it's an eliminated box
                if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                    n_neighbor = assist_box.skel_indices.size();
                }
                // gggggggggggggggggggggg
                // n_neighbor = assist_box.coords.size() / dimension;
            }else{
                throw std::runtime_error("assisting box must be found\n");
            }
        }
        
        if (n_neighbor == 0) {
            throw std::runtime_error(
                "gather_id_workspace: 2-hop neighbor box " + 
                std::to_string(neighbor_morton) + " not found or has no points");
        }
        
        neighbor_point_counts.push_back(n_neighbor);
        total_rows += n_neighbor;
    }
    
    // Add neighbor transposes (if nonsymmetric)
    if (!is_symmetric) {
        for (int64_t n_neighbor : neighbor_point_counts) {
            total_rows += n_neighbor;
        }
    }
    
    // Add proxy points at the end
    total_rows += num_proxy;
    if (!is_symmetric) {
        total_rows += num_proxy;
    }
    
    workspace_rows = total_rows;
    workspace.resize(workspace_rows * workspace_cols, DataType{0.0});
    
    // Step 2: Fill A_NS blocks (use far_field_interaction_map)
    int64_t current_row_offset = 0;
    std::vector<DataType> sliced_block;
    
    
    for (size_t neighbor_idx = 0; neighbor_idx < box->two_hop.size(); ++neighbor_idx) {
        int64_t neighbor_morton = box->two_hop[neighbor_idx];
        int64_t n_neighbor = neighbor_point_counts[neighbor_idx];
        // if(box->far_field_interaction_map.size() != 0){
        //     printf("Box %ld: far_field_interaction_map size: %ld\n", box->morton_index, box->far_field_interaction_map.size());
        //     fflush(stdout);
        //     assert(box->far_field_interaction_map.size() == 0);
        // }
        
        auto it = box->far_field_interaction_map.find(neighbor_morton);
        
        
        if (it != box->far_field_interaction_map.end()) {
            // Found in modified blocks - get with slicing
            int64_t block_idx = it->second;
            const auto& modified_block = box->far_field_modified_interactions[block_idx];
            
            get_sliced_neighbor_block_into(
                neighbor_morton, modified_block, level, 
                true,  // is_A_NS
                workspace_cols,
                sliced_block
            );
            // ggggggggggggggggggg
            // if(sliced_block.size() != modified_block.A_NS.rows * modified_block.A_NS.cols){
            //     printf("box: %ld, neighbor: %ld, sliced_block size: %ld, expected size: %ld\n", box->morton_index, neighbor_morton, sliced_block.size(), modified_block.A_NS.rows * modified_block.A_NS.cols);
            //     fflush(stdout);
            //     assert(sliced_block.size() == modified_block.A_NS.rows * modified_block.A_NS.cols);
            // }
            
            
            // Copy into workspace
            int64_t actual_rows = sliced_block.size() / workspace_cols;
            assert(actual_rows == n_neighbor);  
            for (int64_t col = 0; col < workspace_cols; ++col) {
                for (int64_t row = 0; row < actual_rows; ++row) {
                    workspace[(current_row_offset + row) + col * workspace_rows] = 
                        sliced_block[row + col * actual_rows];
                }
            }
        } else {
            
            // Not in modified blocks - evaluate kernel
            const CoordType* neighbor_coords = nullptr;
            int64_t actual_n_neighbor = n_neighbor;
            
            BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
            
            if (neighbor_box == nullptr) {
                auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
                if (ghost_it != level.ghost_id_to_index.end()) {
                    neighbor_box = &level.ghost_boxes[ghost_it->second];
                }
            }
            
            if (neighbor_box != nullptr) {
                
                
                // Check if eliminated
                bool is_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
                                     level.eliminated_boxes.end());
                
                // Use skeleton ONLY if eliminated AND NOT both on boundary
                // gggggggggggggggggggggg
                
                // if (false) {
                if (is_eliminated) {
                    
                    // Use skeleton points only
                    actual_n_neighbor = neighbor_box->skeleton_indices.size();
                    assert(actual_n_neighbor == n_neighbor);
                    std::vector<CoordType> skeleton_coords(actual_n_neighbor * dimension);
                    for (int64_t i = 0; i < actual_n_neighbor; ++i) {
                        int64_t src_idx = neighbor_box->skeleton_indices[i];
                        for (int d = 0; d < dimension; ++d) {
                            skeleton_coords[i * dimension + d] = 
                                neighbor_box->point_coords[src_idx * dimension + d];
                        }
                    }
                    
                    std::vector<DataType> A_NB(actual_n_neighbor * workspace_cols);
                    kernel->evaluate_block(
                        skeleton_coords.data(), actual_n_neighbor,
                        box->point_coords.data(), workspace_cols,
                        A_NB.data(), actual_n_neighbor
                    );
                    
                    for (int64_t col = 0; col < workspace_cols; ++col) {
                        for (int64_t row = 0; row < actual_n_neighbor; ++row) {
                            workspace[(current_row_offset + row) + col * workspace_rows] = 
                                A_NB[row + col * actual_n_neighbor];
                        }
                    }
                    
                } else {
                    // Use all points (either not eliminated OR both on boundary)
                    neighbor_coords = neighbor_box->point_coords.data();
                    assert(neighbor_box->point_coords.size() == (size_t)n_neighbor * (size_t)dimension);
                    
                    std::vector<DataType> A_NB(n_neighbor * workspace_cols);
                    kernel->evaluate_block(
                        neighbor_coords, n_neighbor,
                        box->point_coords.data(), workspace_cols,
                        A_NB.data(), n_neighbor
                    );
                    
                    for (int64_t col = 0; col < workspace_cols; ++col) {
                        for (int64_t row = 0; row < n_neighbor; ++row) {
                            workspace[(current_row_offset + row) + col * workspace_rows] = 
                                A_NB[row + col * n_neighbor];
                        }
                    }
                    // if(box->morton_index == 64){
                    //     printf("Box %ld: neighbor %ld on_boundary=%d\n", box->morton_index, neighbor_morton, on_boundary);
                        
                    // }
                }
            } else {
                
                // assisting neighbor: use full coords if not eliminated; use skeleton coords if eliminated
                auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
                if (assist_it == level.assisting_box_points_for_kernel_evaluation.end()) {
                    throw std::runtime_error("gather_id_workspace: missing assisting data for neighbor " +
                                            std::to_string(neighbor_morton));
                }

                const int64_t assist_idx = assist_it->second;
                if (assist_idx < 0 || (size_t)assist_idx >= level.assisting_boxes.size()) {
                    throw std::runtime_error("gather_id_workspace: assisting index out of range for neighbor " +
                                            std::to_string(neighbor_morton));
                }

                auto& req = level.assisting_boxes[(size_t)assist_idx];


                // Is the neighbor eliminated? If yes, use its skeleton points only.
                const bool neighbor_eliminated =
                    (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

                const CoordType* neighbor_coords_ptr = nullptr;
                std::vector<CoordType> neighbor_skel_coords; // storage if we need to gather skeleton coords
                // ggggggggggggggggggggg
                // if (false) {
                if (neighbor_eliminated) {
                    if (req.skel_indices.empty()) {
                        throw std::runtime_error("gather_id_workspace: neighbor " + std::to_string(neighbor_morton) +
                                                " is eliminated but assisting skel_indices is empty");
                    }
                    if ((int64_t)req.skel_indices.size() < n_neighbor) {
                        throw std::runtime_error("gather_id_workspace: assisting skel_indices too small for neighbor " +
                                                std::to_string(neighbor_morton) + " (need n_neighbor=" +
                                                std::to_string(n_neighbor) + ", have " +
                                                std::to_string(req.skel_indices.size()) + ")");
                    }
                    if ((int64_t)req.coords.size() < dimension) {
                        throw std::runtime_error("gather_id_workspace: assisting coords empty for neighbor " +
                                                std::to_string(neighbor_morton));
                    }

                    neighbor_skel_coords.resize((size_t)n_neighbor * (size_t)dimension);
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        const int64_t src = req.skel_indices[(size_t)i];
                        if (src < 0 || (size_t)(src * dimension + (dimension - 1)) >= req.coords.size()) {
                            throw std::runtime_error("gather_id_workspace: bad skel index " + std::to_string(src) +
                                                    " for neighbor " + std::to_string(neighbor_morton));
                        }
                        for (int d = 0; d < dimension; ++d) {
                            neighbor_skel_coords[(size_t)i * dimension + (size_t)d] =
                                req.coords[(size_t)src * dimension + (size_t)d];
                        }
                    }
                    neighbor_coords_ptr = neighbor_skel_coords.data();
                    assert(neighbor_skel_coords.size() == (size_t)n_neighbor * (size_t)dimension);

                } else {

                    // full (uneliminated) neighbor coordinates
                    // (Optional sanity check: req.coords.size() == n_neighbor*dim)
                    neighbor_coords_ptr = req.coords.data();
                    assert(req.coords.size() == (size_t)n_neighbor * (size_t)dimension);
                }

                // Evaluate block A(neighbor, box_points)
                std::vector<DataType> A_NB((size_t)n_neighbor * (size_t)workspace_cols);
                kernel->evaluate_block(
                    neighbor_coords_ptr, n_neighbor,
                    box->point_coords.data(), workspace_cols,
                    A_NB.data(), n_neighbor
                );

                // Scatter into workspace (column-major)
                for (int64_t col = 0; col < workspace_cols; ++col) {
                    for (int64_t row = 0; row < n_neighbor; ++row) {
                        workspace[(current_row_offset + row) + col * workspace_rows] =
                            A_NB[(size_t)row + (size_t)col * (size_t)n_neighbor];
                    }
                }
            }
        }
        
        current_row_offset += n_neighbor;
    }
    
    
    
    // Step 4: Add proxy blocks at the end (no changes needed - not 2-hop neighbor interaction)
    if (num_proxy > 0) {
        std::vector<DataType> A_proxy_B(num_proxy * workspace_cols);
        kernel->evaluate_block(
            transformed_proxy.data(), num_proxy,
            box->point_coords.data(), workspace_cols,
            A_proxy_B.data(), num_proxy
        );

        for (int64_t col = 0; col < workspace_cols; ++col) {
            for (int64_t row = 0; row < num_proxy; ++row) {
                workspace[(current_row_offset + row) + col * workspace_rows] =
                    A_proxy_B[row + col * num_proxy];
            }
        }

        current_row_offset += num_proxy;

        if (!is_symmetric) {
            std::vector<DataType> A_B_proxy(workspace_cols * num_proxy);
            kernel->evaluate_block(
                box->point_coords.data(), workspace_cols,
                transformed_proxy.data(), num_proxy,
                A_B_proxy.data(), workspace_cols
            );

            for (int64_t col = 0; col < workspace_cols; ++col) {
                for (int64_t row = 0; row < num_proxy; ++row) {
                    workspace[(current_row_offset + row) + col * workspace_rows] =
                        A_B_proxy[col + row * workspace_cols];
                }
            }

            current_row_offset += num_proxy;
        }
    }
    
    if (current_row_offset != workspace_rows) {
        throw std::runtime_error(
            "gather_id_workspace: Row count mismatch. Expected " + 
            std::to_string(workspace_rows) + ", got " + 
            std::to_string(current_row_offset));
    }

    // if (DEBUG){
    //     print_workspace_segment_stats(workspace, workspace_rows, workspace_cols,
    //             neighbor_point_counts,
    //             nullptr, 0, std::cerr);
    //     for (int64_t neighbor_morton : box->two_hop) {
    //         printf("segment neighbor morton: %ld\n", neighbor_morton);
    //     }
    // }
   
    
    // // DEBUG: Show if any near-field blocks were created
    // if (DEBUG) {
    //     std::cout << "DEBUG gather_id_workspace END: Box " << box->morton_index << std::endl;
    //     std::cout << "  near_field_modified_interactions: " << box->near_field_modified_interactions.size() << std::endl;
    //     for (const auto& [morton, idx] : box->near_field_interaction_map) {
    //         std::cout << "    Neighbor " << morton << std::endl;
    //     }
    // }
}



// Returns pointer to coordinates for `morton` with `n_use` points (n_use must match caller's n_*).
// If the box is eliminated, returns sliced skeleton coords (stored in `tmp`).
template <typename CoordType, typename DataType>
static const CoordType* coords_ptr_maybe_sliced(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton,
    int64_t n_use,
    int dim,
    const BoxData<CoordType, DataType>* box_or_null,                 // non-null if not assisting
    const PointDataRequest<CoordType>* assist_or_null,               // non-null if assisting
    std::vector<CoordType>& tmp)                                     // storage for sliced coords
{
    const bool eliminated = (level.eliminated_boxes.find(morton) != level.eliminated_boxes.end());

    if (!eliminated) {
        if (assist_or_null) {
            if (assist_or_null->coords.empty())
                throw std::runtime_error("coords_ptr_maybe_sliced: assisting coords empty for morton=" +
                                         std::to_string(morton));
            return assist_or_null->coords.data(); // full coords
        } else {
            if (!box_or_null) throw std::runtime_error("coords_ptr_maybe_sliced: box ptr null for morton=" +
                                                       std::to_string(morton));
            if (box_or_null->point_coords.empty())
                throw std::runtime_error("coords_ptr_maybe_sliced: local point_coords empty for morton=" +
                                         std::to_string(morton));
            return box_or_null->point_coords.data(); // full coords
        }
    }

    // eliminated => use skeleton coords
    const std::vector<int64_t>* skel_idx = nullptr;
    const std::vector<CoordType>* coords = nullptr;

    if (assist_or_null) {
        skel_idx = &assist_or_null->skel_indices;
        coords   = &assist_or_null->coords;
    } else {
        if (!box_or_null) throw std::runtime_error("coords_ptr_maybe_sliced: box ptr null for morton=" +
                                                   std::to_string(morton));
        skel_idx = &box_or_null->skeleton_indices;
        coords   = &box_or_null->point_coords;
    }

    if ((int64_t)skel_idx->size() < n_use) {
        throw std::runtime_error("coords_ptr_maybe_sliced: need " + std::to_string(n_use) +
                                 " skeleton indices but have " + std::to_string(skel_idx->size()) +
                                 " for morton=" + std::to_string(morton));
    }
    if ((int64_t)coords->size() < dim) {
        throw std::runtime_error("coords_ptr_maybe_sliced: coords empty for morton=" +
                                 std::to_string(morton));
    }

    tmp.resize((size_t)n_use * (size_t)dim);
    for (int64_t i = 0; i < n_use; ++i) {
        const int64_t src = (*skel_idx)[(size_t)i];
        const int64_t base = src * dim;
        if (src < 0 || (size_t)(base + (dim - 1)) >= coords->size()) {
            throw std::runtime_error("coords_ptr_maybe_sliced: bad skeleton index " +
                                     std::to_string(src) + " for morton=" + std::to_string(morton));
        }
        for (int d = 0; d < dim; ++d) {
            tmp[(size_t)i * (size_t)dim + (size_t)d] =
                (*coords)[(size_t)base + (size_t)d];
        }
    }
    return tmp.data();
}


template<typename CoordType, typename DataType>
void slice_far_field_blocks(
    TreeLevel<CoordType, DataType>& level,
    bool is_symmetric,
    bool is_hermitian = false) {

    
    // (a) Slice current box B's far-field blocks (S dimension: B's skeleton)
    // (b) Slice neighbors' reciprocal far-field blocks (N dimension: B's skeleton)

    for (auto& box : level.local_boxes) {
        int64_t k = box.skeleton_indices.size();
        if (k == 0) continue;

        // ====================================================================
        // Part (a): Slice this box's own far-field blocks
        // A_NS: (n_neighbor × num_points) → (n_neighbor × k)
        // A_SN: (num_points × n_neighbor) → (k × n_neighbor)
        // ====================================================================
        for (auto& block : box.far_field_modified_interactions) {

            if (block.A_NS.is_allocated() && block.A_NS.cols != k) {
                int64_t stored_rows = block.A_NS.rows;
                std::vector<int64_t> all_rows(stored_rows);
                std::iota(all_rows.begin(), all_rows.end(), 0);

                std::vector<DataType> sliced(stored_rows * k);
                extract_submatrix(
                    block.A_NS.data.data(), stored_rows,
                    all_rows, box.skeleton_indices,
                    sliced.data(), stored_rows
                );

                block.A_NS.data = std::move(sliced);
                block.A_NS.cols = k;
            }
        }

        // ====================================================================
        // Part (b): Slice the reciprocal view in each neighbor's far-field
        // i.e. neighbor G has box B in its far-field map; slice the N dimension
        // (which corresponds to B's skeleton, size k)
        // ====================================================================
        for (const auto& block : box.far_field_modified_interactions) {
            int64_t neighbor_morton = block.neighbor_morton;

            BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
            if (neighbor_box == nullptr) {
                neighbor_box = level.find_ghost_box(neighbor_morton);
            }
            if (neighbor_box == nullptr) continue;  // Not on this process

            auto& neighbor_far_map = (is_symmetric || is_hermitian)
                ? neighbor_box->far_field_interaction_map
                : neighbor_box->far_field_interaction_map_nonsymmetry;

            auto it = neighbor_far_map.find(box.morton_index);
            if (it == neighbor_far_map.end()) continue;

            auto& neighbor_block = neighbor_box->far_field_modified_interactions[it->second];

            // Slice A_NS rows: (n_B × k_G) → (k_B × k_G)
            if (neighbor_block.A_NS.is_allocated() && neighbor_block.A_NS.rows != k) {
                int64_t stored_cols = neighbor_block.A_NS.cols;
                std::vector<int64_t> all_cols(stored_cols);
                std::iota(all_cols.begin(), all_cols.end(), 0);

                std::vector<DataType> sliced(k * stored_cols);
                extract_submatrix(
                    neighbor_block.A_NS.data.data(), neighbor_block.A_NS.rows,
                    box.skeleton_indices, all_cols,
                    sliced.data(), k
                );

                neighbor_block.A_NS.data = std::move(sliced);
                neighbor_block.A_NS.rows = k;
                neighbor_block.A_NS.lda  = k;
            }
        }
    }
}




template<typename CoordType, typename DataType, typename KernelType>
void compute_step_two_internal(
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    const std::vector<DataType>& X_BB,
    const std::vector<DataType>& A_NS_all_unmodified,
    const std::vector<DataType>& A_SN_all_unmodified,
    const std::vector<int64_t>& neighbor_point_counts,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    KernelType* kernel,
    std::vector<DataType>& X_NN_full,
    FactorizationThreadScratch<CoordType, DataType>& scratch,
    PendingFactorUpdates<DataType>* pending,
    bool store = false,
    int DEBUG = 0) {
    
    static_assert(std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>,
                  "Only double precision supported currently");
    
    std::unordered_map<int64_t, omp_lock_t*>& box_locks = level.box_locks;
    if (box->skeleton_indices.empty() || box->redundant_indices.empty()) {
        // level.eliminated_boxes.insert(box->morton_index);
        return;
    }
    
    int64_t k = box->skeleton_indices.size();
    int64_t r = box->redundant_indices.size();

    // Mark box as eliminated NOW (before updating neighbors)
    // This allows the slicing function to work correctly in Step 7
    // level.eliminated_boxes.insert(box->morton_index);
    
    int64_t total_neighbor_points = 0;
    if (box->X_NR.is_allocated()) {
        total_neighbor_points = box->X_NR.rows;
    }
    
    // ===== Step 1: Compute temp1 = -X_SR * X_RR^{-1} =====
    
    auto& temp1 = scratch.temp1;
    temp1.resize(static_cast<size_t>(k * r));
    std::copy(box->X_SR.data.begin(), box->X_SR.data.end(), temp1.begin());
    
    apply_right_inverse_in_place(
        box->X_RR, box->X_RR_pivots, factorization_method, temp1, k, r,
        "compute_step_two_internal temp1");
    
    for (int64_t i = 0; i < k * r; ++i) {
        temp1[i] = -temp1[i];
    }
    
    // ===== Step 2: Compute temp2 = -X_NR * X_RR^{-1} =====
    
    auto& temp2 = scratch.temp2;
    temp2.clear();
    if (total_neighbor_points > 0 && box->X_NR.is_allocated()) {
        temp2.resize(static_cast<size_t>(total_neighbor_points * r));
        std::copy(box->X_NR.data.begin(), box->X_NR.data.end(), temp2.begin());
        
        apply_right_inverse_in_place(
            box->X_RR, box->X_RR_pivots, factorization_method, temp2, total_neighbor_points, r,
            "compute_step_two_internal temp2");
        
        for (int64_t i = 0; i < total_neighbor_points * r; ++i) {
            temp2[i] = -temp2[i];
        }
    }

    // ===== Step 3: Get original X_RS =====
    
    auto& X_RS_original = scratch.x_rs_original;
    X_RS_original.resize(static_cast<size_t>(r * k));
    
    if (is_symmetric || is_hermitian) {
        // X_RS = X_SR^T (original X_SR)
        for (int64_t i = 0; i < r; ++i) {
            for (int64_t j = 0; j < k; ++j) {
                X_RS_original[i + j * r] = box->X_SR(j, i);
            }
        }
    } else {
        // X_RS is explicitly stored
        std::copy(box->X_RS.data.begin(), box->X_RS.data.end(), X_RS_original.begin());
    }
    
    // ===== Step 4: Compute Schur complement S = A_SS + temp1 * X_RS =====
    
    std::vector<DataType> A_SS(k * k);
    extract_submatrix(X_BB.data(), box->num_points,
                     box->skeleton_indices, box->skeleton_indices,
                     A_SS.data(), k);
    
    box->schur_complement.set_owned(k, k, std::move(A_SS), MatrixStorage<DataType>::FULL);
    
    int m = k, n = k, kk = r;
    DataType alpha = 1.0, beta = 1.0;
    
    gemm_("N", "N", &m, &n, &kk,
           &alpha, temp1.data(), &m,
           X_RS_original.data(), &kk,
           &beta, box->schur_complement.data.data(), &m);
    
    // ===== Step 5: Update A_NS = A_NS + temp2 * X_RS =====
    
    if (total_neighbor_points > 0 && !temp2.empty()) {
        // Compute full update: temp2 * X_RS
        auto& X_NS_update = scratch.x_ns_update;
        X_NS_update.resize(static_cast<size_t>(total_neighbor_points * k));
        
        m = total_neighbor_points; n = k; kk = r;
        alpha = 1.0; beta = 0.0;
        
        gemm_("N", "N", &m, &n, &kk,
               &alpha, temp2.data(), &m,
               X_RS_original.data(), &kk,
               &beta, X_NS_update.data(), &m);
        
        // Segment and update each neighbor
        int64_t current_row = 0;
        
        for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
            int64_t neighbor_morton = box->one_hop[idx];
            int64_t n_neighbor = neighbor_point_counts[idx];
            
            // Extract update for this neighbor (n_neighbor × k)
            auto& update_G = scratch.update_buffer;
            update_G.resize(static_cast<size_t>(n_neighbor * k));
            for (int64_t j = 0; j < k; ++j) {
                for (int64_t i = 0; i < n_neighbor; ++i) {
                    update_G[i + j * n_neighbor] = 
                        X_NS_update[(current_row + i) + j * total_neighbor_points];
                }
            }
            
            auto it = box->near_field_interaction_map.find(neighbor_morton);
            
            if (it != box->near_field_interaction_map.end()) {
                // Block exists - need double slicing then update
                int64_t block_idx = it->second;
                auto& modified_block = box->near_field_modified_interactions[block_idx];

                // DEBUG: Print dimensions
                if (DEBUG){
                    std::cout << "DEBUG Step 5: Slicing A_NS for neighbor " << neighbor_morton << std::endl;
                    std::cout << "  Stored A_NS dims: " << modified_block.A_NS.rows << " × " << modified_block.A_NS.cols << std::endl;
                    std::cout << "  Current box skeleton size: " << box->skeleton_indices.size() << std::endl;
                    std::cout << "  Expected neighbor size (n_neighbor): " << n_neighbor << std::endl;
                    std::cout << "  Neighbor is eliminated: " << (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end() ? "yes" : "no") << std::endl;
                }
                
                // Get fully sliced block
                std::vector<DataType> sliced = slice_modified_block_both_directions<CoordType, DataType>(
                    modified_block, level, neighbor_morton,
                    box->skeleton_indices, true, n_neighbor
                );
                
                // Size check
                if (sliced.size() != update_G.size()) {
                    throw std::runtime_error(
                        "compute_step_two_internal: A_NS size mismatch. Expected " +
                        std::to_string(update_G.size()) + " got " + std::to_string(sliced.size()));
                }
                
                // Add update
                for (int64_t i = 0; i < sliced.size(); ++i) {
                    sliced[i] += update_G[i];
                }
                
                // Store back (now properly sized: n_neighbor × k)
                modified_block.A_NS.set_owned(
                    n_neighbor, k, std::move(sliced), MatrixStorage<DataType>::FULL);
                
            } else {
               // Create new block - extract from unmodified A_NS and add update
                std::vector<DataType> A_NS_original(n_neighbor * k);
                
                // Extract this neighbor's portion from A_NS_all_unmodified
                for (int64_t j = 0; j < k; ++j) {
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        A_NS_original[i + j * n_neighbor] = 
                            A_NS_all_unmodified[(current_row + i) + j * total_neighbor_points];
                    }
                }

                if (DEBUG) {
                    bool original_has_nan = false;
                    for (auto val : A_NS_original) {
                        if (is_nan(val)) {
                            original_has_nan = true;
                            break;
                        }
                    }
                    std::cout << "  Step 5 (else): A_NS_original (before update) for neighbor " 
                            << neighbor_morton << " contains NaN: " << (original_has_nan ? "YES" : "NO") << std::endl;
                    
                    bool update_has_nan = false;
                    for (auto val : update_G) {
                        if (is_nan(val)) {
                            update_has_nan = true;
                            break;
                        }
                    }
                    std::cout << "  Step 5 (else): update_G contains NaN: " << (update_has_nan ? "YES" : "NO") << std::endl;
                }
                
                // Add update to original
                for (int64_t i = 0; i < n_neighbor * k; ++i) {
                    A_NS_original[i] += update_G[i];
                }

                if (DEBUG) {
                    bool final_has_nan = false;
                    for (auto val : A_NS_original) {
                        if (is_nan(val)) {
                            final_has_nan = true;
                            break;
                        }
                    }
                    std::cout << "  Step 5 (else): A_NS_original (after update) contains NaN: " << (final_has_nan ? "YES" : "NO") << std::endl;
                }
                
                // Store combined result
                ModifiedBlock<DataType> new_block;
                new_block.neighbor_morton = neighbor_morton;
                new_block.A_NS.set_owned(
                    n_neighbor, k, std::move(A_NS_original), MatrixStorage<DataType>::FULL);
                
                int64_t new_idx = box->near_field_modified_interactions.size();
                box->near_field_modified_interactions.push_back(std::move(new_block));
                box->near_field_interaction_map[neighbor_morton] = new_idx;
            }
            
            current_row += n_neighbor;
        }
    }
    
    // ===== Step 6: Nonsymmetric - compute temp3, temp4 and update A_SN =====
    // ===== Declare temp3 and temp4 for nonsymmetric case =====
    
    auto& temp3 = scratch.temp3;
    auto& temp4 = scratch.temp4;
    temp3.clear();
    temp4.clear();
    
    if (!is_symmetric && !is_hermitian) {
        // temp3 = -X_RR^{-1} * X_RS
        // std::vector<DataType> temp3(r * k);
        temp3.resize(r * k);
        std::copy(box->X_RS.data.begin(), box->X_RS.data.end(), temp3.begin());
        
        apply_left_inverse_in_place(
            box->X_RR, box->X_RR_pivots, factorization_method, temp3, r, k,
            "compute_step_two_internal temp3");
        
        for (int64_t i = 0; i < r * k; ++i) {
            temp3[i] = -temp3[i];
        }
        
        // temp4 = -X_RR^{-1} * X_RN
        // std::vector<DataType> temp4;
        if (total_neighbor_points > 0 && box->X_RN.is_allocated()) {
            temp4.resize(r * total_neighbor_points);
            std::copy(box->X_RN.data.begin(), box->X_RN.data.end(), temp4.begin());
            
            apply_left_inverse_in_place(
                box->X_RR, box->X_RR_pivots, factorization_method, temp4, r, total_neighbor_points,
                "compute_step_two_internal temp4");
            
            for (int64_t i = 0; i < r * total_neighbor_points; ++i) {
                temp4[i] = -temp4[i];
            }
            
            // Compute A_SN update = X_SR * temp4
            auto& X_SN_update = scratch.x_sn_update;
            X_SN_update.resize(static_cast<size_t>(k * total_neighbor_points));
            
            int m = k, n = total_neighbor_points, kk = r;
            DataType alpha = 1.0, beta = 0.0;
            
            gemm_("N", "N", &m, &n, &kk,
                   &alpha, box->X_SR.data.data(), &m,
                   temp4.data(), &kk,
                   &beta, X_SN_update.data(), &m);
            
            // Segment and update each neighbor
            int64_t current_col = 0;
            
            for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
                int64_t neighbor_morton = box->one_hop[idx];
                int64_t n_neighbor = neighbor_point_counts[idx];
                
                // Extract update for this neighbor (k × n_neighbor)
                auto& update_G = scratch.update_buffer;
                update_G.resize(static_cast<size_t>(k * n_neighbor));
                for (int64_t j = 0; j < n_neighbor; ++j) {
                    for (int64_t i = 0; i < k; ++i) {
                        update_G[i + j * k] = 
                            X_SN_update[i + (current_col + j) * k];
                    }
                }
                
                auto it = box->near_field_interaction_map_nonsymmetry.find(neighbor_morton);
                
                if (it != box->near_field_interaction_map_nonsymmetry.end()) {
                    // Block exists - need double slicing then update
                    int64_t block_idx = it->second;
                    auto& modified_block = box->near_field_modified_interactions[block_idx];
                    
                    if (!modified_block.A_SN.is_allocated()) {
                        throw std::runtime_error(
                            "compute_step_two_internal: A_SN not allocated");
                    }
                    
                    // Get fully sliced block
                    std::vector<DataType> sliced = slice_modified_block_both_directions<CoordType, DataType>(
                        modified_block, level, neighbor_morton,
                        box->skeleton_indices, false, n_neighbor
                    );
                    
                    // Size check
                    if (sliced.size() != update_G.size()) {
                        throw std::runtime_error(
                            "compute_step_two_internal: A_SN size mismatch. Expected " +
                            std::to_string(update_G.size()) + " got " + std::to_string(sliced.size()));
                    }
                    
                    // Add update
                    for (int64_t i = 0; i < sliced.size(); ++i) {
                        sliced[i] += update_G[i];
                    }
                    
                    // Store back (now properly sized: k × n_neighbor)
                    modified_block.A_SN.set_owned(
                        k, n_neighbor, std::move(sliced), MatrixStorage<DataType>::FULL);
                    
                } else {
                    // Create new block - extract from unmodified A_SN and add update
                    std::vector<DataType> A_SN_original(k * n_neighbor);
                    
                    // Extract this neighbor's portion from A_SN_all_unmodified
                    for (int64_t j = 0; j < n_neighbor; ++j) {
                        for (int64_t i = 0; i < k; ++i) {
                            A_SN_original[i + j * k] = 
                                A_SN_all_unmodified[i + (current_col + j) * k];
                        }
                    }
                    
                    // Add update to original
                    for (int64_t i = 0; i < k * n_neighbor; ++i) {
                        A_SN_original[i] += update_G[i];
                    }
                    
                    // Store combined result
                    ModifiedBlock<DataType> new_block;
                    new_block.neighbor_morton = neighbor_morton;
                    new_block.A_SN.set_owned(
                        k, n_neighbor, std::move(A_SN_original), MatrixStorage<DataType>::FULL);
                    
                    int64_t new_idx = box->near_field_modified_interactions.size();
                    box->near_field_modified_interactions.push_back(std::move(new_block));
                    box->near_field_interaction_map_nonsymmetry[neighbor_morton] = new_idx;
                }
                
                current_col += n_neighbor;
            }
        }
        

    }

    // Helper: find assisting request for a morton (nullptr if not present)
    auto find_assist = [&](int64_t morton) -> PointDataRequest<CoordType>* {
        auto it = level.assisting_box_points_for_kernel_evaluation.find(morton);
        if (it == level.assisting_box_points_for_kernel_evaluation.end()) return nullptr;
        return &level.assisting_boxes[it->second];
    };
    

    
    // ===== Step 7: Pairwise update - neighbors' view of current box B =====
    // Update each neighbor G's modified_interactions to include their view of B.
    // When store=true this entire stage is deferred so the later owner replay can
    // keep the symmetric path lock-free.
    if (!store && total_neighbor_points > 0) {
        int64_t current_row = 0;
        
        for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
            int64_t neighbor_morton = box->one_hop[idx];
            int64_t n_neighbor = neighbor_point_counts[idx];
            
            // Find neighbor G (OWNED?) and (OPTIONAL ghost/assist presence)
            BoxData<CoordType, DataType>* neighbor_local = level.find_local_box(neighbor_morton);
            BoxData<CoordType, DataType>* neighbor_box = neighbor_local;

            if (neighbor_box == nullptr) {
                // keep ghost lookup intact (but we will NOT write into ghosts)
                neighbor_box = level.find_ghost_box(neighbor_morton);
            }

            if (neighbor_local == nullptr) {
                // no-halo model: every remote box must be present as assisting
                if (find_assist(neighbor_morton) == nullptr) {
                    throw std::runtime_error("Step 7: remote neighbor missing assisting entry");
                }
            }

            // If neighbor is NOT owned here, we must cache Step-7 REPLACE (do NOT skip!)
            if (neighbor_local == nullptr) {
                if (pending == nullptr) {
                    throw std::runtime_error("Step 7: pending is null but remote REPLACE is required");
                }

                // Extract X_GS from B's updated A_NS (already has Schur updates)
                auto it = box->near_field_interaction_map.find(neighbor_morton);
                if (it == box->near_field_interaction_map.end()) {
                    throw std::runtime_error("Step 7: neighbor missing from near_field_interaction_map");
                }
                const auto& modified_block = box->near_field_modified_interactions[it->second];
                if (!modified_block.A_NS.is_allocated()) {
                    throw std::runtime_error("Step 7: A_NS should be allocated after step 5");
                }

                if (is_symmetric || is_hermitian) {
                    // Cache REPLACE of neighbor's A_NS for (neighbor_morton <- box->morton_index)
                    // This is exactly what you'd store locally into neighbor_box->A_NS, but remote.
                    std::vector<DataType> X_GS_transpose(k * n_neighbor);
                    for (int64_t i = 0; i < k; ++i) {
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            X_GS_transpose[i + j * k] = modified_block.A_NS.data[j + i * n_neighbor];
                        }
                    }

                    pending->replace_blocks[ReplaceKey{neighbor_morton, box->morton_index}] =
                        DenseBlock<DataType>{k, n_neighbor, std::move(X_GS_transpose)};
                } else {
                    // Nonsymmetric: cache REPLACE into neighbor's A_SN
                    std::vector<DataType> X_GS = modified_block.A_NS.data; // (n_neighbor x k)

                    pending->replace_blocks[ReplaceKey{neighbor_morton, box->morton_index}] =
                        DenseBlock<DataType>{n_neighbor, k, std::move(X_GS)};
                }

                current_row += n_neighbor;
                continue;
            }

            // If we get here, neighbor is owned locally => keep your existing local-update logic
            BoxData<CoordType, DataType>* neighbor_box_owned = neighbor_local;
            neighbor_box = neighbor_box_owned;
            
            // Extract X_GS from B's updated A_NS (already has Schur updates)
            auto it = box->near_field_interaction_map.find(neighbor_morton);
            if (it == box->near_field_interaction_map.end()) {
                throw std::runtime_error(
                    "compute_step_two_internal: Neighbor should be in near_field_interaction_map");
            }
            
            const auto& modified_block = box->near_field_modified_interactions[it->second];
            


            if (!modified_block.A_NS.is_allocated()) {
                throw std::runtime_error(
                    "compute_step_two_internal: A_NS should be allocated after step 5");
            }
            
            

            
            if (is_symmetric || is_hermitian) {
                // ===== Symmetric: Store X_GS^T in G's A_NS =====
                // From B: X_GS is (n_neighbor × k)
                // From G: need (n_neighbor × k) but as A_NS[B]
                
                std::vector<DataType> X_GS_transpose(k * n_neighbor);
                for (int64_t i = 0; i < k; ++i) {
                    for (int64_t j = 0; j < n_neighbor; ++j) {
                        X_GS_transpose[i + j * k] = modified_block.A_NS.data[j + i * n_neighbor];
                    }
                }
                // Write to neighbor INSIDE lock
                omp_set_lock(box_locks.at(neighbor_morton));

                auto it_G = neighbor_box->near_field_interaction_map.find(box->morton_index);
                
                if (it_G != neighbor_box->near_field_interaction_map.end()) {
                    // Block exists - replace with new view from B

                    int64_t block_idx_G = it_G->second;
                    auto& block_G = neighbor_box->near_field_modified_interactions[block_idx_G];
                    
                    // DEBUG: Print info
                    if (DEBUG) {
                        std::cout << "DEBUG Step 7: UPDATING existing block in Box " << neighbor_morton 
                                << "'s view for Box " << box->morton_index << std::endl;
                        std::cout << "  Old A_NS dims: " << block_G.A_NS.rows << " × " << block_G.A_NS.cols << std::endl;
                        std::cout << "  New A_NS dims: " << k << " × " << n_neighbor << std::endl;
                    }
                    block_G.A_NS.set_owned(
                        k, n_neighbor, std::move(X_GS_transpose), MatrixStorage<DataType>::FULL);

                    // DEBUG: Verify allocation
                    if (DEBUG) {
                        bool stored_has_nan = false;
                        for (auto val : block_G.A_NS.data) {
                            if (is_nan(val)) {
                                stored_has_nan = true;
                                break;
                            }
                        }
                        std::cout << "  AFTER STORAGE in Box " << neighbor_morton << "'s block[" << block_idx_G << "]:" << std::endl;
                        std::cout << "    A_NS contains NaN: " << (stored_has_nan ? "YES ⚠️" : "NO") << std::endl;
                        std::cout << "    A_NS dims: " << block_G.A_NS.rows << " × " << block_G.A_NS.cols << std::endl;
                        std::cout << "    A_NS data.size: " << block_G.A_NS.data.size() << std::endl;
                        
                        if (stored_has_nan) {
                            std::cout << "  ⚠️ NaN DETECTED! Box " << box->morton_index 
                                    << " corrupted Box " << neighbor_morton << "'s view!" << std::endl;
                        }
                    }
                } else {
                    // Create new block in G's modified_interactions
                    // DEBUG: Print info
                    if (DEBUG){
                        std::cout << "DEBUG Step 7: Creating NEW block in Box " << neighbor_morton 
                                << "'s view for Box " << box->morton_index << std::endl;
                        std::cout << "  Allocating A_NS: " << k << " × " << n_neighbor << std::endl;
                        bool transpose_has_nan = false;
                        for (auto val : X_GS_transpose) {
                            if (is_nan(val)) {
                                transpose_has_nan = true;
                                break;
                            }
                        }
                        std::cout << "  X_GS_transpose contains NaN: " << (transpose_has_nan ? "YES" : "NO") << std::endl;
                    }
                    ModifiedBlock<DataType> new_block_G;
                    new_block_G.neighbor_morton = box->morton_index;
                    new_block_G.A_NS.set_owned(
                        k, n_neighbor, std::move(X_GS_transpose), MatrixStorage<DataType>::FULL);  // FIX: k × n_neighbor

                    
                    int64_t new_idx = neighbor_box->near_field_modified_interactions.size();
                    neighbor_box->near_field_modified_interactions.push_back(std::move(new_block_G));
                    neighbor_box->near_field_interaction_map[box->morton_index] = new_idx;

                    if (DEBUG) {
                        const auto& stored_block = neighbor_box->near_field_modified_interactions[new_idx];
                        bool stored_has_nan = false;
                        for (auto val : stored_block.A_NS.data) {
                            if (is_nan(val)) {
                                stored_has_nan = true;
                                break;
                            }
                        }
                        std::cout << "  AFTER STORAGE in Box " << neighbor_morton << "'s NEW block[" << new_idx << "]:" << std::endl;
                        std::cout << "    A_NS contains NaN: " << (stored_has_nan ? "YES ⚠️" : "NO") << std::endl;
                        
                        if (stored_has_nan) {
                            std::cout << "  ⚠️ NaN DETECTED! Box " << box->morton_index 
                                    << " corrupted Box " << neighbor_morton << "'s view!" << std::endl;
                        }
                    }
                }
                omp_unset_lock(box_locks.at(neighbor_morton));
            } else {
                // ===== Nonsymmetric: Store X_GS in G's A_SN =====
                // From B: X_GS is (n_neighbor × k) in A_NS
                // From G: X_GS is (n_neighbor × k) in A_SN
                std::vector<DataType> X_GS = modified_block.A_NS.data;
                
                auto it_G = neighbor_box->near_field_interaction_map_nonsymmetry.find(box->morton_index);
                
                if (it_G != neighbor_box->near_field_interaction_map_nonsymmetry.end()) {
                    // Block exists - replace with new view from B
                    int64_t block_idx_G = it_G->second;
                    auto& block_G = neighbor_box->near_field_modified_interactions[block_idx_G];
                    
                    block_G.A_SN.set_owned(
                        n_neighbor, k, std::move(X_GS), MatrixStorage<DataType>::FULL);
                    
                } else {
                    // Create new block
                    ModifiedBlock<DataType> new_block_G;
                    new_block_G.neighbor_morton = box->morton_index;
                    new_block_G.A_SN.set_owned(
                        n_neighbor, k, std::move(X_GS), MatrixStorage<DataType>::FULL);
                    
                    int64_t new_idx = neighbor_box->near_field_modified_interactions.size();
                    neighbor_box->near_field_modified_interactions.push_back(std::move(new_block_G));
                    neighbor_box->near_field_interaction_map_nonsymmetry[box->morton_index] = new_idx;
                }
            }
            
            current_row += n_neighbor;
        }
    }
    
    // ===== Step 8: Compute X_NN = A_NN + temp2 * X_RN =====
    
    // std::vector<DataType> X_NN_full;
    
    // In deferred mode (store=true), Step 8/9 are reconstructed later from the
    // saved source state instead of writing inline under box locks here.
    if (!store) {
        if (total_neighbor_points > 0) {
            auto gather_start = std::chrono::high_resolution_clock::now();
            X_NN_full.resize(total_neighbor_points * total_neighbor_points);
            
            // We will calculate: X_NN = temp2 * X_RN
            // where temp2 is (total_neighbor_points × r)
            // and X_RN is (r × total_neighbor_points).
            // The result X_NN is (total_neighbor_points × total_neighbor_points).
            
            int m = total_neighbor_points; // Rows of op(A) and C
            int n = total_neighbor_points; // Cols of op(B) and C
            int k = r;                     // Cols of op(A) and Rows of op(B)
            DataType alpha = 1.0, beta = 0.0;
            
            if (is_symmetric || is_hermitian) {
                // --- Symmetric Case ---
                // We need X_RN, but we only store X_NR. Since X_RN = X_NR^T,
                // we can compute temp2 * X_NR^T directly using BLAS.
                // Matrix A: temp2 (m x k)
                // Matrix B: box->X_NR (n x k), but we use its transpose (k x n).
                
                // Leading dimension of A (temp2)
                int lda = m; 
                // Leading dimension of B (box->X_NR), which has n rows
                int ldb = n; 
                // Leading dimension of C (X_NN_full)
                int ldc = m;
                
                // dgemm computes C = alpha * op(A) * op(B) + beta * C
                // op(A) = temp2, op(B) = X_NR^T
                // printf("before total thread: %d\n", omp_get_num_threads());
                // omp_set_num_threads(16);
                // printf("after total thread: %d\n", omp_get_num_threads());
                gemm_("N", "T", &m, &n, &k,
                    &alpha, 
                    temp2.data(), &lda,
                    box->X_NR.data.data(), &ldb,
                    &beta, 
                    X_NN_full.data(), &ldc);
                // auto gather_end = std::chrono::high_resolution_clock::now();
                // auto gather_duration = std::chrono::duration_cast<std::chrono::milliseconds>(gather_end - gather_start);
                // std::cout << "  compute time: " << gather_duration.count() << " ms" << std::endl;
                // printf("m: %d, n: %d, k: %d, total neighbor: %d\n", m, n, k, total_neighbor_points);
                    
            } else {
                // --- Nonsymmetric Case ---
                // We use the explicitly stored X_RN matrix.
                // Matrix A: temp2 (m x k)
                // Matrix B: box->X_RN (k x n)
                
                // Leading dimension of A (temp2)
                int lda = m; 
                // Leading dimension of B (box->X_RN), which has k rows
                int ldb = k; 
                // Leading dimension of C (X_NN_full)
                int ldc = m;
                
                // dgemm computes C = alpha * op(A) * op(B) + beta * C
                // op(A) = temp2, op(B) = X_RN
                gemm_("N", "N", &m, &n, &k,
                    &alpha, 
                    temp2.data(), &lda,
                    box->X_RN.data.data(), &ldb,
                    &beta, 
                    X_NN_full.data(), &ldc);
            }
        }
    
        // ===== Step 9: Update X_NN for all neighbor pairs =====
        // Debug
        if (DEBUG) {
            std::cout << "\n=== Step 9: Processing Box " << box->morton_index << " X_NN updates ===" << std::endl;
            std::cout << "total_neighbor_points = " << total_neighbor_points << std::endl;
            std::cout << "X_NN_full.size() = " << X_NN_full.size() << std::endl;
        }

        auto accumulate_dense = [&](DenseBlock<DataType>& dst,
                                int64_t rows, int64_t cols,
                                const std::vector<DataType>& delta)
        {
            if (delta.size() != static_cast<size_t>(rows * cols)) {
                throw std::runtime_error("accumulate_dense: delta size mismatch");
            }
            if (dst.data.empty()) {
                dst.rows = rows;
                dst.cols = cols;
                dst.data = delta; // first contribution
                return;
            }
            if (dst.rows != rows || dst.cols != cols) {
                throw std::runtime_error("accumulate_dense: dimension mismatch");
            }
            for (size_t i = 0; i < dst.data.size(); ++i) dst.data[i] += delta[i];
        };

        
        
        if (total_neighbor_points > 0 && !X_NN_full.empty()) {

            if (DEBUG) {
                // Check for NaN in X_NN_full
                bool has_nan = false;
                for (auto val : X_NN_full) {
                    if (is_nan(val)) {
                        has_nan = true;
                        break;
                    }
                }
                std::cout << "DEBUG Step 9: X_NN_full contains NaN: " << (has_nan ? "YES" : "NO") << std::endl;
                
                // Check temp2
                bool temp2_has_nan = false;
                for (auto val : temp2) {
                    if (is_nan(val)) {
                        temp2_has_nan = true;
                        break;
                    }
                }
                std::cout << "DEBUG Step 9: temp2 contains NaN: " << (temp2_has_nan ? "YES" : "NO") << std::endl;
            }
            
            // Determine dimension for Morton decoding
            int dimension = (box->bounds[4] == box->bounds[5]) ? 2 : 3;
            
            // Column-major traversal
            int64_t col_offset = 0;
            
            for (size_t j_idx = 0; j_idx < box->one_hop.size(); ++j_idx) {
                int64_t RG_morton = box->one_hop[j_idx];
                int64_t n_RG = neighbor_point_counts[j_idx];
                
                // Traverse rows
                size_t i_start = (is_symmetric || is_hermitian) ? j_idx : 0;
                
                // Initialize row_offset to skip rows before i_start
                int64_t row_offset = 0;
                for (size_t skip = 0; skip < i_start; ++skip) {
                    row_offset += neighbor_point_counts[skip];
                }
                
                for (size_t i_idx = i_start; i_idx < box->one_hop.size(); ++i_idx) {
                    int64_t LG_morton = box->one_hop[i_idx];
                    int64_t n_LG = neighbor_point_counts[i_idx];
                    
                    
                    // Extract X(LG, RG) from X_NN_full
                    std::vector<DataType> X_LG_RG(n_LG * n_RG);

                    for (int64_t col = 0; col < n_RG; ++col) {
                        for (int64_t row = 0; row < n_LG; ++row) {
                            X_LG_RG[row + col * n_LG] =
                                X_NN_full[(row_offset + row) + (col_offset + col) * total_neighbor_points];
                        }
                    }
                
                    
                    // Precompute transpose used by LG-side storage (shape n_RG × n_LG, column-major)
                    std::vector<DataType> X_LG_RG_for_LG(n_LG * n_RG);
                    for (int64_t i = 0; i < n_LG; ++i) {
                        for (int64_t j = 0; j < n_RG; ++j) {
                            X_LG_RG_for_LG[j + i * n_RG] = X_LG_RG[i + j * n_LG];
                        }
                    }

                    if (DEBUG) {
                        bool xlg_rg_has_nan = false;
                        for (auto val : X_LG_RG) {
                            if (is_nan(val)) {
                                xlg_rg_has_nan = true;
                                break;
                            }
                        }
                        std::cout << "DEBUG Step 9: Extracted X_LG_RG for LG=" << LG_morton 
                                << " RG=" << RG_morton << std::endl;
                        std::cout << "  X_LG_RG size: " << X_LG_RG.size() 
                                << " (n_LG=" << n_LG << " × n_RG=" << n_RG << ")" << std::endl;
                        std::cout << "  X_LG_RG contains NaN: " << (xlg_rg_has_nan ? "YES" : "NO") << std::endl;
                        std::cout << "  Extracted from X_NN_full at row_offset=" << row_offset 
                                << " col_offset=" << col_offset << std::endl;
                    }
                    
                    // Check if diagonal
                    bool is_diagonal = (LG_morton == RG_morton);

                    if (DEBUG) {
                        std::cout << "Step 9 loop: i_idx=" << i_idx << " j_idx=" << j_idx 
                                << " LG=" << LG_morton << " RG=" << RG_morton 
                                << " is_diagonal=" << is_diagonal << std::endl;
                    }
                    
                    // Find LG and RG endpoints (OWNERSHIP-BASED; no writing to ghosts)
                    BoxData<CoordType, DataType>* LG_local = level.find_local_box(LG_morton);
                    BoxData<CoordType, DataType>* RG_local = level.find_local_box(RG_morton);

                    const bool LG_owned_here = (LG_local != nullptr);
                    const bool RG_owned_here = (RG_local != nullptr);

                    // Keep these names because they are used later in your code
                    BoxData<CoordType, DataType>* LG_box = LG_local;
                    BoxData<CoordType, DataType>* RG_box = RG_local;

                    PointDataRequest<CoordType>* LG_assisting = nullptr;
                    PointDataRequest<CoordType>* RG_assisting = nullptr;

                    // In the no-halo model: if not owned, treat as assisting (must exist)
                    bool LG_is_assisting = !LG_owned_here;
                    bool RG_is_assisting = !RG_owned_here;

                    if (!LG_owned_here) {
                        LG_assisting = find_assist(LG_morton);
                        if (LG_assisting == nullptr) {
                            // (Optional) legacy fallback to ghost pointer (read-only), but this should not happen in no-halo mode
                            LG_box = level.find_ghost_box(LG_morton);
                            throw std::runtime_error("Step 9: nonlocal LG missing assisting data");
                        }
                    } else {
                        // (Optional) legacy ghost lookup (read-only); not used for writes
                        if (LG_box == nullptr) LG_box = level.find_ghost_box(LG_morton);
                    }

                    if (!RG_owned_here) {
                        RG_assisting = find_assist(RG_morton);
                        if (RG_assisting == nullptr) {
                            RG_box = level.find_ghost_box(RG_morton);
                            throw std::runtime_error("Step 9: nonlocal RG missing assisting data");
                        }
                    } else {
                        if (RG_box == nullptr) RG_box = level.find_ghost_box(RG_morton);
                    }
                    
                    // Skip if both are neither local nor ghost nor assisting
                    if ((LG_box == nullptr && !LG_is_assisting) || (RG_box == nullptr && !RG_is_assisting)) {
                        throw std::runtime_error("box missing in level during X_NN update");
                    }

                    // The deferred owner pass will rebuild every local-local
                    // storage update later. Nonlocal pairs still need their
                    // canonical delta recorded now for remote ownership.
                    if (store && LG_owned_here && RG_owned_here) {
                        row_offset += n_LG;
                        continue;
                    }
                    
                    // NOTE: do NOT skip remote-remote here; Step 9 deltas must be accumulated via `pending`.
                    // Local writes below are guarded by LG_owned_here / RG_owned_here.
                    
                    // Determine 1-hop or 2-hop
                    bool is_one_hop = false;
                    
                    if (dimension == 2) {
                        uint32_t LG_x, LG_y, RG_x, RG_y;
                        morton::decode_2d(LG_morton, LG_x, LG_y);
                        morton::decode_2d(RG_morton, RG_x, RG_y);
                        
                        int64_t dx = std::abs(static_cast<int64_t>(LG_x) - static_cast<int64_t>(RG_x));
                        int64_t dy = std::abs(static_cast<int64_t>(LG_y) - static_cast<int64_t>(RG_y));
                        
                        is_one_hop = (dx <= 1 && dy <= 1);
                    } else {
                        uint32_t LG_x, LG_y, LG_z, RG_x, RG_y, RG_z;
                        morton::decode_3d(LG_morton, LG_x, LG_y, LG_z);
                        morton::decode_3d(RG_morton, RG_x, RG_y, RG_z);
                        
                        int64_t dx = std::abs(static_cast<int64_t>(LG_x) - static_cast<int64_t>(RG_x));
                        int64_t dy = std::abs(static_cast<int64_t>(LG_y) - static_cast<int64_t>(RG_y));
                        int64_t dz = std::abs(static_cast<int64_t>(LG_z) - static_cast<int64_t>(RG_z));
                        
                        is_one_hop = (dx <= 1 && dy <= 1 && dz <= 1);
                    }

                    // ---- Remote accumulation hook (Step 9) ----
                    // If a target endpoint is not owned here, you must accumulate/send an ADD update instead of writing locally.
                    // Use:
                    //   - LG-side delta: X_LG_RG_for_LG  (shape n_RG × n_LG)
                    //   - RG-side delta: X_LG_RG         (shape n_LG × n_RG)   (skip if diagonal)
                    //
                    // if (pending) {
                    //     if (!LG_owned_here) pending->accumulate_nn_add(/*target=*/LG_morton, /*other=*/RG_morton,
                    //                                                   /*is_diagonal=*/is_diagonal, /*is_one_hop=*/is_one_hop,
                    //                                                   /*slot=*/(is_symmetric||is_hermitian ? A_NS : A_SN),
                    //                                                   /*rows=*/n_RG, /*cols=*/n_LG, X_LG_RG_for_LG);
                    //     if (!RG_owned_here && !is_diagonal) pending->accumulate_nn_add(/*target=*/RG_morton, /*other=*/LG_morton,
                    //                                                                   /*is_diagonal=*/false, /*is_one_hop=*/is_one_hop,
                    //                                                                   /*slot=*/A_NS,
                    //                                                                   /*rows=*/n_LG, /*cols=*/n_RG, X_LG_RG);
                    // }
                    
                    // ===== Update LG's view of RG (skip if LG is assisting) =====

                    

                    // ---- Step 9: accumulate canonical ADD delta for remote delivery ----
                    // Canonical key: (lo,hi). Stored block is ΔM_{lo→hi} with dims (n_hi × n_lo)
                    // i.e., "lo's view of hi" (rows=hi, cols=lo).
                    if (pending != nullptr) {

                        // Only worth caching if at least one endpoint is nonlocal (otherwise local update already did it)
                        if (!LG_owned_here || !RG_owned_here) {
                            const bool is_diagonal = (LG_morton == RG_morton);

                            EdgeKind kind = EdgeKind::Far;
                            if (is_diagonal) kind = EdgeKind::Diag;
                            else kind = is_one_hop ? EdgeKind::Near : EdgeKind::Far;

                            const int64_t lo = std::min(LG_morton, RG_morton);
                            const int64_t hi = std::max(LG_morton, RG_morton);

                            const int64_t n_lo = (lo == LG_morton) ? n_LG : n_RG;
                            const int64_t n_hi = (hi == LG_morton) ? n_LG : n_RG;

                            // Pick the already-correct orientation for ΔM_{lo→hi} (n_hi × n_lo):
                            // - If LG is lo, then X_LG_RG_for_LG is (n_RG × n_LG) = (n_hi × n_lo).
                            // - If LG is hi, then X_LG_RG is (n_LG × n_RG) = (n_hi × n_lo).
                            std::vector<DataType>& delta_lo_hi =
                                (LG_morton == lo) ? X_LG_RG_for_LG : X_LG_RG;
                            // for(int gg = 0; gg < delta_lo_hi.size(); ++gg) {
                            //     delta_lo_hi[gg] *= 10.0;
                            // }
                            

                            auto& dst = pending->accumulated_deltas[EdgeKey{lo, hi, kind}];
                            accumulate_dense(dst, n_hi, n_lo, delta_lo_hi);
                            
                        }
                    }

                    // ===== Update LG's view of RG (owned only) =====
                    if (LG_owned_here) {
                        // X_LG_RG_for_LG already computed above
                        omp_set_lock(box_locks.at(LG_morton)); 
                        if (is_diagonal) {
                            // Update schur_complement - LG_box must be valid here
                            if (!LG_box->schur_complement.is_allocated()) {
                                // Get LG coordinates
                                std::vector<CoordType> LG_coords;
                                bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                
                                if (LG_eliminated) {
                                    LG_coords.resize(n_LG * dimension);
                                    for (int64_t i = 0; i < n_LG; ++i) {
                                        int64_t src_idx = LG_box->skeleton_indices[i];
                                        for (int d = 0; d < dimension; ++d) {
                                            LG_coords[i * dimension + d] = 
                                                LG_box->point_coords[src_idx * dimension + d];
                                        }
                                    }
                                } else {
                                    LG_coords = LG_box->point_coords;
                                }
                                
                                // Evaluate A(LG, LG)
                                std::vector<DataType> A_LG_LG(n_LG * n_LG);
                                kernel->evaluate_block(
                                    LG_coords.data(), n_LG,
                                    LG_coords.data(), n_LG,
                                    A_LG_LG.data(), n_LG
                                );
                                
                                LG_box->schur_complement.set_owned(
                                    n_LG, n_LG, std::move(A_LG_LG), MatrixStorage<DataType>::FULL);
                            }
                            
                            // Add X_NN update
                            for (int64_t i = 0; i < n_LG * n_LG; ++i) {
                                LG_box->schur_complement.data[i] += X_LG_RG_for_LG[i];
                            }
                            
                        } else {
                            
                            // Update near_field or far_field
                            auto& interaction_map = is_one_hop ? 
                                (is_symmetric || is_hermitian ? 
                                    LG_box->near_field_interaction_map : 
                                    LG_box->near_field_interaction_map_nonsymmetry) :
                                (is_symmetric || is_hermitian ?
                                    LG_box->far_field_interaction_map :
                                    LG_box->far_field_interaction_map_nonsymmetry);
                            
                            auto& modified_interactions = is_one_hop ?
                                LG_box->near_field_modified_interactions :
                                LG_box->far_field_modified_interactions;
                            
                            auto it_LG = interaction_map.find(RG_morton);
                            
                            if (it_LG != interaction_map.end()) {
                                // Block exists - need to add update
                                int64_t block_idx = it_LG->second;
                                auto& block = modified_interactions[block_idx];
                                
                                auto& target_matrix = (is_symmetric || is_hermitian) ?
                                    block.A_NS : block.A_SN;
                                
                                // Slice if needed
                                bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                bool RG_eliminated = (level.eliminated_boxes.find(RG_morton) != level.eliminated_boxes.end());
                                
                                std::vector<DataType> sliced;

                                if (LG_eliminated && RG_eliminated) {
                                    
                                    sliced = slice_modified_block_both_directions<CoordType, DataType>(
                                        block, level, RG_morton, 
                                        LG_box->skeleton_indices,
                                        (is_symmetric || is_hermitian),
                                        n_RG
                                    );
                                    
                                } else if (LG_eliminated) {
                                    int64_t stored_rows = target_matrix.rows;
                                    int64_t stored_cols = target_matrix.cols;
                                    
                                    if (stored_rows == n_RG && stored_cols == n_LG) {
                                        sliced = target_matrix.data;
                                    } else if (stored_rows == n_RG && stored_cols > n_LG) {
                                        sliced.resize(n_RG * n_LG);
                                        for (int64_t j = 0; j < n_LG; ++j) {
                                            int64_t src_col = LG_box->skeleton_indices[j];
                                            for (int64_t i = 0; i < n_RG; ++i) {
                                                sliced[i + j * n_RG] = target_matrix(i, src_col);
                                            }
                                        }
                                    } else {
                                        printf("Step 9 LG view LG_elim: stored_rows=%lld stored_cols=%lld n_RG=%lld n_LG=%lld LG_morton=%lld RG_morton=%lld LG_eliminated=%lld RG_eliminated=%lld\n",
                                            stored_rows, stored_cols, n_RG, n_LG, LG_morton, RG_morton, LG_eliminated, RG_eliminated);
                                        throw std::runtime_error("Step 9 LG view LG_elim: unexpected dims");

                                    }
                                } else if (RG_eliminated) {
                                    int64_t stored_rows = target_matrix.rows;
                                    int64_t stored_cols = target_matrix.cols;
                                    
                                    if (stored_rows == n_RG && stored_cols == n_LG) {
                                        sliced = target_matrix.data;
                                    } else if (stored_rows > n_RG && stored_cols == n_LG) {
                                        sliced.resize(n_RG * n_LG);

                                        // Pick the skeleton index source once (local/ghost RG_box vs assisting RG box)
                                        const std::vector<int64_t>* rg_skel = nullptr;
                                        if (RG_box != nullptr) {
                                            rg_skel = &RG_box->skeleton_indices;
                                        } else {
                                            auto it = level.assisting_box_points_for_kernel_evaluation.find(RG_morton);
                                            if (it == level.assisting_box_points_for_kernel_evaluation.end()) {
                                                throw std::runtime_error("Missing assisting skeleton indices for RG_morton=" +
                                                                        std::to_string(RG_morton));
                                            }
                                            const auto idx = static_cast<size_t>(it->second);
                                            if (idx >= level.assisting_boxes.size()) {
                                                throw std::runtime_error("Assisting index out of range for RG_morton=" +
                                                                        std::to_string(RG_morton));
                                            }
                                            rg_skel = &level.assisting_boxes[idx].skel_indices;
                                        }

                                        if (static_cast<int64_t>(rg_skel->size()) < n_RG) {
                                            throw std::runtime_error("RG skeleton size < n_RG for RG_morton=" +
                                                                    std::to_string(RG_morton));
                                        }

                                        const auto& sk = *rg_skel;

                                        for (int64_t j = 0; j < n_LG; ++j) {
                                            for (int64_t i = 0; i < n_RG; ++i) {
                                                const int64_t src_row = sk[static_cast<size_t>(i)];
                                                sliced[static_cast<size_t>(i + j * n_RG)] = target_matrix(src_row, j);
                                            }
                                        }
                                    } else {
                                        throw std::runtime_error("Step 9 LG update: unexpected dimensions");
                                    }
                                } else {
                                    // Neither eliminated (or RG is assisting)
                                    sliced = target_matrix.data;
                                }
                                
                                if (sliced.size() != X_LG_RG_for_LG.size()) {
                                    // printf("LG eliminated: %lld, RG eliminated: %lld, sliced.size()=%lld, X_LG_RG_for_LG.size()=%lld, LG_morton=%lld, RG_morton=%lld, n_RG=%lld, n_LG=%lld\n", 
                                    //     LG_eliminated,
                                    //     RG_eliminated,
                                    //     sliced.size(),
                                    //     X_LG_RG_for_LG.size(),
                                    //     LG_morton,
                                    //     RG_morton,
                                    //     n_RG,
                                    //     n_LG);
                                    throw std::runtime_error("X_NN size mismatch for LG");
                                }
                                
                                for (int64_t i = 0; i < sliced.size(); ++i) {
                                    sliced[i] += X_LG_RG_for_LG[i];
                                }
                                
                                target_matrix.set_owned(
                                    n_RG, n_LG, std::move(sliced), MatrixStorage<DataType>::FULL);
                                
                            } else {
                                // Block doesn't exist - evaluate A_NN and add update
                                
                                
                                // Get LG coordinates
                                std::vector<CoordType> LG_coords;
                                bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                
                                if (LG_eliminated) {
                                    LG_coords.resize(n_LG * dimension);
                                    for (int64_t i = 0; i < n_LG; ++i) {
                                        int64_t src_idx = LG_box->skeleton_indices[i];
                                        for (int d = 0; d < dimension; ++d) {
                                            LG_coords[i * dimension + d] = 
                                                LG_box->point_coords[src_idx * dimension + d];
                                        }
                                    }
                                } else {
                                    LG_coords = LG_box->point_coords;
                                }
                                
                                // Get RG coordinates
                                std::vector<CoordType> RG_coords_temp;
                                const CoordType* RG_coords_ptr =
                                    coords_ptr_maybe_sliced(level,
                                                            RG_morton,
                                                            n_RG,
                                                            dimension,
                                                            /*box_or_null=*/ RG_is_assisting ? nullptr : RG_box,
                                                            /*assist_or_null=*/ RG_is_assisting ? RG_assisting : nullptr,
                                                            RG_coords_temp);
                                
                                // Kernel evaluate
                                
                                std::vector<DataType> A_LG_RG(n_LG * n_RG);
                                kernel->evaluate_block(
                                    LG_coords.data(), n_LG,
                                    RG_coords_ptr, n_RG,
                                    A_LG_RG.data(), n_LG
                                );
                                
                                // Transpose and add update
                                std::vector<DataType> A_transposed(n_LG * n_RG);
                                for (int64_t i = 0; i < n_LG; ++i) {
                                    for (int64_t j = 0; j < n_RG; ++j) {
                                        A_transposed[j + i * n_RG] = A_LG_RG[i + j * n_LG] + X_LG_RG_for_LG[j + i * n_RG];
                                    }
                                }
                                // if(box->morton_index == 64 && LG_morton == 65 && RG_morton == 21){
                                //     printf("LG morton: %ld, RG morton: %ld, LG assist: %d, RG assist: %d\n", LG_morton, RG_morton, LG_is_assisting, RG_is_assisting);
                                //     print_workspace_stats(A_transposed, "relevant (after gather_id_workspace)");
                                //     fflush(stdout);

                                // }
                                ModifiedBlock<DataType> new_block;
                                new_block.neighbor_morton = RG_morton;
                                auto& target_matrix = (is_symmetric || is_hermitian) ? new_block.A_NS : new_block.A_SN;
                                target_matrix.set_owned(
                                    n_RG, n_LG, std::move(A_transposed), MatrixStorage<DataType>::FULL);
                                
                                int64_t new_idx = modified_interactions.size();
                                modified_interactions.push_back(std::move(new_block));
                                interaction_map[RG_morton] = new_idx;
                            }
                        }
                        omp_unset_lock(box_locks.at(LG_morton));
                    }
                    
                    // ===== Update RG's view of LG (owned only) =====
                    if (RG_owned_here && !is_diagonal) {
                        omp_set_lock(box_locks.at(RG_morton));
                        if (is_symmetric || is_hermitian) {
                            // Store X_LG_RG in RG's A_NS
                            auto& interaction_map = is_one_hop ?
                                RG_box->near_field_interaction_map :
                                RG_box->far_field_interaction_map;
                            
                            auto& modified_interactions = is_one_hop ?
                                RG_box->near_field_modified_interactions :
                                RG_box->far_field_modified_interactions;
                            
                            auto it_RG = interaction_map.find(LG_morton);
                            
                            if (it_RG != interaction_map.end()) {
                                int64_t block_idx = it_RG->second;
                                auto& block = modified_interactions[block_idx];
                                
                                // Slice if needed
                                bool RG_eliminated = (level.eliminated_boxes.find(RG_morton) != level.eliminated_boxes.end());
                                bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                
                                std::vector<DataType> sliced;

                                if (RG_eliminated && LG_eliminated) {
                                    sliced = slice_modified_block_both_directions<CoordType, DataType>(
                                        block, level, LG_morton,
                                        RG_box->skeleton_indices,
                                        true,
                                        n_LG
                                    );
                                } else if (RG_eliminated) {
                                    int64_t stored_rows = block.A_NS.rows;
                                    int64_t stored_cols = block.A_NS.cols;
                                    
                                    if (stored_rows == n_LG && stored_cols == n_RG) {
                                        sliced = block.A_NS.data;
                                    } else if (stored_rows == n_LG && stored_cols > n_RG) {
                                        sliced.resize(n_LG * n_RG);
                                        for (int64_t j = 0; j < n_RG; ++j) {
                                            int64_t src_col = RG_box->skeleton_indices[j];
                                            for (int64_t i = 0; i < n_LG; ++i) {
                                                sliced[i + j * n_LG] = block.A_NS(i, src_col);
                                            }
                                        }
                                    } else {
                                        throw std::runtime_error("Step 9 RG update (sym): unexpected dimensions");
                                    }
                                } else if (LG_eliminated) {
                                    int64_t stored_rows = block.A_NS.rows;
                                    int64_t stored_cols = block.A_NS.cols;
                                    
                                    if (stored_rows == n_LG && stored_cols == n_RG) {
                                        sliced = block.A_NS.data;
                                    } else if (stored_rows > n_LG && stored_cols == n_RG) {
                                        sliced.resize(n_LG * n_RG);

                                        // Pick the skeleton index source once (local/ghost LG_box vs assisting LG box)
                                        const std::vector<int64_t>* lg_skel = nullptr;
                                        if (LG_box != nullptr) {
                                            lg_skel = &LG_box->skeleton_indices;
                                        } else {
                                            auto it = level.assisting_box_points_for_kernel_evaluation.find(LG_morton);
                                            if (it == level.assisting_box_points_for_kernel_evaluation.end()) {
                                                throw std::runtime_error("Missing assisting skeleton indices for LG_morton=" +
                                                                        std::to_string(LG_morton));
                                            }
                                            const auto idx = static_cast<size_t>(it->second);
                                            if (idx >= level.assisting_boxes.size()) {
                                                throw std::runtime_error("Assisting index out of range for LG_morton=" +
                                                                        std::to_string(LG_morton));
                                            }
                                            lg_skel = &level.assisting_boxes[idx].skel_indices;
                                        }

                                        if (static_cast<int64_t>(lg_skel->size()) < n_LG) {
                                            throw std::runtime_error("LG skeleton size < n_LG for LG_morton=" +
                                                                    std::to_string(LG_morton));
                                        }

                                        const auto& sk = *lg_skel;

                                        for (int64_t j = 0; j < n_RG; ++j) {
                                            for (int64_t i = 0; i < n_LG; ++i) {
                                                const int64_t src_row = sk[static_cast<size_t>(i)];
                                                sliced[static_cast<size_t>(i + j * n_LG)] = block.A_NS(src_row, j);
                                            }
                                        }
                                    }else {
                                        throw std::runtime_error("Step 9 RG view sym LG_elim: unexpected dims");
                                    }
                                } else {
                                    // Neither eliminated (or LG is assisting)
                                    sliced = block.A_NS.data;
                                }
                                
                                if (sliced.size() != X_LG_RG.size()) {
                                    throw std::runtime_error("X_NN size mismatch for RG");
                                }
                                
                                for (int64_t i = 0; i < sliced.size(); ++i) {
                                    sliced[i] += X_LG_RG[i];
                                }
                                
                                block.A_NS.set_owned(
                                    n_LG, n_RG, std::move(sliced), MatrixStorage<DataType>::FULL);
                                
                            } else {
                                // Block doesn't exist
                                
                                // Get RG coordinates
                                std::vector<CoordType> RG_coords;
                                bool RG_eliminated = (level.eliminated_boxes.find(RG_morton) != level.eliminated_boxes.end());
                                
                                if (RG_eliminated) {
                                    RG_coords.resize(n_RG * dimension);
                                    for (int64_t i = 0; i < n_RG; ++i) {
                                        int64_t src_idx = RG_box->skeleton_indices[i];
                                        for (int d = 0; d < dimension; ++d) {
                                            RG_coords[i * dimension + d] = 
                                                RG_box->point_coords[src_idx * dimension + d];
                                        }
                                    }
                                } else {
                                    RG_coords = RG_box->point_coords;
                                }
                                
                                // Get LG coordinates
                                std::vector<CoordType> LG_coords_temp;
                                const CoordType* LG_coords_ptr =
                                    coords_ptr_maybe_sliced(level,
                                                            LG_morton,
                                                            n_LG,
                                                            dimension,
                                                            /*box_or_null=*/ LG_is_assisting ? nullptr : LG_box,
                                                            /*assist_or_null=*/ LG_is_assisting ? LG_assisting : nullptr,
                                                            LG_coords_temp);
                                
                                // Kernel evaluate
                                std::vector<DataType> A_LG_RG(n_LG * n_RG);
                                kernel->evaluate_block(
                                    LG_coords_ptr, n_LG,
                                    RG_coords.data(), n_RG,
                                    A_LG_RG.data(), n_LG
                                );
                                
                                // Add update
                                for (int64_t i = 0; i < n_LG * n_RG; ++i) {
                                    A_LG_RG[i] += X_LG_RG[i];
                                }
                                
                                ModifiedBlock<DataType> new_block;
                                new_block.neighbor_morton = LG_morton;
                                new_block.A_NS.set_owned(
                                    n_LG, n_RG, std::move(A_LG_RG), MatrixStorage<DataType>::FULL);
                                
                                int64_t new_idx = modified_interactions.size();
                                modified_interactions.push_back(std::move(new_block));
                                interaction_map[LG_morton] = new_idx;
                            }
                            
                        } else {
                            std::runtime_error("need to fix assist by using proper slicing");
                            // Nonsymmetric: Store in RG's A_NS
                            auto& interaction_map = is_one_hop ?
                                RG_box->near_field_interaction_map_nonsymmetry :
                                RG_box->far_field_interaction_map_nonsymmetry;
                            
                            auto& modified_interactions = is_one_hop ?
                                RG_box->near_field_modified_interactions :
                                RG_box->far_field_modified_interactions;
                            
                            auto it_RG = interaction_map.find(LG_morton);
                            
                            if (it_RG != interaction_map.end()) {
                                int64_t block_idx = it_RG->second;
                                auto& block = modified_interactions[block_idx];
                                
                                // Slice if needed
                                bool RG_eliminated = (level.eliminated_boxes.find(RG_morton) != level.eliminated_boxes.end());
                                bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                
                                std::vector<DataType> sliced;

                                if (RG_eliminated && LG_eliminated) {
                                    sliced = slice_modified_block_both_directions<CoordType, DataType>(
                                        block, level, LG_morton,
                                        RG_box->skeleton_indices,
                                        false,
                                        n_LG
                                    );
                                } else if (RG_eliminated) {
                                    int64_t stored_rows = block.A_NS.rows;
                                    int64_t stored_cols = block.A_NS.cols;
                                    
                                    if (stored_rows == n_LG && stored_cols == n_RG) {
                                        sliced = block.A_NS.data;
                                    } else if (stored_rows == n_LG && stored_cols > n_RG) {
                                        sliced.resize(n_LG * n_RG);
                                        for (int64_t j = 0; j < n_RG; ++j) {
                                            int64_t src_col = RG_box->skeleton_indices[j];
                                            for (int64_t i = 0; i < n_LG; ++i) {
                                                sliced[i + j * n_LG] = block.A_NS(i, src_col);
                                            }
                                        }
                                    } else {
                                        throw std::runtime_error("Step 9 RG view nonsym RG_elim: unexpected dims");
                                    }
                                } else if (LG_eliminated) {
                                    int64_t stored_rows = block.A_NS.rows;
                                    int64_t stored_cols = block.A_NS.cols;
                                    
                                    if (stored_rows == n_LG && stored_cols == n_RG) {
                                        sliced = block.A_NS.data;
                                    } else if (stored_rows > n_LG && stored_cols == n_RG) {
                                        sliced.resize(n_LG * n_RG);
                                        for (int64_t j = 0; j < n_RG; ++j) {
                                            for (int64_t i = 0; i < n_LG; ++i) {
                                                int64_t src_row = LG_box->skeleton_indices[i];
                                                sliced[i + j * n_LG] = block.A_NS(src_row, j);
                                            }
                                        }
                                    } else {
                                        throw std::runtime_error("Step 9 RG view nonsym LG_elim: unexpected dims");
                                    }
                                } else {
                                    // Neither eliminated (or LG is assisting)
                                    sliced = block.A_NS.data;
                                }
                                
                                if (sliced.size() != X_LG_RG.size()) {
                                    throw std::runtime_error("X_NN size mismatch for RG A_NS");
                                }
                                
                                for (int64_t i = 0; i < sliced.size(); ++i) {
                                    sliced[i] += X_LG_RG[i];
                                }
                                
                                block.A_NS.set_owned(
                                    n_LG, n_RG, std::move(sliced), MatrixStorage<DataType>::FULL);
                            } else {
                                // Block doesn't exist
                                
                                // Get RG coordinates
                                std::vector<CoordType> RG_coords;
                                bool RG_eliminated = (level.eliminated_boxes.find(RG_morton) != level.eliminated_boxes.end());
                                
                                if (RG_eliminated) {
                                    RG_coords.resize(n_RG * dimension);
                                    for (int64_t i = 0; i < n_RG; ++i) {
                                        int64_t src_idx = RG_box->skeleton_indices[i];
                                        for (int d = 0; d < dimension; ++d) {
                                            RG_coords[i * dimension + d] = 
                                                RG_box->point_coords[src_idx * dimension + d];
                                        }
                                    }
                                } else {
                                    RG_coords = RG_box->point_coords;
                                }
                                
                                // Get LG coordinates
                                const CoordType* LG_coords_ptr;
                                std::vector<CoordType> LG_coords_temp;
                                
                                if (LG_is_assisting) {
                                    LG_coords_ptr = LG_assisting->coords.data();
                                } else {
                                    bool LG_eliminated = (level.eliminated_boxes.find(LG_morton) != level.eliminated_boxes.end());
                                    
                                    if (LG_eliminated) {
                                        LG_coords_temp.resize(n_LG * dimension);
                                        for (int64_t i = 0; i < n_LG; ++i) {
                                            int64_t src_idx = LG_box->skeleton_indices[i];
                                            for (int d = 0; d < dimension; ++d) {
                                                LG_coords_temp[i * dimension + d] = 
                                                    LG_box->point_coords[src_idx * dimension + d];
                                            }
                                        }
                                        LG_coords_ptr = LG_coords_temp.data();
                                    } else {
                                        LG_coords_ptr = LG_box->point_coords.data();
                                    }
                                }
                                
                                // Kernel evaluate
                                std::vector<DataType> A_LG_RG(n_LG * n_RG);
                                kernel->evaluate_block(
                                    LG_coords_ptr, n_LG,
                                    RG_coords.data(), n_RG,
                                    A_LG_RG.data(), n_LG
                                );
                                
                                // Add update
                                for (int64_t i = 0; i < n_LG * n_RG; ++i) {
                                    A_LG_RG[i] += X_LG_RG[i];
                                }
                                
                                ModifiedBlock<DataType> new_block;
                                new_block.neighbor_morton = LG_morton;
                                new_block.A_NS.set_owned(
                                    n_LG, n_RG, std::move(A_LG_RG), MatrixStorage<DataType>::FULL);
                                
                                int64_t new_idx = modified_interactions.size();
                                modified_interactions.push_back(std::move(new_block));
                                interaction_map[LG_morton] = new_idx;
                            }
                        }
                        omp_unset_lock(box_locks.at(RG_morton)); 
                    }
                    
                    row_offset += n_LG;
                }
                
                col_offset += n_RG;
            }
        }
    }

    if (store) {
        // Keep the original X_NR plus temp2 until the deferred post-wave stage
        // reconstructs Step 7-9 for all three ownership cases:
        //   1. local-local   -> local owner replay + local mirror copy
        //   2. local-remote  -> local replay + transported canonical ADD/REPLACE
        //   3. remote-remote -> transported canonical ADD only
        if (total_neighbor_points > 0 && !temp2.empty()) {
            box->deferred_xnn_temp2 = std::move(temp2);
            box->deferred_xnn_neighbor_point_counts = neighbor_point_counts;
        }
    }

     // ===== Step 9.5: Slice far-field blocks in both directions =====
    // (a) Slice current box B's far-field blocks (S dimension: B's skeleton)
    // (b) Slice neighbors' reciprocal far-field blocks (N dimension: B's skeleton)
    
    // Part (a): Slice B's own far-field blocks
    // for (auto& modified_block : box->far_field_modified_interactions) {
    //     // Slice A_NS: (n_neighbor × num_points) → (n_neighbor × k)
    //     if (modified_block.A_NS.is_allocated()) {
    //         int64_t stored_rows = modified_block.A_NS.rows;
    //         int64_t stored_cols = modified_block.A_NS.cols;
            
    //         if (stored_cols != k) {
    //             std::vector<DataType> sliced(stored_rows * k);
    //             std::vector<int64_t> all_rows(stored_rows);
    //             std::iota(all_rows.begin(), all_rows.end(), 0);
                
    //             extract_submatrix(
    //                 modified_block.A_NS.data.data(), stored_rows,
    //                 all_rows, box->skeleton_indices,
    //                 sliced.data(), stored_rows
    //             );
                
    //             modified_block.A_NS.data = std::move(sliced);
    //             modified_block.A_NS.cols = k;
    //         }
    //     }
        
    //     // Slice A_SN: (num_points × n_neighbor) → (k × n_neighbor)
    //     if (!is_symmetric && !is_hermitian && modified_block.A_SN.is_allocated()) {
    //         int64_t stored_rows = modified_block.A_SN.rows;
    //         int64_t stored_cols = modified_block.A_SN.cols;
            
    //         if (stored_rows != k) {
    //             std::vector<DataType> sliced(k * stored_cols);
    //             std::vector<int64_t> all_cols(stored_cols);
    //             std::iota(all_cols.begin(), all_cols.end(), 0);
                
    //             extract_submatrix(
    //                 modified_block.A_SN.data.data(), stored_rows,
    //                 box->skeleton_indices, all_cols,
    //                 sliced.data(), k
    //             );
                
    //             modified_block.A_SN.data = std::move(sliced);
    //             modified_block.A_SN.rows = k;
    //             modified_block.A_SN.lda = k;
    //         }
    //     }
    // }
    
    // // Part (b): Slice reciprocal views in neighbors' far-field blocks
    // for (const auto& modified_block : box->far_field_modified_interactions) {
    //     int64_t neighbor_morton = modified_block.neighbor_morton;
        
    //     BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
    //     if (neighbor_box == nullptr) {
    //         continue; // do not modify ghosts; remote owner will handle slicing
    //     }
        
    //     // Check if G has B in its far-field map
    //     auto& neighbor_far_map = is_symmetric || is_hermitian ? 
    //         neighbor_box->far_field_interaction_map : 
    //         neighbor_box->far_field_interaction_map_nonsymmetry;
        
    //     auto it = neighbor_far_map.find(box->morton_index);
    //     if (it == neighbor_far_map.end()) {
    //         continue; // G doesn't have B in far-field, skip
    //     }
        
    //     // Found G's far-field block for B - slice N dimension (B's dimension)
    //     int64_t block_idx = it->second;
    //     auto& neighbor_block = neighbor_box->far_field_modified_interactions[block_idx];
        
    //     // Slice A_NS rows: (n_B × k_G) → (k_B × k_G)
    //     if (neighbor_block.A_NS.is_allocated()) {
    //         int64_t stored_rows = neighbor_block.A_NS.rows;
    //         int64_t stored_cols = neighbor_block.A_NS.cols;
            
    //         if (stored_rows != k) {
    //             std::vector<DataType> sliced(k * stored_cols);
    //             std::vector<int64_t> all_cols(stored_cols);
    //             std::iota(all_cols.begin(), all_cols.end(), 0);
                
    //             extract_submatrix(
    //                 neighbor_block.A_NS.data.data(), stored_rows,
    //                 box->skeleton_indices, all_cols,
    //                 sliced.data(), k
    //             );
                
    //             neighbor_block.A_NS.data = std::move(sliced);
    //             neighbor_block.A_NS.rows = k;
    //             neighbor_block.A_NS.lda = k;
    //         }
    //     }
        
    //     // Nonsymmetric: also slice A_SN columns: (k_G × n_B) → (k_G × k_B)
    //     if (!is_symmetric && !is_hermitian && neighbor_block.A_SN.is_allocated()) {
    //         int64_t stored_rows = neighbor_block.A_SN.rows;
    //         int64_t stored_cols = neighbor_block.A_SN.cols;
            
    //         if (stored_cols != k) {
    //             std::vector<DataType> sliced(stored_rows * k);
    //             std::vector<int64_t> all_rows(stored_rows);
    //             std::iota(all_rows.begin(), all_rows.end(), 0);
                
    //             extract_submatrix(
    //                 neighbor_block.A_SN.data.data(), stored_rows,
    //                 all_rows, box->skeleton_indices,
    //                 sliced.data(), stored_rows
    //             );
                
    //             neighbor_block.A_SN.data = std::move(sliced);
    //             neighbor_block.A_SN.cols = k;
    //         }
    //     }
    // }
    
    // ===== Step 10: Store elimination matrices for solve phase =====
    
    box->X_SR.data.assign(temp1.begin(), temp1.end());
    
    if (!temp2.empty() && !store) {
        box->X_NR.data.assign(temp2.begin(), temp2.end());
    }
    
    if (!is_symmetric && !is_hermitian) {
        box->X_RS.data.assign(temp3.begin(), temp3.end());
        if (!temp4.empty()) {
            box->X_RN.data.assign(temp4.begin(), temp4.end());
        }
    }
    

}




/**
 * @brief Perform complete two-step elimination (Sections 2.3-2.4)
 * 
 * Combines first and second elimination steps:
 * - Step 1: ID decomposition, compute X blocks
 * - Step 2: Apply L and U eliminations, compute Schur complement
 * 
 * @param box Current box being factorized
 * @param level Tree level containing neighbors (for workspace gathering)
 * @param kernel Kernel evaluator
 * @param workspace A_FB matrix workspace (will be modified)
 * @param workspace_rows Rows in workspace
 * @param workspace_cols Columns in workspace
 * @param tolerance ID tolerance
 * @param is_symmetric Whether matrix is symmetric
 * @param is_hermitian Whether matrix is Hermitian
 * @param factorization_method How to factorize X_RR
 */
template<typename CoordType, typename DataType, typename KernelType>
void compute_and_modify(
    int dimension,
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    FactorizationThreadScratch<CoordType, DataType>& scratch,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    PendingFactorUpdates<DataType> *pending_updates,
    FactorizationMethod factorization_method, bool store = false,
    int DEBUG = 0, int VERBOSE = 0) {
    
    static_assert(std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>,
                  "Only double precision supported currently");
    auto& workspace = scratch.workspace;
    const int64_t workspace_rows = scratch.workspace_rows;
    const int64_t workspace_cols = scratch.workspace_cols;
    auto& X_NN_full = scratch.x_nn_full;
    auto& sketch_storage = scratch.sketch_storage;
    
    if (box->num_points == 0) {
        workspace.clear();
        // level.eliminated_boxes.insert(box->morton_index);
        return;
    }

    // DEBUG: Show existing blocks at start
    if (DEBUG){
        std::cout << "\n=== START compute_and_modify for Box " << box->morton_index << " ===" << std::endl;
        std::cout << "Existing near_field_modified_interactions: " << box->near_field_modified_interactions.size() << std::endl;
        for (const auto& [morton, idx] : box->near_field_interaction_map) {
            const auto& block = box->near_field_modified_interactions[idx];
            std::cout << "  Neighbor " << morton << ": A_NS = ";
            if (block.A_NS.is_allocated()) {
                std::cout << block.A_NS.rows << " × " << block.A_NS.cols;
            } else {
                std::cout << "not allocated";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Box " << box->morton_index << "'s 1-hop neighbors: ";
        for (auto n : box->one_hop) std::cout << n << " ";
        std::cout << std::endl;
        // ADD THIS:
        std::cout << "near_field_interaction_map contents:" << std::endl;
        for (const auto& [morton, idx] : box->near_field_interaction_map) {
            std::cout << "  morton " << morton << " -> idx " << idx;
            if (idx < box->near_field_modified_interactions.size()) {
                const auto& block = box->near_field_modified_interactions[idx];
                std::cout << " (block " << (block.A_NS.is_allocated() ? "allocated" : "NOT allocated") << ")";
            } else {
                std::cout << " (INDEX OUT OF BOUNDS!)";
            }
            std::cout << std::endl;
        }

    }
   

    


    
    
    // ========================================================================
    // STEP 1: First Elimination (Section 2.3)
    // ========================================================================
    
    // compute ID
    
    // const bool force_id_diagnostics = (box->morton_index == 238736);
    const double sketch_factor = 1.0;
    const int sketch_nonzeros = 4;
    int64_t sketch_rows = static_cast<int64_t>(
        std::ceil(sketch_factor * static_cast<double>(workspace_cols)));
    sketch_rows = std::min(sketch_rows, workspace_rows);
    sketch_rows = std::max(sketch_rows, workspace_cols);
    // if (force_id_diagnostics) {
    //     const auto workspace_summary =
    //         summarize_colmajor_matrix(workspace, workspace_rows, workspace_cols);
    //     std::cerr << "[ID diagnostics] box=" << box->morton_index
    //               << " (pre-RRQR)"
    //               << " num_points=" << box->num_points
    //               << " one_hop=" << box->one_hop.size()
    //               << " two_hop=" << box->two_hop.size()
    //               << " far_field_map=" << box->far_field_interaction_map.size()
    //               << " far_field_blocks=" << box->far_field_modified_interactions.size()
    //               << " near_field_map=" << box->near_field_interaction_map.size()
    //               << " near_field_blocks=" << box->near_field_modified_interactions.size()
    //               << '\n'
    //               << "  "
    //               << format_matrix_diagnostic_summary("workspace", workspace_summary)
    //               << '\n';
    // }
    // auto id_start = std::chrono::high_resolution_clock::now();
    // auto id_result = fmm::compute_id(
    //     workspace.data(),
    //     workspace_rows, workspace_cols, workspace_rows,
    //     tolerance, 0
    // );
    // auto id_result = fmm::compute_id_uniform_sketch(
    //     workspace.data(), sketch_storage,
    //     workspace_rows, workspace_cols, workspace_rows,
    //     tolerance, 1.0, box->morton_index + 1
    // );
    IDResult<DataType> id_result;
    try {
        id_result = fmm::compute_id_sparse_sketch(
            workspace.data(), sketch_storage,
            workspace_rows, workspace_cols, workspace_rows,
            tolerance, sketch_factor, sketch_nonzeros, box->morton_index + 1
        );
    } catch (const std::exception& e) {
        const auto workspace_summary =
            summarize_colmajor_matrix(workspace, workspace_rows, workspace_cols);
        const auto sketch_summary =
            summarize_colmajor_matrix(sketch_storage, sketch_rows, workspace_cols);
        std::ostringstream oss;
        oss << "compute_and_modify: ID failure for box " << box->morton_index
            << " (num_points=" << box->num_points
            << ", one_hop=" << box->one_hop.size()
            << ", two_hop=" << box->two_hop.size()
            << ", far_field_map=" << box->far_field_interaction_map.size()
            << ", far_field_blocks=" << box->far_field_modified_interactions.size()
            << ", near_field_map=" << box->near_field_interaction_map.size()
            << ", near_field_blocks=" << box->near_field_modified_interactions.size()
            << ", workspace_rows=" << workspace_rows
            << ", workspace_cols=" << workspace_cols
            << ", sketch_rows=" << sketch_rows
            << ", sketch_cols=" << workspace_cols
            << ", sketch_factor=" << sketch_factor
            << ", sketch_nonzeros=" << sketch_nonzeros
            << ")\n  "
            << format_matrix_diagnostic_summary("workspace", workspace_summary)
            << "\n  "
            << format_matrix_diagnostic_summary("sketch", sketch_summary)
            << "\n  underlying_error=" << e.what();
        throw std::runtime_error(oss.str());
    }
    // if (force_id_diagnostics) {
    //     const auto sketch_summary =
    //         summarize_colmajor_matrix(sketch_storage, sketch_rows, workspace_cols);
    //     std::cerr << "  "
    //               << format_matrix_diagnostic_summary("sketch", sketch_summary)
    //               << '\n';
    // }
    // auto id_result = fmm::compute_id_with_rrqr(
    //     workspace.data(),
    //     workspace_rows, workspace_cols, workspace_rows,
    //     tolerance, 1.25
    // );
    // auto id_end = std::chrono::high_resolution_clock::now();
    // auto id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
    // std::cout << "  id time: " << id_duration.count() << " ms" << std::endl;
    
    box->skeleton_indices = id_result.skeleton_indices;
    box->redundant_indices = id_result.redundant_indices;
    box->interpolation_matrix = std::move(id_result.interpolation);
    
    int64_t k = box->skeleton_indices.size();
    int64_t r = box->redundant_indices.size();
    
    if (r == 0) {
        workspace.clear();
        // level.eliminated_boxes.insert(box->morton_index);
        return;  // Full rank, nothing to compress
    }

    // **ADD DEBUG HERE** ↓↓↓
    if (DEBUG) {
        std::cout << "DEBUG: After ID decomposition for Box " << box->morton_index << std::endl;
        std::cout << "  Skeleton DOFs (k): " << k << std::endl;
        std::cout << "  Redundant DOFs (r): " << r << std::endl;
        
        // Check interpolation matrix (T)
        bool T_has_nan = false;
        if (box->interpolation_matrix.is_allocated()) {
            for (size_t i = 0; i < box->interpolation_matrix.data.size(); i++) {
                if (is_nan(box->interpolation_matrix.data[i])) {
                    T_has_nan = true;
                    break;
                }
            }
        }
        std::cout << "  Interpolation matrix T contains NaN: " << (T_has_nan ? "YES" : "NO") << std::endl;
        std::cout << "  T allocated: " << (box->interpolation_matrix.is_allocated() ? "yes" : "no") << std::endl;
        if (box->interpolation_matrix.is_allocated()) {
            std::cout << "  T size: " << box->interpolation_matrix.data.size() 
                    << " (expected: " << (k * r) << ")" << std::endl;
        }
        
        // Print first few skeleton and redundant indices
        std::cout << "  First 5 skeleton indices: ";
        for (size_t i = 0; i < std::min(size_t(5), box->skeleton_indices.size()); i++) {
            std::cout << box->skeleton_indices[i] << " ";
        }
        std::cout << std::endl;
        
        std::cout << "  First 5 redundant indices: ";
        for (size_t i = 0; i < std::min(size_t(5), box->redundant_indices.size()); i++) {
            std::cout << box->redundant_indices[i] << " ";
        }
        std::cout << std::endl;
    }


    
    // Form X_BB (self-interaction)
    const std::vector<DataType>* X_BB_ptr = nullptr;

    if (box->schur_complement.is_allocated()) {
        // Box has been modified by previous eliminations - use updated schur complement
        if (box->schur_complement.rows != box->num_points || 
            box->schur_complement.cols != box->num_points) {
            throw std::runtime_error(
                "compute_and_modify: Schur complement size mismatch. Expected " +
                std::to_string(box->num_points) + "×" + std::to_string(box->num_points) +
                ", got " + std::to_string(box->schur_complement.rows) + "×" + 
                std::to_string(box->schur_complement.cols));
        }
        
        // std::cout << "Using updated schur_complement for X_BB (size: " 
        //         << box->schur_complement.rows << " × " << box->schur_complement.cols << ")" << std::endl;
        
        X_BB_ptr = &box->schur_complement.data;
    } else {
        // First time processing this box - evaluate kernel
        // std::cout << "Evaluating fresh kernel for X_BB" << std::endl;
        
        scratch.x_bb.resize(static_cast<size_t>(box->num_points * box->num_points));
        kernel->evaluate_block(
            box->point_coords.data(), box->num_points,
            box->point_coords.data(), box->num_points,
            scratch.x_bb.data(), box->num_points
        );
        X_BB_ptr = &scratch.x_bb;
    }
    const auto& X_BB = *X_BB_ptr;

    if (DEBUG) {
        bool x_bb_has_nan = false;
        for (auto val : X_BB) {
            if (is_nan(val)) {
                x_bb_has_nan = true;
                break;
            }
        }
        std::cout << "DEBUG: After getting X_BB for Box " << box->morton_index << std::endl;
        std::cout << "  X_BB contains NaN: " << (x_bb_has_nan ? "YES" : "NO") << std::endl;
        std::cout << "  X_BB came from: " << (box->schur_complement.is_allocated() ? "schur_complement" : "kernel eval") << std::endl;
    }
    
    // Extract submatrices from X_BB
    std::vector<DataType> A_RR(r * r);
    extract_submatrix(X_BB.data(), box->num_points,
                     box->redundant_indices, box->redundant_indices,
                     A_RR.data(), r);
    
    // Extract A_SR (always needed)
    std::vector<DataType> A_SR(k * r);
    extract_submatrix(X_BB.data(), box->num_points,
                     box->skeleton_indices, box->redundant_indices,
                     A_SR.data(), k);
    
    // Extract A_RS (nonsymmetric only)
    std::vector<DataType> A_RS;
    if (!is_symmetric && !is_hermitian) {
        A_RS.resize(r * k);
        extract_submatrix(X_BB.data(), box->num_points,
                         box->redundant_indices, box->skeleton_indices,
                         A_RS.data(), r);
    }


    
    std::vector<DataType> A_SS(k * k);
    extract_submatrix(X_BB.data(), box->num_points,
                     box->skeleton_indices, box->skeleton_indices,
                     A_SS.data(), k);
    
    const auto& T = box->interpolation_matrix;
    
    // Compute X_RR, X_RS, X_SR using BLAS
    box->X_RR.set_owned(r, r, std::move(A_RR), MatrixStorage<DataType>::FULL);
    
    if (is_symmetric || is_hermitian) {
        std::vector<DataType> temp1(r * r, 0.0);
        int M = r, N = r, K = k;
        DataType alpha = -1.0, beta = 0.0;
        
        // Compute temp1 = -A_SR^T * T (use transpose since A_RS not allocated)
        gemm_("T", "N", &M, &N, &K,
               &alpha, A_SR.data(), &K,  // A_SR is (k × r), transposed to (r × k)
               T.data.data(), &K,
               &beta, temp1.data(), &M);
        
        for (int64_t j = 0; j < r; ++j) {
            for (int64_t i = 0; i < r; ++i) {
                box->X_RR(i, j) += temp1[i + j * r] + temp1[j + i * r];
            }
        }
        
        std::vector<DataType> temp2(k * r, 0.0);
        M = k; N = r; K = k;
        alpha = 1.0; beta = 0.0;
        
        gemm_("N", "N", &M, &N, &K,
               &alpha, A_SS.data(), &M,
               T.data.data(), &K,
               &beta, temp2.data(), &M);
        
        M = r; N = r; K = k;
        alpha = 1.0; beta = 1.0;
        
        gemm_("T", "N", &M, &N, &K,
               &alpha, T.data.data(), &K,
               temp2.data(), &K,
               &beta, box->X_RR.data.data(), &M);
        
    } else {
        int M = r, N = r, K = k;
        DataType alpha = -1.0, beta = 1.0;
        
        gemm_("T", "N", &M, &N, &K,
               &alpha, T.data.data(), &K,
               A_SR.data(), &K,
               &beta, box->X_RR.data.data(), &M);
        
        gemm_("N", "N", &M, &N, &K,
               &alpha, A_RS.data(), &M,
               T.data.data(), &K,
               &beta, box->X_RR.data.data(), &M);
        
        std::vector<DataType> temp(k * r, 0.0);
        M = k; N = r; K = k;
        alpha = 1.0; beta = 0.0;
        
        gemm_("N", "N", &M, &N, &K,
               &alpha, A_SS.data(), &M,
               T.data.data(), &K,
               &beta, temp.data(), &M);
        
        M = r; N = r; K = k;
        alpha = 1.0; beta = 1.0;
        
        gemm_("T", "N", &M, &N, &K,
               &alpha, T.data.data(), &K,
               temp.data(), &K,
               &beta, box->X_RR.data.data(), &M);
    }
    
    // Compute X_SR (always needed)
    box->X_SR.set_owned(k, r, std::move(A_SR), MatrixStorage<DataType>::FULL);
    
    int M = k, N = r, K = k;
    DataType alpha = -1.0, beta = 1.0;
    
    gemm_("N", "N", &M, &N, &K,
           &alpha, A_SS.data(), &M,
           T.data.data(), &K,
           &beta, box->X_SR.data.data(), &M);
    
    // Compute X_RS (nonsymmetric only)
    if (!is_symmetric && !is_hermitian) {
        box->X_RS.set_owned(r, k, std::move(A_RS), MatrixStorage<DataType>::FULL);
        
        M = r; N = k; K = k;
        alpha = -1.0; beta = 1.0;
        
        gemm_("T", "N", &M, &N, &K,
               &alpha, T.data.data(), &K,
               A_SS.data(), &K,
               &beta, box->X_RS.data.data(), &M);
    }
    
    // Factorize X_RR
    if (factorization_method == FactorizationMethod::CHOLESKY) {
        box->X_RR_pivots.clear();

        if (DEBUG){
            // Compute Frobenius norm
            DataType frob_norm = 0.0;
            for (int64_t j = 0; j < r; ++j) {
                for (int64_t i = 0; i < r; ++i) {
                    DataType val = box->X_RR(i, j);
                    frob_norm += std::norm(val);
                }
            }
            
            
            // Check symmetry
            DataType sym_error = 0.0;
            for (int64_t j = 0; j < r; ++j) {
                for (int64_t i = 0; i < j; ++i) {
                    DataType diff = box->X_RR(i, j) - box->X_RR(j, i);
                    sym_error += std::norm(diff);
                }
            }
        
            std::cout << "Frobenius norm: " << std::sqrt(frob_norm) << std::endl;
            std::cout << "Symmetry error: " << std::sqrt(sym_error) << std::endl;
        }
        

        char uplo = 'L';
        int n = r;
        int info = 0;

        if constexpr (std::is_same_v<DataType, double>) {
            // int mpi_rank;
            // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            // if(mpi_rank == 0){
                // auto rep = symmetry_error_colmajor_real(box->X_RR.data.data(),
                //                         box->X_RR.rows, box->X_RR.cols, box->X_RR.lda);
                // print_symmetry_report("X_RR", rep, std::cerr);
            // }
            
            dpotrf_(&uplo, &n, box->X_RR.data.data(), &n, &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zsychol_(&uplo, &n, box->X_RR.data.data(), &n, &info);
        }
        
        if (info != 0) {
            throw std::runtime_error(
                "Cholesky factorization failed with INFO = " + std::to_string(info));
        }
        
        box->X_RR.format = MatrixStorage<DataType>::CHOLESKY_L;
        
    } else if (factorization_method == FactorizationMethod::LU) {
        int m = r, n = r;
        box->X_RR_pivots.resize(static_cast<size_t>(r));
        int info = 0;

        if constexpr (std::is_same_v<DataType, double>) {
            dgetrf_(&m, &n, box->X_RR.data.data(), &m, box->X_RR_pivots.data(), &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgetrf_(&m, &n, box->X_RR.data.data(), &m, box->X_RR_pivots.data(), &info);
        }

        if (info != 0) {
            throw std::runtime_error(
                "LU factorization failed with INFO = " + std::to_string(info));
        }

        box->X_RR.format = MatrixStorage<DataType>::LU_FACTORED;

    } else if (factorization_method == FactorizationMethod::COMPLEX_SYM) {
        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            char uplo = 'L';
            int n = r;
            int lwork = -1;
            int info = 0;
            box->X_RR_pivots.resize(static_cast<size_t>(r));
            std::vector<DataType> work(1);
            zsytrf_(&uplo, &n, box->X_RR.data.data(), &n,
                    box->X_RR_pivots.data(), work.data(), &lwork, &info);
            lwork = static_cast<int>(work[0].real());
            work.resize(static_cast<size_t>(lwork));
            zsytrf_(&uplo, &n, box->X_RR.data.data(), &n,
                    box->X_RR_pivots.data(), work.data(), &lwork, &info);
            if (info != 0)
                throw std::runtime_error(
                    "Bunch-Kaufman factorization failed with INFO = " + std::to_string(info));
            box->X_RR.format = MatrixStorage<DataType>::BUNCH_KAUFMAN;
        } else {
            throw std::runtime_error("COMPLEX_SYM factorization only supported for complex<double>");
        }

    } else if (factorization_method == FactorizationMethod::NONE) {
        int n = r;
        std::vector<int> ipiv(r);
        int info = 0;
        box->X_RR_pivots.clear();

        if constexpr (std::is_same_v<DataType, double>) {
            dgetrf_(&n, &n, box->X_RR.data.data(), &n, ipiv.data(), &info);
            if (info != 0) {
                throw std::runtime_error("LU for inverse failed: " + std::to_string(info));
            }

            int lwork = -1;
            std::vector<double> work(1);
            dgetri_(&n, box->X_RR.data.data(), &n, ipiv.data(),
                    work.data(), &lwork, &info);

            lwork = static_cast<int>(work[0]);
            work.resize(lwork);

            dgetri_(&n, box->X_RR.data.data(), &n, ipiv.data(),
                    work.data(), &lwork, &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgetrf_(&n, &n, box->X_RR.data.data(), &n, ipiv.data(), &info);
            if (info != 0) {
                throw std::runtime_error("LU for inverse failed: " + std::to_string(info));
            }

            int lwork = -1;
            std::vector<std::complex<double>> work(1);
            zgetri_(&n, box->X_RR.data.data(), &n, ipiv.data(),
                    work.data(), &lwork, &info);

            lwork = static_cast<int>(work[0].real());
            work.resize(lwork);

            zgetri_(&n, box->X_RR.data.data(), &n, ipiv.data(),
                    work.data(), &lwork, &info);
        }

        if (info != 0) {
            throw std::runtime_error("Matrix inversion failed: " + std::to_string(info));
        }

        box->X_RR.format = MatrixStorage<DataType>::INVERSE;
    }
    
    // ===== Compute X_RN and X_NR from 1-hop near-field neighbors =====
    // Also preserve unmodified A_NS and A_SN for part 2
    
    int64_t total_neighbor_points = 0;
    auto& neighbor_point_counts = scratch.neighbor_point_counts;
    neighbor_point_counts.clear();
    neighbor_point_counts.reserve(box->one_hop.size());
    
    // First pass: count total neighbor points based on current elimination status
    for (int64_t neighbor_morton : box->one_hop) {
        BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
        
        if (neighbor_box == nullptr) {
            auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
            if (ghost_it != level.ghost_id_to_index.end()) {
                neighbor_box = &level.ghost_boxes[ghost_it->second];
            }
        }
        
        int64_t n_neighbor = 0;
        bool found_assist = false;
        if (neighbor_box != nullptr) {
            // Use current size based on elimination status
            if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                n_neighbor = neighbor_box->skeleton_indices.size();
            } else {
                n_neighbor = neighbor_box->num_points;
            }
        } else {
            // Check assisting boxes
            auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
            if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
                if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                    n_neighbor = level.assisting_boxes[assist_it->second].skel_indices.size();
                    found_assist = true;
                }
                else{
                    int64_t assist_idx = assist_it->second;
                    const auto& assist_box = level.assisting_boxes[assist_idx];
                    n_neighbor = assist_box.coords.size() / dimension;
                    found_assist = true;
                }
                
            }
        }
        
        if ((neighbor_box == nullptr && !found_assist) || n_neighbor == 0) {
            throw std::runtime_error(
                "compute_and_modify: 1-hop neighbor box " + 
                std::to_string(neighbor_morton) + " not found or has no points");
        }
        
        neighbor_point_counts.push_back(n_neighbor);
        total_neighbor_points += n_neighbor;
    }

    // store this information to be used for solve
    level.solve_neighbor_size[box->morton_index - level.local_morton_start] = neighbor_point_counts;
    
    // Unmodified A_NS and A_SN to pass to step 2
    auto& A_NS_all_unmodified = scratch.a_ns_all;
    auto& A_SN_all_unmodified = scratch.a_sn_all;
    A_NS_all_unmodified.clear();
    A_SN_all_unmodified.clear();
    
    if (total_neighbor_points > 0) {
        // Allocate batched arrays
        auto& A_NS_all = scratch.a_ns_all;
        A_NS_all.resize(static_cast<size_t>(total_neighbor_points * k));
        std::vector<DataType> A_NR_all(total_neighbor_points * r);
        
        int64_t current_row = 0;
        
        // Second pass: fill A_NS and A_NR
        for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
            int64_t neighbor_morton = box->one_hop[idx];
            int64_t n_neighbor = neighbor_point_counts[idx];
            
            auto it = box->near_field_interaction_map.find(neighbor_morton);
            
            if (it != box->near_field_interaction_map.end()) {
                int64_t block_idx = it->second;
                const auto& modified_block = box->near_field_modified_interactions[block_idx];

                if (DEBUG) {
                    bool stored_has_nan = false;
                    const auto& check_matrix = modified_block.A_NS;
                    for (auto val : check_matrix.data) {
                        if (is_nan(val)) {
                            stored_has_nan = true;
                            break;
                        }
                    }
                    
                    std::cout << "DEBUG Box " << box->morton_index << " neighbor " << neighbor_morton << std::endl;
                    std::cout << "  Stored block A_NS dims: " << check_matrix.rows << " × " << check_matrix.cols << std::endl;
                    std::cout << "  Stored block ALREADY has NaN: " << (stored_has_nan ? "YES" : "NO") << std::endl;
                    // if(neighbor_morton == 192)
                    // {
                    //     printf("DEBUG: neighbor_morton=%ld, block_idx=%ld, modified_block.A_NS.rows=%ld, modified_block.A_NS.cols=%ld, box->morton_index=%ld\n", neighbor_morton, block_idx, modified_block.A_NS.rows, modified_block.A_NS.cols, box->morton_index);
                    // }
                }
                
                auto& A_NS_block = scratch.eval_buffer;
                get_sliced_neighbor_block_into(
                    neighbor_morton, modified_block, level,
                    true, workspace_cols, A_NS_block
                );
                // std::vector<DataType> A_NS_block = slice_modified_block_both_directions<CoordType, DataType>(
                //     modified_block, level, neighbor_morton,
                //     box->skeleton_indices,  // Current box's skeleton (just computed in step 1)
                //     true,                   // is_A_NS
                //     n_neighbor              // Expected neighbor size
                // );

                if (DEBUG) {
                    bool block_has_nan = false;
                    for (auto val : A_NS_block) {
                        if (is_nan(val)) {
                            block_has_nan = true;
                            break;
                        }
                    }
                    std::cout << "DEBUG Box " << box->morton_index << " filling A_NS from neighbor " << neighbor_morton << std::endl;
                    std::cout << "  Block came from: near_field_modified_interactions[" << block_idx << "]" << std::endl;
                    std::cout << "  A_NS_block size: " << A_NS_block.size() << std::endl;
                    std::cout << "  A_NS_block contains NaN: " << (block_has_nan ? "YES" : "NO") << std::endl;
                    std::cout << "  current_row offset: " << current_row << std::endl;
                }
                
                int64_t actual_rows = A_NS_block.size() / workspace_cols;
                assert(actual_rows == n_neighbor);
                n_neighbor = actual_rows;
                // if(rank == 0 && VERBOSE == 1)
                // {
                //     printf("n_neighbor at 1: %ld\n", n_neighbor);
                //     printf("modified_block row: %d\n", modified_block.A_NS.rows);
                //     printf("modified_block col: %d\n", modified_block.A_NS.cols);
                // }

                for (int64_t j = 0; j < k; ++j) {
                    int64_t src_col = box->skeleton_indices[j];
                    for (int64_t i = 0; i < actual_rows; ++i) {
                        A_NS_all[(current_row + i) + j * total_neighbor_points] = 
                            A_NS_block[i + src_col * actual_rows];
                    }
                }
                
                for (int64_t j = 0; j < r; ++j) {
                    int64_t src_col = box->redundant_indices[j];
                    for (int64_t i = 0; i < actual_rows; ++i) {
                        A_NR_all[(current_row + i) + j * total_neighbor_points] = 
                            A_NS_block[i + src_col * actual_rows];
                    }
                }
                
            } else {
                BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
                PointDataRequest<CoordType>* assisting_neighbor = nullptr;
                bool is_assisting = false;

                if (neighbor_box == nullptr) {
                    auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
                    if (ghost_it != level.ghost_id_to_index.end()) {
                        neighbor_box = &level.ghost_boxes[ghost_it->second];
                    }
                }

                if (neighbor_box == nullptr) {
                    auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
                    if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
                        assisting_neighbor = &level.assisting_boxes[assist_it->second];
                        is_assisting = true;
                    }
                }

                if (neighbor_box == nullptr && assisting_neighbor == nullptr) {
                    throw std::runtime_error(
                        "compute_and_modify: Cannot find neighbor " + 
                        std::to_string(neighbor_morton));
                }

                if (is_assisting) {
                    const bool neighbor_elim =
                        (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end());

                    const int64_t n_full = static_cast<int64_t>(assisting_neighbor->coords.size() / dimension);

                    const CoordType* nb_coords = assisting_neighbor->coords.data();
                        auto& skel_coords = scratch.coord_buffer;
                        skel_coords.clear();
                        int64_t n_use = n_full;

                    if (neighbor_elim) {
                        n_use = static_cast<int64_t>(assisting_neighbor->skel_indices.size());
                        if (n_use <= 0) throw std::runtime_error("assisting neighbor eliminated but skel_indices empty");

                            skel_coords.resize(static_cast<size_t>(n_use * dimension));
                        for (int64_t i = 0; i < n_use; ++i) {
                            const int64_t src = assisting_neighbor->skel_indices[static_cast<size_t>(i)];
                            for (int d = 0; d < dimension; ++d)
                                skel_coords[static_cast<size_t>(i * dimension + d)] =
                                    assisting_neighbor->coords[static_cast<size_t>(src * dimension + d)];
                        }
                        nb_coords = skel_coords.data();
                    }

                    if (n_neighbor != n_use) {
                        printf("ERROR: n_neighbor mismatch for assisting box %ld. Expected %ld, got %ld\n",
                            neighbor_morton, n_neighbor, n_use);
                        assert(false);
                    }
                    n_neighbor = n_use;

                    auto& A_NB = scratch.eval_buffer;
                    A_NB.resize(static_cast<size_t>(n_neighbor * workspace_cols));
                    kernel->evaluate_block(
                        nb_coords, n_neighbor,
                        box->point_coords.data(), workspace_cols,
                        A_NB.data(), n_neighbor
                    );

                    for (int64_t j = 0; j < k; ++j) {
                        const int64_t src_col = box->skeleton_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i)
                            A_NS_all[(current_row + i) + j * total_neighbor_points] =
                                A_NB[i + src_col * n_neighbor];
                    }
                    for (int64_t j = 0; j < r; ++j) {
                        const int64_t src_col = box->redundant_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i)
                            A_NR_all[(current_row + i) + j * total_neighbor_points] =
                                A_NB[i + src_col * n_neighbor];
                    }
                    
                } else if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                    // Local or ghost box that has been eliminated: use skeleton
                    assert(n_neighbor == neighbor_box->skeleton_indices.size());
                    n_neighbor = neighbor_box->skeleton_indices.size();  // UPDATE n_neighbor!
                    auto& skeleton_coords = scratch.coord_buffer;
                    skeleton_coords.resize(static_cast<size_t>(n_neighbor * dimension));
                    
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        int64_t src_idx = neighbor_box->skeleton_indices[i];
                        for (int d = 0; d < dimension; ++d) {
                            skeleton_coords[i * dimension + d] = 
                                neighbor_box->point_coords[src_idx * dimension + d];
                        }
                    }
                    
                    auto& A_NB = scratch.eval_buffer;
                    A_NB.resize(static_cast<size_t>(n_neighbor * workspace_cols));
                    kernel->evaluate_block(
                        skeleton_coords.data(), n_neighbor,
                        box->point_coords.data(), workspace_cols,
                        A_NB.data(), n_neighbor
                    );
                    
                    for (int64_t j = 0; j < k; ++j) {
                        int64_t src_col = box->skeleton_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            A_NS_all[(current_row + i) + j * total_neighbor_points] = 
                                A_NB[i + src_col * n_neighbor];
                        }
                    }
                    
                    for (int64_t j = 0; j < r; ++j) {
                        int64_t src_col = box->redundant_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            A_NR_all[(current_row + i) + j * total_neighbor_points] = 
                                A_NB[i + src_col * n_neighbor];
                        }
                    }
                    
                } else {
                    // Local or ghost box that has NOT been eliminated: use full coords
                    auto& A_NB = scratch.eval_buffer;
                    A_NB.resize(static_cast<size_t>(n_neighbor * workspace_cols));
                    kernel->evaluate_block(
                        neighbor_box->point_coords.data(), n_neighbor,
                        box->point_coords.data(), workspace_cols,
                        A_NB.data(), n_neighbor
                    );
                    
                    for (int64_t j = 0; j < k; ++j) {
                        int64_t src_col = box->skeleton_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            A_NS_all[(current_row + i) + j * total_neighbor_points] = 
                                A_NB[i + src_col * n_neighbor];
                        }
                    }
                    
                    for (int64_t j = 0; j < r; ++j) {
                        int64_t src_col = box->redundant_indices[j];
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            A_NR_all[(current_row + i) + j * total_neighbor_points] = 
                                A_NB[i + src_col * n_neighbor];
                        }
                    }
                }
            }
            
            current_row += n_neighbor;
        }
        
        // Compute X_NR = A_NR - A_NS * T (batched)
        box->X_NR.set_owned(
            total_neighbor_points, r, std::move(A_NR_all), MatrixStorage<DataType>::FULL);
        
        M = total_neighbor_points; N = r; K = k;
        alpha = -1.0; beta = 1.0;
        
        gemm_("N", "N", &M, &N, &K,
               &alpha, A_NS_all.data(), &M,
               T.data.data(), &K,
               &beta, box->X_NR.data.data(), &M);

        if (DEBUG) {
            // Check A_NS_all
            bool ans_has_nan = false;
            for (auto val : A_NS_all) {
                if (is_nan(val)) {
                    ans_has_nan = true;
                    break;
                }
            }
            
            // Check A_NR_all
            bool anr_has_nan = false;
            for (auto val : A_NR_all) {
                if (is_nan(val)) {
                    anr_has_nan = true;
                    break;
                }
            }
            
            // Check X_NR
            bool xnr_has_nan = false;
            for (size_t i = 0; i < box->X_NR.data.size(); i++) {
                if (is_nan(box->X_NR.data[i])) {
                    xnr_has_nan = true;
                    break;
                }
            }
            
            std::cout << "DEBUG: After computing X_NR for Box " << box->morton_index << std::endl;
            std::cout << "  A_NS_all contains NaN: " << (ans_has_nan ? "YES" : "NO") << std::endl;
            std::cout << "  A_NR_all contains NaN: " << (anr_has_nan ? "YES" : "NO") << std::endl;
            std::cout << "  X_NR contains NaN: " << (xnr_has_nan ? "YES" : "NO") << std::endl;
            std::cout << "  X_NR dims: " << box->X_NR.rows << " × " << box->X_NR.cols << std::endl;
            std::cout << "  total_neighbor_points: " << total_neighbor_points << std::endl;
        }

        // For nonsymmetric: compute X_RN
        if (!is_symmetric && !is_hermitian) {
            auto& A_SN_all = scratch.a_sn_all;
            A_SN_all.resize(static_cast<size_t>(k * total_neighbor_points));
            std::vector<DataType> A_RN_all_nonsym(r * total_neighbor_points);
            
            current_row = 0;
            
            for (size_t idx = 0; idx < box->one_hop.size(); ++idx) {
                int64_t neighbor_morton = box->one_hop[idx];
                int64_t n_neighbor = neighbor_point_counts[idx];
                
                auto it = box->near_field_interaction_map_nonsymmetry.find(neighbor_morton);
                
                if (it != box->near_field_interaction_map_nonsymmetry.end()) {
                    int64_t block_idx = it->second;
                    const auto& modified_block = box->near_field_modified_interactions[block_idx];
                    
                    if (!modified_block.A_SN.is_allocated()) {
                        throw std::runtime_error(
                            "compute_and_modify: A_SN not allocated for neighbor " + 
                            std::to_string(neighbor_morton));
                    }
                    
                    std::vector<DataType> A_SN_block = get_sliced_neighbor_block(
                        neighbor_morton, modified_block, level,
                        false, workspace_cols
                    );
                    // std::vector<DataType> A_SN_block = get_sliced_neighbor_block(
                    //     neighbor_morton, modified_block, level,
                    //     false, workspace_cols,
                    //     // &box->skeleton_indices
                    //     nullptr
                    // );
                    
                    int64_t actual_cols = A_SN_block.size() / workspace_cols;
                    n_neighbor = actual_cols;
                    
                    for (int64_t j = 0; j < actual_cols; ++j) {
                        for (int64_t i = 0; i < k; ++i) {
                            int64_t src_row = box->skeleton_indices[i];
                            A_SN_all[i + (current_row + j) * k] = 
                                A_SN_block[src_row + j * workspace_cols];
                        }
                    }
                    
                    for (int64_t j = 0; j < actual_cols; ++j) {
                        for (int64_t i = 0; i < r; ++i) {
                            int64_t src_row = box->redundant_indices[i];
                            A_RN_all_nonsym[i + (current_row + j) * r] = 
                                A_SN_block[src_row + j * workspace_cols];
                        }
                    }
                    
                } else {
                    BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
                    PointDataRequest<CoordType>* assisting_neighbor = nullptr;
                    bool is_assisting = false;

                    if (neighbor_box == nullptr) {
                        auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
                        if (ghost_it != level.ghost_id_to_index.end()) {
                            neighbor_box = &level.ghost_boxes[ghost_it->second];
                        }
                    }

                    if (neighbor_box == nullptr) {
                        auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
                        if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
                            assisting_neighbor = &level.assisting_boxes[assist_it->second];
                            is_assisting = true;
                        }
                    }

                    if (neighbor_box == nullptr && assisting_neighbor == nullptr) {
                        throw std::runtime_error(
                            "compute_and_modify: Cannot find neighbor " + 
                            std::to_string(neighbor_morton));
                    }

                    if (is_assisting) {
                        // Assisting box: use full coords (no slicing, ignore eliminated status)
                        n_neighbor = assisting_neighbor->coords.size() / dimension;
                        
                        auto& A_BN = scratch.eval_buffer;
                        A_BN.resize(static_cast<size_t>(workspace_cols * n_neighbor));
                        kernel->evaluate_block(
                            box->point_coords.data(), workspace_cols,
                            assisting_neighbor->coords.data(), n_neighbor,
                            A_BN.data(), workspace_cols
                        );
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < k; ++i) {
                                int64_t src_row = box->skeleton_indices[i];
                                A_SN_all[i + (current_row + j) * k] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < r; ++i) {
                                int64_t src_row = box->redundant_indices[i];
                                A_RN_all_nonsym[i + (current_row + j) * r] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                        
                    } else if (level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                        // Local or ghost box that has been eliminated: use skeleton
                        n_neighbor = neighbor_box->skeleton_indices.size();  // UPDATE n_neighbor!
                        auto& skeleton_coords = scratch.coord_buffer;
                        skeleton_coords.resize(static_cast<size_t>(n_neighbor * dimension));
                        
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            int64_t src_idx = neighbor_box->skeleton_indices[i];
                            for (int d = 0; d < dimension; ++d) {
                                skeleton_coords[i * dimension + d] = 
                                    neighbor_box->point_coords[src_idx * dimension + d];
                            }
                        }
                        
                        auto& A_BN = scratch.eval_buffer;
                        A_BN.resize(static_cast<size_t>(workspace_cols * n_neighbor));
                        kernel->evaluate_block(
                            box->point_coords.data(), workspace_cols,
                            skeleton_coords.data(), n_neighbor,
                            A_BN.data(), workspace_cols
                        );
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < k; ++i) {
                                int64_t src_row = box->skeleton_indices[i];
                                A_SN_all[i + (current_row + j) * k] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < r; ++i) {
                                int64_t src_row = box->redundant_indices[i];
                                A_RN_all_nonsym[i + (current_row + j) * r] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                        
                    } else {
                        // Local or ghost box that has NOT been eliminated: use full coords
                        auto& A_BN = scratch.eval_buffer;
                        A_BN.resize(static_cast<size_t>(workspace_cols * n_neighbor));
                        kernel->evaluate_block(
                            box->point_coords.data(), workspace_cols,
                            neighbor_box->point_coords.data(), n_neighbor,
                            A_BN.data(), workspace_cols
                        );
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < k; ++i) {
                                int64_t src_row = box->skeleton_indices[i];
                                A_SN_all[i + (current_row + j) * k] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                        
                        for (int64_t j = 0; j < n_neighbor; ++j) {
                            for (int64_t i = 0; i < r; ++i) {
                                int64_t src_row = box->redundant_indices[i];
                                A_RN_all_nonsym[i + (current_row + j) * r] = 
                                    A_BN[src_row + j * workspace_cols];
                            }
                        }
                    }
                }
                
                current_row += n_neighbor;
            }
            
            // Compute X_RN = A_RN - T^T * A_SN
            box->X_RN.set_owned(
                r, total_neighbor_points, std::move(A_RN_all_nonsym), MatrixStorage<DataType>::FULL);
            
            M = r; N = total_neighbor_points; K = k;
            alpha = -1.0; beta = 1.0;
            
            gemm_("T", "N", &M, &N, &K,
                   &alpha, T.data.data(), &K,
                   A_SN_all.data(), &K,
                   &beta, box->X_RN.data.data(), &M);
        }
    }
    
    workspace.clear();
    
    // ========================================================================
    // STEP 2: Second Elimination (Section 2.4) - Call internal helper
    // ========================================================================
    
    compute_step_two_internal(
        box, level, X_BB,
        A_NS_all_unmodified, A_SN_all_unmodified,
        neighbor_point_counts,
        is_symmetric, is_hermitian,
        factorization_method,
        kernel,
        X_NN_full,
        scratch,
        pending_updates,
        store
    );
}

template<typename CoordType, typename DataType, typename KernelType>
void compute_and_modify(
    int dimension,
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    std::vector<DataType>& workspace,
    int64_t workspace_rows,
    int64_t workspace_cols,
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    std::vector<DataType>& X_NN_full,
    std::vector<DataType>& sketch_storage,
    PendingFactorUpdates<DataType>* pending_updates,
    FactorizationMethod factorization_method, bool store = false,
    int DEBUG = 0, int VERBOSE = 0) {
    FactorizationThreadScratch<CoordType, DataType> scratch;
    scratch.workspace.swap(workspace);
    scratch.x_nn_full.swap(X_NN_full);
    scratch.sketch_storage.swap(sketch_storage);
    scratch.workspace_rows = workspace_rows;
    scratch.workspace_cols = workspace_cols;

    compute_and_modify(
        dimension, box, level, kernel, scratch,
        tolerance, is_symmetric, is_hermitian,
        pending_updates, factorization_method, store, DEBUG, VERBOSE
    );

    workspace.swap(scratch.workspace);
    X_NN_full.swap(scratch.x_nn_full);
    sketch_storage.swap(scratch.sketch_storage);
}

/**
 * @brief Update neighbor slicing in modified interactions after level-wide ID
 * 
 * After all boxes at a level have been processed with ID, some neighbors may
 * have been compressed to skeleton DOFs. This function updates all modified
 * interaction blocks to reflect the correct neighbor dimensions.
 * 
 * Key insight: The S direction (current box's skeleton) should already be
 * correct from elimination. Only the N direction (neighbor's DOFs) may need
 * slicing if the neighbor was eliminated after this block was created.
 * 
 * @param level Tree level to process
 * @param is_symmetric Whether matrix is symmetric
 */
template<typename CoordType, typename DataType>
void update_neighbor_slicing_for_level(
    TreeLevel<CoordType, DataType>& level,
    bool is_symmetric) {
    
    if(!level.is_process_active) return;
    
    // Helper lambda to slice a single matrix in the neighbor dimension
    auto slice_neighbor_dimension = [&](
        MatrixStorage<DataType>& matrix,
        int64_t neighbor_morton,
        int64_t current_box_skeleton_size,
        int64_t current_box_morton,
        bool is_A_NS) {
        
        if (!matrix.is_allocated()) {
            return;  // Nothing to slice
        }
        
        int64_t stored_rows = matrix.rows;
        int64_t stored_cols = matrix.cols;
        
        
        
        // Find neighbor box - check local, ghost, then assisting
        BoxData<CoordType, DataType>* neighbor_box = nullptr;
        PointDataRequest<CoordType>* assisting_neighbor = nullptr;
        bool is_assisting = false;
        
        // Check local boxes
        neighbor_box = level.find_local_box(neighbor_morton);
        
        // Check ghost boxes
        if (neighbor_box == nullptr) {
            auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
            if (ghost_it != level.ghost_id_to_index.end()) {
                neighbor_box = &level.ghost_boxes[ghost_it->second];
            }
        }
        
        // Check assisting boxes
        if (neighbor_box == nullptr) {
            auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(neighbor_morton);
            if (assist_it != level.assisting_box_points_for_kernel_evaluation.end()) {
                assisting_neighbor = &level.assisting_boxes[assist_it->second];
                is_assisting = true;
            }
        }
        
        // Not found anywhere - throw error
        if (neighbor_box == nullptr && assisting_neighbor == nullptr) {
            throw std::runtime_error(
                "update_neighbor_slicing: Neighbor box " + std::to_string(neighbor_morton) +
                " not found in local_boxes, ghost_boxes, or assisting_boxes " +
                "(current box: " + std::to_string(current_box_morton) + ")");
        }
        
        // Get skeleton indices based on neighbor type
        std::vector<int64_t> skeleton_indices;
        bool has_skeleton = false;
        
        if (is_assisting) {
            // Assisting box - use skel_indices field
            if(level.eliminated_boxes.find(neighbor_morton) != level.eliminated_boxes.end()) {
                skeleton_indices = assisting_neighbor->skel_indices;
                has_skeleton = true;
            }else{
                has_skeleton = false;
            }
            // bool gg = !(assisting_neighbor->skel_indices.empty());
            // if(gg != has_skeleton){
            //     printf("gggggggggggggggg neighbor morton: %ld, current box morton: %ld, has skel: %d, gg: %d\n", neighbor_morton, current_box_morton, has_skeleton, gg);
            //     fflush(stdout);
            //     assert(gg == has_skeleton);
            // }
            
            
        } else {
            // Local or ghost box
            // Check elimination status and skeleton indices for consistency
            bool neighbor_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
                                       level.eliminated_boxes.end());
            has_skeleton = !neighbor_box->skeleton_indices.empty();
            
            // ERROR CHECKING: Verify consistency
            if (neighbor_eliminated && !has_skeleton) {
                throw std::runtime_error(
                    "update_neighbor_slicing: Box " + std::to_string(neighbor_morton) +
                    " is marked as eliminated but has no skeleton indices!");
            }
            
            if (!neighbor_eliminated && has_skeleton) {
                // This is actually OK - box may have been ID'd but not yet eliminated
                // Just proceed with slicing if needed
            }
            
            if (has_skeleton) {
                skeleton_indices = neighbor_box->skeleton_indices;
            }
        }

        
        if (!has_skeleton) {
            // Neighbor not yet processed - nothing to slice
            return;
        }
        
        int64_t skeleton_size = skeleton_indices.size();
        
        if (is_A_NS) {
            // A_NS: (n_neighbor × k_box)
            // Check rows (neighbor dimension)
            
            if (stored_rows == skeleton_size) {
                // Already correct size
                return;
            } else if (stored_rows > skeleton_size) {
                // Need to slice rows to skeleton
                std::vector<DataType> sliced(skeleton_size * stored_cols);
                
                // Extract skeleton rows, all columns
                std::vector<int64_t> all_cols(stored_cols);
                std::iota(all_cols.begin(), all_cols.end(), 0);
                
                extract_submatrix(
                    matrix.data.data(), stored_rows,
                    skeleton_indices,
                    all_cols,
                    sliced.data(), skeleton_size
                );
                
                // Update matrix in-place
                matrix.data = std::move(sliced);
                matrix.rows = skeleton_size;
                matrix.lda = skeleton_size;
                
            } else {
                // stored_rows < skeleton_size
                throw std::runtime_error(
                    "update_neighbor_slicing: Box " + std::to_string(current_box_morton) +
                    " A_NS has fewer rows (" + std::to_string(stored_rows) + 
                    ") than neighbor's skeleton size (" + std::to_string(skeleton_size) + 
                    ") for neighbor " + std::to_string(neighbor_morton));
            }
            
        } else {
            // A_SN: (k_box × n_neighbor)
            // Check cols (neighbor dimension)
            
            if (stored_cols == skeleton_size) {
                // Already correct size
                return;
            } else if (stored_cols > skeleton_size) {
                // Need to slice columns to skeleton
                std::vector<DataType> sliced(stored_rows * skeleton_size);
                
                // Extract all rows, skeleton columns
                std::vector<int64_t> all_rows(stored_rows);
                std::iota(all_rows.begin(), all_rows.end(), 0);
                
                extract_submatrix(
                    matrix.data.data(), stored_rows,
                    all_rows,
                    skeleton_indices,
                    sliced.data(), stored_rows
                );
                
                // Update matrix in-place
                matrix.data = std::move(sliced);
                matrix.cols = skeleton_size;
                
            } else {
                // stored_cols < skeleton_size
                throw std::runtime_error(
                    "update_neighbor_slicing: Box " + std::to_string(current_box_morton) +
                    " A_SN has fewer cols (" + std::to_string(stored_cols) + 
                    ") than neighbor's skeleton size (" + std::to_string(skeleton_size) + 
                    ") for neighbor " + std::to_string(neighbor_morton));
            }
        }
    };
    
    // Process all local boxes
    for (auto& box : level.local_boxes) {

        int64_t k_box = box.skeleton_indices.size();
        

        // Update near-field blocks
        for (auto& modified_block : box.near_field_modified_interactions) {
            int64_t neighbor_morton = modified_block.neighbor_morton;
            
            // Symmetric: only A_NS needs checking
            if (modified_block.A_NS.is_allocated()) {
                slice_neighbor_dimension(
                    modified_block.A_NS, neighbor_morton, 
                    k_box, box.morton_index, true);
            }
            
            // Nonsymmetric: also check A_SN
            if (!is_symmetric && modified_block.A_SN.is_allocated()) {
                slice_neighbor_dimension(
                    modified_block.A_SN, neighbor_morton, 
                    k_box, box.morton_index, false);
            }
        }
        
        // Update far-field blocks
        for (auto& modified_block : box.far_field_modified_interactions) {
            int64_t neighbor_morton = modified_block.neighbor_morton;
            
            // Symmetric: only A_NS needs checking
            if (modified_block.A_NS.is_allocated()) {
                slice_neighbor_dimension(
                    modified_block.A_NS, neighbor_morton, 
                    k_box, box.morton_index, true);
            }
            
            // Nonsymmetric: also check A_SN
            if (!is_symmetric && modified_block.A_SN.is_allocated()) {
                slice_neighbor_dimension(
                    modified_block.A_SN, neighbor_morton, 
                    k_box, box.morton_index, false);
            }
        }
    }
}


/**
 * @brief Verify symmetry of modified interaction blocks
 * 
 * For symmetric kernels, verifies that for every box B and its neighbor G:
 *   A_NS[B,G]^T = A_NS[G,B]
 * 
 * Where:
 *   A_NS[B,G] is stored in B's modified_interactions (k_G × k_B)
 *   A_NS[G,B] is stored in G's modified_interactions (k_B × k_G)
 * 
 * This checks that the distributed storage maintains exact transpose symmetry.
 * 
 * @param level Tree level to verify
 * @param tolerance Absolute tolerance for floating point comparison
 * @param verbose Print detailed mismatch information
 * @return true if all blocks are symmetric, false otherwise
 */
template<typename CoordType, typename DataType>
bool verify_modified_interaction_symmetry(
    TreeLevel<CoordType, DataType>& level,
    DataType tolerance = 1e-12,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    bool all_symmetric = true;
    int64_t blocks_checked = 0;
    int64_t blocks_failed = 0;
    DataType max_error = 0.0;
    
    if (verbose && rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Verifying Modified Interaction Symmetry" << std::endl;
        std::cout << "Level: " << level.level << std::endl;
        std::cout << "Tolerance: " << std::scientific << tolerance << std::endl;
        std::cout << "========================================\n" << std::endl;
    }
    
    // Process all local boxes
    for (auto& box_B : level.local_boxes) {
        int64_t morton_B = box_B.morton_index;
        int64_t k_B = box_B.skeleton_indices.size();
        
        if (k_B == 0) {
            continue;  // Box not yet processed
        }
        
        // Lambda to verify one interaction list
        auto verify_interaction_list = [&](
            std::vector<ModifiedBlock<DataType>>& interactions,
            const char* list_name) {
            
            for (auto& block_BG : interactions) {
                int64_t morton_G = block_BG.neighbor_morton;
                
                if (!block_BG.A_NS.is_allocated()) {
                    continue;  // Nothing to verify
                }
                
                // A_NS[B,G] stored in box B's modified_interactions
                auto& A_NS_BG = block_BG.A_NS;
                int64_t rows_BG = A_NS_BG.rows;  // k_G (neighbor's skeleton size)
                int64_t cols_BG = A_NS_BG.cols;  // k_B (box's skeleton size)
                
                if (cols_BG != k_B) {
                    if (verbose) {
                        std::cout << "Rank " << rank << " ERROR: Box " << morton_B 
                                  << " " << list_name << " A_NS[" << morton_B << "," << morton_G 
                                  << "] has cols=" << cols_BG << " but k_B=" << k_B << std::endl;
                    }
                    all_symmetric = false;
                    continue;
                }
                
                // Find neighbor box G (could be local, ghost, or assisting)
                BoxData<CoordType, DataType>* box_G = nullptr;
                
                // Check local boxes
                box_G = level.find_local_box(morton_G);
                
                // Check ghost boxes
                if (box_G == nullptr) {
                    auto ghost_it = level.ghost_id_to_index.find(morton_G);
                    if (ghost_it != level.ghost_id_to_index.end()) {
                        box_G = &level.ghost_boxes[ghost_it->second];
                    }
                }
                
                if (box_G == nullptr) {
                    if (verbose) {
                        std::cout << "Rank " << rank << " WARNING: Box " << morton_B 
                                  << " neighbor " << morton_G 
                                  << " not found in local/ghost boxes (skipping verification)" 
                                  << std::endl;
                    }
                    continue;
                }
                
                int64_t k_G = box_G->skeleton_indices.size();
                
                if (k_G == 0) {
                    continue;  // Neighbor not yet processed
                }
                
                if (rows_BG != k_G) {
                    if (verbose) {
                        std::cout << "Rank " << rank << " ERROR: Box " << morton_B 
                                  << " " << list_name << " A_NS[" << morton_B << "," << morton_G 
                                  << "] has rows=" << rows_BG << " but k_G=" << k_G << std::endl;
                    }
                    all_symmetric = false;
                    continue;
                }
                
                // Find A_NS[G,B] in neighbor G's modified_interactions
                ModifiedBlock<DataType>* block_GB = nullptr;
                
                // Check both near and far field lists
                for (auto& candidate : box_G->near_field_modified_interactions) {
                    if (candidate.neighbor_morton == morton_B && candidate.A_NS.is_allocated()) {
                        block_GB = &candidate;
                        break;
                    }
                }
                
                if (block_GB == nullptr) {
                    for (auto& candidate : box_G->far_field_modified_interactions) {
                        if (candidate.neighbor_morton == morton_B && candidate.A_NS.is_allocated()) {
                            block_GB = &candidate;
                            break;
                        }
                    }
                }
                
                if (block_GB == nullptr) {
                    if (verbose) {
                        std::cout << "Rank " << rank << " WARNING: Box " << morton_G 
                                  << " does not have reciprocal interaction with box " 
                                  << morton_B << " (asymmetric graph)" << std::endl;
                    }
                    continue;
                }
                
                auto& A_NS_GB = block_GB->A_NS;
                int64_t rows_GB = A_NS_GB.rows;  // Should be k_B
                int64_t cols_GB = A_NS_GB.cols;  // Should be k_G
                if(block_GB->neighbor_morton == 67)
                {
                    int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    printf("rank: %d, my id: %d, 67 index row: %d, col: %d\n", rank, box_G->morton_index, A_NS_GB.rows, A_NS_GB.cols);
                }
                
                // Dimension check: A_NS[B,G]^T should have same dimensions as A_NS[G,B]
                if (rows_GB != k_B || cols_GB != k_G) {
                    
                    if (verbose) {
                        std::cout << "Rank " << rank << " ERROR: Box " << morton_G 
                                  << " A_NS[" << morton_G << "," << morton_B 
                                  << "] has wrong dimensions: (" << rows_GB << "×" << cols_GB 
                                  << ") expected (" << k_B << "×" << k_G << ")" << std::endl;
                    }
                    all_symmetric = false;
                    continue;
                }
                
                // Entry-wise comparison: A_NS[B,G]^T[i,j] == A_NS[G,B][i,j]
                // A_NS[B,G]^T[i,j] = A_NS[B,G][j,i]
                
                bool block_symmetric = true;
                DataType block_max_error = 0.0;
                int64_t first_error_i = -1, first_error_j = -1;
                DataType first_error_val_BG = 0.0, first_error_val_GB = 0.0;
                
                for (int64_t i = 0; i < k_B; ++i) {
                    for (int64_t j = 0; j < k_G; ++j) {
                        // A_NS[B,G][j,i] (j-th row, i-th col) - column-major indexing
                        DataType val_BG_transpose = A_NS_BG.data[i * rows_BG + j];
                        
                        // A_NS[G,B][i,j] (i-th row, j-th col) - column-major indexing
                        DataType val_GB = A_NS_GB.data[j * rows_GB + i];
                        
                        DataType error = std::abs(val_BG_transpose - val_GB);
                        
                        
                        if (error > tolerance) {
                            block_symmetric = false;
                            if (error > block_max_error) {
                                block_max_error = error;
                                first_error_i = i;
                                first_error_j = j;
                                first_error_val_BG = val_BG_transpose;
                                first_error_val_GB = val_GB;
                            }
                        }
                        
                        if (error > max_error) {
                            max_error = error;
                        }
                    }
                }
                
                blocks_checked++;
                
                if (!block_symmetric) {
                    blocks_failed++;
                    all_symmetric = false;
                    
                    if (verbose) {
                        std::cout << "Rank " << rank << " FAILED: " << list_name 
                                  << " A_NS[" << morton_B << "," << morton_G << "]^T"
                                  << " != A_NS[" << morton_G << "," << morton_B << "]" << std::endl;
                        std::cout << "  Max error: " << std::scientific << block_max_error << std::endl;
                        std::cout << "  First error at (" << first_error_i << "," << first_error_j << "): "
                                  << first_error_val_BG << " vs " << first_error_val_GB << std::endl;
                    }
                }
            }
        };
        
        // Verify near-field interactions
        verify_interaction_list(box_B.near_field_modified_interactions, "near-field");
        
        // Verify far-field interactions
        verify_interaction_list(box_B.far_field_modified_interactions, "far-field");
    }
    
    if (verbose) {
        std::cout << "\nRank " << rank << " Summary:" << std::endl;
        std::cout << "  Blocks checked: " << blocks_checked << std::endl;
        std::cout << "  Blocks failed: " << blocks_failed << std::endl;
        std::cout << "  Max error: " << std::scientific << max_error << std::endl;
        std::cout << "  Result: " << (all_symmetric ? "✓ PASS" : "✗ FAIL") << std::endl;
        std::cout << "========================================\n" << std::endl;
    }
    
    return all_symmetric;
}


/**
 * @brief Extract child-level interaction C(child_i, child_j)
 * 
 * Returns modified interaction from child level, or evaluates kernel for >= 3-hop
 * Returns matrix in (child_i × child_j) orientation for assembly
 */
template<typename CoordType, typename DataType, typename KernelType>
std::vector<DataType> extract_child_interaction(
    BoxData<CoordType, DataType>* child_i,
    BoxData<CoordType, DataType>* child_j,
    TreeLevel<CoordType, DataType>& child_level,
    int dimension,
    KernelType* kernel,
    bool DEBUG = false) {
    
    int64_t n_i = child_i->skeleton_indices.size();
    int64_t n_j = child_j->skeleton_indices.size();
    
    // Diagonal case: use schur complement (if available)
    if (child_i->morton_index == child_j->morton_index) {
        if (!child_i->schur_complement.is_allocated()) {
            if (DEBUG){
                // WARNING: Schur complement not allocated - using direct kernel evaluation
                std::cerr << "WARNING: extract_child_interaction: Child box " 
                        << child_i->morton_index 
                        << " schur complement not allocated (skeleton=" << n_i 
                        << ", redundant=" << child_i->redundant_indices.size()
                        << "). Falling back to kernel evaluation." << std::endl;
            }
           
            
            // Extract skeleton coordinates and evaluate kernel directly
            std::vector<DataType> C_block(n_i * n_i);
            std::vector<CoordType> coords_i(n_i * dimension);
            
            for (int64_t idx = 0; idx < n_i; ++idx) {
                int64_t skel_idx = child_i->skeleton_indices[idx];
                for (int d = 0; d < dimension; ++d) {
                    coords_i[idx * dimension + d] = child_i->point_coords[skel_idx * dimension + d];
                }
            }
            
            // Evaluate kernel: self-interaction
            kernel->evaluate_block(
                coords_i.data(), n_i,
                coords_i.data(), n_i,
                C_block.data(), n_i);
            
            return C_block;
        }
        
        return child_i->schur_complement.data;
    }
    
    // Check hop distance between children
    uint32_t i_x, i_y, i_z = 0, j_x, j_y, j_z = 0;
    if (dimension == 2) {
        morton::decode_2d(child_i->morton_index, i_x, i_y);
        morton::decode_2d(child_j->morton_index, j_x, j_y);
    } else {
        morton::decode_3d(child_i->morton_index, i_x, i_y, i_z);
        morton::decode_3d(child_j->morton_index, j_x, j_y, j_z);
    }
    
    int64_t dx = std::abs(static_cast<int64_t>(i_x) - static_cast<int64_t>(j_x));
    int64_t dy = std::abs(static_cast<int64_t>(i_y) - static_cast<int64_t>(j_y));
    int64_t dz = (dimension == 3) ? std::abs(static_cast<int64_t>(i_z) - static_cast<int64_t>(j_z)) : 0;
    
    bool is_one_hop = (dx <= 1 && dy <= 1 && (dimension == 2 || dz <= 1));
    bool is_two_hop = (dx <= 2 && dy <= 2 && (dimension == 2 || dz <= 2)) &&
                     (dx == 2 || dy == 2 || (dimension == 3 && dz == 2));
    
    // 1-hop or 2-hop: look for modified interaction
    if (is_one_hop || is_two_hop) {
        auto& interaction_map = is_one_hop ? 
            child_i->near_field_interaction_map :
            child_i->far_field_interaction_map;
        
        auto& modified_interactions = is_one_hop ?
            child_i->near_field_modified_interactions :
            child_i->far_field_modified_interactions;
        
        auto it = interaction_map.find(child_j->morton_index);
        if (it == interaction_map.end()) {
            if (DEBUG){
                printf("near field size: %zu, far field size: %zu\n", 
                    child_i->near_field_modified_interactions.size(),
                    child_i->far_field_modified_interactions.size());
                // WARNING: Modified interaction not found - fall back to kernel evaluation
                std::cerr << "WARNING: extract_child_interaction: " 
                        << std::string(is_one_hop ? "1-hop" : "2-hop") 
                        << " interaction not found: child " 
                        << child_i->morton_index << " -> " << child_j->morton_index
                        << ". Falling back to kernel evaluation." << std::endl;
            }
            // Fall through to kernel evaluation below
        } else {
            // *** KEY FIX: A_NS is stored as (neighbor × self) = (child_j × child_i)
            // *** We need to transpose to get C(child_i, child_j) = (child_i × child_j)
            const auto& A_NS = modified_interactions[it->second].A_NS;

            if (A_NS.rows != n_j || A_NS.cols != n_i) {
                if (DEBUG){
                    std::cerr << "WARNING: extract_child_interaction: A_NS dimension mismatch. "
                          << "Expected " << n_j << "×" << n_i
                          << " got " << A_NS.rows << "×" << A_NS.cols
                          << ". Falling back to kernel evaluation." << std::endl;
                }
            } else {
                // Transpose: A_NS is (n_j × n_i), return (n_i × n_j)
                std::vector<DataType> C_block(n_i * n_j);
                for (int64_t i = 0; i < n_i; ++i) {
                    for (int64_t j = 0; j < n_j; ++j) {
                        C_block[i + j * n_i] = A_NS.data[j + i * n_j];
                    }
                }
                
                return C_block;
            }
        }
    }
    
    // >= 3 hops OR fallback case: evaluate kernel directly
    std::vector<DataType> C_block(n_i * n_j);
    
    // Extract skeleton coordinates from children
    std::vector<CoordType> coords_i(n_i * dimension);
    std::vector<CoordType> coords_j(n_j * dimension);
    
    for (int64_t idx = 0; idx < n_i; ++idx) {
        int64_t skel_idx = child_i->skeleton_indices[idx];
        for (int d = 0; d < dimension; ++d) {
            coords_i[idx * dimension + d] = child_i->point_coords[skel_idx * dimension + d];
        }
    }
    
    for (int64_t idx = 0; idx < n_j; ++idx) {
        int64_t skel_idx = child_j->skeleton_indices[idx];
        for (int d = 0; d < dimension; ++d) {
            coords_j[idx * dimension + d] = child_j->point_coords[skel_idx * dimension + d];
        }
    }
    
    // Evaluate kernel: returns (n_i × n_j)
    kernel->evaluate_block(
        coords_i.data(), n_i,
        coords_j.data(), n_j,
        C_block.data(), n_i);
    
    return C_block;
}


// /**
//  * @brief Build parent level boxes and interactions from child level
//  * 
//  * Transitions from child level to parent level by:
//  * 1. Accumulating skeleton points from children to parents
//  * 2. Assembling parent-level interactions from child-level modified interactions
//  * 
//  * @param child_level Current level (being compressed)
//  * @param parent_level Next coarser level (to be built)
//  * @param dimension 2 or 3
//  * @param is_symmetric Whether matrix is symmetric
//  * @param is_hermitian Whether matrix is Hermitian
//  * @param kernel Kernel evaluator for >= 3-hop interactions
//  */
// template<typename CoordType, typename DataType, typename KernelType>
// void build_parent_level_interactions(
//     TreeLevel<CoordType, DataType>& child_level,
//     TreeLevel<CoordType, DataType>& parent_level,
//     int dimension,
//     bool is_symmetric,
//     bool is_hermitian,
//     KernelType* kernel) {
    
//     int num_children = (dimension == 2) ? 4 : 8;
    
//     // Calculate number of parent boxes this process owns
//     int64_t num_parent_boxes = child_level.local_boxes.size() / num_children;
    
//     if (child_level.local_boxes.size() % num_children != 0) {
//         throw std::runtime_error(
//             "build_parent_level: Child boxes not evenly divisible by " + 
//             std::to_string(num_children));
//     }
    
//     parent_level.local_boxes.resize(num_parent_boxes);
//     parent_level.num_boxes_local = num_parent_boxes;
//     parent_level.level = child_level.level - 1;
    
//     // ===== Step 1: Initialize parent boxes and accumulate skeleton points =====
    
//     for (int64_t parent_idx = 0; parent_idx < num_parent_boxes; ++parent_idx) {
//         auto& parent_box = parent_level.local_boxes[parent_idx];
        
//         // Get first child index
//         int64_t first_child_idx = parent_idx * num_children;
//         auto& first_child = child_level.local_boxes[first_child_idx];
        
//         // Calculate parent Morton index
//         parent_box.morton_index = first_child.morton_index / num_children;
//         parent_box.level = parent_level.level;
        
//         // Decode grid coordinates from parent Morton
//         if (dimension == 2) {
//             uint32_t x, y;
//             morton::decode_2d(parent_box.morton_index, x, y);
//             parent_box.grid_coords[0] = static_cast<int32_t>(x);
//             parent_box.grid_coords[1] = static_cast<int32_t>(y);
//             parent_box.grid_coords[2] = 0;
//         } else {
//             uint32_t x, y, z;
//             morton::decode_3d(parent_box.morton_index, x, y, z);
//             parent_box.grid_coords[0] = static_cast<int32_t>(x);
//             parent_box.grid_coords[1] = static_cast<int32_t>(y);
//             parent_box.grid_coords[2] = static_cast<int32_t>(z);
//         }
        
//         // Calculate parent box geometry (2x size of child)
//         parent_box.size = first_child.size * 2.0;
        
//         // Parent bounds encompass all children
//         parent_box.bounds[0] = first_child.bounds[0]; // xmin
//         parent_box.bounds[1] = first_child.bounds[1] + first_child.size; // xmax
//         parent_box.bounds[2] = first_child.bounds[2]; // ymin
//         parent_box.bounds[3] = first_child.bounds[3] + first_child.size; // ymax
//         if (dimension == 3) {
//             parent_box.bounds[4] = first_child.bounds[4]; // zmin
//             parent_box.bounds[5] = first_child.bounds[5] + first_child.size; // zmax
//         } else {
//             parent_box.bounds[4] = 0;
//             parent_box.bounds[5] = 0;
//         }
        
//         // Calculate center
//         for (int d = 0; d < dimension; ++d) {
//             parent_box.center[d] = (parent_box.bounds[2*d] + parent_box.bounds[2*d+1]) / 2.0;
//         }
//         if (dimension == 2) {
//             parent_box.center[2] = 0;
//         }
        
//         // Count total skeleton points from all children
//         int64_t total_skeleton_points = 0;
//         for (int c = 0; c < num_children; ++c) {
//             auto& child = child_level.local_boxes[first_child_idx + c];
//             total_skeleton_points += child.skeleton_indices.size();
//         }
        
//         // Reserve space
//         parent_box.point_indices.reserve(total_skeleton_points);
//         parent_box.point_coords.reserve(total_skeleton_points * dimension);
        
//         // Accumulate skeleton points from all children (in Morton order)
//         for (int c = 0; c < num_children; ++c) {
//             auto& child = child_level.local_boxes[first_child_idx + c];
            
//             // Append skeleton indices (these become parent's point indices)
//             parent_box.point_indices.insert(
//                 parent_box.point_indices.end(),
//                 child.skeleton_indices.begin(),
//                 child.skeleton_indices.end());
            
//             // Append skeleton point coordinates
//             for (int64_t skel_idx : child.skeleton_indices) {
//                 for (int d = 0; d < dimension; ++d) {
//                     parent_box.point_coords.push_back(
//                         child.point_coords[skel_idx * dimension + d]);
//                 }
//             }
//         }
        
//         parent_box.num_points = total_skeleton_points;
        
//         // Initialize parent/children relationships
//         parent_box.num_children = num_children;
//         for (int c = 0; c < num_children; ++c) {
//             parent_box.children_morton[c] = child_level.local_boxes[first_child_idx + c].morton_index;
//         }
//         for (int c = num_children; c < 8; ++c) {
//             parent_box.children_morton[c] = -1;
//         }
//     }
    
//     // ===== Step 2: Build parent-level modified interactions =====
    
//     for (size_t b1_idx = 0; b1_idx < parent_level.local_boxes.size(); ++b1_idx) {
//         auto& B1 = parent_level.local_boxes[b1_idx];
        
//         // Lower triangular for symmetric, full for nonsymmetric
//         size_t b2_start = (is_symmetric || is_hermitian) ? b1_idx : 0;
        
//         for (size_t b2_idx = b2_start; b2_idx < parent_level.local_boxes.size(); ++b2_idx) {
//             auto& B2 = parent_level.local_boxes[b2_idx];
            
//             // Classify relationship between parent boxes
//             bool is_diagonal = (b1_idx == b2_idx);
            
//             int64_t dx = std::abs(B1.grid_coords[0] - B2.grid_coords[0]);
//             int64_t dy = std::abs(B1.grid_coords[1] - B2.grid_coords[1]);
//             int64_t dz = (dimension == 3) ? std::abs(B1.grid_coords[2] - B2.grid_coords[2]) : 0;
            
//             bool is_one_hop = (dx <= 1 && dy <= 1 && (dimension == 2 || dz <= 1));

            
//             if (!is_diagonal && !is_one_hop) {
//                 continue; // >= 3 hops, skip
//             }
            
//             // Calculate dimensions of assembled interaction
//             int64_t total_rows = B1.num_points;
//             int64_t total_cols = B2.num_points;
            
//             // Preallocate I(B1, B2)
//             std::vector<DataType> I_B1_B2(total_rows * total_cols);
            
//             // Assemble I(B1, B2) from child interactions
//             int64_t row_offset = 0;
//             int64_t first_child_b1 = b1_idx * num_children;
            
//             for (int ci = 0; ci < num_children; ++ci) {
//                 auto& child_i = child_level.local_boxes[first_child_b1 + ci];
//                 int64_t n_i = child_i.skeleton_indices.size();
                
//                 int64_t col_offset = 0;
//                 int64_t first_child_b2 = b2_idx * num_children;
                
//                 for (int cj = 0; cj < num_children; ++cj) {
//                     auto& child_j = child_level.local_boxes[first_child_b2 + cj];
//                     int64_t n_j = child_j.skeleton_indices.size();
                    
//                     // Extract C(child_i, child_j)
//                     std::vector<DataType> C_block = extract_child_interaction(
//                         &child_i, &child_j, child_level, dimension, kernel);
                    
//                     // Verify size
//                     if (C_block.size() != n_i * n_j) {
//                         throw std::runtime_error(
//                             "extract_child_interaction size mismatch: expected " +
//                             std::to_string(n_i * n_j) + " got " + std::to_string(C_block.size()));
//                     }
                    
//                     // Copy C_block into I_B1_B2 (column-major)
//                     for (int64_t col = 0; col < n_j; ++col) {
//                         for (int64_t row = 0; row < n_i; ++row) {
//                             I_B1_B2[(row_offset + row) + (col_offset + col) * total_rows] = 
//                                 C_block[row + col * n_i];
//                         }
//                     }
                    
//                     col_offset += n_j;
//                 }
                
//                 row_offset += n_i;
//             }
            
//             // Store I(B1, B2) in appropriate location
//             if (is_diagonal) {
//                 // Goes to B1's schur complement
//                 B1.schur_complement.allocate(total_rows, total_cols, MatrixStorage<DataType>::FULL);
//                 B1.schur_complement.data = std::move(I_B1_B2);
                
//             } else if (is_one_hop) {
//                 // B1's view of B2: transpose to (B2 × B1)
//                 std::vector<DataType> I_transposed(total_rows * total_cols);
//                 for (int64_t i = 0; i < total_rows; ++i) {
//                     for (int64_t j = 0; j < total_cols; ++j) {
//                         I_transposed[j + i * total_cols] = I_B1_B2[i + j * total_rows];
//                     }
//                 }
                
//                 ModifiedBlock<DataType> block_b1;
//                 block_b1.neighbor_morton = B2.morton_index;
//                 block_b1.A_NS.allocate(total_cols, total_rows, MatrixStorage<DataType>::FULL);  // (B2 × B1)
//                 block_b1.A_NS.data = I_transposed;  // Copy (we need it for B2 in nonsymmetric case)
                
//                 int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
//                 B1.near_field_modified_interactions.push_back(std::move(block_b1));
                
//                 if (is_symmetric || is_hermitian) {
//                     B1.near_field_interaction_map[B2.morton_index] = block_idx_b1;
//                 } else {
//                     B1.near_field_interaction_map_nonsymmetry[B2.morton_index] = block_idx_b1;
//                 }
                
//                 // B2's view of B1
//                 if (is_symmetric || is_hermitian) {
//                     // Symmetric: A_NS should be (B1 × B2)
//                     ModifiedBlock<DataType> block_b2;
//                     block_b2.neighbor_morton = B1.morton_index;
//                     block_b2.A_NS.allocate(total_rows, total_cols, MatrixStorage<DataType>::FULL);  // (B1 × B2)
//                     block_b2.A_NS.data = std::move(I_B1_B2);
                    
//                     int64_t block_idx_b2 = B2.near_field_modified_interactions.size();
//                     B2.near_field_modified_interactions.push_back(std::move(block_b2));
//                     B2.near_field_interaction_map[B1.morton_index] = block_idx_b2;
                    
//                 } else {
//                     // ✓ FIX: Nonsymmetric A_SN should be (B2 × B1), same as A_NS!
//                     ModifiedBlock<DataType> block_b2;
//                     block_b2.neighbor_morton = B1.morton_index;
//                     block_b2.A_SN.allocate(total_cols, total_rows, MatrixStorage<DataType>::FULL);  // ✓ (B2 × B1)
//                     block_b2.A_SN.data = std::move(I_transposed);  // ✓ Use transposed version
                    
//                     int64_t block_idx_b2 = B2.near_field_modified_interactions.size();
//                     B2.near_field_modified_interactions.push_back(std::move(block_b2));
//                     B2.near_field_interaction_map_nonsymmetry[B1.morton_index] = block_idx_b2;
//                 }
                
//             } 
//         }
//     }
// }



/**
 * @brief Extract skeleton coordinates from PointDataRequest
 * @param request The PointDataRequest containing coords and skel_indices
 * @param dimension Spatial dimension (2 or 3)
 * @return Sliced coordinates containing only skeleton points (column-major)
 */
template<typename CoordType>
std::vector<CoordType> extract_skeleton_coords(
    const PointDataRequest<CoordType>& request,
    int dimension
) {
    if (request.skel_indices.empty()) {
        // No skeleton indices - return full coords
        return request.coords;
    }
    
    int64_t n_skel = request.skel_indices.size();
    std::vector<CoordType> skeleton_coords(n_skel * dimension);
    
    for (int64_t i = 0; i < n_skel; ++i) {
        int64_t src_idx = request.skel_indices[i];
        for (int d = 0; d < dimension; ++d) {
            skeleton_coords[i * dimension + d] = 
                request.coords[src_idx * dimension + d];
        }
    }
    
    return skeleton_coords;
}




/**
 * @brief Try to extract child interaction from modified blocks, fallback to kernel evaluation
 * 
 * Checks near_field and far_field modified interactions for pre-computed block.
 * If not found, evaluates kernel directly.
 * 
 * Storage convention: A_NS is always (n_neighbor × n_source)
 * - Looking in source's A_NS[target]: need transpose to get (n_source × n_target)
 * - Looking in target's A_NS[source]: already (n_source × n_target), no transpose
 * 
 * @param source_child Child box that owns the interaction (where to look for modified blocks)
 * @param target_morton Morton index of the target child box
 * @param target_coords Skeleton coordinates of target child
 * @param n_source Number of skeleton points in source child
 * @param n_target Number of skeleton points in target child
 * @param child_level Tree level containing children
 * @param dimension 2 or 3
 * @param kernel Kernel evaluator
 * @param transpose_if_found If true, transpose A_NS when found (looking in source's view)
 * @return Interaction block (n_source × n_target, column-major)
 */
template<typename CoordType, typename DataType, typename KernelType>
std::vector<DataType> extract_or_evaluate_child_interaction_for_assisting(
    BoxData<CoordType, DataType>* source_child,
    int64_t target_morton,
    const std::vector<CoordType>& target_coords,
    int64_t n_source,
    int64_t n_target,
    TreeLevel<CoordType, DataType>& child_level,
    int dimension,
    KernelType* kernel,
    bool transpose_if_found = true) {
    

    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(n_target > 50)
    // {
    //     printf("aaaaaaaaaaaaaaaa rank: %d, source: %d, target: %d, n_source: %d, n_target: %d\n", rank, source_child->morton_index, target_morton, n_source, n_target);
    //     fflush(stdout);
    //     assert(n_target < 50);
    // }
    // if(source_child->morton_index == 63 && target_morton == 150){
    //     printf("aaaaaaaaaaaaaaaa rank: %d, source: %d, target: %d, n_source: %d, n_target: %d, one_hop_map_size: %d, two_hop_map_size: %d\n", rank, source_child->morton_index, target_morton, n_source, n_target, source_child->near_field_interaction_map.size(), source_child->far_field_interaction_map.size());
    //     for (auto it = source_child->far_field_interaction_map.begin(); it != source_child->far_field_interaction_map.end(); ++it) {
    //         std::cout << "Key: " << it->first << std::endl;
    //         std::cout << "found: " << (source_child->far_field_interaction_map.find(target_morton) != source_child->far_field_interaction_map.end()) << std::endl;
    //     }
    // }
    
    
    // Try to find in near_field_modified_interactions
    auto near_it = source_child->near_field_interaction_map.find(target_morton);
    if (near_it != source_child->near_field_interaction_map.end()) {
        auto& block = source_child->near_field_modified_interactions[near_it->second];
        
        if (block.A_NS.is_allocated()) {
            // Found A_NS
            if (transpose_if_found) {
                // A_NS is (n_target × n_source), need (n_source × n_target)
                std::vector<DataType> transposed(n_source * n_target);
                for (int64_t i = 0; i < n_target; ++i) {
                    for (int64_t j = 0; j < n_source; ++j) {
                        transposed[j + i * n_source] = block.A_NS.data[i + j * n_target];
                    }
                }
                return transposed;
            } else {
                // Already correct orientation
                return block.A_NS.data;
            }
        }
    }
    
    // Try to find in far_field_modified_interactions
    auto far_it = source_child->far_field_interaction_map.find(target_morton);
    if (far_it != source_child->far_field_interaction_map.end()) {
        auto& block = source_child->far_field_modified_interactions[far_it->second];
        
        if (block.A_NS.is_allocated()) {
            // Found A_NS
            if (transpose_if_found) {
                // A_NS is (n_target × n_source), need (n_source × n_target)
                std::vector<DataType> transposed(n_source * n_target);
                for (int64_t i = 0; i < n_target; ++i) {
                    for (int64_t j = 0; j < n_source; ++j) {
                        transposed[j + i * n_source] = block.A_NS.data[i + j * n_target];
                    }
                }
                return transposed;
            } else {
                // Already correct orientation
                return block.A_NS.data;
            }
        }
    }
    
    // Not found - fallback to kernel evaluation
    std::vector<DataType> C_block(n_source * n_target);
    
    // Extract skeleton coords from source_child
    std::vector<CoordType> source_coords(n_source * dimension);
    for (int64_t idx = 0; idx < n_source; ++idx) {
        int64_t skel_idx = source_child->skeleton_indices[idx];
        for (int d = 0; d < dimension; ++d) {
            source_coords[idx * dimension + d] = 
                source_child->point_coords[skel_idx * dimension + d];
        }
    }
    
    kernel->evaluate_block(
        source_coords.data(), n_source,
        target_coords.data(), n_target,
        C_block.data(), n_source);
    
    return C_block;
}

/**
 * @brief Build parent level boxes and interactions from child level
 * 
 * Transitions from child level to parent level by:
 * 1. Accumulating skeleton points from children to parents
 * 2. Assembling parent-level interactions from child-level modified interactions
 * 
 * Handles both on-process and off-process (ghost/assisting) neighbors:
 * - Symmetric: Store only B1's view for ghosts, mutual views for on-process (lower triangular)
 * - Nonsymmetric: Store both A_NS and A_SN for ghosts, A_NS(B1)+A_SN(B2) for on-process
 * 
 * @param child_level Current level (being compressed)
 * @param dimension 2 or 3
 * @param is_symmetric Whether matrix is symmetric
 * @param is_hermitian Whether matrix is Hermitian
 * @param kernel Kernel evaluator for >= 3-hop interactions
 * @param global_bounds Global domain bounds [xmin, xmax, ymin, ymax, zmin, zmax]
 * @return Vector of parent level local boxes
 */
template<typename CoordType, typename DataType, typename KernelType>
std::vector<BoxData<CoordType, DataType>> build_parent_level_interactions(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    int dimension,
    bool is_symmetric,
    bool is_hermitian,
    KernelType* kernel,
    const CoordType global_bounds[6]) {
    
    int num_children = (dimension == 2) ? 4 : 8;
    
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
    int64_t local_morton_end = child_level.local_boxes.back().morton_index / num_children;
    
    // ===== Step 1: Initialize parent boxes =====
    
    initialize_local_boxes(
        parent_boxes,
        local_morton_start,
        parent_level_num,
        dimension,
        global_bounds);

    compute_neighbor_lists(parent_boxes, dimension, parent_level_num);    
    
    // ===== Step 2: Check boundary conditions =====
    // if no reduction happen, then no multiplier needed, otherwise multiply by 4 (2D) or 8 (3D)
    int transition_multiplier = (parent_level.num_active_processes == child_level.num_active_processes) ? 1 : num_children;
    int64_t start_accounted_for_reduction = parent_level.rank_to_morton[child_level.parent_level_owner] * num_parent_boxes * transition_multiplier;
    int64_t end_accounted_for_reduction = (parent_level.rank_to_morton[child_level.parent_level_owner] + 1) * (num_parent_boxes) * transition_multiplier - 1;
    for (auto& parent_box : parent_boxes) {
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
    for (int64_t parent_idx = 0; parent_idx < num_parent_boxes; ++parent_idx) {
        auto& parent_box = parent_boxes[parent_idx];
        int64_t first_child_idx = parent_idx * num_children;
        
        int64_t total_skeleton_points = 0;
        for (int c = 0; c < num_children; ++c) {
            auto& child = child_level.local_boxes[first_child_idx + c];
            total_skeleton_points += child.skeleton_indices.size();
        }
        
        parent_box.point_indices.reserve(total_skeleton_points);
        parent_box.point_coords.reserve(total_skeleton_points * dimension);
        
        for (int c = 0; c < num_children; ++c) {
            auto& child = child_level.local_boxes[first_child_idx + c];
            
            parent_box.point_indices.insert(
                parent_box.point_indices.end(),
                child.skeleton_indices.begin(),
                child.skeleton_indices.end());
            
            for (int64_t skel_idx : child.skeleton_indices) {
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
    
    // ===== Step 4: Build parent-level modified interactions =====
    
    // Helper struct to store child info (either from ghost or assisting boxes)
    struct ChildInfo {
        BoxData<CoordType, DataType>* box_ptr;  // nullptr if assisting box
        std::vector<CoordType> coords_ptr;
        int64_t num_points;
        bool is_ghost;  // true if from ghost_boxes, false if from assisting_boxes
    };
    
    for (size_t b1_idx = 0; b1_idx < parent_boxes.size(); ++b1_idx) {
        auto& B1 = parent_boxes[b1_idx];
        // printf("Assembling interactions for parent box %lu (morton %lu)\n", 
        //        b1_idx, B1.morton_index);
        // fflush(stdout);
        // Get relevant neighbors: [self, 1-hop]
        std::vector<uint64_t> relevant_neighbors;
        relevant_neighbors.push_back(B1.morton_index);
        
        std::vector<uint64_t> one_hop = (dimension == 2) ?
            morton::neighbors_2d(B1.morton_index, grid_size) :
            morton::neighbors_3d(B1.morton_index, grid_size);
        relevant_neighbors.insert(relevant_neighbors.end(), one_hop.begin(), one_hop.end());
        
        for (uint64_t neighbor_morton : relevant_neighbors) {
            
            
            
            // Check if neighbor is on-process or off-process
            bool is_on_process = (neighbor_morton >= local_morton_start && 
                                 neighbor_morton <= local_morton_end);
            
            // ===== DIAGONAL CASE =====
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
                            &child_i, &child_j, child_level, dimension, kernel);
                        
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
                                &child_i, &child_j, child_level, dimension, kernel);
                            
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
                    
                    // Store B1's view (A_NS)
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = B2.morton_index;
                    block_b1.A_NS.set_owned(
                        total_cols, total_rows, std::move(I_transposed), MatrixStorage<DataType>::FULL);
                    
                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map[B2.morton_index] = block_idx_b1;
                    
                    // Store B2's view (A_NS)
                    ModifiedBlock<DataType> block_b2;
                    block_b2.neighbor_morton = B1.morton_index;
                    block_b2.A_NS.set_owned(
                        total_rows, total_cols, std::move(I_B1_B2), MatrixStorage<DataType>::FULL);
                    
                    int64_t block_idx_b2 = B2.near_field_modified_interactions.size();
                    B2.near_field_modified_interactions.push_back(std::move(block_b2));
                    B2.near_field_interaction_map[B1.morton_index] = block_idx_b2;
                    
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
                            info.coords_ptr = ghost_box.point_coords;
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
                            info.coords_ptr = extract_skeleton_coords(assist_box, dimension);
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
                                    &child_i, child_j_info.box_ptr, child_level, dimension, kernel);
                                
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
                                    child_j_info.coords_ptr,
                                    n_i, 
                                    n_j,
                                    child_level,
                                    dimension,
                                    kernel,
                                    true  // TRANSPOSE: child_i.A_NS[child_j] is (n_j × n_i), need (n_i × n_j)
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
                    block_b1.A_NS.set_owned(
                        total_cols, total_rows, std::move(I_transposed), MatrixStorage<DataType>::FULL);
                    
                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map[neighbor_morton] = block_idx_b1;
                }
                
            } else {
                // ========== NONSYMMETRIC CASE ==========
                
                if (is_on_process) {
                    // Both on process - always process (no triangular optimization)
                    
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
                                &child_i, &child_j, child_level, dimension, kernel);
                            
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
                    
                    // Transpose for A_NS (B2 × B1)
                    std::vector<DataType> I_transposed(total_rows * total_cols);
                    for (int64_t i = 0; i < total_rows; ++i) {
                        for (int64_t j = 0; j < total_cols; ++j) {
                            I_transposed[j + i * total_cols] = I_B1_B2[i + j * total_rows];
                        }
                    }
                    
                    // Store A_NS for B1
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = B2.morton_index;
                    block_b1.A_NS.allocate(total_cols, total_rows, MatrixStorage<DataType>::FULL);
                    block_b1.A_NS.data = I_transposed;
                    
                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map_nonsymmetry[B2.morton_index] = block_idx_b1;
                    
                    // Store A_SN for B2
                    ModifiedBlock<DataType> block_b2;
                    block_b2.neighbor_morton = B1.morton_index;
                    block_b2.A_SN.set_owned(
                        total_cols, total_rows, std::move(I_transposed), MatrixStorage<DataType>::FULL);
                    
                    int64_t block_idx_b2 = B2.near_field_modified_interactions.size();
                    B2.near_field_modified_interactions.push_back(std::move(block_b2));
                    B2.near_field_interaction_map_nonsymmetry[B1.morton_index] = block_idx_b2;
                    
                } else {
                    // B2 is off-process (ghost or assisting) - store both A_NS and A_SN for B1
                    
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
                            info.coords_ptr = ghost_box.point_coords;
                            info.num_points = ghost_box.skeleton_indices.size();
                            info.is_ghost = true;
                            b2_children.push_back(info);
                            continue;
                        }
                        
                        // Check assisting boxes
                        auto assist_it = child_level.assisting_box_points_for_kernel_evaluation.find(child_morton);
                        if (assist_it != child_level.assisting_box_points_for_kernel_evaluation.end()) {
                            auto& assist_box = child_level.assisting_boxes[assist_it->second];
                            ChildInfo info;
                            info.box_ptr = nullptr;
                            info.coords_ptr = extract_skeleton_coords(assist_box, dimension);
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
                    int64_t first_child_b1 = b1_idx * num_children;
                    
                    // ===== Assemble I(B1, B2) for A_SN =====
                    std::vector<DataType> I_B1_B2(total_rows * total_cols);
                    
                    int64_t row_offset = 0;
                    for (int ci = 0; ci < num_children; ++ci) {
                        auto& child_i = child_level.local_boxes[first_child_b1 + ci];
                        int64_t n_i = child_i.skeleton_indices.size();
                        
                        // Extract skeleton coords from child_i
                        std::vector<CoordType> child_i_skel_coords(n_i * dimension);
                        for (int64_t idx = 0; idx < n_i; ++idx) {
                            int64_t skel_idx = child_i.skeleton_indices[idx];
                            for (int d = 0; d < dimension; ++d) {
                                child_i_skel_coords[idx * dimension + d] = 
                                    child_i.point_coords[skel_idx * dimension + d];
                            }
                        }
                        
                        int64_t col_offset = 0;
                        for (int cj = 0; cj < num_children; ++cj) {
                            const auto& child_j_info = b2_children[cj];
                            int64_t n_j = child_j_info.num_points;
                            
                            std::vector<DataType> C_block;
                            
                            if (child_j_info.is_ghost) {
                                C_block = extract_child_interaction(
                                    &child_i, child_j_info.box_ptr, child_level, dimension, kernel);
                            } else {
                                // // Direct kernel evaluation
                                // C_block.resize(n_i * n_j);
                                // kernel->evaluate_block(
                                //     child_i_skel_coords.data(), n_i,
                                //     child_j_info.coords_ptr.data(), n_j,
                                //     C_block.data(), n_i);
                                // Assisting box - try to find in child_i's modified interactions first
                                int64_t child_j_morton = neighbor_morton * num_children + cj;
                                C_block = extract_or_evaluate_child_interaction_for_assisting(
                                    &child_i,           // Looking in child_i's storage
                                    child_j_morton,     // For neighbor child_j
                                    child_j_info.coords_ptr,
                                    n_i, 
                                    n_j,
                                    child_level,
                                    dimension,
                                    kernel,
                                    true  // TRANSPOSE: child_i.A_NS[child_j] is (n_j × n_i), need (n_i × n_j)
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
                    
                    // ===== Assemble I(B2, B1) for A_NS =====
                    std::vector<DataType> I_B2_B1(total_cols * total_rows);
                    
                    row_offset = 0;
                    for (int ci = 0; ci < num_children; ++ci) {
                        const auto& child_i_info = b2_children[ci];
                        int64_t n_i = child_i_info.num_points;
                        
                        int64_t col_offset = 0;
                        for (int cj = 0; cj < num_children; ++cj) {
                            auto& child_j = child_level.local_boxes[first_child_b1 + cj];
                            int64_t n_j = child_j.skeleton_indices.size();
                            
                            // Extract skeleton coords from child_j
                            std::vector<CoordType> child_j_skel_coords(n_j * dimension);
                            for (int64_t idx = 0; idx < n_j; ++idx) {
                                int64_t skel_idx = child_j.skeleton_indices[idx];
                                for (int d = 0; d < dimension; ++d) {
                                    child_j_skel_coords[idx * dimension + d] = 
                                        child_j.point_coords[skel_idx * dimension + d];
                                }
                            }
                            
                            std::vector<DataType> C_block;
                            
                            if (child_i_info.is_ghost) {
                                C_block = extract_child_interaction(
                                    child_i_info.box_ptr, &child_j, child_level, dimension, kernel);
                            } else {
                                // // Direct kernel evaluation
                                // C_block.resize(n_i * n_j);
                                // kernel->evaluate_block(
                                //     child_i_info.coords_ptr.data(), n_i,
                                //     child_j_skel_coords.data(), n_j,
                                //     C_block.data(), n_i);
                                // Assisting box - try to find in child_j's modified interactions first
                                int64_t child_i_morton = neighbor_morton * num_children + ci;
                                
                                // Looking in child_j for child_i
                                // child_j.A_NS[child_i] is (n_i × n_j)
                                // Function call: extract from child_j (source), looking for child_i (target)
                                // Returns: (n_source × n_target) = (n_j × n_i) after transpose
                                C_block = extract_or_evaluate_child_interaction_for_assisting(
                                    &child_j,           // Looking in child_j's storage
                                    child_i_morton,     // For neighbor child_i
                                    child_i_info.coords_ptr,
                                    n_j,                // n_source = child_j skeleton size
                                    n_i,                // n_target = child_i skeleton size  
                                    child_level,
                                    dimension,
                                    kernel,
                                    true  // TRANSPOSE: child_j.A_NS[child_i] is (n_i × n_j), transpose to (n_j × n_i)
                                );
                                
                                // C_block is now (n_j × n_i), but we need (n_i × n_j) for I(B2, B1)
                                // Transpose again
                                std::vector<DataType> C_block_transposed(n_i * n_j);
                                for (int64_t i = 0; i < n_j; ++i) {
                                    for (int64_t j = 0; j < n_i; ++j) {
                                        C_block_transposed[j + i * n_i] = C_block[i + j * n_j];
                                    }
                                }
                                C_block = std::move(C_block_transposed);
                            }
                            
                            for (int64_t col = 0; col < n_j; ++col) {
                                for (int64_t row = 0; row < n_i; ++row) {
                                    I_B2_B1[(row_offset + row) + (col_offset + col) * total_cols] = 
                                        C_block[row + col * n_i];
                                }
                            }
                            col_offset += n_j;
                        }
                        row_offset += n_i;
                    }
                    
                    // Transpose I_B2_B1 for A_NS (B1 × B2)
                    std::vector<DataType> I_B2_B1_transposed(total_cols * total_rows);
                    for (int64_t i = 0; i < total_cols; ++i) {
                        for (int64_t j = 0; j < total_rows; ++j) {
                            I_B2_B1_transposed[j + i * total_rows] = I_B2_B1[i + j * total_cols];
                        }
                    }
                    
                    // Store both blocks for B1
                    ModifiedBlock<DataType> block_b1;
                    block_b1.neighbor_morton = neighbor_morton;
                    
                    // A_NS (B1 × B2)
                    block_b1.A_NS.set_owned(
                        total_rows, total_cols, std::move(I_B2_B1_transposed), MatrixStorage<DataType>::FULL);
                    
                    // A_SN (B1 × B2)
                    block_b1.A_SN.set_owned(
                        total_rows, total_cols, std::move(I_B1_B2), MatrixStorage<DataType>::FULL);
                    
                    int64_t block_idx_b1 = B1.near_field_modified_interactions.size();
                    B1.near_field_modified_interactions.push_back(std::move(block_b1));
                    B1.near_field_interaction_map_nonsymmetry[neighbor_morton] = block_idx_b1;
                }
            }
        }
    }
    
    return parent_boxes;
}


} // namespace fmm

#endif // FACTORIZATION_HPP



// /**
//  * @brief Update neighbor slicing in modified interactions after level-wide ID
//  * 
//  * After all boxes at a level have been processed with ID, some neighbors may
//  * have been compressed to skeleton DOFs. This function updates all modified
//  * interaction blocks to reflect the correct dimensions in BOTH directions.
//  * 
//  * Key insight: Both S direction (current box's skeleton) and N direction 
//  * (neighbor's skeleton) may need slicing, depending on when blocks were created
//  * relative to when boxes were eliminated.
//  * 
//  * Scenarios:
//  * 1. Far-field block created before current box eliminated: S needs slicing
//  * 2. Block created before neighbor eliminated: N needs slicing
//  * 3. Both: Need to slice both dimensions
//  * 
//  * @param level Tree level to process
//  * @param is_symmetric Whether matrix is symmetric
//  */
// template<typename CoordType, typename DataType>
// void update_neighbor_slicing_for_level(
//     TreeLevel<CoordType, DataType>& level,
//     bool is_symmetric) {
    
//     // Helper lambda to slice a matrix in both S and N dimensions if needed
//     auto slice_both_directions = [&](
//         MatrixStorage<DataType>& matrix,
//         int64_t neighbor_morton,
//         int64_t current_box_skeleton_size,
//         int64_t current_box_morton,
//         bool is_A_NS) {
        
//         if (!matrix.is_allocated()) {
//             return;  // Nothing to slice
//         }
        
//         int64_t stored_rows = matrix.rows;
//         int64_t stored_cols = matrix.cols;
        
//         // Find neighbor box to get its skeleton size
//         BoxData<CoordType, DataType>* neighbor_box = level.find_local_box(neighbor_morton);
//         if (neighbor_box == nullptr) {
//             auto ghost_it = level.ghost_id_to_index.find(neighbor_morton);
//             if (ghost_it != level.ghost_id_to_index.end()) {
//                 neighbor_box = &level.ghost_boxes[ghost_it->second];
//             }
//         }
        
//         if (neighbor_box == nullptr) {
//             // Neighbor not on this process - skip
//             return;
//         }
        
//         // Check elimination status
//         bool neighbor_eliminated = (level.eliminated_boxes.find(neighbor_morton) != 
//                                    level.eliminated_boxes.end());
//         bool has_skeleton = !neighbor_box->skeleton_indices.empty();
        
//         // ERROR CHECKING
//         if (neighbor_eliminated && !has_skeleton) {
//             throw std::runtime_error(
//                 "update_neighbor_slicing: Box " + std::to_string(neighbor_morton) +
//                 " is marked as eliminated but has no skeleton indices!");
//         }
        
//         // Determine target dimensions
//         int64_t target_rows, target_cols;
//         bool need_slice_S = false;
//         bool need_slice_N = false;
        
//         if (is_A_NS) {
//             // A_NS: (n_neighbor × k_box)
//             target_rows = has_skeleton ? neighbor_box->skeleton_indices.size() : stored_rows;
//             target_cols = current_box_skeleton_size;
            
//             need_slice_N = (stored_rows != target_rows) && has_skeleton;
//             need_slice_S = (stored_cols != target_cols);
            
//         } else {
//             // A_SN: (k_box × n_neighbor)
//             target_rows = current_box_skeleton_size;
//             target_cols = has_skeleton ? neighbor_box->skeleton_indices.size() : stored_cols;
            
//             need_slice_S = (stored_rows != target_rows);
//             need_slice_N = (stored_cols != target_cols) && has_skeleton;
//         }
        
//         if (!need_slice_S && !need_slice_N) {
//             // Already correct size
//             return;
//         }
        
//         // Perform slicing based on what's needed
//         std::vector<DataType> sliced(target_rows * target_cols);
        
//         if (is_A_NS) {
//             // A_NS: (n_neighbor × k_box)
            
//             // Build row indices (N dimension)
//             std::vector<int64_t> row_indices;
//             if (need_slice_N) {
//                 row_indices = neighbor_box->skeleton_indices;
//             } else {
//                 row_indices.resize(stored_rows);
//                 std::iota(row_indices.begin(), row_indices.end(), 0);
//             }
            
//             // Build column indices (S dimension)
//             std::vector<int64_t> col_indices;
//             if (need_slice_S) {
//                 // Get current box to retrieve its skeleton indices
//                 BoxData<CoordType, DataType>* current_box = level.find_local_box(current_box_morton);
//                 if (current_box == nullptr) {
//                     auto ghost_it = level.ghost_id_to_index.find(current_box_morton);
//                     if (ghost_it != level.ghost_id_to_index.end()) {
//                         current_box = &level.ghost_boxes[ghost_it->second];
//                     }
//                 }
                
//                 if (current_box == nullptr || current_box->skeleton_indices.empty()) {
//                     throw std::runtime_error(
//                         "update_neighbor_slicing: Cannot find current box " + 
//                         std::to_string(current_box_morton) + " or it has no skeleton");
//                 }
                
//                 col_indices = current_box->skeleton_indices;
//             } else {
//                 col_indices.resize(stored_cols);
//                 std::iota(col_indices.begin(), col_indices.end(), 0);
//             }
            
//             // Extract submatrix
//             extract_submatrix(
//                 matrix.data.data(), stored_rows,
//                 row_indices,
//                 col_indices,
//                 sliced.data(), target_rows
//             );
            
//         } else {
//             // A_SN: (k_box × n_neighbor)
            
//             // Build row indices (S dimension)
//             std::vector<int64_t> row_indices;
//             if (need_slice_S) {
//                 BoxData<CoordType, DataType>* current_box = level.find_local_box(current_box_morton);
//                 if (current_box == nullptr) {
//                     auto ghost_it = level.ghost_id_to_index.find(current_box_morton);
//                     if (ghost_it != level.ghost_id_to_index.end()) {
//                         current_box = &level.ghost_boxes[ghost_it->second];
//                     }
//                 }
                
//                 if (current_box == nullptr || current_box->skeleton_indices.empty()) {
//                     throw std::runtime_error(
//                         "update_neighbor_slicing: Cannot find current box for S slicing");
//                 }
                
//                 row_indices = current_box->skeleton_indices;
//             } else {
//                 row_indices.resize(stored_rows);
//                 std::iota(row_indices.begin(), row_indices.end(), 0);
//             }
            
//             // Build column indices (N dimension)
//             std::vector<int64_t> col_indices;
//             if (need_slice_N) {
//                 col_indices = neighbor_box->skeleton_indices;
//             } else {
//                 col_indices.resize(stored_cols);
//                 std::iota(col_indices.begin(), col_indices.end(), 0);
//             }
            
//             // Extract submatrix
//             extract_submatrix(
//                 matrix.data.data(), stored_rows,
//                 row_indices,
//                 col_indices,
//                 sliced.data(), target_rows
//             );
//         }
        
//         // Update matrix in-place
//         matrix.data = std::move(sliced);
//         matrix.rows = target_rows;
//         matrix.cols = target_cols;
//         matrix.lda = target_rows;
//     };
    
//     // Process all local boxes
//     for (auto& box : level.local_boxes) {
        
//         int64_t k_box = box.skeleton_indices.size();
        
//         if (k_box == 0) {
//             // Box not yet processed with ID - skip
//             continue;
//         }
        
//         // Update near-field blocks
//         for (auto& modified_block : box.near_field_modified_interactions) {
//             int64_t neighbor_morton = modified_block.neighbor_morton;
            
//             if (modified_block.A_NS.is_allocated()) {
//                 slice_both_directions(
//                     modified_block.A_NS, neighbor_morton, 
//                     k_box, box.morton_index, true);
//             }
            
//             if (!is_symmetric && modified_block.A_SN.is_allocated()) {
//                 slice_both_directions(
//                     modified_block.A_SN, neighbor_morton, 
//                     k_box, box.morton_index, false);
//             }
//         }
        
//         // Update far-field blocks (this is critical for 2-hop neighbors!)
//         for (auto& modified_block : box.far_field_modified_interactions) {
//             int64_t neighbor_morton = modified_block.neighbor_morton;
            
//             if (modified_block.A_NS.is_allocated()) {
//                 slice_both_directions(
//                     modified_block.A_NS, neighbor_morton, 
//                     k_box, box.morton_index, true);
//             }
            
//             if (!is_symmetric && modified_block.A_SN.is_allocated()) {
//                 slice_both_directions(
//                     modified_block.A_SN, neighbor_morton, 
//                     k_box, box.morton_index, false);
//             }
//         }
//     }
// }
