#ifndef ID_DECOMPOSITION_HPP
#define ID_DECOMPOSITION_HPP

#include <vector>
#include <random>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "tree.hpp"
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <chrono>
#include "blas_declare.hpp"
#include <cassert>

/*used blas functions*/
// dgeqp3rk_, dgemm_, dtrtrs_, daxpy_

namespace fmm {


/**
 * @brief Interpolative decomposition result structure
 */
template<typename DataType>
struct IDResult {
    std::vector<int64_t> skeleton_indices;    ///< Skeleton DOF indices (size k)
    std::vector<int64_t> redundant_indices;   ///< Redundant DOF indices (size n-k)
    MatrixStorage<DataType> interpolation;    ///< Interpolation matrix T (k × n-k)
    int64_t rank;                             ///< Numerical rank k
    
    IDResult() : rank(0) {}
};



/**
 * @brief Perform interpolative decomposition using RRQR
 * 
 * Given matrix A (m × n), finds a partition of columns into:
 * - Skeleton columns S (indices 1:k)
 * - Redundant columns R (indices k+1:n)
 * 
 * Such that: A[:, R] ≈ A[:, S] * T
 * 
 * Implementation uses LAPACK's dgeqp3rk (rank-revealing QR with stopping criteria).
 * 
 * @param A Input matrix (m × n, column-major, will be overwritten)
 * @param m Number of rows
 * @param n Number of columns
 * @param lda Leading dimension of A
 * @param tolerance Relative tolerance for rank determination
 * @param max_rank Maximum rank (0 = no limit, use min(m,n))
 * @return IDResult containing skeleton/redundant partitions and interpolation matrix
 */
template<typename DataType>
IDResult<DataType> compute_id(
    DataType* A,
    int64_t m,
    int64_t n,
    int64_t lda,
    double tolerance = 1e-6,
    int64_t max_rank = 0) {
    
    static_assert(std::is_same_v<DataType, double> || 
                  std::is_same_v<DataType, std::complex<double>>,
                  "Only double and complex<double> supported");
    
    if (m <= 0 || n <= 0 || lda < m) {
        throw std::invalid_argument("Invalid matrix dimensions for ID");
    }
    
    IDResult<DataType> result;
    
    // Setup LAPACK parameters
    int M = static_cast<int>(m);
    int N = static_cast<int>(n);
    int NRHS = 0;
    int LDA = static_cast<int>(lda);
    
    // Stopping criteria
    int KMAX = (max_rank > 0) ? static_cast<int>(max_rank) : std::min(M, N);
    double ABSTOL = -1.0;
    double RELTOL = tolerance;
    
    // Output variables
    int K = 0;
    double MAXC2NRMK = 0.0;
    double RELMAXC2NRMK = 0.0;
    
    // Permutation and Householder reflectors
    std::vector<int> JPIV(N);
    std::vector<DataType> TAU(std::min(M, N));
    
    // Workspace query
    int LWORK = -1;
    std::vector<DataType> WORK(1);
    std::vector<int> IWORK(N);
    int INFO = 0;
    
    if constexpr (std::is_same_v<DataType, double>) {
        // Query optimal workspace size
        dgeqp3rk_(&M, &N, &NRHS, &KMAX, &ABSTOL, &RELTOL,
                  A, &LDA, &K, &MAXC2NRMK, &RELMAXC2NRMK,
                  JPIV.data(), TAU.data(),
                  WORK.data(), &LWORK, IWORK.data(), &INFO);
        
        if (INFO != 0) {
            throw std::runtime_error("RRQR workspace query failed with INFO = " + 
                                   std::to_string(INFO));
        }
        
        LWORK = static_cast<int>(WORK[0]);
        WORK.resize(LWORK);
        
        // Perform factorization
        dgeqp3rk_(&M, &N, &NRHS, &KMAX, &ABSTOL, &RELTOL,
                  A, &LDA, &K, &MAXC2NRMK, &RELMAXC2NRMK,
                  JPIV.data(), TAU.data(),
                  WORK.data(), &LWORK, IWORK.data(), &INFO);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        // Complex version needs separate RWORK array
        std::vector<double> RWORK(2 * N);
        
        // Query optimal workspace size
        zgeqp3rk_(&M, &N, &NRHS, &KMAX, &ABSTOL, &RELTOL,
                  A, &LDA, &K, &MAXC2NRMK, &RELMAXC2NRMK,
                  JPIV.data(), TAU.data(),
                  WORK.data(), &LWORK,
                  RWORK.data(), IWORK.data(), &INFO);
        
        if (INFO != 0) {
            throw std::runtime_error("RRQR workspace query failed with INFO = " + 
                                   std::to_string(INFO));
        }
        
        LWORK = static_cast<int>(WORK[0].real());
        WORK.resize(LWORK);
        
        // Perform factorization
        zgeqp3rk_(&M, &N, &NRHS, &KMAX, &ABSTOL, &RELTOL,
                  A, &LDA, &K, &MAXC2NRMK, &RELMAXC2NRMK,
                  JPIV.data(), TAU.data(),
                  WORK.data(), &LWORK,
                  RWORK.data(), IWORK.data(), &INFO);
    }
    
    if (INFO != 0) {
        throw std::runtime_error("RRQR factorization failed with INFO = " + 
                               std::to_string(INFO));
    }
    
    result.rank = K;
    
    // Handle edge cases
    if (K == 0) {
        for (int64_t i = 0; i < n; ++i) {
            result.redundant_indices.push_back(i);
        }
        result.interpolation.allocate(0, n);
        return result;
    }
    
    if (K == n) {
        for (int64_t i = 0; i < n; ++i) {
            result.skeleton_indices.push_back(i);
        }
        result.interpolation.allocate(n, 0);
        return result;
    }
    
    // Extract skeleton and redundant indices from permutation
    for (int i = 0; i < K; ++i) {
        result.skeleton_indices.push_back(JPIV[i] - 1);
    }
    
    for (int i = K; i < N; ++i) {
        result.redundant_indices.push_back(JPIV[i] - 1);
    }
    
    // Compute interpolation matrix T = R11^{-1} * R12
    int64_t num_redundant = N - K;
    result.interpolation.allocate(K, num_redundant);
    
    // Extract R12 from the factored matrix A
    for (int64_t j = 0; j < num_redundant; ++j) {
        for (int64_t i = 0; i < K; ++i) {
            result.interpolation(i, j) = A[i + (K + j) * LDA];
        }
    }
    
    // Solve R11 * T = R12 for T
    char UPLO = 'U';
    char TRANS = 'N';
    char DIAG = 'N';
    int NK = static_cast<int>(K);
    int NRHS_solve = static_cast<int>(num_redundant);
    int LDB = NK;
    
    if constexpr (std::is_same_v<DataType, double>) {
        dtrtrs_(&UPLO, &TRANS, &DIAG,
                &NK, &NRHS_solve,
                A, &LDA,
                result.interpolation.data.data(), &LDB,
                &INFO);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        ztrtrs_(&UPLO, &TRANS, &DIAG,
                &NK, &NRHS_solve,
                A, &LDA,
                result.interpolation.data.data(), &LDB,
                &INFO);
    }
    
    if (INFO != 0) {
        throw std::runtime_error("Triangular solve for interpolation matrix failed with INFO = " + 
                               std::to_string(INFO));
    }
    
    return result;
}


template<typename DataType>
IDResult<DataType> compute_id_complex(
    DataType* A,
    int64_t m,
    int64_t n,
    int64_t lda,
    double tolerance = 1e-6,
    int64_t max_rank = 0) {
    
    static_assert(std::is_same_v<DataType, double> || 
                  std::is_same_v<DataType, std::complex<double>>,
                  "Only double and complex<double> supported");
    
    if (m <= 0 || n <= 0 || lda < m) {
        throw std::invalid_argument("Invalid matrix dimensions for ID");
    }
    
    IDResult<DataType> result;
    
    // Setup LAPACK parameters
    int M = static_cast<int>(m);
    int N = static_cast<int>(n);
    int LDA = static_cast<int>(lda);
    int minMN = std::min(M, N);
    
    // Permutation and Householder reflectors
    std::vector<int> JPIV(N, 0);  // Initialize to 0: all columns free
    std::vector<DataType> TAU(minMN);
    
    // Workspace query
    int LWORK = -1;
    std::vector<DataType> WORK(1);
    int INFO = 0;
    
    if constexpr (std::is_same_v<DataType, double>) {
        dgeqp3_(&M, &N, A, &LDA,
                JPIV.data(), TAU.data(),
                WORK.data(), &LWORK, &INFO);
        
        if (INFO != 0) {
            throw std::runtime_error("RRQR workspace query failed with INFO = " + 
                                   std::to_string(INFO));
        }
        
        LWORK = static_cast<int>(WORK[0]);
        WORK.resize(LWORK);
        
        dgeqp3_(&M, &N, A, &LDA,
                JPIV.data(), TAU.data(),
                WORK.data(), &LWORK, &INFO);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        std::vector<double> RWORK(2 * N);
        
        zgeqp3_(&M, &N, A, &LDA,
                JPIV.data(), TAU.data(),
                WORK.data(), &LWORK,
                RWORK.data(), &INFO);
        
        if (INFO != 0) {
            throw std::runtime_error("RRQR workspace query failed with INFO = " + 
                                   std::to_string(INFO));
        }
        
        LWORK = static_cast<int>(WORK[0].real());
        WORK.resize(LWORK);
        
        zgeqp3_(&M, &N, A, &LDA,
                JPIV.data(), TAU.data(),
                WORK.data(), &LWORK,
                RWORK.data(), &INFO);
    }
    
    if (INFO != 0) {
        throw std::runtime_error("RRQR factorization failed with INFO = " + 
                               std::to_string(INFO));
    }
    
    // Determine rank from diagonal of R
    // |R(i,i)| <= tolerance * |R(1,1)| means we stop
    double R11 = std::abs(A[0]);  // |R(1,1)|
    int KMAX = (max_rank > 0) ? static_cast<int>(max_rank) : minMN;
    int K = 0;
    
    if (R11 > 0.0) {
        for (int i = 0; i < std::min(KMAX, minMN); ++i) {
            if (std::abs(A[i + i * LDA]) > tolerance * R11) {
                K++;
            } else {
                break;
            }
        }
    }
    
    result.rank = K;
    
    // Handle edge cases
    if (K == 0) {
        for (int64_t i = 0; i < n; ++i) {
            result.redundant_indices.push_back(i);
        }
        result.interpolation.allocate(0, n);
        return result;
    }
    
    if (K == n) {
        for (int64_t i = 0; i < n; ++i) {
            result.skeleton_indices.push_back(i);
        }
        result.interpolation.allocate(n, 0);
        return result;
    }
    
    // Extract skeleton and redundant indices from permutation
    for (int i = 0; i < K; ++i) {
        result.skeleton_indices.push_back(JPIV[i] - 1);
    }
    
    for (int i = K; i < N; ++i) {
        result.redundant_indices.push_back(JPIV[i] - 1);
    }
    
    // Compute interpolation matrix T = R11^{-1} * R12
    int64_t num_redundant = N - K;
    result.interpolation.allocate(K, num_redundant);
    
    // Extract R12 from the factored matrix A
    for (int64_t j = 0; j < num_redundant; ++j) {
        for (int64_t i = 0; i < K; ++i) {
            result.interpolation(i, j) = A[i + (K + j) * LDA];
        }
    }
    
    // Solve R11 * T = R12 for T
    char UPLO = 'U';
    char TRANS = 'N';
    char DIAG = 'N';
    int NK = static_cast<int>(K);
    int NRHS_solve = static_cast<int>(num_redundant);
    int LDB = NK;
    
    if constexpr (std::is_same_v<DataType, double>) {
        dtrtrs_(&UPLO, &TRANS, &DIAG,
                &NK, &NRHS_solve,
                A, &LDA,
                result.interpolation.data.data(), &LDB,
                &INFO);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        ztrtrs_(&UPLO, &TRANS, &DIAG,
                &NK, &NRHS_solve,
                A, &LDA,
                result.interpolation.data.data(), &LDB,
                &INFO);
    }
    
    if (INFO != 0) {
        throw std::runtime_error("Triangular solve for interpolation matrix failed with INFO = " + 
                               std::to_string(INFO));
    }
    
    return result;
}



/**
 * @brief Generate sparse sketch of matrix using stratified random projection
 * 
 * Computes B = S * A where S is a (d × m) sparse random matrix.
 * Each column of S has exactly k nonzeros at stratified random locations.
 * 
 * The k nonzeros per column are placed in disjoint blocks:
 * - Block 0: [0, d/k)
 * - Block 1: [d/k, 2*d/k)
 * - ...
 * - Block k-1: [(k-1)*d/k, d)
 * 
 * This ensures good coverage of the sketch space.
 * 
 * Never forms S explicitly! Uses BLAS AXPY for efficiency.
 * Complexity: O(k * m * n) where typically k = 3 or 4.
 * 
 * @param A Input matrix (m × n, column-major)
 * @param m Number of rows in A
 * @param n Number of columns in A
 * @param lda Leading dimension of A
 * @param B Output sketch matrix (d × n, column-major, must be zero-initialized)
 * @param sketch_size d (number of rows in sketch)
 * @param ldb Leading dimension of B
 * @param k Number of nonzeros per column of S (typically 3 or 4)
 * @param seed Random seed for reproducibility
 */
// inline void sketch_sparse_random(
//     const double* A,
//     int64_t m,
//     int64_t n,
//     int64_t lda,
//     double* B,
//     int64_t sketch_size,
//     int64_t ldb,
//     int k,
//     uint64_t seed) {
    
//     if (m <= 0 || n <= 0 || sketch_size <= 0 || k <= 0) {
//         throw std::invalid_argument("Invalid dimensions for sparse sketching");
//     }
//     if (lda < m) {
//         throw std::invalid_argument("lda must be >= m");
//     }
//     if (ldb < sketch_size) {
//         throw std::invalid_argument("ldb must be >= sketch_size");
//     }
//     if (k > sketch_size) {
//         throw std::invalid_argument("k cannot exceed sketch_size");
//     }
    
//     // Zero-initialize B
//     std::fill(B, B + ldb * n, 0.0);
    
//     int d = static_cast<int>(sketch_size);
//     int block_size = d / k;
    
//     // Single RNG - reuse throughout (no reseeding per column!)
//     std::mt19937_64 rng(seed);
    
//     // Precompute block boundaries to avoid recalculation
//     std::vector<int> block_starts(k);
//     std::vector<int> block_ranges(k);
//     for (int b = 0; b < k; ++b) {
//         block_starts[b] = b * block_size;
//         int end = (b == k - 1) ? d : (b + 1) * block_size;
//         block_ranges[b] = end - block_starts[b];
//     }
    
//     // Normalization factor
//     double norm_factor = 1.0 / std::sqrt(static_cast<double>(k));
    
//     // BLAS parameters (constant throughout)
//     int n_int = static_cast<int>(n);
//     int lda_int = static_cast<int>(lda);
//     int ldb_int = static_cast<int>(ldb);
//     auto id_start = std::chrono::high_resolution_clock::now();
//     // Process each row of A (i.e., each column of S)
//     for (int64_t j = 0; j < m; ++j) {
//         // Generate k stratified indices and signs
//         for (int b = 0; b < k; ++b) {
//             // Generate one random number and extract both index and sign from it
//             uint64_t rand_val = rng();
            
//             // Index in [start, start + range) using fast modulo
//             int idx = block_starts[b] + static_cast<int>(rand_val % block_ranges[b]);
            
//             // Sign from high bit: extract and map to {-norm_factor, +norm_factor}
//             double sign = (rand_val & 1) ? norm_factor : -norm_factor;
            
//             // B[idx, :] += sign * A[j, :]
//             daxpy_(&n_int, &sign, 
//                    const_cast<double*>(&A[j]), &lda_int,
//                    &B[idx], &ldb_int);
//         }
//     }
//     auto id_end = std::chrono::high_resolution_clock::now();
//     auto id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
//     std::cout << "  sketch time: " << id_duration.count() << " ms" << " sketch_size: " << sketch_size << " m: " << m << " n: " << n << std::endl;
// }

template<typename DataType>
inline void sketch_sparse_random(
    const DataType* A,
    int64_t m,
    int64_t n,
    int64_t lda,
    DataType* B,
    int64_t sketch_size,
    int64_t ldb,
    int k,
    uint64_t seed) {
    
    if (m <= 0 || n <= 0 || sketch_size <= 0 || k <= 0) {
        throw std::invalid_argument("Invalid dimensions for sparse sketching");
    }
    if (lda < m) {
        throw std::invalid_argument("lda must be >= m");
    }
    if (ldb < sketch_size) {
        throw std::invalid_argument("ldb must be >= sketch_size");
    }
    if (k > sketch_size) {
        throw std::invalid_argument("k cannot exceed sketch_size");
    }
    
    int d = static_cast<int>(sketch_size);
    
    // Allocate row-major temporary: temp[i, j] = temp[i * n + j]
    // temp is d × n in row-major layout
    std::vector<DataType> temp(d * n, 0.0);
    
    int block_size = d / k;
    
    // Single RNG - reuse throughout
    std::mt19937_64 rng(seed);
    
    // Precompute block boundaries
    std::vector<int> block_starts(k);
    std::vector<int> block_ranges(k);
    for (int b = 0; b < k; ++b) {
        block_starts[b] = b * block_size;
        int end = (b == k - 1) ? d : (b + 1) * block_size;
        block_ranges[b] = end - block_starts[b];
    }

    // Normalization factor
    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<
            std::is_same_v<DataType, std::complex<float>>, float,
            DataType
        >
    >;
    DataType norm_factor = DataType(1.0 / std::sqrt(static_cast<RealType>(k)));

    // BLAS parameters
    int n_int = static_cast<int>(n);
    int lda_int = static_cast<int>(lda);
    int inc_one = 1;

    // Process each row of A (i.e., each column of S)
    for (int64_t j = 0; j < m; ++j) {
        // Generate k stratified indices and signs
        for (int b = 0; b < k; ++b) {
            // Generate random index and sign
            uint64_t rand_val = rng();

            // Index in [start, start + range)
            int idx = block_starts[b] + static_cast<int>(rand_val % block_ranges[b]);

            // Sign from high bit
            DataType sign = ((rand_val >> 63) & 1) ? norm_factor : -norm_factor;

            // temp[idx, :] += sign * A[j, :]
            // Row idx of temp (row-major) starts at temp[idx * n]
            // Row j of A (column-major) starts at A[j] with stride lda
            if constexpr (std::is_same_v<DataType, double>) {
                daxpy_(&n_int, &sign,
                    const_cast<double*>(&A[j]), &lda_int,
                    &temp[idx * n], &inc_one);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zaxpy_(&n_int, &sign,
                    const_cast<std::complex<double>*>(&A[j]), &lda_int,
                    &temp[idx * n], &inc_one);
            }
        }
    }
    
    // Transpose temp (d × n, row-major) into B (d × n, column-major)
    // This is equivalent to copying temp as-is into B^T
    // Or: B = temp^T, but temp is row-major so just copy
    
    // Method 1: Manual transpose (simple and cache-friendly)
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < n_int; ++j) {
            B[i + j * ldb] = temp[i * n_int + j];
        }
    }

    // auto id_end = std::chrono::high_resolution_clock::now();
    // auto id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
    // std::cout << "  sketch time: " << id_duration.count() << " ms" << " sketch_size: " << sketch_size << " m: " << m << " n: " << n << std::endl;
    
    /* Alternative Method 2: Use BLAS dgemm for transpose
    // B = 1.0 * temp^T + 0.0 * B
    // But temp is stored row-major, so treat it as column-major transpose
    char trans_n = 'N', trans_t = 'T';
    int d_int = d;
    double alpha = 1.0, beta = 0.0;
    int ld_temp = n_int;  // leading dim of temp (row-major = # cols)
    int ldb_int = static_cast<int>(ldb);
    
    dgemm_(&trans_t, &trans_n,
           &d_int, &n_int, &n_int,
           &alpha, temp.data(), &ld_temp,
           temp.data(), &ld_temp,  // dummy for identity
           &beta, B, &ldb_int);
    */
}


/**
 * @brief Compute ID using sparse random sketching
 * 
 * For A (m × n) where m >> n:
 * 1. Sparse sketch: Y = S * A where S is (d × m) with k nonzeros per column
 * 2. Compute ID of Y to find skeleton columns
 * 3. Return skeleton/redundant indices (works for A)
 * 
 * Complexity: O(k * m * n + ID_cost) where k << d
 * Much faster than dense sketching when m is very large.
 * 
 * @param A Input matrix (m × n, column-major, NOT overwritten)
 * @param B Output sketch matrix (will be resized to sketch_size × n)
 * @param m Number of rows
 * @param n Number of columns
 * @param lda Leading dimension of A
 * @param rank_or_tol Rank or tolerance for ID
 * @param sketch_factor Sketch size = sketch_factor * n (default: 2.0)
 * @param k Sparsity parameter: nonzeros per column of S (default: 3)
 * @param seed Random seed (0 = auto-generate)
 * @return IDResult
 */
template<typename DataType>
IDResult<DataType> compute_id_sparse_sketch(
    const DataType* A,
    std::vector<DataType>& B,
    int64_t m,
    int64_t n,
    int64_t lda,
    double rank_or_tol = 1e-6,
    double sketch_factor = 2.0,
    int k = 3,
    uint64_t seed = 0) {
    
    static_assert(std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>,
                  "Only double precision supported currently");
    
    if (m <= 0 || n <= 0 || lda < m) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }
    if (sketch_factor < 1.0) {
        throw std::invalid_argument("sketch_factor must be >= 1.0");
    }
    if (k <= 0) {
        throw std::invalid_argument("k must be positive");
    }
    
    // Auto-generate seed if not provided
    if (seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    
    // Sketch size: d = sketch_factor * n
    int64_t sketch_size = static_cast<int64_t>(std::ceil(sketch_factor * n));
    sketch_size = std::min(sketch_size, m);  // Can't exceed m
    sketch_size = std::max(sketch_size, n);  // At least n
    
    // Ensure k doesn't exceed sketch_size
    k = std::min(k, static_cast<int>(sketch_size));
    
    // Allocate sketch matrix Y = S * A (d × n)
    B.resize(sketch_size * n);
    
    // Generate sparse sketch: Y = S * A
    
    sketch_sparse_random(A, m, n, lda, B.data(), sketch_size, sketch_size, k, seed);
   
    
    // Compute ID on sketched matrix Y
    // auto id_start = std::chrono::high_resolution_clock::now();
    if constexpr (std::is_same_v<DataType, double> || std::is_same_v<DataType, float>){
        // IDResult<DataType> result = compute_id(
        //     B.data(),
        //     sketch_size,
        //     n,
        //     sketch_size,
        //     rank_or_tol,
        //     0  // max_rank = 0 means auto-detect from tolerance
        // );
        IDResult<DataType> result = compute_id_complex(
            B.data(),
            sketch_size,
            n,
            sketch_size,
            rank_or_tol,
            0  // max_rank = 0 means auto-detect from tolerance
        );
        return result;
    }
    else if constexpr (std::is_same_v<DataType, std::complex<double>> || std::is_same_v<DataType, std::complex<float>>){
        IDResult<DataType> result = compute_id_complex(
            B.data(),
            sketch_size,
            n,
            sketch_size,
            rank_or_tol,
            0  // max_rank = 0 means auto-detect from tolerance
        );
        return result;
    }
    
    // auto id_end = std::chrono::high_resolution_clock::now();
    // auto id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
    // std::cout << "  id time: " << id_duration.count() << " ms" << std::endl;
    
    
}

/**
 * @brief Generate random sketch of matrix using uniform random projection
 * 
 * Computes B = S * A where S is a (sketch_size × m) random matrix with 
 * entries uniformly sampled from [-1, 1].
 * 
 * Uses BLAS dgemm for efficient multiplication: O(sketch_size * m * n)
 * 
 * @param A Input matrix (m × n, column-major)
 * @param m Number of rows in A
 * @param n Number of columns in A
 * @param lda Leading dimension of A
 * @param B Output sketch matrix (sketch_size × n, column-major, pre-allocated)
 * @param sketch_size Number of rows in sketch
 * @param ldb Leading dimension of B
 * @param seed Random seed for reproducibility
 */
inline void sketch_uniform_random(
    const double* A,
    int64_t m,
    int64_t n,
    int64_t lda,
    double* B,
    int64_t sketch_size,
    int64_t ldb,
    uint64_t seed) {
    
    if (m <= 0 || n <= 0 || sketch_size <= 0) {
        throw std::invalid_argument("Invalid dimensions for sketching");
    }
    if (lda < m) {
        throw std::invalid_argument("lda must be >= m");
    }
    if (ldb < sketch_size) {
        throw std::invalid_argument("ldb must be >= sketch_size");
    }
    if (sketch_size > m) {
        throw std::invalid_argument("Sketch size cannot exceed m");
    }

    throw std::invalid_argument("Mixed precision GEMM not implemented yet");
    
    // Generate random sketch matrix S (sketch_size × m)
    std::vector<double> S(sketch_size * m);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    
    for (int64_t i = 0; i < sketch_size * m; ++i) {
        S[i] = dist(rng);
    }
    
    // Normalize rows for better numerical properties
    double scale = 1.0 / std::sqrt(static_cast<double>(m));
    for (int64_t i = 0; i < sketch_size * m; ++i) {
        S[i] *= scale;
    }
    
    // Compute B = S * A using BLAS
    // B (sketch_size × n) = S (sketch_size × m) * A (m × n)
    char transa = 'N';  // No transpose of S
    char transb = 'N';  // No transpose of A
    int M_blas = static_cast<int>(sketch_size);
    int N_blas = static_cast<int>(n);
    int K_blas = static_cast<int>(m);
    double alpha = 1.0;
    double beta = 0.0;
    int LDA_blas = static_cast<int>(sketch_size);  // Leading dim of S
    int LDB_blas = static_cast<int>(lda);          // Leading dim of A
    int LDC_blas = static_cast<int>(ldb);          // Leading dim of B
    
    dgemm_(&transa, &transb,
           &M_blas, &N_blas, &K_blas,
           &alpha, S.data(), &LDA_blas,
           const_cast<double*>(A), &LDB_blas,
           &beta, B, &LDC_blas);
}

/**
 * @brief Compute ID using random sketching
 * 
 * For A (m × n) where m >> n:
 * 1. Sketch: Y = S * A where S is (d × m) with uniform(-1, 1) entries
 * 2. Compute ID of Y to find skeleton columns
 * 3. Return skeleton/redundant indices (works for A)
 * 
 * The sketch size d is chosen as d = sketch_factor * n, typically:
 * - sketch_factor = 2.0 gives d = 2n (standard)
 * - sketch_factor = 3.0 gives d = 3n (more conservative)
 * 
 * @param A Input matrix (m × n, column-major, NOT overwritten)
 * @param B Output sketch matrix (will be resized to sketch_size × n)
 * @param m Number of rows
 * @param n Number of columns
 * @param lda Leading dimension of A
 * @param rank_or_tol Rank or tolerance for ID
 * @param sketch_factor Sketch size = sketch_factor * n (default: 2.0)
 * @param seed Random seed (0 = auto-generate)
 * @return IDResult
 */
template<typename DataType>
IDResult<DataType> compute_id_uniform_sketch(
    const DataType* A,
    std::vector<DataType>& B,
    int64_t m,
    int64_t n,
    int64_t lda,
    double rank_or_tol = 1e-6,
    double sketch_factor = 2.0,
    uint64_t seed = 0) {
    
    static_assert(std::is_same_v<DataType, double>, 
                  "Only double precision supported");
    
    if (m <= 0 || n <= 0 || lda < m) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }
    if (sketch_factor < 1.0) {
        throw std::invalid_argument("sketch_factor must be >= 1.0");
    }
    
    // // Auto-generate seed if not provided
    // if (seed == 0) {
    //     std::random_device rd;
    //     seed = rd();
    // }
    
    // Sketch size: d = sketch_factor * n
    // This preserves the column space structure of A
    int64_t sketch_size = static_cast<int64_t>(std::ceil(sketch_factor * n));
    sketch_size = std::min(sketch_size, m);  // Can't exceed m
    sketch_size = std::max(sketch_size, n);  // At least n
    
    // Allocate sketch matrix Y = S * A (d × n)
    B.resize(sketch_size * n);
    
    // Generate sketch: Y = S * A where S is (sketch_size × m)
    auto id_start = std::chrono::high_resolution_clock::now();
    sketch_uniform_random(A, m, n, lda, B.data(), sketch_size, sketch_size, seed);
    auto id_end = std::chrono::high_resolution_clock::now();
    auto id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
    std::cout << "  sketch time: " << id_duration.count() << " ms" << std::endl;
    
    id_start = std::chrono::high_resolution_clock::now();
    // Compute ID on sketched matrix Y
    IDResult<DataType> result = compute_id(
        B.data(),
        sketch_size,
        n,
        sketch_size,
        rank_or_tol,
        0  // max_rank = 0 means auto-detect from tolerance
    );
    id_end = std::chrono::high_resolution_clock::now();
    id_duration = std::chrono::duration_cast<std::chrono::milliseconds>(id_end - id_start);
    std::cout << "  id time: " << id_duration.count() << " ms" << std::endl;
    return result;
}

/**
 * @brief Seed generator for deterministic but different seeds
 */
class SeedGenerator {
private:
    std::mt19937_64 rng_;
    
public:
    explicit SeedGenerator(uint64_t master_seed = 12345) : rng_(master_seed) {}
    
    uint64_t next() {
        return rng_();
    }
    
    void reset(uint64_t master_seed) {
        rng_.seed(master_seed);
    }
};

/**
 * @brief FMM compression helper with multiple sketch backends
 */
template<typename DataType>
class FMMCompressor {
private:
    SeedGenerator seed_gen_;
    double tolerance_;
    double sketch_factor_;       // Sketch size = sketch_factor * n
    int64_t sketch_threshold_;   // Use sketching if m > threshold
    int sparse_k_;               // Sparsity for sparse sketching
    
public:
    enum SketchType {
        DIRECT,        // No sketching, direct ID
        UNIFORM,       // Dense uniform random projection
        SPARSE         // Sparse random projection (faster for very tall matrices)
    };
    
    FMMCompressor(
        double tol = 1e-6, 
        double sketch_factor = 2.0,
        int64_t sketch_threshold = 1000,
        int sparse_k = 3,
        uint64_t master_seed = 12345)
        : seed_gen_(master_seed), 
          tolerance_(tol), 
          sketch_factor_(sketch_factor),
          sketch_threshold_(sketch_threshold),
          sparse_k_(sparse_k) {}
    
    /**
     * @brief Compress far-field with automatic backend selection
     */
    IDResult<DataType> compress_far_field(
        const DataType* A, std::vector<DataType>& B, 
        int64_t m, int64_t n, int64_t lda,
        SketchType sketch_type = SPARSE) {
        
        uint64_t seed = seed_gen_.next();
        
        // For smaller matrices, use direct ID
        if (m <= sketch_threshold_ || m <= 10 * n) {
            B.clear();
            return compute_id(const_cast<DataType*>(A), m, n, lda, tolerance_, 0);
        }
        
        // Choose sketch backend
        switch (sketch_type) {
            case UNIFORM:
                return compute_id_uniform_sketch(A, B, m, n, lda, 
                                                tolerance_, sketch_factor_, seed);
            case SPARSE:
                return compute_id_sparse_sketch(A, B, m, n, lda, 
                                               tolerance_, sketch_factor_, sparse_k_, seed);
            default:
                B.clear();
                return compute_id(const_cast<DataType*>(A), m, n, lda, tolerance_, 0);
        }
    }
    
    /**
     * @brief Compress with explicit control over sketch type and seed
     */
    IDResult<DataType> compress_far_field_explicit(
        const DataType* A, std::vector<DataType>& B, 
        int64_t m, int64_t n, int64_t lda,
        SketchType sketch_type,
        uint64_t seed) {
        
        switch (sketch_type) {
            case UNIFORM:
                return compute_id_uniform_sketch(A, B, m, n, lda, 
                                                tolerance_, sketch_factor_, seed);
            case SPARSE:
                return compute_id_sparse_sketch(A, B, m, n, lda, 
                                               tolerance_, sketch_factor_, sparse_k_, seed);
            default:
                B.clear();
                return compute_id(const_cast<DataType*>(A), m, n, lda, tolerance_, 0);
        }
    }
    
    void reset_seeds(uint64_t master_seed = 12345) { seed_gen_.reset(master_seed); }
    void set_sketch_factor(double factor) { sketch_factor_ = factor; }
    void set_sparse_k(int k) { sparse_k_ = k; }
    
    double get_sketch_factor() const { return sketch_factor_; }
    int get_sparse_k() const { return sparse_k_; }
};

} // namespace fmm

#endif // ID_DECOMPOSITION_HPP
