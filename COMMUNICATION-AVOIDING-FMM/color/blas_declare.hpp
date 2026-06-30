#ifndef BLAS_DECLARE_HPP
#define BLAS_DECLARE_HPP

// extern "C" {
//     // =====================
//     // Real double
//     // =====================

//     // Cholesky factorization: A = L*L^T
//     void dpotrf_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO);

//     // LU factorization: A = P*L*U
//     void dgetrf_(const int* M, const int* N, double* A, const int* LDA,
//                  int* IPIV, int* INFO);

//     // Matrix inverse using LU
//     void dgetri_(const int* N, double* A, const int* LDA, const int* IPIV,
//                  double* WORK, const int* LWORK, int* INFO);

//     // =====================
//     // Complex double
//     // =====================

//     // Cholesky factorization: A = L*L^H
//     void zpotrf_(const char* UPLO, const int* N, std::complex<double>* A, const int* LDA, int* INFO);

//     // LU factorization: A = P*L*U
//     void zgetrf_(const int* M, const int* N, std::complex<double>* A, const int* LDA,
//                  int* IPIV, int* INFO);

//     // Matrix inverse using LU
//     void zgetri_(const int* N, std::complex<double>* A, const int* LDA, const int* IPIV,
//                  std::complex<double>* WORK, const int* LWORK, int* INFO);
                 
//     // Double precision rank-revealing QR with column pivoting
//     void dgeqp3rk_(
//         int* M, int* N, int* NRHS,
//         int* KMAX, double* ABSTOL, double* RELTOL,
//         double* A, int* LDA,
//         int* K, double* MAXC2NRMK, double* RELMAXC2NRMK,
//         int* JPIV, double* TAU,
//         double* WORK, int* LWORK,
//         int* IWORK, int* INFO
//     );
    
//     // Solve triangular system
//     void dtrtrs_(
//         char* UPLO, char* TRANS, char* DIAG,
//         int* N, int* NRHS,
//         double* A, int* LDA,
//         double* B, int* LDB,
//         int* INFO
//     );


//     // Additional LAPACK/BLAS declarations
//     // QR factorization without pivoting
//     void dgeqrf_(
//         int* M, int* N,
//         double* A, int* LDA,
//         double* TAU,
//         double* WORK, int* LWORK,
//         int* INFO
//     );
    
//     // QR factorization with column pivoting
//     void dgeqp3_(
//         int* M, int* N,
//         double* A, int* LDA,
//         int* JPIV, double* TAU,
//         double* WORK, int* LWORK,
//         int* INFO
//     );
    
//     // Generate orthogonal matrix Q from QR
//     void dorgqr_(
//         int* M, int* N, int* K,
//         double* A, int* LDA,
//         double* TAU,
//         double* WORK, int* LWORK,
//         int* INFO
//     );
    
//     // Matrix-matrix multiply: C = alpha*A*B + beta*C
//    void dgemm_(const char* TRANSA, const char* TRANSB,
//                 const int* M, const int* N, const int* K,
//                 const double* ALPHA, const double* A, const int* LDA,
//                 const double* B, const int* LDB,
//                 const double* BETA, double* C, const int* LDC);
    
//     // Triangular solve: B = alpha*inv(A)*B
//     void dtrsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
//                 const int* M, const int* N, const double* ALPHA,
//                 const double* A, const int* LDA, double* B, const int* LDB);
    
//     // Triangular matrix-matrix solve
//     void dtrtrs_(
//         char* UPLO, char* TRANS, char* DIAG,
//         int* N, int* NRHS,
//         double* A, int* LDA,
//         double* B, int* LDB,
//         int* INFO
//     );
    
//     // Copy matrix
//     void dlacpy_(
//         char* UPLO,
//         int* M, int* N,
//         double* A, int* LDA,
//         double* B, int* LDB
//     );

//     // BLAS AXPY: y = alpha*x + y
//     void daxpy_(
//         int* N,
//         double* ALPHA,
//         double* X, int* INCX,
//         double* Y, int* INCY
//     );
//     // BLAS Level 2
//     void dgemv_(const char* trans, const int* m, const int* n,
//                 const double* alpha, const double* a, const int* lda,
//                 const double* x, const int* incx,
//                 const double* beta, double* y, const int* incy);
    
//     void zgemv_(const char* trans, const int* m, const int* n,
//                 const void* alpha, const void* a, const int* lda,
//                 const void* x, const int* incx,
//                 const void* beta, void* y, const int* incy);
    

//     void sgemv_(const char* trans, const int* m, const int* n,
//                 const float* alpha, const float* a, const int* lda,
//                 const float* x, const int* incx,
//                 const float* beta, float* y, const int* incy);

//     void cgemv_(const char* trans, const int* m, const int* n,
//                 const std::complex<float>* alpha, const std::complex<float>* a, const int* lda,
//                 const std::complex<float>* x, const int* incx,
//                 const std::complex<float>* beta, std::complex<float>* y, const int* incy);

//     // =====================
//     // NRM2 (Euclidean Norm)
//     // =====================
//     double dnrm2_(const int* n, const double* x, const int* incx);
//     float  snrm2_(const int* n, const float* x, const int* incx);

//     // Note: complex nrm2 returns real type (magnitude)
//     double dznrm2_(const int* n, const std::complex<double>* x, const int* incx);
//     float  scnrm2_(const int* n, const std::complex<float>* x, const int* incx);
    
//     // LAPACK solvers
//     void dpotrs_(const char* uplo, const int* n, const int* nrhs,
//                  const double* a, const int* lda,
//                  double* b, const int* ldb, int* info);
    
//     void zpotrs_(const char* uplo, const int* n, const int* nrhs,
//                  const void* a, const int* lda,
//                  void* b, const int* ldb, int* info);
    
//     void dpotrf_(const char* uplo, const int* n,
//                  double* a, const int* lda, int* info);
    
        
// }


extern "C" {
    // =====================
    // OpenBLAS threading
    // =====================
    void openblas_set_num_threads(int num_threads);
    int openblas_get_num_threads(void);

    // =====================
    // Cholesky factorization
    // =====================
    void dpotrf_(const char* UPLO, const int* N, double* A, const int* LDA, int* INFO);
    void zpotrf_(const char* UPLO, const int* N, std::complex<double>* A, const int* LDA, int* INFO);

    // =====================
    // Cholesky solve
    // =====================
    void dpotrs_(const char* UPLO, const int* N, const int* NRHS,
                 const double* A, const int* LDA,
                 double* B, const int* LDB, int* INFO);
    void zpotrs_(const char* UPLO, const int* N, const int* NRHS,
                 const std::complex<double>* A, const int* LDA,
                 std::complex<double>* B, const int* LDB, int* INFO);

    // =====================
    // LU factorization: A = P*L*U
    // =====================
    void dgetrf_(const int* M, const int* N, double* A, const int* LDA,
                 int* IPIV, int* INFO);
    void zgetrf_(const int* M, const int* N, std::complex<double>* A, const int* LDA,
                 int* IPIV, int* INFO);

    // =====================
    // Matrix inverse using LU
    // =====================
    void dgetri_(const int* N, double* A, const int* LDA, const int* IPIV,
                 double* WORK, const int* LWORK, int* INFO);
    void zgetri_(const int* N, std::complex<double>* A, const int* LDA, const int* IPIV,
                 std::complex<double>* WORK, const int* LWORK, int* INFO);

    // =====================
    // Solve using LU factors
    // =====================
    void dgetrs_(const char* TRANS, const int* N, const int* NRHS,
                 double* A, const int* LDA, int* IPIV,
                 double* B, const int* LDB, int* INFO);
    void zgetrs_(const char* TRANS, const int* N, const int* NRHS,
                 std::complex<double>* A, const int* LDA, int* IPIV,
                 std::complex<double>* B, const int* LDB, int* INFO);

    // =====================
    // Bunch-Kaufman factorization for complex symmetric matrices (A = P*L*D*L^T*P^T)
    // =====================
    void zsytrf_(const char* UPLO, const int* N, std::complex<double>* A, const int* LDA,
                 int* IPIV, std::complex<double>* WORK, const int* LWORK, int* INFO);
    void zsytrs_(const char* UPLO, const int* N, const int* NRHS,
                 const std::complex<double>* A, const int* LDA, const int* IPIV,
                 std::complex<double>* B, const int* LDB, int* INFO);

    // =====================
    // QR factorization without pivoting
    // =====================
    void dgeqrf_(int* M, int* N,
                 double* A, int* LDA,
                 double* TAU,
                 double* WORK, int* LWORK, int* INFO);
    void zgeqrf_(int* M, int* N,
                 std::complex<double>* A, int* LDA,
                 std::complex<double>* TAU,
                 std::complex<double>* WORK, int* LWORK, int* INFO);

    // =====================
    // QR factorization with column pivoting
    // =====================
    void dgeqp3_(int* M, int* N,
                 double* A, int* LDA,
                 int* JPIV, double* TAU,
                 double* WORK, int* LWORK, int* INFO);
    void zgeqp3_(int* M, int* N,
                 std::complex<double>* A, int* LDA,
                 int* JPIV, std::complex<double>* TAU,
                 std::complex<double>* WORK, int* LWORK,
                 double* RWORK, int* INFO);

    // =====================
    // Rank-revealing QR (GEQP3RK)
    // =====================
    void dgeqp3rk_(int* M, int* N, int* NRHS,
                   int* KMAX, double* ABSTOL, double* RELTOL,
                   double* A, int* LDA,
                   int* K, double* MAXC2NRMK, double* RELMAXC2NRMK,
                   int* JPIV, double* TAU,
                   double* WORK, int* LWORK,
                   int* IWORK, int* INFO);
    void zgeqp3rk_(int* M, int* N, int* NRHS,
                   int* KMAX, double* ABSTOL, double* RELTOL,
                   std::complex<double>* A, int* LDA,
                   int* K, double* MAXC2NRMK, double* RELMAXC2NRMK,
                   int* JPIV, std::complex<double>* TAU,
                   std::complex<double>* WORK, int* LWORK,
                   double* RWORK, int* IWORK, int* INFO);


    // =====================
    // GEMM (Matrix-Matrix multiply): C = alpha*A*B + beta*C
    // =====================
    void sgemm_(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K,
            const float* ALPHA, const float* A, const int* LDA,
            const float* B, const int* LDB,
            const float* BETA, float* C, const int* LDC);
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA, const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA, double* C, const int* LDC);
    void cgemm_(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K,
            const std::complex<float>* ALPHA, const std::complex<float>* A, const int* LDA,
            const std::complex<float>* B, const int* LDB,
            const std::complex<float>* BETA, std::complex<float>* C, const int* LDC);
    void zgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const std::complex<double>* ALPHA, const std::complex<double>* A, const int* LDA,
                const std::complex<double>* B, const int* LDB,
                const std::complex<double>* BETA, std::complex<double>* C, const int* LDC);

    // =====================
    // TRSM (Triangular solve, matrix): B = alpha*inv(A)*B
    // =====================
    void dtrsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
                const int* M, const int* N, const double* ALPHA,
                const double* A, const int* LDA, double* B, const int* LDB);
    void ztrsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
                const int* M, const int* N, const std::complex<double>* ALPHA,
                const std::complex<double>* A, const int* LDA,
                std::complex<double>* B, const int* LDB);
    void strsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
            const int* M, const int* N, const float* ALPHA,
            const float* A, const int* LDA, float* B, const int* LDB);
    void ctrsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
            const int* M, const int* N, const std::complex<float>* ALPHA,
            const std::complex<float>* A, const int* LDA,
            std::complex<float>* B, const int* LDB);

    // =====================
    // TRTRS (Triangular solve, LAPACK): A*X = B
    // =====================
    void dtrtrs_(char* UPLO, char* TRANS, char* DIAG,
                 int* N, int* NRHS,
                 double* A, int* LDA,
                 double* B, int* LDB, int* INFO);
    void ztrtrs_(char* UPLO, char* TRANS, char* DIAG,
                 int* N, int* NRHS,
                 std::complex<double>* A, int* LDA,
                 std::complex<double>* B, int* LDB, int* INFO);

    // =====================
    // LACPY (Copy matrix)
    // =====================
    void dlacpy_(char* UPLO, int* M, int* N,
                 double* A, int* LDA,
                 double* B, int* LDB);
    void zlacpy_(char* UPLO, int* M, int* N,
                 std::complex<double>* A, int* LDA,
                 std::complex<double>* B, int* LDB);

    // =====================
    // AXPY: y = alpha*x + y
    // =====================
    void daxpy_(int* N, double* ALPHA,
                double* X, int* INCX,
                double* Y, int* INCY);
    void zaxpy_(int* N, std::complex<double>* ALPHA,
                std::complex<double>* X, int* INCX,
                std::complex<double>* Y, int* INCY);

    // =====================
    // GEMV (Matrix-Vector): y = alpha*A*x + beta*y
    // =====================
    void dgemv_(const char* TRANS, const int* M, const int* N,
                const double* ALPHA, const double* A, const int* LDA,
                const double* X, const int* INCX,
                const double* BETA, double* Y, const int* INCY);
    void sgemv_(const char* TRANS, const int* M, const int* N,
                const float* ALPHA, const float* A, const int* LDA,
                const float* X, const int* INCX,
                const float* BETA, float* Y, const int* INCY);
    void zgemv_(const char* TRANS, const int* M, const int* N,
                const std::complex<double>* ALPHA, const std::complex<double>* A, const int* LDA,
                const std::complex<double>* X, const int* INCX,
                const std::complex<double>* BETA, std::complex<double>* Y, const int* INCY);
    void cgemv_(const char* TRANS, const int* M, const int* N,
                const std::complex<float>* ALPHA, const std::complex<float>* A, const int* LDA,
                const std::complex<float>* X, const int* INCX,
                const std::complex<float>* BETA, std::complex<float>* Y, const int* INCY);

    // =====================
    // NRM2 (Euclidean Norm)
    // =====================
    double dnrm2_(const int* N, const double* X, const int* INCX);
    float  snrm2_(const int* N, const float* X, const int* INCX);
    double dznrm2_(const int* N, const std::complex<double>* X, const int* INCX);
    float  scnrm2_(const int* N, const std::complex<float>* X, const int* INCX);

    // =====================
    // TRMV (Triangular Matrix-Vector multiply): x = A*x
    // =====================
    void dtrmv_(const char* UPLO, const char* TRANS, const char* DIAG,
                const int* N, const double* A, const int* LDA,
                double* X, const int* INCX);
    void ztrmv_(const char* UPLO, const char* TRANS, const char* DIAG,
                const int* N, const std::complex<double>* A, const int* LDA,
                std::complex<double>* X, const int* INCX);

    // =====================
    // LASWP (Row permutations): apply IPIV row swaps to a matrix
    // =====================
    void dlaswp_(const int* N, double* A, const int* LDA,
                 const int* K1, const int* K2, const int* IPIV, const int* INCX);
    void zlaswp_(const int* N, std::complex<double>* A, const int* LDA,
                 const int* K1, const int* K2, const int* IPIV, const int* INCX);
}


template <typename DataType>
void trsm_(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
           const int* M, const int* N, const DataType* ALPHA,
           const DataType* A, const int* LDA,
           DataType* B, const int* LDB) {
    if constexpr (std::is_same_v<DataType, double>) {
        dtrsm_(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
    } else if constexpr (std::is_same_v<DataType, float>) {
        strsm_(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        ztrsm_(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        ctrsm_(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
    } else {
        static_assert(sizeof(DataType) == 0, "Unsupported DataType for trsm_");
    }
}

template <typename DataType>
void gemm_(const char* TRANSA, const char* TRANSB,
           const int* M, const int* N, const int* K,
           const DataType* ALPHA, const DataType* A, const int* LDA,
           const DataType* B, const int* LDB,
           const DataType* BETA, DataType* C, const int* LDC) {
    if constexpr (std::is_same_v<DataType, double>) {
        dgemm_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
    } else if constexpr (std::is_same_v<DataType, float>) {
        sgemm_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgemm_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        cgemm_(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
    } else {
        static_assert(sizeof(DataType) == 0, "Unsupported DataType for gemm_");
    }
}


// =====================
// Complex Symmetric Cholesky: A = L * L^T (no conjugation)
// For complex symmetric matrices (e.g., Helmholtz BEM)
// =====================

// Unblocked kernel for small diagonal blocks
static void zsytrf_chol_unblocked_(
    const int N,
    std::complex<double>* A,
    const int LDA,
    int* INFO
) {
    *INFO = 0;
    for (int k = 0; k < N; ++k) {
        std::complex<double>& Akk = A[k + k * LDA];

        // Akk -= sum_{i=0}^{k-1} L(k,i) * L(k,i)  (no conjugation)
        for (int i = 0; i < k; ++i) {
            Akk -= A[k + i * LDA] * A[k + i * LDA];
        }

        // Check for zero pivot
        if (std::abs(Akk) == 0.0) {
            *INFO = k + 1; // 1-based index
            return;
        }

        Akk = std::sqrt(Akk);
        std::complex<double> AkkInv = 1.0 / Akk;

        // Update column k below diagonal
        // A(j,k) -= sum_{i=0}^{k-1} A(j,i) * A(k,i)  for j = k+1..N-1
        for (int i = 0; i < k; ++i) {
            std::complex<double> Lki = A[k + i * LDA];
            for (int j = k + 1; j < N; ++j) {
                A[j + k * LDA] -= A[j + i * LDA] * Lki;
            }
        }

        // Scale column k
        for (int j = k + 1; j < N; ++j) {
            A[j + k * LDA] *= AkkInv;
        }
    }
}

// Blocked complex symmetric Cholesky: A = L * L^T
// Interface mirrors zpotrf_
//   UPLO: 'L' only (lower triangular stored)
//   N:    order of matrix
//   A:    on entry, complex symmetric matrix; on exit, L in lower triangle
//   LDA:  leading dimension of A
//   INFO: 0 on success, k if pivot k is zero (1-based)
//   NB:   block size (optional tuning parameter, default 64)
void zsychol_(
    const char* UPLO,
    const int* N_,
    std::complex<double>* A,
    const int* LDA_,
    int* INFO,
    int NB = 64
) {
    const int N = *N_;
    const int LDA = *LDA_;
    *INFO = 0;

    if (*UPLO != 'L' && *UPLO != 'l') {
        *INFO = -1;
        return;
    }

    if (N <= 0) return;

    // Fall back to unblocked for small matrices
    if (N <= NB) {
        zsytrf_chol_unblocked_(N, A, LDA, INFO);
        return;
    }

    // Blocked algorithm
    std::complex<double> one(1.0, 0.0);
    std::complex<double> neg_one(-1.0, 0.0);

    for (int k = 0; k < N; k += NB) {
        int kb = std::min(NB, N - k);
        int tail = N - k - kb;

        // Pointer to A(k, k)
        std::complex<double>* Akk = &A[k + k * LDA];

        // Step 1: Factor diagonal block A(k:k+kb-1, k:k+kb-1)
        zsytrf_chol_unblocked_(kb, Akk, LDA, INFO);
        if (*INFO != 0) {
            *INFO += k; // Adjust to global index
            return;
        }

        if (tail > 0) {
            // Pointer to A(k+kb, k)
            std::complex<double>* Aik = &A[(k + kb) + k * LDA];

            // Step 2: Solve panel block
            // A(k+kb:N-1, k:k+kb-1) = A(k+kb:N-1, k:k+kb-1) * L(k:k+kb-1, k:k+kb-1)^{-T}
            // Use TRSM with 'T' (transpose, no conjugation)
            // B := B * inv(L^T)  =>  SIDE='R', UPLO='L', TRANSA='T', DIAG='N'
            ztrsm_("R", "L", "T", "N",
                   &tail, &kb, &one,
                   Akk, &LDA,
                   Aik, &LDA);

            // Step 3: Update trailing submatrix
            // A(k+kb:N-1, k+kb:N-1) -= A(k+kb:N-1, k:k+kb-1) * A(k+kb:N-1, k:k+kb-1)^T
            // Use ZGEMM with 'N', 'T' (no conjugation on transpose)
            std::complex<double>* Aii = &A[(k + kb) + (k + kb) * LDA];
            zgemm_("N", "T",
                   &tail, &tail, &kb,
                   &neg_one, Aik, &LDA,
                   Aik, &LDA,
                   &one, Aii, &LDA);
        }
    }
}


// Complex symmetric Cholesky solve: solves A*X = B given L from zsychol_
// Interface mirrors zpotrs_
//   UPLO:  'L' only
//   N:     order of matrix
//   NRHS:  number of right-hand sides
//   A:     factorized matrix from zsychol_
//   LDA:   leading dimension of A
//   B:     on entry, RHS; on exit, solution
//   LDB:   leading dimension of B
//   INFO:  0 on success
void zsychol_solve_(
    const char* UPLO,
    const int* N,
    const int* NRHS,
    const std::complex<double>* A,
    const int* LDA,
    std::complex<double>* B,
    const int* LDB,
    int* INFO
) {
    *INFO = 0;

    if (*UPLO != 'L' && *UPLO != 'l') {
        *INFO = -1;
        return;
    }

    std::complex<double> one(1.0, 0.0);

    // Forward solve: L * Z = B
    ztrsm_("L", "L", "N", "N", N, NRHS, &one,
           A, LDA, B, LDB);

    // Back solve: L^T * X = Z  (transpose, no conjugation)
    ztrsm_("L", "L", "T", "N", N, NRHS, &one,
           A, LDA, B, LDB);
}


#endif 
