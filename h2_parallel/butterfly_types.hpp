#pragma once

// H2 (format-7) integration: foundational traits, kernel, solver struct, options.

// C++ Standard Library
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <cmath>
#include <iostream>
#include <type_traits>



// MPI
#include <mpi.h>

// FMM Verification
#include <iomanip>
#include <limits>
#include <chrono>
#include <cstdio>

// FMM
#include "factorization.hpp"
#include "runtime_thread_support.hpp"
#include "solver.hpp"
#include "apply_mul.hpp"
#include "tree_impl.hpp"
//#include "kernel.hpp"
#include "id_decomposition.hpp"


// Extra

#include <array>
#include <atomic>
#include <mutex>
#include <omp.h>
#include <random>
#include <complex>
#include <sched.h>
#include <algorithm>
#include <cctype>
#include <thread>
#include <vector>
#include <unordered_set>
#include <unordered_map>


namespace butterfly {
using namespace fmm;

// C_DT (C built-in complex) → the std:: library type the FMM is written for
template<class T> struct fmm_data           { using type = T; };
template<>        struct fmm_data<_Complex double> { using type = std::complex<double>; };
template<>        struct fmm_data<_Complex float>  { using type = std::complex<float>;  };

// inverse: DataType → the C boundary type C_FuncZmn actually uses
template<class T> struct kernel_value_type                        { using type = T; };
template<>        struct kernel_value_type<std::complex<double>>  { using type = _Complex double; };
template<>        struct kernel_value_type<std::complex<float>>   { using type = _Complex float;  };

template <typename T> struct is_complex : std::false_type {};                 // (1) primary
template <typename T> struct is_complex<std::complex<T>> : std::true_type {}; // (2) specialization
template <typename T> inline constexpr bool is_complex_v = is_complex<T>::value; // (3) alias

// // in butterfly_h2/h2_quant.hpp — YOU define this, the user fills it
// struct H2QuantApp {
//     fmm::KernelKind kernel_kind = fmm::KernelKind::LAPLACE;  // default member initializers
//     double wave_divisor    = 32.0;   // Helmholtz only
//     double length_scale    = 0.1;    // Matérn only
//     double nugget          = 1e-6;   // Matérn only
//     double kappa           = 10.0;   // Yukawa only
//     int    num_proxy       = -1;     // -1 = auto (32 for 2D, 256 for 3D)
//     int64_t reduction_threshold = 0;
//     const double* locations = nullptr;
// };

enum class KernelKind {
    LAPLACE,
    HELMHOLTZ,
    MATERN52,
    YUKAWA
};

enum class NumberKind {
    REAL,
    COMPLEX
};

inline const char* kernel_kind_to_string(KernelKind kind) {
    switch (kind) {
        case KernelKind::LAPLACE:
            return "laplace";
        case KernelKind::HELMHOLTZ:
            return "helmholtz";
        case KernelKind::MATERN52:
            return "matern52";
        case KernelKind::YUKAWA:
            return "yukawa";
        default:
            return "unknown";
    }
}

inline const char* number_kind_to_string(NumberKind kind) {
    switch (kind) {
        case NumberKind::REAL:
            return "real";
        case NumberKind::COMPLEX:
            return "complex";
        default:
            return "unknown";
    }
}

struct ProgramOptions {
    int num_levels = 0;
    int64_t N = 0;
    int64_t grid_size = 0;
    double tolerance = 0.0;
    //To DO:NumberKind and KernelKind are not used, probably can delete in future
    KernelKind kernel_kind = KernelKind::LAPLACE;
    NumberKind number_kind = NumberKind::REAL;
    int dimension = 3;
    int64_t reduction_threshold = 0;
    int num_proxy = -1;
    double wave_divisor = 32.0;
    double length_scale = 0.1;   // Matérn length scale ℓ
    double nugget = 1e-6;        // Matérn diagonal nugget σ_n²
    double kappa = 10.0;         // Yukawa screening parameter κ
    int cond_samples = 0;        // Power iteration samples for condition number estimate (0 = skip)
    int precon = 1;              // 1/3: RS-S factorization, 2: compression-only H2
    int verbosity = 0;           // -1: quiet, 0: summaries, 1+: detailed progress
};

enum class H2BuildState {
    UNBUILT,
    H2_COMPRESSED,
    RS_FACTORIZED
};

template<typename CoordType, typename DataType>
struct H2Kernel {
    using kerData = typename kernel_value_type<DataType>::type;
    void (*kernel)(int*, int*, kerData*, void*) = nullptr;
    void* quant = nullptr;
    mutable std::vector<double> entryeval_time_per_thread;


    H2Kernel() = default;

    H2Kernel(void (*kernel_)(int*, int*, kerData*, void*), void* quant_)
        : kernel(kernel_), quant(quant_) {}

    // Evaluate the single (i, j) entry into a column-major matrix A with leading dimension lda:
    //   A[i + j*lda] = K(i, j),  i = row, j = column (0-based global DOF indices).
    // Routes to the same C_FuncZmn callback as evaluate_block_by_index (1-based), so the entry
    // matches exactly what the H2 solver compressed/factored -- including the self/diagonal term
    // when i == j. This is what lets direct-verify test the true operator, not a built-in kernel.
    void evaluate_by_index(const int64_t i, const int64_t j, DataType* A, int64_t lda) const {
        int m = static_cast<int>(i) + 1;   // 1-based row
        int n = static_cast<int>(j) + 1;   // 1-based column
        kernel(&m, &n, reinterpret_cast<kerData*>(&A[i + j * lda]), quant);
    }

    // Fill the x_size-by-y_size block A (column-major, leading dimension lda):
    //   A[i + j*lda] = K( x_indices[i], y_indices[j] )
    // x_indices index the rows, y_indices the columns. Indices are 0-based
    // global DOF indices (as in BoxData::point_indices); C_FuncZmn expects
    // 1-based, so we add 1.
    void evaluate_block_by_index(const int64_t* x_indices, int64_t x_size,
                                 const int64_t* y_indices, int64_t y_size,
                                 DataType* A, int64_t lda) const {
        
        
        double t0 = MPI_Wtime();
        for (int64_t j = 0; j < y_size; ++j) {
            int n = static_cast<int>(y_indices[j]) + 1;   // 1-based column
            for (int64_t i = 0; i < x_size; ++i) {
                int m = static_cast<int>(x_indices[i]) + 1;   // 1-based row
                kernel(&m, &n, reinterpret_cast<kerData*>(&A[i + j * lda]), quant);
            }
        }
        double tf = MPI_Wtime() - t0;
        int tid = omp_get_thread_num();
        if (tid < static_cast<int>(entryeval_time_per_thread.size())) {
            entryeval_time_per_thread[tid] += tf;
        }
    }
};

// construct b_mat
template<typename CoordType, typename DataType>
struct H2 {
    //using CoordType = double;
    MPI_Comm comm;
    int64_t N;
    int dimension;
    //temporary comment
    //fmm::HierarchicalFactorization<CoordType, C_DT, KernelType> factorizer;
    
    H2Kernel<CoordType, DataType> kernel;
    std::unique_ptr<fmm::ParallelTree<CoordType, DataType>> tree;
    ProgramOptions options;

    //temporary comment
    //RedistributionPlan redistribution;
    int64_t last_factor_rankmax = 0;
    size_t factorization_memory = 0;
    H2BuildState build_state = H2BuildState::UNBUILT;
    bool factorized = false;
};

// The sparse test vector x = sum_k weight[k] * e_{idx[k]}  (zero elsewhere).
template<typename DataType>
struct SparseTestVector {
    std::vector<int64_t>  idx;      // support: global column indices of the nonzeros
    std::vector<DataType> weight;   // value of x at each support index
};


} // namespace butterfly
