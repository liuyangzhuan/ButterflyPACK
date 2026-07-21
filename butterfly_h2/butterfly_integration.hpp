#pragma once

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
#include "kernel.hpp"
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
    bool factorized = false;
};

inline int default_num_levels(int64_t grid_size) {
  int k = 0;
  while (grid_size > 1) {  
      grid_size /= 2;
      ++k;
  }
  return k;
}

inline int64_t default_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 256 : 4096;
}

inline int64_t min_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 64 : 512;
}

inline int default_num_proxy_for_dimension(int dimension) {
    return (dimension == 2) ? 32 : 256;
}

inline int64_t is_power_of_base(int64_t value, int64_t base) {
    if (value < 1 || base < 2) {
        return -1;
    }
    double power = 0;
    while (value % base == 0) {
        value /= base;
        power += 1;
    }
    if (value == 1) {
      return power;
    }
    else {
      return -1;
    }
}

inline bool is_valid_reduction_threshold(int64_t reduction_threshold, int dimension) {
    const int64_t base = (dimension == 2) ? 4 : 8;
    return is_power_of_base(reduction_threshold, base) >= 0;
}

inline std::string reduction_threshold_pattern(int dimension) {
    if (dimension == 2) {
        return "4^k (1, 4, 16, 64, 256, ...)";
    }
    return "8^k (1, 8, 64, 512, 4096, ...)";
}

/**
 * @brief Derive global domain bounds from the point coordinates.
 *
 * Locations is point-major (interleaved), stride = dimension:
 *   [x0,y0,z0, x1,y1,z1, ...]     (matches tree_impl.hpp:801-803)
 *
 * Fills bounds as [xmin,xmax, ymin,ymax, zmin,zmax]; z entries are 0 in 2D
 * (compute_box_geometry only reads them when dimension == 3).
 *
 * A small relative pad is applied so points on the upper face fall strictly
 * inside the last box: point_to_morton computes (p - min)/box_size, which at
 * p == max would land exactly on boxes_per_dim and be silently clamped.
 *
 * No MPI reduction: every rank holds the full global point array.
 */
template<typename CoordType>
inline void compute_global_bounds(const CoordType* Locations,
                                  int64_t num_points,
                                  int dimension,
                                  CoordType bounds[6]) {
    if (Locations == nullptr || num_points <= 0) {
        throw std::invalid_argument("compute_global_bounds: no points provided");
    }
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("compute_global_bounds: dimension must be 2 or 3");
    }

    CoordType lo[3] = { std::numeric_limits<CoordType>::max(),
                        std::numeric_limits<CoordType>::max(),
                        std::numeric_limits<CoordType>::max() };
    CoordType hi[3] = { std::numeric_limits<CoordType>::lowest(),
                        std::numeric_limits<CoordType>::lowest(),
                        std::numeric_limits<CoordType>::lowest() };

    for (int64_t i = 0; i < num_points; ++i) {
        for (int d = 0; d < dimension; ++d) {
            const CoordType v = Locations[i * dimension + d];
            if (!std::isfinite(v)) {
                throw std::runtime_error("compute_global_bounds: non-finite coordinate");
            }
            lo[d] = std::min(lo[d], v);
            hi[d] = std::max(hi[d], v);
        }
    }

    // Pad; the (span == 0) case covers a degenerate/planar dimension.
    for (int d = 0; d < dimension; ++d) {
        const CoordType span = hi[d] - lo[d];
        const CoordType pad  = (span > CoordType(0) ? span : CoordType(1)) * CoordType(1e-6);
        lo[d] -= pad;
        hi[d] += pad;
    }

    bounds[0] = lo[0]; bounds[1] = hi[0];
    bounds[2] = lo[1]; bounds[3] = hi[1];
    if (dimension == 3) { bounds[4] = lo[2]; bounds[5] = hi[2]; }
    else                { bounds[4] = CoordType(0);   bounds[5] = CoordType(0);   }
}




// ProgramOptions
  // int num_levels = nlevel;
  // int64_t N = Npo;
  // int64_t grid_size = 0;
  // double tolerance = 0.0;
  // fmm::KernelKind kernel_kind = fmm::KernelKind::LAPLACE;
  // NumberKind number_kind = NumberKind::REAL;
  // int dimension = Ndim;
  // int64_t reduction_threshold = 0;
  // int num_proxy = -1;
  // double wave_divisor = 32.0;
  // double length_scale = 0.1;   // Matérn length scale ℓ
  // double nugget = 1e-6;        // Matérn diagonal nugget σ_n²
  // double kappa = 10.0;         // Yukawa screening parameter κ
  // int cond_samples = 0;        // Power iteration samples for condition number estimate (0 = skip)
inline ProgramOptions parse_program_options(int* Npo, int* Ndim, double* Locations, 
  double tolerance) {

    ProgramOptions h2_options;
    
    h2_options.N = *Npo;
    if (h2_options.N <= 0) {
        throw std::invalid_argument("N must be positive.");
    }

    h2_options.dimension = *Ndim;
    if (h2_options.dimension != 2 && h2_options.dimension != 3) {
        throw std::invalid_argument("H2 solver dimension must be either 2 or 3.");
    }

    //Compute grid_size from Npo and dimension
    // grid_size such that grid_size^dimension == N
    int64_t grid_size = std::llround(std::pow((double)*Npo, 1.0 / *Ndim));

    int64_t check = 1;
    for (int d = 0; d < *Ndim; ++d) check *= grid_size;

    if (check == *Npo) {
        h2_options.grid_size = grid_size;
    } else {
        throw std::invalid_argument(
            "N must equal grid_size^dimension for FFT verification.");
    }

    if (h2_options.grid_size < 2) {
      throw std::invalid_argument(
        "grid size must be >=2 for H2 solver."
      );
    }

    // for h2 just calculate in the package
    h2_options.num_levels = default_num_levels(grid_size);

    // need to check DataType: want to support double and double complex
    h2_options.tolerance = tolerance;


    // H2QuantApp* quant = static_cast<H2QuantApp*> (C_QuantApp);
    // h2_options.wave_divisor = quant->wave_divisor;
    // h2_options.length_scale = quant->length_scale;
    // h2_options.nugget = quant->nugget;
    // h2_options.kappa = quant->kappa;

    // To Do: reduction_threshold not defined in option, we may want to add these properties for butterfly users

    // not using proxy points right now
    h2_options.num_proxy = 0;
    
    if (h2_options.num_levels <= 0) {
        throw std::invalid_argument("num_levels must be positive.");
    }
    if (!(h2_options.tolerance > 0.0)) {
        throw std::invalid_argument("tolerance must be positive.");
    }
    if (h2_options.reduction_threshold == 0) {
        h2_options.reduction_threshold =
            default_reduction_threshold_for_dimension(h2_options.dimension);
    }
    if (h2_options.num_proxy == -1) {
        h2_options.num_proxy = default_num_proxy_for_dimension(h2_options.dimension);
    }

    if (h2_options.reduction_threshold <= 0) {
        throw std::invalid_argument("reduction_threshold must be positive.");
    }
    if (h2_options.reduction_threshold < min_reduction_threshold_for_dimension(h2_options.dimension)) {
        throw std::invalid_argument(
            "reduction_threshold must be at least " +
            std::to_string(min_reduction_threshold_for_dimension(h2_options.dimension)) +
            " for dimension " + std::to_string(h2_options.dimension) + ".");
    }
    if (!is_valid_reduction_threshold(h2_options.reduction_threshold, h2_options.dimension)) {
        throw std::invalid_argument(
            "reduction_threshold must be " +
            reduction_threshold_pattern(h2_options.dimension) +
            " for dimension " + std::to_string(h2_options.dimension) + ".");
    }
    if (h2_options.num_proxy < -1) {
        throw std::invalid_argument(
            "num_proxy must be -1 (default), 0, or a positive integer.");
    }
    if (h2_options.dimension == 3 && h2_options.num_proxy == 1) {
        throw std::invalid_argument(
            "num_proxy=1 is invalid in 3D. Use 0 to disable proxies or >= 2.");
    }
    // if (!(h2_options.wave_divisor > 0.0)) {
    //     throw std::invalid_argument("wave_divisor must be positive.");
    // }

    // if (h2_options.kernel_kind == fmm::KernelKind::LAPLACE &&
    //     h2_options.number_kind != NumberKind::REAL) {
    //     throw std::invalid_argument(
    //         "Unsupported combination: Laplace currently supports only real number type.");
    // }

    // if (h2_options.kernel_kind == fmm::KernelKind::HELMHOLTZ &&
    //     h2_options.number_kind != NumberKind::COMPLEX) {
    //     throw std::invalid_argument(
    //         "Unsupported combination: Helmholtz currently supports only complex number type.");
    // }

    // if (h2_options.kernel_kind == fmm::KernelKind::MATERN52) {
    //     if (h2_options.number_kind != NumberKind::REAL) {
    //         throw std::invalid_argument(
    //             "Unsupported combination: Matern52 currently supports only real number type.");
    //     }
    //     if (!(h2_options.length_scale > 0.0)) {
    //         throw std::invalid_argument("length_scale must be positive.");
    //     }
    //     if (!(h2_options.nugget >= 0.0)) {
    //         throw std::invalid_argument("nugget must be non-negative.");
    //     }
    // }

    // if (h2_options.kernel_kind == fmm::KernelKind::YUKAWA) {
    //     if (h2_options.number_kind != NumberKind::REAL) {
    //         throw std::invalid_argument(
    //             "Unsupported combination: Yukawa currently supports only real number type.");
    //     }
    //     if (!(h2_options.kappa > 0.0)) {
    //         throw std::invalid_argument("kappa must be positive.");
    //     }
    // }

    return h2_options;
}



static void allgather_idx_map_to_new2old(
    MPI_Comm comm,
    const std::vector<int>& idx_map,
    int N_global,
    std::vector<int>& new2old,
    int& idxs,
    int& idxe)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int send_count = static_cast<int>(idx_map.size());

    std::vector<int> recv_counts(size, 0);

    MPI_Allgather(
        &send_count, 1, MPI_INT,
        recv_counts.data(), 1, MPI_INT,
        comm
    );

    std::vector<int> displs(size, 0);
    for (int r = 1; r < size; ++r) {
        displs[r] = displs[r - 1] + recv_counts[r - 1];
    }

    int total_count = displs[size - 1] + recv_counts[size - 1];

    if (total_count != N_global) {
        throw std::runtime_error(
            "allgather_idx_map_to_new2old: total gathered idx_map size != N_global"
        );
    }

    // This rank's global-new-index range, zero-based and inclusive
    idxs = displs[rank];

    if (send_count > 0) {
        idxe = idxs + send_count - 1;
    } else {
        idxe = idxs - 1;   // empty local range
    }

    new2old.resize(total_count);

    MPI_Allgatherv(
        idx_map.data(), send_count, MPI_INT,
        new2old.data(), recv_counts.data(), displs.data(), MPI_INT,
        comm
    );
}




template<typename CoordType, typename DataType>
int h2_initiate(H2<CoordType, DataType>* H2_solver, const ProgramOptions& options, CoordType* Locations, int rank, std::vector<int>& new2old, int& idxs, int& idxe) {
  // from H2 struct, set typename CoordType, DataType
  // from H2 struct, get kernel
  // from H2 struct, get rank, size
  if (rank == 0) {
    std::cout << "=== Hierarchical Factorization Test ("
              << options.dimension << "D "
              << kernel_kind_to_string(options.kernel_kind) << " kernel, "
              << number_kind_to_string(options.number_kind) << ") ===" << std::endl;
  }

  CoordType bounds[6];

  compute_global_bounds(
      Locations,
      static_cast<int64_t>(options.N),
      options.dimension,
      bounds
  );

  std::vector<int> idx_map;
  // Provide Locations to nullptr, check compatibility
  auto tree = std::unique_ptr<fmm::ParallelTree<CoordType, DataType>>(
      fmm::create_uniform_tree<CoordType, DataType>(
          Locations,
          options.N,
          options.num_levels,
          bounds,
          options.dimension,
          H2_solver->comm,
          options.reduction_threshold, 
          idx_map));

  H2_solver->tree = std::move(tree); 

  idxs = 0;
  idxe = -1;
  allgather_idx_map_to_new2old(
    H2_solver->comm,
    idx_map,
    static_cast<int>(options.N),
    new2old,
    idxs,
    idxe
  );
  
  idxs++;// this is to convert from 0-based index to 1-based index for Fortran compatibility
  idxe++;


  // To Do: maybe put this section and below to c_bpack_factor
  return 0;
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
    bool verbose = true) {

    // To Do: NEED TO FIX KERNEL!!!!!
    using clock = std::chrono::high_resolution_clock;
    
    int rank = tree->mpi_rank;
    int size = tree->mpi_size;
    DynamicThreadingContext dynamic_threading =
        make_dynamic_threading_context(tree->comm);
    
    int dimension = tree->dimension;
    int num_levels = tree->num_levels;
    int leaf_level = num_levels - 1;
    int num_children = (dimension == 2) ? 4 : 8;
    const int factorization_header_rank =
        smallest_active_rank(tree->levels[leaf_level]);
    
    if (verbose && rank == factorization_header_rank) {
        const char* factorization_name =
            factorization_method == FactorizationMethod::CHOLESKY ? "Cholesky" :
            factorization_method == FactorizationMethod::LU ? "LU" : "None";
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
        const int level_print_rank = smallest_active_rank(level);
        const LevelThreadPlan level_thread_plan =
            configure_level_thread_plan(dynamic_threading, tree->comm, level);
        
        if (!level.is_process_active) {
            if (verbose && rank == level_print_rank) {
                std::cout << "\n===== Level " << current_level << " (Process " << rank << " inactive) =====" << std::endl;
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
        
        if (verbose && rank == level_print_rank) {
            std::cout << "\n===== Level " << current_level << " =====" << std::endl;
            std::cout << "Active processes: " << level.num_active_processes << " from rank: " << rank << std::endl;
            std::cout << "Boxes per process: " << level.num_boxes_local << std::endl;
        }
        print_level_thread_plan(
            dynamic_threading,
            tree->comm,
            current_level,
            "factorization",
            level_print_rank,
            rank,
            level_thread_plan,
            verbose);
        
        // ===== Step 1: Gather ghost and assisting boxes =====
        

        
        // ===== Step 2: Eliminate boxes in colored order =====
        std::chrono::milliseconds elim_duration{0};
        int64_t total_skeleton = 0;
        int64_t total_redundant = 0;
        
        if (current_level > 1 && level.is_process_active) {
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
                if(!is_interior){
                    const auto comm_duration_raw = transport_and_apply_factor_updates_symmetric_onehop(
                    tree, current_level, kernel, pending_updates, false);
                    level_data_exchange += comm_duration_raw;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                    auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);

                    if (verbose && rank == level_print_rank) {
                        std::cout << "  Comm time: "
                                << comm_duration.count()
                                << " ms" << std::endl;
                    }
                }else if(counter == interior_start_loc){
                    const std::vector<int> neighbor_ranks = compute_one_hop_neighbor_ranks(tree, level, current_level);
                    std::vector<int64_t> need_assist;
                    for (const auto& kv : level.assisting_box_points_for_kernel_evaluation) {
                        if(level.eliminated_boxes.count(kv.first) != 0){
                            need_assist.push_back(kv.first);
                        }
                    }
                    if (!need_assist.empty()) {
                        const auto comm_duration_raw = exchange_assisting_for_mortons_onehop(
                            tree, level, current_level, neighbor_ranks, need_assist);
                        level_data_exchange += comm_duration_raw;
                        auto comm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(comm_duration_raw);
                        if (verbose && rank == level_print_rank) {
                            std::cout << "  Comm time: "
                                    << comm_duration.count()
                                    << " ms" << std::endl;
                        }
                    }
                }
                

              

                const auto& color_list = color_bins[static_cast<size_t>(counter)];

                if (verbose && rank == level_print_rank) {
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

                // Final transport after the very last wave
                if (counter == static_cast<int>(color_bins.size()) - 1) {
                    const auto final_comm_duration = transport_and_apply_factor_updates_symmetric_onehop(
                        tree, current_level, kernel, pending_updates, false);
                    level_data_exchange += final_comm_duration;
                    update_neighbor_slicing_for_level(level, is_symmetric);
                }
            }
            
            auto elim_end = std::chrono::high_resolution_clock::now();
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(elim_end - elim_start);

            for (const auto& box : level.local_boxes) {
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
                local_max_skel = std::max<int64_t>(local_max_skel, box.skeleton_indices.size());
            }
        } else {
            // Level 1: Skip elimination (only 4/8 boxes, no far-field)
            if (verbose && rank == level_print_rank) {
                std::cout << "  Skipping elimination at level 1 (final coarsening step)" << std::endl;
            }
        }

        if (current_level > 1) {
            double min_elim_ms = 0.0;
            double max_elim_ms = 0.0;
            reduce_active_duration_bounds_ms(
                tree->comm,
                level_print_rank,
                level.is_process_active,
                elim_duration,
                min_elim_ms,
                max_elim_ms);
            
            if (verbose && rank == level_print_rank) {
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
        clear_ghosts(level);
        
        

        
        // ===== Step 4: Build parent level interactions =====
        
        // Special handling for level 1: set all DOFs as skeleton (no elimination)
        if (current_level == 1 && level.is_process_active) {
            if (verbose && rank == level_print_rank) {
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
        
        if (verbose && rank == level_print_rank) {
            std::cout << "  Level transition time: " << transition_duration.count() << " ms" << std::endl;
            std::cout << "  Parent boxes created: " << parent_boxes.size() << std::endl;
        }
        
        // ===== Step 5: Handle process reduction =====
        
        auto& parent_level = tree->levels[current_level - 1];
        
        // Error check: Level 2 → 1 should NOT trigger reduction
        if (current_level == 2 && level.parent_level_owner != rank && level.is_process_active) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Level 2 → 1 triggered process reduction! "
                "This is not allowed. parent_level_owner = " + std::to_string(level.parent_level_owner) +
                ", current rank = " + std::to_string(rank));
        }


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
            MPI_Barrier(tree->comm);
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

            if (verbose) {
                std::cout << "  Process " << rank << " sending " << parent_boxes.size() 
                          << " parent boxes to process " << level.parent_level_owner << std::endl;
            }
        }
        
        // Clear modified interaction matrices to free memory
        if (level.is_process_active) {
            clear_modified_interaction_matrices(level);
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
            MPI_Allreduce(local_deltas, max_deltas, 2, MPI_DOUBLE, MPI_MAX, tree->comm);
            total_data_exchange_time += std::chrono::duration_cast<clock::duration>(
                std::chrono::duration<double, std::milli>(max_deltas[0]));
            total_reduction_time += std::chrono::duration_cast<clock::duration>(
                std::chrono::duration<double, std::milli>(max_deltas[1]));
        }

        auto level_end = std::chrono::high_resolution_clock::now();
        auto level_duration = std::chrono::duration_cast<std::chrono::milliseconds>(level_end - level_start);
        
        if (verbose && rank == level_print_rank) {
            std::cout << "  Total level time: " << level_duration.count() << " ms" << std::endl;
        }
        
    }

    int64_t global_max_skel = 0;
    MPI_Allreduce(&local_max_skel, &global_max_skel, 1, MPI_INT64_T, MPI_MAX, tree->comm);
    if (out_rankmax) *out_rankmax = global_max_skel;

    
    // ===== Special handling for level 0 (root) =====
    
    const int root_print_rank = smallest_active_rank(tree->levels[0]);
    if (verbose && rank == root_print_rank) {
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
        
        if (verbose && rank == root_print_rank) {
            std::cout << "  Root box points: " << root_box.num_points << std::endl;
        }
        
        // At level 0, the assembled matrix is just the schur complement
        if (!root_box.schur_complement.is_allocated()) {
            throw std::runtime_error(
                "hierarchical_factorization_parallel: Root box schur complement not allocated");
        }
        
        int64_t n = root_box.schur_complement.rows;
        
        if (verbose && rank == root_print_rank) {
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
            
            if (verbose && rank == root_print_rank) {
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

            if (verbose && rank == root_print_rank) {
                std::cout << "  ✓ Root LU factorization complete" << std::endl;
            }

        } else if (factorization_method == FactorizationMethod::COMPLEX_SYM) {
            if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                root_box.X_RR.data = root_box.schur_complement.data;
                root_box.X_RR_pivots.resize(static_cast<size_t>(n));

                char uplo = 'L';
                int nn = n;
                int lwork = -1;
                int info = 0;
                std::vector<DataType> work(1);
                zsytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), work.data(), &lwork, &info);
                lwork = static_cast<int>(work[0].real());
                work.resize(static_cast<size_t>(lwork));
                zsytrf_(&uplo, &nn, root_box.X_RR.data.data(), &nn,
                        root_box.X_RR_pivots.data(), work.data(), &lwork, &info);

                if (info != 0) {
                    throw std::runtime_error(
                        "hierarchical_factorization_parallel: Bunch-Kaufman factorization of root failed at pivot " +
                        std::to_string(info));
                }

                root_box.X_RR.format = MatrixStorage<DataType>::BUNCH_KAUFMAN;

                if (verbose && rank == root_print_rank) {
                    std::cout << "  ✓ Root Bunch-Kaufman factorization complete" << std::endl;
                }
            } else {
                throw std::runtime_error("COMPLEX_SYM factorization only supported for complex<double>");
            }

        } else {
            // No factorization: just copy schur complement to X_RR
            root_box.X_RR.data = root_box.schur_complement.data;
            root_box.X_RR.format = MatrixStorage<DataType>::FULL;
            root_box.X_RR_pivots.clear();
            
            if (verbose && rank == root_print_rank) {
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

    const double data_exchange_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_data_exchange_time).count();
    const double reduction_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
            total_reduction_time).count();

    if (verbose && rank == root_print_rank) {
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - total_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "✓ Hierarchical Factorization Complete" << std::endl;
        std::cout << "  total time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "  data exchange communication time: " << std::llround(data_exchange_ms) << " ms" << std::endl;
        std::cout << "  process reduction communication time: " << std::llround(reduction_ms) << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    restore_base_process_affinity();
    clear_runtime_fmm_thread_count();
    destroy_dynamic_threading_context(dynamic_threading);


    size_t local_memory_usage = 0;
    for (int current_level = leaf_level; current_level >= 1; current_level--) {
    
        auto& level = tree->levels[current_level];
        for (auto& box : level.local_boxes) {
            local_memory_usage += calculate_box_data_size(box);
        }
    }
    printf("factorization memory usage on rank %d: %.2f GB\n", rank, local_memory_usage / (1024.0 * 1024.0 * 1024.0));
    fflush(stdout);
    
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
    bool verbose = true) {
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
            verbose);   // instantiated ONLY for double types
    } else {
        throw std::runtime_error("H2/FMM only supports double / std::complex<double>");
    }
}

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
void butterfly_factorization_parallel(H2<CoordType,DataType>* solver, double* factorization_time, double* entryeval_time) {

  if (solver->factorized) return;

  int rank = 0;
  MPI_Comm_rank(solver->comm, &rank);

  const auto factorization_method =
    (butterfly::is_complex_v<DataType>)
      ? fmm::FactorizationMethod::COMPLEX_SYM
      : fmm::FactorizationMethod::LU;
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
    true);
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
template<typename CoordType, typename DataType>
void hierarchical_solve_parallel(
    ParallelTree<CoordType, DataType>* tree,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<SolveDataRequest<CoordType, DataType>>> &solve_data,
    bool verbose = true) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
    
    if (verbose && rank == solve_header_rank) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Hierarchical Solve (Parallel MPI)" << std::endl;
        std::cout << "========================================" << std::endl;
        if (dynamic_threading.enabled) {
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

    
    if (verbose && rank == solve_header_rank) {
        std::cout << "Initialized solve data structures" << std::endl;
    }
    
    // ===== Phase 1: Forward Sweep (V^{-1}) with level transitions =====
    
    if (verbose && rank == solve_header_rank) {
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
        if (verbose && rank == level_print_rank) {
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
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Forward sweep time: " << forward_duration.count() << " ms" << std::endl;
    }
    
    // ===== Phase 2: Diagonal Solve (D^{-1}) =====
    
    if (verbose && rank == solve_root_rank) {
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
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Diagonal solves: " << num_diagonal_solves  << " from rank: " << rank << std::endl;
        std::cout << "  Diagonal solve time: " << diagonal_duration.count() << " ms" << std::endl;
    }
    
    
    // ===== Phase 3: Backward Sweep (W^{-1}) with level transitions =====
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "\n===== Phase 3: Backward Sweep (W^{-1}) =====" << std::endl;
    }
    
    auto backward_start = std::chrono::high_resolution_clock::now();
    
    for (int level = 1; level <= leaf_level; level++) {
        auto& tree_level = tree->levels[level];
        const int level_print_rank = smallest_active_rank(tree_level);
        
        if (verbose && rank == level_print_rank) {
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
        
        if (verbose && rank == level_print_rank) {
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
    
    if (verbose && rank == solve_root_rank) {
        std::cout << "  Backward sweep time: " << backward_duration.count() << " ms" << std::endl;
    }
    
    
    if (verbose && rank == solve_root_rank) {
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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


// //provide booleans for whether to do certain checks
// inline int h2_verification(H2<CoordType, DataType>* H2_solver, const ProgramOptions& options, bool verify_solution, bool verify_factorization, bool condition_number) {
//   // pass in tree and h2 options

//   // Step 1: verify solution x from Ax=b is reasonable
//   auto grid_points = fmm::build_grid_points_from_tree(tree.get(), options.grid_size);
//   int can_run_fft_verification = 0;
//   if (rank == 0 && !grid_points.empty()) {
//     can_run_fft_verification = 1;
//   }
//   MPI_Bcast(&can_run_fft_verification, 1, MPI_INT, 0, MPI_COMM_WORLD);
//   if (can_run_fft_verification) {
//     const auto verify_fft_start = std::chrono::high_resolution_clock::now();
//     const auto solve_metrics = fmm::verify_solution_fft(
//       tree.get(),
//       aggregated_rhs,
//       solution,
//       grid_points,
//       options.grid_size,
//       options.kernel_kind,
//       options.dimension,
//       options.wave_divisor,
//       options.length_scale,
//       options.nugget,
//       options.kappa,
//       true);
//     const auto verify_fft_end = std::chrono::high_resolution_clock::now();
//     const auto verify_fft_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//       verify_fft_end - verify_fft_start);

//     if (rank == 0) {
//       std::cout << "FFT verification time: " << verify_fft_duration.count() << " ms" << std::endl;
//       fflush(stdout);
//     }

//     // verify Forward error via multiply: ||F*x - A*x|| / ||A*x|| (no condition number amplification)
//     const double forward_relative_error = fmm::verify_forward_error_fft(
//       tree.get(),
//       grid_points,
//       options.grid_size,
//       options.kernel_kind,
//       options.dimension,
//       options.wave_divisor,
//       options.length_scale,
//       options.nugget,
//       options.kappa,
//       true);
//     if (rank == 0) {
//       const std::streamsize fwd_old_precision = std::cout.precision();
//       const std::ios::fmtflags fwd_old_flags = std::cout.flags();
//       std::cout << "Forward relative error: " << std::scientific
//                 << std::setprecision(std::numeric_limits<double>::max_digits10)
//                 << forward_relative_error
//                 << std::setprecision(fwd_old_precision);
//       std::cout.flags(fwd_old_flags);
//       std::cout << std::endl;
//       fflush(stdout);
//     }

//     // verify solution x from Ax=b is reasonable, this is the ground truth comparison for smaller matrices
//     if (options.N <= 4096) {
//       const double direct_error = fmm::verify_solution_direct(
//         &kernel,
//         aggregated_rhs,
//         solution,
//         grid_points,
//         tree->num_points,
//         options.dimension,
//         true);
//       if (rank == 0) {
//         std::cout << "verifying with direct matrix vector multiply since N <= 4096" << std::endl;
//         std::cout << "Direct error: " << direct_error << std::endl;
//       }
//     }

//     if (rank == 0) {
//       const std::streamsize old_precision = std::cout.precision();
//       const std::ios::fmtflags old_flags = std::cout.flags();
//       std::cout << "\nBackward residual: " << std::scientific
//                 << std::setprecision(std::numeric_limits<double>::max_digits10)
//                 << solve_metrics.relative_residual
//                 << std::setprecision(old_precision);
//       std::cout.flags(old_flags);
//       std::cout << std::endl;
//     }


//     // Estimate Condition Number
//     if (options.cond_samples > 0) {
//       const double kappa_est = fmm::estimate_condition_number(
//         tree.get(),
//         grid_points,
//         options.grid_size,
//         options.kernel_kind,
//         options.dimension,
//         options.cond_samples,
//         options.wave_divisor,
//         options.length_scale,
//         options.nugget,
//         options.kappa);
//       if (rank == 0) {
//         const std::streamsize be_old_precision = std::cout.precision();
//         const std::ios::fmtflags be_old_flags = std::cout.flags();
//         std::cout << "Backward relative error (rough estimate, uses ||A|| approx ||Ax||/||x||): "
//                   << std::scientific
//                   << std::setprecision(std::numeric_limits<double>::max_digits10)
//                   << solve_metrics.backward_error
//                   << std::setprecision(be_old_precision);
//         std::cout.flags(be_old_flags);
//         std::cout << std::endl;
//         std::cout << "Estimated condition number: " << std::scientific
//                   << std::setprecision(6) << kappa_est << std::endl;
//       }
//     }
//   }

//   return 0;
// }
}