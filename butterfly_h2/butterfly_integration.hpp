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

template <typename T> struct is_complex : std::false_type {};                 // (1) primary
template <typename T> struct is_complex<std::complex<T>> : std::true_type {}; // (2) specialization
template <typename T> inline constexpr bool is_complex_v = is_complex<T>::value; // (3) alias


// construct b_mat
template<typename CoordType, typename DataType>
struct H2 {
    //using CoordType = double;
    MPI_Comm comm;
    int64_t N;
    int dimension;
    //temporary comment
    //fmm::HierarchicalFactorization<CoordType, C_DT, KernelType> factorizer;
    
    //temporary comment
    void (*kernel)(int*, int*, DataType*, void*);
    std::unique_ptr<fmm::ParallelTree<CoordType, DataType>> tree;
    ProgramOptions options;
    //std::unique_ptr<ButterflyPACKKernel> kernel;

    //temporary comment
    //RedistributionPlan redistribution;
    bool factorized = false;
};

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

// initiate MPI, first step. requires mpi library
// void MPI_init(int argc, char* argv[]){
//     MPI_Init(&argc, &argv);

//     int rank = 0; // ID of calling process
//     int size = 1; // total number of processes
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
    
//     if (rank == 0) {
//         std::cout << "Thread runtime: fmm_threads=" << fmm::configured_fmm_thread_count()
//                   << ", omp_max_threads=" << omp_get_max_threads()
//                   << ", openblas_threads=" << openblas_get_num_threads()
//                   << ", visible_cpus=" << fmm::visible_process_cpu_count();
//         if (const int dynamic_cpu_cap =
//                 fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
//             dynamic_cpu_cap > 0) {
//             std::cout << ", dynamic_cpu_cap_per_node=" << dynamic_cpu_cap;
//         }
//         std::cout << std::endl;
//     }

//     {
//         const auto& base_cpus = fmm::base_process_cpu_list();
//         const int cap = fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
//         MPI_Comm shared_comm;
//         MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
//                             MPI_INFO_NULL, &shared_comm);
//         int shared_rank = 0, shared_size = 1;
//         MPI_Comm_rank(shared_comm, &shared_rank);
//         MPI_Comm_size(shared_comm, &shared_size);
//         MPI_Comm_free(&shared_comm);

//         int slice_first = -1, slice_last = -1, slice_count = 0;
//         if (!base_cpus.empty()) {
//             const int usable = std::min<int>(cap > 0 ? cap : (int)base_cpus.size(),
//                                               (int)base_cpus.size());
//             const int threads = std::max(1, usable / std::max(1, shared_size));
//             int slice_begin_idx = shared_rank * threads;
//             if (slice_begin_idx >= usable) slice_begin_idx = std::max(0, usable - 1);
//             const int slice_end_idx = std::min(slice_begin_idx + threads, usable);
//             if (slice_end_idx > slice_begin_idx) {
//                 slice_first = base_cpus[slice_begin_idx];
//                 slice_last = base_cpus[slice_end_idx - 1];
//                 slice_count = slice_end_idx - slice_begin_idx;
//             }
//         }

//         constexpr int kSocketCores = 64;
//         constexpr int kCcdCores = 8;
//         const bool socket_ok = (slice_first >= 0) &&
//             ((slice_count >= kSocketCores && slice_first % kSocketCores == 0 && slice_count % kSocketCores == 0) ||
//              (slice_count < kSocketCores && (slice_first / kSocketCores == slice_last / kSocketCores)));
//         const bool ccd_ok = (slice_first >= 0) &&
//             ((slice_count >= kCcdCores && slice_first % kCcdCores == 0 && slice_count % kCcdCores == 0) ||
//              (slice_count < kCcdCores && (slice_first / kCcdCores == slice_last / kCcdCores)));

//         struct AffinityInfo { int rank; int shared_rank; int first; int last; int count; int socket_ok; int ccd_ok; };
//         AffinityInfo local{rank, shared_rank, slice_first, slice_last, slice_count,
//                             socket_ok ? 1 : 0, ccd_ok ? 1 : 0};
//         std::vector<AffinityInfo> all(size);
//         MPI_Gather(&local, sizeof(AffinityInfo), MPI_BYTE,
//                    all.data(), sizeof(AffinityInfo), MPI_BYTE,
//                    0, MPI_COMM_WORLD);
//         if (rank == 0) {
//             std::cout << "Projected per-rank CPU slice (Perlmutter: 64 cores/socket, 8 cores/CCD):" << std::endl;
//             bool any_socket_bad = false, any_ccd_bad = false;
//             for (const auto& a : all) {
//                 std::cout << "  rank " << a.rank
//                           << " (node-local " << a.shared_rank << ")"
//                           << ": cpus=[" << a.first << "," << a.last << "]"
//                           << " count=" << a.count
//                           << (a.socket_ok ? " socket-ok" : " SOCKET-CROSSES")
//                           << (a.ccd_ok ? " ccd-ok" : " ccd-crosses")
//                           << std::endl;
//                 if (!a.socket_ok) any_socket_bad = true;
//                 if (!a.ccd_ok) any_ccd_bad = true;
//             }
//             if (any_socket_bad) {
//                 std::cout << "  WARNING: at least one rank's slice crosses a socket boundary"
//                           << " — NUMA traffic expected." << std::endl;
//             } else if (any_ccd_bad) {
//                 std::cout << "  NOTE: slices cross CCD boundaries (cross-L3 within socket is OK but not optimal)."
//                           << std::endl;
//             }
//         }
//     }
// }

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

template<typename CoordType, typename DataType>
int h2_initiate(H2<CoordType, DataType>* H2_solver, const ProgramOptions& options, CoordType* Locations, int rank) {
  // from H2 struct, set typename CoordType, DataType
  // from H2 struct, get kernel
  // from H2 struct, get rank, size
  if (rank == 0) {
    std::cout << "=== Hierarchical Factorization Test ("
              << options.dimension << "D "
              << kernel_kind_to_string(options.kernel_kind) << " kernel, "
              << number_kind_to_string(options.number_kind) << ") ===" << std::endl;
  }

  double bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, (options.dimension == 2) ? 0.0 : 1.0};

  // Provide Locations to nullptr, check compatibility
  auto tree = std::unique_ptr<fmm::ParallelTree<CoordType, DataType>>(
      fmm::create_uniform_tree<CoordType, DataType>(
          Locations,
          options.N,
          options.num_levels,
          bounds,
          options.dimension,
          H2_solver->comm,
          options.reduction_threshold));

  H2_solver->tree = std::move(tree); 

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
template<typename CoordType, typename DataType>
void hierarchical_factorization_parallel(
    ParallelTree<CoordType, DataType>* tree,
    void (*kernel)(int*, int*, DataType*, void*),
    double tolerance,
    bool is_symmetric,
    bool is_hermitian,
    FactorizationMethod factorization_method,
    const std::vector<CoordType>& unit_proxy_points,
    int num_proxy,
    CoordType proxy_radius,
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
            elim_duration = std::chrono::duration_cast<std::chrono::milliseconds>(elim_end - elim_start);

            for (const auto& box : level.local_boxes) {
                total_skeleton += box.skeleton_indices.size();
                total_redundant += box.redundant_indices.size();
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


// provide booleans for whether to do certain checks
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