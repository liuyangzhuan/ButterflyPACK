#pragma once

// C++ Standard Library
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <cmath>
#include <iostream>

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
    //void (*kernel)(int*, int*, C_DT*,C2Fptr);
    std::unique_ptr<fmm::ParallelTree<CoordType, DataType>> tree;
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
int h2_initiate(H2<CoordType, DataType>* H2_solver, const ProgramOptions& options, double* Locations, int rank) {
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
  // const int leaf_level = options.num_levels - 1;

  // const auto factorization_method =
  //     (options.number_kind == NumberKind::COMPLEX)
  //         ? fmm::FactorizationMethod::COMPLEX_SYM
  //         : fmm::FactorizationMethod::LU;
  // fmm::HierarchicalFactorization<CoordType, DataType, KernelType> factorizer(
  //     options.N,
  //     fmm::MatrixProperty::SYMMETRIC,
  //     &kernel,
  //     options.dimension,
  //     factorization_method,
  //     options.num_proxy);

  // const auto& unit_proxy = factorizer.get_unit_proxy_points();
  // const int num_proxy = factorizer.get_num_proxy_points();

  // const bool is_symmetric = true;
  // const bool is_hermitian = false;

  // MPI_Barrier(MPI_COMM_WORLD);
  return 0;
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