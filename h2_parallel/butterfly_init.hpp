#pragma once

#include "butterfly_types.hpp"

namespace butterfly {
using namespace fmm;

inline int default_num_levels(int64_t grid_size) {
  int k = 0;
  while (grid_size > 1) {  
      grid_size /= 2;
      ++k;
  }
  return k;
}

inline int calc_num_levels(int64_t grid_size, int64_t grid_dim_min_leaf) {
    int k = 0;
    while (grid_dim_min_leaf <= grid_size) {
        grid_dim_min_leaf *= 2;
        ++k;
    }
    return k;
}


inline int64_t default_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 256 : 4096;
}

inline int64_t min_reduction_threshold_for_dimension(int dimension) {
    return (dimension == 2) ? 16 : 64;
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
  double tolerance, int64_t reduction_threshold, int64_t Nmin_leaf) {

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

    // for h2: get nmin_leaf argument and validate that it's possible
    int64_t grid_dim_min_leaf = std::llround(std::pow((double)Nmin_leaf, 1.0 / *Ndim));
    if (grid_dim_min_leaf == 0) {
        h2_options.num_levels = default_num_levels(grid_size);
    } else if (grid_dim_min_leaf >= 2 && *Npo >= Nmin_leaf) {
        h2_options.num_levels = calc_num_levels(grid_size, grid_dim_min_leaf);
    } else {
        throw std::invalid_argument(
            "Nmin_leaf must be 0 (default), >= 2^dimension, or <= Npo for H2 solver."
        );
    }

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
    h2_options.reduction_threshold = reduction_threshold;
    if (h2_options.reduction_threshold <= 0) {
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
  if (rank == 0 && options.verbosity >= 1) {
    std::cout << "=== Hierarchical Factorization Test ("
              << options.dimension << "D "
              << number_kind_to_string(options.number_kind) << ") ===" << std::endl;
  
    if (const int dynamic_cpu_cap =
            fmm::parse_positive_thread_count(std::getenv("FMM_MAX_CPUS_PER_NODE"));
        dynamic_cpu_cap > 0) {
        std::cout << "Dynamic thread cpu cap per node: " << dynamic_cpu_cap << std::endl;
    }
  
  }
  (void) fmm::base_process_cpu_list();

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

  // new2old changed to index-1 based for fortran side
  for (int& v : new2old) v += 1;


  // To Do: maybe put this section and below to c_bpack_factor
  return 0;
}


} // namespace butterfly
