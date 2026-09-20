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


inline int64_t min_reduction_threshold_for_dimension(int dimension) {
    if (dimension < 1 || dimension > 3) {
        throw std::invalid_argument("dimension must be 1, 2, or 3");
    }
    return int64_t{1} << dimension;
}

inline int64_t default_reduction_threshold_for_dimension(int dimension) {
    return min_reduction_threshold_for_dimension(dimension);
}

inline int default_num_proxy_for_dimension(int dimension) {
    switch (dimension) {
        case 1: return 2;
        case 2: return 32;
        case 3: return 256;
        default: throw std::invalid_argument("dimension must be 1, 2, or 3");
    }
}

/**
 * @brief Derive global domain bounds from the point coordinates.
 *
 * Locations is point-major (interleaved), stride = dimension:
 *   [x0,y0,z0, x1,y1,z1, ...]     (matches tree_impl.hpp:801-803)
 *
 * Fills bounds as [xmin,xmax, ymin,ymax, zmin,zmax]. Entries for inactive
 * dimensions are zero.
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
    if (dimension < 1 || dimension > 3) {
        throw std::invalid_argument(
            "compute_global_bounds: dimension must be 1, 2, or 3");
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

    for (int d = 0; d < 3; ++d) {
        if (d < dimension) {
            bounds[2 * d] = lo[d];
            bounds[2 * d + 1] = hi[d];
        } else {
            bounds[2 * d] = CoordType(0);
            bounds[2 * d + 1] = CoordType(0);
        }
    }
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
inline ProgramOptions parse_program_options64(
  int64_t Npo, int Ndim, double tolerance,
  int64_t reduction_threshold, int64_t Nmin_leaf,
  int CA_level = 10000) {

    ProgramOptions h2_options;
    
    h2_options.N = Npo;
    if (h2_options.N <= 0) {
        throw std::invalid_argument("N must be positive.");
    }

    h2_options.dimension = Ndim;
    if (h2_options.dimension < 1 || h2_options.dimension > 3) {
        throw std::invalid_argument("H2 solver dimension must be 1, 2, or 3.");
    }

    // Approximate per-dimension grid size for the existing level-selection logic.
    int64_t grid_size = std::llround(std::pow((double)Npo, 1.0 / Ndim));

    h2_options.grid_size = grid_size;

    if (h2_options.grid_size < 2) {
      throw std::invalid_argument(
        "grid size must be >=2 for H2 solver."
      );
    }

    // for h2: get nmin_leaf argument and validate that it's possible
    int64_t grid_dim_min_leaf = std::llround(std::pow((double)Nmin_leaf, 1.0 / Ndim));
    if (grid_dim_min_leaf == 0) {
        h2_options.num_levels = default_num_levels(grid_size);
    } else if (grid_dim_min_leaf >= 2 && Npo >= Nmin_leaf) {
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
    h2_options.CA_level =
        (h2_options.dimension == 1 || CA_level < 0)
        ? h2_options.num_levels
        : CA_level;
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

inline ProgramOptions parse_program_options(
  int* Npo, int* Ndim, double*, double tolerance,
  int64_t reduction_threshold, int64_t Nmin_leaf,
  int CA_level = 10000) {
    if (Npo == nullptr || Ndim == nullptr) {
        throw std::invalid_argument("parse_program_options: null size argument");
    }
    return parse_program_options64(
        static_cast<int64_t>(*Npo), *Ndim, tolerance,
        reduction_threshold, Nmin_leaf, CA_level);
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


template<typename CoordType>
void compute_distributed_global_bounds(
    MPI_Comm comm,
    const CoordType* local_coordinates,
    int64_t local_count,
    int dimension,
    const CoordType* supplied_bounds,
    bool bounds_provided,
    CoordType bounds[6]) {
  if (dimension < 1 || dimension > 3 || local_count < 0) {
    throw std::invalid_argument(
        "compute_distributed_global_bounds: invalid dimension or local count");
  }
  if (local_count > 0 && local_coordinates == nullptr) {
    throw std::invalid_argument(
        "compute_distributed_global_bounds: local coordinates are null");
  }

  if (bounds_provided) {
    if (supplied_bounds == nullptr) {
      throw std::invalid_argument(
          "compute_distributed_global_bounds: supplied bounds are null");
    }
    for (int d = 0; d < 3; ++d) {
      if (d >= dimension) {
        bounds[2 * d] = CoordType(0);
        bounds[2 * d + 1] = CoordType(0);
        continue;
      }
      bounds[2 * d] = supplied_bounds[2 * d];
      bounds[2 * d + 1] = supplied_bounds[2 * d + 1];
      if (!(bounds[2 * d] < bounds[2 * d + 1]) ||
          !std::isfinite(bounds[2 * d]) ||
          !std::isfinite(bounds[2 * d + 1])) {
        throw std::invalid_argument(
            "compute_distributed_global_bounds: invalid supplied bounds");
      }
    }
    return;
  }

  CoordType local_min[3] = {
      std::numeric_limits<CoordType>::max(),
      std::numeric_limits<CoordType>::max(),
      std::numeric_limits<CoordType>::max()};
  CoordType local_max[3] = {
      std::numeric_limits<CoordType>::lowest(),
      std::numeric_limits<CoordType>::lowest(),
      std::numeric_limits<CoordType>::lowest()};
  for (int64_t point = 0; point < local_count; ++point) {
    for (int d = 0; d < dimension; ++d) {
      const CoordType value =
          local_coordinates[point * dimension + d];
      if (!std::isfinite(value)) {
        throw std::invalid_argument(
            "compute_distributed_global_bounds: non-finite coordinate");
      }
      local_min[d] = std::min(local_min[d], value);
      local_max[d] = std::max(local_max[d], value);
    }
  }

  CoordType global_min[3];
  CoordType global_max[3];
  const MPI_Datatype coordinate_type =
      std::is_same_v<CoordType, float> ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Allreduce(local_min, global_min, 3, coordinate_type, MPI_MIN, comm);
  MPI_Allreduce(local_max, global_max, 3, coordinate_type, MPI_MAX, comm);
  for (int d = 0; d < 3; ++d) {
    if (d >= dimension) {
      bounds[2 * d] = CoordType(0);
      bounds[2 * d + 1] = CoordType(0);
      continue;
    }
    if (!(global_min[d] <= global_max[d])) {
      throw std::runtime_error(
          "compute_distributed_global_bounds: no global coordinates");
    }
    const CoordType span = global_max[d] - global_min[d];
    const CoordType pad =
        (span > CoordType(0) ? span : CoordType(1)) * CoordType(1e-6);
    bounds[2 * d] = global_min[d] - pad;
    bounds[2 * d + 1] = global_max[d] + pad;
  }
}

inline std::vector<int> distributed_displacements(
    const std::vector<int>& counts,
    const char* caller) {
  std::vector<int> displacements(counts.size(), 0);
  int64_t total = 0;
  for (size_t process = 0; process < counts.size(); ++process) {
    if (counts[process] < 0 || total > std::numeric_limits<int>::max()) {
      throw std::overflow_error(
          std::string(caller) + ": MPI displacement exceeds INT_MAX");
    }
    displacements[process] = static_cast<int>(total);
    total += counts[process];
  }
  if (total > std::numeric_limits<int>::max()) {
    throw std::overflow_error(
        std::string(caller) + ": local MPI count exceeds INT_MAX");
  }
  return displacements;
}

template<typename CoordType, typename DataType>
int distributed_leaf_owner(
    const fmm::TreeLevel<CoordType, DataType>& leaf,
    int64_t morton) {
  const int64_t boxes = leaf.num_boxes_global;
  const int processes = leaf.num_active_processes;
  if (morton < 0 || morton >= boxes || processes <= 0 || boxes < processes) {
    throw std::runtime_error(
        "distributed_leaf_owner: unsupported leaf box distribution");
  }
  const int64_t quotient = boxes / processes;
  const int64_t remainder = boxes % processes;
  const int64_t enlarged_span = remainder * (quotient + 1);
  const int region = morton < enlarged_span
      ? static_cast<int>(morton / (quotient + 1))
      : static_cast<int>(remainder +
            (morton - enlarged_span) / quotient);
  const auto owner = leaf.morton_to_rank.find(region);
  if (owner == leaf.morton_to_rank.end()) {
    throw std::runtime_error(
        "distributed_leaf_owner: missing Morton-region owner");
  }
  return owner->second;
}

template<typename CoordType, typename DataType>
void assign_distributed_points(
    fmm::ParallelTree<CoordType, DataType>* tree,
    int64_t global_count,
    int64_t input_local_count,
    const int64_t* input_global_ids,
    const CoordType* input_coordinates,
    DistributedLayout64& layout) {
  if (tree == nullptr || global_count <= 0 || input_local_count < 0) {
    throw std::invalid_argument(
        "assign_distributed_points: invalid point count");
  }
  if (input_local_count > std::numeric_limits<int>::max()) {
    throw std::overflow_error(
        "assign_distributed_points: local point count currently must fit INT_MAX");
  }
  if (input_local_count > 0 &&
      (input_global_ids == nullptr || input_coordinates == nullptr)) {
    throw std::invalid_argument(
        "assign_distributed_points: null local point data");
  }

  int64_t global_input_count = 0;
  MPI_Allreduce(&input_local_count, &global_input_count, 1, MPI_INT64_T,
                MPI_SUM, tree->comm);
  if (global_input_count != global_count) {
    throw std::invalid_argument(
        "assign_distributed_points: sum of local counts does not equal N_global");
  }

  uint64_t local_sum = 0;
  uint64_t local_sum_squares = 0;
  for (int64_t point = 0; point < input_local_count; ++point) {
    const int64_t id = input_global_ids[point];
    if (id < 1 || id > global_count) {
      throw std::invalid_argument(
          "assign_distributed_points: global IDs must be in [1,N_global]");
    }
    const uint64_t value = static_cast<uint64_t>(id);
    local_sum += value;
    local_sum_squares += value * value;
  }
  uint64_t global_sum = 0;
  uint64_t global_sum_squares = 0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_UINT64_T, MPI_SUM,
                tree->comm);
  MPI_Allreduce(&local_sum_squares, &global_sum_squares, 1, MPI_UINT64_T,
                MPI_SUM, tree->comm);
  const unsigned __int128 n = static_cast<uint64_t>(global_count);
  const uint64_t expected_sum = static_cast<uint64_t>(n * (n + 1) / 2);
  const uint64_t expected_sum_squares = static_cast<uint64_t>(
      n * (n + 1) * (2 * n + 1) / 6);
  if (global_sum != expected_sum ||
      global_sum_squares != expected_sum_squares) {
    throw std::invalid_argument(
        "assign_distributed_points: global IDs are not a permutation of [1,N_global]");
  }

  const int leaf_level_number = tree->num_levels - 1;
  auto& leaf = tree->levels[leaf_level_number];
  struct OutgoingPoint {
    int64_t id;
    int64_t input_slot;
    int64_t morton;
    std::array<CoordType, 3> coordinate;
  };
  std::vector<std::vector<OutgoingPoint>> outgoing(
      static_cast<size_t>(tree->mpi_size));
  layout.input_global_ids.clear();
  if (input_local_count > 0) {
    layout.input_global_ids.assign(
        input_global_ids, input_global_ids + input_local_count);
  }
  layout.input_internal_owner.assign(
      static_cast<size_t>(input_local_count), -1);

  for (int64_t point = 0; point < input_local_count; ++point) {
    CoordType coordinate[3] = {CoordType(0), CoordType(0), CoordType(0)};
    for (int d = 0; d < tree->dimension; ++d) {
      coordinate[d] = input_coordinates[point * tree->dimension + d];
    }
    int32_t grid_coordinates[3] = {0, 0, 0};
    const int64_t morton = fmm::point_to_morton(
        coordinate, tree->dimension, tree->global_bounds,
        leaf_level_number, grid_coordinates);
    const int owner = distributed_leaf_owner(leaf, morton);
    OutgoingPoint record;
    record.id = input_global_ids[point] - 1;
    record.input_slot = point;
    record.morton = morton;
    record.coordinate = {coordinate[0], coordinate[1], coordinate[2]};
    outgoing[static_cast<size_t>(owner)].push_back(record);
    layout.input_internal_owner[static_cast<size_t>(point)] = owner;
  }

  layout.forward_send_counts.assign(static_cast<size_t>(tree->mpi_size), 0);
  for (int process = 0; process < tree->mpi_size; ++process) {
    const size_t count = outgoing[static_cast<size_t>(process)].size();
    if (count > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error(
          "assign_distributed_points: per-peer point count exceeds INT_MAX");
    }
    layout.forward_send_counts[static_cast<size_t>(process)] =
        static_cast<int>(count);
  }
  layout.forward_send_displacements = distributed_displacements(
      layout.forward_send_counts, "assign_distributed_points");
  layout.forward_recv_counts.assign(static_cast<size_t>(tree->mpi_size), 0);
  MPI_Alltoall(layout.forward_send_counts.data(), 1, MPI_INT,
               layout.forward_recv_counts.data(), 1, MPI_INT, tree->comm);
  layout.forward_recv_displacements = distributed_displacements(
      layout.forward_recv_counts, "assign_distributed_points");

  const int send_total = input_local_count == 0
      ? 0
      : layout.forward_send_displacements.back() +
            layout.forward_send_counts.back();
  const int recv_total = layout.forward_recv_counts.empty()
      ? 0
      : layout.forward_recv_displacements.back() +
            layout.forward_recv_counts.back();
  std::vector<int64_t> send_metadata(static_cast<size_t>(send_total) * 3);
  std::vector<CoordType> send_coordinates(
      static_cast<size_t>(send_total) * tree->dimension);
  layout.forward_send_input_slots.resize(static_cast<size_t>(send_total));
  for (int process = 0; process < tree->mpi_size; ++process) {
    const int displacement =
        layout.forward_send_displacements[static_cast<size_t>(process)];
    const auto& records = outgoing[static_cast<size_t>(process)];
    for (size_t offset = 0; offset < records.size(); ++offset) {
      const size_t packet = static_cast<size_t>(displacement) + offset;
      send_metadata[3 * packet] = records[offset].id;
      send_metadata[3 * packet + 1] = records[offset].input_slot;
      send_metadata[3 * packet + 2] = records[offset].morton;
      layout.forward_send_input_slots[packet] = records[offset].input_slot;
      for (int d = 0; d < tree->dimension; ++d) {
        send_coordinates[packet * tree->dimension + d] =
            records[offset].coordinate[static_cast<size_t>(d)];
      }
    }
  }
  outgoing.clear();

  auto scaled_mpi_arrays = [](const std::vector<int>& counts,
                              const std::vector<int>& displacements,
                              int scale, std::vector<int>& scaled_counts,
                              std::vector<int>& scaled_displacements) {
    scaled_counts.resize(counts.size());
    scaled_displacements.resize(displacements.size());
    for (size_t process = 0; process < counts.size(); ++process) {
      const int64_t count = static_cast<int64_t>(counts[process]) * scale;
      const int64_t displacement =
          static_cast<int64_t>(displacements[process]) * scale;
      if (count > std::numeric_limits<int>::max() ||
          displacement > std::numeric_limits<int>::max()) {
        throw std::overflow_error(
            "assign_distributed_points: scaled MPI count exceeds INT_MAX");
      }
      scaled_counts[process] = static_cast<int>(count);
      scaled_displacements[process] = static_cast<int>(displacement);
    }
  };

  std::vector<int64_t> recv_metadata(static_cast<size_t>(recv_total) * 3);
  std::vector<CoordType> recv_coordinates(
      static_cast<size_t>(recv_total) * tree->dimension);
  std::vector<int> send_counts_scaled;
  std::vector<int> send_displacements_scaled;
  std::vector<int> recv_counts_scaled;
  std::vector<int> recv_displacements_scaled;
  scaled_mpi_arrays(
      layout.forward_send_counts, layout.forward_send_displacements, 3,
      send_counts_scaled, send_displacements_scaled);
  scaled_mpi_arrays(
      layout.forward_recv_counts, layout.forward_recv_displacements, 3,
      recv_counts_scaled, recv_displacements_scaled);
  MPI_Alltoallv(
      send_metadata.data(), send_counts_scaled.data(),
      send_displacements_scaled.data(), MPI_INT64_T,
      recv_metadata.data(), recv_counts_scaled.data(),
      recv_displacements_scaled.data(), MPI_INT64_T, tree->comm);

  scaled_mpi_arrays(
      layout.forward_send_counts, layout.forward_send_displacements,
      tree->dimension, send_counts_scaled, send_displacements_scaled);
  scaled_mpi_arrays(
      layout.forward_recv_counts, layout.forward_recv_displacements,
      tree->dimension, recv_counts_scaled, recv_displacements_scaled);
  const MPI_Datatype coordinate_type =
      std::is_same_v<CoordType, float> ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Alltoallv(
      send_coordinates.data(), send_counts_scaled.data(),
      send_displacements_scaled.data(), coordinate_type,
      recv_coordinates.data(), recv_counts_scaled.data(),
      recv_displacements_scaled.data(), coordinate_type, tree->comm);

  struct ReceivedPoint {
    int64_t id;
    int64_t input_slot;
    int64_t morton;
    int source_rank;
    int receive_packet;
    std::array<CoordType, 3> coordinate;
  };
  std::vector<ReceivedPoint> received;
  received.reserve(static_cast<size_t>(recv_total));
  for (int source = 0; source < tree->mpi_size; ++source) {
    const int begin =
        layout.forward_recv_displacements[static_cast<size_t>(source)];
    const int end = begin +
        layout.forward_recv_counts[static_cast<size_t>(source)];
    for (int packet = begin; packet < end; ++packet) {
      ReceivedPoint record;
      record.id = recv_metadata[static_cast<size_t>(3 * packet)];
      record.input_slot = recv_metadata[static_cast<size_t>(3 * packet + 1)];
      record.morton = recv_metadata[static_cast<size_t>(3 * packet + 2)];
      record.source_rank = source;
      record.receive_packet = packet;
      record.coordinate = {CoordType(0), CoordType(0), CoordType(0)};
      for (int d = 0; d < tree->dimension; ++d) {
        record.coordinate[static_cast<size_t>(d)] =
            recv_coordinates[static_cast<size_t>(packet) * tree->dimension + d];
      }
      received.push_back(record);
    }
  }
  std::sort(received.begin(), received.end(),
            [](const ReceivedPoint& first, const ReceivedPoint& second) {
              if (first.morton != second.morton) {
                return first.morton < second.morton;
              }
              if (first.id != second.id) return first.id < second.id;
              if (first.source_rank != second.source_rank) {
                return first.source_rank < second.source_rank;
              }
              return first.input_slot < second.input_slot;
            });

  for (auto& box : leaf.local_boxes) {
    box.point_indices.clear();
    box.point_coords.clear();
    box.num_points = 0;
  }
  layout.internal_global_ids.resize(received.size());
  layout.forward_recv_internal_slots.resize(received.size());
  for (size_t internal_slot = 0; internal_slot < received.size();
       ++internal_slot) {
    const auto& record = received[internal_slot];
    if (record.morton < leaf.local_morton_start ||
        record.morton > leaf.local_morton_end) {
      throw std::runtime_error(
          "assign_distributed_points: received point belongs to a nonlocal box");
    }
    const int64_t local_box = record.morton - leaf.local_morton_start;
    auto& box = leaf.local_boxes[static_cast<size_t>(local_box)];
    box.point_indices.push_back(record.id);
    for (int d = 0; d < tree->dimension; ++d) {
      box.point_coords.push_back(
          record.coordinate[static_cast<size_t>(d)]);
    }
    ++box.num_points;
    layout.internal_global_ids[internal_slot] = record.id + 1;
    layout.forward_recv_internal_slots[
        static_cast<size_t>(record.receive_packet)] =
        static_cast<int64_t>(internal_slot);
  }
  for (const auto& box : leaf.local_boxes) {
    if (box.num_points == 0) {
      throw std::runtime_error(
          "assign_distributed_points: empty leaf boxes are not supported yet");
    }
  }

  layout.comm = tree->comm;
  layout.rank = tree->mpi_rank;
  layout.size = tree->mpi_size;
  layout.dimension = tree->dimension;
  layout.global_size = global_count;
  layout.input_local_size = input_local_count;
  layout.internal_local_size = static_cast<int64_t>(received.size());
  int64_t internal_offset = 0;
  MPI_Exscan(&layout.internal_local_size, &internal_offset, 1, MPI_INT64_T,
             MPI_SUM, tree->comm);
  if (tree->mpi_rank == 0) internal_offset = 0;
  layout.internal_global_start = internal_offset + 1;

  std::vector<int64_t> reverse_send(received.size());
  for (size_t packet = 0; packet < reverse_send.size(); ++packet) {
    reverse_send[packet] = layout.internal_global_start +
        layout.forward_recv_internal_slots[packet];
  }
  std::vector<int64_t> reverse_receive(
      layout.forward_send_input_slots.size());
  MPI_Alltoallv(
      reverse_send.data(), layout.forward_recv_counts.data(),
      layout.forward_recv_displacements.data(), MPI_INT64_T,
      reverse_receive.data(), layout.forward_send_counts.data(),
      layout.forward_send_displacements.data(), MPI_INT64_T, tree->comm);
  layout.input_internal_global_index.assign(
      static_cast<size_t>(input_local_count), 0);
  for (size_t packet = 0; packet < reverse_receive.size(); ++packet) {
    const int64_t input_slot = layout.forward_send_input_slots[packet];
    layout.input_internal_global_index[static_cast<size_t>(input_slot)] =
        reverse_receive[packet];
  }
}


template<typename CoordType, typename DataType>
int bpack_initiate_distributed64(
    H2<CoordType, DataType>* solver,
    const ProgramOptions& options,
    int64_t input_local_count,
    const int64_t* input_global_ids,
    const CoordType* input_coordinates,
    bool bounds_provided,
    const CoordType* supplied_bounds) {
  if (solver == nullptr) {
    throw std::invalid_argument(
        "bpack_initiate_distributed64: solver is null");
  }
  if (options.id_neighborhood_radius != 2 || options.id_proxy_mode != 0) {
    throw std::invalid_argument(
        "c_bpack_construct_init_distributed64 currently requires "
        "H2_ID_radius=2 and H2_ID_proxy=0");
  }

  int rank = 0;
  MPI_Comm_rank(solver->comm, &rank);
  if (rank == 0 && options.verbosity >= 1) {
    std::cout << "=== Hierarchical Factorization Test (distributed64, "
              << options.dimension << "D "
              << number_kind_to_string(options.number_kind) << ") ==="
              << std::endl;
  }
  (void)fmm::base_process_cpu_list();

  CoordType bounds[6];
  compute_distributed_global_bounds(
      solver->comm, input_coordinates, input_local_count,
      options.dimension, supplied_bounds, bounds_provided, bounds);

  std::vector<int> unused_index_map;
  auto tree = std::unique_ptr<fmm::ParallelTree<CoordType, DataType>>(
      fmm::create_uniform_tree<CoordType, DataType>(
          nullptr, 0, options.num_levels, bounds, options.dimension,
          solver->comm, options.reduction_threshold, options.CA_level,
          unused_index_map));
  tree->num_points = options.N;
  tree->id_neighborhood_radius = options.id_neighborhood_radius;
  tree->id_proxy_mode = options.id_proxy_mode;
  tree->id_proxy_points = options.id_proxy_points;
  tree->id_adaptive_batch = options.id_adaptive_batch;

  auto layout = std::make_unique<DistributedLayout64>();
  assign_distributed_points(
      tree.get(), options.N, input_local_count, input_global_ids,
      input_coordinates, *layout);

  solver->N = options.N;
  solver->dimension = options.dimension;
  solver->tree = std::move(tree);
  solver->distributed_layout = std::move(layout);
  return 0;
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
          options.CA_level,
          idx_map));

  H2_solver->tree = std::move(tree); 

  H2_solver->tree->id_neighborhood_radius = options.id_neighborhood_radius;
  H2_solver->tree->id_proxy_mode = options.id_proxy_mode;
  H2_solver->tree->id_proxy_points = options.id_proxy_points;
  H2_solver->tree->id_adaptive_batch = options.id_adaptive_batch;

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

  const bool needs_id_source_index =
      options.id_neighborhood_radius > 2 || options.id_proxy_mode != 0;
  if (needs_id_source_index) {
      H2_solver->tree->id_source_point_order.assign(
          new2old.begin(), new2old.end());

      const int leaf_level = H2_solver->tree->num_levels - 1;
      const auto& leaf = H2_solver->tree->levels[leaf_level];
      if (leaf.local_boxes.size() >
          static_cast<size_t>(std::numeric_limits<int>::max())) {
          throw std::runtime_error(
              "H2 ID source index: too many local leaf boxes for MPI_Allgatherv");
      }

      int comm_size = 1;
      MPI_Comm_size(H2_solver->comm, &comm_size);
      const int local_box_count = static_cast<int>(leaf.local_boxes.size());
      std::vector<int> box_counts(static_cast<size_t>(comm_size), 0);
      MPI_Allgather(
          &local_box_count, 1, MPI_INT,
          box_counts.data(), 1, MPI_INT, H2_solver->comm);

      std::vector<int> box_displacements(static_cast<size_t>(comm_size), 0);
      for (int process = 1; process < comm_size; ++process) {
          box_displacements[static_cast<size_t>(process)] =
              box_displacements[static_cast<size_t>(process - 1)] +
              box_counts[static_cast<size_t>(process - 1)];
      }
      const int total_leaf_boxes =
          box_displacements.back() + box_counts.back();
      if (static_cast<int64_t>(total_leaf_boxes) != leaf.num_boxes_global ||
          (local_box_count > 0 &&
           leaf.local_morton_start != box_displacements[static_cast<size_t>(rank)])) {
          throw std::runtime_error(
              "H2 ID source index: leaf boxes are not rank-contiguous Morton ranges");
      }

      std::vector<int64_t> local_point_counts(
          static_cast<size_t>(local_box_count), 0);
      for (int box_index = 0; box_index < local_box_count; ++box_index) {
          local_point_counts[static_cast<size_t>(box_index)] =
              leaf.local_boxes[static_cast<size_t>(box_index)].num_points;
      }
      std::vector<int64_t> global_point_counts(
          static_cast<size_t>(total_leaf_boxes), 0);
      MPI_Allgatherv(
          local_point_counts.data(), local_box_count, MPI_INT64_T,
          global_point_counts.data(), box_counts.data(),
          box_displacements.data(), MPI_INT64_T, H2_solver->comm);

      auto& leaf_offsets = H2_solver->tree->id_source_leaf_offsets;
      leaf_offsets.assign(static_cast<size_t>(total_leaf_boxes + 1), 0);
      for (int box_index = 0; box_index < total_leaf_boxes; ++box_index) {
          leaf_offsets[static_cast<size_t>(box_index + 1)] =
              leaf_offsets[static_cast<size_t>(box_index)] +
              global_point_counts[static_cast<size_t>(box_index)];
      }
      if (leaf_offsets.back() != options.N) {
          throw std::runtime_error(
              "H2 ID source index: leaf point counts do not sum to N");
      }
  }

  if (options.id_proxy_mode == 1) {
      const size_t coordinate_count =
          static_cast<size_t>(options.N) * static_cast<size_t>(options.dimension);
      H2_solver->tree->id_source_point_coords.assign(
          Locations, Locations + coordinate_count);
  }
  
  idxs++;// this is to convert from 0-based index to 1-based index for Fortran compatibility
  idxe++;

  // new2old changed to index-1 based for fortran side
  for (int& v : new2old) v += 1;


  // To Do: maybe put this section and below to c_bpack_factor
  return 0;
}


} // namespace butterfly
