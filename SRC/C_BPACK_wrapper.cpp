#include "BPACK_wrapper.h"
#ifdef HAVE_MPI
#include "butterfly_integration.hpp"
#include <mpi.h>
#endif
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <memory>
#include <cstring>
#include <getopt.h>
#include <iomanip>
#include <stdlib.h>
#include <type_traits>

using namespace std;


static inline int product(int arr[], int n) {
  int out = 1;  // Initialize product to 1
  for (int i = 0; i < n; i++) {
      out *= arr[i];
  }
  return out;
}

#ifndef HAVE_MPI
[[noreturn]] static void format7_requires_mpi(const char* caller) {
  throw std::runtime_error(
    std::string(caller) +
    ": format=7 requires a build with MPI support (enable_mpi=ON)");
}
#endif

#ifdef HAVE_MPI
// option%sym of a format-7 matrix: 1 (transpose-symmetric), or 0 (a general
// matrix) with H2_unstructured=1.
static void check_h2_symmetry_option(
    F2Cptr* option,
    MPI_Comm comm,
    const char* caller) {
  double sym_d = 0.0;
  c_bpack_getoption(option, "sym", &sym_d);
  const int sym = static_cast<int>(std::llround(sym_d));
  if (sym == 1) return;
  double unstructured_d = 0.0;
  c_bpack_getoption(option, "H2_unstructured", &unstructured_d);
  if (sym == 0 && std::llround(unstructured_d) == 1) return;

  const std::string message = std::string(caller) +
    ": ButterflyPACK H2 (format=7) takes option%sym=1 (transpose-symmetric), "
    "or option%sym=0 (nonsymmetric) with H2_unstructured=1; received option%sym=" +
    std::to_string(sym);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) std::cerr << message << std::endl;
  MPI_Abort(comm, 1);
  throw std::runtime_error(message);
}

template<typename CoordType, typename DataType>
int sync_h2_verbosity(
    F2Cptr* option,
    butterfly::H2<CoordType, DataType>* solver) {
  double verbosity_d = 0.0;
  c_bpack_getoption(option, "verbosity", &verbosity_d);
  solver->options.verbosity = static_cast<int>(std::llround(verbosity_d));
  return solver->options.verbosity;
}

// The operation a trans argument asks of format 7 ('N', 'T' or 'C', either
// case).  A transpose-symmetric H2 (option%sym=1) takes 'T' as 'N', and 'C'
// as 'N' on conjugated vectors: A^H x = conj(A conj(x)), and A^{-H} x =
// conj(A^{-1} conj(x)) (for real data 'T' again).  A nonsymmetric one
// (option%sym=0) applies A^T, or its inverse, with A's factors for 'T' (the
// transposed sweeps, h2_transposed), and 'C' as 'T' on conjugated vectors.
static char h2_trans_operation(char const* trans, const char* caller) {
  const char given = (trans != nullptr && trans[0] != '\0') ? trans[0] : 'N';
  const char op = static_cast<char>(std::toupper(static_cast<unsigned char>(given)));
  if (op != 'N' && op != 'T' && op != 'C') {
    throw std::invalid_argument(std::string(caller) +
      ": trans must be 'N', 'T' or 'C'; got '" + std::string(1, given) + "'");
  }
  return op;
}

// Whether operation `op` (h2_trans_operation) runs the transposed sweeps: 'T'
// and 'C' of a nonsymmetric H2.
template<typename DataType>
bool h2_transposed(const butterfly::H2<double, DataType>* solver, char op) {
  return op != 'N' && !solver->options.symmetric;
}

template<typename DataType>
void conjugate_h2_values(std::vector<DataType>& values) {
  if constexpr (butterfly::is_complex_v<DataType>) {
    for (DataType& value : values) value = std::conj(value);
  }
}

// `ncol` columns of the caller's local rows in the solver's (tree) order: a
// distributed64 handle converts from the caller's input ownership (`nloc`
// must be its input-local size), a legacy handle already uses the tree's.
// Returns the internal local row count.
template<typename DataType>
int h2_values_to_internal(
    const butterfly::H2<double, DataType>* solver,
    const DataType* values, int nloc, int ncol,
    std::vector<DataType>& internal, const char* caller) {
  const auto* layout = solver->distributed_layout.get();
  int internal_nloc = nloc;
  if (layout != nullptr) {
    if (layout->input_local_size != nloc) {
      throw std::invalid_argument(std::string(caller) +
        " (distributed64): the local row count must equal the input-local size");
    }
    if (layout->internal_local_size > std::numeric_limits<int>::max()) {
      throw std::overflow_error(std::string(caller) +
        " (distributed64): internal local size exceeds INT_MAX");
    }
    internal_nloc = static_cast<int>(layout->internal_local_size);
  }
  internal.resize(static_cast<size_t>(internal_nloc) * static_cast<size_t>(ncol));
  if (layout != nullptr) {
    layout->input_to_internal(values, internal.data(), ncol);
  } else {
    std::copy_n(values, internal.size(), internal.data());
  }
  return internal_nloc;
}

// ... and back to the caller's order.
template<typename DataType>
void h2_values_from_internal(
    const butterfly::H2<double, DataType>* solver,
    const std::vector<DataType>& internal, DataType* values, int ncol) {
  const auto* layout = solver->distributed_layout.get();
  if (layout != nullptr) {
    layout->internal_to_input(internal.data(), values, ncol);
  } else {
    std::copy(internal.begin(), internal.end(), values);
  }
}

// solution = A^{-1} rhs with the factorization (c_bpack_factor, precon != 2),
// both in the tree's order: the direct solve, and c_bpack_inv_mult.
// `seconds`: the sweeps (with the solve data's set-up) and the gather.
// transposed: A^{-T} rhs (a nonsymmetric H2).
template<typename DataType>
void h2_apply_factorized_inverse(
    butterfly::H2<double, DataType>* solver,
    const std::vector<DataType>& rhs, std::vector<DataType>& solution,
    int internal_nloc, int ncol, int verbosity, double* seconds = nullptr,
    bool transposed = false) {
  const double t0 = MPI_Wtime();
  std::vector<std::vector<fmm::SolveDataRequest<double, DataType>>> solve_data(
    solver->options.num_levels);
  butterfly::dispatch_h2_solve(solver, rhs, solve_data, ncol, verbosity, transposed);
  const double t1 = MPI_Wtime();
  solution.resize(static_cast<size_t>(internal_nloc) * static_cast<size_t>(ncol));
  butterfly::gather_local_solution(
    solver->tree.get(), solve_data, solution.data(), &internal_nloc, ncol);
  if (seconds != nullptr) {
    seconds[0] = t1 - t0;
    seconds[1] = MPI_Wtime() - t1;
  }
}

template<typename CoordType, typename DataType>
void compress_h2_and_update_stats(
    butterfly::H2<CoordType, DataType>* solver,
    F2Cptr* stats) {
  if (solver->build_state == butterfly::H2BuildState::H2_COMPRESSED) {
    return;
  }

  double compression_time = 0.0;
  double entryeval_time = 0.0;
  butterfly::dispatch_h2_compression(
    solver, &compression_time, &entryeval_time);

  c_bpack_setstats(stats, "Time_Fill", &compression_time);
  c_bpack_setstats(stats, "Time_Entry", &entryeval_time);

  double rank_max = static_cast<double>(solver->last_factor_rankmax);
  c_bpack_setstats(stats, "Rank_max_Constr", &rank_max);

  // bytes / 1024e3, the unit the Fortran paths report memory in.  Mem_Fill
  // is the name a caller reads for the size of a compressed matrix (they set
  // it as Mem_Comp_for + Mem_Direct_for); without it a compression-only H2
  // reported no memory at all.
  double compression_memory_mb = solver->factorization_memory / 1024.0e3;
  c_bpack_setstats(stats, "Mem_Comp_for", &compression_memory_mb);
  c_bpack_setstats(stats, "Mem_Fill", &compression_memory_mb);

  if (solver->options.verbosity >= 0) {
    (void)butterfly::dispatch_h2_compression_verification(solver);
  }
}

static butterfly::ProgramOptions read_h2_program_options64(
    F2Cptr* option, int64_t point_count, int dimension) {
  double tolerance = 0.0;
  double reduction_threshold_d = 0.0;
  double Nmin_leaf_d = 0.0;
  double precon_d = 0.0;
  double verbosity_d = 0.0;
  double H2_unstructured_d = 0.0;
  double CA_level_d = 0.0;
  double H2_use_sketch_d = 0.0;
  double H2_lazy_schur_d = 0.0;
  double H2_compress_factors_d = 0.0;
  double H2_GEMM_split_d = 0.0;
  double H2_XRR_factor_d = 0.0;
  double H2_use_gpu_d = 0.0;
  double H2_CA_staged_halo_d = 0.0;
  double H2_CA_owner_component_d = 0.0;
  double H2_CA_pack_shared_d = 0.0;
  double H2_CA_solve_d = 0.0;
  double H2_CA_owner_serial_d = 0.0;
  double H2_ID_radius_d = 0.0;
  double H2_ID_proxy_d = 0.0;
  double H2_ID_proxy_points_d = 0.0;
  double BACA_Batch_d = 0.0;
  double H2_s_op = 0.0;
  c_bpack_getoption(option, "tol_comp", &tolerance);
  c_bpack_getoption(option, "reduction_threshold", &reduction_threshold_d);
  c_bpack_getoption(option, "Nmin_leaf", &Nmin_leaf_d);
  c_bpack_getoption(option, "precon", &precon_d);
  c_bpack_getoption(option, "verbosity", &verbosity_d);
  c_bpack_getoption(option, "H2_unstructured", &H2_unstructured_d);
  c_bpack_getoption(option, "CA_level", &CA_level_d);
  c_bpack_getoption(option, "H2_use_sketch", &H2_use_sketch_d);
  c_bpack_getoption(option, "H2_lazy_schur", &H2_lazy_schur_d);
  c_bpack_getoption(option, "H2_compress_factors", &H2_compress_factors_d);
  c_bpack_getoption(option, "H2_GEMM_split", &H2_GEMM_split_d);
  c_bpack_getoption(option, "H2_XRR_factor", &H2_XRR_factor_d);
  c_bpack_getoption(option, "H2_use_gpu", &H2_use_gpu_d);
  c_bpack_getoption(option, "H2_CA_staged_halo", &H2_CA_staged_halo_d);
  c_bpack_getoption(option, "H2_CA_owner_component", &H2_CA_owner_component_d);
  c_bpack_getoption(option, "H2_CA_pack_shared", &H2_CA_pack_shared_d);
  c_bpack_getoption(option, "H2_CA_solve", &H2_CA_solve_d);
  c_bpack_getoption(option, "H2_CA_owner_serial", &H2_CA_owner_serial_d);
  c_bpack_getoption(option, "H2_ID_radius", &H2_ID_radius_d);
  c_bpack_getoption(option, "H2_ID_proxy", &H2_ID_proxy_d);
  c_bpack_getoption(option, "H2_ID_proxy_points", &H2_ID_proxy_points_d);
  c_bpack_getoption(option, "BACA_Batch", &BACA_Batch_d);
  c_bpack_getoption(option, "H2_s_op", &H2_s_op);

  const int64_t reduction_threshold =
      static_cast<int64_t>(reduction_threshold_d);
  const int64_t Nmin_leaf = static_cast<int64_t>(Nmin_leaf_d);
  const int CA_level = static_cast<int>(std::llround(CA_level_d));
  const int H2_unstructured =
      static_cast<int>(std::llround(H2_unstructured_d));
  const int H2_use_sketch =
      static_cast<int>(std::llround(H2_use_sketch_d));
  const int H2_lazy_schur =
      static_cast<int>(std::llround(H2_lazy_schur_d));
  const int H2_compress_factors =
      static_cast<int>(std::llround(H2_compress_factors_d));
  const int H2_GEMM_split =
      static_cast<int>(std::llround(H2_GEMM_split_d));
  const int H2_XRR_factor =
      static_cast<int>(std::llround(H2_XRR_factor_d));
  const int H2_use_gpu =
      static_cast<int>(std::llround(H2_use_gpu_d));
  const int H2_CA_staged_halo =
      static_cast<int>(std::llround(H2_CA_staged_halo_d));
  const int H2_CA_owner_component =
      static_cast<int>(std::llround(H2_CA_owner_component_d));
  const int H2_CA_pack_shared =
      static_cast<int>(std::llround(H2_CA_pack_shared_d));
  const int H2_CA_solve =
      static_cast<int>(std::llround(H2_CA_solve_d));
  const int H2_CA_owner_serial =
      static_cast<int>(std::llround(H2_CA_owner_serial_d));
  const int H2_ID_radius =
      static_cast<int>(std::llround(H2_ID_radius_d));
  const int H2_ID_proxy =
      static_cast<int>(std::llround(H2_ID_proxy_d));
  const int H2_ID_proxy_points =
      static_cast<int>(std::llround(H2_ID_proxy_points_d));
  const int BACA_Batch =
      static_cast<int>(std::llround(BACA_Batch_d));

  if (H2_use_sketch < 0 || H2_use_sketch > 2) {
    throw std::invalid_argument("H2_use_sketch must be 0, 1, or 2");
  }
  if (H2_lazy_schur < 0 || H2_lazy_schur > 2) {
    throw std::invalid_argument("H2_lazy_schur must be 0, 1, or 2");
  }
  if (H2_compress_factors < 0 || H2_compress_factors > 1) {
    throw std::invalid_argument("H2_compress_factors must be 0 (off) or 1 (on); 2 is not a separate mode: 1 compresses the fill sources too");
  }
  if (H2_GEMM_split < 0) {
    throw std::invalid_argument("H2_GEMM_split must be nonnegative");
  }
  if (H2_XRR_factor != 0 && H2_XRR_factor != 1) {
    throw std::invalid_argument("H2_XRR_factor must be 0 or 1");
  }
  if (H2_use_gpu < 0 || H2_use_gpu > 2) {
    throw std::invalid_argument("H2_use_gpu must be 0 (CPU) or 1 (GPU, FP64 tensor-core GEMMs); 2, the old tensor-core mode, is taken as 1");
  }
#ifndef H2_HAVE_GPU
  if (H2_use_gpu != 0) {
    throw std::invalid_argument("H2_use_gpu=1 requires a build with enable_h2_gpu=ON");
  }
#endif
  if (H2_CA_staged_halo != 0 && H2_CA_staged_halo != 2) {
    throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
  }
  if (H2_CA_owner_component != 0 && H2_CA_owner_component != 3) {
    throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
  }
  if (H2_CA_pack_shared < 0 || H2_CA_pack_shared > 2) {
    throw std::invalid_argument("H2_CA_pack_shared must be 0 (auto), 1 (always) or 2 (never)");
  }
  if (H2_CA_solve < 0 || H2_CA_solve > 1) {
    throw std::invalid_argument("H2_CA_solve must be 0 (replicated) or 1 (home)");
  }
  if (H2_CA_owner_serial != 0 && H2_CA_owner_serial != 1) {
    throw std::invalid_argument("H2_CA_owner_serial must be 0 or 1");
  }
  if (H2_CA_owner_component == 3 && H2_lazy_schur == 0) {
    throw std::invalid_argument(
        "H2_CA_owner_component=3 requires H2_lazy_schur=1 or 2");
  }
  if (H2_lazy_schur != 0 && H2_use_sketch != 2) {
    throw std::invalid_argument("H2_lazy_schur requires H2_use_sketch=2");
  }
  if (H2_ID_radius < 2) {
    throw std::invalid_argument("H2_ID_radius must be at least 2");
  }
  if (H2_ID_proxy < 0 || H2_ID_proxy > 2) {
    throw std::invalid_argument("H2_ID_proxy must be 0, 1, or 2");
  }
  if (H2_ID_proxy == 1 && H2_ID_proxy_points <= 0) {
    throw std::invalid_argument(
        "H2_ID_proxy_points must be positive when H2_ID_proxy=1");
  }
  if (H2_ID_proxy == 2 && BACA_Batch <= 0) {
    throw std::invalid_argument(
        "BACA_Batch must be positive when H2_ID_proxy=2");
  }

  butterfly::ProgramOptions result = butterfly::parse_program_options64(
      point_count, dimension, tolerance, reduction_threshold,
      Nmin_leaf, CA_level);
  result.precon = static_cast<int>(std::llround(precon_d));
  result.verbosity = static_cast<int>(std::llround(verbosity_d));
  result.unstructured = H2_unstructured;
  result.use_sketch = H2_use_sketch;
  result.lazy_schur = H2_lazy_schur;
  result.compress_factors = H2_compress_factors;
  result.gemm_split = H2_GEMM_split;
  result.xrr_factor = H2_XRR_factor;
  result.use_gpu = H2_use_gpu != 0 ? 1 : 0;
  {
    double sym_d = 1.0;
    c_bpack_getoption(option, "sym", &sym_d);
    result.symmetric = std::llround(sym_d) != 0;
  }
  result.ca_staged_halo = H2_CA_staged_halo;
  result.ca_owner_component = H2_CA_owner_component;
  result.ca_pack_shared = H2_CA_pack_shared;
  result.ca_solve = H2_CA_solve;
  result.ca_owner_serial = H2_CA_owner_serial;
  result.id_neighborhood_radius = H2_ID_radius;
  result.id_proxy_mode = H2_ID_proxy;
  if (H2_ID_proxy == 1) {
    result.id_proxy_points = H2_ID_proxy_points;
  } else if (H2_ID_proxy == 2) {
    result.id_adaptive_batch = BACA_Batch;
  }
  result.s_op = H2_s_op;
  if (result.precon == 2) result.CA_level = result.num_levels;
  butterfly::validate_h2_backend_selection(result);
  return result;
}

template<typename H2Data>
static butterfly::H2<double, H2Data>* get_h2_solver(F2Cptr* bmat,
                                                    const char* caller) {
  if (bmat == nullptr || *bmat == nullptr) {
    throw std::invalid_argument(std::string(caller) + ": null matrix handle");
  }
  void* raw = nullptr;
  c_bpack_get_h2(*bmat, &raw);
  if (raw == nullptr) {
    throw std::runtime_error(
        std::string(caller) +
        ": currently implemented only for format=7 (H2)");
  }
  return static_cast<butterfly::H2<double, H2Data>*>(raw);
}

template<typename H2Data>
static butterfly::DistributedLayout64* get_distributed_layout(
    F2Cptr* bmat, const char* caller) {
  auto* solver = get_h2_solver<H2Data>(bmat, caller);
  if (!solver->distributed_layout) {
    throw std::runtime_error(
        std::string(caller) +
        ": matrix was not initialized with the distributed64 API");
  }
  return solver->distributed_layout.get();
}
#endif

// The command line parser for the example related parameters
void c_bpack_set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0) {

	struct OptionHelp {
		const char* name;
		const char* description;
	};

	static const OptionHelp option_help_table[] = {
		{"nmin_leaf",       "leafsize in the hierarchical partitioning"},
		{"tol_comp",        "relative tolerance for matrix construction"},
		{"tol_rand",        "relative tolerance for matrix inversion"},
		{"tol_Rdetect",     "relative tolerance for rank detection during matrix inversion"},
		{"tol_itersol",     "convergence tolerance for the iterative solver (format-7 precon=2 uses H2 BiCGSTAB)"},
		{"n_iter",          "maximum iteration count for the iterative solver"},
		{"IR_HODLR", "maximum HODLR direct-solve refinement steps; 0 disables refinement"},
		{"level_check",     "the level in the hierarchical partitioning where the randomized construction algorithm is tested, set to 10000 by default (no checking)"},
		{"precon",          "the use mode of butterflypack: 1: as a direct solver 2: as an iterative solver (compress the matrix and use it without a preconditioner), 3: as a preconditioned iterative solver (compress and invert the matrix, using the approximate inverse as a preconditioner)"},
		{"xyzsort",         "the hierarchical partitioning algorithm: 0: no permutation 1: permutation based on KD-tree 2: permutation based on cobble-like partitioning"},
		{"lrlevel",         "the level in the hierarchical partitioning (top-down numbered) above which butterfly is used and below which low-rank is used"},
		{"sym",             "matrix symmetry flag; format-7 H2 takes sym=1 (transpose-symmetric) or sym=0 with h2_unstructured=1, and sym=1 selects symmetric HODLR when format=1"},
		{"reduction_threshold", "format-7 H2 boxes-per-process threshold for reducing the active MPI process count"},
		{"h2_unstructured", "format-7 H2 backend: 0 structured Color/CA, 1 color_unstructured (Color only)"},
		{"h2_use_sketch",   "format-7 H2 ID mode: 0 full workspace, 1 materialized sparse sketch, 2 streamed sparse sketch on color levels"},
		{"h2_lazy_schur",   "format-7 H2 color Schur mode: 0 eager, 1 lazy far, 2 lazy far plus generated near"},
		{"h2_gemm_split",   "maximum task split for one format-7 H2 color work item; 0 disables splitting"},
		{"h2_xrr_factor",   "format-7 H2 X_RR pivot factorization: 0 Bunch-Kaufman, 1 LU with partial pivoting"},
		{"h2_use_gpu",      "format-7 H2 Color levels: 0 CPU, 1 GPU backend (FP64 tensor-core GEMMs for real data; 2 is taken as 1) (build with enable_h2_gpu)"},
		{"h2_compress_factors", "format-7 H2 GPU: lossy compression (nvCOMP bitcomp) at tol_comp: 0 off, 1 on - the factors kept on the device for the solve (X_NR, T, X_SR) and the factorization's own data: the fill sources later waves and the level transition read, the near and parent blocks held under memory pressure, and the CA blocks on more than one rank (so the factorization is lossy too, within the tolerance)"},
		{"hodlr_use_gpu",   "format-1 HODLR with lrlevel=0: 0 CPU, 1 GPU construction/multiply/factorization/solve, 2 the same with FP64 tensor-core GEMMs (build with enable_h2_gpu)"},
		{"hodlr_gpu_pieces", "HODLR GPU construction: pieces per rank of a low-rank block shared by several ranks, spread over its ranks by estimated cost (rounded down to a power of 2; 1: the CPU's pieces; default 4)"},
		{"jitter",          "relative size of the pivots and diagonal shifts that keep dense and small factorizations nonsingular"},
		{"h2_ca_staged_halo", "format-7 H2 CA halo mode: 0 legacy, 2 staged/overlapped"},
		{"h2_ca_owner_component", "format-7 H2 CA ownership mode: 0 replicated, 3 asynchronous components"},
		{"h2_ca_solve", "format-7 H2 GPU, replicated CA levels: how the level is solved: replicated (0: the ghosts' factors kept and applied again on every rank holding them) or home (1: each box's factors applied by its home rank only, the vector exchanged by wave)"},
		{"h2_ca_pack_shared", "format-7 H2 GPU, replicated CA levels: near blocks held by more than one rank packed through the ghosts' groups: auto (0, by the level's estimate), always (1), never (2: the card must hold them raw)"},
		{"h2_ca_owner_serial", "serialize format-7 H2 CA component-owner events (0 or 1)"},
		{"h2_id_radius",    "format-7 H2 mandatory ID neighborhood radius; 2 keeps the standard workspace"},
		{"h2_id_proxy",     "format-7 H2 proxy mode: 0 none, 1 geometric surface, 2 adaptive row sampling"},
		{"h2_id_proxy_points", "geometric surface samples when h2_id_proxy=1"},
		{"h2_s_op",         "format-7 H2 IDs: an absolute floor tol_comp*h2_s_op on the rank rule, |R_kk| > max(tol_comp |R_11|, tol_comp h2_s_op); 0 off (default). At most ||A|| for compression-only H2 (precon 2), at most 1/||A^-1|| for the factorization (K + tau^2 I: tau^2)"},
		{"errfillfull",     "errfillfull: a slow (n^2), thorough error checking is performed after the compression of each block"},
		{"baca_batch",      "block size in batched ACA when reclr_leaf=4 or 5; adaptive H2 rows per spatial node when h2_id_proxy=2"},
		{"reclr_leaf",      "low-rank compression algorithms 1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton 7: ACA with naive parallelization"},
		{"nogeo",           "whether there is geometry information provided 1: is no geometry (xyzsort can not be 1 or 2), 0: there is geometry"},
		{"less_adapt",      "1: improved randomized butterfly construction, default to 1"},
		{"errsol",          "1: generating an artificial true solution vector, compute the RHS with compressed matrix, solve the system, and compare with true solution vector"},
		{"lr_blk_num",      "sqrt of #of subblocks in H-BACA, default to 1"},
		{"rank0",           "initial rank guess in the randomized butterfly algorithm, default to 32"},
		{"rankrate",        "increasing ratio of the rank guess in each iteration, default to 2"},
		{"itermax",         "maximum number of iterations in the randomized butterfly algorithm, default to 10"},
		{"powiter",         "order of power iteration in the randomized low-rank construction"},
		{"ilu",             "whether symmetric gauss-seidel is used when format=2"},
		{"nbundle",         "multiply nbundle sets of vectors together in randomized butterfly algorithm for better flop performance, default to 1"},
		{"near_para",       "admissibility parameter when format=2/3/4/5, strong admissibility typically requires near_para>2.0"},
		{"format",          "the hierarchical matrix format: 1: HODLR/HODBF 2: H matrix 3: HSSBF/SHNBF 4: HSSBF_MD/SHNBF_MD 5: block-LR/BF"},
		{"verbosity",       "-1 suppresses output, 0 prints summaries, and 1 or greater prints details"},
		{"rmax",            "preestimate of the maximum rank for allocating buffers, default to 1000"},
		{"sample_para",     "oversampling factor in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
		{"pat_comp",        "pattern of entry-evaluation-based butterfly compression: 1 from right to left, 2 from left to right, 3 from outer to inner"},
		{"knn",             "nearest neighbouring points used in improved BACA and entry-evaluation-based butterfly compression"},
		{"knn_near_para",   "admissibility parameter for guiding the nearest neighbouring points search"},
		{"forwardN15flag",  "whether to use N15 or NlogN algorithm for entry-evaluation-based matrix butterfly compression"},
		{"sample_para_outer","oversampling factor for the outtermost factor matrices in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
		{"elem_extract",    "0: evaluating entries one by one 1: evaluating entries block by block (may requires communication inside the callback function) 2: evaluating entries block by block (no communication allowed inside the callback function)"},
		{"fastsample_tensor","0: uniformly sample each dimension. 1: uniformly sample the rows of the unfolded matrices on top of 0. 2: use translation invariance"},
		{"trans_invariant", "1: reuse HTENSOR blocks for translational-invariant tensor kernels on uniform Cartesian grids"},
		{"htensor_mvp_level_batch", "number of HTENSOR levels grouped per MVP call; 1 keeps level-by-level memory"},
		{"use_fft_circulant","whether to use FFT-based circulant representative blocks"},
		{"fftw_plan_mode",  "FFTW apply-plan mode: 0 estimate, 1 measure, 2 patient, 3 exhaustive"},
		{"use_zfp",         "whether to use zfp compression"},
		{"use_qtt",         "whether to use qtt compression"},
		{"hextralevel",         "HMAT: extra levels for top partitioning of the H matrix based on MPI counts. BLR: Maxlevel-hextralevel is the level for defining B-LR/B-BF blocks"},
		{"iter_solver",         "The choice of iterative solvers. 1: TFQMR, 2: GMRES, 3: IterativeRefinement, or 4: CG"},
		{"help",            "print this help message"}
	};

	double opt_d;
	int opt_i;
	std::vector<std::unique_ptr<char[]>> argv_data(argc);
	std::vector<char*> argv(argc);
	for (int i=0; i<argc; i++) {
		argv_data[i].reset(new char[strlen(cargv[i])+1]);
		argv[i] = argv_data[i].get();
		strcpy(argv[i], cargv[i]);
	}
	option long_options[] =
		{
		{"help",             no_argument,       0, 1000},
		{"nmin_leaf",                     required_argument, 0, 1},
		{"tol_comp",                   required_argument, 0, 2},
		{"tol_rand",                   required_argument, 0, 3},
		{"tol_Rdetect",             required_argument, 0, 4},
		{"tol_itersol",             required_argument, 0, 5},
		{"n_iter",          required_argument, 0, 6},
		{"level_check",         required_argument, 0, 7},
		{"precon",                  required_argument, 0, 8},
		{"xyzsort",                  required_argument, 0, 9},
		{"lrlevel",     required_argument, 0, 10},
		{"errfillfull",       required_argument, 0, 11},
		{"baca_batch",      required_argument, 0, 12},
		{"reclr_leaf",      required_argument, 0, 13},
		{"nogeo",     required_argument, 0, 14},
		{"less_adapt",            required_argument, 0, 15},
		{"errsol",           required_argument, 0, 16},
		{"lr_blk_num",                  required_argument, 0, 17},
		{"rank0",  required_argument, 0, 18},
		{"rankrate", required_argument, 0, 19},
		{"itermax",               required_argument, 0, 20},
		{"powiter",  required_argument, 0, 21},
		{"ilu", required_argument, 0, 22},
		{"nbundle",     required_argument, 0, 23},
		{"near_para",  required_argument, 0, 24},
		{"format",  required_argument, 0, 25},
		{"verbosity", required_argument, 0, 26},
		{"rmax", required_argument, 0, 27},
		{"sample_para", required_argument, 0, 28},
		{"pat_comp",    required_argument, 0, 29},
		{"knn",         required_argument, 0, 30},
		{"knn_near_para",         required_argument, 0, 31},
		{"forwardN15flag",         required_argument, 0, 32},
		{"sample_para_outer",         required_argument, 0, 33},
		{"elem_extract",         required_argument, 0, 34},
		{"fastsample_tensor",         required_argument, 0, 35},
		{"use_zfp",         required_argument, 0, 36},
		{"use_qtt",         required_argument, 0, 37},
		{"hextralevel",         required_argument, 0, 38},
		{"iter_solver",         required_argument, 0, 39},
		{"trans_invariant",         required_argument, 0, 40},
			{"htensor_mvp_level_batch", required_argument, 0, 41},
		{"use_fft_circulant",         required_argument, 0, 42},
		{"fftw_plan_mode",         required_argument, 0, 43},
		{"fft_plan_mode",         required_argument, 0, 44},
		{"sym",         required_argument, 0, 45},
		{"IR_HODLR", required_argument, 0, 46},
		{"h2_use_sketch", required_argument, 0, 47},
		{"H2_use_sketch", required_argument, 0, 47},
		{"h2_id_radius", required_argument, 0, 48},
		{"H2_ID_radius", required_argument, 0, 48},
		{"h2_id_proxy", required_argument, 0, 49},
		{"H2_ID_proxy", required_argument, 0, 49},
		{"h2_id_proxy_points", required_argument, 0, 50},
		{"H2_ID_proxy_points", required_argument, 0, 50},
		{"reduction_threshold", required_argument, 0, 51},
		{"h2_lazy_schur", required_argument, 0, 52},
		{"H2_lazy_schur", required_argument, 0, 52},
		{"h2_gemm_split", required_argument, 0, 53},
		{"H2_GEMM_split", required_argument, 0, 53},
		{"h2_ca_staged_halo", required_argument, 0, 54},
		{"H2_CA_staged_halo", required_argument, 0, 54},
		{"h2_ca_owner_component", required_argument, 0, 55},
		{"H2_CA_owner_component", required_argument, 0, 55},
		{"h2_ca_pack_shared", required_argument, 0, 1002},
		{"H2_CA_pack_shared", required_argument, 0, 1002},
		{"h2_ca_solve", required_argument, 0, 1003},
		{"H2_CA_solve", required_argument, 0, 1003},
		{"h2_ca_owner_serial", required_argument, 0, 56},
		{"H2_CA_owner_serial", required_argument, 0, 56},
		{"h2_unstructured", required_argument, 0, 57},
		{"H2_unstructured", required_argument, 0, 57},
		{"h2_xrr_factor", required_argument, 0, 58},
		{"H2_XRR_factor", required_argument, 0, 58},
		{"h2_use_gpu", required_argument, 0, 59},
		// (new options take codes from 1000 up: getopt returns '?' (63) for
		// an option it does not know, and a driver passes its own options
		// through here too)
		{"h2_compress_factors", required_argument, 0, 1001},
		{"H2_compress_factors", required_argument, 0, 1001},
		{"h2_s_op", required_argument, 0, 1004},
		{"H2_s_op", required_argument, 0, 1004},
		{"H2_use_gpu", required_argument, 0, 59},
		{"hodlr_use_gpu", required_argument, 0, 60},
		{"HODLR_use_gpu", required_argument, 0, 60},
		{"hodlr_gpu_pieces", required_argument, 0, 62},
		{"HODLR_gpu_pieces", required_argument, 0, 62},
		{"jitter", required_argument, 0, 61},
		{NULL, 0, NULL, 0}
		};
	int c, option_index = 0;
	// bool unrecognized_options = false;
	opterr = optind = 0;
	while ((c = getopt_long_only
			(argc, argv.data(), "",
			long_options, &option_index)) != -1) {
		switch (c) {
		case 1000: {
		std::cout << "Available ButterflyPACK Command-Line Options:\n";
		for (const auto& opt : option_help_table) {
			std::cout << "  --" << std::setw(20) << std::left << opt.name
						<< " : " << opt.description << "\n";
		}
		} break;
		case 1: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
		} break;
		case 2: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_comp", opt_d);
		c_bpack_set_D_option(&option0, "tol_rand", opt_d);
		c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
		} break;
		case 3: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_rand", opt_d);
		} break;
		case 4: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
		} break;
		case 5: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
		} break;
		case 6: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "n_iter", opt_i);
		} break;
		case 7: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "level_check", opt_i);
		} break;
		case 8: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "precon", opt_i);
		} break;
		case 9: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "xyzsort", opt_i);
		} break;
		case 10: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "LRlevel", opt_i);
		} break;
		case 11: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
		} break;
		case 12: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
		} break;
		case 13: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
		} break;
		case 14: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "nogeo", opt_i);
		} break;
		case 15: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "less_adapt", opt_i);
		} break;
		case 16: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ErrSol", opt_i);
		} break;
		case 17: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
		} break;
		case 18: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "rank0", opt_i);
		} break;
		case 19: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "rankrate", opt_d);
		} break;
		case 20: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "itermax", opt_i);
		} break;
		case 21: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "powiter", opt_i);
		} break;
		case 22: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ILU", opt_i);
		} break;
		case 23: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "Nbundle", opt_i);
		} break;
		case 24: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "near_para", opt_d);
		} break;
		case 25: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "format", opt_i);
		} break;
		case 26: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "verbosity", opt_i);
		} break;
		case 27: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "rmax", opt_i);
		} break;
		case 28: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "sample_para", opt_d);
		} break;
		case 29: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "pat_comp", opt_i);
		} break;
		case 30: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "knn", opt_i);
		} break;
		case 31: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "knn_near_para", opt_d);
		} break;
		case 32: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "forwardN15flag", opt_i);
		} break;
		case 33: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "sample_para_outer", opt_d);
		} break;
		case 34: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "elem_extract", opt_i);
		} break;
		case 35: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "fastsample_tensor", opt_i);
		} break;
		case 40: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "trans_invariant", opt_i);
		} break;
			case 41: {
			std::istringstream iss(optarg);
			iss >> opt_i;
			c_bpack_set_I_option(&option0, "htensor_mvp_level_batch", opt_i);
			} break;
		case 42: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "use_fft_circulant", opt_i);
		} break;
		case 43:
		case 44: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "fftw_plan_mode", opt_i);
		} break;
		case 45: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "sym", opt_i);
		} break;
		case 46: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "IR_HODLR", opt_i);
		} break;
		case 47: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_use_sketch", opt_i);
		} break;
		case 48: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_ID_radius", opt_i);
		} break;
		case 49: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_ID_proxy", opt_i);
		} break;
		case 50: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_ID_proxy_points", opt_i);
		} break;
		case 51: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "reduction_threshold", opt_i);
		} break;
		case 52: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_lazy_schur", opt_i);
		} break;
		case 53: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_GEMM_split", opt_i);
		} break;
		case 54: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_CA_staged_halo", opt_i);
		} break;
		case 55: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_CA_owner_component", opt_i);
		} break;
		case 56: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_CA_owner_serial", opt_i);
		} break;
		case 57: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_unstructured", opt_i);
		} break;
		case 58: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_XRR_factor", opt_i);
		} break;
		case 59: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_use_gpu", opt_i);
		} break;
		case 60: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "HODLR_use_gpu", opt_i);
		} break;
		case 62: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "HODLR_gpu_pieces", opt_i);
		} break;
		case 1001: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "H2_compress_factors", opt_i);
		} break;
		case 1002: {
		const std::string v(optarg);
		if (v == "auto") opt_i = 0;
		else if (v == "always") opt_i = 1;
		else if (v == "never") opt_i = 2;
		else { std::istringstream iss(v); iss >> opt_i; }
		c_bpack_set_I_option(&option0, "H2_CA_pack_shared", opt_i);
		} break;
		case 1003: {
		const std::string v(optarg);
		if (v == "replicated") opt_i = 0;
		else if (v == "home") opt_i = 1;
		else { std::istringstream iss(v); iss >> opt_i; }
		c_bpack_set_I_option(&option0, "H2_CA_solve", opt_i);
		} break;
		case 1004: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "H2_s_op", opt_d);
		} break;
		case 61: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "jitter", opt_d);
		} break;
		case 36: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "use_zfp", opt_i);
		} break;
		case 37: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "use_qtt", opt_i);
		} break;
		case 38: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "hextralevel", opt_i);
		} break;
		case 39: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "iter_solver", opt_i);
		} break;
		default: break;
		}
	}
	}





void c_bpack_vector_global2local(F2Cptr* ptree, F2Cptr* msh, int* nvec_p, C_DT   *b_global, C_DT   *b)
{
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;
	int nvec = *nvec_p;

    int iam, size;
    int master_rank = 0;


#ifdef HAVE_MPI

	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	fst_row = idxs-1;

	{
	C_DT* b_global_permed=NULL;
	C_DT* btmp= (C_DT*)malloc(m_loc * nvec*sizeof(C_DT));
	int *counts = NULL, *displs = NULL;

	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		b_global_permed= (C_DT*)malloc(m * nvec*sizeof(C_DT));
		for (int j = 0; j < nvec; ++j) {
			for (int i = 0; i < m; ++i) {
				int i_old= i+1;
				int i_new;
				c_bpack_old2new(msh,&i_old,&i_new);
				b_global_permed[j + (i_new-1) * nvec] = b_global[j * m + i];
			}
		}
	}
	m_loc1 = m_loc*nvec;
	fst_row1 = fst_row*nvec;
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, comm);
	MPI_Scatterv(b_global_permed, counts, displs, C_MPI_DT,  btmp, m_loc * nvec, C_MPI_DT, 0, comm);

	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_loc; ++i)
		{
			b[j * m_loc + i] = btmp[j + i * nvec];
		}
	}
	if (iam == 0) {
		free(counts);
		free(displs);
		free(b_global_permed);
	}
	free(btmp);
	}
#else
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			b[j * m + i] = b_global[j * m + i];
		}
	}
#endif

}



void c_bpack_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* nvec_p, C_DT   *b, C_DT   *b_global)
{
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;
	int nvec = *nvec_p;
    int iam, size;
    int master_rank = 0;


#ifdef HAVE_MPI
	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	fst_row = idxs-1;

	{
	C_DT* b_global_permed;
	C_DT* btmp= (C_DT*)malloc(m_loc * nvec*sizeof(C_DT));
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_loc; ++i)
		{
			btmp[j + i * nvec] = b[j * m_loc + i];
		}
	}
	int *counts = NULL, *displs = NULL;
	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		b_global_permed= (C_DT*)malloc(m * nvec*sizeof(C_DT));
	}
	m_loc1 = m_loc*nvec;
	fst_row1 = fst_row*nvec;
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, comm);
	MPI_Gatherv(btmp, m_loc * nvec, C_MPI_DT, b_global_permed, counts, displs, C_MPI_DT, 0, comm);
	if (iam == 0) {
		for (int j = 0; j < nvec; ++j) {
			for (int i = 0; i < m; ++i) {
				int i_old= i+1;
				int i_new;
				c_bpack_old2new(msh,&i_old,&i_new);
				b_global[j * m + i] = b_global_permed[j + (i_new-1) * nvec];
			}
		}
		free(counts);
		free(displs);
		free(b_global_permed);
	}
	free(btmp);
	}
#else

	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			b_global[j * m + i] = b[j * m + i];
		}
	}

#endif

}


void c_bpack_md_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* ndim_p, int* nvec_p, C_DT   *b, C_DT   *b_global)
{

	int nvec = *nvec_p;
	int Ndim= *ndim_p;
    int iam, size;
    int master_rank = 0;


	int *idxs = (int *) malloc(Ndim*sizeof(int));
	int *m_loc = (int *) malloc(Ndim*sizeof(int));
	int *m = (int *) malloc(Ndim*sizeof(int));

#ifdef HAVE_MPI
	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_md_localindices(&Ndim, msh, idxs, m_loc, m);
	size_t m_scalar = product(m,Ndim);
	int m_loc1 = product(m_loc,Ndim)*nvec;

	int *counts = NULL, *displs = NULL, *m_locs=NULL, *idxss=NULL;
	C_DT *b_global_tmp = NULL;
	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		m_locs = (int *) malloc(Ndim*size*sizeof(int));
		idxss = (int *) malloc(Ndim*size*sizeof(int));
		b_global_tmp= (C_DT*)malloc(m_scalar * nvec*sizeof(C_DT));
	}

	MPI_Gather(m_loc, Ndim, MPI_INT, m_locs, Ndim, MPI_INT, 0, comm);
	MPI_Gather(idxs, Ndim, MPI_INT, idxss, Ndim, MPI_INT, 0, comm);
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	if (iam == 0) {
		displs[0] = 0;
		for (int i = 1; i < size; i++) {
			displs[i] = displs[i-1] + counts[i-1];
		}
	}
	MPI_Gatherv(b, m_loc1, C_MPI_DT, b_global, counts, displs, C_MPI_DT, 0, comm);

	if (iam == 0) {
		//Converting from the "stacked" order (rank 0's data, then rank 1's, etc.) to the permuted global ordering
		for (int s=0;s<size;s++){
			int m_loc_1vec = product(&(m_locs[s*Ndim]),Ndim);
			for (int i=0; i<m_loc_1vec; i++){
				int i_new_scalar = i+1;
				int i_new_md[Ndim];
				c_bpack_singleindex_to_multiindex(&Ndim,&(m_locs[s*Ndim]),&i_new_scalar,i_new_md);
				for (int j=0;j<Ndim;j++){
					i_new_md[j]=i_new_md[j]+idxss[s*Ndim+j]-1;
				}
				c_bpack_multiindex_to_singleindex(&Ndim,m,&i_new_scalar,i_new_md);
				for (int nth=0; nth<nvec; nth++){
					b_global_tmp[(i_new_scalar-1)+m_scalar*nth] = b_global[displs[s]+i+m_loc_1vec*nth];
				}
			}
		}

		//Converting from permuted global ordering to the original ordering
		for (int i = 0; i < m_scalar; ++i) {
			int i_old_scalar=i+1;
			int i_new_scalar;
			int i_old_md[Ndim];
			int i_new_md[Ndim];
			c_bpack_singleindex_to_multiindex(&Ndim,m,&i_old_scalar,i_old_md);
			c_bpack_md_old2new(&Ndim, msh, i_old_md, i_new_md);
			c_bpack_multiindex_to_singleindex(&Ndim,m,&i_new_scalar,i_new_md);
			for (int nth=0; nth<nvec; nth++){
				b_global[nth * m_scalar + i] = b_global_tmp[nth * m_scalar + (i_new_scalar-1)];
			}
		}

		free(counts);
		free(displs);
		free(m_locs);
		free(idxss);
		free(b_global_tmp);
	}
#else

	c_bpack_md_localindices(&Ndim, msh, idxs, m_loc, m);
	size_t m_scalar = product(m,Ndim);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_scalar; ++i)
		{
			b_global[j * m_scalar + i] = b[j * m_scalar + i];
		}
	}
#endif

free(idxs);
free(m_loc);
free(m);

}

void c_bpack_construct_init(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms,
	int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree,
	void (*C_FuncDistmn)(int*, int*, double*,C2Fptr),
	void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp){
	// To Do: need to compute msh (idxs, idxe, new2old), Npo_locs, perms
    // C_FuncDistmn: defines distance
	// C_FuncNearFar:

  // correspond to create_uniform_tree
  // arguments in create_uniform_tree:
  //   point_coords: null_ptr (assign uniform grid) or a pointer to an array of points on the grid
  //   num_points: N, the dimension of the matrix
  //   num_levels: number of levels for factorization
  //   global_bounds: bounding min and max of all dimensions
  //   dimension: 2D, or 3D problem
  //   comm: some form of MPI info, gets mpi_rank, mpi_size
  //   reduction_threshold: only uniform reduction pattern is supported, so really there is only one option which is uniform
  //   pattern: only uniform reduction pattern is supported, so really there is only one option which is uniform
  //   returns: tree structure

  //   can also create HierarchicalFactorization object
  //   * @param N Total number of points in the problem
  //   * @param prop Matrix property (symmetric, hermitian, or nonsymmetric)
  //   * @param kernel_func Kernel evaluator
  //   * @param dim Spatial dimension (2 or 3)
  //   * @param factorization_type Method for factorizing/inverting matrices (default: CHOLESKY)
  //   * @param num_proxy Number of proxy points (-1 uses default 32 for 2D, 256 for 3D)
  //   * @param tol Compression tolerance (default: 1e-6)
  //   * @param proxy_factor Proxy surface radius factor (default: 2.5)

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

  // Butterflypack end: need some definition of proxy points

  //   Npo: pass into num_points
  //   Ndim: pass into dimension
  //   Locations: pass into point_coords
  //   nns:
  //   nlevel: pass into num_levels
  //   tree: type difference with tree returned by create_uniform_tree
  //   perms: permutation vector?
  //   Npo_loc:
  //   bmat: this stores the h2 solver struct, h2 tree, kernel, etc.
  //   option:
  //   stats:
  //   msh
  //   ker: kernel types from FMM?
  //   ptree: mpi communicator needed, otherwise not relevant
  //   C_FuncDistmn
  //   C_FuncNearFar
  //   C_QuantApp


  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	MPI_Comm mpi_comm = MPI_Comm_f2c((MPI_Fint)fcomm);
	check_h2_symmetry_option(option, mpi_comm, "c_bpack_construct_init");

	// use datatype C_DT
    // construct H2 solver
	// auto H2_solver = std::make_unique<H2<double, C_DT>>();
	using H2Data = typename butterfly::fmm_data<C_DT>::type;
    butterfly::H2<double, H2Data>* H2_solver = new butterfly::H2<double, H2Data>();

    H2_solver->comm = mpi_comm;
    //*bmat = static_cast<F2Cptr>(H2_solver.release());

	int rank = 0;
    int size = 1;
	MPI_Comm_rank(H2_solver->comm, &rank);
    MPI_Comm_size(H2_solver->comm, &size);

    // // parse options and throw error if any requirements are not satisfied
    // for (int i = 1; i < argc; ++i) {
    //   const std::string arg = argv[i];
    //   if (arg == "--help" || arg == "-h") {
    //     if (rank == 0) {
    //       print_usage(argv[0]);
    //     }
    //     MPI_Finalize();
    //     return 0;
    //   }
    // }

    butterfly::ProgramOptions H2_options;
    try {
      double tolerance;
      double reduction_threshold_d;
      double Nmin_leaf_d;
      double precon_d;
      double verbosity_d;
      double H2_unstructured_d;
      double CA_level_d;
      double H2_use_sketch_d;
      double H2_lazy_schur_d;
      double H2_compress_factors_d;
      double H2_GEMM_split_d;
      double H2_XRR_factor_d;
      double H2_use_gpu_d;
      double H2_CA_staged_halo_d;
      double H2_CA_owner_component_d;
      double H2_CA_pack_shared_d;
      double H2_CA_solve_d;
      double H2_CA_owner_serial_d;
      double H2_ID_radius_d;
      double H2_ID_proxy_d;
      double H2_ID_proxy_points_d;
      double BACA_Batch_d;
      double H2_s_op = 0.0;
      c_bpack_getoption(option, "tol_comp", &tolerance);
      c_bpack_getoption(option, "reduction_threshold", &reduction_threshold_d);
      c_bpack_getoption(option, "Nmin_leaf", &Nmin_leaf_d);
      c_bpack_getoption(option, "precon", &precon_d);
      c_bpack_getoption(option, "verbosity", &verbosity_d);
      c_bpack_getoption(option, "H2_unstructured", &H2_unstructured_d);
      c_bpack_getoption(option, "CA_level", &CA_level_d);
      c_bpack_getoption(option, "H2_use_sketch", &H2_use_sketch_d);
      c_bpack_getoption(option, "H2_lazy_schur", &H2_lazy_schur_d);
      c_bpack_getoption(option, "H2_compress_factors", &H2_compress_factors_d);
      c_bpack_getoption(option, "H2_GEMM_split", &H2_GEMM_split_d);
      c_bpack_getoption(option, "H2_XRR_factor", &H2_XRR_factor_d);
      c_bpack_getoption(option, "H2_use_gpu", &H2_use_gpu_d);
      c_bpack_getoption(option, "H2_CA_staged_halo", &H2_CA_staged_halo_d);
      c_bpack_getoption(option, "H2_CA_owner_component", &H2_CA_owner_component_d);
      c_bpack_getoption(option, "H2_CA_pack_shared", &H2_CA_pack_shared_d);
      c_bpack_getoption(option, "H2_CA_solve", &H2_CA_solve_d);
      c_bpack_getoption(option, "H2_CA_owner_serial", &H2_CA_owner_serial_d);
      c_bpack_getoption(option, "H2_ID_radius", &H2_ID_radius_d);
      c_bpack_getoption(option, "H2_ID_proxy", &H2_ID_proxy_d);
      c_bpack_getoption(option, "H2_ID_proxy_points", &H2_ID_proxy_points_d);
      c_bpack_getoption(option, "BACA_Batch", &BACA_Batch_d);
      c_bpack_getoption(option, "H2_s_op", &H2_s_op);
      int64_t reduction_threshold = (int64_t)reduction_threshold_d;
      int64_t Nmin_leaf = (int64_t)Nmin_leaf_d;
      const int CA_level = static_cast<int>(std::llround(CA_level_d));
      const int H2_unstructured =
          static_cast<int>(std::llround(H2_unstructured_d));
      const int H2_use_sketch = static_cast<int>(std::llround(H2_use_sketch_d));
      const int H2_lazy_schur = static_cast<int>(std::llround(H2_lazy_schur_d));
      const int H2_compress_factors = static_cast<int>(std::llround(H2_compress_factors_d));
      const int H2_GEMM_split = static_cast<int>(std::llround(H2_GEMM_split_d));
      const int H2_XRR_factor = static_cast<int>(std::llround(H2_XRR_factor_d));
      const int H2_use_gpu = static_cast<int>(std::llround(H2_use_gpu_d));
      const int H2_CA_staged_halo =
          static_cast<int>(std::llround(H2_CA_staged_halo_d));
      const int H2_CA_owner_component =
          static_cast<int>(std::llround(H2_CA_owner_component_d));
      const int H2_CA_pack_shared =
          static_cast<int>(std::llround(H2_CA_pack_shared_d));
      const int H2_CA_solve =
          static_cast<int>(std::llround(H2_CA_solve_d));
      const int H2_CA_owner_serial =
          static_cast<int>(std::llround(H2_CA_owner_serial_d));
      const int H2_ID_radius = static_cast<int>(std::llround(H2_ID_radius_d));
      const int H2_ID_proxy = static_cast<int>(std::llround(H2_ID_proxy_d));
      const int H2_ID_proxy_points =
          static_cast<int>(std::llround(H2_ID_proxy_points_d));
      const int BACA_Batch = static_cast<int>(std::llround(BACA_Batch_d));
      if (H2_use_sketch < 0 || H2_use_sketch > 2) {
        throw std::invalid_argument("H2_use_sketch must be 0, 1, or 2");
      }
      if (H2_lazy_schur < 0 || H2_lazy_schur > 2) {
        throw std::invalid_argument("H2_lazy_schur must be 0, 1, or 2");
      }
      if (H2_compress_factors < 0 || H2_compress_factors > 1) {
        throw std::invalid_argument("H2_compress_factors must be 0 (off) or 1 (on); 2 is not a separate mode: 1 compresses the fill sources too");
      }
      if (H2_GEMM_split < 0) {
        throw std::invalid_argument("H2_GEMM_split must be nonnegative");
      }
      if (H2_XRR_factor != 0 && H2_XRR_factor != 1) {
        throw std::invalid_argument("H2_XRR_factor must be 0 or 1");
      }
      if (H2_use_gpu < 0 || H2_use_gpu > 2) {
        throw std::invalid_argument("H2_use_gpu must be 0 (CPU) or 1 (GPU, FP64 tensor-core GEMMs); 2, the old tensor-core mode, is taken as 1");
      }
#ifndef H2_HAVE_GPU
      if (H2_use_gpu != 0) {
        throw std::invalid_argument("H2_use_gpu=1 requires a build with enable_h2_gpu=ON");
      }
#endif
      if (H2_CA_staged_halo != 0 && H2_CA_staged_halo != 2) {
        throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
      }
      if (H2_CA_owner_component != 0 && H2_CA_owner_component != 3) {
        throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
      }
      if (H2_CA_pack_shared < 0 || H2_CA_pack_shared > 2) {
        throw std::invalid_argument("H2_CA_pack_shared must be 0 (auto), 1 (always) or 2 (never)");
      }
      if (H2_CA_solve < 0 || H2_CA_solve > 1) {
        throw std::invalid_argument("H2_CA_solve must be 0 (replicated) or 1 (home)");
      }
      if (H2_CA_owner_serial != 0 && H2_CA_owner_serial != 1) {
        throw std::invalid_argument("H2_CA_owner_serial must be 0 or 1");
      }
      if (H2_CA_owner_component == 3 && H2_lazy_schur == 0) {
        throw std::invalid_argument(
            "H2_CA_owner_component=3 requires H2_lazy_schur=1 or 2");
      }
      if (H2_lazy_schur != 0 && H2_use_sketch != 2) {
        throw std::invalid_argument(
            "H2_lazy_schur requires H2_use_sketch=2");
      }
      if (H2_ID_radius < 2) {
        throw std::invalid_argument("H2_ID_radius must be at least 2");
      }
      if (H2_ID_proxy < 0 || H2_ID_proxy > 2) {
        throw std::invalid_argument("H2_ID_proxy must be 0, 1, or 2");
      }
      if (H2_ID_proxy == 1 && H2_ID_proxy_points <= 0) {
        throw std::invalid_argument(
            "H2_ID_proxy_points must be positive when H2_ID_proxy=1");
      }
      if (H2_ID_proxy == 2 && BACA_Batch <= 0) {
        throw std::invalid_argument(
            "BACA_Batch must be positive when H2_ID_proxy=2");
      }
      H2_options = butterfly::parse_program_options(
          Npo, Ndim, Locations, tolerance, reduction_threshold, Nmin_leaf, CA_level);
      H2_options.precon = static_cast<int>(std::llround(precon_d));
      H2_options.verbosity = static_cast<int>(std::llround(verbosity_d));
      H2_options.unstructured = H2_unstructured;
      {
        double sym_d = 1.0;
        c_bpack_getoption(option, "sym", &sym_d);
        H2_options.symmetric = std::llround(sym_d) != 0;
      }
      H2_options.use_sketch = H2_use_sketch;
      H2_options.lazy_schur = H2_lazy_schur;
      H2_options.compress_factors = H2_compress_factors;
      H2_options.gemm_split = H2_GEMM_split;
      H2_options.xrr_factor = H2_XRR_factor;
      H2_options.use_gpu = H2_use_gpu != 0 ? 1 : 0;
      H2_options.ca_staged_halo = H2_CA_staged_halo;
      H2_options.ca_owner_component = H2_CA_owner_component;
      H2_options.ca_pack_shared = H2_CA_pack_shared;
      H2_options.ca_solve = H2_CA_solve;
      H2_options.ca_owner_serial = H2_CA_owner_serial;
      H2_options.id_neighborhood_radius = H2_ID_radius;
      H2_options.id_proxy_mode = H2_ID_proxy;
      if (H2_ID_proxy == 1) {
        H2_options.id_proxy_points = H2_ID_proxy_points;
      } else if (H2_ID_proxy == 2) {
        H2_options.id_adaptive_batch = BACA_Batch;
      }
      H2_options.s_op = H2_s_op;
      if (H2_options.precon == 2) {
        H2_options.CA_level = H2_options.num_levels;
      }
      butterfly::validate_h2_backend_selection(H2_options);
      H2_solver->options = H2_options;
    } catch (const std::exception& e) {
      if (rank == 0) {
        std::cerr << "Argument error: " << e.what() << std::endl;
      }
      throw;
    }

    try {
      if (rank == 0 && H2_options.verbosity >= 1) {
        std::cout << "ButterflyPACK H2, number_type=" << butterfly::number_kind_to_string(H2_options.number_kind)
                  << ", dimension=" << H2_options.dimension
                  << ", reduction_threshold=" << H2_options.reduction_threshold
                  << ", h2_unstructured=" << H2_options.unstructured
                  << ", CA_level=" << H2_options.CA_level
                  << ", h2_use_sketch=" << H2_options.use_sketch
                  << ", h2_lazy_schur=" << H2_options.lazy_schur
                  << ", h2_gemm_split=" << H2_options.gemm_split
                  << ", h2_xrr_factor=" << H2_options.xrr_factor
                  << ", h2_use_gpu=" << H2_options.use_gpu
                  << ", h2_ca_staged_halo=" << H2_options.ca_staged_halo
                  << ", h2_ca_owner_component=" << H2_options.ca_owner_component
                  << ", h2_ca_pack_shared=" << H2_options.ca_pack_shared
                  << ", h2_ca_solve=" << H2_options.ca_solve
                  << ", h2_ca_owner_serial=" << H2_options.ca_owner_serial
                  << ", h2_id_radius=" << H2_options.id_neighborhood_radius
                  << ", h2_id_proxy=" << H2_options.id_proxy_mode
                  << ", h2_s_op=" << H2_options.s_op;
        if (H2_options.id_proxy_mode == 1) {
          std::cout << ", h2_id_proxy_points=" << H2_options.id_proxy_points;
        } else if (H2_options.id_proxy_mode == 2) {
          std::cout << ", baca_batch=" << H2_options.id_adaptive_batch;
        }
        std::cout << std::endl;
      }

	  // Setting up mesh and permutation variables
      std::vector<int> new2old;
	  int idxs = 0;
	  int idxe = -1;

      butterfly::h2_initiate<double, H2Data>(H2_solver, H2_options, Locations, rank, new2old, idxs, idxe);
	  // convert to perms and Npo_loc.  The ranks' ranges go with it, so the
	  // mesh carries a basis_group reproducing this distribution and a
	  // butterfly block built on it (c_bf_construct_init, for an F12/F21 or a
	  // Schur complement beside this matrix) is distributed the same way.
	  int nproc_h2 = 1;
	  MPI_Comm_size(H2_solver->comm, &nproc_h2);
	  std::vector<int> rank_idxs(nproc_h2, 0), rank_idxe(nproc_h2, -1);
	  MPI_Allgather(&idxs, 1, MPI_INT, rank_idxs.data(), 1, MPI_INT, H2_solver->comm);
	  MPI_Allgather(&idxe, 1, MPI_INT, rank_idxe.data(), 1, MPI_INT, H2_solver->comm);
	  c_bpack_set_mesh_h2(Npo, new2old.data(), &idxs, &idxe, &nproc_h2,
	                      rank_idxs.data(), rank_idxe.data(), msh);
	  *Npo_loc=idxe-idxs+1;
	  if (perms != nullptr) {
		std::copy(new2old.begin(), new2old.end(), perms);
	  }
	  new2old.clear();


    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(H2_solver->comm, 1);
    }

	H2_solver->kernel.entryeval_time_per_thread.assign(omp_get_max_threads(), 0.0);

	double zero = 0.0;
	c_bpack_setstats(stats, "Time_C_Mult_Wrapper", &zero);   // stats bare (already F2Cptr*), &zero
	c_bpack_wrap_h2(bmat, static_cast<C2Fptr>(H2_solver));   // *bmat now = c_loc(Fortran Bmatrix)
#else
	format7_requires_mpi("c_bpack_construct_init");
#endif
  }else{
	  c_bpack_construct_init_fortran(Npo, Ndim, Locations, nns, nlevel, tree, perms, Npo_loc, bmat, option, stats, msh, ker, ptree, C_FuncDistmn, C_FuncNearFar, C_QuantApp);
  }
}

/*
 * Distributed 64-bit ordering API.
 *
 * The caller owns an arbitrary rank-local subset of the global points and
 * supplies one-based int64_t IDs plus point-major coordinates.  This
 * constructor redistributes those points into ButterflyPACK's tree ordering
 * and retains the bidirectional map in the matrix handle.  Consequently,
 * c_bpack_solve and c_bpack_mult still consume and return vectors in the
 * caller's original rank-local order; their format-7 paths perform the
 * required redistribution transparently.
 *
 * These entry-point names are format-neutral so that other ButterflyPACK
 * formats can adopt the same interface later.  The current implementation is
 * intentionally restricted to format=7 and reports that restriction at run
 * time.  Global counts and IDs are 64-bit, while each MPI_Alltoallv peer count
 * and the existing solve/multiply local-size arguments remain INT_MAX-limited.
 */
void c_bpack_construct_init_distributed64(
    const int64_t* N_global, const int64_t* N_input_local,
    const int* Ndim, const int64_t* input_global_ids,
    const double* input_locations, const int* bounds_provided,
    const double* global_bounds, int64_t* N_internal_local,
    F2Cptr* bmat, F2Cptr* option, F2Cptr* stats, F2Cptr* msh,
    F2Cptr* ker, F2Cptr* ptree) {
#ifdef HAVE_MPI
  if (N_global == nullptr || N_input_local == nullptr || Ndim == nullptr ||
      bounds_provided == nullptr || N_internal_local == nullptr ||
      bmat == nullptr || option == nullptr || stats == nullptr ||
      msh == nullptr || ptree == nullptr) {
    throw std::invalid_argument(
        "c_bpack_construct_init_distributed64: null required argument");
  }

  int fcomm = 0;
  c_bpack_get_comm(ptree, &fcomm);
  MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(fcomm));
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  try {
    double format_d = 0.0;
    c_bpack_getoption(option, "format", &format_d);
    const int format = static_cast<int>(std::llround(format_d));
    if (format != 7) {
      throw std::runtime_error(
          "c_bpack_construct_init_distributed64 currently supports only "
          "format=7; support for other ButterflyPACK formats is TODO");
    }
    check_h2_symmetry_option(
        option, comm, "c_bpack_construct_init_distributed64");
    if (*N_global <= 0 || *N_input_local < 0) {
      throw std::invalid_argument(
          "c_bpack_construct_init_distributed64: invalid point count");
    }
    if (*bounds_provided != 0 && *bounds_provided != 1) {
      throw std::invalid_argument(
          "c_bpack_construct_init_distributed64: bounds_provided must be 0 or 1");
    }

    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    auto solver = std::make_unique<butterfly::H2<double, H2Data>>();
    solver->comm = comm;
    solver->options = read_h2_program_options64(
        option, *N_global, *Ndim);
    // global_bounds has 2*Ndim entries in [min0,max0,min1,max1,...]
    // order.  When bounds_provided is zero, the bounds are reduced from the
    // distributed coordinates instead.
    butterfly::bpack_initiate_distributed64(
        solver.get(), solver->options, *N_input_local,
        input_global_ids, input_locations, *bounds_provided != 0,
        global_bounds);
    *N_internal_local =
        solver->distributed_layout->internal_local_size;

    c_bpack_set_mesh_distributed64(
        N_global, N_input_local, N_internal_local,
        &solver->distributed_layout->internal_global_start, msh);
    solver->kernel.entryeval_time_per_thread.assign(
        omp_get_max_threads(), 0.0);
    double zero = 0.0;
    c_bpack_setstats(stats, "Time_C_Mult_Wrapper", &zero);
    if (ker != nullptr) *ker = nullptr;
    c_bpack_wrap_h2(bmat, static_cast<C2Fptr>(solver.release()));
  } catch (const std::exception& error) {
    if (rank == 0) {
      std::cerr << "c_bpack_construct_init_distributed64: "
                << error.what() << std::endl;
    }
    MPI_Abort(comm, 1);
    throw;
  }
#else
  (void)N_global;
  (void)N_input_local;
  (void)Ndim;
  (void)input_global_ids;
  (void)input_locations;
  (void)bounds_provided;
  (void)global_bounds;
  (void)N_internal_local;
  (void)bmat;
  (void)option;
  (void)stats;
  (void)msh;
  (void)ker;
  (void)ptree;
  format7_requires_mpi("c_bpack_construct_init_distributed64");
#endif
}

/*
 * Install coordinate-aware matrix-entry callbacks for a matrix created by
 * c_bpack_construct_init_distributed64.  Scalar callbacks receive one-based
 * global IDs and the two Ndim-coordinate points.  Block callbacks receive
 * point-major row/column coordinates and fill a column-major matrix using the
 * supplied 64-bit leading dimension.  C_QuantApp is passed through unchanged
 * as the callback context.
 */
void c_bpack_construct_element_compute_distributed64(
    F2Cptr* bmat, F2Cptr* option, F2Cptr* stats, F2Cptr* msh,
    F2Cptr* ker, F2Cptr* ptree,
    c_bpack_func_zmn64 C_FuncZmn,
    c_bpack_func_zmn_block64 C_FuncZmnBlock,
    C2Fptr C_QuantApp) {
  (void)stats;
  (void)msh;
  (void)ker;
  (void)ptree;
#ifdef HAVE_MPI
  double format_d = 0.0;
  c_bpack_getoption(option, "format", &format_d);
  if (static_cast<int>(std::llround(format_d)) != 7) {
    throw std::runtime_error(
        "c_bpack_construct_element_compute_distributed64 currently "
        "supports only format=7; support for other formats is TODO");
  }
  if (C_FuncZmn == nullptr) {
    throw std::invalid_argument(
        "c_bpack_construct_element_compute_distributed64: entry callback is null");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  auto* solver = get_h2_solver<H2Data>(
      bmat, "c_bpack_construct_element_compute_distributed64");
  check_h2_symmetry_option(
      option, solver->comm,
      "c_bpack_construct_element_compute_distributed64");
  if (!solver->distributed_layout) {
    throw std::runtime_error(
        "c_bpack_construct_element_compute_distributed64 requires the "
        "distributed64 constructor");
  }
  solver->kernel.kernel = nullptr;
  solver->kernel.block_kernel = nullptr;
  solver->kernel.kernel64 = C_FuncZmn;
  double elem_extract_d = 0.0;
  c_bpack_getoption(option, "elem_extract", &elem_extract_d);
  const int elem_extract =
      static_cast<int>(std::llround(elem_extract_d));
  solver->kernel.block_kernel64 =
      elem_extract == 2 ? C_FuncZmnBlock : nullptr;
  solver->kernel.quant = C_QuantApp;
  solver->kernel.dimension = solver->options.dimension;
  solver->kernel.symmetric = solver->options.symmetric;
  MPI_Comm_rank(solver->comm, &solver->kernel.block_callback_pid);
  solver->kernel.register_level_coordinates(
      solver->tree->levels[solver->tree->num_levels - 1]);
#else
  (void)bmat;
  (void)option;
  (void)C_FuncZmn;
  (void)C_FuncZmnBlock;
  (void)C_QuantApp;
  format7_requires_mpi(
      "c_bpack_construct_element_compute_distributed64");
#endif
}


void c_bpack_construct_element_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,
	F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),
	void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr),
	C2Fptr C_QuantApp){
  // these functions are important to define
  // C_FuncZmn: returns value at (i,j)th element of matrix -- need to update this to work for kernel
  // C_FuncZmnBlock: returns a block, low priority

  // To do: need to define stats and msh or else need to redefine c_bpack_delete for these to pointers
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
	using H2Data = typename butterfly::fmm_data<C_DT>::type;
	static_assert(std::is_same_v<C2Fptr, void*>, "H2::kernel assumes C2Fptr == void*; update butterfly_integration.hpp if this changes");
	void* H2_raw = nullptr;
	c_bpack_get_h2(*bmat, &H2_raw);
	butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(
	  option, H2_solver->comm, "c_bpack_construct_element_compute");
	H2_solver->kernel.kernel = C_FuncZmn;
	H2_solver->kernel.kernel64 = nullptr;
	H2_solver->kernel.block_kernel64 = nullptr;
	H2_solver->kernel.dimension = H2_solver->options.dimension;
	H2_solver->kernel.symmetric = H2_solver->options.symmetric;
	double elem_extract_d = 0.0;
	c_bpack_getoption(option, "elem_extract", &elem_extract_d);
	const int elem_extract = static_cast<int>(std::llround(elem_extract_d));
	H2_solver->kernel.block_kernel =
	  (elem_extract == 2) ? C_FuncZmnBlock : nullptr;
	// the rank within this matrix's own communicator, not the world's: the
	// block callback's pmaps name the owner of an intersection, and a caller
	// that builds an H2 on a sub-communicator (a STRUMPACK front on
	// MPI_COMM_SELF inside a larger job) resolves that id against the
	// matrix's process tree, where a world rank is out of range
	MPI_Comm_rank(H2_solver->comm, &H2_solver->kernel.block_callback_pid);
	H2_solver->kernel.quant = C_QuantApp;
	// elem_extract=1 with a block callback: the callback is collective (it
	// may communicate, as c_bpack_extractelement_bridge does), so the blocks
	// are fetched in rounds of all ranks (collective_entries.hpp); otherwise
	// each rank calls its callbacks when it needs an entry
	// (H2_use_gpu: the device reads the same blocks through the library's
	// evaluator of them, gpu/collective_evaluator.hpp, in place of the
	// application's)
#if defined(H2_HAVE_GPU) && defined(DAT) && (DAT == 0 || DAT == 1)
	if (const auto* e = fmm::gpu::evaluator_of(H2_solver->kernel.gpu_evaluator); e != nullptr && e->collective()) {
	  H2_solver->kernel.gpu_evaluator.reset();
	}
#endif
	H2_solver->kernel.collective.reset();
	if (elem_extract == 1 && C_FuncZmnBlock != nullptr) {
	  // (a nonsymmetric H2: each orientation of a pair of boxes fetched apart)
	  H2_solver->kernel.collective =
	    std::make_shared<fmm::CollectiveEntries<H2Data, typename butterfly::H2Kernel<double, H2Data>::kerData>>(
	      H2_solver->comm, C_FuncZmnBlock, C_QuantApp, H2_solver->kernel.block_callback_pid,
	      H2_solver->options.symmetric);
	  // (c_bpack_set_sketch_callback: the ID sketches of kernel rows fetched
	  // sketched; a nonsymmetric H2 needs c_bpack_set_sketch_transpose_callback
	  // too, else it fetches their blocks, CollectiveEntries::set_sketch)
	  H2_solver->kernel.collective->set_sketch(H2_solver->kernel.sketch_callback);
	  H2_solver->kernel.collective->set_sketch_transposed(H2_solver->kernel.sketch_transpose_callback);
#if defined(H2_HAVE_GPU) && defined(DAT) && (DAT == 0 || DAT == 1)
	  if (H2_solver->options.use_gpu != 0) {
	    H2_solver->kernel.gpu_evaluator =
	      fmm::gpu::make_collective_evaluator<H2Data>(H2_solver->kernel.collective.get());
	  }
#endif
	}
#else
	format7_requires_mpi("c_bpack_construct_element_compute");
#endif
  }else{
	c_bpack_construct_element_compute_fortran(bmat, option, stats, msh, ker, ptree, C_FuncZmn, C_FuncZmnBlock, C_QuantApp);
  }
}


/* Return counts for both orderings and the one-based first internal index. */
void c_bpack_get_distributed_layout64(
    F2Cptr* bmat, int64_t* N_global, int64_t* N_input_local,
    int64_t* N_internal_local, int64_t* internal_global_start) {
#ifdef HAVE_MPI
  if (N_global == nullptr || N_input_local == nullptr ||
      N_internal_local == nullptr || internal_global_start == nullptr) {
    throw std::invalid_argument(
        "c_bpack_get_distributed_layout64: null output argument");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const auto* layout = get_distributed_layout<H2Data>(
      bmat, "c_bpack_get_distributed_layout64");
  *N_global = layout->global_size;
  *N_input_local = layout->input_local_size;
  *N_internal_local = layout->internal_local_size;
  *internal_global_start = layout->internal_global_start;
#else
  (void)bmat;
  (void)N_global;
  (void)N_input_local;
  (void)N_internal_local;
  (void)internal_global_start;
  format7_requires_mpi("c_bpack_get_distributed_layout64");
#endif
}

/*
 * The application's GPU evaluator of the entries (GPU_INTERFACE/bpack_gpu.h,
 * doc/gpu_kernels.md): kept on the H2 solver of a format 7 matrix, else on
 * the matrix's GPU state (HODLR).  Only the double and double complex
 * libraries of a GPU build have a GPU backend; elsewhere it is ignored.
 */
#if defined(HAVE_MPI) && defined(H2_HAVE_GPU) && defined(DAT) && (DAT == 0 || DAT == 1)
#define BPACK_GPU_EVALUATORS 1
#endif

#ifdef BPACK_GPU_EVALUATORS
/* The GPU state of a non-H2 matrix, created on first use. */
static void* bpack_gpu_state(F2Cptr* bmat, const char* caller) {
  void* gpu = nullptr;
  c_bpack_get_gpu(*bmat, &gpu);
  if (gpu == nullptr) {
    c_bpack_gpu_create(&gpu);
    if (gpu == nullptr) throw std::runtime_error(std::string(caller) + ": no GPU state");
    c_bpack_set_gpu(*bmat, gpu);
  }
  return gpu;
}

template<typename Make>
static void bpack_set_gpu_evaluator(F2Cptr* bmat, const char* caller, Make make) {
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  if (bmat == nullptr || *bmat == nullptr) {
    throw std::invalid_argument(std::string(caller) + ": null matrix handle");
  }
  const std::shared_ptr<fmm::gpu::Evaluator> evaluator = make(static_cast<int>(sizeof(H2Data)));
  void* h2 = nullptr;
  c_bpack_get_h2(*bmat, &h2);
  if (h2 != nullptr) {
    get_h2_solver<H2Data>(bmat, caller)->kernel.gpu_evaluator = evaluator;
  } else {
    c_bpack_gpu_set_evaluator(bpack_gpu_state(bmat, caller), &evaluator);
  }
}
#endif

void c_bpack_set_gpu_entry_launchers(F2Cptr* bmat, const bpack_gpu_entry_launchers* launchers) {
#ifdef BPACK_GPU_EVALUATORS
  if (launchers == nullptr) throw std::invalid_argument("c_bpack_set_gpu_entry_launchers: null argument");
  bpack_set_gpu_evaluator(bmat, "c_bpack_set_gpu_entry_launchers", [&](int scalar_bytes) {
    return fmm::gpu::Evaluator::from_launchers(*launchers, scalar_bytes);
  });
#else
  (void)bmat;
  (void)launchers;
#endif
}

void c_bpack_set_gpu_entry_source(F2Cptr* bmat, const char* source, const double* params, const int* nparams,
                                  const int* flags) {
#ifdef BPACK_GPU_EVALUATORS
  if (source == nullptr || flags == nullptr || (nparams != nullptr && *nparams > 0 && params == nullptr)) {
    throw std::invalid_argument("c_bpack_set_gpu_entry_source: null argument");
  }
  const int count = nparams != nullptr ? std::max(*nparams, 0) : 0;
  bpack_set_gpu_evaluator(bmat, "c_bpack_set_gpu_entry_source", [&](int scalar_bytes) {
    return fmm::gpu::Evaluator::from_source(source, params, count, *flags, scalar_bytes);
  });
#else
  (void)bmat;
  (void)source;
  (void)params;
  (void)nparams;
  (void)flags;
#endif
}

void c_bpack_set_gpu_block_evaluator(F2Cptr* bmat, bpack_gpu_block_evaluator evaluate, void* user,
                                     const int* flags) {
#ifdef BPACK_GPU_EVALUATORS
  if (evaluate == nullptr || flags == nullptr) {
    throw std::invalid_argument("c_bpack_set_gpu_block_evaluator: null argument");
  }
  bpack_set_gpu_evaluator(bmat, "c_bpack_set_gpu_block_evaluator", [&](int scalar_bytes) {
    return fmm::gpu::Evaluator::from_block(evaluate, user, *flags, scalar_bytes);
  });
#else
  (void)bmat;
  (void)evaluate;
  (void)user;
  (void)flags;
#endif
}

void c_bpack_set_gpu_list_evaluator(F2Cptr* bmat, bpack_gpu_list_evaluator evaluate, void* user,
                                    const int* flags) {
#ifdef BPACK_GPU_EVALUATORS
  if (evaluate == nullptr || flags == nullptr) {
    throw std::invalid_argument("c_bpack_set_gpu_list_evaluator: null argument");
  }
  bpack_set_gpu_evaluator(bmat, "c_bpack_set_gpu_list_evaluator", [&](int scalar_bytes) {
    return fmm::gpu::Evaluator::from_list(evaluate, user, *flags, scalar_bytes);
  });
#else
  (void)bmat;
  (void)evaluate;
  (void)user;
  (void)flags;
#endif
}

/*
 * The warm-up of the registered GPU evaluator now, on the GPU the matrix will
 * use, rather than before its first level (fmm::gpu::Evaluator::warm_up):
 * for H2 with H2_use_gpu, after the registration.  seconds[0]: its time on
 * this rank, seconds[1]: of which the NVRTC compile of an entry evaluator
 * given as source.  Collective over the matrix's ranks.  HODLR needs none:
 * its construction warms up first.
 */
void c_bpack_gpu_warm_up(F2Cptr* bmat, double* seconds) {
  if (seconds == nullptr) throw std::invalid_argument("c_bpack_gpu_warm_up: null argument");
  seconds[0] = seconds[1] = 0.0;
#ifdef BPACK_GPU_EVALUATORS
  void* h2 = nullptr;
  if (bmat != nullptr && *bmat != nullptr) c_bpack_get_h2(*bmat, &h2);
  if (h2 == nullptr) return;
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  auto* solver = get_h2_solver<H2Data>(bmat, "c_bpack_gpu_warm_up");
  if (solver->options.use_gpu == 0 || !solver->kernel.gpu_evaluator) return;
  // (the device share of the ranks and the operator, as a compression or
  // factorization starts, before the evaluator takes device memory)
  fmm::gpu::begin_operator_build(solver->tree.get(), solver->tree->comm);
  seconds[0] = fmm::gpu::warm_up_registered_evaluator<H2Data>(solver->tree.get(), solver->kernel.gpu_evaluator);
  seconds[1] = fmm::gpu::evaluator_of(solver->kernel.gpu_evaluator)->compile_seconds();
#else
  (void)bmat;
#endif
}

/*
 * Return original one-based point IDs in internal tree order.  The requested
 * rank-local offset is also one-based.
 */
void c_bpack_get_internal_global_ids64(
    F2Cptr* bmat, const int64_t* internal_local_offset,
    const int64_t* count, int64_t* global_ids) {
#ifdef HAVE_MPI
  if (internal_local_offset == nullptr || count == nullptr ||
      (*count > 0 && global_ids == nullptr)) {
    throw std::invalid_argument(
        "c_bpack_get_internal_global_ids64: null argument");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const auto* layout = get_distributed_layout<H2Data>(
      bmat, "c_bpack_get_internal_global_ids64");
  const int64_t offset = *internal_local_offset - 1;
  if (*internal_local_offset < 1 || *count < 0 ||
      offset > layout->internal_local_size ||
      *count > layout->internal_local_size - offset) {
    throw std::out_of_range(
        "c_bpack_get_internal_global_ids64: local range is invalid");
  }
  if (*count == 0) return;
  std::copy_n(layout->internal_global_ids.data() + offset,
              static_cast<size_t>(*count), global_ids);
#else
  (void)bmat;
  (void)internal_local_offset;
  (void)count;
  (void)global_ids;
  format7_requires_mpi("c_bpack_get_internal_global_ids64");
#endif
}

/*
 * Describe where each caller-local point moved.  Owner ranks are zero-based
 * MPI ranks; IDs, input_local_offset, and internal_global_indices are
 * one-based.
 */
void c_bpack_get_input_to_internal_map64(
    F2Cptr* bmat, const int64_t* input_local_offset,
    const int64_t* count, int64_t* global_ids,
    int* internal_owner_ranks, int64_t* internal_global_indices) {
#ifdef HAVE_MPI
  if (input_local_offset == nullptr || count == nullptr ||
      (*count > 0 && (global_ids == nullptr ||
                      internal_owner_ranks == nullptr ||
                      internal_global_indices == nullptr))) {
    throw std::invalid_argument(
        "c_bpack_get_input_to_internal_map64: null argument");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const auto* layout = get_distributed_layout<H2Data>(
      bmat, "c_bpack_get_input_to_internal_map64");
  const int64_t offset = *input_local_offset - 1;
  if (*input_local_offset < 1 || *count < 0 ||
      offset > layout->input_local_size ||
      *count > layout->input_local_size - offset) {
    throw std::out_of_range(
        "c_bpack_get_input_to_internal_map64: local range is invalid");
  }
  for (int64_t index = 0; index < *count; ++index) {
    const size_t source = static_cast<size_t>(offset + index);
    global_ids[index] = layout->input_global_ids[source];
    internal_owner_ranks[index] = layout->input_internal_owner[source];
    internal_global_indices[index] =
        layout->input_internal_global_index[source];
  }
#else
  (void)bmat;
  (void)input_local_offset;
  (void)count;
  (void)global_ids;
  (void)internal_owner_ranks;
  (void)internal_global_indices;
  format7_requires_mpi("c_bpack_get_input_to_internal_map64");
#endif
}

/*
 * Explicitly redistribute a column-major local multivector from caller input
 * ownership to internal tree ownership.  All nrhs columns are packed into one
 * collective rather than issuing one collective per column.
 */
void c_bpack_input_to_internal(
    F2Cptr* bmat, const int* nrhs,
    const C_DT* input_values, C_DT* internal_values) {
#ifdef HAVE_MPI
  if (nrhs == nullptr) {
    throw std::invalid_argument("c_bpack_input_to_internal: nrhs is null");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const auto* layout = get_distributed_layout<H2Data>(
      bmat, "c_bpack_input_to_internal");
  layout->input_to_internal(
      reinterpret_cast<const H2Data*>(input_values),
      reinterpret_cast<H2Data*>(internal_values), *nrhs);
#else
  (void)bmat;
  (void)nrhs;
  (void)input_values;
  (void)internal_values;
  format7_requires_mpi("c_bpack_input_to_internal");
#endif
}

/*
 * The batched-extraction profile of this process (unstructured/extraction.hpp,
 * core/extraction_profile.hpp): seconds of the eight steps (resolve, tasks,
 * vectors up, column side, compute, pack, route, unpack), then the route's
 * MPI seconds, the entries computed, the bytes routed and the calls, summed
 * over every compressed H2 this process served since the start (or since a
 * CollectiveEntries::report, which resets it).  Local, not collective.
 */
void c_bpack_get_extraction_profile(double* vals) {
#ifdef HAVE_MPI
  const fmm::ExtractionProfile& p = fmm::extraction_profile();
  for (int s = 0; s < fmm::ExtractionProfile::kSteps; ++s) vals[s] = p.seconds[s];
  vals[8] = p.route_mpi;
  vals[9] = p.entries_computed;
  vals[10] = p.bytes_routed;
  vals[11] = static_cast<double>(p.calls);
#else
  for (int s = 0; s < 12; ++s) vals[s] = 0.0;
#endif
}

/* Inverse of c_bpack_input_to_internal, with the same batched layout. */
void c_bpack_internal_to_input(
    F2Cptr* bmat, const int* nrhs,
    const C_DT* internal_values, C_DT* input_values) {
#ifdef HAVE_MPI
  if (nrhs == nullptr) {
    throw std::invalid_argument("c_bpack_internal_to_input: nrhs is null");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const auto* layout = get_distributed_layout<H2Data>(
      bmat, "c_bpack_internal_to_input");
  layout->internal_to_input(
      reinterpret_cast<const H2Data*>(internal_values),
      reinterpret_cast<H2Data*>(input_values), *nrhs);
#else
  (void)bmat;
  (void)nrhs;
  (void)internal_values;
  (void)input_values;
  format7_requires_mpi("c_bpack_internal_to_input");
#endif
}


void c_bpack_factor(F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*msh){

  // Correspond to hierarchical_factorization_parallel, arguments
  // tree:
  // kernel: kernel from factorizer
  // tolerance:
  // is_symmetric: bool -- works for general helmholtz and V3D
  // is_hermitian: bool, not supported right now
  // factorization_method: provided by factorizer, factorization_type
  // unit_proxy_points:
  // num_proxy:
  // proxy_radius:

  // bmat: can contain tree, and kernel function
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(option, H2_solver->comm, "c_bpack_factor");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    try {
      sync_h2_verbosity(option, H2_solver);
      double precon_d = 1.0;
	  c_bpack_getoption(option, "precon", &precon_d);
	  H2_solver->options.precon = static_cast<int>(std::llround(precon_d));
	  if (H2_solver->options.precon == 2) {
		compress_h2_and_update_stats(H2_solver, stats);
	  } else {
		double factorization_time = 0.0;
		double entryeval_time = 0.0;
		butterfly::dispatch_h2_factorization(
		  H2_solver, &factorization_time, &entryeval_time);
		c_bpack_setstats(stats, "Time_Factor", &factorization_time);
		c_bpack_setstats(stats, "Time_Entry", &entryeval_time);

		double rank_max = static_cast<double>(H2_solver->last_factor_rankmax);
		c_bpack_setstats(stats, "Rank_max", &rank_max);

		// bytes / 1024e3, as the Fortran paths.  Mem_Fill stays unset: an H2
		// compresses and factors in one pass and keeps no separate compressed
		// form, so its whole size is the factor's.
		double factorization_memory_MB = H2_solver->factorization_memory / 1024.0e3;
		c_bpack_setstats(stats, "Mem_Factor", &factorization_memory_MB);
	  }

    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        throw;
    }
#else
	format7_requires_mpi("c_bpack_factor");
#endif
  }else{
	c_bpack_factor_fortran(bmat, option, stats, ptree, msh);
  }
}

void c_bpack_solve(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree){
  // correspond to hierarchical_solve_parallel, arguments:
  //   tree:
  //   rhs: pass in b
  //   solve_data: SolveDataRequest type, accumulates distributes solution / communication requests during the tree level sweep
  //   verbose: printing

  // and gather_solution_to_root, arguments:
  //   tree
  //   solve_data: pass in from hierarchical_solve_parallel
  //   solution: pass into x
  //   aggregated_rhs:

  // Ax = b
  // x: final solution
  // b: provided rhs
  // Hloc:
  // Nrhs: number of rhs columns
  // bmat: factored matrix, do we need this? because it's in tree
  // option
  // stat
  // ptree

  // note to self: figure out aggregated_rhs, solution, solve_data; and mpi stuff (mpi stuff probably ask Tianyu)
  // need to redistribute x into H2, and then call mul_parallel, then extract mul_data, the nredistribute to Butterfly
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
    if (*Nrhs <= 0) {
      throw std::invalid_argument(
        "c_bpack_solve (format 7): Nrhs must be positive");
    }
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(option, H2_solver->comm, "c_bpack_solve");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    try {
      const int verbosity = sync_h2_verbosity(option, H2_solver);
      double t0 = MPI_Wtime();
	  // A distributed64 handle exposes the caller's input ownership at this
	  // public boundary: the RHS batch goes to tree ownership and the solution
	  // comes back below.
	  std::vector<H2Data> rhs;
	  const int internal_nloc = h2_values_to_internal(
		H2_solver, reinterpret_cast<const H2Data*>(b), *Nloc, *Nrhs, rhs,
		"c_bpack_solve");
	  std::vector<H2Data> internal_solution;
	  double precon_d = 1.0;
	  c_bpack_getoption(option, "precon", &precon_d);
	  H2_solver->options.precon = static_cast<int>(std::llround(precon_d));

	  if (H2_solver->options.precon == 2) {
		if (H2_solver->build_state == butterfly::H2BuildState::UNBUILT) {
		  throw std::runtime_error(
			"c_bpack_solve (format 7): call c_bpack_factor before solving");
		}
		if (H2_solver->build_state != butterfly::H2BuildState::H2_COMPRESSED) {
		  throw std::runtime_error(
			"c_bpack_solve (format 7): precon=2 requires a compression-only H2 representation");
		}

		double tolerance = 0.0;
		double max_iterations_d = 0.0;
		c_bpack_getoption(option, "tol_itersol", &tolerance);
		c_bpack_getoption(option, "n_iter", &max_iterations_d);
		int iterations = 0;
		double residual = 0.0;
		butterfly::dispatch_h2_iterative_solve(
		  H2_solver, rhs, internal_solution, *Nrhs,
		  tolerance, static_cast<int>(std::llround(max_iterations_d)),
		  &iterations, &residual, verbosity >= 1);
	  } else {
		if (H2_solver->build_state != butterfly::H2BuildState::RS_FACTORIZED) {
		  throw std::runtime_error(
			"c_bpack_solve (format 7): direct solve requires c_bpack_factor first");
		}
		h2_apply_factorized_inverse(
		  H2_solver, rhs, internal_solution, internal_nloc, *Nrhs, verbosity);
	  }

	  if (internal_solution.size() !=
		  static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Nrhs)) {
		throw std::runtime_error(
		  "c_bpack_solve (format 7): internal solution length mismatch");
	  }
	  h2_values_from_internal(
		H2_solver, internal_solution, reinterpret_cast<H2Data*>(x), *Nrhs);

	  double t_solve = MPI_Wtime() - t0;
	  MPI_Allreduce(MPI_IN_PLACE, &t_solve, 1, MPI_DOUBLE, MPI_MAX, H2_solver->comm);
	  c_bpack_setstats(stats, "Time_Solve", &t_solve);
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        // Rethrown, not MPI_Abort(MPI_COMM_WORLD): these entries are called
        // from C++ (no Fortran frame to unwind through), and a host solver
        // compressing one front of many decides for itself what a failure
        // means.  c_bpack_h2_inv_mult and c_bpack_h2_mult below catch again
        // for the Fortran callers.
        throw;
    }
#else
	format7_requires_mpi("c_bpack_solve");
#endif
  }else{
	c_bpack_solve_fortran(x, b, Nloc, Nrhs, bmat, option, stats, ptree);
  }
}

void c_bpack_mult(char const * trans, C_DT const * xin,
	C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol,
	F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree){

  // F * xin = xout, where F is the approximated matrix from hierarchical decomposition
  //

  // can call fft_matvec for uniform grid



  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(option, H2_solver->comm, "c_bpack_mult");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    if (*Ncol <= 0) {
      throw std::invalid_argument(
        "c_bpack_mult (format 7): Ncol must be positive");
    }

    try {
      // F x, F^T x (= F x for a transpose-symmetric H2), F^H x = conj(F^T
      // conj(x)) (h2_trans_operation)
      const char op = h2_trans_operation(trans, "c_bpack_mult (format 7)");
      const bool transposed = h2_transposed(H2_solver, op);
      if (*Ninloc != *Noutloc) {
        throw std::invalid_argument(
          "c_bpack_mult (format 7): Ninloc must equal Noutloc");
      }
      const int verbosity = sync_h2_verbosity(option, H2_solver);
      double t0 = MPI_Wtime();
      double precon_d = 1.0;
	  c_bpack_getoption(option, "precon", &precon_d);
      H2_solver->options.precon = static_cast<int>(std::llround(precon_d));

      if (H2_solver->build_state == butterfly::H2BuildState::UNBUILT) {
		throw std::runtime_error(
		  "c_bpack_mult (format 7): call c_bpack_factor before multiplication");
      }

	  // Keep the public multiply contract in caller input order for a
	  // distributed64 handle; the H2 multiply itself operates in tree order.
	  std::vector<H2Data> lhs;
	  int internal_nloc = h2_values_to_internal(
		H2_solver, reinterpret_cast<const H2Data*>(xin), *Ninloc, *Ncol, lhs,
		"c_bpack_mult");
	  if (op == 'C') conjugate_h2_values(lhs);
	  std::vector<H2Data> compressed_output;
	  std::vector<std::vector<fmm::SolveDataRequest<double, H2Data>>> mul_data;

	  if (H2_solver->build_state == butterfly::H2BuildState::H2_COMPRESSED) {
		butterfly::dispatch_h2_compressed_multiply(
		  H2_solver, lhs, compressed_output,
		  *Ncol, verbosity >= 1, transposed);
	  } else {
		mul_data.resize(H2_solver->options.num_levels);
		butterfly::dispatch_h2_multiply(
		  H2_solver, lhs, mul_data,
		  *Ncol, verbosity >= 1, transposed);
	  }

	  std::vector<H2Data> internal_output;
	  if (H2_solver->build_state == butterfly::H2BuildState::H2_COMPRESSED) {
		if (compressed_output.size() !=
		    static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Ncol)) {
		  throw std::runtime_error(
			"c_bpack_mult (format 7): compression-only output length mismatch");
		}
		internal_output = std::move(compressed_output);
	  } else {
		internal_output.resize(
		  static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Ncol));
		butterfly::gather_local_solution(
		  H2_solver->tree.get(), mul_data,
		  internal_output.data(), &internal_nloc, *Ncol);
	  }
	  if (op == 'C') conjugate_h2_values(internal_output);
	  h2_values_from_internal(
		H2_solver, internal_output, reinterpret_cast<H2Data*>(xout), *Ncol);

	  double t_mult = MPI_Wtime() - t0;
	  MPI_Allreduce(MPI_IN_PLACE, &t_mult, 1, MPI_DOUBLE, MPI_MAX,
	                H2_solver->comm);

	  double prev = 0.0;
	  c_bpack_getstats(stats, "Time_C_Mult_Wrapper", &prev);
	  double total = prev + t_mult;
	  c_bpack_setstats(stats, "Time_C_Mult_Wrapper", &total);

    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        // Rethrown, not MPI_Abort(MPI_COMM_WORLD): these entries are called
        // from C++ (no Fortran frame to unwind through), and a host solver
        // compressing one front of many decides for itself what a failure
        // means.  c_bpack_h2_inv_mult and c_bpack_h2_mult below catch again
        // for the Fortran callers.
        throw;
    }
#else
	format7_requires_mpi("c_bpack_mult");
#endif
  }else{
	c_bpack_mult_fortran(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  }
}

void c_bpack_inv_mult(char const * trans, C_DT const * xin,
	C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol,
	F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree){

  // xout = F^{-1} xin ('N'), F^{-T} xin ('T') or F^{-H} xin ('C'), F the
  // factorized matrix; repeated on blocks of vectors (Schur sampling), so
  // quiet below verbosity 1.  Format 7 applies the H2 factorization (the
  // direct solve's sweeps, transposed for 'T' of a nonsymmetric H2; 'C'
  // conjugates, h2_trans_operation).
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(option, H2_solver->comm, "c_bpack_inv_mult");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    try {
      const char op = h2_trans_operation(trans, "c_bpack_inv_mult (format 7)");
      const bool transposed = h2_transposed(H2_solver, op);
      if (*Ncol <= 0) {
        throw std::invalid_argument(
          "c_bpack_inv_mult (format 7): Ncol must be positive");
      }
      if (*Ninloc != *Noutloc) {
        throw std::invalid_argument(
          "c_bpack_inv_mult (format 7): Ninloc must equal Noutloc");
      }
      if (H2_solver->build_state != butterfly::H2BuildState::RS_FACTORIZED) {
        throw std::runtime_error(
          "c_bpack_inv_mult (format 7): applies the factorization; call "
          "c_bpack_factor with precon other than 2 first");
      }
      const int verbosity = sync_h2_verbosity(option, H2_solver);
      double t0 = MPI_Wtime();

      std::vector<H2Data> rhs;
      const int internal_nloc = h2_values_to_internal(
        H2_solver, reinterpret_cast<const H2Data*>(xin), *Ninloc, *Ncol, rhs,
        "c_bpack_inv_mult");
      if (op == 'C') conjugate_h2_values(rhs);
      const double t_in = MPI_Wtime();
      std::vector<H2Data> solution;
      double t_apply[2] = {0.0, 0.0};
      h2_apply_factorized_inverse(
        H2_solver, rhs, solution, internal_nloc, *Ncol,
        verbosity >= 1 ? verbosity : -1, t_apply, transposed);
      const double t_out = MPI_Wtime();
      if (op == 'C') conjugate_h2_values(solution);
      h2_values_from_internal(
        H2_solver, solution, reinterpret_cast<H2Data*>(xout), *Ncol);

      // (as the Fortran inv_mult: this call's time)
      double t_inv_mult = MPI_Wtime() - t0;
      if (verbosity >= 1 && rank == 0) {
        std::cout << "c_bpack_inv_mult (format 7, '" << op << "', " << *Ncol
                  << " columns): " << t_inv_mult << " s on rank 0: to the tree's order "
                  << t_in - t0 << ", sweeps " << t_apply[0] << ", gather "
                  << t_apply[1] << ", back " << MPI_Wtime() - t_out << std::endl;
      }
      MPI_Allreduce(MPI_IN_PLACE, &t_inv_mult, 1, MPI_DOUBLE, MPI_MAX,
                    H2_solver->comm);
      c_bpack_setstats(stats, "Time_C_Mult", &t_inv_mult);
      c_bpack_setstats(stats, "Time_C_Mult_Wrapper", &t_inv_mult);
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        // Rethrown, not MPI_Abort(MPI_COMM_WORLD): these entries are called
        // from C++ (no Fortran frame to unwind through), and a host solver
        // compressing one front of many decides for itself what a failure
        // means.  c_bpack_h2_inv_mult and c_bpack_h2_mult below catch again
        // for the Fortran callers.
        throw;
    }
#else
	format7_requires_mpi("c_bpack_inv_mult");
#endif
  }else{
	c_bpack_inv_mult_fortran(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  }
}

// Fortran-facing wrappers of c_bpack_inv_mult and c_bpack_mult.  BPACK_Inv_Mult
// and BPACK_Mult have no format-7 case, so C_BPACK_Inv_Mult and C_BPACK_Mult
// (BPACK_wrapper.f90) send an H2 matrix here instead; c_bf_schur_deterministic
// reaches its F11 inverse this way.  An exception must not unwind through the
// Fortran frame that called this, so a failure is fatal here.
static void h2_fortran_bridge_failed(const char* who, const std::exception& e) {
  std::cerr << who << ": " << e.what() << std::endl;
  std::cerr.flush();
#ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  std::abort();
}

void c_bpack_h2_inv_mult(char const * trans, C_DT const * xin,
	C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol,
	F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree){
  try {
    c_bpack_inv_mult(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  } catch (const std::exception& e) {
    h2_fortran_bridge_failed("c_bpack_h2_inv_mult", e);
  }
}

void c_bpack_h2_mult(char const * trans, C_DT const * xin,
	C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol,
	F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree){
  try {
    c_bpack_mult(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  } catch (const std::exception& e) {
    h2_fortran_bridge_failed("c_bpack_h2_mult", e);
  }
}

#ifdef HAVE_MPI
// The rank owning each intersection (pmaps of 1 x 1 grids only).
static std::vector<int> single_rank_owners(int ninter, const int* pgidx, int npmap, const int* pmaps,
                                           const char* caller) {
  std::vector<int> owner(static_cast<size_t>(std::max(ninter, 0)));
  for (int n = 0; n < ninter; ++n) {
    const int pg = pgidx[n];
    if (pg < 0 || pg >= npmap) {
      throw std::invalid_argument(std::string(caller) + ": a process group out of range");
    }
    if (pmaps[pg] != 1 || pmaps[npmap + pg] != 1) {
      throw std::invalid_argument(std::string(caller) +
        ": an intersection owned by more than one rank (1 x 1 process grids only)");
    }
    owner[static_cast<size_t>(n)] = pmaps[2 * npmap + pg];
  }
  return owner;
}

template<typename H2Data>
static butterfly::H2<double, H2Data>* h2_of(F2Cptr* bmat) {
  if (bmat == nullptr || *bmat == nullptr) return nullptr;
  void* raw = nullptr;
  c_bpack_get_h2(*bmat, &raw);
  return static_cast<butterfly::H2<double, H2Data>*>(raw);
}

// The blocks of a list of intersections from a compression-only
// unstructured Color H2 (format 7; unstructured/extraction.hpp), each to one
// rank of caller_comm (pmaps of 1 x 1 grids), collective over caller_comm;
// the H2 on a subset of its ranks, bmat null on the others.
static void h2_extract_elements(
    F2Cptr* bmat, MPI_Comm caller_comm, int* Ninter, int* allrows, int* allcols,
    C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap,
    int* pmaps, const char* caller) {
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const std::vector<int> owner =
    single_rank_owners(*Ninter, pgidx, *Npmap, pmaps, (std::string(caller) + " (format 7)").c_str());
  butterfly::color_unstructured::extract_compressed_blocks(
    h2_of<H2Data>(bmat), caller_comm, *Ninter, allrows, allcols, rowidx, colidx, owner,
    reinterpret_cast<H2Data*>(alldat_loc));
}

// The sketched blocks of a list of requests from a compression-only
// unstructured Color H2 (as h2_extract_elements; Omega's entries 1-based).
static void h2_extract_sketches(
    F2Cptr* bmat, MPI_Comm caller_comm, int* Nsketch, int* allrows, int* allcols, int* rowidx, int* colidx,
    int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest, double* omega_val, C_DT* alldat_loc,
    int* pgidx, int* Npmap, int* pmaps, const char* caller, bool transposed) {
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  const std::vector<int> owner =
    single_rank_owners(*Nsketch, pgidx, *Npmap, pmaps, (std::string(caller) + " (format 7)").c_str());
  std::vector<int64_t> omega_off(static_cast<size_t>(std::max(*Nsketch, 0)) + 1, 0);
  for (int n = 0; n < *Nsketch; ++n) {
    if (omegaidx[n] < 0) throw std::invalid_argument(std::string(caller) + ": a negative sketch entry count");
    omega_off[static_cast<size_t>(n) + 1] = omega_off[static_cast<size_t>(n)] + omegaidx[n];
  }
  const int64_t nomega = omega_off.back();
  std::vector<int> orow(static_cast<size_t>(nomega)), odest(static_cast<size_t>(nomega));
  for (int64_t e = 0; e < nomega; ++e) {
    orow[static_cast<size_t>(e)] = omega_row[e] - 1;
    odest[static_cast<size_t>(e)] = omega_dest[e] - 1;
  }
  butterfly::color_unstructured::extract_compressed_sketches(
    h2_of<H2Data>(bmat), caller_comm, *Nsketch, allrows, allcols, rowidx, colidx, sketchdim, omega_off.data(),
    orow.data(), odest.data(), omega_val, owner, reinterpret_cast<H2Data*>(alldat_loc), transposed);
}
#endif

void c_bpack_extractelement(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,
	int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,int* allcols,
	C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps){
  // every rank holds the matrix: its own format
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    MPI_Comm comm = static_cast<butterfly::H2<double, H2Data>*>(H2_raw)->comm;
    try {
      h2_extract_elements(bmat, comm, Ninter, allrows, allcols, alldat_loc,
                          rowidx, colidx, pgidx, Npmap, pmaps, "c_bpack_extractelement");
    } catch (const std::exception& e) {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
      MPI_Abort(comm, 1);
    }
#else
	format7_requires_mpi("c_bpack_extractelement");
#endif
  }else{
	c_bpack_extractelement_fortran(bmat, option, msh, stats, ptree, Ninter, Nallrows, Nallcols,
	  Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps);
  }
}

void c_bpack_extractelement_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,
	F2Cptr* ptree_caller, int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc,
	int* allrows,int* allcols, C_DT* alldat_loc, int* rowidx, int* colidx, int* pgidx,
	int* Npmap, int* pmaps){
#ifdef HAVE_MPI
  // the source's format, from its ranks (the others hold no option)
  int fcomm = 0;
  c_bpack_get_comm(ptree_caller, &fcomm);
  MPI_Comm caller_comm = MPI_Comm_f2c(static_cast<MPI_Fint>(fcomm));
  int format = -1;
  if (option != nullptr && *option != nullptr) {
    double tmp = 0.0;
    c_bpack_getoption(option, "format", &tmp);
    format = static_cast<int>(tmp);
  }
  MPI_Allreduce(MPI_IN_PLACE, &format, 1, MPI_INT, MPI_MAX, caller_comm);
  if (format == 7) {
    try {
      h2_extract_elements(bmat, caller_comm, Ninter, allrows, allcols, alldat_loc,
                          rowidx, colidx, pgidx, Npmap, pmaps, "c_bpack_extractelement_bridge");
    } catch (const std::exception& e) {
      int rank = 0;
      MPI_Comm_rank(caller_comm, &rank);
      std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
      MPI_Abort(caller_comm, 1);
    }
    return;
  }
#endif
  c_bpack_extractelement_bridge_fortran(bmat, option, msh, stats, ptree, ptree_caller, Ninter,
    Nallrows, Nallcols, Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps);
}


void c_bpack_set_sketch_callback(F2Cptr* bmat, c_bpack_func_zmn_sketch C_FuncZmnSketch) {
#ifdef HAVE_MPI
  if (bmat == nullptr || *bmat == nullptr) {
    throw std::invalid_argument("c_bpack_set_sketch_callback: null matrix handle");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  butterfly::H2<double, H2Data>* solver = h2_of<H2Data>(bmat);
  if (solver == nullptr) return;  // (not an H2: its entries are never sketched)
  using Entries = fmm::CollectiveEntries<H2Data, typename butterfly::H2Kernel<double, H2Data>::kerData>;
  solver->kernel.sketch_callback = reinterpret_cast<typename Entries::SketchCallback>(C_FuncZmnSketch);
  if (solver->kernel.collective) solver->kernel.collective->set_sketch(solver->kernel.sketch_callback);
#else
  (void)bmat;
  (void)C_FuncZmnSketch;
#endif
}

void c_bpack_set_sketch_transpose_callback(F2Cptr* bmat, c_bpack_func_zmn_sketch C_FuncZmnSketchT) {
#ifdef HAVE_MPI
  if (bmat == nullptr || *bmat == nullptr) {
    throw std::invalid_argument("c_bpack_set_sketch_transpose_callback: null matrix handle");
  }
  using H2Data = typename butterfly::fmm_data<C_DT>::type;
  butterfly::H2<double, H2Data>* solver = h2_of<H2Data>(bmat);
  if (solver == nullptr) return;  // (not an H2: its entries are never sketched)
  using Entries = fmm::CollectiveEntries<H2Data, typename butterfly::H2Kernel<double, H2Data>::kerData>;
  solver->kernel.sketch_transpose_callback = reinterpret_cast<typename Entries::SketchCallback>(C_FuncZmnSketchT);
  if (solver->kernel.collective) {
    solver->kernel.collective->set_sketch_transposed(solver->kernel.sketch_transpose_callback);
  }
#else
  (void)bmat;
  (void)C_FuncZmnSketchT;
#endif
}

// c_bpack_extractsketch_bridge and its transposed form: Y_n = Omega_n
// A(rows, columns), or with transposed Omega_n A(columns, rows)^T.
static void extractsketch_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,
	F2Cptr* ptree_caller, int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc,
	int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row,
	int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps, bool transposed,
	const char* caller){
#ifdef HAVE_MPI
  (void)Nallomega;
  (void)Nalldat_loc;
  int fcomm = 0;
  c_bpack_get_comm(ptree_caller, &fcomm);
  MPI_Comm caller_comm = MPI_Comm_f2c(static_cast<MPI_Fint>(fcomm));
  int crank = 0;
  MPI_Comm_rank(caller_comm, &crank);
  // the source's format, from its ranks (the others hold no option)
  int format = -1;
  if (option != nullptr && *option != nullptr) {
    double tmp = 0.0;
    c_bpack_getoption(option, "format", &tmp);
    format = static_cast<int>(tmp);
  }
  MPI_Allreduce(MPI_IN_PLACE, &format, 1, MPI_INT, MPI_MAX, caller_comm);
  try {
    if (format == 7) {
      h2_extract_sketches(bmat, caller_comm, Nsketch, allrows, allcols, rowidx, colidx, sketchdim, omegaidx,
                          omega_row, omega_dest, omega_val, alldat_loc, pgidx, Npmap, pmaps, caller, transposed);
      return;
    }
    // other formats: the requests' entries (transposed: A(columns, rows)),
    // sketched on their owners
    using H2Data = typename butterfly::fmm_data<C_DT>::type;
    const std::vector<int> owner = single_rank_owners(*Nsketch, pgidx, *Npmap, pmaps, caller);
    int64_t nentries = 0;
    for (int n = 0; n < *Nsketch; ++n) {
      if (owner[static_cast<size_t>(n)] == crank) nentries += static_cast<int64_t>(rowidx[n]) * colidx[n];
    }
    std::vector<H2Data> entries(static_cast<size_t>(std::max<int64_t>(nentries, 1)));
    if (transposed) {
      c_bpack_extractelement_bridge_fortran(bmat, option, msh, stats, ptree, ptree_caller, Nsketch, Nallcols,
        Nallrows, &nentries, allcols, allrows, reinterpret_cast<C_DT*>(entries.data()), colidx, rowidx, pgidx,
        Npmap, pmaps);
    } else {
      c_bpack_extractelement_bridge_fortran(bmat, option, msh, stats, ptree, ptree_caller, Nsketch, Nallrows,
        Nallcols, &nentries, allrows, allcols, reinterpret_cast<C_DT*>(entries.data()), rowidx, colidx, pgidx,
        Npmap, pmaps);
    }
    H2Data* out = reinterpret_cast<H2Data*>(alldat_loc);
    int64_t at = 0, e0 = 0;
    for (int n = 0; n < *Nsketch; ++n) {
      const int64_t nr = rowidx[n], nc = colidx[n], d = sketchdim[n];
      if (owner[static_cast<size_t>(n)] == crank) {
        std::fill(out, out + d * nc, H2Data{0});
        const H2Data* A = entries.data() + at;  // nr x nc, or transposed nc x nr
        for (int64_t e = e0; e < e0 + omegaidx[n]; ++e) {
          const int64_t r = omega_row[e] - 1, i = omega_dest[e] - 1;
          const H2Data w = static_cast<H2Data>(omega_val[e]);
          if (transposed) {
            for (int64_t j = 0; j < nc; ++j) out[i + j * d] += w * A[j + r * nc];
          } else {
            for (int64_t j = 0; j < nc; ++j) out[i + j * d] += w * A[r + j * nr];
          }
        }
        at += nr * nc;
        out += d * nc;
      }
      e0 += omegaidx[n];
    }
  } catch (const std::exception& e) {
    std::cerr << "Error on rank " << crank << ": " << e.what() << std::endl;
    MPI_Abort(caller_comm, 1);
  }
#else
  (void)bmat; (void)option; (void)msh; (void)stats; (void)ptree; (void)ptree_caller; (void)Nsketch;
  (void)Nallrows; (void)Nallcols; (void)Nallomega; (void)Nalldat_loc; (void)allrows; (void)allcols;
  (void)rowidx; (void)colidx; (void)sketchdim; (void)omegaidx; (void)omega_row; (void)omega_dest;
  (void)omega_val; (void)alldat_loc; (void)pgidx; (void)Npmap; (void)pmaps; (void)transposed;
  throw std::runtime_error(std::string(caller) + " requires MPI");
#endif
}

void c_bpack_extractsketch_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,
	F2Cptr* ptree_caller, int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc,
	int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row,
	int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps){
  extractsketch_bridge(bmat, option, msh, stats, ptree, ptree_caller, Nsketch, Nallrows, Nallcols, Nallomega,
    Nalldat_loc, allrows, allcols, rowidx, colidx, sketchdim, omegaidx, omega_row, omega_dest, omega_val,
    alldat_loc, pgidx, Npmap, pmaps, false, "c_bpack_extractsketch_bridge");
}

void c_bpack_extractsketch_transpose_bridge(F2Cptr* bmat,F2Cptr* option,F2Cptr* msh,F2Cptr* stats,F2Cptr* ptree,
	F2Cptr* ptree_caller, int* Nsketch, int* Nallrows, int* Nallcols, int64_t* Nallomega, int64_t* Nalldat_loc,
	int* allrows, int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row,
	int* omega_dest, double* omega_val, C_DT* alldat_loc, int* pgidx, int* Npmap, int* pmaps){
  extractsketch_bridge(bmat, option, msh, stats, ptree, ptree_caller, Nsketch, Nallrows, Nallcols, Nallomega,
    Nalldat_loc, allrows, allcols, rowidx, colidx, sketchdim, omegaidx, omega_row, omega_dest, omega_val,
    alldat_loc, pgidx, Npmap, pmaps, true, "c_bpack_extractsketch_transpose_bridge");
}


void c_bpack_logdet(C_DT* phase, C_RDT* logabsdet, F2Cptr* option, F2Cptr* bmat){
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
#ifdef HAVE_MPI
	using H2Data = typename butterfly::fmm_data<C_DT>::type;
    void* H2_raw = nullptr;
    c_bpack_get_h2(*bmat, &H2_raw);
    butterfly::H2<double, H2Data>* H2_solver = static_cast<butterfly::H2<double, H2Data>*>(H2_raw);
	check_h2_symmetry_option(option, H2_solver->comm, "c_bpack_logdet");


	if (H2_solver->build_state != butterfly::H2BuildState::RS_FACTORIZED) {
	  throw std::runtime_error(
		"c_bpack_logdet (format 7): log-determinant requires an RS-S factorization");
	}
	double logabs_d = 0.0;
	butterfly::hierarchical_logdet_parallel(H2_solver->tree.get(), &logabs_d, reinterpret_cast<H2Data*>(phase));
	*logabsdet = static_cast<C_RDT>(logabs_d);
#else
	format7_requires_mpi("c_bpack_logdet");
#endif
  }else{
	c_bpack_logdet_fortran(phase, logabsdet, option, bmat);
  }
}

extern "C" void c_bpack_h2_delete(C2Fptr h2_ptr) {
#ifdef HAVE_MPI
	using H2Data = typename butterfly::fmm_data<C_DT>::type;
    auto* H2_solver = static_cast<butterfly::H2<double,H2Data>*>(h2_ptr);
#ifdef H2_HAVE_GPU
    // its device data (kept per operator, h2_gpu/device_heap.hpp)
    if (H2_solver != nullptr) fmm::gpu::release_operator(H2_solver->tree.get());
#endif
    delete H2_solver;
#else
	(void)h2_ptr;
#endif
}


// void c_bpack_delete(F2Cptr* option, F2Cptr*bmat) {
//   double tmp;
//   c_bpack_getoption(option, "format", &tmp);
//   int format=(int)tmp;
//   if(format==7){
// 	butterfly::H2<double, C_DT>* H2_solver = static_cast<butterfly::H2<double, C_DT>*>(*bmat);
// 	delete H2_solver;
// 	H2_solver = nullptr;
// 	*bmat = nullptr;
//   }else{
// 	c_bpack_delete_fortran(bmat);
//   }
// }
