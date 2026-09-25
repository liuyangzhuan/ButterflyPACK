#include "BPACK_wrapper.h"
#ifdef HAVE_MPI
#include "butterfly_integration.hpp"
#include <mpi.h>
#endif
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
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
static void require_symmetric_h2_option(
    F2Cptr* option,
    MPI_Comm comm,
    const char* caller) {
  double sym_d = 0.0;
  c_bpack_getoption(option, "sym", &sym_d);
  const int sym = static_cast<int>(std::llround(sym_d));
  if (sym == 1) return;

  const std::string message = std::string(caller) +
    ": ButterflyPACK H2 (format=7) currently requires option%sym=1 "
    "for a transpose-symmetric matrix; received option%sym=" +
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

  double compression_memory_mb =
    solver->factorization_memory / (1024.0 * 1024.0);
  c_bpack_setstats(stats, "Mem_Comp_for", &compression_memory_mb);

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
  double H2_GEMM_split_d = 0.0;
  double H2_CA_staged_halo_d = 0.0;
  double H2_CA_owner_component_d = 0.0;
  double H2_CA_owner_serial_d = 0.0;
  double H2_ID_radius_d = 0.0;
  double H2_ID_proxy_d = 0.0;
  double H2_ID_proxy_points_d = 0.0;
  double BACA_Batch_d = 0.0;
  c_bpack_getoption(option, "tol_comp", &tolerance);
  c_bpack_getoption(option, "reduction_threshold", &reduction_threshold_d);
  c_bpack_getoption(option, "Nmin_leaf", &Nmin_leaf_d);
  c_bpack_getoption(option, "precon", &precon_d);
  c_bpack_getoption(option, "verbosity", &verbosity_d);
  c_bpack_getoption(option, "H2_unstructured", &H2_unstructured_d);
  c_bpack_getoption(option, "CA_level", &CA_level_d);
  c_bpack_getoption(option, "H2_use_sketch", &H2_use_sketch_d);
  c_bpack_getoption(option, "H2_lazy_schur", &H2_lazy_schur_d);
  c_bpack_getoption(option, "H2_GEMM_split", &H2_GEMM_split_d);
  c_bpack_getoption(option, "H2_CA_staged_halo", &H2_CA_staged_halo_d);
  c_bpack_getoption(option, "H2_CA_owner_component", &H2_CA_owner_component_d);
  c_bpack_getoption(option, "H2_CA_owner_serial", &H2_CA_owner_serial_d);
  c_bpack_getoption(option, "H2_ID_radius", &H2_ID_radius_d);
  c_bpack_getoption(option, "H2_ID_proxy", &H2_ID_proxy_d);
  c_bpack_getoption(option, "H2_ID_proxy_points", &H2_ID_proxy_points_d);
  c_bpack_getoption(option, "BACA_Batch", &BACA_Batch_d);

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
  const int H2_GEMM_split =
      static_cast<int>(std::llround(H2_GEMM_split_d));
  const int H2_CA_staged_halo =
      static_cast<int>(std::llround(H2_CA_staged_halo_d));
  const int H2_CA_owner_component =
      static_cast<int>(std::llround(H2_CA_owner_component_d));
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
  if (H2_GEMM_split < 0) {
    throw std::invalid_argument("H2_GEMM_split must be nonnegative");
  }
  if (H2_CA_staged_halo != 0 && H2_CA_staged_halo != 2) {
    throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
  }
  if (H2_CA_owner_component != 0 && H2_CA_owner_component != 3) {
    throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
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
  if (H2_lazy_schur != 0 && H2_ID_proxy == 2) {
    throw std::invalid_argument(
        "H2_lazy_schur is not yet compatible with H2_ID_proxy=2");
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
  result.gemm_split = H2_GEMM_split;
  result.ca_staged_halo = H2_CA_staged_halo;
  result.ca_owner_component = H2_CA_owner_component;
  result.ca_owner_serial = H2_CA_owner_serial;
  result.id_neighborhood_radius = H2_ID_radius;
  result.id_proxy_mode = H2_ID_proxy;
  if (H2_ID_proxy == 1) {
    result.id_proxy_points = H2_ID_proxy_points;
  } else if (H2_ID_proxy == 2) {
    result.id_adaptive_batch = BACA_Batch;
  }
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
		{"sym",             "matrix symmetry flag; sym=1 is required by format-7 H2 and selects symmetric HODLR when format=1"},
		{"reduction_threshold", "format-7 H2 boxes-per-process threshold for reducing the active MPI process count"},
		{"h2_unstructured", "format-7 H2 backend: 0 structured Color/CA, 1 color_unstructured (Color only)"},
		{"h2_use_sketch",   "format-7 H2 ID mode: 0 full workspace, 1 materialized sparse sketch, 2 streamed sparse sketch on color levels"},
		{"h2_lazy_schur",   "format-7 H2 color Schur mode: 0 eager, 1 lazy far, 2 lazy far plus generated near"},
		{"h2_gemm_split",   "maximum task split for one format-7 H2 color work item; 0 disables splitting"},
		{"h2_ca_staged_halo", "format-7 H2 CA halo mode: 0 legacy, 2 staged/overlapped"},
		{"h2_ca_owner_component", "format-7 H2 CA ownership mode: 0 replicated, 3 asynchronous components"},
		{"h2_ca_owner_serial", "serialize format-7 H2 CA component-owner events (0 or 1)"},
		{"h2_id_radius",    "format-7 H2 mandatory ID neighborhood radius; 2 keeps the standard workspace"},
		{"h2_id_proxy",     "format-7 H2 proxy mode: 0 none, 1 geometric surface, 2 adaptive row sampling"},
		{"h2_id_proxy_points", "geometric surface samples when h2_id_proxy=1"},
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
		{"h2_ca_owner_serial", required_argument, 0, 56},
		{"H2_CA_owner_serial", required_argument, 0, 56},
		{"h2_unstructured", required_argument, 0, 57},
		{"H2_unstructured", required_argument, 0, 57},
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
	require_symmetric_h2_option(option, mpi_comm, "c_bpack_construct_init");

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
      double H2_GEMM_split_d;
      double H2_CA_staged_halo_d;
      double H2_CA_owner_component_d;
      double H2_CA_owner_serial_d;
      double H2_ID_radius_d;
      double H2_ID_proxy_d;
      double H2_ID_proxy_points_d;
      double BACA_Batch_d;
      c_bpack_getoption(option, "tol_comp", &tolerance);
      c_bpack_getoption(option, "reduction_threshold", &reduction_threshold_d);
      c_bpack_getoption(option, "Nmin_leaf", &Nmin_leaf_d);
      c_bpack_getoption(option, "precon", &precon_d);
      c_bpack_getoption(option, "verbosity", &verbosity_d);
      c_bpack_getoption(option, "H2_unstructured", &H2_unstructured_d);
      c_bpack_getoption(option, "CA_level", &CA_level_d);
      c_bpack_getoption(option, "H2_use_sketch", &H2_use_sketch_d);
      c_bpack_getoption(option, "H2_lazy_schur", &H2_lazy_schur_d);
      c_bpack_getoption(option, "H2_GEMM_split", &H2_GEMM_split_d);
      c_bpack_getoption(option, "H2_CA_staged_halo", &H2_CA_staged_halo_d);
      c_bpack_getoption(option, "H2_CA_owner_component", &H2_CA_owner_component_d);
      c_bpack_getoption(option, "H2_CA_owner_serial", &H2_CA_owner_serial_d);
      c_bpack_getoption(option, "H2_ID_radius", &H2_ID_radius_d);
      c_bpack_getoption(option, "H2_ID_proxy", &H2_ID_proxy_d);
      c_bpack_getoption(option, "H2_ID_proxy_points", &H2_ID_proxy_points_d);
      c_bpack_getoption(option, "BACA_Batch", &BACA_Batch_d);
      int64_t reduction_threshold = (int64_t)reduction_threshold_d;
      int64_t Nmin_leaf = (int64_t)Nmin_leaf_d;
      const int CA_level = static_cast<int>(std::llround(CA_level_d));
      const int H2_unstructured =
          static_cast<int>(std::llround(H2_unstructured_d));
      const int H2_use_sketch = static_cast<int>(std::llround(H2_use_sketch_d));
      const int H2_lazy_schur = static_cast<int>(std::llround(H2_lazy_schur_d));
      const int H2_GEMM_split = static_cast<int>(std::llround(H2_GEMM_split_d));
      const int H2_CA_staged_halo =
          static_cast<int>(std::llround(H2_CA_staged_halo_d));
      const int H2_CA_owner_component =
          static_cast<int>(std::llround(H2_CA_owner_component_d));
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
      if (H2_GEMM_split < 0) {
        throw std::invalid_argument("H2_GEMM_split must be nonnegative");
      }
      if (H2_CA_staged_halo != 0 && H2_CA_staged_halo != 2) {
        throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
      }
      if (H2_CA_owner_component != 0 && H2_CA_owner_component != 3) {
        throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
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
      if (H2_lazy_schur != 0 && H2_ID_proxy == 2) {
        throw std::invalid_argument(
            "H2_lazy_schur is not yet compatible with H2_ID_proxy=2");
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
      H2_options.use_sketch = H2_use_sketch;
      H2_options.lazy_schur = H2_lazy_schur;
      H2_options.gemm_split = H2_GEMM_split;
      H2_options.ca_staged_halo = H2_CA_staged_halo;
      H2_options.ca_owner_component = H2_CA_owner_component;
      H2_options.ca_owner_serial = H2_CA_owner_serial;
      H2_options.id_neighborhood_radius = H2_ID_radius;
      H2_options.id_proxy_mode = H2_ID_proxy;
      if (H2_ID_proxy == 1) {
        H2_options.id_proxy_points = H2_ID_proxy_points;
      } else if (H2_ID_proxy == 2) {
        H2_options.id_adaptive_batch = BACA_Batch;
      }
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
                  << ", h2_ca_staged_halo=" << H2_options.ca_staged_halo
                  << ", h2_ca_owner_component=" << H2_options.ca_owner_component
                  << ", h2_ca_owner_serial=" << H2_options.ca_owner_serial
                  << ", h2_id_radius=" << H2_options.id_neighborhood_radius
                  << ", h2_id_proxy=" << H2_options.id_proxy_mode;
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
	  // convert to perms and Npo_loc
	  c_bpack_set_mesh_h2(Npo, new2old.data(), &idxs, &idxe, msh);
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
    require_symmetric_h2_option(
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
  require_symmetric_h2_option(
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
	require_symmetric_h2_option(
	  option, H2_solver->comm, "c_bpack_construct_element_compute");
	H2_solver->kernel.kernel = C_FuncZmn;
	H2_solver->kernel.kernel64 = nullptr;
	H2_solver->kernel.block_kernel64 = nullptr;
	H2_solver->kernel.dimension = H2_solver->options.dimension;
	double elem_extract_d = 0.0;
	c_bpack_getoption(option, "elem_extract", &elem_extract_d);
	const int elem_extract = static_cast<int>(std::llround(elem_extract_d));
	H2_solver->kernel.block_kernel =
	  (elem_extract == 2) ? C_FuncZmnBlock : nullptr;
	MPI_Comm_rank(MPI_COMM_WORLD, &H2_solver->kernel.block_callback_pid);
	H2_solver->kernel.quant = C_QuantApp;
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
	require_symmetric_h2_option(option, H2_solver->comm, "c_bpack_factor");

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

		double factorization_memory_MB = H2_solver->factorization_memory/(1024.0 * 1024.0);
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
	require_symmetric_h2_option(option, H2_solver->comm, "c_bpack_solve");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    try {
      const int verbosity = sync_h2_verbosity(option, H2_solver);
      double t0 = MPI_Wtime();
      const H2Data* b_h2 = reinterpret_cast<const H2Data*>(b);
	  // A distributed64 handle exposes the caller's input ownership at this
	  // public boundary.  Convert the entire RHS batch to tree ownership before
	  // solving and convert the solution back below.  Legacy handles have no
	  // layout object and already use the solver's internal ownership.
	  const auto* distributed_layout = H2_solver->distributed_layout.get();
	  int internal_nloc = *Nloc;
	  if (distributed_layout != nullptr) {
		if (distributed_layout->input_local_size != *Nloc) {
		  throw std::invalid_argument(
			"c_bpack_solve (distributed64): Nloc must equal the input-local size");
		}
		if (distributed_layout->internal_local_size >
			std::numeric_limits<int>::max()) {
		  throw std::overflow_error(
			"c_bpack_solve (distributed64): internal local size exceeds INT_MAX");
		}
		internal_nloc = static_cast<int>(
			distributed_layout->internal_local_size);
	  }
	  std::vector<H2Data> rhs(
		static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Nrhs));
	  if (distributed_layout != nullptr) {
		distributed_layout->input_to_internal(
		  b_h2, rhs.data(), *Nrhs);
	  } else {
		std::copy_n(b_h2, rhs.size(), rhs.data());
	  }
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
		std::vector<std::vector<fmm::SolveDataRequest<double, H2Data>>> solve_data(
		  H2_solver->options.num_levels);
		butterfly::dispatch_h2_solve(
		  H2_solver, rhs, solve_data, *Nrhs, verbosity);
		internal_solution.resize(
		  static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Nrhs));
		butterfly::gather_local_solution(
		  H2_solver->tree.get(), solve_data,
		  internal_solution.data(), &internal_nloc, *Nrhs);
	  }

	  if (internal_solution.size() !=
		  static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Nrhs)) {
		throw std::runtime_error(
		  "c_bpack_solve (format 7): internal solution length mismatch");
	  }
	  if (distributed_layout != nullptr) {
		distributed_layout->internal_to_input(
		  internal_solution.data(), reinterpret_cast<H2Data*>(x), *Nrhs);
	  } else {
		std::copy(internal_solution.begin(), internal_solution.end(),
		          reinterpret_cast<H2Data*>(x));
	  }

	  double t_solve = MPI_Wtime() - t0;
	  MPI_Allreduce(MPI_IN_PLACE, &t_solve, 1, MPI_DOUBLE, MPI_MAX, H2_solver->comm);
	  c_bpack_setstats(stats, "Time_Solve", &t_solve);
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
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
	require_symmetric_h2_option(option, H2_solver->comm, "c_bpack_mult");

    int rank = 0;
    MPI_Comm_rank(H2_solver->comm, &rank);

    if (*Ncol <= 0) {
      throw std::invalid_argument(
        "c_bpack_mult (format 7): Ncol must be positive");
    }

    // Only F·x is implemented for format 7 (transpose/conj-transpose not yet supported).
    const char op = (trans && trans[0]) ? trans[0] : 'N';
    if (op != 'N' && op != 'n') {
      throw std::runtime_error(
          "c_bpack_mult (format 7): only trans == 'N' is currently supported; got '" +
          std::string(1, op) + "'");
    }

    try {
      const int verbosity = sync_h2_verbosity(option, H2_solver);
      double t0 = MPI_Wtime();
      double precon_d = 1.0;
	  c_bpack_getoption(option, "precon", &precon_d);
      H2_solver->options.precon = static_cast<int>(std::llround(precon_d));

      if (H2_solver->build_state == butterfly::H2BuildState::UNBUILT) {
		throw std::runtime_error(
		  "c_bpack_mult (format 7): call c_bpack_factor before multiplication");
      }

      const H2Data* xin_h2 = reinterpret_cast<const H2Data*>(xin);
	  // Keep the public multiply contract in caller input order for a
	  // distributed64 handle; the H2 multiply itself operates in tree order.
	  // Legacy handles bypass both redistributions.
	  const auto* distributed_layout = H2_solver->distributed_layout.get();
	  int internal_nloc = *Ninloc;
	  if (distributed_layout != nullptr) {
		if (distributed_layout->input_local_size != *Ninloc ||
			distributed_layout->input_local_size != *Noutloc) {
		  throw std::invalid_argument(
			"c_bpack_mult (distributed64): Ninloc and Noutloc must equal "
			"the input-local size");
		}
		if (distributed_layout->internal_local_size >
			std::numeric_limits<int>::max()) {
		  throw std::overflow_error(
			"c_bpack_mult (distributed64): internal local size exceeds INT_MAX");
		}
		internal_nloc = static_cast<int>(
			distributed_layout->internal_local_size);
	  }
	  std::vector<H2Data> lhs(
		static_cast<size_t>(internal_nloc) * static_cast<size_t>(*Ncol));
	  if (distributed_layout != nullptr) {
		distributed_layout->input_to_internal(xin_h2, lhs.data(), *Ncol);
	  } else {
		std::copy_n(xin_h2, lhs.size(), lhs.data());
	  }
	  std::vector<H2Data> compressed_output;
	  std::vector<std::vector<fmm::SolveDataRequest<double, H2Data>>> mul_data;

	  if (H2_solver->build_state == butterfly::H2BuildState::H2_COMPRESSED) {
		butterfly::dispatch_h2_compressed_multiply(
		  H2_solver, lhs, compressed_output,
		  *Ncol, verbosity >= 1);
	  } else {
		mul_data.resize(H2_solver->options.num_levels);
		butterfly::dispatch_h2_multiply(
		  H2_solver, lhs, mul_data,
		  *Ncol, verbosity >= 1);
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
	  if (distributed_layout != nullptr) {
		distributed_layout->internal_to_input(
		  internal_output.data(), reinterpret_cast<H2Data*>(xout), *Ncol);
	  } else {
		std::copy(internal_output.begin(), internal_output.end(),
		          reinterpret_cast<H2Data*>(xout));
	  }

	  double t_mult = MPI_Wtime() - t0;
	  MPI_Allreduce(MPI_IN_PLACE, &t_mult, 1, MPI_DOUBLE, MPI_MAX,
	                H2_solver->comm);

	  double prev = 0.0;
	  c_bpack_getstats(stats, "Time_C_Mult_Wrapper", &prev);
	  double total = prev + t_mult;
	  c_bpack_setstats(stats, "Time_C_Mult_Wrapper", &total);

    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#else
	format7_requires_mpi("c_bpack_mult");
#endif
  }else{
	c_bpack_mult_fortran(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  }
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
	require_symmetric_h2_option(option, H2_solver->comm, "c_bpack_logdet");


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
    delete static_cast<butterfly::H2<double,H2Data>*>(h2_ptr);
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
