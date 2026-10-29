#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>
#include <unistd.h>

#include "dBPACK_wrapper.h"
#include "bpack_inv_mult_check.hpp"
#ifdef BPACK_EXAMPLE_GPU
#include "gpu/gpu_evaluators.h"
#endif

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

struct DriverOptions {
  int64_t grid_size = 96;
  int64_t expected_points = -1;
  int num_levels = 0;
  int64_t nmin_leaf = 216;
  bool nmin_leaf_set = false;
  double tolerance = 1e-3;
  int64_t reduction_threshold = 8;
  int h2_unstructured = 0;
  int ca_level = 0;
  int h2_use_sketch = 1;
  int h2_lazy_schur = 0;
  int h2_compress_factors = 0;
  int h2_gemm_split = 16;
  int h2_xrr_factor = 0;
  int h2_use_gpu = 0;
  int hodlr_use_gpu = 0;
  int reclr_leaf = -1;  // -1: the library default
  int hodlr_gpu_pieces = -1;  // -1: the library default
  int baca_batch = -1;
  int h2_id_radius = -1;  // -1: the library default
  int h2_id_proxy = -1;
  int h2_id_proxy_points = -1;
  double h2_s_op = 0.0;  // H2_s_op: the IDs' absolute floor tol_comp * h2_s_op (0 off)
  int h2_ca_staged_halo = 0;
  int h2_ca_owner_component = 0;
  int h2_ca_pack_shared = 0;  // 0 auto, 1 always, 2 never
  int h2_ca_solve = 0;        // 0 replicated, 1 home
  int h2_ca_owner_serial = 0;
  int precon = 1;
  int nrhs = 1;
  int inv_mult_check = 0;  // columns of the inv_mult / transposed mult check (bpack_inv_mult_check.hpp); 0: none
  int extract_check = 0;   // intersections of the extraction check of a compression-only H2 (precon 2); 0: none
  int64_t h2_partial = 0;  // H2_partial: the first h2_partial points (the grid's first planes) form F11
  int h2_partial_extract = 0;  // H2_partial_extract: 0 auto, 1 device, 2 host
  int h2_extract = 0;      // H2_extract: 0 auto, 1 device, 2 host
  int partial_check = 0;   // 1: the partial solves and Schur complement against dense LAPACK on rank 0 (N <= 40000, distributed64 1)
  int partial_extract_bench = 0;  // > 0: that many rounds of parent-like Schur complement blocks through c_bpack_extractelement
  int logdet_check = 0;    // 1: the logdet and a solve again after another H2 is built (format 7, not precon 2)
  int verbosity = 1;
  int elem_extract = -1;  // -1: 2, with hodlr_extract_ranks 1
  int lrlevel = 0;
  int hodlr_extract_ranks = 0;
  int gpu_evaluator = 1;
  int format = 7;
  int sym = 1;
  // a nonsymmetric matrix (--sym 0; H2: --H2_unstructured 1): a dipole term
  // c (d . (x - y)) / r along a fixed direction d, and column weights w(y) =
  // 1 + a sin(2 pi (y0 + 2 y1 + 3 y2)) (exact checks: log det(A W) =
  // log det(A) + sum log w, x = W^{-1} x_sym)
  double skew = 0.0;
  double column_weight = 0.0;
  int distributed64 = 0;
  bool show_help = false;
};

std::string normalize_option_name(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    if (c == '-') return '_';
    return static_cast<char>(std::tolower(c));
  });
  return name;
}

int64_t parse_int64(const std::string& text, const char* name) {
  size_t used = 0;
  long long value = 0;
  try {
    value = std::stoll(text, &used);
  } catch (const std::exception&) {
    throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
  }
  if (used != text.size()) {
    throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
  }
  return static_cast<int64_t>(value);
}

int parse_int(const std::string& text, const char* name) {
  const int64_t value = parse_int64(text, name);
  if (value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::out_of_range(std::string(name) + " is outside the int range");
  }
  return static_cast<int>(value);
}

double parse_double(const std::string& text, const char* name) {
  size_t used = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &used);
  } catch (const std::exception&) {
    throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
  }
  if (used != text.size() || !std::isfinite(value)) {
    throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
  }
  return value;
}

int64_t checked_cube(int64_t value, const char* name) {
  if (value <= 0 || value > std::numeric_limits<int64_t>::max() / value ||
      value * value > std::numeric_limits<int64_t>::max() / value) {
    throw std::overflow_error(std::string(name) + " cubed overflows int64_t");
  }
  return value * value * value;
}

void derive_nmin_leaf_from_levels(DriverOptions& options) {
  if (options.num_levels == 0) return;
  if (options.num_levels < 2 || options.num_levels >= 63) {
    throw std::invalid_argument("num_levels must be between 2 and 62");
  }

  const int64_t boxes_per_dimension = int64_t{1} << (options.num_levels - 1);
  if (options.grid_size % boxes_per_dimension != 0) {
    throw std::invalid_argument(
        "grid_size must be divisible by 2^(num_levels-1)");
  }

  const int64_t leaf_grid_size = options.grid_size / boxes_per_dimension;
  const int64_t derived_nmin_leaf = checked_cube(leaf_grid_size, "leaf grid size");
  if (options.nmin_leaf_set && options.nmin_leaf != derived_nmin_leaf) {
    throw std::invalid_argument(
        "Nmin_leaf conflicts with the requested num_levels and grid_size");
  }
  options.nmin_leaf = derived_nmin_leaf;
}

DriverOptions parse_driver_options(int argc, char** argv) {
  DriverOptions options;
  int first_option = 1;

  // Accept the standalone driver's positional form:
  //   <num_levels> <N> <grid_size> <tolerance>
  if (argc >= 5 && argv[1][0] != '-') {
    options.num_levels = parse_int(argv[1], "num_levels");
    options.expected_points = parse_int64(argv[2], "N");
    options.grid_size = parse_int64(argv[3], "grid_size");
    options.tolerance = parse_double(argv[4], "tolerance");
    first_option = 5;
  }

  for (int i = first_option; i < argc; ++i) {
    std::string argument = argv[i];
    if (argument == "--help" || argument == "-h") {
      options.show_help = true;
      continue;
    }
    if (argument.rfind("--", 0) != 0) {
      throw std::invalid_argument("unexpected positional argument: " + argument);
    }

    argument.erase(0, 2);
    const size_t equals = argument.find('=');
    const std::string name = normalize_option_name(argument.substr(0, equals));
    std::string value;
    if (equals != std::string::npos) {
      value = argument.substr(equals + 1);
    }

    if (value.empty()) {
      if (++i >= argc) {
        throw std::invalid_argument("missing value after --" + name);
      }
      value = argv[i];
    }

    if (name == "grid_size" || name == "grid") {
      options.grid_size = parse_int64(value, "grid_size");
    } else if (name == "n" || name == "points") {
      options.expected_points = parse_int64(value, "N");
    } else if (name == "num_levels" || name == "levels") {
      options.num_levels = parse_int(value, "num_levels");
    } else if (name == "tol_comp" || name == "tolerance") {
      options.tolerance = parse_double(value, "tol_comp");
    } else if (name == "nmin_leaf") {
      options.nmin_leaf = parse_int64(value, "Nmin_leaf");
      options.nmin_leaf_set = true;
    } else if (name == "reduction_threshold") {
      options.reduction_threshold = parse_int64(value, "reduction_threshold");
    } else if (name == "h2_unstructured") {
      options.h2_unstructured = parse_int(value, "H2_unstructured");
    } else if (name == "ca_level") {
      options.ca_level = parse_int(value, "CA_level");
    } else if (name == "h2_use_sketch") {
      options.h2_use_sketch = parse_int(value, "H2_use_sketch");
    } else if (name == "h2_lazy_schur") {
      options.h2_lazy_schur = parse_int(value, "H2_lazy_schur");
    } else if (name == "h2_compress_factors") {
      options.h2_compress_factors = parse_int(value, "H2_compress_factors");
    } else if (name == "h2_gemm_split") {
      options.h2_gemm_split = parse_int(value, "H2_GEMM_split");
    } else if (name == "h2_xrr_factor") {
      options.h2_xrr_factor = parse_int(value, "H2_XRR_factor");
    } else if (name == "h2_use_gpu") {
      options.h2_use_gpu = parse_int(value, "H2_use_gpu");
    } else if (name == "hodlr_use_gpu") {
      options.hodlr_use_gpu = parse_int(value, "HODLR_use_gpu");
    } else if (name == "reclr_leaf") {
      options.reclr_leaf = parse_int(value, "RecLR_leaf");
    } else if (name == "hodlr_gpu_pieces") {
      options.hodlr_gpu_pieces = parse_int(value, "HODLR_gpu_pieces");
    } else if (name == "baca_batch") {
      options.baca_batch = parse_int(value, "BACA_Batch");
    } else if (name == "h2_id_radius") {
      options.h2_id_radius = parse_int(value, "H2_ID_radius");
    } else if (name == "h2_id_proxy") {
      options.h2_id_proxy = parse_int(value, "H2_ID_proxy");
    } else if (name == "h2_id_proxy_points") {
      options.h2_id_proxy_points = parse_int(value, "H2_ID_proxy_points");
    } else if (name == "h2_s_op") {
      options.h2_s_op = parse_double(value, "H2_s_op");
    } else if (name == "h2_ca_staged_halo") {
      options.h2_ca_staged_halo = parse_int(value, "H2_CA_staged_halo");
    } else if (name == "h2_ca_owner_component") {
      options.h2_ca_owner_component =
          parse_int(value, "H2_CA_owner_component");
    } else if (name == "h2_ca_pack_shared") {
      options.h2_ca_pack_shared =
          value == "auto" ? 0 : value == "always" ? 1 : value == "never" ? 2
                                                    : parse_int(value, "H2_CA_pack_shared");
    } else if (name == "h2_ca_solve") {
      options.h2_ca_solve =
          value == "replicated" ? 0 : value == "home" ? 1 : parse_int(value, "H2_CA_solve");
    } else if (name == "h2_ca_owner_serial") {
      options.h2_ca_owner_serial = parse_int(value, "H2_CA_owner_serial");
    } else if (name == "precon") {
      options.precon = parse_int(value, "precon");
    } else if (name == "nrhs") {
      options.nrhs = parse_int(value, "nrhs");
    } else if (name == "inv_mult_check") {
      options.inv_mult_check = parse_int(value, "inv_mult_check");
    } else if (name == "h2_partial") {
      options.h2_partial = parse_int64(value, "H2_partial");
    } else if (name == "h2_partial_extract") {
      options.h2_partial_extract = value == "auto"     ? 0
                                   : value == "device" ? 1
                                   : value == "host"   ? 2
                                                       : parse_int(value, "H2_partial_extract");
    } else if (name == "h2_extract") {
      options.h2_extract = value == "auto"     ? 0
                           : value == "device" ? 1
                           : value == "host"   ? 2
                                               : parse_int(value, "H2_extract");
    } else if (name == "partial_check") {
      options.partial_check = parse_int(value, "partial_check");
    } else if (name == "partial_extract_bench") {
      options.partial_extract_bench = parse_int(value, "partial_extract_bench");
    } else if (name == "extract_check") {
      options.extract_check = parse_int(value, "extract_check");
    } else if (name == "logdet_check") {
      options.logdet_check = parse_int(value, "logdet_check");
    } else if (name == "verbosity") {
      options.verbosity = parse_int(value, "verbosity");
    } else if (name == "elem_extract") {
      options.elem_extract = parse_int(value, "elem_extract");
    } else if (name == "hodlr_extract_ranks" || name == "extract_ranks") {
      options.hodlr_extract_ranks = parse_int(value, "extract_ranks");
    } else if (name == "gpu_evaluator") {
      options.gpu_evaluator = parse_int(value, "gpu_evaluator");
    } else if (name == "lrlevel") {
      options.lrlevel = parse_int(value, "lrlevel");
    } else if (name == "distributed64") {
      options.distributed64 = parse_int(value, "distributed64");
    } else if (name == "kernel") {
      if (normalize_option_name(value) != "laplace") {
        throw std::invalid_argument("this driver supports only --kernel laplace");
      }
    } else if (name == "number_type" || name == "number") {
      if (normalize_option_name(value) != "real") {
        throw std::invalid_argument("this driver supports only --number-type real");
      }
    } else if (name == "dimension") {
      if (parse_int(value, "dimension") != 3) {
        throw std::invalid_argument("this driver supports only --dimension 3");
      }
    } else if (name == "num_proxy") {
      if (parse_int(value, "num_proxy") != 0) {
        throw std::invalid_argument(
            "the ButterflyPACK H2 interface currently supports num_proxy=0 only");
      }
    } else if (name == "format") {
      options.format = parse_int(value, "format");
    } else if (name == "sym") {
      options.sym = parse_int(value, "sym");
    } else if (name == "skew") {
      options.skew = std::stod(value);
    } else if (name == "column_weight") {
      options.column_weight = std::stod(value);
    } else {
      throw std::invalid_argument("unknown option: --" + name);
    }
  }

  if (options.grid_size < 2) {
    throw std::invalid_argument("grid_size must be at least 2");
  }
  if (!(options.tolerance > 0.0)) {
    throw std::invalid_argument("tol_comp must be positive");
  }
  if (options.nmin_leaf <= 0 || options.reduction_threshold <= 0) {
    throw std::invalid_argument("Nmin_leaf and reduction_threshold must be positive");
  }
  if (options.hodlr_extract_ranks < 0) {
    throw std::invalid_argument("hodlr_extract_ranks must be nonnegative");
  }
  if (options.hodlr_extract_ranks > 0 &&
      ((options.format != 1 && options.format != 7) || options.distributed64 != 0)) {
    throw std::invalid_argument(
        "extract_ranks requires format=1 or 7 and distributed64=0");
  }
  if (options.hodlr_extract_ranks > 0 && options.format == 7 &&
      (options.elem_extract == 0 || options.h2_unstructured != 1 || options.precon == 2)) {
    // (the H2's blocks from the compressed H2 through the collective block
    // callback; unstructured Color only; the H2 under test factored)
    throw std::invalid_argument(
        "extract_ranks with format=7 requires elem_extract 1, H2_unstructured 1 and precon 1 or 3");
  }
  if (options.gpu_evaluator != 0 && options.gpu_evaluator != 1) {
    throw std::invalid_argument("gpu_evaluator must be 0 or 1");
  }
  if (options.elem_extract == -1) {
    options.elem_extract = options.hodlr_extract_ranks > 0 ? 1 : 2;
  }
  if (options.elem_extract < 0 || options.elem_extract > 2) {
    throw std::invalid_argument("elem_extract must be 0, 1, or 2");
  }
  if (options.hodlr_extract_ranks > 0 && options.elem_extract == 2) {
    // (the extraction from the source is collective over all ranks)
    throw std::invalid_argument(
        "with extract_ranks, elem_extract must be 1 (or 0)");
  }
  if (options.precon < 1 || options.precon > 3) {
    throw std::invalid_argument("precon must be 1, 2, or 3");
  }
  if (options.nrhs <= 0) {
    throw std::invalid_argument("nrhs must be positive");
  }
  if (options.h2_use_sketch < 0 || options.h2_use_sketch > 2) {
    throw std::invalid_argument("H2_use_sketch must be 0, 1, or 2");
  }
  if (options.h2_unstructured != 0 && options.h2_unstructured != 1) {
    throw std::invalid_argument("H2_unstructured must be 0 or 1");
  }
  if (options.h2_unstructured == 1 && options.h2_lazy_schur == 1) {
    throw std::invalid_argument(
        "color_unstructured supports H2_lazy_schur=0 or 2 only");
  }
  if (options.h2_compress_factors < 0 || options.h2_compress_factors > 1) {
    throw std::invalid_argument("H2_compress_factors must be 0 (off) or 1 (on); 2 is not a separate mode: 1 compresses the fill sources too");
  }
  if (options.h2_lazy_schur < 0 || options.h2_lazy_schur > 2) {
    throw std::invalid_argument("H2_lazy_schur must be 0, 1, or 2");
  }
  if (options.h2_lazy_schur != 0 && options.h2_use_sketch != 2) {
    throw std::invalid_argument("H2_lazy_schur requires H2_use_sketch=2");
  }
  if (options.h2_gemm_split < 0) {
    throw std::invalid_argument("H2_GEMM_split must be nonnegative");
  }
  if (options.h2_xrr_factor != 0 && options.h2_xrr_factor != 1) {
    throw std::invalid_argument("H2_XRR_factor must be 0 or 1");
  }
  if (options.h2_use_gpu < 0 || options.h2_use_gpu > 2) {
    throw std::invalid_argument("H2_use_gpu must be 0 (CPU) or 1 (GPU, FP64 tensor-core GEMMs); 2 is taken as 1");
  }
  if (options.hodlr_use_gpu < 0 || options.hodlr_use_gpu > 2) {
    throw std::invalid_argument("HODLR_use_gpu must be 0, 1 or 2");
  }
  if (options.hodlr_gpu_pieces != -1 && options.hodlr_gpu_pieces < 1) {
    throw std::invalid_argument("HODLR_gpu_pieces must be at least 1");
  }
  if (options.h2_ca_staged_halo != 0 && options.h2_ca_staged_halo != 2) {
    throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
  }
  if (options.h2_ca_owner_component != 0 &&
      options.h2_ca_owner_component != 3) {
    throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
  }
  if (options.h2_ca_pack_shared < 0 || options.h2_ca_pack_shared > 2) {
    throw std::invalid_argument("H2_CA_pack_shared must be auto (0), always (1) or never (2)");
  }
  if (options.h2_partial_extract < 0 || options.h2_partial_extract > 2) {
    throw std::invalid_argument("H2_partial_extract must be auto (0), device (1) or host (2)");
  }
  if (options.h2_extract < 0 || options.h2_extract > 2) {
    throw std::invalid_argument("H2_extract must be auto (0), device (1) or host (2)");
  }
  if (options.h2_ca_solve < 0 || options.h2_ca_solve > 1) {
    throw std::invalid_argument("H2_CA_solve must be replicated (0) or home (1)");
  }
  if (options.h2_ca_owner_serial != 0 && options.h2_ca_owner_serial != 1) {
    throw std::invalid_argument("H2_CA_owner_serial must be 0 or 1");
  }
  if (options.format != 1 && options.format != 7) {
    throw std::invalid_argument("format must be 1 (HODLR) or 7 (H2)");
  }
  if (options.sym != 0 && options.sym != 1) {
    throw std::invalid_argument("sym must be 0 or 1");
  }
  if (options.format == 7 && options.sym != 1 && options.h2_unstructured != 1) {
    throw std::invalid_argument("format=7 (H2) with --sym 0 requires --H2_unstructured 1");
  }
  if ((options.skew != 0.0 || options.column_weight != 0.0) && options.sym != 0) {
    throw std::invalid_argument("--skew and --column_weight make the matrix nonsymmetric: pass --sym 0");
  }
  if (!(options.column_weight >= 0.0 && options.column_weight < 1.0)) {
    throw std::invalid_argument("--column_weight must lie in [0, 1)");
  }
  if (options.distributed64 != 0 && options.distributed64 != 1) {
    throw std::invalid_argument("distributed64 must be 0 or 1");
  }
  if (options.distributed64 == 1 && options.format != 7) {
    throw std::invalid_argument(
        "the distributed64 API currently supports only format=7");
  }

  derive_nmin_leaf_from_levels(options);
  const int64_t points = checked_cube(options.grid_size, "grid_size");
  if (options.expected_points > 0 && options.expected_points != points) {
    throw std::invalid_argument(
        "N must equal grid_size^3; expected " + std::to_string(points));
  }
  if (options.distributed64 == 0 &&
      points > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(
        "the legacy ButterflyPACK C API requires N <= INT_MAX; use "
        "--distributed64 1 for format=7");
  }
  if (options.nmin_leaf > points) {
    throw std::invalid_argument("Nmin_leaf cannot exceed N");
  }
  if (options.nmin_leaf > std::numeric_limits<int>::max() ||
      options.reduction_threshold > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(
        "Nmin_leaf and reduction_threshold must fit in an int");
  }

  return options;
}

void print_usage(const char* executable) {
  std::cout
      << "Usage:\n"
      << "  " << executable
      << " [<num_levels> <N> <grid_size> <tolerance>] [options]\n\n"
      << "Defaults reproduce the active jscript_laplace_large_CA.sh case:\n"
      << "  num_levels=5, N=884736, grid_size=96, tolerance=1e-3\n"
      << "  Nmin_leaf=216, reduction_threshold=8, CA_level=0\n\n"
      << "Options:\n"
      << "  --grid-size <n>\n"
      << "  --num-levels <count>\n"
      << "  --tol-comp <value>\n"
      << "  --Nmin_leaf <count>\n"
      << "  --reduction_threshold <count>\n"
      << "  --H2_unstructured <0|1>\n"
      << "  --CA_level <level>\n"
      << "  --H2_use_sketch <0|1|2>\n"
      << "  --H2_lazy_schur <0|1|2>\n"
      << "  --H2_compress_factors <0|1>  (GPU: lossy compression at tol_comp of\n"
         "                                the kept solve factors and of the\n"
         "                                factorization's fill sources and packed\n"
         "                                blocks; the factorization is lossy too)\n"
      << "  --H2_GEMM_split <count>\n"
      << "  --H2_XRR_factor <0|1>\n"
      << "  --H2_use_gpu <0|1>   (1: GPU, FP64 tensor-core GEMMs; 2 is taken as 1)\n"
      << "  --HODLR_use_gpu <0|1|2>\n"
      << "  --HODLR_gpu_pieces <n>   (pieces per rank of a shared HODLR block on the GPU; default 4)\n"
      << "  --RecLR_leaf <n>   (low-rank compression of HODLR blocks; 5: BACA without overlap)\n"
      << "  --BACA_Batch <n>   (also the rows sampled per spatial node with --H2_ID_proxy 2)\n"
      << "  --H2_ID_radius <n>   (ID neighborhood radius; default 2)\n"
      << "  --H2_ID_proxy <0|1|2>   (0: none, 1: geometric surface, 2: adaptive row sampling; needs --distributed64 0)\n"
      << "  --H2_ID_proxy_points <n>   (surface samples of --H2_ID_proxy 1)\n"
      << "  --H2_s_op <s>   (IDs: absolute floor tol_comp*s on the rank rule; at most ||A|| for --precon 2, at most 1/||A^-1|| factored; 0 off)\n"
      << "  --H2_CA_staged_halo <0|2>\n"
      << "  --H2_CA_owner_component <0|3>\n"
      << "  --H2_CA_pack_shared <auto|always|never>\n"
      << "  --H2_CA_solve <replicated|home>\n"
      << "  --H2_CA_owner_serial <0|1>\n"
      << "  --precon <1|2|3>\n"
      << "  --nrhs <count>\n"
      << "  --inv_mult_check <columns>   (after the solve: inv_mult N/T/C and the\n"
         "                                transposed mult; with precon 2 the\n"
         "                                transposed mult alone)\n"
      << "  --H2_partial <n>   (format=7, H2_unstructured 1, sym 0: eliminate the first n points only;\n"
         "                            the Schur complement on the rest is served by c_bpack_partial_*)\n"
      << "  --H2_partial_extract <auto|device|host>   (with --H2_partial: where the Schur complement's\n"
         "                            entries are computed; auto: the device when the records are on it\n"
         "                            and the call is large enough)\n"
      << "  --partial_extract_bench <R>   (with --H2_partial: R rounds of Schur complement blocks as a\n"
      << "                            parent front would read them through c_bpack_extractelement - near\n"
      << "                            field 64 x 64 blocks and 16 x 4096 panels - timed)\n"
      << "  --partial_check <0|1>   (with --H2_partial, N <= 40000: the partial solves\n"
         "                            and the Schur complement against dense LAPACK on rank 0)\n"
      << "  --extract_check <count>   (format=7, precon 2, distributed64 1: that many\n"
         "                             random blocks extracted from the compressed H2,\n"
         "                             against the kernel and its matvec)\n"
      << "  --H2_extract <auto|device|host>   (format=7, precon 2: where a compressed H2's entries\n"
         "                            are computed; auto: the device when the compression kept its\n"
         "                            blocks there and the call is large enough)\n"
      << "  --logdet_check <0|1>   (format=7, precon 1 or 3: the logdet and a solve of\n"
         "                          the factored H2 again after a compression-only H2\n"
         "                          of the same kernel is built)\n"
      << "  --format <1|7>\n"
      << "  --sym <0|1>          (0: nonsymmetric; format=7 needs --H2_unstructured 1)\n"
      << "  --skew <c>           (a nonsymmetric dipole term c (d . (x - y)) / r; --sym 0)\n"
      << "  --column_weight <a>  (columns times 1 + a sin(...), 0 <= a < 1; --sym 0)\n"
      << "  --distributed64 <0|1>\n"
      << "  --elem_extract <0|1|2>   (default 2; with --hodlr_extract_ranks: of\n"
         "                            the HODLR from the entries, default 1)\n"
      << "  --extract_ranks <P>   (format=1: build the HODLR from the kernel on\n"
         "                            the first P ranks, then the one factored and\n"
         "                            solved on all ranks from its entries; format=7:\n"
         "                            the same with a compression-only H2 (precon 2)\n"
         "                            and an H2 factored from its entries through the\n"
         "                            collective block callback, elem_extract 1;\n"
         "                            alias --hodlr_extract_ranks)\n"
      << "  --gpu_evaluator <0|1>   (0: no GPU evaluator of the kernel; the GPU\n"
         "                            HODLR construction takes the entries through\n"
         "                            the block callback)\n"
      << "  --verbosity <-1|0|1>\n";
}

double laplace_self_cell_integral(int64_t grid_size) {
  const double h = 1.0 / static_cast<double>(grid_size);
  constexpr std::array<double, 5> nodes = {
      -0.9061798459386640, -0.5384693101056831, 0.0,
       0.5384693101056831,  0.9061798459386640};
  constexpr std::array<double, 5> weights = {
      0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
      0.4786286704993665, 0.2369268850561891};

  const double transform = h / 4.0;
  const double shift = h / 4.0;
  double integral = 0.0;
  for (size_t i = 0; i < nodes.size(); ++i) {
    const double x = transform * nodes[i] + shift;
    for (size_t j = 0; j < nodes.size(); ++j) {
      const double y = transform * nodes[j] + shift;
      for (size_t k = 0; k < nodes.size(); ++k) {
        const double z = transform * nodes[k] + shift;
        const double r = std::sqrt(x * x + y * y + z * z);
        integral += weights[i] * weights[j] * weights[k] /
                    (4.0 * kPi * r);
      }
    }
  }

  return 8.0 * integral * transform * transform * transform;
}

class Laplace3DApplication {
 public:
  Laplace3DApplication(int64_t grid_size, bool materialize_locations, double skew = 0.0,
                       double column_weight = 0.0)
      : grid_size_(grid_size),
        point_count_(checked_cube(grid_size, "grid_size")),
        inverse_4pi_n_(1.0 / (4.0 * kPi * static_cast<double>(point_count_))),
        diagonal_(laplace_self_cell_integral(grid_size)),
        skew_(skew),
        column_weight_(column_weight) {
    if (materialize_locations) {
      locations_.resize(static_cast<size_t>(point_count_) * 3);
      for (int64_t index = 0; index < point_count_; ++index) {
        coordinate(index, locations_.data() + static_cast<size_t>(3 * index));
      }
    }
  }

  double* locations() { return locations_.data(); }
  double diagonal() const { return diagonal_; }
  double inverse_4pi_n() const { return inverse_4pi_n_; }
  double skew() const { return skew_; }
  double column_weight() const { return column_weight_; }
  const double* skew_direction() const { return kSkewDirection; }
  bool symmetric() const { return skew_ == 0.0 && column_weight_ == 0.0; }

  // w(y) of --column_weight
  double weight(const double* y) const {
    if (column_weight_ == 0.0) return 1.0;
    return 1.0 + column_weight_ * std::sin(2.0 * kPi * (y[0] + 2.0 * y[1] + 3.0 * y[2]));
  }
  // sum over the points of log w (log det(A W) - log det(A))
  double log_weight_sum() const {
    double sum = 0.0;
    for (int64_t index = 0; index < point_count_; ++index) {
      double y[3];
      coordinate(index, y);
      sum += std::log(weight(y));
    }
    return sum;
  }

  void coordinate(int64_t index, double* point) const {
    const int64_t k = index % grid_size_;
    const int64_t quotient = index / grid_size_;
    const int64_t j = quotient % grid_size_;
    const int64_t i = quotient / grid_size_;
    const double h = 1.0 / static_cast<double>(grid_size_);
    point[0] = (static_cast<double>(i) + 0.5) * h;
    point[1] = (static_cast<double>(j) + 0.5) * h;
    point[2] = (static_cast<double>(k) + 0.5) * h;
  }

  double evaluate_coordinates(int64_t row_id, int64_t column_id,
                              const double* x, const double* y) const {
    if (row_id == column_id) return diagonal_ * weight(y);
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    double value = inverse_4pi_n_ / r;
    if (skew_ != 0.0) {
      value *= 1.0 + skew_ * (kSkewDirection[0] * dx + kSkewDirection[1] * dy + kSkewDirection[2] * dz) / r;
    }
    return value * weight(y);
  }

  double evaluate(int row, int column) const {
    const double* x = locations_.data() + static_cast<size_t>(3 * row);
    const double* y = locations_.data() + static_cast<size_t>(3 * column);
    return evaluate_coordinates(row + 1, column + 1, x, y);
  }

  void evaluate_block(int rows, int columns, const int* row_indices,
                      const int* column_indices, double* output) const {
    for (int j = 0; j < columns; ++j) {
      const int column = column_indices[j] - 1;
      const double* y = locations_.data() + static_cast<size_t>(3 * column);
      for (int i = 0; i < rows; ++i) {
        const int row = row_indices[i] - 1;
        const double* x = locations_.data() + static_cast<size_t>(3 * row);
        output[i + static_cast<int64_t>(j) * rows] = evaluate_coordinates(
            row + 1, column + 1, x, y);
      }
    }
  }

 private:
  // the dipole direction of --skew: (1, 2, 3) / sqrt(14)
  static constexpr double kSkewDirection[3] = {0.2672612419124244, 0.5345224838248488, 0.8017837257372732};
  int64_t grid_size_;
  int64_t point_count_;
  double inverse_4pi_n_;
  double diagonal_;
  double skew_;
  double column_weight_;
  std::vector<double> locations_;
};

void laplace_entry_callback(int* row, int* column, double* value, C2Fptr quant) {
  const auto* application = static_cast<Laplace3DApplication*>(quant);
  *value = application->evaluate(*row - 1, *column - 1);
}

void laplace_block_callback(
    int* ninter, int* nallrows, int* nallcols, int64_t* nalldat_loc,
    int* allrows, int* allcols, double* alldat_loc, int* rowidx,
    int* colidx, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  const auto* application = static_cast<Laplace3DApplication*>(quant);
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int64_t row_offset = 0;
  int64_t column_offset = 0;
  int64_t value_offset = 0;
  for (int interaction = 0; interaction < *ninter; ++interaction) {
    const int process_group = pgidx[interaction];
    const int nprow = pmaps[process_group];
    const int npcol = pmaps[*npmap + process_group];
    const int owner = pmaps[2 * (*npmap) + process_group];
    const int rows = rowidx[interaction];
    const int columns = colidx[interaction];

    if (nprow * npcol != 1) {
      if (world_rank == 0) {
        std::cerr << "laplace_block_callback supports only single-process blocks"
                  << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 2);
    }

    if (*npmap == 1 || world_rank == owner) {
      application->evaluate_block(
          rows, columns, allrows + row_offset, allcols + column_offset,
          alldat_loc + value_offset);
      value_offset += static_cast<int64_t>(rows) * columns;
    }

    row_offset += rows;
    column_offset += columns;
  }

  (void)nallrows;
  (void)nallcols;
  (void)nalldat_loc;
}

void laplace_entry_callback64(
    const int64_t* row, const int64_t* column, const double* row_coordinate,
    const double* column_coordinate, const int* dimension, double* value,
    C2Fptr quant) {
  const auto* application = static_cast<Laplace3DApplication*>(quant);
  if (*dimension != 3) {
    throw std::invalid_argument("laplace_entry_callback64 requires dimension=3");
  }
  *value = application->evaluate_coordinates(
      *row, *column, row_coordinate, column_coordinate);
}

void laplace_block_callback64(
    const int64_t* rows, const int64_t* columns, const int64_t* row_ids,
    const int64_t* column_ids, const double* row_coordinates,
    const double* column_coordinates, const int* dimension, double* output,
    const int64_t* leading_dimension, C2Fptr quant) {
  const auto* application = static_cast<Laplace3DApplication*>(quant);
  if (*dimension != 3 || *leading_dimension < *rows) {
    throw std::invalid_argument("invalid Laplace distributed block metadata");
  }
  for (int64_t column = 0; column < *columns; ++column) {
    const double* y = column_coordinates + 3 * column;
    for (int64_t row = 0; row < *rows; ++row) {
      const double* x = row_coordinates + 3 * row;
      output[row + column * *leading_dimension] =
          application->evaluate_coordinates(
              row_ids[row], column_ids[column], x, y);
    }
  }
}

struct DistributedInput {
  int64_t first = 0;
  int64_t count = 0;
  std::vector<int64_t> global_ids;
  std::vector<double> coordinates;
};

DistributedInput make_distributed_input(
    int64_t global_count, int rank, int mpi_size,
    const Laplace3DApplication& application) {
  DistributedInput input;
  const int64_t base = global_count / mpi_size;
  const int64_t remainder = global_count % mpi_size;
  input.count = base + (rank < remainder ? 1 : 0);
  input.first = static_cast<int64_t>(rank) * base +
      std::min<int64_t>(rank, remainder);
  input.global_ids.resize(static_cast<size_t>(input.count));
  input.coordinates.resize(static_cast<size_t>(input.count) * 3);
  for (int64_t local = 0; local < input.count; ++local) {
    const int64_t global = input.first + local;
    input.global_ids[static_cast<size_t>(local)] = global + 1;
    application.coordinate(
        global, input.coordinates.data() + static_cast<size_t>(3 * local));
  }
  return input;
}

void distance_callback(int*, int*, double* value, C2Fptr) { *value = 0.0; }
void near_far_callback(int*, int*, int* value, C2Fptr) { *value = 0; }

// The source of the entries of the HODLR under test with
// --hodlr_extract_ranks: the HODLR built from the kernel on the first ranks
// (its handles null on the others), and the process tree of the HODLR under
// test, over all ranks
struct ExtractionSource {
  F2Cptr matrix = nullptr;
  F2Cptr option = nullptr;
  F2Cptr mesh = nullptr;
  F2Cptr stats = nullptr;
  F2Cptr process_tree = nullptr;
  F2Cptr caller_tree = nullptr;
};

void extraction_entry_callback(int*, int*, double*, C2Fptr) {
  std::cerr << "the HODLR from another one's entries takes them by blocks "
               "(elem_extract=1)"
            << std::endl;
  MPI_Abort(MPI_COMM_WORLD, 2);
}

// The entries of the HODLR under test (its intersections, in its natural
// order: the source's tree order) from the source, collective over all ranks
// (elem_extract=1)
void extraction_block_callback(
    int* ninter, int* nallrows, int* nallcols, int64_t* nalldat_loc,
    int* allrows, int* allcols, double* alldat_loc, int* rowidx,
    int* colidx, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  auto* source = static_cast<ExtractionSource*>(quant);
  d_c_bpack_extractelement_bridge(
      &source->matrix, &source->option, &source->mesh, &source->stats,
      &source->process_tree, &source->caller_tree, ninter, nallrows,
      nallcols, nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx,
      pgidx, npmap, pmaps);
}

// The sketched blocks of the H2 under test (the ID's sketches of its kernel
// rows) from the source, collective over all ranks (elem_extract=1)
void extraction_sketch_callback(
    int* nsketch, int* nallrows, int* nallcols, int64_t* nallomega, int64_t* nalldat_loc, int* allrows,
    int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest,
    double* omega_val, double* alldat_loc, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  auto* source = static_cast<ExtractionSource*>(quant);
  d_c_bpack_extractsketch_bridge(
      &source->matrix, &source->option, &source->mesh, &source->stats, &source->process_tree,
      &source->caller_tree, nsketch, nallrows, nallcols, nallomega, nalldat_loc, allrows, allcols, rowidx,
      colidx, sketchdim, omegaidx, omega_row, omega_dest, omega_val, alldat_loc, pgidx, npmap, pmaps);
}

// The transposed sketched blocks, Omega A(columns, rows)^T (a nonsymmetric
// H2's ID), from the source, collective over all ranks (elem_extract=1)
void extraction_sketch_transpose_callback(
    int* nsketch, int* nallrows, int* nallcols, int64_t* nallomega, int64_t* nalldat_loc, int* allrows,
    int* allcols, int* rowidx, int* colidx, int* sketchdim, int* omegaidx, int* omega_row, int* omega_dest,
    double* omega_val, double* alldat_loc, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  auto* source = static_cast<ExtractionSource*>(quant);
  d_c_bpack_extractsketch_transpose_bridge(
      &source->matrix, &source->option, &source->mesh, &source->stats, &source->process_tree,
      &source->caller_tree, nsketch, nallrows, nallcols, nallomega, nalldat_loc, allrows, allcols, rowidx,
      colidx, sketchdim, omegaidx, omega_row, omega_dest, omega_val, alldat_loc, pgidx, npmap, pmaps);
}

struct ButterflyResources {
  F2Cptr option = nullptr;
  F2Cptr stats = nullptr;
  F2Cptr process_tree = nullptr;
  F2Cptr matrix = nullptr;
  F2Cptr mesh = nullptr;

  ~ButterflyResources() {
    if (stats != nullptr) d_c_bpack_deletestats(&stats);
    if (mesh != nullptr) d_c_bpack_deletemesh(&mesh);
    if (matrix != nullptr) d_c_bpack_delete(&matrix);
    if (process_tree != nullptr) d_c_bpack_deleteproctree(&process_tree);
    if (option != nullptr) d_c_bpack_deleteoption(&option);
  }
};

}  // namespace

extern "C" void Cblacs_exit(int);

extern "C" {
void dgetrf_(const int*, const int*, double*, const int*, int*, int*);
double dlange_(const char*, const int*, const int*, const double*, const int*, double*);
void dgecon_(const char*, const int*, const double*, const int*, const double*, double*, double*, int*, int*);
void dgetrs_(const char*, const int*, const int*, const double*, const int*, const int*, double*, const int*, int*);
void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*,
            const double*, const int*, const double*, double*, const int*);
}

// --H2_partial: the partial solves and the Schur complement CB = F22 - F21
// F11^-1 F12 through the C API, and with --partial_check against dense
// LAPACK on rank 0 (the kernel over all points).  A local row's input point
// is input.global_ids (1-based): with distributed64 1 the rows
// [first, first + count) in input order, else the rank's range of the tree
// order.
void partial_checks(const DriverOptions& o, const Laplace3DApplication& application, ButterflyResources& resources,
                    int rank, int mpi_size, int local_points, const DistributedInput& input) {
  const int64_t N = checked_cube(o.grid_size, "grid_size");
  const int n11 = static_cast<int>(o.h2_partial), n22 = static_cast<int>(N - o.h2_partial);
  const int nrhs = o.nrhs;
  auto global_of = [&](int r) { return input.global_ids[static_cast<size_t>(r)] - 1; };
  if (n22 == 0) {
    if (rank == 0) std::cout << "H2 partial: no F22 points (H2_partial = N), nothing more to check" << std::endl;
    return;
  }
  // b: a function of the global row (so rank 0 can regenerate it), x2 random
  std::vector<double> b(static_cast<size_t>(local_points) * nrhs);
  for (int c = 0; c < nrhs; ++c)
    for (int r = 0; r < local_points; ++r)
      b[static_cast<size_t>(r) + static_cast<size_t>(c) * local_points] =
          1.0 + 0.125 * c + 0.001 * (global_of(r) % 97);
  std::mt19937_64 rng(20261022);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  std::vector<double> x2(static_cast<size_t>(n22) * nrhs);
  for (double& v : x2) v = uni(rng);
  std::vector<double> out2(x2.size(), 0.0), x(b.size(), 0.0), y2(x2.size(), 0.0);
  int nloc = local_points, N22 = n22, Nrhs = nrhs;
  const double t0 = MPI_Wtime();
  d_c_bpack_partial_forward(out2.data(), b.data(), &nloc, &N22, &Nrhs, &resources.matrix, &resources.option,
                            &resources.stats, &resources.process_tree);
  const double t1 = MPI_Wtime();
  d_c_bpack_partial_backward(x.data(), x2.data(), &nloc, &N22, &Nrhs, &resources.matrix, &resources.option,
                             &resources.stats, &resources.process_tree);
  const double t2 = MPI_Wtime();
  d_c_bpack_partial_schur_mult(y2.data(), x2.data(), &N22, &Nrhs, &resources.matrix, &resources.option);
  const double t3 = MPI_Wtime();
  // the F22 entries of x come back as given
  double x2_diff = 0.0;
  for (int c = 0; c < nrhs; ++c)
    for (int r = 0; r < local_points; ++r) {
      const int64_t g = global_of(r);
      if (g >= n11) x2_diff = std::max(x2_diff, std::abs(x[static_cast<size_t>(r) + static_cast<size_t>(c) * local_points] -
                                                        x2[static_cast<size_t>(g - n11) + static_cast<size_t>(c) * n22]));
    }
  MPI_Allreduce(MPI_IN_PLACE, &x2_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  // three random blocks of the Schur complement, and one over a run of
  // consecutive F22 indices (neighbouring points: the near-field entries)
  constexpr int nblocks = 4;
  std::vector<int> rows, cols;
  std::vector<int> nr(nblocks), nc(nblocks);
  for (int k = 0; k < 3; ++k) {
    nr[k] = 1 + static_cast<int>(rng() % 40);
    nc[k] = 1 + static_cast<int>(rng() % 40);
    for (int i = 0; i < nr[k]; ++i) rows.push_back(static_cast<int>(rng() % static_cast<uint64_t>(n22)));
    for (int j = 0; j < nc[k]; ++j) cols.push_back(static_cast<int>(rng() % static_cast<uint64_t>(n22)));
  }
  nr[3] = nc[3] = std::min(40, n22 - n22 / 3);
  for (int i = 0; i < nr[3]; ++i) rows.push_back(n22 / 3 + i);
  for (int j = 0; j < nc[3]; ++j) cols.push_back(n22 / 3 + j);
  std::vector<std::vector<double>> blocks(nblocks);
  {
    int ra = 0, ca = 0;
    for (int k = 0; k < nblocks; ++k) {
      blocks[k].assign(static_cast<size_t>(nr[k]) * nc[k], 0.0);
      d_c_bpack_partial_schur_extract(blocks[k].data(), rows.data() + ra, &nr[k], cols.data() + ca, &nc[k],
                                      &resources.matrix, &resources.option);
      ra += nr[k];
      ca += nc[k];
    }
  }
  const double t4 = MPI_Wtime();
  // the same blocks through c_bpack_extractelement, as STRUMPACK reads a
  // child's update block: one call, 1-based F22 indices, block k to rank k
  // mod P (1 x 1 grids)
  double ext_diff = 0.0, ext_scale = 0.0;
  {
    std::vector<int> allrows, allcols, rowidx(nblocks), colidx(nblocks), pgidx(nblocks), pmaps(3 * nblocks);
    int64_t nloc_dat = 0;
    int ra = 0, ca = 0;
    for (int k = 0; k < nblocks; ++k) {
      for (int i = 0; i < nr[k]; ++i) allrows.push_back(rows[static_cast<size_t>(ra + i)] + 1);
      for (int j = 0; j < nc[k]; ++j) allcols.push_back(cols[static_cast<size_t>(ca + j)] + 1);
      rowidx[k] = nr[k];
      colidx[k] = nc[k];
      pgidx[k] = k;
      pmaps[k] = 1;
      pmaps[nblocks + k] = 1;
      pmaps[2 * nblocks + k] = k % mpi_size;
      if (k % mpi_size == rank) nloc_dat += static_cast<int64_t>(nr[k]) * nc[k];
      ra += nr[k];
      ca += nc[k];
    }
    std::vector<double> dat(static_cast<size_t>(std::max<int64_t>(nloc_dat, 1)));
    int ninter = nblocks, nallrows = static_cast<int>(allrows.size()), nallcols = static_cast<int>(allcols.size()), npmap = nblocks;
    d_c_bpack_extractelement(&resources.matrix, &resources.option, &resources.mesh, &resources.stats, &resources.process_tree,
                             &ninter, &nallrows, &nallcols, &nloc_dat, allrows.data(), allcols.data(), dat.data(),
                             rowidx.data(), colidx.data(), pgidx.data(), &npmap, pmaps.data());
    int64_t at = 0;
    for (int k = 0; k < nblocks; ++k) {
      if (k % mpi_size != rank) continue;
      for (size_t e = 0; e < blocks[k].size(); ++e) {
        ext_diff = std::max(ext_diff, std::abs(dat[static_cast<size_t>(at) + e] - blocks[k][e]));
        ext_scale = std::max(ext_scale, std::abs(blocks[k][e]));
      }
      at += static_cast<int64_t>(blocks[k].size());
    }
    MPI_Allreduce(MPI_IN_PLACE, &ext_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ext_scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  const double t5 = MPI_Wtime();
  double times[5] = {t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4};
  MPI_Allreduce(MPI_IN_PLACE, times, 5, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << std::setprecision(3) << "H2 partial solves: forward " << times[0] << " s, backward " << times[1]
              << " s, Schur mult " << times[2] << " s, Schur extract " << times[3] << " s; x2 returned to "
              << x2_diff << std::endl;
    std::cout << std::setprecision(3) << "H2 partial extractelement: the same " << nblocks << " blocks in one call "
              << times[4] << " s, against the Schur extract's " << ext_diff / std::max(ext_scale, 1e-300)
              << " (relative)" << std::endl;
  }
  if (o.partial_check == 0) return;
  if (N > 40000) throw std::invalid_argument("--partial_check: N <= 40000");
  // rank 0: the dense kernel, A11 = P L U, W = A11^-1 A12, CB = A22 - A21 W
  const int n = static_cast<int>(N);
  std::vector<double> x1d(static_cast<size_t>(n11) * nrhs, 0.0);
  double err_mult = 0.0, err_extract = 0.0, err_forward = 0.0;
  // kept for the backward solve's residual: A11 before its LU, A12, b1; A11's 1-norm condition estimate
  std::vector<double> A11k, A12k, b1k;
  double cond11 = 0.0;
  if (rank == 0) {
    std::vector<double> A(static_cast<size_t>(n) * n);
    for (int j = 0; j < n; ++j) {
      double y[3];
      application.coordinate(j, y);
      for (int i = 0; i < n; ++i) {
        double xx[3];
        application.coordinate(i, xx);
        A[static_cast<size_t>(i) + static_cast<size_t>(j) * n] = application.evaluate_coordinates(i + 1, j + 1, xx, y);
      }
    }
    auto block = [&](int r0, int nrow, int c0, int ncol) {
      std::vector<double> B(static_cast<size_t>(nrow) * ncol);
      for (int j = 0; j < ncol; ++j)
        for (int i = 0; i < nrow; ++i)
          B[static_cast<size_t>(i) + static_cast<size_t>(j) * nrow] = A[static_cast<size_t>(r0 + i) + static_cast<size_t>(c0 + j) * n];
      return B;
    };
    std::vector<double> A11 = block(0, n11, 0, n11), A12 = block(0, n11, n11, n22), A21 = block(n11, n22, 0, n11),
                        CB = block(n11, n22, n11, n22);
    std::vector<int> piv(static_cast<size_t>(std::max(n11, 1)));
    int info = 0;
    A11k = A11;
    A12k = A12;
    const char one_norm = '1';
    std::vector<double> lwork(static_cast<size_t>(4) * std::max(n11, 1));
    std::vector<int> liwork(static_cast<size_t>(std::max(n11, 1)));
    const double anorm = dlange_(&one_norm, &n11, &n11, A11.data(), &n11, lwork.data());
    dgetrf_(&n11, &n11, A11.data(), &n11, piv.data(), &info);
    if (info != 0) throw std::runtime_error("partial_check: dense LU of A11 failed");
    double rcond = 0.0;
    dgecon_(&one_norm, &n11, A11.data(), &n11, &anorm, &rcond, lwork.data(), liwork.data(), &info);
    cond11 = rcond > 0.0 ? 1.0 / rcond : 0.0;
    const char no = 'N';
    std::vector<double> W = A12;
    dgetrs_(&no, &n11, &n22, A11.data(), &n11, piv.data(), W.data(), &n11, &info);
    const double one = 1.0, minus = -1.0, zero = 0.0;
    dgemm_(&no, &no, &n22, &n22, &n11, &minus, A21.data(), &n22, W.data(), &n11, &one, CB.data(), &n22);
    // y2 = CB x2
    std::vector<double> y2d(static_cast<size_t>(n22) * nrhs, 0.0);
    dgemm_(&no, &no, &n22, &nrhs, &n22, &one, CB.data(), &n22, x2.data(), &n22, &zero, y2d.data(), &n22);
    double d2 = 0.0, r2 = 0.0;
    for (size_t i = 0; i < y2d.size(); ++i) {
      d2 += (y2[i] - y2d[i]) * (y2[i] - y2d[i]);
      r2 += y2d[i] * y2d[i];
    }
    err_mult = std::sqrt(d2 / r2);
    // the blocks
    {
      int ra = 0, ca = 0;
      double bd2 = 0.0, br2 = 0.0;
      for (int k = 0; k < nblocks; ++k) {
        for (int j = 0; j < nc[k]; ++j)
          for (int i = 0; i < nr[k]; ++i) {
            const double ref = CB[static_cast<size_t>(rows[static_cast<size_t>(ra + i)]) +
                                  static_cast<size_t>(cols[static_cast<size_t>(ca + j)]) * n22];
            const double v = blocks[k][static_cast<size_t>(i) + static_cast<size_t>(j) * nr[k]];
            bd2 += (v - ref) * (v - ref);
            br2 += ref * ref;
          }
        ra += nr[k];
        ca += nc[k];
      }
      err_extract = std::sqrt(bd2 / br2);
    }
    // the forward solve: out2 = b2 - A21 A11^-1 b1; the backward: x1 = A11^-1 (b1 - A12 x2)
    std::vector<double> b1(static_cast<size_t>(n11) * nrhs), b2(static_cast<size_t>(n22) * nrhs);
    for (int c = 0; c < nrhs; ++c) {
      for (int r = 0; r < n11; ++r) b1[static_cast<size_t>(r) + static_cast<size_t>(c) * n11] = 1.0 + 0.125 * c + 0.001 * (r % 97);
      for (int r = 0; r < n22; ++r) b2[static_cast<size_t>(r) + static_cast<size_t>(c) * n22] = 1.0 + 0.125 * c + 0.001 * ((r + n11) % 97);
    }
    std::vector<double> y1 = b1;
    dgetrs_(&no, &n11, &nrhs, A11.data(), &n11, piv.data(), y1.data(), &n11, &info);
    dgemm_(&no, &no, &n22, &nrhs, &n11, &minus, A21.data(), &n22, y1.data(), &n11, &one, b2.data(), &n22);
    d2 = r2 = 0.0;
    for (size_t i = 0; i < b2.size(); ++i) {
      d2 += (out2[i] - b2[i]) * (out2[i] - b2[i]);
      r2 += b2[i] * b2[i];
    }
    err_forward = std::sqrt(d2 / r2);
    x1d = b1;
    b1k = b1;
    dgemm_(&no, &no, &n11, &nrhs, &n22, &minus, A12.data(), &n11, x2.data(), &n22, &one, x1d.data(), &n11);
    dgetrs_(&no, &n11, &nrhs, A11.data(), &n11, piv.data(), x1d.data(), &n11, &info);
  }
  MPI_Bcast(x1d.data(), static_cast<int>(x1d.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double d2 = 0.0, r2 = 0.0;
  for (int c = 0; c < nrhs; ++c)
    for (int r = 0; r < local_points; ++r) {
      const int64_t g = global_of(r);
      if (g >= n11) continue;
      const double v = x[static_cast<size_t>(r) + static_cast<size_t>(c) * local_points];
      const double ref = x1d[static_cast<size_t>(g) + static_cast<size_t>(c) * n11];
      d2 += (v - ref) * (v - ref);
      r2 += ref * ref;
    }
  MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &r2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double err_backward = std::sqrt(d2 / r2);
  // the backward solve's residual against the dense blocks, ||b1 - A12 x2 - A11 x1|| / ||b1 - A12 x2||:
  // its error is at most about cond(A11) times this
  std::vector<double> x1h(static_cast<size_t>(n11) * nrhs, 0.0);
  for (int c = 0; c < nrhs; ++c)
    for (int r = 0; r < local_points; ++r) {
      const int64_t g = global_of(r);
      if (g < n11) x1h[static_cast<size_t>(g) + static_cast<size_t>(c) * n11] = x[static_cast<size_t>(r) + static_cast<size_t>(c) * local_points];
    }
  MPI_Allreduce(MPI_IN_PLACE, x1h.data(), static_cast<int>(x1h.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0) {
    const char no = 'N';
    const double one = 1.0, minus = -1.0;
    std::vector<double> rhs = b1k;
    dgemm_(&no, &no, &n11, &nrhs, &n22, &minus, A12k.data(), &n11, x2.data(), &n22, &one, rhs.data(), &n11);
    std::vector<double> res = rhs;
    dgemm_(&no, &no, &n11, &nrhs, &n11, &minus, A11k.data(), &n11, x1h.data(), &n11, &one, res.data(), &n11);
    double rr = 0.0, bb = 0.0;
    for (size_t i = 0; i < res.size(); ++i) {
      rr += res[i] * res[i];
      bb += rhs[i] * rhs[i];
    }
    const double residual = std::sqrt(rr / bb);
    std::cout << std::scientific << std::setprecision(3) << "H2 partial dense check (N " << n << ", F11 " << n11
              << "): schur_mult " << err_mult << "  schur_extract " << err_extract << "  forward " << err_forward
              << "  backward " << err_backward << " (residual " << residual << ", cond1(A11) " << cond11 << ")"
              << std::defaultfloat << std::endl;
  }
  (void)mpi_size;
}

// This process's resident set (GB, /proc/self/statm).
double resident_gb() {
  unsigned long pages = 0, resident = 0;
  if (FILE* f = std::fopen("/proc/self/statm", "r")) {
    if (std::fscanf(f, "%lu %lu", &pages, &resident) != 2) resident = 0;
    std::fclose(f);
  }
  return static_cast<double>(resident) * static_cast<double>(sysconf(_SC_PAGESIZE)) / 1e9;
}

// --partial_extract_bench R: R rounds of Schur complement blocks as a
// STRUMPACK parent front reads a child's update block (extract_2d through
// c_bpack_extractelement, 1-based F22 indices, each block to one rank): per
// round 64 dense blocks between neighbouring 4^3 cubes of the F22 points (a
// compressed front's dense leaves) and 32 panels, 16 points of a cube against
// a 16^3 region not touching it and the transpose (BACA's row and column
// samples).  Prints the time a round and the rate.
void partial_extract_bench(const DriverOptions& o, ButterflyResources& resources, int rank, int mpi_size) {
  const int64_t n = o.grid_size, n22 = n * n * n - o.h2_partial;
  if (n22 <= 0) return;
  const int64_t i0 = (o.h2_partial + n * n - 1) / (n * n);  // the first whole F22 plane
  if (n - i0 < 16 || n < 16) throw std::invalid_argument("--partial_extract_bench: at least 16^3 whole F22 points");
  auto f22 = [&](int64_t i, int64_t j, int64_t k) { return static_cast<int>(k + n * j + n * n * i - o.h2_partial + 1); };
  auto cube = [&](int64_t ci, int64_t cj, int64_t ck, int64_t side, std::vector<int>& out) {
    for (int64_t i = ci; i < ci + side; ++i)
      for (int64_t j = cj; j < cj + side; ++j)
        for (int64_t k = ck; k < ck + side; ++k) out.push_back(f22(i, j, k));
  };
  std::mt19937_64 rng(20261025);
  auto pick = [&](int64_t lo, int64_t hi) { return lo + static_cast<int64_t>(rng() % static_cast<uint64_t>(hi - lo + 1)); };
  double t_total = 0.0, entries = 0.0;
  for (int r = 0; r < o.partial_extract_bench; ++r) {
    std::vector<int> allrows, allcols, rowidx, colidx;
    auto add = [&](const std::vector<int>& rr, const std::vector<int>& cc) {
      allrows.insert(allrows.end(), rr.begin(), rr.end());
      allcols.insert(allcols.end(), cc.begin(), cc.end());
      rowidx.push_back(static_cast<int>(rr.size()));
      colidx.push_back(static_cast<int>(cc.size()));
    };
    for (int b = 0; b < 64; ++b) {  // neighbouring 4^3 cubes
      const int64_t ci = pick(i0, n - 4), cj = pick(0, n - 4), ck = pick(0, n - 4);
      const int64_t di = std::clamp<int64_t>(ci + 4 * pick(-1, 1), i0, n - 4), dj = std::clamp<int64_t>(cj + 4 * pick(-1, 1), 0, n - 4),
                    dk = std::clamp<int64_t>(ck + 4 * pick(-1, 1), 0, n - 4);
      std::vector<int> rr, cc;
      cube(ci, cj, ck, 4, rr);
      cube(di, dj, dk, 4, cc);
      add(rr, cc);
    }
    for (int b = 0; b < 16; ++b) {  // 16 points of a cube against a 16^3 region away from it, and the transpose
      const int64_t ci = pick(i0, n - 4), cj = pick(0, n - 4), ck = pick(0, n - 4);
      std::vector<int> pts, far;
      cube(ci, cj, ck, 4, pts);
      std::shuffle(pts.begin(), pts.end(), rng);
      pts.resize(16);
      int64_t fi = 0, fj = 0, fk = 0;
      for (int tries = 0;; ++tries) {
        fi = pick(i0, n - 16);
        fj = pick(0, n - 16);
        fk = pick(0, n - 16);
        const bool apart = ci + 8 <= fi || fi + 20 <= ci || cj + 8 <= fj || fj + 20 <= cj || ck + 8 <= fk || fk + 20 <= ck;
        if (apart || tries > 100) break;
      }
      cube(fi, fj, fk, 16, far);
      add(pts, far);
      add(far, pts);
    }
    const int ninter = static_cast<int>(rowidx.size());
    std::vector<int> pgidx(ninter), pmaps(3 * ninter);
    int64_t nloc = 0;
    for (int k = 0; k < ninter; ++k) {
      pgidx[k] = k;
      pmaps[k] = pmaps[ninter + k] = 1;
      pmaps[2 * ninter + k] = k % mpi_size;
      if (k % mpi_size == rank) nloc += static_cast<int64_t>(rowidx[k]) * colidx[k];
      entries += static_cast<double>(rowidx[k]) * colidx[k];
    }
    std::vector<double> dat(static_cast<size_t>(std::max<int64_t>(nloc, 1)));
    int nallrows = static_cast<int>(allrows.size()), nallcols = static_cast<int>(allcols.size()), npmap = ninter, ni = ninter;
    MPI_Barrier(MPI_COMM_WORLD);
    const double t0 = MPI_Wtime();
    d_c_bpack_extractelement(&resources.matrix, &resources.option, &resources.mesh, &resources.stats, &resources.process_tree,
                             &ni, &nallrows, &nallcols, &nloc, allrows.data(), allcols.data(), dat.data(), rowidx.data(),
                             colidx.data(), pgidx.data(), &npmap, pmaps.data());
    MPI_Barrier(MPI_COMM_WORLD);
    t_total += MPI_Wtime() - t0;
  }
  if (rank == 0) {
    std::cout << std::setprecision(3) << "H2 partial extract bench: " << o.partial_extract_bench << " rounds of 96 blocks, "
              << entries / o.partial_extract_bench / 1e6 << " M entries a round, " << t_total / o.partial_extract_bench
              << " s a round, " << entries / t_total / 1e6 << " M entries/s" << std::endl;
  }
}

// --H2_partial, after partial_checks: the release of the Schur complement.
// The partial solves give what they gave before it, and the Schur
// complement's entries and products are refused.
void partial_release_check(const DriverOptions& o, ButterflyResources& resources, int rank, int local_points,
                           const DistributedInput& input) {
  const int64_t N = checked_cube(o.grid_size, "grid_size");
  const int n22 = static_cast<int>(N - o.h2_partial), nrhs = o.nrhs;
  if (n22 == 0) return;
  std::vector<double> b(static_cast<size_t>(local_points) * nrhs), x2(static_cast<size_t>(n22) * nrhs);
  for (size_t i = 0; i < b.size(); ++i) b[i] = 1.0 + 0.001 * static_cast<double>((input.first + static_cast<int64_t>(i)) % 89);
  for (size_t i = 0; i < x2.size(); ++i) x2[i] = std::sin(0.37 * static_cast<double>(i));
  int nloc = local_points, N22 = n22, Nrhs = nrhs;
  auto solves = [&](std::vector<double>& out2, std::vector<double>& x) {
    out2.assign(x2.size(), 0.0);
    x.assign(b.size(), 0.0);
    d_c_bpack_partial_forward(out2.data(), b.data(), &nloc, &N22, &Nrhs, &resources.matrix, &resources.option,
                              &resources.stats, &resources.process_tree);
    d_c_bpack_partial_backward(x.data(), x2.data(), &nloc, &N22, &Nrhs, &resources.matrix, &resources.option,
                               &resources.stats, &resources.process_tree);
  };
  std::vector<double> out2a, xa, out2b, xb;
  solves(out2a, xa);
  double rss[2] = {resident_gb(), 0.0};
  d_c_bpack_partial_schur_release(&resources.matrix, &resources.option);
  rss[1] = resident_gb();
  solves(out2b, xb);
  double diff = 0.0, scale = 0.0;
  for (size_t i = 0; i < out2a.size(); ++i) {
    diff = std::max(diff, std::abs(out2a[i] - out2b[i]));
    scale = std::max(scale, std::abs(out2a[i]));
  }
  for (size_t i = 0; i < xa.size(); ++i) {
    diff = std::max(diff, std::abs(xa[i] - xb[i]));
    scale = std::max(scale, std::abs(xa[i]));
  }
  int refused = 0;
  try {
    std::vector<double> y2(x2.size());
    d_c_bpack_partial_schur_mult(y2.data(), x2.data(), &N22, &Nrhs, &resources.matrix, &resources.option);
  } catch (const std::exception&) {
    refused = 1;
  }
  try {
    double v = 0.0;
    int one = 1, r0 = 0;
    d_c_bpack_partial_schur_extract(&v, &r0, &one, &r0, &one, &resources.matrix, &resources.option);
  } catch (const std::exception&) {
    refused += 1;
  }
  double drop[2] = {rss[0] - rss[1], rss[0]};
  MPI_Allreduce(MPI_IN_PLACE, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, drop, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &refused, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << std::setprecision(3) << "H2 partial release: resident set " << drop[1] << " -> " << drop[1] - drop[0]
              << " GB (summed over ranks); the solves after it differ by " << diff / std::max(scale, 1e-300)
              << " (relative); the Schur complement's product and entries " << (refused == 2 ? "refused" : "NOT refused")
              << std::endl;
  }
  if (refused != 2) throw std::runtime_error("partial_release_check: a call on the released Schur complement went through");
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  int mpi_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int return_code = 0;
  MPI_Comm source_comm = MPI_COMM_NULL;
  try {
    const DriverOptions driver_options = parse_driver_options(argc, argv);
    if (driver_options.show_help) {
      if (rank == 0) print_usage(argv[0]);
      MPI_Finalize();
      return 0;
    }

    const int64_t point_count_64 =
        checked_cube(driver_options.grid_size, "grid_size");
    int dimension = 3;
    Laplace3DApplication application(
        driver_options.grid_size, driver_options.distributed64 == 0, driver_options.skew,
        driver_options.column_weight);

    if (rank == 0) {
      std::cout << "=== ButterflyPACK H2 3D Laplace Test ===\n"
                << "Grid size: " << driver_options.grid_size << " x "
                << driver_options.grid_size << " x "
                << driver_options.grid_size << "\n"
                << "Total points: " << point_count_64 << "\n"
                << "Format: " << driver_options.format << "\n"
                << "sym: " << driver_options.sym << "\n"
                << (driver_options.skew != 0.0 ? "Skew (dipole) term: " + std::to_string(driver_options.skew) + "\n" : "")
                << (driver_options.column_weight != 0.0
                        ? "Column weights: amplitude " + std::to_string(driver_options.column_weight) +
                              ", sum of log w " + std::to_string(application.log_weight_sum()) +
                              " (log det(A W) - log det(A))\n"
                        : "")
                << "C API: "
                << (driver_options.distributed64 == 1
                        ? "distributed 64-bit"
                        : "legacy 32-bit")
                << "\n"
                << "Tolerance: " << driver_options.tolerance << "\n"
                << "Nmin_leaf: " << driver_options.nmin_leaf << "\n"
                << "Reduction threshold: "
                << driver_options.reduction_threshold << "\n"
                << "H2_unstructured: "
                << driver_options.h2_unstructured << "\n"
                << "CA_level: " << driver_options.ca_level << "\n"
                << "H2_use_sketch: " << driver_options.h2_use_sketch << "\n"
                << "H2_lazy_schur: " << driver_options.h2_lazy_schur << "\n"
                << "H2_compress_factors: " << driver_options.h2_compress_factors << "\n"
                << "H2_GEMM_split: " << driver_options.h2_gemm_split << "\n"
                << "H2_XRR_factor: " << driver_options.h2_xrr_factor << "\n"
                << "H2_use_gpu: " << driver_options.h2_use_gpu << "\n"
                << "HODLR_use_gpu: " << driver_options.hodlr_use_gpu << "\n"
                << "GPU evaluator: " << driver_options.gpu_evaluator << "\n"
                << "elem_extract: " << driver_options.elem_extract << "\n"
                << "HODLR from the entries of one on ranks: "
                << (driver_options.hodlr_extract_ranks > 0
                        ? std::to_string(driver_options.hodlr_extract_ranks)
                        : std::string("no"))
                << "\n"
                << "H2_CA_staged_halo: "
                << driver_options.h2_ca_staged_halo << "\n"
                << "H2_CA_owner_component: "
                << driver_options.h2_ca_owner_component << "\n"
                << "H2_CA_pack_shared: "
                << driver_options.h2_ca_pack_shared << "\n"
                << "H2_CA_solve: "
                << driver_options.h2_ca_solve << "\n"
                << "H2_CA_owner_serial: "
                << driver_options.h2_ca_owner_serial << "\n"
                << "Number of RHS: " << driver_options.nrhs << "\n"
                << "Gaussian potential: disabled\n"
                << "Proxy points: 0 (current ButterflyPACK H2 interface)\n"
                << "Diagonal self-cell integral: " << std::setprecision(17)
                << application.diagonal() << std::setprecision(6) << std::endl;
    }

    ButterflyResources resources;
    d_c_bpack_createoption(&resources.option);
    d_c_bpack_createstats(&resources.stats);

    std::vector<int> groups(static_cast<size_t>(mpi_size));
    std::iota(groups.begin(), groups.end(), 0);
    MPI_Fint fortran_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
    d_c_bpack_createptree(
        &mpi_size, groups.data(), &fortran_comm, &resources.process_tree);

    // the options of a matrix built here (with --hodlr_extract_ranks the
    // source and the HODLR under test alike)
    auto set_options = [&](F2Cptr* option_handle) {
      d_c_bpack_set_D_option(
          option_handle, "tol_comp", driver_options.tolerance);
      d_c_bpack_set_I_option(
          option_handle, "format", driver_options.format);
      d_c_bpack_set_I_option(option_handle, "sym", driver_options.sym);
      if (driver_options.format == 1) {
        // Keep the forward HODLR blocks after factorization (the unsymmetric
        // factorization frees them otherwise) for the solution check below,
        // and run the library's solve-error test.
        d_c_bpack_set_I_option(option_handle, "ErrSol", 1);
      }
      d_c_bpack_set_I_option(
          option_handle, "Nmin_leaf", static_cast<int>(driver_options.nmin_leaf));
      d_c_bpack_set_I_option(
          option_handle, "reduction_threshold",
          static_cast<int>(driver_options.reduction_threshold));
      d_c_bpack_set_I_option(
          option_handle, "H2_unstructured",
          driver_options.h2_unstructured);
      d_c_bpack_set_I_option(
          option_handle, "CA_level", driver_options.ca_level);
      d_c_bpack_set_I_option(
          option_handle, "H2_use_sketch", driver_options.h2_use_sketch);
      d_c_bpack_set_I_option(
          option_handle, "H2_lazy_schur", driver_options.h2_lazy_schur);
      d_c_bpack_set_I_option(
          option_handle, "H2_compress_factors", driver_options.h2_compress_factors);
      d_c_bpack_set_I_option(
          option_handle, "H2_GEMM_split", driver_options.h2_gemm_split);
      d_c_bpack_set_I_option(
          option_handle, "H2_XRR_factor", driver_options.h2_xrr_factor);
      d_c_bpack_set_I_option(
          option_handle, "H2_use_gpu", driver_options.h2_use_gpu);
      d_c_bpack_set_I_option(
          option_handle, "HODLR_use_gpu", driver_options.hodlr_use_gpu);
      if (driver_options.hodlr_gpu_pieces >= 1) {
        d_c_bpack_set_I_option(
            option_handle, "HODLR_gpu_pieces", driver_options.hodlr_gpu_pieces);
      }
      if (driver_options.reclr_leaf >= 0) {
        d_c_bpack_set_I_option(
            option_handle, "RecLR_leaf", driver_options.reclr_leaf);
      }
      if (driver_options.baca_batch > 0) {
        d_c_bpack_set_I_option(
            option_handle, "BACA_Batch", driver_options.baca_batch);
      }
      if (driver_options.h2_id_radius >= 0) {
        d_c_bpack_set_I_option(
            option_handle, "H2_ID_radius", driver_options.h2_id_radius);
      }
      if (driver_options.h2_id_proxy >= 0) {
        d_c_bpack_set_I_option(
            option_handle, "H2_ID_proxy", driver_options.h2_id_proxy);
      }
      if (driver_options.h2_id_proxy_points >= 0) {
        d_c_bpack_set_I_option(
            option_handle, "H2_ID_proxy_points",
            driver_options.h2_id_proxy_points);
      }
      d_c_bpack_set_D_option(option_handle, "H2_s_op", driver_options.h2_s_op);
      d_c_bpack_set_I_option(
          option_handle, "H2_CA_staged_halo",
          driver_options.h2_ca_staged_halo);
      d_c_bpack_set_I_option(
          option_handle, "H2_CA_owner_component",
          driver_options.h2_ca_owner_component);
      d_c_bpack_set_I_option(
          option_handle, "H2_CA_pack_shared",
          driver_options.h2_ca_pack_shared);
      d_c_bpack_set_I_option(
          option_handle, "H2_CA_solve",
          driver_options.h2_ca_solve);
      d_c_bpack_set_I_option(
          option_handle, "H2_CA_owner_serial",
          driver_options.h2_ca_owner_serial);
      d_c_bpack_set_I_option(option_handle, "precon", driver_options.precon);
      d_c_bpack_set_I_option(option_handle, "H2_partial", static_cast<int>(driver_options.h2_partial));
      d_c_bpack_set_I_option(option_handle, "H2_partial_extract", driver_options.h2_partial_extract);
      d_c_bpack_set_I_option(option_handle, "H2_extract", driver_options.h2_extract);
      d_c_bpack_set_I_option(
          option_handle, "verbosity", driver_options.verbosity);
      d_c_bpack_set_I_option(
          option_handle, "elem_extract", driver_options.elem_extract);
      d_c_bpack_set_I_option(
          option_handle, "LRlevel", driver_options.lrlevel);
      d_c_bpack_set_I_option(option_handle, "cpp", 1);
      d_c_bpack_set_I_option(option_handle, "nogeo", 0);
    };

    // with --hodlr_extract_ranks: the source of the entries, a HODLR of the
    // kernel on the first ranks
    ExtractionSource extraction;
    ButterflyResources source;
    F2Cptr source_kernel_quantities = nullptr;
    const int source_ranks = driver_options.hodlr_extract_ranks;
    if (source_ranks > mpi_size) {
      throw std::invalid_argument("hodlr_extract_ranks exceeds the number of ranks");
    }
    if (source_ranks > 0) {
      MPI_Comm_split(MPI_COMM_WORLD, rank < source_ranks ? 0 : MPI_UNDEFINED,
                     rank, &source_comm);
    }
    if (source_ranks > 0 && rank < source_ranks) {
      d_c_bpack_createoption(&source.option);
      d_c_bpack_createstats(&source.stats);
      std::vector<int> source_groups(static_cast<size_t>(source_ranks));
      std::iota(source_groups.begin(), source_groups.end(), 0);
      MPI_Fint source_fortran_comm = MPI_Comm_c2f(source_comm);
      int source_size = source_ranks;
      d_c_bpack_createptree(&source_size, source_groups.data(),
                            &source_fortran_comm, &source.process_tree);
      set_options(&source.option);
      // (not factored, or for an H2 compressed only; the Laplace block
      // callback on its own ranks)
      d_c_bpack_set_I_option(&source.option, "ErrSol", 0);
      d_c_bpack_set_I_option(&source.option, "elem_extract", 2);
      if (driver_options.format == 7) {
        // (with --H2_partial the H2 under test is the partial one)
        d_c_bpack_set_I_option(&source.option, "precon", 2);
        d_c_bpack_set_I_option(&source.option, "H2_partial", 0);
      }
      int point_count = static_cast<int>(point_count_64);
      int nlevel = 0;
      int user_tree = point_count;
      int source_local_points = 0;
      std::vector<int> permutation(static_cast<size_t>(point_count));
      d_c_bpack_construct_init(
          &point_count, &dimension, application.locations(), nullptr,
          &nlevel, &user_tree, permutation.data(), &source_local_points,
          &source.matrix, &source.option, &source.stats, &source.mesh,
          &source_kernel_quantities, &source.process_tree,
          &distance_callback, &near_far_callback, &application);
#ifdef BPACK_EXAMPLE_GPU
      if ((driver_options.format == 7 ? driver_options.h2_use_gpu : driver_options.hodlr_use_gpu) > 0 &&
          driver_options.gpu_evaluator == 1) {
        laplace3d_gpu_register_general(&source.matrix, application.inverse_4pi_n(), application.diagonal(),
                                       application.skew(), application.skew_direction(),
                                       application.column_weight());
      }
#endif
      d_c_bpack_construct_element_compute(
          &source.matrix, &source.option, &source.stats, &source.mesh,
          &source_kernel_quantities, &source.process_tree,
          &laplace_entry_callback, &laplace_block_callback, &application);
      if (driver_options.format == 7) {
        // the compression-only H2 (precon 2): the source of the entries
        d_c_bpack_factor(&source.matrix, &source.option, &source.stats,
                         &source.process_tree, &source.mesh);
      }
      extraction.matrix = source.matrix;
      extraction.option = source.option;
      extraction.mesh = source.mesh;
      extraction.stats = source.stats;
      extraction.process_tree = source.process_tree;
      if (rank == 0) {
        const char* kind = driver_options.format == 7 ? "H2" : "HODLR";
        std::cout << "\nBuilt the source " << kind << " on " << source_ranks
                  << " ranks; the " << kind << " under test from its entries on "
                  << mpi_size << " ranks:" << std::endl;
      }
    }

    set_options(&resources.option);
    if (source_ranks > 0) {
      // the HODLR under test in the source's tree order (natural order, no
      // geometry), its entries by blocks from the source; an H2 keeps its
      // geometry (the source's indices are the points', as its own)
      if (driver_options.format == 1) {
        d_c_bpack_set_I_option(&resources.option, "nogeo", 1);
        d_c_bpack_set_I_option(&resources.option, "xyzsort", 0);
        d_c_bpack_set_I_option(&resources.option, "knn", 0);
      }
      extraction.caller_tree = resources.process_tree;
    }

    int local_points = 0;
    F2Cptr kernel_quantities = nullptr;
    DistributedInput distributed_input;
    if (driver_options.distributed64 == 1) {
      distributed_input = make_distributed_input(
          point_count_64, rank, mpi_size, application);
      if (distributed_input.count > std::numeric_limits<int>::max()) {
        throw std::overflow_error(
            "the local input count exceeds the current solve API limit");
      }
      int bounds_provided = 1;
      const double global_bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
      int64_t internal_local_points = 0;
      d_c_bpack_construct_init_distributed64(
          &point_count_64, &distributed_input.count, &dimension,
          distributed_input.global_ids.data(),
          distributed_input.coordinates.data(), &bounds_provided,
          global_bounds, &internal_local_points, &resources.matrix,
          &resources.option, &resources.stats, &resources.mesh,
          &kernel_quantities, &resources.process_tree);
      local_points = static_cast<int>(distributed_input.count);
    } else {
      int point_count = static_cast<int>(point_count_64);
      int nlevel = 0;
      int user_tree = point_count;
      std::vector<int> permutation(static_cast<size_t>(point_count));
      d_c_bpack_construct_init(
          &point_count, &dimension, application.locations(), nullptr,
          &nlevel, &user_tree, permutation.data(), &local_points,
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &distance_callback, &near_far_callback, &application);
      // the local rows are the rank's range of the tree order: their input
      // points (1-based), for the partial checks
      int64_t count = local_points, first = 0;
      MPI_Exscan(&count, &first, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0) first = 0;
      distributed_input.first = first;
      distributed_input.count = count;
      distributed_input.global_ids.assign(permutation.begin() + first, permutation.begin() + first + count);
    }

    d_c_bpack_printoption(&resources.option, &resources.process_tree);
#ifdef BPACK_EXAMPLE_GPU
    if (source_ranks == 0 && driver_options.gpu_evaluator == 1) {
      // GPU evaluator of the entries for H2_use_gpu / HODLR_use_gpu
      // (gpu/laplace3d_gpu.cu): 1 / (4 pi N r), with the self-cell integral
      // on the diagonal
      double h2_use_gpu = 0.0, hodlr_use_gpu = 0.0;
      d_c_bpack_getoption(&resources.option, "H2_use_gpu", &h2_use_gpu);
      d_c_bpack_getoption(&resources.option, "HODLR_use_gpu", &hodlr_use_gpu);
      if (h2_use_gpu > 0 || hodlr_use_gpu > 0) {
        laplace3d_gpu_register_general(&resources.matrix, application.inverse_4pi_n(), application.diagonal(),
                                       application.skew(), application.skew_direction(),
                                       application.column_weight());
      }
    }
#endif
    if (driver_options.distributed64 == 1) {
      d_c_bpack_construct_element_compute_distributed64(
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &laplace_entry_callback64, &laplace_block_callback64, &application);
    } else if (source_ranks > 0) {
      // an H2 asks for the ID's sketches of its kernel rows sketched (a
      // nonsymmetric one also for the transposed sketches)
      if (driver_options.format == 7) {
        d_c_bpack_set_sketch_callback(&resources.matrix, &extraction_sketch_callback);
        d_c_bpack_set_sketch_transpose_callback(&resources.matrix, &extraction_sketch_transpose_callback);
      }
      d_c_bpack_construct_element_compute(
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &extraction_entry_callback, &extraction_block_callback,
          &extraction);
    } else {
      d_c_bpack_construct_element_compute(
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &laplace_entry_callback, &laplace_block_callback, &application);
    }

    if (rank == 0) {
      std::cout << "\nFactoring the 3D Laplace operator:" << std::endl;
    }
    d_c_bpack_factor(
        &resources.matrix, &resources.option, &resources.stats,
        &resources.process_tree, &resources.mesh);

    // --extract_check: random blocks of the compression-only H2 through the
    // extraction bridge (each to a random rank), against the kernel; then
    // this rank's rows of two columns against the compressed matvec of their
    // unit vectors (the same representation: rounding)
    if (driver_options.extract_check > 0 && driver_options.format == 7 &&
        driver_options.precon == 2 && driver_options.distributed64 == 1) {
      const int64_t n_global = point_count_64;
      std::mt19937_64 rng(20261016);
      auto pick = [&](int count, std::vector<int>& out) {
        // a random box-sized cluster of consecutive ids, or scattered ids
        const bool cluster = rng() % 2 == 0;
        const int64_t first = static_cast<int64_t>(rng() % static_cast<uint64_t>(n_global));
        for (int i = 0; i < count; ++i) {
          const int64_t id = cluster ? (first + i * 7) % n_global
                                     : static_cast<int64_t>(rng() % static_cast<uint64_t>(n_global));
          out.push_back(static_cast<int>(id + 1));
        }
      };
      int ninter = driver_options.extract_check;
      std::vector<int> allrows, allcols, rowidx, colidx, pgidx;
      for (int n = 0; n < ninter; ++n) {
        const int nr = 1 + static_cast<int>(rng() % 60), nc = 1 + static_cast<int>(rng() % 60);
        pick(nr, allrows);
        pick(nc, allcols);
        rowidx.push_back(nr);
        colidx.push_back(nc);
        pgidx.push_back(static_cast<int>(rng() % static_cast<uint64_t>(mpi_size)));
      }
      // this rank's rows of two random columns, owned here
      std::vector<int> probe(2);
      for (int& j : probe) j = static_cast<int>(rng() % static_cast<uint64_t>(n_global)) + 1;
      for (int p = 0; p < mpi_size; ++p) {
        for (int j : probe) {
          // (every rank must pass the same list: rank p's rows)
          const int64_t base = n_global / mpi_size, extra = n_global % mpi_size;
          const int64_t first = p * base + std::min<int64_t>(p, extra);
          const int64_t count = base + (p < extra ? 1 : 0);
          for (int64_t i = 0; i < count; ++i) allrows.push_back(static_cast<int>(first + i + 1));
          allcols.push_back(j);
          rowidx.push_back(static_cast<int>(count));
          colidx.push_back(1);
          pgidx.push_back(p);
        }
      }
      ninter = static_cast<int>(rowidx.size());
      std::vector<int> pmaps(3 * mpi_size);
      for (int p = 0; p < mpi_size; ++p) {
        pmaps[p] = 1;
        pmaps[mpi_size + p] = 1;
        pmaps[2 * mpi_size + p] = p;
      }
      int64_t mine = 0;
      for (int n = 0; n < ninter; ++n) {
        if (pgidx[n] == rank) mine += static_cast<int64_t>(rowidx[n]) * colidx[n];
      }
      std::vector<double> values(static_cast<size_t>(std::max<int64_t>(mine, 1)));
      int nallrows = static_cast<int>(allrows.size()), nallcols = static_cast<int>(allcols.size());
      int npmap = mpi_size;
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      d_c_bpack_extractelement_bridge(
          &resources.matrix, &resources.option, &resources.mesh, &resources.stats,
          &resources.process_tree, &resources.process_tree, &ninter, &nallrows,
          &nallcols, &mine, allrows.data(), allcols.data(), values.data(),
          rowidx.data(), colidx.data(), pgidx.data(), &npmap, pmaps.data());
      double t_extract = MPI_Wtime() - t0;
      MPI_Allreduce(MPI_IN_PLACE, &t_extract, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      {  // where that call's time went (core/extraction_profile.hpp)
        double pr[12] = {0.0};
        d_c_bpack_get_extraction_profile(pr);
        MPI_Allreduce(MPI_IN_PLACE, pr, 12, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) {
          static const char* names[8] = {"resolve", "tasks",   "vectors up", "column side",
                                         "compute", "pack",    "route",      "unpack"};
          std::printf("  extraction profile (max over ranks, %lld calls): ", static_cast<long long>(pr[11]));
          for (int k = 0; k < 8; ++k) std::printf("%s %.4f  ", names[k], pr[k]);
          std::printf("| entries computed %.0f\n", pr[9]);
          std::fflush(stdout);
        }
      }

      // against the kernel (the random blocks), and the matvec (the probes)
      std::vector<double> unit(static_cast<size_t>(local_points) * 2, 0.0);
      for (int c = 0; c < 2; ++c) {
        const int64_t local = probe[static_cast<size_t>(c)] - 1 - distributed_input.first;
        if (local >= 0 && local < local_points) unit[static_cast<size_t>(local + c * local_points)] = 1.0;
      }
      std::vector<double> product(unit.size(), 0.0);
      {
        int two = 2, n_in = local_points, n_out = local_points;
        const char trans = 'N';
        d_c_bpack_mult(&trans, unit.data(), product.data(), &n_in, &n_out, &two,
                       &resources.matrix, &resources.option, &resources.stats,
                       &resources.process_tree);
      }
      double kernel_error2 = 0.0, kernel_norm2 = 0.0, matvec_error = 0.0, matvec_scale = 0.0;
      int64_t entries = 0;
      int64_t at = 0, row_at = 0, col_at = 0;
      const int random_blocks = driver_options.extract_check;
      for (int n = 0; n < ninter; ++n) {
        const int nr = rowidx[n], nc = colidx[n];
        if (pgidx[n] == rank) {
          for (int j = 0; j < nc; ++j) {
            for (int i = 0; i < nr; ++i) {
              const double v = values[static_cast<size_t>(at + i + static_cast<int64_t>(j) * nr)];
              if (n < random_blocks) {
                // (distributed64: the coordinates from the ids, none stored)
                const int64_t a = allrows[static_cast<size_t>(row_at + i)], b = allcols[static_cast<size_t>(col_at + j)];
                double x[3], y[3];
                application.coordinate(a - 1, x);
                application.coordinate(b - 1, y);
                const double k = application.evaluate_coordinates(a, b, x, y);
                kernel_error2 += (v - k) * (v - k);
                kernel_norm2 += k * k;
              } else {
                const int c = ((n - random_blocks) % 2);
                const double y = product[static_cast<size_t>(i + c * local_points)];
                matvec_error = std::max(matvec_error, std::abs(v - y));
                matvec_scale = std::max(matvec_scale, std::abs(y));
              }
            }
          }
          at += static_cast<int64_t>(nr) * nc;
          entries += static_cast<int64_t>(nr) * nc;
        }
        row_at += nr;
        col_at += nc;
      }
      MPI_Allreduce(MPI_IN_PLACE, &kernel_error2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &kernel_norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &matvec_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &matvec_scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &entries, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

      // the random blocks sketched (c_bpack_extractsketch_bridge), against
      // the same sketch of their extracted entries
      int nsketch = random_blocks;
      std::vector<int> sketchdim, omegaidx, omega_row, omega_dest;
      std::vector<double> omega_val;
      int nsrows = 0, nscols = 0;
      for (int n = 0; n < nsketch; ++n) {
        const int d = 1 + static_cast<int>(rng() % 40);
        sketchdim.push_back(d);
        for (int i = 0; i < rowidx[n]; ++i) {
          for (int e = 0; e < 4; ++e) {
            omega_row.push_back(i + 1);
            omega_dest.push_back(static_cast<int>(rng() % static_cast<uint64_t>(d)) + 1);
            omega_val.push_back(rng() % 2 == 0 ? 0.5 : -0.5);
          }
        }
        omegaidx.push_back(4 * rowidx[n]);
        nsrows += rowidx[n];
        nscols += colidx[n];
      }
      int64_t sketch_mine = 0;
      for (int n = 0; n < nsketch; ++n) {
        if (pgidx[n] == rank) sketch_mine += static_cast<int64_t>(sketchdim[n]) * colidx[n];
      }
      std::vector<double> sketched(static_cast<size_t>(std::max<int64_t>(sketch_mine, 1)));
      int64_t nallomega = static_cast<int64_t>(omega_val.size());
      MPI_Barrier(MPI_COMM_WORLD);
      const double t_sketch0 = MPI_Wtime();
      d_c_bpack_extractsketch_bridge(
          &resources.matrix, &resources.option, &resources.mesh, &resources.stats,
          &resources.process_tree, &resources.process_tree, &nsketch, &nsrows, &nscols, &nallomega,
          &sketch_mine, allrows.data(), allcols.data(), rowidx.data(), colidx.data(), sketchdim.data(),
          omegaidx.data(), omega_row.data(), omega_dest.data(), omega_val.data(), sketched.data(),
          pgidx.data(), &npmap, pmaps.data());
      double t_sketch = MPI_Wtime() - t_sketch0;
      MPI_Allreduce(MPI_IN_PLACE, &t_sketch, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      double sketch_error = 0.0, sketch_scale = 0.0;
      {
        int64_t at_entries = 0, at_sketch = 0, e0 = 0;
        for (int n = 0; n < nsketch; ++n) {
          const int nr = rowidx[n], nc = colidx[n], d = sketchdim[n];
          if (pgidx[n] == rank) {
            std::vector<double> reference(static_cast<size_t>(d) * nc, 0.0);
            for (int64_t e = e0; e < e0 + omegaidx[n]; ++e) {
              const int r = omega_row[static_cast<size_t>(e)] - 1, i = omega_dest[static_cast<size_t>(e)] - 1;
              for (int j = 0; j < nc; ++j) {
                reference[static_cast<size_t>(i + j * d)] +=
                    omega_val[static_cast<size_t>(e)] * values[static_cast<size_t>(at_entries + r + static_cast<int64_t>(j) * nr)];
              }
            }
            for (size_t q = 0; q < reference.size(); ++q) {
              sketch_error = std::max(sketch_error, std::abs(sketched[static_cast<size_t>(at_sketch) + q] - reference[q]));
              sketch_scale = std::max(sketch_scale, std::abs(reference[q]));
            }
            at_entries += static_cast<int64_t>(nr) * nc;
            at_sketch += static_cast<int64_t>(d) * nc;
          }
          e0 += omegaidx[n];
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &sketch_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sketch_scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      // the same transposed (c_bpack_extractsketch_transpose_bridge): Omega
      // A(J, I)^T against that sketch of the extracted A(J, I)
      std::vector<double> values_t;
      {
        int64_t mine_t = 0;
        for (int n = 0; n < nsketch; ++n) {
          if (pgidx[n] == rank) mine_t += static_cast<int64_t>(rowidx[n]) * colidx[n];
        }
        values_t.assign(static_cast<size_t>(std::max<int64_t>(mine_t, 1)), 0.0);
        d_c_bpack_extractelement_bridge(
            &resources.matrix, &resources.option, &resources.mesh, &resources.stats,
            &resources.process_tree, &resources.process_tree, &nsketch, &nscols, &nsrows, &mine_t,
            allcols.data(), allrows.data(), values_t.data(), colidx.data(), rowidx.data(), pgidx.data(), &npmap,
            pmaps.data());
      }
      std::vector<double> sketched_t(static_cast<size_t>(std::max<int64_t>(sketch_mine, 1)));
      MPI_Barrier(MPI_COMM_WORLD);
      const double t_sketch_t0 = MPI_Wtime();
      d_c_bpack_extractsketch_transpose_bridge(
          &resources.matrix, &resources.option, &resources.mesh, &resources.stats,
          &resources.process_tree, &resources.process_tree, &nsketch, &nsrows, &nscols, &nallomega,
          &sketch_mine, allrows.data(), allcols.data(), rowidx.data(), colidx.data(), sketchdim.data(),
          omegaidx.data(), omega_row.data(), omega_dest.data(), omega_val.data(), sketched_t.data(),
          pgidx.data(), &npmap, pmaps.data());
      double t_sketch_t = MPI_Wtime() - t_sketch_t0;
      MPI_Allreduce(MPI_IN_PLACE, &t_sketch_t, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      double sketch_t_error = 0.0, sketch_t_scale = 0.0;
      {
        int64_t at_entries = 0, at_sketch = 0, e0 = 0;
        for (int n = 0; n < nsketch; ++n) {
          const int nr = rowidx[n], nc = colidx[n], d = sketchdim[n];
          if (pgidx[n] == rank) {
            std::vector<double> reference(static_cast<size_t>(d) * nc, 0.0);
            for (int64_t e = e0; e < e0 + omegaidx[n]; ++e) {
              const int r = omega_row[static_cast<size_t>(e)] - 1, i = omega_dest[static_cast<size_t>(e)] - 1;
              for (int j = 0; j < nc; ++j) {  // (A(J, I): nc x nr)
                reference[static_cast<size_t>(i + j * d)] +=
                    omega_val[static_cast<size_t>(e)] * values_t[static_cast<size_t>(at_entries + j + static_cast<int64_t>(r) * nc)];
              }
            }
            for (size_t q = 0; q < reference.size(); ++q) {
              sketch_t_error = std::max(sketch_t_error, std::abs(sketched_t[static_cast<size_t>(at_sketch) + q] - reference[q]));
              sketch_t_scale = std::max(sketch_t_scale, std::abs(reference[q]));
            }
            at_entries += static_cast<int64_t>(nr) * nc;
            at_sketch += static_cast<int64_t>(d) * nc;
          }
          e0 += omegaidx[n];
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &sketch_t_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sketch_t_scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) {
        std::cout << std::setprecision(4) << "H2 extraction check: " << ninter << " blocks, " << entries
                  << " entries in " << t_extract << " s; random blocks against the kernel "
                  << std::sqrt(kernel_error2 / kernel_norm2) << " (Frobenius, relative); "
                  << "rows of 2 columns against the matvec " << matvec_error / matvec_scale
                  << " (max, relative); the random blocks sketched (" << t_sketch
                  << " s) against their entries sketched " << sketch_error / sketch_scale << " (max, relative)"
                  << "; transposed (" << t_sketch_t << " s) " << sketch_t_error / sketch_t_scale
                  << std::setprecision(6) << std::endl;
      }
    }

    if (driver_options.h2_partial == 0) {
    int number_of_rhs = driver_options.nrhs;
    const size_t value_count = static_cast<size_t>(local_points) *
        static_cast<size_t>(number_of_rhs);
    std::vector<double> rhs(value_count);
    for (int column = 0; column < number_of_rhs; ++column) {
      for (int row = 0; row < local_points; ++row) {
        const int64_t global_row = driver_options.distributed64 == 1
            ? distributed_input.first + row
            : static_cast<int64_t>(row) + 17 * rank;
        rhs[static_cast<size_t>(row) +
            static_cast<size_t>(column) * local_points] =
            1.0 + 0.125 * column +
            0.001 * (global_row % 97);
      }
    }

    if (rank == 0) {
      std::cout << "\nSolving the 3D Laplace system:" << std::endl;
    }
    std::vector<double> solution(value_count, 0.0);
    d_c_bpack_solve(
        solution.data(), rhs.data(), &local_points,
        &number_of_rhs,
        &resources.matrix, &resources.option, &resources.stats,
        &resources.process_tree);

    {
      std::vector<double> product(value_count, 0.0);
      int output_local_points = local_points;
      const char trans = 'N';
      d_c_bpack_mult(
          &trans, solution.data(), product.data(),
          &local_points, &output_local_points, &number_of_rhs,
          &resources.matrix, &resources.option, &resources.stats,
          &resources.process_tree);

      std::vector<double> solution_sum(static_cast<size_t>(number_of_rhs), 0.0);
      std::vector<double> solution_norm2(static_cast<size_t>(number_of_rhs), 0.0);
      std::vector<double> residual_norm2(static_cast<size_t>(number_of_rhs), 0.0);
      std::vector<double> rhs_norm2(static_cast<size_t>(number_of_rhs), 0.0);
      for (int column = 0; column < number_of_rhs; ++column) {
        for (int row = 0; row < local_points; ++row) {
          const size_t index = static_cast<size_t>(row) +
              static_cast<size_t>(column) * local_points;
          const double residual = product[index] - rhs[index];
          solution_sum[static_cast<size_t>(column)] += solution[index];
          solution_norm2[static_cast<size_t>(column)] +=
              solution[index] * solution[index];
          residual_norm2[static_cast<size_t>(column)] += residual * residual;
          rhs_norm2[static_cast<size_t>(column)] += rhs[index] * rhs[index];
        }
      }
      MPI_Allreduce(
          MPI_IN_PLACE, solution_sum.data(), number_of_rhs,
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(
          MPI_IN_PLACE, solution_norm2.data(), number_of_rhs,
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(
          MPI_IN_PLACE, residual_norm2.data(), number_of_rhs,
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(
          MPI_IN_PLACE, rhs_norm2.data(), number_of_rhs,
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if (rank == 0) {
        for (int column = 0; column < number_of_rhs; ++column) {
          const size_t index = static_cast<size_t>(column);
          const double relative_residual = std::sqrt(
              residual_norm2[index] / rhs_norm2[index]);
          std::cout << std::setprecision(17)
                    << (driver_options.format == 7 ? "H2" : "HODLR")
                    << " solution check RHS " << column
                    << ": sum=" << solution_sum[index]
                    << ", norm=" << std::sqrt(solution_norm2[index])
                    << ", relative residual=" << relative_residual
                    << std::setprecision(6) << std::endl;
        }
      }
    }

    auto mult = [&](char trans, const double* in, double* out, int ncol) {
      int n_in = local_points, n_out = local_points;
      d_c_bpack_mult(&trans, in, out, &n_in, &n_out, &ncol,
                     &resources.matrix, &resources.option, &resources.stats,
                     &resources.process_tree);
    };
    if (driver_options.inv_mult_check > 0 && driver_options.precon == 2) {
      bpack_example::check_mult_transpose<double>(
          MPI_COMM_WORLD, local_points, driver_options.inv_mult_check, mult,
          driver_options.format == 7 ? "H2" : "HODLR");
    } else if (driver_options.inv_mult_check > 0) {
      auto inv_mult = [&](char trans, const double* in, double* out, int ncol) {
        int n_in = local_points, n_out = local_points;
        d_c_bpack_inv_mult(&trans, in, out, &n_in, &n_out, &ncol,
                           &resources.matrix, &resources.option, &resources.stats,
                           &resources.process_tree);
      };
      bpack_example::check_inv_mult<double>(
          MPI_COMM_WORLD, local_points, driver_options.inv_mult_check, 5,
          rhs, solution, number_of_rhs, inv_mult, mult,
          driver_options.format == 7 ? "H2" : "HODLR");
    }

    // The Fortran entry points must reach the same H2 code as the C++ ones:
    // BPACK_Inv_Mult and BPACK_Mult have no format-7 case, so C_BPACK_Inv_Mult
    // and C_BPACK_Mult hand an H2 straight to the C++ wrapper.  This is the
    // path c_bf_schur_deterministic takes to F11's inverse, and before it
    // existed an H2 there returned its output array untouched.
    if (driver_options.inv_mult_check > 0 && driver_options.format == 7 &&
        driver_options.precon != 2) {
      const int ncol = 2;
      const size_t n = static_cast<size_t>(local_points) * ncol;
      std::vector<double> u(n), a(n, 0.0), b(n, 0.0);
      for (size_t i = 0; i < n; ++i)
        u[i] = std::sin(0.37 * static_cast<double>(i + 1) + 0.11 * rank);
      int n_in = local_points, n_out = local_points, nc = ncol;
      double worst = 0.0;
      for (char trans : {'N', 'T'}) {
        std::fill(a.begin(), a.end(), 0.0);
        std::fill(b.begin(), b.end(), 0.0);
        d_c_bpack_inv_mult(&trans, u.data(), a.data(), &n_in, &n_out, &nc,
                           &resources.matrix, &resources.option,
                           &resources.stats, &resources.process_tree);
        d_c_bpack_inv_mult_fortran(&trans, u.data(), b.data(), &n_in, &n_out, &nc,
                                   &resources.matrix, &resources.option,
                                   &resources.stats, &resources.process_tree);
        double d2 = 0.0, r2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
          d2 += (a[i] - b[i]) * (a[i] - b[i]);
          r2 += a[i] * a[i];
        }
        MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &r2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double rel = r2 > 0.0 ? std::sqrt(d2 / r2) : std::sqrt(d2);
        worst = std::max(worst, rel);
        std::fill(a.begin(), a.end(), 0.0);
        std::fill(b.begin(), b.end(), 0.0);
        d_c_bpack_mult(&trans, u.data(), a.data(), &n_in, &n_out, &nc,
                       &resources.matrix, &resources.option,
                       &resources.stats, &resources.process_tree);
        d_c_bpack_mult_fortran(&trans, u.data(), b.data(), &n_in, &n_out, &nc,
                               &resources.matrix, &resources.option,
                               &resources.stats, &resources.process_tree);
        d2 = 0.0; r2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
          d2 += (a[i] - b[i]) * (a[i] - b[i]);
          r2 += a[i] * a[i];
        }
        MPI_Allreduce(MPI_IN_PLACE, &d2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &r2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        worst = std::max(worst, r2 > 0.0 ? std::sqrt(d2 / r2) : std::sqrt(d2));
      }
      if (rank == 0) {
        std::cout << "H2 Fortran entry check (inv_mult and mult, N and T): "
                  << "worst relative difference from the C++ entry = "
                  << std::setprecision(3) << worst << std::setprecision(6)
                  << std::endl;
      }
    }

    // --logdet_check: the log-determinant and a solve of the factored H2
    // again after another H2 operator was built in the process (a
    // compression-only H2 of the same kernel on all ranks): the device data
    // is kept per operator, and each call must find this one's
    if (driver_options.logdet_check > 0 && driver_options.format == 7 && driver_options.precon != 2) {
      double phase_before = 0.0, logabs_before = 0.0;
      d_c_bpack_logdet(&phase_before, &logabs_before, &resources.option, &resources.matrix);
      ButterflyResources other;
      d_c_bpack_createoption(&other.option);
      d_c_bpack_createstats(&other.stats);
      d_c_bpack_createptree(&mpi_size, groups.data(), &fortran_comm, &other.process_tree);
      set_options(&other.option);
      d_c_bpack_set_I_option(&other.option, "precon", 2);
      d_c_bpack_set_I_option(&other.option, "ErrSol", 0);
      d_c_bpack_set_I_option(&other.option, "elem_extract", 2);
      F2Cptr other_quantities = nullptr;
      if (driver_options.distributed64 == 1) {
        int bounds_provided = 1;
        const double global_bounds[6] = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
        int64_t internal_local_points = 0;
        d_c_bpack_construct_init_distributed64(
            &point_count_64, &distributed_input.count, &dimension, distributed_input.global_ids.data(),
            distributed_input.coordinates.data(), &bounds_provided, global_bounds, &internal_local_points,
            &other.matrix, &other.option, &other.stats, &other.mesh, &other_quantities, &other.process_tree);
      } else {
        int point_count = static_cast<int>(point_count_64);
        int nlevel = 0;
        int user_tree = point_count;
        int other_local_points = 0;
        std::vector<int> permutation(static_cast<size_t>(point_count));
        d_c_bpack_construct_init(
            &point_count, &dimension, application.locations(), nullptr, &nlevel, &user_tree, permutation.data(),
            &other_local_points, &other.matrix, &other.option, &other.stats, &other.mesh, &other_quantities,
            &other.process_tree, &distance_callback, &near_far_callback, &application);
      }
#ifdef BPACK_EXAMPLE_GPU
      if (driver_options.h2_use_gpu > 0 && driver_options.gpu_evaluator == 1) {
        laplace3d_gpu_register_general(&other.matrix, application.inverse_4pi_n(), application.diagonal(),
                                       application.skew(), application.skew_direction(),
                                       application.column_weight());
      }
#endif
      if (driver_options.distributed64 == 1) {
        d_c_bpack_construct_element_compute_distributed64(
            &other.matrix, &other.option, &other.stats, &other.mesh, &other_quantities, &other.process_tree,
            &laplace_entry_callback64, &laplace_block_callback64, &application);
      } else {
        d_c_bpack_construct_element_compute(
            &other.matrix, &other.option, &other.stats, &other.mesh, &other_quantities, &other.process_tree,
            &laplace_entry_callback, &laplace_block_callback, &application);
      }
      d_c_bpack_factor(&other.matrix, &other.option, &other.stats, &other.process_tree, &other.mesh);
      double phase_after = 0.0, logabs_after = 0.0;
      d_c_bpack_logdet(&phase_after, &logabs_after, &resources.option, &resources.matrix);
      std::vector<double> again(value_count, 0.0);
      d_c_bpack_solve(again.data(), rhs.data(), &local_points, &number_of_rhs, &resources.matrix, &resources.option,
                      &resources.stats, &resources.process_tree);
      double difference = 0.0, scale = 0.0;
      for (size_t i = 0; i < again.size(); ++i) {
        difference = std::max(difference, std::abs(again[i] - solution[i]));
        scale = std::max(scale, std::abs(solution[i]));
      }
      MPI_Allreduce(MPI_IN_PLACE, &difference, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) {
        std::cout << std::setprecision(17) << "H2 logdet check: log|det| " << logabs_before << " (phase "
                  << phase_before << ") before and " << logabs_after << " (phase " << phase_after
                  << ") after another H2 was built, difference " << std::abs(logabs_after - logabs_before)
                  << "; the solve again: max difference " << difference / scale << " (relative)"
                  << std::setprecision(6) << std::endl;
      }
    }
    } else {
      partial_checks(driver_options, application, resources, rank, mpi_size, local_points, distributed_input);
      if (driver_options.partial_extract_bench > 0) partial_extract_bench(driver_options, resources, rank, mpi_size);
      partial_release_check(driver_options, resources, rank, local_points, distributed_input);
    }

    d_c_bpack_printstats(&resources.stats, &resources.process_tree);
  } catch (const std::exception& error) {
    std::cerr << "Laplace3D_H2Test_Driver error on rank " << rank << ": "
              << error.what() << std::endl;
    return_code = 1;
    MPI_Abort(MPI_COMM_WORLD, return_code);
  }

  if (source_comm != MPI_COMM_NULL) MPI_Comm_free(&source_comm);
  Cblacs_exit(1);
  MPI_Finalize();
  return return_code;
}
