#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
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
  int h2_ca_staged_halo = 0;
  int h2_ca_owner_component = 0;
  int h2_ca_pack_shared = 0;  // 0 auto, 1 always, 2 never
  int h2_ca_solve = 0;        // 0 replicated, 1 home
  int h2_ca_owner_serial = 0;
  int precon = 1;
  int nrhs = 1;
  int inv_mult_check = 0;  // columns of the inv_mult / transposed mult check (bpack_inv_mult_check.hpp); 0: none
  int extract_check = 0;   // intersections of the extraction check of a compression-only H2 (precon 2); 0: none
  int logdet_check = 0;    // 1: the logdet and a solve again after another H2 is built (format 7, not precon 2)
  int verbosity = 1;
  int elem_extract = -1;  // -1: 2, with hodlr_extract_ranks 1
  int lrlevel = 0;
  int hodlr_extract_ranks = 0;
  int gpu_evaluator = 1;
  int format = 7;
  int sym = 1;
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
  if (options.format == 7 && options.sym != 1) {
    throw std::invalid_argument("format=7 (H2) requires --sym 1");
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
      << "  --H2_CA_staged_halo <0|2>\n"
      << "  --H2_CA_owner_component <0|3>\n"
      << "  --H2_CA_pack_shared <auto|always|never>\n"
      << "  --H2_CA_solve <replicated|home>\n"
      << "  --H2_CA_owner_serial <0|1>\n"
      << "  --precon <1|2|3>\n"
      << "  --nrhs <count>\n"
      << "  --inv_mult_check <columns>   (after the solve: inv_mult N/T/C and the\n"
         "                                transposed mult; not with precon 2)\n"
      << "  --extract_check <count>   (format=7, precon 2, distributed64 1: that many\n"
         "                             random blocks extracted from the compressed H2,\n"
         "                             against the kernel and its matvec)\n"
      << "  --logdet_check <0|1>   (format=7, precon 1 or 3: the logdet and a solve of\n"
         "                          the factored H2 again after a compression-only H2\n"
         "                          of the same kernel is built)\n"
      << "  --format <1|7>\n"
      << "  --sym <0|1>          (0: unsymmetric HODLR, format=1 only)\n"
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
  Laplace3DApplication(int64_t grid_size, bool materialize_locations)
      : grid_size_(grid_size),
        point_count_(checked_cube(grid_size, "grid_size")),
        inverse_4pi_n_(1.0 / (4.0 * kPi * static_cast<double>(point_count_))),
        diagonal_(laplace_self_cell_integral(grid_size)) {
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
    if (row_id == column_id) return diagonal_;
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    return inverse_4pi_n_ / r;
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
  int64_t grid_size_;
  int64_t point_count_;
  double inverse_4pi_n_;
  double diagonal_;
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
        driver_options.grid_size, driver_options.distributed64 == 0);

    if (rank == 0) {
      std::cout << "=== ButterflyPACK H2 3D Laplace Test ===\n"
                << "Grid size: " << driver_options.grid_size << " x "
                << driver_options.grid_size << " x "
                << driver_options.grid_size << "\n"
                << "Total points: " << point_count_64 << "\n"
                << "Format: " << driver_options.format << "\n"
                << "sym: " << driver_options.sym << "\n"
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
      if (driver_options.format == 7) d_c_bpack_set_I_option(&source.option, "precon", 2);
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
        laplace3d_gpu_register(&source.matrix, application.inverse_4pi_n(),
                               application.diagonal());
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
        laplace3d_gpu_register(&resources.matrix, application.inverse_4pi_n(),
                               application.diagonal());
      }
    }
#endif
    if (driver_options.distributed64 == 1) {
      d_c_bpack_construct_element_compute_distributed64(
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &laplace_entry_callback64, &laplace_block_callback64, &application);
    } else if (source_ranks > 0) {
      // an H2 asks for the ID's sketches of its kernel rows sketched
      if (driver_options.format == 7) {
        d_c_bpack_set_sketch_callback(&resources.matrix, &extraction_sketch_callback);
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
      if (rank == 0) {
        std::cout << std::setprecision(4) << "H2 extraction check: " << ninter << " blocks, " << entries
                  << " entries in " << t_extract << " s; random blocks against the kernel "
                  << std::sqrt(kernel_error2 / kernel_norm2) << " (Frobenius, relative); "
                  << "rows of 2 columns against the matvec " << matvec_error / matvec_scale
                  << " (max, relative); the random blocks sketched (" << t_sketch
                  << " s) against their entries sketched " << sketch_error / sketch_scale << " (max, relative)"
                  << std::setprecision(6) << std::endl;
      }
    }

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

    if (driver_options.inv_mult_check > 0 && driver_options.precon != 2) {
      auto inv_mult = [&](char trans, const double* in, double* out, int ncol) {
        int n_in = local_points, n_out = local_points;
        d_c_bpack_inv_mult(&trans, in, out, &n_in, &n_out, &ncol,
                           &resources.matrix, &resources.option, &resources.stats,
                           &resources.process_tree);
      };
      auto mult = [&](char trans, const double* in, double* out, int ncol) {
        int n_in = local_points, n_out = local_points;
        d_c_bpack_mult(&trans, in, out, &n_in, &n_out, &ncol,
                       &resources.matrix, &resources.option, &resources.stats,
                       &resources.process_tree);
      };
      bpack_example::check_inv_mult<double>(
          MPI_COMM_WORLD, local_points, driver_options.inv_mult_check, 5,
          rhs, solution, number_of_rhs, inv_mult, mult,
          driver_options.format == 7 ? "H2" : "HODLR");
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
        laplace3d_gpu_register(&other.matrix, application.inverse_4pi_n(), application.diagonal());
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
