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
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "dBPACK_wrapper.h"

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
  int h2_gemm_split = 16;
  int h2_ca_staged_halo = 0;
  int h2_ca_owner_component = 0;
  int h2_ca_owner_serial = 0;
  int precon = 1;
  int nrhs = 1;
  int verbosity = 1;
  int elem_extract = 2;
  int lrlevel = 0;
  int format = 7;
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
    } else if (name == "h2_gemm_split") {
      options.h2_gemm_split = parse_int(value, "H2_GEMM_split");
    } else if (name == "h2_ca_staged_halo") {
      options.h2_ca_staged_halo = parse_int(value, "H2_CA_staged_halo");
    } else if (name == "h2_ca_owner_component") {
      options.h2_ca_owner_component =
          parse_int(value, "H2_CA_owner_component");
    } else if (name == "h2_ca_owner_serial") {
      options.h2_ca_owner_serial = parse_int(value, "H2_CA_owner_serial");
    } else if (name == "precon") {
      options.precon = parse_int(value, "precon");
    } else if (name == "nrhs") {
      options.nrhs = parse_int(value, "nrhs");
    } else if (name == "verbosity") {
      options.verbosity = parse_int(value, "verbosity");
    } else if (name == "elem_extract") {
      options.elem_extract = parse_int(value, "elem_extract");
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
      if (parse_int(value, "sym") != 1) {
        throw std::invalid_argument("this driver requires --sym 1");
      }
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
  if (options.elem_extract != 0 && options.elem_extract != 2) {
    throw std::invalid_argument("elem_extract must be 0 or 2 for this driver");
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
  if (options.h2_lazy_schur < 0 || options.h2_lazy_schur > 2) {
    throw std::invalid_argument("H2_lazy_schur must be 0, 1, or 2");
  }
  if (options.h2_lazy_schur != 0 && options.h2_use_sketch != 2) {
    throw std::invalid_argument("H2_lazy_schur requires H2_use_sketch=2");
  }
  if (options.h2_gemm_split < 0) {
    throw std::invalid_argument("H2_GEMM_split must be nonnegative");
  }
  if (options.h2_ca_staged_halo != 0 && options.h2_ca_staged_halo != 2) {
    throw std::invalid_argument("H2_CA_staged_halo must be 0 or 2");
  }
  if (options.h2_ca_owner_component != 0 &&
      options.h2_ca_owner_component != 3) {
    throw std::invalid_argument("H2_CA_owner_component must be 0 or 3");
  }
  if (options.h2_ca_owner_serial != 0 && options.h2_ca_owner_serial != 1) {
    throw std::invalid_argument("H2_CA_owner_serial must be 0 or 1");
  }
  if (options.format != 1 && options.format != 7) {
    throw std::invalid_argument("format must be 1 (HODLR) or 7 (H2)");
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
      << "  --H2_GEMM_split <count>\n"
      << "  --H2_CA_staged_halo <0|2>\n"
      << "  --H2_CA_owner_component <0|3>\n"
      << "  --H2_CA_owner_serial <0|1>\n"
      << "  --precon <1|2|3>\n"
      << "  --nrhs <count>\n"
      << "  --format <1|7>\n"
      << "  --distributed64 <0|1>\n"
      << "  --elem_extract <0|2>\n"
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
                << "H2_GEMM_split: " << driver_options.h2_gemm_split << "\n"
                << "H2_CA_staged_halo: "
                << driver_options.h2_ca_staged_halo << "\n"
                << "H2_CA_owner_component: "
                << driver_options.h2_ca_owner_component << "\n"
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

    d_c_bpack_set_D_option(
        &resources.option, "tol_comp", driver_options.tolerance);
    d_c_bpack_set_I_option(
        &resources.option, "format", driver_options.format);
    d_c_bpack_set_I_option(&resources.option, "sym", 1);
    d_c_bpack_set_I_option(
        &resources.option, "Nmin_leaf", static_cast<int>(driver_options.nmin_leaf));
    d_c_bpack_set_I_option(
        &resources.option, "reduction_threshold",
        static_cast<int>(driver_options.reduction_threshold));
    d_c_bpack_set_I_option(
        &resources.option, "H2_unstructured",
        driver_options.h2_unstructured);
    d_c_bpack_set_I_option(
        &resources.option, "CA_level", driver_options.ca_level);
    d_c_bpack_set_I_option(
        &resources.option, "H2_use_sketch", driver_options.h2_use_sketch);
    d_c_bpack_set_I_option(
        &resources.option, "H2_lazy_schur", driver_options.h2_lazy_schur);
    d_c_bpack_set_I_option(
        &resources.option, "H2_GEMM_split", driver_options.h2_gemm_split);
    d_c_bpack_set_I_option(
        &resources.option, "H2_CA_staged_halo",
        driver_options.h2_ca_staged_halo);
    d_c_bpack_set_I_option(
        &resources.option, "H2_CA_owner_component",
        driver_options.h2_ca_owner_component);
    d_c_bpack_set_I_option(
        &resources.option, "H2_CA_owner_serial",
        driver_options.h2_ca_owner_serial);
    d_c_bpack_set_I_option(&resources.option, "precon", driver_options.precon);
    d_c_bpack_set_I_option(
        &resources.option, "verbosity", driver_options.verbosity);
    d_c_bpack_set_I_option(
        &resources.option, "elem_extract", driver_options.elem_extract);
    d_c_bpack_set_I_option(
        &resources.option, "LRlevel", driver_options.lrlevel);
    d_c_bpack_set_I_option(&resources.option, "cpp", 1);
    d_c_bpack_set_I_option(&resources.option, "nogeo", 0);

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
    if (driver_options.distributed64 == 1) {
      d_c_bpack_construct_element_compute_distributed64(
          &resources.matrix, &resources.option, &resources.stats,
          &resources.mesh, &kernel_quantities, &resources.process_tree,
          &laplace_entry_callback64, &laplace_block_callback64, &application);
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

    if (driver_options.format == 7) {
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
                    << "H2 solution check RHS " << column
                    << ": sum=" << solution_sum[index]
                    << ", norm=" << std::sqrt(solution_norm2[index])
                    << ", relative residual=" << relative_residual
                    << std::setprecision(6) << std::endl;
        }
      }
    }

    d_c_bpack_printstats(&resources.stats, &resources.process_tree);
  } catch (const std::exception& error) {
    std::cerr << "Laplace3D_H2Test_Driver error on rank " << rank << ": "
              << error.what() << std::endl;
    return_code = 1;
    MPI_Abort(MPI_COMM_WORLD, return_code);
  }

  Cblacs_exit(1);
  MPI_Finalize();
  return return_code;
}
