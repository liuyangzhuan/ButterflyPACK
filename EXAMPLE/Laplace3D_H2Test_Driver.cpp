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
  int64_t reduction_threshold = 4096;
  int ca_level = 0;
  int precon = 1;
  int verbosity = 1;
  int elem_extract = 2;
  int lrlevel = 0;
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
    } else if (name == "ca_level") {
      options.ca_level = parse_int(value, "CA_level");
    } else if (name == "precon") {
      options.precon = parse_int(value, "precon");
    } else if (name == "verbosity") {
      options.verbosity = parse_int(value, "verbosity");
    } else if (name == "elem_extract") {
      options.elem_extract = parse_int(value, "elem_extract");
    } else if (name == "lrlevel") {
      options.lrlevel = parse_int(value, "lrlevel");
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
      if (parse_int(value, "format") != 7) {
        throw std::invalid_argument("this driver requires --format 7");
      }
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

  derive_nmin_leaf_from_levels(options);
  const int64_t points = checked_cube(options.grid_size, "grid_size");
  if (options.expected_points > 0 && options.expected_points != points) {
    throw std::invalid_argument(
        "N must equal grid_size^3; expected " + std::to_string(points));
  }
  if (points > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(
        "this ButterflyPACK C interface requires N <= INT_MAX");
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
      << "  Nmin_leaf=216, reduction_threshold=4096, CA_level=0\n\n"
      << "Options:\n"
      << "  --grid-size <n>\n"
      << "  --num-levels <count>\n"
      << "  --tol-comp <value>\n"
      << "  --Nmin_leaf <count>\n"
      << "  --reduction_threshold <count>\n"
      << "  --CA_level <level>\n"
      << "  --precon <1|2|3>\n"
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
  explicit Laplace3DApplication(int64_t grid_size)
      : grid_size_(grid_size),
        point_count_(checked_cube(grid_size, "grid_size")),
        inverse_4pi_n_(1.0 / (4.0 * kPi * static_cast<double>(point_count_))),
        diagonal_(laplace_self_cell_integral(grid_size)),
        locations_(static_cast<size_t>(point_count_) * 3) {
    const double h = 1.0 / static_cast<double>(grid_size_);
    for (int64_t i = 0; i < grid_size_; ++i) {
      const double x = (static_cast<double>(i) + 0.5) * h;
      for (int64_t j = 0; j < grid_size_; ++j) {
        const double y = (static_cast<double>(j) + 0.5) * h;
        for (int64_t k = 0; k < grid_size_; ++k) {
          const double z = (static_cast<double>(k) + 0.5) * h;
          const int64_t index = (i * grid_size_ + j) * grid_size_ + k;
          locations_[static_cast<size_t>(3 * index)] = x;
          locations_[static_cast<size_t>(3 * index + 1)] = y;
          locations_[static_cast<size_t>(3 * index + 2)] = z;
        }
      }
    }
  }

  double* locations() { return locations_.data(); }
  double diagonal() const { return diagonal_; }

  double evaluate(int row, int column) const {
    if (row == column) return diagonal_;

    const double* x = locations_.data() + static_cast<size_t>(3 * row);
    const double* y = locations_.data() + static_cast<size_t>(3 * column);
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    return inverse_4pi_n_ / r;
  }

  void evaluate_block(int rows, int columns, const int* row_indices,
                      const int* column_indices, double* output) const {
    for (int j = 0; j < columns; ++j) {
      const int column = column_indices[j] - 1;
      const double* y = locations_.data() + static_cast<size_t>(3 * column);
      for (int i = 0; i < rows; ++i) {
        const int row = row_indices[i] - 1;
        if (row == column) {
          output[i + static_cast<int64_t>(j) * rows] = diagonal_;
          continue;
        }

        const double* x = locations_.data() + static_cast<size_t>(3 * row);
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        const double dz = x[2] - y[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        output[i + static_cast<int64_t>(j) * rows] = inverse_4pi_n_ / r;
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
    int point_count = static_cast<int>(point_count_64);
    int dimension = 3;
    Laplace3DApplication application(driver_options.grid_size);

    if (rank == 0) {
      std::cout << "=== ButterflyPACK H2 3D Laplace Test ===\n"
                << "Grid size: " << driver_options.grid_size << " x "
                << driver_options.grid_size << " x "
                << driver_options.grid_size << "\n"
                << "Total points: " << point_count_64 << "\n"
                << "Tolerance: " << driver_options.tolerance << "\n"
                << "Nmin_leaf: " << driver_options.nmin_leaf << "\n"
                << "Reduction threshold: "
                << driver_options.reduction_threshold << "\n"
                << "CA_level: " << driver_options.ca_level << "\n"
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
    d_c_bpack_set_I_option(&resources.option, "format", 7);
    d_c_bpack_set_I_option(&resources.option, "sym", 1);
    d_c_bpack_set_I_option(
        &resources.option, "Nmin_leaf", static_cast<int>(driver_options.nmin_leaf));
    d_c_bpack_set_I_option(
        &resources.option, "reduction_threshold",
        static_cast<int>(driver_options.reduction_threshold));
    d_c_bpack_set_I_option(
        &resources.option, "CA_level", driver_options.ca_level);
    d_c_bpack_set_I_option(&resources.option, "precon", driver_options.precon);
    d_c_bpack_set_I_option(
        &resources.option, "verbosity", driver_options.verbosity);
    d_c_bpack_set_I_option(
        &resources.option, "elem_extract", driver_options.elem_extract);
    d_c_bpack_set_I_option(
        &resources.option, "LRlevel", driver_options.lrlevel);
    d_c_bpack_set_I_option(&resources.option, "cpp", 1);
    d_c_bpack_set_I_option(&resources.option, "nogeo", 0);

    int nlevel = 0;
    int user_tree = point_count;
    int local_points = 0;
    F2Cptr kernel_quantities = nullptr;
    d_c_bpack_construct_init(
        &point_count, &dimension, application.locations(), nullptr,
        &nlevel, &user_tree, nullptr, &local_points,
        &resources.matrix, &resources.option, &resources.stats,
        &resources.mesh, &kernel_quantities, &resources.process_tree,
        &distance_callback, &near_far_callback, &application);

    d_c_bpack_printoption(&resources.option, &resources.process_tree);
    d_c_bpack_construct_element_compute(
        &resources.matrix, &resources.option, &resources.stats,
        &resources.mesh, &kernel_quantities, &resources.process_tree,
        &laplace_entry_callback, &laplace_block_callback, &application);

    if (rank == 0) {
      std::cout << "\nFactoring the 3D Laplace operator:" << std::endl;
    }
    d_c_bpack_factor(
        &resources.matrix, &resources.option, &resources.stats,
        &resources.process_tree, &resources.mesh);
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
