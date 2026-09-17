#include <algorithm>
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

constexpr double kSqrt5 = 2.2360679774997896964091736687747632;

struct DriverOptions {
  int64_t grid_size = 65536;
  int64_t expected_points = -1;
  int num_levels = 0;
  int64_t nmin_leaf = 64;
  bool nmin_leaf_set = false;
  double tolerance = 1e-8;
  int64_t reduction_threshold = 2;
  int ca_level = 10000;
  int precon = 1;
  int nrhs = 1;
  double length_scale = 0.1;
  double nugget = 1e-3;
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

void derive_nmin_leaf_from_levels(DriverOptions& options) {
  if (options.num_levels == 0) return;
  if (options.num_levels < 2 || options.num_levels >= 63) {
    throw std::invalid_argument("num_levels must be between 2 and 62");
  }

  const int64_t leaf_boxes = int64_t{1} << (options.num_levels - 1);
  if (options.grid_size % leaf_boxes != 0) {
    throw std::invalid_argument(
        "grid_size must be divisible by 2^(num_levels-1)");
  }

  const int64_t derived_nmin_leaf = options.grid_size / leaf_boxes;
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
    if (equals != std::string::npos) value = argument.substr(equals + 1);
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
    } else if (name == "nrhs") {
      options.nrhs = parse_int(value, "nrhs");
    } else if (name == "length_scale") {
      options.length_scale = parse_double(value, "length_scale");
    } else if (name == "nugget") {
      options.nugget = parse_double(value, "nugget");
    } else if (name == "verbosity") {
      options.verbosity = parse_int(value, "verbosity");
    } else if (name == "elem_extract") {
      options.elem_extract = parse_int(value, "elem_extract");
    } else if (name == "lrlevel") {
      options.lrlevel = parse_int(value, "lrlevel");
    } else if (name == "kernel") {
      const std::string kernel = normalize_option_name(value);
      if (kernel != "matern52" && kernel != "matern5/2" &&
          kernel != "matern_52") {
        throw std::invalid_argument(
            "this driver supports only --kernel matern52");
      }
    } else if (name == "number_type" || name == "number") {
      if (normalize_option_name(value) != "real") {
        throw std::invalid_argument(
            "this driver supports only --number-type real");
      }
    } else if (name == "dimension") {
      if (parse_int(value, "dimension") != 1) {
        throw std::invalid_argument("this driver supports only --dimension 1");
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
  derive_nmin_leaf_from_levels(options);
  if (!(options.tolerance > 0.0)) {
    throw std::invalid_argument("tol_comp must be positive");
  }
  if (!(options.length_scale > 0.0) || options.nugget < 0.0) {
    throw std::invalid_argument(
        "length_scale must be positive and nugget must be nonnegative");
  }
  if (options.nmin_leaf < 2) {
    throw std::invalid_argument("Nmin_leaf must be at least 2");
  }
  if (options.reduction_threshold < 2) {
    throw std::invalid_argument(
        "reduction_threshold must be at least 2 for a 1D H2 tree");
  }
  if (options.precon != 1 && options.precon != 3) {
    throw std::invalid_argument(
        "this driver supports factorized modes precon=1 or precon=3");
  }
  if (options.nrhs <= 0) {
    throw std::invalid_argument("nrhs must be positive");
  }
  if (options.elem_extract != 0 && options.elem_extract != 2) {
    throw std::invalid_argument("elem_extract must be 0 or 2 for this driver");
  }

  const int64_t points = options.grid_size;
  if (options.expected_points > 0 && options.expected_points != points) {
    throw std::invalid_argument(
        "N must equal grid_size for the 1D driver; expected " +
        std::to_string(points));
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
      << "Defaults:\n"
      << "  N=grid_size=65536, tolerance=1e-8, Nmin_leaf=64\n"
      << "  reduction_threshold=2, CA_level=10000 (fully color)\n"
      << "  length_scale=0.1, nugget=1e-3, precon=1\n\n"
      << "Options:\n"
      << "  --grid-size <n>\n"
      << "  --num-levels <count>\n"
      << "  --tol-comp <value>\n"
      << "  --Nmin_leaf <count>\n"
      << "  --reduction_threshold <count>\n"
      << "  --CA_level <level>  (accepted for compatibility; 1D stays color)\n"
      << "  --length-scale <value>\n"
      << "  --nugget <value>\n"
      << "  --precon <1|3>\n"
      << "  --nrhs <count>\n"
      << "  --elem_extract <0|2>\n"
      << "  --verbosity <-1|0|1>\n";
}

class Matern1DApplication {
 public:
  Matern1DApplication(int64_t point_count, double length_scale, double nugget)
      : point_count_(point_count),
        length_scale_(length_scale),
        nugget_(nugget),
        locations_(static_cast<size_t>(point_count)) {
    const double spacing = 1.0 / static_cast<double>(point_count_);
    for (int64_t i = 0; i < point_count_; ++i) {
      locations_[static_cast<size_t>(i)] =
          (static_cast<double>(i) + 0.5) * spacing;
    }
  }

  double* locations() { return locations_.data(); }

  double evaluate(int row, int column) const {
    if (row == column) return 1.0 + nugget_;
    const double distance = std::abs(
        locations_[static_cast<size_t>(row)] -
        locations_[static_cast<size_t>(column)]);
    return kernel_value(distance);
  }

  void evaluate_block(int rows, int columns, const int* row_indices,
                      const int* column_indices, double* output) const {
    for (int j = 0; j < columns; ++j) {
      const int column = column_indices[j] - 1;
      for (int i = 0; i < rows; ++i) {
        const int row = row_indices[i] - 1;
        output[i + static_cast<int64_t>(j) * rows] = evaluate(row, column);
      }
    }
  }

 private:
  double kernel_value(double distance) const {
    const double scaled = kSqrt5 * distance / length_scale_;
    return (1.0 + scaled + scaled * scaled / 3.0) * std::exp(-scaled);
  }

  int64_t point_count_;
  double length_scale_;
  double nugget_;
  std::vector<double> locations_;
};

void matern_entry_callback(int* row, int* column, double* value,
                           C2Fptr quant) {
  const auto* application = static_cast<Matern1DApplication*>(quant);
  *value = application->evaluate(*row - 1, *column - 1);
}

void matern_block_callback(
    int* ninter, int* nallrows, int* nallcols, int64_t* nalldat_loc,
    int* allrows, int* allcols, double* alldat_loc, int* rowidx,
    int* colidx, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  const auto* application = static_cast<Matern1DApplication*>(quant);
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
        std::cerr << "matern_block_callback supports only single-process blocks"
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
  F2Cptr kernel_quantities = nullptr;

  ~ButterflyResources() {
    if (stats != nullptr) d_c_bpack_deletestats(&stats);
    if (process_tree != nullptr) d_c_bpack_deleteproctree(&process_tree);
    if (mesh != nullptr) d_c_bpack_deletemesh(&mesh);
    if (kernel_quantities != nullptr) {
      d_c_bpack_deletekernelquant(&kernel_quantities);
    }
    if (matrix != nullptr) d_c_bpack_delete(&matrix);
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

    const int64_t point_count_64 = driver_options.grid_size;
    int point_count = static_cast<int>(point_count_64);
    int dimension = 1;
    Matern1DApplication application(
        point_count_64, driver_options.length_scale, driver_options.nugget);

    if (rank == 0) {
      std::cout << "=== ButterflyPACK H2 1D Matern 5/2 Test ===\n"
                << "Grid size: " << driver_options.grid_size << "\n"
                << "Total points: " << point_count_64 << "\n"
                << "Tolerance: " << driver_options.tolerance << "\n"
                << "Length scale: " << driver_options.length_scale << "\n"
                << "Nugget: " << driver_options.nugget << "\n"
                << "Nmin_leaf: " << driver_options.nmin_leaf << "\n"
                << "Reduction threshold: "
                << driver_options.reduction_threshold << "\n"
                << "Requested CA_level: " << driver_options.ca_level
                << " (1D H2 is forced to fully color)\n"
                << "Preconditioner mode: " << driver_options.precon << "\n"
                << "Number of RHS: " << driver_options.nrhs << "\n"
                << "Proxy points: 0 (current ButterflyPACK H2 interface)"
                << std::endl;
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
        &resources.option, "Nmin_leaf",
        static_cast<int>(driver_options.nmin_leaf));
    d_c_bpack_set_I_option(
        &resources.option, "reduction_threshold",
        static_cast<int>(driver_options.reduction_threshold));
    d_c_bpack_set_I_option(
        &resources.option, "CA_level", driver_options.ca_level);
    d_c_bpack_set_I_option(
        &resources.option, "precon", driver_options.precon);
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
    d_c_bpack_construct_init(
        &point_count, &dimension, application.locations(), nullptr,
        &nlevel, &user_tree, nullptr, &local_points,
        &resources.matrix, &resources.option, &resources.stats,
        &resources.mesh, &resources.kernel_quantities,
        &resources.process_tree, &distance_callback, &near_far_callback,
        &application);

    d_c_bpack_printoption(&resources.option, &resources.process_tree);
    d_c_bpack_construct_element_compute(
        &resources.matrix, &resources.option, &resources.stats,
        &resources.mesh, &resources.kernel_quantities,
        &resources.process_tree, &matern_entry_callback,
        &matern_block_callback, &application);

    if (rank == 0) {
      std::cout << "\nFactoring the 1D Matern operator:" << std::endl;
    }
    d_c_bpack_factor(
        &resources.matrix, &resources.option, &resources.stats,
        &resources.process_tree, &resources.mesh);

    const int number_of_rhs = driver_options.nrhs;
    const size_t value_count = static_cast<size_t>(local_points) *
        static_cast<size_t>(number_of_rhs);
    std::vector<int> old_global_indices(static_cast<size_t>(local_points));
    for (int row = 0; row < local_points; ++row) {
      int new_local_index = row + 1;
      int old_global_index = 0;
      d_c_bpack_new2old(
          &resources.mesh, &new_local_index, &old_global_index);
      old_global_indices[static_cast<size_t>(row)] = old_global_index - 1;
    }

    std::vector<double> rhs(value_count);
    for (int column = 0; column < number_of_rhs; ++column) {
      for (int row = 0; row < local_points; ++row) {
        const int old_index = old_global_indices[static_cast<size_t>(row)];
        rhs[static_cast<size_t>(row) +
            static_cast<size_t>(column) * local_points] =
            1.0 + 0.125 * column + 0.001 * (old_index % 97);
      }
    }

    if (rank == 0) {
      std::cout << "\nSolving the 1D Matern system:" << std::endl;
    }
    std::vector<double> solution(value_count, 0.0);
    int solve_rhs = number_of_rhs;
    d_c_bpack_solve(
        solution.data(), rhs.data(), &local_points, &solve_rhs,
        &resources.matrix, &resources.option, &resources.stats,
        &resources.process_tree);

    int local_nonfinite = std::any_of(
        solution.begin(), solution.end(),
        [](double value) { return !std::isfinite(value); }) ? 1 : 0;
    int global_nonfinite = 0;
    MPI_Allreduce(
        &local_nonfinite, &global_nonfinite, 1, MPI_INT, MPI_MAX,
        MPI_COMM_WORLD);
    if (global_nonfinite != 0) {
      throw std::runtime_error("the H2 solve returned a non-finite value");
    }

    d_c_bpack_printstats(&resources.stats, &resources.process_tree);
  } catch (const std::exception& error) {
    std::cerr << "Matern1D_H2Test_Driver error on rank " << rank << ": "
              << error.what() << std::endl;
    return_code = 1;
    MPI_Abort(MPI_COMM_WORLD, return_code);
  }

  Cblacs_exit(1);
  MPI_Finalize();
  return return_code;
}
