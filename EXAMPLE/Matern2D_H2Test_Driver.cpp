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
constexpr uint64_t kTrueSolutionSeed = 0x123456789abcdef0ULL;

struct DriverOptions {
  int64_t grid_size = 2048;
  int64_t expected_points = -1;
  int num_levels = 0;
  int64_t nmin_leaf = 64;
  bool nmin_leaf_set = false;
  double tolerance = 1e-9;
  int64_t reduction_threshold = 256;
  int ca_level = 0;
  int precon = 3;
  int iter_solver = 4;
  double cg_tolerance = 1e-10;
  int cg_max_iterations = 100;
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

int64_t checked_square(int64_t value, const char* name) {
  if (value <= 0 || value > std::numeric_limits<int64_t>::max() / value) {
    throw std::overflow_error(std::string(name) + " squared overflows int64_t");
  }
  return value * value;
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
  const int64_t derived_nmin_leaf =
      checked_square(leaf_grid_size, "leaf grid size");
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
    if (argument == "--cg") {
      options.iter_solver = 4;
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
    } else if (name == "iter_solver") {
      options.iter_solver = parse_int(value, "iter_solver");
    } else if (name == "cg_tol" || name == "tol_itersol") {
      options.cg_tolerance = parse_double(value, "cg_tol");
    } else if (name == "cg_max_iterations" || name == "n_iter") {
      options.cg_max_iterations = parse_int(value, "cg_max_iterations");
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
      if (parse_int(value, "dimension") != 2) {
        throw std::invalid_argument("this driver supports only --dimension 2");
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
  if (!(options.tolerance > 0.0) || !(options.cg_tolerance > 0.0)) {
    throw std::invalid_argument("tol_comp and cg_tol must be positive");
  }
  if (!(options.length_scale > 0.0) || options.nugget < 0.0) {
    throw std::invalid_argument(
        "length_scale must be positive and nugget must be nonnegative");
  }
  if (options.nmin_leaf <= 0 || options.reduction_threshold <= 0) {
    throw std::invalid_argument(
        "Nmin_leaf and reduction_threshold must be positive");
  }
  if (options.cg_max_iterations <= 0) {
    throw std::invalid_argument("cg_max_iterations must be positive");
  }
  if (options.elem_extract != 0 && options.elem_extract != 2) {
    throw std::invalid_argument("elem_extract must be 0 or 2 for this driver");
  }
  if (options.precon != 3) {
    throw std::invalid_argument(
        "this PCG driver requires precon=3 (factorized H2 preconditioner)");
  }
  // if (options.iter_solver != 4) {
  //   throw std::invalid_argument("this driver requires iter_solver=4 (CG)");
  // }

  derive_nmin_leaf_from_levels(options);
  const int64_t points = checked_square(options.grid_size, "grid_size");
  if (options.expected_points > 0 && options.expected_points != points) {
    throw std::invalid_argument(
        "N must equal grid_size^2; expected " + std::to_string(points));
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
      << "Defaults reproduce the active CA_matern.sh case:\n"
      << "  num_levels=9, N=4194304, grid_size=2048, tolerance=1e-9\n"
      << "  length_scale=0.1, nugget=1e-3, cg_tol=1e-10\n"
      << "  Nmin_leaf=64, reduction_threshold=256, CA_level=0\n\n"
      << "Options:\n"
      << "  --grid-size <n>\n"
      << "  --num-levels <count>\n"
      << "  --tol-comp <value>\n"
      << "  --Nmin_leaf <count>\n"
      << "  --reduction_threshold <count>\n"
      << "  --CA_level <level>\n"
      << "  --length-scale <value>\n"
      << "  --nugget <value>\n"
      << "  --cg-tol <value>\n"
      << "  --cg-max-iterations <count>\n"
      << "  --precon 3\n"
      << "  --iter_solver 4\n"
      << "  --elem_extract <0|2>\n"
      << "  --verbosity <-1|0|1>\n";
}

class Matern2DApplication {
 public:
  Matern2DApplication(int64_t grid_size, double length_scale, double nugget,
                      MPI_Comm communicator)
      : grid_size_(grid_size),
        point_count_(checked_square(grid_size, "grid_size")),
        length_scale_(length_scale),
        nugget_(nugget),
        communicator_(communicator),
        locations_(static_cast<size_t>(point_count_) * 2) {
    const double spacing = 1.0 / static_cast<double>(grid_size_);
    for (int64_t i = 0; i < grid_size_; ++i) {
      const double x = (static_cast<double>(i) + 0.5) * spacing;
      for (int64_t j = 0; j < grid_size_; ++j) {
        const double y = (static_cast<double>(j) + 0.5) * spacing;
        const int64_t index = i * grid_size_ + j;
        locations_[static_cast<size_t>(2 * index)] = x;
        locations_[static_cast<size_t>(2 * index + 1)] = y;
      }
    }
  }

  double* locations() { return locations_.data(); }
  int64_t point_count() const { return point_count_; }
  const std::vector<int>& local_old_indices() const {
    return local_old_indices_;
  }

  double evaluate(int row, int column) const {
    if (row == column) return 1.0 + nugget_;

    const double* x = locations_.data() + static_cast<size_t>(2 * row);
    const double* y = locations_.data() + static_cast<size_t>(2 * column);
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    return kernel_value(std::sqrt(dx * dx + dy * dy));
  }

  void evaluate_block(int rows, int columns, const int* row_indices,
                      const int* column_indices, double* output) const {
    for (int j = 0; j < columns; ++j) {
      const int column = column_indices[j] - 1;
      const double* y =
          locations_.data() + static_cast<size_t>(2 * column);
      for (int i = 0; i < rows; ++i) {
        const int row = row_indices[i] - 1;
        if (row == column) {
          output[i + static_cast<int64_t>(j) * rows] = 1.0 + nugget_;
          continue;
        }

        const double* x = locations_.data() + static_cast<size_t>(2 * row);
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        output[i + static_cast<int64_t>(j) * rows] =
            kernel_value(std::sqrt(dx * dx + dy * dy));
      }
    }
  }

  void configure_distribution(F2Cptr mesh, int local_points) {
    int local_start = 0;
    int reported_local_points = 0;
    int reported_global_points = 0;
    d_c_bpack_localindices(
        &mesh, &local_start, &reported_local_points, &reported_global_points);
    if (reported_local_points != local_points ||
        reported_global_points != point_count_) {
      throw std::runtime_error(
          "ButterflyPACK local-index metadata is inconsistent");
    }

    local_points_ = local_points;
    local_start_ = local_start;
    local_old_indices_.resize(static_cast<size_t>(local_points_));
    for (int i = 0; i < local_points_; ++i) {
      int new_local_index = i + 1;
      int old_global_index = 0;
      d_c_bpack_new2old(&mesh, &new_local_index, &old_global_index);
      if (old_global_index < 1 || old_global_index > point_count_) {
        throw std::runtime_error("ButterflyPACK permutation is invalid");
      }
      local_old_indices_[static_cast<size_t>(i)] = old_global_index - 1;
    }
  }

  void verify_matching_distribution(F2Cptr mesh, int local_points) const {
    int local_start = 0;
    int reported_local_points = 0;
    int reported_global_points = 0;
    d_c_bpack_localindices(
        &mesh, &local_start, &reported_local_points, &reported_global_points);

    int local_mismatch =
        reported_local_points != local_points ||
        reported_local_points != local_points_ ||
        reported_global_points != point_count_ || local_start != local_start_;
    if (local_mismatch == 0) {
      for (int i = 0; i < local_points_; ++i) {
        int new_local_index = i + 1;
        int old_global_index = 0;
        d_c_bpack_new2old(&mesh, &new_local_index, &old_global_index);
        if (old_global_index - 1 !=
            local_old_indices_[static_cast<size_t>(i)]) {
          local_mismatch = 1;
          break;
        }
      }
    }

    int global_mismatch = 0;
    MPI_Allreduce(&local_mismatch, &global_mismatch, 1, MPI_INT, MPI_MAX,
                  communicator_);
    if (global_mismatch != 0) {
      throw std::runtime_error(
          "operator and preconditioner H2 distributions do not match");
    }
  }

  void attach_operator(F2Cptr matrix, F2Cptr option, F2Cptr stats,
                       F2Cptr process_tree, int verbosity) {
    operator_matrix_ = matrix;
    operator_option_ = option;
    operator_stats_ = stats;
    operator_process_tree_ = process_tree;
    verbosity_ = verbosity;
  }

  void attach_preconditioner(F2Cptr matrix, F2Cptr option, F2Cptr stats,
                             F2Cptr process_tree) {
    preconditioner_matrix_ = matrix;
    preconditioner_option_ = option;
    preconditioner_stats_ = stats;
    preconditioner_process_tree_ = process_tree;
  }

  void apply_operator(const char* trans, int input_local_size,
                      int output_local_size, int number_of_vectors,
                      const double* input, double* output) {
    if (operator_matrix_ == nullptr || operator_option_ == nullptr ||
        operator_stats_ == nullptr || operator_process_tree_ == nullptr) {
      throw std::runtime_error("compression-only H2 operator is not attached");
    }
    if ((trans[0] != 'N' && trans[0] != 'n') ||
        input_local_size != local_points_ ||
        output_local_size != local_points_ || number_of_vectors != 1) {
      throw std::runtime_error(
          "the compression-only H2 operator requires one non-transposed RHS");
    }

    d_c_bpack_set_I_option(&operator_option_, "verbosity", -1);
    try {
      d_c_bpack_mult(trans, input, output, &input_local_size,
                     &output_local_size, &number_of_vectors,
                     &operator_matrix_, &operator_option_, &operator_stats_,
                     &operator_process_tree_);
    } catch (...) {
      d_c_bpack_set_I_option(&operator_option_, "verbosity", verbosity_);
      throw;
    }
    d_c_bpack_set_I_option(&operator_option_, "verbosity", verbosity_);
  }

  void apply_preconditioner(const char* trans, int input_local_size,
                            int output_local_size, int number_of_vectors,
                            const double* input, double* output) {
    if (preconditioner_matrix_ == nullptr ||
        preconditioner_option_ == nullptr || preconditioner_stats_ == nullptr ||
        preconditioner_process_tree_ == nullptr) {
      throw std::runtime_error("H2 preconditioner is not attached");
    }
    if ((trans[0] != 'N' && trans[0] != 'n') ||
        input_local_size != local_points_ ||
        output_local_size != local_points_ || number_of_vectors != 1) {
      throw std::runtime_error(
          "the H2 PCG preconditioner requires one non-transposed RHS");
    }

    // c_bpack_inv_mult does not dispatch format 7. Use the H2-aware solve
    // entry point and silence only its nested per-iteration banners.
    d_c_bpack_set_I_option(&preconditioner_option_, "verbosity", -1);
    try {
      d_c_bpack_solve(
          output, const_cast<double*>(input), &input_local_size,
          &number_of_vectors, &preconditioner_matrix_, &preconditioner_option_,
          &preconditioner_stats_, &preconditioner_process_tree_);
    } catch (...) {
      d_c_bpack_set_I_option(
          &preconditioner_option_, "verbosity", verbosity_);
      throw;
    }
    d_c_bpack_set_I_option(
        &preconditioner_option_, "verbosity", verbosity_);
  }

 private:
  double kernel_value(double distance) const {
    const double scaled = kSqrt5 * distance / length_scale_;
    return (1.0 + scaled + scaled * scaled / 3.0) * std::exp(-scaled);
  }

  int64_t grid_size_;
  int64_t point_count_;
  double length_scale_;
  double nugget_;
  MPI_Comm communicator_;
  int local_points_ = 0;
  int local_start_ = 0;
  std::vector<double> locations_;
  std::vector<int> local_old_indices_;
  F2Cptr operator_matrix_ = nullptr;
  F2Cptr operator_option_ = nullptr;
  F2Cptr operator_stats_ = nullptr;
  F2Cptr operator_process_tree_ = nullptr;
  F2Cptr preconditioner_matrix_ = nullptr;
  F2Cptr preconditioner_option_ = nullptr;
  F2Cptr preconditioner_stats_ = nullptr;
  F2Cptr preconditioner_process_tree_ = nullptr;
  int verbosity_ = 1;
};

void matern_entry_callback(int* row, int* column, double* value,
                           C2Fptr quant) {
  const auto* application = static_cast<Matern2DApplication*>(quant);
  *value = application->evaluate(*row - 1, *column - 1);
}

void matern_block_callback(
    int* ninter, int* nallrows, int* nallcols, int64_t* nalldat_loc,
    int* allrows, int* allcols, double* alldat_loc, int* rowidx,
    int* colidx, int* pgidx, int* npmap, int* pmaps, C2Fptr quant) {
  const auto* application = static_cast<Matern2DApplication*>(quant);
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

void h2_operator_callback(const char* trans, int* input_local_size,
                          int* output_local_size, int* number_of_vectors,
                          const double* input, double* output, C2Fptr quant) {
  auto* application = static_cast<Matern2DApplication*>(quant);
  try {
    application->apply_operator(
        trans, *input_local_size, *output_local_size, *number_of_vectors,
        input, output);
  } catch (const std::exception& error) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cerr << "compression-only H2 Matvec callback failed on rank " << rank
              << ": " << error.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
}

void h2_preconditioner_callback(
    const char* trans, int* input_local_size, int* output_local_size,
    int* number_of_vectors, const double* input, double* output,
    C2Fptr quant) {
  auto* application = static_cast<Matern2DApplication*>(quant);
  try {
    application->apply_preconditioner(
        trans, *input_local_size, *output_local_size, *number_of_vectors,
        input, output);
  } catch (const std::exception& error) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cerr << "H2 preconditioner callback failed on rank " << rank << ": "
              << error.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 4);
  }
}

void distance_callback(int*, int*, double* value, C2Fptr) { *value = 0.0; }
void near_far_callback(int*, int*, int* value, C2Fptr) { *value = 0; }

uint64_t splitmix64(uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

double centered_uniform_from_hash(uint64_t key) {
  constexpr long double scale =
      1.0L / static_cast<long double>(std::numeric_limits<uint64_t>::max());
  const long double uniform =
      static_cast<long double>(splitmix64(key)) * scale;
  return static_cast<double>(uniform - 0.5L);
}

double make_true_solution_entry(int global_index) {
  const uint64_t base =
      splitmix64(kTrueSolutionSeed ^ static_cast<uint64_t>(global_index));
  return centered_uniform_from_hash(base);
}

double relative_forward_error(const std::vector<double>& solution,
                              const std::vector<double>& true_solution,
                              MPI_Comm communicator) {
  if (solution.size() != true_solution.size()) {
    throw std::invalid_argument("forward-error vectors have different sizes");
  }

  double local_values[2] = {0.0, 0.0};
  for (size_t i = 0; i < solution.size(); ++i) {
    const double difference = solution[i] - true_solution[i];
    local_values[0] += difference * difference;
    local_values[1] += true_solution[i] * true_solution[i];
  }

  double global_values[2] = {0.0, 0.0};
  MPI_Allreduce(local_values, global_values, 2, MPI_DOUBLE, MPI_SUM,
                communicator);
  if (global_values[1] == 0.0) {
    throw std::runtime_error("x_true has zero norm");
  }
  return std::sqrt(global_values[0] / global_values[1]);
}

struct ButterflyResources {
  F2Cptr option = nullptr;
  F2Cptr stats = nullptr;
  F2Cptr process_tree = nullptr;
  F2Cptr matrix = nullptr;
  F2Cptr mesh = nullptr;
  F2Cptr construction_kernel = nullptr;
  F2Cptr iterative_kernel = nullptr;

  ~ButterflyResources() {
    if (stats != nullptr) d_c_bpack_deletestats(&stats);
    if (process_tree != nullptr) d_c_bpack_deleteproctree(&process_tree);
    if (mesh != nullptr) d_c_bpack_deletemesh(&mesh);
    if (iterative_kernel != nullptr) {
      d_c_bpack_deletekernelquant(&iterative_kernel);
    }
    if (construction_kernel != nullptr) {
      d_c_bpack_deletekernelquant(&construction_kernel);
    }
    if (matrix != nullptr) d_c_bpack_delete(&matrix);
    if (option != nullptr) d_c_bpack_deleteoption(&option);
  }
};

void initialize_h2_resources(ButterflyResources& resources,
                             const DriverOptions& driver_options,
                             double compression_tolerance, int precon,
                             MPI_Comm communicator) {
  d_c_bpack_createoption(&resources.option);
  d_c_bpack_createstats(&resources.stats);

  int mpi_size = 1;
  MPI_Comm_size(communicator, &mpi_size);
  std::vector<int> groups(static_cast<size_t>(mpi_size));
  std::iota(groups.begin(), groups.end(), 0);
  MPI_Fint fortran_comm = MPI_Comm_c2f(communicator);
  d_c_bpack_createptree(
      &mpi_size, groups.data(), &fortran_comm, &resources.process_tree);

  d_c_bpack_set_D_option(
      &resources.option, "tol_comp", compression_tolerance);
  d_c_bpack_set_D_option(
      &resources.option, "tol_itersol", driver_options.cg_tolerance);
  d_c_bpack_set_I_option(
      &resources.option, "n_iter", driver_options.cg_max_iterations);
  d_c_bpack_set_I_option(
      &resources.option, "iter_solver", driver_options.iter_solver);
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
  d_c_bpack_set_I_option(&resources.option, "precon", precon);
  d_c_bpack_set_I_option(
      &resources.option, "verbosity", driver_options.verbosity);
  d_c_bpack_set_I_option(
      &resources.option, "elem_extract", driver_options.elem_extract);
  d_c_bpack_set_I_option(
      &resources.option, "LRlevel", driver_options.lrlevel);
  d_c_bpack_set_I_option(&resources.option, "cpp", 1);
  d_c_bpack_set_I_option(&resources.option, "nogeo", 0);
}

}  // namespace

extern "C" void Cblacs_exit(int);

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int return_code = 0;
  try {
    const DriverOptions driver_options = parse_driver_options(argc, argv);
    if (driver_options.show_help) {
      if (rank == 0) print_usage(argv[0]);
      MPI_Finalize();
      return 0;
    }

    const int64_t point_count_64 =
        checked_square(driver_options.grid_size, "grid_size");
    const double operator_tolerance = 1e-13; //driver_options.tolerance * 0.01;
    if (!(operator_tolerance > 0.0) ||
        !std::isfinite(operator_tolerance)) {
      throw std::invalid_argument("0.01*tol_comp is outside the double range");
    }
    int point_count = static_cast<int>(point_count_64);
    int dimension = 2;
    Matern2DApplication application(
        driver_options.grid_size, driver_options.length_scale,
        driver_options.nugget, MPI_COMM_WORLD);

    if (rank == 0) {
      std::cout << "=== ButterflyPACK H2 2D Matern 5/2 PCG Test ===\n"
                << "Grid size: " << driver_options.grid_size << " x "
                << driver_options.grid_size << "\n"
                << "Total points: " << point_count_64 << "\n"
                << "H2 operator tolerance: " << operator_tolerance << "\n"
                << "H2 preconditioner tolerance: "
                << driver_options.tolerance << "\n"
                << "Length scale: " << driver_options.length_scale << "\n"
                << "Nugget: " << driver_options.nugget << "\n"
                << "Nmin_leaf: " << driver_options.nmin_leaf << "\n"
                << "Reduction threshold: "
                << driver_options.reduction_threshold << "\n"
                << "CA_level: " << driver_options.ca_level << "\n"
                << "CG tolerance: " << driver_options.cg_tolerance << "\n"
                << "CG maximum iterations: "
                << driver_options.cg_max_iterations << "\n"
                << "Operator: compression-only H2 matvec (precon=2)\n"
                << "Preconditioner: factorized H2 inverse (precon=3)\n"
                << "Proxy points: 0 (current ButterflyPACK H2 interface)"
                << std::endl;
    }

    ButterflyResources operator_resources;
    ButterflyResources preconditioner_resources;
    initialize_h2_resources(
        operator_resources, driver_options, operator_tolerance, 2,
        MPI_COMM_WORLD);
    initialize_h2_resources(
        preconditioner_resources, driver_options, driver_options.tolerance,
        driver_options.precon, MPI_COMM_WORLD);

    int operator_nlevel = 0;
    int operator_user_tree = point_count;
    int operator_local_points = 0;
    d_c_bpack_construct_init(
        &point_count, &dimension, application.locations(), nullptr,
        &operator_nlevel, &operator_user_tree, nullptr, &operator_local_points,
        &operator_resources.matrix, &operator_resources.option,
        &operator_resources.stats, &operator_resources.mesh,
        &operator_resources.construction_kernel,
        &operator_resources.process_tree,
        &distance_callback, &near_far_callback, &application);

    application.configure_distribution(
        operator_resources.mesh, operator_local_points);
    application.attach_operator(
        operator_resources.matrix, operator_resources.option,
        operator_resources.stats, operator_resources.process_tree,
        driver_options.verbosity);
    if (rank == 0) {
      std::cout << "\nBuilding compression-only H2 operator at tolerance "
                << std::scientific << operator_tolerance << std::defaultfloat
                << ":" << std::endl;
    }
    d_c_bpack_printoption(
        &operator_resources.option, &operator_resources.process_tree);
    d_c_bpack_construct_element_compute(
        &operator_resources.matrix, &operator_resources.option,
        &operator_resources.stats, &operator_resources.mesh,
        &operator_resources.construction_kernel,
        &operator_resources.process_tree, &matern_entry_callback,
        &matern_block_callback, &application);
    d_c_bpack_factor(
        &operator_resources.matrix, &operator_resources.option,
        &operator_resources.stats, &operator_resources.process_tree,
        &operator_resources.mesh);

    int preconditioner_nlevel = 0;
    int preconditioner_user_tree = point_count;
    int preconditioner_local_points = 0;
    d_c_bpack_construct_init(
        &point_count, &dimension, application.locations(), nullptr,
        &preconditioner_nlevel, &preconditioner_user_tree, nullptr,
        &preconditioner_local_points, &preconditioner_resources.matrix,
        &preconditioner_resources.option, &preconditioner_resources.stats,
        &preconditioner_resources.mesh,
        &preconditioner_resources.construction_kernel,
        &preconditioner_resources.process_tree, &distance_callback,
        &near_far_callback, &application);

    application.verify_matching_distribution(
        preconditioner_resources.mesh, preconditioner_local_points);
    application.attach_preconditioner(
        preconditioner_resources.matrix, preconditioner_resources.option,
        preconditioner_resources.stats,
        preconditioner_resources.process_tree);
    if (rank == 0) {
      std::cout << "\nBuilding factorized H2 preconditioner at tolerance "
                << std::scientific << driver_options.tolerance
                << std::defaultfloat << ":" << std::endl;
    }
    d_c_bpack_printoption(
        &preconditioner_resources.option,
        &preconditioner_resources.process_tree);
    d_c_bpack_construct_element_compute(
        &preconditioner_resources.matrix, &preconditioner_resources.option,
        &preconditioner_resources.stats, &preconditioner_resources.mesh,
        &preconditioner_resources.construction_kernel,
        &preconditioner_resources.process_tree, &matern_entry_callback,
        &matern_block_callback, &application);
    d_c_bpack_factor(
        &preconditioner_resources.matrix, &preconditioner_resources.option,
        &preconditioner_resources.stats,
        &preconditioner_resources.process_tree,
        &preconditioner_resources.mesh);

    int local_points = operator_local_points;

    std::vector<double> true_solution(static_cast<size_t>(local_points));
    std::vector<double> rhs(static_cast<size_t>(local_points));
    std::vector<double> solution(static_cast<size_t>(local_points));
    const std::vector<int>& local_old_indices =
        application.local_old_indices();
    for (int i = 0; i < local_points; ++i) {
      true_solution[static_cast<size_t>(i)] =
          make_true_solution_entry(local_old_indices[static_cast<size_t>(i)]);
    }
    application.apply_operator(
        "N", local_points, local_points, 1, true_solution.data(), rhs.data());

    int number_of_rhs = 1;
    application.apply_preconditioner(
        "N", local_points, local_points, number_of_rhs, rhs.data(),
        solution.data());
    const double forward_error_before_cg =
        relative_forward_error(solution, true_solution, MPI_COMM_WORLD);
    if (rank == 0) {
      std::cout << "Forward error before CG "
                   "(x=H2_preconditioner^{-1}*b): "
                << std::scientific << std::setprecision(16)
                << forward_error_before_cg << std::defaultfloat << std::endl;
    }

    if (rank == 0) {
      std::cout << "\nPCG: b=H2_operator*x_true, initialized with "
                   "x=H2_preconditioner^{-1}*b"
                << std::endl;
    }
    d_c_bpack_iter_usermatvec_precon(
        solution.data(), rhs.data(), &local_points, &number_of_rhs,
        &preconditioner_resources.option, &preconditioner_resources.stats,
        &preconditioner_resources.process_tree,
        &preconditioner_resources.iterative_kernel, &h2_operator_callback,
        &h2_preconditioner_callback, &application);

    const double forward_error_after_cg =
        relative_forward_error(solution, true_solution, MPI_COMM_WORLD);
    if (rank == 0) {
      std::cout << "Forward error after CG: " << std::scientific
                << std::setprecision(16) << forward_error_after_cg
                << std::defaultfloat << std::endl;
    }

    std::vector<double> operator_product(static_cast<size_t>(local_points));
    application.apply_operator(
        "N", local_points, local_points, 1, solution.data(),
        operator_product.data());
    double local_residual_squared = 0.0;
    double local_rhs_squared = 0.0;
    for (int i = 0; i < local_points; ++i) {
      const double residual =
          operator_product[static_cast<size_t>(i)] - rhs[static_cast<size_t>(i)];
      local_residual_squared += residual * residual;
      local_rhs_squared +=
          rhs[static_cast<size_t>(i)] * rhs[static_cast<size_t>(i)];
    }
    double global_residual_squared = 0.0;
    double global_rhs_squared = 0.0;
    MPI_Allreduce(&local_residual_squared, &global_residual_squared, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_rhs_squared, &global_rhs_squared, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    const double operator_relative_residual =
        std::sqrt(global_residual_squared / global_rhs_squared);
    if (rank == 0) {
      std::cout << "CG ||H2_operator*x-b||_2/||b||_2: " << std::scientific
                << std::setprecision(16) << operator_relative_residual
                << std::defaultfloat << std::endl;
    }

    if (rank == 0) {
      std::cout << "\nCompression-only H2 operator statistics:" << std::endl;
    }
    d_c_bpack_printstats(
        &operator_resources.stats, &operator_resources.process_tree);
    if (rank == 0) {
      std::cout << "\nFactorized H2 preconditioner statistics:" << std::endl;
    }
    d_c_bpack_printstats(
        &preconditioner_resources.stats,
        &preconditioner_resources.process_tree);
  } catch (const std::exception& error) {
    std::cerr << "Matern2D_H2Test_Driver error on rank " << rank << ": "
              << error.what() << std::endl;
    return_code = 1;
    MPI_Abort(MPI_COMM_WORLD, return_code);
  }

  Cblacs_exit(1);
  MPI_Finalize();
  return return_code;
}
