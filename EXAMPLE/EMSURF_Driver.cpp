/*
 * C++ EMSURF driver using ButterflyPACK's legacy 32-bit global-array API.
 * The mesh and RWG quadrature routines are shared with EMSURF_Module.f90
 * through the C bindings in EMSURF_C_Bindings.f90.
 */
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex.h>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "zBPACK_wrapper.h"

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

struct EMSurfOptions {
  std::string data_dir = "../EXAMPLE/EM3D_DATA/sphere_2300";
  double wavelength = 2.0;
  double frequency = 0.0;
  double cfie_alpha = 1.0;
  double scaling = 1.0;
  int use_frequency = 0;
  int rcs_static = 2;
  int rcs_nsample = 1000;
  int mesh_normal = 1;
};

double real_part(_Complex double value) { return __real__ value; }
double imag_part(_Complex double value) { return __imag__ value; }
double abs_squared(_Complex double value) {
  const double re = real_part(value);
  const double im = imag_part(value);
  return re * re + im * im;
}

double parse_double(const char* option, const char* value) {
  char* end = nullptr;
  const double parsed = std::strtod(value, &end);
  if (end == value || *end != '\0') {
    throw std::runtime_error(std::string("invalid value for ") + option + ": " + value);
  }
  return parsed;
}

int parse_int(const char* option, const char* value) {
  char* end = nullptr;
  const long parsed = std::strtol(value, &end, 10);
  if (end == value || *end != '\0') {
    throw std::runtime_error(std::string("invalid value for ") + option + ": " + value);
  }
  return static_cast<int>(parsed);
}

EMSurfOptions parse_emsurf_options(int argc, char** argv) {
  EMSurfOptions options;
  for (int arg = 1; arg < argc; ++arg) {
    const std::string name(argv[arg]);
    auto next = [&]() -> const char* {
      if (arg + 1 >= argc) throw std::runtime_error("missing value after " + name);
      return argv[++arg];
    };

    if (name == "--data_dir") {
      options.data_dir = next();
    } else if (name == "--wavelength") {
      options.wavelength = parse_double(name.c_str(), next());
      options.use_frequency = 0;
    } else if (name == "--freq") {
      options.frequency = parse_double(name.c_str(), next());
      options.use_frequency = 1;
    } else if (name == "--cfie_alpha") {
      options.cfie_alpha = parse_double(name.c_str(), next());
    } else if (name == "--rcs_static") {
      options.rcs_static = parse_int(name.c_str(), next());
    } else if (name == "--rcs_nsample") {
      options.rcs_nsample = parse_int(name.c_str(), next());
    } else if (name == "--mesh_normal") {
      options.mesh_normal = parse_int(name.c_str(), next());
    } else if (name == "--scaling") {
      options.scaling = parse_double(name.c_str(), next());
    }
  }
  if (options.rcs_nsample <= 0) throw std::runtime_error("--rcs_nsample must be positive");
  if (options.rcs_static != 1 && options.rcs_static != 2) {
    throw std::runtime_error("--rcs_static must be 1 (monostatic) or 2 (bistatic)");
  }
  if (options.cfie_alpha < 0.0 || options.cfie_alpha > 1.0) {
    throw std::runtime_error("--cfie_alpha must lie in [0,1]");
  }
  return options;
}

extern "C" {
void emsurf_initialize_c(const char*, double, double, int, double, int, int, int,
                         double, int);
void emsurf_get_problem_size_c(int*);
void emsurf_get_coordinates_c(double*);
void emsurf_get_minedgelength_c(double*);
void emsurf_get_wavenumber_c(double*);
void emsurf_entry_c(int*, int*, _Complex double*, C2Fptr);
void emsurf_block_c(int*, int*, int*, int64_t*, int*, int*, _Complex double*, int*,
                    int*, int*, int*, int*, C2Fptr);
void emsurf_incident_c(int, int, double, double, _Complex double*);
void emsurf_rcs_contribution_c(int, int, double, double, _Complex double,
                               _Complex double*);
void emsurf_finalize_c();

#ifdef HAVE_MPI
void Cblacs_exit(int);
#endif
}

void distance_callback(int*, int*, double* value, C2Fptr) { *value = 0.0; }
void near_far_callback(int*, int*, int* value, C2Fptr) { *value = 0; }

void print_symmetry_probe(int rank, int nunk) {
  if (rank != 0) return;
  const int count = std::min(nunk, 12);
  double maximum_relative_difference = 0.0;
  for (int sample = 1; sample <= count; ++sample) {
    int row = sample;
    int col = nunk - sample + 1;
    _Complex double a_row_col = 0.0;
    _Complex double a_col_row = 0.0;
    emsurf_entry_c(&row, &col, &a_row_col, nullptr);
    emsurf_entry_c(&col, &row, &a_col_row, nullptr);
    const double difference = std::sqrt(abs_squared(a_row_col - a_col_row));
    const double scale = std::max({1.0, std::sqrt(abs_squared(a_row_col)),
                                   std::sqrt(abs_squared(a_col_row))});
    maximum_relative_difference = std::max(maximum_relative_difference,
                                           difference / scale);
  }
  std::cout << std::scientific << std::setprecision(6)
            << "EMSURF complex-symmetry probe: max relative |A(i,j)-A(j,i)| = "
            << maximum_relative_difference << '\n';
}

void print_solution_checksums(int rank, int local_size, int nrhs,
                              const std::vector<_Complex double>& solution) {
  for (int rhs = 0; rhs < nrhs; ++rhs) {
    double local[3] = {0.0, 0.0, 0.0};
    for (int local_index = 0; local_index < local_size; ++local_index) {
      const auto value = solution[local_index + rhs * local_size];
      local[0] += real_part(value);
      local[1] += imag_part(value);
      local[2] += abs_squared(value);
    }
    double global[3] = {0.0, 0.0, 0.0};
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
      std::cout << std::scientific << std::setprecision(16)
                << "EMSURF solution checksum rhs " << rhs + 1 << ": sum = "
                << global[0] << " " << global[1]
                << ", norm2 = " << std::sqrt(global[2]) << '\n';
    }
  }
}

_Complex double far_field_sum(int local_size, F2Cptr* mesh,
                              const _Complex double* current, int polarization,
                              double theta, double phi) {
  double local_re = 0.0;
  double local_im = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : local_re, local_im)
#endif
  for (int local_index = 0; local_index < local_size; ++local_index) {
    int new_index = local_index + 1;
    int old_index = 0;
    z_c_bpack_new2old(mesh, &new_index, &old_index);
    _Complex double contribution = 0.0;
    emsurf_rcs_contribution_c(old_index, polarization, theta, phi,
                              current[local_index], &contribution);
    local_re += real_part(contribution);
    local_im += imag_part(contribution);
  }
  double local[2] = {local_re, local_im};
  double global[2] = {0.0, 0.0};
  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  _Complex double result = 0.0;
  __real__ result = global[0];
  __imag__ result = global[1];
  return result;
}

void write_bistatic_rcs(int rank, int local_size, F2Cptr* mesh, int samples,
                        double wavenumber,
                        const std::vector<_Complex double>& solution) {
  std::ofstream vv;
  std::ofstream hh;
  if (rank == 0) {
    vv.open("VV_bistatic.txt");
    hh.open("HH_bistatic.txt");
    vv << std::scientific << std::setprecision(16);
    hh << std::scientific << std::setprecision(16);
  }
  const double theta = 90.0;
  const double dphi = 180.0 / samples;
  for (int sample = 0; sample <= samples; ++sample) {
    const double phi = sample * dphi;
    const _Complex double vv_sum =
        far_field_sum(local_size, mesh, solution.data(), 0, theta, phi) * wavenumber;
    const _Complex double hh_sum =
        far_field_sum(local_size, mesh, solution.data() + local_size, 1, theta, phi) *
        wavenumber;
    if (rank == 0) {
      vv << phi << ' ' << real_part(vv_sum) << ' ' << imag_part(vv_sum) << '\n';
      hh << phi << ' ' << real_part(hh_sum) << ' ' << imag_part(hh_sum) << '\n';
    }
  }
}

void write_monostatic_rcs(int rank, int local_size, F2Cptr* mesh, int samples,
                          double wavenumber,
                          const std::vector<_Complex double>& solution) {
  std::ofstream output;
  if (rank == 0) {
    output.open("monostaticH.out");
    output << std::scientific << std::setprecision(16);
  }
  const double theta = 90.0;
  const double dphi = 180.0 / samples;
  for (int sample = 0; sample <= samples; ++sample) {
    const double phi = sample * dphi;
    const _Complex double field = far_field_sum(
        local_size, mesh, solution.data() + sample * local_size, 1, theta, phi);
    const double rcs = 10.0 * std::log10(wavenumber * wavenumber *
                                         abs_squared(field) / (4.0 * pi));
    if (rank == 0) output << phi << ' ' << rcs << '\n';
  }
}

}  // namespace

int main(int argc, char** argv) {
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  int rank = 0;
  int process_count = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);

  int return_code = 0;
  try {
    const EMSurfOptions em_options = parse_emsurf_options(argc, argv);
    const MPI_Fint fortran_communicator = MPI_Comm_c2f(MPI_COMM_WORLD);
    emsurf_initialize_c(em_options.data_dir.c_str(), em_options.wavelength,
                        em_options.frequency, em_options.use_frequency,
                        em_options.cfie_alpha, em_options.rcs_static,
                        em_options.rcs_nsample, em_options.mesh_normal,
                        em_options.scaling, fortran_communicator);

    int nunk = 0;
    emsurf_get_problem_size_c(&nunk);
    std::vector<double> coordinates(3 * static_cast<std::size_t>(nunk));
    emsurf_get_coordinates_c(coordinates.data());
    print_symmetry_probe(rank, nunk);

    std::vector<int> group_members(process_count);
    for (int process = 0; process < process_count; ++process) group_members[process] = process;

    F2Cptr matrix = nullptr;
    F2Cptr option = nullptr;
    F2Cptr statistics = nullptr;
    F2Cptr mesh = nullptr;
    F2Cptr kernel = nullptr;
    F2Cptr process_tree = nullptr;
    MPI_Fint mutable_communicator = fortran_communicator;
    z_c_bpack_createptree(&process_count, group_members.data(), &mutable_communicator,
                          &process_tree);
    z_c_bpack_createoption(&option);
    z_c_bpack_createstats(&statistics);

    z_c_bpack_set_I_option(&option, "ErrSol", 1);
    z_c_bpack_set_I_option(&option, "format", 1);
    z_c_bpack_set_D_option(&option, "near_para", 2.01);
    z_c_bpack_set_I_option(&option, "verbosity", 1);
    z_c_bpack_set_I_option(&option, "ILU", 0);
    z_c_bpack_set_I_option(&option, "forwardN15flag", 0);
    z_c_bpack_set_I_option(&option, "LRlevel", 100);
    z_c_bpack_set_D_option(&option, "tol_itersol", 1e-5);
    z_c_bpack_set_D_option(&option, "sample_para", 4.0);
    z_c_bpack_set_I_option(&option, "knn", 50);
    z_c_bpack_set_I_option(&option, "cpp", 1);
    z_c_bpack_set_option_from_command_line(argc, argv, option);
    double format_value = 0.0;
    z_c_bpack_getoption(&option, "format", &format_value);
    const int matrix_format = static_cast<int>(std::llround(format_value));

    double minimum_edge_length = 0.0;
    emsurf_get_minedgelength_c(&minimum_edge_length);
    z_c_bpack_set_D_option(&option, "touch_para", 3.0 * minimum_edge_length);
    z_c_bpack_printoption(&option, &process_tree);

    int dimension = 3;
    int nlevel = 0;
    int tree[1] = {nunk};
    int dummy_nearest_neighbor = 0;
    int local_size = 0;
    int context = 0;
    std::vector<int> permutation(nunk);
    z_c_bpack_construct_init(
        &nunk, &dimension, coordinates.data(), &dummy_nearest_neighbor, &nlevel, tree,
        permutation.data(), &local_size, &matrix, &option, &statistics, &mesh, &kernel,
        &process_tree, distance_callback, near_far_callback, &context);
    z_c_bpack_construct_element_compute(&matrix, &option, &statistics, &mesh, &kernel,
                                        &process_tree, emsurf_entry_c, emsurf_block_c,
                                        &context);
    z_c_bpack_factor(&matrix, &option, &statistics, &process_tree, &mesh);

    const int nrhs = em_options.rcs_static == 1 ? em_options.rcs_nsample + 1 : 2;
    std::vector<_Complex double> right_hand_side(
        static_cast<std::size_t>(local_size) * nrhs);
    std::vector<_Complex double> solution(
        static_cast<std::size_t>(local_size) * nrhs, 0.0);
    const double theta = 90.0;

    for (int rhs = 0; rhs < nrhs; ++rhs) {
      const int polarization = em_options.rcs_static == 1 ? 1 : rhs;
      const double phi = em_options.rcs_static == 1
                             ? rhs * (180.0 / em_options.rcs_nsample)
                             : 0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int local_index = 0; local_index < local_size; ++local_index) {
        int new_index = local_index + 1;
        int old_index = 0;
        z_c_bpack_new2old(&mesh, &new_index, &old_index);
        emsurf_incident_c(old_index, polarization, theta, phi,
                          &right_hand_side[local_index + rhs * local_size]);
      }
    }

    int mutable_nrhs = nrhs;
    z_c_bpack_solve(solution.data(), right_hand_side.data(), &local_size, &mutable_nrhs,
                    &matrix, &option, &statistics, &process_tree);
    print_solution_checksums(rank, local_size, nrhs, solution);

    double wavenumber = 0.0;
    emsurf_get_wavenumber_c(&wavenumber);
    if (em_options.rcs_static == 1) {
      write_monostatic_rcs(rank, local_size, &mesh, em_options.rcs_nsample,
                           wavenumber, solution);
    } else {
      write_bistatic_rcs(rank, local_size, &mesh, em_options.rcs_nsample,
                         wavenumber, solution);
    }

    if (rank == 0) std::cout << "Printing EMSURF ButterflyPACK statistics\n";
    z_c_bpack_printstats(&statistics, &process_tree);

    z_c_bpack_deletestats(&statistics);
    z_c_bpack_deletemesh(&mesh);
    if (matrix_format != 7 && kernel != nullptr) z_c_bpack_deletekernelquant(&kernel);
    z_c_bpack_delete(&matrix);
    z_c_bpack_deleteoption(&option);
    z_c_bpack_deleteproctree(&process_tree);
    emsurf_finalize_c();
  } catch (const std::exception& error) {
    if (rank == 0) std::cerr << "cie3d: " << error.what() << '\n';
    return_code = 1;
  }

#ifdef HAVE_MPI
  Cblacs_exit(1);
#endif
  MPI_Finalize();
  return return_code;
}
