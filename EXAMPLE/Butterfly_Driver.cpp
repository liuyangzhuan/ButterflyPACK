/* “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
  Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
  U.S. Dept. of Energy). All rights reserved.

  If you have questions about your rights to use or distribute this software, please contact
  Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

  NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
  U.S. Government consequently retains certain rights. As such, the U.S. Government has been
  granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
  worldwide license in the Software to reproduce, distribute copies to the public, prepare
  derivative works, and perform publicly and display publicly, and to permit other to do so.

  Developers: Yang Liu
             (Lawrence Berkeley National Lab, Computational Research Division).
*/
/*! @file
 * @brief One butterfly block through the C interface (c_bf_*), checked and timed on the CPU or,
 * with --hodlr_use_gpu 1, on the GPU:
 *
 *   1. BF1 from the entries of a kernel     (c_bf_construct_element_compute)
 *   2. BF2 from the entries of BF1          (its block callback is c_bf_extractelement on BF1)
 *   3. BF3 from the matvec of BF1           (its matvec callback is c_bf_mult on BF1; on the GPU a
 *                                            device black box, c_bf_gpu_mult on BF1)
 *
 * and BF1's own matvec and entry extraction.  These are the four operations STRUMPACK uses for
 * the F12, F21 and Schur blocks of a front.  Every line of results starts with "BFTEST".
 * This file works on the double data type.
 */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "dBPACK_wrapper.h"
#ifdef BPACK_EXAMPLE_GPU
#include "gpu/gpu_evaluators.h"
#endif

extern "C" void Cblacs_exit(int);

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

struct DriverOptions {
  std::string kernel = "dft";  // dft: cos(2 pi s x xi); laplace: 1/(4 pi r) between two squares
  int rows = 4096;             // M
  int cols = 4096;             // N
  double scale = 0.25;         // dft: the oscillation s (1: the real part of the DFT)
  double jitter = 1.0;         // dft: the points' random offsets, in grid cells (0: the DFT's grid)
  double gap = 0.5;            // laplace: the distance between the two unit squares
  double tolerance = 1e-4;     // tol_comp (and, as the library sets them, tol_rand, tol_Rdetect)
  std::string checks = "123";
  int nvec = 4;                // columns of the accuracy checks
  int sample_rows = 200;       // rows of BF1's check against the kernel
  std::vector<int> mult_nrhs{1, 16, 64};
  int mult_reps = 5;
  int extract_blocks = 64;     // intersections of the extraction timing
  int extract_size = 64;       // their rows and columns
  int extract_reps = 3;
  int gpu_evaluator = 1;       // with the GPU on, BF1's entries from a GPU evaluator of the kernel
  int gpu_matvec = 1;          // with the GPU on, BF3's black box BF1's matvec on device vectors
  bool show_help = false;
};

std::vector<int> parse_list(const std::string& value) {
  std::vector<int> out;
  std::stringstream in(value);
  std::string item;
  while (std::getline(in, item, ',')) {
    if (!item.empty()) out.push_back(std::stoi(item));
  }
  return out;
}

// The driver's own options; everything else (--nmin_leaf, --lrlevel, --forwardN15flag,
// --hodlr_use_gpu, ...) is ButterflyPACK's and goes to c_bpack_set_option_from_command_line,
// which in turn skips these.
DriverOptions parse_driver_options(int argc, char** argv) {
  DriverOptions o;
  for (int i = 1; i < argc; ++i) {
    std::string argument = argv[i];
    if (argument == "--help" || argument == "-h") {
      o.show_help = true;
      continue;
    }
    if (argument.rfind("--", 0) != 0) throw std::invalid_argument("unexpected argument: " + argument);
    argument.erase(0, 2);
    std::string name = argument, value;
    const size_t equals = argument.find('=');
    if (equals != std::string::npos) {
      name = argument.substr(0, equals);
      value = argument.substr(equals + 1);
    } else if (i + 1 < argc) {
      value = argv[++i];
    }
    if (name == "kernel") {
      o.kernel = value;
    } else if (name == "rows") {
      o.rows = std::stoi(value);
    } else if (name == "cols") {
      o.cols = std::stoi(value);
    } else if (name == "scale") {
      o.scale = std::stod(value);
    } else if (name == "jitter") {
      o.jitter = std::stod(value);
    } else if (name == "gap") {
      o.gap = std::stod(value);
    } else if (name == "tol_comp") {
      o.tolerance = std::stod(value);
    } else if (name == "checks") {
      o.checks = value;
    } else if (name == "nvec") {
      o.nvec = std::stoi(value);
    } else if (name == "sample_rows") {
      o.sample_rows = std::stoi(value);
    } else if (name == "mult_nrhs") {
      o.mult_nrhs = parse_list(value);
    } else if (name == "mult_reps") {
      o.mult_reps = std::stoi(value);
    } else if (name == "extract_blocks") {
      o.extract_blocks = std::stoi(value);
    } else if (name == "extract_size") {
      o.extract_size = std::stoi(value);
    } else if (name == "extract_reps") {
      o.extract_reps = std::stoi(value);
    } else if (name == "gpu_evaluator") {
      o.gpu_evaluator = std::stoi(value);
    } else if (name == "gpu_matvec") {
      o.gpu_matvec = std::stoi(value);
    }
  }
  if (o.kernel != "dft" && o.kernel != "laplace") throw std::invalid_argument("--kernel is dft or laplace");
  if (o.rows < 2 || o.cols < 2) throw std::invalid_argument("--rows and --cols must be at least 2");
  return o;
}

void print_usage(const char* program) {
  std::cout
      << "usage: " << program << " [driver options] [ButterflyPACK options]\n"
      << "  --kernel dft|laplace    dft: K(i,j) = cos(2 pi s x_i xi_j), x_i = (i + u_i)/M,\n"
      << "                          xi_j = j + v_j - N/2, u, v uniform in [0, jitter)\n"
      << "                          laplace: K = 1/(4 pi r), rows on the unit square z=0,\n"
      << "                          columns on the unit square z=gap\n"
      << "  --rows M --cols N       the block's size (4096 x 4096)\n"
      << "  --scale s               dft's oscillation (0.25; the ranks grow as s * Nmin_leaf)\n"
      << "  --jitter j              dft's random offsets (1; 0 is the DFT's grid, whose rows the\n"
      << "                          evenly spaced row sampling aliases)\n"
      << "  --gap g                 laplace's distance between the squares (0.5)\n"
      << "  --tol_comp t            compression tolerance (1e-4)\n"
      << "  --checks 123            which of BF1 (always built), BF2 from BF1's entries, BF3 from\n"
      << "                          BF1's matvec to build; 1 alone times only BF1\n"
      << "  --nvec k                columns of the accuracy checks (4)\n"
      << "  --sample_rows k         rows of BF1's check against the kernel (200)\n"
      << "  --mult_nrhs 1,16,64     right-hand sides of BF1's matvec timing (empty: none)\n"
      << "  --mult_reps r           repetitions of each (5)\n"
      << "  --extract_blocks b      intersections of BF1's extraction timing (64; 0: none)\n"
      << "  --extract_size s        rows and columns of each (64)\n"
      << "  --extract_reps r        repetitions (3)\n"
      << "  --gpu_evaluator 0|1     with the GPU on, BF1's entries from a GPU evaluator of the\n"
      << "                          kernel (1) or from the CPU callbacks (0)\n"
      << "  --gpu_matvec 0|1        with the GPU on, BF3's black box on the device (1: BF1's\n"
      << "                          matvec on device vectors, c_bf_set_gpu_matvec) or the CPU\n"
      << "                          callback (0)\n"
      << "The butterfly needs --lrlevel >= 1 (the driver sets 100); --hodlr_use_gpu 1 asks for the\n"
      << "GPU paths.  Results: lines starting with BFTEST.\n";
}

// The kernel, by the original 0-based indices of its rows and columns.
class Kernel {
 public:
  explicit Kernel(const DriverOptions& o) : laplace_(o.kernel == "laplace"), M_(o.rows), N_(o.cols),
                                            scale_(o.scale), gap_(o.gap) {
    side_r_ = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(M_))));
    side_c_ = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(N_))));
    if (!laplace_) {
      // x_i = (i + u_i)/M and xi_j = j + v_j - N/2, u and v uniform in [0, jitter): still sorted, so
      // the natural order is the tree's.  Without the jitter the matrix is the real part of a
      // DFT, whose rows the construction's evenly spaced row sampling aliases.
      std::mt19937_64 gen(2026);
      std::uniform_real_distribution<double> u(0.0, o.jitter);
      x_.resize(M_);
      xi_.resize(N_);
      for (int i = 0; i < M_; ++i) x_[i] = (i + u(gen)) / M_;
      for (int j = 0; j < N_; ++j) xi_[j] = j + u(gen) - 0.5 * N_;
    }
  }

  bool geometric() const { return laplace_; }
  double x(int64_t i) const { return x_[i]; }
  double xi(int64_t j) const { return xi_[j]; }
  int dim() const { return laplace_ ? 3 : 1; }

  double operator()(int64_t i, int64_t j) const {
    if (!laplace_) return std::cos(2.0 * kPi * scale_ * x_[i] * xi_[j]);
    double p[3], q[3];
    row_point(i, p);
    col_point(j, q);
    const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
    return 1.0 / (4.0 * kPi * std::sqrt(dx * dx + dy * dy + dz * dz));
  }

  // The points the row and column trees are built on: the squares' grid points for laplace;
  // for dft (built without geometry) x_i and xi_j, which the library does not read.
  std::vector<double> row_points() const {
    std::vector<double> xyz(static_cast<size_t>(M_) * dim());
    for (int i = 0; i < M_; ++i) {
      if (laplace_) {
        row_point(i, &xyz[3 * static_cast<size_t>(i)]);
      } else {
        xyz[i] = x_[i];
      }
    }
    return xyz;
  }
  std::vector<double> col_points() const {
    std::vector<double> xyz(static_cast<size_t>(N_) * dim());
    for (int j = 0; j < N_; ++j) {
      if (laplace_) {
        col_point(j, &xyz[3 * static_cast<size_t>(j)]);
      } else {
        xyz[j] = xi_[j];
      }
    }
    return xyz;
  }

  void row_point(int64_t i, double* p) const {
    const double h = 1.0 / side_r_;
    p[0] = (static_cast<double>(i % side_r_) + 0.5) * h;
    p[1] = (static_cast<double>(i / side_r_) + 0.5) * h;
    p[2] = 0.0;
  }
  void col_point(int64_t j, double* q) const {
    const double h = 1.0 / side_c_;
    q[0] = (static_cast<double>(j % side_c_) + 0.5) * h;
    q[1] = (static_cast<double>(j / side_c_) + 0.5) * h;
    q[2] = gap_;
  }

 private:

  bool laplace_;
  int M_, N_;
  double scale_, gap_;
  int side_r_ = 1, side_c_ = 1;
  std::vector<double> x_, xi_;  // dft's points
};

// ---- BF1's entries: the kernel, by the indices of the row and column trees ----

struct KernelSource {
  const Kernel* kernel = nullptr;
  const std::vector<int>* row_perm = nullptr;  // tree order (1-based) -> original (1-based)
  const std::vector<int>* col_perm = nullptr;
  int myid = 0;
};

// A positive index is a row and a negative one a column, in either order.
void kernel_entry(int* m, int* n, double* value, C2Fptr quant) {
  const KernelSource* s = static_cast<const KernelSource*>(quant);
  int r = *m, c = *n;
  if (r < 0) std::swap(r, c);
  *value = (*s->kernel)((*s->row_perm)[r - 1] - 1, (*s->col_perm)[-c - 1] - 1);
}

// The intersections are listed on every rank; a rank fills those of its own (1 x 1 process
// grids, which is what a butterfly's construction asks for), stacked in order, column major.
void kernel_block(int* Ninter, int* /*Nallrows*/, int* /*Nallcols*/, int64_t* /*Nalldat_loc*/, int* allrows,
                  int* allcols, double* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap, int* pmaps,
                  C2Fptr quant) {
  const KernelSource* s = static_cast<const KernelSource*>(quant);
  int64_t at = 0, r0 = 0, c0 = 0;
  for (int k = 0; k < *Ninter; ++k) {
    const int m = rowidx[k], n = colidx[k];
    const int pg = pgidx[k];
    const int nprow = pmaps[pg], npcol = pmaps[pg + *Npmap], head = pmaps[pg + 2 * *Npmap];
    if (nprow * npcol != 1) throw std::runtime_error("kernel_block: only 1 x 1 process grids are handled");
    if (head == s->myid) {
      double* out = alldat_loc + at;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int j = 0; j < n; ++j) {
        const int64_t jo = (*s->col_perm)[allcols[c0 + j] - 1] - 1;
        for (int i = 0; i < m; ++i) {
          out[i + static_cast<int64_t>(j) * m] = (*s->kernel)((*s->row_perm)[allrows[r0 + i] - 1] - 1, jo);
        }
      }
      at += static_cast<int64_t>(m) * n;
    }
    r0 += m;
    c0 += n;
  }
}

// ---- BF2's entries: extracted from BF1 ----

struct Handles {
  F2Cptr bf = nullptr, option = nullptr, stats = nullptr, msh = nullptr, ker = nullptr;
};

struct ExtractionSource {
  Handles* bf1 = nullptr;
  F2Cptr* ptree = nullptr;
  double seconds = 0.0;
  int64_t calls = 0;
};

void extraction_entry(int*, int*, double*, C2Fptr) {
  throw std::runtime_error("BF2 is built with elem_extract 1: its entry callback is never called");
}

void extraction_block(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows,
                      int* allcols, double* alldat_loc, int* rowidx, int* colidx, int* pgidx, int* Npmap,
                      int* pmaps, C2Fptr quant) {
  ExtractionSource* s = static_cast<ExtractionSource*>(quant);
  const double t0 = MPI_Wtime();
  d_c_bf_extractelement(&s->bf1->bf, &s->bf1->option, &s->bf1->msh, &s->bf1->stats, s->ptree, Ninter, Nallrows,
                        Nallcols, Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps);
  s->seconds += MPI_Wtime() - t0;
  ++s->calls;
}

// ---- BF3's matvec: BF1's ----

struct MatvecSource {
  Handles* bf1 = nullptr;
  F2Cptr* ptree = nullptr;
  std::vector<double> buffer;
  double seconds = 0.0;
  int64_t calls = 0, columns = 0;
};

// xout = b xout + a op(BF1) xin, with nin and nout the local lengths of xin and xout
void matvec_bf1(char const* trans, int* nin, int* nout, int* nvec, double const* xin, double* xout, C2Fptr quant,
                double* a, double* b) {
  MatvecSource* s = static_cast<MatvecSource*>(quant);
  const double t0 = MPI_Wtime();
  const size_t count = static_cast<size_t>(*nout) * static_cast<size_t>(*nvec);
  s->buffer.resize(std::max<size_t>(count, 1));
  d_c_bf_mult(trans, xin, s->buffer.data(), nin, nout, nvec, &s->bf1->bf, &s->bf1->option, &s->bf1->stats,
              s->ptree);
  for (size_t k = 0; k < count; ++k) xout[k] = *b * xout[k] + *a * s->buffer[k];
  s->seconds += MPI_Wtime() - t0;
  ++s->calls;
  s->columns += *nvec;
}

#ifdef BPACK_EXAMPLE_GPU
// The same on device vectors (bpack_gpu_bf_matvec): y = op(BF1) x, complete on return
void matvec_bf1_device(char trans, int nin, int nout, int nvec, const void* x, void* y, void*, void* user) {
  MatvecSource* s = static_cast<MatvecSource*>(user);
  const double t0 = MPI_Wtime();
  const double one = 1.0, zero = 0.0;
  double flops = 0.0;
  int info = 0;
  d_c_bf_gpu_mult(s->bf1->bf, *s->ptree, &trans, &nvec, x, &nin, y, &nout, &one, &zero, &flops, &info);
  if (info != 0) throw std::runtime_error("BF1's device matvec declined");
  s->seconds += MPI_Wtime() - t0;
  ++s->calls;
  s->columns += nvec;
}
#endif

void no_distance(int*, int*, double*, C2Fptr) {}
void no_admissibility(int*, int*, int* value, C2Fptr) { *value = 1; }

double stat(F2Cptr& stats, const char* name) {
  double value = 0.0;
  d_c_bpack_getstats(&stats, name, &value);
  return value;
}

// Timers a call accumulates into a block's stats: read before, then the per-call difference
// (the maximum over the ranks) after, printed as name=seconds.
class StatDelta {
 public:
  StatDelta(F2Cptr& stats, std::vector<const char*> names) : stats_(stats), names_(std::move(names)) {
    for (const char* n : names_) before_.push_back(stat(stats_, n));
  }
  std::string report(int calls) {
    std::vector<double> d(names_.size());
    for (size_t k = 0; k < names_.size(); ++k) d[k] = (stat(stats_, names_[k]) - before_[k]) / std::max(calls, 1);
    MPI_Allreduce(MPI_IN_PLACE, d.data(), static_cast<int>(d.size()), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    std::ostringstream out;
    for (size_t k = 0; k < names_.size(); ++k) out << " " << names_[k] << "=" << d[k];
    return out.str();
  }

 private:
  F2Cptr& stats_;
  std::vector<const char*> names_;
  std::vector<double> before_;
};

// ---- helpers over the distributed vectors (each rank its rows / columns in tree order) ----

std::vector<double> random_vector(size_t count, unsigned seed) {
  std::mt19937_64 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::vector<double> v(std::max<size_t>(count, 1));
  for (auto& x : v) x = dist(gen);
  return v;
}

// ||a - b||_F / ||a||_F over all ranks
double relative_difference(const std::vector<double>& a, const std::vector<double>& b, size_t count) {
  double sums[2] = {0.0, 0.0};
  for (size_t k = 0; k < count; ++k) {
    sums[0] += (a[k] - b[k]) * (a[k] - b[k]);
    sums[1] += a[k] * a[k];
  }
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sums[1] > 0.0 ? std::sqrt(sums[0] / sums[1]) : std::sqrt(sums[0]);
}

void apply(Handles& h, F2Cptr& ptree, char trans, const std::vector<double>& x, std::vector<double>& y, int nin,
           int nout, int nvec) {
  y.assign(std::max<size_t>(static_cast<size_t>(nout) * nvec, 1), 0.0);
  d_c_bf_mult(&trans, x.data(), y.data(), &nin, &nout, &nvec, &h.bf, &h.option, &h.stats, &ptree);
}

void print_layout(int myid, const char* name, Handles& h) {
  int L = 0, half = 0, M = 0, N = 0, Ml = 0, Nl = 0, nlev = 0, nblk = 0;
  d_c_bf_get_layout(h.bf, &L, &half, &M, &N, &Ml, &Nl, &nlev, &nblk);
  int total = nblk;
  MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid == 0) {
    std::cout << "BFTEST layout " << name << " M=" << M << " N=" << N << " level_butterfly=" << L
              << " level_half=" << half << " blocks=" << total << " rank_max=" << stat(h.stats, "Rank_max")
              << std::endl;
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int myid = 0, nproc = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  int return_code = 0;

  try {
    const DriverOptions opt = parse_driver_options(argc, argv);
    if (opt.show_help) {
      if (myid == 0) print_usage(argv[0]);
      MPI_Finalize();
      return 0;
    }
    const Kernel kernel(opt);
    int M = opt.rows, N = opt.cols;

    std::vector<int> groups(nproc);
    std::iota(groups.begin(), groups.end(), 0);
    MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
    F2Cptr ptree = nullptr, option = nullptr;
    d_c_bpack_createptree(&nproc, groups.data(), &fcomm, &ptree);
    d_c_bpack_createoption(&option);

    // the driver's defaults, then the command line over them
    d_c_bpack_set_D_option(&option, "tol_comp", opt.tolerance);
    d_c_bpack_set_D_option(&option, "tol_rand", opt.tolerance);
    d_c_bpack_set_D_option(&option, "tol_Rdetect", 0.1 * opt.tolerance);
    d_c_bpack_set_I_option(&option, "nogeo", kernel.geometric() ? 0 : 1);
    d_c_bpack_set_I_option(&option, "xyzsort", kernel.geometric() ? 1 : 0);
    d_c_bpack_set_I_option(&option, "Nmin_leaf", 64);
    d_c_bpack_set_I_option(&option, "LRlevel", 100);  // a butterfly, not a low-rank block
    d_c_bpack_set_I_option(&option, "cpp", 1);
    d_c_bpack_set_I_option(&option, "verbosity", -1);
    d_c_bpack_set_I_option(&option, "less_adapt", 1);  // as STRUMPACK builds its Schur blocks
    d_c_bpack_set_I_option(&option, "elem_extract", 0);
    d_c_bpack_set_option_from_command_line(argc, argv, option);

    double gpu = 0.0, n15 = 0.0, lrlevel = 0.0, leaf = 0.0, knn = 0.0, elem = 0.0, less_adapt = 0.0;
    d_c_bpack_getoption(&option, "HODLR_use_gpu", &gpu);
    d_c_bpack_getoption(&option, "forwardN15flag", &n15);
    d_c_bpack_getoption(&option, "LRlevel", &lrlevel);
    d_c_bpack_getoption(&option, "Nmin_leaf", &leaf);
    d_c_bpack_getoption(&option, "knn", &knn);
    d_c_bpack_getoption(&option, "elem_extract", &elem);
    d_c_bpack_getoption(&option, "less_adapt", &less_adapt);
    int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    nthreads = omp_get_num_threads();
#endif
    if (myid == 0) {
      std::cout << "BFTEST config kernel=" << opt.kernel << " M=" << M << " N=" << N << " scale=" << opt.scale
                << " jitter=" << opt.jitter << " gap=" << opt.gap << " tol=" << opt.tolerance << " nproc=" << nproc
                << " threads=" << nthreads << " hodlr_use_gpu=" << gpu << " forwardN15flag=" << n15
                << " lrlevel=" << lrlevel << " Nmin_leaf=" << leaf << " knn=" << knn
                << " elem_extract=" << elem << " less_adapt=" << less_adapt << std::endl;
    }

    // ---- the row and column trees, each from a mesh of its own points ----
    // c_bpack_construct_init lays out a HODLR on each, which with the GPU on requires LRlevel 0;
    // the meshes are all that is wanted of it, so it gets the options without either
    const int dim = kernel.dim();
    std::vector<double> row_xyz = kernel.row_points(), col_xyz = kernel.col_points();
    std::vector<int> row_perm(M), col_perm(N);
    std::vector<int> nns(1, 0);  // no neighbour lists (knn from the geometry, if asked for)
    int Ml = 0, Nl = 0, nlevel = 0, Mtree = M, Ntree = N;
    int dim_r = dim, dim_c = dim;
    F2Cptr mesh_option = nullptr;
    d_c_bpack_copyoption(&option, &mesh_option);
    d_c_bpack_set_I_option(&mesh_option, "HODLR_use_gpu", 0);
    d_c_bpack_set_I_option(&mesh_option, "LRlevel", 0);
    F2Cptr bmat_r = nullptr, stats_r = nullptr, msh_r = nullptr, ker_r = nullptr;
    F2Cptr bmat_c = nullptr, stats_c = nullptr, msh_c = nullptr, ker_c = nullptr;
    d_c_bpack_createstats(&stats_r);
    d_c_bpack_createstats(&stats_c);
    d_c_bpack_construct_init(&M, &dim_r, row_xyz.data(), nns.data(), &nlevel, &Mtree, row_perm.data(), &Ml,
                             &bmat_r, &mesh_option, &stats_r, &msh_r, &ker_r, &ptree, &no_distance,
                             &no_admissibility, nullptr);
    d_c_bpack_construct_init(&N, &dim_c, col_xyz.data(), nns.data(), &nlevel, &Ntree, col_perm.data(), &Nl,
                             &bmat_c, &mesh_option, &stats_c, &msh_c, &ker_c, &ptree, &no_distance,
                             &no_admissibility, nullptr);

    KernelSource kernel_source;
    kernel_source.kernel = &kernel;
    kernel_source.row_perm = &row_perm;
    kernel_source.col_perm = &col_perm;
    kernel_source.myid = myid;

    // ---- 1. BF1 from the kernel's entries ----
    Handles bf1;
    d_c_bpack_copyoption(&option, &bf1.option);
    d_c_bpack_createstats(&bf1.stats);
    int Ml1 = Ml, Nl1 = Nl;
    d_c_bf_construct_init(&M, &N, &Ml1, &Nl1, nns.data(), nns.data(), &msh_r, &msh_c, &bf1.bf, &bf1.option,
                          &bf1.stats, &bf1.msh, &bf1.ker, &ptree, &no_distance, &no_admissibility, &kernel_source);
#ifdef BPACK_EXAMPLE_GPU
    // the kernel's GPU evaluator: its tables by the ids it receives, the trees' order (the CPU
    // callbacks' indices minus one), so through the trees' permutations
    if (gpu > 0 && opt.gpu_evaluator) {
      if (kernel.geometric()) {
        std::vector<double> p(static_cast<size_t>(M) * 3), q(static_cast<size_t>(N) * 3);
        for (int r = 0; r < M; ++r) kernel.row_point(row_perm[r] - 1, &p[3 * static_cast<size_t>(r)]);
        for (int c = 0; c < N; ++c) kernel.col_point(col_perm[c] - 1, &q[3 * static_cast<size_t>(c)]);
        butterfly_gpu_register_laplace(&bf1.bf, p.data(), M, q.data(), N);
      } else {
        std::vector<double> xr(M), xc(N);
        for (int r = 0; r < M; ++r) xr[r] = kernel.x(row_perm[r] - 1);
        for (int c = 0; c < N; ++c) xc[c] = kernel.xi(col_perm[c] - 1);
        butterfly_gpu_register_dft(&bf1.bf, xr.data(), M, xc.data(), N, opt.scale);
      }
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    d_c_bf_construct_element_compute(&bf1.bf, &bf1.option, &bf1.stats, &bf1.msh, &bf1.ker, &ptree, &kernel_entry,
                                     &kernel_block, &kernel_source);
    MPI_Barrier(MPI_COMM_WORLD);
    const double t_bf1 = MPI_Wtime() - t0;
    {
      double entry = stat(bf1.stats, "Time_Entry"), flops = stat(bf1.stats, "Flop_Fill");
      double mem = stat(bf1.stats, "Mem_Fill"), warm = stat(bf1.stats, "Time_GPU_Warmup");
      double entry_max = entry;
      MPI_Allreduce(MPI_IN_PLACE, &entry_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &flops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (myid == 0) {
        std::cout << "BFTEST construct bf=1 from=kernel seconds=" << t_bf1 - warm << " entry_max=" << entry_max
                  << " flops=" << flops << " mem_MB=" << mem << " gpu_warmup=" << warm << std::endl;
      }
    }
    print_layout(myid, "bf1", bf1);
    if (Ml1 != Ml || Nl1 != Nl) throw std::runtime_error("BF1's local sizes differ from its trees'");

    // BF1 against the kernel, on sampled rows: y = BF1 x, and the same rows of K x summed directly
    {
      std::vector<double> x = random_vector(static_cast<size_t>(Nl) * opt.nvec, 1234u + myid), y;
      apply(bf1, ptree, 'N', x, y, Nl, Ml, opt.nvec);
      const int S = std::min(opt.sample_rows, M);
      std::vector<int> sample(S);
      if (myid == 0) {
        std::mt19937 gen(99);
        std::vector<int> all(M);
        std::iota(all.begin(), all.end(), 0);
        std::shuffle(all.begin(), all.end(), gen);
        std::copy(all.begin(), all.begin() + S, sample.begin());
      }
      MPI_Bcast(sample.data(), S, MPI_INT, 0, MPI_COMM_WORLD);
      std::unordered_map<int, int> local_row;  // original row -> local row
      for (int l = 1; l <= Ml; ++l) {
        int l1 = l, old = 0;
        d_c_bf_new2old_row(&msh_r, &l1, &old);
        local_row[old - 1] = l - 1;
      }
      std::vector<int> col_old(Nl);
      for (int l = 1; l <= Nl; ++l) {
        int l1 = l, old = 0;
        d_c_bf_new2old_col(&msh_c, &l1, &old);
        col_old[l - 1] = old - 1;
      }
      std::vector<double> exact(static_cast<size_t>(S) * opt.nvec, 0.0), approx(exact.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int s = 0; s < S; ++s) {
        for (int l = 0; l < Nl; ++l) {
          const double k = kernel(sample[s], col_old[l]);
          for (int v = 0; v < opt.nvec; ++v) exact[s + static_cast<size_t>(v) * S] += k * x[l + static_cast<size_t>(v) * Nl];
        }
        auto it = local_row.find(sample[s]);
        if (it != local_row.end()) {
          for (int v = 0; v < opt.nvec; ++v) approx[s + static_cast<size_t>(v) * S] = y[it->second + static_cast<size_t>(v) * Ml];
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, exact.data(), static_cast<int>(exact.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, approx.data(), static_cast<int>(approx.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      double num = 0.0, den = 0.0;
      for (size_t k = 0; k < exact.size(); ++k) {
        num += (exact[k] - approx[k]) * (exact[k] - approx[k]);
        den += exact[k] * exact[k];
      }
      if (myid == 0) {
        std::cout << "BFTEST check bf=1 against=kernel rows=" << S << " nvec=" << opt.nvec
                  << " rel=" << std::sqrt(num / den) << std::endl;
      }
    }

    // ---- BF1's matvec: on the GPU, first against the host's.  Every device kernel sums in k
    // order, as OpenBLAS does at an even column count: there the two must agree to the bit.  At an
    // odd count OpenBLAS changes order on blocks of real rank, so there only to rounding.
    if (gpu > 0) {
      std::vector<double> xn = random_vector(static_cast<size_t>(Nl) * opt.nvec, 555u + myid);
      std::vector<double> xt = random_vector(static_cast<size_t>(Ml) * opt.nvec, 666u + myid);
      std::vector<double> yn_gpu, yt_gpu, yn_cpu, yt_cpu;
      apply(bf1, ptree, 'N', xn, yn_gpu, Nl, Ml, opt.nvec);
      apply(bf1, ptree, 'T', xt, yt_gpu, Ml, Nl, opt.nvec);
      d_c_bpack_set_I_option(&bf1.option, "HODLR_use_gpu", 0);
      apply(bf1, ptree, 'N', xn, yn_cpu, Nl, Ml, opt.nvec);
      apply(bf1, ptree, 'T', xt, yt_cpu, Ml, Nl, opt.nvec);
      d_c_bpack_set_I_option(&bf1.option, "HODLR_use_gpu", static_cast<int>(gpu));
      const double rn = relative_difference(yn_cpu, yn_gpu, static_cast<size_t>(Ml) * opt.nvec);
      const double rt = relative_difference(yt_cpu, yt_gpu, static_cast<size_t>(Nl) * opt.nvec);
      const bool exact = opt.nvec % 2 == 0;
      const bool pass = exact ? (rn == 0.0 && rt == 0.0) : (rn < 1e-13 && rt < 1e-13);
      if (myid == 0) {
        std::cout << "BFTEST check bf=1 mult gpu_vs_cpu rel_N=" << rn << " rel_T=" << rt << " expect="
                  << (exact ? "exact" : "1e-13") << " " << (pass ? "PASS" : "FAIL") << std::endl;
      }
    }
    for (int nrhs : opt.mult_nrhs) {
      for (char trans : {'N', 'T'}) {
        const int nin = trans == 'N' ? Nl : Ml, nout = trans == 'N' ? Ml : Nl;
        std::vector<double> x = random_vector(static_cast<size_t>(nin) * nrhs, 77u + myid), y;
        apply(bf1, ptree, trans, x, y, nin, nout, nrhs);  // warm-up
        double best = 1e300, sum = 0.0;
        StatDelta parts(bf1.stats, {"Time_BF_MVP_Gemm", "Time_BF_MVP_Exchange", "Time_BF_MVP_All2All",
                                    "Time_BF_MVP_Other"});
        for (int r = 0; r < opt.mult_reps; ++r) {
          MPI_Barrier(MPI_COMM_WORLD);
          const double s0 = MPI_Wtime();
          apply(bf1, ptree, trans, x, y, nin, nout, nrhs);
          MPI_Barrier(MPI_COMM_WORLD);
          const double dt = MPI_Wtime() - s0;
          best = std::min(best, dt);
          sum += dt;
        }
        double flops = stat(bf1.stats, "Flop_C_Mult");
        MPI_Allreduce(MPI_IN_PLACE, &flops, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const std::string split = parts.report(opt.mult_reps);
        if (myid == 0) {
          std::cout << "BFTEST mult bf=1 op=" << trans << " nrhs=" << nrhs << " seconds_min=" << best
                    << " seconds_avg=" << sum / std::max(opt.mult_reps, 1) << " gflops="
                    << flops / best * 1e-9 << split << std::endl;
        }
        // the same product on vectors already on the device (c_bf_gpu_mult), the form a caller
        // chaining operators uses: what is left is the sweep alone
        if (gpu > 0) {
          void *dx = nullptr, *dy = nullptr;
          int64_t ein = static_cast<int64_t>(std::max(nin, 1)) * nrhs, eout = static_cast<int64_t>(std::max(nout, 1)) * nrhs;
          d_c_bpack_gpu_workspace_alloc(&ein, &dx);
          d_c_bpack_gpu_workspace_alloc(&eout, &dy);
          int to_device = 1;
          d_c_bpack_gpu_workspace_copy(dx, x.data(), &ein, &to_device);
          int ldx = std::max(nin, 1), ldy = std::max(nout, 1), info = 0;
          double one = 1.0, zero = 0.0, fl = 0.0;
          char t = trans;
          d_c_bf_gpu_mult(bf1.bf, ptree, &t, &nrhs, dx, &ldx, dy, &ldy, &one, &zero, &fl, &info);  // warm-up
          double dbest = 1e300;
          for (int r = 0; r < opt.mult_reps; ++r) {
            MPI_Barrier(MPI_COMM_WORLD);
            const double s0 = MPI_Wtime();
            d_c_bf_gpu_mult(bf1.bf, ptree, &t, &nrhs, dx, &ldx, dy, &ldy, &one, &zero, &fl, &info);
            MPI_Barrier(MPI_COMM_WORLD);
            dbest = std::min(dbest, MPI_Wtime() - s0);
          }
          MPI_Allreduce(MPI_IN_PLACE, &fl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          if (myid == 0) {
            std::cout << "BFTEST mult_device bf=1 op=" << trans << " nrhs=" << nrhs << " seconds_min=" << dbest
                      << " gflops=" << fl / dbest * 1e-9 << " info=" << info << std::endl;
          }
          d_c_bpack_gpu_workspace_free(dy);
          d_c_bpack_gpu_workspace_free(dx);
        }
      }
    }

    // ---- BF1's entry extraction: random intersections, each owned by one rank, against the kernel ----
    if (opt.extract_blocks > 0) {
      const int B = opt.extract_blocks, sz = std::min({opt.extract_size, M, N});
      std::vector<int> rows(static_cast<size_t>(B) * sz), cols(static_cast<size_t>(B) * sz);
      if (myid == 0) {
        std::mt19937 gen(7);
        std::uniform_int_distribution<int> pick_r(1, M), pick_c(1, N);
        for (auto& r : rows) r = pick_r(gen);
        for (auto& c : cols) c = pick_c(gen);
      }
      MPI_Bcast(rows.data(), static_cast<int>(rows.size()), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(cols.data(), static_cast<int>(cols.size()), MPI_INT, 0, MPI_COMM_WORLD);
      std::vector<int> rowidx(B, sz), colidx(B, sz), pgidx(B);
      std::vector<int> pmaps(static_cast<size_t>(3) * nproc);
      for (int p = 0; p < nproc; ++p) {
        pmaps[p] = 1;
        pmaps[p + nproc] = 1;
        pmaps[p + 2 * nproc] = p;
      }
      int64_t mine = 0;
      for (int b = 0; b < B; ++b) {
        pgidx[b] = b % nproc;
        if (pgidx[b] == myid) mine += static_cast<int64_t>(sz) * sz;
      }
      std::vector<double> data(std::max<int64_t>(mine, 1));
      int Ninter = B, Nallrows = B * sz, Nallcols = B * sz, Npmap = nproc;
      double best = 1e300, sum = 0.0;
      StatDelta parts(bf1.stats, {"Time_Entry_Traverse", "Time_Entry_BF", "Time_Extract_Gemm", "Time_Entry_Comm"});
      for (int r = 0; r < std::max(opt.extract_reps, 1); ++r) {
        MPI_Barrier(MPI_COMM_WORLD);
        const double s0 = MPI_Wtime();
        d_c_bf_extractelement(&bf1.bf, &bf1.option, &bf1.msh, &bf1.stats, &ptree, &Ninter, &Nallrows, &Nallcols,
                              &mine, rows.data(), cols.data(), data.data(), rowidx.data(), colidx.data(),
                              pgidx.data(), &Npmap, pmaps.data());
        MPI_Barrier(MPI_COMM_WORLD);
        const double dt = MPI_Wtime() - s0;
        best = std::min(best, dt);
        sum += dt;
      }
      double sums[2] = {0.0, 0.0};
      int64_t at = 0;
      for (int b = 0; b < B; ++b) {
        if (pgidx[b] != myid) continue;
        for (int j = 0; j < sz; ++j) {
          for (int i = 0; i < sz; ++i) {
            const double k = kernel(row_perm[rows[static_cast<size_t>(b) * sz + i] - 1] - 1,
                                    col_perm[cols[static_cast<size_t>(b) * sz + j] - 1] - 1);
            const double d = data[at + i + static_cast<int64_t>(j) * sz] - k;
            sums[0] += d * d;
            sums[1] += k * k;
          }
        }
        at += static_cast<int64_t>(sz) * sz;
      }
      MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      const std::string split = parts.report(std::max(opt.extract_reps, 1));
      if (myid == 0) {
        std::cout << "BFTEST extract bf=1 blocks=" << B << " size=" << sz << " seconds_min=" << best
                  << " seconds_avg=" << sum / std::max(opt.extract_reps, 1)
                  << " rel=" << std::sqrt(sums[0] / sums[1]) << split << std::endl;
      }
      // on the GPU, against the host's extraction of the same butterfly: the products are the
      // same up to their order, so they agree to rounding
      if (gpu > 0) {
        std::vector<double> host(data.size());
        d_c_bpack_set_I_option(&bf1.option, "HODLR_use_gpu", 0);
        MPI_Barrier(MPI_COMM_WORLD);
        const double s0 = MPI_Wtime();
        d_c_bf_extractelement(&bf1.bf, &bf1.option, &bf1.msh, &bf1.stats, &ptree, &Ninter, &Nallrows, &Nallcols,
                              &mine, rows.data(), cols.data(), host.data(), rowidx.data(), colidx.data(),
                              pgidx.data(), &Npmap, pmaps.data());
        MPI_Barrier(MPI_COMM_WORLD);
        const double t_host = MPI_Wtime() - s0;
        d_c_bpack_set_I_option(&bf1.option, "HODLR_use_gpu", static_cast<int>(gpu));
        const double rel = relative_difference(host, data, static_cast<size_t>(mine));
        if (myid == 0) {
          std::cout << "BFTEST check bf=1 extract gpu_vs_cpu rel=" << rel << " seconds_cpu=" << t_host << std::endl;
        }
      }
    }

    // the vectors BF2 and BF3 are compared with BF1 on
    std::vector<double> xcheck = random_vector(static_cast<size_t>(Nl) * opt.nvec, 4321u + myid), y1;
    apply(bf1, ptree, 'N', xcheck, y1, Nl, Ml, opt.nvec);
    std::vector<double> xcheck_t = random_vector(static_cast<size_t>(Ml) * opt.nvec, 8765u + myid), y1t;
    apply(bf1, ptree, 'T', xcheck_t, y1t, Ml, Nl, opt.nvec);

    // ---- 2. BF2 from BF1's entries ----
    Handles bf2;
    if (opt.checks.find('2') != std::string::npos) {
      ExtractionSource source;
      source.bf1 = &bf1;
      source.ptree = &ptree;
      d_c_bpack_copyoption(&option, &bf2.option);
      d_c_bpack_set_I_option(&bf2.option, "elem_extract", 1);  // BF1's extraction is collective
      d_c_bpack_createstats(&bf2.stats);
      int Ml2 = Ml, Nl2 = Nl;
      d_c_bf_construct_init(&M, &N, &Ml2, &Nl2, nns.data(), nns.data(), &msh_r, &msh_c, &bf2.bf, &bf2.option,
                            &bf2.stats, &bf2.msh, &bf2.ker, &ptree, &no_distance, &no_admissibility, &source);
      MPI_Barrier(MPI_COMM_WORLD);
      // where BF1's extraction spends the time BF2 waits on it
      StatDelta served(bf1.stats, {"Time_Entry_Traverse", "Time_Entry_BF", "Time_Extract_Gemm", "Time_Entry_Comm"});
      t0 = MPI_Wtime();
      d_c_bf_construct_element_compute(&bf2.bf, &bf2.option, &bf2.stats, &bf2.msh, &bf2.ker, &ptree,
                                       &extraction_entry, &extraction_block, &source);
      MPI_Barrier(MPI_COMM_WORLD);
      const double t_bf2 = MPI_Wtime() - t0;
      const std::string split = served.report(1);
      double extract_max = source.seconds, mem = stat(bf2.stats, "Mem_Fill");
      MPI_Allreduce(MPI_IN_PLACE, &extract_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std::vector<double> y2, y2t;
      apply(bf2, ptree, 'N', xcheck, y2, Nl, Ml, opt.nvec);
      apply(bf2, ptree, 'T', xcheck_t, y2t, Ml, Nl, opt.nvec);
      const double rel = relative_difference(y1, y2, static_cast<size_t>(Ml) * opt.nvec);
      const double relt = relative_difference(y1t, y2t, static_cast<size_t>(Nl) * opt.nvec);
      if (myid == 0) {
        std::cout << "BFTEST construct bf=2 from=bf1_entries seconds=" << t_bf2 << " extract_max=" << extract_max
                  << " callbacks=" << source.calls << " mem_MB=" << mem << split << std::endl;
      }
      print_layout(myid, "bf2", bf2);
      if (myid == 0) std::cout << "BFTEST check bf=2 against=bf1 rel_N=" << rel << " rel_T=" << relt << std::endl;
    }

    // ---- 3. BF3 from BF1's matvec ----
    Handles bf3;
    if (opt.checks.find('3') != std::string::npos) {
      MatvecSource source;
      source.bf1 = &bf1;
      source.ptree = &ptree;
      d_c_bpack_copyoption(&option, &bf3.option);
      d_c_bpack_createstats(&bf3.stats);
      int Ml3 = Ml, Nl3 = Nl;
      d_c_bf_construct_init(&M, &N, &Ml3, &Nl3, nns.data(), nns.data(), &msh_r, &msh_c, &bf3.bf, &bf3.option,
                            &bf3.stats, &bf3.msh, &bf3.ker, &ptree, &no_distance, &no_admissibility, &source);
#ifdef BPACK_EXAMPLE_GPU
      if (gpu > 0 && opt.gpu_matvec) d_c_bf_set_gpu_matvec(&bf3.bf, &matvec_bf1_device, &source);
#endif
      MPI_Barrier(MPI_COMM_WORLD);
      StatDelta parts(bf3.stats, {"Time_Random_Sample", "Time_Random_Resolve", "Time_Random_OneBlock"});
      StatDelta served(bf1.stats, {"Time_BF_MVP_Gemm", "Time_BF_MVP_Exchange", "Time_BF_MVP_All2All"});
      t0 = MPI_Wtime();
      d_c_bf_construct_matvec_compute(&bf3.bf, &bf3.option, &bf3.stats, &bf3.msh, &bf3.ker, &ptree, &matvec_bf1,
                                      &source);
      MPI_Barrier(MPI_COMM_WORLD);
      const double t_bf3 = MPI_Wtime() - t0;
      const std::string split = parts.report(1) + served.report(1);
      double matvec_max = source.seconds, mem = stat(bf3.stats, "Mem_Fill");
      MPI_Allreduce(MPI_IN_PLACE, &matvec_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mem, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      std::vector<double> y3, y3t;
      apply(bf3, ptree, 'N', xcheck, y3, Nl, Ml, opt.nvec);
      apply(bf3, ptree, 'T', xcheck_t, y3t, Ml, Nl, opt.nvec);
      const double rel = relative_difference(y1, y3, static_cast<size_t>(Ml) * opt.nvec);
      const double relt = relative_difference(y1t, y3t, static_cast<size_t>(Nl) * opt.nvec);
      if (myid == 0) {
        std::cout << "BFTEST construct bf=3 from=bf1_matvec seconds=" << t_bf3 << " matvec_max=" << matvec_max
                  << " callbacks=" << source.calls << " columns=" << source.columns << " mem_MB=" << mem << split
                  << std::endl;
      }
      print_layout(myid, "bf3", bf3);
      if (myid == 0) std::cout << "BFTEST check bf=3 against=bf1 rel_N=" << rel << " rel_T=" << relt << std::endl;
    }

    for (Handles* h : {&bf3, &bf2, &bf1}) {
      if (h->bf == nullptr) continue;
      d_c_bf_deletebf(&h->bf);
      d_c_bpack_deletemesh(&h->msh);
      d_c_bpack_deletekernelquant(&h->ker);
      d_c_bpack_deletestats(&h->stats);
      d_c_bpack_deleteoption(&h->option);
    }
    d_c_bpack_deletestats(&stats_r);
    d_c_bpack_deletestats(&stats_c);
    d_c_bpack_deleteoption(&mesh_option);
    d_c_bpack_deleteoption(&option);
    d_c_bpack_deleteproctree(&ptree);
#ifdef BPACK_EXAMPLE_GPU
    butterfly_gpu_release();
#endif
  } catch (const std::exception& error) {
    std::cerr << "Butterfly_Driver error on rank " << myid << ": " << error.what() << std::endl;
    return_code = 1;
    MPI_Abort(MPI_COMM_WORLD, return_code);
  }

  Cblacs_exit(1);
  MPI_Finalize();
  return return_code;
}
