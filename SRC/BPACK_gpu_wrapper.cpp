// C entry points of the GPU backend of the formats other than H2 (HODLR for
// now; hodlr_gpu/), called from Fortran (BPACK_gpu.f90) and from
// C_BPACK_wrapper.cpp.  PrecisionPreprocessing.sh copies this file into each
// precision; the double and double complex libraries of a GPU build
// (enable_h2_gpu) get the backend, the others stubs that report its absence.

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bpack_env.hpp"

#if defined(H2_HAVE_GPU) && defined(DAT) && (DAT == 0 || DAT == 1)
#define BPACK_GPU_ENABLED 1
#include "hodlr_gpu/gpu_state.hpp"
#endif

namespace {

[[noreturn]] void gpu_fail(const char* where, const char* what) {
  std::fprintf(stderr, "ButterflyPACK GPU backend error in %s: %s\n", where, what);
  std::fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort();
}

#ifdef BPACK_GPU_ENABLED
#if DAT == 1
using Scalar = double;
#else
using Scalar = fmm::gpu::dcomplex;
#endif
using State = bpack::gpu::GpuState<Scalar>;

State* state_of(void* gpu, const char* where) {
  if (gpu == nullptr) gpu_fail(where, "the matrix has no GPU state");
  return static_cast<State*>(gpu);
}

bpack::gpu::HodlrConstruct<Scalar>& construct_of(void* gpu, const char* where) {
  return state_of(gpu, where)->hodlr_construct();
}

// 0-based block indices from the 1-based Fortran ones
std::vector<int> block_list(const int* na, const int* bl) {
  std::vector<int> out(static_cast<size_t>(std::max(*na, 0)));
  for (int a = 0; a < *na; ++a) out[static_cast<size_t>(a)] = bl[a] - 1;
  return out;
}
#endif

}  // namespace

extern "C" {

// 1 if this library has a GPU backend, else 0
void c_bpack_gpu_available(int* available) {
#ifdef BPACK_GPU_ENABLED
  *available = 1;
#else
  *available = 0;
#endif
}

// A new, empty GPU state (nullptr without a GPU backend).
void c_bpack_gpu_create(void** gpu) {
#ifdef BPACK_GPU_ENABLED
  try {
    *gpu = new State();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_create", e.what());
  }
#else
  *gpu = nullptr;
#endif
}

void c_bpack_gpu_delete(void* gpu) {
#ifdef BPACK_GPU_ENABLED
  try {
    delete static_cast<State*>(gpu);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_delete", e.what());
  }
#else
  (void)gpu;
#endif
}

// The application's evaluator of the entries (color_gpu/evaluator.hpp):
// `evaluator` points to a std::shared_ptr<fmm::gpu::Evaluator> (an empty one
// clears it), made by the registrations of C_BPACK_wrapper.cpp.
void c_bpack_gpu_set_evaluator(void* gpu, const void* evaluator) {
#ifdef BPACK_GPU_ENABLED
  if (evaluator == nullptr) gpu_fail("c_bpack_gpu_set_evaluator", "null argument");
  state_of(gpu, "c_bpack_gpu_set_evaluator")->evaluator =
      *static_cast<const std::shared_ptr<fmm::gpu::Evaluator>*>(evaluator);
#else
  (void)gpu;
  (void)evaluator;
#endif
}

// ---- HODLR (hodlr_gpu/hodlr_device.hpp); rows are local and 0-based ----

void c_bpack_hodlr_gpu_reset(void* gpu, const int64_t* n_loc, const int* maxlevel) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_reset");
    st->reset_factors();
    st->hodlr_device().reset(*n_loc, *maxlevel);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_reset", e.what());
  }
#else
  (void)gpu;
  (void)n_loc;
  (void)maxlevel;
  gpu_fail("c_bpack_hodlr_gpu_reset", "this library has no GPU backend");
#endif
}

void c_bpack_hodlr_gpu_add_leaf(void* gpu, const int64_t* row0, const int* m, const void* d, const int* ldd) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_hodlr_gpu_add_leaf")
        ->hodlr_device()
        .add_leaf(*row0, *m, static_cast<const Scalar*>(d), *ldd);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_add_leaf", e.what());
  }
#else
  (void)gpu;
  (void)row0;
  (void)m;
  (void)d;
  (void)ldd;
  gpu_fail("c_bpack_hodlr_gpu_add_leaf", "this library has no GPU backend");
#endif
}

void c_bpack_hodlr_gpu_add_lowrank(void* gpu, const int* level, const int64_t* row0, const int* m,
                                   const int64_t* col0, const int* n, const int* k, const void* u, const int* ldu,
                                   const void* v, const int* ldv, const int* sym, const void* dev_v) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_hodlr_gpu_add_lowrank")
        ->hodlr_device()
        .add_lowrank(*level, *row0, *m, *col0, *n, *k, static_cast<const Scalar*>(u), *ldu,
                     static_cast<const Scalar*>(v), *ldv, *sym != 0, static_cast<const Scalar*>(dev_v));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_add_lowrank", e.what());
  }
#else
  (void)gpu;
  (void)level;
  (void)row0;
  (void)m;
  (void)col0;
  (void)n;
  (void)k;
  (void)u;
  (void)ldu;
  (void)v;
  (void)ldv;
  (void)sym;
  (void)dev_v;
  gpu_fail("c_bpack_hodlr_gpu_add_lowrank", "this library has no GPU backend");
#endif
}

// Upload the recorded blocks; *mbytes returns their device size in MB.
void c_bpack_hodlr_gpu_commit_forward(void* gpu, double* mbytes, double* mirrored) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto& dev = state_of(gpu, "c_bpack_hodlr_gpu_commit_forward")->hodlr_device();
    dev.commit_forward();
    *mbytes = static_cast<double>(dev.forward_bytes()) / 1.0e6;
    *mirrored = static_cast<double>(dev.mirrored_bytes()) / 1.0e6;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_commit_forward", e.what());
  }
#else
  (void)gpu;
  (void)mbytes;
  (void)mirrored;
  gpu_fail("c_bpack_hodlr_gpu_commit_forward", "this library has no GPU backend");
#endif
}

// Block idx (from 0, in the order of add_lowrank) of the forward blocks on the
// device: its shape m, n, k, and its U, V copied to the host arrays hu, hv of
// leading dimensions ldu, ldv (null: not copied)
void c_bpack_hodlr_gpu_fetch_block(void* gpu, const int* idx, void* hu, const int* ldu, void* hv, const int* ldv,
                                   int* m, int* n, int* k) {
#ifdef BPACK_GPU_ENABLED
  try {
    const auto& dev = state_of(gpu, "c_bpack_hodlr_gpu_fetch_block")->hodlr_device();
    const auto& b = dev.committed_block(*idx);
    *m = b.m;
    *n = b.n;
    *k = b.k;
    dev.fetch_block(*idx, static_cast<Scalar*>(hu), *ldu, static_cast<Scalar*>(hv), *ldv);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_fetch_block", e.what());
  }
#else
  (void)gpu;
  (void)idx;
  (void)hu;
  (void)ldu;
  (void)hv;
  (void)ldv;
  (void)m;
  (void)n;
  (void)k;
  gpu_fail("c_bpack_hodlr_gpu_fetch_block", "this library has no GPU backend");
#endif
}

// Rows of U (role 0) or V (role 1) of forward block idx on the device to the
// host: out (n x k, leading dimension n) = the rows rows[0..n), numbered as
// this rank's rows from 0; k: the block's rank
void c_bpack_hodlr_gpu_gather_rows(void* gpu, const int* idx, const int* role, const int* n, const int* rows, void* out,
                                   int* k) {
#ifdef BPACK_GPU_ENABLED
  try {
    *k = state_of(gpu, "c_bpack_hodlr_gpu_gather_rows")
             ->hodlr_device()
             .gather_rows(*idx, *role, *n, rows, static_cast<Scalar*>(out));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_gather_rows", e.what());
  }
#else
  (void)gpu;
  (void)idx;
  (void)role;
  (void)n;
  (void)rows;
  (void)out;
  (void)k;
  gpu_fail("c_bpack_hodlr_gpu_gather_rows", "this library has no GPU backend");
#endif
}

// a device matrix id (DistQr::dm_wrap) viewing V of forward block idx on the
// device (n x k), 0 when it has no entries
void c_bpack_hodlr_gpu_block_v_dm(void* gpu, const int* idx, int* id) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto* st = state_of(gpu, "c_bpack_hodlr_gpu_block_v_dm");
    const auto& b = st->hodlr_device().committed_block(*idx);
    *id = b.n > 0 && b.k > 0 ? st->dist_qr().dm_wrap(b.v, b.n, b.k) : 0;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_block_v_dm", e.what());
  }
#else
  (void)gpu;
  (void)idx;
  (void)id;
  gpu_fail("c_bpack_hodlr_gpu_block_v_dm", "this library has no GPU backend");
#endif
}

// y = op(A) x, op 'N' or 'T', for the local n_loc x nrhs host arrays x, y.
// (flops: this rank's)
void c_bpack_hodlr_gpu_mult(void* gpu, const char* trans, const int* nrhs, const void* x, void* y,
                            const int* mode, double* flops) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto& dev = state_of(gpu, "c_bpack_hodlr_gpu_mult")->hodlr_device();
    dev.mult(*trans, *nrhs, static_cast<const Scalar*>(x), static_cast<Scalar*>(y), *mode);
    *flops = dev.mult_flops();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_mult", e.what());
  }
#else
  (void)gpu;
  (void)trans;
  (void)nrhs;
  (void)x;
  (void)y;
  (void)mode;
  (void)flops;
  gpu_fail("c_bpack_hodlr_gpu_mult", "this library has no GPU backend");
#endif
}

// Symmetric factorization (option%sym=1) of the uploaded HODLR.  out[0..11]:
// log|det|, phase (re, im), flops, factor MB on the device, seconds of the
// leaves, leaf propagation and nodes, largest jitter, leaves and nodes that
// needed jitter, largest rank.
void c_bpack_hodlr_gpu_factor_sym(void* gpu, const double* jitter, const int* mode, double* out) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_factor_sym");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    if (!st->sym) st->sym = std::make_unique<bpack::gpu::HodlrSymFactor<Scalar>>();
    st->sym->factor(st->hodlr_device(), *jitter);
    const auto& s = st->sym->stats();
    const double values[12] = {s.logabsdet, s.phase.real(), s.phase.imag(), s.flops, s.mbytes, s.t_leaf,
                               s.t_prop, s.t_nodes, s.max_jitter, static_cast<double>(s.jitter_leaves),
                               static_cast<double>(s.jitter_nodes), static_cast<double>(s.maxrank)};
    for (int i = 0; i < 12; ++i) out[i] = values[i];
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_factor_sym", e.what());
  }
#else
  (void)gpu;
  (void)jitter;
  (void)mode;
  (void)out;
  gpu_fail("c_bpack_hodlr_gpu_factor_sym", "this library has no GPU backend");
#endif
}

// 1 if the symmetric factors are on the device, else 0
void c_bpack_hodlr_gpu_sym_ready(void* gpu, int* ready) {
#ifdef BPACK_GPU_ENABLED
  State* st = static_cast<State*>(gpu);
  *ready = (st != nullptr && st->sym && st->sym->ready()) ? 1 : 0;
#else
  (void)gpu;
  *ready = 0;
#endif
}

// y = A^{-1} x with the symmetric factors, for the local n_loc x nrhs host arrays
void c_bpack_hodlr_gpu_solve_sym(void* gpu, const int* nrhs, const void* x, void* y, const int* mode, double* flops) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_solve_sym");
    if (!st->sym) gpu_fail("c_bpack_hodlr_gpu_solve_sym", "no symmetric factorization on the device");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    st->sym->solve(*nrhs, static_cast<const Scalar*>(x), static_cast<Scalar*>(y));
    *flops = st->sym->solve_flops();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_solve_sym", e.what());
  }
#else
  (void)gpu;
  (void)nrhs;
  (void)x;
  (void)y;
  (void)mode;
  (void)flops;
  gpu_fail("c_bpack_hodlr_gpu_solve_sym", "this library has no GPU backend");
#endif
}

// Copies of the symmetric factor of the node of `level` whose child 0 starts
// at local row c0 (0-based): Q0, Q1, G0, G1, the LU of S and its pivots;
// *k its rank, *found 0 if there is no such node.
void c_bpack_hodlr_gpu_download_symnode(void* gpu, const int* level, const int64_t* c0, const int* n0,
                                        const int* n1, void* q0, void* q1, void* g0, void* g1, void* s, int* ipiv,
                                        int* k, int* found, void* z0, void* z1) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_download_symnode");
    *found = 0;
    if (st->sym && st->sym->download_node(*level, *c0, *n0, *n1, static_cast<Scalar*>(q0), static_cast<Scalar*>(q1),
                                          static_cast<Scalar*>(g0), static_cast<Scalar*>(g1),
                                          static_cast<Scalar*>(s), ipiv, k, static_cast<Scalar*>(z0),
                                          static_cast<Scalar*>(z1))) {
      *found = 1;
    }
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_download_symnode", e.what());
  }
#else
  (void)gpu;
  (void)level;
  (void)c0;
  (void)n0;
  (void)n1;
  (void)q0;
  (void)q1;
  (void)g0;
  (void)g1;
  (void)s;
  (void)ipiv;
  (void)k;
  (void)found;
  (void)z0;
  (void)z1;
  gpu_fail("c_bpack_hodlr_gpu_download_symnode", "this library has no GPU backend");
#endif
}

// Unsymmetric factorization (option%sym=0) of the uploaded HODLR.  out[0..11]:
// log|det|, phase (re, im), flops, factor MB on the device, seconds of the
// leaves, of the updates of the blocks above each level (Sblock) and of the
// node inverses, 0, leaves and nodes whose LU had pivots raised to the
// threshold, largest rank.
void c_bpack_hodlr_gpu_factor_unsym(void* gpu, const double* jitter, const int* mode, double* out) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_factor_unsym");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    if (!st->unsym) st->unsym = std::make_unique<bpack::gpu::HodlrUnsymFactor<Scalar>>();
    st->unsym->factor(st->hodlr_device(), *jitter);
    const auto& s = st->unsym->stats();
    const double values[12] = {s.logabsdet, s.phase.real(), s.phase.imag(), s.flops, s.mbytes, s.t_leaf,
                               s.t_sblock, s.t_nodes, 0.0, static_cast<double>(s.raised_leaves),
                               static_cast<double>(s.raised_nodes), static_cast<double>(s.maxrank)};
    for (int i = 0; i < 12; ++i) out[i] = values[i];
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_factor_unsym", e.what());
  }
#else
  (void)gpu;
  (void)jitter;
  (void)mode;
  (void)out;
  gpu_fail("c_bpack_hodlr_gpu_factor_unsym", "this library has no GPU backend");
#endif
}

// 1 if the unsymmetric factors are on the device, else 0
void c_bpack_hodlr_gpu_unsym_ready(void* gpu, int* ready) {
#ifdef BPACK_GPU_ENABLED
  State* st = static_cast<State*>(gpu);
  *ready = (st != nullptr && st->unsym && st->unsym->ready()) ? 1 : 0;
#else
  (void)gpu;
  *ready = 0;
#endif
}

// y = op(A)^{-1} x, op 'N' or 'T', with the unsymmetric factors, for the
// local n_loc x nrhs host arrays
void c_bpack_hodlr_gpu_solve_unsym(void* gpu, const char* trans, const int* nrhs, const void* x, void* y,
                                   const int* mode, double* flops) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_solve_unsym");
    if (!st->unsym) gpu_fail("c_bpack_hodlr_gpu_solve_unsym", "no unsymmetric factorization on the device");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    st->unsym->solve(*trans, *nrhs, static_cast<const Scalar*>(x), static_cast<Scalar*>(y));
    *flops = st->unsym->solve_flops();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_solve_unsym", e.what());
  }
#else
  (void)gpu;
  (void)trans;
  (void)nrhs;
  (void)x;
  (void)y;
  (void)mode;
  (void)flops;
  gpu_fail("c_bpack_hodlr_gpu_solve_unsym", "this library has no GPU backend");
#endif
}

// Copies of the unsymmetric factors for the comparison with the CPU (local,
// 0-based rows; *found 0 if there is no such block): what = 0: U' of the
// block of `level` at rows row0, columns col0 (m x k); 1: Uo of the node of
// `level` whose child 0 starts at row0 (m x k); 2: the inverse of the leaf
// at row0 (m x m).
void c_bpack_hodlr_gpu_download_unsym(void* gpu, const int* what, const int* level, const int64_t* row0,
                                      const int64_t* col0, const int* m, const int* k, void* out, int* found) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_download_unsym");
    Scalar* o = static_cast<Scalar*>(out);
    bool ok = false;
    if (st->unsym) {
      if (*what == 0) ok = st->unsym->download_block(*level, *row0, *col0, *m, *k, o);
      if (*what == 1) ok = st->unsym->download_schur(*level, *row0, *m, *k, o);
      if (*what == 2) ok = st->unsym->download_leaf_inverse(*row0, *m, o);
    }
    *found = ok ? 1 : 0;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_download_unsym", e.what());
  }
#else
  (void)gpu;
  (void)what;
  (void)level;
  (void)row0;
  (void)col0;
  (void)m;
  (void)k;
  (void)out;
  (void)found;
  gpu_fail("c_bpack_hodlr_gpu_download_unsym", "this library has no GPU backend");
#endif
}

// Level `level` of the uploaded HODLR is shared by the ranks of the Fortran
// communicator comm (the node of this rank at that level)
void c_bpack_hodlr_gpu_add_shared_level(void* gpu, const int* level, const MPI_Fint* comm) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_hodlr_gpu_add_shared_level")->hodlr_device().add_shared_level(*level, MPI_Comm_f2c(*comm));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_add_shared_level", e.what());
  }
#else
  (void)gpu;
  (void)level;
  (void)comm;
#endif
}

// Collective over the Fortran communicator comm (once per run): the ranks
// of a node that use the same GPU split its memory (fmm::gpu::share_device);
// *share: how many ranks use this rank's GPU.  Also the report of the
// environment variables (h2_parallel/bpack_env.hpp).
void c_bpack_gpu_share_device(const MPI_Fint* comm_f, int* share) {
  MPI_Comm comm = MPI_Comm_f2c(*comm_f);
  fmm::env::report_environment(comm);
#ifdef BPACK_GPU_ENABLED
  try {
    *share = fmm::gpu::share_device(comm);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_share_device", e.what());
  }
#else
  *share = 1;
#endif
}

// Multi-rank runs with GPU-aware MPI: the device buffers MPI reads and
// writes come from the exchange arena, reserved before the main device heap
// (fmm::gpu::reserve_exchange_arena: BPACK_GPU_EXCHANGE_MB; the merges of
// the HODLR construction move up to a few hundred MB per message)
void c_bpack_gpu_init_exchange(const int* nproc, const int* share) {
#ifdef BPACK_GPU_ENABLED
  try {
    fmm::gpu::reserve_exchange_arena(*nproc, *share);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_init_exchange", e.what());
  }
#else
  (void)nproc;
  (void)share;
#endif
}

// BPACK_CHECK and BPACK_TRACE for the Fortran code: 1 if the list names
// `value` (a null-terminated string; doc/environment_variables.md)
int c_bpack_env_check(const char* value) { return fmm::env::check(value) ? 1 : 0; }
int c_bpack_env_trace(const char* value) { return fmm::env::trace(value) ? 1 : 0; }

// ---- HODLR construction (hodlr_gpu/hodlr_construct.hpp); point slots are
// tree indices - 1, block indices and rows / columns 1-based ----

// Before the matrix's first GPU work, not in the construction time
// (GpuState::warm_up): the factorization's kernels, with *construct the GPU
// construction's routines and handles, with *evaluate the evaluator;
// seconds[0] the first two, seconds[1] the evaluator
void c_bpack_hodlr_gpu_warm_up(void* gpu, const int* construct, const int* evaluate, double* seconds) {
  seconds[0] = seconds[1] = 0.0;
#ifdef BPACK_GPU_ENABLED
  try {
    const std::pair<double, double> s =
        state_of(gpu, "c_bpack_hodlr_gpu_warm_up")->warm_up(*construct != 0, *evaluate != 0);
    seconds[0] = s.first;
    seconds[1] = s.second;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_warm_up", e.what());
  }
#else
  (void)gpu;
  (void)construct;
  (void)evaluate;
#endif
}

// 1 if the matrix has a GPU evaluator of the entries and its points
void c_bpack_hodlr_gpu_construct_ready(void* gpu, int* ready) {
  *ready = 0;
#ifdef BPACK_GPU_ENABLED
  if (gpu == nullptr) return;
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_construct_ready");
    if (st->evaluator && st->construct && st->construct->has_points()) {
      if (st->evaluator->scalar_bytes() != static_cast<int>(sizeof(Scalar))) {
        gpu_fail("c_bpack_hodlr_gpu_construct_ready", "the GPU evaluator's value type does not match the matrix");
      }
      if (st->evaluator->needs_coordinates() && st->construct->dim() == 0) {
        gpu_fail("c_bpack_hodlr_gpu_construct_ready",
                 "the GPU evaluator reads coordinates (BPACK_GPU_COORDINATES) but the matrix has none");
      }
      *ready = 1;
    }
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_construct_ready", e.what());
  }
#else
  (void)gpu;
#endif
}

// The points in tree order: xyz (dim x n; dim 0: none) and the 0-based
// original indices
void c_bpack_hodlr_gpu_set_points(void* gpu, const int64_t* n, const int* dim, const double* xyz,
                                  const int64_t* ids) {
#ifdef BPACK_GPU_ENABLED
  try {
    construct_of(gpu, "c_bpack_hodlr_gpu_set_points").set_points(*n, *dim, xyz, ids);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_set_points", e.what());
  }
#else
  (void)gpu;
  (void)n;
  (void)dim;
  (void)xyz;
  (void)ids;
#endif
}

// out[b] (m[b] x n[b], host) = scale * A(r0[b] .., c0[b] ..), r0 and c0 0-based slots
void c_bpack_hodlr_gpu_eval_dense(void* gpu, const double* scale, const int* count, const int64_t* r0,
                                  const int* m, const int64_t* c0, const int* n, void** out) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_eval_dense");
    if (!st->evaluator) throw std::logic_error("no GPU evaluator registered");
    st->hodlr_construct().eval_dense(*st->evaluator, *scale, *count, r0, m, c0, n,
                                     reinterpret_cast<Scalar* const*>(out));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_eval_dense", e.what());
  }
#else
  (void)gpu;
  (void)scale;
  (void)count;
  (void)r0;
  (void)m;
  (void)c0;
  (void)n;
  (void)out;
  gpu_fail("c_bpack_hodlr_gpu_eval_dense", "this library has no GPU backend");
#endif
}

// variant: option%RecLR_leaf, 4 (BACA) or 5 (BACA without overlap)
void c_bpack_hodlr_gpu_baca_begin(void* gpu, const double* scale, const int* mode, const int* variant,
                                  const int* nb, const int64_t* r0, const int* m, const int64_t* c0, const int* n,
                                  const int* r_est) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_baca_begin");
    if (!st->evaluator) throw std::logic_error("no GPU evaluator registered");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    st->hodlr_construct().begin(*st->evaluator, *scale, *variant, *nb, r0, m, c0, n,
                                r_est);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_begin", e.what());
  }
#else
  (void)gpu;
  (void)scale;
  (void)mode;
  (void)variant;
  (void)nb;
  (void)r0;
  (void)m;
  (void)c0;
  (void)n;
  (void)r_est;
  gpu_fail("c_bpack_hodlr_gpu_baca_begin", "this library has no GPU backend");
#endif
}

void c_bpack_hodlr_gpu_baca_set_columns(void* gpu, const int* b, const int* cols) {
#ifdef BPACK_GPU_ENABLED
  try {
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_set_columns").set_columns(*b - 1, cols);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_set_columns", e.what());
  }
#else
  (void)gpu;
  (void)b;
  (void)cols;
#endif
}

void c_bpack_hodlr_gpu_baca_panels(void* gpu, const int* na, const int* bl, void* core) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_panels").panels(*na, list.data(), static_cast<Scalar*>(core));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_panels", e.what());
  }
#else
  (void)gpu;
  (void)na;
  (void)bl;
  (void)core;
#endif
}

void c_bpack_hodlr_gpu_baca_knn_panels(void* gpu, const int* na, const int* bl, const int* nc, const int* cols,
                                       const int* nr, const int* rows, void* core) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_knn_panels")
        .knn_panels(*na, list.data(), nc, cols, nr, rows, static_cast<Scalar*>(core));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_knn_panels", e.what());
  }
#else
  (void)gpu;
  (void)na;
  (void)bl;
  (void)nc;
  (void)cols;
  (void)nr;
  (void)rows;
  (void)core;
#endif
}

void c_bpack_hodlr_gpu_baca_append(void* gpu, const int* knn, const int* na, const int* bl, const int* ru,
                                   const int* jpvt, const void* w, const int* rskip, void* grams) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_append")
        .append(*knn, *na, list.data(), ru, jpvt, static_cast<const Scalar*>(w), rskip, static_cast<Scalar*>(grams));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_append", e.what());
  }
#else
  (void)gpu;
  (void)knn;
  (void)na;
  (void)bl;
  (void)ru;
  (void)jpvt;
  (void)w;
  (void)rskip;
  (void)grams;
#endif
}

// The recompression of the listed blocks after their BACA, all on the GPU
// (HodlrConstruct::recompress): rn the new ranks; the new factors stay for
// c_bpack_hodlr_gpu_baca_download
void c_bpack_hodlr_gpu_baca_recompress(void* gpu, const int* na, const int* bl, const double* tol, const double* underflow,
                                       int* rn) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    auto* st = state_of(gpu, "c_bpack_hodlr_gpu_baca_recompress");
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_recompress")
        .recompress(*na, list.data(), *tol, *underflow, st->device_svd(), rn);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_recompress", e.what());
  }
#else
  (void)gpu;
  (void)na;
  (void)bl;
  (void)tol;
  (void)underflow;
  (void)rn;
  gpu_fail("c_bpack_hodlr_gpu_baca_recompress", "this library has no GPU backend");
#endif
}

// the new factors of c_bpack_hodlr_gpu_baca_recompress as device matrices of
// the merges (DistQr::dm_*): ids[2a], ids[2a + 1] those of block a's U, V
void c_bpack_hodlr_gpu_baca_to_dm(void* gpu, const int* na, const int* bl, const int* rn, int* ids) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    auto* st = state_of(gpu, "c_bpack_hodlr_gpu_baca_to_dm");
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_to_dm").to_dm(*na, list.data(), rn, st->dist_qr(), ids);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_to_dm", e.what());
  }
#else
  (void)gpu;
  (void)na;
  (void)bl;
  (void)rn;
  (void)ids;
  gpu_fail("c_bpack_hodlr_gpu_baca_to_dm", "this library has no GPU backend");
#endif
}

// the new factors of c_bpack_hodlr_gpu_baca_recompress to the host: u[a] (m x rn), v[a] (n x rn);
// keep: their device copies stay as the device HODLR's mirrors of these arrays
void c_bpack_hodlr_gpu_baca_download(void* gpu, const int* na, const int* bl, const int* rn, void** u, void** v,
                                     const int* keep) {
#ifdef BPACK_GPU_ENABLED
  try {
    const std::vector<int> list = block_list(na, bl);
    auto* st = state_of(gpu, "c_bpack_hodlr_gpu_baca_download");
    construct_of(gpu, "c_bpack_hodlr_gpu_baca_download")
        .download(*na, list.data(), rn, reinterpret_cast<Scalar* const*>(u), reinterpret_cast<Scalar* const*>(v),
                  *keep > 0 ? &st->hodlr_device() : nullptr, *keep < 2);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_download", e.what());
  }
#else
  (void)gpu;
  (void)na;
  (void)bl;
  (void)rn;
  (void)u;
  (void)v;
  (void)keep;
  gpu_fail("c_bpack_hodlr_gpu_baca_download", "this library has no GPU backend");
#endif
}

// ---- the TSQR of the merges of shared blocks (hodlr_gpu/hodlr_distqr.hpp) ----

// ---- the device matrices of the merges of shared blocks (DistQr::dm_*) ----

// a device matrix id viewing the device HODLR's mirror of host (rows x cols), or 0 without one
void c_bpack_gpu_dm_from_mirror(void* gpu, const void* host, const int* rows, const int* cols, int* id) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto* st = state_of(gpu, "c_bpack_gpu_dm_from_mirror");
    const size_t bytes = static_cast<size_t>(std::max(*rows, 0)) * std::max(*cols, 0) * sizeof(Scalar);
    const void* d = st->hodlr_device().mirror_ptr(host, bytes);
    *id = d != nullptr ? st->dist_qr().dm_wrap(static_cast<const Scalar*>(d), *rows, *cols) : 0;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_from_mirror", e.what());
  }
#else
  (void)gpu;
  (void)host;
  (void)rows;
  (void)cols;
  (void)id;
  gpu_fail("c_bpack_gpu_dm_from_mirror", "this library has no GPU backend");
#endif
}

// the flops of the merges (DistQr: TSQR and products) since the last call, this rank's
void c_bpack_gpu_dm_flops(void* gpu, double* flops) {
#ifdef BPACK_GPU_ENABLED
  try {
    *flops = state_of(gpu, "c_bpack_gpu_dm_flops")->dist_qr().take_flops();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_flops", e.what());
  }
#else
  (void)gpu;
  *flops = 0.0;
#endif
}

// device matrix id's memory to the device HODLR (freed after its commit); dev: its address
void c_bpack_gpu_dm_give(void* gpu, const int* id, void** dev) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto* st = state_of(gpu, "c_bpack_gpu_dm_give");
    Scalar* d = st->dist_qr().dm_give(*id);
    st->hodlr_device().own_mirror_memory(d);
    *dev = d;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_give", e.what());
  }
#else
  (void)gpu;
  (void)id;
  (void)dev;
  gpu_fail("c_bpack_gpu_dm_give", "this library has no GPU backend");
#endif
}

// device matrix id, whose copy on the host is host: kept as the device HODLR's mirror of it
// (stale != 0: host was left unfilled)
void c_bpack_gpu_dm_keep(void* gpu, const int* id, const void* host, const int* stale) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto* st = state_of(gpu, "c_bpack_gpu_dm_keep");
    st->dist_qr().dm_keep(*id, static_cast<const Scalar*>(host), st->hodlr_device(), *stale != 0);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_keep", e.what());
  }
#else
  (void)gpu;
  (void)id;
  (void)host;
  (void)stale;
  gpu_fail("c_bpack_gpu_dm_keep", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_alloc(void* gpu, const int* rows, const int* cols, const void* host, int* id) {
#ifdef BPACK_GPU_ENABLED
  try {
    *id = state_of(gpu, "c_bpack_gpu_dm_alloc")->dist_qr().dm_alloc(*rows, *cols, static_cast<const Scalar*>(host));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_alloc", e.what());
  }
#else
  (void)gpu;
  (void)rows;
  (void)cols;
  (void)host;
  (void)id;
  gpu_fail("c_bpack_gpu_dm_alloc", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_download(void* gpu, const int* id, void* host) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_gpu_dm_download")->dist_qr().dm_download(*id, static_cast<Scalar*>(host));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_download", e.what());
  }
#else
  (void)gpu;
  (void)id;
  (void)host;
  gpu_fail("c_bpack_gpu_dm_download", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_free(void* gpu, const int* id) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_gpu_dm_free")->dist_qr().dm_free(*id);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_free", e.what());
  }
#else
  (void)gpu;
  (void)id;
  gpu_fail("c_bpack_gpu_dm_free", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_redistribute(void* gpu, const int* comm, const int* src, const int* ncols, const int* dst, const int* c0, const int* nsend, const int* s_rank, const int* s_off, const int* s_rows, const int* nrecv, const int* r_rank, const int* r_off, const int* r_rows, const int* tag) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_gpu_dm_redistribute")
        ->dist_qr()
        .dm_redistribute(MPI_Comm_f2c(*comm), *src, *ncols, *dst, *c0, *nsend, s_rank, s_off, s_rows, *nrecv, r_rank, r_off,
                         r_rows, *tag);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_redistribute", e.what());
  }
#else
  (void)gpu;
  (void)comm;
  (void)src;
  (void)ncols;
  (void)dst;
  (void)c0;
  (void)nsend;
  (void)s_rank;
  (void)s_off;
  (void)s_rows;
  (void)nrecv;
  (void)r_rank;
  (void)r_off;
  (void)r_rows;
  (void)tag;
  gpu_fail("c_bpack_gpu_dm_redistribute", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_tsqr(void* gpu, const int* comm, const int* x, const double* tol, const double* underflow, int* rn, double* tsec, int* xnew) {
#ifdef BPACK_GPU_ENABLED
  try {
    *rn = state_of(gpu, "c_bpack_gpu_dm_tsqr")
              ->dist_qr()
              .dm_tsqr(MPI_Comm_f2c(*comm), *x, *tol, *underflow, tsec, xnew);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_tsqr", e.what());
  }
#else
  (void)gpu;
  (void)comm;
  (void)x;
  (void)tol;
  (void)underflow;
  (void)rn;
  (void)tsec;
  (void)xnew;
  gpu_fail("c_bpack_gpu_dm_tsqr", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_fetch_sw(void* gpu, void* sw) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_gpu_dm_fetch_sw")->dist_qr().dm_fetch_sw(static_cast<Scalar*>(sw));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_fetch_sw", e.what());
  }
#else
  (void)gpu;
  (void)sw;
  gpu_fail("c_bpack_gpu_dm_fetch_sw", "this library has no GPU backend");
#endif
}

void c_bpack_gpu_dm_gemm_nt(void* gpu, const int* a, const int* m, const int* k, const void* b, const int* n, const int* c, const int* c_row0) {
#ifdef BPACK_GPU_ENABLED
  try {
    state_of(gpu, "c_bpack_gpu_dm_gemm_nt")
        ->dist_qr()
        .dm_gemm_nt(*a, *m, *k, static_cast<const Scalar*>(b), *n, *c, *c_row0);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_dm_gemm_nt", e.what());
  }
#else
  (void)gpu;
  (void)a;
  (void)m;
  (void)k;
  (void)b;
  (void)n;
  (void)c;
  (void)c_row0;
  gpu_fail("c_bpack_gpu_dm_gemm_nt", "this library has no GPU backend");
#endif
}

// End the BACA of a level; out[0..4]: seconds of the dense leaves, panels,
// appends and recompression, and device flops, since the last call
void c_bpack_hodlr_gpu_baca_end(void* gpu, double* out) {
#ifdef BPACK_GPU_ENABLED
  try {
    auto& c = construct_of(gpu, "c_bpack_hodlr_gpu_baca_end");
    c.end();
    const auto& s = c.stats();
    out[0] = s.t_eval;
    out[1] = s.t_panels;
    out[2] = s.t_append;
    out[3] = s.t_recomp;
    out[4] = s.flops;
    out[5] = s.t_entry;
    c.reset_stats();
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_baca_end", e.what());
  }
#else
  (void)gpu;
  for (int i = 0; i < 6; ++i) out[i] = 0.0;
#endif
}

}  // extern "C"
