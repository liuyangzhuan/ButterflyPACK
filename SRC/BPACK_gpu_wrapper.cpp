// C entry points of the GPU backend of the formats other than H2 (HODLR for
// now; hodlr_gpu/), called from Fortran (BPACK_gpu.f90) and from
// C_BPACK_wrapper.cpp.  PrecisionPreprocessing.sh copies this file into each
// precision; the double and double complex libraries of a GPU build
// (enable_h2_gpu) get the backend, the others stubs that report its absence.

#ifdef HAVE_MPI
#include <mpi.h>
#else
// A build without MPI (enable_mpi=OFF) only has the stubs below; the Fortran
// communicator arguments are plain integers there, as in BPACK_wrapper.h.
#ifndef MPI_Fint
#define MPI_Fint int
#endif
#endif

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "gpu_common/bpack_env.hpp"

#if defined(H2_HAVE_GPU) && defined(DAT) && (DAT == 0 || DAT == 1)
#ifndef HAVE_MPI
#error "ButterflyPACK's GPU backend needs a build with MPI support (enable_mpi=ON)"
#endif
#define BPACK_GPU_ENABLED 1
#include "hodlr_gpu/gpu_state.hpp"
#if DAT == 1
// the butterfly backend is real double for now (bf_gpu's DataType)
#define BPACK_BF_GPU_ENABLED 1
#include "bf_gpu/bf_id.hpp"
#include "bf_gpu/bf_randomized.hpp"
#include "bf_gpu/bf_store.hpp"
#endif
#endif

#include "bpack_gpu.h"  // bpack_gpu_bf_matvec (GPU_INTERFACE)

#ifdef BPACK_BF_GPU_ENABLED
extern "C" {
void c_bf_get_layout(void* bf, int* level_butterfly, int* level_half, int* M, int* N, int* M_loc,
                     int* N_loc, int* nlevels, int* nblocks);
void c_bf_get_blocks(void* bf, int* nr, int* nc, int* idx_r, int* idx_c, int* inc_r, int* inc_c,
                     int* rows, int* cols, void** ptrs);
void c_bf_get_vec_ranges(void* bf, void* ptree, char dir, int* idx_r, int* inc_r, int* nr, int* idx_c,
                         int* inc_c, int* nc);
void c_bf_get_mvp_vec_ranges(void* bf, void* ptree, char chara, int* idx_r, int* inc_r, int* nr,
                             int* idx_c, int* inc_c, int* nc);
void c_bf_get_sub_head(void* bf, void* ptree, int* is_head, int* comm, int* nproc);
void c_bf_get_exchange_dest(void* bf, void* ptree, int level, char mode, char collect, int nr, int nc,
                            int idx_r, int inc_r, int idx_c, int inc_c, int* pid);
void c_bf_get_all2all_dest(void* bf, void* ptree, int level, char mode_new, int nr, int nc, int idx_r,
                           int inc_r, int idx_c, int inc_c, int* pid);
void c_bf_get_all2all_range(void* bf, void* ptree, int level, char mode_new, int* idx_r, int* inc_r,
                            int* nr, int* idx_c, int* inc_c, int* nc);
void c_bf_get_block_group(void* bf, void* ptree, int level, int index_i, int index_j, char mode,
                          int* nproc, int* comm, int* is_head);
}
#endif

namespace {

[[noreturn]] void gpu_fail(const char* where, const char* what) {
  std::fprintf(stderr, "ButterflyPACK GPU backend error in %s: %s\n", where, what);
  std::fflush(stderr);
#ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
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

// The application's evaluator of the entries (h2_gpu/evaluator.hpp):
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

// The rank of forward block idx on the device (after the commit)
void c_bpack_hodlr_gpu_block_rank(void* gpu, const int* idx, int* k) {
#ifdef BPACK_GPU_ENABLED
  try {
    *k = state_of(gpu, "c_bpack_hodlr_gpu_block_rank")->hodlr_device().committed_block(*idx).k;
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_block_rank", e.what());
  }
#else
  (void)gpu;
  (void)idx;
  (void)k;
  gpu_fail("c_bpack_hodlr_gpu_block_rank", "this library has no GPU backend");
#endif
}

// The products of entry extraction on the device (HodlrDevice::extract_products):
// np pieces, each rows of U (role 0) or V (role 1) of forward block idx times a
// smaller factor block from the host; flops: their flops
void c_bpack_hodlr_gpu_extract_products(void* gpu, const int* np, const int* idx, const int* role, const int* nbig,
                                        const int* nsmall, const int* small_first, const int* big_rows,
                                        const void* small, void* out, double* flops) {
#ifdef BPACK_GPU_ENABLED
  try {
    *flops = state_of(gpu, "c_bpack_hodlr_gpu_extract_products")
                 ->hodlr_device()
                 .extract_products(*np, idx, role, nbig, nsmall, small_first, big_rows,
                                   static_cast<const Scalar*>(small), static_cast<Scalar*>(out));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_extract_products", e.what());
  }
#else
  (void)gpu;
  (void)np;
  (void)idx;
  (void)role;
  (void)nbig;
  (void)nsmall;
  (void)small_first;
  (void)big_rows;
  (void)small;
  (void)out;
  (void)flops;
  gpu_fail("c_bpack_hodlr_gpu_extract_products", "this library has no GPU backend");
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
// environment variables (GPU_BACKEND/gpu_common/bpack_env.hpp).
void c_bpack_gpu_share_device(const MPI_Fint* comm_f, int* share) {
#ifdef HAVE_MPI
  MPI_Comm comm = MPI_Comm_f2c(*comm_f);
  fmm::env::report_environment(comm);
#else
  (void)comm_f;
  fmm::env::report_environment_as(0);  // a build without MPI is one rank
#endif
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

// 1 if the matrix has a GPU evaluator of the entries and its points; 2 if it
// has its points but no GPU evaluator (the GPU construction then takes the
// entries through the host: c_bpack_hodlr_gpu_baca_begin with host_fn)
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
    } else if (!st->evaluator && st->construct && st->construct->has_points()) {
      *ready = 2;
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

// The low-rank term U V^T (n x n, tree order) the construction adds to the entries it takes from the host
// (c_bpack_set_entry_addend): the nu rows of U and nv of V this rank holds (r columns, column-major), and the
// 1-based slot of every row there (umap, vmap, n each, 0: not held); n = 0 frees it
void c_bpack_hodlr_gpu_set_addend(void* gpu, const int64_t* n, const int* nu, const int* nv, const int* r,
                                  const void* u, const void* v, const int* umap, const int* vmap) {
#ifdef BPACK_GPU_ENABLED
  try {
    construct_of(gpu, "c_bpack_hodlr_gpu_set_addend")
        .set_addend(*n, *nu, *nv, *r, reinterpret_cast<const Scalar*>(u), reinterpret_cast<const Scalar*>(v), umap,
                    vmap);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_hodlr_gpu_set_addend", e.what());
  }
#else
  (void)gpu;
  (void)n;
  (void)nu;
  (void)nv;
  (void)r;
  (void)u;
  (void)v;
  (void)umap;
  (void)vmap;
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

// variant: option%RecLR_leaf, 4 (BACA) or 5 (BACA without overlap).  The
// entries from the matrix's GPU evaluator, or without one through the host:
// host_fn (c_bpack_hodlr_gpu_host_entries, HodlrConstruct::HostEntryFn) with
// host_ctx
void c_bpack_hodlr_gpu_baca_begin(void* gpu, const double* scale, const int* mode, const int* variant,
                                  const int* nb, const int64_t* r0, const int* m, const int64_t* c0, const int* n,
                                  const int* r_est, void* host_fn, void* host_ctx) {
#ifdef BPACK_GPU_ENABLED
  try {
    State* st = state_of(gpu, "c_bpack_hodlr_gpu_baca_begin");
    const auto fn = reinterpret_cast<bpack::gpu::HodlrConstruct<Scalar>::HostEntryFn>(host_fn);
    if (!st->evaluator && fn == nullptr) throw std::logic_error("no GPU evaluator registered and no host entries");
    fmm::gpu::tensor_core_gemm() = (*mode == 2);
    st->hodlr_construct().begin(st->evaluator ? st->evaluator.get() : nullptr, fn, host_ctx, *scale, *variant, *nb, r0,
                                m, c0, n, r_est);
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
  (void)host_fn;
  (void)host_ctx;
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

#ifdef BPACK_BF_GPU_ENABLED
namespace {

// The device's setup for a butterfly's ranks, before the device pool is
// first reserved: ranks of a node that share a GPU split its memory
// (fmm::gpu::share_device), and a run of several ranks takes its exchange
// arena.  Collective over comm.  Both take effect once per process, as for an
// H2 (begin_operator_build) and a HODLR (BPACK_gpu_bind).
void bf_device_setup(MPI_Comm comm, int nproc) {
  const int share = fmm::gpu::share_device(comm);
  fmm::gpu::reserve_exchange_arena(nproc, share);
}

// A butterfly's GPU evaluator of its entries (c_bf_set_gpu_*): the
// application's evaluator, and the points it is evaluated over - the
// butterfly's own index order, rows 0..M-1 then columns M..M+N-1, with ids
// the CPU callbacks' indices minus one (a row's in the row mesh's order, a
// column's in the column mesh's) and the coordinates of a mesh with geometry.
struct BfEntries {
  std::shared_ptr<fmm::gpu::Evaluator> evaluator;
  int M = -1, N = -1, dim = 0;
  int64_t* d_ids = nullptr;
  double* d_xyz = nullptr;
  fmm::gpu::MetaBuilder meta;
  fmm::gpu::DeviceBuffer meta_device;
  double warm_seconds = 0.0;
  // a level's entries an evaluation left on the device (c_bf_gpu_eval_blocks
  // with keep), for the rank-reveal that follows (c_bf_gpu_level_id_kept)
  double* kept = nullptr;
  int64_t kept_elements = 0;

  ~BfEntries() {
    drop_kept();
    drop_points();
  }
  void drop_kept() {
    if (kept != nullptr) fmm::gpu::DeviceHeap::instance().free(kept);
    kept = nullptr;
    kept_elements = 0;
  }
  void drop_points() {
    auto& heap = fmm::gpu::DeviceHeap::instance();
    if (d_ids != nullptr) heap.free(d_ids);
    if (d_xyz != nullptr) heap.free(d_xyz);
    d_ids = nullptr;
    d_xyz = nullptr;
    M = N = -1;
    dim = 0;
  }
  fmm::gpu::PointTable points() const {
    fmm::gpu::PointTable p;
    p.xyz = d_xyz;
    p.ids = d_ids;
    p.dim = std::max(dim, 1);
    return p;
  }
};

std::map<const void*, std::unique_ptr<BfEntries>>& bf_entries() {
  static std::map<const void*, std::unique_ptr<BfEntries>> registry;
  return registry;
}

// A butterfly's black box on the device (c_bf_set_gpu_matvec), for its
// construction from matvecs (c_bf_gpu_randomized).
struct BfMatvec {
  bpack_gpu_bf_matvec fn = nullptr;
  void* user = nullptr;
};

std::map<const void*, BfMatvec>& bf_matvecs() {
  static std::map<const void*, BfMatvec> registry;
  return registry;
}

// The C callback matvec of C_BF_Construct_Matvec_Compute, called as
// Bmatvec_user_C calls it: (trans, Nin, Nout, nvec, Vin, Vout, quant, a, b).
using BfHostMatvec = void (*)(const char*, int*, int*, int*, const double*, double*, void*, double*, double*);

// A construction from matvecs between c_bf_gpu_randomized and its fetch, and
// the pinned staging of a host black box.
struct BfRandomizedRun {
  std::unique_ptr<bf::gpu::Randomized> run;
  fmm::gpu::PinnedBuffer hx, hy;
};

// Make a butterfly resident: upload its factors, read the grids the sweep
// needs, and wire the closures through which bf_gpu asks the Fortran for its
// distribution decisions (bf_ranges.hpp explains why they are closures).
// Returns nullptr when the device has no room for the factors.  Collective
// over the butterfly's ranks the first time, as its matvec is.
bf::gpu::BfResident* bf_resident(void* bf, void* ptree) {
  bf::gpu::BfFingerprint print;
  c_bf_get_layout(bf, &print.level_butterfly, &print.level_half, &print.M, &print.N, &print.M_loc,
                  &print.N_loc, &print.nlevels, &print.nblocks);
  if (bf::gpu::BfResident* hit = bf::gpu::bf_store_find(bf, print)) return hit;
  if (print.nlevels <= 0 || print.nblocks <= 0) return nullptr;

  // the factor layout and the host pointers of every block
  const int nlev = print.nlevels, nblk = print.nblocks;
  std::vector<int> nr(nlev), nc(nlev), ir(nlev), ic(nlev), er(nlev), ec(nlev);
  std::vector<int> brows(nblk), bcols(nblk);
  std::vector<void*> bptr(nblk);
  c_bf_get_blocks(bf, nr.data(), nc.data(), ir.data(), ic.data(), er.data(), ec.data(), brows.data(),
                  bcols.data(), bptr.data());

  std::vector<bf::gpu::LevelDesc> lv(nlev);
  std::vector<std::vector<const double*>> up(nlev);
  std::vector<std::vector<std::pair<int, int>>> shp(nlev);
  int at = 0;
  for (int l = 0; l < nlev; ++l) {
    lv[l].nr = nr[l]; lv[l].nc = nc[l];
    lv[l].idx_r = ir[l]; lv[l].idx_c = ic[l];
    lv[l].inc_r = er[l]; lv[l].inc_c = ec[l];
    const int want = nr[l] * nc[l];
    for (int b = 0; b < want; ++b) {
      up[l].push_back(reinterpret_cast<const double*>(bptr[at + b]));
      shp[l].push_back({brows[at + b], bcols[at + b]});
    }
    at += want;
  }

  int is_head = 0, comm_f = 0, nproc = 1, myid = 0;
  c_bf_get_sub_head(bf, ptree, &is_head, &comm_f, &nproc);
  MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));
  if (nproc > 1) MPI_Comm_rank(comm, &myid);
  bf_device_setup(comm, nproc);

  bf::gpu::BfResident& r = bf::gpu::bf_store_add(bf, print);
  int uploaded = bf::gpu::upload(lv, up, shp, print.level_butterfly, print.level_half, print.M, print.N,
                                 print.M_loc, print.N_loc, r.dev) ? 1 : 0;
  // every rank declines together, or those that went on would wait in the
  // plan's collectives for those that did not
  if (nproc > 1) MPI_Allreduce(MPI_IN_PLACE, &uploaded, 1, MPI_INT, MPI_MIN, comm);
  if (!uploaded) {
    bf::gpu::bf_store_release(bf);     // no room: the caller stays on the host
    return nullptr;
  }

  auto ranges = [&](char which, int levels) {
    bf::gpu::SweepRanges g;
    g.resize(levels);
    if (which == 'R' || which == 'C') {
      c_bf_get_vec_ranges(bf, ptree, which, g.idx_r.data(), g.inc_r.data(), g.nr.data(),
                          g.idx_c.data(), g.inc_c.data(), g.nc.data());
    } else {
      c_bf_get_mvp_vec_ranges(bf, ptree, which, g.idx_r.data(), g.inc_r.data(), g.nr.data(),
                              g.idx_c.data(), g.inc_c.data(), g.nc.data());
    }
    return g;
  };
  const int L = print.level_butterfly;
  bf::gpu::set_ranges(r.dev, ranges('R', L + 2), ranges('C', L + 2), ranges('N', L + 3),
                      ranges('T', L + 3), nproc, myid);

  r.bridge.comm = comm;
  r.bridge.nproc = nproc;
  r.bridge.exchange_dest = [bf, ptree](int level, char mode, char collect,
                                       const bf::gpu::VecLevel& src, int* dest) {
    c_bf_get_exchange_dest(bf, ptree, level, mode, collect, src.nr, src.nc, src.idx_r, src.inc_r,
                           src.idx_c, src.inc_c, dest);
  };
  r.bridge.all2all_dest = [bf, ptree](int level, char mode_new, const bf::gpu::VecLevel& src,
                                      int* dest) {
    c_bf_get_all2all_dest(bf, ptree, level, mode_new, src.nr, src.nc, src.idx_r, src.inc_r,
                          src.idx_c, src.inc_c, dest);
  };
  r.bridge.all2all_range = [bf, ptree](int level, char mode_new) {
    bf::gpu::Range2 g;
    c_bf_get_all2all_range(bf, ptree, level, mode_new, &g.idx_r, &g.inc_r, &g.nr, &g.idx_c,
                           &g.inc_c, &g.nc);
    return g;
  };
  r.bridge.block_group = [bf, ptree](int level, int index_i, int index_j, char mode) {
    bf::gpu::BfBridge::Group g;
    int np = 1, cf = 0, head = 0;
    c_bf_get_block_group(bf, ptree, level, index_i, index_j, mode, &np, &cf, &head);
    g.nproc = np;
    g.comm = MPI_Comm_f2c(static_cast<MPI_Fint>(cf));
    g.is_head = head != 0;
    return g;
  };
  return &r;
}

}  // namespace
#endif

#ifdef BPACK_BF_GPU_ENABLED
namespace {

// The batched ID of a level's panels, d_entries on the device laid out as the
// host's buffer at host_base (panel b at blk[b]).  False when the device has
// no room for the batch.
bool level_id_on_device(const double* d_entries, int count, const int* rows, const int* cols, void* const* blk,
                        int transposed, const double* host_base, int64_t nnz, double rtol, double atol, int* rank,
                        int* jpvt, double* host_x) {
  bf::gpu::IdBatch batch;
  batch.panels.resize(static_cast<size_t>(count));
  for (int b = 0; b < count; ++b) {
    bf::gpu::IdPanel& p = batch.panels[static_cast<size_t>(b)];
    p.rows = rows[b];
    p.cols = cols[b];
    p.transposed = transposed;
    p.entry_offset = static_cast<const double*>(blk[b]) - host_base;
    if (p.entry_offset < 0 || p.entry_offset + static_cast<int64_t>(p.rows) * p.cols > nnz) {
      gpu_fail("c_bf_gpu_level_id", "a panel lies outside the level's entry buffer");
    }
  }
  if (!bf::gpu::plan_id_batch(batch)) return false;
  bf::gpu::run_id_batch(d_entries, batch, rtol, atol, nullptr);
  fmm::gpu::check_cuda(cudaMemcpy(rank, batch.d_rank, static_cast<size_t>(count) * sizeof(int),
                                  cudaMemcpyDeviceToHost), "c_bf_gpu_level_id rank");
  fmm::gpu::check_cuda(cudaMemcpy(jpvt, batch.d_jpvt, static_cast<size_t>(batch.jpvt_elements) * sizeof(int),
                                  cudaMemcpyDeviceToHost), "c_bf_gpu_level_id jpvt");
  fmm::gpu::check_cuda(cudaMemcpy(host_x, batch.d_x, static_cast<size_t>(batch.mat_elements) * sizeof(double),
                                  cudaMemcpyDeviceToHost), "c_bf_gpu_level_id X");
  // the Fortran's pivots are one based
  for (int64_t t = 0; t < batch.jpvt_elements; ++t) jpvt[t] += 1;
  bf::gpu::release_id_batch(batch);
  return true;
}

}  // namespace
#endif

// The rank-reveal of one level of a butterfly construction, on the device
// (GPU_BACKEND/bf_gpu).  BF_compress_NlogN does it per block under OpenMP -
// a truncated QRCP, Q^H applied to the panel, a triangular solve - and the
// level's panels already sit in one contiguous buffer, because
// element_Zmn_blocklist_user filled alldat_loc_in for the whole level.  So
// the level goes to the device as one copy and comes back as one batch.
//
// `entries` is the level's buffer and `blk` each panel's pointer into it; the
// offsets are taken by subtraction here rather than recomputed in the Fortran,
// because the fill compacts masked rows and skips empty blocks, so its block
// order is not a reliable cumulative sum.  `transposed` is 1 for the 'C' half,
// which factorizes the transpose.
//
// The results come back in the layout the batch lays out: `rank` one per
// block, `jpvt` the cumulative sum of each panel's column count, `xout` the
// cumulative sum of m * n, both in block order.  info 0 on success, 1 when the
// device had no room, so the caller runs the CPU path instead.
void c_bf_gpu_level_id(const int* nblk, const int* rows, const int* cols, void* const* blk,
                       const int* transposed, const void* entries, const int64_t* nnz,
                       const double* rtol, const double* atol, int* rank, int* jpvt, void* xout,
                       int* info) {
#ifdef BPACK_BF_GPU_ENABLED
  *info = 0;
  const int count = *nblk;
  if (count <= 0) return;
  // real double only, which the guard above already ensures
  const double* host_entries = static_cast<const double*>(entries);
  double* host_x = static_cast<double*>(xout);
  try {
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    double* d_entries = reinterpret_cast<double*>(
        heap.try_alloc(static_cast<size_t>(*nnz) * sizeof(double), fmm::gpu::ResidentKind::Adaptive));
    if (d_entries == nullptr) { *info = 1; return; }
    fmm::gpu::check_cuda(cudaMemcpy(d_entries, host_entries, static_cast<size_t>(*nnz) * sizeof(double),
                                    cudaMemcpyHostToDevice), "c_bf_gpu_level_id entries");
    if (!level_id_on_device(d_entries, count, rows, cols, blk, *transposed, host_entries, *nnz, *rtol, *atol, rank,
                            jpvt, host_x)) {
      *info = 1;
    }
    heap.free(reinterpret_cast<char*>(d_entries));
  } catch (const std::exception& e) {
    gpu_fail("c_bf_gpu_level_id", e.what());
  }
#else
  (void)nblk; (void)rows; (void)cols; (void)blk; (void)transposed; (void)entries; (void)nnz;
  (void)rtol; (void)atol; (void)rank; (void)jpvt; (void)xout;
  *info = 1;   // no butterfly backend in this precision: the caller stays on the host
#endif
}

// c_bf_gpu_level_id on the level's entries the GPU evaluator `handle` left on
// the device (c_bf_gpu_eval_blocks with keep): no copy of the entries either
// way.  `entries` is only the host buffer they mirror, for the panels' offsets.
// info 1: no room for the batch - the entries are then copied into `entries`,
// so the host path that follows finds them there.
void c_bf_gpu_level_id_kept(void* handle, const int* nblk, const int* rows, const int* cols, void* const* blk,
                            const int* transposed, void* entries, const int64_t* nnz, const double* rtol,
                            const double* atol, int* rank, int* jpvt, void* xout, int* info) {
#ifdef BPACK_BF_GPU_ENABLED
  *info = 0;
  try {
    BfEntries& e = *static_cast<BfEntries*>(handle);
    if (e.kept == nullptr || e.kept_elements != *nnz) gpu_fail("c_bf_gpu_level_id_kept", "no entries kept on the device");
    if (*nblk > 0 && !level_id_on_device(e.kept, *nblk, rows, cols, blk, *transposed,
                                         static_cast<const double*>(entries), *nnz, *rtol, *atol, rank, jpvt,
                                         static_cast<double*>(xout))) {
      fmm::gpu::check_cuda(cudaMemcpy(entries, e.kept, static_cast<size_t>(*nnz) * sizeof(double),
                                      cudaMemcpyDeviceToHost), "c_bf_gpu_level_id_kept entries");
      *info = 1;
    }
    e.drop_kept();
  } catch (const std::exception& ex) {
    gpu_fail("c_bf_gpu_level_id_kept", ex.what());
  }
#else
  (void)handle; (void)nblk; (void)rows; (void)cols; (void)blk; (void)transposed; (void)entries; (void)nnz;
  (void)rtol; (void)atol; (void)rank; (void)jpvt; (void)xout;
  *info = 1;
#endif
}

// A butterfly matvec on device-resident vectors (GPU_BACKEND/bf_gpu):
// y = alpha op(bf) x + beta y, with x and y already in device memory.  This
// is the device-pointer overload STRUMPACK's Schur sampling needs: it applies
// F12, then F11's inverse, then F21, and with host intermediates every stage
// pays a transfer - which is what made enabling the GPU for an H2 front a
// regression before this existed.
//
// The factors and the sweep plans stay resident between calls (bf_store.hpp),
// so the tens of thousands of samples per front upload nothing and build no
// plan.  info 1 means the device declined - no room for the factors, or a
// case the sweep does not handle yet, such as the closing factor's sub-group
// spanning several ranks - and the caller runs its host path instead.
void c_bf_gpu_mult(void* bf, void* ptree, const char* op, const int* nrhs, const void* x,
                   const int* ldx, void* y, const int* ldy, const double* alpha,
                   const double* beta, double* flops, int* info) {
#ifdef BPACK_BF_GPU_ENABLED
  *info = 1;
  if (flops != nullptr) *flops = 0.0;
  try {
    fmm::gpu::HostRange range("bf matvec, device vectors");
    bf::gpu::BfResident* r = bf_resident(bf, ptree);
    if (r == nullptr) return;
    const char o = (*op == 'N' || *op == 'n') ? 'N' : 'T';
    bf::gpu::SweepPlan& plan = r->plan(o);
    if (!plan.built) {
      // two small collectives per level, once per block and orientation
      bf::gpu::build_sweep_plan(r->bridge, r->dev, o, *nrhs, plan);
    }
    // everything on the backend's stream, and the product complete on return
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    const double f = bf::gpu::sweep(r->dev, plan, *nrhs, static_cast<const double*>(x), *ldx,
                                    static_cast<double*>(y), *ldy, *alpha, *beta, r->scratch, stream);
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "c_bf_gpu_mult");
    if (flops != nullptr) *flops = f;
    *info = 0;
  } catch (const std::exception&) {
    // a case the sweep refuses: leave info 1 and let the caller use the host
    *info = 1;
  }
#else
  (void)bf; (void)ptree; (void)op; (void)nrhs; (void)x; (void)ldx; (void)y; (void)ldy;
  (void)alpha; (void)beta;
  if (flops != nullptr) *flops = 0.0;
  *info = 1;
#endif
}

// The same product on the caller's host arrays (C_BF_Mult): x goes to the
// device, the resident sweep runs, y = op(bf) x comes back.  For a caller
// that keeps its vectors on the host; one that can keep them on the device
// saves both copies with c_bf_gpu_mult.  The device copies of x and y are
// taken per call from the arena, which is a bump of a pointer.  info 1: the
// device declined, and the caller runs its host path.
void c_bf_gpu_mult_host(void* bf, void* ptree, const char* op, const int* nrhs, const void* x,
                        const int* ldx, void* y, const int* ldy, double* flops, int* info) {
#ifdef BPACK_BF_GPU_ENABLED
  *info = 1;
  *flops = 0.0;
  try {
    fmm::gpu::HostRange range("bf matvec, host vectors");
    bf::gpu::BfResident* r = bf_resident(bf, ptree);
    if (r == nullptr) return;
    const char o = (*op == 'N' || *op == 'n') ? 'N' : 'T';
    const int nin = o == 'N' ? r->dev.N_loc : r->dev.M_loc;
    const int nout = o == 'N' ? r->dev.M_loc : r->dev.N_loc;
    bf::gpu::SweepPlan& plan = r->plan(o);
    if (!plan.built) bf::gpu::build_sweep_plan(r->bridge, r->dev, o, *nrhs, plan);

    auto& heap = fmm::gpu::DeviceHeap::instance();
    const size_t bin = static_cast<size_t>(std::max(nin, 1)) * *nrhs * sizeof(double);
    const size_t bout = static_cast<size_t>(std::max(nout, 1)) * *nrhs * sizeof(double);
    char* dx = heap.try_alloc(bin, fmm::gpu::ResidentKind::Matvec);
    char* dy = heap.try_alloc(bout, fmm::gpu::ResidentKind::Matvec);
    // a rank without room must not leave the others in the sweep's exchanges
    int room = (dx != nullptr && dy != nullptr) ? 1 : 0;
    if (r->bridge.nproc > 1) MPI_Allreduce(MPI_IN_PLACE, &room, 1, MPI_INT, MPI_MIN, r->bridge.comm);
    if (!room) {
      if (dx != nullptr) heap.free(dx);
      if (dy != nullptr) heap.free(dy);
      return;
    }
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    if (nin > 0) {
      fmm::gpu::check_cuda(cudaMemcpy2DAsync(dx, static_cast<size_t>(nin) * sizeof(double), x,
                                             static_cast<size_t>(*ldx) * sizeof(double),
                                             static_cast<size_t>(nin) * sizeof(double), *nrhs,
                                             cudaMemcpyHostToDevice, stream),
                           "c_bf_gpu_mult_host x");
    }
    const double f = bf::gpu::sweep(r->dev, plan, *nrhs, reinterpret_cast<const double*>(dx), std::max(nin, 1),
                                    reinterpret_cast<double*>(dy), std::max(nout, 1), 1.0, 0.0, r->scratch,
                                    stream);
    if (nout > 0) {
      fmm::gpu::check_cuda(cudaMemcpy2DAsync(y, static_cast<size_t>(*ldy) * sizeof(double), dy,
                                             static_cast<size_t>(nout) * sizeof(double),
                                             static_cast<size_t>(nout) * sizeof(double), *nrhs,
                                             cudaMemcpyDeviceToHost, stream),
                           "c_bf_gpu_mult_host y");
    }
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "c_bf_gpu_mult_host");
    heap.free(dx);
    heap.free(dy);
    *flops = f;
    *info = 0;
  } catch (const std::exception&) {
    *info = 1;
  }
#else
  (void)bf; (void)ptree; (void)op; (void)nrhs; (void)x; (void)ldx; (void)y; (void)ldy;
  *flops = 0.0;
  *info = 1;
#endif
}

// The entries of a list of intersections of a butterfly, on the device
// (bf_extract.hpp), for C_BF_ExtractElement: the intersections as it receives
// them (rows 1..M and columns 1..N in the butterfly's order, stacked; pgidx,
// 0 based, and pmaps the process grid of each), every block into alldat_loc
// in order, column major.  A butterfly on one process whose every block this
// rank owns in a 1 x 1 grid; otherwise info 1, and the host path runs.  Not
// collective, except for the device's setup at a butterfly's first use, which
// on one process involves no other rank.
void c_bf_gpu_extract(void* bf, void* ptree, const int* ninter, const int* allrows, const int* allcols,
                      const int* rowidx, const int* colidx, const int* pgidx, const int* npmap,
                      const int* pmaps, void* alldat_loc, const int64_t* nalldat_loc, double* flops,
                      int* info) {
#ifdef BPACK_BF_GPU_ENABLED
  *info = 1;
  *flops = 0.0;
  try {
    const int count = *ninter;
    int is_head = 0, comm_f = 0, nproc = 1;
    c_bf_get_sub_head(bf, ptree, &is_head, &comm_f, &nproc);
    if (nproc != 1) return;
    int64_t want = 0, nr = 0, nc = 0;
    for (int k = 0; k < count; ++k) {
      const int pg = pgidx[k];
      if (pmaps[pg] != 1 || pmaps[pg + *npmap] != 1 || pmaps[pg + 2 * *npmap] != 0) return;
      want += static_cast<int64_t>(rowidx[k]) * colidx[k];
      nr += rowidx[k];
      nc += colidx[k];
    }
    if (want != *nalldat_loc) return;
    bf::gpu::BfResident* r = bf_resident(bf, ptree);
    if (r == nullptr) return;
    std::vector<int> rows(allrows, allrows + nr), cols(allcols, allcols + nc);
    for (int& i : rows) i -= 1;
    for (int& j : cols) j -= 1;
    bf::gpu::ExtractProfile prof;
    if (!bf::gpu::extract(r->dev, r->leaves, count, rowidx, colidx, rows.data(), cols.data(),
                          static_cast<double*>(alldat_loc), prof, r->extract_meta, r->extract_meta_device)) {
      return;  // no room on the device for the partial products
    }
    *flops = prof.flops;
    *info = 0;
    if (fmm::env::trace("phase")) {
      std::printf("BF GPU extract: %d intersections, %lld entries, %lld GEMMs, %.3e flops: plan %.4f s, "
                  "run %.4f s\n",
                  count, static_cast<long long>(prof.entries), static_cast<long long>(prof.items), prof.flops,
                  prof.plan, prof.run);
      std::fflush(stdout);
    }
  } catch (const std::exception& e) {
    gpu_fail("c_bf_gpu_extract", e.what());
  }
#else
  (void)bf; (void)ptree; (void)ninter; (void)allrows; (void)allcols; (void)rowidx; (void)colidx;
  (void)pgidx; (void)npmap; (void)pmaps; (void)alldat_loc; (void)nalldat_loc;
  *flops = 0.0;
  *info = 1;
#endif
}

// The registration of a butterfly's GPU evaluator (C_BPACK_wrapper.cpp's
// c_bf_set_gpu_*), before its construction: `evaluator` points to a
// std::shared_ptr<fmm::gpu::Evaluator>, an empty one clearing it.
void c_bf_gpu_set_evaluator(void* bf, const void* evaluator) {
#ifdef BPACK_BF_GPU_ENABLED
  if (bf == nullptr || evaluator == nullptr) gpu_fail("c_bf_gpu_set_evaluator", "null argument");
  const auto& ev = *static_cast<const std::shared_ptr<fmm::gpu::Evaluator>*>(evaluator);
  if (!ev) {
    bf_entries().erase(bf);
    return;
  }
  auto& slot = bf_entries()[bf];
  if (!slot) slot = std::make_unique<BfEntries>();
  slot->evaluator = ev;
#else
  (void)bf;
  (void)evaluator;
#endif
}

// For C_BF_Construct_Element_Compute: the handle of the butterfly's GPU
// evaluator, null when it has none, with its points laid out for M rows and
// N columns.  xyz: the butterfly mesh's coordinates as the Fortran holds them,
// xyz(1:dim, -N:M) - row i at i, column j at -j - or null for a mesh without
// geometry.  The evaluator's first-use costs (an NVRTC compile, a first
// launch) are paid here, outside the construction's own work, and returned.
void c_bf_gpu_entries(void* bf, const int* M, const int* N, const int* dim, const double* xyz, void** handle,
                      double* warm_seconds) {
  *handle = nullptr;
  *warm_seconds = 0.0;
#ifdef BPACK_BF_GPU_ENABLED
  try {
    auto it = bf_entries().find(bf);
    if (it == bf_entries().end() || !it->second->evaluator) return;
    BfEntries& e = *it->second;
    if (e.evaluator->needs_coordinates() && xyz == nullptr) {
      gpu_fail("c_bf_gpu_entries",
               "the GPU evaluator reads coordinates, but the butterfly has none (built with nogeo 0 or 4 it has)");
    }
    if (e.M != *M || e.N != *N) {
      e.drop_points();
      auto& heap = fmm::gpu::DeviceHeap::instance();
      heap.ensure_initialized();
      const int64_t np = static_cast<int64_t>(*M) + *N;
      std::vector<int64_t> ids(static_cast<size_t>(np));
      for (int64_t p = 0; p < np; ++p) ids[p] = p < *M ? p : p - *M;
      e.d_ids = heap.alloc<int64_t>(ids.size() * sizeof(int64_t), fmm::gpu::ResidentKind::Matvec);
      fmm::gpu::check_cuda(cudaMemcpy(e.d_ids, ids.data(), ids.size() * sizeof(int64_t), cudaMemcpyHostToDevice),
                           "c_bf_gpu_entries ids");
      if (xyz != nullptr && *dim > 0) {
        const int d = *dim;
        std::vector<double> table(static_cast<size_t>(np) * d);
        for (int64_t r = 0; r < *M; ++r)
          for (int q = 0; q < d; ++q) table[r * d + q] = xyz[(r + 1 + *N) * d + q];
        for (int64_t c = 0; c < *N; ++c)
          for (int q = 0; q < d; ++q) table[(*M + c) * d + q] = xyz[(*N - c - 1) * d + q];
        e.d_xyz = heap.alloc<double>(table.size() * sizeof(double), fmm::gpu::ResidentKind::Matvec);
        fmm::gpu::check_cuda(cudaMemcpy(e.d_xyz, table.data(), table.size() * sizeof(double),
                                        cudaMemcpyHostToDevice),
                             "c_bf_gpu_entries coordinates");
        e.dim = d;
      }
      e.M = *M;
      e.N = *N;
    }
    const int count = std::max(1, std::min(64, std::min(*M, *N)));
    e.warm_seconds = e.evaluator->warm_up<double>(e.points(), count, fmm::gpu::Context::instance().stream());
    *warm_seconds = e.warm_seconds;
    *handle = &e;
  } catch (const std::exception& ex) {
    gpu_fail("c_bf_gpu_entries", ex.what());
  }
#else
  (void)bf; (void)M; (void)N; (void)dim; (void)xyz;
#endif
}

// The entries of nblk blocks by the GPU evaluator `handle` (c_bf_gpu_entries),
// for element_Zmn_blocklist_user: block b's rowidx[b] rows and colidx[b]
// columns, stacked in rows and cols as the butterfly's mesh numbers them
// (rows 1..M, columns M+1..M+N), into out (host), each block column major,
// stacked in order.  keep 1: left on the device instead, for the level's
// rank-reveal (c_bf_gpu_level_id_kept), and out is not written.  info 1: no
// room on the device; the caller computes them on the host instead.
void c_bf_gpu_eval_blocks(void* handle, const int* nblk, const int* rows, const int* cols, const int* rowidx,
                          const int* colidx, void* out, const int64_t* nout, const int* keep, int* info) {
  *info = 1;
#ifdef BPACK_BF_GPU_ENABLED
  try {
    BfEntries& e = *static_cast<BfEntries*>(handle);
    const int count = *nblk;
    std::vector<int64_t> r_at(static_cast<size_t>(count) + 1, 0), c_at(r_at.size(), 0), o_at(r_at.size(), 0);
    int max_m = 0, max_n = 0;
    for (int k = 0; k < count; ++k) {
      r_at[k + 1] = r_at[k] + rowidx[k];
      c_at[k + 1] = c_at[k] + colidx[k];
      o_at[k + 1] = o_at[k] + static_cast<int64_t>(rowidx[k]) * colidx[k];
      max_m = std::max(max_m, rowidx[k]);
      max_n = std::max(max_n, colidx[k]);
    }
    if (o_at[count] != *nout) gpu_fail("c_bf_gpu_eval_blocks", "the blocks do not fill the output");
    if (o_at[count] == 0) {
      *info = 0;
      return;
    }
    const int64_t nr = r_at[count], nc = c_at[count];
    auto& heap = fmm::gpu::DeviceHeap::instance();
    char* d_out = heap.try_alloc(static_cast<size_t>(o_at[count]) * sizeof(double), fmm::gpu::ResidentKind::Matvec);
    if (d_out == nullptr) return;
    // the points (the mesh's index minus one, rows then columns) and the
    // blocks, in one image; the blocks name their points by offset into it
    e.meta.clear();
    const size_t at_pts = e.meta.alloc(static_cast<size_t>(nr + nc) * sizeof(int));
    int* pts = reinterpret_cast<int*>(e.meta.host(at_pts));
    for (int64_t i = 0; i < nr; ++i) pts[i] = rows[i] - 1;
    for (int64_t j = 0; j < nc; ++j) pts[nr + j] = cols[j] - 1;
    double* dout = reinterpret_cast<double*>(d_out);
    std::vector<fmm::gpu::EvalItemT<double>> items(static_cast<size_t>(count));
    for (int k = 0; k < count; ++k) {
      fmm::gpu::EvalItemT<double>& it = items[k];
      it.out = dout + o_at[k];
      it.ld = std::max(rowidx[k], 1);
      it.m = rowidx[k];
      it.n = colidx[k];
      it.rows = fmm::gpu::IndexList{static_cast<int64_t>(at_pts + r_at[k] * sizeof(int)), 0, nullptr};
      it.cols = fmm::gpu::IndexList{static_cast<int64_t>(at_pts + (nr + c_at[k]) * sizeof(int)), 0, nullptr};
    }
    const size_t at_items = e.meta.append(items);
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    char* dm = e.meta.upload(e.meta_device, stream);
    e.evaluator->eval<double>(items, reinterpret_cast<const fmm::gpu::EvalItemT<double>*>(dm + at_items), max_m,
                              max_n, dm, e.points(), stream);
    if (*keep != 0) {
      // for the rank-reveal that follows: no copy to the host
      fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "c_bf_gpu_eval_blocks");
      e.drop_kept();
      e.kept = dout;
      e.kept_elements = o_at[count];
      *info = 0;
      return;
    }
    fmm::gpu::check_cuda(cudaMemcpyAsync(out, d_out, static_cast<size_t>(o_at[count]) * sizeof(double),
                                         cudaMemcpyDeviceToHost, stream),
                         "c_bf_gpu_eval_blocks entries");
    fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "c_bf_gpu_eval_blocks");
    heap.free(d_out);
    *info = 0;
  } catch (const std::exception& ex) {
    gpu_fail("c_bf_gpu_eval_blocks", ex.what());
  }
#else
  (void)handle; (void)nblk; (void)rows; (void)cols; (void)rowidx; (void)colidx; (void)out; (void)nout;
#endif
}

// The registration of a butterfly's black box on the device (bpack_gpu.h),
// before its construction from matvecs; a null matvec clears it.
void c_bf_set_gpu_matvec(void** bf, bpack_gpu_bf_matvec matvec, void* user) {
#ifdef BPACK_BF_GPU_ENABLED
  if (bf == nullptr || *bf == nullptr) gpu_fail("c_bf_set_gpu_matvec", "null butterfly handle");
  if (matvec == nullptr) {
    bf_matvecs().erase(*bf);
    return;
  }
  bf_matvecs()[*bf] = BfMatvec{matvec, user};
#else
  (void)bf; (void)matvec; (void)user;
#endif
}

// BF_randomized (less_adapt 1) on the device, for BF_randomized_gpu
// (bf_gpu/bf_randomized.hpp): a butterfly on one process with L levels, its
// leaves' ends ms and ns (one based, cumulative, 2^L each).  The black box is
// the device one registered for bf, else `matvec`, the C callback, on pinned
// host copies of the vectors.  On success the shapes of every block - V's,
// Kerl(1..L)'s row by row, U's - are in rows and cols, the blocks wait on the
// device under `handle` for c_bf_gpu_randomized_fetch, and prof holds the
// seconds of all, the black box, the sweeps, the range finding and solves,
// the error test; the black box's calls and columns, the tries of the rank
// guess and the test's error.  info 1: the device declined (a butterfly on
// several processes, no room, a rank guess that never settled).
void c_bf_gpu_randomized(void* bf, void* ptree, const int* L, const int* level_half, const int* M, const int* N,
                         const int* ms, const int* ns, const int* rank0, const double* rankrate, const int* itermax,
                         const int* nbundle, const double* tol_rdetect, const double* tol_ls, void* matvec,
                         void* quant, void** handle, int* rows, int* cols, double* prof, int* info) {
  *info = 1;
  *handle = nullptr;
#ifdef BPACK_BF_GPU_ENABLED
  try {
    int is_head = 0, comm_f = 0, nproc = 1;
    c_bf_get_sub_head(bf, ptree, &is_head, &comm_f, &nproc);
    if (nproc != 1 || *L < 1) return;
    bf_device_setup(MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f)), nproc);

    bf::gpu::RandomizedConfig c;
    c.L = *L;
    c.level_half = *level_half;
    c.M = *M;
    c.N = *N;
    const int leaves = 1 << *L;
    c.row_off.assign(static_cast<size_t>(leaves) + 1, 0);
    c.col_off.assign(static_cast<size_t>(leaves) + 1, 0);
    for (int k = 0; k < leaves; ++k) {
      c.row_off[k + 1] = ms[k];
      c.col_off[k + 1] = ns[k];
    }
    if (c.row_off.back() != c.M || c.col_off.back() != c.N) gpu_fail("c_bf_gpu_randomized", "the leaves do not cover the block");
    c.rank0 = *rank0;
    c.rankrate = *rankrate;
    c.itermax = *itermax;
    c.nbundle = std::max(*nbundle, 1);
    c.tol_rdetect = *tol_rdetect;
    c.tol_ls = *tol_ls;

    auto state = std::make_unique<BfRandomizedRun>();
    BfRandomizedRun* s = state.get();
    cudaStream_t stream = fmm::gpu::Context::instance().stream();
    bf::gpu::BlackBox box;
    auto dev = bf_matvecs().find(bf);
    if (dev != bf_matvecs().end()) {
      const BfMatvec m = dev->second;
      box = [m, c, stream](char op, int nvec, const double* x, double* y) {
        const int nin = static_cast<int>(op == 'N' ? c.N : c.M), nout = static_cast<int>(op == 'N' ? c.M : c.N);
        // x is ready for a black box on a stream of its own too
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "the butterfly's device black box");
        m.fn(op, nin, nout, nvec, x, y, stream, m.user);
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "the butterfly's device black box");
      };
    } else {
      if (matvec == nullptr) return;
      const BfHostMatvec fn = reinterpret_cast<BfHostMatvec>(matvec);
      box = [s, fn, quant, c, stream](char op, int nvec, const double* x, double* y) {
        int nin = static_cast<int>(op == 'N' ? c.N : c.M), nout = static_cast<int>(op == 'N' ? c.M : c.N);
        const size_t bin = static_cast<size_t>(nin) * nvec * sizeof(double);
        const size_t bout = static_cast<size_t>(nout) * nvec * sizeof(double);
        double* hx = static_cast<double*>(s->hx.reserve(std::max<size_t>(bin, 1)));
        double* hy = static_cast<double*>(s->hy.reserve(std::max<size_t>(bout, 1)));
        fmm::gpu::check_cuda(cudaMemcpyAsync(hx, x, bin, cudaMemcpyDeviceToHost, stream), "the butterfly's black box x");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "the butterfly's black box x");
        std::memset(hy, 0, bout);
        const char trans[2] = {op, '\0'};
        double one = 1.0, zero = 0.0;
        fn(trans, &nin, &nout, &nvec, hx, hy, quant, &one, &zero);
        fmm::gpu::check_cuda(cudaMemcpyAsync(y, hy, bout, cudaMemcpyHostToDevice, stream), "the butterfly's black box y");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "the butterfly's black box y");
      };
    }
    s->run = std::make_unique<bf::gpu::Randomized>(c, box);
    if (!s->run->run()) return;
    std::vector<int> r, k;
    s->run->shapes(r, k);
    std::copy(r.begin(), r.end(), rows);
    std::copy(k.begin(), k.end(), cols);
    const bf::gpu::RandomizedProfile& p = s->run->prof;
    prof[0] = p.total;
    prof[1] = p.matvec;
    prof[2] = p.sweep;
    prof[3] = p.range;
    prof[4] = p.test;
    prof[5] = static_cast<double>(p.calls);
    prof[6] = static_cast<double>(p.columns);
    prof[7] = p.tries;
    prof[8] = p.error;
    // the pinned staging goes now; the factors stay for the fetch
    s->hx.release();
    s->hy.release();
    *handle = state.release();
    *info = 0;
  } catch (const std::exception& e) {
    gpu_fail("c_bf_gpu_randomized", e.what());
  }
#else
  (void)bf; (void)ptree; (void)L; (void)level_half; (void)M; (void)N; (void)ms; (void)ns; (void)rank0;
  (void)rankrate; (void)itermax; (void)nbundle; (void)tol_rdetect; (void)tol_ls; (void)matvec; (void)quant;
  (void)rows; (void)cols; (void)prof;
#endif
}

// The blocks of a c_bf_gpu_randomized into ptrs, in its order (null: not
// wanted), and the handle released.  seconds: the copy's time.
void c_bf_gpu_randomized_fetch(void* handle, void* const* ptrs, double* seconds) {
  *seconds = 0.0;
#ifdef BPACK_BF_GPU_ENABLED
  if (handle == nullptr) return;
  std::unique_ptr<BfRandomizedRun> state(static_cast<BfRandomizedRun*>(handle));
  try {
    const double t0 = MPI_Wtime();
    std::vector<int> r, k;
    state->run->shapes(r, k);
    std::vector<double*> dst(r.size());
    for (size_t b = 0; b < dst.size(); ++b) dst[b] = static_cast<double*>(ptrs[b]);
    state->run->download(dst);
    state->run->release();
    *seconds = MPI_Wtime() - t0;
  } catch (const std::exception& e) {
    gpu_fail("c_bf_gpu_randomized_fetch", e.what());
  }
#else
  (void)handle; (void)ptrs;
#endif
}

// Drop what a butterfly holds on the device.  The caller must do this before
// the block is freed: the store is keyed on its address, and although every
// lookup re-checks the block's layout, releasing is what keeps the device
// memory from being held by a block that no longer exists.
void c_bf_gpu_release(void* bf) {
#ifdef BPACK_BF_GPU_ENABLED
  try {
    bf::gpu::bf_store_release(bf);
    bf_entries().erase(bf);  // and its GPU evaluator, if it had one
    bf_matvecs().erase(bf);  // and its device black box
  } catch (const std::exception& e) {
    gpu_fail("c_bf_gpu_release", e.what());
  }
#else
  (void)bf;
#endif
}

// Device memory for a caller that holds vectors across several backend calls
// (STRUMPACK's Schur sampling chain: F12's product, F11's solve, F21's
// product, without a host hop between them).
//
// It comes from the backend's own arena rather than a cudaMalloc of the
// caller's, because the arena takes most of what the device had free when it
// was reserved (BPACK_GPU_HEAP_FRACTION, 0.85 by default): a separate
// allocation would be competing with it for what little is left, and would
// fail unpredictably.  `elements` counts matrix entries, not bytes, so the
// caller does not need to know the precision.  ptr is null when the arena has
// no room, and the caller keeps its host path.
void c_bpack_gpu_workspace_alloc(const int64_t* elements, void** ptr) {
#ifdef BPACK_GPU_ENABLED
  *ptr = nullptr;
  try {
    if (*elements <= 0) return;
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    *ptr = heap.try_alloc(static_cast<size_t>(*elements) * sizeof(Scalar),
                          fmm::gpu::ResidentKind::Matvec);
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_workspace_alloc", e.what());
  }
#else
  (void)elements;
  *ptr = nullptr;
#endif
}

void c_bpack_gpu_workspace_free(void* ptr) {
#ifdef BPACK_GPU_ENABLED
  try {
    if (ptr != nullptr) fmm::gpu::DeviceHeap::instance().free(reinterpret_cast<char*>(ptr));
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_workspace_free", e.what());
  }
#else
  (void)ptr;
#endif
}

// Copy between the caller's host memory and such a workspace.  `elements`
// counts matrix entries; `to_device` selects the direction.
void c_bpack_gpu_workspace_copy(void* device, void* host, const int64_t* elements,
                                const int* to_device) {
#ifdef BPACK_GPU_ENABLED
  try {
    if (device == nullptr || host == nullptr || *elements <= 0) return;
    const size_t bytes = static_cast<size_t>(*elements) * sizeof(Scalar);
    fmm::gpu::check_cuda(*to_device ? cudaMemcpy(device, host, bytes, cudaMemcpyHostToDevice)
                                    : cudaMemcpy(host, device, bytes, cudaMemcpyDeviceToHost),
                         "c_bpack_gpu_workspace_copy");
  } catch (const std::exception& e) {
    gpu_fail("c_bpack_gpu_workspace_copy", e.what());
  }
#else
  (void)device; (void)host; (void)elements; (void)to_device;
#endif
}

}  // extern "C"
