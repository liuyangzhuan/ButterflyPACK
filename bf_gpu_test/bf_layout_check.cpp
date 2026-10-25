// Check bf_gpu's sweep against a real butterfly that ButterflyPACK built, at
// any rank count: upload through the bridge, build the sweep plan from the
// Fortran's own grids, and compare both orientations against
// BF_block_MVP_dat.  This is what catches a wrong index map - a wrong one
// gives plausible numbers, not an obvious failure.
//
// At one process it also checks the structural assumptions a single-process
// sweep could have been written against, and which do NOT hold once the
// blocks are distributed:
//   (1) a level's grid is 2^(l-1) x 2^(L-l+1);
//   (2) idx_r/idx_c are 1 and inc_r/inc_c are 1;
//   (3) every kernel block in a row shares a row count.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <mpi.h>
#include <cstdlib>
#include <random>
#include <vector>

#include "dBPACK_wrapper.h"
#include "bf_gpu/bf_upload.hpp"
#include "bf_gpu/bf_sweep_n.hpp"
#include <cuda_runtime.h>

namespace {
// a smooth kernel: 1/(1 + |x - y|), so the block compresses but still has
// real butterfly levels when LRlevel asks for them
double gx(int i, int n) { return double(i) / double(n); }
void elem(int* m, int* n, double* v, C2Fptr) {
  *v = 1.0 / (1.0 + std::fabs(gx(*m, 4096) - gx(*n, 4096) + 0.5));
}
// The block callback, which BF_compress_NlogN calls through
// element_Zmn_blocklist_user even at elem_extract 0: the intersections are
// stacked, each block column major, and here every one is owned by this rank
// (1 x 1 process groups), so alldat_loc holds them all in order.
void elem_block(int* Ninter, int* /*Nallrows*/, int* /*Nallcols*/,
                int64_t* /*Nalldat_loc*/, int* allrows, int* allcols,
                double* alldat_loc, int* rowidx, int* colidx, int* /*pgidx*/,
                int* /*Npmap*/, int* /*pmaps*/, C2Fptr q) {
  int64_t at = 0;
  int r0 = 0, c0 = 0;
  for (int k = 0; k < *Ninter; k++) {
    const int m = rowidx[k], n = colidx[k];
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        int gi = allrows[r0 + i], gj = allcols[c0 + j];
        elem(&gi, &gj, &alldat_loc[at + i + int64_t(j) * m], q);
      }
    }
    at += int64_t(m) * n;
    r0 += m; c0 += n;
  }
}

void distmn(int*, int*, double*, C2Fptr) {}
void nearfar(int*, int*, int* r, C2Fptr) { *r = 1; }  // everything far
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  // before the device pool: CUDA-aware MPI registers the whole allocation
  // that holds a message buffer, so the messages get their own small arena
  fmm::gpu::reserve_exchange_arena(size, 1);
  std::vector<int> groups(size);
  for (int i = 0; i < size; i++) groups[i] = i;

  F2Cptr ptree, option, stats, msh, kerquant, bmat;
  d_c_bpack_createptree(&size, groups.data(), &fcomm, &ptree);
  d_c_bpack_createoption(&option);
  d_c_bpack_createstats(&stats);
  // without cpp=1 the entry path takes ker%FuncZmn, the Fortran callback a
  // C++ caller never sets, and jumps to a null pointer
  d_c_bpack_set_I_option(&option, "cpp", 1);
  d_c_bpack_set_I_option(&option, "nogeo", 1);
  d_c_bpack_set_I_option(&option, "xyzsort", 0);
  d_c_bpack_set_I_option(&option, "verbosity", -2);
  int nmin_leaf = argc > 1 ? std::atoi(argv[1]) : 64;
  d_c_bpack_set_I_option(&option, "Nmin_leaf", nmin_leaf);
  d_c_bpack_set_D_option(&option, "tol_comp", 1e-4);
  d_c_bpack_set_I_option(&option, "LRlevel", 100);   // force butterfly levels
  d_c_bpack_set_I_option(&option, "elem_extract", 0);
  d_c_bpack_set_I_option(&option, "forwardN15flag", 0);

  // real arrays even with nogeo=1: the Fortran dummies are assumed size, so a
  // null actual argument is undefined even where the values go unused
  int Npo = argc > 2 ? std::atoi(argv[2]) : 1024;
  int Ndim = 1, nlevel = 0, myseg = 0;
  std::vector<int> tree(2, Npo), perms(Npo);
  if (rank == 0) std::printf("config: Npo=%d Nmin_leaf=%d P=%d\n", Npo, nmin_leaf, size);
  std::vector<double> coords(static_cast<size_t>(Npo) * Ndim);
  for (int i = 0; i < Npo; i++) coords[i] = gx(i, Npo);
  std::vector<int> nns(Npo, 0);
  d_c_bpack_construct_init(&Npo, &Ndim, coords.data(), nns.data(), &nlevel,
                           tree.data(), perms.data(), &myseg, &bmat, &option,
                           &stats, &msh, &kerquant, &ptree, &distmn, &nearfar,
                           nullptr);

  F2Cptr bf, option2, stats2, msh1, kerquant2, ptree2;
  d_c_bpack_createptree(&size, groups.data(), &fcomm, &ptree2);
  d_c_bpack_copyoption(&option, &option2);
  d_c_bpack_createstats(&stats2);
  int M = Npo, N = Npo, myrow = 0, mycol = 0;
  d_c_bf_construct_init(&M, &N, &myrow, &mycol, nns.data(), nns.data(), &msh, &msh, &bf,
                        &option2, &stats2, &msh1, &kerquant2, &ptree2,
                        &distmn, &nearfar, nullptr);
  d_c_bf_construct_element_compute(&bf, &option2, &stats2, &msh1, &kerquant2,
                                   &ptree2, &elem, &elem_block, nullptr);

  // ---- the bridge
  int L = 0, LH = 0, bM = 0, bN = 0, bMl = 0, bNl = 0, nlev = 0, nblk = 0;
  d_c_bf_get_layout(bf, &L, &LH, &bM, &bN, &bMl, &bNl, &nlev, &nblk);
  std::printf("butterfly: L=%d level_half=%d M=%d N=%d M_loc=%d N_loc=%d "
              "levels=%d blocks=%d\n", L, LH, bM, bN, bMl, bNl, nlev, nblk);
  std::vector<int> nr(nlev), nc(nlev), ir(nlev), ic(nlev), er(nlev), ec(nlev);
  std::vector<int> brows(nblk), bcols(nblk);
  std::vector<void*> bptr(nblk);
  d_c_bf_get_blocks(bf, nr.data(), nc.data(), ir.data(), ic.data(), er.data(),
                    ec.data(), brows.data(), bcols.data(), bptr.data());

  int bad = 0;
  int fails_construct = 0;
  for (int l = 0; l < nlev; l++) {
    std::printf("  level %2d: nr=%-5d nc=%-5d idx_r=%-4d idx_c=%-4d inc_r=%-3d inc_c=%-3d\n",
                l, nr[l], nc[l], ir[l], ic[l], er[l], ec[l]);
    if (size == 1 && (er[l] != 1 || ec[l] != 1 || ir[l] != 1 || ic[l] != 1)) {
      std::printf("    (2) FAIL: one process but idx/inc are not 1\n"); bad++;
    }
  }
  // (1)-(3) hold at one process only: with several ranks a rank owns a
  // strided subset of each level, so the grids are neither the full powers of
  // two nor contiguous.
  if (size == 1) {
  if (nr[0] != 1 || nc[0] != (1 << L)) {
    std::printf("    (1) FAIL: V is %dx%d, expected 1x%d\n", nr[0], nc[0], 1 << L);
    bad++;
  }
  for (int l = 1; l <= L; l++) {
    if (nr[l] != (1 << l) || nc[l] != (1 << (L - l + 1))) {
      std::printf("    (1) FAIL: Kerl(%d) is %dx%d, expected %dx%d\n",
                  l, nr[l], nc[l], 1 << l, 1 << (L - l + 1));
      bad++;
    }
  }
  if (nr[L + 1] != (1 << L) || nc[L + 1] != 1) {
    std::printf("    (1) FAIL: U is %dx%d, expected %dx1\n", nr[L+1], nc[L+1], 1 << L);
    bad++;
  }
  // (3) do all kernel blocks in a row share a row count?
  int at = nr[0] * nc[0];
  for (int l = 1; l <= L; l++) {
    for (int i = 0; i < nr[l]; i++) {
      int seen = -1;
      for (int j = 0; j < nc[l]; j++) {
        const int r = brows[at + i * nc[l] + j];
        if (r == 0) continue;
        if (seen < 0) seen = r;
        else if (seen != r) {
          std::printf("    (3) FAIL: Kerl(%d) row %d has rows %d and %d\n", l, i, seen, r);
          bad++;
        }
      }
    }
    at += nr[l] * nc[l];
  }
  }

  if (bad) {
    std::printf("LAYOUT CHECKS FAILED\n");
    MPI_Finalize();
    return 1;
  }
  std::printf("LAYOUT CHECKS PASS\n");

  // ---- a second construction of the same block with the device
  // ---- rank-reveal, compared against this one.  The two are independent
  // ---- compressions of one operator and may pick different pivots, so they
  // ---- agree to the compression tolerance rather than bitwise; what would
  // ---- show a broken rank-reveal is a difference far above it.
  if (argc > 3 && std::atoi(argv[3]) != 0) {
    F2Cptr bfg, optiong, statsg, mshg, kerquantg, ptreeg;
    d_c_bpack_createptree(&size, groups.data(), &fcomm, &ptreeg);
    d_c_bpack_copyoption(&option, &optiong);
    d_c_bpack_createstats(&statsg);
    // argv[3] == 2 is the control: build the second one on the host too, so
    // that an exactly-zero difference confirms the device path really ran when
    // it is 1 (a gate that silently declined would also "agree")
    const int use_gpu = std::atoi(argv[3]) == 2 ? 0 : 1;
    d_c_bpack_set_I_option(&optiong, "HODLR_use_gpu", use_gpu);
    int Mg = Npo, Ng = Npo, myrowg = 0, mycolg = 0;
    d_c_bf_construct_init(&Mg, &Ng, &myrowg, &mycolg, nns.data(), nns.data(), &msh, &msh, &bfg,
                          &optiong, &statsg, &mshg, &kerquantg, &ptreeg, &distmn, &nearfar, nullptr);
    d_c_bf_construct_element_compute(&bfg, &optiong, &statsg, &mshg, &kerquantg, &ptreeg, &elem,
                                     &elem_block, nullptr);
    int Lg = 0, LHg = 0, bMg = 0, bNg = 0, bMlg = 0, bNlg = 0, nlevg = 0, nblkg = 0;
    d_c_bf_get_layout(bfg, &Lg, &LHg, &bMg, &bNg, &bMlg, &bNlg, &nlevg, &nblkg);
    const int nrhs_c = 2;
    std::vector<double> xc(size_t(bNl) * nrhs_c), y0(size_t(bMl) * nrhs_c, 0.0),
                        y1(size_t(bMl) * nrhs_c, 0.0);
    std::mt19937_64 crng(4242);
    std::uniform_real_distribution<double> cuni(-1.0, 1.0);
    for (auto& v : xc) v = cuni(crng);
    int nin_c = bNl, nout_c = bMl, nv_c = nrhs_c;
    char opn = 'N';
    d_c_bf_mult(&opn, xc.data(), y0.data(), &nin_c, &nout_c, &nv_c, &bf, &option2, &stats2, &ptree2);
    d_c_bf_mult(&opn, xc.data(), y1.data(), &nin_c, &nout_c, &nv_c, &bfg, &optiong, &statsg, &ptreeg);
    double num = 0, den = 0;
    for (size_t i = 0; i < y0.size(); i++) { num += (y1[i]-y0[i])*(y1[i]-y0[i]); den += y0[i]*y0[i]; }
    const double rel = den > 0 ? std::sqrt(num/den) : std::sqrt(num);
    const double tol = 1e-4;   // tol_comp set above
    std::printf("  [rank %d] construction: host L=%d blocks=%d, second(use_gpu=%d) L=%d blocks=%d, "
                "||y2 - y_host||/||y_host|| = %.3e  %s\n", rank, L, nblk, use_gpu, Lg, nblkg, rel,
                rel < 50.0 * tol ? "PASS" : "FAIL");
    std::fflush(stdout);
    if (!(rel < 50.0 * tol)) fails_construct++;
  }

  // ---- upload through the bridge, then compare the device sweep against
  // ---- ButterflyPACK's own matvec for both orientations
  std::vector<bf::gpu::LevelDesc> lv(nlev);
  std::vector<std::vector<const double*>> up(nlev);
  std::vector<std::vector<std::pair<int,int>>> shp(nlev);
  {
    int at2 = 0;
    for (int l = 0; l < nlev; l++) {
      lv[l].nr = nr[l]; lv[l].nc = nc[l];
      lv[l].idx_r = ir[l]; lv[l].idx_c = ic[l];
      lv[l].inc_r = er[l]; lv[l].inc_c = ec[l];
      const int want = nr[l] * nc[l];
      for (int b = 0; b < want; b++) {
        up[l].push_back(reinterpret_cast<const double*>(bptr[at2 + b]));
        shp[l].push_back({brows[at2 + b], bcols[at2 + b]});
      }
      at2 += want;
    }
  }
  bf::gpu::BfDevice dev;
  if (!bf::gpu::upload(lv, up, shp, L, LH, bM, bN, bMl, bNl, dev)) {
    std::printf("FAIL: upload refused\n"); MPI_Finalize(); return 1;
  }
  std::printf("uploaded %lld elements over %zu blocks\n",
              (long long)dev.elements, dev.blocks.size());

  // The grids the sweep runs over, and the rules it needs, from the Fortran.
  // GPU_BACKEND is compiled once while these entry points carry a precision
  // prefix, so bf_gpu takes them as closures rather than calling them.
  auto kernel_ranges = [&](char dir) {
    bf::gpu::SweepRanges r;
    r.resize(L + 2);
    d_c_bf_get_vec_ranges(bf, ptree2, dir, r.idx_r.data(), r.inc_r.data(), r.nr.data(),
                          r.idx_c.data(), r.inc_c.data(), r.nc.data());
    return r;
  };
  auto mvp_ranges = [&](char chara) {
    bf::gpu::SweepRanges r;
    r.resize(L + 3);
    d_c_bf_get_mvp_vec_ranges(bf, ptree2, chara, r.idx_r.data(), r.inc_r.data(), r.nr.data(),
                              r.idx_c.data(), r.inc_c.data(), r.nc.data());
    return r;
  };
  int is_head = 0, comm_f = 0, bnproc = 1;
  d_c_bf_get_sub_head(bf, ptree2, &is_head, &comm_f, &bnproc);
  bf::gpu::set_ranges(dev, kernel_ranges('R'), kernel_ranges('C'), mvp_ranges('N'), mvp_ranges('T'),
                      bnproc, rank);

  bf::gpu::BfBridge bridge;
  bridge.comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));
  bridge.nproc = bnproc;
  bridge.exchange_dest = [&](int level, char mode, char collect, const bf::gpu::VecLevel& src, int* dest) {
    d_c_bf_get_exchange_dest(bf, ptree2, level, mode, collect, src.nr, src.nc, src.idx_r, src.inc_r,
                             src.idx_c, src.inc_c, dest);
  };
  bridge.all2all_dest = [&](int level, char mode_new, const bf::gpu::VecLevel& src, int* dest) {
    d_c_bf_get_all2all_dest(bf, ptree2, level, mode_new, src.nr, src.nc, src.idx_r, src.inc_r,
                            src.idx_c, src.inc_c, dest);
  };
  bridge.block_group = [&](int level, int index_i, int index_j, char mode) {
    bf::gpu::BfBridge::Group g;
    int np = 1, cf = 0, head = 0;
    d_c_bf_get_block_group(bf, ptree2, level, index_i, index_j, mode, &np, &cf, &head);
    g.nproc = np;
    g.comm = MPI_Comm_f2c(static_cast<MPI_Fint>(cf));
    g.is_head = head != 0;
    return g;
  };
  bridge.all2all_range = [&](int level, char mode_new) {
    bf::gpu::Range2 r;
    d_c_bf_get_all2all_range(bf, ptree2, level, mode_new, &r.idx_r, &r.inc_r, &r.nr, &r.idx_c,
                             &r.inc_c, &r.nc);
    return r;
  };
  if (rank == 0) {
    for (int t = 0; t < 2; t++) {
      const bf::gpu::SweepRanges& v = t == 0 ? dev.vec_n : dev.vec_t;
      std::printf("  vec grid '%c':", t == 0 ? 'N' : 'T');
      for (int l = 0; l <= L + 2; l++)
        std::printf(" [%d](%dx%d @%d,%d +%d,%d)", l, v.nr[l], v.nc[l], v.idx_r[l],
                    v.idx_c[l], v.inc_r[l], v.inc_c[l]);
      std::printf("\n");
    }
  }

  // The plan is built once per orientation and reused, which is the point of
  // caching it: a matvec runs no collective of its own.  So each orientation
  // is checked several times through one plan, including at a different nrhs,
  // where only the element offsets and message sizes change.
  const int nrhs_cases[] = {3, 7, 3};
  int fails = 0;
  for (int t = 0; t < 2; t++) {
    const char op = t == 0 ? 'N' : 'T';
    const int nin  = (op == 'N') ? bNl : bMl;
    const int nout = (op == 'N') ? bMl : bNl;
    bf::gpu::SweepScratch scratch;
    bf::gpu::SweepPlan plan;
    bf::gpu::build_sweep_plan(bridge, dev, op, nrhs_cases[0], plan);
    const bf::gpu::SweepPlan::Traffic tr = plan.traffic();
    for (int c = 0; c < 3; c++) {
      const int nrhs = nrhs_cases[c];
      std::vector<double> x(size_t(nin) * nrhs), yref(size_t(nout) * nrhs, 0.0),
                          ygpu(size_t(nout) * nrhs, 0.0);
      std::mt19937_64 rng(777 + t * 10 + c);
      std::uniform_real_distribution<double> uni(-1.0, 1.0);
      for (auto& v : x) v = uni(rng);

      // reference: ButterflyPACK's own matvec (alpha 1, beta 0 semantics)
      int nin_c = nin, nout_c = nout, nv = nrhs;
      char opc = op;
      d_c_bf_mult(&opc, x.data(), yref.data(), &nin_c, &nout_c, &nv,
                  &bf, &option2, &stats2, &ptree2);

      auto& heap = fmm::gpu::DeviceHeap::instance();
      heap.ensure_initialized();
      double* dx = heap.alloc<double>(x.size() * sizeof(double), fmm::gpu::ResidentKind::Matvec);
      double* dy = heap.alloc<double>(ygpu.size() * sizeof(double), fmm::gpu::ResidentKind::Matvec);
      cudaMemcpy(dx, x.data(), x.size() * sizeof(double), cudaMemcpyHostToDevice);
      cudaMemset(dy, 0, ygpu.size() * sizeof(double));
      // argv[4] != 0: go through the library's device-pointer entry point
      // (c_bf_gpu_mult), which also exercises the resident store, the cached
      // sweep plan and the bridge closures the wrapper builds.  Otherwise call
      // the sweep directly, as before.
      double fl = 0.0;
      if (argc > 4 && std::atoi(argv[4]) != 0) {
        int nrhs_c2 = nrhs, ldx = nin, ldy = nout, gi = 1;
        double al = 1.0, be = 0.0;
        char opc2 = op;
        d_c_bf_gpu_mult(bf, ptree2, &opc2, &nrhs_c2, dx, &ldx, dy, &ldy, &al, &be, &fl, &gi);
        if (gi != 0) {
          std::printf("  [rank %d] op=%c use=%d: c_bf_gpu_mult declined (info=%d)\n", rank, op, c, gi);
          std::fflush(stdout);
          fails++;
          heap.free(dy); heap.free(dx);
          continue;
        }
      } else {
        fl = bf::gpu::sweep(dev, plan, nrhs, dx, nin, dy, nout, 1.0, 0.0, scratch);
      }
      cudaDeviceSynchronize();
      cudaMemcpy(ygpu.data(), dy, ygpu.size() * sizeof(double), cudaMemcpyDeviceToHost);
      heap.free(dy); heap.free(dx);

      double num = 0.0, den = 0.0;
      for (size_t i = 0; i < yref.size(); i++) {
        num += (ygpu[i] - yref[i]) * (ygpu[i] - yref[i]);
        den += yref[i] * yref[i];
      }
      const double rel = den > 0 ? std::sqrt(num / den) : std::sqrt(num);
      if (rank == 0 || !(rel < 1e-12)) {
        std::printf("  [rank %d] op=%c use=%d nrhs=%d flops=%.3e rel=%.3e  %s   "
                    "exchange: %d levels, pair %lld->/%lld<- (%lld rows), mid %lld->/%lld<- (%lld rows)\n",
                    rank, op, c, nrhs, fl, rel, rel < 1e-12 ? "PASS" : "FAIL", tr.levels_exchanging,
                    (long long)tr.pair_send, (long long)tr.pair_recv, (long long)tr.pair_rows,
                    (long long)tr.mid_send, (long long)tr.mid_recv, (long long)tr.mid_rows);
        std::fflush(stdout);
      }
      if (!(rel < 1e-12)) fails++;
    }
    scratch.release();
    plan.release();
  }
  bf::gpu::release(dev);
  int total = 0;
  fails += fails_construct;
  MPI_Allreduce(&fails, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    std::printf("P=%d: %s\n", size, total == 0 ? "MATVEC MATCHES CPU" : "MATVEC DIFFERS");
  MPI_Finalize();
  return total == 0 ? 0 : 1;
}
#if 0
  std::printf("%s\n", bad == 0 ? "LAYOUT CHECKS PASS" : "LAYOUT CHECKS FAILED");
  MPI_Finalize();
  return bad == 0 ? 0 : 1;
}
#endif
