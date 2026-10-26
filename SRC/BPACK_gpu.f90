! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so.

!> @file BPACK_gpu.f90
!> @brief Fortran side of the GPU backend of the formats other than H2: HODLR
!> (format 1 with LRlevel 0, option%HODLR_use_gpu > 0).  The Fortran code keeps
!> the tree and host copies of all blocks and hands them to the device backend
!> (hodlr_gpu/, through SRC/BPACK_gpu_wrapper.cpp).  The CPU routines stay the
!> numerical reference.

#include "ButterflyPACK_config.fi"
module BPACK_GPU
   use BPACK_DEFS
   use MISC_Utilities
   use iso_c_binding
   implicit none

   !> seconds in the merges' TSQR since the last reset: the local QR, up the
   !> tree, the root's SVD, down the tree, the local products (read and reset
   !> by the shared-block construction)
   real(kind=8) :: hodlr_gpu_tsqr_time(5) = 0

   interface
      subroutine c_bpack_gpu_available(available) bind(c, name="c_bpack_gpu_available")
         import :: c_int
         integer(c_int) :: available
      end subroutine c_bpack_gpu_available

      subroutine c_bpack_gpu_create(gpu) bind(c, name="c_bpack_gpu_create")
         import :: c_ptr
         type(c_ptr) :: gpu
      end subroutine c_bpack_gpu_create

      subroutine c_bpack_gpu_delete(gpu) bind(c, name="c_bpack_gpu_delete")
         import :: c_ptr
         type(c_ptr), value :: gpu
      end subroutine c_bpack_gpu_delete

      subroutine c_bpack_hodlr_gpu_reset(gpu, n_loc, maxlevel) bind(c, name="c_bpack_hodlr_gpu_reset")
         import :: c_ptr, c_int, c_int64_t
         type(c_ptr), value :: gpu
         integer(c_int64_t) :: n_loc
         integer(c_int) :: maxlevel
      end subroutine c_bpack_hodlr_gpu_reset

      subroutine c_bpack_hodlr_gpu_add_shared_level(gpu, level, comm) bind(c, name="c_bpack_hodlr_gpu_add_shared_level")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: level
         integer :: comm
      end subroutine c_bpack_hodlr_gpu_add_shared_level

      subroutine c_bpack_gpu_dm_alloc(gpu, rows, cols, host, id) bind(c, name="c_bpack_gpu_dm_alloc")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, host
         integer(c_int) :: rows, cols, id
      end subroutine c_bpack_gpu_dm_alloc

      subroutine c_bpack_gpu_dm_download(gpu, id, host) bind(c, name="c_bpack_gpu_dm_download")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, host
         integer(c_int) :: id
      end subroutine c_bpack_gpu_dm_download

      subroutine c_bpack_gpu_dm_free(gpu, id) bind(c, name="c_bpack_gpu_dm_free")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: id
      end subroutine c_bpack_gpu_dm_free

      subroutine c_bpack_gpu_dm_redistribute(gpu, comm, src, ncols, dst, c0, nsend, s_rank, s_off, s_rows, &
                                             nrecv, r_rank, r_off, r_rows, tag) bind(c, name="c_bpack_gpu_dm_redistribute")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer :: comm
         integer(c_int) :: src, ncols, dst, c0, nsend, nrecv, tag
         integer(c_int) :: s_rank(*), s_off(*), s_rows(*), r_rank(*), r_off(*), r_rows(*)
      end subroutine c_bpack_gpu_dm_redistribute

      subroutine c_bpack_gpu_dm_tsqr(gpu, comm, x, tol, underflow, rn, tsec, xnew) bind(c, name="c_bpack_gpu_dm_tsqr")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu
         integer :: comm
         integer(c_int) :: x, rn, xnew
         real(c_double) :: tol, underflow
         real(c_double) :: tsec(5)
      end subroutine c_bpack_gpu_dm_tsqr

      subroutine c_bpack_gpu_dm_fetch_sw(gpu, sw) bind(c, name="c_bpack_gpu_dm_fetch_sw")
         import :: c_ptr
         type(c_ptr), value :: gpu, sw
      end subroutine c_bpack_gpu_dm_fetch_sw

      subroutine c_bpack_gpu_dm_gemm_nt(gpu, a, m, k, b, n, c, c_row0) bind(c, name="c_bpack_gpu_dm_gemm_nt")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, b
         integer(c_int) :: a, m, k, n, c, c_row0
      end subroutine c_bpack_gpu_dm_gemm_nt

      subroutine c_bpack_gpu_share_device(comm, share) bind(c, name="c_bpack_gpu_share_device")
         import :: c_int
         integer :: comm
         integer(c_int) :: share
      end subroutine c_bpack_gpu_share_device

      subroutine c_bpack_gpu_init_exchange(nproc, share) bind(c, name="c_bpack_gpu_init_exchange")
         import :: c_int
         integer(c_int) :: nproc, share
      end subroutine c_bpack_gpu_init_exchange

      integer(c_int) function c_bpack_env_check(value) bind(c, name="c_bpack_env_check")
         import :: c_int, c_char
         character(kind=c_char), dimension(*) :: value
      end function c_bpack_env_check

      integer(c_int) function c_bpack_env_trace(value) bind(c, name="c_bpack_env_trace")
         import :: c_int, c_char
         character(kind=c_char), dimension(*) :: value
      end function c_bpack_env_trace

      subroutine c_bpack_hodlr_gpu_add_leaf(gpu, row0, m, d, ldd) bind(c, name="c_bpack_hodlr_gpu_add_leaf")
         import :: c_ptr, c_int, c_int64_t
         type(c_ptr), value :: gpu, d
         integer(c_int64_t) :: row0
         integer(c_int) :: m, ldd
      end subroutine c_bpack_hodlr_gpu_add_leaf

      subroutine c_bpack_hodlr_gpu_add_lowrank(gpu, level, row0, m, col0, n, k, u, ldu, v, ldv, sym, dev_v) &
         bind(c, name="c_bpack_hodlr_gpu_add_lowrank")
         import :: c_ptr, c_int, c_int64_t
         type(c_ptr), value :: gpu, u, v, dev_v
         integer(c_int) :: level, m, n, k, ldu, ldv, sym
         integer(c_int64_t) :: row0, col0
      end subroutine c_bpack_hodlr_gpu_add_lowrank

      subroutine c_bpack_hodlr_gpu_commit_forward(gpu, mbytes, mirrored) bind(c, name="c_bpack_hodlr_gpu_commit_forward")
         import :: c_ptr, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: mbytes, mirrored
      end subroutine c_bpack_hodlr_gpu_commit_forward

      subroutine c_bpack_hodlr_gpu_fetch_block(gpu, idx, hu, ldu, hv, ldv, m, n, k) &
         bind(c, name="c_bpack_hodlr_gpu_fetch_block")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, hu, hv
         integer(c_int) :: idx, ldu, ldv, m, n, k
      end subroutine c_bpack_hodlr_gpu_fetch_block

      subroutine c_bpack_hodlr_gpu_gather_rows(gpu, idx, role, n, rows, out, k) &
         bind(c, name="c_bpack_hodlr_gpu_gather_rows")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, out
         integer(c_int) :: idx, role, n, k
         integer(c_int) :: rows(*)
      end subroutine c_bpack_hodlr_gpu_gather_rows

      subroutine c_bpack_hodlr_gpu_block_rank(gpu, idx, k) bind(c, name="c_bpack_hodlr_gpu_block_rank")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: idx, k
      end subroutine c_bpack_hodlr_gpu_block_rank

      subroutine c_bpack_hodlr_gpu_extract_products(gpu, np, idx, role, nbig, nsmall, small_first, big_rows, small, &
                                                    out, flops) bind(c, name="c_bpack_hodlr_gpu_extract_products")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu, small, out
         integer(c_int) :: np
         integer(c_int) :: idx(*), role(*), nbig(*), nsmall(*), small_first(*), big_rows(*)
         real(c_double) :: flops
      end subroutine c_bpack_hodlr_gpu_extract_products

      subroutine c_bpack_hodlr_gpu_block_v_dm(gpu, idx, id) bind(c, name="c_bpack_hodlr_gpu_block_v_dm")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: idx, id
      end subroutine c_bpack_hodlr_gpu_block_v_dm

      subroutine c_bpack_hodlr_gpu_mult(gpu, trans, nrhs, x, y, mode, flops) bind(c, name="c_bpack_hodlr_gpu_mult")
         import :: c_ptr, c_int, c_char, c_double
         type(c_ptr), value :: gpu, x, y
         character(kind=c_char) :: trans
         integer(c_int) :: nrhs, mode
         real(c_double) :: flops
      end subroutine c_bpack_hodlr_gpu_mult

      subroutine c_bpack_hodlr_gpu_factor_sym(gpu, jitter, mode, out) bind(c, name="c_bpack_hodlr_gpu_factor_sym")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: jitter
         integer(c_int) :: mode
         real(c_double) :: out(12)
      end subroutine c_bpack_hodlr_gpu_factor_sym

      subroutine c_bpack_hodlr_gpu_sym_ready(gpu, ready) bind(c, name="c_bpack_hodlr_gpu_sym_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: ready
      end subroutine c_bpack_hodlr_gpu_sym_ready

      subroutine c_bpack_hodlr_gpu_solve_sym(gpu, nrhs, x, y, mode, flops) bind(c, name="c_bpack_hodlr_gpu_solve_sym")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu, x, y
         integer(c_int) :: nrhs, mode
         real(c_double) :: flops
      end subroutine c_bpack_hodlr_gpu_solve_sym

      subroutine c_bpack_hodlr_gpu_factor_unsym(gpu, jitter, mode, out) bind(c, name="c_bpack_hodlr_gpu_factor_unsym")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: jitter
         integer(c_int) :: mode
         real(c_double) :: out(12)
      end subroutine c_bpack_hodlr_gpu_factor_unsym

      subroutine c_bpack_hodlr_gpu_unsym_ready(gpu, ready) bind(c, name="c_bpack_hodlr_gpu_unsym_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: ready
      end subroutine c_bpack_hodlr_gpu_unsym_ready

      subroutine c_bpack_hodlr_gpu_solve_unsym(gpu, trans, nrhs, x, y, mode, flops) &
         bind(c, name="c_bpack_hodlr_gpu_solve_unsym")
         import :: c_ptr, c_int, c_char, c_double
         type(c_ptr), value :: gpu, x, y
         character(kind=c_char) :: trans
         integer(c_int) :: nrhs, mode
         real(c_double) :: flops
      end subroutine c_bpack_hodlr_gpu_solve_unsym

      subroutine c_bpack_hodlr_gpu_download_unsym(gpu, what, level, row0, col0, m, k, out, found) &
         bind(c, name="c_bpack_hodlr_gpu_download_unsym")
         import :: c_ptr, c_int, c_int64_t
         type(c_ptr), value :: gpu, out
         integer(c_int) :: what, level, m, k, found
         integer(c_int64_t) :: row0, col0
      end subroutine c_bpack_hodlr_gpu_download_unsym

      subroutine c_bpack_hodlr_gpu_download_symnode(gpu, level, c0, n0, n1, q0, q1, g0, g1, s, ipiv, k, found, &
         z0, z1) bind(c, name="c_bpack_hodlr_gpu_download_symnode")
         import :: c_ptr, c_int, c_int64_t
         type(c_ptr), value :: gpu, q0, q1, g0, g1, s, ipiv, z0, z1
         integer(c_int) :: level, n0, n1, k, found
         integer(c_int64_t) :: c0
      end subroutine c_bpack_hodlr_gpu_download_symnode
      subroutine c_bpack_hodlr_gpu_warm_up(gpu, construct, evaluate, seconds) bind(c, name="c_bpack_hodlr_gpu_warm_up")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu
         integer(c_int) :: construct, evaluate
         real(c_double) :: seconds(2)
      end subroutine c_bpack_hodlr_gpu_warm_up
      subroutine c_bpack_hodlr_gpu_construct_ready(gpu, ready) bind(c, name="c_bpack_hodlr_gpu_construct_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: ready
      end subroutine c_bpack_hodlr_gpu_construct_ready

      subroutine c_bpack_hodlr_gpu_set_points(gpu, n, dim, xyz, ids) bind(c, name="c_bpack_hodlr_gpu_set_points")
         import :: c_ptr, c_int64_t, c_int, c_double
         type(c_ptr), value :: gpu
         integer(c_int64_t) :: n
         integer(c_int) :: dim
         real(c_double) :: xyz(*)
         integer(c_int64_t) :: ids(*)
      end subroutine c_bpack_hodlr_gpu_set_points

      subroutine c_bpack_hodlr_gpu_set_addend(gpu, n, nu, nv, r, u, v, umap, vmap) bind(c, name="c_bpack_hodlr_gpu_set_addend")
         import :: c_ptr, c_int64_t, c_int
         type(c_ptr), value :: gpu, u, v, umap, vmap
         integer(c_int64_t) :: n
         integer(c_int) :: nu, nv, r
      end subroutine c_bpack_hodlr_gpu_set_addend

      subroutine c_bpack_hodlr_gpu_eval_dense(gpu, scale, count, r0, m, c0, n, out) &
         bind(c, name="c_bpack_hodlr_gpu_eval_dense")
         import :: c_ptr, c_int, c_int64_t, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: scale
         integer(c_int) :: count, m(*), n(*)
         integer(c_int64_t) :: r0(*), c0(*)
         type(c_ptr) :: out(*)
      end subroutine c_bpack_hodlr_gpu_eval_dense

      ! (host_fn, host_ctx: the entries through the host when the matrix has no GPU evaluator,
      ! c_bpack_hodlr_gpu_host_entries; c_null_funptr otherwise)
      subroutine c_bpack_hodlr_gpu_baca_begin(gpu, scale, mode, variant, nb, r0, m, c0, n, r_est, host_fn, host_ctx) &
         bind(c, name="c_bpack_hodlr_gpu_baca_begin")
         import :: c_ptr, c_funptr, c_int, c_int64_t, c_double
         type(c_ptr), value :: gpu, host_ctx
         type(c_funptr), value :: host_fn
         real(c_double) :: scale
         integer(c_int) :: mode, variant, nb, m(*), n(*), r_est(*)
         integer(c_int64_t) :: r0(*), c0(*)
      end subroutine c_bpack_hodlr_gpu_baca_begin

      subroutine c_bpack_hodlr_gpu_baca_set_columns(gpu, b, cols) bind(c, name="c_bpack_hodlr_gpu_baca_set_columns")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: b, cols(*)
      end subroutine c_bpack_hodlr_gpu_baca_set_columns

      subroutine c_bpack_hodlr_gpu_baca_panels(gpu, na, bl, core) bind(c, name="c_bpack_hodlr_gpu_baca_panels")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, core
         integer(c_int) :: na, bl(*)
      end subroutine c_bpack_hodlr_gpu_baca_panels

      subroutine c_bpack_hodlr_gpu_baca_knn_panels(gpu, na, bl, nc, cols, nr, rows, core) &
         bind(c, name="c_bpack_hodlr_gpu_baca_knn_panels")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, core
         integer(c_int) :: na, bl(*), nc(*), cols(*), nr(*), rows(*)
      end subroutine c_bpack_hodlr_gpu_baca_knn_panels

      subroutine c_bpack_hodlr_gpu_baca_append(gpu, knn, na, bl, ru, jpvt, w, rskip, grams) &
         bind(c, name="c_bpack_hodlr_gpu_baca_append")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, w, grams
         integer(c_int) :: knn, na, bl(*), ru(*), jpvt(*), rskip(*)
      end subroutine c_bpack_hodlr_gpu_baca_append

      subroutine c_bpack_hodlr_gpu_baca_recompress(gpu, na, bl, tol, underflow, rn) &
         bind(c, name="c_bpack_hodlr_gpu_baca_recompress")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: gpu
         integer(c_int) :: na
         integer(c_int) :: bl(*), rn(*)
         real(c_double) :: tol, underflow
      end subroutine c_bpack_hodlr_gpu_baca_recompress

      subroutine c_bpack_hodlr_gpu_baca_download(gpu, na, bl, rn, u, v, keep) bind(c, name="c_bpack_hodlr_gpu_baca_download")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: na, keep
         integer(c_int) :: bl(*), rn(*)
         type(c_ptr) :: u(*), v(*)
      end subroutine c_bpack_hodlr_gpu_baca_download

      subroutine c_bpack_hodlr_gpu_baca_to_dm(gpu, na, bl, rn, ids) bind(c, name="c_bpack_hodlr_gpu_baca_to_dm")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: na
         integer(c_int) :: bl(*), rn(*), ids(*)
      end subroutine c_bpack_hodlr_gpu_baca_to_dm

      subroutine c_bpack_gpu_dm_from_mirror(gpu, host, rows, cols, id) bind(c, name="c_bpack_gpu_dm_from_mirror")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, host
         integer(c_int) :: rows, cols, id
      end subroutine c_bpack_gpu_dm_from_mirror

      subroutine c_bpack_gpu_dm_give(gpu, id, dev) bind(c, name="c_bpack_gpu_dm_give")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: id
         type(c_ptr) :: dev
      end subroutine c_bpack_gpu_dm_give

      subroutine c_bpack_gpu_dm_keep(gpu, id, host, stale) bind(c, name="c_bpack_gpu_dm_keep")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu, host
         integer(c_int) :: id, stale
      end subroutine c_bpack_gpu_dm_keep

      subroutine c_bpack_hodlr_gpu_baca_end(gpu, out) bind(c, name="c_bpack_hodlr_gpu_baca_end")
         import :: c_ptr, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: out(6)
      end subroutine c_bpack_hodlr_gpu_baca_end

      subroutine c_bpack_gpu_dm_flops(gpu, flops) bind(c, name="c_bpack_gpu_dm_flops")
         import :: c_ptr, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: flops
      end subroutine c_bpack_gpu_dm_flops
      !>**** The rank-reveal of one level of a butterfly construction on the
      !> device (GPU_BACKEND/bf_gpu).  `entry_off` is each panel's offset into
      !> `entries`, in elements; `transposed` is 1 for the 'C' half.  The
      !> results come back in block order: `rank` one per block, `jpvt` the
      !> cumulative sum of the panels' column counts, `xout` the cumulative sum
      !> of m*n.  info 1 means the device declined (no room, or a precision
      !> without the butterfly backend) and the caller runs the CPU path.
      subroutine c_bf_gpu_level_id(nblk, rows, cols, blk, transposed, entries, nnz, rtol, atol, &
                                   rank, jpvt, xout, info) bind(c, name="c_bf_gpu_level_id")
         import :: c_int, c_int64_t, c_double, c_ptr
         integer(c_int) :: nblk, transposed
         integer(c_int) :: rows(*), cols(*)
         type(c_ptr) :: blk(*)
         integer(c_int64_t) :: nnz
         !>**** the entry buffer and the output X cross as addresses, so this
         !> binding is the same in every precision; the backend is real double
         !> for now and declines (info 1) in the others
         type(c_ptr), value :: entries, xout
         real(c_double) :: rtol, atol
         integer(c_int) :: rank(*), jpvt(*), info
      end subroutine c_bf_gpu_level_id

      !>**** c_bf_gpu_level_id on the entries the butterfly's GPU evaluator
      !> (handle) left on the device; `entries` is the host buffer they mirror,
      !> filled only when the device declines (info 1)
      subroutine c_bf_gpu_level_id_kept(handle, nblk, rows, cols, blk, transposed, entries, nnz, rtol, atol, &
                                        rank, jpvt, xout, info) bind(c, name="c_bf_gpu_level_id_kept")
         import :: c_int, c_int64_t, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int) :: nblk, transposed
         integer(c_int) :: rows(*), cols(*)
         type(c_ptr) :: blk(*)
         integer(c_int64_t) :: nnz
         type(c_ptr), value :: entries, xout
         real(c_double) :: rtol, atol
         integer(c_int) :: rank(*), jpvt(*), info
      end subroutine c_bf_gpu_level_id_kept

      !>**** A butterfly matvec on the caller's host arrays, run on the device
      !> (GPU_BACKEND/bf_gpu): y = op(bf) x, the factors resident between
      !> calls.  Collective over the butterfly's ranks.  info 1 means the
      !> device declined and the caller runs the CPU path.
      subroutine c_bf_gpu_mult_host(bf, ptree, op, nrhs, x, ldx, y, ldy, flops, info) &
         bind(c, name="c_bf_gpu_mult_host")
         import :: c_int, c_double, c_ptr, c_char
         type(c_ptr), value :: bf, ptree, x, y
         character(kind=c_char) :: op
         integer(c_int) :: nrhs, ldx, ldy, info
         real(c_double) :: flops
      end subroutine c_bf_gpu_mult_host

      !>**** The entries of a list of intersections of a butterfly, on the
      !> device (GPU_BACKEND/bf_gpu/bf_extract.hpp), in C_BF_ExtractElement's
      !> arguments (columns 1..N, pgidx 0 based).  info 1: the device declined
      !> (a butterfly on several ranks, a block another rank owns, no room)
      !> and the caller runs the host path.
      subroutine c_bf_gpu_extract(bf, ptree, ninter, allrows, allcols, rowidx, colidx, pgidx, npmap, pmaps, &
                                  alldat_loc, nalldat_loc, flops, info) bind(c, name="c_bf_gpu_extract")
         import :: c_int, c_int64_t, c_double, c_ptr
         type(c_ptr), value :: bf, ptree, alldat_loc
         integer(c_int) :: ninter, npmap, info
         integer(c_int) :: allrows(*), allcols(*), rowidx(*), colidx(*), pgidx(*), pmaps(*)
         integer(c_int64_t) :: nalldat_loc
         real(c_double) :: flops
      end subroutine c_bf_gpu_extract

      !>**** The handle of a butterfly's GPU evaluator (c_bf_set_gpu_*), null
      !> when it has none, for its M rows and N columns; xyz: the butterfly
      !> mesh's coordinates xyz(1:dim, -N:M), or null.  warm: the evaluator's
      !> first-use costs, paid here.
      subroutine c_bf_gpu_entries(bf, M, N, dim, xyz, handle, warm) bind(c, name="c_bf_gpu_entries")
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: bf, xyz
         integer(c_int) :: M, N, dim
         type(c_ptr) :: handle
         real(c_double) :: warm
      end subroutine c_bf_gpu_entries

      !>**** BF_randomized (less_adapt 1) on the device, for BF_randomized_gpu
      !> (GPU_BACKEND/bf_gpu/bf_randomized.hpp): a butterfly of L >= 1 levels
      !> on its process group, laid out as bf_rand (BF_Init_randomized's
      !> structure, no factors yet): this rank's blocks of every level, its
      !> leaves' ends ms and ns as BF_Init_randomized counts them, over its U
      !> and its V blocks, and the global index (zero based) of its first row
      !> and column.  The black box is the device one registered for bf
      !> (c_bf_set_gpu_matvec), else the C callback matvec with its quant.
      !> rows/cols: the shapes of this rank's blocks - V's, Kerl(1:L)'s row by
      !> row, U's - for c_bf_gpu_randomized_fetch, which takes the handle.
      !> prof: seconds of all, the black box, the sweeps, the range finding
      !> and solves, the error test; then the black box's calls and columns,
      !> the tries of the rank guess and the test's error.  info 1: the device
      !> declined, on every rank together.  Collective over the group.
      subroutine c_bf_gpu_randomized(bf, bf_rand, ptree, L, level_half, M, N, ms, ns, row0, col0, rank0, rankrate, &
                                     itermax, nbundle, tol_rdetect, tol_ls, matvec, quant, handle, rows, cols, prof, &
                                     info) bind(c, name="c_bf_gpu_randomized")
         import :: c_int, c_double, c_ptr, c_funptr
         type(c_ptr), value :: bf, bf_rand, ptree, quant
         type(c_funptr), value :: matvec
         integer(c_int) :: L, level_half, M, N, row0, col0, rank0, itermax, nbundle, info
         integer(c_int) :: ms(*), ns(*), rows(*), cols(*)
         real(c_double) :: rankrate, tol_rdetect, tol_ls, prof(*)
         type(c_ptr) :: handle
      end subroutine c_bf_gpu_randomized

      !>**** The blocks of c_bf_gpu_randomized into ptrs (null: skipped), and
      !> the handle released; seconds: the time of the copy
      subroutine c_bf_gpu_randomized_fetch(handle, ptrs, seconds) bind(c, name="c_bf_gpu_randomized_fetch")
         import :: c_double, c_ptr
         type(c_ptr), value :: handle
         type(c_ptr) :: ptrs(*)
         real(c_double) :: seconds
      end subroutine c_bf_gpu_randomized_fetch

      !>**** Drop what a butterfly holds on the device; before it is freed
      subroutine c_bf_gpu_release(bf) bind(c, name="c_bf_gpu_release")
         import :: c_ptr
         type(c_ptr), value :: bf
      end subroutine c_bf_gpu_release

   end interface

   !> the device indices of the forward blocks of a rank (HODLR_gpu_block_table)
   type hodlr_gpu_itab
      integer, allocatable :: i(:)
   end type hodlr_gpu_itab

contains

   !>**** Whether a butterfly level's rank-reveal runs on the device
   !> (GPU_BACKEND/bf_gpu): the options ask for the GPU, this library has the
   !> backend, the block really is a butterfly (a low-rank block goes through
   !> BF_gpu_construct instead, see BF_gpu_construct_ok), and the level's
   !> entries arrive in one contiguous buffer - which they do not when there
   !> are masks, so a masked block stays on the host.
   !>
   !> EITHER option asks for it, because a butterfly block reaches here from
   !> two places: coupling two HODLRs, where it inherits HODLR_use_gpu, or as
   !> the F12/F21 of an H2 front, where the format is 7 and the caller sets
   !> H2_use_gpu instead.  A STRUMPACK H2 front does the latter and
   !> deliberately leaves HODLR_use_gpu at 0, because HODLR_gpu_check_options
   !> requires format 1 and LRlevel 0 - neither true of an H2 with butterfly
   !> couplings - so reading HODLR_use_gpu alone would never fire there.
   logical function BF_gpu_rankreveal_ok(blocks, option, Nboundall)
      implicit none
      type(matrixblock)::blocks
      type(Hoption)::option
      integer::Nboundall

      BF_gpu_rankreveal_ok = .false.
      if (option%HODLR_use_gpu <= 0 .and. option%H2_use_gpu <= 0) return
      if (blocks%level_butterfly == 0) return
      if (Nboundall /= 0) return
      if (.not. BPACK_gpu_available()) return
      BF_gpu_rankreveal_ok = .true.
   end function BF_gpu_rankreveal_ok

   !>**** Whether a butterfly's matvec runs on the device (GPU_BACKEND/bf_gpu):
   !> the options ask for the GPU and this library has the backend.  A
   !> low-rank block (level_butterfly 0) is one too: on several ranks they all
   !> share its one leaf, which the device sums and broadcasts as the host
   !> does.  The device copy of the factors is made at the first call and kept.
   !> Either option asks for it, as for BF_gpu_rankreveal_ok: a STRUMPACK H2
   !> front sets H2_use_gpu and leaves HODLR_use_gpu at 0.
   logical function BF_gpu_mult_ok(blocks, option)
      implicit none
      type(matrixblock)::blocks
      type(Hoption)::option

      BF_gpu_mult_ok = .false.
      if (option%HODLR_use_gpu <= 0 .and. option%H2_use_gpu <= 0) return
      if (.not. BPACK_gpu_available()) return
      BF_gpu_mult_ok = .true.
   end function BF_gpu_mult_ok

   !>**** Whether a butterfly's entry extraction is tried on the device
   !> (GPU_BACKEND/bf_gpu/bf_extract.hpp): as for the matvec, a butterfly
   !> rather than a low-rank block.  The device itself declines what it does
   !> not handle yet - a butterfly on several ranks, a block another rank owns
   !> - and the host path runs.
   logical function BF_gpu_extract_ok(blocks, option)
      implicit none
      type(matrixblock)::blocks
      type(Hoption)::option

      BF_gpu_extract_ok = BF_gpu_mult_ok(blocks, option) .and. blocks%level_butterfly > 0
   end function BF_gpu_extract_ok

   !>**** Whether a butterfly's construction from matvecs runs on the device
   !> (BF_randomized_gpu): as for the matvec, a butterfly rather than a
   !> low-rank block, with the less adaptive scheme (less_adapt 1; the other
   !> one is BF_randomized_old's, as is a low-rank block's).
   logical function BF_gpu_sketch_ok(blocks, option, ptree)
      implicit none
      type(matrixblock)::blocks
      type(Hoption)::option
      type(proctree)::ptree

      BF_gpu_sketch_ok = .false.
      if (.not. BF_gpu_mult_ok(blocks, option)) return
      if (blocks%level_butterfly == 0) return
      if (option%less_adapt /= 1) return
      if (.not. IOwnPgrp(ptree, blocks%pgno)) return
      BF_gpu_sketch_ok = .true.
   end function BF_gpu_sketch_ok

   !> Whether this library has the GPU backend (the double or double complex
   !> library of a build with enable_h2_gpu=ON)
   logical function BPACK_gpu_available()
      integer(c_int) :: available
      call c_bpack_gpu_available(available)
      BPACK_gpu_available = (available == 1)
   end function BPACK_gpu_available

   !> Stop with a clear message unless the options allow the HODLR GPU backend
   subroutine HODLR_gpu_check_options(option, ptree)
      type(Hoption)::option
      type(proctree)::ptree

      call BPACK_require(option%HODLR_use_gpu >= 0 .and. option%HODLR_use_gpu <= 2, 'HODLR_use_gpu must be 0, 1 or 2')
      if (option%HODLR_use_gpu == 0) return
      call BPACK_require(BPACK_gpu_available(), &
         'HODLR_use_gpu>0 requires the double or double complex library of a build with enable_h2_gpu=ON')
      call BPACK_require(option%format == HODLR, 'HODLR_use_gpu>0 requires format=1 (HODLR)')
      call BPACK_require(option%LRlevel == 0, 'HODLR_use_gpu>0 requires LRlevel=0 (HODBF has no GPU path)')
      call BPACK_require(option%use_zfp /= 1, 'HODLR_use_gpu>0 does not support ZFP-compressed dense blocks')
      call BPACK_require(option%HODLR_gpu_pieces >= 1, 'HODLR_gpu_pieces must be at least 1')
   end subroutine HODLR_gpu_check_options

   !> Give bmat a GPU state for a HODLR GPU run (HODLR_use_gpu > 0) and let
   !> its HODLR borrow it.  Called at the entry points that receive bmat.
   subroutine BPACK_gpu_bind(bmat, option, ptree)
      type(Bmatrix)::bmat
      type(Hoption)::option
      type(proctree)::ptree

      integer(c_int) :: nproc, share
      integer :: share_max, ierr

      if (option%HODLR_use_gpu <= 0 .or. option%format /= HODLR) return
      call HODLR_gpu_check_options(option, ptree)
      nproc = ptree%nproc
      ! ranks sharing a GPU split its memory (collective, once), then the
      ! exchange arena (before the first device allocation; its default size
      ! split among them too)
      call c_bpack_gpu_share_device(ptree%Comm, share)
      call c_bpack_gpu_init_exchange(nproc, share)
      share_max = share
      call MPI_ALLREDUCE(MPI_IN_PLACE, share_max, 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (share_max > 1 .and. ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, '(A,I3,A)') ' HODLR GPU: up to', share_max, ' ranks share a GPU; each takes an equal part of its memory'
      endif
      if (.not. c_associated(bmat%gpu)) call c_bpack_gpu_create(bmat%gpu)
      if (associated(bmat%ho_bf)) bmat%ho_bf%gpu = bmat%gpu
   end subroutine BPACK_gpu_bind

   !> Free the GPU state of bmat (device copies and registered kernel)
   !> The low-rank entry term of ker (c_bpack_set_entry_addend) to the device for the GPU construction of
   !> bmat, which adds it to the entries it takes from the host (HodlrConstruct::host_fill; a construction
   !> by a GPU evaluator cannot take one): this rank's rows, and while a shared block is built the rows its
   !> pieces need (HODLR_gpu_addend_prefetch).  Else the host adds it (element_Zmn_blocklist_user)
   subroutine HODLR_gpu_addend_begin(bmat, option, ker)
      type(Bmatrix)::bmat
      type(Hoption)::option
      type(kernelquant)::ker
      integer(c_int) :: ready

      if (.not. allocated(ker%addend)) return
      ker%addend%on_device = .false.
      if (option%format /= HODLR .or. option%HODLR_use_gpu <= 0 .or. .not. c_associated(bmat%gpu)) return
      call c_bpack_hodlr_gpu_construct_ready(bmat%gpu, ready)
      call BPACK_require(ready /= 1, 'a low-rank entry term (c_bpack_set_entry_addend) needs the entries from the host, '// &
                         'not a GPU evaluator')
      if (ready /= 2) return
      call BPACK_require(size(ker%addend%vmap) == size(ker%addend%umap), 'a low-rank entry term of a square matrix')
      ker%addend%on_device = .true.
      call HODLR_gpu_addend_sync(bmat%gpu, ker%addend, option)
   end subroutine HODLR_gpu_addend_begin

   !> The rows of the low-rank entry term ad held now (with their slots) to the device, U scaled as the host
   !> entries (scale_factor)
   subroutine HODLR_gpu_addend_sync(gpu, ad, option)
      type(c_ptr) :: gpu
      type(entry_addend)::ad
      type(Hoption)::option
      DT, allocatable, target :: ut(:, :), vt(:, :)
      integer(c_int), allocatable, target :: um(:), vm(:)
      integer(c_int) :: nu, nv, r
      integer(c_int64_t) :: n

      if (.not. ad%on_device) return
      allocate (ut(max(ad%nu, 1), ad%r), vt(max(ad%nv, 1), ad%r), um(size(ad%umap)), vm(size(ad%vmap)))
      if (ad%nu > 0) ut(1:ad%nu, :) = ad%u(1:ad%nu, 1:ad%r)*option%scale_factor
      if (ad%nv > 0) vt(1:ad%nv, :) = ad%v(1:ad%nv, 1:ad%r)
      um = ad%umap
      vm = ad%vmap
      n = size(ad%umap)
      nu = ad%nu
      nv = ad%nv
      r = ad%r
      call c_bpack_hodlr_gpu_set_addend(gpu, n, nu, nv, r, c_loc(ut), c_loc(vt), c_loc(um), c_loc(vm))
      deallocate (ut, vt, um, vm)
   end subroutine HODLR_gpu_addend_sync

   !> The device copy of the low-rank entry term freed (after the construction of bmat)
   subroutine HODLR_gpu_addend_end(bmat, ker)
      type(Bmatrix)::bmat
      type(kernelquant)::ker
      integer(c_int) :: nu, nv, r
      integer(c_int64_t) :: n

      if (.not. allocated(ker%addend)) return
      if (.not. ker%addend%on_device) return
      n = 0
      nu = 0
      nv = 0
      r = 0
      call c_bpack_hodlr_gpu_set_addend(bmat%gpu, n, nu, nv, r, c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr)
      ker%addend%on_device = .false.
   end subroutine HODLR_gpu_addend_end

   subroutine BPACK_gpu_delete(bmat)
      type(Bmatrix)::bmat

      if (c_associated(bmat%gpu)) call c_bpack_gpu_delete(bmat%gpu)
      bmat%gpu = c_null_ptr
      if (associated(bmat%ho_bf)) bmat%ho_bf%gpu = c_null_ptr
   end subroutine BPACK_gpu_delete

   !> First local row (global, 1-based) and local row count of this rank
   subroutine HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      type(hobf)::ho_bf1
      type(proctree)::ptree
      integer idx_start_glo, pp
      integer(c_int64_t) :: n_loc
      type(matrixblock), pointer::root

      root => ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)
      pp = ptree%MyID - ptree%pgrp(root%pgno)%head + 1
      idx_start_glo = root%headn + root%N_p(pp, 1) - 1
      n_loc = root%N_p(pp, 2) - root%N_p(pp, 1) + 1
   end subroutine HODLR_gpu_local_rows

   !> Upload the forward blocks of ho_bf1 (dense leaves and low-rank blocks)
   !> to the GPU; the symmetric HODLR stores A21 only (A12 is its transpose).
   !> With several ranks, a block of a shared node (one whose process group
   !> has more than one rank) goes up as the rows of U this rank owns and the
   !> rows of V it owns once V is moved to the ranks that own the rows it
   !> multiplies (the child of the other side); the node's communicator comes
   !> with its level.
   subroutine HODLR_gpu_upload_forward(ho_bf1, option, stats, ptree)
      type(hobf), target::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type dtbuf
         DT, allocatable :: a(:, :)
      end type dtbuf
      type(matrixblock), pointer::blk, colchild
      type(dtbuf), allocatable, target :: vstore(:)
      DT, target :: dummy(1, 1)
      DT, allocatable :: vsrc(:, :)
      integer level, ii, nd, pp, idx_start_glo, ierr, pgno, nsrc, nv, kk, crank
      logical shared, devv, check
      integer(c_int) :: src_id, dst_id, rows_c, cols_c
      type(c_ptr) :: pvd
      integer(c_int) :: level_c, m, n, k, ldu, ldv, sym, maxlevel_c
      integer(c_int64_t) :: n_loc, row0, col0
      type(c_ptr) :: pu, pv
      real(c_double) :: mbytes, mirrored
      real(kind=8) :: t1, t2

      call assert(c_associated(ho_bf1%gpu), 'HODLR_gpu_upload_forward: the HODLR has no GPU state')
      t1 = MPI_Wtime()
      call HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      maxlevel_c = ho_bf1%Maxlevel
      call c_bpack_hodlr_gpu_reset(ho_bf1%gpu, n_loc, maxlevel_c)

      level = ho_bf1%Maxlevel + 1
      do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
         blk => ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1)
         if (.not. IOwnPgrp(ptree, blk%pgno)) cycle
         call assert(associated(blk%fullmat), 'HODLR_gpu_upload_forward: a dense leaf has no entries')
         row0 = blk%headm - idx_start_glo
         m = blk%M
         ldu = size(blk%fullmat, 1)
         call c_bpack_hodlr_gpu_add_leaf(ho_bf1%gpu, row0, m, c_loc(blk%fullmat(1, 1)), ldu)
      enddo

      sym = 0
      if (option%sym > 0) sym = 1
      check = HODLR_gpu_check_level() > 0
      allocate (vstore(2*max(1, ho_bf1%Maxlevel)))
      nv = 0
      dummy = 0
      do level = 1, ho_bf1%Maxlevel
         do nd = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            pgno = ho_bf1%levels(level)%BP_inverse(nd)%pgno
            if (.not. IOwnPgrp(ptree, pgno)) cycle
            shared = ptree%pgrp(pgno)%nproc > 1
            level_c = level
            if (shared) then
               ! (the device's head of a node, rank 0 of its communicator, is the CPU's)
               call MPI_COMM_RANK(ptree%pgrp(pgno)%Comm, crank, ierr)
               call assert((crank == 0) .eqv. (ptree%MyID == ptree%pgrp(pgno)%head), &
                  'HODLR_gpu_upload_forward: rank 0 of a node communicator is not the group head')
               call c_bpack_hodlr_gpu_add_shared_level(ho_bf1%gpu, level_c, ptree%pgrp(pgno)%Comm)
            endif
            do ii = nd*2 - 1, nd*2
               if (sym == 1 .and. mod(ii, 2) == 1) cycle
               blk => ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1)
               ! the child whose rows the block's columns are
               if (mod(ii, 2) == 1) then
                  colchild => ho_bf1%levels(level + 1)%BP_inverse(nd*2)%LL(1)%matrices_block(1)
               else
                  colchild => ho_bf1%levels(level + 1)%BP_inverse(nd*2 - 1)%LL(1)%matrices_block(1)
               endif
               kk = 0
               if (IOwnPgrp(ptree, blk%pgno)) then
                  call assert(blk%level_butterfly == 0, 'HODLR_gpu_upload_forward: a block is not low rank (LRlevel>0)')
                  kk = size(blk%ButterflyU%blocks(1)%matrix, 2)
                  call assert(size(blk%ButterflyV%blocks(1)%matrix, 2) == kk, &
                     'HODLR_gpu_upload_forward: U and V ranks differ')
               endif
               if (shared) call MPI_ALLREDUCE(MPI_IN_PLACE, kk, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
               k = kk
               ! U: the rows this rank owns
               m = 0
               row0 = 0
               ldu = 1
               pu = c_loc(dummy(1, 1))
               if (IOwnPgrp(ptree, blk%pgno)) then
                  pp = ptree%MyID - ptree%pgrp(blk%pgno)%head + 1
                  row0 = blk%headm + blk%M_p(pp, 1) - 1 - idx_start_glo
                  m = blk%M_loc
                  ldu = max(1, size(blk%ButterflyU%blocks(1)%matrix, 1))
                  if (m > 0 .and. k > 0) pu = c_loc(blk%ButterflyU%blocks(1)%matrix(1, 1))
               endif
               ! V: the rows this rank owns, after moving V to the column child's ranks
               n = 0
               col0 = 0
               ldv = 1
               pv = c_loc(dummy(1, 1))
               pvd = c_null_ptr
               ! a shared block's V moves on the device when all ranks of the node have its device copy
               ! (the construction's, kept for the upload) and the block's columns are the column child's rows
               devv = .false.
               if (shared .and. k > 0) then
                  devv = blk%headn == colchild%headm
                  src_id = 0
                  if (devv .and. IOwnPgrp(ptree, blk%pgno)) then
                     if (blk%N_loc > 0) then
                        rows_c = blk%N_loc
                        cols_c = k
                        call c_bpack_gpu_dm_from_mirror(ho_bf1%gpu, c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1)), rows_c, &
                                                        cols_c, src_id)
                        devv = src_id /= 0
                     endif
                  endif
                  call MPI_ALLREDUCE(MPI_IN_PLACE, devv, 1, MPI_LOGICAL, MPI_LAND, ptree%pgrp(pgno)%Comm, ierr)
                  call assert(devv .or. .not. ho_bf1%gpu_host_stale, &
                     'HODLR_gpu_upload_forward: a shared block has neither its host V nor its device copy')
                  if (devv) then
                     dst_id = 0
                     if (IOwnPgrp(ptree, colchild%pgno)) then
                        rows_c = colchild%M_loc
                        cols_c = k
                        call c_bpack_gpu_dm_alloc(ho_bf1%gpu, rows_c, cols_c, c_null_ptr, dst_id)
                     endif
                     call HODLR_gpu_dm_redist(ho_bf1%gpu, src_id, blk%N_p, blk%pgno, dst_id, 0, colchild%M_p, &
                                              colchild%pgno, k, ptree)
                     if (dst_id /= 0) call c_bpack_gpu_dm_give(ho_bf1%gpu, dst_id, pvd)
                  endif
                  if (src_id /= 0) call c_bpack_gpu_dm_free(ho_bf1%gpu, src_id)
               endif
               if (.not. shared) then
                  pp = ptree%MyID - ptree%pgrp(blk%pgno)%head + 1
                  col0 = blk%headn + blk%N_p(pp, 1) - 1 - idx_start_glo
                  n = blk%N_loc
                  ldv = max(1, size(blk%ButterflyV%blocks(1)%matrix, 1))
                  if (n > 0 .and. k > 0) pv = c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1))
               elseif (devv .and. .not. check) then
                  if (IOwnPgrp(ptree, colchild%pgno)) then
                     pp = ptree%MyID - ptree%pgrp(colchild%pgno)%head + 1
                     col0 = colchild%headm + colchild%M_p(pp, 1) - 1 - idx_start_glo
                     n = colchild%M_loc
                     ldv = max(1, colchild%M_loc)
                  endif
               elseif (k > 0) then  ! (on the host; with BPACK_CHECK=hodlr also next to the device's, compared at the upload)
                  nsrc = 0
                  if (IOwnPgrp(ptree, blk%pgno)) nsrc = blk%N_loc
                  allocate (vsrc(max(1, nsrc), k))
                  vsrc = 0
                  if (nsrc > 0) vsrc(1:nsrc, 1:k) = blk%ButterflyV%blocks(1)%matrix(1:nsrc, 1:k)
                  nv = nv + 1
                  allocate (vstore(nv)%a(max(1, colchild%M_loc), k))
                  vstore(nv)%a = 0
                  call Redistribute1Dto1D(vsrc, max(1, nsrc), blk%N_p, blk%headn, blk%pgno, &
                     vstore(nv)%a, max(1, colchild%M_loc), colchild%M_p, colchild%headm, colchild%pgno, k, ptree)
                  deallocate (vsrc)
                  if (IOwnPgrp(ptree, colchild%pgno)) then
                     pp = ptree%MyID - ptree%pgrp(colchild%pgno)%head + 1
                     col0 = colchild%headm + colchild%M_p(pp, 1) - 1 - idx_start_glo
                     n = colchild%M_loc
                     ldv = max(1, colchild%M_loc)
                     if (n > 0) pv = c_loc(vstore(nv)%a(1, 1))
                  endif
               endif
               call c_bpack_hodlr_gpu_add_lowrank(ho_bf1%gpu, level_c, row0, m, col0, n, k, pu, ldu, pv, ldv, sym, pvd)
            enddo
         enddo
      enddo

      call c_bpack_hodlr_gpu_commit_forward(ho_bf1%gpu, mbytes, mirrored)
      do ii = 1, nv
         deallocate (vstore(ii)%a)
      enddo
      deallocate (vstore)
      t2 = MPI_Wtime()
      call MPI_ALLREDUCE(MPI_IN_PLACE, mbytes, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, mirrored, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, '(A,F10.1,A,F10.1,A,F8.3,A)') ' HODLR GPU: forward blocks on the device, ', mbytes, ' MB (', &
            mirrored, ' MB kept from the construction) in ', t2 - t1, ' s'
      endif
      ! (BPACK_CHECK=hodlr: the way back to the host, compared with the host copies; then taken, the host
      ! copies overwritten by the device's, which the CPU comparisons that follow use)
      if (check .and. .not. ho_bf1%gpu_host_stale) then
         call HODLR_gpu_fetch_host(ho_bf1, option, ptree, compare=.true.)
         ho_bf1%gpu_host_stale = .true.
         call HODLR_gpu_fetch_host(ho_bf1, option, ptree)
      endif
   end subroutine HODLR_gpu_upload_forward

   !> The host copies of the low-rank factors of ho_bf1 marked unfilled
   !> (gpu_host_stale), all of them, from the forward blocks on the device:
   !> the blocks in the order of HODLR_gpu_upload_forward, a shared block's V
   !> moved back from its column child's layout to the block's (collective
   !> over the ranks of ptree).  compare: into temporaries, compared with the
   !> host copies (which must be filled).  For BPACK_CHECK=hodlr
   subroutine HODLR_gpu_fetch_host(ho_bf1, option, ptree, compare)
      type(hobf), target::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      logical, optional :: compare
      type(matrixblock), pointer::blk, colchild
      DT, allocatable, target :: tu(:, :), tv(:, :)
      DT :: one
      integer level, nd, ii, pgno, ierr, nb, nbad
      logical shared, cmp, own
      integer(c_int) :: idx, m, n, k, ldu, ldv, src_id, dst_id, rows_c, cols_c
      type(c_ptr) :: pu, pv
      real(kind=8) :: t1, t2, mbytes

      cmp = .false.
      if (present(compare)) cmp = compare
      if (.not. (ho_bf1%gpu_host_stale .or. cmp)) return
      call assert(c_associated(ho_bf1%gpu), 'HODLR_gpu_fetch_host: the HODLR has no GPU state')
      t1 = MPI_Wtime()
      mbytes = 0
      nb = 0
      nbad = 0
      idx = 0
      do level = 1, ho_bf1%Maxlevel
         do nd = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            pgno = ho_bf1%levels(level)%BP_inverse(nd)%pgno
            if (.not. IOwnPgrp(ptree, pgno)) cycle
            shared = ptree%pgrp(pgno)%nproc > 1
            do ii = nd*2 - 1, nd*2
               if (option%sym > 0 .and. mod(ii, 2) == 1) cycle  ! (A12 is a view of A21)
               blk => ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1)
               if (mod(ii, 2) == 1) then
                  colchild => ho_bf1%levels(level + 1)%BP_inverse(nd*2)%LL(1)%matrices_block(1)
               else
                  colchild => ho_bf1%levels(level + 1)%BP_inverse(nd*2 - 1)%LL(1)%matrices_block(1)
               endif
               own = IOwnPgrp(ptree, blk%pgno)
               ! U, and V of a block of one rank, straight into the host arrays (or temporaries)
               pu = c_null_ptr
               pv = c_null_ptr
               ldu = 1
               ldv = 1
               if (own) then
                  k = size(blk%ButterflyU%blocks(1)%matrix, 2)
                  ldu = max(1, size(blk%ButterflyU%blocks(1)%matrix, 1))
                  ldv = max(1, size(blk%ButterflyV%blocks(1)%matrix, 1))
                  call assert(ldu >= blk%M_loc .and. ldv >= blk%N_loc .and. size(blk%ButterflyV%blocks(1)%matrix, 2) == k, &
                              'HODLR_gpu_fetch_host: unexpected host factors')
                  ! (a shared block's V comes down whole: N_loc rows, leading dimension N_loc)
                  if (shared .and. blk%N_loc > 0) call assert(ldv == blk%N_loc, 'HODLR_gpu_fetch_host: unexpected host V')
                  allocate (tu(ldu, k), tv(ldv, k))
                  if (blk%M_loc > 0 .and. k > 0) then
                     pu = c_loc(blk%ButterflyU%blocks(1)%matrix(1, 1))
                     if (cmp) pu = c_loc(tu(1, 1))
                  endif
                  if (.not. shared .and. blk%N_loc > 0 .and. k > 0) then
                     pv = c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1))
                     if (cmp) pv = c_loc(tv(1, 1))
                  endif
               endif
               call c_bpack_hodlr_gpu_fetch_block(ho_bf1%gpu, idx, pu, ldu, pv, ldv, m, n, k)
               if (own) then
                  call assert(m == blk%M_loc .and. k == size(blk%ButterflyU%blocks(1)%matrix, 2), &
                     'HODLR_gpu_fetch_host: a block differs in shape from its device copy')
                  mbytes = mbytes + dble(blk%M_loc + blk%N_loc)*k*(storage_size(one)/8)/1d6
               endif
               ! V of a shared block: back from the column child's layout (the device's) to the block's
               if (shared .and. k > 0) then
                  src_id = 0
                  call c_bpack_hodlr_gpu_block_v_dm(ho_bf1%gpu, idx, src_id)
                  dst_id = 0
                  if (own) then
                     rows_c = blk%N_loc
                     cols_c = k
                     call c_bpack_gpu_dm_alloc(ho_bf1%gpu, rows_c, cols_c, c_null_ptr, dst_id)
                  endif
                  call HODLR_gpu_dm_redist(ho_bf1%gpu, src_id, colchild%M_p, colchild%pgno, dst_id, 0, blk%N_p, blk%pgno, &
                                           int(k), ptree)
                  if (dst_id /= 0) then
                     if (blk%N_loc > 0) then
                        pv = c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1))
                        if (cmp) pv = c_loc(tv(1, 1))
                        call c_bpack_gpu_dm_download(ho_bf1%gpu, dst_id, pv)
                     endif
                     call c_bpack_gpu_dm_free(ho_bf1%gpu, dst_id)
                  endif
                  if (src_id /= 0) call c_bpack_gpu_dm_free(ho_bf1%gpu, src_id)
               endif
               if (own) then
                  if (cmp) then
                     nb = nb + 1
                     if (any(tu(1:blk%M_loc, 1:k) /= blk%ButterflyU%blocks(1)%matrix(1:blk%M_loc, 1:k)) .or. &
                         any(tv(1:blk%N_loc, 1:k) /= blk%ButterflyV%blocks(1)%matrix(1:blk%N_loc, 1:k))) nbad = nbad + 1
                  endif
                  deallocate (tu, tv)
               endif
               idx = idx + 1
            enddo
         enddo
      enddo
      ho_bf1%gpu_host_stale = .false.
      t2 = MPI_Wtime()
      call MPI_ALLREDUCE(MPI_IN_PLACE, mbytes, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (cmp) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, nb, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, nbad, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
         if (ptree%MyID == Main_ID) write (*, '(A,I6,A,I6,A,F10.1,A,F8.3,A)') ' HODLR GPU check (fetch): ', nbad, ' of', nb, &
            ' low-rank blocks brought back from the device differ from the host (', mbytes, ' MB in', t2 - t1, ' s)'
      elseif (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, '(A,F10.1,A,F8.3,A)') ' HODLR GPU: host copies of the low-rank factors from the device, ', mbytes, &
            ' MB in ', t2 - t1, ' s'
      endif
   end subroutine HODLR_gpu_fetch_host

   !> Vout = op(A) Vin with the forward blocks on the GPU (the whole HODLR,
   !> levels 1 to Maxlevel+1), as HODLR_Mult does on the CPU
   subroutine HODLR_gpu_Mult(trans, Ns, num_vectors, Vin, Vout, ho_bf1, ptree, option, stats)
      character trans
      integer Ns, num_vectors
      DT, target::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hoption)::option
      type(Hstat)::stats
      character(kind=c_char) :: op
      integer(c_int) :: nrhs, mode
      real(c_double) :: flops

      call assert(c_associated(ho_bf1%gpu), 'HODLR_gpu_Mult: the HODLR has no GPU state')
      op = 'N'
      if (option%sym <= 0 .and. trans /= 'N') op = 'T'
#if DAT==0 || DAT==2
      if (trans == 'C') Vin = conjg(Vin)
#endif
      nrhs = num_vectors
      mode = option%HODLR_use_gpu
      call c_bpack_hodlr_gpu_mult(ho_bf1%gpu, op, nrhs, c_loc(Vin(1, 1)), c_loc(Vout(1, 1)), mode, flops)
#if DAT==0 || DAT==2
      if (trans == 'C') then
         Vout = conjg(Vout)
         Vin = conjg(Vin)
      endif
#endif
      Vout = Vout/option%scale_factor
      stats%Flop_Tmp = flops  ! (this rank's, as HODLR_Mult counts them)
   end subroutine HODLR_gpu_Mult

   !> Symmetric HODLR factorization on the GPU (the device form of
   !> HODLR_factorization_sym, with LU leaves); the factors stay on the device
   subroutine HODLR_gpu_factor_sym(ho_bf1, option, stats, ptree)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      real(c_double) :: out(12), jitter
      integer(c_int) :: mode
      integer level, ii, ierr
      real(kind=8) :: t0, t1
      DTR :: logdet_loc
      DT :: phase_loc

      call assert(c_associated(ho_bf1%gpu), 'HODLR_gpu_factor_sym: the HODLR has no GPU state')
      t0 = MPI_Wtime()
      jitter = option%jitter
      mode = option%HODLR_use_gpu
      call c_bpack_hodlr_gpu_factor_sym(ho_bf1%gpu, jitter, mode, out)
      t1 = MPI_Wtime()

      ! (this rank's part of the log-determinant, summed as the CPU does)
      logdet_loc = out(1)
#if DAT==0 || DAT==2
      phase_loc = cmplx(out(2), out(3), kind=8)
#else
      phase_loc = out(2)
#endif
      call MPI_ALLREDUCE(logdet_loc, ho_bf1%logabsdet, 1, MPI_DTR, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(phase_loc, ho_bf1%phase, 1, MPI_DT, MPI_PROD, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      stats%Time_Factor = t1 - t0
      stats%Flop_Factor = stats%Flop_Factor + out(4)
      stats%Mem_Factor = out(5)*1d6/1024d3  ! (resident factor data, in the CPU's MB: bytes/1024/1000)
      if (.not. allocated(stats%rankmax_of_level_global_factor)) then
         allocate (stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel))
      endif
      stats%rankmax_of_level_global_factor = 0
      do level = 1, ho_bf1%Maxlevel
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            if (.not. IOwnPgrp(ptree, ho_bf1%levels(level)%BP(2*ii)%pgno)) cycle
            stats%rankmax_of_level_global_factor(level) = max(stats%rankmax_of_level_global_factor(level), &
               ho_bf1%levels(level)%BP(2*ii)%LL(1)%matrices_block(1)%rankmax)
         enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel), &
         ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) 'symmetric HODLR logdet:', ho_bf1%phase, ho_bf1%logabsdet
         write (*, '(A,F9.3,A,F8.3,A,F8.3,A,F8.3,A,F9.1,A,F10.1,A)') ' HODLR GPU factor (sym): ', t1 - t0, &
            ' s (leaves ', out(6), ', leaf solves ', out(7), ', nodes ', out(8), '), ', &
            out(4)/max(t1 - t0, 1d-12)/1d9, ' GFLOP/s, ', out(5), ' MB on the device'
         if (out(10) + out(11) > 0) write (*, '(A,I8,A,I8,A,Es11.3)') ' HODLR GPU factor (sym): jitter on ', &
            nint(out(10)), ' leaves and ', nint(out(11)), ' nodes, largest ', out(9)
      endif
   end subroutine HODLR_gpu_factor_sym

   !> Whether the symmetric factors of ho_bf1 are on the GPU
   logical function HODLR_gpu_sym_ready(ho_bf1)
      type(hobf)::ho_bf1
      integer(c_int) :: ready

      HODLR_gpu_sym_ready = .false.
      if (.not. c_associated(ho_bf1%gpu)) return
      call c_bpack_hodlr_gpu_sym_ready(ho_bf1%gpu, ready)
      HODLR_gpu_sym_ready = (ready == 1)
   end function HODLR_gpu_sym_ready

   !> Vout = A^{-1} Vin (op(A) for trans 'C', the symmetric A being its own
   !> transpose) with the symmetric factors on the GPU, as HODLR_Sym_Inv_Apply
   subroutine HODLR_gpu_Sym_Inv_Apply(trans, Ns, num_vectors, Vin, Vout, ho_bf1, mode, flops)
      character trans
      integer Ns, num_vectors, mode
      DT, target::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hobf)::ho_bf1
      real(kind=8) :: flops  ! (out: this rank's)
      integer(c_int) :: nrhs, mode_c
      real(c_double) :: flops_c

      nrhs = num_vectors
      mode_c = mode
#if DAT==0 || DAT==2
      if (trans == 'C') Vin = conjg(Vin)
#endif
      call c_bpack_hodlr_gpu_solve_sym(ho_bf1%gpu, nrhs, c_loc(Vin(1, 1)), c_loc(Vout(1, 1)), mode_c, flops_c)
      flops = flops_c
#if DAT==0 || DAT==2
      if (trans == 'C') then
         Vout = conjg(Vout)
         Vin = conjg(Vin)
      endif
#endif
   end subroutine HODLR_gpu_Sym_Inv_Apply

   !> Compare the symmetric factors on the GPU with those the CPU left in
   !> ho_bf1%levels(:)%SymFactor (BPACK_CHECK=hodlr): the largest relative
   !> differences of Q, G and the LU of S on each level
   subroutine HODLR_gpu_compare_symfactor(ho_bf1, option, ptree)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      DT, allocatable, target::q0(:, :), q1(:, :), g0(:, :), g1(:, :), s(:, :), z0(:, :), z1(:, :)
      integer, allocatable, target::ipiv(:)
      integer level, ii, idx_start_glo, ierr, npiv, nl0, nl1
      logical has0, has1, head
      integer(c_int) :: level_c, k, found, n0_c, n1_c
      integer(c_int64_t) :: c0, n_loc
      real(kind=8) :: err(4), res(2), res_gpu, res_cpu, num, den

      call HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      res = 0
      do level = ho_bf1%Maxlevel, 1, -1
         err = 0
         npiv = 0
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            associate (fac => ho_bf1%levels(level)%SymFactor(ii))
               if (.not. IOwnPgrp(ptree, fac%pgno) .or. fac%rank == 0) cycle
               ! the rows of the children this rank holds (one child on a shared level), G and S on the head
               has0 = allocated(fac%Q0)
               has1 = allocated(fac%Q1)
               head = allocated(fac%S)
               nl0 = 0
               nl1 = 0
               if (has0) nl0 = fac%nloc0
               if (has1) nl1 = fac%nloc1
               k = fac%rank
               allocate (q0(max(1, nl0), k), q1(max(1, nl1), k), z0(max(1, nl0), k), z1(max(1, nl1), k))
               allocate (g0(k, k), g1(k, k), s(k, k), ipiv(k))
               level_c = level
               if (has0) then
                  c0 = fac%head0 - idx_start_glo
               else
                  c0 = fac%head1 - idx_start_glo
               endif
               n0_c = nl0
               n1_c = nl1
               call c_bpack_hodlr_gpu_download_symnode(ho_bf1%gpu, level_c, c0, n0_c, n1_c, c_loc(q0), c_loc(q1), &
                  c_loc(g0), c_loc(g1), c_loc(s), c_loc(ipiv), k, found, c_loc(z0), c_loc(z1))
               call assert(found == 1 .and. k == fac%rank, 'HODLR_gpu_compare_symfactor: node not found on the GPU')
               ! Z: the device bases against the CPU's; Q: the factors
               num = 0
               den = 0
               if (has0) then
                  num = num + fnorm(z0(1:nl0, :) - fac%Z0, nl0, k)**2
                  den = den + fnorm(fac%Z0, nl0, k)**2
               endif
               if (has1) then
                  num = num + fnorm(z1(1:nl1, :) - fac%Z1, nl1, k)**2
                  den = den + fnorm(fac%Z1, nl1, k)**2
               endif
               err(4) = max(err(4), sqrt(num/max(den, tiny(1d0))))
               num = 0
               den = 0
               if (has0) then
                  num = num + fnorm(q0(1:nl0, :) - fac%Q0, nl0, k)**2
                  den = den + fnorm(fac%Q0, nl0, k)**2
               endif
               if (has1) then
                  num = num + fnorm(q1(1:nl1, :) - fac%Q1, nl1, k)**2
                  den = den + fnorm(fac%Q1, nl1, k)**2
               endif
               err(1) = max(err(1), sqrt(num/max(den, tiny(1d0))))
               if (head) then
                  err(2) = max(err(2), sqrt((fnorm(g0 - fac%G0, k, k)**2 + fnorm(g1 - fac%G1, k, k)**2) &
                     /max(fnorm(fac%G0, k, k)**2 + fnorm(fac%G1, k, k)**2, tiny(1d0))))
                  err(3) = max(err(3), fnorm(s - fac%S, k, k)/max(fnorm(fac%S, k, k), tiny(1d0)))
                  if (any(ipiv /= fac%ipiv)) npiv = npiv + 1
               endif
               if (level == ho_bf1%Maxlevel .and. has0 .and. has1) then
                  ! the children are dense leaves: residuals of D Q = Z, GPU and CPU
                  call leaf_residual(ho_bf1%levels(level + 1)%BP(2*ii - 1)%LL(1)%matrices_block(1)%fullmat, &
                     q0, fac%Q0, fac%Z0, res_gpu, res_cpu)
                  res(1) = max(res(1), res_gpu)
                  res(2) = max(res(2), res_cpu)
                  call leaf_residual(ho_bf1%levels(level + 1)%BP(2*ii)%LL(1)%matrices_block(1)%fullmat, &
                     q1, fac%Q1, fac%Z1, res_gpu, res_cpu)
                  res(1) = max(res(1), res_gpu)
                  res(2) = max(res(2), res_cpu)
               endif
               deallocate (q0, q1, z0, z1, g0, g1, s, ipiv)
            end associate
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, err, 4, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, npiv, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
         if (ptree%MyID == Main_ID) write (*, '(A,I3,A,Es10.2,A,Es10.2,A,Es10.2,A,Es10.2,A,I5)') &
            ' HODLR GPU check (sym factor) level', level, ': Q', err(1), ', G', err(2), ', LU(S)', err(3), &
            ', Z vs CPU', err(4), ', nodes with other pivots', npiv
         if (level == ho_bf1%Maxlevel) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, res, 2, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
            if (ptree%MyID == Main_ID) write (*, '(A,I3,A,Es10.2,A,Es10.2)') &
               ' HODLR GPU check (sym factor) level', level, ': largest |D Q - Z|/|Z| on the GPU', res(1), &
               ', on the CPU', res(2)
         endif
      enddo

   contains

      subroutine leaf_residual(d, q_gpu, q_cpu, z, res_gpu, res_cpu)
         DT, pointer::d(:, :)
         DT::q_gpu(:, :), q_cpu(:, :), z(:, :)
         real(kind=8)::res_gpu, res_cpu, zn
         DT, allocatable::r(:, :)

         allocate (r(size(z, 1), size(z, 2)))
         zn = max(fnorm(z, size(z, 1), size(z, 2)), tiny(1d0))
         r = matmul(d, q_gpu) - z
         res_gpu = fnorm(r, size(z, 1), size(z, 2))/zn
         r = matmul(d, q_cpu) - z
         res_cpu = fnorm(r, size(z, 1), size(z, 2))/zn
         deallocate (r)
      end subroutine leaf_residual
   end subroutine HODLR_gpu_compare_symfactor

   !> Unsymmetric HODLR factorization on the GPU (the device form of
   !> HODLR_factorization with LRlevel 0); the factors stay on the device
   subroutine HODLR_gpu_factor_unsym(ho_bf1, option, stats, ptree)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      real(c_double) :: out(12), jitter
      integer(c_int) :: mode
      integer level, ii, ierr
      real(kind=8) :: t0, t1
      DTR :: logdet_loc
      DT :: phase_loc

      call assert(c_associated(ho_bf1%gpu), 'HODLR_gpu_factor_unsym: the HODLR has no GPU state')
      t0 = MPI_Wtime()
      jitter = option%jitter
      mode = option%HODLR_use_gpu
      call c_bpack_hodlr_gpu_factor_unsym(ho_bf1%gpu, jitter, mode, out)
      t1 = MPI_Wtime()

      ! (this rank's part of the log-determinant, summed as the CPU does)
      logdet_loc = out(1)
#if DAT==0 || DAT==2
      phase_loc = cmplx(out(2), out(3), kind=8)
#else
      phase_loc = out(2)
#endif
      call MPI_ALLREDUCE(logdet_loc, ho_bf1%logabsdet, 1, MPI_DTR, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(phase_loc, ho_bf1%phase, 1, MPI_DT, MPI_PROD, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      stats%Time_Factor = t1 - t0
      stats%Flop_Factor = stats%Flop_Factor + out(4)
      stats%Mem_Factor = out(5)*1d6/1024d3  ! (resident factor data, in the CPU's MB: bytes/1024/1000)
      if (.not. allocated(stats%rankmax_of_level_global_factor)) then
         allocate (stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel))
      endif
      stats%rankmax_of_level_global_factor = 0
      do level = 1, ho_bf1%Maxlevel
         do ii = ho_bf1%levels(level)%Bidxs*2 - 1, ho_bf1%levels(level)%Bidxe*2
            if (.not. associated(ho_bf1%levels(level)%BP(ii)%LL)) cycle  ! (freed by a CPU factorization)
            if (.not. IOwnPgrp(ptree, ho_bf1%levels(level)%BP(ii)%pgno)) cycle
            stats%rankmax_of_level_global_factor(level) = max(stats%rankmax_of_level_global_factor(level), &
               ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1)%rankmax)
         enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel), &
         ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) 'logdet:', ho_bf1%phase, ho_bf1%logabsdet
         write (*, '(A,F9.3,A,F8.3,A,F8.3,A,F8.3,A,F9.1,A,F10.1,A)') ' HODLR GPU factor (unsym): ', t1 - t0, &
            ' s (leaves ', out(6), ', block updates ', out(7), ', nodes ', out(8), '), ', &
            out(4)/max(t1 - t0, 1d-12)/1d9, ' GFLOP/s, ', out(5), ' MB on the device'
         if (out(10) + out(11) > 0) write (*, '(A,I8,A,I8,A)') ' HODLR GPU factor (unsym): pivots raised to the '// &
            'jitter threshold in ', nint(out(10)), ' leaves and ', nint(out(11)), ' nodes'
      endif
   end subroutine HODLR_gpu_factor_unsym

   !> Whether the unsymmetric factors of ho_bf1 are on the GPU
   logical function HODLR_gpu_unsym_ready(ho_bf1)
      type(hobf)::ho_bf1
      integer(c_int) :: ready

      HODLR_gpu_unsym_ready = .false.
      if (.not. c_associated(ho_bf1%gpu)) return
      call c_bpack_hodlr_gpu_unsym_ready(ho_bf1%gpu, ready)
      HODLR_gpu_unsym_ready = (ready == 1)
   end function HODLR_gpu_unsym_ready

   !> Vout = op(A)^{-1} Vin with the unsymmetric factors on the GPU, as HODLR_Inv_Apply
   subroutine HODLR_gpu_Inv_Apply_unsym(trans, Ns, num_vectors, Vin, Vout, ho_bf1, mode, flops)
      character trans
      integer Ns, num_vectors, mode
      DT, target::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hobf)::ho_bf1
      real(kind=8) :: flops  ! (out: this rank's)
      character(kind=c_char) :: op
      integer(c_int) :: nrhs, mode_c
      real(c_double) :: flops_c

      nrhs = num_vectors
      mode_c = mode
      op = 'N'
      if (trans /= 'N') op = 'T'
#if DAT==0 || DAT==2
      if (trans == 'C') Vin = conjg(Vin)
#endif
      call c_bpack_hodlr_gpu_solve_unsym(ho_bf1%gpu, op, nrhs, c_loc(Vin(1, 1)), c_loc(Vout(1, 1)), mode_c, flops_c)
      flops = flops_c
#if DAT==0 || DAT==2
      if (trans == 'C') then
         Vout = conjg(Vout)
         Vin = conjg(Vin)
      endif
#endif
   end subroutine HODLR_gpu_Inv_Apply_unsym

   !> Compare the unsymmetric factors on the GPU with those the CPU left in
   !> ho_bf1 (BPACK_CHECK=hodlr): per level, the largest relative differences
   !> of the updated blocks (BP_inverse_update) and of the Schur corrections
   !> (BP_inverse_schur); the leaf inverses (BP_inverse)
   subroutine HODLR_gpu_compare_unsymfactor(ho_bf1, option, ptree)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      DT, allocatable, target::buf(:, :)
      type(matrixblock), pointer::blk
      integer level, ii, idx_start_glo, ierr
      logical shared
      integer(c_int) :: what, level_c, m, k, found
      integer(c_int64_t) :: row0, col0, n_loc
      real(kind=8) :: err(2), d

      call HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      err = 0
      level = ho_bf1%Maxlevel + 1
      do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
         blk => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
         if (.not. associated(blk%fullmat)) cycle
         m = blk%M
         allocate (buf(m, m))
         what = 2
         level_c = level
         row0 = blk%headm - idx_start_glo
         col0 = row0
         call c_bpack_hodlr_gpu_download_unsym(ho_bf1%gpu, what, level_c, row0, col0, m, m, c_loc(buf), found)
         call assert(found == 1, 'HODLR_gpu_compare_unsymfactor: leaf not found on the GPU')
         d = fnorm(buf - blk%fullmat, m, m)/max(fnorm(blk%fullmat, m, m), tiny(1d0))
         err(1) = max(err(1), d)
         deallocate (buf)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A,Es10.2)') ' HODLR GPU check (unsym factor) leaf inverses:', err(1)
      do level = ho_bf1%Maxlevel, 1, -1
         err = 0
         ! (the CPU keeps the updated blocks of a shared level in its node's
         ! layout, not in the device's: those levels are compared through the
         ! log-determinant and the solves)
         shared = .false.
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            if (ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%nproc > 1) shared = .true.
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, shared, 1, MPI_LOGICAL, MPI_LOR, ptree%Comm, ierr)
         if (shared) then
            if (ptree%MyID == Main_ID) write (*, '(A,I3,A)') ' HODLR GPU check (unsym factor) level', level, &
               ': shared by several ranks (checked through the log-determinant and the solves)'
            cycle
         endif
         do ii = ho_bf1%levels(level)%Bidxs*2 - 1, ho_bf1%levels(level)%Bidxe*2
            blk => ho_bf1%levels(level)%BP_inverse_update(ii)%LL(1)%matrices_block(1)
            m = blk%M_loc
            k = size(blk%ButterflyU%blocks(1)%matrix, 2)
            allocate (buf(max(m, 1), max(k, 1)))
            what = 0
            level_c = level
            row0 = blk%headm - idx_start_glo
            col0 = blk%headn - idx_start_glo
            call c_bpack_hodlr_gpu_download_unsym(ho_bf1%gpu, what, level_c, row0, col0, m, k, c_loc(buf), found)
            call assert(found == 1, 'HODLR_gpu_compare_unsymfactor: block not found on the GPU')
            d = fnorm(buf(1:m, 1:k) - blk%ButterflyU%blocks(1)%matrix, m, k) &
               /max(fnorm(blk%ButterflyU%blocks(1)%matrix, m, k), tiny(1d0))
            err(1) = max(err(1), d)
            deallocate (buf)
         enddo
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            blk => ho_bf1%levels(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)
            m = blk%M_loc
            k = size(blk%ButterflyU%blocks(1)%matrix, 2)
            allocate (buf(max(m, 1), max(k, 1)))
            what = 1
            level_c = level
            row0 = ho_bf1%levels(level)%BP_inverse_update(2*ii - 1)%LL(1)%matrices_block(1)%headm - idx_start_glo
            col0 = row0
            call c_bpack_hodlr_gpu_download_unsym(ho_bf1%gpu, what, level_c, row0, col0, m, k, c_loc(buf), found)
            call assert(found == 1, 'HODLR_gpu_compare_unsymfactor: node not found on the GPU')
            d = fnorm(buf(1:m, 1:k) - blk%ButterflyU%blocks(1)%matrix, m, k) &
               /max(fnorm(blk%ButterflyU%blocks(1)%matrix, m, k), tiny(1d0))
            err(2) = max(err(2), d)
            deallocate (buf)
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, err, 2, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
         if (ptree%MyID == Main_ID) write (*, '(A,I3,A,Es10.2,A,Es10.2)') ' HODLR GPU check (unsym factor) level', level, &
            ': updated blocks', err(1), ', Schur corrections', err(2)
      enddo
   end subroutine HODLR_gpu_compare_unsymfactor

   !> Print the relative difference of a GPU result to the CPU one (BPACK_CHECK=hodlr)
   subroutine HODLR_gpu_report(what, trans, Ns, num_vectors, Vgpu, Vref, ptree)
      character(len=*) what
      character trans
      integer Ns, num_vectors, ierr
      DT::Vgpu(Ns, num_vectors), Vref(Ns, num_vectors)
      type(proctree)::ptree
      real(kind=8)::norms(2)

      norms(1) = fnorm(Vgpu - Vref, Ns, num_vectors)**2
      norms(2) = fnorm(Vref, Ns, num_vectors)**2
      call MPI_ALLREDUCE(MPI_IN_PLACE, norms, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A,A,A,A,A,I6,A,Es11.3)') ' HODLR GPU check (', what, ' ', trans, &
         ', nrhs', num_vectors, '): relative difference to the CPU', sqrt(norms(1)/max(norms(2), tiny(norms(2))))
   end subroutine HODLR_gpu_report

   !> Give the GPU backend the points of the matrix in tree order (before
   !> msh%xyz is freed), for the HODLR construction with a GPU evaluator: their
   !> 0-based original indices, and their coordinates when msh%xyz has them.
   subroutine HODLR_gpu_set_points(bmat, option, msh)
      type(Bmatrix)::bmat
      type(Hoption)::option
      type(mesh)::msh
      real(c_double), allocatable :: xyz(:, :)
      integer(c_int64_t), allocatable :: ids(:)
      integer(c_int64_t) :: n
      integer(c_int) :: dim
      integer i

      if (option%HODLR_use_gpu <= 0 .or. option%format /= HODLR .or. .not. c_associated(bmat%gpu)) return
      if (.not. allocated(msh%new2old)) return
      dim = 0
      if (allocated(msh%xyz)) then
         if (lbound(msh%xyz, 2) <= 1 .and. ubound(msh%xyz, 2) >= msh%Nunk) dim = size(msh%xyz, 1)
      endif
      n = msh%Nunk
      allocate (xyz(max(1, dim), max(1, msh%Nunk)), ids(max(1, msh%Nunk)))
      xyz = 0
      do i = 1, msh%Nunk
         if (dim > 0) xyz(:, i) = msh%xyz(:, msh%new2old(i))
         ids(i) = msh%new2old(i) - 1
      enddo
      call c_bpack_hodlr_gpu_set_points(bmat%gpu, n, dim, xyz, ids)
      deallocate (xyz, ids)
   end subroutine HODLR_gpu_set_points

   !> Rows of device matrix src (ncols columns; its rows lay_i over group
   !> pgno_i, 0 on the ranks outside) into columns c0 .. c0 + ncols - 1 of
   !> device matrix dst (rows lay_o over group pgno_o, 0 outside): the plan
   !> of Redistribute1Dto1D, the data between the GPUs
   !> (DistQr::dm_redistribute)
   subroutine HODLR_gpu_dm_redist(gpu, src, lay_i, pgno_i, dst, c0, lay_o, pgno_o, ncols, ptree)
      implicit none
      type(c_ptr) :: gpu
      integer(c_int) :: src, dst
      integer c0, pgno_i, pgno_o, ncols
      integer :: lay_i(:, :), lay_o(:, :)
      type(proctree)::ptree
      integer(c_int), allocatable :: srank(:), soff(:), srows(:), rrank(:), roff(:), rrows(:)
      integer(c_int) :: src_c, dst_c, ncols_c, c0_c, ns_c, nr_c, tag_c
      integer nproc_i, nproc_o, head_i, head_o, ii, jj, is, ie, os, oe, lo, hi, ns, nr

      if (ncols <= 0) return
      nproc_i = ptree%pgrp(pgno_i)%nproc
      nproc_o = ptree%pgrp(pgno_o)%nproc
      head_i = ptree%pgrp(pgno_i)%head
      head_o = ptree%pgrp(pgno_o)%head
      allocate (srank(max(1, nproc_o)), soff(max(1, nproc_o)), srows(max(1, nproc_o)))
      allocate (rrank(max(1, nproc_i)), roff(max(1, nproc_i)), rrows(max(1, nproc_i)))
      ns = 0
      nr = 0
      src_c = 0
      dst_c = 0
      if (IOwnPgrp(ptree, pgno_i)) then
         src_c = src
         ii = ptree%MyID - head_i + 1
         is = lay_i(ii, 1)
         ie = lay_i(ii, 2)
         do jj = 1, nproc_o
            lo = max(is, lay_o(jj, 1))
            hi = min(ie, lay_o(jj, 2))
            if (hi >= lo) then
               ns = ns + 1
               srank(ns) = jj + head_o - 1
               soff(ns) = lo - is
               srows(ns) = hi - lo + 1
            endif
         enddo
      endif
      if (IOwnPgrp(ptree, pgno_o)) then
         dst_c = dst
         jj = ptree%MyID - head_o + 1
         os = lay_o(jj, 1)
         oe = lay_o(jj, 2)
         do ii = 1, nproc_i
            lo = max(os, lay_i(ii, 1))
            hi = min(oe, lay_i(ii, 2))
            if (hi >= lo) then
               nr = nr + 1
               rrank(nr) = ii + head_i - 1
               roff(nr) = lo - os
               rrows(nr) = hi - lo + 1
            endif
         enddo
      endif
      ncols_c = ncols
      c0_c = c0
      ns_c = ns
      nr_c = nr
      tag_c = pgno_o
      call c_bpack_gpu_dm_redistribute(gpu, ptree%Comm, src_c, ncols_c, dst_c, c0_c, ns_c, srank, soff, srows, &
                                       nr_c, rrank, roff, rrows, tag_c)
      deallocate (srank, soff, srows, rrank, roff, rrows)
   end subroutine HODLR_gpu_dm_redist

   !> tab(level)%i(ii): the device index (from 0, the order of
   !> HODLR_gpu_upload_forward) of the forward block BP(ii) of level on this
   !> rank, -1 without one
   subroutine HODLR_gpu_block_table(ho_bf1, option, ptree, tab)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      type(hodlr_gpu_itab), allocatable :: tab(:)
      integer level, nd, ii, idx

      allocate (tab(max(1, ho_bf1%Maxlevel)))
      idx = 0
      do level = 1, ho_bf1%Maxlevel
         allocate (tab(level)%i(2*ho_bf1%levels(level)%Bidxs - 1:2*ho_bf1%levels(level)%Bidxe))
         tab(level)%i = -1
         do nd = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            if (.not. IOwnPgrp(ptree, ho_bf1%levels(level)%BP_inverse(nd)%pgno)) cycle
            do ii = nd*2 - 1, nd*2
               if (option%sym > 0 .and. mod(ii, 2) == 1) cycle
               tab(level)%i(ii) = idx
               idx = idx + 1
            enddo
         enddo
      enddo
   end subroutine HODLR_gpu_block_table

   !> The rows of the factors that BPACK_ExtractElement reads, for a HODLR
   !> built on the GPU whose host copies are unfilled (HODLR_gpu_defer_host):
   !> for each low-rank block of lstblk (this rank's blocks of the
   !> intersections), the rows of U of its requested rows and the rows of V
   !> of its requested columns in this rank's parts (M_p, N_p), into the
   !> block's gpu_u, gpu_v (gpu_urow, gpu_vrow: where each local row is), which
   !> LR_block_extraction and LR_all2all_extraction read instead of the host
   !> copies (writing scattered rows into those would touch most of their
   !> pages).  A row of U or V of a block is on the device of the rank that
   !> owns that row of the matrix (a shared block's V in its column child's
   !> layout; the symmetric A12 = A21^T reads A21's V and U), so each is
   !> gathered there and sent here.  Collective over ptree%Comm
   subroutine HODLR_gpu_fill_rows(ho_bf1, lstblk, inters, option, ptree)
      type(hobf), target::ho_bf1
      type(list)::lstblk
      type(intersect)::inters(:)
      type(Hoption)::option
      type(proctree)::ptree
      type(matrixblock), pointer::blocks
      type(nod), pointer::cur, curr, curc
      class(*), pointer::ptrr, ptrc
      type(block_ptr), allocatable :: blist(:)
      type(hodlr_gpu_itab), allocatable :: tab(:)
      integer, allocatable :: q_owner(:), q_int(:, :), q_blk(:), q_hrole(:), q_row(:)
      integer, allocatable :: rmap(:)
      DT, allocatable :: vals(:)
      integer(kind=8), allocatable :: voff(:)
      integer nblk, nq, nqmax, b, nn, t, g, r, lc, ii, iis, role, hrole, pp, k
      integer ncache
      integer nloc
      logical view

      ! this rank's low-rank blocks and an upper bound of the rows they need
      allocate (blist(max(1, lstblk%num_nods)))
      nblk = 0
      nqmax = 0
      cur => lstblk%head
      do b = 1, lstblk%num_nods
         select type (ptr => cur%item)
         type is (block_ptr)
            if (ptr%ptr%style /= 1) then
               call assert(ptr%ptr%level_butterfly == 0, 'HODLR_gpu_fill_rows: a block is not low rank')
               nblk = nblk + 1
               blist(nblk)%ptr => ptr%ptr
               curr => ptr%ptr%lstr%head
               curc => ptr%ptr%lstc%head
               do nn = 1, ptr%ptr%lstr%num_nods
                  ptrr => curr%item
                  ptrc => curc%item
                  select type (ptrr)
                  type is (iarray)
                     nqmax = nqmax + ptrr%num_nods
                  end select
                  select type (ptrc)
                  type is (iarray)
                     nqmax = nqmax + ptrc%num_nods
                  end select
                  curr => curr%next
                  curc => curc%next
               enddo
            endif
         end select
         cur => cur%next
      enddo
      ! the requests, one per distinct host row: (level, stored block, role on the device, row of the matrix, rank)
      allocate (q_owner(max(1, nqmax)), q_int(5, max(1, nqmax)), q_blk(max(1, nqmax)), q_hrole(max(1, nqmax)), &
                q_row(max(1, nqmax)))
      nq = 0
      do b = 1, nblk
         blocks => blist(b)%ptr
         lc = GetTreelevel(blocks%row_group) - 1
         ii = blocks%row_group - 2**lc + 1
         call assert(lc >= 1 .and. lc <= ho_bf1%Maxlevel, 'HODLR_gpu_fill_rows: a block outside the low-rank levels')
         call assert(ho_bf1%levels(lc)%BP(ii)%LL(1)%matrices_block(1)%row_group == blocks%row_group, &
                     'HODLR_gpu_fill_rows: a block not where its row group puts it')
         view = blocks%is_transpose_view == 1
         iis = ii
         if (view) iis = ii + 1  ! (A12 = A21^T: its U is A21's V, its V A21's U)
         pp = ptree%MyID - ptree%pgrp(blocks%pgno)%head + 1
         k = size(blocks%ButterflyU%blocks(1)%matrix, 2)
         do hrole = 0, 1  ! the host rows of U (of the requested rows), then of V (of the requested columns)
            if (hrole == 0) then
               nloc = blocks%M_p(pp, 2) - blocks%M_p(pp, 1) + 1
            else
               nloc = blocks%N_p(pp, 2) - blocks%N_p(pp, 1) + 1
            endif
            if (nloc <= 0) cycle
            role = hrole
            if (view) role = 1 - hrole
            allocate (rmap(nloc))
            rmap = 0
            ncache = 0
            curr => blocks%lstr%head
            curc => blocks%lstc%head
            do nn = 1, blocks%lstr%num_nods
               if (hrole == 0) then
                  ptrr => curr%item
               else
                  ptrr => curc%item
               endif
               select type (ptrr)
               type is (iarray)
                  do t = 1, ptrr%num_nods
                     if (hrole == 0) then
                        g = inters(ptrr%idx)%rows(ptrr%dat(t))
                        r = g - blocks%headm + 1 - blocks%M_p(pp, 1) + 1
                     else
                        g = inters(ptrr%idx)%cols(ptrr%dat(t))
                        r = g - blocks%headn + 1 - blocks%N_p(pp, 1) + 1
                     endif
                     if (r < 1 .or. r > nloc) cycle
                     if (rmap(r) > 0) cycle
                     ncache = ncache + 1
                     rmap(r) = ncache
                     nq = nq + 1
                     q_owner(nq) = HODLR_gpu_row_owner(ho_bf1, ptree, g)
                     q_int(:, nq) = (/lc, iis, role, g, k/)
                     q_blk(nq) = b
                     q_hrole(nq) = hrole
                     q_row(nq) = ncache
                  enddo
               end select
               curr => curr%next
               curc => curc%next
            enddo
            ! (the block's cache of these rows)
            if (hrole == 0) then
               if (allocated(blocks%gpu_urow)) deallocate (blocks%gpu_urow)
               if (allocated(blocks%gpu_u)) deallocate (blocks%gpu_u)
               call move_alloc(rmap, blocks%gpu_urow)
               allocate (blocks%gpu_u(max(1, ncache), max(1, k)))
            else
               if (allocated(blocks%gpu_vrow)) deallocate (blocks%gpu_vrow)
               if (allocated(blocks%gpu_v)) deallocate (blocks%gpu_v)
               call move_alloc(rmap, blocks%gpu_vrow)
               allocate (blocks%gpu_v(max(1, ncache), max(1, k)))
            endif
         enddo
      enddo
      ! the rows from the devices of their owners
      call HODLR_gpu_fetch_rows(ho_bf1, option, ptree, nq, q_owner, q_int, vals, voff)
      do t = 1, nq
         blocks => blist(q_blk(t))%ptr
         k = q_int(5, t)
         if (q_hrole(t) == 0) then
            blocks%gpu_u(q_row(t), 1:k) = vals(voff(t) + 1:voff(t) + k)
         else
            blocks%gpu_v(q_row(t), 1:k) = vals(voff(t) + 1:voff(t) + k)
         endif
      enddo
      deallocate (blist, q_owner, q_int, q_blk, q_hrole, q_row, vals, voff)
   end subroutine HODLR_gpu_fill_rows

   !> The rank that owns row g of the matrix (the layout of the vectors, where
   !> the rows of the forward blocks are on the devices)
   integer function HODLR_gpu_row_owner(ho_bf1, ptree, g)
      type(hobf), target::ho_bf1
      type(proctree)::ptree
      integer g
      type(matrixblock), pointer::root
      integer lo, hi, mid

      root => ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)
      lo = 1
      hi = size(root%N_p, 1)
      do while (lo < hi)
         mid = (lo + hi)/2
         if (root%N_p(mid, 2) < g - root%headn + 1) then
            lo = mid + 1
         else
            hi = mid
         endif
      enddo
      HODLR_gpu_row_owner = ptree%pgrp(root%pgno)%head + lo - 1
   end function HODLR_gpu_row_owner

   !> Rows of the forward blocks on the GPUs, from the ranks that hold them:
   !> request t (q_int(:, t): the level, the block in the level, its role on
   !> the device (0: U, 1: V), the row of the matrix, the block's rank k) to
   !> q_owner(t), the rank that owns that row (HODLR_gpu_row_owner), which
   !> gathers it from its device; vals(voff(t) + 1:voff(t) + k) its k values.
   !> Collective over ptree%Comm
   subroutine HODLR_gpu_fetch_rows(ho_bf1, option, ptree, nq, q_owner, q_int, vals, voff)
      type(hobf), target::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      integer nq
      integer :: q_owner(:), q_int(:, :)
      DT, allocatable :: vals(:)
      integer(kind=8), allocatable :: voff(:)
      type(hodlr_gpu_itab), allocatable :: tab(:)
      integer, allocatable :: order(:)
      integer, allocatable :: scount(:), rcount(:), sdisp(:), rdisp(:), svcount(:), rvcount(:), svdisp(:), rvdisp(:)
      integer, allocatable :: r_int(:, :), sbuf(:), rbuf(:)
      integer(c_int), allocatable :: grows(:)
      DT, allocatable, target :: vsend(:), vrecv(:), gout(:)
      integer t, lc, iis, nproc, ierr, q, k, nrecv
      integer idx_start_glo, run0, run1, pos, c, m
      integer(kind=8) :: ntot
      integer(c_int) :: idx_c, role_c, n_c, k_c
      integer(c_int64_t) :: n_loc

      nproc = ptree%nproc
      ! the requests to their owners (in the order of the owners, stable)
      allocate (scount(nproc), rcount(nproc), sdisp(nproc), rdisp(nproc), order(max(1, nq)))
      scount = 0
      do t = 1, nq
         scount(q_owner(t) + 1) = scount(q_owner(t) + 1) + 1
      enddo
      sdisp(1) = 0
      do q = 2, nproc
         sdisp(q) = sdisp(q - 1) + scount(q - 1)
      enddo
      rdisp = sdisp
      do t = 1, nq
         rdisp(q_owner(t) + 1) = rdisp(q_owner(t) + 1) + 1
         order(rdisp(q_owner(t) + 1)) = t
      enddo
      call MPI_ALLTOALL(scount, 1, MPI_INTEGER, rcount, 1, MPI_INTEGER, ptree%Comm, ierr)
      rdisp(1) = 0
      do q = 2, nproc
         rdisp(q) = rdisp(q - 1) + rcount(q - 1)
      enddo
      nrecv = sum(rcount)
      allocate (sbuf(5*max(1, nq)), rbuf(5*max(1, nrecv)))
      do t = 1, nq
         sbuf(5*t - 4:5*t) = q_int(:, order(t))
      enddo
      call MPI_ALLTOALLV(sbuf, 5*scount, 5*sdisp, MPI_INTEGER, rbuf, 5*rcount, 5*rdisp, MPI_INTEGER, ptree%Comm, ierr)
      allocate (r_int(5, max(1, nrecv)))
      r_int(:, 1:nrecv) = reshape(rbuf(1:5*nrecv), (/5, nrecv/))
      ! the owner: the requested rows from its device, a run of requests of one block and role at a time
      call HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      call HODLR_gpu_block_table(ho_bf1, option, ptree, tab)
      allocate (svcount(nproc), rvcount(nproc), svdisp(nproc), rvdisp(nproc))
      rvcount = 0  ! (the values this rank answers to each requester)
      do q = 1, nproc
         do t = rdisp(q) + 1, rdisp(q) + rcount(q)
            rvcount(q) = rvcount(q) + r_int(5, t)
         enddo
      enddo
      rvdisp(1) = 0
      do q = 2, nproc
         rvdisp(q) = rvdisp(q - 1) + rvcount(q - 1)
      enddo
      allocate (vsend(max(1, sum(rvcount))))
      run0 = 1
      pos = 0
      do while (run0 <= nrecv)
         run1 = run0
         do while (run1 < nrecv)
            if (any(r_int(1:3, run1 + 1) /= r_int(1:3, run0))) exit
            run1 = run1 + 1
         enddo
         lc = r_int(1, run0)
         iis = r_int(2, run0)
         call assert(lc >= 1 .and. lc <= ho_bf1%Maxlevel, 'HODLR_gpu_fetch_rows: a request outside the low-rank levels')
         call assert(iis >= lbound(tab(lc)%i, 1) .and. iis <= ubound(tab(lc)%i, 1), &
                     'HODLR_gpu_fetch_rows: a request for a block of another rank')
         idx_c = tab(lc)%i(iis)
         call assert(idx_c >= 0, 'HODLR_gpu_fetch_rows: a request for a block not on this device')
         role_c = r_int(3, run0)
         n_c = run1 - run0 + 1
         allocate (grows(n_c))
         do t = run0, run1
            grows(t - run0 + 1) = r_int(4, t) - idx_start_glo
         enddo
         k = r_int(5, run0)
         allocate (gout(max(1, n_c*k)))
         call c_bpack_hodlr_gpu_gather_rows(ho_bf1%gpu, idx_c, role_c, n_c, grows, c_loc(gout(1)), k_c)
         call assert(k_c == k, 'HODLR_gpu_fetch_rows: a block of another rank on the host than on the device')
         do t = 1, n_c  ! (row t of gout, n_c x k: the k values of request run0 + t - 1)
            do c = 1, k
               vsend(pos + c) = gout(t + (c - 1)*n_c)
            enddo
            pos = pos + k
         enddo
         deallocate (grows, gout)
         run0 = run1 + 1
      enddo
      ! the rows back to the requesters
      svcount = 0
      do t = 1, nq
         svcount(q_owner(t) + 1) = svcount(q_owner(t) + 1) + q_int(5, t)
      enddo
      svdisp(1) = 0
      do q = 2, nproc
         svdisp(q) = svdisp(q - 1) + svcount(q - 1)
      enddo
      allocate (vrecv(max(1, sum(svcount))))
      call MPI_ALLTOALLV(vsend, rvcount, rvdisp, MPI_DT, vrecv, svcount, svdisp, MPI_DT, ptree%Comm, ierr)
      ! (in the order of the requests)
      allocate (voff(max(1, nq)))
      ntot = 0
      do t = 1, nq
         voff(t) = ntot
         ntot = ntot + q_int(5, t)
      enddo
      allocate (vals(max(1_8, ntot)))
      pos = 0
      do t = 1, nq
         m = order(t)
         k = q_int(5, m)
         vals(voff(m) + 1:voff(m) + k) = vrecv(pos + 1:pos + k)
         pos = pos + k
      enddo
      deallocate (order, scount, rcount, sdisp, rdisp, sbuf, rbuf, r_int)
      deallocate (svcount, rvcount, svdisp, rvdisp, vsend, vrecv, tab)
   end subroutine HODLR_gpu_fetch_rows

   !> BPACK_ExtractElement of a HODLR built on the GPU whose host copies of the
   !> factors are unfilled (HODLR_gpu_defer_host): the products where the
   !> factors are.  lstblk
   !> holds the blocks of the intersections of every node this rank belongs to
   !> (HODLR_MapIntersec2Block with members), with the positions of their rows
   !> and columns in the intersections.  The piece of an intersection in a
   !> low-rank block U V^T, rows R and columns C, is computed by the ranks
   !> holding the rows of its larger side on their devices: with R >= C the
   !> owners of the rows of R, U(R_d, :) V(C, :)^T; else the owners of C,
   !> U(R, :) V(C_d, :)^T.  Only the rows of the smaller side move (from their
   !> owners' devices, HODLR_gpu_fetch_rows), the products run on the devices
   !> (c_bpack_hodlr_gpu_extract_products), and the entries go to the ranks
   !> holding them: the intersection's rank, or for a 2D grid (pmaps) each rank
   !> its block-cyclic part.  A piece in a dense leaf is computed on the host by
   !> the leaf's rank.  Collective over ptree%Comm
   subroutine HODLR_gpu_extract_blocks(ho_bf1, lstblk, inters, Npmap, pmaps, option, stats, ptree)
      type(hobf), target::ho_bf1
      type(list)::lstblk
      type(intersect)::inters(:)
      integer Npmap, pmaps(Npmap, 3)
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(matrixblock), pointer::blocks
      type(nod), pointer::cur, curr, curc
      class(*), pointer::ptrr, ptrc
      type(hodlr_gpu_itab), allocatable :: tab(:)
      ! the pieces this rank computes, the low-rank ones (on the device) first, then the leaves'
      integer, allocatable :: p_nn(:), p_nrow(:), p_ncol(:), p_q0(:), pos(:)
      integer(kind=8), allocatable :: p_pos(:), p_val(:)
      integer(c_int), allocatable :: d_idx(:), d_role(:), d_nbig(:), d_nsmall(:), d_first(:), d_rows(:)
      integer, allocatable :: q_owner(:), q_int(:, :)
      integer, allocatable :: scount(:), rcount(:), sdisp(:), rdisp(:), svcount(:), rvcount(:), svdisp(:), rvdisp(:)
      integer, allocatable :: sbuf(:), rbuf(:)
      DT, allocatable, target :: small(:), outv(:), vals(:), svals(:), rvals(:)
      integer(kind=8), allocatable :: voff(:), ipos_d(:), ival_d(:)
      integer pass, b, it, nn, nR, nC, nb, ns, t, s, c, g, lc, ii, iis, urole, k, nlr, nlf, ilr, ilf, iq, nq, p, np
      integer dest, nproc, ierr, i, j, idx_start_glo, ge  ! (nr, nc below: nR, nC)
      integer pass2, npr, npc, gh, pr, pc, nrd, ncd, rr, cc, lr_, lc_
      integer, allocatable :: rsel(:), csel(:)
      integer(kind=8) :: npos, nbig_tot, nout, nout_lr, nsmallv, ipos, ibig, ival, ism, nis, nvs, ip, iv
      integer(c_int) :: idx_c, k_c, np_c
      integer(c_int64_t) :: n_loc
      real(c_double) :: flops
      logical view, rowmode, leaf

      nproc = ptree%nproc
      call HODLR_gpu_local_rows(ho_bf1, ptree, idx_start_glo, n_loc)
      ge = idx_start_glo + int(n_loc) - 1
      call HODLR_gpu_block_table(ho_bf1, option, ptree, tab)
      nlr = 0
      nlf = 0
      nq = 0
      npos = 0
      nbig_tot = 0
      nout_lr = 0
      nout = 0
      nsmallv = 0
      do pass = 1, 2  ! (count, then record)
         if (pass == 2) then
            np = nlr + nlf
            allocate (p_nn(max(1, np)), p_nrow(max(1, np)), p_ncol(max(1, np)), p_q0(max(1, np)), p_pos(max(1, np)), &
                      p_val(max(1, np)), pos(max(1_8, npos)))
            allocate (d_idx(max(1, nlr)), d_role(max(1, nlr)), d_nbig(max(1, nlr)), d_nsmall(max(1, nlr)), &
                      d_first(max(1, nlr)), d_rows(max(1_8, nbig_tot)))
            allocate (q_owner(max(1, nq)), q_int(5, max(1, nq)), outv(max(1_8, nout)), small(max(1_8, nsmallv)))
         endif
         ilr = 0
         ilf = nlr
         iq = 0
         ipos = 0
         ibig = 0
         ival = nout_lr  ! (the leaves' values after the products')
         cur => lstblk%head
         do b = 1, lstblk%num_nods
            select type (ptr => cur%item)
            type is (block_ptr)
               blocks => ptr%ptr
               leaf = blocks%style == 1
               if (.not. leaf) then
                  call assert(blocks%level_butterfly == 0, 'HODLR_gpu_extract_blocks: a block is not low rank')
                  lc = GetTreelevel(blocks%row_group) - 1
                  ii = blocks%row_group - 2**lc + 1
                  view = blocks%is_transpose_view == 1
                  iis = ii
                  if (view) iis = ii + 1  ! (A12 = A21^T: its U is A21's V, its V A21's U)
                  idx_c = tab(lc)%i(iis)
                  call assert(idx_c >= 0, 'HODLR_gpu_extract_blocks: a block not on this device')
                  call c_bpack_hodlr_gpu_block_rank(ho_bf1%gpu, idx_c, k_c)
                  k = k_c
                  urole = 0
                  if (view) urole = 1
               endif
               curr => blocks%lstr%head
               curc => blocks%lstc%head
               do it = 1, blocks%lstr%num_nods
                  ptrr => curr%item
                  ptrc => curc%item
                  select type (ptrr)
                  type is (iarray)
                     select type (ptrc)
                     type is (iarray)
                        nn = ptrr%idx
                        nR = ptrr%num_nods
                        nC = ptrc%num_nods
                        if (leaf) then
                           if (ptree%MyID == ptree%pgrp(blocks%pgno)%head .and. nR > 0 .and. nC > 0) then
                              if (pass == 1) then
                                 nlf = nlf + 1
                                 npos = npos + nR + nC
                                 nout = nout + int(nR, 8)*nC
                              else
                                 ilf = ilf + 1
                                 p_nn(ilf) = nn
                                 p_nrow(ilf) = nR
                                 p_ncol(ilf) = nC
                                 p_pos(ilf) = ipos
                                 pos(ipos + 1:ipos + nR) = ptrr%dat(1:nR)
                                 pos(ipos + nR + 1:ipos + nR + nC) = ptrc%dat(1:nC)
                                 ipos = ipos + nR + nC
                                 p_val(ilf) = ival
                                 do j = 1, nC
                                    c = inters(nn)%cols(ptrc%dat(j)) - blocks%headn + 1
                                    do i = 1, nR
                                       outv(ival + i + int(j - 1, 8)*nR) = &
                                          blocks%fullmat(inters(nn)%rows(ptrr%dat(i)) - blocks%headm + 1, c)
                                    enddo
                                 enddo
                                 ival = ival + int(nR, 8)*nC
                              endif
                           endif
                        else
                           ! the rows of the larger side on this rank's device
                           rowmode = nR >= nC
                           nb = 0
                           if (rowmode) then
                              do t = 1, nR
                                 g = inters(nn)%rows(ptrr%dat(t))
                                 if (g >= idx_start_glo .and. g <= ge) nb = nb + 1
                              enddo
                              ns = nC
                           else
                              do t = 1, nC
                                 g = inters(nn)%cols(ptrc%dat(t))
                                 if (g >= idx_start_glo .and. g <= ge) nb = nb + 1
                              enddo
                              ns = nR
                           endif
                           if (nb > 0 .and. ns > 0) then
                              if (pass == 1) then
                                 nlr = nlr + 1
                                 npos = npos + nb + ns
                                 nbig_tot = nbig_tot + nb
                                 nq = nq + ns
                                 nout_lr = nout_lr + int(nb, 8)*ns
                                 nout = nout + int(nb, 8)*ns
                                 nsmallv = nsmallv + int(ns, 8)*k
                              else
                                 ilr = ilr + 1
                                 p_nn(ilr) = nn
                                 p_pos(ilr) = ipos
                                 p_q0(ilr) = iq
                                 d_idx(ilr) = idx_c
                                 d_nbig(ilr) = nb
                                 d_nsmall(ilr) = ns
                                 if (rowmode) then  ! U(R_d, :) V(C, :)^T: rows R_d, columns C
                                    p_nrow(ilr) = nb
                                    p_ncol(ilr) = ns
                                    d_role(ilr) = urole
                                    d_first(ilr) = 0
                                    do t = 1, nR
                                       g = inters(nn)%rows(ptrr%dat(t))
                                       if (g < idx_start_glo .or. g > ge) cycle
                                       ipos = ipos + 1
                                       pos(ipos) = ptrr%dat(t)
                                       ibig = ibig + 1
                                       d_rows(ibig) = g - idx_start_glo
                                    enddo
                                    do t = 1, nC
                                       g = inters(nn)%cols(ptrc%dat(t))
                                       ipos = ipos + 1
                                       pos(ipos) = ptrc%dat(t)
                                       iq = iq + 1
                                       q_owner(iq) = HODLR_gpu_row_owner(ho_bf1, ptree, g)
                                       q_int(:, iq) = (/lc, iis, 1 - urole, g, k/)
                                    enddo
                                 else  ! U(R, :) V(C_d, :)^T: rows R, columns C_d
                                    p_nrow(ilr) = ns
                                    p_ncol(ilr) = nb
                                    d_role(ilr) = 1 - urole
                                    d_first(ilr) = 1
                                    do t = 1, nR
                                       g = inters(nn)%rows(ptrr%dat(t))
                                       ipos = ipos + 1
                                       pos(ipos) = ptrr%dat(t)
                                       iq = iq + 1
                                       q_owner(iq) = HODLR_gpu_row_owner(ho_bf1, ptree, g)
                                       q_int(:, iq) = (/lc, iis, urole, g, k/)
                                    enddo
                                    do t = 1, nC
                                       g = inters(nn)%cols(ptrc%dat(t))
                                       if (g < idx_start_glo .or. g > ge) cycle
                                       ipos = ipos + 1
                                       pos(ipos) = ptrc%dat(t)
                                       ibig = ibig + 1
                                       d_rows(ibig) = g - idx_start_glo
                                    enddo
                                 endif
                              endif
                           endif
                        endif
                     end select
                  end select
                  curr => curr%next
                  curc => curc%next
               enddo
            end select
            cur => cur%next
         enddo
      enddo

      ! the smaller sides from the devices that hold them, then the products on this rank's device
      call HODLR_gpu_fetch_rows(ho_bf1, option, ptree, nq, q_owner, q_int, vals, voff)
      ism = 0
      ival = 0
      do p = 1, nlr
         ns = d_nsmall(p)
         k = q_int(5, p_q0(p) + 1)
         do s = 1, ns
            iq = p_q0(p) + s
            do c = 1, k
               small(ism + s + int(c - 1, 8)*ns) = vals(voff(iq) + c)
            enddo
         enddo
         ism = ism + int(ns, 8)*k
         p_val(p) = ival
         ival = ival + int(d_nbig(p), 8)*ns
      enddo
      flops = 0
      if (nlr > 0) then
         np_c = nlr
         call c_bpack_hodlr_gpu_extract_products(ho_bf1%gpu, np_c, d_idx, d_role, d_nbig, d_nsmall, d_first, d_rows, &
                                                 c_loc(small(1)), c_loc(outv(1)), flops)
      endif
      stats%Flop_Tmp = stats%Flop_Tmp + flops

      ! the entries to the ranks holding them: per piece and rank of its intersection's grid (pmaps: nprow x
      ! npcol, block cyclic by nbslpk from process (0, 0), column-major ranks from pmaps(:, 3) as Gridinfo_2D; a
      ! 1 x 1 grid takes the whole piece) its intersection, its numbers of rows and columns there and their
      ! positions in the intersection, then its values (column-major); the receiver puts them at their local
      ! positions
      allocate (scount(nproc), rcount(nproc), sdisp(nproc), rdisp(nproc), svcount(nproc), rvcount(nproc), &
                svdisp(nproc), rvdisp(nproc), ipos_d(nproc), ival_d(nproc))
      nrd = 1
      ncd = 1
      if (np > 0) nrd = max(1, maxval(p_nrow(1:np)))
      if (np > 0) ncd = max(1, maxval(p_ncol(1:np)))
      allocate (rsel(nrd), csel(ncd))
      scount = 0
      svcount = 0
      nis = 0
      nvs = 0
      do pass2 = 1, 2  ! (count, then fill)
         if (pass2 == 2) then
            call assert(nis < 2_8**31 .and. nvs < 2_8**31, 'HODLR_gpu_extract_blocks: more than 2^31 entries to send')
            sdisp(1) = 0
            svdisp(1) = 0
            do dest = 2, nproc
               sdisp(dest) = sdisp(dest - 1) + scount(dest - 1)
               svdisp(dest) = svdisp(dest - 1) + svcount(dest - 1)
            enddo
            allocate (sbuf(max(1_8, nis)), svals(max(1_8, nvs)))
            ipos_d = sdisp
            ival_d = svdisp
         endif
         do p = 1, np
            nr = p_nrow(p)
            nc = p_ncol(p)
            npr = pmaps(inters(p_nn(p))%pg, 1)
            npc = pmaps(inters(p_nn(p))%pg, 2)
            gh = pmaps(inters(p_nn(p))%pg, 3)
            do pc = 0, npc - 1
               ncd = 0
               do j = 1, nc
                  if (mod((pos(p_pos(p) + nr + j) - 1)/nbslpk, npc) /= pc) cycle
                  ncd = ncd + 1
                  csel(ncd) = j
               enddo
               if (ncd == 0) cycle
               do pr = 0, npr - 1
                  nrd = 0
                  do i = 1, nr
                     if (mod((pos(p_pos(p) + i) - 1)/nbslpk, npr) /= pr) cycle
                     nrd = nrd + 1
                     rsel(nrd) = i
                  enddo
                  if (nrd == 0) cycle
                  dest = gh + pr + pc*npr + 1
                  if (pass2 == 1) then
                     scount(dest) = scount(dest) + 3 + nrd + ncd
                     svcount(dest) = svcount(dest) + nrd*ncd
                     nis = nis + 3 + nrd + ncd
                     nvs = nvs + int(nrd, 8)*ncd
                  else
                     sbuf(ipos_d(dest) + 1:ipos_d(dest) + 3) = (/p_nn(p), nrd, ncd/)
                     do i = 1, nrd
                        sbuf(ipos_d(dest) + 3 + i) = pos(p_pos(p) + rsel(i))
                     enddo
                     do j = 1, ncd
                        sbuf(ipos_d(dest) + 3 + nrd + j) = pos(p_pos(p) + nr + csel(j))
                     enddo
                     ipos_d(dest) = ipos_d(dest) + 3 + nrd + ncd
                     do j = 1, ncd
                        do i = 1, nrd
                           svals(ival_d(dest) + i + int(j - 1, 8)*nrd) = outv(p_val(p) + rsel(i) + int(csel(j) - 1, 8)*nr)
                        enddo
                     enddo
                     ival_d(dest) = ival_d(dest) + int(nrd, 8)*ncd
                  endif
               enddo
            enddo
         enddo
      enddo
      call MPI_ALLTOALL(scount, 1, MPI_INTEGER, rcount, 1, MPI_INTEGER, ptree%Comm, ierr)
      call MPI_ALLTOALL(svcount, 1, MPI_INTEGER, rvcount, 1, MPI_INTEGER, ptree%Comm, ierr)
      rdisp(1) = 0
      rvdisp(1) = 0
      do dest = 2, nproc
         rdisp(dest) = rdisp(dest - 1) + rcount(dest - 1)
         rvdisp(dest) = rvdisp(dest - 1) + rvcount(dest - 1)
      enddo
      allocate (rbuf(max(1, sum(rcount))), rvals(max(1, sum(rvcount))))
      call MPI_ALLTOALLV(sbuf, scount, sdisp, MPI_INTEGER, rbuf, rcount, rdisp, MPI_INTEGER, ptree%Comm, ierr)
      call MPI_ALLTOALLV(svals, svcount, svdisp, MPI_DT, rvals, rvcount, rvdisp, MPI_DT, ptree%Comm, ierr)
      ip = 0
      iv = 0
      do while (ip < sum(rcount))
         nn = rbuf(ip + 1)
         nr = rbuf(ip + 2)
         nc = rbuf(ip + 3)
         call assert(associated(inters(nn)%dat_loc), 'HODLR_gpu_extract_blocks: entries of an intersection of another rank')
         npr = pmaps(inters(nn)%pg, 1)
         npc = pmaps(inters(nn)%pg, 2)
         do j = 1, nc
            cc = rbuf(ip + 3 + nr + j)
            lc_ = ((cc - 1)/(nbslpk*npc))*nbslpk + mod(cc - 1, nbslpk) + 1
            do i = 1, nr
               rr = rbuf(ip + 3 + i)
               lr_ = ((rr - 1)/(nbslpk*npr))*nbslpk + mod(rr - 1, nbslpk) + 1
               inters(nn)%dat_loc(lr_, lc_) = rvals(iv + i + int(j - 1, 8)*nr)
            enddo
         enddo
         ip = ip + 3 + nr + nc
         iv = iv + int(nr, 8)*nc
      enddo
      deallocate (rsel, csel)
      deallocate (p_nn, p_nrow, p_ncol, p_q0, p_pos, p_val, pos, d_idx, d_role, d_nbig, d_nsmall, d_first, d_rows)
      deallocate (q_owner, q_int, outv, small, vals, voff, tab)
      deallocate (scount, rcount, sdisp, rdisp, svcount, rvcount, svdisp, rvdisp, ipos_d, ival_d, sbuf, svals, rbuf, rvals)
   end subroutine HODLR_gpu_extract_blocks

   !> BPACK_CHECK=hodlr: the host factors of the low-rank blocks of this rank
   !> overwritten by a sentinel and marked unfilled, so that an extraction
   !> must take every row it reads from the GPUs (HODLR_gpu_fill_rows);
   !> HODLR_gpu_fetch_host restores them
   subroutine HODLR_gpu_poison_host(ho_bf1, option, ptree)
      type(hobf), target::ho_bf1
      type(Hoption)::option
      type(proctree)::ptree
      type(matrixblock), pointer::blk
      integer level, nd, ii

      do level = 1, ho_bf1%Maxlevel
         do nd = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            if (.not. IOwnPgrp(ptree, ho_bf1%levels(level)%BP_inverse(nd)%pgno)) cycle
            do ii = nd*2 - 1, nd*2
               if (option%sym > 0 .and. mod(ii, 2) == 1) cycle  ! (A12 is a view of A21)
               blk => ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1)
               if (.not. IOwnPgrp(ptree, blk%pgno)) cycle
               blk%ButterflyU%blocks(1)%matrix = 1d30
               blk%ButterflyV%blocks(1)%matrix = 1d30
            enddo
         enddo
      enddo
      ho_bf1%gpu_host_stale = .true.
   end subroutine HODLR_gpu_poison_host

   !> Stop when a CPU routine needs the host copies of the factors of a HODLR
   !> built on the GPU, which stay on the GPU (HODLR_gpu_defer_host)
   subroutine HODLR_gpu_require_host(ho_bf1, what)
      type(hobf)::ho_bf1
      character(len=*) :: what

      if (.not. ho_bf1%gpu_host_stale) return
      call BPACK_require(.false., 'HODLR GPU: '//what//' reads the host copies of the factors, which stay on the GPU. '// &
                  'TODO: not supported yet; BPACK_CHECK=hodlr or ErrFillFull=1 keep the host copies')
   end subroutine HODLR_gpu_require_host

   !> The host copies of the low-rank factors built on the GPU are left
   !> unfilled: the device keeps the only copy (HODLR_gpu_upload_forward
   !> moves it into the forward blocks), BPACK_ExtractElement takes the rows it
   !> reads from the GPUs (HODLR_gpu_fill_rows), and the CPU multiply and
   !> factorization stop (HODLR_gpu_require_host).  Not with BPACK_CHECK=hodlr
   !> (the CPU comparisons read them) or ErrFillFull=1, which download them
   !> during the construction.
   logical function HODLR_gpu_defer_host(option)
      type(Hoption)::option
      HODLR_gpu_defer_host = HODLR_gpu_check_level() <= 0 .and. option%ErrFillFull /= 1
   end function HODLR_gpu_defer_host

   !> The OpenMP threads of the host loops over the small blocks of the GPU
   !> construction (the BACA core QRs): OMP_NUM_THREADS, at most 8 (many
   !> threads making small LAPACK calls to an OpenBLAS built with
   !> USE_LOCKING mostly wait for its lock: 4 threads beat 64)
   integer function HODLR_gpu_host_threads()
#ifdef HAVE_OPENMP
      use omp_lib, only: omp_get_max_threads
      HODLR_gpu_host_threads = max(1, min(omp_get_max_threads(), 8))
#else
      HODLR_gpu_host_threads = 1
#endif
   end function HODLR_gpu_host_threads

   !> How the GPU construction splits a block compressed by several ranks
   !> into pieces.  1 (default): a block of the symmetric HODLR over its
   !> node's group (both children's ranks, the other half idle otherwise);
   !> with the option HODLR_gpu_pieces > 1 each rank's part split further and
   !> the pieces spread over the group by their estimated cost.  0
   !> (BPACK_CHECK=hodlr-exact): as the CPU (LR_HBACA_Leaflevel over the
   !> block's group, with the same random numbers, for comparisons).
   integer function HODLR_gpu_split_mode()
      integer, save :: mode = -1
      if (mode < 0) then
         mode = 1
         if (c_bpack_env_check('hodlr-exact'//c_null_char) /= 0) mode = 0
      endif
      HODLR_gpu_split_mode = mode
   end function HODLR_gpu_split_mode

   !> BPACK_CHECK=hodlr (1): compare each GPU result with the CPU routine it
   !> replaces and print the relative difference; hodlr-transpose (2): also
   !> the transposed products of each multiply (doc/environment_variables.md)
   integer function HODLR_gpu_check_level()
      integer, save :: mode = -1
      if (mode < 0) then
         mode = 0
         if (c_bpack_env_check('hodlr'//c_null_char) /= 0) mode = 1
         if (c_bpack_env_check('hodlr-transpose'//c_null_char) /= 0) mode = 2
      endif
      HODLR_gpu_check_level = mode
   end function HODLR_gpu_check_level

end module BPACK_GPU
