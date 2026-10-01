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
      subroutine c_bpack_hodlr_gpu_construct_ready(gpu, ready) bind(c, name="c_bpack_hodlr_gpu_construct_ready")
         import :: c_ptr, c_int
         type(c_ptr), value :: gpu
         integer(c_int) :: ready
      end subroutine c_bpack_hodlr_gpu_construct_ready

      subroutine c_bpack_hodlr_gpu_set_points(gpu, n, xyz, ids) bind(c, name="c_bpack_hodlr_gpu_set_points")
         import :: c_ptr, c_int64_t, c_double
         type(c_ptr), value :: gpu
         integer(c_int64_t) :: n
         real(c_double) :: xyz(*)
         integer(c_int64_t) :: ids(*)
      end subroutine c_bpack_hodlr_gpu_set_points

      subroutine c_bpack_hodlr_gpu_eval_dense(gpu, scale, count, r0, m, c0, n, out) &
         bind(c, name="c_bpack_hodlr_gpu_eval_dense")
         import :: c_ptr, c_int, c_int64_t, c_double
         type(c_ptr), value :: gpu
         real(c_double) :: scale
         integer(c_int) :: count, m(*), n(*)
         integer(c_int64_t) :: r0(*), c0(*)
         type(c_ptr) :: out(*)
      end subroutine c_bpack_hodlr_gpu_eval_dense

      subroutine c_bpack_hodlr_gpu_baca_begin(gpu, scale, mode, variant, nb, r0, m, c0, n, r_est) &
         bind(c, name="c_bpack_hodlr_gpu_baca_begin")
         import :: c_ptr, c_int, c_int64_t, c_double
         type(c_ptr), value :: gpu
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
   end interface

   !> the device indices of the forward blocks of a rank (HODLR_gpu_block_table)
   type hodlr_gpu_itab
      integer, allocatable :: i(:)
   end type hodlr_gpu_itab

contains

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

      call assert(option%HODLR_use_gpu >= 0 .and. option%HODLR_use_gpu <= 2, 'HODLR_use_gpu must be 0, 1 or 2')
      if (option%HODLR_use_gpu == 0) return
      call assert(BPACK_gpu_available(), &
         'HODLR_use_gpu>0 requires the double or double complex library of a build with enable_h2_gpu=ON')
      call assert(option%format == HODLR, 'HODLR_use_gpu>0 requires format=1 (HODLR)')
      call assert(option%LRlevel == 0, 'HODLR_use_gpu>0 requires LRlevel=0 (HODBF has no GPU path)')
      call assert(option%use_zfp /= 1, 'HODLR_use_gpu>0 does not support ZFP-compressed dense blocks')
      call assert(option%HODLR_gpu_pieces >= 1, 'HODLR_gpu_pieces must be at least 1')
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
      character(len=16) :: value
      integer :: status
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
      value = ''
      call get_environment_variable('HODLR_GPU_CHECK', value, status=status)
      check = status == 0 .and. len_trim(value) > 0 .and. trim(value) /= '0'
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
               ! (the construction's, HODLR_GPU_KEEP_FACTORS) and the block's columns are the column child's rows
               devv = .false.
               if (shared .and. k > 0) then
                  devv = HODLR_gpu_keep_factors() .and. blk%headn == colchild%headm
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
               elseif (k > 0) then  ! (on the host; with HODLR_GPU_CHECK also next to the device's, compared at the upload)
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
      ! (HODLR_GPU_CHECK: the way back to the host, compared with the host copies; then taken, the host
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
   !> host copies (which must be filled).  For HODLR_GPU_CHECK
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
   !> ho_bf1%levels(:)%SymFactor (HODLR_GPU_CHECK): the largest relative
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
   !> ho_bf1 (HODLR_GPU_CHECK): per level, the largest relative differences
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

   !> Print the relative difference of a GPU result to the CPU one (HODLR_GPU_CHECK)
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
   !> msh%xyz is freed), for the HODLR construction with a device kernel;
   !> nothing without 3D coordinates in msh%xyz.
   subroutine HODLR_gpu_set_points(bmat, option, msh)
      type(Bmatrix)::bmat
      type(Hoption)::option
      type(mesh)::msh
      real(c_double), allocatable :: xyz(:, :)
      integer(c_int64_t), allocatable :: ids(:)
      integer(c_int64_t) :: n
      integer i

      if (option%HODLR_use_gpu <= 0 .or. option%format /= HODLR .or. .not. c_associated(bmat%gpu)) return
      if (.not. allocated(msh%xyz)) return
      if (size(msh%xyz, 1) /= 3 .or. lbound(msh%xyz, 2) > 1 .or. ubound(msh%xyz, 2) < msh%Nunk) return
      if (.not. allocated(msh%new2old)) return
      n = msh%Nunk
      allocate (xyz(3, max(1, msh%Nunk)), ids(max(1, msh%Nunk)))
      do i = 1, msh%Nunk
         xyz(:, i) = msh%xyz(:, msh%new2old(i))
         ids(i) = msh%new2old(i) - 1
      enddo
      call c_bpack_hodlr_gpu_set_points(bmat%gpu, n, xyz, ids)
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
   !> built on the GPU whose host copies are unfilled (HODLR_GPU_DEFER_HOST):
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
      type(matrixblock), pointer::blocks, root
      type(nod), pointer::cur, curr, curc
      class(*), pointer::ptrr, ptrc
      type(block_ptr), allocatable :: blist(:)
      type(hodlr_gpu_itab), allocatable :: tab(:)
      integer, allocatable :: q_owner(:), q_int(:, :), q_blk(:), q_hrole(:), q_row(:), order(:)
      integer, allocatable :: scount(:), rcount(:), sdisp(:), rdisp(:), svcount(:), rvcount(:), svdisp(:), rvdisp(:)
      integer, allocatable :: r_int(:, :), sbuf(:), rbuf(:)
      integer(c_int), allocatable :: grows(:)
      integer, allocatable :: rmap(:)
      DT, allocatable, target :: vsend(:), vrecv(:), gout(:)
      integer nblk, nq, nqmax, b, nn, t, g, r, lc, ii, iis, role, hrole, pp, lo, hi, mid, np, nproc, ierr, q, k, nrecv
      integer ncache
      integer nloc, idx_start_glo, run0, run1, pos, c, m
      integer(c_int) :: idx_c, role_c, n_c, k_c
      integer(c_int64_t) :: n_loc
      logical view

      nproc = ptree%nproc
      root => ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)
      np = size(root%N_p, 1)
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
                     ! the rank owning row g of the matrix
                     lo = 1
                     hi = np
                     do while (lo < hi)
                        mid = (lo + hi)/2
                        if (root%N_p(mid, 2) < g - root%headn + 1) then
                           lo = mid + 1
                        else
                           hi = mid
                        endif
                     enddo
                     nq = nq + 1
                     q_owner(nq) = ptree%pgrp(root%pgno)%head + lo - 1
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
         call assert(lc >= 1 .and. lc <= ho_bf1%Maxlevel, 'HODLR_gpu_fill_rows: a request outside the low-rank levels')
         call assert(iis >= lbound(tab(lc)%i, 1) .and. iis <= ubound(tab(lc)%i, 1), &
                     'HODLR_gpu_fill_rows: a request for a block of another rank')
         idx_c = tab(lc)%i(iis)
         call assert(idx_c >= 0, 'HODLR_gpu_fill_rows: a request for a block not on this device')
         role_c = r_int(3, run0)
         n_c = run1 - run0 + 1
         allocate (grows(n_c))
         do t = run0, run1
            grows(t - run0 + 1) = r_int(4, t) - idx_start_glo
         enddo
         k = r_int(5, run0)
         allocate (gout(max(1, n_c*k)))
         call c_bpack_hodlr_gpu_gather_rows(ho_bf1%gpu, idx_c, role_c, n_c, grows, c_loc(gout(1)), k_c)
         call assert(k_c == k, 'HODLR_gpu_fill_rows: a block of another rank on the host than on the device')
         do t = 1, n_c  ! (row t of gout, n_c x k: the k values of request run0 + t - 1)
            do c = 1, k
               vsend(pos + c) = gout(t + (c - 1)*n_c)
            enddo
            pos = pos + k
         enddo
         deallocate (grows, gout)
         run0 = run1 + 1
      enddo
      ! the rows back to the requesters, into the host arrays
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
      pos = 0
      do t = 1, nq
         m = order(t)
         blocks => blist(q_blk(m))%ptr
         k = q_int(5, m)
         if (q_hrole(m) == 0) then
            blocks%gpu_u(q_row(m), 1:k) = vrecv(pos + 1:pos + k)
         else
            blocks%gpu_v(q_row(m), 1:k) = vrecv(pos + 1:pos + k)
         endif
         pos = pos + k
      enddo
      deallocate (blist, q_owner, q_int, q_blk, q_hrole, q_row, order, scount, rcount, sdisp, rdisp, sbuf, rbuf, r_int)
      deallocate (svcount, rvcount, svdisp, rvdisp, vsend, vrecv, tab)
   end subroutine HODLR_gpu_fill_rows

   !> HODLR_GPU_CHECK: the host factors of the low-rank blocks of this rank
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
   !> built on the GPU, which stay on the GPU (HODLR_GPU_DEFER_HOST)
   subroutine HODLR_gpu_require_host(ho_bf1, what)
      type(hobf)::ho_bf1
      character(len=*) :: what

      if (.not. ho_bf1%gpu_host_stale) return
      call assert(.false., 'HODLR GPU: '//what//' reads the host copies of the factors, which stay on the GPU '// &
                  '(HODLR_GPU_DEFER_HOST=1). TODO: not supported yet; run with HODLR_GPU_DEFER_HOST=0')
   end subroutine HODLR_gpu_require_host

   !> HODLR_GPU_KEEP_FACTORS (environment): 1 (default) keeps the device
   !> copies of the low-rank factors built on the GPU until the forward
   !> blocks go to the device (HODLR_gpu_upload_forward), which then copies
   !> them on the device instead of uploading the host arrays; 0 frees them
   !> (less device memory during the construction)
   logical function HODLR_gpu_keep_factors()
      integer, save :: mode = -1
      character(len=16) :: value
      integer :: status

      if (mode < 0) then
         value = ''
         call get_environment_variable('HODLR_GPU_KEEP_FACTORS', value, status=status)
         mode = 1
         if (status == 0 .and. len_trim(value) > 0) read (value, *, iostat=status) mode
         if (status /= 0 .or. mode < 0) mode = 1
      endif
      HODLR_gpu_keep_factors = mode > 0
   end function HODLR_gpu_keep_factors

   !> HODLR_GPU_DEFER_HOST (environment): 1 (default) leaves the host copies
   !> of the low-rank factors built on the GPU unfilled: the device keeps
   !> the only copy (HODLR_gpu_upload_forward moves it into the forward
   !> blocks), BPACK_ExtractElement takes the rows it reads from the GPUs
   !> (HODLR_gpu_fill_rows), and the CPU multiply and factorization stop
   !> (HODLR_gpu_require_host); 0 downloads them during the construction, as
   !> always with HODLR_GPU_CHECK (the CPU comparisons read them),
   !> ErrFillFull=1 or HODLR_GPU_KEEP_FACTORS=0
   logical function HODLR_gpu_defer_host(option)
      type(Hoption)::option
      integer, save :: mode = -1
      character(len=16) :: value
      integer :: status

      if (mode < 0) then
         value = ''
         call get_environment_variable('HODLR_GPU_DEFER_HOST', value, status=status)
         mode = 1
         if (status == 0 .and. len_trim(value) > 0) read (value, *, iostat=status) mode
         if (status /= 0 .or. mode < 0) mode = 1
      endif
      HODLR_gpu_defer_host = mode > 0 .and. HODLR_gpu_keep_factors() .and. HODLR_gpu_check_level() <= 0 &
                             .and. option%ErrFillFull /= 1
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

   !> HODLR_GPU_SPLIT (environment): how the GPU construction splits a block
   !> compressed by several ranks into pieces.  0: as the CPU
   !> (LR_HBACA_Leaflevel over the block's group, with the same random
   !> numbers, for comparisons).  1 (default): a block of the symmetric HODLR
   !> over its node's group (both children's ranks, the other half idle
   !> otherwise); with the option HODLR_gpu_pieces > 1 each rank's part split
   !> further and the pieces spread over the group by their estimated cost
   integer function HODLR_gpu_split_mode()
      integer, save :: mode = -1
      character(len=16) :: value
      integer :: status

      if (mode < 0) then
         value = ''
         call get_environment_variable('HODLR_GPU_SPLIT', value, status=status)
         mode = 1
         if (status == 0 .and. len_trim(value) > 0) read (value, *, iostat=status) mode
         if (status /= 0 .or. mode < 0) mode = 1
      endif
      HODLR_gpu_split_mode = mode
   end function HODLR_gpu_split_mode

   !> HODLR_GPU_CHECK (environment): 1 compares each GPU result with the CPU
   !> routine it replaces and prints the relative difference; 2 also checks
   !> the transposed products of each multiply
   integer function HODLR_gpu_check_level()
      integer, save :: mode = -1
      character(len=16) :: value
      integer :: status

      if (mode < 0) then
         value = ''
         call get_environment_variable('HODLR_GPU_CHECK', value, status=status)
         mode = 0
         if (status == 0 .and. len_trim(value) > 0) read (value, *, iostat=status) mode
         if (status /= 0) mode = 0
      endif
      HODLR_gpu_check_level = mode
   end function HODLR_gpu_check_level

end module BPACK_GPU
