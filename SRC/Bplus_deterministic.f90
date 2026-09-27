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

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).
!> @file Bplus_deterministic.f90
!> @brief Deterministic butterfly algebra for the H-BF LU (option%bf_algebra=1), following A. Heldring, E. Ubeda and
!> J. M. Rius, "Fast direct high-frequency MoM solver using butterfly algebra", IEEE TAP 74(1), 2026: butterfly
!> recompression, sum, product (BRowSwap, with a truncation after every swap in a canonical form), partitioned block
!> times butterfly, butterfly back-substitution, and split/merge between a butterfly and its 2x2 children at their native
!> levels. Only butterflies owned by a single process are supported (the H-matrix blocks below the distribution level).
!>
!> Storage conventions (level_butterfly = K, all indices global and 1-based):
!> ButterflyU%blocks(t): m_t x r for row leaf t (stage-K node (t,1))
!> ButterflyV%blocks(j): n_j x r for column leaf j (stage-0 node (1,j)), the transpose of the row basis
!> ButterflyKerl(l)%blocks(i,j), 1<=i<=2^l, 1<=j<=2^(K-l+1): maps the stage-(l-1) node ((i+1)/2, j) to the stage-l node (i, (j+1)/2)

#include "ButterflyPACK_config.fi"
module Bplus_deterministic
   use BPACK_DEFS
   use MISC_Utilities
   use BPACK_Utilities
   use Bplus_randomizedop
   use Bplus_factor
   implicit none

   !>**** one interior factor of a butterfly product chain (BMult). Its nodes are triples (t,u,v) of an X-row, middle and
   !> Y-column tree node at levels s(1:3), sum(s)=K. The factor maps input state s to s+e(inc)-e(dec): output node o gets
   !> the two children (in tree dec) of its parent (in tree inc)
   type bfd_fac
      integer :: s(3) = 0 !< input state (levels of the X-row, middle and Y-column trees)
      integer :: inc = 0 !< the tree refined by this factor
      integer :: dec = 0 !< the tree coarsened by this factor
      type(butterflymatrix), allocatable :: blk(:, :) !< blk(d,o): block from input slot d (1:2) to output node o
   end type bfd_fac

   !>**** a product chain U*f(1)*...*f(nf)*V of 3-tree factors (BMult)
   type bfd_chain
      integer :: K = 0 !< number of butterfly levels
      integer :: nf = 0 !< number of interior factors
      real(kind=8) :: scale = 0d0 !< truncation scale of the swaps (BF_TruncRank)
      type(butterflymatrix), allocatable :: U(:) !< U(t): m_t x r, state (K,0,0)
      type(butterflymatrix), allocatable :: V(:) !< V(v): n_v x r (transposed), state (0,0,K)
      type(bfd_fac), allocatable :: f(:) !< interior factors, f(1) next to U
   end type bfd_chain

   !>**** a product alpha*F(1)*F(2)*...*F(n) of butterflies (F(k) is replaced by I+F(k) if iplus(k)), applied to vectors
   !> by BFD_ChainRand_MVP (used by the diagnostics). All factors live on one process group with conforming layouts.
   type bfd_chainop
      integer :: n = 0
      type(block_ptr) :: f(8)
      logical :: iplus(8) = .false.
      DT :: alpha
   end type bfd_chainop

   !>**** distribution of the nodes of a 3-tree state s = (s_a, s_b, s_c) over the 2^p chain ranks of a distributed BMult:
   !> the owner (0-based) of node q is the concatenation of the top nb(k) bits of q(tr(k)), k = 1..nseg (trees 1 = a,
   !> 2 = b, 3 = c; ButterflyPACK's 'C' layout of a butterfly is [(rows),(columns)], its 'R' layout [(columns),(rows)]).
   !> The local index of a node is formed by its remaining bits, tree a first, so that no index tables are needed.
   type bfd_layout
      integer :: nseg = 0
      integer :: tr(3) = 0
      integer :: nb(3) = 0
   end type bfd_layout

   !>**** one block of a distributed-BMult exchange: destination (rank in the process group's communicator), integer keys
   type bfd_item
      integer :: dst = 0
      integer :: hdr(4) = 0
      DT, allocatable :: mat(:, :)
   end type bfd_item

   type bfd_items
      integer :: n = 0
      type(bfd_item), allocatable :: it(:)
   end type bfd_items

   !>**** a chain factor distributed by output node: blk(d, li) for the output node with local index li in layout lay
   type bfd_dfac
      integer :: s(3) = 0
      integer :: inc = 0
      integer :: dec = 0
      type(bfd_layout) :: lay
      type(butterflymatrix), allocatable :: blk(:, :)
   end type bfd_dfac

   !>**** a BMult chain distributed over the 2^p chain ranks of a process group; chain rank v is the group rank v*stride
   !> (myv = -1 on the other ranks). U and V are the local leaves (states (K,0,0) and (0,0,K)) in layouts Ulay, Vlay.
   type bfd_dchain
      integer :: K = 0
      integer :: p = 0
      integer :: nf = 0
      integer :: stride = 1
      integer :: myv = -1
      integer :: comm = 0
      integer :: pgno = 0
      type(bfd_layout) :: Ulay, Vlay
      type(butterflymatrix), allocatable :: U(:), V(:)
      type(bfd_dfac), allocatable :: f(:)
   end type bfd_dchain

   !>**** wall-clock breakdown of the deterministic algebra (reset and printed by HODLR_factorization with verbosity >= 1);
   !> BFD_rankgrow = largest rank of an updated forward block before and after the final recompression of BFD_Sblock
   integer, parameter :: BFD_NT = 15
   real(kind=8) :: BFD_time(BFD_NT) = 0d0
   integer :: BFD_count(BFD_NT) = 0
   integer :: BFD_rankgrow(2) = 0
   !> recursion statistics of BFD_HxBF (1) and BFD_Lsolve (2) per nesting depth d (depth of the combined HxBF/Lsolve
   !> recursion below the calling H-LU operation): calls, stored entries and largest block dimension of the butterfly
   !> argument, and inclusive flops and time (H-BF factorization, printed by BFD_recstats)
   integer, parameter :: BFD_RD = 40
   integer :: BFD_rec_depth = 0
   integer(kind=8) :: BFD_rec_calls(0:BFD_RD - 1, 2) = 0
   integer :: BFD_rec_rank(0:BFD_RD - 1, 2) = 0
   real(kind=8) :: BFD_rec_size(0:BFD_RD - 1, 2) = 0d0, BFD_rec_flop(0:BFD_RD - 1, 2) = 0d0, BFD_rec_time(0:BFD_RD - 1, 2) = 0d0
   !> breakdown of BFD_Multiply (H-BF products) per H-level of the target (BFD_mul_level, -1 outside BFD_Multiply):
   !> calls, flops and time of its parts and of the BMult steps (called directly, or nested in BFD_HxBF), and the
   !> largest block dimension of the BMult chain after the swaps (direct, nested); printed by BFD_mulstats
   integer, parameter :: BFD_MC = 22
   integer :: BFD_mul_level = -1
   integer(kind=8) :: BFD_mul_cnt(0:BFD_RD - 1, BFD_MC) = 0
   real(kind=8) :: BFD_mul_flop(0:BFD_RD - 1, BFD_MC) = 0d0, BFD_mul_time(0:BFD_RD - 1, BFD_MC) = 0d0
   integer :: BFD_mul_rank(0:BFD_RD - 1, 2) = 0
   !> per H-level of the target (verbosity >= 2): ratio ||P||/||C|| of the norms of a product and of its target block (sum
   !> of log10, largest, smallest), the largest and summed rank of the products, and the largest and smallest ratio of the
   !> stored norm of the target (BFD_Hmat_normest, before the LU) to a fresh estimate (BFD_mulstats)
   real(kind=8) :: BFD_mul_nlog(0:BFD_RD - 1) = 0d0, BFD_mul_nmax(0:BFD_RD - 1) = 0d0, BFD_mul_nmin(0:BFD_RD - 1) = 1d300
   integer :: BFD_mul_prank(0:BFD_RD - 1) = 0
   real(kind=8) :: BFD_mul_prsum(0:BFD_RD - 1) = 0d0
   real(kind=8) :: BFD_mul_cmax(0:BFD_RD - 1) = 0d0, BFD_mul_cmin(0:BFD_RD - 1) = 1d300
   !> time of BFD_Hmat_normest
   real(kind=8) :: BFD_norm_time = 0d0

contains

!======================================================================================== dense helpers

   subroutine BFD_matnew(bm, m, n)
      implicit none
      type(butterflymatrix)::bm
      integer m, n
      if (associated(bm%matrix)) deallocate (bm%matrix)
      allocate (bm%matrix(m, n))
      if (m > 0 .and. n > 0) bm%matrix = 0
   end subroutine BFD_matnew

   subroutine BFD_matfree(bm)
      implicit none
      type(butterflymatrix)::bm
      if (associated(bm%matrix)) deallocate (bm%matrix)
      nullify (bm%matrix)
   end subroutine BFD_matfree

   subroutine BFD_matset(bm, A)
      implicit none
      type(butterflymatrix)::bm
      DT::A(:, :)
      call BFD_matnew(bm, size(A, 1), size(A, 2))
      if (size(A, 1) > 0 .and. size(A, 2) > 0) bm%matrix = A
   end subroutine BFD_matset

   subroutine BFD_matmove(src, dst)
      implicit none
      type(butterflymatrix)::src, dst
      call BFD_matfree(dst)
      dst%matrix => src%matrix
      nullify (src%matrix)
   end subroutine BFD_matmove

   subroutine BFD_eye(bm, n)
      implicit none
      type(butterflymatrix)::bm
      integer n, i
      call BFD_matnew(bm, n, n)
      do i = 1, n
         bm%matrix(i, i) = BPACK_cone
      enddo
   end subroutine BFD_eye

   !>**** C = alpha*op(A)*op(B) + beta*C
   subroutine BFD_gemm(transa, transb, A, B, C, alpha, beta, stats)
      implicit none
      character transa, transb
      DT::A(:, :), B(:, :), C(:, :)
      DT alpha, beta
      type(Hstat)::stats
      integer m, n, k
      real(kind=8)::flop
      if (transa == 'N') then
         m = size(A, 1)
         k = size(A, 2)
      else
         m = size(A, 2)
         k = size(A, 1)
      endif
      if (transb == 'N') then
         n = size(B, 2)
      else
         n = size(B, 1)
      endif
      if (m == 0 .or. n == 0) return
      if (k == 0) then
         if (beta == BPACK_czero) then
            C(1:m, 1:n) = 0
         else
            C(1:m, 1:n) = beta*C(1:m, 1:n)
         endif
         return
      endif
      flop = 0
      call gemmf90(A, size(A, 1), B, size(B, 1), C, size(C, 1), transa, transb, m, n, k, alpha, beta, flop=flop)
      !$omp atomic
      stats%Flop_Tmp = stats%Flop_Tmp + flop
   end subroutine BFD_gemm

   !>**** move the flops accumulated in stats%Flop_Tmp into stats%Flop_Factor. Used around LR_Sblock and LR_minusBC,
   !> which reset stats%Flop_Tmp on entry and add their own flops to stats%Flop_Factor.
   subroutine BFD_flops_flush(stats)
      implicit none
      type(Hstat)::stats
      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
      stats%Flop_Tmp = 0
   end subroutine BFD_flops_flush

   !>**** bm = op(A)*op(B)
   subroutine BFD_mulnew(transa, A, transb, B, bm, stats)
      implicit none
      character transa, transb
      DT::A(:, :), B(:, :)
      type(butterflymatrix)::bm
      type(Hstat)::stats
      integer m, n
      if (transa == 'N') then
         m = size(A, 1)
      else
         m = size(A, 2)
      endif
      if (transb == 'N') then
         n = size(B, 2)
      else
         n = size(B, 1)
      endif
      call BFD_matnew(bm, m, n)
      call BFD_gemm(transa, transb, A, B, bm%matrix, BPACK_cone, BPACK_czero, stats)
   end subroutine BFD_mulnew

   !>**** bm = T*bm
   subroutine BFD_lmul(T, bm, stats)
      implicit none
      DT::T(:, :)
      type(butterflymatrix)::bm, tmp
      type(Hstat)::stats
      call BFD_matmove(bm, tmp)
      call BFD_mulnew('N', T, 'N', tmp%matrix, bm, stats)
      call BFD_matfree(tmp)
   end subroutine BFD_lmul

   !>**** bm = bm*T
   subroutine BFD_rmul(bm, T, stats)
      implicit none
      DT::T(:, :)
      type(butterflymatrix)::bm, tmp
      type(Hstat)::stats
      call BFD_matmove(bm, tmp)
      call BFD_mulnew('N', tmp%matrix, 'N', T, bm, stats)
      call BFD_matfree(tmp)
   end subroutine BFD_rmul

   !>**** bm = [A, beta*B]
   subroutine BFD_hcat(A, B, beta, bm)
      implicit none
      DT::A(:, :), B(:, :)
      DT beta
      type(butterflymatrix)::bm
      integer n1
      call assert(size(A, 1) == size(B, 1), 'BFD_hcat: row mismatch')
      n1 = size(A, 2)
      call BFD_matnew(bm, size(A, 1), n1 + size(B, 2))
      if (size(A, 1) > 0) then
         bm%matrix(:, 1:n1) = A
         bm%matrix(:, n1 + 1:) = beta*B
      endif
   end subroutine BFD_hcat

   !>**** bm = [A; B]
   subroutine BFD_vcat(A, B, bm)
      implicit none
      DT::A(:, :), B(:, :)
      type(butterflymatrix)::bm
      integer m1
      call assert(size(A, 2) == size(B, 2), 'BFD_vcat: column mismatch')
      m1 = size(A, 1)
      call BFD_matnew(bm, m1 + size(B, 1), size(A, 2))
      if (size(A, 2) > 0) then
         bm%matrix(1:m1, :) = A
         bm%matrix(m1 + 1:, :) = B
      endif
   end subroutine BFD_vcat

   !>**** bm = blkdiag(A, B)
   subroutine BFD_blkdiag(A, B, bm)
      implicit none
      DT::A(:, :), B(:, :)
      type(butterflymatrix)::bm
      integer m1, n1
      m1 = size(A, 1)
      n1 = size(A, 2)
      call BFD_matnew(bm, m1 + size(B, 1), n1 + size(B, 2))
      bm%matrix(1:m1, 1:n1) = A
      bm%matrix(m1 + 1:, n1 + 1:) = B
   end subroutine BFD_blkdiag

   !>**** SVD A = W*diag(S)*Z (A is not modified); r = min(m,n), or the truncation rank if trunc (BF_TruncRank with the
   !> optional scale)
   subroutine BFD_svd(A, W, S, Z, r, tol, trunc, stats, scale)
      implicit none
      DT::A(:, :)
      DT, allocatable::W(:, :), Z(:, :)
      DTR, allocatable::S(:)
      integer r, m, n, mn
      real(kind=8)::tol, flop
      real(kind=8), optional::scale
      logical trunc
      type(Hstat)::stats
      m = size(A, 1)
      n = size(A, 2)
      mn = min(m, n)
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
      allocate (W(m, max(mn, 1)), Z(max(mn, 1), n), S(max(mn, 1)))
      W = 0
      Z = 0
      S = 0
      r = mn
      if (mn == 0) return
      flop = 0
      call gesvd_robust(A, S, W, Z, m, n, mn, flop=flop)
      !$omp atomic
      stats%Flop_Tmp = stats%Flop_Tmp + flop
      if (trunc) r = BF_TruncRank(S, mn, tol, scale)
   end subroutine BFD_svd

   !>**** thin QR A = Q*R (A is not modified), Q: m x k, R: k x n, k = min(m,n)
   subroutine BFD_qr(A, Q, R, stats)
      implicit none
      DT::A(:, :)
      DT, allocatable::Q(:, :), R(:, :), W(:, :), tau(:)
      integer m, n, k, i
      real(kind=8)::flop, flop1
      type(Hstat)::stats
      m = size(A, 1)
      n = size(A, 2)
      k = min(m, n)
      if (allocated(Q)) deallocate (Q)
      if (allocated(R)) deallocate (R)
      allocate (Q(m, k), R(k, n))
      Q = 0
      R = 0
      if (k == 0) return
      allocate (W(m, n), tau(k))
      W = A
      flop = 0
      call geqrff90(W, tau, flop=flop)
      do i = 1, k
         R(i, i:n) = W(i, i:n)
      enddo
      flop1 = 0
      call un_or_gqrf90(W, tau, m, k, k, flop=flop1)
      !$omp atomic
      stats%Flop_Tmp = stats%Flop_Tmp + flop + flop1
      Q = W(:, 1:k)
      deallocate (W, tau)
   end subroutine BFD_qr

!======================================================================================== tolerances and diagnostics (verbosity >= 3)

   !>**** internal truncation tolerance of a deterministic operation: its intermediate stages are truncated with
   !> tol_rand*tolfac (the result is recompressed once more with tol_rand). The recursions of BFD_HxBF, BFD_Lsolve and
   !> BFD_SubHH pass this tolerance unchanged to their children (a fixed threshold, as in Heldring et al., IEEE TAP 2026):
   !> tightening it per recursion level lets the ranks of the pieces grow with the recursion depth, i.e. with log N.
   real(kind=8) function BFD_tolfac()
      implicit none
      BFD_tolfac = 0.1d0
   end function BFD_tolfac

   !>**** scale of the truncation of a product P added to an H-block C (BFD_Multiply, BFD_SubHH): P is truncated relative to
   !> max(largest singular value of the truncated block, BFD_scalefac()*||C||) (BF_TruncRank), not relative to P itself.
   !> The two-hop products of the Schur updates have butterfly ranks that grow with N at a truncation relative to P, but
   !> the growing part lies far below ||C||. The stored blocks are truncated relative to the largest singular value of
   !> every butterfly node, which for the far interactions within C is much smaller than ||C||: with the internal
   !> tolerance tol_rand*BFD_tolfac(), 0.01 keeps the parts of P that matter for the nodes of C down to 1d-3*||C||
   !> (semi-circle H-BF, tol_rand=1d-2: the accuracy of the truncation relative to P, product ranks independent of N)
   real(kind=8) function BFD_scalefac()
      implicit none
      BFD_scalefac = 0.01d0
   end function BFD_scalefac

   !>**** expensive diagnostics of the deterministic algebra (dense or random-vector checks of every operation)
   logical function BFD_debug(option)
      implicit none
      type(Hoption)::option
      BFD_debug = option%verbosity >= 3
   end function BFD_debug

   subroutine BFD_tadd(i, t0)
      implicit none
      integer i
      real(kind=8)::t0, dt
      dt = MPI_Wtime() - t0
      !$omp atomic
      BFD_time(i) = BFD_time(i) + dt
      !$omp atomic
      BFD_count(i) = BFD_count(i) + 1
   end subroutine BFD_tadd

   !>**** mode 0: reset the breakdown timers; mode 1: print them (max over processes, verbosity >= 1)
   subroutine BFD_timers(ptree, option, mode)
      implicit none
      type(proctree)::ptree
      type(Hoption)::option
      integer mode, i, ierr
      character(len=24), parameter :: names(BFD_NT) = [character(len=24) :: 'Sblock total', 'Sblock leaf level', &
         'Sblock 1-level nodes', 'Sblock split/merge', '(I+C)^-1 total', 'BMult build', 'BMult sweeps', 'BMult swaps', &
         'BMult merges', 'BMult extract', 'recompress', 'sum', ' recompress: patterns', ' recompress: sweep 1', ' recompress: sweep 2']
      if (mode == 0) then
         BFD_time = 0
         BFD_count = 0
         BFD_rankgrow = 0
         return
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_time, BFD_NT, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_count, BFD_NT, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rankgrow, 2, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) then
         write (*, *) 'deterministic butterfly algebra, time breakdown (max over processes):'
         do i = 1, BFD_NT
            write (*, '(A26,F10.3,A,I8)') names(i), BFD_time(i), ' s   calls ', BFD_count(i)
         enddo
         write (*, '(A,I6,A,I6)') '   Sblock: largest rank before / after the final recompression ', BFD_rankgrow(1), ' / ', BFD_rankgrow(2)
      endif
   end subroutine BFD_timers

   !>**** stored entries sz and largest block dimension rk of the U, V and kernel blocks of a butterfly
   subroutine BFD_bfstat(X, sz, rk)
      implicit none
      type(matrixblock)::X
      real(kind=8)::sz
      integer rk, i, j, l
      sz = 0
      rk = 0
      if (allocated(X%ButterflyU%blocks)) then
         do i = 1, size(X%ButterflyU%blocks)
            if (.not. associated(X%ButterflyU%blocks(i)%matrix)) cycle
            sz = sz + size(X%ButterflyU%blocks(i)%matrix)
            rk = max(rk, size(X%ButterflyU%blocks(i)%matrix, 2))
         enddo
      endif
      if (allocated(X%ButterflyV%blocks)) then
         do i = 1, size(X%ButterflyV%blocks)
            if (.not. associated(X%ButterflyV%blocks(i)%matrix)) cycle
            sz = sz + size(X%ButterflyV%blocks(i)%matrix)
            rk = max(rk, size(X%ButterflyV%blocks(i)%matrix, 2))
         enddo
      endif
      if (allocated(X%ButterflyKerl)) then
         do l = 1, size(X%ButterflyKerl)
            if (.not. allocated(X%ButterflyKerl(l)%blocks)) cycle
            do j = 1, size(X%ButterflyKerl(l)%blocks, 2)
               do i = 1, size(X%ButterflyKerl(l)%blocks, 1)
                  if (.not. associated(X%ButterflyKerl(l)%blocks(i, j)%matrix)) cycle
                  sz = sz + size(X%ButterflyKerl(l)%blocks(i, j)%matrix)
                  rk = max(rk, size(X%ButterflyKerl(l)%blocks(i, j)%matrix, 1), size(X%ButterflyKerl(l)%blocks(i, j)%matrix, 2))
               enddo
            enddo
         enddo
      endif
   end subroutine BFD_bfstat

   !>**** record a call of BFD_HxBF (kind 1) or BFD_Lsolve (kind 2) with butterfly argument X at the current depth d;
   !> f0 and t0 are the flop counter and time at entry, passed back to BFD_rec_exit
   subroutine BFD_rec_enter(kind, X, stats, d, f0, t0)
      implicit none
      integer kind, d, rk
      type(matrixblock)::X
      type(Hstat)::stats
      real(kind=8)::f0, t0, sz
      d = min(BFD_rec_depth, BFD_RD - 1)
      call BFD_bfstat(X, sz, rk)
      BFD_rec_calls(d, kind) = BFD_rec_calls(d, kind) + 1
      BFD_rec_size(d, kind) = BFD_rec_size(d, kind) + sz
      BFD_rec_rank(d, kind) = max(BFD_rec_rank(d, kind), rk)
      BFD_rec_depth = BFD_rec_depth + 1
      f0 = stats%Flop_Tmp
      t0 = MPI_Wtime()
   end subroutine BFD_rec_enter

   subroutine BFD_rec_exit(kind, stats, d, f0, t0)
      implicit none
      integer kind, d
      type(Hstat)::stats
      real(kind=8)::f0, t0
      BFD_rec_depth = BFD_rec_depth - 1
      BFD_rec_flop(d, kind) = BFD_rec_flop(d, kind) + stats%Flop_Tmp - f0
      BFD_rec_time(d, kind) = BFD_rec_time(d, kind) + MPI_Wtime() - t0
   end subroutine BFD_rec_exit

   !>**** mode 0: reset the HxBF/Lsolve recursion statistics; mode 1: print them (calls, flops summed over processes,
   !> largest block dimension and time max over processes). Flops and time are inclusive of deeper calls; "self" is
   !> the depth-d total minus the depth-(d+1) total of both routines.
   subroutine BFD_recstats(ptree, option, mode)
      implicit none
      type(proctree)::ptree
      type(Hoption)::option
      integer mode, d, ierr, dmax
      real(kind=8)::fs(0:BFD_RD), ts(0:BFD_RD)
      if (mode == 0) then
         BFD_rec_depth = 0
         BFD_rec_calls = 0
         BFD_rec_rank = 0
         BFD_rec_size = 0
         BFD_rec_flop = 0
         BFD_rec_time = 0
         return
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rec_calls, 2*BFD_RD, MPI_INTEGER8, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rec_size, 2*BFD_RD, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rec_flop, 2*BFD_RD, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rec_rank, 2*BFD_RD, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_rec_time, 2*BFD_RD, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID /= Main_ID .or. sum(BFD_rec_calls) == 0) return
      dmax = 0
      do d = 0, BFD_RD - 1
         if (sum(BFD_rec_calls(d, :)) > 0) dmax = d
      enddo
      fs = 0
      ts = 0
      do d = 0, dmax
         fs(d) = sum(BFD_rec_flop(d, :))
         ts(d) = sum(BFD_rec_time(d, :))
      enddo
      write (*, *) 'HxBF / Lsolve recursion per depth: calls, mean stored entries of the butterfly argument, largest block'
      write (*, *) 'dimension, inclusive flops (sum over processes); self flops and time of both routines at this depth'
      write (*, '(A6,2(A12,A12,A6,A11),A11,A10)') 'depth', 'HxBF calls', 'mean size', 'rank', 'incl flops', &
         'Lsol calls', 'mean size', 'rank', 'incl flops', 'self flops', 'self time'
      do d = 0, dmax
         write (*, '(I6,2(I12,Es12.3,I6,Es11.2),Es11.2,F10.1)') d, &
            BFD_rec_calls(d, 1), BFD_rec_size(d, 1)/max(1_8, BFD_rec_calls(d, 1)), BFD_rec_rank(d, 1), BFD_rec_flop(d, 1), &
            BFD_rec_calls(d, 2), BFD_rec_size(d, 2)/max(1_8, BFD_rec_calls(d, 2)), BFD_rec_rank(d, 2), BFD_rec_flop(d, 2), &
            fs(d) - fs(d + 1), ts(d) - ts(d + 1)
      enddo
   end subroutine BFD_recstats

   subroutine BFD_mstart(stats, f0, t0)
      implicit none
      type(Hstat)::stats
      real(kind=8)::f0, t0
      f0 = stats%Flop_Tmp
      t0 = MPI_Wtime()
   end subroutine BFD_mstart

   !>**** add one call of category cat of the BFD_Multiply breakdown, started at (f0, t0), at the current H-level
   subroutine BFD_mstat(cat, stats, f0, t0)
      implicit none
      integer cat, l
      type(Hstat)::stats
      real(kind=8)::f0, t0
      if (BFD_mul_level < 0) return
      l = min(BFD_mul_level, BFD_RD - 1)
      BFD_mul_cnt(l, cat) = BFD_mul_cnt(l, cat) + 1
      BFD_mul_flop(l, cat) = BFD_mul_flop(l, cat) + stats%Flop_Tmp - f0
      BFD_mul_time(l, cat) = BFD_mul_time(l, cat) + MPI_Wtime() - t0
   end subroutine BFD_mstat

   !>**** largest block dimension of the interior factors of a BMult chain
   integer function BFD_ch_maxdim(ch)
      implicit none
      type(bfd_chain)::ch
      integer q, i, j
      BFD_ch_maxdim = 0
      if (.not. allocated(ch%f)) return
      do q = 1, size(ch%f)
         if (.not. allocated(ch%f(q)%blk)) cycle
         do j = 1, size(ch%f(q)%blk, 2)
            do i = 1, size(ch%f(q)%blk, 1)
               if (.not. associated(ch%f(q)%blk(i, j)%matrix)) cycle
               BFD_ch_maxdim = max(BFD_ch_maxdim, size(ch%f(q)%blk(i, j)%matrix, 1), size(ch%f(q)%blk(i, j)%matrix, 2))
            enddo
         enddo
      enddo
   end function BFD_ch_maxdim

   !>**** mode 0: reset the BFD_Multiply breakdown; mode 1: print it per H-level of the target (calls and flops summed
   !> over processes, time max over processes)
   subroutine BFD_mulstats(ptree, option, mode)
      implicit none
      type(proctree)::ptree
      type(Hoption)::option
      integer mode, l, c, ierr
      character(len=30), parameter :: names(BFD_MC) = [character(len=30) :: 'Multiply total', ' Product: BMult', &
         ' Product: H x BF', ' Product: BF x H', ' Product: dense x dense', ' adjust level', ' recompress', &
         '  BMult: build', '  BMult: left sweep', '  BMult: right sweeps', '  BMult: swaps', '  BMult: merges', &
         '  BMult: extract', '  BMult: recompress', '  BMult in HxBF: build', '  BMult in HxBF: left sweep', &
         '  BMult in HxBF: right sweeps', '  BMult in HxBF: swaps', '  BMult in HxBF: merges', '  BMult in HxBF: extract', &
         '  BMult in HxBF: recompress', ' norm of the target (missing)']
      if (mode == 0) then
         BFD_mul_level = -1
         BFD_mul_cnt = 0
         BFD_mul_flop = 0
         BFD_mul_time = 0
         BFD_mul_rank = 0
         BFD_mul_nlog = 0
         BFD_mul_nmax = 0
         BFD_mul_nmin = 1d300
         BFD_mul_prank = 0
         BFD_mul_prsum = 0
         BFD_mul_cmax = 0
         BFD_mul_cmin = 1d300
         return
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_cnt, BFD_RD*BFD_MC, MPI_INTEGER8, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_flop, BFD_RD*BFD_MC, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_time, BFD_RD*BFD_MC, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_rank, BFD_RD*2, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_nlog, BFD_RD, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_nmax, BFD_RD, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_nmin, BFD_RD, MPI_DOUBLE_PRECISION, MPI_MIN, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_prank, BFD_RD, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_prsum, BFD_RD, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_cmax, BFD_RD, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_mul_cmin, BFD_RD, MPI_DOUBLE_PRECISION, MPI_MIN, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, BFD_norm_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID /= Main_ID .or. sum(BFD_mul_cnt(:, 1)) == 0) return
      write (*, *) 'Multiply (BFD_Multiply) breakdown per H-level of the target: calls, flops (sum over processes), time (max)'
      write (*, '(A,F11.2)') ' norm estimates of the H-blocks before the LU (BFD_Hmat_normest), time:', BFD_norm_time
      do l = 0, BFD_RD - 1
         if (BFD_mul_cnt(l, 1) == 0) cycle
         write (*, '(A,I3,A,2I6)') ' level', l, '   largest BMult chain block dimension (direct, in HxBF):', BFD_mul_rank(l, 1), BFD_mul_rank(l, 2)
         if (option%verbosity >= 2) then
            write (*, '(A,3Es10.2,A,I6,F9.1)') '           ||P||/||C|| (geometric mean, max, min):', 10d0**(BFD_mul_nlog(l)/BFD_mul_cnt(l, 1)), &
               BFD_mul_nmax(l), BFD_mul_nmin(l), '   product rank (max, mean):', BFD_mul_prank(l), BFD_mul_prsum(l)/BFD_mul_cnt(l, 1)
            write (*, '(A,2Es10.2)') '           stored ||C|| / fresh estimate (max, min):', BFD_mul_cmax(l), BFD_mul_cmin(l)
         endif
         do c = 1, BFD_MC
            if (BFD_mul_cnt(l, c) == 0) cycle
            write (*, '(A32,I10,Es11.2,F11.2)') names(c), BFD_mul_cnt(l, c), BFD_mul_flop(l, c), BFD_mul_time(l, c)
         enddo
      enddo
   end subroutine BFD_mulstats

   !>**** dense reconstruction of a single-process butterfly (debugging only)
   subroutine BFD_todense(blk, D, stats)
      implicit none
      type(matrixblock)::blk
      DT, allocatable::D(:, :)
      type(Hstat)::stats
      type(butterflymatrix), allocatable::cur(:), nxt(:)
      integer K, l, i, j, jj, nb, off, M, N, n0
      K = blk%level_butterfly
      nb = 2**K
      M = 0
      N = 0
      do i = 1, nb
         M = M + size(blk%ButterflyU%blocks(i)%matrix, 1)
         N = N + size(blk%ButterflyV%blocks(i)%matrix, 1)
      enddo
      if (allocated(D)) deallocate (D)
      allocate (D(M, N))
      D = 0
      allocate (cur(nb))
      off = 0
      do j = 1, nb
         n0 = size(blk%ButterflyV%blocks(j)%matrix, 1)
         call BFD_matnew(cur(j), size(blk%ButterflyV%blocks(j)%matrix, 2), N)
         cur(j)%matrix(:, off + 1:off + n0) = transpose(blk%ButterflyV%blocks(j)%matrix)
         off = off + n0
      enddo
      do l = 1, K
         allocate (nxt(nb))
         do i = 1, 2**l
            do jj = 1, 2**(K - l)
               call BFD_mulnew('N', blk%ButterflyKerl(l)%blocks(i, 2*jj - 1)%matrix, 'N', cur(((i + 1)/2 - 1)*2**(K - l + 1) + 2*jj - 1)%matrix, nxt((i - 1)*2**(K - l) + jj), stats)
               call BFD_gemm('N', 'N', blk%ButterflyKerl(l)%blocks(i, 2*jj)%matrix, cur(((i + 1)/2 - 1)*2**(K - l + 1) + 2*jj)%matrix, nxt((i - 1)*2**(K - l) + jj)%matrix, BPACK_cone, BPACK_cone, stats)
            enddo
         enddo
         call BFD_matarray_free(cur)
         call move_alloc(nxt, cur)
      enddo
      off = 0
      do i = 1, nb
         n0 = size(blk%ButterflyU%blocks(i)%matrix, 1)
         call BFD_gemm('N', 'N', blk%ButterflyU%blocks(i)%matrix, cur(i)%matrix, D(off + 1:off + n0, :), BPACK_cone, BPACK_czero, stats)
         off = off + n0
      enddo
      call BFD_matarray_free(cur)
   end subroutine BFD_todense

   !>**** dense op(A) of an H-block (debugging only)
   subroutine BFD_hbdense(A, trans, D, ptree, stats)
      implicit none
      type(matrixblock)::A
      character trans
      DT, allocatable::D(:, :), Id(:, :)
      type(proctree)::ptree
      type(Hstat)::stats
      integer i, m, n
      m = A%M
      n = A%N
      if (trans == 'N') then
         allocate (Id(n, n), D(m, n))
         Id = 0
         do i = 1, n
            Id(i, i) = BPACK_cone
         enddo
         D = 0
         call Hmat_block_MVP_dat(A, 'N', A%headm, A%headn, n, Id, n, D, m, BPACK_cone, ptree, stats)
      else
         allocate (Id(m, m), D(n, m))
         Id = 0
         do i = 1, m
            Id(i, i) = BPACK_cone
         enddo
         D = 0
         call Hmat_block_MVP_dat(A, 'T', A%headm, A%headn, m, Id, m, D, n, BPACK_cone, ptree, stats)
      endif
      deallocate (Id)
   end subroutine BFD_hbdense

   subroutine BFD_report(name, D, Dref, K)
      implicit none
      character(*) name
      DT::D(:, :), Dref(:, :)
      integer K
      real(kind=8)::e, nr
      nr = fnorm(Dref, size(Dref, 1), size(Dref, 2))
      e = fnorm(D - Dref, size(Dref, 1), size(Dref, 2))
      if (nr > 0) e = e/nr
      write (*, '(A,A,A,I3,A,Es12.4,A,I6,I6)') ' BFD_DEBUG ', name, ' K=', K, ' relerr=', e, ' size', size(Dref, 1), size(Dref, 2)
   end subroutine BFD_report

   !>**** check that the leaves of a butterfly with cluster metadata are the clusters K levels below (debugging only)
   subroutine BFD_checkleaves(blk, msh, name)
      implicit none
      type(matrixblock)::blk
      type(mesh)::msh
      character(*) name
      integer K, i, g, nbad
      K = blk%level_butterfly
      nbad = 0
      do i = 1, 2**K
         g = blk%row_group*2**K + i - 1
         if (size(blk%ButterflyU%blocks(i)%matrix, 1) /= msh%basis_group(g)%tail - msh%basis_group(g)%head + 1) nbad = nbad + 1
         g = blk%col_group*2**K + i - 1
         if (size(blk%ButterflyV%blocks(i)%matrix, 1) /= msh%basis_group(g)%tail - msh%basis_group(g)%head + 1) nbad = nbad + 1
      enddo
      if (nbad > 0) write (*, *) 'BFD_DEBUG leaves mismatch in ', name, ' K', K, ' groups', blk%row_group, blk%col_group, ' level', blk%level, ' nbad', nbad
   end subroutine BFD_checkleaves

!======================================================================================== butterfly bookkeeping

   !>**** check that a butterfly is owned by one process with global (untruncated) local index ranges
   subroutine BFD_check_seq(blk, ptree)
      implicit none
      type(matrixblock)::blk
      type(proctree)::ptree
      integer level, K
      call assert(blk%style == 2, 'BFD: expecting a butterfly/low-rank block')
      call assert(ptree%pgrp(blk%pgno)%nproc == 1, 'BFD: deterministic butterfly algebra requires blocks owned by one process')
      K = blk%level_butterfly
      call assert(blk%ButterflyU%nblk_loc == 2**K .and. blk%ButterflyV%nblk_loc == 2**K, 'BFD: butterfly U/V are not fully local')
      do level = 1, K
         call assert(blk%ButterflyKerl(level)%nr == 2**level .and. blk%ButterflyKerl(level)%nc == 2**(K - level + 1) &
                     .and. blk%ButterflyKerl(level)%idx_r == 1 .and. blk%ButterflyKerl(level)%idx_c == 1 &
                     .and. blk%ButterflyKerl(level)%inc_r == 1 .and. blk%ButterflyKerl(level)%inc_c == 1, &
                     'BFD: butterfly kernel is not fully local')
      enddo
   end subroutine BFD_check_seq

   !>**** set the metadata of a single-process K-level butterfly and allocate its (empty) block arrays
   subroutine BFD_alloc(blk, K, pgno)
      implicit none
      type(matrixblock)::blk
      integer K, pgno, level, nb
      blk%style = 2
      blk%level_butterfly = K
      blk%level_half = floor_safe(dble(K)/2d0)
      blk%pgno = pgno
      blk%pgno_db = pgno
      nb = 2**K
      blk%ButterflyU%num_blk = nb
      blk%ButterflyU%nblk_loc = nb
      blk%ButterflyU%idx = 1
      blk%ButterflyU%inc = 1
      blk%ButterflyV%num_blk = nb
      blk%ButterflyV%nblk_loc = nb
      blk%ButterflyV%idx = 1
      blk%ButterflyV%inc = 1
      allocate (blk%ButterflyU%blocks(nb))
      allocate (blk%ButterflyV%blocks(nb))
      if (K > 0) then
         allocate (blk%ButterflyKerl(K))
         do level = 1, K
            blk%ButterflyKerl(level)%num_row = 2**level
            blk%ButterflyKerl(level)%num_col = 2**(K - level + 1)
            blk%ButterflyKerl(level)%nr = 2**level
            blk%ButterflyKerl(level)%nc = 2**(K - level + 1)
            blk%ButterflyKerl(level)%idx_r = 1
            blk%ButterflyKerl(level)%inc_r = 1
            blk%ButterflyKerl(level)%idx_c = 1
            blk%ButterflyKerl(level)%inc_c = 1
            allocate (blk%ButterflyKerl(level)%blocks(2**level, 2**(K - level + 1)))
         enddo
      endif
   end subroutine BFD_alloc

   !>**** dimensions, local layouts and accumulated leaf sizes from the U/V blocks
   subroutine BFD_setdims(blk)
      implicit none
      type(matrixblock)::blk
      integer i, nb
      nb = 2**blk%level_butterfly
      if (associated(blk%ms)) deallocate (blk%ms)
      if (associated(blk%ns)) deallocate (blk%ns)
      allocate (blk%ms(nb), blk%ns(nb))
      do i = 1, nb
         blk%ms(i) = size(blk%ButterflyU%blocks(i)%matrix, 1)
         blk%ns(i) = size(blk%ButterflyV%blocks(i)%matrix, 1)
         if (i > 1) then
            blk%ms(i) = blk%ms(i) + blk%ms(i - 1)
            blk%ns(i) = blk%ns(i) + blk%ns(i - 1)
         endif
      enddo
      blk%M = blk%ms(nb)
      blk%N = blk%ns(nb)
      blk%M_loc = blk%M
      blk%N_loc = blk%N
      if (associated(blk%M_p)) deallocate (blk%M_p)
      if (associated(blk%N_p)) deallocate (blk%N_p)
      allocate (blk%M_p(1, 2), blk%N_p(1, 2))
      blk%M_p(1, 1) = 1
      blk%M_p(1, 2) = blk%M
      blk%N_p(1, 1) = 1
      blk%N_p(1, 2) = blk%N
      blk%rankmax = 0
      blk%rankmin = 0
      do i = 1, nb
         blk%rankmax = max(blk%rankmax, size(blk%ButterflyU%blocks(i)%matrix, 2), size(blk%ButterflyV%blocks(i)%matrix, 2))
      enddo
   end subroutine BFD_setdims

   !>**** copy the cluster metadata of a target block into a result block
   subroutine BFD_copymeta(res, tmpl)
      implicit none
      type(matrixblock)::res, tmpl
      res%level = tmpl%level
      res%row_group = tmpl%row_group
      res%col_group = tmpl%col_group
      res%headm = tmpl%headm
      res%headn = tmpl%headn
      res%pgno = tmpl%pgno
      res%pgno_db = tmpl%pgno_db
      call assert(res%M == tmpl%M .and. res%N == tmpl%N, 'BFD: result dimensions do not match the target block')
   end subroutine BFD_copymeta

   !>**** replace the content of the butterfly C by res (res is emptied)
   subroutine BFD_Install(C, res, ptree)
      implicit none
      type(matrixblock)::C, res
      type(proctree)::ptree
      call BFD_copymeta(res, C)
      call BF_delete(C, 1)
      call BF_copy_delete(res, C)
      call BF_delete(res, 1)
      call BF_get_rank(C, ptree)
   end subroutine BFD_Install

   !>**** deep copy of a single-process butterfly
   subroutine BFD_copy(A, B)
      implicit none
      type(matrixblock)::A, B
      integer K, l, i, j
      K = A%level_butterfly
      call BFD_alloc(B, K, A%pgno)
      do i = 1, 2**K
         call BFD_matset(B%ButterflyU%blocks(i), A%ButterflyU%blocks(i)%matrix)
         call BFD_matset(B%ButterflyV%blocks(i), A%ButterflyV%blocks(i)%matrix)
      enddo
      do l = 1, K
         do i = 1, 2**l
            do j = 1, 2**(K - l + 1)
               call BFD_matset(B%ButterflyKerl(l)%blocks(i, j), A%ButterflyKerl(l)%blocks(i, j)%matrix)
            enddo
         enddo
      enddo
      call BFD_setdims(B)
      call BFD_copymeta(B, A)
   end subroutine BFD_copy

   integer function BFD_pattern(blk)
      implicit none
      type(matrixblock)::blk
      if (blk%level_half == blk%level_butterfly) then
         BFD_pattern = 1
      elseif (blk%level_half == 0) then
         BFD_pattern = 2
      else
         BFD_pattern = 3
      endif
   end function BFD_pattern

   !>**** leaf sizes of a cluster at K levels below it
   subroutine BFD_leafsizes(group, K, msh, sz)
      implicit none
      integer group, K, i, g
      integer sz(:)
      type(mesh)::msh
      do i = 1, 2**K
         g = group*2**K + i - 1
         sz(i) = msh%basis_group(g)%tail - msh%basis_group(g)%head + 1
      enddo
   end subroutine BFD_leafsizes

!======================================================================================== BRecompress

   !>**** low-rank recompression U*V^T with a truncated SVD
   subroutine BFD_recompress_lr(blk, tol, stats, scale)
      implicit none
      type(matrixblock)::blk
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(Hstat)::stats
      DT, allocatable::W(:, :), Z(:, :), T(:, :), W2(:, :), Z2(:, :)
      DTR, allocatable::S(:), S2(:)
      integer r, r2, i, M, N
      M = size(blk%ButterflyU%blocks(1)%matrix, 1)
      N = size(blk%ButterflyV%blocks(1)%matrix, 1)
      call BFD_svd(blk%ButterflyU%blocks(1)%matrix, W, S, Z, r, tol, .false., stats)
      do i = 1, r
         Z(i, :) = Z(i, :)*S(i)
      enddo
      allocate (T(r, N))
      call BFD_gemm('N', 'T', Z(1:r, :), blk%ButterflyV%blocks(1)%matrix, T, BPACK_cone, BPACK_czero, stats)
      call BFD_svd(T, W2, S2, Z2, r2, tol, .true., stats, scale)
      do i = 1, r2
         W2(:, i) = W2(:, i)*S2(i)
      enddo
      call BFD_matnew(blk%ButterflyU%blocks(1), M, r2)
      call BFD_gemm('N', 'N', W(:, 1:r), W2(1:r, 1:r2), blk%ButterflyU%blocks(1)%matrix, BPACK_cone, BPACK_czero, stats)
      call BFD_matnew(blk%ButterflyV%blocks(1), N, r2)
      blk%ButterflyV%blocks(1)%matrix = transpose(Z2(1:r2, :))
      deallocate (W, Z, S, T, W2, Z2, S2)
   end subroutine BFD_recompress_lr

   !>**** recompression of a low-rank block U*V^T whose U and V are row-distributed over the process group of blk:
   !> U = Qu*Mu^T and V = Qv*Mv^T (BF_LeafOrth), then a truncated SVD of the small core Mu^T*Mv on the head
   subroutine BFD_recompress_lr_dist(blk, tol, stats, ptree, scale)
      implicit none
      type(matrixblock)::blk
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(Hstat)::stats
      type(proctree)::ptree
      DT, allocatable::Mu(:, :), Mv(:, :), Cc(:, :), UU(:, :), VV(:, :), Wb(:, :), Zb(:, :)
      DTR, allocatable::S(:)
      integer comm, nproc, myid, ierr, ru, rv, mn, r, i, dims(1)
      real(kind=8)::flop

      comm = ptree%pgrp(blk%pgno)%Comm
      nproc = ptree%pgrp(blk%pgno)%nproc
      call MPI_Comm_rank(comm, myid, ierr)
      call BF_LeafOrth(blk%ButterflyU%blocks(1), Mu, .false., tol, comm, nproc, stats)
      call BF_LeafOrth(blk%ButterflyV%blocks(1), Mv, .false., tol, comm, nproc, stats)
      ru = size(Mu, 2)
      rv = size(Mv, 2)
      if (myid == 0) then
         allocate (Cc(ru, rv))
         call BFD_gemm('T', 'N', Mu, Mv, Cc, BPACK_cone, BPACK_czero, stats)
         mn = min(ru, rv)
         allocate (UU(ru, mn), VV(mn, rv), S(mn))
         flop = 0
         call gesvd_robust(Cc, S, UU, VV, ru, rv, mn, flop=flop)
         !$omp atomic
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         dims(1) = BF_TruncRank(S, mn, tol, scale)
      endif
      call MPI_Bcast(dims, 1, MPI_INTEGER, 0, comm, ierr)
      r = dims(1)
      allocate (Wb(ru, r), Zb(rv, r))
      if (myid == 0) then
         do i = 1, r
            Wb(:, i) = UU(:, i)*S(i)
            Zb(:, i) = VV(i, :)
         enddo
         deallocate (Cc, UU, VV, S)
      endif
      call MPI_Bcast(Wb, ru*r, MPI_DT, 0, comm, ierr)
      call MPI_Bcast(Zb, rv*r, MPI_DT, 0, comm, ierr)
      flop = 0
      call BF_LeafRmul(blk%ButterflyU%blocks(1), Wb, flop)
      call BF_LeafRmul(blk%ButterflyV%blocks(1), Zb, flop)
      !$omp atomic
      stats%Flop_Tmp = stats%Flop_Tmp + flop
      deallocate (Mu, Mv, Wb, Zb)
   end subroutine BFD_recompress_lr_dist

   !>**** BRecompress: sweep A (V->U) makes V and the kernels row-orthonormal without truncation; sweep B (U->V) truncates
   !> every node with the relative tolerance tol (the part already swept is column-orthonormal and the rest row-orthonormal,
   !> so each local truncation is sharp). Both sweeps are BF_MoveSingular_Ker over all levels 0..K+1 and work for
   !> distributed butterflies.
   subroutine BFD_recompress(blk, tol, option, stats, ptree, scale)
      implicit none
      type(matrixblock)::blk
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      integer K, pat0
      real(kind=8)::t0, t1

      t0 = MPI_Wtime()
      K = blk%level_butterfly
      if (K == 0) then
         if (ptree%pgrp(blk%pgno)%nproc == 1) then
            call BFD_recompress_lr(blk, tol, stats, scale)
            call BFD_setdims(blk)
         else
            call BFD_recompress_lr_dist(blk, tol, stats, ptree, scale)
         endif
         call BF_get_rank(blk, ptree)
         call BFD_tadd(11, t0)
         return
      endif

      pat0 = BFD_pattern(blk)
      t1 = MPI_Wtime()
      call BF_ChangePattern(blk, pat0, 1, stats, ptree)
      call BFD_tadd(13, t1)
      t1 = MPI_Wtime()
      call BF_MoveSingular_Ker(blk, 'N', 0, K + 1, ptree, stats, tol, useqr=.true.)
      call BFD_tadd(14, t1)
      t1 = MPI_Wtime()
      call BF_ChangePattern(blk, 1, 2, stats, ptree)
      call BFD_tadd(13, t1)
      t1 = MPI_Wtime()
      call BF_MoveSingular_Ker(blk, 'T', K + 1, 0, ptree, stats, tol, truncate=.true., scale=scale)
      call BFD_tadd(15, t1)
      t1 = MPI_Wtime()
      call BF_ChangePattern(blk, 2, pat0, stats, ptree)
      call BFD_tadd(13, t1)
      if (ptree%pgrp(blk%pgno)%nproc == 1) call BFD_setdims(blk)
      call BF_get_rank(blk, ptree)
      call BFD_tadd(11, t0)
   end subroutine BFD_recompress

!======================================================================================== exact structural operations

   !>**** S = A + beta*B (same level and leaves): U side by side, kernels block diagonal, V stacked
   subroutine BFD_sum(A, B, beta, S, pgno)
      implicit none
      type(matrixblock)::A, B, S
      DT beta
      integer pgno, K, l, i, j
      K = A%level_butterfly
      call assert(B%level_butterfly == K, 'BFD_sum: level mismatch')
      call BFD_alloc(S, K, pgno)
      do i = 1, 2**K
         call BFD_hcat(A%ButterflyU%blocks(i)%matrix, B%ButterflyU%blocks(i)%matrix, beta, S%ButterflyU%blocks(i))
         call assert(size(A%ButterflyV%blocks(i)%matrix, 1) == size(B%ButterflyV%blocks(i)%matrix, 1), 'BFD_sum: leaf mismatch')
         call BFD_hcat(A%ButterflyV%blocks(i)%matrix, B%ButterflyV%blocks(i)%matrix, BPACK_cone, S%ButterflyV%blocks(i))
      enddo
      do l = 1, K
         do i = 1, 2**l
            do j = 1, 2**(K - l + 1)
               call BFD_blkdiag(A%ButterflyKerl(l)%blocks(i, j)%matrix, B%ButterflyKerl(l)%blocks(i, j)%matrix, S%ButterflyKerl(l)%blocks(i, j))
            enddo
         enddo
      enddo
      call BFD_setdims(S)
   end subroutine BFD_sum

   !>**** S = [A, beta*B] (same rows and level); column leaf j of S is [leaf j of A; leaf j of B]
   subroutine BFD_hconcat(A, B, beta, S, pgno)
      implicit none
      type(matrixblock)::A, B, S
      DT beta
      integer pgno, K, l, i, j
      K = A%level_butterfly
      call assert(B%level_butterfly == K, 'BFD_hconcat: level mismatch')
      call BFD_alloc(S, K, pgno)
      do i = 1, 2**K
         call BFD_hcat(A%ButterflyU%blocks(i)%matrix, B%ButterflyU%blocks(i)%matrix, beta, S%ButterflyU%blocks(i))
         call BFD_blkdiag(A%ButterflyV%blocks(i)%matrix, B%ButterflyV%blocks(i)%matrix, S%ButterflyV%blocks(i))
      enddo
      do l = 1, K
         do i = 1, 2**l
            do j = 1, 2**(K - l + 1)
               call BFD_blkdiag(A%ButterflyKerl(l)%blocks(i, j)%matrix, B%ButterflyKerl(l)%blocks(i, j)%matrix, S%ButterflyKerl(l)%blocks(i, j))
            enddo
         enddo
      enddo
      call BFD_setdims(S)
   end subroutine BFD_hconcat

   !>**** T = A^T (plain transpose)
   subroutine BFD_transpose(A, T, pgno)
      implicit none
      type(matrixblock)::A, T
      integer pgno, K, l, i, j
      K = A%level_butterfly
      call BFD_alloc(T, K, pgno)
      do i = 1, 2**K
         call BFD_matset(T%ButterflyU%blocks(i), A%ButterflyV%blocks(i)%matrix)
         call BFD_matset(T%ButterflyV%blocks(i), A%ButterflyU%blocks(i)%matrix)
      enddo
      do l = 1, K
         do j = 1, 2**l
            do i = 1, 2**(K - l + 1)
               call BFD_matset(T%ButterflyKerl(l)%blocks(j, i), transpose(A%ButterflyKerl(K - l + 1)%blocks(i, j)%matrix))
            enddo
         enddo
      enddo
      call BFD_setdims(T)
   end subroutine BFD_transpose

   !>**** exact K -> K+1 on leaves refined once (fine leaf sizes rsz, csz): V restricted to the fine leaves, a new first kernel
   !> of identities, every old kernel block copied to both row children, U restricted to the fine leaves
   subroutine BFD_levelup(A, B, rsz, csz, pgno)
      implicit none
      type(matrixblock)::A, B
      integer rsz(:), csz(:)
      integer pgno, K, l, i, j, jc, off, r
      K = A%level_butterfly
      call BFD_alloc(B, K + 1, pgno)
      do j = 1, 2**(K + 1)
         jc = (j + 1)/2
         off = 0
         if (mod(j, 2) == 0) off = csz(j - 1)
         if (mod(j, 2) == 1) call assert(csz(j) + csz(j + 1) == size(A%ButterflyV%blocks(jc)%matrix, 1), 'BFD_levelup: column leaf mismatch')
         call BFD_matset(B%ButterflyV%blocks(j), A%ButterflyV%blocks(jc)%matrix(off + 1:off + csz(j), :))
         r = size(A%ButterflyV%blocks(jc)%matrix, 2)
         call BFD_eye(B%ButterflyKerl(1)%blocks(1, j), r)
         call BFD_eye(B%ButterflyKerl(1)%blocks(2, j), r)
      enddo
      do l = 1, K
         do i = 1, 2**(l + 1)
            do j = 1, 2**(K - l + 1)
               call BFD_matset(B%ButterflyKerl(l + 1)%blocks(i, j), A%ButterflyKerl(l)%blocks((i + 1)/2, j)%matrix)
            enddo
         enddo
      enddo
      do i = 1, 2**(K + 1)
         jc = (i + 1)/2
         off = 0
         if (mod(i, 2) == 0) off = rsz(i - 1)
         if (mod(i, 2) == 1) call assert(rsz(i) + rsz(i + 1) == size(A%ButterflyU%blocks(jc)%matrix, 1), 'BFD_levelup: row leaf mismatch')
         call BFD_matset(B%ButterflyU%blocks(i), A%ButterflyU%blocks(jc)%matrix(off + 1:off + rsz(i), :))
      enddo
      call BFD_setdims(B)
   end subroutine BFD_levelup

   !>**** exact K -> K-1 on leaves merged pairwise: every stage-l node (i,j) stacks the sibling nodes (i,2j-1), (i,2j)
   subroutine BFD_leveldown(A, B, pgno, stats)
      implicit none
      type(matrixblock)::A, B
      type(Hstat)::stats
      integer pgno, K, l, i, jj, J1, h0, h1, c1, c2, r0
      type(butterflymatrix)::T1, T2, Ttop, Tbot
      K = A%level_butterfly
      call assert(K >= 1, 'BFD_leveldown: needs at least one level')
      call BFD_alloc(B, K - 1, pgno)
      do jj = 1, 2**(K - 1)
         call BFD_blkdiag(A%ButterflyV%blocks(2*jj - 1)%matrix, A%ButterflyV%blocks(2*jj)%matrix, B%ButterflyV%blocks(jj))
      enddo
      do l = 1, K - 1
         do i = 1, 2**l
            do jj = 1, 2**(K - l)
               J1 = (jj + 1)/2
               h0 = size(A%ButterflyKerl(l)%blocks(i, 4*J1 - 3)%matrix, 1)
               h1 = size(A%ButterflyKerl(l)%blocks(i, 4*J1 - 1)%matrix, 1)
               c1 = size(A%ButterflyKerl(l)%blocks(i, 2*jj - 1)%matrix, 2)
               c2 = size(A%ButterflyKerl(l)%blocks(i, 2*jj)%matrix, 2)
               call BFD_matnew(B%ButterflyKerl(l)%blocks(i, jj), h0 + h1, c1 + c2)
               r0 = 0
               if (mod(jj, 2) == 0) r0 = h0
               B%ButterflyKerl(l)%blocks(i, jj)%matrix(r0 + 1:r0 + size(A%ButterflyKerl(l)%blocks(i, 2*jj - 1)%matrix, 1), 1:c1) = A%ButterflyKerl(l)%blocks(i, 2*jj - 1)%matrix
               B%ButterflyKerl(l)%blocks(i, jj)%matrix(r0 + 1:r0 + size(A%ButterflyKerl(l)%blocks(i, 2*jj)%matrix, 1), c1 + 1:c1 + c2) = A%ButterflyKerl(l)%blocks(i, 2*jj)%matrix
            enddo
         enddo
      enddo
      do i = 1, 2**(K - 1)
         call BFD_mulnew('N', A%ButterflyU%blocks(2*i - 1)%matrix, 'N', A%ButterflyKerl(K)%blocks(2*i - 1, 1)%matrix, T1, stats)
         call BFD_mulnew('N', A%ButterflyU%blocks(2*i - 1)%matrix, 'N', A%ButterflyKerl(K)%blocks(2*i - 1, 2)%matrix, T2, stats)
         call BFD_hcat(T1%matrix, T2%matrix, BPACK_cone, Ttop)
         call BFD_mulnew('N', A%ButterflyU%blocks(2*i)%matrix, 'N', A%ButterflyKerl(K)%blocks(2*i, 1)%matrix, T1, stats)
         call BFD_mulnew('N', A%ButterflyU%blocks(2*i)%matrix, 'N', A%ButterflyKerl(K)%blocks(2*i, 2)%matrix, T2, stats)
         call BFD_hcat(T1%matrix, T2%matrix, BPACK_cone, Tbot)
         call BFD_vcat(Ttop%matrix, Tbot%matrix, B%ButterflyU%blocks(i))
      enddo
      call BFD_matfree(T1)
      call BFD_matfree(T2)
      call BFD_matfree(Ttop)
      call BFD_matfree(Tbot)
      call BFD_setdims(B)
   end subroutine BFD_leveldown

   !>**** bring a butterfly with cluster metadata to level Kt (exact level changes, then recompression)
   subroutine BFD_adjustlevel(P, Kt, tol, option, stats, ptree, msh, scale)
      implicit none
      type(matrixblock)::P
      integer Kt
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(matrixblock)::T
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      integer, allocatable::rsz(:), csz(:)
      if (P%level_butterfly == Kt) return
      do while (P%level_butterfly < Kt)
         allocate (rsz(2**(P%level_butterfly + 1)), csz(2**(P%level_butterfly + 1)))
         call BFD_leafsizes(P%row_group, P%level_butterfly + 1, msh, rsz)
         call BFD_leafsizes(P%col_group, P%level_butterfly + 1, msh, csz)
         call BFD_levelup(P, T, rsz, csz, P%pgno)
         call BFD_Install(P, T, ptree)
         deallocate (rsz, csz)
      enddo
      do while (P%level_butterfly > Kt)
         call BFD_leveldown(P, T, P%pgno, stats)
         call BFD_Install(P, T, ptree)
      enddo
      call BFD_recompress(P, tol, option, stats, ptree, scale)
   end subroutine BFD_adjustlevel

!======================================================================================== row split / merge

   !>**** H = rows of half a (1 or 2) of the K-level butterfly B as a (K-1)-level butterfly whose column leaves are the
   !> stage-1 nodes (a,j) of B (V = identity); B = blkdiag(H_1, H_2)*W with W = (kernel 1 of B)*V
   subroutine BFD_rowhalf(B, a, H, pgno)
      implicit none
      type(matrixblock)::B, H
      integer a, pgno, K, hh, t, j, lp, ip
      K = B%level_butterfly
      hh = 2**(K - 1)
      call BFD_alloc(H, K - 1, pgno)
      do t = 1, hh
         call BFD_matset(H%ButterflyU%blocks(t), B%ButterflyU%blocks((a - 1)*hh + t)%matrix)
         call BFD_eye(H%ButterflyV%blocks(t), size(B%ButterflyKerl(1)%blocks(a, 2*t - 1)%matrix, 1))
      enddo
      do lp = 1, K - 1
         do ip = 1, 2**lp
            do j = 1, 2**(K - lp)
               call BFD_matset(H%ButterflyKerl(lp)%blocks(ip, j), B%ButterflyKerl(lp + 1)%blocks((a - 1)*2**lp + ip, j)%matrix)
            enddo
         enddo
      enddo
      call BFD_setdims(H)
   end subroutine BFD_rowhalf

   !>**** inverse of the row split: rows [P1; P2], kernel 1 (a,jf) = P_a%V(ceil(jf/2))^T * W_a(jf), V copied from B
   subroutine BFD_rowmerge(P1, P2, W1, W2, B, O, pgno, stats)
      implicit none
      type(matrixblock)::P1, P2, B, O
      type(butterflymatrix)::W1(:), W2(:)
      type(Hstat)::stats
      integer pgno, K, h, t, jf, lp, ip, j
      K = P1%level_butterfly + 1
      h = 2**(K - 1)
      call BFD_alloc(O, K, pgno)
      do t = 1, h
         call BFD_matset(O%ButterflyU%blocks(t), P1%ButterflyU%blocks(t)%matrix)
         call BFD_matset(O%ButterflyU%blocks(h + t), P2%ButterflyU%blocks(t)%matrix)
      enddo
      do jf = 1, 2**K
         call BFD_mulnew('T', P1%ButterflyV%blocks((jf + 1)/2)%matrix, 'N', W1(jf)%matrix, O%ButterflyKerl(1)%blocks(1, jf), stats)
         call BFD_mulnew('T', P2%ButterflyV%blocks((jf + 1)/2)%matrix, 'N', W2(jf)%matrix, O%ButterflyKerl(1)%blocks(2, jf), stats)
         call BFD_matset(O%ButterflyV%blocks(jf), B%ButterflyV%blocks(jf)%matrix)
      enddo
      do lp = 1, K - 1
         do ip = 1, 2**lp
            do j = 1, 2**(K - lp)
               call BFD_matset(O%ButterflyKerl(lp + 1)%blocks(ip, j), P1%ButterflyKerl(lp)%blocks(ip, j)%matrix)
               call BFD_matset(O%ButterflyKerl(lp + 1)%blocks(2**lp + ip, j), P2%ButterflyKerl(lp)%blocks(ip, j)%matrix)
            enddo
         enddo
      enddo
      call BFD_setdims(O)
   end subroutine BFD_rowmerge

   subroutine BFD_matarray_free(W)
      implicit none
      type(butterflymatrix), allocatable::W(:)
      integer i
      if (.not. allocated(W)) return
      do i = 1, size(W)
         call BFD_matfree(W(i))
      enddo
      deallocate (W)
   end subroutine BFD_matarray_free

!======================================================================================== layouts of the distributed BMult

   !>**** number of owner bits of each tree (1 = a, 2 = b, 3 = c)
   subroutine BFD_lay_bits(lay, nbt)
      implicit none
      type(bfd_layout)::lay
      integer nbt(3), k
      nbt = 0
      do k = 1, lay%nseg
         nbt(lay%tr(k)) = lay%nb(k)
      enddo
   end subroutine BFD_lay_bits

   integer function BFD_lay_p(lay)
      implicit none
      type(bfd_layout)::lay
      BFD_lay_p = sum(lay%nb(1:lay%nseg))
   end function BFD_lay_p

   logical function BFD_lay_same(l1, l2)
      implicit none
      type(bfd_layout)::l1, l2
      BFD_lay_same = l1%nseg == l2%nseg
      if (BFD_lay_same .and. l1%nseg > 0) BFD_lay_same = all(l1%tr(1:l1%nseg) == l2%tr(1:l1%nseg)) .and. all(l1%nb(1:l1%nseg) == l2%nb(1:l1%nseg))
   end function BFD_lay_same

   !>**** owner (0-based chain rank) of node q (1-based tree coordinates) of state s
   integer function BFD_lay_owner(lay, s, q)
      implicit none
      type(bfd_layout)::lay
      integer s(3), q(3), k, t
      BFD_lay_owner = 0
      do k = 1, lay%nseg
         t = lay%tr(k)
         call assert(lay%nb(k) <= s(t), 'BFD_lay_owner: owner bits exceed the tree level')
         BFD_lay_owner = BFD_lay_owner*2**lay%nb(k) + ishft(q(t) - 1, -(s(t) - lay%nb(k)))
      enddo
   end function BFD_lay_owner

   !>**** local index (1-based) of node q of state s on its owner: the non-owner bits, tree a first
   integer function BFD_lay_lidx(lay, s, q)
      implicit none
      type(bfd_layout)::lay
      integer s(3), q(3), nbt(3), t, fr
      call BFD_lay_bits(lay, nbt)
      BFD_lay_lidx = 0
      do t = 1, 3
         fr = s(t) - nbt(t)
         BFD_lay_lidx = BFD_lay_lidx*2**fr + iand(q(t) - 1, 2**fr - 1)
      enddo
      BFD_lay_lidx = BFD_lay_lidx + 1
   end function BFD_lay_lidx

   !>**** coordinates q (1-based) of the node with local index li on chain rank v
   subroutine BFD_lay_node(lay, s, v, li, q)
      implicit none
      type(bfd_layout)::lay
      integer s(3), v, li, q(3), nbt(3), fixed(3), k, t, sh, x, fr
      call BFD_lay_bits(lay, nbt)
      fixed = 0
      sh = BFD_lay_p(lay)
      do k = 1, lay%nseg
         sh = sh - lay%nb(k)
         fixed(lay%tr(k)) = iand(ishft(v, -sh), 2**lay%nb(k) - 1)
      enddo
      x = li - 1
      do t = 3, 1, -1
         fr = s(t) - nbt(t)
         q(t) = fixed(t)*2**fr + mod(x, 2**fr) + 1
         x = x/2**fr
      enddo
   end subroutine BFD_lay_node

   !>**** number of local nodes of a state (every chain rank owns the same number)
   integer function BFD_lay_nloc(lay, s)
      implicit none
      type(bfd_layout)::lay
      integer s(3)
      BFD_lay_nloc = 2**(sum(s) - BFD_lay_p(lay))
   end function BFD_lay_nloc

   !>**** a factor with input state s refining tree inc and coarsening tree dec is local in the layout iff no owner bit is
   !> the bit its inc tree gains or its dec tree loses
   logical function BFD_lay_valid(lay, s, inc, dec)
      implicit none
      type(bfd_layout)::lay
      integer s(3), inc, dec, k, t
      BFD_lay_valid = .true.
      do k = 1, lay%nseg
         t = lay%tr(k)
         if (t == dec) then
            if (lay%nb(k) > s(t) - 1) BFD_lay_valid = .false.
         else
            if (lay%nb(k) > s(t)) BFD_lay_valid = .false.
         endif
      enddo
   end function BFD_lay_valid

   !>**** layout filling p owner bits in the tree order prio(1:3), each tree up to caps(tree)
   subroutine BFD_lay_canon(caps, p, prio, lay)
      implicit none
      integer caps(3), p, prio(3), k, left, nb
      type(bfd_layout)::lay
      lay%nseg = 0
      left = p
      do k = 1, 3
         nb = min(max(caps(prio(k)), 0), left)
         if (nb > 0) then
            lay%nseg = lay%nseg + 1
            lay%tr(lay%nseg) = prio(k)
            lay%nb(lay%nseg) = nb
            left = left - nb
         endif
      enddo
      call assert(left == 0, 'BFD_lay_canon: no valid layout')
   end subroutine BFD_lay_canon

!======================================================================================== BMult (3-tree chains)

   integer function BFD_nid(s, q)
      implicit none
      integer s(3), q(3)
      BFD_nid = ((q(1) - 1)*2**s(2) + (q(2) - 1))*2**s(3) + q(3)
   end function BFD_nid

   subroutine BFD_ncoord(s, id, q)
      implicit none
      integer s(3), id, q(3), r
      r = id - 1
      q(3) = mod(r, 2**s(3)) + 1
      r = r/2**s(3)
      q(2) = mod(r, 2**s(2)) + 1
      r = r/2**s(2)
      q(1) = r + 1
   end subroutine BFD_ncoord

   subroutine BFD_sout(f, so)
      implicit none
      type(bfd_fac)::f
      integer so(3)
      so = f%s
      so(f%inc) = so(f%inc) + 1
      so(f%dec) = so(f%dec) - 1
   end subroutine BFD_sout

   !>**** input node of f feeding output node o through slot d
   integer function BFD_fin(f, o, d)
      implicit none
      type(bfd_fac)::f
      integer o, d, so(3), q(3)
      call BFD_sout(f, so)
      call BFD_ncoord(so, o, q)
      q(f%inc) = (q(f%inc) + 1)/2
      q(f%dec) = 2*q(f%dec) - 2 + d
      BFD_fin = BFD_nid(f%s, q)
   end function BFD_fin

   !>**** the two output nodes of f fed by input node i, and the slot d they use
   subroutine BFD_fusers(f, i, o1, o2, d)
      implicit none
      type(bfd_fac)::f
      integer i, o1, o2, d, so(3), q(3), qo(3)
      call BFD_sout(f, so)
      call BFD_ncoord(f%s, i, q)
      qo = q
      qo(f%dec) = (q(f%dec) + 1)/2
      d = q(f%dec) - 2*(qo(f%dec) - 1)
      qo(f%inc) = 2*q(f%inc) - 1
      o1 = BFD_nid(so, qo)
      qo(f%inc) = 2*q(f%inc)
      o2 = BFD_nid(so, qo)
   end subroutine BFD_fusers

   subroutine BFD_facnew(f, s, inc, dec, nn)
      implicit none
      type(bfd_fac)::f
      integer s(3), inc, dec, nn
      f%s = s
      f%inc = inc
      f%dec = dec
      allocate (f%blk(2, nn))
   end subroutine BFD_facnew

   subroutine BFD_facfree(f)
      implicit none
      type(bfd_fac)::f
      integer o, d
      if (.not. allocated(f%blk)) return
      do o = 1, size(f%blk, 2)
         do d = 1, 2
            call BFD_matfree(f%blk(d, o))
         enddo
      enddo
      deallocate (f%blk)
   end subroutine BFD_facfree

   subroutine BFD_facmove(src, dst)
      implicit none
      type(bfd_fac)::src, dst
      call BFD_facfree(dst)
      dst%s = src%s
      dst%inc = src%inc
      dst%dec = src%dec
      call move_alloc(src%blk, dst%blk)
   end subroutine BFD_facmove

   subroutine BFD_ch_free(ch)
      implicit none
      type(bfd_chain)::ch
      integer i
      if (allocated(ch%U)) then
         do i = 1, size(ch%U)
            call BFD_matfree(ch%U(i))
         enddo
         deallocate (ch%U)
      endif
      if (allocated(ch%V)) then
         do i = 1, size(ch%V)
            call BFD_matfree(ch%V(i))
         enddo
         deallocate (ch%V)
      endif
      if (allocated(ch%f)) then
         do i = 1, size(ch%f)
            call BFD_facfree(ch%f(i))
         enddo
         deallocate (ch%f)
      endif
   end subroutine BFD_ch_free

   !>**** chain U^X * X-kernels ('ba') * M ('ca') * Y-kernels ('cb') * V^Y, with M = R^X_1 * V^X * U^Y * R^Y_K (paper step 2)
   subroutine BFD_ch_build(X, Y, ch, stats)
      implicit none
      type(matrixblock)::X, Y
      type(bfd_chain)::ch
      type(Hstat)::stats
      type(butterflymatrix)::T1, T2
      integer K, nn, p, l, q, i, j, d, e, o, a1, u, jm
      K = X%level_butterfly
      nn = 2**K
      ch%K = K
      ch%nf = 2*K - 1
      allocate (ch%U(nn), ch%V(nn), ch%f(2*K - 1))
      do i = 1, nn
         call BFD_matset(ch%U(i), X%ButterflyU%blocks(i)%matrix)
         call BFD_matset(ch%V(i), Y%ButterflyV%blocks(i)%matrix)
      enddo
      do p = 1, K - 1
         l = K + 1 - p
         call BFD_facnew(ch%f(p), [l - 1, K - l + 1, 0], 1, 2, nn)
         do i = 1, 2**l
            do j = 1, 2**(K - l)
               o = BFD_nid([l, K - l, 0], [i, j, 1])
               do d = 1, 2
                  call BFD_matset(ch%f(p)%blk(d, o), X%ButterflyKerl(l)%blocks(i, 2*j - 2 + d)%matrix)
               enddo
            enddo
         enddo
      enddo
      call BFD_facnew(ch%f(K), [0, K - 1, 1], 1, 3, nn)
      do a1 = 1, 2
         do u = 1, 2**(K - 1)
            o = BFD_nid([1, K - 1, 0], [a1, u, 1])
            do d = 1, 2
               do e = 1, 2
                  jm = 2*u - 2 + e
                  call BFD_mulnew('T', X%ButterflyV%blocks(jm)%matrix, 'N', Y%ButterflyU%blocks(jm)%matrix, T1, stats)
                  call BFD_mulnew('N', T1%matrix, 'N', Y%ButterflyKerl(K)%blocks(jm, d)%matrix, T2, stats)
                  if (e == 1) then
                     call BFD_mulnew('N', X%ButterflyKerl(1)%blocks(a1, jm)%matrix, 'N', T2%matrix, ch%f(K)%blk(d, o), stats)
                  else
                     call BFD_gemm('N', 'N', X%ButterflyKerl(1)%blocks(a1, jm)%matrix, T2%matrix, ch%f(K)%blk(d, o)%matrix, BPACK_cone, BPACK_cone, stats)
                  endif
               enddo
            enddo
         enddo
      enddo
      call BFD_matfree(T1)
      call BFD_matfree(T2)
      do q = 1, K - 1
         l = K - q
         call BFD_facnew(ch%f(K + q), [0, l - 1, K - l + 1], 2, 3, nn)
         do i = 1, 2**l
            do j = 1, 2**(K - l)
               o = BFD_nid([0, l, K - l], [1, i, j])
               do d = 1, 2
                  call BFD_matset(ch%f(K + q)%blk(d, o), Y%ButterflyKerl(l)%blocks(i, 2*j - 2 + d)%matrix)
               enddo
            enddo
         enddo
      enddo
   end subroutine BFD_ch_build

   !>**** positions lo..hi (0 = U, 1..nf = factors): column-orthonormalize (truncated SVD of every column block) and push
   !> S*Z into the next position
   subroutine BFD_ch_rsweep(ch, lo, hi, tol, stats)
      implicit none
      type(bfd_chain)::ch
      integer lo, hi
      real(kind=8)::tol
      type(Hstat)::stats
      DT, allocatable::W(:, :), Z(:, :), T(:, :)
      DTR, allocatable::S(:)
      type(butterflymatrix), allocatable::Rc(:)
      integer nn, p, it, i, j, r, o1, o2, d, h1, h2, dd
      nn = 2**ch%K
      allocate (Rc(nn))
      do p = lo, hi
         if (p == 0) then
            do it = 1, nn
               call BFD_qr(ch%U(it)%matrix, W, Z, stats)
               call BFD_matset(ch%U(it), W)
               do dd = 1, 2
                  call BFD_lmul(Z, ch%f(1)%blk(dd, it), stats)
               enddo
            enddo
         else
            do i = 1, nn
               call BFD_fusers(ch%f(p), i, o1, o2, d)
               h1 = size(ch%f(p)%blk(d, o1)%matrix, 1)
               h2 = size(ch%f(p)%blk(d, o2)%matrix, 1)
               allocate (T(h1 + h2, size(ch%f(p)%blk(d, o1)%matrix, 2)))
               T(1:h1, :) = ch%f(p)%blk(d, o1)%matrix
               T(h1 + 1:, :) = ch%f(p)%blk(d, o2)%matrix
               call BFD_qr(T, W, Z, stats)
               call BFD_matset(ch%f(p)%blk(d, o1), W(1:h1, :))
               call BFD_matset(ch%f(p)%blk(d, o2), W(h1 + 1:h1 + h2, :))
               call BFD_matset(Rc(i), Z)
               deallocate (T)
            enddo
            if (p < ch%nf) then
               do i = 1, nn
                  do dd = 1, 2
                     call BFD_lmul(Rc(i)%matrix, ch%f(p + 1)%blk(dd, i), stats)
                  enddo
               enddo
            else
               do i = 1, nn
                  call BFD_rmul(ch%V(i), transpose(Rc(i)%matrix), stats)
               enddo
            endif
         endif
      enddo
      call BFD_matarray_free(Rc)
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
   end subroutine BFD_ch_rsweep

   !>**** positions hi down to lo (1..nf = factors, nf+1 = V): row-orthonormalize (truncated SVD of every row block) and push
   !> W*S into the previous position
   subroutine BFD_ch_lsweep(ch, lo, hi, tol, stats)
      implicit none
      type(bfd_chain)::ch
      integer lo, hi
      real(kind=8)::tol
      type(Hstat)::stats
      DT, allocatable::W(:, :), Z(:, :), T(:, :)
      DTR, allocatable::S(:)
      type(butterflymatrix), allocatable::L(:)
      integer nn, p, v, o, j, r, o1, o2, d, n1, n2
      nn = 2**ch%K
      allocate (L(nn))
      do p = hi, lo, -1
         if (p == ch%nf + 1) then
            do v = 1, nn
               call BFD_qr(ch%V(v)%matrix, W, Z, stats)
               call BFD_matset(ch%V(v), W)
               call BFD_fusers(ch%f(ch%nf), v, o1, o2, d)
               call BFD_rmul(ch%f(ch%nf)%blk(d, o1), transpose(Z), stats)
               call BFD_rmul(ch%f(ch%nf)%blk(d, o2), transpose(Z), stats)
            enddo
         else
            do o = 1, nn
               n1 = size(ch%f(p)%blk(1, o)%matrix, 2)
               n2 = size(ch%f(p)%blk(2, o)%matrix, 2)
               allocate (T(size(ch%f(p)%blk(1, o)%matrix, 1), n1 + n2))
               T(:, 1:n1) = ch%f(p)%blk(1, o)%matrix
               T(:, n1 + 1:) = ch%f(p)%blk(2, o)%matrix
               call BFD_qr(transpose(T), W, Z, stats)
               call BFD_matset(ch%f(p)%blk(1, o), transpose(W(1:n1, :)))
               call BFD_matset(ch%f(p)%blk(2, o), transpose(W(n1 + 1:n1 + n2, :)))
               call BFD_matset(L(o), transpose(Z))
               deallocate (T)
            enddo
            if (p > 1) then
               do o = 1, nn
                  call BFD_fusers(ch%f(p - 1), o, o1, o2, d)
                  call BFD_rmul(ch%f(p - 1)%blk(d, o1), L(o)%matrix, stats)
                  call BFD_rmul(ch%f(p - 1)%blk(d, o2), L(o)%matrix, stats)
               enddo
            else
               do o = 1, nn
                  call BFD_rmul(ch%U(o), L(o)%matrix, stats)
               enddo
            endif
         endif
      enddo
      call BFD_matarray_free(L)
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
   end subroutine BFD_ch_lsweep

   !>**** truncate the input state of f(p) (between f(p) and f(p+1)): SVD of the column blocks of f(p), keep W*S in f(p)
   !> and multiply Z into the rows of f(p+1)
   subroutine BFD_ch_ctrunc(ch, p, tol, stats)
      implicit none
      type(bfd_chain)::ch
      integer p
      real(kind=8)::tol
      type(Hstat)::stats
      DT, allocatable::W(:, :), Z(:, :), T(:, :)
      DTR, allocatable::S(:)
      integer nn, m, o1, o2, d, h1, h2, r, j, e
      nn = 2**ch%K
      do m = 1, nn
         call BFD_fusers(ch%f(p), m, o1, o2, d)
         h1 = size(ch%f(p)%blk(d, o1)%matrix, 1)
         h2 = size(ch%f(p)%blk(d, o2)%matrix, 1)
         allocate (T(h1 + h2, size(ch%f(p)%blk(d, o1)%matrix, 2)))
         T(1:h1, :) = ch%f(p)%blk(d, o1)%matrix
         T(h1 + 1:, :) = ch%f(p)%blk(d, o2)%matrix
         call BFD_svd(T, W, S, Z, r, tol, .true., stats, ch%scale)
         do j = 1, r
            W(:, j) = W(:, j)*S(j)
         enddo
         call BFD_matset(ch%f(p)%blk(d, o1), W(1:h1, 1:r))
         call BFD_matset(ch%f(p)%blk(d, o2), W(h1 + 1:h1 + h2, 1:r))
         do e = 1, 2
            call BFD_lmul(Z(1:r, :), ch%f(p + 1)%blk(e, m), stats)
         enddo
         deallocate (T)
      enddo
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
   end subroutine BFD_ch_ctrunc

   !>**** BRowSwap: f(p) ('ca') * f(p+1) ('cb') = B ('cb') * C ('ca') exactly, C = stacked identities (copies of the two
   !> sibling inputs); then the new state is truncated
   subroutine BFD_ch_swap(ch, p, tol, stats)
      implicit none
      type(bfd_chain)::ch
      integer p
      real(kind=8)::tol
      type(Hstat)::stats
      type(bfd_fac)::Bf, Cf
      type(butterflymatrix)::T
      integer nn, s0(3), smid(3), o, d, mo, m, i1, i2, n1, n2, o1, o2, dd
      nn = 2**ch%K
      s0 = ch%f(p + 1)%s
      smid = s0
      smid(1) = smid(1) + 1
      smid(3) = smid(3) - 1
      call BFD_facnew(Bf, smid, 2, 3, nn)
      call BFD_facnew(Cf, s0, 1, 3, nn)
      do o = 1, nn
         do d = 1, 2
            mo = BFD_fin(ch%f(p), o, d)
            call BFD_hcat(ch%f(p + 1)%blk(1, mo)%matrix, ch%f(p + 1)%blk(2, mo)%matrix, BPACK_cone, T)
            call BFD_mulnew('N', ch%f(p)%blk(d, o)%matrix, 'N', T%matrix, Bf%blk(d, o), stats)
         enddo
      enddo
      call BFD_matfree(T)
      do m = 1, nn
         i1 = BFD_fin(Cf, m, 1)
         i2 = BFD_fin(Cf, m, 2)
         call BFD_fusers(ch%f(p + 1), i1, o1, o2, dd)
         n1 = size(ch%f(p + 1)%blk(dd, o1)%matrix, 2)
         call BFD_fusers(ch%f(p + 1), i2, o1, o2, dd)
         n2 = size(ch%f(p + 1)%blk(dd, o1)%matrix, 2)
         call BFD_matnew(Cf%blk(1, m), n1 + n2, n1)
         call BFD_matnew(Cf%blk(2, m), n1 + n2, n2)
         do dd = 1, n1
            Cf%blk(1, m)%matrix(dd, dd) = BPACK_cone
         enddo
         do dd = 1, n2
            Cf%blk(2, m)%matrix(n1 + dd, dd) = BPACK_cone
         enddo
      enddo
      call BFD_facmove(Bf, ch%f(p))
      call BFD_facmove(Cf, ch%f(p + 1))
      call BFD_ch_ctrunc(ch, p, tol, stats)
   end subroutine BFD_ch_swap

   !>**** f(p) ('ba') * f(p+1) ('cb') = one 'ca' factor (the middle tree is contracted); the chain gets one factor shorter
   subroutine BFD_ch_merge(ch, p, stats)
      implicit none
      type(bfd_chain)::ch
      integer p
      type(Hstat)::stats
      type(bfd_fac)::Pf
      integer nn, o, e, m1, m2, q
      nn = 2**ch%K
      call BFD_facnew(Pf, ch%f(p + 1)%s, 1, 3, nn)
      do o = 1, nn
         m1 = BFD_fin(ch%f(p), o, 1)
         m2 = BFD_fin(ch%f(p), o, 2)
         do e = 1, 2
            call BFD_mulnew('N', ch%f(p)%blk(1, o)%matrix, 'N', ch%f(p + 1)%blk(e, m1)%matrix, Pf%blk(e, o), stats)
            call BFD_gemm('N', 'N', ch%f(p)%blk(2, o)%matrix, ch%f(p + 1)%blk(e, m2)%matrix, Pf%blk(e, o)%matrix, BPACK_cone, BPACK_cone, stats)
         enddo
      enddo
      call BFD_facmove(Pf, ch%f(p))
      call BFD_facfree(ch%f(p + 1))
      do q = p + 1, ch%nf - 1
         call BFD_facmove(ch%f(q + 1), ch%f(q))
      enddo
      ch%nf = ch%nf - 1
   end subroutine BFD_ch_merge

   !>**** chain U * K 'ca' factors * V -> K-level butterfly
   subroutine BFD_ch_extract(ch, Z, pgno)
      implicit none
      type(bfd_chain)::ch
      type(matrixblock)::Z
      integer pgno, K, p, l, i, jo, o, d
      K = ch%K
      call assert(ch%nf == K, 'BFD_ch_extract: chain is not a butterfly')
      call BFD_alloc(Z, K, pgno)
      do i = 1, 2**K
         call BFD_matmove(ch%U(i), Z%ButterflyU%blocks(i))
         call BFD_matmove(ch%V(i), Z%ButterflyV%blocks(i))
      enddo
      do p = 1, K
         l = K + 1 - p
         call assert(ch%f(p)%s(1) == l - 1 .and. ch%f(p)%s(2) == 0 .and. ch%f(p)%inc == 1 .and. ch%f(p)%dec == 3, 'BFD_ch_extract: unexpected factor')
         do i = 1, 2**l
            do jo = 1, 2**(K - l)
               o = BFD_nid([l, 0, K - l], [i, 1, jo])
               do d = 1, 2
                  call BFD_matmove(ch%f(p)%blk(d, o), Z%ButterflyKerl(l)%blocks(i, 2*jo - 2 + d))
               enddo
            enddo
         enddo
      enddo
      call BFD_setdims(Z)
   end subroutine BFD_ch_extract

   !>**** Z = X*Y for two K-level butterflies sharing the middle leaves (paper IV-C). After every BRowSwap the new state is
   !> truncated while the chain is kept orthonormal around the swap position, so all intermediate ranks stay at their true values.
   !> All truncations are relative to max(largest singular value, scale) (BF_TruncRank)
   subroutine BFD_bmult(X, Y, Z, pgno, tol, option, stats, ptree, scale)
      implicit none
      type(matrixblock)::X, Y, Z
      integer pgno
      real(kind=8), optional::scale
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(bfd_chain)::ch
      type(butterflymatrix)::T
      integer K, m, it, oc, l
      real(kind=8)::tol, f0, t0
      DT, allocatable::Dx(:, :), Dy(:, :), Dz(:, :)
      K = X%level_butterfly
      call assert(Y%level_butterfly == K, 'BFD_bmult: the two butterflies must have the same level')
      call assert(X%N == Y%M, 'BFD_bmult: inner dimensions do not match')
      if (K == 0) then
         call BFD_alloc(Z, 0, pgno)
         call BFD_mulnew('T', X%ButterflyV%blocks(1)%matrix, 'N', Y%ButterflyU%blocks(1)%matrix, T, stats)
         call BFD_mulnew('N', X%ButterflyU%blocks(1)%matrix, 'N', T%matrix, Z%ButterflyU%blocks(1), stats)
         call BFD_matset(Z%ButterflyV%blocks(1), Y%ButterflyV%blocks(1)%matrix)
         call BFD_matfree(T)
         call BFD_setdims(Z)
         call BFD_recompress(Z, tol, option, stats, ptree, scale)
         return
      endif
      oc = 7
      if (BFD_rec_depth > 0) oc = 14
      call BFD_mstart(stats, f0, t0)
      call BFD_ch_build(X, Y, ch, stats)
      if (present(scale)) ch%scale = scale
      call BFD_mstat(oc + 1, stats, f0, t0)
      call BFD_mstart(stats, f0, t0)
      call BFD_ch_lsweep(ch, K + 1, 2*K, tol, stats)
      call BFD_mstat(oc + 2, stats, f0, t0)
      do m = 1, K - 1
         call BFD_mstart(stats, f0, t0)
         call BFD_ch_rsweep(ch, 0, K - 1, tol, stats)
         call BFD_mstat(oc + 3, stats, f0, t0)
         call BFD_mstart(stats, f0, t0)
         do it = 1, m
            call BFD_ch_swap(ch, K + 1 - it, tol, stats)
         enddo
         call BFD_mstat(oc + 4, stats, f0, t0)
         if (BFD_mul_level >= 0) then
            l = min(BFD_mul_level, BFD_RD - 1)
            BFD_mul_rank(l, oc/7) = max(BFD_mul_rank(l, oc/7), BFD_ch_maxdim(ch))
         endif
         call BFD_mstart(stats, f0, t0)
         call BFD_ch_merge(ch, K - m, stats)
         call BFD_mstat(oc + 5, stats, f0, t0)
      enddo
      call BFD_mstart(stats, f0, t0)
      call BFD_ch_extract(ch, Z, pgno)
      call BFD_ch_free(ch)
      call BFD_mstat(oc + 6, stats, f0, t0)
      if (BFD_debug(option)) then
         call BFD_todense(X, Dx, stats)
         call BFD_todense(Y, Dy, stats)
         call BFD_todense(Z, Dz, stats)
         call BFD_report('bmult(before recompress)', Dz, matmul(Dx, Dy), K)
      endif
      call BFD_mstart(stats, f0, t0)
      call BFD_recompress(Z, tol, option, stats, ptree, scale)
      call BFD_mstat(oc + 7, stats, f0, t0)
      if (BFD_debug(option)) then
         call BFD_todense(Z, Dz, stats)
         call BFD_report('bmult', Dz, matmul(Dx, Dy), K)
      endif
   end subroutine BFD_bmult

!======================================================================================== distributed BMult

   !>**** append an item (a copy of A) to the list L
   subroutine BFD_items_add(L, dst, h1, h2, h3, h4, A)
      implicit none
      type(bfd_items)::L
      integer dst, h1, h2, h3, h4, k
      DT::A(:, :)
      type(bfd_item), allocatable::tmp(:)
      if (.not. allocated(L%it)) allocate (L%it(64))
      if (L%n == size(L%it)) then
         allocate (tmp(2*size(L%it)))
         do k = 1, L%n
            tmp(k)%dst = L%it(k)%dst
            tmp(k)%hdr = L%it(k)%hdr
            call move_alloc(L%it(k)%mat, tmp(k)%mat)
         enddo
         call move_alloc(tmp, L%it)
      endif
      L%n = L%n + 1
      L%it(L%n)%dst = dst
      L%it(L%n)%hdr = [h1, h2, h3, h4]
      allocate (L%it(L%n)%mat(size(A, 1), size(A, 2)))
      if (size(A) > 0) L%it(L%n)%mat = A
   end subroutine BFD_items_add

   !>**** append an item to L by moving its matrix (the item is emptied)
   subroutine BFD_items_take(L, item)
      implicit none
      type(bfd_items)::L
      type(bfd_item)::item
      type(bfd_item), allocatable::tmp(:)
      integer k
      if (.not. allocated(L%it)) allocate (L%it(64))
      if (L%n == size(L%it)) then
         allocate (tmp(2*size(L%it)))
         do k = 1, L%n
            tmp(k)%dst = L%it(k)%dst
            tmp(k)%hdr = L%it(k)%hdr
            call move_alloc(L%it(k)%mat, tmp(k)%mat)
         enddo
         call move_alloc(tmp, L%it)
      endif
      L%n = L%n + 1
      L%it(L%n)%dst = -1
      L%it(L%n)%hdr = item%hdr
      call move_alloc(item%mat, L%it(L%n)%mat)
   end subroutine BFD_items_take

   subroutine BFD_items_free(L)
      implicit none
      type(bfd_items)::L
      if (allocated(L%it)) deallocate (L%it)
      L%n = 0
   end subroutine BFD_items_free

   !>**** sparse all-to-all over the communicator comm: the items of Lout go to their destinations and are appended to Lin
   !> (collective: every rank of comm calls it, possibly with no items)
   subroutine BFD_exchange(Lout, Lin, comm)
      implicit none
      type(bfd_items)::Lout, Lin
      integer comm, np, ierr, k, dst, m, n, pi, pd, ni, nd, me
      integer, allocatable::sci(:), rci(:), sdi(:), rdi(:), scd(:), rcd(:), sdd(:), rdd(:), ibuf(:), ribuf(:), posi(:), posd(:)
      DT, allocatable::dbuf(:), rdbuf(:)

      call MPI_Comm_size(comm, np, ierr)
      call MPI_Comm_rank(comm, me, ierr)
      if (np == 1) then ! everything is local: move the items without communication
         do k = 1, Lout%n
            call BFD_items_take(Lin, Lout%it(k))
         enddo
         call BFD_items_free(Lout)
         return
      endif
      allocate (sci(0:np - 1), rci(0:np - 1), sdi(0:np - 1), rdi(0:np - 1), scd(0:np - 1), rcd(0:np - 1), sdd(0:np - 1), rdd(0:np - 1), posi(0:np - 1), posd(0:np - 1))
      sci = 0
      scd = 0
      do k = 1, Lout%n
         dst = Lout%it(k)%dst
         if (dst == me) cycle
         sci(dst) = sci(dst) + 6
         scd(dst) = scd(dst) + size(Lout%it(k)%mat)
      enddo
      call MPI_Alltoall(sci, 1, MPI_INTEGER, rci, 1, MPI_INTEGER, comm, ierr)
      call MPI_Alltoall(scd, 1, MPI_INTEGER, rcd, 1, MPI_INTEGER, comm, ierr)
      sdi(0) = 0
      sdd(0) = 0
      rdi(0) = 0
      rdd(0) = 0
      do k = 1, np - 1
         sdi(k) = sdi(k - 1) + sci(k - 1)
         sdd(k) = sdd(k - 1) + scd(k - 1)
         rdi(k) = rdi(k - 1) + rci(k - 1)
         rdd(k) = rdd(k - 1) + rcd(k - 1)
      enddo
      allocate (ibuf(max(1, sum(sci))), dbuf(max(1, sum(scd))), ribuf(max(1, sum(rci))), rdbuf(max(1, sum(rcd))))
      posi = sdi
      posd = sdd
      do k = 1, Lout%n
         dst = Lout%it(k)%dst
         if (dst == me) then ! kept local: moved to Lin without packing
            call BFD_items_take(Lin, Lout%it(k))
            cycle
         endif
         m = size(Lout%it(k)%mat, 1)
         n = size(Lout%it(k)%mat, 2)
         ibuf(posi(dst) + 1:posi(dst) + 4) = Lout%it(k)%hdr
         ibuf(posi(dst) + 5) = m
         ibuf(posi(dst) + 6) = n
         posi(dst) = posi(dst) + 6
         if (m*n > 0) dbuf(posd(dst) + 1:posd(dst) + m*n) = reshape(Lout%it(k)%mat, [m*n])
         posd(dst) = posd(dst) + m*n
      enddo
      call BFD_items_free(Lout)
      call MPI_Alltoallv(ibuf, sci, sdi, MPI_INTEGER, ribuf, rci, rdi, MPI_INTEGER, comm, ierr)
      call MPI_Alltoallv(dbuf, scd, sdd, MPI_DT, rdbuf, rcd, rdd, MPI_DT, comm, ierr)
      ni = sum(rci)
      nd = sum(rcd)
      pi = 0
      pd = 0
      do while (pi < ni)
         m = ribuf(pi + 5)
         n = ribuf(pi + 6)
         call BFD_items_add(Lin, -1, ribuf(pi + 1), ribuf(pi + 2), ribuf(pi + 3), ribuf(pi + 4), reshape(rdbuf(pd + 1:pd + m*n), [m, n]))
         pi = pi + 6
         pd = pd + m*n
      enddo
      deallocate (sci, rci, sdi, rdi, scd, rcd, sdd, rdd, posi, posd, ibuf, dbuf, ribuf, rdbuf)
   end subroutine BFD_exchange

   !>**** output state of a factor (input state s, refining inc, coarsening dec)
   subroutine BFD_q_sout(s, inc, dec, so)
      implicit none
      integer s(3), inc, dec, so(3)
      so = s
      so(inc) = so(inc) + 1
      so(dec) = so(dec) - 1
   end subroutine BFD_q_sout

   !>**** input node qi of a factor feeding its output node qo through slot d
   subroutine BFD_q_fin(inc, dec, qo, d, qi)
      implicit none
      integer inc, dec, qo(3), d, qi(3)
      qi = qo
      qi(inc) = (qo(inc) + 1)/2
      qi(dec) = 2*qo(dec) - 2 + d
   end subroutine BFD_q_fin

   !>**** the two output nodes of a factor fed by its input node qi, and the slot they use
   subroutine BFD_q_fusers(inc, dec, qi, qo1, qo2, d)
      implicit none
      integer inc, dec, qi(3), qo1(3), qo2(3), d
      qo1 = qi
      qo1(dec) = (qi(dec) + 1)/2
      d = qi(dec) - 2*(qo1(dec) - 1)
      qo2 = qo1
      qo1(inc) = 2*qi(inc) - 1
      qo2(inc) = 2*qi(inc)
   end subroutine BFD_q_fusers

   !>**** canonical layout of a factor (or of a pair: caps are the minimum over both): 'ca' factors keep their a and c
   !> levels while the chain evolves, so they fill a and c first (larger first); 'ba' fill a first, 'cb' c first
   subroutine BFD_dch_canon(caps, p, kind, lay)
      implicit none
      integer caps(3), p, kind
      type(bfd_layout)::lay
      if (kind == 1) then ! 'ba'
         call BFD_lay_canon(caps, p, [1, 2, 3], lay)
      elseif (kind == 3) then ! 'cb'
         call BFD_lay_canon(caps, p, [3, 2, 1], lay)
      elseif (caps(1) >= caps(3)) then ! 'ca'
         call BFD_lay_canon(caps, p, [1, 3, 2], lay)
      else
         call BFD_lay_canon(caps, p, [3, 1, 2], lay)
      endif
   end subroutine BFD_dch_canon

   subroutine BFD_dfac_caps(F, caps)
      implicit none
      type(bfd_dfac)::F
      integer caps(3), t
      do t = 1, 3
         caps(t) = F%s(t)
         if (t == F%dec) caps(t) = caps(t) - 1
      enddo
   end subroutine BFD_dfac_caps

   integer function BFD_dfac_kind(F)
      implicit none
      type(bfd_dfac)::F
      if (F%dec == 2) then
         BFD_dfac_kind = 1 ! 'ba'
      elseif (F%inc == 2) then
         BFD_dfac_kind = 3 ! 'cb'
      else
         BFD_dfac_kind = 2 ! 'ca'
      endif
   end function BFD_dfac_kind

   subroutine BFD_dfac_alloc(dc, F)
      implicit none
      type(bfd_dchain)::dc
      type(bfd_dfac)::F
      if (allocated(F%blk)) deallocate (F%blk)
      if (dc%myv >= 0) allocate (F%blk(2, BFD_lay_nloc(F%lay, F%s)))
   end subroutine BFD_dfac_alloc

   subroutine BFD_dfac_free(F)
      implicit none
      type(bfd_dfac)::F
      integer i, d
      if (allocated(F%blk)) then
         do i = 1, size(F%blk, 2)
            do d = 1, 2
               call BFD_matfree(F%blk(d, i))
            enddo
         enddo
         deallocate (F%blk)
      endif
   end subroutine BFD_dfac_free

   subroutine BFD_dfac_move(src, dst)
      implicit none
      type(bfd_dfac)::src, dst
      call BFD_dfac_free(dst)
      dst%s = src%s
      dst%inc = src%inc
      dst%dec = src%dec
      dst%lay = src%lay
      if (allocated(src%blk)) call move_alloc(src%blk, dst%blk)
   end subroutine BFD_dfac_move

   subroutine BFD_dch_free(dc)
      implicit none
      type(bfd_dchain)::dc
      integer i
      if (allocated(dc%f)) then
         do i = 1, size(dc%f)
            call BFD_dfac_free(dc%f(i))
         enddo
         deallocate (dc%f)
      endif
      if (allocated(dc%U)) call BFD_matarray_free(dc%U)
      if (allocated(dc%V)) call BFD_matarray_free(dc%V)
   end subroutine BFD_dch_free

   !>**** the chain U^X * ba^(K-1) * M * cb^(K-1) * V^Y (as BFD_ch_build) distributed over the chain ranks: layouts are
   !> inherited from the left neighbour when valid, canonical otherwise; the blocks of X and Y are sent from their
   !> ButterflyPACK layouts in one exchange, and M is formed where its middle leaves arrive.
   subroutine BFD_dch_build(X, Y, dc, ptree, stats)
      implicit none
      type(matrixblock)::X, Y
      type(bfd_dchain)::dc
      type(proctree)::ptree
      type(Hstat)::stats
      type(bfd_items)::Lout, Lin
      type(butterflymatrix), allocatable::Xk1(:, :), VX(:), UY(:), YkK(:, :)
      type(butterflymatrix)::T1, T2
      integer K, nf, pphys, nproc, me, ierr, q, l, ii, jj, i, j, d, e, jm, li, kk, t, caps(3), so(3), qo(3), nbt(3), soM(3), dst, kind
      logical valid

      K = X%level_butterfly
      call assert(Y%level_butterfly == K .and. Y%pgno == X%pgno, 'BFD_dch_build: X and Y must have the same levels and process group')
      dc%K = K
      dc%pgno = X%pgno
      dc%comm = ptree%pgrp(X%pgno)%Comm
      call MPI_Comm_rank(dc%comm, me, ierr)
      call MPI_Comm_size(dc%comm, nproc, ierr)
      pphys = nint(log(dble(nproc))/log(2d0))
      call assert(2**pphys == nproc, 'BFD_dch_build: the process group must have a power-of-2 size')
      call assert(pphys <= K, 'BFD_dch_build: more processes than leaves')
      dc%p = max(0, min(pphys, K - 2))
      dc%stride = 2**(pphys - dc%p)
      dc%myv = -1
      if (mod(me, dc%stride) == 0) dc%myv = me/dc%stride
      nf = 2*K - 1
      dc%nf = nf
      allocate (dc%f(nf))
      do q = 1, K - 1
         l = K + 1 - q
         dc%f(q)%s = [l - 1, K - l + 1, 0]
         dc%f(q)%inc = 1
         dc%f(q)%dec = 2
      enddo
      dc%f(K)%s = [0, K - 1, 1]
      dc%f(K)%inc = 1
      dc%f(K)%dec = 3
      do q = 1, K - 1
         l = K - q
         dc%f(K + q)%s = [0, l - 1, K - l + 1]
         dc%f(K + q)%inc = 2
         dc%f(K + q)%dec = 3
      enddo
      do q = 1, nf
         valid = .false.
         if (q > 1) valid = BFD_lay_valid(dc%f(q - 1)%lay, dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec)
         if (valid) then
            dc%f(q)%lay = dc%f(q - 1)%lay
         else
            call BFD_dfac_caps(dc%f(q), caps)
            call BFD_dch_canon(caps, dc%p, BFD_dfac_kind(dc%f(q)), dc%f(q)%lay)
         endif
         call BFD_dfac_alloc(dc, dc%f(q))
      enddo
      dc%Ulay = dc%f(1)%lay
      dc%Vlay = dc%f(nf)%lay
      if (dc%myv >= 0) allocate (dc%U(2**(K - dc%p)), dc%V(2**(K - dc%p)))
      soM = [1, K - 1, 0]
      call BFD_lay_bits(dc%f(K)%lay, nbt)
      call assert(nbt(1) == 0, 'BFD_dch_build: M must not be distributed by the X-row tree')

      !>**** kernels of X (levels 2..K) and Y (levels 1..K-1)
      do l = 1, K
         do kk = 1, 2
            if (kk == 1 .and. l == 1) cycle
            if (kk == 2 .and. l == K) cycle
            if (kk == 1) then
               if (.not. allocated(X%ButterflyKerl(l)%blocks)) cycle
               do jj = 1, X%ButterflyKerl(l)%nc
                  do ii = 1, X%ButterflyKerl(l)%nr
                     if (.not. associated(X%ButterflyKerl(l)%blocks(ii, jj)%matrix)) cycle
                     i = (ii - 1)*X%ButterflyKerl(l)%inc_r + X%ButterflyKerl(l)%idx_r
                     j = (jj - 1)*X%ButterflyKerl(l)%inc_c + X%ButterflyKerl(l)%idx_c
                     q = K + 1 - l
                     so = [l, K - l, 0]
                     qo = [i, (j + 1)/2, 1]
                     d = j - 2*((j + 1)/2) + 2
                     dst = BFD_lay_owner(dc%f(q)%lay, so, qo)*dc%stride
                     call BFD_items_add(Lout, dst, q, BFD_nid(so, qo), d, 0, X%ButterflyKerl(l)%blocks(ii, jj)%matrix)
                  enddo
               enddo
            else
               if (.not. allocated(Y%ButterflyKerl(l)%blocks)) cycle
               do jj = 1, Y%ButterflyKerl(l)%nc
                  do ii = 1, Y%ButterflyKerl(l)%nr
                     if (.not. associated(Y%ButterflyKerl(l)%blocks(ii, jj)%matrix)) cycle
                     i = (ii - 1)*Y%ButterflyKerl(l)%inc_r + Y%ButterflyKerl(l)%idx_r
                     j = (jj - 1)*Y%ButterflyKerl(l)%inc_c + Y%ButterflyKerl(l)%idx_c
                     q = 2*K - l
                     so = [0, l, K - l]
                     qo = [1, i, (j + 1)/2]
                     d = j - 2*((j + 1)/2) + 2
                     dst = BFD_lay_owner(dc%f(q)%lay, so, qo)*dc%stride
                     call BFD_items_add(Lout, dst, q, BFD_nid(so, qo), d, 0, Y%ButterflyKerl(l)%blocks(ii, jj)%matrix)
                  enddo
               enddo
            endif
         enddo
      enddo
      !>**** the pieces of M = X ker 1 * V^X U^Y * Y ker K, sent to the owner of the middle node (jm+1)/2
      if (allocated(X%ButterflyKerl(1)%blocks)) then
         do jj = 1, X%ButterflyKerl(1)%nc
            do ii = 1, X%ButterflyKerl(1)%nr
               if (.not. associated(X%ButterflyKerl(1)%blocks(ii, jj)%matrix)) cycle
               i = (ii - 1)*X%ButterflyKerl(1)%inc_r + X%ButterflyKerl(1)%idx_r
               jm = (jj - 1)*X%ButterflyKerl(1)%inc_c + X%ButterflyKerl(1)%idx_c
               dst = BFD_lay_owner(dc%f(K)%lay, soM, [i, (jm + 1)/2, 1])*dc%stride
               call BFD_items_add(Lout, dst, -1, i, jm, 0, X%ButterflyKerl(1)%blocks(ii, jj)%matrix)
            enddo
         enddo
      endif
      if (allocated(Y%ButterflyKerl(K)%blocks)) then
         do jj = 1, Y%ButterflyKerl(K)%nc
            do ii = 1, Y%ButterflyKerl(K)%nr
               if (.not. associated(Y%ButterflyKerl(K)%blocks(ii, jj)%matrix)) cycle
               jm = (ii - 1)*Y%ButterflyKerl(K)%inc_r + Y%ButterflyKerl(K)%idx_r
               d = (jj - 1)*Y%ButterflyKerl(K)%inc_c + Y%ButterflyKerl(K)%idx_c
               dst = BFD_lay_owner(dc%f(K)%lay, soM, [1, (jm + 1)/2, 1])*dc%stride
               call BFD_items_add(Lout, dst, -4, jm, d, 0, Y%ButterflyKerl(K)%blocks(ii, jj)%matrix)
            enddo
         enddo
      endif
      do kk = 1, X%ButterflyV%nblk_loc
         jm = (kk - 1)*X%ButterflyV%inc + X%ButterflyV%idx
         dst = BFD_lay_owner(dc%f(K)%lay, soM, [1, (jm + 1)/2, 1])*dc%stride
         call BFD_items_add(Lout, dst, -2, jm, 0, 0, X%ButterflyV%blocks(kk)%matrix)
      enddo
      do kk = 1, Y%ButterflyU%nblk_loc
         jm = (kk - 1)*Y%ButterflyU%inc + Y%ButterflyU%idx
         dst = BFD_lay_owner(dc%f(K)%lay, soM, [1, (jm + 1)/2, 1])*dc%stride
         call BFD_items_add(Lout, dst, -3, jm, 0, 0, Y%ButterflyU%blocks(kk)%matrix)
      enddo
      !>**** the leaves U^X and V^Y
      do kk = 1, X%ButterflyU%nblk_loc
         t = (kk - 1)*X%ButterflyU%inc + X%ButterflyU%idx
         dst = BFD_lay_owner(dc%Ulay, [K, 0, 0], [t, 1, 1])*dc%stride
         call BFD_items_add(Lout, dst, 0, t, 0, 0, X%ButterflyU%blocks(kk)%matrix)
      enddo
      do kk = 1, Y%ButterflyV%nblk_loc
         t = (kk - 1)*Y%ButterflyV%inc + Y%ButterflyV%idx
         dst = BFD_lay_owner(dc%Vlay, [0, 0, K], [1, 1, t])*dc%stride
         call BFD_items_add(Lout, dst, -5, t, 0, 0, Y%ButterflyV%blocks(kk)%matrix)
      enddo

      call BFD_exchange(Lout, Lin, dc%comm)

      if (dc%myv >= 0) then
         allocate (Xk1(2, 2**K), VX(2**K), UY(2**K), YkK(2**K, 2))
         do kk = 1, Lin%n
            q = Lin%it(kk)%hdr(1)
            if (q >= 1) then
               call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, so)
               call BFD_ncoord(so, Lin%it(kk)%hdr(2), qo)
               li = BFD_lay_lidx(dc%f(q)%lay, so, qo)
               call BFD_matset(dc%f(q)%blk(Lin%it(kk)%hdr(3), li), Lin%it(kk)%mat)
            elseif (q == 0) then
               li = BFD_lay_lidx(dc%Ulay, [K, 0, 0], [Lin%it(kk)%hdr(2), 1, 1])
               call BFD_matset(dc%U(li), Lin%it(kk)%mat)
            elseif (q == -5) then
               li = BFD_lay_lidx(dc%Vlay, [0, 0, K], [1, 1, Lin%it(kk)%hdr(2)])
               call BFD_matset(dc%V(li), Lin%it(kk)%mat)
            elseif (q == -1) then
               call BFD_matset(Xk1(Lin%it(kk)%hdr(2), Lin%it(kk)%hdr(3)), Lin%it(kk)%mat)
            elseif (q == -2) then
               call BFD_matset(VX(Lin%it(kk)%hdr(2)), Lin%it(kk)%mat)
            elseif (q == -3) then
               call BFD_matset(UY(Lin%it(kk)%hdr(2)), Lin%it(kk)%mat)
            elseif (q == -4) then
               call BFD_matset(YkK(Lin%it(kk)%hdr(2), Lin%it(kk)%hdr(3)), Lin%it(kk)%mat)
            endif
         enddo
         !>**** M(d, (a1,u)) = sum_e X ker1(a1, jm) * V^X(jm)^T U^Y(jm) * Y kerK(jm, d), jm = 2u-2+e
         !$omp parallel do default(shared) private(li, qo, d, e, jm, T1, T2) schedule(dynamic)
         do li = 1, BFD_lay_nloc(dc%f(K)%lay, soM)
            nullify (T1%matrix, T2%matrix)
            call BFD_lay_node(dc%f(K)%lay, soM, dc%myv, li, qo)
            do d = 1, 2
               do e = 1, 2
                  jm = 2*qo(2) - 2 + e
                  call BFD_mulnew('T', VX(jm)%matrix, 'N', UY(jm)%matrix, T1, stats)
                  call BFD_mulnew('N', T1%matrix, 'N', YkK(jm, d)%matrix, T2, stats)
                  if (e == 1) then
                     call BFD_mulnew('N', Xk1(qo(1), jm)%matrix, 'N', T2%matrix, dc%f(K)%blk(d, li), stats)
                  else
                     call BFD_gemm('N', 'N', Xk1(qo(1), jm)%matrix, T2%matrix, dc%f(K)%blk(d, li)%matrix, BPACK_cone, BPACK_cone, stats)
                  endif
               enddo
            enddo
            call BFD_matfree(T1)
            call BFD_matfree(T2)
         enddo
         !$omp end parallel do
         do jm = 1, 2**K
            call BFD_matfree(Xk1(1, jm))
            call BFD_matfree(Xk1(2, jm))
            call BFD_matfree(VX(jm))
            call BFD_matfree(UY(jm))
            call BFD_matfree(YkK(jm, 1))
            call BFD_matfree(YkK(jm, 2))
         enddo
         deallocate (Xk1, VX, UY, YkK)
      endif
      call BFD_items_free(Lin)
   end subroutine BFD_dch_build

   !>**** move the leaves U (state (K,0,0)) or V (state (0,0,K)) to a new layout
   subroutine BFD_dch_moveleaves(dc, which, lay)
      implicit none
      type(bfd_dchain)::dc
      character which
      type(bfd_layout)::lay
      type(bfd_items)::Lout, Lin
      type(butterflymatrix), allocatable::new(:)
      integer s(3), li, q(3), dst, k
      if (which == 'U') then
         s = [dc%K, 0, 0]
      else
         s = [0, 0, dc%K]
      endif
      if (dc%myv >= 0) then
         allocate (new(2**(dc%K - dc%p)))
         do li = 1, size(new)
            if (which == 'U') then
               call BFD_lay_node(dc%Ulay, s, dc%myv, li, q)
            else
               call BFD_lay_node(dc%Vlay, s, dc%myv, li, q)
            endif
            dst = BFD_lay_owner(lay, s, q)
            if (dst == dc%myv) then
               if (which == 'U') then
                  call BFD_matmove(dc%U(li), new(BFD_lay_lidx(lay, s, q)))
               else
                  call BFD_matmove(dc%V(li), new(BFD_lay_lidx(lay, s, q)))
               endif
            elseif (which == 'U') then
               call BFD_items_add(Lout, dst*dc%stride, BFD_nid(s, q), 0, 0, 0, dc%U(li)%matrix)
            else
               call BFD_items_add(Lout, dst*dc%stride, BFD_nid(s, q), 0, 0, 0, dc%V(li)%matrix)
            endif
         enddo
      endif
      call BFD_exchange(Lout, Lin, dc%comm)
      if (dc%myv >= 0) then
         do k = 1, Lin%n
            call BFD_ncoord(s, Lin%it(k)%hdr(1), q)
            call BFD_matset(new(BFD_lay_lidx(lay, s, q)), Lin%it(k)%mat)
         enddo
         if (which == 'U') then
            call BFD_matarray_free(dc%U)
            call move_alloc(new, dc%U)
         else
            call BFD_matarray_free(dc%V)
            call move_alloc(new, dc%V)
         endif
      endif
      call BFD_items_free(Lin)
      if (which == 'U') then
         dc%Ulay = lay
      else
         dc%Vlay = lay
      endif
   end subroutine BFD_dch_moveleaves

   !>**** move factor q to a new layout (U follows f(1), V follows f(nf))
   subroutine BFD_dch_redist(dc, q, lay)
      implicit none
      type(bfd_dchain)::dc
      integer q
      type(bfd_layout)::lay
      type(bfd_items)::Lout, Lin
      type(butterflymatrix), allocatable::new(:, :)
      integer so(3), li, qo(3), dst, d, k, lnew
      call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, so)
      if (dc%myv >= 0) then
         allocate (new(2, BFD_lay_nloc(lay, so)))
         do li = 1, size(dc%f(q)%blk, 2)
            call BFD_lay_node(dc%f(q)%lay, so, dc%myv, li, qo)
            dst = BFD_lay_owner(lay, so, qo)
            do d = 1, 2
               if (dst == dc%myv) then
                  call BFD_matmove(dc%f(q)%blk(d, li), new(d, BFD_lay_lidx(lay, so, qo)))
               else
                  call BFD_items_add(Lout, dst*dc%stride, BFD_nid(so, qo), d, 0, 0, dc%f(q)%blk(d, li)%matrix)
               endif
            enddo
         enddo
      endif
      call BFD_exchange(Lout, Lin, dc%comm)
      if (dc%myv >= 0) then
         do k = 1, Lin%n
            call BFD_ncoord(so, Lin%it(k)%hdr(1), qo)
            lnew = BFD_lay_lidx(lay, so, qo)
            call BFD_matset(new(Lin%it(k)%hdr(2), lnew), Lin%it(k)%mat)
         enddo
         call BFD_dfac_free(dc%f(q))
         call move_alloc(new, dc%f(q)%blk)
      endif
      call BFD_items_free(Lin)
      dc%f(q)%lay = lay
      if (q == 1) call BFD_dch_moveleaves(dc, 'U', lay)
      if (q == dc%nf) call BFD_dch_moveleaves(dc, 'V', lay)
   end subroutine BFD_dch_redist

   !>**** make f(q) and f(q+1) share one layout valid for both, so that a swap or merge at q is local
   subroutine BFD_dch_pair(dc, q)
      implicit none
      type(bfd_dchain)::dc
      integer q, c1(3), c2(3)
      type(bfd_layout)::lay
      logical v11, v12, v21, v22
      v11 = BFD_lay_valid(dc%f(q)%lay, dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec)
      v12 = BFD_lay_valid(dc%f(q)%lay, dc%f(q + 1)%s, dc%f(q + 1)%inc, dc%f(q + 1)%dec)
      v21 = BFD_lay_valid(dc%f(q + 1)%lay, dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec)
      v22 = BFD_lay_valid(dc%f(q + 1)%lay, dc%f(q + 1)%s, dc%f(q + 1)%inc, dc%f(q + 1)%dec)
      if (BFD_lay_same(dc%f(q)%lay, dc%f(q + 1)%lay) .and. v11 .and. v12) return
      if (v11 .and. v12) then
         lay = dc%f(q)%lay
         call BFD_dch_redist(dc, q + 1, lay)
      elseif (v21 .and. v22) then
         lay = dc%f(q + 1)%lay
         call BFD_dch_redist(dc, q, lay)
      else
         call BFD_dfac_caps(dc%f(q), c1)
         call BFD_dfac_caps(dc%f(q + 1), c2)
         call BFD_dch_canon(min(c1, c2), dc%p, 2, lay)
         call BFD_dch_redist(dc, q, lay)
         call BFD_dch_redist(dc, q + 1, lay)
      endif
   end subroutine BFD_dch_pair

   !>**** positions lo..hi (0 = U): column-orthonormalize (QR of every column block) and push R into the next position;
   !> the pushes cross ranks only where two neighbouring factors have different layouts
   subroutine BFD_dch_rsweep(dc, lo, hi, stats)
      implicit none
      type(bfd_dchain)::dc
      integer lo, hi
      type(Hstat)::stats
      type(bfd_items)::Lout, Lin
      DT, allocatable::W(:, :), Z(:, :), T(:, :)
      integer q, li, qi(3), qo1(3), qo2(3), d, l1, l2, h1, h2, dst, dd, lg, kk, so(3), K, nl
      integer, allocatable::dsts(:), qis(:, :)
      type(butterflymatrix), allocatable::Rs(:)
      logical cross
      K = dc%K
      do q = lo, hi
         if (q == 0) then
            call assert(BFD_lay_same(dc%Ulay, dc%f(1)%lay), 'BFD_dch_rsweep: U and f(1) layouts differ')
            if (dc%myv >= 0) then
               do li = 1, size(dc%U)
                  call BFD_qr(dc%U(li)%matrix, W, Z, stats)
                  call BFD_matset(dc%U(li), W)
                  do dd = 1, 2
                     call BFD_lmul(Z, dc%f(1)%blk(dd, li), stats)
                  enddo
               enddo
            endif
            cycle
         endif
         if (q < dc%nf) then
            cross = .not. BFD_lay_same(dc%f(q)%lay, dc%f(q + 1)%lay)
         else
            cross = .not. BFD_lay_same(dc%f(q)%lay, dc%Vlay)
         endif
         call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, so)
         if (dc%myv >= 0) then
            nl = BFD_lay_nloc(dc%f(q)%lay, dc%f(q)%s)
            allocate (Rs(nl), dsts(nl), qis(3, nl))
            !$omp parallel do default(shared) private(li, qi, qo1, qo2, d, l1, l2, h1, h2, T, W, Z) schedule(dynamic)
            do li = 1, nl
               call BFD_lay_node(dc%f(q)%lay, dc%f(q)%s, dc%myv, li, qi)
               qis(:, li) = qi
               call BFD_q_fusers(dc%f(q)%inc, dc%f(q)%dec, qi, qo1, qo2, d)
               l1 = BFD_lay_lidx(dc%f(q)%lay, so, qo1)
               l2 = BFD_lay_lidx(dc%f(q)%lay, so, qo2)
               h1 = size(dc%f(q)%blk(d, l1)%matrix, 1)
               h2 = size(dc%f(q)%blk(d, l2)%matrix, 1)
               allocate (T(h1 + h2, size(dc%f(q)%blk(d, l1)%matrix, 2)))
               T(1:h1, :) = dc%f(q)%blk(d, l1)%matrix
               T(h1 + 1:, :) = dc%f(q)%blk(d, l2)%matrix
               call BFD_qr(T, W, Z, stats)
               call BFD_matset(dc%f(q)%blk(d, l1), W(1:h1, :))
               call BFD_matset(dc%f(q)%blk(d, l2), W(h1 + 1:h1 + h2, :))
               call BFD_matset(Rs(li), Z)
               deallocate (T)
               if (q < dc%nf) then
                  dsts(li) = BFD_lay_owner(dc%f(q + 1)%lay, dc%f(q)%s, qi)
               else
                  dsts(li) = BFD_lay_owner(dc%Vlay, dc%f(q)%s, qi)
               endif
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(li) schedule(dynamic)
            do li = 1, nl
               if (dsts(li) == dc%myv) call BFD_dch_rpush(dc, q, qis(:, li), Rs(li)%matrix, stats)
            enddo
            !$omp end parallel do
            do li = 1, nl
               if (dsts(li) /= dc%myv) call BFD_items_add(Lout, dsts(li)*dc%stride, BFD_nid(dc%f(q)%s, qis(:, li)), 0, 0, 0, Rs(li)%matrix)
            enddo
            call BFD_matarray_free(Rs)
            deallocate (dsts, qis)
         endif
         if (cross) then
            call BFD_exchange(Lout, Lin, dc%comm)
            do kk = 1, Lin%n
               call BFD_ncoord(dc%f(q)%s, Lin%it(kk)%hdr(1), qi)
               call BFD_dch_rpush(dc, q, qi, Lin%it(kk)%mat, stats)
            enddo
            call BFD_items_free(Lin)
         endif
      enddo
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
   end subroutine BFD_dch_rsweep

   !>**** R of the column block of f(q) at its input node qi goes into the next position (f(q+1) or V)
   subroutine BFD_dch_rpush(dc, q, qi, R, stats)
      implicit none
      type(bfd_dchain)::dc
      integer q, qi(3), lg, dd
      DT::R(:, :)
      type(Hstat)::stats
      if (q < dc%nf) then
         lg = BFD_lay_lidx(dc%f(q + 1)%lay, dc%f(q)%s, qi)
         do dd = 1, 2
            call BFD_lmul(R, dc%f(q + 1)%blk(dd, lg), stats)
         enddo
      else
         lg = BFD_lay_lidx(dc%Vlay, dc%f(q)%s, qi)
         call BFD_rmul(dc%V(lg), transpose(R), stats)
      endif
   end subroutine BFD_dch_rpush

   !>**** positions hi down to lo (nf+1 = V): row-orthonormalize (LQ of every row block) and push L into the previous
   !> position; the pushes cross ranks only where two neighbouring factors have different layouts
   subroutine BFD_dch_lsweep(dc, lo, hi, stats)
      implicit none
      type(bfd_dchain)::dc
      integer lo, hi
      type(Hstat)::stats
      type(bfd_items)::Lout, Lin
      DT, allocatable::W(:, :), Z(:, :), T(:, :)
      integer q, li, qv(3), qo(3), qo1(3), qo2(3), d, l1, l2, n1, n2, dst, kk, so(3), gso(3), K, nl
      integer, allocatable::dsts(:), qos(:, :)
      type(butterflymatrix), allocatable::Ls(:)
      logical cross
      K = dc%K
      do q = hi, lo, -1
         if (q == dc%nf + 1) then
            call assert(BFD_lay_same(dc%Vlay, dc%f(dc%nf)%lay), 'BFD_dch_lsweep: V and f(nf) layouts differ')
            if (dc%myv >= 0) then
               call BFD_q_sout(dc%f(dc%nf)%s, dc%f(dc%nf)%inc, dc%f(dc%nf)%dec, so)
               do li = 1, size(dc%V)
                  call BFD_lay_node(dc%Vlay, [0, 0, K], dc%myv, li, qv)
                  call BFD_qr(dc%V(li)%matrix, W, Z, stats)
                  call BFD_matset(dc%V(li), W)
                  call BFD_q_fusers(dc%f(dc%nf)%inc, dc%f(dc%nf)%dec, qv, qo1, qo2, d)
                  call BFD_rmul(dc%f(dc%nf)%blk(d, BFD_lay_lidx(dc%f(dc%nf)%lay, so, qo1)), transpose(Z), stats)
                  call BFD_rmul(dc%f(dc%nf)%blk(d, BFD_lay_lidx(dc%f(dc%nf)%lay, so, qo2)), transpose(Z), stats)
               enddo
            endif
            cycle
         endif
         if (q > 1) then
            cross = .not. BFD_lay_same(dc%f(q)%lay, dc%f(q - 1)%lay)
         else
            cross = .not. BFD_lay_same(dc%f(q)%lay, dc%Ulay)
         endif
         call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, so)
         if (dc%myv >= 0) then
            nl = size(dc%f(q)%blk, 2)
            allocate (Ls(nl), dsts(nl), qos(3, nl))
            !$omp parallel do default(shared) private(li, qo, n1, n2, T, W, Z, qo1, qo2, d, gso) schedule(dynamic)
            do li = 1, nl
               call BFD_lay_node(dc%f(q)%lay, so, dc%myv, li, qo)
               qos(:, li) = qo
               n1 = size(dc%f(q)%blk(1, li)%matrix, 2)
               n2 = size(dc%f(q)%blk(2, li)%matrix, 2)
               allocate (T(size(dc%f(q)%blk(1, li)%matrix, 1), n1 + n2))
               T(:, 1:n1) = dc%f(q)%blk(1, li)%matrix
               T(:, n1 + 1:) = dc%f(q)%blk(2, li)%matrix
               call BFD_qr(transpose(T), W, Z, stats)
               call BFD_matset(dc%f(q)%blk(1, li), transpose(W(1:n1, :)))
               call BFD_matset(dc%f(q)%blk(2, li), transpose(W(n1 + 1:n1 + n2, :)))
               call BFD_matset(Ls(li), transpose(Z))
               deallocate (T)
               if (q > 1) then
                  call BFD_q_fusers(dc%f(q - 1)%inc, dc%f(q - 1)%dec, qo, qo1, qo2, d)
                  call BFD_q_sout(dc%f(q - 1)%s, dc%f(q - 1)%inc, dc%f(q - 1)%dec, gso)
                  dsts(li) = BFD_lay_owner(dc%f(q - 1)%lay, gso, qo1)
               else
                  dsts(li) = BFD_lay_owner(dc%Ulay, [K, 0, 0], qo)
               endif
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(li) schedule(dynamic)
            do li = 1, nl
               if (dsts(li) == dc%myv) call BFD_dch_lpush(dc, q, qos(:, li), Ls(li)%matrix, stats)
            enddo
            !$omp end parallel do
            do li = 1, nl
               if (dsts(li) /= dc%myv) call BFD_items_add(Lout, dsts(li)*dc%stride, BFD_nid(so, qos(:, li)), 0, 0, 0, Ls(li)%matrix)
            enddo
            call BFD_matarray_free(Ls)
            deallocate (dsts, qos)
         endif
         if (cross) then
            call BFD_exchange(Lout, Lin, dc%comm)
            do kk = 1, Lin%n
               call BFD_ncoord(so, Lin%it(kk)%hdr(1), qo)
               call BFD_dch_lpush(dc, q, qo, Lin%it(kk)%mat, stats)
            enddo
            call BFD_items_free(Lin)
         endif
      enddo
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
   end subroutine BFD_dch_lsweep

   !>**** L of the row block of f(q) at its output node qo goes into the previous position (f(q-1) or U)
   subroutine BFD_dch_lpush(dc, q, qo, L, stats)
      implicit none
      type(bfd_dchain)::dc
      integer q, qo(3), qo1(3), qo2(3), d, gso(3)
      DT::L(:, :)
      type(Hstat)::stats
      if (q > 1) then
         call BFD_q_fusers(dc%f(q - 1)%inc, dc%f(q - 1)%dec, qo, qo1, qo2, d)
         call BFD_q_sout(dc%f(q - 1)%s, dc%f(q - 1)%inc, dc%f(q - 1)%dec, gso)
         call BFD_rmul(dc%f(q - 1)%blk(d, BFD_lay_lidx(dc%f(q - 1)%lay, gso, qo1)), L, stats)
         call BFD_rmul(dc%f(q - 1)%blk(d, BFD_lay_lidx(dc%f(q - 1)%lay, gso, qo2)), L, stats)
      else
         call BFD_rmul(dc%U(BFD_lay_lidx(dc%Ulay, [dc%K, 0, 0], qo)), L, stats)
      endif
   end subroutine BFD_dch_lpush

   !>**** BRowSwap f(q) ('ca') * f(q+1) ('cb') = Bf ('cb') * Cf ('ca') followed by the truncation of the new state (as
   !> BFD_ch_swap + BFD_ch_ctrunc); Cf's stacked identities are not formed: its blocks are the column blocks of the
   !> truncation's right factor. Local once both factors share a layout.
   subroutine BFD_dch_swap(dc, q, tol, stats)
      implicit none
      type(bfd_dchain)::dc
      integer q
      real(kind=8)::tol
      type(Hstat)::stats
      type(bfd_dfac)::Bf, Cf
      type(butterflymatrix)::T
      DT, allocatable::W(:, :), Z(:, :), TT(:, :)
      DTR, allocatable::S(:)
      integer, allocatable::split(:, :)
      integer li, qo(3), mo(3), qm(3), qo1(3), qo2(3), d, lm, l1, l2, h1, h2, r, j, n1, s1(3), s2(3), nloc
      call BFD_dch_pair(dc, q)
      s1 = dc%f(q)%s
      call BFD_q_sout(s1, dc%f(q)%inc, dc%f(q)%dec, s2)
      Cf%s = dc%f(q + 1)%s
      Cf%inc = 1
      Cf%dec = 3
      Cf%lay = dc%f(q)%lay
      Bf%s = [Cf%s(1) + 1, Cf%s(2), Cf%s(3) - 1]
      Bf%inc = 2
      Bf%dec = 3
      Bf%lay = dc%f(q)%lay
      call assert(BFD_lay_valid(Bf%lay, Bf%s, Bf%inc, Bf%dec) .and. BFD_lay_valid(Cf%lay, Cf%s, Cf%inc, Cf%dec), 'BFD_dch_swap: layout not valid')
      call BFD_dfac_alloc(dc, Bf)
      call BFD_dfac_alloc(dc, Cf)
      if (dc%myv >= 0) then
         nloc = size(Bf%blk, 2)
         allocate (split(2, nloc))
         !$omp parallel do default(shared) private(li, qo, d, mo, lm, T) schedule(dynamic)
         do li = 1, nloc
            nullify (T%matrix)
            call BFD_lay_node(Bf%lay, s2, dc%myv, li, qo)
            do d = 1, 2
               call BFD_q_fin(dc%f(q)%inc, dc%f(q)%dec, qo, d, mo)
               lm = BFD_lay_lidx(dc%f(q + 1)%lay, s1, mo)
               call BFD_hcat(dc%f(q + 1)%blk(1, lm)%matrix, dc%f(q + 1)%blk(2, lm)%matrix, BPACK_cone, T)
               call BFD_mulnew('N', dc%f(q)%blk(d, li)%matrix, 'N', T%matrix, Bf%blk(d, li), stats)
               split(d, li) = size(dc%f(q + 1)%blk(1, lm)%matrix, 2)
            enddo
            call BFD_matfree(T)
         enddo
         !$omp end parallel do
         !$omp parallel do default(shared) private(li, qm, qo1, qo2, d, l1, l2, h1, h2, TT, W, S, Z, r, j, n1) schedule(dynamic)
         do li = 1, nloc
            call BFD_lay_node(Cf%lay, Bf%s, dc%myv, li, qm)
            call BFD_q_fusers(Bf%inc, Bf%dec, qm, qo1, qo2, d)
            l1 = BFD_lay_lidx(Bf%lay, s2, qo1)
            l2 = BFD_lay_lidx(Bf%lay, s2, qo2)
            h1 = size(Bf%blk(d, l1)%matrix, 1)
            h2 = size(Bf%blk(d, l2)%matrix, 1)
            allocate (TT(h1 + h2, size(Bf%blk(d, l1)%matrix, 2)))
            TT(1:h1, :) = Bf%blk(d, l1)%matrix
            TT(h1 + 1:, :) = Bf%blk(d, l2)%matrix
            call BFD_svd(TT, W, S, Z, r, tol, .true., stats)
            do j = 1, r
               W(:, j) = W(:, j)*S(j)
            enddo
            call BFD_matset(Bf%blk(d, l1), W(1:h1, 1:r))
            call BFD_matset(Bf%blk(d, l2), W(h1 + 1:h1 + h2, 1:r))
            n1 = split(d, l1)
            call BFD_matset(Cf%blk(1, li), Z(1:r, 1:n1))
            call BFD_matset(Cf%blk(2, li), Z(1:r, n1 + 1:))
            deallocate (TT)
         enddo
         !$omp end parallel do
         deallocate (split)
      endif
      call BFD_dfac_move(Bf, dc%f(q))
      call BFD_dfac_move(Cf, dc%f(q + 1))
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
   end subroutine BFD_dch_swap

   !>**** f(q) ('ba') * f(q+1) ('cb') = one 'ca' factor (as BFD_ch_merge); local once both share a layout
   subroutine BFD_dch_merge(dc, q, stats)
      implicit none
      type(bfd_dchain)::dc
      integer q
      type(Hstat)::stats
      type(bfd_dfac)::Pf
      integer li, qo(3), m1(3), m2(3), l1, l2, e, s2(3), k
      call BFD_dch_pair(dc, q)
      call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, s2)
      Pf%s = dc%f(q + 1)%s
      Pf%inc = 1
      Pf%dec = 3
      Pf%lay = dc%f(q)%lay
      call assert(BFD_lay_valid(Pf%lay, Pf%s, Pf%inc, Pf%dec), 'BFD_dch_merge: layout not valid')
      call BFD_dfac_alloc(dc, Pf)
      if (dc%myv >= 0) then
         !$omp parallel do default(shared) private(li, qo, m1, m2, l1, l2, e) schedule(dynamic)
         do li = 1, size(Pf%blk, 2)
            call BFD_lay_node(Pf%lay, s2, dc%myv, li, qo)
            call BFD_q_fin(dc%f(q)%inc, dc%f(q)%dec, qo, 1, m1)
            call BFD_q_fin(dc%f(q)%inc, dc%f(q)%dec, qo, 2, m2)
            l1 = BFD_lay_lidx(dc%f(q + 1)%lay, dc%f(q)%s, m1)
            l2 = BFD_lay_lidx(dc%f(q + 1)%lay, dc%f(q)%s, m2)
            do e = 1, 2
               call BFD_mulnew('N', dc%f(q)%blk(1, li)%matrix, 'N', dc%f(q + 1)%blk(e, l1)%matrix, Pf%blk(e, li), stats)
               call BFD_gemm('N', 'N', dc%f(q)%blk(2, li)%matrix, dc%f(q + 1)%blk(e, l2)%matrix, Pf%blk(e, li)%matrix, BPACK_cone, BPACK_cone, stats)
            enddo
         enddo
         !$omp end parallel do
      endif
      call BFD_dfac_move(Pf, dc%f(q))
      call BFD_dfac_free(dc%f(q + 1))
      do k = q + 1, dc%nf - 1
         call BFD_dfac_move(dc%f(k + 1), dc%f(k))
      enddo
      dc%nf = dc%nf - 1
   end subroutine BFD_dch_merge

   !>**** the final chain U * K 'ca' factors * V sent to the ButterflyPACK layout of Z (level_half = K/2), whose structure
   !> is set up from X's rows and Y's columns
   subroutine BFD_dch_extract(dc, X, Y, Z, option, ptree, msh)
      implicit none
      type(bfd_dchain)::dc
      type(matrixblock)::X, Y, Z
      type(Hoption)::option
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::tmpl
      type(bfd_items)::Lout, Lin
      integer K, q, l, li, qo(3), d, i, j, pgno_sub, dst, head, kk, ii, jj, so(3)
      K = dc%K
      call assert(dc%nf == K, 'BFD_dch_extract: chain is not a butterfly')
      tmpl%M = X%M
      tmpl%N = Y%N
      tmpl%headm = X%headm
      tmpl%headn = Y%headn
      tmpl%level = X%level
      tmpl%M_loc = X%M_loc
      tmpl%N_loc = Y%N_loc
      tmpl%pgno = X%pgno
      tmpl%pgno_db = X%pgno_db
      tmpl%M_p => X%M_p
      tmpl%N_p => Y%N_p
      if (associated(X%ms)) tmpl%ms => X%ms
      if (associated(Y%ns)) tmpl%ns => Y%ns
      call BF_Init_randomized(K, 0, X%row_group, Y%col_group, tmpl, Z, msh, ptree, option, 0)
      head = ptree%pgrp(dc%pgno)%head
      if (dc%myv >= 0) then
         do q = 1, K
            l = K + 1 - q
            call BFD_q_sout(dc%f(q)%s, dc%f(q)%inc, dc%f(q)%dec, so)
            call assert(so(1) == l .and. so(2) == 0 .and. dc%f(q)%inc == 1 .and. dc%f(q)%dec == 3, 'BFD_dch_extract: unexpected factor')
            do li = 1, size(dc%f(q)%blk, 2)
               call BFD_lay_node(dc%f(q)%lay, so, dc%myv, li, qo)
               do d = 1, 2
                  i = qo(1)
                  j = 2*qo(3) - 2 + d
                  if (l <= Z%level_half) then
                     call GetBlockPID(ptree, dc%pgno, l, K, i, qo(3), 'R', pgno_sub)
                  else
                     call GetBlockPID(ptree, dc%pgno, l, K, (i + 1)/2, j, 'C', pgno_sub)
                  endif
                  dst = ptree%pgrp(pgno_sub)%head - head
                  call BFD_items_add(Lout, dst, l, i, j, 0, dc%f(q)%blk(d, li)%matrix)
               enddo
            enddo
         enddo
         do li = 1, size(dc%U)
            call BFD_lay_node(dc%Ulay, [K, 0, 0], dc%myv, li, qo)
            call GetBlockPID(ptree, dc%pgno, K + 1, K, qo(1), 1, 'C', pgno_sub)
            call BFD_items_add(Lout, ptree%pgrp(pgno_sub)%head - head, K + 1, qo(1), 0, 0, dc%U(li)%matrix)
         enddo
         do li = 1, size(dc%V)
            call BFD_lay_node(dc%Vlay, [0, 0, K], dc%myv, li, qo)
            call GetBlockPID(ptree, dc%pgno, 0, K, 1, qo(3), 'R', pgno_sub)
            call BFD_items_add(Lout, ptree%pgrp(pgno_sub)%head - head, 0, qo(3), 0, 0, dc%V(li)%matrix)
         enddo
      endif
      call BFD_exchange(Lout, Lin, dc%comm)
      do kk = 1, Lin%n
         l = Lin%it(kk)%hdr(1)
         i = Lin%it(kk)%hdr(2)
         j = Lin%it(kk)%hdr(3)
         if (l == K + 1) then
            call BFD_matset(Z%ButterflyU%blocks((i - Z%ButterflyU%idx)/Z%ButterflyU%inc + 1), Lin%it(kk)%mat)
         elseif (l == 0) then
            call BFD_matset(Z%ButterflyV%blocks((i - Z%ButterflyV%idx)/Z%ButterflyV%inc + 1), Lin%it(kk)%mat)
         else
            ii = (i - Z%ButterflyKerl(l)%idx_r)/Z%ButterflyKerl(l)%inc_r + 1
            jj = (j - Z%ButterflyKerl(l)%idx_c)/Z%ButterflyKerl(l)%inc_c + 1
            call BFD_matset(Z%ButterflyKerl(l)%blocks(ii, jj), Lin%it(kk)%mat)
         endif
      enddo
      call BFD_items_free(Lin)
      call BF_get_rank(Z, ptree)
   end subroutine BFD_dch_extract

   !>**** Z = X*Y for two K-level butterflies on the same process group sharing the middle leaves: the distributed
   !> version of BFD_bmult (the chain lives on 2^min(p, K-2) ranks of the 2^p ranks of the group). Truncations use tol,
   !> the result is recompressed with tol unless norecomp (the caller sums and recompresses it right away).
   !> K = 0: Z = U^X (V^X^T U^Y) V^Y^T with an allreduced core.
   subroutine BFD_bmult_dist(X, Y, Z, tol, option, stats, ptree, msh, norecomp)
      implicit none
      type(matrixblock)::X, Y, Z
      logical, optional::norecomp
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(bfd_dchain)::dc
      type(matrixblock)::tmpl
      DT, allocatable::core(:, :)
      integer K, m, it, lo, ierr, rx, ry
      real(kind=8)::t0
      logical skip
      skip = .false.
      if (present(norecomp)) skip = norecomp
      K = X%level_butterfly
      call assert(Y%level_butterfly == K .and. Y%pgno == X%pgno, 'BFD_bmult_dist: X and Y must have the same levels and process group')
      if (K == 0) then
         tmpl%M = X%M
         tmpl%N = Y%N
         tmpl%headm = X%headm
         tmpl%headn = Y%headn
         tmpl%level = X%level
         tmpl%M_loc = X%M_loc
         tmpl%N_loc = Y%N_loc
         tmpl%pgno = X%pgno
         tmpl%pgno_db = X%pgno_db
         tmpl%M_p => X%M_p
         tmpl%N_p => Y%N_p
         call BF_Init_randomized(0, 0, X%row_group, Y%col_group, tmpl, Z, msh, ptree, option, 0)
         rx = size(X%ButterflyV%blocks(1)%matrix, 2)
         ry = size(Y%ButterflyU%blocks(1)%matrix, 2)
         allocate (core(rx, ry))
         core = 0
         if (X%N_loc > 0) call BFD_gemm('T', 'N', X%ButterflyV%blocks(1)%matrix, Y%ButterflyU%blocks(1)%matrix, core, BPACK_cone, BPACK_czero, stats)
         if (ptree%pgrp(X%pgno)%nproc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, core, rx*ry, MPI_DT, MPI_SUM, ptree%pgrp(X%pgno)%Comm, ierr)
         call BFD_mulnew('N', X%ButterflyU%blocks(1)%matrix, 'N', core, Z%ButterflyU%blocks(1), stats)
         call BFD_matset(Z%ButterflyV%blocks(1), Y%ButterflyV%blocks(1)%matrix)
         deallocate (core)
         call BF_get_rank(Z, ptree)
         if (.not. skip) call BFD_recompress(Z, tol, option, stats, ptree)
         return
      endif
      t0 = MPI_Wtime()
      call BFD_dch_build(X, Y, dc, ptree, stats)
      call BFD_tadd(6, t0)
      t0 = MPI_Wtime()
      call BFD_dch_lsweep(dc, K + 1, 2*K, stats)
      call BFD_tadd(7, t0)
      do m = 1, K - 1
         lo = K - m + 1
         if (m == 1) lo = 0
         t0 = MPI_Wtime()
         call BFD_dch_rsweep(dc, lo, K - 1, stats)
         call BFD_tadd(7, t0)
         do it = 1, m
            t0 = MPI_Wtime()
            call BFD_dch_swap(dc, K + 1 - it, tol, stats)
            call BFD_tadd(8, t0)
         enddo
         t0 = MPI_Wtime()
         call BFD_dch_merge(dc, K - m, stats)
         call BFD_tadd(9, t0)
      enddo
      t0 = MPI_Wtime()
      call BFD_dch_extract(dc, X, Y, Z, option, ptree, msh)
      call BFD_dch_free(dc)
      call BFD_tadd(10, t0)
      if (.not. skip) call BFD_recompress(Z, tol, option, stats, ptree)
   end subroutine BFD_bmult_dist

!======================================================================================== products and solves with H-blocks

   !>**** P = op(A)*B, op(A) = A (trans='N') or A^T (trans='T'), A an H-block (dense, butterfly or 2x2 partitioned),
   !> B a butterfly whose rows are the columns of op(A) (paper IV-D for a partitioned A)
   recursive subroutine BFD_HxBF(A, trans, B, P, pgno, tol, option, stats, ptree, scale)
      implicit none
      type(matrixblock), target::A
      type(matrixblock)::B, P
      character trans
      integer pgno
      real(kind=8), optional::scale
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(matrixblock)::X, H1, H2, P11, P22, P12, P21, D, E
      type(matrixblock), pointer::S11, S12, S21, S22
      type(butterflymatrix), allocatable::W1(:), W2(:)
      integer K, r, mA, jf
      real(kind=8)::tol, rf0, rt0
      DT, allocatable::Dp(:, :), Db(:, :), Da(:, :)
      integer rd

      call BFD_rec_enter(1, B, stats, rd, rf0, rt0)
      K = B%level_butterfly
      if (K == 0) then
         r = size(B%ButterflyU%blocks(1)%matrix, 2)
         if (trans == 'N') then
            mA = A%M
         else
            mA = A%N
         endif
         call BFD_alloc(P, 0, pgno)
         call BFD_matnew(P%ButterflyU%blocks(1), mA, r)
         call Hmat_block_MVP_dat(A, trans, A%headm, A%headn, r, B%ButterflyU%blocks(1)%matrix, size(B%ButterflyU%blocks(1)%matrix, 1), &
                                 P%ButterflyU%blocks(1)%matrix, mA, BPACK_cone, ptree, stats)
         call BFD_matset(P%ButterflyV%blocks(1), B%ButterflyV%blocks(1)%matrix)
         call BFD_setdims(P)
      elseif (A%style == 2) then
         if (trans == 'N') then
            call BFD_bmult(A, B, P, pgno, tol, option, stats, ptree, scale)
         else
            call BFD_transpose(A, X, pgno)
            call BFD_bmult(X, B, P, pgno, tol, option, stats, ptree, scale)
            call BF_delete(X, 1)
         endif
      elseif (A%style == 4) then
         S11 => A%sons(1, 1)
         S22 => A%sons(2, 2)
         if (trans == 'N') then
            S12 => A%sons(1, 2)
            S21 => A%sons(2, 1)
         else
            S12 => A%sons(2, 1)
            S21 => A%sons(1, 2)
         endif
         call BFD_rowhalf(B, 1, H1, pgno)
         call BFD_rowhalf(B, 2, H2, pgno)
         call BFD_HxBF(S11, trans, H1, P11, pgno, tol, option, stats, ptree, scale)
         call BFD_HxBF(S22, trans, H2, P22, pgno, tol, option, stats, ptree, scale)
         call BFD_HxBF(S12, trans, H2, P12, pgno, tol, option, stats, ptree, scale)
         call BFD_HxBF(S21, trans, H1, P21, pgno, tol, option, stats, ptree, scale)
         allocate (W1(2**K), W2(2**K))
         do jf = 1, 2**K
            call BFD_matset(W1(jf), B%ButterflyKerl(1)%blocks(1, jf)%matrix)
            call BFD_matset(W2(jf), B%ButterflyKerl(1)%blocks(2, jf)%matrix)
         enddo
         if (BFD_debug(option)) then
            call BFD_rowmerge(H1, H2, W1, W2, B, D, pgno, stats)
            call BFD_todense(D, Dp, stats)
            call BFD_todense(B, Db, stats)
            call BFD_report('rowsplit/rowmerge', Dp, Db, K)
            call BF_delete(D, 1)
         endif
         call BFD_rowmerge(P11, P22, W1, W2, B, D, pgno, stats) ! diagonal part
         call BFD_rowmerge(P12, P21, W2, W1, B, E, pgno, stats) ! antidiagonal part, halves of W swapped
         call BFD_matarray_free(W1)
         call BFD_matarray_free(W2)
         call BF_delete(H1, 1)
         call BF_delete(H2, 1)
         call BF_delete(P11, 1)
         call BF_delete(P22, 1)
         call BF_delete(P12, 1)
         call BF_delete(P21, 1)
         call BFD_sum(D, E, BPACK_cone, P, pgno)
         call BF_delete(D, 1)
         call BF_delete(E, 1)
         if (BFD_debug(option)) then
            call BFD_todense(P, Dp, stats)
            call BFD_todense(B, Db, stats)
            call BFD_hbdense(A, trans, Da, ptree, stats)
            call BFD_report('HxBF(before recompress)', Dp, matmul(Da, Db), K)
            if (fnorm(Dp - matmul(Da, Db), size(Dp, 1), size(Dp, 2)) > 1d-3*fnorm(matmul(Da, Db), size(Dp, 1), size(Dp, 2))) then
               write (*, *) 'BFD_DEBUG HxBF bad: trans ', trans, ' A style', A%style, ' sons', S11%style, S12%style, S21%style, S22%style, &
                  ' A%M,N', A%M, A%N, ' B%M,N', B%M, B%N, ' S11 M,N', S11%M, S11%N, ' B rows half', B%ms(2**(K - 1))
            endif
         endif
         call BFD_recompress(P, tol, option, stats, ptree, scale)
         if (BFD_debug(option)) then
            call BFD_todense(P, Dp, stats)
            call BFD_report('HxBF', Dp, matmul(Da, Db), K)
         endif
      else
         write (*, *) 'BFD_HxBF: a dense block times a butterfly with level_butterfly>0 is not supported', A%row_group, A%col_group, K
         stop
      endif
      call BFD_rec_exit(1, stats, rd, rf0, rt0)
   end subroutine BFD_HxBF

   !>**** P = B*A, B a butterfly, A an H-block: (A^T B^T)^T
   subroutine BFD_BFxH(B, A, P, pgno, tol, option, stats, ptree, scale)
      implicit none
      type(matrixblock)::B, A, P
      integer pgno
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(matrixblock)::Bt, Pt
      call BFD_transpose(B, Bt, pgno)
      call BFD_HxBF(A, 'T', Bt, Pt, pgno, tol, option, stats, ptree, scale)
      call BFD_transpose(Pt, P, pgno)
      call BF_delete(Bt, 1)
      call BF_delete(Pt, 1)
   end subroutine BFD_BFxH

   !>**** Y = L^{-1} X (uplo='L', the unit lower triangle of the LU-factored block) or Y = U^{-T} X (uplo='U', the upper
   !> triangle), X a butterfly (paper IV-G: row split, two recursive calls)
   recursive subroutine BFD_Lsolve(Lb, uplo, X, Y, pgno, tol, option, stats, ptree)
      implicit none
      type(matrixblock), target::Lb
      type(matrixblock)::X, Y
      character uplo, tr
      integer pgno
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(matrixblock)::H1, H2, Q, T, R, Y2
      type(matrixblock), pointer::S11, S21, S22
      type(butterflymatrix), allocatable::Wt(:), Wb(:)
      integer K, r0, jf
      real(kind=8)::tol, rf0, rt0
      DT, allocatable::Dx(:, :), Dy(:, :)
      integer rd

      call BFD_rec_enter(2, X, stats, rd, rf0, rt0)
      K = X%level_butterfly
      if (K == 0) then
         call BFD_copy(X, Y)
         Y%pgno = pgno
         r0 = size(Y%ButterflyU%blocks(1)%matrix, 2)
         if (uplo == 'L') then
            call Hmat_Lsolve(Lb, 'N', Lb%headm, r0, Y%ButterflyU%blocks(1)%matrix, Y%M, ptree, stats)
         else
            call Hmat_Usolve(Lb, 'T', Lb%headm, r0, Y%ButterflyU%blocks(1)%matrix, Y%M, ptree, stats)
         endif
      elseif (Lb%style == 4) then
         S11 => Lb%sons(1, 1)
         S22 => Lb%sons(2, 2)
         if (uplo == 'L') then
            S21 => Lb%sons(2, 1)
            tr = 'N'
         else
            S21 => Lb%sons(1, 2)
            tr = 'T'
         endif
         call BFD_rowhalf(X, 1, H1, pgno)
         call BFD_rowhalf(X, 2, H2, pgno)
         call BFD_Lsolve(S11, uplo, H1, Q, pgno, tol, option, stats, ptree)
         call BFD_HxBF(S21, tr, Q, T, pgno, tol, option, stats, ptree)
         call BFD_hconcat(H2, T, -BPACK_cone, R, pgno)
         !>**** R = [H2, -T] stacks two half-size butterflies (twice the ranks, stacked column leaves): recompress it to the
         !> size of one half before the second call (paper IV-G: BF_R x V_R), otherwise the ranks of the right-hand sides
         !> grow along the recursion and the solve is no longer O(n log^3 n)
         call BFD_recompress(R, tol, option, stats, ptree)
         call BF_delete(H1, 1)
         call BF_delete(H2, 1)
         call BF_delete(T, 1)
         call BFD_Lsolve(S22, uplo, R, Y2, pgno, tol, option, stats, ptree)
         call BF_delete(R, 1)
         allocate (Wt(2**K), Wb(2**K))
         do jf = 1, 2**K
            call BFD_matset(Wt(jf), X%ButterflyKerl(1)%blocks(1, jf)%matrix)
            call BFD_vcat(X%ButterflyKerl(1)%blocks(2, jf)%matrix, X%ButterflyKerl(1)%blocks(1, jf)%matrix, Wb(jf))
         enddo
         call BFD_rowmerge(Q, Y2, Wt, Wb, X, Y, pgno, stats)
         call BFD_matarray_free(Wt)
         call BFD_matarray_free(Wb)
         call BF_delete(Q, 1)
         call BF_delete(Y2, 1)
         if (BFD_debug(option)) then
            call BFD_todense(X, Dx, stats)
            if (uplo == 'L') then
               call Hmat_Lsolve(Lb, 'N', Lb%headm, size(Dx, 2), Dx, size(Dx, 1), ptree, stats)
            else
               call Hmat_Usolve(Lb, 'T', Lb%headm, size(Dx, 2), Dx, size(Dx, 1), ptree, stats)
            endif
            call BFD_todense(Y, Dy, stats)
            call BFD_report('Lsolve(before recompress) '//uplo, Dy, Dx, K)
         endif
         call BFD_recompress(Y, tol, option, stats, ptree)
         if (BFD_debug(option)) then
            call BFD_todense(Y, Dy, stats)
            call BFD_report('Lsolve '//uplo, Dy, Dx, K)
         endif
      else
         write (*, *) 'BFD_Lsolve: a dense diagonal block with a butterfly of level_butterfly>0 is not supported', X%level_butterfly
         stop
      endif
      call BFD_rec_exit(2, stats, rd, rf0, rt0)
   end subroutine BFD_Lsolve

   !>**** estimate of the spectral norm of a single-process H-block (dense, butterfly or partitioned): largest |A*x| over
   !> a few steps of the power iteration on A^H*A with a small block of random unit vectors (a lower bound)
   real(kind=8) function BFD_normest(A, ptree, stats)
      implicit none
      type(matrixblock)::A
      type(proctree)::ptree
      type(Hstat)::stats
      integer, parameter :: nv = 2, nit = 2
      DT, allocatable::X(:, :), Y(:, :)
      integer m, n, k, j, it
      real(kind=8)::nx
      m = A%M
      n = A%N
      BFD_normest = 0d0
      if (m == 0 .or. n == 0) return
      k = min(nv, n)
      allocate (X(n, k), Y(m, k))
      call RandomMat(n, k, k, X, 0)
      do it = 1, nit
         do j = 1, k
            nx = sqrt(sum(abs(X(:, j))**2))
            if (nx > 0d0) X(:, j) = X(:, j)/nx
         enddo
         Y = 0
         call Hmat_block_MVP_dat(A, 'N', A%headm, A%headn, k, X, n, Y, m, BPACK_cone, ptree, stats)
         do j = 1, k
            BFD_normest = max(BFD_normest, sqrt(sum(abs(Y(:, j))**2)))
         enddo
         if (it == nit) exit
#if DAT==0 || DAT==2
         Y = conjg(Y)
#endif
         X = 0
         call Hmat_block_MVP_dat(A, 'T', A%headm, A%headn, k, Y, m, X, n, BPACK_cone, ptree, stats)
#if DAT==0 || DAT==2
         X = conjg(X)
#endif
      enddo
      deallocate (X, Y)
   end function BFD_normest

   !>**** the norm estimate of the H-block C stored by BFD_Hmat_normest, estimated and stored now if it is missing
   real(kind=8) function BFD_blknorm(C, ptree, stats)
      implicit none
      type(matrixblock)::C
      type(proctree)::ptree
      type(Hstat)::stats
      real(kind=8)::f0, t0
      if (C%normest < 0d0) then
         call BFD_mstart(stats, f0, t0)
         C%normest = BFD_normest(C, ptree, stats)
         call BFD_mstat(22, stats, f0, t0)
      endif
      BFD_blknorm = C%normest
   end function BFD_blknorm

   !>**** count (blks not allocated) or collect the nodes of the block tree of blk
   recursive subroutine BFD_blkcollect(blk, blks, n)
      implicit none
      type(matrixblock), target::blk
      type(block_ptr), allocatable::blks(:)
      integer n, i, j
      n = n + 1
      if (allocated(blks)) blks(n)%ptr => blk
      if (blk%style == 4) then
         do j = 1, 2
            do i = 1, 2
               call BFD_blkcollect(blk%sons(i, j), blks, n)
            enddo
         enddo
      endif
   end subroutine BFD_blkcollect

   !>**** norm estimates (BFD_normest) of all H-blocks of h_mat, i.e. of every node of the block trees of the local blocks,
   !> stored in blocks%normest before the H-LU (bf_algebra=1). They set the scale of the truncation of the products added
   !> to the blocks during the factorization (BFD_Multiply, BFD_SubHH); these updates are small relative to the blocks, so
   !> the norms are not updated. The blocks of one level are disjoint and estimated in parallel.
   subroutine BFD_Hmat_normest(h_mat, stats, ptree)
      implicit none
      type(Hmat)::h_mat
      type(Hstat)::stats
      type(proctree)::ptree
      type(block_ptr), allocatable::blks(:)
      integer i, j, k, l, n, lmin, lmax
      real(kind=8)::t0, f0
      t0 = MPI_Wtime()
      f0 = stats%Flop_Tmp
      do k = 1, 2
         n = 0
         do i = 1, h_mat%myArows
            do j = 1, h_mat%myAcols
               call BFD_blkcollect(h_mat%Local_blocks(j, i), blks, n)
            enddo
         enddo
         if (k == 1) allocate (blks(n))
      enddo
      if (n > 0) then
         lmin = minval([(blks(k)%ptr%level, k=1, n)])
         lmax = maxval([(blks(k)%ptr%level, k=1, n)])
         do l = lmin, lmax
            !$omp parallel do default(shared) private(k) schedule(dynamic)
            do k = 1, n
               if (blks(k)%ptr%level == l) blks(k)%ptr%normest = BFD_normest(blks(k)%ptr, ptree, stats)
            enddo
            !$omp end parallel do
         enddo
      endif
      deallocate (blks)
      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp - f0
      stats%Flop_Tmp = f0
      BFD_norm_time = MPI_Wtime() - t0
   end subroutine BFD_Hmat_normest

   !>**** P = A*B for two H-blocks, at least one of them a butterfly, or two dense leaves (then P is low rank). With scale
   !> (the norm of the block P is added to), all truncations are relative to max(largest singular value, scale)
   subroutine BFD_Product(A, B, P, pgno, tol, option, stats, ptree, scale)
      implicit none
      type(matrixblock)::A, B, P
      integer pgno
      real(kind=8)::tol
      real(kind=8), optional::scale
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      DT, allocatable::Id(:, :), Bd(:, :), Pd(:, :), W(:, :), Z(:, :)
      DTR, allocatable::S(:)
      integer r, j
      if (A%style == 2 .and. B%style == 2) then
         call BFD_bmult(A, B, P, pgno, tol, option, stats, ptree, scale)
      elseif (B%style == 2) then
         call BFD_HxBF(A, 'N', B, P, pgno, tol, option, stats, ptree, scale)
      elseif (A%style == 2) then
         call BFD_BFxH(A, B, P, pgno, tol, option, stats, ptree, scale)
      elseif (A%style == 1 .and. B%style == 1) then
         allocate (Id(B%N, B%N), Bd(B%M, B%N), Pd(A%M, B%N))
         Id = 0
         Bd = 0
         Pd = 0
         do j = 1, B%N
            Id(j, j) = BPACK_cone
         enddo
         call Hmat_block_MVP_dat(B, 'N', B%headm, B%headn, B%N, Id, B%N, Bd, B%M, BPACK_cone, ptree, stats)
         call Hmat_block_MVP_dat(A, 'N', A%headm, A%headn, B%N, Bd, B%M, Pd, A%M, BPACK_cone, ptree, stats)
         deallocate (Id)
         call BFD_svd(Pd, W, S, Z, r, tol, .true., stats, scale)
         do j = 1, r
            W(:, j) = W(:, j)*S(j)
         enddo
         call BFD_alloc(P, 0, pgno)
         call BFD_matset(P%ButterflyU%blocks(1), W(:, 1:r))
         call BFD_matset(P%ButterflyV%blocks(1), transpose(Z(1:r, :)))
         call BFD_setdims(P)
         deallocate (Bd, Pd, W, Z, S)
      else
         write (*, *) 'BFD_Product: unsupported block styles', A%style, B%style
         stop
      endif
   end subroutine BFD_Product

!======================================================================================== split / merge at native levels

   !>**** P%sons(2,2) = the children of the butterfly P at their native level max(K-1,0): BF_split (exact restriction to
   !> max(K-2,0) levels) followed by the exact level-up
   subroutine BFD_SplitNative(P, option, stats, ptree, msh)
      implicit none
      type(matrixblock), target::P
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::T
      integer K, i, j, Kc
      integer, allocatable::rsz(:), csz(:)
      K = P%level_butterfly
      allocate (P%sons(2, 2))
      call BF_split(P, P, ptree, stats, msh, option)
      if (K >= 2) then
         do i = 1, 2
            do j = 1, 2
               Kc = P%sons(i, j)%level_butterfly
               call assert(Kc == K - 2, 'BFD_SplitNative: unexpected level from BF_split')
               allocate (rsz(2**(Kc + 1)), csz(2**(Kc + 1)))
               call BFD_leafsizes(P%sons(i, j)%row_group, Kc + 1, msh, rsz)
               call BFD_leafsizes(P%sons(i, j)%col_group, Kc + 1, msh, csz)
               call BFD_levelup(P%sons(i, j), T, rsz, csz, P%sons(i, j)%pgno)
               call BFD_Install(P%sons(i, j), T, ptree)
               deallocate (rsz, csz)
               if (BFD_debug(option)) call BFD_checkleaves(P%sons(i, j), msh, 'SplitNative child')
            enddo
         enddo
      endif
   end subroutine BFD_SplitNative

   !>**** four low-rank children -> 1-level butterfly (BF_Aggregate needs level_butterfly >= 2)
   subroutine BFD_aggregate_lr1(part, O, pgno, tol, option, stats)
      implicit none
      type(matrixblock)::part, O
      integer pgno
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(butterflymatrix)::T, T1, T3
      DT, allocatable::W(:, :), Z(:, :)
      DTR, allocatable::S(:)
      integer a, b, r
      call BFD_alloc(O, 1, pgno)
      do a = 1, 2
         call BFD_hcat(part%sons(a, 1)%ButterflyU%blocks(1)%matrix, part%sons(a, 2)%ButterflyU%blocks(1)%matrix, BPACK_cone, T)
         call BFD_svd(T%matrix, W, S, Z, r, tol, .true., stats)
         call BFD_matset(O%ButterflyU%blocks(a), W(:, 1:r))
      enddo
      do b = 1, 2
         call BFD_hcat(part%sons(1, b)%ButterflyV%blocks(1)%matrix, part%sons(2, b)%ButterflyV%blocks(1)%matrix, BPACK_cone, T)
         call BFD_svd(T%matrix, W, S, Z, r, tol, .true., stats)
         call BFD_matset(O%ButterflyV%blocks(b), W(:, 1:r))
      enddo
      do a = 1, 2
         do b = 1, 2
            call BFD_mulnew('C', O%ButterflyU%blocks(a)%matrix, 'N', part%sons(a, b)%ButterflyU%blocks(1)%matrix, T1, stats)
            call BFD_mulnew('C', O%ButterflyV%blocks(b)%matrix, 'N', part%sons(a, b)%ButterflyV%blocks(1)%matrix, T3, stats)
            call BFD_mulnew('N', T1%matrix, 'T', T3%matrix, O%ButterflyKerl(1)%blocks(a, b), stats)
         enddo
      enddo
      call BFD_matfree(T)
      call BFD_matfree(T1)
      call BFD_matfree(T3)
      if (allocated(W)) deallocate (W)
      if (allocated(Z)) deallocate (Z)
      if (allocated(S)) deallocate (S)
      call BFD_setdims(O)
   end subroutine BFD_aggregate_lr1

   !>**** four low-rank children -> low rank
   subroutine BFD_aggregate_lr0(part, O, pgno, tol, option, stats)
      implicit none
      type(matrixblock)::part, O
      integer pgno
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      integer a, b, M, N, m1, n1, rs, r, roff, coff
      m1 = part%sons(1, 1)%M
      n1 = part%sons(1, 1)%N
      M = m1 + part%sons(2, 1)%M
      N = n1 + part%sons(1, 2)%N
      rs = 0
      do a = 1, 2
         do b = 1, 2
            rs = rs + size(part%sons(a, b)%ButterflyU%blocks(1)%matrix, 2)
         enddo
      enddo
      call BFD_alloc(O, 0, pgno)
      call BFD_matnew(O%ButterflyU%blocks(1), M, rs)
      call BFD_matnew(O%ButterflyV%blocks(1), N, rs)
      rs = 0
      do a = 1, 2
         do b = 1, 2
            r = size(part%sons(a, b)%ButterflyU%blocks(1)%matrix, 2)
            roff = 0
            if (a == 2) roff = m1
            coff = 0
            if (b == 2) coff = n1
            O%ButterflyU%blocks(1)%matrix(roff + 1:roff + part%sons(a, b)%M, rs + 1:rs + r) = part%sons(a, b)%ButterflyU%blocks(1)%matrix
            O%ButterflyV%blocks(1)%matrix(coff + 1:coff + part%sons(a, b)%N, rs + 1:rs + r) = part%sons(a, b)%ButterflyV%blocks(1)%matrix
            rs = rs + r
         enddo
      enddo
      call BFD_setdims(O)
      call BFD_recompress_lr(O, tol, stats)
      call BFD_setdims(O)
   end subroutine BFD_aggregate_lr0

   !>**** the butterfly O (with the cluster metadata of the target tmpl, level K) from the children part%sons at their native
   !> level max(K-1,0): level-down + BF_Aggregate for K>=2
   subroutine BFD_MergeNative(part, tmpl, O, tol, option, stats, ptree, msh)
      implicit none
      type(matrixblock), target::part
      type(matrixblock)::tmpl, O
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::T
      integer K, i, j
      K = tmpl%level_butterfly
      if (K >= 2) then
         do i = 1, 2
            do j = 1, 2
               call BFD_leveldown(part%sons(i, j), T, part%sons(i, j)%pgno, stats)
               call BFD_Install(part%sons(i, j), T, ptree)
               call BFD_recompress(part%sons(i, j), tol, option, stats, ptree)
            enddo
         enddo
         call BF_Init_randomized(K, tmpl%rankmax, tmpl%row_group, tmpl%col_group, tmpl, O, msh, ptree, option, 1)
         call BF_Aggregate(part, O, ptree, stats, option, msh)
         call BFD_recompress(O, tol, option, stats, ptree)
      elseif (K == 1) then
         call BFD_aggregate_lr1(part, O, tmpl%pgno, tol, option, stats)
         call BFD_recompress(O, tol, option, stats, ptree)
      else
         call BFD_aggregate_lr0(part, O, tmpl%pgno, tol, option, stats)
      endif
      call BFD_copymeta(O, tmpl)
      call BF_get_rank(O, ptree)
   end subroutine BFD_MergeNative

!======================================================================================== entry points for the H-BF LU

   !>**** res = L^{-1} X (Hmat_LXM): Lb is the LU-factored diagonal block, X a butterfly in its block row
   subroutine BFD_LXM(Lb, X, res, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::Lb, X, res
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      call BFD_check_seq(X, ptree)
      if (BFD_debug(option)) call BFD_checkleaves(X, msh, 'LXM input')
      call BFD_Lsolve(Lb, 'L', X, res, X%pgno, option%tol_rand*BFD_tolfac(), option, stats, ptree)
      if (BFD_tolfac() < 1d0) call BFD_recompress(res, option%tol_rand, option, stats, ptree)
      call BFD_copymeta(res, X)
      call BF_get_rank(res, ptree)
      if (BFD_debug(option)) call BFD_checkleaves(res, msh, 'LXM output')
   end subroutine BFD_LXM

   !>**** res = X U^{-1} (Hmat_XUM): Ub is the LU-factored diagonal block, X a butterfly in its block column
   subroutine BFD_XUM(Ub, X, res, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::Ub, X, res
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::Xt, Yt
      call BFD_check_seq(X, ptree)
      if (BFD_debug(option)) call BFD_checkleaves(X, msh, 'XUM input')
      call BFD_transpose(X, Xt, X%pgno)
      call BFD_Lsolve(Ub, 'U', Xt, Yt, X%pgno, option%tol_rand*BFD_tolfac(), option, stats, ptree)
      call BFD_transpose(Yt, res, X%pgno)
      if (BFD_tolfac() < 1d0) call BFD_recompress(res, option%tol_rand, option, stats, ptree)
      call BF_delete(Xt, 1)
      call BF_delete(Yt, 1)
      call BFD_copymeta(res, X)
      call BF_get_rank(res, ptree)
      if (BFD_debug(option)) call BFD_checkleaves(res, msh, 'XUM output')
   end subroutine BFD_XUM

   !>**** res = A*B at the level of the target tmpl (Hmat_add_multiply_Hblock3, at least one operand is a butterfly)
   !>**** res = A*B at the level of tmpl, to be added to the H-block C (Hmat_add_multiply_Hblock3). The product is
   !> truncated relative to the norm of C (BF_TruncRank with scale BFD_scalefac()*||C||): its parts that are negligible
   !> relative to the block it is added to are dropped, however large they are relative to the product itself
   subroutine BFD_Multiply(A, B, tmpl, res, C, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::A, B, tmpl, res, C
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      real(kind=8)::f0, t0, f1, t1, cn, pn, cf
      integer cat, l
      BFD_mul_level = tmpl%level
      call BFD_mstart(stats, f0, t0)
      if (BFD_debug(option) .and. A%style == 2) call BFD_checkleaves(A, msh, 'Multiply input A')
      if (BFD_debug(option) .and. B%style == 2) call BFD_checkleaves(B, msh, 'Multiply input B')
      cn = BFD_scalefac()*BFD_blknorm(C, ptree, stats)
      call BFD_mstart(stats, f1, t1)
      call BFD_Product(A, B, res, tmpl%pgno, option%tol_rand*BFD_tolfac(), option, stats, ptree, cn)
      if (A%style == 2 .and. B%style == 2) then
         cat = 2
      elseif (B%style == 2) then
         cat = 3
      elseif (A%style == 2) then
         cat = 4
      else
         cat = 5
      endif
      call BFD_mstat(cat, stats, f1, t1)
      call BFD_copymeta(res, tmpl)
      call BFD_mstart(stats, f1, t1)
      call BFD_adjustlevel(res, tmpl%level_butterfly, option%tol_rand*BFD_tolfac(), option, stats, ptree, msh, cn)
      call BFD_mstat(6, stats, f1, t1)
      call BFD_mstart(stats, f1, t1)
      if (BFD_tolfac() < 1d0) call BFD_recompress(res, option%tol_rand, option, stats, ptree, cn)
      call BFD_mstat(7, stats, f1, t1)
      call BF_get_rank(res, ptree)
      if (option%verbosity >= 2) then
         f1 = stats%Flop_Tmp
         pn = BFD_normest(res, ptree, stats)
         cf = BFD_normest(C, ptree, stats)
         stats%Flop_Tmp = f1
         l = min(BFD_mul_level, BFD_RD - 1)
         if (C%normest > 0d0 .and. pn > 0d0) then
            BFD_mul_nlog(l) = BFD_mul_nlog(l) + log10(pn/C%normest)
            BFD_mul_nmax(l) = max(BFD_mul_nmax(l), pn/C%normest)
            BFD_mul_nmin(l) = min(BFD_mul_nmin(l), pn/C%normest)
         endif
         if (cf > 0d0) then
            BFD_mul_cmax(l) = max(BFD_mul_cmax(l), C%normest/cf)
            BFD_mul_cmin(l) = min(BFD_mul_cmin(l), C%normest/cf)
         endif
         BFD_mul_prank(l) = max(BFD_mul_prank(l), res%rankmax)
         BFD_mul_prsum(l) = BFD_mul_prsum(l) + res%rankmax
      endif
      if (BFD_debug(option)) call BFD_checkleaves(res, msh, 'Multiply output')
      call BFD_mstat(1, stats, f0, t0)
      BFD_mul_level = -1
   end subroutine BFD_Multiply

   !>**** res = C + P (chara='+') or C - P (chara='-'), C a butterfly (Hmat_BF_add)
   subroutine BFD_AddBF(C, chara, P, res, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::C, P, res
      character chara
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::Pa
      DT beta
      call BFD_check_seq(C, ptree)
      beta = BPACK_cone
      if (chara == '-') beta = -BPACK_cone
      if (P%level_butterfly == C%level_butterfly) then
         call BFD_sum(C, P, beta, res, C%pgno)
      else
         call BFD_copy(P, Pa)
         Pa%pgno = C%pgno
         call BFD_copymeta(Pa, C)
         call BFD_adjustlevel(Pa, C%level_butterfly, option%tol_rand*BFD_tolfac(), option, stats, ptree, msh)
         call BFD_sum(C, Pa, beta, res, C%pgno)
         call BF_delete(Pa, 1)
      endif
      call BFD_recompress(res, option%tol_rand, option, stats, ptree)
      call BFD_copymeta(res, C)
      if (BFD_debug(option)) then
         call BFD_checkleaves(C, msh, 'AddBF input C')
         call BFD_checkleaves(P, msh, 'AddBF input P')
         call BFD_checkleaves(res, msh, 'AddBF output')
      endif
   end subroutine BFD_AddBF

   !>**** res = C - A*B (chara='-') or C + A*B, C a butterfly and A, B partitioned H-blocks (Hmat_add_multiply): C is split
   !> into its children at their native level, updated recursively and merged back
   subroutine BFD_SubHH(C, chara, A, B, res, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::C, res
      type(matrixblock)::A, B
      character chara
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      call BFD_SubHH_rec(C, chara, A, B, res, option%tol_rand*BFD_tolfac(), BFD_scalefac()*BFD_blknorm(C, ptree, stats), option, stats, ptree, msh)
      if (BFD_tolfac() < 1d0) call BFD_recompress(res, option%tol_rand, option, stats, ptree)
      if (BFD_debug(option)) then
         call BFD_checkleaves(C, msh, 'SubHH input C')
         call BFD_checkleaves(res, msh, 'SubHH output')
      endif
   end subroutine BFD_SubHH

   !>**** the products are truncated relative to cn, the scaled norm of the top-level block C (see BFD_Multiply)
   recursive subroutine BFD_SubHH_rec(C, chara, A, B, res, tol, cn, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::C, res
      type(matrixblock), target::A, B
      character chara
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock), target::Cw
      type(matrixblock)::T, Pr
      type(matrixblock), pointer::Aik, Bkj, Cij
      DT beta
      integer i, j, k
      real(kind=8)::tol, cn
      call BFD_check_seq(C, ptree)
      beta = BPACK_cone
      if (chara == '-') beta = -BPACK_cone
      call BF_copy('N', C, Cw)
      call BFD_SplitNative(Cw, option, stats, ptree, msh)
      do i = 1, 2
         do j = 1, 2
            Cij => Cw%sons(i, j)
            do k = 1, 2
               Aik => A%sons(i, k)
               Bkj => B%sons(k, j)
               if (Aik%style == 4 .and. Bkj%style == 4) then
                  call BFD_SubHH_rec(Cij, chara, Aik, Bkj, T, tol, cn, option, stats, ptree, msh)
                  call BFD_Install(Cij, T, ptree)
               else
                  call BFD_Product(Aik, Bkj, Pr, Cij%pgno, tol, option, stats, ptree, cn)
                  call BFD_copymeta(Pr, Cij)
                  call BFD_adjustlevel(Pr, Cij%level_butterfly, tol, option, stats, ptree, msh, cn)
                  call BFD_sum(Cij, Pr, beta, T, Cij%pgno)
                  call BF_delete(Pr, 1)
                  call BFD_recompress(T, tol, option, stats, ptree)
                  call BFD_Install(Cij, T, ptree)
               endif
            enddo
         enddo
      enddo
      call BFD_MergeNative(Cw, C, res, tol, option, stats, ptree, msh)
      do i = 1, 2
         do j = 1, 2
            call BF_delete(Cw%sons(i, j), 1)
         enddo
      enddo
      deallocate (Cw%sons)
      call BF_delete(Cw, 1)
   end subroutine BFD_SubHH_rec

!======================================================================================== distributed HODBF (bf_algebra=1)

   !>**** A = A + beta*B in place, for two (possibly distributed) butterflies on the same process group with the same levels
   !> and leaves: U side by side, kernels block diagonal, V side by side. Every process only touches its local blocks.
   subroutine BFD_sumdist(A, B, beta, stats, ptree)
      implicit none
      type(matrixblock)::A, B
      DT beta
      type(Hstat)::stats
      type(proctree)::ptree
      type(butterflymatrix)::T
      integer K, l, i, j
      real(kind=8)::t0

      t0 = MPI_Wtime()
      K = A%level_butterfly
      call assert(B%level_butterfly == K .and. B%pgno == A%pgno, 'BFD_sumdist: level or process group mismatch')
      if (K > 0 .and. B%level_half /= A%level_half) call BF_ChangePattern(B, BFD_pattern(B), BFD_pattern(A), stats, ptree)
      if (IOwnPgrp(ptree, A%pgno)) then
         call assert(A%ButterflyU%nblk_loc == B%ButterflyU%nblk_loc .and. A%ButterflyV%nblk_loc == B%ButterflyV%nblk_loc, 'BFD_sumdist: layout mismatch')
         do i = 1, A%ButterflyU%nblk_loc
            call BFD_hcat(A%ButterflyU%blocks(i)%matrix, B%ButterflyU%blocks(i)%matrix, beta, T)
            call BFD_matmove(T, A%ButterflyU%blocks(i))
         enddo
         do j = 1, A%ButterflyV%nblk_loc
            call BFD_hcat(A%ButterflyV%blocks(j)%matrix, B%ButterflyV%blocks(j)%matrix, BPACK_cone, T)
            call BFD_matmove(T, A%ButterflyV%blocks(j))
         enddo
         do l = 1, K
            if (A%ButterflyKerl(l)%nr > 0 .and. A%ButterflyKerl(l)%nc > 0) then
               call assert(B%ButterflyKerl(l)%nr == A%ButterflyKerl(l)%nr .and. B%ButterflyKerl(l)%nc == A%ButterflyKerl(l)%nc &
                           .and. B%ButterflyKerl(l)%idx_r == A%ButterflyKerl(l)%idx_r .and. B%ButterflyKerl(l)%idx_c == A%ButterflyKerl(l)%idx_c, &
                           'BFD_sumdist: kernel layout mismatch')
               do j = 1, A%ButterflyKerl(l)%nc
                  do i = 1, A%ButterflyKerl(l)%nr
                     call BFD_blkdiag(A%ButterflyKerl(l)%blocks(i, j)%matrix, B%ButterflyKerl(l)%blocks(i, j)%matrix, T)
                     call BFD_matmove(T, A%ButterflyKerl(l)%blocks(i, j))
                  enddo
               enddo
            endif
         enddo
      endif
      call BF_get_rank(A, ptree)
      call BFD_tadd(12, t0)
   end subroutine BFD_sumdist

   !>**** products with the chain alpha*F(1)*...*F(n) (BMatVec interface, operand = bfd_chainop)
   subroutine BFD_ChainRand_MVP(operand, block_o, trans, M, N, num_vect_sub, Vin, ldi, Vout, ldo, a, b, ptree, stats, operand1)
      implicit none
      class(*)::operand
      class(*), optional::operand1
      type(matrixblock)::block_o
      character trans
      integer M, N, num_vect_sub, ldi, ldo
      type(proctree)::ptree
      type(Hstat)::stats
      DT :: Vin(ldi, *), Vout(ldo, *), a, b
      DT, allocatable::v(:, :), w(:, :)
      type(matrixblock), pointer::F
      integer k, nv, mi, mo

      select type (operand)
      type is (bfd_chainop)
         nv = num_vect_sub
         if (trans == 'N') then
            mi = operand%f(operand%n)%ptr%N_loc
            call assert(mi == N, 'BFD_ChainRand_MVP: input size mismatch')
         else
            mi = operand%f(1)%ptr%M_loc
            call assert(mi == M, 'BFD_ChainRand_MVP: input size mismatch')
         endif
         allocate (v(max(mi, 1), nv))
         v = 0
         if (mi > 0) v(1:mi, :) = Vin(1:mi, 1:nv)
         do k = 1, operand%n
            if (trans == 'N') then
               F => operand%f(operand%n - k + 1)%ptr
               call assert(F%N_loc == mi, 'BFD_ChainRand_MVP: nonconforming factors')
               mo = F%M_loc
            else
               F => operand%f(k)%ptr
               call assert(F%M_loc == mi, 'BFD_ChainRand_MVP: nonconforming factors')
               mo = F%N_loc
            endif
            allocate (w(max(mo, 1), nv))
            w = 0
            call BF_block_MVP_dat(F, trans, F%M_loc, F%N_loc, nv, v, max(mi, 1), w, max(mo, 1), BPACK_cone, BPACK_czero, ptree, stats)
            if ((trans == 'N' .and. operand%iplus(operand%n - k + 1)) .or. (trans /= 'N' .and. operand%iplus(k))) then
               call assert(mo == mi, 'BFD_ChainRand_MVP: I+F needs a square factor')
               if (mo > 0) w(1:mo, :) = w(1:mo, :) + v(1:mo, :)
            endif
            call move_alloc(w, v)
            mi = mo
         enddo
         if (trans == 'N') then
            call assert(mi == M, 'BFD_ChainRand_MVP: output size mismatch')
         else
            call assert(mi == N, 'BFD_ChainRand_MVP: output size mismatch')
         endif
         if (mi > 0) Vout(1:mi, 1:nv) = a*operand%alpha*v(1:mi, :) + b*Vout(1:mi, 1:nv)
         deallocate (v)
      class default
         write (*, *) "unexpected type in BFD_ChainRand_MVP"
         stop
      end select
   end subroutine BFD_ChainRand_MVP

   !>**** set op to alpha*F1*F2*... (up to five factors, iplus marks the factors that enter as I+F)
   subroutine BFD_chainset(op, alpha, iplus, F1, F2, F3, F4, F5)
      implicit none
      type(bfd_chainop)::op
      DT alpha
      logical iplus(:)
      type(matrixblock), target::F1
      type(matrixblock), target, optional::F2, F3, F4, F5
      op%alpha = alpha
      op%n = 1
      op%f(1)%ptr => F1
      if (present(F2)) then
         op%n = 2
         op%f(2)%ptr => F2
      endif
      if (present(F3)) then
         op%n = 3
         op%f(3)%ptr => F3
      endif
      if (present(F4)) then
         op%n = 4
         op%f(4)%ptr => F4
      endif
      if (present(F5)) then
         op%n = 5
         op%f(5)%ptr => F5
      endif
      call assert(size(iplus) == op%n, 'BFD_chainset: iplus size mismatch')
      op%iplus(1:op%n) = iplus
   end subroutine BFD_chainset

   !>**** Z = X*Y (+ W if present) by the distributed deterministic BMult (and a sum), truncated with tol
   subroutine BFD_mul(X, Y, Z, tol, option, stats, ptree, msh, W)
      implicit none
      type(matrixblock)::X, Y, Z
      type(matrixblock), optional::W
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      if (present(W)) then
         call BFD_bmult_dist(X, Y, Z, tol, option, stats, ptree, msh, norecomp=.true.)
         call BFD_sumdist(Z, W, BPACK_cone, stats, ptree)
         call BFD_recompress(Z, tol, option, stats, ptree)
      else
         call BFD_bmult_dist(X, Y, Z, tol, option, stats, ptree, msh)
      endif
   end subroutine BFD_mul

   subroutine BFD_scale(Z, alpha)
      implicit none
      type(matrixblock)::Z
      DT alpha
      integer k
      do k = 1, Z%ButterflyU%nblk_loc
         Z%ButterflyU%blocks(k)%matrix = alpha*Z%ButterflyU%blocks(k)%matrix
      enddo
   end subroutine BFD_scale

   !>**** replace the content of blk by Z (Z is emptied); blk keeps its cluster level and groups
   subroutine BFD_install_dist(Z, blk, ptree)
      implicit none
      type(matrixblock)::Z, blk
      type(proctree)::ptree
      integer lvl, rg, cg
      lvl = blk%level
      rg = blk%row_group
      cg = blk%col_group
      call BF_delete(blk, 1)
      call BF_copy_delete(Z, blk)
      blk%level = lvl
      blk%row_group = rg
      blk%col_group = cg
      call BF_get_rank(blk, ptree)
   end subroutine BFD_install_dist

   !>**** gather one leaf (U or V) block that is row-distributed over the sharing group pgno_sub on its head; the other
   !> processes drop it (nblk_loc=0). Same data movement as the input preparation in BF_split.
   subroutine BFD_leaf_to_head(leaves, gg, pgno_sub, ptree, msh)
      implicit none
      type(butterfly_UV)::leaves
      integer gg, pgno_sub
      type(proctree)::ptree
      type(mesh)::msh
      integer nsub, mloc, rank
      integer, allocatable::M_p_sub(:, :), M_p_sub1(:, :)
      DT, allocatable::t1(:, :), t2(:, :)

      nsub = ptree%pgrp(pgno_sub)%nproc
      allocate (M_p_sub(nsub, 2), M_p_sub1(nsub, 2))
      call ComputeParallelIndicesSub(gg, pgno_sub, ptree, msh, M_p_sub)
      M_p_sub1 = -1
      M_p_sub1(1, 1) = 1
      M_p_sub1(1, 2) = msh%basis_group(gg)%tail - msh%basis_group(gg)%head + 1
      mloc = size(leaves%blocks(1)%matrix, 1)
      rank = size(leaves%blocks(1)%matrix, 2)
      allocate (t1(max(mloc, 1), rank), t2(M_p_sub1(1, 2), rank))
      t1 = 0
      if (mloc > 0) t1(1:mloc, :) = leaves%blocks(1)%matrix
      call Redistribute1Dto1D(t1, max(mloc, 1), M_p_sub, 0, pgno_sub, t2, M_p_sub1(1, 2), M_p_sub1, 0, pgno_sub, rank, ptree)
      deallocate (leaves%blocks(1)%matrix)
      if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
         allocate (leaves%blocks(1)%matrix(M_p_sub1(1, 2), rank))
         leaves%blocks(1)%matrix = t2
      else
         deallocate (leaves%blocks)
         leaves%nblk_loc = 0
      endif
      deallocate (t1, t2, M_p_sub, M_p_sub1)
   end subroutine BFD_leaf_to_head

   !>**** for a butterfly with more processes than leaves, gather every shared leaf on the head of its sharing group, the
   !> layout from which BF_all2all_UV (hence BF_ReDistribute_Inplace) can move it
   subroutine BFD_leaves_to_head(blk, ptree, msh)
      implicit none
      type(matrixblock)::blk
      type(proctree)::ptree
      type(mesh)::msh
      integer pgno_sub, gg
      if (blk%level_butterfly == 0 .or. .not. IOwnPgrp(ptree, blk%pgno)) return
      call GetPgno_Sub(ptree, blk%pgno, blk%level_butterfly, pgno_sub)
      if (ptree%pgrp(pgno_sub)%nproc == 1) return
      gg = blk%row_group*2**blk%level_butterfly + blk%ButterflyU%idx - 1
      call BFD_leaf_to_head(blk%ButterflyU, gg, pgno_sub, ptree, msh)
      gg = blk%col_group*2**blk%level_butterfly + blk%ButterflyV%idx - 1
      call BFD_leaf_to_head(blk%ButterflyV, gg, pgno_sub, ptree, msh)
   end subroutine BFD_leaves_to_head

   !>**** one off-diagonal chain of the partitioned inverse: T = X(I+P), then Y = -(I+Q)T; T is deleted if freeT.
   !> X, P and Q are only read (X is the added term, whose pattern may change).
   subroutine BFD_offdiag_chain(X, P, Q, T, Y, freeT, tol, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::X, P, Q, T, Y
      logical freeT
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      call BFD_mul(X, P, T, tol, option, stats, ptree, msh, X)
      call BFD_mul(Q, T, Y, tol, option, stats, ptree, msh, T)
      call BFD_scale(Y, -BPACK_cone)
      if (freeT) call BF_delete(T, 1)
   end subroutine BFD_offdiag_chain

   !>**** (I+C)^{-1} = I+C' in place by partitioned inversion with the deterministic butterfly algebra: split into the four
   !> (K-2)-level children (moved to a group of at most 2^(K-2) processes, so that every block owns whole leaves), D' and
   !> A' recursively with the Schur complement A = A - B(C + D'C) in between, then T = B(I+D'), Y12 = -(I+A')T,
   !> T2 = C(I+A'), Y21 = -(I+D')T2, Y22 = D' - Y21*T, and the aggregation of [A' Y12; Y21 Y22] (BF_Aggregate for K >= 2,
   !> BFD_aggregate_lr1_dist for K = 1) followed by recompression. Products are distributed BMults; intermediate results
   !> are truncated with tol_rand*tolfac, the result with tol_rand. K = 0 uses LR_SMW.
   recursive subroutine BFD_IplusInverse(blocks_io, recurlevel, option, stats, ptree, msh, pgno)
      implicit none
      type(matrixblock)::blocks_io
      integer recurlevel, pgno
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::partitioned_block, Z, Y12, Y21, Y22, T1, T2
      type(matrixblock), target::Cin
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      type(bfd_chainop)::op
      logical dbg, conc
      real(kind=8)::Memory, error, n1, n2, tol, tolin
      integer ii, jj, ierr, K, Kc, pgno_agg, nth, lvl
      DT::phaseA, phaseD
      DTR::ldA, ldD

      tol = option%tol_rand
      tolin = tol*BFD_tolfac()
      K = blocks_io%level_butterfly
      if (K == 0) then
         call LR_SMW(blocks_io, Memory, ptree, option, stats, pgno)
         return
      endif
      dbg = BFD_debug(option) .and. recurlevel > 0
      if (dbg) call BF_copy('N', blocks_io, Cin)

      allocate (partitioned_block%sons(2, 2))
      blocks_A => partitioned_block%sons(1, 1)
      blocks_B => partitioned_block%sons(1, 2)
      blocks_C => partitioned_block%sons(2, 1)
      blocks_D => partitioned_block%sons(2, 2)
      n1 = MPI_Wtime()
      call BF_split(blocks_io, partitioned_block, ptree, stats, msh, option)
      n2 = MPI_Wtime()
      stats%Time_split = stats%Time_split + n2 - n1
      Kc = blocks_A%level_butterfly

      !>**** BF_split may leave the children on more processes than leaves: move them to a group of at most 2^Kc processes
      pgno_agg = blocks_A%pgno
      do while (ptree%pgrp(pgno_agg)%nproc > 2**Kc)
         pgno_agg = pgno_agg*2
      enddo
      if (pgno_agg /= blocks_A%pgno) then
         if (IOwnPgrp(ptree, blocks_A%pgno)) then
            do jj = 1, 2
               do ii = 1, 2
                  call BFD_leaves_to_head(partitioned_block%sons(ii, jj), ptree, msh)
                  call BF_ReDistribute_Inplace(partitioned_block%sons(ii, jj), pgno_agg, stats, ptree, msh)
               enddo
            enddo
         endif
         do jj = 1, 2 ! every process of blocks_io%pgno must see the new process group and layouts (as after BF_split)
            do ii = 1, 2
               partitioned_block%sons(ii, jj)%pgno = pgno_agg
               if (.not. IOwnPgrp(ptree, pgno_agg)) call ComputeParallelIndices(partitioned_block%sons(ii, jj), pgno_agg, ptree, msh)
            enddo
         enddo
      endif

      if (IOwnPgrp(ptree, blocks_D%pgno)) then
         call BFD_IplusInverse(blocks_D, recurlevel + 1, option, stats, ptree, msh, blocks_D%pgno)
         if (Kc == 0) then
            call LR_A_minusBDinvC(partitioned_block, ptree, option, stats)
         else
            !>**** A = A - B(C + D'C)
            call BFD_mul(blocks_D, blocks_C, T1, tolin, option, stats, ptree, msh, blocks_C)
            call BFD_bmult_dist(blocks_B, T1, Z, tolin, option, stats, ptree, msh, norecomp=.true.)
            call BF_delete(T1, 1)
            call BFD_sumdist(blocks_A, Z, -BPACK_cone, stats, ptree)
            call BF_delete(Z, 1)
            call BFD_recompress(blocks_A, tolin, option, stats, ptree)
         endif
         call BFD_IplusInverse(blocks_A, recurlevel + 1, option, stats, ptree, msh, blocks_D%pgno)
         phaseA = blocks_A%phase
         phaseD = blocks_D%phase
         ldA = blocks_A%logabsdet
         ldD = blocks_D%logabsdet
         !>**** T = B(I+D'), Y12 = -(I+A')T, T2 = C(I+A'), Y21 = -(I+D')T2, Y22 = D' - Y21*T. The chains T -> Y12 and
         !> T2 -> Y21 are independent: on one process they run as two concurrent OpenMP tasks, each of which uses half of
         !> the threads for its own parallel loops (one more active nesting level).
         conc = .false.
#ifdef HAVE_OPENMP
         nth = omp_get_max_threads()
         conc = ptree%pgrp(blocks_D%pgno)%nproc == 1 .and. nth > 1 .and. .not. omp_in_parallel()
#endif
         if (conc) then
#ifdef HAVE_OPENMP
            lvl = omp_get_max_active_levels()
            call omp_set_max_active_levels(max(lvl, 2))
            !$omp parallel num_threads(2) default(shared)
            !$omp single
            !$omp task default(shared)
            call omp_set_num_threads(max(1, nth/2))
            call BFD_offdiag_chain(blocks_B, blocks_D, blocks_A, T1, Y12, .false., tolin, option, stats, ptree, msh)
            !$omp end task
            !$omp task default(shared)
            call omp_set_num_threads(max(1, nth - nth/2))
            call BFD_offdiag_chain(blocks_C, blocks_A, blocks_D, T2, Y21, .true., tolin, option, stats, ptree, msh)
            !$omp end task
            !$omp end single
            !$omp end parallel
            call omp_set_max_active_levels(lvl)
#endif
         else
            call BFD_offdiag_chain(blocks_B, blocks_D, blocks_A, T1, Y12, .false., tolin, option, stats, ptree, msh)
            call BFD_offdiag_chain(blocks_C, blocks_A, blocks_D, T2, Y21, .true., tolin, option, stats, ptree, msh)
         endif
         call BFD_bmult_dist(Y21, T1, Y22, tolin, option, stats, ptree, msh, norecomp=.true.)
         call BF_delete(T1, 1)
         call BFD_scale(Y22, -BPACK_cone)
         call BFD_sumdist(Y22, blocks_D, BPACK_cone, stats, ptree)
         call BFD_recompress(Y22, tolin, option, stats, ptree)
         call BFD_install_dist(Y12, blocks_B, ptree)
         call BFD_install_dist(Y21, blocks_C, ptree)
         call BFD_install_dist(Y22, blocks_D, ptree)
      else
         phaseA = 1
         phaseD = 1
         ldA = 0
         ldD = 0
      endif
      call MPI_Bcast(phaseA, 1, MPI_DT, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
      call MPI_Bcast(ldA, 1, MPI_DTR, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
      call MPI_Bcast(phaseD, 1, MPI_DT, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
      call MPI_Bcast(ldD, 1, MPI_DTR, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)

      if (K == 1) then
         call BFD_aggregate_lr1_dist(partitioned_block, blocks_io, pgno_agg, tol, option, stats, ptree, msh)
      else
         call BF_Aggregate(partitioned_block, blocks_io, ptree, stats, option, msh)
      endif
      if (IOwnPgrp(ptree, blocks_io%pgno)) call BFD_recompress(blocks_io, tol, option, stats, ptree)

      blocks_io%phase = phaseA*phaseD
      blocks_io%logabsdet = ldA + ldD
      if (dbg) then
         call BFD_chainset(op, BPACK_cone, [.true., .true.], Cin, blocks_io)
         error = BFD_chainerr_identity(op, blocks_io, ptree, stats)
         if (ptree%MyID == ptree%pgrp(blocks_io%pgno)%head) write (*, '(A,I3,A,I3,A,I4,A,I4,A,ES10.3)') '   BFD (I+C)^-1 recurlevel ', recurlevel, ' K ', K, &
            ' nproc ', ptree%pgrp(blocks_io%pgno)%nproc, ' nproc_children ', ptree%pgrp(pgno_agg)%nproc, ' err ', error
         call BF_delete(Cin, 1)
      endif
      if (option%verbosity >= 2 .and. recurlevel == 0 .and. ptree%MyID == Main_ID) write (*, '(A23,A6,I3,A8,I3)') ' RecursiveI (det) ', ' rank:', blocks_io%rankmax, ' L_butt:', blocks_io%level_butterfly

      do jj = 1, 2
         do ii = 1, 2
            call BF_delete(partitioned_block%sons(ii, jj), 1)
         enddo
      enddo
      deallocate (partitioned_block%sons)
   end subroutine BFD_IplusInverse

   !>**** four low-rank children on the single process pgno_c -> 1-level butterfly (BFD_aggregate_lr1), moved to the
   !> process group of blocks_io and installed there
   subroutine BFD_aggregate_lr1_dist(part, blocks_io, pgno_c, tol, option, stats, ptree, msh)
      implicit none
      type(matrixblock)::part, blocks_io
      integer pgno_c
      real(kind=8)::tol
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock)::O
      if (IOwnPgrp(ptree, pgno_c)) then
         call assert(ptree%pgrp(pgno_c)%nproc == 1, 'BFD_aggregate_lr1_dist: the low-rank children must be on one process')
         call BFD_aggregate_lr1(part, O, pgno_c, tol, option, stats)
      endif
      O%style = 2
      O%level_butterfly = 1
      O%level_half = 0
      O%pgno = pgno_c
      O%pgno_db = pgno_c
      O%level = blocks_io%level
      O%row_group = blocks_io%row_group
      O%col_group = blocks_io%col_group
      O%headm = blocks_io%headm
      O%headn = blocks_io%headn
      O%M = blocks_io%M
      O%N = blocks_io%N
      if (pgno_c /= blocks_io%pgno) then
         call BF_ReDistribute_Inplace(O, blocks_io%pgno, stats, ptree, msh)
      else
         call ComputeParallelIndices(O, O%pgno, ptree, msh)
      endif
      call BFD_install_dist(O, blocks_io, ptree)
   end subroutine BFD_aggregate_lr1_dist

   !>**** HODBF with bf_algebra=1: Schur complement I+C, C = -B12'*B21' (distributed BMult), of the node rowblock at
   !> level_c and its inverse I+C' (replaces BF_inverse_schur_partitionedinverse). With verbosity >= 3 the result is
   !> checked by ||(I+C)(I+C')x - x||/||x|| on random vectors.
   subroutine BFD_inverse_schur_partitionedinverse(ho_bf1, level_c, rowblock, option, stats, ptree, msh)
      implicit none
      type(hobf)::ho_bf1
      integer level_c, rowblock
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock), pointer::block_o, block_off1, block_off2
      type(matrixblock), target::Csave
      type(matrixblock)::Zd
      type(bfd_chainop)::op
      real(kind=8)::n1, n2, err
      logical dbg

      block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)
      block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
      block_o%level_butterfly = block_off1%level_butterfly

      if (block_off1%level_butterfly == 0 .or. block_off2%level_butterfly == 0) then
         call BFD_flops_flush(stats)
         call LR_minusBC(ho_bf1, level_c, rowblock, ptree, stats)
         stats%Flop_Tmp = 0 ! LR_minusBC already added its flops to stats%Flop_Factor
      else
         call BFD_bmult_dist(block_off1, block_off2, Zd, option%tol_rand, option, stats, ptree, msh)
         call BFD_scale(Zd, -BPACK_cone)
         call BFD_install_dist(Zd, block_o, ptree)
      endif
      call BF_get_rank(block_o, ptree)

      dbg = BFD_debug(option)
      if (dbg) call BF_copy('N', block_o, Csave)
      n1 = MPI_Wtime()
      call BFD_IplusInverse(block_o, 0, option, stats, ptree, msh, block_o%pgno)
      n2 = MPI_Wtime()
      stats%Time_SMW = stats%Time_SMW + n2 - n1
      call BFD_tadd(5, n1)
      if (dbg) then
         call BFD_chainset(op, BPACK_cone, [.true., .true.], Csave, block_o)
         err = BFD_chainerr_identity(op, block_o, ptree, stats)
         if (ptree%MyID == ptree%pgrp(block_o%pgno)%head) write (*, '(A,I3,A,I5,A,I3,A,I4,A,ES10.3)') ' BFD (I+C)^-1 level ', level_c, ' blk ', rowblock, &
            ' K ', block_o%level_butterfly, ' rank ', block_o%rankmax, ' ||(I+C)(I+C^-1-I)x-x||/||x|| ', err
         call BF_delete(Csave, 1)
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) write (*, '(A10,I5,A6,I3,A8,I3)') 'OneL No. ', rowblock, ' rank:', block_o%rankmax, ' L_butt:', block_o%level_butterfly
   end subroutine BFD_inverse_schur_partitionedinverse

   !>**** ||op*x - x||/||x|| for random x (op must be square)
   real(kind=8) function BFD_chainerr_identity(op, tmpl, ptree, stats)
      implicit none
      type(bfd_chainop)::op
      type(matrixblock)::tmpl
      type(proctree)::ptree
      type(Hstat)::stats
      DT, allocatable::x(:, :), y(:, :)
      real(kind=8), allocatable::xr(:, :)
      real(kind=8)::nrm(2)
      integer m, nvec, ierr
      nvec = 4
      m = tmpl%M_loc
      allocate (x(max(m, 1), nvec), y(max(m, 1), nvec), xr(max(m, 1), nvec))
      call random_number(xr)
      x = xr - 0.5d0
      y = 0
      call BFD_ChainRand_MVP(op, tmpl, 'N', m, m, nvec, x, max(m, 1), y, max(m, 1), BPACK_cone, BPACK_czero, ptree, stats)
      nrm = 0
      if (m > 0) then
         nrm(1) = sum(abs(y(1:m, :) - x(1:m, :))**2)
         nrm(2) = sum(abs(x(1:m, :))**2)
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, nrm, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(tmpl%pgno)%Comm, ierr)
      BFD_chainerr_identity = sqrt(nrm(1)/max(nrm(2), BPACK_SafeUnderflow))
      deallocate (x, y, xr)
   end function BFD_chainerr_identity

!======================================================================================== deterministic HODBF Sblock (bf_algebra=1)

   !>**** layout of a butterfly whose columns are not a cluster group (a partial): M_loc, N_loc from the local leaves,
   !> M_p/N_p by an allgather over its process group, ms/ns local accumulations (as at the end of BF_extract_partial)
   subroutine BFD_setlayout_nc(blk, ptree)
      implicit none
      type(matrixblock)::blk
      type(proctree)::ptree
      integer np, ierr, pp, k
      integer, allocatable::mloc(:), nloc(:)
      np = ptree%pgrp(blk%pgno)%nproc
      blk%M_loc = 0
      blk%N_loc = 0
      if (associated(blk%ms)) deallocate (blk%ms)
      if (associated(blk%ns)) deallocate (blk%ns)
      allocate (blk%ms(blk%ButterflyU%nblk_loc), blk%ns(blk%ButterflyV%nblk_loc))
      do k = 1, blk%ButterflyU%nblk_loc
         blk%M_loc = blk%M_loc + size(blk%ButterflyU%blocks(k)%matrix, 1)
         blk%ms(k) = blk%M_loc
      enddo
      do k = 1, blk%ButterflyV%nblk_loc
         blk%N_loc = blk%N_loc + size(blk%ButterflyV%blocks(k)%matrix, 1)
         blk%ns(k) = blk%N_loc
      enddo
      allocate (mloc(np), nloc(np))
      if (np == 1) then ! no collective on one process
         mloc(1) = blk%M_loc
         nloc(1) = blk%N_loc
      else
         call MPI_ALLGATHER(blk%M_loc, 1, MPI_INTEGER, mloc, 1, MPI_INTEGER, ptree%pgrp(blk%pgno)%Comm, ierr)
         call MPI_ALLGATHER(blk%N_loc, 1, MPI_INTEGER, nloc, 1, MPI_INTEGER, ptree%pgrp(blk%pgno)%Comm, ierr)
      endif
      if (associated(blk%M_p)) deallocate (blk%M_p)
      if (associated(blk%N_p)) deallocate (blk%N_p)
      allocate (blk%M_p(np, 2), blk%N_p(np, 2))
      do pp = 1, np
         if (pp == 1) then
            blk%M_p(pp, 1) = 1
            blk%N_p(pp, 1) = 1
         else
            blk%M_p(pp, 1) = blk%M_p(pp - 1, 2) + 1
            blk%N_p(pp, 1) = blk%N_p(pp - 1, 2) + 1
         endif
         blk%M_p(pp, 2) = blk%M_p(pp, 1) + mloc(pp) - 1
         blk%N_p(pp, 2) = blk%N_p(pp, 1) + nloc(pp) - 1
      enddo
      blk%M = sum(mloc)
      blk%N = sum(nloc)
      deallocate (mloc, nloc)
   end subroutine BFD_setlayout_nc

   !>**** move a butterfly with non-cluster columns to another (nested) process group, as BF_ReDistribute_Inplace
   !> but with the layout recomputed from the local leaves. Called by every process of the larger group; the processes
   !> without data must know level_butterfly, level_half and the old pgno.
   subroutine BFD_redist_nc(blk, pgno_new, stats, ptree)
      implicit none
      type(matrixblock)::blk
      integer pgno_new, level, K
      type(Hstat)::stats
      type(proctree)::ptree
      if (blk%pgno == pgno_new) return
      K = blk%level_butterfly
      if (.not. allocated(blk%ButterflyKerl)) allocate (blk%ButterflyKerl(K))
      do level = 0, K + 1
         if (level == 0) then
            call BF_all2all_UV(blk, blk%pgno, blk%ButterflyV, level, 0, blk, pgno_new, blk%ButterflyV, level, stats, ptree)
         elseif (level == K + 1) then
            call BF_all2all_UV(blk, blk%pgno, blk%ButterflyU, level, 0, blk, pgno_new, blk%ButterflyU, level, stats, ptree)
         else
            call BF_all2all_ker(blk, blk%pgno, blk%ButterflyKerl(level), level, 0, 0, blk, pgno_new, blk%ButterflyKerl(level), level, stats, ptree)
         endif
      enddo
      blk%pgno = pgno_new
      blk%pgno_db = pgno_new
      if (IOwnPgrp(ptree, pgno_new)) then
         call BFD_setlayout_nc(blk, ptree)
         call BF_get_rank(blk, ptree)
      else
         deallocate (blk%ButterflyKerl)
         blk%ButterflyU%nblk_loc = 0
         blk%ButterflyV%nblk_loc = 0
         if (allocated(blk%ButterflyU%blocks)) deallocate (blk%ButterflyU%blocks)
         if (allocated(blk%ButterflyV%blocks)) deallocate (blk%ButterflyV%blocks)
      endif
   end subroutine BFD_redist_nc

   !>**** P <- F^{-1} P for the partial P (rows of the HODBF node ii at level, V = identity), F^{-1} the node's factor:
   !> y1 = (I+C')(x1 - B12' x2), y2 = x2 - B21' y1. With one level, F^{-1} is applied to the multi-vector blkdiag(U1,U2)
   !> by BF_block_MVP_inverse_dat and the level-1 kernel is stacked; otherwise the row halves X_a = X1, X2 (extracted on
   !> the child groups, columns padded to [nodes of half 1; nodes of half 2]) are updated by distributed BMults on the
   !> node's group and merged back (U and kernels 2..K by BF_copyback_partial, kernel 1 = V_a^T * [W1; W2]).
   subroutine BFD_apply_factor(ho_bf1, level, ii, P, option, stats, ptree, msh)
      implicit none
      type(hobf)::ho_bf1
      integer level, ii
      type(matrixblock)::P
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock), pointer::B12, B21, Cs
      type(matrixblock)::X(2), Tm, Y1, Pold
      type(bfd_items)::Lout, Lin
      type(butterflymatrix), allocatable::Vr(:, :)
      type(butterflymatrix)::Wb
      DT, allocatable::Xv(:, :), Yv(:, :)
      integer pgno, KP, Kh, pg(2), grp(2), a, k, j, jf, e, r, n, ierr, rr(2), mm, nv, off, m, t, jj, pgno_sub, dst, nc
      integer, allocatable::rsz(:, :)
      real(kind=8)::tolin, err, t0
      logical dbg

      B12 => ho_bf1%levels(level)%BP_inverse_update(2*ii - 1)%LL(1)%matrices_block(1)
      B21 => ho_bf1%levels(level)%BP_inverse_update(2*ii)%LL(1)%matrices_block(1)
      Cs => ho_bf1%levels(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)
      pgno = P%pgno
      KP = P%level_butterfly
      Kh = KP - 1
      tolin = option%tol_rand*BFD_tolfac()
      dbg = BFD_debug(option)
      if (dbg) call BF_copy('N', P, Pold)

      t0 = MPI_Wtime()
      if (KP == 1) then
         !>**** one level: F^{-1} applied to blkdiag(U1, U2), both halves' columns side by side
         rr = 0
         do k = 1, P%ButterflyU%nblk_loc
            t = (k - 1)*P%ButterflyU%inc + P%ButterflyU%idx
            rr(t) = size(P%ButterflyU%blocks(k)%matrix, 2)
         enddo
         if (ptree%pgrp(pgno)%nproc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, rr, 2, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
         nv = rr(1) + rr(2)
         mm = P%M_loc
         allocate (Xv(max(mm, 1), nv), Yv(max(mm, 1), nv))
         Xv = 0
         Yv = 0
         off = 0
         do k = 1, P%ButterflyU%nblk_loc
            t = (k - 1)*P%ButterflyU%inc + P%ButterflyU%idx
            m = size(P%ButterflyU%blocks(k)%matrix, 1)
            Xv(off + 1:off + m, (t - 1)*rr(1) + 1:(t - 1)*rr(1) + rr(t)) = P%ButterflyU%blocks(k)%matrix
            off = off + m
         enddo
         call BF_block_MVP_inverse_dat(ho_bf1, level, ii, 'N', mm, nv, Xv, max(mm, 1), Yv, max(mm, 1), ptree, stats)
         off = 0
         do k = 1, P%ButterflyU%nblk_loc
            m = size(P%ButterflyU%blocks(k)%matrix, 1)
            call BFD_matset(P%ButterflyU%blocks(k), Yv(off + 1:off + m, :))
            off = off + m
         enddo
         deallocate (Xv, Yv)
         if (allocated(P%ButterflyKerl(1)%blocks)) then
            do jj = 1, P%ButterflyKerl(1)%nc
               call BFD_vcat2(P%ButterflyKerl(1)%blocks(1, jj)%matrix, P%ButterflyKerl(1)%blocks(2, jj)%matrix, Wb)
               call BFD_matset(P%ButterflyKerl(1)%blocks(1, jj), Wb%matrix)
               call BFD_matset(P%ButterflyKerl(1)%blocks(2, jj), Wb%matrix)
            enddo
            call BFD_matfree(Wb)
         endif
         call BF_get_rank(P, ptree)
         call BFD_tadd(3, t0)
      else
         !>**** 1. row halves on the child groups
         if (ptree%pgrp(pgno)%nproc > 1) then
            pg = [pgno*2, pgno*2 + 1]
         else
            pg = [pgno, pgno]
         endif
         do a = 1, 2
            grp(a) = P%row_group*2 + a - 1
            X(a)%level_butterfly = Kh
            X(a)%level_half = 0
            X(a)%pgno = pg(a)
            if (IOwnPgrp(ptree, pg(a))) call BF_extract_partial(P, Kh, a, msh%basis_group(grp(a))%head, grp(a), 'L', X(a), pg(a), ptree)
            X(a)%style = 2
            X(a)%level = -1
            X(a)%row_group = grp(a)
            X(a)%col_group = -1
            X(a)%headm = msh%basis_group(grp(a))%head
            X(a)%headn = 0
         enddo
         !>**** 2. pad the columns: leaf j of both halves = [stage-1 node (1,j); stage-1 node (2,j)] of P
         allocate (rsz(2, 2**Kh))
         rsz = 0
         do a = 1, 2
            if (.not. IOwnPgrp(ptree, pg(a))) cycle
            do k = 1, X(a)%ButterflyV%nblk_loc
               j = (k - 1)*X(a)%ButterflyV%inc + X(a)%ButterflyV%idx
               rsz(a, j) = size(X(a)%ButterflyV%blocks(k)%matrix, 1)
            enddo
         enddo
         if (ptree%pgrp(pgno)%nproc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, rsz, 2*2**Kh, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
         do a = 1, 2
            if (.not. IOwnPgrp(ptree, pg(a))) cycle
            do k = 1, X(a)%ButterflyV%nblk_loc
               j = (k - 1)*X(a)%ButterflyV%inc + X(a)%ButterflyV%idx
               r = rsz(a, j)
               n = rsz(1, j) + rsz(2, j)
               call BFD_matnew(X(a)%ButterflyV%blocks(k), n, r)
               do e = 1, r
                  if (a == 1) then
                     X(a)%ButterflyV%blocks(k)%matrix(e, e) = BPACK_cone
                  else
                     X(a)%ButterflyV%blocks(k)%matrix(rsz(1, j) + e, e) = BPACK_cone
                  endif
               enddo
            enddo
            call BFD_setlayout_nc(X(a), ptree)
         enddo
         deallocate (rsz)
         !>**** 3. both halves to the node's group
         do a = 1, 2
            call BFD_redist_nc(X(a), pgno, stats, ptree)
         enddo
         call BFD_tadd(4, t0)
         !>**** 4. Z = X1 - B12'*X2, Y1 = Z + C'*Z, Y2 = X2 - B21'*Y1
         call BFD_bmult_dist(B12, X(2), Tm, tolin, option, stats, ptree, msh, norecomp=.true.)
         call BFD_sumdist(X(1), Tm, -BPACK_cone, stats, ptree)
         call BF_delete(Tm, 1)
         call BFD_recompress(X(1), tolin, option, stats, ptree)
         call BFD_mul(Cs, X(1), Y1, tolin, option, stats, ptree, msh, X(1))
         call BF_delete(X(1), 1)
         call BF_copy_delete(Y1, X(1))
         call BFD_bmult_dist(B21, X(1), Tm, tolin, option, stats, ptree, msh, norecomp=.true.)
         call BFD_sumdist(X(2), Tm, -BPACK_cone, stats, ptree)
         call BF_delete(Tm, 1)
         call BFD_recompress(X(2), tolin, option, stats, ptree)
         !>**** 5. back to the child groups and into P: U and kernels 2..K by BF_copyback_partial
         t0 = MPI_Wtime()
         do a = 1, 2
            X(a)%row_group = grp(a)
            X(a)%col_group = -1
            call BF_ChangePattern(X(a), BFD_pattern(X(a)), 2, stats, ptree)
            call BFD_redist_nc(X(a), pg(a), stats, ptree)
            if (IOwnPgrp(ptree, pg(a))) call BF_copyback_partial(P, Kh, a, 'L', X(a), pg(a), ptree)
         enddo
         !>**** kernel 1 of P: block (a, jf) = V_a(j)^T * [W1(jf); W2(jf)], j = (jf+1)/2
         do a = 1, 2
            if (.not. IOwnPgrp(ptree, pg(a))) cycle
            do k = 1, X(a)%ButterflyV%nblk_loc
               j = (k - 1)*X(a)%ButterflyV%inc + X(a)%ButterflyV%idx
               do e = 1, 2
                  jf = 2*j - 2 + e
                  call GetBlockPID(ptree, pgno, 1, KP, 1, jf, 'C', pgno_sub)
                  dst = ptree%pgrp(pgno_sub)%head - ptree%pgrp(pgno)%head
                  call BFD_items_add(Lout, dst, a, jf, 0, 0, X(a)%ButterflyV%blocks(k)%matrix)
               enddo
            enddo
         enddo
         call BFD_exchange(Lout, Lin, ptree%pgrp(pgno)%Comm)
         nc = 0
         if (allocated(P%ButterflyKerl(1)%blocks)) nc = P%ButterflyKerl(1)%nc
         allocate (Vr(2, max(nc, 1)))
         do k = 1, Lin%n
            jj = (Lin%it(k)%hdr(2) - P%ButterflyKerl(1)%idx_c)/P%ButterflyKerl(1)%inc_c + 1
            call BFD_matset(Vr(Lin%it(k)%hdr(1), jj), Lin%it(k)%mat)
         enddo
         call BFD_items_free(Lin)
         do jj = 1, nc
            call BFD_vcat2(P%ButterflyKerl(1)%blocks(1, jj)%matrix, P%ButterflyKerl(1)%blocks(2, jj)%matrix, Wb)
            do a = 1, 2
               call BFD_mulnew('T', Vr(a, jj)%matrix, 'N', Wb%matrix, P%ButterflyKerl(1)%blocks(a, jj), stats)
               call BFD_matfree(Vr(a, jj))
            enddo
         enddo
         call BFD_matfree(Wb)
         deallocate (Vr)
         do a = 1, 2
            call BF_delete(X(a), 1)
         enddo
         call BF_get_rank(P, ptree)
         call BFD_tadd(4, t0)
      endif

      if (dbg) then
         err = BFD_factor_err(ho_bf1, level, ii, Pold, P, ptree, stats)
         if (ptree%MyID == ptree%pgrp(pgno)%head) write (*, '(A,I3,A,I5,A,I3,A,I4,A,ES10.3)') '   BFD Sblock level ', level, ' node ', ii, ' K ', KP, &
            ' nproc ', ptree%pgrp(pgno)%nproc, ' ||P_new x - F^-1 P_old x||/||F^-1 P_old x|| ', err
         call BF_delete(Pold, 1)
      endif
   end subroutine BFD_apply_factor

   !>**** bm = [A; B]
   subroutine BFD_vcat2(A, B, bm)
      implicit none
      DT::A(:, :), B(:, :)
      type(butterflymatrix)::bm
      call assert(size(A, 2) == size(B, 2), 'BFD_vcat2: column mismatch')
      call BFD_matnew(bm, size(A, 1) + size(B, 1), size(A, 2))
      if (size(A) > 0) bm%matrix(1:size(A, 1), :) = A
      if (size(B) > 0) bm%matrix(size(A, 1) + 1:, :) = B
   end subroutine BFD_vcat2

   !>**** ||Pnew x - F^{-1}(Pold x)|| / ||F^{-1}(Pold x)|| for random x (debugging of BFD_apply_factor)
   real(kind=8) function BFD_factor_err(ho_bf1, level, ii, Pold, Pnew, ptree, stats)
      implicit none
      type(hobf)::ho_bf1
      integer level, ii
      type(matrixblock)::Pold, Pnew
      type(proctree)::ptree
      type(Hstat)::stats
      DT, allocatable::x(:, :), y1(:, :), y2(:, :), z(:, :)
      real(kind=8), allocatable::xr(:, :)
      real(kind=8)::nrm(2)
      integer m, n, nvec, ierr
      nvec = 4
      m = Pold%M_loc
      n = Pold%N_loc
      allocate (x(max(n, 1), nvec), xr(max(n, 1), nvec), y1(max(m, 1), nvec), y2(max(m, 1), nvec), z(max(m, 1), nvec))
      call random_number(xr)
      x = xr - 0.5d0
      y1 = 0
      y2 = 0
      z = 0
      call BF_block_MVP_dat(Pold, 'N', m, n, nvec, x, max(n, 1), y1, max(m, 1), BPACK_cone, BPACK_czero, ptree, stats)
      call BF_block_MVP_inverse_dat(ho_bf1, level, ii, 'N', m, nvec, y1, max(m, 1), z, max(m, 1), ptree, stats)
      call BF_block_MVP_dat(Pnew, 'N', m, n, nvec, x, max(n, 1), y2, max(m, 1), BPACK_cone, BPACK_czero, ptree, stats)
      nrm = 0
      if (m > 0) then
         nrm(1) = sum(abs(y2(1:m, :) - z(1:m, :))**2)
         nrm(2) = sum(abs(z(1:m, :))**2)
      endif
      if (ptree%pgrp(Pold%pgno)%nproc > 1) call MPI_ALLREDUCE(MPI_IN_PLACE, nrm, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(Pold%pgno)%Comm, ierr)
      BFD_factor_err = sqrt(nrm(1)/max(nrm(2), BPACK_SafeUnderflow))
      deallocate (x, xr, y1, y2, z)
   end function BFD_factor_err

   !>**** node ii (at level, below the node rowblock of BFD_Sblock) of the Sblock update: extract the node's row part of
   !> block_o, apply the node's factor inverse to it (BFD_apply_factor) and copy it back
   subroutine BFD_Sblock_node(ho_bf1, level, rowblock, ii, N_diag, block_o, option, stats, ptree, msh)
      implicit none
      type(hobf)::ho_bf1
      integer level, rowblock, ii, N_diag
      type(matrixblock)::block_o
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock), pointer::blocks
      type(matrixblock)::agent
      integer ii_loc, level_butterfly_loc
      ii_loc = ii - ((rowblock - 1)*N_diag + 1) + 1
      level_butterfly_loc = ho_bf1%Maxlevel + 1 - level
      blocks => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
      call assert(IOwnPgrp(ptree, blocks%pgno), 'BFD_Sblock: I do not own this pgno')
      call BF_extract_partial(block_o, level_butterfly_loc, ii_loc, blocks%headm, blocks%row_group, 'L', agent, blocks%pgno, ptree)
      call BFD_apply_factor(ho_bf1, level, ii, agent, option, stats, ptree, msh)
      call BF_copyback_partial(block_o, level_butterfly_loc, ii_loc, 'L', agent, blocks%pgno, ptree)
      call BF_delete(agent, 1)
   end subroutine BFD_Sblock_node

   !>**** HODBF with bf_algebra=1: the forward block Z_ij^l of rowblock at level_c is multiplied by the inverse of its
   !> diagonal block, D^{-1} Z (replaces Bplus_Sblock_randomized_memfree): dense leaf inverses on U, then from the finest
   !> node level to the coarsest the partial of every node is updated by BFD_apply_factor; one recompression at the end.
   subroutine BFD_Sblock(ho_bf1, level_c, rowblock, option, stats, ptree, msh)
      implicit none
      type(hobf)::ho_bf1
      integer level_c, rowblock
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock), pointer::block_o, blocks
      integer K, level, N_diag, idx_start_diag, idx_end_diag, ii, index_i, index_i_loc_k, num_vect_sub, pat0
      integer nth, nnode, nout, nin, lvl
      logical conc
      real(kind=8)::t0, t1

      t0 = MPI_Wtime()
      call Bplus_copy(ho_bf1%levels(level_c)%BP(rowblock), ho_bf1%levels(level_c)%BP_inverse_update(rowblock))
      block_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)
      K = block_o%level_butterfly
      if (K == 0) then
         call BFD_flops_flush(stats)
         call LR_Sblock(ho_bf1, level_c, rowblock, ptree, stats)
         stats%Flop_Tmp = 0 ! LR_Sblock already added its flops to stats%Flop_Factor
         return
      endif
      pat0 = BFD_pattern(block_o)
      call BF_ChangePattern(block_o, pat0, 2, stats, ptree)
      do level = ho_bf1%Maxlevel + 1, level_c + 1, -1
         N_diag = 2**(level - level_c - 1)
         idx_start_diag = max((rowblock - 1)*N_diag + 1, ho_bf1%levels(level)%Bidxs)
         idx_end_diag = min(rowblock*N_diag, ho_bf1%levels(level)%Bidxe)
         if (level == ho_bf1%Maxlevel + 1) then
            t1 = MPI_Wtime()
            do ii = idx_start_diag, idx_end_diag
               blocks => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
               index_i = ii - ((rowblock - 1)*N_diag + 1) + 1
               index_i_loc_k = (index_i - block_o%ButterflyU%idx)/block_o%ButterflyU%inc + 1
               num_vect_sub = size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix, 2)
               call Full_block_MVP_dat(blocks, 'N', blocks%M, num_vect_sub, block_o%ButterflyU%blocks(index_i_loc_k)%matrix, &
                  size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix, 1), block_o%ButterflyU%blocks(index_i_loc_k)%matrix, &
                  size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix, 1), BPACK_cone, BPACK_czero)
            enddo
            call BFD_tadd(2, t1)
         else
            !>**** the nodes of this level update disjoint row parts of block_o: when every node is on one process
            !> they run as concurrent OpenMP tasks; with fewer nodes than threads each task gets nth/nodes threads for
            !> its own parallel loops (one more active nesting level)
            conc = .false.
#ifdef HAVE_OPENMP
            nth = omp_get_max_threads()
            nnode = idx_end_diag - idx_start_diag + 1
            conc = nth > 1 .and. nnode > 1 .and. .not. omp_in_parallel()
            do ii = idx_start_diag, idx_end_diag
               if (ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%pgno)%nproc /= 1) conc = .false.
            enddo
#endif
            if (conc) then
#ifdef HAVE_OPENMP
               nout = min(nnode, nth)
               nin = max(1, nth/nout)
               lvl = omp_get_max_active_levels()
               if (nin > 1) call omp_set_max_active_levels(max(lvl, 2))
               !$omp parallel num_threads(nout) default(shared) private(ii)
               !$omp single
               do ii = idx_start_diag, idx_end_diag
                  !$omp task default(shared) firstprivate(ii)
                  call omp_set_num_threads(nin)
                  call BFD_Sblock_node(ho_bf1, level, rowblock, ii, N_diag, block_o, option, stats, ptree, msh)
                  !$omp end task
               enddo
               !$omp end single
               !$omp end parallel
               call omp_set_max_active_levels(lvl)
#endif
            else
               do ii = idx_start_diag, idx_end_diag
                  call BFD_Sblock_node(ho_bf1, level, rowblock, ii, N_diag, block_o, option, stats, ptree, msh)
               enddo
            endif
         endif
      enddo
      call BF_get_rank(block_o, ptree)
      BFD_rankgrow(1) = max(BFD_rankgrow(1), block_o%rankmax)
      call BFD_recompress(block_o, option%tol_rand, option, stats, ptree)
      BFD_rankgrow(2) = max(BFD_rankgrow(2), block_o%rankmax)
      call BFD_tadd(1, t0)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) write (*, '(A10,I5,A6,I3,A8,I3)') 'OneL No. ', rowblock, ' rank:', block_o%rankmax, ' L_butt:', block_o%level_butterfly
   end subroutine BFD_Sblock

end module Bplus_deterministic
