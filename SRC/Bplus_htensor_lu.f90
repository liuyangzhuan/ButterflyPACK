! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National
! Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual
! Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government
! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting
! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to
! the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).
!> @file Bplus_htensor_lu.f90
!> @brief LU factorization, inverse and solution of a forest of H-TENSOR matrices with Tucker (tensor low-rank) admissible blocks
!> @details A P x P forest of HTENSOR matrices (block (p,q) between trees p and q, all trees of the same structure) is copied
!> into recursive block trees. A node is a dense leaf (style 1), a Tucker block (style 2) or a hierarchical block (style 4)
!> whose sons are the 2^Ndim x 2^Ndim pairs of the children of its row and column clusters. Tensors are in natural order
!> (dimension 1 fastest), and a Tucker block is B = (U_Ndim x ... x U_1) C (V_Ndim x ... x V_1)^T.
!> The factorization is a block LU with whole-cluster pivots, over the trees and recursively inside them (htlu_BLU): for
!> the diagonal blocks i = 1..K of a block, factor the current diagonal block S_ii, set the lower blocks X_ji = A_ji S_ii^{-1}
!> (XSM) and update the trailing blocks A_jl = A_jl - X_ji A_il (Schur updates C = C - A*B). The upper blocks keep their
!> Schur-updated values, and a solve goes through the block factors (htlu_bsolve). Pivoting is only inside the dense
!> diagonal leaves. Tucker arithmetic, truncated at option%tol_rand:
!>   - a solve with a Tucker block applies the solve to its column Kronecker factor and recompresses it by a truncated
!>     HOSVD of the result;
!>   - a Schur update of a Tucker block by a product of two Tucker blocks is formed exactly in Tucker form and added by a
!>     QR of the stacked factors; other updates of a Tucker block rebuild it from products: per-mode bases from a few
!>     random products of the updated block (every product gives all fibers of its output) and the core by an exact
!>     projection;
!>   - a product with a Tucker factor that updates a hierarchical block is formed in Tucker form and added block by block.
!> The inverse (HTLU_dist_inverse) is an in-place block Gauss-Jordan sweep on the same trees with the same kernels.
!> Singular values are dropped below option%tol_rand relative to each block and below option%tol_rand*normest, with
!> normest an estimate of ||A||_2 (LU) or ||A^{-1}||_2 (inverse) given by the caller.
!> Distributed memory (htlu_dist): the blocks between the clusters of a level Ld of all trees (K = P*2^(Ndim*Ld) clusters, at
!> least the level at which every cluster lies on one process in the HTENSOR layout) are in a 2D block-cyclic layout over the MPI
!> processes. Every process runs the same recursion over the coarse levels (the trees, then the children inside them); each of its
!> steps is a set of independent operations on level-Ld blocks (factor or invert a diagonal block, a solve or a product with a
!> diagonal block, Schur updates), done by the owner of the updated block after the remote operands are fetched. This is the
!> factorization of the one-process code with the coarse blocks split at level Ld. OpenMP tasks run the operations of a step.
!> Vectors of the solves and products are full on every process (trees one after another, natural order in every tree).
#include "ButterflyPACK_config.fi"
module Bplus_htensor_lu
   use BPACK_DEFS
!$ use omp_lib
   use MISC_Utilities
   use MISC_DenseLA
   implicit none

   !>**** a matrix (factor of a Tucker block, or a child piece of a vector set)
   type htlu_fmat
      DT, allocatable :: mat(:, :)
   end type htlu_fmat

   !>**** a block of the recursive tree
   type htlu_node
      integer :: style = 0 !< 1: dense, 2: Tucker, 4: hierarchical
      integer :: level = 0 !< tree level of the row and column clusters
      integer :: Ndim = 0
      integer :: nch = 0 !< children per cluster (2^Ndim) of a hierarchical block
      integer, allocatable :: m(:), n(:) !< row and column sizes per dimension
      integer, allocatable :: rg(:), cg(:) !< row and column groups per dimension
      type(htlu_node), pointer :: sons(:, :) => null() !< sons(i,j): row child i, column child j; child i is c_d = mod((i-1)/2^(d-1),2)+1
      DT, allocatable :: D(:, :) !< dense block, prod(m) x prod(n); for a diagonal leaf its LU after the factorization
      integer, allocatable :: ipiv(:) !< pivots of the LU of a dense diagonal leaf
      type(htlu_fmat), allocatable :: U(:), V(:) !< Tucker factors per dimension, m(d) x r(d) and n(d) x rc(d)
      DT, allocatable :: C(:, :) !< Tucker core, prod(r) x prod(rc)
   end type htlu_node

   !>**** one term coef*a*b (b absent: coef*a) of an operator whose Tucker form is rebuilt
   type htlu_term
      DT :: coef = 1
      type(htlu_node), pointer :: a => null(), b => null()
   end type htlu_term

   !>**** the HTENSOR blocks of one level, keyed by their row and column clusters
   type htlu_index_level
      integer :: nb = 0
      integer, allocatable :: head(:), next(:), cidx(:), ll(:), bb(:)
   end type htlu_index_level

   !>**** a pointer to a block
   type htlu_nptr
      type(htlu_node), pointer :: p => null()
   end type htlu_nptr

   !>**** the forest distributed over the MPI processes: the blocks between the clusters of level Ld of all trees (K = P*nc
   !>**** clusters, numbered tree by tree and in the recursive child order inside a tree) in a 2D block-cyclic layout; block (a,b)
   !>**** is on the process at row mod(a-1,nprow) and column mod(b-1+skew*((b-1)/npcol),npcol) of the process grid (rank =
   !>**** row*npcol + column): the skew spreads the diagonal blocks over all processes
   type htlu_dist
      integer :: Ndim = 0, P = 0, Ld = 0, nch = 0, nc = 0, K = 0, nt = 0 !< nc clusters per tree, nt entries per tree
      integer :: comm = 0, nproc = 1, myid = 0, nprow = 1, npcol = 1, myrow = 0, mycol = 0, skew = 0
      integer :: nr = 0, nl = 0 !< numbers of local block rows and columns
      integer, allocatable :: Ns(:) !< grid size per dimension of every tree
      integer, allocatable :: coff(:, :), csz(:, :), cgrp(:, :) !< offset, size and group per dimension of cluster m (1..nc) of a tree
      integer, allocatable :: boxes(:, :) !< boxes(1:2*Ndim, r): offsets (from 0) and sizes per dimension of the vector box of rank r in the HTENSOR layout
      type(htlu_node), allocatable :: blk(:, :) !< blk(i, j): block (a, b) with (a-1)/nprow = i-1 and (b-1)/npcol = j-1
      integer, allocatable :: rmap(:, :) !< scratch: index of a fetched block by its position
      integer :: inverse = 0 !< 1: blk holds the inverse
      real(kind=8) :: tcomm = 0, tcomp = 0, vrecv = 0 !< seconds in the exchanges and in the block operations, MB received
      integer :: nstep = 0 !< number of exchanges
      integer :: rr = 0 !< offset of the round-robin assignment of spread steps
      !>**** profile of the steps by kind (1 factor, 2 invert a diagonal block, 3 solve, 4 right and 5 left product with a
      !>**** diagonal block, 6 Schur updates, 7 products into new blocks) and by the number c of summed clusters (class
      !>**** ceiling(log_{2^Ndim} c)): wall seconds, block-operation seconds of this process, steps, local targets
      real(kind=8) :: pt0 = 0, pc0 = 0
      real(kind=8) :: pwall(7, 0:9) = 0, pwork(7, 0:9) = 0
      integer :: pcnt(7, 0:9) = 0, pntg(7, 0:9) = 0
   end type htlu_dist

   !>**** blocks above this size (rows) create OpenMP tasks for their children
   integer, parameter :: htlu_task_rows = 1024
   !>**** products with at most this many multiply-adds use matmul instead of BLAS (whose per-call overhead dominates for the
   !>**** small Tucker factors and cores, and whose buffer allocation contends between threads)
   real(kind=8), parameter :: htlu_small_gemm = 32768d0

   !>**** thread-seconds and counts of the leaf operations of the factorization (diagnostics): 1 Tucker rebuilds of Schur
   !>**** updates, 2 block-LU solves of Tucker blocks (XSM), 3 Tucker products, 4 Tucker multiplications of the inverse, 5 dense leaf
   !>**** operations, 6 Tucker rebuilds of sums of Tucker blocks, 7 deterministic Tucker sums (QR of the concatenated factors)
   real(kind=8) :: htlu_tsec(7) = 0d0
   integer(kind=8) :: htlu_tcnt(7) = 0
   !>**** deterministic Tucker sums whose new core would exceed this many entries fall back to the randomized rebuild
   integer(kind=8), parameter :: htlu_sum_cap = 4000000
   !>**** rebuilds by level of the target and kind (1: Schur product terms, 2: sums of Tucker blocks)
   real(kind=8) :: htlu_rlev_t(0:15, 2) = 0d0
   integer(kind=8) :: htlu_rlev_n(0:15, 2) = 0

contains

   !>**** add thread-seconds dt to category c
   subroutine htlu_tadd(c, dt)
      integer c
      real(kind=8) dt
      !$omp atomic
      htlu_tsec(c) = htlu_tsec(c) + dt
      !$omp atomic
      htlu_tcnt(c) = htlu_tcnt(c) + 1
   end subroutine htlu_tadd

   !>**** print and reset the timers of the leaf operations; wall: wall time of the factorization, nth: threads
   subroutine HTLU_print_timers(wall, nth)
      real(kind=8) wall
      integer nth, c, l
      character(len=24), parameter :: cname(7) = (/'Tucker rebuilds (Schur) ', 'Tucker solves           ', &
         'Tucker products         ', 'Tucker inverse products ', 'dense leaf operations   ', 'Tucker rebuilds (sums)  ', &
         'Tucker sums (QR)        '/)
      write (*, '(A,F10.1,A,I4,A)') ' tensor LU leaf operations (thread-seconds; wall', wall, ' s x', nth, ' threads):'
      do c = 1, 7
         if (htlu_tcnt(c) > 0) write (*, '(4X,A24,F12.1,A,I12,A)') cname(c), htlu_tsec(c), ' s in', htlu_tcnt(c), ' calls'
      enddo
      write (*, '(4X,A24,F12.1,A)') 'sum / (wall x threads)', 100d0*sum(htlu_tsec)/max(1d-30, wall*nth), ' %'
      write (*, '(A)') '    rebuilds by level: level, Schur (thread-s, calls), sums (thread-s, calls)'
      do l = 0, 15
         if (htlu_rlev_n(l, 1) + htlu_rlev_n(l, 2) > 0) write (*, '(I10,2(F12.1,I12))') l, htlu_rlev_t(l, 1), htlu_rlev_n(l, 1), &
            htlu_rlev_t(l, 2), htlu_rlev_n(l, 2)
      enddo
      htlu_tsec = 0
      htlu_tcnt = 0
      htlu_rlev_t = 0
      htlu_rlev_n = 0
   end subroutine HTLU_print_timers

   !>**** complex conjugate (identity for real data)
   elemental function htlu_conj(x) result(y)
      DT, intent(in) :: x
      DT :: y
#if DAT==0 || DAT==2
      y = conjg(x)
#else
      y = x
#endif
   end function htlu_conj

   !>**** uniform random entries in [-0.5,0.5] (real and imaginary parts)
   subroutine htlu_random(A)
      DT :: A(:, :)
      real(kind=8), allocatable :: r1(:, :), r2(:, :)
      allocate (r1(size(A, 1), size(A, 2)), r2(size(A, 1), size(A, 2)))
      call random_number(r1)
      call random_number(r2)
#if DAT==0 || DAT==2
      A = cmplx(r1 - 0.5d0, r2 - 0.5d0, kind=kind(A))
#else
      A = r1 - 0.5d0
#endif
      deallocate (r1, r2)
   end subroutine htlu_random

   !>**** C = op(A)*op(B), C allocated here
   subroutine htlu_gemm(ta, tb, A, B, C)
      character ta, tb
      DT :: A(:, :), B(:, :)
      DT, allocatable :: C(:, :)
      integer m, n, k
      if (ta == 'N') then
         m = size(A, 1); k = size(A, 2)
      else
         m = size(A, 2); k = size(A, 1)
      endif
      if (tb == 'N') then
         n = size(B, 2)
      else
         n = size(B, 1)
      endif
      if (allocated(C)) deallocate (C)
      allocate (C(m, n))
      if (m*n == 0) return
      if (k == 0) then
         C = 0
         return
      endif
      if (dble(m)*dble(n)*dble(k) <= htlu_small_gemm) then
         if (ta == 'N' .and. tb == 'N') then
            C = matmul(A, B)
         elseif (ta == 'N') then
            if (tb == 'T') then
               C = matmul(A, transpose(B))
            else
               C = matmul(A, htlu_conj(transpose(B)))
            endif
         elseif (tb == 'N') then
            if (ta == 'T') then
               C = matmul(transpose(A), B)
            else
               C = matmul(htlu_conj(transpose(A)), B)
            endif
         else
            call gemmf90(A, size(A, 1), B, size(B, 1), C, m, ta, tb, m, n, k, BPACK_cone, BPACK_czero)
         endif
         return
      endif
      call gemmf90(A, size(A, 1), B, size(B, 1), C, m, ta, tb, m, n, k, BPACK_cone, BPACK_czero)
   end subroutine htlu_gemm

   !>**** child multi-index c(1:Ndim) in {1,2} of child number i
   subroutine htlu_child(Ndim, i, c)
      integer Ndim, i, c(Ndim), d
      do d = 1, Ndim
         c(d) = mod((i - 1)/2**(d - 1), 2) + 1
      enddo
   end subroutine htlu_child

   !>**** offsets and sizes per dimension of row child i (side 'r') or column child i (side 'c') of a hierarchical block
   subroutine htlu_child_box(A, side, i, offs, sz)
      type(htlu_node) :: A
      character side
      integer i, offs(:), sz(:), c(A%Ndim), d
      call htlu_child(A%Ndim, i, c)
      do d = 1, A%Ndim
         if (side == 'r') then
            sz(d) = A%sons(i, 1)%m(d)
            offs(d) = 0
            if (c(d) == 2) offs(d) = A%sons(1, 1)%m(d)
         else
            sz(d) = A%sons(1, i)%n(d)
            offs(d) = 0
            if (c(d) == 2) offs(d) = A%sons(1, 1)%n(d)
         endif
      enddo
   end subroutine htlu_child_box

   !>**** box of a natural-order tensor with columns. mode 0: B = P(box), 1: P(box) = B, 2: P(box) = P(box) + B
   subroutine htlu_box(Ndim, dP, offs, sz, P, B, mode)
      integer Ndim, dP(Ndim), offs(Ndim), sz(Ndim), mode
      DT :: P(:, :), B(:, :)
      integer k, t, d, id, col, nrun, l, ip, stride
      nrun = 1
      do d = 2, Ndim
         nrun = nrun*sz(d)
      enddo
      do col = 1, size(B, 2)
         do k = 0, nrun - 1
            ip = offs(1) + 1
            stride = dP(1)
            t = k
            do d = 2, Ndim
               id = mod(t, sz(d))
               t = t/sz(d)
               ip = ip + (offs(d) + id)*stride
               stride = stride*dP(d)
            enddo
            l = k*sz(1)
            select case (mode)
            case (0)
               B(l + 1:l + sz(1), col) = P(ip:ip + sz(1) - 1, col)
            case (1)
               P(ip:ip + sz(1) - 1, col) = B(l + 1:l + sz(1), col)
            case (2)
               P(ip:ip + sz(1) - 1, col) = P(ip:ip + sz(1) - 1, col) + B(l + 1:l + sz(1), col)
            end select
         enddo
      enddo
   end subroutine htlu_box

   !>**** the pieces Xc(i) of X (rows on the row clusters, side 'r', or the column clusters, side 'c', of block A) on the children
   subroutine htlu_gather(A, side, X, Xc)
      type(htlu_node) :: A
      character side
      DT :: X(:, :)
      type(htlu_fmat) :: Xc(:)
      integer i, offs(A%Ndim), sz(A%Ndim), dP(A%Ndim)
      if (side == 'r') then
         dP = A%m
      else
         dP = A%n
      endif
      do i = 1, A%nch
         call htlu_child_box(A, side, i, offs, sz)
         allocate (Xc(i)%mat(product(sz), size(X, 2)))
         call htlu_box(A%Ndim, dP, offs, sz, X, Xc(i)%mat, 0)
      enddo
   end subroutine htlu_gather

   !>**** the inverse of htlu_gather; mode 1: copy, 2: add. The pieces are freed
   subroutine htlu_scatter(A, side, Xc, X, mode)
      type(htlu_node) :: A
      character side
      DT :: X(:, :)
      type(htlu_fmat) :: Xc(:)
      integer mode, i, offs(A%Ndim), sz(A%Ndim), dP(A%Ndim)
      if (side == 'r') then
         dP = A%m
      else
         dP = A%n
      endif
      do i = 1, A%nch
         call htlu_child_box(A, side, i, offs, sz)
         call htlu_box(A%Ndim, dP, offs, sz, X, Xc(i)%mat, mode)
         deallocate (Xc(i)%mat)
      enddo
   end subroutine htlu_scatter

   !>**** K = mats(Ndim) x ... x mats(1) (Kronecker product, dimension 1 fastest), conjugated with conjflag = 1
   subroutine htlu_kron(Ndim, mats, K, conjflag)
      integer Ndim, conjflag, d, id, jd, pm, pr, md, rd
      type(htlu_fmat) :: mats(:)
      DT, allocatable :: K(:, :), K0(:, :)
      DT f
      if (allocated(K)) deallocate (K)
      allocate (K(size(mats(1)%mat, 1), size(mats(1)%mat, 2)))
      if (conjflag == 1) then
         K = htlu_conj(mats(1)%mat)
      else
         K = mats(1)%mat
      endif
      do d = 2, Ndim
         pm = size(K, 1)
         pr = size(K, 2)
         md = size(mats(d)%mat, 1)
         rd = size(mats(d)%mat, 2)
         call move_alloc(K, K0)
         allocate (K(pm*md, pr*rd))
         do jd = 1, rd
            do id = 1, md
               f = mats(d)%mat(id, jd)
               if (conjflag == 1) f = htlu_conj(f)
               K((id - 1)*pm + 1:id*pm, (jd - 1)*pr + 1:jd*pr) = f*K0
            enddo
         enddo
         deallocate (K0)
      enddo
   end subroutine htlu_kron

   !>**** one mode product: B = op(M) applied along dimension d of the natural-order tensor A (dims, R columns); op(M) is sout x dims(d)
   subroutine htlu_mode1(Ndim, dims, d, R, A, M, op, sout, B)
      integer Ndim, dims(Ndim), d, R, sout, pre, sin, k, suf, i
      DT :: A(*), B(*), M(:, :)
      character op
      DT, allocatable :: At(:, :), Bt(:, :)
      pre = 1
      do k = 1, d - 1
         pre = pre*dims(k)
      enddo
      sin = dims(d)
      suf = R
      do k = d + 1, Ndim
         suf = suf*dims(k)
      enddo
      if (sout*suf == 0) return
      if (sin == 0) then
         B(1:pre*sout*suf) = 0
         return
      endif
      if (pre == 1) then
         if (dble(sout)*dble(suf)*dble(sin) <= htlu_small_gemm) then
            allocate (At(sin, suf))
            At = reshape(A(1:sin*suf), (/sin, suf/))
            call htlu_gemm(op, 'N', M, At, Bt)
            B(1:sout*suf) = reshape(Bt, (/sout*suf/))
            deallocate (At, Bt)
         else
            call gemmf90(M, size(M, 1), A, sin, B, sout, op, 'N', sout, suf, sin, BPACK_cone, BPACK_czero)
         endif
      else
         ! slice k: B_k (pre x sout) = A_k (pre x sin) * op(M)^T, i.e. B_k^T = op(M) A_k^T: all slices in one product
         allocate (At(sin, pre*suf))
         do k = 1, suf
            do i = 1, sin
               At(i, (k - 1)*pre + 1:k*pre) = A((k - 1)*pre*sin + (i - 1)*pre + 1:(k - 1)*pre*sin + i*pre)
            enddo
         enddo
         call htlu_gemm(op, 'N', M, At, Bt)
         do k = 1, suf
            do i = 1, sout
               B((k - 1)*pre*sout + (i - 1)*pre + 1:(k - 1)*pre*sout + i*pre) = Bt(i, (k - 1)*pre + 1:k*pre)
            enddo
         enddo
         deallocate (At, Bt)
      endif
   end subroutine htlu_mode1

   !>**** Y = (op(mats(Ndim)) x ... x op(mats(1))) X for the natural-order tensor X of sizes dims_in with columns; Y allocated here
   subroutine htlu_modes(Ndim, dims_in, X, mats, op, Y)
      integer Ndim, dims_in(Ndim), dims(Ndim), d, R, sout
      DT :: X(:, :)
      type(htlu_fmat) :: mats(:)
      character op
      DT, allocatable :: Y(:, :), w1(:), w2(:)
      R = size(X, 2)
      dims = dims_in
      allocate (w1(max(1, size(X))))
      if (size(X) > 0) w1(1:size(X)) = reshape(X, (/size(X)/))
      do d = 1, Ndim
         if (op == 'N') then
            sout = size(mats(d)%mat, 1)
         else
            sout = size(mats(d)%mat, 2)
         endif
         allocate (w2(max(1, (product(dims)/max(1, dims(d)))*sout*R)))
         call htlu_mode1(Ndim, dims, d, R, w1, mats(d)%mat, op, sout, w2)
         dims(d) = sout
         call move_alloc(w2, w1)
      enddo
      if (allocated(Y)) deallocate (Y)
      allocate (Y(product(dims), R))
      if (size(Y) > 0) Y = reshape(w1(1:size(Y)), (/product(dims), R/))
      deallocate (w1)
   end subroutine htlu_modes

   !>**** orthonormal bases Q(d) of the mode-d unfoldings of the natural-order tensor Y (dims, columns): singular values below
   !>**** tol(1) times the largest one of the unfolding or below tol(2) are dropped (eigenvalues of the Gram matrix)
   subroutine htlu_mode_bases(Ndim, dims, Y, tol, Q)
      integer Ndim, dims(Ndim), d, md, pre, suf, k, rank, R, i
      DT :: Y(:, :)
      real(kind=8) tol(2)
      type(htlu_fmat) :: Q(:)
      DT, allocatable :: G(:, :), UU(:, :), VV(:, :), w(:), Wd(:, :)
      DTR, allocatable :: S(:)
      R = size(Y, 2)
      allocate (w(max(1, size(Y))))
      if (size(Y) > 0) w(1:size(Y)) = reshape(Y, (/size(Y)/))
      do d = 1, Ndim
         md = dims(d)
         pre = 1
         do k = 1, d - 1
            pre = pre*dims(k)
         enddo
         suf = R
         do k = d + 1, Ndim
            suf = suf*dims(k)
         enddo
         allocate (UU(md, md), VV(md, md), S(md))
         ! the mode-d unfolding Wd (md x pre*suf), Wd(i, p + pre*(k-1)) = Y(p, i, k), and its Gram matrix Wd Wd^H
         allocate (Wd(md, pre*suf))
         do k = 1, suf
            do i = 1, md
               Wd(i, (k - 1)*pre + 1:k*pre) = w((k - 1)*pre*md + (i - 1)*pre + 1:(k - 1)*pre*md + i*pre)
            enddo
         enddo
         call htlu_gemm('N', 'C', Wd, Wd, G)
         deallocate (Wd)
         call SVD_Truncate(G, md, md, md, UU, VV, S, tol(1)**2, max(tol(2)**2, BPACK_SafeUnderflow), rank)
         if (allocated(Q(d)%mat)) deallocate (Q(d)%mat)
         allocate (Q(d)%mat(md, rank))
         if (S(1) < BPACK_SafeUnderflow) then
            Q(d)%mat = 0
            Q(d)%mat(1, 1) = 1
         else
            Q(d)%mat = UU(:, 1:rank)
         endif
         deallocate (G, UU, VV, S)
      enddo
      deallocate (w)
   end subroutine htlu_mode_bases

   !>**** W ~ kron(Q) G: bases Q(d) of the modes of W (dims, columns) and G = kron(Q)^H W
   subroutine htlu_tuckerize(Ndim, dims, W, tol, Q, G)
      integer Ndim, dims(Ndim)
      DT :: W(:, :)
      real(kind=8) tol(2)
      type(htlu_fmat) :: Q(:)
      DT, allocatable :: G(:, :)
      call htlu_mode_bases(Ndim, dims, W, tol, Q)
      call htlu_modes(Ndim, dims, W, Q, 'C', G)
   end subroutine htlu_tuckerize

   !>**** truncated HOSVD of the core of the Tucker block A, whose factors have orthonormal columns: every mode of the core
   !>**** (row and column dimensions) is truncated at tol, and the factors are rotated accordingly
   subroutine htlu_core_trim(A, tol)
      type(htlu_node) :: A
      real(kind=8) tol(2)
      integer Nd, d, pr, pc, dims(2*A%Ndim)
      type(htlu_fmat), allocatable :: Q(:)
      DT, allocatable :: G(:, :), core(:, :), T(:, :)
      Nd = A%Ndim
      do d = 1, Nd
         dims(d) = size(A%U(d)%mat, 2)
         dims(Nd + d) = size(A%V(d)%mat, 2)
      enddo
      allocate (core(size(A%C), 1), Q(2*Nd))
      core(:, 1) = reshape(A%C, (/size(A%C)/))
      call htlu_tuckerize(2*Nd, dims, core, tol, Q, G)
      pr = 1
      pc = 1
      do d = 1, Nd
         call htlu_gemm('N', 'N', A%U(d)%mat, Q(d)%mat, T)
         call move_alloc(T, A%U(d)%mat)
         call htlu_gemm('N', 'N', A%V(d)%mat, Q(Nd + d)%mat, T)
         call move_alloc(T, A%V(d)%mat)
         pr = pr*size(Q(d)%mat, 2)
         pc = pc*size(Q(Nd + d)%mat, 2)
      enddo
      deallocate (A%C)
      allocate (A%C(pr, pc))
      A%C = reshape(G(:, 1), (/pr, pc/))
      do d = 1, 2*Nd
         deallocate (Q(d)%mat)
      enddo
      deallocate (Q, G, core)
   end subroutine htlu_core_trim

   !>**** thin QR per dimension, mats(d) = Q(d) R(d), and Rk = R(Ndim) x ... x R(1)
   subroutine htlu_orth(Ndim, mats, Q, Rk)
      integer Ndim, d, md, rd, kd, i, j
      type(htlu_fmat) :: mats(:), Q(:)
      DT, allocatable :: Rk(:, :), A(:, :), tau(:)
      type(htlu_fmat), allocatable :: R(:)
      allocate (R(Ndim))
      do d = 1, Ndim
         md = size(mats(d)%mat, 1)
         rd = size(mats(d)%mat, 2)
         kd = min(md, rd)
         allocate (A(md, rd), tau(kd))
         A = mats(d)%mat
         call geqrff90(A, tau)
         allocate (R(d)%mat(kd, rd))
         R(d)%mat = 0
         do j = 1, rd
            do i = 1, min(j, kd)
               R(d)%mat(i, j) = A(i, j)
            enddo
         enddo
         call un_or_gqrf90(A, tau, md, kd, kd)
         if (allocated(Q(d)%mat)) deallocate (Q(d)%mat)
         allocate (Q(d)%mat(md, kd))
         Q(d)%mat = A(:, 1:kd)
         deallocate (A, tau)
      enddo
      call htlu_kron(Ndim, R, Rk, 0)
      do d = 1, Ndim
         deallocate (R(d)%mat)
      enddo
      deallocate (R)
   end subroutine htlu_orth

   !>**** ranks of a Tucker block: r(d) (row side) and rc(d) (column side)
   subroutine htlu_ranks(A, r, rc)
      type(htlu_node) :: A
      integer r(:), rc(:), d
      do d = 1, A%Ndim
         r(d) = size(A%U(d)%mat, 2)
         rc(d) = size(A%V(d)%mat, 2)
      enddo
   end subroutine htlu_ranks

   !>**** Y = Y + alpha*op(A)*X for a Tucker block A
   subroutine htlu_tucker_apply(A, op, X, Y, alpha)
      type(htlu_node) :: A
      character op
      DT :: X(:, :), Y(:, :), alpha
      DT, allocatable :: T1(:, :), T2(:, :), T3(:, :)
      integer r(A%Ndim), rc(A%Ndim)
      call htlu_ranks(A, r, rc)
      if (op == 'N') then
         call htlu_modes(A%Ndim, A%n, X, A%V, 'T', T1)
         call htlu_gemm('N', 'N', A%C, T1, T2)
         call htlu_modes(A%Ndim, r, T2, A%U, 'N', T3)
      else
         call htlu_modes(A%Ndim, A%m, X, A%U, 'T', T1)
         call htlu_gemm('T', 'N', A%C, T1, T2)
         call htlu_modes(A%Ndim, rc, T2, A%V, 'N', T3)
      endif
      Y = Y + alpha*T3
      deallocate (T1, T2, T3)
   end subroutine htlu_tucker_apply

   !>**** the dense matrix of a leaf-level block (dense or Tucker)
   subroutine htlu_dense_of(A, D)
      type(htlu_node) :: A
      DT, allocatable :: D(:, :), KU(:, :), KV(:, :), T(:, :)
      if (allocated(D)) deallocate (D)
      if (A%style == 1) then
         allocate (D(size(A%D, 1), size(A%D, 2)))
         D = A%D
      else
         call assert(A%style == 2, 'htlu_dense_of: a hierarchical block')
         call htlu_kron(A%Ndim, A%U, KU, 0)
         call htlu_kron(A%Ndim, A%V, KV, 0)
         call htlu_gemm('N', 'T', A%C, KV, T)
         call htlu_gemm('N', 'N', KU, T, D)
         deallocate (KU, KV, T)
      endif
   end subroutine htlu_dense_of

   !>**** Y = Y + alpha*op(A)*X for a block A of any style; X and Y hold natural-order tensors on the column (row) clusters of A
   recursive subroutine htlu_mult(A, op, X, Y, alpha)
      type(htlu_node) :: A
      character op
      DT :: X(:, :), Y(:, :), alpha
      type(htlu_fmat), allocatable :: Xc(:), Yc(:)
      integer i, j, R, mm, nn
      R = size(X, 2)
      if (R == 0) return
      select case (A%style)
      case (1)
         mm = size(A%D, 1)
         nn = size(A%D, 2)
         if (op == 'N') then
            call gemmf90(A%D, mm, X, nn, Y, mm, 'N', 'N', mm, R, nn, alpha, BPACK_cone)
         else
            call gemmf90(A%D, mm, X, mm, Y, nn, 'T', 'N', nn, R, mm, alpha, BPACK_cone)
         endif
      case (2)
         call htlu_tucker_apply(A, op, X, Y, alpha)
      case (4)
         allocate (Xc(A%nch), Yc(A%nch))
         if (op == 'N') then
            call htlu_gather(A, 'c', X, Xc)
            do i = 1, A%nch
               allocate (Yc(i)%mat(product(A%sons(i, 1)%m), R))
               Yc(i)%mat = 0
            enddo
         else
            call htlu_gather(A, 'r', X, Xc)
            do i = 1, A%nch
               allocate (Yc(i)%mat(product(A%sons(1, i)%n), R))
               Yc(i)%mat = 0
            enddo
         endif
         do i = 1, A%nch
            !$omp task default(shared) firstprivate(i) private(j) if(product(A%m) >= htlu_task_rows)
            do j = 1, A%nch
               if (op == 'N') then
                  call htlu_mult(A%sons(i, j), 'N', Xc(j)%mat, Yc(i)%mat, alpha)
               else
                  call htlu_mult(A%sons(j, i), 'T', Xc(j)%mat, Yc(i)%mat, alpha)
               endif
            enddo
            !$omp end task
         enddo
         !$omp taskwait
         if (op == 'N') then
            call htlu_scatter(A, 'r', Yc, Y, 2)
         else
            call htlu_scatter(A, 'c', Yc, Y, 2)
         endif
         do i = 1, A%nch
            if (allocated(Xc(i)%mat)) deallocate (Xc(i)%mat)
         enddo
         deallocate (Xc, Yc)
      end select
   end subroutine htlu_mult

   !>**** Y = sum_t coef_t*op(a_t*b_t)*X (Y zeroed here)
   subroutine htlu_terms_apply(nt, terms, op, X, Y)
      integer nt, t
      type(htlu_term) :: terms(:)
      character op
      DT :: X(:, :), Y(:, :)
      DT, allocatable :: W(:, :)
      Y = 0
      do t = 1, nt
         if (associated(terms(t)%b)) then
            if (op == 'N') then
               allocate (W(product(terms(t)%b%m), size(X, 2)))
               W = 0
               call htlu_mult(terms(t)%b, 'N', X, W, BPACK_cone)
               call htlu_mult(terms(t)%a, 'N', W, Y, terms(t)%coef)
            else
               allocate (W(product(terms(t)%a%n), size(X, 2)))
               W = 0
               call htlu_mult(terms(t)%a, 'T', X, W, BPACK_cone)
               call htlu_mult(terms(t)%b, 'T', W, Y, terms(t)%coef)
            endif
            deallocate (W)
         else
            call htlu_mult(terms(t)%a, op, X, Y, terms(t)%coef)
         endif
      enddo
   end subroutine htlu_terms_apply

   !>**** per-mode bases of the range (op = 'N') or of the transposed range (op = 'T') of the operator sum_t terms_t, from
   !>**** random products; every product gives all the fibers of its output tensor, so a few vectors suffice. The number of
   !>**** vectors doubles until each basis is clearly smaller than its sample
   subroutine htlu_range_bases(Ndim, nt, terms, op, dims_out, dims_in, rguess, tol, Q)
      integer Ndim, nt, dims_out(Ndim), dims_in(Ndim), rguess, d, p, fib, rk
      type(htlu_term) :: terms(:)
      character op
      real(kind=8) tol(2)
      type(htlu_fmat) :: Q(:)
      DT, allocatable :: Om(:, :), Y(:, :)
      logical ok
      fib = huge(1)
      do d = 1, Ndim
         fib = min(fib, product(dims_out)/dims_out(d))
      enddo
      p = max(2, ceiling(dble(2*rguess + 8)/dble(max(1, fib))))
      do
         allocate (Om(product(dims_in), p), Y(product(dims_out), p))
         call htlu_random(Om)
         call htlu_terms_apply(nt, terms, op, Om, Y)
         call htlu_mode_bases(Ndim, dims_out, Y, (/tol(1), 0d0/), Q)
         deallocate (Om, Y)
         ok = .true.
         do d = 1, Ndim
            rk = size(Q(d)%mat, 2)
            if (rk < dims_out(d) .and. 2*rk > (product(dims_out)/dims_out(d))*p) ok = .false.
         enddo
         if (ok) exit
         p = 2*p
      enddo
   end subroutine htlu_range_bases

   !>**** the Tucker block C (style 2) becomes the truncated Tucker form of sum_t terms_t (whose terms may refer to C itself):
   !>**** bases from random products of the operator and its transpose, core = Q^H * op * conj(Qc) exactly
   subroutine htlu_tucker_rebuild(C, nt, terms, tol)
      type(htlu_node) :: C
      integer nt
      type(htlu_term) :: terms(:)
      real(kind=8) tol(2)
      type(htlu_fmat), allocatable :: Q(:), Qc(:)
      DT, allocatable :: K(:, :), T(:, :), core(:, :), coreT(:, :)
      integer Ndim, d, rg, pr, pc, kind
      real(kind=8) t0
      t0 = MPI_Wtime()
      Ndim = C%Ndim
      rg = 1
      do d = 1, Ndim
         if (allocated(C%U)) rg = max(rg, size(C%U(d)%mat, 2), size(C%V(d)%mat, 2))
      enddo
      allocate (Q(Ndim), Qc(Ndim))
      call htlu_range_bases(Ndim, nt, terms, 'N', C%m, C%n, rg, tol, Q)
      call htlu_range_bases(Ndim, nt, terms, 'T', C%n, C%m, rg, tol, Qc)
      pr = 1
      pc = 1
      do d = 1, Ndim
         pr = pr*size(Q(d)%mat, 2)
         pc = pc*size(Qc(d)%mat, 2)
      enddo
      if (pc <= pr) then
         call htlu_kron(Ndim, Qc, K, 1)
         allocate (T(product(C%m), pc))
         call htlu_terms_apply(nt, terms, 'N', K, T)
         call htlu_modes(Ndim, C%m, T, Q, 'C', core)
      else
         call htlu_kron(Ndim, Q, K, 1)
         allocate (T(product(C%n), pr))
         call htlu_terms_apply(nt, terms, 'T', K, T)
         call htlu_modes(Ndim, C%n, T, Qc, 'C', coreT)
         allocate (core(pr, pc))
         core = transpose(coreT)
         deallocate (coreT)
      endif
      deallocate (K, T)
      if (allocated(C%U)) deallocate (C%U)
      if (allocated(C%V)) deallocate (C%V)
      call move_alloc(Q, C%U)
      call move_alloc(Qc, C%V)
      call move_alloc(core, C%C)
      call htlu_core_trim(C, tol)
      kind = 2
      do d = 1, nt
         if (associated(terms(d)%b)) kind = 1
      enddo
      t0 = MPI_Wtime() - t0
      if (kind == 1) then
         call htlu_tadd(1, t0)
      else
         call htlu_tadd(6, t0)
      endif
      d = min(15, C%level)
      !$omp atomic
      htlu_rlev_t(d, kind) = htlu_rlev_t(d, kind) + t0
      !$omp atomic
      htlu_rlev_n(d, kind) = htlu_rlev_n(d, kind) + 1
   end subroutine htlu_tucker_rebuild

   !>**** P = A*B in Tucker form, for A or B a Tucker block (P is a new block on the rows of A and the columns of B)
   subroutine htlu_tucker_product(A, B, P, tol)
      type(htlu_node) :: A, B, P
      real(kind=8) tol(2)
      type(htlu_fmat), allocatable :: S(:), Q(:)
      DT, allocatable :: K(:, :), W(:, :), G(:, :), T(:, :)
      integer Ndim, d
      real(kind=8) t0
      t0 = MPI_Wtime()
      Ndim = A%Ndim
      P%style = 2
      P%Ndim = Ndim
      P%level = A%level
      allocate (P%m(Ndim), P%n(Ndim), P%U(Ndim), P%V(Ndim))
      P%m = A%m
      P%n = B%n
      if (A%style == 2 .and. B%style == 2) then
         allocate (S(Ndim))
         do d = 1, Ndim
            call htlu_gemm('T', 'N', A%V(d)%mat, B%U(d)%mat, S(d)%mat)
         enddo
         call htlu_kron(Ndim, S, K, 0)
         call htlu_gemm('N', 'N', K, B%C, T)
         call htlu_gemm('N', 'N', A%C, T, P%C)
         do d = 1, Ndim
            allocate (P%U(d)%mat(size(A%U(d)%mat, 1), size(A%U(d)%mat, 2)))
            P%U(d)%mat = A%U(d)%mat
            allocate (P%V(d)%mat(size(B%V(d)%mat, 1), size(B%V(d)%mat, 2)))
            P%V(d)%mat = B%V(d)%mat
            deallocate (S(d)%mat)
         enddo
         deallocate (S, K, T)
      elseif (A%style == 2) then
         ! P = kron(QU) Z^T with U_A = QU RU and Z = (B^T kron(V_A)) (kron(RU) C_A)^T ~ kron(Q) G: the truncation sees the block
         allocate (Q(Ndim), S(Ndim))
         call htlu_orth(Ndim, A%U, S, T)
         call htlu_kron(Ndim, A%V, K, 0)
         allocate (W(product(B%n), size(K, 2)))
         W = 0
         call htlu_mult(B, 'T', K, W, BPACK_cone)
         call htlu_gemm('N', 'N', T, A%C, K)
         call htlu_gemm('N', 'T', W, K, T)
         call htlu_tuckerize(Ndim, B%n, T, tol, Q, G)
         allocate (P%C(size(G, 2), size(G, 1)))
         P%C = transpose(G)
         do d = 1, Ndim
            call move_alloc(S(d)%mat, P%U(d)%mat)
            call move_alloc(Q(d)%mat, P%V(d)%mat)
         enddo
         deallocate (Q, S, K, W, G, T)
      else
         call assert(B%style == 2, 'htlu_tucker_product: no Tucker operand')
         ! P = Z kron(QV)^T with V_B = QV RV and Z = (A kron(U_B)) C_B kron(RV)^T ~ kron(Q) G
         allocate (Q(Ndim), S(Ndim))
         call htlu_orth(Ndim, B%V, S, T)
         call htlu_kron(Ndim, B%U, K, 0)
         allocate (W(product(A%m), size(K, 2)))
         W = 0
         call htlu_mult(A, 'N', K, W, BPACK_cone)
         call htlu_gemm('N', 'T', B%C, T, K)
         call htlu_gemm('N', 'N', W, K, T)
         call htlu_tuckerize(Ndim, A%m, T, tol, Q, G)
         call move_alloc(G, P%C)
         do d = 1, Ndim
            call move_alloc(Q(d)%mat, P%U(d)%mat)
            call move_alloc(S(d)%mat, P%V(d)%mat)
         enddo
         deallocate (Q, S, K, W, T)
      endif
      if (A%style /= 2 .or. B%style /= 2) call htlu_core_trim(P, tol)
      call htlu_tadd(3, MPI_Wtime() - t0)
   end subroutine htlu_tucker_product

   !>**** Pc = the restriction of the Tucker block P to row child i and column child j of the hierarchical block A
   subroutine htlu_restrict(P, A, i, j, Pc)
      type(htlu_node) :: P, A, Pc
      integer i, j
      integer ro(A%Ndim), rs(A%Ndim), co(A%Ndim), cs(A%Ndim)
      call htlu_child_box(A, 'r', i, ro, rs)
      call htlu_child_box(A, 'c', j, co, cs)
      call htlu_restrict_box(P, ro, rs, co, cs, A%level + 1, Pc)
   end subroutine htlu_restrict

   !>**** Pc (tree level lev) = the restriction of the Tucker block P to the rows ro+1..ro+rs and the columns co+1..co+cs per dimension
   subroutine htlu_restrict_box(P, ro, rs, co, cs, lev, Pc)
      type(htlu_node) :: P, Pc
      integer ro(:), rs(:), co(:), cs(:), lev, d, Ndim
      Ndim = P%Ndim
      Pc%style = 2
      Pc%Ndim = Ndim
      Pc%level = lev
      allocate (Pc%m(Ndim), Pc%n(Ndim), Pc%U(Ndim), Pc%V(Ndim))
      Pc%m = rs
      Pc%n = cs
      do d = 1, Ndim
         allocate (Pc%U(d)%mat(rs(d), size(P%U(d)%mat, 2)))
         Pc%U(d)%mat = P%U(d)%mat(ro(d) + 1:ro(d) + rs(d), :)
         allocate (Pc%V(d)%mat(cs(d), size(P%V(d)%mat, 2)))
         Pc%V(d)%mat = P%V(d)%mat(co(d) + 1:co(d) + cs(d), :)
      enddo
      allocate (Pc%C(size(P%C, 1), size(P%C, 2)))
      Pc%C = P%C
   end subroutine htlu_restrict_box

   !>**** thin QR A = Q R (Q: m x k orthonormal columns, R: k x n, k = min(m, n))
   subroutine htlu_qr(A, Q, R)
      DT :: A(:, :)
      DT, allocatable :: Q(:, :), R(:, :), W(:, :), tau(:)
      integer m, n, k, i, j
      m = size(A, 1)
      n = size(A, 2)
      k = min(m, n)
      allocate (W(m, n), tau(max(1, k)))
      W = A
      if (allocated(Q)) deallocate (Q)
      if (allocated(R)) deallocate (R)
      allocate (Q(m, k), R(k, n))
      R = 0
      if (k > 0) then
         call geqrff90(W, tau)
         do j = 1, n
            do i = 1, min(j, k)
               R(i, j) = W(i, j)
            enddo
         enddo
         call un_or_gqrf90(W, tau, m, k, k)
         Q = W(:, 1:k)
      endif
      deallocate (W, tau)
   end subroutine htlu_qr

   !>**** C = C + alpha*P for Tucker blocks C and P on the same clusters, without random products: per dimension the thin QR
   !>**** [U_C U_P] = Q [RA RB] and [V_C V_P] = Qc [RcA RcB], the core G = kron(RA) C_C kron(RcA)^T + alpha kron(RB) C_P
   !>**** kron(RcB)^T, then the truncated HOSVD of G. Returns ok = .false. (C unchanged) if G would exceed htlu_sum_cap entries
   subroutine htlu_tucker_sum(C, P, alpha, tol, ok)
      type(htlu_node) :: C, P
      DT alpha
      real(kind=8) tol(2)
      logical ok
      type(htlu_fmat), allocatable :: Q(:), Qc(:), RA(:), RB(:), RcA(:), RcB(:)
      DT, allocatable :: Cat(:, :), R(:, :), KA(:, :), KB(:, :), T(:, :), G(:, :)
      integer Ndim, d, na, nb, k, kc
      integer(kind=8) gsize
      real(kind=8) t0
      t0 = MPI_Wtime()
      Ndim = C%Ndim
      gsize = 1
      do d = 1, Ndim
         k = min(C%m(d), size(C%U(d)%mat, 2) + size(P%U(d)%mat, 2))
         kc = min(C%n(d), size(C%V(d)%mat, 2) + size(P%V(d)%mat, 2))
         gsize = gsize*k*kc
      enddo
      ok = gsize <= htlu_sum_cap
      if (.not. ok) return
      allocate (Q(Ndim), Qc(Ndim), RA(Ndim), RB(Ndim), RcA(Ndim), RcB(Ndim))
      do d = 1, Ndim
         na = size(C%U(d)%mat, 2)
         nb = size(P%U(d)%mat, 2)
         allocate (Cat(C%m(d), na + nb))
         Cat(:, 1:na) = C%U(d)%mat
         Cat(:, na + 1:na + nb) = P%U(d)%mat
         call htlu_qr(Cat, Q(d)%mat, R)
         allocate (RA(d)%mat(size(R, 1), na), RB(d)%mat(size(R, 1), nb))
         RA(d)%mat = R(:, 1:na)
         RB(d)%mat = R(:, na + 1:na + nb)
         deallocate (Cat, R)
         na = size(C%V(d)%mat, 2)
         nb = size(P%V(d)%mat, 2)
         allocate (Cat(C%n(d), na + nb))
         Cat(:, 1:na) = C%V(d)%mat
         Cat(:, na + 1:na + nb) = P%V(d)%mat
         call htlu_qr(Cat, Qc(d)%mat, R)
         allocate (RcA(d)%mat(size(R, 1), na), RcB(d)%mat(size(R, 1), nb))
         RcA(d)%mat = R(:, 1:na)
         RcB(d)%mat = R(:, na + 1:na + nb)
         deallocate (Cat, R)
      enddo
      ! G = kron(RA) C_C kron(RcA)^T + alpha kron(RB) C_P kron(RcB)^T
      call htlu_kron(Ndim, RA, KA, 0)
      call htlu_kron(Ndim, RcA, KB, 0)
      call htlu_gemm('N', 'T', C%C, KB, T)
      call htlu_gemm('N', 'N', KA, T, G)
      deallocate (KA, KB, T)
      call htlu_kron(Ndim, RB, KA, 0)
      call htlu_kron(Ndim, RcB, KB, 0)
      call htlu_gemm('N', 'T', P%C, KB, T)
      call htlu_gemm('N', 'N', KA, T, R)
      G = G + alpha*R
      deallocate (KA, KB, T, R)
      do d = 1, Ndim
         call move_alloc(Q(d)%mat, C%U(d)%mat)
         call move_alloc(Qc(d)%mat, C%V(d)%mat)
         deallocate (RA(d)%mat, RB(d)%mat, RcA(d)%mat, RcB(d)%mat)
      enddo
      call move_alloc(G, C%C)
      deallocate (Q, Qc, RA, RB, RcA, RcB)
      call htlu_core_trim(C, tol)
      call htlu_tadd(7, MPI_Wtime() - t0)
   end subroutine htlu_tucker_sum

   !>**** C = C + alpha*P for a block C of any style and a Tucker block P on the same clusters
   recursive subroutine htlu_add_tucker(C, P, alpha, tol)
      type(htlu_node), target :: C, P
      DT alpha
      real(kind=8) tol(2)
      type(htlu_term) :: terms(2)
      type(htlu_node), allocatable :: Pc(:, :)
      DT, allocatable :: D(:, :)
      integer i, j
      logical ok
      select case (C%style)
      case (1)
         call htlu_dense_of(P, D)
         C%D = C%D + alpha*D
         deallocate (D)
      case (2)
         call htlu_tucker_sum(C, P, alpha, tol, ok)
         if (.not. ok) then
            terms(1)%coef = BPACK_cone
            terms(1)%a => C
            terms(2)%coef = alpha
            terms(2)%a => P
            call htlu_tucker_rebuild(C, 2, terms, tol)
         endif
      case (4)
         allocate (Pc(C%nch, C%nch))
         do j = 1, C%nch
            do i = 1, C%nch
               !$omp task default(shared) firstprivate(i,j) if(product(C%m) >= htlu_task_rows)
               call htlu_restrict(P, C, i, j, Pc(i, j))
               call htlu_add_tucker(C%sons(i, j), Pc(i, j), alpha, tol)
               call HTLU_delete(Pc(i, j))
               !$omp end task
            enddo
         enddo
         !$omp taskwait
         deallocate (Pc)
      end select
   end subroutine htlu_add_tucker

   !>**** C = C - A*B (Schur update)
   recursive subroutine htlu_subhh(C, A, B, tol)
      type(htlu_node), target :: C, A, B
      real(kind=8) tol(2)
      type(htlu_term) :: terms(2)
      type(htlu_node) :: P
      DT, allocatable :: Bd(:, :)
      integer i, j, k
      real(kind=8) t0
      logical ok
      select case (C%style)
      case (1)
         t0 = MPI_Wtime()
         call htlu_dense_of(B, Bd)
         call htlu_mult(A, 'N', Bd, C%D, -BPACK_cone)
         deallocate (Bd)
         call htlu_tadd(5, MPI_Wtime() - t0)
      case (2)
         ! both operands Tucker: closed-form product and deterministic sum; otherwise the randomized rebuild, which is cheaper
         ! than the weighted product with a hierarchical operand followed by the sum
         ok = .false.
         if (A%style == 2 .and. B%style == 2) then
            call htlu_tucker_product(A, B, P, tol)
            call htlu_tucker_sum(C, P, -BPACK_cone, tol, ok)
            call HTLU_delete(P)
         endif
         if (.not. ok) then
            terms(1)%coef = BPACK_cone
            terms(1)%a => C
            terms(2)%coef = -BPACK_cone
            terms(2)%a => A
            terms(2)%b => B
            call htlu_tucker_rebuild(C, 2, terms, tol)
         endif
      case (4)
         if (A%style == 4 .and. B%style == 4) then
            do j = 1, C%nch
               do i = 1, C%nch
                  !$omp task default(shared) firstprivate(i,j) private(k) if(product(C%m) >= htlu_task_rows)
                  do k = 1, C%nch
                     call htlu_subhh(C%sons(i, j), A%sons(i, k), B%sons(k, j), tol)
                  enddo
                  !$omp end task
               enddo
            enddo
            !$omp taskwait
         else
            call htlu_tucker_product(A, B, P, tol)
            call htlu_add_tucker(C, P, -BPACK_cone, tol)
            call HTLU_delete(P)
         endif
      end select
   end subroutine htlu_subhh

   !>**** block LU with whole-cluster pivots (htlu_BLU) of a diagonal block D = Lb Ub: the lower sons hold X_ji = At_ji S_ii^{-1}
   !>**** (Lb block unit lower triangular), the upper sons the Schur-updated blocks At_ij and the diagonal sons the factored
   !>**** current diagonal blocks S_ii (Ub block upper triangular); dense leaves hold their partially pivoted LU.
   !>**** X = D^{-1} X: y_j = x_j - sum_{i<j} X_ji y_i, then x_i = S_ii^{-1} (y_i - sum_{l>i} At_il x_l) for i = K..1
   recursive subroutine htlu_bsolve(D, X)
      type(htlu_node) :: D
      DT :: X(:, :)
      type(htlu_fmat), allocatable :: Xc(:)
      integer i, j
      if (size(X, 2) == 0) return
      if (D%style == 1) then
         call getrsf90(D%D, D%ipiv, X, 'N')
      else
         call assert(D%style == 4, 'htlu_bsolve: a Tucker diagonal block')
         allocate (Xc(D%nch))
         call htlu_gather(D, 'r', X, Xc)
         do j = 2, D%nch
            do i = 1, j - 1
               call htlu_mult(D%sons(j, i), 'N', Xc(i)%mat, Xc(j)%mat, -BPACK_cone)
            enddo
         enddo
         do i = D%nch, 1, -1
            do j = i + 1, D%nch
               call htlu_mult(D%sons(i, j), 'N', Xc(j)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
            call htlu_bsolve(D%sons(i, i), Xc(i)%mat)
         enddo
         call htlu_scatter(D, 'r', Xc, X, 1)
         deallocate (Xc)
      endif
   end subroutine htlu_bsolve

   !>**** X = D^{-T} X for the block-LU factored diagonal block D: D^{-T} = Lb^{-T} Ub^{-T}, i.e. z_i = S_ii^{-T} (x_i - sum_{l<i}
   !>**** At_li^T z_l) for i = 1..K, then y_i = z_i - sum_{j>i} X_ji^T y_j for i = K..1
   recursive subroutine htlu_bsolve_T(D, X)
      type(htlu_node) :: D
      DT :: X(:, :)
      type(htlu_fmat), allocatable :: Xc(:)
      integer i, j
      if (size(X, 2) == 0) return
      if (D%style == 1) then
         call getrsf90(D%D, D%ipiv, X, 'T')
      else
         call assert(D%style == 4, 'htlu_bsolve_T: a Tucker diagonal block')
         allocate (Xc(D%nch))
         call htlu_gather(D, 'c', X, Xc)
         do i = 1, D%nch
            do j = 1, i - 1
               call htlu_mult(D%sons(j, i), 'T', Xc(j)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
            call htlu_bsolve_T(D%sons(i, i), Xc(i)%mat)
         enddo
         do i = D%nch - 1, 1, -1
            do j = i + 1, D%nch
               call htlu_mult(D%sons(j, i), 'T', Xc(j)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
         enddo
         call htlu_scatter(D, 'c', Xc, X, 1)
         deallocate (Xc)
      endif
   end subroutine htlu_bsolve_T

   !>**** M = M D^{-1} for the block-LU factored diagonal block D of the column clusters of M
   recursive subroutine htlu_XSM(D, M, tol)
      type(htlu_node), target :: D, M
      real(kind=8) tol(2)
      type(htlu_fmat), allocatable :: Q(:), S(:)
      DT, allocatable :: W(:, :), G(:, :), Z0(:, :), Rk(:, :), T(:, :)
      integer a, b, c, k
      real(kind=8) t0
      t0 = MPI_Wtime()
      select case (M%style)
      case (1)
         call assert(D%style == 1, 'htlu_XSM: a dense block with a hierarchical factor')
         allocate (T(size(M%D, 2), size(M%D, 1)))
         T = transpose(M%D)
         call getrsf90(D%D, D%ipiv, T, 'T')
         M%D = transpose(T)
         deallocate (T)
         call htlu_tadd(5, MPI_Wtime() - t0)
      case (2)
         ! M D^{-1} = kron(QU) Z^T with U_M = QU RU and Z = (D^{-T} kron(V_M)) (kron(RU) C_M)^T ~ kron(Q) G
         allocate (Q(M%Ndim), S(M%Ndim))
         call htlu_orth(M%Ndim, M%U, S, Rk)
         call htlu_gemm('N', 'N', Rk, M%C, Z0)
         call htlu_kron(M%Ndim, M%V, W, 0)
         call htlu_bsolve_T(D, W)
         call htlu_gemm('N', 'T', W, Z0, Rk)
         call htlu_tuckerize(M%Ndim, M%n, Rk, tol, Q, G)
         do k = 1, M%Ndim
            call move_alloc(S(k)%mat, M%U(k)%mat)
            call move_alloc(Q(k)%mat, M%V(k)%mat)
         enddo
         deallocate (M%C)
         allocate (M%C(size(G, 2), size(G, 1)))
         M%C = transpose(G)
         deallocate (Q, S, W, Z0, Rk, G)
         call htlu_core_trim(M, tol)
         call htlu_tadd(2, MPI_Wtime() - t0)
      case (4)
         call assert(D%style == 4, 'htlu_XSM: a hierarchical block with a leaf factor')
         ! M D^{-1} = (M Ub^{-1}) Lb^{-1}, row son by row son
         do a = 1, M%nch
            !$omp task default(shared) firstprivate(a) private(b,c) if(product(M%m) >= htlu_task_rows)
            do b = 1, M%nch
               do c = 1, b - 1
                  call htlu_subhh(M%sons(a, b), M%sons(a, c), D%sons(c, b), tol)
               enddo
               call htlu_XSM(D%sons(b, b), M%sons(a, b), tol)
            enddo
            do b = M%nch - 1, 1, -1
               do c = b + 1, M%nch
                  call htlu_subhh(M%sons(a, b), M%sons(a, c), D%sons(c, b), tol)
               enddo
            enddo
            !$omp end task
         enddo
         !$omp taskwait
      end select
   end subroutine htlu_XSM

   !>**** block LU with whole-cluster pivots of the diagonal block A in place (see htlu_bsolve): for i = 1..K, factor A_ii,
   !>**** A_ji = A_ji A_ii^{-1} (j > i), A_jl = A_jl - A_ji A_il (j, l > i); the upper blocks A_ij keep their Schur-updated values
   recursive subroutine htlu_BLU(A, tol)
      type(htlu_node), target :: A
      real(kind=8) tol(2)
      integer i, j, l
      real(kind=8) t0
      if (A%style == 1) then
         if (allocated(A%ipiv)) deallocate (A%ipiv)
         allocate (A%ipiv(size(A%D, 1)))
         t0 = MPI_Wtime()
         call getrff90(A%D, A%ipiv)
         call htlu_tadd(5, MPI_Wtime() - t0)
         return
      endif
      call assert(A%style == 4, 'htlu_BLU: a Tucker diagonal block')
      do i = 1, A%nch
         call htlu_BLU(A%sons(i, i), tol)
         do j = i + 1, A%nch
            !$omp task default(shared) firstprivate(i,j) if(product(A%m) >= htlu_task_rows)
            call htlu_XSM(A%sons(i, i), A%sons(j, i), tol)
            !$omp end task
         enddo
         !$omp taskwait
         do l = i + 1, A%nch
            do j = i + 1, A%nch
               !$omp task default(shared) firstprivate(i,j,l) if(product(A%m) >= htlu_task_rows)
               call htlu_subhh(A%sons(j, l), A%sons(j, i), A%sons(i, l), tol)
               !$omp end task
            enddo
         enddo
         !$omp taskwait
      enddo
   end subroutine htlu_BLU

   !>**** A = alpha*A for a block of any style
   recursive subroutine htlu_scale(A, alpha)
      type(htlu_node) :: A
      DT alpha
      integer i, j
      select case (A%style)
      case (1)
         A%D = alpha*A%D
      case (2)
         A%C = alpha*A%C
      case (4)
         do j = 1, A%nch
            do i = 1, A%nch
               call htlu_scale(A%sons(i, j), alpha)
            enddo
         enddo
      end select
   end subroutine htlu_scale

   !>**** Y = a block with the structure of A and zero entries
   recursive subroutine htlu_zero_clone(A, Y)
      type(htlu_node) :: A, Y
      integer i, j, d
      Y%style = A%style
      Y%level = A%level
      Y%Ndim = A%Ndim
      Y%nch = A%nch
      allocate (Y%m(A%Ndim), Y%n(A%Ndim), Y%rg(A%Ndim), Y%cg(A%Ndim))
      Y%m = A%m
      Y%n = A%n
      Y%rg = A%rg
      Y%cg = A%cg
      select case (A%style)
      case (1)
         allocate (Y%D(size(A%D, 1), size(A%D, 2)))
         Y%D = 0
      case (2)
         allocate (Y%U(A%Ndim), Y%V(A%Ndim), Y%C(1, 1))
         do d = 1, A%Ndim
            allocate (Y%U(d)%mat(A%m(d), 1), Y%V(d)%mat(A%n(d), 1))
            Y%U(d)%mat = 0
            Y%V(d)%mat = 0
         enddo
         Y%C = 0
      case (4)
         allocate (Y%sons(A%nch, A%nch))
         do j = 1, A%nch
            do i = 1, A%nch
               call htlu_zero_clone(A%sons(i, j), Y%sons(i, j))
            enddo
         enddo
      end select
   end subroutine htlu_zero_clone

   !>**** X = X*P (side 'R') or X = P*X (side 'L') for the square block P on the column (row) clusters of X
   recursive subroutine htlu_inv_mult(X, P, side, tol)
      type(htlu_node), target :: X, P
      character side
      real(kind=8) tol(2)
      type(htlu_node), target :: Y
      type(htlu_fmat), allocatable :: Q(:), S(:)
      DT, allocatable :: W(:, :), W2(:, :), G(:, :), Z0(:, :), Rk(:, :)
      integer i, j, k, d
      real(kind=8) t0
      t0 = MPI_Wtime()
      select case (X%style)
      case (1)
         call assert(P%style == 1, 'htlu_inv_mult: a dense block with a hierarchical factor')
         if (side == 'R') then
            call htlu_gemm('N', 'N', X%D, P%D, W)
         else
            call htlu_gemm('N', 'N', P%D, X%D, W)
         endif
         call move_alloc(W, X%D)
         call htlu_tadd(5, MPI_Wtime() - t0)
      case (2)
         allocate (Q(X%Ndim), S(X%Ndim))
         if (side == 'R') then
            ! X P = kron(QU) Z^T with U_X = QU RU and Z = (P^T kron(V_X)) (kron(RU) C_X)^T ~ kron(Q) G
            call htlu_orth(X%Ndim, X%U, S, Rk)
            call htlu_gemm('N', 'N', Rk, X%C, Z0)
            call htlu_kron(X%Ndim, X%V, W, 0)
            allocate (W2(product(X%n), size(W, 2)))
            W2 = 0
            call htlu_mult(P, 'T', W, W2, BPACK_cone)
            call htlu_gemm('N', 'T', W2, Z0, Rk)
            call htlu_tuckerize(X%Ndim, X%n, Rk, tol, Q, G)
            do d = 1, X%Ndim
               call move_alloc(S(d)%mat, X%U(d)%mat)
               call move_alloc(Q(d)%mat, X%V(d)%mat)
            enddo
            deallocate (X%C)
            allocate (X%C(size(G, 2), size(G, 1)))
            X%C = transpose(G)
         else
            ! P X = Z kron(QV)^T with V_X = QV RV and Z = (P kron(U_X)) C_X kron(RV)^T ~ kron(Q) G
            call htlu_orth(X%Ndim, X%V, S, Rk)
            call htlu_gemm('N', 'T', X%C, Rk, Z0)
            call htlu_kron(X%Ndim, X%U, W, 0)
            allocate (W2(product(X%m), size(W, 2)))
            W2 = 0
            call htlu_mult(P, 'N', W, W2, BPACK_cone)
            call htlu_gemm('N', 'N', W2, Z0, Rk)
            call htlu_tuckerize(X%Ndim, X%m, Rk, tol, Q, G)
            do d = 1, X%Ndim
               call move_alloc(Q(d)%mat, X%U(d)%mat)
               call move_alloc(S(d)%mat, X%V(d)%mat)
            enddo
            call move_alloc(G, X%C)
         endif
         deallocate (Q, S, W, W2, Z0, Rk)
         if (allocated(G)) deallocate (G)
         call htlu_core_trim(X, tol)
         call htlu_tadd(4, MPI_Wtime() - t0)
      case (4)
         call assert(P%style == 4, 'htlu_inv_mult: a hierarchical block with a leaf factor')
         ! Y = -X*P (or -P*X) accumulated block by block in a zero copy of X, which then replaces X
         call htlu_zero_clone(X, Y)
         do j = 1, X%nch
            do i = 1, X%nch
               !$omp task default(shared) firstprivate(i,j) private(k) if(product(X%m) >= htlu_task_rows)
               do k = 1, X%nch
                  if (side == 'R') then
                     call htlu_subhh(Y%sons(i, j), X%sons(i, k), P%sons(k, j), tol)
                  else
                     call htlu_subhh(Y%sons(i, j), P%sons(i, k), X%sons(k, j), tol)
                  endif
               enddo
               !$omp end task
            enddo
         enddo
         !$omp taskwait
         call htlu_scale(Y, -BPACK_cone)
         do j = 1, X%nch
            do i = 1, X%nch
               call HTLU_delete(X%sons(i, j))
            enddo
         enddo
         deallocate (X%sons)
         X%sons => Y%sons
         nullify (Y%sons)
         call HTLU_delete(Y)
      end select
   end subroutine htlu_inv_mult

   !>**** inverse of a dense leaf in place
   subroutine htlu_dense_inv(A)
      type(htlu_node) :: A
      DT, allocatable :: LU(:, :)
      integer, allocatable :: ipiv(:)
      integer i, n
      real(kind=8) t0
      t0 = MPI_Wtime()
      n = size(A%D, 1)
      allocate (LU(n, n), ipiv(n))
      LU = A%D
      call getrff90(LU, ipiv)
      A%D = 0
      do i = 1, n
         A%D(i, i) = 1
      enddo
      call getrsf90(LU, ipiv, A%D, 'N')
      deallocate (LU, ipiv)
      call htlu_tadd(5, MPI_Wtime() - t0)
   end subroutine htlu_dense_inv

   !>**** block Gauss-Jordan sweep of the K x K array of blocks B (square diagonal blocks): B = -B^{-1}. Step k:
   !>**** B_kk = inv(B_kk), B_jk = B_jk B_kk, B_jl = B_jl - B_jk B_kl, B_kl = B_kk B_kl, B_kk = -B_kk (j, l /= k)
   recursive subroutine htlu_sweep(K, B, tol)
      integer K
      type(htlu_node), target :: B(:, :)
      real(kind=8) tol(2)
      integer i, j, l
      do i = 1, K
         call htlu_inv(B(i, i), tol)
         do j = 1, K
            if (j == i) cycle
            !$omp task default(shared) firstprivate(i,j) if(product(B(1, 1)%m) >= htlu_task_rows)
            call htlu_inv_mult(B(j, i), B(i, i), 'R', tol)
            !$omp end task
         enddo
         !$omp taskwait
         do l = 1, K
            if (l == i) cycle
            do j = 1, K
               if (j == i) cycle
               !$omp task default(shared) firstprivate(i,j,l) if(product(B(1, 1)%m) >= htlu_task_rows)
               call htlu_subhh(B(j, l), B(j, i), B(i, l), tol)
               !$omp end task
            enddo
         enddo
         !$omp taskwait
         do l = 1, K
            if (l == i) cycle
            !$omp task default(shared) firstprivate(i,l) if(product(B(1, 1)%m) >= htlu_task_rows)
            call htlu_inv_mult(B(i, l), B(i, i), 'L', tol)
            !$omp end task
         enddo
         !$omp taskwait
         call htlu_scale(B(i, i), -BPACK_cone)
      enddo
   end subroutine htlu_sweep

   !>**** A = A^{-1} for a square block (dense inverse of a leaf, block Gauss-Jordan sweep over the sons otherwise)
   recursive subroutine htlu_inv(A, tol)
      type(htlu_node), target :: A
      real(kind=8) tol(2)
      if (A%style == 1) then
         call htlu_dense_inv(A)
         return
      endif
      call assert(A%style == 4, 'htlu_inv: a Tucker diagonal block')
      call htlu_sweep(A%nch, A%sons, tol)
      call htlu_scale(A, -BPACK_cone)
   end subroutine htlu_inv

   !>**** free a block and its sons
   recursive subroutine HTLU_delete(A)
      type(htlu_node) :: A
      integer i, j
      if (associated(A%sons)) then
         do j = 1, size(A%sons, 2)
            do i = 1, size(A%sons, 1)
               call HTLU_delete(A%sons(i, j))
            enddo
         enddo
         deallocate (A%sons)
      endif
      if (allocated(A%D)) deallocate (A%D)
      if (allocated(A%ipiv)) deallocate (A%ipiv)
      if (allocated(A%U)) deallocate (A%U)
      if (allocated(A%V)) deallocate (A%V)
      if (allocated(A%C)) deallocate (A%C)
      if (allocated(A%m)) deallocate (A%m)
      if (allocated(A%n)) deallocate (A%n)
      if (allocated(A%rg)) deallocate (A%rg)
      if (allocated(A%cg)) deallocate (A%cg)
      A%style = 0
   end subroutine HTLU_delete

   !>**** memory (MB) of the dense and the Tucker blocks of a tree, maximum Tucker rank per level, numbers of blocks
   recursive subroutine HTLU_stats(A, memd, memt, rk, nd, nt)
      type(htlu_node) :: A
      real(kind=8) memd, memt
      integer rk(0:), nd, nt, i, j, d
      real(kind=8) b
      b = dble(storage_size(BPACK_cone))/8d0/1024d0/1024d0
      select case (A%style)
      case (1)
         memd = memd + b*size(A%D)
         nd = nd + 1
      case (2)
         memt = memt + b*size(A%C)
         do d = 1, A%Ndim
            memt = memt + b*(size(A%U(d)%mat) + size(A%V(d)%mat))
            if (A%level <= ubound(rk, 1)) rk(A%level) = max(rk(A%level), size(A%U(d)%mat, 2), size(A%V(d)%mat, 2))
         enddo
         nt = nt + 1
      case (4)
         do j = 1, A%nch
            do i = 1, A%nch
               call HTLU_stats(A%sons(i, j), memd, memt, rk, nd, nt)
            enddo
         enddo
      end select
   end subroutine HTLU_stats

   !>**** Tucker blocks by level and by part of the matrix (part 1: strictly lower, 2: strictly upper, 3: on the diagonal path):
   !>**** count, sum and maximum of the largest row-side (side 1) and column-side (side 2) rank, memory (MB). part = 0 on entry
   !>**** marks a block on the diagonal path
   recursive subroutine HTLU_rank_profile(A, part, cnt, rsum, rmax, mem)
      type(htlu_node) :: A
      integer part, cnt(0:, :), rmax(0:, :, :), i, j, p2, l, d, rr, rc
      real(kind=8) rsum(0:, :, :), mem(0:, :), b
      b = dble(storage_size(BPACK_cone))/8d0/1024d0/1024d0
      select case (A%style)
      case (2)
         l = A%level
         if (l > ubound(cnt, 1)) return
         p2 = part
         if (p2 == 0) p2 = 3
         rr = 0
         rc = 0
         do d = 1, A%Ndim
            rr = max(rr, size(A%U(d)%mat, 2))
            rc = max(rc, size(A%V(d)%mat, 2))
            mem(l, p2) = mem(l, p2) + b*(size(A%U(d)%mat) + size(A%V(d)%mat))
         enddo
         mem(l, p2) = mem(l, p2) + b*size(A%C)
         cnt(l, p2) = cnt(l, p2) + 1
         rsum(l, p2, 1) = rsum(l, p2, 1) + rr
         rsum(l, p2, 2) = rsum(l, p2, 2) + rc
         rmax(l, p2, 1) = max(rmax(l, p2, 1), rr)
         rmax(l, p2, 2) = max(rmax(l, p2, 2), rc)
      case (4)
         do j = 1, A%nch
            do i = 1, A%nch
               p2 = part
               if (part == 0 .and. i > j) p2 = 1
               if (part == 0 .and. i < j) p2 = 2
               call HTLU_rank_profile(A%sons(i, j), p2, cnt, rsum, rmax, mem)
            enddo
         enddo
      end select
   end subroutine HTLU_rank_profile

   !>**** linear index (0-based) of the cluster of groups g(1:Ndim) at level l
   integer function htlu_cluster_index(Ndim, l, g)
      integer Ndim, l, g(Ndim), d, s
      htlu_cluster_index = 0
      s = 1
      do d = 1, Ndim
         htlu_cluster_index = htlu_cluster_index + (g(d) - 2**l)*s
         s = s*2**l
      enddo
   end function htlu_cluster_index

   !>**** index of the local blocks (dense or Tucker) of an HTENSOR by level and by row and column clusters
   subroutine htlu_index_build(h_mat_md, Ndim, idx, Lmax)
      type(Hmat_md), target :: h_mat_md
      integer Ndim, Lmax
      type(htlu_index_level), allocatable :: idx(:)
      type(matrixblock_MD), pointer :: blk
      integer ll, bb, l, ri, k
      Lmax = 0
      do ll = 1, h_mat_md%BP%Lplus
         do bb = 1, h_mat_md%BP%LL(ll)%Nbound_loc
            blk => h_mat_md%BP%LL(ll)%matrices_block(bb)
            if (blk%style /= 4) Lmax = max(Lmax, blk%level)
         enddo
      enddo
      allocate (idx(0:Lmax))
      do ll = 1, h_mat_md%BP%Lplus
         do bb = 1, h_mat_md%BP%LL(ll)%Nbound_loc
            blk => h_mat_md%BP%LL(ll)%matrices_block(bb)
            if (blk%style /= 4) idx(blk%level)%nb = idx(blk%level)%nb + 1
         enddo
      enddo
      do l = 0, Lmax
         allocate (idx(l)%head(0:2**(l*Ndim) - 1), idx(l)%next(max(1, idx(l)%nb)), idx(l)%cidx(max(1, idx(l)%nb)), &
            idx(l)%ll(max(1, idx(l)%nb)), idx(l)%bb(max(1, idx(l)%nb)))
         idx(l)%head = 0
         idx(l)%nb = 0
      enddo
      do ll = 1, h_mat_md%BP%Lplus
         do bb = 1, h_mat_md%BP%LL(ll)%Nbound_loc
            blk => h_mat_md%BP%LL(ll)%matrices_block(bb)
            if (blk%style == 4) cycle
            l = blk%level
            ri = htlu_cluster_index(Ndim, l, blk%row_group)
            idx(l)%nb = idx(l)%nb + 1
            k = idx(l)%nb
            idx(l)%cidx(k) = htlu_cluster_index(Ndim, l, blk%col_group)
            idx(l)%ll(k) = ll
            idx(l)%bb(k) = bb
            idx(l)%next(k) = idx(l)%head(ri)
            idx(l)%head(ri) = k
         enddo
      enddo
   end subroutine htlu_index_build

   !>**** free the index of htlu_index_build
   subroutine htlu_index_free(idx)
      type(htlu_index_level), allocatable :: idx(:)
      integer l
      do l = lbound(idx, 1), ubound(idx, 1)
         deallocate (idx(l)%head, idx(l)%next, idx(l)%cidx, idx(l)%ll, idx(l)%bb)
      enddo
      deallocate (idx)
   end subroutine htlu_index_free

   !>**** groups g(1:Ndim) of the cluster of linear index t (from 0, dimension 1 fastest) at level l
   subroutine htlu_cluster_groups(Ndim, l, t, g)
      integer Ndim, l, t, g(Ndim), d, s
      s = t
      do d = 1, Ndim
         g(d) = 2**l + mod(s, 2**l)
         s = s/2**l
      enddo
   end subroutine htlu_cluster_groups

   recursive subroutine htlu_build(A, Ndim, l, rg, cg, msh, h_mat_md, idx, Lmax)
      type(htlu_node) :: A
      integer Ndim, l, rg(Ndim), cg(Ndim), Lmax
      type(mesh) :: msh(Ndim)
      type(Hmat_md), target :: h_mat_md
      type(htlu_index_level) :: idx(0:)
      type(matrixblock_MD), pointer :: blk
      integer d, k, ri, ci, i, j, rgc(Ndim), cgc(Ndim), c1(Ndim), c2(Ndim)
      A%Ndim = Ndim
      A%level = l
      allocate (A%m(Ndim), A%n(Ndim), A%rg(Ndim), A%cg(Ndim))
      A%rg = rg
      A%cg = cg
      do d = 1, Ndim
         A%m(d) = msh(d)%basis_group(rg(d))%tail - msh(d)%basis_group(rg(d))%head + 1
         A%n(d) = msh(d)%basis_group(cg(d))%tail - msh(d)%basis_group(cg(d))%head + 1
      enddo
      k = 0
      if (l <= Lmax) then
         ri = htlu_cluster_index(Ndim, l, rg)
         ci = htlu_cluster_index(Ndim, l, cg)
         k = idx(l)%head(ri)
         do while (k /= 0)
            if (idx(l)%cidx(k) == ci) exit
            k = idx(l)%next(k)
         enddo
      endif
      if (k /= 0) then
         blk => h_mat_md%BP%LL(idx(l)%ll(k))%matrices_block(idx(l)%bb(k))
         call assert(all(blk%M == A%m) .and. all(blk%N == A%n), 'htlu_build: block sizes do not match the clusters')
         if (blk%style == 1) then
            A%style = 1
            allocate (A%D(size(blk%fullmat, 1), size(blk%fullmat, 2)))
            A%D = blk%fullmat
         else
            call assert(blk%style == 2 .and. blk%level_butterfly == 0, 'htlu_build: only Tucker (lrlevel 0) admissible blocks')
            call assert(.not. allocated(blk%MiddleZFP) .and. .not. allocated(blk%MiddleQTT), 'htlu_build: ZFP/QTT cores')
            call assert(size(blk%ButterflyU(1)%blocks, 1) == 1 .and. size(blk%ButterflyV(1)%blocks, 1) == 1, 'htlu_build: Tucker factors')
            A%style = 2
            allocate (A%U(Ndim), A%V(Ndim))
            do d = 1, Ndim
               allocate (A%U(d)%mat(size(blk%ButterflyU(1)%blocks(1, d)%matrix, 1), size(blk%ButterflyU(1)%blocks(1, d)%matrix, 2)))
               A%U(d)%mat = blk%ButterflyU(1)%blocks(1, d)%matrix
               allocate (A%V(d)%mat(size(blk%ButterflyV(1)%blocks(1, d)%matrix, 1), size(blk%ButterflyV(1)%blocks(1, d)%matrix, 2)))
               A%V(d)%mat = blk%ButterflyV(1)%blocks(1, d)%matrix
            enddo
            allocate (A%C(size(blk%ButterflyMiddle(1)%matrix, 1), size(blk%ButterflyMiddle(1)%matrix, 2)))
            A%C = blk%ButterflyMiddle(1)%matrix
         endif
      else
         call assert(2*rg(1) <= msh(1)%Maxgroup, 'htlu_build: a leaf-level block is missing')
         A%style = 4
         A%nch = 2**Ndim
         allocate (A%sons(A%nch, A%nch))
         do j = 1, A%nch
            call htlu_child(Ndim, j, c2)
            cgc = 2*cg + c2 - 1
            do i = 1, A%nch
               call htlu_child(Ndim, i, c1)
               rgc = 2*rg + c1 - 1
               call htlu_build(A%sons(i, j), Ndim, l + 1, rgc, cgc, msh, h_mat_md, idx, Lmax)
            enddo
         enddo
      endif
   end subroutine htlu_build

   !>**** the block A between the clusters of groups rg and cg at level Ld from the local HTENSOR blocks indexed by idx: the
   !>**** restriction of a coarser Tucker block that contains it, or the copy of the blocks inside it (htlu_build)
   subroutine htlu_build_ld(A, Ndim, Ld, rg, cg, msh, h_mat_md, idx, Lmax)
      type(htlu_node) :: A
      integer Ndim, Ld, rg(Ndim), cg(Ndim), Lmax
      type(mesh) :: msh(Ndim)
      type(Hmat_md), target :: h_mat_md
      type(htlu_index_level) :: idx(0:)
      type(htlu_node) :: T
      integer l, d, k, ri, ci, rga(Ndim), cga(Ndim), ro(Ndim), rs(Ndim), co(Ndim), cs(Ndim)
      do l = 0, min(Ld - 1, Lmax)
         rga = rg/2**(Ld - l)
         cga = cg/2**(Ld - l)
         ri = htlu_cluster_index(Ndim, l, rga)
         ci = htlu_cluster_index(Ndim, l, cga)
         k = idx(l)%head(ri)
         do while (k /= 0)
            if (idx(l)%cidx(k) == ci) exit
            k = idx(l)%next(k)
         enddo
         if (k == 0) cycle
         call htlu_build(T, Ndim, l, rga, cga, msh, h_mat_md, idx, Lmax)
         call assert(T%style == 2, 'htlu_build_ld: a dense block above the distribution level')
         do d = 1, Ndim
            ro(d) = msh(d)%basis_group(rg(d))%head - msh(d)%basis_group(rga(d))%head
            co(d) = msh(d)%basis_group(cg(d))%head - msh(d)%basis_group(cga(d))%head
            rs(d) = msh(d)%basis_group(rg(d))%tail - msh(d)%basis_group(rg(d))%head + 1
            cs(d) = msh(d)%basis_group(cg(d))%tail - msh(d)%basis_group(cg(d))%head + 1
         enddo
         call htlu_restrict_box(T, ro, rs, co, cs, Ld, A)
         allocate (A%rg(Ndim), A%cg(Ndim))
         A%rg = rg
         A%cg = cg
         call HTLU_delete(T)
         return
      enddo
      call htlu_build(A, Ndim, Ld, rg, cg, msh, h_mat_md, idx, Lmax)
   end subroutine htlu_build_ld

   !>**** T = S, moving the storage (S becomes empty)
   subroutine htlu_move(S, T)
      type(htlu_node) :: S, T
      call assert(.not. associated(T%sons), 'htlu_move: the target has sons')
      T%style = S%style
      T%level = S%level
      T%Ndim = S%Ndim
      T%nch = S%nch
      if (allocated(S%m)) call move_alloc(S%m, T%m)
      if (allocated(S%n)) call move_alloc(S%n, T%n)
      if (allocated(S%rg)) call move_alloc(S%rg, T%rg)
      if (allocated(S%cg)) call move_alloc(S%cg, T%cg)
      if (allocated(S%D)) call move_alloc(S%D, T%D)
      if (allocated(S%ipiv)) call move_alloc(S%ipiv, T%ipiv)
      if (allocated(S%U)) call move_alloc(S%U, T%U)
      if (allocated(S%V)) call move_alloc(S%V, T%V)
      if (allocated(S%C)) call move_alloc(S%C, T%C)
      T%sons => S%sons
      nullify (S%sons)
      S%style = 0
   end subroutine htlu_move

   !>**** sizes of the integer and the data parts of the serialized block A (added to ni and nd)
   recursive subroutine htlu_pack_count(A, ni, nd)
      type(htlu_node) :: A
      integer(kind=8) ni, nd
      integer i, j, d
      call assert(A%style /= 0, 'htlu_pack_count: an empty block')
      ni = ni + 5 + 4*A%Ndim
      select case (A%style)
      case (1)
         ni = ni + 3
         if (allocated(A%ipiv)) ni = ni + size(A%ipiv)
         nd = nd + size(A%D)
      case (2)
         ni = ni + 4*A%Ndim + 2
         do d = 1, A%Ndim
            nd = nd + size(A%U(d)%mat) + size(A%V(d)%mat)
         enddo
         nd = nd + size(A%C)
      case (4)
         do j = 1, A%nch
            do i = 1, A%nch
               call htlu_pack_count(A%sons(i, j), ni, nd)
            enddo
         enddo
      end select
   end subroutine htlu_pack_count

   subroutine htlu_pack_mat(M, db, dp)
      DT :: M(:, :), db(:)
      integer dp, j, n1
      n1 = size(M, 1)
      do j = 1, size(M, 2)
         db(dp + 1:dp + n1) = M(:, j)
         dp = dp + n1
      enddo
   end subroutine htlu_pack_mat

   subroutine htlu_unpack_mat(M, n1, n2, db, dp)
      DT, allocatable :: M(:, :)
      DT :: db(:)
      integer n1, n2, dp, j
      allocate (M(n1, n2))
      do j = 1, n2
         M(:, j) = db(dp + 1:dp + n1)
         dp = dp + n1
      enddo
   end subroutine htlu_unpack_mat

   !>**** serialize the block A into ib(ip+1:) and db(dp+1:) (ip and dp advanced)
   recursive subroutine htlu_pack(A, ib, ip, db, dp)
      type(htlu_node) :: A
      integer ib(:), ip, dp
      DT :: db(:)
      integer i, j, d, Nd
      Nd = A%Ndim
      ib(ip + 1) = A%style
      ib(ip + 2) = A%level
      ib(ip + 3) = Nd
      ib(ip + 4) = A%nch
      ib(ip + 5) = 0
      if (allocated(A%rg)) ib(ip + 5) = 1
      ip = ip + 5
      ib(ip + 1:ip + Nd) = A%m
      ib(ip + Nd + 1:ip + 2*Nd) = A%n
      ib(ip + 2*Nd + 1:ip + 4*Nd) = 0
      if (allocated(A%rg)) then
         ib(ip + 2*Nd + 1:ip + 3*Nd) = A%rg
         ib(ip + 3*Nd + 1:ip + 4*Nd) = A%cg
      endif
      ip = ip + 4*Nd
      select case (A%style)
      case (1)
         ib(ip + 1) = size(A%D, 1)
         ib(ip + 2) = size(A%D, 2)
         ib(ip + 3) = 0
         if (allocated(A%ipiv)) ib(ip + 3) = size(A%ipiv)
         ip = ip + 3
         if (allocated(A%ipiv)) then
            ib(ip + 1:ip + size(A%ipiv)) = A%ipiv
            ip = ip + size(A%ipiv)
         endif
         call htlu_pack_mat(A%D, db, dp)
      case (2)
         do d = 1, Nd
            ib(ip + 1) = size(A%U(d)%mat, 1)
            ib(ip + 2) = size(A%U(d)%mat, 2)
            ib(ip + 3) = size(A%V(d)%mat, 1)
            ib(ip + 4) = size(A%V(d)%mat, 2)
            ip = ip + 4
         enddo
         ib(ip + 1) = size(A%C, 1)
         ib(ip + 2) = size(A%C, 2)
         ip = ip + 2
         do d = 1, Nd
            call htlu_pack_mat(A%U(d)%mat, db, dp)
            call htlu_pack_mat(A%V(d)%mat, db, dp)
         enddo
         call htlu_pack_mat(A%C, db, dp)
      case (4)
         do j = 1, A%nch
            do i = 1, A%nch
               call htlu_pack(A%sons(i, j), ib, ip, db, dp)
            enddo
         enddo
      end select
   end subroutine htlu_pack

   !>**** the block A from its serialization at ib(ip+1:) and db(dp+1:) (ip and dp advanced)
   recursive subroutine htlu_unpack(A, ib, ip, db, dp)
      type(htlu_node) :: A
      integer ib(:), ip, dp
      DT :: db(:)
      integer i, j, d, Nd, hasrg, n1, n2, np
      integer, allocatable :: su(:, :)
      A%style = ib(ip + 1)
      A%level = ib(ip + 2)
      Nd = ib(ip + 3)
      A%Ndim = Nd
      A%nch = ib(ip + 4)
      hasrg = ib(ip + 5)
      ip = ip + 5
      allocate (A%m(Nd), A%n(Nd))
      A%m = ib(ip + 1:ip + Nd)
      A%n = ib(ip + Nd + 1:ip + 2*Nd)
      if (hasrg == 1) then
         allocate (A%rg(Nd), A%cg(Nd))
         A%rg = ib(ip + 2*Nd + 1:ip + 3*Nd)
         A%cg = ib(ip + 3*Nd + 1:ip + 4*Nd)
      endif
      ip = ip + 4*Nd
      select case (A%style)
      case (1)
         n1 = ib(ip + 1)
         n2 = ib(ip + 2)
         np = ib(ip + 3)
         ip = ip + 3
         if (np > 0) then
            allocate (A%ipiv(np))
            A%ipiv = ib(ip + 1:ip + np)
            ip = ip + np
         endif
         call htlu_unpack_mat(A%D, n1, n2, db, dp)
      case (2)
         allocate (su(4, Nd), A%U(Nd), A%V(Nd))
         do d = 1, Nd
            su(:, d) = ib(ip + 1:ip + 4)
            ip = ip + 4
         enddo
         n1 = ib(ip + 1)
         n2 = ib(ip + 2)
         ip = ip + 2
         do d = 1, Nd
            call htlu_unpack_mat(A%U(d)%mat, su(1, d), su(2, d), db, dp)
            call htlu_unpack_mat(A%V(d)%mat, su(3, d), su(4, d), db, dp)
         enddo
         call htlu_unpack_mat(A%C, n1, n2, db, dp)
         deallocate (su)
      case (4)
         allocate (A%sons(A%nch, A%nch))
         do j = 1, A%nch
            do i = 1, A%nch
               call htlu_unpack(A%sons(i, j), ib, ip, db, dp)
            enddo
         enddo
      end select
   end subroutine htlu_unpack

   !>**** rank of the owner of block (a, b)
   integer function htlu_dist_owner(D, a, b)
      type(htlu_dist) :: D
      integer a, b
      htlu_dist_owner = mod(a - 1, D%nprow)*D%npcol + htlu_dist_pc(D, b)
   end function htlu_dist_owner

   !>**** process column of the blocks of column b
   integer function htlu_dist_pc(D, b)
      type(htlu_dist) :: D
      integer b
      htlu_dist_pc = mod(b - 1 + D%skew*((b - 1)/D%npcol), D%npcol)
   end function htlu_dist_pc

   !>**** global row (side 'r') or column (side 'c') of local block row or column i
   integer function htlu_dist_glob(D, side, i)
      type(htlu_dist) :: D
      character side
      integer i
      if (side == 'r') then
         htlu_dist_glob = D%myrow + 1 + (i - 1)*D%nprow
      else
         htlu_dist_glob = (i - 1)*D%npcol + mod(D%mycol - D%skew*(i - 1) + D%npcol*(i + 1), D%npcol) + 1
      endif
   end function htlu_dist_glob

   !>**** block (a, b): local, or fetched into recv (htlu_dist_fetch)
   function htlu_dist_ptr(D, recv, a, b) result(p)
      type(htlu_dist), target :: D
      type(htlu_node), target :: recv(:)
      integer a, b
      type(htlu_node), pointer :: p
      if (htlu_dist_owner(D, a, b) == D%myid) then
         p => D%blk((a - 1)/D%nprow + 1, (b - 1)/D%npcol + 1)
      else
         call assert(D%rmap(a, b) > 0, 'htlu_dist_ptr: a remote block was not fetched')
         p => recv(D%rmap(a, b))
      endif
   end function htlu_dist_ptr

   !>**** set up the distributed forest of P trees on the processes of ptree (the blocks are filled by HTLU_dist_add_pair):
   !>**** msh: the per-dimension cluster trees (the same for all trees); Ns_loc and idxs: the box of this process in the vector
   !>**** layout of the HTENSOR products; Ldreq: requested distribution level, raised to the coarsest level at which every cluster
   !>**** lies on one process in the HTENSOR layout (collective)
   subroutine HTLU_dist_init(D, Ndim, P, msh, ptree, Ns_loc, idxs, Ldreq)
      type(htlu_dist) :: D
      integer Ndim, P, Ns_loc(Ndim), idxs(Ndim), Ldreq
      type(mesh) :: msh(Ndim)
      type(proctree) :: ptree
      integer i, dd, l, m, t, dig, Lr, Maxlevel, g(Ndim), c(Ndim), ierr
      integer, allocatable :: mybox(:)
      logical ok
      D%Ndim = Ndim
      D%P = P
      D%comm = ptree%Comm
      D%nproc = ptree%nproc
      D%myid = ptree%MyID
      do i = int(sqrt(dble(D%nproc)) + 1d-8), 1, -1
         if (mod(D%nproc, i) == 0) exit
      enddo
      D%nprow = i
      D%npcol = D%nproc/i
      D%myrow = D%myid/D%npcol
      D%mycol = mod(D%myid, D%npcol)
      Maxlevel = GetTreelevel(msh(1)%Maxgroup) - 1
      Lr = 0
      do
         ok = .true.
         do t = 0, 2**(Lr*Ndim) - 1
            call htlu_cluster_groups(Ndim, Lr, t, g)
            if (ptree%pgrp(GetMshGroup_Pgno(ptree, Ndim, g))%nproc /= 1) ok = .false.
         enddo
         if (ok) exit
         Lr = Lr + 1
         call assert(Lr <= Maxlevel, 'HTLU_dist_init: more processes than leaf clusters')
      enddo
      D%Ld = max(Lr, Ldreq)
      call assert(D%Ld <= Maxlevel, 'HTLU_dist_init: the distribution level is below the leaf level')
      D%nch = 2**Ndim
      D%nc = D%nch**D%Ld
      D%K = P*D%nc
      allocate (D%Ns(Ndim), D%coff(Ndim, D%nc), D%csz(Ndim, D%nc), D%cgrp(Ndim, D%nc))
      do dd = 1, Ndim
         call assert(msh(dd)%basis_group(1)%head == 1, 'HTLU_dist_init: the cluster tree does not start at 1')
         D%Ns(dd) = msh(dd)%basis_group(1)%tail
      enddo
      D%nt = product(D%Ns)
      do m = 1, D%nc
         g = 1
         do l = 1, D%Ld
            dig = mod((m - 1)/D%nch**(D%Ld - l), D%nch) + 1
            call htlu_child(Ndim, dig, c)
            g = 2*g + c - 1
         enddo
         D%cgrp(:, m) = g
         do dd = 1, Ndim
            D%coff(dd, m) = msh(dd)%basis_group(g(dd))%head - 1
            D%csz(dd, m) = msh(dd)%basis_group(g(dd))%tail - msh(dd)%basis_group(g(dd))%head + 1
         enddo
      enddo
      allocate (D%boxes(2*Ndim, 0:D%nproc - 1), mybox(2*Ndim))
      mybox(1:Ndim) = idxs - 1
      mybox(Ndim + 1:2*Ndim) = Ns_loc
      call MPI_Allgather(mybox, 2*Ndim, MPI_INTEGER, D%boxes, 2*Ndim, MPI_INTEGER, D%comm, ierr)
      deallocate (mybox)
      D%nr = 0
      do i = 1, D%K
         if (mod(i - 1, D%nprow) == D%myrow) D%nr = D%nr + 1
      enddo
      D%skew = 0
      if (mod(D%K, D%npcol) == 0) D%skew = 1
      D%nl = 0
      do i = 1, D%K
         if (htlu_dist_pc(D, i) == D%mycol) D%nl = D%nl + 1
      enddo
      allocate (D%blk(D%nr, D%nl), D%rmap(D%K, D%K))
      D%rmap = 0
   end subroutine HTLU_dist_init

   !>**** fill the blocks between tree p (rows) and tree q (columns) from the HTENSOR h_mat_md of the pair (collective): every
   !>**** process builds the level-Ld blocks of its rows in the HTENSOR layout and sends them to their owners
   subroutine HTLU_dist_add_pair(D, p, q, h_mat_md, msh, ptree)
      type(htlu_dist), target :: D
      integer p, q
      type(Hmat_md), target :: h_mat_md
      type(mesh) :: msh(:)
      type(proctree) :: ptree
      type(htlu_index_level), allocatable :: idx(:)
      type(htlu_node), allocatable, target :: tmp(:)
      type(htlu_node), allocatable :: rnd(:)
      type(htlu_nptr), allocatable :: sptr(:)
      integer, allocatable :: sdest(:), sid(:, :), rid(:, :)
      integer Lmax, m, n, a, b, ns, nrcv, k, fo, ow
      call htlu_index_build(h_mat_md, D%Ndim, idx, Lmax)
      allocate (tmp(D%nc*D%nc), sptr(D%nc*D%nc), sdest(D%nc*D%nc), sid(2, D%nc*D%nc))
      ns = 0
      do m = 1, D%nc
         fo = ptree%pgrp(GetMshGroup_Pgno(ptree, D%Ndim, D%cgrp(:, m)))%head
         if (fo /= D%myid) cycle
         a = (p - 1)*D%nc + m
         do n = 1, D%nc
            b = (q - 1)*D%nc + n
            ow = htlu_dist_owner(D, a, b)
            if (ow == D%myid) then
               call htlu_build_ld(D%blk((a - 1)/D%nprow + 1, (b - 1)/D%npcol + 1), D%Ndim, D%Ld, D%cgrp(:, m), D%cgrp(:, n), &
                  msh, h_mat_md, idx, Lmax)
            else
               ns = ns + 1
               call htlu_build_ld(tmp(ns), D%Ndim, D%Ld, D%cgrp(:, m), D%cgrp(:, n), msh, h_mat_md, idx, Lmax)
               sptr(ns)%p => tmp(ns)
               sdest(ns) = ow
               sid(:, ns) = (/a, b/)
            endif
         enddo
      enddo
      call htlu_dist_xchg(D, ns, sdest, sid, sptr, nrcv, rid, rnd)
      do k = 1, ns
         call HTLU_delete(tmp(k))
      enddo
      do k = 1, nrcv
         call htlu_move(rnd(k), D%blk((rid(1, k) - 1)/D%nprow + 1, (rid(2, k) - 1)/D%npcol + 1))
      enddo
      deallocate (tmp, sptr, sdest, sid, rid, rnd)
      call htlu_index_free(idx)
   end subroutine HTLU_dist_add_pair

   !>**** send the blocks sptr(i) (positions sid(:, i)) to the ranks sdest(i), i = 1..ns, and receive nrcv blocks rnd with
   !>**** positions rid (collective)
   subroutine htlu_dist_xchg(D, ns, sdest, sid, sptr, nrcv, rid, rnd)
      type(htlu_dist) :: D
      integer ns, sdest(:), sid(:, :), nrcv
      type(htlu_nptr) :: sptr(:)
      integer, allocatable :: rid(:, :)
      type(htlu_node), allocatable :: rnd(:)
      integer, allocatable :: scnt(:, :), rcnt(:, :), sdi(:), sdd(:), rdi(:), rdd(:), ip(:), dp(:), ci(:), cd(:), rci(:), rcd(:)
      integer, allocatable :: ibuf(:), ribuf(:)
      DT, allocatable :: dbuf(:), rdbuf(:)
      integer(kind=8) ni, nd, tsi, tsd, tri, trd
      integer i, r, k, ierr, np, ipos, dpos
      np = D%nproc
      allocate (scnt(3, 0:np - 1), rcnt(3, 0:np - 1), sdi(0:np - 1), sdd(0:np - 1), rdi(0:np - 1), rdd(0:np - 1), ip(0:np - 1), &
         dp(0:np - 1), ci(0:np - 1), cd(0:np - 1), rci(0:np - 1), rcd(0:np - 1))
      scnt = 0
      do i = 1, ns
         ni = 2
         nd = 0
         call htlu_pack_count(sptr(i)%p, ni, nd)
         r = sdest(i)
         call assert(dble(scnt(1, r)) + dble(ni) < dble(huge(1)) .and. dble(scnt(2, r)) + dble(nd) < dble(huge(1)), &
            'htlu_dist_xchg: message too large')
         scnt(1, r) = scnt(1, r) + int(ni)
         scnt(2, r) = scnt(2, r) + int(nd)
         scnt(3, r) = scnt(3, r) + 1
      enddo
      call MPI_Alltoall(scnt, 3, MPI_INTEGER, rcnt, 3, MPI_INTEGER, D%comm, ierr)
      tsi = 0
      tsd = 0
      tri = 0
      trd = 0
      do r = 0, np - 1
         sdi(r) = int(tsi)
         sdd(r) = int(tsd)
         rdi(r) = int(tri)
         rdd(r) = int(trd)
         tsi = tsi + scnt(1, r)
         tsd = tsd + scnt(2, r)
         tri = tri + rcnt(1, r)
         trd = trd + rcnt(2, r)
      enddo
      call assert(max(tsi, tsd, tri, trd) < huge(1), 'htlu_dist_xchg: messages too large')
      allocate (ibuf(max(1_8, tsi)), dbuf(max(1_8, tsd)), ribuf(max(1_8, tri)), rdbuf(max(1_8, trd)))
      ip = sdi
      dp = sdd
      do i = 1, ns
         r = sdest(i)
         ibuf(ip(r) + 1) = sid(1, i)
         ibuf(ip(r) + 2) = sid(2, i)
         ip(r) = ip(r) + 2
         call htlu_pack(sptr(i)%p, ibuf, ip(r), dbuf, dp(r))
      enddo
      ci = scnt(1, :)
      cd = scnt(2, :)
      rci = rcnt(1, :)
      rcd = rcnt(2, :)
      call MPI_Alltoallv(ibuf, ci, sdi, MPI_INTEGER, ribuf, rci, rdi, MPI_INTEGER, D%comm, ierr)
      call MPI_Alltoallv(dbuf, cd, sdd, MPI_DT, rdbuf, rcd, rdd, MPI_DT, D%comm, ierr)
      deallocate (ibuf, dbuf)
      nrcv = sum(rcnt(3, :))
      allocate (rid(2, max(1, nrcv)), rnd(nrcv))
      ipos = 0
      dpos = 0
      do k = 1, nrcv
         rid(:, k) = ribuf(ipos + 1:ipos + 2)
         ipos = ipos + 2
         call htlu_unpack(rnd(k), ribuf, ipos, rdbuf, dpos)
      enddo
      D%vrecv = D%vrecv + dble(trd)*dble(storage_size(BPACK_cone))/8d0/1024d0/1024d0
      deallocate (ribuf, rdbuf, scnt, rcnt, sdi, sdd, rdi, rdd, ip, dp, ci, cd, rci, rcd)
   end subroutine htlu_dist_xchg

   !>**** mark block (a, b) as needed (once)
   subroutine htlu_dist_addneed(D, a, b, nneed, need)
      type(htlu_dist) :: D
      integer a, b, nneed, need(:, :)
      if (D%rmap(a, b) /= 0) return
      D%rmap(a, b) = -1
      nneed = nneed + 1
      need(1, nneed) = a
      need(2, nneed) = b
   end subroutine htlu_dist_addneed

   !>**** fetch the remote blocks at the positions need(:, 1:nneed) into recv, with D%rmap(a, b) its index in recv (collective)
   subroutine htlu_dist_fetch(D, nneed, need, recv)
      type(htlu_dist), target :: D
      integer nneed, need(:, :)
      type(htlu_node), allocatable :: recv(:)
      integer, allocatable :: scnt(:), rcnt(:), sdsp(:), rdsp(:), pos(:), sreq(:), rreq(:), sdest(:), sid(:, :), rid(:, :)
      type(htlu_nptr), allocatable :: sptr(:)
      integer i, k, r, a, b, ns, nrcv, ierr
      real(kind=8) t0
      t0 = MPI_Wtime()
      allocate (scnt(0:D%nproc - 1), rcnt(0:D%nproc - 1), sdsp(0:D%nproc - 1), rdsp(0:D%nproc - 1), pos(0:D%nproc - 1))
      scnt = 0
      do i = 1, nneed
         r = htlu_dist_owner(D, need(1, i), need(2, i))
         scnt(r) = scnt(r) + 2
      enddo
      call MPI_Alltoall(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, D%comm, ierr)
      sdsp(0) = 0
      rdsp(0) = 0
      do r = 1, D%nproc - 1
         sdsp(r) = sdsp(r - 1) + scnt(r - 1)
         rdsp(r) = rdsp(r - 1) + rcnt(r - 1)
      enddo
      allocate (sreq(max(1, sum(scnt))), rreq(max(1, sum(rcnt))))
      pos = sdsp
      do i = 1, nneed
         r = htlu_dist_owner(D, need(1, i), need(2, i))
         sreq(pos(r) + 1:pos(r) + 2) = need(1:2, i)
         pos(r) = pos(r) + 2
      enddo
      call MPI_Alltoallv(sreq, scnt, sdsp, MPI_INTEGER, rreq, rcnt, rdsp, MPI_INTEGER, D%comm, ierr)
      ns = sum(rcnt)/2
      allocate (sdest(max(1, ns)), sid(2, max(1, ns)), sptr(max(1, ns)))
      i = 0
      do r = 0, D%nproc - 1
         do k = 1, rcnt(r)/2
            i = i + 1
            a = rreq(rdsp(r) + 2*k - 1)
            b = rreq(rdsp(r) + 2*k)
            call assert(htlu_dist_owner(D, a, b) == D%myid, 'htlu_dist_fetch: a request for a block of another process')
            sdest(i) = r
            sid(:, i) = (/a, b/)
            sptr(i)%p => D%blk((a - 1)/D%nprow + 1, (b - 1)/D%npcol + 1)
         enddo
      enddo
      call htlu_dist_xchg(D, ns, sdest, sid, sptr, nrcv, rid, recv)
      call assert(nrcv == nneed, 'htlu_dist_fetch: wrong number of blocks received')
      do i = 1, nrcv
         D%rmap(rid(1, i), rid(2, i)) = i
      enddo
      deallocate (scnt, rcnt, sdsp, rdsp, pos, sreq, rreq, sdest, sid, sptr, rid)
      D%tcomm = D%tcomm + MPI_Wtime() - t0
      D%nstep = D%nstep + 1
   end subroutine htlu_dist_fetch

   !>**** free the fetched blocks
   subroutine htlu_dist_release(D, nneed, need, recv)
      type(htlu_dist) :: D
      integer nneed, need(:, :), i
      type(htlu_node), allocatable :: recv(:)
      do i = 1, nneed
         D%rmap(need(1, i), need(2, i)) = 0
      enddo
      do i = 1, size(recv)
         call HTLU_delete(recv(i))
      enddo
      deallocate (recv)
   end subroutine htlu_dist_release

   !>**** lo..hi
   function htlu_irange(lo, hi) result(v)
      integer lo, hi, i
      integer, allocatable :: v(:)
      allocate (v(max(0, hi - lo + 1)))
      do i = lo, hi
         v(i - lo + 1) = i
      enddo
   end function htlu_irange

   !>**** lo..hi without xlo..xhi
   function htlu_irange_ex(lo, hi, xlo, xhi) result(v)
      integer lo, hi, xlo, xhi, i, k
      integer, allocatable :: v(:)
      allocate (v(max(0, hi - lo + 1 - max(0, min(hi, xhi) - max(lo, xlo) + 1))))
      k = 0
      do i = lo, hi
         if (i >= xlo .and. i <= xhi) cycle
         k = k + 1
         v(k) = i
      enddo
   end function htlu_irange_ex

   !>**** the sons of the cluster range lo..hi: the trees for the whole forest, the 2^Ndim children inside a tree
   subroutine htlu_dist_sons(D, lo, hi, ns, slo, shi)
      type(htlu_dist) :: D
      integer lo, hi, ns, i, sz
      integer, allocatable :: slo(:), shi(:)
      if (hi - lo + 1 > D%nc) then
         ns = D%P
         sz = D%nc
      else
         ns = D%nch
         sz = (hi - lo + 1)/D%nch
      endif
      allocate (slo(ns), shi(ns))
      do i = 1, ns
         slo(i) = lo + (i - 1)*sz
         shi(i) = slo(i) + sz - 1
      enddo
   end subroutine htlu_dist_sons

   !>**** start a step of the profile (collective)
   subroutine htlu_dist_pbeg(D)
      type(htlu_dist) :: D
      integer ierr
      call MPI_Barrier(D%comm, ierr)
      D%pt0 = MPI_Wtime()
      D%pc0 = D%tcomp
   end subroutine htlu_dist_pbeg

   !>**** end a step of kind typ summing c clusters with ntg local targets (collective)
   subroutine htlu_dist_pend(D, typ, c, ntg)
      type(htlu_dist) :: D
      integer typ, c, ntg, cls, t, ierr
      call MPI_Barrier(D%comm, ierr)
      cls = 0
      t = 1
      do while (t < c .and. cls < 9)
         t = t*D%nch
         cls = cls + 1
      enddo
      D%pwall(typ, cls) = D%pwall(typ, cls) + MPI_Wtime() - D%pt0
      D%pwork(typ, cls) = D%pwork(typ, cls) + D%tcomp - D%pc0
      D%pcnt(typ, cls) = D%pcnt(typ, cls) + 1
      D%pntg(typ, cls) = D%pntg(typ, cls) + ntg
   end subroutine htlu_dist_pend

   !>**** the targets (a, b), a in rows and b in cols (rows fastest), done by this process: every target by its owner, or round
   !>**** robin over all processes when the owners of the targets are fewer than the processes and than the targets (a step on
   !>**** a few block rows or columns is then spread over all processes, with the targets moved to them and back)
   subroutine htlu_dist_myops(D, rows, cols, nmy, ta, tb)
      type(htlu_dist) :: D
      integer rows(:), cols(:), nmy, n, nown, i, j, t, ex
      integer, allocatable :: ta(:), tb(:)
      logical, allocatable :: used(:)
      logical spread
      n = size(rows)*size(cols)
      allocate (used(0:D%nproc - 1))
      used = .false.
      do j = 1, size(cols)
         do i = 1, size(rows)
            used(htlu_dist_owner(D, rows(i), cols(j))) = .true.
         enddo
      enddo
      nown = count(used)
      spread = nown < D%nproc .and. n > nown
      allocate (ta(max(1, n)), tb(max(1, n)))
      nmy = 0
      t = 0
      do j = 1, size(cols)
         do i = 1, size(rows)
            t = t + 1
            if (spread) then
               ex = mod(t - 1 + D%rr, D%nproc)
            else
               ex = htlu_dist_owner(D, rows(i), cols(j))
            endif
            if (ex == D%myid) then
               nmy = nmy + 1
               ta(nmy) = rows(i)
               tb(nmy) = cols(j)
            endif
         enddo
      enddo
      if (spread) D%rr = mod(D%rr + n, D%nproc)
      deallocate (used)
   end subroutine htlu_dist_myops

   !>**** send the targets updated here for other owners (fetched into recv, or the new blocks Y) back to them, and replace the
   !>**** local targets updated elsewhere or by new blocks (collective)
   subroutine htlu_dist_return(D, nmy, ta, tb, recv, Y, newblk)
      type(htlu_dist), target :: D
      integer nmy, ta(:), tb(:)
      type(htlu_node), target :: recv(:), Y(:)
      logical newblk
      type(htlu_node), allocatable :: rnd(:)
      type(htlu_nptr), allocatable :: sptr(:)
      type(htlu_node), pointer :: Tp
      integer, allocatable :: sdest(:), sid(:, :), rid(:, :)
      integer t, ns, nrcv, k, ow
      real(kind=8) t0
      t0 = MPI_Wtime()
      allocate (sdest(max(1, nmy)), sid(2, max(1, nmy)), sptr(max(1, nmy)))
      ns = 0
      do t = 1, nmy
         ow = htlu_dist_owner(D, ta(t), tb(t))
         Tp => D%blk((ta(t) - 1)/D%nprow + 1, (tb(t) - 1)/D%npcol + 1)
         if (ow == D%myid) then
            if (newblk) then
               call HTLU_delete(Tp)
               call htlu_move(Y(t), Tp)
            endif
         else
            ns = ns + 1
            sdest(ns) = ow
            sid(:, ns) = (/ta(t), tb(t)/)
            if (newblk) then
               sptr(ns)%p => Y(t)
            else
               sptr(ns)%p => recv(D%rmap(ta(t), tb(t)))
            endif
         endif
      enddo
      call htlu_dist_xchg(D, ns, sdest, sid, sptr, nrcv, rid, rnd)
      do k = 1, nrcv
         Tp => D%blk((rid(1, k) - 1)/D%nprow + 1, (rid(2, k) - 1)/D%npcol + 1)
         call HTLU_delete(Tp)
         call htlu_move(rnd(k), Tp)
      enddo
      deallocate (sdest, sid, sptr, rid, rnd)
      D%tcomm = D%tcomm + MPI_Wtime() - t0
   end subroutine htlu_dist_return

   !>**** G(a, b) = G(a, b) - sum_{k = k1..k2} G(a, k) G(k, b) for a in rows and b in cols, G the distributed blocks (newblk: G(a, b) =
   !>**** sum_k G(a, k) G(k, b), accumulated in new blocks that replace the targets at the end, so the sources may include the
   !>**** targets). The targets are updated by their owners or spread (htlu_dist_myops), after the remote operands are fetched
   !>**** (collective)
   subroutine htlu_dist_gemm(D, rows, cols, k1, k2, newblk, tol)
      type(htlu_dist), target :: D
      integer rows(:), cols(:), k1, k2
      logical newblk
      real(kind=8) tol(2)
      integer, allocatable :: ta(:), tb(:), need(:, :)
      type(htlu_node), allocatable, target :: recv(:), Y(:)
      type(htlu_node), pointer :: Ap, Bp, Tp
      integer nneed, k, t, a, b, nmy
      real(kind=8) t0
      call htlu_dist_pbeg(D)
      call htlu_dist_myops(D, rows, cols, nmy, ta, tb)
      nneed = 0
      allocate (need(2, max(1, nmy*(1 + 2*(k2 - k1 + 1)))))
      do t = 1, nmy
         a = ta(t)
         b = tb(t)
         if (htlu_dist_owner(D, a, b) /= D%myid) call htlu_dist_addneed(D, a, b, nneed, need)
         do k = k1, k2
            if (htlu_dist_owner(D, a, k) /= D%myid) call htlu_dist_addneed(D, a, k, nneed, need)
            if (htlu_dist_owner(D, k, b) /= D%myid) call htlu_dist_addneed(D, k, b, nneed, need)
         enddo
      enddo
      call htlu_dist_fetch(D, nneed, need, recv)
      t0 = MPI_Wtime()
      allocate (Y(nmy))
      if (nmy > 0) then
         !$omp parallel
         !$omp single
         do t = 1, nmy
            !$omp task default(shared) firstprivate(t) private(a, b, k, Ap, Bp, Tp)
            a = ta(t)
            b = tb(t)
            Tp => htlu_dist_ptr(D, recv, a, b)
            if (newblk) then
               call htlu_zero_clone(Tp, Y(t))
               Tp => Y(t)
            endif
            do k = k1, k2
               Ap => htlu_dist_ptr(D, recv, a, k)
               Bp => htlu_dist_ptr(D, recv, k, b)
               call htlu_subhh(Tp, Ap, Bp, tol)
            enddo
            if (newblk) call htlu_scale(Tp, -BPACK_cone)
            !$omp end task
         enddo
         !$omp end single
         !$omp end parallel
      endif
      D%tcomp = D%tcomp + MPI_Wtime() - t0
      call htlu_dist_return(D, nmy, ta, tb, recv, Y, newblk)
      call htlu_dist_release(D, nneed, need, recv)
      deallocate (need, ta, tb, Y)
      if (newblk) then
         call htlu_dist_pend(D, 7, k2 - k1 + 1, nmy)
      else
         call htlu_dist_pend(D, 6, k2 - k1 + 1, nmy)
      endif
   end subroutine htlu_dist_gemm

   !>**** the operation with the diagonal block S = G(d, d) on G(a, b) for a in rows and b in cols: op 1: G = G S^{-1} with the
   !>**** block-LU factored S (htlu_XSM), 2: G = G S, 3: G = S G (htlu_inv_mult); targets as in htlu_dist_gemm (collective)
   subroutine htlu_dist_diag(D, rows, cols, dd, op, tol)
      type(htlu_dist), target :: D
      integer rows(:), cols(:), dd, op
      real(kind=8) tol(2)
      integer, allocatable :: ta(:), tb(:), need(:, :)
      type(htlu_node), allocatable, target :: recv(:), Y(:)
      type(htlu_node), pointer :: Sp, Tp
      integer nneed, t, a, b, nmy
      real(kind=8) t0
      call htlu_dist_pbeg(D)
      call htlu_dist_myops(D, rows, cols, nmy, ta, tb)
      nneed = 0
      allocate (need(2, max(1, nmy + 1)))
      if (nmy > 0 .and. htlu_dist_owner(D, dd, dd) /= D%myid) call htlu_dist_addneed(D, dd, dd, nneed, need)
      do t = 1, nmy
         if (htlu_dist_owner(D, ta(t), tb(t)) /= D%myid) call htlu_dist_addneed(D, ta(t), tb(t), nneed, need)
      enddo
      call htlu_dist_fetch(D, nneed, need, recv)
      t0 = MPI_Wtime()
      if (nmy > 0) then
         !$omp parallel
         !$omp single
         do t = 1, nmy
            !$omp task default(shared) firstprivate(t) private(a, b, Sp, Tp)
            a = ta(t)
            b = tb(t)
            Tp => htlu_dist_ptr(D, recv, a, b)
            Sp => htlu_dist_ptr(D, recv, dd, dd)
            select case (op)
            case (1)
               call htlu_XSM(Sp, Tp, tol)
            case (2)
               call htlu_inv_mult(Tp, Sp, 'R', tol)
            case (3)
               call htlu_inv_mult(Tp, Sp, 'L', tol)
            end select
            !$omp end task
         enddo
         !$omp end single
         !$omp end parallel
      endif
      D%tcomp = D%tcomp + MPI_Wtime() - t0
      allocate (Y(0))
      call htlu_dist_return(D, nmy, ta, tb, recv, Y, .false.)
      call htlu_dist_release(D, nneed, need, recv)
      deallocate (need, ta, tb, Y)
      call htlu_dist_pend(D, 2 + op, 1, nmy)
   end subroutine htlu_dist_diag

   !>**** factor (op 1, htlu_BLU) or invert (op 2, htlu_inv) the diagonal block G(d, d) on its owner
   subroutine htlu_dist_local(D, dd, op, tol)
      type(htlu_dist), target :: D
      integer dd, op
      real(kind=8) tol(2), t0
      integer ntg
      call htlu_dist_pbeg(D)
      ntg = 0
      if (htlu_dist_owner(D, dd, dd) == D%myid) then
         ntg = 1
         t0 = MPI_Wtime()
         !$omp parallel
         !$omp single
         if (op == 1) then
            call htlu_BLU(D%blk((dd - 1)/D%nprow + 1, (dd - 1)/D%npcol + 1), tol)
         else
            call htlu_inv(D%blk((dd - 1)/D%nprow + 1, (dd - 1)/D%npcol + 1), tol)
         endif
         !$omp end single
         !$omp end parallel
         D%tcomp = D%tcomp + MPI_Wtime() - t0
      endif
      call htlu_dist_pend(D, op, 1, ntg)
   end subroutine htlu_dist_local

   !>**** G(a, b) = -G(a, b) for a and b in lo..hi
   subroutine htlu_dist_scale(D, lo, hi)
      type(htlu_dist) :: D
      integer lo, hi, i, j, a, b
      do j = 1, D%nl
         b = htlu_dist_glob(D, 'c', j)
         if (b < lo .or. b > hi) cycle
         do i = 1, D%nr
            a = htlu_dist_glob(D, 'r', i)
            if (a < lo .or. a > hi) cycle
            call htlu_scale(D%blk(i, j), -BPACK_cone)
         enddo
      enddo
   end subroutine htlu_dist_scale

   !>**** block LU with whole-cluster pivots (htlu_BLU) of the diagonal block on the clusters lo..hi
   recursive subroutine htlu_dist_blu(D, lo, hi, tol)
      type(htlu_dist), target :: D
      integer lo, hi, i, ns
      integer, allocatable :: slo(:), shi(:)
      real(kind=8) tol(2)
      if (lo == hi) then
         call htlu_dist_local(D, lo, 1, tol)
         return
      endif
      call htlu_dist_sons(D, lo, hi, ns, slo, shi)
      do i = 1, ns
         call htlu_dist_blu(D, slo(i), shi(i), tol)
         if (i < ns) then
            call htlu_dist_xsm(D, slo(i), shi(i), shi(i) + 1, hi, tol)
            call htlu_dist_gemm(D, htlu_irange(shi(i) + 1, hi), htlu_irange(shi(i) + 1, hi), slo(i), shi(i), .false., tol)
         endif
      enddo
      deallocate (slo, shi)
   end subroutine htlu_dist_blu

   !>**** M = M S^{-1} (htlu_XSM) for the blocks M = G(mlo..mhi, dlo..dhi) and the factored diagonal block S on dlo..dhi
   recursive subroutine htlu_dist_xsm(D, dlo, dhi, mlo, mhi, tol)
      type(htlu_dist), target :: D
      integer dlo, dhi, mlo, mhi, b, ns
      integer, allocatable :: slo(:), shi(:)
      real(kind=8) tol(2)
      if (dlo == dhi) then
         call htlu_dist_diag(D, htlu_irange(mlo, mhi), (/dlo/), dlo, 1, tol)
         return
      endif
      call htlu_dist_sons(D, dlo, dhi, ns, slo, shi)
      ! M Ub^{-1}, then the product with Lb^{-1}
      do b = 1, ns
         if (b > 1) call htlu_dist_gemm(D, htlu_irange(mlo, mhi), htlu_irange(slo(b), shi(b)), dlo, slo(b) - 1, .false., tol)
         call htlu_dist_xsm(D, slo(b), shi(b), mlo, mhi, tol)
      enddo
      do b = ns - 1, 1, -1
         call htlu_dist_gemm(D, htlu_irange(mlo, mhi), htlu_irange(slo(b), shi(b)), shi(b) + 1, dhi, .false., tol)
      enddo
      deallocate (slo, shi)
   end subroutine htlu_dist_xsm

   !>**** block Gauss-Jordan sweep (htlu_sweep) of the diagonal block on the clusters lo..hi: it becomes minus its inverse
   recursive subroutine htlu_dist_sweep(D, lo, hi, tol)
      type(htlu_dist), target :: D
      integer lo, hi, i, ns
      integer, allocatable :: slo(:), shi(:), ex(:), si(:)
      real(kind=8) tol(2)
      call htlu_dist_sons(D, lo, hi, ns, slo, shi)
      do i = 1, ns
         call htlu_dist_inv(D, slo(i), shi(i), tol)
         ex = htlu_irange_ex(lo, hi, slo(i), shi(i))
         si = htlu_irange(slo(i), shi(i))
         if (slo(i) == shi(i)) then
            call htlu_dist_diag(D, ex, si, slo(i), 2, tol)
         else
            call htlu_dist_gemm(D, ex, si, slo(i), shi(i), .true., tol)
         endif
         call htlu_dist_gemm(D, ex, ex, slo(i), shi(i), .false., tol)
         if (slo(i) == shi(i)) then
            call htlu_dist_diag(D, si, ex, slo(i), 3, tol)
         else
            call htlu_dist_gemm(D, si, ex, slo(i), shi(i), .true., tol)
         endif
         call htlu_dist_scale(D, slo(i), shi(i))
      enddo
      deallocate (slo, shi)
   end subroutine htlu_dist_sweep

   !>**** the diagonal block on the clusters lo..hi becomes its inverse (htlu_inv)
   recursive subroutine htlu_dist_inv(D, lo, hi, tol)
      type(htlu_dist), target :: D
      integer lo, hi
      real(kind=8) tol(2)
      if (lo == hi) then
         call htlu_dist_local(D, lo, 2, tol)
         return
      endif
      call htlu_dist_sweep(D, lo, hi, tol)
      call htlu_dist_scale(D, lo, hi)
   end subroutine htlu_dist_inv

   !>**** sum the leaf-operation timers and the exchange statistics over the processes and print them (rank 0)
   subroutine htlu_dist_timers(D, wall, nth)
      type(htlu_dist) :: D
      real(kind=8) wall, ts(7), tr(0:15, 2), v(3), vmax(3), vsum(3), pw(7, 0:9), pwx(7, 0:9)
      integer(kind=8) tc(7), nr(0:15, 2)
      integer nth, ierr, pn(7, 0:9), k, c
      character(len=6), parameter :: kname(7) = (/'factor', 'invert', 'solve ', 'mul R ', 'mul L ', 'Schur ', 'mul   '/)
      call MPI_Allreduce(htlu_tsec, ts, 7, MPI_DOUBLE_PRECISION, MPI_SUM, D%comm, ierr)
      call MPI_Allreduce(htlu_tcnt, tc, 7, MPI_INTEGER8, MPI_SUM, D%comm, ierr)
      call MPI_Allreduce(htlu_rlev_t, tr, 32, MPI_DOUBLE_PRECISION, MPI_SUM, D%comm, ierr)
      call MPI_Allreduce(htlu_rlev_n, nr, 32, MPI_INTEGER8, MPI_SUM, D%comm, ierr)
      v = (/D%tcomm, D%tcomp, D%vrecv/)
      call MPI_Allreduce(v, vmax, 3, MPI_DOUBLE_PRECISION, MPI_MAX, D%comm, ierr)
      call MPI_Allreduce(v, vsum, 3, MPI_DOUBLE_PRECISION, MPI_SUM, D%comm, ierr)
      call MPI_Allreduce(D%pwork, pw, 70, MPI_DOUBLE_PRECISION, MPI_SUM, D%comm, ierr)
      call MPI_Allreduce(D%pwork, pwx, 70, MPI_DOUBLE_PRECISION, MPI_MAX, D%comm, ierr)
      call MPI_Allreduce(D%pntg, pn, 70, MPI_INTEGER, MPI_SUM, D%comm, ierr)
      if (D%myid == 0) then
         htlu_tsec = ts
         htlu_tcnt = tc
         htlu_rlev_t = tr
         htlu_rlev_n = nr
         call HTLU_print_timers(wall, nth*D%nproc)
         write (*, '(A,I6,A,I3,A,I4,A,I4,A,I8,A)') ' tensor LU distribution:', D%K, ' clusters (level', D%Ld, ') on', D%nprow, ' x', &
            D%npcol, ' processes,', D%nstep, ' exchange steps'
         write (*, '(A,2F12.3,A,2F12.3,A,F12.1,F14.1)') '    exchanges (s, max mean)', vmax(1), vsum(1)/D%nproc, &
            '   block operations (s, max mean)', vmax(2), vsum(2)/D%nproc, '   received MB (max total)', vmax(3), vsum(3)
         write (*, '(A)') '    steps: kind, summed clusters <=, steps, targets, wall (s), block operations (s, sum and max over processes), efficiency'
         do k = 1, 7
            do c = 0, 9
               if (D%pcnt(k, c) == 0) cycle
               write (*, '(6X,A6,I8,I8,I10,F10.2,2F12.2,F8.1,A)') kname(k), D%nch**c, D%pcnt(k, c), pn(k, c), D%pwall(k, c), pw(k, c), &
                  pwx(k, c), 100d0*pw(k, c)/max(1d-30, D%pwall(k, c)*D%nproc), ' %'
            enddo
         enddo
      else
         htlu_tsec = 0
         htlu_tcnt = 0
         htlu_rlev_t = 0
         htlu_rlev_n = 0
      endif
      D%tcomm = 0
      D%tcomp = 0
      D%vrecv = 0
      D%nstep = 0
      D%pwall = 0
      D%pwork = 0
      D%pcnt = 0
      D%pntg = 0
   end subroutine htlu_dist_timers

   !>**** block LU factorization with whole-cluster pivots of the distributed forest in place (see htlu_BLU and htlu_bsolve; the
   !>**** same block LU over the trees as inside each tree) (normest: estimate of ||A||_2; Tucker truncation at option%tol_rand
   !>**** relative to each block and to normest) (collective)
   subroutine HTLU_dist_factor(D, option, stats, normest)
      type(htlu_dist), target :: D
      type(Hoption) :: option
      type(Hstat) :: stats
      real(kind=8) normest, t1, t2, tol(2)
      integer nth, ierr
      tol = (/option%tol_rand, option%tol_rand*normest/)
      nth = 1
!$    nth = omp_get_max_threads()
      call MPI_Barrier(D%comm, ierr)
      t1 = MPI_Wtime()
      call htlu_dist_blu(D, 1, D%K, tol)
      call MPI_Barrier(D%comm, ierr)
      t2 = MPI_Wtime()
      if (option%verbosity >= 0) call htlu_dist_timers(D, t2 - t1, nth)
      stats%Time_Factor = stats%Time_Factor + t2 - t1
   end subroutine HTLU_dist_factor

   !>**** the distributed forest becomes its inverse (H-TENSOR inverse with Tucker blocks, block Gauss-Jordan; apply it with
   !>**** HTLU_dist_mult) (normest: estimate of ||A^{-1}||_2) (collective)
   subroutine HTLU_dist_inverse(D, option, stats, normest)
      type(htlu_dist), target :: D
      type(Hoption) :: option
      type(Hstat) :: stats
      real(kind=8) normest, t1, t2, tol(2)
      integer nth, ierr
      tol = (/option%tol_rand, option%tol_rand*normest/)
      nth = 1
!$    nth = omp_get_max_threads()
      call MPI_Barrier(D%comm, ierr)
      t1 = MPI_Wtime()
      call htlu_dist_inv(D, 1, D%K, tol)
      call MPI_Barrier(D%comm, ierr)
      t2 = MPI_Wtime()
      D%inverse = 1
      if (option%verbosity >= 0) call htlu_dist_timers(D, t2 - t1, nth)
      stats%Time_Factor = stats%Time_Factor + t2 - t1
   end subroutine HTLU_dist_inverse

   !>**** xs = the entries of cluster a of the full vector set X (natural order in the cluster box)
   subroutine htlu_dist_seg(D, a, X, xs)
      type(htlu_dist) :: D
      integer a, f, m
      DT :: X(:, :)
      DT, allocatable :: xs(:, :)
      f = (a - 1)/D%nc + 1
      m = mod(a - 1, D%nc) + 1
      if (allocated(xs)) deallocate (xs)
      allocate (xs(product(D%csz(:, m)), size(X, 2)))
      call htlu_box(D%Ndim, D%Ns, D%coff(:, m), D%csz(:, m), X((f - 1)*D%nt + 1:f*D%nt, :), xs, 0)
   end subroutine htlu_dist_seg

   !>**** the entries of cluster a of the full vector set X become xs (mode 1) or are incremented by xs (mode 2)
   subroutine htlu_dist_put(D, a, xs, X, mode)
      type(htlu_dist) :: D
      integer a, mode, f, m
      DT :: X(:, :), xs(:, :)
      f = (a - 1)/D%nc + 1
      m = mod(a - 1, D%nc) + 1
      call htlu_box(D%Ndim, D%Ns, D%coff(:, m), D%csz(:, m), X((f - 1)*D%nt + 1:f*D%nt, :), xs, mode)
   end subroutine htlu_dist_put

   !>**** Y = op(A) X (op 'N' or 'T') for the distributed forest; X and Y are full vector sets (trees one after another, natural
   !>**** order in every tree), the same on all processes (collective)
   subroutine HTLU_dist_mult(D, op, X, Y)
      type(htlu_dist), target :: D
      character op
      DT :: X(:, :), Y(:, :)
      DT, allocatable :: xs(:, :), ys(:, :)
      integer R, i, j, a, b, ierr
      R = size(X, 2)
      call assert(size(X, 1) == D%P*D%nt .and. size(Y, 1) == D%P*D%nt, 'HTLU_dist_mult: not a full vector set')
      Y(1:D%P*D%nt, 1:R) = 0
      !$omp parallel
      !$omp single
      if (op == 'N') then
         do i = 1, D%nr
            !$omp task default(shared) firstprivate(i) private(j, a, b, xs, ys)
            a = htlu_dist_glob(D, 'r', i)
            allocate (ys(product(D%csz(:, mod(a - 1, D%nc) + 1)), R))
            ys = 0
            do j = 1, D%nl
               b = htlu_dist_glob(D, 'c', j)
               call htlu_dist_seg(D, b, X, xs)
               call htlu_mult(D%blk(i, j), 'N', xs, ys, BPACK_cone)
            enddo
            call htlu_dist_put(D, a, ys, Y, 2)
            deallocate (ys)
            if (allocated(xs)) deallocate (xs)
            !$omp end task
         enddo
      else
         do j = 1, D%nl
            !$omp task default(shared) firstprivate(j) private(i, a, b, xs, ys)
            b = htlu_dist_glob(D, 'c', j)
            allocate (ys(product(D%csz(:, mod(b - 1, D%nc) + 1)), R))
            ys = 0
            do i = 1, D%nr
               a = htlu_dist_glob(D, 'r', i)
               call htlu_dist_seg(D, a, X, xs)
               call htlu_mult(D%blk(i, j), 'T', xs, ys, BPACK_cone)
            enddo
            call htlu_dist_put(D, b, ys, Y, 2)
            deallocate (ys)
            if (allocated(xs)) deallocate (xs)
            !$omp end task
         enddo
      endif
      !$omp end single
      !$omp end parallel
      call MPI_Allreduce(MPI_IN_PLACE, Y, D%P*D%nt*R, MPI_DT, MPI_SUM, D%comm, ierr)
   end subroutine HTLU_dist_mult

   !>**** estimate of ||A||_2 for the distributed forest: nit power iterations with A^H A (collective)
   real(kind=8) function HTLU_dist_norm2(D, nit)
      type(htlu_dist), target :: D
      integer nit, it, ierr
      DT, allocatable :: x(:, :), y(:, :)
      real(kind=8) s
      allocate (x(D%P*D%nt, 1), y(D%P*D%nt, 1))
      if (D%myid == 0) call htlu_random(x)
      call MPI_Bcast(x, D%P*D%nt, MPI_DT, 0, D%comm, ierr)
      x = x/sqrt(sum(abs(x)**2))
      s = 0
      do it = 1, nit
         call HTLU_dist_mult(D, 'N', x, y)
         y = htlu_conj(y)
         call HTLU_dist_mult(D, 'T', y, x)
         x = htlu_conj(x)
         s = sqrt(sum(abs(x)**2))
         x = x/s
      enddo
      HTLU_dist_norm2 = sqrt(s)
      deallocate (x, y)
   end function HTLU_dist_norm2

   !>**** X(clusters r1..r2) = X(r1..r2) - sum_{b = c1..c2} G(a, b) X(b) for the full vector set X (collective)
   subroutine htlu_dist_matvec(D, r1, r2, c1, c2, X)
      type(htlu_dist), target :: D
      integer r1, r2, c1, c2
      DT :: X(:, :)
      DT, allocatable :: W(:, :), xs(:, :), ys(:, :)
      integer, allocatable :: woff(:)
      integer R, i, j, a, b, ierr
      R = size(X, 2)
      allocate (woff(r1:r2 + 1))
      woff(r1) = 0
      do a = r1, r2
         woff(a + 1) = woff(a) + product(D%csz(:, mod(a - 1, D%nc) + 1))
      enddo
      allocate (W(woff(r2 + 1), R))
      W = 0
      !$omp parallel
      !$omp single
      do i = 1, D%nr
         a = htlu_dist_glob(D, 'r', i)
         if (a < r1 .or. a > r2) cycle
         !$omp task default(shared) firstprivate(i, a) private(j, b, xs, ys)
         allocate (ys(woff(a + 1) - woff(a), R))
         ys = 0
         do j = 1, D%nl
            b = htlu_dist_glob(D, 'c', j)
            if (b < c1 .or. b > c2) cycle
            call htlu_dist_seg(D, b, X, xs)
            call htlu_mult(D%blk(i, j), 'N', xs, ys, BPACK_cone)
         enddo
         W(woff(a) + 1:woff(a + 1), :) = ys
         deallocate (ys)
         if (allocated(xs)) deallocate (xs)
         !$omp end task
      enddo
      !$omp end single
      !$omp end parallel
      call MPI_Allreduce(MPI_IN_PLACE, W, woff(r2 + 1)*R, MPI_DT, MPI_SUM, D%comm, ierr)
      do a = r1, r2
         allocate (ys(woff(a + 1) - woff(a), R))
         ys = -W(woff(a) + 1:woff(a + 1), :)
         call htlu_dist_put(D, a, ys, X, 2)
         deallocate (ys)
      enddo
      deallocate (W, woff)
   end subroutine htlu_dist_matvec

   !>**** X = S^{-1} X (htlu_bsolve) with the block-LU factored diagonal block S on the clusters lo..hi (collective)
   recursive subroutine htlu_dist_bsolve(D, lo, hi, X)
      type(htlu_dist), target :: D
      integer lo, hi
      DT :: X(:, :)
      DT, allocatable :: xs(:, :)
      integer, allocatable :: slo(:), shi(:)
      integer i, j, ns, ow, ierr
      if (lo == hi) then
         ow = htlu_dist_owner(D, lo, lo)
         call htlu_dist_seg(D, lo, X, xs)
         if (ow == D%myid) then
            !$omp parallel
            !$omp single
            call htlu_bsolve(D%blk((lo - 1)/D%nprow + 1, (lo - 1)/D%npcol + 1), xs)
            !$omp end single
            !$omp end parallel
         endif
         call MPI_Bcast(xs, size(xs), MPI_DT, ow, D%comm, ierr)
         call htlu_dist_put(D, lo, xs, X, 1)
         deallocate (xs)
         return
      endif
      call htlu_dist_sons(D, lo, hi, ns, slo, shi)
      do j = 2, ns
         call htlu_dist_matvec(D, slo(j), shi(j), lo, shi(j - 1), X)
      enddo
      do i = ns, 1, -1
         if (i < ns) call htlu_dist_matvec(D, slo(i), shi(i), slo(i + 1), hi, X)
         call htlu_dist_bsolve(D, slo(i), shi(i), X)
      enddo
      deallocate (slo, shi)
   end subroutine htlu_dist_bsolve

   !>**** X = A^{-1} X with the distributed block LU factors (HTLU_dist_factor); X is a full vector set, the same on all processes
   !>**** (collective)
   subroutine HTLU_dist_solve(D, X)
      type(htlu_dist), target :: D
      DT :: X(:, :)
      call htlu_dist_bsolve(D, 1, D%K, X)
   end subroutine HTLU_dist_solve

   !>**** the full vector set Xf from the local vectors Xl of the HTENSOR layout (the box D%boxes(:, rank) of every tree, trees one
   !>**** after another) (collective)
   subroutine HTLU_dist_l2g(D, Xl, Xf)
      type(htlu_dist) :: D
      DT :: Xl(:, :), Xf(:, :)
      integer R, ir, f, c, nl, ierr, Nd
      integer, allocatable :: cnt(:), dsp(:)
      DT, allocatable :: buf(:), sb(:), piece(:, :)
      Nd = D%Ndim
      R = size(Xl, 2)
      allocate (cnt(0:D%nproc - 1), dsp(0:D%nproc - 1))
      do ir = 0, D%nproc - 1
         cnt(ir) = D%P*product(D%boxes(Nd + 1:2*Nd, ir))*R
      enddo
      dsp(0) = 0
      do ir = 1, D%nproc - 1
         dsp(ir) = dsp(ir - 1) + cnt(ir - 1)
      enddo
      nl = product(D%boxes(Nd + 1:2*Nd, D%myid))
      allocate (sb(max(1, D%P*nl*R)), buf(max(1, sum(cnt))))
      do c = 1, R
         sb((c - 1)*D%P*nl + 1:c*D%P*nl) = Xl(1:D%P*nl, c)
      enddo
      call MPI_Allgatherv(sb, cnt(D%myid), MPI_DT, buf, cnt, dsp, MPI_DT, D%comm, ierr)
      do ir = 0, D%nproc - 1
         nl = product(D%boxes(Nd + 1:2*Nd, ir))
         allocate (piece(nl, 1))
         do c = 1, R
            do f = 1, D%P
               piece(:, 1) = buf(dsp(ir) + (c - 1)*D%P*nl + (f - 1)*nl + 1:dsp(ir) + (c - 1)*D%P*nl + f*nl)
               call htlu_box(Nd, D%Ns, D%boxes(1:Nd, ir), D%boxes(Nd + 1:2*Nd, ir), Xf((f - 1)*D%nt + 1:f*D%nt, c:c), piece, 1)
            enddo
         enddo
         deallocate (piece)
      enddo
      deallocate (cnt, dsp, sb, buf)
   end subroutine HTLU_dist_l2g

   !>**** the local vectors Xl of the HTENSOR layout from the full vector set Xf
   subroutine HTLU_dist_g2l(D, Xf, Xl)
      type(htlu_dist) :: D
      DT :: Xf(:, :), Xl(:, :)
      integer f, nl, Nd
      Nd = D%Ndim
      nl = product(D%boxes(Nd + 1:2*Nd, D%myid))
      do f = 1, D%P
         call htlu_box(Nd, D%Ns, D%boxes(1:Nd, D%myid), D%boxes(Nd + 1:2*Nd, D%myid), Xf((f - 1)*D%nt + 1:f*D%nt, :), &
            Xl((f - 1)*nl + 1:f*nl, :), 0)
      enddo
   end subroutine HTLU_dist_g2l

   !>**** memory (MB) of the dense and the Tucker blocks summed over the processes, the largest memory of one process, the maximum
   !>**** Tucker rank per level and the numbers of dense and Tucker blocks (collective)
   subroutine HTLU_dist_stats(D, memd, memt, memmax, rk, nd, nt)
      type(htlu_dist) :: D
      real(kind=8) memd, memt, memmax, v(2), vs(2)
      integer rk(0:), nd, nt, i, j, ierr, n2(2), n2s(2)
      integer, allocatable :: rks(:)
      memd = 0
      memt = 0
      rk = 0
      nd = 0
      nt = 0
      do j = 1, D%nl
         do i = 1, D%nr
            call HTLU_stats(D%blk(i, j), memd, memt, rk, nd, nt)
         enddo
      enddo
      v(1) = memd + memt
      call MPI_Allreduce(v(1), memmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, D%comm, ierr)
      v = (/memd, memt/)
      call MPI_Allreduce(v, vs, 2, MPI_DOUBLE_PRECISION, MPI_SUM, D%comm, ierr)
      memd = vs(1)
      memt = vs(2)
      n2 = (/nd, nt/)
      call MPI_Allreduce(n2, n2s, 2, MPI_INTEGER, MPI_SUM, D%comm, ierr)
      nd = n2s(1)
      nt = n2s(2)
      allocate (rks(size(rk)))
      call MPI_Allreduce(rk, rks, size(rk), MPI_INTEGER, MPI_MAX, D%comm, ierr)
      rk = rks
      deallocate (rks)
   end subroutine HTLU_dist_stats

   !>**** print the Tucker profile of the distributed forest (by level and part: strictly lower, strictly upper, diagonal path)
   !>**** on rank 0 (collective)
   subroutine HTLU_dist_print_profile(D, lmax, title)
      type(htlu_dist) :: D
      integer lmax, part, l, i, j, a, b, ierr
      character(*) title
      integer cnt(0:lmax, 3), rmax(0:lmax, 3, 2), cnts(0:lmax, 3), rmaxs(0:lmax, 3, 2)
      real(kind=8) rsum(0:lmax, 3, 2), mem(0:lmax, 3), rsums(0:lmax, 3, 2), mems(0:lmax, 3)
      character(len=5), parameter :: pname(3) = (/'lower', 'upper', 'diag '/)
      cnt = 0
      rmax = 0
      rsum = 0
      mem = 0
      do j = 1, D%nl
         b = htlu_dist_glob(D, 'c', j)
         do i = 1, D%nr
            a = htlu_dist_glob(D, 'r', i)
            part = 0
            if (a > b) part = 1
            if (a < b) part = 2
            call HTLU_rank_profile(D%blk(i, j), part, cnt, rsum, rmax, mem)
         enddo
      enddo
      call MPI_Reduce(cnt, cnts, size(cnt), MPI_INTEGER, MPI_SUM, 0, D%comm, ierr)
      call MPI_Reduce(rmax, rmaxs, size(rmax), MPI_INTEGER, MPI_MAX, 0, D%comm, ierr)
      call MPI_Reduce(rsum, rsums, size(rsum), MPI_DOUBLE_PRECISION, MPI_SUM, 0, D%comm, ierr)
      call MPI_Reduce(mem, mems, size(mem), MPI_DOUBLE_PRECISION, MPI_SUM, 0, D%comm, ierr)
      if (D%myid /= 0) return
      write (*, '(A,A)') ' Tucker profile, ', title
      write (*, '(A)') '   level  part      blocks   row rank (mean, max)   col rank (mean, max)   memory (MB)'
      do l = 0, lmax
         do part = 1, 3
            if (cnts(l, part) == 0) cycle
            write (*, '(I8,2X,A5,I12,F14.1,I8,F15.1,I8,F14.1)') l, pname(part), cnts(l, part), rsums(l, part, 1)/cnts(l, part), &
               rmaxs(l, part, 1), rsums(l, part, 2)/cnts(l, part), rmaxs(l, part, 2), mems(l, part)
         enddo
      enddo
   end subroutine HTLU_dist_print_profile

   !>**** free the distributed forest
   subroutine HTLU_dist_delete(D)
      type(htlu_dist) :: D
      integer i, j
      if (allocated(D%blk)) then
         do j = 1, size(D%blk, 2)
            do i = 1, size(D%blk, 1)
               call HTLU_delete(D%blk(i, j))
            enddo
         enddo
         deallocate (D%blk)
      endif
      if (allocated(D%rmap)) deallocate (D%rmap)
      if (allocated(D%Ns)) deallocate (D%Ns)
      if (allocated(D%coff)) deallocate (D%coff)
      if (allocated(D%csz)) deallocate (D%csz)
      if (allocated(D%cgrp)) deallocate (D%cgrp)
      if (allocated(D%boxes)) deallocate (D%boxes)
      D%K = 0
      D%inverse = 0
   end subroutine HTLU_dist_delete

end module Bplus_htensor_lu
