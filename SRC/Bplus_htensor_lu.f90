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
!> The factorization is the H-matrix LU (as Hmat_LU) on these trees, with a block LU over the trees on top: LU of the
!> diagonal blocks, triangular solves (LXM, XUM) and Schur updates (C = C - A*B). Pivoting is only inside the dense
!> diagonal leaves. Tucker arithmetic, truncated at option%tol_rand:
!>   - a triangular solve with a Tucker block applies the solve to its row (column) Kronecker factor and recompresses it
!>     by a truncated HOSVD of the result;
!>   - a Schur update of a Tucker block rebuilds the block from products: per-mode bases from a few random products of
!>     the updated block (every product gives all fibers of its output) and the core by an exact projection;
!>   - a product with a Tucker factor that updates a hierarchical block is formed in Tucker form and added block by block.
!> A variant is the block LU with whole-cluster pivots (htlu_BLU): the lower blocks are A_ji S_ii^{-1} with S_ii the whole
!> current diagonal block, the upper blocks stay the Schur-updated blocks, and solves go through the block factors.
!> The inverse (HTLU_forest_inverse) is an in-place block Gauss-Jordan sweep on the same trees with the same kernels.
!> Singular values are dropped below option%tol_rand relative to each block and below option%tol_rand*normest, with
!> normest an estimate of ||A||_2 (LU) or ||A^{-1}||_2 (inverse) given by the caller.
!> The current implementation is for one MPI process (OpenMP tasks over independent blocks).
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

   !>**** blocks above this size (rows) create OpenMP tasks for their children
   integer, parameter :: htlu_task_rows = 1024
   !>**** products with at most this many multiply-adds use matmul instead of BLAS (whose per-call overhead dominates for the
   !>**** small Tucker factors and cores, and whose buffer allocation contends between threads)
   real(kind=8), parameter :: htlu_small_gemm = 32768d0

   !>**** thread-seconds and counts of the leaf operations of the factorization (diagnostics): 1 Tucker rebuilds of Schur
   !>**** updates, 2 triangular solves of Tucker blocks, 3 Tucker products, 4 Tucker multiplications of the inverse, 5 dense leaf
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

   !>**** X = L^{-1} X with the (unit) lower factor and the leaf pivots of the diagonal block L
   recursive subroutine htlu_lsolve(L, X)
      type(htlu_node) :: L
      DT :: X(:, :)
      type(htlu_fmat), allocatable :: Xc(:)
      DT, allocatable :: row(:)
      integer i, k, n
      if (size(X, 2) == 0) return
      if (L%style == 1) then
         n = size(L%D, 1)
         allocate (row(size(X, 2)))
         do i = 1, n
            if (L%ipiv(i) /= i) then
               row = X(i, :)
               X(i, :) = X(L%ipiv(i), :)
               X(L%ipiv(i), :) = row
            endif
         enddo
         deallocate (row)
         call trsmf90(L%D, X, 'L', 'L', 'N', 'U', n, size(X, 2))
      else
         call assert(L%style == 4, 'htlu_lsolve: a Tucker diagonal block')
         allocate (Xc(L%nch))
         call htlu_gather(L, 'r', X, Xc)
         do i = 1, L%nch
            do k = 1, i - 1
               call htlu_mult(L%sons(i, k), 'N', Xc(k)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
            call htlu_lsolve(L%sons(i, i), Xc(i)%mat)
         enddo
         call htlu_scatter(L, 'r', Xc, X, 1)
         deallocate (Xc)
      endif
   end subroutine htlu_lsolve

   !>**** X = U^{-1} X with the upper factor of the diagonal block U
   recursive subroutine htlu_usolve(U, X)
      type(htlu_node) :: U
      DT :: X(:, :)
      type(htlu_fmat), allocatable :: Xc(:)
      integer i, k, n
      if (size(X, 2) == 0) return
      if (U%style == 1) then
         n = size(U%D, 1)
         call trsmf90(U%D, X, 'L', 'U', 'N', 'N', n, size(X, 2))
      else
         call assert(U%style == 4, 'htlu_usolve: a Tucker diagonal block')
         allocate (Xc(U%nch))
         call htlu_gather(U, 'r', X, Xc)
         do i = U%nch, 1, -1
            do k = i + 1, U%nch
               call htlu_mult(U%sons(i, k), 'N', Xc(k)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
            call htlu_usolve(U%sons(i, i), Xc(i)%mat)
         enddo
         call htlu_scatter(U, 'r', Xc, X, 1)
         deallocate (Xc)
      endif
   end subroutine htlu_usolve

   !>**** X = U^{-T} X with the upper factor of the diagonal block U
   recursive subroutine htlu_usolve_T(U, X)
      type(htlu_node) :: U
      DT :: X(:, :)
      type(htlu_fmat), allocatable :: Xc(:)
      integer i, k, n
      if (size(X, 2) == 0) return
      if (U%style == 1) then
         n = size(U%D, 1)
         call trsmf90(U%D, X, 'L', 'U', 'T', 'N', n, size(X, 2))
      else
         call assert(U%style == 4, 'htlu_usolve_T: a Tucker diagonal block')
         allocate (Xc(U%nch))
         call htlu_gather(U, 'c', X, Xc)
         do i = 1, U%nch
            do k = 1, i - 1
               call htlu_mult(U%sons(k, i), 'T', Xc(k)%mat, Xc(i)%mat, -BPACK_cone)
            enddo
            call htlu_usolve_T(U%sons(i, i), Xc(i)%mat)
         enddo
         call htlu_scatter(U, 'c', Xc, X, 1)
         deallocate (Xc)
      endif
   end subroutine htlu_usolve_T

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
      integer i, j, d, Ndim
      integer ro(A%Ndim), rs(A%Ndim), co(A%Ndim), cs(A%Ndim)
      Ndim = A%Ndim
      call htlu_child_box(A, 'r', i, ro, rs)
      call htlu_child_box(A, 'c', j, co, cs)
      Pc%style = 2
      Pc%Ndim = Ndim
      Pc%level = A%level + 1
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
   end subroutine htlu_restrict

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

   !>**** M = L^{-1} M for the diagonal block L of the row clusters of M
   recursive subroutine htlu_LXM(L, M, tol)
      type(htlu_node), target :: L, M
      real(kind=8) tol(2)
      type(htlu_fmat), allocatable :: Q(:), S(:)
      DT, allocatable :: W(:, :), G(:, :), Z0(:, :), Rk(:, :)
      integer i, j, k, d
      real(kind=8) t0
      t0 = MPI_Wtime()
      select case (M%style)
      case (1)
         call htlu_lsolve(L, M%D)
         call htlu_tadd(5, MPI_Wtime() - t0)
      case (2)
         ! L^{-1} M = Z kron(QV)^T with V_M = QV RV and Z = (L^{-1} kron(U_M)) C_M kron(RV)^T ~ kron(Q) G
         allocate (Q(M%Ndim), S(M%Ndim))
         call htlu_orth(M%Ndim, M%V, S, Rk)
         call htlu_gemm('N', 'T', M%C, Rk, Z0)
         call htlu_kron(M%Ndim, M%U, W, 0)
         call htlu_lsolve(L, W)
         call htlu_gemm('N', 'N', W, Z0, Rk)
         call htlu_tuckerize(M%Ndim, M%m, Rk, tol, Q, G)
         do d = 1, M%Ndim
            call move_alloc(Q(d)%mat, M%U(d)%mat)
            call move_alloc(S(d)%mat, M%V(d)%mat)
         enddo
         call move_alloc(G, M%C)
         deallocate (Q, S, W, Z0, Rk)
         call htlu_core_trim(M, tol)
         call htlu_tadd(2, MPI_Wtime() - t0)
      case (4)
         call assert(L%style == 4, 'htlu_LXM: a hierarchical block under a leaf')
         do j = 1, M%nch
            !$omp task default(shared) firstprivate(j) private(i,k) if(product(M%m) >= htlu_task_rows)
            do i = 1, M%nch
               do k = 1, i - 1
                  call htlu_subhh(M%sons(i, j), L%sons(i, k), M%sons(k, j), tol)
               enddo
               call htlu_LXM(L%sons(i, i), M%sons(i, j), tol)
            enddo
            !$omp end task
         enddo
         !$omp taskwait
      end select
   end subroutine htlu_LXM

   !>**** M = M U^{-1} for the diagonal block U of the column clusters of M
   recursive subroutine htlu_XUM(U, M, tol)
      type(htlu_node), target :: U, M
      real(kind=8) tol(2)
      type(htlu_fmat), allocatable :: Q(:), S(:)
      DT, allocatable :: W(:, :), G(:, :), Z0(:, :), Rk(:, :)
      integer i, j, k, d
      real(kind=8) t0
      t0 = MPI_Wtime()
      select case (M%style)
      case (1)
         call trsmf90(U%D, M%D, 'R', 'U', 'N', 'N', size(M%D, 1), size(M%D, 2))
         call htlu_tadd(5, MPI_Wtime() - t0)
      case (2)
         ! M U^{-1} = kron(QU) Z^T with U_M = QU RU and Z = (U^{-T} kron(V_M)) (kron(RU) C_M)^T ~ kron(Q) G
         allocate (Q(M%Ndim), S(M%Ndim))
         call htlu_orth(M%Ndim, M%U, S, Rk)
         call htlu_gemm('N', 'N', Rk, M%C, Z0)
         call htlu_kron(M%Ndim, M%V, W, 0)
         call htlu_usolve_T(U, W)
         call htlu_gemm('N', 'T', W, Z0, Rk)
         call htlu_tuckerize(M%Ndim, M%n, Rk, tol, Q, G)
         do d = 1, M%Ndim
            call move_alloc(S(d)%mat, M%U(d)%mat)
            call move_alloc(Q(d)%mat, M%V(d)%mat)
         enddo
         deallocate (M%C)
         allocate (M%C(size(G, 2), size(G, 1)))
         M%C = transpose(G)
         deallocate (Q, S, W, Z0, Rk, G)
         call htlu_core_trim(M, tol)
         call htlu_tadd(2, MPI_Wtime() - t0)
      case (4)
         call assert(U%style == 4, 'htlu_XUM: a hierarchical block under a leaf')
         do i = 1, M%nch
            !$omp task default(shared) firstprivate(i) private(j,k) if(product(M%m) >= htlu_task_rows)
            do j = 1, M%nch
               do k = 1, j - 1
                  call htlu_subhh(M%sons(i, j), M%sons(i, k), U%sons(k, j), tol)
               enddo
               call htlu_XUM(U%sons(j, j), M%sons(i, j), tol)
            enddo
            !$omp end task
         enddo
         !$omp taskwait
      end select
   end subroutine htlu_XUM

   !>**** LU factorization of the diagonal block A in place
   recursive subroutine htlu_LU(A, tol)
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
      call assert(A%style == 4, 'htlu_LU: a Tucker diagonal block')
      do i = 1, A%nch
         call htlu_LU(A%sons(i, i), tol)
         do j = i + 1, A%nch
            !$omp task default(shared) firstprivate(i,j) if(product(A%m) >= htlu_task_rows)
            call htlu_LXM(A%sons(i, i), A%sons(i, j), tol)
            !$omp end task
            !$omp task default(shared) firstprivate(i,j) if(product(A%m) >= htlu_task_rows)
            call htlu_XUM(A%sons(i, i), A%sons(j, i), tol)
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
   end subroutine htlu_LU

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

   !>**** the P x P forest A becomes its inverse (H-TENSOR inverse with Tucker blocks; apply it with HTLU_forest_mult)
   !>**** (normest: estimate of ||A^{-1}||_2; Tucker truncation at option%tol_rand relative to each block and to normest)
   subroutine HTLU_forest_inverse(P, A, option, stats, ptree, normest)
      integer P, p1, q1
      type(htlu_node), target :: A(:, :)
      type(Hoption) :: option
      type(Hstat) :: stats
      type(proctree) :: ptree
      real(kind=8) normest
      integer nth
      real(kind=8) t1, t2, tol(2)
      call assert(ptree%nproc == 1, 'HTLU_forest_inverse: one MPI process only')
      tol = (/option%tol_rand, option%tol_rand*normest/)
      nth = 1
!$    nth = omp_get_max_threads()
      t1 = MPI_Wtime()
      !$omp parallel
      !$omp single
      call htlu_sweep(P, A, tol)
      !$omp end single
      !$omp end parallel
      do q1 = 1, P
         do p1 = 1, P
            call htlu_scale(A(p1, q1), -BPACK_cone)
         enddo
      enddo
      t2 = MPI_Wtime()
      if (option%verbosity >= 0) call HTLU_print_timers(t2 - t1, nth)
      stats%Time_Factor = stats%Time_Factor + t2 - t1
   end subroutine HTLU_forest_inverse

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

   !>**** print the Tucker profile of the P x P forest A (rank 0 only)
   subroutine HTLU_print_profile(P, A, lmax, title)
      integer P, lmax, p1, q1, part, l
      type(htlu_node) :: A(:, :)
      character(*) title
      integer cnt(0:lmax, 3), rmax(0:lmax, 3, 2)
      real(kind=8) rsum(0:lmax, 3, 2), mem(0:lmax, 3)
      character(len=5), parameter :: pname(3) = (/'lower', 'upper', 'diag '/)
      cnt = 0
      rmax = 0
      rsum = 0
      mem = 0
      do q1 = 1, P
         do p1 = 1, P
            part = 0
            if (p1 > q1) part = 1
            if (p1 < q1) part = 2
            call HTLU_rank_profile(A(p1, q1), part, cnt, rsum, rmax, mem)
         enddo
      enddo
      write (*, '(A,A)') ' Tucker profile, ', title
      write (*, '(A)') '   level  part      blocks   row rank (mean, max)   col rank (mean, max)   memory (MB)'
      do l = 0, lmax
         do part = 1, 3
            if (cnt(l, part) == 0) cycle
            write (*, '(I8,2X,A5,I12,F14.1,I8,F15.1,I8,F14.1)') l, pname(part), cnt(l, part), rsum(l, part, 1)/cnt(l, part), &
               rmax(l, part, 1), rsum(l, part, 2)/cnt(l, part), rmax(l, part, 2), mem(l, part)
         enddo
      enddo
   end subroutine HTLU_print_profile

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

   !>**** copy the blocks of an HTENSOR (one process) into the recursive tree A
   subroutine HTLU_from_HTENSOR(A, h_mat_md, msh, Ndim)
      type(htlu_node) :: A
      type(Hmat_md), target :: h_mat_md
      integer Ndim
      type(mesh) :: msh(Ndim)
      type(htlu_index_level), allocatable :: idx(:)
      type(matrixblock_MD), pointer :: blk
      integer ll, bb, l, Lmax, ri, k, g1(Ndim)
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
      g1 = 1
      call htlu_build(A, Ndim, 0, g1, g1, msh, h_mat_md, idx, Lmax)
      do l = 0, Lmax
         deallocate (idx(l)%head, idx(l)%next, idx(l)%cidx, idx(l)%ll, idx(l)%bb)
      enddo
      deallocate (idx)
   end subroutine HTLU_from_HTENSOR

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
            call assert(blk%style == 2 .and. blk%level_butterfly == 0, 'HTLU_from_HTENSOR: only Tucker (lrlevel 0) admissible blocks')
            call assert(.not. allocated(blk%MiddleZFP) .and. .not. allocated(blk%MiddleQTT), 'HTLU_from_HTENSOR: ZFP/QTT cores')
            call assert(size(blk%ButterflyU(1)%blocks, 1) == 1 .and. size(blk%ButterflyV(1)%blocks, 1) == 1, 'HTLU_from_HTENSOR: Tucker factors')
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
         call assert(2*rg(1) <= msh(1)%Maxgroup, 'HTLU_from_HTENSOR: a leaf-level block is missing')
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

   !>**** Y = A*X for the P x P forest A; X and Y hold the trees one after another (nloc entries each)
   subroutine HTLU_forest_mult(P, A, nloc, X, Y)
      integer P, nloc, p1, q1, R
      type(htlu_node) :: A(:, :)
      DT :: X(:, :), Y(:, :)
      DT, allocatable :: xq(:, :), yp(:, :)
      R = size(X, 2)
      allocate (xq(nloc, R), yp(nloc, R))
      Y(1:P*nloc, 1:R) = 0
      !$omp parallel
      !$omp single
      do p1 = 1, P
         yp = 0
         do q1 = 1, P
            xq = X((q1 - 1)*nloc + 1:q1*nloc, 1:R)
            call htlu_mult(A(p1, q1), 'N', xq, yp, BPACK_cone)
         enddo
         Y((p1 - 1)*nloc + 1:p1*nloc, 1:R) = yp
      enddo
      !$omp end single
      !$omp end parallel
      deallocate (xq, yp)
   end subroutine HTLU_forest_mult

   !>**** Y = A^T X for the P x P forest A
   subroutine HTLU_forest_multT(P, A, nloc, X, Y)
      integer P, nloc, p1, q1, R
      type(htlu_node) :: A(:, :)
      DT :: X(:, :), Y(:, :)
      DT, allocatable :: xp(:, :), yq(:, :)
      R = size(X, 2)
      allocate (xp(nloc, R), yq(nloc, R))
      !$omp parallel
      !$omp single
      do q1 = 1, P
         yq = 0
         do p1 = 1, P
            xp = X((p1 - 1)*nloc + 1:p1*nloc, 1:R)
            call htlu_mult(A(p1, q1), 'T', xp, yq, BPACK_cone)
         enddo
         Y((q1 - 1)*nloc + 1:q1*nloc, 1:R) = yq
      enddo
      !$omp end single
      !$omp end parallel
      deallocate (xp, yq)
   end subroutine HTLU_forest_multT

   !>**** estimate of ||A||_2 for the P x P forest A: nit power iterations with A^H A
   real(kind=8) function HTLU_forest_norm2(P, A, nloc, nit)
      integer P, nloc, nit, it
      type(htlu_node) :: A(:, :)
      DT, allocatable :: x(:, :), y(:, :)
      real(kind=8) s
      allocate (x(P*nloc, 1), y(P*nloc, 1))
      call htlu_random(x)
      x = x/sqrt(sum(abs(x)**2))
      s = 0
      do it = 1, nit
         call HTLU_forest_mult(P, A, nloc, x, y)
         y = htlu_conj(y)
         call HTLU_forest_multT(P, A, nloc, y, x)
         x = htlu_conj(x)
         s = sqrt(sum(abs(x)**2))
         x = x/s
      enddo
      HTLU_forest_norm2 = sqrt(s)
      deallocate (x, y)
   end function HTLU_forest_norm2

   !>**** LU factorization of the P x P forest A in place: block LU over the trees, H-matrix LU inside (blocklu = 0) or block LU
   !>**** with whole-cluster pivots inside (blocklu = 1, see htlu_bsolve)
   !>**** (normest: estimate of ||A||_2; Tucker truncation at option%tol_rand relative to each block and to normest)
   subroutine HTLU_forest_factor(P, A, option, stats, ptree, normest, blocklu)
      integer P, p1, q1, r1, blocklu
      type(htlu_node), target :: A(:, :)
      type(Hoption) :: option
      type(Hstat) :: stats
      type(proctree) :: ptree
      real(kind=8) normest
      integer nth
      real(kind=8) t1, t2, tol(2)
      call assert(ptree%nproc == 1, 'HTLU_forest_factor: one MPI process only')
      tol = (/option%tol_rand, option%tol_rand*normest/)
      nth = 1
!$    nth = omp_get_max_threads()
      t1 = MPI_Wtime()
      !$omp parallel
      !$omp single
      do p1 = 1, P
         if (blocklu == 1) then
            call htlu_BLU(A(p1, p1), tol)
            do q1 = p1 + 1, P
               !$omp task default(shared) firstprivate(p1,q1)
               call htlu_XSM(A(p1, p1), A(q1, p1), tol)
               !$omp end task
            enddo
         else
            call htlu_LU(A(p1, p1), tol)
            do q1 = p1 + 1, P
               !$omp task default(shared) firstprivate(p1,q1)
               call htlu_LXM(A(p1, p1), A(p1, q1), tol)
               !$omp end task
               !$omp task default(shared) firstprivate(p1,q1)
               call htlu_XUM(A(p1, p1), A(q1, p1), tol)
               !$omp end task
            enddo
         endif
         !$omp taskwait
         do r1 = p1 + 1, P
            do q1 = p1 + 1, P
               !$omp task default(shared) firstprivate(p1,q1,r1)
               call htlu_subhh(A(q1, r1), A(q1, p1), A(p1, r1), tol)
               !$omp end task
            enddo
         enddo
         !$omp taskwait
      enddo
      !$omp end single
      !$omp end parallel
      t2 = MPI_Wtime()
      if (option%verbosity >= 0) call HTLU_print_timers(t2 - t1, nth)
      stats%Time_Factor = stats%Time_Factor + t2 - t1
   end subroutine HTLU_forest_factor

   !>**** X = A^{-1} X with the LU (blocklu = 0) or the block LU (blocklu = 1) of the P x P forest A (trees one after another,
   !>**** nloc entries each)
   subroutine HTLU_forest_solve(P, A, nloc, X, blocklu)
      integer P, nloc, p1, q1, R, blocklu
      type(htlu_node) :: A(:, :)
      DT :: X(:, :)
      DT, allocatable :: xq(:, :), xp(:, :)
      R = size(X, 2)
      allocate (xq(nloc, R), xp(nloc, R))
      !$omp parallel
      !$omp single
      do p1 = 1, P
         xp = X((p1 - 1)*nloc + 1:p1*nloc, 1:R)
         do q1 = 1, p1 - 1
            xq = X((q1 - 1)*nloc + 1:q1*nloc, 1:R)
            call htlu_mult(A(p1, q1), 'N', xq, xp, -BPACK_cone)
         enddo
         if (blocklu == 0) call htlu_lsolve(A(p1, p1), xp)
         X((p1 - 1)*nloc + 1:p1*nloc, 1:R) = xp
      enddo
      do p1 = P, 1, -1
         xp = X((p1 - 1)*nloc + 1:p1*nloc, 1:R)
         do q1 = p1 + 1, P
            xq = X((q1 - 1)*nloc + 1:q1*nloc, 1:R)
            call htlu_mult(A(p1, q1), 'N', xq, xp, -BPACK_cone)
         enddo
         if (blocklu == 1) then
            call htlu_bsolve(A(p1, p1), xp)
         else
            call htlu_usolve(A(p1, p1), xp)
         endif
         X((p1 - 1)*nloc + 1:p1*nloc, 1:R) = xp
      enddo
      !$omp end single
      !$omp end parallel
      deallocate (xq, xp)
   end subroutine HTLU_forest_solve

end module Bplus_htensor_lu
