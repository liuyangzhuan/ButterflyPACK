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

!> @file BPACK_solve_mul.f90
!> @brief Top-level subroutines for multiplying a BPACK (H/HODBF/HODLR/HSS-BF) matrix (or its inverse or triangular factors) with vectors and associated communication routines

#include "ButterflyPACK_config.fi"
module BPACK_Solve_Mul
   use BPACK_DEFS
   use Bplus_compress

contains

!>**** eigen solver using ARPACK
   !bmat_shift,option_sh,ptree_sh,stats_sh: matrix, option, process tree and statistics for A-sigmaI or A-sigma*B and its inverse, not referenced if SI=0
   !bmat_A,option_A,ptree_A,stats_A: matrix, option, process tree and statistics for A
   !bmat_B,option_B,ptree_B,stats_B: matrix, option, process tree and statistics for B, not referenced if CMmode=0
   !CMmode: 0: solve the eigen mode with AV=VD; 1: solve the characteristic mode with AV=BVD, B=real(A)
   !SI: 0: regular mode 1: shift-invert mode
   !Nunk: size of A
   !Nunk_loc: local size of A
   !nev: number of eigen value required
   !nconv: number of eigen value converged
   !tol: tolerance in ARPACK
   !shift: shift value, not referenced if SI=0
   !eigval: returned eigenvalues of size(nev)
   !eigvec: returned eigenvectors  of size(Nunk_loc x nev)
   !which: which eigenvlaues are computed: 'LM', 'SM', LR', 'SR', LI', 'SI'
   subroutine BPACK_Eigen(bmat_A, option_A, ptree_A, stats_A, bmat_B, option_B, ptree_B, stats_B, bmat_sh, option_sh, ptree_sh, stats_sh, Nunk, Nunk_loc, nev, tol, CMmode, SI, shift, which, nconv, eigval, eigvec)


      implicit none

      integer Nunk, Nunk_loc, nev, nconv, CMmode, SI
      character(len=2) which
      complex(kind=8) shift
      complex(kind=8) eigval(nev)
      complex(kind=8) eigvec(Nunk_loc, nev)
      integer i, j, ii, jj, iii, jjj
      real(kind=8):: rel_error, t1, t2
      type(Hoption)::option_sh, option_A, option_B
      type(proctree)::ptree_sh, ptree_A, ptree_B
      type(Hstat)::stats_sh, stats_A, stats_B
      type(Bmatrix)::bmat_sh, bmat_A, bmat_B
      real(kind=8) n1, n2, rtemp

      integer maxn, maxnev, maxncv, ldv, ierr
      integer iparam(11), ipntr(14)
      logical, allocatable:: select(:)
      complex(kind=8), allocatable:: ax(:), mx(:), d(:), v(:, :), workd(:), workev(:), resid(:), workl(:)
      real(kind=8), allocatable:: rwork(:), rd(:, :)

      character bmattype
      integer ido, n, nx, ncv, lworkl, info, maxitr, ishfts, mode
      complex(kind=8) sigma
      real(kind=8) tol
      logical rvec
      real(kind=8), external :: pdznorm2, dlapy2

      nconv = 0

#ifdef HAVE_ARPACK

#if DAT==0

      t1 = MPI_Wtime()
      maxn = Nunk
      ldv = Nunk_loc
      maxnev = nev
      ! maxncv = nev*2
      maxncv = min(nev*4,Nunk)

      !!!! the following if test has been disabled. Make sure "else if (ncv .le. nev+1 .or.  ncv .gt. n) then" in pzneupd.f is changed to "else if (ncv .le. nev+1) then"
      ! if(maxncv>Nunk_loc)then
      !    if(ptree_A%MyID==Main_ID)print *, ' PARPACK requires ncv<=Nunk_loc. Please reduce MPI count or nev'
      !    stop
      ! endif

      n = maxn
      ncv = maxncv
      lworkl = 3*ncv**2 + 5*ncv
      ido = 0
      info = 0

      ishfts = 1
      maxitr = 3000
      iparam(1) = ishfts
      iparam(3) = maxitr

      allocate (select(maxncv))
      allocate (ax(Nunk_loc))
      allocate (mx(Nunk_loc))
      allocate (d(maxn))
      allocate (v(ldv, maxncv))
      allocate (workd(3*Nunk_loc))
      workd=0
      allocate (workev(3*maxncv))
      allocate (resid(maxn))
      allocate (workl(3*maxncv*maxncv + 5*maxncv))
      allocate (rwork(maxncv))
      allocate (rd(maxncv, 3))

      if (CMmode == 0) then ! solve the eigen mode

         do while (ido /= 99)

            if (SI == 0) then  ! regular mode
               bmattype = 'I'
               ! which='LM'  ! largest eigenvalues
               iparam(7) = 1
               sigma = 0d0
               call pznaupd(ptree_A%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
               if (ido == -1 .or. ido == 1) then !Perform  y <--- OP*x = A*x
                  call BPACK_Mult('N', Nunk_loc, 1, workd(ipntr(1)), workd(ipntr(2)), bmat_A, ptree_A, option_A, stats_A)
               endif
            else                ! shift-invert mode
               bmattype = 'I'
               ! which='LM' ! eigenvalues closest to the shifts
               iparam(7) = 3
               sigma = shift
               call pznaupd(ptree_A%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
               if (ido == -1 .or. ido == 1) then !Perform  y <--- OP*x = inv[A-SIGMA*I]*x
                  call BPACK_Solution(bmat_sh, workd(ipntr(2)), workd(ipntr(1)), Nunk_loc, 1, option_sh, ptree_sh, stats_sh)
               endif
            endif

         enddo

         if (info < 0) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd.'
            print *, ' '
         else
            rvec = .true.
            call pzneupd(ptree_A%Comm, rvec, 'A', select, d, v, ldv, sigma, workev, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
            if (ierr /= 0) then
               print *, ' '
               print *, ' Error with _neupd, info = ', ierr
               print *, ' Check the documentation of _neupd. '
               print *, ' '
               if(ierr==-3)then
                  print *, ' Consider modifying pzneupd.f in PARPACK to disable this error'
               endif
            else
               nconv = iparam(5)
               do j = 1, nconv

                  call BPACK_Mult('N', Nunk_loc, 1, v(1, j), ax, bmat_A, ptree_A, option_A, stats_A)
                  mx(1:Nunk_loc) = v(1:Nunk_loc, j)
                  call zaxpy(Nunk_loc, -d(j), mx, 1, ax, 1)
                  rd(j, 1) = dble(d(j))
                  rd(j, 2) = dimag(d(j))
                  rd(j, 3) = pdznorm2(ptree_A%Comm, Nunk_loc, ax, 1)
                  rd(j, 3) = rd(j, 3)/dlapy2(rd(j, 1), rd(j, 2))
               enddo
               if (option_A%verbosity>=0) call pdmout(ptree_A%Comm, 6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and direct residuals')

               eigval(1:nconv) = d(1:nconv)
               eigvec(1:Nunk_loc, 1:nconv) = v(1:Nunk_loc, 1:nconv)

            end if

            if (ptree_A%MyID == Main_ID .and. option_A%verbosity>=0) then
               ! write(*,*)resid(1:nconv)
               if (info .eq. 1) then
                  print *, ' '
                  print *, ' Maximum number of iterations reached.'
                  print *, ' '
               else if (info .eq. 3) then
                  print *, ' '
                  print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
                  print *, ' '
               end if
               print *, ' '
               print *, '_NDRV1'
               print *, '====== '
               print *, ' '
               print *, ' Size of the matrix is ', n
               print *, ' The number of processors is ', ptree_A%nproc
               print *, ' The number of Ritz values requested is ', nev
               print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
               print *, ' What portion of the spectrum: ', which
               print *, ' The number of converged Ritz values is ', nconv
               print *, ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
               print *, ' The number of OP*x is ', iparam(9)
               print *, ' The convergence criterion is ', tol
               print *, ' '
            endif
         end if

      else if (CMmode == 1) then ! solve the characteristic mode
         do while (ido /= 99)

            if (SI == 0) then ! regular mode
               bmattype = 'G'
               ! which='LM'         ! largest eigenvalues
               iparam(7) = 2
               sigma = 0
               call pznaupd(ptree_A%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
               if (ido == -1 .or. ido == -1) then !Perform  y <--- OP*x = inv[M]*A*x
                  call BPACK_Mult('N', Nunk_loc, 1, workd(ipntr(1)), workd(ipntr(2) + Nunk_loc), bmat_A, ptree_A, option_A, stats_A)
                  call BPACK_Solution(bmat_B, workd(ipntr(2)), workd(ipntr(2) + Nunk_loc), Nunk_loc, 1, option_B, ptree_B, stats_B)
               else if (ido == 2) then !Perform  y <--- M*x
                  call BPACK_Mult('N', Nunk_loc, 1, workd(ipntr(1)), workd(ipntr(2)), bmat_B, ptree_B, option_B, stats_B)
               endif
            else ! shift-invert mode
               bmattype = 'G'
               ! which='LM'        ! eigenvalues closest to the shifts
               iparam(7) = 3
               sigma = shift
               call pznaupd(ptree_A%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)
               if (ido == -1) then !Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x
                  call BPACK_Mult('N', Nunk_loc, 1, workd(ipntr(1)), workd(ipntr(2) + Nunk_loc), bmat_B, ptree_B, option_B, stats_B)
                  call BPACK_Solution(bmat_sh, workd(ipntr(2)), workd(ipntr(2) + Nunk_loc), Nunk_loc, 1, option_sh, ptree_sh, stats_sh)
               else if (ido == 1) then !Perform y <-- OP*x = inv[A-sigma*M]*M*x, M*x has been saved in workd(ipntr(3))
                  call BPACK_Solution(bmat_sh, workd(ipntr(2)), workd(ipntr(3)), Nunk_loc, 1, option_sh, ptree_sh, stats_sh)
               else if (ido == 2) then !Perform  y <--- M*x
                  call BPACK_Mult('N', Nunk_loc, 1, workd(ipntr(1)), workd(ipntr(2)), bmat_B, ptree_B, option_B, stats_B)
               endif
            endif
         enddo

         if (info < 0) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd.'
            print *, ' '
         else
            rvec = .true.
            call pzneupd(ptree_A%Comm, rvec, 'A', select, d, v, ldv, sigma, workev, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
            if (ierr /= 0) then
               print *, ' '
               print *, ' Error with _neupd, info = ', ierr
               print *, ' Check the documentation of _neupd. '
               print *, ' '
            else
               nconv = iparam(5)
               do j = 1, nconv

                  call BPACK_Mult('N', Nunk_loc, 1, v(1, j), ax, bmat_A, ptree_A, option_A, stats_A)
                  call BPACK_Mult('N', Nunk_loc, 1, v(1, j), mx, bmat_B, ptree_B, option_B, stats_B)
                  call zaxpy(Nunk_loc, -d(j), mx, 1, ax, 1)
                  rd(j, 1) = dble(d(j))
                  rd(j, 2) = dimag(d(j))
                  rd(j, 3) = pdznorm2(ptree_A%Comm, Nunk_loc, ax, 1)
                  rd(j, 3) = rd(j, 3)/dlapy2(rd(j, 1), rd(j, 2))
               enddo
               if (option_A%verbosity>=0) call pdmout(ptree_A%Comm, 6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and direct residuals')
               eigval(1:nconv) = d(1:nconv)
               eigvec(1:Nunk_loc, 1:nconv) = v(1:Nunk_loc, 1:nconv)
            end if

            if (ptree_A%MyID == Main_ID .and. option_A%verbosity>=0) then
               if (info .eq. 1) then
                  print *, ' '
                  print *, ' Maximum number of iterations reached.'
                  print *, ' '
               else if (info .eq. 3) then
                  print *, ' '
                  print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
                  print *, ' '
               end if
               print *, ' '
               print *, '_NDRV1'
               print *, '====== '
               print *, ' '
               print *, ' Size of the matrix is ', n
               print *, ' The number of processors is ', ptree_A%nproc
               print *, ' The number of Ritz values requested is ', nev
               print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
               print *, ' What portion of the spectrum: ', which
               print *, ' The number of converged Ritz values is ', nconv
               print *, ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
               print *, ' The number of OP*x is ', iparam(9)
               print *, ' The convergence criterion is ', tol
               print *, ' '
            endif
         end if

      endif
      t2 = MPI_Wtime()
      if (ptree_A%MyID == Main_ID .and. option_A%verbosity >= 0) write (*, *) 'Eigen Solve time: ', t2 - t1

      deallocate (select)
      deallocate (ax)
      deallocate (mx)
      deallocate (d)
      deallocate (v)
      deallocate (workd)
      deallocate (workev)
      deallocate (resid)
      deallocate (workl)
      deallocate (rwork)
      deallocate (rd)

#else
      if (ptree_A%MyID == Main_ID) write (*, *) 'eigen driver for real matrices are not yet implemented'
#endif

#else
      if (ptree_A%MyID == Main_ID) write (*, *) 'arpack not found, returning nconv=0'
#endif

   end subroutine BPACK_Eigen



   subroutine BPACK_Convert2Dense(bmat,option,stats,msh,ker,ptree)
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer num_vect,Nloc,ii
      DT,allocatable :: RandVectIn(:, :), RandVectOut(:, :), mat2D(:,:)
      type(matrixblock)::block_dummy
      integer pgno
      integer tempi, ctxt, info, iproc, jproc, myi, myj, myArows, myAcols, myrow, mycol, nprow, npcol, M, N, mnmin
      integer::descsMat1D(9), descsMat2D(9)


      ! construct a dummy block for auxiliary purposes
      pgno=1
      block_dummy%level = 0
      block_dummy%row_group = 1
      block_dummy%col_group = 1
      block_dummy%pgno = 1
      block_dummy%M = msh%Nunk
      block_dummy%N = msh%Nunk
      block_dummy%headm = 1
      block_dummy%headn = 1
      call ComputeParallelIndices(block_dummy, block_dummy%pgno, ptree, msh)


      num_vect = msh%Nunk
      Nloc = msh%idxe - msh%idxs + 1
      allocate (RandVectIn(Nloc, num_vect))
      RandVectIn=0
      allocate (RandVectOut(Nloc, num_vect))
      RandVectOut=0
      do ii=1,msh%Nunk
         if(ii>=msh%idxs .and. ii<=msh%idxe)then
            RandVectIn(ii-msh%idxs+1,ii)=1d0
         endif
      enddo
      call BPACK_Mult('N', Nloc, num_vect, RandVectIn, RandVectOut, bmat, ptree, option, stats)


      !!!!>**** generate 2D grid blacs quantities
      ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(block_dummy%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(block_dummy%N, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMat2D, block_dummy%M, block_dummy%N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (mat2D(max(1,myArows), max(1,myAcols)))
         mat2D = 0
      else
         descsMat2D(2) = -1
         allocate (mat2D(1, 1))
         mat2D = 0
      endif

      !!!!>**** redistribution of input matrix
      call Redistribute1Dto2D(RandVectOut, block_dummy%M_p, 0, pgno, mat2D, block_dummy%M, 0, pgno, block_dummy%N, ptree)
      deallocate(RandVectIn)
      deallocate(RandVectOut)
      call BF_delete(block_dummy, 1)

      if (associated(bmat%h_mat)) then
         allocate(bmat%h_mat%fullmat2D(size(mat2D,1),size(mat2D,2)))
         bmat%h_mat%fullmat2D=mat2D
      endif
      if (associated(bmat%ho_bf)) then
         allocate(bmat%ho_bf%fullmat2D(size(mat2D,1),size(mat2D,2)))
         bmat%ho_bf%fullmat2D=mat2D
      endif
      if (associated(bmat%hss_bf)) then
         allocate(bmat%hss_bf%fullmat2D(size(mat2D,1),size(mat2D,2)))
         bmat%hss_bf%fullmat2D=mat2D
      endif
      deallocate(mat2D)

   end subroutine BPACK_Convert2Dense


   subroutine BPACK_Eigen_Dense(bmat,option,stats,msh,ker,ptree,eigval,eigvec)
      implicit none
      DTC eigval(:)
      DT eigvec(:, :)
      DT,allocatable:: eigvec2d(:, :)
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer num_vect,Nloc,ii
      DT,allocatable :: RandVectIn(:, :), RandVectOut(:, :), mat2D(:,:)
      type(matrixblock)::block_dummy
      integer pgno
      integer tempi, ctxt, info, iproc, jproc, myi, myj, myArows, myAcols, myrow, mycol, nprow, npcol, M, N, mnmin
      integer::descsMat1D(9), descsMat2D(9)

      ! construct a dummy block for auxiliary purposes
      pgno=1
      block_dummy%level = 0
      block_dummy%row_group = 1
      block_dummy%col_group = 1
      block_dummy%pgno = 1
      block_dummy%M = msh%Nunk
      block_dummy%N = msh%Nunk
      block_dummy%headm = 1
      block_dummy%headn = 1
      call ComputeParallelIndices(block_dummy, block_dummy%pgno, ptree, msh)

      !!!!>**** generate 2D grid blacs quantities
      ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(block_dummy%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(block_dummy%N, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMat2D, block_dummy%M, block_dummy%N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (mat2D(max(1,myArows), max(1,myAcols)))
         if (associated(bmat%h_mat)) then
            mat2D=bmat%h_mat%fullmat2D
         endif
         if (associated(bmat%ho_bf)) then
            mat2D=bmat%ho_bf%fullmat2D
         endif
         if (associated(bmat%hss_bf)) then
            mat2D=bmat%hss_bf%fullmat2D
         endif
         allocate(eigvec2d(max(1,myArows), max(1,myAcols)))
      else
         descsMat2D(2) = -1
         allocate (mat2D(1, 1))
         mat2D = 0
         allocate(eigvec2d(1,1))
      endif

      call pgeeigf90(mat2D, block_dummy%M, descsMat2D, eigval, eigvec2d, ptree%pgrp(pgno)%ctxt, ptree%pgrp(pgno)%ctxt_head)

      !!!!>**** redistribution of input matrix
      call Redistribute2Dto1D(eigvec2d, block_dummy%M, 0, pgno, eigvec, block_dummy%M_p, 0, pgno, block_dummy%N, ptree)
	   deallocate(eigvec2d)
      deallocate(mat2D)
      call BF_delete(block_dummy, 1)


   end subroutine BPACK_Eigen_Dense



   subroutine BPACK_Solution(bmat, x, b, Ns_loc, num_vectors, option, ptree, stats)


      implicit none

      integer i, j, ii, jj, iii, jjj
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, Ns_loc, num_vectors
      real(kind=8) theta, phi, dphi, rcs_V, rcs_H
      real T0
      real(kind=8):: rel_error
      type(Hoption)::option
      type(proctree)::ptree
      type(Hstat)::stats
      ! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
      type(Bmatrix)::bmat
      DT::x(Ns_loc, num_vectors), b(Ns_loc, num_vectors)
      DT, allocatable::r0_initial(:)
      real(kind=8) n1, n2, rtemp

      n1 = MPI_Wtime()

      if (option%precon /= DIRECT) then
         allocate (r0_initial(1:Ns_loc))
         do ii = 1, Ns_loc
            call random_dp_number(r0_initial(ii))
         end do

         do ii = 1, num_vectors
            iter = 0
            rel_error = option%tol_itersol
            call BPACK_Ztfqmr(option%precon, option%n_iter, Ns_loc, b(:, ii), x(:, ii), rel_error, iter, r0_initial, bmat, ptree, option, stats)
         end do

         deallocate (r0_initial)
      else
         call BPACK_Inv_Mult('N', Ns_loc, num_vectors, b, x, bmat, ptree, option, stats)
         stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      end if

      n2 = MPI_Wtime()
      stats%Time_Sol = stats%Time_Sol + n2 - n1

      return

   end subroutine BPACK_Solution


   !!!!>**** expose the TFQMR as an API using user-supplied matvec, the TFQMR is the same as BPACK_Ztfqmr
   ! subroutine BPACK_Ztfqmr_usermatvec_noprecon(ntotal, nn_loc, b, x, err, iter, r0_initial, blackbox_MVP, ptree, option, stats, ker)
   !    implicit none
   !    integer level_c, rowblock, ierr
   !    integer, intent(in)::ntotal
   !    integer::iter, itmax, it, nn_loc
   !    DT, dimension(1:nn_loc,1)::x, bb, b, ytmp
   !    real(kind=8)::err, rerr
   !    DT, dimension(1:nn_loc,1)::w, yo, ayo, ye, aye, r, d, v
   !    real(kind=8)::ta, we, cm
   !    DT::we_local, we_sum, rerr_local, rerr_sum, err_local, err_sum
   !    DT::ta_local, ta_sum, bmag_local, bmag_sum1, dumb_ali(6)
   !    DT::etha, rho, rho_local, amgis, amgis_local, ahpla, dum, dum_local, beta
   !    real(kind=8)::bmag
   !    real(kind=8)::tim1, tim2
   !    integer::kk, ll
   !    ! Variables for storing current
   !    integer::count_unk, srcbox_id
   !    DT, dimension(:), allocatable::curr_coefs_dum_local, curr_coefs_dum_global
   !    real(kind=8)::mem_est
   !    character:: trans
   !    type(Hstat)::stats
   !    type(kernelquant)::ker
   !    procedure(HMatVec)::blackbox_MVP

   !    DT::r0_initial(:,:)
   !    type(proctree)::ptree
   !    type(Hoption)::option


   !    if(myisnan(sum(abs(x))))then
   !       write(*,*)'In BPACK_Ztfqmr, an initial guess of x is needed'
   !       stop
   !    endif

   !    itmax = iter

   !    ! call BPACK_ApplyPrecon(precond, nn_loc, b, bb, ptree, bmat, option, stats)
   !    bb=b

   !    ! ! ! if (myid == main_id) then
   !    ! ! ! call cpu_time(tim1)
   !    ! ! ! open(unit=32,file='iterations.out',status='unknown')
   !    ! ! ! end if

   !    if (iter .eq. 0) itmax = ntotal

   !    !  set initial values
   !    !
   !    d = 0d0
   !    ! write(*,*)'1'
   !    ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
   !    call blackbox_MVP('N',nn_loc,nn_loc,1,x,ytmp,ker)
   !    stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !    ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)
   !    r=ytmp

   !    r = bb - r !residual from the initial guess
   !    w = r
   !    yo = r
   !    ! ! write(*,*)'2'
   !    ! ! if(myisnan(sum(abs(yo)**2)))then
   !    ! ! write(*,*)'shitddd'
   !    ! ! stop
   !    ! ! end if
   !    ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
   !    call blackbox_MVP('N',nn_loc,nn_loc,1,yo,ytmp,ker)
   !    stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !    ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, ayo, ptree, bmat, option, stats)
   !    ayo = ytmp

   !    v = ayo
   !    we = 0d0
   !    etha = 0d0

   !    ta_local = dot_product(r(:,1), r(:,1))
   !    rho_local = dot_product(r0_initial(:,1), r(:,1))
   !    bmag_local = dot_product(bb(:,1), bb(:,1))

   !    dumb_ali(1:3) = (/ta_local, rho_local, bmag_local/)
   !    call MPI_ALLREDUCE(dumb_ali(1:3), dumb_ali(4:6), 3, MPI_DT, &
   !                       MPI_SUM, ptree%Comm, ierr)
   !    ta_sum = dumb_ali(4); rho = dumb_ali(5); bmag_sum1 = dumb_ali(6)
   !    ta = sqrt(abs(ta_sum))
   !    bmag = sqrt(abs(bmag_sum1))
   !    rerr = ta/bmag

   !    iters: do it = 1, itmax
   !       amgis_local = dot_product(r0_initial(:,1), v(:,1))
   !       call MPI_ALLREDUCE(amgis_local, amgis, 1, MPI_DT, MPI_SUM, &
   !                          ptree%Comm, ierr)
   !       ahpla = rho/amgis
   !       ye = yo - ahpla*v
   !       ! write(*,*)'3'
   !       ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,ye,aye)
   !       call blackbox_MVP('N',nn_loc,nn_loc,1,ye,ytmp,ker)
   !       stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !       ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, aye, ptree, bmat, option, stats)
   !       aye=ytmp

   !       !  start odd (2n-1) m loop
   !       d = yo + (we*we*etha/ahpla)*d
   !       w = w - ahpla*ayo
   !       we_local = dot_product(w(:,1), w(:,1))
   !       call MPI_ALLREDUCE(we_local, we_sum, 1, MPI_DT, MPI_SUM, &
   !                          ptree%Comm, ierr)
   !       we = sqrt(abs(we_sum))/ta
   !       cm = 1.0d0/sqrt(1.0d0 + we*we)
   !       ta = ta*we*cm
   !       etha = ahpla*cm*cm
   !       x = x + etha*d
   !       !  check if the result has converged.
   !       !a        if (err*bmag .gt. ta*sqrt(2.*it)) then
   !       !
   !       !  start even (2n)  m loop
   !       d = ye + (we*we*etha/ahpla)*d
   !       w = w - ahpla*aye
   !       we_local = dot_product(w(:,1), w(:,1))
   !       call MPI_ALLREDUCE(we_local, we_sum, 1, MPI_DT, MPI_SUM, &
   !                          ptree%Comm, ierr)
   !       we = sqrt(abs(we_sum))/ta
   !       cm = 1.0d0/sqrt(1.0d0 + we*we)
   !       ta = ta*we*cm
   !       etha = ahpla*cm*cm
   !       x = x + etha*d

   !       !  check if the result has converged.
   !       if (mod(it, 1) == 0 .or. rerr < 1d0*err) then
   !          ! write(*,*)'4'
   !          ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
   !          call blackbox_MVP('N',nn_loc,nn_loc,1,x,ytmp,ker)
   !          stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !          ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)
   !          r=ytmp

   !          r = bb - r
   !          rerr_local = dot_product(r(:,1), r(:,1))
   !          call MPI_ALLREDUCE(rerr_local, rerr_sum, 1, MPI_DT, MPI_SUM, &
   !                             ptree%Comm, ierr)
   !          rerr = sqrt(abs(rerr_sum))/bmag

   !          if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
   !             print *, '# ofiter,error:', it, rerr
   !             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
   !          end if

   !          if (err > rerr) then
   !             err = rerr
   !             iter = it

   !             ! ! ! if (myid == main_id) then
   !             ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
   !             ! ! ! end if

   !             return
   !          endif
   !       end if
   !       !  make preparations for next iteration
   !       dum_local = dot_product(r0_initial(:,1), w(:,1))
   !       call MPI_ALLREDUCE(dum_local, dum, 1, MPI_DT, MPI_SUM, &
   !                          ptree%Comm, ierr)
   !       beta = dum/rho
   !       rho = dum
   !       yo = w + beta*ye
   !       ! write(*,*)'5'
   !       ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
   !       call blackbox_MVP('N',nn_loc,nn_loc,1,yo,ytmp,ker)
   !       stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !       ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, ayo, ptree, bmat, option, stats)
   !       ayo=ytmp

   !       !MAGIC
   !       v = ayo + beta*(aye + beta*v)
   !    enddo iters
   !    ! write(*,*)'6'
   !    ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
   !    call blackbox_MVP('N',nn_loc,nn_loc,1,x,ytmp,ker)
   !    stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
   !    ! call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)
   !    r=ytmp

   !    !MAGIC
   !    r = bb - r
   !    err_local = dot_product(r(:,1), r(:,1))
   !    call MPI_ALLREDUCE(err_local, err_sum, 1, MPI_DT, MPI_SUM, ptree%Comm, &
   !                       ierr)
   !    err = sqrt(abs(err_sum))/bmag
   !    iter = itmax

   !    print *, 'Iterative solver is terminated without convergence!!!', it, err
   !    stop

   !    return
   ! end subroutine BPACK_Ztfqmr_usermatvec_noprecon


   !!!!>**** expose the TFQMR as an API using user-supplied matvec, the TFQMR is translated from matlab's tfqmr implementaion and slightly different compared to BPACK_Ztfqmr
   subroutine BPACK_Ztfqmr_usermatvec_noprecon(ntotal, nn_loc, b, x, err, iter, r0_initial, blackbox_MVP, ptree, option, stats, ker)
      implicit none
      integer ierr
      integer, intent(in)::ntotal
      integer::iter, itmax, it, nn_loc
      DT, dimension(1:nn_loc,1)::x, b, ytmp, r, r0, u_m, w, pu_m, v, Au, d, u_mp1, Ad, Au_new
      DT, allocatable:: eye(:,:),fullmat(:,:),UU(:,:),VV(:,:)
      DTR, allocatable::Singular(:)
      real(kind=8)::err
      integer::mm
      type(Hstat)::stats
      type(kernelquant)::ker
      procedure(HMatVec)::blackbox_MVP
      DT::r0_initial(:,:) ! not used
      type(proctree)::ptree
      type(Hoption)::option
      DT::nb_local,nr_local,nw_local,eta,alpha,r0v,r0w,sigma,beta,rho,rho_old
      integer flag,even
      real(kind=8)::tolb,nb,nr,nw,nr_act,tau,theta,c_mp1


      !!!!!!!!! check condition number of the operator
      ! if (ptree%nproc==1)then
      ! allocate(eye(nn_loc,nn_loc))
      ! allocate(fullmat(nn_loc,nn_loc))
      ! allocate(UU(nn_loc,nn_loc))
      ! allocate(VV(nn_loc,nn_loc))
      ! allocate(Singular(nn_loc))
      ! eye=0
      ! do mm=1,nn_loc
      !    eye(mm,mm)=1
      ! enddo
      ! fullmat=0
      ! call blackbox_MVP('N',nn_loc,nn_loc,nn_loc,eye,fullmat,ker)

      ! UU = 0
      ! VV = 0
      ! Singular = 0
      ! call gesvd_robust(fullmat, Singular, UU, VV, nn_loc, nn_loc, nn_loc)
      ! write(*,*)'condition number:', Singular(1)/Singular(nn_loc)
      ! deallocate(eye)
      ! deallocate(fullmat)
      ! deallocate(UU)
      ! deallocate(VV)
      ! deallocate(Singular)
      ! endif

      if(myisnan(sum(abs(x))))then
         write(*,*)'In BPACK_Ztfqmr, an initial guess of x is needed'
         stop
      endif

      itmax = ntotal


      nb_local = dot_product(b(:,1), b(:,1))
      call MPI_ALLREDUCE(MPI_IN_PLACE, nb_local, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
      nb = sqrt(abs(nb_local))

      flag=1
      tolb=err*nb

      call blackbox_MVP('N',nn_loc,nn_loc,1,x,ytmp,ker)
      r = b-ytmp

      nr_local = dot_product(r(:,1), r(:,1))
      call MPI_ALLREDUCE(MPI_IN_PLACE, nr_local, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
      nr = sqrt(abs(nr_local))
      nr_act = nr

      if (nr <= tolb)then  ! Initial guess is a good enough solution
         flag = 0
         iter = 0
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            print *, '# ofiter,error:', iter, nr
         end if
         return
      endif

      r0 = r
      u_m = r
      w = r
      pu_m = u_m
      call blackbox_MVP('N',nn_loc,nn_loc,1,pu_m,v,ker)
      Au = v
      d = 0
      Ad = d

      tau = nr
      theta = 0
      eta = 0

      rho = dot_product(r(:,1), r(:,1))
      call MPI_ALLREDUCE(MPI_IN_PLACE, rho, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
      rho_old = rho

      even=1

      do mm=1,itmax
         if(even==1)then
            r0v = dot_product(r0(:,1), v(:,1))
            call MPI_ALLREDUCE(MPI_IN_PLACE, r0v, 1, MPI_DT, MPI_SUM, &
                                 ptree%Comm, ierr)
            alpha = rho/r0v
            u_mp1 = u_m - alpha*v
         endif
         w = w - alpha*Au
         sigma = (theta**2/alpha)*eta
         d = pu_m + sigma*d
         Ad = Au + sigma*Ad

         nw_local = dot_product(w(:,1), w(:,1))
         call MPI_ALLREDUCE(MPI_IN_PLACE, nw_local, 1, MPI_DT, MPI_SUM, &
                              ptree%Comm, ierr)
         nw = sqrt(abs(nw_local))
         theta = nw/tau

         c_mp1 = 1d0/sqrt(1+theta**2);
         tau = tau*theta*c_mp1;
         eta = c_mp1**2*alpha;

         x = x + eta*d
         r = r - eta*Ad
         nr_local = dot_product(r(:,1), r(:,1))
         call MPI_ALLREDUCE(MPI_IN_PLACE, nr_local, 1, MPI_DT, MPI_SUM, &
                              ptree%Comm, ierr)
         nr = sqrt(abs(nr_local))
         nr_act = nr;

         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            print *, '# ofiter,error:', mm, nr_act
         end if
         iter = mm

         ! check for convergence
         if (nr <= tolb)then
            call blackbox_MVP('N',nn_loc,nn_loc,1,x,ytmp,ker)
            r = b - ytmp

            nr_local = dot_product(r(:,1), r(:,1))
            call MPI_ALLREDUCE(MPI_IN_PLACE, nr_local, 1, MPI_DT, MPI_SUM, &
                                 ptree%Comm, ierr)
            nr_act = sqrt(abs(nr_local))
            if (nr_act <= tolb)then
                  flag = 0
                  exit
            endif
         endif

         if(even==0)then
            r0w = dot_product(r0(:,1), w(:,1))
            call MPI_ALLREDUCE(MPI_IN_PLACE, r0w, 1, MPI_DT, MPI_SUM, &
                              ptree%Comm, ierr)
            rho = r0w;
            beta = rho/rho_old
            rho_old = rho
            u_mp1 = w + beta*u_m
         endif

         pu_m = u_mp1;
         call blackbox_MVP('N',nn_loc,nn_loc,1,pu_m,Au_new,ker)

         if(even==0)then
            v = Au_new + beta*(Au+beta*v)
         endif

         Au = Au_new
         u_m = u_mp1
         even=1-even
      enddo

      if (ptree%MyID == Main_ID )then
      if(iter>=itmax)then
         print *, 'Iterative solver is terminated without convergence!!!', iter, nr_act
         stop
      endif
      endif

      return
   end subroutine BPACK_Ztfqmr_usermatvec_noprecon

   subroutine BPACK_Ztfqmr(precond, ntotal, nn_loc, b, x, err, iter, r0_initial, bmat, ptree, option, stats)
      implicit none
      integer level_c, rowblock, ierr
      integer, intent(in)::ntotal
      integer::iter, itmax, it, nn_loc
      DT, dimension(1:nn_loc)::x, bb, b, ytmp
      real(kind=8)::err, rerr
      DT, dimension(1:nn_loc)::w, yo, ayo, ye, aye, r, d, v
      real(kind=8)::ta, we, cm
      DT::we_local, we_sum, rerr_local, rerr_sum, err_local, err_sum
      DT::ta_local, ta_sum, bmag_local, bmag_sum1, dumb_ali(6)
      DT::etha, rho, rho_local, amgis, amgis_local, ahpla, dum, dum_local, beta
      real(kind=8)::bmag
      real(kind=8)::tim1, tim2
      integer::kk, ll
      ! Variables for storing current
      integer::count_unk, srcbox_id
      DT, dimension(:), allocatable::curr_coefs_dum_local, curr_coefs_dum_global
      real(kind=8)::mem_est
      character:: trans
      type(Hstat)::stats

      ! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
      type(Bmatrix)::bmat
      DT::r0_initial(:)
      integer precond
      type(proctree)::ptree
      type(Hoption)::option


      if(myisnan(sum(abs(x))))then
         write(*,*)'In BPACK_Ztfqmr, an initial guess of x is needed'
         stop
      endif

      itmax = iter

      call BPACK_ApplyPrecon(precond, nn_loc, b, bb, ptree, bmat, option, stats)

      ! ! ! if (myid == main_id) then
      ! ! ! call cpu_time(tim1)
      ! ! ! open(unit=32,file='iterations.out',status='unknown')
      ! ! ! end if

      if (iter .eq. 0) itmax = ntotal

      !  set initial values
      !
      d = 0d0
      ! write(*,*)'1'
      ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
      call BPACK_Mult('N', nn_loc, 1, x, ytmp, bmat, ptree, option, stats)
      stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)

      r = bb - r !residual from the initial guess
      w = r
      yo = r
      ! ! write(*,*)'2'
      ! ! if(myisnan(sum(abs(yo)**2)))then
      ! ! write(*,*)'shitddd'
      ! ! stop
      ! ! end if
      ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
      call BPACK_Mult('N', nn_loc, 1, yo, ytmp, bmat, ptree, option, stats)
      stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      call BPACK_ApplyPrecon(precond, nn_loc, ytmp, ayo, ptree, bmat, option, stats)

      v = ayo
      we = 0d0
      etha = 0d0

      ta_local = dot_product(r, r)
      rho_local = dot_product(r0_initial, r)
      bmag_local = dot_product(bb, bb)

      dumb_ali(1:3) = (/ta_local, rho_local, bmag_local/)
      call MPI_ALLREDUCE(dumb_ali(1:3), dumb_ali(4:6), 3, MPI_DT, &
                         MPI_SUM, ptree%Comm, ierr)
      ta_sum = dumb_ali(4); rho = dumb_ali(5); bmag_sum1 = dumb_ali(6)
      ta = sqrt(abs(ta_sum))
      bmag = sqrt(abs(bmag_sum1))
      rerr = ta/bmag

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         print *, '# ofiter,error:', 0, rerr
         ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
      end if
      if (err > rerr) then
         err = rerr
         iter = 0
         return
      endif


      iters: do it = 1, itmax
         amgis_local = dot_product(r0_initial, v)
         call MPI_ALLREDUCE(amgis_local, amgis, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
         ahpla = rho/amgis
         ye = yo - ahpla*v
         ! write(*,*)'3'
         ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,ye,aye)
         call BPACK_Mult('N', nn_loc, 1, ye, ytmp, bmat, ptree, option, stats)
         stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
         call BPACK_ApplyPrecon(precond, nn_loc, ytmp, aye, ptree, bmat, option, stats)

         !  start odd (2n-1) m loop
         d = yo + (we*we*etha/ahpla)*d
         w = w - ahpla*ayo
         we_local = dot_product(w, w)
         call MPI_ALLREDUCE(we_local, we_sum, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
         we = sqrt(abs(we_sum))/ta
         cm = 1.0d0/sqrt(1.0d0 + we*we)
         ta = ta*we*cm
         etha = ahpla*cm*cm
         x = x + etha*d
         !  check if the result has converged.
         !a        if (err*bmag .gt. ta*sqrt(2.*it)) then
         !
         !  start even (2n)  m loop
         d = ye + (we*we*etha/ahpla)*d
         w = w - ahpla*aye
         we_local = dot_product(w, w)
         call MPI_ALLREDUCE(we_local, we_sum, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
         we = sqrt(abs(we_sum))/ta
         cm = 1.0d0/sqrt(1.0d0 + we*we)
         ta = ta*we*cm
         etha = ahpla*cm*cm
         x = x + etha*d

         !  check if the result has converged.
         if (mod(it, 1) == 0 .or. rerr < 1d0*err) then
            ! write(*,*)'4'
            ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
            call BPACK_Mult('N', nn_loc, 1, x, ytmp, bmat, ptree, option, stats)
            stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
            call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)

            r = bb - r
            rerr_local = dot_product(r, r)
            call MPI_ALLREDUCE(rerr_local, rerr_sum, 1, MPI_DT, MPI_SUM, &
                               ptree%Comm, ierr)
            rerr = sqrt(abs(rerr_sum))/bmag

            if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
               print *, '# ofiter,error:', it, rerr
               ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
            end if

            if (err > rerr) then
               err = rerr
               iter = it

               ! ! ! if (myid == main_id) then
               ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
               ! ! ! end if

               return
            endif
         end if
         !  make preparations for next iteration
         dum_local = dot_product(r0_initial, w)
         call MPI_ALLREDUCE(dum_local, dum, 1, MPI_DT, MPI_SUM, &
                            ptree%Comm, ierr)
         beta = dum/rho
         rho = dum
         yo = w + beta*ye
         ! write(*,*)'5'
         ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
         call BPACK_Mult('N', nn_loc, 1, yo, ytmp, bmat, ptree, option, stats)
         stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
         call BPACK_ApplyPrecon(precond, nn_loc, ytmp, ayo, ptree, bmat, option, stats)

         !MAGIC
         v = ayo + beta*(aye + beta*v)
      enddo iters
      ! write(*,*)'6'
      ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
      call BPACK_Mult('N', nn_loc, 1, x, ytmp, bmat, ptree, option, stats)
      stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      call BPACK_ApplyPrecon(precond, nn_loc, ytmp, r, ptree, bmat, option, stats)

      !MAGIC
      r = bb - r
      err_local = dot_product(r, r)
      call MPI_ALLREDUCE(err_local, err_sum, 1, MPI_DT, MPI_SUM, ptree%Comm, &
                         ierr)
      err = sqrt(abs(err_sum))/bmag
      iter = itmax

      print *, 'Iterative solver is terminated without convergence!!!', it, err
      stop

      return
   end subroutine BPACK_Ztfqmr

   subroutine BPACK_ApplyPrecon(precond, nn_loc, x, y, ptree, bmat, option, stats)
      implicit none
      integer nn_loc
      DT, dimension(1:nn_loc)::x, y
      integer precond
      type(Bmatrix)::bmat
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option

      if (precond == NOPRECON) then
         y = x
      else if (precond == HODLRPRECON) then
         call BPACK_Inv_Mult('N', nn_loc, 1, x, y, bmat, ptree, option, stats)
         stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      endif
   end subroutine BPACK_ApplyPrecon

   subroutine BPACK_Test_Solve_error(bmat, N_unk_loc, option, ptree, stats)



      implicit none

      integer i, j, ii, jj, iii, jjj, ierr
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
      real(kind=8) theta, phi, dphi, rcs_V, rcs_H
      real T0
      real(kind=8) n1, n2, rtemp
      DT value_Z
      DT, allocatable:: Voltage_pre(:), x(:, :), xtrue(:, :), b(:, :), btrue(:, :)
      real(kind=8):: rel_error, rtemp1, rtemp2, rtemp3, rtemp4, norm1, norm2, norm3, norm4
      type(Hoption)::option
      ! type(mesh)::msh
      ! type(kernelquant)::ker
      type(proctree)::ptree
      type(Bmatrix)::bmat
      type(Hstat)::stats
      DT, allocatable:: current(:), voltage(:)
      integer idxs, idxe

      allocate (x(N_unk_loc, 1))
      x = 0
      allocate (xtrue(N_unk_loc, 1))
      xtrue = 0
      call RandomMat(N_unk_loc, 1, 1, xtrue, 0)
      allocate (btrue(N_unk_loc, 1))
      btrue = 0
      allocate (b(N_unk_loc, 1))
      b = 0
      call BPACK_Mult('N', N_unk_loc, 1, xtrue, btrue, bmat, ptree, option, stats)
      stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      call BPACK_Solution(bmat, x, btrue, N_unk_loc, 1, option, ptree, stats)
      call BPACK_Mult('N', N_unk_loc, 1, x, b, bmat, ptree, option, stats)
      stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp

      rtemp1 = fnorm(xtrue - x, N_unk_loc, 1)**2d0;
      rtemp2 = fnorm(xtrue, N_unk_loc, 1)**2d0;
      call MPI_ALLREDUCE(rtemp1, norm1, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(rtemp2, norm2, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) '||X_t-H\(H*X_t)||_F/||X_t||_F: ', sqrt(norm1/norm2)
      endif

      rtemp3 = fnorm(btrue - b, N_unk_loc, 1)**2d0;
      rtemp4 = fnorm(btrue, N_unk_loc, 1)**2d0;
      call MPI_ALLREDUCE(rtemp3, norm3, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(rtemp4, norm4, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) '||B-H*(H\B)||_F/||B||_F: ', sqrt(norm3/norm4)
      endif

      deallocate (x)
      deallocate (xtrue)
      deallocate (btrue)
      deallocate (b)

   end subroutine BPACK_Test_Solve_error

   subroutine BPACK_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, bmat, ptree, option, stats)
      implicit none

      integer Ns
      character trans
      integer num_vectors
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(Bmatrix)::bmat
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option

      select case (option%format)
      case (HODLR)
         call HODLR_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, bmat%ho_bf, ptree, option, stats)
      case (HSS)
         call HSS_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, bmat%hss_bf, ptree, option, stats)
      case (HMAT)
         call Hmat_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, bmat%h_mat, ptree, option, stats)
      end select

   end subroutine BPACK_Inv_Mult

   subroutine Test_BPACK_Mult(Ns, bmat, ptree, option, stats)

      implicit none
      integer Ns
      DT::Vin(Ns, 1), Vout(Ns, 1)
      type(Bmatrix)::bmat
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::rtemp0, fnorm0
      integer ierr

      Vin = 1d0
      call BPACK_Mult('N', Ns, 1, Vin, Vout, bmat, ptree, option, stats)
      rtemp0 = fnorm(Vout, Ns, 1)**2d0
      call MPI_ALLREDUCE(rtemp0, fnorm0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, *) 'Test_BPACK_Mult: |Ax|_F: ', fnorm0

      Vin = 1d0
      call BPACK_Mult('T', Ns, 1, Vin, Vout, bmat, ptree, option, stats)
      rtemp0 = fnorm(Vout, Ns, 1)**2d0
      call MPI_ALLREDUCE(rtemp0, fnorm0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, *) 'Test_BPACK_Mult: |A^Tx|_F: ', fnorm0

   end subroutine Test_BPACK_Mult

   subroutine BPACK_Mult(trans, Ns, num_vectors, Vin, Vout, bmat, ptree, option, stats)

      implicit none
      character trans
      integer Ns
      integer num_vectors
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(Bmatrix)::bmat
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option

      select case (option%format)
      case (HODLR)
         call HODLR_Mult(trans, Ns, num_vectors, 1, bmat%ho_bf%Maxlevel + 1, Vin, Vout, bmat%ho_bf, ptree, option, stats)
      case (HMAT)
         call Hmat_Mult(trans, Ns, num_vectors, 1, bmat%h_mat%Maxlevel + 1, Vin, Vout, bmat%h_mat, ptree, option, stats,1)
      case (HSS)
         call HSS_Mult(trans, Ns, num_vectors, Vin, Vout, bmat%hss_bf, ptree, option, stats)
      end select

   end subroutine BPACK_Mult

   subroutine HODLR_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, ho_bf1, ptree, option, stats)



      implicit none

      integer Ns
      integer level_c, rowblock, head, tail
      integer i, j, k, level, ii, jj, kk, test, num_vectors, pp
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character trans, trans_tmp
      real(kind=8) a, b, c, d
      DT ctemp, ctemp1, ctemp2
      type(matrixblock), pointer::block_o

      ! type(vectorsblock), pointer :: random1, random2

      integer idx_start_glo, N_diag, idx_start_diag, idx_start_loc, idx_end_loc

      DT, allocatable::vec_old(:, :), vec_new(:, :)
      ! complex(kind=8)::Vin(:),Vout(:)
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer istart, iend, iinc

      idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1)

      stats%Flop_Tmp = 0

      trans_tmp = trans
      if (trans == 'C') then
         trans_tmp = 'T'
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      ! get the right multiplied vectors
      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      ! allocate(vec_old(Ns,num_vectors))
      allocate (vec_new(Ns, num_vectors))
      Vout = Vin
      ! write(*,*)'ddddd',Ns,num_vectors
      ! write(*,*)'begin'

      if (trans == 'N') then
         istart = ho_bf1%Maxlevel + 1
         iend = 1
         iinc = -1
      else
         istart = 1
         iend = ho_bf1%Maxlevel + 1
         iinc = 1
      endif

      do level = istart, iend, iinc
         vec_new = 0
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe
            pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%head + 1
            head = ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%headn - 1
            tail = head + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_loc - 1
            idx_start_loc = head - idx_start_glo + 1
            idx_end_loc = tail - idx_start_glo + 1

            if (level == ho_bf1%Maxlevel + 1) then
               call Full_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1), trans_tmp, idx_end_loc - idx_start_loc + 1, num_vectors,&
&Vout(idx_start_loc, 1), Ns, vec_new(idx_start_loc, 1), Ns, ctemp1, ctemp2)
               stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(idx_end_loc - idx_start_loc + 1, num_vectors, idx_end_loc - idx_start_loc + 1)
            else
               call Bplus_block_MVP_inverse_dat(ho_bf1, level, ii, trans_tmp, idx_end_loc - idx_start_loc + 1, num_vectors, Vout(idx_start_loc, 1), Ns, vec_new(idx_start_loc, 1), Ns, ptree, stats)

            endif
         end do
         Vout = vec_new
      end do
      ! Vout = vec_new(1:Ns,1)
      ! deallocate(vec_old)
      deallocate (vec_new)

      ! do ii=1,Ns
      ! write(131,*)abs(Vout(ii,1))
      ! enddo

      if (trans == 'C') then
         Vout = conjg(cmplx(Vout, kind=8))
         Vin = conjg(cmplx(Vin, kind=8))
      endif
      Vout = Vout*option%scale_factor

      return

   end subroutine HODLR_Inv_Mult

   subroutine HSS_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, hss_bf1, ptree, option, stats)



      implicit none

      integer Ns
      integer level_c, rowblock, head, tail
      integer i, j, k, level, ii, jj, kk, test, num_vectors, pp
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character trans, trans_tmp
      real(kind=8) a, b, c, d
      DT ctemp, ctemp1, ctemp2
      type(matrixblock), pointer::block_o

      ! type(vectorsblock), pointer :: random1, random2

      integer idx_start_glo, N_diag, idx_start_diag, idx_start_loc, idx_end_loc

      DT, allocatable::vec_old(:, :), vec_new(:, :)
      ! complex(kind=8)::Vin(:),Vout(:)
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hssbf)::hss_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer istart, iend, iinc

      stats%Flop_Tmp = 0

      trans_tmp = trans
      if (trans == 'C') then
         trans_tmp = 'T'
         Vin = conjg(cmplx(Vin, kind=8))
      endif
      Vout=0
      call Bplus_block_MVP_dat(hss_bf1%BP_inverse, trans, Ns, Ns, num_vectors, Vin, Ns, Vout, Ns, BPACK_cone, BPACK_czero, ptree, stats)

      if (trans == 'C') then
         Vout = conjg(cmplx(Vout, kind=8))
         Vin = conjg(cmplx(Vin, kind=8))
      endif
      Vout = Vout*option%scale_factor

      return

   end subroutine HSS_Inv_Mult

   subroutine HODLR_Mult(trans, Ns, num_vectors, level_start, level_end, Vin, Vout, ho_bf1, ptree, option, stats)

      implicit none

      character trans, trans_tmp
      integer Ns, level_start, level_end
      integer level_c, rowblock
      integer i, j, k, level, ii, jj, kk, test, num_vectors
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character chara
      real(kind=8) a, b, c, d
      DT ctemp, ctemp1, ctemp2
      ! type(matrixblock),pointer::block_o
      type(blockplus), pointer::bplus_o
      type(proctree)::ptree
      ! type(vectorsblock), pointer :: random1, random2
      type(Hstat)::stats
      type(Hoption)::option

      integer idx_start_glo, N_diag, idx_start_diag, idx_start_m, idx_end_m, idx_start_n, idx_end_n, pp, head, tail, idx_start_loc, idx_end_loc

      DT, allocatable::vec_old(:, :), vec_new(:, :)
      ! complex(kind=8)::Vin(:,:),Vout(:,:)
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hobf)::ho_bf1

      idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1)

      trans_tmp = trans
      if (trans == 'C') then
         trans_tmp = 'T'
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      ! get the right multiplied vectors
      ctemp1 = 1.0d0; ctemp2 = 1.0d0
      ! allocate(vec_old(Ns,num_vectors))
      allocate (vec_new(Ns, num_vectors))
      ! vec_old(1:Ns,1:num_vectors) = Vin
      vec_new = 0
      stats%Flop_Tmp = 0

      do level = level_start, level_end !ho_bf1%Maxlevel+1
         do ii = ho_bf1%levels(level)%Bidxs, ho_bf1%levels(level)%Bidxe

            pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%head + 1
            head = ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%headn - 1
            tail = head + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_loc - 1
            idx_start_loc = head - idx_start_glo + 1
            idx_end_loc = tail - idx_start_glo + 1

            if (level == ho_bf1%Maxlevel + 1) then
               call Full_block_MVP_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1), trans_tmp, idx_end_loc - idx_start_loc + 1, num_vectors,&
&Vin(idx_start_loc, 1), Ns, vec_new(idx_start_loc, 1), Ns, ctemp1, ctemp2)
               stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(idx_end_loc - idx_start_loc + 1, num_vectors, idx_end_loc - idx_start_loc + 1)
            else
               call Bplus_block_MVP_twoforward_dat(ho_bf1, level, ii, trans_tmp, idx_end_loc - idx_start_loc + 1, num_vectors, Vin(idx_start_loc, 1), Ns, vec_new(idx_start_loc, 1), Ns, ctemp1, ctemp2, ptree, stats)
            endif

         end do
      end do

      Vout = vec_new(1:Ns, 1:num_vectors)
      ! deallocate(vec_old)
      deallocate (vec_new)

      if (trans == 'C') then
         Vout = conjg(cmplx(Vout, kind=8))
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      Vout = Vout/option%scale_factor

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)"output norm: ",fnorm(Vout,Ns,num_vectors)**2d0

      return

   end subroutine HODLR_Mult



   !>***** redistribute the vector fed to Hmat from 1D to 2D layouts
   subroutine Hmat_Redistribute1Dto2D_Vector(Vin, Ns, num_vectors, vector2D, h_mat, ptree, nproc, stats, mode)

      implicit none
      integer Ns, num_vectors, i,i1, j, j1, gg, iproc, jproc, myrow, mycol, nprow, npcol, myi, myj, offr, offs, num_blocks
      integer ii, jj, ij, pp, tt
      integer level
      DT::Vin(Ns, num_vectors)
      type(vectorsblock):: vector2D(:)
      type(Hmat)::h_mat
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_sub, pgno_sub_mine, pid, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      integer::sendIDactive(nproc), recvIDactive(nproc)
      character::mode
      real(kind=8)::n1, n2
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist, kerflag
      integer::idxs_i,idxe_i, idxs_o,idxe_o

      n1 = MPI_Wtime()
      Ncol = num_vectors
      tag = 1

      ! allocation of communication quantities
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo

      Nsendactive = 0
      Nrecvactive = 0

      call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
      nprow = ptree%pgrp(1)%nprow
      npcol = ptree%pgrp(1)%npcol
      num_blocks = 2**h_mat%Dist_level
      do ii = 1, nproc
         idxs_i = h_mat%N_p(ii, 1)
         idxe_i = h_mat%N_p(ii, 2)
         do gg = 1, num_blocks
            idxs_o = h_mat%basis_group(gg)%head
            idxe_o = h_mat%basis_group(gg)%tail
            if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
               if(mode=='C')then
                  call g2l(gg, num_blocks, npcol, 1, jproc, myj)
                  do i1=1,nprow
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, i1-1, jproc)
                        pp = pid+1
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        sendquant(pp)%size = sendquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                     if(i1-1==myrow .and. jproc==mycol)then
                        pp=ii
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        recvquant(pp)%size = recvquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                  enddo
               else
                  call g2l(gg, num_blocks, nprow, 1, iproc, myi)
                  do j1=1,npcol
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, iproc, j1-1)
                        pp = pid+1
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        sendquant(pp)%size = sendquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                     if(iproc==myrow .and. j1-1==mycol)then
                        pp=ii
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        recvquant(pp)%size = recvquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo


      ! pack the send buffer in the second pass
      do ii = 1, nproc
         idxs_i = h_mat%N_p(ii, 1)
         idxe_i = h_mat%N_p(ii, 2)
         do gg = 1, num_blocks
            idxs_o = h_mat%basis_group(gg)%head
            idxe_o = h_mat%basis_group(gg)%tail
            if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
               if(mode=='C')then
                  call g2l(gg, num_blocks, npcol, 1, jproc, myj)
                  do i1=1,nprow
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, i1-1, jproc)
                        pp = pid+1

                        Nrow = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        offs = max(idxs_i, idxs_o) - idxs_i
                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = gg
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = max(idxs_i, idxs_o) - idxs_o
                        sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
                        sendquant(pp)%size = sendquant(pp)%size + 3
                        do i = 1, Nrow*Ncol
                           rr = mod(i - 1, Nrow) + 1
                           cc = (i - 1)/Nrow + 1
                           sendquant(pp)%dat(sendquant(pp)%size + i, 1) = Vin(offs+rr,cc)
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
                     endif
                  enddo
               else
                  call g2l(gg, num_blocks, nprow, 1, iproc, myi)
                  do j1=1,npcol
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, iproc, j1-1)
                        pp = pid+1

                        Nrow = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        offs = max(idxs_i, idxs_o) - idxs_i
                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = gg
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = max(idxs_i, idxs_o) - idxs_o
                        sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
                        sendquant(pp)%size = sendquant(pp)%size + 3
                        do i = 1, Nrow*Ncol
                           rr = mod(i - 1, Nrow) + 1
                           cc = (i - 1)/Nrow + 1
                           sendquant(pp)%dat(sendquant(pp)%size + i, 1) = Vin(offs+rr,cc)
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag, ptree%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag, ptree%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            gg = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            offr = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))

            if(mode=='C')then
               call g2l(gg, num_blocks, npcol, 1, jproc, myj)
               do cc = 1, Ncol
                  do rr = 1, Nrow
                     i = i + 1
                     vector2D(myj)%vector(offr+rr,cc) = recvquant(pp)%dat(i, 1)
                  enddo
               enddo
            else
               call g2l(gg, num_blocks, nprow, 1, iproc, myi)
               do cc = 1, Ncol
                  do rr = 1, Nrow
                     i = i + 1
                     vector2D(myi)%vector(offr+rr,cc) = recvquant(pp)%dat(i, 1)
                  enddo
               enddo
            endif
         enddo
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1


   end subroutine Hmat_Redistribute1Dto2D_Vector



   !>***** redistribute the vector fed to Hmat from 2D to 1D layouts
   subroutine Hmat_Redistribute2Dto1D_Vector(Vin, Ns, num_vectors, vector2D, h_mat, ptree, nproc, stats, mode)

      implicit none
      integer Ns, num_vectors, i,i1, j, j1, gg, iproc, jproc, myrow, mycol, nprow, npcol, myi, myj, offr, offs, num_blocks
      integer ii, jj, ij, pp, tt
      integer level
      DT::Vin(Ns, num_vectors)
      type(vectorsblock):: vector2D(:)
      type(Hmat)::h_mat
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_sub, pgno_sub_mine, pid, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      integer::sendIDactive(nproc), recvIDactive(nproc)
      character::mode
      real(kind=8)::n1, n2
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist, kerflag
      integer::idxs_i,idxe_i, idxs_o,idxe_o

      n1 = MPI_Wtime()
      Vin=0
      tag = 1
      Ncol = num_vectors

      ! allocation of communication quantities
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo

      Nsendactive = 0
      Nrecvactive = 0

      call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
      nprow = ptree%pgrp(1)%nprow
      npcol = ptree%pgrp(1)%npcol
      num_blocks = 2**h_mat%Dist_level
      do ii = 1, nproc
         idxs_o = h_mat%N_p(ii, 1)
         idxe_o = h_mat%N_p(ii, 2)
         do gg = 1, num_blocks
            idxs_i = h_mat%basis_group(gg)%head
            idxe_i = h_mat%basis_group(gg)%tail
            if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
               if(mode=='C')then
                  call g2l(gg, num_blocks, npcol, 1, jproc, myj)
                  do i1=1,nprow
                     if(i1-1==myrow .and. jproc==mycol)then
                        pp=ii
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        sendquant(pp)%size = sendquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, i1-1, jproc)
                        pp = pid+1
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        recvquant(pp)%size = recvquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                  enddo
               else
                  call g2l(gg, num_blocks, nprow, 1, iproc, myi)
                  do j1=1,npcol
                     if(iproc==myrow .and. j1-1==mycol)then
                        pp=ii
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        sendquant(pp)%size = sendquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                     if(ii==ptree%MyID+1)then
                        pid = blacs_pnum_wp(nprow,npcol, iproc, j1-1)
                        pp = pid+1
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        recvquant(pp)%size = recvquant(pp)%size + 3 + (min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1)*num_vectors
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo


      ! pack the send buffer in the second pass
      do ii = 1, nproc
         idxs_o = h_mat%N_p(ii, 1)
         idxe_o = h_mat%N_p(ii, 2)
         do gg = 1, num_blocks
            idxs_i = h_mat%basis_group(gg)%head
            idxe_i = h_mat%basis_group(gg)%tail
            if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
               if(mode=='C')then
                  call g2l(gg, num_blocks, npcol, 1, jproc, myj)
                  do i1=1,nprow
                     if(i1-1==myrow .and. jproc==mycol)then
                        pp=ii

                        Nrow = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        offs = max(idxs_i, idxs_o) - idxs_i
                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = gg
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = max(idxs_i, idxs_o) - idxs_o
                        sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
                        sendquant(pp)%size = sendquant(pp)%size + 3
                        do i = 1, Nrow*Ncol
                           rr = mod(i - 1, Nrow) + 1
                           cc = (i - 1)/Nrow + 1
                           sendquant(pp)%dat(sendquant(pp)%size + i, 1) = vector2D(myj)%vector(offs+rr,cc)
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
                     endif
                  enddo
               else
                  call g2l(gg, num_blocks, nprow, 1, iproc, myi)
                  do j1=1,npcol
                     if(iproc==myrow .and. j1-1==mycol)then
                        pp=ii

                        Nrow = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        offs = max(idxs_i, idxs_o) - idxs_i
                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = gg
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = max(idxs_i, idxs_o) - idxs_o
                        sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
                        sendquant(pp)%size = sendquant(pp)%size + 3
                        do i = 1, Nrow*Ncol
                           rr = mod(i - 1, Nrow) + 1
                           cc = (i - 1)/Nrow + 1
                           sendquant(pp)%dat(sendquant(pp)%size + i, 1) = vector2D(myi)%vector(offs+rr,cc)
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
                     endif
                  enddo
               endif
            endif
         enddo
      enddo

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag, ptree%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag, ptree%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            gg = NINT(dble(recvquant(pp)%dat(i, 1)))  ! this information is not needed for 2D-to-1D
            i = i + 1
            offr = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))

            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  Vin(offr+rr,cc) = Vin(offr+rr,cc) + recvquant(pp)%dat(i, 1)
               enddo
            enddo
         enddo
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1


   end subroutine Hmat_Redistribute2Dto1D_Vector



   subroutine HSS_Mult(trans, Ns, num_vectors, Vin, Vout, hss_bf1, ptree, option, stats)

      implicit none

      character trans, trans_tmp
      integer Ns
      integer level_c, rowblock
      integer i, j, k, level, ii, jj, kk, test, num_vectors
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character chara
      real(kind=8) a, b, c, d
      DT ctemp, ctemp1, ctemp2
      ! type(matrixblock),pointer::block_o
      type(blockplus), pointer::bplus_o
      type(proctree)::ptree
      ! type(vectorsblock), pointer :: random1, random2
      type(Hstat)::stats
      type(Hoption)::option

      integer idx_start_glo, N_diag, idx_start_diag, idx_start_m, idx_end_m, idx_start_n, idx_end_n, pp, head, tail, idx_start_loc, idx_end_loc

      DT, allocatable::vec_old(:, :), vec_new(:, :)
      ! complex(kind=8)::Vin(:,:),Vout(:,:)
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(hssbf)::hss_bf1

      trans_tmp = trans
      if (trans == 'C') then
         trans_tmp = 'T'
         Vin = conjg(cmplx(Vin, kind=8))
      endif
      Vout=0
      stats%Flop_Tmp = 0
      call Bplus_block_MVP_dat(hss_bf1%BP, trans, Ns, Ns, num_vectors, Vin, Ns, Vout, Ns, BPACK_cone, BPACK_czero, ptree, stats)

      if (trans == 'C') then
         Vout = conjg(cmplx(Vout, kind=8))
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      Vout = Vout/option%scale_factor

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)"output norm: ",fnorm(Vout,Ns,num_vectors)**2d0

      return

   end subroutine HSS_Mult

   subroutine Hmat_Inv_Mult(trans, Ns, num_vectors, Vin, Vout, h_mat, ptree, option, stats)
      implicit none

      integer Ns, ii
      character trans, trans_tmp
      integer num_vectors
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(Hmat)::h_mat
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option

      stats%Flop_Tmp = 0

      trans_tmp = trans
      if (trans == 'C') then
         trans_tmp = 'T'
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      Vout = Vin

      if (trans == 'N') then
         ! write(*,*)fnorm(Vout,Ns,num_vectors),'before L solve'
         call Hmat_Lsolve_Toplevel(h_mat, trans_tmp, Vout, Ns, num_vectors, ptree, stats)
         ! write(*,*)fnorm(Vout,Ns,num_vectors),'before U solve',abs(Vout)
         call Hmat_Usolve_Toplevel(h_mat, trans_tmp, Vout, Ns, num_vectors, ptree, stats)
         ! write(*,*)fnorm(Vout,Ns,num_vectors),'after LU solve'
         ! do ii=1,Ns
         ! write(130,*)abs(Vout(ii,1))
         ! enddo
      else
         call Hmat_Usolve_Toplevel(h_mat, trans_tmp, Vout, Ns, num_vectors, ptree, stats)
         call Hmat_Lsolve_Toplevel(h_mat, trans_tmp, Vout, Ns, num_vectors, ptree, stats)
      endif

      if (trans == 'C') then
         Vout = conjg(cmplx(Vout, kind=8))
         Vin = conjg(cmplx(Vin, kind=8))
      endif

      Vout = Vout*option%scale_factor

   end subroutine Hmat_Inv_Mult

   subroutine Hmat_Mult(trans, Ns, num_vectors, level_start, level_end, Vin, Vout, h_mat, ptree, option, stats, use_blockcopy)

      implicit none
      integer level_start, level_end
      character trans, trans_tmp
      integer Ns
      integer level_c, rowblock
      integer i, j, k, level, num_blocks, ii, jj, kk, test, num_vectors
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character chara
      real(kind=8) vecnorm, n1, n2
      DT ctemp, ctemp1, ctemp2
      ! type(matrixblock),pointer::block_o
      type(blockplus), pointer::bplus_o
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option

      integer use_blockcopy, idx_start_glo, N_diag, idx_start_diag, idx_start_m, idx_end_m, idx_start_n, idx_end_n, pp, head, tail, idx_start_loc, idx_end_loc, Nmax, Nmsg,N_glo
      type(matrixblock), pointer :: blocks_i, blocks_j
      type(matrixblock) :: blocks_dummy
      DT, allocatable::vec_old(:, :), vec_new(:, :), vin_tmp(:, :), vout_tmp(:, :), vec_buffer(:, :), Vout_glo(:, :), Vin_glo(:, :)

      ! DT::Vin(:,:),Vout(:,:)
      DT::Vin(Ns, num_vectors), Vout(Ns, num_vectors)
      type(Hmat)::h_mat
      integer ierr, m_size
      integer, allocatable::status_all(:, :), srequest_all(:)
      integer :: status(MPI_Status_size)
      type(vectorsblock),allocatable:: vector2D_i(:),vector2D_o(:)
      integer:: nprow, npcol, myrow, mycol
      character mode_i,mode_o

      if (ptree%Comm /= MPI_COMM_NULL) then

         trans_tmp = trans
         if (trans == 'C') then
            trans_tmp = 'T'
            Vin = conjg(cmplx(Vin, kind=8))
         endif

         call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
         num_blocks = 2**h_mat%Dist_level
         if (trans == 'N') then
            mode_i='C'
            mode_o='R'
            allocate(vector2D_i(max(h_mat%myAcols,1)))
            allocate(vector2D_o(max(h_mat%myArows,1)))
            if(h_mat%myAcols>0 .and. h_mat%myArows>0)then
            do j = 1, h_mat%myAcols
               call l2g(j, mycol, num_blocks, npcol, 1, jj)
               blocks_i => h_mat%Local_blocks(j, 1)
               allocate(vector2D_i(j)%vector(blocks_i%N,num_vectors))
               vector2D_i(j)%vector=0
            enddo
            do i = 1, h_mat%myArows
               call l2g(i, myrow, num_blocks, nprow, 1, ii)
               blocks_i => h_mat%Local_blocks(1, i)
               allocate(vector2D_o(i)%vector(blocks_i%M,num_vectors))
               vector2D_o(i)%vector=0
            enddo
            endif
         else
            mode_i='R'
            mode_o='C'
            allocate(vector2D_i(max(h_mat%myArows,1)))
            allocate(vector2D_o(max(h_mat%myAcols,1)))
            if(h_mat%myAcols>0 .and. h_mat%myArows>0)then
            do i = 1, h_mat%myArows
               call l2g(i, myrow, num_blocks, nprow, 1, ii)
               blocks_i => h_mat%Local_blocks(1, i)
               allocate(vector2D_i(i)%vector(blocks_i%M,num_vectors))
               vector2D_i(i)%vector=0
            enddo
            do j = 1, h_mat%myAcols
               call l2g(j, mycol, num_blocks, npcol, 1, jj)
               blocks_i => h_mat%Local_blocks(j, 1)
               allocate(vector2D_o(j)%vector(blocks_i%N,num_vectors))
               vector2D_o(j)%vector=0
            enddo
            endif
         endif


         call Hmat_Redistribute1Dto2D_Vector(Vin, Ns, num_vectors, vector2D_i, h_mat, ptree, ptree%nproc, stats, mode_i)

         ! call MPI_barrier(ptree%Comm, ierr)

         n1 = MPI_Wtime()

         do i = 1, h_mat%myArows
            do j = 1, h_mat%myAcols
               if(use_blockcopy==1)then
                  blocks_i => h_mat%Local_blocks_copy(j, i)
               else
                  blocks_i => h_mat%Local_blocks(j, i)
               endif
               if (trans == 'N') then
                  call Hmat_block_MVP_dat(blocks_i, trans_tmp, blocks_i%headm, blocks_i%headn, num_vectors, vector2D_i(j)%vector, blocks_i%N, vector2D_o(i)%vector, blocks_i%M, BPACK_cone, ptree, stats,level_start, level_end)
               else
                  call Hmat_block_MVP_dat(blocks_i, trans_tmp, blocks_i%headm, blocks_i%headn, num_vectors, vector2D_i(i)%vector, blocks_i%M, vector2D_o(j)%vector, blocks_i%N, BPACK_cone, ptree, stats,level_start, level_end)
               endif
            enddo
         enddo

         call Hmat_Redistribute2Dto1D_Vector(Vout, Ns, num_vectors, vector2D_o, h_mat, ptree, ptree%nproc, stats, mode_o)

         if (trans == 'N') then
            if(h_mat%myAcols>0 .and. h_mat%myArows>0)then
            do j = 1, h_mat%myAcols
               deallocate(vector2D_i(j)%vector)
            enddo
            do i = 1, h_mat%myArows
               deallocate(vector2D_o(i)%vector)
            enddo
            endif
            deallocate(vector2D_i)
            deallocate(vector2D_o)
         else
            if(h_mat%myAcols>0 .and. h_mat%myArows>0)then
            do i = 1, h_mat%myArows
               deallocate(vector2D_i(i)%vector)
            enddo
            do j = 1, h_mat%myAcols
               deallocate(vector2D_o(j)%vector)
            enddo
            endif
            deallocate(vector2D_i)
            deallocate(vector2D_o)
         endif


         if (trans == 'C') then
            Vout = conjg(cmplx(Vout, kind=8))
            Vin = conjg(cmplx(Vin, kind=8))
         endif

         Vout = Vout/option%scale_factor

         call MPI_barrier(ptree%Comm, ierr)
         n2 = MPI_Wtime()

         ! vecnorm = fnorm(Vout, Ns, num_vectors)**2d0
         ! call MPI_AllREDUCE(MPI_IN_PLACE, vecnorm, 1, MPI_double, MPI_SUM, ptree%Comm, ierr)
         ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)"output norm: ",sqrt(vecnorm)
      endif

      return

   end subroutine Hmat_Mult

   subroutine Hmat_Lsolve_Toplevel(h_mat, trans, xloc, nloc, nvec, ptree, stats)
      implicit none

      type(Hmat)::h_mat
      character::trans
      type(proctree)::ptree
      type(Hstat)::stats

      integer i,i1, j, k, ii, jj, kk, iii, jjj, pp, num_blocks, mm, nn, groupm, groupn
      integer vectors_start, vectors_x, vectors_y, id_l, nvec, nloc, tag

      type(matrixblock), pointer :: blocks_l
      integer Nreq, Nmod, Bufsize
      DT, allocatable::recv_buf(:), vin(:, :)
      DT::xloc(:, :)
      DT,allocatable::xglo(:,:)
      integer :: status(MPI_Status_size)
      integer, allocatable::status_all(:, :), request_all(:)
      integer idx_start,m_size
      integer ierr,num_vectors,selflag,offflag,nrecvx,nrecvmod,nprow, npcol, myrow, mycol,iproc,myi,jproc,myj,jproc1,myj1,receiver, nrecv, N_glo
      integer,allocatable:: fmod(:), frecv(:),sendflagarray(:), receiverlists(:,:)
      type(vectorsblock),allocatable:: sendbufx(:),sendbufmod(:)

      if (trans == 'N') then
         Nreq=0
         num_blocks = 2**h_mat%Dist_level
         call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
         nrecvx=0
         nrecvmod=0
         N_glo = h_mat%N
         allocate(xglo(N_glo,nvec))
         xglo=0
         xglo(h_mat%idxs:h_mat%idxe,:) = xloc(1:nloc,:)
         call MPI_ALLREDUCE(MPI_IN_PLACE, xglo, N_glo*nvec, MPI_DT, MPI_SUM, ptree%Comm, ierr)

         allocate(receiverlists(max(1,nprow),max(1,h_mat%myAcols)))
         receiverlists=0

         allocate(sendbufx(max(h_mat%myAcols,1)))
         allocate(sendbufmod(max(h_mat%myArows,1)))

         do j = 1, h_mat%myAcols
            selflag=0
            call l2g(j, mycol, num_blocks, npcol, 1, jj)
            do i = 1, h_mat%myArows
               blocks_l => h_mat%Local_blocks(j, i)
               if(blocks_l%row_group==blocks_l%col_group)then
                  selflag=1
                  allocate(sendbufx(j)%vector(blocks_l%N,nvec))
                  sendbufx(j)%vector = xglo(blocks_l%headn:blocks_l%headn+blocks_l%N-1,:)
                  do ii=1,num_blocks
                     if(ii>jj)then
                        call g2l(ii, num_blocks, nprow, 1, iproc, myi)
                        receiverlists(iproc+1,j)=1
                     endif
                  enddo
                  do ii=1,nprow
                     if(receiverlists(ii,j)==1)then
                        Nreq = Nreq+1
                     endif
                  enddo
               endif
            enddo

            offflag=0
            do i = 1, h_mat%myArows
               blocks_l => h_mat%Local_blocks(j, i)
               if(blocks_l%row_group>blocks_l%col_group)then
                  offflag=1
                  exit
               endif
            enddo
            if(offflag==1)then
               nrecvx = nrecvx + 1
            endif
         enddo

         allocate(fmod(max(h_mat%myArows,1)))
         fmod=0
         allocate(frecv(max(h_mat%myArows,1)))
         frecv=0
         allocate(sendflagarray(num_blocks))
         sendflagarray=0
         do i = 1, h_mat%myArows
            call l2g(i, myrow, num_blocks, nprow, 1, ii)
            do j = 1, h_mat%myAcols
               blocks_l => h_mat%Local_blocks(j, i)
               if(blocks_l%row_group>blocks_l%col_group)then
                  fmod(i) = fmod(i) + 1
                  sendflagarray(ii)=1
               endif
            enddo
            if(sendflagarray(ii)==1)then
               Nreq = Nreq+1
               blocks_l => h_mat%Local_blocks(1, i)
               allocate(sendbufmod(i)%vector(blocks_l%M,nvec))
               sendbufmod(i)%vector=0
            endif
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, sendflagarray, num_blocks, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
         do i = 1, h_mat%myArows
            call l2g(i, myrow, num_blocks, nprow, 1, ii)
            do j = 1, h_mat%myAcols
               blocks_l => h_mat%Local_blocks(j, i)
               if(blocks_l%row_group==blocks_l%col_group)then
                  nrecvmod = nrecvmod + sendflagarray(ii)
                  frecv(i) = frecv(i) + sendflagarray(ii)
               endif
            enddo
         enddo
         deallocate(sendflagarray)

         if (Nreq > 0) then
            allocate (status_all(MPI_status_size, Nreq))
            allocate (request_all(Nreq))
         end if
         Nreq=0

         ii=1
         jj=1
         call g2l(ii, num_blocks, nprow, 1, iproc, myi)
         call g2l(jj, num_blocks, npcol, 1, jproc, myj)

         if(iproc==myrow .and. jproc==mycol)then
            blocks_l => h_mat%Local_blocks(myj, myi)

            idx_start = blocks_l%headn
            call Hmat_Lsolve(blocks_l, 'N', idx_start, nvec, sendbufx(myj)%vector, blocks_l%N, ptree, stats)

            do i=1,nprow
               if(receiverlists(i,myj)==1)then
                  receiver = blacs_pnum_wp(nprow,npcol, i-1, jproc)
                  Nreq = Nreq + 1
                  call MPI_Isend(sendbufx(myj)%vector, blocks_l%N*nvec, MPI_DT, receiver, jj, ptree%Comm, request_all(Nreq), ierr)
               endif
            enddo
         endif
         do nrecv=1, nrecvx+nrecvmod
            call MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, ptree%Comm, status,ierr)
            pp = status(MPI_SOURCE)
            tag = status(MPI_TAG)
            call MPI_Get_count(status, MPI_DT, m_size,ierr)
            nn = m_size/nvec
            allocate (vin(nn, nvec))
            call MPI_Recv(vin, m_size, MPI_DT, pp, tag, ptree%Comm, status, ierr)

            if(tag<=num_blocks)then
               jj=tag
               call g2l(jj, num_blocks, npcol, 1, jproc, myj)
               do i = 1, h_mat%myArows
                  call l2g(i, myrow, num_blocks, nprow, 1, ii)
                  if(ii>jj)then
                     blocks_l => h_mat%Local_blocks(myj, i)
                     nn = blocks_l%N

                     call Hmat_block_MVP_dat(blocks_l, 'N', blocks_l%headm, blocks_l%headn, nvec, vin, nn, sendbufmod(i)%vector, blocks_l%M,-BPACK_cone, ptree, stats)
                     fmod(i) = fmod(i)-1

                     if(fmod(i)==0)then
                        call g2l(ii, num_blocks, npcol, 1, jproc1, myj1)

                        ! offdiagonal block: send right to diagonal block
                        receiver = blacs_pnum_wp(nprow,npcol, myrow, jproc1)
                        Nreq = Nreq + 1
                        call MPI_Isend(sendbufmod(i)%vector, blocks_l%M*nvec, MPI_DT, receiver, ii+num_blocks, ptree%Comm, request_all(Nreq), ierr)

                     endif
                  endif
               enddo
            else ! diagonal block: solve the diagonal block, send down to off-digonal blocks
               ii=tag-num_blocks
               call g2l(ii, num_blocks, nprow, 1, iproc, myi)
               call g2l(ii, num_blocks, npcol, 1, jproc1, myj1)
               sendbufx(myj1)%vector = sendbufx(myj1)%vector + vin
               frecv(myi) = frecv(myi)-1
               if(frecv(myi)==0 .and. fmod(myi)==0)then
                  blocks_l => h_mat%Local_blocks(myj1, myi)
                  idx_start = blocks_l%headn
                  call Hmat_Lsolve(blocks_l, 'N', idx_start, nvec, sendbufx(myj1)%vector, blocks_l%N, ptree, stats)

                  do i1=1,nprow
                     if(receiverlists(i1,myj1)==1)then
                        receiver = blacs_pnum_wp(nprow,npcol, i1-1, jproc1)
                        Nreq = Nreq + 1
                        call MPI_Isend(sendbufx(myj1)%vector, blocks_l%N*nvec, MPI_DT, receiver, ii, ptree%Comm, request_all(Nreq), ierr)
                     endif
                  enddo
               endif
            endif
            deallocate (vin)
         enddo

         if (Nreq > 0) then
            call MPI_waitall(Nreq, request_all, status_all, ierr)
            deallocate (status_all)
            deallocate (request_all)
         endif

         xglo=0
         do j = 1, h_mat%myAcols
            do i = 1, h_mat%myArows
               blocks_l => h_mat%Local_blocks(j, i)
               if(blocks_l%row_group==blocks_l%col_group)then
                  xglo(blocks_l%headn:blocks_l%headn+blocks_l%N-1,:) = sendbufx(j)%vector
               endif
            enddo
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, xglo, N_glo*nvec, MPI_DT, MPI_SUM, ptree%Comm, ierr)
         xloc(1:nloc,:) = xglo(h_mat%idxs:h_mat%idxe,:)
         deallocate(xglo)
         deallocate(fmod)
         deallocate(frecv)
         deallocate(receiverlists)

         do j=1,h_mat%myAcols
            if(allocated(sendbufx(j)%vector))deallocate(sendbufx(j)%vector)
         enddo
         deallocate(sendbufx)

         do i=1,h_mat%myArows
            if(allocated(sendbufmod(i)%vector))deallocate(sendbufmod(i)%vector)
         enddo
         deallocate(sendbufmod)
      else
         write (*, *) 'XxL^-1 with MPI is not yet implemented'
         stop
      endif

      call MPI_barrier(ptree%Comm, ierr)
      return

   end subroutine Hmat_Lsolve_Toplevel

   subroutine Hmat_Usolve_Toplevel(h_mat, trans, xloc, nloc, nvec, ptree, stats)

      implicit none

      type(Hmat)::h_mat
      character::trans
      type(proctree)::ptree
      type(Hstat)::stats

      integer i,i1, j, k, ii, jj, kk, pp, iii, jjj, num_blocks, mm, nn, idx_start, nvec, nloc
      integer vectors_start, vectors_x, vectors_y, id_u, groupm, groupn, ierr

      type(matrixblock), pointer :: blocks_u
      integer :: status(MPI_Status_size)
      integer Nreq, Nmod, Bufsize, m_size
      DT, allocatable::recv_buf(:), vin(:, :)
      DT::xloc(:, :)
      DT,allocatable::xglo(:,:)
      integer, allocatable::status_all(:, :), request_all(:)

      integer selflag,offflag,nrecvx,nrecvmod,nprow, npcol, myrow, mycol,iproc,myi,jproc,myj,jproc1,myj1,receiver, nrecv, N_glo,tag
      integer,allocatable:: bmod(:), brecv(:),sendflagarray(:), receiverlists(:,:)
      type(vectorsblock),allocatable:: sendbufx(:),sendbufmod(:)



      if (trans == 'N') then

         Nreq=0
         num_blocks = 2**h_mat%Dist_level
         call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
         nrecvx=0
         nrecvmod=0
         N_glo = h_mat%N
         allocate(xglo(N_glo,nvec))
         xglo=0
         xglo(h_mat%idxs:h_mat%idxe,:) = xloc(1:nloc,:)
         call MPI_ALLREDUCE(MPI_IN_PLACE, xglo, N_glo*nvec, MPI_DT, MPI_SUM, ptree%Comm, ierr)

         allocate(receiverlists(max(1,nprow),max(1,h_mat%myAcols)))
         receiverlists=0

         allocate(sendbufx(max(h_mat%myAcols,1)))
         allocate(sendbufmod(max(h_mat%myArows,1)))

         do j = 1, h_mat%myAcols
            selflag=0
            call l2g(j, mycol, num_blocks, npcol, 1, jj)
            do i = 1, h_mat%myArows
               blocks_u => h_mat%Local_blocks(j, i)
               if(blocks_u%row_group==blocks_u%col_group)then
                  selflag=1
                  allocate(sendbufx(j)%vector(blocks_u%N,nvec))
                  sendbufx(j)%vector = xglo(blocks_u%headn:blocks_u%headn+blocks_u%N-1,:)
                  do ii=1,num_blocks
                     if(ii<jj)then
                        call g2l(ii, num_blocks, nprow, 1, iproc, myi)
                        receiverlists(iproc+1,j)=1
                     endif
                  enddo
                  do ii=1,nprow
                     if(receiverlists(ii,j)==1)then
                        Nreq = Nreq+1
                     endif
                  enddo
               endif
            enddo

            offflag=0
            do i = 1, h_mat%myArows
               blocks_u => h_mat%Local_blocks(j, i)
               if(blocks_u%row_group<blocks_u%col_group)then
                  offflag=1
                  exit
               endif
            enddo
            if(offflag==1)then
               nrecvx = nrecvx + 1
            endif
         enddo

         allocate(bmod(max(h_mat%myArows,1)))
         bmod=0
         allocate(brecv(max(h_mat%myArows,1)))
         brecv=0
         allocate(sendflagarray(num_blocks))
         sendflagarray=0
         do i = 1, h_mat%myArows
            call l2g(i, myrow, num_blocks, nprow, 1, ii)
            do j = 1, h_mat%myAcols
               blocks_u => h_mat%Local_blocks(j, i)
               if(blocks_u%row_group<blocks_u%col_group)then
                  bmod(i) = bmod(i) + 1
                  sendflagarray(ii)=1
               endif
            enddo
            if(sendflagarray(ii)==1)then
               Nreq = Nreq+1
               blocks_u => h_mat%Local_blocks(1, i)
               allocate(sendbufmod(i)%vector(blocks_u%M,nvec))
               sendbufmod(i)%vector=0
            endif
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, sendflagarray, num_blocks, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
         do i = 1, h_mat%myArows
            call l2g(i, myrow, num_blocks, nprow, 1, ii)
            do j = 1, h_mat%myAcols
               blocks_u => h_mat%Local_blocks(j, i)
               if(blocks_u%row_group==blocks_u%col_group)then
                  nrecvmod = nrecvmod + sendflagarray(ii)
                  brecv(i) = brecv(i) + sendflagarray(ii)
               endif
            enddo
         enddo
         deallocate(sendflagarray)

         if (Nreq > 0) then
            allocate (status_all(MPI_status_size, Nreq))
            allocate (request_all(Nreq))
         end if
         Nreq=0

         ii=num_blocks
         jj=num_blocks
         call g2l(ii, num_blocks, nprow, 1, iproc, myi)
         call g2l(jj, num_blocks, npcol, 1, jproc, myj)

         if(iproc==myrow .and. jproc==mycol)then
            blocks_u => h_mat%Local_blocks(myj, myi)

            idx_start = blocks_u%headn
            call Hmat_Usolve(blocks_u, 'N', idx_start, nvec, sendbufx(myj)%vector, blocks_u%N, ptree, stats)

            do i=1,nprow
               if(receiverlists(i,myj)==1)then
                  receiver = blacs_pnum_wp(nprow,npcol, i-1, jproc)
                  Nreq = Nreq + 1
                  call MPI_Isend(sendbufx(myj)%vector, blocks_u%N*nvec, MPI_DT, receiver, jj, ptree%Comm, request_all(Nreq), ierr)
               endif
            enddo
         endif

         do nrecv=1, nrecvx+nrecvmod
            call MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, ptree%Comm, status,ierr)
            pp = status(MPI_SOURCE)
            tag = status(MPI_TAG)
            call MPI_Get_count(status, MPI_DT, m_size,ierr)
            nn = m_size/nvec
            allocate (vin(nn, nvec))
            call MPI_Recv(vin, m_size, MPI_DT, pp, tag, ptree%Comm, status, ierr)

            if(tag<=num_blocks)then
               jj=tag
               call g2l(jj, num_blocks, npcol, 1, jproc, myj)
               do i = 1, h_mat%myArows
                  call l2g(i, myrow, num_blocks, nprow, 1, ii)
                  if(ii<jj)then
                     blocks_u => h_mat%Local_blocks(myj, i)
                     nn = blocks_u%N


                     call Hmat_block_MVP_dat(blocks_u, 'N', blocks_u%headm, blocks_u%headn, nvec, vin, nn, sendbufmod(i)%vector, blocks_u%M,-BPACK_cone, ptree, stats)
                     bmod(i) = bmod(i)-1

                     if(bmod(i)==0)then
                        call g2l(ii, num_blocks, npcol, 1, jproc1, myj1)
                        ! offdiagonal block: send left to diagonal block
                        receiver = blacs_pnum_wp(nprow,npcol, myrow, jproc1)
                        Nreq = Nreq + 1
                        call MPI_Isend(sendbufmod(i)%vector, blocks_u%M*nvec, MPI_DT, receiver, ii+num_blocks, ptree%Comm, request_all(Nreq), ierr)
                     endif
                  endif
               enddo
            else ! diagonal block: solve the diagonal block, send left to off-digonal blocks
               ii=tag-num_blocks
               call g2l(ii, num_blocks, nprow, 1, iproc, myi)
               call g2l(ii, num_blocks, npcol, 1, jproc1, myj1)
               sendbufx(myj1)%vector = sendbufx(myj1)%vector + vin
               brecv(myi) = brecv(myi)-1
               if(brecv(myi)==0 .and. bmod(myi)==0)then
                  blocks_u => h_mat%Local_blocks(myj1, myi)
                  idx_start = blocks_u%headn
                  call Hmat_Usolve(blocks_u, 'N', idx_start, nvec, sendbufx(myj1)%vector, blocks_u%N, ptree, stats)

                  do i1=1,nprow
                     if(receiverlists(i1,myj1)==1)then
                        receiver = blacs_pnum_wp(nprow,npcol, i1-1, jproc1)
                        Nreq = Nreq + 1
                        call MPI_Isend(sendbufx(myj1)%vector, blocks_u%N*nvec, MPI_DT, receiver, ii, ptree%Comm, request_all(Nreq), ierr)
                     endif
                  enddo
               endif
            endif
            deallocate (vin)
         enddo

         if (Nreq > 0) then
            call MPI_waitall(Nreq, request_all, status_all, ierr)
            deallocate (status_all)
            deallocate (request_all)
         endif


         xglo=0
         do j = 1, h_mat%myAcols
            do i = 1, h_mat%myArows
               blocks_u => h_mat%Local_blocks(j, i)
               if(blocks_u%row_group==blocks_u%col_group)then
                  xglo(blocks_u%headn:blocks_u%headn+blocks_u%N-1,:) = sendbufx(j)%vector
               endif
            enddo
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, xglo, N_glo*nvec, MPI_DT, MPI_SUM, ptree%Comm, ierr)
         xloc(1:nloc,:) = xglo(h_mat%idxs:h_mat%idxe,:)
         deallocate(xglo)
         deallocate(bmod)
         deallocate(brecv)
         deallocate(receiverlists)

         do j=1,h_mat%myAcols
            if(allocated(sendbufx(j)%vector))deallocate(sendbufx(j)%vector)
         enddo
         deallocate(sendbufx)

         do i=1,h_mat%myArows
            if(allocated(sendbufmod(i)%vector))deallocate(sendbufmod(i)%vector)
         enddo
         deallocate(sendbufmod)

      else
         write (*, *) 'XxU^-1 with MPI is not yet implemented'
         stop
      endif

      call MPI_barrier(ptree%Comm, ierr)
      return
   end subroutine Hmat_Usolve_Toplevel

end module BPACK_Solve_Mul
