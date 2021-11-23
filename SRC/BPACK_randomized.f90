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

#include "ButterflyPACK_config.fi"
module BPACK_randomMVP
   use BPACK_DEFS
   use MISC_Utilities
   use Bplus_randomizedop
   use BPACK_Solve_Mul
   use Bplus_compress

contains

   subroutine matvec_user(trans, M, N, num_vect, Vin, Vout, ker)

      class(*), pointer :: quant
      integer, INTENT(IN):: M, N, num_vect
      DT::Vin(:, :), Vout(:, :)
      ! type(mesh)::msh
      ! type(proctree)::ptree
      type(kernelquant)::ker
      ! type(Hstat)::stats
      character trans

      procedure(F_HMatVec), POINTER :: proc
      proc => ker%FuncHMatVec
      call proc(trans, M, N, num_vect, Vin, Vout, ker%QuantApp)

      return

   end subroutine matvec_user

!**** Computation of the construction phase with matrix-vector multiplication
   subroutine BPACK_construction_Matvec(bmat, blackbox_BMAT_MVP, Memory, error, option, stats, ker, ptree, msh)
      implicit none

      real(kind=8):: Memory, error, t1, t2
      procedure(HMatVec)::blackbox_BMAT_MVP
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree

      if (allocated(msh%xyz)) deallocate (msh%xyz)

      t1 = OMP_get_wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "FastMATVEC-based Matrix construction......"
      select case (option%format)
      case (HODLR)
         call HODLR_randomized(bmat%ho_bf, blackbox_BMAT_MVP, Memory, error, option, stats, ker, ptree, msh)
      case (HMAT)
         write (*, *) 'FastMATVEC-based Matrix construction is not yet supported for Hmatrix'
         stop
      end select
      t2 = OMP_get_wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "FastMATVEC-based Matrix construction finished", t2 - t1, 'secnds. Error: ', error

   end subroutine BPACK_construction_Matvec

   subroutine HODLR_randomized(ho_bf1, blackbox_HODLR_MVP, Memory, error, option, stats, ker, ptree, msh)


      implicit none
      real(kind=8):: n1, n2, n3, n4, Memory, error_inout, error_lastiter, Memtmp, tmpfact, error, tmp1, tmp2, norm1, norm2
      integer level_c, level_butterfly, bb, rank_new_max, ii, groupm, groupn, Nloc, rank_max_lastiter, rank_max_lastlevel, rank_pre_max, converged
      type(matrixblock), pointer::block_o, block_ref
      DT, allocatable::Vin(:, :), Vout1(:, :), Vout2(:, :)
      type(matrixblock), allocatable::block_rand(:)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      real(kind=8):: time_gemm1, tolerance_abs
      type(kernelquant)::ker
      procedure(HMatVec)::blackbox_HODLR_MVP
      type(proctree)::ptree
      type(mesh)::msh
      integer Bidxs, Bidxe, ierr, tt
      integer vecCNT,num_vect,nn,mm,ranktmp,rank,mn
      DT, allocatable:: RandVectIn(:, :),RandVectOut(:, :)
      DTR, allocatable:: Singular(:)

      if (.not. allocated(stats%rankmax_of_level)) allocate (stats%rankmax_of_level(0:ho_bf1%Maxlevel))
      stats%rankmax_of_level = 0
      if (.not. allocated(stats%rankmax_of_level_global)) allocate (stats%rankmax_of_level_global(0:ho_bf1%Maxlevel))
      stats%rankmax_of_level_global = 0

      rank_max_lastlevel = option%rank0
      Nloc = msh%idxe - msh%idxs + 1

      call assert(option%less_adapt == 0, 'HODLR_randomized does not support less_adapt=1 currently')

      n3 = OMP_get_wtime()


      !!!!!!!!!!!!!!!!!!  get the 2-norm of the HODLR and compute an absolute tolerance
      nn=msh%Nunk
      mm=msh%Nunk
      num_vect = min(10, msh%Nunk)
      allocate (RandVectIn(Nloc, num_vect))
      allocate (RandVectOut(Nloc, num_vect))
      RandVectOut = 0
      call RandomMat(Nloc, num_vect, min(Nloc, num_vect), RandVectIn, 1)
      ! computation of AR
      call blackbox_HODLR_MVP('N', Nloc, Nloc, num_vect, RandVectIn, RandVectOut, ker)

      tmp1 = fnorm(RandVectOut, Nloc, num_vect)**2d0
      call MPI_ALLREDUCE(tmp1, norm1, 1, MPI_double_precision, MPI_SUM, ptree%pgrp(ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      norm1 = sqrt(norm1/num_vect)

      ! computation of range Q of AR
      call PComputeRange(ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%M_p, num_vect, RandVectOut, ranktmp, BPACK_SafeEps, ptree, ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%pgno)
      ! computation of B = Q^c*A
      RandVectOut = conjg(cmplx(RandVectOut, kind=8))
      call blackbox_HODLR_MVP('T', Nloc, Nloc, num_vect, RandVectOut, RandVectIn, ker)
      RandVectOut = conjg(cmplx(RandVectOut, kind=8))
      ! computation of singular values of B
      mn = min(mm, ranktmp)
      allocate (Singular(mn))
      Singular = 0
      call PSVDTruncateSigma(ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1), RandVectIn, ranktmp, rank, Singular, option, stats, ptree)
      tolerance_abs = min(Singular(1)*sqrt(dble(mm)),norm1)*option%tol_Rdetect
      ! if(ptree%MyID==0)write(*,*)Singular(1)*sqrt(dble(mm)),norm1

      deallocate (Singular)
      deallocate (RandVectIn)
      deallocate (RandVectOut)


      Memory = 0
      do level_c = 1, ho_bf1%Maxlevel + 1
         if (level_c == ho_bf1%Maxlevel + 1) then
            call HODLR_randomized_OneL_Fullmat(ho_bf1, blackbox_HODLR_MVP, Nloc, level_c, Memtmp, ker, ptree, option, stats, msh)
            stats%Mem_Direct_for = stats%Mem_Direct_for + Memtmp
         else
            if (level_c > option%LRlevel) then
               level_butterfly = 0
            else
               ! level_butterfly=int((ho_bf1%Maxlevel-level_c)/2)*2
               level_butterfly = ho_bf1%Maxlevel - level_c
            endif

            if (level_c /= ho_bf1%Maxlevel + 1) then
               Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
               Bidxe = ho_bf1%levels(level_c)%Bidxe*2
            else
               Bidxs = ho_bf1%levels(level_c)%Bidxs
               Bidxe = ho_bf1%levels(level_c)%Bidxe
            endif

            converged = 0
            rank_max_lastiter = rank_max_lastlevel
            error_lastiter = BPACK_Bigvalue
            do tt = 1, option%itermax
               rank_pre_max = ceiling_safe(rank_max_lastlevel*option%rankrate**(tt - 1)) + 1

               if (level_butterfly == 0) then
                  n1 = OMP_get_wtime()
                  allocate (block_rand(Bidxe - Bidxs + 1))
                  do bb = Bidxs, Bidxe
                     ! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
                     groupm = ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here
                     groupn = ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here
                     call BF_Init_randomized(level_butterfly, rank_pre_max, groupm, groupn, ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1), block_rand(bb - Bidxs + 1), msh, ptree, option, 1)
                     ! endif
                  enddo
                  n2 = OMP_get_wtime()
                  stats%Time_random(1) = stats%Time_random(1) + n2 - n1

                  call HODLR_randomized_OneL_Lowrank(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, rank_pre_max, option, ker, ptree, stats, msh, tolerance_abs)
               else

                  n1 = OMP_get_wtime()
                  allocate (block_rand(Bidxe - Bidxs + 1))
                  do bb = Bidxs, Bidxe
                     ! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
                     groupm = ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here
                     groupn = ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here
                     call BF_Init_randomized(level_butterfly, rank_pre_max, groupm, groupn, ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1), block_rand(bb - Bidxs + 1), msh, ptree, option, 0)
                     ! endif
                  enddo
                  n2 = OMP_get_wtime()
                  stats%Time_random(1) = stats%Time_random(1) + n2 - n1

                  n1 = OMP_get_wtime()
                  call HODLR_Reconstruction_LL(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, vecCNT, option, stats, ker, ptree, msh)
                  n2 = OMP_get_wtime()
                  if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'reconstructLL: ', n2 - n1, 'vecCNT', vecCNT

                  n1 = OMP_get_wtime()
                  call HODLR_Reconstruction_RR(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, vecCNT, option, stats, ker, ptree, msh)
                  n2 = OMP_get_wtime()
                  if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'reconstructRR: ', n2 - n1, 'vecCNT', vecCNT
               end if

               call HODLR_Test_Error_RR(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, error_inout, ker, ptree, stats, msh, option)

               rank_new_max = 0
               do bb = Bidxs, Bidxe
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                  call BF_get_rank(block_rand(bb - Bidxs + 1), ptree)
                  rank_new_max = max(rank_new_max, block_rand(bb - Bidxs + 1)%rankmax)
               endif
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new_max, 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)

               if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A10,I5,A6,I5,A8,I3, A8,I3,A7,Es14.7,A9,I5)') ' Level ', level_c, ' rank:', rank_new_max, ' Ntrial:', tt, ' L_butt:', level_butterfly, ' error:', error_inout, ' #sample:', rank_pre_max

               ! !!!!*** terminate if 1. error small enough or 2. error not decreasing or 3. rank not increasing
               ! if(error_inout>option%tol_rand .and. error_inout<error_lastiter .and. ((rank_new_max>rank_max_lastiter .and. tt>1).or.tt==1))then

               ! !!!!*** terminate if 1. error small enough or 2. error not decreasing or 3. rank smaller than num_vec
               ! if(error_inout>option%tol_rand .and. error_inout<error_lastiter .and. rank_new_max==rank_pre_max)then

               !!!!*** terminate if 1. error small enough or 2. rank smaller than num_vec
               if (error_inout > option%tol_rand .and. rank_new_max == rank_pre_max) then
                  do bb = Bidxs, Bidxe
                     call BF_delete(block_rand(bb - Bidxs + 1), 1)
                  end do
                  deallocate (block_rand)
                  error_lastiter = error_inout
                  rank_max_lastiter = rank_new_max
               else
                  do bb = Bidxs, Bidxe
                  if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                     block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
                     call BF_copy('N', block_rand(bb - Bidxs + 1), block_o, Memtmp)
                     Memory = Memory + Memtmp
                  endif
                  end do

                  do bb = Bidxs, Bidxe
                     call BF_delete(block_rand(bb - Bidxs + 1), 1)
                  end do
                  deallocate (block_rand)
                  stats%rankmax_of_level(level_c) = rank_new_max
                  rank_max_lastlevel = rank_new_max
                  converged = 1
                  exit
               endif

            end do
            if (converged == 0) then
               write (*, *) 'randomized scheme not converged. level: ', level_c, ' rank:', rank_new_max, ' L_butt:', level_butterfly, ' error:', error_inout
               stop
            end if
         end if
      end do

      n4 = OMP_get_wtime()
      stats%Time_Fill = stats%Time_Fill + n4 - n3

      stats%Mem_Comp_for = stats%Mem_Comp_for + Memory
      call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel), stats%rankmax_of_level_global(0:ho_bf1%Maxlevel), ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level

      allocate (Vin(Nloc, 1))
      allocate (Vout1(Nloc, 1))
      allocate (Vout2(Nloc, 1))
      do ii = 1, Nloc
         call random_dp_number(Vin(ii, 1))
      end do

      call blackbox_HODLR_MVP('N', Nloc, Nloc, 1, Vin, Vout1, ker)
      call HODLR_Mult('N', Nloc, 1, 1, ho_bf1%Maxlevel + 1, Vin, Vout2, ho_bf1, ptree, option, stats)

      tmp1 = fnorm(Vout2 - Vout1, Nloc, 1)**2d0
      call MPI_ALLREDUCE(tmp1, norm1, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      tmp2 = fnorm(Vout1, Nloc, 1)**2d0
      call MPI_ALLREDUCE(tmp2, norm2, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      error = sqrt(norm1)/sqrt(norm2)

      deallocate (Vin, Vout1, Vout2)

   end subroutine HODLR_randomized

   subroutine HODLR_randomized_OneL_Lowrank(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, rmax, option, ker, ptree, stats, msh, tolerance_abs)


      implicit none
      real(kind=8):: n1, n2, Memory, error_inout, Memtmp, tolerance_abs
      integer mn, rankref, level_c, rmax, rmaxloc, level_butterfly, bb, bb1, bb_inv, rank_new_max, rank, num_vect, groupn, groupm, header_n, header_m, tailer_m, tailer_n, ii, jj, k, mm, nn
      type(matrixblock), pointer::block_o, block_ref, block_inv
      DT, allocatable::RandVectTmp(:, :)
      DT, allocatable :: matRcol(:, :), matZRcol(:, :), matRrow(:, :), matZcRrow(:, :), mattmp1(:, :), mattmp2(:, :), matrixtemp(:, :), matrixtemp1(:, :), matrixtempQ(:, :), matrixtempin(:, :), matrixtempout(:, :)
      DTR, allocatable:: Singular(:)
      DT::ctemp1, ctemp2
      DT, allocatable::UU(:, :), VV(:, :)
      integer q, qq, Nloc, pp
      integer, allocatable::perms(:), ranks(:)
      type(hobf)::ho_bf1
      type(matrixblock)::block_rand(:)
      type(Hoption)::option
      DT, allocatable:: RandVectInR(:, :), RandVectOutR(:, :), RandVectInL(:, :), RandVectOutL(:, :)
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(mesh)::msh
      procedure(HMatVec)::blackbox_HODLR_MVP
      integer Bidxs, Bidxe, head, tail, idx_start_loc, idx_end_loc, ierr

      Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
      Bidxe = ho_bf1%levels(level_c)%Bidxe*2

      level_butterfly = 0

      ctemp1 = 1.0d0; ctemp2 = 0.0d0

      rank_new_max = 0

      num_vect = rmax
      allocate (RandVectInR(Nloc, num_vect))
      RandVectInR = 0
      allocate (RandVectOutR(Nloc, num_vect))

      allocate (ranks(Bidxe - Bidxs + 1))
      ranks = rmax

      ! write(*,*)Bidxs,Bidxe,ptree%MyID,'wocao'

      do bb = Bidxs, Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
            block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

            mm = block_o%M_loc

            pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
            header_m = block_o%M_p(pp, 1) + block_o%headm - 1

            k = header_m - msh%idxs
            ranks(bb - Bidxs + 1) = min(min(block_o%M, block_o%N), num_vect)

            allocate (matrixtemp(mm, num_vect))
            call RandomMat(mm, num_vect, min(mm, num_vect), matrixtemp, 1)

            ! write(*,*)1+k,mm+k,header_m,msh%idxs,msh%idxe,block_o%row_group,block_o%col_group,ptree%MyID

            RandVectInR(1 + k:mm + k, 1:num_vect) = matrixtemp
            deallocate (matrixtemp)
         endif
      end do

      call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'N', RandVectInR, RandVectOutR, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)
      ! computation of range Q
      do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
         head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
         tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
         idx_start_loc = head - msh%idxs + 1
         idx_end_loc = tail - msh%idxs + 1
         call PComputeRange_twoforward(ho_bf1, level_c, Bidxs, bb_inv, ranks, RandVectOutR(idx_start_loc:idx_end_loc, 1:num_vect), BPACK_SafeEps, ptree, stats)
      end do


      ! power iteration of order q, orthogonalize each multiplication results , see algorithm 4.4 Halko 2010
      do qq = 1, option%powiter
         RandVectOutR = conjg(cmplx(RandVectOutR, kind=8))
         call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'T', RandVectOutR, RandVectInR, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)
         RandVectInR = conjg(cmplx(RandVectInR, kind=8))
         ! computation of range Q
         do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
            pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
            head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
            tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
            idx_start_loc = head - msh%idxs + 1
            idx_end_loc = tail - msh%idxs + 1
            call PComputeRange_twoforward(ho_bf1, level_c, Bidxs, bb_inv, ranks, RandVectInR(idx_start_loc:idx_end_loc, 1:num_vect), BPACK_SafeEps, ptree, stats)
         end do


         call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'N', RandVectInR, RandVectOutR, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)
         ! computation of range Q
         do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
            pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
            head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
            tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
            idx_start_loc = head - msh%idxs + 1
            idx_end_loc = tail - msh%idxs + 1
            call PComputeRange_twoforward(ho_bf1, level_c, Bidxs, bb_inv, ranks, RandVectOutR(idx_start_loc:idx_end_loc, 1:num_vect), BPACK_SafeEps, ptree, stats)
         end do

      enddo



      !!!!!!!!!!!!!!!!!!!! need add a redistrubtion here

      ! computation of B = Q^c*A
      RandVectOutR = conjg(cmplx(RandVectOutR, kind=8))
      call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'T', RandVectOutR, RandVectInR, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)
      RandVectOutR = conjg(cmplx(RandVectOutR, kind=8))

      ! computation of SVD of B and LR of A

      do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe

         pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
         head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
         tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
         idx_start_loc = head - msh%idxs + 1
         idx_end_loc = tail - msh%idxs + 1

         call PQxSVDTruncate_twoforward(ho_bf1, level_c, Bidxs, bb_inv, ranks, RandVectOutR(idx_start_loc:idx_end_loc, 1:num_vect), RandVectInR(idx_start_loc:idx_end_loc, 1:num_vect), block_rand, option, ptree, stats,tolerance_abs)

      end do

      deallocate (RandVectOutR, RandVectInR, ranks)

   end subroutine HODLR_randomized_OneL_Lowrank

   subroutine HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, trans, VectIn, VectOut, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)


      implicit none
      real(kind=8):: n1, n2, Memory, error_inout, Memtmp
      integer Nloc, mn, rankref, level_c, rmax, rmaxloc, bb, rank_new_max, rank, num_vect, groupn, groupm, header_n, header_m, tailer_m, tailer_n, ii, jj, k, mm, nn
      type(matrixblock), pointer::block_o, block_ref
      DT, allocatable::RandVectTmp(:, :)
      DT, allocatable :: matRcol(:, :), matZRcol(:, :), matRrow(:, :), matZcRrow(:, :), mattmp1(:, :), mattmp2(:, :), matrixtemp(:, :)
      DTR, allocatable:: Singular(:)
      DT::ctemp1, ctemp2
      complex(kind=8), allocatable::UU(:, :), VV(:, :)
      integer, allocatable::perms(:)
      character trans
      DT::VectIn(:, :), VectOut(:, :)
      DT, allocatable:: RandVectIn(:, :), RandVectOut(:, :)
      type(hobf)::ho_bf1
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh
      integer Bidxs, Bidxe, pp

      procedure(HMatVec)::blackbox_HODLR_MVP

      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      Memory = 0
      rank_new_max = 0

      if (level_c /= ho_bf1%Maxlevel + 1) then
         Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
         Bidxe = ho_bf1%levels(level_c)%Bidxe*2

         if (trans == 'N') then

            VectOut = 0

            allocate (RandVectIn(Nloc, num_vect))
            allocate (RandVectOut(Nloc, num_vect))
            allocate (RandVectTmp(Nloc, num_vect))

            ! Compute the odd block MVP first
            RandVectIn = 0
            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 1) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb + 1)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb + 1)%LL(1)%matrices_block(1)
                  mm = block_o%M_loc

                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  RandVectIn(1 + k:mm + k, 1:num_vect) = VectIn(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            call blackbox_HODLR_MVP('N', Nloc, Nloc, num_vect, RandVectIn, RandVectTmp, ker)
            call HODLR_Mult('N', Nloc, num_vect, 1, level_c - 1, RandVectIn, RandVectOut, ho_bf1, ptree, option, stats)
            RandVectOut = RandVectTmp - RandVectOut
            stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 1) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
                  groupm = block_o%row_group  ! Note: row_group and col_group interchanged here

                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  VectOut(1 + k:mm + k, 1:num_vect) = RandVectOut(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            ! Compute the even block MVP next
            RandVectIn = 0
            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 0) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb - 1)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb - 1)%LL(1)%matrices_block(1)
                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  RandVectIn(1 + k:mm + k, 1:num_vect) = VectIn(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            call blackbox_HODLR_MVP('N', Nloc, Nloc, num_vect, RandVectIn, RandVectTmp, ker)

            call HODLR_Mult('N', Nloc, num_vect, 1, level_c - 1, RandVectIn, RandVectOut, ho_bf1, ptree, option, stats)
            RandVectOut = RandVectTmp - RandVectOut
            stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 0) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
                  groupm = block_o%row_group  ! Note: row_group and col_group interchanged here

                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  VectOut(1 + k:mm + k, 1:num_vect) = RandVectOut(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            deallocate (RandVectIn)
            deallocate (RandVectOut)
            deallocate (RandVectTmp)

         else if (trans == 'T') then
            VectOut = 0

            allocate (RandVectIn(Nloc, num_vect))
            allocate (RandVectOut(Nloc, num_vect))

            allocate (RandVectTmp(Nloc, num_vect))

            ! Compute the odd block MVP first
            RandVectIn = 0

            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 1) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
                  groupm = block_o%row_group  ! Note: row_group and col_group interchanged here
                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  RandVectIn(1 + k:mm + k, 1:num_vect) = VectIn(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            call blackbox_HODLR_MVP('T', Nloc, Nloc, num_vect, RandVectIn, RandVectTmp, ker)

            call HODLR_Mult('T', Nloc, num_vect, 1, level_c - 1, RandVectIn, RandVectOut, ho_bf1, ptree, option, stats)
            RandVectOut = RandVectTmp - RandVectOut
            stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 1) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb + 1)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb + 1)%LL(1)%matrices_block(1)
                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  VectOut(1 + k:mm + k, 1:num_vect) = RandVectOut(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            ! Compute the even block MVP next
            RandVectIn = 0
            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 0) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
                  groupm = block_o%row_group  ! Note: row_group and col_group interchanged here
                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  RandVectIn(1 + k:mm + k, 1:num_vect) = VectIn(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            call blackbox_HODLR_MVP('T', Nloc, Nloc, num_vect, RandVectIn, RandVectTmp, ker)
            call HODLR_Mult('T', Nloc, num_vect, 1, level_c - 1, RandVectIn, RandVectOut, ho_bf1, ptree, option, stats)
            stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
            RandVectOut = RandVectTmp - RandVectOut

            do bb = Bidxs, Bidxe
               if (mod(bb, 2) == 0) then
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb - 1)%pgno)) then
                  block_o => ho_bf1%levels(level_c)%BP(bb - 1)%LL(1)%matrices_block(1)
                  mm = block_o%M_loc
                  pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
                  header_m = block_o%M_p(pp, 1) + block_o%headm - 1
                  k = header_m - msh%idxs
                  VectOut(1 + k:mm + k, 1:num_vect) = RandVectOut(1 + k:mm + k, 1:num_vect)
               endif
               endif
            end do

            deallocate (RandVectIn)
            deallocate (RandVectOut)
            deallocate (RandVectTmp)

         endif

      else

         VectOut = 0
         allocate (RandVectTmp(Nloc, num_vect))
         ! Compute the odd block MVP first
         call blackbox_HODLR_MVP('N', Nloc, Nloc, num_vect, VectIn, RandVectTmp, ker)
         call HODLR_Mult('N', Nloc, num_vect, 1, level_c - 1, VectIn, VectOut, ho_bf1, ptree, option, stats)
         VectOut = RandVectTmp - VectOut
         stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
         deallocate (RandVectTmp)
      endif

   end subroutine HODLR_MVP_randomized_OneL

   subroutine HODLR_randomized_OneL_Fullmat(ho_bf1, blackbox_HODLR_MVP, N, level_c, Memory, ker, ptree, option, stats, msh)


      implicit none
      real(kind=8):: n1, n2, Memory, error_inout, Memtmp
      integer N, rankref, level_c, rmaxloc, level_butterfly, bb, rank_new_max, rank, num_vect, groupn, groupm, header_n, header_m, tailer_m, tailer_n, ii, jj, k, mm, nn
      type(matrixblock), pointer::block_o, block_ref
      DT, allocatable::RandVectTmp(:, :)
      DT, allocatable :: matRcol(:, :), matZRcol(:, :), matRrow(:, :), matZcRrow(:, :), mattmp1(:, :), mattmp2(:, :)
      DTR, allocatable:: Singular(:)
      DT::ctemp1, ctemp2
      type(hobf)::ho_bf1
      DT, allocatable:: RandVectInR(:, :), RandVectOutR(:, :), RandVectInL(:, :), RandVectOutL(:, :)
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh
      integer ierr, tempi
      integer Bidxs, Bidxe, N_unk_loc

      procedure(HMatVec)::blackbox_HODLR_MVP

      if (level_c /= ho_bf1%Maxlevel + 1) then
         Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
         Bidxe = ho_bf1%levels(level_c)%Bidxe*2
      else
         Bidxs = ho_bf1%levels(level_c)%Bidxs
         Bidxe = ho_bf1%levels(level_c)%Bidxe
      endif

      ctemp1 = 1.0d0; ctemp2 = 0.0d0

      Memory = 0
      rank_new_max = 0

      num_vect = 0

      do bb = Bidxs, Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
            block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
            nn = block_o%N
            num_vect = max(num_vect, nn)
         endif
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, num_vect, 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)

      N_unk_loc = msh%idxe - msh%idxs + 1

      allocate (RandVectInR(N_unk_loc, num_vect))
      RandVectInR = 0
      allocate (RandVectTmp(N_unk_loc, num_vect))
      allocate (RandVectOutR(N_unk_loc, num_vect))

      do bb = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
         mm = block_o%M_loc
         header_m = block_o%headm
         k = header_m - msh%idxs
         do ii = 1, mm
            RandVectInR(ii + k, ii) = 1d0
         enddo
      end do

      call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'N', RandVectInR, RandVectOutR, N_unk_loc, level_c, num_vect, ker, ptree, stats, msh, option)

      do bb = Bidxs, Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
            block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
            groupm = block_o%row_group  ! Note: row_group and col_group interchanged here
            mm = block_o%M_loc
            header_m = block_o%headm
            k = header_m - msh%idxs

            nn = block_o%N_loc

            block_o%style = 1
            allocate (block_o%fullmat(mm, nn))
            ! call copymatN(RandVectOutR(k+1:k+mm,1:nn),block_o%fullmat,mm,nn)
            block_o%fullmat = RandVectOutR(k + 1:k + mm, 1:nn)
            Memory = Memory + SIZEOF(block_o%fullmat)/1024.0d3
         endif
      end do

      deallocate (RandVectInR, RandVectOutR, RandVectTmp)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A10,I5,A13)') ' Level ', level_c, ' fullmat done'

   end subroutine HODLR_randomized_OneL_Fullmat

   subroutine HODLR_Reconstruction_LL(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, vecCNT, option, stats, ker, ptree, msh)


      implicit none

      integer level_c, rowblock, Nloc
      integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
      integer i, j, ii, jj, bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb, Bidxs
      integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, num_vect_sub, num_vect_subsub
      DT ctemp, a, b
      character chara
      integer level_right_start

      type(RandomBlock), allocatable :: vec_rand(:)
      integer Nsub, Ng, nth, nth_s, nth_e
      integer Nbind
      real(kind=8)::n1, n2

      integer blocks1, blocks2, blocks3, level_butterfly
      integer tt
      type(matrixblock), pointer::block_o

      integer::rank_new_max, dimension_rank
      real(kind=8)::rank_new_avr, error
      DT, allocatable::matrixtmp(:, :)
      integer niter, unique_nth  ! level# where each block is touched only once
      real(kind=8):: error_inout
      integer, allocatable::perms(:)

      type(matrixblock)::block_rand(:)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(kernelquant)::ker
      type(proctree)::ptree
      type(mesh)::msh
      DT, allocatable:: RandVectIn(:, :), RandVectOut(:, :)
      integer bb_inv, idx_start_loc, idx_end_loc, pp, head, tail
      integer vecCNT

      procedure(HMatVec)::blackbox_HODLR_MVP

      vecCNT = 0

      dimension_rank = block_rand(1)%dimension_rank   ! be careful here
      num_blocks = 2**level_butterfly
      num_vect_subsub = dimension_rank + vec_oversample ! be careful with the oversampling factor here
      level_right_start = block_rand(1)%level_half
      Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1

      do level = 0, level_right_start
         Nsub = NINT(2**ceiling_safe((level_butterfly - 1)/2d0)/dble(2**(level_right_start - level)))   !  check here later
         Ng = 2**level_butterfly/Nsub

         Nbind = min(option%Nbundle, Nsub)
         num_vect_sub = num_vect_subsub*Nbind

         allocate (RandVectIn(Nloc, num_vect_sub))
         RandVectIn = 0
         allocate (RandVectOut(Nloc, num_vect_sub))
         RandVectOut = 0

         do ii = 1, Nsub/Nbind
            nth_s = (ii - 1)*Nbind + 1
            nth_e = ii*Nbind

            n1 = OMP_get_wtime()
            call HODLR_Randomized_Vectors('L', ho_bf1, block_rand, RandVectIn, RandVectOut, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, nth_s, nth_e, num_vect_sub, level, ker, ptree, stats, msh, option)
            vecCNT = vecCNT + num_vect_sub*2

            n2 = OMP_get_wtime()
            ! time_getvec = time_getvec + n2-n1
            stats%Time_Random(2) = stats%Time_Random(2) + n2 - n1

            n1 = OMP_get_wtime()

            do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
               pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
               head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
               tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
               idx_start_loc = head - msh%idxs + 1
               idx_end_loc = tail - msh%idxs + 1

               call BF_Resolving_Butterfly_LL_dat_twoforward(ho_bf1, level_c, num_vect_sub, nth_s, nth_e, Ng, level, Bidxs, bb_inv, block_rand, RandVectIn(idx_start_loc:idx_end_loc, 1:num_vect_sub), RandVectOut(idx_start_loc:idx_end_loc, 1:num_vect_sub), option, ptree, msh, stats)

            enddo
            n2 = OMP_get_wtime()
            stats%Time_Random(3) = stats%Time_Random(3) + n2 - n1
         end do

         deallocate (RandVectIn)
         deallocate (RandVectOut)

      end do

      return

   end subroutine HODLR_Reconstruction_LL

   subroutine HODLR_Reconstruction_RR(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, vecCNT, option, stats, ker, ptree, msh)


      implicit none

      integer level_c, rowblock, Nloc
      integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
      integer i, j, ii, jj, bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb, Bidxs
      integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, num_vect_sub, num_vect_subsub
      DT ctemp, a, b
      character chara
      integer level_left_start

      type(RandomBlock), allocatable :: vec_rand(:)
      integer Nsub, Ng, nth, nth_s, nth_e
      integer Nbind
      real(kind=8)::n1, n2

      integer blocks1, blocks2, blocks3, level_butterfly
      integer tt
      type(matrixblock), pointer::block_o

      integer::rank_new_max, dimension_rank
      real(kind=8)::rank_new_avr, error
      DT, allocatable::matrixtmp(:, :)
      integer niter, unique_nth  ! level# where each block is touched only once
      real(kind=8):: error_inout
      integer, allocatable::perms(:)

      type(matrixblock)::block_rand(:)
      type(hobf)::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(kernelquant)::ker
      type(proctree)::ptree
      type(mesh)::msh
      DT, allocatable:: RandVectIn(:, :), RandVectOut(:, :)
      integer bb_inv, idx_start_loc, idx_end_loc, pp, head, tail
      integer vecCNT
      procedure(HMatVec)::blackbox_HODLR_MVP

      vecCNT = 0

      dimension_rank = block_rand(1)%dimension_rank   ! be careful here
      num_blocks = 2**level_butterfly
      num_vect_subsub = dimension_rank + vec_oversample ! be careful with the oversampling factor here
      level_left_start = block_rand(1)%level_half + 1 !  check here later
      Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1

      do level = level_butterfly + 1, level_left_start, -1
         Nsub = NINT(2**ceiling_safe((level_butterfly)/2d0)/dble(2**(level - level_left_start)))    !  check here later
         Ng = 2**level_butterfly/Nsub

         Nbind = min(option%Nbundle, Nsub)
         num_vect_sub = num_vect_subsub*Nbind

         allocate (RandVectIn(Nloc, num_vect_sub))
         RandVectIn = 0
         allocate (RandVectOut(Nloc, num_vect_sub))
         RandVectOut = 0

         do ii = 1, Nsub/Nbind
            nth_s = (ii - 1)*Nbind + 1
            nth_e = ii*Nbind

            n1 = OMP_get_wtime()
            call HODLR_Randomized_Vectors('R', ho_bf1, block_rand, RandVectIn, RandVectOut, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, nth_s, nth_e, num_vect_sub, level, ker, ptree, stats, msh, option)
            vecCNT = vecCNT + num_vect_sub*2

            n2 = OMP_get_wtime()
            ! time_getvec = time_getvec + n2-n1
            stats%Time_Random(2) = stats%Time_Random(2) + n2 - n1

            n1 = OMP_get_wtime()

            do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
               pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
               head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
               tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
               idx_start_loc = head - msh%idxs + 1
               idx_end_loc = tail - msh%idxs + 1

               call BF_Resolving_Butterfly_RR_dat_twoforward(ho_bf1, level_c, num_vect_sub, nth_s, nth_e, Ng, level, Bidxs, bb_inv, block_rand, RandVectIn(idx_start_loc:idx_end_loc, 1:num_vect_sub), RandVectOut(idx_start_loc:idx_end_loc, 1:num_vect_sub), option, ptree, msh, stats)

            enddo
            n2 = OMP_get_wtime()
            stats%Time_Random(3) = stats%Time_Random(3) + n2 - n1
         end do

         deallocate (RandVectIn)
         deallocate (RandVectOut)

      end do

      return

   end subroutine HODLR_Reconstruction_RR

   subroutine HODLR_Test_Error_RR(ho_bf1, block_rand, blackbox_HODLR_MVP, Nloc, level_c, error, ker, ptree, stats, msh, option)


      implicit none

      integer nth
      integer i, j, k, level, ii, jj, kk, test, groupm, groupn, bb, pp
      integer Nloc, mm, nn
      real(kind=8) a, b, c, d, condition_number, norm1_R, norm2_R, norm3_R, norm4_R
      DT ctemp, ctemp1, ctemp2

      ! type(matricesblock), pointer :: blocks
      type(RandomBlock), pointer :: random
      integer Nsub, Ng, num_vect, nth_s, nth_e, level_butterfly, ierr
      integer*8 idx_start
      real(kind=8)::error, tmp1, tmp2, norm1, norm2
      integer level_c, rowblock, dimension_m, header_m, tailer_m, header_n, tailer_n
      DT, allocatable::RandomVectors_Output_ref(:, :)
      type(matrixblock), pointer::block_o

      type(hobf)::ho_bf1
      type(matrixblock)::block_rand(:)
      type(vectorsblock), pointer:: RandomVectors_InOutput(:)
      type(kernelquant)::ker
      type(proctree)::ptree
      procedure(HMatVec)::blackbox_HODLR_MVP
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh
      integer Bidxs, Bidxe, bb_inv, head, tail, idx_start_loc, idx_end_loc

      if (level_c /= ho_bf1%Maxlevel + 1) then
         Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
         Bidxe = ho_bf1%levels(level_c)%Bidxe*2
      else
         Bidxs = ho_bf1%levels(level_c)%Bidxs
         Bidxe = ho_bf1%levels(level_c)%Bidxe
      endif

      num_vect = 1
      allocate (RandomVectors_InOutput(3))
      do ii = 1, 3
         allocate (RandomVectors_InOutput(ii)%vector(Nloc, num_vect))
         RandomVectors_InOutput(ii)%vector = 0d0
      end do

      do bb = Bidxs, Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
            block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
            mm = block_o%M_loc
            pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
            header_m = block_o%M_p(pp, 1) + block_o%headm - 1
            k = header_m - msh%idxs
            do jj = 1, num_vect
               do ii = 1, mm
                  call random_dp_number(RandomVectors_InOutput(1)%vector(ii + k, jj))        ! matrixtemp1(jj,ii) !
               enddo
            enddo

            ! !!! this is needed to use HODLR_Mult, otherwise I need to write a separate subroutine for matvec of block_rand at one level
            ! write(*,*)block_rand(bb-Bidxs+1)%row_group,block_rand(bb-Bidxs+1)%col_group,'b copy',bb-Bidxs+1,ptree%MyID
            ! call BF_copy('N',block_rand(bb-Bidxs+1),block_o)
            ! write(*,*)block_o%row_group,block_o%col_group,'a copy'
         endif
      end do

      call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, 'N', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, Nloc, level_c, num_vect, ker, ptree, stats, msh, option)

      do bb_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
         head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp, 1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn - 1
         tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc - 1
         idx_start_loc = head - msh%idxs + 1
         idx_end_loc = tail - msh%idxs + 1
         if (level_c == ho_bf1%Maxlevel + 1) then
            call Full_block_MVP_dat(block_rand(bb_inv - Bidxs + 1), 'N', idx_end_loc - idx_start_loc + 1, num_vect,&
                        &RandomVectors_InOutput(1)%vector(idx_start_loc, 1),Nloc, RandomVectors_InOutput(2)%vector(idx_start_loc, 1),Nloc, BPACK_cone, BPACK_czero)
         else
            call BF_block_MVP_twoforward_dat(ho_bf1, level_c, bb_inv, block_rand, 'N', idx_end_loc - idx_start_loc + 1, num_vect, RandomVectors_InOutput(1)%vector(idx_start_loc, 1), Nloc, RandomVectors_InOutput(2)%vector(idx_start_loc, 1), Nloc, BPACK_cone, BPACK_czero, ptree, stats)
         endif
      end do

      ! call HODLR_Mult('N',Nloc,num_vect,level_c,level_c,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ho_bf1,ptree,option,stats)
      ! do bb =Bidxs,Bidxe
      ! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
      ! block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
      ! !!! this is needed to use HODLR_Mult, otherwise I need to write a seperate subroutine for matvec of block_rand at one level
      ! ! write(*,*)block_o%row_group,block_o%col_group,'b delete'
      ! call BF_delete(block_o,1)
      ! ! write(*,*)block_o%row_group,block_o%col_group,'a delete'
      ! endif
      ! end do

      tmp1 = fnorm(RandomVectors_InOutput(2)%vector - RandomVectors_InOutput(3)%vector, Nloc, 1)**2d0
      call MPI_ALLREDUCE(tmp1, norm1, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      tmp2 = fnorm(RandomVectors_InOutput(3)%vector, Nloc, 1)**2d0
      call MPI_ALLREDUCE(tmp2, norm2, 1, MPI_double_precision, MPI_SUM, ptree%Comm, ierr)
      error = sqrt(norm1)/sqrt(norm2)
      ! error=0
      do ii = 1, 3
         deallocate (RandomVectors_InOutput(ii)%vector)
      end do
      deallocate (RandomVectors_InOutput)

      return

   end subroutine HODLR_Test_Error_RR

   subroutine HODLR_Randomized_Vectors(side, ho_bf1, block_rand, RandVectIn, RandVectOut, blackbox_HODLR_MVP, Nloc, level_c, level_butterfly, nth_s, nth_e, num_vect_sub, unique_nth, ker, ptree, stats, msh, option)




      implicit none

      integer level_c, rowblock, unique_nth
      integer Nloc, i, j, k, level, ii, jj, kk, test, bb
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character chara
      real(kind=8) a, b, c, d
      ! DT ctemp, ctemp1, ctemp2
      type(matrixblock), pointer::block_o

      type(vectorsblock), pointer :: random1, random2

      DTR, allocatable :: Singular(:)
      integer idx_start_glo, N_diag, idx_start_diag, idx_start_loc, idx_end_loc
      DT, allocatable::vec_old(:, :), vec_new(:, :), matrixtemp1(:, :)

      integer Nsub, Ng
      integer*8 idx_start
      integer level_blocks
      integer groupm_start, groupn_start
      integer header_mm, header_nn
      integer header_m, header_n, tailer_m, tailer_n

      integer nth_s, nth_e, num_vect_sub, nth, num_vect_subsub
      real(kind=8)::n1, n2

      type(hobf)::ho_bf1
      type(matrixblock)::block_rand(:)
      DT:: RandVectIn(:, :), RandVectOut(:, :)
      type(vectorsblock), pointer:: RandomVectors_InOutput(:)
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh
      integer Bidxs, Bidxe
      integer idx_r, inc_r, nr, idx_c, inc_c, nc, idx_m_s, idx_m_e, idxs, idxe
      character::trans, side
      integer level_left_start, level_right_start

      procedure(HMatVec)::blackbox_HODLR_MVP

      RandVectIn = 0
      RandVectOut = 0
      num_vect_subsub = num_vect_sub/(nth_e - nth_s + 1)

      if (side == 'L') then
         trans = 'T'
         level_right_start = block_rand(1)%level_half
         Nsub = NINT(2**ceiling_safe((level_butterfly - 1)/2d0)/dble(2**(level_right_start - unique_nth)))   !  check here later
      else if (side == 'R') then
         trans = 'N'
         level_left_start = block_rand(1)%level_half + 1 !  check here later
         Nsub = NINT(2**ceiling_safe((level_butterfly)/2d0)/dble(2**(unique_nth - level_left_start)))    !  check here later
      endif

      Ng = 2**level_butterfly/Nsub
      Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
      Bidxe = ho_bf1%levels(level_c)%Bidxe*2

      do bb = Bidxs, Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(bb)%pgno)) then
            block_o => ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
            ! mm=block_o%M_loc
            ! pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
            ! header_m = block_o%M_p(pp,1) + block_o%headm -1
            ! k=header_m-msh%idxs

            groupm = block_o%row_group  ! Note: row_group and col_group interchanged here
            groupm_start = groupm*2**(level_butterfly)

            call GetLocalBlockRange(ptree, block_rand(bb - Bidxs + 1)%pgno, level_butterfly + 1, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
            idx_m_s = idx_r
            idx_m_e = idx_r + nr - 1

            do nth = nth_s, nth_e
               !$omp parallel do default(shared) private(i,mm,idxs,idxe)
               do i = (nth - 1)*Ng + 1, nth*Ng
                  if (i >= idx_m_s .and. i <= idx_m_e) then
                     idxs = msh%basis_group(groupm_start + i - 1)%head - msh%idxs + 1
                     idxe = msh%basis_group(groupm_start + i - 1)%tail - msh%idxs + 1
                     mm = idxe - idxs + 1
                     call RandomSubMat(idxs, idxe, (nth - nth_s)*num_vect_subsub + 1, (nth - nth_s)*num_vect_subsub + num_vect_subsub, min(mm, num_vect_subsub), RandVectIn, 0)
                  endif
               enddo
               !$omp end parallel do
            enddo
         endif
      enddo

      call HODLR_MVP_randomized_OneL(ho_bf1, blackbox_HODLR_MVP, trans, RandVectIn, RandVectOut, Nloc, level_c, num_vect_sub, ker, ptree, stats, msh, option)

      return

   end subroutine HODLR_Randomized_Vectors

   subroutine PComputeRange_twoforward(ho_bf1, level, Bidxs, ii, ranks, AR, eps, ptree, stats)

      implicit none
      integer ranks(:)
      integer level, ii, bb
      DT :: AR(:, :)
      DT, pointer :: matrixtemp1(:, :), matrixtemp2(:, :), matrixtemp(:, :)
      type(matrixblock), pointer::block_o, block_inv, block_schur, block_off1, block_off2, block_off
      integer tempi, groupn, groupm, mm(2), mm1, nn1, mm2, nn2, ierr, nin1, nout1, nin2, nout2, offout(2), offout1, offout2, rank
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      real(kind=8)::n1, n2, eps, flop
      integer, pointer::M_p(:, :)
      integer Bidxs

      block_off1 => ho_bf1%levels(level)%BP(ii*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level)%BP(ii*2)%LL(1)%matrices_block(1)
      block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)

      mm(1) = block_off1%M_loc
      mm(2) = block_off2%M_loc

      offout(1) = 0
      offout(2) = block_off1%M

      !!!*** redistribute AR from process layout of hodlr to the process layout of block_off1 and block_off2
      call MPI_ALLREDUCE(MPI_IN_PLACE, ranks(ii*2 - 1 - Bidxs + 1), 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(block_inv%pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ranks(ii*2 - Bidxs + 1), 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(block_inv%pgno)%Comm, ierr)

      do bb = 1, 2
         block_off => ho_bf1%levels(level)%BP(ii*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         ! if (mm(bb) > 0) then
            if (bb == 1) then
               allocate (matrixtemp1(max(1,mm(bb)), ranks(ii*2 - 1 + bb - 1 - Bidxs + 1)))
               matrixtemp => matrixtemp1
            endif
            if (bb == 2) then
               allocate (matrixtemp2(max(1,mm(bb)), ranks(ii*2 - 1 + bb - 1 - Bidxs + 1)))
               matrixtemp => matrixtemp2
            endif
         ! endif
         n1 = OMP_get_wtime()
         call Redistribute1Dto1D(AR, size(AR,1), block_inv%M_p, 0, block_inv%pgno, matrixtemp, max(1,mm(bb)), M_p, offout(bb), block_off%pgno, ranks(ii*2 - 1 + bb - 1 - Bidxs + 1), ptree)
         n2 = OMP_get_wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2 - n1
      enddo

      !!!*** compute range of AR from QR for block_off1 and block_off2
      do bb = 1, 2
         block_off => ho_bf1%levels(level)%BP(ii*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         if (bb == 1) matrixtemp => matrixtemp1
         if (bb == 2) matrixtemp => matrixtemp2

         if (mm(bb) > 0) then
            call PComputeRange(M_p, ranks(ii*2 - 1 + bb - 1 - Bidxs + 1), matrixtemp, rank, eps, ptree, block_off%pgno, flop)

            ranks(ii*2 - 1 + bb - 1 - Bidxs + 1) = rank
            stats%Flop_Fill = stats%Flop_Fill + flop
         endif
      enddo

      call MPI_ALLREDUCE(MPI_IN_PLACE, ranks(ii*2 - 1 - Bidxs + 1), 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(block_inv%pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ranks(ii*2 - Bidxs + 1), 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(block_inv%pgno)%Comm, ierr)
      ! write(*,*)'wonima',ii*2-1-Bidxs+1,ranks(ii*2-1-Bidxs+1),ii*2-Bidxs+1,ranks(ii*2-Bidxs+1),ptree%MyID,mm(1)

      !!!*** redistribute AR from process layout of block_off1 and block_off2  to the process layout of hodlr
      do bb = 1, 2
         block_off => ho_bf1%levels(level)%BP(ii*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         if (bb == 1) matrixtemp => matrixtemp1
         if (bb == 2) matrixtemp => matrixtemp2

         n1 = OMP_get_wtime()
         call Redistribute1Dto1D(matrixtemp, mm(bb), M_p, offout(bb), block_off%pgno, AR, size(AR,1),block_inv%M_p, 0, block_inv%pgno, ranks(ii*2 - 1 + bb - 1 - Bidxs + 1), ptree)
         n2 = OMP_get_wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2 - n1
         deallocate (matrixtemp)
      enddo

   end subroutine PComputeRange_twoforward

!!!!!***** this subroutine is part of the randomized SVD.
! Given B^T = (Q^cA)^T (N_loc x ranks(bb)) and Q (M_loc x ranks(bb)) in the process layout of hodlr, it computes SVD B=USV and output A = (QU)*(SV)
   subroutine PQxSVDTruncate_twoforward(ho_bf1, level, Bidxs, bb_inv, ranks, Q, QcA_trans, block_rand, option, ptree, stats,tolerance_abs)

      implicit none
      integer ranks(:)
      integer level, ii, bb, bb_inv
      DT :: Q(:, :), QcA_trans(:, :)
      DT, pointer :: mat1(:, :), mat2(:, :), mat(:, :), matQ1(:, :), matQ2(:, :), matQ(:, :), matQ2D(:, :), matQcA_trans1(:, :), matQcA_trans2(:, :), matQcA_trans(:, :), matQcA_trans2D(:, :), matQUt2D(:, :), UU(:, :), VV(:, :), mattemp(:, :)
      type(matrixblock), pointer::block_o, block_inv, block_schur, block_off1, block_off2, block_off
      integer groupn, groupm, mm(2), nn(2), ierr, nin1, nout1, nin2, nout2, offM(2), offN(2), offout1, offout2, rank
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::n1, n2, eps, flop, tolerance_abs
      integer, pointer::M_p(:, :), N_p(:, :)
      DTR, pointer::Singular(:)
      integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9), descQUt2D(9)
      integer tempi, ctxt, info, iproc, jproc, myi, myj, myArows, myAcols, myrow, mycol, nprow, npcol, M, N, mnmin
      type(matrixblock)::block_rand(:)
      integer Bidxs
      stats%Flop_Tmp = 0
      block_off1 => ho_bf1%levels(level)%BP(bb_inv*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level)%BP(bb_inv*2)%LL(1)%matrices_block(1)
      block_inv => ho_bf1%levels(level)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)

      mm(1) = block_off1%M_loc
      mm(2) = block_off2%M_loc
      nn(1) = block_off1%N_loc
      nn(2) = block_off2%N_loc

      offM(1) = 0
      offM(2) = block_off1%M
      offN(1) = block_off1%M
      offN(2) = 0

      !!!*** redistribute Q and B^T from process layout of hodlr to the process layout of block_off1 and block_off2
      n1 = OMP_get_wtime()
      do bb = 1, 2
         block_off => ho_bf1%levels(level)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         ! if (mm(bb) > 0) then
            if (bb == 1) then
               allocate (matQ1(max(mm(bb),1), ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1)))
               matQ => matQ1
            endif
            if (bb == 2) then
               allocate (matQ2(max(mm(bb),1), ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1)))
               matQ => matQ2
            endif
         ! endif
         call Redistribute1Dto1D(Q, size(Q,1), block_inv%M_p, 0, block_inv%pgno, matQ, max(mm(bb),1), M_p, offM(bb), block_off%pgno, ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), ptree)

         ! if (nn(bb) > 0) then
            if (bb == 1) then
               allocate (matQcA_trans1(max(nn(bb),1), ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1)))
               matQcA_trans => matQcA_trans1
            endif
            if (bb == 2) then
               allocate (matQcA_trans2(max(nn(bb),1), ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1)))
               matQcA_trans => matQcA_trans2
            endif
         ! endif
         call Redistribute1Dto1D(QcA_trans, size(QcA_trans,1), block_inv%N_p, 0, block_inv%pgno, matQcA_trans, max(1,nn(bb)), N_p, offN(bb), block_off%pgno, ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), ptree)
      enddo
      n2 = OMP_get_wtime()
      stats%Time_RedistV = stats%Time_RedistV + n2 - n1

      !!!*** compute B^T=V^TS^TU^T and A = (QU)*(SV)
      do bb = 1, 2
         block_off => ho_bf1%levels(level)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         M = block_off%M
         N = block_off%N
         if (bb == 1) matQcA_trans => matQcA_trans1
         if (bb == 2) matQcA_trans => matQcA_trans2
         if (bb == 1) matQ => matQ1
         if (bb == 2) matQ => matQ2
         if (mm(bb) > 0) then
            call PQxSVDTruncate(block_rand(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), matQ, matQcA_trans, ranks(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), rank, option, stats, ptree, tolerance_abs, flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
         endif
         deallocate (matQcA_trans)
         deallocate (matQ)
      enddo

   end subroutine PQxSVDTruncate_twoforward

!!!!!***** this subroutine is part of the randomized HODLR_BF.
! The difference between this subroutine and BF_Resolving_Butterfly_LL_dat is that this subroutine requires redistribution of RandVectIn and RandVectOut to match the data layout of block_rand(bb_inv*2-1-Bidxs+1) and block_rand(bb_inv*2-Bidxs+1). Therefore this subroutine reconstructs two neighbouring butterflies together.
   subroutine BF_Resolving_Butterfly_LL_dat_twoforward(ho_bf1, level_c, num_vect_sub, nth_s, nth_e, Ng, level, Bidxs, bb_inv, block_rand, RandVectIn, RandVectOut, option, ptree, msh, stats)

      implicit none
      integer level, level_c, ii, bb, bb_inv
      DT :: RandVectIn(:, :), RandVectOut(:, :)
      DT, pointer :: mat1(:, :), mat2(:, :), mat(:, :), matQ1(:, :), matQ2(:, :), matQ(:, :), matQ2D(:, :), matQcA_trans1(:, :), matQcA_trans2(:, :), matQcA_trans(:, :), matQcA_trans2D(:, :), matQUt2D(:, :), UU(:, :), VV(:, :), mattemp(:, :), matOut1(:, :), matOut2(:, :), matOut(:, :), matIn1(:, :), matIn2(:, :), matIn(:, :)
      type(matrixblock), pointer::block_o, block_inv, block_schur, block_off1, block_off2, block_off
      integer groupn, groupm, mm(2), nn(2), ierr, nin1, nout1, nin2, nout2, offM(2), offN(2), offout1, offout2, rank
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::n1, n2, eps, flop
      integer, pointer::M_p(:, :), N_p(:, :)
      DTR, pointer::Singular(:)
      integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9), descQUt2D(9)
      integer tempi, ctxt, info, iproc, jproc, myi, myj, myArows, myAcols, myrow, mycol, nprow, npcol, M, N, mnmin
      type(matrixblock)::block_rand(:)
      integer Bidxs, num_vect_sub, nth_s, nth_e, Ng
      type(mesh)::msh

      stats%Flop_Tmp = 0
      block_off1 => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level_c)%BP(bb_inv*2)%LL(1)%matrices_block(1)
      block_inv => ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)

      mm(1) = block_off1%M_loc
      mm(2) = block_off2%M_loc
      nn(1) = block_off1%N_loc
      nn(2) = block_off2%N_loc

      offM(1) = 0
      offM(2) = block_off1%M
      offN(1) = block_off1%M
      offN(2) = 0

      !!!*** redistribute Q and B^T from process layout of hodlr to the process layout of block_off1 and block_off2
      n1 = OMP_get_wtime()
      do bb = 1, 2
         block_off => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         ! if (nn(bb) > 0) then
            if (bb == 1) then
               allocate (matOut1(max(1,nn(bb)), num_vect_sub))
               matOut => matOut1
            endif
            if (bb == 2) then
               allocate (matOut2(max(1,nn(bb)), num_vect_sub))
               matOut => matOut2
            endif
         ! endif
         call Redistribute1Dto1D(RandVectOut, block_inv%N_loc, block_inv%N_p, 0, block_inv%pgno, matOut, max(1,nn(bb)), N_p, offN(bb), block_off%pgno, num_vect_sub, ptree)

         ! if (mm(bb) > 0) then
            if (bb == 1) then
               allocate (matIn1(max(1,mm(bb)), num_vect_sub))
               matIn => matIn1
            endif
            if (bb == 2) then
               allocate (matIn2(max(1,mm(bb)), num_vect_sub))
               matIn => matIn2
            endif
         ! endif
         call Redistribute1Dto1D(RandVectIn, block_inv%M_loc, block_inv%M_p, 0, block_inv%pgno, matIn, max(1,mm(bb)),M_p, offM(bb), block_off%pgno, num_vect_sub, ptree)
      enddo
      n2 = OMP_get_wtime()
      stats%Time_RedistV = stats%Time_RedistV + n2 - n1

      !!!*** call BF_Resolving_Butterfly_LL_dat
      do bb = 1, 2
         block_off => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         M = block_off%M
         N = block_off%N
         if (bb == 1) matOut => matOut1
         if (bb == 2) matOut => matOut2
         if (bb == 1) matIn => matIn1
         if (bb == 2) matIn => matIn2
         if (mm(bb) > 0) then
            call BF_Resolving_Butterfly_LL_dat(num_vect_sub, nth_s, nth_e, Ng, level, block_rand(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), matIn, matOut, option, ptree, msh, stats)
         endif
         deallocate (matOut)
         deallocate (matIn)
      enddo
      stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

   end subroutine BF_Resolving_Butterfly_LL_dat_twoforward

!!!!!***** this subroutine is part of the randomized HODLR_BF.
! The difference between this subroutine and BF_Resolving_Butterfly_RR_dat is that this subroutine requires redistribution of RandVectIn and RandVectOut to match the data layout of block_rand(bb_inv*2-1-Bidxs+1) and block_rand(bb_inv*2-Bidxs+1). Therefore this subroutine reconstructs two neighbouring butterflies together.
   subroutine BF_Resolving_Butterfly_RR_dat_twoforward(ho_bf1, level_c, num_vect_sub, nth_s, nth_e, Ng, level, Bidxs, bb_inv, block_rand, RandVectIn, RandVectOut, option, ptree, msh, stats)

      implicit none
      integer level, level_c, ii, bb, bb_inv
      DT :: RandVectIn(:, :), RandVectOut(:, :)
      DT, pointer :: mat1(:, :), mat2(:, :), mat(:, :), matQ1(:, :), matQ2(:, :), matQ(:, :), matQ2D(:, :), matQcA_trans1(:, :), matQcA_trans2(:, :), matQcA_trans(:, :), matQcA_trans2D(:, :), matQUt2D(:, :), UU(:, :), VV(:, :), mattemp(:, :), matOut1(:, :), matOut2(:, :), matOut(:, :), matIn1(:, :), matIn2(:, :), matIn(:, :)
      type(matrixblock), pointer::block_o, block_inv, block_schur, block_off1, block_off2, block_off
      integer groupn, groupm, mm(2), nn(2), ierr, nin1, nout1, nin2, nout2, offM(2), offN(2), offout1, offout2, rank
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::n1, n2, eps, flop
      integer, pointer::M_p(:, :), N_p(:, :)
      DTR, pointer::Singular(:)
      integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9), descQUt2D(9)
      integer tempi, ctxt, info, iproc, jproc, myi, myj, myArows, myAcols, myrow, mycol, nprow, npcol, M, N, mnmin
      type(matrixblock)::block_rand(:)
      integer Bidxs, num_vect_sub, nth_s, nth_e, Ng
      type(mesh)::msh

      stats%Flop_Tmp = 0
      block_off1 => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level_c)%BP(bb_inv*2)%LL(1)%matrices_block(1)
      block_inv => ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)

      mm(1) = block_off1%M_loc
      mm(2) = block_off2%M_loc
      nn(1) = block_off1%N_loc
      nn(2) = block_off2%N_loc

      offM(1) = 0
      offM(2) = block_off1%M
      offN(1) = block_off1%M
      offN(2) = 0

      !!!*** redistribute Q and B^T from process layout of hodlr to the process layout of block_off1 and block_off2
      n1 = OMP_get_wtime()
      do bb = 1, 2
         block_off => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         ! if (mm(bb) > 0) then
            if (bb == 1) then
               allocate (matOut1(max(1,mm(bb)), num_vect_sub))
               matOut => matOut1
            endif
            if (bb == 2) then
               allocate (matOut2(max(1,mm(bb)), num_vect_sub))
               matOut => matOut2
            endif
         ! endif
         call Redistribute1Dto1D(RandVectOut, block_inv%M_loc, block_inv%M_p, 0, block_inv%pgno, matOut, max(1,mm(bb)), M_p, offM(bb), block_off%pgno, num_vect_sub, ptree)

         ! if (nn(bb) > 0) then
            if (bb == 1) then
               allocate (matIn1(max(nn(bb),1), num_vect_sub))
               matIn => matIn1
            endif
            if (bb == 2) then
               allocate (matIn2(max(nn(bb),1), num_vect_sub))
               matIn => matIn2
            endif
         ! endif
         call Redistribute1Dto1D(RandVectIn, block_inv%N_loc, block_inv%N_p, 0, block_inv%pgno, matIn, max(1,nn(bb)),N_p, offN(bb), block_off%pgno, num_vect_sub, ptree)
      enddo
      n2 = OMP_get_wtime()
      stats%Time_RedistV = stats%Time_RedistV + n2 - n1

      !!!*** call BF_Resolving_Butterfly_RR_dat
      do bb = 1, 2
         block_off => ho_bf1%levels(level_c)%BP(bb_inv*2 - 1 + bb - 1)%LL(1)%matrices_block(1)
         M_p => block_off%M_p
         N_p => block_off%N_p
         M = block_off%M
         N = block_off%N
         if (bb == 1) matOut => matOut1
         if (bb == 2) matOut => matOut2
         if (bb == 1) matIn => matIn1
         if (bb == 2) matIn => matIn2
         if (mm(bb) > 0) then
            call BF_Resolving_Butterfly_RR_dat(num_vect_sub, nth_s, nth_e, Ng, level, block_rand(bb_inv*2 - 1 + bb - 1 - Bidxs + 1), matIn, matOut, option, ptree, msh, stats)
         endif
         deallocate (matOut)
         deallocate (matIn)
      enddo
      stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

   end subroutine BF_Resolving_Butterfly_RR_dat_twoforward

end module BPACK_randomMVP
