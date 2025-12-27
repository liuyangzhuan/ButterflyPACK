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
!> @file Bplus_factor.f90
!> @brief Low-level subroutines for factorizing a BPACK (HMAT-LR/HMAT-BF/BLR/HODBF/HODLR/HSS-BF) matrix



#include "ButterflyPACK_config.fi"
module Bplus_factor
   use BPACK_DEFS
   use MISC_Utilities
   use Bplus_compress
   use Bplus_randomizedop


contains

   subroutine Full_LU(blocks, option, stats)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats

      integer size_m, size_n
      integer i, j, k, ii, jj, kk, zfpflag
      real*8 T0, T1, tol_used
      type(matrixblock) :: blocks
      real(kind=8) flop

      T0 = MPI_Wtime()
      size_m = size(blocks%fullmat, 1)
      if (option%ILU == 0) then

      zfpflag=0
      if(allocated(blocks%FullmatZFP%buffer_r))zfpflag=1
#if HAVE_ZFP
      if(zfpflag==1)call ZFP_Decompress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,0)
#endif
         ! do ii=1,size_m
         ! do jj=1,size_m
         ! write(777,*)dble(blocks%fullmat(ii,jj)),aimag(blocks%fullmat(ii,jj))
         ! enddo
         ! enddo
         ! call getrff90(blocks%fullmat, blocks%ipiv, flop=flop)
         call getrfmodf90(blocks%fullmat, option%jitter, blocks%ipiv, flop=flop,phase=blocks%phase,logabsdet=blocks%logabsdet)

         stats%Flop_Factor = stats%Flop_Factor + flop
         ! do ii=1,size_m
         ! do jj=1,size_m
         ! write(778,*)dble(blocks%fullmat(ii,jj)),aimag(blocks%fullmat(ii,jj))
         ! enddo
         ! enddo
#if HAVE_ZFP
         if(zfpflag==1)call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,option%tol_comp,0)
#endif
      else
         do ii = 1, size_m
            blocks%ipiv(ii) = ii
         enddo
      endif
      T1 = MPI_Wtime()
      stats%Time_Direct_LU = stats%Time_Direct_LU + T1 - T0
      return

   end subroutine Full_LU

   subroutine Full_add_multiply(block3, chara, block1, block2, h_mat, option, stats, ptree, msh)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Hmat)::h_mat
      type(proctree)::ptree
      type(mesh)::msh

      integer level_butterfly, flag
      integer i, j, k, ii, level, mm, nn, kk, rank, level_blocks, mn, group_k, rmax, rank_new
      integer style(3), data_type(3), id1, id2, id3, zfpflag
      character chara
      DT, allocatable::Vin(:, :), Vin1(:, :), Vin2(:, :), fullmat(:, :), fullmatrix(:, :), matrix_U(:,:), matrix_V(:,:)
      real*8 T0, T1,T2,T3, tol_used, flop_tmp
      type(matrixblock) :: block1, block2, block3

      stats%Flop_Tmp = 0

      T0 = MPI_Wtime()
      style(3) = block3%style
      style(1) = block1%style
      style(2) = block2%style
      level_blocks = block3%level

      group_k = block1%col_group
      kk = msh%basis_group(group_k)%tail - msh%basis_group(group_k)%head + 1


      if(style(3) == 1)then
         zfpflag=0
         if(allocated(block3%FullmatZFP%buffer_r))zfpflag=1
#if HAVE_ZFP
         if(zfpflag==1)call ZFP_Decompress(block3%fullmat,block3%FullmatZFP,block3%M,block3%N,tol_used,0)
#endif
      else
         call assert(style(1) == 1 .and. style(2) == 1, 'both block1 and block2 supposed to be style 1')
         zfpflag=0
         if(allocated(block1%FullmatZFP%buffer_r))zfpflag=1
      endif

         mm = block3%m
         nn = block3%n

if(style(3) == 1)then

         allocate (Vin(nn, nn))
         Vin = 0d0
         do ii = 1, nn
            Vin(ii, ii) = 1d0
         enddo
         allocate (Vin1(kk, nn))
         Vin1 = 0d0
         allocate (fullmatrix(mm, nn))
         fullmatrix = 0d0

         call Hmat_block_MVP_dat(block2, 'N', msh%basis_group(block2%row_group)%head, msh%basis_group(block2%col_group)%head, nn, Vin, nn, Vin1, kk, BPACK_cone, ptree, stats)
         call Hmat_block_MVP_dat(block1, 'N', msh%basis_group(block1%row_group)%head, msh%basis_group(block1%col_group)%head, nn, Vin1, kk, fullmatrix, mm, BPACK_cone, ptree, stats)

         if (chara == '-') fullmatrix = -fullmatrix
         block3%fullmat = block3%fullmat + fullmatrix
         deallocate (Vin)
         deallocate (Vin1)
         deallocate (fullmatrix)

#if HAVE_ZFP
         if(zfpflag==1)call ZFP_Compress(block3%fullmat,block3%FullmatZFP,block3%M,block3%N,option%tol_comp,0)
#endif

else
         T2 = MPI_Wtime()
#if HAVE_ZFP
         if(allocated(block1%FullmatZFP%buffer_r))call ZFP_Decompress(block1%fullmat,block1%FullmatZFP,block1%M,block1%N,tol_used,1)
         if(allocated(block2%FullmatZFP%buffer_r))call ZFP_Decompress(block2%fullmat,block2%FullmatZFP,block2%M,block2%N,tol_used,1)
#endif

         allocate (fullmatrix(mm, nn))
         fullmatrix = 0d0
         call gemmf90(block1%fullmat, block1%M, block2%fullmat, block2%M, fullmatrix, mm, 'N', 'N', mm, nn, kk, BPACK_cone, BPACK_czero, flop=flop_tmp)
         stats%Flop_Tmp = stats%Flop_Tmp + flop_tmp
         if (chara == '-') fullmatrix = -fullmatrix

#if HAVE_ZFP
         if(allocated(block1%FullmatZFP%buffer_r))call ZFP_Compress(block1%fullmat,block1%FullmatZFP,block1%M,block1%N,option%tol_comp,1)
         if(allocated(block2%FullmatZFP%buffer_r))call ZFP_Compress(block2%fullmat,block2%FullmatZFP,block2%M,block2%N,option%tol_comp,1)
#endif
         allocate (Vin2(nn, nn))
         Vin2 = 0d0
         do ii = 1, nn
            Vin2(ii, ii) = 1d0
         enddo
         allocate (fullmat(mm, nn))
         fullmat = 0d0
         call Hmat_block_MVP_dat(block3, 'N', msh%basis_group(block3%row_group)%head, msh%basis_group(block3%col_group)%head, nn, Vin2, nn, fullmat, mm, BPACK_cone, ptree, stats)
         fullmat = fullmat + fullmatrix


         rmax = min(mm, nn)
         allocate (matrix_U(mm, rmax))
         allocate (matrix_V(rmax, nn))
         call ACA_CompressionFull(fullmat, matrix_U, matrix_V, mm, nn, rmax, rank_new, option%tol_comp*0.3d0, option%tol_comp)

         deallocate(block3%ButterflyU%blocks(1)%matrix)
         allocate(block3%ButterflyU%blocks(1)%matrix(mm,rank_new))
         block3%ButterflyU%blocks(1)%matrix=matrix_U(:,1:rank_new)

         deallocate(block3%ButterflyV%blocks(1)%matrix)
         allocate(block3%ButterflyV%blocks(1)%matrix(nn,rank_new))
         call copymatT(matrix_V(1:rank_new,:), block3%ButterflyV%blocks(1)%matrix, rank_new, nn)

         block3%rankmax=rank_new
         block3%rankmin=rank_new
         deallocate(matrix_U)
         deallocate(matrix_V)

         deallocate(Vin2)
         deallocate(fullmat)
         deallocate (fullmatrix)
         T3 = MPI_Wtime()
         ! time_tmp = time_tmp + T3-T2
endif


      T1 = MPI_Wtime()
      stats%Time_Add_Multiply = stats%Time_Add_Multiply + T1 - T0
      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp

      return

   end subroutine Full_add_multiply

   subroutine Full_add(block3, chara, block1, ptree, stats, option)

      implicit none

      type(Hoption)::option
      integer level_butterfly, flag, group_n, group_m, zfpflag
      integer i, j, k, level, mm, nn, rank, level_blocks, mn, ii, jj
      integer style(3), data_type(3), id1, id2, id3
      character chara
      DT, allocatable:: Vin(:, :)
      real*8 T0, T1, tol_used
      type(matrixblock) :: block1, block3
      type(Hstat):: stats
      type(proctree):: ptree
      DT, allocatable::fullmatrix(:, :)

      stats%Flop_Tmp = 0

      style(3) = block3%style
      style(1) = block1%style
      level_blocks = block3%level

      T0 = MPI_Wtime()
      call assert(style(1) /= 1, 'block1 not supposed to be full')
      call assert(style(3) == 1, 'block3 supposed to be full')

      group_m = block3%row_group
      group_n = block3%col_group

      mm = block3%M
      nn = block3%N
      allocate (Vin(nn, nn))

      Vin = 0d0
      do ii = 1, nn
         Vin(ii, ii) = 1d0
      enddo

      allocate (fullmatrix(mm, nn))
      fullmatrix = 0d0

      call BF_block_MVP_dat(block1, 'N', mm, nn, nn, Vin, nn, fullmatrix, mm, BPACK_cone, BPACK_czero, ptree, stats)

      if (chara == '-') fullmatrix = -fullmatrix

      zfpflag=0
      if(allocated(block3%FullmatZFP%buffer_r))zfpflag=1

#if HAVE_ZFP
      if(zfpflag==1)call ZFP_Decompress(block3%fullmat,block3%FullmatZFP,block3%M,block3%N,tol_used,0)
#endif
      block3%fullmat = block3%fullmat + fullmatrix
#if HAVE_ZFP
      if(zfpflag==1)call ZFP_Compress(block3%fullmat,block3%FullmatZFP,block3%M,block3%N,option%tol_comp,0)
#endif
      deallocate (fullmatrix)
      deallocate (Vin)

      T1 = MPI_Wtime()
      stats%Time_Add_Multiply = stats%Time_Add_Multiply + T1 - T0
      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
      return

   end subroutine Full_add

   subroutine LR_minusBC(ho_bf1, level_c, rowblock, ptree, stats)




      implicit none

      integer level_c, rowblock
      integer i, j, k, level, num_blocks, num_row, num_col, ii, jj, kk, test
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm_diag
      character chara
      real(kind=8) a, b, c, d
      DT ctemp1, ctemp2
      type(matrixblock), pointer::block_o

      ! type(vectorsblock), pointer :: random1, random2

      DTR, allocatable :: Singular(:)
      integer idx_start_glo, N_diag, idx_start_diag, idx_start_loc, idx_end_loc
      ! DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),myA(:,:),BUold(:,:),BVold(:,:),CUold(:,:),CVold(:,:),BU(:,:),BV(:,:),CU(:,:),CV(:,:),BVCU(:,:),BUBVCU(:,:)

      integer Nsub, Ng, unique_nth, level_left_start, ll
      integer*8 idx_start
      integer level_blocks
      integer groupm_start, groupn_start, dimension_rank, rank1, rank
      integer header_mm, header_nn
      integer header_m, header_n, tailer_m, tailer_n

      real(kind=8)::n2, n1
      type(hobf)::ho_bf1
      type(matrixblock), pointer::block_off1, block_off2
      type(proctree)::ptree
      integer pgno, pgno1, pgno2
      integer descBUold(9), descBVold(9), descCUold(9), descCVold(9), descBU(9), descBV(9), descCU(9), descCV(9), descBVCU(9), descBUBVCU(9)
      integer ctxt1, ctxt2, ctxt, ctxtall, info, myrow, mycol, myArows, myAcols
      type(Hstat)::stats

      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)

      block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
      call BF_delete(block_o, 1)

      block_o%level_butterfly = 0
      block_o%ButterflyU%nblk_loc = 1
      block_o%ButterflyU%inc = 1
      block_o%ButterflyU%idx = 1
      block_o%ButterflyV%nblk_loc = 1
      block_o%ButterflyV%inc = 1
      block_o%ButterflyV%idx = 1

      allocate (block_o%ButterflyU%blocks(1))
      allocate (block_o%ButterflyV%blocks(1))

      pgno = block_o%pgno
      pgno1 = block_off1%pgno
      pgno2 = block_off2%pgno
      ! if(ptree%MyID==2)write(*,*)pgno,pgno1,pgno2,'nana'
      call assert(pgno == pgno1, 'block_o and block_off1 should be on the same process group')
      call assert(pgno == pgno2, 'block_o and block_off2 should be on the same process group')

      mm = block_off1%M
      nn = block_off1%N

      rank = block_off2%rankmax
      block_o%rankmax = rank
      block_o%M_loc = block_off1%M_loc
      block_o%N_loc = block_o%M_loc
      allocate (block_o%M_p(size(block_off1%M_p, 1), 2))
      block_o%M_p = block_off1%M_p
      allocate (block_o%N_p(size(block_off1%M_p, 1), 2))
      block_o%N_p = block_off1%M_p

      allocate (block_o%ButterflyU%blocks(1)%matrix(block_o%M_loc, rank))
      block_o%ButterflyU%blocks(1)%matrix = 0
      allocate (block_o%ButterflyV%blocks(1)%matrix(block_o%M_loc, rank))
      block_o%ButterflyV%blocks(1)%matrix = block_off2%ButterflyV%blocks(1)%matrix

      stats%Flop_Tmp = 0
      call BF_block_MVP_dat(block_off1, 'N', block_off1%M_loc, block_off1%N_loc, rank, block_off2%ButterflyU%blocks(1)%matrix, block_off1%N_loc, block_o%ButterflyU%blocks(1)%matrix, block_off1%M_loc, ctemp1, ctemp2, ptree, stats)
      block_o%ButterflyU%blocks(1)%matrix = -block_o%ButterflyU%blocks(1)%matrix
      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
      stats%Flop_Tmp = 0

      return

   end subroutine LR_minusBC

   subroutine LR_SMW(block_o, Memory, ptree, option, stats, pgno)



      implicit none

      integer level_c, rowblock, kover, rank, kk1, kk2, nprow, npcol
      integer i, j, k, level, num_blocks, blocks3, num_row, num_col, ii, jj, kk, level_butterfly, mm, nn
      integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn, index_j, index_i
      real(kind=8) a, b, c, d, Memory, flop, norm1
      DT ctemp, TEMP(1)
      type(matrixblock)::block_o
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), matrixtemp2(:, :), matrixtemp3(:, :), UU(:, :), VV(:, :), matrix_small(:, :), matrix_small1(:, :), vin(:, :), vout1(:, :), vout2(:, :), vout3(:, :), matU(:, :), matU2D(:, :), matU2D1(:, :), matV2D(:, :)
      DTR, allocatable:: Singular(:)
      integer, allocatable :: ipiv(:), iwork(:)
      type(proctree)::ptree
      integer pgno, ctxt, ctxt_head, myrow, mycol, myArows, myAcols, iproc, myi, jproc, myj, info
      integer descUV(9), descU(9), descV(9), descsmall(9), desctemp(9), TEMPI(1)

      integer lwork, liwork, lcmrc, ierr
      DT, allocatable:: work(:)
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::dlamch

      ctxt = ptree%pgrp(pgno)%ctxt
      ctxt_head = ptree%pgrp(pgno)%ctxt_head

      rank = size(block_o%ButterflyU%blocks(1)%matrix, 2)

      if (rank <= nbslpk) then
         allocate (matrixtemp(rank, rank))
         matrixtemp = 0
         allocate (matrixtemp1(rank, rank))
         matrixtemp1 = 0
         allocate (matU(block_o%M_loc, rank))
         matU = block_o%ButterflyU%blocks(1)%matrix

         call gemmf90(block_o%ButterflyV%blocks(1)%matrix, block_o%M_loc, block_o%ButterflyU%blocks(1)%matrix, block_o%M_loc, matrixtemp, rank, 'T', 'N', rank, rank, block_o%M_loc, BPACK_cone, BPACK_czero, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop

         ! write(*,*)'goog1'
         call assert(MPI_COMM_NULL /= ptree%pgrp(pgno)%Comm, 'communicator should not be null 1')
         call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*rank, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
         ! write(*,*)'goog2'
         do ii = 1, rank
            matrixtemp1(ii, ii) = matrixtemp1(ii, ii) + 1
         enddo

         ! write(*,*)abs(matrixtemp1),rank,'gggddd'
#if 1
         allocate (ipiv(rank))
         ipiv = 0
         ! call getrff90(matrixtemp1, ipiv, flop=flop)
         call getrfmodf90(matrixtemp1, option%jitter, ipiv, flop=flop,phase=block_o%phase,logabsdet=block_o%logabsdet)
         stats%Flop_Factor = stats%Flop_Factor + flop
         call getrif90(matrixtemp1, ipiv, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
         deallocate (ipiv)
#else
         matrixtemp = matrixtemp1
         call GeneralInverse(rank, rank, matrixtemp, matrixtemp1, BPACK_SafeEps, Flops=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
#endif

         call gemmf90(matU, block_o%M_loc, matrixtemp1, rank, block_o%ButterflyU%blocks(1)%matrix, block_o%M_loc, 'N', 'N', block_o%M_loc, rank, rank, BPACK_cone, BPACK_czero, flop=flop)
         block_o%ButterflyU%blocks(1)%matrix = -block_o%ButterflyU%blocks(1)%matrix
         stats%Flop_Factor = stats%Flop_Factor + flop

         deallocate (matrixtemp, matrixtemp1, matU)

      else
         call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(block_o%M, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (matU2D(max(1,myArows), max(1,myAcols)))
            matU2D = 0
            allocate (matU2D1(max(1,myArows), max(1,myAcols)))
            matU2D1 = 0
            call descinit_wp(descU, block_o%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descU')
            myArows = numroc_wp(block_o%N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (matV2D(max(1,myArows), max(1,myAcols)))
            matV2D = 0
            call descinit_wp(descV, block_o%N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descV')
         else
            allocate (matU2D(1, 1))
            allocate (matU2D1(1, 1))
            descU(2) = -1
            allocate (matV2D(1, 1))
            descV(2) = -1
         endif
         call Redistribute1Dto2D(block_o%ButterflyU%blocks(1)%matrix, block_o%M_p, 0, pgno, matU2D, block_o%M, 0, pgno, rank, ptree)
         call Redistribute1Dto2D(block_o%ButterflyV%blocks(1)%matrix, block_o%N_p, 0, pgno, matV2D, block_o%N, 0, pgno, rank, ptree)


         block_o%phase=1
         block_o%logabsdet=0
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (matrix_small(max(1,myArows), max(1,myAcols)))
            matrix_small = 0
            call descinit_wp(descsmall, rank, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descsmall')

            call pgemmf90('T', 'N', rank, rank, block_o%N, BPACK_cone, matV2D, 1, 1, descV, matU2D, 1, 1, descU, BPACK_czero, matrix_small, 1, 1, descsmall, flop=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
            do ii = 1, rank
               call g2l(ii, rank, nprow, nbslpk, iproc, myi)
               call g2l(ii, rank, npcol, nbslpk, jproc, myj)
               if (iproc == myrow .and. jproc == mycol) then
                  matrix_small(myi, myj) = matrix_small(myi, myj) + 1
               endif
            enddo
#if 1
            allocate (ipiv(max(1,myArows) + nbslpk))
            ipiv = 0
            ! call pgetrff90(rank, rank, matrix_small, 1, 1, descsmall, ipiv, info, flop=flop)
            call pgetrfmodf90(rank, rank, matrix_small, 1, 1, descsmall, option%jitter, ipiv, info, flop=flop, phase_loc=block_o%phase,logabsdet_loc=block_o%logabsdet)
            stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
            ! norm1 = pfnorm(rank, rank, matrix_small, 1, 1, descsmall, '1')
            ! do ii = 1, rank
            !    call g2l(ii, rank, nprow, nbslpk, iproc, myi)
            !    call g2l(ii, rank, npcol, nbslpk, jproc, myj)
            !    if (iproc == myrow .and. jproc == mycol) then
            !       if(abs(matrix_small(myi, myj))<dlamch('E'))then
            !          matrix_small(myi, myj) = sign(1d0,dble(matrix_small(myi, myj)))*sqrt(dlamch('E'))*norm1
            !       endif
            !    endif
            ! enddo
            call pgetrif90(rank, matrix_small, 1, 1, descsmall, ipiv, flop=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
            deallocate (ipiv)
#else
            allocate (matrix_small1(max(1,myArows), max(1,myAcols)))
            matrix_small1=0
            call PGeneralInverse(rank, rank, matrix_small, matrix_small1, BPACK_SafeEps, ctxt, Flops=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
            matrix_small=matrix_small1
            deallocate(matrix_small1)

#endif
            call pgemmf90('N', 'N', block_o%M, rank, rank, BPACK_cone, matU2D, 1, 1, descU, matrix_small, 1, 1, descsmall, BPACK_czero, matU2D1, 1, 1, descU, flop=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
            matU2D1 = -matU2D1
            deallocate (matrix_small)
         endif

         call MPI_ALLREDUCE(MPI_IN_PLACE, block_o%phase, 1, MPI_DT, MPI_PROD, ptree%pgrp(pgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, block_o%logabsdet, 1, MPI_DTR, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         call Redistribute2Dto1D(matU2D1, block_o%M, 0, pgno, block_o%ButterflyU%blocks(1)%matrix, block_o%M_p, 0, pgno, rank, ptree)

         deallocate (matU2D)
         deallocate (matU2D1)
         deallocate (matV2D)

      endif

      Memory = 0
      Memory = Memory + SIZEOF(block_o%ButterflyV%blocks(1)%matrix)/1024.0d3
      Memory = Memory + SIZEOF(block_o%ButterflyU%blocks(1)%matrix)/1024.0d3

      return

   end subroutine LR_SMW

   subroutine LR_Sblock(ho_bf1, level_c, rowblock, ptree, stats)




      implicit none

      integer level_c, rowblock
      integer i, j, k, level, num_blocks, num_row, num_col, ii, jj, kk, test, pp, qq
      integer mm, nn, mn, blocks1, blocks2, blocks3, level_butterfly, groupm, groupn, groupm_diag
      character chara
      real(kind=8) a, b, c, d
      type(matrixblock), pointer::block_o, blocks

      type(vectorsblock), pointer :: random1, random2

      DTR, allocatable :: Singular(:)
      integer idx_start_glo, N_diag, idx_start_diag, idx_end_diag, idx_start_loc, idx_end_loc
      DT, allocatable::vec_old(:, :), vec_new(:, :)

      integer Nsub, Ng, unique_nth, level_left_start
      integer*8 idx_start
      integer level_blocks, head, tail
      integer groupm_start, groupn_start, dimension_rank
      integer header_mm, header_nn
      integer header_m, header_n, tailer_m, tailer_n

      integer nth_s, nth_e, num_vect_sub, nth
      real(kind=8)::n2, n1
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(Hstat)::stats

      block_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

      ! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2)),'dfdU1',ptree%MyID
      ! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),'dfdV1',ptree%MyID

      level_butterfly = block_o%level_butterfly
      call assert(level_butterfly == 0, 'Butterfly_Sblock_LowRank only works with LowRank blocks')

      num_blocks = 2**level_butterfly

      num_vect_sub = size(block_o%ButterflyU%blocks(1)%matrix, 2)
      ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here

      ! get the right multiplied vectors
      pp = ptree%myid - ptree%pgrp(block_o%pgno)%head + 1
      idx_start_glo = block_o%headm + block_o%M_p(pp, 1) - 1

      ! mm=block_o%M
      mm = block_o%M_loc
      allocate (vec_old(mm, num_vect_sub))
      allocate (vec_new(mm, num_vect_sub))
      vec_old = block_o%ButterflyU%blocks(1)%matrix
      stats%Flop_Tmp = 0
      do level = ho_bf1%Maxlevel + 1, level_c + 1, -1
         N_diag = 2**(level - level_c - 1)
         idx_start_diag = max((rowblock - 1)*N_diag + 1, ho_bf1%levels(level)%Bidxs)
         idx_end_diag = min(rowblock*N_diag, ho_bf1%levels(level)%Bidxe)
         vec_new = 0

         n1 = MPI_Wtime()
         do ii = idx_start_diag, idx_end_diag

            if (associated(ho_bf1%levels(level)%BP_inverse(ii)%LL)) then
               blocks => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
               if (IOwnPgrp(ptree, blocks%pgno)) then

                  qq = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  head = blocks%headm + blocks%M_p(qq, 1) - 1
                  tail = head + blocks%M_loc - 1
                  idx_start_loc = head - idx_start_glo + 1
                  idx_end_loc = tail - idx_start_glo + 1
                  if (level == ho_bf1%Maxlevel + 1) then
                     call Full_block_MVP_dat(blocks, 'N', idx_end_loc - idx_start_loc + 1, num_vect_sub,&
               &vec_old(idx_start_loc, 1),mm, vec_new(idx_start_loc, 1),mm, BPACK_cone, BPACK_czero)
                  else
                     call BF_block_MVP_inverse_dat(ho_bf1, level, ii, 'N', idx_end_loc - idx_start_loc + 1, num_vect_sub, vec_old(idx_start_loc, 1), mm, vec_new(idx_start_loc, 1), mm, ptree, stats)
                  endif

               endif
            endif
         end do
         n2 = MPI_Wtime()
         ! time_tmp = time_tmp + n2 - n1

         vec_old = vec_new
      end do
      ! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
      block_o%ButterflyU%blocks(1)%matrix = vec_new
      ! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2)),'dfdU',ptree%MyID
      ! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),'dfdV',ptree%MyID
      deallocate (vec_old)
      deallocate (vec_new)

      stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp

      return

   end subroutine LR_Sblock

   subroutine LR_A_minusBDinvC(partitioned_block, ptree, option, stats)


      implicit none
      integer level, ii, num_vect_sub, mv, nv
      DT, allocatable :: V_tmp(:, :), V_tmp2(:, :), Vin_tmp(:, :), Vinter_tmp(:, :), Vout_tmp(:, :), U_tmp(:, :)
      DT :: ctemp1, ctemp2, a, b
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      integer groupn, groupm, mm, nn, rank, rank1, rank2
      type(matrixblock)::partitioned_block
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer ctxt, nprow, npcol, myrow, mycol, pgno, myArows, myAcols, iproc, myi, jproc, myj, info
      DT, allocatable:: matU(:, :), matV(:, :), matU2D(:, :), matV2D(:, :), matU2Dnew(:, :), matV2Dnew(:, :)
      integer::descsMatU2D(9), descsMatV2D(9), descsMatSml(9), descsMatSmlRR1(9), descsMatSmlRR2(9), descsMatU2Dnew(9), descsMatV2Dnew(9), descsUUSml(9), descsVVSml(9), descsUU_u(9), descsVV_u(9), descsUU_v(9), descsVV_v(9)
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), UUu(:, :), UUv(:, :), VVu(:, :), VVv(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:), Singularuv(:)
      real(kind=8)::flop
      integer ierr, mn1, mn2, jj, ranknew

      blocks_A => partitioned_block%sons(1, 1)
      blocks_B => partitioned_block%sons(1, 2)
      blocks_C => partitioned_block%sons(2, 1)
      blocks_D => partitioned_block%sons(2, 2)

      groupn = blocks_B%col_group    ! Note: row_group and col_group interchanged here
      nn = blocks_B%N_loc
      groupm = blocks_B%row_group    ! Note: row_group and col_group interchanged here
      mm = blocks_B%M_loc
      num_vect_sub = blocks_B%rankmax
      allocate (Vin_tmp(blocks_B%N_loc, num_vect_sub))
      Vin_tmp = blocks_B%ButterflyV%blocks(1)%matrix
      allocate (Vout_tmp(blocks_B%M_loc, num_vect_sub))
      Vout_tmp = 0
      allocate (Vinter_tmp(blocks_B%N_loc, num_vect_sub))
      Vinter_tmp = 0

      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      call BF_block_MVP_dat(blocks_D, 'T', nn, nn, num_vect_sub, Vin_tmp, blocks_B%N_loc, Vinter_tmp, blocks_B%N_loc, ctemp1, ctemp2, ptree, stats)
      Vinter_tmp = Vinter_tmp + Vin_tmp
      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      call BF_block_MVP_dat(blocks_C, 'T', nn, mm, num_vect_sub, Vinter_tmp, blocks_B%N_loc, Vout_tmp, blocks_B%M_loc, ctemp1, ctemp2, ptree, stats)
      rank = blocks_B%rankmax + blocks_A%rankmax
      allocate (matU(blocks_A%M_loc, rank))
      matU(:, 1:blocks_A%rankmax) = blocks_A%ButterflyU%blocks(1)%matrix
      matU(:, 1 + blocks_A%rankmax:rank) = blocks_B%ButterflyU%blocks(1)%matrix
      allocate (matV(blocks_A%N_loc, rank))
      matV(:, 1:blocks_A%rankmax) = blocks_A%ButterflyV%blocks(1)%matrix
      matV(:, 1 + blocks_A%rankmax:rank) = -Vout_tmp

      ! deallocate(blocks_A%ButterflyU%blocks(1)%matrix)
      ! deallocate(blocks_A%ButterflyV%blocks(1)%matrix)
      ! blocks_A%rankmax=rank
      ! blocks_A%rankmin=rank
      ! allocate(blocks_A%ButterflyU%blocks(1)%matrix(blocks_A%M_loc,rank))
      ! allocate(blocks_A%ButterflyV%blocks(1)%matrix(blocks_A%N_loc,rank))
      ! blocks_A%ButterflyU%blocks(1)%matrix=matU
      ! blocks_A%ButterflyV%blocks(1)%matrix=matV

! SVD recompression



      !!!!>**** generate 2D grid blacs quantities for matU
      pgno = blocks_A%pgno
      ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatU2D, blocks_A%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matU2D(max(1,myArows), max(1,myAcols)))
         matU2D = 0

      else
         descsMatU2D(2) = -1
         allocate (matU2D(1, 1))
         matU2D = 0
      endif
      !!!!>**** redistribution of input matrix
      call Redistribute1Dto2D(matU, blocks_A%M_p, 0, pgno, matU2D, blocks_A%M, 0, pgno, rank, ptree)

      !!!!>**** generate 2D grid blacs quantities for matV transpose
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatV2D, blocks_A%N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matV2D(max(1,myArows), max(1,myAcols)))
         matV2D = 0
      else
         descsMatV2D(2) = -1
         allocate (matV2D(1, 1))
         matV2D = 0
      endif
      !!!!>**** redistribution of input matrix

      call Redistribute1Dto2D(matV, blocks_A%N_p, 0, pgno, matV2D, blocks_A%N, 0, pgno, rank, ptree)

      if (myrow /= -1 .and. mycol /= -1) then

         ! GCC 9 was segfaulting when using PXGEQRF when it calls PXGEQR2 ->PXLARFG->PXSCAL
#if __GNUC__ < 9
! #if 0
         mn1 = min(blocks_A%M, rank)
         myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (tau_Q(max(1,myAcols)))
         call pgeqrff90(blocks_A%M, rank, matU2D, 1, 1, descsMatU2D, tau_Q, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatSmlRR1, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (RR1(max(1,myArows), max(1,myAcols)))
         RR1 = 0d0
         do myj = 1, myAcols
            call l2g(myj, mycol, mn1, npcol, nbslpk, jj)
            do myi = 1, myArows
               call l2g(myi, myrow, rank, nprow, nbslpk, ii)
               if (ii <= jj) RR1(myi, myj) = matU2D(myi, myj)
            enddo
         enddo
         call pun_or_gqrf90(ctxt, matU2D, tau_Q, blocks_A%M, mn1, mn1, descsMatU2D, 1, 1, flop=flop)
         deallocate (tau_Q)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

         mn2 = min(blocks_A%N, rank)
         myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (tau_Q(max(1,myAcols)))
         call pgeqrff90(blocks_A%N, rank, matV2D, 1, 1, descsMatV2D, tau_Q, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

         myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatSmlRR2, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (RR2(max(1,myArows), max(1,myAcols)))
         RR2 = 0d0
         do myj = 1, myAcols
            call l2g(myj, mycol, mn2, npcol, nbslpk, jj)
            do myi = 1, myArows
               call l2g(myi, myrow, rank, nprow, nbslpk, ii)
               if (ii <= jj) RR2(myi, myj) = matV2D(myi, myj)
            enddo
         enddo
         call pun_or_gqrf90(ctxt, matV2D, tau_Q, blocks_A%N, mn2, mn2, descsMatV2D, 1, 1, flop=flop)
         deallocate (tau_Q)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
         allocate (mattemp(max(1,myArows), max(1,myAcols)))
         mattemp = 0
         call descinit_wp(descsMatSml, mn1, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', mn1, mn2, rank, BPACK_cone, RR1, 1, 1, descsMatSmlRR1, RR2, 1, 1, descsMatSmlRR2, BPACK_czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(min(mn1, mn2), nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUUSml, mn1, min(mn1, mn2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUsml(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(min(mn1, mn2), nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVVSml, min(mn1, mn2), mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVsml(max(1,myArows), max(1,myAcols)))

         allocate (Singularsml(min(mn1, mn2)))
         call PSVD_Truncate(mn1, mn2, mattemp, descsMatSml, UUsml, VVsml, descsUUSml, descsVVSml, Singularsml, option%tol_rand, ranknew, ctxt, BPACK_SafeUnderflow, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (matU2Dnew(max(1,myArows), max(1,myAcols)))
         matU2Dnew = 0
         call descinit_wp(descsMatU2Dnew, blocks_A%M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', blocks_A%M, ranknew, mn1, BPACK_cone, matU2D, 1, 1, descsMatU2D, UUsml, 1, 1, descsUUSml, BPACK_czero, matU2Dnew, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            matU2Dnew(:, myj) = matU2Dnew(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (matV2Dnew(max(1,myArows), max(1,myAcols)))
         matV2Dnew = 0
         call descinit_wp(descsMatV2Dnew, blocks_A%N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', blocks_A%N, ranknew, mn2, BPACK_cone, matV2D, 1, 1, descsMatV2D, VVsml, 1, 1, descsVVSml, BPACK_czero, matV2Dnew, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         rank = ranknew

         deallocate (mattemp, RR1, UUsml, VVsml, Singularsml)
         deallocate (RR2)

#else
         mn1 = min(blocks_A%M, rank)
         myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn1, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUU_u, blocks_A%M, mn1, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUu(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVV_u, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVu(max(1,myArows), max(1,myAcols)))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(blocks_A%M, rank, matU2D, descsMatU2D, UUu, VVu, descsUU_u, descsVV_u, Singularuv, option%tol_rand, rank1, ctxt, BPACK_SafeUnderflow, flop=flop)
         do ii = 1, rank1
            call g2l(ii, rank1, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VVu(myi, :) = VVu(myi, :)*Singularuv(ii)
            endif
         enddo
         deallocate (Singularuv)

         mn2 = min(blocks_A%N, rank)
         myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUU_v, blocks_A%N, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUv(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVV_v, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVv(max(1,myArows), max(1,myAcols)))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(blocks_A%N, rank, matV2D, descsMatV2D, UUv, VVv, descsUU_v, descsVV_v, Singularuv, option%tol_rand, rank2, ctxt, BPACK_SafeUnderflow, flop=flop)
         do ii = 1, rank2
            call g2l(ii, rank2, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VVv(myi, :) = VVv(myi, :)*Singularuv(ii)
            endif
         enddo
         deallocate (Singularuv)

         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         allocate (mattemp(max(1,myArows), max(1,myAcols)))
         mattemp = 0
         call descinit_wp(descsMatSml, rank1, rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', rank1, rank2, rank, BPACK_cone, VVu, 1, 1, descsVV_u, VVv, 1, 1, descsVV_v, BPACK_czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(min(rank1, rank2), nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUUSml, rank1, min(rank1, rank2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUsml(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(min(rank1, rank2), nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVVSml, min(rank1, rank2), rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVsml(max(1,myArows), max(1,myAcols)))

         allocate (Singularsml(min(rank1, rank2)))
         call PSVD_Truncate(rank1, rank2, mattemp, descsMatSml, UUsml, VVsml, descsUUSml, descsVVSml, Singularsml, option%tol_rand, ranknew, ctxt, BPACK_SafeUnderflow, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (matU2Dnew(max(1,myArows), max(1,myAcols)))
         matU2Dnew = 0
         call descinit_wp(descsMatU2Dnew, blocks_A%M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', blocks_A%M, ranknew, rank1, BPACK_cone, UUu, 1, 1, descsUU_u, UUsml, 1, 1, descsUUSml, BPACK_czero, matU2Dnew, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            matU2Dnew(:, myj) = matU2Dnew(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (matV2Dnew(max(1,myArows), max(1,myAcols)))
         matV2Dnew = 0
         call descinit_wp(descsMatV2Dnew, blocks_A%N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', blocks_A%N, ranknew, rank2, BPACK_cone, UUv, 1, 1, descsUU_v, VVsml, 1, 1, descsVVSml, BPACK_czero, matV2Dnew, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         rank = ranknew

         deallocate (mattemp, UUsml, VVsml, Singularsml)
         deallocate (UUu, VVu, UUv, VVv)

#endif

      else
         allocate (matU2Dnew(1, 1))
         allocate (matV2Dnew(1, 1))
      endif
      call MPI_Bcast(rank, 1, MPI_INTEGER, Main_ID, ptree%pgrp(pgno)%Comm, ierr)

      deallocate (blocks_A%ButterflyU%blocks(1)%matrix)
      deallocate (blocks_A%ButterflyV%blocks(1)%matrix)
      blocks_A%rankmax = rank
      blocks_A%rankmin = rank
      allocate (blocks_A%ButterflyU%blocks(1)%matrix(blocks_A%M_loc, rank))
      allocate (blocks_A%ButterflyV%blocks(1)%matrix(blocks_A%N_loc, rank))

      call Redistribute2Dto1D(matU2Dnew, blocks_A%M, 0, pgno, blocks_A%ButterflyU%blocks(1)%matrix, blocks_A%M_p, 0, pgno, rank, ptree)
      call Redistribute2Dto1D(matV2Dnew, blocks_A%N, 0, pgno, blocks_A%ButterflyV%blocks(1)%matrix, blocks_A%N_p, 0, pgno, rank, ptree)

      deallocate (matU2D)
      deallocate (matV2D)
      deallocate (Vin_tmp, Vout_tmp, Vinter_tmp)
      deallocate (matU, matV)
      deallocate (matU2Dnew, matV2Dnew)
   end subroutine LR_A_minusBDinvC

   subroutine BF_inverse_schur_partitionedinverse(ho_bf1, level_c, rowblock, error_inout, option, stats, ptree, msh)




#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_c, rowblock
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
      integer num_col, num_row, level, mm, nn, ii, jj, tt, ll
      character chara
      real(kind=8) T0
      type(matrixblock), pointer::block_o, block_off1, block_off2
      integer rank_new_max, rank0
      real(kind=8):: rank_new_avr, error
      integer niter
      real(kind=8):: error_inout, rate, err_avr
      integer itermax, ntry
      real(kind=8):: n1, n2, Memory
      type(Hoption)::option
      type(Hstat)::stats
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(mesh)::msh
      integer pgno

      error_inout = 0

      block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2 - 1)%LL(1)%matrices_block(1)
      block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)

      block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
      block_o%level_butterfly = block_off1%level_butterfly
      level_butterfly = block_o%level_butterfly

      Memory = 0

      if (block_off1%level_butterfly == 0 .or. block_off2%level_butterfly == 0) then
         call LR_minusBC(ho_bf1, level_c, rowblock, ptree, stats)
      else
         ho_bf1%ind_lv = level_c
         ho_bf1%ind_bk = rowblock
         rank0 = max(block_off1%rankmax, block_off2%rankmax)
         rate = option%rankrate !1.2d0
         call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, ho_bf1, BF_block_MVP_inverse_minusBC_dat, error, 'minusBC', option, stats, ptree, msh)
         stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
         error_inout = max(error_inout, error)
      endif

      pgno = block_o%pgno

      n1 = MPI_Wtime()
      ! if(block_o%level==3)then
      call BF_get_rank(block_o, ptree)
      if (level_butterfly >= option%schulzlevel) then
         call BF_inverse_schulziteration_IplusButter(block_o, error, option, stats, ptree, msh)
      else
         call BF_inverse_partitionedinverse_IplusButter(block_o, level_butterfly, 0, option, error, stats, ptree, msh, pgno)
      endif

      error_inout = max(error_inout, error)

      n2 = MPI_Wtime()
      stats%Time_SMW = stats%Time_SMW + n2 - n1
      ! write(*,*)'I+B Inversion Time:',n2-n1

      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) write (*, '(A10,I5,A6,I3,A8,I3,A11,Es14.7)') 'OneL No. ', rowblock, ' rank:', block_o%rankmax, ' L_butt:', block_o%level_butterfly, ' error:', error_inout

      return

   end subroutine BF_inverse_schur_partitionedinverse

   subroutine BF_inverse_schulziteration_IplusButter(block_o, error_inout, option, stats, ptree, msh)




#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_c, rowblock
      integer groupm, blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
      integer num_col, num_row, level, mm, nn, ii, jj, tt, ll, pp
      character chara
      real(kind=8) T0
      type(matrixblock)::block_o, block_Xn
      integer rank_new_max, rank0, num_vect
      real(kind=8):: rank_new_avr, error
      integer niter
      real(kind=8):: error_inout
      DTR:: error_inout_tmp
      real(kind=8):: rate, err_avr
      integer itermax, ntry, converged
      real(kind=8):: n1, n2, Memory, memory_temp, norm1, norm2, scale_new, rr
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(schulz_operand)::schulz_op
      DT, allocatable::VecIn(:, :), VecOut(:, :), VecBuff(:, :)
      DT::ctemp1, ctemp2
      character(len=10)::iternumber
      integer ierr

      error_inout = 0
      level_butterfly = block_o%level_butterfly

      Memory = 0

      mm = block_o%M_loc

      num_vect = 1
      allocate (VecIn(mm, num_vect))
      VecIn = 0
      allocate (VecOut(mm, num_vect))
      VecOut = 0
      allocate (VecBuff(mm, num_vect))
      VecBuff = 0

      call BF_copy('N', block_o, schulz_op%matrices_block, memory_temp)
      call BF_copy('N', block_o, block_Xn, memory_temp)

      call BF_compute_schulz_init(schulz_op, option, ptree, stats, msh)

      itermax = 100
      converged = 0
      ! n1 = MPI_Wtime()
      do ii = 1, itermax

         write (iternumber, "(I4)") ii

         rank0 = block_Xn%rankmax

         rate = option%rankrate !1.2d0
         call BF_randomized(block_Xn%pgno, level_butterfly, rank0, rate, block_Xn, schulz_op, BF_block_MVP_schulz_dat, error, 'schulz iter'//TRIM(iternumber), option, stats, ptree, msh, operand1=ii)
         stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp


         scale_new = 1d0
         rr=1d0
         do pp=1,schulz_op%order-1
            rr = (1-schulz_op%scale)*rr
            scale_new = scale_new + rr
         enddo
         scale_new = schulz_op%scale*scale_new
         schulz_op%scale = scale_new


         ! test error

         ctemp1 = 1.0d0; ctemp2 = 0.0d0
         call RandomMat(mm, num_vect, min(mm, num_vect), VecIn, 1)
         ! XnR
         call BF_block_MVP_schulz_Xn_dat(schulz_op, block_Xn, 'N', mm, mm, num_vect, VecIn, mm, VecBuff, mm, ctemp1, ctemp2, ptree, stats, ii + 1)

         ! AXnR
         call BF_block_MVP_dat(schulz_op%matrices_block, 'N', mm, mm, num_vect, VecBuff, mm, VecOut, mm, ctemp1, ctemp2, ptree, stats)
         VecOut = VecBuff + VecOut

         norm1 = fnorm(VecOut - VecIn, mm, num_vect)**2d0
         norm2 = fnorm(VecIn, mm, num_vect)**2d0
         call MPI_ALLREDUCE(MPI_IN_PLACE, norm1, 1, MPI_double_precision, MPI_SUM, ptree%pgrp(schulz_op%matrices_block%pgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, norm2, 1, MPI_double_precision, MPI_SUM, ptree%pgrp(schulz_op%matrices_block%pgno)%Comm, ierr)
         error_inout = sqrt(norm1)/sqrt(norm2)

         if (ptree%MyID == Main_ID .and. option%verbosity >= 1) write (*, '(A22,A6,I3,A8,I2,A8,I3,A7,Es14.7)') ' Schultz ', ' rank:', block_Xn%rankmax, ' Iter:', ii, ' L_butt:', block_Xn%level_butterfly, ' error:', error_inout

         if (error_inout < option%tol_rand) then
            converged = 1
            exit
         endif
         error_inout_tmp=error_inout
         if (myisnan(error_inout_tmp)) then
            converged = 0
            exit
         endif

      enddo
      ! n2 = MPI_Wtime()

      if (converged == 0) then
         write (*, *) 'Schulz Iteration does not converge, consider:'
         write (*, *) '  1: decrease tol_Rdetect such that the randomized contruction in each iteration is accurate'
         write (*, *) '  2: decrease schulzsplitlevel such that initial guess is good enough'
         stop
      else
         ! write(*,*)'Schulz Iteration Time:',n2-n1
         call BF_delete(block_o, 1)
         call BF_get_rank(block_Xn, ptree)
         rank_new_max = block_Xn%rankmax
         call BF_copy_delete(block_Xn, block_o, Memory)
         call BF_delete(schulz_op%matrices_block, 1)
         if(option%schulzhardstart==0)then
         if (allocated(schulz_op%diags)) deallocate (schulz_op%diags)
            do ii=schulz_op%bdiags%Bidxs,schulz_op%bdiags%Bidxe
               call BF_delete(schulz_op%bdiags%BF_inverse(ii),1)
            enddo
            deallocate(schulz_op%bdiags%BF_inverse)
         endif
      endif

      deallocate (VecIn)
      deallocate (VecOut)
      deallocate (VecBuff)

      return

   end subroutine BF_inverse_schulziteration_IplusButter



   recursive subroutine BF_bdiag_approximation_precompute(recurlevel, bidx,bdiags, option, stats, ptree, msh, pgno)
#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_c, rowblock,bidx, Maxgrp
      type(bdiag)::bdiags
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
      integer num_col, num_row, recurlevel, mm, nn, ii, jj, tt, kk1, kk2, rank, err_cnt
      character chara
      real(kind=8) T0, err_avr
      integer rank_new_max, rank0
      real(kind=8):: rank_new_avr, error, rate, norm
      integer niter
      real(kind=8):: error_inout
      integer itermax, ntry
      real(kind=8):: n1, n2, Memory
      DT, allocatable::matrix_small(:, :),vecin(:,:),vecout(:,:)
      type(matrixblock)::block_rand
      type(Hoption)::option
      type(Hstat)::stats
      integer level_butterfly_target, pgno, pgno1,pgnoA,pgnoD
      type(proctree)::ptree
      type(mesh)::msh
      integer ierr
      type(matrixblock)::partitioned_block

      if(recurlevel==bdiags%splitlevel)then
         bdiags%bidxs = min(bdiags%bidxs,bidx)
         bdiags%bidxe = max(bdiags%bidxe,bidx)
      else
         !>*** split process groups row-wise
         Maxgrp = 2**(ptree%nlevel) - 1
         pgnoA = pgno
         pgnoD = pgno
         if (pgnoA*2 <= Maxgrp) then
            pgnoA = pgno*2
            pgnoD = pgno*2+1
         endif

         if (IOwnPgrp(ptree, pgnoA)) then
            call BF_bdiag_approximation_precompute(recurlevel+1, bidx*2-1, bdiags, option, stats, ptree, msh, pgnoA)
         endif
         if (IOwnPgrp(ptree, pgnoD)) then
            call BF_bdiag_approximation_precompute(recurlevel+1, bidx*2, bdiags, option, stats, ptree, msh, pgnoD)
         endif
      endif

   end subroutine BF_bdiag_approximation_precompute




   recursive subroutine BF_bdiag_approximation(recurlevel, bidx, blocks_io,bdiags, option, stats, ptree, msh, pgno)
#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_c, rowblock,bidx
      type(bdiag)::bdiags
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
      integer num_col, num_row, recurlevel, mm, nn, ii, jj, tt, kk1, kk2, rank, err_cnt
      character chara
      real(kind=8) T0, err_avr
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      type(matrixblock)::blocks_io
      type(matrixblock)::blocks_schur
      integer rank_new_max, rank0
      real(kind=8):: rank_new_avr, error, rate, norm
      integer niter
      real(kind=8):: error_inout
      integer itermax, ntry
      real(kind=8):: n1, n2, Memory
      DT, allocatable::matrix_small(:, :),vecin(:,:),vecout(:,:)
      type(matrixblock)::block_rand
      type(Hoption)::option
      type(Hstat)::stats
      integer level_butterfly_target, pgno, pgno1
      type(proctree)::ptree
      type(mesh)::msh
      integer ierr
      type(matrixblock)::partitioned_block

      if(recurlevel==bdiags%splitlevel)then
         call BF_copy('N',blocks_io,bdiags%BF_inverse(bidx))
         call BF_inverse_partitionedinverse_IplusButter(bdiags%BF_inverse(bidx), blocks_io%level_butterfly, 0, option, error_inout, stats, ptree, msh, pgno)
      else

         allocate (partitioned_block%sons(2, 2))

         blocks_A => partitioned_block%sons(1, 1)
         blocks_B => partitioned_block%sons(1, 2)
         blocks_C => partitioned_block%sons(2, 1)
         blocks_D => partitioned_block%sons(2, 2)

         ! split into four smaller butterflies
         n1 = MPI_Wtime()
         call BF_split(blocks_io, partitioned_block, ptree, stats, msh, option,1)
         n2 = MPI_Wtime()
         stats%Time_split = stats%Time_split + n2 - n1

         if (IOwnPgrp(ptree, blocks_A%pgno)) then
            call BF_bdiag_approximation(recurlevel+1, bidx*2-1, blocks_A,bdiags, option, stats, ptree, msh, blocks_A%pgno)
         endif
         if (IOwnPgrp(ptree, blocks_D%pgno)) then
            call BF_bdiag_approximation(recurlevel+1, bidx*2, blocks_D,bdiags, option, stats, ptree, msh, blocks_D%pgno)
         endif

         do ii = 1, 2
         do jj = 1, 2
            call BF_delete(partitioned_block%sons(ii, jj), 1)
         enddo
         enddo
         deallocate (partitioned_block%sons)
      endif

   end subroutine BF_bdiag_approximation







   subroutine BF_compute_schulz_init(schulz_op, option, ptree, stats, msh)
#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_butterfly
      integer mm, nn, mn, ii
      real(kind=8) T0

      real(kind=8):: error
      integer niter, groupm, groupn, Nblock
      real(kind=8):: error_inout
      integer num_vect, rank, ranktmp, q, qq
      real(kind=8):: n1, n2, memory_temp, flop
      type(Hoption)::option
      type(schulz_operand)::schulz_op
      DTR, allocatable:: Singular(:)
      DT, allocatable::UU(:, :), VV(:, :), RandVectIn(:, :), RandVectOut(:, :), matrixtmp(:, :), matrixtmp1(:, :)
      type(proctree)::ptree
      type(Hstat)::stats
      type(mesh)::msh

      stats%Flop_tmp = 0

      schulz_op%order = option%schulzorder
      schulz_op%hardstart = option%schulzhardstart

      error_inout = 0

      level_butterfly = schulz_op%matrices_block%level_butterfly

      if(option%schulzhardstart==0)then
         schulz_op%bdiags%splitlevel = min(level_butterfly,option%schulzsplitlevel)
         Nblock = 2**schulz_op%bdiags%splitlevel
         schulz_op%bdiags%bidxs=10000000
         schulz_op%bdiags%bidxe=-10000000
         call BF_bdiag_approximation_precompute(0, 1, schulz_op%bdiags, option, stats, ptree, msh, schulz_op%matrices_block%pgno)
         allocate(schulz_op%bdiags%BF_inverse(schulz_op%bdiags%bidxs:schulz_op%bdiags%bidxe))
         call BF_bdiag_approximation(0, 1, schulz_op%matrices_block, schulz_op%bdiags, option, stats, ptree, msh, schulz_op%matrices_block%pgno)
      else
         mm = schulz_op%matrices_block%M_loc
         nn = mm
         num_vect = min(10, nn)

         allocate (RandVectIn(nn, num_vect))
         allocate (RandVectOut(mm, num_vect))
         RandVectOut = 0
         call RandomMat(nn, num_vect, min(nn, num_vect), RandVectIn, 1)

         ! computation of AR
         call BF_block_MVP_dat(schulz_op%matrices_block, 'N', mm, nn, num_vect, RandVectIn, nn,RandVectOut, mm, BPACK_cone, BPACK_czero, ptree, stats)
         RandVectOut = RandVectIn + RandVectOut
         ! computation of range Q of AR
         call PComputeRange(schulz_op%matrices_block%M_p, num_vect, RandVectOut, ranktmp, BPACK_SafeEps, ptree, schulz_op%matrices_block%pgno, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop

         ! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
         q = 6
         do qq = 1, q
            RandVectOut = conjg(cmplx(RandVectOut, kind=8))
            call BF_block_MVP_dat(schulz_op%matrices_block, 'T', mm, nn, num_vect, RandVectOut, mm, RandVectIn, nn, BPACK_cone, BPACK_czero, ptree, stats)
            RandVectIn = RandVectOut + RandVectIn
            RandVectIn = conjg(cmplx(RandVectIn, kind=8))
            call PComputeRange(schulz_op%matrices_block%M_p, num_vect, RandVectIn, ranktmp, BPACK_SafeEps, ptree, schulz_op%matrices_block%pgno, flop)

            call BF_block_MVP_dat(schulz_op%matrices_block, 'N', mm, nn, num_vect, RandVectIn, nn, RandVectOut, mm, BPACK_cone, BPACK_czero, ptree, stats)
            RandVectOut = RandVectIn + RandVectOut
            call PComputeRange(schulz_op%matrices_block%M_p, num_vect, RandVectOut, ranktmp, BPACK_SafeEps, ptree, schulz_op%matrices_block%pgno, flop)
         enddo


         ! computation of B = Q^c*A
         RandVectOut = conjg(cmplx(RandVectOut, kind=8))
         call BF_block_MVP_dat(schulz_op%matrices_block, 'T', mm, nn, num_vect, RandVectOut, mm,RandVectIn, nn,BPACK_cone, BPACK_czero, ptree, stats)
         RandVectIn = RandVectOut + RandVectIn
         RandVectOut = conjg(cmplx(RandVectOut, kind=8))

         ! computation of singular values of B
         mn = min(schulz_op%matrices_block%M, ranktmp)
         allocate (Singular(mn))
         Singular = 0
         call PSVDTruncateSigma(schulz_op%matrices_block, RandVectIn, ranktmp, rank, Singular, option, stats, ptree, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         schulz_op%A2norm = Singular(1)
         deallocate (Singular)

         deallocate (RandVectIn)
         deallocate (RandVectOut)
         stats%Flop_factor = stats%Flop_tmp
      endif



   end subroutine BF_compute_schulz_init

   recursive subroutine BF_inverse_partitionedinverse_IplusButter(blocks_io, level_butterfly_target, recurlevel, option, error_inout, stats, ptree, msh, pgno)




#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none

      integer level_c, rowblock
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
      integer num_col, num_row, recurlevel, mm, nn, ii, jj, tt, kk1, kk2, rank, err_cnt
      character chara
      real(kind=8) T0, err_avr
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      type(matrixblock)::blocks_io
      type(matrixblock)::blocks_schur
      integer rank_new_max, rank0
      real(kind=8):: rank_new_avr, error, rate, norm
      integer niter
      real(kind=8):: error_inout
      integer itermax, ntry
      real(kind=8):: n1, n2, Memory
      DT, allocatable::matrix_small(:, :),vecin(:,:),vecout(:,:)
      type(matrixblock)::block_rand
      type(Hoption)::option
      type(Hstat)::stats
      integer level_butterfly_target, pgno, pgno1
      type(proctree)::ptree
      type(mesh)::msh
      integer ierr
      type(matrixblock)::partitioned_block

      error_inout = 0

      if (blocks_io%level_butterfly == 0) then
         call LR_SMW(blocks_io, Memory, ptree, option, stats, pgno)
         return
      else
         allocate (partitioned_block%sons(2, 2))

         blocks_A => partitioned_block%sons(1, 1)
         blocks_B => partitioned_block%sons(1, 2)
         blocks_C => partitioned_block%sons(2, 1)
         blocks_D => partitioned_block%sons(2, 2)

         ! split into four smaller butterflies
         n1 = MPI_Wtime()
         call BF_split(blocks_io, partitioned_block, ptree, stats, msh, option)
         n2 = MPI_Wtime()
         stats%Time_split = stats%Time_split + n2 - n1

         if (IOwnPgrp(ptree, blocks_D%pgno)) then

            ! partitioned inverse of D
            ! level_butterfly=level_butterfly_target-1
            level_butterfly = blocks_D%level_butterfly
            pgno1 = blocks_D%pgno
            call BF_inverse_partitionedinverse_IplusButter(blocks_D, level_butterfly, recurlevel + 1, option, error, stats, ptree, msh, pgno1)
            error_inout = max(error_inout, error)

            ! construct the schur complement A-BD^-1C
            ! level_butterfly = level_butterfly_target-1
            level_butterfly = blocks_A%level_butterfly

            ! write(*,*)'A-BDC',level_butterfly,level
            call BF_get_rank_ABCD(partitioned_block, rank0)
            if (level_butterfly == 0) then
               ! n1 = MPI_Wtime()
               call LR_A_minusBDinvC(partitioned_block, ptree, option, stats)
               ! n2 = MPI_Wtime()
               ! time_tmp1 = time_tmp1 + n2-n1
            else
               rate = option%rankrate !1.2d0
               ! if(option%format==3)option%tol_Rdetect = option%tol_Rdetect/max(1,blocks_A%level_butterfly/2)
               call BF_randomized(blocks_A%pgno, level_butterfly, rank0, rate, blocks_A, partitioned_block, BF_block_MVP_inverse_A_minusBDinvC_dat, error, 'A-BD^-1C', option, stats, ptree, msh)
               ! if(option%format==3)option%tol_Rdetect = option%tol_Rdetect*max(1,blocks_A%level_butterfly/2)
               stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
               error_inout = max(error_inout, error)
            endif

            ! write(*,*)'ddd1'
            ! partitioned inverse of the schur complement
            ! level_butterfly=level_butterfly_target-1
            level_butterfly = blocks_A%level_butterfly
            pgno1 = blocks_D%pgno
            call BF_inverse_partitionedinverse_IplusButter(blocks_A, level_butterfly, recurlevel + 1, option, error, stats, ptree, msh, pgno1)
            error_inout = max(error_inout, error)
            call BF_get_rank_ABCD(partitioned_block, rank0)
         else
            rank0 = 0
            blocks_A%phase=1
            blocks_A%logabsdet=0
            blocks_D%phase=1
            blocks_D%logabsdet=0
         endif
         call MPI_ALLREDUCE(MPI_IN_PLACE, rank0, 1, MPI_integer, MPI_MAX, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, error_inout, 1, MPI_double_precision, MPI_MAX, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
         call MPI_Bcast(blocks_A%phase, 1, MPI_DT, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
         call MPI_Bcast(blocks_A%logabsdet, 1, MPI_DTR, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
         call MPI_Bcast(blocks_D%phase, 1, MPI_DT, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
         call MPI_Bcast(blocks_D%logabsdet, 1, MPI_DTR, Main_ID, ptree%pgrp(blocks_io%pgno)%Comm, ierr)


         if (level_butterfly_target <= 0) then ! try to use deterministic algorithms for merging four LRs into a bigger LR
            pgno1 = blocks_io%pgno  ! recording pgno as blocks_io%pgno will be set to partitioned_block%sons(1,1)%pgno in LR_ABCDinverse
            blocks_io%level_butterfly = 0
            call LR_ABCDinverse(partitioned_block, blocks_io, ptree, stats, option, msh)
            call BF_ReDistribute_Inplace(blocks_io, pgno1, stats, ptree, msh)
         else
            level_butterfly = level_butterfly_target
            rate = option%rankrate !1.2d0

            !>**** Check the estimated norm of the operator, if too small skip the randomized construction
            allocate(vecin(blocks_io%N_loc,1))
            vecin=1/sqrt(dble(blocks_io%N))
            allocate(vecout(blocks_io%M_loc,1))
            vecout=0
            call BF_block_MVP_inverse_ABCD_dat(partitioned_block, blocks_io, 'N', blocks_io%M_loc, blocks_io%N_loc, 1, vecin, blocks_io%N_loc, vecout, blocks_io%M_loc, BPACK_cone, BPACK_czero, ptree, stats)
            norm = fnorm(vecout,blocks_io%M_loc,1)**2d0
            call MPI_ALLREDUCE(MPI_IN_PLACE, norm, 1, MPI_double_precision, MPI_SUM, ptree%pgrp(blocks_io%pgno)%Comm, ierr)
            norm = sqrt(norm)
            deallocate(vecin)
            deallocate(vecout)

            if(norm<option%tol_Rdetect*1d-3)then
            ! if(norm<1d-12)then
               call BF_Zero(level_butterfly, blocks_io%row_group, blocks_io%col_group, blocks_io, block_rand, msh, ptree, option)

               call BF_delete(blocks_io, 1)
               call BF_get_rank(block_rand, ptree)
               ! rank_new_max = block_rand%rankmax
               call BF_copy_delete(block_rand, blocks_io)

            else
               ! if(option%format==3)option%tol_Rdetect = option%tol_Rdetect/max(1,blocks_io%level_butterfly/2)
               call BF_randomized(blocks_io%pgno, level_butterfly, rank0, rate, blocks_io, partitioned_block, BF_block_MVP_inverse_ABCD_dat, error, 'ABCDinverse', option, stats, ptree, msh)
               ! if(option%format==3)option%tol_Rdetect = option%tol_Rdetect*max(1,blocks_io%level_butterfly/2)

            endif

            stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
            error_inout = max(error_inout, error)
         endif

         ! stop

         blocks_io%phase = blocks_A%phase*blocks_D%phase
         blocks_io%logabsdet = blocks_A%logabsdet+blocks_D%logabsdet

         if (option%verbosity >= 2 .and. recurlevel == 0 .and. ptree%MyID == Main_ID) write (*, '(A23,A6,I3,A8,I3,A11,Es14.7)') ' RecursiveI ', ' rank:', blocks_io%rankmax, ' L_butt:', blocks_io%level_butterfly, ' error:', error_inout

         do ii = 1, 2
         do jj = 1, 2
            call BF_delete(partitioned_block%sons(ii, jj), 1)
         enddo
         enddo
         deallocate (partitioned_block%sons)

         return

      end if
   end subroutine BF_inverse_partitionedinverse_IplusButter

!>**** Merge four child LR into a bigger one
   !partitioned_block: (input) partitioned_block%sons store the four smaller BF
   !ptree: (input) process tree
   !stats: (input) statistics
   !option: (input) containing compression options
   !msh: (input) containing mesh
   !blocks_o: (inout) the parent BF
   subroutine LR_ABCDinverse(partitioned_block, blocks_o, ptree, stats, option, msh)

      implicit none
      type(matrixblock)::partitioned_block
      type(matrixblock)::blocks_o
      integer rank, ranktmp, leafsize
      integer header_m, header_n
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, mn
      type(mesh)::msh
      type(Hoption)::option
      type(proctree)::ptree
      integer pgno
      type(grid), pointer::gd
      integer:: cridx, info
      integer::mrange_dummy(1), nrange_dummy(1)
      type(Hstat)::stats
      integer::passflag = 0
      integer::frow, rmax
      real(kind=8)::error
      DT:: mat_dummy(1, 1)

      call BF_delete(blocks_o, 1)

      blocks_o%pgno = partitioned_block%sons(1, 1)%pgno
      call ComputeParallelIndices(blocks_o, blocks_o%pgno, ptree, msh)

      call LR_BuildABCD(blocks_o, partitioned_block, option, msh, stats, ptree, partitioned_block%sons(1, 1)%pgno, 0)

      stats%Flop_tmp=0
      call LR_HMerge(blocks_o, rank, option, msh, stats, ptree, partitioned_block%sons(1, 1)%pgno, 0, 0)
      stats%Flop_Factor=stats%Flop_Factor+stats%Flop_tmp
   end subroutine LR_ABCDinverse

!>**** Use low-rank arithmetic to form the four child LRs in inverse_ABCD
   !partitioned_block: (input) partitioned_block%sons store the four smaller LR
   !ptree: (input) process tree
   !stats: (input) statistics
   !option: (input) containing compression options
   !msh: (input) containing mesh
   !blocks: (inout) the current LR in the recursion
   !pgno: (in) the process group used for the four children
   !gd: (in) the process grid from process group pgno
   recursive subroutine LR_BuildABCD(blocks, partitioned_block, option, msh, stats, ptree, pgno, cridx)

      implicit none
      type(matrixblock)::partitioned_block
      integer rank, ranktmp
      integer header_m, header_n
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, rmax, mn
      type(mesh)::msh
      type(Hoption)::option
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno
      integer:: cridx, info
      DT, allocatable:: UU(:, :), UU1(:, :), UU2(:, :), VV(:, :), VV1(:, :), VV2(:, :), SS1(:, :), TT1(:, :), matU(:, :), matV(:, :), matU1(:, :), matV1(:, :), matU2(:, :), matV2(:, :), tmp(:, :), matU1D(:, :), matV1D(:, :), Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Fullmat(:, :), QQ1(:, :), matU2D(:, :), matV2D(:, :)
      DTR, allocatable::Singular(:)
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, M1, N1, M2, N2, rank1, rank2, ierr
      integer::descsmatU(9), descsmatV(9), descsmatU1(9), descsmatV1(9), descsmatU2(9), descsmatV2(9), descUU(9), descVV(9), descsmatU1c(9), descsmatU2c(9), descsmatV1c(9), descsmatV2c(9), descButterflyV(9), descButterflyU(9), descButterU1D(9), descButterV1D(9), descVin(9), descVout(9), descVinter(9), descFull(9)
      integer dims(6), dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
      DT:: TEMP(1)
      integer LWORK, mnmax, mnmin, rank_new
      DT, allocatable:: WORK(:)
      real(kind=8), allocatable::RWORK(:), center(:)
      real(kind=8):: rtemp, dist, error, rtemp1, rtemp0, fnorm1, fnorm0, flop
      integer :: nb1Dc, nb1Dr, frow, Dimn, edge_n, edge_m, MyID, Ntest, rmaxc, rmaxr
      integer, allocatable::M_p(:, :), N_p(:, :), mrange(:), nrange(:)
      type(Hstat)::stats
      integer::passflag = 0
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      integer ctxt

      blocks_A => partitioned_block%sons(1, 1)
      blocks_B => partitioned_block%sons(1, 2)
      blocks_C => partitioned_block%sons(2, 1)
      blocks_D => partitioned_block%sons(2, 2)

      rank = 0
      blocks%ButterflyU%idx = 1
      blocks%ButterflyU%inc = 1
      blocks%ButterflyU%nblk_loc = 1
      blocks%ButterflyV%idx = 1
      blocks%ButterflyV%inc = 1
      blocks%ButterflyV%nblk_loc = 1
      allocate (blocks%ButterflyU%blocks(1))
      allocate (blocks%ButterflyV%blocks(1))

      if (cridx == 2) then ! reach bottom level
         blocks%level = blocks_A%level

         if (blocks%row_group == blocks_A%row_group .and. blocks%col_group == blocks_A%col_group) then
            rank = blocks_A%rankmax
            allocate (UU(blocks_A%M_loc, rank))
            UU = blocks_A%ButterflyU%blocks(1)%matrix
            allocate (VV(blocks_A%N_loc, rank))
            VV = blocks_A%ButterflyV%blocks(1)%matrix
         elseif (blocks%row_group == blocks_B%row_group .and. blocks%col_group == blocks_B%col_group) then
            rank = blocks_B%rankmax
            allocate (UU(blocks_B%M_loc, rank))
            UU = 0
            call BF_block_MVP_dat(blocks_A, 'N', blocks_A%M_loc, blocks_A%N_loc, rank, blocks_B%ButterflyU%blocks(1)%matrix, blocks_A%N_loc, UU, blocks_A%M_loc, BPACK_cone, BPACK_czero, ptree, stats)
            UU = -UU - blocks_B%ButterflyU%blocks(1)%matrix
            allocate (VV(blocks_B%N_loc, rank))
            VV = 0
            call BF_block_MVP_dat(blocks_D, 'T', blocks_D%M_loc, blocks_D%N_loc, rank, blocks_B%ButterflyV%blocks(1)%matrix, blocks_D%M_loc, VV, blocks_D%N_loc, BPACK_cone, BPACK_czero, ptree, stats)
            VV = VV + blocks_B%ButterflyV%blocks(1)%matrix

         elseif (blocks%row_group == blocks_C%row_group .and. blocks%col_group == blocks_C%col_group) then
            rank = blocks_C%rankmax
            allocate (UU(blocks_C%M_loc, rank))
            UU = 0
            call BF_block_MVP_dat(blocks_D, 'N', blocks_D%M_loc, blocks_D%N_loc, rank, blocks_C%ButterflyU%blocks(1)%matrix, blocks_D%N_loc, UU, blocks_D%M_loc, BPACK_cone, BPACK_czero, ptree, stats)
            UU = -UU - blocks_C%ButterflyU%blocks(1)%matrix
            allocate (VV(blocks_C%N_loc, rank))
            VV = 0
            call BF_block_MVP_dat(blocks_A, 'T', blocks_A%M_loc, blocks_A%N_loc, rank, blocks_C%ButterflyV%blocks(1)%matrix,blocks_A%M_loc, VV, blocks_A%N_loc, BPACK_cone, BPACK_czero, ptree, stats)
            VV = VV + blocks_C%ButterflyV%blocks(1)%matrix

         elseif (blocks%row_group == blocks_D%row_group .and. blocks%col_group == blocks_D%col_group) then
            rank1 = blocks_C%rankmax
            allocate (UU1(blocks_C%M_loc, rank1))
            UU1 = 0
            call BF_block_MVP_dat(blocks_D, 'N', blocks_D%M_loc, blocks_D%N_loc, rank1, blocks_C%ButterflyU%blocks(1)%matrix,blocks_D%N_loc, UU1, blocks_D%M_loc, BPACK_cone, BPACK_czero, ptree, stats)
            UU1 = UU1 + blocks_C%ButterflyU%blocks(1)%matrix
            allocate (TT1(blocks_A%N_loc, rank1))
            TT1 = 0
            call BF_block_MVP_dat(blocks_A, 'T', blocks_A%M_loc, blocks_A%N_loc, rank1, blocks_C%ButterflyV%blocks(1)%matrix, blocks_A%M_loc, TT1, blocks_A%N_loc, BPACK_cone, BPACK_czero, ptree, stats)
            TT1 = TT1 + blocks_C%ButterflyV%blocks(1)%matrix
            allocate (SS1(blocks_B%N_loc, rank1))
            SS1 = 0
            call BF_block_MVP_dat(blocks_B, 'T', blocks_B%M_loc, blocks_B%N_loc, rank1, TT1, blocks_B%M_loc, SS1, blocks_B%N_loc, BPACK_cone, BPACK_czero, ptree, stats)
            allocate (VV1(blocks_B%N_loc, rank1))
            VV1 = 0
            call BF_block_MVP_dat(blocks_D, 'T', blocks_D%M_loc, blocks_D%N_loc, rank1, SS1, blocks_D%M_loc, VV1, blocks_D%N_loc, BPACK_cone, BPACK_czero, ptree, stats)
            VV1 = VV1 + SS1

            rank2 = blocks_D%rankmax
            allocate (UU2(blocks_D%M_loc, rank2))
            UU2 = blocks_D%ButterflyU%blocks(1)%matrix
            allocate (VV2(blocks_D%N_loc, rank2))
            VV2 = blocks_D%ButterflyV%blocks(1)%matrix
            rank = rank1 + rank2

            allocate (UU(blocks_D%M_loc, rank))
            UU(:, 1:rank1) = UU1
            UU(:, 1 + rank1:rank) = UU2
            allocate (VV(blocks_D%N_loc, rank))
            VV(:, 1:rank1) = VV1
            VV(:, 1 + rank1:rank) = VV2

            deallocate (UU1, UU2, VV1, VV2, TT1, SS1)

         endif

         blocks%rankmax = rank
         call ComputeParallelIndices(blocks, pgno, ptree, msh)

         !!!!>**** generate 2D grid blacs quantities
         ctxt = ptree%pgrp(pgno)%ctxt
         call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            call descinit_wp(descsmatU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            allocate (blocks%ButterflyU%blocks(1)%matrix(max(myArows,1), max(myAcols,1)))
            blocks%ButterflyU%blocks(1)%matrix = 0

            myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            call descinit_wp(descsmatV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            allocate (blocks%ButterflyV%blocks(1)%matrix(max(myArows,1), max(myAcols,1)))
            blocks%ButterflyV%blocks(1)%matrix = 0
         else
            descsmatU(2) = -1
            descsmatV(2) = -1
            allocate (blocks%ButterflyU%blocks(1)%matrix(1, 1))
            allocate (blocks%ButterflyV%blocks(1)%matrix(1, 1))
         endif

         !!!!>**** redistribution of UV factors
         call Redistribute1Dto2D(UU, blocks%M_p, 0, pgno, blocks%ButterflyU%blocks(1)%matrix, blocks%M, 0, pgno, rank, ptree)
         call Redistribute1Dto2D(VV, blocks%N_p, 0, pgno, blocks%ButterflyV%blocks(1)%matrix, blocks%N, 0, pgno, rank, ptree)

         deallocate (UU, VV)

      else
         if (IOwnPgrp(ptree, pgno)) then
            allocate (blocks%sons(2, 1))
            if (mod(cridx + 1, 2) == 0) then  ! split along column dimension
               blocks%sons(1, 1)%row_group = blocks%row_group
               blocks%sons(1, 1)%headm = blocks%headm
               blocks%sons(1, 1)%col_group = blocks%col_group*2
               blocks%sons(1, 1)%headn = blocks%headn
            else  ! split along row dimension
               blocks%sons(1, 1)%col_group = blocks%col_group
               blocks%sons(1, 1)%headn = blocks%headn
               blocks%sons(1, 1)%row_group = blocks%row_group*2
               blocks%sons(1, 1)%headm = blocks%headm
            endif

            blocks%sons(1, 1)%M = msh%basis_group(blocks%sons(1, 1)%row_group)%tail - msh%basis_group(blocks%sons(1, 1)%row_group)%head + 1
            blocks%sons(1, 1)%N = msh%basis_group(blocks%sons(1, 1)%col_group)%tail - msh%basis_group(blocks%sons(1, 1)%col_group)%head + 1
            call LR_BuildABCD(blocks%sons(1, 1), partitioned_block, option, msh, stats, ptree, pgno, cridx + 1)

            if (mod(cridx + 1, 2) == 0) then  ! split along column dimension
               blocks%sons(2, 1)%row_group = blocks%row_group
               blocks%sons(2, 1)%headm = blocks%headm
               blocks%sons(2, 1)%col_group = blocks%col_group*2 + 1
               blocks%sons(2, 1)%headn = blocks%headn + blocks%sons(1, 1)%N
            else  ! split along row dimension
               blocks%sons(2, 1)%col_group = blocks%col_group
               blocks%sons(2, 1)%headn = blocks%headn
               blocks%sons(2, 1)%row_group = blocks%row_group*2 + 1
               blocks%sons(2, 1)%headm = blocks%headm + blocks%sons(1, 1)%M
            endif
            blocks%sons(2, 1)%M = msh%basis_group(blocks%sons(2, 1)%row_group)%tail - msh%basis_group(blocks%sons(2, 1)%row_group)%head + 1
            blocks%sons(2, 1)%N = msh%basis_group(blocks%sons(2, 1)%col_group)%tail - msh%basis_group(blocks%sons(2, 1)%col_group)%head + 1
            call LR_BuildABCD(blocks%sons(2, 1), partitioned_block, option, msh, stats, ptree, pgno, cridx + 1)
         endif

      endif

   end subroutine LR_BuildABCD

!>**** Merge four child butterflies (of the same level) into a bigger one
   !partitioned_block: (input) partitioned_block%sons store the four smaller BF
   !ptree: (input) process tree
   !stats: (input) statistics
   !option: (input) containing compression options
   !msh: (input) containing mesh
   !blocks_o: (inout) the parent BF
   subroutine BF_Aggregate(partitioned_block, blocks_o, ptree, stats, option, msh)


      implicit none
      integer level_p, ADflag, iii, jjj
      integer mm1, mm2, nn1, nn2, M1, M2, N1, N2, ii, jj, ss, kk, j, i, mm, nn
      integer level_butterfly, num_blocks, level_butterfly_c, level_butterfly_o, level_butterfly_dummy, level_final, num_blocks_c, level, num_col, num_row, num_rowson, num_colson

      type(matrixblock)::partitioned_block
      type(matrixblock)::blocks_o, blocks_dummyL(2), blocks_dummyR(2)
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D, blocks
      DT, allocatable:: matrixtemp1(:, :), matrixtemp2(:, :), vin(:, :), vout1(:, :), vout2(:, :)
      DT::ctemp1, ctemp2
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh
      integer pgno_i, pgno_o, pgno
      DTR, allocatable:: Singular(:)
      DT, allocatable::UU(:, :), VV(:, :)
      integer idx_r, inc_r, nr, idx_c, inc_c, nc, idx, inc, index_i, index_j, mn, nblk_loc, num_blk, rank, row_group, col_group, level_o
      integer ierr
      character mode

      blocks_A => partitioned_block%sons(1, 1)
      blocks_B => partitioned_block%sons(1, 2)
      blocks_C => partitioned_block%sons(2, 1)
      blocks_D => partitioned_block%sons(2, 2)

      pgno_i = blocks_A%pgno
      pgno_o = blocks_o%pgno
      pgno = min(pgno_i, pgno_o)

      if (IOwnPgrp(ptree, pgno_o)) then
         level_butterfly_o = blocks_o%level_butterfly
      else
         level_butterfly_o = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_c = blocks_A%level_butterfly
      else
         level_butterfly_c = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_c, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      call assert(level_butterfly_o == 2 + level_butterfly_c, 'BF_ABCD only supports merging L-2-level BFs into a L-level BF')

      call BF_delete(blocks_o, 0)
      if (.not. allocated(blocks_o%ButterflyKerl)) then
         allocate (blocks_o%ButterflyKerl(level_butterfly_o))
      endif
      blocks_o%pgno = pgno_i
      blocks_o%level_half = floor_safe(dble(level_butterfly_o)/2d0) ! from outer to inner

      num_blocks = 2**level_butterfly_o
      !>****** row-wise ordering from right side
      do level = 0, blocks_o%level_half

         call GetLocalBlockRange(ptree, blocks_o%pgno, level, level_butterfly_o, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
         if (level == 0) then
            blocks_o%ButterflyV%idx = idx_c
            blocks_o%ButterflyV%inc = inc_c
            blocks_o%ButterflyV%nblk_loc = nc
            blocks_o%ButterflyV%num_blk = num_blocks
         elseif (level == level_butterfly_o + 1) then
            blocks_o%ButterflyU%idx = idx_r
            blocks_o%ButterflyU%inc = inc_r
            blocks_o%ButterflyU%nblk_loc = nr
            blocks_o%ButterflyU%num_blk = num_blocks
         else
            blocks_o%ButterflyKerl(level)%num_row = 2**level
            blocks_o%ButterflyKerl(level)%num_col = 2**(level_butterfly_o - level + 1)
            blocks_o%ButterflyKerl(level)%idx_r = idx_r
            blocks_o%ButterflyKerl(level)%inc_r = inc_r
            blocks_o%ButterflyKerl(level)%nr = nr
            blocks_o%ButterflyKerl(level)%idx_c = idx_c*2 - 1
            blocks_o%ButterflyKerl(level)%inc_c = inc_c
            blocks_o%ButterflyKerl(level)%nc = nc*2
         endif
      enddo

      !>****** column-wise ordering from left side
      level_final = blocks_o%level_half + 1
      do level = level_butterfly_o + 1, level_final, -1

         call GetLocalBlockRange(ptree, blocks_o%pgno, level, level_butterfly_o, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

         if (level == 0) then
            blocks_o%ButterflyV%idx = idx_c
            blocks_o%ButterflyV%inc = inc_c
            blocks_o%ButterflyV%nblk_loc = nc
            blocks_o%ButterflyV%num_blk = num_blocks
         elseif (level == level_butterfly_o + 1) then
            blocks_o%ButterflyU%idx = idx_r
            blocks_o%ButterflyU%inc = inc_r
            blocks_o%ButterflyU%nblk_loc = nr
            blocks_o%ButterflyU%num_blk = num_blocks
         else
            blocks_o%ButterflyKerl(level)%num_row = 2**level
            blocks_o%ButterflyKerl(level)%num_col = 2**(level_butterfly_o - level + 1)
            blocks_o%ButterflyKerl(level)%idx_r = idx_r*2 - 1
            blocks_o%ButterflyKerl(level)%inc_r = inc_r
            blocks_o%ButterflyKerl(level)%nr = nr*2
            blocks_o%ButterflyKerl(level)%idx_c = idx_c
            blocks_o%ButterflyKerl(level)%inc_c = inc_c
            blocks_o%ButterflyKerl(level)%nc = nc
         endif
      enddo

      !!!! level level_butterfly_o+1 and level_butterfly_o only requires flop operations and communication
      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_dummy = max(level_butterfly_o - 1, 0)
         write (*, *) 1, 's', ptree%MyID, ptree%pgrp(pgno_i)%nproc, ptree%pgrp(pgno_o)%nproc
         do iii = 1, 2
            blocks_dummyL(iii)%level_butterfly = level_butterfly_dummy
            blocks_dummyL(iii)%level_half = level_butterfly_dummy - 1 ! the value of level_half is only used to guarantee column-wise all2all
            mode = 'C'
            call GetLocalBlockRange(ptree, pgno_i, level_butterfly_dummy + 1, level_butterfly_dummy, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)
            num_blk = 2**level_butterfly_dummy
            idx = idx_r
            inc = inc_r
            nblk_loc = nr
            if (nblk_loc > 0) then
               ! blocks_dummyL(iii)%ButterflyU%idx=idx+(iii-1)*num_blk/2
               blocks_dummyL(iii)%ButterflyU%idx = idx
               blocks_dummyL(iii)%ButterflyU%inc = inc
               blocks_dummyL(iii)%ButterflyU%nblk_loc = nblk_loc
               blocks_dummyL(iii)%ButterflyU%num_blk = num_blk
               blocks_dummyL(iii)%pgno = pgno_i
               allocate (blocks_dummyL(iii)%ButterflyU%blocks(blocks_dummyL(iii)%ButterflyU%nblk_loc))
            endif

            allocate (blocks_dummyL(iii)%ButterflyKerl(1:level_butterfly_dummy)) ! only need one kernel level
            ! level_o = level_butterfly_o
            ! num_row=2**level_o
            ! num_col=2**(level_butterfly_o-level_o+1)
            num_row = 2**level_butterfly_dummy
            num_col = 2

            call GetLocalBlockRange(ptree, pgno_i, level_butterfly_dummy, level_butterfly_dummy, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)
            if (mode == 'R') then
               idx_c = idx_c*2 - 1
               inc_c = inc_c
               nc = nc*2
            elseif (mode == 'C') then
               idx_r = idx_r*2 - 1
               inc_r = inc_r
               nr = nr*2
            endif

            if (nr > 0 .and. nc > 0) then
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_r = idx_r
               ! blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_r=idx_r+(iii-1)*num_row/2
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_c = idx_c
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%inc_r = inc_r
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%inc_c = inc_c
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nr = nr
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nc = nc
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%num_row = num_row
               blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%num_col = num_col
               allocate (blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nr, blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nc))
            endif

            if (nblk_loc > 0 .and. nr > 0 .and. nc > 0) then
               do ii = 1, partitioned_block%sons(iii, 1)%ButterflyU%nblk_loc
                  nn1 = size(partitioned_block%sons(iii, 1)%ButterflyU%blocks(ii)%matrix, 2)
                  nn2 = size(partitioned_block%sons(iii, 2)%ButterflyU%blocks(ii)%matrix, 2)
                  mm = size(partitioned_block%sons(iii, 1)%ButterflyU%blocks(ii)%matrix, 1)
                  allocate (matrixtemp1(mm, nn1 + nn2))
                  index_i = (ii - 1)*blocks_dummyL(iii)%ButterflyU%inc + blocks_dummyL(iii)%ButterflyU%idx
                  row_group = blocks_o%row_group*2**level_butterfly_o + (index_i*2 - 1) - 1
                  mm1 = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
                  mm2 = mm - mm1
                  matrixtemp1(:, 1:nn1) = partitioned_block%sons(iii, 1)%ButterflyU%blocks(ii)%matrix
                  matrixtemp1(:, 1 + nn1:nn1 + nn2) = partitioned_block%sons(iii, 2)%ButterflyU%blocks(ii)%matrix

                  !!!! low rank for the first half rows
                  allocate (matrixtemp2(mm1, nn1 + nn2))
                  matrixtemp2 = matrixtemp1(1:mm1, :)
                  mn = min(mm1, nn1 + nn2)
                  allocate (UU(mm1, mn), VV(mn, nn1 + nn2), Singular(mn))
                  call SVD_Truncate(matrixtemp2, mm1, nn1 + nn2, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
                  do ss = 1, rank
                     UU(:, ss) = UU(:, ss)*Singular(ss)
                  enddo
                  allocate (blocks_dummyL(iii)%ButterflyU%blocks(ii*2 - 1)%matrix(mm1, rank))
                  blocks_dummyL(iii)%ButterflyU%blocks(ii*2 - 1)%matrix = UU(:, 1:rank)

                  allocate (blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2 - 1, 1)%matrix(rank, nn1))
                  blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2 - 1, 1)%matrix = VV(1:rank, 1:nn1)
                  allocate (blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2 - 1, 2)%matrix(rank, nn2))
                  blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2 - 1, 2)%matrix = VV(1:rank, 1 + nn1:nn1 + nn2)
                  deallocate (UU, VV, Singular, matrixtemp2)

                  !!!! low rank for the second half rows
                  allocate (matrixtemp2(mm2, nn1 + nn2))
                  matrixtemp2 = matrixtemp1(1 + mm1:mm, :)
                  mn = min(mm2, nn1 + nn2)
                  allocate (UU(mm2, mn), VV(mn, nn1 + nn2), Singular(mn))
                  call SVD_Truncate(matrixtemp2, mm2, nn1 + nn2, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
                  do ss = 1, rank
                     UU(:, ss) = UU(:, ss)*Singular(ss)
                  enddo
                  allocate (blocks_dummyL(iii)%ButterflyU%blocks(ii*2)%matrix(mm2, rank))
                  blocks_dummyL(iii)%ButterflyU%blocks(ii*2)%matrix = UU(:, 1:rank)

                  allocate (blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2, 1)%matrix(rank, nn1))
                  blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2, 1)%matrix = VV(1:rank, 1:nn1)
                  allocate (blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2, 2)%matrix(rank, nn2))
                  blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2, 2)%matrix = VV(1:rank, 1 + nn1:nn1 + nn2)
                  deallocate (UU, VV, Singular, matrixtemp2)
                  deallocate (matrixtemp1)
               enddo
            endif
            num_blk = 2**level_butterfly_dummy
            call BF_all2all_UV(blocks_dummyL(iii), pgno_i, blocks_dummyL(iii)%ButterflyU, level_butterfly_dummy + 1, (iii - 1)*num_blk, blocks_o, pgno_i, blocks_o%ButterflyU, level_butterfly_o + 1, stats, ptree)
            call BF_all2all_ker(blocks_dummyL(iii), pgno_i, blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy), level_butterfly_dummy, (iii - 1)*num_blk, 0, blocks_o, pgno_i, blocks_o%ButterflyKerl(level_butterfly_o), level_butterfly_o, stats, ptree)

         enddo

         do jjj = 1, 2
            blocks_dummyR(jjj)%level_butterfly = level_butterfly_dummy
            blocks_dummyR(jjj)%level_half = level_butterfly_dummy + 1  ! the value of level_half is only used to guarantee row-wise all2all
            mode = 'R'
            call GetLocalBlockRange(ptree, pgno_i, 0, level_butterfly_dummy, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)
            num_blk = 2**level_butterfly_dummy
            ! num_blk=2**level_butterfly_o
            idx = idx_c
            inc = inc_c
            nblk_loc = nc
            if (nblk_loc > 0) then
               ! blocks_dummyR(jjj)%ButterflyV%idx=idx+(jjj-1)*num_blk/2
               blocks_dummyR(jjj)%ButterflyV%idx = idx
               blocks_dummyR(jjj)%ButterflyV%inc = inc
               blocks_dummyR(jjj)%ButterflyV%nblk_loc = nblk_loc
               blocks_dummyR(jjj)%ButterflyV%num_blk = num_blk
               blocks_dummyR(jjj)%pgno = pgno_i
               allocate (blocks_dummyR(jjj)%ButterflyV%blocks(blocks_dummyR(jjj)%ButterflyV%nblk_loc))
            endif

            allocate (blocks_dummyR(jjj)%ButterflyKerl(1:level_butterfly_dummy)) ! only need one kernel level
            ! level_o = 1
            ! num_row=2**level_o
            ! num_col=2**(level_butterfly_o-level_o+1)
            num_row = 2
            num_col = 2**level_butterfly_dummy
            call GetLocalBlockRange(ptree, pgno_i, 1, level_butterfly_dummy, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)
            if (mode == 'R') then
               idx_c = idx_c*2 - 1
               inc_c = inc_c
               nc = nc*2
            elseif (mode == 'C') then
               idx_r = idx_r*2 - 1
               inc_r = inc_r
               nr = nr*2
            endif

            if (nr > 0 .and. nc > 0) then
               blocks_dummyR(jjj)%ButterflyKerl(1)%idx_r = idx_r
               ! blocks_dummyR(jjj)%ButterflyKerl(1)%idx_c=idx_c+(jjj-1)*num_row/2
               blocks_dummyR(jjj)%ButterflyKerl(1)%idx_c = idx_c
               blocks_dummyR(jjj)%ButterflyKerl(1)%inc_r = inc_r
               blocks_dummyR(jjj)%ButterflyKerl(1)%inc_c = inc_c
               blocks_dummyR(jjj)%ButterflyKerl(1)%nr = nr
               blocks_dummyR(jjj)%ButterflyKerl(1)%nc = nc
               blocks_dummyR(jjj)%ButterflyKerl(1)%num_row = num_row
               blocks_dummyR(jjj)%ButterflyKerl(1)%num_col = num_col
               allocate (blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(blocks_dummyR(jjj)%ButterflyKerl(1)%nr, blocks_dummyR(jjj)%ButterflyKerl(1)%nc))
            endif

            if (nblk_loc > 0 .and. nr > 0 .and. nc > 0) then
               do ii = 1, partitioned_block%sons(1, jjj)%ButterflyV%nblk_loc
                  mm1 = size(partitioned_block%sons(1, jjj)%ButterflyV%blocks(ii)%matrix, 2)
                  mm2 = size(partitioned_block%sons(2, jjj)%ButterflyV%blocks(ii)%matrix, 2)
                  nn = size(partitioned_block%sons(1, jjj)%ButterflyV%blocks(ii)%matrix, 1)
                  allocate (matrixtemp1(mm1 + mm2, nn))

                  index_j = (ii - 1)*blocks_dummyR(jjj)%ButterflyV%inc + blocks_dummyR(jjj)%ButterflyV%idx
                  col_group = blocks_o%col_group*2**level_butterfly_o + (index_j*2 - 1) - 1
                  nn1 = msh%basis_group(col_group)%tail - msh%basis_group(col_group)%head + 1
                  nn2 = nn - nn1
                  call copymatT(partitioned_block%sons(1, jjj)%ButterflyV%blocks(ii)%matrix, matrixtemp1(1:mm1, :), nn, mm1)
                  call copymatT(partitioned_block%sons(2, jjj)%ButterflyV%blocks(ii)%matrix, matrixtemp1(1 + mm1:mm1 + mm2, :), nn, mm2)

                  !!!! low rank for the first half rows
                  allocate (matrixtemp2(mm1 + mm2, nn1))
                  matrixtemp2 = matrixtemp1(:, 1:nn1)
                  mn = min(mm1 + mm2, nn1)
                  allocate (UU(mm1 + mm2, mn), VV(mn, nn1), Singular(mn))
                  call SVD_Truncate(matrixtemp2, mm1 + mm2, nn1, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
                  do ss = 1, rank
                     VV(ss, :) = VV(ss, :)*Singular(ss)
                  enddo
                  allocate (blocks_dummyR(jjj)%ButterflyV%blocks(ii*2 - 1)%matrix(nn1, rank))
                  call copymatT(VV(1:rank, :), blocks_dummyR(jjj)%ButterflyV%blocks(ii*2 - 1)%matrix, rank, nn1)
                  allocate (blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1, ii*2 - 1)%matrix(mm1, rank))
                  blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1, ii*2 - 1)%matrix = UU(1:mm1, 1:rank)
                  allocate (blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2, ii*2 - 1)%matrix(mm2, rank))
                  blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2, ii*2 - 1)%matrix = UU(1 + mm1:mm1 + mm2, 1:rank)
                  deallocate (UU, VV, Singular, matrixtemp2)

                  !!!! low rank for the second half rows
                  allocate (matrixtemp2(mm1 + mm2, nn2))
                  matrixtemp2 = matrixtemp1(:, 1 + nn1:nn1 + nn2)
                  mn = min(mm1 + mm2, nn2)
                  allocate (UU(mm1 + mm2, mn), VV(mn, nn2), Singular(mn))
                  call SVD_Truncate(matrixtemp2, mm1 + mm2, nn2, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
                  do ss = 1, rank
                     VV(ss, :) = VV(ss, :)*Singular(ss)
                  enddo
                  allocate (blocks_dummyR(jjj)%ButterflyV%blocks(ii*2)%matrix(nn2, rank))
                  call copymatT(VV(1:rank, :), blocks_dummyR(jjj)%ButterflyV%blocks(ii*2)%matrix, rank, nn2)
                  allocate (blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1, ii*2)%matrix(mm1, rank))
                  blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1, ii*2)%matrix = UU(1:mm1, 1:rank)
                  allocate (blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2, ii*2)%matrix(mm2, rank))
                  blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2, ii*2)%matrix = UU(1 + mm1:mm1 + mm2, 1:rank)
                  deallocate (UU, VV, Singular, matrixtemp2)
                  deallocate (matrixtemp1)
               enddo
            endif
            num_blk = 2**level_butterfly_dummy
            call BF_all2all_UV(blocks_dummyR(jjj), pgno_i, blocks_dummyR(jjj)%ButterflyV, 0, (jjj - 1)*num_blk, blocks_o, pgno_i, blocks_o%ButterflyV, 0, stats, ptree)
            call BF_all2all_ker(blocks_dummyR(jjj), pgno_i, blocks_dummyR(jjj)%ButterflyKerl(1), 1, 0, (jjj - 1)*num_blk, blocks_o, pgno_i, blocks_o%ButterflyKerl(1), 1, stats, ptree)
         enddo

         do iii = 1, 2
         do jjj = 1, 2
         do level = 2, blocks_o%level_butterfly - 1
            num_row = 2**level/2
            num_col = 2**(level_butterfly_o - level + 1)/2
            call BF_all2all_ker(partitioned_block%sons(iii, jjj), pgno_i, partitioned_block%sons(iii, jjj)%ButterflyKerl(level - 1), level - 1, (iii - 1)*num_row, (jjj - 1)*num_col, blocks_o, pgno_i, blocks_o%ButterflyKerl(level), level, stats, ptree)
         enddo
         enddo
         enddo

         call BF_delete(blocks_dummyL(1), 1)
         call BF_delete(blocks_dummyL(2), 1)
         call BF_delete(blocks_dummyR(1), 1)
         call BF_delete(blocks_dummyR(2), 1)
      endif

      !!!! redistribute the parent BF to its desired process group
      if (IOwnPgrp(ptree, pgno)) then
         call BF_ReDistribute_Inplace(blocks_o, pgno_o, stats, ptree, msh)
      endif

      if (IOwnPgrp(ptree, pgno_o)) then
         call BF_get_rank(blocks_o, ptree)
      endif

      ! mm1=msh%basis_group(blocks_A%row_group)%tail-msh%basis_group(blocks_A%row_group)%head+1
      ! nn1=msh%basis_group(blocks_A%col_group)%tail-msh%basis_group(blocks_A%col_group)%head+1
      ! mm2=msh%basis_group(blocks_D%row_group)%tail-msh%basis_group(blocks_D%row_group)%head+1
      ! nn2=msh%basis_group(blocks_D%col_group)%tail-msh%basis_group(blocks_D%col_group)%head+1

      ! allocate(vin(nn1+nn2,1))
      ! vin = 1
      ! allocate(vout1(mm1+mm2,1))
      ! vout1 = 0
      ! allocate(vout2(mm1+mm2,1))
      ! vout2 = 0

      ! ctemp1 = 1d0; ctemp2 = 0d0
      ! call BF_block_MVP_dat(blocks_i,'N',mm1+mm2,nn1+nn2,1,vin,vout1,ctemp1,ctemp2)

      ! ctemp1 = 1d0; ctemp2 = 0d0
      ! call BF_block_MVP_dat(blocks_A,'N',mm1,nn1,1,vin(1:nn1,:),vout2(1:mm1,:),ctemp1,ctemp2)
      ! ctemp1 = 1d0; ctemp2 = 1d0
      ! call BF_block_MVP_dat(blocks_B,'N',mm1,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1:mm1,:),ctemp1,ctemp2)

      ! ctemp1 = 1d0; ctemp2 = 0d0
      ! call BF_block_MVP_dat(blocks_C,'N',mm2,nn1,1,vin(1:nn1,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
      ! ctemp1 = 1d0; ctemp2 = 1d0
      ! call BF_block_MVP_dat(blocks_D,'N',mm2,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)

      ! write(*,*)'spliting error:',fnorm(vout1-vout2,mm1+mm2,1)/fnorm(vout1,mm1+mm2,1)
      ! deallocate(vin,vout1,vout2)

   end subroutine BF_Aggregate

   subroutine BF_split(blocks_i, blocks_o, ptree, stats, msh, option, splitpg)



#ifdef HAVE_OPENMP
     use omp_lib
#endif

      implicit none
      integer level_p, ADflag, iii, jjj
      integer mm1, mm2, nn1, nn2, M1, M2, N1, N2, ii, jj, kk, j, i, mm, nn, pgno_sub_mine, gg
      integer level_butterfly, num_blocks, level_butterfly_c, num_blocks_c, level, num_col, num_row, num_rowson, num_colson

      type(matrixblock), target::blocks_i
      type(matrixblock)::blocks_o, blocks_dummy, blocks_i_copy
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D, blocks
      DT, allocatable:: matrixtemp1(:, :), matrixtemp2(:, :), vin(:, :), vout1(:, :), vout2(:, :)
      DT::ctemp1, ctemp2
      type(mesh)::msh
      type(proctree)::ptree
      type(Hstat)::stats
      integer pgno, pgnos(2,2), Maxlevel, ld
      integer, allocatable::M_p_sub(:, :), M_p_sub1(:, :), N_p_sub(:, :), N_p_sub1(:, :)
      type(matrixblock), pointer::block_c_o
      integer rank, ierr, Maxgrp
      type(Hoption)::option
      real(kind=8)::n3,n4
      integer,optional:: splitpg

      blocks_A => blocks_o%sons(1, 1)
      blocks_B => blocks_o%sons(1, 2)
      blocks_C => blocks_o%sons(2, 1)
      blocks_D => blocks_o%sons(2, 2)

      call BF_copy('N', blocks_i, blocks_i_copy)

      if (blocks_i_copy%level_butterfly == 0) then
         level_butterfly = 0
      else
         level_butterfly = max(blocks_i_copy%level_butterfly - 2, 0)
      endif

      Maxlevel = GetTreelevel(msh%Maxgroup) - 1

      do iii = 1, 2
      do jjj = 1, 2

         if(present(splitpg))then
            !>*** split process groups row-wise
            Maxgrp = 2**(ptree%nlevel) - 1
            pgno = blocks_i_copy%pgno
            if (pgno*2 <= Maxgrp) then
               pgno = pgno*2 + iii-1
            endif
         else
            !>*** try to use the same process group as blocks_i
            pgno = blocks_i_copy%pgno
            do while (Maxlevel - blocks_i%level - 1 < ptree%nlevel - GetTreelevel(pgno))
               ! do while(level_butterfly<ptree%nlevel-GetTreelevel(pgno))
               pgno = pgno*2
            enddo
         endif
         pgnos(iii,jjj)=pgno

         blocks => blocks_o%sons(iii, jjj)
         blocks%level = blocks_i_copy%level + 1
         blocks%row_group = blocks_i_copy%row_group*2 + iii - 1
         blocks%col_group = blocks_i_copy%col_group*2 + jjj - 1
         blocks%level_butterfly = max(blocks_i_copy%level_butterfly - 2, 0)
         blocks%level_half = floor_safe(dble(blocks%level_butterfly)/2d0) ! from outer to inner
         blocks%style = 2
         blocks%headm = msh%basis_group(blocks%row_group)%head
         blocks%M = msh%basis_group(blocks%row_group)%tail - msh%basis_group(blocks%row_group)%head + 1
         blocks%headn = msh%basis_group(blocks%col_group)%head
         blocks%N = msh%basis_group(blocks%col_group)%tail - msh%basis_group(blocks%col_group)%head + 1
         blocks%pgno = pgno
         call ComputeParallelIndices(blocks, pgno, ptree, msh)
      enddo
      enddo

      if (blocks_i_copy%level_butterfly == 0) then

         kk = size(blocks_i_copy%ButterflyU%blocks(1)%matrix, 2)
         do iii = 1, 2
         do jjj = 1, 2
            blocks => blocks_o%sons(iii, jjj)
            blocks%level_butterfly = 0
            blocks%level_half = 0
            allocate (blocks%ButterflyU%blocks(1))
            allocate (blocks%ButterflyV%blocks(1))
            if (IOwnPgrp(ptree, blocks%pgno)) then
               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc, kk))
               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc, kk))
               blocks%rankmax = kk
               blocks%rankmin = kk
               blocks%ButterflyU%nblk_loc = 1
               blocks%ButterflyU%inc = 1
               blocks%ButterflyU%idx = 1
               blocks%ButterflyV%nblk_loc = 1
               blocks%ButterflyV%inc = 1
               blocks%ButterflyV%idx = 1
            endif
            n3 = MPI_Wtime()
            call Redistribute1Dto1D(blocks_i_copy%ButterflyU%blocks(1)%matrix, blocks_i_copy%M_loc, blocks_i_copy%M_p, blocks_i_copy%headm, blocks_i_copy%pgno, blocks%ButterflyU%blocks(1)%matrix, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, kk, ptree)
            call Redistribute1Dto1D(blocks_i_copy%ButterflyV%blocks(1)%matrix, blocks_i_copy%N_loc, blocks_i_copy%N_p, blocks_i_copy%headn, blocks_i_copy%pgno, blocks%ButterflyV%blocks(1)%matrix, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, kk, ptree)
            n4 = MPI_Wtime()
            stats%Time_RedistB = stats%Time_RedistB + n4-n3
         enddo
         enddo

      else

         if (IOwnPgrp(ptree, blocks_i_copy%pgno)) then
            call GetPgno_Sub(ptree, blocks_i_copy%pgno, blocks_i_copy%level_butterfly, pgno_sub_mine)
            if (ptree%pgrp(pgno_sub_mine)%nproc > 1) then
               allocate (M_p_sub(ptree%pgrp(pgno_sub_mine)%nproc, 2))
               gg = blocks_i_copy%row_group*2**blocks_i_copy%level_butterfly + blocks_i_copy%ButterflyU%idx - 1
               call ComputeParallelIndicesSub(gg, pgno_sub_mine, ptree, msh, M_p_sub)
               allocate (M_p_sub1(ptree%pgrp(pgno_sub_mine)%nproc, 2))
               M_p_sub1 = -1
               M_p_sub1(1, 1) = 1
               M_p_sub1(1, 2) = msh%basis_group(gg)%tail - msh%basis_group(gg)%head + 1
               rank = size(blocks_i_copy%ButterflyU%blocks(1)%matrix, 2)
               allocate (matrixtemp1(blocks_i_copy%M_loc, rank))
               matrixtemp1 = blocks_i_copy%ButterflyU%blocks(1)%matrix
               deallocate (blocks_i_copy%ButterflyU%blocks(1)%matrix)
               allocate (matrixtemp2(M_p_sub1(1, 2), rank))
               n3 = MPI_Wtime()
               call Redistribute1Dto1D(matrixtemp1, blocks_i_copy%M_loc, M_p_sub, 0, pgno_sub_mine, matrixtemp2, M_p_sub1(1, 2), M_p_sub1, 0, pgno_sub_mine, rank, ptree)
               n4 = MPI_Wtime()
               stats%Time_RedistB = stats%Time_RedistB + n4-n3
               if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                  allocate (blocks_i_copy%ButterflyU%blocks(1)%matrix(M_p_sub1(1, 2), rank))
                  blocks_i_copy%ButterflyU%blocks(1)%matrix = matrixtemp2
               endif
               deallocate (matrixtemp1)
               deallocate (matrixtemp2)
               deallocate (M_p_sub)
               deallocate (M_p_sub1)

               allocate (N_p_sub(ptree%pgrp(pgno_sub_mine)%nproc, 2))
               gg = blocks_i_copy%col_group*2**blocks_i_copy%level_butterfly + blocks_i_copy%ButterflyV%idx - 1
               call ComputeParallelIndicesSub(gg, pgno_sub_mine, ptree, msh, N_p_sub)
               allocate (N_p_sub1(ptree%pgrp(pgno_sub_mine)%nproc, 2))
               N_p_sub1 = -1
               N_p_sub1(1, 1) = 1
               N_p_sub1(1, 2) = msh%basis_group(gg)%tail - msh%basis_group(gg)%head + 1
               rank = size(blocks_i_copy%ButterflyV%blocks(1)%matrix, 2)
               allocate (matrixtemp1(blocks_i_copy%N_loc, rank))
               matrixtemp1 = blocks_i_copy%ButterflyV%blocks(1)%matrix
               deallocate (blocks_i_copy%ButterflyV%blocks(1)%matrix)
               allocate (matrixtemp2(N_p_sub1(1, 2), rank))
               n3 = MPI_Wtime()
               call Redistribute1Dto1D(matrixtemp1, blocks_i_copy%N_loc, N_p_sub, 0, pgno_sub_mine, matrixtemp2, N_p_sub1(1, 2), N_p_sub1, 0, pgno_sub_mine, rank, ptree)
               n4 = MPI_Wtime()
               stats%Time_RedistB = stats%Time_RedistB + n4-n3
               if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                  allocate (blocks_i_copy%ButterflyV%blocks(1)%matrix(N_p_sub1(1, 2), rank))
                  blocks_i_copy%ButterflyV%blocks(1)%matrix = matrixtemp2
               endif
               deallocate (matrixtemp1)
               deallocate (matrixtemp2)
               deallocate (N_p_sub)
               deallocate (N_p_sub1)

            endif
         endif

         !>**** first redistribute blocks_i into blocks_dummy%sons of the same butterfly levels
         blocks_dummy%level_butterfly = blocks_i_copy%level_butterfly
         blocks_dummy%level_half = blocks_i_copy%level_half
         allocate (blocks_dummy%sons(2, 2))
         do ii = 1, 2
         do jj = 1, 2
            allocate (blocks_dummy%sons(ii, jj)%ButterflyKerl(blocks_dummy%level_butterfly))
            blocks_dummy%sons(ii, jj)%level_butterfly = blocks_dummy%level_butterfly
            blocks_dummy%sons(ii, jj)%level_half = blocks_dummy%level_half
         enddo
         enddo

         do level = 0, blocks_dummy%level_butterfly + 1
            if (level == 0) then
               call BF_all2all_V_split(blocks_i_copy, blocks_i_copy%pgno, level, blocks_dummy, pgnos, level, stats, ptree)
            elseif (level == blocks_i_copy%level_butterfly + 1) then
               call BF_all2all_U_split(blocks_i_copy, blocks_i_copy%pgno, level, blocks_dummy, pgnos, level, stats, ptree)
            else
               call BF_all2all_ker_split(blocks_i_copy, blocks_i_copy%pgno, level, blocks_dummy, pgnos, level, stats, ptree)
            endif
         enddo
         !>**** next convert blocks_dummy%sons into  blocks_o%sons
         call BF_convert_to_smallBF(blocks_dummy, blocks_o, stats, ptree)

         do iii = 1, 2
            do jjj = 1, 2
               pgno = pgnos(iii,jjj)
               if (IOwnPgrp(ptree, pgno)) then
                  level_butterfly = max(blocks_i_copy%level_butterfly - 2, 0)
                  call GetPgno_Sub(ptree, pgno, level_butterfly, pgno_sub_mine)
                  if (ptree%pgrp(pgno_sub_mine)%nproc > 1) then

                     block_c_o => blocks_o%sons(iii, jjj)
                     block_c_o%level_butterfly = level_butterfly
                     block_c_o%level_half = floor_safe(dble(block_c_o%level_butterfly)/2d0) ! from outer to inner

                     if (ptree%pgrp(pgno_sub_mine)%head /= ptree%MyID) then
                        if (block_c_o%level_butterfly > 0) then
                           allocate (block_c_o%ButterflyKerl(block_c_o%level_butterfly))
                        endif
                     endif
                     do level = 0, block_c_o%level_butterfly + 1
                        if (level == 0) then
                           block_c_o%ButterflyV%num_blk = 1
                           block_c_o%ButterflyV%nblk_loc = 1
                           block_c_o%ButterflyV%inc = 1
                           if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                              rank = size(block_c_o%ButterflyV%blocks(1)%matrix, 2)
                           endif
                           call MPI_Bcast(block_c_o%ButterflyV%idx, 1, MPI_integer, Main_ID, ptree%pgrp(pgno_sub_mine)%Comm, ierr)
                           call MPI_Bcast(rank, 1, MPI_integer, Main_ID, ptree%pgrp(pgno_sub_mine)%Comm, ierr)

                           allocate (N_p_sub(ptree%pgrp(pgno_sub_mine)%nproc, 2))
                           gg = block_c_o%col_group*2**block_c_o%level_butterfly + block_c_o%ButterflyV%idx - 1
                           call ComputeParallelIndicesSub(gg, pgno_sub_mine, ptree, msh, N_p_sub)
                           allocate (N_p_sub1(ptree%pgrp(pgno_sub_mine)%nproc, 2))
                           N_p_sub1 = -1
                           N_p_sub1(1, 1) = 1
                           N_p_sub1(1, 2) = msh%basis_group(gg)%tail - msh%basis_group(gg)%head + 1

                           if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                              allocate (matrixtemp1(N_p_sub1(1, 2), rank))
                              ld = N_p_sub1(1, 2)
                              matrixtemp1 = block_c_o%ButterflyV%blocks(1)%matrix
                              deallocate (block_c_o%ButterflyV%blocks(1)%matrix)
                           else
                              allocate (matrixtemp1(1, 1))
                              ld = 1
                              matrixtemp1 = 0
                              allocate (block_c_o%ButterflyV%blocks(1))
                           endif
                           allocate (block_c_o%ButterflyV%blocks(1)%matrix(block_c_o%N_loc, rank))
                           n3 = MPI_Wtime()
                           call Redistribute1Dto1D(matrixtemp1, ld, N_p_sub1, 0, pgno_sub_mine, block_c_o%ButterflyV%blocks(1)%matrix, block_c_o%N_loc, N_p_sub, 0, pgno_sub_mine, rank, ptree)
                           n4 = MPI_Wtime()
                           stats%Time_RedistB = stats%Time_RedistB + n4-n3
                           deallocate (N_p_sub)
                           deallocate (N_p_sub1)
                           deallocate (matrixtemp1)
                        elseif (level == block_c_o%level_butterfly + 1) then

                           block_c_o%ButterflyU%num_blk = 1
                           block_c_o%ButterflyU%nblk_loc = 1
                           block_c_o%ButterflyU%inc = 1
                           if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                              rank = size(block_c_o%ButterflyU%blocks(1)%matrix, 2)
                           endif
                           call MPI_Bcast(block_c_o%ButterflyU%idx, 1, MPI_integer, Main_ID, ptree%pgrp(pgno_sub_mine)%Comm, ierr)
                           call MPI_Bcast(rank, 1, MPI_integer, Main_ID, ptree%pgrp(pgno_sub_mine)%Comm, ierr)

                           allocate (M_p_sub(ptree%pgrp(pgno_sub_mine)%nproc, 2))
                           gg = block_c_o%row_group*2**block_c_o%level_butterfly + block_c_o%ButterflyU%idx - 1
                           call ComputeParallelIndicesSub(gg, pgno_sub_mine, ptree, msh, M_p_sub)
                           allocate (M_p_sub1(ptree%pgrp(pgno_sub_mine)%nproc, 2))
                           M_p_sub1 = -1
                           M_p_sub1(1, 1) = 1
                           M_p_sub1(1, 2) = msh%basis_group(gg)%tail - msh%basis_group(gg)%head + 1

                           if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
                              allocate (matrixtemp1(M_p_sub1(1, 2), rank))
                              ld=M_p_sub1(1, 2)
                              matrixtemp1 = block_c_o%ButterflyU%blocks(1)%matrix
                              deallocate (block_c_o%ButterflyU%blocks(1)%matrix)
                           else
                              allocate (matrixtemp1(1, 1))
                              ld=1
                              matrixtemp1 = 0
                              allocate (block_c_o%ButterflyU%blocks(1))
                           endif
                           allocate (block_c_o%ButterflyU%blocks(1)%matrix(block_c_o%M_loc, rank))
                           n3 = MPI_Wtime()
                           call Redistribute1Dto1D(matrixtemp1, ld, M_p_sub1, 0, pgno_sub_mine, block_c_o%ButterflyU%blocks(1)%matrix, block_c_o%M_loc, M_p_sub, 0, pgno_sub_mine, rank, ptree)
                           n4 = MPI_Wtime()
                           stats%Time_RedistB = stats%Time_RedistB + n4-n3
                           deallocate (M_p_sub)
                           deallocate (M_p_sub1)
                           deallocate (matrixtemp1)

                        else
                           if (ptree%pgrp(pgno_sub_mine)%head /= ptree%MyID) then
                              block_c_o%ButterflyKerl(level)%nr = 0
                              block_c_o%ButterflyKerl(level)%inc_r = 0
                              block_c_o%ButterflyKerl(level)%idx_r = 0
                              block_c_o%ButterflyKerl(level)%nc = 0
                              block_c_o%ButterflyKerl(level)%inc_c = 0
                              block_c_o%ButterflyKerl(level)%idx_c = 0
                           endif
                        endif
                     enddo
                  endif
               endif
            enddo
         enddo


         do ii = 1, 2
         do jj = 1, 2
            call BF_delete(blocks_dummy%sons(ii, jj), 1)
         enddo
         enddo
         deallocate (blocks_dummy%sons)
      endif

      do ii = 1, 2
      do jj = 1, 2
         call BF_get_rank(blocks_o%sons(ii, jj), ptree)
      enddo
      enddo

      call BF_delete(blocks_i_copy, 1)
      if(.not. present(splitpg))then
      call BF_split_checkerror(blocks_i, blocks_o, ptree, stats, option)
      endif
   end subroutine BF_split

! Compare a block with its children
   subroutine BF_split_checkerror(blocks_i, blocks_o, ptree, stats, option)


      implicit none
      integer level, ii, M, N, M_loc, N_loc, num_vect_sub, mv, nv
      character trans
      DT, allocatable :: Vin(:, :), Vout(:, :)
      DT, allocatable :: Vin_tmp(:, :), Vout_tmp(:, :), Vbuff(:, :), V1(:, :), V2(:, :), Vout1(:, :), Vout2(:, :), Vo1(:, :), Vo2(:, :)
      DT :: ctemp1, ctemp2, a, b
      type(matrixblock), pointer::blocks_A, blocks_B, blocks_C, blocks_D
      integer groupn, groupm, mm, nn
      type(matrixblock)::blocks_i, blocks_o
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      real(kind=8)::error,n1,n2

      if (IOwnPgrp(ptree, blocks_i%pgno)) then

         blocks_A => blocks_o%sons(1, 1)
         blocks_B => blocks_o%sons(1, 2)
         blocks_C => blocks_o%sons(2, 1)
         blocks_D => blocks_o%sons(2, 2)

         num_vect_sub = 1
         allocate (vin(blocks_i%N_loc, num_vect_sub))
         vin = 1
         allocate (vout1(blocks_i%M_loc, num_vect_sub))
         vout1 = 0
         allocate (vout2(blocks_i%M_loc, num_vect_sub))
         vout2 = 0
         call BF_block_MVP_dat(blocks_i, 'N', blocks_i%M_loc, blocks_i%N_loc, num_vect_sub, vin, blocks_i%N_loc, vout1, blocks_i%M_loc, BPACK_cone, BPACK_czero, ptree, stats)

         allocate (V1(max(1, blocks_A%N_loc), num_vect_sub))
         allocate (V2(max(1, blocks_B%N_loc), num_vect_sub))

         allocate (Vo1(max(1, blocks_A%M_loc), num_vect_sub))
         Vo1 = 0
         allocate (Vo2(max(1, blocks_C%M_loc), num_vect_sub))
         Vo2 = 0

         n1 = MPI_Wtime()
         ! call Redistribute1Dto1D(vin, blocks_i%N_p, 0, blocks_i%pgno, V1, blocks_A%N_p, 0, blocks_A%pgno, num_vect_sub, ptree)
         ! call Redistribute1Dto1D(vin, blocks_i%N_p, 0, blocks_i%pgno, V2, blocks_B%N_p, blocks_A%N, blocks_B%pgno, num_vect_sub, ptree)

         call Redistribute1Dto1D_OnetoTwo(vin, blocks_i%N_loc, blocks_i%N_p, 0, blocks_i%pgno, V1, max(1, blocks_A%N_loc),blocks_A%N_p, 0, blocks_A%pgno,V2, max(1, blocks_B%N_loc), blocks_B%N_p, blocks_A%N, blocks_B%pgno, num_vect_sub, ptree)
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1




         if (blocks_A%M_loc > 0) then
            call BF_block_MVP_dat(blocks_A, 'N', blocks_A%M_loc, blocks_A%N_loc, num_vect_sub, V1, max(1, blocks_A%N_loc), Vo1, max(1, blocks_A%M_loc), BPACK_cone, BPACK_cone, ptree, stats)
            call BF_block_MVP_dat(blocks_B, 'N', blocks_B%M_loc, blocks_B%N_loc, num_vect_sub, V2, max(1, blocks_B%N_loc), Vo1, max(1, blocks_A%M_loc), BPACK_cone, BPACK_cone, ptree, stats)
            call BF_block_MVP_dat(blocks_C, 'N', blocks_C%M_loc, blocks_C%N_loc, num_vect_sub, V1, max(1, blocks_A%N_loc), Vo2, max(1, blocks_C%M_loc), BPACK_cone, BPACK_cone, ptree, stats)
            call BF_block_MVP_dat(blocks_D, 'N', blocks_D%M_loc, blocks_D%N_loc, num_vect_sub, V2, max(1, blocks_B%N_loc), Vo2, max(1, blocks_C%M_loc), BPACK_cone, BPACK_cone, ptree, stats)
         endif

         n1 = MPI_Wtime()
         ! call Redistribute1Dto1D(Vo1, blocks_A%M_p, 0, blocks_A%pgno, vout2, blocks_i%M_p, 0, blocks_i%pgno, num_vect_sub, ptree)
         ! call Redistribute1Dto1D(Vo2, blocks_C%M_p, blocks_A%M, blocks_C%pgno, vout2, blocks_i%M_p, 0, blocks_i%pgno, num_vect_sub, ptree)

         call Redistribute1Dto1D_TwotoOne(Vo1, max(1, blocks_A%M_loc), blocks_A%M_p, 0, blocks_A%pgno,Vo2, max(1, blocks_C%M_loc), blocks_C%M_p, blocks_A%M, blocks_C%pgno, vout2, blocks_i%M_loc, blocks_i%M_p, 0, blocks_i%pgno, num_vect_sub, ptree)
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1

         if(fnorm(vout1, blocks_i%M_loc, 1)<BPACK_SafeUnderflow .and. fnorm(vout2, blocks_i%M_loc, 1)<BPACK_SafeUnderflow)then
            error = 0d0
         else
            error = fnorm(vout1 - vout2, blocks_i%M_loc, 1)/fnorm(vout1, blocks_i%M_loc, 1)
         endif

         if (ptree%MyID == ptree%pgrp(blocks_i%pgno)%head .and. option%verbosity >= 2) write (*, '(A38,I5,A8,I5,A8,Es14.7,A8,I5)') 'Split L_in:', blocks_i%level_butterfly, 'L_out:', blocks_A%level_butterfly, ' error:', error, ' #nproc:', ptree%pgrp(blocks_A%pgno)%nproc

         deallocate (vin, vout1, vout2)
         deallocate (V1)
         deallocate (V2)
         deallocate (Vo1)
         deallocate (Vo2)
      endif

   end subroutine BF_split_checkerror

   subroutine BF_get_rank_ABCD(partitioned_block, rankmax)


      implicit none

      integer rankmax, ii, jj
      type(matrixblock)::partitioned_block

      rankmax = -1000
      do ii = 1, 2
      do jj = 1, 2
         rankmax = max(rankmax, partitioned_block%sons(ii, jj)%rankmax)
      enddo
      enddo
   end subroutine BF_get_rank_ABCD

!>**** Update one off-diagonal block in HODLR/HODBF compressed as
! Bplus/Butterfly/LR by multiplying on it left the inverse of diagonal block
! If LR, call LR_Sblock; if butterfly, call BF_randomized; if Bplus, call Bplus_randomized_constr
   !ho_bf1: working HODLR/HODBF
   !level_c: level# of the block in HODLR/HODBF
   !rowblock: block# of the block at this level in HODLR/HODBF
   !option: containing compression options
   !stats: statistics
   !ptree: process tree
   subroutine Bplus_Sblock_randomized_memfree(ho_bf1, level_c, rowblock, option, stats, ptree, msh)


#ifdef HAVE_OPENMP
     use omp_lib
#endif
      implicit none

      integer level_c, rowblock
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks, idx_start_diag, idx_end_diag, N_diag, num_vect_sub, level_butterfly_loc, index_i_loc_k, index_i, ii_loc
      integer num_col, num_row, level, mm, nn, ii, jj, tt
      character chara
      real(kind=8) T0
      type(blockplus), pointer::bplus
      type(matrixblock)::block_old
      type(matrixblock), pointer::block_o, blocks
      integer::rank_new_max
      real(kind=8)::rank_new_avr, error, rate, rankrate_inner, rankrate_outter
      integer niter, rank, ntry, rank0, rank0_inner, rank0_outter
      real(kind=8):: error_inout
      real(kind=8):: n1, n2
      type(Hoption)::option
      type(Hstat)::stats
      type(hobf)::ho_bf1
      type(proctree)::ptree
      type(mesh)::msh
      type(matrixblock):: agent_block

      error_inout = 0


      call Bplus_copy(ho_bf1%levels(level_c)%BP(rowblock), ho_bf1%levels(level_c)%BP_inverse_update(rowblock))
      !!!!!!! the forward block BP can be deleted if not used in solution phase

      bplus => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
      if (bplus%Lplus == 1) then

         block_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)
         level_butterfly = block_o%level_butterfly




      if (level_butterfly == 0) then
         call LR_Sblock(ho_bf1, level_c, rowblock, ptree, stats)
      else
#if 0
         ho_bf1%ind_lv = level_c
         ho_bf1%ind_bk = rowblock
         rank0 = block_o%rankmax
         rate = option%rankrate !1.2d0
         call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, ho_bf1, BF_block_MVP_Sblock_dat, error_inout, 'Sblock', option, stats, ptree, msh, operand1=msh)
         stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp

#else




         error_inout=0
         call assert(option%pat_comp/=2,'pat_comp==2 not yet supported in BF_MoveSingular_Ker')
         if(option%pat_comp==3 .and. block_o%level_butterfly>0)then
            call BF_ChangePattern(block_o, option%pat_comp, 1, stats, ptree)
            call BF_MoveSingular_Ker(block_o, 'N', floor_safe(dble(block_o%level_butterfly)/2d0) +1, block_o%level_butterfly, ptree, stats, option%tol_rand)
         endif
         call BF_ChangePattern(block_o, 1, 2, stats, ptree)



         do level = ho_bf1%Maxlevel + 1,level_c+1,-1
            N_diag = 2**(level-level_c-1)
            idx_start_diag = max((rowblock - 1)*N_diag + 1, ho_bf1%levels(level)%Bidxs)
            idx_end_diag = min(rowblock*N_diag, ho_bf1%levels(level)%Bidxe)

            if(level==ho_bf1%Maxlevel + 1)then

               !!$omp parallel do default(shared) private(ii,blocks,index_i,index_i_loc_k,num_vect_sub)
               do ii = idx_start_diag,idx_end_diag
                  blocks => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
                  index_i = ii - ((rowblock - 1)*N_diag + 1) + 1
                  index_i_loc_k = (index_i - block_o%ButterflyU%idx)/block_o%ButterflyU%inc + 1
                  num_vect_sub = size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix,2)
                  call Full_block_MVP_dat(blocks, 'N', blocks%M, num_vect_sub,&
                     &block_o%ButterflyU%blocks(index_i_loc_k)%matrix, size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix,1),block_o%ButterflyU%blocks(index_i_loc_k)%matrix, size(block_o%ButterflyU%blocks(index_i_loc_k)%matrix,1), BPACK_cone, BPACK_czero)
               end do
               !!$omp end parallel do

            else

               do ii = idx_start_diag,idx_end_diag
                  ii_loc = ii - ((rowblock - 1)*N_diag + 1) + 1
                  level_butterfly_loc = ho_bf1%Maxlevel+1-level    !!!! all even level butterfly could be problematic here

                  blocks => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
                  call assert(IOwnPgrp(ptree, blocks%pgno),'I do not own this pgno')
                  n1=MPI_Wtime()
                  call BF_extract_partial(block_o, level_butterfly_loc, ii_loc,blocks%headm,blocks%row_group, 'L', agent_block,blocks%pgno,ptree)
                  rank0 = agent_block%rankmax
                  rate = option%rankrate !1.2d0
                  n2=MPI_Wtime()
                  ! time_tmp = time_tmp + n2-n1

                  ho_bf1%ind_lv = level
                  ho_bf1%ind_bk = ii
                  rank0 = agent_block%rankmax
                  rate = option%rankrate !1.2d0
                  call BF_randomized(agent_block%pgno, level_butterfly_loc, rank0, rate, agent_block, ho_bf1, BF_block_MVP_Sblock_Sml_dat, error, 'Sblock_sml', option, stats, ptree, msh, operand1=msh, vskip=.true.)
                  stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                  error_inout = max(error_inout, error)

                  call BF_ChangePattern(agent_block, 3, 2, stats, ptree)

                  n1=MPI_Wtime()
                  call BF_copyback_partial(block_o, level_butterfly_loc, ii_loc, 'L', agent_block,blocks%pgno,ptree)
                  n2=MPI_Wtime()
                  ! time_tmp = time_tmp + n2-n1

                  call BF_delete(agent_block,1)

               end do


            end if

         end do

#endif

      end if


         if (ptree%MyID == Main_ID .and. option%verbosity >= 1) write (*, '(A10,I5,A6,I3,A8,I3,A11,Es14.7)') 'OneL No. ', rowblock, ' rank:', block_o%rankmax, ' L_butt:', block_o%level_butterfly, ' error:', error_inout

      else

         ho_bf1%ind_lv = level_c
         ho_bf1%ind_bk = rowblock
         Bplus => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
         block_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

         rank0_inner = Bplus%LL(2)%rankmax
         rankrate_inner = 2.0d0

         rank0_outter = block_o%rankmax
         rankrate_outter = option%rankrate !1.2d0
         level_butterfly = block_o%level_butterfly
         call Bplus_randomized_constr(level_butterfly, Bplus, ho_bf1, rank0_inner, rankrate_inner, Bplus_block_MVP_Sblock_dat, rank0_outter, rankrate_outter, Bplus_block_MVP_Outter_Sblock_dat, error_inout, 'Sblock+', option, stats, ptree, msh)

         block_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

         if (option%verbosity >= 1 .and. ptree%myid == ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%head) write (*, '(A10,I5,A6,I3,A8,I3,A11,Es14.7)') 'Mult No. ', rowblock, ' rank:', block_o%rankmax, ' L_butt:', block_o%level_butterfly, ' error:', error_inout

      end if

      return

   end subroutine Bplus_Sblock_randomized_memfree

   subroutine Bplus_inverse_schur_partitionedinverse(ho_bf1, level_c, rowblock, option, stats, ptree, msh)




#ifdef HAVE_OPENMP
     use omp_lib
#endif


      implicit none

      integer level_c, rowblock, ierr
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks, level_butterfly_loc
      integer num_col, num_row, level, mm, nn, ii, jj, tt, ll, llplus, bb, mmb
      character chara
      real(kind=8) T0
      type(matrixblock), pointer::block_o, block_off1, block_off2
      type(matrixblock), pointer::blocks_o_D
      type(matrixblock)::block_tmp
      type(blockplus), pointer::Bplus, Bplus_schur
      integer rank_new_max
      real(kind=8):: rank_new_avr, error, err_avr, err_max
      integer niter
      real(kind=8):: error_inout, rate, rankrate_inner, rankrate_outter
      integer itermax, ntry, cnt, cnt_partial
      real(kind=8):: n1, n2, n3, n4, Memory
      integer rank0, rank0_inner, rank0_outter, Lplus, level_BP, levelm, groupm_start, ij_loc, edge_s, edge_e, edge_first, idx_end_m_ref, idx_start_m_ref, idx_start_b, idx_end_b
      DT, allocatable:: matin(:, :), matout(:, :), matin_tmp(:, :), matout_tmp(:, :)
      DT:: ctemp1, ctemp2
      integer, allocatable :: ipiv(:)
      type(Hoption)::option
      type(Hstat)::stats
      type(hobf)::ho_bf1
      type(matrixblock):: agent_block
      type(blockplus):: agent_bplus
      type(proctree) :: ptree
      type(mesh) :: msh
      real(kind=8) flop

      bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)

      if (bplus%Lplus == 1) then
         call BF_inverse_schur_partitionedinverse(ho_bf1, level_c, rowblock, error_inout, option, stats, ptree, msh)
      else
         ctemp1 = 1d0
         ctemp2 = 0d0

         block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2 - 1)%LL(1)%matrices_block(1)
         block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)
         block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
         ! write(*,*)block_o%row_group,block_o%col_group,level_c,rowblock,block_o%level,'diao'
         block_o%level_butterfly = block_off1%level_butterfly

         Memory = 0

         error_inout = 0

         ho_bf1%ind_lv = level_c
         ho_bf1%ind_bk = rowblock

         rank0_inner = ho_bf1%levels(level_c)%BP_inverse_update(2*rowblock - 1)%LL(2)%rankmax
         rankrate_inner = 2.0d0

         rank0_outter = max(block_off1%rankmax, block_off2%rankmax)
         rankrate_outter = option%rankrate !1.2d0

         level_butterfly = block_o%level_butterfly

         call Bplus_randomized_constr(level_butterfly, Bplus, ho_bf1, rank0_inner, rankrate_inner, Bplus_block_MVP_minusBC_dat, rank0_outter, rankrate_outter, Bplus_block_MVP_Outter_minusBC_dat, error, 'mBC+', option, stats, ptree, msh)
         error_inout = max(error_inout, error)

         ! write(*,*)'good!!!!'
         ! stop
         Bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)
         Lplus = ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%Lplus
         do llplus = Lplus, 1, -1
            do bb = 1, Bplus%LL(llplus)%Nbound
               block_o => Bplus%LL(llplus)%matrices_block(bb)
               if (IOwnPgrp(ptree, block_o%pgno)) then
                  n1 = MPI_Wtime()
                  !!!!! partial update butterflies at level llplus from left B1 = D^-1xB
                  if (llplus /= Lplus) then
                     rank0 = block_o%rankmax
                     rate = option%rankrate !1.2d0
                     level_butterfly = block_o%level_butterfly
                     Bplus%ind_ll = llplus
                     Bplus%ind_bk = bb
                     call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, Bplus, Bplus_block_MVP_diagBinvB_dat, error, 'L update', option, stats, ptree, msh,operand1= msh)
                     stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                     error_inout = max(error_inout, error)
                  endif
                  n2 = MPI_Wtime()
                  stats%Time_PartialUpdate = stats%Time_PartialUpdate + n2 - n1

                  n1 = MPI_Wtime()
                  !!!!! invert I+B1 to be I+B2
                  level_butterfly = block_o%level_butterfly
                  call BF_inverse_partitionedinverse_IplusButter(block_o, level_butterfly, 0, option, error, stats, ptree, msh, block_o%pgno)
                  error_inout = max(error_inout, error)
                  n2 = MPI_Wtime()
                  stats%Time_SMW = stats%Time_SMW + n2 - n1

                  n1 = MPI_Wtime()
                  if (llplus /= Lplus) then
                     rank0 = block_o%rankmax
                     rate = option%rankrate !1.2d0
                     level_butterfly = block_o%level_butterfly
                     Bplus%ind_ll = llplus
                     Bplus%ind_bk = bb
                     call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, Bplus, Bplus_block_MVP_BdiagBinv_dat, error, 'R update', option, stats, ptree, msh, operand1=msh)
                     stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                     error_inout = max(error_inout, error)
                  endif
                  n2 = MPI_Wtime()
                  stats%Time_PartialUpdate = stats%Time_PartialUpdate + n2 - n1
               endif
            end do
         end do


         do ll = 1, Bplus%Lplus
            Bplus%LL(ll)%rankmax = 0
            do bb = 1, Bplus%LL(ll)%Nbound
               Bplus%LL(ll)%rankmax = max(Bplus%LL(ll)%rankmax, Bplus%LL(ll)%matrices_block(bb)%rankmax)
            enddo
            call MPI_ALLREDUCE(MPI_IN_PLACE, Bplus%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
         end do

         rank_new_max = 0
         do ll = 1, Lplus
            rank_new_max = max(rank_new_max, Bplus%LL(ll)%rankmax)
         end do

         if (option%verbosity >= 1 .and. ptree%myid == ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%head) write (*, '(A10,I5,A6,I3,A8,I3,A11,Es14.7)') 'Mult No. ', rowblock, ' rank:', rank_new_max, ' L_butt:', Bplus%LL(1)%matrices_block(1)%level_butterfly, ' error:', error_inout

      endif

      return

   end subroutine Bplus_inverse_schur_partitionedinverse

   subroutine Bplus_inverse_schur_partitionedinverse_hss(bplus, option, stats, ptree, msh)




#ifdef HAVE_OPENMP
     use omp_lib
#endif


      implicit none

      integer ierr
      integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks, level_butterfly_loc
      integer num_col, num_row, level, mm, nn, ii, jj, tt, ll, llplus, bb, mmb
      character chara
      real(kind=8) T0
      type(matrixblock), pointer::block_o, block_off1, block_off2
      type(matrixblock), pointer::blocks_o_D
      type(matrixblock)::block_tmp
      type(blockplus)::Bplus
      integer rank_new_max
      real(kind=8):: rank_new_avr, error, err_avr, err_max, tol_used
      integer niter
      real(kind=8):: error_inout, rate, rankrate_inner, rankrate_outter
      integer itermax, ntry, cnt, cnt_partial
      real(kind=8):: n1, n2, n3, n4, Memory
      integer rank0, rank0_inner, rank0_outter, Lplus, level_BP, levelm, groupm_start, ij_loc, edge_s, edge_e, edge_first, idx_end_m_ref, idx_start_m_ref, idx_start_b, idx_end_b, idxs,idxe,groupm
      DT, allocatable:: matin(:, :), matout(:, :), matin_tmp(:, :), matout_tmp(:, :), matrixtemp(:,:)
      DT:: ctemp1, ctemp2
      integer, allocatable :: ipiv(:)
      type(Hoption)::option
      type(Hstat)::stats
      type(matrixblock):: agent_block
      type(blockplus):: agent_bplus
      type(proctree) :: ptree
      type(mesh) :: msh
      real(kind=8) flop

      ctemp1 = 1d0
      ctemp2 = 0d0
      Memory = 0
      error_inout = 0


      Lplus = Bplus%Lplus
      do llplus = Lplus, 1, -1
         call MPI_Barrier(ptree%Comm,ierr)
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'HSS inverse at level:', llplus
         do bb = 1, Bplus%LL(llplus)%Nbound
            block_o => Bplus%LL(llplus)%matrices_block(bb)
            if (IOwnPgrp(ptree, block_o%pgno)) then

               n1 = MPI_Wtime()
               !!!!! partial update butterflies at level llplus from left B1 = D^-1xB
               if (llplus /= Lplus) then


                  level_butterfly = block_o%level_butterfly
                  level_BP = Bplus%level
                  levelm = ceiling_safe(dble(level_butterfly)/2d0)
                  level_butterfly_loc = levelm
                  groupm_start=block_o%row_group*2**levelm
                  edge_s =msh%basis_group(block_o%row_group)%head
                  edge_e =msh%basis_group(block_o%row_group)%tail

#if 1
                  if(Bplus%LL(llplus+1)%Nbound>0)then
                  groupm = findgroup(edge_s, msh, levelm, block_o%row_group)
                  idxs =  groupm - Bplus%LL(llplus+1)%matrices_block(1)%row_group+1
                  groupm = findgroup(edge_e, msh, levelm, block_o%row_group)
                  idxe =  groupm - Bplus%LL(llplus+1)%matrices_block(1)%row_group+1
                  ! !call BF_MoveSingulartoLeft(block_o)
                  do ii=idxs,idxe
                     ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1
                     if(level_butterfly_loc==0)then
                        write(*,*)'level_butterfly_loc==0 not done'
                        stop
                     else
                        if (IOwnPgrp(ptree, Bplus%LL(llplus+1)%matrices_block(ii)%pgno)) then
                        call BF_extract_partial(block_o, level_butterfly_loc, ij_loc,Bplus%LL(llplus+1)%matrices_block(ii)%headm,Bplus%LL(llplus+1)%matrices_block(ii)%row_group, 'L', agent_block,Bplus%LL(llplus+1)%matrices_block(ii)%pgno,ptree)
                        rank0 = agent_block%rankmax
                        rate = option%rankrate !1.2d0
                        Bplus%ind_ll = llplus
                        Bplus%ind_bk = bb
                        if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect/max(1,block_o%level_butterfly/2)
                        call BF_randomized(agent_block%pgno, level_butterfly_loc, rank0, rate, agent_block, Bplus, Bplus_block_MVP_diagBinvBHSS_dat, error, 'L update', option, stats, ptree, msh, operand1=msh,vskip=.true.)
                        if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect*max(1,block_o%level_butterfly/2)
                        stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                        error_inout = max(error_inout, error)
                        call BF_ChangePattern(agent_block, 3, 2, stats, ptree)
                        ! call BF_MoveSingulartoRight(agent_block)

                        call BF_copyback_partial(block_o, level_butterfly_loc, ij_loc, 'L', agent_block,Bplus%LL(llplus+1)%matrices_block(ii)%pgno,ptree)

                        call BF_delete(agent_block,1)
                        endif
                     end if
                  end do
                  endif
#else
                  rank0 = block_o%rankmax
                  rate = option%rankrate !1.2d0
                  level_butterfly = block_o%level_butterfly
                  Bplus%ind_ll = llplus
                  Bplus%ind_bk = bb
                  if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect/max(1,block_o%level_butterfly/2)
                  call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, Bplus, Bplus_block_MVP_diagBinvBHSS_dat, error, 'L update', option, stats, ptree, msh, operand1=msh)
                  if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect*max(1,block_o%level_butterfly/2)
                  stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                  error_inout = max(error_inout, error)
#endif
               endif
               n2 = MPI_Wtime()
               stats%Time_PartialUpdate = stats%Time_PartialUpdate + n2 - n1


               n1 = MPI_Wtime()
               if (block_o%style == 1) then
#if HAVE_ZFP
                  if(option%use_zfp==1 .or. (option%use_zfp==2 .and. block_o%row_group /=block_o%col_group))call ZFP_Decompress(block_o%fullmat,block_o%FullmatZFP,block_o%M,block_o%N,tol_used,0)
#endif
#if 0
                  allocate (ipiv(block_o%M))
                  call getrff90(block_o%fullmat, ipiv, flop=flop)
                  stats%Flop_Factor = stats%Flop_Factor + flop
                  call getrif90(block_o%fullmat, ipiv, flop=flop)
                     ! do ii=1,block_o%M
                     !          block_o%fullmat(ii,ii) = block_o%fullmat(ii,ii)-1d0  ! this is needed for the later multplication as we assume the inverse is I+ something
                     ! enddo
                  stats%Flop_Factor = stats%Flop_Factor + flop
                  deallocate (ipiv)
#else
                  allocate(matrixtemp(block_o%M,block_o%M))
                  matrixtemp = block_o%fullmat
                  call GeneralInverse(block_o%M, block_o%M, matrixtemp, block_o%fullmat, BPACK_SafeEps, Flops=flop)
                  stats%Flop_Factor = stats%Flop_Factor + flop
                  deallocate(matrixtemp)
#endif
#if HAVE_ZFP
                  if(option%use_zfp==1 .or. (option%use_zfp==2 .and. block_o%row_group /=block_o%col_group))call ZFP_Compress(block_o%fullmat,block_o%FullmatZFP,block_o%M,block_o%N,option%tol_comp,0)
#endif
               else
                  !!!!! invert I+B1 to be I+B2
                  level_butterfly = block_o%level_butterfly
                  call BF_get_rank(block_o, ptree)
                  if (level_butterfly >= option%schulzlevel) then
                     call BF_inverse_schulziteration_IplusButter(block_o, error, option, stats, ptree, msh)
                  else
                     call BF_inverse_partitionedinverse_IplusButter(block_o, level_butterfly, 0, option, error, stats, ptree, msh, block_o%pgno)
                  endif
                  error_inout = max(error_inout, error)
               endif
               n2 = MPI_Wtime()
               stats%Time_SMW = stats%Time_SMW + n2 - n1

               n1 = MPI_Wtime()
               if (llplus /= Lplus) then
                  level_butterfly = block_o%level_butterfly
                  levelm = ceiling_safe(dble(level_butterfly)/2d0)
                  level_butterfly_loc = levelm
                  groupm_start=block_o%row_group*2**levelm
                  edge_s =msh%basis_group(block_o%row_group)%head
                  edge_e =msh%basis_group(block_o%row_group)%tail

#if 1
                  ! call BF_MoveSingulartoRight(block_o)


                  call BF_MoveSingular_Ker(block_o, 'T', floor_safe(dble(level_butterfly_loc)/2d0)+ block_o%level_half +1, block_o%level_half, ptree, stats, option%tol_rand)

                  if(Bplus%LL(llplus+1)%Nbound>0)then
                  groupm = findgroup(edge_s, msh, levelm, block_o%row_group)
                  idxs =  groupm - Bplus%LL(llplus+1)%matrices_block(1)%row_group+1
                  groupm = findgroup(edge_e, msh, levelm, block_o%row_group)
                  idxe =  groupm - Bplus%LL(llplus+1)%matrices_block(1)%row_group+1

                  do ii=idxs,idxe
                     ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1
                     if(level_butterfly_loc==0)then
                        write(*,*)'level_butterfly_loc==0 not done'
                        stop
                     else
                        if (IOwnPgrp(ptree,Bplus%LL(llplus+1)%matrices_block(ii)%pgno)) then
                        call BF_extract_partial(block_o, level_butterfly_loc, ij_loc, Bplus%LL(llplus+1)%matrices_block(ii)%headm,Bplus%LL(llplus+1)%matrices_block(ii)%row_group, 'R', agent_block,Bplus%LL(llplus+1)%matrices_block(ii)%pgno,ptree)

                        rank0 = agent_block%rankmax
                        rate = option%rankrate !1.2d0
                        Bplus%ind_ll = llplus
                        Bplus%ind_bk = bb
                        if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect/max(1,block_o%level_butterfly/2)
                        call BF_randomized(agent_block%pgno, level_butterfly_loc, rank0, rate, agent_block, Bplus, Bplus_block_MVP_BdiagBinvHSS_dat, error, 'R update', option, stats, ptree, msh, operand1=msh,uskip=.true.)
                        if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect*max(1,block_o%level_butterfly/2)
                        stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                        error_inout = max(error_inout, error)
                        call BF_ChangePattern(agent_block, 3, 1, stats, ptree)
                        ! call BF_MoveSingulartoLeft(agent_block)
                        call BF_copyback_partial(block_o, level_butterfly_loc, ij_loc, 'R', agent_block,Bplus%LL(llplus+1)%matrices_block(ii)%pgno,ptree)
                        call BF_delete(agent_block,1)
                        endif
                     end if
                  end do
                  endif
#else
                  rank0 = block_o%rankmax
                  rate = option%rankrate !1.2d0
                  level_butterfly = block_o%level_butterfly
                  Bplus%ind_ll = llplus
                  Bplus%ind_bk = bb
                  if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect/max(1,block_o%level_butterfly/2)
                  call BF_randomized(block_o%pgno, level_butterfly, rank0, rate, block_o, Bplus, Bplus_block_MVP_BdiagBinvHSS_dat, error, 'R update', option, stats, ptree, msh, operand1=msh)
                  if(option%format==3 .and. option%near_para<=0.1d0)option%tol_Rdetect = option%tol_Rdetect*max(1,block_o%level_butterfly/2)
                  stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
                  error_inout = max(error_inout, error)
#endif
               endif
               n2 = MPI_Wtime()
               stats%Time_PartialUpdate = stats%Time_PartialUpdate + n2 - n1
            endif
         end do
      end do


      do ll = 1, Bplus%Lplus
         Bplus%LL(ll)%rankmax = 0
         do bb = 1, Bplus%LL(ll)%Nbound
            if(IOwnPgrp(ptree,Bplus%LL(ll)%matrices_block(bb)%pgno))Bplus%LL(ll)%rankmax = max(Bplus%LL(ll)%rankmax, Bplus%LL(ll)%matrices_block(bb)%rankmax)
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, Bplus%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      end do

      rank_new_max = 0
      do ll = 1, Lplus
         rank_new_max = max(rank_new_max, Bplus%LL(ll)%rankmax)
      end do

      if (option%verbosity >= 1 .and. ptree%myid == ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%head) write (*, '(A14,A6,I6,A8,I3,A11,Es14.7)') 'HSS inverse: ', ' rank:', rank_new_max, ' L_butt:', Bplus%LL(1)%matrices_block(1)%level_butterfly, ' error:', error_inout

      return

   end subroutine Bplus_inverse_schur_partitionedinverse_hss

end module Bplus_factor
