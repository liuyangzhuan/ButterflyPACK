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
module Bplus_randomizedop
! use Utilites_randomized

use MISC_Utilities
use BPACK_Utilities

contains



subroutine BF_block_MVP_inverse_dat(ho_bf1,level,ii,trans,N,num_vect_sub,Vin,Vout,ptree,stats)
   use BPACK_DEFS


   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vin1(:,:),Vin2(:,:),Vout1(:,:),Vout2(:,:)
   DT :: ctemp1,ctemp2
   type(matrixblock),pointer::block_inv,block_schur,block_off1,block_off2
   integer groupn,groupm,mm,nn,ierr
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2

   ctemp1=1.0d0
   ctemp2=0.0d0


	block_off1 => ho_bf1%levels(level)%BP_inverse_update(ii*2-1)%LL(1)%matrices_block(1)
	block_off2 => ho_bf1%levels(level)%BP_inverse_update(ii*2)%LL(1)%matrices_block(1)
	block_schur => ho_bf1%levels(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)
	block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)


	nn=block_off1%N_loc
	mm=block_off1%M_loc

   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin

   ! call MPI_barrier(ptree%pgrp(block_inv%pgno)%Comm,ierr)
   n1 = OMP_get_wtime()
   allocate(Vin1(mm,num_vect_sub))
   call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin1,block_off1%M_p,0,block_off1%pgno,num_vect_sub,ptree)
   ! Vin1 = Vin(1:mm,1:num_vect_sub)

   allocate(Vin2(nn,num_vect_sub))
   call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin2,block_off1%N_p,block_off1%M,block_off1%pgno,num_vect_sub,ptree)
   ! Vin2 = Vin(1+mm:N,1:num_vect_sub)
   n2 = OMP_get_wtime()
   stats%Time_RedistV = stats%Time_RedistV + n2-n1

   allocate(Vout1(mm,num_vect_sub))
   Vout1=0
   allocate(Vout2(nn,num_vect_sub))
   Vout2=0


	if(trans=='N')then
		call BF_block_MVP_dat(block_off1,trans,mm,nn,num_vect_sub,&
		&Vin2,Vout1,ctemp1,ctemp2,ptree,stats)
		Vout1 = Vin1- Vout1
		Vout2 = Vin2

		! write(2111,*)abs(Vout)

		call BF_block_MVP_dat(block_schur,trans,mm,mm,num_vect_sub,&
		&Vout1,Vin1,ctemp1,ctemp2,ptree,stats)
		Vin1 = Vout1 + Vin1
		Vin2 = Vout2

		! write(2112,*)abs(Vin)

		call BF_block_MVP_dat(block_off2,trans,nn,mm,num_vect_sub,&
		&Vin1,Vout2,ctemp1,ctemp2,ptree,stats)
		Vout2 = Vin2 - Vout2
		Vout1 = Vin1

		! write(2113,*)abs(Vout)
		! stop

	else if(trans=='T')then
		call BF_block_MVP_dat(block_off2,trans,nn,mm,num_vect_sub,&
		&Vin2,Vout1,ctemp1,ctemp2,ptree,stats)
		Vout1 = Vin1 - Vout1
		Vout2 = Vin2

		call BF_block_MVP_dat(block_schur,trans,mm,mm,num_vect_sub,&
		&Vout1,Vin1,ctemp1,ctemp2,ptree,stats)
		Vin1 = Vout1 + Vin1
		Vin2 = Vout2

		call BF_block_MVP_dat(block_off1,trans,mm,nn,num_vect_sub,&
		&Vin1,Vout2,ctemp1,ctemp2,ptree,stats)
		Vout2 = Vin2 - Vout2
		Vout1 = Vin1

	end if

	! call MPI_barrier(ptree%pgrp(block_inv%pgno)%Comm,ierr)
	n1 = OMP_get_wtime()
	! Vout(1:mm,1:num_vect_sub) = Vout1
	call Redistribute1Dto1D(Vout1,block_off1%M_p,0,block_off1%pgno,Vout,block_inv%M_p,0,block_inv%pgno,num_vect_sub,ptree)
    ! Vout(1+mm:N,1:num_vect_sub) = Vout2
	call Redistribute1Dto1D(Vout2,block_off1%N_p,block_off1%M,block_off1%pgno,Vout,block_inv%M_p,0,block_inv%pgno,num_vect_sub,ptree)
   n2 = OMP_get_wtime()
   stats%Time_RedistV = stats%Time_RedistV + n2-n1

   Vin = Vin_tmp

   deallocate(Vin_tmp)
   deallocate(Vin1)
   deallocate(Vin2)
   deallocate(Vout1)
   deallocate(Vout2)

end subroutine BF_block_MVP_inverse_dat




subroutine BF_Delete_RandVect(chara,random,level_butterfly)
	use BPACK_DEFS
	implicit none
    character chara
	integer num_col,num_row,level,i,j
	type(RandomBlock):: random
	integer level_butterfly

	if(chara=='T')then
		do level=0, level_butterfly+2
			num_row=random%RandomVectorLL(level)%num_row
			num_col=random%RandomVectorLL(level)%num_col
			do j=1, num_col
				do i=1, num_row
					if(allocated(random%RandomVectorLL(level)%blocks(i,j)%matrix))deallocate (random%RandomVectorLL(level)%blocks(i,j)%matrix)
				enddo
			enddo
			deallocate (random%RandomVectorLL(level)%blocks)
		enddo
		deallocate (random%RandomVectorLL)
	else if(chara=='N')then
		do level=0, level_butterfly+2
			num_row=random%RandomVectorRR(level)%num_row
			num_col=random%RandomVectorRR(level)%num_col
			do j=1, num_col
				do i=1, num_row
					if(allocated(random%RandomVectorRR(level)%blocks(i,j)%matrix))deallocate (random%RandomVectorRR(level)%blocks(i,j)%matrix)
				enddo
			enddo
			deallocate (random%RandomVectorRR(level)%blocks)
		enddo
		deallocate (random%RandomVectorRR)
	end if
end subroutine BF_Delete_RandVect



subroutine BF_Init_randomized(level_butterfly,rankmax,groupm,groupn,block,block_rand,msh,ptree,option,nodataflag)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max,dimension_rank, dimension_m, dimension_n, blocks, groupm, groupm_start,groupn_start,groupn,index_j,index_i
    real(kind=8) a,b,c,d
    DT ctemp
	type(matrixblock)::block,block_rand
	DT, allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real(kind=8), allocatable:: Singular(:)
	integer rankmax
	type(mesh)::msh
	type(Hoption)::option
	type(proctree)::ptree
	integer nodataflag,level_final,level_half
	integer idx_r,inc_r,nr,idx_c,inc_c,nc

	! allocate (butterfly_block_randomized(1))

    block_rand%level_butterfly=level_butterfly
    num_blocks=2**level_butterfly
	dimension_rank= rankmax

	! level_half = BF_Switchlevel(level_butterfly,option)
	level_half = floor_safe(dble(level_butterfly)/2d0) ! from outer to inner


	block_rand%level_half=level_half


	! write(*,*)dimension_rank


    block_rand%dimension_rank=dimension_rank

    block_rand%style=2
    block_rand%row_group=groupm
    block_rand%col_group=groupn

	block_rand%M = block%M
	block_rand%N = block%N
	block_rand%headm = block%headm
	block_rand%headn = block%headn
	block_rand%level = block%level


	block_rand%M_loc = block%M_loc
	block_rand%N_loc = block%N_loc
	block_rand%pgno = block%pgno
	block_rand%pgno_db = block%pgno_db


	if(associated(block%N_p))then
		allocate(block_rand%N_p(size(block%N_p,1),2))
		block_rand%N_p = block%N_p
	endif
	if(associated(block%M_p))then
		allocate(block_rand%M_p(size(block%M_p,1),2))
		block_rand%M_p = block%M_p
	endif


	if(IOwnPgrp(ptree,block%pgno))then

		groupm_start=groupm*2**level_butterfly
		groupn_start=groupn*2**level_butterfly
		if (level_butterfly/=0) then
			allocate (block_rand%ButterflyKerl(level_butterfly))
		endif

		!****** row-wise ordering from right side
		do level=0, level_half
			if(level_butterfly==0)then
				if(level==0)then
					block_rand%ButterflyV%idx=1
					block_rand%ButterflyV%inc=1
					block_rand%ButterflyV%nblk_loc=1
					block_rand%ButterflyV%num_blk=num_blocks
				elseif(level==level_butterfly+1)then
					block_rand%ButterflyU%idx=1
					block_rand%ButterflyU%inc=1
					block_rand%ButterflyU%nblk_loc=1
					block_rand%ButterflyU%num_blk=num_blocks
				endif
			else
				call GetLocalBlockRange(ptree,block_rand%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')
				if(level==0)then
					block_rand%ButterflyV%idx=idx_c
					block_rand%ButterflyV%inc=inc_c
					block_rand%ButterflyV%nblk_loc=nc
					block_rand%ButterflyV%num_blk=num_blocks
				elseif(level==level_butterfly+1)then
					block_rand%ButterflyU%idx=idx_r
					block_rand%ButterflyU%inc=inc_r
					block_rand%ButterflyU%nblk_loc=nr
					block_rand%ButterflyU%num_blk=num_blocks
				else
					block_rand%ButterflyKerl(level)%num_row=2**level
					block_rand%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
					block_rand%ButterflyKerl(level)%idx_r=idx_r
					block_rand%ButterflyKerl(level)%inc_r=inc_r
					block_rand%ButterflyKerl(level)%nr=nr
					block_rand%ButterflyKerl(level)%idx_c=idx_c*2-1
					block_rand%ButterflyKerl(level)%inc_c=inc_c
					block_rand%ButterflyKerl(level)%nc=nc*2
				endif
			endif
		enddo

		!****** column-wise ordering from left side
		level_final=level_half+1
		do level=level_butterfly+1,level_final, -1
			if(level_butterfly==0)then
				if(level==0)then
					block_rand%ButterflyV%idx=1
					block_rand%ButterflyV%inc=1
					block_rand%ButterflyV%nblk_loc=1
					block_rand%ButterflyV%num_blk=num_blocks
				elseif(level==level_butterfly+1)then
					block_rand%ButterflyU%idx=1
					block_rand%ButterflyU%inc=1
					block_rand%ButterflyU%nblk_loc=1
					block_rand%ButterflyU%num_blk=num_blocks
				endif
			else
				call GetLocalBlockRange(ptree,block_rand%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')

				if(level==0)then
					block_rand%ButterflyV%idx=idx_c
					block_rand%ButterflyV%inc=inc_c
					block_rand%ButterflyV%nblk_loc=nc
					block_rand%ButterflyV%num_blk=num_blocks
				elseif(level==level_butterfly+1)then
					block_rand%ButterflyU%idx=idx_r
					block_rand%ButterflyU%inc=inc_r
					block_rand%ButterflyU%nblk_loc=nr
					block_rand%ButterflyU%num_blk=num_blocks
				else
					block_rand%ButterflyKerl(level)%num_row=2**level
					block_rand%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
					block_rand%ButterflyKerl(level)%idx_r=idx_r*2-1
					block_rand%ButterflyKerl(level)%inc_r=inc_r
					block_rand%ButterflyKerl(level)%nr=nr*2
					block_rand%ButterflyKerl(level)%idx_c=idx_c
					block_rand%ButterflyKerl(level)%inc_c=inc_c
					block_rand%ButterflyKerl(level)%nc=nc
				endif
			endif
		enddo


		if(nodataflag==0)then
			!****** row-wise ordering from right side
			do level=0, level_half
				if(level==0)then
					allocate(block_rand%ButterflyV%blocks(block_rand%ButterflyV%nblk_loc))
				elseif(level==level_butterfly+1)then
					allocate(block_rand%ButterflyU%blocks(block_rand%ButterflyU%nblk_loc))
				else
					allocate(block_rand%ButterflyKerl(level)%blocks(block_rand%ButterflyKerl(level)%nr,block_rand%ButterflyKerl(level)%nc))
				endif
			enddo

			!****** column-wise ordering from left side
			level_final=level_half+1
			do level=level_butterfly+1,level_final, -1
				if(level==0)then
					allocate(block_rand%ButterflyV%blocks(block_rand%ButterflyV%nblk_loc))
				elseif(level==level_butterfly+1)then
					allocate(block_rand%ButterflyU%blocks(block_rand%ButterflyU%nblk_loc))
				else
					allocate(block_rand%ButterflyKerl(level)%blocks(block_rand%ButterflyKerl(level)%nr,block_rand%ButterflyKerl(level)%nc))
				endif
			enddo
		endif
	endif

    return

end subroutine BF_Init_randomized



subroutine BF_Resolving_Butterfly_LL_dat(num_vect_sub,nth_s,nth_e,Ng,level,blocks,RandVectIn,RandVectOut,option,ptree,msh,stats)

   use BPACK_DEFS


   implicit none

   integer nth_s,nth_e,level
   integer i,j,k,num_blocks,num_row,num_col,ii,jj,mm,kk, rs,re,rank
   integer index_i, index_j, iter, vector1, vector2, direction, round, flag
   real(kind=8) a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   DT ctemp
   integer kmax
   type(Hoption)::option
   type(proctree)::ptree
   type(mesh)::msh
   type(Hstat)::stats
   DT:: RandVectIn(:,:),RandVectOut(:,:)

   real(kind=8), allocatable :: Singular(:)
   DT, allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   DT, allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,noe,Ng,dimension_nn,nn1,nn2,ieo,level_butterfly
   real(kind=8)::n1,n2
   type(matrixblock) :: blocks
   type(butterfly_vec) :: BFvec
   integer idx_r,inc_r,nr,idx_c,inc_c,nc

   level_butterfly=blocks%level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
   mm=num_vect_subsub !rank

	!********* multiply BF^C with vectors
	RandVectOut = conjg(cmplx(RandVectOut,kind=8))
	call BF_block_MVP_partial(blocks,'N',num_vect_sub,RandVectOut,BFvec,level-1,ptree,msh,stats)
	do j=1, BFvec%vec(level)%nc
		do i=1, BFvec%vec(level)%nr
			BFvec%vec(level)%blocks(i,j)%matrix=conjg(cmplx(BFvec%vec(level)%blocks(i,j)%matrix,kind=8))
		enddo
	enddo


	!********* compute row spaces and reconstruct blocks at level level
	do nth = nth_s,nth_e
		if (level==0) then
			index_i = 1
		else
			num_col=blocks%ButterflyKerl(level)%num_col
			index_i = ceiling_safe(nth*Ng*2d0/num_col)
		endif

		call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')
		if(mod(index_i-idx_r,inc_r)==0 .and. index_i>=idx_r .and. index_i<=idx_r+(nr-1)*inc_r)then ! I own this block row number
			!$omp parallel do default(shared) private(j,index_j)
			do j=1,nc
				index_j=idx_c+(j-1)*inc_c
				call BF_OneBlock_LL(index_i, index_j, level,num_vect_sub,num_vect_subsub,nth,nth_s,blocks,BFvec,option,stats)
			enddo
			!$omp end parallel do
		endif
	enddo

	!********* delete BFvec%vec(level), note that all the other levels have already been deleted in BF_block_MVP_partial
	do j=1, BFvec%vec(level)%nc
		do i=1, BFvec%vec(level)%nr
			deallocate(BFvec%vec(level)%blocks(i,j)%matrix)
		enddo
	enddo
	deallocate(BFvec%vec(level)%blocks)
	deallocate(BFvec%vec)

   return

end subroutine BF_Resolving_Butterfly_LL_dat



subroutine BF_OneBlock_LL(index_i, index_j,level,num_vect_sub,mm,nth,nth_s,blocks,BFvec,option,stats)
   use BPACK_DEFS


   implicit none
   type(matrixblock) :: blocks
   DT, allocatable :: matB(:,:)
   integer index_i,index_j,i,j,level,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,level_butterfly
   type(Hoption) :: option
   type(Hstat) :: stats
   real(kind=8)::flop,n1,n2
   real(kind=8)::Flops
   type(butterfly_vec) :: BFvec
   integer index_i_loc_k,index_j_loc_k,index_ii_loc,index_jj_loc,index_ii,index_jj

	Flops=0
	n1 = OMP_get_wtime()
   ! blocks => butterfly_block_randomized(1)
   level_butterfly=blocks%level_butterfly

	if (level==0) then
		index_jj_loc=(index_j-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1
		index_j_loc_k=(index_j-blocks%ButterflyV%idx)/blocks%ButterflyV%inc+1
		dimension_nn=size(BFvec%vec(level)%blocks(1,index_jj_loc)%matrix,1)
		allocate(matB(dimension_nn,mm))
		matB = BFvec%vec(level)%blocks(1,index_jj_loc)%matrix(1:dimension_nn,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
		call ComputeRange(dimension_nn,mm,matB,rank,1,option%tol_Rdetect,Flops=flop)
		Flops = Flops + flop
		if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
		if(allocated(blocks%ButterflyV%blocks(index_j_loc_k)%matrix))deallocate(blocks%ButterflyV%blocks(index_j_loc_k)%matrix)
		allocate(blocks%ButterflyV%blocks(index_j_loc_k)%matrix(dimension_nn,rank))
		blocks%ButterflyV%blocks(index_j_loc_k)%matrix = matB(1:dimension_nn,1:rank)
		deallocate(matB)


	elseif (level==level_butterfly+1) then
		write(*,*)'should not come here in BF_OneBlock_LL'
		stop
	else

		index_ii=int((index_i+1)/2) ; index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level)

		index_ii_loc=(index_ii-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level)
		index_jj_loc=(index_jj-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1

		index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
		index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

		nn1 = size(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,1)
		nn2 = size(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix,1)

		allocate (matB(nn1+nn2,mm))
		matB(1:nn1,1:mm) = BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
		matB(1+nn1:nn2+nn1,1:mm) = BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
		call ComputeRange(nn1+nn2,mm,matB,rank,1,option%tol_Rdetect,Flops=flop)
		Flops = Flops + flop
		if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
		allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix(rank,nn1))
		allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix(rank,nn2))

		call copymatT(matB(1:nn1,1:rank),blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,nn1,rank)
		call copymatT(matB(1+nn1:nn2+nn1,1:rank),blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,nn2,rank)
		deallocate(matB)
	endif

	n2 = OMP_get_wtime()
	stats%Time_random(5) = stats%Time_random(5) + n2-n1

	stats%Flop_Tmp = stats%Flop_Tmp + flops

end subroutine BF_OneBlock_LL




subroutine BF_Resolving_Butterfly_RR_dat(num_vect_sub,nth_s,nth_e,Ng,level,blocks,RandVectIn,RandVectOut,option,ptree,msh,stats)

   use BPACK_DEFS


   implicit none

   integer nth_s,nth_e,level
   integer i,j,k,num_blocks,num_row,num_col,ii,jj,mm,kk, rs,re,rank
   integer index_i, index_j, iter, vector1, vector2, direction, round, flag
   real(kind=8) a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   DT ctemp
   integer kmax
   type(Hoption)::option
   type(proctree)::ptree
   type(mesh)::msh
   type(Hstat)::stats
   DT:: RandVectIn(:,:),RandVectOut(:,:)

   real(kind=8), allocatable :: Singular(:)
   DT, allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   DT, allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,noe,Ng,dimension_nn,nn1,nn2,ieo,level_butterfly,level_left_start
   real(kind=8)::n1,n2
   type(matrixblock) :: blocks
   type(butterfly_vec) :: BFvec
   type(butterfly_vec) :: BFvec1
   integer idx_r,inc_r,nr,idx_c,inc_c,nc

   level_butterfly=blocks%level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
   mm=num_vect_subsub !rank

   level_left_start= blocks%level_half+1    !  check here later
   if(level_left_start>0 .and. level_left_start==level)then
		n1 = OMP_get_wtime()
		call BF_block_MVP_partial(blocks,'N',num_vect_sub,RandVectIn,BFvec1,level-1,ptree,msh,stats)
		n2 = OMP_get_wtime()
		call BF_all2all_matvec(blocks,BFvec1%vec(level),num_vect_sub,stats,ptree,level-1,'R','C')
		! time_halfbuttermul = time_halfbuttermul + n2-n1
	endif

	!********* multiply BF^C with vectors
	RandVectOut = conjg(cmplx(RandVectOut,kind=8))
	call BF_block_MVP_partial(blocks,'T',num_vect_sub,RandVectOut,BFvec,level+1,ptree,msh,stats)
	do j=1, BFvec%vec(level_butterfly-level+1)%nc
		do i=1, BFvec%vec(level_butterfly-level+1)%nr
			BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix=conjg(cmplx(BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix,kind=8))
		enddo
	enddo

	!********* compute column spaces and reconstruct blocks at level level
	do nth = nth_s,nth_e
		if (level==level_butterfly+1) then
			index_j = 1
		else
			num_row=blocks%ButterflyKerl(level)%num_row
			index_j = ceiling_safe(nth*Ng*2d0/num_row)
		endif

		call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')
		if(mod(index_j-idx_c,inc_c)==0 .and. index_j>=idx_c .and. index_j<=idx_c+(nc-1)*inc_c)then ! I own this block column number
			!$omp parallel do default(shared) private(i,index_i)
			do i=1,nr
				index_i=idx_r+(i-1)*inc_r
				call BF_OneBlock_RR(index_i, index_j, level,num_vect_sub,num_vect_subsub,nth,nth_s,blocks,BFvec,BFvec1,option,stats)
			enddo
			!$omp end parallel do
		endif
	enddo

	!********* delete BFvec%vec(level_butterfly-level+1), note that all the other levels have already been deleted in BF_block_MVP_partial
	do j=1, BFvec%vec(level_butterfly-level+1)%nc
		do i=1, BFvec%vec(level_butterfly-level+1)%nr
			deallocate(BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix)
		enddo
	enddo
	deallocate(BFvec%vec(level_butterfly-level+1)%blocks)
	deallocate(BFvec%vec)

	!********* delete BFvec1%vec(level), note that all the other levels have already been deleted in BF_block_MVP_partial
	if(level_left_start>0 .and. level_left_start==level)then
	do j=1, BFvec1%vec(level)%nc
		do i=1, BFvec1%vec(level)%nr
			deallocate(BFvec1%vec(level)%blocks(i,j)%matrix)
		enddo
	enddo
	deallocate(BFvec1%vec(level)%blocks)
	deallocate(BFvec1%vec)
	endif

   return

end subroutine BF_Resolving_Butterfly_RR_dat





subroutine BF_OneBlock_RR(index_i, index_j,level,num_vect_sub,mm,nth,nth_s,blocks,BFvec,BFvec1,option,stats)
   use BPACK_DEFS


   implicit none
   type(matrixblock) :: blocks
   DT, allocatable :: matA(:,:),matB(:,:),matC(:,:)
   integer index_i,index_j,i,j,level,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,level_butterfly,level_left_start
   type(Hoption) :: option
   type(Hstat) :: stats
   real(kind=8)::flop,n1,n2
   real(kind=8)::Flops
   type(butterfly_vec) :: BFvec,BFvec1
   integer index_i_loc_k,index_j_loc_k,index_i_loc_s,index_j_loc_s,index_ii_loc,index_jj_loc,index_ii,index_jj
   integer dimension_mm
	Flops=0

   ! blocks => butterfly_block_randomized(1)
   level_butterfly=blocks%level_butterfly
   level_left_start= blocks%level_half+1    !  check here later



	n1 = OMP_get_wtime()
	if (level==level_butterfly+1) then
		if(level==level_left_start)then
			index_ii_loc=(index_i-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1
			index_i_loc_k=(index_i-blocks%ButterflyU%idx)/blocks%ButterflyU%inc+1
			index_i_loc_s=(index_i-BFvec1%vec(level)%idx_r)/BFvec1%vec(level)%inc_r+1

			dimension_mm=size(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,1)%matrix,1)
			allocate(matB(mm,dimension_mm))

			call copymatT(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)

			rank = size(BFvec1%vec(level)%blocks(index_i_loc_s,1)%matrix,1)
			allocate(blocks%ButterflyU%blocks(index_i_loc_k)%matrix(dimension_mm,rank))
			allocate(matC(rank,dimension_mm),matA(mm,rank))
			call copymatT(BFvec1%vec(level)%blocks(index_i_loc_s,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,dimension_mm,matA,matB,matC,option%tol_LS,Flops=flop)
			Flops = Flops + flop
			! write(*,*)fnorm(matC,rank,dimension_mm),'U',level_left,level_butterfly

			call copymatT(matC,blocks%ButterflyU%blocks(index_i_loc_k)%matrix,rank,dimension_mm)
			deallocate(matB,matC,matA)
		else
			index_ii_loc=(index_i-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1
			index_i_loc_k=(index_i-blocks%ButterflyU%idx)/blocks%ButterflyU%inc+1
			dimension_mm=size(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,1)%matrix,1)
			allocate(matB(dimension_mm,mm))
			matB = BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
			! write(*,*)fnorm(matB,dimension_mm,mm),'Ubeforerange',dimension_mm,mm
			call ComputeRange(dimension_mm,mm,matB,rank,1,option%tol_Rdetect,Flops=flop)
			Flops = Flops + flop
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank

			allocate(blocks%ButterflyU%blocks(index_i_loc_k)%matrix(dimension_mm,rank))
			blocks%ButterflyU%blocks(index_i_loc_k)%matrix=matB(1:dimension_mm,1:rank)
			! write(*,*)fnorm(matB(1:dimension_mm,1:rank),dimension_mm,rank),'U'
			deallocate(matB)
		endif
	elseif (level==0) then
		write(*,*)'should not come here in BF_OneBlock_RR'
		stop
	else

		!!!!!! index_ii is global index in blocks%ButterflyKerl(level) and BFvec1%vec(level)

		index_ii=2*index_i-1; index_jj=int((index_j+1)/2) ;  !index_ii is global index in BFvec%vec(level_butterfly-level+1)

		index_ii_loc=(index_ii-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
		index_jj_loc=(index_jj-BFvec%vec(level_butterfly-level+1)%idx_c)/BFvec%vec(level_butterfly-level+1)%inc_c+1

		index_i_loc_k=(2*index_i-1-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
		index_j_loc_k=(index_j-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

		nn1 = size(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,1)
		nn2 = size(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc+1,index_jj_loc)%matrix,1)

		if(level==level_left_start)then
			index_i_loc_s=(index_i-BFvec1%vec(level)%idx_r)/BFvec1%vec(level)%inc_r+1 !index_i_loc_s is local index in BFvec1%vec(level)
			index_j_loc_s=(index_j-BFvec1%vec(level)%idx_c)/BFvec1%vec(level)%inc_c+1

			allocate (matB(mm,nn1+nn2))
			call copymatT(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
			call copymatT(BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc+1,index_jj_loc)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)

			rank = size(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix,1)
			allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix(nn2,rank))

			allocate(matC(rank,nn1+nn2),matA(mm,rank))
			call copymatT(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,option%tol_LS,Flops=flop)
			! write(*,*)fnorm(matA,mm,rank),fnorm(matB,mm,nn1+nn2),fnorm(matC,rank,nn1+nn2),'Rker'

			Flops = Flops + flop
			call copymatT(matC(1:rank,1:nn1),blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,rank,nn1)
			call copymatT(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,rank,nn2)
			deallocate(matB,matC,matA)

		else
			allocate (matB(nn1+nn2,mm))
			matB(1:nn1,1:mm)=BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
			matB(1+nn1:nn2+nn1,1:mm)=BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc+1,index_jj_loc)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm)
			call ComputeRange(nn1+nn2,mm,matB,rank,1,option%tol_Rdetect,Flops=flop)
			Flops = Flops + flop
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank

			allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix(nn2,rank))

			blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix = matB(1:nn1,1:rank)
			blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix = matB(1+nn1:nn2+nn1,1:rank)
			deallocate(matB)
		end if
	endif

	n2 = OMP_get_wtime()
	stats%Time_random(5) = stats%Time_random(5) + n2-n1

	stats%Flop_Tmp = stats%Flop_Tmp + Flops

end subroutine BF_OneBlock_RR


subroutine BF_randomized(pgno_large,level_butterfly,rank0,rankrate,blocks_o,operand,blackbox_MVP_dat,error_inout,strings,option,stats,ptree,msh,operand1)

    use BPACK_DEFS
	use MISC_Utilities
    use omp_lib

    implicit none

	integer level_c,rowblock,ierr
    integer blocks1, blocks2, blocks3, level_butterfly, rank0, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,r1,r2,r3,r3tmp,mn,rank
    character chara
    real(kind=8) T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock)::blocks_o
    type(matrixblock),pointer::blocks_A, blocks_B, blocks_C, blocks_D,block_tmp
    integer rank_new_max,rank_pre_max
	real(kind=8):: rank_new_avr,error,rankrate
	integer niter,groupm,groupn
	real(kind=8):: error_inout
	integer itermax,ntry
	real(kind=8):: n1,n2,Memory
	DT, allocatable::matrix_small(:,:),U1(:,:),V1(:,:),U2(:,:),V2(:,:),U3(:,:),V3(:,:),U3tmp(:,:),V3tmp(:,:),UUtmp(:,:),VVtmp(:,:),UU(:,:),VV(:,:),UUr(:,:),VVr(:,:)
	real(kind=8),allocatable :: Singular(:)
	DT, allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:),Vout3(:,:),Vout4(:,:),Vout(:,:),Vinter(:,:)
	DT::ctemp1,ctemp2
	DT,allocatable:: matin(:,:),matout(:,:),matsub_tmp(:,:)
	integer idx_start_m_ref
	class(*):: operand
	class(*),optional:: operand1
	character(*)  :: strings
	type(matrixblock),allocatable::block_rand(:)
	type(Hoption)::option
	type(Hstat)::stats
	procedure(BMatVec)::blackbox_MVP_dat
	type(proctree)::ptree
	type(mesh)::msh
	integer converged
	integer pgno_large


	ctemp1 = 1d0; ctemp2 = 0d0
	Memory = 0

	stats%Flop_Tmp=0
	converged=0

	do tt = 1,option%itermax

		rank_pre_max = ceiling_safe(rank0*option%rankrate**(tt-1))+1


		groupm=blocks_o%row_group
		groupn=blocks_o%col_group

		if(level_butterfly==0)then
			allocate (block_rand(1))
			call BF_Init_randomized(level_butterfly,rank_pre_max,groupm,groupn,blocks_o,block_rand(1),msh,ptree,option,1)
			call BF_Reconstruction_Lowrank(block_rand(1),blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)
		else


			allocate (block_rand(1))
			n1 = OMP_get_wtime()
			call BF_Init_randomized(level_butterfly,rank_pre_max,groupm,groupn,blocks_o,block_rand(1),msh,ptree,option,0)
			n2 = OMP_get_wtime()
			stats%Time_random(1) = stats%Time_random(1) + n2-n1
			n1 = OMP_get_wtime()
			call BF_Reconstruction_LL(block_rand(1),blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)
			n2 = OMP_get_wtime()
			time_tmp = time_tmp + n2-n1
			call BF_Reconstruction_RR(block_rand(1),blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)

		endif

		call BF_Test_Reconstruction_Error(block_rand(1),blocks_o,operand,blackbox_MVP_dat,error_inout,ptree,stats,operand1)
		n2 = OMP_get_wtime()

		call BF_get_rank(block_rand(1),ptree)
		call MPI_ALLREDUCE(MPI_IN_PLACE, error_inout, 1,MPI_double_precision, MPI_MAX, ptree%pgrp(pgno_large)%Comm,ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE, block_rand(1)%rankmax, 1,MPI_integer, MPI_MAX, ptree%pgrp(pgno_large)%Comm,ierr)

		if(ptree%MyID==ptree%pgrp(blocks_o%pgno)%head .and. option%verbosity>=2)write(*,'(A38,A6,I3,A8,I2,A8,I3,A7,Es14.7,A9,I5,A8,I5)')' '//TRIM(strings)//' ',' rank:',block_rand(1)%rankmax,' Ntrial:',tt,' L_butt:',block_rand(1)%level_butterfly,' error:',error_inout,' #sample:',rank_pre_max,' #nproc:',ptree%pgrp(block_rand(1)%pgno)%nproc

		!!!!*** terminate if 1. error small enough or 2. rank smaller than num_vec
		if(error_inout>option%tol_rand .and. block_rand(1)%rankmax==rank_pre_max)then
			call BF_get_rank(block_rand(1),ptree)
			rank_new_max = block_rand(1)%rankmax
			call BF_delete(block_rand(1),1)

			if(rank_pre_max>=min(block_rand(1)%M,block_rand(1)%N))then
				write(*,*)'rank_pre_max: ',rank_pre_max,' exceeding min(M,N): ',min(block_rand(1)%M,block_rand(1)%N)
				stop
			endif
			deallocate(block_rand)
        else

			call BF_delete(blocks_o,1)
			call BF_get_rank(block_rand(1),ptree)
			rank_new_max = block_rand(1)%rankmax
			call BF_copy_delete(block_rand(1),blocks_o,Memory)
			deallocate(block_rand)
			converged=1
			exit
		end if
	end do

	if(converged==0)then
		write(*,*)'randomized scheme not converged in '//TRIM(strings)//'. level: ',blocks_o%level_butterfly,error_inout,rank_new_max
		stop
	endif

    return

end subroutine BF_randomized



subroutine BF_Reconstruction_Lowrank(block_rand,blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock
    integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer ranks,rank,i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2

    real(kind=8)::n1,n2,flop


	class(*):: operand
	class(*),optional:: operand1

	type(proctree)::ptree
	type(mesh)::msh

	type(matrixblock)::blocks_o,block_rand


	type(Hoption)::option
	type(Hstat)::stats
	procedure(BMatVec)::blackbox_MVP_dat
	DT::ctemp1,ctemp2
	integer num_vect,level_butterfly,rmax
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:)
	real(kind=8), allocatable:: Singular(:)
	integer mm,nn,q,qq,Nloc,pp
	DT,pointer :: matQ2D(:,:),matQcA_trans2D(:,:),matQUt2D(:,:),UU(:,:),VV(:,:)
	integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9),descQUt2D(9)
	integer tempi,ctxt,info,iproc,jproc,myi,myj,myArows,myAcols,myrow,mycol,nprow,npcol,M,N,mnmin

	stats%Flop_Tmp=0

	level_butterfly=0
	num_vect = block_rand%dimension_rank
	rmax = num_vect
	mm=blocks_o%M_loc
	nn=blocks_o%N_loc


	if(mm>0)then
		allocate(RandVectInR(block_rand%N_loc,num_vect))
		RandVectInR=0
		allocate(RandVectOutR(block_rand%M_loc,num_vect))
		RandVectOutR=0
		call RandomMat(nn,num_vect,min(nn,num_vect),RandVectInR,1)
	endif

	call blackbox_MVP_dat(operand,blocks_o,'N',mm,nn,num_vect,RandVectInR,RandVectOutR,cone,czero,ptree,stats,operand1)

	! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
	do qq=1,option%powiter
		if(mm>0)RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
		call blackbox_MVP_dat(operand,blocks_o,'T',mm,nn,num_vect,RandVectOutR,RandVectInR,cone,czero,ptree,stats,operand1)
		if(mm>0)RandVectInR=conjg(cmplx(RandVectInR,kind=8))
		call blackbox_MVP_dat(operand,blocks_o,'N',mm,nn,num_vect,RandVectInR,RandVectOutR,cone,czero,ptree,stats,operand1)
	enddo

	! computation of range Q
	call PComputeRange(block_rand%M_p,num_vect,RandVectOutR,ranks,option%tol_Rdetect,ptree,block_rand%pgno,flop)
	stats%Flop_Tmp = stats%Flop_Tmp + flop

	! computation of B^T = (Q^c*A)^T
	if(mm>0)RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
	call blackbox_MVP_dat(operand,blocks_o,'T',mm,nn,num_vect,RandVectOutR,RandVectInR,cone,czero,ptree,stats,operand1)
	if(mm>0)RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))

	! computation of SVD B=USV and output A = (QU)*(SV)
	call PQxSVDTruncate(block_rand,RandVectOutR,RandVectInR,ranks,rank,option,stats,ptree,flop)
	stats%Flop_Tmp = stats%Flop_Tmp + flop

	if(mm>0)deallocate(RandVectOutR,RandVectInR)

end subroutine BF_Reconstruction_Lowrank






subroutine PQxSVDTruncate(block_rand,matQ,matQcA_trans,rmax,rank,option,stats,ptree,flops)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock
    integer rank,rmax,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2

    real(kind=8)::n1,n2,flop
    real(kind=8),optional::flops
	type(proctree)::ptree
	type(matrixblock)::block_rand
	type(Hoption)::option
	type(Hstat)::stats
	DT::matQ(:,:),matQcA_trans(:,:)
	integer num_vect,level_butterfly
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:)
	real(kind=8), allocatable:: Singular(:)
	integer q,qq,Nloc,pp
	DT,pointer :: matQ2D(:,:),matQcA_trans2D(:,:),matQUt2D(:,:),UU(:,:),VV(:,:)
	integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9),descQUt2D(9)
	integer tempi,ctxt,info,iproc,jproc,myi,myj,myArows,myAcols,myrow,mycol,nprow,npcol,M,N,mnmin
	integer ierr

	if(present(flops))flops = 0

	if(IOwnPgrp(ptree,block_rand%pgno))then

	!!!!**** generate 2D grid blacs quantities
	ctxt = ptree%pgrp(block_rand%pgno)%ctxt
	call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
	if(myrow/=-1 .and. mycol/=-1)then
		myArows = numroc_wp(block_rand%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rmax, nbslpk, mycol, 0, npcol)
		! write(*,*)ptree%MyID,'descQ2D',M, ranks(bb_inv*2-1+bb-1-Bidxs+1)
		call descinit( descQ2D, block_rand%M, rmax, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descQ2D')
		allocate(matQ2D(myArows,myAcols))
		matQ2D=0

		myArows = numroc_wp(block_rand%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rmax, nbslpk, mycol, 0, npcol)
		! write(*,*)ptree%MyID,'descQcA_trans2D',N, ranks(bb_inv*2-1+bb-1-Bidxs+1)
		call descinit( descQcA_trans2D, block_rand%N, rmax, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descQcA_trans2D')
		allocate(MatQcA_trans2D(myArows,myAcols))
		MatQcA_trans2D=0

		mnmin=min(block_rand%N,rmax)

		myArows = numroc_wp(block_rand%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
		allocate(UU(myArows,myAcols))
		! write(*,*)ptree%MyID,'descUU',N, mnmin
		call descinit( descUU, block_rand%N, mnmin, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descUU')
		UU=0

		myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rmax, nbslpk, mycol, 0, npcol)
		allocate(VV(myArows,myAcols))
		! write(*,*)ptree%MyID,'descVV', mnmin, ranks(bb_inv*2-1+bb-1-Bidxs+1)
		call descinit( descVV, mnmin, rmax, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descVV')
		VV=0

		allocate(Singular(mnmin))
		Singular=0

	else
		descQ2D(2)=-1
		descQcA_trans2D(2)=-1
		descUU(2)=-1
		descVV(2)=-1
		allocate(matQ2D(1,1))   ! required for Redistribute1Dto2D
		matQ2D=0
		allocate(matQcA_trans2D(1,1)) ! required for Redistribute1Dto2D
		matQcA_trans2D=0
		allocate(UU(1,1))  ! required for Redistribute2Dto1D
		UU=0
		allocate(VV(1,1))
		VV=0
	endif


!!!!**** redistribution into 2D grid
	call Redistribute1Dto2D(matQ,block_rand%M_p,0,block_rand%pgno,matQ2D,block_rand%M,0,block_rand%pgno,rmax,ptree)
	call Redistribute1Dto2D(matQcA_trans,block_rand%N_p,0,block_rand%pgno,matQcA_trans2D,block_rand%N,0,block_rand%pgno,rmax,ptree)


!!!!**** compute B^T=V^TS^TU^T
	rank=0
	if(myrow/=-1 .and. mycol/=-1)then
		call PSVD_Truncate(block_rand%N, rmax,matQcA_trans2D,descQcA_trans2D,UU,VV,descUU,descVV,Singular,option%tol_Rdetect,rank,ctxt,flop=flop)
		if(present(flops))flops = flops + flop/dble(nprow*npcol)
		do ii=1,rank
			call g2l(ii,rank,npcol,nbslpk,jproc,myj)
			if(jproc==mycol)then
				UU(:,myj) = UU(:,myj)*Singular(ii)
			endif
		enddo


		myArows = numroc_wp(block_rand%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		allocate(matQUt2D(myArows,myAcols))
		! write(*,*)'descQUt2D', M, rank
		call descinit( descQUt2D, block_rand%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descQUt2D')
		matQUt2D=0

		call pgemmf90('N','T',block_rand%M,rank,rmax,cone, matQ2D,1,1,descQ2D,VV,1,1,descVV,czero,matQUt2D,1,1,descQUt2D,flop=flop)
		if(present(flops))flops = flops + flop/dble(nprow*npcol)
	else
		allocate(matQUt2D(1,1)) ! required for Redistribute2Dto1D
	endif

	call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1,MPI_integer, MPI_MAX, ptree%pgrp(block_rand%pgno)%Comm,ierr)

	block_rand%rankmax = rank
	block_rand%rankmin = rank
	allocate(block_rand%ButterflyU%blocks(1))
	allocate(block_rand%ButterflyV%blocks(1))
	allocate(block_rand%ButterflyU%blocks(1)%matrix(block_rand%M_loc,rank))
	allocate(block_rand%ButterflyV%blocks(1)%matrix(block_rand%N_loc,rank))

	!!!!**** redistribution into 1D grid conformal to leaf sizes
	call Redistribute2Dto1D(matQUt2D,block_rand%M,0,block_rand%pgno,block_rand%ButterflyU%blocks(1)%matrix,block_rand%M_p,0,block_rand%pgno,rank,ptree)
	call Redistribute2Dto1D(UU,block_rand%N,0,block_rand%pgno,block_rand%ButterflyV%blocks(1)%matrix,block_rand%N_p,0,block_rand%pgno,rank,ptree)


	if(myrow/=-1 .and. mycol/=-1)then
		deallocate(Singular)
	endif
	deallocate(matQ2D)
	deallocate(MatQcA_trans2D)
	deallocate(UU)
	deallocate(VV)
	deallocate(matQUt2D)
	endif

end subroutine PQxSVDTruncate






subroutine PSVDTruncateSigma(block_rand,matQcA_trans,rmax,rank,Singular,option,stats,ptree,flops)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock
    integer rank,rmax,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2

    real(kind=8)::n1,n2,flop
    real(kind=8),optional::flops
	type(proctree)::ptree
	type(matrixblock)::block_rand
	type(Hoption)::option
	type(Hstat)::stats
	DT::matQcA_trans(:,:)
	integer num_vect,level_butterfly
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:)
	real(kind=8):: Singular(:)
	integer q,qq,Nloc,pp
	DT,pointer :: matQ2D(:,:),matQcA_trans2D(:,:),matQUt2D(:,:),UU(:,:),VV(:,:)
	integer descQcA_trans2D(9), descUU(9), descVV(9),descQUt2D(9)
	integer tempi,ctxt,info,iproc,jproc,myi,myj,myArows,myAcols,myrow,mycol,nprow,npcol,M,N,mnmin
	integer ierr

	if(present(flops))flops = 0
	Singular=0

	!!!!**** generate 2D grid blacs quantities
	ctxt = ptree%pgrp(block_rand%pgno)%ctxt
	call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
	if(myrow/=-1 .and. mycol/=-1)then

		myArows = numroc_wp(block_rand%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rmax, nbslpk, mycol, 0, npcol)
		! write(*,*)ptree%MyID,'descQcA_trans2D',N, ranks(bb_inv*2-1+bb-1-Bidxs+1)
		call descinit( descQcA_trans2D, block_rand%N, rmax, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descQcA_trans2D')
		allocate(MatQcA_trans2D(myArows,myAcols))
		MatQcA_trans2D=0

		mnmin=min(block_rand%N,rmax)

		myArows = numroc_wp(block_rand%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
		allocate(UU(myArows,myAcols))
		! write(*,*)ptree%MyID,'descUU',N, mnmin
		call descinit( descUU, block_rand%N, mnmin, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descUU')
		UU=0

		myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rmax, nbslpk, mycol, 0, npcol)
		allocate(VV(myArows,myAcols))
		! write(*,*)ptree%MyID,'descVV', mnmin, ranks(bb_inv*2-1+bb-1-Bidxs+1)
		call descinit( descVV, mnmin, rmax, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call assert(info==0,'descinit fail for descVV')
		VV=0
	else
		descQcA_trans2D(2)=-1
		descUU(2)=-1
		descVV(2)=-1
		allocate(matQcA_trans2D(1,1)) ! required for Redistribute1Dto2D
		matQcA_trans2D=0
		allocate(UU(1,1))  ! required for Redistribute2Dto1D
		UU=0
		allocate(VV(1,1))
		VV=0
	endif


!!!!**** redistribution into 2D grid
	call Redistribute1Dto2D(matQcA_trans,block_rand%N_p,0,block_rand%pgno,matQcA_trans2D,block_rand%N,0,block_rand%pgno,rmax,ptree)

!!!!**** compute singular values
	rank=0
	if(myrow/=-1 .and. mycol/=-1)then
		call PSVD_Truncate(block_rand%N, rmax,matQcA_trans2D,descQcA_trans2D,UU,VV,descUU,descVV,Singular,option%tol_Rdetect,rank,ctxt,flop=flop)
		if(present(flops))flops = flops + flop/dble(nprow*npcol)
	else
	endif

	call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1,MPI_integer, MPI_MAX, ptree%pgrp(block_rand%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE, Singular, rank,MPI_double_precision, MPI_MAX, ptree%pgrp(block_rand%pgno)%Comm,ierr)

	deallocate(MatQcA_trans2D)
	deallocate(UU)
	deallocate(VV)

end subroutine PSVDTruncateSigma




subroutine BF_Reconstruction_LL(block_rand,blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    DT ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row

    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real(kind=8)::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real(kind=8)::rank_new_avr,error
	DT,allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	real(kind=8):: error_inout
	integer,allocatable::perms(:)
	class(*):: operand
	class(*),optional:: operand1
	type(proctree)::ptree
	type(mesh)::msh
	type(matrixblock)::blocks_o,block_rand
	DT,allocatable:: RandVectIn(:,:),RandVectOut(:,:)
	type(Hoption)::option
	type(Hstat)::stats
	procedure(BMatVec)::blackbox_MVP_dat
	integer idx_r,inc_r,nr,idx_c,inc_c,nc

	mm=block_rand%M_loc
	nn=block_rand%N_loc

	level_butterfly=block_rand%level_butterfly
    num_blocks=2**level_butterfly
	dimension_rank =block_rand%dimension_rank
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here

	level_right_start = block_rand%level_half !  check here later

	do level = 0,level_right_start

		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-level)))   !  check here later
		Ng = 2**level_butterfly/Nsub


		Nbind = min(option%Nbundle,Nsub)
		num_vect_sub = num_vect_subsub*Nbind

		allocate(RandVectIn(block_rand%M_loc,num_vect_sub))
		RandVectIn=0
		allocate(RandVectOut(block_rand%N_loc,num_vect_sub))
		RandVectOut=0

		do ii = 1,Nsub/Nbind
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind

			n1 = OMP_get_wtime()

			call BF_Randomized_Vectors_dat('L',block_rand,RandVectIn,RandVectOut,blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,level,ptree,msh,stats,operand1)
			n2 = OMP_get_wtime()
			stats%Time_random(2) = stats%Time_random(2) + n2-n1
			! Time_Vector_inverse = Time_Vector_inverse + n2-n1

			n1 = OMP_get_wtime()
			call BF_Resolving_Butterfly_LL_dat(num_vect_sub,nth_s,nth_e,Ng,level,block_rand,RandVectIn,RandVectOut,option,ptree,msh,stats)
			n2 = OMP_get_wtime()
			stats%Time_random(3) = stats%Time_random(3) + n2-n1
		end do

		deallocate(RandVectIn)
		deallocate(RandVectOut)

	end do


    return

end subroutine BF_Reconstruction_LL


subroutine BF_Reconstruction_RR(block_rand,blocks_o,operand,blackbox_MVP_dat,operand1,option,stats,ptree,msh)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    DT ctemp, a, b
    character chara
	integer level_left_start,num_row,num_col
    real(kind=8)::n1,n2
	type(Hoption)::option
	type(Hstat)::stats

    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	! type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real(kind=8)::rank_new_avr
	DT,allocatable:: RandVectIn(:,:),RandVectOut(:,:)
	integer niter,level

	type(matrixblock)::blocks_o,block_rand
	class(*):: operand
	class(*),optional::operand1
	procedure(BMatVec)::blackbox_MVP_dat
	type(proctree)::ptree
	type(mesh)::msh

	level_butterfly=block_rand%level_butterfly
	mm=block_rand%M_loc
	nn=block_rand%N_loc

    num_blocks=2**level_butterfly
	dimension_rank =block_rand%dimension_rank
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here

    level_left_start= block_rand%level_half+1   !  check here later

	do level=level_butterfly+1,level_left_start,-1
		Nsub = NINT(2**ceiling_safe((level_butterfly)/2d0)/dble(2**(level-level_left_start)))    !  check here later
		Ng = 2**level_butterfly/Nsub

		Nbind = min(option%Nbundle,Nsub)
		num_vect_sub = num_vect_subsub*Nbind

		allocate(RandVectIn(block_rand%N_loc,num_vect_sub))
		RandVectIn=0
		allocate(RandVectOut(block_rand%M_loc,num_vect_sub))
		RandVectOut=0

		do ii = 1,Nsub/Nbind
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind

			n1 = OMP_get_wtime()
			call BF_Randomized_Vectors_dat('R',block_rand,RandVectIn,RandVectOut,blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,level,ptree,msh,stats,operand1)

			n2 = OMP_get_wtime()
			stats%Time_random(2) = stats%Time_random(2) + n2-n1
			! Time_Vector_inverse = Time_Vector_inverse + n2-n1

			n1 = OMP_get_wtime()
			call BF_Resolving_Butterfly_RR_dat(num_vect_sub,nth_s,nth_e,Ng,level,block_rand,RandVectIn,RandVectOut,option,ptree,msh,stats)
			n2 = OMP_get_wtime()
			stats%Time_random(3) = stats%Time_random(3) + n2-n1
		end do

		deallocate(RandVectIn)
		deallocate(RandVectOut)

	end do

    return

end subroutine BF_Reconstruction_RR



subroutine BF_Test_Reconstruction_Error(block_rand,block_o,operand,blackbox_MVP_dat,error,ptree,stats,operand1)

    use BPACK_DEFS
    implicit none

	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm
    integer mm,nn
    real(kind=8) a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    DT ctemp

    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real(kind=8)::error,tmp1,tmp2,norm1,norm2
	integer level_c,rowblock,dimension_m
	DT,allocatable::Vdref(:,:),Id(:,:),Vd(:,:)
	type(proctree)::ptree
	type(Hstat)::stats
	type(matrixblock)::block_o,block_rand
	class(*)::operand
	class(*),optional::operand1
	procedure(BMatVec)::blackbox_MVP_dat
	integer ierr
	error=0

		level_butterfly=block_rand%level_butterfly
		num_blocks=2**level_butterfly

		num_vect = 1

		mm=block_rand%M_loc
		nn=block_rand%N_loc

	if(IOwnPgrp(ptree,block_rand%pgno))then
		allocate (Vdref(mm,num_vect))
		Vdref=0
		allocate (Id(nn,num_vect))
		Id=0
		allocate (Vd(mm,num_vect))
		Vd=0

		call RandomMat(nn,num_vect,min(nn,num_vect),Id,0)
	endif

		call blackbox_MVP_dat(operand,block_o,'N',mm,nn,num_vect,Id,Vdref,cone,czero,ptree,stats,operand1)

	if(IOwnPgrp(ptree,block_rand%pgno))then
		call BF_block_MVP_dat(block_rand,'N',mm,nn,num_vect,Id,Vd,cone,czero,ptree,stats)

		tmp1 = fnorm(Vd-Vdref,mm,num_vect)**2d0
		call MPI_ALLREDUCE(tmp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%pgrp(block_rand%pgno)%Comm,ierr)
		tmp2 = fnorm(Vdref,mm,num_vect)**2d0
		call MPI_ALLREDUCE(tmp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%pgrp(block_rand%pgno)%Comm,ierr)
		error = sqrt(norm1)/sqrt(norm2)

		deallocate(Vdref)
		deallocate(Vd)
		deallocate(Id)
	endif
    return

end subroutine BF_Test_Reconstruction_Error




subroutine BF_Randomized_Vectors_dat(side,block_rand,RandVectIn,RandVectOut,blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth,ptree,msh,stats,operand1)

    use BPACK_DEFS

	use MISC_Utilities
    implicit none

    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk
    integer mm,nn,mm1,nn1,mn,level_butterfly
	character chara
	integer*8 idx_start

    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n

	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_right_start,level_left_start
	DT:: RandVectIn(:,:),RandVectOut(:,:)
	integer Nsub,Ng
	integer,allocatable:: mm_end(:),nn_end(:)

	DT,allocatable::matrix_small(:,:)
	class(*):: operand
	class(*),optional::operand1
	type(matrixblock)::blocks_o,block_rand
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	type(proctree)::ptree
	type(mesh)::msh
	type(Hstat)::stats
	procedure(BMatVec)::blackbox_MVP_dat
	integer idx_r,inc_r,nr,idx_c,inc_c,nc,idx_n_s,idx_n_e,idx_m_s,idx_m_e,idxs,idxe
	character side, trans
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)

    level_butterfly=block_rand%level_butterfly
	groupm_start=blocks_o%row_group*2**level_butterfly
	groupn_start=blocks_o%col_group*2**level_butterfly

	call GetLocalBlockRange(ptree,block_rand%pgno,0,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')
	idx_n_s=idx_c
	idx_n_e=idx_c+nc-1
	idxs=msh%basis_group(groupn_start+idx_n_s-1)%head
	idxe=msh%basis_group(groupn_start+idx_n_e-1)%tail
	nn = idxe-idxs+1


	call GetLocalBlockRange(ptree,block_rand%pgno,level_butterfly+1,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')
	idx_m_s=idx_r
	idx_m_e=idx_r+nr-1
	idxs=msh%basis_group(groupm_start+idx_m_s-1)%head
	idxe=msh%basis_group(groupm_start+idx_m_e-1)%tail
	mm = idxe-idxs+1



	RandVectIn=0
	RandVectOut=0
	dimension_rank =block_rand%dimension_rank

	if(side=='L')then
		trans='T'
		level_right_start = block_rand%level_half !  check here later
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later
		Ng = 2**level_butterfly/Nsub

		do nth= nth_s,nth_e
			!$omp parallel do default(shared) private(i,mm1,idxs,idxe)
			do i=(nth-1)*Ng+1, nth*Ng
				if(i>=idx_m_s .and. i<=idx_m_e)then
					idxs=msh%basis_group(groupm_start+i-1)%head-msh%basis_group(groupm_start+idx_m_s-1)%head+1
					idxe=msh%basis_group(groupm_start+i-1)%tail-msh%basis_group(groupm_start+idx_m_s-1)%head+1
					mm1 =idxe-idxs+1
					call RandomSubMat(idxs,idxe,(nth-nth_s)*num_vect_subsub+1,(nth-nth_s)*num_vect_subsub+num_vect_subsub,min(mm1,num_vect_subsub),RandVectIn,0)
				end if
			end do
			!$omp end parallel do
		end do

	elseif(side=='R')then
		trans='N'
		level_left_start = block_rand%level_half+1 !  check here later
		Nsub = NINT(2**ceiling_safe((level_butterfly)/2d0)/dble(2**(unique_nth-level_left_start)))    !  check here later
		Ng = 2**level_butterfly/Nsub

		do nth= nth_s,nth_e
			!$omp parallel do default(shared) private(i,nn1,idxs,idxe)
			do i=(nth-1)*Ng+1, nth*Ng
				if(i>=idx_n_s .and. i<=idx_n_e)then
					idxs=msh%basis_group(groupn_start+i-1)%head-msh%basis_group(groupn_start+idx_n_s-1)%head+1
					idxe=msh%basis_group(groupn_start+i-1)%tail-msh%basis_group(groupn_start+idx_n_s-1)%head+1
					nn1 =idxe-idxs+1
					call RandomSubMat(idxs,idxe,(nth-nth_s)*num_vect_subsub+1,(nth-nth_s)*num_vect_subsub+num_vect_subsub,min(nn1,num_vect_subsub),RandVectIn,0)
				end if
			end do
			!$omp end parallel do
		end do

	endif

	! get the left multiplied vectors
	call blackbox_MVP_dat(operand,blocks_o,trans,mm,nn,num_vect_sub,RandVectIn,RandVectOut,cone,czero,ptree,stats,operand1)

    return

end subroutine BF_Randomized_Vectors_dat



! blocks_D: D^-1 - I
! blocks_B: B
! blocks_C: C
! blocks_A: (A-BD^-1C)^-1 - I
subroutine BF_block_MVP_inverse_ABCD_dat(partitioned_block,block_o,trans,M_loc,N_loc,num_vect_sub,Vinin,Voutout,a,b,ptree,stats,operand1)
   use BPACK_DEFS


   implicit none
   integer level, ii, M, N, M_loc, N_loc, num_vect_sub, mv,nv
   character trans
   DT :: Vinin(:,:), Voutout(:,:)  ! the leading dimensions are block_o%M_loc or block_o%N_loc
   DT,allocatable :: Vin(:,:), Vout(:,:) ! the leading dimensions are blocks_B%M_loc+blocks_B%N_loc
   DT,allocatable :: Vin_tmp(:,:),Vout_tmp(:,:),Vbuff(:,:),V1(:,:),V2(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
   class(*),optional::operand1
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(partitioned_block)

   type is (matrixblock)
		call assert(M_loc==N_loc,'M/=N in BF_block_MVP_inverse_ABCD_dat')

		blocks_A => partitioned_block%sons(1,1)
		blocks_B => partitioned_block%sons(1,2)
		blocks_C => partitioned_block%sons(2,1)
		blocks_D => partitioned_block%sons(2,2)

		ctemp1=1.0d0
		ctemp2=0.0d0

		groupn=blocks_B%col_group    ! Note: row_group and col_group interchanged here
		nn=blocks_B%N_loc
		groupm=blocks_B%row_group    ! Note: row_group and col_group interchanged here
		mm=blocks_B%M_loc

		N=mm+nn

		!** convert the layout of input and output to those of child groups
		allocate(V1(mm,num_vect_sub))
		allocate(V2(nn,num_vect_sub))
		allocate(Vout(N,num_vect_sub))
		allocate(Vin(N,num_vect_sub))
		call Redistribute1Dto1D(Vinin,block_o%M_p,0,block_o%pgno,V1,blocks_A%M_p,0,blocks_A%pgno,num_vect_sub,ptree)
		if(mm>0)Vin(1:blocks_A%M_loc,:)=V1
		call Redistribute1Dto1D(Vinin,block_o%M_p,0,block_o%pgno,V2,blocks_D%M_p,blocks_A%M,blocks_D%pgno,num_vect_sub,ptree)
		if(nn>0)Vin(1+blocks_A%M_loc:N,:)=V2
		call Redistribute1Dto1D(Voutout,block_o%M_p,0,block_o%pgno,V1,blocks_A%M_p,0,blocks_A%pgno,num_vect_sub,ptree)
		if(mm>0)Vout(1:blocks_A%M_loc,:)=V1
		call Redistribute1Dto1D(Voutout,block_o%M_p,0,block_o%pgno,V2,blocks_D%M_p,blocks_A%M,blocks_D%pgno,num_vect_sub,ptree)
		if(nn>0)Vout(1+blocks_A%M_loc:N,:)=V2


		if(IOwnPgrp(ptree,blocks_A%pgno))then

			allocate(Vout_tmp(N,num_vect_sub))
			Vout_tmp = Vout
			allocate(Vin_tmp(N,num_vect_sub))
			Vin_tmp = Vin
			Vout = 0
			allocate(Vbuff(nn,num_vect_sub))
			Vbuff = 0

			if(trans=='N')then
				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vin(1+mm:N,1:num_vect_sub),Vbuff,ctemp1,ctemp2,ptree,stats)
				Vbuff = Vbuff + Vin(1+mm:N,1:num_vect_sub)
				call BF_block_MVP_dat(blocks_B,trans,mm,nn,num_vect_sub,&
				&Vbuff,Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)- Vout(1:mm,1:num_vect_sub)
				Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)

				if(isnan(fnorm(Vout,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD11N'
					stop
				end if

				! write(2111,*)abs(Vout)

				call BF_block_MVP_dat(blocks_A,trans,mm,mm,num_vect_sub,&
				&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vout(1+mm:mm+nn,1:num_vect_sub),Vin(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vin = Vout + Vin

				if(isnan(fnorm(Vin,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD22N'
					stop
				end if
				! write(2112,*)abs(Vin)

				call BF_block_MVP_dat(blocks_C,trans,nn,mm,num_vect_sub,&
				&Vin(1:mm,1:num_vect_sub),Vbuff,ctemp1,ctemp2,ptree,stats)
				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vbuff, Vout(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vout(1+mm:mm+nn,1:num_vect_sub) = Vout(1+mm:mm+nn,1:num_vect_sub) + Vbuff

				if(isnan(fnorm(Vout,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD33N'
					stop
				end if
				Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
				Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)

				! write(2113,*)abs(Vout)
				! stop

			else if(trans=='T')then

				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vin(1+mm:N,1:num_vect_sub),Vbuff,ctemp1,ctemp2,ptree,stats)
				Vbuff = Vbuff + Vin(1+mm:N,1:num_vect_sub)
				call BF_block_MVP_dat(blocks_C,trans,nn,mm,num_vect_sub,&
				&Vbuff,Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub) - Vout(1:mm,1:num_vect_sub)
				Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)

				if(isnan(fnorm(Vout,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD11T'
					stop
				end if

				call BF_block_MVP_dat(blocks_A,trans,mm,mm,num_vect_sub,&
				&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vout(1+mm:mm+nn,1:num_vect_sub),Vin(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vin = Vout + Vin

				if(isnan(fnorm(Vin,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD22T'
					stop
				end if

				call BF_block_MVP_dat(blocks_B,trans,mm,nn,num_vect_sub,&
				&Vin(1:mm,1:num_vect_sub),Vbuff,ctemp1,ctemp2,ptree,stats)
				call BF_block_MVP_dat(blocks_D,trans,nn,nn,num_vect_sub,&
				&Vbuff,Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2,ptree,stats)
				Vout(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub)+Vbuff

				if(isnan(fnorm(Vout,N,num_vect_sub)))then
					write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD33T'
					stop
				end if
				Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
				Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
			end if

		   Vin = Vin_tmp
		   Vout = Vout-Vin
		   Vout = a*Vout + b*Vout_tmp
		   deallocate(Vout_tmp)
		   deallocate(Vin_tmp)
		   deallocate(Vbuff)
	    endif

	    if(mm>0)V1=Vin(1:blocks_A%M_loc,:)
		call Redistribute1Dto1D(V1,blocks_A%M_p,0,blocks_A%pgno,Vinin,block_o%M_p,0,block_o%pgno,num_vect_sub,ptree)
		if(nn>0)V2=Vin(1+blocks_A%M_loc:N,:)
		call Redistribute1Dto1D(V2,blocks_D%M_p,blocks_A%M,blocks_D%pgno,Vinin,block_o%M_p,0,block_o%pgno,num_vect_sub,ptree)
		if(mm>0)V1=Vout(1:blocks_A%M_loc,:)
		call Redistribute1Dto1D(V1,blocks_A%M_p,0,blocks_A%pgno,Voutout,block_o%M_p,0,block_o%pgno,num_vect_sub,ptree)
		if(nn>0)V2=Vout(1+blocks_A%M_loc:N,:)
		call Redistribute1Dto1D(V2,blocks_D%M_p,blocks_A%M,blocks_D%pgno,Voutout,block_o%M_p,0,block_o%pgno,num_vect_sub,ptree)

	   deallocate(Vin)
	   deallocate(V1)
	   deallocate(V2)
	   deallocate(Vout)


   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_inverse_ABCD_dat






! blocks_D: D^-1 - I
! blocks_B: B
! blocks_C: C
! blocks_A: (A-BD^-1C)^-1 - I
subroutine BF_block_MVP_inverse_A_minusBDinvC_dat(partitioned_block,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS


   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: V_tmp1(:,:),V_tmp2(:,:),Vin_tmp(:,:),Vout_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
   class(*),optional::operand1
	type(matrixblock)::block_o
	type(proctree)::ptree
	type(Hstat)::stats

   select TYPE(partitioned_block)

   type is (matrixblock)
		call assert(M==N,'M/=N in BF_block_MVP_inverse_A_minusBDinvC_dat')

		blocks_A => partitioned_block%sons(1,1)
		blocks_B => partitioned_block%sons(1,2)
		blocks_C => partitioned_block%sons(2,1)
		blocks_D => partitioned_block%sons(2,2)


		groupn=blocks_B%col_group    ! Note: row_group and col_group interchanged here
		nn=blocks_B%N
		groupm=blocks_B%row_group    ! Note: row_group and col_group interchanged here
		mm=blocks_B%M


		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout


		Vout=0
		allocate(Vin_tmp(M,num_vect_sub))
		Vin_tmp=Vin

		allocate(V_tmp1(nn,num_vect_sub))
		V_tmp1 = 0
		allocate(V_tmp2(nn,num_vect_sub))
		V_tmp2 = 0

		if(trans=='N')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_C,'N',nn,mm,num_vect_sub,Vin_tmp,V_tmp1,ctemp1,ctemp2,ptree,stats)
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_D,'N',nn,nn,num_vect_sub,V_tmp1,V_tmp2,ctemp1,ctemp2,ptree,stats)
			V_tmp2 = V_tmp2 + V_tmp1
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_B,'N',mm,nn,num_vect_sub,V_tmp2,Vout,ctemp1,ctemp2,ptree,stats)
			ctemp1=1.0d0 ; ctemp2=-1.0d0
			call BF_block_MVP_dat(blocks_A,'N',mm,mm,num_vect_sub,Vin_tmp,Vout,ctemp1,ctemp2,ptree,stats)


		else if(trans=='T')then

			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_B,'T',mm,nn,num_vect_sub,Vin_tmp,V_tmp1,ctemp1,ctemp2,ptree,stats)
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_D,'T',nn,nn,num_vect_sub,V_tmp1,V_tmp2,ctemp1,ctemp2,ptree,stats)
			V_tmp2 = V_tmp2 + V_tmp1
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(blocks_C,'T',nn,mm,num_vect_sub,V_tmp2,Vout,ctemp1,ctemp2,ptree,stats)
			ctemp1=1.0d0 ; ctemp2=-1.0d0
			call BF_block_MVP_dat(blocks_A,'T',mm,mm,num_vect_sub,Vin_tmp,Vout,ctemp1,ctemp2,ptree,stats)
		end if



	   Vin = Vin_tmp


	   deallocate(Vin_tmp)
	   deallocate(V_tmp1)
	   deallocate(V_tmp2)

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_inverse_A_minusBDinvC_dat




subroutine BF_block_MVP_inverse_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS


   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:),tmpU(:,:),tmpV(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::block_off1,block_off2
   integer groupn,groupm,mm,nn
   class(*)::ho_bf1
   class(*),optional::operand1
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats
   integer rank
   select TYPE(ho_bf1)

   type is (hobf)

		block_off1 => ho_bf1%levels(ho_bf1%ind_lv)%BP_inverse_update(ho_bf1%ind_bk*2-1)%LL(1)%matrices_block(1)
		block_off2 => ho_bf1%levels(ho_bf1%ind_lv)%BP_inverse_update(ho_bf1%ind_bk*2)%LL(1)%matrices_block(1)

		groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here
		nn=block_off1%N_loc
		groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here
		mm=block_off1%M_loc

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout






		! rank=1
		! allocate(tmpU(mm,rank))
		! tmpU=1d0
		! allocate(tmpV(rank,mm))
		! tmpV=1d0
		! allocate(Vbuff(rank,num_vect_sub))
		! Vbuff=0
		! if(trans=='N')then
			! call gemmf90(tmpV,rank,Vin,mm,Vbuff,rank,'N','N',rank,num_vect_sub,mm,cone,czero)
			! call gemmf90(tmpU,mm,Vbuff,rank,Vout,M,'N','N',mm,num_vect_sub,rank,cone,czero)
		! else
			! call gemmf90(tmpU,mm,Vin,mm,Vbuff,rank,'T','N',rank,num_vect_sub,mm,cone,czero)
			! call gemmf90(tmpV,rank,Vbuff,rank,Vout,N,'T','N',mm,num_vect_sub,rank,cone,czero)
		! endif

		! deallocate(tmpU)
		! deallocate(tmpV)
		! deallocate(Vbuff)







		allocate(Vin_tmp(N,num_vect_sub))
		Vin_tmp = Vin
		Vout = 0

		allocate(Vbuff(nn,num_vect_sub))
		Vbuff = 0

		if(trans=='N')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(block_off2,'N',nn,mm,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,ptree,stats)
			call BF_block_MVP_dat(block_off1,'N',mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats)
			Vout = -Vout

		else if(trans=='T')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call BF_block_MVP_dat(block_off1,'T',mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,ptree,stats)
			call BF_block_MVP_dat(block_off2,'T',nn,mm,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats)
			Vout = -Vout
		end if

	   Vin = Vin_tmp

	   deallocate(Vin_tmp)
	   deallocate(Vbuff)

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_inverse_minusBC_dat




subroutine BF_block_MVP_schulz_dat(schulz_op,block_Xn,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS


   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vbuff1(:,:),Vout_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock)::block_Xn
   integer groupn,groupm,mm,nn
   class(*)::schulz_op
   class(*),optional::operand1
   type(matrixblock)::block_o
   real(kind=8)::scale_new
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(schulz_op)
   type is (schulz_operand)
	   select TYPE(operand1)
	   type is (integer)

			groupn=block_Xn%col_group    ! Note: row_group and col_group interchanged here
			nn=block_Xn%N_loc
			groupm=block_Xn%row_group    ! Note: row_group and col_group interchanged here
			mm=block_Xn%M_loc
			call assert(mm==nn,'block nonsquare')

			if(schulz_op%order==2)then
				mv=size(Vout,1)
				nv=size(Vout,2)
				allocate(Vout_tmp(mv,nv))
				Vout_tmp = Vout

				allocate(Vin_tmp(N,num_vect_sub))
				Vin_tmp = Vin
				Vout = 0

				allocate(Vbuff(nn,num_vect_sub))
				Vbuff = 0


				if(trans=='N')then
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! XnR
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,ptree,stats,operand1)

					! AXnR
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats)
					Vout = 	Vbuff+Vout

					! (2-AXn)R
					Vbuff = 2*Vin-Vout

					! Xn(2-AXn)R
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats,operand1)

				else if(trans=='T')then
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! RXn
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2,ptree,stats,operand1)

					! RXnA
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vout,Vbuff,ctemp1,ctemp2,ptree,stats)
					Vbuff=	Vout+Vbuff

					! RXnAXn
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vin,ctemp1,ctemp2,ptree,stats,operand1)

					! 2RXn - RXnAXn
					Vout = 2*Vout - Vin
				end if


				Vin = Vin_tmp

				scale_new=schulz_op%scale*(2-schulz_op%scale)

				Vout = Vout-Vin*scale_new ! Xn(2-AXn)-I

				deallocate(Vin_tmp)
				deallocate(Vbuff)

			else if(schulz_op%order==3)then

				mv=size(Vout,1)
				nv=size(Vout,2)
				allocate(Vout_tmp(mv,nv))
				Vout_tmp = Vout

				allocate(Vin_tmp(nn,num_vect_sub))
				Vin_tmp = Vin
				Vout = 0

				allocate(Vbuff(nn,num_vect_sub))
				Vbuff = 0
				allocate(Vbuff1(nn,num_vect_sub))
				Vbuff1 = 0



				if(trans=='N')then
					ctemp1=1.0d0 ; ctemp2=0.0d0

					! AXnR
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,ptree,stats,operand1)
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vbuff1,ctemp1,ctemp2,ptree,stats)
					Vbuff1 = 	Vbuff+Vbuff1


					! (AXn)^2R
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff1,Vbuff,ctemp1,ctemp2,ptree,stats,operand1)
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats)
					Vout = 	Vout+Vbuff

					! (3-3AXn+(AXn)^2)R
					Vbuff1 = 3*Vin-3*Vbuff1+Vout

					! Xn(3-3AXn+(AXn)^2)R
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff1,Vout,ctemp1,ctemp2,ptree,stats,operand1)

				else if(trans=='T')then
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! RXn
					Vout=Vin
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vout,Vin,ctemp1,ctemp2,ptree,stats,operand1)

					! RXn*AXn
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,ptree,stats)
					Vbuff=	Vin+Vbuff
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vbuff1,ctemp1,ctemp2,ptree,stats,operand1)

					! RXn*(AXn)^2
					call BF_block_MVP_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff1,Vbuff,ctemp1,ctemp2,ptree,stats)
					Vbuff=	Vbuff1+Vbuff
					call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,ptree,stats,operand1)


					! RXn*(3-3AXn+(AXn)^2)
					Vout = 3*Vin - 3*Vbuff1+Vout
				end if


				Vin = Vin_tmp

				scale_new=schulz_op%scale*(3 - 3*schulz_op%scale + schulz_op%scale**2d0)

				Vout = Vout-Vin*scale_new ! Xn(2-AXn)-I

				deallocate(Vin_tmp)
				deallocate(Vbuff)
				deallocate(Vbuff1)

			endif


			Vout = a*Vout + b*Vout_tmp
			deallocate(Vout_tmp)

	   end select
   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_schulz_dat




subroutine BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS


   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans,trans_new
   real(kind=8)::eps,memory
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:),matrixtmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock)::block_Xn
   integer groupn,groupm,mm,nn
   class(*)::schulz_op
   class(*),optional::operand1
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(schulz_op)
   type is (schulz_operand)
		select TYPE(operand1)
		type is (integer)
			eps=0.8d0
			ctemp1=1d0
			ctemp2=0d0
			if(operand1==1)then ! X0
				Vin=conjg(cmplx(Vin,kind=8))
				if(trans=='N')trans_new='T'
				if(trans=='T')trans_new='N'
				call BF_block_MVP_dat(schulz_op%matrices_block,trans_new,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2,ptree,stats)
				Vout = Vout + Vin
				Vin=conjg(cmplx(Vin,kind=8))
				Vout=conjg(cmplx(Vout,kind=8))
				schulz_op%scale=(2d0-eps)/schulz_op%A2norm**2d0
				Vout = Vout*schulz_op%scale



				! if(trans=='N')trans_new='T'
				! if(trans=='T')trans_new='N'
				! call BF_block_MVP_dat(schulz_op%matrices_block,trans_new,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				! Vout = Vin - Vout/schulz_op%A2norm



				! call BF_copy('N',schulz_op%matrices_block,block_o,memory)
				! call LR_SMW(block_o,memory)

				! call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				! Vout = Vout + Vin


			else ! Xn
				call BF_block_MVP_dat(block_Xn,trans,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2,ptree,stats)
				Vout=	Vout+Vin*schulz_op%scale
			endif

		end select
   end select

end subroutine BF_block_MVP_schulz_Xn_dat



subroutine BF_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,msh)

    use BPACK_DEFS

	use MISC_Utilities
    implicit none

	integer level_c,rowblock,unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,M,N,mv,nv
    integer mm,nn,mn,qq,pp,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,head,tail
    character trans
    ! real(kind=8) a,b,c,d
    DT ctemp,a,b


    ! type(vectorsblock), pointer :: random1, random2

    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_end_diag,idx_start_loc,idx_end_loc,rank
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vbuff(:,:),Vout_tmp(:,:),tmpU(:,:),tmpV(:,:)

	integer Nsub,Ng
	integer*8 idx_start
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n

	integer nth_s,nth_e,num_vect_sub,nth,level_right_start
	real(kind=8)::n1,n2
	DT :: Vin(:,:), Vout(:,:)
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)

	class(*)::ho_bf1
    class(*),optional::msh
	type(matrixblock)::block_o
	type(matrixblock),pointer::blocks
	type(proctree)::ptree
	type(Hstat)::stats

	call assert(present(msh),'operand1 cannot be skipped')

   select TYPE(msh)
   type is (mesh)
   select TYPE(ho_bf1)
   type is (hobf)

		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		if(trans=='N')then

			level_butterfly=block_o%level_butterfly

			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here
			nn=block_o%N_loc

			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here
			mm=block_o%M_loc

			allocate(Vbuff(mm,num_vect_sub))
			Vbuff=0

			! get the right multiplied vectors
			pp = ptree%myid-ptree%pgrp(block_o%pgno)%head+1
			idx_start_glo = block_o%headm + block_o%M_p(pp,1) -1


		n1 = OMP_get_wtime()
		  call BF_block_MVP_dat(block_o,'N',mm,nn,num_vect_sub,Vin,Vbuff,cone,czero,ptree,stats)
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1
			mm=block_o%M_loc
			allocate(vec_new(mm,num_vect_sub))

			do level = ho_bf1%Maxlevel+1,level_c+1,-1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = max((rowblock-1)*N_diag+1,ho_bf1%levels(level)%Bidxs)
				idx_end_diag = min(rowblock*N_diag,ho_bf1%levels(level)%Bidxe)
				vec_new = 0

				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_end_diag
					if(associated(ho_bf1%levels(level)%BP_inverse(ii)%LL))then
					blocks=>ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
					if(IOwnPgrp(ptree,blocks%pgno))then
					qq = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
					head = blocks%headm + blocks%M_p(qq,1) -1
					tail = head + blocks%M_loc - 1
					idx_start_loc = head-idx_start_glo+1
					idx_end_loc = tail-idx_start_glo+1

					if(level==ho_bf1%Maxlevel+1)then
						call Full_block_MVP_dat(blocks,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),cone,czero)
					else
						call BF_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
					endif
					endif
					endif
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1

				Vbuff = vec_new
			end do

			Vout = vec_new

			deallocate(Vbuff)
			deallocate(vec_new)

		else

			level_butterfly=block_o%level_butterfly

			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here
			nn=block_o%N_loc

			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here
			mm=block_o%M_loc

			allocate(Vbuff(mm,num_vect_sub))
			Vbuff=0


			! get the left multiplied vectors
			pp = ptree%myid-ptree%pgrp(block_o%pgno)%head+1
			idx_start_glo = block_o%headm + block_o%M_p(pp,1) -1

			allocate(vec_new(mm,num_vect_sub))
			Vbuff = Vin
			do level = level_c+1,ho_bf1%Maxlevel+1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = max((rowblock-1)*N_diag+1,ho_bf1%levels(level)%Bidxs)
				idx_end_diag = min(rowblock*N_diag,ho_bf1%levels(level)%Bidxe)
				vec_new = 0

				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_end_diag
					if(associated(ho_bf1%levels(level)%BP_inverse(ii)%LL))then
					blocks=>ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
					if(IOwnPgrp(ptree,blocks%pgno))then
					qq = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
					head = blocks%headm + blocks%M_p(qq,1) -1
					tail = head + blocks%M_loc - 1
					idx_start_loc = head-idx_start_glo+1
					idx_end_loc = tail-idx_start_glo+1

					if(level==ho_bf1%Maxlevel+1)then
						call Full_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),cone,czero)
					else
						call BF_block_MVP_inverse_dat(ho_bf1,level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
					endif
					endif
					endif
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1

				Vbuff = vec_new
			end do

			Vbuff = vec_new

			deallocate(vec_new)

			mm=block_o%M_loc
			nn=block_o%N_loc
			n1 = OMP_get_wtime()
			call BF_block_MVP_dat(block_o,'T',mm,nn,num_vect_sub,Vbuff,Vout,cone,czero,ptree,stats)
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1

			deallocate(Vbuff)

		end if

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop
   end select
   class default
		write(*,*)"unexpected type"
		stop
   end select
    return

end subroutine BF_block_MVP_Sblock_dat





! chara='m': block_1 x block_2
! chara='a': block_o + block_1
! chara='s': block_o - block_1
! chara='+': block_o + block_1 x block_2
! chara='-': block_o - block_1 x block_2
subroutine BF_block_MVP_Add_Multiply_dat(h_mat,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,chara)
   use BPACK_DEFS

   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::block_off1,block_off2
   integer groupn,groupm,groupk,mm,nn,kk
   class(*)::h_mat
   class(*),optional::chara
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(h_mat)
   type is (Hmat)
   select TYPE(chara)
   type is(character(*))

		block_off1 => h_mat%blocks_1
		block_off2 => h_mat%blocks_2

		mm=block_o%M
		nn=block_o%N
		kk=block_off1%N

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout
		allocate(Vin_tmp(N,num_vect_sub))
		Vin_tmp = Vin
		Vout = 0

		allocate(Vbuff(kk,num_vect_sub))
		Vbuff = 0


		if(trans=='N')then
			if(chara=='a' .or. chara=='s')then  ! block_o +- block_1
				call assert(kk==nn,'block dimensions do not match')
				Vbuff=Vin
			else
				call Hmat_block_MVP_dat(block_off2,trans,block_off2%headm,block_off2%headn,num_vect_sub,Vin,Vbuff,cone,ptree,stats)
			endif
			call Hmat_block_MVP_dat(block_off1,trans,block_off1%headm,block_off1%headn,num_vect_sub,Vbuff,Vout,cone,ptree,stats)
		else if(trans=='T')then
			call Hmat_block_MVP_dat(block_off1,trans,block_off1%headm,block_off1%headn,num_vect_sub,Vin,Vbuff,cone,ptree,stats)
			if(chara=='a' .or. chara=='s')then  ! block_o +- block_1
				call assert(kk==nn,'block dimensions do not match')
				Vout = Vbuff
			else
				call Hmat_block_MVP_dat(block_off2,trans,block_off2%headm,block_off2%headn,num_vect_sub,Vbuff,Vout,cone,ptree,stats)
			endif
	    endif

		if(chara=='+' .or. chara=='a')then ! block_o + block_1 x block_2 or block_o + block_1
			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vin,Vout,cone,cone,ptree,stats)
		else if(chara=='-' .or. chara=='s')then ! block_o - block_1 x block_2 or block_o - block_1
			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vin,Vout,cone,-cone,ptree,stats)
		else if(chara=='m')then ! block_1 x block_2
			!!!! nothing needs to be done here
		endif

	   Vin = Vin_tmp
	   deallocate(Vin_tmp)
	   deallocate(Vbuff)
	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   end select
   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_Add_Multiply_dat





subroutine BF_block_MVP_XLM_dat(blocks_l,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS

   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   integer groupn,groupm,groupk,mm,nn
   class(*)::blocks_l
   class(*),optional::operand1
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(blocks_l)
   type is (matrixblock)

		mm=block_o%M
		nn=block_o%N

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout
		allocate(Vin_tmp(N,num_vect_sub))
		Vin_tmp = Vin
		Vout = 0

		allocate(Vbuff(mm,num_vect_sub))
		Vbuff = 0


		if(trans=='N')then

			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vin,Vbuff,cone,czero,ptree,stats)
			Vout = Vbuff
			call Hmat_Lsolve(blocks_l,trans,blocks_l%headm,num_vect_sub,Vout,ptree,stats)


		else if(trans=='T')then
			Vbuff = Vin
			call Hmat_Lsolve(blocks_l,trans,blocks_l%headm,num_vect_sub,Vbuff,ptree,stats)
			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vbuff,Vout,cone,czero,ptree,stats)
	    endif


	   Vin = Vin_tmp
	   deallocate(Vin_tmp)
	   deallocate(Vbuff)
	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_XLM_dat



subroutine BF_block_MVP_XUM_dat(blocks_u,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
   use BPACK_DEFS

   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   integer groupn,groupm,groupk,mm,nn
   class(*)::blocks_u
   class(*),optional::operand1
   type(matrixblock)::block_o
   type(proctree)::ptree
   type(Hstat)::stats

   select TYPE(blocks_u)
   type is (matrixblock)


		mm=block_o%M
		nn=block_o%N


		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout
		allocate(Vin_tmp(N,num_vect_sub))
		Vin_tmp = Vin
		Vout = 0

		allocate(Vbuff(nn,num_vect_sub))
		Vbuff = 0


		if(trans=='N')then
			Vbuff = Vin
			call Hmat_Usolve(blocks_u,trans,blocks_u%headm,num_vect_sub,Vbuff,ptree,stats)
			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vbuff,Vout,cone,czero,ptree,stats)

		else if(trans=='T')then
			call BF_block_MVP_dat(block_o,trans,M,N,num_vect_sub,Vin,Vbuff,cone,czero,ptree,stats)
			Vout = Vbuff
			call Hmat_Usolve(blocks_u,trans,blocks_u%headm,num_vect_sub,Vout,ptree,stats)
	    endif


	   Vin = Vin_tmp
	   deallocate(Vin_tmp)
	   deallocate(Vbuff)
	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select


end subroutine BF_block_MVP_XUM_dat



subroutine Bplus_block_MVP_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real(kind=8) a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	type(proctree)::ptree
	type(Hstat)::stats
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1


	real(kind=8)::n2,n1

   select TYPE(bplus)

   type is (blockplus)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout
		Vout=0

		ctemp1=1.0d0 ; ctemp2=0.0d0

		! groupn=bplus%col_group  ! Note: row_group and col_group interchanged here
		! nn=bplus%LL(1)%matrices_block(1)%N
		! groupm=bplus%row_group  ! Note: row_group and col_group interchanged here
		! mm=bplus%LL(1)%matrices_block(1)%M


		call Bplus_block_MVP_dat(bplus,trans,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2,ptree,stats)


	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

    class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_Exact_dat




subroutine Bplus_block_MVP_Outter_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT,allocatable :: Vout_tmp(:,:)
	DT :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv

	type(Hstat)::stats
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1
	type(proctree)::ptree
	real(kind=8)::n2,n1

	ctemp3=-1.0d0 ; ctemp4=1.0d0


   select TYPE(bplus)

   type is (blockplus)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		call assert(present(operand1),'operand1 cannot be skipped')

		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,cone,czero,ptree,stats,operand1)

			call Bplus_block_MVP_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,ptree,stats,2,operand1%Lplus)
		end select

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_Outter_Exact_dat




subroutine Bplus_block_MVP_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn
	type(proctree)::ptree
	type(Hstat)::stats

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real(kind=8) a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2

	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1


	real(kind=8)::n2,n1

   select TYPE(ho_bf1)

   type is (hobf)
		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk

		bplus_off1 => ho_bf1%levels(level_c)%BP_inverse_update(2*rowblock-1)
		bplus_off2 => ho_bf1%levels(level_c)%BP_inverse_update(2*rowblock)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		if(trans=='N')then
			groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here
			nn=bplus_off1%LL(1)%matrices_block(1)%N_loc
			groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here
			mm=bplus_off1%LL(1)%matrices_block(1)%M_loc
			allocate(vec_new(nn,num_vect_sub))
			vec_new = 0

			! get the right multiplied vectors
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_dat(bplus_off2,'N',nn,mm,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)
			call Bplus_block_MVP_dat(bplus_off1,'N',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
			Vout = -Vout
			deallocate(vec_new)

	   else if(trans=='T')then
			groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here
			nn=bplus_off1%LL(1)%matrices_block(1)%N_loc
			groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here
			mm=bplus_off1%LL(1)%matrices_block(1)%M_loc

			allocate(vec_new(nn,num_vect_sub))
			vec_new=0

			! get the right multiplied vectors
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_dat(bplus_off1,'T',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)
			call Bplus_block_MVP_dat(bplus_off2,'T',nn,mm,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
			Vout = -Vout
			deallocate(vec_new)

	   end if

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

    class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_minusBC_dat




subroutine Bplus_block_MVP_Outter_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT,allocatable :: Vout_tmp(:,:)
	DT :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv
	type(proctree)::ptree
	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1
	type(Hstat)::stats

	real(kind=8)::n2,n1

	ctemp3=-1.0d0 ; ctemp4=1.0d0


   select TYPE(ho_bf1)

   type is (hobf)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		call assert(present(operand1),'operand1 cannot be skipped')

		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,cone,czero,ptree,stats,operand1)
			call Bplus_block_MVP_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,ptree,stats,2,operand1%Lplus)
		end select

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_Outter_minusBC_dat


subroutine Bplus_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,Ctemp,a,b
	type(blockplus),pointer::bplus_o
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,M,N,mv,nv,qq,head,tail,pp
	integer level_butterfly,groupm_diag
	! real(kind=8) a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_end_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)

	class(*):: ho_bf1
	type(matrixblock)::block_o
	type(matrixblock),pointer::blocks
	class(*),optional::operand1
	type(proctree)::ptree
	type(Hstat)::stats

	real(kind=8)::n2,n1


   select TYPE(ho_bf1)
   type is (hobf)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk

		if(trans=='N')then
			bplus_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)

			groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here
			nn=bplus_o%LL(1)%matrices_block(1)%N_loc

			groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here
			mm=bplus_o%LL(1)%matrices_block(1)%M_loc
			! get the right multiplied vectors
			pp = ptree%myid-ptree%pgrp(bplus_o%LL(1)%matrices_block(1)%pgno)%head+1
			idx_start_glo = bplus_o%LL(1)%matrices_block(1)%headm + bplus_o%LL(1)%matrices_block(1)%M_p(pp,1) -1

			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_dat(bplus_o,'N',mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2,ptree,stats)
			allocate(vec_new(mm,num_vect_sub))

			do level = ho_bf1%Maxlevel+1,level_c+1,-1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = max((rowblock-1)*N_diag+1,ho_bf1%levels(level)%Bidxs)
				idx_end_diag = min(rowblock*N_diag,ho_bf1%levels(level)%Bidxe)
				vec_new = 0


				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_end_diag
					if(associated(ho_bf1%levels(level)%BP_inverse(ii)%LL))then
					blocks=>ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
					if(IOwnPgrp(ptree,blocks%pgno))then
					qq = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
					head = blocks%headm + blocks%M_p(qq,1) -1
					tail = head + blocks%M_loc - 1
					idx_start_loc = head-idx_start_glo+1
					idx_end_loc = tail-idx_start_glo+1

					if(level==ho_bf1%Maxlevel+1)then
						call Full_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%	matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					else
						call Bplus_block_MVP_inverse_dat(ho_bf1, level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
					endif
					endif
					endif
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1


				Vout = vec_new
			end do
			deallocate(vec_new)


	   else if(trans=='T')then
			bplus_o => ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
			groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here
			nn=bplus_o%LL(1)%matrices_block(1)%N_loc
			groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here
			mm=bplus_o%LL(1)%matrices_block(1)%M_loc

			ctemp1=1.0d0 ; ctemp2=0.0d0
			! get the left multiplied vectors
			pp = ptree%myid-ptree%pgrp(bplus_o%LL(1)%matrices_block(1)%pgno)%head+1
			idx_start_glo = bplus_o%LL(1)%matrices_block(1)%headm + bplus_o%LL(1)%matrices_block(1)%M_p(pp,1) -1
			allocate(vec_old(mm,num_vect_sub))
			allocate(vec_new(mm,num_vect_sub))

			vec_old = Vin
			do level = level_c+1,ho_bf1%Maxlevel+1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = max((rowblock-1)*N_diag+1,ho_bf1%levels(level)%Bidxs)
				idx_end_diag = min(rowblock*N_diag,ho_bf1%levels(level)%Bidxe)
				vec_new = 0

				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_end_diag
					if(associated(ho_bf1%levels(level)%BP_inverse(ii)%LL))then
					blocks=>ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
					if(IOwnPgrp(ptree,blocks%pgno))then
					qq = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
					head = blocks%headm + blocks%M_p(qq,1) -1
					tail = head + blocks%M_loc - 1
					idx_start_loc = head-idx_start_glo+1
					idx_end_loc = tail-idx_start_glo+1
					if(level==ho_bf1%Maxlevel+1)then
						call Full_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					else
						call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
					endif
					endif
					endif
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1

				vec_old = vec_new
			end do
			deallocate(vec_new)
			n1 = OMP_get_wtime()
			call Bplus_block_MVP_dat(bplus_o,'T',mm,nn,num_vect_sub,vec_old,Vout,ctemp1,ctemp2,ptree,stats)
			n2 = OMP_get_wtime()
			deallocate(vec_old)
	   end if

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

    class default
		write(*,*)"unexpected type"
		stop
   end select

end subroutine Bplus_block_MVP_Sblock_dat



subroutine Bplus_block_MVP_Outter_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT,allocatable :: Vout_tmp(:,:)
	DT :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv
	type(Hstat)::stats

	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1
	type(proctree)::ptree
	real(kind=8)::n2,n1

	ctemp3=-1.0d0 ; ctemp4=1.0d0


   select TYPE(ho_bf1)

   type is (hobf)

		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout

		call assert(present(operand1),'operand1 cannot be skipped')

		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,cone,czero,ptree,stats,operand1)
			call Bplus_block_MVP_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,ptree,stats,2,operand1%Lplus)
		end select

	   Vout = a*Vout + b*Vout_tmp
	   deallocate(Vout_tmp)

   class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_Outter_Sblock_dat







subroutine Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,trans,N,num_vect_sub,Vin,Vout,ptree,stats)
   use BPACK_DEFS


   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vin1(:,:),Vin2(:,:),Vout1(:,:),Vout2(:,:)
   DT :: ctemp1,ctemp2
   type(matrixblock),pointer::block_o,block_inv,block_schur,block_off1,block_off2
   type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
   integer groupn,groupm,mm,nn,ierr
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2
   ctemp1=1.0d0
   ctemp2=0.0d0

	block_off1 => ho_bf1%levels(level)%BP_inverse_update(ii*2-1)%LL(1)%matrices_block(1)
	block_off2 => ho_bf1%levels(level)%BP_inverse_update(ii*2)%LL(1)%matrices_block(1)
	block_schur => ho_bf1%levels(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)
	block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)


	nn=block_off1%N_loc
	mm=block_off1%M_loc
   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin

   ! call MPI_barrier(ptree%pgrp(block_inv%pgno)%Comm,ierr)
   n1 = OMP_get_wtime()
   allocate(Vin1(mm,num_vect_sub))
   call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin1,block_off1%M_p,0,block_off1%pgno,num_vect_sub,ptree)
   ! Vin1 = Vin(1:mm,1:num_vect_sub)

   allocate(Vin2(nn,num_vect_sub))
   call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin2,block_off1%N_p,block_off1%M,block_off1%pgno,num_vect_sub,ptree)
   ! Vin2 = Vin(1+mm:N,1:num_vect_sub)
   n2 = OMP_get_wtime()
   stats%Time_RedistV = stats%Time_RedistV + n2-n1

   allocate(Vout1(mm,num_vect_sub))
   Vout1=0
   allocate(Vout2(nn,num_vect_sub))
   Vout2=0

	bplus_off1 => ho_bf1%levels(level)%BP_inverse_update(ii*2-1)
	bplus_off2 => ho_bf1%levels(level)%BP_inverse_update(ii*2)
	bplus_o => ho_bf1%levels(level)%BP_inverse_schur(ii)
	if(trans=='N')then
		call Bplus_block_MVP_dat(bplus_off1,trans,mm,nn,num_vect_sub,&
		&Vin2,Vout1,ctemp1,ctemp2,ptree,stats)
		Vout1 = Vin1- Vout1
		Vout2 = Vin2

		! write(2111,*)abs(Vout)

		call Bplus_block_MVP_dat(bplus_o,trans,mm,mm,num_vect_sub,&
		&Vout1,Vin1,ctemp1,ctemp2,ptree,stats)
		Vin1 = Vout1 + Vin1
		Vin2 = Vout2

		! write(2112,*)abs(Vin)

		call Bplus_block_MVP_dat(bplus_off2,trans,nn,mm,num_vect_sub,&
		&Vin1,Vout2,ctemp1,ctemp2,ptree,stats)
		Vout2 = Vin2 - Vout2
		Vout1 = Vin1

		! write(2113,*)abs(Vout)
		! stop

	else if(trans=='T')then
	! write(*,*)'good1'
		call Bplus_block_MVP_dat(bplus_off2,trans,nn,mm,num_vect_sub,&
		&Vin2,Vout1,ctemp1,ctemp2,ptree,stats)
		Vout1 = Vin1 - Vout1
		Vout2 = Vin2
	! write(*,*)'good2'
		call Bplus_block_MVP_dat(bplus_o,trans,mm,mm,num_vect_sub,&
		&Vout1,Vin1,ctemp1,ctemp2,ptree,stats)
		Vin1 = Vout1 + Vin1
		Vin2 = Vout2

	! write(*,*)'good3'
		call Bplus_block_MVP_dat(bplus_off1,trans,mm,nn,num_vect_sub,&
		&Vin1,Vout2,ctemp1,ctemp2,ptree,stats)
		Vout2 = Vin2 - Vout2
		Vout1 = Vin1
	! write(*,*)'good4'
	end if

	n1 = OMP_get_wtime()
	! Vout(1:mm,1:num_vect_sub) = Vout1
	call Redistribute1Dto1D(Vout1,block_off1%M_p,0,block_off1%pgno,Vout,block_inv%M_p,0,block_inv%pgno,num_vect_sub,ptree)
    ! Vout(1+mm:N,1:num_vect_sub) = Vout2
	call Redistribute1Dto1D(Vout2,block_off1%N_p,block_off1%M,block_off1%pgno,Vout,block_inv%M_p,0,block_inv%pgno,num_vect_sub,ptree)
   n2 = OMP_get_wtime()
   stats%Time_RedistV = stats%Time_RedistV + n2-n1

   Vin = Vin_tmp

   deallocate(Vin_tmp)
   deallocate(Vin1)
   deallocate(Vin2)
   deallocate(Vout1)
   deallocate(Vout2)

end subroutine Bplus_block_MVP_inverse_dat








subroutine Bplus_block_MVP_twoforward_dat(ho_bf1,level,ii,trans,N,num_vect_sub,Vin,Vout,a,b,ptree,stats)
   use BPACK_DEFS


   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vin1(:,:),Vin2(:,:),Vout1(:,:),Vout2(:,:)
   DT :: ctemp1,ctemp2, a, b
   type(matrixblock),pointer::block_o,block_inv,block_schur,block_off1,block_off2
   type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
   integer groupn,groupm,mm1,nn1,mm2,nn2,ierr,nin1,nout1,nin2,nout2,offin1,offout1,offin2,offout2
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2
   integer,pointer::Nin_p1(:,:),Nin_p2(:,:),Nout_p1(:,:),Nout_p2(:,:)
   ! ctemp1=1.0d0
   ! ctemp2=0.0d0

	block_off1 => ho_bf1%levels(level)%BP(ii*2-1)%LL(1)%matrices_block(1)
	block_off2 => ho_bf1%levels(level)%BP(ii*2)%LL(1)%matrices_block(1)
	block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
	bplus_off1 => ho_bf1%levels(level)%BP(ii*2-1)
	bplus_off2 => ho_bf1%levels(level)%BP(ii*2)

	mm1=block_off1%M_loc
	nn1=block_off1%N_loc

	mm2=block_off2%M_loc
	nn2=block_off2%N_loc



	if(trans=='N')then
		nin1 = nn1
		nout1 = mm1
		nin2 = nn2
		nout2 = mm2
		Nin_p1 => block_off1%N_p
		Nin_p2 => block_off2%N_p
		Nout_p1 => block_off1%M_p
		Nout_p2 => block_off2%M_p
		offin1 = block_off1%M
		offout1 = 0
		offin2 = 0
		offout2 = block_off1%M
	else
		nin1 = mm1
		nout1 = nn1
		nin2 = mm2
		nout2 = nn2
		Nin_p1 => block_off1%M_p
		Nin_p2 => block_off2%M_p
		Nout_p1 => block_off1%N_p
		Nout_p2 => block_off2%N_p
		offin1 = 0
		offout1 = block_off1%M
		offin2 = block_off1%M
		offout2 = 0
	endif

	! allocate(Vin_tmp(N,num_vect_sub))
	! Vin_tmp = Vin


	if(mm1>0)then
		allocate(Vin1(nin1,num_vect_sub))
		allocate(Vout1(nout1,num_vect_sub))
	endif
	if(mm2>0)then
		allocate(Vin2(nin2,num_vect_sub))
		allocate(Vout2(nout2,num_vect_sub))
	endif


	n1 = OMP_get_wtime()
	call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin1,Nin_p1,offin1,block_off1%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin2,Nin_p2,offin2,block_off2%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout,block_inv%N_p,0,block_inv%pgno,Vout1,Nout_p1,offout1,block_off1%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout,block_inv%N_p,0,block_inv%pgno,Vout2,Nout_p2,offout2,block_off2%pgno,num_vect_sub,ptree)

    n2 = OMP_get_wtime()
    stats%Time_RedistV = stats%Time_RedistV + n2-n1

	if(mm1>0)then
		call Bplus_block_MVP_dat(bplus_off1,trans,mm1,nn1,num_vect_sub,Vin1,Vout1,a,b,ptree,stats)
	endif
	if(mm2>0)then
		call Bplus_block_MVP_dat(bplus_off2,trans,mm2,nn2,num_vect_sub,Vin2,Vout2,a,b,ptree,stats)
	endif

	n1 = OMP_get_wtime()
	call Redistribute1Dto1D(Vout1,Nout_p1,offout1,block_off1%pgno,Vout,block_inv%N_p,0,block_inv%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout2,Nout_p2,offout2,block_off2%pgno,Vout,block_inv%N_p,0,block_inv%pgno,num_vect_sub,ptree)
    n2 = OMP_get_wtime()
    stats%Time_RedistV = stats%Time_RedistV + n2-n1


	if(mm1>0)then
		deallocate(Vin1)
		deallocate(Vout1)
	endif
	if(mm2>0)then
		deallocate(Vin2)
		deallocate(Vout2)
	endif


   ! Vin = Vin_tmp
   ! deallocate(Vin_tmp)

end subroutine Bplus_block_MVP_twoforward_dat






subroutine BF_block_MVP_twoforward_dat(ho_bf1,level,ii,block_rand,trans,N,num_vect_sub,Vin,Vout,a,b,ptree,stats)
   use BPACK_DEFS


   implicit none
   type(matrixblock)::block_rand(:)
   integer level, ii, N, num_vect_sub
   character trans
   DT :: Vin(:,:), Vout(:,:)
   DT,allocatable :: Vin_tmp(:,:),Vin1(:,:),Vin2(:,:),Vout1(:,:),Vout2(:,:)
   DT :: ctemp1,ctemp2, a, b
   type(matrixblock),pointer::block_o,block_inv,block_schur
   type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
   integer groupn,groupm,mm1,nn1,mm2,nn2,ierr,nin1,nout1,nin2,nout2,offin1,offout1,offin2,offout2,Bidxs
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2
   integer,pointer::Nin_p1(:,:),Nin_p2(:,:),Nout_p1(:,:),Nout_p2(:,:)

    Bidxs = ho_bf1%levels(level)%Bidxs*2-1

	! block_off1 => block_rand(ii*2-1-Bidxs+1)
	! block_off2 => block_rand(ii*2-Bidxs+1)
	block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)

	mm1=block_rand(ii*2-1-Bidxs+1)%M_loc
	nn1=block_rand(ii*2-1-Bidxs+1)%N_loc

	mm2=block_rand(ii*2-Bidxs+1)%M_loc
	nn2=block_rand(ii*2-Bidxs+1)%N_loc



	if(trans=='N')then
		nin1 = nn1
		nout1 = mm1
		nin2 = nn2
		nout2 = mm2
		Nin_p1 => block_rand(ii*2-1-Bidxs+1)%N_p
		Nin_p2 => block_rand(ii*2-Bidxs+1)%N_p
		Nout_p1 => block_rand(ii*2-1-Bidxs+1)%M_p
		Nout_p2 => block_rand(ii*2-Bidxs+1)%M_p
		offin1 = block_rand(ii*2-1-Bidxs+1)%M
		offout1 = 0
		offin2 = 0
		offout2 = block_rand(ii*2-1-Bidxs+1)%M
	else
		nin1 = mm1
		nout1 = nn1
		nin2 = mm2
		nout2 = nn2
		Nin_p1 => block_rand(ii*2-1-Bidxs+1)%M_p
		Nin_p2 => block_rand(ii*2-Bidxs+1)%M_p
		Nout_p1 => block_rand(ii*2-1-Bidxs+1)%N_p
		Nout_p2 => block_rand(ii*2-Bidxs+1)%N_p
		offin1 = 0
		offout1 = block_rand(ii*2-1-Bidxs+1)%M
		offin2 = block_rand(ii*2-1-Bidxs+1)%M
		offout2 = 0
	endif


	if(mm1>0)then
		allocate(Vin1(nin1,num_vect_sub))
		allocate(Vout1(nout1,num_vect_sub))
	endif
	if(mm2>0)then
		allocate(Vin2(nin2,num_vect_sub))
		allocate(Vout2(nout2,num_vect_sub))
	endif


	n1 = OMP_get_wtime()
	call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin1,Nin_p1,offin1,block_rand(ii*2-1-Bidxs+1)%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vin,block_inv%N_p,0,block_inv%pgno,Vin2,Nin_p2,offin2,block_rand(ii*2-Bidxs+1)%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout,block_inv%N_p,0,block_inv%pgno,Vout1,Nout_p1,offout1,block_rand(ii*2-1-Bidxs+1)%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout,block_inv%N_p,0,block_inv%pgno,Vout2,Nout_p2,offout2,block_rand(ii*2-Bidxs+1)%pgno,num_vect_sub,ptree)

    n2 = OMP_get_wtime()
    stats%Time_RedistV = stats%Time_RedistV + n2-n1

	if(mm1>0)then
		call BF_block_MVP_dat(block_rand(ii*2-1-Bidxs+1),trans,mm1,nn1,num_vect_sub,Vin1,Vout1,a,b,ptree,stats)
	endif
	if(mm2>0)then
		call BF_block_MVP_dat(block_rand(ii*2-Bidxs+1),trans,mm2,nn2,num_vect_sub,Vin2,Vout2,a,b,ptree,stats)
	endif

	n1 = OMP_get_wtime()
	call Redistribute1Dto1D(Vout1,Nout_p1,offout1,block_rand(ii*2-1-Bidxs+1)%pgno,Vout,block_inv%N_p,0,block_inv%pgno,num_vect_sub,ptree)
	call Redistribute1Dto1D(Vout2,Nout_p2,offout2,block_rand(ii*2-Bidxs+1)%pgno,Vout,block_inv%N_p,0,block_inv%pgno,num_vect_sub,ptree)
    n2 = OMP_get_wtime()
    stats%Time_RedistV = stats%Time_RedistV + n2-n1


	if(mm1>0)then
		deallocate(Vin1)
		deallocate(Vout1)
	endif
	if(mm2>0)then
		deallocate(Vin2)
		deallocate(Vout2)
	endif

end subroutine BF_block_MVP_twoforward_dat







subroutine Bplus_block_MVP_BplusB_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real(kind=8) a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	! type(vectorsblock), pointer :: random1, random2

	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1
	type(proctree)::ptree
	type(Hstat)::stats
	real(kind=8)::n2,n1

   select TYPE(bplus)

   type is (blockplus)

	 select TYPE(operand1)
	   type is (matrixblock)

			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout


			level_butterfly = block_o%level_butterfly
			num_blocks=2**level_butterfly
			mm=0
			nn=0
			do i=1, num_blocks
				mm = mm+size(block_o%ButterflyU%blocks(i)%matrix,1)
				nn = nn+size(block_o%ButterflyV%blocks(i)%matrix,1)
			enddo



			allocate(vec_new(mm,num_vect_sub))
			vec_new = 0

			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call BF_block_MVP_dat(operand1,'N',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)

				call Bplus_block_MVP_dat(bplus,'N',mm,mm,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				Vout = Vout + vec_new

				deallocate(vec_new)

		   else if(trans=='T')then

				! get the left multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call Bplus_block_MVP_dat(bplus,'T',mm,mm,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)
				vec_new = vec_new + Vin

				call BF_block_MVP_dat(operand1,'T',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				deallocate(vec_new)

		   end if

		   Vout = a*Vout + b*Vout_tmp
		   deallocate(Vout_tmp)

		class default
			write(*,*)"unexpected type"
			stop

	   end select

    class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_BplusB_dat



subroutine Bplus_block_MVP_diagBinvB_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,msh)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,pp,test
	integer level_butterfly,groupm_diag,edge_first
	! real(kind=8) a,b,c,d
	integer idx_start_glo_sml,idx_start_glo_tmp,idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc,Maxgrp,pgno,row_group
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:),Vin_sml(:,:),Vout_sml(:,:),Vin_tmp1(:,:),Vout_tmp1(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	! type(vectorsblock), pointer :: random1, random2

	class(*):: bplus
	type(matrixblock)::block_o,block_tmp
	type(matrixblock),pointer::blocks,blocks_sml
	class(*),optional::msh
	type(proctree)::ptree
	type(Hstat)::stats
	real(kind=8)::n2,n1

	call assert(present(msh),'operand1 cannot be skipped')
   select TYPE(msh)
   type is (mesh)
   select TYPE(bplus)
   type is (blockplus)
			blocks=> bplus%LL(1)%matrices_block(1)

			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout
			Vout=0
			allocate(vec_new(mv,nv))
			vec_new = 0

			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call BF_block_MVP_dat(blocks,trans,M,N,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)

				do ii=1,Bplus%LL(2)%Nbound
					blocks_sml=>Bplus%LL(2)%matrices_block(ii)
					edge_first = blocks_sml%headm
					if(edge_first>=blocks%headm .and. edge_first<=blocks%headm+blocks%M-1)then


						! generate a tmporay block block_tmp whose layout is conformal to blocks, then redistribute data from block_tmp to blocks_sml
						block_tmp%M=blocks_sml%M
						block_tmp%headm=blocks_sml%headm
						block_tmp%N=blocks_sml%N
						block_tmp%headn=blocks_sml%headn
						block_tmp%row_group=blocks_sml%row_group
						block_tmp%col_group=blocks_sml%col_group
						block_tmp%level=blocks_sml%level


						Maxgrp=2**(ptree%nlevel)-1

						row_group=blocks%row_group
						pgno=blocks%pgno
						do while(row_group/=block_tmp%row_group)
							if(msh%basis_group(2*row_group)%head<=msh%basis_group(block_tmp%row_group)%head .and. msh%basis_group(2*row_group)%tail>=msh%basis_group(block_tmp%row_group)%tail)then
								row_group = 2*row_group
								if(pgno*2<=Maxgrp)pgno=pgno*2
							else
								row_group = 2*row_group+1
								if(pgno*2+1<=Maxgrp)pgno=pgno*2+1
							endif
						enddo
						block_tmp%pgno=pgno
						call ComputeParallelIndices(block_tmp,block_tmp%pgno,ptree,msh)

						if(blocks_sml%M_loc>0)then
						allocate(Vin_sml(blocks_sml%M_loc,num_vect_sub))
						allocate(Vout_sml(blocks_sml%M_loc,num_vect_sub))
						Vout_sml=0
						endif
						if(block_tmp%M_loc>0)then
							allocate(Vin_tmp1(block_tmp%M_loc,num_vect_sub))
							allocate(Vout_tmp1(block_tmp%M_loc,num_vect_sub))
							pp = ptree%myid-ptree%pgrp(block_tmp%pgno)%head+1
							idx_start_glo_tmp = block_tmp%headm + block_tmp%M_p(pp,1) -1
							pp = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
							idx_start_glo = blocks%headm + blocks%M_p(pp,1) -1
							idx_start_loc = idx_start_glo_tmp-idx_start_glo+1
							idx_end_loc = idx_start_loc + block_tmp%M_loc-1
							Vin_tmp1 = vec_new(idx_start_loc:idx_end_loc,:)
						endif
						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vin_tmp1,block_tmp%M_p,0,block_tmp%pgno,Vin_sml,blocks_sml%M_p,0,blocks_sml%pgno,num_vect_sub,ptree)



						if(IOwnPgrp(ptree,blocks_sml%pgno))then
							call BF_block_MVP_dat(blocks_sml,trans,blocks_sml%M_loc,blocks_sml%N_loc,num_vect_sub,Vin_sml,Vout_sml,ctemp1,ctemp2,ptree,stats)
						endif

						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vout_sml,blocks_sml%M_p,0,blocks_sml%pgno,Vout_tmp1,block_tmp%M_p,0,block_tmp%pgno,num_vect_sub,ptree)
						if(block_tmp%M_loc>0)then
							Vout(idx_start_loc:idx_end_loc,:) = Vout_tmp1
							deallocate(Vin_tmp1)
							deallocate(Vout_tmp1)
						endif

						if(blocks_sml%M_loc>0)then
							deallocate(Vout_sml)
							deallocate(Vin_sml)
						endif
						call BF_delete(block_tmp,1)


					endif
				enddo
				Vout = Vout + vec_new

				deallocate(vec_new)

		   else if(trans=='T')then
				! get the left multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				do ii=1,Bplus%LL(2)%Nbound
					blocks_sml=>Bplus%LL(2)%matrices_block(ii)
					edge_first = blocks_sml%headm
					if(edge_first>=blocks%headm .and. edge_first<=blocks%headm+blocks%M-1)then




						! generate a tmporay block block_tmp whose layout is conformal to blocks, then redistribute data from block_tmp to blocks_sml
						block_tmp%M=blocks_sml%M
						block_tmp%headm=blocks_sml%headm
						block_tmp%N=blocks_sml%N
						block_tmp%headn=blocks_sml%headn
						block_tmp%row_group=blocks_sml%row_group
						block_tmp%col_group=blocks_sml%col_group
						block_tmp%level=blocks_sml%level


						Maxgrp=2**(ptree%nlevel)-1

						row_group=blocks%row_group
						pgno=blocks%pgno
						do while(row_group/=block_tmp%row_group)
							if(msh%basis_group(2*row_group)%head<=msh%basis_group(block_tmp%row_group)%head .and. msh%basis_group(2*row_group)%tail>=msh%basis_group(block_tmp%row_group)%tail)then
								row_group = 2*row_group
								if(pgno*2<=Maxgrp)pgno=pgno*2
							else
								row_group = 2*row_group+1
								if(pgno*2+1<=Maxgrp)pgno=pgno*2+1
							endif
						enddo
						block_tmp%pgno=pgno
						call ComputeParallelIndices(block_tmp,block_tmp%pgno,ptree,msh)

						if(blocks_sml%M_loc>0)then
						allocate(Vin_sml(blocks_sml%M_loc,num_vect_sub))
						allocate(Vout_sml(blocks_sml%M_loc,num_vect_sub))
						Vout_sml=0
						endif
						if(block_tmp%M_loc>0)then
							allocate(Vin_tmp1(block_tmp%M_loc,num_vect_sub))
							allocate(Vout_tmp1(block_tmp%M_loc,num_vect_sub))
							pp = ptree%myid-ptree%pgrp(block_tmp%pgno)%head+1
							idx_start_glo_tmp = block_tmp%headm + block_tmp%M_p(pp,1) -1
							pp = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
							idx_start_glo = blocks%headm + blocks%M_p(pp,1) -1
							idx_start_loc = idx_start_glo_tmp-idx_start_glo+1
							idx_end_loc = idx_start_loc + block_tmp%M_loc-1
							Vin_tmp1 = Vin(idx_start_loc:idx_end_loc,:)
						endif
						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vin_tmp1,block_tmp%M_p,0,block_tmp%pgno,Vin_sml,blocks_sml%M_p,0,blocks_sml%pgno,num_vect_sub,ptree)

						if(IOwnPgrp(ptree,blocks_sml%pgno))then
							call BF_block_MVP_dat(blocks_sml,trans,blocks_sml%M_loc,blocks_sml%N_loc,num_vect_sub,Vin_sml,Vout_sml,ctemp1,ctemp2,ptree,stats)
						endif

						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vout_sml,blocks_sml%M_p,0,blocks_sml%pgno,Vout_tmp1,block_tmp%M_p,0,block_tmp%pgno,num_vect_sub,ptree)
						if(block_tmp%M_loc>0)then
							vec_new(idx_start_loc:idx_end_loc,:) = Vout_tmp1
							deallocate(Vin_tmp1)
							deallocate(Vout_tmp1)
						endif

						if(blocks_sml%M_loc>0)then
							deallocate(Vout_sml)
							deallocate(Vin_sml)
						endif
						call BF_delete(block_tmp,1)

					endif
				enddo
				vec_new = vec_new + Vin

				call BF_block_MVP_dat(blocks,trans,M,N,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				deallocate(vec_new)
		   end if

		   Vout = a*Vout + b*Vout_tmp
		   deallocate(Vout_tmp)

    class default
		write(*,*)"unexpected type"
		stop

   end select
   end select

end subroutine Bplus_block_MVP_diagBinvB_dat







subroutine Bplus_block_MVP_BBplus_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
	use BPACK_DEFS


	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real(kind=8) a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	type(proctree)::ptree
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1
	type(Hstat)::stats

	real(kind=8)::n2,n1

   select TYPE(bplus)

   type is (blockplus)

	 select TYPE(operand1)
	   type is (matrixblock)

			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout

			level_butterfly = block_o%level_butterfly
			num_blocks=2**level_butterfly
			mm=0
			nn=0
			do i=1, num_blocks
				mm = mm+size(block_o%ButterflyU%blocks(i)%matrix,1)
				nn = nn+size(block_o%ButterflyV%blocks(i)%matrix,1)
			enddo


			allocate(vec_new(nn,num_vect_sub))
			vec_new = 0

			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0

				call Bplus_block_MVP_dat(bplus,'N',nn,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)
				vec_new = vec_new + Vin

				call BF_block_MVP_dat(operand1,'N',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				deallocate(vec_new)

		   else if(trans=='T')then

				! get the left multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call BF_block_MVP_dat(operand1,'T',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)

				call Bplus_block_MVP_dat(bplus,'T',nn,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				Vout = Vout + vec_new

				deallocate(vec_new)

		   end if

		   Vout = a*Vout + b*Vout_tmp
		   deallocate(Vout_tmp)

		class default
			write(*,*)"unexpected type"
			stop

	   end select

    class default
		write(*,*)"unexpected type"
		stop

   end select

end subroutine Bplus_block_MVP_BBplus_dat




subroutine Bplus_block_MVP_BdiagBinv_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,msh)
	use BPACK_DEFS

	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	DT :: Vin(:,:), Vout(:,:)
	DT :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,pp,test
	integer level_butterfly,groupm_diag,edge_first
	! real(kind=8) a,b,c,d
	integer idx_start_glo_sml,idx_start_glo_tmp,idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc,Maxgrp,pgno,row_group
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:),Vin_sml(:,:),Vout_sml(:,:),Vin_tmp1(:,:),Vout_tmp1(:,:)

	integer*8 idx_start
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	! type(vectorsblock), pointer :: random1, random2

	class(*):: bplus
	type(matrixblock)::block_o,block_tmp
	type(matrixblock),pointer::blocks,blocks_sml
	class(*),optional::msh
	type(proctree)::ptree
	type(Hstat)::stats
	real(kind=8)::n2,n1


   call assert(present(msh),'operand1 cannot be skipped')
   select TYPE(msh)
   type is (mesh)
   select TYPE(bplus)
   type is (blockplus)
			blocks=> bplus%LL(1)%matrices_block(1)

			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout
			Vout=0
			allocate(vec_new(mv,nv))
			vec_new = 0

			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				do ii=1,Bplus%LL(2)%Nbound
					blocks_sml=>Bplus%LL(2)%matrices_block(ii)
					edge_first = blocks_sml%headm
					if(edge_first>=blocks%headm .and. edge_first<=blocks%headm+blocks%M-1)then




						! generate a tmporay block block_tmp whose layout is conformal to blocks, then redistribute data from block_tmp to blocks_sml
						block_tmp%M=blocks_sml%M
						block_tmp%headm=blocks_sml%headm
						block_tmp%N=blocks_sml%N
						block_tmp%headn=blocks_sml%headn
						block_tmp%row_group=blocks_sml%row_group
						block_tmp%col_group=blocks_sml%col_group
						block_tmp%level=blocks_sml%level


						Maxgrp=2**(ptree%nlevel)-1

						row_group=blocks%row_group
						pgno=blocks%pgno
						do while(row_group/=block_tmp%row_group)
							if(msh%basis_group(2*row_group)%head<=msh%basis_group(block_tmp%row_group)%head .and. msh%basis_group(2*row_group)%tail>=msh%basis_group(block_tmp%row_group)%tail)then
								row_group = 2*row_group
								if(pgno*2<=Maxgrp)pgno=pgno*2
							else
								row_group = 2*row_group+1
								if(pgno*2+1<=Maxgrp)pgno=pgno*2+1
							endif
						enddo
						block_tmp%pgno=pgno
						call ComputeParallelIndices(block_tmp,block_tmp%pgno,ptree,msh)

						if(blocks_sml%M_loc>0)then
						allocate(Vin_sml(blocks_sml%M_loc,num_vect_sub))
						allocate(Vout_sml(blocks_sml%M_loc,num_vect_sub))
						Vout_sml=0
						endif
						if(block_tmp%M_loc>0)then
							allocate(Vin_tmp1(block_tmp%M_loc,num_vect_sub))
							allocate(Vout_tmp1(block_tmp%M_loc,num_vect_sub))
							pp = ptree%myid-ptree%pgrp(block_tmp%pgno)%head+1
							idx_start_glo_tmp = block_tmp%headm + block_tmp%M_p(pp,1) -1
							pp = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
							idx_start_glo = blocks%headm + blocks%M_p(pp,1) -1
							idx_start_loc = idx_start_glo_tmp-idx_start_glo+1
							idx_end_loc = idx_start_loc + block_tmp%M_loc-1
							Vin_tmp1 = Vin(idx_start_loc:idx_end_loc,:)
						endif
						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vin_tmp1,block_tmp%M_p,0,block_tmp%pgno,Vin_sml,blocks_sml%M_p,0,blocks_sml%pgno,num_vect_sub,ptree)

						if(IOwnPgrp(ptree,blocks_sml%pgno))then
							call BF_block_MVP_dat(blocks_sml,trans,blocks_sml%M_loc,blocks_sml%N_loc,num_vect_sub,Vin_sml,Vout_sml,ctemp1,ctemp2,ptree,stats)
						endif

						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vout_sml,blocks_sml%M_p,0,blocks_sml%pgno,Vout_tmp1,block_tmp%M_p,0,block_tmp%pgno,num_vect_sub,ptree)
						if(block_tmp%M_loc>0)then
							vec_new(idx_start_loc:idx_end_loc,:) = Vout_tmp1
							deallocate(Vin_tmp1)
							deallocate(Vout_tmp1)
						endif

						if(blocks_sml%M_loc>0)then
							deallocate(Vout_sml)
							deallocate(Vin_sml)
						endif
						call BF_delete(block_tmp,1)

					endif
				enddo
				vec_new = vec_new + Vin

				call BF_block_MVP_dat(blocks,trans,M,N,num_vect_sub,vec_new,Vout,ctemp1,ctemp2,ptree,stats)
				deallocate(vec_new)

		   else if(trans=='T')then
				! get the left multiplied vectors

				ctemp1=1.0d0 ; ctemp2=0.0d0
				call BF_block_MVP_dat(blocks,trans,M,N,num_vect_sub,Vin,vec_new,ctemp1,ctemp2,ptree,stats)

				do ii=1,Bplus%LL(2)%Nbound
					blocks_sml=>Bplus%LL(2)%matrices_block(ii)
					edge_first = blocks_sml%headm
					if(edge_first>=blocks%headm .and. edge_first<=blocks%headm+blocks%M-1)then


						! generate a tmporay block block_tmp whose layout is conformal to blocks, then redistribute data from block_tmp to blocks_sml
						block_tmp%M=blocks_sml%M
						block_tmp%headm=blocks_sml%headm
						block_tmp%N=blocks_sml%N
						block_tmp%headn=blocks_sml%headn
						block_tmp%row_group=blocks_sml%row_group
						block_tmp%col_group=blocks_sml%col_group
						block_tmp%level=blocks_sml%level


						Maxgrp=2**(ptree%nlevel)-1

						row_group=blocks%row_group
						pgno=blocks%pgno
						do while(row_group/=block_tmp%row_group)
							if(msh%basis_group(2*row_group)%head<=msh%basis_group(block_tmp%row_group)%head .and. msh%basis_group(2*row_group)%tail>=msh%basis_group(block_tmp%row_group)%tail)then
								row_group = 2*row_group
								if(pgno*2<=Maxgrp)pgno=pgno*2
							else
								row_group = 2*row_group+1
								if(pgno*2+1<=Maxgrp)pgno=pgno*2+1
							endif
						enddo
						block_tmp%pgno=pgno
						call ComputeParallelIndices(block_tmp,block_tmp%pgno,ptree,msh)

						if(blocks_sml%M_loc>0)then
						allocate(Vin_sml(blocks_sml%M_loc,num_vect_sub))
						allocate(Vout_sml(blocks_sml%M_loc,num_vect_sub))
						Vout_sml=0
						endif
						if(block_tmp%M_loc>0)then
							allocate(Vin_tmp1(block_tmp%M_loc,num_vect_sub))
							allocate(Vout_tmp1(block_tmp%M_loc,num_vect_sub))
							pp = ptree%myid-ptree%pgrp(block_tmp%pgno)%head+1
							idx_start_glo_tmp = block_tmp%headm + block_tmp%M_p(pp,1) -1
							pp = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
							idx_start_glo = blocks%headm + blocks%M_p(pp,1) -1
							idx_start_loc = idx_start_glo_tmp-idx_start_glo+1
							idx_end_loc = idx_start_loc + block_tmp%M_loc-1
							Vin_tmp1 = vec_new(idx_start_loc:idx_end_loc,:)
						endif
						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vin_tmp1,block_tmp%M_p,0,block_tmp%pgno,Vin_sml,blocks_sml%M_p,0,blocks_sml%pgno,num_vect_sub,ptree)



						if(IOwnPgrp(ptree,blocks_sml%pgno))then
							call BF_block_MVP_dat(blocks_sml,trans,blocks_sml%M_loc,blocks_sml%N_loc,num_vect_sub,Vin_sml,Vout_sml,ctemp1,ctemp2,ptree,stats)
						endif

						if(blocks_sml%M_loc>0 .or. block_tmp%M_loc>0)	call Redistribute1Dto1D(Vout_sml,blocks_sml%M_p,0,blocks_sml%pgno,Vout_tmp1,block_tmp%M_p,0,block_tmp%pgno,num_vect_sub,ptree)
						if(block_tmp%M_loc>0)then
							Vout(idx_start_loc:idx_end_loc,:) = Vout_tmp1
							deallocate(Vin_tmp1)
							deallocate(Vout_tmp1)
						endif

						if(blocks_sml%M_loc>0)then
							deallocate(Vout_sml)
							deallocate(Vin_sml)
						endif
						call BF_delete(block_tmp,1)


					endif
				enddo
				Vout = Vout + vec_new

				deallocate(vec_new)

		   end if

		   Vout = a*Vout + b*Vout_tmp
		   deallocate(Vout_tmp)

    class default
		write(*,*)"unexpected type"
		stop

   end select
   end select
end subroutine Bplus_block_MVP_BdiagBinv_dat





subroutine Bplus_MultiLrandomized_Onesubblock(rank0,rankrate,blocks,operand,blackbox_MVP_dat,error_inout,strings,option,stats,ptree,msh,operand1)

   use BPACK_DEFS


   use MISC_Utilities
   implicit none

    type(blockplus),pointer::bplus
	integer:: ii,ll,bb,jj,bb_o,tt,rank0
    real(kind=8) Memory,rtemp,error_inout,n2,n1,mem_vec,rankrate
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	DT,allocatable::Vout1(:,:),Vout2(:,:),Vout3(:,:),Vin(:,:)
	integer M,N,idx_start_n,idx_start_m,idx_start_n_loc,idx_end_n_loc,idx_start_m_loc,idx_end_m_loc,mm,nn,rmax,rank,idx_start_n_ref,idx_start_m_ref,idx_end_n_ref,idx_end_m_ref,head,tail
	DT::ctemp1,ctemp2,Ctemp
	type(matrixblock)::blocks
	type(matrixblock)::block_dummy
	DT, allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),matrix_small(:,:)
	real(kind=8), allocatable :: Singular(:)
	integer level_c,rowblock,Nactive
	integer,allocatable::boxindex(:)
	integer Chunksize, Nchunk, Nidx, idx_s,cc
	class(*):: operand
	class(*),optional:: operand1
	character(*)  :: strings
	type(Hoption)::option
	type(Hstat)::stats
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	DT, allocatable :: matU_glo(:,:), matV_glo(:,:)
	procedure(BMatVec)::blackbox_MVP_dat
	type(proctree)::ptree
	type(mesh)::msh
	real(kind=8)flop
	type(matrixblock),pointer::block_large

	select TYPE(operand1)
	type is (blockplus)
		block_large=>operand1%LL(1)%matrices_block(1)

		stats%Flop_tmp=0
		option%tol_rand=option%tol_rand*1d-1;
		call BF_randomized(block_large%pgno,blocks%level_butterfly,rank0,rankrate,blocks,operand,Bplus_block_MVP_Onesubblock_dat,error_inout,'Sub',option,stats,ptree,msh,operand1)
		option%tol_rand=option%tol_rand*1d1
		stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
		stats%Flop_Tmp=0
	end select

contains

	subroutine Bplus_block_MVP_Onesubblock_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
		use BPACK_DEFS


		implicit none
		integer level_c, rowblock, num_vect_sub
		character trans
		DT :: Vin(:,:), Vout(:,:)
		DT,allocatable :: Vin_large(:,:),Vout_large(:,:)
		DT :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
		integer M,N,mv,nv,M_loc_large,N_loc_large
		type(Hstat)::stats

		class(*):: operand
		type(matrixblock)::block_o
		type(matrixblock),pointer::block_large
		class(*),optional::operand1
		type(proctree)::ptree
		real(kind=8)::n2,n1

		integer,pointer::Nin_p(:,:),Nin_pl(:,:),Nout_p(:,:),Nout_pl(:,:)
		integer nin,ninl,nout,noutl,offin,offinl,offout,offoutl

		call assert(present(operand1),'operand1 cannot be skipped')

		select TYPE(operand1)
		type is (blockplus)
			block_large=>operand1%LL(1)%matrices_block(1)

			if(trans=='N')then
				nin = block_o%N_loc
				ninl = block_large%N_loc
				Nin_p => block_o%N_p
				Nin_pl => block_large%N_p
				offin = block_o%headn
				offinl = block_large%headn
				noutl = block_large%M_loc
				nout = block_o%M_loc
				Nout_pl => block_large%M_p
				Nout_p => block_o%M_p
				offoutl = block_large%headm
				offout = block_o%headm
			else
				nin = block_o%M_loc
				ninl = block_large%M_loc
				Nin_p => block_o%M_p
				Nin_pl => block_large%M_p
				offin = block_o%headm
				offinl = block_large%headm
				noutl = block_large%N_loc
				nout = block_o%N_loc
				Nout_pl => block_large%N_p
				Nout_p => block_o%N_p
				offoutl = block_large%headn
				offout = block_o%headn
			endif

			allocate(Vin_large(ninl,num_vect_sub))
			Vin_large=0
			allocate(Vout_large(noutl,num_vect_sub))
			Vout_large=0

			call Redistribute1Dto1D(Vin,Nin_p,offin,block_o%pgno,Vin_large,Nin_pl,offinl,block_large%pgno,num_vect_sub,ptree)
			call Redistribute1Dto1D(Vout,Nout_p,offout,block_o%pgno,Vout_large,Nout_pl,offoutl,block_large%pgno,num_vect_sub,ptree)

			call blackbox_MVP_dat(operand,block_large,trans,block_large%M_loc,block_large%N_loc,num_vect_sub,Vin_large,Vout_large,a,b,ptree,stats,msh)

			call Redistribute1Dto1D(Vout_large,Nout_pl,offoutl,block_large%pgno,Vout,Nout_p,offout,block_o%pgno,num_vect_sub,ptree)

		end select

	end subroutine Bplus_block_MVP_Onesubblock_dat


end subroutine Bplus_MultiLrandomized_Onesubblock



subroutine Bplus_randomized_constr(level_butterfly,bplus_o,operand,rank0_inner,rankrate_inner,blackbox_MVP_dat_inner,rank0_outter,rankrate_outter,blackbox_MVP_dat_outter,error_inout,strings,option,stats,ptree,msh)

   use BPACK_DEFS


   use MISC_Utilities
   implicit none

    ! type(blockplus),pointer::bplus
    type(blockplus)::bplus_o
	integer:: ii,ll,bb
    real(kind=8) rtemp,error,Memory,n2,n1,rate,error_inout,err_avr,rankrate_inner,rankrate_outter
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall,M,N,err_cnt
	DT,allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	DT ctemp, ctemp1, ctemp2
	integer level_c,rowblock,rank_new_max,rank0_inner,rank0_outter
	type(matrixblock),pointer::block_off1,block_off2,block_o
	class(*)::operand
	character(*)  :: strings
	integer rankthusfar
	type(Hoption)::option
	type(Hstat)::stats
	type(blockplus) :: Bplus_randomized
	procedure(BMatVec)::blackbox_MVP_dat_inner,blackbox_MVP_dat_outter
	type(proctree)::ptree
	type(mesh)::msh
	integer ierr,rank0

	error_inout=0
	! Memory = 0
	call assert(bplus_o%Lplus>=2,'this is not a multi Bplus in Bplus_randomized')

	call Bplus_Init_FromInput(bplus_o,Bplus_randomized,msh,ptree,option)

	! write(*,*)'hhhhh1',Bplus_randomized%LL(1)%matrices_block(1)%level
	n1 = OMP_get_wtime()

	rankthusfar = 0
	do bb =1,Bplus_randomized%LL(2)%Nbound
		! ! write(*,*)bb,Bplus_randomized%LL(2)%Nbound,'dddd'

		! rank0_inner = ho_bf1%levels(level_c)%BP(2*rowblock-1)%LL(2)%rankmax
		! rankrate_inner = 2.0d0
		block_o => Bplus_randomized%LL(2)%matrices_block(bb)
		! ho_bf1%ind_lv = level_c
		! ho_bf1%ind_bk = rowblock
		rank0 = min(max(rank0_inner+3,rankthusfar),min(block_o%M,block_o%N))
		call Bplus_MultiLrandomized_Onesubblock(rank0,rankrate_inner,block_o,operand,blackbox_MVP_dat_inner,error,strings,option,stats,ptree,msh,Bplus_randomized)
		error_inout = max(error_inout, error)
		! write(*,*)'go'

		call MPI_ALLREDUCE(MPI_IN_PLACE,rankthusfar,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(Bplus_randomized%LL(1)%matrices_block(1)%pgno)%Comm,ierr)

	end do
	! call Test_Error_RR_Inner_Exact(bplus_o)
	n2 = OMP_get_wtime()
	stats%Time_random(4) = stats%Time_random(4) + n2-n1


	! level_c = ho_bf1%ind_lv
	! rowblock = ho_bf1%ind_bk

	! block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	! ! block_off1 => ho_bf1%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)
	! block_off2 => ho_bf1%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)


	! level_butterfly=int((maxlevel_for_blocks-bplus_o%level)/2)*2
	! rank0_outter = max(block_off1%rankmax,block_off2%rankmax)
	! rate_outter=1.2d0
	call BF_randomized(Bplus_randomized%LL(1)%matrices_block(1)%pgno,level_butterfly,rank0_outter,rankrate_outter,Bplus_randomized%LL(1)%matrices_block(1),operand,blackbox_MVP_dat_outter,error,'Outter',option,stats,ptree,msh,Bplus_randomized)
	error_inout = max(error_inout, error)
	stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp

	do ll=1,Bplus_randomized%Lplus
	Bplus_randomized%LL(ll)%rankmax=0
	do bb=1,Bplus_randomized%LL(ll)%Nbound
		Bplus_randomized%LL(ll)%rankmax=max(Bplus_randomized%LL(ll)%rankmax,Bplus_randomized%LL(ll)%matrices_block(bb)%rankmax)
	enddo
	call MPI_ALLREDUCE(MPI_IN_PLACE,Bplus_randomized%LL(ll)%rankmax,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(Bplus_randomized%LL(1)%matrices_block(1)%pgno)%Comm,ierr)
	end do


	call Bplus_delete(bplus_o)
	call Bplus_copy_delete(Bplus_randomized,bplus_o,Memory)
	! deallocate(Bplus_randomized)


	rank_new_max = 0
	do ll=1,bplus_o%Lplus
		if(ll==1 .or. option%bp_cnt_lr==1)rank_new_max = max(rank_new_max,bplus_o%LL(ll)%rankmax)
	end do



	if(option%verbosity>=2 .and. ptree%myid==ptree%pgrp(bplus_o%LL(1)%matrices_block(1)%pgno)%head)write(*,'(A20,A8,I3,A8,I3,A11,Es14.7)')strings,' rank:',rank_new_max,' L_butt:',bplus_o%LL(1)%matrices_block(1)%level_butterfly,' error:',error_inout
    return

end subroutine Bplus_randomized_constr




subroutine Bplus_Init_FromInput(Bplus,Bplus_randomized,msh,ptree,option)
	use MISC_Utilities
    use BPACK_DEFS


    implicit none

    integer level_c,rowblock
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real(kind=8) a,b,c,d
    DT ctemp
	DT, allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real(kind=8), allocatable:: Singular(:)
	integer mn_min,index_i,index_j,blocks_idx
	type(matrixblock)::block
	type(blockplus)::Bplus,Bplus_randomized
	integer seed_myid(2)
	integer time(8)
	integer ll,bb
	integer:: level_BP,levelm,groupm_start,Nboundall
	type(mesh)::msh
	type(proctree)::ptree
	type(Hoption)::option


	! allocate(Bplus_randomized)

	Bplus_randomized%level = Bplus%level
	Bplus_randomized%col_group = Bplus%col_group
	Bplus_randomized%row_group = Bplus%row_group
	Bplus_randomized%Lplus = Bplus%Lplus
	Bplus_randomized%pgno = Bplus%pgno
	Bplus_randomized%boundary = Bplus%boundary

	allocate(Bplus_randomized%LL(LplusMax))


	! Bplus_randomized%LL(1)%Nbound = 1
	! allocate(Bplus_randomized%LL(1)%matrices_block(1))
	! Bplus_randomized%LL(1)%matrices_block(1)%level = Bplus_randomized%level
	! Bplus_randomized%LL(1)%matrices_block(1)%col_group = Bplus_randomized%col_group
	! Bplus_randomized%LL(1)%matrices_block(1)%row_group = Bplus_randomized%row_group
	! Bplus_randomized%LL(1)%matrices_block(1)%style = Bplus%LL(1)%matrices_block(1)%style
	allocate(Bplus_randomized%LL(1)%boundary_map(1))
	Bplus_randomized%LL(1)%boundary_map(1) = Bplus%LL(1)%boundary_map(1)

	do ll=1,LplusMax
	Bplus_randomized%LL(ll)%Nbound = 0
	end do

	do ll=1,LplusMax-1
		if(Bplus%LL(ll)%Nbound>0)then

			Bplus_randomized%LL(ll)%rankmax = Bplus%LL(ll)%rankmax
			Bplus_randomized%LL(ll)%Nbound = Bplus%LL(ll)%Nbound
			allocate(Bplus_randomized%LL(ll)%matrices_block(Bplus_randomized%LL(ll)%Nbound))

			do bb =1,Bplus_randomized%LL(ll)%Nbound
				call BF_Init_randomized(Bplus%LL(ll)%matrices_block(bb)%level_butterfly,Bplus%LL(ll)%rankmax,Bplus%LL(ll)%matrices_block(bb)%row_group,Bplus%LL(ll)%matrices_block(bb)%col_group,Bplus%LL(ll)%matrices_block(bb),Bplus_randomized%LL(ll)%matrices_block(bb),msh,ptree,option,1)
			end do


			if(Bplus%LL(ll+1)%Nbound==0)then
				Bplus_randomized%LL(ll+1)%Nbound=0
			else
				level_butterfly = Bplus%LL(ll)%matrices_block(1)%level_butterfly
				! level_butterfly = 0 !!! only lowrank
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)
				groupm_start=Bplus%LL(ll)%matrices_block(1)%row_group*2**levelm
				! Nboundall = 2**(Bplus%LL(ll)%matrices_block(1)%level+levelm-level_BP)
				Nboundall=size(Bplus%LL(ll+1)%boundary_map)

				allocate(Bplus_randomized%LL(ll+1)%boundary_map(Nboundall))
				! write(*,*)shape(Bplus%LL(ll+1)%boundary_map),shape(Bplus_randomized%LL(ll+1)%boundary_map),'didi',ll

				! write(*,*)'gali',ll,Nboundall,shape(Bplus%LL(ll+1)%boundary_map)
				Bplus_randomized%LL(ll+1)%boundary_map = Bplus%LL(ll+1)%boundary_map
			end if
		else
			exit
		end if

	end do

    return

end subroutine Bplus_Init_FromInput




end module Bplus_randomizedop
