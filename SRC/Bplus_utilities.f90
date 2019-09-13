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
module Bplus_Utilities
use MISC_Utilities
contains


subroutine Bplus_delete(bplus)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block
type(blockplus)::bplus

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8)::rtemp

if(associated(bplus%LL))then
do ll=1,LplusMax
	if(bplus%LL(ll)%Nbound>0)then
		if(associated(bplus%LL(ll)%matrices_block))then
		do bb=1,bplus%LL(ll)%Nbound
			! write(*,*)ll,bplus%Lplus,bb,bplus%LL(ll)%Nbound,'fff'
			call BF_delete(bplus%LL(ll)%matrices_block(bb),1)
		end do
		deallocate(bplus%LL(ll)%matrices_block)
		endif
		if(allocated(bplus%LL(ll)%boundary_map))deallocate(bplus%LL(ll)%boundary_map)
	end if
end do
deallocate(bplus%LL)
endif

end subroutine Bplus_delete


subroutine Bplus_copy(bplus_i,bplus_o,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8),optional::memory
real(kind=8)::rtemp

call Bplus_delete(bplus_o)

if(present(memory))memory=0

allocate(bplus_o%LL(LplusMax))
bplus_o%Lplus = bplus_i%Lplus
bplus_o%boundary = bplus_i%boundary
bplus_o%level = bplus_i%level
bplus_o%col_group = bplus_i%col_group
bplus_o%row_group = bplus_i%row_group
bplus_o%pgno = bplus_i%pgno


do ll=1,LplusMax
	bplus_o%LL(ll)%Nbound=bplus_i%LL(ll)%Nbound
	bplus_o%LL(ll)%rankmax=bplus_i%LL(ll)%rankmax



	if(bplus_i%LL(ll)%Nbound>0)then
		allocate(bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
		do bb=1,bplus_i%LL(ll)%Nbound
			call BF_copy('N',bplus_i%LL(ll)%matrices_block(bb),bplus_o%LL(ll)%matrices_block(bb),rtemp)
			if(present(memory))memory=memory+rtemp
		end do
		if(allocated(bplus_i%LL(ll)%boundary_map))then
			Nboundall=size(bplus_i%LL(ll)%boundary_map)
			allocate(bplus_o%LL(ll)%boundary_map(Nboundall))
			if(present(memory))memory=memory+ SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
			bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
		endif
	end if
end do

end subroutine Bplus_copy


subroutine Bplus_copy_delete(bplus_i,bplus_o,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8),optional::memory
real(kind=8)::rtemp

if(present(memory))memory=0

allocate(bplus_o%LL(LplusMax))
bplus_o%Lplus = bplus_i%Lplus
bplus_o%boundary = bplus_i%boundary
bplus_o%level = bplus_i%level
bplus_o%col_group = bplus_i%col_group
bplus_o%row_group = bplus_i%row_group


do ll=1,LplusMax
	bplus_o%LL(ll)%Nbound=bplus_i%LL(ll)%Nbound
	bplus_o%LL(ll)%rankmax=bplus_i%LL(ll)%rankmax
	if(bplus_i%LL(ll)%Nbound>0)then
		allocate(bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
		do bb=1,bplus_i%LL(ll)%Nbound
			call BF_copy_delete(bplus_i%LL(ll)%matrices_block(bb),bplus_o%LL(ll)%matrices_block(bb),rtemp)
			if(present(memory))memory=memory+rtemp
		end do
		deallocate(bplus_i%LL(ll)%matrices_block)
		Nboundall=size(bplus_i%LL(ll)%boundary_map)
		allocate(bplus_o%LL(ll)%boundary_map(Nboundall))
		if(present(memory))memory=memory+ SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
		bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
		deallocate(bplus_i%LL(ll)%boundary_map)
	end if
end do

deallocate(bplus_i%LL)
end subroutine Bplus_copy_delete



subroutine Bplus_extract_partial(bplus_i,ll_s,row_group,agent_bplus,msh)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,agent_bplus

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb,bb_o
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8)::rtemp
integer row_group,ll_s,idx_s,idx_e
type(mesh)::msh

call assert(bplus_i%row_group==bplus_i%col_group,'only works for square matrix')

idx_s = msh%basis_group(row_group)%head
idx_e = msh%basis_group(row_group)%tail

! allocate(agent_bplus)
allocate(agent_bplus%LL(LplusMax))
do ll=1,LplusMax
agent_bplus%LL(ll)%Nbound = 0
end do

agent_bplus%Lplus = bplus_i%Lplus - ll_s + 1
agent_bplus%row_group = 	row_group
agent_bplus%col_group = 	row_group
agent_bplus%level = GetTreelevel(row_group)-1

do ll=1,agent_bplus%Lplus
	agent_bplus%LL(ll)%Nbound = 0
	agent_bplus%LL(ll)%rankmax=bplus_i%LL(ll+ll_s-1)%rankmax
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(msh%basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. msh%basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			agent_bplus%LL(ll)%Nbound = agent_bplus%LL(ll)%Nbound + 1
		end if
	end do
	if(agent_bplus%LL(ll)%Nbound>0)then
		allocate(agent_bplus%LL(ll)%matrices_block(agent_bplus%LL(ll)%Nbound))
	end if
end do


do ll=1,agent_bplus%Lplus
	bb_o = 0
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(msh%basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. msh%basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			bb_o = bb_o + 1
			call BF_copy('N',bplus_i%LL(ll+ll_s-1)%matrices_block(bb),agent_bplus%LL(ll)%matrices_block(bb_o))
		end if
	end do
end do



end subroutine Bplus_extract_partial



subroutine Bplus_ComputeMemory(bplus_i,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8)::memory
real(kind=8)::rtemp

memory=0

do ll=1,LplusMax
	if(bplus_i%LL(ll)%Nbound>0)then
		do bb=1,bplus_i%LL(ll)%Nbound
			call BF_ComputeMemory(bplus_i%LL(ll)%matrices_block(bb),rtemp)
			memory=memory+rtemp
		end do
	end if
end do

end subroutine Bplus_ComputeMemory




logical function Bplus_checkNAN(bplus_i)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall
real(kind=8)::rtemp
Bplus_checkNAN = .false.

do ll=1,LplusMax
	if(bplus_i%LL(ll)%Nbound>0)then
		do bb=1,bplus_i%LL(ll)%Nbound
			if(BF_checkNAN(bplus_i%LL(ll)%matrices_block(bb)))then
				Bplus_checkNAN = .true.
				return
			end if
		end do
	end if
end do

end function Bplus_checkNAN

subroutine Bplus_block_MVP_dat(bplus,chara,M,N,Nrnd,random1,random2,a,b,ptree,stats,level_start,level_end)

    use BPACK_DEFS
	use MISC_Utilities
    implicit none

    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    DT ctemp, a, b,ctemp1,ctemp2
    character chara
	type(matrixblock),pointer::blocks,blocks_1
	type(blockplus)::bplus
	integer ll,bb
	integer,optional:: level_start,level_end
	integer:: level_s,level_e
	type(proctree)::ptree
	type(Hstat)::stats

    DT :: random1(:,:), random2(:,:)
    DT,allocatable :: Vout(:,:),Vin_loc(:,:),Vout_loc(:,:)
	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:)

	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)

	integer idx_start_m,idx_start_n, idx_start_n_loc,idx_start_m_loc, idx_end_n_loc,idx_end_m_loc,idx_start_i_loc,idx_start_o_loc,idx_end_i_loc,idx_end_o_loc

	level_s=1
	level_e=bplus%Lplus
	if(present(level_start))level_s=level_start
	if(present(level_end))level_e=level_end


	if (chara=='N')allocate(Vout(M,Nrnd))
	if (chara=='T')allocate(Vout(N,Nrnd))
	Vout = 0


	ctemp1=1.0d0 ; ctemp2=1.0d0

	blocks_1 => bplus%LL(1)%matrices_block(1)

	do ll=level_s,level_e
		do bb = 1,bplus%LL(ll)%Nbound
			blocks => bplus%LL(ll)%matrices_block(bb)

			if (chara=='N')then

				if(blocks%M_loc>0)allocate(Vout_loc(blocks%M_loc,Nrnd))
				if(blocks%N_loc>0)allocate(Vin_loc(blocks%N_loc,Nrnd))
				call Redistribute1Dto1D(random1,blocks_1%N_p,blocks_1%headn,blocks_1%pgno,Vin_loc,blocks%N_p,blocks%headn,blocks%pgno,Nrnd,ptree)
				call Redistribute1Dto1D(Vout,blocks_1%M_p,blocks_1%headm,blocks_1%pgno,Vout_loc,blocks%M_p,blocks%headm,blocks%pgno,Nrnd,ptree)
			else

				if(blocks%N_loc>0)allocate(Vout_loc(blocks%N_loc,Nrnd))
				if(blocks%M_loc>0)allocate(Vin_loc(blocks%M_loc,Nrnd))
				call Redistribute1Dto1D(random1,blocks_1%M_p,blocks_1%headm,blocks_1%pgno,Vin_loc,blocks%M_p,blocks%headm,blocks%pgno,Nrnd,ptree)
				call Redistribute1Dto1D(Vout,blocks_1%N_p,blocks_1%headn,blocks_1%pgno,Vout_loc,blocks%N_p,blocks%headn,blocks%pgno,Nrnd,ptree)

			endif

			if(blocks%N_loc>0 .or. blocks%M_loc>0)then
				if(blocks%style==1)then
					write(*,*)'style 1 not implemented'
					stop
				else
					! write(*,*)'ddd1',ll,bb
					call BF_block_MVP_dat(blocks,chara,blocks%M_loc,blocks%N_loc,Nrnd,&
					&Vin_loc,Vout_loc,ctemp1,ctemp2,ptree,stats)
					! write(*,*)'ddd2'
				end if
			endif

			if (chara=='N')then
				call Redistribute1Dto1D(Vout_loc,blocks%M_p,blocks%headm,blocks%pgno,Vout,blocks_1%M_p,blocks_1%headm,blocks_1%pgno,Nrnd,ptree)
				if(blocks%M_loc>0)deallocate(Vout_loc)
				if(blocks%N_loc>0)deallocate(Vin_loc)
			else
				call Redistribute1Dto1D(Vout_loc,blocks%N_p,blocks%headn,blocks%pgno,Vout,blocks_1%N_p,blocks_1%headn,blocks_1%pgno,Nrnd,ptree)
				if(blocks%N_loc>0)deallocate(Vout_loc)
				if(blocks%M_loc>0)deallocate(Vin_loc)
			endif

		end do
	end do

	random2 = random2*b + Vout*a
	deallocate(Vout)

end subroutine Bplus_block_MVP_dat



! redistribute Bplus
subroutine Bplus_DoubleDistribute(bplus_o,stats,ptree)
implicit none

integer nproc_i, nproc_o,idxs_i,idxs_o,idxe_i,idxe_o,ii,jj,iii,jjj,level
type(proctree)::ptree
integer,allocatable::S_req(:),R_req(:)
integer,allocatable:: statuss(:,:),statusr(:,:)
integer tag,Nreqs,Nreqr,recvid,sendid,ierr,head_i,head_o,rank,rankmax
type(blockplus)::bplus_o
type(matrixblock),pointer::blocks
DT,pointer::dat_new(:,:),dat_old(:,:)
real(kind=8)::n1,n2
type(Hstat)::stats
integer pgno,ll,bb


if(associated(bplus_o%LL))then
do ll=1,LplusMax
	if(bplus_o%LL(ll)%Nbound>0)then
		if(associated(bplus_o%LL(ll)%matrices_block))then
		do bb=1,bplus_o%LL(ll)%Nbound
			pgno = bplus_o%LL(ll)%matrices_block(bb)%pgno_db
			if(IOwnPgrp(ptree,pgno))call BF_DoubleDistribute(bplus_o%LL(ll)%matrices_block(bb),stats,ptree)
		end do
		endif
		call MPI_ALLREDUCE(MPI_IN_PLACE,bplus_o%LL(ll)%rankmax,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(bplus_o%LL(1)%matrices_block(1)%pgno)%Comm,ierr)
	end if
end do
endif

end subroutine Bplus_DoubleDistribute



! redistribute Butterfly
subroutine BF_DoubleDistribute(blocks,stats,ptree)
implicit none

integer nproc_i, nproc_o,idxs_i,idxs_o,idxe_i,idxe_o,ii,jj,iii,jjj,level
type(proctree)::ptree
integer,allocatable::S_req(:),R_req(:)
integer,allocatable:: statuss(:,:),statusr(:,:)
integer tag,Nreqs,Nreqr,recvid,sendid,ierr,head_i,head_o,rank,rankmax
type(matrixblock)::blocks
DT,pointer::dat_new(:,:),dat_old(:,:)
real(kind=8)::n1,n2
type(Hstat)::stats

dat_new=>null()
dat_old=>null()


	! call MPI_barrier(ptree%pgrp(blocks%pgno_db)%Comm,ierr)
	n1 = OMP_get_wtime()

	if(blocks%level_butterfly==0)then

		if(blocks%pgno/=blocks%pgno_db)then
			! communicate block sizes first
			if(blocks%M_loc>0)then
				rank = blocks%rankmax
			else
				rank = 0
			endif
			call assert(MPI_COMM_NULL/=ptree%pgrp(blocks%pgno_db)%Comm,'communicator should not be null 4')
			rankmax=0
			call MPI_ALLREDUCE(rank,rankmax,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(blocks%pgno_db)%Comm,ierr)
			rank=rankmax

			! redistribute U
			if(blocks%M_loc>0)then
				allocate(dat_old(blocks%M_loc,rank))
				dat_old=blocks%ButterflyU%blocks(1)%matrix
			endif
			if(blocks%M_loc_db>0)then
				allocate(dat_new(blocks%M_loc_db,rank))
				dat_new=0
			endif
			call Redistribute1Dto1D(dat_old,blocks%M_p,0,blocks%pgno,dat_new,blocks%M_p_db,0,blocks%pgno_db,rank,ptree)
			if(blocks%M_loc>0)then
				deallocate(blocks%ButterflyU%blocks(1)%matrix)
				deallocate(dat_old)
				deallocate(blocks%M_p)
			endif
			if(blocks%M_loc_db>0)then
				if(.not.allocated(blocks%ButterflyU%blocks))allocate(blocks%ButterflyU%blocks(1))
				if(.not.allocated(blocks%ButterflyU%blocks(1)%matrix))allocate(blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc_db,rank))
				blocks%ButterflyU%blocks(1)%matrix = dat_new
				deallocate(dat_new)
				blocks%M_loc =blocks%M_loc_db
				blocks%M_p=>blocks%M_p_db
				blocks%M_p_db=>NULL()
			endif



			! redistribute V
			if(blocks%N_loc>0)then
				allocate(dat_old(blocks%N_loc,rank))
				dat_old=blocks%ButterflyV%blocks(1)%matrix
			endif
			if(blocks%N_loc_db>0)then
				allocate(dat_new(blocks%N_loc_db,rank))
				dat_new=0
			endif
			call Redistribute1Dto1D(dat_old,blocks%N_p,0,blocks%pgno,dat_new,blocks%N_p_db,0,blocks%pgno_db,rank,ptree)
			if(blocks%N_loc>0)then
				deallocate(blocks%ButterflyV%blocks(1)%matrix)
				deallocate(dat_old)
				deallocate(blocks%N_p)
			endif


			if(blocks%N_loc_db>0)then
				if(.not.allocated(blocks%ButterflyV%blocks))allocate(blocks%ButterflyV%blocks(1))
				if(.not.allocated(blocks%ButterflyV%blocks(1)%matrix))allocate(blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc_db,rank))
				blocks%ButterflyV%blocks(1)%matrix = dat_new
				deallocate(dat_new)
				blocks%N_loc =blocks%N_loc_db
				blocks%N_p=>blocks%N_p_db
				blocks%N_p_db=>NULL()
				blocks%rankmax = rank
				blocks%pgno=blocks%pgno_db
			endif
		endif
	else
		if(blocks%pgno/=blocks%pgno_db)then

			!*** make sure every process has blocks%ButterflyKerl allocated
			if(.not. allocated(blocks%ButterflyKerl))then
				allocate(blocks%ButterflyKerl(blocks%level_butterfly))
			endif

			do level=0,blocks%level_butterfly+1
				if(level==0)then
					call BF_all2all_UV(blocks,blocks%pgno,blocks%ButterflyV,level,blocks,blocks%pgno_db,blocks%ButterflyV,level,stats,ptree)
				elseif(level==blocks%level_butterfly+1)then
					call BF_all2all_UV(blocks,blocks%pgno,blocks%ButterflyU,level,blocks,blocks%pgno_db,blocks%ButterflyU,level,stats,ptree)
				else
					call BF_all2all_ker(blocks,blocks%pgno,blocks%ButterflyKerl(level),level,blocks,blocks%pgno_db,blocks%ButterflyKerl(level),level,stats,ptree)
				endif
			enddo


			!*** delete dummy blocks%ButterflyKerl if I don't share the output butterfly
			if(.not. IOwnPgrp(ptree,blocks%pgno_db))then
				deallocate(blocks%ButterflyKerl)
			endif


			if(blocks%M_loc>0)then
				deallocate(blocks%M_p)
			endif
			if(blocks%M_loc_db>0)then
				blocks%M_loc =blocks%M_loc_db
				blocks%M_p=>blocks%M_p_db
				blocks%M_p_db=>NULL()
				blocks%pgno=blocks%pgno_db
			endif

			if(blocks%N_loc>0)then
				deallocate(blocks%N_p)
			endif
			if(blocks%N_loc_db>0)then
				blocks%N_loc =blocks%N_loc_db
				blocks%N_p=>blocks%N_p_db
				blocks%N_p_db=>NULL()
				blocks%pgno=blocks%pgno_db
			endif

			call BF_get_rank(blocks,ptree)

		endif
	endif

	n2 = OMP_get_wtime()
	stats%Time_RedistB=stats%Time_RedistB + n2-n1



end subroutine BF_DoubleDistribute



subroutine BF_delete(blocks,allflag)

    use BPACK_DEFS
    implicit none

    integer butterflyB_inuse, num_col, num_row
    integer i, j, mm, nn, rank, num_blocks, level, level_butterfly,index_i_m,index_j_m,levelm
    real(kind=8) memory_butterfly, rtemp
    type(matrixblock)::blocks
	integer allflag

        level_butterfly=blocks%level_butterfly

		if(allocated(blocks%ButterflyU%blocks))then
        ! !$omp parallel do default(shared) private(i)
        do i=1, blocks%ButterflyU%nblk_loc
            if(allocated(blocks%ButterflyU%blocks(i)%matrix))deallocate (blocks%ButterflyU%blocks(i)%matrix)
        enddo
        ! !$omp end parallel do
        deallocate (blocks%ButterflyU%blocks)
        end if

		if(allocated(blocks%ButterflyV%blocks))then
        ! !$omp parallel do default(shared) private(i)
        do i=1, blocks%ButterflyV%nblk_loc
            if(allocated(blocks%ButterflyV%blocks(i)%matrix))deallocate (blocks%ButterflyV%blocks(i)%matrix)
        enddo
        ! !$omp end parallel do
        deallocate (blocks%ButterflyV%blocks)
        end if


		if(allocated(blocks%ButterflyKerl))then
        if (level_butterfly/=0) then
            ! !$omp parallel do default(shared) private(level,i,j,num_col,num_row)
            do level=1, level_butterfly
                if(allocated(blocks%ButterflyKerl(level)%blocks))then
                do j=1, blocks%ButterflyKerl(level)%nc
                    do i=1, blocks%ButterflyKerl(level)%nr
                        if(allocated(blocks%ButterflyKerl(level)%blocks(i,j)%matrix))deallocate (blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
                    enddo
                enddo
                deallocate (blocks%ButterflyKerl(level)%blocks)
				endif
            enddo
            ! !$omp end parallel do
            deallocate (blocks%ButterflyKerl)
        endif
		end if

		if(allocated(blocks%ButterflyMiddle))then
		write(*,*)'warning: this part has not been parallelized in BF_delete'
			levelm = ceiling_safe(dble(level_butterfly)/2d0)
			do index_i_m=1, 2**levelm
				do index_j_m=1, 2**(level_butterfly-levelm)
					if(allocated(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix))deallocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix)
				end do
			end do
			deallocate(blocks%ButterflyMiddle)
		end if

        ! blocks%level_butterfly=0
		blocks%rankmax = -1000
		blocks%rankmin = 1000

		if(allocated(blocks%fullmat))deallocate (blocks%fullmat)
		if(allocated(blocks%fullmat_MPI))deallocate (blocks%fullmat_MPI)
		if(allocated(blocks%ipiv))deallocate (blocks%ipiv)
        if (allocated(blocks%Butterfly_data_MPI))deallocate (blocks%Butterfly_data_MPI)
        if (allocated(blocks%Butterfly_index_MPI))deallocate (blocks%Butterfly_index_MPI)

		if(allflag==1)then
			if(associated(blocks%N_p))deallocate(blocks%N_p)
			if(associated(blocks%M_p))deallocate(blocks%M_p)
			if(associated(blocks%M_p_db))deallocate(blocks%M_p_db)
			if(associated(blocks%N_p_db))deallocate(blocks%N_p_db)
		endif
    return

end subroutine BF_delete


subroutine BF_copy(trans,block_i,block_o,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i,block_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
character::trans
real(kind=8),optional::memory

if(present(memory))memory=0

if(trans=='N')then

	block_o%level = block_i%level
	block_o%col_group = block_i%col_group
	block_o%row_group = block_i%row_group
	block_o%style = block_i%style
	block_o%level_butterfly = block_i%level_butterfly
	block_o%level_half = block_i%level_half
	block_o%rankmax = block_i%rankmax
	block_o%rankmin = block_i%rankmin
	block_o%dimension_rank = block_i%dimension_rank
	block_o%M = block_i%M
	block_o%N = block_i%N
	block_o%headm = block_i%headm
	block_o%headn = block_i%headn

	block_o%M_loc = block_i%M_loc
	block_o%M_loc_db = block_i%M_loc_db
	block_o%N_loc = block_i%N_loc
	block_o%N_loc_db = block_i%N_loc_db
	block_o%pgno = block_i%pgno
	block_o%pgno_db = block_i%pgno_db



	if(associated(block_i%N_p))then
		if(associated(block_o%N_p))deallocate(block_o%N_p)
		allocate(block_o%N_p(size(block_i%N_p,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%N_p)/1024.0d3
		block_o%N_p = block_i%N_p
	endif
	if(associated(block_i%M_p))then
		if(associated(block_o%M_p))deallocate(block_o%M_p)
		allocate(block_o%M_p(size(block_i%M_p,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%M_p)/1024.0d3
		block_o%M_p = block_i%M_p
	endif
	if(associated(block_i%N_p_db))then
		if(associated(block_o%N_p_db))deallocate(block_o%N_p_db)
		allocate(block_o%N_p_db(size(block_i%N_p_db,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%N_p_db)/1024.0d3
		block_o%N_p_db = block_i%N_p_db
	endif
	if(associated(block_i%M_p_db))then
		if(associated(block_o%M_p_db))deallocate(block_o%M_p_db)
		allocate(block_o%M_p_db(size(block_i%M_p_db,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%M_p_db)/1024.0d3
		block_o%M_p_db = block_i%M_p_db
	endif



	level_butterfly = block_i%level_butterfly
	num_blocks=2**level_butterfly

	if(block_i%style==2)then
		if(allocated(block_i%ButterflyU%blocks))then
			if (level_butterfly/=0) then
				allocate(block_o%ButterflyKerl(level_butterfly))
			end if

			do level=0, level_butterfly+1
				if(level==0)then
					block_o%ButterflyV%num_blk=block_i%ButterflyV%num_blk
					block_o%ButterflyV%nblk_loc=block_i%ButterflyV%nblk_loc
					block_o%ButterflyV%idx=block_i%ButterflyV%idx
					block_o%ButterflyV%inc=block_i%ButterflyV%inc
					allocate(block_o%ButterflyV%blocks(block_o%ButterflyV%nblk_loc))

					do jj=1,block_o%ButterflyV%nblk_loc
						nn=size(block_i%ButterflyV%blocks(jj)%matrix,1)
						rank=size(block_i%ButterflyV%blocks(jj)%matrix,2)
						allocate(block_o%ButterflyV%blocks(jj)%matrix(nn,rank))
						block_o%ButterflyV%blocks(jj)%matrix = block_i%ButterflyV%blocks(jj)%matrix
						if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV%blocks(jj)%matrix)/1024.0d3
					enddo

				else if(level==level_butterfly+1)then
					block_o%ButterflyU%num_blk=block_i%ButterflyU%num_blk
					block_o%ButterflyU%nblk_loc=block_i%ButterflyU%nblk_loc
					block_o%ButterflyU%idx=block_i%ButterflyU%idx
					block_o%ButterflyU%inc=block_i%ButterflyU%inc
					allocate(block_o%ButterflyU%blocks(block_o%ButterflyU%nblk_loc))

					do ii=1,block_o%ButterflyU%nblk_loc
						nn=size(block_i%ButterflyU%blocks(ii)%matrix,1)
						rank=size(block_i%ButterflyU%blocks(ii)%matrix,2)
						allocate(block_o%ButterflyU%blocks(ii)%matrix(nn,rank))
						block_o%ButterflyU%blocks(ii)%matrix = block_i%ButterflyU%blocks(ii)%matrix
						if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU%blocks(ii)%matrix)/1024.0d3
					enddo
				else
					block_o%ButterflyKerl(level)%num_row=block_i%ButterflyKerl(level)%num_row
					block_o%ButterflyKerl(level)%num_col=block_i%ButterflyKerl(level)%num_col
					block_o%ButterflyKerl(level)%nc=block_i%ButterflyKerl(level)%nc
					block_o%ButterflyKerl(level)%nr=block_i%ButterflyKerl(level)%nr
					block_o%ButterflyKerl(level)%idx_c=block_i%ButterflyKerl(level)%idx_c
					block_o%ButterflyKerl(level)%idx_r=block_i%ButterflyKerl(level)%idx_r
					block_o%ButterflyKerl(level)%inc_c=block_i%ButterflyKerl(level)%inc_c
					block_o%ButterflyKerl(level)%inc_r=block_i%ButterflyKerl(level)%inc_r

					allocate(block_o%ButterflyKerl(level)%blocks(block_o%ButterflyKerl(level)%nr,block_o%ButterflyKerl(level)%nc))

					do ii=1,block_o%ButterflyKerl(level)%nr
						do jj=1,block_o%ButterflyKerl(level)%nc
							nn=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,2)
							rank=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,1)
							allocate(block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix(rank,nn))
							block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix=block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix
							if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix)/1024.0d3
						enddo
					enddo
				endif
			enddo
		endif

		if(allocated(block_i%ButterflyMiddle))then
			write(*,*)'ButterflyMiddle is not yet handled when the butterfly is distributed'
			levelm = ceiling_safe(dble(level_butterfly)/2d0)
			allocate(block_o%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))

			do index_i_m=1, 2**levelm
				do index_j_m=1, 2**(level_butterfly-levelm)
					rank = size(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,1)
					allocate(block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
					block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix
				end do
			end do
		end if


	else if(block_i%style==1)then
		if(allocated(block_i%fullmat))then
			mm = size(block_i%fullmat,1)
			nn = size(block_i%fullmat,2)
			allocate (block_o%fullmat(mm,nn))
			block_o%fullmat = block_i%fullmat
			if(present(memory))memory = memory + SIZEOF(block_o%fullmat)/1024.0d3
		endif
	else
		! write(*,*)'block style not implemented'
		! stop
	end if
else if(trans=='T')then
write(*,*)'transposed copy is not well tested if the butterfly is distributed'
stop

	block_o%level = block_i%level
	block_o%col_group = block_i%row_group
	block_o%row_group = block_i%col_group
	block_o%style = block_i%style
	block_o%level_butterfly = block_i%level_butterfly
	block_o%rankmax = block_i%rankmax
	block_o%rankmin = block_i%rankmin
	block_o%dimension_rank = block_i%dimension_rank
	block_o%M = block_i%N
	block_o%N = block_i%M
	block_o%headm = block_i%headn
	block_o%headn = block_i%headm

	block_o%M_loc = block_i%N_loc
	block_o%M_loc_db = block_i%N_loc_db
	block_o%N_loc = block_i%M_loc
	block_o%N_loc_db = block_i%M_loc_db
	block_o%pgno = block_i%pgno
	block_o%pgno_db = block_i%pgno_db



	if(associated(block_i%N_p))then
		if(associated(block_o%M_p))deallocate(block_o%M_p)
		allocate(block_o%M_p(size(block_i%N_p,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%M_p)/1024.0d3
		block_o%M_p = block_i%N_p
	endif
	if(associated(block_i%M_p))then
		if(associated(block_o%N_p))deallocate(block_o%N_p)
		allocate(block_o%N_p(size(block_i%M_p,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%N_p)/1024.0d3
		block_o%N_p = block_i%M_p
	endif
	if(associated(block_i%N_p_db))then
		if(associated(block_o%M_p_db))deallocate(block_o%M_p_db)
		allocate(block_o%M_p_db(size(block_i%N_p_db,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%N_p_db)/1024.0d3
		block_o%M_p_db = block_i%N_p_db
	endif
	if(associated(block_i%M_p_db))then
		if(associated(block_o%N_p_db))deallocate(block_o%N_p_db)
		allocate(block_o%N_p_db(size(block_i%M_p_db,1),2))
		if(present(memory))memory = memory + SIZEOF(block_o%M_p_db)/1024.0d3
		block_o%N_p_db = block_i%M_p_db
	endif



	level_butterfly = block_i%level_butterfly
	num_blocks=2**level_butterfly

	if(block_i%style==2)then
		if(allocated(block_i%ButterflyU%blocks))then
			allocate(block_o%ButterflyU%blocks(num_blocks))
			allocate(block_o%ButterflyV%blocks(num_blocks))
			if (level_butterfly/=0) then
				allocate(block_o%ButterflyKerl(level_butterfly))
			end if

			do level=0, level_butterfly
				index_ij=0
				if (level>0) then
					block_o%ButterflyKerl(level)%num_col=2**level
					block_o%ButterflyKerl(level)%num_row=2**(level_butterfly-level+1)
					allocate(block_o%ButterflyKerl(level)%blocks(2**(level_butterfly-level+1),2**level))
				endif
				do index_i=1, 2**level
					do index_j=1, 2**(level_butterfly-level)
						index_ij=index_ij+1
						if (level==0) then
							nn=size(block_i%ButterflyV%blocks(index_ij)%matrix,1)
							rank=size(block_i%ButterflyV%blocks(index_ij)%matrix,2)
							allocate(block_o%ButterflyU%blocks(index_ij)%matrix(nn,rank))
							block_o%ButterflyU%blocks(index_ij)%matrix = block_i%ButterflyV%blocks(index_ij)%matrix
							if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU%blocks(index_ij)%matrix)/1024.0d3
						else
							nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
							rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
							allocate(block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j-1,index_i)%matrix(nn,rank))
							call copymatT(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j-1,index_i)%matrix,rank,nn)
							if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j-1,index_i)%matrix)/1024.0d3
							nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
							allocate(block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j,index_i)%matrix(nn,rank))
							call copymatT(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j,index_i)%matrix,rank,nn)
							if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level+1)%blocks(2*index_j,index_i)%matrix)/1024.0d3
						endif
						if (level==level_butterfly) then
							mm=size(block_i%ButterflyU%blocks(index_ij)%matrix,1)
							rank=size(block_i%ButterflyU%blocks(index_ij)%matrix,2)
							allocate(block_o%ButterflyV%blocks(index_ij)%matrix(mm,rank))
							block_o%ButterflyV%blocks(index_ij)%matrix = block_i%ButterflyU%blocks(index_ij)%matrix
							if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV%blocks(index_ij)%matrix)/1024.0d3
						endif
					enddo
				enddo
			enddo
		endif

		if(allocated(block_i%ButterflyMiddle))then
			levelm = ceiling_safe(dble(level_butterfly)/2d0)
			allocate(block_o%ButterflyMiddle(2**(level_butterfly-levelm),2**levelm))

			do index_i_m=1, 2**levelm
				do index_j_m=1, 2**(level_butterfly-levelm)
					rank = size(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,1)
					allocate(block_o%ButterflyMiddle(index_j_m,index_i_m)%matrix(rank,rank))
					call copymatT(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,block_o%ButterflyMiddle(index_j_m,index_i_m)%matrix,rank,rank)
				end do
			end do
		end if
	else if(block_i%style==1)then
		if(allocated(block_i%fullmat))then
			mm = size(block_i%fullmat,1)
			nn = size(block_i%fullmat,2)
			allocate (block_o%fullmat(nn,mm))
			call copymatT(block_i%fullmat,block_o%fullmat,mm,nn)
			if(present(memory))memory = memory + SIZEOF(block_o%fullmat)/1024.0d3
		endif
	else
		! write(*,*)'block style not implemented'
		! stop
	end if
endif


end subroutine BF_copy

subroutine BF_copy_delete(block_i,block_o,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i,block_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
real(kind=8),optional::memory
if(present(memory))memory=0
block_o%level = block_i%level
block_o%col_group = block_i%col_group
block_o%row_group = block_i%row_group
block_o%style = block_i%style
block_o%level_butterfly = block_i%level_butterfly
block_o%level_half = block_i%level_half
block_o%rankmax = block_i%rankmax
block_o%rankmin = block_i%rankmin
block_o%M = block_i%M
block_o%N = block_i%N
block_o%headm = block_i%headm
block_o%headn = block_i%headn

block_o%M_loc = block_i%M_loc
block_o%M_loc_db = block_i%M_loc_db
block_o%N_loc = block_i%N_loc
block_o%N_loc_db = block_i%N_loc_db
block_o%pgno = block_i%pgno
block_o%pgno_db = block_i%pgno_db

if(associated(block_i%N_p))then
	allocate(block_o%N_p(size(block_i%N_p,1),2))
	block_o%N_p = block_i%N_p
endif
if(associated(block_i%M_p))then
	allocate(block_o%M_p(size(block_i%M_p,1),2))
	block_o%M_p = block_i%M_p
endif
if(associated(block_i%N_p_db))then
	allocate(block_o%N_p_db(size(block_i%N_p_db,1),2))
	block_o%N_p_db = block_i%N_p_db
endif
if(associated(block_i%M_p_db))then
	allocate(block_o%M_p_db(size(block_i%M_p_db,1),2))
	block_o%M_p_db = block_i%M_p_db
endif


level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

if(block_i%style==2)then

	if (level_butterfly/=0) then
		if(allocated(block_i%ButterflyKerl))allocate(block_o%ButterflyKerl(level_butterfly))
	end if

	do level=0, level_butterfly+1
		if(level==0)then
			block_o%ButterflyV%num_blk=block_i%ButterflyV%num_blk
			block_o%ButterflyV%nblk_loc=block_i%ButterflyV%nblk_loc
			block_o%ButterflyV%idx=block_i%ButterflyV%idx
			block_o%ButterflyV%inc=block_i%ButterflyV%inc
			if(allocated(block_i%ButterflyV%blocks))allocate(block_o%ButterflyV%blocks(block_o%ButterflyV%nblk_loc))

			do jj=1,block_o%ButterflyV%nblk_loc
			if(allocated(block_i%ButterflyV%blocks(jj)%matrix))then
				nn=size(block_i%ButterflyV%blocks(jj)%matrix,1)
				rank=size(block_i%ButterflyV%blocks(jj)%matrix,2)
				allocate(block_o%ButterflyV%blocks(jj)%matrix(nn,rank))
				block_o%ButterflyV%blocks(jj)%matrix = block_i%ButterflyV%blocks(jj)%matrix
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV%blocks(jj)%matrix)/1024.0d3
				deallocate(block_i%ButterflyV%blocks(jj)%matrix)
			endif
			enddo

		else if(level==level_butterfly+1)then
			block_o%ButterflyU%num_blk=block_i%ButterflyU%num_blk
			block_o%ButterflyU%nblk_loc=block_i%ButterflyU%nblk_loc
			block_o%ButterflyU%idx=block_i%ButterflyU%idx
			block_o%ButterflyU%inc=block_i%ButterflyU%inc
			if(allocated(block_i%ButterflyU%blocks))allocate(block_o%ButterflyU%blocks(block_o%ButterflyU%nblk_loc))

			do ii=1,block_o%ButterflyU%nblk_loc
				if(allocated(block_i%ButterflyU%blocks(ii)%matrix))then
				nn=size(block_i%ButterflyU%blocks(ii)%matrix,1)
				rank=size(block_i%ButterflyU%blocks(ii)%matrix,2)
				allocate(block_o%ButterflyU%blocks(ii)%matrix(nn,rank))
				block_o%ButterflyU%blocks(ii)%matrix = block_i%ButterflyU%blocks(ii)%matrix
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU%blocks(ii)%matrix)/1024.0d3
				deallocate(block_i%ButterflyU%blocks(ii)%matrix)
				endif
			enddo
		else
			block_o%ButterflyKerl(level)%num_row=block_i%ButterflyKerl(level)%num_row
			block_o%ButterflyKerl(level)%num_col=block_i%ButterflyKerl(level)%num_col
			block_o%ButterflyKerl(level)%nc=block_i%ButterflyKerl(level)%nc
			block_o%ButterflyKerl(level)%nr=block_i%ButterflyKerl(level)%nr
			block_o%ButterflyKerl(level)%idx_c=block_i%ButterflyKerl(level)%idx_c
			block_o%ButterflyKerl(level)%idx_r=block_i%ButterflyKerl(level)%idx_r
			block_o%ButterflyKerl(level)%inc_c=block_i%ButterflyKerl(level)%inc_c
			block_o%ButterflyKerl(level)%inc_r=block_i%ButterflyKerl(level)%inc_r

			if(allocated(block_i%ButterflyKerl(level)%blocks))then
			allocate(block_o%ButterflyKerl(level)%blocks(block_o%ButterflyKerl(level)%nr,block_o%ButterflyKerl(level)%nc))
			do ii=1,block_o%ButterflyKerl(level)%nr
				do jj=1,block_o%ButterflyKerl(level)%nc
					if(allocated(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix))then
					nn=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,1)
					allocate(block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix(rank,nn))
					block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix=block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(ii,jj)%matrix)/1024.0d3
					deallocate(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix)
					endif
				enddo
			enddo
			deallocate(block_i%ButterflyKerl(level)%blocks)
			endif
		endif
	enddo
	if(allocated(block_i%ButterflyU%blocks))deallocate(block_i%ButterflyU%blocks)
	if(allocated(block_i%ButterflyV%blocks))deallocate(block_i%ButterflyV%blocks)
	if(allocated(block_i%ButterflyKerl))then
	if(level_butterfly/=0)deallocate(block_i%ButterflyKerl)
	endif
	if(allocated(block_i%ButterflyMiddle))then
		write(*,*)'ButterflyMiddle is not yet handled when the butterfly is distributed'
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		allocate(block_o%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))

		do index_i_m=1, 2**levelm
			do index_j_m=1, 2**(level_butterfly-levelm)
				rank = size(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,1)
				allocate(block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix
				deallocate(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix)
			end do
		end do
		deallocate(block_i%ButterflyMiddle)
	end if

else if(block_i%style==1)then
	mm = size(block_i%fullmat,1)
	nn = size(block_i%fullmat,2)
	allocate (block_o%fullmat(mm,nn))
	block_o%fullmat = block_i%fullmat
	deallocate(block_i%fullmat)
else
	write(*,*)'block style not implemented'
	stop
end if

end subroutine BF_copy_delete

subroutine BF_ComputeMemory(block_i,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
real(kind=8)::memory
memory=0

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

if(block_i%style==2)then

	do level=0, level_butterfly+1
		if(level==0)then
			do jj=1,block_i%ButterflyV%nblk_loc
				memory = memory + SIZEOF(block_i%ButterflyV%blocks(jj)%matrix)/1024.0d3
			enddo
		elseif(level==level_butterfly+1)then
			do jj=1,block_i%ButterflyU%nblk_loc
				memory = memory + SIZEOF(block_i%ButterflyU%blocks(jj)%matrix)/1024.0d3
			enddo
		else
			do ii=1, block_i%ButterflyKerl(level)%nr
				do jj=1, block_i%ButterflyKerl(level)%nc
					memory = memory + SIZEOF(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix)/1024.0d3
				enddo
			enddo
		endif
	enddo

else if(block_i%style==1)then
	memory = memory + SIZEOF(block_i%fullmat)/1024.0d3
else
	write(*,*)'block style not implemented'
	stop
end if

end subroutine BF_ComputeMemory



integer function BF_Switchlevel(level_butterfly,option)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(Hoption)::option
integer::level_butterfly

if(option%pat_comp==1)BF_Switchlevel = level_butterfly ! from right to left until the second last level
if(option%pat_comp==2)BF_Switchlevel = 0 ! from left to right until the second last level
if(option%pat_comp==3)BF_Switchlevel = floor_safe(dble(level_butterfly)/2d0)  ! from outer to inner

end function BF_Switchlevel


logical function BF_checkNAN(block_i)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
real(kind=8):: temp

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly
temp = 0

if(block_i%style==2)then
	do level=0, level_butterfly+1
		if (level==0) then
			if(allocated(block_i%ButterflyV%blocks))then
			do index_j=1,block_i%ButterflyV%nblk_loc
				if(allocated(block_i%ButterflyV%blocks(index_j)%matrix))then
				mm=size(block_i%ButterflyV%blocks(index_j)%matrix,1)
				nn=size(block_i%ButterflyV%blocks(index_j)%matrix,2)
				temp = temp + fnorm(block_i%ButterflyV%blocks(index_j)%matrix,mm,nn)
				! write(*,*)'V',level_butterfly,index_j,fnorm(block_i%ButterflyV%blocks(index_j)%matrix,mm,nn),mm,nn
				endif
			enddo
			endif
		elseif (level==level_butterfly+1) then
			if(allocated(block_i%ButterflyU%blocks))then
			do index_i=1,block_i%ButterflyU%nblk_loc
				if(allocated(block_i%ButterflyU%blocks(index_i)%matrix))then
				mm=size(block_i%ButterflyU%blocks(index_i)%matrix,1)
				nn=size(block_i%ButterflyU%blocks(index_i)%matrix,2)
				temp = temp + fnorm(block_i%ButterflyU%blocks(index_i)%matrix,mm,nn)
				! write(*,*)'U',level_butterfly,index_i,fnorm(block_i%ButterflyU%blocks(index_i)%matrix,mm,nn),mm,nn
				endif
			enddo
			endif
		else
			if(allocated(block_i%ButterflyKerl))then
			if(allocated(block_i%ButterflyKerl(level)%blocks))then
			do index_i=1, block_i%ButterflyKerl(level)%nr
				do index_j=1, block_i%ButterflyKerl(level)%nc
					if(allocated(block_i%ButterflyKerl(level)%blocks(index_i,index_j)%matrix))then
					mm=size(block_i%ButterflyKerl(level)%blocks(index_i,index_j)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,index_j)%matrix,2)
					temp = temp + fnorm(block_i%ButterflyKerl(level)%blocks(index_i,index_j)%matrix,mm,nn)
					! write(*,*)'Ker',level_butterfly,level,index_i,index_j,fnorm(block_i%ButterflyKerl(level)%blocks(index_i,index_j)%matrix,mm,nn),mm,nn
					endif
				enddo
			enddo
			endif
			endif
		endif
	enddo
else if(block_i%style==1)then
	if(allocated(block_i%fullmat))then
	mm=size(block_i%fullmat,1)
	nn=size(block_i%fullmat,2)
	temp = temp + fnorm(block_i%fullmat,mm,nn)
	endif
else
	write(*,*)'block style not implemented'
	stop
end if

BF_checkNAN = isnan(temp)

end function BF_checkNAN



subroutine BF_print_size_rank(block_i,tolerance)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,truerank,index_i,index_j,levelm,index_i_m,index_j_m,mm1,mm2,nn1,nn2
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
DT,allocatable::matrixtemp(:,:),mat11(:,:),mat12(:,:),mat21(:,:),mat22(:,:)
real(kind=8)::tolerance

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

do level=0, level_butterfly+1
	! write(*,*)level
	if (level==0) then
		do index_ij=1, 2**level_butterfly
			nn=size(block_i%ButterflyV%blocks(index_ij)%matrix,1)
			rank=size(block_i%ButterflyV%blocks(index_ij)%matrix,2)
			allocate(matrixtemp(nn,rank))
			matrixtemp = block_i%ButterflyV%blocks(index_ij)%matrix
			call GetRank(nn,rank,matrixtemp,truerank,tolerance)
			write(*,*)level,index_ij,nn,rank,truerank
			deallocate(matrixtemp)
		end do
	else if (level==level_butterfly+1) then
		do index_ij=1, 2**level_butterfly
			mm=size(block_i%ButterflyU%blocks(index_ij)%matrix,1)
			rank=size(block_i%ButterflyU%blocks(index_ij)%matrix,2)
			allocate(matrixtemp(mm,rank))
			matrixtemp = block_i%ButterflyU%blocks(index_ij)%matrix
			call GetRank(mm,rank,matrixtemp,truerank,tolerance)
			write(*,*)level,index_ij,mm,rank,truerank
			deallocate(matrixtemp)
		end do
	else
		do index_i=1, 2**(level-1)
			do index_j=1, 2**(level_butterfly-level)

				mm1=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix,1)
				nn1=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix,2)

				mm2=size(block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix,1)
				nn2=size(block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix,2)

				allocate(mat11(mm1,nn1))
				mat11 = block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix
				allocate(mat12(mm1,nn2))
				mat12 = block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j)%matrix
				allocate(mat21(mm2,nn1))
				mat21 = block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j-1)%matrix
				allocate(mat22(mm2,nn2))
				mat22 = block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix
				allocate(matrixtemp(mm1+mm2,nn1+nn2))
				matrixtemp(1:mm1,1:nn1) = mat11
				matrixtemp(1:mm1,1+nn1:nn2+nn1) = mat12
				matrixtemp(1+mm1:mm2+mm1,1:nn1) = mat21
				matrixtemp(1+mm1:mm2+mm1,1+nn1:nn2+nn1) = mat22
				call GetRank(mm1+mm2,nn1+nn2,matrixtemp,truerank,tolerance)
				write(*,*)level,index_i,index_j,(mm1+mm2),(nn1+nn2),truerank
				deallocate(mat11)
				deallocate(mat12)
				deallocate(mat21)
				deallocate(mat22)
				deallocate(matrixtemp)

			enddo
		enddo

	end if

enddo

end subroutine BF_print_size_rank



subroutine BF_extract_partial(block_o,level_butterfly_loc,ij_loc,LR,agent_block)
	use MISC_Utilities
    use BPACK_DEFS
    implicit none

	type(matrixblock)::block_o,agent_block
	integer level_butterfly,level_butterfly_loc, ij_loc,index_i,index_i_start,index_j_start,index_j,level,ii,nn,mm,num_blocks,rank
	character LR

	! allocate(agent_block)


	call assert(level_butterfly_loc>=1,'level_butterfly_loc cannot be zero')

	agent_block%row_group=-1
	agent_block%col_group=-1

	agent_block%style = block_o%style
	agent_block%level_butterfly = level_butterfly_loc
	agent_block%rankmax = block_o%rankmax
	agent_block%rankmin = block_o%rankmin
	level_butterfly = block_o%level_butterfly



	num_blocks=2**level_butterfly_loc




	allocate(agent_block%ButterflyU%blocks(num_blocks))
	allocate(agent_block%ButterflyV%blocks(num_blocks))

	allocate(agent_block%ButterflyKerl(level_butterfly_loc))


	if(LR=='L')then
		do level=1, level_butterfly_loc
			agent_block%ButterflyKerl(level)%num_row=2**level
			agent_block%ButterflyKerl(level)%num_col=2**(level_butterfly_loc-level+1)
			allocate(agent_block%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly_loc-level+1)))
			do index_i=1, 2**level
				do index_j=1, 2**(level_butterfly_loc-level)

					index_i_start = (ij_loc-1)*2**level
					nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,2)
					rank=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,1)
					allocate(agent_block%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn))
					agent_block%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix = block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix

					nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,2)
					allocate(agent_block%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn))
					agent_block%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix = block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix

					if (level==level_butterfly_loc) then
						index_i_start = (ij_loc-1)*2**level

						mm=size(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix,1)
						rank=size(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix,2)
						allocate(agent_block%ButterflyU%blocks(index_i)%matrix(mm,rank))
						agent_block%ButterflyU%blocks(index_i)%matrix = block_o%ButterflyU%blocks(index_i+index_i_start)%matrix
					endif
				enddo
			enddo

			if(level==1)then
				do index_i=1, 1
					do index_j=1, 2**(level_butterfly_loc-level)
						index_i_start = (ij_loc-1)*2**level
						nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,2)
						allocate(agent_block%ButterflyV%blocks(2*index_j-1)%matrix(nn,nn))
						agent_block%ButterflyV%blocks(2*index_j-1)%matrix = 0
						do ii=1,nn
							agent_block%ButterflyV%blocks(2*index_j-1)%matrix(ii,ii)=1
						end do
						nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,2)
						allocate(agent_block%ButterflyV%blocks(2*index_j)%matrix(nn,nn))
						agent_block%ButterflyV%blocks(2*index_j)%matrix = 0
						do ii=1,nn
							agent_block%ButterflyV%blocks(2*index_j)%matrix(ii,ii)=1
						end do
					end do
				end do
			end if

		enddo
	else if(LR=='R')then
		do level=1, level_butterfly_loc
			agent_block%ButterflyKerl(level)%num_row=2**level
			agent_block%ButterflyKerl(level)%num_col=2**(level_butterfly_loc-level+1)
			allocate(agent_block%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly_loc-level+1)))
			do index_i=1, 2**(level-1)
				do index_j=1, 2**(level_butterfly_loc-level+1)

					index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)

					mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,1)
					rank=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,2)
					allocate(agent_block%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm,rank))
					agent_block%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix = block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix

					mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,1)
					allocate(agent_block%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm,rank))
					agent_block%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix = block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix

					if (level==1) then
						index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)

						nn=size(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix,1)
						rank=size(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix,2)
						allocate(agent_block%ButterflyV%blocks(index_j)%matrix(nn,rank))
						agent_block%ButterflyV%blocks(index_j)%matrix = block_o%ButterflyV%blocks(index_j+index_j_start)%matrix
					endif
				enddo
			enddo

			if(level==level_butterfly_loc)then
				do index_i=1, 2**(level_butterfly_loc-1)
					do index_j=1, 1
						index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)
						mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,1)
						allocate(agent_block%ButterflyU%blocks(2*index_i-1)%matrix(mm,mm))
						agent_block%ButterflyU%blocks(2*index_i-1)%matrix = 0
						do ii=1,mm
							agent_block%ButterflyU%blocks(2*index_i-1)%matrix(ii,ii)=1
						end do
						mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,1)
						allocate(agent_block%ButterflyU%blocks(2*index_i)%matrix(mm,mm))
						agent_block%ButterflyU%blocks(2*index_i)%matrix = 0
						do ii=1,mm
							agent_block%ButterflyU%blocks(2*index_i)%matrix(ii,ii)=1
						end do
					end do
				end do
			end if

		enddo

	end if
end subroutine BF_extract_partial


subroutine BF_copy_partial(block_i,block_o,level_butterfly_loc,ij_loc,LR,memory)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_o,block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,dimension_m,dimension_n,index_i_start,index_j_start
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,level_butterfly_loc,ij_loc
character LR
real(kind=8),optional::memory
if(present(memory))memory=0

!!!!! be careful here, may need changes later
block_o%rankmax = max(block_o%rankmax,block_i%rankmax)
block_o%rankmin = max(block_o%rankmin,block_i%rankmin)


call assert(level_butterfly_loc>=1,'level_butterfly_loc cannot be zero')
call assert(level_butterfly_loc==block_i%level_butterfly,'level_butterfly_loc/=block_i%level_butterfly')

level_butterfly = block_o%level_butterfly
num_blocks=2**level_butterfly_loc


if(LR=='L')then

	do level=1, level_butterfly_loc
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly_loc-level)
				index_i_start = (ij_loc-1)*2**level

				if(level==1)then
					dimension_n = size(block_i%ButterflyV%blocks(2*index_j-1)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix(rank,dimension_n))
					! call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix, block_i%ButterflyV%blocks(2*index_j-1)%matrix, &
					! &block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,dimension_n,nn)
					call gemmf90(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,rank, block_i%ButterflyV%blocks(2*index_j-1)%matrix,dimension_n, block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank, 'N','T',rank,dimension_n,nn,cone,czero)



					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,dimension_n)))then
						write(*,*)'NAN in L 1'
					end if


					dimension_n = size(block_i%ButterflyV%blocks(2*index_j)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix(rank,dimension_n))
					! call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix, block_i%ButterflyV%blocks(2*index_j)%matrix, &
					! &block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,dimension_n,nn)
					call gemmf90(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,rank, block_i%ButterflyV%blocks(2*index_j)%matrix,dimension_n, block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank, 'N','T',rank,dimension_n,nn,cone,czero)

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,dimension_n)))then
						write(*,*)'NAN in L 2'
					end if


				else
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix(rank,nn))
					block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,nn)))then
						write(*,*)'NAN in L 3'
					end if

					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix(rank,nn))
					block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,nn)))then
						write(*,*)'NAN in L 4'
					end if

				end if

				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)/1024.0d3
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)/1024.0d3

				if (level==level_butterfly_loc) then
					index_i_start = (ij_loc-1)*2**level
					mm=size(block_i%ButterflyU%blocks(index_i)%matrix,1)
					rank=size(block_i%ButterflyU%blocks(index_i)%matrix,2)
					deallocate(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix)
					allocate(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix(mm,rank))
					block_o%ButterflyU%blocks(index_i+index_i_start)%matrix = block_i%ButterflyU%blocks(index_i)%matrix
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix)/1024.0d3
					if(isnan(fnorm(block_o%ButterflyU%blocks(index_i+index_i_start)%matrix,mm,rank)))then
						write(*,*)'NAN in L 5'
					end if
				endif
			enddo
		enddo
	enddo

else if(LR=='R')then


	do level=1, level_butterfly_loc
		do index_i=1, 2**(level-1)
			do index_j=1, 2**(level_butterfly_loc-level+1)
			! write(*,*)level,index_i,index_j
				index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)
				if(level==level_butterfly_loc)then
				! write(*,*)'good 1'
					dimension_m = size(block_i%ButterflyU%blocks(2*index_i-1)%matrix,1)
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)
					! write(*,*)dimension_m,mm,rank,'d'
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix(dimension_m,rank))
					! call gemm_omp(block_i%ButterflyU%blocks(2*index_i-1)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,&
					! &block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,rank,mm)

					call gemmf90(block_i%ButterflyU%blocks(2*index_i-1)%matrix,dimension_m,block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,mm,block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,'N','N',dimension_m,rank,mm,cone,czero)

! write(*,*)'good 1.1'

					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,rank)))then
						write(*,*)'NAN in R 1'
					end if

					dimension_m = size(block_i%ButterflyU%blocks(2*index_i)%matrix,1)
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix(dimension_m,rank))
					! call gemm_omp(block_i%ButterflyU%blocks(2*index_i)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,&
					! &block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,rank,mm)

					call gemmf90(block_i%ButterflyU%blocks(2*index_i)%matrix,dimension_m,block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,mm,block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,'N','N',dimension_m,rank,mm,cone,czero)

! write(*,*)'good 2'
					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,rank)))then
						write(*,*)'NAN in R 2'
					end if
				else
				! write(*,*)'good 3'
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix(mm,rank))
					block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix

					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,mm,rank)))then
						write(*,*)'NAN in R 3'
					end if


					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix(mm,rank))
					block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix
				! write(*,*)'good 4'
					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,mm,rank)))then
						write(*,*)'NAN in R 4'
					end if
				end if

				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)/1024.0d3
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)/1024.0d3

				if (level==1) then
					index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)
					nn=size(block_i%ButterflyV%blocks(index_j)%matrix,1)
					rank=size(block_i%ButterflyV%blocks(index_j)%matrix,2)
					deallocate(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix(nn,rank))
					block_o%ButterflyV%blocks(index_j+index_j_start)%matrix = block_i%ButterflyV%blocks(index_j)%matrix
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix)/1024.0d3

					if(isnan(fnorm(block_o%ButterflyV%blocks(index_j+index_j_start)%matrix,nn,rank)))then
						write(*,*)'NAN in R 5'
					end if
				endif
			end do
		end do
	end do

end if

end subroutine BF_copy_partial




subroutine BF_Partial_MVP_Half(block_rand,chara,level_start,level_end,random,num_vect_sub,nth_s,nth_e,Ng)

    use BPACK_DEFS
    implicit none

    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, ij, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, level_end, level_start
    DT ctemp, a, b
    character chara
	integer num_vect_sub,num_vect_subsub,nth_s,nth_e,Ng,nth,dimension_rank,level_butterfly

    type(RandomBlock) :: random

    type(matrixblock)::block_rand
   ! write(*,*)'in '

	level_butterfly=block_rand%level_butterfly
	dimension_rank = block_rand%dimension_rank
    num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)

    if (chara=='N') then

        num_blocks=2**level_butterfly

        do level=level_start, level_end
            if (level==0) then
                num_groupn=num_blocks
                do nth=nth_s, nth_e
					! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
					do j = (nth-1)*Ng+1,nth*Ng
						rank=size(block_rand%ButterflyV%blocks(j)%matrix,2)
						nn=size(block_rand%ButterflyV%blocks(j)%matrix,1)
						if(.not. allocated(random%RandomVectorRR(1)%blocks(1,j)%matrix))allocate(random%RandomVectorRR(1)%blocks(1,j)%matrix(rank,num_vect_subsub))
						random%RandomVectorRR(1)%blocks(1,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vect_subsub
							do ii=1, rank
								ctemp=0d0
								do kk=1, nn
									ctemp=ctemp+block_rand%ButterflyV%blocks(j)%matrix(kk,ii)*random%RandomVectorRR(0)%blocks(1,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorRR(1)%blocks(1,j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						!$omp end parallel do

					enddo
					! !$omp end parallel do
                enddo
            elseif (level==level_butterfly+1) then

            else
                num_groupm=block_rand%ButterflyKerl(level)%num_row
                num_groupn=block_rand%ButterflyKerl(level)%num_col

				! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm,nth)
				do ij=1,(num_groupm/2)*(num_groupn/2)
					index_i = (ij-1)/(num_groupn/2)+1
					index_j = mod(ij-1,(num_groupn/2)) + 1
					i = index_i*2-1
					j = index_j*2-1

					do nth = nth_s,nth_e

						if((j>=(nth-1)*Ng/2**(level-1)+1 .and. j<=nth*Ng/2**(level-1)) .or. &
						& (j+1>=(nth-1)*Ng/2**(level-1)+1 .and. j+1<=nth*Ng/2**(level-1)))then
							nn1=size(block_rand%ButterflyKerl(level)%blocks(i,j)%matrix,2)
							nn2=size(block_rand%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
							mm=size(block_rand%ButterflyKerl(level)%blocks(i,j)%matrix,1)
							! write(*,*)ij,i,j,level,'ha',index_i
							if(.not. allocated(random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix))allocate(random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(mm,num_vect_subsub))
							random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do jj=1, num_vect_subsub
								do ii=1, mm
									ctemp=0d0
									do kk=1, nn1
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									do kk=1, nn2
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do

							nn1=size(block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix,2)
							nn2=size(block_rand%ButterflyKerl(level)%blocks(i+1,j+1)%matrix,2)
							mm=size(block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix,1)
							! write(*,*)ij,i,j,level,'ha',index_i
							if(.not. allocated(random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix))allocate(random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(mm,num_vect_subsub))
							random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do jj=1, num_vect_subsub
								do ii=1, mm
									ctemp=0d0
									do kk=1, nn1
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									do kk=1, nn2
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do

							! write(*,*)ij,i,j,level,'ha done0',index_i
							deallocate(random%RandomVectorRR(level)%blocks(index_i,j)%matrix)
							deallocate(random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix)
							! write(*,*)ij,i,j,level,'ha done',index_i
						end if
					end do
				enddo
				! !$omp end parallel do
            endif
        enddo


    elseif (chara=='T') then

        num_blocks=2**level_butterfly

        do level=level_start, level_end
            if (level==0) then
                num_groupm=num_blocks
                do nth=nth_s, nth_e
					! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
					do i = (nth-1)*Ng+1,nth*Ng
						rank=size(block_rand%ButterflyU%blocks(i)%matrix,2)
						mm=size(block_rand%ButterflyU%blocks(i)%matrix,1)
						if(.not. allocated(random%RandomVectorLL(1)%blocks(i,1)%matrix))allocate(random%RandomVectorLL(1)%blocks(i,1)%matrix(rank,num_vect_subsub))
						random%RandomVectorLL(1)%blocks(i,1)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vect_subsub
							do ii=1, rank
								ctemp=0d0
								do kk=1, mm
									ctemp=ctemp+block_rand%ButterflyU%blocks(i)%matrix(kk,ii)*random%RandomVectorLL(0)%blocks(i,1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorLL(1)%blocks(i,1)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						!$omp end parallel do
					end do
					! !$omp end parallel do
                enddo
            elseif (level==level_butterfly+1) then
            else
                num_groupm=block_rand%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=block_rand%ButterflyKerl(level_butterfly-level+1)%num_col

				! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn,nth)
				do ij=1,(num_groupn/2)*(num_groupm/2)
					index_j = (ij-1)/(num_groupm/2)+1
					index_i = mod(ij-1,(num_groupm/2)) + 1
					j = 2*index_j-1
					i = 2*index_i-1




					do nth = nth_s,nth_e

						if((i>=(nth-1)*Ng/2**(level-1)+1 .and. i<=nth*Ng/2**(level-1)) .or. &
						& (i+1>=(nth-1)*Ng/2**(level-1)+1 .and. i+1<=nth*Ng/2**(level-1)))then
							mm1=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
							mm2=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
							nn=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
							if(.not. allocated(random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix))allocate(random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(nn,num_vect_subsub))
							random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do ii=1, num_vect_subsub
								do jj=1, nn
									ctemp=0d0
									do kk=1, mm1
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
									enddo
									do kk=1, mm2
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
									enddo
									random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(jj,ii+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do
							mm1=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix,1)
							mm2=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j+1)%matrix,1)
							nn=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix,2)
							if(.not. allocated(random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix))allocate(random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(nn,num_vect_subsub))
							random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do ii=1, num_vect_subsub
								do jj=1, nn
									ctemp=0d0
									do kk=1, mm1
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix(kk,jj)
									enddo
									do kk=1, mm2
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j+1)%matrix(kk,jj)
									enddo
									random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(jj,ii+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do
							deallocate(random%RandomVectorLL(level)%blocks(i,index_j)%matrix)
							deallocate(random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix)
						end if
					end do
				enddo
				! !$omp end parallel do

            endif
        enddo

    endif
       ! write(*,*)'out '
    return

end subroutine BF_Partial_MVP_Half




subroutine BF_exchange_extraction(blocks,kerls,stats,ptree,level,collect)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i,index_i0,index_i1, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j1,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::blocks
	type(Hstat)::stats
	type(proctree)::ptree
	type(butterfly_kerl)::kerls

    integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    integer,allocatable:: select_row_rr(:), select_column_rr(:)
    DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:),core(:,:),tau(:)

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pids,pidr,pid,pid0,tag,nproc,Ncol,Ncol1,Nrow,Nreqr,Nreqs,recvid,sendid

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive
	logical::sendflag,recvflag
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	integer rr,cc
	character mode,modetrans,collect
	DT,allocatable::mat1(:,:),mat2(:,:),mat(:,:)

	real(kind=8)::n1,n2

	n1 = OMP_get_wtime()

	mode='R'
	modetrans='C'


	level_butterfly=blocks%level_butterfly
	nproc = ptree%pgrp(blocks%pgno)%nproc
	tag = blocks%pgno


	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0

	! calculate send buffer sizes in the first pass

	do nn=1,size(kerls%index,1)
		ii = kerls%index(nn,1)
		jj = kerls%index(nn,2)
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c

		sendflag=.false.
		recvflag=.false.
		if(collect=='R')then ! pair-wise reduction
			if(mode=='R')then
				index_i0=floor_safe((index_i-1)/2d0)+1
				index_j0=2*index_j-mod(index_i,2)
				index_i1=floor_safe((index_i-1)/2d0)+1
				index_j1=2*index_j-mod(index_i-1,2)
			elseif(mode=='C')then
				write(*,*)'mode=C not needed in BF_exchange_extraction'
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,modetrans,pids)
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i1,index_j1,modetrans,pidr)
		elseif(collect=='B')then ! pair-wise broadcast
			if(mode=='R')then
				index_j0 = index_j+2*mod(index_j,2)-1
				index_i0=index_i
				index_j1 = index_j
				index_i1=index_i
			elseif(mode=='C')then
				write(*,*)'mode=C not needed in BF_exchange_extraction'
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,mode,pids)
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i1,index_j1,mode,pidr)
		endif
		sendflag = pids/=ptree%MyID
		recvflag = pidr/=ptree%MyID

		if(recvflag)then
			pp=pidr-ptree%pgrp(blocks%pgno)%head+1
			if(recvquant(pp)%active==0)then
				recvquant(pp)%active=1
				Nrecvactive=Nrecvactive+1
				recvIDactive(Nrecvactive)=pp
			endif
		endif

		if(sendflag)then
			pp=pids-ptree%pgrp(blocks%pgno)%head+1
			if(sendquant(pp)%active==0)then
				sendquant(pp)%active=1
				Nsendactive=Nsendactive+1
				sendIDactive(Nsendactive)=pp
			endif
			if(allocated(kerls%blocks(ii,jj)%matrix))then
				sendquant(pp)%size=sendquant(pp)%size+4+size(kerls%blocks(ii,jj)%matrix,1)*size(kerls%blocks(ii,jj)%matrix,2)
			endif
		endif
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif


	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do nn=1,size(kerls%index,1)
		ii = kerls%index(nn,1)
		jj = kerls%index(nn,2)
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c

		sendflag=.false.
		if(collect=='R')then ! pair-wise reduction
			if(mode=='R')then
				index_i0=floor_safe((index_i-1)/2d0)+1
				index_j0=2*index_j-mod(index_i,2)
			elseif(mode=='C')then
				write(*,*)'mode=C not needed in BF_exchange_extraction'
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,modetrans,pid)
		elseif(collect=='B')then ! pair-wise broadcast

			if(mode=='R')then
				index_j0 = index_j+2*mod(index_j,2)-1
				index_i0=index_i
			elseif(mode=='C')then
				write(*,*)'mode=C not needed in BF_exchange_extraction'
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,mode,pid)
		endif
		sendflag = pid/=ptree%MyID

		if(sendflag)then
			pp=pid-ptree%pgrp(blocks%pgno)%head+1
			if(allocated(kerls%blocks(ii,jj)%matrix))then
				Nrow=size(kerls%blocks(ii,jj)%matrix,1)
				Ncol=size(kerls%blocks(ii,jj)%matrix,2)
				sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
				sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
				sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
				sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
				sendquant(pp)%size=sendquant(pp)%size+4
				do i=1,Nrow*Ncol
					rr = mod(i-1,Nrow)+1
					cc = (i-1)/Nrow+1
					sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls%blocks(ii,jj)%matrix(rr,cc)
				enddo
				sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
			endif
		endif
	enddo


	! communicate the data buffer
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo

	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			ii=(index_i-kerls%idx_r)/kerls%inc_r+1
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jj=(index_j-kerls%idx_c)/kerls%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))

			if(allocated(kerls%blocks(ii,jj)%matrix))then
				Ncol1=size(kerls%blocks(ii,jj)%matrix,2)
				allocate(mat(Nrow,Ncol1+Ncol))
				mat(1:Nrow,1:Ncol1)=kerls%blocks(ii,jj)%matrix
				deallocate(kerls%blocks(ii,jj)%matrix)
			else
				Ncol1=0
				allocate(mat(Nrow,Ncol1+Ncol))
			endif
			allocate(kerls%blocks(ii,jj)%matrix(Nrow,Ncol1+Ncol))

			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				mat(rr,cc+Ncol1) = recvquant(pp)%dat(i+j,1)
			enddo
			kerls%blocks(ii,jj)%matrix = mat
			deallocate(mat)
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_exchange_extraction





subroutine BF_exchange_matvec(blocks,kerls,Ncol,stats,ptree,level,mode,collect)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i,index_i0, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::blocks
	type(Hstat)::stats
	type(proctree)::ptree
	type(butterfly_kerl)::kerls

    integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    integer,allocatable:: select_row_rr(:), select_column_rr(:)
    DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:),core(:,:),tau(:)

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,pid0,tag,nproc,Ncol,Nrow,Nreqr,Nreqs,recvid,sendid

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive
	logical::sendflag
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	integer rr,cc
	character mode,modetrans,collect

	real(kind=8)::n1,n2

	n1 = OMP_get_wtime()

	if(mode=='R')modetrans='C'
	if(mode=='C')modetrans='R'

	level_butterfly=blocks%level_butterfly
	nproc = ptree%pgrp(blocks%pgno)%nproc
	tag = blocks%pgno

	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0

	! calculate send buffer sizes in the first pass
	do ii=1,kerls%nr
	do jj=1,kerls%nc
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c

		sendflag=.false.
		if(collect=='R')then ! pair-wise reduction
			if(mode=='R')then
				index_i0=floor_safe((index_i-1)/2d0)+1
				index_j0=2*index_j-mod(index_i,2)
			elseif(mode=='C')then
				index_i0=2*index_i-mod(index_j,2)
				index_j0=floor_safe((index_j-1)/2d0)+1
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,modetrans,pid)
		elseif(collect=='B')then ! pair-wise broadcast
			if(mode=='R')then
				index_j0 = index_j+2*mod(index_j,2)-1
				index_i0=index_i
			elseif(mode=='C')then
				index_i0 = index_i+2*mod(index_i,2)-1
				index_j0=index_j
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,mode,pid)
		endif
		sendflag = pid/=ptree%MyID

		if(sendflag)then
			pp=pid-ptree%pgrp(blocks%pgno)%head+1
			if(recvquant(pp)%active==0)then
				recvquant(pp)%active=1
				Nrecvactive=Nrecvactive+1
				recvIDactive(Nrecvactive)=pp
			endif

			if(sendquant(pp)%active==0)then
				sendquant(pp)%active=1
				Nsendactive=Nsendactive+1
				sendIDactive(Nsendactive)=pp
			endif
			sendquant(pp)%size=sendquant(pp)%size+3+size(kerls%blocks(ii,jj)%matrix,1)*Ncol
		endif
	enddo
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif


	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,kerls%nr
	do jj=1,kerls%nc

		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c

		sendflag=.false.
		if(collect=='R')then ! pair-wise reduction
			if(mode=='R')then
				index_i0=floor_safe((index_i-1)/2d0)+1
				index_j0=2*index_j-mod(index_i,2)
			elseif(mode=='C')then
				index_i0=2*index_i-mod(index_j,2)
				index_j0=floor_safe((index_j-1)/2d0)+1
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,modetrans,pid)
		elseif(collect=='B')then ! pair-wise broadcast
			if(mode=='R')then
				index_j0 = index_j+2*mod(index_j,2)-1
				index_i0=index_i
			elseif(mode=='C')then
				index_i0 = index_i+2*mod(index_i,2)-1
				index_j0=index_j
			endif
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i0,index_j0,mode,pid)
		endif
		sendflag = pid/=ptree%MyID

		if(sendflag)then
			pp=pid-ptree%pgrp(blocks%pgno)%head+1
			Nrow=size(kerls%blocks(ii,jj)%matrix,1)
			sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
			sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
			sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
			sendquant(pp)%size=sendquant(pp)%size+3
			do i=1,Nrow*Ncol
				rr = mod(i-1,Nrow)+1
				cc = (i-1)/Nrow+1
				sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls%blocks(ii,jj)%matrix(rr,cc)
			enddo
			sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
		endif
	enddo
	enddo


	! communicate the data buffer
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo

	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			ii=(index_i-kerls%idx_r)/kerls%inc_r+1
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jj=(index_j-kerls%idx_c)/kerls%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			if(.not. allocated(kerls%blocks(ii,jj)%matrix))then
				allocate(kerls%blocks(ii,jj)%matrix(Nrow,Ncol))
				kerls%blocks(ii,jj)%matrix=0
			endif
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				kerls%blocks(ii,jj)%matrix(rr,cc) = kerls%blocks(ii,jj)%matrix(rr,cc) + recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_exchange_matvec



!*********** all to all communication of extraction results of one butterfly level from row-wise ordering to column-wise ordering or the reverse
subroutine BF_all2all_extraction(blocks,kerls,kerls1,stats,ptree,level,mode,mode_new)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i_loc_k,index_i_loc_s, index_j,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::blocks
	type(Hstat)::stats
	type(proctree)::ptree

    integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    integer,allocatable:: select_row_rr(:), select_column_rr(:)
    DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:),core(:,:),tau(:)

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Ncol,Nrow,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	character::mode,mode_new
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	type(butterfly_kerl)::kerls,kerls1
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist



	n1 = OMP_get_wtime()

	call assert(mode/=mode_new,'only row2col or col2row is supported')

	level_butterfly=blocks%level_butterfly
	nproc = ptree%pgrp(blocks%pgno)%nproc
	tag = blocks%pgno

	! mode_new and level_new determine the block range in the new mode
	if(mode_new=='R')then
		level_new=max(level-1,0)
	elseif(mode_new=='C')then
		level_new=min(level+1,level_butterfly+1)
	endif


	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	! calculate send buffer sizes in the first pass

	do nn=1,size(kerls1%index,1)
		ii = kerls1%index(nn,1)
		jj = kerls1%index(nn,2)
		index_i = (ii-1)*kerls1%inc_r+kerls1%idx_r
		index_j = (jj-1)*kerls1%inc_c+kerls1%idx_c
		call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(blocks%pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo


	do nn=1,size(kerls%index,1)
		ii = kerls%index(nn,1)
		jj = kerls%index(nn,2)
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c
		call GetBlockPID(ptree,blocks%pgno,level_new,level_butterfly,index_i,index_j,mode_new,pid)
		pp=pid-ptree%pgrp(blocks%pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(kerls%blocks(ii,jj)%matrix,1)*size(kerls%blocks(ii,jj)%matrix,2)
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do nn=1,size(kerls%index,1)
		ii = kerls%index(nn,1)
		jj = kerls%index(nn,2)
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c
		call GetBlockPID(ptree,blocks%pgno,level_new,level_butterfly,index_i,index_j,mode_new,pid)

		pp=pid-ptree%pgrp(blocks%pgno)%head+1
		Nrow=size(kerls%blocks(ii,jj)%matrix,1)
		Ncol=size(kerls%blocks(ii,jj)%matrix,2)

		sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
		sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
		sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
		sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
		sendquant(pp)%size=sendquant(pp)%size+4
		do i=1,Nrow*Ncol
			rr = mod(i-1,Nrow)+1
			cc = (i-1)/Nrow+1
			sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls%blocks(ii,jj)%matrix(rr,cc)
		enddo
		sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
		deallocate(kerls%blocks(ii,jj)%matrix)
		if(allocated(kerls%blocks(ii,jj)%index))deallocate(kerls%blocks(ii,jj)%index)
	enddo

	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif
	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(blocks%pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			ii=(index_i-kerls1%idx_r)/kerls1%inc_r+1
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jj=(index_j-kerls1%idx_c)/kerls1%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			call assert(.not. allocated(kerls1%blocks(ii,jj)%matrix),'receiving dat alreay exists locally')
			allocate(kerls1%blocks(ii,jj)%matrix(Nrow,Ncol))
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				kerls1%blocks(ii,jj)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_extraction





!*********** all to all communication of matvec results of one butterfly level from row-wise ordering to column-wise ordering or the reverse
subroutine BF_all2all_matvec(blocks,kerls,Ncol,stats,ptree,level,mode,mode_new)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i_loc_k,index_i_loc_s, index_j,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::blocks
	type(Hstat)::stats
	type(proctree)::ptree

    integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    integer,allocatable:: select_row_rr(:), select_column_rr(:)
    DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:),core(:,:),tau(:)

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Ncol,Nrow,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	character::mode,mode_new
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	type(butterfly_kerl)::kerls
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist



	n1 = OMP_get_wtime()

	call assert(mode/=mode_new,'only row2col or col2row is supported')

	level_butterfly=blocks%level_butterfly
	nproc = ptree%pgrp(blocks%pgno)%nproc
	tag = blocks%pgno

	! mode_new and level_new determine the block range in the new mode
	if(mode_new=='R')then
		level_new=max(level-1,0)
	elseif(mode_new=='C')then
		level_new=min(level+1,level_butterfly+1)
	endif
	call GetLocalBlockRange(ptree,blocks%pgno,level_new,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,mode_new)


	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	! calculate send buffer sizes in the first pass
	do ii=1,nr
	do jj=1,nc
		index_i = (ii-1)*inc_r+idx_r
		index_j = (jj-1)*inc_c+idx_c
		call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(blocks%pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo
	enddo


	do ii=1,kerls%nr
	do jj=1,kerls%nc
		index_i = (ii-1)*kerls%inc_r+kerls%idx_r
		index_j = (jj-1)*kerls%inc_c+kerls%idx_c
		call GetBlockPID(ptree,blocks%pgno,level_new,level_butterfly,index_i,index_j,mode_new,pid)
		pp=pid-ptree%pgrp(blocks%pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+3+size(kerls%blocks(ii,jj)%matrix,1)*Ncol
	enddo
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(blocks%pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(blocks%pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,kerls%nr
	do jj=1,kerls%nc
			index_i = (ii-1)*kerls%inc_r+kerls%idx_r
			index_j = (jj-1)*kerls%inc_c+kerls%idx_c
			call GetBlockPID(ptree,blocks%pgno,level_new,level_butterfly,index_i,index_j,mode_new,pid)

			pp=pid-ptree%pgrp(blocks%pgno)%head+1
			Nrow=size(kerls%blocks(ii,jj)%matrix,1)

			sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
			sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
			sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
			sendquant(pp)%size=sendquant(pp)%size+3
			do i=1,Nrow*Ncol
				rr = mod(i-1,Nrow)+1
				cc = (i-1)/Nrow+1
				sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls%blocks(ii,jj)%matrix(rr,cc)
			enddo
			sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
			deallocate(kerls%blocks(ii,jj)%matrix)
	enddo
	enddo
	deallocate(kerls%blocks)

	kerls%idx_r=idx_r
	kerls%idx_c=idx_c
	kerls%inc_r=inc_r
	kerls%inc_c=inc_c
	kerls%nr=nr
	kerls%nc=nc

	allocate(kerls%blocks(kerls%nr,kerls%nc))



	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(blocks%pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			ii=(index_i-kerls%idx_r)/kerls%inc_r+1
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jj=(index_j-kerls%idx_c)/kerls%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			call assert(.not. allocated(kerls%blocks(ii,jj)%matrix),'receiving dat alreay exists locally')
			allocate(kerls%blocks(ii,jj)%matrix(Nrow,Ncol))
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				kerls%blocks(ii,jj)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_matvec








!*********** all to all communication of one level of a butterfly from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
subroutine BF_all2all_ker(block_i,pgno_i,kerls_i,level_i,block_o,pgno_o,kerls_o,level_o,stats,ptree)

   use BPACK_DEFS
   implicit none
	integer pgno_i,pgno_o,pgno,level_i,level_o
    integer i, j, level_butterfly_i,level_butterfly_o, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,num_row,num_col,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	type(butterfly_kerl)::kerls_i,kerls_o
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist
	character::mode



	n1 = OMP_get_wtime()

	nproc = max(ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc)
	pgno = min(pgno_i,pgno_o)
	tag = pgno

	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_i=block_i%level_butterfly
	else
		level_butterfly_i=-1
		block_i%level_half=-1
	endif
	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o=block_o%level_butterfly
	else
		level_butterfly_o=-1
		block_o%level_half=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_i,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_o%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)


	if(level_i<=block_i%level_half)then
		call assert(level_o<=block_o%level_half,'row-wise ordering is only redistributed to row-wise ordering')
		mode='R'
	endif

	if(level_i>block_i%level_half)then
		call assert(level_o>block_o%level_half,'column-wise ordering is only redistributed to column-wise ordering')
		mode='C'
	endif

	call assert((ptree%pgrp(pgno_i)%head<=ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail>=ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head<=ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail>=ptree%pgrp(pgno_i)%tail),'pgno_i or pgno_o should be contained in the other')


	call GetLocalBlockRange(ptree,pgno_o,level_o,level_butterfly_o,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)

	num_row=2**level_o
	num_col=2**(level_butterfly_o-level_o+1)

	if(mode=='R')then
		idx_c=idx_c*2-1
		inc_c=inc_c
		nc=nc*2
	elseif(mode=='C')then
		idx_r=idx_r*2-1
		inc_r=inc_r
		nr=nr*2
	endif

	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	! calculate send buffer sizes in the first pass
	do ii=1,nr
	do jj=1,nc
		index_i = (ii-1)*inc_r+idx_r
		index_j = (jj-1)*inc_c+idx_c
		if(mode=='R')then
			index_j0=floor_safe((index_j-1)/2d0)+1
			index_i0=index_i
		endif
		if(mode=='C')then
			index_i0=floor_safe((index_i-1)/2d0)+1
			index_j0=index_j
		endif
		call GetBlockPID(ptree,pgno_i,level_i,level_butterfly_i,index_i0,index_j0,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo
	enddo


	do ii=1,kerls_i%nr
	do jj=1,kerls_i%nc
		index_i = (ii-1)*kerls_i%inc_r+kerls_i%idx_r
		index_j = (jj-1)*kerls_i%inc_c+kerls_i%idx_c
		if(mode=='R')then
			index_j0=floor_safe((index_j-1)/2d0)+1
			index_i0=index_i
		endif
		if(mode=='C')then
			index_i0=floor_safe((index_i-1)/2d0)+1
			index_j0=index_j
		endif
		call GetBlockPID(ptree,pgno_o,level_o,level_butterfly_o,index_i0,index_j0,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(kerls_i%blocks(ii,jj)%matrix,1)*size(kerls_i%blocks(ii,jj)%matrix,2)
	enddo
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,kerls_i%nr
	do jj=1,kerls_i%nc
			index_i = (ii-1)*kerls_i%inc_r+kerls_i%idx_r
			index_j = (jj-1)*kerls_i%inc_c+kerls_i%idx_c

			if(mode=='R')then
				index_j0=floor_safe((index_j-1)/2d0)+1
				index_i0=index_i
			endif
			if(mode=='C')then
				index_i0=floor_safe((index_i-1)/2d0)+1
				index_j0=index_j
			endif

			call GetBlockPID(ptree,pgno_o,level_o,level_butterfly_o,index_i0,index_j0,mode,pid)

			pp=pid-ptree%pgrp(pgno)%head+1
			Nrow=size(kerls_i%blocks(ii,jj)%matrix,1)
			Ncol=size(kerls_i%blocks(ii,jj)%matrix,2)

			sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
			sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
			sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
			sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
			sendquant(pp)%size=sendquant(pp)%size+4
			do i=1,Nrow*Ncol
				rr = mod(i-1,Nrow)+1
				cc = (i-1)/Nrow+1
				sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls_i%blocks(ii,jj)%matrix(rr,cc)
			enddo
			sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
			deallocate(kerls_i%blocks(ii,jj)%matrix)
	enddo
	enddo
	if(allocated(kerls_i%blocks))deallocate(kerls_i%blocks)

	if(nr>0 .and. nc>0)then
	kerls_o%idx_r=idx_r
	kerls_o%idx_c=idx_c
	kerls_o%inc_r=inc_r
	kerls_o%inc_c=inc_c
	kerls_o%nr=nr
	kerls_o%nc=nc
	kerls_o%num_row=num_row
	kerls_o%num_col=num_col
	allocate(kerls_o%blocks(kerls_o%nr,kerls_o%nc))
	endif


	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			ii=(index_i-kerls_o%idx_r)/kerls_o%inc_r+1
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jj=(index_j-kerls_o%idx_c)/kerls_o%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			call assert(.not. allocated(kerls_o%blocks(ii,jj)%matrix),'receiving dat alreay exists locally')
			allocate(kerls_o%blocks(ii,jj)%matrix(Nrow,Ncol))
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				kerls_o%blocks(ii,jj)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_ker



!*********** convert blocks in block_i%sons to block_o%sons, this is a local function without MPI communication, it is assumed block_i%sons has L levels, and block_o%sons will have max(L-2,0) levels
subroutine BF_convert_to_smallBF(block_i,block_o,stats,ptree)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly_i,level_butterfly_o, level_butterfly_c_o,num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_ic, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_jc,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt, iii, jjj
	integer mm1,nn1,mm2,nn2,M1,N1,kk
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(matrixblock),pointer::block_c_i,block_c_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,num_row,num_col,level_new,idx

	real(kind=8)::t1,t2
	DT,allocatable::matrixtemp1(:,:),matrixtemp2(:,:)

	if(IOwnPgrp(ptree,block_o%sons(1,1)%pgno))then


	t1 = OMP_get_wtime()

	do iii=1,2
	do jjj=1,2
		block_c_o=>block_o%sons(iii,jjj)
		block_c_i=>block_i%sons(iii,jjj)
		if(block_i%level_butterfly==1)then  ! l-level butterfly becomes 0-level butterfly
			block_c_o%level_butterfly=0
			block_c_o%level_half = 0
			allocate(block_c_o%ButterflyU%blocks(1))
			allocate(block_c_o%ButterflyV%blocks(1))
			block_c_o%ButterflyU%nblk_loc=1
			block_c_o%ButterflyU%inc=1
			block_c_o%ButterflyU%idx=1
			block_c_o%ButterflyV%nblk_loc=1
			block_c_o%ButterflyV%inc=1
			block_c_o%ButterflyV%idx=1

			call assert(block_c_i%ButterflyU%nblk_loc==1,'parent block has more than one ButterflyU block')
			call assert(block_c_i%ButterflyV%nblk_loc==1,'parent block has more than one ButterflyV block')
			call assert(block_c_i%ButterflyKerl(1)%nr==1 .and. block_c_i%ButterflyKerl(1)%nc==1,'parent block has more than one ButterflyKerl block')

			mm1 = size(block_c_i%ButterflyKerl(1)%blocks(1,1)%matrix,1)
			nn1 = size(block_c_i%ButterflyKerl(1)%blocks(1,1)%matrix,2)
			M1 = size(block_c_i%ButterflyU%blocks(1)%matrix,1)
			N1 = size(block_c_i%ButterflyV%blocks(1)%matrix,1)

			allocate(block_c_o%ButterflyU%blocks(1)%matrix(M1,nn1))
			allocate(block_c_o%ButterflyV%blocks(1)%matrix(N1,nn1))
			call gemmf90(block_c_i%ButterflyU%blocks(1)%matrix,M1,block_c_i%ButterflyKerl(1)%blocks(1,1)%matrix,mm1,block_c_o%ButterflyU%blocks(1)%matrix,M1,'N','N',M1, nn1, mm1,cone,czero)
			block_c_o%ButterflyV%blocks(1)%matrix = block_c_i%ButterflyV%blocks(1)%matrix
		else ! L-level butterfly becomes (L-2)-level butterfly
			block_c_o%level_butterfly=block_i%level_butterfly-2
			block_c_o%level_half = floor_safe(dble(block_c_o%level_butterfly)/2d0) ! from outer to inner
			call assert(block_c_o%level_butterfly>=0,'negative level_butterfly!')
			if(block_c_o%level_butterfly>0)then
				allocate(block_c_o%ButterflyKerl(block_c_o%level_butterfly))
			endif
			do level=0,block_c_o%level_butterfly+1
				if(level==0)then
					block_c_o%ButterflyV%num_blk=2**block_c_o%level_butterfly
					call assert(mod(block_c_i%ButterflyV%nblk_loc,2)==0,'parent block should have even number of ButterflyV blocks')
					block_c_o%ButterflyV%nblk_loc=block_c_i%ButterflyV%nblk_loc/2
					block_c_o%ButterflyV%inc=block_c_i%ButterflyV%inc
					idx =block_c_i%ButterflyV%idx
					if(idx>block_c_o%ButterflyV%num_blk*2)idx=idx-block_c_o%ButterflyV%num_blk*2
					block_c_o%ButterflyV%idx=ceiling_safe(idx/2d0)
					allocate(block_c_o%ButterflyV%blocks(block_c_o%ButterflyV%nblk_loc))
					do ii =1,block_c_o%ButterflyV%nblk_loc
						mm1 = size(block_c_i%ButterflyV%blocks(2*ii-1)%matrix,1)
						nn1 = size(block_c_i%ButterflyV%blocks(2*ii-1)%matrix,2)
						mm2 = size(block_c_i%ButterflyV%blocks(2*ii)%matrix,1)
						nn2 = size(block_c_i%ButterflyV%blocks(2*ii)%matrix,2)
						kk = size(block_c_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,1)
						allocate(block_c_o%ButterflyV%blocks(ii)%matrix(mm1+mm2,kk))
						N1=N1+mm1+mm2
						allocate(matrixtemp1(mm1,kk))
						allocate(matrixtemp2(mm2,kk))
						call gemmf90(block_c_i%ButterflyV%blocks(2*ii-1)%matrix,mm1, block_c_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,kk, matrixtemp1,mm1, 'N','T',mm1,kk,nn1,cone,czero)
						call gemmf90(block_c_i%ButterflyV%blocks(2*ii)%matrix,mm2, block_c_i%ButterflyKerl(1)%blocks(1,2*ii)%matrix,kk, matrixtemp2,mm2, 'N','T',mm2,kk,nn2,cone,czero)
						block_c_o%ButterflyV%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
						block_c_o%ButterflyV%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
						deallocate(matrixtemp1)
						deallocate(matrixtemp2)
					enddo
				elseif(level==block_c_o%level_butterfly+1)then
					block_c_o%ButterflyU%num_blk=2**block_c_o%level_butterfly
					call assert(mod(block_c_i%ButterflyU%nblk_loc,2)==0,'parent block should have even number of ButterflyU blocks')
					block_c_o%ButterflyU%nblk_loc=block_c_i%ButterflyU%nblk_loc/2
					block_c_o%ButterflyU%inc=block_c_i%ButterflyU%inc
					idx =block_c_i%ButterflyU%idx
					if(idx>block_c_o%ButterflyU%num_blk*2)idx=idx-block_c_o%ButterflyU%num_blk*2
					block_c_o%ButterflyU%idx=ceiling_safe(idx/2d0)
					allocate(block_c_o%ButterflyU%blocks(block_c_o%ButterflyU%nblk_loc))
					do ii =1,block_c_o%ButterflyU%nblk_loc
						mm1 = size(block_c_i%ButterflyU%blocks(2*ii-1)%matrix,1)
						nn1 = size(block_c_i%ButterflyU%blocks(2*ii-1)%matrix,2)
						mm2 = size(block_c_i%ButterflyU%blocks(2*ii)%matrix,1)
						nn2 = size(block_c_i%ButterflyU%blocks(2*ii)%matrix,2)
						kk = size(block_c_i%ButterflyKerl(block_c_o%level_butterfly+2)%blocks(2*ii-1,1)%matrix,2)
						allocate(block_c_o%ButterflyU%blocks(ii)%matrix(mm1+mm2,kk))
						M1=M1+mm1+mm2
						allocate(matrixtemp1(mm1,kk))
						allocate(matrixtemp2(mm2,kk))
						call gemmf90(block_c_i%ButterflyU%blocks(2*ii-1)%matrix,mm1,block_c_i%ButterflyKerl(block_c_i%level_butterfly)%blocks(2*ii-1,1)%matrix,nn1,matrixtemp1,mm1,'N','N',mm1,kk,nn1,cone,czero)
						call gemmf90(block_c_i%ButterflyU%blocks(2*ii)%matrix,mm2,block_c_i%ButterflyKerl(block_c_i%level_butterfly)%blocks(2*ii,1)%matrix,nn2,matrixtemp2,mm2,'N','N',mm2,kk,nn2,cone,czero)
						block_c_o%ButterflyU%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
						block_c_o%ButterflyU%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
						deallocate(matrixtemp1)
						deallocate(matrixtemp2)
					end do
				else
					num_col=block_c_i%ButterflyKerl(level+1)%num_col
					num_row=block_c_i%ButterflyKerl(level+1)%num_row
					block_c_o%ButterflyKerl(level)%num_row=num_row/2
					block_c_o%ButterflyKerl(level)%num_col=num_col/2
					block_c_o%ButterflyKerl(level)%nr=block_c_i%ButterflyKerl(level+1)%nr
					block_c_o%ButterflyKerl(level)%inc_r=block_c_i%ButterflyKerl(level+1)%inc_r
					block_c_o%ButterflyKerl(level)%idx_r=block_c_i%ButterflyKerl(level+1)%idx_r
					if( block_c_o%ButterflyKerl(level)%idx_r> block_c_o%ButterflyKerl(level)%num_row)block_c_o%ButterflyKerl(level)%idx_r=block_c_o%ButterflyKerl(level)%idx_r- block_c_o%ButterflyKerl(level)%num_row
					block_c_o%ButterflyKerl(level)%nc=block_c_i%ButterflyKerl(level+1)%nc
					block_c_o%ButterflyKerl(level)%inc_c=block_c_i%ButterflyKerl(level+1)%inc_c
					block_c_o%ButterflyKerl(level)%idx_c=block_c_i%ButterflyKerl(level+1)%idx_c
					if( block_c_o%ButterflyKerl(level)%idx_c> block_c_o%ButterflyKerl(level)%num_col)block_c_o%ButterflyKerl(level)%idx_c=block_c_o%ButterflyKerl(level)%idx_c- block_c_o%ButterflyKerl(level)%num_col
					allocate(block_c_o%ButterflyKerl(level)%blocks(block_c_o%ButterflyKerl(level)%nr,block_c_o%ButterflyKerl(level)%nc))

					do ii=1,block_c_o%ButterflyKerl(level)%nr
					do jj=1,block_c_o%ButterflyKerl(level)%nc
						mm=size(block_c_i%ButterflyKerl(level+1)%blocks(ii,jj)%matrix,1)
						nn=size(block_c_i%ButterflyKerl(level+1)%blocks(ii,jj)%matrix,2)
						allocate(block_c_o%ButterflyKerl(level)%blocks(ii,jj)%matrix(mm,nn))
						block_c_o%ButterflyKerl(level)%blocks(ii,jj)%matrix = block_c_i%ButterflyKerl(level+1)%blocks(ii,jj)%matrix
					enddo
					enddo
				endif
			enddo
		endif
	enddo
	enddo

	t2 = OMP_get_wtime()
	! time_tmp = time_tmp + t2 - t1
	endif
end subroutine BF_convert_to_smallBF




















!*********** all to all communication of one level of a butterfly into four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
subroutine BF_all2all_ker_split(block_i,pgno_i,level_i,block_o,pgno_o,level_o,stats,ptree)

   use BPACK_DEFS
   implicit none
	integer pgno_i,pgno_o,pgno,level_i,level_o
    integer i, j, level_butterfly_i,level_butterfly_o, level_butterfly_c_o,num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_ic, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_jc,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt, iii, jjj
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,num_row,num_col,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	! type(butterfly_kerl)::kerls_i,kerls_o
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist
	character::mode



	n1 = OMP_get_wtime()

	nproc = max(ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc)
	pgno = min(pgno_i,pgno_o)
	tag = pgno

	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_i=block_i%level_butterfly
	else
		level_butterfly_i=-1
		block_i%level_half=-1
	endif
	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o=block_o%level_butterfly
	else
		level_butterfly_o=-1
		block_o%level_half=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_i,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_o%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)

	if(level_i<=block_i%level_half)then
		call assert(level_o<=block_o%level_half,'row-wise ordering is only redistributed to row-wise ordering')
		mode='R'
	endif

	if(level_i>block_i%level_half)then
		call assert(level_o>block_o%level_half,'column-wise ordering is only redistributed to column-wise ordering')
		mode='C'
	endif

	call assert((ptree%pgrp(pgno_i)%head<=ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail>=ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head<=ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail>=ptree%pgrp(pgno_i)%tail),'pgno_i or pgno_o should be contained in the other')


	num_row=2**level_o
	num_col=2**(level_butterfly_o-level_o+1)
	level_butterfly_c_o = max(level_butterfly_o-2,0)
	call GetLocalBlockRange(ptree,pgno_o,level_o-1,level_butterfly_c_o,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)
	if(mode=='R')then
		idx_c=idx_c*2-1
		inc_c=inc_c
		if(level_butterfly_o>1)nc=nc*2
	elseif(mode=='C')then
		idx_r=idx_r*2-1
		inc_r=inc_r
		if(level_butterfly_o>1)nr=nr*2
	endif

	do iii=1,2
	do jjj=1,2
		if(nr>0 .and. nc>0)then
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_r=idx_r+(iii-1)*num_row/2
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_c=idx_c+(jjj-1)*num_col/2
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_r=inc_r
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_c=inc_c
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nr=nr
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nc=nc
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%num_row=num_row
		block_o%sons(iii,jjj)%ButterflyKerl(level_o)%num_col=num_col
		allocate(block_o%sons(iii,jjj)%ButterflyKerl(level_o)%blocks(block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nr,block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nc))
		endif
	enddo
	enddo


	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0

	! calculate send buffer sizes in the first pass
	do iii=1,2
	do jjj=1,2
	do ii=1,block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nr
	do jj=1,block_o%sons(iii,jjj)%ButterflyKerl(level_o)%nc
		index_i = (ii-1)*block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_r+block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_r
		index_j = (jj-1)*block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_c+block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_c
		if(mode=='R')then
			index_j0=floor_safe((index_j-1)/2d0)+1
			index_i0=index_i
		endif
		if(mode=='C')then
			index_i0=floor_safe((index_i-1)/2d0)+1
			index_j0=index_j
		endif
		call GetBlockPID(ptree,pgno_i,level_i,level_butterfly_i,index_i0,index_j0,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo
	enddo
	enddo
	enddo

	do ii=1,block_i%ButterflyKerl(level_i)%nr
	do jj=1,block_i%ButterflyKerl(level_i)%nc
		index_i = (ii-1)*block_i%ButterflyKerl(level_i)%inc_r+block_i%ButterflyKerl(level_i)%idx_r
		index_j = (jj-1)*block_i%ButterflyKerl(level_i)%inc_c+block_i%ButterflyKerl(level_i)%idx_c
		index_ic = index_i
		index_jc = index_j
		if(index_ic>num_row/2)index_ic=index_ic-num_row/2
		if(index_jc>num_col/2)index_jc=index_jc-num_col/2


		if(mode=='R')then
			index_j0=floor_safe((index_jc-1)/2d0)+1
			index_i0=index_ic
		endif
		if(mode=='C')then
			index_i0=floor_safe((index_ic-1)/2d0)+1
			index_j0=index_jc
		endif
		call GetBlockPID(ptree,pgno_o,level_o-1,level_butterfly_c_o,index_i0,index_j0,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix,1)*size(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix,2)
	enddo
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,block_i%ButterflyKerl(level_i)%nr
	do jj=1,block_i%ButterflyKerl(level_i)%nc
			index_i = (ii-1)*block_i%ButterflyKerl(level_i)%inc_r+block_i%ButterflyKerl(level_i)%idx_r
			index_j = (jj-1)*block_i%ButterflyKerl(level_i)%inc_c+block_i%ButterflyKerl(level_i)%idx_c

			index_ic = index_i
			index_jc = index_j
			if(index_ic>num_row/2)index_ic=index_ic-num_row/2
			if(index_jc>num_col/2)index_jc=index_jc-num_col/2


			if(mode=='R')then
				index_j0=floor_safe((index_jc-1)/2d0)+1
				index_i0=index_ic
			endif
			if(mode=='C')then
				index_i0=floor_safe((index_ic-1)/2d0)+1
				index_j0=index_jc
			endif
			call GetBlockPID(ptree,pgno_o,level_o-1,level_butterfly_c_o,index_i0,index_j0,mode,pid)

			pp=pid-ptree%pgrp(pgno)%head+1
			Nrow=size(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix,1)
			Ncol=size(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix,2)

			sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
			sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
			sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
			sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
			sendquant(pp)%size=sendquant(pp)%size+4
			do i=1,Nrow*Ncol
				rr = mod(i-1,Nrow)+1
				cc = (i-1)/Nrow+1
				sendquant(pp)%dat(sendquant(pp)%size+i,1) = block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix(rr,cc)
			enddo
			sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
			! deallocate(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix)
	enddo
	enddo
	! if(allocated(block_i%ButterflyKerl(level_i)%blocks))deallocate(block_i%ButterflyKerl(level_i)%blocks)


	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			iii=1
			if(index_i>num_row/2)iii=2
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))
			jjj=1
			if(index_j>num_col/2)jjj=2

			ii=(index_i-block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_r)/block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_r+1
			jj=(index_j-block_o%sons(iii,jjj)%ButterflyKerl(level_o)%idx_c)/block_o%sons(iii,jjj)%ButterflyKerl(level_o)%inc_c+1
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			call assert(.not. allocated(block_o%sons(iii,jjj)%ButterflyKerl(level_o)%blocks(ii,jj)%matrix),'receiving dat alreay exists locally')
			allocate(block_o%sons(iii,jjj)%ButterflyKerl(level_o)%blocks(ii,jj)%matrix(Nrow,Ncol))
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				block_o%sons(iii,jjj)%ButterflyKerl(level_o)%blocks(ii,jj)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_ker_split





!*********** all to all communication of one level of a butterfly from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
subroutine BF_all2all_UV(block_i,pgno_i,kerls_i,level_i,block_o,pgno_o,kerls_o,level_o,stats,ptree)

   use BPACK_DEFS
   implicit none
	integer pgno_i,pgno_o,pgno,level_i,level_o
    integer i, j, level_butterfly_i,level_butterfly_o, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,idx,inc_r,inc_c,inc,nr,nc,nblk_loc,num_blk,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	type(butterfly_UV)::kerls_i,kerls_o
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist
	character::mode

	n1 = OMP_get_wtime()

	nproc = max(ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc)
	pgno = min(pgno_i,pgno_o)
	tag = pgno


	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_i=block_i%level_butterfly
	else
		level_butterfly_i=-1
		block_i%level_half=-1
	endif
	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o=block_o%level_butterfly
	else
		level_butterfly_o=-1
		block_o%level_half=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_i,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_o%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)

	if(level_i<=block_i%level_half)then
		call assert(level_o<=block_o%level_half,'row-wise ordering is only redistributed to row-wise ordering')
		mode='R'
	endif

	if(level_i>block_i%level_half)then
		call assert(level_o>block_o%level_half,'column-wise ordering is only redistributed to column-wise ordering')
		mode='C'
	endif

	call assert((ptree%pgrp(pgno_i)%head<=ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail>=ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head<=ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail>=ptree%pgrp(pgno_i)%tail),'pgno_i or pgno_o should be contained in the other')


	call GetLocalBlockRange(ptree,pgno_o,level_o,level_butterfly_o,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)

	if(mode=='R')then
		num_blk=2**level_butterfly_o
		idx=idx_c
		inc=inc_c
		nblk_loc=nc
	elseif(mode=='C')then
		num_blk=2**level_butterfly_o
		idx=idx_r
		inc=inc_r
		nblk_loc=nr
	endif

	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	! calculate send buffer sizes in the first pass
	do ii=1,nblk_loc
		! convert indices from output to input
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*inc+idx
		elseif(mode=='C')then
			index_i = (ii-1)*inc+idx
			index_j = 1
		endif
		call GetBlockPID(ptree,pgno_i,level_i,level_butterfly_i,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo

	do ii=1,kerls_i%nblk_loc
		! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*kerls_i%inc+kerls_i%idx
		elseif(mode=='C')then
			index_i = (ii-1)*kerls_i%inc+kerls_i%idx
			index_j = 1
		endif
		call GetBlockPID(ptree,pgno_o,level_o,level_butterfly_o,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(kerls_i%blocks(ii)%matrix,1)*size(kerls_i%blocks(ii)%matrix,2)
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,kerls_i%nblk_loc
	   ! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*kerls_i%inc+kerls_i%idx
		elseif(mode=='C')then
			index_i = (ii-1)*kerls_i%inc+kerls_i%idx
			index_j = 1
		endif

		call GetBlockPID(ptree,pgno_o,level_o,level_butterfly_o,index_i,index_j,mode,pid)

		pp=pid-ptree%pgrp(pgno)%head+1
		Nrow=size(kerls_i%blocks(ii)%matrix,1)
		Ncol=size(kerls_i%blocks(ii)%matrix,2)

		sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
		sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
		sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
		sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
		sendquant(pp)%size=sendquant(pp)%size+4
		do i=1,Nrow*Ncol
			rr = mod(i-1,Nrow)+1
			cc = (i-1)/Nrow+1
			sendquant(pp)%dat(sendquant(pp)%size+i,1) = kerls_i%blocks(ii)%matrix(rr,cc)
		enddo
		sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
		deallocate(kerls_i%blocks(ii)%matrix)
	enddo
	if(allocated(kerls_i%blocks))deallocate(kerls_i%blocks)

	if(nblk_loc>0)then
	kerls_o%idx=idx
	kerls_o%inc=inc
	kerls_o%nblk_loc=nblk_loc
	kerls_o%num_blk=num_blk
	allocate(kerls_o%blocks(kerls_o%nblk_loc))
	endif

	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))

			if(mode=='R')then
				ii=(index_j-kerls_o%idx)/kerls_o%inc+1
			elseif(mode=='C')then
				ii=(index_i-kerls_o%idx)/kerls_o%inc+1
			endif
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			call assert(.not. allocated(kerls_o%blocks(ii)%matrix),'receiving dat alreay exists locally')
			allocate(kerls_o%blocks(ii)%matrix(Nrow,Ncol))
			do j=1,Nrow*Ncol
				rr = mod(j-1,Nrow)+1
				cc = (j-1)/Nrow+1
				kerls_o%blocks(ii)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_UV





!*********** all to all communication of one level of a butterfly to four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
subroutine BF_all2all_U_split(block_i,pgno_i,level_i,block_o,pgno_o,level_o,stats,ptree)

   use BPACK_DEFS
   implicit none
	integer pgno_i,pgno_o,pgno,level_i,level_o,level_c_o
    integer i, j, level_butterfly_i,level_butterfly_o, level_butterfly_c_o,num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1,iii,jjj,si,sj,ni,nj
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,idx,inc_r,inc_c,inc,nr,nc,nblk_loc,num_blk,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	! type(butterfly_UV)::kerls_i,kerls_o
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist
	character::mode

	n1 = OMP_get_wtime()

	nproc = max(ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc)
	pgno = min(pgno_i,pgno_o)
	tag = pgno


	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_i=block_i%level_butterfly
	else
		level_butterfly_i=-1
		block_i%level_half=-1
	endif
	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o=block_i%level_butterfly
	else
		level_butterfly_o=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_i,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)


	level_butterfly_c_o = max(level_butterfly_o-1,0)
	if(level_i<=block_i%level_half)then
		mode='R'
		level_c_o=0
	endif

	if(level_i>block_i%level_half)then
		mode='C'
		level_c_o=level_butterfly_c_o+1
	endif

	call assert((ptree%pgrp(pgno_i)%head<=ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail>=ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head<=ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail>=ptree%pgrp(pgno_i)%tail),'pgno_i or pgno_o should be contained in the other')


	call GetLocalBlockRange(ptree,pgno_o,level_c_o,level_butterfly_c_o,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)

	if(mode=='R')then
		num_blk=2**level_butterfly_o
		idx=idx_c
		inc=inc_c
		nblk_loc=nc
	elseif(mode=='C')then
		num_blk=2**level_butterfly_o
		idx=idx_r
		inc=inc_r
		nblk_loc=nr
	endif


	do iii=1,2
	do jjj=1,2
		if(nblk_loc>0)then
		block_o%sons(iii,jjj)%ButterflyU%idx=idx+(iii-1)*num_blk/2
		block_o%sons(iii,jjj)%ButterflyU%inc=inc
		block_o%sons(iii,jjj)%ButterflyU%nblk_loc=nblk_loc
		block_o%sons(iii,jjj)%ButterflyU%num_blk=num_blk
		allocate(block_o%sons(iii,jjj)%ButterflyU%blocks(block_o%sons(iii,jjj)%ButterflyU%nblk_loc))
		endif
	enddo
	enddo


	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	if(mode=='R')then
		nj=2
		ni=1
	elseif(mode=='C')then
		nj=1
		ni=2
	endif

	! calculate send buffer sizes in the first pass
	do iii=1,ni
	do jjj=1,nj
	do ii=1,nblk_loc
		! convert indices from output to input
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_o%sons(iii,jjj)%ButterflyU%inc+block_o%sons(iii,jjj)%ButterflyU%idx
		elseif(mode=='C')then
			index_i = (ii-1)*block_o%sons(iii,jjj)%ButterflyU%inc+block_o%sons(iii,jjj)%ButterflyU%idx
			index_j = 1
		endif
		call GetBlockPID(ptree,pgno_i,level_i,level_butterfly_i,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo
	enddo
	enddo

	do ii=1,block_i%ButterflyU%nblk_loc
		! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
			index_i0 = 1
			index_j0 = index_j
			if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
		elseif(mode=='C')then
			index_i = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
			index_j = 1
			index_j0 = 1
			index_i0 = index_i
			if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
		endif

		call GetBlockPID(ptree,pgno_o,level_c_o,level_butterfly_c_o,index_i0,index_j0,mode,pid)

		pp=pid-ptree%pgrp(pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(block_i%ButterflyU%blocks(ii)%matrix,1)*size(block_i%ButterflyU%blocks(ii)%matrix,2)
	enddo

	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo


	! pack the send buffer in the second pass
	do ii=1,block_i%ButterflyU%nblk_loc
	   ! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
			index_i0 = 1
			index_j0 = index_j
			if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
		elseif(mode=='C')then
			index_i = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
			index_j = 1
			index_j0 = 1
			index_i0 = index_i
			if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
		endif

		call GetBlockPID(ptree,pgno_o,level_c_o,level_butterfly_c_o,index_i0,index_j0,mode,pid)

		pp=pid-ptree%pgrp(pgno)%head+1
		Nrow=size(block_i%ButterflyU%blocks(ii)%matrix,1)
		Ncol=size(block_i%ButterflyU%blocks(ii)%matrix,2)

		sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
		sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
		sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
		sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
		sendquant(pp)%size=sendquant(pp)%size+4
		do i=1,Nrow*Ncol
			rr = mod(i-1,Nrow)+1
			cc = (i-1)/Nrow+1
			sendquant(pp)%dat(sendquant(pp)%size+i,1) = block_i%ButterflyU%blocks(ii)%matrix(rr,cc)
		enddo
		sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
		! deallocate(block_i%ButterflyU%blocks(ii)%matrix)
	enddo
	! if(allocated(block_i%ButterflyU%blocks))deallocate(block_i%ButterflyU%blocks)



	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))

			if(mode=='R')then
				nj=1
				sj=1
				if(index_j>num_blk/2)sj=2
				ni=2
				si=1
				ii=(index_j-block_o%sons(si,sj)%ButterflyU%idx)/block_o%sons(si,sj)%ButterflyU%inc+1
			elseif(mode=='C')then
				ni=1
				si=1
				if(index_i>num_blk/2)si=2
				nj=2
				sj=1
				ii=(index_i-block_o%sons(si,sj)%ButterflyU%idx)/block_o%sons(si,sj)%ButterflyU%inc+1
			endif
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			do iii=si,si+ni-1
			do jjj=sj,sj+nj-1
				call assert(.not. allocated(block_o%sons(iii,jjj)%ButterflyU%blocks(ii)%matrix),'receiving dat alreay exists locally')
				allocate(block_o%sons(iii,jjj)%ButterflyU%blocks(ii)%matrix(Nrow,Ncol))
				do j=1,Nrow*Ncol
					rr = mod(j-1,Nrow)+1
					cc = (j-1)/Nrow+1
					block_o%sons(iii,jjj)%ButterflyU%blocks(ii)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
				enddo
			enddo
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_U_split



!*********** all to all communication of one level of a butterfly to four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
subroutine BF_all2all_V_split(block_i,pgno_i,level_i,block_o,pgno_o,level_o,stats,ptree)

   use BPACK_DEFS
   implicit none
	integer pgno_i,pgno_o,pgno,level_i,level_o,level_c_o
    integer i, j, level_butterfly_i,level_butterfly_o, level_butterfly_c_o,num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k,index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1,iii,jjj,si,sj,ni,nj
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::block_i,block_o
	type(Hstat)::stats
	type(proctree)::ptree

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Ncol,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,idx,inc_r,inc_c,inc,nr,nc,nblk_loc,num_blk,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	real(kind=8)::n1,n2
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	! type(butterfly_UV)::kerls_i,kerls_o
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist
	character::mode

	n1 = OMP_get_wtime()

	nproc = max(ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc)
	pgno = min(pgno_i,pgno_o)
	tag = pgno


	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_i=block_i%level_butterfly
	else
		level_butterfly_i=-1
		block_i%level_half=-1
	endif
	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o=block_i%level_butterfly
	else
		level_butterfly_o=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_i,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%level_half,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)


	level_butterfly_c_o = max(level_butterfly_o-1,0)
	if(level_i<=block_i%level_half)then
		mode='R'
		level_c_o=0
	endif

	if(level_i>block_i%level_half)then
		mode='C'
		level_c_o=level_butterfly_c_o+1
	endif

	call assert((ptree%pgrp(pgno_i)%head<=ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail>=ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head<=ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail>=ptree%pgrp(pgno_i)%tail),'pgno_i or pgno_o should be contained in the other')


	call GetLocalBlockRange(ptree,pgno_o,level_c_o,level_butterfly_c_o,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)

	if(mode=='R')then
		num_blk=2**level_butterfly_o
		idx=idx_c
		inc=inc_c
		nblk_loc=nc
	elseif(mode=='C')then
		num_blk=2**level_butterfly_o
		idx=idx_r
		inc=inc_r
		nblk_loc=nr
	endif

	do iii=1,2
	do jjj=1,2
		if(nblk_loc>0)then
		block_o%sons(iii,jjj)%ButterflyV%idx=idx+(jjj-1)*num_blk/2
		block_o%sons(iii,jjj)%ButterflyV%inc=inc
		block_o%sons(iii,jjj)%ButterflyV%nblk_loc=nblk_loc
		block_o%sons(iii,jjj)%ButterflyV%num_blk=num_blk
		allocate(block_o%sons(iii,jjj)%ButterflyV%blocks(block_o%sons(iii,jjj)%ButterflyV%nblk_loc))
		endif
	enddo
	enddo

	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	if(mode=='R')then
		nj=2
		ni=1
	elseif(mode=='C')then
		nj=1
		ni=2
	endif
	! calculate send buffer sizes in the first pass
	do iii=1,ni
	do jjj=1,nj
	do ii=1,nblk_loc
		! convert indices from output to input
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_o%sons(iii,jjj)%ButterflyV%inc+block_o%sons(iii,jjj)%ButterflyV%idx
		elseif(mode=='C')then
			index_i = (ii-1)*block_o%sons(iii,jjj)%ButterflyV%inc+block_o%sons(iii,jjj)%ButterflyV%idx
			index_j = 1
		endif
		call GetBlockPID(ptree,pgno_i,level_i,level_butterfly_i,index_i,index_j,mode,pid)
		pp=pid-ptree%pgrp(pgno)%head+1
		if(recvquant(pp)%active==0)then
			recvquant(pp)%active=1
			Nrecvactive=Nrecvactive+1
			recvIDactive(Nrecvactive)=pp
		endif
	enddo
	enddo
	enddo
	do ii=1,block_i%ButterflyV%nblk_loc
		! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
			index_i0 = 1
			index_j0 = index_j
			if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
		elseif(mode=='C')then
			index_i = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
			index_j = 1
			index_j0 = 1
			index_i0 = index_i
			if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
		endif
		call GetBlockPID(ptree,pgno_o,level_c_o,level_butterfly_c_o,index_i0,index_j0,mode,pid)

		pp=pid-ptree%pgrp(pgno)%head+1
		if(sendquant(pp)%active==0)then
			sendquant(pp)%active=1
			Nsendactive=Nsendactive+1
			sendIDactive(Nsendactive)=pp
		endif
		sendquant(pp)%size=sendquant(pp)%size+4+size(block_i%ButterflyV%blocks(ii)%matrix,1)*size(block_i%ButterflyV%blocks(ii)%matrix,2)
	enddo
	! communicate receive buffer sizes
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		recvid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Isend(sendquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,S_req(tt),ierr)
	enddo

	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		sendid=pp-1+ptree%pgrp(pgno)%head
		call MPI_Irecv(recvquant(pp)%size,1,MPI_INTEGER,pp-1,tag,ptree%pgrp(pgno)%Comm,R_req(tt),ierr)
	enddo
	if(Nsendactive>0)then
		call MPI_waitall(Nsendactive,S_req,statuss,ierr)
	endif
	if(Nrecvactive>0)then
		call MPI_waitall(Nrecvactive,R_req,statusr,ierr)
	endif

	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo

	! pack the send buffer in the second pass
	do ii=1,block_i%ButterflyV%nblk_loc
	   ! convert indices from input to output
		if(mode=='R')then
			index_i = 1
			index_j = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
			index_i0 = 1
			index_j0 = index_j
			if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
		elseif(mode=='C')then
			index_i = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
			index_j = 1
			index_j0 = 1
			index_i0 = index_i
			if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
		endif

		call GetBlockPID(ptree,pgno_o,level_c_o,level_butterfly_c_o,index_i0,index_j0,mode,pid)

		pp=pid-ptree%pgrp(pgno)%head+1
		Nrow=size(block_i%ButterflyV%blocks(ii)%matrix,1)
		Ncol=size(block_i%ButterflyV%blocks(ii)%matrix,2)

		sendquant(pp)%dat(sendquant(pp)%size+1,1)=index_i
		sendquant(pp)%dat(sendquant(pp)%size+2,1)=index_j
		sendquant(pp)%dat(sendquant(pp)%size+3,1)=Nrow
		sendquant(pp)%dat(sendquant(pp)%size+4,1)=Ncol
		sendquant(pp)%size=sendquant(pp)%size+4
		do i=1,Nrow*Ncol
			rr = mod(i-1,Nrow)+1
			cc = (i-1)/Nrow+1
			sendquant(pp)%dat(sendquant(pp)%size+i,1) = block_i%ButterflyV%blocks(ii)%matrix(rr,cc)
		enddo
		sendquant(pp)%size=sendquant(pp)%size+Nrow*Ncol
		! deallocate(block_i%ButterflyV%blocks(ii)%matrix)
	enddo
	! if(allocated(block_i%ButterflyV%blocks))deallocate(block_i%ButterflyV%blocks)


	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif

	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(pgno)%Comm, ierr)

		do pp=1,nproc
			recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)

	else

		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			index_i=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			index_j=NINT(dble(recvquant(pp)%dat(i,1)))

			if(mode=='R')then
				nj=1
				sj=1
				if(index_j>num_blk/2)sj=2
				ni=2
				si=1
				ii=(index_j-block_o%sons(si,sj)%ButterflyV%idx)/block_o%sons(si,sj)%ButterflyV%inc+1
			elseif(mode=='C')then
				ni=1
				si=1
				if(index_i>num_blk/2)si=2
				nj=2
				sj=1
				ii=(index_i-block_o%sons(si,sj)%ButterflyV%idx)/block_o%sons(si,sj)%ButterflyV%inc+1
			endif
			i=i+1
			Nrow=NINT(dble(recvquant(pp)%dat(i,1)))
			i=i+1
			Ncol=NINT(dble(recvquant(pp)%dat(i,1)))
			do iii=si,si+ni-1
			do jjj=sj,sj+nj-1
				call assert(.not. allocated(block_o%sons(iii,jjj)%ButterflyV%blocks(ii)%matrix),'receiving dat alreay exists locally')
				allocate(block_o%sons(iii,jjj)%ButterflyV%blocks(ii)%matrix(Nrow,Ncol))
				do j=1,Nrow*Ncol
					rr = mod(j-1,Nrow)+1
					cc = (j-1)/Nrow+1
					block_o%sons(iii,jjj)%ButterflyV%blocks(ii)%matrix(rr,cc) = recvquant(pp)%dat(i+j,1)
				enddo
			enddo
			enddo
			i=i+Nrow*Ncol
		enddo
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)

	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1

end subroutine BF_all2all_V_split




subroutine BF_block_MVP_dat(blocks,chara,M,N,Nrnd,random1,random2,a,b,ptree,stats)

    use BPACK_DEFS
	use MISC_Utilities
    implicit none

    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,level_half,level_final
	integer idx_r,inc_r,nr,idx_c,inc_c,nc
	integer idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0
    DT ctemp, a, b
    character chara
	type(matrixblock)::blocks
	type(proctree)::ptree
	integer pgno,comm,ierr
	type(Hstat)::stats
	real(kind=8)::flop,flops
	integer index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc0,index_i_loc0, index_j_loc_s,index_j_loc_k

    type(butterfly_vec) :: BFvec
    DT :: random1(:,:), random2(:,:)
	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)

	level_butterfly=blocks%level_butterfly
	pgno = blocks%pgno
	! write(*,*)blocks%col_group,blocks%row_group,blocks%pgno,level_butterfly,'dd'
	comm = ptree%pgrp(pgno)%comm
	if(comm==MPI_COMM_NULL)then
		write(*,*)'ninin',pgno,comm==MPI_COMM_NULL,ptree%MyID
	endif

	call assert(IOwnPgrp(ptree,pgno),'I do not share this block!')

	if(level_butterfly==0)then
		rank = size(blocks%ButterflyU%blocks(1)%matrix,2)
		call assert(rank>0,'rank incorrect in blocks%ButterflyU')
		allocate(matrixtemp(rank,Nrnd))
		matrixtemp=0
		allocate(matrixtemp1(rank,Nrnd))
		matrixtemp1=0
		allocate(Vout_tmp(size(random2,1),size(random2,2)))
		Vout_tmp = 0
		! for implementation simplicity, MPI_ALLREDUCE is used even when nproc==1
		if (chara=='N') then !Vout=U*V^T*Vin
			call gemmf90(blocks%ButterflyV%blocks(1)%matrix,size(blocks%ButterflyV%blocks(1)%matrix,1),random1,size(random1,1),matrixtemp,rank,'T','N',rank,Nrnd,size(blocks%ButterflyV%blocks(1)%matrix,1),cone,czero,flop)
			stats%Flop_Tmp = stats%Flop_Tmp + flop
			call assert(MPI_COMM_NULL/=comm,'communicator should not be null 2')
			call MPI_ALLREDUCE(matrixtemp,matrixtemp1,rank*Nrnd,MPI_DT,MPI_SUM,comm,ierr)
			call gemmf90(blocks%ButterflyU%blocks(1)%matrix,size(blocks%ButterflyU%blocks(1)%matrix,1),matrixtemp1,rank,Vout_tmp,size(random2,1),'N','N',size(blocks%ButterflyU%blocks(1)%matrix,1),Nrnd,rank,cone,czero,flop)
			stats%Flop_Tmp = stats%Flop_Tmp + flop
			random2 = b*random2+a*Vout_tmp
		else if(chara=='T')then !Vout=V*U^T*Vin
			call gemmf90(blocks%ButterflyU%blocks(1)%matrix,size(blocks%ButterflyU%blocks(1)%matrix,1),random1,size(random1,1),matrixtemp,rank,'T','N',rank,Nrnd,size(blocks%ButterflyU%blocks(1)%matrix,1),cone,czero,flop)
			stats%Flop_Tmp = stats%Flop_Tmp + flop
			call assert(MPI_COMM_NULL/=comm,'communicator should not be null 3')
			call MPI_ALLREDUCE(matrixtemp,matrixtemp1,rank*Nrnd,MPI_DT,MPI_SUM,comm,ierr)
			call gemmf90(blocks%ButterflyV%blocks(1)%matrix,size(blocks%ButterflyV%blocks(1)%matrix,1),matrixtemp1,rank,Vout_tmp,size(random2,1),'N','N',size(blocks%ButterflyV%blocks(1)%matrix,1),Nrnd,rank,cone,czero,flop)
			stats%Flop_Tmp = stats%Flop_Tmp + flop
			random2 = b*random2+a*Vout_tmp
		endif

		deallocate(matrixtemp)
		deallocate(matrixtemp1)
		deallocate(Vout_tmp)

	else
		allocate(arr_acc_n(blocks%ButterflyV%nblk_loc))
		allocate(arr_acc_m(blocks%ButterflyU%nblk_loc))
		k1=0
		do i=1, blocks%ButterflyV%nblk_loc
			arr_acc_n(i) = k1
			nn=size(blocks%ButterflyV%blocks(i)%matrix,1)
			k1 =k1 +nn
		enddo

		k2=0
		do i=1, blocks%ButterflyU%nblk_loc
			arr_acc_m(i) = k2
			mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
			k2 =k2 +mm
		enddo

		num_vectors=Nrnd

		! if(BF_checkNAN(blocks))then
			! write(*,*)'NAN in 0 BF_block_MVP_dat'
			! stop
		! end if

		if (chara=='N') then

			if(isnan(sum(abs(random1(:,1))**2)))then
				write(*,*)'NAN in 1 BF_block_MVP_dat'
				stop
			end if

			level_butterfly=blocks%level_butterfly
			num_blocks=2**level_butterfly
			level_half = blocks%level_half


			allocate(BFvec%vec(0:level_butterfly+2))


			allocate (BFvec%vec(0)%blocks(1,blocks%ButterflyV%nblk_loc))
			BFvec%vec(0)%num_row=1
			BFvec%vec(0)%num_col=num_blocks
			BFvec%vec(0)%idx_r=1
			BFvec%vec(0)%inc_r=1
			BFvec%vec(0)%nr=1
			BFvec%vec(0)%idx_c=blocks%ButterflyV%idx
			BFvec%vec(0)%inc_c=blocks%ButterflyV%inc
			BFvec%vec(0)%nc=blocks%ButterflyV%nblk_loc

			!$omp parallel do default(shared) private(i,nn,ii,jj)
			do i=1, BFvec%vec(0)%nc
				nn=size(blocks%ButterflyV%blocks(i)%matrix,1)
				allocate (BFvec%vec(0)%blocks(1,i)%matrix(nn,num_vectors))
				do ii=1, nn
					do jj=1, num_vectors
						BFvec%vec(0)%blocks(1,i)%matrix(ii,jj)=random1(ii+arr_acc_n(i),jj)
					enddo
				enddo
			enddo
			!$omp end parallel do

			do level=0, level_half
				call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')

				BFvec%vec(level+1)%idx_r=idx_r
				BFvec%vec(level+1)%inc_r=inc_r
				BFvec%vec(level+1)%nr=nr
				BFvec%vec(level+1)%idx_c=idx_c
				BFvec%vec(level+1)%inc_c=inc_c
				BFvec%vec(level+1)%nc=nc
				if (level/=level_butterfly+1) then
					BFvec%vec(level+1)%num_row=2**level
					BFvec%vec(level+1)%num_col=2**(level_butterfly-level)
				else
					BFvec%vec(level+1)%num_row=2**level_butterfly
					BFvec%vec(level+1)%num_col=1
				endif
				if(level_half/=level)then ! the last level doesn't require doubling block columns
				if(nc==1 .and. 2**(level_butterfly-level)>1)then ! double the number of local block columns used for MPI communication
					BFvec%vec(level+1)%nc=2
					BFvec%vec(level+1)%idx_c=BFvec%vec(level+1)%idx_c-1+mod(BFvec%vec(level+1)%idx_c,2)
				endif
				endif
				allocate(BFvec%vec(level+1)%blocks(BFvec%vec(level+1)%nr,BFvec%vec(level+1)%nc))

				if (level==0) then
					flops=0
					!$omp parallel do default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s) reduction(+:flops)
					do j=1, blocks%ButterflyV%nblk_loc
						index_j=(j-1)*inc_c+idx_c
						index_j_loc_s=(index_j-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
						rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
						nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
						allocate (BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix(rank,num_vectors))
						BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix=0

						call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(0)%blocks(1,j)%matrix,nn,BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix,rank,'T','N',rank,num_vectors,nn,cone,czero,flop=flop)
						flops=flops+flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				elseif (level==level_butterfly+1) then
					flops=0
					!$omp parallel do default(shared) private(i,rank,mm,flop) reduction(+:flops)
					do i=1, blocks%ButterflyU%nblk_loc
						rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
						mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
						allocate (BFvec%vec(level+1)%blocks(i,1)%matrix(mm,num_vectors))
						BFvec%vec(level+1)%blocks(i,1)%matrix=0

						call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(level)%blocks(i,1)%matrix,rank,BFvec%vec(level+1)%blocks(i,1)%matrix,mm,'N','N',mm,num_vectors,rank,cone,czero,flop=flop)
						flops = flops + flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				else
					flops=0
					!$omp parallel do default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)

					do index_ij=1, nr*nc
						index_j_loc = (index_ij-1)/nr+1
						index_i_loc= mod(index_ij-1,nr) + 1  !index_i_loc is local index of row-wise ordering at current level
						index_i=(index_i_loc-1)*inc_r+idx_r  !index_i is global index of row-wise ordering at current level
						index_j=(index_j_loc-1)*inc_c+idx_c

						index_ii=int((index_i+1)/2) ; index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level)

						index_ii_loc=(index_ii-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level)
						index_jj_loc=(index_jj-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1


						index_i_loc_s=(index_i-BFvec%vec(level+1)%idx_r)/BFvec%vec(level+1)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level+1)
						index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
						index_j_loc_s=(index_j-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
						index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

						nn1=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
						nn2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,2)
						mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)

						allocate (BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm,num_vectors))
						BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0

						call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn1,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn1,cone,cone,flop=flop)
						flops = flops + flop

						call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix,nn2,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn2,cone,cone,flop=flop)
						flops = flops + flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				endif

				do j=1, BFvec%vec(level)%nc
					do i=1, BFvec%vec(level)%nr
						if(allocated(BFvec%vec(level)%blocks(i,j)%matrix))deallocate (BFvec%vec(level)%blocks(i,j)%matrix)
					enddo
				enddo
				if(level_half/=level)then
					call BF_exchange_matvec(blocks,BFvec%vec(level+1),num_vectors,stats,ptree,level,'R','B')
				endif
			enddo

			if(level_half+1/=0)then
				call BF_all2all_matvec(blocks,BFvec%vec(level_half+1),num_vectors,stats,ptree,level_half,'R','C')
			else
				call BF_all2all_matvec(blocks,BFvec%vec(level_half+1),num_vectors,stats,ptree,level_half+1,'R','C')
			endif

			do level=level_half+1,level_butterfly+1
				call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0,'C')

				! convert the local column-wise kernel block ranges to local row-wise output vector ranges
				if(level/=0 .and. level/=level_butterfly+1)then
					idx_r = idx_r0*2-1
					nr = nr0*2
					inc_r=inc_r0
					idx_c = ceiling_safe(idx_c0/2d0)
					if(inc_c0>1)then
					nc=nc0
					else
					nc = ceiling_safe(nc0/2d0)
					endif
					inc_c = ceiling_safe(inc_c0/2d0)
				else
					idx_r=idx_r0
					nr=nr0
					inc_r=inc_r0
					idx_c=idx_c0
					nc=nc0
					inc_c=inc_c0
				endif

				BFvec%vec(level+1)%idx_r=idx_r
				BFvec%vec(level+1)%inc_r=inc_r
				BFvec%vec(level+1)%nr=nr
				BFvec%vec(level+1)%idx_c=idx_c
				BFvec%vec(level+1)%inc_c=inc_c
				BFvec%vec(level+1)%nc=nc
				if (level/=level_butterfly+1) then
					BFvec%vec(level+1)%num_row=2**level
					BFvec%vec(level+1)%num_col=2**(level_butterfly-level)
				else
					BFvec%vec(level+1)%num_row=2**level_butterfly
					BFvec%vec(level+1)%num_col=1
				endif

				allocate(BFvec%vec(level+1)%blocks(BFvec%vec(level+1)%nr,BFvec%vec(level+1)%nc))


				if (level==0) then
					flops=0
					!$omp parallel do default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s) reduction(+:flops)
					do j=1, nc0
						index_j=(j-1)*inc_c0+idx_c0
						index_j_loc_s=(index_j-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
						rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
						nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
						allocate (BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix(rank,num_vectors))
						BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix=0
						call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(0)%blocks(1,j)%matrix,nn,BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix,rank,'T','N',rank,num_vectors,nn,cone,czero,flop=flop)
						flops=flops+flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				elseif (level==level_butterfly+1) then
					flops=0
					!$omp parallel do default(shared) private(i,index_i,index_i_loc_s,rank,mm,flop) reduction(+:flops)
					do i=1, nr0
						index_i=(i-1)*inc_r0+idx_r0
						index_i_loc_s=(index_i-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1

						rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
						mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
						allocate (BFvec%vec(level+1)%blocks(i,1)%matrix(mm,num_vectors))
						BFvec%vec(level+1)%blocks(i,1)%matrix=0

						call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(level)%blocks(index_i_loc_s,1)%matrix,rank,BFvec%vec(level+1)%blocks(i,1)%matrix,mm,'N','N',mm,num_vectors,rank,cone,czero,flop=flop)
						flops = flops + flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				else
					flops=0

					if(nc0>1 .and. inc_c0==1)then  ! this special treatment makes sure two threads do not write to the same address simultaneously
						!$omp parallel do default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,index_j_loc0,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
						do index_ij=1, nr0*nc0/2
							index_j_loc0=(index_ij-1)/nr0+1
							do jj=1,2
							index_j_loc = 2*(index_j_loc0-1)+jj       !index_i_loc is local index of column-wise ordering at current level
							index_i_loc= mod(index_ij-1,nr0) + 1
							index_i=(index_i_loc-1)*inc_r0+idx_r0  !index_i is global index of column-wise ordering at current level
							index_j=(index_j_loc-1)*inc_c0+idx_c0

							index_ii_loc=(index_i-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1  !index_ii_loc is local index in BFvec%vec(level)
							index_jj_loc=(index_j-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1

							index_ii=2*index_i-1; index_jj=int((index_j+1)/2)  !index_ii is global index in BFvec%vec(level+1)

							index_i_loc_s=(index_ii-BFvec%vec(level+1)%idx_r)/BFvec%vec(level+1)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level+1)
							index_i_loc_k=(2*index_i-1-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
							index_j_loc_s=(index_jj-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
							index_j_loc_k=(index_j-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm,num_vectors))
								BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn,cone,cone,flop=flop)
							flops = flops + flop

							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix(mm,num_vectors))
								BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn,BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn,cone,cone,flop=flop)
							flops = flops + flop
							enddo
						enddo
						!$omp end parallel do

					else
						!$omp parallel do default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
						do index_ij=1, nr0*nc0
							index_j_loc = (index_ij-1)/nr0+1       !index_i_loc is local index of column-wise ordering at current level
							index_i_loc= mod(index_ij-1,nr0) + 1
							index_i=(index_i_loc-1)*inc_r0+idx_r0  !index_i is global index of column-wise ordering at current level
							index_j=(index_j_loc-1)*inc_c0+idx_c0

							index_ii_loc=(index_i-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1  !index_ii_loc is local index in BFvec%vec(level)
							index_jj_loc=(index_j-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1

							index_ii=2*index_i-1; index_jj=int((index_j+1)/2)  !index_ii is global index in BFvec%vec(level+1)

							index_i_loc_s=(index_ii-BFvec%vec(level+1)%idx_r)/BFvec%vec(level+1)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level+1)
							index_i_loc_k=(2*index_i-1-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
							index_j_loc_s=(index_jj-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
							index_j_loc_k=(index_j-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm,num_vectors))
								BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn,cone,cone,flop=flop)
							flops = flops + flop

							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix(mm,num_vectors))
								BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn,BFvec%vec(level+1)%blocks(index_i_loc_s+1,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn,cone,cone,flop=flop)
							flops = flops + flop

						enddo
						!$omp end parallel do
					endif
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				endif


				do j=1, BFvec%vec(level)%nc
					do i=1, BFvec%vec(level)%nr
						if(allocated(BFvec%vec(level)%blocks(i,j)%matrix))deallocate (BFvec%vec(level)%blocks(i,j)%matrix)
					enddo
				enddo

				if(level/=level_butterfly+1)then
					call BF_exchange_matvec(blocks,BFvec%vec(level+1),num_vectors,stats,ptree,level,'R','R')
				endif
			enddo

			!$omp parallel do default(shared) private(i,mm,ii,jj)
			do jj=1, num_vectors
				do i=1, blocks%ButterflyU%nblk_loc
					mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
					do ii=1, mm
						random2(ii+arr_acc_m(i),jj)=b*random2(ii+arr_acc_m(i),jj)+a*BFvec%vec(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)
					enddo
				enddo
			enddo
			!$omp end parallel do

			if(isnan(sum(abs(random2(:,1))**2)))then
				write(*,*)'NAN in 2 BF_block_MVP_dat',blocks%row_group,blocks%col_group,blocks%level,blocks%level_butterfly
				stop
			end if
			!deallocate (BFvec%vec)

		elseif (chara=='T') then

			level_butterfly=blocks%level_butterfly
			num_blocks=2**level_butterfly
			level_half = blocks%level_half

			allocate (BFvec%vec(0:level_butterfly+2))
			allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc,1))
			BFvec%vec(0)%num_row=num_blocks
			BFvec%vec(0)%num_col=1
			BFvec%vec(0)%idx_r=blocks%ButterflyU%idx
			BFvec%vec(0)%inc_r=blocks%ButterflyU%inc
			BFvec%vec(0)%nr=blocks%ButterflyU%nblk_loc
			BFvec%vec(0)%idx_c=1
			BFvec%vec(0)%inc_c=1
			BFvec%vec(0)%nc=1

			!$omp parallel do default(shared) private(i,mm,ii,jj)
			do i=1, BFvec%vec(0)%nr
				mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
				allocate (BFvec%vec(0)%blocks(i,1)%matrix(mm,num_vectors))
				do ii=1, mm
					do jj=1, num_vectors
						BFvec%vec(0)%blocks(i,1)%matrix(ii,jj)=random1(ii+arr_acc_m(i),jj)
					enddo
				enddo
			enddo
			!$omp end parallel do


			do level=level_butterfly+1, level_half+1,-1
				call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')

				BFvec%vec(level_butterfly-level+2)%idx_r=idx_r
				BFvec%vec(level_butterfly-level+2)%inc_r=inc_r
				BFvec%vec(level_butterfly-level+2)%nr=nr
				BFvec%vec(level_butterfly-level+2)%idx_c=idx_c
				BFvec%vec(level_butterfly-level+2)%inc_c=inc_c
				BFvec%vec(level_butterfly-level+2)%nc=nc
				if (level/=0) then
					BFvec%vec(level_butterfly-level+2)%num_row=2**(level-1)
					BFvec%vec(level_butterfly-level+2)%num_col=2**(level_butterfly-level+1)
				else
					BFvec%vec(level_butterfly-level+2)%num_row=1
					BFvec%vec(level_butterfly-level+2)%num_col=2**level_butterfly
				endif
				if(level_half+1/=level)then ! the last level doesn't require doubling block rows
				if(nr==1 .and. 2**(level-1)>1)then ! double the number of local block rows used for MPI communication
					BFvec%vec(level_butterfly-level+2)%nr=2
					BFvec%vec(level_butterfly-level+2)%idx_r=BFvec%vec(level_butterfly-level+2)%idx_r-1+mod(BFvec%vec(level_butterfly-level+2)%idx_r,2)
				endif
				endif
				allocate(BFvec%vec(level_butterfly-level+2)%blocks(BFvec%vec(level_butterfly-level+2)%nr,BFvec%vec(level_butterfly-level+2)%nc))

				if (level==level_butterfly+1) then
					flops=0
					!$omp parallel do default(shared) private(i,rank,mm,flop,index_i,index_i_loc_s) reduction(+:flops)
					do i=1, blocks%ButterflyU%nblk_loc
						index_i=(i-1)*blocks%ButterflyU%inc+blocks%ButterflyU%idx
						index_i_loc_s=(index_i-BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r+1

						rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
						mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
						allocate (BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix(rank,num_vectors))
						BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix=0

						call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(0)%blocks(i,1)%matrix,mm,BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix,rank,'T','N',rank,num_vectors,mm,cone,czero,flop=flop)
						flops = flops + flop

					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				elseif (level==0) then
					flops=0
					!$omp parallel do default(shared) private(j,rank,nn,flop) reduction(+:flops)
					do j=1, blocks%ButterflyV%nblk_loc
						nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
						rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
						allocate (BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix(nn,num_vectors))
						BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix=0
						call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(level_butterfly+1)%blocks(1,j)%matrix,rank,BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix,nn,'N','N',nn,num_vectors,rank,cone,czero,flop=flop)
						flops = flops + flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				else
					flops=0
					!$omp parallel do default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
					do index_ij=1, nr*nc
						index_j_loc = (index_ij-1)/nr+1
						index_i_loc= mod(index_ij-1,nr) + 1  !index_i_loc is local index of column-wise ordering at current level
						index_i=(index_i_loc-1)*inc_r+idx_r  !index_i is global index of column-wise ordering at current level
						index_j=(index_j_loc-1)*inc_c+idx_c

						index_ii=2*index_i-1; index_jj=int((index_j+1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

						index_ii_loc=(index_ii-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
						index_jj_loc=(index_jj-BFvec%vec(level_butterfly-level+1)%idx_c)/BFvec%vec(level_butterfly-level+1)%inc_c+1

						index_i_loc_s=(index_i-BFvec%vec(level_butterfly-level+2)%idx_r)/BFvec%vec(level_butterfly-level+2)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
						index_i_loc_k=(2*index_i-1-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
						index_j_loc_s=(index_j-BFvec%vec(level_butterfly-level+2)%idx_c)/BFvec%vec(level_butterfly-level+2)%inc_c+1
						index_j_loc_k=(index_j-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

						mm1=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
						mm2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,1)
						nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
						allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix(nn,num_vectors))
						BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0


						call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm1,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm1,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm1,cone,cone,flop=flop)
						flops = flops + flop
						call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,mm2,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc+1,index_jj_loc)%matrix,mm2,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm2,cone,cone,flop=flop)
						flops = flops + flop

					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				endif

				do j=1, BFvec%vec(level_butterfly-level+1)%nc
					do i=1, BFvec%vec(level_butterfly-level+1)%nr
						deallocate (BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix)
					enddo
				enddo

				if(level_half+1/=level)then
					call BF_exchange_matvec(blocks,BFvec%vec(level_butterfly-level+2),num_vectors,stats,ptree,level,'C','B')
				endif
			enddo

			if(level_half/=level_butterfly+1)then
				call BF_all2all_matvec(blocks,BFvec%vec(level_butterfly-level_half+1),num_vectors,stats,ptree,level_half+1,'C','R')
			else
				call BF_all2all_matvec(blocks,BFvec%vec(level_butterfly-level_half+1),num_vectors,stats,ptree,level_half,'C','R')
			endif

			do level=level_half,0,-1
				call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0,'R')

				! convert the local row-wise kernel block ranges to local column-wise output vector ranges
				if(level/=0 .and. level/=level_butterfly+1)then
					idx_r = ceiling_safe(idx_r0/2d0)
					if(inc_r0>1)then
					nr=nr0
					else
					nr = ceiling_safe(nr0/2d0)
					endif
					inc_r = ceiling_safe(inc_r0/2d0)
					idx_c = idx_c0*2-1
					nc = nc0*2
					inc_c=inc_c0
				else
					idx_r=idx_r0
					nr=nr0
					inc_r=inc_r0
					idx_c=idx_c0
					nc=nc0
					inc_c=inc_c0
				endif

				BFvec%vec(level_butterfly-level+2)%idx_r=idx_r
				BFvec%vec(level_butterfly-level+2)%inc_r=inc_r
				BFvec%vec(level_butterfly-level+2)%nr=nr
				BFvec%vec(level_butterfly-level+2)%idx_c=idx_c
				BFvec%vec(level_butterfly-level+2)%inc_c=inc_c
				BFvec%vec(level_butterfly-level+2)%nc=nc
				if (level/=0) then
					BFvec%vec(level+1)%num_row=2**(level-1)
					BFvec%vec(level+1)%num_col=2**(level_butterfly-level+1)
				else
					BFvec%vec(level+1)%num_row=1
					BFvec%vec(level+1)%num_col=2**level_butterfly
				endif

				allocate(BFvec%vec(level_butterfly-level+2)%blocks(BFvec%vec(level_butterfly-level+2)%nr,BFvec%vec(level_butterfly-level+2)%nc))


				if (level==level_butterfly+1) then
					flops=0
					!$omp parallel do default(shared) private(i,rank,mm,flop,index_i,index_i_loc_s) reduction(+:flops)
					do i=1, blocks%ButterflyU%nblk_loc
						index_i=(i-1)*blocks%ButterflyU%inc+blocks%ButterflyU%idx
						index_i_loc_s=(index_i-BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r+1

						rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
						mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
						allocate (BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix(rank,num_vectors))
						BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix=0

						call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(0)%blocks(i,1)%matrix,mm,BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix,rank,'T','N',rank,num_vectors,mm,cone,czero,flop=flop)
						flops = flops + flop

					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				elseif (level==0) then
					flops=0
					!$omp parallel do default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s) reduction(+:flops)
					do j=1, blocks%ButterflyV%nblk_loc
						index_j=(j-1)*blocks%ButterflyV%inc+blocks%ButterflyV%idx
						index_j_loc_s=(index_j-BFvec%vec(level_butterfly+1)%idx_c)/BFvec%vec(level_butterfly+1)%inc_c+1

						nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
						rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
						allocate (BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix(nn,num_vectors))
						BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix=0
						call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(level_butterfly+1)%blocks(1,index_j_loc_s)%matrix,rank,BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix,nn,'N','N',nn,num_vectors,rank,cone,czero,flop=flop)
						flops = flops + flop
					enddo
					!$omp end parallel do
					stats%Flop_Tmp = stats%Flop_Tmp + flops
				else

					flops=0
					if(nr0>1 .and. inc_r0==1)then ! this special treatment makes sure two threads do not write to the same address simultaneously
	!$omp parallel do default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,index_i_loc0,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
						do index_ij=1, nr0*nc0/2
						index_i_loc0 = (index_ij-1)/nc0+1
						do ii=1,2
							index_i_loc=(index_i_loc0-1)*2+ii
							index_j_loc= mod(index_ij-1,nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
							index_i=(index_i_loc-1)*inc_r0+idx_r0  !index_i is global index of row-wise ordering at current level
							index_j=(index_j_loc-1)*inc_c0+idx_c0

							index_ii=int((index_i+1)/2); index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

							index_ii_loc=(index_i-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
							index_jj_loc=(index_j-BFvec%vec(level_butterfly-level+1)%idx_c)/BFvec%vec(level_butterfly-level+1)%inc_c+1

							index_i_loc_s=(index_ii-BFvec%vec(level_butterfly-level+2)%idx_r)/BFvec%vec(level_butterfly-level+2)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
							index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
							index_j_loc_s=(index_jj-BFvec%vec(level_butterfly-level+2)%idx_c)/BFvec%vec(level_butterfly-level+2)%inc_c+1
							index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1


							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix(nn,num_vectors))
								BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							! write(*,*)index_ii_loc,index_jj_loc,shape(BFvec%vec(level_butterfly-level+1)%blocks),index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level_butterfly-level+2)%blocks),'lv:',level,shape(blocks%ButterflyKerl(level)%blocks)
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm,cone,cone,flop=flop)
							flops = flops + flop


							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix))then
								allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix(nn,num_vectors))
								BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,mm,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix,nn,'T','N',nn,num_vectors,mm,cone,cone,flop=flop)
							flops = flops + flop
						enddo
						enddo
						!$omp end parallel do
					else
	!$omp parallel do default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
						do index_ij=1, nr0*nc0
							index_j_loc = (index_ij-1)/nr0+1
							index_i_loc= mod(index_ij-1,nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
							index_i=(index_i_loc-1)*inc_r0+idx_r0  !index_i is global index of row-wise ordering at current level
							index_j=(index_j_loc-1)*inc_c0+idx_c0

							index_ii=int((index_i+1)/2); index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

							index_ii_loc=(index_i-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
							index_jj_loc=(index_j-BFvec%vec(level_butterfly-level+1)%idx_c)/BFvec%vec(level_butterfly-level+1)%inc_c+1

							index_i_loc_s=(index_ii-BFvec%vec(level_butterfly-level+2)%idx_r)/BFvec%vec(level_butterfly-level+2)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
							index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
							index_j_loc_s=(index_jj-BFvec%vec(level_butterfly-level+2)%idx_c)/BFvec%vec(level_butterfly-level+2)%inc_c+1
							index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1


							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then
								allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix(nn,num_vectors))
								BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0
							endif
							! !$omp end critical
							! write(*,*)index_ii_loc,index_jj_loc,shape(BFvec%vec(level_butterfly-level+1)%blocks),index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level_butterfly-level+2)%blocks),'lv:',level,shape(blocks%ButterflyKerl(level)%blocks)
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm,cone,cone,flop=flop)
							flops = flops + flop


							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,1)
							nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,2)
							! !$omp critical
							if(.not. allocated(BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix))then
								allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix(nn,num_vectors))
								BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix=0
							endif
							! !$omp end critical
							call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,mm,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s+1)%matrix,nn,'T','N',nn,num_vectors,mm,cone,cone,flop=flop)
							flops = flops + flop
						enddo
						!$omp end parallel do

					endif

					stats%Flop_Tmp = stats%Flop_Tmp + flops
				endif


				do j=1, BFvec%vec(level_butterfly-level+1)%nc
					do i=1, BFvec%vec(level_butterfly-level+1)%nr
						deallocate (BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix)
					enddo
				enddo

				if(level/=0)then
					call BF_exchange_matvec(blocks,BFvec%vec(level_butterfly-level+2),num_vectors,stats,ptree,level,'C','R')
				endif
			enddo

			!$omp parallel do default(shared) private(j,nn,ii,jj)
			do jj=1, num_vectors
				do j=1, blocks%ButterflyV%nblk_loc
					nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
					! write(*,*)nn,arr_acc_n(j)
					do ii=1, nn
						random2(ii+arr_acc_n(j),jj)=b*random2(ii+arr_acc_n(j),jj)+a*BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix(ii,jj)
					enddo
				enddo
			enddo
			!$omp end parallel do
		endif


		do level=0, level_butterfly+2
			do j=1, BFvec%vec(level)%nc
				do i=1, BFvec%vec(level)%nr
					if(allocated(BFvec%vec(level)%blocks(i,j)%matrix))deallocate (BFvec%vec(level)%blocks(i,j)%matrix)
				enddo
			enddo
			deallocate (BFvec%vec(level)%blocks)
		enddo

		deallocate (BFvec%vec)
		deallocate(arr_acc_m,arr_acc_n)


	endif

    return

end subroutine BF_block_MVP_dat



!**** Matvec of partial levels of BF with vectors
! if chara=='N', out=BF(level_end:0)*vec, if chara=='T', out=vec*BF(level_butterfly+1:level_end).
	!blocks: working BF
	!chara: 'N' or 'T'
	!num_vectors: number of vectors
	!VectIn: dimension (mnloc,num_vectors) the local input vectors
	!BFvec: storing the result of partial matvec
	!level_end: the last level of the operator
	!ptree: process tree
	!stats: statistics
subroutine BF_block_MVP_partial(blocks,chara,num_vectors,VectIn,BFvec,level_end,ptree,msh,stats)

    use BPACK_DEFS
	use MISC_Utilities
    implicit none

    integer index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,level_half,level_final
	integer idx_r,inc_r,nr,idx_c,inc_c,nc
	integer idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0
    DT ctemp
    character chara
	type(matrixblock)::blocks
	type(proctree)::ptree
	integer pgno,comm,ierr
	type(Hstat)::stats
	type(mesh)::msh
	real(kind=8)::flop,flops
	integer index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,level_end

    type(butterfly_vec) :: BFvec
    DT :: VectIn(:,:)
	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)

	integer idxs,groupn_start,groupm_start

	level_butterfly=blocks%level_butterfly
	pgno = blocks%pgno
	comm = ptree%pgrp(pgno)%comm
	if(comm==MPI_COMM_NULL)then
		write(*,*)'ninin',pgno,comm==MPI_COMM_NULL,ptree%MyID
	endif

	call assert(IOwnPgrp(ptree,pgno),'I do not share this block!')



	if(BF_checkNAN(blocks))then
		write(*,*)'NAN in 0 BF_block_MVP_partial'
		stop
	end if

	if (chara=='N') then

		if(isnan(sum(abs(VectIn(:,1))**2)))then
			write(*,*)'NAN in 1 BF_block_MVP_partial'
			stop
		end if

		level_butterfly=blocks%level_butterfly
		num_blocks=2**level_butterfly
		level_half = blocks%level_half
		call assert(level_half>=level_end,'partial matvec with chara=N requires row-wise ordering')

		if(.not. allocated(BFvec%vec))allocate(BFvec%vec(0:level_butterfly+2))

		allocate (BFvec%vec(0)%blocks(1,blocks%ButterflyV%nblk_loc))
		BFvec%vec(0)%num_row=1
		BFvec%vec(0)%num_col=num_blocks
		BFvec%vec(0)%idx_r=1
		BFvec%vec(0)%inc_r=1
		BFvec%vec(0)%nr=1
		BFvec%vec(0)%idx_c=blocks%ButterflyV%idx
		BFvec%vec(0)%inc_c=blocks%ButterflyV%inc
		BFvec%vec(0)%nc=blocks%ButterflyV%nblk_loc



		groupn_start=blocks%col_group*2**level_butterfly
		call GetLocalBlockRange(ptree,blocks%pgno,0,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')
		!$omp parallel do default(shared) private(i,nn,ii,jj,idxs)
		do i=1, BFvec%vec(0)%nc
			nn = msh%basis_group(groupn_start+(i-1)*inc_c+idx_c-1)%tail-msh%basis_group(groupn_start+(i-1)*inc_c+idx_c-1)%head+1
			idxs = msh%basis_group(groupn_start+(i-1)*inc_c+idx_c-1)%head-msh%basis_group(groupn_start+idx_c-1)%head
			allocate (BFvec%vec(0)%blocks(1,i)%matrix(nn,num_vectors))
			do ii=1, nn
				do jj=1, num_vectors
					BFvec%vec(0)%blocks(1,i)%matrix(ii,jj)=VectIn(ii+idxs,jj)
				enddo
			enddo
		enddo
		!$omp end parallel do




		do level=0, level_end
			call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')

			BFvec%vec(level+1)%idx_r=idx_r
			BFvec%vec(level+1)%inc_r=inc_r
			BFvec%vec(level+1)%nr=nr
			BFvec%vec(level+1)%idx_c=idx_c
			BFvec%vec(level+1)%inc_c=inc_c
			BFvec%vec(level+1)%nc=nc
			if (level/=level_butterfly+1) then
				BFvec%vec(level+1)%num_row=2**level
				BFvec%vec(level+1)%num_col=2**(level_butterfly-level)
			else
				BFvec%vec(level+1)%num_row=2**level_butterfly
				BFvec%vec(level+1)%num_col=1
			endif
			if(level_half/=level)then ! the last level doesn't require doubling block columns
			if(nc==1 .and. 2**(level_butterfly-level)>1)then ! double the number of local block columns used for MPI communication
				BFvec%vec(level+1)%nc=2
				BFvec%vec(level+1)%idx_c=BFvec%vec(level+1)%idx_c-1+mod(BFvec%vec(level+1)%idx_c,2)
			endif
			endif
			allocate(BFvec%vec(level+1)%blocks(BFvec%vec(level+1)%nr,BFvec%vec(level+1)%nc))

			if (level==0) then
				flops=0
				!$omp parallel do default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s) reduction(+:flops)
				do j=1, blocks%ButterflyV%nblk_loc
					index_j=(j-1)*inc_c+idx_c
					index_j_loc_s=(index_j-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
					rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
					nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
					allocate (BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix(rank,num_vectors))
					BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix=0

					call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(0)%blocks(1,j)%matrix,nn,BFvec%vec(1)%blocks(1,index_j_loc_s)%matrix,rank,'T','N',rank,num_vectors,nn,cone,czero,flop=flop)
					flops=flops+flop
				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			elseif (level==level_butterfly+1) then
				flops=0
				!$omp parallel do default(shared) private(i,rank,mm,flop) reduction(+:flops)
				do i=1, blocks%ButterflyU%nblk_loc
					rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
					mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
					allocate (BFvec%vec(level+1)%blocks(i,1)%matrix(mm,num_vectors))
					BFvec%vec(level+1)%blocks(i,1)%matrix=0

					call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(level)%blocks(i,1)%matrix,rank,BFvec%vec(level+1)%blocks(i,1)%matrix,mm,'N','N',mm,num_vectors,rank,cone,czero,flop=flop)
					flops = flops + flop
				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			else
				flops=0
				!$omp parallel do default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)

				do index_ij=1, nr*nc
					index_j_loc = (index_ij-1)/nr+1
					index_i_loc= mod(index_ij-1,nr) + 1  !index_i_loc is local index of row-wise ordering at current level
					index_i=(index_i_loc-1)*inc_r+idx_r  !index_i is global index of row-wise ordering at current level
					index_j=(index_j_loc-1)*inc_c+idx_c

					index_ii=int((index_i+1)/2) ; index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level)

					index_ii_loc=(index_ii-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level)
					index_jj_loc=(index_jj-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1


					index_i_loc_s=(index_i-BFvec%vec(level+1)%idx_r)/BFvec%vec(level+1)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level+1)
					index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
					index_j_loc_s=(index_j-BFvec%vec(level+1)%idx_c)/BFvec%vec(level+1)%inc_c+1
					index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

					nn1=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
					nn2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,2)
					mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)

					allocate (BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm,num_vectors))
					BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0

					call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,nn1,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn1,cone,cone,flop=flop)
					flops = flops + flop

					call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix,nn2,BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,mm,'N','N',mm,num_vectors,nn2,cone,cone,flop=flop)
					flops = flops + flop
				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			endif

			do j=1, BFvec%vec(level)%nc
				do i=1, BFvec%vec(level)%nr
					if(allocated(BFvec%vec(level)%blocks(i,j)%matrix))deallocate (BFvec%vec(level)%blocks(i,j)%matrix)
				enddo
			enddo
			if(level_half/=level)then
				call BF_exchange_matvec(blocks,BFvec%vec(level+1),num_vectors,stats,ptree,level,'R','B')
			endif
		enddo
	elseif (chara=='T') then
		level_butterfly=blocks%level_butterfly
		num_blocks=2**level_butterfly
		level_half = blocks%level_half
		call assert(level_half+1<=level_end,'partial matvec with chara=T requires column-wise ordering')

		if(.not. allocated(BFvec%vec))allocate(BFvec%vec(0:level_butterfly+2))
		allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc,1))
		BFvec%vec(0)%num_row=num_blocks
		BFvec%vec(0)%num_col=1
		BFvec%vec(0)%idx_r=blocks%ButterflyU%idx
		BFvec%vec(0)%inc_r=blocks%ButterflyU%inc
		BFvec%vec(0)%nr=blocks%ButterflyU%nblk_loc
		BFvec%vec(0)%idx_c=1
		BFvec%vec(0)%inc_c=1
		BFvec%vec(0)%nc=1


		groupm_start=blocks%row_group*2**level_butterfly
		call GetLocalBlockRange(ptree,blocks%pgno,level_butterfly+1,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')
		!$omp parallel do default(shared) private(i,mm,ii,jj,idxs)
		do i=1, BFvec%vec(0)%nr
			mm = msh%basis_group(groupm_start+(i-1)*inc_r+idx_r-1)%tail-msh%basis_group(groupm_start+(i-1)*inc_r+idx_r-1)%head+1
			idxs = msh%basis_group(groupm_start+(i-1)*inc_r+idx_r-1)%head-msh%basis_group(groupm_start+idx_r-1)%head
			allocate (BFvec%vec(0)%blocks(i,1)%matrix(mm,num_vectors))
			do ii=1, mm
				do jj=1, num_vectors
					BFvec%vec(0)%blocks(i,1)%matrix(ii,jj)=VectIn(ii+idxs,jj)
				enddo
			enddo
		enddo
		!$omp end parallel do


		do level=level_butterfly+1, level_end,-1
			call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')

			BFvec%vec(level_butterfly-level+2)%idx_r=idx_r
			BFvec%vec(level_butterfly-level+2)%inc_r=inc_r
			BFvec%vec(level_butterfly-level+2)%nr=nr
			BFvec%vec(level_butterfly-level+2)%idx_c=idx_c
			BFvec%vec(level_butterfly-level+2)%inc_c=inc_c
			BFvec%vec(level_butterfly-level+2)%nc=nc
			if (level/=0) then
				BFvec%vec(level_butterfly-level+2)%num_row=2**(level-1)
				BFvec%vec(level_butterfly-level+2)%num_col=2**(level_butterfly-level+1)
			else
				BFvec%vec(level_butterfly-level+2)%num_row=1
				BFvec%vec(level_butterfly-level+2)%num_col=2**level_butterfly
			endif
			if(level_half+1/=level)then ! the last level doesn't require doubling block rows
			if(nr==1 .and. 2**(level-1)>1)then ! double the number of local block rows used for MPI communication
				BFvec%vec(level_butterfly-level+2)%nr=2
				BFvec%vec(level_butterfly-level+2)%idx_r=BFvec%vec(level_butterfly-level+2)%idx_r-1+mod(BFvec%vec(level_butterfly-level+2)%idx_r,2)
			endif
			endif
			allocate(BFvec%vec(level_butterfly-level+2)%blocks(BFvec%vec(level_butterfly-level+2)%nr,BFvec%vec(level_butterfly-level+2)%nc))

			if (level==level_butterfly+1) then
				flops=0
				!$omp parallel do default(shared) private(i,rank,mm,flop,index_i,index_i_loc_s) reduction(+:flops)
				do i=1, blocks%ButterflyU%nblk_loc
					index_i=(i-1)*blocks%ButterflyU%inc+blocks%ButterflyU%idx
					index_i_loc_s=(index_i-BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r+1

					rank=size(blocks%ButterflyU%blocks(i)%matrix,2)
					mm=size(blocks%ButterflyU%blocks(i)%matrix,1)
					allocate (BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix(rank,num_vectors))
					BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix=0

					call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm,BFvec%vec(0)%blocks(i,1)%matrix,mm,BFvec%vec(1)%blocks(index_i_loc_s,1)%matrix,rank,'T','N',rank,num_vectors,mm,cone,czero,flop=flop)
					flops = flops + flop

				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			elseif (level==0) then
				flops=0
				!$omp parallel do default(shared) private(j,rank,nn,flop) reduction(+:flops)
				do j=1, blocks%ButterflyV%nblk_loc
					nn=size(blocks%ButterflyV%blocks(j)%matrix,1)
					rank=size(blocks%ButterflyV%blocks(j)%matrix,2)
					allocate (BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix(nn,num_vectors))
					BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix=0
					call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn,BFvec%vec(level_butterfly+1)%blocks(1,j)%matrix,rank,BFvec%vec(level_butterfly+2)%blocks(1,j)%matrix,nn,'N','N',nn,num_vectors,rank,cone,czero,flop=flop)
					flops = flops + flop
				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			else
				flops=0
				!$omp parallel do default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
				do index_ij=1, nr*nc
					index_j_loc = (index_ij-1)/nr+1
					index_i_loc= mod(index_ij-1,nr) + 1  !index_i_loc is local index of column-wise ordering at current level
					index_i=(index_i_loc-1)*inc_r+idx_r  !index_i is global index of column-wise ordering at current level
					index_j=(index_j_loc-1)*inc_c+idx_c

					index_ii=2*index_i-1; index_jj=int((index_j+1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

					index_ii_loc=(index_ii-BFvec%vec(level_butterfly-level+1)%idx_r)/BFvec%vec(level_butterfly-level+1)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
					index_jj_loc=(index_jj-BFvec%vec(level_butterfly-level+1)%idx_c)/BFvec%vec(level_butterfly-level+1)%inc_c+1

					index_i_loc_s=(index_i-BFvec%vec(level_butterfly-level+2)%idx_r)/BFvec%vec(level_butterfly-level+2)%inc_r+1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
					index_i_loc_k=(2*index_i-1-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
					index_j_loc_s=(index_j-BFvec%vec(level_butterfly-level+2)%idx_c)/BFvec%vec(level_butterfly-level+2)%inc_c+1
					index_j_loc_k=(index_j-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

					mm1=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
					mm2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,1)
					nn=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
					allocate (BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix(nn,num_vectors))
					BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix=0


					call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,mm1,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc,index_jj_loc)%matrix,mm1,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm1,cone,cone,flop=flop)
					flops = flops + flop
					call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix,mm2,BFvec%vec(level_butterfly-level+1)%blocks(index_ii_loc+1,index_jj_loc)%matrix,mm2,BFvec%vec(level_butterfly-level+2)%blocks(index_i_loc_s,index_j_loc_s)%matrix,nn,'T','N',nn,num_vectors,mm2,cone,cone,flop=flop)
					flops = flops + flop

				enddo
				!$omp end parallel do
				stats%Flop_Tmp = stats%Flop_Tmp + flops
			endif

			do j=1, BFvec%vec(level_butterfly-level+1)%nc
				do i=1, BFvec%vec(level_butterfly-level+1)%nr
					deallocate (BFvec%vec(level_butterfly-level+1)%blocks(i,j)%matrix)
				enddo
			enddo

			if(level_half+1/=level)then
				call BF_exchange_matvec(blocks,BFvec%vec(level_butterfly-level+2),num_vectors,stats,ptree,level,'C','B')
			endif
		enddo
	endif

    return

end subroutine BF_block_MVP_partial





subroutine Full_block_extraction(blocks,inters,ptree,msh,stats)
    use BPACK_DEFS
    implicit none
	type(matrixblock)::blocks
	type(proctree)::ptree
	type(mesh)::msh
	type(Hstat)::stats
	integer nn,ri,ci,nng,headm,headn,pp,row_group,col_group
	type(intersect)::inters(:)
	integer ii,jj

	headm = msh%basis_group(blocks%row_group)%head
	headn = msh%basis_group(blocks%col_group)%head
	do nn=1,size(blocks%inters,1)
	nng=blocks%inters(nn)%idx
	do ii=1,blocks%inters(nn)%nr_loc
		ri=inters(nng)%rows(blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii)))-headm+1
		do jj=1,blocks%inters(nn)%nc
			ci=inters(nng)%cols(blocks%inters(nn)%cols(jj))-headn+1
			blocks%inters(nn)%dat_loc(ii,jj)=blocks%fullmat(ri,ci)
		enddo
	enddo
	enddo

end subroutine Full_block_extraction



subroutine LR_block_extraction(blocks,inters,ptree,msh,stats)
    use BPACK_DEFS
    implicit none
	type(matrixblock)::blocks
	type(proctree)::ptree
	type(mesh)::msh
	type(Hstat)::stats
	integer nn,ri,ci,nng,headm,headn,head,tail,pp,row_group,col_group
	type(intersect)::inters(:)
	integer ii,jj,rank,ncol,nrow,iidx,pgno,ierr,nr_loc
	DT,allocatable::Vpartial(:,:),matU(:,:)
	real(kind=8)::t1,t2,t3,t4

	headm = blocks%headm
	headn = blocks%headn
	pgno = blocks%pgno

	pp = ptree%myid-ptree%pgrp(pgno)%head+1
	head=blocks%N_p(pp,1)
	tail=blocks%N_p(pp,2)

	t1 = OMP_get_wtime()
	rank=size(blocks%ButterflyU%blocks(1)%matrix,2)
	ncol=0
	do nn=1,size(blocks%inters,1)
		ncol = ncol + blocks%inters(nn)%nc
	enddo
	allocate(Vpartial(rank,ncol))
	t2 = OMP_get_wtime()
	call LR_all2all_extraction(blocks,inters,Vpartial,rank,ncol,stats,ptree,msh)
	t3 = OMP_get_wtime()

	nr_loc=0
	do nn=1,size(blocks%inters,1)
	nr_loc = max(nr_loc,blocks%inters(nn)%nr_loc)
	enddo
	allocate(matU(nr_loc,rank))

	iidx=0
	do nn=1,size(blocks%inters,1)
		nng=blocks%inters(nn)%idx
		if(blocks%inters(nn)%nr_loc>0)then
			do ii=1,blocks%inters(nn)%nr_loc
				ri=inters(nng)%rows(blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii)))-headm+1-blocks%M_p(pp,1)+1
				matU(ii,:)=blocks%ButterflyU%blocks(1)%matrix(ri,:)
			enddo

			call gemmf77('N','N',blocks%inters(nn)%nr_loc,blocks%inters(nn)%nc,rank, cone, matU, nr_loc,Vpartial(1,iidx+1),rank,czero,blocks%inters(nn)%dat_loc(1,1),blocks%inters(nn)%nr_loc)
			stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(blocks%inters(nn)%nr_loc,blocks%inters(nn)%nc,rank)
		endif
		iidx = iidx + blocks%inters(nn)%nc
	enddo

	deallocate(Vpartial)
	deallocate(matU)

	t4 = OMP_get_wtime()
	! time_tmp = time_tmp + t3 - t2

end subroutine LR_block_extraction







!*********** all to all communication of columns in the V factor from the 1D block column layout to that needed by the 1D block row layout
subroutine LR_all2all_extraction(blocks,inters,Vpartial,rank,ncol,stats,ptree,msh)

   use BPACK_DEFS
   implicit none
    integer i, j, level_butterfly, num_blocks, k, attempt,edge_m,edge_n,header_m,header_n,leafsize,nn_start,rankmax_r,rankmax_c,rankmax_min,rank_new
    integer group_m, group_n, mm, nn, index_i, index_i_loc_k,index_i_loc_s, index_j,index_j_loc_k,index_j_loc_s, ii, jj,ij,pp,tt
    integer level, length_1, length_2, level_blocks
    integer level_p, ncol, rank, rankmax, butterflyB_inuse, rank1, rank2
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
    real(kind=8) flop
    DT ctemp
	type(matrixblock)::blocks
	DT::Vpartial(rank,ncol)
	type(Hstat)::stats
	type(proctree)::ptree
	type(mesh)::msh
	type(intersect)::inters(:)
	integer ri,ci,nng,head,tail,row_group,col_group
    integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    integer,allocatable:: select_row_rr(:), select_column_rr(:)
    DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:),core(:,:),tau(:)

	integer,allocatable::jpvt(:)
	integer ierr,nsendrecv,pid,tag,nproc,Nrow,Nreqr,Nreqs,recvid,sendid,tmpi
	integer idx_r,idx_c,inc_r,inc_c,nr,nc,level_new

	type(commquant1D),allocatable::sendquant(:),recvquant(:)
	integer,allocatable::S_req(:),R_req(:)
	integer,allocatable:: statuss(:,:),statusr(:,:)
	character::mode,mode_new
	real(kind=8)::n1,n2,n3,n4
	integer,allocatable::sendIDactive(:),recvIDactive(:)
	integer Nsendactive,Nrecvactive,Nsendactive_min,Nrecvactive_min
	type(butterfly_kerl)::kerls,kerls1
	integer rr,cc
	logical all2all
	integer,allocatable::sdispls(:),sendcounts(:),rdispls(:),recvcounts(:),col_idx_loc(:),activeproc(:,:),row_idx(:,:)
	DT,allocatable::sendbufall2all(:),recvbufall2all(:)
	integer::dist,pgno,ncol_loc,iidx
	integer::headm,headn,nc_loc,ncmax,nrmax,head1,tail1


	n1 = OMP_get_wtime()

	nproc = ptree%pgrp(blocks%pgno)%nproc
	tag = blocks%pgno
	level_p = ptree%nlevel-GetTreelevel(blocks%pgno)
	headm = blocks%headm
	headn = blocks%headn
	pp = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
	head=blocks%N_p(pp,1)
	tail=blocks%N_p(pp,2)

	! allocation of communication quantities
	allocate(statuss(MPI_status_size,nproc))
	allocate(statusr(MPI_status_size,nproc))
	allocate(S_req(nproc))
	allocate(R_req(nproc))
	allocate(sendquant(nproc))
	do ii=1,nproc
		sendquant(ii)%size=0
		sendquant(ii)%active=0
	enddo
	allocate(recvquant(nproc))
	do ii=1,nproc
		recvquant(ii)%size=0
		recvquant(ii)%active=0
	enddo
	allocate(sendIDactive(nproc))
	allocate(recvIDactive(nproc))
	Nsendactive=0
	Nrecvactive=0


	! calculate send buffer sizes in the first pass

	ncmax=0
	nrmax=0
	do nn=1,size(blocks%inters,1)
		ncmax = max(blocks%inters(nn)%nc,ncmax)
		nrmax = max(blocks%inters(nn)%nr,nrmax)
	enddo
	allocate(col_idx_loc(ncmax))
	allocate(row_idx(nrmax,1))
	allocate(activeproc(nproc,size(blocks%inters,1)))
	activeproc=0

	! iidx=0
	do nn=1,size(blocks%inters,1)

		nng=blocks%inters(nn)%idx
		nc_loc=0
		do jj=1,blocks%inters(nn)%nc
			ci=inters(nng)%cols(blocks%inters(nn)%cols(jj))-headn+1
			if(ci>=head .and. ci<=tail)then
				nc_loc=nc_loc+1
				col_idx_loc(nc_loc)=jj
			endif
		enddo

		row_idx(1:blocks%inters(nn)%nr,1)=0
		do ii=1,blocks%inters(nn)%nr
			ri=inters(nng)%rows(blocks%inters(nn)%rows(ii))-headm+1
			row_idx(ii,1)=ri
		enddo
		call PIKSRT_INT_Multi(blocks%inters(nn)%nr,1,row_idx)


		ii=1
		do pp=1,nproc
			head1=blocks%M_p(pp,1)
			tail1=blocks%M_p(pp,2)
			do while(row_idx(ii,1)<head1)
			ii=ii+1
			if(ii>blocks%inters(nn)%nr)exit
			enddo
			if(ii<=blocks%inters(nn)%nr)then
			if(tail1>=row_idx(ii,1))activeproc(pp,nn)=1
			else
			exit
			endif
		enddo

		do pp=1,nproc
			if(activeproc(pp,nn)==1)then
			do jj=1,nc_loc
				ci=inters(nng)%cols(blocks%inters(nn)%cols(col_idx_loc(jj)))-headn+1
				if(sendquant(pp)%active==0)then
					sendquant(pp)%active=1
					Nsendactive=Nsendactive+1
					sendIDactive(Nsendactive)=pp
				endif
				sendquant(pp)%size=sendquant(pp)%size+(1+rank)
			enddo
			endif
		enddo

	enddo
	deallocate(row_idx)


	do nn=1,size(blocks%inters,1)
		nng=blocks%inters(nn)%idx
		if(blocks%inters(nn)%nr_loc>0)then
			do jj=1,blocks%inters(nn)%nc
				ci=inters(nng)%cols(blocks%inters(nn)%cols(jj))
				pgno = findpggroup(ci,msh,ptree,blocks%col_group,blocks%pgno)
				pp = ptree%pgrp(pgno)%head -ptree%pgrp(blocks%pgno)%head+1
				if(recvquant(pp)%active==0)then
					recvquant(pp)%active=1
					Nrecvactive=Nrecvactive+1
					recvIDactive(Nrecvactive)=pp
				endif
				recvquant(pp)%size=recvquant(pp)%size+(1+rank)
			enddo
		endif
	enddo

n2 = OMP_get_wtime()


	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		allocate(sendquant(pp)%dat(sendquant(pp)%size,1))
		sendquant(pp)%size=0
	enddo
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		allocate(recvquant(pp)%dat(recvquant(pp)%size,1))
	enddo

	! pack the send buffer in the second pass
	iidx=0
	do nn=1,size(blocks%inters,1)
		nng=blocks%inters(nn)%idx
		nc_loc=0
		do jj=1,blocks%inters(nn)%nc
			ci=inters(nng)%cols(blocks%inters(nn)%cols(jj))-headn+1
			if(ci>=head .and. ci<=tail)then
				nc_loc=nc_loc+1
				col_idx_loc(nc_loc)=jj
			endif
		enddo

		do pp=1,nproc
			if(activeproc(pp,nn)==1)then
			do jj=1,nc_loc
				ci=inters(nng)%cols(blocks%inters(nn)%cols(col_idx_loc(jj)))-headn+1
				sendquant(pp)%dat(sendquant(pp)%size+1,1)=iidx+col_idx_loc(jj)
				sendquant(pp)%size=sendquant(pp)%size+1
				sendquant(pp)%dat(sendquant(pp)%size+1:sendquant(pp)%size+rank,1)=blocks%ButterflyV%blocks(1)%matrix(ci-head+1,:)
				sendquant(pp)%size=sendquant(pp)%size+rank
			enddo
			endif
		enddo
		iidx = iidx + blocks%inters(nn)%nc
	enddo
	deallocate(col_idx_loc)
	deallocate(activeproc)



n3 = OMP_get_wtime()


	call MPI_ALLREDUCE(Nsendactive,Nsendactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(Nrecvactive,Nrecvactive_min,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(blocks%pgno)%Comm,ierr)
#if 0
	all2all=(nproc==Nsendactive_min .and. Nsendactive_min==Nrecvactive_min)
#else
	all2all=.false.
#endif
	if(all2all)then ! if truly all-to-all, use MPI_ALLTOALLV
		allocate(sdispls(nproc))
		allocate(sendcounts(nproc))
		dist=0
		do pp=1,nproc
			sendcounts(pp)=sendquant(pp)%size
			sdispls(pp)=dist
			dist = dist + sendquant(pp)%size
		enddo
		allocate(sendbufall2all(dist))
		do pp=1,nproc
			if(sendquant(pp)%size>0)sendbufall2all(sdispls(pp)+1:sdispls(pp)+sendcounts(pp))=sendquant(pp)%dat(:,1)
		enddo

		allocate(rdispls(nproc))
		allocate(recvcounts(nproc))
		dist=0
		do pp=1,nproc
			recvcounts(pp)=recvquant(pp)%size
			rdispls(pp)=dist
			dist = dist + recvquant(pp)%size
		enddo
		allocate(recvbufall2all(dist))

		call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_DT, recvbufall2all, recvcounts,rdispls, MPI_DT, ptree%pgrp(blocks%pgno)%Comm, ierr)

		do pp=1,nproc
			if(recvquant(pp)%size>0)recvquant(pp)%dat(:,1) = recvbufall2all(rdispls(pp)+1:rdispls(pp)+recvcounts(pp))
		enddo

		deallocate(sdispls)
		deallocate(sendcounts)
		deallocate(sendbufall2all)
		deallocate(rdispls)
		deallocate(recvcounts)
		deallocate(recvbufall2all)
	else
		Nreqs=0
		do tt=1,Nsendactive
			pp=sendIDactive(tt)
			recvid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(recvid/=ptree%MyID)then
				Nreqs=Nreqs+1
				call MPI_Isend(sendquant(pp)%dat,sendquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,S_req(Nreqs),ierr)
			else
				if(sendquant(pp)%size>0)recvquant(pp)%dat=sendquant(pp)%dat
			endif
		enddo

		Nreqr=0
		do tt=1,Nrecvactive
			pp=recvIDactive(tt)
			sendid=pp-1+ptree%pgrp(blocks%pgno)%head
			if(sendid/=ptree%MyID)then
				Nreqr=Nreqr+1
				call MPI_Irecv(recvquant(pp)%dat,recvquant(pp)%size,MPI_DT,pp-1,tag+1,ptree%pgrp(blocks%pgno)%Comm,R_req(Nreqr),ierr)
			endif
		enddo

		if(Nreqs>0)then
			call MPI_waitall(Nreqs,S_req,statuss,ierr)
		endif
		if(Nreqr>0)then
			call MPI_waitall(Nreqr,R_req,statusr,ierr)
		endif
	endif

	! copy data from buffer to target
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		i=0
		do while(i<recvquant(pp)%size)
			i=i+1
			iidx=NINT(dble(recvquant(pp)%dat(i,1)))
			Vpartial(:,iidx) = recvquant(pp)%dat(i+1:i+rank,1)
			i=i+rank
		enddo
	enddo

n4 = OMP_get_wtime()
	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do tt=1,Nsendactive
		pp=sendIDactive(tt)
		if(allocated(sendquant(pp)%dat))deallocate(sendquant(pp)%dat)
	enddo
	deallocate(sendquant)
	do tt=1,Nrecvactive
		pp=recvIDactive(tt)
		if(allocated(recvquant(pp)%dat))deallocate(recvquant(pp)%dat)
	enddo
	deallocate(recvquant)
	deallocate(sendIDactive)
	deallocate(recvIDactive)



	! time_tmp = time_tmp + n4 - n1

end subroutine LR_all2all_extraction









subroutine BF_block_extraction(blocks,inters,ptree,msh,stats)
    use BPACK_DEFS
    implicit none
	type(matrixblock)::blocks
	type(proctree)::ptree
	type(mesh)::msh
	type(Hstat)::stats
	integer ri,ci,nng,headm,headn,head,tail,pp,row_group,col_group
	type(intersect)::inters(:)
	integer ii,iii, jj,jjj,nnn,gg,gg1,gg2,rank,ncol,nrow,iidx,jidx,pgno,comm,ierr,idx1,idx2,idx,idx_r,idx_c,idx_r0,idx_c0,nc,nc0,nr,nr0,inc_c,inc_c0,inc_r,inc_r0
    integer level_butterfly, num_blocks, level_half
    integer group_m, group_n, mm, nn, index_i,index_i0, index_i_k, index_j_k, index_i_loc_k,index_i_s, index_j_s, index_i_loc_s, index_j,index_j0,index_j_loc_k,index_j_loc_s,ij,tt,nvec1,nvec2
    integer level
    real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
	integer header_n1,header_n2,nn1,nn2,mmm,index_ii,index_jj,index_ii_loc,index_jj_loc,nnn1
	integer nsendrecv,pid,pid0,tag,nproc,Nreqr,Nreqs,recvid,sendid
	type(butterfly_vec) :: BFvec,BFvec1
	type(butterfly_kerl),allocatable::g_idx_m(:),g_idx_n(:)
	integer,allocatable::group_ms(:),group_ns(:),group_ms1(:),group_ns1(:)
	type(nod),pointer::cur
	class(*),pointer::ptr
	type(ipair):: p
	DT,allocatable::mat1(:,:),mat2(:,:),mat(:,:)
	DT,allocatable::Vpartial(:,:)
	DT::val
	integer idxr,idxc
	integer,allocatable::num_nods_i(:),num_nods_j(:)
	real(kind=8)::n1,n2,n3,n4,n5


	n1 = OMP_get_wtime()

	headm = blocks%headm
	headn = blocks%headn
	pgno = blocks%pgno
	comm = ptree%pgrp(pgno)%comm

	pp = ptree%myid-ptree%pgrp(pgno)%head+1
	head=blocks%N_p(pp,1)
	tail=blocks%N_p(pp,2)

	level_butterfly=blocks%level_butterfly
	num_blocks=2**level_butterfly
	level_half = blocks%level_half

	! write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,'inin'

	!******* preallocate BFvec and BFvec1, #of blocks are exactly those in BF_block_MVP_dat, BFvec for the first 0:level_half+1 levels and BFvec for the next level_half+1 to level_butterfly+2 levels,note that level level_half+1 is duplicated for all2all communication
	allocate(BFvec%vec(0:level_half+1))
	allocate (BFvec%vec(0)%blocks(1,blocks%ButterflyV%nblk_loc))
	BFvec%vec(0)%num_row=1
	BFvec%vec(0)%num_col=num_blocks
	BFvec%vec(0)%idx_r=1
	BFvec%vec(0)%inc_r=1
	BFvec%vec(0)%nr=1
	BFvec%vec(0)%idx_c=blocks%ButterflyV%idx
	BFvec%vec(0)%inc_c=blocks%ButterflyV%inc
	BFvec%vec(0)%nc=blocks%ButterflyV%nblk_loc

	do level=0, level_half
		call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')

		BFvec%vec(level+1)%idx_r=idx_r
		BFvec%vec(level+1)%inc_r=inc_r
		BFvec%vec(level+1)%nr=nr
		BFvec%vec(level+1)%idx_c=idx_c
		BFvec%vec(level+1)%inc_c=inc_c
		BFvec%vec(level+1)%nc=nc
		if (level/=level_butterfly+1) then
			BFvec%vec(level+1)%num_row=2**level
			BFvec%vec(level+1)%num_col=2**(level_butterfly-level)
		else
			BFvec%vec(level+1)%num_row=2**level_butterfly
			BFvec%vec(level+1)%num_col=1
		endif
		if(level_half/=level)then ! the last level doesn't require doubling block columns
		if(nc==1 .and. 2**(level_butterfly-level)>1)then ! double the number of local block columns used for MPI communication
			BFvec%vec(level+1)%nc=2
			BFvec%vec(level+1)%idx_c=BFvec%vec(level+1)%idx_c-1+mod(BFvec%vec(level+1)%idx_c,2)
		endif
		endif
		allocate(BFvec%vec(level+1)%blocks(BFvec%vec(level+1)%nr,BFvec%vec(level+1)%nc))
	enddo


	allocate(BFvec1%vec(level_half+1:level_butterfly+2))
	do level=level_half,level_butterfly+1
		if(level==level_half)then
			call GetLocalBlockRange(ptree,blocks%pgno,level_half+1,level_butterfly,idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0,'C') ! this is the same as the target block range used in BF_all2all_matvec
			idx_r=idx_r0
			nr=nr0
			inc_r=inc_r0
			idx_c=idx_c0
			nc=nc0
			inc_c=inc_c0
		else
			call GetLocalBlockRange(ptree,blocks%pgno,level,level_butterfly,idx_r0,inc_r0,nr0,idx_c0,inc_c0,nc0,'C')
			! convert the local column-wise kernel block ranges to local row-wise output vector ranges
			if(level/=0 .and. level/=level_butterfly+1)then
				idx_r = idx_r0*2-1
				nr = nr0*2
				inc_r=inc_r0
				idx_c = ceiling_safe(idx_c0/2d0)
				if(inc_c0>1)then
				nc=nc0
				else
				nc = ceiling_safe(nc0/2d0)
				endif
				inc_c = ceiling_safe(inc_c0/2d0)
			else
				idx_r=idx_r0
				nr=nr0
				inc_r=inc_r0
				idx_c=idx_c0
				nc=nc0
				inc_c=inc_c0
			endif
		endif

		BFvec1%vec(level+1)%idx_r=idx_r
		BFvec1%vec(level+1)%inc_r=inc_r
		BFvec1%vec(level+1)%nr=nr
		BFvec1%vec(level+1)%idx_c=idx_c
		BFvec1%vec(level+1)%inc_c=inc_c
		BFvec1%vec(level+1)%nc=nc
		if (level/=level_butterfly+1) then
			BFvec1%vec(level+1)%num_row=2**level
			BFvec1%vec(level+1)%num_col=2**(level_butterfly-level)
		else
			BFvec1%vec(level+1)%num_row=2**level_butterfly
			BFvec1%vec(level+1)%num_col=1
		endif
		allocate(BFvec1%vec(level+1)%blocks(BFvec1%vec(level+1)%nr,BFvec1%vec(level+1)%nc))
	enddo


	!**** compute group_ms and group_ns which stores the leaf block number (from 1 to 2^L) of each index. group_ms1 and group_ns1 are used to track the block number in each butterfly level.
	nrow=0
	ncol=0
	do nn=1,size(blocks%inters,1)
	nng=blocks%inters(nn)%idx
	nrow = nrow + blocks%inters(nn)%nr
	ncol = ncol + blocks%inters(nn)%nc
	enddo
	allocate(group_ms(nrow))
	allocate(group_ns(ncol))
	iidx=0
	jidx=0
	do nn=1,size(blocks%inters,1)
	nng=blocks%inters(nn)%idx
	do ii=1,blocks%inters(nn)%nr
		iidx = iidx + 1
		group_ms(iidx)=findgroup(inters(nng)%rows(blocks%inters(nn)%rows(ii)),msh,level_butterfly,blocks%row_group)-blocks%row_group*2**level_butterfly+1
	enddo
	do jj=1,blocks%inters(nn)%nc
		jidx = jidx + 1
		group_ns(jidx)=findgroup(inters(nng)%cols(blocks%inters(nn)%cols(jj)),msh,level_butterfly,blocks%col_group)-blocks%col_group*2**level_butterfly+1
	enddo
	enddo
	allocate(group_ms1(nrow))
	group_ms1=group_ms
	allocate(group_ns1(ncol))
	group_ns1=group_ns


	!**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
	allocate(g_idx_m(0:level_half+1))
	allocate(g_idx_n(0:level_half+1))
	do level=0, level_half+1
		g_idx_m(level)%idx_r=BFvec%vec(level)%idx_r
		g_idx_m(level)%inc_r=BFvec%vec(level)%inc_r
		g_idx_m(level)%nr=BFvec%vec(level)%nr
		g_idx_m(level)%idx_c=1
		g_idx_m(level)%inc_c=1
		g_idx_m(level)%nc=1
		allocate(g_idx_m(level)%blocks(g_idx_m(level)%nr,g_idx_m(level)%nc))
		allocate(g_idx_m(level)%index(iidx,1))

		g_idx_n(level)%idx_c=BFvec%vec(level)%idx_c
		g_idx_n(level)%inc_c=BFvec%vec(level)%inc_c
		g_idx_n(level)%nc=BFvec%vec(level)%nc
		g_idx_n(level)%idx_r=1
		g_idx_n(level)%inc_r=1
		g_idx_n(level)%nr=1
		allocate(g_idx_n(level)%blocks(g_idx_n(level)%nr,g_idx_n(level)%nc))
		allocate(g_idx_n(level)%index(jidx,1))
	enddo

	iidx=0
	jidx=0
	allocate(num_nods_i(0:level_half+1))
	allocate(num_nods_j(0:level_half+1))
	do nn=1,size(blocks%inters,1)
		!*** for each index in each intersection, traverse the row and column tree
		num_nods_i=0
		num_nods_j=0
		do level=level_butterfly+2,0,-1
			do ii=1,blocks%inters(nn)%nr

				if(level<=level_half+1)then
				index_i_s = group_ms1(iidx+ii)
				index_i_loc_s = (index_i_s-g_idx_m(level)%idx_r)/g_idx_m(level)%inc_r+1
				if(index_i_s>=g_idx_m(level)%idx_r .and. mod(index_i_s-g_idx_m(level)%idx_r,g_idx_m(level)%inc_r)==0 .and. index_i_loc_s<=g_idx_m(level)%nr)then
					if(g_idx_m(level)%blocks(index_i_loc_s,1)%lst%idx/=nn)then ! not seen this index_i_loc_s before
						num_nods_i(level) = num_nods_i(level)+1
						g_idx_m(level)%index(num_nods_i(level),1)=index_i_loc_s
						g_idx_m(level)%blocks(index_i_loc_s,1)%lst%idx=nn
					endif
				endif
				endif

			enddo
			if(level>0 .and. level<level_butterfly+2)group_ms1(iidx+1:iidx+blocks%inters(nn)%nr) = floor((group_ms1(iidx+1:iidx+blocks%inters(nn)%nr)+1)/2d0)
		enddo
		iidx = iidx + blocks%inters(nn)%nr

		do level=0, level_butterfly+2
			do jj=1,blocks%inters(nn)%nc
				if(level<=level_half+1)then
				index_j_s = group_ns1(jidx+jj)
				index_j_loc_s = (index_j_s-g_idx_n(level)%idx_c)/g_idx_n(level)%inc_c+1
				if(index_j_s>=g_idx_n(level)%idx_c .and. mod(index_j_s-g_idx_n(level)%idx_c,g_idx_n(level)%inc_c)==0 .and. index_j_loc_s<=g_idx_n(level)%nc)then
					if(g_idx_n(level)%blocks(1,index_j_loc_s)%lst%idx/=nn)then ! not seen this index_j_loc_s before
						num_nods_j(level) = num_nods_j(level)+1
						g_idx_n(level)%index(num_nods_j(level),1)=index_j_loc_s
						g_idx_n(level)%blocks(1,index_j_loc_s)%lst%idx=nn
					endif
					if(level==0)BFvec%vec(level)%blocks(1,index_j_loc_s)%ndim = BFvec%vec(level)%blocks(1,index_j_loc_s)%ndim+1
				endif
				endif

			enddo
			if(level>0 .and. level<level_butterfly+2)group_ns1(jidx+1:jidx+blocks%inters(nn)%nc) = floor((group_ns1(jidx+1:jidx+blocks%inters(nn)%nc)+1)/2d0)

		enddo
		jidx = jidx + blocks%inters(nn)%nc


		!**** construct BFvec%vec(level)%lst for active BF blocks and BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst for list of intersection#s
		do level=0, level_half+1
			do ii=1,num_nods_i(level)
			index_i_loc_s=g_idx_m(level)%index(ii,1)
				do jj=1,num_nods_j(level)
				index_j_loc_s=g_idx_n(level)%index(jj,1)
					if(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods==0)then ! first see this block
						p%i=index_i_loc_s
						p%j=index_j_loc_s
						! write(*,*)p%i,p%j,BFvec%vec(level)%lst%num_nods
						call append(BFvec%vec(level)%lst,p)
					endif
					call append(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst,nn)
				enddo
			enddo
		enddo

	enddo

	do level=0, level_half+1
		deallocate(g_idx_m(level)%index)
		deallocate(g_idx_n(level)%index)
	enddo
	deallocate(g_idx_m)
	deallocate(g_idx_n)
	deallocate(group_ms1)
	deallocate(group_ns1)
	deallocate(num_nods_i)
	deallocate(num_nods_j)


	!**** copy *%lst to *%index
	do level=0, level_half+1
		allocate(BFvec%vec(level)%index(BFvec%vec(level)%lst%num_nods,2))
		cur=>BFvec%vec(level)%lst%head
		do nn=1,BFvec%vec(level)%lst%num_nods
		select type(ptr=>cur%item)
		type is(ipair)
			BFvec%vec(level)%index(nn,1)=ptr%i
			BFvec%vec(level)%index(nn,2)=ptr%j
		end select
		cur=>cur%next
		enddo
		call list_finalizer(BFvec%vec(level)%lst)

		do nn=1,size(BFvec%vec(level)%index,1)
			index_i_loc_s=BFvec%vec(level)%index(nn,1)
			index_j_loc_s=BFvec%vec(level)%index(nn,2)
			! if(ptree%MyID>=2)write(*,*)ptree%MyID,'my',blocks%row_group,blocks%col_group,'BFvec range',level,index_i_loc_s,index_j_loc_s,(index_i_loc_s-1)*BFvec%vec(level)%inc_r+BFvec%vec(level)%idx_r,(index_j_loc_s-1)*BFvec%vec(level)%inc_c+BFvec%vec(level)%idx_c
			allocate(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods,1))
			cur=>BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%head
			do ii=1,BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods
			select type(ptr=>cur%item)
			type is(integer)
				BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,1)=ptr
			end select
			cur=>cur%next
			enddo
			call list_finalizer(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst)
		enddo
	enddo


	allocate(group_ms1(nrow))
	group_ms1=group_ms
	allocate(group_ns1(ncol))
	group_ns1=group_ns

	!**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
	allocate(g_idx_m(level_half+1:level_butterfly+2))
	allocate(g_idx_n(level_half+1:level_butterfly+2))
	do level=level_half+1, level_butterfly+2
		g_idx_m(level)%idx_r=BFvec1%vec(level)%idx_r
		g_idx_m(level)%inc_r=BFvec1%vec(level)%inc_r
		g_idx_m(level)%nr=BFvec1%vec(level)%nr
		g_idx_m(level)%idx_c=1
		g_idx_m(level)%inc_c=1
		g_idx_m(level)%nc=1
		allocate(g_idx_m(level)%blocks(g_idx_m(level)%nr,g_idx_m(level)%nc))
		allocate(g_idx_m(level)%index(iidx,1))

		g_idx_n(level)%idx_c=BFvec1%vec(level)%idx_c
		g_idx_n(level)%inc_c=BFvec1%vec(level)%inc_c
		g_idx_n(level)%nc=BFvec1%vec(level)%nc
		g_idx_n(level)%idx_r=1
		g_idx_n(level)%inc_r=1
		g_idx_n(level)%nr=1
		allocate(g_idx_n(level)%blocks(g_idx_n(level)%nr,g_idx_n(level)%nc))
		allocate(g_idx_n(level)%index(jidx,1))
	enddo


	!**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
	iidx=0
	jidx=0
	allocate(num_nods_i(level_half+1:level_butterfly+2))
	allocate(num_nods_j(level_half+1:level_butterfly+2))
	do nn=1,size(blocks%inters,1)
		num_nods_i=0
		num_nods_j=0
		!*** for each index in each intersection, traverse the row and column tree
		do level=level_butterfly+2,0,-1
			do ii=1,blocks%inters(nn)%nr


				if(level>=level_half+1)then
				index_i_s = group_ms1(iidx+ii)
				index_i_loc_s = (index_i_s-g_idx_m(level)%idx_r)/g_idx_m(level)%inc_r+1
				if(index_i_s>=g_idx_m(level)%idx_r .and. mod(index_i_s-g_idx_m(level)%idx_r,g_idx_m(level)%inc_r)==0 .and. index_i_loc_s<=g_idx_m(level)%nr)then
					if(g_idx_m(level)%blocks(index_i_loc_s,1)%lst%idx/=nn)then ! not seen this index_i_loc_s before
						num_nods_i(level) = num_nods_i(level)+1
						g_idx_m(level)%index(num_nods_i(level),1)=index_i_loc_s
						g_idx_m(level)%blocks(index_i_loc_s,1)%lst%idx=nn
					endif
				endif
				endif

			enddo
			if(level>0 .and. level<level_butterfly+2)group_ms1(iidx+1:iidx+blocks%inters(nn)%nr) = floor((group_ms1(iidx+1:iidx+blocks%inters(nn)%nr)+1)/2d0)
		enddo
		iidx = iidx + blocks%inters(nn)%nr

		do level=0, level_butterfly+2
			do jj=1,blocks%inters(nn)%nc
				if(level>=level_half+1)then
				index_j_s = group_ns1(jidx+jj)
				index_j_loc_s = (index_j_s-g_idx_n(level)%idx_c)/g_idx_n(level)%inc_c+1
				if(index_j_s>=g_idx_n(level)%idx_c .and. mod(index_j_s-g_idx_n(level)%idx_c,g_idx_n(level)%inc_c)==0 .and. index_j_loc_s<=g_idx_n(level)%nc)then
					if(g_idx_n(level)%blocks(1,index_j_loc_s)%lst%idx/=nn)then ! not seen this index_j_loc_s before
						num_nods_j(level) = num_nods_j(level)+1
						g_idx_n(level)%index(num_nods_j(level),1)=index_j_loc_s
						g_idx_n(level)%blocks(1,index_j_loc_s)%lst%idx=nn

					endif
					if(level==0)BFvec1%vec(level)%blocks(1,index_j_loc_s)%ndim = BFvec1%vec(level)%blocks(1,index_j_loc_s)%ndim+1
				endif
				endif

			enddo
			if(level>0 .and. level<level_butterfly+2)group_ns1(jidx+1:jidx+blocks%inters(nn)%nc) = floor((group_ns1(jidx+1:jidx+blocks%inters(nn)%nc)+1)/2d0)
		enddo
		jidx = jidx + blocks%inters(nn)%nc

		!**** construct BFvec1%vec(level)%lst for active BF blocks and BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst for list of intersection#s
		do level=level_half+1, level_butterfly+2
			do ii=1,num_nods_i(level)
			index_i_loc_s=g_idx_m(level)%index(ii,1)
				do jj=1,num_nods_j(level)
				index_j_loc_s=g_idx_n(level)%index(jj,1)
					if(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods==0)then ! first see this block
						p%i=index_i_loc_s
						p%j=index_j_loc_s
						call append(BFvec1%vec(level)%lst,p)
					endif
					call append(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst,nn)
				enddo
			enddo
		enddo
	enddo
	do level=level_half+1, level_butterfly+2
		deallocate(g_idx_m(level)%index)
		deallocate(g_idx_n(level)%index)
	enddo
	deallocate(g_idx_m)
	deallocate(g_idx_n)
	deallocate(group_ms1)
	deallocate(group_ns1)
	deallocate(num_nods_i)
	deallocate(num_nods_j)

	!**** copy *%lst to *%index
	do level=level_half+1, level_butterfly+2
		allocate(BFvec1%vec(level)%index(BFvec1%vec(level)%lst%num_nods,2))
		cur=>BFvec1%vec(level)%lst%head
		do nn=1,BFvec1%vec(level)%lst%num_nods
		select type(ptr=>cur%item)
		type is(ipair)
			BFvec1%vec(level)%index(nn,1)=ptr%i
			BFvec1%vec(level)%index(nn,2)=ptr%j
		end select
		cur=>cur%next
		enddo
		call list_finalizer(BFvec1%vec(level)%lst)

		do nn=1,size(BFvec1%vec(level)%index,1)
			index_i_loc_s=BFvec1%vec(level)%index(nn,1)
			index_j_loc_s=BFvec1%vec(level)%index(nn,2)

! if(ptree%MyID>=2)write(*,*)ptree%MyID,'my',blocks%row_group,blocks%col_group,'BFvec1 range',level,index_i_loc_s,index_j_loc_s,(index_i_loc_s-1)*BFvec1%vec(level)%inc_r+BFvec1%vec(level)%idx_r,(index_j_loc_s-1)*BFvec1%vec(level)%inc_c+BFvec1%vec(level)%idx_c


			if(level/=level_butterfly+2)then ! last level doesn't need a list of intersection#
			allocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods,1))
			cur=>BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%head
			do ii=1,BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods
			select type(ptr=>cur%item)
			type is(integer)
				BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,1)=ptr
			end select
			cur=>cur%next
			enddo
			endif
			call list_finalizer(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst)
		enddo
	enddo



	!**** create a list of row_loc indices for each block of BFvec1%vec(level_butterfly+2),note that BFvec1%vec(level)%blocks(index_i_loc_s,1)%lst and BFvec1%vec(level)%blocks(index_i_loc_s,1)%index are used at this level compared to other levels
	iidx=0
	do nn=1,size(blocks%inters,1)
		!*** for each index in each intersection, traverse the row and column tree
		do ii=1,blocks%inters(nn)%nr
			iidx = iidx + 1
			level=level_butterfly+2
			index_i_s = group_ms(iidx)
			index_i_loc_s = (index_i_s-BFvec1%vec(level)%idx_r)/BFvec1%vec(level)%inc_r+1
			if(index_i_s>=BFvec1%vec(level)%idx_r .and. mod(index_i_s-BFvec1%vec(level)%idx_r,BFvec1%vec(level)%inc_r)==0 .and. index_i_loc_s<=BFvec1%vec(level)%nr)then
				p%i=nn
				p%j=blocks%inters(nn)%glo2loc(ii)
				call append(BFvec1%vec(level)%blocks(index_i_loc_s,1)%lst,p)
				! write(*,*)ptree%MyID,'idtarg',index_i_s,index_i_loc_s
			endif
		enddo
	enddo

	!**** copy *%lst to *%index
	level=level_butterfly+2
	do nn=1,size(BFvec1%vec(level)%index,1)
		index_i_loc_s=BFvec1%vec(level)%index(nn,1)
		index_j_loc_s=BFvec1%vec(level)%index(nn,2)

		allocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods,2))
		cur=>BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%head
		do ii=1,BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst%num_nods
		select type(ptr=>cur%item)
		type is(ipair)
			BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,1)=ptr%i
			BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,2)=ptr%j
		end select
		cur=>cur%next
		enddo
		call list_finalizer(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst)
	enddo

	n2 = OMP_get_wtime()


	!**** generate data in BFvec%vec(0)
	do nn=1,size(BFvec%vec(0)%index,1)
		index_j_loc_s = BFvec%vec(0)%index(nn,2)
		allocate(BFvec%vec(0)%blocks(1,index_j_loc_s)%matrix(3,BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim))
		BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim=0
	enddo
	jidx=0
	do nn=1,size(blocks%inters,1)
		nng=blocks%inters(nn)%idx
		do jj=1,blocks%inters(nn)%nc
			jidx = jidx + 1
			index_j_s = group_ns(jidx)
			index_j_loc_s = (index_j_s-BFvec%vec(0)%idx_c)/BFvec%vec(0)%inc_c+1
			if(index_j_s>=BFvec%vec(0)%idx_c .and. mod(index_j_s-BFvec%vec(0)%idx_c,BFvec%vec(0)%inc_c)==0 .and. index_j_loc_s<=BFvec%vec(0)%nc)then
				BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim = BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim + 1
				BFvec%vec(0)%blocks(1,index_j_loc_s)%matrix(1,BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim)=nn
				BFvec%vec(0)%blocks(1,index_j_loc_s)%matrix(2,BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim)=jj
				BFvec%vec(0)%blocks(1,index_j_loc_s)%matrix(3,BFvec%vec(0)%blocks(1,index_j_loc_s)%ndim)=inters(nng)%cols(blocks%inters(nn)%cols(jj))-msh%basis_group(blocks%col_group*2**level_butterfly+index_j_s-1)%head+1
			endif
		enddo
	enddo


	!**** multiply BF with BFvec
	do level=0, level_half
		do nn=1,size(BFvec%vec(level+1)%index,1)
			index_i_loc_s = BFvec%vec(level+1)%index(nn,1)
			index_i_s = (index_i_loc_s-1)*BFvec%vec(level+1)%inc_r+BFvec%vec(level+1)%idx_r
			index_j_loc_s = BFvec%vec(level+1)%index(nn,2)
			index_j_s = (index_j_loc_s-1)*BFvec%vec(level+1)%inc_c+BFvec%vec(level+1)%idx_c
			call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i_s,index_j_s,'R',pid)
			if(pid==ptree%MyID)then
				if(level==0)then
					index_j = index_j_s
					index_j_loc_k = (index_j-blocks%ButterflyV%idx)/blocks%ButterflyV%inc+1
					index_jj_loc = index_j_loc_k
					rank=size(blocks%ButterflyV%blocks(index_j_loc_k)%matrix,2)
					allocate(BFvec%vec(level+1)%blocks(1,index_j_loc_s)%matrix(rank+2,BFvec%vec(level)%blocks(1,index_jj_loc)%ndim))
					BFvec%vec(level+1)%blocks(1,index_j_loc_s)%matrix(1:2,:) = BFvec%vec(level)%blocks(1,index_jj_loc)%matrix(1:2,:)
					do nnn=1,BFvec%vec(level)%blocks(1,index_jj_loc)%ndim
						jjj =NINT(dble(BFvec%vec(level)%blocks(1,index_jj_loc)%matrix(3,nnn)))
						BFvec%vec(level+1)%blocks(1,index_j_loc_s)%matrix(3:rank+2,nnn) =  blocks%ButterflyV%blocks(index_j_loc_k)%matrix(jjj,:)
					enddo
				else if(level==level_butterfly+1)then
					write(*,*)'should not come here as level_half<=level_butterfly'
					stop
				else
					index_i = index_i_s
					index_j = index_j_s
					index_ii=floor_safe((index_i+1)/2d0) ; index_jj=2*index_j-1 !index_ii is global index in BFvec%vec(level)

					index_ii_loc=(index_ii-BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r+1 !index_ii_loc is local index in BFvec%vec(level)
					index_jj_loc=(index_jj-BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c+1

					index_i_loc_k=(index_i-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
					index_j_loc_k=(2*index_j-1-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

					nn1=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
					nn2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix,2)
					mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
					if(allocated(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix))then
					nvec1=size(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,2)
					else
					nvec1=0
					endif

					if(allocated(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix))then
					nvec2=size(BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix,2)
					else
					nvec2=0
					endif

					allocate(mat1(mm+2,nvec1))
					if(nvec1>0)then
					mat1=0
					mat1(1:2,:)=BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix(1:2,:)
					call gemmf77('N','N',mm,nvec1,nn1, cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix, mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix(3,1),nn1+2,czero,mat1(3,1),mm+2)
					stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm,nvec1,nn1)
					endif

					allocate(mat2(mm+2,nvec2))
					if(nvec2>0)then
					mat2=0
					mat2(1:2,:)=BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix(1:2,:)
					call gemmf77('N','N',mm,nvec2,nn2, cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k+1)%matrix, mm,BFvec%vec(level)%blocks(index_ii_loc,index_jj_loc+1)%matrix(3,1),nn2+2,czero,mat2(3,1),mm+2)
					stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm,nvec2,nn2)
					endif

					!**** filter out the columns of mat1 and mat2 that are not specified by BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index
					allocate(mat(mm+2,nvec1+nvec2))
					idx1=0
					idx2=0
					idx=0
					gg1=0
					gg2=0
					do ii = 1, size(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index,1)
						gg=NINT(dble(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,1)))
						do while(gg1<=gg)
							if(gg1==gg)then
								idx=idx+1
								mat(:,idx)=mat1(:,idx1)
							endif
							idx1=idx1+1
							if(idx1>nvec1)exit
							gg1=NINT(dble(mat1(1,idx1)))
						enddo
						do while(gg2<=gg)
							if(gg2==gg)then
								idx=idx+1
								mat(:,idx)=mat2(:,idx2)
							endif
							idx2=idx2+1
							if(idx2>nvec2)exit
							gg2=NINT(dble(mat2(1,idx2)))
						enddo
					enddo
					allocate(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm+2,idx))
					BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix = mat(:,1:idx)
					deallocate(mat1)
					deallocate(mat2)
					deallocate(mat)

				endif
			endif

		enddo

		do nn=1,size(BFvec%vec(level)%index,1)
			index_i_loc_s = BFvec%vec(level)%index(nn,1)
			index_j_loc_s = BFvec%vec(level)%index(nn,2)
			if(allocated(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index))deallocate(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index)
			if(allocated(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix))deallocate(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix)
		enddo
		if(level_half/=level)then
		call BF_exchange_extraction(blocks,BFvec%vec(level+1),stats,ptree,level,'B')
		endif
	enddo

	n3 = OMP_get_wtime()

	! write(*,*)blocks%row_group,blocks%col_group,'myid',ptree%MyID,'level',level,'before all2all'
	! all2all communication from BFvec%vec(level_half+1) to BFvec1%vec(level_half+1)
	call BF_all2all_extraction(blocks,BFvec%vec(level_half+1),BFvec1%vec(level_half+1),stats,ptree,level_half,'R','C')
	! write(*,*)blocks%row_group,blocks%col_group,'myid',ptree%MyID,'level',level,'after all2all'

	n4 = OMP_get_wtime()

	!**** multiply BF with BFvec1
	do level=level_half+1,level_butterfly+1
		if (level==0) then
			write(*,*)'should not come here as level_half>=0'
			stop
		elseif (level==level_butterfly+1) then

			do nn=1,size(BFvec1%vec(level+1)%index,1)
				index_i_loc_s = BFvec1%vec(level+1)%index(nn,1)
				index_j_loc_s = BFvec1%vec(level+1)%index(nn,2)

				idxr=0
				idxc=0

				do while(idxr<size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index,1))
					idx=BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+1,1)
					nr=0
					do while(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+1+nr,1)==idx)
					nr=nr+1
					if(idxr+1+nr>size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index,1))exit
					enddo
					! if(ptree%MyID>=2)write(*,*)'id',ptree%MyID,index_i_loc_s,index_j_loc_s,idx,'shape',shape(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix),idxc
					idx1=NINT(dble(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(1,idxc+1)))
					call assert(idx==idx1,'row and column intersection# not match')
					nc=0
					do while(NINT(dble(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(1,idxc+1+nc)))==idx)
					nc=nc+1
					if(idxc+1+nc>size(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix,2))exit
					enddo


					rank = size(blocks%ButterflyU%blocks(index_i_loc_s)%matrix,2)
					allocate(Vpartial(nr,nc))
					Vpartial=0
					allocate(mat(nr,rank))

					nng=blocks%inters(idx)%idx
					do ii=1,nr
						iii = BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+ii,2)
						mmm = inters(nng)%rows(blocks%inters(idx)%rows(blocks%inters(idx)%rows_loc(iii)))
						ri=mmm-msh%basis_group(findgroup(mmm,msh,level_butterfly,blocks%row_group))%head+1
						mat(ii,:) = blocks%ButterflyU%blocks(index_i_loc_s)%matrix(ri,:)
					enddo

					call gemmf77('N','N',nr,nc,rank, cone, mat, nr,BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(3,idxc+1),rank+2,czero,Vpartial,nr)
					stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(nr,nc,rank)

					do ii=1,nr
						iii = BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+ii,2)
						do jj=1,nc
							jjj = NINT(dble(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(2,idxc+jj)))
							blocks%inters(idx)%dat_loc(iii,jjj) = Vpartial(ii,jj)
						enddo
					enddo

					! do ii=1,nr
						! iii = BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+ii,2)
						! iii=inters(nng)%rows(blocks%inters(idx)%rows(blocks%inters(idx)%rows_loc(iii)))-msh%basis_group(blocks%row_group)%head+1
						! do jj=1,nc
							! jjj = BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(2,idxc+jj)
							! jjj=inters(nng)%cols(blocks%inters(idx)%cols(jjj))-msh%basis_group(blocks%col_group)%head+1
							! call BF_value(iii,jjj,blocks,val)
							! write(*,*)'seq:',abs(val)
						! enddo
					! enddo

					deallocate(Vpartial)
					deallocate(mat)

					idxr = idxr + nr
					idxc = idxc + nc
				enddo
			enddo


			do nn=1,size(BFvec1%vec(level)%index,1)
				index_i_loc_s = BFvec1%vec(level)%index(nn,1)
				index_j_loc_s = BFvec1%vec(level)%index(nn,2)
				if(allocated(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index))deallocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index)
				if(allocated(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix))deallocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix)
			enddo

			do nn=1,size(BFvec1%vec(level+1)%index,1)
				index_i_loc_s = BFvec1%vec(level+1)%index(nn,1)
				index_j_loc_s = BFvec1%vec(level+1)%index(nn,2)
				if(allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index))deallocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index)
				if(allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix))deallocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)
			enddo


		else
			do nn=1,size(BFvec1%vec(level+1)%index,1)
				index_i_loc_s = BFvec1%vec(level+1)%index(nn,1)
				index_i_s = (index_i_loc_s-1)*BFvec1%vec(level+1)%inc_r+BFvec1%vec(level+1)%idx_r
				index_j_loc_s = BFvec1%vec(level+1)%index(nn,2)
				index_j_s = (index_j_loc_s-1)*BFvec1%vec(level+1)%inc_c+BFvec1%vec(level+1)%idx_c
				do jjj=0,1
					index_i=floor_safe((index_i_s-1)/2d0)+1
					index_j=2*index_j_s-jjj
					call GetBlockPID(ptree,blocks%pgno,level,level_butterfly,index_i,index_j,'C',pid)
					if(pid==ptree%MyID)then
						index_i_k = 2*index_i-mod(index_i_s,2)
						index_j_k = index_j

						index_i_loc_k=(index_i_k-blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r+1 !index_i_loc_k is local index of kernels at current level
						index_j_loc_k=(index_j_k-blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c+1

						index_ii_loc=(index_i-BFvec1%vec(level)%idx_r)/BFvec1%vec(level)%inc_r+1 !index_ii_loc is local index in BFvec1%vec(level)
						index_jj_loc=(index_j-BFvec1%vec(level)%idx_c)/BFvec1%vec(level)%inc_c+1


						if(allocated(BFvec1%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix))then
							mm=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,1)
							nn2=size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix,2)
							nvec2=size(BFvec1%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix,2)

							allocate(mat2(mm+2,nvec2))
							mat2=0
							mat2(1:2,:)=BFvec1%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix(1:2,:)
							call gemmf77('N','N',mm,nvec2,nn2, cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix, mm,BFvec1%vec(level)%blocks(index_ii_loc,index_jj_loc)%matrix(3,1),nn2+2,czero,mat2(3,1),mm+2)
							stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm,nvec2,nn2)

							if(allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then
								nvec1=size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,2)
								allocate(mat1(mm+2,nvec1))
								mat1=BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix
								deallocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)
							else
								nvec1=0
								allocate(mat1(mm+2,nvec1))
							endif

							allocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm+2,nvec1+nvec2))
							if(nvec1>0)BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(:,1:nvec1)=mat1
							if(nvec2>0)BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(:,1+nvec1:nvec1+nvec2)=mat2

							deallocate(mat1)
							deallocate(mat2)

						endif


					endif
				enddo
			enddo


			do nn=1,size(BFvec1%vec(level)%index,1)
				index_i_loc_s = BFvec1%vec(level)%index(nn,1)
				index_j_loc_s = BFvec1%vec(level)%index(nn,2)
				if(allocated(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index))deallocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index)
				if(allocated(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix))deallocate(BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix)
			enddo



			call BF_exchange_extraction(blocks,BFvec1%vec(level+1),stats,ptree,level,'R')



		    !!!!!! sort the intersection# here
			do nn=1,size(BFvec1%vec(level+1)%index,1)
				index_i_loc_s = BFvec1%vec(level+1)%index(nn,1)
				index_j_loc_s = BFvec1%vec(level+1)%index(nn,2)

				if(allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix))then


! if(ptree%MyID>=2)write(*,*)ptree%MyID,'fanibef',level+1,index_i_loc_s,index_j_loc_s,shape(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)

					gg1=0
					nvec1=0
					do jj = 1, size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,2)
					gg=NINT(dble(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(1,jj)))

					if(gg<gg1)then
						nvec1=jj-1
						exit
					endif
					gg1=gg
					enddo
					nvec2=size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,2)-nvec1

					mm = size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix,1)
					allocate(mat1(mm,nvec1))
					if(nvec1>0)mat1=BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(:,1:nvec1)
					allocate(mat2(mm,nvec2))
					if(nvec2>0)mat2=BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(:,1+nvec1:nvec1+nvec2)
					deallocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)

					!**** filter out the columns of mat1 and mat2 that are not specified by BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index
					allocate(mat(mm,nvec1+nvec2))
					idx1=0
					idx2=0
					idx=0
					gg1=0
					gg2=0
					do ii = 1, size(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index,1)
						gg=NINT(dble(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(ii,1)))
						do while(gg1<=gg)
							if(gg1==gg)then
								idx=idx+1
								mat(:,idx)=mat1(:,idx1)
							endif
							idx1=idx1+1
							if(idx1>nvec1)exit
							gg1=NINT(dble(mat1(1,idx1)))
						enddo
						do while(gg2<=gg)
							if(gg2==gg)then
								idx=idx+1
								mat(:,idx)=mat2(:,idx2)
							endif
							idx2=idx2+1
							if(idx2>nvec2)exit
							gg2=NINT(dble(mat2(1,idx2)))
						enddo
					enddo

					allocate(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix(mm,idx))
					BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix = mat(:,1:idx)

					deallocate(mat)
					deallocate(mat1)
					deallocate(mat2)
				endif

				! if(ptree%MyID>=2)write(*,*)ptree%MyID,'fani',level+1,index_i_loc_s,index_j_loc_s,shape(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec1%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)

			enddo
		endif
	enddo


	deallocate(group_ms)
	deallocate(group_ns)

	do level=level_half+1,level_butterfly+2
		if(allocated(BFvec1%vec(level)%index))deallocate(BFvec1%vec(level)%index)
	enddo
	deallocate(BFvec1%vec)

	do level=0, level_half+1
		if(allocated(BFvec%vec(level)%index))deallocate(BFvec%vec(level)%index)
	enddo
	deallocate(BFvec%vec)

	n5 = OMP_get_wtime()
	time_tmp = time_tmp + n5 - n1

end subroutine BF_block_extraction



!*** Find the group index of point idx at the (group%level+level) level
integer function findgroup(idx,msh,level,group)
    use BPACK_DEFS
    implicit none
	integer idx,level,ll,group,group1
	type(mesh)::msh
	group1=group
	if(idx<msh%basis_group(group1)%head .or. idx>msh%basis_group(group1)%tail)then
		findgroup=-1
	else
		do ll=1,level
			if(idx<=msh%basis_group(2*group1)%tail)then
				group1=group1*2
			else
				group1=group1*2+1
			endif
		enddo
		findgroup=group1
	endif

end function findgroup



!*** Find the process group index of point idx in a group
integer function findpggroup(idx,msh,ptree,group,pgno)
    use BPACK_DEFS
    implicit none
	integer idx,ll,group,group1
	type(mesh)::msh
	type(proctree)::ptree
	integer level_p,pgno,pgno1

	level_p = ptree%nlevel-GetTreelevel(pgno)

	group1=group
	pgno1=pgno
	if(idx<msh%basis_group(group1)%head .or. idx>msh%basis_group(group1)%tail)then
		findpggroup=-1
	else
		do ll=1,level_p
			if(idx<=msh%basis_group(2*group1)%tail)then
				group1=group1*2
				pgno1=pgno1*2
			else
				group1=group1*2+1
				pgno1=pgno1*2+1
			endif
		enddo
		findpggroup=pgno1
	endif

end function findpggroup


subroutine BF_value(mi,nj,blocks,value)

    use BPACK_DEFS
    implicit none

    integer mm, nn, mi, nj, groupm_start, groupn_start, level_butterfly, flag
    integer i, j, ii, jj, rank, group_m, group_n, header_mm, header_nn, k, kk
    integer group, level, mii, njj, rank1, rank2, index_ij, level_blocks, flag1
    DT ctemp, value

    type(matrixblock) :: blocks
    type(vectorset),allocatable:: vectors_set(:)
    integer,allocatable :: group_index_mm(:), group_index_nn(:)





    level_butterfly=blocks%level_butterfly




    allocate (group_index_mm(0:level_butterfly),group_index_nn(0:level_butterfly))

    flag=0; i=0; k=0
    do while (flag==0)
        i=i+1
        if (size(blocks%ButterflyU%blocks(i)%matrix,1)+k>=mi) then
            flag=1
        endif
        k=k+size(blocks%ButterflyU%blocks(i)%matrix,1)
    enddo
    group_index_mm(0)=i
    mii=mi-k+size(blocks%ButterflyU%blocks(i)%matrix,1)

    flag=0; j=0; k=0
    do while (flag==0)
        j=j+1
        if (size(blocks%ButterflyV%blocks(j)%matrix,1)+k>=nj) then
            flag=1
        endif
        k=k+size(blocks%ButterflyV%blocks(j)%matrix,1)
    enddo
    group_index_nn(0)=j
    njj=nj-k+size(blocks%ButterflyV%blocks(j)%matrix,1)

    if (level_butterfly>0) then
        group_index_mm(1)=group_index_mm(0)
        group_index_nn(1)= group_index_nn(0)
        do level=1, level_butterfly-1
            group_index_mm(level+1)=int((group_index_mm(level)+1)/2)
            group_index_nn(level+1)=int((group_index_nn(level)+1)/2)
        enddo
    endif

!     if (group_index_mm(0)/=group_m .or. group_index_nn(0)/=group_n) then
!         write (*,*) 'BF_value_func error1!'
!         pause
!         continue
!     endif

!     do level=0, level_butterfly
!         group_index_mm(level)=group_index_mm(level)-group_m*2**level+1
!         group_index_nn(level)=group_index_nn(level)-group_n*2**level+1
!     enddo

    allocate (vectors_set(0:level_butterfly))
    do level=0, level_butterfly
        if (level==0) then
            rank=size(blocks%ButterflyV%blocks(group_index_nn(0))%matrix,2)
            allocate (vectors_set(level)%vector(rank))
            !!$omp parallel do default(shared) private(i)
            do i=1, rank
                vectors_set(level)%vector(i)=blocks%ButterflyV%blocks(group_index_nn(0))%matrix(njj,i)
            enddo
            !!$omp end parallel do
			! write(*,*)'seq: ',level, abs(sum(vectors_set(level)%vector))
        else
            rank1=size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,2)
            rank2=size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,1)
            allocate (vectors_set(level)%vector(rank2))
            !!$omp parallel do default(shared) private(i,j,ctemp)
            do i=1, rank2
                ctemp=0
                do j=1, rank1
                    ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix(i,j)*vectors_set(level-1)%vector(j)
                enddo
                vectors_set(level)%vector(i)=ctemp
            enddo
            !!$omp end parallel do
            deallocate (vectors_set(level-1)%vector)
			! write(*,*)'seq: ',level, abs(sum(vectors_set(level)%vector))
        endif
        if (level==level_butterfly) then
            rank=size(vectors_set(level)%vector,1)
            ctemp=0
            !!$omp parallel do default(shared) private(i) reduction(+:ctemp)
			! write(*,*)'seq: ', level, abs(sum(blocks%ButterflyU%blocks(group_index_mm(0))%matrix(mii,1:rank))),mii
            do i=1, rank
                ctemp=ctemp+blocks%ButterflyU%blocks(group_index_mm(0))%matrix(mii,i)*vectors_set(level)%vector(i)
            enddo
            !!$omp end parallel do
            value=ctemp
            deallocate (vectors_set(level)%vector)
        endif
    enddo
    deallocate (vectors_set)

    return

end subroutine BF_value


subroutine BF_get_rank(block_i,ptree)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly
integer::ierr
type(proctree)::ptree

block_i%rankmin = 100000
block_i%rankmax = -100000

if(IOwnPgrp(ptree,block_i%pgno))then

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly



do level=0, level_butterfly+1
	if(level==0)then
		do jj=1,block_i%ButterflyV%nblk_loc
			nn=size(block_i%ButterflyV%blocks(jj)%matrix,1)
			rank=size(block_i%ButterflyV%blocks(jj)%matrix,2)
			block_i%rankmin = min(block_i%rankmin,rank)
			block_i%rankmax = max(block_i%rankmax,rank)
		enddo
	elseif(level==level_butterfly+1)then
		do jj=1,block_i%ButterflyU%nblk_loc
			mm=size(block_i%ButterflyU%blocks(jj)%matrix,1)
			rank=size(block_i%ButterflyU%blocks(jj)%matrix,2)
			block_i%rankmin = min(block_i%rankmin,rank)
			block_i%rankmax = max(block_i%rankmax,rank)
		enddo
	else
		do ii=1, block_i%ButterflyKerl(level)%nr
			do jj=1, block_i%ButterflyKerl(level)%nc
				nn=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,2)
				rank=size(block_i%ButterflyKerl(level)%blocks(ii,jj)%matrix,1)
				block_i%rankmin = min(block_i%rankmin,rank)
				block_i%rankmax = max(block_i%rankmax,rank)
			enddo
		enddo
	endif
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%rankmax,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(block_i%pgno)%Comm,ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE,block_i%rankmin,1,MPI_INTEGER,MPI_MIN,ptree%pgrp(block_i%pgno)%Comm,ierr)

endif


end subroutine BF_get_rank







subroutine BF_sym2asym(blocks)

    use BPACK_DEFS
	use MISC_Utilities


    implicit none

    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    DT ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: dimension_n,num_row,num_col,mn_min

	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:)

    real(kind=8), allocatable :: Singular(:)
    DT, allocatable :: UU(:,:),VV(:,:)


	if(allocated(blocks%ButterflyMiddle))then



        group_m=blocks%row_group ! Note: row_group and col_group interchanged here
        group_n=blocks%col_group
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
	    levelm = floor_safe(dble(level_butterfly)/2d0)

		call assert(level_butterfly>=2,'level_butterfly not correct')

		level = levelm
		num_groupm=blocks%ButterflyKerl(level)%num_row
		num_groupn=blocks%ButterflyKerl(level)%num_col


		! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
		do ij=1,num_groupm*(num_groupn/2)
			i = (ij-1)/(num_groupn/2)+1
			j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
			index_i=int((i+1)/2)
			index_j=int((j+1)/2)

			nn1=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
			nn2=size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
			mm=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)

			allocate(matrixtemp(mm,nn1))
			matrixtemp = blocks%ButterflyKerl(level)%blocks(i,j)%matrix
			! call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm,nn1,mm)
			call gemmf90(blocks%ButterflyMiddle(i,index_j)%matrix,mm,matrixtemp,mm,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm,'N','N',mm,nn1,mm,cone,czero)
			deallocate(matrixtemp)

			allocate(matrixtemp(mm,nn2))
			matrixtemp = blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix
			! call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,nn2,mm)
			call gemmf90(blocks%ButterflyMiddle(i,index_j)%matrix,mm,matrixtemp,mm,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,'N','N',mm,nn2,mm,cone,czero)
			deallocate(matrixtemp)

			deallocate(blocks%ButterflyMiddle(i,index_j)%matrix)
		enddo
		! !$omp end parallel do

		deallocate(blocks%ButterflyMiddle)

		do level=0, levelm-1
			if(level==0)then

				iijj=0
				do j=1, num_blocks
					iijj = iijj + 1
					dimension_n=size(blocks%ButterflyV%blocks(j)%matrix,1)
					rank = size(blocks%ButterflyV%blocks(j)%matrix,2)
					mn_min = min(dimension_n,rank)

					allocate(matrixtemp(rank,dimension_n))
					allocate(UU(rank,mn_min))
					allocate(VV(mn_min,dimension_n))
					allocate(Singular(mn_min))

					call copymatT(blocks%ButterflyV%blocks(j)%matrix,matrixtemp,dimension_n,rank)

					call gesvd_robust(matrixtemp,Singular,UU,VV,rank,dimension_n,mn_min)
					do ii=1,mn_min
						UU(:,ii) = UU(:,ii)*Singular(ii)
					end do

					deallocate(blocks%ButterflyV%blocks(j)%matrix)
					allocate(blocks%ButterflyV%blocks(j)%matrix(dimension_n,mn_min))
					call copymatT(VV,blocks%ButterflyV%blocks(j)%matrix,mn_min,dimension_n)

					index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
					index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))
					mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
					allocate(matrixtemp1(mm1,mn_min))
					matrixtemp1=0
					! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
					call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,mm1,UU,rank,matrixtemp1,mm1,'N','N',mm1,mn_min,rank,cone,czero)

					deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
					allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
					blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
					deallocate(matrixtemp1)

					mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
					allocate(matrixtemp1(mm2,mn_min))
					! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
					call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,mm2,UU,rank,matrixtemp1,mm2,'N','N',mm2,mn_min,rank,cone,czero)
					deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
					allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
					blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
					deallocate(matrixtemp1)

					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)

				enddo
			else
				num_row=blocks%ButterflyKerl(level)%num_row
				num_col=blocks%ButterflyKerl(level)%num_col

				iijj=0
				do i=1,	num_row
					do j =1, num_col, 2
						iijj = iijj + 1
						rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						nn1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2 = size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mn_min = min(nn1+nn2,rank)

						allocate(matrixtemp(rank,nn1+nn2))
						allocate(UU(rank,mn_min))
						allocate(VV(mn_min,nn1+nn2))
						allocate(Singular(mn_min))

						! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
						matrixtemp(1:rank,1:nn1) = blocks%ButterflyKerl(level)%blocks(i,j)%matrix
						! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
						matrixtemp(1:rank,1+nn1:nn2+nn1) = blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix

						call gesvd_robust(matrixtemp,Singular,UU,VV,rank,nn1+nn2,mn_min)
						do ii=1,mn_min
							UU(:,ii) = UU(:,ii)*Singular(ii)
						end do

						deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
						allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mn_min,nn1))
						! call copymatN(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
						blocks%ButterflyKerl(level)%blocks(i,j)%matrix = VV(1:mn_min,1:nn1)
						deallocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix)
						allocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(mn_min,nn2))
						! call copymatN(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
						blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix = VV(1:mn_min,1+nn1:nn2+nn1)


						index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
						index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))


						mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
						call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,mm1,UU,rank,matrixtemp1,mm1,'N','N',mm1,mn_min,rank,cone,czero)
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)

						mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
						allocate(matrixtemp1(mm2,mn_min))
						! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
						call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,mm2,UU,rank,matrixtemp1,mm2,'N','N',mm2,mn_min,rank,cone,czero)
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)

						deallocate(matrixtemp)
						deallocate(UU)
						deallocate(VV)
						deallocate(Singular)

					end do
				end do
			end if
		end do

	end if

end subroutine BF_sym2asym




subroutine BF_MoveSingulartoLeft(blocks)

    use BPACK_DEFS
	use MISC_Utilities


    implicit none

    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    DT ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: dimension_n,dimension_m,num_row,num_col,mn_min

	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:)

    real(kind=8), allocatable :: Singular(:)
    DT, allocatable :: UU(:,:),VV(:,:)


	group_m=blocks%row_group ! Note: row_group and col_group interchanged here
	group_n=blocks%col_group
	level_butterfly=blocks%level_butterfly
	num_blocks=2**level_butterfly

	do level=0, level_butterfly
		if(level==0)then
			iijj=0
			do j=1, num_blocks
				iijj = iijj + 1
				dimension_n=size(blocks%ButterflyV%blocks(j)%matrix,1)
				rank = size(blocks%ButterflyV%blocks(j)%matrix,2)
				mn_min = min(dimension_n,rank)

				allocate(matrixtemp(rank,dimension_n))
				allocate(UU(rank,mn_min))
				allocate(VV(mn_min,dimension_n))
				allocate(Singular(mn_min))

				call copymatT(blocks%ButterflyV%blocks(j)%matrix,matrixtemp,dimension_n,rank)
				call assert(.not. isnan(fnorm(matrixtemp,rank,dimension_n)),'matrixtemp NAN at 3')

				call gesvd_robust(matrixtemp,Singular,UU,VV,rank,dimension_n,mn_min)
				call assert(.not. isnan(sum(Singular)),'Singular NAN at 3')

				do ii=1,mn_min
					UU(:,ii) = UU(:,ii)*Singular(ii)
				end do


				deallocate(blocks%ButterflyV%blocks(j)%matrix)
				allocate(blocks%ButterflyV%blocks(j)%matrix(dimension_n,mn_min))
				call copymatT(VV,blocks%ButterflyV%blocks(j)%matrix,mn_min,dimension_n)


				index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
				index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))

				mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
				allocate(matrixtemp1(mm1,mn_min))
				! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
				call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,mm1,UU,rank,matrixtemp1,mm1,'N','N',mm1,mn_min,rank,cone,czero)

				deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
				allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
				blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
				deallocate(matrixtemp1)

				mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
				allocate(matrixtemp1(mm2,mn_min))
				! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
				call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,mm2,UU,rank,matrixtemp1,mm2,'N','N',mm2,mn_min,rank,cone,czero)
				deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
				allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
				blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
				deallocate(matrixtemp1)

				deallocate(matrixtemp)
				deallocate(UU)
				deallocate(VV)
				deallocate(Singular)

			enddo
		else
			num_row=blocks%ButterflyKerl(level)%num_row
			num_col=blocks%ButterflyKerl(level)%num_col

			iijj=0
			do i=1,	num_row
				do j =1, num_col, 2
					iijj = iijj + 1
					rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
					nn1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
					nn2 = size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
					mn_min = min(nn1+nn2,rank)

					allocate(matrixtemp(rank,nn1+nn2))
					allocate(UU(rank,mn_min))
					allocate(VV(mn_min,nn1+nn2))
					allocate(Singular(mn_min))

					! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
					matrixtemp(1:rank,1:nn1) = blocks%ButterflyKerl(level)%blocks(i,j)%matrix
					! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
					matrixtemp(1:rank,1+nn1:nn2+nn1) = blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix
					call assert(.not. isnan(fnorm(matrixtemp,rank,nn1+nn2)),'matrixtemp NAN at 4')
					call gesvd_robust(matrixtemp,Singular,UU,VV,rank,nn1+nn2,mn_min)
					call assert(.not. isnan(sum(Singular)),'Singular NAN at 4')

					do ii=1,mn_min
						UU(:,ii) = UU(:,ii)*Singular(ii)
					end do

					deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mn_min,nn1))
					! call copymatN(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
					blocks%ButterflyKerl(level)%blocks(i,j)%matrix = VV(1:mn_min,1:nn1)
					deallocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(mn_min,nn2))
					! call copymatN(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
					blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix = VV(1:mn_min,1+nn1:nn2+nn1)

					if(level/=level_butterfly)then
						index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
						index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))

						mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)

						call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,mm1,UU,rank,matrixtemp1,mm1,'N','N',mm1,mn_min,rank,cone,czero)

						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)

						mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
						allocate(matrixtemp1(mm2,mn_min))
						! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
						call gemmf90(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,mm2,UU,rank,matrixtemp1,mm2,'N','N',mm2,mn_min,rank,cone,czero)

						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)
					else
						mm1 = size(blocks%ButterflyU%blocks(i)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						! call gemm_omp(blocks%ButterflyU%blocks(i)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
						call gemmf90(blocks%ButterflyU%blocks(i)%matrix,mm1,UU,rank,matrixtemp1,mm1,'N','N',mm1,mn_min,rank,cone,czero)
						deallocate(blocks%ButterflyU%blocks(i)%matrix)
						allocate(blocks%ButterflyU%blocks(i)%matrix(mm1,mn_min))
						blocks%ButterflyU%blocks(i)%matrix = matrixtemp1
						deallocate(matrixtemp1)
					end if

					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)

				end do
			end do
		end if
	end do


end subroutine BF_MoveSingulartoLeft





subroutine BF_MoveSingulartoRight(blocks)

    use BPACK_DEFS
	use MISC_Utilities


    implicit none

    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    DT ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: dimension_n,dimension_m,num_row,num_col,mn_min

	DT,allocatable::matrixtemp(:,:),matrixtemp1(:,:)

    real(kind=8), allocatable :: Singular(:)
    DT, allocatable :: UU(:,:),VV(:,:)


	group_m=blocks%row_group ! Note: row_group and col_group interchanged here
	group_n=blocks%col_group
	level_butterfly=blocks%level_butterfly
	num_blocks=2**level_butterfly

	do level=level_butterfly+1, 1,-1
		if(level==level_butterfly+1)then
			iijj=0
			do i=1, num_blocks
				iijj = iijj + 1
				dimension_m=size(blocks%ButterflyU%blocks(i)%matrix,1)
				rank = size(blocks%ButterflyU%blocks(i)%matrix,2)
				mn_min = min(dimension_m,rank)

				allocate(matrixtemp(dimension_m,rank))
				allocate(UU(dimension_m,mn_min))
				allocate(VV(mn_min,rank))
				allocate(Singular(mn_min))

				! call copymatN(blocks%ButterflyU%blocks(i)%matrix,matrixtemp,dimension_m,rank)
				matrixtemp = blocks%ButterflyU%blocks(i)%matrix
				call assert(.not. isnan(fnorm(matrixtemp,dimension_m,rank)),'matrixtemp NAN at 1')

				call gesvd_robust(matrixtemp,Singular,UU,VV,dimension_m,rank,mn_min)
				call assert(.not. isnan(sum(Singular)),'Singular NAN at 1')

				do ii=1,mn_min
					VV(ii,:) = VV(ii,:)*Singular(ii)
				end do

				deallocate(blocks%ButterflyU%blocks(i)%matrix)
				allocate(blocks%ButterflyU%blocks(i)%matrix(dimension_m,mn_min))
				! call copymatN(UU,blocks%ButterflyU%blocks(i)%matrix,dimension_m,mn_min)
				blocks%ButterflyU%blocks(i)%matrix = UU

				index_i = mod(iijj-1,blocks%ButterflyKerl(level-1)%num_row)+1
				index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level-1)%num_row))

				nn1=size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,2)
				allocate(matrixtemp1(mn_min,nn1))
				! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,nn1,rank)
				call gemmf90(VV,mn_min,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,rank,matrixtemp1,mn_min,'N','N',mn_min,nn1,rank,cone,czero)

				deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix)
				allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix(mn_min,nn1))
				blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix = matrixtemp1
				deallocate(matrixtemp1)

				nn2=size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,2)
				allocate(matrixtemp1(mn_min,nn2))
				! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,nn2,rank)
				call gemmf90(VV,mn_min,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,rank,matrixtemp1,mn_min,'N','N',mn_min,nn2,rank,cone,czero)
				deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix)
				allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix(mn_min,nn2))
				blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix = matrixtemp1
				deallocate(matrixtemp1)

				deallocate(matrixtemp)
				deallocate(UU)
				deallocate(VV)
				deallocate(Singular)

			enddo
		else
			num_row=blocks%ButterflyKerl(level)%num_row
			num_col=blocks%ButterflyKerl(level)%num_col

			iijj=0
			do j=1,	num_col
				do i =1, num_row, 2
					iijj = iijj + 1
					rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)

					mm1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
					mm2 = size(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,1)
					mn_min = min(mm1+mm2,rank)

					allocate(matrixtemp(mm1+mm2,rank))
					allocate(UU(mm1+mm2,mn_min))
					allocate(VV(mn_min,rank))
					allocate(Singular(mn_min))

					! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:mm1,1:rank),mm1,rank)
					matrixtemp(1:mm1,1:rank) = blocks%ButterflyKerl(level)%blocks(i,j)%matrix
					! call copymatN(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,matrixtemp(1+mm1:mm2+mm1,1:rank),mm2,rank)
					matrixtemp(1+mm1:mm2+mm1,1:rank) = blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix
					call assert(.not. isnan(fnorm(matrixtemp,mm1+mm2,rank)),'matrixtemp NAN at 2')

					call gesvd_robust(matrixtemp,Singular,UU,VV,mm1+mm2,rank,mn_min)
					! if(isnan(sum(Singular)).and. mm1+mm2<rank)then
						! write(*,*)mm1+mm2,rank,mm1+mm2>=rank,'rank too large?'
					! end if

					! call assert(.not. isnan(sum(Singular)),'Singular NAN at 2')
					if(isnan(sum(Singular)))then
						write(*,*)'Singular NAN at 2',mm1+mm2,rank
						do ii=1,mm1+mm2
							do jj=1,rank
								write(777,*)dble(matrixtemp(ii,jj)),aimag(cmplx(matrixtemp(ii,jj),kind=8)),abs(matrixtemp(ii,jj))
							end do
						end do
						stop
					end if



					do ii=1,mn_min
						VV(ii,:) = VV(ii,:)*Singular(ii)
					end do

					deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mm1,mn_min))
					! call copymatN(UU(1:mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm1,mn_min)
					blocks%ButterflyKerl(level)%blocks(i,j)%matrix = UU(1:mm1,1:mn_min)
					deallocate(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix(mm2,mn_min))
					! call copymatN(UU(1+mm1:mm2+mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,mm2,mn_min)
					blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix = UU(1+mm1:mm2+mm1,1:mn_min)

					if(level/=1)then
						index_i = mod(iijj-1,blocks%ButterflyKerl(level-1)%num_row)+1
						index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level-1)%num_row))
						nn1 = size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,2)

						allocate(matrixtemp1(mn_min,nn1))
						! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,nn1,rank)
						call gemmf90(VV,mn_min,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,rank,matrixtemp1,mn_min,'N','N',mn_min,nn1,rank,cone,czero)

						deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix)
						allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix(mn_min,nn1))
						blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix = matrixtemp1
						deallocate(matrixtemp1)

						nn2 = size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,2)
						allocate(matrixtemp1(mn_min,nn2))
						! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,nn2,rank)
						call gemmf90(VV,mn_min,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,rank,matrixtemp1,mn_min,'N','N',mn_min,nn2,rank,cone,czero)


						deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix)
						allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix(mn_min,nn2))
						blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix = matrixtemp1
						deallocate(matrixtemp1)
					else
						nn1 = size(blocks%ButterflyV%blocks(j)%matrix,1)
						allocate(matrixtemp1(nn1,mn_min))
						! call gemmNT_omp(blocks%ButterflyV%blocks(j)%matrix,VV,matrixtemp1,nn1,mn_min,rank)
						call gemmf90(blocks%ButterflyV%blocks(j)%matrix,nn1, VV,mn_min, matrixtemp1,nn1, 'N','T',nn1,mn_min,rank,cone,czero)
						deallocate(blocks%ButterflyV%blocks(j)%matrix)
						allocate(blocks%ButterflyV%blocks(j)%matrix(nn1,mn_min))
						blocks%ButterflyV%blocks(j)%matrix = matrixtemp1
						deallocate(matrixtemp1)
					end if

					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)

				end do
			end do
		end if
	end do


end subroutine BF_MoveSingulartoRight





subroutine BF_Init_blocks(level_butterfly,groupm,groupn,pgno,block_rand,msh,ptree)

    use BPACK_DEFS
    implicit none

    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max,dimension_m, dimension_n, blocks, groupm, groupm_start,groupn_start,groupn,index_j,index_i
    real(kind=8) a,b,c,d
    DT ctemp
	type(matrixblock)::block,block_rand
	DT, allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real(kind=8), allocatable:: Singular(:)
	type(mesh)::msh
	! type(Hoption)::option
	type(proctree)::ptree
	integer level_final,level_half
	integer idx_r,inc_r,nr,idx_c,inc_c,nc
	integer pgno

    block_rand%level_butterfly=level_butterfly
    num_blocks=2**level_butterfly

	! level_half = BF_Switchlevel(level_butterfly,option)
	level_half = floor_safe(dble(level_butterfly)/2d0) ! from outer to inner
	block_rand%level_half=level_half

    block_rand%style=2
    block_rand%row_group=groupm
    block_rand%col_group=groupn

	block_rand%M = msh%basis_group(block_rand%row_group)%tail - msh%basis_group(block_rand%row_group)%head + 1
	block_rand%N = msh%basis_group(block_rand%col_group)%tail - msh%basis_group(block_rand%col_group)%head + 1
	block_rand%headm = msh%basis_group(block_rand%row_group)%head
	block_rand%headn = msh%basis_group(block_rand%col_group)%head


	block_rand%pgno = pgno

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

    return

end subroutine BF_Init_blocks


recursive subroutine Hmat_block_copy(trans,block2,block1,memory)

    use BPACK_DEFS
    implicit none

    integer blocks, flag_recv, count1, count2, recv_count, mm, nn, length
    integer i, ii, j, jj, style, send_ID, group_m, group_n, indices, requests
    character chara

    type(matrixblock), pointer :: block1, block2, blocks_son1, blocks_son2
	character::trans
	real(kind=8),optional::memory
	real(kind=8)::memory_tmp

	block2%style = block1%style

	block2%level=block1%level
	block2%row_group=block1%row_group
	block2%col_group=block1%col_group
	block2%level_butterfly=0
	group_m = block2%row_group
	group_n = block2%col_group
	block2%pgno = block1%pgno
	block2%M = block1%M
	block2%N = block1%N
	block2%headm = block1%headm
	block2%headn = block1%headn

	if(associated(block1%N_p))then
		if(associated(block2%N_p))deallocate(block2%N_p)
		allocate(block2%N_p(size(block1%N_p,1),2))
		block2%N_p = block1%N_p
	endif
	if(associated(block1%M_p))then
		if(associated(block2%M_p))deallocate(block2%M_p)
		allocate(block2%M_p(size(block1%M_p,1),2))
		block2%M_p = block1%M_p
	endif


    style=block2%style
    if (style==4) then
        allocate(block2%sons(2,2))
        do j=1,2
            do i=1,2
                block2%sons(i,j)%father=>block2
            enddo
        enddo

		blocks_son1=>block1%sons(1,1)
		blocks_son2=>block2%sons(1,1)
		call Hmat_block_copy(trans,blocks_son2,blocks_son1,memory)
		blocks_son1=>block1%sons(2,1)
		blocks_son2=>block2%sons(2,1)
		call Hmat_block_copy(trans,blocks_son2,blocks_son1,memory)
		 blocks_son1=>block1%sons(1,2)
		blocks_son2=>block2%sons(1,2)
		call Hmat_block_copy(trans,blocks_son2,blocks_son1,memory)
		 blocks_son1=>block1%sons(2,2)
		blocks_son2=>block2%sons(2,2)
		call Hmat_block_copy(trans,blocks_son2,blocks_son1,memory)

    else
		call BF_copy(trans,block1,block2,memory_tmp)
		if(present(memory))memory = memory + memory_tmp
    endif

    return

end subroutine Hmat_block_copy


recursive subroutine Hmat_block_delete(blocks)


    implicit none

    integer level_actual, num_col, num_row
    integer i, j, mm, nn, rank, num_blocks, level, level_butterfly
    real*8 memory_butterfly, rtemp
    type(matrixblock) :: blocks
    type(matrixblock), pointer :: blocks_son

    if (blocks%style==4) then

		blocks_son=>blocks%sons(1,1)
		call Hmat_block_delete(blocks_son)
		blocks_son=>blocks%sons(2,1)
		call Hmat_block_delete(blocks_son)
		blocks_son=>blocks%sons(1,2)
		call Hmat_block_delete(blocks_son)
		blocks_son=>blocks%sons(2,2)
		call Hmat_block_delete(blocks_son)

        deallocate (blocks%sons)

    else
        call BF_delete(blocks,1)
    endif

    return

end subroutine Hmat_block_delete




recursive subroutine Hmat_block_ComputeMemory(blocks,memory)


    implicit none

    integer level_actual, num_col, num_row
    integer i, j, mm, nn, rank, num_blocks, level, level_butterfly
    real*8 memory_butterfly, rtemp,memory
    type(matrixblock) :: blocks
    type(matrixblock), pointer :: blocks_son

    if (blocks%style==4) then

		blocks_son=>blocks%sons(1,1)
		call Hmat_block_ComputeMemory(blocks_son,memory)
		blocks_son=>blocks%sons(2,1)
		call Hmat_block_ComputeMemory(blocks_son,memory)
		blocks_son=>blocks%sons(1,2)
		call Hmat_block_ComputeMemory(blocks_son,memory)
		blocks_son=>blocks%sons(2,2)
		call Hmat_block_ComputeMemory(blocks_son,memory)
    else
        call BF_ComputeMemory(blocks,rtemp)
		memory = memory + rtemp
    endif

    return

end subroutine Hmat_block_ComputeMemory




recursive subroutine Hmat_Lsolve(blocks_l,trans,idx_start,nvec,Vinout,ptree,stats)
    implicit none

    ! integer vectors_y
    integer style(3)
    integer i, j, k, ii
    integer mm, nn, nvec,idxs_m, idx_start ! idx_start means the global indice of the first element of Vinout
    integer head, tail
    DT ctemp
    DT:: Vinout(:,:)
    type(matrixblock) :: blocks_l !!!! modified by Yang Liu. passing pointer is dangerous, blocks_u row/row_group becomes different once in this subroutine
	character trans ! 'N' means multiple L^-1 from left, 'T' means multiple L^-1 from right
	type(proctree)::ptree
	type(Hstat)::stats

    if (blocks_l%style==4) then
		if(trans=='N')then
			call Hmat_Lsolve(blocks_l%sons(1,1),trans,idx_start,nvec,Vinout,ptree,stats)
			call Hmat_block_MVP_dat(blocks_l%sons(2,1),trans,idx_start,idx_start,nvec,Vinout,Vinout,-cone,ptree,stats)
			call Hmat_Lsolve(blocks_l%sons(2,2),trans,idx_start,nvec,Vinout,ptree,stats)
		else
			call Hmat_Lsolve(blocks_l%sons(2,2),trans,idx_start,nvec,Vinout,ptree,stats)
			call Hmat_block_MVP_dat(blocks_l%sons(2,1),trans,idx_start,idx_start,nvec,Vinout,Vinout,-cone,ptree,stats)
			call Hmat_Lsolve(blocks_l%sons(1,1),trans,idx_start,nvec,Vinout,ptree,stats)
		end if
    else
		mm = blocks_l%M
		idxs_m = blocks_l%headm - idx_start + 1

		if(trans=='N')then
			do i=1, mm
				ii=blocks_l%ipiv(i)
				if (ii/=i) then
					!$omp parallel do default(shared) private(j,ctemp)
					do j=1, nvec
						ctemp=Vinout(idxs_m+i-1,j)
						Vinout(idxs_m+i-1,j)=Vinout(idxs_m+ii-1,j)
						Vinout(idxs_m+ii-1,j)=ctemp
					enddo
					!$omp end parallel do
				endif
			enddo
		endif
		! write(*,*)blocks_l%level,blocks_l%pgno,ptree%MyID,blocks_l%headm,mm,idx_start,'daha'
		call trsmf90(blocks_l%fullmat,Vinout(idxs_m:idxs_m+mm-1,1:nvec),'L','L',trans,'U',mm,nvec)
		if(trans/='N')then
			do i=mm,1,-1
				ii=blocks_l%ipiv(i)
				if (ii/=i) then
					!$omp parallel do default(shared) private(j,ctemp)
					do j=1, nvec
						ctemp=Vinout(idxs_m+i-1,j)
						Vinout(idxs_m+i-1,j)=Vinout(idxs_m+ii-1,j)
						Vinout(idxs_m+ii-1,j)=ctemp
					enddo
					!$omp end parallel do
				endif
			enddo
		end if
    endif

    return

end subroutine Hmat_Lsolve


recursive subroutine Hmat_Usolve(blocks_u,trans,idx_start,nvec,Vinout,ptree,stats)
    implicit none


	type(proctree)::ptree
	type(Hstat)::stats

	integer vectors_x, vectors_y
    integer style(3), mark
    integer i, j, k,ii
    integer mm, nn, nvec
    integer head, tail
    DT Vinout(:,:)
    type(matrixblock) :: blocks_u, blocks !!!! modified by Yang Liu. passing pointer is dangerous, blocks_u row/row_group becomes different once in this subroutine
	character trans
    integer idx_start,idxs_m

    mark=0
    if (blocks_u%style==4) then
		if(trans=='N')then
			call Hmat_Usolve(blocks_u%sons(2,2),trans,idx_start,nvec,Vinout,ptree,stats)
			call Hmat_block_MVP_dat(blocks_u%sons(1,2),trans,idx_start,idx_start,nvec,Vinout,Vinout,-cone,ptree,stats)
			call Hmat_Usolve(blocks_u%sons(1,1),trans,idx_start,nvec,Vinout,ptree,stats)
		else
			call Hmat_Usolve(blocks_u%sons(1,1),trans,idx_start,nvec,Vinout,ptree,stats)
			call Hmat_block_MVP_dat(blocks_u%sons(1,2),trans,idx_start,idx_start,nvec,Vinout,Vinout,-cone,ptree,stats)
			call Hmat_Usolve(blocks_u%sons(2,2),trans,idx_start,nvec,Vinout,ptree,stats)
		end if

    else
		mm = blocks_u%M
		idxs_m = blocks_u%headm - idx_start + 1
		call trsmf90(blocks_u%fullmat,Vinout(idxs_m:idxs_m+mm-1,1:nvec),'L','U',trans,'N',mm,nvec)
    endif

    return

end subroutine Hmat_Usolve

recursive subroutine Hmat_block_MVP_dat(blocks,trans,idx_start_m,idx_start_n,Nrnd,Vin,Vout,a,ptree,stats)

    implicit none
	integer idx_start_m,idx_start_n
    integer Nrnd
    integer mm, nn, idxs_m,idxs_n
    DT a
    character trans
	type(matrixblock)::blocks
	type(matrixblock),pointer::blocks_son
	integer:: style
	DT,allocatable::Vintmp(:,:),Vouttmp(:,:)
	DT::Vin(:,:),Vout(:,:)
	type(proctree)::ptree
	type(Hstat)::stats


	style = blocks%style
	mm = blocks%M
	idxs_m = blocks%headm - idx_start_m + 1
	nn = blocks%N
	idxs_n = blocks%headn - idx_start_n + 1


    if (style==4) then
		blocks_son=>blocks%sons(1,1)
		call Hmat_block_MVP_dat(blocks_son,trans,idx_start_m,idx_start_n,Nrnd,Vin,Vout,a,ptree,stats)
		blocks_son=>blocks%sons(1,2)
		call Hmat_block_MVP_dat(blocks_son,trans,idx_start_m,idx_start_n,Nrnd,Vin,Vout,a,ptree,stats)
		blocks_son=>blocks%sons(2,1)
		call Hmat_block_MVP_dat(blocks_son,trans,idx_start_m,idx_start_n,Nrnd,Vin,Vout,a,ptree,stats)
		blocks_son=>blocks%sons(2,2)
		call Hmat_block_MVP_dat(blocks_son,trans,idx_start_m,idx_start_n,Nrnd,Vin,Vout,a,ptree,stats)
    else
        if (style==1) then
            if (trans=='N') then
				allocate(Vintmp(nn,Nrnd))
				Vintmp = Vin(idxs_n:idxs_n+nn-1,1:Nrnd)
				allocate(Vouttmp(mm,Nrnd))
				Vouttmp = 0
				call gemmf90(blocks%fullmat,mm,Vintmp,nn,Vouttmp,mm,trans,'N',mm,Nrnd,nn,a,czero)
				Vout(idxs_m:idxs_m+mm-1,1:Nrnd) = Vout(idxs_m:idxs_m+mm-1,1:Nrnd)+Vouttmp
				deallocate(Vintmp)
				deallocate(Vouttmp)
			else
				allocate(Vintmp(mm,Nrnd))
				Vintmp = Vin(idxs_m:idxs_m+mm-1,1:Nrnd)
				allocate(Vouttmp(nn,Nrnd))
				Vouttmp = 0
				call gemmf90(blocks%fullmat,mm,Vintmp,mm,Vouttmp,nn,trans,'N',nn,Nrnd,mm,a,czero)
				Vout(idxs_n:idxs_n+nn-1,1:Nrnd) = Vout(idxs_n:idxs_n+nn-1,1:Nrnd)+Vouttmp
				deallocate(Vintmp)
				deallocate(Vouttmp)
		   endif
        else
			if (trans=='N') then
				call BF_block_MVP_dat(blocks,trans,mm,nn,Nrnd,Vin(idxs_n:idxs_n+nn-1,1:Nrnd),Vout(idxs_m:idxs_m+mm-1,1:Nrnd),a,cone,ptree,stats)
            else
				call BF_block_MVP_dat(blocks,trans,mm,nn,Nrnd,Vin(idxs_m:idxs_m+mm-1,1:Nrnd),Vout(idxs_n:idxs_n+nn-1,1:Nrnd),a,cone,ptree,stats)
            endif
        endif
    endif

end subroutine Hmat_block_MVP_dat



subroutine Full_block_MVP_dat(blocks,chara,M,N,random1,random2,a,b)
	use BPACK_DEFS


	use MISC_Utilities
    implicit none

    integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2
    DT ctemp, a, b
    character chara
	type(matrixblock)::blocks
	integer M,N
    DT :: random1(M,N), random2(M,N)
	DT:: al,be
	DT,allocatable :: random2tmp(:,:)

	allocate(random2tmp(M,N))



	al=1d0
	be=0d0

	num_vectors=size(random1,2)


	random2tmp = random2
	call assert(size(blocks%fullmat,1)==size(blocks%fullmat,2) ,'M not square')
	if(size(blocks%fullmat,1)/=M)write(*,*)M,N,shape(blocks%fullmat),blocks%row_group,blocks%col_group,'niao'
	call assert(size(blocks%fullmat,1)==M,'M not equal fullmat dim')

	if (chara=='N') then
        group_m=blocks%row_group  ! Note: row_group and col_group interchanged here
        group_n=blocks%col_group
		call assert(group_m==group_n,'fullmat not square')
        ! level_blocks=blocks%level
		! write(*,*)shape(blocks%fullmat),shape(random1),shape(random2),num_vectors

		! call gemm_omp(blocks%fullmat, random1, random2,M,N,M)
		call gemmf90(blocks%fullmat,M, random1,M, random2,M,'N','N',M,N,M,cone,czero)
    elseif (chara=='T') then
        group_m=blocks%row_group  ! Note: row_group and col_group interchanged here
        group_n=blocks%col_group
		call assert(group_m==group_n,'fullmat not square')
        ! level_blocks=blocks%level
		! call gemmTN_omp(blocks%fullmat, random1, random2,M,N,M)
		call gemmf90(blocks%fullmat,M, random1,M, random2,M, 'T','N',M,N,M,al,be)
	end if

	random2 = a*random2+ b*random2tmp
	! write(*,*)'wo cao ni ma'
	deallocate(random2tmp)
end subroutine Full_block_MVP_dat




! compute arrays M_p(1:P+1) and N_p(1:P+1) the holds the start and end column/row of each process sharing this block
subroutine ComputeParallelIndices(block,pgno,ptree,msh,flag)
implicit none
	type(matrixblock)::block
	integer pgno,level,level_p,level_butterfly,flag,nproc,num_blocks,proc,gg,ii,ii_new,Maxlevel
	type(proctree)::ptree
	integer,pointer::M_p(:,:),N_p(:,:)
	type(mesh)::msh

	if(flag==0)block%M_loc = 0
	if(flag==1)block%M_loc_db = 0
	if(flag==0)block%N_loc = 0
	if(flag==1)block%N_loc_db = 0

	Maxlevel = GetTreelevel(msh%Maxgroup)-1
	! write(*,*)msh%Maxgroup,GetTreelevel(msh%Maxgroup),Maxlevel-block%level,block%level,ptree%nlevel-GetTreelevel(pgno),pgno
	call assert(Maxlevel-block%level>=ptree%nlevel-GetTreelevel(pgno),'too many process sharing this group')

	! if(IOwnPgrp(ptree,pgno))then

		! level_butterfly = block%level_butterfly
		level_p = ptree%nlevel-GetTreelevel(pgno)
		nproc = ptree%pgrp(pgno)%nproc
		num_blocks = 2**level_p

		if(flag==0)then
			if(associated(block%M_p))deallocate(block%M_p)
			if(associated(block%N_p))deallocate(block%N_p)
			allocate(block%M_p(nproc,2))
			allocate(block%N_p(nproc,2))
			M_p => block%M_p
			N_p => block%N_p
		else
			if(associated(block%M_p_db))deallocate(block%M_p_db)
			if(associated(block%N_p_db))deallocate(block%N_p_db)
			allocate(block%M_p_db(nproc,2))
			allocate(block%N_p_db(nproc,2))
			M_p => block%M_p_db
			N_p => block%N_p_db
		endif

		M_p(:,1) = block%M+1
		N_p(:,1) = block%N+1
		M_p(:,2) = -block%M-1
		N_p(:,2) = -block%N-1

		do ii=1,num_blocks

			! if(flag==1)then  ! compute optimal renumbering of data pieces among the twice many processes
				! if(mod(ii,2)==1)then
					! ii_new=ceiling_safe(ii/2d0)
				! else
					! ii_new=ii/2+num_blocks/2
				! endif
			! else
				ii_new=ii
			! endif

			gg = block%row_group*2**level_p+ii_new-1
			proc = ptree%pgrp(pgno*2**level_p+ii-1)%head - ptree%pgrp(pgno)%head
			M_p(proc+1,1) = min(M_p(proc+1,1),msh%basis_group(gg)%head-msh%basis_group(block%row_group)%head+1)
			M_p(proc+1,2) = max(M_p(proc+1,2),msh%basis_group(gg)%tail-msh%basis_group(block%row_group)%head+1)
			gg = block%col_group*2**level_p+ii_new-1
			N_p(proc+1,1) = min(N_p(proc+1,1),msh%basis_group(gg)%head-msh%basis_group(block%col_group)%head+1)
			N_p(proc+1,2) = max(N_p(proc+1,2),msh%basis_group(gg)%tail-msh%basis_group(block%col_group)%head+1)
		enddo

		if(IOwnPgrp(ptree,pgno))then
			ii = ptree%myid-ptree%pgrp(pgno)%head+1
			if(flag==0)block%M_loc = M_p(ii,2)-M_p(ii,1)+1
			if(flag==1)block%M_loc_db = M_p(ii,2)-M_p(ii,1)+1
			if(flag==0)block%N_loc = N_p(ii,2)-N_p(ii,1)+1
			if(flag==1)block%N_loc_db = N_p(ii,2)-N_p(ii,1)+1
		endif
		! write(*,*)level_butterfly,level_p,block%M_loc,block%N_loc,'nima',M_p,N_p,block%M,block%N,block%row_group,block%col_group
	! endif
end subroutine ComputeParallelIndices


function node_score_block_ptr_row(this) result(score)
	implicit none
	type(nod)::this
	real(kind=8)::score
	class(*),pointer::ptr

	select TYPE(ptr=>this%item)
		type is (block_ptr)
			score=dble(ptr%ptr%row_group)
		class default
			write(*,*)'unexpected item type in node_score_dble'
			stop
	end select
end function node_score_block_ptr_row





subroutine element_Zmn_block_user(nrow,ncol,mrange,nrange,values,msh,option,ker,myflag,passflag,ptree,stats)

	use BPACK_DEFS
	implicit none

	integer ii, jj,nn,pp,ij,i,j,nrow,ncol,passflag,myflag,Ninter,idx,nc,nr,pgno,ctxt,nprow,npcol,myrow,mycol
	integer mrange(nrow)
	integer nrange(ncol)
	DT:: value_e,values(nrow,ncol)
	type(mesh)::msh
	type(proctree)::ptree
	type(Hoption)::option
	type(Hstat)::stats
	type(kernelquant)::ker
	integer ierr,idx_row,idx_col,idx_dat
	integer,allocatable:: flags(:),dests(:),colidx1(:),rowidx1(:),colidx(:),rowidx(:),allrows(:),allcols(:),disps(:),pgidx(:),pmaps(:,:)
	procedure(F_Zelem_block), POINTER :: proc
	procedure(C_Zelem_block), POINTER :: proc_c
	procedure(F_Zelem), POINTER :: proc1
	procedure(C_Zelem), POINTER :: proc1_c
	DT,allocatable::alldat_loc(:)
	type(intersect),allocatable::inters(:)
	integer myArows,myAcols,Npmap
	real(kind=8)::t1,t2,t3,t4
	integer reqm,reqn
	integer statusm(MPI_status_size),statusn(MPI_status_size)

	if(option%elem_extract==0)then

		t1 = OMP_get_wtime()

		if(option%cpp==1)then
			call c_f_procpointer(ker%C_FuncZmn, proc1_C)
#ifdef HAVE_TASKLOOP
			!$omp parallel
			!$omp single
			!$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
			!$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
			do ij=1,ncol*nrow
				jj = (ij-1)/nrow+1
				ii = mod(ij-1,nrow) + 1
				value_e=0
				call proc1_C(msh%new2old(mrange(ii)),msh%new2old(nrange(jj)),value_e,ker%C_QuantApp)
				value_e =value_e*option%scale_factor
				values(ii,jj) = value_e
			enddo
#ifdef HAVE_TASKLOOP
			!$omp end taskloop
			!$omp end single
			!$omp end parallel
#else
			!$omp end parallel do
#endif
		else
			proc1 => ker%FuncZmn
			if(nrow*ncol>0)then
#ifdef HAVE_TASKLOOP
				!$omp parallel
				!$omp single
				!$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
				!$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
				do ij=1,ncol*nrow
					jj = (ij-1)/nrow+1
					ii = mod(ij-1,nrow) + 1
					value_e=0
					call proc1(msh%new2old(mrange(ii)),msh%new2old(nrange(jj)),value_e,ker%QuantApp)
					value_e =value_e*option%scale_factor
					values(ii,jj)=value_e
				enddo
#ifdef HAVE_TASKLOOP
			!$omp end taskloop
			!$omp end single
			!$omp end parallel
#else
			!$omp end parallel do
#endif

			endif
		endif

		t2 = OMP_get_wtime()
		stats%Time_Entry=stats%Time_Entry+t2-t1

		passflag=2
	else if(option%elem_extract==1)then

		allocate(flags(ptree%nproc))


#ifdef HAVE_MPI3
		call MPI_IALLGATHER(myflag, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
		call MPI_Wait(reqm,statusm,ierr)
#else
		call MPI_ALLGATHER(myflag, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif




		passflag=minval(flags)

		t1 = OMP_get_wtime()

		if(passflag==0)then

			allocate(colidx1(ptree%nproc))
			allocate(rowidx1(ptree%nproc))
			allocate(disps(ptree%nproc))


#ifdef HAVE_MPI3
			call MPI_IALLGATHER(nrow, 1, MPI_INTEGER, rowidx1, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
			call MPI_IALLGATHER(ncol, 1, MPI_INTEGER, colidx1, 1, MPI_INTEGER, ptree%Comm, reqn, ierr)
#else
			call MPI_ALLGATHER(nrow, 1, MPI_INTEGER, rowidx1, 1, MPI_INTEGER, ptree%Comm, ierr)
			call MPI_ALLGATHER(ncol, 1, MPI_INTEGER, colidx1, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif



			Npmap=ptree%nproc
			allocate(pmaps(Npmap,3))
			do pp=1,Npmap
				pmaps(pp,1)=1
				pmaps(pp,2)=1
				pmaps(pp,3)=pp-1
			enddo

			Ninter=0
			do pp=1,ptree%nproc
			if(flags(pp)==0)then
			Ninter=Ninter+1
			endif
			enddo

			allocate(colidx(Ninter))
			allocate(rowidx(Ninter))
			allocate(pgidx(Ninter))

#ifdef HAVE_MPI3
			call MPI_Wait(reqm,statusm,ierr)
			call MPI_Wait(reqn,statusn,ierr)
#endif
			!***** Count number of active intersections Ninter
			Ninter=0
			do pp=1,ptree%nproc
			if(flags(pp)==0)then
			Ninter=Ninter+1
			pgidx(Ninter)=pp
			rowidx(Ninter)=rowidx1(pp)
			colidx(Ninter)=colidx1(pp)
			endif
			enddo

			!***** count number of local data
			idx_dat=0
			do nn=1,Ninter
				nr=rowidx(nn)
				nc=colidx(nn)
				! datidx(nn)=ntot_loc
				nprow=pmaps(pgidx(nn),1)
				npcol=pmaps(pgidx(nn),2)
				call Gridinfo_2D(pmaps(pgidx(nn),:),ptree%MyID,myrow,mycol)
				if(myrow/=-1 .and. mycol/=-1)then
					myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
					myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
					idx_dat = idx_dat + myArows*myAcols
				endif
			enddo
			allocate(alldat_loc(idx_dat))
			if(idx_dat>0)alldat_loc=0


			!***** Broadcast mrange and nrange for each intersection
			idx_row=sum(rowidx)
			allocate(allrows(idx_row))
			idx=0
			do pp=1,ptree%nproc
				disps(pp)=idx
				idx=idx+rowidx1(pp)
			enddo

#ifdef HAVE_MPI3
			call MPI_IALLGATHERV(mrange, nrow, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, reqm, ierr)
#else
			call MPI_ALLGATHERV(mrange, nrow, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif

			idx_col=sum(colidx)
			allocate(allcols(idx_col))
			idx=0
			do pp=1,ptree%nproc
				disps(pp)=idx
				idx=idx+colidx1(pp)
			enddo
#ifdef HAVE_MPI3
			call MPI_IALLGATHERV(nrange, ncol, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, reqn, ierr)
			call MPI_Wait(reqm,statusm,ierr)
			call MPI_Wait(reqn,statusn,ierr)
#else
			call MPI_ALLGATHERV(nrange, ncol, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif
			if(option%cpp==1)then
				call c_f_procpointer(ker%C_FuncZmnBlock, proc_C)
				! !***** parallel extraction of the data
				do ii=1,idx_row
					allrows(ii)=abs(msh%new2old(allrows(ii)))
				enddo
				do jj=1,idx_col
					allcols(jj)=abs(msh%new2old(allcols(jj)))
				enddo
				pgidx=pgidx-1
				call proc_C(Ninter,idx_row,idx_col,idx_dat,allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps,ker%C_QuantApp)
			else
				proc => ker%FuncZmnBlock
				! !***** parallel extraction of the data
				do ii=1,idx_row
					allrows(ii)=abs(msh%new2old(allrows(ii)))
				enddo
				do jj=1,idx_col
					allcols(jj)=abs(msh%new2old(allcols(jj)))
				enddo
				call proc(Ninter,allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps,ker%QuantApp)
			endif


			idx=0
			do jj=1,ncol ! note that alldat_loc has column major
			do ii=1,nrow
				idx=idx+1
				values(ii,jj)=alldat_loc(idx)
			enddo
			enddo

			deallocate(allrows)
			deallocate(allcols)
			deallocate(colidx1)
			deallocate(colidx)
			deallocate(rowidx1)
			deallocate(rowidx)
			deallocate(disps)
			deallocate(alldat_loc)
			deallocate(pgidx)
			deallocate(pmaps)

		endif
		t2 = OMP_get_wtime()
		stats%Time_Entry=stats%Time_Entry+t2-t1
		deallocate(flags)
	endif

	return

end subroutine element_Zmn_block_user


end module Bplus_Utilities