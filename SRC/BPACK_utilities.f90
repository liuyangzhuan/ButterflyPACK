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
module BPACK_Utilities
use MISC_Utilities
use Bplus_Utilities

contains



subroutine copy_HOBF(ho_bf_i,ho_bf_o)
use BPACK_DEFS
use MISC_Utilities
implicit none

type(hobf)::ho_bf_i,ho_bf_o

integer ii
integer level_c

! real(kind=8),optional::memory
! real(kind=8)::rtemp


ho_bf_o%Maxlevel = ho_bf_i%Maxlevel
ho_bf_o%N = ho_bf_i%N


allocate(ho_bf_o%levels(ho_bf_o%Maxlevel+1))
do level_c = 1,ho_bf_o%Maxlevel+1
	ho_bf_o%levels(level_c)%level = ho_bf_i%levels(level_c)%level
	ho_bf_o%levels(level_c)%N_block_forward = ho_bf_i%levels(level_c)%N_block_forward
	ho_bf_o%levels(level_c)%N_block_inverse = ho_bf_i%levels(level_c)%N_block_inverse
	ho_bf_o%levels(level_c)%Bidxs = ho_bf_i%levels(level_c)%Bidxs
	ho_bf_o%levels(level_c)%Bidxe = ho_bf_i%levels(level_c)%Bidxe

	allocate(ho_bf_o%levels(level_c)%BP(ho_bf_o%levels(level_c)%N_block_forward))
	! write(*,*)ho_bf_o%levels(level_c)%N_block_inverse,'g'
	allocate(ho_bf_o%levels(level_c)%BP_inverse(ho_bf_o%levels(level_c)%N_block_inverse))
	! write(*,*)ho_bf_o%levels(level_c)%N_block_inverse,'g1'
	do ii = 1, ho_bf_o%levels(level_c)%N_block_forward
		call Bplus_copy(ho_bf_i%levels(level_c)%BP(ii),ho_bf_o%levels(level_c)%BP(ii))
	end do
	! if(level_c/=ho_bf_o%Maxlevel+1)then
		do ii = 1, ho_bf_o%levels(level_c)%N_block_inverse
			! write(*,*)ii,'6642'
			call Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse(ii),ho_bf_o%levels(level_c)%BP_inverse(ii))
		end do
	! end if
end do

end subroutine copy_HOBF



recursive subroutine Hmat_delete_global_tree(blocks)
    implicit none
    type(global_matricesblock), pointer :: blocks, blocks_son
    integer group_m, group_n, i, j, k, level


    group_m=blocks%row_group
    group_n=blocks%col_group
    level=blocks%level

    if (associated(blocks%sons)) then
       	blocks_son=>blocks%sons(1,1)
		call Hmat_delete_global_tree(blocks_son)
		blocks_son=>blocks%sons(2,1)
		call Hmat_delete_global_tree(blocks_son)
		blocks_son=>blocks%sons(1,2)
		call Hmat_delete_global_tree(blocks_son)
		blocks_son=>blocks%sons(2,2)
		call Hmat_delete_global_tree(blocks_son)
		deallocate(blocks%sons)
    endif

    return

end subroutine Hmat_delete_global_tree


subroutine Hmat_delete(h_mat)
use BPACK_DEFS
use MISC_Utilities
implicit none

type(Hmat)::h_mat
integer bm,bn,ii,jj,level

call Hmat_delete_global_tree(h_mat%blocks_root)
deallocate(h_mat%blocks_root)


if(associated(h_mat%First_block_eachlevel))deallocate(h_mat%First_block_eachlevel)
if(associated(h_mat%Local_blocks))then
	bm = size(h_mat%Local_blocks,1)
	bn = size(h_mat%Local_blocks,2)
	do ii=1,bm
	do jj=1,bn
		call Hmat_block_delete(h_mat%Local_blocks(ii,jj))
	enddo
	enddo
	deallocate(h_mat%Local_blocks)
endif

if(associated(h_mat%Local_blocks_copy))then
	bm = size(h_mat%Local_blocks_copy,1)
	bn = size(h_mat%Local_blocks_copy,2)
	do ii=1,bm
	do jj=1,bn
		call Hmat_block_delete(h_mat%Local_blocks_copy(ii,jj))
	enddo
	enddo
	deallocate(h_mat%Local_blocks_copy)
endif

if(allocated(h_mat%lstblks))then
	do level=0,h_mat%Maxlevel
		call list_finalizer(h_mat%lstblks(level))
	enddo
	deallocate(h_mat%lstblks)
endif

end subroutine Hmat_delete



subroutine HODLR_delete(ho_bf_o)
use BPACK_DEFS
use MISC_Utilities
implicit none

type(hobf)::ho_bf_o

integer ii
integer level_c

do level_c = 1,ho_bf_o%Maxlevel+1
	do ii = 1, ho_bf_o%levels(level_c)%N_block_forward
		call Bplus_delete(ho_bf_o%levels(level_c)%BP(ii))
		call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse_update(ii))
	end do
	do ii = 1, ho_bf_o%levels(level_c)%N_block_inverse
		call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse(ii))
		call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse_schur(ii))
	end do
	deallocate(ho_bf_o%levels(level_c)%BP)
	deallocate(ho_bf_o%levels(level_c)%BP_inverse_update)
	deallocate(ho_bf_o%levels(level_c)%BP_inverse)
	deallocate(ho_bf_o%levels(level_c)%BP_inverse_schur)
end do
deallocate(ho_bf_o%levels)

end subroutine HODLR_delete

subroutine BPACK_delete(bmat)
use BPACK_DEFS
implicit none
type(Bmatrix)::bmat
if(associated(bmat%ho_bf))then
	call HODLR_delete(bmat%ho_bf)
	deallocate(bmat%ho_bf)
	bmat%ho_bf=>null()
endif
if(associated(bmat%h_mat))then
	call Hmat_delete(bmat%h_mat)
	deallocate(bmat%h_mat)
	bmat%h_mat=>null()
endif
end subroutine BPACK_delete


subroutine delete_kernelquant(ker)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(kernelquant)::ker
if(allocated(ker%matZ_glo))deallocate(ker%matZ_glo)
end subroutine delete_kernelquant


subroutine delete_mesh(msh)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(mesh)::msh
integer ii

if(allocated(msh%xyz))deallocate(msh%xyz)
if(allocated(msh%nns))deallocate(msh%nns)
if(allocated(msh%new2old))deallocate(msh%new2old)
if(allocated(msh%old2new))deallocate(msh%old2new)
if(allocated(msh%pretree))deallocate(msh%pretree)
if(allocated(msh%basis_group))then
do ii=1,msh%Maxgroup
	if(allocated(msh%basis_group(ii)%center))deallocate(msh%basis_group(ii)%center)
	if(allocated(msh%basis_group(ii)%nlist))deallocate(msh%basis_group(ii)%nlist)
	msh%basis_group(ii)%nn=0
enddo
deallocate(msh%basis_group)
endif

end subroutine delete_mesh




subroutine delete_proctree(ptree)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(proctree)::ptree
integer ii,Maxgrp
integer ierr

if(allocated(ptree%pgrp))then
Maxgrp=2**(ptree%nlevel)-1
do ii=1,Maxgrp
	! if(associated(ptree%pgrp(ii)%gd))then
		! call delete_grid(ptree%pgrp(ii)%gd)
		! deallocate(ptree%pgrp(ii)%gd)
		! ptree%pgrp(ii)%gd=>null()
	! endif
	if(ptree%pgrp(ii)%ctxt/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt)
	if(ptree%pgrp(ii)%ctxt1D/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt1D)
	if(ptree%pgrp(ii)%ctxt1DCol/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt1DCol)
	if(ptree%pgrp(ii)%ctxt_head/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt_head)
	if(ptree%pgrp(ii)%Comm/=MPI_COMM_NULL)call MPI_Comm_free(ptree%pgrp(ii)%Comm,ierr)
enddo
deallocate(ptree%pgrp)
endif
if(ptree%Comm/=MPI_COMM_NULL)call MPI_Comm_free(ptree%Comm,ierr)

end subroutine delete_proctree


! recursive subroutine delete_grid(gd)
! use BPACK_DEFS
! use MISC_Utilities
! implicit none
! type(grid)::gd
! integer ierr

! if(.not. associated(gd%gdc))then
	! if(gd%ctxt/=-1)call blacs_gridexit(gd%ctxt)
	! if(gd%Comm/=MPI_COMM_NULL)call MPI_Comm_free(gd%Comm,ierr)
	! return
! else
	! call delete_grid(gd%gdc(1))
	! call delete_grid(gd%gdc(2))
	! deallocate(gd%gdc)
	! gd%gdc=>null()
! endif
! end subroutine delete_grid


subroutine delete_Hstat(stats)
use BPACK_DEFS
use MISC_Utilities
implicit none
type(Hstat)::stats

if(allocated(stats%rankmax_of_level))deallocate(stats%rankmax_of_level)
if(allocated(stats%rankmin_of_level))deallocate(stats%rankmin_of_level)
if(allocated(stats%rankmax_of_level_global))deallocate(stats%rankmax_of_level_global)
if(allocated(stats%Add_random_CNT))deallocate(stats%Add_random_CNT)
if(allocated(stats%Mul_random_CNT))deallocate(stats%Mul_random_CNT)
if(allocated(stats%XLUM_random_CNT))deallocate(stats%XLUM_random_CNT)
if(allocated(stats%Add_random_Time))deallocate(stats%Add_random_Time)
if(allocated(stats%Mul_random_Time))deallocate(stats%Mul_random_Time)
if(allocated(stats%XLUM_random_Time))deallocate(stats%XLUM_random_Time)


end subroutine delete_Hstat

recursive subroutine copy_basis_group(basis_group1,node1,Maxgroup1,basis_group2,node2,Maxgroup2,offset)
implicit none
type(basisgroup):: basis_group1(:),basis_group2(:)
integer node1,node2,Maxgroup1,Maxgroup2,offset
if(node2<=Maxgroup2 .and. node1<=Maxgroup1)then

	basis_group2(node2)%head =basis_group1(node1)%head+offset
	basis_group2(node2)%tail =basis_group1(node1)%tail+offset
	basis_group2(node2)%pgno =basis_group1(node1)%pgno

	call copy_basis_group(basis_group1,node1*2,Maxgroup1,basis_group2,node2*2,Maxgroup2,offset)
	call copy_basis_group(basis_group1,node1*2+1,Maxgroup1,basis_group2,node2*2+1,Maxgroup2,offset)
endif

end subroutine copy_basis_group





subroutine InitStat(stats)
	implicit none
	type(Hstat)::stats

	stats%Time_random=0  ! Intialization, MVP, Reconstruction
	stats%Time_Sblock=0
	stats%Time_Sol=0
	stats%Time_C_Mult=0
	stats%Time_C_Extract=0
	stats%Time_Inv=0
	stats%Time_RedistB=0
	stats%Time_RedistV=0
	stats%Time_SMW=0
	stats%Time_Fill=0
	stats%Time_Entry=0
	stats%Mem_peak=0
	stats%Mem_Sblock=0
	stats%Mem_SMW=0
	stats%Mem_Direct_for=0
	stats%Mem_Direct_inv=0
	stats%Mem_int_vec=0
	stats%Mem_Comp_for=0
	stats%Mem_Fill=0
	stats%Mem_Factor=0
	stats%Flop_Fill=0
	stats%Flop_Factor=0
	stats%Flop_Sol=0
	stats%Flop_C_Mult=0
	stats%Flop_C_Extract=0

	stats%Time_Direct_LU=0
	stats%Time_Add_Multiply=0
	stats%Time_Multiply=0
	stats%Time_XLUM=0
	stats%Time_Split=0
	stats%Time_Comm=0
	stats%Time_Idle=0
	stats%Time_Factor=0

	time_tmp = 0
end subroutine InitStat



subroutine PrintStat(stats,ptree)
	implicit none
	type(Hstat)::stats
	type(proctree)::ptree
	real(kind=8)::rtemp,rtemp1,rtemp2
	integer ierr


	! stats%Time_random=0  ! Intialization, MVP, Reconstruction
	! stats%Time_Sblock=0
	! stats%Time_Sol=0
	! stats%Time_Inv=0
	! stats%Time_RedistB=0
	! stats%Time_RedistV=0
	! stats%Time_SMW=0
	! stats%Time_Fill=0
	! stats%Mem_peak=0
	! stats%Mem_Sblock=0
	! stats%Mem_SMW=0
	! stats%Mem_Direct_for=0
	! stats%Mem_Direct_inv=0
	! stats%Mem_int_vec=0
	! stats%Mem_Comp_for=0
	! stats%Flop_Fill=0
	! stats%Flop_Factor=0
	! stats%Flop_Sol=0



	call MPI_ALLREDUCE(stats%Time_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Constr time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Time_Entry,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'EntryEval time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A3)') 'Tot constr mem:',rtemp+rtemp1,'MB'
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Constr flops:',rtemp
	if(ptree%MyID==Main_ID)write (*,'(A21,I14)') 'Rank before factor:', maxval(stats%rankmax_of_level_global)


	call MPI_ALLREDUCE(stats%Time_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Factor time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Mem_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A3)') 'Tot factor mem:',rtemp,'MB'
	call MPI_ALLREDUCE(stats%Flop_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Factor flops:',rtemp

	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Solve time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Solve flops:',rtemp


	call MPI_ALLREDUCE(stats%Time_C_Mult,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'C_mult time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_C_Mult,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'C_mult flops:',rtemp

	call MPI_ALLREDUCE(stats%Time_C_Extract,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'C_extract time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_C_Extract,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'C_extract flops:',rtemp

	call MPI_ALLREDUCE(stats%Mem_peak,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A3)') 'Peak mem:',rtemp,'MB'

end subroutine PrintStat



subroutine SetDefaultOptions(option)
	implicit none
	type(Hoption)::option

	option%Nmin_leaf=200
	option%tol_comp=1d-4
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%n_iter=1000
	option%tol_rand=option%tol_comp
	option%tol_Rdetect=option%tol_comp*1d-1
	option%level_check=10000
	option%precon=DIRECT
	option%xyzsort=TM
	option%lnoBP=40000
	option%TwoLayerOnly=1
	option%touch_para=0d0
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=0
	option%BACA_Batch=16
	option%RecLR_leaf=BACA
	option%nogeo=0
	option%ErrSol=0
	option%LR_BLK_NUM=1
	option%rank0=32
	option%rankrate=1.5d0
	option%itermax=10
	option%powiter=0
	option%ILU=0
	option%Nbundle=1
	option%near_para=SafeEps
	option%format=HODLR
	option%verbosity=0
	option%scale_factor=1d0
	option%rmax=3000
	option%forwardN15flag=0
	option%sample_para=1.2d0
	option%sample_heuristic=1
	option%pat_comp=3
	option%elem_extract=0
	option%knn=0
	option%cpp=0
	option%bp_cnt_lr=0
	option%less_adapt=0

end subroutine SetDefaultOptions


subroutine ReadOption(option,ptree,ii)
	implicit none
	type(Hoption)::option
	type(proctree)::ptree
	integer ii
	integer nargs,flag
	character(len=1024)  :: strings,strings1

	nargs = iargc()
	flag=1
	do while(flag==1)
		ii=ii+1
		if(ii<=nargs)then
			call getarg(ii,strings)
			if(strings(1:2)=='--')then
				ii=ii+1
				call getarg(ii,strings1)

				if(trim(strings)=='--nmin_leaf')then
					read(strings1,*)option%Nmin_leaf
				else if	(trim(strings)=='--tol_comp')then
					read(strings1,*)option%tol_comp
					option%tol_rand=option%tol_comp
					option%tol_Rdetect=option%tol_comp*1d-1
				else if	(trim(strings)=='--tol_itersol')then
					read(strings1,*)option%tol_itersol
				else if	(trim(strings)=='--n_iter')then
					read(strings1,*)option%n_iter
				else if	(trim(strings)=='--level_check')then
					read(strings1,*)option%level_check
				else if	(trim(strings)=='--precon')then
					read(strings1,*)option%precon
				else if	(trim(strings)=='--xyzsort')then
					read(strings1,*)option%xyzsort
				else if	(trim(strings)=='--schulzorder')then
					read(strings1,*)option%schulzorder
				else if	(trim(strings)=='--schulzlevel')then
					read(strings1,*)option%schulzlevel
				else if	(trim(strings)=='--lrlevel')then
					read(strings1,*)option%LRlevel
				else if	(trim(strings)=='--errfillfull')then
					read(strings1,*)option%ErrFillFull
				else if	(trim(strings)=='--baca_batch')then
					read(strings1,*)option%BACA_Batch
				else if	(trim(strings)=='--reclr_leaf')then
					read(strings1,*)option%RecLR_leaf
				else if	(trim(strings)=='--nogeo')then
					read(strings1,*)option%nogeo
				else if	(trim(strings)=='--less_adapt')then
					read(strings1,*)option%less_adapt
				else if	(trim(strings)=='--errsol')then
					read(strings1,*)option%ErrSol
				else if	(trim(strings)=='--lr_blk_num')then
					read(strings1,*)option%LR_BLK_NUM
				else if	(trim(strings)=='--rank0')then
					read(strings1,*)option%rank0
				else if	(trim(strings)=='--rankrate')then
					read(strings1,*)option%rankrate
				else if	(trim(strings)=='--itermax')then
					read(strings1,*)option%itermax
				else if	(trim(strings)=='--powiter')then
					read(strings1,*)option%powiter
				else if	(trim(strings)=='--ilu')then
					read(strings1,*)option%ILU
				else if	(trim(strings)=='--nbundle')then
					read(strings1,*)option%Nbundle
				else if	(trim(strings)=='--near_para')then
					read(strings1,*)option%near_para
				else if	(trim(strings)=='--format')then
					read(strings1,*)option%format
				else if	(trim(strings)=='--verbosity')then
					read(strings1,*)option%verbosity
				else if	(trim(strings)=='--rmax')then
					read(strings1,*)option%rmax
				else if	(trim(strings)=='--sample_para')then
					read(strings1,*)option%sample_para
				else if	(trim(strings)=='--sample_heuristic')then
					read(strings1,*)option%sample_heuristic
				else if	(trim(strings)=='--pat_comp')then
					read(strings1,*)option%pat_comp
				else if	(trim(strings)=='--elem_extract')then
					read(strings1,*)option%elem_extract
				else if	(trim(strings)=='--knn')then
					read(strings1,*)option%knn
				else if	(trim(strings)=='--cpp')then
					read(strings1,*)option%cpp
				else if	(trim(strings)=='--lnobp')then
					read(strings1,*)option%lnoBP
				else if	(trim(strings)=='--bp_cnt_lr')then
					read(strings1,*)option%bp_cnt_lr
				else if	(trim(strings)=='--touch_para')then
					read(strings1,*)option%touch_para
				else
					if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown option: ', trim(strings)
				endif
			else
				flag=0
			endif
		else
			flag=0
		endif
	enddo

end subroutine ReadOption


subroutine CopyOptions(option,option1)
	implicit none
	type(Hoption)::option,option1

	option1%Nmin_leaf = option%Nmin_leaf
	option1%tol_comp = option%tol_comp
	option1%tol_Rdetect = option%tol_Rdetect
	option1%tol_LS = option%tol_LS
	option1%tol_itersol = option%tol_itersol
	option1%n_iter = option%n_iter
	option1%tol_rand = option%tol_rand
	option1%level_check = option%level_check
	option1%precon = option%precon
	option1%xyzsort = option%xyzsort
	option1%lnoBP = option%lnoBP
	option1%bp_cnt_lr = option%bp_cnt_lr
	option1%TwoLayerOnly = option%TwoLayerOnly
	option1%touch_para = option%touch_para
    option1%schulzorder = option%schulzorder
    option1%schulzlevel = option%schulzlevel
	option1%LRlevel = option%LRlevel
	option1%ErrFillFull = option%ErrFillFull
	option1%BACA_Batch = option%BACA_Batch
	option1%RecLR_leaf = option%RecLR_leaf
	option1%nogeo = option%nogeo
	option1%ErrSol = option%ErrSol
	option1%LR_BLK_NUM = option%LR_BLK_NUM
	option1%rank0 = option%rank0
	option1%rankrate = option%rankrate
	option1%itermax = option%itermax
	option1%powiter = option%powiter
	option1%ILU = option%ILU
	option1%Nbundle = option%Nbundle
	option1%near_para = option%near_para
	option1%format = option%format
	option1%verbosity = option%verbosity
	option1%scale_factor = option%scale_factor
	option1%rmax = option%rmax
	option1%forwardN15flag = option%forwardN15flag
	option1%sample_para = option%sample_para
	option1%sample_heuristic = option%sample_heuristic
	option1%pat_comp = option%pat_comp
	option1%elem_extract = option%elem_extract
	option1%cpp = option%cpp
	option1%knn = option%knn
	option1%less_adapt = option%less_adapt

end subroutine CopyOptions



subroutine PrintOptions(option,ptree)
	implicit none
	type(Hoption)::option
	type(proctree)::ptree

	if(ptree%MyID==Main_ID)then
		write(*,*) ' '
		write(*,*) '***************************'
		write(*,'(A25)') 'Printing Solver Options:'
		write(*,'(A18,I8)') 'Nmin_leaf',option%Nmin_leaf
		write(*,'(A18,I8)') 'n_iter', option%n_iter
		write(*,'(A18,I8)') 'level_check', option%level_check
		write(*,'(A18,I8)') 'precon', option%precon
		write(*,'(A18,I8)') 'xyzsort', option%xyzsort
		write(*,'(A18,I8)') 'lnoBP', option%lnoBP
		write(*,'(A18,I8)') 'bp_cnt_lr', option%bp_cnt_lr
		write(*,'(A18,I8)') 'TwoLayerOnly', option%TwoLayerOnly
		write(*,'(A18,I8)') 'schulzorder', option%schulzorder
		write(*,'(A18,I8)') 'schulzlevel', option%schulzlevel
		write(*,'(A18,I8)') 'LRlevel', option%LRlevel
		write(*,'(A18,I8)') 'BACA_Batch', option%BACA_Batch
		write(*,'(A18,I8)') 'RecLR_leaf', option%RecLR_leaf
		write(*,'(A18,I8)') 'nogeo', option%nogeo
		write(*,'(A18,I8)') 'LR_BLK_NUM', option%LR_BLK_NUM
		write(*,'(A18,I8)') 'rank0', option%rank0
		write(*,'(A18,I8)') 'itermax', option%itermax
		write(*,'(A18,I8)') 'powiter', option%powiter
		write(*,'(A18,I8)') 'ILU', option%ILU
		write(*,'(A18,I8)') 'Nbundle', option%Nbundle
		write(*,'(A18,I8)') 'verbosity', option%verbosity
		write(*,'(A18,I8)') 'rmax', option%rmax
		write(*,'(A18,I8)') 'forwardN15flag', option%forwardN15flag
		write(*,'(A18,I8)') 'pat_comp', option%pat_comp
		write(*,'(A18,I8)') 'elem_extract', option%elem_extract
		write(*,'(A18,I8)') 'cpp', option%cpp
		write(*,'(A18,I8)') 'knn', option%knn
		write(*,'(A18,I8)') 'ErrFillFull', option%ErrFillFull
		write(*,'(A18,I8)') 'ErrSol', option%ErrSol
		write(*,'(A18,I8)') 'sample_heuristic', option%sample_heuristic
		write(*,'(A18,I8)') 'less_adapt', option%less_adapt

		write(*,'(A18,Es14.7)') 'rankrate', option%rankrate
		write(*,'(A18,Es14.7)') 'tol_comp', option%tol_comp
		write(*,'(A18,Es14.7)') 'tol_Rdetect', option%tol_Rdetect
		write(*,'(A18,Es14.7)') 'tol_LS', option%tol_LS
		write(*,'(A18,Es14.7)') 'tol_itersol', option%tol_itersol
		write(*,'(A18,Es14.7)') 'tol_rand', option%tol_rand
		write(*,'(A18,Es14.7)') 'touch_para', option%touch_para
		write(*,'(A18,Es14.7)') 'near_para', option%near_para
		write(*,'(A18,Es14.7)') 'scale_factor', option%scale_factor
		write(*,'(A18,Es14.7)') 'sample_para', option%sample_para
		write(*,*) '***************************'
		write(*,*) ' '
	endif
end subroutine PrintOptions


subroutine BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
use BPACK_DEFS
implicit none
integer v_major,v_minor,v_bugfix
v_major=BPACK_MAJOR_VERSION
v_minor=BPACK_MINOR_VERSION
v_bugfix=BPACK_PATCH_VERSION

end subroutine BPACK_GetVersionNumber

end module BPACK_Utilities
