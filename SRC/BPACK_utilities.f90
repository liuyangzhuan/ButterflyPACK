#include "HODLR_config.fi"
module BPACK_Utilities
use misc
use Bplus_Utilities

contains
 


subroutine copy_HOBF(ho_bf_i,ho_bf_o)
use BPACK_DEFS
use misc
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
	else 
		deallocate(blocks)
    endif
    
    return    

end subroutine Hmat_delete_global_tree


subroutine Hmat_delete(h_mat)
use BPACK_DEFS
use misc
implicit none 

type(Hmat)::h_mat
integer bm,bn,ii,jj

call Hmat_delete_global_tree(h_mat%blocks_root)

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

end subroutine Hmat_delete



subroutine HODLR_delete(ho_bf_o)
use BPACK_DEFS
use misc
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



subroutine delete_kernelquant(ker)
use BPACK_DEFS
use misc
implicit none 
type(kernelquant)::ker
if(allocated(ker%matZ_glo))deallocate(ker%matZ_glo)
end subroutine delete_kernelquant


subroutine delete_mesh(msh)
use BPACK_DEFS
use misc
implicit none 
type(mesh)::msh
integer ii

if(allocated(msh%xyz))deallocate(msh%xyz)
if(allocated(msh%new2old))deallocate(msh%new2old)
if(allocated(msh%old2new))deallocate(msh%old2new)
if(allocated(msh%pretree))deallocate(msh%pretree)
if(allocated(msh%basis_group))then
! do ii=1,msh%Maxgroup
	! if(allocated(msh%basis_group(ii)%center))deallocate(msh%basis_group(ii)%center)
! enddo
deallocate(msh%basis_group)
endif

end subroutine delete_mesh

subroutine delete_proctree(ptree)
use BPACK_DEFS
use misc
implicit none 
type(proctree)::ptree
integer ii,Maxgrp
integer ierr

if(allocated(ptree%pgrp))then
Maxgrp=2**(ptree%nlevel)-1		
do ii=1,Maxgrp
	if(associated(ptree%pgrp(ii)%gd))then
		call delete_grid(ptree%pgrp(ii)%gd)
		deallocate(ptree%pgrp(ii)%gd)
		ptree%pgrp(ii)%gd=>null()
	endif
	if(ptree%pgrp(ii)%ctxt/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt)
	if(ptree%pgrp(ii)%ctxt1D/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt1D)
	if(ptree%pgrp(ii)%ctxt_head/=-1)call blacs_gridexit(ptree%pgrp(ii)%ctxt_head)
	if(ptree%pgrp(ii)%Comm/=MPI_COMM_NULL)call MPI_Comm_free(ptree%pgrp(ii)%Comm,ierr)
enddo
deallocate(ptree%pgrp)
endif
if(ptree%Comm/=MPI_COMM_NULL)call MPI_Comm_free(ptree%Comm,ierr)

end subroutine delete_proctree


recursive subroutine delete_grid(gd)
use BPACK_DEFS
use misc
implicit none 
type(grid)::gd
integer ierr

if(.not. associated(gd%gdc))then
	if(gd%ctxt/=-1)call blacs_gridexit(gd%ctxt)
	if(gd%Comm/=MPI_COMM_NULL)call MPI_Comm_free(gd%Comm,ierr)
	return
else
	call delete_grid(gd%gdc(1))
	call delete_grid(gd%gdc(2))
	deallocate(gd%gdc)
	gd%gdc=>null()
endif
end subroutine delete_grid


subroutine delete_Hstat(stats)
use BPACK_DEFS
use misc
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
	stats%Time_Inv=0
	stats%Time_RedistB=0
	stats%Time_RedistV=0
	stats%Time_SMW=0
	stats%Time_Fill=0
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
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Construction time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A3)') 'Construction mem:',rtemp+rtemp1,'MB'	
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Construction flops:',rtemp



	call MPI_ALLREDUCE(stats%Time_Inv,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Factorization time:',rtemp,'Seconds'	
	call MPI_ALLREDUCE(stats%Mem_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	call MPI_ALLREDUCE(stats%Mem_SMW,rtemp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	call MPI_ALLREDUCE(stats%Mem_Direct_inv,rtemp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)	
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A3)') 'Factorization mem:',rtemp+rtemp1+rtemp2,'MB'	
	call MPI_ALLREDUCE(stats%Flop_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Factorization flops:',rtemp	

	
	
	
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'Solve time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'Solve flops:',rtemp


	call MPI_ALLREDUCE(stats%Time_C_Mult,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2,A8)') 'C_mult time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_C_Mult,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A21,Es14.2)') 'C_mult flops:',rtemp	
	
end subroutine PrintStat



subroutine SetDefaultOptions(option)
	implicit none 
	type(Hoption)::option	

	option%Nmin_leaf=200
	option%tol_comp=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%n_iter=1000
	option%tol_rand=1d-3
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
	option%BACA_Batch=64
	option%RecLR_leaf=BACA
	option%nogeo=0
	option%ErrSol=0
	option%LR_BLK_NUM=1
	option%rank0=32
	option%rankrate=1.2d0
	option%itermax=10
	option%powiter=0
	option%near_para=SafeEps
	option%format=HODLR
	option%verbosity=0

end subroutine SetDefaultOptions	



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
	option1%near_para = option%near_para
	option1%format = option%format
	option1%verbosity = option%verbosity

end subroutine CopyOptions	



end module BPACK_Utilities
