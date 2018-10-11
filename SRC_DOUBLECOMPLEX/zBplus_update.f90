#include "HODLR_config.fi"
module z_Bplus_update

use z_Bplus_randomized
integer rankthusfarS

		
contains



!**** Update one off-diagonal block in HODLR compressed as 
! Bplus/Butterfly/LR by multiplying on it left the inverse of diagonal block
! If LR, call LR_Sblock; if butterfly, call BF_randomized; if Bplus, call z_Bplus_randomized_constr  
	!ho_bf1: working HODLR     
	!level_c: level# of the block in HODLR   
	!rowblock: block# of the block at this level in HODLR    
	!option: containing compression options
	!stats: statistics 
	!ptree: process tree
subroutine Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats,ptree)

    use z_HODLR_DEFS
	use omp_lib
    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt
    character chara
    real(kind=8) T0
    type(blockplus),pointer::bplus
	type(matrixblock)::block_old
	type(matrixblock),pointer::block_o
    integer::rank_new_max
	real(kind=8)::rank_new_avr,error,rate,rankrate_inner,rankrate_outter 
	integer niter,rank,ntry,rank0,rank0_inner,rank0_outter
	real(kind=8):: error_inout
	real(kind=8):: n1,n2
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree
	
	call copy_Bplus(ho_bf1%levels(level_c)%BP(rowblock),ho_bf1%levels(level_c)%BP_inverse_update(rowblock))
	!!!!!!! the forward block BP can be deleted if not used in solution phase
	
	
    bplus =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
	if(bplus%Lplus==1)then
	
		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1) 									 
		level_butterfly=block_o%level_butterfly

		if(level_butterfly==0)then
			call LR_Sblock(ho_bf1,level_c,rowblock,ptree,stats)
		else 
			ho_bf1%ind_lv=level_c
			ho_bf1%ind_bk=rowblock
			rank0 = block_o%rankmax
			rate = 1.2d0
			call BF_randomized(level_butterfly,rank0,rate,block_o,ho_bf1,BF_block_MVP_Sblock_dat,error_inout,'Sblock',option,stats,ptree)			
#if PRNTlevel >= 1
			write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout	
#endif
		end if
	else 

		ho_bf1%ind_lv=level_c
		ho_bf1%ind_bk=rowblock
		Bplus =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1) 		
		
		rank0_inner = Bplus%LL(2)%rankmax
		rankrate_inner = 2.0d0
		
		rank0_outter = block_o%rankmax
		rankrate_outter=1.2d0	
		level_butterfly = block_o%level_butterfly
		call z_Bplus_randomized_constr(level_butterfly,Bplus,ho_bf1,rank0_inner,rankrate_inner,Bplus_block_MVP_Sblock_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Sblock_dat,error_inout,'Sblock+',option,stats,ptree)
		
		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1) 	
#if PRNTlevel >= 1
		write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout	
#endif		
		
	end if
	
    return

end subroutine Bplus_Sblock_randomized_memfree




subroutine LR_Sblock(ho_bf1,level_c,rowblock,ptree,stats)

    use z_HODLR_DEFS
    
	use z_misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,pp,qq
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
	type(matrixblock),pointer::block_o,blocks
	
    type(vectorsblock), pointer :: random1, random2
    
    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_end_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng,unique_nth,level_left_start
	integer*8 idx_start   
    integer level_blocks,head,tail
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	real(kind=8)::n2,n1
	type(hobf)::ho_bf1
	type(proctree)::ptree
	type(Hstat)::stats

	
	
	block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1) 
	  
    level_butterfly=block_o%level_butterfly
    call assert(level_butterfly==0,'Butterfly_Sblock_LowRank only works with LowRank blocks')
		
	
	
	num_blocks=2**level_butterfly

	
	num_vect_sub = size(block_o%ButterflyU%blocks(1)%matrix,2)
    ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   

	! get the right multiplied vectors
	pp = ptree%myid-ptree%pgrp(block_o%pgno)%head+1
	idx_start_glo = block_o%headm + block_o%M_p(pp,1) -1

	
	
	! mm=block_o%M
	mm=block_o%M_loc
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))
	vec_old = block_o%ButterflyU%blocks(1)%matrix
	stats%Flop_Tmp=0
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
				call fullmat_block_MVP_dat(blocks,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),cone,czero)
			else 
				call BF_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
			endif	
			endif
			endif
		end do		
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1			

		
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	block_o%ButterflyU%blocks(1)%matrix = vec_new
	deallocate(vec_old)
	deallocate(vec_new)

	stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
	
    return                

end subroutine LR_Sblock



end module z_Bplus_update
