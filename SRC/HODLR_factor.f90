#include "HODLR_config.fi"
module HODLR_factor
use Bplus_update
use Bplus_inversion


contains 

subroutine HODLR_Factorization(ho_bf1,option,stats,ptree)

    use HODLR_DEFS
    
    	
	use misc
	use omp_lib
    implicit none

    integer i, j, ii, jj, iii, jjj,index_ij,mm,nn
    integer level, blocks, edge, patch, node, group,level_c,groupm_diag
    integer rank, index_near, m, n, length, flag, itemp
    real T0
	real(kind=8) rtemp,tmpfact
    real(kind=8) Memory, Memory_near
	integer,allocatable:: index_old(:),index_new(:) 
	integer::block_num,block_num_new,num_blocks,level_butterfly	
	integer, allocatable :: ipiv(:)
	integer rowblock,pgno1,pgno2,pgno,ierr,rowblock_inv
	type(matrixblock),pointer::block_o,block_off,block_off1,block_off2
	type(matrixblock)::block_tmp
	real(kind=8) n1,n2,flop
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree


    if(ptree%MyID==Main_ID)write (*,*) ''
	
	call MPI_barrier(ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,*) 'Computing block inverse at level Maxlevel_for_blocks+1...'	
	level_c = ho_bf1%Maxlevel+1
	do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
			
		ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat
		nn = size(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,1)
		allocate(ipiv(nn))
		call getrff90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		call getrif90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv,flop=flop)		
		stats%Flop_Factor = stats%Flop_Factor + flop
		!!!!!!! the forward block BP can be deleted if not used in solution phase
		
		stats%Mem_Direct_inv=stats%Mem_Direct_inv+SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3	

		deallocate(ipiv)
	end do		


	call MPI_barrier(ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) 'Computing block inverse at higher levels...'
	do level_c = ho_bf1%Maxlevel,1,-1

		!!!***** update the forward off-diagonal block by left multiplication of inverse of diagonal blocks in Z: Z_ij^l -> Z_ii^-1*Z_ij^l
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write(*,*)'update forward blocks at level:',level_c
					
		n1 = OMP_get_wtime()
		do rowblock_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		do rowblock=rowblock_inv*2-1,rowblock_inv*2
		
		if(ptree%MyID >=ptree%pgrp(ho_bf1%levels(level_c)%BP(rowblock)%pgno)%head .and. ptree%MyID <=ptree%pgrp(ho_bf1%levels(level_c)%BP(rowblock)%pgno)%tail)then			
			
			call Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats,ptree) 
			
			call ComputeMemory_Bplus(ho_bf1%levels(level_c)%BP_inverse_update(rowblock),rtemp)
			stats%Mem_Sblock = stats%Mem_Sblock + rtemp
			
			
			! if(level_c==6)then			
				! call print_butterfly_size_rank(ho_bf1%levels(level_c)%matrices_block(rowblock),option%tol_comp)
				! stop
			! end if	
			
		end if	
		end do
		end do
		n2 = OMP_get_wtime()
		stats%Time_Sblock=stats%Time_Sblock+n2-n1 
	

		!!!***** compute the inverse of each block 2x2 submatrices whose two off-diagonal blocks are butterflies 
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write(*,*)'compute block inverse at level:',level_c
		n1 = OMP_get_wtime()
		do rowblock = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
			
			pgno =  ho_bf1%levels(level_c)%BP_inverse(rowblock)%pgno		
			if((ptree%MyID >=ptree%pgrp(pgno)%head .and. ptree%MyID <=ptree%pgrp(pgno)%tail))then	
				
				call DoubleDistributeBplus(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2-1),stats,ptree)
				call DoubleDistributeBplus(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2),stats,ptree)
				
				call Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,stats,ptree)
				call ComputeMemory_Bplus(ho_bf1%levels(level_c)%BP_inverse_schur(rowblock),rtemp)
				stats%Mem_SMW=stats%Mem_SMW+rtemp

			endif
		end do
		n2 = OMP_get_wtime()		
		stats%Time_Inv=stats%Time_Inv + n2-n1	
	end do
	
	call MPI_ALLREDUCE(stats%Time_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,*) 'computing updated forward block time:',rtemp,'Seconds'		
	call MPI_ALLREDUCE(stats%Time_Inv,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,*) 'computing inverse block time:',rtemp,'Seconds'	
	call MPI_ALLREDUCE(stats%Time_random(1),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_Init:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(2),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_MVP:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(3),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_Reconstruct:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(4),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_Onesub:', rtemp
	call MPI_ALLREDUCE(stats%Time_SMW,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_SMW:', rtemp
	call MPI_ALLREDUCE(stats%Time_RedistB,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistB:', rtemp	
	call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistV:', rtemp		
	call MPI_ALLREDUCE(time_tmp,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*)'time_tmp',time_tmp
	call MPI_ALLREDUCE(stats%Flop_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21Es14.2)') 'Factorization flops:',rtemp	
	
    if(ptree%MyID==Main_ID)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_SMW,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for butterfly inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for butterfly Sblocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_inv,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for direct inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_int_vec,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for storing intermidiate vectors'
    if(ptree%MyID==Main_ID)write(*,*)''	

	
    return

end subroutine HODLR_Factorization


end module HODLR_factor