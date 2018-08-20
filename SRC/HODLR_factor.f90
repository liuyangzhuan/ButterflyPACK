module cascading_factorization
! use Butterfly_inversion_schur_iterative
! use Butterfly_inversion_schur_partition
use Bplus_rightmultiply
use Bplus_inversion_schur_partition
contains 

subroutine cascading_factorizing(ho_bf1,option,stats,ptree)

    use MODULE_FILE
    ! use lapack95
    ! use blas95	
	use misc
	use omp_lib
    implicit none

    integer i, j, ii, jj, iii, jjj,index_ij,mm,nn
    integer level, blocks, edge, patch, node, group,level_c,groupm_diag
    integer rank, index_near, m, n, length, flag, itemp
    real T0
	real*8 rtemp,tmpfact
    real*8 Memory, Memory_near
	integer,allocatable:: index_old(:),index_new(:) 
	integer::block_num,block_num_new,num_blocks,level_butterfly	
	integer, allocatable :: ipiv(:)
	integer rowblock,pgno1,pgno2,pgno,ierr,rowblock_inv
	type(matrixblock),pointer::block_o,block_off,block_off1,block_off2
	type(matrixblock)::block_tmp
	real*8 n1,n2
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree
	! Memory_int_vec = 0
	! Memory_tfqmr_vec = 0

    level_c=0
    flag=0

    if(ptree%MyID==Main_ID)write (*,*) ''
	
	
	! Time_Init_forward=0
	! Time_Vector_forward=0
	! Time_Oneblock_forward=0					
	! Time_Reconstruct_forward=0
	! Time_Init_inverse=0
	! Time_Vector_inverse=0
	! Time_Reconstruct_inverse=0
	! Time_Oneblock_forward=0					
	! Time_InvertBlock = 0
	! time_tmp = 0
	
	! ! call Butterfly_forward_acc_check
	! ! stop
	call MPI_barrier(ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,*) 'Computing block inverse at level Maxlevel_for_blocks+1...'	
	level_c = ho_bf1%Maxlevel+1
	! allocate(ho_bf1%levels(level_c)%matrices_block_inverse(ho_bf1%levels(level_c)%N_block_inverse))
	do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
		! if(ptree%MyID >=ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)%head .and. ptree%MyID <=ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)%tail)then
			
			ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat
			nn = size(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,1)
			allocate(ipiv(nn))
			call getrff90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv)
			stats%Flop_Factor = stats%Flop_Factor + flops_zgetrf(nn,nn)
			call getrif90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv)		
			stats%Flop_Factor = stats%Flop_Factor + flops_zgetri(nn)
			! stats%Mem_Direct=stats%Mem_Direct+SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3		

			deallocate(ipiv)
		! endif	
	end do		


	call MPI_barrier(ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) 'Computing block inverse at higher levels...'
	do level_c = ho_bf1%Maxlevel,1,-1

		! update the forward butterfly after left multiplication of inverse 
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write(*,*)'update forward blocks at level:',level_c
					
		n1 = OMP_get_wtime()
		do rowblock_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		do rowblock=rowblock_inv*2-1,rowblock_inv*2
		
		if(ptree%MyID >=ptree%pgrp(ho_bf1%levels(level_c)%BP(rowblock)%pgno)%head .and. ptree%MyID <=ptree%pgrp(ho_bf1%levels(level_c)%BP(rowblock)%pgno)%tail)then			
			
			! if(level_c/=0)then	
			! if(level_c<ho_bf1%Maxlevel-8)then
				
				! call Butterfly_Sblock_randomized_memfree(level_c,rowblock,rtemp)
				! call Butterfly_Sblock_randomized_partialupdate_memfree(level_c,rowblock,rtemp)

				call Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats,ptree) 
				
				
					! if(ptree%MyID==2)write(*,*)fnorm(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyU%blocks(1)%matrix,size(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyU%blocks(1)%matrix,1),size(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyU%blocks(1)%matrix,2)),fnorm(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyV%blocks(1)%matrix,size(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyV%blocks(1)%matrix,1),size(ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1)%ButterflyV%blocks(1)%matrix,2)),'nddanaerfeffe',ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe,rowblock				
				
				
				
				! call Bplus_Sblock_randomized_symmetric(level_c,rowblock) 
				
				! call Butterfly_Sblock_randomized_symmetric(level_c,rowblock) 
				! call Butterfly_sym2asym(ho_bf1%levels(level_c)%matrices_block(rowblock))
				
				
				! call Butterfly_Sblock_randomized_logN(level_c,rowblock)
				
				! ! call Butterfly_Sblock_randomized(level_c,rowblock)
				! ! call Butterfly_Sblock_randomized_memfree_exact(level_c,rowblock)
				! ! stop
			! else
				! ! call Butterfly_Sblock_randomized_memfree(level_c,rowblock)
				! ! call Butterfly_Sblock_randomized(level_c,rowblock)
				! ! call Butterfly_Sblock_randomized_symmetric(level_c,rowblock)
			! end if
			
				call ComputeMemory_Bplus(ho_bf1%levels(level_c)%BP(rowblock),rtemp)
				stats%Mem_Sblock = stats%Mem_Sblock + rtemp
			
			
			! if(level_c==6)then			
				! call print_butterfly_size_rank(ho_bf1%levels(level_c)%matrices_block(rowblock),option%tol_comp)
				! stop
			! end if	
				

			
			! ! call Checkrank_Sblock(level_c,rowblock)
			! ! call Checkrank_Sblock_exact(level_c,rowblock)
			! if(level_c==ho_bf1%Maxlevel-2)stop
			
			
		end if	
		end do
		end do
		n2 = OMP_get_wtime()
		stats%Time_Sblock=stats%Time_Sblock+n2-n1 
	
		! write(*,*)'aha'
		! stop
		
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write(*,*)'compute block inverse at level:',level_c
		! compute the inverse butterfly
		n1 = OMP_get_wtime()
		! allocate(ho_bf1%levels(level_c)%matrices_block_inverse(ho_bf1%levels(level_c)%N_block_inverse))
		do rowblock = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		
			pgno =  ho_bf1%levels(level_c)%BP_inverse(rowblock)%pgno			
			if((ptree%MyID >=ptree%pgrp(pgno)%head .and. ptree%MyID <=ptree%pgrp(pgno)%tail))then	
				
				

				
				call DoubleDistributeBplus(ho_bf1%levels(level_c)%BP(rowblock*2-1),stats,ptree)
				call DoubleDistributeBplus(ho_bf1%levels(level_c)%BP(rowblock*2),stats,ptree)
				
				! write(*,*)level_c,rowblock,ptree%MyID,'aha'
				call Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,stats,ptree)
				! write(*,*)level_c,rowblock,ptree%MyID,'ahadone'
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
	! write (*,*) '     Time_Sblock_local:', stats%Time_Sblock	
	if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistB:', rtemp	
	call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistV:', rtemp		
	call MPI_ALLREDUCE(time_tmp,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*)'time_tmp',time_tmp
	! write (*,*)'time_resolve',time_resolve
	! write (*,*)'time_halfbuttermul',time_halfbuttermul
	call MPI_ALLREDUCE(stats%Flop_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write (*,'(A21Es14.2)') 'Factorization flops:',rtemp	
	
    if(ptree%MyID==Main_ID)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_SMW,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for butterfly inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for butterfly Sblocks'
	call MPI_ALLREDUCE(stats%Mem_Direct,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for direct inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_int_vec,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for storing intermidiate vectors'
    ! write(*,*)Memory_tfqmr_vec,'MB costed for storing intermidiate vectors in tfqmr'	
    if(ptree%MyID==Main_ID)write(*,*)''	

	
    return

end subroutine cascading_factorizing


end module cascading_factorization