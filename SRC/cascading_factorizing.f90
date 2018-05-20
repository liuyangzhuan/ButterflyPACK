module cascading_factorization
! use Butterfly_inversion_schur_iterative
use Butterfly_inversion_schur_partition
use Bplus_rightmultiply
use Bplus_inversion_schur_partition
contains 

subroutine cascading_factorizing(ho_bf1,option,stats)

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
	complex(kind=8), allocatable :: matrixtemp1(:,:)
	integer, allocatable :: ipiv(:)
	integer rowblock
	type(matrixblock),pointer::block_o,block_off,block_off1,block_off2
	type(matrixblock)::block_tmp
	real*8 n1,n2
	complex(kind=8),allocatable::Vin(:,:),Vout(:,:)
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	
	! Memory_int_vec = 0
	! Memory_tfqmr_vec = 0

    level_c=0
    flag=0

    write (*,*) ''
	
	call InitStat(stats)
	
	! Time_Init_forward=0
	! Time_Vector_forward=0
	! Time_Oneblock_forward=0					
	! Time_Reconstruct_forward=0
	! Time_Init_inverse=0
	! Time_Vector_inverse=0
	! Time_Reconstruct_inverse=0
	! Time_Oneblock_forward=0					
	! Time_InvertBlock = 0
	time_tmp = 0
	
	! ! call Butterfly_forward_acc_check
	! ! stop
	
    write (*,*) 'Computing block inverse at level Maxlevel_for_blocks+1...'	
	level_c = ho_bf1%Maxlevel+1
	! allocate(ho_bf1%levels(level_c)%matrices_block_inverse(ho_bf1%levels(level_c)%N_block_inverse))
	do ii = 1, 2**(level_c-1)

		nn = size(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,1)

		allocate(ipiv(nn))
	    call getrff90(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,ipiv)
        call getrif90(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,ipiv)		

		stats%Mem_Direct=stats%Mem_Direct+SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3		

		deallocate(ipiv)
	end do		


	
	write (*,*) 'Computing block inverse at higher levels...'
	do level_c = ho_bf1%Maxlevel,1,-1

		! update the forward butterfly after left multiplication of inverse 
		write(*,*)'update forward blocks at level:',level_c
					
		n1 = OMP_get_wtime()
		do rowblock = 1,2**level_c
			
			! if(level_c/=0)then	
			! if(level_c<ho_bf1%Maxlevel-8)then
				
				! call Butterfly_Sblock_randomized_memfree(level_c,rowblock,rtemp)
				! call Butterfly_Sblock_randomized_partialupdate_memfree(level_c,rowblock,rtemp)
				
				call Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats) 

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
				! call print_butterfly_size_rank(ho_bf1%levels(level_c)%matrices_block(rowblock),option%tol_SVD)
				! stop
			! end if	
				

			
			! ! call Checkrank_Sblock(level_c,rowblock)
			! ! call Checkrank_Sblock_exact(level_c,rowblock)
			! if(level_c==ho_bf1%Maxlevel-2)stop
			
			
			
		end do
		n2 = OMP_get_wtime()
		stats%Time_Sblock=stats%Time_Sblock+n2-n1 
	
		! write(*,*)'aha'
		! stop
		
		write(*,*)'compute block inverse at level:',level_c
		! compute the inverse butterfly
		n1 = OMP_get_wtime()
		! allocate(ho_bf1%levels(level_c)%matrices_block_inverse(ho_bf1%levels(level_c)%N_block_inverse))
		do rowblock = 1,2**(level_c-1)
			
			call Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,stats)
			call ComputeMemory_Bplus(ho_bf1%levels(level_c)%BP_inverse_schur(rowblock),rtemp)
			stats%Mem_SMW=stats%Mem_SMW+rtemp

		end do
		n2 = OMP_get_wtime()		
		stats%Time_SMW=stats%Time_SMW + n2-n1	
	end do
	
	
    write (*,*) 'computing updated forward block time:',stats%Time_Sblock,'Seconds'	
    write (*,*) 'computing inverse block time:',stats%Time_SMW,'Seconds'	
	write (*,*) '     Time_Init:', stats%Time_random(1)
	write (*,*) '     Time_MVP:', stats%Time_random(2)
	write (*,*) '     Time_Reconstruct:', stats%Time_random(3)
	write (*,*) '     Time_Onesub:', stats%Time_random(4)															 
	! write (*,*) '     Time_InvertBlock:', Time_InvertBlock
	write (*,*)'time_tmp',time_tmp
	! write (*,*)'time_resolve',time_resolve
	! write (*,*)'time_halfbuttermul',time_halfbuttermul
	
    write(*,*)''
    write(*,*)stats%Mem_SMW,'MB costed for butterfly inverse blocks'
    write(*,*)stats%Mem_Sblock,'MB costed for butterfly Sblocks'
    write(*,*)stats%Mem_Direct,'MB costed for direct inverse blocks'
    write(*,*)stats%Mem_int_vec,'MB costed for storing intermidiate vectors'
    ! write(*,*)Memory_tfqmr_vec,'MB costed for storing intermidiate vectors in tfqmr'	
    write(*,*)''	

	
    return

end subroutine cascading_factorizing


end module cascading_factorization