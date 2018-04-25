module cascading_factorization
! use Butterfly_inversion_schur_iterative
use Butterfly_inversion_schur_partition
use Bplus_rightmultiply
use Bplus_inversion_schur_partition
contains 

subroutine cascading_factorizing(tolerance)

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
	real*8 tolerance, rtemp,tmpfact
    real*8 Memory, Memory_near,Memory_butterfly_inverse,Memory_direct_inverse,Memory_butterfly_Sblock
	integer,allocatable:: index_old(:),index_new(:) 
	integer::block_num,block_num_new,num_blocks,level_butterfly	
	complex(kind=8), allocatable :: matrixtemp1(:,:)
	integer, allocatable :: ipiv(:)
	integer rowblock
	type(matrixblock),pointer::block_o,block_off,block_off1,block_off2
	type(matrixblock)::block_tmp
	real*8 n1,n2
	complex(kind=8),allocatable::Vin(:,:),Vout(:,:)
	
	
    Memory_direct_inverse=0
    Memory_butterfly_inverse=0
    Memory_butterfly_Sblock=0
	Memory_int_vec = 0
	Memory_tfqmr_vec = 0

    level_c=0
    flag=0

    write (*,*) ''
	Time_Rightmultiply_inverse_randomized = 0
	Time_Inversion_diagnal_randomized = 0
	Time_Init_forward=0
	Time_Vector_forward=0
	Time_Oneblock_forward=0					
	Time_Reconstruct_forward=0
	Time_Init_inverse=0
	Time_Vector_inverse=0
	Time_Reconstruct_inverse=0
	Time_Oneblock_forward=0					
	Time_InvertBlock = 0
	time_tmp = 0
	
	! ! call Butterfly_forward_acc_check
	! ! stop
	
	
	
	
    write (*,*) 'Computing block inverse at level Maxlevel_for_blocks+1...'	
	level_c = Maxlevel_for_blocks+1
	! allocate(cascading_factors(level_c)%matrices_block_inverse(cascading_factors(level_c)%N_block_inverse))
	do ii = 1, 2**(level_c-1)
		! cascading_factors(level_c)%matrices_block_inverse(ii)%level = cascading_factors(level_c)%matrices_block(ii)%level
		! cascading_factors(level_c)%matrices_block_inverse(ii)%col_group = cascading_factors(level_c)%matrices_block(ii)%col_group
		! cascading_factors(level_c)%matrices_block_inverse(ii)%row_group = cascading_factors(level_c)%matrices_block(ii)%row_group
		! cascading_factors(level_c)%matrices_block_inverse(ii)%nested_num = cascading_factors(level_c)%matrices_block(ii)%nested_num
		
		! cascading_factors(level_c)%matrices_block_inverse(ii)%style = cascading_factors(level_c)%matrices_block(ii)%style
		! cascading_factors(level_c)%matrices_block_inverse(ii)%data_type = cascading_factors(level_c)%matrices_block(ii)%data_type
		nn = size(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,1)
		allocate(cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat(nn,nn))
		allocate(matrixtemp1(nn,nn))
		allocate(ipiv(nn))
		matrixtemp1 = cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat
	    call getrff90(matrixtemp1,ipiv)
		! write(*,*)shape(matrixtemp1)
        call getrif90(matrixtemp1,ipiv)		
		cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat = matrixtemp1
		
		Memory_direct_inverse=Memory_direct_inverse+SIZEOF(cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3		
		
		call delete_blocks(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1))
		
		deallocate(matrixtemp1)
		deallocate(ipiv)
	end do		


	
	write (*,*) 'Computing block inverse at higher levels...'
	do level_c = Maxlevel_for_blocks,1,-1

		! update the forward butterfly after left multiplication of inverse 
		write(*,*)'update forward blocks at level:',level_c
		tmpfact = Rank_detection_factor
		! Rank_detection_factor = Rank_detection_factor * 0.1														
		n1 = OMP_get_wtime()
		do rowblock = 1,2**level_c
			
			adaptive_flag=1   ! use adaptive rank in S block construction
			if(level_c/=0)then	
			! if(level_c<Maxlevel_for_blocks-8)then
				
				! call Butterfly_Sblock_randomized_memfree(level_c,rowblock,rtemp)
				! call Butterfly_Sblock_randomized_partialupdate_memfree(level_c,rowblock,rtemp)
				
				call Bplus_Sblock_randomized_memfree(level_c,rowblock,rtemp) 
				
				! call Bplus_Sblock_randomized_symmetric(level_c,rowblock) 
				
				! call Butterfly_Sblock_randomized_symmetric(level_c,rowblock) 
				! call Butterfly_sym2asym(cascading_factors(level_c)%matrices_block(rowblock))
				
				
				! call Butterfly_Sblock_randomized_logN(level_c,rowblock)
				
				! ! call Butterfly_Sblock_randomized(level_c,rowblock)
				! ! call Butterfly_Sblock_randomized_memfree_exact(level_c,rowblock)
				! ! stop
			else
				! call Butterfly_Sblock_randomized_memfree(level_c,rowblock)
				! call Butterfly_Sblock_randomized(level_c,rowblock)
				! call Butterfly_Sblock_randomized_symmetric(level_c,rowblock)
			end if
			
			
			! if(level_c==6)then			
				! call print_butterfly_size_rank(cascading_factors(level_c)%matrices_block(rowblock))
				! stop
			! end if	
				

			
			! ! call Checkrank_Sblock(level_c,rowblock)
			! ! call Checkrank_Sblock_exact(level_c,rowblock)
			! if(level_c==Maxlevel_for_blocks-2)stop
			
			Memory_butterfly_Sblock = Memory_butterfly_Sblock + rtemp
			
		end do
		n2 = OMP_get_wtime()
		Time_Rightmultiply_inverse_randomized=Time_Rightmultiply_inverse_randomized+n2-n1 
		Rank_detection_factor = tmpfact

		! write(*,*)'aha'
		! stop
		
		write(*,*)'compute block inverse at level:',level_c
		! compute the inverse butterfly
		n1 = OMP_get_wtime()
		! allocate(cascading_factors(level_c)%matrices_block_inverse(cascading_factors(level_c)%N_block_inverse))
		do rowblock = 1,2**(level_c-1)
			
			if(schurinv==0)then
				write(*,*)'schurinv=0 removed'
				stop
			else 
				! ! ! block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
				! ! block_off => cascading_factors(level_c)%matrices_block(rowblock*2-1)			
				! ! ! block_o%level_butterfly = block_off%level_butterfly+1
				! ! block_o => cascading_factors(level_c)%matrices_block_inverse_schur(rowblock)
				! ! block_o%level_butterfly = block_off%level_butterfly
				

				
				if(directschur==1)then
					! use fixed rank in block inverse
					adaptive_flag=1
					tmpfact = Rank_detection_factor
					! Rank_detection_factor = Rank_detection_factor * 0.001						
					! call Butterfly_inverse_schur_partitionedinverse(level_c,rowblock,rtemp)
					call Bplus_inverse_schur_partitionedinverse(level_c,rowblock,rtemp)
					Rank_detection_factor = tmpfact
				else				
					write(*,*)'directshur==0 removed'
				end if

				
			end if
			
			Memory_butterfly_inverse=Memory_butterfly_inverse+rtemp
			
	
		end do
		n2 = OMP_get_wtime()		
		Time_Inversion_diagnal_randomized=Time_Inversion_diagnal_randomized + n2-n1	
	end do
	
	
    write (*,*) 'computing updated forward block time:',Time_Rightmultiply_inverse_randomized,'Seconds'	
	write (*,*) '     Time_Init_forward:', Time_Init_forward
	write (*,*) '     Time_Vector_forward:', Time_Vector_forward
	write (*,*) '     Time_Reconstruct_forward:', Time_Reconstruct_forward
	write (*,*) '     Time_Oneblock_forward:', Time_Oneblock_forward															 
	
    write (*,*) 'computing inverse block time:',Time_Inversion_diagnal_randomized,'Seconds'	
	write (*,*) '     Time_Init_inverse:', Time_Init_inverse
	write (*,*) '     Time_Vector_inverse:', Time_Vector_inverse
	write (*,*) '     Time_Reconstruct_inverse:', Time_Reconstruct_inverse
	write (*,*) '     Time_Oneblock_inverse:', Time_Oneblock_inverse															 
	write (*,*) '     Time_InvertBlock:', Time_InvertBlock
	write (*,*)'time_tmp',time_tmp
	write (*,*)'time_resolve',time_resolve
	write (*,*)'time_halfbuttermul',time_halfbuttermul
	
    write(*,*)''
    write(*,*)Memory_butterfly_inverse,'MB costed for butterfly inverse blocks'
    write(*,*)Memory_butterfly_Sblock,'MB costed for butterfly Sblocks'
    write(*,*)Memory_direct_inverse,'MB costed for direct inverse blocks'
    write(*,*)Memory_int_vec,'MB costed for storing intermidiate vectors'
    write(*,*)Memory_tfqmr_vec,'MB costed for storing intermidiate vectors in tfqmr'	
    write(*,*)''	

	
    return

end subroutine cascading_factorizing


! ! ! subroutine testfactorization
	
    ! ! ! use MODULE_FILE
    ! ! ! ! use lapack95
    ! ! ! implicit none
    
	! ! ! integer level_c,rowblock
    ! ! ! integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    ! ! ! integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    ! ! ! character chara
    ! ! ! real*8 a,b,c,d
    ! ! ! complex(kind=8) ctemp, ctemp1, ctemp2
	! ! ! type(matrixblock),pointer::block_o
	
    ! ! ! type(vectorsblock), pointer :: random1, random2
    
    ! ! ! real*8,allocatable :: Singular(:)
	! ! ! integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	! ! ! complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
   
    
	! ! ! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	
	
    ! ! ! level_butterfly=maxlevel_for_blocks-block_o%level
    ! ! ! num_blocks=2**level_butterfly
    ! ! ! allocate (RandomVectors_InOutput(6))

    ! ! ! groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    ! ! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
    
	! ! ! ! num_vectors=(level_butterfly+1)*(int(nn/2**level_butterfly)+1)+10
    ! ! ! num_vectors=int(8*butterfly_block_randomized(1)%dimension_rank)     !!!!!!!!!!!!! tune this 
    
	
    ! ! ! groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    ! ! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	! ! ! allocate (RandomVectors_InOutput(1)%vector(nn,num_vectors))
	! ! ! allocate (RandomVectors_InOutput(6)%vector(nn,num_vectors))
    
    ! ! ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    ! ! ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    ! ! ! allocate (RandomVectors_InOutput(2)%vector(mm,num_vectors))
    ! ! ! allocate (RandomVectors_InOutput(3)%vector(mm,num_vectors))    
	! ! ! allocate (RandomVectors_InOutput(4)%vector(mm,num_vectors))
    ! ! ! allocate (RandomVectors_InOutput(5)%vector(mm,num_vectors))
	! ! ! do ii =1,6
		! ! ! RandomVectors_InOutput(ii)%vector = 0
	! ! ! end do
	 
	! ! ! ! write(*,*)groupn,basis_group(groupn)%head,groupm,basis_group(groupm)%head
	 
	 
	! ! ! idx_start_glo = basis_group(groupm)%head	 
	   
    ! ! ! !$omp parallel do default(shared) private(ii,jj,a,b,c,d)
    ! ! ! do jj=1, num_vectors
        ! ! ! do ii=1, nn
            ! ! ! call random_number(a)
            ! ! ! call random_number(b)
            ! ! ! call random_number(c)
            ! ! ! call random_number(d)
            ! ! ! if (c<0.5) then
                ! ! ! a=-a
            ! ! ! endif
            ! ! ! if (d<0.5) then
                ! ! ! b=-b
            ! ! ! endif
            ! ! ! RandomVectors_InOutput(1)%vector(ii,jj)=a+junit*b
        ! ! ! enddo
        ! ! ! do ii=1, mm
            ! ! ! call random_number(a)
            ! ! ! call random_number(b)
            ! ! ! call random_number(c)
            ! ! ! call random_number(d)
            ! ! ! if (c<0.5) then
                ! ! ! a=-a
            ! ! ! endif
            ! ! ! if (d<0.5) then
                ! ! ! b=-b
            ! ! ! endif
            ! ! ! RandomVectors_InOutput(4)%vector(ii,jj)=a+junit*b
        ! ! ! enddo
    ! ! ! enddo
    ! ! ! !$omp end parallel do
    
	
	
	! ! ! ! get the right multiplied vectors
	
    ! ! ! random1=>RandomVectors_InOutput(1)
    ! ! ! random2=>RandomVectors_InOutput(2)
    ! ! ! ctemp1=1.0 ; ctemp2=0.0
    ! ! ! call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	
	! ! ! allocate(vec_old(mm,num_vectors))
	! ! ! allocate(vec_new(mm,num_vectors))
	! ! ! vec_old = RandomVectors_InOutput(2)%vector

	! ! ! ! write(*,*)'begin'
	
	! ! ! do level = Maxlevel_for_blocks+1,level_c+1,-1
		! ! ! N_diag = 2**(level-level_c-1)
		! ! ! idx_start_diag = (rowblock-1)*N_diag+1
		! ! ! vec_new = 0
		! ! ! do ii = idx_start_diag,idx_start_diag+N_diag-1
			! ! ! ! write(*,*)level,ii
			! ! ! groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			! ! ! idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			! ! ! idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			
			! ! ! ! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
			
			! ! ! ! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
			! ! ! if(level==Maxlevel_for_blocks+1)then
				! ! ! call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
				! ! ! &vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! ! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
			! ! ! else 
				! ! ! call butterfly_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,idx_end_loc-idx_start_loc+1,num_vectors,&
				! ! ! &vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors) + vec_new(idx_start_loc:idx_end_loc,1:num_vectors)
				! ! ! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
			! ! ! end if
		! ! ! end do
		! ! ! vec_old = vec_new
	! ! ! end do
	! ! ! ! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	! ! ! RandomVectors_InOutput(3)%vector = vec_new
	! ! ! deallocate(vec_old)
	! ! ! deallocate(vec_new)
	
	! ! ! ! ! ! ! ! write(*,*)'end'
	
     
    ! ! ! return                

! ! ! end subroutine testfactorization




end module cascading_factorization