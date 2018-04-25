module Butterfly_inversion
use Utilites_randomized
! use Butterfly_compression_givenfullmat
use Bplus_rightmultiply

integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
complex(kind=dp),allocatable::r0_initial(:)

contains 

subroutine Butterfly_inverse_diagonal_randomized(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax
	real*8:: n1,n2,Memory
		
    block_o =>  cascading_factors(level_c)%matrices_block_inverse(rowblock) 
	Memory = 0
	
	do tt = 1,1
		!  write(*,*)'10'
		n1 = OMP_get_wtime()
		call Initialize_Butterfly_inverse(level_c,rowblock,tt-1)
		n2 = OMP_get_wtime()
		Time_Init_inverse = Time_Init_inverse + n2-n1
		!  write(*,*)'11'
		n1 = OMP_get_wtime()
		itermax = 0	
		call Get_Randomized_Vectors_inverse(level_c,rowblock,itermax)        
		n2 = OMP_get_wtime()
		Time_Vector_inverse = Time_Vector_inverse + n2-n1
		!  write(*,*)'12'
		n1 = OMP_get_wtime()
		call Resolving_Butterfly_v8(niter,error_inout,block_o)
		n2 = OMP_get_wtime()	
		Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
		!  write(*,*)'13'
		if(error_inout>iter_tolerance)then
			call Delete_randomized_butterfly()
		else 
			call delete_blocks(block_o)
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			rank_new_max = butterfly_block_randomized(1)%rankmax
			call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
			deallocate(butterfly_block_randomized)
			! call copy_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
			! call Delete_randomized_butterfly()
				
			write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7,A7,I3)')'     No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',block_o%level_butterfly,' error:',error_inout,' Niter:',itermax
			
			return			
		end if		
	end do
	
	write(*,*)'randomized scheme not converged gaga'
	stop

    return

end subroutine Butterfly_inverse_diagonal_randomized


subroutine Butterfly_inverse_diagonal_randomized_memfree(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
		
    block_o =>  cascading_factors(level_c)%matrices_block_inverse(rowblock) 
	Memory = 0
	
	do tt = 1,10
		do ntry=1,1
		itermax = 0
		n1 = OMP_get_wtime()
		call Initialize_Butterfly_inverse(level_c,rowblock,tt-1)
		n2 = OMP_get_wtime()
		Time_Init_inverse = Time_Init_inverse + n2-n1
		
		n1 = OMP_get_wtime()
		call Reconstruction_LL_inverse(level_c,rowblock,itermax)	
		call Reconstruction_RR_inverse(level_c,rowblock,error_inout,itermax)
		n2 = OMP_get_wtime()	
		Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
		
		! write(*,*)tt,error_inout
		
		
		if(error_inout>iter_tolerance)then
		! if(0)then
			call Delete_randomized_butterfly()
		else 
			call delete_blocks(block_o)
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			rank_new_max = butterfly_block_randomized(1)%rankmax
			call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
			deallocate(butterfly_block_randomized)
			! call copy_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
			! call Delete_randomized_butterfly()
			write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7,A7,I3)')'     No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',block_o%level_butterfly,' error:',error_inout,' Niter:',itermax
				
			
			return			
		end if		
		end do
	end do
	
	write(*,*)'randomized scheme not converged gan'
	stop

    return

end subroutine Butterfly_inverse_diagonal_randomized_memfree



subroutine Butterfly_inverse_diagonal_randomized_symmetric(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax
	real*8:: n1,n2,Memory
		
    block_o =>  cascading_factors(level_c)%matrices_block_inverse(rowblock) 
	block_o%level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	level_butterfly = int((maxlevel_for_blocks-block_o%level)/2)*2
	
	Memory = 0
	if(level_butterfly<2)then
		call Butterfly_inverse_diagonal_randomized(level_c,rowblock,Memory)
	else 
		do tt = 1,1
			!  write(*,*)'10'
			! n1 = OMP_get_wtime()
			! call Initialize_Butterfly_inverse(level_c,rowblock,tt-1)
			! n2 = OMP_get_wtime()
			! Time_Init_inverse = Time_Init_inverse + n2-n1
			!  write(*,*)'11'
			n1 = OMP_get_wtime()
			itermax = 0	
			call Get_Randomized_Vectors_inverse_Memeff(level_c,rowblock,tt-1,itermax)        
			n2 = OMP_get_wtime()
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			!  write(*,*)'12'
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_Symmetric(block_o%row_group,block_o%col_group)
			call Test_Error_RR_inverse(level_c,rowblock,error_inout)			
			n2 = OMP_get_wtime()	
			Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
			
			!  write(*,*)'13'
			if(error_inout>iter_tolerance)then
				call Delete_randomized_butterfly()
			else 
				call delete_blocks(block_o)
				call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
				rank_new_max = butterfly_block_randomized(1)%rankmax
				call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
				deallocate(butterfly_block_randomized)
				! call copy_randomizedbutterfly(butterfly_block_randomized(1),block_o,Memory)
				! call Delete_randomized_butterfly()
					
				write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7,A7,I3)')'     No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',block_o%level_butterfly,' error:',error_inout,' Niter:',itermax
				
				return			
			end if		
		end do
		
		write(*,*)'randomized scheme not converged',error_inout
		stop	
	end if

    return

end subroutine Butterfly_inverse_diagonal_randomized_symmetric


subroutine Reconstruction_LL_inverse(level_c,rowblock,itermax)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock),pointer::block_o,block_off1,block_off2
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,itermax,unique_nth
	real*8:: error_inout
	integer,allocatable::perms(:)

	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
	level_butterfly=butterfly_block_randomized(1)%level_butterfly
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
    allocate (Random_Block(1))
    
    allocate (Random_Block(1)%RandomVectorLL(0:level_butterfly+2))    
	Nbind = 1
	
	num_vect_sub = num_vect_subsub*Nbind
	
    random=>Random_Block(1)
	call Init_RandVect_Empty('T',random,num_vect_sub)	
	
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	! ! level_right_start = level_butterfly+1
	
	! call Zero_Butterfly(0,level_right_start)
 
    ! allocate(perms(Nsub))
	! call rperm(Nsub, perms)
	
	
	do unique_nth = 0,level_right_start
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
		Ng = 2**level_butterfly/Nsub	
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
		! do ii = 1,Nsub		
			! nth_s = perms(ii)
			! nth_e = perms(ii)
			
			n1 = OMP_get_wtime()
			call Get_Randomized_Vectors_LL_inverse(level_c,rowblock,nth_s,nth_e,num_vect_sub,itermax,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1		
		end do
	end do
	
	! deallocate(perms)
	
	
	random=>Random_Block(1)
	call Delete_RandVect('T',random,level_butterfly)

    return
    
end subroutine Reconstruction_LL_inverse




subroutine Reconstruction_RR_inverse(level_c,rowblock,error,itermax)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,rowblock    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_left_start,num_row,num_col
    real*8::n1,n2,error
	integer itermax
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock),pointer::block_o,block_off1,block_off2
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	
	level_butterfly=butterfly_block_randomized(1)%level_butterfly
		
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
	! ! call assert(num_vectors==2**level_butterfly*dimension_rank,'incorrect num_vectors')
	
	! call assert(num_vectors==Nsub*dimension_rank,'incorrect num_vectors') !  check here later
	
    ! ! allocate (Random_Block(1))   !  check here later 
    
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    
	
	Nbind = 1
	
	num_vect_sub = num_vect_subsub*Nbind
	
    random=>Random_Block(1)
	call Init_RandVect_Empty('N',random,num_vect_sub)

    level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    ! level_left_start = 0
	
	! call Zero_Butterfly(level_left_start,level_butterfly+1)

	do unique_nth=level_butterfly+1,level_left_start,-1
		if(mod(level_butterfly,2)==0)then
			Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
		else 
			Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
		end if	
		Ng = 2**level_butterfly/Nsub	
		
		do ii = 1,Nsub/Nbind
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			n1 = OMP_get_wtime()
			call Get_Randomized_Vectors_RR_inverse(level_c,rowblock,nth_s,nth_e,num_vect_sub,itermax,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_inverse(level_c,rowblock,error)


	
    return
    
end subroutine Reconstruction_RR_inverse

subroutine Test_Error_RR_inverse(level_c,rowblock,error)

    use MODULE_FILE
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm
    integer mm,nn
    real*8 a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    complex(kind=8) ctemp
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real*8::error
	integer level_c,rowblock,dimension_m 
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock),pointer::block_o
	
	
	block_o =>  cascading_factors(level_c)%matrices_block_inverse(rowblock) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	num_blocks=2**level_butterfly
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call Get_Randomized_Vectors_RR_Test_inverse(level_c,rowblock,num_vect)

	k=0
	do i=1, num_blocks
		dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, dimension_m
			do jj=1, num_vect
				RandomVectors_Output_ref(ii+k,jj)=random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+dimension_m
	enddo 	

	call Butterfly_partial_MVP('N',0,level_butterfly+1,random)

	k=0
	norm3_R=0 ; norm4_R=0
	do i=1, num_blocks
		 dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		 norm1_R=0 ; norm2_R=0
		 do ii=1, dimension_m
			do jj =1,num_vect
				 norm1_R=norm1_R+abs(random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj))**2
				 norm2_R=norm2_R+abs(random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)-RandomVectors_Output_ref(ii+k,jj))**2
			enddo
 		 enddo
		 norm3_R=norm3_R+norm1_R
		 norm4_R=norm4_R+norm2_R
		 k=k+dimension_m
	enddo 
	error = sqrt(norm4_R/norm3_R)
	
	
	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	deallocate(Random_Block)

	deallocate(RandomVectors_Output_ref)
	
    return                

end subroutine Test_Error_RR_inverse


subroutine Initialize_Butterfly_inverse(level_c,rowblock,kover)

    use MODULE_FILE
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,index_j,index_i
    real*8 a,b,c,d
    complex (kind=8) ctemp
	type(matrixblock),pointer::block_o,block_off1,block_off2
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	! write(*,*)level_c+1,rowblock*2-1
	
	allocate (butterfly_block_randomized(1))
    
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    
    num_blocks=2**level_butterfly
 	
	dimension_rank= max(block_off1%rankmax,block_off2%rankmax)+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	
	
	! if(level_c==2)dimension_rank=11
	if(level_c==1)dimension_rank=9+kover
	
	! write(*,*)dimension_rank
	
    groupm=block_o%row_group   ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    butterfly_block_randomized(1)%dimension_rank=dimension_rank
    !write (*,*) dimension_rank, int(real(kk)/num_blocks/2.0), mm
    
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))
    ! allocate (butterfly_block_randomized(1)%ButterflyU_old(2**level_butterfly))
    ! allocate (butterfly_block_randomized(1)%ButterflyV_old(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyUInv(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyVInv(2**level_butterfly))
    

	dimension_max = 2*dimension_rank
	do blocks=1, num_blocks	
		dimension_m=size(block_off1%ButterflyU(blocks)%matrix,1)
		dimension_n=size(block_off1%ButterflyV(blocks)%matrix,1)
		dimension_m=size(block_off2%ButterflyU(blocks)%matrix,1)
		dimension_n=size(block_off2%ButterflyV(blocks)%matrix,1)		
		dimension_max = max(dimension_max,dimension_m)	
		dimension_max = max(dimension_max,dimension_n)	
	end do	
	allocate(butterfly_block_randomized(1)%KerInv(dimension_max,dimension_max))
	call RandomMat(dimension_max,dimension_max,dimension_max,butterfly_block_randomized(1)%KerInv,3)	
	
	
    do blocks=1, num_blocks
		
		if(blocks<=num_blocks/2)then
			dimension_m=size(block_off1%ButterflyU(blocks)%matrix,1)
		else
			dimension_m=size(block_off2%ButterflyU(blocks-num_blocks/2)%matrix,1)
		end if		
			
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyUInv(blocks)%matrix(dimension_rank,dimension_m))
        
		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        
		!$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		!$omp end parallel do
		
		! butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyU(blocks)%matrix
		deallocate(matrixtemp1)		
		
		if(blocks<=num_blocks/2)then
			dimension_n=size(block_off2%ButterflyV(blocks)%matrix,1)
		else
			dimension_n=size(block_off1%ButterflyV(blocks-num_blocks/2)%matrix,1)
		end if	

        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        
		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
		
		!$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		!$omp end parallel do
		
		! butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyV(blocks)%matrix
		deallocate(matrixtemp1)		
		
    enddo
    ! call invert_Butterfly_U()
    ! call invert_Butterfly_V()
    
    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
        allocate (butterfly_block_randomized(1)%ButterflyKerl(level_butterfly))
        ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level_butterfly))
        allocate (butterfly_block_randomized(1)%ButterflyInv(level_butterfly))
        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%num_row=num_row
            ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))
            ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(num_row,num_col))
            allocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(num_row/2,num_col/2))
            ! do j=1, num_col, 2
                ! index_j=int((j+1)/2)
                ! do i=1, num_row, 2
                    ! index_i=int((i+1)/2)
                    ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(dimension_rank,dimension_rank))
                    ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i+1,j)%matrix(dimension_rank,dimension_rank))
                    ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j+1)%matrix(dimension_rank,dimension_rank))
                    ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(dimension_rank,dimension_rank))
                    ! ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i,j)%matrix(dimension_rank,dimension_rank))
                    ! ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i+1,j)%matrix(dimension_rank,dimension_rank))
                    ! ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i,j+1)%matrix(dimension_rank,dimension_rank))
                    ! ! ! ! allocate (butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i+1,j+1)%matrix(dimension_rank,dimension_rank))
                   
				    ! ! ! call RandomMat(2*dimension_rank,2*dimension_rank,2*dimension_rank,matrixtemp1,0)	
					
                    ! ! ! !$omp parallel do default(shared) private(ii,jj)
                    ! ! ! do jj=1, dimension_rank
                        ! ! ! do ii=1, dimension_rank
                            ! ! ! butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=matrixtemp1(ii,jj)
                            ! ! ! ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i,j)%matrix(ii,jj)=matrixtemp1(ii,jj)
                            ! ! ! butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i+1,j)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj)
                            ! ! ! ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i+1,j)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj)
                            ! ! ! butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,jj)=matrixtemp1(ii,jj+dimension_rank)
                            ! ! ! ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i,j+1)%matrix(ii,jj)=matrixtemp1(ii,jj+dimension_rank)
                            ! ! ! butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj+dimension_rank)
                            ! ! ! ! butterfly_block_randomized(1)%ButterflyKerl_old(level)%blocks(i+1,j+1)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj+dimension_rank)
                        ! ! ! enddo
                    ! ! ! enddo
                    ! ! ! !$omp end parallel do
                    ! allocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(index_i,index_j)%matrix(2*dimension_rank,2*dimension_rank))
                ! enddo
            ! enddo
            ! call invert_Butterfly_Kernel(level)
        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine Initialize_Butterfly_inverse


subroutine Get_Randomized_Vectors_LL_inverse(level_c,rowblock,nth_s,nth_e,num_vect_sub,itermax,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    integer itermax
	character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_right_start
	type(RandomBlock), pointer :: random
	
	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 


    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i,header_m,tailer_m,mm,k)	
		do i=(nth-1)*Ng+1, nth*Ng	
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_m=basis_group(groupm_start+i-1)%head
				tailer_m=basis_group(groupm_start+i-1)%tail
				mm=tailer_m-header_m+1
				k=header_m-header_mm	

				! allocate(matrixtemp1(num_vect_subsub,mm))
				call RandomMat(mm,num_vect_subsub,min(mm,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:mm+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
				
				
				
				
				
				! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) !random_complex_number()	
					 ! enddo
				 ! enddo
				 ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	call GetInputVector_Inverse_tfqmr_batch_LL(level_c,rowblock,num_vect_sub,mm,RandomVectors_InOutput(1)%vector&
	&,RandomVectors_InOutput(2)%vector,itermax)		
		
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector - RandomVectors_InOutput(1)%vector
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+mm
	enddo 
	
	k=0
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+nn
	enddo 	

    !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_LL_inverse



subroutine Get_Randomized_Vectors_RR_inverse(level_c,rowblock,nth_s,nth_e,num_vect_sub,itermax,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks,itermax
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	type(RandomBlock), pointer :: random
	
	
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub	
		
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i,header_m,tailer_m,mm,k)	
		do i=(nth-1)*Ng+1, nth*Ng	
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_m=basis_group(groupm_start+i-1)%head
				tailer_m=basis_group(groupm_start+i-1)%tail
				mm=tailer_m-header_m+1
				k=header_m-header_mm	

				! allocate(matrixtemp1(num_vect_subsub,mm))
				call RandomMat(mm,num_vect_subsub,min(mm,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:mm+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
				
				
				
				
				! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! random_complex_number()	
					 ! enddo
				 ! enddo
				 ! !$omp end parallel do
				! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	call GetInputVector_Inverse_tfqmr_batch_RR(level_c,rowblock,num_vect_sub,mm,RandomVectors_InOutput(1)%vector&
	&,RandomVectors_InOutput(2)%vector,itermax)		
		
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector - RandomVectors_InOutput(1)%vector
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+mm
	enddo 
	
    !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_inverse




subroutine Get_Randomized_Vectors_RR_Test_inverse(level_c,rowblock,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks,itermax
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer num_vect_sub
	type(RandomBlock), pointer :: random
	
	
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	
	do i=1, num_blocks			
		header_m=basis_group(groupm_start+i-1)%head
		tailer_m=basis_group(groupm_start+i-1)%tail
		mm=tailer_m-header_m+1
		k=header_m-header_mm	
		!$omp parallel do default(shared) private(ii,jj)
		 do jj=1,num_vect_sub
			 do ii=1, mm
				 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	
			 enddo
		 enddo
		 !$omp end parallel do
	end do

	! get the right multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	call GetInputVector_Inverse_tfqmr_batch_RR(level_c,rowblock,num_vect_sub,mm,RandomVectors_InOutput(1)%vector&
	&,RandomVectors_InOutput(2)%vector,itermax)		
		
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector - RandomVectors_InOutput(1)%vector
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+mm
	enddo 
	
    !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_Test_inverse


subroutine Get_Randomized_Vectors_inverse(level_c,rowblock,itermax)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	integer itermax


   
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(6))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)
    Ng = 2**level_butterfly/Nsub
	dimension_rank = butterfly_block_randomized(1)%dimension_rank
	num_vectors = Nsub*dimension_rank	

	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(3)%vector(nn,num_vectors))    
	allocate (RandomVectors_InOutput(4)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(5)%vector(nn,num_vectors))
	allocate (RandomVectors_InOutput(6)%vector(nn,num_vectors))	 
	do ii =1,6
		RandomVectors_InOutput(ii)%vector = 0
	end do

	idx_start_glo = basis_group(groupn)%head		
	level_blocks = block_o%level



	groupn_start=groupn*2**(level_butterfly)
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	do i=1, num_blocks
		header_n=basis_group(groupn_start+i-1)%head
		tailer_n=basis_group(groupn_start+i-1)%tail
		nn=tailer_n-header_n+1
		k=header_n-header_nn
		
		
		 !$omp parallel do default(shared) private(ii,jj)
		 do jj=idx_start,idx_start+dimension_rank-1
			 do ii=1, nn
				 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()
				 RandomVectors_InOutput(4)%vector(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo 	 
	   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	if(level_c==6)then
	call GetInputVector_Inverse_direct(level_c,rowblock,num_vectors,nn,RandomVectors_InOutput(1)%vector&
	&,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(5)%vector)	
	end if
	
	! ! call GetInputVector_Inverse_tfqmr(level_c,rowblock,num_vectors,nn,RandomVectors_InOutput(1)%vector&
	! ! &,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(5)%vector,itermax)		
	!  write(*,*)'heihei'
	call GetInputVector_Inverse_tfqmr_batch(level_c,rowblock,num_vectors,nn,RandomVectors_InOutput(1)%vector&
	&,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(5)%vector,itermax)		
		
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector - RandomVectors_InOutput(1)%vector
	RandomVectors_InOutput(6)%vector = RandomVectors_InOutput(5)%vector - RandomVectors_InOutput(4)%vector	
		!  write(*,*)'heihei111'
    return                

end subroutine Get_Randomized_Vectors_inverse



subroutine Get_Randomized_Vectors_inverse_Memeff(level_c,rowblock,kover,itermax)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d,mem_vec
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng,num_vectorsR,num_vectorsL
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	integer itermax,kover,DimMax,levelm
	integer Nbatch,batch_size,ll,idxs,idxe

   
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	! dimension_rank= max(block_off1%rankmax,block_off2%rankmax)+3+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= max(block_off1%rankmax,block_off2%rankmax)*2+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
    
	allocate (butterfly_block_randomized(1))
	butterfly_block_randomized(1)%dimension_rank=dimension_rank
	butterfly_block_randomized(1)%level_butterfly = level_butterfly    
	levelm = ceiling_safe(dble(level_butterfly)/2d0)
	num_blocks=2**level_butterfly
    
! first consider the right multiplied vectors 
	Nsub = 2**(level_butterfly-levelm)
    Ng = 2**level_butterfly/Nsub
	num_vectorsR=Nsub*dimension_rank
	
	groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	DimMax = nn
	
	allocate(RandVectInR(DimMax,num_vectorsR))
	allocate(RandVectOutR(DimMax,num_vectorsR))
	RandVectInR=0
	groupn_start=groupn*2**level_butterfly
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	do i=1, num_blocks
		header_n=basis_group(groupn_start+i-1)%head
		tailer_n=basis_group(groupn_start+i-1)%tail
		nn=tailer_n-header_n+1
		k=header_n-header_nn
		 !$omp parallel do default(shared) private(ii,jj)
		 do jj=idx_start,idx_start+dimension_rank-1
			 do ii=1, nn
				 RandVectInR(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo 	
	
	batch_size = 10
	Nbatch = ceiling_safe(num_vectorsR/dble(batch_size))
	do ll=1,Nbatch
		idxs = (ll-1)*batch_size+1
		idxe = min(ll*batch_size,num_vectorsR)
		call GetInputVector_Inverse_tfqmr_batch_RR(level_c,rowblock,idxe-idxs+1,DimMax,RandVectInR(1:DimMax,idxs:idxe),RandVectOutR(1:DimMax,idxs:idxe),itermax)		
	end do

	RandVectOutR = RandVectOutR - RandVectInR
	
	
	
! next consider the left multiplied vectors 
	Nsub = 2**levelm
    Ng = 2**level_butterfly/Nsub
	num_vectorsL=Nsub*dimension_rank
	
	groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	DimMax = nn
	
	allocate(RandVectInL(DimMax,num_vectorsL))
	allocate(RandVectOutL(DimMax,num_vectorsL))
	RandVectInL=0
	groupn_start=groupn*2**level_butterfly
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	do i=1, num_blocks
		header_n=basis_group(groupn_start+i-1)%head
		tailer_n=basis_group(groupn_start+i-1)%tail
		nn=tailer_n-header_n+1
		k=header_n-header_nn
		 !$omp parallel do default(shared) private(ii,jj)
		 do jj=idx_start,idx_start+dimension_rank-1
			 do ii=1, nn
				 RandVectInL(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo 	

	batch_size = 10
	Nbatch = ceiling_safe(num_vectorsL/dble(batch_size))
	do ll=1,Nbatch
		idxs = (ll-1)*batch_size+1
		idxe = min(ll*batch_size,num_vectorsL)
		call GetInputVector_Inverse_tfqmr_batch_LL(level_c,rowblock,idxe-idxs+1,DimMax,RandVectInL(1:DimMax,idxs:idxe),RandVectOutL(1:DimMax,idxs:idxe),itermax)		
	end do			
	RandVectOutL = RandVectOutL - RandVectInL	
	
	mem_vec=0
	mem_vec =mem_vec + SIZEOF(RandVectInL)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectInR)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectOutL)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectOutR)/1024.0d3
	Memory_int_vec = max(Memory_int_vec,mem_vec)
	
    return                

end subroutine Get_Randomized_Vectors_inverse_Memeff






subroutine GetInputVector_Inverse_direct(level_c,rowblock,num_vectors,N,V_R,I_R,V_L,I_L)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	complex(kind=8), allocatable :: Vin(:,:),Vout(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::V_R(N,num_vectors),I_R(N,num_vectors),V_L(N,num_vectors),I_L(N,num_vectors)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)
	complex(kind=8):: al,be
	al=1d0
	be=0d0
	
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput_tmp(2))

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	allocate (RandomVectors_InOutput_tmp(1)%vector(nn,nn))
    allocate (RandomVectors_InOutput_tmp(2)%vector(nn,nn))

	do ii =1,2
		RandomVectors_InOutput_tmp(ii)%vector = 0
	end do
	
	idx_start_glo = basis_group(groupn)%head	 
	   
    do ii=1, nn
        RandomVectors_InOutput_tmp(1)%vector(ii,ii)=1
	end do
	
	! get the full forward matrix
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1

	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! compute the matrix to be inverted from the updated forward blocks 
	
	
    ctemp1=1.0d0 ; ctemp2=0.0d0
	! write(*,*)shape(RandomVectors_InOutput_tmp(2)%vector),idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	call butterfly_block_MVP_randomized_dat(block_off1,'N',idx_end_loc1-idx_start_loc1+1,idx_end_loc2-idx_start_loc2+1,nn,&
	&RandomVectors_InOutput_tmp(1)%vector(idx_start_loc2:idx_end_loc2,1:nn),RandomVectors_InOutput_tmp(2)%vector(idx_start_loc1:idx_end_loc1,1:nn),ctemp1,ctemp2)	
	RandomVectors_InOutput_tmp(2)%vector(idx_start_loc1:idx_end_loc1,1:nn) = RandomVectors_InOutput_tmp(2)%vector(idx_start_loc1:idx_end_loc1,1:nn) + &
																			&RandomVectors_InOutput_tmp(1)%vector(idx_start_loc1:idx_end_loc1,1:nn)
	call butterfly_block_MVP_randomized_dat(block_off2,'N',idx_end_loc2-idx_start_loc2+1,idx_end_loc1-idx_start_loc1+1,nn,&
	&RandomVectors_InOutput_tmp(1)%vector(idx_start_loc1:idx_end_loc1,1:nn),RandomVectors_InOutput_tmp(2)%vector(idx_start_loc2:idx_end_loc2,1:nn),ctemp1,ctemp2)	
	RandomVectors_InOutput_tmp(2)%vector(idx_start_loc2:idx_end_loc2,1:nn) = RandomVectors_InOutput_tmp(2)%vector(idx_start_loc2:idx_end_loc2,1:nn) + &
																			&RandomVectors_InOutput_tmp(1)%vector(idx_start_loc2:idx_end_loc2,1:nn)	
																																		
	allocate(matrixtemp0(nn,nn))
	matrixtemp0 = RandomVectors_InOutput_tmp(2)%vector
	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! compute the matrix to be inverted exactly 	
	
	! group_m1 = groupm_off1
	! group_n1 = groupm_off2
	! mm=basis_group(group_m1)%tail-basis_group(group_m1)%head+1
	! kk=basis_group(group_n1)%tail-basis_group(group_n1)%head+1
	! header_m=basis_group(group_m1)%head
	! header_n=basis_group(group_n1)%head
	
	! allocate(fullmat1(mm,kk))
	! allocate(fullmat2(kk,mm))
	! do m = 1,mm
		! do k = 1,kk
			! edge_m = header_m + m-1
			! edge_n = header_n + k-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! fullmat1(m,k) = ctemp
			! call element_Zmn(edge_n,edge_m,ctemp)
			! fullmat2(k,m) = ctemp
		! end do
	! end do	
	! allocate(matrixtmp1(mm,mm))
	! do m=1,mm
		! do k =1,mm
			! edge_m = header_m + m-1
			! edge_n = header_m + k-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! matrixtmp1(m,k) = ctemp	
		! end do
	! end do
	! allocate(ipiv(mm))
	! call getrff90(matrixtmp1,ipiv)
	! ! write(*,*)shape(matrixtemp1)
	! call getrif90(matrixtmp1,ipiv)	
	! deallocate(ipiv)
	! allocate(matrixtmp2(kk,kk))
	! do m=1,kk
		! do k =1,kk
			! edge_m = header_n + m-1
			! edge_n = header_n + k-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! matrixtmp2(m,k) = ctemp			
		! end do
	! end do
	! allocate(ipiv(kk))
	! call getrff90(matrixtmp2,ipiv)
	! ! write(*,*)shape(matrixtemp1)
	! call getrif90(matrixtmp2,ipiv)	
	! deallocate(ipiv)				
	! allocate(fullmat3(mm,kk))
	! call gemmf90(matrixtmp1, fullmat1, fullmat3,'N','N') 
	! allocate(fullmat4(kk,mm))
	! call gemmf90(matrixtmp2, fullmat2, fullmat4,'N','N')
	! allocate(fullmat_eye(mm+kk,mm+kk))
	! fullmat_eye = 0d0
	! do ii =1, mm+kk								
		! fullmat_eye(ii,ii) = 1d0									
	! end do
	! allocate(fullmat(mm+kk,mm+kk))
	! fullmat = fullmat_eye								
	! fullmat(1:mm,1+mm:kk+mm) = fullmat3
	! fullmat(mm+1:mm+kk,1:mm) = fullmat4	
	! allocate(matrixtemp0(nn,nn))
	! matrixtemp0 = fullmat
	! deallocate(fullmat1,fullmat2,fullmat,fullmat3,fullmat4,fullmat_eye,matrixtmp2,matrixtmp1)
	

	
	! ! if(level_butterfly==4)then
		! ! do ii=1,nn
		! ! do jj=1,nn
			! ! write(333,*) real(matrixtemp0(ii,jj)),imag(matrixtemp0(ii,jj))
		! ! end do
		! ! end do
		
		! ! do ii=1,nn
		! ! do jj=1,num_vectors
			! ! write(334,*) real(V_R(ii,jj)),imag(V_R(ii,jj))
		! ! end do
		! ! end do	
		! ! stop
	! ! end if
	
	
	allocate(ipiv(nn))
	call getrff90(matrixtemp0,ipiv)
	! write(*,*)shape(matrixtemp0)
	call getrif90(matrixtemp0,ipiv)																				

	! write(*,*)shape(matrixtemp0),shape(V_R),shape(I_R)

	call gemmf90(matrixtemp0, V_R, I_R,'N','N',al,be)      

	call gemmf90(matrixtemp0, V_L, I_L,'T','N',al,be)    
	
	
	
	deallocate(matrixtemp0)																		
	deallocate(ipiv)
	deallocate(RandomVectors_InOutput_tmp(1)%vector)	
	deallocate(RandomVectors_InOutput_tmp(2)%vector)	
	deallocate(RandomVectors_InOutput_tmp)	

	
    return                

end subroutine GetInputVector_Inverse_direct




subroutine GetInputVector_Inverse_tfqmr(level_c,rowblock,num_vectors,N,V_R,I_R,V_L,I_L,itermax)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::V_R(N,num_vectors),I_R(N,num_vectors),V_L(N,num_vectors),I_L(N,num_vectors)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

	integer N_iter_max,iter,itermax
	real*8:: rel_error
	
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	idx_start_glo = basis_group(groupn)%head	 
	
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1


	! initialize the iterative solver 
	call initialize_r0(nn)
	I_L = 0
	I_R = 0
	itermax = 0
	
	! write(*,*)'ddd'
	! get the right multiplied input vectors
	do ii=1,num_vectors
	! write(*,*)'dddddd',ii
		N_iter_max = 100
		iter = 0
		rel_error = tfqmr_tolerance
		call ztfqmr_parallel('N',level_c,rowblock,N_iter_max,nn,V_R(:,ii),I_R(:,ii),rel_error,iter)
		itermax = max(itermax,iter)
	end do
	
	! get the left multiplied input vectors
	do ii=1,num_vectors
		N_iter_max = 100
		iter = 0
		rel_error = tfqmr_tolerance
		call ztfqmr_parallel('T',level_c,rowblock,N_iter_max,nn,V_L(:,ii),I_L(:,ii),rel_error,iter)
		itermax = max(itermax,iter)
	end do	
	
	deallocate(r0_initial)
		
    return                

end subroutine GetInputVector_Inverse_tfqmr



subroutine GetInputVector_Inverse_tfqmr_batch(level_c,rowblock,num_vectors,N,V_R,I_R,V_L,I_L,itermax)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::V_R(N,num_vectors),I_R(N,num_vectors),V_L(N,num_vectors),I_L(N,num_vectors)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

	integer N_iter_max,iter,itermax
	real*8:: rel_error
	
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	idx_start_glo = basis_group(groupn)%head	 
	
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1


	! initialize the iterative solver 
	call initialize_r0(nn)
	I_L = 0
	I_R = 0
	itermax = 0
	
	! write(*,*)'ddd'
	! get the right multiplied input vectors

	N_iter_max = 100
	iter = 0
	rel_error = tfqmr_tolerance
	! write(*,*)'in1'
	call ztfqmr_batch('N',level_c,rowblock,N_iter_max,nn,num_vectors,V_R,I_R,rel_error,iter)
	itermax = max(itermax,iter)
	! write(*,*)'out1'
	
	! get the left multiplied input vectors
	N_iter_max = 100
	iter = 0
	rel_error = tfqmr_tolerance
		! write(*,*)'in2'
	call ztfqmr_batch('T',level_c,rowblock,N_iter_max,nn,num_vectors,V_L,I_L,rel_error,iter)
	itermax = max(itermax,iter)
		! write(*,*)'out2'
	
	deallocate(r0_initial)
		
    return                

end subroutine GetInputVector_Inverse_tfqmr_batch




subroutine GetInputVector_Inverse_tfqmr_batch_LL(level_c,rowblock,num_vectors,N,V_L,I_L,itermax)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::V_L(N,num_vectors),I_L(N,num_vectors)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

	integer N_iter_max,iter,itermax
	real*8:: rel_error
	
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	idx_start_glo = basis_group(groupn)%head	 
	
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1


	! initialize the iterative solver 
	call initialize_r0(nn)
	I_L = 0
	itermax = 0
	
	! get the left multiplied input vectors
	N_iter_max = 100
	iter = 0
	rel_error = tfqmr_tolerance*1d3
		! write(*,*)'in2'
	call ztfqmr_batch('T',level_c,rowblock,N_iter_max,nn,num_vectors,V_L,I_L,rel_error,iter)
	itermax = max(itermax,iter)
		! write(*,*)'out2'
	
	deallocate(r0_initial)
		
    return                

end subroutine GetInputVector_Inverse_tfqmr_batch_LL


subroutine GetInputVector_Inverse_tfqmr_batch_RR(level_c,rowblock,num_vectors,N,V_R,I_R,itermax)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::V_R(N,num_vectors),I_R(N,num_vectors)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

	integer N_iter_max,iter,itermax
	real*8:: rel_error
	
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	idx_start_glo = basis_group(groupn)%head	 
	
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1


	! initialize the iterative solver 
	call initialize_r0(nn)
	I_R = 0
	itermax = 0
	
	! write(*,*)'ddd'
	! get the right multiplied input vectors

	N_iter_max = 100
	iter = 0
	rel_error = tfqmr_tolerance*1d3
	! write(*,*)'in1'
	call ztfqmr_batch('N',level_c,rowblock,N_iter_max,nn,num_vectors,V_R,I_R,rel_error,iter)
	itermax = max(itermax,iter)
	! write(*,*)'out1'
	
	deallocate(r0_initial)
		
    return                

end subroutine GetInputVector_Inverse_tfqmr_batch_RR

! subroutine SmartMultifly(trans,N,level_c,rowblock,num_vectors,Id,Vd)

    ! use MODULE_FILE
    ! ! use lapack95
	! ! use blas95	
    ! implicit none
    
	! character trans
	! integer level_c,rowblock
    ! integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N
    ! integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    ! character chara
    ! real*8 a,b,c,d
    ! complex(kind=8) ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o,block_off1,block_off2
	! type(matrixblock):: block_tmp	
    ! type(vectorsblock), pointer :: random1, random2
    
    ! real*8,allocatable :: Singular(:)
	! integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	! complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	! complex(kind=8), allocatable :: matrixtemp0(:,:)
	! integer, allocatable :: ipiv(:)
	! real*8::rtemp
	! complex(kind=8)::Vd(N,num_vectors),Id(N,num_vectors)
	
	! integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	! complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

    
	! block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	! block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	! block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    ! level_butterfly=maxlevel_for_blocks-block_o%level
    ! num_blocks=2**level_butterfly

    ! groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	! idx_start_glo = basis_group(groupn)%head	 
	
	! groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	! idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	! idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	! groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	! idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	! idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1

	! Vd = 0
	
		! !  write(*,*)'a',trans
	! if(trans=='N')then
		! ctemp1=1.0d0
		! ctemp2=0.0d0
		! call butterfly_block_MVP_randomized_dat(block_off1,'N',idx_end_loc1-idx_start_loc1+1,idx_end_loc2-idx_start_loc2+1,num_vectors,&
		! &Id(idx_start_loc2:idx_end_loc2,1:num_vectors),Vd(idx_start_loc1:idx_end_loc1,1:num_vectors),ctemp1,ctemp2)	
		! !  write(*,*)'aha'
		! Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) = Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) + Id(idx_start_loc1:idx_end_loc1,1:num_vectors)
		! call butterfly_block_MVP_randomized_dat(block_off2,'N',idx_end_loc2-idx_start_loc2+1,idx_end_loc1-idx_start_loc1+1,num_vectors,&
		! &Id(idx_start_loc1:idx_end_loc1,1:num_vectors),Vd(idx_start_loc2:idx_end_loc2,1:num_vectors),ctemp1,ctemp2)	
		! Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) = Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) + Id(idx_start_loc2:idx_end_loc2,1:num_vectors)							
	! else 
		! ctemp1=1.0d0
		! ctemp2=0.0d0
		! call butterfly_block_MVP_randomized_dat(block_off2,'T',idx_end_loc1-idx_start_loc1+1,idx_end_loc2-idx_start_loc2+1,num_vectors,&
		! &Id(idx_start_loc2:idx_end_loc2,1:num_vectors),Vd(idx_start_loc1:idx_end_loc1,1:num_vectors),ctemp1,ctemp2)	
		! Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) = Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) + Id(idx_start_loc1:idx_end_loc1,1:num_vectors)	
		! call butterfly_block_MVP_randomized_dat(block_off1,'T',idx_end_loc2-idx_start_loc2+1,idx_end_loc1-idx_start_loc1+1,num_vectors,&
		! &Id(idx_start_loc1:idx_end_loc1,1:num_vectors),Vd(idx_start_loc2:idx_end_loc2,1:num_vectors),ctemp1,ctemp2)	
		! Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) = Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) + Id(idx_start_loc2:idx_end_loc2,1:num_vectors)		
	! end if

	
		! !  write(*,*)'b'
! end subroutine SmartMultifly



subroutine SmartMultifly(trans,N,level_c,rowblock,num_vectors,Id,Vd)

    use MODULE_FILE
    ! use lapack95
	! use blas95	
    implicit none
    
	character trans
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,N,Nmax,Nmax1,Nmax2
    integer mm,nn,mn,m,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2
	type(matrixblock):: block_tmp	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8), allocatable :: matrixtemp0(:,:)
	integer, allocatable :: ipiv(:)
	real*8::rtemp
	complex(kind=8)::Vd(N,num_vectors),Id(N,num_vectors)
	complex(kind=8),allocatable::Id_tmp(:,:),Vd_tmp(:,:)
	
	integer header_m,header_n,group_m1,group_n1,group_m2,group_n2,edge_m,edge_n
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:)

    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	call CountMaxIntermidiateVector(block_off1,Nmax1)
	call CountMaxIntermidiateVector(block_off2,Nmax2)
	Nmax = max(Nmax1,Nmax2)
	allocate(Id_tmp(Nmax,num_vectors))
	allocate(Vd_tmp(Nmax,num_vectors))
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly

    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	idx_start_glo = basis_group(groupn)%head	 
	
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1

	Vd = 0
	
		!  write(*,*)'a',trans
	if(trans=='N')then
		Id_tmp(1:idx_end_loc2-idx_start_loc2+1,1:num_vectors) = Id(idx_start_loc2:idx_end_loc2,1:num_vectors)
		call Butterfly_Partial_MVP_dat_memeff(block_off1,'N',0,block_off1%level_butterfly+1,N,num_vectors,Id_tmp,Vd_tmp)
		Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) = Vd_tmp(1:idx_end_loc1-idx_start_loc1+1,1:num_vectors) + Id(idx_start_loc1:idx_end_loc1,1:num_vectors)
		
		Id_tmp(1:idx_end_loc1-idx_start_loc1+1,1:num_vectors)=Id(idx_start_loc1:idx_end_loc1,1:num_vectors)
		call Butterfly_Partial_MVP_dat_memeff(block_off2,'N',0,block_off2%level_butterfly+1,N,num_vectors,Id_tmp,Vd_tmp)
		Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) = Vd_tmp(1:idx_end_loc2-idx_start_loc2+1,1:num_vectors) + Id(idx_start_loc2:idx_end_loc2,1:num_vectors)
	else 
		Id_tmp(1:idx_end_loc2-idx_start_loc2+1,1:num_vectors) = Id(idx_start_loc2:idx_end_loc2,1:num_vectors)
		call Butterfly_Partial_MVP_dat_memeff(block_off2,'T',0,block_off2%level_butterfly+1,N,num_vectors,Id_tmp,Vd_tmp)
		Vd(idx_start_loc1:idx_end_loc1,1:num_vectors) = Vd_tmp(1:idx_end_loc1-idx_start_loc1+1,1:num_vectors) + Id(idx_start_loc1:idx_end_loc1,1:num_vectors)
		
		Id_tmp(1:idx_end_loc1-idx_start_loc1+1,1:num_vectors)=Id(idx_start_loc1:idx_end_loc1,1:num_vectors)
		call Butterfly_Partial_MVP_dat_memeff(block_off1,'T',0,block_off1%level_butterfly+1,N,num_vectors,Id_tmp,Vd_tmp)
		Vd(idx_start_loc2:idx_end_loc2,1:num_vectors) = Vd_tmp(1:idx_end_loc2-idx_start_loc2+1,1:num_vectors) + Id(idx_start_loc2:idx_end_loc2,1:num_vectors)
	end if

	
	deallocate(Id_tmp)
	deallocate(Vd_tmp)
	
		!  write(*,*)'b'
end subroutine SmartMultifly













subroutine Test_Randomized_Inversion(level_c,rowblock,error)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,groupm_off1,groupm_off2
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_off1,block_off2	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc1,idx_end_loc1,idx_start_loc2,idx_end_loc2
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	real*8::error
   
    
	block_o => cascading_factors(level_c)%matrices_block_inverse(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(6))

    num_vectors=1  
	
    groupn=block_o%col_group    ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(3)%vector(nn,num_vectors))    
	allocate (RandomVectors_InOutput(4)%vector(nn,num_vectors))
    allocate (RandomVectors_InOutput(5)%vector(nn,num_vectors))
	allocate (RandomVectors_InOutput(6)%vector(nn,num_vectors))	 
	do ii =1,6
		RandomVectors_InOutput(ii)%vector = 0
	end do
	idx_start_glo = basis_group(groupn)%head	 
	   
    ! !$omp parallel do default(shared) private(ii,jj,a,b,c,d)
    do jj=1, num_vectors
        do ii=1, nn
            call random_number(a)
            call random_number(b)
            call random_number(c)
            call random_number(d)
            if (c<0.5d0) then
                a=-a
            endif
            if (d<0.5d0) then
                b=-b
            endif
            RandomVectors_InOutput(2)%vector(ii,jj)=1
        enddo
        do ii=1, nn
            call random_number(a)
            call random_number(b)
            call random_number(c)
            call random_number(d)
            if (c<0.5) then
                a=-a
            endif
            if (d<0.5) then
                b=-b
            endif
            RandomVectors_InOutput(5)%vector(ii,jj)=a+junit*b
        enddo
    enddo
    ! !$omp end parallel do
    
	
	
	! get the right multiplied vectors
	groupm_off1 = block_off1%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc1 = basis_group(groupm_off1)%head-idx_start_glo+1
	idx_end_loc1 = basis_group(groupm_off1)%tail-idx_start_glo+1
	groupm_off2 = block_off2%row_group  ! Note: row_group and col_group interchanged here   
	idx_start_loc2 = basis_group(groupm_off2)%head-idx_start_glo+1
	idx_end_loc2 = basis_group(groupm_off2)%tail-idx_start_glo+1
	
    ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(block_off1,'N',idx_end_loc1-idx_start_loc1+1,idx_end_loc2-idx_start_loc2+1,num_vectors,&
	&RandomVectors_InOutput(2)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors),RandomVectors_InOutput(1)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors),ctemp1,ctemp2)	
	RandomVectors_InOutput(1)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors) = RandomVectors_InOutput(1)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors) + &
																			&RandomVectors_InOutput(2)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors)
	call butterfly_block_MVP_randomized_dat(block_off2,'N',idx_end_loc2-idx_start_loc2+1,idx_end_loc1-idx_start_loc1+1,num_vectors,&
	&RandomVectors_InOutput(2)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors),RandomVectors_InOutput(1)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors),ctemp1,ctemp2)	
	RandomVectors_InOutput(1)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors) = RandomVectors_InOutput(1)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors) + &
																			&RandomVectors_InOutput(2)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors)	

	! write(*,*)level_c,rowblock
	! write(*,*)block_o%ButterflyV(1)%matrix
	! write(*,*)'haha'
	call butterfly_block_MVP_randomized_dat(block_o,'N',nn,nn,num_vectors,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	
	
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(1)%vector + RandomVectors_InOutput(3)%vector	
	
	error = sqrt(sum(abs(RandomVectors_InOutput(2)%vector(:,1)-RandomVectors_InOutput(3)%vector(:,1))**2)/sum(abs(RandomVectors_InOutput(2)%vector(:,1))**2))
	! write(*,*)error,sum(abs(RandomVectors_InOutput(2)%vector(:,1))**2),sum(abs(RandomVectors_InOutput(3)%vector(:,1))**2)
	
	
	
	! get the left multiplied vectors	
	call butterfly_block_MVP_randomized_dat(block_off2,'T',idx_end_loc2-idx_start_loc2+1,idx_end_loc1-idx_start_loc1+1,num_vectors,&
	&RandomVectors_InOutput(5)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors),RandomVectors_InOutput(4)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors),ctemp1,ctemp2)	
	RandomVectors_InOutput(4)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors) = RandomVectors_InOutput(4)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors) + &
																			&RandomVectors_InOutput(5)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors)	
	call butterfly_block_MVP_randomized_dat(block_off1,'T',idx_end_loc1-idx_start_loc1+1,idx_end_loc2-idx_start_loc2+1,num_vectors,&
	&RandomVectors_InOutput(5)%vector(idx_start_loc1:idx_end_loc1,1:num_vectors),RandomVectors_InOutput(4)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors),ctemp1,ctemp2)	
	RandomVectors_InOutput(4)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors) = RandomVectors_InOutput(4)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors) + &
																			&RandomVectors_InOutput(5)%vector(idx_start_loc2:idx_end_loc2,1:num_vectors)		

	call butterfly_block_MVP_randomized_dat(block_o,'T',nn,nn,num_vectors,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(6)%vector,ctemp1,ctemp2)	
	
	
	RandomVectors_InOutput(6)%vector = RandomVectors_InOutput(4)%vector + RandomVectors_InOutput(6)%vector	
	
	error = max(error,sqrt(sum(abs(RandomVectors_InOutput(5)%vector(:,1)-RandomVectors_InOutput(6)%vector(:,1))**2)/sum(abs(RandomVectors_InOutput(5)%vector(:,1))**2)))																			
	
	! write(*,*)error
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 6
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)	
	
	
    return                

end subroutine Test_Randomized_Inversion



  !*********************************************************************
  !  this subroutine solves the system ax=b using transpose free quasi-
  !  minimal residual (tfqmr) algorithm. reference:
  !    siam j. sci. compt. vol.14, no.2, pp. 470-482, march 93
  !       by roland. w. freund
  !
  !  the program terminates when the required precision is reached.
  !  if the required precision is not established only n iterations
  !  are carried out.     a.a.ergin may 1995
  !
  ! Modified by a.e. yilmaz 2003
  ! needs mat_vec_mult and initialize_r0 routine
  ! that computes an initial random vector r0_initial
  subroutine ztfqmr_parallel(trans,level_c,rowblock,ntotal,nn,b,x,err,iter)
    implicit none
	integer level_c,rowblock
    integer,intent(in)::ntotal
    integer::iter,itmax,it,nn	
    complex(kind=dp),dimension(1:nn)::x,bb,b
    real(kind=dp)::err,rerr
    complex(kind=dp),dimension(1:nn)::w,yo,ayo,ye,aye,r,d,v
    real(kind=dp)::ta,we,cm
    complex(kind=dp)::we_local,we_sum,rerr_local,rerr_sum,err_local,err_sum
    complex(kind=dp)::ta_local,ta_sum,bmag_local,bmag_sum1,dumb_ali(6)
    complex(kind=dp)::etha,rho,rho_local,amgis,amgis_local,ahpla,dum,dum_local,beta
    real(kind=dp)::bmag
    real(kind=dp)::tim1,tim2
    integer::kk,ll
    ! Variables for storing current
    integer::count_unk,srcbox_id
    complex(kind=dp),dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    real(kind=dp)::mem_est
	character:: trans

    itmax=iter

    ! ! ! if (myid == main_id) then
       ! ! ! call cpu_time(tim1)
       ! ! ! open(unit=32,file='iterations.out',status='unknown')
    ! ! ! end if
    
    if (iter.eq.0) itmax=ntotal
    bb=b 
    
    !  set initial values
    !
    d=cmplx(0.0_dp,0.0_dp,dp)
    ! write(*,*)'1'
	call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)    
    
    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if		
	call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)    

    v=ayo
    we=0.0_dp
    etha=cmplx(0.0_dp,0.0_dp,dp)
    
    ta_local=dot_product(r,r)
    rho_local=dot_product(r0_initial,r)
    bmag_local=dot_product(bb,bb)
    
    dumb_ali(1:3)=(/ta_local,rho_local,bmag_local/)
    ! call MPI_ALLREDUCE(dumb_ali(1:3),dumb_ali(4:6),3,MPI_DOUBLE_COMPLEX,&
         ! MPI_SUM,MPI_Comm_world,ierr)
	dumb_ali(4:6) = dumb_ali(1:3)	 
    ta_sum=dumb_ali(4);rho=dumb_ali(5);bmag_sum1=dumb_ali(6)
    ta=sqrt(abs(ta_sum))
    bmag=sqrt(abs(bmag_sum1))
    rerr=ta/bmag
    
    iters: do it=1,itmax
       amgis_local=dot_product(r0_initial,v)
       ! ! call MPI_ALLREDUCE(amgis_local,amgis,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! ! MPI_Comm_world,ierr)
	   amgis = 	amgis_local
       ahpla=rho/amgis
       ye=yo-ahpla*v
           ! write(*,*)'3'
       call SmartMultifly(trans,nn,level_c,rowblock,1,ye,aye)
	
       !  start odd (2n-1) m loop
       d=yo+(we*we*etha/ahpla)*d
       w=w-ahpla*ayo
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
       !  check if the result has converged.
       !a        if (err*bmag .gt. ta*sqrt(2.*it)) then
       !
       !  start even (2n)  m loop
       d=ye+(we*we*etha/ahpla)*d
       w=w-ahpla*aye
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
      
       
       !  check if the result has converged.
       if (mod(it,1)==0 .or. rerr<1.0_dp*err) then
    ! write(*,*)'4'
		  call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
	   
          r=bb-r
          rerr_local=dot_product(r,r)
          ! call MPI_ALLREDUCE(rerr_local,rerr_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
               ! MPI_Comm_world,ierr)
			rerr_sum = rerr_local   
          rerr=sqrt(abs(rerr_sum))/bmag
          
          ! ! if (myid==main_id) then
             ! print*,'# ofiter,error:',it,rerr
             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          ! ! end if
          
          if (err > rerr) then
             err=rerr
             iter=it

             ! ! ! if (myid == main_id) then
                ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
             ! ! ! end if
             
             return
          endif
       end if
       !  make preparations for next iteration
       dum_local=dot_product( r0_initial,w)
       ! call MPI_ALLREDUCE(dum_local,dum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		dum = dum_local	
       beta=dum/rho
       rho=dum
       yo=w+beta*ye
           ! write(*,*)'5'
       call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)
       
       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
    call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
	   
    !MAGIC
    r=bb-r
    err_local=dot_product(r,r)
    ! call MPI_ALLREDUCE(err_local,err_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_Comm_world,&
         ! ierr)
	err_sum = err_local	 
    err=sqrt(abs(err_sum))/bmag
    iter=itmax



   print*,'Iterative solver is terminated without convergence!!!',it,err
   stop
	
    return
  end subroutine ztfqmr_parallel


  
  
    subroutine ztfqmr_fullmat(ntotal,nn,A,b,x,err,iter)
    implicit none
	integer level_c,rowblock
    integer,intent(in)::ntotal
    integer::iter,itmax,it,nn,mm,ii,jj	
    complex(kind=dp),dimension(1:nn)::x,bb,b
    complex(kind=dp),dimension(1:nn,1:nn)::A
    real(kind=dp)::err,rerr
    complex(kind=dp),dimension(1:nn)::w,yo,ayo,ye,aye,r,d,v
    real(kind=dp)::ta,we,cm
    complex(kind=dp)::we_local,we_sum,rerr_local,rerr_sum,err_local,err_sum,ctemp
    complex(kind=dp)::ta_local,ta_sum,bmag_local,bmag_sum1,dumb_ali(6)
    complex(kind=dp)::etha,rho,rho_local,amgis,amgis_local,ahpla,dum,dum_local,beta
    real(kind=dp)::bmag
    real(kind=dp)::tim1,tim2
    integer::kk,ll
    ! Variables for storing current
    integer::count_unk,srcbox_id
    complex(kind=dp),dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    real(kind=dp)::mem_est
	character:: trans

    itmax=iter

	
       allocate(r0_initial(1:nn))

		do ii=1,nn
		   r0_initial(ii)= random_complex_number()
		end do		
	
	
    ! ! ! if (myid == main_id) then
       ! ! ! call cpu_time(tim1)
       ! ! ! open(unit=32,file='iterations.out',status='unknown')
    ! ! ! end if
    
    if (iter.eq.0) itmax=ntotal
    bb=b 
    
    !  set initial values
    !
    d=cmplx(0.0_dp,0.0_dp,dp)
    ! write(*,*)'1'
	
	
	do ii =1,nn
	ctemp = 0
	do kk=1,nn
		ctemp = ctemp + A(ii,kk)*x(kk)	
	end do
	r(ii) = ctemp
	end do

	! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)    
    
    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if	


	do ii =1,nn
	ctemp = 0
	do kk=1,nn
		ctemp = ctemp + A(ii,kk)*yo(kk)	
	end do
	ayo(ii) = ctemp
	end do
	
	! call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)    

    v=ayo
    we=0.0_dp
    etha=cmplx(0.0_dp,0.0_dp,dp)
    
    ta_local=dot_product(r,r)
    rho_local=dot_product(r0_initial,r)
    bmag_local=dot_product(bb,bb)
    
    dumb_ali(1:3)=(/ta_local,rho_local,bmag_local/)
    ! call MPI_ALLREDUCE(dumb_ali(1:3),dumb_ali(4:6),3,MPI_DOUBLE_COMPLEX,&
         ! MPI_SUM,MPI_Comm_world,ierr)
	dumb_ali(4:6) = dumb_ali(1:3)	 
    ta_sum=dumb_ali(4);rho=dumb_ali(5);bmag_sum1=dumb_ali(6)
    ta=sqrt(abs(ta_sum))
    bmag=sqrt(abs(bmag_sum1))
    rerr=ta/bmag
    
    iters: do it=1,itmax
       amgis_local=dot_product(r0_initial,v)
       ! ! call MPI_ALLREDUCE(amgis_local,amgis,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! ! MPI_Comm_world,ierr)
	   amgis = 	amgis_local
       ahpla=rho/amgis
       ye=yo-ahpla*v
           ! write(*,*)'3'
		   
		   
		do ii =1,nn
		ctemp = 0
		do kk=1,nn
			ctemp = ctemp + A(ii,kk)*ye(kk)	
		end do
		aye(ii) = ctemp
		end do		   
       ! call SmartMultifly(trans,nn,level_c,rowblock,1,ye,aye)
	
       !  start odd (2n-1) m loop
       d=yo+(we*we*etha/ahpla)*d
       w=w-ahpla*ayo
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
       !  check if the result has converged.
       !a        if (err*bmag .gt. ta*sqrt(2.*it)) then
       !
       !  start even (2n)  m loop
       d=ye+(we*we*etha/ahpla)*d
       w=w-ahpla*aye
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
      
       
       !  check if the result has converged.
       if (mod(it,1)==0 .or. rerr<1.0_dp*err) then
    ! write(*,*)'4'
	
		do ii =1,nn
		ctemp = 0
		do kk=1,nn
			ctemp = ctemp + A(ii,kk)*x(kk)	
		end do
		r(ii) = ctemp
		end do		
		  ! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
	   
          r=bb-r
          rerr_local=dot_product(r,r)
          ! call MPI_ALLREDUCE(rerr_local,rerr_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
               ! MPI_Comm_world,ierr)
			rerr_sum = rerr_local   
          rerr=sqrt(abs(rerr_sum))/bmag
          
          ! ! if (myid==main_id) then
             print*,'# ofiter,error:',it,rerr
             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          ! ! end if
          
          if (err > rerr) then
             err=rerr
             iter=it

             ! ! ! if (myid == main_id) then
                ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
             ! ! ! end if
             deallocate(r0_initial)
             return
          endif
       end if
       !  make preparations for next iteration
       dum_local=dot_product( r0_initial,w)
       ! call MPI_ALLREDUCE(dum_local,dum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		dum = dum_local	
       beta=dum/rho
       rho=dum
       yo=w+beta*ye
           ! write(*,*)'5'
		   
		do ii =1,nn
		ctemp = 0
		do kk=1,nn
			ctemp = ctemp + A(ii,kk)*yo(kk)	
		end do
		ayo(ii) = ctemp
		end do			   
       ! call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)
       
       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
		do ii =1,nn
		ctemp = 0
		do kk=1,nn
			ctemp = ctemp + A(ii,kk)*x(kk)	
		end do
		r(ii) = ctemp
		end do		
    ! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
	   
    !MAGIC
    r=bb-r
    err_local=dot_product(r,r)
    ! call MPI_ALLREDUCE(err_local,err_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_Comm_world,&
         ! ierr)
	err_sum = err_local	 
    err=sqrt(abs(err_sum))/bmag
    iter=itmax

	deallocate(r0_initial)

   print*,'Iterative solver is terminated without convergence!!!',it,err
   stop
	
    return
  end subroutine ztfqmr_fullmat

  subroutine ztfqmr_batch(trans,level_c,rowblock,ntotal,nn,nvec,b,x,err,iter)
    implicit none
	integer level_c,rowblock
    integer,intent(in)::ntotal
    integer::iter,itmax,it,nn,nvec,ii	
    complex(kind=dp),dimension(1:nn,nvec)::x,b
    real(kind=dp)::err,rerr(nvec),rerr_old(nvec)
    ! complex(kind=dp),dimension(1:nn,nvec)::w,yo,ayo,ye,aye,r,d,v,bb
    complex(kind=dp),dimension(:,:),allocatable::w,yo,ayo,ye,aye,r,d,v
	real(kind=dp)::ta(nvec),we(nvec),cm(nvec)
    complex(kind=dp)::we_local,we_sum(nvec),rerr_local,rerr_sum(nvec),err_local,err_sum(nvec)
    complex(kind=dp)::ta_local,ta_sum(nvec),bmag_local(nvec),bmag_sum1(nvec),dumb_ali(6)
    complex(kind=dp)::etha(nvec),rho(nvec),rho_local(nvec),amgis(nvec),amgis_local(nvec),ahpla(nvec),dum(nvec),dum_local(nvec),beta(nvec)
    real(kind=dp)::bmag(nvec)
    real(kind=dp)::tim1,tim2
    integer::kk,ll
    ! Variables for storing current
    integer::count_unk,srcbox_id
    complex(kind=dp),dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    real(kind=dp)::mem_est
	character:: trans
	integer::Ngood,Nbad
	real(kind=dp)::mem_tmp
	mem_tmp = SIZEOF(x)/1024.0d3
	mem_tmp = mem_tmp*11d0
	Memory_tfqmr_vec = max(Memory_tfqmr_vec,mem_tmp)
	
	allocate(w(nn,Nvec))
	allocate(yo(nn,Nvec))
	allocate(ayo(nn,Nvec))
	allocate(ye(nn,Nvec))
	allocate(aye(nn,Nvec))
	allocate(r(nn,Nvec))
	allocate(d(nn,Nvec))
	allocate(v(nn,Nvec))
	
	
    itmax=iter

    if (iter.eq.0) itmax=ntotal
    ! bb=b 
    
    !  set initial values
    !
    d=cmplx(0.0_dp,0.0_dp,dp)
	! ! write(*,*)'1'
	call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)    
    
    r=b-r !residual from the initial guess
    
	! ! do ii =1,nvec
		! ! write(*,*)sum(abs(r(:,ii))**2)
	! ! end do
	
	
	w=r
    yo=r	
	! write(*,*)'2'
	call SmartMultifly(trans,nn,level_c,rowblock,nvec,yo,ayo)    
	       		  ! ! write(*,*)abs(ayo)
					! ! pause
    v=ayo
    we=0.0_dp
    etha=cmplx(0.0_dp,0.0_dp,dp)
    
	call dot_product_batch(nn,nvec,r,r,ta_sum)
	call dot_product_one_time_batch(nn,nvec,r0_initial,r,rho)	
	call dot_product_batch(nn,nvec,b,b,bmag_sum1)	
	! ! write(*,*)'d',abs(rho),'dd'
	! ! pause

    ta=sqrt(abs(ta_sum))
    bmag=sqrt(abs(bmag_sum1))
	do ii =1,nvec
		rerr(ii)=ta(ii)/bmag(ii)
	end do
	
    iters: do it=1,itmax
	   call dot_product_one_time_batch(nn,nvec,r0_initial,v,amgis)
	
	do ii =1,nvec
		ahpla(ii)=rho(ii)/amgis(ii)
		ye(:,ii)=yo(:,ii)-ahpla(ii)*v(:,ii)
	end do	   
! write(*,*)'3'
       call SmartMultifly(trans,nn,level_c,rowblock,nvec,ye,aye)
				  
	do ii =1,nvec
		d(:,ii)=yo(:,ii)+(we(ii)*we(ii)*etha(ii)/ahpla(ii))*d(:,ii)
		w(:,ii)=w(:,ii)-ahpla(ii)*ayo(:,ii)
	end do	
	call dot_product_batch(nn,nvec,w,w,we_sum)
	do ii =1,nvec
		we(ii)=sqrt(abs(we_sum(ii)))/ta(ii)
		cm(ii)=1.0d0/sqrt(1.0d0+we(ii)*we(ii))
		ta(ii)=ta(ii)*we(ii)*cm(ii)
		etha(ii)=ahpla(ii)*cm(ii)*cm(ii)
		x(:,ii)=x(:,ii)+etha(ii)*d(:,ii)
		d(:,ii)=ye(:,ii)+(we(ii)*we(ii)*etha(ii)/ahpla(ii))*d(:,ii)
		w(:,ii)=w(:,ii)-ahpla(ii)*aye(:,ii)
	end do
	

	
	call dot_product_batch(nn,nvec,w,w,we_sum)
	do ii =1,nvec
		we(ii)=sqrt(abs(we_sum(ii)))/ta(ii)
		cm(ii)=1.0d0/sqrt(1.0d0+we(ii)*we(ii))
		ta(ii)=ta(ii)*we(ii)*cm(ii)
		etha(ii)=ahpla(ii)*cm(ii)*cm(ii)
		x(:,ii)=x(:,ii)+etha(ii)*d(:,ii)		
	end do   
	   


       !  check if the result has converged.
       if (it>=6 .or. maxval(rerr)<1.0_dp*err) then
		  ! ! write(*,*)'4'
		  call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)		
          r=b-r

		  call dot_product_batch(nn,nvec,r,r,rerr_sum)
		  do ii =1,nvec
			  rerr(ii)=sqrt(abs(rerr_sum(ii)))/bmag(ii)
		  end do
		  
		  if(isnan(sum(rerr)))then
			write(*,*)it,rerr_old
			stop
		  end if
		  
          ! ! if (myid==main_id) then
		  ! if(level_c==2 .or. level_c==1)then
			 Ngood = count(rerr<err)
			 Nbad = count(rerr>err)
             ! ! print*,'# ofiter,error:',it,maxval(rerr),Ngood,Nbad
		  ! end if
			 ! ! pause
             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          ! ! end if
          
          if (err > maxval(rerr) .or. Ngood>Nbad*10) then
             err=maxval(rerr)
             iter=it

             ! ! ! if (myid == main_id) then
                ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
             ! ! ! end if
			deallocate(w,yo,ayo,ye,aye,r,d,v)
             
             return
          endif
       end if
       !  make prepdeallocate(w,yo,ayo,ye,aye,r,d,v)
		call dot_product_one_time_batch(nn,nvec,r0_initial,w,dum)	   
		do ii =1,nvec
			beta(ii)=dum(ii)/rho(ii)
			rho(ii)=dum(ii)
			yo(:,ii)=w(:,ii)+beta(ii)*ye(:,ii)
		end do		
		
       call SmartMultifly(trans,nn,level_c,rowblock,nvec,yo,ayo)
       
       !MAGIC
	   do ii =1,nvec
		  v(:,ii)=ayo(:,ii)+beta(ii)*( aye(:,ii)+beta(ii)*v(:,ii) )
	   end do
	   rerr_old = rerr
	   
    enddo iters
    ! write(*,*)'6'
    call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)
	   
    !MAGIC
    r=b-r
	call dot_product_batch(nn,nvec,r,r,err_sum)	
	err = 0
	do ii =1,nvec	 
		err=max(err,sqrt(abs(err_sum(ii)))/bmag(ii))
    end do
	
	iter=itmax
	deallocate(w,yo,ayo,ye,aye,r,d,v)


   print*,'Iterative solver is terminated without convergence!!!',it,err,rerr
   stop
	
    return
  end subroutine ztfqmr_batch  

  
  ! subroutine ztfqmr_batch(trans,level_c,rowblock,ntotal,nn,nvec,b,x,err,iter)
    ! implicit none
	! integer level_c,rowblock
    ! integer,intent(in)::ntotal
    ! integer::iter,itmax,it,nn,nvec,ii	
    ! complex(kind=dp),dimension(1:nn,nvec)::x,bb,b
    ! real(kind=dp)::err,rerr(nvec),rerr_old(nvec)
    ! complex(kind=dp),dimension(1:nn,nvec)::w,yo,ayo,ye,aye,r,d,v
    ! real(kind=dp)::ta(nvec),we(nvec),cm(nvec)
    ! complex(kind=dp)::we_local,we_sum(nvec),rerr_local,rerr_sum(nvec),err_local,err_sum(nvec)
    ! complex(kind=dp)::ta_local,ta_sum(nvec),bmag_local(nvec),bmag_sum1(nvec),dumb_ali(6)
    ! complex(kind=dp)::etha(nvec),rho(nvec),rho_local(nvec),amgis(nvec),amgis_local(nvec),ahpla(nvec),dum(nvec),dum_local(nvec),beta(nvec)
    ! real(kind=dp)::bmag(nvec)
    ! real(kind=dp)::tim1,tim2
    ! integer::kk,ll
    ! ! Variables for storing current
    ! integer::count_unk,srcbox_id
    ! complex(kind=dp),dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    ! real(kind=dp)::mem_est
	! character:: trans
	! integer::Ngood,Nbad
	! real(kind=dp)::mem_tmp
	! mem_tmp = SIZEOF(x)/1024.0d3
	! mem_tmp = mem_tmp*11d0
	! Memory_int_vec = max(Memory_int_vec,mem_tmp)
	
    ! itmax=iter

    ! if (iter.eq.0) itmax=ntotal
    ! bb=b 
    
    ! !  set initial values
    ! !
    ! d=cmplx(0.0_dp,0.0_dp,dp)
	! ! ! write(*,*)'1'
	! call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)    
    
    ! r=bb-r !residual from the initial guess
    
	! ! ! do ii =1,nvec
		! ! ! write(*,*)sum(abs(r(:,ii))**2)
	! ! ! end do
	
	
	! w=r
    ! yo=r	
	! ! write(*,*)'2'
	! call SmartMultifly(trans,nn,level_c,rowblock,nvec,yo,ayo)    
	       		  ! ! ! write(*,*)abs(ayo)
					! ! ! pause
    ! v=ayo
    ! we=0.0_dp
    ! etha=cmplx(0.0_dp,0.0_dp,dp)
    
	! call dot_product_batch(nn,nvec,r,r,ta_sum)
	! call dot_product_one_time_batch(nn,nvec,r0_initial,r,rho)	
	! call dot_product_batch(nn,nvec,bb,bb,bmag_sum1)	
	! ! ! write(*,*)'d',abs(rho),'dd'
	! ! ! pause

    ! ta=sqrt(abs(ta_sum))
    ! bmag=sqrt(abs(bmag_sum1))
	! do ii =1,nvec
		! rerr(ii)=ta(ii)/bmag(ii)
	! end do
	
    ! iters: do it=1,itmax
	   ! call dot_product_one_time_batch(nn,nvec,r0_initial,v,amgis)
	
	! do ii =1,nvec
		! ahpla(ii)=rho(ii)/amgis(ii)
		! ye(:,ii)=yo(:,ii)-ahpla(ii)*v(:,ii)
	! end do	   
! ! write(*,*)'3'
       ! call SmartMultifly(trans,nn,level_c,rowblock,nvec,ye,aye)
				  
	! do ii =1,nvec
		! d(:,ii)=yo(:,ii)+(we(ii)*we(ii)*etha(ii)/ahpla(ii))*d(:,ii)
		! w(:,ii)=w(:,ii)-ahpla(ii)*ayo(:,ii)
	! end do	
	! call dot_product_batch(nn,nvec,w,w,we_sum)
	! do ii =1,nvec
		! we(ii)=sqrt(abs(we_sum(ii)))/ta(ii)
		! cm(ii)=1.0d0/sqrt(1.0d0+we(ii)*we(ii))
		! ta(ii)=ta(ii)*we(ii)*cm(ii)
		! etha(ii)=ahpla(ii)*cm(ii)*cm(ii)
		! x(:,ii)=x(:,ii)+etha(ii)*d(:,ii)
		! d(:,ii)=ye(:,ii)+(we(ii)*we(ii)*etha(ii)/ahpla(ii))*d(:,ii)
		! w(:,ii)=w(:,ii)-ahpla(ii)*aye(:,ii)
	! end do
	

	
	! call dot_product_batch(nn,nvec,w,w,we_sum)
	! do ii =1,nvec
		! we(ii)=sqrt(abs(we_sum(ii)))/ta(ii)
		! cm(ii)=1.0d0/sqrt(1.0d0+we(ii)*we(ii))
		! ta(ii)=ta(ii)*we(ii)*cm(ii)
		! etha(ii)=ahpla(ii)*cm(ii)*cm(ii)
		! x(:,ii)=x(:,ii)+etha(ii)*d(:,ii)		
	! end do   
	   


       ! !  check if the result has converged.
       ! if (it>=6 .or. maxval(rerr)<1.0_dp*err) then
		  ! ! ! write(*,*)'4'
		  ! call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)		
          ! r=bb-r

		  ! call dot_product_batch(nn,nvec,r,r,rerr_sum)
		  ! do ii =1,nvec
			  ! rerr(ii)=sqrt(abs(rerr_sum(ii)))/bmag(ii)
		  ! end do
		  
		  ! if(isnan(sum(rerr)))then
			! write(*,*)it,rerr_old
			! stop
		  ! end if
		  
          ! ! ! if (myid==main_id) then
		  ! ! if(level_c==2 .or. level_c==1)then
			 ! Ngood = count(rerr<err)
			 ! Nbad = count(rerr>err)
             ! ! ! print*,'# ofiter,error:',it,maxval(rerr),Ngood,Nbad
		  ! ! end if
			 ! ! ! pause
             ! ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          ! ! ! end if
          
          ! if (err > maxval(rerr) .or. Ngood>Nbad*10) then
             ! err=maxval(rerr)
             ! iter=it

             ! ! ! ! if (myid == main_id) then
                ! ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
             ! ! ! ! end if
             
             ! return
          ! endif
       ! end if
       ! !  make preparations for next iteration
	   
		! call dot_product_one_time_batch(nn,nvec,r0_initial,w,dum)	   
		! do ii =1,nvec
			! beta(ii)=dum(ii)/rho(ii)
			! rho(ii)=dum(ii)
			! yo(:,ii)=w(:,ii)+beta(ii)*ye(:,ii)
		! end do		
		
       ! call SmartMultifly(trans,nn,level_c,rowblock,nvec,yo,ayo)
       
       ! !MAGIC
	   ! do ii =1,nvec
		  ! v(:,ii)=ayo(:,ii)+beta(ii)*( aye(:,ii)+beta(ii)*v(:,ii) )
	   ! end do
	   ! rerr_old = rerr
	   
    ! enddo iters
    ! ! write(*,*)'6'
    ! call SmartMultifly(trans,nn,level_c,rowblock,nvec,x,r)
	   
    ! !MAGIC
    ! r=bb-r
	! call dot_product_batch(nn,nvec,r,r,err_sum)	
	! err = 0
	! do ii =1,nvec	 
		! err=max(err,sqrt(abs(err_sum(ii)))/bmag(ii))
    ! end do
	
	! iter=itmax



   ! print*,'Iterative solver is terminated without convergence!!!',it,err
   ! stop
	
    ! return
  ! end subroutine ztfqmr_batch  
  
  
  
  
  subroutine dot_product_batch(N,Nvec,x,y,dots)
	implicit none 
	integer N,Nvec,ii
	complex(kind=8):: x(N,Nvec),y(N,Nvec),dots(Nvec),x0(N),y0(N)
	do ii =1,Nvec
		x0 = x(:,ii)
		y0 = y(:,ii)
		dots(ii) = dot_product(x0,y0)
	end do
	
  end subroutine dot_product_batch

  subroutine dot_product_one_time_batch(N,Nvec,x,y,dots)
	implicit none 
	integer N,Nvec,ii
	complex(kind=8):: x(N,1),y(N,Nvec),dots(Nvec),x0(N),y0(N)
	x0 = x(:,1)
	do ii =1,Nvec
		y0 = y(:,ii)
		dots(ii) = dot_product(x0,y0)
	end do
	
  end subroutine dot_product_one_time_batch  
  
  subroutine initialize_r0(nn)
    implicit none
    integer::jran,i,nn
    real(kind=dp)::r0_dummy(1:2*nn)
    !if (first_frequency) then
       allocate(r0_initial(1:nn))

       !mem_est=nsuinf(3)*complex_mem/1024.0_dp/1024.0_dp
       !if (myid == main_id) then
       !   print*, 'r0_initial requires (MB):',mem_est
       !   write(20,*) 'r0_initial requires (MB):',mem_est        
       !end if

    !end if
!    jran=1211
    ! ! ! ! ! jran=-3
    ! ! ! ! ! do i=1,2*nn
       ! ! ! ! ! r0_dummy(i)=2.0_dp*real(ran1(jran),dp)-1.0_dp
    ! ! ! ! ! end do
    ! ! ! ! ! do i=1,nn
       ! ! ! ! ! r0_initial(i)= cmplx(r0_dummy(i),r0_dummy(nn+i),dp)
    ! ! ! ! ! end do
	
	do i=1,nn
       r0_initial(i)= random_complex_number()
    end do	
	
    return
  end subroutine initialize_r0


  function ran1(idum)
    implicit none
    integer,intent(inout)::idum
    integer,PARAMETER::IM=2147483647,IQ=127773,IR=2836,IA=16807,NTAB=32
    integer::NDIV
    real,parameter::EPS=1.2E-7,RNMX=1.-EPS
    real::am,ran1
    integer::j,k,iv(NTAB)=0,iy=0
    SAVE iv,iy

    NDIV=1+(IM-1)/NTAB
    AM=1./IM

    if (idum<= 0 .or. iy== 0) then
       idum=max(-idum,1)
       do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum .lt. 0) idum=idum+IM
          if (j .le. NTAB) iv(j)=idum
       end do
       iy=iv(1)
    end if
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum .lt. 0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)
    return
  end function ran1


  
 
subroutine MVM_Z_factorized(Ns,Vin,Vout)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer Ns
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8)::Vin(:),Vout(:)
	! complex(kind=8)::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
   
 
    
	level_c = 1
	rowblock = 1
	
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	

	idx_start_glo = 1	 
	   

	! if(num_vectors==1)then
		
		
		num_vectors = 1
		! get the right multiplied vectors
		ctemp1=1.0d0 ; ctemp2=0.0d0
		allocate(vec_old(Ns,num_vectors))
		allocate(vec_new(Ns,num_vectors))
		vec_old(1:Ns,1) = Vin

		! write(*,*)'begin'
		
		do level = Maxlevel_for_blocks+1,1,-1
			N_diag = 2**(level-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			 
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
				if(level==Maxlevel_for_blocks+1)then
					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
					&vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
					! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
				else 
					call Bplus_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
				end if
			end do
			vec_old = vec_new
		end do
		Vout = vec_new(1:Ns,1)
		deallocate(vec_old)
		deallocate(vec_new)	
	! else 
		! write(*,*)'multiple RHS not implemented yet'
	! end if	 
    return                

end subroutine MVM_Z_factorized 
  
  


subroutine MVM_Z_factorized_MRHS(Ns,num_vectors,Vin,Vout)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer Ns
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	! complex(kind=8)::Vin(:),Vout(:)
	complex(kind=8)::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
   
 
    
	level_c = 1
	rowblock = 1
	
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	

	idx_start_glo = 1	 
	   

	! if(num_vectors==1)then
		
		
		! num_vectors = 1
		! get the right multiplied vectors
		ctemp1=1.0d0 ; ctemp2=0.0d0
		! allocate(vec_old(Ns,num_vectors))
		allocate(vec_new(Ns,num_vectors))
		Vout = Vin

		! write(*,*)'begin'
		
		do level = Maxlevel_for_blocks+1,1,-1
			N_diag = 2**(level-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
				if(level==Maxlevel_for_blocks+1)then
					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
					&Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
					! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
				else 
					call Bplus_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
				end if
			end do
			Vout = vec_new
		end do
		! Vout = vec_new(1:Ns,1)
		! deallocate(vec_old)
		deallocate(vec_new)	
	! else 
		! write(*,*)'multiple RHS not implemented yet'
	! end if	 
    return                

end subroutine MVM_Z_factorized_MRHS  

end module Butterfly_inversion
