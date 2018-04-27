module Butterfly_inversion_schur_partition
use Butterfly_inversion

contains 

subroutine Butterfly_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off
    type(matrixblock),pointer::blocks_minusBC
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory

	block_off => cascading_factors(level_c)%matrices_block(rowblock*2-1)			
	block_o => cascading_factors(level_c)%matrices_block_inverse_schur(rowblock)
	block_o%level_butterfly = block_off%level_butterfly	
		
    block_o =>  cascading_factors(level_c)%matrices_block_inverse_schur(rowblock) 
	Memory = 0
	
	allocate(partitioned_blocks(0:block_o%level_butterfly+1))
	do ll=0,block_o%level_butterfly+1
		allocate(partitioned_blocks(ll)%blocks_A)
		allocate(partitioned_blocks(ll)%blocks_B)
		allocate(partitioned_blocks(ll)%blocks_C)
		allocate(partitioned_blocks(ll)%blocks_D)
	end do
	blocks_minusBC => partitioned_blocks(0)%blocks_D
	
	
	call Butterfly_diagonal_minusBC(level_c,rowblock,blocks_minusBC)
	
	call Butterfly_inverse_partitionedinverse_IplusButter(0,1)	


	call copy_delete_butterfly(blocks_minusBC,block_o,Memory)
	! call delete_blocks(blocks_minusBC)		


	do ll=0,block_o%level_butterfly+1
		deallocate(partitioned_blocks(ll)%blocks_A)
		deallocate(partitioned_blocks(ll)%blocks_B)
		deallocate(partitioned_blocks(ll)%blocks_C)
		deallocate(partitioned_blocks(ll)%blocks_D)
	end do
	deallocate(partitioned_blocks)

	
    return

end subroutine Butterfly_inverse_schur_partitionedinverse






recursive subroutine Butterfly_inverse_partitionedinverse_IplusButter(level,ADflag)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock,ADflag
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::blocks_io,blocks_A,blocks_B,blocks_C,blocks_D
    type(matrixblock)::blocks_schur
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:)
	
	if(ADflag==1)then
		blocks_io=> partitioned_blocks(level)%blocks_D
	else 
		blocks_io=> partitioned_blocks(level)%blocks_A
	end if
	
	! write(*,*)'inverse ABDC',blocks_io%row_group,blocks_io%col_group,blocks_io%level,blocks_io%level_butterfly		
	
	if(blocks_io%level_butterfly==0)then
		call Butterfly_inverse_IplusButter_woodbury(blocks_io,Memory)
		return
    else 
		blocks_A => partitioned_blocks(level+1)%blocks_A
		blocks_B => partitioned_blocks(level+1)%blocks_B
		blocks_C => partitioned_blocks(level+1)%blocks_C
		blocks_D => partitioned_blocks(level+1)%blocks_D	
		! split into four smaller butterflies
		call Butterfly_split(blocks_io, blocks_A, blocks_B, blocks_C, blocks_D)
		
		! partitioned inverse of D
		call Butterfly_inverse_partitionedinverse_IplusButter(level+1,1)
		
		! construct the schur complement A-BD^-1C
		if(reducelevel_flag==1)then
			level_butterfly = blocks_A%level_butterfly
		else 
			! level_butterfly = int((Maxlevel-blocks_A%level)/2)*2
			level_butterfly = Maxlevel-blocks_A%level
		end if
		! write(*,*)'A-BDC',level_butterfly,level
		call Butterfly_A_minusBDinvC(level_butterfly,blocks_A, blocks_B, blocks_C, blocks_D,error_inout) 
		error_cnt = error_cnt + 1
		error_avr_glo = error_avr_glo + error_inout
		
		! partitioned inverse of the schur complement 
		call Butterfly_inverse_partitionedinverse_IplusButter(level+1,2)
		
		
		! construct the inverse
		if(reducelevel_flag==1)then
			level_butterfly = blocks_io%level_butterfly
		else 
			level_butterfly = Maxlevel-blocks_io%level
			if(level==0)level_butterfly = int((Maxlevel-blocks_io%level)/2)*2
		end if
		! write(*,*)'inverse ABDC',blocks_io%row_group,blocks_io%col_group,blocks_io%level,blocks_io%level_butterfly		
		call Butterfly_inverse_ABCD(level_butterfly,blocks_io,blocks_A, blocks_B, blocks_C, blocks_D,error_inout) 
		error_cnt = error_cnt + 1
		error_avr_glo = error_avr_glo + error_inout
		! stop
		
		if(level==0 .and. verboselevel>=1)write(*,'(A23,A6,I3,A8,I3,A7,Es14.7)')' SchurI ',' rank:',blocks_io%rankmax,' L_butt:',blocks_io%level_butterfly,' error:',error_inout
			
		call delete_blocks(blocks_A)
		call delete_blocks(blocks_B)
		call delete_blocks(blocks_C)
		call delete_blocks(blocks_D)
		
		return
	
	end if	
end subroutine Butterfly_inverse_partitionedinverse_IplusButter



subroutine Butterfly_diagonal_minusBC(level_c,rowblock,blocks_minusBC)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock)::blocks_minusBC
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:)
	
	Memory = 0
	
    block_o =>  cascading_factors(level_c)%matrices_block_inverse_schur(rowblock) 
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
	blocks_minusBC%level = block_o%level
	blocks_minusBC%level_butterfly = block_o%level_butterfly
	blocks_minusBC%col_group = block_o%col_group
	blocks_minusBC%row_group = block_o%row_group
	blocks_minusBC%style = block_o%style

	if(blocks_minusBC%level_butterfly==0)then
		allocate(blocks_minusBC%ButterflyU(1))
		allocate(blocks_minusBC%ButterflyV(1))

		nn = size(block_off1%ButterflyV(1)%matrix,1)
		mm = size(block_off1%ButterflyU(1)%matrix,1)

		kk1 = size(block_off1%ButterflyV(1)%matrix,2)
		kk2 = size(block_off2%ButterflyV(1)%matrix,2)

		rank = min(kk1,kk2)
		blocks_minusBC%rankmax = rank
		blocks_minusBC%rankmin = rank

		allocate(blocks_minusBC%ButterflyU(1)%matrix(mm,rank))
		allocate(blocks_minusBC%ButterflyV(1)%matrix(mm,rank))
		allocate(matrix_small(kk1,kk2))
		call gemmTN_omp(block_off1%ButterflyV(1)%matrix,block_off2%ButterflyU(1)%matrix,matrix_small,kk1,nn,kk2)
		if(kk1>kk2)then
			call gemm_omp(block_off1%ButterflyU(1)%matrix,matrix_small,blocks_minusBC%ButterflyU(1)%matrix,mm,kk1,rank)
			blocks_minusBC%ButterflyV(1)%matrix = -block_off2%ButterflyV(1)%matrix
		else 
			call gemmNT_omp(block_off2%ButterflyV(1)%matrix,matrix_small,blocks_minusBC%ButterflyV(1)%matrix,mm,kk2,rank)
			blocks_minusBC%ButterflyU(1)%matrix = -block_off1%ButterflyU(1)%matrix
		end if
		deallocate(matrix_small)
	else 
		do tt = 1,10
			do ntry=1,1
			itermax = 0
			n1 = OMP_get_wtime()
			call Initialize_Butterfly_inverse_BC(level_c,rowblock,tt-1)
			n2 = OMP_get_wtime()
			Time_Init_inverse = Time_Init_inverse + n2-n1
			
			n1 = OMP_get_wtime()
			call Reconstruction_LL_BC(level_c,rowblock)	
			call Reconstruction_RR_BC(level_c,rowblock,error_inout)
			n2 = OMP_get_wtime()	
			Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
			
			! write(*,*)tt,error_inout
			
			
			if(error_inout>iter_tolerance)then
			! if(0)then
				call Delete_randomized_butterfly()
			else 
				call delete_blocks(blocks_minusBC)
				call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
				rank_new_max = butterfly_block_randomized(1)%rankmax
				call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),blocks_minusBC,Memory) 
				deallocate(butterfly_block_randomized)
				! call copy_randomizedbutterfly(butterfly_block_randomized(1),blocks_minusBC,Memory) 
				! call Delete_randomized_butterfly()
				if(verboselevel>=1)write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' -BC No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',blocks_minusBC%level_butterfly,' error:',error_inout
					
				
				return			
			end if		
			end do
		end do
		write(*,*)'randomized scheme not converged in BC'
		stop
	end if	
    return

end subroutine Butterfly_diagonal_minusBC



subroutine Reconstruction_LL_BC(level_c,rowblock)
    
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
	integer niter,unique_nth
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
			call Get_Randomized_Vectors_LL_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
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
    
end subroutine Reconstruction_LL_BC



subroutine Reconstruction_RR_BC(level_c,rowblock,error)
    
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
			call Get_Randomized_Vectors_RR_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_BC(level_c,rowblock,error)


	
    return
    
end subroutine Reconstruction_RR_BC




subroutine Test_Error_RR_BC(level_c,rowblock,error)

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
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:),Id(:,:),Vd(:,:)
	type(matrixblock),pointer::block_o
	
	
	block_o =>  cascading_factors(level_c)%matrices_block_inverse_schur(rowblock) 
	level_butterfly=block_o%level_butterfly
	num_blocks=2**level_butterfly
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call Get_Randomized_Vectors_RR_Test_BC(level_c,rowblock,num_vect)
	
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


	! allocate (Id(mm,num_vect))
	! allocate (Vd(mm,num_vect))
	! k=0
	! do i=1, num_blocks
		! dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		! do ii=1, dimension_m
			! do jj=1, num_vect
				! Id(ii+k,jj)=random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)
			! enddo
		! enddo
		! !$omp end parallel do
		! k=k+dimension_m
	! enddo 		
	
	! call SmartMultiflySchur('N',mm,level_c,rowblock,1,Id,Vd)
	! Vd = Vd - Id
	! write(888,*)abs(Vd)
	! write(889,*)abs(RandomVectors_Output_ref)
	! stop
	
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

end subroutine Test_Error_RR_BC




subroutine Get_Randomized_Vectors_LL_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	
	block_o => cascading_factors(level_c)%matrices_block_inverse_schur(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
    level_butterfly=block_o%level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_off1%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
    groupm=block_off1%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 


	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(block_o%level_butterfly)
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
				
				
				
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) !random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(block_off1,'T',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)	
	call butterfly_block_MVP_randomized_dat(block_off2,'T',nn,mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	RandomVectors_InOutput(3)%vector = -RandomVectors_InOutput(3)%vector
	

	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
	k=0
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 	

    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_LL_BC




subroutine Get_Randomized_Vectors_RR_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	type(RandomBlock), pointer :: random
	
	
	block_o => cascading_factors(level_c)%matrices_block_inverse_schur(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
    level_butterfly=block_o%level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_off1%col_group  ! Note: row_group and col_group interchanged here   
    groupm=block_off1%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	

	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(block_o%level_butterfly)
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
				
				
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(block_off2,'N',nn,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)	
	call butterfly_block_MVP_randomized_dat(block_off1,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	RandomVectors_InOutput(3)%vector = -RandomVectors_InOutput(3)%vector	
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_BC



subroutine Get_Randomized_Vectors_RR_Test_BC(level_c,rowblock,num_vect_sub)

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
	type(matrixblock),pointer::block_o,block_off1,block_off2
	
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
	
	
	block_o => cascading_factors(level_c)%matrices_block_inverse_schur(rowblock)
	block_off1 => cascading_factors(level_c)%matrices_block(rowblock*2-1)	
	block_off2 => cascading_factors(level_c)%matrices_block(rowblock*2)	
	
    level_butterfly=block_o%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_off1%col_group  ! Note: row_group and col_group interchanged here   
    groupm=block_off1%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(block_o%level_butterfly)
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
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(block_off2,'N',nn,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)	
	call butterfly_block_MVP_randomized_dat(block_off1,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	RandomVectors_InOutput(3)%vector = -RandomVectors_InOutput(3)%vector
	
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

end subroutine Get_Randomized_Vectors_RR_Test_BC



subroutine Initialize_Butterfly_inverse_BC(level_c,rowblock,kover)

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
    
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	! write(*,*)level_c+1,rowblock*2-1
	
	allocate (butterfly_block_randomized(1))
    
    level_butterfly=block_o%level_butterfly
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    
    num_blocks=2**level_butterfly
 	
	! dimension_rank= max(block_off1%rankmax,block_off2%rankmax)+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= max(block_off1%rankmax,block_off2%rankmax) *1.2d0**(kover) 
	
	! if(level_c==2)dimension_rank=11
	! if(level_c==1)dimension_rank=9+kover
	
	! write(*,*)dimension_rank
	
    groupm=block_o%row_group   ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    butterfly_block_randomized(1)%dimension_rank=dimension_rank

    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))


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
		
		dimension_m=size(block_off1%ButterflyU(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        
		! !$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
		dimension_n=size(block_off2%ButterflyV(blocks)%matrix,1)

        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        
		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
		
		! !$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
    enddo

    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
        allocate (butterfly_block_randomized(1)%ButterflyKerl(level_butterfly))

        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))

        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine Initialize_Butterfly_inverse_BC




subroutine Butterfly_A_minusBDinvC(level_butterfly,blocks_A, blocks_B, blocks_C, blocks_D,error_inout) 

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,r1,r2,r3,r3tmp,mn,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock)::blocks_A, blocks_B, blocks_C, blocks_D, blocks_Schur
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:),U1(:,:),V1(:,:),U2(:,:),V2(:,:),U3(:,:),V3(:,:),U3tmp(:,:),V3tmp(:,:),UUtmp(:,:),VVtmp(:,:),UU(:,:),VV(:,:),UUr(:,:),VVr(:,:)
	real*8,allocatable :: Singular(:)
	complex (kind=8), allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:),Vout3(:,:),Vout4(:,:),Vout(:,:)
	complex (kind=8)::ctemp1,ctemp2
	
	ctemp1 = 1d0; ctemp2 = 0d0
	Memory = 0
	
	if(blocks_A%level_butterfly==0)then
		
		n1 = OMP_get_wtime()					  
		mm=size(blocks_B%ButterflyU(1)%matrix,1)
		nn=size(blocks_B%ButterflyV(1)%matrix,1)
		
		!  copy A
		r1 = size(blocks_A%ButterflyU(1)%matrix,2)
		allocate(U1(mm,r1))
		allocate(V1(mm,r1))
		U1 = blocks_A%ButterflyU(1)%matrix
		V1 = blocks_A%ButterflyV(1)%matrix
		
		
		!  construct BC
		kk1 = size(blocks_B%ButterflyV(1)%matrix,2)
		kk2 = size(blocks_C%ButterflyV(1)%matrix,2)		
		r2 = kk2
		allocate(U2(mm,r2))
		allocate(V2(mm,r2))		
		allocate(matrix_small(kk1,kk2))
		call gemmTN_omp(blocks_B%ButterflyV(1)%matrix,blocks_C%ButterflyU(1)%matrix,matrix_small,kk1,nn,kk2)
		call gemm_omp(blocks_B%ButterflyU(1)%matrix,matrix_small,U2,mm,kk1,r2)
		V2 = -blocks_C%ButterflyV(1)%matrix
		deallocate(matrix_small)
		

		!  construct B(Dinv-I)C
		kk1 = size(blocks_D%ButterflyV(1)%matrix,2)
		kk2 = size(blocks_C%ButterflyV(1)%matrix,2)		
		r3tmp = kk2
		allocate(U3tmp(nn,r3tmp))
		allocate(V3tmp(mm,r3tmp))		
		allocate(matrix_small(kk1,kk2))
		call gemmTN_omp(blocks_D%ButterflyV(1)%matrix,blocks_C%ButterflyU(1)%matrix,matrix_small,kk1,nn,kk2)
		call gemm_omp(blocks_D%ButterflyU(1)%matrix,matrix_small,U3tmp,nn,kk1,r3tmp)
		V3tmp = blocks_C%ButterflyV(1)%matrix
		deallocate(matrix_small)
		
		kk1 = size(blocks_B%ButterflyV(1)%matrix,2)
		kk2 = size(V3tmp,2)		
		r3 = kk2
		allocate(U3(mm,r3))
		allocate(V3(mm,r3))		
		allocate(matrix_small(kk1,kk2))
		call gemmTN_omp(blocks_B%ButterflyV(1)%matrix,U3tmp,matrix_small,kk1,nn,kk2)
		call gemm_omp(blocks_B%ButterflyU(1)%matrix,matrix_small,U3,mm,kk1,r3)
		V3 = -V3tmp
		deallocate(matrix_small)		
		
		! constrcut the schur complement 
		allocate(UUtmp(mm,r1+r2+r3))
		UUtmp(1:mm,1:r1) = U1
		UUtmp(1:mm,1+r1:r1+r2) = U2
		UUtmp(1:mm,1+r1+r2:r1+r2+r3) = U3

		allocate(VVtmp(mm,r1+r2+r3))
		VVtmp(1:mm,1:r1) = V1
		VVtmp(1:mm,1+r1:r1+r2) = V2
		VVtmp(1:mm,1+r1+r2:r1+r2+r3) = V3		
		

		mn=min(mm,r1+r2+r3)
		allocate (UU(mm,mn),VV(mn,r1+r2+r3),Singular(mn))
		call SVD_Truncate(UUtmp,mm,r1+r2+r3,mn,UU,VV,Singular,SVD_tolerance_forward*1d0,rank)			
		! rank = 18
		! write(*,*)'hahahahahahaha',rank,blocks_A%rankmax,blocks_B%rankmax,blocks_C%rankmax,blocks_D%rankmax,r1+r2+r3
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1						 
		
		! if the low rank scheme increases the rank, use butterfly instead 
		if(rank>r1 .and. rank>r2 .and. rank>r3 .and. level_butterfly/=0)then
			deallocate(U1,V1,U2,V2,U3,V3,U3tmp,V3tmp,UUtmp,VVtmp,UU,VV,Singular)			
			
			do tt = 1,10
				do ntry=1,1
				itermax = 0
				n1 = OMP_get_wtime()
				call Initialize_Butterfly_A_minusBDinvC(level_butterfly,blocks_A, blocks_B, blocks_C, blocks_D,tt-1)
				n2 = OMP_get_wtime()
				Time_Init_inverse = Time_Init_inverse + n2-n1
				
				n1 = OMP_get_wtime()
				call Reconstruction_LL_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D)	
				call Reconstruction_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,error_inout)
				n2 = OMP_get_wtime()	
				Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
				
				! write(*,*)tt,error_inout
				
				
				if(error_inout>iter_tolerance)then
				! if(0)then
					call Delete_randomized_butterfly()
				else 
					call delete_blocks(blocks_A)
					call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
					rank_new_max = butterfly_block_randomized(1)%rankmax
					call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),blocks_A,Memory) 
					deallocate(butterfly_block_randomized)
					! call copy_randomizedbutterfly(butterfly_block_randomized(1),blocks_A,Memory) 
					! call Delete_randomized_butterfly()
					if(verboselevel>=2)write(*,'(A35,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' A-BD^-1C ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',blocks_A%level_butterfly,' error:',error_inout
					
					return			
				end if		
				end do
			end do
			write(*,*)'randomized scheme not converged in A_minusBDinvC'
			stop			
			
		else 		
			! ! if(level_butterfly/=0)write(*,*)'aha this is a bug man!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'																					  
			allocate(UUr(mm,rank))
			allocate(VVr(rank,r1+r2+r3))
			do ii =1,rank
				UUr(:,ii) = UU(:,ii)*Singular(ii)
			end do
			VVr = VV(1:rank,:)
			
			
			! ! allocate(Vin(mm,1))		
			! ! Vin = 1
			! ! allocate(Vout1(nn,1))	
			! ! allocate(Vout2(nn,1))	
			! ! allocate(Vout3(mm,1))	
			! ! allocate(Vout4(mm,1))	
			! ! call butterfly_block_MVP_randomized_dat(blocks_C,'N',nn,mm,1,vin,vout1,ctemp1,ctemp2) 
			! ! call butterfly_block_MVP_randomized_dat(blocks_D,'N',nn,nn,1,vout1,vout2,ctemp1,ctemp2) 
			! ! vout2 = vout1 + vout2
			! ! call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm,nn,1,vout2,vout3,ctemp1,ctemp2) 
			! ! call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,1,vin,vout4,ctemp1,ctemp2) 		
			! ! vout4 = vout4 - vout3
			
			
			call delete_blocks(blocks_A)
			allocate(blocks_A%ButterflyU(1))	
			allocate(blocks_A%ButterflyU(1)%matrix(mm,rank))	
			blocks_A%ButterflyU(1)%matrix = UUr
			allocate(blocks_A%ButterflyV(1))	
			allocate(blocks_A%ButterflyV(1)%matrix(mm,rank))
			call gemmNT_omp(VVtmp,VVr,blocks_A%ButterflyV(1)%matrix,mm,r1+r2+r3,rank)
			blocks_A%rankmax = rank
			blocks_A%rankmin = rank

			deallocate(U1,V1,U2,V2,U3,V3,U3tmp,V3tmp,UUtmp,VVtmp,UU,VV,UUr,VVr,Singular)
			
			
			! ! allocate(Vout(mm,1))	
			! ! call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,1,vin,vout,ctemp1,ctemp2) 				
			! ! write(*,*)r1,r2,r3,rank,'LR schur error: ',fnorm(vout-vout4,mm,1)/fnorm(vout4,mm,1)		
			error_inout = SVD_tolerance_forward
		end if		

	else 
		do tt = 1,10
			do ntry=1,1
			itermax = 0
			n1 = OMP_get_wtime()
			call Initialize_Butterfly_A_minusBDinvC(level_butterfly,blocks_A, blocks_B, blocks_C, blocks_D,tt-1)
			n2 = OMP_get_wtime()
			Time_Init_inverse = Time_Init_inverse + n2-n1
			
			n1 = OMP_get_wtime()
			call Reconstruction_LL_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D)	
			call Reconstruction_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,error_inout)
			n2 = OMP_get_wtime()	
			Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
			
			! write(*,*)tt,error_inout
			
			
			if(error_inout>iter_tolerance)then
			! if(0)then
				call Delete_randomized_butterfly()
			else 
				call delete_blocks(blocks_A)
				call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
				rank_new_max = butterfly_block_randomized(1)%rankmax
				call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),blocks_A,Memory) 
				deallocate(butterfly_block_randomized)
				! call copy_randomizedbutterfly(butterfly_block_randomized(1),blocks_A,Memory) 
				! call Delete_randomized_butterfly()
				if(verboselevel>=2)write(*,'(A35,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' A-BD^-1C ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',blocks_A%level_butterfly,' error:',error_inout
				
				return			
			end if		
			end do
		end do
		write(*,*)'randomized scheme not converged in A_minusBDinvC'
		stop
	end if	
    return

end subroutine Butterfly_A_minusBDinvC


subroutine Reconstruction_LL_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D)
    
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
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	real*8:: error_inout
	integer,allocatable::perms(:)
	
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
			call Get_Randomized_Vectors_LL_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
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
    
end subroutine Reconstruction_LL_A_minusBDinvC





subroutine Reconstruction_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,error)
    
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
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	
	level_butterfly=butterfly_block_randomized(1)%level_butterfly
		
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
			call Get_Randomized_Vectors_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,error)


	
    return
    
end subroutine Reconstruction_RR_A_minusBDinvC




subroutine Test_Error_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,error)

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
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:),Id(:,:),Vd(:,:)
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
	

	level_butterfly=butterfly_block_randomized(1)%level_butterfly
	num_blocks=2**level_butterfly
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=blocks_A%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call Get_Randomized_Vectors_RR_Test_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,num_vect)
	
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


	! allocate (Id(mm,num_vect))
	! allocate (Vd(mm,num_vect))
	! k=0
	! do i=1, num_blocks
		! dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		! do ii=1, dimension_m
			! do jj=1, num_vect
				! Id(ii+k,jj)=random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)
			! enddo
		! enddo
		! !$omp end parallel do
		! k=k+dimension_m
	! enddo 		
	
	! call SmartMultiflySchur('N',mm,level_c,rowblock,1,Id,Vd)
	! Vd = Vd - Id
	! write(888,*)abs(Vd)
	! write(889,*)abs(RandomVectors_Output_ref)
	! stop
	
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

end subroutine Test_Error_RR_A_minusBDinvC



subroutine Get_Randomized_Vectors_LL_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)

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
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
	
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
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(4))

    groupn=blocks_B%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
    groupm=blocks_B%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 


	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(4)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used for output 
	
	do ii =1,4
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
				
				
				
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) !random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_B,'T',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)
	if(isnan(fnorm(RandomVectors_InOutput(2)%vector,nn,num_vect_sub)))then
		write(*,*)fnorm(RandomVectors_InOutput(1)%vector,mm,num_vect_sub),fnorm(RandomVectors_InOutput(2)%vector,nn,num_vect_sub),'1111'
		stop
	end if
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_D,'T',nn,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,ctemp1,ctemp2)
	if(isnan(fnorm(RandomVectors_InOutput(4)%vector,nn,num_vect_sub)))then
		write(*,*)fnorm(RandomVectors_InOutput(2)%vector,nn,num_vect_sub),fnorm(RandomVectors_InOutput(4)%vector,nn,num_vect_sub),'2222'
		stop
	end if	
	RandomVectors_InOutput(4)%vector = RandomVectors_InOutput(4)%vector + RandomVectors_InOutput(2)%vector
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_C,'T',nn,mm,num_vect_sub,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
	if(isnan(fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub)))then
		write(*,*)fnorm(RandomVectors_InOutput(4)%vector,nn,num_vect_sub),fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub),'3333'
		stop
	end if		
	ctemp1=1.0d0 ; ctemp2=-1.0d0
	call butterfly_block_MVP_randomized_dat(blocks_A,'T',mm,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
	if(isnan(fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub)))then
		write(*,*)fnorm(RandomVectors_InOutput(1)%vector,mm,num_vect_sub),fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub),'4444'
		stop
	end if


	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
	k=0
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! write(*,*)num_blocks,i,nn,shape(random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix),blocks_A%level_butterfly,blocks_B%level_butterfly,butterfly_block_randomized(1)%level_butterfly
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 	

    ! !$omp parallel do default(shared) private(i)
    do i=1, 4
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_LL_A_minusBDinvC


subroutine Get_Randomized_Vectors_RR_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)

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
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
	
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
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	type(RandomBlock), pointer :: random
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(4))

    groupn=blocks_B%col_group  ! Note: row_group and col_group interchanged here   
    groupm=blocks_B%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	

	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(4)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used as output 
	do ii =1,4
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


				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_C,'N',nn,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_D,'N',nn,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,ctemp1,ctemp2)		
	RandomVectors_InOutput(4)%vector = RandomVectors_InOutput(4)%vector + RandomVectors_InOutput(2)%vector
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	ctemp1=1.0d0 ; ctemp2=-1.0d0
	call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 4
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_A_minusBDinvC



subroutine Get_Randomized_Vectors_RR_Test_A_minusBDinvC(blocks_A,blocks_B,blocks_C,blocks_D,num_vect_sub)

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
	type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
	
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
	

    level_butterfly=butterfly_block_randomized(1)%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(4))

    groupn=blocks_B%col_group  ! Note: row_group and col_group interchanged here   
    groupm=blocks_B%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
    allocate (RandomVectors_InOutput(4)%vector(nn,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,4
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
		! !$omp parallel do default(shared) private(ii,jj)
		 ! do jj=1,num_vect_sub
			 ! do ii=1, mm
				 ! RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	
			 ! enddo
		 ! enddo
		 ! !$omp end parallel do
		 
		 call RandomMat(mm,num_vect_sub,min(mm,num_vect_sub),RandomVectors_InOutput(1)%vector(1+k:mm+k,1:num_vect_sub),0)						
	end do

	! get the right multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_C,'N',nn,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_D,'N',nn,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(4)%vector,ctemp1,ctemp2)		
	RandomVectors_InOutput(4)%vector = RandomVectors_InOutput(4)%vector + RandomVectors_InOutput(2)%vector
	ctemp1=1.0d0 ; ctemp2=0.0d0
	call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(4)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)	
	ctemp1=1.0d0 ; ctemp2=-1.0d0
	call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
	
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
    do i=1, 4
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_Test_A_minusBDinvC


subroutine Initialize_Butterfly_A_minusBDinvC(level_butterfly,blocks_A, blocks_B, blocks_C, blocks_D,kover)

    use MODULE_FILE
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max,dimension_rank, dimension_m, dimension_n, blocks, groupm, groupm_start,groupn,index_j,index_i
    real*8 a,b,c,d
    complex (kind=8) ctemp
	type(matrixblock)::blocks_A, blocks_B, blocks_C, blocks_D
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer rankmax
    
	allocate (butterfly_block_randomized(1))
    
    ! level_butterfly=blocks_A%level_butterfly
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    
    num_blocks=2**level_butterfly
 	
	rankmax = -1000
	rankmax = max(rankmax,blocks_A%rankmax)
	rankmax = max(rankmax,blocks_B%rankmax)
	rankmax = max(rankmax,blocks_C%rankmax)
	rankmax = max(rankmax,blocks_D%rankmax)
	
																					  
 
	! dimension_rank= rankmax+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= rankmax *1.2d0**(kover) 
	
	! if(level_c==2)dimension_rank=11
	! if(level_c==1)dimension_rank=9+kover
	
	! write(*,*)dimension_rank
	
    groupm=blocks_A%row_group   ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    butterfly_block_randomized(1)%dimension_rank=dimension_rank

	groupm_start=groupm*2**level_butterfly

	
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))

	dimension_max = 2*dimension_rank
	do blocks=1, num_blocks	
		dimension_m=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
		dimension_n=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
		dimension_max = max(dimension_max,dimension_m)	
		dimension_max = max(dimension_max,dimension_n)	
	end do	
	allocate(butterfly_block_randomized(1)%KerInv(dimension_max,dimension_max))
	call RandomMat(dimension_max,dimension_max,dimension_max,butterfly_block_randomized(1)%KerInv,3)	
    
    do blocks=1, num_blocks
		
		dimension_m= basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        
		! !$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
		dimension_n=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
	! write(*,*)dimension_m,dimension_n,'hahah'
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        
		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
		
		! !$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
    enddo

    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
        allocate (butterfly_block_randomized(1)%ButterflyKerl(level_butterfly))

        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))

        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine Initialize_Butterfly_A_minusBDinvC


subroutine Butterfly_inverse_ABCD(level_butterfly,blocks_o,blocks_A, blocks_B, blocks_C, blocks_D,error_inout) 

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	use Butterfly_compress_forward
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,r1,r2,r3,r3tmp,mn,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock)::blocks_o,blocks_A, blocks_B, blocks_C, blocks_D,block_tmp
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:),U1(:,:),V1(:,:),U2(:,:),V2(:,:),U3(:,:),V3(:,:),U3tmp(:,:),V3tmp(:,:),UUtmp(:,:),VVtmp(:,:),UU(:,:),VV(:,:),UUr(:,:),VVr(:,:)
	real*8,allocatable :: Singular(:)
	complex (kind=8), allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:),Vout3(:,:),Vout4(:,:),Vout(:,:)
	complex (kind=8)::ctemp1,ctemp2
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matsub_tmp(:,:)
	integer idx_start_m_ref
	
	
	
	ctemp1 = 1d0; ctemp2 = 0d0
	Memory = 0
	
	! write(*,*)blocks_A%level_butterfly,blocks_A%rankmax,blocks_B%level_butterfly,blocks_B%rankmax,blocks_C%level_butterfly,blocks_C%rankmax,blocks_D%level_butterfly,blocks_D%rankmax
	! write(*,*)'AD group ',blocks_A%row_group,blocks_D%row_group
		
		
	! mm = basis_group(blocks_A%row_group)%tail - basis_group(blocks_A%row_group)%head + 1
		! ! write(*,*)'e',mm
	! allocate(matin(mm,mm))
	! allocate(matout(mm,mm))
		! ! write(*,*)'f'
	! matin = 0
	! do ii=1,mm
		! matin(ii,ii) = 1
	! end do
	! call Butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,mm,matin,matout,ctemp1,ctemp2)
	! write(*,*)'A',fnorm(matout,mm,mm)
	! deallocate(matin)
	! deallocate(matout)	

		
	! mm = basis_group(blocks_D%row_group)%tail - basis_group(blocks_D%row_group)%head + 1
		! ! write(*,*)'e',mm
	! allocate(matin(mm,mm))
	! allocate(matout(mm,mm))
		! ! write(*,*)'f'
	! matin = 0
	! do ii=1,mm
		! matin(ii,ii) = 1
	! end do
	! call Butterfly_block_MVP_randomized_dat(blocks_D,'N',mm,mm,mm,matin,matout,ctemp1,ctemp2)
	! write(*,*)'D',fnorm(matout,mm,mm)
	! deallocate(matin)
	! deallocate(matout)			
		
	

	! mm = basis_group(blocks_o%row_group)%tail - basis_group(blocks_o%row_group)%head + 1
		! ! write(*,*)'e',mm
	! allocate(matin(mm,mm))
	! ! allocate(matout(mm,mm))
	! allocate(matsub_glo(mm,mm))
	! matin = 0
	! do ii=1,mm
		! matin(ii,ii) = 1
	! end do
	
	! call butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,'N',mm,mm,matin,matsub_glo)
	! matsub_glo = matsub_glo - matin		
		
	! block_tmp%row_group = blocks_o%row_group
	! block_tmp%col_group = blocks_o%col_group
	! block_tmp%level = blocks_o%level
	! block_tmp%style = blocks_o%style
	
	! idx_start_m_ref = basis_group(blocks_o%row_group)%head

	! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_m_ref)
	
	! write(*,*)'wocaonimaABCD',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,mm)

	! deallocate(matin)
	! deallocate(matsub_glo)
	! call delete_blocks(block_tmp)	
	
	
	
	! mm = basis_group(blocks_o%row_group)%tail - basis_group(blocks_o%row_group)%head + 1
		! ! write(*,*)'e',mm
	! allocate(matin(mm,mm))
	! ! allocate(matout(mm,mm))
	! allocate(matsub_glo(mm,mm))
	! allocate(matsub_tmp(mm,mm))
	! matin = 0
	! do ii=1,mm
		! matin(ii,ii) = 1
	! end do
	
	! call butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,'T',mm,mm,matin,matsub_tmp)
	! matsub_tmp = matsub_tmp - matin	
	! call copymatT_omp(matsub_tmp,matsub_glo,mm,mm)
	! deallocate(matsub_tmp)
	
	! block_tmp%row_group = blocks_o%row_group
	! block_tmp%col_group = blocks_o%col_group
	! block_tmp%level = blocks_o%level
	! block_tmp%style = blocks_o%style
	
	! idx_start_m_ref = basis_group(blocks_o%row_group)%head

	! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_m_ref)
	
	! write(*,*)'wocaonimaABCD',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,mm)

	! deallocate(matin)
	! deallocate(matsub_glo)
	! call delete_blocks(block_tmp)		
	
	
	do tt = 1,10
		do ntry=1,1
		itermax = 0
		n1 = OMP_get_wtime()
		call Initialize_Butterfly_inverse_ABCD(level_butterfly,blocks_o,blocks_A, blocks_B, blocks_C, blocks_D,(tt-1))
		n2 = OMP_get_wtime()
		Time_Init_inverse = Time_Init_inverse + n2-n1
		
		n1 = OMP_get_wtime()

		call Reconstruction_LL_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D)	
		call Reconstruction_RR_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,error_inout)
		n2 = OMP_get_wtime()	
		Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
		
		call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
		write(*,*)tt,error_inout,butterfly_block_randomized(1)%dimension_rank,butterfly_block_randomized(1)%rankmax,'ABCD'
		if(error_inout>iter_tolerance)then
		! if(0)then
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			rank_new_max = butterfly_block_randomized(1)%rankmax
			call Delete_randomized_butterfly()
		else 
			call delete_blocks(blocks_o)
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			rank_new_max = butterfly_block_randomized(1)%rankmax
			call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),blocks_o,Memory) 
			deallocate(butterfly_block_randomized)
			! call copy_randomizedbutterfly(butterfly_block_randomized(1),blocks_o,Memory) 
			! call Delete_randomized_butterfly()
			if(verboselevel>=2)write(*,'(A38,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' ABCDinverse ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',blocks_o%level_butterfly,' error:',error_inout
			
			
			return			
		end if		
		end do
	end do
	write(*,*)'randomized scheme not converged in inverse_ABCD. level: ',blocks_o%level_butterfly,error_inout,rank_new_max
	stop

    return

end subroutine Butterfly_inverse_ABCD


subroutine Reconstruction_LL_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D)
    
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
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	real*8:: error_inout
	integer,allocatable::perms(:)
	
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
			call Get_Randomized_Vectors_LL_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
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
    
end subroutine Reconstruction_LL_inverse_ABCD






subroutine Reconstruction_RR_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,error)
    
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
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	
	level_butterfly=butterfly_block_randomized(1)%level_butterfly
		
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
			call Get_Randomized_Vectors_RR_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,error)


	
    return
    
end subroutine Reconstruction_RR_inverse_ABCD






subroutine Test_Error_RR_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,error)

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
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:),Id(:,:),Vd(:,:)
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
	

	level_butterfly=butterfly_block_randomized(1)%level_butterfly
	num_blocks=2**level_butterfly
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=blocks_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call Get_Randomized_Vectors_RR_Test_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,num_vect)
	
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


	! allocate (Id(mm,num_vect))
	! allocate (Vd(mm,num_vect))
	! k=0
	! do i=1, num_blocks
		! dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		! do ii=1, dimension_m
			! do jj=1, num_vect
				! Id(ii+k,jj)=random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)
			! enddo
		! enddo
		! !$omp end parallel do
		! k=k+dimension_m
	! enddo 		
	
	! call SmartMultiflySchur('N',mm,level_c,rowblock,1,Id,Vd)
	! Vd = Vd - Id
	! write(888,*)abs(Vd)
	! write(889,*)abs(RandomVectors_Output_ref)
	! stop
	
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

end subroutine Test_Error_RR_inverseABCD






subroutine Get_Randomized_Vectors_LL_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)

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
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
	
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
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=blocks_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
    groupm=blocks_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 


	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used for output 
	
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
				
				
				
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) !random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	call butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,'T',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(3)%vector - RandomVectors_InOutput(1)%vector
	
	! write(*,*)'aha',fnorm(RandomVectors_InOutput(1)%vector,mm,num_vect_sub),fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub)
	
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
	k=0
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 	

    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_LL_inverseABCD





subroutine Get_Randomized_Vectors_RR_inverseABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,nth_s,nth_e,num_vect_sub,unique_nth)

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
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
	
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
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	type(RandomBlock), pointer :: random
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=blocks_o%col_group  ! Note: row_group and col_group interchanged here   
    groupm=blocks_o%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	

	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used as output 
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(butterfly_block_randomized(1)%level_butterfly)
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
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! random_complex_number()	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	
	call butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,'N',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(3)%vector - RandomVectors_InOutput(1)%vector	
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_inverseABCD







subroutine Get_Randomized_Vectors_RR_Test_inverse_ABCD(blocks_o,blocks_A,blocks_B,blocks_C,blocks_D,num_vect_sub)

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
	type(matrixblock)::blocks_o,blocks_A,blocks_B,blocks_C,blocks_D
	
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
	

    level_butterfly=butterfly_block_randomized(1)%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=blocks_o%col_group  ! Note: row_group and col_group interchanged here   
    groupm=blocks_o%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

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
		! !$omp parallel do default(shared) private(ii,jj)
		 ! do jj=1,num_vect_sub
			 ! do ii=1, mm
				 ! RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	
			 ! enddo
		 ! enddo
		 ! !$omp end parallel do
		 call RandomMat(mm,num_vect_sub,min(mm,num_vect_sub),RandomVectors_InOutput(1)%vector(1+k:mm+k,1:num_vect_sub),0)
	end do

	! get the right multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	call butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,'N',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(3)%vector - RandomVectors_InOutput(1)%vector
	
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

end subroutine Get_Randomized_Vectors_RR_Test_inverse_ABCD






subroutine Initialize_Butterfly_inverse_ABCD(level_butterfly,blocks_o,blocks_A, blocks_B, blocks_C, blocks_D,kover)

    use MODULE_FILE
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,groupm_start,index_j,index_i
    real*8 a,b,c,d
    complex (kind=8) ctemp
	type(matrixblock)::blocks_o,blocks_A, blocks_B, blocks_C, blocks_D
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer rankmax
    
	allocate (butterfly_block_randomized(1))
    
    ! level_butterfly=blocks_o%level_butterfly
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    
    num_blocks=2**level_butterfly
 	
	rankmax = -1000
	rankmax = max(rankmax,blocks_A%rankmax)
	rankmax = max(rankmax,blocks_B%rankmax)
	rankmax = max(rankmax,blocks_C%rankmax)
	rankmax = max(rankmax,blocks_D%rankmax)
	
																					  
 
	! dimension_rank= rankmax+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= rankmax *1.2d0**(kover)
	 
	! if(level_c==2)dimension_rank=11
	! if(level_c==1)dimension_rank=9+kover
	
	! write(*,*)dimension_rank
	
    groupm=blocks_o%row_group   ! Note: row_group and col_group interchanged here   
	groupm_start=groupm*2**level_butterfly
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    butterfly_block_randomized(1)%dimension_rank=dimension_rank

	
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))

	dimension_max = 2*dimension_rank
	do blocks=1, num_blocks	
		dimension_m=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1        
		dimension_n=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1     
		dimension_max = max(dimension_max,dimension_m)	
		dimension_max = max(dimension_max,dimension_n)	
	end do	
	allocate(butterfly_block_randomized(1)%KerInv(dimension_max,dimension_max))
	call RandomMat(dimension_max,dimension_max,dimension_max,butterfly_block_randomized(1)%KerInv,3)
	
    do blocks=1, num_blocks
		
		! dimension_m=size(blocks_o%ButterflyU(blocks)%matrix,1)
		dimension_m= basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1        
		allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        
		! !$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
		! dimension_n=size(blocks_o%ButterflyV(blocks)%matrix,1)
		dimension_n= basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1     
        
		allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        
		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
		
		! !$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
    enddo

    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
        allocate (butterfly_block_randomized(1)%ButterflyKerl(level_butterfly))
		
        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))
        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine Initialize_Butterfly_inverse_ABCD


subroutine Butterfly_inverse_IplusButter_woodbury(block_o,Memory)

    use MODULE_FILE
    ! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock,kover,rank,kk1,kk2
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,index_j,index_i
    real*8 a,b,c,d,Memory
    complex (kind=8) ctemp
	type(matrixblock)::block_o
	complex (kind=8), allocatable::matrixtemp1(:,:),matrixtemp2(:,:),matrixtemp3(:,:),UU(:,:),VV(:,:),matrix_small(:,:),matrix_eye(:,:),vin(:,:),vout1(:,:),vout2(:,:),vout3(:,:)
	real*8, allocatable:: Singular(:)
    integer, allocatable :: ipiv(:)
	

	nn = size(block_o%ButterflyV(1)%matrix,1)
	rank = size(block_o%ButterflyV(1)%matrix,2)

	call assert(size(block_o%ButterflyV(1)%matrix,2)==size(block_o%ButterflyU(1)%matrix,2),'rank not correct in woodbury')

	allocate(matrix_small(rank,rank))
	allocate(matrix_eye(rank,rank))
	matrix_eye = 0
	do ii=1,rank
	matrix_eye(ii,ii) = 1
	end do	
	call gemmTN_omp(block_o%ButterflyV(1)%matrix,block_o%ButterflyU(1)%matrix,matrix_small,rank,nn,rank)
	matrix_small = matrix_eye+matrix_small
	allocate(ipiv(rank))
	call getrff90(matrix_small,ipiv)
	call getrif90(matrix_small,ipiv)	
	deallocate(ipiv)		
		
	allocate(matrixtemp1(nn,rank))	
		
		
	call gemm_omp(block_o%ButterflyU(1)%matrix,matrix_small,matrixtemp1,nn,rank,rank)	
	block_o%ButterflyU(1)%matrix = 	-matrixtemp1
		
	deallocate(matrixtemp1,matrix_small,matrix_eye)
		
	Memory = 0	
	Memory = Memory + SIZEOF(block_o%ButterflyV(1)%matrix)/1024.0d3
	Memory = Memory + SIZEOF(block_o%ButterflyU(1)%matrix)/1024.0d3
	
	! error_cnt = error_cnt+1
	
    return

end subroutine Butterfly_inverse_IplusButter_woodbury




subroutine Butterfly_split(blocks_i,blocks_A,blocks_B,blocks_C,blocks_D)
    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
	integer level_p,ADflag
	integer mm1,mm2,nn1,nn2,M1,M2,N1,N2,ii,jj,kk,j,i,mm,nn
	integer level_butterfly, num_blocks, level_butterfly_c, num_blocks_c,level,num_col,num_row,num_rowson,num_colson
	
    type(matrixblock)::blocks_i,blocks_A,blocks_B,blocks_C,blocks_D
	complex(kind=8),allocatable:: matrixtemp1(:,:),matrixtemp2(:,:),vin(:,:),vout1(:,:),vout2(:,:)
	complex(kind=8)::ctemp1,ctemp2

	blocks_A%level = blocks_i%level+1
	blocks_A%row_group = blocks_i%row_group*2	
	blocks_A%col_group = blocks_i%col_group*2
	blocks_A%style = blocks_i%style

	blocks_B%level = blocks_i%level+1
	blocks_B%row_group = blocks_i%row_group*2	
	blocks_B%col_group = blocks_i%col_group*2+1
	blocks_B%style = blocks_i%style
	
	blocks_C%level = blocks_i%level+1
	blocks_C%row_group = blocks_i%row_group*2+1	
	blocks_C%col_group = blocks_i%col_group*2
	blocks_C%style = blocks_i%style
	
	blocks_D%level = blocks_i%level+1
	blocks_D%row_group = blocks_i%row_group*2+1	
	blocks_D%col_group = blocks_i%col_group*2+1
	blocks_D%style = blocks_i%style	
	
	
	if(blocks_i%level_butterfly==0)then
		blocks_A%level_butterfly = 0
		blocks_B%level_butterfly = 0
		blocks_C%level_butterfly = 0
		blocks_D%level_butterfly = 0
		
		allocate(blocks_A%ButterflyU(1))
		allocate(blocks_B%ButterflyU(1))
		allocate(blocks_C%ButterflyU(1))
		allocate(blocks_D%ButterflyU(1))
		
		allocate(blocks_A%ButterflyV(1))
		allocate(blocks_B%ButterflyV(1))
		allocate(blocks_C%ButterflyV(1))
		allocate(blocks_D%ButterflyV(1))	

		mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
		nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
		mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
		nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1
		kk = size(blocks_i%ButterflyU(1)%matrix,2)
		
		
		
		allocate(blocks_A%ButterflyU(1)%matrix(mm1,kk))
		allocate(blocks_A%ButterflyV(1)%matrix(nn1,kk))		
		blocks_A%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1:mm1,1:kk)
		blocks_A%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1:nn1,1:kk)
		blocks_A%rankmax = kk
		blocks_A%rankmin = kk		
		
		allocate(blocks_B%ButterflyU(1)%matrix(mm1,kk))
		allocate(blocks_B%ButterflyV(1)%matrix(nn2,kk))		
		blocks_B%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1:mm1,1:kk)
		blocks_B%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1+nn1:nn1+nn2,1:kk)
		blocks_B%rankmax = kk
		blocks_B%rankmin = kk	
		
		allocate(blocks_C%ButterflyU(1)%matrix(mm2,kk))
		allocate(blocks_C%ButterflyV(1)%matrix(nn1,kk))		
		blocks_C%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1+mm1:mm1+mm2,1:kk)
		blocks_C%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1:nn1,1:kk)
		blocks_C%rankmax = kk
		blocks_C%rankmin = kk
		
		allocate(blocks_D%ButterflyU(1)%matrix(mm2,kk))
		allocate(blocks_D%ButterflyV(1)%matrix(nn2,kk))		
		blocks_D%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1+mm1:mm1+mm2,1:kk)
		blocks_D%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1+nn1:nn1+nn2,1:kk)
		blocks_D%rankmax = kk
		blocks_D%rankmin = kk		
		
	
	else if(blocks_i%level_butterfly==1)then
		blocks_A%level_butterfly = 0
		blocks_B%level_butterfly = 0
		blocks_C%level_butterfly = 0
		blocks_D%level_butterfly = 0
		
		allocate(blocks_A%ButterflyU(1))
		allocate(blocks_B%ButterflyU(1))
		allocate(blocks_C%ButterflyU(1))
		allocate(blocks_D%ButterflyU(1))
		
		allocate(blocks_A%ButterflyV(1))
		allocate(blocks_B%ButterflyV(1))
		allocate(blocks_C%ButterflyV(1))
		allocate(blocks_D%ButterflyV(1))
		
		mm1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,1)
		nn1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,2)
		mm2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,1)
		nn2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,2)
		
		M1 = size(blocks_i%ButterflyU(1)%matrix,1)
		M2 = size(blocks_i%ButterflyU(2)%matrix,1)
		N1 = size(blocks_i%ButterflyV(1)%matrix,1)
		N2 = size(blocks_i%ButterflyV(2)%matrix,1)
		
		allocate(blocks_A%ButterflyU(1)%matrix(M1,nn1))
		allocate(blocks_A%ButterflyV(1)%matrix(N1,nn1))		
		call gemm_omp(blocks_i%ButterflyU(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,blocks_A%ButterflyU(1)%matrix, M1, mm1, nn1)
		blocks_A%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix
		blocks_A%rankmax = nn1
		blocks_A%rankmin = nn1
		
		allocate(blocks_B%ButterflyU(1)%matrix(M1,nn2))
		allocate(blocks_B%ButterflyV(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2)%matrix,blocks_B%ButterflyU(1)%matrix, M1, mm1, nn2)
		blocks_B%ButterflyV(1)%matrix = blocks_i%ButterflyV(2)%matrix		
		blocks_B%rankmax = nn2
		blocks_B%rankmin = nn2
		
		allocate(blocks_C%ButterflyU(1)%matrix(M2,nn1))
		allocate(blocks_C%ButterflyV(1)%matrix(N1,nn1))
		call gemm_omp(blocks_i%ButterflyU(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,1)%matrix,blocks_C%ButterflyU(1)%matrix, M2, mm2, nn1)
		blocks_C%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix
		blocks_C%rankmax = nn1
		blocks_C%rankmin = nn1				
		
		allocate(blocks_D%ButterflyU(1)%matrix(M2,nn2))
		allocate(blocks_D%ButterflyV(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,blocks_D%ButterflyU(1)%matrix, M2, mm2, nn2)
		blocks_D%ButterflyV(1)%matrix = blocks_i%ButterflyV(2)%matrix			
		blocks_D%rankmax = nn2
		blocks_D%rankmin = nn2
		
		
	else 
		blocks_A%level_butterfly = blocks_i%level_butterfly-2
		blocks_B%level_butterfly = blocks_i%level_butterfly-2
		blocks_C%level_butterfly = blocks_i%level_butterfly-2
		blocks_D%level_butterfly = blocks_i%level_butterfly-2
		
		level_butterfly_c = blocks_i%level_butterfly-2
		num_blocks_c = 2**level_butterfly_c

		allocate(blocks_A%ButterflyU(num_blocks_c))
		allocate(blocks_A%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,2)
			allocate(blocks_A%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_A%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
		
			mm1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,1)
			allocate(blocks_A%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_A%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
		end do
		
		allocate(blocks_B%ButterflyU(num_blocks_c))
		allocate(blocks_B%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,2)
			allocate(blocks_B%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_B%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_B%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))	
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_B%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)			
		end do

		allocate(blocks_C%ButterflyU(num_blocks_c))
		allocate(blocks_C%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,2)
			allocate(blocks_C%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_C%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,1)
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			allocate(blocks_C%ButterflyV(ii)%matrix(mm1+mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_C%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)				
		end do
		allocate(blocks_D%ButterflyU(num_blocks_c))
		allocate(blocks_D%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,2)
			allocate(blocks_D%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_D%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_D%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_D%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)		
		end do		
	
	
		if(level_butterfly_c/=0)then
			allocate(blocks_A%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_B%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_C%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_D%ButterflyKerl(level_butterfly_c))	
		end if  
		do level=1, level_butterfly_c
             num_col=blocks_i%ButterflyKerl(level+1)%num_col
             num_row=blocks_i%ButterflyKerl(level+1)%num_row
             num_colson=num_col/2
             num_rowson=num_row/2
			 blocks_A%ButterflyKerl(level)%num_row=num_rowson
			 blocks_A%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_A%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_B%ButterflyKerl(level)%num_row=num_rowson
			 blocks_B%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_B%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_C%ButterflyKerl(level)%num_row=num_rowson
			 blocks_C%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_C%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_D%ButterflyKerl(level)%num_row=num_rowson
			 blocks_D%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_D%ButterflyKerl(level)%blocks(num_rowson,num_colson))			 

			do j=1, num_col
				 do i=1, num_row
					 mm=size(blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix,1)
					 nn=size(blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix,2)
					 if (i<=num_rowson .and. j<=num_colson) then
						 allocate (blocks_A%ButterflyKerl(level)%blocks(i,j)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_A%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j<=num_colson) then
						 allocate (blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)                      
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i<=num_rowson .and. j>num_colson) then
						 allocate (blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j>num_colson) then
						 allocate (blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 endif
				 enddo
			 enddo
		enddo
		
		call get_butterfly_minmaxrank(blocks_A)
		call get_butterfly_minmaxrank(blocks_B)
		call get_butterfly_minmaxrank(blocks_C)
		call get_butterfly_minmaxrank(blocks_D)
		
	end if
	
	! if(verboselevel>=1)write(*,'(A33,I5,A7,I3,I3,I3,I3)')'  In split. L_butt: ',blocks_i%level_butterfly,' rank: ', blocks_A%rankmax,blocks_B%rankmax,blocks_C%rankmax,blocks_D%rankmax
	

	
	mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
	nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
	mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
	nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1

	
	! allocate(vin(nn1+nn2,1))
	! vin = 1
	! allocate(vout1(mm1+mm2,1))
	! vout1 = 0
	! allocate(vout2(mm1+mm2,1))
	! vout2 = 0
	
	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_i,'N',mm1+mm2,nn1+nn2,1,vin,vout1,ctemp1,ctemp2)
	
	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm1,nn1,1,vin(1:nn1,:),vout2(1:mm1,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm1,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1:mm1,:),ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_C,'N',mm2,nn1,1,vin(1:nn1,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call butterfly_block_MVP_randomized_dat(blocks_D,'N',mm2,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	
	! write(*,*)'spliting error:',fnorm(vout1-vout2,mm1+mm2,1)/fnorm(vout1,mm1+mm2,1)
	! deallocate(vin,vout1,vout2)
	
	
end subroutine Butterfly_split


! blocks_D: D^-1 - I
! blocks_B: B
! blocks_C: C
! blocks_A: (A-BD^-1C)^-1 - I
subroutine butterfly_block_MVP_inverse_ABCD_dat(blocks_A,blocks_B,blocks_C,blocks_D,trans,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock)::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn

    ctemp1=1.0d0
    ctemp2=0.0d0

	groupn=blocks_B%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=blocks_B%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	call assert(mm+nn==N,'mm+nn/=N')  
	
	
    allocate(Vin_tmp(N,num_vect_sub))
    Vin_tmp = Vin
 	Vout = 0  

    allocate(Vbuff(nn,num_vect_sub))
    Vbuff = 0
   
	if(trans=='N')then
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vbuff,ctemp1,ctemp2)
		Vbuff = Vbuff + Vin(1+mm:N,1:num_vect_sub)
		call butterfly_block_MVP_randomized_dat(blocks_B,trans,mm,nn,num_vect_sub,&
		&Vbuff,Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)- Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)

		if(isnan(fnorm(Vout,N,num_vect_sub)))then
			write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD11N'
			stop
		end if
		
		! write(2111,*)abs(Vout)
		
		call butterfly_block_MVP_randomized_dat(blocks_A,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vout(1+mm:mm+nn,1:num_vect_sub),Vin(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2)
		Vin = Vout + Vin

		if(isnan(fnorm(Vin,N,num_vect_sub)))then
			write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD22N'
			stop
		end if		
		! write(2112,*)abs(Vin)			
		
		call butterfly_block_MVP_randomized_dat(blocks_C,trans,nn,mm,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vbuff,ctemp1,ctemp2)
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vbuff, Vout(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2)
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
	
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vbuff,ctemp1,ctemp2)
		Vbuff = Vbuff + Vin(1+mm:N,1:num_vect_sub)
		call butterfly_block_MVP_randomized_dat(blocks_C,trans,nn,mm,num_vect_sub,&
		&Vbuff,Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub) - Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) 

		if(isnan(fnorm(Vout,N,num_vect_sub)))then
			write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD11T'
			stop
		end if		
		
		call butterfly_block_MVP_randomized_dat(blocks_A,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)				
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vout(1+mm:mm+nn,1:num_vect_sub),Vin(1+mm:mm+nn,1:num_vect_sub),ctemp1,ctemp2)
		Vin = Vout + Vin

		if(isnan(fnorm(Vin,N,num_vect_sub)))then
			write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD22T'
			stop
		end if		
		
		call butterfly_block_MVP_randomized_dat(blocks_B,trans,mm,nn,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vbuff,ctemp1,ctemp2)
		call butterfly_block_MVP_randomized_dat(blocks_D,trans,nn,nn,num_vect_sub,&
		&Vbuff,Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub)+Vbuff
		
		if(isnan(fnorm(Vout,N,num_vect_sub)))then
			write(*,*)fnorm(Vin,N,num_vect_sub),fnorm(Vout,N,num_vect_sub),'ABCD33T'
			stop
		end if	
		  

		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
		
	end if



   Vin = Vin_tmp
   deallocate(Vin_tmp)
   deallocate(Vbuff)
   
end subroutine butterfly_block_MVP_inverse_ABCD_dat



end module Butterfly_inversion_schur_partition
