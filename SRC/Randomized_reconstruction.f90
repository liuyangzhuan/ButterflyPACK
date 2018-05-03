module Randomized_reconstruction
use Utilites_randomized
use Butterfly_rightmultiply
contains 


subroutine Butterfly_randomized(level_butterfly,rank0,rankrate,blocks_o,operand,blackbox_MVP_dat,error_inout) 

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	use Butterfly_compress_forward
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, rank0, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,r1,r2,r3,r3tmp,mn,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock)::blocks_o
    type(matrixblock),pointer::blocks_A, blocks_B, blocks_C, blocks_D,block_tmp
    integer rank_new_max,rank_pre_max
	real*8:: rank_new_avr,error,rankrate
	integer niter,groupm,groupn
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:),U1(:,:),V1(:,:),U2(:,:),V2(:,:),U3(:,:),V3(:,:),U3tmp(:,:),V3tmp(:,:),UUtmp(:,:),VVtmp(:,:),UU(:,:),VV(:,:),UUr(:,:),VVr(:,:)
	real*8,allocatable :: Singular(:)
	complex (kind=8), allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:),Vout3(:,:),Vout4(:,:),Vout(:,:)
	complex (kind=8)::ctemp1,ctemp2
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matsub_tmp(:,:)
	integer idx_start_m_ref,option
	class(*):: operand
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface
	
	ctemp1 = 1d0; ctemp2 = 0d0
	Memory = 0
	
	groupm=blocks_o%row_group
	groupn=blocks_o%col_group
	
	do tt = 1,10
	
		rank_pre_max = rank0*rankrate**(tt-1)
	
		n1 = OMP_get_wtime()
		call Initialize_Butterfly_randomized(level_butterfly,rank_pre_max,groupm,groupn)
		n2 = OMP_get_wtime()
		Time_Init_inverse = Time_Init_inverse + n2-n1
		
		n1 = OMP_get_wtime()

		call Reconstruction_LL(blocks_o,operand,blackbox_MVP_dat)	
		call Reconstruction_RR(blocks_o,operand,blackbox_MVP_dat)
		call Test_Reconstruction_Error(blocks_o,operand,blackbox_MVP_dat,error_inout)
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
	write(*,*)'randomized scheme not converged in inverse_ABCD. level: ',blocks_o%level_butterfly,error_inout,rank_new_max
	stop

    return

end subroutine Butterfly_randomized





subroutine Reconstruction_LL(blocks_o,operand,blackbox_MVP_dat)
    
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
	integer level_right_start,num_col,num_row,option
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	real*8:: error_inout
	integer,allocatable::perms(:)
	type(partitionedblocks)::partitioned_block
	class(*):: operand	
	type(matrixblock)::blocks_o
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface	
	
	
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
			call Get_Randomized_Vectors_LL(blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth)
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
    
end subroutine Reconstruction_LL



subroutine Reconstruction_RR(blocks_o,operand,blackbox_MVP_dat)
    
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
    real*8::n1,n2
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
	! type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	type(partitionedblocks)::partitioned_block
	
	type(matrixblock)::blocks_o
	class(*):: operand
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface	
	
	
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
			call Get_Randomized_Vectors_RR(blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth)
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
	
    return
    
end subroutine Reconstruction_RR




subroutine Test_Reconstruction_Error(block_o,operand,blackbox_MVP_dat,error)

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

	
	
	type(matrixblock)::block_o
	class(*)::operand
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface
	
	
	level_butterfly=butterfly_block_randomized(1)%level_butterfly
	num_blocks=2**level_butterfly
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call Get_Randomized_Vectors_RR_Test(block_o,operand,blackbox_MVP_dat,num_vect)
	
	
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

end subroutine Test_Reconstruction_Error




subroutine Get_Randomized_Vectors_LL(blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk
    integer mm,nn,mn,level_butterfly,groupm,groupn,groupm_diag
	character chara
	integer*8 idx_start   
   
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_right_start
	type(RandomBlock), pointer :: random
	integer Nsub,Ng
	
	class(*):: operand	
	type(matrixblock)::blocks_o	
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface		
	
	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

	groupm=blocks_o%row_group ! Note: row_group and col_group interchanged here
	groupn=blocks_o%col_group	
	
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 


	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub)) ! #3 is used for output 
	
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
				
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	

	call blackbox_MVP_dat(operand,'T',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)		


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

end subroutine Get_Randomized_Vectors_LL




subroutine Get_Randomized_Vectors_RR(blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,level_butterfly,groupm,groupn
    character chara
    
	integer*8 idx_start   
    integer groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	type(RandomBlock), pointer :: random
	integer Nsub,Ng
	
	class(*):: operand	
	type(matrixblock)::blocks_o	
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface		
	


	
	
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

	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used as output 
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupn_start=groupn*2**(butterfly_block_randomized(1)%level_butterfly)
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i,header_m,tailer_m,nn,k)	
		do i=(nth-1)*Ng+1, nth*Ng
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_n=basis_group(groupn_start+i-1)%head
				tailer_n=basis_group(groupn_start+i-1)%tail
				nn=tailer_n-header_n+1
				k=header_n-header_nn	
				call RandomMat(nn,num_vect_subsub,min(nn,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:nn+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	call blackbox_MVP_dat(operand,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)		

	
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

end subroutine Get_Randomized_Vectors_RR


subroutine Get_Randomized_Vectors_RR_Test(block_o,operand,blackbox_MVP_dat,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,level_butterfly,groupm,groupn
    character chara
    complex(kind=8) ctemp, ctemp1, ctemp2

    real*8,allocatable :: Singular(:)

	integer*8 idx_start   
    integer groupn_start,dimension_rank
    integer header_nn
	integer header_n, tailer_n
	
	type(RandomBlock), pointer :: random
	
	type(matrixblock)::block_o
	class(*)::operand
	integer num_vect_sub
	
	interface
		subroutine blackbox_MVP_dat(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout)
			implicit none
			class(*)::operand	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
		end subroutine blackbox_MVP_dat
	end interface
	
    level_butterfly=butterfly_block_randomized(1)%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    groupm=block_o%row_group
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupn_start=groupm*2**(level_butterfly)
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	
	do i=1, num_blocks			
		header_n=basis_group(groupn_start+i-1)%head
		tailer_n=basis_group(groupn_start+i-1)%tail
		nn=tailer_n-header_n+1
		k=header_n-header_nn	
		call RandomMat(nn,num_vect_sub,min(nn,num_vect_sub),RandomVectors_InOutput(1)%vector(1+k:nn+k,1:num_vect_sub),0)
	end do

	! get the right multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
	call blackbox_MVP_dat(operand,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)

	
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

end subroutine Get_Randomized_Vectors_RR_Test



subroutine get_minmaxrank_ABCD(partitioned_block,rankmax)

    use MODULE_FILE
    implicit none
    
	integer rankmax
    type(partitionedblocks)::partitioned_block
	
	rankmax = -1000
	rankmax = max(rankmax,partitioned_block%blocks_A%rankmax)
	rankmax = max(rankmax,partitioned_block%blocks_B%rankmax)
	rankmax = max(rankmax,partitioned_block%blocks_C%rankmax)
	rankmax = max(rankmax,partitioned_block%blocks_D%rankmax)
	
end subroutine get_minmaxrank_ABCD



! blocks_D: D^-1 - I
! blocks_B: B
! blocks_C: C
! blocks_A: (A-BD^-1C)^-1 - I
subroutine butterfly_block_MVP_inverse_ABCD_dat(partitioned_block,block_o,trans,M,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
   type(matrixblock)::block_o

   select TYPE(partitioned_block)
   
   type is (partitionedblocks)
		call assert(M==N,'M/=N in butterfly_block_MVP_inverse_ABCD_dat')
	   
		blocks_A => partitioned_block%blocks_A
		blocks_B => partitioned_block%blocks_B
		blocks_C => partitioned_block%blocks_C
		blocks_D => partitioned_block%blocks_D	
	   
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
	   Vout = Vout-Vin 
	   
	   deallocate(Vin_tmp)
	   deallocate(Vbuff)  

   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_inverse_ABCD_dat





! blocks_D: D^-1 - I
! blocks_B: B 
! blocks_C: C
! blocks_A: (A-BD^-1C)^-1 - I
subroutine butterfly_block_MVP_inverse_A_minusBDinvC_dat(partitioned_block,block_o,trans,M,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: V_tmp1(:,:),V_tmp2(:,:),Vin_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
	type(matrixblock)::block_o

   select TYPE(partitioned_block)
   
   type is (partitionedblocks)
		call assert(M==N,'M/=N in butterfly_block_MVP_inverse_A_minusBDinvC_dat')
	   
		blocks_A => partitioned_block%blocks_A
		blocks_B => partitioned_block%blocks_B
		blocks_C => partitioned_block%blocks_C
		blocks_D => partitioned_block%blocks_D	


		groupn=blocks_B%col_group    ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
		groupm=blocks_B%row_group    ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	

		Vout=0
		allocate(Vin_tmp(M,num_vect_sub))
		Vin_tmp=Vin
		
		allocate(V_tmp1(nn,num_vect_sub))
		V_tmp1 = 0
		allocate(V_tmp2(nn,num_vect_sub))
		V_tmp2 = 0
	   
		if(trans=='N')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_C,'N',nn,mm,num_vect_sub,Vin_tmp,V_tmp1,ctemp1,ctemp2)	
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_D,'N',nn,nn,num_vect_sub,V_tmp1,V_tmp2,ctemp1,ctemp2)		
			V_tmp2 = V_tmp2 + V_tmp1
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm,nn,num_vect_sub,V_tmp2,Vout,ctemp1,ctemp2)	
			ctemp1=1.0d0 ; ctemp2=-1.0d0
			call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm,mm,num_vect_sub,Vin_tmp,Vout,ctemp1,ctemp2)
	
			
		else if(trans=='T')then
		
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_B,'T',mm,nn,num_vect_sub,Vin_tmp,V_tmp1,ctemp1,ctemp2)
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_D,'T',nn,nn,num_vect_sub,V_tmp1,V_tmp2,ctemp1,ctemp2)
			V_tmp2 = V_tmp2 + V_tmp1
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(blocks_C,'T',nn,mm,num_vect_sub,V_tmp2,Vout,ctemp1,ctemp2)	
			ctemp1=1.0d0 ; ctemp2=-1.0d0
			call butterfly_block_MVP_randomized_dat(blocks_A,'T',mm,mm,num_vect_sub,Vin_tmp,Vout,ctemp1,ctemp2)		
		end if



	   Vin = Vin_tmp

	   
	   deallocate(Vin_tmp)
	   deallocate(V_tmp1)  
	   deallocate(V_tmp2)  

   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_inverse_A_minusBDinvC_dat



subroutine butterfly_block_MVP_inverse_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock),pointer::block_off1,block_off2
   integer groupn,groupm,mm,nn
   class(*)::ho_bf1
   type(matrixblock)::block_o

   select TYPE(ho_bf1)
   
   type is (hobf)

		block_off1 => ho_bf%levels(ho_bf1%ind_lv)%BP(ho_bf1%ind_bk*2-1)%LL(1)%matrices_block(1)	
		block_off2 => ho_bf%levels(ho_bf1%ind_lv)%BP(ho_bf1%ind_bk*2)%LL(1)%matrices_block(1)			
		
		groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
		groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
		
		allocate(Vin_tmp(N,num_vect_sub))
		Vin_tmp = Vin
		Vout = 0  

		allocate(Vbuff(nn,num_vect_sub))
		Vbuff = 0
	   
		if(trans=='N')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(block_off2,'N',nn,mm,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2)	
			call butterfly_block_MVP_randomized_dat(block_off1,'N',mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2)	
			Vout = -Vout		
			
		else if(trans=='T')then
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call butterfly_block_MVP_randomized_dat(block_off1,'T',mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2)	
			call butterfly_block_MVP_randomized_dat(block_off2,'T',nn,mm,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2)	
			Vout = -Vout			
		end if

	   Vin = Vin_tmp

	   deallocate(Vin_tmp)
	   deallocate(Vbuff)  

   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_inverse_minusBC_dat



subroutine butterfly_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,M,N
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character trans
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vbuff(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,level_right_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)
	
	class(*)::ho_bf1
	type(matrixblock)::block_o
	

   select TYPE(ho_bf1)
   
   type is (hobf)
		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk
	
		if(trans=='N')then

			level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
			num_blocks=2**level_butterfly

			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	

			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 

			allocate(Vbuff(mm,num_vect_sub))
			Vbuff=0
			
			groupn_start=groupn*2**(level_butterfly)
			header_nn=basis_group(groupn_start)%head
			idx_start = 1

			! get the right multiplied vectors
			idx_start_glo = basis_group(groupm)%head
			ctemp1=1.0d0 ; ctemp2=0.0d0
		n1 = OMP_get_wtime()  
		  call butterfly_block_MVP_randomized_dat(block_o,'N',Vin,Vbuff,ctemp1,ctemp2)
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1	
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
			allocate(vec_new(mm,num_vect_sub))

			do level = Maxlevel_for_blocks+1,level_c+1,-1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==Maxlevel_for_blocks+1)then
					n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						! write(*,*)level,ii
						groupm_diag = ho_bf%levels(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1

						call fullmat_block_MVP_randomized_dat(ho_bf%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
		
					end do
					! !$omp end parallel do
					
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						! write(*,*)level,ii
						groupm_diag = ho_bf%levels(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						
						call OneL_block_MVP_inverse_dat(ho_bf,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))

					end do		
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1			
				end if
				
				Vbuff = vec_new
			end do

			Vout = vec_new
			
			deallocate(Vbuff)
			deallocate(vec_new)
	
		else 
			ctemp1=1.0d0 ; ctemp2=0.0d0	

			level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
			num_blocks=2**level_butterfly

			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	

			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			allocate(Vbuff(mm,num_vect_sub))
			Vbuff=0 
			 
			groupm_start=groupm*2**(level_butterfly)
			header_mm=basis_group(groupm_start)%head
			idx_start = 1
			
			! get the left multiplied vectors
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
			idx_start_glo = basis_group(groupm)%head		

			allocate(vec_new(mm,num_vect_sub))	
			Vbuff = Vin
			do level = level_c+1,Maxlevel_for_blocks+1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==Maxlevel_for_blocks+1)then
					n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf%levels(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						call fullmat_block_MVP_randomized_dat(ho_bf%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)

					end do
					! !$omp end parallel do
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf%levels(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
						call OneL_block_MVP_inverse_dat(ho_bf,level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))					
					end do
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1	
				end if
				Vbuff = vec_new
			end do	
			
			Vbuff = vec_new
			
			deallocate(vec_new)	

			n1 = OMP_get_wtime()
			call butterfly_block_MVP_randomized_dat(block_o,'T',Vbuff,Vout,ctemp1,ctemp2)	
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1		

			deallocate(Vbuff)
	
		end if	
	
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
    return                

end subroutine butterfly_block_MVP_Sblock_dat




end module Randomized_reconstruction