module Bplus_inversion_schur_partition
use Butterfly_inversion_schur_partition

integer rankthusfarBC
contains 



subroutine Bplus_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
	type(blockplus),pointer::bplus
	real*8:: n1,n2,Memory

    bplus => cascading_factors(level_c)%BP_inverse_schur(rowblock)	
	
	if(bplus%Lplus==1)then
		call OneL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)
	else 
		call MultiL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)
	end if	
	
	

end subroutine Bplus_inverse_schur_partitionedinverse




subroutine OneL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

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

	block_off => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)			
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)	
	block_o%level_butterfly = block_off%level_butterfly	
		
	Memory = 0
	
	allocate(partitioned_blocks(0:block_o%level_butterfly+1))
	do ll=0,block_o%level_butterfly+1
		allocate(partitioned_blocks(ll)%blocks_A)
		allocate(partitioned_blocks(ll)%blocks_B)
		allocate(partitioned_blocks(ll)%blocks_C)
		allocate(partitioned_blocks(ll)%blocks_D)
	end do
	blocks_minusBC => partitioned_blocks(0)%blocks_D
	
	
	error_avr_glo = 0
	error_cnt = 0	
	call OneL_diagonal_minusBC(level_c,rowblock,blocks_minusBC)
	call Butterfly_inverse_partitionedinverse_IplusButter(0,1)	

	if(error_cnt/=0)write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',blocks_minusBC%rankmax,' L_butt:',blocks_minusBC%level_butterfly,' error_avr:',error_avr_glo/error_cnt	

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

end subroutine OneL_inverse_schur_partitionedinverse




subroutine OneL_diagonal_minusBC(level_c,rowblock,blocks_minusBC)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	use Butterfly_compress_forward
	
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
	integer itermax,ntry,idx_start_m_ref,idx_start_n_ref
	real*8:: n1,n2,Memory,tmpfact
	complex (kind=8), allocatable::matrix_small(:,:)
	complex (kind=8)::ctemp1,ctemp2
	complex(kind=8),allocatable:: matin(:,:),matinter(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:),matsub_tmp(:,:)
	type(matrixblock)::block_tmp
	
	Memory = 0
	
    block_o =>  cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1) 
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
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
	
	



	! if(level_c<=2)then
		! tmpfact = SVD_tolerance_factor
		! ! SVD_tolerance_factor = SVD_tolerance_factor* 1d-2
		
		! write(*,*)'explicit computing butterfly for BC '

		! mm = basis_group(blocks_minusBC%row_group)%tail - basis_group(blocks_minusBC%row_group)%head + 1
		! nn = basis_group(blocks_minusBC%col_group)%tail - basis_group(blocks_minusBC%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(nn,nn))
		! ! allocate(matout(mm,mm))
		! allocate(matsub_glo(nn,nn))
		! allocate(matinter(mm,nn))
		! matin = 0
		! do ii=1,nn
			! matin(ii,ii) = 1
		! end do
		
		! ctemp1=1.0d0 ; ctemp2=0.0d0
		! call butterfly_block_MVP_randomized_dat(block_off1,'N',mm,nn,nn,matin,matinter,ctemp1,ctemp2)	
		! call butterfly_block_MVP_randomized_dat(block_off2,'N',nn,mm,nn,matinter,matsub_glo,ctemp1,ctemp2)	
		! matsub_glo = -matsub_glo
					
		! block_tmp%row_group = blocks_minusBC%row_group
		! block_tmp%col_group = blocks_minusBC%col_group
		! block_tmp%level = blocks_minusBC%level
		! block_tmp%style = blocks_minusBC%style
		
		! idx_start_m_ref = basis_group(blocks_minusBC%row_group)%head
		! idx_start_n_ref = basis_group(blocks_minusBC%col_group)%head
		! ! write(*,*)'ddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		! write(*,*)'wocaonimaBC',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! deallocate(matin)
		! deallocate(matsub_glo)
		! deallocate(matinter)
		! call delete_blocks(block_tmp)	



		! mm = basis_group(blocks_minusBC%row_group)%tail - basis_group(blocks_minusBC%row_group)%head + 1
		! nn = basis_group(blocks_minusBC%col_group)%tail - basis_group(blocks_minusBC%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(nn,nn))
		! allocate(matinter(mm,nn))
		! allocate(matsub_glo(nn,nn))
		! allocate(matsub_tmp(nn,nn))
		! matin = 0
		! do ii=1,nn
			! matin(ii,ii) = 1
		! end do
		
		
		! ctemp1=1.0d0 ; ctemp2=0.0d0
		! call butterfly_block_MVP_randomized_dat(block_off2,'T',nn,mm,nn,matin,matinter,ctemp1,ctemp2)	
		! call butterfly_block_MVP_randomized_dat(block_off1,'T',mm,nn,nn,matinter,matsub_tmp,ctemp1,ctemp2)				
		! call copymatT_OMP(matsub_tmp,matsub_glo,nn,nn)
		! matsub_glo = -matsub_glo		
		! deallocate(matsub_tmp)
			
		! block_tmp%row_group = blocks_minusBC%row_group
		! block_tmp%col_group = blocks_minusBC%col_group
		! block_tmp%level = blocks_minusBC%level
		! block_tmp%style = blocks_minusBC%style
		
		! idx_start_m_ref = basis_group(blocks_minusBC%row_group)%head
		! idx_start_n_ref = basis_group(blocks_minusBC%col_group)%head
! ! write(*,*)'dddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		
		! write(*,*)'wocaonimaBC',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! ! ranktmp_glo = block_tmp%rankmax*1.2
		
		! deallocate(matin)
		! deallocate(matsub_glo)
		! deallocate(matinter)
		! call delete_blocks(block_tmp)	
		! write(*,*)'Done computing ' 
		! SVD_tolerance_factor = tmpfact
	! end if	
	
		do tt = 1,10
			do ntry=1,1
			itermax = 0
			n1 = OMP_get_wtime()
			call Initialize_Butterfly_inverse_BC(level_c,rowblock,tt-1)
			n2 = OMP_get_wtime()
			Time_Init_inverse = Time_Init_inverse + n2-n1
			
			n1 = OMP_get_wtime()
			call OneL_Reconstruction_LL_BC(level_c,rowblock)	
			call OneL_Reconstruction_RR_BC(level_c,rowblock,error_inout)
			n2 = OMP_get_wtime()	
			Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2 - n1
			
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			! write(*,*)tt,error_inout,butterfly_block_randomized(1)%dimension_rank,butterfly_block_randomized(1)%rankmax,'OneL_mBC'
			
			
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
				if(verboselevel>=1)write(*,'(A20,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' mBC ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',blocks_minusBC%level_butterfly,' error:',error_inout	
				error_avr_glo = error_avr_glo + error_inout
				error_cnt = error_cnt+1
				
				return			
			end if		
			end do
		end do
		write(*,*)'randomized scheme not converged in BC'
		stop
	end if	
    return

end subroutine OneL_diagonal_minusBC



subroutine OneL_Reconstruction_LL_BC(level_c,rowblock)
    
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

	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
	
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
			call OneL_Get_Randomized_Vectors_LL_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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
    
end subroutine OneL_Reconstruction_LL_BC



subroutine OneL_Reconstruction_RR_BC(level_c,rowblock,error)
    
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
		
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
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
			call OneL_Get_Randomized_Vectors_RR_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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

	call OneL_Test_Error_RR_BC(level_c,rowblock,error)


	
    return
    
end subroutine OneL_Reconstruction_RR_BC




subroutine OneL_Test_Error_RR_BC(level_c,rowblock,error)

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
	
	
	block_o =>  cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1) 
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

	call OneL_Get_Randomized_Vectors_RR_Test_BC(level_c,rowblock,num_vect)
	
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

end subroutine OneL_Test_Error_RR_BC




subroutine OneL_Get_Randomized_Vectors_LL_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
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

end subroutine OneL_Get_Randomized_Vectors_LL_BC




subroutine OneL_Get_Randomized_Vectors_RR_BC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	
	
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
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
		! do i=1, num_blocks
		!$omp parallel do default(shared) private(i,header_m,tailer_m,mm,k)	
		do i=(nth-1)*Ng+1, nth*Ng		
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

end subroutine OneL_Get_Randomized_Vectors_RR_BC



subroutine OneL_Get_Randomized_Vectors_RR_Test_BC(level_c,rowblock,num_vect_sub)

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
	
	
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	block_off1 => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	
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

end subroutine OneL_Get_Randomized_Vectors_RR_Test_BC

subroutine MultiL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
    use Butterfly_compress_forward
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks,level_butterfly_loc
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll,llplus,bb,mmb
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off
    type(matrixblock),pointer::blocks_o_D
    type(matrixblock)::block_tmp
	type(blockplus),pointer::Bplus,Bplus_schur
    integer rank_new_max
	real*8:: rank_new_avr,error,error_avr
	integer niter
	real*8:: error_inout
	integer itermax,ntry,cnt,cnt_partial
	real*8:: n1,n2,Memory
	integer Lplus,level_BP,levelm,groupm_start,ij_loc,edge_s,edge_e,edge_first,idx_end_m_ref,idx_start_m_ref,idx_start_b,idx_end_b
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:)
	complex(kind=8):: ctemp1,ctemp2 
	integer, allocatable :: ipiv(:)
	
	ctemp1 = 1d0
	ctemp2 = 0d0
	
	block_off => cascading_factors(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)			
	block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
! write(*,*)block_o%row_group,block_o%col_group,level_c,rowblock,block_o%level,'diao'	
	block_o%level_butterfly = block_off%level_butterfly	
		
	Memory = 0
	

	call MultiL_diagonal_minusBC(level_c,rowblock) 
	
	
	! ! Bplus_schur => cascading_factors(level_c)%BP_inverse_schur(rowblock)
	! ! idx_start_m_ref = basis_group(Bplus_schur%row_group)%head
	! ! idx_end_m_ref = basis_group(Bplus_schur%row_group)%tail
	! ! mm = idx_end_m_ref - idx_start_m_ref + 1
		! ! ! write(*,*)'e',mm
	! ! allocate(matin(mm,mm))
	! ! allocate(matsub_glo(mm,mm))
		! ! ! write(*,*)'f'
	! ! matin = 0
	! ! do ii=1,mm
		! ! matin(ii,ii) = 1
	! ! end do
	! ! ! write(*,*)'g'
	! ! ! call Bplus_block_MVP_randomized_dat_partial(Bplus_schur,'N',mm,mm,mm,matin,matsub_glo,ctemp1,ctemp2,1,Bplus_schur%Lplus)
	! ! call Bplus_block_MVP_randomized_dat_partial(Bplus_schur,'N',mm,mm,mm,matin,matsub_glo,ctemp1,ctemp2,1,1)
	! ! do ii=1,mm
		! ! matsub_glo(ii,ii) = 1 + matsub_glo(ii,ii)
	! ! end do	
		! ! ! write(*,*)'h'
	
	! ! allocate(ipiv(mm))
	! ! call getrf(matsub_glo,ipiv)
	! ! ! write(*,*)shape(matrixtemp1)
	! ! call getri(matsub_glo,ipiv)	
	! ! deallocate(ipiv)	
	
	! ! do ii=1,mm
		! ! matsub_glo(ii,ii) =  matsub_glo(ii,ii) - 1
	! ! end do		
	
	! ! block_o => cascading_factors(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	! ! block_tmp%row_group = block_o%row_group
	! ! block_tmp%col_group = block_o%col_group
	! ! block_tmp%level = block_o%level
	! ! block_tmp%style = block_o%style
	
	! ! ! write(*,*)'wocaonidddma',block_tmp%level,block_o%level
	
	! ! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_m_ref)
	
	! ! write(*,*)'wocaonima',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,block_o%level
	! stop
	
	error_avr_glo = 0
	error_cnt = 0
	
	! write(*,*)'good!!!!'
	! stop
	Bplus => cascading_factors(level_c)%BP_inverse_schur(rowblock)
	Lplus = cascading_factors(level_c)%BP_inverse_schur(rowblock)%Lplus
	do llplus =Lplus,1,-1
		do bb=1,Bplus%LL(llplus)%Nbound
			block_o => Bplus%LL(llplus)%matrices_block(bb)
		
			!!!!! partial update butterflies at level llplus from left B1 = D^-1xB
			if(llplus/=Lplus)then
			
				n1 = OMP_get_wtime()
				call Butterfly_MoveSingulartoLeft(block_o)
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
				
				level_butterfly = int((Maxlevel_for_blocks - Bplus%LL(llplus)%matrices_block(1)%level)/2)*2 
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				level_butterfly_loc = levelm
				groupm_start=block_o%row_group*2**levelm
				edge_s =basis_group(block_o%row_group)%head 
				edge_e =basis_group(block_o%row_group)%tail 
				error_avr = 0
				cnt = 0
				cnt_partial = 0
				
				do ii=1,Bplus%LL(llplus+1)%Nbound
					edge_first = basis_group(Bplus%LL(llplus+1)%matrices_block(ii)%row_group)%head
					if(edge_first>=edge_s .and. edge_first<=edge_e)then
						ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1					
						if(level_butterfly_loc==0)then
							write(*,*)'level_butterfly_loc==0 not done'
							stop
						else 

							cnt_partial = cnt_partial + 1
							! allocate(agent_block(1))

							call Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,'L')
							call Extract_partial_Bplus(Bplus,llplus+1,Bplus%LL(llplus+1)%matrices_block(ii)%row_group)							
															

							do tt=1,10
								n1 = OMP_get_wtime()							
								call Initialize_Butterfly_SblockSmall(agent_block(1),tt-1)
								n2 = OMP_get_wtime()
								Time_Init_inverse = Time_Init_inverse + n2 -n1 				

					
								n1 = OMP_get_wtime()
								call MultiL_Reconstruction_LL_SchurSmall(agent_bplus(1),agent_block(1),'L')	
								call MultiL_Reconstruction_RR_SchurSmall(agent_bplus(1),agent_block(1),error_inout,'L')	 
								n2 = OMP_get_wtime()
								Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2-n1

								call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
								rank_new_max = butterfly_block_randomized(1)%rankmax
								
								if(error_inout>iter_tolerance)then
									call Delete_randomized_butterfly()
								else 
									call copy_randomizedbutterfly_partial(butterfly_block_randomized(1),block_o,level_butterfly_loc,ij_loc,'L',Memory)
									
									if(verboselevel>=2)write(*,'(A37,I7,A6,I3,A8,I2,A8,I3,A7,Es14.7)')'L small cnt: ',cnt_partial,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly_loc,' error:',error_inout
									error_avr = error_avr + error_inout
									cnt = cnt + 1
									error_avr_glo = error_avr_glo + error_inout
									error_cnt = error_cnt + 1
									
									call Delete_randomized_butterfly()
									call delete_blocks(agent_block(1))
									deallocate(agent_block)	
									call delete_bplus(agent_bplus(1))
									deallocate(agent_bplus)		
									exit
								end if					
							end do
							
							if(error_inout>iter_tolerance)then
								write(*,*)'randomized scheme not converged in left partial updates',error_inout,rank_new_max
								stop
							end if
									
						end if
					end if
				end do
				
				
				
				if(cnt>0)error_avr = error_avr/cnt
				if(verboselevel>=1)write(*,'(A30,I7,A6,I3,A11,Es14.7)')' L partial: ll ',llplus,' bb:',bb,' error_avr:',error_avr
			end if
			
			
   
			! write(*,*)block_o%level_butterfly,'ahaha'
			
			!!!!! invert I+B1 to be I+B2		
			allocate(partitioned_blocks(0:block_o%level_butterfly+1))
			do ll=0,block_o%level_butterfly+1
				allocate(partitioned_blocks(ll)%blocks_A)
				allocate(partitioned_blocks(ll)%blocks_B)
				allocate(partitioned_blocks(ll)%blocks_C)
				allocate(partitioned_blocks(ll)%blocks_D)
			end do
			blocks_o_D => partitioned_blocks(0)%blocks_D
			call copy_delete_butterfly(block_o,blocks_o_D)
			! call delete_blocks(block_o)
			
! write(*,*)blocks_o_D%level_butterfly,blocks_o_D%rankmax,'heihei'
! stop			
			call Butterfly_inverse_partitionedinverse_IplusButter(0,1)		
							   
			call copy_delete_butterfly(blocks_o_D,block_o,Memory)
			! call delete_blocks(blocks_o_D)	
			do ll=0,block_o%level_butterfly+1
				deallocate(partitioned_blocks(ll)%blocks_A)
				deallocate(partitioned_blocks(ll)%blocks_B)
				deallocate(partitioned_blocks(ll)%blocks_C)
				deallocate(partitioned_blocks(ll)%blocks_D)
			end do
			deallocate(partitioned_blocks)			
			
			
		
			!!!!! partial update butterflies at level llplus from right B2xD^-1
			if(llplus/=Lplus)then
				! write(*,*)'hhe'
				n1 = OMP_get_wtime()	
				call Butterfly_MoveSingulartoRight(block_o)
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1		
				! write(*,*)'hhe1'
				level_butterfly = int((Maxlevel_for_blocks - Bplus%LL(llplus)%matrices_block(1)%level)/2)*2 
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				level_butterfly_loc = levelm
				groupm_start=block_o%row_group*2**levelm
				edge_s =basis_group(block_o%row_group)%head 
				edge_e =basis_group(block_o%row_group)%tail 
				error_avr = 0
				cnt = 0
				cnt_partial = 0
				
				do ii=1,Bplus%LL(llplus+1)%Nbound
					edge_first = basis_group(Bplus%LL(llplus+1)%matrices_block(ii)%row_group)%head
					if(edge_first>=edge_s .and. edge_first<=edge_e)then
						ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1					
						if(level_butterfly_loc==0)then
							write(*,*)'level_butterfly_loc==0 not done'
							stop
						else 
! write(*,*)'hhe2',level_butterfly_loc,ij_loc,block_o%level_butterfly
							cnt_partial = cnt_partial + 1
							! allocate(agent_block(1))
						   
							call Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,'R')
							call Extract_partial_Bplus(Bplus,llplus+1,Bplus%LL(llplus+1)%matrices_block(ii)%row_group)							
							
							do tt=1,10
							 
								n1 = OMP_get_wtime()
								call Initialize_Butterfly_SblockSmall(agent_block(1),tt-1)
								n2 = OMP_get_wtime()
								Time_Init_inverse = Time_Init_inverse + n2 -n1 				

			! write(*,*)'hhe3'		
								n1 = OMP_get_wtime()
								call MultiL_Reconstruction_LL_SchurSmall(agent_bplus(1),agent_block(1),'R')	
			! write(*,*)'hhe3.1'
								call MultiL_Reconstruction_RR_SchurSmall(agent_bplus(1),agent_block(1),error_inout,'R')	 
								n2 = OMP_get_wtime()
								Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2-n1
			! write(*,*)'hhe4'
								call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
								rank_new_max = butterfly_block_randomized(1)%rankmax
			! write(*,*)'hhe5'
			
								if(error_inout>iter_tolerance)then
									call Delete_randomized_butterfly()
								else 			
									call copy_randomizedbutterfly_partial(butterfly_block_randomized(1),block_o,level_butterfly_loc,ij_loc,'R',Memory)
				! write(*,*)'hhe6'					
									if(verboselevel>=2)write(*,'(A37,I7,A6,I3,A8,I2,A8,I3,A7,Es14.7)')'R small cnt: ',cnt_partial,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly_loc,' error:',error_inout
									error_avr = error_avr + error_inout
									cnt = cnt + 1
									error_avr_glo = error_avr_glo + error_inout
									error_cnt = error_cnt + 1
									call Delete_randomized_butterfly()
									call delete_blocks(agent_block(1))
									deallocate(agent_block)	
									call delete_bplus(agent_bplus(1))
									deallocate(agent_bplus)	
									exit									
								end if
							end do
							if(error_inout>iter_tolerance)then
								write(*,*)'randomized scheme not converged in right partial updates',error_inout,rank_new_max
								stop
							end if							
						end if
					end if
				end do
				
				! call Butterfly_sym2asym(block_o)
				
				
				if(cnt>0)error_avr = error_avr/cnt
				if(verboselevel>=1)write(*,'(A30,I7,A6,I3,A11,Es14.7)')' R partial: ll ',llplus,' bb:',bb,' error_avr:',error_avr
			end if
		end do
	end do
	
	call ComputeMemory_Bplus(Bplus,Memory)	
	
	rank_new_max = 0
	do ll=1,Lplus
		rank_new_max = max(rank_new_max,Bplus%LL(ll)%rankmax)
	end do
	
	write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',rank_new_max,' L_butt:',Bplus%LL(1)%matrices_block(1)%level_butterfly,' error_avr:',error_avr_glo/error_cnt		
	
	
    return

end subroutine MultiL_inverse_schur_partitionedinverse






subroutine MultiL_diagonal_minusBC(level_c,rowblock)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    type(blockplus),pointer::bplus
	integer:: ii,ll,bb
    real*8 rtemp,error,Memory,n2,n1	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall,M,N
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	complex(kind=8) ctemp, ctemp1, ctemp2	
	integer level_c,rowblock,rank_new_max

	bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock)	
	! Memory = 0	
	call assert(bplus%Lplus>=2,'this is not a multi Bplus in MultiL_diagonal_minusBC')
	
	call Initialize_Bplus_FromInput(bplus)
	
	! write(*,*)'hhhhh1',Bplus_randomized(1)%LL(1)%matrices_block(1)%level
	n1 = OMP_get_wtime()
	error_cnt = 0
	error_avr_glo = 0
	rankthusfarBC = 0	
	do bb =1,Bplus_randomized(1)%LL(2)%Nbound
		! write(*,*)bb,Bplus_randomized(1)%LL(2)%Nbound,'dddd'
		call MultiLrandomized_Onesubblock_minusBC(level_c,rowblock,bb)
		! write(*,*)'go'
	end do
	! call Test_Error_RR_Inner_Exact(Bplus)
	n2 = OMP_get_wtime()
	Time_Oneblock_inverse = Time_Oneblock_inverse + n2-n1	
	
	call MultiLrandomized_Outter_minusBC_memfree(level_c,rowblock,rtemp)
		
	! write(*,*)'hhhhh2',Bplus_randomized(1)%LL(1)%matrices_block(1)%level,bplus%LL(1)%matrices_block(1)%level	
		
	call delete_Bplus(bplus)
	call copy_delete_Bplus(Bplus_randomized(1),bplus,Memory)
	! call copy_Bplus(Bplus_randomized(1),bplus,Memory)
	! call delete_Bplus(Bplus_randomized(1))
	deallocate(Bplus_randomized)
	
	! write(*,*)'hhhhh3',bplus%LL(1)%matrices_block(1)%level	
	
	rank_new_max = 0
	do ll=1,bplus%Lplus
		rank_new_max = max(rank_new_max,bplus%LL(ll)%rankmax)
	end do	
	
	if(verboselevel>=1)write(*,'(A20,A8,I3,A8,I3,A11,Es14.7)')' mBC ',' rank:',rank_new_max,' L_butt:',bplus%LL(1)%matrices_block(1)%level_butterfly,' error_avr:',error_avr_glo/error_cnt
	
    return

end subroutine MultiL_diagonal_minusBC





subroutine MultiLrandomized_Onesubblock_minusBC(level_c,rowblock,bb_o)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   use Butterfly_compress_forward
   implicit none

    type(blockplus),pointer::bplus_off1
	integer:: ii,ll,bb,jj,bb_o,tt
    real*8 Memory,rtemp,error,n1,n2, mem_vec	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vout3(:,:),Vin(:,:)
	integer M,N,idx_start_n,idx_start_m,idx_start_n_loc,idx_end_n_loc,idx_start_m_loc,idx_end_m_loc,mm,nn,rmax,rank,idx_start_n_ref,idx_start_m_ref,idx_end_n_ref,idx_end_m_ref
	complex(kind=8)::ctemp1,ctemp2
	type(matrixblock),pointer::blocks
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:)
	real*8, allocatable :: Singular(:)
	integer level_c,rowblock,Nactive
	integer,allocatable::boxindex(:)
	integer Chunksize, Nchunk, Nidx, idx_s,cc
		
	bplus_off1 =>  cascading_factors(level_c)%BP(2*rowblock-1)
	
	N = basis_group(Bplus_randomized(1)%col_group)%tail - basis_group(Bplus_randomized(1)%col_group)%head + 1	
	M = basis_group(Bplus_randomized(1)%row_group)%tail - basis_group(Bplus_randomized(1)%row_group)%head + 1	
	
	
	ctemp1 = 1.0d0
	ctemp2 = 0.0d0
	
	idx_start_n = basis_group(Bplus_randomized(1)%col_group)%head
	idx_start_m = basis_group(Bplus_randomized(1)%row_group)%head
		
	blocks => Bplus_randomized(1)%LL(2)%matrices_block(bb_o)
	idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
	idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
	idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
	idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
	
	mm = idx_end_m_loc - idx_start_m_loc + 1
	nn = idx_end_n_loc - idx_start_n_loc + 1
	
	
	do tt=1,10
	
		! rmax = bplus_off1%LL(2)%rankmax*2**(tt) !+ (tt-1)*10  !!!!! be careful here
		rmax = max(bplus_off1%LL(2)%rankmax*2,rankthusfarBC)*2**(tt-1) !+ (tt-1)*10  !!!!! be careful here
		
		
		
		! rmax = bplus_off1%LL(2)%rankmax*2 + (tt-1)*10  !!!!! be careful here
		if(rmax>min(mm,nn))rmax = min(mm,nn)
													   
		
		
		n1 = OMP_get_wtime()
		
		
		
		
		Chunksize = 100 !rmax
		Nchunk = ceiling_safe(dble(rmax)/dble(Chunksize))		
 
		allocate(matRrow(mm,rmax))
		matRrow=0	
		allocate(matZcRrow(nn,rmax))	
		matZcRrow=0
		allocate(matRcol(nn,rmax))
		matRcol=0
		allocate(matZRcol(mm,rmax))		
		matZRcol=0
		do cc=1,Nchunk
			idx_s = (cc-1)*Chunksize+1
			if(cc == Nchunk)then
				Nidx = rmax - (cc-1)*Chunksize
			else 
				Nidx = Chunksize
			end if
			
			!!!!!!!!!!!!!!!!!!!!!!!! get right multiplied results			
			allocate(RandVectInR(N,Nidx))
			RandVectInR=0
			allocate(RandVectOutR(M,Nidx))
			RandVectOutR=0
			
			mem_vec=0
			mem_vec =mem_vec + SIZEOF(RandVectInR)/1024.0d3
			mem_vec =mem_vec + SIZEOF(RandVectOutR)/1024.0d3
			Memory_int_vec = max(Memory_int_vec,mem_vec)			
			
			! do ii=idx_start_n_loc,idx_end_n_loc
			! do jj=1,rmax
				! RandVectInR(ii,jj)=random_complex_number()
			! end do
			! end do
			
			call RandomMat(nn,Nidx,min(nn,Nidx),RandVectInR(idx_start_n_loc:idx_end_n_loc,1:Nidx),0)	

			! write(*,*)'yani 1'
			! call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,rmax,RandVectInR,RandVectOutR,ctemp1,ctemp2)
			! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',Nidx,RandVectInR,RandVectOutR)
			call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'N',Nidx,RandVectInR,RandVectOutR)
			
			! write(*,*)'yani 2'				
			matRcol(1:nn,idx_s:idx_s+Nidx-1) = RandVectInR(idx_start_n_loc:idx_start_n_loc+nn-1,1:Nidx)	
			deallocate(RandVectInR)
			matZRcol(1:mm,idx_s:idx_s+Nidx-1) = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:Nidx)
			deallocate(RandVectOutR)	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
			
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!! get left multiplied results
			allocate(RandVectInL(M,Nidx))
			RandVectInL=0
			allocate(RandVectOutL(N,Nidx))
			RandVectOutL=0

			mem_vec=0
			mem_vec =mem_vec + SIZEOF(RandVectInL)/1024.0d3
			mem_vec =mem_vec + SIZEOF(RandVectOutL)/1024.0d3
			Memory_int_vec = max(Memory_int_vec,mem_vec)			
			! write(*,*)mem_vec,'dd2',M,N,Nidx
			! do ii=idx_start_m_loc,idx_end_m_loc
			! do jj=1,rmax
				! RandVectInL(ii,jj)=random_complex_number()
			! end do
			! end do
			
			call RandomMat(mm,Nidx,min(mm,Nidx),RandVectInL(idx_start_m_loc:idx_end_m_loc,1:Nidx),0)		
					
			! call Bplus_block_MVP_randomized_dat(bplus,'T',M,N,rmax,RandVectInL,RandVectOutL,ctemp1,ctemp2)
			! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',Nidx,RandVectInL,RandVectOutL)
			call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'T',Nidx,RandVectInL,RandVectOutL)
			
			! write(*,*)'yani 3'	
			matRrow(1:mm,idx_s:idx_s+Nidx-1) = RandVectInL(idx_start_m_loc:idx_start_m_loc+mm-1,1:Nidx)
			deallocate(RandVectInL)
			matZcRrow(1:nn,idx_s:idx_s+Nidx-1) = RandVectOutL(idx_start_n_loc:idx_start_n_loc+nn-1,1:Nidx)	
			deallocate(RandVectOutL)
		end do
		matRrow = conjg(matRrow)
		matZcRrow = conjg(matZcRrow)
					
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

		! !!!!!!!!!!!!!!!!!!!!!!!! get right multiplied results

		! allocate(RandVectInR(N,rmax))
		! RandVectInR=0
		! allocate(RandVectOutR(M,rmax))
		
		! mem_vec=0
		! mem_vec =mem_vec + SIZEOF(RandVectInR)/1024.0d3
		! mem_vec =mem_vec + SIZEOF(RandVectOutR)/1024.0d3
		! Memory_int_vec = max(Memory_int_vec,mem_vec)	
		! write(*,*)mem_vec,'ze1'
		! ! do ii=idx_start_n_loc,idx_end_n_loc
		! ! do jj=1,rmax
			! ! RandVectInR(ii,jj)=random_complex_number()
		! ! end do
		! ! end do

		! call RandomMat(nn,rmax,min(nn,rmax),RandVectInR(idx_start_n_loc:idx_end_n_loc,1:rmax),0)
		
		! ! write(*,*)'yani 1'
		! ! call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,rmax,RandVectInR,RandVectOutR,ctemp1,ctemp2)
		! call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'N',rmax,RandVectInR,RandVectOutR)
		
		! ! write(*,*)'yani 2'		
		! allocate(matRcol(nn,rmax))
		! matRcol = RandVectInR(idx_start_n_loc:idx_start_n_loc+nn-1,1:rmax)	
		! deallocate(RandVectInR)
		! allocate(matZRcol(mm,rmax))
		! matZRcol = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:rmax)
		! deallocate(RandVectOutR)

		! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
				
		! !!!!!!!!!!!!!!!!!!!!!!!!!!!! get left multiplied results
		
		! allocate(RandVectInL(M,rmax))
		! RandVectInL=0
		! allocate(RandVectOutL(N,rmax))
		
		! mem_vec=0
		! mem_vec =mem_vec + SIZEOF(RandVectInL)/1024.0d3
		! mem_vec =mem_vec + SIZEOF(RandVectOutL)/1024.0d3
		! Memory_int_vec = max(Memory_int_vec,mem_vec)		
				! write(*,*)mem_vec,'ze2',M,N,rmax
		! ! do ii=idx_start_m_loc,idx_end_m_loc
		! ! do jj=1,rmax
			! ! RandVectInL(ii,jj)=random_complex_number()
		! ! end do
		! ! end do
		
		! call RandomMat(mm,rmax,min(mm,rmax),RandVectInL(idx_start_m_loc:idx_end_m_loc,1:rmax),0)		
		
		! ! call Bplus_block_MVP_randomized_dat(bplus,'T',M,N,rmax,RandVectInL,RandVectOutL,ctemp1,ctemp2)
		! call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'T',rmax,RandVectInL,RandVectOutL)
		! ! write(*,*)'yani 3'	
		! allocate(matRrow(mm,rmax))
		! matRrow = RandVectInL(idx_start_m_loc:idx_start_m_loc+mm-1,1:rmax)
		! matRrow = conjg(matRrow)
		! deallocate(RandVectInL)
		! allocate(matZcRrow(nn,rmax))
		! matZcRrow = RandVectOutL(idx_start_n_loc:idx_start_n_loc+nn-1,1:rmax)
		! matZcRrow = conjg(matZcRrow)	
		! deallocate(RandVectOutL)

		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		allocate(matU_glo(mm,rmax))
		allocate(matV_glo(rmax,nn))
		allocate(Singular(rmax))

		call RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU_glo,matV_glo,Singular,mm,nn,rmax,rank,LS_tolerance,SVD_tolerance_factor)				
		rankthusfarBC = max(rankthusfarBC,rank)
		! write(*,*)'yani 4'
		do ii=1,rank
			matV_glo(ii,:) = matV_glo(ii,:) * Singular(ii)
		end do
		
		deallocate(matRcol,matZRcol,matRrow,matZcRrow,Singular)



		if(mm*nn*16/1e9>5)then
			write(*,*)'warning: full storage of matSub_glo is too costly'
		end if
		
																				  
							   

		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! test the results
		allocate(Vin(nn,1))
		allocate(Vout1(mm,1))
		allocate(Vout2(mm,1))
		allocate(Vout3(rank,1))
		do ii=1,nn
			Vin(ii,1)=random_complex_number()
		end do
		
		
		allocate(RandVectInR(N,1))
		allocate(RandVectOutR(M,1))
		RandVectInR=0
		RandVectOutR=0
		RandVectInR(idx_start_n_loc:idx_end_n_loc,1:1) = Vin
		! call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,1,RandVectInR,RandVectOutR,ctemp1,ctemp2)
		! write(*,*)'yani 4.1'
		call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'N',1,RandVectInR,RandVectOutR)
		! write(*,*)'yani 5'
		Vout1 = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:1)
		deallocate(RandVectInR)
		deallocate(RandVectOutR)


		call gemm_omp(matV_glo(1:rank,1:nn),Vin,Vout3,rank,nn,1)
		call gemm_omp(matU_glo(1:mm,1:rank),Vout3,Vout2,mm,rank,1)

		
		error = fnorm(Vout2-Vout1,mm,1)/fnorm(Vout1,mm,1)
		! write(*,*)error,bb_o,'ninini'		
		
		deallocate(Vin)
		deallocate(Vout1)
		deallocate(Vout2)
		deallocate(Vout3)
		
		! write(*,*)tt,rmax,'minusBC',error,rank
		if(error>iter_tolerance*1d-1)then
			deallocate(matU_glo)
			deallocate(matV_glo)
			write(*,*)tt,rmax,'minusBC',error,rank
			if(min(mm,nn)==rmax)then
			!	write(*,*)tt,rmax,'minusBC',error,rank
				write(*,*)'no need to increase rmax, try increase RandomizedSVD tolerance'
				stop
			end if
		else
			if(verboselevel>=2)write(*,'(A33,A15,I3,A8,I4,A8,I2,A7,Es14.7)')' Onesub ',' bb_o:',bb_o,' rank:',rank,'Ntrial:',tt,' error:',error	
			! write(*,*)bb_o,error,rmax,rank,'good'
			error_avr_glo = error_avr_glo + error
			error_cnt = error_cnt + 1	
			exit
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	end do

	if(error>iter_tolerance*1d-1)then
		write(*,*)'matSub_glo no correct, needs more work'
		stop
	end if
  
	!!!! construct butterfly at all subsequent levels including level 2 (this part is generic and doesn't require accuracy test)
	idx_start_n_ref = idx_start_n_loc + idx_start_n - 1
	idx_end_n_ref = idx_end_n_loc + idx_start_n - 1
	idx_start_m_ref = idx_start_m_loc + idx_start_m - 1
	idx_end_m_ref = idx_end_m_loc + idx_start_m - 1

	if(TwoLayerOnly==1)then
		ll=2
		allocate(boxindex(Bplus_randomized(1)%LL(ll)%Nbound))
		Nactive = 0
		do bb = 1,Bplus_randomized(1)%LL(ll)%Nbound
			if(basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
				Nactive = Nactive + 1
				boxindex(Nactive) = bb
			end if
		end do
		call assert(Nactive==1,'Nactive should be one')
		

		bb = boxindex(1)
 
																										
									  
														
																				 
		blocks => Bplus_randomized(1)%LL(ll)%matrices_block(bb)
		
		blocks%level_butterfly=0
		allocate(blocks%ButterflyU(1))
		allocate(blocks%ButterflyV(1))
		blocks%rankmax = rank
		blocks%rankmin = rank
		
		allocate (blocks%ButterflyV(1)%matrix(nn,rank))
		call copymatT_OMP(matV_glo(1:rank,1:nn),blocks%ButterflyV(1)%matrix,rank,nn)
		allocate (blocks%ButterflyU(1)%matrix(mm,rank))
		call copymatN_OMP(matU_glo(1:mm,1:rank),blocks%ButterflyU(1)%matrix,mm,rank)
		deallocate(matV_glo)
		deallocate(matU_glo)

		do ii = 1,Nactive
			bb = boxindex(ii)
			Bplus_randomized(1)%LL(ll)%rankmax = max(Bplus_randomized(1)%LL(ll)%rankmax,Bplus_randomized(1)%LL(ll)%matrices_block(bb)%rankmax)
		end do
		
		deallocate(boxindex)		
			
	else 
	
		allocate(matSub_glo(mm,nn))
		n1 = OMP_get_wtime()
		call gemm_omp(matU_glo(1:mm,1:rank),matV_glo(1:rank,1:nn),matSub_glo,mm,rank,nn)
		deallocate(matU_glo,matV_glo)
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1
			
		n1 = OMP_get_wtime()
		do ll=2,Bplus_randomized(1)%Lplus
		
			! level_butterfly = int((Maxlevel_for_blocks - Bplus_randomized(1)%LL(ll)%matrices_block(1)%level)/2)*2 
			! level_BP = Bplus_randomized(1)%level
			! levelm = ceiling_safe(dble(level_butterfly)/2d0)						
			! groupm_start=Bplus_randomized(1)%LL(ll)%matrices_block(1)%row_group*2**levelm		
			! Nboundall = 2**(Bplus_randomized(1)%LL(ll)%matrices_block(1)%level+levelm-level_BP)
			
			
			allocate(boxindex(Bplus_randomized(1)%LL(ll)%Nbound))
			Nactive = 0
			do bb = 1,Bplus_randomized(1)%LL(ll)%Nbound
				if(basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
					Nactive = Nactive + 1
					boxindex(Nactive) = bb
				end if
			end do
			
			!$omp parallel do default(shared) private(bb,ii)
			do ii = 1,Nactive
				! if(basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
				bb = boxindex(ii)		
				if(Bplus_randomized(1)%LL(ll+1)%Nbound==0)then
					! write(*,*)'666',ll
					call Butterfly_compress_N15_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),idx_start_m_ref,idx_start_n_ref)
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
				else
					! write(*,*)'777'				
					call Butterfly_compress_N15_withoutBoundary_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),Bplus_randomized(1)%LL(ll+1)%boundary_map,Nboundall,groupm_start, rtemp, idx_start_m_ref,idx_start_n_ref)
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
				end if				
			end do
			!$omp end parallel do
			
			do ii = 1,Nactive
				bb = boxindex(ii)
				Bplus_randomized(1)%LL(ll)%rankmax = max(Bplus_randomized(1)%LL(ll)%rankmax,Bplus_randomized(1)%LL(ll)%matrices_block(bb)%rankmax)
			end do
			
			deallocate(boxindex)			
			
		! write(*,*)Bplus_randomized(1)%LL(ll)%rankmax,'qunima'
		end do
														
	   

		deallocate(matSub_glo)
		
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1
	end if
	
    return

end subroutine MultiLrandomized_Onesubblock_minusBC






subroutine Bplus_block_MVP_minusBC_dat(level_c,rowblock,trans,num_vect_sub,Vin,Vout)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2
	type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	
	real*8::n2,n1 	
	
	bplus_off1 => cascading_factors(level_c)%BP(2*rowblock-1)	
	bplus_off2 => cascading_factors(level_c)%BP(2*rowblock)	
	
	
	if(trans=='N')then

		! allocate (RandomVectors_InOutput(3))
		groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		! allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
		! allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
		! allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		! do ii =1,3
			! RandomVectors_InOutput(ii)%vector = 0
		! end do	 
		! RandomVectors_InOutput(1)%vector = Vin
   
		allocate(vec_new(nn,num_vect_sub))
		vec_new = 0

		! get the right multiplied vectors
		! random1=>RandomVectors_InOutput(1)
		! random2=>RandomVectors_InOutput(2)
		ctemp1=1.0d0 ; ctemp2=0.0d0
		! write(*,*)'hah1',CheckNAN_Bplus(bplus_off2)
		call Bplus_block_MVP_randomized_dat(bplus_off2,'N',nn,mm,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
		! random1=>RandomVectors_InOutput(2)
		! random2=>RandomVectors_InOutput(3)		
		! write(*,*)'hah2',CheckNAN_Bplus(bplus_off1)
		call Bplus_block_MVP_randomized_dat(bplus_off1,'N',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
		Vout = -Vout
		! write(*,*)'hah3'
		deallocate(vec_new)
		! Vout = RandomVectors_InOutput(3)%vector
	   
		! ! !$omp parallel do default(shared) private(i)
		! do i=1, 3
			! deallocate (RandomVectors_InOutput(i)%vector)
		! enddo
		! ! !$omp end parallel do
		! deallocate (RandomVectors_InOutput)		   
   
   else if(trans=='T')then
		! allocate (RandomVectors_InOutput(3))
		groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		! allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
		! allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
		! allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		! do ii =1,3
			! RandomVectors_InOutput(ii)%vector = 0
		! end do	 
		! RandomVectors_InOutput(1)%vector = Vin
		
		allocate(vec_new(nn,num_vect_sub))
		vec_new=0

		! get the right multiplied vectors
		! random1=>RandomVectors_InOutput(1)
		! random2=>RandomVectors_InOutput(2)
		ctemp1=1.0d0 ; ctemp2=0.0d0
		call Bplus_block_MVP_randomized_dat(bplus_off1,'T',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
		! random1=>RandomVectors_InOutput(2)
		! random2=>RandomVectors_InOutput(3)		
		call Bplus_block_MVP_randomized_dat(bplus_off2,'T',nn,mm,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
		Vout = -Vout
		deallocate(vec_new)
		! Vout = RandomVectors_InOutput(3)%vector
	   
		! ! !$omp parallel do default(shared) private(i)
		! do i=1, 3
			! deallocate (RandomVectors_InOutput(i)%vector)
		! enddo
		! ! !$omp end parallel do
		! deallocate (RandomVectors_InOutput)		
 
   end if
   
 
end subroutine Bplus_block_MVP_minusBC_dat







subroutine MultiLrandomized_Outter_minusBC_memfree(level_c,rowblock,Memory)

    use MODULE_FILE
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
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,rank,ntry
	real*8:: error_inout
	real*8:: n1,n2,Memory
	type(blockplus),pointer::Bplus
	
	Memory = 0
	
	Bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock)
    ! block_o =>  Bplus%LL(1)%matrices_block(1)
	
	level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2

	do tt =1,10
		do ntry=1,1
		! write(*,*)'haha1'
		n1 = OMP_get_wtime()
		call MultiLInitialize_Butterfly_inverse_BC(level_c,rowblock,tt-1)
		n2 = OMP_get_wtime()
		! Time_Init_forward = Time_Init_forward + n2 -n1 
		
! write(*,*)'haha2'
		n1 = OMP_get_wtime()
		call MultiLReconstruction_LL_Outter_minusBC(level_c,rowblock)	
		! write(*,*)'haha3'
		call MultiLReconstruction_RR_Outter_minusBC(level_c,rowblock,error_inout)
		n2 = OMP_get_wtime()
		Time_Reconstruct_inverse = Time_Reconstruct_inverse + n2-n1
		! time_tmp = time_tmp + n2 - n1
! write(*,*)'haha4'
		! write(*,*)tt,error_inout
		
		if(error_inout>iter_tolerance)then
		! if(0)then			
			call Delete_randomized_butterfly()
		else 
			call delete_blocks(Bplus_randomized(1)%LL(1)%matrices_block(1))
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			rank_new_max = butterfly_block_randomized(1)%rankmax				
			call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),Bplus_randomized(1)%LL(1)%matrices_block(1),Memory)
			deallocate(butterfly_block_randomized)
			! call copy_randomizedbutterfly(butterfly_block_randomized(1),Bplus_randomized(1)%LL(1)%matrices_block(1),Memory)
			! call Delete_randomized_butterfly()
			Bplus_randomized(1)%LL(1)%rankmax =rank_new_max
									 
			if(verboselevel>=2)write(*,'(A36,A15,I3,A8,I2,A8,I3,A7,Es14.7)')' OutterBC: ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
			error_avr_glo = error_avr_glo + error_inout
			error_cnt = error_cnt + 1	
			return
		end if
		end do
	end do
	write(*,*)'randomized scheme not converged in MultiLrandomized_Outter_minusBC_memfree',error_inout
	stop
	
    return

end subroutine MultiLrandomized_Outter_minusBC_memfree

subroutine MultiLInitialize_Butterfly_inverse_BC(level_c,rowblock,kover)

    use MODULE_FILE
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn,ll
    integer dimension_rank, dimension_rank_dummy, dimension_max,dimension_m, dimension_n, blocks, groupm, groupn,index_j,index_i,rankmax1,rankmax2
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
 	
	
	rankmax1 = 0
	do ll=1,1!cascading_factors(level_c)%BP(rowblock*2-1)%Lplus
		rankmax1 = max(rankmax1, cascading_factors(level_c)%BP(rowblock*2-1)%LL(ll)%rankmax)
		! write(*,*)level_c,rowblock*2-1,ll,cascading_factors(level_c)%BP(rowblock*2-1)%LL(ll)%rankmax,'d1'
	end do
	rankmax2 = 0
	do ll=1,1!cascading_factors(level_c)%BP(rowblock*2)%Lplus
		rankmax2 = max(rankmax2, cascading_factors(level_c)%BP(rowblock*2)%LL(ll)%rankmax)
		! write(*,*)level_c,rowblock*2,ll,cascading_factors(level_c)%BP(rowblock*2)%LL(ll)%rankmax,'d2'
	end do	
	
																									 
 
 
	! dimension_rank= max(rankmax1,rankmax2)+kover  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= max(rankmax1,rankmax2)*1.2d0**(kover)   !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank_dummy = 1
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
		! write(*,*)dimension_rank,dimension_m,'1a'
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank_dummy))

		allocate(matrixtemp1(dimension_rank_dummy,dimension_m))
		call RandomMat(dimension_rank_dummy,dimension_m,min(dimension_m,dimension_rank_dummy),matrixtemp1,0)
        
		! !$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank_dummy
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
		dimension_n=size(block_off2%ButterflyV(blocks)%matrix,1)
		! write(*,*)dimension_rank,dimension_n,'2a'
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank_dummy))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        
		allocate(matrixtemp1(dimension_rank_dummy,dimension_n))
		call RandomMat(dimension_rank_dummy,dimension_n,min(dimension_n,dimension_rank_dummy),matrixtemp1,0)
		
		! !$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank_dummy
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
		! write(*,*)dimension_rank,level,level_butterfly,'1c'
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))

        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine MultiLInitialize_Butterfly_inverse_BC
subroutine MultiLReconstruction_LL_Outter_minusBC(level_c,rowblock)
    
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
    type(matrixblock),pointer::block_off1
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)

    block_off1 =>  cascading_factors(level_c)%BP(2*rowblock-1)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_off1%level)/2)*2
	
	
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
 
    ! ! allocate(perms(Nsub))
	! ! call rperm(Nsub, perms)
	! ! do ii = 1,Nsub		
		! ! nth_s = perms(ii)
		! ! nth_e = perms(ii)
	
	do unique_nth = 0,level_right_start
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()
			call MultiLGet_Randomized_Vectors_LL_Outter_minusBC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1		
		end do
	end do
	
	! pause
	! call Resolving_Butterfly_LL_rankcompletion()
	
	! deallocate(perms)
	
	
	random=>Random_Block(1)
	call Delete_RandVect('T',random,level_butterfly)

    return
    
end subroutine MultiLReconstruction_LL_Outter_minusBC




subroutine MultiLReconstruction_RR_Outter_minusBC(level_c,rowblock,error)
    
    use MODULE_FILE
    implicit none
	   
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_left_start,num_row,num_col
    real*8::n1,n2,error
	integer level_c,rowblock
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    ! type(blockplus)::Bplus
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	integer,allocatable::perms(:)
	type(matrixblock),pointer::block_off1
		
    block_off1 =>  cascading_factors(level_c)%BP(2*rowblock-1)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_off1%level)/2)*2
		
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

    ! ! allocate(perms(Nsub))
	! ! call rperm(Nsub, perms)
	! ! do ii = 1,Nsub		
		! ! nth_s = perms(ii)
		! ! nth_e = perms(ii)
	
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
			call MultiLGet_Randomized_Vectors_RR_Outter_minusBC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do
	
	! deallocate(perms)

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call MultiLTest_Error_RR_Outter_minusBC(level_c,rowblock,error)


	
    return
    
end subroutine MultiLReconstruction_RR_Outter_minusBC




subroutine MultiLTest_Error_RR_Outter_minusBC(level_c,rowblock,error)

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
	integer dimension_m 
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(blockplus),pointer::Bplus
	integer level_c,rowblock

	Bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock) 
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
	num_blocks=2**level_butterfly
	! write(*,*)level_butterfly,'heiyou',maxlevel_for_blocks,block_o%level
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call MultiLGet_Randomized_Vectors_RR_Test_Outter_minusBC(level_c,rowblock,num_vect)

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

end subroutine MultiLTest_Error_RR_Outter_minusBC




subroutine MultiLGet_Randomized_Vectors_LL_Outter_minusBC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
	integer unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2, ctemp3, ctemp4
	type(blockplus),pointer::Bplus,Bplus_off1,Bplus_off2
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	complex(kind=8),allocatable::Vin(:,:),Vout(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_right_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	
    ctemp1=1.0d0 ; ctemp2=0.0d0	
    ctemp3=-1.0d0 ; ctemp4=1.0d0	
	
	Bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    ! allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	allocate (Vin(mm,num_vect_sub))
    allocate (Vout(mm,num_vect_sub))
	Vin = 0
	Vout = 0
	
 
	 
	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		! do i=1, num_blocks
		!$omp parallel do default(shared) private(i,header_m,tailer_m,mm,k)			
		do i=(nth-1)*Ng+1, nth*Ng		
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_m=basis_group(groupm_start+i-1)%head
				tailer_m=basis_group(groupm_start+i-1)%tail
				mm=tailer_m-header_m+1
				k=header_m-header_mm	

				! allocate(matrixtemp1(num_vect_subsub,mm))
				call RandomMat(mm,num_vect_subsub,min(mm,num_vect_subsub),Vin(1+k:mm+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, mm
						 ! Vin(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 ! deallocate(matrixtemp1)
				 
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	idx_start_glo = basis_group(groupm)%head		

	
	
	n1 = OMP_get_wtime()
    ! call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
	
	call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'T',num_vect_sub,Vin,Vout)	
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'T',mm,mm,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1		
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=Vin(ii+k,jj)
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
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=Vout(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo
	deallocate(Vin)
	deallocate(Vout)

	
    return                

end subroutine MultiLGet_Randomized_Vectors_LL_Outter_minusBC






subroutine MultiLGet_Randomized_Vectors_RR_Outter_minusBC(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2, ctemp3, ctemp4
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    type(blockplus),pointer::bplus
	
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	complex(kind=8),allocatable::Vin(:,:),Vout(:,:)
		
	integer Nsub,Ng,unique_nth,level_left_start
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub
	type(RandomBlock), pointer :: random
	real*8::n2,n1

	Bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock) 	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    ! allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	level_left_start= floor_safe(level_butterfly/2d0)+1
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub
	
	
	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	! allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	! allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
    ! allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	! do ii =1,3
		! RandomVectors_InOutput(ii)%vector = 0
	! end do	 
	 
	allocate(Vin(mm,num_vect_sub)) 
	allocate(Vout(mm,num_vect_sub)) 
	Vin = 0
	Vout = 0
	 
	groupn_start=groupn*2**(level_butterfly)
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i,header_n,tailer_n,nn,k)	
		do i=(nth-1)*Ng+1, nth*Ng	
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_n=basis_group(groupn_start+i-1)%head
				tailer_n=basis_group(groupn_start+i-1)%tail
				nn=tailer_n-header_n+1
				k=header_n-header_nn

				! allocate(matrixtemp1(num_vect_subsub,nn))
				call RandomMat(nn,num_vect_subsub,min(nn,num_vect_subsub),Vin(1+k:nn+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
				
				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, nn
						 ! Vin(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	
	! get the right multiplied vectors
	idx_start_glo = basis_group(groupm)%head
    ! random1=>RandomVectors_InOutput(1)
    ! random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    ctemp3=-1.0d0 ; ctemp4=1.0d0
n1 = OMP_get_wtime()  
  
  ! call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)
  call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'N',num_vect_sub,Vin,Vout)
  call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,mm,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
  n2 = OMP_get_wtime()
! time_tmp = time_tmp + n2 - n1	
   


	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=Vin(ii+k,jj)
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
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=Vout(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	deallocate(Vin)
	deallocate(Vout)
	
    return             

end subroutine MultiLGet_Randomized_Vectors_RR_Outter_minusBC
subroutine MultiLGet_Randomized_Vectors_RR_Test_Outter_minusBC(level_c,rowblock,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2,  ctemp3, ctemp4
	type(blockplus),pointer::Bplus
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	complex(kind=8),allocatable::Vin(:,:),Vout(:,:)
		
	
	Bplus =>  cascading_factors(level_c)%BP_inverse_schur(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	! allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	! allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
    ! allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	! do ii =1,3
		! RandomVectors_InOutput(ii)%vector = 0
	! end do	 
	allocate(Vin(mm,num_vect_sub)) 
	allocate(Vout(mm,num_vect_sub)) 
	Vin = 0
	Vout = 0
	
	groupn_start=groupn*2**(level_butterfly)
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	
	do i=1, num_blocks
		header_n=basis_group(groupn_start+i-1)%head
		tailer_n=basis_group(groupn_start+i-1)%tail
		nn=tailer_n-header_n+1
		k=header_n-header_nn

		! !$omp parallel do default(shared) private(ii,jj)
		 do jj=1,num_vect_sub
			 do ii=1, nn
				 Vin(ii+k,jj)=random_complex_number()	
			 enddo
		 enddo
		 ! !$omp end parallel do
	end do
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	
	! get the right multiplied vectors
	
    ! random1=>RandomVectors_InOutput(1)
    ! random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    ctemp3=-1.0d0 ; ctemp4=1.0d0
	call Bplus_block_MVP_minusBC_dat(level_c,rowblock,'N',num_vect_sub,Vin,Vout)
	! call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)		
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,mm,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	 


	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=Vin(ii+k,jj)
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
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=Vout(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
	deallocate(Vin)
	deallocate(Vout)
	
	
    return                

end subroutine MultiLGet_Randomized_Vectors_RR_Test_Outter_minusBC






subroutine MultiL_Reconstruction_LL_SchurSmall(bplus,block_o,LR)
    
    use MODULE_FILE
    implicit none
	
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara,LR
	integer level_right_start,num_col,num_row
    integer level_c,rowblock,ll_o,ij_loc
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock)::block_o
	type(blockplus)::bplus
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)

	level_butterfly=block_o%level_butterfly
	
	
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
 
    ! ! allocate(perms(Nsub))
	! ! call rperm(Nsub, perms)
	! ! do ii = 1,Nsub		
		! ! nth_s = perms(ii)
		! ! nth_e = perms(ii)
	
	do unique_nth = 0,level_right_start
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()
			! call Get_Randomized_Vectors_LL_SblockSmall(level_c,ii_inv,block_o,nth_s,nth_e,num_vect_sub,unique_nth)
			call MultiLGet_Randomized_Vectors_LL_SchurSmall(bplus,block_o,nth_s,nth_e,num_vect_sub,unique_nth,LR)

			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1		
		end do
	end do
	
	! pause
	! call Resolving_Butterfly_LL_rankcompletion()
	
	! deallocate(perms)
	
	
	random=>Random_Block(1)
	call Delete_RandVect('T',random,level_butterfly)

    return
    
end subroutine MultiL_Reconstruction_LL_SchurSmall








subroutine MultiL_Reconstruction_RR_SchurSmall(bplus,block_o,error,LR)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,ii_inv    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara,LR
	integer level_left_start,num_row,num_col
    real*8::n1,n2,error
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock)::block_o
    type(blockplus)::bplus
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	integer,allocatable::perms(:)
		
	level_butterfly=block_o%level_butterfly
		
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

    ! ! allocate(perms(Nsub))
	! ! call rperm(Nsub, perms)
	! ! do ii = 1,Nsub		
		! ! nth_s = perms(ii)
		! ! nth_e = perms(ii)
	
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
			call MultiLGet_Randomized_Vectors_RR_SchurSmall(bplus,block_o,nth_s,nth_e,num_vect_sub,unique_nth,LR)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do
	
	! deallocate(perms)

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call MultiLTest_Error_RR_SchurSmall(bplus,block_o,error,LR)


	
    return
    
end subroutine MultiL_Reconstruction_RR_SchurSmall



! LR: L: bplus*block_o, R: block_o*bplus
subroutine MultiLGet_Randomized_Vectors_LL_SchurSmall(bplus,block_o,nth_s,nth_e,num_vect_sub,unique_nth,LR)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,mm_loc,nn_loc
    character chara,LR
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	type(blockplus)::bplus
	integer level_c, rowblock, ll_o
	
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
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_right_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	
    ctemp1=1.0d0 ; ctemp2=0.0d0	

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    level_butterfly=block_o%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    mm=0
	nn=0	
	do i=1, num_blocks
		mm = mm + size(block_o%butterflyU(i)%matrix,1)
		nn = nn + size(block_o%butterflyV(i)%matrix,1)
	end do
 	
	! call assert(mm==nn,'this should be a square matrix in MultiLGet_Randomized_Vectors_LL_SchurSmall')
	
	
	
	
	if(LR=='L')then
	
	
		allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub))    
		allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		
		
		do nth= nth_s,nth_e
			k = 0
			do i=1, num_blocks
				mm_loc=size(block_o%butterflyU(i)%matrix,1)	
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					! allocate(matrixtemp1(num_vect_subsub,mm_loc))
					call RandomMat(mm_loc,num_vect_subsub,min(mm_loc,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:mm_loc+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
	 
					! ! !$omp parallel do default(shared) private(ii,jj)
					 ! do jj=1,num_vect_subsub
						 ! do ii=1, mm_loc
							 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
						 ! enddo
					 ! enddo
					 ! ! !$omp end parallel do
					 ! deallocate(matrixtemp1)
				 end if
				 k=k+mm_loc
			end do
		end do
		
		n1 = OMP_get_wtime()
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'T',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'T',mm,mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(2)%vector = RandomVectors_InOutput(1)%vector + RandomVectors_InOutput(2)%vector
		
		random1=>RandomVectors_InOutput(2)
		random2=>RandomVectors_InOutput(3)
		call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1		
	
	else if(LR=='R')then
		allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub))    
		allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		
		
		do nth= nth_s,nth_e
			k = 0
			do i=1, num_blocks
				mm_loc=size(block_o%butterflyU(i)%matrix,1)	
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					! allocate(matrixtemp1(num_vect_subsub,mm_loc))
					call RandomMat(mm_loc,num_vect_subsub,min(mm_loc,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:mm_loc+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
	 
					! ! !$omp parallel do default(shared) private(ii,jj)
					 ! do jj=1,num_vect_subsub
						 ! do ii=1, mm_loc
							 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
						 ! enddo
					 ! enddo
					 ! ! !$omp end parallel do
					 ! deallocate(matrixtemp1)
				 end if
				 k=k+mm_loc
			end do
		end do
		
		n1 = OMP_get_wtime()
		random1=>RandomVectors_InOutput(1)
		random2=>RandomVectors_InOutput(2)
		call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
		
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'T',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'T',nn,nn,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector + RandomVectors_InOutput(3)%vector
		

		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1			
	end if
	

	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		mm_loc=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm_loc
			do jj=1, num_vect_sub
				random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm_loc
	enddo 
	
	k=0
	do i=1, num_blocks
		nn_loc=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn_loc
			do jj=1, num_vect_sub
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn_loc
	enddo 	

    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine MultiLGet_Randomized_Vectors_LL_SchurSmall





! LR: L: bplus*block_o, R: block_o,bplus
subroutine MultiLGet_Randomized_Vectors_RR_SchurSmall(bplus,block_o,nth_s,nth_e,num_vect_sub,unique_nth,LR)


    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,ii_inv
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,mm_loc,nn_loc
    character chara,LR
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	type(blockplus)::bplus
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng,unique_nth,level_left_start
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub
	type(RandomBlock), pointer :: random
	real*8::n2,n1
	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	  
    level_butterfly=block_o%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    mm=0
	nn=0	
	do i=1, num_blocks
		mm = mm + size(block_o%butterflyU(i)%matrix,1)
		nn = nn + size(block_o%butterflyV(i)%matrix,1)
	end do

	level_left_start= floor_safe(level_butterfly/2d0)+1
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	
	if(LR=='L')then
	
		allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		 
		do nth= nth_s,nth_e
			k=0
			do i=1, num_blocks
				nn_loc=size(block_o%butterflyV(i)%matrix,1)	
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	

					! allocate(matrixtemp1(num_vect_subsub,nn_loc))
					call RandomMat(nn_loc,num_vect_subsub,min(nn_loc,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:nn_loc+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
					
					! ! !$omp parallel do default(shared) private(ii,jj)
					 ! do jj=1,num_vect_subsub
						 ! do ii=1, nn_loc
							 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
						 ! enddo
					 ! enddo
					 ! ! !$omp end parallel do
					 
					 ! deallocate(matrixtemp1)
				 end if
				 k=k+nn_loc
			end do		
		end do
		
		n1 = OMP_get_wtime() 
		ctemp1=1.0d0 ; ctemp2=0.0d0
		random1=>RandomVectors_InOutput(1)
		random2=>RandomVectors_InOutput(2)
		call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
		
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'N',mm,mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(3)%vector + RandomVectors_InOutput(2)%vector
		n2 = OMP_get_wtime()	
		! time_tmp = time_tmp + n2 - n1		
	
	else if(LR=='R')then
		allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		 
		do nth= nth_s,nth_e
			k=0
			do i=1, num_blocks
				nn_loc=size(block_o%butterflyV(i)%matrix,1)	
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	

					! allocate(matrixtemp1(num_vect_subsub,nn_loc))
					call RandomMat(nn_loc,num_vect_subsub,min(nn_loc,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:nn_loc+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
	 
					! ! !$omp parallel do default(shared) private(ii,jj)
					 ! do jj=1,num_vect_subsub
						 ! do ii=1, nn_loc
							 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
						 ! enddo
					 ! enddo
					 ! ! !$omp end parallel do
					 
					 ! deallocate(matrixtemp1)
				 end if
				 k=k+nn_loc
			end do		
		end do
		
		n1 = OMP_get_wtime() 
		ctemp1=1.0d0 ; ctemp2=0.0d0
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'N',nn,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(2)%vector = RandomVectors_InOutput(2)%vector + RandomVectors_InOutput(1)%vector
		
		random1=>RandomVectors_InOutput(2)
		random2=>RandomVectors_InOutput(3)
		call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
		

		n2 = OMP_get_wtime()	
		! time_tmp = time_tmp + n2 - n1		
	end if
	
	
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn_loc=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn_loc
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn_loc
	enddo 

	k=0
	do i=1, num_blocks
		mm_loc=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm_loc
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm_loc
	enddo 
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine MultiLGet_Randomized_Vectors_RR_SchurSmall


! LR: L: bplus*block_o, R: block_o,bplus
subroutine MultiLGet_Randomized_Vectors_RR_Test_SchurSmall(bplus,block_o,num_vect_sub,LR)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,mm_loc,nn_loc
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara,LR
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	type(blockplus)::bplus
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	

    level_butterfly=block_o%level_butterfly
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    mm=0
	nn=0	
	do i=1, num_blocks
		mm = mm + size(block_o%butterflyU(i)%matrix,1)
		nn = nn + size(block_o%butterflyV(i)%matrix,1)
	end do
	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	if(LR=='L')then
	
		allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		 
		k=0
		do i=1, num_blocks
			nn_loc=size(block_o%butterflyV(i)%matrix,1)	
			! !$omp parallel do default(shared) private(ii,jj)
			 do jj=1,num_vect_sub
				 do ii=1, nn_loc
					 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number() 
				 enddo
			 enddo
			 ! !$omp end parallel do
			 k=k+nn_loc
		end do	


		ctemp1=1.0d0 ; ctemp2=0.0d0
		random1=>RandomVectors_InOutput(1)
		random2=>RandomVectors_InOutput(2)
		call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'N',mm,mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(3)%vector + RandomVectors_InOutput(2)%vector
		
	else if (LR=='R')then
		allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
		allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput(ii)%vector = 0
		end do	 
		 
		k=0
		do i=1, num_blocks
			nn_loc=size(block_o%butterflyV(i)%matrix,1)	
			! !$omp parallel do default(shared) private(ii,jj)
			 do jj=1,num_vect_sub
				 do ii=1, nn_loc
					 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number() 
				 enddo
			 enddo
			 ! !$omp end parallel do
			 k=k+nn_loc
		end do	


		ctemp1=1.0d0 ; ctemp2=0.0d0
		! call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
		call Bplus_block_MVP_randomized_dat(bplus,'N',nn,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ctemp1,ctemp2)
		RandomVectors_InOutput(2)%vector = RandomVectors_InOutput(2)%vector + RandomVectors_InOutput(1)%vector
	
		random1=>RandomVectors_InOutput(2)
		random2=>RandomVectors_InOutput(3)	
		call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)

	end if	
		
	k=0
	random=>random_Block(1)
	do i=1, num_blocks
		nn_loc=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn_loc
			do jj=1, num_vect_sub
				random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn_loc
	enddo 

	k=0
	do i=1, num_blocks
		mm_loc=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm_loc
			do jj=1, num_vect_sub
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm_loc
	enddo 
	
    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine MultiLGet_Randomized_Vectors_RR_Test_SchurSmall





subroutine Extract_partial_Bplus(bplus_i,ll_s,row_group)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb,bb_o
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::rtemp
integer row_group,ll_s,idx_s,idx_e

call assert(bplus_i%row_group==bplus_i%col_group,'only works for square matrix')

idx_s = basis_group(row_group)%head
idx_e = basis_group(row_group)%tail

allocate(agent_bplus(1))
allocate(agent_bplus(1)%LL(LplusMax))
do ll=1,LplusMax
agent_bplus(1)%LL(ll)%Nbound = 0
end do

agent_bplus(1)%Lplus = bplus_i%Lplus - ll_s + 1
agent_bplus(1)%row_group = 	row_group
agent_bplus(1)%col_group = 	row_group
agent_bplus(1)%level = basis_group(row_group)%level

do ll=1,agent_bplus(1)%Lplus
	agent_bplus(1)%LL(ll)%Nbound = 0
	agent_bplus(1)%LL(ll)%rankmax=bplus_i%LL(ll+ll_s-1)%rankmax
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			agent_bplus(1)%LL(ll)%Nbound = agent_bplus(1)%LL(ll)%Nbound + 1
		end if
	end do
	if(agent_bplus(1)%LL(ll)%Nbound>0)then
		allocate(agent_bplus(1)%LL(ll)%matrices_block(agent_bplus(1)%LL(ll)%Nbound))
	end if
end do


do ll=1,agent_bplus(1)%Lplus
	bb_o = 0
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			bb_o = bb_o + 1
			call copy_butterfly(bplus_i%LL(ll+ll_s-1)%matrices_block(bb),agent_bplus(1)%LL(ll)%matrices_block(bb_o))
		end if
	end do
end do



end subroutine Extract_partial_Bplus




subroutine MultiLTest_Error_RR_SchurSmall(bplus,block_o,error,LR)

    use MODULE_FILE
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn
    real*8 a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    complex(kind=8) ctemp
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real*8::error
	integer rowblock,dimension_m 
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock)::block_o
	type(blockplus)::bplus
	character LR
	
	level_butterfly=block_o%level_butterfly
	num_blocks=2**level_butterfly
	! write(*,*)level_butterfly,'heiyou',maxlevel_for_blocks,block_o%level
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	

    mm=0
	nn=0	
	do i=1, num_blocks
		mm = mm + size(block_o%butterflyU(i)%matrix,1)
		nn = nn + size(block_o%butterflyV(i)%matrix,1)
	end do	
	
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call MultiLGet_Randomized_Vectors_RR_Test_SchurSmall(bplus,block_o,num_vect,LR)

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

end subroutine MultiLTest_Error_RR_SchurSmall




end module Bplus_inversion_schur_partition
