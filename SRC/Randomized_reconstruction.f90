module Randomized_reconstruction
use Utilites_randomized
use Butterfly_rightmultiply
contains 

subroutine Butterfly_randomized(level_butterfly,rank0,rankrate,blocks_o,operand,blackbox_MVP_dat,error_inout,strings,option,stats,operand1) 

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
	integer idx_start_m_ref
	class(*):: operand
	class(*),optional:: operand1
	character(*)  :: strings	
	type(matrixblock),allocatable::block_rand(:)
	type(Hoption)::option
	type(Hstat)::stats
	procedure(BF_MVP_blk)::blackbox_MVP_dat
	
	ctemp1 = 1d0; ctemp2 = 0d0
	Memory = 0
	
	do tt = 1,10
	
		rank_pre_max = ceiling_safe(rank0*rankrate**(tt-1))
	
		n1 = OMP_get_wtime()
		groupm=blocks_o%row_group
		groupn=blocks_o%col_group
		
		allocate (block_rand(1))
		
		call Initialize_Butterfly_randomized(level_butterfly,rank_pre_max,groupm,groupn,blocks_o,block_rand(1))
		n2 = OMP_get_wtime()
		stats%Time_random(1) = stats%Time_random(1) + n2-n1
		
		n1 = OMP_get_wtime()

		call Reconstruction_LL(block_rand(1),blocks_o,operand,blackbox_MVP_dat,operand1,option,stats)	
		call Reconstruction_RR(block_rand(1),blocks_o,operand,blackbox_MVP_dat,operand1,option,stats)
		call Test_Reconstruction_Error(block_rand(1),blocks_o,operand,blackbox_MVP_dat,error_inout,operand1)
		n2 = OMP_get_wtime()	
		
		
		call get_butterfly_minmaxrank(block_rand(1))
		
#if PRNTlevel >= 2			
		write(*,'(A38,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' '//TRIM(strings)//' ',' rank:',block_rand(1)%rankmax,' Ntrial:',tt,' L_butt:',block_rand(1)%level_butterfly,' error:',error_inout
#endif	

		if(error_inout>option%tol_rand)then
		! if(0)then
			call get_butterfly_minmaxrank(block_rand(1))
			rank_new_max = block_rand(1)%rankmax
			call delete_blocks(block_rand(1))
			deallocate(block_rand)
		else 
			call delete_blocks(blocks_o)
			call get_butterfly_minmaxrank(block_rand(1))
			rank_new_max = block_rand(1)%rankmax
			call copy_delete_randomizedbutterfly(block_rand(1),blocks_o,Memory) 
			deallocate(block_rand)

			return			
		end if		
	end do
	write(*,*)'randomized scheme not converged in '//TRIM(strings)//'. level: ',blocks_o%level_butterfly,error_inout,rank_new_max
	stop

    return

end subroutine Butterfly_randomized





subroutine Reconstruction_LL(block_rand,blocks_o,operand,blackbox_MVP_dat,operand1,option,stats)
    
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
	type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth
	real*8:: error_inout
	integer,allocatable::perms(:)
	type(partitionedblocks)::partitioned_block
	class(*):: operand	
	class(*),optional:: operand1
	
	type(matrixblock)::blocks_o,block_rand
	
	type(RandomBlock),allocatable :: vec_rand(:)
	type(Hoption)::option
	type(Hstat)::stats
	procedure(BF_MVP_blk)::blackbox_MVP_dat
	
	level_butterfly=block_rand%level_butterfly
    num_blocks=2**level_butterfly
	dimension_rank =block_rand%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
    allocate (vec_rand(1))
    
    allocate (vec_rand(1)%RandomVectorLL(0:level_butterfly+2))    
	Nbind = 1
	
	num_vect_sub = num_vect_subsub*Nbind
	
    ! random=>vec_rand(1)
	call Init_RandVect_Empty('T',vec_rand(1),num_vect_sub,block_rand,stats)	
	
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
			call Get_Randomized_Vectors_LL(block_rand,vec_rand(1),blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth,operand1)
			n2 = OMP_get_wtime()
			stats%Time_random(2) = stats%Time_random(2) + n2-n1	
			! Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,block_rand,vec_rand(1),option,stats)
			n2 = OMP_get_wtime()
			stats%Time_random(3) = stats%Time_random(3) + n2-n1		
		end do
	end do
	
	! deallocate(perms)
	
	
	! random=>vec_rand(1)
	call Delete_RandVect('T',vec_rand(1),level_butterfly)
	deallocate(vec_rand)

    return
    
end subroutine Reconstruction_LL



subroutine Reconstruction_RR(block_rand,blocks_o,operand,blackbox_MVP_dat,operand1,option,stats)
    
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
	type(Hoption)::option
	type(Hstat)::stats
	
    ! type(matricesblock), pointer :: blocks
    ! type(RandomBlock), pointer :: random
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
	
	type(matrixblock)::blocks_o,block_rand
	class(*):: operand
	class(*),optional::operand1
	type(RandomBlock),allocatable :: vec_rand(:)
	procedure(BF_MVP_blk)::blackbox_MVP_dat
	
	
	level_butterfly=block_rand%level_butterfly
		
    num_blocks=2**level_butterfly
	dimension_rank =block_rand%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
	! ! call assert(num_vectors==2**level_butterfly*dimension_rank,'incorrect num_vectors')
	
	! call assert(num_vectors==Nsub*dimension_rank,'incorrect num_vectors') !  check here later
	
    ! ! allocate (Random_Block(1))   !  check here later 
    
	allocate (vec_rand(1))
    allocate (vec_rand(1)%RandomVectorRR(0:level_butterfly+2))    
	
	Nbind = 1
	
	num_vect_sub = num_vect_subsub*Nbind
	
    ! random=>Random_Block(1)
	call Init_RandVect_Empty('N',vec_rand(1),num_vect_sub,block_rand,stats)

    level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    ! level_left_start = 0
	
	! call Zero_Butterfly(level_left_start,level_butterfly+1)

	do unique_nth=level_butterfly+1,level_left_start,-1
		if(mod(level_butterfly,2)==0)then
			Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))    !  check here later
		else 
			Nsub = NINT(2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))
		end if	
		Ng = 2**level_butterfly/Nsub	
		
		do ii = 1,Nsub/Nbind
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			n1 = OMP_get_wtime()
			call Get_Randomized_Vectors_RR(block_rand,vec_rand(1),blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth,operand1)
			n2 = OMP_get_wtime()
			stats%Time_random(2) = stats%Time_random(2) + n2-n1	
			! Time_Vector_inverse = Time_Vector_inverse + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,block_rand,vec_rand(1),option,stats)
			n2 = OMP_get_wtime()
			stats%Time_random(3) = stats%Time_random(3) + n2-n1		
		end do
	end do

	! random=>Random_Block(1)
	call Delete_RandVect('N',vec_rand(1),level_butterfly)
	
	deallocate(vec_rand)
	
    return
    
end subroutine Reconstruction_RR




subroutine Test_Reconstruction_Error(block_rand,block_o,operand,blackbox_MVP_dat,error,operand1)

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
	complex(kind=8),allocatable::Vdref(:,:),Id(:,:),Vd(:,:)

	
	
	type(matrixblock)::block_o,block_rand
	class(*)::operand
	class(*),optional::operand1
	procedure(BF_MVP_blk)::blackbox_MVP_dat
	
	
	level_butterfly=block_rand%level_butterfly
	num_blocks=2**level_butterfly
	
	num_vect = 1

	mm=0
	nn=0
	do i=1, num_blocks
		mm= mm + size(block_rand%butterflyU(i)%matrix,1)
		nn= nn + size(block_rand%butterflyV(i)%matrix,1)
	enddo
	
	
	allocate (Vdref(mm,num_vect))		
	Vdref=0
	allocate (Id(nn,num_vect))		
	Id=0
	allocate (Vd(mm,num_vect))		
	Vd=0
	
	call RandomMat(nn,num_vect,min(nn,num_vect),Id,0)

	call blackbox_MVP_dat(operand,block_o,'N',mm,nn,num_vect,Id,Vdref,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)

	call butterfly_block_MVP_randomized_dat(block_rand,'N',mm,nn,num_vect,Id,Vd,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))

	error = fnorm(Vd-Vdref,mm,num_vect)/fnorm(Vdref,mm,num_vect)
	
	
	deallocate(Vdref)
	deallocate(Vd)
	deallocate(Id)
	
    return                

end subroutine Test_Reconstruction_Error




subroutine Get_Randomized_Vectors_LL(block_rand,vec_rand,blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth,operand1)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk
    integer mm,nn,mm1,nn1,mn,level_butterfly
	character chara
	integer*8 idx_start   
   
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_right_start
	! type(RandomBlock), pointer :: random
	type(RandomBlock) :: vec_rand
	integer Nsub,Ng
	integer,allocatable:: mm_end(:),nn_end(:)
	
	class(*):: operand
	class(*),optional::operand1		
	type(matrixblock)::blocks_o,block_rand	
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	
	procedure(BF_MVP_blk)::blackbox_MVP_dat	
	
	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	
    level_butterfly=block_rand%level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

	groupm_start=blocks_o%row_group*2**level_butterfly
	groupn_start=blocks_o%col_group*2**level_butterfly	

	allocate(mm_end(0:num_blocks))
	mm_end=0
	allocate(nn_end(0:num_blocks))
	nn_end=0
    mm1=0
	nn1=0	
	do i=1, num_blocks
		if(allocated(blocks_o%ButterflyU) .and. size(blocks_o%ButterflyU)==num_blocks)then
			mm1 = size(blocks_o%butterflyU(i)%matrix,1)
			nn1=size(blocks_o%ButterflyV(i)%matrix,1)
		else 	
			mm1=basis_group(groupm_start+i-1)%tail-basis_group(groupm_start+i-1)%head+1
			nn1=basis_group(groupn_start+i-1)%tail-basis_group(groupn_start+i-1)%head+1
		endif	
		mm_end(i)=mm_end(i-1)+mm1
		nn_end(i)=nn_end(i-1)+nn1		
	end do
	mm = mm_end(num_blocks)
	nn = nn_end(num_blocks)	
	
	
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
	Ng = 2**level_butterfly/Nsub
	dimension_rank =block_rand%dimension_rank 


	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub)) ! #3 is used for output 
	
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	! groupm_start=groupm*2**(level_butterfly)
	! header_mm=basis_group(groupm_start)%head
	! idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i)	
		do i=(nth-1)*Ng+1, nth*Ng
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				! header_m=basis_group(groupm_start+i-1)%head
				! tailer_m=basis_group(groupm_start+i-1)%tail
				! mm=tailer_m-header_m+1
				! k=header_m-header_mm	

				! allocate(matrixtemp1(num_vect_subsub,mm))
				call RandomMat(mm_end(i)-mm_end(i-1),num_vect_subsub,min(mm_end(i)-mm_end(i-1),num_vect_subsub),RandomVectors_InOutput(1)%vector(mm_end(i-1)+1:mm_end(i),(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
				
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the left multiplied vectors
	
	call blackbox_MVP_dat(operand,blocks_o,'T',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)		


	! write(*,*)'aha',fnorm(RandomVectors_InOutput(1)%vector,mm,num_vect_sub),fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect_sub)
	
	
	k=0
	! random=>random_Block(1)
	do i=1, num_blocks
		mm=size(block_rand%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				vec_rand%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+mm
	enddo 
	
	k=0
	do i=1, num_blocks
		nn=size(block_rand%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				vec_rand%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
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
	
	deallocate(mm_end)
	deallocate(nn_end)
	
    return                

end subroutine Get_Randomized_Vectors_LL




subroutine Get_Randomized_Vectors_RR(block_rand,vec_rand,blocks_o,operand,blackbox_MVP_dat,nth_s,nth_e,num_vect_sub,unique_nth,operand1)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm1,nn1,mm,nn,mn,level_butterfly
    character chara
    
	integer*8 idx_start   
    integer dimension_rank
    ! integer header_mm, header_nn
	! integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,unique_nth,level_left_start
	! type(RandomBlock), pointer :: random
	type(RandomBlock) :: vec_rand
	integer Nsub,Ng,groupm_start,groupn_start
	integer,allocatable::mm_end(:),nn_end(:)
	
	class(*):: operand
	class(*),optional:: operand1
	
	type(matrixblock)::blocks_o,block_rand	
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	
	procedure(BF_MVP_blk)::blackbox_MVP_dat
	
	
    level_butterfly=block_rand%level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

	groupm_start=blocks_o%row_group*2**level_butterfly
	groupn_start=blocks_o%col_group*2**level_butterfly	

	allocate(mm_end(0:num_blocks))
	mm_end=0
	allocate(nn_end(0:num_blocks))
	nn_end=0
    mm1=0
	nn1=0	
	do i=1, num_blocks
		if(allocated(blocks_o%ButterflyU) .and. size(blocks_o%ButterflyU)==num_blocks)then
			mm1 = size(blocks_o%butterflyU(i)%matrix,1)
			nn1=size(blocks_o%ButterflyV(i)%matrix,1)
		else 	
			mm1=basis_group(groupm_start+i-1)%tail-basis_group(groupm_start+i-1)%head+1
			nn1=basis_group(groupn_start+i-1)%tail-basis_group(groupn_start+i-1)%head+1
		endif	
		mm_end(i)=mm_end(i-1)+mm1
		nn_end(i)=nn_end(i-1)+nn1		
	end do	
	mm = mm_end(num_blocks)
	nn = nn_end(num_blocks)	
	
	if(mod(level_butterfly,2)==0)then
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))    !  check here later
	else 
		Nsub = NINT(2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))
	end if	
	Ng = 2**level_butterfly/Nsub	
	dimension_rank =block_rand%dimension_rank 
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	

	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub)) ! #3 is used as output 
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	! groupn_start=groupn*2**(block_rand%level_butterfly)
	! header_nn=basis_group(groupn_start)%head
	! idx_start = 1
	
	do nth= nth_s,nth_e
		!$omp parallel do default(shared) private(i)	
		do i=(nth-1)*Ng+1, nth*Ng
		! do i=1, num_blocks
			! if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				! header_n=basis_group(groupn_start+i-1)%head
				! tailer_n=basis_group(groupn_start+i-1)%tail
				! nn=tailer_n-header_n+1
				! k=header_n-header_nn
				
				call RandomMat(nn_end(i)-nn_end(i-1),num_vect_subsub,min(nn_end(i)-nn_end(i-1),num_vect_subsub),RandomVectors_InOutput(1)%vector(nn_end(i-1)+1:nn_end(i),(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
	! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	call blackbox_MVP_dat(operand,blocks_o,'N',mm,nn,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)		

	
	k=0
	! random=>random_Block(1)
	do i=1, num_blocks
		nn=size(block_rand%butterflyV(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, nn
			do jj=1, num_vect_sub
				vec_rand%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
			enddo
		enddo
		! !$omp end parallel do
		k=k+nn
	enddo 

	k=0
	do i=1, num_blocks
		mm=size(block_rand%butterflyU(i)%matrix,1)
		! !$omp parallel do default(shared) private(ii,jj)
		do ii=1, mm
			do jj=1, num_vect_sub
				vec_rand%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
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
	
	deallocate(mm_end)
	deallocate(nn_end)
	
    return                

end subroutine Get_Randomized_Vectors_RR




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
subroutine butterfly_block_MVP_inverse_ABCD_dat(partitioned_block,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub, mv,nv
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vout_tmp(:,:),Vbuff(:,:)
   complex(kind=8) :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
   class(*),optional::operand1
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
		
		if(mm+nn/=N)write(*,*)'d3d',mm,nn,N
		call assert(mm+nn==N,'mm+nn/=N')  
		
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout
		
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

	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)
	   
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
subroutine butterfly_block_MVP_inverse_A_minusBDinvC_dat(partitioned_block,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: V_tmp1(:,:),V_tmp2(:,:),Vin_tmp(:,:),Vout_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn
   class(*)::partitioned_block
   class(*),optional::operand1
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

		
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout		
		
		
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
	   
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	   

   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_inverse_A_minusBDinvC_dat



subroutine butterfly_block_MVP_inverse_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::block_off1,block_off2
   integer groupn,groupm,mm,nn
   class(*)::ho_bf1
   class(*),optional::operand1
   type(matrixblock)::block_o

   select TYPE(ho_bf1)
   
   type is (hobf)

		block_off1 => ho_bf1%levels(ho_bf1%ind_lv)%BP(ho_bf1%ind_bk*2-1)%LL(1)%matrices_block(1)	
		block_off2 => ho_bf1%levels(ho_bf1%ind_lv)%BP(ho_bf1%ind_bk*2)%LL(1)%matrices_block(1)			
		
		groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
		groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
		
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout		
		
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

	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)
	   
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_inverse_minusBC_dat




subroutine butterfly_block_MVP_schulz_dat(schulz_op,block_Xn,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vbuff1(:,:),Vout_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2,a,b
   type(matrixblock)::block_Xn
   integer groupn,groupm,mm,nn
   class(*)::schulz_op
   class(*),optional::operand1
   type(matrixblock)::block_o
   real*8::scale_new

   select TYPE(schulz_op)   
   type is (schulz_operand)
	   select TYPE(operand1)
	   type is (integer)
			
			groupn=block_Xn%col_group    ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
			groupm=block_Xn%row_group    ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
			call assert(mm==nn,'block nonsquare')
			
			if(schulz_op%order==2)then
				mv=size(Vout,1)
				nv=size(Vout,2)
				allocate(Vout_tmp(mv,nv))
				Vout_tmp = Vout		
				
				allocate(Vin_tmp(N,num_vect_sub))
				Vin_tmp = Vin
				Vout = 0  

				allocate(Vbuff(nn,num_vect_sub))
				Vbuff = 0			
				
				
				if(trans=='N')then			
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! XnR
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,operand1)
					
					! AXnR
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2)
					Vout = 	Vbuff+Vout	
					
					! (2-AXn)R
					Vbuff = 2*Vin-Vout
					
					! Xn(2-AXn)R
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,operand1)					

				else if(trans=='T')then
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! RXn
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2,operand1)				
					
					! RXnA
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vout,Vbuff,ctemp1,ctemp2)
					Vbuff=	Vout+Vbuff
					
					! RXnAXn
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vin,ctemp1,ctemp2,operand1)
					
					! 2RXn - RXnAXn
					Vout = 2*Vout - Vin
				end if			

				
				Vin = Vin_tmp

				scale_new=schulz_op%scale*(2-schulz_op%scale)
				
				Vout = Vout-Vin*scale_new ! Xn(2-AXn)-I
				
				deallocate(Vin_tmp)
				deallocate(Vbuff)  
			
			else if(schulz_op%order==3)then
				
				mv=size(Vout,1)
				nv=size(Vout,2)
				allocate(Vout_tmp(mv,nv))
				Vout_tmp = Vout		
				
				allocate(Vin_tmp(nn,num_vect_sub))
				Vin_tmp = Vin
				Vout = 0  

				allocate(Vbuff(nn,num_vect_sub))
				Vbuff = 0				
				allocate(Vbuff1(nn,num_vect_sub))
				Vbuff1 = 0				
			
			
			
				if(trans=='N')then			
					ctemp1=1.0d0 ; ctemp2=0.0d0
					
					! AXnR
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2,operand1)
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vbuff1,ctemp1,ctemp2)
					Vbuff1 = 	Vbuff+Vbuff1	
					
					
					! (AXn)^2R
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff1,Vbuff,ctemp1,ctemp2,operand1)
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2)
					Vout = 	Vout+Vbuff
					
					! (3-3AXn+(AXn)^2)R
					Vbuff1 = 3*Vin-3*Vbuff1+Vout
					
					! Xn(3-3AXn+(AXn)^2)R
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff1,Vout,ctemp1,ctemp2,operand1)					

				else if(trans=='T')then
					ctemp1=1.0d0 ; ctemp2=0.0d0
					! RXn
					Vout=Vin
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vout,Vin,ctemp1,ctemp2,operand1)				
					
					! RXn*AXn
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2)
					Vbuff=	Vin+Vbuff
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vbuff1,ctemp1,ctemp2,operand1)

					! RXn*(AXn)^2
					call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans,mm,nn,num_vect_sub,Vbuff1,Vbuff,ctemp1,ctemp2)
					Vbuff=	Vbuff1+Vbuff
					call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2,operand1)

					
					! RXn*(3-3AXn+(AXn)^2)
					Vout = 3*Vin - 3*Vbuff1+Vout
				end if				
			
			
				Vin = Vin_tmp

				scale_new=schulz_op%scale*(3 - 3*schulz_op%scale + schulz_op%scale**2d0)
				
				Vout = Vout-Vin*scale_new ! Xn(2-AXn)-I
				
				deallocate(Vin_tmp)
				deallocate(Vbuff)  			
				deallocate(Vbuff1)  			
						
			endif


			Vout = a*Vout + b*Vout_tmp	
			deallocate(Vout_tmp)			

	   end select
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
   
end subroutine butterfly_block_MVP_schulz_dat




subroutine butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, M, N, num_vect_sub,mv,nv
   character trans,trans_new
   real*8::eps,memory
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:),Vbuff(:,:),Vout_tmp(:,:),matrixtmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2,a,b
   type(matrixblock)::block_Xn
   integer groupn,groupm,mm,nn
   class(*)::schulz_op
   class(*),optional::operand1
   type(matrixblock)::block_o

   select TYPE(schulz_op)   
   type is (schulz_operand)
		select TYPE(operand1)
		type is (integer)
			eps=0.8d0
			ctemp1=1d0
			ctemp2=0d0
			if(operand1==1)then ! X0
				Vin=conjg(Vin)
				if(trans=='N')trans_new='T'
				if(trans=='T')trans_new='N'				
				call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans_new,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				Vout = Vout + Vin
				Vin=conjg(Vin)
				Vout=conjg(Vout)
				schulz_op%scale=(2d0-eps)/schulz_op%A2norm**2d0
				Vout = Vout*schulz_op%scale
				

				
				! if(trans=='N')trans_new='T'
				! if(trans=='T')trans_new='N'				
				! call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,trans_new,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				! Vout = Vin - Vout/schulz_op%A2norm	

				
				
				! call copy_butterfly(schulz_op%matrices_block,block_o,memory)
				! call Butterfly_inverse_IplusButter_woodbury(block_o,memory)
				
				! call butterfly_block_MVP_randomized_dat(block_o,trans,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				! Vout = Vout + Vin
	
				
			else ! Xn
				call butterfly_block_MVP_randomized_dat(block_Xn,trans,M,N,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
				Vout=	Vout+Vin*schulz_op%scale					
			endif
			
		end select	   
   end select	   
			
end subroutine butterfly_block_MVP_schulz_Xn_dat



subroutine butterfly_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,M,N,mv,nv
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character trans
    ! real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2,a,b
	
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vbuff(:,:),Vout_tmp(:,:)
	
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
    class(*),optional::operand1	
	type(matrixblock)::block_o
	

   select TYPE(ho_bf1)
   
   type is (hobf)
		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk
	
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	
			
		if(trans=='N')then

			level_butterfly=block_o%level_butterfly
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
		  call butterfly_block_MVP_randomized_dat(block_o,'N',mm,nn,num_vect_sub,Vin,Vbuff,ctemp1,ctemp2)
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1	
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
			allocate(vec_new(mm,num_vect_sub))

			do level = ho_bf1%Maxlevel+1,level_c+1,-1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==ho_bf1%Maxlevel+1)then
					n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						! write(*,*)level,ii
						groupm_diag = ho_bf1%levels(level)%BP(ii)%row_group ! Note: row_group and col_group interchanged here   

						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1

						call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
		
					end do
					! !$omp end parallel do
					
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						! write(*,*)level,ii
						groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   

						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						
						call OneL_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))

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

			level_butterfly=block_o%level_butterfly
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
			do level = level_c+1,ho_bf1%Maxlevel+1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==ho_bf1%Maxlevel+1)then
					n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf1%levels(level)%BP(ii)%row_group ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)

					end do
					! !$omp end parallel do
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
						call OneL_block_MVP_inverse_dat(ho_bf1,level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,Vbuff(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))					
					end do
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1	
				end if
				Vbuff = vec_new
			end do	
			
			Vbuff = vec_new
			
			deallocate(vec_new)	

			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			n1 = OMP_get_wtime()
			call butterfly_block_MVP_randomized_dat(block_o,'T',mm,nn,num_vect_sub,Vbuff,Vout,ctemp1,ctemp2)	
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1		

			deallocate(Vbuff)
	
		end if	
	
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	
	
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select
   
    return                

end subroutine butterfly_block_MVP_Sblock_dat




subroutine Bplus_block_MVP_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,a,b
	type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1	
	
	
	real*8::n2,n1 	
	
   select TYPE(bplus)
   
   type is (blockplus)	
		
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	
		Vout=0	
		
		ctemp1=1.0d0 ; ctemp2=0.0d0
		
		groupn=bplus%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		groupm=bplus%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 		
		
		if(trans=='N')then
			call Bplus_block_MVP_randomized_dat(bplus,'N',mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2)			
		else if(trans=='T')then
			call Bplus_block_MVP_randomized_dat(bplus,'T',mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2)				
		endif
	   
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	   
	   
    class default
		write(*,*)"unexpected type"
		stop
	   
   end select     
 
end subroutine Bplus_block_MVP_Exact_dat




subroutine Bplus_block_MVP_Outter_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8),allocatable :: Vout_tmp(:,:)
	complex(kind=8) :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv

	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1
	
	real*8::n2,n1 	
	
	ctemp3=-1.0d0 ; ctemp4=1.0d0
	

   select TYPE(bplus)
   
   type is (blockplus)
   
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	   
   
		call assert(present(operand1),'operand1 cannot be skipped')
		
		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_Exact_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)
			
			call Bplus_block_MVP_randomized_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,operand1%Lplus)					
		end select
	
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)
	   
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select

end subroutine Bplus_block_MVP_Outter_Exact_dat




subroutine Bplus_block_MVP_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,a,b
	type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	
	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1	
	
	
	real*8::n2,n1 	
	
   select TYPE(ho_bf1)
   
   type is (hobf)	
		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk
		
		bplus_off1 => ho_bf1%levels(level_c)%BP(2*rowblock-1)	
		bplus_off2 => ho_bf1%levels(level_c)%BP(2*rowblock)	
		
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	
			
		if(trans=='N')then
			groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
			groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			allocate(vec_new(nn,num_vect_sub))
			vec_new = 0

			! get the right multiplied vectors
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_randomized_dat(bplus_off2,'N',nn,mm,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
			call Bplus_block_MVP_randomized_dat(bplus_off1,'N',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
			Vout = -Vout
			deallocate(vec_new)

	   else if(trans=='T')then
			groupn=bplus_off1%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
			groupm=bplus_off1%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			
			allocate(vec_new(nn,num_vect_sub))
			vec_new=0

			! get the right multiplied vectors
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_randomized_dat(bplus_off1,'T',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)		
			call Bplus_block_MVP_randomized_dat(bplus_off2,'T',nn,mm,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
			Vout = -Vout
			deallocate(vec_new)
			
	   end if
	   
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	   
	   
    class default
		write(*,*)"unexpected type"
		stop
	   
   end select     
 
end subroutine Bplus_block_MVP_minusBC_dat




subroutine Bplus_block_MVP_Outter_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8),allocatable :: Vout_tmp(:,:)
	complex(kind=8) :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv

	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1
	
	real*8::n2,n1 	
	
	ctemp3=-1.0d0 ; ctemp4=1.0d0
	

   select TYPE(ho_bf1)
   
   type is (hobf)
   
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	   
   
		call assert(present(operand1),'operand1 cannot be skipped')
		
		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_minusBC_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)
			call Bplus_block_MVP_randomized_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,operand1%Lplus)					
		end select
	
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)
	   
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select

end subroutine Bplus_block_MVP_Outter_minusBC_dat


subroutine Bplus_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,Ctemp,a,b
	type(blockplus),pointer::bplus_o
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,M,N,mv,nv
	integer level_butterfly,groupm_diag
	! real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)
	
	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1
	
	
	real*8::n2,n1 	

	
	
   select TYPE(ho_bf1)
   
   type is (hobf)
   
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	
			
		level_c = ho_bf1%ind_lv
		rowblock = ho_bf1%ind_bk
		if(trans=='N')then
			bplus_o => ho_bf1%levels(level_c)%BP(rowblock)	
			
			groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	

			groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 

			! get the right multiplied vectors
			idx_start_glo = basis_group(groupm)%head
			ctemp1=1.0d0 ; ctemp2=0.0d0
			call Bplus_block_MVP_randomized_dat(bplus_o,'N',mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
			allocate(vec_new(mm,num_vect_sub))

			do level = ho_bf1%Maxlevel+1,level_c+1,-1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==ho_bf1%Maxlevel+1)then
					n1 = OMP_get_wtime()
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf1%levels(level)%BP(ii)%row_group ! Note: row_group and col_group interchanged here   

						
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						
						call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)		
					end do
					! !$omp end parallel do
					
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						! write(*,*)level,ii
						groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   

						
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						
						call Bplus_block_MVP_inverse_dat(ho_bf1, level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))

					end do		
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1			
				end if
				
				Vout = vec_new
			end do
			deallocate(vec_new)
   
	   
	   else if(trans=='T')then
			bplus_o => ho_bf1%levels(level_c)%BP(rowblock)  
			groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
			groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1   
	   
			ctemp1=1.0d0 ; ctemp2=0.0d0
			! get the left multiplied vectors 
			idx_start_glo = basis_group(groupm)%head		
			allocate(vec_old(mm,num_vect_sub))
			allocate(vec_new(mm,num_vect_sub))	
			
			vec_old = Vin
			do level = level_c+1,ho_bf1%Maxlevel+1
				N_diag = 2**(level-level_c-1)
				idx_start_diag = (rowblock-1)*N_diag+1
				vec_new = 0
				
				if(level==ho_bf1%Maxlevel+1)then
					n1 = OMP_get_wtime()
					! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf1%levels(level)%BP(ii)%row_group ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
						call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
						&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					end do
					! !$omp end parallel do
					n2 = OMP_get_wtime()
					! ! time_tmp = time_tmp + n2 - n1
				else 
					n1 = OMP_get_wtime()
					do ii = idx_start_diag,idx_start_diag+N_diag-1
						groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   				
						idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
						idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
						call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))					
					end do
					n2 = OMP_get_wtime()
					! time_tmp = time_tmp + n2 - n1	
				end if
				vec_old = vec_new
			end do	
			deallocate(vec_new)
			n1 = OMP_get_wtime()
			call Bplus_block_MVP_randomized_dat(bplus_o,'T',mm,nn,num_vect_sub,vec_old,Vout,ctemp1,ctemp2)	
			n2 = OMP_get_wtime()
			deallocate(vec_old)	
	   end if
	   
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	   
	   
    class default
		write(*,*)"unexpected type"
		stop
	   
   end select  
 
end subroutine Bplus_block_MVP_Sblock_dat



subroutine Bplus_block_MVP_Outter_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8),allocatable :: Vout_tmp(:,:)
	complex(kind=8) :: ctemp1,ctemp2,ctemp3,ctemp4,a,b
	integer M,N,mv,nv

	class(*):: ho_bf1
	type(matrixblock)::block_o
	class(*),optional::operand1
	
	real*8::n2,n1 	
	
	ctemp3=-1.0d0 ; ctemp4=1.0d0
	

   select TYPE(ho_bf1)
   
   type is (hobf)
   
		mv=size(Vout,1)
		nv=size(Vout,2)
		allocate(Vout_tmp(mv,nv))
		Vout_tmp = Vout	   
   
		call assert(present(operand1),'operand1 cannot be skipped')
		
		select TYPE(operand1)
		type is (blockplus)
			call Bplus_block_MVP_Sblock_dat(ho_bf1,block_o,trans,M,N,num_vect_sub,Vin,Vout,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8),operand1)
			call Bplus_block_MVP_randomized_dat(operand1,trans,M,N,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,operand1%Lplus)					
		end select
	
	   Vout = a*Vout + b*Vout_tmp	
	   deallocate(Vout_tmp)	
	
   class default
		write(*,*)"unexpected type"
		stop
	   
   end select

end subroutine Bplus_block_MVP_Outter_Sblock_dat




subroutine OneL_block_MVP_inverse_dat(ho_bf1,level,ii,trans,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock),pointer::block_schur,block_off1,block_off2
   integer groupn,groupm,mm,nn
	type(hobf)::ho_bf1
   
   ctemp1=1.0d0
   ctemp2=0.0d0
   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin
   
	block_off1 => ho_bf1%levels(level)%BP(ii*2-1)%LL(1)%matrices_block(1)	
	block_off2 => ho_bf1%levels(level)%BP(ii*2)%LL(1)%matrices_block(1)
	block_schur => ho_bf1%levels(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)	
	
	groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	if(mm+nn/=N)write(*,*)'d33d',mm,nn,N
	
	call assert(mm+nn==N,'mm+nn/=N')  
		

		
	if(trans=='N')then
		call butterfly_block_MVP_randomized_dat(block_off1,trans,mm,nn,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)- Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)
		
		! write(2111,*)abs(Vout)
		
		call butterfly_block_MVP_randomized_dat(block_schur,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)			
		Vin(1:mm,1:num_vect_sub) = Vout(1:mm,1:num_vect_sub) + Vin(1:mm,1:num_vect_sub)
		Vin(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub) 

		! write(2112,*)abs(Vin)			
		
		call butterfly_block_MVP_randomized_dat(block_off2,trans,nn,mm,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)			
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
		
		! write(2113,*)abs(Vout)
		! stop
		
	else if(trans=='T')then
		call butterfly_block_MVP_randomized_dat(block_off2,trans,nn,mm,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub) - Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) 
		
		call butterfly_block_MVP_randomized_dat(block_schur,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)				
		Vin(1:mm,1:num_vect_sub) = Vout(1:mm,1:num_vect_sub) + Vin(1:mm,1:num_vect_sub)
		Vin(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub) 
		
		call butterfly_block_MVP_randomized_dat(block_off1,trans,mm,nn,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
		
	end if



   Vin = Vin_tmp
   deallocate(Vin_tmp)
   
end subroutine OneL_block_MVP_inverse_dat




subroutine Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,trans,N,num_vect_sub,Vin,Vout)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer level, ii, N, num_vect_sub
   character trans
   complex(kind=8) :: Vin(:,:), Vout(:,:)
   complex(kind=8),allocatable :: Vin_tmp(:,:)
   complex(kind=8) :: ctemp1,ctemp2
   type(matrixblock),pointer::block_o,block_off1,block_off2
   type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
   integer groupn,groupm,mm,nn
   type(hobf)::ho_bf1

   ctemp1=1.0d0
   ctemp2=0.0d0
   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin
   
	bplus_off1 => ho_bf1%levels(level)%BP(ii*2-1)	
	bplus_off2 => ho_bf1%levels(level)%BP(ii*2)
	
	groupn=bplus_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=bplus_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	if(mm+nn/=N)write(*,*)'dd',mm,nn,N
	
	call assert(mm+nn==N,'mm+nn/=N')  
		
  
	bplus_o => ho_bf1%levels(level)%BP_inverse_schur(ii)
	if(trans=='N')then
		call Bplus_block_MVP_randomized_dat(bplus_off1,trans,mm,nn,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)- Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)
		
		! write(2111,*)abs(Vout)
		
		
		
		call Bplus_block_MVP_randomized_dat(bplus_o,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)			
		Vin(1:mm,1:num_vect_sub) = Vout(1:mm,1:num_vect_sub) + Vin(1:mm,1:num_vect_sub)
		Vin(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub) 

		! write(2112,*)abs(Vin)			
		
		call Bplus_block_MVP_randomized_dat(bplus_off2,trans,nn,mm,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)			
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
		
		! write(2113,*)abs(Vout)
		! stop
		
	else if(trans=='T')then
	! write(*,*)'good1'
		call Bplus_block_MVP_randomized_dat(bplus_off2,trans,nn,mm,num_vect_sub,&
		&Vin(1+mm:N,1:num_vect_sub),Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub) - Vout(1:mm,1:num_vect_sub)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) 
	! write(*,*)'good2'	
		call Bplus_block_MVP_randomized_dat(bplus_o,trans,mm,mm,num_vect_sub,&
		&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)				
		Vin(1:mm,1:num_vect_sub) = Vout(1:mm,1:num_vect_sub) + Vin(1:mm,1:num_vect_sub)
		Vin(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub) 
	! write(*,*)'good3'	
		call Bplus_block_MVP_randomized_dat(bplus_off1,trans,mm,nn,num_vect_sub,&
		&Vin(1:mm,1:num_vect_sub),Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
		Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
	! write(*,*)'good4'	
	end if
  


   Vin = Vin_tmp
   deallocate(Vin_tmp)
   
end subroutine Bplus_block_MVP_inverse_dat



subroutine Bplus_block_MVP_BplusB_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1	
	
	
	real*8::n2,n1 	
	
   select TYPE(bplus)
   
   type is (blockplus)	
   
	 select TYPE(operand1)
	   type is (matrixblock)	
		
			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout	
		
			! groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
			! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			
			level_butterfly = block_o%level_butterfly
			num_blocks=2**level_butterfly			
			mm=0
			nn=0	
			do i=1, num_blocks
				mm = mm+size(block_o%butterflyU(i)%matrix,1)
				nn = nn+size(block_o%ButterflyV(i)%matrix,1)
			enddo
		
			
			
			allocate(vec_new(mm,num_vect_sub))
			vec_new = 0		
		
			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call butterfly_block_MVP_randomized_dat(operand1,'N',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
				
				call Bplus_block_MVP_randomized_dat(bplus,'N',mm,mm,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
				Vout = Vout + vec_new
				
				deallocate(vec_new)

		   else if(trans=='T')then

				! get the left multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call Bplus_block_MVP_randomized_dat(bplus,'T',mm,mm,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
				vec_new = vec_new + Vin
				
				call butterfly_block_MVP_randomized_dat(operand1,'T',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
				deallocate(vec_new)
				
		   end if
	   
		   Vout = a*Vout + b*Vout_tmp	
		   deallocate(Vout_tmp)	   
	   
		class default
			write(*,*)"unexpected type"
			stop
		   
	   end select 	   
	   
    class default
		write(*,*)"unexpected type"
		stop
	   
   end select     
 
end subroutine Bplus_block_MVP_BplusB_dat



subroutine Bplus_block_MVP_BBplus_dat(bplus,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub,M,N,mv,nv
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,a,b
	! type(blockplus),pointer::bplus_o,bplus_off1,bplus_off2
	integer groupn,groupm,mm,nn

	integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
	integer level_butterfly,groupm_diag
	! real*8 a,b,c,d
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),Vout_tmp(:,:)

	integer*8 idx_start   
	integer level_blocks
	integer groupm_start, groupn_start,dimension_rank
	integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	type(vectorsblock), pointer :: random1, random2
	
	class(*):: bplus
	type(matrixblock)::block_o
	class(*),optional::operand1	
	
	
	real*8::n2,n1 	
	
   select TYPE(bplus)
   
   type is (blockplus)	
   
	 select TYPE(operand1)
	   type is (matrixblock)	
		
			mv=size(Vout,1)
			nv=size(Vout,2)
			allocate(Vout_tmp(mv,nv))
			Vout_tmp = Vout			
		
			! groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
			! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 

			level_butterfly = block_o%level_butterfly
			num_blocks=2**level_butterfly			
			mm=0
			nn=0	
			do i=1, num_blocks
				mm = mm+size(block_o%butterflyU(i)%matrix,1)
				nn = nn+size(block_o%ButterflyV(i)%matrix,1)
			enddo
				
			
			allocate(vec_new(nn,num_vect_sub))
			vec_new = 0		
		
			if(trans=='N')then
				! get the right multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				
				call Bplus_block_MVP_randomized_dat(bplus,'N',nn,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
				vec_new = vec_new + Vin
				
				call butterfly_block_MVP_randomized_dat(operand1,'N',mm,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
				deallocate(vec_new)

		   else if(trans=='T')then

				! get the left multiplied vectors
				ctemp1=1.0d0 ; ctemp2=0.0d0
				call butterfly_block_MVP_randomized_dat(operand1,'T',mm,nn,num_vect_sub,Vin,vec_new,ctemp1,ctemp2)
				
				call Bplus_block_MVP_randomized_dat(bplus,'T',nn,nn,num_vect_sub,vec_new,Vout,ctemp1,ctemp2)
				Vout = Vout + vec_new
				
				deallocate(vec_new)
				
		   end if
		   
		   Vout = a*Vout + b*Vout_tmp	
		   deallocate(Vout_tmp)		   
	   
		class default
			write(*,*)"unexpected type"
			stop
		   
	   end select 	   
	   
    class default
		write(*,*)"unexpected type"
		stop
	   
   end select     
 
end subroutine Bplus_block_MVP_BBplus_dat


	

subroutine MultiLrandomized_Onesubblock(rank0,rankrate,rankthusfar,blocks,operand,blackbox_MVP_dat,error_inout,strings,option,stats,operand1)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   use Butterfly_compress_forward
   implicit none

    type(blockplus),pointer::bplus
	integer:: ii,ll,bb,jj,bb_o,tt,rank0,rankthusfar
    real*8 Memory,rtemp,error_inout,n2,n1,mem_vec,rankrate	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vout3(:,:),Vin(:,:)
	integer M,N,idx_start_n,idx_start_m,idx_start_n_loc,idx_end_n_loc,idx_start_m_loc,idx_end_m_loc,mm,nn,rmax,rank,idx_start_n_ref,idx_start_m_ref,idx_end_n_ref,idx_end_m_ref
	complex(kind=8)::ctemp1,ctemp2,Ctemp
	type(matrixblock)::blocks
	type(matrixblock)::block_dummy
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:)
	real*8, allocatable :: Singular(:)
	integer level_c,rowblock,Nactive
	integer,allocatable::boxindex(:)
	integer Chunksize, Nchunk, Nidx, idx_s,cc
	class(*):: operand
	class(*),optional:: operand1
	character(*)  :: strings
	type(Hoption)::option
	type(Hstat)::stats
	complex(kind=8),allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	complex(kind=8), allocatable :: matU_glo(:,:), matV_glo(:,:)
	procedure(BF_MVP_blk)::blackbox_MVP_dat

	select TYPE(operand1)
	type is (blockplus)
		! select TYPE(operand)
		! type is (hobf)	
			
			! level_c = operand%ind_lv
			! rowblock = operand%ind_bk
		
			N = basis_group(operand1%col_group)%tail - basis_group(operand1%col_group)%head + 1	
			M = basis_group(operand1%row_group)%tail - basis_group(operand1%row_group)%head + 1	
			
			
			ctemp1 = 1.0d0
			ctemp2 = 0.0d0
			
			idx_start_n = basis_group(operand1%col_group)%head
			idx_start_m = basis_group(operand1%row_group)%head
				
			! blocks => operand1%LL(2)%matrices_block(bb_o)
			
			idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
			idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
			idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
			idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
			
			mm = idx_end_m_loc - idx_start_m_loc + 1
			nn = idx_end_n_loc - idx_start_n_loc + 1
			
			
			do tt=1,10
			
				rmax = ceiling_safe(max(rank0*2,rankthusfar)*rankrate**(tt-1)) !+ (tt-1)*10  !!!!! be careful here
				
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
					stats%Mem_int_vec = max(stats%Mem_int_vec,mem_vec)			
					
					! do ii=idx_start_n_loc,idx_end_n_loc
					! do jj=1,rmax
						! RandVectInR(ii,jj)=random_complex_number()
					! end do
					! end do
					
					call RandomMat(nn,Nidx,min(nn,Nidx),RandVectInR(idx_start_n_loc:idx_end_n_loc,1:Nidx),0)	


					call blackbox_MVP_dat(operand,block_dummy,'N',M,N,Nidx,RandVectInR,RandVectOutR,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
					
					
					
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
					stats%Mem_int_vec = max(stats%Mem_int_vec,mem_vec)			
					! write(*,*)mem_vec,'dd2',M,N,Nidx
					! do ii=idx_start_m_loc,idx_end_m_loc
					! do jj=1,rmax
						! RandVectInL(ii,jj)=random_complex_number()
					! end do
					! end do
					
					call RandomMat(mm,Nidx,min(mm,Nidx),RandVectInL(idx_start_m_loc:idx_end_m_loc,1:Nidx),0)		
							
		
					call blackbox_MVP_dat(operand,block_dummy,'T',M,N,Nidx,RandVectInL,RandVectOutL,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
					! write(*,*)'yani 3'	
					matRrow(1:mm,idx_s:idx_s+Nidx-1) = RandVectInL(idx_start_m_loc:idx_start_m_loc+mm-1,1:Nidx)
					deallocate(RandVectInL)
					matZcRrow(1:nn,idx_s:idx_s+Nidx-1) = RandVectOutL(idx_start_n_loc:idx_start_n_loc+nn-1,1:Nidx)	
					deallocate(RandVectOutL)
				end do
																	
				matRrow = conjg(matRrow)															   
				matZcRrow = conjg(matZcRrow)
							

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
				
				n1 = OMP_get_wtime()
				allocate(matU_glo(mm,rmax))
				allocate(matV_glo(rmax,nn))
				allocate(Singular(rmax))

				
				! write(*,*)mm,nn,rmax,'didi'
				
				call RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU_glo,matV_glo,Singular,mm,nn,rmax,rank,option%tol_LS,option%tol_SVD)				
				rankthusfar = max(rankthusfar,rank)
				! write(*,*)'yani 4'
				do ii=1,rank
					matV_glo(ii,:) = matV_glo(ii,:) * Singular(ii)
				end do
				
				deallocate(matRcol,matZRcol,matRrow,matZcRrow,Singular)
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1


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

				call blackbox_MVP_dat(operand,block_dummy,'N',M,N,1,RandVectInR,RandVectOutR,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
				
				Vout1 = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:1)
				deallocate(RandVectInR)
				deallocate(RandVectOutR)
				! write(*,*)'yani 5'

				call gemm_omp(matV_glo(1:rank,1:nn),Vin,Vout3,rank,nn,1)
				call gemm_omp(matU_glo(1:mm,1:rank),Vout3,Vout2,mm,rank,1)
				
				
				error_inout = fnorm(Vout2-Vout1,mm,1)/fnorm(Vout1,mm,1)
				! write(*,*)error_inout,bb_o,'ninini'				
				
				deallocate(Vin)
				deallocate(Vout1)
				deallocate(Vout2)
				deallocate(Vout3)
				
				! write(*,*)tt,rmax,rank,'Sblock',error_inout
				
				if(error_inout>option%tol_rand*1d-1)then
					if(min(mm,nn)==rmax)then
						write(*,*)tt,rmax,strings,error_inout,rank
						write(*,*)'no need to increase rmax, try increase RandomizedSVD tolerance'
						stop
					end if		
					deallocate(matU_glo)
					deallocate(matV_glo)
				else
#if PRNTlevel >= 2				
					write(*,'(A33,A8,I4,A8,I2,A7,Es14.7)')' Onesub ',' rank:',rank,'Ntrial:',tt,' error:',error_inout
#endif					

					! write(*,*)bb_o,error_inout,rmax,rank,'good'
					exit
				end if
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
			end do

			if(error_inout>option%tol_rand*1d-1)then
				write(*,*)'matSub_glo no correct, needs more work'
				stop
			end if

			
			!!!! construct butterfly at all subsequent levels including level 2 (this part is generic and doesn't require accuracy test)
			idx_start_n_ref = idx_start_n_loc + idx_start_n - 1
			idx_end_n_ref = idx_end_n_loc + idx_start_n - 1
			idx_start_m_ref = idx_start_m_loc + idx_start_m - 1
			idx_end_m_ref = idx_end_m_loc + idx_start_m - 1	
			
			if(option%TwoLayerOnly==1)then
				ll=2
				Nactive = 0
				do bb = 1,operand1%LL(ll)%Nbound
					if(basis_group(operand1%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(operand1%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
						Nactive = Nactive + 1
						! boxindex(Nactive) = bb
					end if
				end do
				call assert(Nactive==1,'Nactive should be one')
				

				! bb = boxindex(1)
				! blocks => Bplus_randomized(1)%LL(ll)%matrices_block(bb)
				
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

				
																				 
				operand1%LL(ll)%rankmax = max(operand1%LL(ll)%rankmax,blocks%rankmax)
						
			else 
				write(*,*)"twolayer only"
				stop
			end if

			return		
		! end select		
	class default
		write(*,*)"unexpected type"
		stop					
	end select	

end subroutine MultiLrandomized_Onesubblock



subroutine Buplus_randomized(level_butterfly,bplus_o,operand,rank0_inner,rankrate_inner,blackbox_MVP_dat_inner,rank0_outter,rankrate_outter,blackbox_MVP_dat_outter,error_inout,strings,option,stats)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    ! type(blockplus),pointer::bplus
    type(blockplus)::bplus_o
	integer:: ii,ll,bb
    real*8 rtemp,error,Memory,n2,n1,rate,error_inout,err_avr,rankrate_inner,rankrate_outter	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall,M,N,err_cnt
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	complex(kind=8) ctemp, ctemp1, ctemp2	
	integer level_c,rowblock,rank_new_max,rank0_inner,rank0_outter
	type(matrixblock),pointer::block_off1,block_off2,block_o
	class(*)::operand
	character(*)  :: strings
	integer rankthusfar
	type(Hoption)::option
	type(Hstat)::stats
	type(blockplus) :: Bplus_randomized
	procedure(BF_MVP_blk)::blackbox_MVP_dat_inner,blackbox_MVP_dat_outter
	
	error_inout=0
	! Memory = 0	
	call assert(bplus_o%Lplus>=2,'this is not a multi Bplus in Buplus_randomized')
	
	call Initialize_Bplus_FromInput(bplus_o,Bplus_randomized)
	
	! write(*,*)'hhhhh1',Bplus_randomized%LL(1)%matrices_block(1)%level
	n1 = OMP_get_wtime()

	rankthusfar = 0	
	do bb =1,Bplus_randomized%LL(2)%Nbound
		! ! write(*,*)bb,Bplus_randomized%LL(2)%Nbound,'dddd'

		
		! rank0_inner = ho_bf1%levels(level_c)%BP(2*rowblock-1)%LL(2)%rankmax
		! rankrate_inner = 2.0d0
		block_o => Bplus_randomized%LL(2)%matrices_block(bb)
		! ho_bf1%ind_lv = level_c
		! ho_bf1%ind_bk = rowblock
		
		call MultiLrandomized_Onesubblock(rank0_inner,rankrate_inner,rankthusfar,block_o,operand,blackbox_MVP_dat_inner,error,strings,option,stats,Bplus_randomized)		
		error_inout = max(error_inout, error)
		! write(*,*)'go'
	end do
	! call Test_Error_RR_Inner_Exact(bplus_o)
	n2 = OMP_get_wtime()
	stats%Time_random(4) = stats%Time_random(4) + n2-n1	
	
	
	! level_c = ho_bf1%ind_lv
	! rowblock = ho_bf1%ind_bk
	
	! block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	! ! block_off1 => ho_bf1%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	! block_off2 => ho_bf1%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	
	

	! level_butterfly=int((maxlevel_for_blocks-bplus_o%level)/2)*2
	! rank0_outter = max(block_off1%rankmax,block_off2%rankmax)
	! rate_outter=1.2d0
	call Butterfly_randomized(level_butterfly,rank0_outter,rankrate_outter,Bplus_randomized%LL(1)%matrices_block(1),operand,blackbox_MVP_dat_outter,error,'Outter',option,stats,Bplus_randomized)
	error_inout = max(error_inout, error)

	
	call delete_Bplus(bplus_o)
	call copy_delete_Bplus(Bplus_randomized,bplus_o,Memory)
	! deallocate(Bplus_randomized)
	

	rank_new_max = 0
	do ll=1,bplus_o%Lplus
		rank_new_max = max(rank_new_max,bplus_o%LL(ll)%rankmax)
	end do	
	
	
	
#if PRNTlevel >= 2
	write(*,'(A20,A8,I3,A8,I3,A11,Es14.7)')strings,' rank:',rank_new_max,' L_butt:',bplus_o%LL(1)%matrices_block(1)%level_butterfly,' error:',error_inout
#endif	
    return

end subroutine Buplus_randomized



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

    return

end subroutine Butterfly_inverse_IplusButter_woodbury



end module Randomized_reconstruction