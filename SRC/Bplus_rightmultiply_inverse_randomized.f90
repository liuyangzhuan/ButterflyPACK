module Bplus_rightmultiply
use Butterfly_rightmultiply
integer rankthusfarS					
contains





subroutine Bplus_Sblock_randomized_memfree(level_c,rowblock,Memory)

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
    type(blockplus),pointer::bplus
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,rank,ntry
	real*8:: error_inout
	real*8:: n1,n2,Memory
	
	Memory = 0
	! write(*,*)'caca',level_c,rowblock
    bplus =>  cascading_factors(level_c)%BP(rowblock)
	! write(*,*)'caca1',level_c,rowblock
	if(bplus%Lplus==1)then
		call OneL_Sblock_randomized_memfree(level_c,rowblock,Memory)
	else 
	! write(*,*)'cacadd1',level_c,rowblock
		call MultiLSblock_randomized_memfree(level_c,rowblock,Memory)
	end if
	
    return

end subroutine Bplus_Sblock_randomized_memfree


subroutine Bplus_Sblock_randomized_symmetric(level_c,rowblock)

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
    type(blockplus),pointer::bplus
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,rank,ntry
	real*8:: error_inout
	real*8:: n1,n2
	
    bplus =>  cascading_factors(level_c)%BP(rowblock)
	
	if(bplus%Lplus==1)then
		call OneL_Sblock_randomized_symmetric(level_c,rowblock)
		call Butterfly_sym2asym(cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1))
	else 
		write(*,*)'not done yet Bplus_Sblock_randomized_symmetric'
		stop
	end if
	
    return

end subroutine Bplus_Sblock_randomized_symmetric




subroutine OneL_Sblock_randomized_memfree(level_c,rowblock,Memory)

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
	real*8::rank_new_avr,error,tmpfact
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,rank,ntry,idx_start_m_ref,idx_start_n_ref
	real*8:: error_inout
	real*8:: n1,n2,Memory
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:),matsub_tmp(:,:)
	type(blockplus),pointer::bplus
	type(matrixblock)::block_tmp					 
	
	Memory = 0
	
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	bplus =>  cascading_factors(level_c)%BP(rowblock)												 
	
	! ! ! call copy_butterfly(block_o,block_old)
	
	
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2

	if(level_butterfly==0)then
		if(level_c>=maxlevel_for_blocks-1)then
		
			call OneL_Sblock_LowRank(level_c,rowblock)
		
		
		
			! mm=size(block_o%ButterflyU(1)%matrix,1)
			! nn=size(block_o%ButterflyU(1)%matrix,2)
			! allocate(matrixtmp(mm,nn))
			! call gemmf90(cascading_factors(Maxlevel_for_blocks+1)%matrices_block_inverse(rowblock)%fullmat, block_o%ButterflyU(1)%matrix, matrixtmp,'N','N')
			! block_o%ButterflyU(1)%matrix = matrixtmp
			! deallocate(matrixtmp)


		else 		
			write(*,*)'unexpected level_c'
			stop
		end if		
	else 
	! if(level_c<=2)then
		! tmpfact = SVD_tolerance_factor
		! ! SVD_tolerance_factor = SVD_tolerance_factor* 1d-2
		
		! write(*,*)'explicit computing butterfly for Sblock '

		! mm = basis_group(bplus%row_group)%tail - basis_group(bplus%row_group)%head + 1
		! nn = basis_group(bplus%col_group)%tail - basis_group(bplus%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(nn,nn))
		! ! allocate(matout(mm,mm))
		! allocate(matsub_glo(mm,nn))
		! matin = 0
		! do ii=1,nn
			! matin(ii,ii) = 1
		! end do
! ! write(*,*)'ddd1'
		! ! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',nn,matin,matsub_glo)
		! call OneL_block_MVP_Sblock_dat(level_c,rowblock,'N',nn,matin,matsub_glo)
		! ! write(*,*)'ddd2'		
		! ! call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,nn,matin,matsub_glo,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)		

			
		! block_tmp%row_group = bplus%LL(1)%matrices_block(1)%row_group
		! block_tmp%col_group = bplus%LL(1)%matrices_block(1)%col_group
		! block_tmp%level = bplus%LL(1)%matrices_block(1)%level
		! block_tmp%style = bplus%LL(1)%matrices_block(1)%style
		
		! idx_start_m_ref = basis_group(bplus%row_group)%head
		! idx_start_n_ref = basis_group(bplus%col_group)%head
		! write(*,*)'ddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		! write(*,*)'wocaonimaSblockOutter',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! deallocate(matin)
		! deallocate(matsub_glo)
		! call delete_blocks(block_tmp)	



		! mm = basis_group(bplus%row_group)%tail - basis_group(bplus%row_group)%head + 1
		! nn = basis_group(bplus%col_group)%tail - basis_group(bplus%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(mm,mm))
		! ! allocate(matout(mm,mm))
		! allocate(matsub_glo(mm,nn))
		! allocate(matsub_tmp(nn,mm))
		! matin = 0
		! do ii=1,mm
			! matin(ii,ii) = 1
		! end do
! write(*,*)'dddd1'
		! ! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',mm,matin,matsub_tmp)
		! call OneL_block_MVP_Sblock_dat(level_c,rowblock,'T',mm,matin,matsub_tmp)
! write(*,*)'dddd2'		
		! ! call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'T',mm,nn,mm,matin,matsub_tmp,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)		
		! call copymatT_OMP(matsub_tmp,matsub_glo,nn,mm)
		! deallocate(matsub_tmp)
			
		! block_tmp%row_group = bplus%LL(1)%matrices_block(1)%row_group
		! block_tmp%col_group = bplus%LL(1)%matrices_block(1)%col_group
		! block_tmp%level = bplus%LL(1)%matrices_block(1)%level
		! block_tmp%style = bplus%LL(1)%matrices_block(1)%style
		
		! idx_start_m_ref = basis_group(bplus%row_group)%head
		! idx_start_n_ref = basis_group(bplus%col_group)%head
! write(*,*)'dddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		
		! write(*,*)'wocaonimaSblockOutter',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! ranktmp_glo = block_tmp%rankmax*1.2
		
		! deallocate(matin)
		! deallocate(matsub_glo)
		! call delete_blocks(block_tmp)	
		! write(*,*)'Done computing ' 
		! SVD_tolerance_factor = tmpfact
	! end if
		! ! T0=secnds(0.0)
	
		do tt =1,10
			do ntry=1,1
			n1 = OMP_get_wtime()
			call Initialize_Butterfly_Sblock(block_o,level_c,tt-1)
			n2 = OMP_get_wtime()
			Time_Init_forward = Time_Init_forward + n2 -n1 
			
	
			n1 = OMP_get_wtime()
			call OneL_Reconstruction_LL_Sblock(level_c,rowblock)	
			call OneL_Reconstruction_RR_Sblock(level_c,rowblock,error_inout)
			n2 = OMP_get_wtime()
			Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1
				
			call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
			! rank_new_max = butterfly_block_randomized(1)%rankmax
				
			! write(*,*)tt,error_inout,butterfly_block_randomized(1)%dimension_rank,butterfly_block_randomized(1)%rankmax,'OneL_Sb'
			
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
				write(*,'(A10,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7)')'OneL No. ',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
				return
			end if
			end do
		end do
		write(*,*)'randomized scheme not converged',error_inout
		stop
		
	end if
	
    return

end subroutine OneL_Sblock_randomized_memfree




subroutine OneL_Sblock_randomized_symmetric(level_c,rowblock)

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
	integer niter,rank
	real*8:: error_inout,rtmp
	real*8:: n1,n2
	
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	if(level_butterfly==0)then
		

		! mm=size(block_o%ButterflyU(1)%matrix,1)
		! nn=size(block_o%ButterflyU(1)%matrix,2)
		! allocate(matrixtmp(mm,nn))
		! call gemmf90(cascading_factors(Maxlevel_for_blocks+1)%matrices_block_inverse(rowblock)%fullmat, block_o%ButterflyU(1)%matrix, matrixtmp,'N','N')
		! block_o%ButterflyU(1)%matrix = matrixtmp
		! deallocate(matrixtmp)
		! rank_new_max = nn
		! rank_new_avr = nn
		
		call OneL_Sblock_LowRank(level_c,rowblock)
		



	else if(level_butterfly<=2)then
	! write(*,*)'here?'

		! call Butterfly_Sblock_randomized(level_c,rowblock)
		call OneL_Sblock_randomized_memfree(level_c,rowblock,rtmp)
		! write(*,*)'herddde?'


	else 
		! ! T0=secnds(0.0)
	
		do tt =1,1
			
			n1 = OMP_get_wtime()      
			call OneL_Get_Randomized_Vectors_Sblock_Memeff(level_c,rowblock,tt-1)
			n2 = OMP_get_wtime()
			Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_Symmetric(block_o%row_group,block_o%col_group)
			call OneL_Test_Error_RR_Sblock(level_c,rowblock,error_inout)
			n2 = OMP_get_wtime()
			Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1
			
			
			if(error_inout>iter_tolerance)then
				call Delete_randomized_butterfly()
			else 
				call delete_blocks(block_o)
				call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
				rank_new_max = butterfly_block_randomized(1)%rankmax				
				call copy_delete_randomizedbutterfly(butterfly_block_randomized(1),block_o)
				deallocate(butterfly_block_randomized)
				! call copy_randomizedbutterfly(butterfly_block_randomized(1),block_o)
				! call Delete_randomized_butterfly()

				write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7)')'     No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
				! stop


				return
			end if
			
		end do
		write(*,*)'randomized scheme not converged',error_inout
		stop
		
	end if
	
	


    return

end subroutine OneL_Sblock_randomized_symmetric




subroutine OneL_Get_Randomized_Vectors_Sblock_Memeff(level_c,rowblock,kover)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d,mem_vec
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),vec_old_loc(:,:),vec_new_loc(:,:),matrixtmp1(:,:),fullmat(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank,num_vectorsR,num_vectorsL,DimMax,DimMax_loc
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n,kover,levelm
	integer, allocatable :: ipiv(:)
	integer m,n,edge_m,edge_n
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	dimension_rank= block_o%rankmax*2 + kover !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	! dimension_rank= block_o%rankmax +3 + kover !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	! if(level_c==2)dimension_rank=9
	! if(level_c==1)dimension_rank=18
	allocate (butterfly_block_randomized(1))
    butterfly_block_randomized(1)%dimension_rank=dimension_rank
	butterfly_block_randomized(1)%level_butterfly = level_butterfly
	
	levelm = ceiling_safe(dble(level_butterfly)/2d0)
	num_blocks=2**level_butterfly

! first consider the right multiplied vectors
	Nsub = 2**(level_butterfly-levelm)
    Ng = 2**level_butterfly/Nsub
	num_vectorsR=Nsub*dimension_rank
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    groupn=block_o%col_group  
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	DimMax = max(mm,nn)

		! ! allocate(fullmat(mm,nn))
		! ! do m=1,mm
			! ! do n =1,nn
				! ! edge_m = basis_group(groupm)%head + m-1
				! ! edge_n = basis_group(groupn)%head + n-1
				! ! call element_Zmn(edge_m,edge_n,ctemp)
				! ! fullmat(m,n) = ctemp	
			! ! end do
		! ! end do		
	
		! ! allocate(matrixtmp1(mm,mm))
		! ! do m=1,mm
			! ! do n =1,mm
				! ! edge_m = basis_group(groupm)%head + m-1
				! ! edge_n = basis_group(groupm)%head + n-1
				! ! call element_Zmn(edge_m,edge_n,ctemp)
				! ! matrixtmp1(m,n) = ctemp	
			! ! end do
		! ! end do


		! ! allocate(ipiv(mm))
		! ! call getrff90(matrixtmp1,ipiv)
		! ! ! write(*,*)shape(matrixtemp1)
		! ! call getrif90(matrixtmp1,ipiv)	
		! ! deallocate(ipiv)	
	
	
	
	
	
	allocate(RandVectInR(DimMax,num_vectorsR))
	RandVectInR=0
	allocate(vec_old(DimMax,num_vectorsR))
	allocate(vec_new(DimMax,num_vectorsR))
	
	idx_start_glo = basis_group(groupm)%head	
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
				 RandVectInR(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo 	

	! get the right multiplied vectors
    vec_old = RandVectInR
	call Butterfly_Partial_MVP_dat_memeff(block_o,'N',0,level_butterfly+1,DimMax,num_vectorsR,vec_old,vec_new)
	vec_old = vec_new
	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1

			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectorsR,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vectorsR),vec_new(idx_start_loc:idx_end_loc,1:num_vectorsR),ctemp1,ctemp2)		
			else 
				DimMax_loc=basis_group(groupm_diag)%tail-basis_group(groupm_diag)%head+1 
				allocate(vec_old_loc(idx_end_loc-idx_start_loc+1,num_vectorsR))
				allocate(vec_new_loc(idx_end_loc-idx_start_loc+1,num_vectorsR))
				vec_old_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsR) = vec_old(idx_start_loc:idx_end_loc,1:num_vectorsR)
				vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsR) = vec_new(idx_start_loc:idx_end_loc,1:num_vectorsR)				
				! call Butterfly_Partial_MVP_dat_memeff(cascading_factors(level)%matrices_block_inverse(ii),'N',0,cascading_factors(level)%matrices_block_inverse(ii)%level_butterfly+1,&
				! &DimMax_loc,num_vectorsR,vec_old_loc,vec_new_loc)
				
				
				call OneL_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectorsR,vec_old_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsR),vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsR))
				! write(*,*)'R:really?',fnorm(vec_old_loc,idx_end_loc-idx_start_loc+1,num_vectorsR),fnorm(vec_new_loc,idx_end_loc-idx_start_loc+1,num_vectorsR),fnorm(vec_old_loc-vec_new_loc,idx_end_loc-idx_start_loc+1,num_vectorsR)
				
				vec_new(idx_start_loc:idx_end_loc,1:num_vectorsR) = vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsR)  !!!!+vec_old(idx_start_loc:idx_end_loc,1:num_vectorsR)
				deallocate(vec_old_loc)
				deallocate(vec_new_loc)
			end if
		end do
		vec_old = vec_new
	end do
	
	! ! vec_old = RandVectInR
    ! ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    ! ! groupn=block_o%col_group  
    ! ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    ! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	! ! call gemm_omp(fullmat,vec_old(1:nn,1:num_vectorsR),vec_new(1:mm,1:num_vectorsR),mm,nn,num_vectorsR)
	! ! vec_old = vec_new
	! ! call gemm_omp(matrixtmp1,vec_old(1:mm,1:num_vectorsR),vec_new(1:mm,1:num_vectorsR),mm,mm,num_vectorsR)
	
	
	
	deallocate(vec_old)
	allocate(RandVectOutR(DimMax,num_vectorsR))
	RandVectOutR = vec_new
	deallocate(vec_new)
	

! next consider the left multiplied vectors
	Nsub = 2**levelm
    Ng = 2**level_butterfly/Nsub
	num_vectorsL=Nsub*dimension_rank
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    groupn=block_o%col_group  
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	DimMax = max(mm,nn)

	allocate(RandVectInL(DimMax,num_vectorsL))
	RandVectInL=0
	allocate(vec_old(DimMax,num_vectorsL))
	allocate(vec_new(DimMax,num_vectorsL))

	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	do i=1, num_blocks
		header_m=basis_group(groupm_start+i-1)%head
		tailer_m=basis_group(groupm_start+i-1)%tail
		mm=tailer_m-header_m+1
		k=header_m-header_mm
		
		
		 !$omp parallel do default(shared) private(ii,jj)
		 do jj=idx_start,idx_start+dimension_rank-1
			 do ii=1, mm
				 RandVectInL(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo

	! get the left multiplied vectors
  
	! vec_old = RandVectInL
    ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    ! groupn=block_o%col_group  
    ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	! ! write(*,*)'dddd',mm,num_vectorsL
	! call gemmTN_omp(matrixtmp1,vec_old(1:mm,1:num_vectorsL),vec_new(1:mm,1:num_vectorsL),mm,mm,num_vectorsL)
	! ! write(*,*)'eeee'
	
	vec_old = RandVectInL
	do level = level_c+1,Maxlevel_for_blocks+1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vectorsL,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vectorsL),vec_new(idx_start_loc:idx_end_loc,1:num_vectorsL),ctemp1,ctemp2)
			else 
				DimMax_loc=basis_group(groupm_diag)%tail-basis_group(groupm_diag)%head+1 
				allocate(vec_old_loc(idx_end_loc-idx_start_loc+1,num_vectorsL))
				allocate(vec_new_loc(idx_end_loc-idx_start_loc+1,num_vectorsL))
				vec_old_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsL) = vec_old(idx_start_loc:idx_end_loc,1:num_vectorsL)
				vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsL) = vec_new(idx_start_loc:idx_end_loc,1:num_vectorsL)				
				
				
				
				! call Butterfly_Partial_MVP_dat_memeff(cascading_factors(level)%matrices_block_inverse(ii),'T',0,cascading_factors(level)%matrices_block_inverse(ii)%level_butterfly+1,&
				! &DimMax_loc,num_vectorsL,vec_old_loc,vec_new_loc)
				
				call OneL_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vectorsL,vec_old_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsL),vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsL))
				
				! write(*,*)'L:really?',fnorm(vec_new_loc,idx_end_loc-idx_start_loc+1,num_vectorsL)
				
				vec_new(idx_start_loc:idx_end_loc,1:num_vectorsL) = vec_new_loc(1:idx_end_loc-idx_start_loc+1,1:num_vectorsL) !!!!+vec_old(idx_start_loc:idx_end_loc,1:num_vectorsL)
				deallocate(vec_old_loc)
				deallocate(vec_new_loc)		
			end if
		end do
		vec_old = vec_new
	end do	
	
	vec_old = vec_new
    call Butterfly_Partial_MVP_dat_memeff(block_o,'T',0,level_butterfly+1,DimMax,num_vectorsL,vec_old,vec_new)      



    ! ! vec_old = RandVectInL
    ! ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    ! ! groupn=block_o%col_group  
    ! ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    ! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	! ! call gemmTN_omp(matrixtmp1,vec_old(1:mm,1:num_vectorsL),vec_new(1:mm,1:num_vectorsL),mm,mm,num_vectorsL)
	! ! vec_old = vec_new
	! ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    ! ! groupn=block_o%col_group  
    ! ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
    ! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	! ! call gemmTN_omp(fullmat,vec_old(1:mm,1:num_vectorsL),vec_new(1:nn,1:num_vectorsL),nn,mm,num_vectorsL)

	
	deallocate(vec_old)
	allocate(RandVectOutL(DimMax,num_vectorsL))
	RandVectOutL = vec_new
	deallocate(vec_new)	

	mem_vec=0
	mem_vec =mem_vec + SIZEOF(RandVectInL)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectInR)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectOutL)/1024.0d3
	mem_vec =mem_vec + SIZEOF(RandVectOutR)/1024.0d3
	Memory_int_vec = max(Memory_int_vec,mem_vec)
	
	
    return                

end subroutine OneL_Get_Randomized_Vectors_Sblock_Memeff


subroutine OneL_block_MVP_Sblock_dat(level_c,rowblock,trans,num_vect_sub,Vin,Vout)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character trans
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
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
	
	if(trans=='N')then
		block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
		  
		level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
		num_blocks=2**level_butterfly
		allocate (RandomVectors_InOutput_tmp(3))
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		allocate (RandomVectors_InOutput_tmp(1)%vector(nn,num_vect_sub))
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		allocate (RandomVectors_InOutput_tmp(2)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput_tmp(3)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput_tmp(ii)%vector = 0
		end do	 
		 
		groupn_start=groupn*2**(level_butterfly)
		header_nn=basis_group(groupn_start)%head
		idx_start = 1
		
		RandomVectors_InOutput_tmp(1)%vector = Vin
		
		! get the right multiplied vectors
		idx_start_glo = basis_group(groupm)%head
		random1=>RandomVectors_InOutput_tmp(1)
		random2=>RandomVectors_InOutput_tmp(2)
		ctemp1=1.0d0 ; ctemp2=0.0d0
	n1 = OMP_get_wtime()  
	  call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1	
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
		allocate(vec_old(mm,num_vect_sub))
		allocate(vec_new(mm,num_vect_sub))
		vec_old = RandomVectors_InOutput_tmp(2)%vector

		do level = Maxlevel_for_blocks+1,level_c+1,-1
			N_diag = 2**(level-level_c-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			
			if(level==Maxlevel_for_blocks+1)then
				n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
				! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

					
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					
					! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
					! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
					&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
				end do
				! !$omp end parallel do
				
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
			else 
				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

					
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					
					! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
					
					! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

					call OneL_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)	
				end do		
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1			
			end if
			
			vec_old = vec_new
		end do
		! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
		RandomVectors_InOutput_tmp(3)%vector = vec_new
		deallocate(vec_old)
		deallocate(vec_new)
		
		Vout = RandomVectors_InOutput_tmp(3)%vector
		
		! !$omp parallel do default(shared) private(i)
		do i=1, 3
			deallocate (RandomVectors_InOutput_tmp(i)%vector)
		enddo
		! !$omp end parallel do
		deallocate (RandomVectors_InOutput_tmp)		
		
	else 
		ctemp1=1.0d0 ; ctemp2=0.0d0	
		block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

		level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
		num_blocks=2**level_butterfly
		allocate (RandomVectors_InOutput_tmp(3))

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		allocate (RandomVectors_InOutput_tmp(3)%vector(nn,num_vect_sub))
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		allocate (RandomVectors_InOutput_tmp(1)%vector(mm,num_vect_sub))
		allocate (RandomVectors_InOutput_tmp(2)%vector(mm,num_vect_sub))
		do ii =1,3
			RandomVectors_InOutput_tmp(ii)%vector = 0
		end do	 
		 
		groupm_start=groupm*2**(level_butterfly)
		header_mm=basis_group(groupm_start)%head
		idx_start = 1
		
		RandomVectors_InOutput_tmp(1)%vector = Vin
		
		! get the left multiplied vectors
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
		idx_start_glo = basis_group(groupm)%head		
		allocate(vec_old(mm,num_vect_sub))
		allocate(vec_new(mm,num_vect_sub))	
		vec_old = RandomVectors_InOutput_tmp(1)%vector
		do level = level_c+1,Maxlevel_for_blocks+1
			N_diag = 2**(level-level_c-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			
			if(level==Maxlevel_for_blocks+1)then
				n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
				! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
					&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
				end do
				! !$omp end parallel do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
			else 
				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
					call OneL_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))			
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1	
			end if
			vec_old = vec_new
		end do	
		! ! write(*,*)vec_new(1,1),RandomVectors_InOutput_tmp(4)%vector(1,1)
		RandomVectors_InOutput_tmp(2)%vector = vec_new
		deallocate(vec_old)
		deallocate(vec_new)	
		
		
		random1=>RandomVectors_InOutput_tmp(2)
		random2=>RandomVectors_InOutput_tmp(3)
		
		n1 = OMP_get_wtime()
		call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1		
		
		
		Vout = RandomVectors_InOutput_tmp(3)%vector
		

		! !$omp parallel do default(shared) private(i)
		do i=1, 3
			deallocate (RandomVectors_InOutput_tmp(i)%vector)
		enddo
		! !$omp end parallel do
		deallocate (RandomVectors_InOutput_tmp)			
	end if	
	
    return                

end subroutine OneL_block_MVP_Sblock_dat

subroutine OneL_Get_Randomized_Vectors_LL_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
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
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_right_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	
	! write(*,*)'1a'
    ctemp1=1.0d0 ; ctemp2=0.0d0	
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**(level_butterfly)
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	! write(*,*)'1b'
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
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! random_complex_number()	! 
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 ! deallocate(matrixtemp1)
				 
			 ! end if
		end do
		!$omp end parallel do
	end do
	! write(*,*)'1c'
	! get the left multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	idx_start_glo = basis_group(groupm)%head		
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))	
	vec_old = RandomVectors_InOutput(1)%vector
	do level = level_c+1,Maxlevel_for_blocks+1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		
		if(level==Maxlevel_for_blocks+1)then
			n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
			end do
			! !$omp end parallel do
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1
		else 
			n1 = OMP_get_wtime()
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
				call OneL_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))			
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
			end do
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1	
		end if
		vec_old = vec_new
	end do	
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(4)%vector(1,1)
	RandomVectors_InOutput(2)%vector = vec_new
	deallocate(vec_old)
	deallocate(vec_new)	
	
 
	! write(*,*)'1d'
    random1=>RandomVectors_InOutput(2)
    random2=>RandomVectors_InOutput(3)
	
	n1 = OMP_get_wtime()
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1		
	
	! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)
	
		! write(*,*)'1e'
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

end subroutine OneL_Get_Randomized_Vectors_LL_Sblock



subroutine OneL_Get_Randomized_Vectors_RR_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	level_left_start= floor_safe(level_butterfly/2d0)+1
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub
	
	
	
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
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
				call RandomMat(nn,num_vect_subsub,min(nn,num_vect_subsub),RandomVectors_InOutput(1)%vector(1+k:nn+k,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s)*num_vect_subsub+num_vect_subsub),0)

				! ! !$omp parallel do default(shared) private(ii,jj)
				 ! do jj=1,num_vect_subsub
					 ! do ii=1, nn
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=matrixtemp1(jj,ii) ! 	random_complex_number() ! 
					 ! enddo
				 ! enddo
				 ! ! !$omp end parallel do
				 
				 ! deallocate(matrixtemp1)
			 ! end if
		end do
		!$omp end parallel do
	end do
	
	! get the right multiplied vectors
	idx_start_glo = basis_group(groupm)%head
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
	
	
n1 = OMP_get_wtime()  
  call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
n2 = OMP_get_wtime()
! time_tmp = time_tmp + n2 - n1	
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))
	vec_old = RandomVectors_InOutput(2)%vector

	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		
		if(level==Maxlevel_for_blocks+1)then
			n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
			end do
			! !$omp end parallel do
			
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1
		else 
			n1 = OMP_get_wtime()
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call OneL_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)	
			end do		
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1			
		end if
		
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	RandomVectors_InOutput(3)%vector = vec_new

	deallocate(vec_old)
	deallocate(vec_new)


	! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector)
	
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

end subroutine OneL_Get_Randomized_Vectors_RR_Sblock
subroutine OneL_Sblock_LowRank(level_c,rowblock)

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
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    call assert(level_butterfly==0,'Butterfly_Sblock_LowRank only works with LowRank blocks')
	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	
	
	
	
	num_blocks=2**level_butterfly

	
	num_vect_sub = size(block_o%ButterflyU(1)%matrix,2)
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   

	! get the right multiplied vectors
	idx_start_glo = basis_group(groupm)%head
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))
	vec_old = block_o%ButterflyU(1)%matrix

	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		
		if(level==Maxlevel_for_blocks+1)then
			n1 = OMP_get_wtime() ! comment: will this omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
			end do
			! !$omp end parallel do
			
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1
		else 
			n1 = OMP_get_wtime()
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call OneL_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)	
			end do		
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1			
		end if
		
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	block_o%ButterflyU(1)%matrix = vec_new
	deallocate(vec_old)
	deallocate(vec_new)


	
	
    return                

end subroutine OneL_Sblock_LowRank



subroutine OneL_Get_Randomized_Vectors_RR_Test_Sblock(level_c,rowblock,num_vect_sub)

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
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	
	
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(3)%vector(mm,num_vect_sub))
	do ii =1,3
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
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
				 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	
			 enddo
		 enddo
		 ! !$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	
	! RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector
	
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
	idx_start_glo = basis_group(groupm)%head
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))
	vec_old = RandomVectors_InOutput(2)%vector

	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			! write(*,*)level,ii
			groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			
			! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
			
			! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
			else 
				call OneL_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)	
			end if
		end do
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	RandomVectors_InOutput(3)%vector = vec_new
	deallocate(vec_old)
	deallocate(vec_new)


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

end subroutine OneL_Get_Randomized_Vectors_RR_Test_Sblock



subroutine OneL_Reconstruction_LL_Sblock(level_c,rowblock)
    
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
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)

	n1 = OMP_get_wtime()									  
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	
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
 
	n2 = OMP_get_wtime()						   
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
			call OneL_Get_Randomized_Vectors_LL_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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
	
	n1 = OMP_get_wtime()
	random=>Random_Block(1)													
	call Delete_RandVect('T',random,level_butterfly)
 	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1			
    return
    
end subroutine OneL_Reconstruction_LL_Sblock





subroutine OneL_Reconstruction_RR_Sblock(level_c,rowblock,error)
    
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
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	integer,allocatable::perms(:)
		
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
		
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
	! ! call assert(num_vectors==2**level_butterfly*dimension_rank,'incorrect num_vectors')
	
	! call assert(num_vectors==Nsub*dimension_rank,'incorrect num_vectors') !  check here later
	
    ! ! allocate (Random_Block(1))   !  check here later 
    
	n1 = OMP_get_wtime()									
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    
	
	Nbind = 1
	num_vect_sub = num_vect_subsub*Nbind
	
    random=>Random_Block(1)
	call Init_RandVect_Empty('N',random,num_vect_sub)

    level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    ! level_left_start = 0
	
	! call Zero_Butterfly(level_left_start,level_butterfly+1)

  	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1																								
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
			call OneL_Get_Randomized_Vectors_RR_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do
	
	! deallocate(perms)

	n1 = OMP_get_wtime()

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	deallocate(Random_Block)

	call OneL_Test_Error_RR_Sblock(level_c,rowblock,error)

 	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1		
    return
    
end subroutine OneL_Reconstruction_RR_Sblock






subroutine OneL_Test_Error_RR_Sblock(level_c,rowblock,error)

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
	

	
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	num_blocks=2**level_butterfly
	! write(*,*)level_butterfly,'heiyou',maxlevel_for_blocks,block_o%level
	
	allocate (Random_Block(1))
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    

	num_vect = 1
    random=>Random_Block(1)	
	
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_Output_ref(mm,num_vect))		
	
	call Init_RandVect_Empty('N',random,num_vect)	

	call OneL_Get_Randomized_Vectors_RR_Test_Sblock(level_c,rowblock,num_vect)

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

end subroutine OneL_Test_Error_RR_Sblock




subroutine OneL_block_MVP_inverse_dat(level,ii,trans,N,num_vect_sub,Vin,Vout)
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
   integer groupn,groupm,mm,nn

   ctemp1=1.0d0
   ctemp2=0.0d0
   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin
   
	block_off1 => cascading_factors(level)%BP(ii*2-1)%LL(1)%matrices_block(1)	
	block_off2 => cascading_factors(level)%BP(ii*2)%LL(1)%matrices_block(1)
	
	groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	call assert(mm+nn==N,'mm+nn/=N')  
		
   if(schurinv==0)then
		write(*,*)'schurinv=0 removed '
		stop
		! block_o => cascading_factors(level)%matrices_block_inverse(ii)
		! call butterfly_block_MVP_randomized_dat(block_o,trans,N,N,num_vect_sub,&
		! &Vin(1:N,1:num_vect_sub),Vout(1:N,1:num_vect_sub),ctemp1,ctemp2)
		! Vout(1:N,1:num_vect_sub) = Vin(1:N,1:num_vect_sub) + Vout(1:N,1:num_vect_sub)
   else 
		block_o => cascading_factors(level)%BP_inverse_schur(ii)%LL(1)%matrices_block(1)	
		if(trans=='N')then
			call butterfly_block_MVP_randomized_dat(block_off1,trans,mm,nn,num_vect_sub,&
			&Vin(1+mm:N,1:num_vect_sub),Vout(1:mm,1:num_vect_sub),ctemp1,ctemp2)
			Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)- Vout(1:mm,1:num_vect_sub)
			Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub)
			
			! write(2111,*)abs(Vout)
			
			call butterfly_block_MVP_randomized_dat(block_o,trans,mm,mm,num_vect_sub,&
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
			
			call butterfly_block_MVP_randomized_dat(block_o,trans,mm,mm,num_vect_sub,&
			&Vout(1:mm,1:num_vect_sub),Vin(1:mm,1:num_vect_sub),ctemp1,ctemp2)				
			Vin(1:mm,1:num_vect_sub) = Vout(1:mm,1:num_vect_sub) + Vin(1:mm,1:num_vect_sub)
			Vin(1+mm:N,1:num_vect_sub) = Vout(1+mm:N,1:num_vect_sub) 
			
			call butterfly_block_MVP_randomized_dat(block_off1,trans,mm,nn,num_vect_sub,&
			&Vin(1:mm,1:num_vect_sub),Vout(1+mm:N,1:num_vect_sub),ctemp1,ctemp2)
			Vout(1+mm:N,1:num_vect_sub) = Vin(1+mm:N,1:num_vect_sub) - Vout(1+mm:N,1:num_vect_sub)
			Vout(1:mm,1:num_vect_sub) = Vin(1:mm,1:num_vect_sub)
			
		end if
   end if


   Vin = Vin_tmp
   deallocate(Vin_tmp)
   
end subroutine OneL_block_MVP_inverse_dat



subroutine MultiLInitialize_Butterfly_Sblock(level_c,rowblock,kover)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,ll,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_rank_dummy, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj,rankmax
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j
	type(matrixblock),pointer::block
	
	
    block =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1)
    allocate (butterfly_block_randomized(1))
    
    level_butterfly=int((maxlevel_for_blocks-block%level)/2)*2
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    

    num_blocks=2**level_butterfly

	
	if(TwoLayerOnly==1)then
		rankmax = 0
		do ll=1,1 !cascading_factors(level_c)%BP(rowblock)%Lplus
			rankmax = max(rankmax, cascading_factors(level_c)%BP(rowblock)%LL(ll)%rankmax)
			! write(*,*)level_c,rowblock*2-1,ll,cascading_factors(level_c)%BP(rowblock*2-1)%LL(ll)%rankmax,'d1'
		end do	
	else 
		rankmax = 0
		do ll=1,cascading_factors(level_c)%BP(rowblock)%Lplus
			rankmax = max(rankmax, cascading_factors(level_c)%BP(rowblock)%LL(ll)%rankmax)
			! write(*,*)level_c,rowblock*2-1,ll,cascading_factors(level_c)%BP(rowblock*2-1)%LL(ll)%rankmax,'d1'
		end do		
	end if
	
	
	! write(*,*)'wocaca',block%level,block%row_group,block%col_group,block%rankmax
	! dimension_rank= max(ranktmp_glo,block%rankmax + kover) !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 

	dimension_rank= rankmax *1.2d0**(kover) !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank_dummy = 1
	
	! write(*,*)dimension_rank,ranktmp_glo,block%rankmax + kover,'zao'
	
	! if(level_c==2)dimension_rank=9
	! if(level_c==1)dimension_rank=max(dimension_rank,maxlevel_for_blocks)+kover

	 ! if(level_c==1)dimension_rank = 63
	
	! write(*,*)dimension_rank,'ha1'
	
    groupm=block%row_group         ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    ! ! if (int(mm/num_blocks)<dimension_rank) then
        ! ! dimension_rank=int(mm/num_blocks)
    ! ! endif
    butterfly_block_randomized(1)%dimension_rank=dimension_rank
    !write (*,*) dimension_rank, int(real(kk)/num_blocks/2.0), mm
    
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))
    ! allocate (butterfly_block_randomized(1)%ButterflyU_old(2**level_butterfly))
    ! allocate (butterfly_block_randomized(1)%ButterflyV_old(2**level_butterfly))
	! allocate (butterfly_block_randomized(1)%ButterflyV_qr(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyUInv(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyVInv(2**level_butterfly))
    
	dimension_max = 2*dimension_rank
	do blocks=1, num_blocks	
		dimension_m=size(block%ButterflyU(blocks)%matrix,1)
		dimension_n=size(block%ButterflyV(blocks)%matrix,1)
		dimension_max = max(dimension_max,dimension_m)	
		dimension_max = max(dimension_max,dimension_n)	
	end do	
	allocate(butterfly_block_randomized(1)%KerInv(dimension_max,dimension_max))
	call RandomMat(dimension_max,dimension_max,dimension_max,butterfly_block_randomized(1)%KerInv,3)	
	
    do blocks=1, num_blocks
        dimension_m=size(block%ButterflyU(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank_dummy))
        ! allocate (butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyUInv(blocks)%matrix(dimension_rank,dimension_m))
        ! write(*,*)shape(butterfly_block_randomized(1)%ButterflyUInv(blocks)%matrix),'hi1'
		
		allocate(matrixtemp1(dimension_rank_dummy,dimension_m))
		call RandomMat(dimension_rank_dummy,dimension_m,min(dimension_m,dimension_rank_dummy),matrixtemp1,0)
        do j=1, dimension_rank_dummy
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyU(blocks)%matrix
		deallocate(matrixtemp1)

		
        dimension_n=size(block%ButterflyV(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank_dummy))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_qr(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        ! allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%list(dimension_n))
		! butterfly_block_randomized(1)%ButterflyV(blocks)%list = block%ButterflyV(blocks)%list
		

		allocate(matrixtemp1(dimension_rank_dummy,dimension_n))
		call RandomMat(dimension_rank_dummy,dimension_n,min(dimension_n,dimension_rank_dummy),matrixtemp1,0)
        do j=1, dimension_rank_dummy
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyV(blocks)%matrix
		deallocate(matrixtemp1)

    enddo
	
	
	! ! write(*,*)'helloo'
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
					
                    ! ! ! ! !$omp parallel do default(shared) private(ii,jj)
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
                    ! ! ! ! !$omp end parallel do
                    ! allocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(index_i,index_j)%matrix(2*dimension_rank,2*dimension_rank))
                ! enddo
            ! enddo
            ! call invert_Butterfly_Kernel(level)
        enddo
        deallocate (matrixtemp1)
    endif	
	
    return

end subroutine MultiLInitialize_Butterfly_Sblock

subroutine MultiLSblock_randomized_memfree(level_c,rowblock,Memory)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   use Butterfly_compress_forward
   implicit none

    type(blockplus),pointer::bplus
	integer:: ii,ll,bb,rank_new_max
    real*8 Memory,rtemp,error,n2,n1	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall,M,N
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	complex(kind=8) ctemp, ctemp1, ctemp2, ctemp3, ctemp4		
	integer level_c,rowblock,idx_start_n_ref,idx_start_m_ref
	type(matrixblock)::block_tmp
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:),matsub_tmp(:,:)
	integer mm,nn
	! write(*,*)'ddd111'
	ctemp3=-1.0d0 ; ctemp4=1.0d0

	bplus =>  cascading_factors(level_c)%BP(rowblock)	
	Memory = 0	
	call assert(cascading_factors(level_c)%BP(rowblock)%Lplus>=2,'this is not a multi Bplus in MultiLSblock_randomized_memfree')
	! write(*,*)'ddd'
	call Initialize_Bplus_FromInput(bplus)
	! write(*,*)'dddddd'
	
	n1 = OMP_get_wtime()
	error_cnt = 0
	error_avr_glo = 0
	rankthusfarS = 0
	do bb =1,Bplus_randomized(1)%LL(2)%Nbound
		! write(*,*)bb,Bplus_randomized(1)%LL(2)%Nbound,'dddd'
		call MultiLrandomized_Onesubblock_Sblock(level_c,rowblock,bb)
		! write(*,*)'go'
	end do
	! call Test_Error_RR_Inner_Exact(Bplus)
	n2 = OMP_get_wtime()
	Time_Oneblock_forward = Time_Oneblock_forward + n2-n1

	! if(level_c<=3)then
		! write(*,*)'explicit computing butterfly for Sblock '

		! mm = basis_group(bplus%row_group)%tail - basis_group(bplus%row_group)%head + 1
		! nn = basis_group(bplus%col_group)%tail - basis_group(bplus%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(nn,nn))
		! ! allocate(matout(mm,mm))
		! allocate(matsub_glo(mm,nn))
		! matin = 0
		! do ii=1,nn
			! matin(ii,ii) = 1
		! end do
! ! write(*,*)'ddd1'
		! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',nn,matin,matsub_glo)
! ! write(*,*)'ddd2'		
		! call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,nn,matin,matsub_glo,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)		

			
		! block_tmp%row_group = bplus%LL(1)%matrices_block(1)%row_group
		! block_tmp%col_group = bplus%LL(1)%matrices_block(1)%col_group
		! block_tmp%level = bplus%LL(1)%matrices_block(1)%level
		! block_tmp%style = bplus%LL(1)%matrices_block(1)%style
		
		! idx_start_m_ref = basis_group(bplus%row_group)%head
		! idx_start_n_ref = basis_group(bplus%col_group)%head
		! write(*,*)'ddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		! write(*,*)'wocaonimaSblockOutter',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! deallocate(matin)
		! deallocate(matsub_glo)
		! call delete_blocks(block_tmp)	



		! mm = basis_group(bplus%row_group)%tail - basis_group(bplus%row_group)%head + 1
		! nn = basis_group(bplus%col_group)%tail - basis_group(bplus%col_group)%head + 1
			! ! write(*,*)'e',mm
		! allocate(matin(mm,mm))
		! ! allocate(matout(mm,mm))
		! allocate(matsub_glo(mm,nn))
		! allocate(matsub_tmp(nn,mm))
		! matin = 0
		! do ii=1,mm
			! matin(ii,ii) = 1
		! end do
! write(*,*)'dddd1'
		! call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',mm,matin,matsub_tmp)
! write(*,*)'dddd2'		
		! call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'T',mm,nn,mm,matin,matsub_tmp,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)		
		! call copymatT_OMP(matsub_tmp,matsub_glo,nn,mm)
		! deallocate(matsub_tmp)
			
		! block_tmp%row_group = bplus%LL(1)%matrices_block(1)%row_group
		! block_tmp%col_group = bplus%LL(1)%matrices_block(1)%col_group
		! block_tmp%level = bplus%LL(1)%matrices_block(1)%level
		! block_tmp%style = bplus%LL(1)%matrices_block(1)%style
		
		! idx_start_m_ref = basis_group(bplus%row_group)%head
		! idx_start_n_ref = basis_group(bplus%col_group)%head
! write(*,*)'dddd3'
		! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_n_ref)
		
		! write(*,*)'wocaonimaSblockOutter',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,fnorm(matsub_glo,mm,nn)

		! deallocate(matin)
		! deallocate(matsub_glo)
		! call delete_blocks(block_tmp)	
		! write(*,*)'Done computing ' 
	! end if
	
	


	
	
	call MultiLrandomized_Outter_Sblock_memfree(level_c,rowblock,rtemp)
	
	rank_new_max = 0
	do ll=1,Bplus_randomized(1)%Lplus
		rank_new_max = max(rank_new_max,Bplus_randomized(1)%LL(ll)%rankmax)
	end do
	
	write(*,'(A10,I5,A6,I3,A8,I3,A7,Es14.7)')'Mult No. ',rowblock,' rank:',rank_new_max,' L_butt:',bplus%LL(1)%matrices_block(1)%level_butterfly,' error:',error_avr_glo/error_cnt
			
	! write(*,*)'dd1'		
	call delete_Bplus(bplus)
	! write(*,*)'dd11'
	call copy_delete_Bplus(Bplus_randomized(1),bplus,Memory)
	! call copy_Bplus(Bplus_randomized(1),bplus,Memory)
	! ! write(*,*)'dd2'
	! call delete_Bplus(Bplus_randomized(1))
	! write(*,*)'dd3'
	deallocate(Bplus_randomized)
			
	
    return

end subroutine MultiLSblock_randomized_memfree







subroutine MultiLrandomized_Onesubblock_Sblock(level_c,rowblock,bb_o)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   use Butterfly_compress_forward
   implicit none

    type(blockplus),pointer::bplus
	integer:: ii,ll,bb,jj,bb_o,tt
    real*8 Memory,rtemp,error,n2,n1,mem_vec	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vout3(:,:),Vin(:,:)
	integer M,N,idx_start_n,idx_start_m,idx_start_n_loc,idx_end_n_loc,idx_start_m_loc,idx_end_m_loc,mm,nn,rmax,rank,idx_start_n_ref,idx_start_m_ref,idx_end_n_ref,idx_end_m_ref
	complex(kind=8)::ctemp1,ctemp2,Ctemp
	type(matrixblock),pointer::blocks
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:)
	real*8, allocatable :: Singular(:)
	integer level_c,rowblock,Nactive
	integer,allocatable::boxindex(:)
	integer Chunksize, Nchunk, Nidx, idx_s,cc
	
	bplus =>  cascading_factors(level_c)%BP(rowblock)
	
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
	
		rmax = max(bplus%LL(2)%rankmax*2,rankthusfarS)*2**(tt-1) !+ (tt-1)*10  !!!!! be careful here
		
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
			call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',Nidx,RandVectInR,RandVectOutR)
			
			
			
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
			call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',Nidx,RandVectInL,RandVectOutL)
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
		
		call RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU_glo,matV_glo,Singular,mm,nn,rmax,rank,LS_tolerance,SVD_tolerance_factor)				
		rankthusfarS = max(rankthusfarS,rank)
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
		! call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,1,RandVectInR,RandVectOutR,ctemp1,ctemp2)
		call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',1,RandVectInR,RandVectOutR)
		
		Vout1 = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:1)
		deallocate(RandVectInR)
		deallocate(RandVectOutR)
		! write(*,*)'yani 5'

		call gemm_omp(matV_glo(1:rank,1:nn),Vin,Vout3,rank,nn,1)
		call gemm_omp(matU_glo(1:mm,1:rank),Vout3,Vout2,mm,rank,1)
		
		
		error = fnorm(Vout2-Vout1,mm,1)/fnorm(Vout1,mm,1)
		! write(*,*)error,bb_o,'ninini'				
		
		deallocate(Vin)
		deallocate(Vout1)
		deallocate(Vout2)
		deallocate(Vout3)
		
		! write(*,*)tt,rmax,rank,'Sblock',error
		
		if(error>iter_tolerance*1d-1)then
			if(min(mm,nn)==rmax)then
				write(*,*)tt,rmax,'Sblock',error,rank
				write(*,*)'no need to increase rmax, try increase RandomizedSVD tolerance'
				stop
			end if		
			deallocate(matU_glo)
			deallocate(matV_glo)
		else
			if(verboselevel>=2)write(*,'(A33,A15,I3,A8,I4,A8,I2,A7,Es14.7)')' Onesub ',' bb_o:',bb_o,' rank:',rank,'Ntrial:',tt,' error:',error		
			error_avr_glo = error_avr_glo + error
			error_cnt = error_cnt + 1
			! write(*,*)bb_o,error,rmax,rank,'good'
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
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1		
		deallocate(matU_glo,matV_glo)			
		
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
			! write(*,*)Nactive,'jiba'
			!$omp parallel do default(shared) private(bb,ii)
			do ii = 1,Nactive
				! if(basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
				bb = boxindex(ii)
				if(Bplus_randomized(1)%LL(ll+1)%Nbound==0)then
					! write(*,*)'666',ll
					
					call Butterfly_compress_N15_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),idx_start_m_ref,idx_start_n_ref)
					
					! blocks=>Bplus_randomized(1)%LL(ll)%matrices_block(bb)
					! M = basis_group(blocks%row_group)%tail - basis_group(blocks%row_group)%head + 1
					! N = basis_group(blocks%col_group)%tail - basis_group(blocks%col_group)%head + 1
					! allocate(Vin(N,1))
					! Vin = 0
					! allocate(Vout1(M,1))
					! Vout1 = 0
					! call butterfly_block_MVP_randomized_dat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),'N',M,N,1,Vin,Vout1,ctemp1,ctemp2)
					! deallocate(Vin,Vout1)
					
					! write(*,*)'667'	
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
					
					! blocks=>Bplus_randomized(1)%LL(ll)%matrices_block(bb)
					! M = basis_group(blocks%row_group)%tail - basis_group(blocks%row_group)%head + 1
					! N = basis_group(blocks%col_group)%tail - basis_group(blocks%col_group)%head + 1
					! allocate(Vin(N,1))
					! Vin = 0
					! allocate(Vout1(M,1))
					! Vout1 = 0
					! call butterfly_block_MVP_randomized_dat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),'N',M,N,1,Vin,Vout1,ctemp1,ctemp2)
					! deallocate(Vin,Vout1)
					
				else
					! write(*,*)'777'				
					call Butterfly_compress_N15_withoutBoundary_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),Bplus_randomized(1)%LL(ll+1)%boundary_map,Nboundall,groupm_start, rtemp, idx_start_m_ref,idx_start_n_ref)
					
					! blocks=>Bplus_randomized(1)%LL(ll)%matrices_block(bb)
					! M = basis_group(blocks%row_group)%tail - basis_group(blocks%row_group)%head + 1
					! N = basis_group(blocks%col_group)%tail - basis_group(blocks%col_group)%head + 1
					! allocate(Vin(N,1))
					! Vin = 0
					! allocate(Vout1(M,1))
					! Vout1 = 0
					! call butterfly_block_MVP_randomized_dat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),'N',M,N,1,Vin,Vout1,ctemp1,ctemp2)
					! deallocate(Vin,Vout1)
					
					! write(*,*)'778'	
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
					
					! blocks=>Bplus_randomized(1)%LL(ll)%matrices_block(bb)
					! M = basis_group(blocks%row_group)%tail - basis_group(blocks%row_group)%head + 1
					! N = basis_group(blocks%col_group)%tail - basis_group(blocks%col_group)%head + 1
					! allocate(Vin(N,1))
					! Vin = 0
					! allocate(Vout1(M,1))
					! Vout1 = 0
					! call butterfly_block_MVP_randomized_dat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),'N',M,N,1,Vin,Vout1,ctemp1,ctemp2)
					! deallocate(Vin,Vout1)					
					
				end if				
				! end if	
			end do
			!$omp end parallel do
			
			do ii = 1,Nactive
				bb = boxindex(ii)
				Bplus_randomized(1)%LL(ll)%rankmax = max(Bplus_randomized(1)%LL(ll)%rankmax,Bplus_randomized(1)%LL(ll)%matrices_block(bb)%rankmax)
			end do
			
			deallocate(boxindex)
				
		! write(*,*)Bplus_randomized(1)%LL(ll)%rankmax,'wonima'
		end do
		
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1

		deallocate(matSub_glo)		
	end if

    return

end subroutine MultiLrandomized_Onesubblock_Sblock




subroutine MultiLrandomized_Outter_Sblock_memfree(level_c,rowblock,Memory)

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
	
	! Bplus =>  cascading_factors(level_c)%BP(rowblock)
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2

	do tt =1,10
		do ntry=1,1
		! write(*,*)'haha1'
		n1 = OMP_get_wtime()
		call MultiLInitialize_Butterfly_Sblock(level_c,rowblock,(tt-1))
		n2 = OMP_get_wtime()
		Time_Init_forward = Time_Init_forward + n2 -n1 
		
! write(*,*)'haha2'
		n1 = OMP_get_wtime()
		call MultiLReconstruction_LL_Outter_Sblock(level_c,rowblock)	
		! write(*,*)'haha3'
		call MultiLReconstruction_RR_Outter_Sblock(level_c,rowblock,error_inout)
		n2 = OMP_get_wtime()
		Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1
		! time_tmp = time_tmp + n2-n1
		
! write(*,*)'haha4'
		
		call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
		rank_new_max = butterfly_block_randomized(1)%rankmax
		
		write(*,*)tt,error_inout,butterfly_block_randomized(1)%rankmax,butterfly_block_randomized(1)%dimension_rank
		
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
									 
			if(verboselevel>=2)write(*,'(A34,A8,I3,A8,I2,A8,I3,A7,Es14.7)')' Outter: ',' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
			error_avr_glo = error_avr_glo + error_inout
			error_cnt = error_cnt + 1
			return
		end if
		end do
	end do
	write(*,*)'randomized scheme not converged in MultiLrandomized_Outter_Sblock_memfree',error_inout,rank_new_max
	stop
	
    return

end subroutine MultiLrandomized_Outter_Sblock_memfree


subroutine MultiLReconstruction_LL_Outter_Sblock(level_c,rowblock)
    
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
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)

	n1 = OMP_get_wtime()
	
    block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	
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
 
 	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1		 
 
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
			call MultiLGet_Randomized_Vectors_LL_Outter_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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
	
	n1 =OMP_get_wtime()
	random=>Random_Block(1)														
	call Delete_RandVect('T',random,level_butterfly)
 	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1			
	
    return
    
end subroutine MultiLReconstruction_LL_Outter_Sblock

subroutine MultiLReconstruction_RR_Outter_Sblock(level_c,rowblock,error)
    
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
	type(matrixblock),pointer::block_o
	
	n1 = OMP_get_wtime()
	
	block_o =>  cascading_factors(level_c)%BP(rowblock)%LL(1)%matrices_block(1)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
		
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	
	! ! call assert(num_vectors==2**level_butterfly*dimension_rank,'incorrect num_vectors')
	
	! call assert(num_vectors==Nsub*dimension_rank,'incorrect num_vectors') !  check here later
	
    ! ! allocate (Random_Block(1))   !  check here later 
    
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    
	
	Nbind = 1
	num_vect_sub = num_vect_subsub*Nbind
		! write(*,*)'ggod 1'
    random=>Random_Block(1)
	call Init_RandVect_Empty('N',random,num_vect_sub)

    level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    ! level_left_start = 0
	
	! call Zero_Butterfly(level_left_start,level_butterfly+1)
! write(*,*)'ggod 2'
    ! ! allocate(perms(Nsub))
	! ! call rperm(Nsub, perms)
	! ! do ii = 1,Nsub		
		! ! nth_s = perms(ii)
		! ! nth_e = perms(ii)

 	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1		
		
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
			
			! write(*,*)'ggod 3',ii
			n1 = OMP_get_wtime()
			call MultiLGet_Randomized_Vectors_RR_Outter_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			Time_Vector_forward = Time_Vector_forward + n2-n1
			! write(*,*)'ggod 4',ii
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1			
		end do
	end do
	
	! deallocate(perms)

	n1 = OMP_get_wtime()

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	deallocate(Random_Block)

	! write(*,*)'heij'
	call MultiLTest_Error_RR_Outter_Sblock(level_c,rowblock,error)
	! write(*,*)'heijd'

	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1	
	
    return
    
end subroutine MultiLReconstruction_RR_Outter_Sblock
subroutine MultiLTest_Error_RR_Outter_Sblock(level_c,rowblock,error)

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

	Bplus =>  cascading_factors(level_c)%BP(rowblock) 
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

	call MultiLGet_Randomized_Vectors_RR_Test_Outter_Sblock(level_c,rowblock,num_vect)

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

end subroutine MultiLTest_Error_RR_Outter_Sblock









subroutine MultiLGet_Randomized_Vectors_LL_Outter_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	type(blockplus),pointer::Bplus
	
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
	
	Bplus =>  cascading_factors(level_c)%BP(rowblock) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    ! allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (Vin(mm,num_vect_sub))
    allocate (Vout(nn,num_vect_sub))
	Vin = 0
	Vout = 0
	
 
	 
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
	
	call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'T',num_vect_sub,Vin,Vout)
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'T',mm,nn,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	
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

end subroutine MultiLGet_Randomized_Vectors_LL_Outter_Sblock








subroutine MultiLGet_Randomized_Vectors_RR_Outter_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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

	Bplus =>  cascading_factors(level_c)%BP(rowblock) 	
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
	 
	allocate(Vin(nn,num_vect_sub)) 
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
  ! write(*,*)'p1'
  ! call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)
  call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',num_vect_sub,Vin,Vout)
  ! write(*,*)'p2'
  call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
  n2 = OMP_get_wtime()
! time_tmp = time_tmp + n2 - n1	
   ! write(*,*)'p3'


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

end subroutine MultiLGet_Randomized_Vectors_RR_Outter_Sblock









subroutine MultiLGet_Randomized_Vectors_RR_Test_Outter_Sblock(level_c,rowblock,num_vect_sub)

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
		
	
	Bplus =>  cascading_factors(level_c)%BP(rowblock) 
	  
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
	allocate(Vin(nn,num_vect_sub)) 
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
	
	! write(*,*)'d1'
	
	call Bplus_block_MVP_Sblock_dat(level_c,rowblock,'N',num_vect_sub,Vin,Vout)
	! call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)		
	! write(*,*)'d2'
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,num_vect_sub,Vin,Vout,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	! write(*,*)'d3'


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

end subroutine MultiLGet_Randomized_Vectors_RR_Test_Outter_Sblock




subroutine Bplus_block_MVP_inverse_dat(level,ii,trans,N,num_vect_sub,Vin,Vout)
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

   ctemp1=1.0d0
   ctemp2=0.0d0
   allocate(Vin_tmp(N,num_vect_sub))
   Vin_tmp = Vin
   
	bplus_off1 => cascading_factors(level)%BP(ii*2-1)	
	bplus_off2 => cascading_factors(level)%BP(ii*2)
	
	groupn=bplus_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=bplus_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	call assert(mm+nn==N,'mm+nn/=N')  
		
   if(schurinv==0)then
		write(*,*)'schurinv=0 removed '
		stop
		! block_o => cascading_factors(level)%matrices_block_inverse(ii)
		! call butterfly_block_MVP_randomized_dat(block_o,trans,N,N,num_vect_sub,&
		! &Vin(1:N,1:num_vect_sub),Vout(1:N,1:num_vect_sub),ctemp1,ctemp2)
		! Vout(1:N,1:num_vect_sub) = Vin(1:N,1:num_vect_sub) + Vout(1:N,1:num_vect_sub)
   else 
		bplus_o => cascading_factors(level)%BP_inverse_schur(ii)
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
   end if


   Vin = Vin_tmp
   deallocate(Vin_tmp)
   
end subroutine Bplus_block_MVP_inverse_dat




subroutine Bplus_block_MVP_Sblock_dat(level_c,rowblock,trans,num_vect_sub,Vin,Vout)
	use MODULE_FILE
	! use lapack95
	! use blas95
	implicit none
	integer level_c, rowblock, num_vect_sub
	character trans
	complex(kind=8) :: Vin(:,:), Vout(:,:)
	complex(kind=8) :: ctemp1,ctemp2,Ctemp
	type(blockplus),pointer::bplus_o
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
	type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)
	
	real*8::n2,n1 	
	
	! write(*,*)'I am here'
	if(trans=='N')then
		bplus_o => cascading_factors(level_c)%BP(rowblock)	
		! allocate (RandomVectors_InOutput_tmp(3))
		groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		! allocate (RandomVectors_InOutput_tmp(1)%vector(nn,num_vect_sub))
		
		groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		! allocate (RandomVectors_InOutput_tmp(2)%vector(mm,num_vect_sub))
		! allocate (RandomVectors_InOutput_tmp(3)%vector(mm,num_vect_sub))
		! do ii =1,3
			! RandomVectors_InOutput_tmp(ii)%vector = 0
		! end do	 
		! RandomVectors_InOutput_tmp(1)%vector = Vin
   

		! get the right multiplied vectors
		idx_start_glo = basis_group(groupm)%head
		! random1=>RandomVectors_InOutput_tmp(1)
		! random2=>RandomVectors_InOutput_tmp(2)
		ctemp1=1.0d0 ; ctemp2=0.0d0
! write(*,*)'1 I am'
		call Bplus_block_MVP_randomized_dat(bplus_o,'N',mm,nn,num_vect_sub,Vin,Vout,ctemp1,ctemp2)
! write(*,*)'2 I am'
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
		! allocate(vec_old(mm,num_vect_sub))
		allocate(vec_new(mm,num_vect_sub))
		! vec_old = RandomVectors_InOutput_tmp(2)%vector
		
		do level = Maxlevel_for_blocks+1,level_c+1,-1
			N_diag = 2**(level-level_c-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			
			if(level==Maxlevel_for_blocks+1)then
				n1 = OMP_get_wtime()
				! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

					
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					
					! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
					! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head
! write(*,*)'3 I am'
					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
					&Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
! write(*,*)'4 I am'					
					! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
				end do
				! !$omp end parallel do
				
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
			else 
				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

					
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					
					! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
					
					! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head
! write(*,*)'5 I am'
					call Bplus_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,Vout(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
! write(*,*)'6 I am'					
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)	
				end do		
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1			
			end if
			
			Vout = vec_new
		end do
		! ! write(*,*)vec_new(1,1),RandomVectors_InOutput_tmp(2)%vector(1,1)
		! RandomVectors_InOutput_tmp(3)%vector = vec_new
		! deallocate(vec_old)
		deallocate(vec_new)

		! Vout = RandomVectors_InOutput_tmp(3)%vector
	   
		! ! !$omp parallel do default(shared) private(i)
		! do i=1, 3
			! deallocate (RandomVectors_InOutput_tmp(i)%vector)
		! enddo
		! ! !$omp end parallel do
		! deallocate (RandomVectors_InOutput_tmp)		   
   
   else if(trans=='T')then
		bplus_o => cascading_factors(level_c)%BP(rowblock)  
 		! allocate (RandomVectors_InOutput_tmp(3))
		groupn=bplus_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
		groupm=bplus_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		
		! allocate (RandomVectors_InOutput_tmp(1)%vector(mm,num_vect_sub))
		! allocate (RandomVectors_InOutput_tmp(2)%vector(mm,num_vect_sub))
		! allocate (RandomVectors_InOutput_tmp(3)%vector(nn,num_vect_sub))
		! do ii =1,3
			! RandomVectors_InOutput_tmp(ii)%vector = 0
		! end do	 
		! RandomVectors_InOutput_tmp(1)%vector = Vin  
   
   		ctemp1=1.0d0 ; ctemp2=0.0d0
   ! write(*,*)'1 I amm'
		! get the left multiplied vectors 
		idx_start_glo = basis_group(groupm)%head		
		allocate(vec_old(mm,num_vect_sub))
		allocate(vec_new(mm,num_vect_sub))	

			! if(level_c==1)then
			! do while(.true.)
				! Ctemp = 1d0
			! end do
			! end if
		
		vec_old = Vin
		! Vout = Vin
		do level = level_c+1,Maxlevel_for_blocks+1
			N_diag = 2**(level-level_c-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			
			if(level==Maxlevel_for_blocks+1)then
				n1 = OMP_get_wtime()
				! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii,'ggg'
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
					call fullmat_block_MVP_randomized_dat(cascading_factors(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
					&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
				end do
				! !$omp end parallel do
				n2 = OMP_get_wtime()
				! ! time_tmp = time_tmp + n2 - n1
			else 
				n1 = OMP_get_wtime()
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					! write(*,*)level,ii,'gggddd',idx_start_diag,idx_start_diag+N_diag-1
					groupm_diag = cascading_factors(level)%BP_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
					call Bplus_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))			
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
				end do
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1	
			end if
			vec_old = vec_new
		end do	
		! ! write(*,*)vec_new(1,1),RandomVectors_InOutput_tmp(4)%vector(1,1)
		! RandomVectors_InOutput_tmp(2)%vector = vec_new
		! deallocate(vec_old)
					  
		! write(*,*)'dddgegegege'
		deallocate(vec_new)
		! random1=>RandomVectors_InOutput_tmp(2)
		! random2=>RandomVectors_InOutput_tmp(3)
		   ! write(*,*)'2 I amm'
		n1 = OMP_get_wtime()
		call Bplus_block_MVP_randomized_dat(bplus_o,'T',mm,nn,num_vect_sub,vec_old,Vout,ctemp1,ctemp2)	
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1			
			   ! write(*,*)'3 I amm'
		! Vout = RandomVectors_InOutput_tmp(3)%vector
												
		   
												
	   
						 
										 
		
		deallocate(vec_old)	
		
		! ! !$omp parallel do default(shared) private(i)
		! do i=1, 3
			! deallocate (RandomVectors_InOutput_tmp(i)%vector)
		! enddo
		! ! !$omp end parallel do
		! deallocate (RandomVectors_InOutput_tmp)				
		
 
   end if
   
 
end subroutine Bplus_block_MVP_Sblock_dat



end module Bplus_rightmultiply
