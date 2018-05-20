module Bplus_rightmultiply
use Butterfly_rightmultiply
use Randomized_reconstruction
integer rankthusfarS					
contains





subroutine Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats)

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
	type(matrixblock),pointer::block_o
    integer::rank_new_max
	real*8::rank_new_avr,error,rate,rankrate_inner,rankrate_outter 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,rank,ntry,rank0,rank0_inner,rank0_outter
	real*8:: error_inout
	real*8:: n1,n2,Memory
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	
	Memory = 0
	! write(*,*)'caca',level_c,rowblock
    bplus =>  ho_bf1%levels(level_c)%BP(rowblock)
	! write(*,*)'caca1',level_c,rowblock
	if(bplus%Lplus==1)then
	
		block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
												 
		level_butterfly=block_o%level_butterfly

		if(level_butterfly==0)then
			if(level_c>=ho_bf1%Maxlevel-1)then
				call OneL_Sblock_LowRank(ho_bf1,level_c,rowblock)
			else 		
				write(*,*)'unexpected level_c'
				stop
			end if	
		else 
			ho_bf1%ind_lv=level_c
			ho_bf1%ind_bk=rowblock
			rank0 = block_o%rankmax
			rate = 1.2d0
			call Butterfly_randomized(level_butterfly,rank0,rate,block_o,ho_bf1,butterfly_block_MVP_Sblock_dat,error_inout,'Sblock',option,stats)
							
#if PRNTlevel >= 1
				write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout	
#endif
			
		end if
	else 
	! write(*,*)'cacadd1',level_c,rowblock
		! call MultiLSblock_randomized_memfree(level_c,rowblock,Memory)
		
		ho_bf1%ind_lv=level_c
		ho_bf1%ind_bk=rowblock
		Bplus =>  ho_bf1%levels(level_c)%BP(rowblock)
		block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 		
		
		rank0_inner = Bplus%LL(2)%rankmax
		rankrate_inner = 2.0d0
		
		rank0_outter = block_o%rankmax
		rankrate_outter=1.2d0	
		level_butterfly = block_o%level_butterfly
		call Buplus_randomized(level_butterfly,Bplus,ho_bf1,rank0_inner,rankrate_inner,Bplus_block_MVP_Sblock_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Sblock_dat,error_inout,'Sblock+',option,stats)
		
		block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 	
#if PRNTlevel >= 1
		write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout	
#endif		
		
	end if
	
    return

end subroutine Bplus_Sblock_randomized_memfree




subroutine OneL_Sblock_LowRank(ho_bf1,level_c,rowblock)

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
	type(hobf)::ho_bf1
	
	
	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
	block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 
	  
    level_butterfly=block_o%level_butterfly
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
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
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
				groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call OneL_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
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



end module Bplus_rightmultiply
