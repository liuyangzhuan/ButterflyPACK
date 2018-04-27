module Butterfly_rightmultiply
use Utilites_randomized
! use Butterfly_compression_givenfullmat
use omp_lib
contains


subroutine Butterfly_Sblock_randomized_memfree(level_c,rowblock,Memory)

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
	
	Memory = 0
	
    block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	
	! ! ! call copy_butterfly(block_o,block_old)
	
	
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2

	if(level_butterfly==0)then
		if(level_c>=maxlevel_for_blocks-1)then
		
			call Butterfly_Sblock_LowRank(level_c,rowblock)
		
		
		
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
		! ! T0=secnds(0.0)
	
		do tt =1,10
			do ntry=1,1
			n1 = OMP_get_wtime()
			call Initialize_Butterfly_Sblock(block_o,level_c,tt-1)
			n2 = OMP_get_wtime()
			Time_Init_forward = Time_Init_forward + n2 -n1 
			
	
			n1 = OMP_get_wtime()
			call Reconstruction_LL_Sblock(level_c,rowblock)	
			call Reconstruction_RR_Sblock(level_c,rowblock,error_inout)
			n2 = OMP_get_wtime()
			Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1

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
				write(*,'(A8,I5,A6,I3,A8,I2,A8,I3,A7,Es14.7)')'     No.',rowblock,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
				return
			end if
			end do
		end do
		write(*,*)'randomized scheme not converged',error_inout
		stop
		
	end if
	
    return

end subroutine Butterfly_Sblock_randomized_memfree




subroutine Butterfly_Sblock_randomized_partialupdate_memfree(level_c,rowblock,Memory)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,ii_loc,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,level_butterfly_loc,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d,error_inout
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o,block_inverse
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng,unique_nth,level_left_start
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n, dimension_m
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,rank_new_max
	type(RandomBlock), pointer :: random
	real*8::n2,n1,Memory
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	if(level_butterfly==0)then
		call Butterfly_Sblock_randomized_memfree(level_c,rowblock,Memory)
	else 
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		 

		
		do level = Maxlevel_for_blocks+1,level_c+1,-1
			N_diag = 2**(level-level_c-1)
			idx_start_diag = (rowblock-1)*N_diag+1

			
			if(level==Maxlevel_for_blocks+1)then

				!$omp parallel do default(shared) private(ii,ii_loc,dimension_m,dimension_rank)
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					ii_loc = ii - idx_start_diag + 1
					dimension_m = size(block_o%ButterflyU(ii_loc)%matrix,1)
					dimension_rank = size(block_o%ButterflyU(ii_loc)%matrix,2)
					call gemm_omp(cascading_factors(level)%matrices_block_inverse(ii)%fullmat,block_o%ButterflyU(ii_loc)%matrix, block_o%ButterflyU(ii_loc)%matrix,dimension_m,dimension_m,dimension_rank)
				end do
				!$omp end parallel do
				
			else
				do ii = idx_start_diag,idx_start_diag+N_diag-1
					ii_loc = ii - idx_start_diag + 1
					level_butterfly_loc = Maxlevel_for_blocks+1-level    !!!! all even level butterfly could be problematic here
					
					! allocate(agent_block(1))
					n1 = OMP_get_wtime()
					call Extract_partial_butterfly(block_o,level_butterfly_loc,ii_loc,'L')
					call Initialize_Butterfly_SblockSmall(agent_block(1),0)
					n2 = OMP_get_wtime()
					Time_Init_forward = Time_Init_forward + n2 -n1 				

		
					n1 = OMP_get_wtime()
					call Reconstruction_LL_SblockSmall(level,ii,agent_block(1))	
					call Reconstruction_RR_SblockSmall(level,ii,agent_block(1),error_inout)
					n2 = OMP_get_wtime()
					Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1

				


					call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(1))
					rank_new_max = butterfly_block_randomized(1)%rankmax

					call copy_randomizedbutterfly_partial(butterfly_block_randomized(1),block_o,level_butterfly_loc,ii_loc,'L',Memory)
					
					if(level_butterfly_loc==level_butterfly)write(*,'(A8,I5,A6,I3,A8,I3,A7,Es14.7)')'     No.',rowblock,' rank:',rank_new_max,' L_butt:',level_butterfly_loc,' error:',error_inout
					call Delete_randomized_butterfly()
					call delete_blocks(agent_block(1))
					deallocate(agent_block)					
				end do

			end if

		end do
	end if
    return                

end subroutine Butterfly_Sblock_randomized_partialupdate_memfree




subroutine Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,LR)
	use misc
    use MODULE_FILE
    implicit none
	
	type(matrixblock)::block_o
	integer level_butterfly,level_butterfly_loc, ij_loc,index_i,index_i_start,index_j_start,index_j,level,ii,nn,mm,num_blocks,rank
	character LR
	
	allocate(agent_block(1))
	

	call assert(level_butterfly_loc>=1,'level_butterfly_loc cannot be zero')

	agent_block(1)%style = block_o%style
	agent_block(1)%level_butterfly = level_butterfly_loc
	agent_block(1)%rankmax = block_o%rankmax
	agent_block(1)%rankmin = block_o%rankmin
	level_butterfly = block_o%level_butterfly
	
	
	
	num_blocks=2**level_butterfly_loc




	allocate(agent_block(1)%ButterflyU(num_blocks))
	allocate(agent_block(1)%ButterflyV(num_blocks))
	
	allocate(agent_block(1)%ButterflyKerl(level_butterfly_loc))

	
	if(LR=='L')then
		do level=1, level_butterfly_loc
			agent_block(1)%ButterflyKerl(level)%num_row=2**level
			agent_block(1)%ButterflyKerl(level)%num_col=2**(level_butterfly_loc-level+1)
			allocate(agent_block(1)%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly_loc-level+1)))
			do index_i=1, 2**level
				do index_j=1, 2**(level_butterfly_loc-level)

					index_i_start = (ij_loc-1)*2**level	
					nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,2)
					rank=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,1)
					allocate(agent_block(1)%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn))
					agent_block(1)%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix = block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix

					nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,2)
					allocate(agent_block(1)%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn))
					agent_block(1)%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix = block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix

					if (level==level_butterfly_loc) then
						index_i_start = (ij_loc-1)*2**level	
						
						mm=size(block_o%ButterflyU(index_i+index_i_start)%matrix,1)
						rank=size(block_o%ButterflyU(index_i+index_i_start)%matrix,2)
						allocate(agent_block(1)%ButterflyU(index_i)%matrix(mm,rank))
						agent_block(1)%ButterflyU(index_i)%matrix = block_o%ButterflyU(index_i+index_i_start)%matrix					
					endif
				enddo
			enddo
			
			if(level==1)then
				do index_i=1, 1
					do index_j=1, 2**(level_butterfly_loc-level)
						index_i_start = (ij_loc-1)*2**level	
						nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,2)
						allocate(agent_block(1)%ButterflyV(2*index_j-1)%matrix(nn,nn))
						agent_block(1)%ButterflyV(2*index_j-1)%matrix = 0
						do ii=1,nn
							agent_block(1)%ButterflyV(2*index_j-1)%matrix(ii,ii)=1
						end do
						nn=size(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,2)
						allocate(agent_block(1)%ButterflyV(2*index_j)%matrix(nn,nn))
						agent_block(1)%ButterflyV(2*index_j)%matrix = 0
						do ii=1,nn
							agent_block(1)%ButterflyV(2*index_j)%matrix(ii,ii)=1
						end do					
					end do
				end do
			end if
			
		enddo
	else if(LR=='R')then
		do level=1, level_butterfly_loc
			agent_block(1)%ButterflyKerl(level)%num_row=2**level
			agent_block(1)%ButterflyKerl(level)%num_col=2**(level_butterfly_loc-level+1)
			allocate(agent_block(1)%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly_loc-level+1)))
			do index_i=1, 2**(level-1)
				do index_j=1, 2**(level_butterfly_loc-level+1)

					index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)	
					
					mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,1)
					rank=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,2)
					allocate(agent_block(1)%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm,rank))
					agent_block(1)%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix = block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix

					mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,1)
					allocate(agent_block(1)%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm,rank))
					agent_block(1)%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix = block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix

					if (level==1) then
						index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)	
						
						nn=size(block_o%ButterflyV(index_j+index_j_start)%matrix,1)
						rank=size(block_o%ButterflyV(index_j+index_j_start)%matrix,2)
						allocate(agent_block(1)%ButterflyV(index_j)%matrix(nn,rank))
						agent_block(1)%ButterflyV(index_j)%matrix = block_o%ButterflyV(index_j+index_j_start)%matrix					
					endif
				enddo
			enddo
			
			if(level==level_butterfly_loc)then
				do index_i=1, 2**(level_butterfly_loc-1)
					do index_j=1, 1
						index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)	
						mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,1)
						allocate(agent_block(1)%ButterflyU(2*index_i-1)%matrix(mm,mm))
						agent_block(1)%ButterflyU(2*index_i-1)%matrix = 0
						do ii=1,mm
							agent_block(1)%ButterflyU(2*index_i-1)%matrix(ii,ii)=1
						end do
						mm=size(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,1)
						allocate(agent_block(1)%ButterflyU(2*index_i)%matrix(mm,mm))
						agent_block(1)%ButterflyU(2*index_i)%matrix = 0
						do ii=1,mm
							agent_block(1)%ButterflyU(2*index_i)%matrix(ii,ii)=1
						end do				
					end do
				end do
			end if
			
		enddo
	
	end if	
end subroutine Extract_partial_butterfly






subroutine copy_randomizedbutterfly_partial(block_i,block_o,level_butterfly_loc,ij_loc,LR,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_o
type(butterflyblock_randomized)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,dimension_m,dimension_n,index_i_start,index_j_start
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,level_butterfly_loc,ij_loc	
character LR
real*8,optional::memory
if(present(memory))memory=0

!!!!! be careful here, may need changes later 
block_o%rankmax = max(block_o%rankmax,block_i%rankmax)
block_o%rankmin = max(block_o%rankmin,block_i%rankmin)


call assert(level_butterfly_loc>=1,'level_butterfly_loc cannot be zero')
call assert(level_butterfly_loc==block_i%level_butterfly,'level_butterfly_loc/=block_i%level_butterfly')

level_butterfly = block_o%level_butterfly
num_blocks=2**level_butterfly_loc


if(LR=='L')then

	do level=1, level_butterfly_loc
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly_loc-level)
				index_i_start = (ij_loc-1)*2**level	

				if(level==1)then
					dimension_n = size(block_i%ButterflyV(2*index_j-1)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)					
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix(rank,dimension_n))
					call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix, block_i%ButterflyV(2*index_j-1)%matrix, &
					&block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,nn,dimension_n)

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,dimension_n)))then
						write(*,*)'NAN in L 1'
					end if
					
					
					dimension_n = size(block_i%ButterflyV(2*index_j)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)					
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix(rank,dimension_n))
					call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix, block_i%ButterflyV(2*index_j)%matrix, &
					&block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,nn,dimension_n)

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,dimension_n)))then
						write(*,*)'NAN in L 2'
					end if					
					
					
				else 
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix(rank,nn))
					block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix				

					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,nn)))then
						write(*,*)'NAN in L 3'
					end if	
					
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)
					allocate(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix(rank,nn))
					block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix							
					
					if(isnan(fnorm(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,nn)))then
						write(*,*)'NAN in L 4'
					end if					
				
				end if
				
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix)/1024.0d3
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix)/1024.0d3

				if (level==level_butterfly_loc) then
					index_i_start = (ij_loc-1)*2**level	
					mm=size(block_i%ButterflyU(index_i)%matrix,1)
					rank=size(block_i%ButterflyU(index_i)%matrix,2)
					deallocate(block_o%ButterflyU(index_i+index_i_start)%matrix)
					allocate(block_o%ButterflyU(index_i+index_i_start)%matrix(mm,rank))
					block_o%ButterflyU(index_i+index_i_start)%matrix = block_i%ButterflyU(index_i)%matrix 					
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU(index_i+index_i_start)%matrix)/1024.0d3
					if(isnan(fnorm(block_o%ButterflyU(index_i+index_i_start)%matrix,mm,rank)))then
						write(*,*)'NAN in L 5'
					end if	
				endif
			enddo
		enddo
	enddo

else if(LR=='R')then


	do level=1, level_butterfly_loc
		do index_i=1, 2**(level-1)
			do index_j=1, 2**(level_butterfly_loc-level+1)
			! write(*,*)level,index_i,index_j
				index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)	
				if(level==level_butterfly_loc)then
				! write(*,*)'good 1'
					dimension_m = size(block_i%ButterflyU(2*index_i-1)%matrix,1)
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)
					! write(*,*)dimension_m,mm,rank,'d'
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)					
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix(dimension_m,rank))
					call gemm_omp(block_i%ButterflyU(2*index_i-1)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,&
					&block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,mm,rank)
! write(*,*)'good 1.1'

					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,rank)))then
						write(*,*)'NAN in R 1'
					end if

					dimension_m = size(block_i%ButterflyU(2*index_i)%matrix,1)
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)					
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix(dimension_m,rank))
					call gemm_omp(block_i%ButterflyU(2*index_i)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,&
					&block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,mm,rank)
! write(*,*)'good 2'
					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,rank)))then
						write(*,*)'NAN in R 2'
					end if
				else 
				! write(*,*)'good 3'
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix(mm,rank))
					block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix				

					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,mm,rank)))then
						write(*,*)'NAN in R 3'
					end if
					
	
					mm=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,1)
					rank=size(block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,2)
					deallocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix(mm,rank))
					block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix						
				! write(*,*)'good 4'
					if(isnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,mm,rank)))then
						write(*,*)'NAN in R 4'
					end if				
				end if
				
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix)/1024.0d3
				if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix)/1024.0d3

				if (level==1) then
					index_j_start = (ij_loc-1)*2**(level_butterfly_loc-level+1)	
					nn=size(block_i%ButterflyV(index_j)%matrix,1)
					rank=size(block_i%ButterflyV(index_j)%matrix,2)
					deallocate(block_o%ButterflyV(index_j+index_j_start)%matrix)
					allocate(block_o%ButterflyV(index_j+index_j_start)%matrix(nn,rank))
					block_o%ButterflyV(index_j+index_j_start)%matrix = block_i%ButterflyV(index_j)%matrix 					
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV(index_j+index_j_start)%matrix)/1024.0d3
					
					if(isnan(fnorm(block_o%ButterflyV(index_j+index_j_start)%matrix,nn,rank)))then
						write(*,*)'NAN in R 5'
					end if					
				endif					
			end do
		end do
	end do
					
end if	

end subroutine copy_randomizedbutterfly_partial


subroutine Initialize_Butterfly_SblockSmall(block,kover)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer rowblock
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j,kover
	type(matrixblock)::block
	
	
    allocate (butterfly_block_randomized(1))
    
    level_butterfly=block%level_butterfly
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    

    num_blocks=2**level_butterfly

    kk=0
    do j=1, num_blocks
        kk=max(kk,size(block%butterflyU(j)%matrix,2))
        kk=max(kk,size(block%butterflyV(j)%matrix,2))
	enddo
	

	! dimension_rank= block%rankmax !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	dimension_rank= block%rankmax *1.2d0**(kover)
	
	
    butterfly_block_randomized(1)%dimension_rank=dimension_rank
    !write (*,*) dimension_rank, int(real(kk)/num_blocks/2.0), mm
    
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))

	
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
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

		
		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyU(blocks)%matrix
		deallocate(matrixtemp1)

		
        dimension_n=size(block%ButterflyV(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_qr(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        ! allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%list(dimension_n))
		! butterfly_block_randomized(1)%ButterflyV(blocks)%list = block%ButterflyV(blocks)%list
		

		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyV(blocks)%matrix
		deallocate(matrixtemp1)

    enddo
	
	
	! ! write(*,*)'helloo'

	
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

end subroutine Initialize_Butterfly_SblockSmall

subroutine Initialize_Butterfly_Sblock(block,level_c,kover)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j
	type(matrixblock)::block
	
	
    ! ! block =>  cascading_factors(level_c)%matrices_block(rowblock) 
    allocate (butterfly_block_randomized(1))
    
    level_butterfly=int((maxlevel_for_blocks-block%level)/2)*2
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    

    num_blocks=2**level_butterfly

    kk=0
    do j=1, num_blocks
        kk=max(kk,size(block%butterflyU(j)%matrix,2))
        kk=max(kk,size(block%butterflyV(j)%matrix,2))
    enddo
	
	! write(*,*)'wocaca',block%level,block%row_group,block%col_group,block%rankmax
	! dimension_rank= max(ranktmp_glo,block%rankmax + kover) !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 

	dimension_rank= block%rankmax *1.2d0**(kover+1) !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	
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
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyU(blocks)%matrix
		deallocate(matrixtemp1)

		
        dimension_n=size(block%ButterflyV(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyV_qr(blocks)%matrix(dimension_n,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))
        ! allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%list(dimension_n))
		! butterfly_block_randomized(1)%ButterflyV(blocks)%list = block%ButterflyV(blocks)%list
		

		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyV(blocks)%matrix
		deallocate(matrixtemp1)

    enddo
	
	
	! write(*,*)'helloo'
    ! call invert_Butterfly_U()
    ! call invert_Butterfly_V()
	
	
	
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

end subroutine Initialize_Butterfly_Sblock



subroutine Initialize_Butterfly_Sblock_Empty(block,level_c,kover)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j
	type(matrixblock)::block
	
    ! ! block =>  cascading_factors(level_c)%matrices_block(rowblock) 
    allocate (butterfly_block_randomized(1))
    
    level_butterfly=int((maxlevel_for_blocks-block%level)/2)*2
    butterfly_block_randomized(1)%level_butterfly=level_butterfly
    

    num_blocks=2**level_butterfly

    kk=0
    do j=1, num_blocks
        kk=max(kk,size(block%butterflyU(j)%matrix,2))
        kk=max(kk,size(block%butterflyV(j)%matrix,2))
    enddo
	

	dimension_rank= block%rankmax + kover !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 
	! if(level_c==2)dimension_rank=9
	if(level_c==1)dimension_rank=max(dimension_rank,maxlevel_for_blocks)+kover

    groupm=block%row_group         ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    ! ! if (int(mm/num_blocks)<dimension_rank) then
        ! ! dimension_rank=int(mm/num_blocks)
    ! ! endif
    butterfly_block_randomized(1)%dimension_rank=dimension_rank
    !write (*,*) dimension_rank, int(real(kk)/num_blocks/2.0), mm
    
    allocate (butterfly_block_randomized(1)%ButterflyU(2**level_butterfly))
    allocate (butterfly_block_randomized(1)%ButterflyV(2**level_butterfly))

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
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))
		
		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_m
				butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyU(blocks)%matrix
		deallocate(matrixtemp1)

		
        dimension_n=size(block%ButterflyV(blocks)%matrix,1)
        allocate (butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
        do j=1, dimension_rank
            do i=1, dimension_n
				butterfly_block_randomized(1)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		! butterfly_block_randomized(1)%ButterflyV_old(blocks)%matrix=butterfly_block_randomized(1)%ButterflyV(blocks)%matrix
		deallocate(matrixtemp1)

    enddo
	
    if (level_butterfly/=0) then
        allocate (butterfly_block_randomized(1)%ButterflyKerl(level_butterfly))

        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_row=num_row
            butterfly_block_randomized(1)%ButterflyKerl(level)%num_col=num_col
            allocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(num_row,num_col))
        enddo
    endif	
	

    return

end subroutine Initialize_Butterfly_Sblock_Empty

subroutine Get_Randomized_Vectors_Sblock(level_c,rowblock)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
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
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(6))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	! if(mod(level_butterfly,2)==0)then
		! Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)
    ! else
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)
	! end if
	Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	num_vectors = Nsub*dimension_rank
	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vectors))
	allocate (RandomVectors_InOutput(6)%vector(nn,num_vectors))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vectors))
    allocate (RandomVectors_InOutput(3)%vector(mm,num_vectors))    
	allocate (RandomVectors_InOutput(4)%vector(mm,num_vectors))
    allocate (RandomVectors_InOutput(5)%vector(mm,num_vectors))
	do ii =1,6
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
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
				 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo 	


	
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
				 RandomVectors_InOutput(4)%vector(ii+k,jj)=random_complex_number()
			 enddo
		 enddo
		 !$omp end parallel do
		 if(mod(i,Ng)==0)idx_start = idx_start + dimension_rank			
	enddo	 
	 

	! get the right multiplied vectors
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1	
	allocate(vec_old(mm,num_vectors))
	allocate(vec_new(mm,num_vectors))
	vec_old = RandomVectors_InOutput(2)%vector

	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			! write(*,*)level,ii
			groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   	

			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			
			! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
			
			! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
			else 
				call butterfly_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
			end if
		end do
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	RandomVectors_InOutput(3)%vector = vec_new
	deallocate(vec_old)
	deallocate(vec_new)
	

	
	! get the left multiplied vectors	
	allocate(vec_old(mm,num_vectors))
	allocate(vec_new(mm,num_vectors))	
	vec_old = RandomVectors_InOutput(4)%vector
	do level = level_c+1,Maxlevel_for_blocks+1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'T',idx_end_loc-idx_start_loc+1,num_vectors,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
			else 				
				call butterfly_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vectors,vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
			end if
		end do
		vec_old = vec_new
	end do	
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(4)%vector(1,1)
	RandomVectors_InOutput(5)%vector = vec_new
	deallocate(vec_old)
	deallocate(vec_new)	
	
    random1=>RandomVectors_InOutput(5)
    random2=>RandomVectors_InOutput(6)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	          
    return                

end subroutine Get_Randomized_Vectors_Sblock





subroutine Get_Randomized_Vectors_LL_Sblock_exact(level_c,rowblock,nth_s,nth_e,num_vect_sub)

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
	
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)
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
	
	do nth= nth_s,nth_e
		do i=1, num_blocks
			if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_m=basis_group(groupm_start+i-1)%head
				tailer_m=basis_group(groupm_start+i-1)%tail
				mm=tailer_m-header_m+1
				k=header_m-header_mm	


				!$omp parallel do default(shared) private(ii,jj)
				 do jj=1,dimension_rank
					 do ii=1, mm
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*dimension_rank+jj)=random_complex_number()	
					 enddo
				 enddo
				 !$omp end parallel do
			 end if
		end do
	end do
	
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	

	
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

end subroutine Get_Randomized_Vectors_LL_Sblock_exact





subroutine Get_Randomized_Vectors_LL_SblockSmall(level_c,ii_inv,block_o,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,mm_loc,nn_loc
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	integer level_c,ii_inv
	
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
				allocate(matrixtemp1(num_vect_subsub,mm_loc))
				call RandomMat(num_vect_subsub,mm_loc,min(mm_loc,num_vect_subsub),matrixtemp1,0)
				
				! !$omp parallel do default(shared) private(ii,jj)
				 do jj=1,num_vect_subsub
					 do ii=1, mm_loc
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
					 enddo
				 enddo
				 ! !$omp end parallel do
				 deallocate(matrixtemp1)
			 end if
			 k=k+mm_loc
		end do
	end do
	
	n1 = OMP_get_wtime()
	call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'T',mm,num_vect_sub,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector)
	
    random1=>RandomVectors_InOutput(2)
    random2=>RandomVectors_InOutput(3)
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1		
	
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

end subroutine Get_Randomized_Vectors_LL_SblockSmall

subroutine Get_Randomized_Vectors_LL_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	
    ctemp1=1.0d0 ; ctemp2=0.0d0	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 

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
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
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
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))	
	vec_old = RandomVectors_InOutput(1)%vector
	do level = level_c+1,Maxlevel_for_blocks+1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		
		if(level==Maxlevel_for_blocks+1)then
			n1 = OMP_get_wtime()   !!!! comment: will the omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'T',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
			end do
			! !$omp end parallel do
			n2 = OMP_get_wtime()
			! time_tmp = time_tmp + n2 - n1
		else 
			n1 = OMP_get_wtime()
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
				call butterfly_block_MVP_inverse_dat(level,ii,'T',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))			
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
	
	
    random1=>RandomVectors_InOutput(2)
    random2=>RandomVectors_InOutput(3)
	
	n1 = OMP_get_wtime()
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
	n2 = OMP_get_wtime()
	! time_tmp = time_tmp + n2 - n1			
	
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

end subroutine Get_Randomized_Vectors_LL_Sblock


subroutine Get_Randomized_Vectors_RR_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)

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
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
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
						 ! RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
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
			n1 = OMP_get_wtime() !!!! comment: will the omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
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
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call butterfly_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
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

end subroutine Get_Randomized_Vectors_RR_Sblock
subroutine Butterfly_Sblock_LowRank(level_c,rowblock)

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
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
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
			n1 = OMP_get_wtime() !!!! comment: will the omp cause segment fault?
			! !$omp parallel do default(shared) private(ii,groupm_diag,idx_start_loc,idx_end_loc)
			do ii = idx_start_diag,idx_start_diag+N_diag-1
				! write(*,*)level,ii
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
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
				groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call butterfly_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
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

end subroutine Butterfly_Sblock_LowRank




subroutine Get_Randomized_Vectors_RR_Sblock_exact(level_c,rowblock,nth_s,nth_e,num_vect_sub)

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
	
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)    
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)
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
		do i=1, num_blocks
			if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_n=basis_group(groupn_start+i-1)%head
				tailer_n=basis_group(groupn_start+i-1)%tail
				nn=tailer_n-header_n+1
				k=header_n-header_nn

				!$omp parallel do default(shared) private(ii,jj)
				 do jj=1,dimension_rank
					 do ii=1, nn
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*dimension_rank+jj)=random_complex_number()	
					 enddo
				 enddo
				 !$omp end parallel do
			 end if
		end do
	end do
	
	! get the right multiplied vectors
	idx_start_glo = basis_group(groupm)%head
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	
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

end subroutine Get_Randomized_Vectors_RR_Sblock_exact


subroutine Get_Randomized_Vectors_RR_Test_Sblock(level_c,rowblock,num_vect_sub)

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
	
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
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
			groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			
			! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
			
			! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub)			
			else 
				call butterfly_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub))
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

end subroutine Get_Randomized_Vectors_RR_Test_Sblock


subroutine Get_Randomized_Vectors_RR_Test_Sblock_exact(level_c,rowblock,num_vect_sub)

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
	
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
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

		!$omp parallel do default(shared) private(ii,jj)
		 do jj=1,num_vect_sub
			 do ii=1, nn
				 RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	
			 enddo
		 enddo
		 !$omp end parallel do
	end do
	
	! get the right multiplied vectors
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	

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

end subroutine Get_Randomized_Vectors_RR_Test_Sblock_exact




subroutine Get_Randomized_Vectors_RR_SblockSmall(level_c,ii_inv,block_o,nth_s,nth_e,num_vect_sub,unique_nth)


    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,ii_inv
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag,mm_loc,nn_loc
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	
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

				allocate(matrixtemp1(num_vect_subsub,nn_loc))
				call RandomMat(num_vect_subsub,nn_loc,min(nn_loc,num_vect_subsub),matrixtemp1,0)
				
				! !$omp parallel do default(shared) private(ii,jj)
				 do jj=1,num_vect_subsub
					 do ii=1, nn_loc
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
					 enddo
				 enddo
				 ! !$omp end parallel do
				 
				 deallocate(matrixtemp1)
			 end if
			 k=k+nn_loc
		end do		
	end do
	
	n1 = OMP_get_wtime() 
	ctemp1=1.0d0 ; ctemp2=0.0d0
	random1=>RandomVectors_InOutput(1)
	random2=>RandomVectors_InOutput(2)
	call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
	call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
	n2 = OMP_get_wtime()	
	! time_tmp = time_tmp + n2 - n1		
	
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

end subroutine Get_Randomized_Vectors_RR_SblockSmall



subroutine Get_Randomized_Vectors_RR_Test_SblockSmall(level_c,ii_inv,block_o,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,ii_inv
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,mm_loc,nn_loc
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	
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
	call butterfly_block_MVP_inverse_dat(level_c,ii_inv,'N',mm,num_vect_sub,RandomVectors_InOutput(2)%vector,RandomVectors_InOutput(3)%vector)
	

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

end subroutine Get_Randomized_Vectors_RR_Test_SblockSmall





subroutine Test_Randomized_Sblock(level_c,rowblock,block_old,error)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	type(matrixblock)::block_old
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	real*8::error
	
    
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(6))


    num_vectors=1
    
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vectors))
	allocate (RandomVectors_InOutput(6)%vector(nn,num_vectors))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
    allocate (RandomVectors_InOutput(2)%vector(mm,num_vectors))
    allocate (RandomVectors_InOutput(3)%vector(mm,num_vectors))    
	allocate (RandomVectors_InOutput(4)%vector(mm,num_vectors))
    allocate (RandomVectors_InOutput(5)%vector(mm,num_vectors))
	do ii =1,6
		RandomVectors_InOutput(ii)%vector = 0
	end do 
	 
	! write(*,*)groupn,basis_group(groupn)%head,groupm,basis_group(groupm)%head
	 
	 
	idx_start_glo = basis_group(groupm)%head	 
	   
    ! !$omp parallel do default(shared) private(ii,jj,a,b,c,d)
    do jj=1, num_vectors
        do ii=1, nn
            ! call random_number(a)
            ! call random_number(b)
            ! call random_number(c)
            ! call random_number(d)
            ! if (c<0.5d0) then
                ! a=-a
            ! endif
            ! if (d<0.5d0) then
                ! b=-b
            ! endif
            RandomVectors_InOutput(1)%vector(ii,jj)=random_complex_number()
        enddo
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
    enddo
    ! !$omp end parallel do
    
	
	
	! get the right multiplied vectors
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_old,'N',random1,random2,ctemp1,ctemp2)
	
	allocate(vec_old(mm,num_vectors))
	allocate(vec_new(mm,num_vectors))
	vec_old = RandomVectors_InOutput(2)%vector

	! write(*,*)'begin'
	
	do level = Maxlevel_for_blocks+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = (rowblock-1)*N_diag+1
		vec_new = 0
		do ii = idx_start_diag,idx_start_diag+N_diag-1
			! write(*,*)level,ii
			groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

			
			idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			
			! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
			
			! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
			else 
				call butterfly_block_MVP_inverse_dat(level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
				! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
			end if
		end do
		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	RandomVectors_InOutput(3)%vector = vec_new
	deallocate(vec_old)
	deallocate(vec_new)
	
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(4)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)	

	error = sqrt(sum(abs(RandomVectors_InOutput(3)%vector(:,1)-RandomVectors_InOutput(4)%vector(:,1))**2)/sum(abs(RandomVectors_InOutput(3)%vector(:,1))**2))
	! ! write(*,*)error
	
	
	! ! ! ! ! ! ! ! write(*,*)'end'
	
	! ! ! ! get the left multiplied vectors	
	! ! ! allocate(vec_old(mm,num_vectors))
	! ! ! allocate(vec_new(mm,num_vectors))	
	! ! ! vec_old = RandomVectors_InOutput(4)%vector
	! ! ! do level = level_c+1,Maxlevel_for_blocks+1
		! ! ! N_diag = 2**(level-level_c-1)
		! ! ! idx_start_diag = (rowblock-1)*N_diag+1
		! ! ! do ii = idx_start_diag,idx_start_diag+N_diag-1
			! ! ! groupm_diag = cascading_factors(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   
			! ! ! idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
			! ! ! idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
			! ! ! if(level==Maxlevel_for_blocks+1)then
				! ! ! call fullmat_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'T',idx_end_loc-idx_start_loc+1,num_vectors,&
				! ! ! &vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! ! ! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors)
			! ! ! else 				
				! ! ! call butterfly_block_MVP_randomized_dat(cascading_factors(level)%matrices_block_inverse(ii),'T',idx_end_loc-idx_start_loc+1,idx_end_loc-idx_start_loc+1,num_vectors,&
				! ! ! &vec_old(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = vec_old(idx_start_loc:idx_end_loc,1:num_vectors) + vec_new(idx_start_loc:idx_end_loc,1:num_vectors)
				! ! ! ! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
			! ! ! end if
		! ! ! end do
		! ! ! vec_old = vec_new
	! ! ! end do	
	! ! ! ! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(4)%vector(1,1)
	! ! ! RandomVectors_InOutput(5)%vector = vec_new
	! ! ! deallocate(vec_old)
	! ! ! deallocate(vec_new)	
	
    ! ! ! random1=>RandomVectors_InOutput(5)
    ! ! ! random2=>RandomVectors_InOutput(6)
    ! ! ! ctemp1=1.0 ; ctemp2=0.0
    ! ! ! call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	          
   


    ! !$omp parallel do default(shared) private(i)
    do i=1, 6
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)



   return                

end subroutine Test_Randomized_Sblock




subroutine Reconstruction_LL_Sblock(level_c,rowblock)
    
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

    block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
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
			call Get_Randomized_Vectors_LL_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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
    
end subroutine Reconstruction_LL_Sblock


subroutine Reconstruction_LL_SblockSmall(level_c,ii_inv,block_o)
    
    use MODULE_FILE
    implicit none
	
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
    integer level_c,ii_inv
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock)::block_o
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
			call Get_Randomized_Vectors_LL_SblockSmall(level_c,ii_inv,block_o,nth_s,nth_e,num_vect_sub,unique_nth)
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
    
end subroutine Reconstruction_LL_SblockSmall



subroutine Reconstruction_RR_Sblock(level_c,rowblock,error)
    
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
		
    block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
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
			call Get_Randomized_Vectors_RR_Sblock(level_c,rowblock,nth_s,nth_e,num_vect_sub,unique_nth)
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

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_Sblock(level_c,rowblock,error)


	
    return
    
end subroutine Reconstruction_RR_Sblock




subroutine Reconstruction_RR_SblockSmall(level_c,ii_inv,block_o,error)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,ii_inv    
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
    type(matrixblock)::block_o
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
			call Get_Randomized_Vectors_RR_SblockSmall(level_c,ii_inv,block_o,nth_s,nth_e,num_vect_sub,unique_nth)
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

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)

	deallocate(Random_Block)

	call Test_Error_RR_SblockSmall(level_c,ii_inv,block_o,error)


	
    return
    
end subroutine Reconstruction_RR_SblockSmall


subroutine Test_Error_RR_Sblock_exact(level_c,rowblock,error)

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
	
	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
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

	call Get_Randomized_Vectors_RR_Test_Sblock_exact(level_c,rowblock,num_vect)

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

end subroutine Test_Error_RR_Sblock_exact


subroutine Test_Error_RR_Sblock(level_c,rowblock,error)

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
	

	
	block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
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

	call Get_Randomized_Vectors_RR_Test_Sblock(level_c,rowblock,num_vect)

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

end subroutine Test_Error_RR_Sblock





subroutine Test_Error_RR_SblockSmall(level_c,ii_inv,block_o,error)

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
	integer level_c,rowblock,dimension_m,ii_inv 
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock)::block_o
	
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

	call Get_Randomized_Vectors_RR_Test_SblockSmall(level_c,ii_inv,block_o,num_vect)

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

end subroutine Test_Error_RR_SblockSmall


subroutine Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer nth_s,nth_e,unique_nth
   integer num_vect_sub,Ng

	if(adaptive_flag==1)then
		call Resolving_Butterfly_LL_new_adaptive(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
	else 
		call Resolving_Butterfly_LL_new_uniform(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
	end if
   
end subroutine Resolving_Butterfly_LL_new
   
   
subroutine Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   integer nth_s,nth_e,unique_nth
   integer num_vect_sub,Ng

	if(adaptive_flag==1)then
		call Resolving_Butterfly_RR_new_adaptive(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
	else 
		call Resolving_Butterfly_RR_new_uniform(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
	end if
   
end subroutine Resolving_Butterfly_RR_new   
   








subroutine Resolving_Butterfly_LL_new_adaptive(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer nth_s,nth_e,unique_nth
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,mm,kk,level_left,level_right, rs,re,rank,level_right_start,level_left_start
   integer index_i, index_j, iter, vector1, vector2, direction, round, flag
   real*8 a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   integer kmax
   
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real*8, allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,noe,Ng,dimension_nn,nn1,nn2,ieo,level_butterfly
   real*8::n1,n2
   type(butterflyblock_randomized), pointer :: blocks
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
	
   if (level_butterfly/=0) then
       mm=num_vect_subsub !rank

	   
	   level_right_start = floor_safe(level_butterfly/2d0)	!  check here later		   
	   ! ! level_right_start = level_butterfly+1
	   
	   
	   do level_right=0,unique_nth !level_right_start
		   ! kmax = ceiling_safe(rank/dble(2**(level_right_start-level_right)))+1
		   ! if(level_butterfly==9)write(*,*)level_right,kmax
			! write(*,*)level_right,'haha'
		   if (level_right==0) then 
			   do nth = nth_s,nth_e
				   !$omp parallel do default(shared) private(j)
				   do j=1, num_blocks
						call OneV_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s)	   
					end do
					!$omp end parallel do
				end do
	
		   elseif (level_right==level_butterfly+1) then
				write(*,*)'the right half scheme should not touch leftmost matrix'
				stop
		   else

			   num_row=blocks%ButterflyKerl(level_right)%num_row
			   num_col=blocks%ButterflyKerl(level_right)%num_col
			   
			   do nth = nth_s,nth_e
				   noe = ceiling_safe(nth*Ng*2d0/num_col)
				   index_i = ceiling_safe(noe/2d0)
				   !$omp parallel do default(shared) private(j,index_j)
				   do j=1, num_col, 2
						index_j=int((j+1)/2)
						call OneKernel_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s)
				   enddo
				   !$omp end parallel do
			   enddo
		   endif	   
	   end do
	   
   endif
   
   
   return

end subroutine Resolving_Butterfly_LL_new_adaptive

subroutine OneV_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,level_butterfly
	
   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
   if(level_right==unique_nth)then
	   dimension_nn=size(blocks%butterflyV(j)%matrix,1)
	   allocate(matB(mm,dimension_nn))
	   call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix(1:dimension_nn,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_nn,mm)
	   call GetRank(mm,dimension_nn,matB,rank,Rank_detection_factor)
	   if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
							   
	   if(allocated(blocks%butterflyV(j)%matrix))deallocate(blocks%butterflyV(j)%matrix)
	   ! if(allocated(blocks%ButterflyVInv(j)%matrix))deallocate(blocks%ButterflyVInv(j)%matrix)
	   if(allocated(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix))deallocate(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix)
	   allocate(blocks%butterflyV(j)%matrix(dimension_nn,rank))
	   ! allocate(blocks%ButterflyVInv(j)%matrix(rank,dimension_nn))
	   allocate(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(rank,num_vect_sub))
	   ! call RandomMat(rank,dimension_nn,min(rank,dimension_nn),blocks%ButterflyVInv(j)%matrix,0)
	   
	   allocate(matC(rank,dimension_nn),matA(mm,rank),matinv(dimension_nn,rank))
	   ! call copymatT_omp(blocks%ButterflyVInv(j)%matrix,matinv,rank,dimension_nn)
	   call copymatT_omp(blocks%KerInv(1:rank,1:dimension_nn),matinv,rank,dimension_nn)																		   
	   if(isnan(fnorm(matB,mm,dimension_nn)) .or. isnan(fnorm(matinv,dimension_nn,rank)))then
		write(*,*)shape(matB),fnorm(matB,mm,dimension_nn),fnorm(matinv,dimension_nn,rank),j,'hei',fnorm(random_Block(1)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix,dimension_nn,mm)
		stop
	   end if
	   call gemm_omp(matB,matinv,matA,mm,dimension_nn,rank)
	   call LeastSquare(mm,rank,dimension_nn,matA,matB,matC,LS_tolerance)
	   call copymatT_omp(matC,blocks%ButterflyV(j)%matrix,rank,dimension_nn)						   
	   deallocate(matB,matC,matA,matinv)						   
   else 
	   rank=size(blocks%butterflyV(j)%matrix,2)
	   dimension_nn=size(blocks%butterflyV(j)%matrix,1)									
	   allocate(matB(mm,dimension_nn),matA(mm,rank),matinv(dimension_nn,rank))
	   call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix(1:dimension_nn,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_nn,mm)
	   ! call copymatT_omp(blocks%ButterflyVInv(j)%matrix,matinv,rank,dimension_nn)
	   call copymatT_omp(blocks%KerInv(1:rank,1:dimension_nn),matinv,rank,dimension_nn)																					   
	   if(isnan(fnorm(matB,mm,dimension_nn)) .or. isnan(fnorm(matinv,dimension_nn,rank)))then
		write(*,*)fnorm(matB,mm,dimension_nn),fnorm(matinv,dimension_nn,rank),j,'hei1'
		stop
	   end if
	   call gemm_omp(matB,matinv,matA,mm,dimension_nn,rank)
	   if(.not. allocated(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix))allocate(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(rank,num_vect_sub))
	   call copymatT_omp(matA,random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)					   
	   deallocate(matB,matA,matinv)	
   end if   
   
end subroutine OneV_LL 	


subroutine OneKernel_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer index_i,index_j,i,j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,ieo,noe,rs,re,level_butterfly

   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
   i = index_i*2-1
   j = index_j*2-1
   ieo = i + 1 - mod(noe,2)

	nn1 = size(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix,1)
	nn2 = size(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix,1)

	if(level_right==unique_nth)then
		allocate (matB(mm,nn1+nn2))
		call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)
		if(mod(noe,2)==1)then
			call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
			
			rs = 1
			re = rank
		else 
			call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
									   
			rs = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)+1
			re = rs+rank-1		
		end if


		if(allocated(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix))deallocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix)
		if(allocated(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix))deallocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix)
		allocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix(rank,nn1))
		allocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix(rank,nn2))
		if(allocated(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix))deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix)
		allocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix(rank,mm))
		

		allocate (matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatN_omp(blocks%KerInv(1:nn1+nn2,rs:re),matinv,nn1+nn2,rank)																  
		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'ho'
		 stop
	    end if
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
		call copymatN_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix,rank,nn1)
		call copymatN_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix,rank,nn2)							
		deallocate(matB,matC,matA,matinv)
		deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix)
		deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix)
		deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix)																									
	else 
		if(mod(noe,2)==1)then
			rank = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)
			rs = 1
			re = rank
		else 
			rank = size(blocks%ButterflyKerl(level_right)%blocks(i+1,j)%matrix,1)
			rs = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)+1
			re = rs+rank-1
		end if
		allocate (matB(mm,nn1+nn2),matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)
										

		call copymatN_omp(blocks%KerInv(1:nn1+nn2,rs:re),matinv,nn1+nn2,rank)																		  
		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'ho1'
		 stop
	    end if		
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		if(.not. allocated(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix))allocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix(rank,num_vect_sub))
		call copymatT_omp(matA,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matC,matA,matinv)	
		deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix)
		deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix)
	end if
		! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_right,nth,i,j,error0,'L' 


   
end subroutine OneKernel_LL   



subroutine Resolving_Butterfly_RR_new_adaptive(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank,level_right_start,level_left_start,nn1,nn2,rs,re
   integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, flag
   real*8 a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   
   type(butterflyblock_randomized), pointer :: blocks
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real*8, allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,ind_c,noe,Ng,nth_s,nth_e,dimension_mm,dimension_n,jeo,level_butterfly
   real*8::n1,n2
   integer::kmax,unique_nth
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 

   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
   
   if (level_butterfly/=0) then
       mm=num_vect_subsub !rank

	   level_left_start= floor_safe(level_butterfly/2d0)+1    !  check here later		   
	   ! ! level_left_start = 0
	   
	   random=>random_Block(1)
	   if(level_left_start>0 .and. level_left_start==unique_nth)then
			n1 = OMP_get_wtime()
			call Butterfly_partial_MVP_Half('N',0,level_left_start-1,random,num_vect_sub,nth_s,nth_e,Ng)
			! call Butterfly_partial_MVP(blocks,'N',0,level_left_start-1,random,num_vect_sub)
			n2 = OMP_get_wtime()
			time_halfbuttermul = time_halfbuttermul + n2-n1		 
		endif 
	   
	   
	   do level_left = level_butterfly+1,unique_nth,-1
			if (level_left==level_butterfly+1) then
				do nth=nth_s,nth_e
					!$omp parallel do default(shared) private(i)
					do i=1, num_blocks   
						call OneU_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s)
					end do
					!$omp end parallel do					
				end do
			elseif (level_left==0) then
				write(*,*)'the left half scheme should not touch rightmost matrix'
				stop
			else 
			   num_row=blocks%ButterflyKerl(level_left)%num_row
			   num_col=blocks%ButterflyKerl(level_left)%num_col
			   do nth = nth_s,nth_e
				   noe = ceiling_safe(nth*Ng*2d0/num_row)
				   index_j = ceiling_safe(noe/2d0)
				   !$omp parallel do default(shared) private(i,index_i)
				   do i=1, num_row, 2
					   index_i=int((i+1)/2)
					   call OneKernel_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s)									
				   enddo
				   !$omp end parallel do	
			   enddo	   
			end if
	   end do
   endif
   
   
   return

end subroutine Resolving_Butterfly_RR_new_adaptive


subroutine OneU_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer i,level_left,unique_nth,dimension_mm,mm,rank,num_vect_sub,nth,nth_s,level_butterfly
	
   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
	if(level_left==unique_nth)then
		dimension_mm=size(blocks%butterflyU(i)%matrix,1)	
		allocate(matB(mm,dimension_mm))
		call copymatT_omp(random_Block(1)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)							
		call GetRank(mm,dimension_mm,matB,rank,Rank_detection_factor)
		if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
		
		if(allocated(blocks%butterflyU(i)%matrix))deallocate(blocks%butterflyU(i)%matrix)
		if(allocated(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix))deallocate(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix)
		allocate(blocks%butterflyU(i)%matrix(dimension_mm,rank))
		allocate(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(rank,num_vect_sub))
		allocate(matC(rank,dimension_mm),matA(mm,rank),matinv(dimension_mm,rank))	
		call copymatT_omp(blocks%KerInv(1:rank,1:dimension_mm),matinv,rank,dimension_mm)																					 
		if(isnan(fnorm(matB,mm,dimension_mm)) .or. isnan(fnorm(matinv,dimension_mm,rank)))then
		 write(*,*)fnorm(matB,mm,dimension_mm),fnorm(matinv,dimension_mm,rank),i,'heee'
		 stop
	    end if		
		call gemm_omp(matB,matinv,matA,mm,dimension_mm,rank)							
		call LeastSquare(mm,rank,dimension_mm,matA,matB,matC,LS_tolerance)
		call copymatT_omp(matC,blocks%ButterflyU(i)%matrix,rank,dimension_mm)							
		deallocate(matB,matC,matA,matinv)
	else 
		dimension_mm=size(blocks%butterflyU(i)%matrix,1)						
		rank=size(blocks%butterflyU(i)%matrix,2)						
		allocate(matB(mm,dimension_mm),matA(mm,rank),matinv(dimension_mm,rank))	
		call copymatT_omp(random_Block(1)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)
		call copymatT_omp(blocks%KerInv(1:rank,1:dimension_mm),matinv,rank,dimension_mm)																					 
		if(isnan(fnorm(matB,mm,dimension_mm)) .or. isnan(fnorm(matinv,dimension_mm,rank)))then
		 write(*,*)fnorm(matB,mm,dimension_mm),fnorm(matinv,dimension_mm,rank),i,'heee1'
		 stop
	    end if			
		call gemm_omp(matB,matinv,matA,mm,dimension_mm,rank)
		if(.not. allocated(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix))allocate(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(rank,num_vect_sub))
		call copymatT_OMP(matA,random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matA,matinv)					
	end if	   
   

end subroutine OneU_RR  



subroutine OneKernel_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer index_i,index_j,i,j,level_left,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,jeo,noe,rs,re,level_left_start,level_butterfly

   blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
	
	i = index_i*2-1
	j = index_j*2-1
	jeo = j + 1 - mod(noe,2)					   
	
	nn1 = size(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix,1)
	nn2 = size(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix,1)

	if(level_left==unique_nth)then
		if(level_left==level_left_start)then
			allocate (matB(mm,nn1+nn2))
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
			
			if(mod(noe,2)==1)then
				rank = size(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix,1)
				rs = 1
				re = rank
			else 
				rank = size(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix,1)
				rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
				re = rs+rank-1		
			end if																		

			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix)
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix)
			allocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix(nn2,rank))
			allocate(matC(rank,nn1+nn2),matA(mm,rank))
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
			call copymatT_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix,rank,nn1)
			call copymatT_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix,rank,nn2)										
			deallocate(matB,matC,matA)

		else 
			allocate (matB(mm,nn1+nn2))
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
			if(mod(noe,2)==1)then
				call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
				if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
				
				rs = 1
				re = rank
			else 
				call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
				if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
				
				rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
				re = rs+rank-1		
			end if																		
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix)
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix)
			allocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix(nn2,rank))
			if(allocated(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix))deallocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix)
			allocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(rank,mm))
			
			allocate(matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
			call copymatT_OMP(blocks%KerInv(rs:re,1:nn1+nn2),matinv,rank,nn1+nn2)															   

			if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
			 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'helo'
			 stop
			end if	
			call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
			
			call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
			call copymatT_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix,rank,nn1)
			call copymatT_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix,rank,nn2)										
			deallocate(matB,matC,matA,matinv)
			
		end if
		deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix)		
		deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix)		
		deallocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix)																			   
	else 
		if(mod(noe,2)==1)then
			rank = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)
			rs = 1
			re = rank
		else 
			rank = size(blocks%ButterflyKerl(level_left)%blocks(i,j+1)%matrix,2)
			rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
			re = rs+rank-1
		end if							
		allocate (matB(mm,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
		
		call copymatT_omp(blocks%KerInv(rs:re,1:nn1+nn2),matinv,rank,nn1+nn2)																		  

		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'helo'
		 stop
		end if			
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		if(.not. allocated(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix))allocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(rank,num_vect_sub))
		call copymatT_omp(matA,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matA,matinv)
		deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix)		
		deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix)				
	end if

	! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_left,nth,i,j,error0,'R'
end subroutine OneKernel_RR  

 

subroutine Resolving_Butterfly_LL_new_uniform(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer nth_s,nth_e,unique_nth
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank,level_right_start,level_left_start
   integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, flag
   real(kind=8) a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   integer kmax
   
   ! type(matricesblock), pointer :: blocks
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real(kind=8), allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer num_vect_sub,num_vect_subsub,nth,nulldim,ind_r,noe,Ng,level_butterfly,dimension_n
   real*8::n1,n2
   
   
   norm3R=0. ; norm4R=0.
   level_butterfly=butterfly_block_randomized(1)%level_butterfly 
   num_blocks=2**level_butterfly
   rank=butterfly_block_randomized(1)%dimension_rank
   ! write(*,*)rank,'how come'
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
	
   if (level_butterfly/=0) then
       mm=num_vect_subsub !rank
       nn=2*rank
       error=1000.
	   iter = 1	
	   
	   level_right_start = floor_safe(level_butterfly/2d0)	!  check here later		   
	   ! ! level_right_start = level_butterfly+1
	   
	   random=>random_Block(1)
	   ! ! if(level_right_start<level_butterfly+1 .and. unique_nth==level_right_start)then
			! ! n1 = OMP_get_wtime()
			! ! call Butterfly_partial_MVP_Half('T',0,level_butterfly-level_right_start,random,num_vect_sub,nth_s,nth_e,Ng)
			! ! ! call Butterfly_partial_MVP('T',0,level_butterfly-level_right_start,random)
			! ! n2 = OMP_get_wtime()
			! ! time_halfbuttermul = time_halfbuttermul + n2-n1		
	   ! ! endif 	  
	   
	   
	   do level_right=0,unique_nth !level_right_start
		   ! kmax = ceiling_safe(rank/dble(2**(level_right_start-level_right)))+1
		   ! if(level_butterfly==9)write(*,*)level_right,kmax
		   kmax = rank
			! write(*,*)level_right,'haha'
		   if (level_right==0) then 
			   do nth = nth_s,nth_e
				   do j=1, num_blocks
					   dimension_n=size(butterfly_block_randomized(1)%butterflyV(j)%matrix,1)
					   allocate(matB(mm,dimension_n),matC(rank,dimension_n),matA(mm,rank),matinv(dimension_n,rank))
					   call copymatT_OMP(random_Block(1)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix(1:dimension_n,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_n,mm)
					   ! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyVInv(j)%matrix,matinv,rank,dimension_n)
					   ! write(*,*)'haa'
					   call copymatT_OMP(butterfly_block_randomized(1)%KerInv(1:rank,1:dimension_n),matinv,rank,dimension_n)
					   ! write(*,*)'hab'
					   call gemm_omp(matB,matinv,matA,mm,dimension_n,rank)
					   if(.not. allocated(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix))allocate(random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(rank,mm))
					   call copymatT_OMP(matA,random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
					   
					   if(level_right==unique_nth)then
						   nulldim = random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%nulldim
						   ! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyV(j)%matrix,matC,dimension_n,rank)
							! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%null,mm,rank,&
						    ! &dimension_n,random_Block(1)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%nulldim,up_tolerance,LS_tolerance,error0,kmax)
							call LeastSquare(mm,rank,dimension_n,matA,matB,matC,LS_tolerance)
							if(.not. allocated(butterfly_block_randomized(1)%ButterflyV(j)%matrix))allocate(butterfly_block_randomized(1)%ButterflyV(j)%matrix(dimension_n,rank))
							call copymatT_OMP(matC,butterfly_block_randomized(1)%ButterflyV(j)%matrix,rank,dimension_n)
							! write(*,'(I5,I5,I5,Es16.7E3,A2)')unique_nth,nth,j,error0,'L'						   
					    end if
						deallocate(matB,matC,matA,matinv)
					end do
				
				end do
				
				
				
		   elseif (level_right==level_butterfly+1) then
				write(*,*)'the right half scheme should not touch leftmost matrix'
				stop
		   else
			   allocate (matB(mm,rank*2),matC(rank,rank*2),matA(mm,rank),matinv(rank*2,rank))
			   
			   num_row=butterfly_block_randomized(1)%ButterflyKerl(level_right)%num_row
			   num_col=butterfly_block_randomized(1)%ButterflyKerl(level_right)%num_col
			   
			   do nth = nth_s,nth_e
			   
				   noe = ceiling_safe(nth*Ng*2d0/num_col)
				   ind_r = ceiling_safe(noe/2d0)
				   
				   do j=1, num_col, 2
					   index_j=int((j+1)/2)
					   do i=1, num_row, 2
						   index_i=int((i+1)/2)
						   if(index_i==ind_r)then
								
																																													 
																																															  
		
								call copymatT_OMP(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:rank),rank,mm)
								call copymatT_OMP(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+rank:rank*2),rank,mm)
								deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix)
								deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix)

								if(mod(noe,2)==1)then
									! ! if(level_right==level_right_start)then
										! ! call copymatT_OMP(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
									! ! else 
										call copymatN_omp(butterfly_block_randomized(1)%KerInv(1:2*rank,1:rank),matinv,2*rank,rank)

										call gemm_omp(matB,matinv,matA,mm,2*rank,rank)
										if(.not. allocated(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%matrix))allocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%matrix(rank,mm))
										call copymatT_OMP(matA,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
									! ! end if
		 
									if(level_right==unique_nth)then
										deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%matrix)
										! call copymatN_omp(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j)%matrix,matC(1:rank,1:rank),rank,rank)
										! call copymatN_omp(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j+1)%matrix,matC(1:rank,rank+1:2*rank),rank,rank)	
										! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%null,mm,rank,&
										! &2*rank,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i,index_j)%nulldim,up_tolerance,LS_tolerance,error0,kmax)									
										call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
										if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j)%matrix(rank,rank))
										if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j+1)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j+1)%matrix(rank,rank))
										call copymatN_omp(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j)%matrix,rank,rank)
										call copymatN_omp(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i,j+1)%matrix,rank,rank)
									end if
								else
									! ! if(level_right==level_right_start)then
										! ! call copymatT_OMP(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
									! ! else
										call copymatN_omp(butterfly_block_randomized(1)%KerInv(1:2*rank,1+rank:2*rank),matinv,2*rank,rank)								
									
										call gemm_omp(matB,matinv,matA,mm,2*rank,rank)
										if(.not. allocated(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%matrix))allocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%matrix(rank,mm))
										call copymatT_OMP(matA,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)											
									! ! end if
									
									if(level_right==unique_nth)then
										deallocate(random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%matrix)
										! call copymatN_omp(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j)%matrix,matC(1:rank,1:rank),rank,rank)
										! call copymatN_omp(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j+1)%matrix,matC(1:rank,rank+1:2*rank),rank,rank)	
										! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%null,mm,rank,&
										! &2*rank,random_Block(1)%RandomVectorLL(level_butterfly-level_right+1)%blocks(i+1,index_j)%nulldim,up_tolerance,LS_tolerance,error0,kmax)
										call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
										if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j)%matrix(rank,rank))
										if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j+1)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j+1)%matrix(rank,rank))										
										call copymatN_omp(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j)%matrix,rank,rank)
										call copymatN_omp(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_right)%blocks(i+1,j+1)%matrix,rank,rank)									
									end if
								end if

								! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_right,nth,i,j,error0,'L'
						   end if   
					   enddo
				   enddo
			   enddo
			   deallocate(matB,matC,matA,matinv)
		   endif	   
	   end do
	   
   endif
   
   
   return

end subroutine Resolving_Butterfly_LL_new_uniform


subroutine Resolving_Butterfly_RR_new_uniform(num_vect_sub,nth_s,nth_e,Ng,unique_nth)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank,level_right_start,level_left_start,level_butterfly
   integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, flag
   real(kind=8) a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   
   ! type(matricesblock), pointer :: blocks
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real(kind=8), allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer num_vect_sub,num_vect_subsub,nth,nulldim,ind_r,ind_c,noe,Ng,nth_s,nth_e,dimension_m,dimension_n
   real*8::n1,n2
   integer::kmax,unique_nth
   
   norm3R=0. ; norm4R=0.
   level_butterfly=butterfly_block_randomized(1)%level_butterfly 
   num_blocks=2**level_butterfly
   rank=butterfly_block_randomized(1)%dimension_rank
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
   
   if (level_butterfly/=0) then
       mm=num_vect_subsub !rank
       nn=2*rank
       error=1000.
	   iter = 1	
	   	   
	   level_left_start= floor_safe(level_butterfly/2d0)+1    !  check here later		   
	   ! ! level_left_start = 0
	   
	   random=>random_Block(1)
	   if(level_left_start>0 .and. level_left_start==unique_nth)then
			n1 = OMP_get_wtime()
			call Butterfly_partial_MVP_Half('N',0,level_left_start-1,random,num_vect_sub,nth_s,nth_e,Ng)
			! call Butterfly_partial_MVP('N',0,level_left_start-1,random)
			n2 = OMP_get_wtime()
			time_halfbuttermul = time_halfbuttermul + n2-n1		 
		endif 
	   
	   
	   do level_left = level_butterfly+1,unique_nth,-1
			! kmax = ceiling_safe(rank/dble(2**(level_left-level_left_start)))+1
			kmax = rank
			if (level_left==level_butterfly+1) then
			   do nth=nth_s,nth_e
				   do i=1, num_blocks
					   dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)	
					   allocate(matB(mm,dimension_m),matC(rank,dimension_m),matA(mm,rank),matinv(dimension_m,rank))	
					   call copymatT_OMP(random_Block(1)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_m,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_m,mm)
					   call copymatT_OMP(butterfly_block_randomized(1)%KerInv(1:rank,1:dimension_m),matinv,rank,dimension_m)
					   call gemm_omp(matB,matinv,matA,mm,dimension_m,rank)
					   if(.not. allocated(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix))allocate(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(rank,mm))
					   
					   call copymatT_OMP(matA,random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
					   
					   if(level_left==unique_nth)then
							! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyU(i)%matrix,matC,dimension_m,rank)
							! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%null,mm,rank,&
							! &dimension_m,random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%nulldim,up_tolerance,LS_tolerance,error0,kmax)
							call LeastSquare(mm,rank,dimension_m,matA,matB,matC,LS_tolerance)
							if(.not. allocated(butterfly_block_randomized(1)%ButterflyU(i)%matrix))allocate(butterfly_block_randomized(1)%ButterflyU(i)%matrix(dimension_m,rank))
							call copymatT_OMP(matC,butterfly_block_randomized(1)%ButterflyU(i)%matrix,rank,dimension_m)
					   end if
					   ! write(*,'(I5,I5,I5,Es16.7E3,A2)')unique_nth,nth,i,error0,'R'
					   deallocate(matB,matC,matA,matinv)
					end do
				end do
			elseif (level_left==0) then
				write(*,*)'the left half scheme should not touch rightmost matrix'
				stop
			else 
			   allocate (matB(mm,rank*2),matC(rank,rank*2),matA(mm,rank),matinv(rank*2,rank))
			   
			   num_row=butterfly_block_randomized(1)%ButterflyKerl(level_left)%num_row
			   num_col=butterfly_block_randomized(1)%ButterflyKerl(level_left)%num_col
			   do nth = nth_s,nth_e
				   noe = ceiling_safe(nth*Ng*2d0/num_row)
				   ind_c = ceiling_safe(noe/2d0)
				   
				   do j=1, num_col, 2
					   index_j=int((j+1)/2)
					   do i=1, num_row, 2
						   index_i=int((i+1)/2)
						   if(index_j==ind_c)then
								
								call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:rank),rank,mm)
								call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+rank:rank*2),rank,mm)
								
								deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix)
								deallocate(random_Block(1)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix)
								
								! ! ! ! ! if(level_left==level_left_start)then
									! ! ! ! ! if(mod(noe,2)==1)then
										! ! ! ! ! call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
										! ! ! ! ! call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
										! ! ! ! ! call copymatT_OMP(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j)%matrix,rank,rank)
										! ! ! ! ! call copymatT_OMP(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j)%matrix,rank,rank)								
									! ! ! ! ! else 
										! ! ! ! ! call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
										! ! ! ! ! call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
										! ! ! ! ! call copymatT_OMP(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j+1)%matrix,rank,rank)
										! ! ! ! ! call copymatT_OMP(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j+1)%matrix,rank,rank)									
									! ! ! ! ! end if
								! ! ! ! ! else 
									if(mod(noe,2)==1)then
										if(level_left==level_left_start)then
											call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
										else 
											
											call copymatT_OMP(butterfly_block_randomized(1)%KerInv(1:rank,1:2*rank),matinv,rank,2*rank)
											call gemm_omp(matB,matinv,matA,mm,2*rank,rank)
											if(.not. allocated(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix))allocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix(rank,mm))
											call copymatT_OMP(matA,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
										end if
										
										if(level_left==unique_nth)then
											! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j)%matrix,matC(1:rank,1:rank),rank,rank)
											! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j)%matrix,matC(1:rank,rank+1:2*rank),rank,rank)	
											! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%null,mm,rank,&
											! &2*rank,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%nulldim,up_tolerance,LS_tolerance,error0,kmax)
											call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
											if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j)%matrix(rank,rank))
											if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j)%matrix(rank,rank))
											call copymatT_OMP(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j)%matrix,rank,rank)
											call copymatT_OMP(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j)%matrix,rank,rank)
											deallocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j)%matrix)
										end if
									else
										if(level_left==level_left_start)then
											call copymatT_OMP(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
										else 
											call copymatT_OMP(butterfly_block_randomized(1)%KerInv(1+rank:2*rank,1:2*rank),matinv,rank,2*rank)
											
											call gemm_omp(matB,matinv,matA,mm,2*rank,rank)
											if(.not. allocated(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix))allocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix(rank,mm))
											call copymatT_OMP(matA,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
										end if
										if(level_left==unique_nth)then
											! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j+1)%matrix,matC(1:rank,1:rank),rank,rank)
											! call copymatT_OMP(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j+1)%matrix,matC(1:rank,rank+1:2*rank),rank,rank)	
											! call KernelUpdate(matA,matB,matC,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%null,mm,rank,&
											! &2*rank,random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%nulldim,up_tolerance,LS_tolerance,error0,kmax)
											call LeastSquare(mm,rank,2*rank,matA,matB,matC,LS_tolerance)
											if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j+1)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j+1)%matrix(rank,rank))
											if(.not. allocated(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j+1)%matrix))allocate(butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j+1)%matrix(rank,rank))											
											call copymatT_omp(matC(1:rank,1:rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i,j+1)%matrix,rank,rank)
											call copymatT_omp(matC(1:rank,rank+1:2*rank),butterfly_block_randomized(1)%ButterflyKerl(level_left)%blocks(i+1,j+1)%matrix,rank,rank)
											deallocate(random_Block(1)%RandomVectorRR(level_left)%blocks(index_i,j+1)%matrix)
										end if							
									end if
								! ! ! ! ! end if
								! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_left,nth,i,j,error0,'R'
																
						   end if   
					   enddo
				   enddo
			   enddo	   
			   deallocate(matB,matC,matA,matinv)
			end if
	   end do
   endif
   
   
   return

end subroutine Resolving_Butterfly_RR_new_uniform




subroutine butterfly_block_MVP_inverse_dat(level,ii,trans,N,num_vect_sub,Vin,Vout)
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
   
	block_off1 => cascading_factors(level)%matrices_block(ii*2-1)	
	block_off2 => cascading_factors(level)%matrices_block(ii*2)
	
	groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	call assert(mm+nn==N,'mm+nn/=N')  
		
   if(schurinv==0)then
		block_o => cascading_factors(level)%matrices_block_inverse(ii)
		call butterfly_block_MVP_randomized_dat(block_o,trans,N,N,num_vect_sub,&
		&Vin(1:N,1:num_vect_sub),Vout(1:N,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:N,1:num_vect_sub) = Vin(1:N,1:num_vect_sub) + Vout(1:N,1:num_vect_sub)
   else 
		block_o => cascading_factors(level)%matrices_block_inverse_schur(ii)
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
   
end subroutine butterfly_block_MVP_inverse_dat









end module Butterfly_rightmultiply
