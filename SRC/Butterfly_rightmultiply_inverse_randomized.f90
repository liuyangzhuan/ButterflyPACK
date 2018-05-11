module Butterfly_rightmultiply
use Utilites_randomized
! use Butterfly_compression_givenfullmat
use omp_lib
contains






subroutine Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,LR)
	use misc
    use MODULE_FILE
    implicit none
	
	type(matrixblock)::block_o
	integer level_butterfly,level_butterfly_loc, ij_loc,index_i,index_i_start,index_j_start,index_j,level,ii,nn,mm,num_blocks,rank
	character LR
	
	allocate(agent_block(1))
	

	call assert(level_butterfly_loc>=1,'level_butterfly_loc cannot be zero')

	agent_block(1)%row_group=-1
	agent_block(1)%col_group=-1
	
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
type(matrixblock)::block_i

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



subroutine Copy_butterfly_partial(block_i,block_o,level_butterfly_loc,ij_loc,LR,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_o,block_i

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

end subroutine Copy_butterfly_partial



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
	block_o =>  ho_bf%levels(level_c)%matrices_block(rowblock) 
	  
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
				groupm_diag = ho_bf%levels(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
				idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
				idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1
				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vect_sub,mm !,block_o%col_group,basis_group(block_o%col_group)%head

				call fullmat_block_MVP_randomized_dat(ho_bf%levels(level)%matrices_block_inverse(ii),'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
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
				groupm_diag = ho_bf%levels(level)%matrices_block_inverse(ii)%row_group ! Note: row_group and col_group interchanged here   

				
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











subroutine Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,blocks)

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
   type(matrixblock) :: blocks
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   ! blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
	
   ! if (level_butterfly/=0) then
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
						call OneV_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)	   
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
						call OneKernel_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
				   enddo
				   !$omp end parallel do
			   enddo
		   endif	   
	   end do
	   
   ! endif
   
   
   return

end subroutine Resolving_Butterfly_LL_new

subroutine OneV_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(matrixblock) :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,level_butterfly
	
   ! blocks => butterfly_block_randomized(1)   
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


subroutine OneKernel_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(matrixblock) :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer index_i,index_j,i,j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,ieo,noe,rs,re,level_butterfly

   ! blocks => butterfly_block_randomized(1)   
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



subroutine Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,blocks)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank,level_right_start,level_left_start,nn1,nn2,rs,re
   integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, flag
   real*8 a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   
   type(matrixblock) :: blocks
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real*8, allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,ind_c,noe,Ng,nth_s,nth_e,dimension_mm,dimension_n,jeo,level_butterfly
   real*8::n1,n2
   integer::kmax,unique_nth
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   ! blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 

   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
   
   ! if (level_butterfly/=0) then
       mm=num_vect_subsub !rank

	   level_left_start= floor_safe(level_butterfly/2d0)+1    !  check here later		   
	   ! ! level_left_start = 0
	   
	   random=>random_Block(1)
	   if(level_left_start>0 .and. level_left_start==unique_nth)then
			n1 = OMP_get_wtime()
			call Butterfly_partial_MVP_Half(blocks,'N',0,level_left_start-1,random,num_vect_sub,nth_s,nth_e,Ng)
			n2 = OMP_get_wtime()
			time_halfbuttermul = time_halfbuttermul + n2-n1		 
		endif 
	   
	   
	   do level_left = level_butterfly+1,unique_nth,-1
			if (level_left==level_butterfly+1) then
				do nth=nth_s,nth_e
					!$omp parallel do default(shared) private(i)
					do i=1, num_blocks   
						call OneU_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
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
					   call OneKernel_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)									
				   enddo
				   !$omp end parallel do	
			   enddo	   
			end if
	   end do
   ! endif
   
   
   return

end subroutine Resolving_Butterfly_RR_new


subroutine OneU_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(matrixblock) :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer i,level_left,unique_nth,dimension_mm,mm,rank,num_vect_sub,nth,nth_s,level_butterfly
	
   ! blocks => butterfly_block_randomized(1)   
   level_butterfly=blocks%level_butterfly 
   
	if(level_left==unique_nth)then
		if(level_butterfly>0)then
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
			allocate(matB(mm,dimension_mm))
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)									
			rank = size(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix,1)
			if(allocated(blocks%butterflyU(i)%matrix))deallocate(blocks%butterflyU(i)%matrix)
			allocate(blocks%butterflyU(i)%matrix(dimension_mm,rank))
			allocate(matC(rank,dimension_mm),matA(mm,rank))	
			call copymatT_omp(random_Block(1)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,dimension_mm,matA,matB,matC,LS_tolerance)
			! write(*,*)fnorm(matC,rank,dimension_mm),'U',level_left,level_butterfly
			call copymatT_omp(matC,blocks%ButterflyU(i)%matrix,rank,dimension_mm)							
			deallocate(matB,matC,matA)			
		endif			
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



subroutine OneKernel_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s,blocks)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(matrixblock) :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer index_i,index_j,i,j,level_left,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,jeo,noe,rs,re,level_left_start,level_butterfly

   ! blocks => butterfly_block_randomized(1)   
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
   
	block_off1 => ho_bf%levels(level)%matrices_block(ii*2-1)	
	block_off2 => ho_bf%levels(level)%matrices_block(ii*2)
	
	groupn=block_off1%col_group    ! Note: row_group and col_group interchanged here   
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1     
	groupm=block_off1%row_group    ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	call assert(mm+nn==N,'mm+nn/=N')  
		
   if(schurinv==0)then
		block_o => ho_bf%levels(level)%matrices_block_inverse(ii)
		call butterfly_block_MVP_randomized_dat(block_o,trans,N,N,num_vect_sub,&
		&Vin(1:N,1:num_vect_sub),Vout(1:N,1:num_vect_sub),ctemp1,ctemp2)
		Vout(1:N,1:num_vect_sub) = Vin(1:N,1:num_vect_sub) + Vout(1:N,1:num_vect_sub)
   else 
		block_o => ho_bf%levels(level)%matrices_block_inverse_schur(ii)
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
