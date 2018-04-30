module Utilities
use misc
contains
 

real*8 function memory_of_blocks_func(blocks)
    
    use MODULE_FILE
    implicit none
    
    integer blocks, butterflyB_inuse, level_actual, num_col, num_row
    integer i, j, mm, nn, rank, num_blocks, level, level_butterfly
    real*8 memory_butterfly, rtemp
    
    level_butterfly=matrices_block(blocks,0)%level_butterfly
    level_actual=Maxlevel_for_blocks-matrices_block(blocks,0)%level
    
    memory_butterfly=0.
    
    if (matrices_block(blocks,0)%style==2) then
    
        rtemp=0
        !$omp parallel do default(shared) private(i,mm,nn,rank) reduction(+:rtemp)
        do i=1, 2**level_actual
            mm=size(matrices_block(blocks,0)%ButterflyU(i)%matrix,1)
            rank=size(matrices_block(blocks,0)%ButterflyU(i)%matrix,2)
            rtemp=rtemp+real(mm*rank)
            nn=size(matrices_block(blocks,0)%ButterflyV(i)%matrix,1)
            rank=size(matrices_block(blocks,0)%ButterflyV(i)%matrix,2)
            rtemp=rtemp+real(nn*rank)
        enddo
        !$omp end parallel do
        memory_butterfly=memory_butterfly+rtemp
                
        if (level_butterfly/=0) then
            do level=1, level_butterfly
                num_col=matrices_block(blocks,0)%ButterflyKerl(level)%num_col
                num_row=matrices_block(blocks,0)%ButterflyKerl(level)%num_row
                rtemp=0
                !$omp parallel do default(shared) private(i,j,mm,nn) reduction(+:rtemp)
                do j=1, num_col
                    do i=1, num_row
                        mm=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
                        nn=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
                        rtemp=rtemp+real(mm*nn)
                    enddo
                enddo
                !$omp end parallel do
                memory_butterfly=memory_butterfly+rtemp
            enddo
        endif
        
    elseif (matrices_block(blocks,0)%style==1) then
        
        mm=size(matrices_block(blocks,0)%fullmat,1)
        nn=size(matrices_block(blocks,0)%fullmat,2)
        memory_butterfly=real(mm*nn)
    
    endif
    
    memory_of_blocks_func=memory_butterfly
   
    return

end function memory_of_blocks_func



subroutine delete_blocks(blocks)

    use MODULE_FILE
    implicit none
    
    integer butterflyB_inuse, level_actual, num_col, num_row
    integer i, j, mm, nn, rank, num_blocks, level, level_butterfly,index_i_m,index_j_m,levelm
    real*8 memory_butterfly, rtemp
    type(matrixblock)::blocks
	
    if (blocks%style/=1) then
    
        level_butterfly=blocks%level_butterfly
		
        ! level_actual=Maxlevel_for_blocks-blocks%level
        level_actual=level_butterfly
		
		if(allocated(blocks%ButterflyU))then
        !$omp parallel do default(shared) private(i)
        do i=1, 2**level_actual
													
            if(allocated(blocks%ButterflyU(i)%matrix))deallocate (blocks%ButterflyU(i)%matrix)
        enddo
        !$omp end parallel do
        deallocate (blocks%ButterflyU)
        end if

		if(allocated(blocks%ButterflyV))then
        !$omp parallel do default(shared) private(i)
        do i=1, 2**level_actual
            if(allocated(blocks%ButterflyV(i)%matrix))deallocate (blocks%ButterflyV(i)%matrix)
        enddo
        !$omp end parallel do
        deallocate (blocks%ButterflyV)
        end if

		
		if(allocated(blocks%ButterflyKerl))then
        if (level_butterfly/=0) then
            !$omp parallel do default(shared) private(level,i,j,num_col,num_row)
            do level=1, level_butterfly
                num_col=blocks%ButterflyKerl(level)%num_col
                num_row=blocks%ButterflyKerl(level)%num_row
                do j=1, num_col
                    do i=1, num_row
                        if(allocated(blocks%ButterflyKerl(level)%blocks(i,j)%matrix))deallocate (blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
                    enddo
                enddo
                deallocate (blocks%ButterflyKerl(level)%blocks)
            enddo
            !$omp end parallel do
            deallocate (blocks%ButterflyKerl)
        endif
		end if
		
		if(allocated(blocks%ButterflyMiddle))then
			levelm = ceiling_safe(dble(level_butterfly)/2d0)
			do index_i_m=1, 2**levelm
				do index_j_m=1, 2**(level_butterfly-levelm)	
					if(allocated(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix))deallocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix)
				end do
			end do
			deallocate(blocks%ButterflyMiddle)
		end if	
	
        blocks%level_butterfly=0
		blocks%rankmax = -1000
		blocks%rankmin = 1000
        
    elseif (blocks%style==1) then
        
        if(allocated(blocks%fullmat))deallocate (blocks%fullmat)
    
    endif
    
    return

end subroutine delete_blocks




subroutine delete_Bplus(bplus)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block
type(blockplus)::bplus

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::rtemp


do ll=1,LplusMax
	if(bplus%LL(ll)%Nbound>0)then
		do bb=1,bplus%LL(ll)%Nbound
			! write(*,*)ll,bplus%Lplus,bb,bplus%LL(ll)%Nbound,'fff'
			call delete_blocks(bplus%LL(ll)%matrices_block(bb))
		end do		
		deallocate(bplus%LL(ll)%matrices_block)
		if(allocated(bplus%LL(ll)%boundary_map))deallocate(bplus%LL(ll)%boundary_map)
	end if
end do
deallocate(bplus%LL)

end subroutine delete_Bplus



subroutine delete_Bplus_onlyButterflydata(bplus)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block
type(blockplus)::bplus

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::rtemp


do ll=1,LplusMax
	if(bplus%LL(ll)%Nbound>0)then
		do bb=1,bplus%LL(ll)%Nbound
			! write(*,*)ll,bplus%Lplus,bb,bplus%LL(ll)%Nbound,'fff'
			call delete_blocks(bplus%LL(ll)%matrices_block(bb))
		end do		
		! deallocate(bplus%LL(ll)%matrices_block)
		! if(allocated(bplus%LL(ll)%boundary_map))deallocate(bplus%LL(ll)%boundary_map)
	end if
end do
! deallocate(bplus%LL)

end subroutine delete_Bplus_onlyButterflydata





subroutine copy_butterfly(block_i,block_o,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i,block_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	
real*8,optional::memory
if(present(memory))memory=0
block_o%level = block_i%level
block_o%col_group = block_i%col_group
block_o%row_group = block_i%row_group
block_o%nested_num = block_i%nested_num
block_o%style = block_i%style
block_o%data_type = block_i%data_type
block_o%level_butterfly = block_i%level_butterfly
block_o%rankmax = block_i%rankmax
block_o%rankmin = block_i%rankmin
level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

if(block_i%style==2)then

	allocate(block_o%ButterflyU(num_blocks))
	allocate(block_o%ButterflyV(num_blocks))
	if (level_butterfly/=0) then
		allocate(block_o%ButterflyKerl(level_butterfly))
	end if

	do level=0, level_butterfly
		index_ij=0
		if (level>0) then
			block_o%ButterflyKerl(level)%num_row=2**level
			block_o%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
			allocate(block_o%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
		endif
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly-level)
				index_ij=index_ij+1
				if (level==0) then
					nn=size(block_i%ButterflyV(index_ij)%matrix,1)
					rank=size(block_i%ButterflyV(index_ij)%matrix,2)
					allocate(block_o%ButterflyV(index_ij)%matrix(nn,rank))
					block_o%ButterflyV(index_ij)%matrix = block_i%ButterflyV(index_ij)%matrix
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV(index_ij)%matrix)/1024.0d3
				else                    
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					allocate(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn))
					block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					allocate(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn))
					block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix                    
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
				endif
				if (level==level_butterfly) then
					mm=size(block_i%ButterflyU(index_ij)%matrix,1)
					rank=size(block_i%ButterflyU(index_ij)%matrix,2)
					allocate(block_o%ButterflyU(index_ij)%matrix(mm,rank))
					block_o%ButterflyU(index_ij)%matrix = block_i%ButterflyU(index_ij)%matrix					
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU(index_ij)%matrix)/1024.0d3
				endif
			enddo
		enddo
	enddo
			
	if(allocated(block_i%ButterflyMiddle))then
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		allocate(block_o%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))

		do index_i_m=1, 2**levelm
			do index_j_m=1, 2**(level_butterfly-levelm)	
				rank = size(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,1)
				allocate(block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix
			end do
		end do
	end if
else if(block_i%style==1)then
	mm = size(block_i%fullmat,1)
	nn = size(block_i%fullmat,2)
	allocate (block_o%fullmat(mm,nn))
	block_o%fullmat = block_i%fullmat
else 
	write(*,*)'block style not implemented'
	stop
end if

end subroutine copy_butterfly

subroutine copy_delete_butterfly(block_i,block_o,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i,block_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	
real*8,optional::memory
if(present(memory))memory=0
block_o%level = block_i%level
block_o%col_group = block_i%col_group
block_o%row_group = block_i%row_group
block_o%nested_num = block_i%nested_num
block_o%style = block_i%style
block_o%data_type = block_i%data_type
block_o%level_butterfly = block_i%level_butterfly
block_o%rankmax = block_i%rankmax
block_o%rankmin = block_i%rankmin
level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

if(block_i%style==2)then

	allocate(block_o%ButterflyU(num_blocks))
	allocate(block_o%ButterflyV(num_blocks))
	if (level_butterfly/=0) then
		allocate(block_o%ButterflyKerl(level_butterfly))
	end if

	do level=0, level_butterfly
		index_ij=0
		if (level>0) then
			block_o%ButterflyKerl(level)%num_row=2**level
			block_o%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
			allocate(block_o%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
		endif
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly-level)
				index_ij=index_ij+1
				if (level==0) then
					nn=size(block_i%ButterflyV(index_ij)%matrix,1)
					rank=size(block_i%ButterflyV(index_ij)%matrix,2)
					allocate(block_o%ButterflyV(index_ij)%matrix(nn,rank))
					block_o%ButterflyV(index_ij)%matrix = block_i%ButterflyV(index_ij)%matrix
					deallocate(block_i%ButterflyV(index_ij)%matrix)
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyV(index_ij)%matrix)/1024.0d3
				else                    
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					allocate(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn))
					block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix
					deallocate(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					allocate(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn))
					block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix = block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix                    
					deallocate(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
				endif
				if (level==level_butterfly) then
					mm=size(block_i%ButterflyU(index_ij)%matrix,1)
					rank=size(block_i%ButterflyU(index_ij)%matrix,2)
					allocate(block_o%ButterflyU(index_ij)%matrix(mm,rank))
					block_o%ButterflyU(index_ij)%matrix = block_i%ButterflyU(index_ij)%matrix					
					deallocate(block_i%ButterflyU(index_ij)%matrix)
					if(present(memory))memory = memory + SIZEOF(block_o%ButterflyU(index_ij)%matrix)/1024.0d3
				endif
			enddo
		enddo
		if (level>0) then
			deallocate(block_i%ButterflyKerl(level)%blocks)
		end if
	enddo
	deallocate(block_i%ButterflyU)
	deallocate(block_i%ButterflyV)
	if(level_butterfly/=0)deallocate(block_i%ButterflyKerl)
	
	if(allocated(block_i%ButterflyMiddle))then
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		allocate(block_o%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))

		do index_i_m=1, 2**levelm
			do index_j_m=1, 2**(level_butterfly-levelm)	
				rank = size(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix,1)
				allocate(block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				block_o%ButterflyMiddle(index_i_m,index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix
				deallocate(block_i%ButterflyMiddle(index_i_m,index_j_m)%matrix)
			end do
		end do
		deallocate(block_i%ButterflyMiddle)
	end if
else if(block_i%style==1)then
	mm = size(block_i%fullmat,1)
	nn = size(block_i%fullmat,2)
	allocate (block_o%fullmat(mm,nn))
	block_o%fullmat = block_i%fullmat
	deallocate(block_i%fullmat)
else 
	write(*,*)'block style not implemented'
	stop
end if

							 

end subroutine copy_delete_butterfly

subroutine ComputeMemory_butterfly(block_i,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	
real*8::memory
memory=0

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

if(block_i%style==2)then

	do level=0, level_butterfly
		index_ij=0
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly-level)
				index_ij=index_ij+1
				if (level==0) then
					memory = memory + SIZEOF(block_i%ButterflyV(index_ij)%matrix)/1024.0d3
				else                    
					memory = memory + SIZEOF(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3                  
					memory = memory + SIZEOF(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
				endif
				if (level==level_butterfly) then				
					memory = memory + SIZEOF(block_i%ButterflyU(index_ij)%matrix)/1024.0d3
				endif
			enddo
		enddo
	enddo
			
else if(block_i%style==1)then
	memory = memory + SIZEOF(block_i%fullmat)/1024.0d3
else 
	write(*,*)'block style not implemented'
	stop
end if

end subroutine ComputeMemory_butterfly






logical function CheckNAN_butterfly(block_i)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	
real*8:: temp

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly
temp = 0

if(block_i%style==2)then

	do level=0, level_butterfly
		index_ij=0
		do index_i=1, 2**level
			do index_j=1, 2**(level_butterfly-level)
				index_ij=index_ij+1
				if (level==0) then
					mm=size(block_i%ButterflyV(index_ij)%matrix,1)
					nn=size(block_i%ButterflyV(index_ij)%matrix,2)
					temp = temp + fnorm(block_i%ButterflyV(index_ij)%matrix,mm,nn)
				else
					mm=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
					temp = temp + fnorm(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,mm,nn)
					
					mm=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,1)
					nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,2)
					temp = temp + fnorm(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix,mm,nn)
				endif
				if (level==level_butterfly) then				
					mm=size(block_i%ButterflyU(index_ij)%matrix,1)
					nn=size(block_i%ButterflyU(index_ij)%matrix,2)
					temp = temp + fnorm(block_i%ButterflyU(index_ij)%matrix,mm,nn)
				endif
			enddo
		enddo
	enddo
			
else if(block_i%style==1)then
	mm=size(block_i%fullmat,1)
	nn=size(block_i%fullmat,2)
	temp = temp + fnorm(block_i%fullmat,mm,nn)
else 
	write(*,*)'block style not implemented'
	stop
end if

CheckNAN_butterfly = isnan(temp)

end function CheckNAN_butterfly







subroutine copy_Bplus(bplus_i,bplus_o,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8,optional::memory
real*8::rtemp

if(present(memory))memory=0

allocate(bplus_o%LL(LplusMax))
bplus_o%Lplus = bplus_i%Lplus
bplus_o%boundary = bplus_i%boundary
bplus_o%level = bplus_i%level
bplus_o%col_group = bplus_i%col_group
bplus_o%row_group = bplus_i%row_group


do ll=1,LplusMax
	bplus_o%LL(ll)%Nbound=bplus_i%LL(ll)%Nbound
	bplus_o%LL(ll)%rankmax=bplus_i%LL(ll)%rankmax
	if(bplus_i%LL(ll)%Nbound>0)then
		allocate(bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
		do bb=1,bplus_i%LL(ll)%Nbound
			call copy_butterfly(bplus_i%LL(ll)%matrices_block(bb),bplus_o%LL(ll)%matrices_block(bb),rtemp)
			if(present(memory))memory=memory+rtemp
		end do
		Nboundall=size(bplus_i%LL(ll)%boundary_map)
		allocate(bplus_o%LL(ll)%boundary_map(Nboundall))
		if(present(memory))memory=memory+ SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
		bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
	end if
end do

end subroutine copy_Bplus


subroutine copy_delete_Bplus(bplus_i,bplus_o,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8,optional::memory
real*8::rtemp

if(present(memory))memory=0

allocate(bplus_o%LL(LplusMax))
bplus_o%Lplus = bplus_i%Lplus
bplus_o%boundary = bplus_i%boundary
bplus_o%level = bplus_i%level
bplus_o%col_group = bplus_i%col_group
bplus_o%row_group = bplus_i%row_group


do ll=1,LplusMax
	bplus_o%LL(ll)%Nbound=bplus_i%LL(ll)%Nbound
	bplus_o%LL(ll)%rankmax=bplus_i%LL(ll)%rankmax
	if(bplus_i%LL(ll)%Nbound>0)then
		allocate(bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
		do bb=1,bplus_i%LL(ll)%Nbound
			call copy_delete_butterfly(bplus_i%LL(ll)%matrices_block(bb),bplus_o%LL(ll)%matrices_block(bb),rtemp)
			if(present(memory))memory=memory+rtemp
		end do
		deallocate(bplus_i%LL(ll)%matrices_block)
		Nboundall=size(bplus_i%LL(ll)%boundary_map)
		allocate(bplus_o%LL(ll)%boundary_map(Nboundall))
		if(present(memory))memory=memory+ SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
		bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
		deallocate(bplus_i%LL(ll)%boundary_map)
	end if
end do

deallocate(bplus_i%LL)
end subroutine copy_delete_Bplus


subroutine ComputeMemory_Bplus(bplus_i,memory)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::memory
real*8::rtemp

memory=0

do ll=1,LplusMax
	if(bplus_i%LL(ll)%Nbound>0)then
		do bb=1,bplus_i%LL(ll)%Nbound
			call ComputeMemory_butterfly(bplus_i%LL(ll)%matrices_block(bb),rtemp)
			memory=memory+rtemp
		end do
	end if
end do

end subroutine ComputeMemory_Bplus




logical function CheckNAN_Bplus(bplus_i)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i,bplus_o

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::rtemp
CheckNAN_Bplus = .false.

do ll=1,LplusMax
	if(bplus_i%LL(ll)%Nbound>0)then
		do bb=1,bplus_i%LL(ll)%Nbound
			if(CheckNAN_butterfly(bplus_i%LL(ll)%matrices_block(bb)))then
				CheckNAN_Bplus = .true.
				return
			end if
		end do
	end if
end do

end function CheckNAN_Bplus





subroutine print_butterfly_size_rank(block_i)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,truerank,index_i,index_j,levelm,index_i_m,index_j_m,mm1,mm2,nn1,nn2
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	
complex(kind=8),allocatable::matrixtemp(:,:),mat11(:,:),mat12(:,:),mat21(:,:),mat22(:,:)

level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

do level=0, level_butterfly+1
	! write(*,*)level
	if (level==0) then
		do index_ij=1, 2**level_butterfly
			nn=size(block_i%ButterflyV(index_ij)%matrix,1)
			rank=size(block_i%ButterflyV(index_ij)%matrix,2)
			allocate(matrixtemp(nn,rank))
			matrixtemp = block_i%ButterflyV(index_ij)%matrix
			call GetRank(nn,rank,matrixtemp,truerank,SVD_tolerance_forward)
			write(*,*)level,index_ij,nn,rank,truerank
			deallocate(matrixtemp)
		end do
	else if (level==level_butterfly+1) then
		do index_ij=1, 2**level_butterfly
			mm=size(block_i%ButterflyU(index_ij)%matrix,1)
			rank=size(block_i%ButterflyU(index_ij)%matrix,2)	
			allocate(matrixtemp(mm,rank))
			matrixtemp = block_i%ButterflyU(index_ij)%matrix
			call GetRank(mm,rank,matrixtemp,truerank,SVD_tolerance_forward)			
			write(*,*)level,index_ij,mm,rank,truerank	
			deallocate(matrixtemp)
		end do
	else
		do index_i=1, 2**(level-1)
			do index_j=1, 2**(level_butterfly-level)
                    
				mm1=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix,1)
				nn1=size(block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix,2)
				
				mm2=size(block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix,1)
				nn2=size(block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix,2)
				
				allocate(mat11(mm1,nn1))
				mat11 = block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j-1)%matrix
				allocate(mat12(mm1,nn2))
				mat12 = block_i%ButterflyKerl(level)%blocks(2*index_i-1,2*index_j)%matrix
				allocate(mat21(mm2,nn1))
				mat21 = block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j-1)%matrix				
				allocate(mat22(mm2,nn2))
				mat22 = block_i%ButterflyKerl(level)%blocks(2*index_i,2*index_j)%matrix
				allocate(matrixtemp(mm1+mm2,nn1+nn2))
				matrixtemp(1:mm1,1:nn1) = mat11
				matrixtemp(1:mm1,1+nn1:nn2+nn1) = mat12
				matrixtemp(1+mm1:mm2+mm1,1:nn1) = mat21				
				matrixtemp(1+mm1:mm2+mm1,1+nn1:nn2+nn1) = mat22				
				call GetRank(mm1+mm2,nn1+nn2,matrixtemp,truerank,SVD_tolerance_forward)		
				write(*,*)level,index_i,index_j,(mm1+mm2),(nn1+nn2),truerank					
				deallocate(mat11)
				deallocate(mat12)
				deallocate(mat21)
				deallocate(mat22)
				deallocate(matrixtemp)
				
			enddo
		enddo		
	
	end if

enddo

end subroutine print_butterfly_size_rank




subroutine butterfly_block_MVP_randomized(blocks,chara,random1,random2,a,b)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1,k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,num_col,num_row
    complex(kind=8) ctemp, a, b
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
    character chara
	type(matrixblock)::blocks
	integer middleflag,levelm
	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)
    
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    type(vectorsblock), pointer :: random1, random2
    
    num_vectors=size(random1%vector,2)
    middleflag = 0
	if(allocated(blocks%ButterflyMiddle))middleflag=1
	
    level_butterfly=blocks%level_butterfly
    num_blocks=2**level_butterfly	
	allocate(arr_acc_m(num_blocks))
	allocate(arr_acc_n(num_blocks))
	
	k1=0
	k2=0
	do i=1, num_blocks
		arr_acc_n(i) = k1
		arr_acc_m(i) = k2
		nn=size(blocks%ButterflyV(i)%matrix,1)
		k1 =k1 +nn
		mm=size(blocks%ButterflyU(i)%matrix,1)
		k2 =k2 +mm		
	enddo 
	
	
    if (chara=='N') then
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
        levelm = ceiling_safe(dble(level_butterfly)/2d0)
		
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(1,num_blocks))
        ButterflyVector(0)%num_row=1
        ButterflyVector(0)%num_col=num_blocks
        
		
		!$omp parallel do default(shared) private(i,nn,ii,jj)
		do i=1, num_blocks
            nn=size(blocks%ButterflyV(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(1,i)%matrix(nn,num_vectors))
            do ii=1, nn
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(1,i)%matrix(ii,jj)=random1%vector(ii+arr_acc_n(i),jj)
                enddo
            enddo
		enddo 
		!$omp end parallel do
        
        do level=0, level_butterfly
            if (level==0) then
                num_groupn=num_blocks
                allocate (ButterflyVector(1)%blocks(1,num_groupn))
                ButterflyVector(1)%num_row=1
                ButterflyVector(1)%num_col=num_groupn
                ! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_groupn
                    rank=size(blocks%ButterflyV(j)%matrix,2)
                    nn=size(blocks%ButterflyV(j)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(1,j)%matrix(rank,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
					do jj=1, num_vectors						
						do ii=1, rank
                            ctemp=0
                            do kk=1, nn
                                ctemp=ctemp+blocks%ButterflyV(j)%matrix(kk,ii)*ButterflyVector(0)%blocks(1,j)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo                    
				! !$omp end parallel do
            else
                num_groupm=blocks%ButterflyKerl(level)%num_row
                num_groupn=blocks%ButterflyKerl(level)%num_col
                if (num_groupn/=1) then
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,int(num_groupn/2)))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=int(num_groupn/2)                    
                    
                    ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
					do ij=1,num_groupm*(num_groupn/2)
						i = (ij-1)/(num_groupn/2)+1
						j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
						index_i=int((i+1)/2)
						index_j=int((j+1)/2)	
						
						nn1=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2=size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mm=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						allocate (ButterflyVector(level+1)%blocks(i,index_j)%matrix(mm,num_vectors))
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vectors
							do ii=1, mm
								ctemp=0
								do kk=1, nn1
									ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j)%matrix(kk,jj)
								enddo
								do kk=1, nn2
									ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j+1)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(i,index_j)%matrix(ii,jj)=ctemp
							enddo
						enddo
						!$omp end parallel do
						if(level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							! ! allocate(matrixtemp(mm,num_vectors))
							! ! call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,matrixtemp,mm,mm,num_vectors)
							! ! ButterflyVector(level+1)%blocks(i,index_j)%matrix = matrixtemp
							! ! deallocate(matrixtemp)
							
							call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,mm,mm,num_vectors)
						end if								
                    enddo
					! !$omp end parallel do
                else
					write(*,*)'this part is not optimzed yet'
					stop
					
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,1))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=1
                    do i=1, num_groupm
                        index_i=int((i+1)/2)
                        nn=size(blocks%ButterflyKerl(level)%blocks(i,1)%matrix,2)
                        mm=size(blocks%ButterflyKerl(level)%blocks(i,1)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(i,1)%matrix(mm,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                        do jj=1, num_vectors						
							do ii=1, mm
                                ctemp=0                            
                                do kk=1, nn
                                    ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,1)%matrix(kk,jj)
                                enddo
                                ButterflyVector(level+1)%blocks(i,1)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(num_blocks,1))
                ButterflyVector(level+2)%num_row=num_blocks
                ButterflyVector(level+2)%num_col=1
                ! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_blocks
                    rank=size(blocks%ButterflyU(i)%matrix,2)
                    mm=size(blocks%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(level+2)%blocks(i,1)%matrix(mm,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do jj=1, num_vectors
						do ii=1, mm
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+blocks%ButterflyU(i)%matrix(ii,kk)*ButterflyVector(level+1)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(level+2)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
				! !$omp end parallel do
            endif
			
			if (level/=0) then
				do j=1, ButterflyVector(level)%num_col
					do i=1, ButterflyVector(level)%num_row
						deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
					enddo
				enddo
			end if	
			if (level==level_butterfly) then			
				do j=1, ButterflyVector(level+1)%num_col
					do i=1, ButterflyVector(level+1)%num_row
						deallocate (ButterflyVector(level+1)%blocks(i,j)%matrix)
					enddo
				enddo			
			end if
        enddo
                                  
		!$omp parallel do default(shared) private(index_i,mm,ii,jj)
		do jj=1, num_vectors
			do index_i=1, num_blocks
				mm=size(blocks%ButterflyU(index_i)%matrix,1)
				do ii=1, mm
									
                    random2%vector(ii+arr_acc_m(index_i),jj)=b*random2%vector(ii+arr_acc_m(index_i),jj)+a*ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj)
                enddo
            enddo
		enddo      
        !$omp end parallel do
					
        !deallocate (butterflyvector)
                    
    elseif (chara=='T') then
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
        levelm = ceiling_safe(dble(level_butterfly)/2d0)
        
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(num_blocks,1))
        ButterflyVector(0)%num_row=num_blocks
        ButterflyVector(0)%num_col=1
        
		!$omp parallel do default(shared) private(i,mm,ii,jj)
		do i=1, num_blocks
            mm=size(blocks%ButterflyU(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(i,1)%matrix(mm,num_vectors))
            do ii=1, mm
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(i,1)%matrix(ii,jj)=random1%vector(ii+arr_acc_m(i),jj)
                enddo
            enddo
		enddo 
		!$omp end parallel do
        
        do level=0, level_butterfly
            if (level==0) then
                num_groupm=num_blocks
                allocate (ButterflyVector(1)%blocks(num_groupm,1))
                ButterflyVector(1)%num_row=num_groupm
                ButterflyVector(1)%num_col=1
                ! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_groupm
                    rank=size(blocks%ButterflyU(i)%matrix,2)
                    mm=size(blocks%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(i,1)%matrix(rank,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do jj=1, num_vectors
						do ii=1, rank
                            ctemp=0
                            do kk=1, mm
                                ctemp=ctemp+blocks%ButterflyU(i)%matrix(kk,ii)*ButterflyVector(0)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo 
				! !$omp end parallel do				
            else
                num_groupm=blocks%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=blocks%ButterflyKerl(level_butterfly-level+1)%num_col
                if (num_groupm/=1) then
                    allocate (ButterflyVector(level+1)%blocks(int(num_groupm/2),num_groupn))
                    ButterflyVector(level+1)%num_row=int(num_groupm/2)
                    ButterflyVector(level+1)%num_col=num_groupn                    
                    ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn)
					do ij=1,num_groupn*(num_groupm/2)
						j = (ij-1)/(num_groupm/2)+1
						i = (mod(ij-1,(num_groupm/2)) + 1)*2-1	
						index_j=int((j+1)/2)
						index_i=int((i+1)/2)	
						
						mm1=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
						mm2=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
						nn=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
						allocate (ButterflyVector(level+1)%blocks(index_i,j)%matrix(nn,num_vectors))
						
						! ! ! ! ! write(*,*)fnorm(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,mm1,nn),fnorm(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,mm2,nn),&
						! ! ! ! ! &fnorm(ButterflyVector(level)%blocks(i,index_j)%matrix(:,1),mm1,1),fnorm(ButterflyVector(level)%blocks(i+1,index_j)%matrix(:,1),mm2,1),'nimei'
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do ii=1, num_vectors
							do jj=1, nn
								ctemp=0
								do kk=1, mm1
									ctemp=ctemp+ButterflyVector(level)%blocks(i,index_j)%matrix(kk,ii)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
								enddo
								do kk=1, mm2
									ctemp=ctemp+ButterflyVector(level)%blocks(i+1,index_j)%matrix(kk,ii)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(index_i,j)%matrix(jj,ii)=ctemp
							enddo
						enddo
						!$omp end parallel do
						if(level_butterfly-level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							! allocate(matrixtemp(nn,num_vectors))
							! allocate(matrixtemp1(nn,nn))
							! call copymatT_OMP(blocks%ButterflyMiddle(index_i,j)%matrix,matrixtemp1,nn,nn)
							! call gemm_omp(matrixtemp1,ButterflyVector(level+1)%blocks(index_i,j)%matrix,matrixtemp,nn,nn,num_vectors)
							! ButterflyVector(level+1)%blocks(index_i,j)%matrix = matrixtemp
							! deallocate(matrixtemp)
							! deallocate(matrixtemp1)	
							
							call gemmTN_omp(blocks%ButterflyMiddle(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,nn,nn,num_vectors)															
						end if							
                    enddo
					! !$omp end parallel do
                else
					write(*,*)'this part is not optimzed yet'
					stop
                    
					allocate (ButterflyVector(level+1)%blocks(1,num_groupn))
                    ButterflyVector(level+1)%num_row=1
                    ButterflyVector(level+1)%num_col=num_groupn
                    do j=1, num_groupn
                        index_j=int((j+1)/2)
                        nn=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,2)
                        mm=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(1,j)%matrix(nn,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                        do jj=1, num_vectors						
							do ii=1, nn
                                ctemp=0                                
                                do kk=1, mm
                                    ctemp=ctemp+ButterflyVector(level)%blocks(1,index_j)%matrix(kk,jj)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix(kk,ii)
                                enddo
                                ButterflyVector(level+1)%blocks(1,j)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(1,num_blocks))
                ButterflyVector(level+2)%num_row=1
                ButterflyVector(level+2)%num_col=num_blocks
                ! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_blocks
                    nn=size(blocks%ButterflyV(j)%matrix,1)
                    rank=size(blocks%ButterflyV(j)%matrix,2)
                    allocate (ButterflyVector(level+2)%blocks(1,j)%matrix(nn,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do jj=1, num_vectors
						do ii=1, nn
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+ButterflyVector(level+1)%blocks(1,j)%matrix(kk,jj)*blocks%ButterflyV(j)%matrix(ii,kk)
                            enddo
                            ButterflyVector(level+2)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
				! !$omp end parallel do
            endif
			if (level/=0) then
				do j=1, ButterflyVector(level)%num_col
					do i=1, ButterflyVector(level)%num_row
						deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
					enddo
				enddo
			end if	
			if (level==level_butterfly) then			
				do j=1, ButterflyVector(level+1)%num_col
					do i=1, ButterflyVector(level+1)%num_row
						deallocate (ButterflyVector(level+1)%blocks(i,j)%matrix)
					enddo
				enddo			
			end if			
        enddo

  
  
		
	! ! ! do level=0, level_butterfly+2
		! ! ! num_row=ButterflyVector(level)%num_row
		! ! ! num_col=ButterflyVector(level)%num_col
		! ! ! write(*,*)'level:', level
		! ! ! do j=1, num_col
			! ! ! do i=1, num_row
				! ! ! mm = size(ButterflyVector(level)%blocks(i,j)%matrix,1)
				! ! ! nn = size(ButterflyVector(level)%blocks(i,j)%matrix,2)
				! ! ! write(*,*)fnorm(ButterflyVector(level)%blocks(i,j)%matrix(:,1),mm,1)
			! ! ! enddo
		! ! ! enddo
	! ! ! enddo		
		

		!$omp parallel do default(shared) private(index_j,nn,ii,jj)
        do jj=1, num_vectors
			do index_j=1, num_blocks
				nn=size(blocks%ButterflyV(index_j)%matrix,1)
				do ii=1, nn
									
                    random2%vector(ii+arr_acc_n(index_j),jj)=b*random2%vector(ii+arr_acc_n(index_j),jj)+a*ButterflyVector(level_butterfly+2)%blocks(1,index_j)%matrix(ii,jj)
                enddo
            enddo
        enddo 
		!$omp end parallel do
    
    endif
    
    ! !$omp parallel do default(shared) private(level,i,j)
    do level=0, level_butterfly+2
        do j=1, ButterflyVector(level)%num_col
            do i=1, ButterflyVector(level)%num_row
                if(allocated(ButterflyVector(level)%blocks(i,j)%matrix))deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
            enddo
        enddo
        deallocate (ButterflyVector(level)%blocks)
    enddo
    ! !$omp end parallel do
    deallocate (ButterflyVector)   
    
	deallocate(arr_acc_m,arr_acc_n)
	
    return
    
end subroutine butterfly_block_MVP_randomized


subroutine butterfly_block_MVP_randomized_dat(blocks,chara,M,N,Nrnd,random1,random2,a,b)
    
    use MODULE_FILE
	use misc
    implicit none
    
    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: middleflag
	
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    !  complex(kind=8) :: random1(N,Nrnd), random2(M,Nrnd)
        complex(kind=8) :: random1(:,:), random2(:,:)
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'
	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)
	
	middleflag = 0
	if(allocated(blocks%ButterflyMiddle))middleflag=1
	
    level_butterfly=blocks%level_butterfly
    num_blocks=2**level_butterfly	
	allocate(arr_acc_m(num_blocks))
	allocate(arr_acc_n(num_blocks))
	
	k1=0
	k2=0
	do i=1, num_blocks
		arr_acc_n(i) = k1
		arr_acc_m(i) = k2
		nn=size(blocks%ButterflyV(i)%matrix,1)
		k1 =k1 +nn
		mm=size(blocks%ButterflyU(i)%matrix,1)
		k2 =k2 +mm		
	enddo 
	
    num_vectors=Nrnd
    ! write(*,*)num_vectors
	! stop
	
	if(CheckNAN_Butterfly(blocks))then
			write(*,*)'NAN in 0 butterfly_block_MVP_randomized_dat'
			stop		
	end if
	
    if (chara=='N') then
    
		if(isnan(sum(abs(random1(:,1))**2)))then
			write(*,*)'NAN in 1 butterfly_block_MVP_randomized_dat'
			stop
		end if	
	
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
	    levelm = ceiling_safe(dble(level_butterfly)/2d0)        
      
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(1,num_blocks))
        ButterflyVector(0)%num_row=1
        ButterflyVector(0)%num_col=num_blocks
                !  write(*,*)'nima0'
		!$omp parallel do default(shared) private(i,nn,ii,jj)
		do i=1, num_blocks
            nn=size(blocks%ButterflyV(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(1,i)%matrix(nn,num_vectors))
            do ii=1, nn
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(1,i)%matrix(ii,jj)=random1(ii+arr_acc_n(i),jj)
                enddo
            enddo
		enddo 
		!$omp end parallel do
		
        !  write(*,*)'nima1'
        do level=0, level_butterfly
	        !  write(*,*)'nima1',level
            if (level==0) then
                num_groupn=num_blocks
                allocate (ButterflyVector(1)%blocks(1,num_groupn))
                ButterflyVector(1)%num_row=1
                ButterflyVector(1)%num_col=num_groupn
				
                ! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_groupn
				! write(*,*)num_groupn
                    rank=size(blocks%ButterflyV(j)%matrix,2)
                    nn=size(blocks%ButterflyV(j)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(1,j)%matrix(rank,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
					do jj=1, num_vectors                    
						do ii=1, rank
                            ctemp=0
                            do kk=1, nn
                                ctemp=ctemp+blocks%ButterflyV(j)%matrix(kk,ii)*ButterflyVector(0)%blocks(1,j)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
					
					! write(*,*)shape(blocks%ButterflyV(j)%matrix),'da'
					! ! if(is_nan_mat_c(blocks%ButterflyV(j)%matrix,nn,rank))then
						! ! write(*,*)'V'
						! ! stop
					! ! end if					
                enddo  
				! !$omp end parallel do
				
            else
                num_groupm=blocks%ButterflyKerl(level)%num_row
                num_groupn=blocks%ButterflyKerl(level)%num_col
                if (num_groupn/=1) then
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,int(num_groupn/2)))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=int(num_groupn/2)                    
                    
					
                    ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
					do ij=1,num_groupm*(num_groupn/2)
						i = (ij-1)/(num_groupn/2)+1
						j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
						index_i=int((i+1)/2)
						index_j=int((j+1)/2)					

						nn1=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2=size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mm=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						allocate (ButterflyVector(level+1)%blocks(i,index_j)%matrix(mm,num_vectors))
						if(size(ButterflyVector(level)%blocks(index_i,j)%matrix,1)<nn1)then
							write(*,*)blocks%row_group,blocks%col_group,blocks%level_butterfly,level,'nimade'
							stop
						end if
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vectors						
							do ii=1, mm
								ctemp=0
								do kk=1, nn1
									ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j)%matrix(kk,jj)
								enddo
								do kk=1, nn2
									ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j+1)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(i,index_j)%matrix(ii,jj)=ctemp
							enddo
						enddo
						!$omp end parallel do
				! ! if(is_nan_mat_c(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm,nn1) .or. is_nan_mat_c(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,nn2))then
					! ! write(*,*)'kernel'
					! ! stop
				! ! end if
						if(level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							! allocate(matrixtemp(mm,num_vectors))
							call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,mm,mm,num_vectors)
							! ButterflyVector(level+1)%blocks(i,index_j)%matrix = matrixtemp
							! deallocate(matrixtemp)
						end if					
                    enddo
					! !$omp end parallel do
                else
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,1))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=1
                    do i=1, num_groupm
                        index_i=int((i+1)/2)
                        nn=size(blocks%ButterflyKerl(level)%blocks(i,1)%matrix,2)
                        mm=size(blocks%ButterflyKerl(level)%blocks(i,1)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(i,1)%matrix(mm,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vectors	
							do ii=1, mm
                                ctemp=0                       
                                do kk=1, nn
                                    ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(i,1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,1)%matrix(kk,jj)
                                enddo
                                ButterflyVector(level+1)%blocks(i,1)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
					! ! if(is_nan_mat_c(blocks%ButterflyKerl(level)%blocks(i,1)%matrix,mm,nn))then
						! ! write(*,*)'kernel2'
						! ! stop
					! ! end if
													
						
						
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(num_blocks,1))
                ButterflyVector(level+2)%num_row=num_blocks
                ButterflyVector(level+2)%num_col=1
                ! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_blocks
                    rank=size(blocks%ButterflyU(i)%matrix,2)
                    mm=size(blocks%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(level+2)%blocks(i,1)%matrix(mm,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
					do jj=1, num_vectors
						do ii=1, mm
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+blocks%ButterflyU(i)%matrix(ii,kk)*ButterflyVector(level+1)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(level+2)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
					
					! ! if(is_nan_mat_c(blocks%ButterflyU(i)%matrix,mm,rank))then
						! ! write(*,*)'U'
						! ! stop
					! ! end if
					
					
                enddo
				! !$omp end parallel do
            endif
			
			if (level/=0) then
				do j=1, ButterflyVector(level)%num_col
					do i=1, ButterflyVector(level)%num_row
						deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
					enddo
				enddo
			end if	
			if (level==level_butterfly) then			
				do j=1, ButterflyVector(level+1)%num_col
					do i=1, ButterflyVector(level+1)%num_row
						deallocate (ButterflyVector(level+1)%blocks(i,j)%matrix)
					enddo
				enddo			
			end if			
			
        enddo
                                  

        !$omp parallel do default(shared) private(index_i,mm,ii,jj)
		do jj=1, num_vectors
			do index_i=1, num_blocks
				mm=size(blocks%ButterflyU(index_i)%matrix,1)
				do ii=1, mm
				   
						random2(ii+arr_acc_m(index_i),jj)=b*random2(ii+arr_acc_m(index_i),jj)+a*ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj)
						! if(isnan(abs(b*random2(ii+k,jj)+a*ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj))))write(*,*)index_i,ii,k,jj,ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj),random2(ii+k,jj),a,b
					
				enddo
			enddo   
		enddo		
		!$omp end parallel do		
 
		if(isnan(sum(abs(random2(:,1))**2)))then
			write(*,*)'NAN in 2 butterfly_block_MVP_randomized_dat',blocks%row_group,blocks%col_group,blocks%level,blocks%level_butterfly
			stop
		end if
        !deallocate (butterflyvector)
                    
    elseif (chara=='T') then
    
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
        levelm = ceiling_safe(dble(level_butterfly)/2d0)        
        
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(num_blocks,1))
        ButterflyVector(0)%num_row=num_blocks
        ButterflyVector(0)%num_col=1
        
		!$omp parallel do default(shared) private(i,mm,ii,jj)
		do i=1, num_blocks
			mm=size(blocks%ButterflyU(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(i,1)%matrix(mm,num_vectors))
            do ii=1, mm
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(i,1)%matrix(ii,jj)=random1(ii+arr_acc_m(i),jj)
                enddo
            enddo
        enddo 
        !$omp end parallel do		
        
        do level=0, level_butterfly
            if (level==0) then
                num_groupm=num_blocks
                allocate (ButterflyVector(1)%blocks(num_groupm,1))
                ButterflyVector(1)%num_row=num_groupm
                ButterflyVector(1)%num_col=1
                ! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_groupm
                    rank=size(blocks%ButterflyU(i)%matrix,2)
                    mm=size(blocks%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(i,1)%matrix(rank,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do jj=1, num_vectors					
						do ii=1, rank
                            ctemp=0
                            do kk=1, mm
                                ctemp=ctemp+blocks%ButterflyU(i)%matrix(kk,ii)*ButterflyVector(0)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo                   
				! !$omp end parallel do				
            else
                num_groupm=blocks%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=blocks%ButterflyKerl(level_butterfly-level+1)%num_col
                if (num_groupm/=1) then
                    allocate (ButterflyVector(level+1)%blocks(int(num_groupm/2),num_groupn))
                    ButterflyVector(level+1)%num_row=int(num_groupm/2)
                    ButterflyVector(level+1)%num_col=num_groupn                    
                    
					
                    ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn)
					do ij=1,num_groupn*(num_groupm/2)
						j = (ij-1)/(num_groupm/2)+1
						i = (mod(ij-1,(num_groupm/2)) + 1)*2-1	
						index_j=int((j+1)/2)
						index_i=int((i+1)/2)					

						mm1=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
						mm2=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
						nn=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
						allocate (ButterflyVector(level+1)%blocks(index_i,j)%matrix(nn,num_vectors))
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do ii=1, num_vectors
							do jj=1, nn
								ctemp=0
								do kk=1, mm1
									ctemp=ctemp+ButterflyVector(level)%blocks(i,index_j)%matrix(kk,ii)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
								enddo
								do kk=1, mm2
									ctemp=ctemp+ButterflyVector(level)%blocks(i+1,index_j)%matrix(kk,ii)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(index_i,j)%matrix(jj,ii)=ctemp
							enddo
						enddo
						!$omp end parallel do
						
						if(level_butterfly-level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							
							! allocate(matrixtemp(nn,num_vectors))
							! allocate(matrixtemp1(nn,nn))
							! call copymatT_OMP(blocks%ButterflyMiddle(index_i,j)%matrix,matrixtemp1,nn,nn)
							! call gemm_omp(matrixtemp1,ButterflyVector(level+1)%blocks(index_i,j)%matrix,matrixtemp,nn,nn,num_vectors)
							! ButterflyVector(level+1)%blocks(index_i,j)%matrix = matrixtemp
							! deallocate(matrixtemp)
							! deallocate(matrixtemp1)	
							
							call gemmTN_omp(blocks%ButterflyMiddle(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,nn,nn,num_vectors)
							
						end if
                    enddo
					! !$omp end parallel do
                else
                    allocate (ButterflyVector(level+1)%blocks(1,num_groupn))
                    ButterflyVector(level+1)%num_row=1
                    ButterflyVector(level+1)%num_col=num_groupn
                    do j=1, num_groupn
                        index_j=int((j+1)/2)
                        nn=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,2)
                        mm=size(blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(1,j)%matrix(nn,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                        do jj=1, num_vectors
							do ii=1, nn
                                ctemp=0                           
                                do kk=1, mm
                                    ctemp=ctemp+ButterflyVector(level)%blocks(1,index_j)%matrix(kk,jj)*blocks%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix(kk,ii)
                                enddo
                                ButterflyVector(level+1)%blocks(1,j)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(1,num_blocks))
                ButterflyVector(level+2)%num_row=1
                ButterflyVector(level+2)%num_col=num_blocks
                ! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_blocks
                    nn=size(blocks%ButterflyV(j)%matrix,1)
                    rank=size(blocks%ButterflyV(j)%matrix,2)
                    allocate (ButterflyVector(level+2)%blocks(1,j)%matrix(nn,num_vectors))
                    !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do jj=1, num_vectors					
						do ii=1, nn
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+ButterflyVector(level+1)%blocks(1,j)%matrix(kk,jj)*blocks%ButterflyV(j)%matrix(ii,kk)
                            enddo
                            ButterflyVector(level+2)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
				! !$omp end parallel do
            endif
			if (level/=0) then
				do j=1, ButterflyVector(level)%num_col
					do i=1, ButterflyVector(level)%num_row
						deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
					enddo
				enddo
			end if	
			if (level==level_butterfly) then			
				do j=1, ButterflyVector(level+1)%num_col
					do i=1, ButterflyVector(level+1)%num_row
						deallocate (ButterflyVector(level+1)%blocks(i,j)%matrix)
					enddo
				enddo			
			end if			
			
        enddo
        
		!$omp parallel do default(shared) private(index_j,nn,ii,jj)
		do jj=1, num_vectors
			do index_j=1, num_blocks
				nn=size(blocks%ButterflyV(index_j)%matrix,1)
				do ii=1, nn    
									
                    random2(ii+arr_acc_n(index_j),jj)=b*random2(ii+arr_acc_n(index_j),jj)+a*ButterflyVector(level_butterfly+2)%blocks(1,index_j)%matrix(ii,jj)
                enddo
            enddo
        enddo 
		!$omp end parallel do
    
    endif
    
    ! !$omp parallel do default(shared) private(level,i,j)
    do level=0, level_butterfly+2
        do j=1, ButterflyVector(level)%num_col
            do i=1, ButterflyVector(level)%num_row
                if(allocated(ButterflyVector(level)%blocks(i,j)%matrix))deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
            enddo
        enddo
        deallocate (ButterflyVector(level)%blocks)
    enddo
    ! !$omp end parallel do
    deallocate (ButterflyVector)   
    deallocate(arr_acc_m,arr_acc_n)
	
    return
    
end subroutine butterfly_block_MVP_randomized_dat
subroutine Bplus_block_MVP_randomized_dat(bplus,chara,M,N,Nrnd,random1,random2,a,b)
    
    use MODULE_FILE
	use misc
    implicit none
    
    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b,ctemp1,ctemp2
    character chara
	type(matrixblock),pointer::blocks
    integer:: middleflag
	type(blockplus)::bplus
	integer ll,bb
	
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    !  complex(kind=8) :: random1(N,Nrnd), random2(M,Nrnd)
        complex(kind=8) :: random1(:,:), random2(:,:)
        complex(kind=8),allocatable :: Vout(:,:),Vin_loc(:,:),Vout_loc(:,:)
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'
	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)
	
	integer idx_start_m,idx_start_n, idx_start_n_loc,idx_start_m_loc, idx_end_n_loc,idx_end_m_loc
	
	
	ctemp1=1.0d0 ; ctemp2=1.0d0
	
	if (chara=='N') then
		
		allocate(Vout(M,Nrnd))
		Vout = 0
		idx_start_n = basis_group(bplus%col_group)%head
		idx_start_m = basis_group(bplus%row_group)%head
		
		do ll=1,bplus%Lplus
			do bb = 1,bplus%LL(ll)%Nbound
				blocks => bplus%LL(ll)%matrices_block(bb)
				idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
				idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
				idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
				idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
				
				if(blocks%style==1)then
					write(*,*)'style 1 not implemented'
					stop
				else 
					call butterfly_block_MVP_randomized_dat(blocks,'N',idx_end_m_loc-idx_start_m_loc+1,idx_end_n_loc-idx_start_n_loc+1,Nrnd,&
					&random1(idx_start_n_loc:idx_end_n_loc,1:Nrnd),Vout(idx_start_m_loc:idx_end_m_loc,1:Nrnd),ctemp1,ctemp2)					
				end if
			end do
		end do			
		
		random2 = random2*b + Vout*a
		deallocate(Vout)
		
	else 
		! write(*,*)'yo'
		allocate(Vout(N,Nrnd))
		Vout = 0
		idx_start_n = basis_group(bplus%col_group)%head
		idx_start_m = basis_group(bplus%row_group)%head
		
		do ll=1,bplus%Lplus
			do bb = 1,bplus%LL(ll)%Nbound
			! write(*,*)ll,bplus%Lplus,bb,bplus%LL(ll)%Nbound,'dfdfdfddfd'
				blocks => bplus%LL(ll)%matrices_block(bb)
				! write(*,*)'1dfdfdfddfd'
				idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
				idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
				idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
				idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
				
				if(blocks%style==1)then
					write(*,*)'style 1 not implemented'
					stop
				else 
					! write(*,*)'dudud',blocks%level,blocks%level_butterfly,blocks%style
					call butterfly_block_MVP_randomized_dat(blocks,'T',idx_end_m_loc-idx_start_m_loc+1,idx_end_n_loc-idx_start_n_loc+1,Nrnd,&
					&random1(idx_start_m_loc:idx_end_m_loc,1:Nrnd),Vout(idx_start_n_loc:idx_end_n_loc,1:Nrnd),ctemp1,ctemp2)					
					! write(*,*)'dududdddd'
				end if
			end do
		end do			
		! write(*,*)'nimabia'
		
		random2 = random2*b + Vout*a
		deallocate(Vout)
		! write(*,*)'yohi'
	
	end if
	
	
	

end subroutine Bplus_block_MVP_randomized_dat




subroutine Bplus_block_MVP_randomized_dat_partial(bplus,chara,M,N,Nrnd,random1,random2,a,b,level_start,level_end)
    
    use MODULE_FILE
	use misc
    implicit none
    
    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b,ctemp1,ctemp2
    character chara
	type(matrixblock),pointer::blocks
    integer:: middleflag
	type(blockplus)::bplus
	integer ll,bb,level_start,level_end
	
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    !  complex(kind=8) :: random1(N,Nrnd), random2(M,Nrnd)
        complex(kind=8) :: random1(:,:), random2(:,:)
        complex(kind=8),allocatable :: Vout(:,:),Vin_loc(:,:),Vout_loc(:,:)
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'
	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)
	
	integer idx_start_m,idx_start_n, idx_start_n_loc,idx_start_m_loc, idx_end_n_loc,idx_end_m_loc
	
	
	ctemp1=1.0d0 ; ctemp2=1.0d0
	
	if (chara=='N') then
		
		allocate(Vout(M,Nrnd))
		Vout = 0
		idx_start_n = basis_group(bplus%col_group)%head
		idx_start_m = basis_group(bplus%row_group)%head
		
		do ll=level_start,level_end
			do bb = 1,bplus%LL(ll)%Nbound
				blocks => bplus%LL(ll)%matrices_block(bb)
				idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
				idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
				idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
				idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
				
				if(blocks%style==1)then
					write(*,*)'style 1 not implemented'
					stop
				else 
					! write(*,*)'ddd1',ll,bb
					call butterfly_block_MVP_randomized_dat(blocks,'N',idx_end_m_loc-idx_start_m_loc+1,idx_end_n_loc-idx_start_n_loc+1,Nrnd,&
					&random1(idx_start_n_loc:idx_end_n_loc,1:Nrnd),Vout(idx_start_m_loc:idx_end_m_loc,1:Nrnd),ctemp1,ctemp2)					
					! write(*,*)'ddd2'
				end if
			end do
		end do			
		
		random2 = random2*b + Vout*a
		deallocate(Vout)
		
	else 

		allocate(Vout(N,Nrnd))
		Vout = 0
		idx_start_n = basis_group(bplus%col_group)%head
		idx_start_m = basis_group(bplus%row_group)%head
		 
		do ll=level_start,level_end
			do bb = 1,bplus%LL(ll)%Nbound
				blocks => bplus%LL(ll)%matrices_block(bb)
				! write(*,*)'diao',ll,bb,blocks%row_group,blocks%col_group
				idx_start_n_loc = basis_group(blocks%col_group)%head - idx_start_n + 1
				idx_end_n_loc = basis_group(blocks%col_group)%tail - idx_start_n + 1
				idx_start_m_loc = basis_group(blocks%row_group)%head - idx_start_m + 1
				idx_end_m_loc = basis_group(blocks%row_group)%tail	- idx_start_m + 1
				
				if(blocks%style==1)then
					write(*,*)'style 1 not implemented'
					stop
				else 
					call butterfly_block_MVP_randomized_dat(blocks,'T',idx_end_m_loc-idx_start_m_loc+1,idx_end_n_loc-idx_start_n_loc+1,Nrnd,&
					&random1(idx_start_m_loc:idx_end_m_loc,1:Nrnd),Vout(idx_start_n_loc:idx_end_n_loc,1:Nrnd),ctemp1,ctemp2)					
				end if
			end do
		end do			
		
		random2 = random2*b + Vout*a
		deallocate(Vout)

	
	end if
	
	
	

end subroutine Bplus_block_MVP_randomized_dat_partial




subroutine ComputeRandVectIndexArray(blocks,chara,level,IndexArray)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, Ni,Nj, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,level_start,level_end,DimMax,Dimtmp
    integer::IndexArray(:,:,:)
	integer::level_butterfly
	complex(kind=8) ctemp, a, b
    character chara

	type(matrixblock)::blocks
    
	level_butterfly = blocks%level_butterfly
	
    if (chara=='N') then
        num_blocks=2**level_butterfly
		
		Dimtmp = 0
		if(level==0)then
			num_groupn=num_blocks
			do j=1, num_groupn
				nn=size(blocks%ButterflyV(j)%matrix,1)
				Dimtmp = Dimtmp + nn
				IndexArray(1,j,1) = Dimtmp-nn+1
				IndexArray(1,j,2) = Dimtmp
			end do
		else if(level==1)then
			num_groupn=num_blocks
			do j=1, num_groupn
				rank=size(blocks%ButterflyV(j)%matrix,2)
				Dimtmp = Dimtmp + rank
				IndexArray(1,j,1) = Dimtmp-rank+1
				IndexArray(1,j,2) = Dimtmp
			end do
		else if(level<=level_butterfly+1)then
			num_groupm=blocks%ButterflyKerl(level-1)%num_row
			num_groupn=blocks%ButterflyKerl(level-1)%num_col			
			if(num_groupn/=1)then
				do i=1, num_groupm
					index_i=int((i+1)/2)
					do j=1, num_groupn, 2
						index_j=int((j+1)/2)
						mm=size(blocks%ButterflyKerl(level-1)%blocks(i,j)%matrix,1)
						Dimtmp = Dimtmp + mm
						IndexArray(i,index_j,1) = Dimtmp-mm+1
						IndexArray(i,index_j,2) = Dimtmp
					end do
				end do
			else 
				write(*,*)'not implemented'
				stop
			end if
		else if(level==level_butterfly+2)then
			num_groupm = num_blocks
			do i=1, num_groupm
				mm=size(blocks%ButterflyU(i)%matrix,1)
				Dimtmp = Dimtmp + mm
				IndexArray(i,1,1) = Dimtmp-mm+1
				IndexArray(i,1,2) = Dimtmp
			end do	
		end if

			
	else if (chara=='T') then
        num_blocks=2**level_butterfly

		Dimtmp = 0
		if(level==0)then
			num_groupm=num_blocks
			do i=1, num_groupm
				mm=size(blocks%ButterflyU(i)%matrix,1)
				Dimtmp = Dimtmp + mm
				IndexArray(i,1,1) = Dimtmp - mm +1
				IndexArray(i,1,2) = Dimtmp
			end do		
		else if(level==1)then
			num_groupm=num_blocks
			do i=1, num_groupm
				rank=size(blocks%ButterflyU(i)%matrix,2)
				Dimtmp = Dimtmp + rank
				IndexArray(i,1,1) = Dimtmp - rank +1
				IndexArray(i,1,2) = Dimtmp
			end do	
		else if(level<=level_butterfly+1)then
			num_groupm=blocks%ButterflyKerl(level_butterfly-level+2)%num_row
			num_groupn=blocks%ButterflyKerl(level_butterfly-level+2)%num_col
			if (num_groupm/=1) then    
				do j=1, num_groupn
					index_j=int((j+1)/2)
					do i=1, num_groupm, 2
						index_i=int((i+1)/2)
						nn=size(blocks%ButterflyKerl(level_butterfly-level+2)%blocks(i,j)%matrix,2)
						Dimtmp = Dimtmp + nn
						IndexArray(index_i,j,1)=Dimtmp-nn+1
						IndexArray(index_i,j,2)=Dimtmp
					end do
				end do
			else 
				write(*,*)'not implemented'
				stop
			end if 
		else if(level==level_butterfly+2)then
			num_groupn = num_blocks
			do j=1, num_groupn
				nn=size(blocks%ButterflyV(j)%matrix,1)
				Dimtmp = Dimtmp + nn
				IndexArray(1,j,1) = Dimtmp-nn+1					
				IndexArray(1,j,2) = Dimtmp					
			end do		
		end if

	end if		
end subroutine ComputeRandVectIndexArray




subroutine CountMaxIntermidiateVector(blocks,Nmax)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start,Nmax
    integer i, j, ii, jj, level, Ni,Nj, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,level_start,level_end,DimMax,Dimtmp
	integer::level_butterfly
	complex(kind=8) ctemp, a, b
    character chara

	type(matrixblock)::blocks
    
	level_butterfly = blocks%level_butterfly
	
		Nmax = 0
        num_blocks=2**level_butterfly
		
		do level = 0,level_butterfly+2
			Dimtmp = 0
			if(level==0)then
				num_groupn=num_blocks
				do j=1, num_groupn
					nn=size(blocks%ButterflyV(j)%matrix,1)
					Dimtmp = Dimtmp + nn
				end do
			else if(level==1)then
				num_groupn=num_blocks
				do j=1, num_groupn
					rank=size(blocks%ButterflyV(j)%matrix,2)
					Dimtmp = Dimtmp + rank
				end do
			else if(level<=level_butterfly+1)then
				num_groupm=blocks%ButterflyKerl(level-1)%num_row
				num_groupn=blocks%ButterflyKerl(level-1)%num_col			
				if(num_groupn/=1)then
					do i=1, num_groupm
						index_i=int((i+1)/2)
						do j=1, num_groupn, 2
							index_j=int((j+1)/2)
							mm=size(blocks%ButterflyKerl(level-1)%blocks(i,j)%matrix,1)
							Dimtmp = Dimtmp + mm
						end do
					end do
				else 
					write(*,*)'not implemented'
					stop
				end if
			else if(level==level_butterfly+2)then
				num_groupm = num_blocks
				do i=1, num_groupm
					mm=size(blocks%ButterflyU(i)%matrix,1)
					Dimtmp = Dimtmp + mm
				end do	
			end if
			Nmax = max(Nmax,Dimtmp)
		end do
		
		
			
		
end subroutine CountMaxIntermidiateVector



subroutine CountIndexArrayShape(blocks,chara,level,Ni,Nj)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, Ni,Nj, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,level_butterfly
    complex(kind=8) ctemp, a, b
    character chara
    
	type(matrixblock)::blocks


	level_butterfly = blocks%level_butterfly
    
    if (chara=='N') then
        num_blocks=2**level_butterfly
		if(level==0 .or. level==1)then
			Ni=1
			Nj=num_blocks
		else if(level<=level_butterfly+1)then
            num_groupm=blocks%ButterflyKerl(level-1)%num_row
            num_groupn=blocks%ButterflyKerl(level-1)%num_col			
			if(num_groupn/=1)then
				Ni = num_groupm
				Nj = int(num_groupn/2)
			else 
				Ni = num_groupm
				Nj = 1			
			end if
		else if(level==level_butterfly+2)then
			Ni=num_blocks
			Nj=1			
		end if
		
	else if (chara=='T') then
        num_blocks=2**level_butterfly
		! write(*,*)level,level_butterfly+2
		if(level==0 .or. level==1)then
			Ni=num_blocks
			Nj=1
		else if(level<=level_butterfly+1)then
            num_groupm=blocks%ButterflyKerl(level_butterfly-level+2)%num_row
            num_groupn=blocks%ButterflyKerl(level_butterfly-level+2)%num_col			
			if(num_groupm/=1)then
				Ni = int(num_groupm/2)
				Nj = num_groupn
			else 
				Ni = 1
				Nj = num_groupn			
			end if
		else if(level==level_butterfly+2)then
			Ni=1
			Nj=num_blocks			
		end if	
	end if	
end subroutine CountIndexArrayShape
 


subroutine Butterfly_value(mi,nj,blocks,value)

    use MODULE_FILE
    implicit none
    
    integer mm, nn, mi, nj, butterflyB_inuse, groupm_start, groupn_start
    integer i, j, ii, jj, rank, level_butterfly, group_m, group_n, header_mm, header_nn
    integer group, level, flag, mii, njj, rank1, rank2, index_ij, level_blocks
    complex(kind=8) ctemp, value
    
    type(matrixblock) :: blocks
    type(vectorset),allocatable:: vectors_set(:)    
    integer,allocatable :: group_index_mm(:), group_index_nn(:)
    
    group_m=blocks%row_group ! Note: row_group and col_group interchanged here  
    group_n=blocks%col_group
    header_mm=basis_group(group_m)%head
    header_nn=basis_group(group_n)%head
    level_butterfly=blocks%level_butterfly
    level_blocks=blocks%level
    groupm_start=group_m*2**min(level_butterfly,Maxlevel_for_blocks-level_blocks)    !!! modified by Yang Liu
    groupn_start=group_n*2**min(level_butterfly,Maxlevel_for_blocks-level_blocks) 
    
    allocate (group_index_mm(0:level_butterfly),group_index_nn(0:level_butterfly))
    
    flag=0; i=-1
    do while (flag==0)
        i=i+1
        if (basis_group(groupm_start+i)%tail-header_mm+1>=mi) then
            flag=1
        endif
    enddo
    group_index_mm(0)=1+i
    mii=mi-(basis_group(groupm_start+i)%head-header_mm)
    
    flag=0; j=-1
    do while (flag==0)
        j=j+1
        if (basis_group(groupn_start+j)%tail-header_nn+1>=nj) then
            flag=1
        endif
    enddo
    group_index_nn(0)=1+j
    njj=nj-(basis_group(groupn_start+j)%head-header_nn)
    
    if (level_butterfly>0) then
        group_index_mm(1)=group_index_mm(0)
        group_index_nn(1)= group_index_nn(0)
        do level=1, level_butterfly-1
            group_index_mm(level+1)=int((group_index_mm(level)+1)/2)
            group_index_nn(level+1)=int((group_index_nn(level)+1)/2)
        enddo
    endif
    
!     if (group_index_mm(0)/=group_m .or. group_index_nn(0)/=group_n) then
!         write (*,*) 'Butterfly_value_func error1!'
!         pause
!         continue
!     endif
    
!     do level=0, level_butterfly
!         group_index_mm(level)=group_index_mm(level)-group_m*2**level+1
!         group_index_nn(level)=group_index_nn(level)-group_n*2**level+1
!     enddo
    
    allocate (vectors_set(0:level_butterfly))
    do level=0, level_butterfly
        if (level==0) then
            rank=size(blocks%ButterflyV(group_index_nn(0))%matrix,2)
            allocate (vectors_set(level)%vector(rank))
            !!$omp parallel do default(shared) private(i)
            do i=1, rank
                vectors_set(level)%vector(i)=blocks%ButterflyV(group_index_nn(0))%matrix(njj,i)
            enddo
            !!$omp end parallel do
        else
            rank1=size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,2)
            rank2=size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,1)
            allocate (vectors_set(level)%vector(rank2))
            !!$omp parallel do default(shared) private(i,j,ctemp)
            do i=1, rank2
                ctemp=0
                do j=1, rank1
                    ctemp=ctemp+blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix(i,j)*vectors_set(level-1)%vector(j)
                enddo
                vectors_set(level)%vector(i)=ctemp
            enddo
            !!$omp end parallel do
            deallocate (vectors_set(level-1)%vector)
        endif
        if (level==level_butterfly) then
            rank=size(vectors_set(level)%vector,1)
            ctemp=0
            !!$omp parallel do default(shared) private(i) reduction(+:ctemp)
            do i=1, rank
                ctemp=ctemp+blocks%butterflyU(group_index_mm(0))%matrix(mii,i)*vectors_set(level)%vector(i)
            enddo
            !!$omp end parallel do
            value=ctemp
            deallocate (vectors_set(level)%vector)
        endif
    enddo
    deallocate (vectors_set)        
     
    return

end subroutine Butterfly_value

subroutine Butterfly_get_vector(ij,chara,mn,blocks,vector_ij)

    use MODULE_FILE
    implicit none
    
    integer mm, nn, blocks, mi, nj, butterflyB_inuse, groupm_start, groupn_start, ij, mn, num
    integer i, j, ii, jj, rank, level_butterfly, group_m, group_n, header_mm, header_nn, k, kk
    integer group, level, flag, mii, njj, rank1, rank2, index_ij, level_blocks, num_row, num_col
    complex(kind=8) ctemp, value, vector_ij(mn)
    character chara
    
    type(butterfly_multivectors),allocatable:: vectors_set(:)    
    integer,allocatable :: group_index_mm(:), group_index_nn(:)
    
    if (chara=='C') then
    
        !group_m=matrices_block(blocks,0)%col_group
        group_n=matrices_block(blocks,0)%row_group
        !header_mm=basis_group(group_m)%head
        header_nn=basis_group(group_n)%head
        level_butterfly=matrices_block(blocks,0)%level_butterfly
        level_blocks=matrices_block(blocks,0)%level
        !groupm_start=group_m*2**(Maxlevel_for_blocks-level_blocks)
        groupn_start=group_n*2**(Maxlevel_for_blocks-level_blocks)
        
        !allocate (group_index_mm(0:level_butterfly))
        allocate (group_index_nn(0:level_butterfly))
        
        !flag=0; i=-1
        !do while (flag==0)
        !    i=i+1
        !    if (basis_group(groupm_start+i)%tail-header_mm+1>=ij) then
        !        flag=1
        !    endif
        !enddo
        !group_index_mm(0)=1+i
        !mii=ij-(basis_group(groupm_start+i)%head-header_mm)
        
        flag=0; j=-1
        do while (flag==0)
            j=j+1
            if (basis_group(groupn_start+j)%tail-header_nn+1>=ij) then
                flag=1
            endif
        enddo
        group_index_nn(0)=1+j
        njj=ij-(basis_group(groupn_start+j)%head-header_nn)
        
        if (level_butterfly>0) then
            !group_index_mm(1)=group_index_mm(0)
            group_index_nn(1)= group_index_nn(0)
            do level=1, level_butterfly-1
                !group_index_mm(level+1)=int((group_index_mm(level)+1)/2)
                group_index_nn(level+1)=int((group_index_nn(level)+1)/2)
            enddo
        endif
        
        allocate (vectors_set(0:level_butterfly))
        do level=0, level_butterfly
            if (level==0) then
                rank=size(matrices_block(blocks,0)%ButterflyV(group_index_nn(0))%matrix,2)
                vectors_set(level)%num=1
                allocate (vectors_set(level)%blocks(1))
                allocate (vectors_set(level)%blocks(1)%vector(rank))
                !!$omp parallel do default(shared) private(ii)
                do ii=1, rank
                    vectors_set(level)%blocks(1)%vector(ii)=matrices_block(blocks,0)%ButterflyV(group_index_nn(0))%matrix(njj,ii)
                enddo
                !!$omp end parallel do
            else
                !rank1=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,2)
                !rank2=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly-level+1),group_index_nn(level))%matrix,1)
                num_row=matrices_block(blocks,0)%ButterflyKerl(level)%num_row
                vectors_set(level)%num=num_row
                allocate (vectors_set(level)%blocks(num_row))
                if (vectors_set(level)%num==vectors_set(level-1)%num) then
                    do i=1, num_row
                        rank1=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix,2)
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix,1)
                        allocate (vectors_set(level)%blocks(i)%vector(rank2))
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do ii=1, rank2
                            ctemp=0
                            do jj=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix(ii,jj)*vectors_set(level-1)%blocks(i)%vector(jj)
                            enddo
                            vectors_set(level)%blocks(i)%vector(ii)=ctemp
                        enddo
                        !!$omp end parallel do
                        deallocate (vectors_set(level-1)%blocks(i)%vector)
                    enddo
                elseif (vectors_set(level)%num==2*vectors_set(level-1)%num) then
                    do i=1, num_row, 2
                        k=int((i+1)/2)
                        rank1=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix,2)
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix,1)
                        allocate (vectors_set(level)%blocks(i)%vector(rank2))   
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do ii=1, rank2
                            ctemp=0
                            do jj=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,group_index_nn(level))%matrix(ii,jj)*vectors_set(level-1)%blocks(k)%vector(jj)
                            enddo
                            vectors_set(level)%blocks(i)%vector(ii)=ctemp
                        enddo
                        !!$omp end parallel do
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i+1,group_index_nn(level))%matrix,1)
                        allocate (vectors_set(level)%blocks(i+1)%vector(rank2))   
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do ii=1, rank2
                            ctemp=0
                            do jj=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i+1,group_index_nn(level))%matrix(ii,jj)*vectors_set(level-1)%blocks(k)%vector(jj)
                            enddo
                            vectors_set(level)%blocks(i+1)%vector(ii)=ctemp
                        enddo
                        !!$omp end parallel do
                        deallocate (vectors_set(level-1)%blocks(k)%vector)
                    enddo
                endif                    
            endif
            if (level==level_butterfly) then
                num=vectors_set(level)%num
                kk=0
                do i=1, num
                    mm=size(matrices_block(blocks,0)%ButterflyU(i)%matrix,1)            
                    rank=size(matrices_block(blocks,0)%ButterflyU(i)%matrix,2)
                    !!$omp parallel do default(shared) private(ii,jj,ctemp)
                    do ii=1, mm
                        ctemp=0
                        do jj=1, rank
                            ctemp=ctemp+matrices_block(blocks,0)%ButterflyU(i)%matrix(ii,jj)*vectors_set(level)%blocks(i)%vector(jj)
                        enddo
                        vector_ij(ii+kk)=ctemp
                    enddo
                    !!$omp end parallel do
                    deallocate (vectors_set(level)%blocks(i)%vector)
                    kk=kk+mm
                enddo
            endif
        enddo
        
        do level=0, level_butterfly
            deallocate (vectors_set(level)%blocks)
        enddo
        deallocate (vectors_set)
        
    elseif (chara=='R') then
        
        group_m=matrices_block(blocks,0)%col_group
        !group_n=matrices_block(blocks,0)%row_group
        header_mm=basis_group(group_m)%head
        !header_nn=basis_group(group_n)%head
        level_butterfly=matrices_block(blocks,0)%level_butterfly
        level_blocks=matrices_block(blocks,0)%level
        groupm_start=group_m*2**(Maxlevel_for_blocks-level_blocks)
        !groupn_start=group_n*2**(Maxlevel_for_blocks-level_blocks)
        
        allocate (group_index_mm(0:level_butterfly))
        !allocate (group_index_nn(0:level_butterfly))
        
        flag=0; i=-1
        do while (flag==0)
            i=i+1
            if (basis_group(groupm_start+i)%tail-header_mm+1>=ij) then
                flag=1
            endif
        enddo
        group_index_mm(0)=1+i
        mii=ij-(basis_group(groupm_start+i)%head-header_mm)
        
        !flag=0; j=-1
        !do while (flag==0)
        !    j=j+1
        !    if (basis_group(groupn_start+j)%tail-header_nn+1>=ij) then
        !        flag=1
        !    endif
        !enddo
        !group_index_nn(0)=1+j
        !njj=ij-(basis_group(groupn_start+j)%head-header_nn)
        
        if (level_butterfly>0) then
            group_index_mm(1)=group_index_mm(0)
            !group_index_nn(1)= group_index_nn(0)
            do level=1, level_butterfly-1
                group_index_mm(level+1)=int((group_index_mm(level)+1)/2)
                !group_index_nn(level+1)=int((group_index_nn(level)+1)/2)
            enddo
        endif
        
        allocate (vectors_set(0:level_butterfly))
        do level=0, level_butterfly
            if (level==0) then
                rank=size(matrices_block(blocks,0)%ButterflyU(group_index_mm(0))%matrix,2)
                vectors_set(level)%num=1
                allocate (vectors_set(level)%blocks(1))
                allocate (vectors_set(level)%blocks(1)%vector(rank))
                !!$omp parallel do default(shared) private(jj)
                do jj=1, rank
                    vectors_set(level)%blocks(1)%vector(jj)=matrices_block(blocks,0)%ButterflyU(group_index_mm(0))%matrix(mii,jj)
                enddo
                !!$omp end parallel do
            else
                num_col=matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%num_col
                vectors_set(level)%num=num_col
                allocate (vectors_set(level)%blocks(num_col))
                if (vectors_set(level)%num==vectors_set(level-1)%num) then
                    do j=1, num_col
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix,2)
                        rank1=size(matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix,1)
                        allocate (vectors_set(level)%blocks(j)%vector(rank2))
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do jj=1, rank2
                            ctemp=0
                            do ii=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix(ii,jj)*vectors_set(level-1)%blocks(j)%vector(ii)
                            enddo
                            vectors_set(level)%blocks(j)%vector(jj)=ctemp
                        enddo
                        !!$omp end parallel do
                        deallocate (vectors_set(level-1)%blocks(j)%vector)
                    enddo
                elseif (vectors_set(level)%num==2*vectors_set(level-1)%num) then
                    do j=1, num_col, 2
                        k=int((j+1)/2)
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix,2)
                        rank1=size(matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix,1)
                        allocate (vectors_set(level)%blocks(j)%vector(rank2))   
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do jj=1, rank2
                            ctemp=0
                            do ii=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j)%matrix(ii,jj)*vectors_set(level-1)%blocks(k)%vector(ii)
                            enddo
                            vectors_set(level)%blocks(j)%vector(jj)=ctemp
                        enddo
                        !!$omp end parallel do
                        rank2=size(matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j+1)%matrix,2)
                        allocate (vectors_set(level)%blocks(j+1)%vector(rank2))   
                        !!$omp parallel do default(shared) private(ii,jj,ctemp)
                        do jj=1, rank2
                            ctemp=0
                            do ii=1, rank1
                                ctemp=ctemp+matrices_block(blocks,0)%ButterflyKerl(level_butterfly-level+1)%blocks(group_index_mm(level),j+1)%matrix(ii,jj)*vectors_set(level-1)%blocks(k)%vector(ii)
                            enddo
                            vectors_set(level)%blocks(j+1)%vector(jj)=ctemp
                        enddo
                        !!$omp end parallel do
                        deallocate (vectors_set(level-1)%blocks(k)%vector)
                    enddo
                endif                    
            endif
            if (level==level_butterfly) then
                num=vectors_set(level)%num
                kk=0
                do j=1, num
                    nn=size(matrices_block(blocks,0)%ButterflyV(j)%matrix,1)            
                    rank=size(matrices_block(blocks,0)%ButterflyV(j)%matrix,2)
                    !!$omp parallel do default(shared) private(ii,jj,ctemp)
                    do jj=1, nn
                        ctemp=0
                        do ii=1, rank
                            ctemp=ctemp+matrices_block(blocks,0)%ButterflyV(j)%matrix(jj,ii)*vectors_set(level)%blocks(j)%vector(ii)
                        enddo
                        vector_ij(jj+kk)=ctemp
                    enddo
                    !!$omp end parallel do
                    deallocate (vectors_set(level)%blocks(j)%vector)
                    kk=kk+nn
                enddo
            endif
        enddo
        
        do level=0, level_butterfly
            deallocate (vectors_set(level)%blocks)
        enddo
        deallocate (vectors_set)
       
    endif    
     
    return

end subroutine Butterfly_get_vector

recursive subroutine element_extraction(mi,nj,blocks,value)
    
    use MODULE_FILE
    implicit none
    
    integer mi, nj, blocks, i, j, k, ii, jj, kk, nested_num
    integer group_m, group_n
    integer level, level_butterfly, mm, nn
    complex(kind=8) value, ctemp
    
    if (matrices_block(blocks,0)%style==3) then
        group_m=matrices_block(4*blocks+1,0)%col_group
        group_n=matrices_block(4*blocks+1,0)%row_group
        mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
        nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
        if (mi<=mm .and. nj<=nn) then
            call element_extraction(mi,nj,4*blocks+1,ctemp)
            value=ctemp
        elseif (mi>mm .and. nj<=nn) then
            call element_extraction(mi-mm,nj,4*blocks+2,ctemp)
            value=ctemp
        elseif (mi<=mm .and. nj>nn) then
            call element_extraction(mi,nj-nn,4*blocks+3,ctemp)
            value=ctemp
        elseif (mi>mm .and. nj>nn) then
            call element_extraction(mi-mm,nj-nn,4*blocks+4,ctemp)
            value=ctemp
        endif
    elseif (matrices_block(blocks,0)%style==4) then
        nested_num=matrices_block(blocks,0)%nested_num
        if (nested_num==1) then
            call element_extraction(mi,nj,int((blocks-1)/4),ctemp)
            value=ctemp
        elseif (nested_num==2) then
            group_m=matrices_block(blocks-1,0)%col_group
            mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
            call element_extraction(mi+mm,nj,int((blocks-1)/4),ctemp)
            value=ctemp
        elseif (nested_num==3) then
            group_n=matrices_block(blocks-2,0)%col_group
            nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
            call element_extraction(mi,nj+nn,int((blocks-1)/4),ctemp)
            value=ctemp
        elseif (nested_num==4) then
            group_m=matrices_block(blocks-3,0)%col_group
            group_n=matrices_block(blocks-3,0)%row_group
            mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
            nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
            call element_extraction(mi+mm,nj+nn,int((blocks-1)/4),ctemp)
            value=ctemp
        endif
    elseif (matrices_block(blocks,0)%style==2) then
        call Butterfly_value(mi,nj,matrices_block(blocks,0),ctemp)
        value=ctemp
    elseif (matrices_block(blocks,0)%style==1) then
        value=matrices_block(blocks,0)%fullmat(mi,nj)
    endif
    
    return
    
end subroutine element_extraction





subroutine fullmat_block_MVP_randomized_dat(blocks,chara,M,N,random1,random2,a,b)		
	use MODULE_FILE
    ! use lapack95
    ! use blas95	
	use misc
    implicit none
    
    integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2
    complex(kind=8) ctemp, a, b
    character chara
	type(matrixblock)::blocks
	integer M,N
    complex(kind=8) :: random1(M,N), random2(M,N)
	complex(kind=8):: al,be
	complex(kind=8),allocatable :: random2tmp(:,:)
	allocate(random2tmp(M,N))
	
	al=1d0
	be=0d0
	
	num_vectors=size(random1,2)

	
	random2tmp = random2
	call assert(size(blocks%fullmat,1)==size(blocks%fullmat,2) ,'M not square')
	if(size(blocks%fullmat,1)/=M)write(*,*)M,N,shape(blocks%fullmat),blocks%row_group,blocks%col_group,'niao'
	call assert(size(blocks%fullmat,1)==M,'M not equal fullmat dim')
	
	if (chara=='N') then
        group_m=blocks%row_group  ! Note: row_group and col_group interchanged here   
        group_n=blocks%col_group
		call assert(group_m==group_n,'fullmat not square')
        level_blocks=blocks%level
		! write(*,*)shape(blocks%fullmat),shape(random1),shape(random2),num_vectors
		! call gemmf90(blocks%fullmat, random1, random2,'N','N',al,be)                   
		call gemm_omp(blocks%fullmat, random1, random2,M,M,N)                   
    elseif (chara=='T') then
        group_m=blocks%row_group  ! Note: row_group and col_group interchanged here   
        group_n=blocks%col_group
		call assert(group_m==group_n,'fullmat not square')
        level_blocks=blocks%level
		! call gemmf90(blocks%fullmat, random1, random2,'T','N',al,be)    
		call gemmTN_omp(blocks%fullmat, random1, random2,M,M,N)
	end if
	
	random2 = a*random2+ b*random2tmp
	! write(*,*)'wo cao ni ma'
	deallocate(random2tmp)
end subroutine fullmat_block_MVP_randomized_dat


subroutine get_butterfly_minmaxrank(block_i)
use MODULE_FILE
use misc
implicit none 
type(matrixblock)::block_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly	


block_i%rankmin = 100000
block_i%rankmax = -100000


level_butterfly = block_i%level_butterfly
num_blocks=2**level_butterfly

do level=0, level_butterfly
	index_ij=0
	do index_i=1, 2**level
		do index_j=1, 2**(level_butterfly-level)
			index_ij=index_ij+1
			if (level==0) then
				nn=size(block_i%ButterflyV(index_ij)%matrix,1)
				rank=size(block_i%ButterflyV(index_ij)%matrix,2)
				block_i%rankmin = min(block_i%rankmin,rank)
				block_i%rankmax = max(block_i%rankmax,rank)
			else                    
				nn=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,2)
				rank=size(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
				block_i%rankmin = min(block_i%rankmin,rank)
				block_i%rankmax = max(block_i%rankmax,rank)				
			endif
			if (level==level_butterfly) then
				mm=size(block_i%ButterflyU(index_ij)%matrix,1)
				rank=size(block_i%ButterflyU(index_ij)%matrix,2)
				block_i%rankmin = min(block_i%rankmin,rank)
				block_i%rankmax = max(block_i%rankmax,rank)		
			endif
		enddo
	enddo
enddo

end subroutine get_butterfly_minmaxrank







subroutine Butterfly_sym2asym(blocks)
    
    use MODULE_FILE
	use misc
	! use lapack95
    ! use blas95	
    implicit none
    
    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: middleflag,dimension_n,num_row,num_col,mn_min
	
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'

    real*8, allocatable :: Singular(:)
    complex(kind=8), allocatable :: UU(:,:),VV(:,:)

	
	if(allocated(blocks%ButterflyMiddle))then

		
		
        group_m=blocks%row_group ! Note: row_group and col_group interchanged here   
        group_n=blocks%col_group
        level_butterfly=blocks%level_butterfly
        num_blocks=2**level_butterfly
	    levelm = ceiling_safe(dble(level_butterfly)/2d0)       
		
		call assert(level_butterfly>=2,'level_butterfly not correct')

		level = levelm
		num_groupm=blocks%ButterflyKerl(level)%num_row
		num_groupn=blocks%ButterflyKerl(level)%num_col

			
		! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
		do ij=1,num_groupm*(num_groupn/2)
			i = (ij-1)/(num_groupn/2)+1
			j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
			index_i=int((i+1)/2)
			index_j=int((j+1)/2)					

			nn1=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
			nn2=size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
			mm=size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)

			call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm,mm,nn1)
			call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,mm,nn2)			
			deallocate(blocks%ButterflyMiddle(i,index_j)%matrix)
		enddo
		! !$omp end parallel do

		deallocate(blocks%ButterflyMiddle)
		
		do level=0, levelm-1
			if(level==0)then

				iijj=0
				do j=1, num_blocks
					iijj = iijj + 1
					dimension_n=size(blocks%butterflyV(j)%matrix,1)
					rank = size(blocks%butterflyV(j)%matrix,2)
					mn_min = min(dimension_n,rank)
					
					allocate(matrixtemp(rank,dimension_n))
					allocate(UU(rank,mn_min))
					allocate(VV(mn_min,dimension_n))
					allocate(Singular(mn_min))
					
					call copymatT_omp(blocks%butterflyV(j)%matrix,matrixtemp,dimension_n,rank)
					
					call gesvd_robust(matrixtemp,Singular,UU,VV,rank,dimension_n,mn_min)
					do ii=1,mn_min
						UU(:,ii) = UU(:,ii)*Singular(ii)
					end do
					
					deallocate(blocks%butterflyV(j)%matrix)
					allocate(blocks%butterflyV(j)%matrix(dimension_n,mn_min))
					call copymatT_omp(VV,blocks%butterflyV(j)%matrix,mn_min,dimension_n)
						
					index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
					index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))
					mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
					allocate(matrixtemp1(mm1,mn_min))
					call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,rank,mn_min)	
					deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
					allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
					blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
					deallocate(matrixtemp1)
					
					mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
					allocate(matrixtemp1(mm2,mn_min))
					call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,rank,mn_min)	
					deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
					allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
					blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
					deallocate(matrixtemp1)
					
					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)

				enddo
			else 
				num_row=blocks%ButterflyKerl(level)%num_row
				num_col=blocks%ButterflyKerl(level)%num_col
				
				iijj=0	
				do i=1,	num_row
					do j =1, num_col, 2
						iijj = iijj + 1
						rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						nn1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2 = size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mn_min = min(nn1+nn2,rank)
						
						allocate(matrixtemp(rank,nn1+nn2))
						allocate(UU(rank,mn_min))
						allocate(VV(mn_min,nn1+nn2))
						allocate(Singular(mn_min))
						
						call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
						call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
						
						call gesvd_robust(matrixtemp,Singular,UU,VV,rank,nn1+nn2,mn_min)
						do ii=1,mn_min
							UU(:,ii) = UU(:,ii)*Singular(ii)
						end do								

						deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
						allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mn_min,nn1))						
						call copymatN_omp(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
						deallocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix)
						allocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(mn_min,nn2))												
						call copymatN_omp(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
												
						
						index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
						index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))
						

						mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,rank,mn_min)	
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)						
						
						mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
						allocate(matrixtemp1(mm2,mn_min))
						call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,rank,mn_min)	
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)						
						
						deallocate(matrixtemp)
						deallocate(UU)
						deallocate(VV)
						deallocate(Singular)
					
					end do						
				end do		
			end if
		end do		

	end if        
    
end subroutine Butterfly_sym2asym




subroutine Butterfly_MoveSingulartoLeft(blocks)
    
    use MODULE_FILE
	use misc
	! use lapack95
    ! use blas95	
    implicit none
    
    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: middleflag,dimension_n,dimension_m,num_row,num_col,mn_min
	
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'

    real*8, allocatable :: Singular(:)
    complex(kind=8), allocatable :: UU(:,:),VV(:,:)

	
	group_m=blocks%row_group ! Note: row_group and col_group interchanged here   
	group_n=blocks%col_group
	level_butterfly=blocks%level_butterfly
	num_blocks=2**level_butterfly

	do level=0, level_butterfly
		if(level==0)then
			iijj=0
			do j=1, num_blocks
				iijj = iijj + 1
				dimension_n=size(blocks%butterflyV(j)%matrix,1)
				rank = size(blocks%butterflyV(j)%matrix,2)
				mn_min = min(dimension_n,rank)
				
				allocate(matrixtemp(rank,dimension_n))
				allocate(UU(rank,mn_min))
				allocate(VV(mn_min,dimension_n))
				allocate(Singular(mn_min))
				
				call copymatT_omp(blocks%butterflyV(j)%matrix,matrixtemp,dimension_n,rank)
				call assert(.not. isnan(fnorm(matrixtemp,rank,dimension_n)),'matrixtemp NAN at 3')
				
				call gesvd_robust(matrixtemp,Singular,UU,VV,rank,dimension_n,mn_min)
				call assert(.not. isnan(sum(Singular)),'Singular NAN at 3')
				
				do ii=1,mn_min
					UU(:,ii) = UU(:,ii)*Singular(ii)
				end do
				
				
				deallocate(blocks%butterflyV(j)%matrix)
				allocate(blocks%butterflyV(j)%matrix(dimension_n,mn_min))				
				call copymatT_omp(VV,blocks%butterflyV(j)%matrix,mn_min,dimension_n)
					
	
				index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
				index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))

				mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
				allocate(matrixtemp1(mm1,mn_min))
				call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,rank,mn_min)	
				deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
				allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
				blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
				deallocate(matrixtemp1)						
				
				mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
				allocate(matrixtemp1(mm2,mn_min))
				call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,rank,mn_min)	
				deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
				allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
				blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
				deallocate(matrixtemp1)						
				
				deallocate(matrixtemp)
				deallocate(UU)
				deallocate(VV)
				deallocate(Singular)

			enddo
		else 
			num_row=blocks%ButterflyKerl(level)%num_row
			num_col=blocks%ButterflyKerl(level)%num_col
			
			iijj=0	
			do i=1,	num_row
				do j =1, num_col, 2
					iijj = iijj + 1
					rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
					nn1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
					nn2 = size(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
					mn_min = min(nn1+nn2,rank)
					
					allocate(matrixtemp(rank,nn1+nn2))
					allocate(UU(rank,mn_min))
					allocate(VV(mn_min,nn1+nn2))
					allocate(Singular(mn_min))
					
					call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
					call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
					
					call assert(.not. isnan(fnorm(matrixtemp,rank,nn1+nn2)),'matrixtemp NAN at 4')
					call gesvd_robust(matrixtemp,Singular,UU,VV,rank,nn1+nn2,mn_min)
					call assert(.not. isnan(sum(Singular)),'Singular NAN at 4')
					
					do ii=1,mn_min
						UU(:,ii) = UU(:,ii)*Singular(ii)
					end do								
					
					deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mn_min,nn1))						
					call copymatN_omp(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
					deallocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix(mn_min,nn2))												
					call copymatN_omp(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
											
					if(level/=level_butterfly)then
						index_j = mod(iijj-1,blocks%ButterflyKerl(level+1)%num_col)+1
						index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level+1)%num_col))					

						mm1=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,rank,mn_min)	
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix(mm1,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)						
						
						mm2=size(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,1)
						allocate(matrixtemp1(mm2,mn_min))
						call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,rank,mn_min)	
						deallocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix)
						allocate(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix(mm2,mn_min))
						blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix = matrixtemp1
						deallocate(matrixtemp1)															
					else 
						mm1 = size(blocks%ButterflyU(i)%matrix,1)
						allocate(matrixtemp1(mm1,mn_min))
						call gemm_omp(blocks%ButterflyU(i)%matrix,UU,matrixtemp1,mm1,rank,mn_min)
						deallocate(blocks%ButterflyU(i)%matrix)
						allocate(blocks%ButterflyU(i)%matrix(mm1,mn_min))
						blocks%ButterflyU(i)%matrix = matrixtemp1
						deallocate(matrixtemp1)									
					end if

					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)
				
				end do						
			end do		
		end if
	end do		
       
    
end subroutine Butterfly_MoveSingulartoLeft





subroutine Butterfly_MoveSingulartoRight(blocks)
    
    use MODULE_FILE
	use misc
	! use lapack95
    ! use blas95	
    implicit none
    
    integer M,N, Nrnd,group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b
    character chara
	type(matrixblock)::blocks
    integer:: middleflag,dimension_n,dimension_m,num_row,num_col,mn_min
	
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'

    real*8, allocatable :: Singular(:)
    complex(kind=8), allocatable :: UU(:,:),VV(:,:)

	
	group_m=blocks%row_group ! Note: row_group and col_group interchanged here   
	group_n=blocks%col_group
	level_butterfly=blocks%level_butterfly
	num_blocks=2**level_butterfly

	do level=level_butterfly+1, 1,-1
		if(level==level_butterfly+1)then
			iijj=0
			do i=1, num_blocks
				iijj = iijj + 1
				dimension_m=size(blocks%butterflyU(i)%matrix,1)
				rank = size(blocks%butterflyU(i)%matrix,2)
				mn_min = min(dimension_m,rank)
				
				allocate(matrixtemp(dimension_m,rank))
				allocate(UU(dimension_m,mn_min))
				allocate(VV(mn_min,rank))
				allocate(Singular(mn_min))
				
				call copymatN_omp(blocks%butterflyU(i)%matrix,matrixtemp,dimension_m,rank)
				call assert(.not. isnan(fnorm(matrixtemp,dimension_m,rank)),'matrixtemp NAN at 1')
				
				call gesvd_robust(matrixtemp,Singular,UU,VV,dimension_m,rank,mn_min)
				call assert(.not. isnan(sum(Singular)),'Singular NAN at 1')
				
				do ii=1,mn_min
					VV(ii,:) = VV(ii,:)*Singular(ii)
				end do
				
				deallocate(blocks%butterflyU(i)%matrix)
				allocate(blocks%butterflyU(i)%matrix(dimension_m,mn_min))				
				call copymatN_omp(UU,blocks%butterflyU(i)%matrix,dimension_m,mn_min)
					
				index_i = mod(iijj-1,blocks%ButterflyKerl(level-1)%num_row)+1
				index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level-1)%num_row))
				
				nn1=size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,2)
				allocate(matrixtemp1(mn_min,nn1))
				call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,rank,nn1)			
				deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix)
				allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix(mn_min,nn1))
				blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix = matrixtemp1
				deallocate(matrixtemp1)

				nn2=size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,2)
				allocate(matrixtemp1(mn_min,nn2))
				call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,rank,nn2)			
				deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix)
				allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix(mn_min,nn2))
				blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix = matrixtemp1
				deallocate(matrixtemp1)
				
				deallocate(matrixtemp)
				deallocate(UU)
				deallocate(VV)
				deallocate(Singular)  

			enddo
		else 
			num_row=blocks%ButterflyKerl(level)%num_row
			num_col=blocks%ButterflyKerl(level)%num_col
			
			iijj=0	
			do j=1,	num_col
				do i =1, num_row, 2
					iijj = iijj + 1
					rank = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,2)
					
					mm1 = size(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,1)
					mm2 = size(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,1)
					mn_min = min(mm1+mm2,rank)
					
					allocate(matrixtemp(mm1+mm2,rank))
					allocate(UU(mm1+mm2,mn_min))
					allocate(VV(mn_min,rank))
					allocate(Singular(mn_min))
					
					call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:mm1,1:rank),mm1,rank)
					call copymatN_omp(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,matrixtemp(1+mm1:mm2+mm1,1:rank),mm2,rank)
					call assert(.not. isnan(fnorm(matrixtemp,mm1+mm2,rank)),'matrixtemp NAN at 2')
					
					call gesvd_robust(matrixtemp,Singular,UU,VV,mm1+mm2,rank,mn_min)
					! if(isnan(sum(Singular)).and. mm1+mm2<rank)then
						! write(*,*)mm1+mm2,rank,mm1+mm2>=rank,'rank too large?'
					! end if
					
					! call assert(.not. isnan(sum(Singular)),'Singular NAN at 2')
					if(isnan(sum(Singular)))then
						write(*,*)'Singular NAN at 2',mm1+mm2,rank
						do ii=1,mm1+mm2
							do jj=1,rank
								write(777,*)dble(matrixtemp(ii,jj)),aimag(matrixtemp(ii,jj)),abs(matrixtemp(ii,jj))
							end do
						end do
						stop
					end if
					
					
					
					do ii=1,mn_min
						VV(ii,:) = VV(ii,:)*Singular(ii)
					end do								
					
					deallocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i,j)%matrix(mm1,mn_min))
					call copymatN_omp(UU(1:mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm1,mn_min)
					deallocate(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix)
					allocate(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix(mm2,mn_min))					
					call copymatN_omp(UU(1+mm1:mm2+mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,mm2,mn_min)
											
											
					if(level/=1)then
						index_i = mod(iijj-1,blocks%ButterflyKerl(level-1)%num_row)+1
						index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level-1)%num_row))					
						nn1 = size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,2)
						
						allocate(matrixtemp1(mn_min,nn1))
						call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,rank,nn1)
						deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix)
						allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix(mn_min,nn1))
						blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix = matrixtemp1
						deallocate(matrixtemp1)
						
						nn2 = size(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,2)
						allocate(matrixtemp1(mn_min,nn2))
						call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,rank,nn2)
						deallocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix)
						allocate(blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix(mn_min,nn2))
						blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix = matrixtemp1
						deallocate(matrixtemp1)						
					else 
						nn1 = size(blocks%ButterflyV(j)%matrix,1)
						allocate(matrixtemp1(nn1,mn_min))
						call gemmNT_omp(blocks%ButterflyV(j)%matrix,VV,matrixtemp1,nn1,rank,mn_min)
						deallocate(blocks%ButterflyV(j)%matrix)
						allocate(blocks%ButterflyV(j)%matrix(nn1,mn_min))
						blocks%ButterflyV(j)%matrix = matrixtemp1
						deallocate(matrixtemp1)
					end if

					deallocate(matrixtemp)
					deallocate(UU)
					deallocate(VV)
					deallocate(Singular)
				
				end do						
			end do		
		end if
	end do		
       
    
end subroutine Butterfly_MoveSingulartoRight




end module Utilities
