PROGRAM RANDOMIZED_BUTTERFLY_CONSTRUCTION

    use MODULE_FILE
	use omp_lib
	use Utilities
	
    implicit none
    
	real*8:: n1,n2,error
	integer:: i,j,k,level,num_row,num_col,mm,nn, dimension_m,dimension_n,dimension_rank,inner_rank,num_blocks,level_butterfly,threads_num,num_vect
	type(matrixblock),pointer::block
	complex(kind=8) ctemp, ctemp1, ctemp2
	complex(kind=8),allocatable:: vec_in(:,:),vec_out(:,:)
	
	time_indexarray = 0
	time_leastsquare = 0
	time_buttermul = 0
	time_buttermulinv = 0
	
    ! open (90,file='input.txt')
	! read (90,*)
	! read (90,*) dimension_m
	! read (90,*) dimension_n
	! read (90,*) dimension_rank
	! read (90,*) inner_rank
	! read (90,*) rank_reconstruction   
    ! read (90,*) level_butterfly
	! read (90,*) iter_tolerance
	! read (90,*) iter_max
	! read (90,*) relax_lamda
	! read (90,*) SolvingMethod
    ! read (90,*) threads_num
    ! read (90,*) LS_tolerance
    ! read (90,*) SVD_tolerance
    ! close (90)
	
	dimension_m = 8
	dimension_n = 8
	dimension_rank = 8
	inner_rank = 8
	level_butterfly = 2
	threads_num = 1 
	num_vect = 100
	
    call OMP_set_num_threads(threads_num)
    call OMP_set_dynamic(.true.)

	allocate(block)
	call create_butterfly_simple(block,level_butterfly,dimension_m,dimension_n,dimension_rank,inner_rank)

	num_blocks=2**level_butterfly	
	mm = num_blocks*dimension_m
	nn = num_blocks*dimension_n		
	ctemp1=1.0d0 ; ctemp2=0.0d0
	
	allocate (vec_in(nn,num_vect))
	allocate (vec_out(mm,num_vect))
	
	vec_in=1d0
	
	call butterfly_block_MVP_randomized_dat(block,'N',mm,nn,num_vect,vec_in,vec_out,ctemp1,ctemp2)
	
	write(*,*)'output fnorm: ',fnorm(vec_out,mm,1)
	
	call delete_blocks(block)
	deallocate(block)
	deallocate(vec_in)
	deallocate(vec_out)

END PROGRAM RANDOMIZED_BUTTERFLY_CONSTRUCTION




subroutine create_butterfly_simple(block,level_butterfly,dimension_m,dimension_n,dimension_rank,inner_rank)
	use misc
    use MODULE_FILE
	! ! use lapack95
	! ! use blas95
    implicit none
    
    integer level_c,rowblock
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_rank, inner_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j,blocks_idx
	type(matrixblock)::block
	integer seed_myid(12)
	integer times(8)
	
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	call RANDOM_SEED(PUT=seed_myid)
	
	
	block%level = 0
	block%col_group = 1
	block%row_group = 1
	! block%nested_num = 1
	block%style = 2
	! block%data_type = 1
	block%level_butterfly = level_butterfly
	block%rankmax = dimension_rank
	block%rankmin = dimension_rank
	
	num_blocks=2**level_butterfly	

	mm = num_blocks*dimension_m
	nn = num_blocks*dimension_n	
	
    allocate (block%ButterflyU(2**level_butterfly))
    allocate (block%ButterflyV(2**level_butterfly))

    do blocks=1, num_blocks
        allocate (block%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))
		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,inner_rank,matrixtemp1,1)
        do j=1, dimension_rank
            do i=1, dimension_m
				block%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		deallocate(matrixtemp1)

		
        allocate (block%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
		
		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,inner_rank,matrixtemp1,1)
        do j=1, dimension_rank
            do i=1, dimension_n
				block%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do	
		deallocate(matrixtemp1)

    enddo
		
    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*inner_rank))
        allocate (block%ButterflyKerl(level_butterfly))
        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            block%ButterflyKerl(level)%num_row=num_row
            block%ButterflyKerl(level)%num_col=num_col
            allocate (block%ButterflyKerl(level)%blocks(num_row,num_col))
            do j=1, num_col, 2
                index_j=int((j+1)/2)
                do i=1, num_row, 2
                    index_i=int((i+1)/2)
                    allocate (block%ButterflyKerl(level)%blocks(i,j)%matrix(dimension_rank,dimension_rank))
                    allocate (block%ButterflyKerl(level)%blocks(i+1,j)%matrix(dimension_rank,dimension_rank))
                    allocate (block%ButterflyKerl(level)%blocks(i,j+1)%matrix(dimension_rank,dimension_rank))
                    allocate (block%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(dimension_rank,dimension_rank))

				    call RandomMat(2*dimension_rank,2*dimension_rank,2*inner_rank,matrixtemp1,1)	
					
                    !$omp parallel do default(shared) private(ii,jj)
                    do jj=1, dimension_rank
                        do ii=1, dimension_rank
                            block%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=matrixtemp1(ii,jj)
                            block%ButterflyKerl(level)%blocks(i+1,j)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj)
                            block%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,jj)=matrixtemp1(ii,jj+dimension_rank)
                            block%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj+dimension_rank)
                        enddo
                    enddo
                    !$omp end parallel do
                enddo
            enddo
        enddo
        deallocate (matrixtemp1)
    endif	  
    return

end subroutine create_butterfly_simple
