module Butterfly_exact
use Utilites_randomized
use Butterfly_rightmultiply
! use Butterfly_compression_givenfullmat
use omp_lib
contains


subroutine Initialize_Butterfly_Exact(block)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max, dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j,blocks_idx
	type(matrixblock)::block
	integer seed_myid(12)
	integer time(8)
	
	call DATE_AND_TIME(values=time)     ! Get the current time 
	seed_myid(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
	call RANDOM_SEED(PUT=seed_myid)
	
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
	

	dimension_rank= block%rankmax  !rank_tmp    !!!!!!!!!!!!!!! be careful with the rank 


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
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyUInv(blocks)%matrix(dimension_rank,dimension_m))
        
		
		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,1)
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
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,1)
        do j=1, dimension_rank
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
                   
				    ! ! ! call RandomMat(2*dimension_rank,2*dimension_rank,2*dimension_rank,matrixtemp1,1)	
					
                    ! ! ! !$omp parallel do default(shared) private(ii,jj)
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
                    ! ! ! !$omp end parallel do
                    ! allocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(index_i,index_j)%matrix(2*dimension_rank,2*dimension_rank))
                ! enddo
            ! enddo
            ! call invert_Butterfly_Kernel(level)
        enddo
        deallocate (matrixtemp1)
    endif	
  
  
  
  
  
  
		! ! open(211,file='fort.211_i')
		! ! open(212,file='fort.212_i')
		! ! open(213,file='fort.213_i')
	
	
		! ! do blocks_idx=1, num_blocks
			! ! dimension_m=size(butterfly_block_randomized(1)%ButterflyU(blocks_idx)%matrix,1)
			! ! write(*,*)'U',dimension_m
		! ! end do
		
		! ! do blocks_idx=1, num_blocks
			! ! dimension_n=size(butterfly_block_randomized(1)%ButterflyV(blocks_idx)%matrix,1)
			! ! write(*,*)'V',dimension_n
		! ! end do		
	
		! ! do blocks_idx=1, num_blocks
			! ! dimension_m=size(butterfly_block_randomized(1)%ButterflyU(blocks_idx)%matrix,1)
			! ! do i=1, dimension_m
				! ! do j=1, dimension_rank
					! ! ! !read(211,*)tmpi,a,b
					! ! ! ! matrices_block(blocks,0)%ButterflyU(blocks_idx)%matrix(i,j)=a+junit*b
					
					! ! tmpi = 1
					! ! write(211,*)tmpi,real(butterfly_block_randomized(1)%ButterflyU(blocks_idx)%matrix(i,j)),aimag(butterfly_block_randomized(1)%ButterflyU(blocks_idx)%matrix(i,j))
					
				! ! enddo
			! ! enddo	
			! ! dimension_n=size(butterfly_block_randomized(1)%ButterflyV(blocks_idx)%matrix,1)
			! ! do j=1, dimension_rank
				! ! do i=1, dimension_n
					! ! ! read(213,*)tmpi,a,b
					! ! ! matrices_block(blocks,0)%ButterflyV(blocks_idx)%matrix(i,j)=a+junit*b
					! ! tmpi = 1
					! ! write(213,*)tmpi,real(butterfly_block_randomized(1)%ButterflyV(blocks_idx)%matrix(i,j)),aimag(butterfly_block_randomized(1)%ButterflyV(blocks_idx)%matrix(i,j))
					
				! ! enddo
			! ! enddo
			! ! ! write(*,*)matrices_block(blocks,0)%ButterflyV(1)%matrix
		! ! enddo	

		! ! ! write(*,*)matrices_block(blocks,0)%ButterflyV(1)%matrix
		! ! ! stop
		
	   ! ! if (level_butterfly/=0) then
			! ! do level=1, level_butterfly
				! ! num_row=2**level
				! ! num_col=2**(level_butterfly-level+1)
				! ! do i=1, num_row
					! ! do j=1, num_col
						! ! do jj=1, dimension_rank
							! ! do ii=1, dimension_rank                           
								! ! ! read(212,*)tmpi,tmpj,a,b
								! ! ! matrices_block(blocks,0)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=a+junit*b
								! ! tmpi = 1
								! ! tmpj = 1
								! ! write(212,*)tmpi,tmpj,real(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)),aimag(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj))
								
							! ! enddo
						! ! enddo
					! ! enddo
				! ! enddo
			! ! enddo
		! ! endif
		! ! close(211)
		! ! close(212)
		! ! close(213)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    return

end subroutine Initialize_Butterfly_Exact



subroutine Initialize_Butterfly_Outter_Exact(block,kover)
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
	! if(level_c==1)dimension_rank=max(dimension_rank,maxlevel_for_blocks)+kover

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
        allocate (butterfly_block_randomized(1)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyU_old(blocks)%matrix(dimension_m,dimension_rank))
        ! allocate (butterfly_block_randomized(1)%ButterflyUInv(blocks)%matrix(dimension_rank,dimension_m))
        
		
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

end subroutine Initialize_Butterfly_Outter_Exact




subroutine Get_Randomized_Vectors_LL_Exact(block_o,nth_s,nth_e,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    type(matrixblock)::block_o
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	
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
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(2))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	

    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(1)%vector(mm,num_vect_sub))
    allocate (RandomVectors_InOutput(2)%vector(nn,num_vect_sub))
	do ii =1,2
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupm_start=groupm*2**level_butterfly
	header_mm=basis_group(groupm_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		do i=1, num_blocks
			if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_m=basis_group(groupm_start+i-1)%head
				tailer_m=basis_group(groupm_start+i-1)%tail
				mm=tailer_m-header_m+1
				k=header_m-header_mm	

				allocate(matrixtemp1(dimension_rank,mm))
				call RandomMat(dimension_rank,mm,min(mm,dimension_rank),matrixtemp1,1)

				!$omp parallel do default(shared) private(ii,jj)
				 do jj=1,dimension_rank
					 do ii=1, mm
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*dimension_rank+jj)=matrixtemp1(jj,ii) !random_complex_number()	
					 enddo
				 enddo
				 !$omp end parallel do
				 
				 deallocate(matrixtemp1)
			 end if
		end do
	end do
	
	! ! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	! ! do jj=1, num_vect_sub
		 ! ! do ii=1, mm
			 ! ! write(115,*)dble(RandomVectors_InOutput(1)%vector(ii,jj)),dimag(RandomVectors_InOutput(1)%vector(ii,jj))
		! ! end do
	! ! end do		
	
	

	! get the left multiplied vectors	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
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
				random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(2)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+nn
	enddo 

	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	nn=basis_group(groupn)%tail-basis_group(groupn)%head+1
	! ! write(*,*)fnorm(RandomVectors_InOutput(1)%vector,mm,num_vect_sub),fnorm(RandomVectors_InOutput(2)%vector,nn,num_vect_sub),'gogo'
	! ! write(*,*)fnorm(RandomVectors_InOutput(1)%vector(:,1:1),mm,1),fnorm(RandomVectors_InOutput(2)%vector(:,1:1),nn,1),'gogo'
	
	
	

	
	
	
	
	
	
	
	
	
	
	! ! ! do ii=1,mm
		! ! ! write(555,*)abs(RandomVectors_InOutput(1)%vector(ii,1))
	! ! ! end do
	! ! ! write(555,*)'gog'
	
	! ! ! do ii=1,nn
		! ! ! write(556,*)abs(RandomVectors_InOutput(2)%vector(ii,1))
	! ! ! end do
	! ! ! write(556,*)'gog'	
	

	! ! k=0
	! ! do i=1, num_blocks
		! ! nn=size(butterfly_block_randomized(1)%butterflyV(i)%matrix,1)
		! ! write(*,*)fnorm(random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(:,1:1),nn,1),'heihei'
		! ! k=k+nn
	! ! enddo 

	! ! k=0
	! ! do i=1, num_blocks
		! ! mm=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		! ! write(*,*)fnorm(random%RandomVectorLL(0)%blocks(i,1)%matrix(:,1:1),mm,1),'ddheihei'
		! ! k=k+mm
	! ! enddo 
	
    !$omp parallel do default(shared) private(i)
    do i=1, 2
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_LL_Exact









subroutine Get_Randomized_Vectors_RR_Exact(block_o,nth_s,nth_e,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock)::block_o
	
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
	
	integer nth_s,nth_e,num_vect_sub,nth
	type(RandomBlock), pointer :: random
	
	 
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(2))

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
	do ii =1,2
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupn_start=groupn*2**level_butterfly
	header_nn=basis_group(groupn_start)%head
	idx_start = 1
	
	do nth= nth_s,nth_e
		do i=1, num_blocks
			if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
				header_n=basis_group(groupn_start+i-1)%head
				tailer_n=basis_group(groupn_start+i-1)%tail
				nn=tailer_n-header_n+1
				k=header_n-header_nn

				allocate(matrixtemp1(dimension_rank,nn))
				call RandomMat(dimension_rank,nn,min(nn,dimension_rank),matrixtemp1,1)
				
				!$omp parallel do default(shared) private(ii,jj)
				 do jj=1,dimension_rank
					 do ii=1, nn
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*dimension_rank+jj)=matrixtemp1(jj,ii)  !random_complex_number()	
					 enddo
				 enddo
				 !$omp end parallel do
				 
				 deallocate(matrixtemp1)
				 
			 end if
		end do
	end do


	! ! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	! ! do jj=1, num_vect_sub
		 ! ! do ii=1, nn
			 ! ! write(116,*)dble(RandomVectors_InOutput(1)%vector(ii,jj)),dimag(RandomVectors_InOutput(1)%vector(ii,jj))
		! ! end do
	! ! end do	
	
	! get the right multiplied vectors
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
    ! call butterfly_block_MVP_randomized_printCond(block_o,'N',random1,random2,ctemp1,ctemp2)
	


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
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(2)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+mm
	enddo 
	
    !$omp parallel do default(shared) private(i)
    do i=1, 2
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_Exact




subroutine Get_Randomized_Vectors_RR_Test_Exact(block_o,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
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
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(2))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
    mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	allocate (RandomVectors_InOutput(2)%vector(mm,num_vect_sub))

	do ii =1,2
		RandomVectors_InOutput(ii)%vector = 0
	end do	 
	 
	groupn_start=groupn*2**level_butterfly
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
    random2=>RandomVectors_InOutput(2)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
    ! call butterfly_block_MVP_randomized_printCond(block_o,'N',random1,random2,ctemp1,ctemp2)


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
				random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(2)%vector(ii+k,jj)
			enddo
		enddo
		!$omp end parallel do
		k=k+mm
	enddo 
	
    !$omp parallel do default(shared) private(i)
    do i=1, 2
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine Get_Randomized_Vectors_RR_Test_Exact




! subroutine Checkrank_Sblock_exact(level_c,rowblock)

    ! use MODULE_FILE
    ! ! use lapack95
	! ! use blas95
    ! implicit none
    
	! integer level_c,rowblock
    ! integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,m,n
    ! integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    ! character chara
    ! real*8 a,b,c,d
    ! complex(kind=8) ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o
	! type(matrixblock):: block_tmp	
    ! type(vectorsblock), pointer :: random1, random2
    
    ! real*8,allocatable :: Singular(:)
	! integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	! complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	! real*8::rtemp

	! integer header_m,header_n,edge_m,edge_n,group_m,group_n
	! complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat1(:,:) 
	! integer, allocatable :: ipiv(:)
	
    
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	
	! block_tmp%col_group = block_o%col_group
	! block_tmp%row_group = block_o%row_group
	! block_tmp%level = block_o%level

	! group_m = block_tmp%row_group
	! group_n = block_tmp%col_group
	! mm=basis_group(block_tmp%row_group)%tail-basis_group(block_tmp%row_group)%head+1
	! nn=basis_group(block_tmp%col_group)%tail-basis_group(block_tmp%col_group)%head+1
	! header_m=basis_group(group_m)%head
	! header_n=basis_group(group_n)%head
	! allocate(block_tmp%fullmat(mm,nn))
	! ! write(*,*)shape(block_tmp%fullmat)
	
	
	! do m = 1,mm
		! do n = 1,nn
			! edge_m = header_m + m-1
			! edge_n = header_n + n-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! block_tmp%fullmat(m,n) = ctemp
		! end do
	! end do	

	! allocate(fullmat(mm,nn))
	! fullmat = block_tmp%fullmat
	
	! allocate(matrixtmp1(mm,mm))
	! do m=1,mm
		! do n =1,mm
			! edge_m = header_m + m-1
			! edge_n = header_m + n-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! matrixtmp1(m,n) = ctemp	
		! end do
	! end do


	! allocate(ipiv(mm))
	! call getrff90(matrixtmp1,ipiv)
	! ! write(*,*)shape(matrixtemp1)
	! call getrif90(matrixtmp1,ipiv)	
	! deallocate(ipiv)

	
	! allocate(matrixtmp2(nn,nn))
	! do m=1,nn
		! do n =1,nn
			! edge_m = header_n + m-1
			! edge_n = header_n + n-1
			! call element_Zmn(edge_m,edge_n,ctemp)
			! matrixtmp2(m,n) = ctemp			
		! end do
	! end do
	
	! allocate(ipiv(nn))
	! call getrff90(matrixtmp2,ipiv)
	! ! write(*,*)shape(matrixtemp1)
	! call getrif90(matrixtmp2,ipiv)	
	! deallocate(ipiv)				

	! allocate(fullmat1(mm,nn))
	
	! ! call gemmf90(matrixtmp1, fullmat, fullmat1,'N','N') 
	! ! call gemmf90(fullmat1, matrixtmp2, block_tmp%fullmat,'N','N') 

	! call gemmf90(matrixtmp1, fullmat, block_tmp%fullmat,'N','N') 

	! ! call gemmf90(fullmat, matrixtmp2, block_tmp%fullmat,'N','N') 
			
	! ! call Butterfly_compress_givenfullmat(block_tmp,rtemp)
	! ! call Butterfly_compress_QRSVD_givenfullmat(block_tmp,rtemp)		
	! call Butterfly_compress_ID_givenfullmat(block_tmp,rtemp)		
	
	
	! if(block_o%level==2)then
	! write(*,*)mm,nn
	! do ii=1,mm
	! do jj=1,nn
		! write(112,*)abs(block_tmp%fullmat(ii,jj))
	! end do
	! end do	
	! stop
	! end if	
	
    ! return                

! end subroutine Checkrank_Sblock_exact


subroutine Reconstruction_LL_Exact(block_o)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub
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
    type(matrixblock)::block_o
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter
	real*8:: error_inout
	integer,allocatable:: perms(:)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	
	Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)   !  check here later
	! ! Nsub = 2**level_butterfly  
	
    Ng = 2**level_butterfly/Nsub
	
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
	
	! write(*,*)Nsub,Ng,dimension_rank
	
    allocate (Random_Block(1))
    
    allocate (Random_Block(1)%RandomVectorLL(0:level_butterfly+2))    
	Nbind = 1
	if(Nsub<Nbind)Nbind = Nsub
	
	num_vect_sub = dimension_rank*Nbind
	
    random=>Random_Block(1)
	call Init_RandVect('T',random,num_vect_sub)
	
	level_right_start = floor_safe(level_butterfly/2d0) !  check here later
	! ! level_right_start = level_butterfly+1
	
	call Zero_Butterfly(0,level_right_start)
 

	do ii = 1,Nsub/Nbind
		
		nth_s = (ii-1)*Nbind+1
		nth_e = ii*Nbind

		n1 = OMP_get_wtime()
		call Get_Randomized_Vectors_LL_Exact(block_o,nth_s,nth_e,num_vect_sub)
		n2 = OMP_get_wtime()
		time_getvec = time_getvec + n2-n1
		
		n1 = OMP_get_wtime()
		call Resolving_Butterfly_LL(num_vect_sub,nth_s,nth_e,Ng)
		n2 = OMP_get_wtime()
		time_resolve = time_resolve + n2-n1		
	end do


	
	
	random=>Random_Block(1)
	do level=0, level_butterfly+2
		num_row=random%RandomVectorLL(level)%num_row
		num_col=random%RandomVectorLL(level)%num_col
		do j=1, num_col
			do i=1, num_row
				deallocate (random%RandomVectorLL(level)%blocks(i,j)%matrix)
				if(level/=0 .and. level/=level_butterfly+2)deallocate (random%RandomVectorLL(level)%blocks(i,j)%null)
			enddo
		enddo
		deallocate (random%RandomVectorLL(level)%blocks)
	enddo
	deallocate (random%RandomVectorLL)

    return
    
end subroutine Reconstruction_LL_Exact





subroutine Reconstruction_LL_Outter_Exact(Bplus)
    
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
    type(blockplus)::Bplus
	
	
    ! type(matrix), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)

	level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
	
	
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
			call Get_Randomized_Vectors_LL_Outter_Exact(Bplus,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			! Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()
			call Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			! time_resolve = time_resolve + n2-n1		
		end do
	end do
	
	! pause
	! call Resolving_Butterfly_LL_rankcompletion()
	
	! deallocate(perms)
	
	
	random=>Random_Block(1)
	call Delete_RandVect('T',random,level_butterfly)

    return
    
end subroutine Reconstruction_LL_Outter_Exact





subroutine Reconstruction_RR_Exact(block_o,error)
    
    use MODULE_FILE
    implicit none
	    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub
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
	integer niter
	integer,allocatable::perms(:)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)   !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)
	end if
	
	! ! Nsub = 2**level_butterfly  
	
    Ng = 2**level_butterfly/Nsub
    num_blocks=2**level_butterfly
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 
		
	! ! call assert(num_vectors==2**level_butterfly*dimension_rank,'incorrect num_vectors')
	
	! call assert(num_vectors==Nsub*dimension_rank,'incorrect num_vectors') !  check here later
	
    ! ! allocate (Random_Block(1))   !  check here later 
    
    allocate (Random_Block(1)%RandomVectorRR(0:level_butterfly+2))    
	
	Nbind = 1
	if(Nsub<Nbind)Nbind = Nsub
	
	num_vect_sub = dimension_rank*Nbind
	
    random=>Random_Block(1)
	call Init_RandVect('N',random,num_vect_sub)

    level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
    ! level_left_start = 0
	
	call Zero_Butterfly(level_left_start,level_butterfly+1)
		
	do ii = 1,Nsub/Nbind
		nth_s = (ii-1)*Nbind+1
		nth_e = ii*Nbind
		n1 = OMP_get_wtime()
		call Get_Randomized_Vectors_RR_Exact(block_o,nth_s,nth_e,num_vect_sub)
		n2 = OMP_get_wtime()
		time_getvec = time_getvec + n2-n1	
		
		n1 = OMP_get_wtime()		
		call Resolving_Butterfly_RR(num_vect_sub,nth_s,nth_e,Ng)
		n2 = OMP_get_wtime()
		time_resolve = time_resolve + n2-n1			
	end do

	random=>Random_Block(1)
	do level=0, level_butterfly+2
		num_row=random%RandomVectorRR(level)%num_row
		num_col=random%RandomVectorRR(level)%num_col
		do j=1, num_col
			do i=1, num_row
				deallocate (random%RandomVectorRR(level)%blocks(i,j)%matrix)
				if(level/=0 .and. level/=level_butterfly+2)deallocate (random%RandomVectorRR(level)%blocks(i,j)%null)
			enddo
		enddo
		deallocate (random%RandomVectorRR(level)%blocks)
	enddo
	deallocate (random%RandomVectorRR)
	
	deallocate(Random_Block)

	call Test_Error_RR_Exact(block_o,error)


	
    return
    
end subroutine Reconstruction_RR_Exact




subroutine Reconstruction_RR_Outter_Exact(Bplus,error)
    
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
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
	
    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(blockplus)::Bplus
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	integer,allocatable::perms(:)
		
	level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
		
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
			call Get_Randomized_Vectors_RR_Outter_Exact(Bplus,nth_s,nth_e,num_vect_sub,unique_nth)
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1	
			! Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()		
			call Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth)
			n2 = OMP_get_wtime()
			! time_resolve = time_resolve + n2-n1			
		end do
	end do
	
	! deallocate(perms)

	random=>Random_Block(1)
	call Delete_RandVect('N',random,level_butterfly)
	
	deallocate(Random_Block)

	call Test_Error_RR_Outter_Exact(Bplus,error)


	
    return
    
end subroutine Reconstruction_RR_Outter_Exact





subroutine Test_Error_RR_Outter_Exact(Bplus,error)

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
	type(blockplus)::Bplus
	

	
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

	call Get_Randomized_Vectors_RR_Test_Outter_Exact(Bplus,num_vect)

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

end subroutine Test_Error_RR_Outter_Exact






subroutine Test_Error_RR_Exact(block_o,error)

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
	type(matrixblock)::block_o
	
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

	call Get_Randomized_Vectors_RR_Test_Exact(block_o,num_vect)

	k=0
	do i=1, num_blocks
		dimension_m=size(butterfly_block_randomized(1)%butterflyU(i)%matrix,1)
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, dimension_m
			do jj=1, num_vect
				RandomVectors_Output_ref(ii+k,jj)=random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)
			enddo
		enddo
		!$omp end parallel do
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
				 norm1_R=norm1_R+abs(RandomVectors_Output_ref(ii+k,jj))**2
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

end subroutine Test_Error_RR_Exact



subroutine Butterfly_randomized_Exact_test_memfree(block_o)

    use MODULE_FILE
	! use lapack95
    ! use blas95	
	use misc
    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,nn1,nn2
    character chara
    real*8 T0
    type(matrixblock)::block_o
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter
	real*8:: error_inout
	real*8:: maxvalue(1:20)
	complex(kind=8),allocatable:: tmp_vec(:),tmp_mat1(:,:),tmp_mat2(:,:)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	if(level_butterfly==0)then
		write(*,*)'not implemented '
		stop
	else 
		call Initialize_Butterfly_Exact(block_o)


		call Reconstruction_LL_Exact(block_o)	
		call Reconstruction_RR_Exact(block_o,error_inout)		
		

		
		level_butterfly=butterfly_block_randomized(1)%level_butterfly
		num_blocks=2**level_butterfly
		
		maxvalue = 0
		rank_new_max = 0
		rank_new_avr = 0
		do i=1, num_blocks
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)
			
	
			mm=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(1) = max(maxvalue(1),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyV(i)%matrix(:,jj),mm)))
			end do
			
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyU(i)%matrix(:,jj),mm)))	
			end do			
			
			
			rank_new_max = max(rank_new_max,nn)
			rank_new_avr = rank_new_avr + nn
			
			deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
		enddo
		rank_new_avr = rank_new_avr/num_blocks
		deallocate (butterfly_block_randomized(1)%ButterflyU, butterfly_block_randomized(1)%ButterflyV)
		if (level_butterfly/=0) then
			do level=1, level_butterfly
				num_row=butterfly_block_randomized(1)%ButterflyKerl(level)%num_row
				num_col=butterfly_block_randomized(1)%ButterflyKerl(level)%num_col
				do i=1, num_row
					do j=1, num_col
						mm=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						nn=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						
						! do ii = 1,mm
						! do jj = 1,nn
							! maxvalue(level+1) = max(maxvalue(level+1),abs(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)))	
							! ! write(*,*)butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)
						! end do
						! end do
						
						
						
						if(mod(j,2)==1)then
							allocate(tmp_mat1(mm,nn))
							tmp_mat1 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix
						else
							allocate(tmp_mat2(mm,nn))
							tmp_mat2 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix

							nn2=size(tmp_mat2,2)
							nn1=size(tmp_mat1,2)						
								
							allocate(tmp_vec(nn1+nn2))
							do ii = 1, mm
								tmp_vec(1:nn1) = tmp_mat1(ii,:)
								tmp_vec(1+nn1:nn2+nn1) = tmp_mat2(ii,:)
								maxvalue(level+1) = max(maxvalue(level+1),sqrt(norm_vector(tmp_vec,nn2+nn1)))
							end do
							deallocate(tmp_vec)
							
							deallocate(tmp_mat1)
							deallocate(tmp_mat2)
						end if
						
					
						deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix)
					enddo
				enddo
				do j=1, num_col/2
					do i=1, num_row/2
						deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(i,j)%matrix)
					enddo
				enddo
				deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks)
				deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks)
			enddo
			deallocate (butterfly_block_randomized(1)%ButterflyKerl)
			deallocate (butterfly_block_randomized(1)%ButterflyInv)
		endif
		
		!$omp parallel do default(shared) private(i)
		do i=1, num_blocks
			!deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			!deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyUInv(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyVInv(i)%matrix)
		enddo
		!$omp end parallel do
		!deallocate (butterfly_block_randomized(1)%ButterflyU)
		!deallocate (butterfly_block_randomized(1)%ButterflyV)
		deallocate (butterfly_block_randomized(1)%ButterflyUInv)
		deallocate (butterfly_block_randomized(1)%ButterflyVInv)
		
		deallocate (butterfly_block_randomized)
	end if
    ! ! Time_Rightmultiply_inverse_randomized=Time_Rightmultiply_inverse_randomized+secnds(T0) 
    
	
	! ! call Test_Randomized_Sblock(level_c,rowblock,block_old,error)
	! ! call delete_blocks(block_old)
	! ! write(*,*)level_c,rowblock,rank_new_max,error,level_butterfly
	if(level_butterfly>0)then
		write(*,*)'randomized exact memfree',level_butterfly,error_inout
	end if
    return

end subroutine Butterfly_randomized_Exact_test_memfree



subroutine Butterfly_randomized_Exact_test(block_o)

    use MODULE_FILE
	! use lapack95
    ! use blas95	
	use misc
    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,nn1,nn2
    character chara
    real*8 T0
    type(matrixblock)::block_o
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter
	real*8:: error_inout
	real*8:: maxvalue(1:20)
	complex(kind=8),allocatable:: tmp_vec(:),tmp_mat1(:,:),tmp_mat2(:,:)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	if(level_butterfly==0)then
		write(*,*)'not implemented '
		stop
	else 
		call Initialize_Butterfly_Exact(block_o)
		call Get_Randomized_Vectors_Exact(block_o)        
		call Resolving_Butterfly_v8(niter,error_inout,block_o)

		
		level_butterfly=butterfly_block_randomized(1)%level_butterfly
		num_blocks=2**level_butterfly
		
		maxvalue = 0
		rank_new_max = 0
		rank_new_avr = 0
		do i=1, num_blocks
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)
			
			
			! do ii = 1,mm
			! do jj = 1,nn
				! maxvalue(1) = max(maxvalue(1),abs(butterfly_block_randomized(1)%ButterflyV(i)%matrix(ii,jj)))	
				! maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),abs(butterfly_block_randomized(1)%ButterflyU(i)%matrix(ii,jj)))	
			! end do
			! end do
	
			mm=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(1) = max(maxvalue(1),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyV(i)%matrix(:,jj),mm)))
			end do
			
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyU(i)%matrix(:,jj),mm)))	
			end do			
			
			
			rank_new_max = max(rank_new_max,nn)
			rank_new_avr = rank_new_avr + nn
			
			deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
		enddo
		rank_new_avr = rank_new_avr/num_blocks
		deallocate (butterfly_block_randomized(1)%ButterflyU, butterfly_block_randomized(1)%ButterflyV)
		if (level_butterfly/=0) then
			do level=1, level_butterfly
				num_row=butterfly_block_randomized(1)%ButterflyKerl(level)%num_row
				num_col=butterfly_block_randomized(1)%ButterflyKerl(level)%num_col
				do i=1, num_row
					do j=1, num_col
						mm=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						nn=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						
						! do ii = 1,mm
						! do jj = 1,nn
							! maxvalue(level+1) = max(maxvalue(level+1),abs(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)))	
							! ! write(*,*)butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)
						! end do
						! end do
						
						
						
						if(mod(j,2)==1)then
							allocate(tmp_mat1(mm,nn))
							tmp_mat1 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix
						else
							allocate(tmp_mat2(mm,nn))
							tmp_mat2 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix

							nn2=size(tmp_mat2,2)
							nn1=size(tmp_mat1,2)						
								
							allocate(tmp_vec(nn1+nn2))
							do ii = 1, mm
								tmp_vec(1:nn1) = tmp_mat1(ii,:)
								tmp_vec(1+nn1:nn2+nn1) = tmp_mat2(ii,:)
								maxvalue(level+1) = max(maxvalue(level+1),sqrt(norm_vector(tmp_vec,nn2+nn1)))
							end do
							deallocate(tmp_vec)
							
							deallocate(tmp_mat1)
							deallocate(tmp_mat2)
						end if
						
					
						deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix)
					enddo
				enddo
				do j=1, num_col/2
					do i=1, num_row/2
						deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(i,j)%matrix)
					enddo
				enddo
				deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks)
				deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks)
			enddo
			deallocate (butterfly_block_randomized(1)%ButterflyKerl)
			deallocate (butterfly_block_randomized(1)%ButterflyInv)
		endif
		
		write(*,*)'max value: ',maxvalue(1:9)
		
		
		! ! ! !$omp parallel do default(shared) private(i)
		! ! ! do i=1, 6
			! ! ! deallocate (RandomVectors_InOutput(i)%vector)
		! ! ! enddo
		! ! ! !$omp end parallel do
		! ! ! deallocate (RandomVectors_InOutput)

		! ! ! do k=1, 2
			! ! ! do level=0, level_butterfly+2
				! ! ! num_row=Random_Block(k)%RandomVectorRR(level)%num_row
				! ! ! num_col=Random_Block(k)%RandomVectorRR(level)%num_col
				! ! ! do j=1, num_col
					! ! ! do i=1, num_row
						! ! ! deallocate (Random_Block(k)%RandomVectorRR(level)%blocks(i,j)%matrix)
					! ! ! enddo
				! ! ! enddo
				! ! ! deallocate (Random_Block(k)%RandomVectorRR(level)%blocks)
			! ! ! enddo
			! ! ! deallocate (Random_Block(k)%RandomVectorRR)
			! ! ! do level=0, level_butterfly+2
				! ! ! num_row=Random_Block(k)%RandomVectorLL(level)%num_row
				! ! ! num_col=Random_Block(k)%RandomVectorLL(level)%num_col
				! ! ! do j=1, num_col
					! ! ! do i=1, num_row
						! ! ! deallocate (Random_Block(k)%RandomVectorLL(level)%blocks(i,j)%matrix)
					! ! ! enddo
				! ! ! enddo
				! ! ! deallocate (Random_Block(k)%RandomVectorLL(level)%blocks)
			! ! ! enddo
			! ! ! deallocate (Random_Block(k)%RandomVectorLL)
		! ! ! enddo
		! ! ! deallocate (Random_Block)
		
		!$omp parallel do default(shared) private(i)
		do i=1, num_blocks
			!deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			!deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyUInv(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyVInv(i)%matrix)
		enddo
		!$omp end parallel do
		!deallocate (butterfly_block_randomized(1)%ButterflyU)
		!deallocate (butterfly_block_randomized(1)%ButterflyV)
		deallocate (butterfly_block_randomized(1)%ButterflyUInv)
		deallocate (butterfly_block_randomized(1)%ButterflyVInv)
		
		deallocate (butterfly_block_randomized)
	end if
    ! ! Time_Rightmultiply_inverse_randomized=Time_Rightmultiply_inverse_randomized+secnds(T0) 
    
	
	! ! call Test_Randomized_Sblock(level_c,rowblock,block_old,error)
	! ! call delete_blocks(block_old)
	! ! write(*,*)level_c,rowblock,rank_new_max,error,level_butterfly
	if(level_butterfly>0)then
		write(*,*)'randomized exact',rank_new_max,level_butterfly,niter,error_inout
	end if
    return

end subroutine Butterfly_randomized_Exact_test













subroutine Butterfly_randomized_Exact_test_symmetric(block_o)

    use MODULE_FILE
	! use lapack95
    ! use blas95	
	use misc
    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,nn1,nn2
    character chara
    real*8 T0
    type(matrixblock)::block_o
	type(matrixblock)::block_old
    integer::rank_new_max
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter
	real*8:: error_inout
	real*8:: maxvalue(1:20)
	complex(kind=8),allocatable:: tmp_vec(:),tmp_mat1(:,:),tmp_mat2(:,:)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	if(level_butterfly==0)then
		write(*,*)'not implemented '
		stop
	else 
		call Get_Randomized_Vectors_Exact_Memeff(block_o)
		call Resolving_Butterfly_Symmetric(block_o%row_group,block_o%col_group)
			
		call Test_Error_RR_Exact(block_o,error_inout)
		
		
		level_butterfly=butterfly_block_randomized(1)%level_butterfly
		num_blocks=2**level_butterfly
		
		maxvalue = 0
		rank_new_max = 0
		rank_new_avr = 0
		do i=1, num_blocks
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)
			
			
			! do ii = 1,mm
			! do jj = 1,nn
				! maxvalue(1) = max(maxvalue(1),abs(butterfly_block_randomized(1)%ButterflyV(i)%matrix(ii,jj)))	
				! maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),abs(butterfly_block_randomized(1)%ButterflyU(i)%matrix(ii,jj)))	
			! end do
			! end do
	
			mm=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyV(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(1) = max(maxvalue(1),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyV(i)%matrix(:,jj),mm)))
			end do
			
			mm=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,1)
			nn=size(butterfly_block_randomized(1)%ButterflyU(i)%matrix,2)	
			do jj = 1,nn
				maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),sqrt(norm_vector(butterfly_block_randomized(1)%ButterflyU(i)%matrix(:,jj),mm)))	
			end do			
			
			
			rank_new_max = max(rank_new_max,nn)
			rank_new_avr = rank_new_avr + nn
			
			deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
		enddo
		rank_new_avr = rank_new_avr/num_blocks
		deallocate (butterfly_block_randomized(1)%ButterflyU, butterfly_block_randomized(1)%ButterflyV)
		if (level_butterfly/=0) then
			do level=1, level_butterfly
				num_row=butterfly_block_randomized(1)%ButterflyKerl(level)%num_row
				num_col=butterfly_block_randomized(1)%ButterflyKerl(level)%num_col
				do i=1, num_row
					do j=1, num_col
						mm=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						nn=size(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						
						! do ii = 1,mm
						! do jj = 1,nn
							! maxvalue(level+1) = max(maxvalue(level+1),abs(butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)))	
							! ! write(*,*)butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)
						! end do
						! end do
						
						
						
						if(mod(j,2)==1)then
							allocate(tmp_mat1(mm,nn))
							tmp_mat1 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix
						else
							allocate(tmp_mat2(mm,nn))
							tmp_mat2 = butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix

							nn2=size(tmp_mat2,2)
							nn1=size(tmp_mat1,2)						
								
							allocate(tmp_vec(nn1+nn2))
							do ii = 1, mm
								tmp_vec(1:nn1) = tmp_mat1(ii,:)
								tmp_vec(1+nn1:nn2+nn1) = tmp_mat2(ii,:)
								maxvalue(level+1) = max(maxvalue(level+1),sqrt(norm_vector(tmp_vec,nn2+nn1)))
							end do
							deallocate(tmp_vec)
							
							deallocate(tmp_mat1)
							deallocate(tmp_mat2)
						end if
						
					
						deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks(i,j)%matrix)
					enddo
				enddo
				! do j=1, num_col/2
					! do i=1, num_row/2
						! deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks(i,j)%matrix)
					! enddo
				! enddo
				! deallocate (butterfly_block_randomized(1)%ButterflyInv(level)%blocks)
				deallocate (butterfly_block_randomized(1)%ButterflyKerl(level)%blocks)
			enddo
			deallocate (butterfly_block_randomized(1)%ButterflyKerl)
			! deallocate (butterfly_block_randomized(1)%ButterflyInv)
		endif
		
		write(*,*)'max value: ',maxvalue(1:9)

		!$omp parallel do default(shared) private(i)
		do i=1, num_blocks
			!deallocate (butterfly_block_randomized(1)%ButterflyU(i)%matrix)
			!deallocate (butterfly_block_randomized(1)%ButterflyV(i)%matrix)
			! deallocate (butterfly_block_randomized(1)%ButterflyUInv(i)%matrix)
			! deallocate (butterfly_block_randomized(1)%ButterflyVInv(i)%matrix)
		enddo
		!$omp end parallel do
		!deallocate (butterfly_block_randomized(1)%ButterflyU)
		!deallocate (butterfly_block_randomized(1)%ButterflyV)
		! deallocate (butterfly_block_randomized(1)%ButterflyUInv)
		! deallocate (butterfly_block_randomized(1)%ButterflyVInv)
		
		deallocate (butterfly_block_randomized)
	end if
    ! ! Time_Rightmultiply_inverse_randomized=Time_Rightmultiply_inverse_randomized+secnds(T0) 
    
	
	! ! call Test_Randomized_Sblock(level_c,rowblock,block_old,error)
	! ! call delete_blocks(block_old)
	! ! write(*,*)level_c,rowblock,rank_new_max,error,level_butterfly
	if(level_butterfly>0)then
		write(*,*)'randomized exact',rank_new_max,level_butterfly,error_inout
	end if
    return

end subroutine Butterfly_randomized_Exact_test_symmetric



subroutine Get_Randomized_Vectors_Exact(block_o)

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
	type(matrixblock)::block_o
	
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
	
	
    level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(6))

    groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)
    Ng = 2**level_butterfly/Nsub
	dimension_rank = butterfly_block_randomized(1)%dimension_rank
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
	 
 
	level_blocks = block_o%level



	groupn_start=groupn*2**level_butterfly
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


	
	groupm_start=groupm*2**level_butterfly
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

	! ! ! open(115,file='fort.115')
		! ! ! do jj=1, num_vectors
			! ! ! do ii=1, nn
				! ! ! read(115,*)a,b
				! ! ! RandomVectors_InOutput(1)%vector(ii,jj)=a+junit*b
			! ! ! end do
		! ! ! end do
	! ! ! close(115)
	! ! ! open(116,file='fort.116')
		! ! ! do jj=1, num_vectors
			! ! ! do ii=1, mm
				! ! ! read(116,*)a,b
				! ! ! RandomVectors_InOutput(4)%vector(ii,jj)=a+junit*b
			! ! ! end do
		! ! ! end do
	! ! ! close(116)	
	
	
	! get the right multiplied vectors
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'N',random1,random2,ctemp1,ctemp2)
    ! call butterfly_block_MVP_randomized_printCond(block_o,'N',random1,random2,ctemp1,ctemp2)
	

	
	
	! get the left multiplied vectors	
    random1=>RandomVectors_InOutput(4)
    random2=>RandomVectors_InOutput(6)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	          
    ! call butterfly_block_MVP_randomized_printCond(block_o,'T',random1,random2,ctemp1,ctemp2)	          
    return                

end subroutine Get_Randomized_Vectors_Exact






subroutine Get_Randomized_Vectors_Exact_Memeff(block_o)

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
	type(matrixblock)::block_o
	
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
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2
	
	dimension_rank= block_o%rankmax
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

    call Butterfly_Partial_MVP_dat_memeff(block_o,'T',0,level_butterfly+1,DimMax,num_vectorsL,vec_old,vec_new)      

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

end subroutine Get_Randomized_Vectors_Exact_Memeff





subroutine Get_Randomized_Vectors_LL_Outter_Exact(Bplus,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer unique_nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2, ctemp3, ctemp4
	type(blockplus)::Bplus
	
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
    ctemp3=-1.0d0 ; ctemp4=1.0d0	
	
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(3)%vector(nn,num_vect_sub))
    
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
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

				allocate(matrixtemp1(num_vect_subsub,mm))
				call RandomMat(num_vect_subsub,mm,min(mm,num_vect_subsub),matrixtemp1,0)
				
				! !$omp parallel do default(shared) private(ii,jj)
				 do jj=1,num_vect_subsub
					 do ii=1, mm
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
					 enddo
				 enddo
				 ! !$omp end parallel do
				 deallocate(matrixtemp1)
				 
			 end if
		end do
	end do
	
	! get the left multiplied vectors
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	idx_start_glo = basis_group(groupm)%head		

	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
	
	n1 = OMP_get_wtime()
    ! call butterfly_block_MVP_randomized(block_o,'T',random1,random2,ctemp1,ctemp2)	
	
	call Bplus_block_MVP_randomized_dat(Bplus,'T',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)	
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'T',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	
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

end subroutine Get_Randomized_Vectors_LL_Outter_Exact






subroutine Get_Randomized_Vectors_RR_Outter_Exact(Bplus,nth_s,nth_e,num_vect_sub,unique_nth)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2, ctemp3, ctemp4
	type(blockplus)::Bplus
	
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
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

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
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
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

				allocate(matrixtemp1(num_vect_subsub,nn))
				call RandomMat(num_vect_subsub,nn,min(nn,num_vect_subsub),matrixtemp1,0)
				
				! !$omp parallel do default(shared) private(ii,jj)
				 do jj=1,num_vect_subsub
					 do ii=1, nn
						 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number() ! matrixtemp1(jj,ii) ! 	
					 enddo
				 enddo
				 ! !$omp end parallel do
				 
				 deallocate(matrixtemp1)
			 end if
		end do
	end do
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 	
	
	
	! get the right multiplied vectors
	idx_start_glo = basis_group(groupm)%head
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    ctemp3=-1.0d0 ; ctemp4=1.0d0
n1 = OMP_get_wtime()  
  
  call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)
  call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
  n2 = OMP_get_wtime()
! time_tmp = time_tmp + n2 - n1	
   


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

end subroutine Get_Randomized_Vectors_RR_Outter_Exact




subroutine Get_Randomized_Vectors_RR_Test_Outter_Exact(Bplus,num_vect_sub)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	! integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2,  ctemp3, ctemp4
	type(blockplus)::Bplus
	
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
	
	
	! block_o =>  cascading_factors(level_c)%matrices_block(rowblock) 
	  
    level_butterfly=int((maxlevel_for_blocks-Bplus%level)/2)*2
    num_blocks=2**level_butterfly
    allocate (RandomVectors_InOutput(3))

    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 

	dimension_rank =butterfly_block_randomized(1)%dimension_rank 

	
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	allocate (RandomVectors_InOutput(1)%vector(nn,num_vect_sub))
    
    groupm=Bplus%row_group  ! Note: row_group and col_group interchanged here   
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
	
    groupn=Bplus%col_group  ! Note: row_group and col_group interchanged here   
    nn=basis_group(groupn)%tail-basis_group(groupn)%head+1  	
	
	! get the right multiplied vectors
	
    random1=>RandomVectors_InOutput(1)
    random2=>RandomVectors_InOutput(3)
    ctemp1=1.0d0 ; ctemp2=0.0d0
    ctemp3=-1.0d0 ; ctemp4=1.0d0

	call Bplus_block_MVP_randomized_dat(Bplus,'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp1,ctemp2)		
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',mm,nn,num_vect_sub,random1%vector,random2%vector,ctemp3,ctemp4,2,Bplus_randomized(1)%Lplus)	
	


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

end subroutine Get_Randomized_Vectors_RR_Test_Outter_Exact



subroutine Bplus_randomized_Exact_test(bplus)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    type(blockplus)::bplus
	integer:: ii,ll,bb
    real*8 Memory,rtemp,error	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall,M,N
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	complex(kind=8) ctemp, ctemp1, ctemp2
	
	call assert(bplus%Lplus>=2,'this is not a Bplus in Bplus_randomized_Exact_test')
	
	call Initialize_Bplus_FromInput(bplus)
	
	
	do bb =1,Bplus_randomized(1)%LL(2)%Nbound
		! write(*,*)bb,Bplus_randomized(1)%LL(2)%Nbound,'dddd'
		call Bplus_randomized_Exact_onesubblock(bplus,bb)
		! write(*,*)'go'
	end do
	call Test_Error_RR_Inner_Exact(Bplus)
		
	call Bplus_randomized_outter_Exact_memfree(bplus,Memory)
	
    return

end subroutine Bplus_randomized_Exact_test







subroutine Test_Error_RR_Inner_Exact(Bplus)

    use MODULE_FILE
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm,ll
    integer mm,nn
    real*8 a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    complex(kind=8) ctemp,ctemp1,ctemp2
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly,rank_new_max
	integer*8 idx_start
	real*8::error
	integer dimension_m, M,N
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)	
	type(blockplus)::Bplus
	

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! test the results of the inner factorization
	
	M = basis_group(Bplus_randomized(1)%row_group)%tail - basis_group(Bplus_randomized(1)%row_group)%head + 1
	N = basis_group(Bplus_randomized(1)%col_group)%tail - basis_group(Bplus_randomized(1)%col_group)%head + 1
	ctemp1=1.0d0 ; ctemp2=0.0d0
	
	allocate(RandVectInR(N,1))
	allocate(RandVectOutR(M,1))
	RandVectInR=0
	do ii=1,N
		RandVectInR(ii,1) = random_complex_number()
	end do

	call Bplus_block_MVP_randomized_dat_partial(Bplus,'N',M,N,1,RandVectInR,RandVectOutR,ctemp1,ctemp2,2,Bplus%Lplus)
	allocate(Vout1(M,1))
	Vout1 = RandVectOutR

	
	call Bplus_block_MVP_randomized_dat_partial(Bplus_randomized(1),'N',M,N,1,RandVectInR,RandVectOutR,ctemp1,ctemp2,2,Bplus_randomized(1)%Lplus)
	allocate(Vout2(M,1))
	Vout2 = RandVectOutR
	deallocate(RandVectInR)
	deallocate(RandVectOutR)	
	
	
	error = fnorm(Vout2-Vout1,M,1)/fnorm(Vout1,M,1)
	! write(*,*)error,'ninini after factorization'	
	deallocate(Vout1)
	deallocate(Vout2)	

	
	if(error>iter_tolerance)then
		write(*,*)'inner factorization not correct, needs more work'
		stop
	end if		

	rank_new_max = 0
	do ll=2,Bplus_randomized(1)%Lplus
		rank_new_max = max(rank_new_max,Bplus_randomized(1)%LL(ll)%rankmax)
	end do
	
	write(*,'(A15,I3,A9,I3,A7,Es14.7)')'Inner: rank:',rank_new_max,' Lplus-1:',Bplus%Lplus,' error:',error	
	
	
    return                

end subroutine Test_Error_RR_Inner_Exact





 

subroutine Bplus_randomized_Exact_onesubblock(bplus,bb_o)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   use Butterfly_compress_forward
   implicit none

    type(blockplus)::bplus
	integer:: ii,ll,bb,jj,bb_o,tt
    real*8 Memory,rtemp,error	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	integer M,N,idx_start_n,idx_start_m,idx_start_n_loc,idx_end_n_loc,idx_start_m_loc,idx_end_m_loc,mm,nn,rmax,rank,idx_start_n_ref,idx_start_m_ref,idx_end_n_ref,idx_end_m_ref
	complex(kind=8)::ctemp1,ctemp2
	type(matrixblock),pointer::blocks
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:)
	real*8, allocatable :: Singular(:)
	
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
	
		rmax = bplus%LL(2)%rankmax*2 + (tt-1)*10  !!!!! be careful here
	
		!!!!!!!!!!!!!!!!!!!!!!!! get right multiplied results
		
		allocate(RandVectInR(N,rmax))
		RandVectInR=0
		allocate(RandVectOutR(M,rmax))
		do ii=idx_start_n_loc,idx_end_n_loc
		do jj=1,rmax
			RandVectInR(ii,jj)=random_complex_number()
		end do
		end do

		
		call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,rmax,RandVectInR,RandVectOutR,ctemp1,ctemp2)
		
		allocate(matRcol(nn,rmax))
		matRcol = RandVectInR(idx_start_n_loc:idx_start_n_loc+nn-1,1:rmax)	
		deallocate(RandVectInR)
		allocate(matZRcol(mm,rmax))
		matZRcol = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:rmax)
		deallocate(RandVectOutR)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		
				
		!!!!!!!!!!!!!!!!!!!!!!!!!!!! get left multiplied results
		
		allocate(RandVectInL(M,rmax))
		RandVectInL=0
		allocate(RandVectOutL(N,rmax))
		do ii=idx_start_m_loc,idx_end_m_loc
		do jj=1,rmax
			RandVectInL(ii,jj)=random_complex_number()
		end do
		end do
				
		call Bplus_block_MVP_randomized_dat(bplus,'T',M,N,rmax,RandVectInL,RandVectOutL,ctemp1,ctemp2)
		
				
		allocate(matRrow(mm,rmax))
		matRrow = RandVectInL(idx_start_m_loc:idx_start_m_loc+mm-1,1:rmax)
		matRrow = conjg(matRrow)
		deallocate(RandVectInL)
		allocate(matZcRrow(nn,rmax))
		matZcRrow = RandVectOutL(idx_start_n_loc:idx_start_n_loc+nn-1,1:rmax)
		matZcRrow = conjg(matZcRrow)	
		deallocate(RandVectOutL)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		allocate(matU_glo(mm,rmax))
		allocate(matV_glo(rmax,nn))
		allocate(Singular(rmax))

		call RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU_glo,matV_glo,Singular,mm,nn,rmax,rank,LS_tolerance,SVD_tolerance_factor)				
		
		do ii=1,rank
			matV_glo(ii,:) = matV_glo(ii,:) * Singular(ii)
		end do
		
		deallocate(matRcol,matZRcol,matRrow,matZcRrow,Singular)



		if(mm*nn*16/1e9>5)then
			write(*,*)'warning: full storage of matSub_glo is too costly'
		end if
		allocate(matSub_glo(mm,nn))
		call gemm_omp(matU_glo(1:mm,1:rank),matV_glo(1:rank,1:nn),matSub_glo,mm,rank,nn)
		deallocate(matU_glo,matV_glo)

		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! test the results
		allocate(Vin(nn,1))
		allocate(Vout1(mm,1))
		allocate(Vout2(mm,1))
		do ii=1,nn
			Vin(ii,1)=random_complex_number()
		end do
		
		
		allocate(RandVectInR(N,1))
		allocate(RandVectOutR(M,1))
		RandVectInR=0
		RandVectInR(idx_start_n_loc:idx_end_n_loc,1:1) = Vin
		call Bplus_block_MVP_randomized_dat(bplus,'N',M,N,1,RandVectInR,RandVectOutR,ctemp1,ctemp2)
		Vout1 = RandVectOutR(idx_start_m_loc:idx_start_m_loc+mm-1,1:1)
		deallocate(RandVectInR)
		deallocate(RandVectOutR)

		call gemm_omp(matSub_glo,Vin,Vout2,mm,nn,1)
		
		error = fnorm(Vout2-Vout1,mm,1)/fnorm(Vout1,mm,1)
		! write(*,*)error,bb_o,'ninini'		

		deallocate(Vin)
		deallocate(Vout1)
		deallocate(Vout2)
		
		if(error>iter_tolerance)then
			deallocate(matSub_glo)
			if(min(mm,nn)==rmax)then
				write(*,*)tt,rmax,'Exact_onesubblock',error,rank
				write(*,*)'no need to increase rmax, try increase RandomizedSVD tolerance'
			end if
		else
			! write(*,*)error,bb_o,'good in Bplus_randomized_Exact_onesubblock'
			exit
		end if
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	end do

	if(error>iter_tolerance)then
		write(*,*)'matSub_glo no correct, needs more work'
		stop
	end if

	!!!! construct butterfly at all subsequent levels including level 2 (this part is generic and doesn't require accuracy test)
	idx_start_n_ref = idx_start_n_loc + idx_start_n - 1
	idx_end_n_ref = idx_end_n_loc + idx_start_n - 1
	idx_start_m_ref = idx_start_m_loc + idx_start_m - 1
	idx_end_m_ref = idx_end_m_loc + idx_start_m - 1


	do ll=2,Bplus_randomized(1)%Lplus
	
		level_butterfly = int((Maxlevel_for_blocks - Bplus_randomized(1)%LL(ll)%matrices_block(1)%level)/2)*2 
		level_BP = Bplus_randomized(1)%level
		levelm = ceiling_safe(dble(level_butterfly)/2d0)						
		groupm_start=Bplus_randomized(1)%LL(ll)%matrices_block(1)%row_group*2**levelm		
		Nboundall = 2**(Bplus_randomized(1)%LL(ll)%matrices_block(1)%level+levelm-level_BP)
		
		
		do bb = 1,Bplus_randomized(1)%LL(ll)%Nbound
			if(basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%head>=idx_start_m_ref .and. basis_group(Bplus_randomized(1)%LL(ll)%matrices_block(bb)%row_group)%tail<=idx_end_m_ref)then
			
				if(Bplus_randomized(1)%LL(ll+1)%Nbound==0)then
					write(*,*)'666',ll
					call Butterfly_compress_N15_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),idx_start_m_ref,idx_start_n_ref)
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
				else
					write(*,*)'777'				
					call Butterfly_compress_N15_withoutBoundary_givenfullmat(Bplus_randomized(1)%LL(ll)%matrices_block(bb),Bplus_randomized(1)%LL(ll+1)%boundary_map,Nboundall,groupm_start, rtemp, idx_start_m_ref,idx_start_n_ref)
					call Butterfly_sym2asym(Bplus_randomized(1)%LL(ll)%matrices_block(bb))
				end if				
			
			end if	
			Bplus_randomized(1)%LL(ll)%rankmax = max(Bplus_randomized(1)%LL(ll)%rankmax,Bplus_randomized(1)%LL(ll)%matrices_block(bb)%rankmax)			
		end do
	end do

	deallocate(matSub_glo)
		
    return

end subroutine Bplus_randomized_Exact_onesubblock






subroutine Bplus_randomized_outter_Exact_memfree(Bplus,Memory)

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
	type(blockplus)::Bplus
	
	Memory = 0
	
    block_o =>  Bplus%LL(1)%matrices_block(1)
	
	level_butterfly=int((maxlevel_for_blocks-block_o%level)/2)*2

	do tt =1,10
		do ntry=1,1
		n1 = OMP_get_wtime()
		call Initialize_Butterfly_Outter_Exact(block_o,tt-1)
		n2 = OMP_get_wtime()
		! Time_Init_forward = Time_Init_forward + n2 -n1 
		

		n1 = OMP_get_wtime()
		call Reconstruction_LL_Outter_Exact(Bplus)	
		call Reconstruction_RR_Outter_Exact(Bplus,error_inout)
		n2 = OMP_get_wtime()
		! Time_Reconstruct_forward = Time_Reconstruct_forward + n2-n1

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
			write(*,'(A15,I3,A8,I2,A8,I3,A7,Es14.7)')'Outter: rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout
			return
		end if
		end do
	end do
	write(*,*)'randomized scheme not converged in Bplus_randomized_outter_Exact_memfree',error_inout
	stop
	
    return

end subroutine Bplus_randomized_outter_Exact_memfree








end module Butterfly_exact