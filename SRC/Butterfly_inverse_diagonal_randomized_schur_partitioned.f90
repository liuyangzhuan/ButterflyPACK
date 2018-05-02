module Butterfly_inversion_schur_partition
use Butterfly_inversion
use Randomized_reconstruction
contains 





recursive subroutine Butterfly_inverse_partitionedinverse_IplusButter(level,ADflag)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock,ADflag
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,rank
    character chara
    real*8 T0
    type(matrixblock),pointer::blocks_io,blocks_A,blocks_B,blocks_C,blocks_D
    type(matrixblock)::blocks_schur
    integer rank_new_max,rank0
	real*8:: rank_new_avr,error,rate
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:)
	
	if(ADflag==1)then
		blocks_io=> partitioned_blocks(level)%blocks_D
	else 
		blocks_io=> partitioned_blocks(level)%blocks_A
	end if
	
	! write(*,*)'inverse ABDC',blocks_io%row_group,blocks_io%col_group,blocks_io%level,blocks_io%level_butterfly		
	
	if(blocks_io%level_butterfly==0)then
		call Butterfly_inverse_IplusButter_woodbury(blocks_io,Memory)
		return
    else 
		blocks_A => partitioned_blocks(level+1)%blocks_A
		blocks_B => partitioned_blocks(level+1)%blocks_B
		blocks_C => partitioned_blocks(level+1)%blocks_C
		blocks_D => partitioned_blocks(level+1)%blocks_D	
		! split into four smaller butterflies
		call Butterfly_split(blocks_io, blocks_A, blocks_B, blocks_C, blocks_D)
		
		! partitioned inverse of D
		call Butterfly_inverse_partitionedinverse_IplusButter(level+1,1)
		
		! construct the schur complement A-BD^-1C
		if(reducelevel_flag==1)then
			level_butterfly = blocks_A%level_butterfly
		else 
			! level_butterfly = int((Maxlevel-blocks_A%level)/2)*2
			level_butterfly = Maxlevel-blocks_A%level
		end if
		! write(*,*)'A-BDC',level_butterfly,level
		

		call get_minmaxrank_ABCD(partitioned_blocks(level+1),rank0)
		rate=1.2d0
		call Butterfly_randomized(level_butterfly,rank0,rate,blocks_A,partitioned_blocks(level+1),butterfly_block_MVP_inverse_A_minusBDinvC_dat,error_inout) 
					
		
		error_cnt = error_cnt + 1
		error_avr_glo = error_avr_glo + error_inout
		
		! partitioned inverse of the schur complement 
		call Butterfly_inverse_partitionedinverse_IplusButter(level+1,2)
		
		
		! construct the inverse
		if(reducelevel_flag==1)then
			level_butterfly = blocks_io%level_butterfly
		else 
			level_butterfly = Maxlevel-blocks_io%level
			if(level==0)level_butterfly = int((Maxlevel-blocks_io%level)/2)*2
		end if
		! write(*,*)'inverse ABDC',blocks_io%row_group,blocks_io%col_group,blocks_io%level,blocks_io%level_butterfly		
		
		! call Butterfly_inverse_ABCD(level_butterfly,blocks_io,partitioned_blocks(level+1),error_inout) 
		
		call get_minmaxrank_ABCD(partitioned_blocks(level+1),rank0)
		rate=1.2d0
		call Butterfly_randomized(level_butterfly,rank0,rate,blocks_io,partitioned_blocks(level+1),butterfly_block_MVP_inverse_ABCD_dat,error_inout) 
		
		
		error_cnt = error_cnt + 1
		error_avr_glo = error_avr_glo + error_inout
		! stop
		
		if(level==0 .and. verboselevel>=1)write(*,'(A23,A6,I3,A8,I3,A7,Es14.7)')' SchurI ',' rank:',blocks_io%rankmax,' L_butt:',blocks_io%level_butterfly,' error:',error_inout
			
		call delete_blocks(blocks_A)
		call delete_blocks(blocks_B)
		call delete_blocks(blocks_C)
		call delete_blocks(blocks_D)
		
		return
	
	end if	
end subroutine Butterfly_inverse_partitionedinverse_IplusButter



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
	
	! error_cnt = error_cnt+1
	
    return

end subroutine Butterfly_inverse_IplusButter_woodbury




subroutine Butterfly_split(blocks_i,blocks_A,blocks_B,blocks_C,blocks_D)
    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
	integer level_p,ADflag
	integer mm1,mm2,nn1,nn2,M1,M2,N1,N2,ii,jj,kk,j,i,mm,nn
	integer level_butterfly, num_blocks, level_butterfly_c, num_blocks_c,level,num_col,num_row,num_rowson,num_colson
	
    type(matrixblock)::blocks_i,blocks_A,blocks_B,blocks_C,blocks_D
	complex(kind=8),allocatable:: matrixtemp1(:,:),matrixtemp2(:,:),vin(:,:),vout1(:,:),vout2(:,:)
	complex(kind=8)::ctemp1,ctemp2

	blocks_A%level = blocks_i%level+1
	blocks_A%row_group = blocks_i%row_group*2	
	blocks_A%col_group = blocks_i%col_group*2
	blocks_A%style = blocks_i%style

	blocks_B%level = blocks_i%level+1
	blocks_B%row_group = blocks_i%row_group*2	
	blocks_B%col_group = blocks_i%col_group*2+1
	blocks_B%style = blocks_i%style
	
	blocks_C%level = blocks_i%level+1
	blocks_C%row_group = blocks_i%row_group*2+1	
	blocks_C%col_group = blocks_i%col_group*2
	blocks_C%style = blocks_i%style
	
	blocks_D%level = blocks_i%level+1
	blocks_D%row_group = blocks_i%row_group*2+1	
	blocks_D%col_group = blocks_i%col_group*2+1
	blocks_D%style = blocks_i%style	
	
	
	if(blocks_i%level_butterfly==0)then
		blocks_A%level_butterfly = 0
		blocks_B%level_butterfly = 0
		blocks_C%level_butterfly = 0
		blocks_D%level_butterfly = 0
		
		allocate(blocks_A%ButterflyU(1))
		allocate(blocks_B%ButterflyU(1))
		allocate(blocks_C%ButterflyU(1))
		allocate(blocks_D%ButterflyU(1))
		
		allocate(blocks_A%ButterflyV(1))
		allocate(blocks_B%ButterflyV(1))
		allocate(blocks_C%ButterflyV(1))
		allocate(blocks_D%ButterflyV(1))	

		mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
		nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
		mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
		nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1
		kk = size(blocks_i%ButterflyU(1)%matrix,2)
		
		
		
		allocate(blocks_A%ButterflyU(1)%matrix(mm1,kk))
		allocate(blocks_A%ButterflyV(1)%matrix(nn1,kk))		
		blocks_A%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1:mm1,1:kk)
		blocks_A%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1:nn1,1:kk)
		blocks_A%rankmax = kk
		blocks_A%rankmin = kk		
		
		allocate(blocks_B%ButterflyU(1)%matrix(mm1,kk))
		allocate(blocks_B%ButterflyV(1)%matrix(nn2,kk))		
		blocks_B%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1:mm1,1:kk)
		blocks_B%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1+nn1:nn1+nn2,1:kk)
		blocks_B%rankmax = kk
		blocks_B%rankmin = kk	
		
		allocate(blocks_C%ButterflyU(1)%matrix(mm2,kk))
		allocate(blocks_C%ButterflyV(1)%matrix(nn1,kk))		
		blocks_C%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1+mm1:mm1+mm2,1:kk)
		blocks_C%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1:nn1,1:kk)
		blocks_C%rankmax = kk
		blocks_C%rankmin = kk
		
		allocate(blocks_D%ButterflyU(1)%matrix(mm2,kk))
		allocate(blocks_D%ButterflyV(1)%matrix(nn2,kk))		
		blocks_D%ButterflyU(1)%matrix = blocks_i%ButterflyU(1)%matrix(1+mm1:mm1+mm2,1:kk)
		blocks_D%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix(1+nn1:nn1+nn2,1:kk)
		blocks_D%rankmax = kk
		blocks_D%rankmin = kk		
		
	
	else if(blocks_i%level_butterfly==1)then
		blocks_A%level_butterfly = 0
		blocks_B%level_butterfly = 0
		blocks_C%level_butterfly = 0
		blocks_D%level_butterfly = 0
		
		allocate(blocks_A%ButterflyU(1))
		allocate(blocks_B%ButterflyU(1))
		allocate(blocks_C%ButterflyU(1))
		allocate(blocks_D%ButterflyU(1))
		
		allocate(blocks_A%ButterflyV(1))
		allocate(blocks_B%ButterflyV(1))
		allocate(blocks_C%ButterflyV(1))
		allocate(blocks_D%ButterflyV(1))
		
		mm1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,1)
		nn1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,2)
		mm2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,1)
		nn2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,2)
		
		M1 = size(blocks_i%ButterflyU(1)%matrix,1)
		M2 = size(blocks_i%ButterflyU(2)%matrix,1)
		N1 = size(blocks_i%ButterflyV(1)%matrix,1)
		N2 = size(blocks_i%ButterflyV(2)%matrix,1)
		
		allocate(blocks_A%ButterflyU(1)%matrix(M1,nn1))
		allocate(blocks_A%ButterflyV(1)%matrix(N1,nn1))		
		call gemm_omp(blocks_i%ButterflyU(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,blocks_A%ButterflyU(1)%matrix, M1, mm1, nn1)
		blocks_A%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix
		blocks_A%rankmax = nn1
		blocks_A%rankmin = nn1
		
		allocate(blocks_B%ButterflyU(1)%matrix(M1,nn2))
		allocate(blocks_B%ButterflyV(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2)%matrix,blocks_B%ButterflyU(1)%matrix, M1, mm1, nn2)
		blocks_B%ButterflyV(1)%matrix = blocks_i%ButterflyV(2)%matrix		
		blocks_B%rankmax = nn2
		blocks_B%rankmin = nn2
		
		allocate(blocks_C%ButterflyU(1)%matrix(M2,nn1))
		allocate(blocks_C%ButterflyV(1)%matrix(N1,nn1))
		call gemm_omp(blocks_i%ButterflyU(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,1)%matrix,blocks_C%ButterflyU(1)%matrix, M2, mm2, nn1)
		blocks_C%ButterflyV(1)%matrix = blocks_i%ButterflyV(1)%matrix
		blocks_C%rankmax = nn1
		blocks_C%rankmin = nn1				
		
		allocate(blocks_D%ButterflyU(1)%matrix(M2,nn2))
		allocate(blocks_D%ButterflyV(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,blocks_D%ButterflyU(1)%matrix, M2, mm2, nn2)
		blocks_D%ButterflyV(1)%matrix = blocks_i%ButterflyV(2)%matrix			
		blocks_D%rankmax = nn2
		blocks_D%rankmin = nn2
		
		
	else 
		blocks_A%level_butterfly = blocks_i%level_butterfly-2
		blocks_B%level_butterfly = blocks_i%level_butterfly-2
		blocks_C%level_butterfly = blocks_i%level_butterfly-2
		blocks_D%level_butterfly = blocks_i%level_butterfly-2
		
		level_butterfly_c = blocks_i%level_butterfly-2
		num_blocks_c = 2**level_butterfly_c

		allocate(blocks_A%ButterflyU(num_blocks_c))
		allocate(blocks_A%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,2)
			allocate(blocks_A%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_A%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
		
			mm1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,1)
			allocate(blocks_A%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_A%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
		end do
		
		allocate(blocks_B%ButterflyU(num_blocks_c))
		allocate(blocks_B%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,2)
			allocate(blocks_B%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_B%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_B%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))	
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_B%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)			
		end do

		allocate(blocks_C%ButterflyU(num_blocks_c))
		allocate(blocks_C%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,2)
			allocate(blocks_C%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_C%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,1)
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			allocate(blocks_C%ButterflyV(ii)%matrix(mm1+mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_C%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)				
		end do
		allocate(blocks_D%ButterflyU(num_blocks_c))
		allocate(blocks_D%ButterflyV(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,2)
			allocate(blocks_D%ButterflyU(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_D%ButterflyU(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyU(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_D%ButterflyV(ii)%matrix(mm1+mm2,kk))
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_D%ButterflyV(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyV(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)		
		end do		
	
	
		if(level_butterfly_c/=0)then
			allocate(blocks_A%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_B%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_C%ButterflyKerl(level_butterfly_c))	
			allocate(blocks_D%ButterflyKerl(level_butterfly_c))	
		end if  
		do level=1, level_butterfly_c
             num_col=blocks_i%ButterflyKerl(level+1)%num_col
             num_row=blocks_i%ButterflyKerl(level+1)%num_row
             num_colson=num_col/2
             num_rowson=num_row/2
			 blocks_A%ButterflyKerl(level)%num_row=num_rowson
			 blocks_A%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_A%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_B%ButterflyKerl(level)%num_row=num_rowson
			 blocks_B%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_B%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_C%ButterflyKerl(level)%num_row=num_rowson
			 blocks_C%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_C%ButterflyKerl(level)%blocks(num_rowson,num_colson))
			 blocks_D%ButterflyKerl(level)%num_row=num_rowson
			 blocks_D%ButterflyKerl(level)%num_col=num_colson
			 allocate (blocks_D%ButterflyKerl(level)%blocks(num_rowson,num_colson))			 

			do j=1, num_col
				 do i=1, num_row
					 mm=size(blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix,1)
					 nn=size(blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix,2)
					 if (i<=num_rowson .and. j<=num_colson) then
						 allocate (blocks_A%ButterflyKerl(level)%blocks(i,j)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_A%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j<=num_colson) then
						 allocate (blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)                      
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i<=num_rowson .and. j>num_colson) then
						 allocate (blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j>num_colson) then
						 allocate (blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%matrix(mm,nn))
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 endif
				 enddo
			 enddo
		enddo
		
		call get_butterfly_minmaxrank(blocks_A)
		call get_butterfly_minmaxrank(blocks_B)
		call get_butterfly_minmaxrank(blocks_C)
		call get_butterfly_minmaxrank(blocks_D)
		
	end if
	
	! if(verboselevel>=1)write(*,'(A33,I5,A7,I3,I3,I3,I3)')'  In split. L_butt: ',blocks_i%level_butterfly,' rank: ', blocks_A%rankmax,blocks_B%rankmax,blocks_C%rankmax,blocks_D%rankmax
	

	
	mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
	nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
	mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
	nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1

	
	! allocate(vin(nn1+nn2,1))
	! vin = 1
	! allocate(vout1(mm1+mm2,1))
	! vout1 = 0
	! allocate(vout2(mm1+mm2,1))
	! vout2 = 0
	
	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_i,'N',mm1+mm2,nn1+nn2,1,vin,vout1,ctemp1,ctemp2)
	
	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_A,'N',mm1,nn1,1,vin(1:nn1,:),vout2(1:mm1,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call butterfly_block_MVP_randomized_dat(blocks_B,'N',mm1,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1:mm1,:),ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call butterfly_block_MVP_randomized_dat(blocks_C,'N',mm2,nn1,1,vin(1:nn1,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call butterfly_block_MVP_randomized_dat(blocks_D,'N',mm2,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	
	! write(*,*)'spliting error:',fnorm(vout1-vout2,mm1+mm2,1)/fnorm(vout1,mm1+mm2,1)
	! deallocate(vin,vout1,vout2)
	
	
end subroutine Butterfly_split


end module Butterfly_inversion_schur_partition
