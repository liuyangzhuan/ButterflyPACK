module Butterfly_inversion_schur_partition
use Butterfly_inversion
use Randomized_reconstruction
contains 

recursive subroutine Butterfly_inverse_partitionedinverse_IplusButter(blocks_io,level_butterfly_target,option,error_inout,stats,ptree,pgno)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,kk1,kk2,rank,err_cnt
    character chara
    real*8 T0,err_avr
    type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    type(matrixblock)::blocks_io
    type(matrixblock)::blocks_schur
    integer rank_new_max,rank0
	real*8:: rank_new_avr,error,rate
	integer niter
	real*8:: error_inout
	integer itermax,ntry
	real*8:: n1,n2,Memory
	complex (kind=8), allocatable::matrix_small(:,:)
	type(Hoption)::option
	type(Hstat)::stats
	integer level_butterfly_target,pgno,pgno1
	type(proctree)::ptree
	
	type(partitionedblocks)::partitioned_block
	
	error_inout=0
	
	! write(*,*)'inverse ABDC',blocks_io%row_group,blocks_io%col_group,blocks_io%level,blocks_io%level_butterfly		
	
	if(blocks_io%level_butterfly==0)then
		call Butterfly_inverse_IplusButter_woodbury(blocks_io,Memory,ptree,stats,pgno)
		return
    else
		allocate(partitioned_block%blocks_A)
		allocate(partitioned_block%blocks_B)
		allocate(partitioned_block%blocks_C)
		allocate(partitioned_block%blocks_D)
	
		blocks_A => partitioned_block%blocks_A
		blocks_B => partitioned_block%blocks_B
		blocks_C => partitioned_block%blocks_C
		blocks_D => partitioned_block%blocks_D	
		
		! split into four smaller butterflies
		call Butterfly_split(blocks_io, blocks_A, blocks_B, blocks_C, blocks_D)
		
		! partitioned inverse of D
		! level_butterfly=level_butterfly_target-1
		level_butterfly=blocks_D%level_butterfly
		if(GetTreelevel(pgno)==ptree%nlevel)then  
			pgno1=pgno
		else
			pgno1=pgno*2
		endif
		call Butterfly_inverse_partitionedinverse_IplusButter(blocks_D,level_butterfly,option,error,stats,ptree,pgno1)
		error_inout = max(error_inout, error)
		
		! construct the schur complement A-BD^-1C

		! level_butterfly = level_butterfly_target-1
		level_butterfly = blocks_A%level_butterfly

		! write(*,*)'A-BDC',level_butterfly,level
		

		call get_minmaxrank_ABCD(partitioned_block,rank0)
		rate=1.2d0
		call Butterfly_randomized(level_butterfly,rank0,rate,blocks_A,partitioned_block,butterfly_block_MVP_inverse_A_minusBDinvC_dat,error,'A-BD^-1C',option,stats,ptree) 
		error_inout = max(error_inout, error)
		
		! write(*,*)'ddd1'
		! partitioned inverse of the schur complement 
		! level_butterfly=level_butterfly_target-1
		level_butterfly=blocks_A%level_butterfly
		if(GetTreelevel(pgno)==ptree%nlevel)then
			pgno1=pgno
		else
			pgno1=pgno*2+1
		endif		
		call Butterfly_inverse_partitionedinverse_IplusButter(blocks_A,level_butterfly,option,error,stats,ptree,pgno1)
		error_inout = max(error_inout, error)		

		level_butterfly = level_butterfly_target
		call get_minmaxrank_ABCD(partitioned_block,rank0)
		rate=1.2d0
		call Butterfly_randomized(level_butterfly,rank0,rate,blocks_io,partitioned_block,butterfly_block_MVP_inverse_ABCD_dat,error,'ABCDinverse',option,stats,ptree) 
		error_inout = max(error_inout, error)		
		
		
		
		! stop
#if PRNTlevel >= 2		
		if(level==0)write(*,'(A23,A6,I3,A8,I3,A11,Es14.7)')' SchurI ',' rank:',blocks_io%rankmax,' L_butt:',blocks_io%level_butterfly,' error:',error_inout
#endif			
		call delete_blocks(blocks_A)
		call delete_blocks(blocks_B)
		call delete_blocks(blocks_C)
		call delete_blocks(blocks_D)
		
		return
	
	end if	
end subroutine Butterfly_inverse_partitionedinverse_IplusButter


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
	
		write(*,*)'the following is very unlikely to be used'
		stop
		
		! blocks_A%level_butterfly = 0
		! blocks_B%level_butterfly = 0
		! blocks_C%level_butterfly = 0
		! blocks_D%level_butterfly = 0
		
		! allocate(blocks_A%ButterflyU%blocks(1))
		! allocate(blocks_B%ButterflyU%blocks(1))
		! allocate(blocks_C%ButterflyU%blocks(1))
		! allocate(blocks_D%ButterflyU%blocks(1))
		
		! allocate(blocks_A%ButterflyV%blocks(1))
		! allocate(blocks_B%ButterflyV%blocks(1))
		! allocate(blocks_C%ButterflyV%blocks(1))
		! allocate(blocks_D%ButterflyV%blocks(1))	

		! mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
		! nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
		! mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
		! nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1
		! kk = size(blocks_i%ButterflyU%blocks(1)%matrix,2)
		
		
		
		! allocate(blocks_A%ButterflyU%blocks(1)%matrix(mm1,kk))
		! allocate(blocks_A%ButterflyV%blocks(1)%matrix(nn1,kk))		
		! blocks_A%ButterflyU%blocks(1)%matrix = blocks_i%ButterflyU%blocks(1)%matrix(1:mm1,1:kk)
		! blocks_A%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix(1:nn1,1:kk)
		! blocks_A%rankmax = kk
		! blocks_A%rankmin = kk		
		
		! allocate(blocks_B%ButterflyU%blocks(1)%matrix(mm1,kk))
		! allocate(blocks_B%ButterflyV%blocks(1)%matrix(nn2,kk))		
		! blocks_B%ButterflyU%blocks(1)%matrix = blocks_i%ButterflyU%blocks(1)%matrix(1:mm1,1:kk)
		! blocks_B%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix(1+nn1:nn1+nn2,1:kk)
		! blocks_B%rankmax = kk
		! blocks_B%rankmin = kk	
		
		! allocate(blocks_C%ButterflyU%blocks(1)%matrix(mm2,kk))
		! allocate(blocks_C%ButterflyV%blocks(1)%matrix(nn1,kk))		
		! blocks_C%ButterflyU%blocks(1)%matrix = blocks_i%ButterflyU%blocks(1)%matrix(1+mm1:mm1+mm2,1:kk)
		! blocks_C%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix(1:nn1,1:kk)
		! blocks_C%rankmax = kk
		! blocks_C%rankmin = kk
		
		! allocate(blocks_D%ButterflyU%blocks(1)%matrix(mm2,kk))
		! allocate(blocks_D%ButterflyV%blocks(1)%matrix(nn2,kk))		
		! blocks_D%ButterflyU%blocks(1)%matrix = blocks_i%ButterflyU%blocks(1)%matrix(1+mm1:mm1+mm2,1:kk)
		! blocks_D%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix(1+nn1:nn1+nn2,1:kk)
		! blocks_D%rankmax = kk
		! blocks_D%rankmin = kk		
		
	
	else if(blocks_i%level_butterfly==1)then
		blocks_A%level_butterfly = 0
		blocks_B%level_butterfly = 0
		blocks_C%level_butterfly = 0
		blocks_D%level_butterfly = 0
		
		allocate(blocks_A%ButterflyU%blocks(1))
		allocate(blocks_B%ButterflyU%blocks(1))
		allocate(blocks_C%ButterflyU%blocks(1))
		allocate(blocks_D%ButterflyU%blocks(1))
		
		allocate(blocks_A%ButterflyV%blocks(1))
		allocate(blocks_B%ButterflyV%blocks(1))
		allocate(blocks_C%ButterflyV%blocks(1))
		allocate(blocks_D%ButterflyV%blocks(1))
		
		mm1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,1)
		nn1 = size(blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,2)
		mm2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,1)
		nn2 = size(blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,2)
		
		M1 = size(blocks_i%ButterflyU%blocks(1)%matrix,1)
		M2 = size(blocks_i%ButterflyU%blocks(2)%matrix,1)
		N1 = size(blocks_i%ButterflyV%blocks(1)%matrix,1)
		N2 = size(blocks_i%ButterflyV%blocks(2)%matrix,1)
		
		
		
		blocks_A%ButterflyU%blocks(1)%mdim=M1
		blocks_A%ButterflyU%blocks(1)%ndim=nn1
		blocks_A%ButterflyV%blocks(1)%mdim=N1
		blocks_A%ButterflyV%blocks(1)%ndim=nn1		
		allocate(blocks_A%ButterflyU%blocks(1)%matrix(M1,nn1))
		allocate(blocks_A%ButterflyV%blocks(1)%matrix(N1,nn1))		
		call gemm_omp(blocks_i%ButterflyU%blocks(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,1)%matrix,blocks_A%ButterflyU%blocks(1)%matrix, M1, mm1, nn1)
		blocks_A%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix
		blocks_A%rankmax = nn1
		blocks_A%rankmin = nn1
		
			
		blocks_B%ButterflyU%blocks(1)%mdim=M1
		blocks_B%ButterflyU%blocks(1)%ndim=nn2
		blocks_B%ButterflyV%blocks(1)%mdim=N2
		blocks_B%ButterflyV%blocks(1)%ndim=nn2		
		allocate(blocks_B%ButterflyU%blocks(1)%matrix(M1,nn2))
		allocate(blocks_B%ButterflyV%blocks(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU%blocks(1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2)%matrix,blocks_B%ButterflyU%blocks(1)%matrix, M1, mm1, nn2)
		blocks_B%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(2)%matrix		
		blocks_B%rankmax = nn2
		blocks_B%rankmin = nn2
		
		
		blocks_C%ButterflyU%blocks(1)%mdim=M2
		blocks_C%ButterflyU%blocks(1)%ndim=nn1
		blocks_C%ButterflyV%blocks(1)%mdim=N1
		blocks_C%ButterflyV%blocks(1)%ndim=nn1			
		allocate(blocks_C%ButterflyU%blocks(1)%matrix(M2,nn1))
		allocate(blocks_C%ButterflyV%blocks(1)%matrix(N1,nn1))
		call gemm_omp(blocks_i%ButterflyU%blocks(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,1)%matrix,blocks_C%ButterflyU%blocks(1)%matrix, M2, mm2, nn1)
		blocks_C%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(1)%matrix
		blocks_C%rankmax = nn1
		blocks_C%rankmin = nn1				
		
				
		blocks_D%ButterflyU%blocks(1)%mdim=M2
		blocks_D%ButterflyU%blocks(1)%ndim=nn2
		blocks_D%ButterflyV%blocks(1)%mdim=N2
		blocks_D%ButterflyV%blocks(1)%ndim=nn2			
		allocate(blocks_D%ButterflyU%blocks(1)%matrix(M2,nn2))
		allocate(blocks_D%ButterflyV%blocks(1)%matrix(N2,nn2))
		call gemm_omp(blocks_i%ButterflyU%blocks(2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2)%matrix,blocks_D%ButterflyU%blocks(1)%matrix, M2, mm2, nn2)
		blocks_D%ButterflyV%blocks(1)%matrix = blocks_i%ButterflyV%blocks(2)%matrix			
		blocks_D%rankmax = nn2
		blocks_D%rankmin = nn2
		
		
		blocks_A%headm=blocks_i%headm
		blocks_A%M=M1
		blocks_A%headn=blocks_i%headn
		blocks_A%N=N1	
		blocks_B%headm=blocks_i%headm
		blocks_B%M=M1	
		blocks_B%headn=blocks_i%headn+N1
		blocks_B%N=N2
		blocks_C%headm=blocks_i%headm+M1
		blocks_C%M=M2			
		blocks_C%headn=blocks_i%headn
		blocks_C%N=N1			
		blocks_D%headm=blocks_i%headm+M1
		blocks_D%M=M2			
		blocks_D%headn=blocks_i%headn+N1
		blocks_D%N=N2
		
		write(*,*)'M_loc,N_loc,N_p,M_p needs to be defined for ABCD'
		stop		
		
		
	else 
		blocks_A%level_butterfly = blocks_i%level_butterfly-2
		blocks_B%level_butterfly = blocks_i%level_butterfly-2
		blocks_C%level_butterfly = blocks_i%level_butterfly-2
		blocks_D%level_butterfly = blocks_i%level_butterfly-2
		
		level_butterfly_c = blocks_i%level_butterfly-2
		num_blocks_c = 2**level_butterfly_c
		
		allocate(blocks_A%ButterflyU%blocks(num_blocks_c))
		allocate(blocks_A%ButterflyV%blocks(num_blocks_c))
		M1=0
		N1=0
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU%blocks(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU%blocks(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,2)
			allocate(blocks_A%ButterflyU%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_A%ButterflyU%blocks(ii)%mdim=mm1+mm2
			blocks_A%ButterflyU%blocks(ii)%ndim=kk
			M1=M1+blocks_A%ButterflyU%blocks(ii)%mdim
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_A%ButterflyU%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyU%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
		
			mm1 = size(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV%blocks(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV%blocks(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,1)
			allocate(blocks_A%ButterflyV%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_A%ButterflyV%blocks(ii)%mdim=mm1+mm2
			N1=N1+blocks_A%ButterflyV%blocks(ii)%mdim
			blocks_A%ButterflyV%blocks(ii)%ndim=kk
			
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_A%ButterflyV%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_A%ButterflyV%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
		end do
		
		allocate(blocks_B%ButterflyU%blocks(num_blocks_c))
		allocate(blocks_B%ButterflyV%blocks(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyU%blocks(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyU%blocks(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,2)
			allocate(blocks_B%ButterflyU%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_B%ButterflyU%blocks(ii)%mdim=mm1+mm2
			blocks_B%ButterflyU%blocks(ii)%ndim=kk
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii-1)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_B%ButterflyU%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyU%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_B%ButterflyV%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_B%ButterflyV%blocks(ii)%mdim=mm1+mm2
			blocks_B%ButterflyV%blocks(ii)%ndim=kk
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))	
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(1,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_B%ButterflyV%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_B%ButterflyV%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)			
		end do

		allocate(blocks_C%ButterflyU%blocks(num_blocks_c))
		allocate(blocks_C%ButterflyV%blocks(num_blocks_c))
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,2)
			allocate(blocks_C%ButterflyU%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_C%ButterflyU%blocks(ii)%mdim=mm1+mm2
			blocks_C%ButterflyU%blocks(ii)%ndim=kk
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,1)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,1)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_C%ButterflyU%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyU%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,1)
			nn1 = size(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,2)
			mm2 = size(blocks_i%ButterflyV%blocks(2*ii)%matrix,1)
			nn2 = size(blocks_i%ButterflyV%blocks(2*ii)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,1)
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))			
			allocate(blocks_C%ButterflyV%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_C%ButterflyV%blocks(ii)%mdim=mm1+mm2
			blocks_C%ButterflyV%blocks(ii)%ndim=kk			
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii-1)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_C%ButterflyV%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_C%ButterflyV%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)				
		end do
		
		allocate(blocks_D%ButterflyU%blocks(num_blocks_c))
		allocate(blocks_D%ButterflyV%blocks(num_blocks_c))
		M2=0
		N2=0		
		do ii =1,num_blocks_c
			mm1 = size(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,2)
			allocate(blocks_D%ButterflyU%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_D%ButterflyU%blocks(ii)%mdim=mm1+mm2
			M2=M2+blocks_D%ButterflyU%blocks(ii)%mdim
			blocks_D%ButterflyU%blocks(ii)%ndim=kk			
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii-1+num_blocks_c*2,2)%matrix,matrixtemp1,mm1,nn1,kk)
			call gemm_omp(blocks_i%ButterflyU%blocks(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(level_butterfly_c+2)%blocks(2*ii+num_blocks_c*2,2)%matrix,matrixtemp2,mm2,nn2,kk)
			blocks_D%ButterflyU%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyU%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2
			deallocate(matrixtemp1)
			deallocate(matrixtemp2)
			
			mm1 = size(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,1)
			nn1 = size(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,2)
			mm2 = size(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,1)
			nn2 = size(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,2)
			kk = size(blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,1)
			allocate(blocks_D%ButterflyV%blocks(ii)%matrix(mm1+mm2,kk))
			blocks_D%ButterflyV%blocks(ii)%mdim=mm1+mm2
			N2=N2+blocks_D%ButterflyV%blocks(ii)%mdim
			blocks_D%ButterflyV%blocks(ii)%ndim=kk					
			allocate(matrixtemp1(mm1,kk))
			allocate(matrixtemp2(mm2,kk))
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii-1+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii-1+num_blocks_c*2)%matrix,matrixtemp1, mm1,nn1,kk)
			call gemmNT_omp(blocks_i%ButterflyV%blocks(2*ii+num_blocks_c*2)%matrix,blocks_i%ButterflyKerl(1)%blocks(2,2*ii+num_blocks_c*2)%matrix,matrixtemp2, mm2,nn2,kk)
			blocks_D%ButterflyV%blocks(ii)%matrix(1:mm1,1:kk) = matrixtemp1
			blocks_D%ButterflyV%blocks(ii)%matrix(1+mm1:mm1+mm2,1:kk) = matrixtemp2			
			deallocate(matrixtemp1) 
			deallocate(matrixtemp2)		
		end do		

		blocks_A%headm=blocks_i%headm
		blocks_A%M=M1
		blocks_A%headn=blocks_i%headn
		blocks_A%N=N1	
		blocks_B%headm=blocks_i%headm
		blocks_B%M=M1	
		blocks_B%headn=blocks_i%headn+N1
		blocks_B%N=N2
		blocks_C%headm=blocks_i%headm+M1
		blocks_C%M=M2			
		blocks_C%headn=blocks_i%headn
		blocks_C%N=N1			
		blocks_D%headm=blocks_i%headm+M1
		blocks_D%M=M2			
		blocks_D%headn=blocks_i%headn+N1
		blocks_D%N=N2
		
		write(*,*)'M_loc,N_loc,N_p,M_p needs to be defined for ABCD'
		stop
	
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
						 blocks_A%ButterflyKerl(level)%blocks(i,j)%mdim=mm
						 blocks_A%ButterflyKerl(level)%blocks(i,j)%ndim=nn
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_A%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j<=num_colson) then
						 allocate (blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(mm,nn))
						 blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%mdim=mm
						 blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%ndim=nn
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_C%ButterflyKerl(level)%blocks(i-num_rowson,j)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)                      
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i<=num_rowson .and. j>num_colson) then
						 allocate (blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(mm,nn))
						 blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%mdim=mm
						 blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%ndim=nn						 
						 !$omp parallel do default(shared) private(ii,jj)
						 do jj=1, nn
							 do ii=1, mm
								 blocks_B%ButterflyKerl(level)%blocks(i,j-num_colson)%matrix(ii,jj)=blocks_i%ButterflyKerl(level+1)%blocks(i,j)%matrix(ii,jj)
							 enddo
						 enddo
						 !$omp end parallel do
					 elseif (i>num_rowson .and. j>num_colson) then
						 allocate (blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%matrix(mm,nn))
						 blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%mdim=mm
						 blocks_D%ButterflyKerl(level)%blocks(i-num_rowson,j-num_colson)%ndim=nn	 
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
	
	
	! mm1=basis_group(blocks_A%row_group)%tail-basis_group(blocks_A%row_group)%head+1
	! nn1=basis_group(blocks_A%col_group)%tail-basis_group(blocks_A%col_group)%head+1
	! mm2=basis_group(blocks_D%row_group)%tail-basis_group(blocks_D%row_group)%head+1
	! nn2=basis_group(blocks_D%col_group)%tail-basis_group(blocks_D%col_group)%head+1

	
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
