module Butterfly_compress_forward 
use Utilites_randomized 
use H_structure
use element_Z 
contains 
! subroutine Butterfly_compress(blocks,Memory)

   ! use MODULE_FILE
   ! implicit none

    ! integer i, j, level_butterfly, num_blocks, k, attempt
    ! integer group_m, group_n, mm, nn, index_i, index_j, ii, jj
    ! integer level, length_1, length_2, level_blocks, index_ij
    ! integer rank, rankmax, butterflyB_inuse, rank1, rank2
    ! real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    ! real*8 Memory
    ! complex(kind=8) ctemp
	! type(matrixblock)::blocks
	
    ! Memory=0.
	
	! blocks%rankmax = -100000
	! blocks%rankmin = 100000
	
    ! group_m=blocks%row_group ! Note: row_group and col_group interchanged here   
    ! group_n=blocks%col_group
    ! level_blocks=blocks%level
    ! !level_butterfly=Maxlevel-level_blocks
    ! Preset_level_butterfly=Maxlevel_for_blocks
    ! level_butterfly=Preset_level_butterfly-level_blocks
! !     if (Maxlevel-level_blocks<8) then
! !         level_butterfly=Maxlevel-level_blocks
! !     endif
    ! blocks%level_butterfly=level_butterfly
    ! allocate (rankmax_for_butterfly(0:level_butterfly))
    ! rankmax_for_butterfly=0

    ! num_blocks=2**level_butterfly

    ! allocate(blocks%ButterflyU(num_blocks))
    ! allocate(blocks%ButterflyV(num_blocks))
    ! if (level_butterfly/=0) then
        ! allocate(blocks%ButterflyKerl(level_butterfly))
        ! allocate(blocks%ButterflyColSelect(num_blocks,0:level_butterfly))
    ! endif
    
    ! memory_butterfly=0.
    ! do level=0, level_butterfly
        ! index_ij=0
        ! if (level>0) then
            ! blocks%ButterflyKerl(level)%num_row=2**level
            ! blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
            ! allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
        ! endif
        ! do index_i=1, 2**level
            ! do index_j=1, 2**(level_butterfly-level)
                ! call decomposition_UkerlV(index_i,index_j,level,blocks,tolerance)
                ! index_ij=index_ij+1
                ! if (level==0) then
                    ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_ij)%matrix)/1024.0d3
                ! else                    
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
                ! endif
                ! if (level==level_butterfly) then
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_ij)%matrix)/1024.0d3
                ! endif				
            ! enddo
        ! enddo
    ! enddo

	! if(level_blocks==1)then
	! write(*,*)'ha',rankmax_for_butterfly
	! ! stop
	! end if
	! rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),rankmax_of_level(level_blocks))
	
    ! deallocate (rankmax_for_butterfly)
    
    ! if (level_butterfly/=0) then
        ! !$omp parallel do default(shared) private(i,level)
        ! do i=1, num_blocks
            ! do level=0, level_butterfly
                ! deallocate (blocks%ButterflyColSelect(i,level)%select_columns)
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (blocks%ButterflyColSelect)
    ! endif            

    ! Memory=memory_butterfly
    ! !write (*,*) memory_butterfly
    ! !pause

    ! return

! end subroutine butterfly_compress

! subroutine decomposition_UkerlV(index_i,index_j,level,blocks,tolerance)

    ! use MODULE_FILE
    ! implicit none

    ! integer mm, nn, level, level_butterfly
    ! integer i, j, ii, jj, index_i, index_j
    ! real*8 tolerance
	! type(matrixblock)::blocks
	
    ! if (level==0) then
        ! call butterfly_recomposition_FastSampling_initial(index_j,blocks)
    ! else
        ! call butterfly_recomposition_FastSampling(index_i,index_j,level,blocks)
    ! endif

    ! return

! end subroutine decomposition_UkerlV

! subroutine butterfly_recomposition_FastSampling_initial(index_j,blocks)

    ! use MODULE_FILE
    ! ! use lapack95
    ! implicit none

    ! integer i, j, k, level_butterfly, level_blocks, ii, jj
    ! integer group_m, group_n, mm, nn, nn_start, edge_m, edge_n
    ! integer level, length_1, length_2, flag
    ! integer rank, rankmax,rankmax_c,rankmax_r, rankmax_min,index_i, index_j, flag0, rank_new, header_m, header_n
    ! real*8 rate, tolerance, rtemp1, rtemp2, rtemp
    ! complex(kind=8) ctemp

    ! integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    ! integer,allocatable:: select_row_rr(:), select_column_rr(:)
    ! complex(kind=8),allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:)
    ! real*8,allocatable:: Singular(:)
	! integer Nlayer
	! type(matrixblock)::blocks
	
    ! group_m=blocks%row_group  ! Note: row_group and col_group interchanged here    
    ! group_n=blocks%col_group
    ! level_blocks=blocks%level
    ! level_butterfly=blocks%level_butterfly
    ! group_n=group_n*2**level_butterfly-1+index_j

    ! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
    ! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
    ! nn_start=basis_group(group_n)%head-basis_group(blocks%col_group)%head
    ! header_m=basis_group(group_m)%head
    ! header_n=basis_group(group_n)%head

    ! rankmax=rank_approximate_func(group_m, group_n, 1)

    ! if (rankmax>nn) then
        ! rankmax=nn
    ! endif

	! ! modified by Yang, this may improve accuracy as one dimension is constant
	! rankmax=min(mm,nn)	

	! rankmax_r = min(rank_approximate_func(group_m, group_n, 1)*2,mm)
	! rankmax_c = nn
	! rankmax_min = min(rankmax_r,rankmax_c)
	
	! ! write(*,*)rankmax,rank_approximate_func(group_m, group_n, 1),'i'
		
    ! allocate(select_row(rankmax_r),select_column(rankmax_c))


    ! select_row(1)=int(real(mm)/real(rankmax_r)/2.)
    ! if (select_row(1)==0) then
        ! select_row(1)=1
    ! endif
    ! !$omp parallel do default(shared) private(i)
    ! do i=1, rankmax_c
        ! select_column(i)=i
    ! enddo
    ! !$omp end parallel do

	! call linspaceI(1,mm,rankmax_r,select_row(1:mm))	
	
	! Nlayer = 3
	! ! modified by Yang, for 2D problem this can include the touching unknowns
	! do ii =1,Nlayer
		! select_row(ii)=ii
	! end do	
	! do ii =1,Nlayer
		! select_row(rankmax_r-ii+1)=mm-ii+1
	! end do	
	
	! call linspaceI(Nlayer+1,mm-Nlayer,rankmax_r-Nlayer*2,select_row(Nlayer+1:mm-Nlayer))
	
    
    ! allocate (MatrixSubselection(rankmax_r,rankmax_c))

    ! !$omp parallel do default(shared) private(i,j,k,edge_m,edge_n,ctemp)
    ! do j=1, rankmax_c
        ! do i=1, rankmax_r
            ! edge_m=header_m+select_row(i)-1
            ! edge_n=header_n+select_column(j)-1
            ! call element_Zmn(edge_m,edge_n,ctemp)
            ! MatrixSubselection(i,j)=ctemp
        ! enddo
    ! enddo
    ! !$omp end parallel do
    
    ! if (level_butterfly==0) then
        ! allocate (matrix_U(mm,rankmax_c),matrix_V(nn,rankmax_r))
        ! !$omp parallel do default(shared) private(i,j,k,ii,edge_m,edge_n,ctemp)
        ! do j=1, rankmax_c
            ! ii=1
            ! do i=1, mm
                ! if (i==select_row(ii)) then
                    ! matrix_U(i,j)=MatrixSubselection(ii,j)
                    ! if (ii<rankmax_c) then
                        ! ii=ii+1
                    ! endif
                ! else
                    ! edge_m=header_m+i-1
                    ! edge_n=header_n+select_column(j)-1
                    ! call element_Zmn(edge_m,edge_n,ctemp)
                    ! matrix_U(i,j)=ctemp
                ! endif
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(i,j,k,jj,edge_m,edge_n,ctemp)
        ! do i=1, rankmax_r
            ! jj=1
            ! do j=1, nn
                ! if (j==select_column(jj)) then
                    ! matrix_V(j,i)=MatrixSubselection(i,jj)
                    ! if (jj<rankmax_r) then
                        ! jj=jj+1
                    ! endif
                ! else
                    ! edge_m=header_m+select_row(i)-1
                    ! edge_n=header_n+j-1
                    ! call element_Zmn(edge_m,edge_n,ctemp)
                    ! matrix_V(j,i)=ctemp
                ! endif
            ! enddo
        ! enddo
        ! !$omp end parallel do

        ! deallocate (select_column,select_row)
        
        ! allocate (UU(rankmax_r,rankmax_min),VV(rankmax_min,rankmax_c),Singular(rankmax_min))
        ! call gesvd_robust(MatrixSubselection,Singular,UU,VV,rankmax_r,rankmax_c,rankmax_min)
        ! deallocate (MatrixSubselection)
        
        ! flag0=0 ; i=0
        ! do while (flag0==0 .and. i<rankmax_min)
            ! i=i+1
            ! if (Singular(i)/Singular(1)<=SVD_tolerance_forward) then
                ! flag0=1
            ! endif
        ! enddo
        ! if (Singular(i)<1.0e-16) then
            ! i=i-1
        ! endif
        ! rank_new=i
        
		
		! ! ! ! write(*,*)rankmax,rank_new,rankmax_for_butterfly(0)
        ! if (rank_new>rankmax_for_butterfly(0)) then
            ! rankmax_for_butterfly(0)=rank_new
        ! endif
		
		! blocks%rankmax = max(blocks%rankmax,rank_new)
		! blocks%rankmin = min(blocks%rankmin,rank_new)

						
        ! allocate (blocks%ButterflyU(index_j)%matrix(mm,rank_new))
        ! !$omp parallel do default(shared) private(i,j,k,ctemp)
        ! do j=1, rank_new
            ! do i=1, mm
                ! ctemp=0
                ! do k=1, rankmax_c
                    ! ctemp=ctemp+matrix_U(i,k)*conjg(VV(j,k))
                ! enddo
                ! blocks%ButterflyU(index_j)%matrix(i,j)=ctemp
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_U,VV)
        
        ! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank_new))
        ! !$omp parallel do default(shared) private(i,j,k,ctemp)
         ! do j=1, rank_new
            ! do i=1, nn
                ! ctemp=0
                ! do k=1, rankmax_r
                    ! ctemp=ctemp+conjg(UU(k,j))*matrix_V(i,k)
                ! enddo
                ! blocks%ButterflyV(index_j)%matrix(i,j)=ctemp/Singular(j)
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_V,UU,Singular)
    
    ! else
        
        ! allocate (matrixtemp_U(rankmax_r,rankmax_c),matrixtemp_V(rankmax_c,rankmax_r))
        ! allocate (select_column_rr(rankmax_c),select_row_rr(rankmax_r))
        ! ! write(*,*)rankmax_r,rankmax_c
		! call ACA_SubsetSelection(select_column_rr,select_row_rr,rankmax_r,rankmax_c,rank_new)
         ! ! write(*,*)rankmax_r,rankmax_c,'dddd'
		! ! write(*,*)rank_new
		
        ! !write (*,*) index_j,rankmax,rank_new
        
        ! !if (rank_new<rankmax) then
        ! !    !$omp parallel do default(shared) private(i)
        ! !    do i=rank_new+1, rankmax
        ! !        column_pivot(i)=rankmax+1
        ! !        row_pivot(i)=rankmax+1
        ! !    enddo
        ! !    !$omp end parallel do
        ! !endif
        ! !
        ! !allocate (select_column_rr(rank_new),select_row_rr(rank_new))
        ! !
        ! !do j=1, rank_new
        ! !    jj=minloc(column_pivot,1)
        ! !    select_column_rr(j)=column_pivot(jj)
        ! !    column_pivot(jj)=rankmax+1
        ! !    jj=minloc(row_pivot,1)
        ! !    select_row_rr(j)=row_pivot(jj)
        ! !    row_pivot(jj)=rankmax+1
        ! !enddo
        ! !deallocate (column_pivot,row_pivot)
            
        ! deallocate (matrixtemp_U,matrixtemp_V)
        
		! ! ! write(*,*)rankmax,rank_new,rankmax_for_butterfly(0)
		! ! if(level_blocks==1)write(*,*)rank_new
        ! if (rank_new>rankmax_for_butterfly(0)) then
            ! rankmax_for_butterfly(0)=rank_new
        ! endif
		
		! blocks%rankmax = max(blocks%rankmax,rank_new)
		! blocks%rankmin = min(blocks%rankmin,rank_new)

        ! allocate (blocks%ButterflyColSelect(index_j,0)%select_columns(rank_new))
        ! !$omp parallel do default(shared) private(j)
        ! do j=1, rank_new
            ! blocks%ButterflyColSelect(index_j,0)%select_columns(j)=select_column(select_column_rr(j))
        ! enddo
        ! !$omp end parallel do
                
        ! allocate (matrix_little(rank_new,rank_new))
        
        ! !$omp parallel do default(shared) private(ii,jj)
        ! do jj=1, rank_new
            ! do ii=1, rank_new
                ! matrix_little(ii,jj)=MatrixSubselection(select_row_rr(ii),select_column_rr(jj))
            ! enddo
        ! enddo
        ! !$omp end parallel do
        
        ! allocate (matrix_V(rank_new,nn))
        ! !$omp parallel do default(shared) private(i,j,jj,edge_m,edge_n,ctemp)
        ! do i=1, rank_new
            ! jj=1
            ! do j=1, nn
                ! if (j==select_column(jj)) then
                    ! matrix_V(i,j)=MatrixSubselection(select_row_rr(i),jj)
                    ! if (jj<rankmax_c) then
                        ! jj=jj+1
                    ! endif
                ! else
                    ! edge_m=header_m+select_row(select_row_rr(i))-1
                    ! edge_n=header_n+j-1
                    ! call element_Zmn(edge_m,edge_n,ctemp)
                    ! matrix_V(i,j)=ctemp
                ! endif
            ! enddo
        ! enddo
        ! !$omp end parallel do
        
        ! deallocate (MatrixSubselection)
        

		! allocate(matrix_little_inv(rank_new,rank_new))
		! allocate (matrix_V_tmp(rank_new,nn))
		! call GeneralInverse(rank_new,rank_new,matrix_little,matrix_little_inv,ACA_tolerance_forward)
		! call gemm_omp(matrix_little_inv,matrix_V,matrix_V_tmp,rank_new,rank_new,nn)
		! matrix_V = matrix_V_tmp
		! deallocate(matrix_little_inv)
		! deallocate(matrix_V_tmp)		
		
		
        ! ! ! allocate (column_pivot(rank_new))
        ! ! ! call getrff90(matrix_little,column_pivot)
        ! ! ! call getrsf90(matrix_little,column_pivot,matrix_V,'N')
        ! ! ! deallocate (column_pivot,matrix_little)
        
		
        ! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank_new))
        ! !$omp parallel do default(shared) private(i,j)
         ! do j=1, rank_new
            ! do i=1, nn
                ! blocks%ButterflyV(index_j)%matrix(i,j)=matrix_v(j,i)
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_V)
        
        ! deallocate (select_row, select_row_rr, select_column, select_column_rr)    
        
    ! endif

    ! return

! end subroutine butterfly_recomposition_FastSampling_initial
 
! subroutine butterfly_recomposition_FastSampling(index_i,index_j,level,blocks)

    ! use MODULE_FILE
    ! ! use lapack95
    ! implicit none

    ! integer i, j, level_butterfly, level_blocks, num_groupm, num_groupn, k, ii, jj, mmm, nnn, nnn1
    ! integer group_m, group_n, mm, nn, index_m, index_n, option, mm_start, nn1, nn2, header_m, header_n
    ! integer level, length_1, length_2, butterflyB_inuse, index_ii, index_jj, header_n1, header_n2
    ! integer rank, rankmax,rankmax_c,rankmax_r,rankmax_min,index_i, index_j, index_ij, index_iijj, flag, flag0, rank_new, edge_m, edge_n
    ! real*8 rate, tolerance, rtemp1, rtemp2, rtemp
    ! complex(kind=8) ctemp

    ! integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    ! integer,allocatable:: select_row_rr(:), select_column_rr(:)
    ! complex(kind=8),allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:)
    ! real*8,allocatable:: Singular(:)
	! type(matrixblock)::blocks
		
    ! level_butterfly=blocks%level_butterfly
    ! group_m=blocks%row_group    ! Note: row_group and col_group interchanged here   
    ! group_n=blocks%col_group    
    ! ii=basis_group(group_m)%head
    ! jj=basis_group(group_n)%head
    ! group_m=group_m*2**level-1+index_i
    ! group_n=group_n*2**(level_butterfly-level)-1+index_j
    ! header_m=basis_group(group_m)%head
    ! header_n1=basis_group(group_n)%head
    ! header_n2=basis_group(2*group_n+1)%head
    ! nnn1=basis_group(2*group_n)%tail-basis_group(2*group_n)%head+1

    ! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
    ! mmm=basis_group(int(group_m/2))%tail-basis_group(int(group_m/2))%head+1
    ! index_ii=int((index_i+1)/2) ; index_jj=2*index_j-1
    ! num_groupn=2**(level_butterfly-level+1)
    ! index_iijj=(index_ii-1)*num_groupn+index_jj
    ! nn1=size(blocks%ButterflyColSelect(index_iijj,level-1)%select_columns,1)
    ! nn2=size(blocks%ButterflyColSelect(index_iijj+1,level-1)%select_columns,1)
    ! nn=nn1+nn2
    ! num_groupn=2**(level_butterfly-level)
    ! index_ij=(index_i-1)*num_groupn+index_j

    ! rankmax=rank_approximate_func(group_m, group_n, 1)

    ! if (rankmax>nn) then
        ! rankmax=nn
    ! endif

	! ! modified by Yang, this may improve accuracy as one dimension is constant
	! rankmax=min(mm,nn)	

	! rankmax_r = min(rank_approximate_func(group_m, group_n, 1)*2,mm)
	! rankmax_c = nn
	! rankmax_min = min(rankmax_r,rankmax_c)
	
    ! ! if (rank_control_forward/=0) then
        ! ! if (rankmax_for_butterfly(level-1)<rankmax) then
            ! ! rankmax=rankmax_for_butterfly(level-1)
        ! ! endif
    ! ! endif

	! if(rankmax==1)write(*,*)group_m,group_n,rank_approximate_func(group_m, group_n, 1),mm,nn
    ! allocate(select_row(rankmax_r),select_column(rankmax_c))

    ! !$omp parallel do default(shared) private(i)
    ! do i=1, rankmax_c
        ! select_column(i)=i
    ! enddo
    ! !$omp end parallel do	
	
	
	! call linspaceI(1,mm,rankmax_r,select_row(1:mm))
	
	
	! ! modified by Yang, for 2D problem this can include the touching unknowns
	! select_row(1)=1
	! select_row(2)=2
	! select_row(3)=3	
	
	! select_row(rankmax_r)=mm
	! select_row(rankmax_r-1)=mm-1
	! select_row(rankmax_r-2)=mm-2
	! call linspaceI(4,mm-3,rankmax_r-6,select_row(4:mm-3))
	
	
	
    ! if (mod(index_i,2)/=0) then
        ! mm_start=0
    ! else
        ! mm_start=mmm-mm
    ! endif

    ! allocate (MatrixSubselection(rankmax_r,rankmax_c))

    ! !$omp parallel do default(shared) private(i,j,k,edge_m,edge_n,ctemp)
    ! do j=1, rankmax_c
        ! do i=1, rankmax_r
            ! if (select_column(j)<=nn1) then
                ! edge_m=select_row(i)+header_m-1
                ! edge_n=blocks%ButterflyColSelect(index_iijj,level-1)%select_columns(select_column(j))+header_n1-1
                ! call element_Zmn(edge_m,edge_n,ctemp)
                ! MatrixSubselection(i,j)=ctemp
            ! else
                ! edge_m=select_row(i)+header_m-1
                ! edge_n=blocks%ButterflyColSelect(index_iijj+1,level-1)%select_columns(select_column(j)-nn1)+header_n2-1
                ! call element_Zmn(edge_m,edge_n,ctemp)
                ! MatrixSubselection(i,j)=ctemp
            ! endif
			! ! if(level==1)then
				! ! write(*,*)nn1,nn2,j,i,edge_m-basis_group(blocks%row_group)%head+1,edge_n-basis_group(blocks%col_group)%head+1
			! ! end if
        ! enddo
    ! enddo
    ! !$omp end parallel do
    
    ! allocate (matrixtemp_U(rankmax_r,rankmax_c),matrixtemp_V(rankmax_c,rankmax_r))
    ! allocate (select_column_rr(rankmax_c),select_row_rr(rankmax_r))
    ! call ACA_SubsetSelection(select_column_rr,select_row_rr,rankmax_r,rankmax_c,rank_new)
    
    ! ! ! if(level==1)then
		! ! ! write(*,*)rank_new,rankmax
		! ! ! ! write(*,*)MatrixSubselection
	! ! ! end if
	
    ! !write (*,*) index_j,index_i,rankmax,rank_new
    
    ! !if (rank_new<rankmax) then
    ! !    !$omp parallel do default(shared) private(i)
    ! !    do i=rank_new+1, rankmax
    ! !        column_pivot(i)=rankmax+1
    ! !        row_pivot(i)=rankmax+1
    ! !    enddo
    ! !    !$omp end parallel do
    ! !endif
    ! !
    ! !allocate (select_column_rr(rank_new),select_row_rr(rank_new))
    ! !
    ! !do j=1, rank_new
    ! !    jj=minloc(column_pivot,1)
    ! !    select_column_rr(j)=column_pivot(jj)
    ! !    column_pivot(jj)=rankmax+1
    ! !    jj=minloc(row_pivot,1)
    ! !    select_row_rr(j)=row_pivot(jj)
    ! !    row_pivot(jj)=rankmax+1
    ! !enddo
    ! !deallocate (column_pivot,row_pivot)
        
    ! deallocate (matrixtemp_U,matrixtemp_V)

    ! ! if (rank_control_forward/=0) then
        ! if (rank_new>rankmax_for_butterfly(level)) then
            ! rankmax_for_butterfly(level)=rank_new
        ! endif
		
		! blocks%rankmax = max(blocks%rankmax,rank_new)
		! blocks%rankmin = min(blocks%rankmin,rank_new)
		
    ! ! endif
    
    ! allocate (blocks%ButterflyColSelect(index_ij,level)%select_columns(rank_new))
    ! !$omp parallel do default(shared) private(j)
    ! do j=1, rank_new
        ! if (select_column(select_column_rr(j))<=nn1) then
            ! blocks%ButterflyColSelect(index_ij,level)%select_columns(j)=blocks%ButterflyColSelect(index_iijj,level-1)%select_columns(select_column(select_column_rr(j)))
        ! else
            ! blocks%ButterflyColSelect(index_ij,level)%select_columns(j)=blocks%ButterflyColSelect(index_iijj+1,level-1)%select_columns(select_column(select_column_rr(j))-nn1)+nnn1
        ! endif
    ! enddo
    ! !$omp end parallel do                
            
    ! allocate (matrix_little(rank_new,rank_new))
    
    ! !$omp parallel do default(shared) private(ii,jj)
    ! do jj=1, rank_new
        ! do ii=1, rank_new
            ! matrix_little(ii,jj)=MatrixSubselection(select_row_rr(ii),select_column_rr(jj))
        ! enddo
    ! enddo
    ! !$omp end parallel do
    
    ! if (level==level_butterfly) then
        ! allocate (blocks%ButterflyU(index_ij)%matrix(mm,rank_new)) 
        ! !$omp parallel do default(shared) private(i,j,ii,edge_m,edge_n,ctemp)
        ! do j=1, rank_new
            ! ii=1
            ! do i=1, mm
                ! if (i==select_row(ii)) then
                    ! blocks%ButterflyU(index_ij)%matrix(i,j)=MatrixSubselection(ii,select_column_rr(j))
                    ! if (ii<rankmax_r) then
                        ! ii=ii+1
                    ! endif
                ! else
                    ! edge_m=i+header_m-1
                    ! edge_n=blocks%ButterflyColSelect(index_ij,level)%select_columns(j)+header_n1-1 
                    ! call element_Zmn(edge_m,edge_n,ctemp)                    
                    ! blocks%ButterflyU(index_ij)%matrix(i,j)=ctemp
                ! endif
            ! enddo
        ! enddo
        ! !$omp end parallel do
    ! endif
    
    ! allocate (matrix_V(rank_new,nn))
    ! !$omp parallel do default(shared) private(i,j,jj,edge_m,edge_n,ctemp)
    ! do i=1, rank_new
        ! jj=1
        ! do j=1, nn
            ! if (j==select_column(jj)) then
                ! matrix_V(i,j)=MatrixSubselection(select_row_rr(i),jj)
                ! if (jj<rankmax_c) then
                    ! jj=jj+1
                ! endif            
            ! else
                ! if (j<=nn1) then
                    ! edge_m=select_row(select_row_rr(i))+header_m-1
                    ! edge_n=blocks%ButterflyColSelect(index_iijj,level-1)%select_columns(j)+header_n1-1
                ! else
                    ! edge_m=select_row(select_row_rr(i))+header_m-1
                    ! edge_n=blocks%ButterflyColSelect(index_iijj+1,level-1)%select_columns(j-nn1)+header_n2-1
                ! endif
                ! call element_Zmn(edge_m,edge_n,ctemp)                    
                ! matrix_V(i,j)=ctemp
            ! endif
        ! enddo
    ! enddo
    ! !$omp end parallel do
    
    ! deallocate (MatrixSubselection)
    
	! allocate(matrix_little_inv(rank_new,rank_new))
	! allocate (matrix_V_tmp(rank_new,nn))
	! call GeneralInverse(rank_new,rank_new,matrix_little,matrix_little_inv,ACA_tolerance_forward)
	! call gemm_omp(matrix_little_inv,matrix_V,matrix_V_tmp,rank_new,rank_new,nn)
	! matrix_V = matrix_V_tmp
	! deallocate(matrix_little_inv)
	! deallocate(matrix_V_tmp)
	
	
    ! ! allocate (column_pivot(rank_new))
    ! ! call getrff90(matrix_little,column_pivot)
    ! ! call getrsf90(matrix_little,column_pivot,matrix_V,'N')
    ! ! deallocate (column_pivot,matrix_little)
	
    
    ! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank_new,nn1))
    ! !$omp parallel do default(shared) private(i,j)
     ! do j=1, nn1
        ! do i=1, rank_new
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=matrix_V(i,j)
        ! enddo
    ! enddo
    ! !$omp end parallel do
    ! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank_new,nn2))
    ! !$omp parallel do default(shared) private(i,j)
     ! do j=1, nn2
        ! do i=1, rank_new
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=matrix_V(i,j+nn1)
        ! enddo
    ! enddo
    ! !$omp end parallel do
    ! deallocate (matrix_V)                    
    
    ! deallocate (select_row, select_row_rr, select_column, select_column_rr)

    ! return

! end subroutine butterfly_recomposition_FastSampling



subroutine Bplus_compress_N15(bplus,Memory)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    type(blockplus)::bplus
	integer:: ii,ll,bb
    real*8 Memory,rtemp	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	
	Memory = 0
	
	do ll=1,bplus%Lplus
		bplus%LL(ll)%rankmax=0
																						  
						
														
																   
																	   
		do bb = 1,bplus%LL(ll)%Nbound
			! write(*,*)bplus%level,ll,bb
			if(bplus%LL(ll+1)%Nbound==0)then
				
				bplus%LL(ll)%matrices_block(bb)%level_butterfly = int((Maxlevel_for_blocks-bplus%LL(ll)%matrices_block(bb)%level)/2)*2
				if(TwoLayerOnly==1 .and. bplus%Lplus==2)bplus%LL(ll)%matrices_block(bb)%level_butterfly = 0
				call Butterfly_compress_N15(bplus%LL(ll)%matrices_block(bb),rtemp)
				call Butterfly_sym2asym(bplus%LL(ll)%matrices_block(bb))
				Memory = Memory + rtemp
			else 		
				level_butterfly = int((Maxlevel_for_blocks - bplus%LL(ll)%matrices_block(1)%level)/2)*2 
				level_BP = bplus%level			
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				groupm_start=bplus%LL(ll)%matrices_block(1)%row_group*2**levelm		
				Nboundall = 2**(bplus%LL(ll)%matrices_block(1)%level+levelm-level_BP)			
				call Butterfly_compress_N15_withoutBoundary(bplus%LL(ll)%matrices_block(bb),bplus%LL(ll+1)%boundary_map,Nboundall,groupm_start, rtemp)
				call Butterfly_sym2asym(bplus%LL(ll)%matrices_block(bb))
				Memory = Memory + rtemp
			end if	
			bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax,bplus%LL(ll)%matrices_block(bb)%rankmax)			
		end do
	end do
	! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'
	
    return

end subroutine Bplus_compress_N15



subroutine Butterfly_compress_N15_withoutBoundary(blocks,boundary_map,Nboundall, groupm_start,Memory)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none
   
	type(matrixblock)::blocks
	integer Nboundall
	integer boundary_map(Nboundall)
	integer groupm_start
	
    integer blocks_idx, i, j, level_butterfly,level_butterflyL,level_butterflyR, num_blocks, k, attempt
    integer group_m, group_n, mm, nn,mn, index_i, index_j,index_i_m, index_j_m,index_i_loc, index_j_loc, index_ij_loc,ii, jj,nn1,nn2,mm1,mm2,idxs_m,idxs_n
    integer level,levelm,level_loc, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    integer rank, rankmax, butterflyB_inuse, rank1, rank2,rmax
    real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real*8 Memory
    complex(kind=8) ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	complex(kind=8), allocatable :: tau_Q(:)
	real*8,allocatable :: Singular(:)
	integer flag,tmpi,tmpj
	real*8:: a,b,n1,n2
	real*8:: maxvalue(1:20)
	real*8:: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row
	
	integer cnt_tmp,rankFar,rankNear
	
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	ForwardSymmetricFlag = 1	
	maxvalue = 0
	minvalue = 10000
    Memory=0

    level_blocks=blocks%level
    !level_butterfly=Maxlevel-level_blocks
    level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2
																  
	
	call assert(level_butterfly>=2,'level_butterfly should be at least 2')
	

!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
	blocks%rankmax = -100000
	blocks%rankmin = 100000
	
    blocks%level_butterfly=level_butterfly


    num_blocks=2**level_butterfly

    allocate(blocks%ButterflyU(num_blocks))
    allocate(blocks%ButterflyV(num_blocks))
    if (level_butterfly/=0) then
        allocate(blocks%ButterflyKerl(level_butterfly))
    endif
    
    memory_butterfly=0.
    

	do level=1, level_butterfly
		blocks%ButterflyKerl(level)%num_row=2**level
		blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
		allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
	end do
	
	levelm = ceiling_safe(dble(level_butterfly)/2d0)
	level_butterflyL = level_butterfly-levelm
	level_butterflyR = level_butterfly-level_butterflyL
	
	allocate (rankmax_for_butterfly(-level_butterflyR:level_butterflyL))
	rankmax_for_butterfly=0		
	allocate (rankmin_for_butterfly(-level_butterflyR:level_butterflyL))
	rankmin_for_butterfly=100000	
	
	allocate(blocks%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))
	
	
	! construct the middle level and the left half
	
	do index_i_m=1, 2**levelm
		level_loc = 0
		index_i_loc = 1
		allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
		
		do index_j_m=1, 2**(level_butterfly-levelm)	
			index_j_loc = index_j_m
			group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
			group_n=blocks%col_group    
			group_m=group_m*2**levelm-1+index_i_m
			group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
		
			mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
			nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
			idxs_m = basis_group(group_m)%head
			idxs_n = basis_group(group_n)%head
			
			rmax = min(500,min(mm,nn))
			allocate(matU(mm,rmax))
			allocate(matV(rmax,nn))
			allocate(Singular(rmax))
			! call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)					
			! rank = min(rank,37)
			
			
			
			! if(near_or_far(group_m,group_n,2d0))write(*,*)near_or_far(group_m,group_n,2d0),rank
			
			! ! write(*,*)near_or_far(group_m,group_n,1.5d0),rank
			! if(.not. near_or_far(group_m,group_n,1.5d0))then
				! cnt_tmp = cnt_tmp+1
				! rankNear = max(rankNear,rank)
			! else 
				! rankFar= max(rankFar,rank)
			! end if
			

			
			
			if(boundary_map(group_m-groupm_start+1)==group_n)then			
			! if(.not. near_or_far(group_m,group_n,2d0))then
				rank = 1
				Singular(1:rank) = 0
				! if(blocks==342)write(111,*)Singular(1:rank)
				rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
				rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
				
				allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
				do ii=1,rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 0d0
				end do	
			else 
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)	
				! if(rank==53)then
					! write(*,*)group_m,group_n,boundary_map(group_m-groupm_start+1)
					! stop
				! end if
						
				! if(blocks==342)write(111,*)Singular(1:rank)
				rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
				rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
				
				allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
				do ii=1,rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 1d0/Singular(ii)
				end do
			end if
			

			
			allocate(mat_tmp(mm,rank))
			!$omp parallel do default(shared) private(i,j,k,ctemp)
			do j=1, rank
				do i=1, mm
					mat_tmp(i,j)=matU(i,j)*Singular(j)
				enddo
			enddo
			!$omp end parallel do
			
			mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
			allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
			allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
			
			ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
			ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)				
			
			deallocate(matU,matV,Singular,mat_tmp)
		end do
		
		n1 = OMP_get_wtime()
		
		do level_loc = 1,level_butterflyL
			level = level_loc+levelm
			
			if(level_loc/=level_butterflyL)then
				allocate(ButterflyP%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
			end if		
			!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)	
			do index_ij_loc = 1, 2**level_butterflyL
				index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,SVD_tolerance_forward,ButterflyP_old,ButterflyP)				
			enddo
			!$omp end parallel do	
			
			do index_ij_loc = 1, 2**level_butterflyL
				index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				index_j = index_j_loc
				index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
				
				rank = size(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
	
				rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
				rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
				
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
				if (level_loc==level_butterflyL) then
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
				endif
			enddo						
			
			if(level_loc/=level_butterflyL)then
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_loc)
						do jj =1,2**(level_butterflyL-level_loc+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if
				allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
				do ii = 1,2**(level_loc+1)
					do jj =1,2**(level_butterflyL-level_loc)
						mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
						nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
						ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
						deallocate(ButterflyP%blocks(ii,jj)%matrix)
					end do
				end do
				deallocate(ButterflyP%blocks)
			else 
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_loc)
						do jj =1,2**(level_butterflyL-level_loc+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if					
			end if				
			
		end do	
		n2 = OMP_get_wtime()
		time_tmp = time_tmp + n2 - n1
		
	enddo

	
	! write(*,*)'stat:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	


	! construct the the right half
	do index_j_m=1, 2**(level_butterfly-levelm)	
		level_loc = 0
		index_j_loc = 1
		allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
		do index_i_m=1, 2**levelm			
			index_i_loc = index_i_m
			group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
			group_n=blocks%col_group    
			group_m=group_m*2**levelm-1+index_i_m
			group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
		
			mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
			nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
			idxs_m = basis_group(group_m)%head
			idxs_n = basis_group(group_n)%head
			
			rmax = min(500,min(mm,nn))
			allocate(matU(mm,rmax))
			allocate(matV(rmax,nn))
			allocate(Singular(rmax))
			
			! rank = min(rank,37)
			
			if(boundary_map(group_m-groupm_start+1)==group_n)then	
			! if(.not. near_or_far(group_m,group_n,2d0))then
				rank = 1
				Singular(1:rank) = 0
				matV(1:rank,1:nn)=0
				! if(blocks==342)write(111,*)Singular(1:rank)
			else 
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)	
			end if
			
			
			
			allocate(mat_tmp(rank,nn))
			!$omp parallel do default(shared) private(i,j,k,ctemp)
			do j=1, nn
				do i=1, rank
					mat_tmp(i,j)=matV(i,j)*Singular(i)
				enddo
			enddo
			!$omp end parallel do
			
			nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
			
			ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
			ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)				
			
			deallocate(matU,matV,Singular,mat_tmp)
		end do
		
		n1 = OMP_get_wtime()
		do level_loc = 1,level_butterflyR
			level = levelm+1-level_loc
			
			if(level_loc/=level_butterflyR)then
				allocate(ButterflyP%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
			end if

			
			!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
			do index_ij_loc = 1, 2**level_butterflyR
				index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
				call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,SVD_tolerance_forward,ButterflyP_old,ButterflyP)				
			enddo
			!$omp end parallel do	
			
			do index_ij_loc = 1, 2**level_butterflyR
				index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
				index_i = index_i_loc
				index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
				rank = size(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)				
				rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
				rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
								
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
				if (level_loc==level_butterflyR) then
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
				endif					
			enddo				
			
			if(level_loc/=level_butterflyR)then
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_butterflyR-level_loc+1)
						do jj =1,2**(level_loc)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if
				allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
				do ii = 1,2**(level_butterflyR-level_loc)
					do jj =1,2**(level_loc+1)
						mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
						nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
						ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
						deallocate(ButterflyP%blocks(ii,jj)%matrix)
					end do
				end do
				deallocate(ButterflyP%blocks)
			else 
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_butterflyR-level_loc+1)
						do jj =1,2**(level_loc)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if					
			end if				
			
		end do	
		n2 = OMP_get_wtime()
		time_tmp = time_tmp + n2 - n1

	enddo		

	! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
	rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),rankmax_of_level(level_blocks))
	! if(blocks==342)then
	! write(*,*)rankmax_for_butterfly
	! write(*,*)rankmin_for_butterfly
	! end if
	
	! write(*,*)'max value: ',maxvalue(1:9)
	! write(*,*)'min value: ',minvalue(1:9)
	
	
    deallocate (rankmax_for_butterfly)
    deallocate (rankmin_for_butterfly)
    
    Memory=memory_butterfly
    !write (*,*) memory_butterfly
    !pause

    return

end subroutine Butterfly_compress_N15_withoutBoundary





subroutine Butterfly_compress_N15_withoutBoundary_givenfullmat(blocks,boundary_map,Nboundall, groupm_start,Memory,idx_m_ref,idx_n_ref)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none
   
	type(matrixblock)::blocks
	integer Nboundall
	integer boundary_map(Nboundall)
	integer groupm_start
	
    integer blocks_idx, i, j, level_butterfly,level_butterflyL,level_butterflyR, num_blocks, k, attempt
    integer group_m, group_n, mm, nn,mn, index_i, index_j,index_i_m, index_j_m,index_i_loc, index_j_loc,index_ij_loc, ii, jj,nn1,nn2,mm1,mm2,idxs_m,idxs_n
    integer level,levelm,level_loc, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    integer rank, rankmax, butterflyB_inuse, rank1, rank2,rmax, ranktotL,ranktotR
    real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real*8 Memory
    complex(kind=8) ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	complex(kind=8), allocatable :: tau_Q(:)
	real*8,allocatable :: Singular(:)
	integer flag,tmpi,tmpj,mt,nt
	real*8:: a,b,n1,n2
	real*8:: maxvalue(1:20)
	real*8:: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row
	
	integer cnt_tmp,rankFar,rankNear,idx_m_ref,idx_n_ref
	
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	ForwardSymmetricFlag = 1	
	maxvalue = 0
	minvalue = 10000
    Memory=0

    level_blocks=blocks%level
    !level_butterfly=Maxlevel-level_blocks
    level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2
																  
	
	call assert(level_butterfly>=2,'level_butterfly should be at least 2')
	

!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
	blocks%rankmax = -100000
	blocks%rankmin = 100000
	
    blocks%level_butterfly=level_butterfly


    num_blocks=2**level_butterfly

    allocate(blocks%ButterflyU(num_blocks))
    allocate(blocks%ButterflyV(num_blocks))
    if (level_butterfly/=0) then
        allocate(blocks%ButterflyKerl(level_butterfly))
    endif
    
    memory_butterfly=0.
    

	do level=1, level_butterfly
		blocks%ButterflyKerl(level)%num_row=2**level
		blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
		allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
	end do
	
	levelm = ceiling_safe(dble(level_butterfly)/2d0)
	level_butterflyL = level_butterfly-levelm
	level_butterflyR = level_butterfly-level_butterflyL
	
	
	allocate(blocks%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))
	
	ranktotL = 0
	! construct the middle level and the left half
	
	do index_i_m=1, 2**levelm
		level_loc = 0
		index_i_loc = 1
		allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
		
		do index_j_m=1, 2**(level_butterfly-levelm)	
			index_j_loc = index_j_m
			group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
			group_n=blocks%col_group    
			group_m=group_m*2**levelm-1+index_i_m
			group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
		
			mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
			nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
			idxs_m = basis_group(group_m)%head
			idxs_n = basis_group(group_n)%head
			
			rmax = min(500,min(mm,nn))
			allocate(matU(mm,rmax))
			allocate(matV(rmax,nn))
			allocate(Singular(rmax))
			! call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)					
			! rank = min(rank,37)
			
			
			
			! if(near_or_far(group_m,group_n,2d0))write(*,*)near_or_far(group_m,group_n,2d0),rank
			
			! ! write(*,*)near_or_far(group_m,group_n,1.5d0),rank
			! if(.not. near_or_far(group_m,group_n,1.5d0))then
				! cnt_tmp = cnt_tmp+1
				! rankNear = max(rankNear,rank)
			! else 
				! rankFar= max(rankFar,rank)
			! end if
			

			
			
			if(boundary_map(group_m-groupm_start+1)==group_n)then			
			! if(.not. near_or_far(group_m,group_n,2d0))then
				rank = 1
				Singular(1:rank) = 0
				! if(blocks==342)write(111,*)Singular(1:rank)
				
				allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
				do ii=1,rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 0d0
				end do	
			else 
				call ACA_CompressionForward_givenfullmat(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_factor*0.1,SVD_tolerance_factor,idx_m_ref,idx_n_ref)	

				
				allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
				do ii=1,rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 1d0/Singular(ii)
				end do
			end if
			ranktotL = ranktotL + rank
			blocks%rankmax = min(blocks%rankmax,rank)
			blocks%rankmin = min(blocks%rankmin,rank)
			
			! if(blocks%row_group==34 .and. blocks%col_group==40 .and. index_i_m==1 .and. index_j_m==1)then
				! mt = size(matsub_glo,1)
				! nt = size(matsub_glo,1)
				! write(*,*)index_i_m,index_j_m,rank,'L',idxs_m,idxs_n,mm,nn,rmax,idx_m_ref,idx_n_ref,mt,nt,fnorm(matsub_glo,mt,nt)
			! end if
			
			allocate(mat_tmp(mm,rank))
			!$omp parallel do default(shared) private(i,j,k,ctemp)
			do j=1, rank
				do i=1, mm
					mat_tmp(i,j)=matU(i,j)*Singular(j)
				enddo
			enddo
			!$omp end parallel do
			
			mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
			allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
			allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
			
			ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
			ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)				
			
			deallocate(matU,matV,Singular,mat_tmp)
		end do
		n1 = OMP_get_wtime()
		do level_loc = 1,level_butterflyL
			level = level_loc+levelm
			
			if(level_loc/=level_butterflyL)then
				allocate(ButterflyP%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
			end if

			
			!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)	
			do index_ij_loc = 1, 2**level_butterflyL
				index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,SVD_tolerance_factor,ButterflyP_old,ButterflyP)				
			enddo
			!$omp end parallel do	
			
			do index_ij_loc = 1, 2**level_butterflyL
				index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				index_j = index_j_loc
				index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
				
				rank = size(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
	
				! rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
				! rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
				
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
				if (level_loc==level_butterflyL) then
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
				endif
			enddo					
			
			
			! do index_i_loc=1, 2**level_loc
				! do index_j_loc=1, 2**(level_butterflyL-level_loc)
					! index_j = index_j_loc
					! index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
					! group_m=blocks%row_group   ! Note: row_group and col_group interchanged here    
					! group_m=group_m*2**level-1+index_i
					! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
								
					! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
					! if(size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,1)/=mm)then
						! write(*,*)'mm incorrect'
						! stop
					! end if
					! nn1=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,2)
					! nn2=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix,2)
					! nn = nn1+nn2
					

					! allocate(QQ(mm,nn))
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn1
						! do i=1, mm
							! QQ(i,j)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn2
						! do i=1, mm
							! QQ(i,j+nn1)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					
					
					
					! mn=min(mm,nn)
					! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
					! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_factor,rank)			
					! ! rank = min(rank,37)

					! blocks%rankmax = max(blocks%rankmax,rank)
					! blocks%rankmin = min(blocks%rankmin,rank)

					! allocate(mat_tmp(mm,rank))
					! !$omp parallel do default(shared) private(i,j,k,ctemp)
					! do j=1, rank
						! do i=1, mm
							! mat_tmp(i,j)=UU(i,j)*Singular(j)
						! enddo
					! enddo
					! !$omp end parallel do							
					


					! if(level_loc/=level_butterflyL)then
						! mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
						! allocate(ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
						! allocate(ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
						
						! ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
						! ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)						
					! else 
						! allocate (blocks%ButterflyU(index_i)%matrix(mm,rank))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, mm
								! blocks%ButterflyU(index_i)%matrix(i,j)=mat_tmp(i,j)
								! ! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do

					! end if		

					! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1))
					! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2))
					
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, rank
						! do j=1, nn1
							! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=VV(i,j)
							! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, rank
						! do j=1, nn2
							! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
							! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
						! enddo
					! enddo
					! !$omp end parallel do						
					
					! deallocate(QQ,UU,VV,Singular,mat_tmp)
				  
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
					! if (level_loc==level_butterflyL) then
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
					! endif
				! enddo
			! enddo
			
			if(level_loc/=level_butterflyL)then
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_loc)
						do jj =1,2**(level_butterflyL-level_loc+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if
				allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
				do ii = 1,2**(level_loc+1)
					do jj =1,2**(level_butterflyL-level_loc)
						mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
						nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
						ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
						deallocate(ButterflyP%blocks(ii,jj)%matrix)
					end do
				end do
				deallocate(ButterflyP%blocks)
			else 
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_loc)
						do jj =1,2**(level_butterflyL-level_loc+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if					
			end if				
			
		end do	
		n2 = OMP_get_wtime()
		time_tmp = time_tmp + n2 - n1

	enddo

	
	! write(*,*)'stat:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	

	ranktotR=0
	! construct the the right half
	do index_j_m=1, 2**(level_butterfly-levelm)	
		level_loc = 0
		index_j_loc = 1
		allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
		do index_i_m=1, 2**levelm			
			index_i_loc = index_i_m
			group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
			group_n=blocks%col_group    
			group_m=group_m*2**levelm-1+index_i_m
			group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
		
			mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
			nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
			idxs_m = basis_group(group_m)%head
			idxs_n = basis_group(group_n)%head
			
			rmax = min(500,min(mm,nn))
			allocate(matU(mm,rmax))
			allocate(matV(rmax,nn))
			allocate(Singular(rmax))
			
			! rank = min(rank,37)
			
			if(boundary_map(group_m-groupm_start+1)==group_n)then	
			! if(.not. near_or_far(group_m,group_n,2d0))then
				rank = 1
				Singular(1:rank) = 0
				matV(1:rank,1:nn)=0
				! if(blocks==342)write(111,*)Singular(1:rank)
			else 
				call ACA_CompressionForward_givenfullmat(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_factor*0.1,SVD_tolerance_factor,idx_m_ref,idx_n_ref)	
			end if
			ranktotR = ranktotR + rank
			blocks%rankmax = min(blocks%rankmax,rank)
			blocks%rankmin = min(blocks%rankmin,rank)			
			! if(blocks%row_group==34 .and. blocks%col_group==40 .and. index_i_m==1 .and. index_j_m==1)then
				! mt = size(matsub_glo,1)
				! nt = size(matsub_glo,1)
				! write(*,*)index_i_m,index_j_m,rank,'R',idxs_m,idxs_n,mm,nn,rmax,idx_m_ref,idx_n_ref,mt,nt,fnorm(matsub_glo,mt,nt)
			! end if
			
			allocate(mat_tmp(rank,nn))
			!$omp parallel do default(shared) private(i,j,k,ctemp)
			do j=1, nn
				do i=1, rank
					mat_tmp(i,j)=matV(i,j)*Singular(i)
				enddo
			enddo
			!$omp end parallel do
			
			nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
			
			ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
			ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)				
			
			deallocate(matU,matV,Singular,mat_tmp)
		end do
		n1 = OMP_get_wtime()
		do level_loc = 1,level_butterflyR
			level = levelm+1-level_loc
			
			if(level_loc/=level_butterflyR)then
				allocate(ButterflyP%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
			end if

			!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
			do index_ij_loc = 1, 2**level_butterflyR
				index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
				call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,SVD_tolerance_factor,ButterflyP_old,ButterflyP)				
			enddo
			!$omp end parallel do	
			
			do index_ij_loc = 1, 2**level_butterflyR
				index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
				index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))

				rank = size(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)				
				! rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
				! rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
								
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
				memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
				if (level_loc==level_butterflyR) then
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
				endif					
			enddo					
			
			
			! do index_i_loc=1, 2**(level_butterflyR-level_loc)
				! do index_j_loc=1,2**level_loc 
					! index_i = index_i_loc
					! index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
					! group_n=blocks%col_group   ! Note: row_group and col_group interchanged here    
					! group_n=group_n*2**(level_butterfly-level+1)-1+index_j
					! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
								
					! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
					! if(size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,2)/=nn)then
						! write(*,*)'nn incorrect'
						! stop
					! end if
					! mm1=size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,1)
					! mm2=size(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix,1)
					! mm = mm1+mm2
					
					! !!!!!!!!!!!!!!!!!!
					
					
					! allocate(QQ(mm,nn))
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn
						! do i=1, mm1
							! QQ(i,j)=ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn
						! do i=1, mm2
							! QQ(i+mm1,j)=ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					
					
					
					! mn=min(mm,nn)
					! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
					! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			

					! blocks%rankmax = max(blocks%rankmax,rank)
					! blocks%rankmin = min(blocks%rankmin,rank)

					! allocate(mat_tmp(rank,nn))
					! !$omp parallel do default(shared) private(i,j,k,ctemp)
					! do j=1, nn
						! do i=1, rank
							! mat_tmp(i,j)=VV(i,j)*Singular(i)
						! enddo
					! enddo
					! !$omp end parallel do							
					


					! if(level_loc/=level_butterflyR)then
						! nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
						! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
						! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
						
						! ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
						! ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)						
					! else 
						! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, nn
								! blocks%ButterflyV(index_j)%matrix(i,j)=mat_tmp(j,i)
								! ! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do

					! end if		

					! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
					! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm2,rank))
					
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, mm1
						! do j=1, rank
							! blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(i,j)=UU(i,j)
							! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, mm2
						! do j=1, rank
							! blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
							! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
						! enddo
					! enddo
					! !$omp end parallel do						
					
					! deallocate(QQ,UU,VV,Singular,mat_tmp)
				  
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
					! if (level_loc==level_butterflyR) then
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
					! endif
				! enddo
			! enddo
			
			if(level_loc/=level_butterflyR)then
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_butterflyR-level_loc+1)
						do jj =1,2**(level_loc)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if
				allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
				do ii = 1,2**(level_butterflyR-level_loc)
					do jj =1,2**(level_loc+1)
						mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
						nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
						ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
						deallocate(ButterflyP%blocks(ii,jj)%matrix)
					end do
				end do
				deallocate(ButterflyP%blocks)
			else 
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level_butterflyR-level_loc+1)
						do jj =1,2**(level_loc)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if					
			end if				
			
		end do	
		n2 = OMP_get_wtime()
		time_tmp = time_tmp + n2 - n1

	enddo		

	if(ranktotL/=ranktotR)then
		write(*,*)ranktotL,ranktotR,'ACA is not unique for applying the same matrix twice'
		stop
	end if
    
    Memory=memory_butterfly
    !write (*,*) memory_butterfly
    !pause

    return

end subroutine Butterfly_compress_N15_withoutBoundary_givenfullmat



subroutine Butterfly_compress_N15(blocks,Memory)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    integer blocks_idx, i, j, level_butterfly,level_butterflyL,level_butterflyR, num_blocks, k, attempt
    integer group_m, group_n, mm, nn,mn, index_i, index_j,index_ij_loc,index_i_m, index_j_m,index_i_loc, index_j_loc, ii, jj,nn1,nn2,mm1,mm2,idxs_m,idxs_n
    integer level,levelm,level_loc, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    integer rank, rankmax, butterflyB_inuse, rank1, rank2,rmax
    real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real*8 Memory,n1,n2
    complex(kind=8) ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	complex(kind=8), allocatable :: tau_Q(:)
	real*8,allocatable :: Singular(:)
	integer flag,tmpi,tmpj
	real*8:: a,b
	real*8:: maxvalue(1:20)
	real*8:: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row
	type(matrixblock)::blocks
	integer cnt_tmp,rankFar,rankNear
	
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	ForwardSymmetricFlag = 1	
	maxvalue = 0
	minvalue = 10000
    Memory=0

	
	! write(*,*)blocks%row_group,blocks%col_group,'In Butterfly_compress_N15'
	
    level_blocks=blocks%level
    !level_butterfly=Maxlevel-level_blocks
    ! level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2

!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
	blocks%rankmax = -100000
	blocks%rankmin = 100000
	
    ! blocks%level_butterfly=level_butterfly

	level_butterfly = blocks%level_butterfly
		
    num_blocks=2**level_butterfly

    allocate(blocks%ButterflyU(num_blocks))
    allocate(blocks%ButterflyV(num_blocks))
    if (level_butterfly/=0) then
        allocate(blocks%ButterflyKerl(level_butterfly))
    endif
    
    memory_butterfly=0.
    
	if(level_butterfly==0)then
		allocate (rankmax_for_butterfly(0:level_butterfly))
		rankmax_for_butterfly=0
		allocate (rankmin_for_butterfly(0:level_butterfly))
		rankmin_for_butterfly=10000
		
		group_m=blocks%row_group  ! Note: row_group and col_group interchanged here  
		group_n=blocks%col_group

		mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
		nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
		
		
		! !!!!! SVD
		! mn=min(mm,nn)
		! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))		
		! allocate(QQ(mm,nn))
		! do ii=1,mm
			! do jj =1,nn
				! edge_m = basis_group(group_m)%head + ii - 1
				! edge_n = basis_group(group_n)%head + jj - 1
				! call element_Zmn(edge_m,edge_n,ctemp)
				! QQ(ii,jj) = ctemp
			! end do
		! end do
		! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			
		! deallocate(QQ)

		
		!!!!! ACA-SVD
		idxs_m = basis_group(group_m)%head
		idxs_n = basis_group(group_n)%head
		
		rmax = min(500,min(mm,nn))
		allocate(UU(mm,rmax))
		allocate(VV(rmax,nn))
		allocate(Singular(rmax))
		call ACA_CompressionForward(UU,VV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)		
		
		
		
		
		! rank = 7
		rankmax_for_butterfly(0)=max(rank,rankmax_for_butterfly(0))
		rankmin_for_butterfly(0)=min(rank,rankmin_for_butterfly(0))

		blocks%rankmax = max(blocks%rankmax,rank)
		blocks%rankmin = min(blocks%rankmin,rank)
		
		allocate (blocks%ButterflyV(1)%matrix(nn,rank))
		
		!$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, nn
				blocks%ButterflyV(1)%matrix(i,j)=VV(j,i)
				! blocks%ButterflyV(1)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		!$omp end parallel do	
		
		
		allocate (blocks%ButterflyU(1)%matrix(mm,rank))
		
		!$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, mm
				blocks%ButterflyU(1)%matrix(i,j)=UU(i,j)*Singular(j)
				! blocks%ButterflyU(1)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		!$omp end parallel do							
		deallocate (UU,VV,Singular)			
			
		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(1)%matrix)/1024.0d3
		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(1)%matrix)/1024.0d3
		
	else if(level_butterfly==1)then
		allocate (rankmax_for_butterfly(0:level_butterfly))
		rankmax_for_butterfly=0
		allocate (rankmin_for_butterfly(0:level_butterfly))
		rankmin_for_butterfly=10000
		
		do level=0, level_butterfly
			index_ij=0
			if (level>0) then
				blocks%ButterflyKerl(level)%num_row=2**level
				blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
				allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
			endif
			if(level/=level_butterfly)then
				allocate(ButterflyP%blocks(2**(level+1),2**(level_butterfly-level)))
			end if
			
			do index_i=1, 2**level
				do index_j=1, 2**(level_butterfly-level)
					
					if(level==0)then
						group_m=blocks%row_group  ! Note: row_group and col_group interchanged here  
						group_n=blocks%col_group
						group_n=group_n*2**level_butterfly-1+index_j

						mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
						nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
						allocate(QQ(mm,nn))
						do ii=1,mm
							do jj =1,nn
								edge_m = basis_group(group_m)%head + ii - 1
								edge_n = basis_group(group_n)%head + jj - 1
								call element_Zmn(edge_m,edge_n,ctemp)
								QQ(ii,jj) = ctemp
							end do
						end do

						mn=min(mm,nn)
						allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			

						! rank = 7
						! write(*,*)rank
						rankmax_for_butterfly(0)=max(rank,rankmax_for_butterfly(0))
						rankmin_for_butterfly(0)=min(rank,rankmin_for_butterfly(0))
						blocks%rankmax = max(blocks%rankmax,rank)
						blocks%rankmin = min(blocks%rankmin,rank)

						allocate(mat_tmp(mm,rank))
						!$omp parallel do default(shared) private(i,j,k,ctemp)
						do j=1, rank
							do i=1, mm
								mat_tmp(i,j)=UU(i,j)*Singular(j)
							enddo
						enddo
						!$omp end parallel do
						
						mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
						allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
						allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank))
						
						ButterflyP%blocks(2*index_i-1,index_j)%matrix=mat_tmp(1:mm1,1:rank)
						ButterflyP%blocks(2*index_i,index_j)%matrix=mat_tmp(1+mm1:mm,1:rank)
					
						allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
						
						!$omp parallel do default(shared) private(i,j)
						do j=1, rank
							do i=1, nn
								blocks%ButterflyV(index_j)%matrix(i,j)=VV(j,i)
								! blocks%ButterflyV(index_j)%matrix(i,j)=random_complex_number()
							enddo
						enddo
						!$omp end parallel do					
						deallocate (QQ,UU,VV,Singular,mat_tmp)	

					else

						group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
						group_n=blocks%col_group    
						group_m=group_m*2**level-1+index_i
						group_n=group_n*2**(level_butterfly-level)-1+index_j
					
						mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
									
						! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
						if(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)/=mm)then
							write(*,*)'mm incorrect'
							stop
						end if
						nn1=size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,2)
						nn2=size(ButterflyP_old%blocks(index_i,2*index_j)%matrix,2)
						nn = nn1+nn2
						

						allocate(QQ(mm,nn))
						!$omp parallel do default(shared) private(i,j)
						do j=1, nn1
							do i=1, mm
								QQ(i,j)=ButterflyP_old%blocks(index_i,2*index_j-1)%matrix(i,j)
							enddo
						enddo
						!$omp end parallel do						
						!$omp parallel do default(shared) private(i,j)
						do j=1, nn2
							do i=1, mm
								QQ(i,j+nn1)=ButterflyP_old%blocks(index_i,2*index_j)%matrix(i,j)
							enddo
						enddo
						!$omp end parallel do						
						
						
						
						mn=min(mm,nn)
						allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			

						! rank = 7
						! write(*,*)rank
						rankmax_for_butterfly(level)=max(rank,rankmax_for_butterfly(level))
						rankmin_for_butterfly(level)=min(rank,rankmin_for_butterfly(level))
						blocks%rankmax = max(blocks%rankmax,rank)
						blocks%rankmin = min(blocks%rankmin,rank)

						allocate(mat_tmp(mm,rank))
						!$omp parallel do default(shared) private(i,j,k,ctemp)
						do j=1, rank
							do i=1, mm
								mat_tmp(i,j)=UU(i,j)*Singular(j)
							enddo
						enddo
						!$omp end parallel do							
						


						if(level/=level_butterfly)then
							mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
							allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
							allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank))
							
							ButterflyP%blocks(2*index_i-1,index_j)%matrix=mat_tmp(1:mm1,1:rank)
							ButterflyP%blocks(2*index_i,index_j)%matrix=mat_tmp(1+mm1:mm,1:rank)						
						else 
							allocate (blocks%ButterflyU(index_i)%matrix(mm,rank))
							!$omp parallel do default(shared) private(i,j)
							do j=1, rank
								do i=1, mm
									blocks%ButterflyU(index_i)%matrix(i,j)=mat_tmp(i,j)
									! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
								enddo
							enddo
							!$omp end parallel do

						end if		

						allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1))
						allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2))
						
						!$omp parallel do default(shared) private(i,j)
						do i=1, rank
							do j=1, nn1
								blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=VV(i,j)
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
							enddo
						enddo
						!$omp end parallel do						
						!$omp parallel do default(shared) private(i,j)
						do i=1, rank
							do j=1, nn2
								blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
							enddo
						enddo
						!$omp end parallel do						
						
						deallocate(QQ,UU,VV,Singular,mat_tmp)

					end if
					
					index_ij=index_ij+1
					if (level==0) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_ij)%matrix)/1024.0d3
					else                    
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
					endif
					if (level==level_butterfly) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_ij)%matrix)/1024.0d3
					endif
				enddo
			enddo
			
			if(level/=level_butterfly)then
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level)
						do jj =1,2**(level_butterfly-level+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if
				allocate(ButterflyP_old%blocks(2**(level+1),2**(level_butterfly-level)))
				do ii = 1,2**(level+1)
					do jj =1,2**(level_butterfly-level)
						mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
						nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
						ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
						deallocate(ButterflyP%blocks(ii,jj)%matrix)
					end do
				end do
				deallocate(ButterflyP%blocks)
			else 
				if(allocated(ButterflyP_old%blocks))then
					do ii = 1,2**(level)
						do jj =1,2**(level_butterfly-level+1)
							deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						end do
					end do	
					deallocate(ButterflyP_old%blocks)
				end if					
			end if
		enddo		
	else 
		
		do level=1, level_butterfly
			blocks%ButterflyKerl(level)%num_row=2**level
			blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
			allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
		end do
		
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		level_butterflyL = level_butterfly-levelm
		level_butterflyR = level_butterfly-level_butterflyL
		
		allocate (rankmax_for_butterfly(-level_butterflyR:level_butterflyL))
		rankmax_for_butterfly=0		
		allocate (rankmin_for_butterfly(-level_butterflyR:level_butterflyL))
		rankmin_for_butterfly=100000	
		
		allocate(blocks%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))
		
		
		! construct the middle level and the left half
		
		do index_i_m=1, 2**levelm
			
			write(*,*)'ddc1'
			
			level_loc = 0
			index_i_loc = 1
			allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
			
			do index_j_m=1, 2**(level_butterfly-levelm)	
				
				write(*,*)index_i_m,2**levelm,index_j_m,2**(level_butterfly-levelm),'c1'
				
				index_j_loc = index_j_m
				group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
				group_n=blocks%col_group    
				group_m=group_m*2**levelm-1+index_i_m
				group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
			
				mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
				nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
				idxs_m = basis_group(group_m)%head
				idxs_n = basis_group(group_n)%head
				
				rmax = min(500,min(mm,nn))
				allocate(matU(mm,rmax))
				allocate(matV(rmax,nn))
				allocate(Singular(rmax))
				
				write(*,*)index_i_m,2**levelm,index_j_m,2**(level_butterfly-levelm),'cao'
				
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)					
				! rank = min(rank,37)
				
				write(*,*)index_i_m,2**levelm,index_j_m,2**(level_butterfly-levelm),'cao1',Singular(1:rank)
				
				! if(rank==53)then
					! write(*,*)group_m,group_n,blocks%row_group,blocks%col_group,blocks%level_butterfly,blocks%level
					! stop
				! end if
				
				
				! if(near_or_far(group_m,group_n,2d0))write(*,*)near_or_far(group_m,group_n,2d0),rank
				
				! write(*,*)near_or_far(group_m,group_n,1.5d0),rank
				! if(.not. near_or_far(group_m,group_n,1.5d0))then
					! cnt_tmp = cnt_tmp+1
					! rankNear = max(rankNear,rank)
				! else 
					! rankFar= max(rankFar,rank)
				! end if
				
				
				! if(.not. near_or_far(group_m,group_n,2d0))then
					! rank = 1
					! Singular(1:rank) = 0
					! ! if(blocks==342)write(111,*)Singular(1:rank)
					! rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
					! rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
					
					! allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
					! blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
					! do ii=1,rank
						! blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 0d0
					! end do					
				! else 
					! if(blocks==342)write(111,*)Singular(1:rank)
					rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
					rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
					blocks%rankmax = max(blocks%rankmax,rank)
					blocks%rankmin = min(blocks%rankmin,rank)
					
					
					allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
					do ii=1,rank
						blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 1d0/Singular(ii)
					end do
				! end if
				
				write(*,*)'aaa'
				
				allocate(mat_tmp(mm,rank))
				!$omp parallel do default(shared) private(i,j,k,ctemp)
				do j=1, rank
					do i=1, mm
						mat_tmp(i,j)=matU(i,j)*Singular(j)
					enddo
				enddo
				!$omp end parallel do
				
				write(*,*)'aaass'
				
				mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
				allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
				allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
				
				ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
				ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)				
				
				deallocate(matU,matV,Singular,mat_tmp)
				
				write(*,*)'addaa'
				
			end do
			
			n1 = OMP_get_wtime()
			
			do level_loc = 1,level_butterflyL
				level = level_loc+levelm
				
				if(level_loc/=level_butterflyL)then
					allocate(ButterflyP%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
				end if

				! do index_i_loc=1, 2**level_loc
					! do index_j_loc=1, 2**(level_butterflyL-level_loc)
					
				!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
				do index_ij_loc = 1, 2**level_butterflyL
					index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				
					call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,SVD_tolerance_forward,ButterflyP_old,ButterflyP)				
				enddo
				!$omp end parallel do
				
				do index_ij_loc = 1, 2**level_butterflyL
					index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
					index_j = index_j_loc
					index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
					
					rank = size(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
		
					rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
					rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
					
					blocks%rankmax = max(blocks%rankmax,rank)
					blocks%rankmin = min(blocks%rankmin,rank)
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
					if (level_loc==level_butterflyL) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
					endif
				enddo				
				
				
				if(level_loc/=level_butterflyL)then
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_loc)
							do jj =1,2**(level_butterflyL-level_loc+1)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if
					allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
					do ii = 1,2**(level_loc+1)
						do jj =1,2**(level_butterflyL-level_loc)
							mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
							nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
							ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
							deallocate(ButterflyP%blocks(ii,jj)%matrix)
						end do
					end do
					deallocate(ButterflyP%blocks)
				else 
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_loc)
							do jj =1,2**(level_butterflyL-level_loc+1)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if					
				end if				
				
			end do	
			n2 = OMP_get_wtime()
			time_tmp = time_tmp + n2 - n1

		enddo

		
		! write(*,*)'stat:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	


		! construct the the right half
		do index_j_m=1, 2**(level_butterfly-levelm)	
			level_loc = 0
			index_j_loc = 1
			allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
			do index_i_m=1, 2**levelm			
				index_i_loc = index_i_m
				group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
				group_n=blocks%col_group    
				group_m=group_m*2**levelm-1+index_i_m
				group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
			
				mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
				nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
				idxs_m = basis_group(group_m)%head
				idxs_n = basis_group(group_n)%head
				
				rmax = min(500,min(mm,nn))
				allocate(matU(mm,rmax))
				allocate(matV(rmax,nn))
				allocate(Singular(rmax))
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)	
				! rank = min(rank,37)
				

				! if(.not. near_or_far(group_m,group_n,2d0))then
					! rank = 1
					! Singular(1:rank) = 0
					! ! if(blocks==342)write(111,*)Singular(1:rank)

				! end if
				
				
				
				allocate(mat_tmp(rank,nn))
				!$omp parallel do default(shared) private(i,j,k,ctemp)
				do j=1, nn
					do i=1, rank
						mat_tmp(i,j)=matV(i,j)*Singular(i)
					enddo
				enddo
				!$omp end parallel do
				
				nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
				
				ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
				ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)				
				
				deallocate(matU,matV,Singular,mat_tmp)
			end do
			
			n1 = OMP_get_wtime()
			do level_loc = 1,level_butterflyR
				level = levelm+1-level_loc
				
				if(level_loc/=level_butterflyR)then
					allocate(ButterflyP%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
				end if

				
				!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
				do index_ij_loc = 1, 2**level_butterflyR
					index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
					call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,SVD_tolerance_forward,ButterflyP_old,ButterflyP)				
				enddo
				!$omp end parallel do	
				
				do index_ij_loc = 1, 2**level_butterflyR
					index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
					
					index_i = index_i_loc
					index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
						
					rank = size(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)				
					rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
					rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
					blocks%rankmax = max(blocks%rankmax,rank)
					blocks%rankmin = min(blocks%rankmin,rank)
									
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
					if (level_loc==level_butterflyR) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
					endif					
				enddo				
				
				! do index_i_loc=1, 2**(level_butterflyR-level_loc)
					! do index_j_loc=1,2**level_loc 
					
						! index_i = index_i_loc
						! index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
						! group_n=blocks%col_group   ! Note: row_group and col_group interchanged here    
						! group_n=group_n*2**(level_butterfly-level+1)-1+index_j
						! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
									
						! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
						! if(size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,2)/=nn)then
							! write(*,*)'nn incorrect'
							! stop
						! end if
						! mm1=size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,1)
						! mm2=size(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix,1)
						! mm = mm1+mm2
						
						! !!!!!!!!!!!!!!!!!!
						
						
						! allocate(QQ(mm,nn))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn
							! do i=1, mm1
								! QQ(i,j)=ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn
							! do i=1, mm2
								! QQ(i+mm1,j)=ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						
						
						
						! mn=min(mm,nn)
						! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			
						! ! rank = min(rank,37)

						! ! rank = 7
						! ! write(*,*)rank
						! rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
						! rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
						! blocks%rankmax = max(blocks%rankmax,rank)
						! blocks%rankmin = min(blocks%rankmin,rank)

						! allocate(mat_tmp(rank,nn))
						! !$omp parallel do default(shared) private(i,j,k,ctemp)
						! do j=1, nn
							! do i=1, rank
								! mat_tmp(i,j)=VV(i,j)*Singular(i)
							! enddo
						! enddo
						! !$omp end parallel do							
						


						! if(level_loc/=level_butterflyR)then
							! nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
							! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
							! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
							
							! ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
							! ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)						
						! else 
							! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
							! !$omp parallel do default(shared) private(i,j)
							! do j=1, rank
								! do i=1, nn
									! blocks%ButterflyV(index_j)%matrix(i,j)=mat_tmp(j,i)
									! ! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
								! enddo
							! enddo
							! !$omp end parallel do

						! end if		

						! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
						! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm2,rank))
						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, mm1
							! do j=1, rank
								! blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(i,j)=UU(i,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, mm2
							! do j=1, rank
								! blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						
						! deallocate(QQ,UU,VV,Singular,mat_tmp)
					  
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
						! if (level_loc==level_butterflyR) then
							! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
						! endif
					! enddo
				! enddo
				
				if(level_loc/=level_butterflyR)then
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_butterflyR-level_loc+1)
							do jj =1,2**(level_loc)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if
					allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
					do ii = 1,2**(level_butterflyR-level_loc)
						do jj =1,2**(level_loc+1)
							mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
							nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
							ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
							deallocate(ButterflyP%blocks(ii,jj)%matrix)
						end do
					end do
					deallocate(ButterflyP%blocks)
				else 
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_butterflyR-level_loc+1)
							do jj =1,2**(level_loc)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if					
				end if				
				
			end do	

			n2 = OMP_get_wtime()
			time_tmp = time_tmp + n2 - n1
		enddo		

			


	end if	
	 
	! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
	rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),rankmax_of_level(level_blocks))
	
	! write(*,*)rankmax_of_level,'nitaima',rankmax_for_butterfly
	! stop
	
	! if(blocks==342)then
	! write(*,*)rankmax_for_butterfly
	! write(*,*)rankmin_for_butterfly
	! end if
	
	! write(*,*)'max value: ',maxvalue(1:9)
	! write(*,*)'min value: ',minvalue(1:9)
	
	
    deallocate (rankmax_for_butterfly)
    deallocate (rankmin_for_butterfly)
    
    Memory=memory_butterfly
    !write (*,*) memory_butterfly
    !pause

    return

end subroutine Butterfly_compress_N15







subroutine Butterfly_compress_N15_givenfullmat(blocks,idx_m_ref,idx_n_ref)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   use misc
   implicit none

    integer index_ij_loc,blocks_idx, i, j, level_butterfly,level_butterflyL,level_butterflyR, num_blocks, k, attempt
    integer group_m, group_n, mm, nn,mn, index_i, index_j,index_i_m, index_j_m,index_i_loc, index_j_loc, ii, jj,nn1,nn2,mm1,mm2,idxs_m,idxs_n
    integer level,levelm,level_loc, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    integer rank, rankmax, butterflyB_inuse, rank1, rank2,rmax,ranktotL,ranktotR
    real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real*8 Memory
    complex(kind=8) ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	complex(kind=8), allocatable :: tau_Q(:)
	real*8,allocatable :: Singular(:)
	integer flag,tmpi,tmpj
	real*8:: a,b,n1,n2
	real*8:: maxvalue(1:20)
	real*8:: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row
	type(matrixblock)::blocks
	integer cnt_tmp,rankFar,rankNear,idx_m_ref,idx_n_ref
	
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	ForwardSymmetricFlag = 1	
	maxvalue = 0
	minvalue = 10000
    Memory=0

    level_blocks=blocks%level
    !level_butterfly=Maxlevel-level_blocks
    level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2
	
!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
	blocks%rankmax = -100000
	blocks%rankmin = 100000
	
    blocks%level_butterfly=level_butterfly


	! write(*,*)blocks%level,level_butterfly,'nimaib'
	
	
	
    num_blocks=2**level_butterfly

    allocate(blocks%ButterflyU(num_blocks))
    allocate(blocks%ButterflyV(num_blocks))
    if (level_butterfly/=0) then
        allocate(blocks%ButterflyKerl(level_butterfly))
    endif
    
    memory_butterfly=0.
    
	if(level_butterfly==0)then
		
		group_m=blocks%row_group  ! Note: row_group and col_group interchanged here  
		group_n=blocks%col_group

		mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
		nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
		
		!!!!! ACA-SVD
		idxs_m = basis_group(group_m)%head
		idxs_n = basis_group(group_n)%head
		
		rmax = min(500,min(mm,nn))
		allocate(UU(mm,rmax))
		allocate(VV(rmax,nn))
		allocate(Singular(rmax))
		! call ACA_CompressionForward(UU,VV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_forward*0.1,SVD_tolerance_forward)		
		call ACA_CompressionForward_givenfullmat(UU,VV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_factor*0.1,SVD_tolerance_factor,idx_m_ref,idx_n_ref)	
		
		
		blocks%rankmax = max(blocks%rankmax,rank)
		blocks%rankmin = min(blocks%rankmin,rank)
		
		allocate (blocks%ButterflyV(1)%matrix(nn,rank))
		
		!$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, nn
				blocks%ButterflyV(1)%matrix(i,j)=VV(j,i)
				! blocks%ButterflyV(1)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		!$omp end parallel do	
		
		
		allocate (blocks%ButterflyU(1)%matrix(mm,rank))
		
		!$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, mm
				blocks%ButterflyU(1)%matrix(i,j)=UU(i,j)*Singular(j)
				! blocks%ButterflyU(1)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		!$omp end parallel do							
		deallocate (UU,VV,Singular)			
			
		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(1)%matrix)/1024.0d3
		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(1)%matrix)/1024.0d3
				
	else if(level_butterfly==1)then
		write(*,*)'1 level butterfly removed'	
	else 
		
		do level=1, level_butterfly
			blocks%ButterflyKerl(level)%num_row=2**level
			blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
			allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
		end do
		
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		level_butterflyL = level_butterfly-levelm
		level_butterflyR = level_butterfly-level_butterflyL

		
		allocate(blocks%ButterflyMiddle(2**levelm,2**(level_butterfly-levelm)))
		
		ranktotL = 0
		! construct the middle level and the left half
		
		do index_i_m=1, 2**levelm
			level_loc = 0
			index_i_loc = 1
			allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
			
			do index_j_m=1, 2**(level_butterfly-levelm)	
				index_j_loc = index_j_m
				group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
				group_n=blocks%col_group    
				group_m=group_m*2**levelm-1+index_i_m
				group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
			
				mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
				nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
				idxs_m = basis_group(group_m)%head
				idxs_n = basis_group(group_n)%head
				
				rmax = min(500,min(mm,nn))
				allocate(matU(mm,rmax))
				allocate(matV(rmax,nn))
				allocate(Singular(rmax))
				call ACA_CompressionForward_givenfullmat(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_factor*0.1,SVD_tolerance_factor,idx_m_ref,idx_n_ref)					
				ranktotL = ranktotL + rank
				
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
				
				allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank))
				blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
				do ii=1,rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 1d0/Singular(ii)
				end do

				

				
				allocate(mat_tmp(mm,rank))
				!$omp parallel do default(shared) private(i,j,k,ctemp)
				do j=1, rank
					do i=1, mm
						mat_tmp(i,j)=matU(i,j)*Singular(j)
					enddo
				enddo
				!$omp end parallel do
				
				mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
				allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
				allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
				
				ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
				ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)				
				
				deallocate(matU,matV,Singular,mat_tmp)
			end do
			
			n1 = OMP_get_wtime()
			do level_loc = 1,level_butterflyL
				level = level_loc+levelm
				
				if(level_loc/=level_butterflyL)then
					allocate(ButterflyP%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
				end if

				!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)	
				do index_ij_loc = 1, 2**level_butterflyL
					index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
					call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,SVD_tolerance_factor,ButterflyP_old,ButterflyP)				
				enddo
				!$omp end parallel do	
				
				do index_ij_loc = 1, 2**level_butterflyL
					index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
					index_j = index_j_loc
					index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
					
					rank = size(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix,1)
		
					! rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
					! rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))
					
					blocks%rankmax = max(blocks%rankmax,rank)
					blocks%rankmin = min(blocks%rankmin,rank)
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
					if (level_loc==level_butterflyL) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
					endif
				enddo					
				
				
				! do index_i_loc=1, 2**level_loc
					! do index_j_loc=1, 2**(level_butterflyL-level_loc)
						! index_j = index_j_loc
						! index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
						! group_m=blocks%row_group   ! Note: row_group and col_group interchanged here    
						! group_m=group_m*2**level-1+index_i
						! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
									
						! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
						! if(size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,1)/=mm)then
							! write(*,*)'mm incorrect',size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,1),mm,basis_group(group_m)%tail,basis_group(group_m)%head,group_m
							! stop
						! end if
						! nn1=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,2)
						! nn2=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix,2)
						! nn = nn1+nn2
						

						! allocate(QQ(mm,nn))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn1
							! do i=1, mm
								! QQ(i,j)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn2
							! do i=1, mm
								! QQ(i,j+nn1)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						
						
						
						! mn=min(mm,nn)
						! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_factor,rank)			

						
						! blocks%rankmax = max(blocks%rankmax,rank)
						! blocks%rankmin = min(blocks%rankmin,rank)

						! allocate(mat_tmp(mm,rank))
						! !$omp parallel do default(shared) private(i,j,k,ctemp)
						! do j=1, rank
							! do i=1, mm
								! mat_tmp(i,j)=UU(i,j)*Singular(j)
							! enddo
						! enddo
						! !$omp end parallel do							
						


						! if(level_loc/=level_butterflyL)then
							! mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
							! allocate(ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
							! allocate(ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
							
							! ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
							! ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)						
						! else 
							! allocate (blocks%ButterflyU(index_i)%matrix(mm,rank))
							! !$omp parallel do default(shared) private(i,j)
							! do j=1, rank
								! do i=1, mm
									! blocks%ButterflyU(index_i)%matrix(i,j)=mat_tmp(i,j)
									! ! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
								! enddo
							! enddo
							! !$omp end parallel do

						! end if		

						! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1))
						! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2))
						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, rank
							! do j=1, nn1
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=VV(i,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, rank
							! do j=1, nn2
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						
						! deallocate(QQ,UU,VV,Singular,mat_tmp)
					  
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
						! if (level_loc==level_butterflyL) then
							! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
						! endif
					! enddo
				! enddo
				
				if(level_loc/=level_butterflyL)then
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_loc)
							do jj =1,2**(level_butterflyL-level_loc+1)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if
					allocate(ButterflyP_old%blocks(2**(level_loc+1),2**(level_butterflyL-level_loc)))
					do ii = 1,2**(level_loc+1)
						do jj =1,2**(level_butterflyL-level_loc)
							mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
							nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
							ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
							deallocate(ButterflyP%blocks(ii,jj)%matrix)
						end do
					end do
					deallocate(ButterflyP%blocks)
				else 
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_loc)
							do jj =1,2**(level_butterflyL-level_loc+1)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if					
				end if				
				
			end do	
			n2 = OMP_get_wtime()
			time_tmp = time_tmp + n2 - n1

		enddo

		
		! write(*,*)'stat:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	

		ranktotR = 0
		! construct the the right half
		do index_j_m=1, 2**(level_butterfly-levelm)	
			level_loc = 0
			index_j_loc = 1
			allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
			do index_i_m=1, 2**levelm			
				index_i_loc = index_i_m
				group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
				group_n=blocks%col_group    
				group_m=group_m*2**levelm-1+index_i_m
				group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
			
				mm=basis_group(group_m)%tail-basis_group(group_m)%head+1	
				nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
				idxs_m = basis_group(group_m)%head
				idxs_n = basis_group(group_n)%head
				
				rmax = min(500,min(mm,nn))
				allocate(matU(mm,rmax))
				allocate(matV(rmax,nn))
				allocate(Singular(rmax))
				call ACA_CompressionForward_givenfullmat(matU,matV,Singular,idxs_m,idxs_n,mm,nn,rmax,rank,SVD_tolerance_factor*0.1,SVD_tolerance_factor,idx_m_ref,idx_n_ref)	
				! rank = min(rank,37)
				ranktotR = ranktotR + rank
				
				blocks%rankmax = max(blocks%rankmax,rank)
				blocks%rankmin = min(blocks%rankmin,rank)
				
				allocate(mat_tmp(rank,nn))
				!$omp parallel do default(shared) private(i,j,k,ctemp)
				do j=1, nn
					do i=1, rank
						mat_tmp(i,j)=matV(i,j)*Singular(i)
					enddo
				enddo
				!$omp end parallel do
				
				nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
				
				ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
				ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)				
				
				deallocate(matU,matV,Singular,mat_tmp)
			end do
			
			n1 = OMP_get_wtime()
			do level_loc = 1,level_butterflyR
				level = levelm+1-level_loc
				
				if(level_loc/=level_butterflyR)then
					allocate(ButterflyP%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
				end if

				
				!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
				do index_ij_loc = 1, 2**level_butterflyR
					index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
					call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,SVD_tolerance_factor,ButterflyP_old,ButterflyP)				
				enddo
				!$omp end parallel do	
				
				do index_ij_loc = 1, 2**level_butterflyR
					index_j_loc = mod(index_ij_loc-1,2**level_loc)+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
					index_i = index_i_loc
					index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
					
					rank = size(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,2)				
					! rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
					! rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
					blocks%rankmax = max(blocks%rankmax,rank)
					blocks%rankmin = min(blocks%rankmin,rank)
									
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
					if (level_loc==level_butterflyR) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
					endif					
				enddo					
				
				! do index_i_loc=1, 2**(level_butterflyR-level_loc)
					! do index_j_loc=1,2**level_loc 
						! index_i = index_i_loc
						! index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
						! group_n=blocks%col_group   ! Note: row_group and col_group interchanged here    
						! group_n=group_n*2**(level_butterfly-level+1)-1+index_j
						! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
									
						! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
						! if(size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,2)/=nn)then
							! write(*,*)'nn incorrect'
							! stop
						! end if
						! mm1=size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,1)
						! mm2=size(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix,1)
						! mm = mm1+mm2
						
						! !!!!!!!!!!!!!!!!!!
						
						
						! allocate(QQ(mm,nn))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn
							! do i=1, mm1
								! QQ(i,j)=ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, nn
							! do i=1, mm2
								! QQ(i+mm1,j)=ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(i,j)
							! enddo
						! enddo
						! !$omp end parallel do						
						
						
						
						! mn=min(mm,nn)
						! allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_factor,rank)			
						! ! rank = min(rank,37)

						! blocks%rankmax = max(blocks%rankmax,rank)
						! blocks%rankmin = min(blocks%rankmin,rank)

						! allocate(mat_tmp(rank,nn))
						! !$omp parallel do default(shared) private(i,j,k,ctemp)
						! do j=1, nn
							! do i=1, rank
								! mat_tmp(i,j)=VV(i,j)*Singular(i)
							! enddo
						! enddo
						! !$omp end parallel do							
						


						! if(level_loc/=level_butterflyR)then
							! nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
							! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
							! allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
							
							! ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
							! ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)						
						! else 
							! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
							! !$omp parallel do default(shared) private(i,j)
							! do j=1, rank
								! do i=1, nn
									! blocks%ButterflyV(index_j)%matrix(i,j)=mat_tmp(j,i)
									! ! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
								! enddo
							! enddo
							! !$omp end parallel do

						! end if		

						! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
						! allocate (blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm2,rank))
						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, mm1
							! do j=1, rank
								! blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(i,j)=UU(i,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, mm2
							! do j=1, rank
								! blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						
						! deallocate(QQ,UU,VV,Singular,mat_tmp)
					  
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
						! if (level_loc==level_butterflyR) then
							! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
						! endif
					! enddo
				! enddo
				
				if(level_loc/=level_butterflyR)then
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_butterflyR-level_loc+1)
							do jj =1,2**(level_loc)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if
					allocate(ButterflyP_old%blocks(2**(level_butterflyR-level_loc),2**(level_loc+1)))
					do ii = 1,2**(level_butterflyR-level_loc)
						do jj =1,2**(level_loc+1)
							mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
							nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
							ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
							deallocate(ButterflyP%blocks(ii,jj)%matrix)
						end do
					end do
					deallocate(ButterflyP%blocks)
				else 
					if(allocated(ButterflyP_old%blocks))then
						do ii = 1,2**(level_butterflyR-level_loc+1)
							do jj =1,2**(level_loc)
								deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
							end do
						end do	
						deallocate(ButterflyP_old%blocks)
					end if					
				end if				
				
			end do	

			n2 = OMP_get_wtime()
			time_tmp = time_tmp + n2 - n1
		enddo		

		if(ranktotL/=ranktotR)then
			write(*,*)ranktotL,ranktotR,'ACA is not unique for applying the same matrix twice in Butterfly_compress_N15_givenfullmat'
			stop
		end if			


	end if	
	 

	
    Memory=memory_butterfly
    !write (*,*) memory_butterfly
    !pause

    return

end subroutine Butterfly_compress_N15_givenfullmat





! subroutine Butterfly_compress_ID(blocks,Memory)

   ! use MODULE_FILE
   ! ! use lapack95
   ! ! use blas95
   ! use misc
   ! implicit none

    ! integer i, j, level_butterfly, num_blocks, k, attempt
    ! integer group_m, group_n, mm, nn,mn, index_i, index_j, ii, jj,mm1,nn1,nn2
    ! integer level, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    ! integer rank, rankmax, butterflyB_inuse, rank1, rank2
    ! real*8 rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    ! real*8 Memory
    ! complex(kind=8) ctemp
	! type(butterfly_Kerl)ButterflyP_old,ButterflyP
	! complex(kind=8), allocatable :: mu_mat_sub(:,:), mu_mat_sub_0(:,:),LL(:,:),RR(:,:),RR_permute(:,:),proj(:,:),identity(:,:),UU(:,:), VV(:,:),mat_tmp(:,:)
	! complex(kind=8), allocatable :: tau_Q(:)
	! real*8,allocatable :: Singular(:)
	! integer flag
	! integer,parameter:: len = 100000000
	! integer list(len),list_tmp(len)
	! real*8 work(len),approx(len)
	! integer,allocatable::ind_map(:)
	! real*8:: maxvalue(1:20)
	! type(matrixblock)::blocks
	
	! maxvalue = 0
    ! Memory=0

    ! level_blocks=blocks%level
    ! !level_butterfly=Maxlevel-level_blocks
    ! Preset_level_butterfly=Maxlevel_for_blocks
    ! level_butterfly=Preset_level_butterfly-level_blocks
	
	! blocks%rankmax = -100000
	! blocks%rankmin = 100000
	
! !     if (Maxlevel-level_blocks<8) then
! !         level_butterfly=Maxlevel-level_blocks
! !     endif
    ! blocks%level_butterfly=level_butterfly
    ! allocate (rankmax_for_butterfly(0:level_butterfly))
    ! rankmax_for_butterfly=0

    ! num_blocks=2**level_butterfly

    ! allocate(blocks%ButterflyU(num_blocks))
    ! allocate(blocks%ButterflyV(num_blocks))
    ! if (level_butterfly/=0) then
        ! allocate(blocks%ButterflyKerl(level_butterfly))
    ! endif
    
    ! memory_butterfly=0.
    ! do level=0, level_butterfly
        ! index_ij=0
        ! if (level>0) then
            ! blocks%ButterflyKerl(level)%num_row=2**level
            ! blocks%ButterflyKerl(level)%num_col=2**(level_butterfly-level+1)
            ! allocate (blocks%ButterflyKerl(level)%blocks(2**level,2**(level_butterfly-level+1)))
        ! endif
		! if(level/=level_butterfly)then
			! allocate(ButterflyP%blocks(2**(level+1),2**(level_butterfly-level)))
		! end if
		
		
		
		! ! ! if(level==0)then
			! ! ! group_n=blocks%col_group
			! ! ! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
			! ! ! allocate(ind_map(nn))
			! ! ! call rperm(nn,ind_map)
			! ! ! ! write(*,*)nn,'ha',ind_map
		! ! ! end if
		
		
		
        ! do index_i=1, 2**level
            ! do index_j=1, 2**(level_butterfly-level)
                
				! if(level==0)then
					! group_m=blocks%row_group  ! Note: row_group and col_group interchanged here  
					! group_n=blocks%col_group
					! group_n=group_n*2**level_butterfly-1+index_j

					! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
					! nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
					
					! allocate(mu_mat_sub(mm,nn))
					! allocate(mu_mat_sub_0(mm,nn))
					! do ii=1,mm
						! do jj =1,nn
							! edge_m = basis_group(group_m)%head + ii - 1
							! ! edge_n = basis_group(blocks%col_group)%head  +  ind_map(basis_group(group_n)%head - basis_group(blocks%col_group)%head + jj) - 1
							! edge_n = basis_group(group_n)%head + jj - 1 
							! call element_Zmn(edge_m,edge_n,ctemp)
							! mu_mat_sub(ii,jj) = ctemp
						! end do
					! end do
					
					! mu_mat_sub_0 = mu_mat_sub
				
					! call idzp_id(SVD_tolerance_forward,mm,nn,mu_mat_sub,rank,list,work)

					! ! rank = 7
					! ! call idzr_id(mm,nn,mu_mat_sub,rank,list,work)
					! ! write(*,*)level,list(1:rank)
					
					! call assert(rank<nn,'rank==nn is not handled properly')
					
					! if (rank>rankmax_for_butterfly(0)) then
						! rankmax_for_butterfly(0)=rank
					! endif
					! blocks%rankmax = max(blocks%rankmax,rank)
					! blocks%rankmin = min(blocks%rankmin,rank)

						
					! ! write(*,*)level,rank
						
					! allocate(LL(mm,rank))
					! call idz_copycols(mm,nn,mu_mat_sub_0,rank,list,LL)
					! allocate(proj(rank,nn-rank))
					! call DirectCopy(mu_mat_sub,proj,rank,nn-rank)
					! allocate(identity(rank,rank))
					! identity = 0
					! do ii=1,rank
						! identity(ii,ii)=1
					! end do
					! allocate(RR_permute(rank,nn))
					! RR_permute(1:rank,1:rank) = identity
					! RR_permute(1:rank,1+rank:nn) = proj
					! deallocate(proj,identity,mu_mat_sub_0,mu_mat_sub)
					! allocate(RR(rank,nn))
					! do ii=1,nn
						! RR(:,list(ii)) = RR_permute(:,ii)
					! end do
					! deallocate(RR_permute)
					
					
					! if(level/=level_butterfly)then
						
						! mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
						! allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
						! allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank))
						
						! ButterflyP%blocks(2*index_i-1,index_j)%matrix=LL(1:mm1,1:rank)
						! ButterflyP%blocks(2*index_i,index_j)%matrix=LL(1+mm1:mm,1:rank)
                    
						! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, nn
								! blocks%ButterflyV(index_j)%matrix(i,j)=RR(j,i)
							! enddo
						! enddo
						! !$omp end parallel do	
						
						! allocate(blocks%ButterflyV(index_j)%list(nn))
						! blocks%ButterflyV(index_j)%list = 0
						! do ii = 1,rank
							! blocks%ButterflyV(index_j)%list(list(ii)) = ii
						! end do
						

						
						! ! do ii = 1,nn
						! ! do jj = 1,rank
							! ! maxvalue(1) = max(maxvalue(1),abs(blocks%ButterflyV(index_j)%matrix(ii,jj)))	
						! ! end do
						! ! end do
						

						! do jj = 1,rank
							! maxvalue(1) = max(maxvalue(1),sqrt(norm_vector(blocks%ButterflyV(index_j)%matrix(:,jj),nn)))	
						! end do
					
						
					! else
						
						! allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, nn
								! blocks%ButterflyV(index_j)%matrix(i,j)=RR(j,i)
							! enddo
						! enddo
						! !$omp end parallel do	
						
						
						! allocate (blocks%ButterflyU(index_j)%matrix(mm,rank))
						
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, mm
								! blocks%ButterflyU(index_j)%matrix(i,j)=LL(i,j)
							! enddo
						! enddo
						! !$omp end parallel do										
					! end if
					! deallocate(LL,RR)
					
				! else
					! group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
					! group_n=blocks%col_group    
					! group_m=group_m*2**level-1+index_i
					! group_n=group_n*2**(level_butterfly-level)-1+index_j
				
					! mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
								
					! ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
					! if(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)/=mm)then
						! write(*,*)'mm incorrect'
						! stop
					! end if
					! nn1=size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,2)
					! nn2=size(ButterflyP_old%blocks(index_i,2*index_j)%matrix,2)
					! nn = nn1+nn2
					
					! allocate(mu_mat_sub(mm,nn))
					! allocate(mu_mat_sub_0(mm,nn))
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn1
						! do i=1, mm
							! mu_mat_sub(i,j)=ButterflyP_old%blocks(index_i,2*index_j-1)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do j=1, nn2
						! do i=1, mm
							! mu_mat_sub(i,j+nn1)=ButterflyP_old%blocks(index_i,2*index_j)%matrix(i,j)
						! enddo
					! enddo
					! !$omp end parallel do	
					
					! mu_mat_sub_0 = mu_mat_sub
					! call idzp_id(SVD_tolerance_forward,mm,nn,mu_mat_sub,rank,list,work)
					
					! ! rank = 7
					! ! call idzr_id(mm,nn,mu_mat_sub,rank,list,work)
					! ! write(*,*)level,list(1:rank)
					! call assert(rank<nn,'rank==nn is not handled properly')
					
					! ! ! write(*,*)level,rank
					
					! if (rank>rankmax_for_butterfly(level)) then
						! rankmax_for_butterfly(level)=rank
					! endif
					! blocks%rankmax = max(blocks%rankmax,rank)
					! blocks%rankmin = min(blocks%rankmin,rank)
					
					! ! write(*,*)level,rank
					
					! allocate(LL(mm,rank))
					! call idz_copycols(mm,nn,mu_mat_sub_0,rank,list,LL)
					! allocate(proj(rank,nn-rank))
					! call DirectCopy(mu_mat_sub,proj,rank,nn-rank)
					! allocate(identity(rank,rank))
					! identity = 0
					! do ii=1,rank
						! identity(ii,ii)=1
					! end do
					! allocate(RR_permute(rank,nn))
					! RR_permute(1:rank,1:rank) = identity
					! RR_permute(1:rank,1+rank:nn) = proj
					! deallocate(proj,identity,mu_mat_sub_0,mu_mat_sub)
					! allocate(RR(rank,nn))
					! do ii=1,nn
						! RR(:,list(ii)) = RR_permute(:,ii)
					! end do
					! deallocate(RR_permute)
					

					! if(level/=level_butterfly)then
						! mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
						! allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
						! allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank))
						
						! ButterflyP%blocks(2*index_i-1,index_j)%matrix=LL(1:mm1,1:rank)
						! ButterflyP%blocks(2*index_i,index_j)%matrix=LL(1+mm1:mm,1:rank)						
					! else 
						! allocate (blocks%ButterflyU(index_i)%matrix(mm,rank))
						! !$omp parallel do default(shared) private(i,j)
						! do j=1, rank
							! do i=1, mm
								! blocks%ButterflyU(index_i)%matrix(i,j)=LL(i,j)
							! enddo
						! enddo
						! !$omp end parallel do		

						! do jj =1,rank							
						! maxvalue(level_butterfly+2) = max(maxvalue(level_butterfly+2),sqrt(norm_vector(blocks%ButterflyU(index_i)%matrix(:,jj),nn)))							
						! end do	
						
					! end if		

					! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1))
					! allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2))
					
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, rank
						! do j=1, nn1
							! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=RR(i,j)
						! enddo
					! enddo
					! !$omp end parallel do						
					! !$omp parallel do default(shared) private(i,j)
					! do i=1, rank
						! do j=1, nn2
							! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=RR(i,j+nn1)
						! enddo
					! enddo
					! !$omp end parallel do

					! list_tmp(1:nn) = 0
					! do ii = 1,rank
						! list_tmp(list(ii)) = ii
					! end do
					
					! ! write(*,*)level,'level',list_tmp(1:nn)
					
					! allocate(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%list(nn1))
					! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%list = list_tmp(1:nn1)
					! allocate(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%list(nn2))
					! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%list = list_tmp(nn1+1:nn)
					

					! ! do ii = 1,rank
					! ! do jj = 1,nn1
						! ! maxvalue(level+1) = max(maxvalue(level+1),abs(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(ii,jj)))	
					! ! end do
					! ! end do
					! ! do ii = 1,rank
					! ! do jj = 1,nn2
						! ! maxvalue(level+1) = max(maxvalue(level+1),abs(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(ii,jj)))	
					! ! end do
					! ! end do
					
					! do jj = 1,rank
						! maxvalue(level+1) = max(maxvalue(level+1),sqrt(norm_vector(RR(jj,:),nn)))	
					! end do					
					
					! deallocate(LL,RR)
					
				! end if
									
				! index_ij=index_ij+1
                ! if (level==0) then
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_ij)%matrix)/1024.0d3
                ! else                    
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3					
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3					
				! endif
                ! if (level==level_butterfly) then
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_ij)%matrix)/1024.0d3
                ! endif
            ! enddo
        ! enddo
		
		! if(level/=level_butterfly)then
			! if(allocated(ButterflyP_old%blocks))then
				! do ii = 1,2**(level)
					! do jj =1,2**(level_butterfly-level+1)
						! deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
					! end do
				! end do	
				! deallocate(ButterflyP_old%blocks)
			! end if
			! allocate(ButterflyP_old%blocks(2**(level+1),2**(level_butterfly-level)))
			! do ii = 1,2**(level+1)
				! do jj =1,2**(level_butterfly-level)
					! mm = size(ButterflyP%blocks(ii,jj)%matrix,1)
					! nn = size(ButterflyP%blocks(ii,jj)%matrix,2)
					! allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn))
					! ButterflyP_old%blocks(ii,jj)%matrix = ButterflyP%blocks(ii,jj)%matrix
					! deallocate(ButterflyP%blocks(ii,jj)%matrix)
				! end do
			! end do
			! deallocate(ButterflyP%blocks)
		! else 
			! if(level_butterfly/=0)then
				! if(allocated(ButterflyP_old%blocks))then
					! do ii = 1,2**(level)
						! do jj =1,2**(level_butterfly-level+1)
							! deallocate(ButterflyP_old%blocks(ii,jj)%matrix)
						! end do
					! end do	
					! deallocate(ButterflyP_old%blocks)
				! end if					
			! end if
		! end if
		
		
		
    ! enddo

	! rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),rankmax_of_level(level_blocks))
	
	! ! write(*,*)rankmax_for_butterfly
	! ! write(*,*)'max value: ',maxvalue(1:9)
	
    ! deallocate (rankmax_for_butterfly)
    
    ! Memory=memory_butterfly
    ! !write (*,*) memory_butterfly
    ! !pause

    ! return

! end subroutine Butterfly_compress_ID



subroutine ACA_CompressionForward(matU,matV,Singular,header_m,header_n,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real*8 norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax,idxs_r,idxs_c
    complex(kind=8) value_Z,value_UV,maxvalue
    complex(kind=8) inner_U,inner_V,ctemp
    real*8 inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	complex(kind=8)::matU(rankmax_r,rmax),matV(rmax,rankmax_c)
	real*8::Singular(rmax)
    complex(kind=8),allocatable:: row_R(:),column_R(:)
    real*8,allocatable:: norm_row_R(:),norm_column_R(:)
	
	complex(kind=8), allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:),QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)	
	real*8, allocatable :: Singularsml(:)
	
	
	allocate(select_column(rankmax_c))
	allocate(select_row(rankmax_r))
	
		
	rankmax_min = min(rankmax_r,rankmax_c)
    norm_Z=0
	select_column = 0
	select_row = 0
	
    allocate(row_R(rankmax_c),column_R(rankmax_r))
    allocate(norm_row_R(rankmax_c),norm_column_R(rankmax_r))

    select_row(1)=1

    !$omp parallel do default(shared) private(j,value_Z,edge_m,edge_n)
    do j=1, rankmax_c
        ! value_Z=mat(select_row(1),j)
		edge_m = header_m + select_row(1) - 1
		edge_n = header_n + j - 1
		call element_Zmn(edge_m,edge_n,value_Z)							 
        row_R(j)=value_Z
        norm_row_R(j)=value_Z*conjg(value_Z)
    enddo
    !$omp end parallel do

    select_column(1)=maxloc(norm_row_R,1)
    maxvalue=row_R(select_column(1))

    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        row_R(j)=row_R(j)/maxvalue
    enddo
    ! !$omp end parallel do
    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        matV(1,j)=row_R(j)
    enddo
    ! !$omp end parallel do

    !$omp parallel do default(shared) private(i,value_Z,edge_m,edge_n)
    do i=1,rankmax_r
		edge_m = header_m + i - 1
		edge_n = header_n + select_column(1) - 1
		call element_Zmn(edge_m,edge_n,value_Z)
        ! value_Z=mat(i,select_column(1))
        column_R(i)=value_Z
        norm_column_R(i)=value_Z*conjg(value_Z)
    enddo
    !$omp end parallel do

    norm_column_R(select_row(1))=0

    ! !$omp parallel do default(shared) private(i)
    do i=1,rankmax_r
        matU(i,1)=column_R(i)
    enddo
    ! !$omp end parallel do

    norm_U=norm_vector(column_R,rankmax_r)
    norm_V=norm_vector(row_R,rankmax_c)
    norm_Z=norm_Z+norm_U*norm_V

	! if(rankmax<2)write(*,*)'rankmax'
    select_row(2)=maxloc(norm_column_R,1)

    rank=1
	! write(*,*)column_R,row_R
	! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
    do while (norm_Z*tolerance**2<norm_U*norm_V .and. rank<rankmax_min)

        !$omp parallel do default(shared) private(j,i,value_Z,value_UV,edge_m,edge_n)
        do j=1,rankmax_c		
			edge_m = header_m + select_row(rank+1) - 1
			edge_n = header_n + j - 1
			call element_Zmn(edge_m,edge_n,value_Z)	
            ! value_Z=mat(select_row(rank+1),j)
            value_UV=0
            do i=1,rank
                value_UV=value_UV+matU(select_row(rank+1),i)*matV(i,j)
            enddo
            row_R(j)=value_Z-value_UV
            norm_row_R(j)=row_R(j)*conjg(row_R(j))
        enddo
        !$omp end parallel do

        do i=1,rank
        	norm_row_R(select_column(i))=0
        enddo

        select_column(rank+1)=maxloc(norm_row_R,1)
        maxvalue=row_R(select_column(rank+1))

        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
    	row_R(j)=row_R(j)/maxvalue
        enddo
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
            matV(rank+1,j)=row_R(j)
        enddo
        ! !$omp end parallel do

        !$omp parallel do default(shared) private(i,j,value_Z,value_UV,edge_m,edge_n)
        do i=1,rankmax_r
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(rank+1) - 1
			call element_Zmn(edge_m,edge_n,value_Z)		
            ! value_Z=mat(i,select_column(rank+1))
            value_UV=0
            do j=1,rank
                value_UV=value_UV+matU(i,j)*matV(j,select_column(rank+1))
            enddo
            column_R(i)=value_Z-value_UV
            norm_column_R(i)=column_R(i)*conjg(column_R(i))
        enddo
        !$omp end parallel do

        do i=1,rank+1
            norm_column_R(select_row(i))=0
        enddo

        ! !$omp parallel do default(shared) private(i)
        do i=1,rankmax_r
            matU(i,rank+1)=column_R(i)
        enddo
        ! !$omp end parallel do

        norm_U=norm_vector(column_R,rankmax_r)
        norm_V=norm_vector(row_R,rankmax_c)
        inner_UV=0
        do j=1,rank
            inner_U=0
            inner_V=0
	        ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i=1,rankmax_r
                ctemp=matU(i,rank+1)*matU(i,j)
                inner_U=inner_U+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i=1,rankmax_c
                ctemp=matV(rank+1,i)*matV(j,i)
                inner_V=inner_V+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            inner_UV=inner_UV+2*sqrt(inner_U*inner_V)
        enddo
        norm_Z=norm_Z+inner_UV+norm_U*norm_V

        rank=rank+1
		if(rank+1>rmax)then
			write(*,*)'increase rmax',rank,rmax
			stop
		end if
        if (rank<rankmax_min) then
            select_row(rank+1)=maxloc(norm_column_R,1)
        endif

    enddo

	! write(*,*)select_row(1:rank),select_column(1:rank)
	
    deallocate(row_R,column_R)
    deallocate(norm_row_R,norm_column_R)

	
! ACA followed by SVD	
	
	allocate(QQ1(rankmax_r,rank))
	call copymatN_omp(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
	allocate (tau_Q(rank))
	call geqrff90(QQ1,tau_Q)
	
	allocate (RR1(rank,rank))
	RR1=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ1,tau_Q)
	deallocate(tau_Q)


	allocate(QQ2tmp(rankmax_c,rank))
	call copymatT_omp(matV(1:rank,1:rankmax_c),QQ2tmp,rank,rankmax_c)
	allocate (tau_Q(rank))
	call geqrff90(QQ2tmp,tau_Q)
	
	allocate (RR2tmp(rank,rank))
	RR2tmp=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR2tmp(i,j)=QQ2tmp(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ2tmp,tau_Q)
	deallocate(tau_Q)
	
	allocate(QQ2(rank,rankmax_c))
	call copymatT_omp(QQ2tmp,QQ2,rankmax_c,rank)	
	allocate(RR2(rank,rank))
	call copymatT_omp(RR2tmp,RR2,rank,rank)	
	
	
	
	! allocate(matU1(rankmax_r,rank))
	! allocate(matV1(rank,rankmax_c))
	! call gemm_omp(QQ1,RR1,matU1,rankmax_r,rank,rank)
	
	! call gemm_omp(RR2,QQ2,matV1,rank,rank,rankmax_c)
	
	! write(*,*)fnorm(matU1-matU(1:rankmax_r,1:rank),rankmax_r,rank),fnorm(matV1-matV(1:rank,1:rankmax_c),rank,rankmax_c)
	
	
	
	deallocate(QQ2tmp,RR2tmp)
	allocate(mattemp(rank,rank))
	call gemm_omp(RR1,RR2,mattemp,rank,rank,rank)
	
	
	
	
	allocate(UUsml(rank,rank),VVsml(rank,rank),Singularsml(rank))
	call SVD_Truncate(mattemp,rank,rank,rank,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew)
	
	call gemm_omp(QQ1,UUsml(1:rank,1:ranknew),matU(1:rankmax_r,1:ranknew),rankmax_r,rank,ranknew)
	
	call gemm_omp(VVsml(1:ranknew,1:rank),QQ2,matV(1:ranknew,1:rankmax_c),ranknew,rank,rankmax_c)
	! write(*,*)'aca rank:',rank,'after svd',ranknew
	
	rank = ranknew
	Singular(1:ranknew) = Singularsml(1:ranknew)
	
	deallocate(mattemp,RR1,RR2,QQ1,QQ2,UUsml,VVsml,Singularsml)

	deallocate(select_column)
	deallocate(select_row)
	
    return

end subroutine ACA_CompressionForward



subroutine ACA_CompressionForward_givenfullmat(matU,matV,Singular,header_m,header_n,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance,idx_m_ref,idx_n_ref)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real*8 norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax,idxs_r,idxs_c
    complex(kind=8) value_Z,value_UV,maxvalue
    complex(kind=8) inner_U,inner_V,ctemp
    real*8 inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	complex(kind=8)::matU(rankmax_r,rmax),matV(rmax,rankmax_c)
	real*8::Singular(rmax)
    complex(kind=8),allocatable:: row_R(:),column_R(:)
    real*8,allocatable:: norm_row_R(:),norm_column_R(:)
	
	complex(kind=8), allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:),QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)	
	real*8, allocatable :: Singularsml(:)
	
	integer idx_m_ref,idx_n_ref
		
	allocate(select_column(rankmax_c))
	allocate(select_row(rankmax_r))		
		
	rankmax_min = min(rankmax_r,rankmax_c)
    norm_Z=0
	select_column = 0
	select_row = 0
	
    allocate(row_R(rankmax_c),column_R(rankmax_r))
    allocate(norm_row_R(rankmax_c),norm_column_R(rankmax_r))

    select_row(1)=1

	! !$omp parallel do default(shared) private(j,edge_m,edge_n,value_Z)
    do j=1, rankmax_c
        ! value_Z=mat(select_row(1),j)
		edge_m = header_m + select_row(1) - 1
		edge_n = header_n + j - 1
		! call element_Zmn(edge_m,edge_n,value_Z)
		value_Z = matSub_glo(edge_m-idx_m_ref+1,edge_n-idx_n_ref+1)
        row_R(j)=value_Z
        norm_row_R(j)=value_Z*conjg(value_Z)
    enddo
    ! !$omp end parallel do

    select_column(1)=maxloc(norm_row_R,1)
    maxvalue=row_R(select_column(1))

    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        row_R(j)=row_R(j)/maxvalue
    enddo
    ! !$omp end parallel do
    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        matV(1,j)=row_R(j)
    enddo
    ! !$omp end parallel do

    ! !$omp parallel do default(shared) private(i,edge_m,edge_n,value_Z)
    do i=1,rankmax_r
		edge_m = header_m + i - 1
		edge_n = header_n + select_column(1) - 1
		! call element_Zmn(edge_m,edge_n,value_Z)	
		value_Z = matSub_glo(edge_m-idx_m_ref+1,edge_n-idx_n_ref+1)
        ! value_Z=mat(i,select_column(1))
        column_R(i)=value_Z
        norm_column_R(i)=value_Z*conjg(value_Z)
    enddo
    ! !$omp end parallel do

    norm_column_R(select_row(1))=0

    ! !$omp parallel do default(shared) private(i)
    do i=1,rankmax_r
        matU(i,1)=column_R(i)
    enddo
    ! !$omp end parallel do

    norm_U=norm_vector(column_R,rankmax_r)
    norm_V=norm_vector(row_R,rankmax_c)
    norm_Z=norm_Z+norm_U*norm_V

	! if(rankmax<2)write(*,*)'rankmax'
    select_row(2)=maxloc(norm_column_R,1)

    rank=1
	! write(*,*)column_R,row_R
	! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
    do while (norm_Z*tolerance**2<norm_U*norm_V .and. rank<rankmax_min)

		!$omp parallel do default(shared) private(j,i,edge_m,edge_n,value_Z,value_UV)
        do j=1,rankmax_c		
			edge_m = header_m + select_row(rank+1) - 1
			edge_n = header_n + j - 1
			! call element_Zmn(edge_m,edge_n,value_Z)
			value_Z = matSub_glo(edge_m-idx_m_ref+1,edge_n-idx_n_ref+1)			
            ! value_Z=mat(select_row(rank+1),j)
            value_UV=0
            do i=1,rank
                value_UV=value_UV+matU(select_row(rank+1),i)*matV(i,j)
            enddo
            row_R(j)=value_Z-value_UV
            norm_row_R(j)=row_R(j)*conjg(row_R(j))
        enddo
        !$omp end parallel do

        do i=1,rank
        	norm_row_R(select_column(i))=0
        enddo

        select_column(rank+1)=maxloc(norm_row_R,1)
        maxvalue=row_R(select_column(rank+1))

        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
    	row_R(j)=row_R(j)/maxvalue
        enddo
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
            matV(rank+1,j)=row_R(j)
        enddo
        ! !$omp end parallel do

        !$omp parallel do default(shared) private(j,i,edge_m,edge_n,value_Z,value_UV)
        do i=1,rankmax_r
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(rank+1) - 1
			! call element_Zmn(edge_m,edge_n,value_Z)	
			value_Z = matSub_glo(edge_m-idx_m_ref+1,edge_n-idx_n_ref+1)			
            ! value_Z=mat(i,select_column(rank+1))
            value_UV=0
            do j=1,rank
                value_UV=value_UV+matU(i,j)*matV(j,select_column(rank+1))
            enddo
            column_R(i)=value_Z-value_UV
            norm_column_R(i)=column_R(i)*conjg(column_R(i))
        enddo
        !$omp end parallel do

        do i=1,rank+1
            norm_column_R(select_row(i))=0
        enddo

        ! !$omp parallel do default(shared) private(i)
        do i=1,rankmax_r
            matU(i,rank+1)=column_R(i)
        enddo
        ! !$omp end parallel do

        norm_U=norm_vector(column_R,rankmax_r)
        norm_V=norm_vector(row_R,rankmax_c)
        inner_UV=0
        do j=1,rank
            inner_U=0
            inner_V=0
	        ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i=1,rankmax_r
                ctemp=matU(i,rank+1)*matU(i,j)
                inner_U=inner_U+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i=1,rankmax_c
                ctemp=matV(rank+1,i)*matV(j,i)
                inner_V=inner_V+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            inner_UV=inner_UV+2*sqrt(inner_U*inner_V)
        enddo
        norm_Z=norm_Z+inner_UV+norm_U*norm_V

        rank=rank+1
		
		! if(header_m==879 .and. header_n==3513)then
			! write(*,*)rank,norm_Z,inner_UV,norm_U*norm_V
		! end if
		
		
		if(rank+1>rmax)then
			write(*,*)'increase rmax',rank,rmax
			stop
		end if
        if (rank<rankmax_min) then
            select_row(rank+1)=maxloc(norm_column_R,1)
        endif

    enddo

	! write(*,*)select_row(1:rank),select_column(1:rank)
	
    deallocate(row_R,column_R)
    deallocate(norm_row_R,norm_column_R)

	
! ACA followed by SVD	
	
	allocate(QQ1(rankmax_r,rank))
	call copymatN_omp(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
	allocate (tau_Q(rank))
	call geqrff90(QQ1,tau_Q)
	
	allocate (RR1(rank,rank))
	RR1=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ1,tau_Q)
	deallocate(tau_Q)


	allocate(QQ2tmp(rankmax_c,rank))
	call copymatT_omp(matV(1:rank,1:rankmax_c),QQ2tmp,rank,rankmax_c)
	allocate (tau_Q(rank))
	call geqrff90(QQ2tmp,tau_Q)
	
	allocate (RR2tmp(rank,rank))
	RR2tmp=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR2tmp(i,j)=QQ2tmp(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ2tmp,tau_Q)
	deallocate(tau_Q)
	
	allocate(QQ2(rank,rankmax_c))
	call copymatT_omp(QQ2tmp,QQ2,rankmax_c,rank)	
	allocate(RR2(rank,rank))
	call copymatT_omp(RR2tmp,RR2,rank,rank)	
	
	
	
	! allocate(matU1(rankmax_r,rank))
	! allocate(matV1(rank,rankmax_c))
	! call gemm_omp(QQ1,RR1,matU1,rankmax_r,rank,rank)
	
	! call gemm_omp(RR2,QQ2,matV1,rank,rank,rankmax_c)
	
	! write(*,*)fnorm(matU1-matU(1:rankmax_r,1:rank),rankmax_r,rank),fnorm(matV1-matV(1:rank,1:rankmax_c),rank,rankmax_c)
	
	
	
	deallocate(QQ2tmp,RR2tmp)
	allocate(mattemp(rank,rank))
	call gemm_omp(RR1,RR2,mattemp,rank,rank,rank)
	
	
	
	
	allocate(UUsml(rank,rank),VVsml(rank,rank),Singularsml(rank))
	call SVD_Truncate(mattemp,rank,rank,rank,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew)
	
		! if(header_m==879 .and. header_n==3513)then
			! write(*,*)rank,ranknew,fnorm(mattemp,rank,rank),SVD_tolerance,'zhong',Singularsml,minloc(abs(Singularsml-Singularsml(1)*SVD_tolerance))
		! end if	
	
	
	call gemm_omp(QQ1,UUsml(1:rank,1:ranknew),matU(1:rankmax_r,1:ranknew),rankmax_r,rank,ranknew)
	
	call gemm_omp(VVsml(1:ranknew,1:rank),QQ2,matV(1:ranknew,1:rankmax_c),ranknew,rank,rankmax_c)
	! write(*,*)'aca rank:',rank,'after svd',ranknew
	
	rank = ranknew
	Singular(1:ranknew) = Singularsml(1:ranknew)
	
	deallocate(mattemp,RR1,RR2,QQ1,QQ2,UUsml,VVsml,Singularsml)

	deallocate(select_column)
	deallocate(select_row)
	
    return

end subroutine ACA_CompressionForward_givenfullmat


subroutine MVM_Z_forward(Ns,Vin,Vout,cascading_factors1)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer Ns
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o
	type(blockplus),pointer::bplus_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_m,idx_end_m,idx_start_n,idx_end_n
	
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8)::Vin(:),Vout(:)
	! complex(kind=8)::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(cascadingfactors)::cascading_factors1(:)
 
	num_vectors = 1   
	
	
	! get the right multiplied vectors
	ctemp1=1.0d0 ; ctemp2=1.0d0
	allocate(vec_old(Ns,num_vectors))
	allocate(vec_new(Ns,num_vectors))
	vec_old(1:Ns,1) = Vin
	vec_new = 0


	
	do level = 1,Maxlevel_for_blocks+1
		do ii =1, cascading_factors1(level)%N_block_forward
			bplus_o =>  cascading_factors1(level)%BP(ii)
			groupm = bplus_o%row_group
			groupn = bplus_o%col_group
			idx_start_m = basis_group(groupm)%head
			idx_end_m = basis_group(groupm)%tail				
			idx_start_n = basis_group(groupn)%head
			idx_end_n = basis_group(groupn)%tail
			
			if(level==Maxlevel_for_blocks+1)then
				call fullmat_block_MVP_randomized_dat(bplus_o%LL(1)%matrices_block(1),'N',idx_end_m-idx_start_m+1,num_vectors,&
				&vec_old(idx_start_n:idx_end_n,1:num_vectors),vec_new(idx_start_m:idx_end_m,1:num_vectors),ctemp1,ctemp2)
			else 
				call Bplus_block_MVP_randomized_dat(bplus_o,'N',idx_end_m-idx_start_m+1,idx_end_n-idx_start_n+1,num_vectors,&
				&vec_old(idx_start_n:idx_end_n,1:num_vectors),vec_new(idx_start_m:idx_end_m,1:num_vectors),ctemp1,ctemp2)
			end if
		end do				
	end do
	

	Vout = vec_new(1:Ns,1)
	deallocate(vec_old)
	deallocate(vec_new)	
	 
    return                

end subroutine MVM_Z_forward 



subroutine LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,SVD_tolerance,ButterflyP_old,ButterflyP)
use MISC
implicit none 
integer index_i_loc,index_j_loc,level_loc,level_butterflyL,index_i_m,index_i,index_j,level,group_m,mm,nn,nn1,nn2,j,i,mn,rank,mm1
type(butterfly_Kerl)ButterflyP_old,ButterflyP
complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
complex(kind=8), allocatable :: tau_Q(:)
real*8,allocatable :: Singular(:)
type(matrixblock)::blocks
real*8:: SVD_tolerance


	index_j = index_j_loc
	index_i = (index_i_m-1)*(2**level_loc)+index_i_loc
	group_m=blocks%row_group   ! Note: row_group and col_group interchanged here    
	group_m=group_m*2**level-1+index_i
	mm=basis_group(group_m)%tail-basis_group(group_m)%head+1
				
	! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
	if(size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,1)/=mm)then
		write(*,*)'mm incorrect'
		stop
	end if
	nn1=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix,2)
	nn2=size(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix,2)
	nn = nn1+nn2


	allocate(QQ(mm,nn))
	! !$omp parallel do default(shared) private(i,j)
	do j=1, nn1
		do i=1, mm
			QQ(i,j)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(i,j)
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do j=1, nn2
		do i=1, mm
			QQ(i,j+nn1)=ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(i,j)
		enddo
	enddo
	! !$omp end parallel do						



	mn=min(mm,nn)
	allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
	call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			
	! rank = min(rank,37)

	! rank = 7
	! write(*,*)rank
	
	! rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
	! rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))

	! blocks%rankmax = max(blocks%rankmax,rank)
	! blocks%rankmin = min(blocks%rankmin,rank)

	allocate(mat_tmp(mm,rank))
	! !$omp parallel do default(shared) private(i,j,k,ctemp)
	do j=1, rank
		do i=1, mm
			mat_tmp(i,j)=UU(i,j)*Singular(j)
		enddo
	enddo
	! !$omp end parallel do							



	if(level_loc/=level_butterflyL)then
		mm1 = basis_group(group_m*2)%tail-basis_group(group_m*2)%head+1
		allocate(ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank))
		allocate(ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank))
		
		ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
		ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)						
	else 
		allocate (blocks%ButterflyU(index_i)%matrix(mm,rank))
		! !$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, mm
				blocks%ButterflyU(index_i)%matrix(i,j)=mat_tmp(i,j)
				! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		! !$omp end parallel do

	end if		

	allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1))
	allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2))

	! !$omp parallel do default(shared) private(i,j)
	do i=1, rank
		do j=1, nn1
			blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=VV(i,j)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do i=1, rank
		do j=1, nn2
			blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
		enddo
	enddo
	! !$omp end parallel do						

	deallocate(QQ,UU,VV,Singular,mat_tmp)

	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
	! if (level_loc==level_butterflyL) then
		! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU(index_i)%matrix)/1024.0d3
	! endif


end subroutine  LocalButterflySVD_Left







subroutine LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,SVD_tolerance,ButterflyP_old,ButterflyP)
use MISC
implicit none 
integer index_i_loc,index_j_loc,level_loc,level_butterflyR,level_butterfly,index_j_m,index_i,index_j,level,group_n,mm,nn,nn1,mm1,mm2,j,i,mn,rank
type(butterfly_Kerl)ButterflyP_old,ButterflyP
complex(kind=8), allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
complex(kind=8), allocatable :: tau_Q(:)
real*8,allocatable :: Singular(:)
type(matrixblock)::blocks
real*8:: SVD_tolerance


	index_i = index_i_loc
	index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
	group_n=blocks%col_group   ! Note: row_group and col_group interchanged here    
	group_n=group_n*2**(level_butterfly-level+1)-1+index_j
	nn=basis_group(group_n)%tail-basis_group(group_n)%head+1
				
	! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
	if(size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,2)/=nn)then
		write(*,*)'nn incorrect'
		stop
	end if
	mm1=size(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix,1)
	mm2=size(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix,1)
	mm = mm1+mm2
	
	!!!!!!!!!!!!!!!!!!
	
	
	allocate(QQ(mm,nn))
	! !$omp parallel do default(shared) private(i,j)
	do j=1, nn
		do i=1, mm1
			QQ(i,j)=ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(i,j)
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do j=1, nn
		do i=1, mm2
			QQ(i+mm1,j)=ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(i,j)
		enddo
	enddo
	! !$omp end parallel do						
	
	
	
	mn=min(mm,nn)
	allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
	call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,SVD_tolerance_forward,rank)			
	! rank = min(rank,37)

	! rank = 7
	! write(*,*)rank
	
	! rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
	! rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
	! blocks%rankmax = max(blocks%rankmax,rank)
	! blocks%rankmin = min(blocks%rankmin,rank)

	allocate(mat_tmp(rank,nn))
	! !$omp parallel do default(shared) private(i,j)
	do j=1, nn
		do i=1, rank
			mat_tmp(i,j)=VV(i,j)*Singular(i)
		enddo
	enddo
	! !$omp end parallel do							
	


	if(level_loc/=level_butterflyR)then
		nn1 = basis_group(group_n*2)%tail-basis_group(group_n*2)%head+1
		allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1))
		allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1))
		
		ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
		ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)						
	else 
		allocate (blocks%ButterflyV(index_j)%matrix(nn,rank))
		! !$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, nn
				blocks%ButterflyV(index_j)%matrix(i,j)=mat_tmp(j,i)
				! blocks%ButterflyU(index_i)%matrix(i,j)=random_complex_number()
			enddo
		enddo
		! !$omp end parallel do

	end if		

	allocate (blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm1,rank))
	allocate (blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm2,rank))
	
	! !$omp parallel do default(shared) private(i,j)
	do i=1, mm1
		do j=1, rank
			blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(i,j)=UU(i,j)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_complex_number()
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do i=1, mm2
		do j=1, rank
			blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_complex_number()
		enddo
	enddo
	! !$omp end parallel do						
	
	deallocate(QQ,UU,VV,Singular,mat_tmp)
  
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
	! if (level_loc==level_butterflyR) then
		! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV(index_j)%matrix)/1024.0d3
	! endif

end subroutine LocalButterflySVD_Right



end module Butterfly_compress_forward
