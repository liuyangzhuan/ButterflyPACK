#include "HODLR_config.fi"
module Bplus_compress 
use Bplus_randomized 
use HODLR_structure
! use element_Z 


contains 
! subroutine Butterfly_compress(blocks,Memory,msh,element_Zmn,ker)

   ! use HODLR_DEFS
   ! implicit none
	! type(mesh)::msh
	! type(kernelquant)::ker
    ! integer i, j, level_butterfly, num_blocks, k, attempt
    ! integer group_m, group_n, mm, nn, index_i, index_j, ii, jj
    ! integer level, length_1, length_2, level_blocks, index_ij
    ! integer rank, rankmax, butterflyB_inuse, rank1, rank2
    ! real(kind=8) rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    ! real(kind=8) Memory
    ! DT ctemp
	! type(matrixblock)::blocks
	! procedure(Z_elem)::element_Zmn
	
	! integer, allocatable :: rankmax_for_butterfly(:),rankmin_for_butterfly(:)
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

    ! allocate(blocks%ButterflyU%blocks(num_blocks))
    ! allocate(blocks%ButterflyV%blocks(num_blocks))
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
                ! call decomposition_UkerlV(index_i,index_j,level,blocks,tolerance,msh,element_Zmn)
                ! index_ij=index_ij+1
                ! if (level==0) then
                    ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_ij)%matrix)/1024.0d3
                ! else                    
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
                ! endif
                ! if (level==level_butterfly) then
					! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_ij)%matrix)/1024.0d3
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

! subroutine decomposition_UkerlV(index_i,index_j,level,blocks,tolerance,msh,element_Zmn,ker)

    ! use HODLR_DEFS
    ! implicit none

    ! integer mm, nn, level, level_butterfly
    ! integer i, j, ii, jj, index_i, index_j
    ! real(kind=8) tolerance
	! type(matrixblock)::blocks
	! type(mesh)::msh
	! type(kernelquant)::ker
	! procedure(Z_elem)::element_Zmn
	
    ! if (level==0) then
        ! call butterfly_recomposition_FastSampling_initial(index_j,blocks,msh,element_Zmn)
    ! else
        ! call butterfly_recomposition_FastSampling(index_i,index_j,level,blocks,msh,element_Zmn)
    ! endif

    ! return

! end subroutine decomposition_UkerlV

! subroutine butterfly_recomposition_FastSampling_initial(index_j,blocks,msh,element_Zmn,ker)

    ! use HODLR_DEFS
    ! 
    ! implicit none

    ! integer i, j, k, level_butterfly, level_blocks, ii, jj
    ! integer group_m, group_n, mm, nn, nn_start, edge_m, edge_n
    ! integer level, length_1, length_2, flag
    ! integer rank, rankmax,rankmax_c,rankmax_r, rankmax_min,index_i, index_j, flag0, rank_new, header_m, header_n
    ! real(kind=8) rate, tolerance, rtemp1, rtemp2, rtemp
    ! DT ctemp

    ! integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    ! integer,allocatable:: select_row_rr(:), select_column_rr(:)
    ! DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:)
    ! real(kind=8),allocatable:: Singular(:)
	! integer Nlayer
	! type(matrixblock)::blocks
	! type(mesh)::msh
	! type(kernelquant)::ker
    ! procedure(Z_elem)::element_Zmn
	
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

    ! rankmax=rank_approximate_func(group_m, group_n, 1,ker)

    ! if (rankmax>nn) then
        ! rankmax=nn
    ! endif

	! ! modified by Yang, this may improve accuracy as one dimension is constant
	! rankmax=min(mm,nn)	

	! rankmax_r = min(rank_approximate_func(group_m, group_n, 1)*2,mm,ker)
	! rankmax_c = nn
	! rankmax_min = min(rankmax_r,rankmax_c)
	
	! ! write(*,*)rankmax,rank_approximate_func(group_m, group_n, 1,ker),'i'
		
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
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
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
                    ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
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
                    ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
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

						
        ! allocate (blocks%ButterflyU%blocks(index_j)%matrix(mm,rank_new))
        ! !$omp parallel do default(shared) private(i,j,k,ctemp)
        ! do j=1, rank_new
            ! do i=1, mm
                ! ctemp=0
                ! do k=1, rankmax_c
                    ! ctemp=ctemp+matrix_U(i,k)*conjg(VV(j,k))
                ! enddo
                ! blocks%ButterflyU%blocks(index_j)%matrix(i,j)=ctemp
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_U,VV)
        
        ! allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn,rank_new))
        ! !$omp parallel do default(shared) private(i,j,k,ctemp)
         ! do j=1, rank_new
            ! do i=1, nn
                ! ctemp=0
                ! do k=1, rankmax_r
                    ! ctemp=ctemp+conjg(UU(k,j))*matrix_V(i,k)
                ! enddo
                ! blocks%ButterflyV%blocks(index_j)%matrix(i,j)=ctemp/Singular(j)
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_V,UU,Singular)
    
    ! else
        
        ! allocate (matrixtemp_U(rankmax_r,rankmax_c),matrixtemp_V(rankmax_c,rankmax_r))
        ! allocate (select_column_rr(rankmax_c),select_row_rr(rankmax_r))
        ! ! write(*,*)rankmax_r,rankmax_c
		! call ACA_SubsetSelection(MatrixSubselection,select_column_rr,select_row_rr,rankmax_r,rankmax_c,rank_new,option%tol_comp)
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
                    ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
                    ! matrix_V(i,j)=ctemp
                ! endif
            ! enddo
        ! enddo
        ! !$omp end parallel do
        
        ! deallocate (MatrixSubselection)
        

		! allocate(matrix_little_inv(rank_new,rank_new))
		! allocate (matrix_V_tmp(rank_new,nn))
		! call GeneralInverse(rank_new,rank_new,matrix_little,matrix_little_inv,ACA_tolerance_forward)
		! call gemm_omp(matrix_little_inv,matrix_V,matrix_V_tmp,rank_new,nn,rank_new)
		! matrix_V = matrix_V_tmp
		! deallocate(matrix_little_inv)
		! deallocate(matrix_V_tmp)		
		
		
        ! ! ! allocate (column_pivot(rank_new))
        ! ! ! call getrff90(matrix_little,column_pivot)
        ! ! ! call getrsf90(matrix_little,column_pivot,matrix_V,'N')
        ! ! ! deallocate (column_pivot,matrix_little)
        
		
        ! allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn,rank_new))
        ! !$omp parallel do default(shared) private(i,j)
         ! do j=1, rank_new
            ! do i=1, nn
                ! blocks%ButterflyV%blocks(index_j)%matrix(i,j)=matrix_v(j,i)
            ! enddo
        ! enddo
        ! !$omp end parallel do
        ! deallocate (matrix_V)
        
        ! deallocate (select_row, select_row_rr, select_column, select_column_rr)    
        
    ! endif

    ! return

! end subroutine butterfly_recomposition_FastSampling_initial
 
! subroutine butterfly_recomposition_FastSampling(index_i,index_j,level,blocks,msh,element_Zmn,ker)

    ! use HODLR_DEFS
    ! 
    ! implicit none

    ! integer i, j, level_butterfly, level_blocks, num_groupm, num_groupn, k, ii, jj, mmm, nnn, nnn1
    ! integer group_m, group_n, mm, nn, index_m, index_n, option, mm_start, nn1, nn2, header_m, header_n
    ! integer level, length_1, length_2, butterflyB_inuse, index_ii, index_jj, header_n1, header_n2
    ! integer rank, rankmax,rankmax_c,rankmax_r,rankmax_min,index_i, index_j, index_ij, index_iijj, flag, flag0, rank_new, edge_m, edge_n
    ! real(kind=8) rate, tolerance, rtemp1, rtemp2, rtemp
    ! DT ctemp

    ! integer,allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
    ! integer,allocatable:: select_row_rr(:), select_column_rr(:)
    ! DT,allocatable:: UU(:,:), VV(:,:), matrix_little(:,:),matrix_little_inv(:,:), matrix_U(:,:), matrix_V(:,:),matrix_V_tmp(:,:), matrix_little_cc(:,:)
    ! real(kind=8),allocatable:: Singular(:)
	! type(matrixblock)::blocks
	! type(mesh)::msh
	! type(kernelquant)::ker
	! procedure(Z_elem)::element_Zmn
	
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

    ! rankmax=rank_approximate_func(group_m, group_n, 1,ker)

    ! if (rankmax>nn) then
        ! rankmax=nn
    ! endif

	! ! modified by Yang, this may improve accuracy as one dimension is constant
	! rankmax=min(mm,nn)	

	! rankmax_r = min(rank_approximate_func(group_m, group_n, 1)*2,mm,ker)
	! rankmax_c = nn
	! rankmax_min = min(rankmax_r,rankmax_c)
	
    ! ! if (rank_control_forward/=0) then
        ! ! if (rankmax_for_butterfly(level-1)<rankmax) then
            ! ! rankmax=rankmax_for_butterfly(level-1)
        ! ! endif
    ! ! endif

	! if(rankmax==1)write(*,*)group_m,group_n,rank_approximate_func(group_m, group_n, 1,ker),mm,nn
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
                ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
                ! MatrixSubselection(i,j)=ctemp
            ! else
                ! edge_m=select_row(i)+header_m-1
                ! edge_n=blocks%ButterflyColSelect(index_iijj+1,level-1)%select_columns(select_column(j)-nn1)+header_n2-1
                ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
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
    ! call ACA_SubsetSelection(MatrixSubselection,select_column_rr,select_row_rr,rankmax_r,rankmax_c,rank_new,option%tol_comp)
    
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
        ! allocate (blocks%ButterflyU%blocks(index_ij)%matrix(mm,rank_new)) 
        ! !$omp parallel do default(shared) private(i,j,ii,edge_m,edge_n,ctemp)
        ! do j=1, rank_new
            ! ii=1
            ! do i=1, mm
                ! if (i==select_row(ii)) then
                    ! blocks%ButterflyU%blocks(index_ij)%matrix(i,j)=MatrixSubselection(ii,select_column_rr(j))
                    ! if (ii<rankmax_r) then
                        ! ii=ii+1
                    ! endif
                ! else
                    ! edge_m=i+header_m-1
                    ! edge_n=blocks%ButterflyColSelect(index_ij,level)%select_columns(j)+header_n1-1 
                    ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)                    
                    ! blocks%ButterflyU%blocks(index_ij)%matrix(i,j)=ctemp
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
                ! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)                    
                ! matrix_V(i,j)=ctemp
            ! endif
        ! enddo
    ! enddo
    ! !$omp end parallel do
    
    ! deallocate (MatrixSubselection)
    
	! allocate(matrix_little_inv(rank_new,rank_new))
	! allocate (matrix_V_tmp(rank_new,nn))
	! call GeneralInverse(rank_new,rank_new,matrix_little,matrix_little_inv,ACA_tolerance_forward)
	! call gemm_omp(matrix_little_inv,matrix_V,matrix_V_tmp,rank_new,nn,rank_new)
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



subroutine Bplus_compress_N15(bplus,option,Memory,stats,msh,ker,element_Zmn,ptree)

   use HODLR_DEFS
   
   
   use misc
   implicit none

    type(blockplus)::bplus
	integer:: ii,ll,bb
    real(kind=8) Memory,rtemp	
	integer:: level_butterfly,level_BP,levelm,groupm_start,Nboundall
	type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	
	Memory = 0
	
	do ll=1,bplus%Lplus
		bplus%LL(ll)%rankmax=0
															   
		do bb = 1,bplus%LL(ll)%Nbound
			! write(*,*)bplus%level,ll,bb
			if(bplus%LL(ll+1)%Nbound==0)then
				
				! bplus%LL(ll)%matrices_block(bb)%level_butterfly = int((Maxlevel_for_blocks-bplus%LL(ll)%matrices_block(bb)%level)/2)*2
				! if(option%TwoLayerOnly==1 .and. bplus%Lplus==2)bplus%LL(ll)%matrices_block(bb)%level_butterfly = 0
				call Butterfly_compress_N15(bplus%LL(ll)%matrices_block(bb),option,rtemp,stats,msh,ker,element_Zmn,ptree)
				call Butterfly_sym2asym(bplus%LL(ll)%matrices_block(bb))
				Memory = Memory + rtemp
			else 		
				level_butterfly = bplus%LL(ll)%matrices_block(1)%level_butterfly
				level_BP = bplus%level			
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				groupm_start=bplus%LL(ll)%matrices_block(1)%row_group*2**levelm		
				Nboundall = 2**(bplus%LL(ll)%matrices_block(1)%level+levelm-level_BP)			
				call Butterfly_compress_N15_withoutBoundary(bplus%LL(ll)%matrices_block(bb),bplus%LL(ll+1)%boundary_map,Nboundall,groupm_start, option, rtemp,stats,msh,ker,element_Zmn,ptree)
				call Butterfly_sym2asym(bplus%LL(ll)%matrices_block(bb))
				Memory = Memory + rtemp
			end if	
			bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax,bplus%LL(ll)%matrices_block(bb)%rankmax)			
		end do
	end do
	! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'
	
    return

end subroutine Bplus_compress_N15



subroutine Butterfly_compress_N15_withoutBoundary(blocks,boundary_map,Nboundall, groupm_start,option,Memory,stats,msh,ker,element_Zmn,ptree)

   use HODLR_DEFS
   
   
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
    real(kind=8) rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real(kind=8) Memory
    DT ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	DT, allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	DT, allocatable :: tau_Q(:)
	real(kind=8),allocatable :: Singular(:)
	integer flag,tmpi,tmpj
	real(kind=8):: a,b,n1,n2,error
	real(kind=8):: maxvalue(1:20)
	real(kind=8):: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row,frow
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	
	integer cnt_tmp,rankFar,rankNear
	type(Hoption)::option
	type(Hstat)::stats
	integer, allocatable :: rankmax_for_butterfly(:),rankmin_for_butterfly(:)
	
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	! ForwardSymmetricFlag = 1	
	maxvalue = 0
	minvalue = 10000
    Memory=0

    level_blocks=blocks%level
    !level_butterfly=Maxlevel-level_blocks
    level_butterfly=blocks%level_butterfly
																  
	
	call assert(level_butterfly>=2,'level_butterfly should be at least 2')
	

!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
	blocks%rankmax = -100000
	blocks%rankmin = 100000
	
    ! blocks%level_butterfly=level_butterfly


    num_blocks=2**level_butterfly

    allocate(blocks%ButterflyU%blocks(num_blocks))
    allocate(blocks%ButterflyV%blocks(num_blocks))
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
			frow=1
			! call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ptree)					
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
				frow=1
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ker,stats,element_Zmn,ptree,error)	
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
			allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank));ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%mdim=mm1;ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%ndim=rank
			allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank));ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%mdim=mm-mm1;ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%ndim=rank
			
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
				call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,option,ButterflyP_old,ButterflyP)				
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
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_i)%matrix)/1024.0d3
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
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn));ButterflyP_old%blocks(ii,jj)%mdim=mm;ButterflyP_old%blocks(ii,jj)%ndim=nn
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
		! time_tmp = time_tmp + n2 - n1
		
	enddo

	
	! write(*,*)'stats:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	


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
				frow=1
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ker,stats,element_Zmn,ptree,error)	
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
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1));ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%mdim=rank;ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%ndim=nn1
			allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1));ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%mdim=rank;ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%ndim=nn-nn1
			
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
				call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,option,ButterflyP_old,ButterflyP)				
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
					memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
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
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn));ButterflyP_old%blocks(ii,jj)%mdim=mm;ButterflyP_old%blocks(ii,jj)%ndim=nn
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
		! time_tmp = time_tmp + n2 - n1

	enddo		

	! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
	stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),stats%rankmax_of_level(level_blocks))
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


subroutine Butterfly_compress_N15(blocks,option,Memory,stats,msh,ker,element_Zmn,ptree)

   use HODLR_DEFS
   
   
   use misc
   implicit none

    integer blocks_idx, i, j, level_butterfly,level_butterflyL,level_butterflyR, num_blocks, k, attempt
    integer group_m, group_n, mm, nn,mn, index_i, index_j,index_ij_loc,index_i_m, index_j_m,index_i_loc, index_j_loc, ii, jj,nn1,nn2,mm1,mm2,idxs_m,idxs_n
    integer level,levelm,level_loc, length_1, length_2, level_blocks, index_ij,edge_m,edge_n
    integer rank, rankmax, butterflyB_inuse, rank1, rank2,rmax
    real(kind=8) rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
    real(kind=8) Memory,n1,n2
    DT ctemp
	type(butterfly_Kerl)ButterflyP_old,ButterflyP
	DT, allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
	DT, allocatable :: tau_Q(:)
	real(kind=8),allocatable :: Singular(:)
	integer flag,tmpi,tmpj
	real(kind=8):: a,b,error
	real(kind=8):: maxvalue(1:20)
	real(kind=8):: minvalue(1:20)
	integer dimension_m,dimension_n,dimension_rank,num_col,num_row
	type(matrixblock)::blocks
	integer cnt_tmp,rankFar,rankNear,leafsize,frow
	type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	integer, allocatable :: rankmax_for_butterfly(:),rankmin_for_butterfly(:)
	cnt_tmp	= 0
	rankFar = 0
	rankNear = 0
	! ForwardSymmetricFlag = 1	
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

    allocate(blocks%ButterflyU%blocks(num_blocks))
    allocate(blocks%ButterflyV%blocks(num_blocks))
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
				! call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
				! QQ(ii,jj) = ctemp
			! end do
		! end do
		! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			
		! deallocate(QQ)

		
		
		
		
		! !!!!! ACA-SVD
		! idxs_m = basis_group(group_m)%head
		! idxs_n = basis_group(group_n)%head
		
		! rmax = min(500,min(mm,nn))
		! allocate(UU(mm,rmax))
		! allocate(VV(rmax,nn))
		! allocate(Singular(rmax))
		! frow=1
		! call ACA_CompressionForward(UU,VV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ker,element_Zmn,ptree)		
		
		
		
		
		! ! rank = 7
		! rankmax_for_butterfly(0)=max(rank,rankmax_for_butterfly(0))
		! rankmin_for_butterfly(0)=min(rank,rankmin_for_butterfly(0))

		! blocks%rankmax = max(blocks%rankmax,rank)
		! blocks%rankmin = min(blocks%rankmin,rank)
		
		! allocate (blocks%ButterflyV%blocks(1)%matrix(nn,rank));blocks%ButterflyV%blocks(1)%mdim=nn;blocks%ButterflyV%blocks(1)%ndim=rank
		
		! !$omp parallel do default(shared) private(i,j)
		! do j=1, rank
			! do i=1, nn
				! blocks%ButterflyV%blocks(1)%matrix(i,j)=VV(j,i)
				! ! blocks%ButterflyV%blocks(1)%matrix(i,j)=random_dp_number()
			! enddo
		! enddo
		! !$omp end parallel do	
		
		
		! allocate (blocks%ButterflyU%blocks(1)%matrix(mm,rank));blocks%ButterflyU%blocks(1)%mdim=mm;blocks%ButterflyU%blocks(1)%ndim=rank
		
		! !$omp parallel do default(shared) private(i,j)
		! do j=1, rank
			! do i=1, mm
				! blocks%ButterflyU%blocks(1)%matrix(i,j)=UU(i,j)*Singular(j)
				! ! blocks%ButterflyU%blocks(1)%matrix(i,j)=random_dp_number()
			! enddo
		! enddo
		! !$omp end parallel do							
		! deallocate (UU,VV,Singular)		

		 
		

		!!!! parallel ACA
		
		! do ii=1,5
		! leafsize = max(blocks%M,blocks%N)/2**ii
		leafsize = max(blocks%M,blocks%N)
		
		! leafsize = 2502
		if(allocated(blocks%ButterflyU%blocks(1)%matrix))deallocate(blocks%ButterflyU%blocks(1)%matrix)
		if(allocated(blocks%ButterflyV%blocks(1)%matrix))deallocate(blocks%ButterflyV%blocks(1)%matrix)
		call BlockLR(blocks,leafsize,rank,option,msh,ker,stats,element_Zmn,ptree,blocks%pgno,ptree%pgrp(blocks%pgno)%gd,0)		
		! enddo
		
		rankmax_for_butterfly(0)=max(blocks%rankmax,rankmax_for_butterfly(0))
		rankmin_for_butterfly(0)=min(blocks%rankmin,rankmin_for_butterfly(0))		
		
		

		!!!! pseudo skeleton

		! ! ! call SeudoSkeleton_CompressionForward(blocks,blocks%headm,blocks%headn,mm,nn,min(min(mm,nn),1000),min(min(mm,nn),1000),rank,option%tol_comp,option%tol_comp,msh,ker,stats,element_Zmn,ptree,ptree%pgrp(blocks%pgno)%ctxt,blocks%pgno)
		! call SeudoSkeleton_CompressionForward(blocks,blocks%headm,blocks%headn,mm,nn,min(mm,nn),min(mm,nn),rank,option%tol_comp,option%tol_comp,msh,ker,stats,element_Zmn,ptree,ptree%pgrp(blocks%pgno)%ctxt,blocks%pgno)
		! rankmax_for_butterfly(0)=max(blocks%rankmax,rankmax_for_butterfly(0))
		! rankmin_for_butterfly(0)=min(blocks%rankmin,rankmin_for_butterfly(0))

		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(1)%matrix)/1024.0d3
		memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(1)%matrix)/1024.0d3
		
	
		! write(*,*)fnorm(blocks%ButterflyU%blocks(1)%matrix,mm,rank),fnorm(blocks%ButterflyV%blocks(1)%matrix,nn,rank),'ganga'
		! stop
				
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
								call element_Zmn(edge_m,edge_n,ctemp,msh,ker)
								QQ(ii,jj) = ctemp
							end do
						end do

						mn=min(mm,nn)
						allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
						call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			

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
						allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank));ButterflyP%blocks(2*index_i-1,index_j)%mdim=mm1;ButterflyP%blocks(2*index_i-1,index_j)%ndim=rank
						allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank));ButterflyP%blocks(2*index_i,index_j)%mdim=mm-mm1;ButterflyP%blocks(2*index_i,index_j)%ndim=rank
						
						ButterflyP%blocks(2*index_i-1,index_j)%matrix=mat_tmp(1:mm1,1:rank)
						ButterflyP%blocks(2*index_i,index_j)%matrix=mat_tmp(1+mm1:mm,1:rank)
					
						allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn,rank));blocks%ButterflyV%blocks(index_j)%mdim=nn;blocks%ButterflyV%blocks(index_j)%ndim=rank
						
						!$omp parallel do default(shared) private(i,j)
						do j=1, rank
							do i=1, nn
								blocks%ButterflyV%blocks(index_j)%matrix(i,j)=VV(j,i)
								! blocks%ButterflyV%blocks(index_j)%matrix(i,j)=random_dp_number()
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
						call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			

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
							allocate(ButterflyP%blocks(2*index_i-1,index_j)%matrix(mm1,rank));ButterflyP%blocks(2*index_i-1,index_j)%mdim=mm1;ButterflyP%blocks(2*index_i-1,index_j)%ndim=rank
							allocate(ButterflyP%blocks(2*index_i,index_j)%matrix(mm-mm1,rank));ButterflyP%blocks(2*index_i,index_j)%mdim=mm-mm1;ButterflyP%blocks(2*index_i,index_j)%ndim=rank
							
							ButterflyP%blocks(2*index_i-1,index_j)%matrix=mat_tmp(1:mm1,1:rank)
							ButterflyP%blocks(2*index_i,index_j)%matrix=mat_tmp(1+mm1:mm,1:rank)						
						else 
							allocate (blocks%ButterflyU%blocks(index_i)%matrix(mm,rank));blocks%ButterflyU%blocks(index_i)%mdim=mm;blocks%ButterflyU%blocks(index_i)%ndim=rank
							!$omp parallel do default(shared) private(i,j)
							do j=1, rank
								do i=1, mm
									blocks%ButterflyU%blocks(index_i)%matrix(i,j)=mat_tmp(i,j)
									! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
								enddo
							enddo
							!$omp end parallel do

						end if		

						allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(rank,nn1));blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%mdim=rank;blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%ndim=nn1
						allocate (blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(rank,nn2));blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%mdim=rank;blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%ndim=nn2
						
						!$omp parallel do default(shared) private(i,j)
						do i=1, rank
							do j=1, nn1
								blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=VV(i,j)
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
							enddo
						enddo
						!$omp end parallel do						
						!$omp parallel do default(shared) private(i,j)
						do i=1, rank
							do j=1, nn2
								blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
								! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
							enddo
						enddo
						!$omp end parallel do						
						
						deallocate(QQ,UU,VV,Singular,mat_tmp)

					end if
					
					index_ij=index_ij+1
					if (level==0) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_ij)%matrix)/1024.0d3
					else                    
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
					endif
					if (level==level_butterfly) then
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_ij)%matrix)/1024.0d3
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
						allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn));ButterflyP_old%blocks(ii,jj)%mdim=mm;ButterflyP_old%blocks(ii,jj)%ndim=nn
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
				
				frow=1
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ker,stats,element_Zmn,ptree,error)					
				! rank = min(rank,37)
				
				
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
					
					
					allocate(blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(rank,rank));blocks%ButterflyMiddle(index_i_m,index_j_m)%mdim=rank;blocks%ButterflyMiddle(index_i_m,index_j_m)%ndim=rank
					blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix = 0
					do ii=1,rank
						blocks%ButterflyMiddle(index_i_m,index_j_m)%matrix(ii,ii) = 1d0/Singular(ii)
					end do
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
				allocate(ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank));ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%mdim=mm1;ButterflyP_old%blocks(2*index_i_loc-1,index_j_loc)%ndim=rank
				allocate(ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank));ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%mdim=mm-mm1;ButterflyP_old%blocks(2*index_i_loc,index_j_loc)%ndim=rank
				
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

				! do index_i_loc=1, 2**level_loc
					! do index_j_loc=1, 2**(level_butterflyL-level_loc)
				
				! write(*,*)'addaa111111'

				
				!$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
				do index_ij_loc = 1, 2**level_butterflyL
					index_j_loc = mod(index_ij_loc-1,2**(level_butterflyL-level_loc))+1
					index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL-level_loc)))
				
					call LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,option,ButterflyP_old,ButterflyP)				
				enddo
				!$omp end parallel do
				
				
				! write(*,*)'addaa1111112222'
				
				
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
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_i)%matrix)/1024.0d3
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
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn));ButterflyP_old%blocks(ii,jj)%mdim=mm;ButterflyP_old%blocks(ii,jj)%ndim=nn
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
			! time_tmp = time_tmp + n2 - n1

		enddo

		
		! write(*,*)'stats:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear	


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
				frow=1
				call ACA_CompressionForward(matU,matV,Singular,idxs_m,idxs_n,mm,nn,frow,rmax,rank,option%tol_comp*0.1,option%tol_comp,msh,ker,stats,element_Zmn,ptree,error)	
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
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1));ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%mdim=rank;ButterflyP_old%blocks(index_i_loc,2*index_j_loc-1)%ndim=nn1
				allocate(ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1));ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%mdim=rank;ButterflyP_old%blocks(index_i_loc,2*index_j_loc)%ndim=nn-nn1
				
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
					call LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,option,ButterflyP_old,ButterflyP)				
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
						memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
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
						! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			
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
							! allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn,rank))
							! !$omp parallel do default(shared) private(i,j)
							! do j=1, rank
								! do i=1, nn
									! blocks%ButterflyV%blocks(index_j)%matrix(i,j)=mat_tmp(j,i)
									! ! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
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
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						! !$omp parallel do default(shared) private(i,j)
						! do i=1, mm2
							! do j=1, rank
								! blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
								! ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
							! enddo
						! enddo
						! !$omp end parallel do						
						
						! deallocate(QQ,UU,VV,Singular,mat_tmp)
					  
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
						! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
						! if (level_loc==level_butterflyR) then
							! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
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
							allocate(ButterflyP_old%blocks(ii,jj)%matrix(mm,nn));ButterflyP_old%blocks(ii,jj)%mdim=mm;ButterflyP_old%blocks(ii,jj)%ndim=nn
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
			! time_tmp = time_tmp + n2 - n1
		enddo		

			


	end if	
	 
	! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
	stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly),stats%rankmax_of_level(level_blocks))
	
	! write(*,*)stats%rankmax_of_level,'nitaima',rankmax_for_butterfly
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



recursive subroutine BlockLR(blocks,leafsize,rank,option,msh,ker,stats,element_Zmn,ptree,pgno,gd,cridx)
use HODLR_DEFS
implicit none 
    integer rank,ranktmp,leafsize
    integer header_m, header_n
    integer N,M,i,j,ii,jj,myi,myj,iproc,jproc,rmax,mn
	type(mesh)::msh
	type(Hoption)::option
	type(kernelquant)::ker
	type(matrixblock)::blocks,blockc(2)
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	integer pgno
	type(grid),pointer::gd
	type(grid),pointer::gdc1,gdc2
	integer:: cridx,info
	DT,allocatable:: UU(:,:), VV(:,:),matU(:,:),matV(:,:),matU1(:,:),matV1(:,:),matU2(:,:),matV2(:,:),tmp(:,:),matU1D(:,:),matV1D(:,:),Vin(:,:),Vout1(:,:),Vout2(:,:),Vinter(:,:),Fullmat(:,:),QQ1(:,:)
	real(kind=8),allocatable::Singular(:)
	integer nsproc1,nsproc2,nprow,npcol,nprow1D,npcol1D,myrow,mycol,nprow1,npcol1,myrow1,mycol1,nprow2,npcol2,myrow2,mycol2,myArows,myAcols,M1,N1,M2,N2,rank1,rank2,ierr
	integer::descsmatU(9),descsmatV(9),descsmatU1(9),descsmatV1(9),descsmatU2(9),descsmatV2(9),descUU(9),descVV(9),descsmatU1c(9),descsmatU2c(9),descsmatV1c(9),descsmatV2c(9),descButterflyV(9),descButterflyU(9),descButterU1D(9),descButterV1D(9),descVin(9),descVout(9),descVinter(9),descFull(9)
	integer dims(6),dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
	DT:: TEMP(1)
	integer LWORK,mnmax,mnmin,rank_new
	DT,allocatable:: WORK(:)	
	real(kind=8),allocatable::RWORK(:),center(:)
	real(kind=8):: rtemp,dist,error,rtemp1,rtemp0,fnorm1,fnorm0,flop
	integer :: numroc   ! blacs routine
	integer :: nb1Dc, nb1Dr, ctxt1D,frow,Dimn,edge_n,edge_m,MyID,Ntest,rmaxc,rmaxr
	integer,allocatable::M_p(:,:),N_p(:,:)
	type(Hstat)::stats
	
	rank=0
	
	if(.not. (min(blocks%M,blocks%N)>leafsize .or. (associated(gd%gdc))))then ! reach bottom level, call sequential aca

		if(option%RecLR_leaf==RRQR)then
			! !!!!! RRQR
			
			! rmaxc = min(blocks%M,blocks%N)
			! rmaxr = min(blocks%M,blocks%N)
			rmaxc = blocks%N
			rmaxr = blocks%M
			
			call SeudoSkeleton_CompressionForward(blocks,blocks%headm,blocks%headn,blocks%M,blocks%N,rmaxc,rmaxr,rank,option%tol_comp,option%tol_comp,msh,ker,stats,element_Zmn,ptree,gd%ctxt)
			
		else if(option%RecLR_leaf==ACA)then 
			!!!!! ACA-SVD
			rmax = min(3000,min(blocks%M,blocks%N))
			allocate(UU(blocks%M,rmax))
			allocate(VV(rmax,blocks%N))
			allocate(Singular(rmax))
			frow = 1
		   ! !!!!!!!!!!!! picking a good first row may requires some geometry information. The following can be commented out if geometry info is not available	
		   ! Dimn = size(msh%xyz,1)
		   ! dist = 1D300
		   ! allocate(center(Dimn))
		   ! center = 0
		   ! header_n = blocks%headn
		   ! do j=1,blocks%N
			  ! edge_n = header_n-1+j
			  ! center = center + msh%xyz(1:Dimn,msh%info_unk(0,edge_n))
		   ! enddo
		   ! center = center/blocks%N
		   
		   ! header_m = blocks%headm
		   ! do i=1,blocks%M
				! edge_m=header_m-1+i
				! rtemp = sum((msh%xyz(1:Dimn,msh%info_unk(0,edge_m))-center(1:Dimn))**2d0)
				! if(rtemp<dist)then
					! dist = rtemp
					! frow = i
				! endif
		   ! enddo		
			! deallocate(center)
			! !!!!!!!!!!!! 
			
			call ACA_CompressionForward(UU,VV,Singular,blocks%headm,blocks%headn,blocks%M,blocks%N,frow,rmax,rank,option%tol_comp,option%tol_comp,msh,ker,stats,element_Zmn,ptree,error)	

			! if(error>option%tol_comp)then
				! write(*,*)'niam',error
				! deallocate (UU,VV,Singular)	
				! goto 100
			! endif
			
			blocks%rankmax = rank
			blocks%rankmin = rank
			
			allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N,rank));blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%N
					blocks%ButterflyV%blocks(1)%matrix(i,j)=VV(j,i)
				enddo
			enddo
			!$omp end parallel do	
			
			allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M,rank));blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyU%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%M
					blocks%ButterflyU%blocks(1)%matrix(i,j)=UU(i,j)*Singular(j)
				enddo
			enddo
			!$omp end parallel do							
			deallocate (UU,VV,Singular)	
			
		else if(option%RecLR_leaf==BACA)then	
			!!! blocked ACA
			
			rmax = min(5000,min(blocks%M,blocks%N))
			allocate(UU(blocks%M,rmax))
			allocate(VV(rmax,blocks%N))
		
			call BatchACA_CompressionForward(UU,VV,blocks%headm,blocks%headn,blocks%M,blocks%N,rmax,rank,option%tol_comp,option%tol_comp,option%BACA_Batch,msh,ker,stats,element_Zmn,ptree,error)	
					
			blocks%rankmax = rank
			blocks%rankmin = rank
			
			allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N,rank));blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%N
					blocks%ButterflyV%blocks(1)%matrix(i,j)=VV(j,i)
				enddo
			enddo
			!$omp end parallel do	
			
			allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M,rank));blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyU%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%M
					blocks%ButterflyU%blocks(1)%matrix(i,j)=UU(i,j)
				enddo
			enddo
			!$omp end parallel do							
			deallocate (UU,VV)	

		else if(option%RecLR_leaf==SVD)then
			!!!!! SVD
			mn=min(blocks%M,blocks%N)
			allocate (UU(blocks%M,mn),VV(mn,blocks%N),Singular(mn))		
			allocate(QQ1(blocks%M,blocks%N))
			do ii=1,blocks%M
				do jj =1,blocks%N
					edge_m = blocks%headm +ii - 1
					edge_n = blocks%headn +jj - 1
					call element_Zmn(edge_m,edge_n,QQ1(ii,jj),msh,ker)
				end do
			end do
			call SVD_Truncate(QQ1,blocks%M,blocks%N,mn,UU,VV,Singular,option%tol_comp,rank,flop=flop)		
			stats%Flop_Fill = stats%Flop_Fill + flop
			blocks%rankmax = rank
			blocks%rankmin = rank
			
			allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N,rank));blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%N
					blocks%ButterflyV%blocks(1)%matrix(i,j)=VV(j,i)
				enddo
			enddo
			!$omp end parallel do	
			
			allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M,rank));blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyU%blocks(1)%ndim=rank
			
			!$omp parallel do default(shared) private(i,j)
			do j=1, rank
				do i=1, blocks%M
					blocks%ButterflyU%blocks(1)%matrix(i,j)=UU(i,j)*Singular(j)
				enddo
			enddo
			!$omp end parallel do				
			
			deallocate(QQ1,UU,VV,Singular)
				
		endif
		
		! !!!!!!! check error
		if(cridx==0)then
		call CheckLRError(blocks,option,msh,ker,stats,element_Zmn,ptree,pgno,gd)
		endif
		! !!!!!!! check error
	else
100		if(.not. associated(gd%gdc))then
			gdc1=>gd
			gdc2=>gd
		else
			gdc1=>gd%gdc(1)
			gdc2=>gd%gdc(2)
		endif
	
		call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
			nsproc1 = gdc1%nsprow*gdc1%nspcol
			nsproc2 = gdc2%nsprow*gdc2%nspcol
			
			dims_tmp(1:3)=0
			call blacs_gridinfo(gdc1%ctxt, nprow1, npcol1, myrow1, mycol1)
			
			! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow1,npcol1,myrow1,mycol1,'dddd'
			if(nprow1/=-1 .and. npcol1/=-1)then
				
				! proportional mapping along row or column dimensions
				if(mod(cridx+1,2)==1)then  ! split along column dimension
					blockc(1)%headm = blocks%headm
					blockc(1)%M = blocks%M
					blockc(1)%headn = blocks%headn
					blockc(1)%N = INT(blocks%N*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					call assert(blockc(1)%N>0 .and. blockc(1)%N<blocks%N,'column of blockc(1) or blockc(2) cannot be empty')
				else  ! split along row dimension
					blockc(1)%headn = blocks%headn
					blockc(1)%N = blocks%N
					blockc(1)%headm = blocks%headm
					blockc(1)%M = INT(blocks%M*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					call assert(blockc(1)%M>0 .and. blockc(1)%M<blocks%M,'row of blockc(1) or blockc(2) cannot be empty')			
				endif
				! write(*,*)blockc(1)%M,blockc(1)%N,'ha1',nsproc1,nsproc2
				allocate (blockc(1)%ButterflyU%blocks(1))
				allocate (blockc(1)%ButterflyV%blocks(1))

				call BlockLR(blockc(1),leafsize,rank,option,msh,ker,stats,element_Zmn,ptree,pgno,gdc1,cridx+1)
				dims_tmp(1)=blockc(1)%M
				dims_tmp(2)=blockc(1)%N
				dims_tmp(3)=blockc(1)%rankmax
			endif
			
			dims_tmp(4:6)=0
			call blacs_gridinfo(gdc2%ctxt, nprow2, npcol2, myrow2, mycol2)
			! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow2,npcol2,myrow2,mycol2,'dddd2'
			if(nprow2/=-1 .and. npcol2/=-1)then
				
				! proportional mapping along row or column dimensions
				if(mod(cridx+1,2)==1)then  ! split along column dimension
					blockc(2)%headm = blocks%headm
					blockc(2)%M = blocks%M
					blockc(2)%headn = blocks%headn + INT(blocks%N*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					blockc(2)%N = blocks%N - INT(blocks%N*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					call assert(blockc(2)%N>0 .and. blockc(2)%N<blocks%N,'column of blockc(1) or blockc(2) cannot be empty')
				else  ! split along row dimension
					blockc(2)%headn = blocks%headn
					blockc(2)%N = blocks%N
					blockc(2)%headm = blocks%headm + INT(blocks%M*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					blockc(2)%M = blocks%M - INT(blocks%M*dble(nsproc1)/(dble(nsproc1+nsproc2)))
					call assert(blockc(2)%M>0 .and. blockc(2)%M<blocks%M,'row of blockc(1) or blockc(2) cannot be empty')			
				endif
				! write(*,*)blockc(2)%M,blockc(2)%N,'ha2'
				allocate (blockc(2)%ButterflyU%blocks(1))
				allocate (blockc(2)%ButterflyV%blocks(1))							
				call BlockLR(blockc(2),leafsize,rank,option,msh,ker,stats,element_Zmn,ptree,pgno,gdc2,cridx+1)
				dims_tmp(4)=blockc(2)%M
				dims_tmp(5)=blockc(2)%N
				dims_tmp(6)=blockc(2)%rankmax
			endif		
			! write(*,*)ptree%MyID,cridx+1
			call MPI_ALLREDUCE(dims_tmp,dims,6,MPI_INTEGER,&
			 MPI_MAX,gd%Comm,ierr)

			M1=dims(1)
			N1=dims(2)
			rank1=dims(3)
			M2=dims(4)
			N2=dims(5)
			rank2=dims(6)
			
			if(mod(cridx+1,2)==1)then  ! merge along column dimension
				
				call assert(M1==M2,'M1/=M2 in column merge')
				
				myArows = numroc(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1+rank2, nbslpk, mycol, 0, npcol)	
				allocate(matU(myArows,myAcols))
				call descinit( descsmatU, M1, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descsmatU')

				myArows = numroc(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1, nbslpk, mycol, 0, npcol)	
				allocate(matV1(myArows,myAcols))
				call descinit( descsmatV1, N1, rank1, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descsmatV1')
				
				myArows = numroc(N2, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank2, nbslpk, mycol, 0, npcol)	
				allocate(matV2(myArows,myAcols))
				call descinit( descsmatV2, N2, rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )				
				call assert(info==0,'descinit fail for descsmatV2')
				
				
				
				if(nprow1/=-1 .and. npcol1/=-1)then
					! redistribute U1
					myArows = numroc(M1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc(rank1, nbslpk, mycol1, 0, npcol1)	
					call descinit( descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )	
					call assert(info==0,'descinit fail for descsmatU1c')

					call pgemr2df90(M1, rank1, blockc(1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, gd%ctxt)
					
					! redistribute V1
					myArows = numroc(N1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc(rank1, nbslpk, mycol1, 0, npcol1)	
					call descinit( descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV1c')					
					call pgemr2df90(N1, rank1, blockc(1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, gd%ctxt)				
				else
					descsmatU1c(2)=-1					
					call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, gd%ctxt)
					descsmatV1c(2)=-1					
					call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, gd%ctxt)
				endif
				
				if(nprow2/=-1 .and. npcol2/=-1)then
					! redistribute U2
					myArows = numroc(M2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc(rank2, nbslpk, mycol2, 0, npcol2)	
					call descinit( descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )	
					call assert(info==0,'descinit fail for descsmatU2c')					
					call pgemr2df90(M2, rank2, blockc(2)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU, 1, 1+rank1, descsmatU, gd%ctxt)
					
					! redistribute V2
					myArows = numroc(N2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc(rank2, nbslpk, mycol2, 0, npcol2)	
					call descinit( descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV2c')					
					call pgemr2df90(N2, rank2, blockc(2)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, gd%ctxt)				
				else
					descsmatU2c(2)=-1
					call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU, 1, 1+rank1, descsmatU, gd%ctxt)
					descsmatV2c(2)=-1
					call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, gd%ctxt)
				endif			
				
				! compute truncated SVD on matU
				mnmax=max(M1,rank1+rank2)
				mnmin=min(M1,rank1+rank2)
				
				myArows = numroc(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(mnmin, nbslpk, mycol, 0, npcol)		
				allocate(UU(myArows,myAcols))
				call descinit( descUU, M1, mnmin, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descUU')
				UU=0			
				myArows = numroc(mnmin, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1+rank2, nbslpk, mycol, 0, npcol)		
				allocate(VV(myArows,myAcols))
				call descinit( descVV, mnmin, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descVV')
				VV=0				
				
				allocate(Singular(mnmin))
				Singular=0

				call pgesvdf90('V', 'V', M1, rank1+rank2, matU, 1, 1, descsmatU, Singular, UU, 1, 1, descUU, VV, 1, 1, descVV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
			
				
				deallocate (matU)
				
				rank_new = mnmin
				rank = rank_new
				if(Singular(1)>SafeUnderflow)then	
					rank_new = mnmin
					do i=1,mnmin
						if (Singular(i)/Singular(1)<=option%tol_comp) then
							rank_new=i
							if(Singular(i)<Singular(1)*option%tol_comp/10)rank_new = i -1
							exit
						end if
					end do	
					rank=rank_new
					
					do ii=1,rank_new
						call g2l(ii,rank_new,nprow,nbslpk,iproc,myi)
						if(iproc==myrow)then
							VV(myi,:) = VV(myi,:)*Singular(ii) 		
						endif
					enddo			
				else
					rank_new=1
					rank=1
					Singular(1)=0
					VV=0
					UU=0
				endif
				
				! compute butterfly U and V	
				blocks%rankmax = rank
				blocks%rankmin = rank
				
				myArows = numroc(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
				allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyV%blocks(1)%matrix=0
				blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
				call descinit( descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descButterflyV')		
				
				myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
				allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
				call descinit( descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descButterflyU')			
				blocks%ButterflyU%blocks(1)%matrix=UU(1:myArows,1:myAcols)
				
				blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyU%blocks(1)%ndim=rank		

				call pgemmf90('N','T',N1,rank,rank1,cone, matV1,1,1,descsmatV1,VV,1,1,descVV,czero,blocks%ButterflyV%blocks(1)%matrix,1,1,descButterflyV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
				call pgemmf90('N','T',N2,rank,rank2,cone, matV2,1,1,descsmatV2,VV,1,1+rank1,descVV,czero,blocks%ButterflyV%blocks(1)%matrix,1+N1,1,descButterflyV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)			
				deallocate(UU,VV,Singular,matV1,matV2)
		
			else
				call assert(N1==N2,'N1/=N2 in row merge')
				myArows = numroc(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1+rank2, nbslpk, mycol, 0, npcol)	
				allocate(matV(myArows,myAcols))
				call descinit( descsmatV, N1, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descsmatV')

				myArows = numroc(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1, nbslpk, mycol, 0, npcol)	
				allocate(matU1(myArows,myAcols))
				call descinit( descsmatU1, M1, rank1, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descsmatU1')
				
				myArows = numroc(M2, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank2, nbslpk, mycol, 0, npcol)	
				allocate(matU2(myArows,myAcols))
				call descinit( descsmatU2, M2, rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )				
				call assert(info==0,'descinit fail for descsmatU2')
				
				
				
				if(nprow1/=-1 .and. npcol1/=-1)then
					! redistribute U1
					myArows = numroc(M1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc(rank1, nbslpk, mycol1, 0, npcol1)	
					call descinit( descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )	
					call assert(info==0,'descinit fail for descsmatU1c')				
					! write(*,*)shape(blockc(1)%ButterflyU%blocks(1)%matrix),shape(matU1),rank1,M1,blockc(1)%rankmax
					call pgemr2df90(M1, rank1, blockc(1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, gd%ctxt)
					
					! redistribute V1
					myArows = numroc(N1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc(rank1, nbslpk, mycol1, 0, npcol1)	
					call descinit( descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV1c')						
					call pgemr2df90(N1, rank1, blockc(1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, gd%ctxt)				
				else
					descsmatU1c(2)=-1
					call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, gd%ctxt)
					descsmatV1c(2)=-1
					call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, gd%ctxt)
				endif
				
				if(nprow2/=-1 .and. npcol2/=-1)then
					! redistribute U2
					myArows = numroc(M2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc(rank2, nbslpk, mycol2, 0, npcol2)	
					call descinit( descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )	
					call assert(info==0,'descinit fail for descsmatU2c')				
					call pgemr2df90(M2, rank2, blockc(2)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, gd%ctxt)
					
					! redistribute V2
					myArows = numroc(N2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc(rank2, nbslpk, mycol2, 0, npcol2)	
					call descinit( descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV2c')		
					call pgemr2df90(N2, rank2, blockc(2)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV, 1, 1+rank1, descsmatV, gd%ctxt)				
				else
					descsmatU2c(2)=-1
					call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, gd%ctxt)
					descsmatV2c(2)=-1
					call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV, 1, 1+rank1, descsmatV, gd%ctxt)
				endif			
				
				! compute truncated SVD on matV
				mnmax=max(N1,rank1+rank2)
				mnmin=min(N1,rank1+rank2)
				
				myArows = numroc(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc(mnmin, nbslpk, mycol, 0, npcol)		
				allocate(UU(myArows,myAcols))
				call descinit( descUU, N1, mnmin, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descUU')
				UU=0			
				myArows = numroc(mnmin, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank1+rank2, nbslpk, mycol, 0, npcol)		
				allocate(VV(myArows,myAcols))
				call descinit( descVV, mnmin, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descVV')
				VV=0				
				
				allocate(Singular(mnmin))
				Singular=0	
				call pgesvdf90('V', 'V', N1, rank1+rank2, matV, 1, 1, descsmatV, Singular, UU, 1, 1, descUU, VV, 1, 1, descVV,flop=flop)			
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

				deallocate (matV)
					
				rank_new = mnmin
				rank = rank_new
				if(Singular(1)>SafeUnderflow)then	
					rank_new = mnmin
					do i=1,mnmin
						if (Singular(i)/Singular(1)<=option%tol_comp) then
							rank_new=i
							if(Singular(i)<Singular(1)*option%tol_comp/10)rank_new = i -1
							exit
						end if
					end do	
					rank=rank_new
					
					do ii=1,rank_new
						call g2l(ii,rank_new,nprow,nbslpk,iproc,myi)
						if(iproc==myrow)then
							VV(myi,:) = VV(myi,:)*Singular(ii) 		
						endif
					enddo	

				else
					rank_new=1
					rank=1
					Singular(1)=0
					VV=0
					UU=0
				endif
								

				
				! compute butterfly U and V	
				blocks%rankmax = rank
				blocks%rankmin = rank
				
				myArows = numroc(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
				allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyV%blocks(1)%matrix=UU(1:myArows,1:myAcols)
				blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
				call descinit( descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descButterflyV')		
				
				myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
				allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
				call descinit( descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descButterflyU')			
				blocks%ButterflyU%blocks(1)%matrix=0
				
				blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyU%blocks(1)%ndim=rank		

				call pgemmf90('N','T',M1,rank,rank1,cone, matU1,1,1,descsmatU1,VV,1,1,descVV,czero,blocks%ButterflyU%blocks(1)%matrix,1,1,descButterflyU,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
			
				call pgemmf90('N','T',M2,rank,rank2,cone, matU2,1,1,descsmatU2,VV,1,1+rank1,descVV,czero,blocks%ButterflyU%blocks(1)%matrix,1+M1,1,descButterflyU,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
			
			
				deallocate(UU,VV,Singular,matU1,matU2)		
			endif
		endif

		! write(*,*)ptree%MyID,cridx,rank,'hei',blocks%M,blocks%N
		
		if(cridx==0)then
		
			! !!!!!!! check error
			call CheckLRError(blocks,option,msh,ker,stats,element_Zmn,ptree,pgno,gd)
			! !!!!!!! check error
				
			! distribute UV factor into 1D grid
			ranktmp=rank
			call MPI_ALLREDUCE(ranktmp,rank,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
			
			ctxt1D = ptree%pgrp(pgno)%ctxt1D
		
		
			call blacs_gridinfo(ctxt1D, nprow, npcol, myrow, mycol)
			nb1Dc=rank
			nb1Dr=ceiling_safe(blocks%M/dble(nprow))			
			myArows = numroc(blocks%M, nb1Dr, myrow, 0, nprow)
			myAcols = numroc(rank, nb1Dr, mycol, 0, npcol)	
			allocate(matU1D(myArows,myAcols))
			call descinit(descButterU1D, blocks%M, rank, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
			call assert(info==0,'descinit fail for descButterU1D')
			matU1D=0	
			allocate(M_p(nprow,2))
			do ii=1,nprow
				M_p(ii,1) = (ii-1)*nb1Dr+1
				M_p(ii,2) = ii*nb1Dr
			enddo	
			M_p(nprow,2) = blocks%M
			call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
			if(myrow/=-1 .and. mycol/=-1)then
				call pgemr2df90(blocks%M, rank, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descButterflyU, matU1D, 1, 1, descButterU1D, ctxt1D)
			else 
				descButterflyU(2)=-1
				call pgemr2df90(blocks%M, rank, tmp, 1, 1, descButterflyU, matU1D, 1, 1, descButterU1D, ctxt1D)
			endif
			

			call blacs_gridinfo(ctxt1D, nprow, npcol, myrow, mycol)
			nb1Dc=rank
			nb1Dr=ceiling_safe(blocks%N/dble(nprow))
			myArows = numroc(blocks%N, nb1Dr, myrow, 0, nprow)
			myAcols = numroc(rank, nb1Dr, mycol, 0, npcol)	
			allocate(matV1D(myArows,myAcols))
			call descinit(descButterV1D, blocks%N, rank, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
			call assert(info==0,'descinit fail for descButterV1D')
			matV1D=0	
			allocate(N_p(nprow,2))
			do ii=1,nprow
				N_p(ii,1) = (ii-1)*nb1Dr+1
				N_p(ii,2) = ii*nb1Dr
			enddo	
			N_p(nprow,2) = blocks%N			
			call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
			if(myrow/=-1 .and. mycol/=-1)then			
				call pgemr2df90(blocks%N, rank, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descButterflyV, matV1D, 1, 1, descButterV1D, ctxt1D)	
			else
				descButterflyV(2)=-1
				call pgemr2df90(blocks%N, rank, tmp, 1, 1, descButterflyV, matV1D, 1, 1, descButterV1D, ctxt1D)				
			endif


			
			! redistribute from blacs 1D grid to 1D grid conformal to leaf sizes
			if(allocated(blocks%ButterflyU%blocks(1)%matrix))deallocate(blocks%ButterflyU%blocks(1)%matrix)
			allocate(blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc,rank))
			blocks%ButterflyU%blocks(1)%matrix=0

			if(allocated(blocks%ButterflyV%blocks(1)%matrix))deallocate(blocks%ButterflyV%blocks(1)%matrix)
			allocate(blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc,rank))
			blocks%ButterflyV%blocks(1)%matrix=0	

			call Redistribute1Dto1D(matU1D,M_p,0,pgno,blocks%ButterflyU%blocks(1)%matrix,blocks%M_p,0,pgno,rank,ptree)

			call Redistribute1Dto1D(matV1D,N_p,0,pgno,blocks%ButterflyV%blocks(1)%matrix,blocks%N_p,0,pgno,rank,ptree)

		
			deallocate(matU1D)
			deallocate(matV1D)
			deallocate(N_p)
			deallocate(M_p)
		endif
	
	endif
	

end subroutine BlockLR




subroutine ACA_CompressionForward(matU,matV,Singular,header_m,header_n,rankmax_r,rankmax_c,frow,rmax,rank,tolerance,SVD_tolerance,msh,ker,stats,element_Zmn,ptree,error)
	
	
    use HODLR_DEFS
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real(kind=8) norm_Z,norm_U,norm_V,tolerance,SVD_tolerance,dist
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n,Dimn,mn,Navr,itr
    integer rank, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax,idxs_r,idxs_c,frow
    DT value_Z,maxvalue
    DT inner_U,inner_V,ctemp,value_UVs
    real(kind=8) inner_UV,n1,n2,a,error,flop
    integer,allocatable:: select_column(:), select_row(:)
	DT::matU(rankmax_r,rmax),matV(rmax,rankmax_c)
	real(kind=8)::Singular(rmax)
    DT,allocatable:: row_R(:),column_R(:),value_UV(:)
    real(kind=8),allocatable:: norm_row_R(:),norm_column_R(:),norm_UVavrbynorm_Z(:)
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)	
	real(kind=8), allocatable :: Singularsml(:)
	type(Hstat)::stats
	

	Navr=3 !5 !10
	itr=1
	allocate(norm_UVavrbynorm_Z(Navr))
	norm_UVavrbynorm_Z=0

	
	n1 = OMP_get_wtime()
	
	allocate(select_column(rankmax_c))
	allocate(select_row(rankmax_r))
	allocate(value_UV(max(rankmax_c,rankmax_r)))
	value_UV=0
	
	rankmax_min = min(rankmax_r,rankmax_c)
	norm_Z=0
	select_column = 0
	select_row = 0
	
	allocate(row_R(rankmax_c),column_R(rankmax_r))
	allocate(norm_row_R(rankmax_c),norm_column_R(rankmax_r))

	select_row(1)=frow	

	!$omp parallel do default(shared) private(j,value_Z,edge_m,edge_n)
	do j=1, rankmax_c
		! value_Z=mat(select_row(1),j)
		edge_m = header_m + select_row(1) - 1
		edge_n = header_n + j - 1
		call element_Zmn(edge_m,edge_n,value_Z,msh,ker)							 
		row_R(j)=value_Z
		norm_row_R(j)=dble(value_Z*conjg(cmplx(value_Z,kind = 8)))
	enddo
	!$omp end parallel do
	
	select_column(1)=maxloc(norm_row_R,1)
	maxvalue=row_R(select_column(1))	
	
	
	
	if(abs(maxvalue)<SafeUnderflow)then
	
		do ii=1,100
			call random_number(a)
			select_row(1)=floor_safe(a*(rankmax_r-1))+1
			!$omp parallel do default(shared) private(j,value_Z,edge_m,edge_n)
			do j=1, rankmax_c
				! value_Z=mat(select_row(1),j)
				edge_m = header_m + select_row(1) - 1
				edge_n = header_n + j - 1
				call element_Zmn(edge_m,edge_n,value_Z,msh,ker)							 
				row_R(j)=value_Z
				norm_row_R(j)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
			enddo
			!$omp end parallel do
			
			select_column(1)=maxloc(norm_row_R,1)
			maxvalue=row_R(select_column(1))
			if(abs(maxvalue)>SafeUnderflow)exit
		end do
		if(abs(maxvalue)<SafeUnderflow)then
			rank = 1
			matU(:,1)=0
			matV(1,:)=0
			Singular(1)=0
			
			deallocate(select_column)
			deallocate(select_row)
			deallocate(value_UV)
			deallocate(row_R,column_R)
			deallocate(norm_row_R,norm_column_R)	
			deallocate(norm_UVavrbynorm_Z)	
			return	
		endif
	end if
	

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
		call element_Zmn(edge_m,edge_n,value_Z,msh,ker)
		! value_Z=mat(i,select_column(1))
		column_R(i)=value_Z
		norm_column_R(i)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
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
	if(norm_Z>SafeUnderflow)then
		norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
	else
		norm_UVavrbynorm_Z(itr) = 0
	endif
	itr=mod(itr,Navr)+1
	
	! if(rankmax<2)write(*,*)'rankmax'
	select_row(2)=maxloc(norm_column_R,1)

	rank=1
	! write(*,*)column_R,row_R
	! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
	do while (tolerance**2<(sum(norm_UVavrbynorm_Z)/min(Navr,rank)) .and. rank<rankmax_min)

	
		!$omp parallel do default(shared) private(j,i,value_Z,edge_m,edge_n)
		do j=1,rankmax_c
			edge_m = header_m + select_row(rank+1) - 1
			edge_n = header_n + j - 1
			call element_Zmn(edge_m,edge_n,row_R(j),msh,ker)	
		enddo
		!$omp end parallel do

		call gemmf77('N','N',1,rankmax_c,rank, cone, matU(select_row(rank+1),1), rankmax_r,matV,rmax,czero,value_UV,1)

		! call gemmf90(matU(select_row(rank+1),1), rankmax_r,matV,rmax,value_UV,1,'N','N',1,rankmax_c,rank,cone,czero)	

		
		row_R=row_R-value_UV(1:rankmax_c)
		norm_row_R=dble(row_R*conjg(cmplx(row_R,kind=8)))
		stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c
		
		! write(*,*)'aha',rank
		do i=1,rank
			norm_row_R(select_column(i))=0
		enddo

		select_column(rank+1)=maxloc(norm_row_R,1)
		maxvalue=row_R(select_column(rank+1))
		
		if(abs(maxvalue)<SafeUnderflow)then
			! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
			exit  
			row_R=0
		else 
			do j=1,rankmax_c
			row_R(j)=row_R(j)/maxvalue
			enddo			
		endif
		! !$omp parallel do default(shared) private(j)
		

		
		! !$omp end parallel do
		! !$omp parallel do default(shared) private(j)
		do j=1,rankmax_c
			matV(rank+1,j)=row_R(j)
		enddo
		! !$omp end parallel do

		
		!$omp parallel do default(shared) private(i,j,value_Z,value_UVs,edge_m,edge_n)
		do i=1,rankmax_r
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(rank+1) - 1
			call element_Zmn(edge_m,edge_n,column_R(i),msh,ker)		
		enddo
		!$omp end parallel do		



		call gemmf77('N','N',rankmax_r,1,rank, cone, matU, rankmax_r,matV(1,select_column(rank+1)),rmax,czero,value_UV,rankmax_r)

	! call gemmf90(matU, rankmax_r,matV(1,select_column(rank+1)),rmax,value_UV,rankmax_r,'N','N',rankmax_r,1,rank,cone,czero)			
		
		
		column_R=column_R-value_UV(1:rankmax_r)
		norm_column_R=dble(column_R*conjg(cmplx(column_R,kind=8)))		
		stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_r
		
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
		!$omp parallel do default(shared) private(j,i,ctemp,inner_V,inner_U) reduction(+:inner_UV)
		do j=1,rank
			inner_U=0
			inner_V=0
			! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
			do i=1,rankmax_r
				ctemp=matU(i,rank+1)*conjg(cmplx(matU(i,j),kind=8))
				inner_U=inner_U+ctemp
			enddo
			! !$omp end parallel do
			! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
			do i=1,rankmax_c
				ctemp=matV(rank+1,i)*conjg(cmplx(matV(j,i),kind=8))
				inner_V=inner_V+ctemp
			enddo
			! !$omp end parallel do
			inner_UV=inner_UV+2*dble(inner_U*inner_V)
		enddo
		!$omp end parallel do
		
		norm_Z=norm_Z+inner_UV+norm_U*norm_V
		
		
		! ! write(*,*)norm_Z,inner_UV,norm_U,norm_V,maxvalue,rank,'gan'
		! if(isnan(sqrt(norm_Z)))then
			! write(*,*)inner_UV,norm_U,norm_V,maxvalue
			! stop
		! endif
		
		if(norm_Z>SafeUnderflow)then
			norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
		else
			norm_UVavrbynorm_Z(itr) = 0
		endif
		itr=mod(itr,Navr)+1
		
		stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c + rank*rankmax_r

		rank=rank+1
		if(rank>rmax)then
			write(*,*)'increase rmax',rank,rmax
			stop
		end if
		if (rank<rankmax_min) then
			select_row(rank+1)=maxloc(norm_column_R,1)
		endif

		if(norm_Z<0)exit
		
		! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)		
		
	enddo
	! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)
	
	error = sqrt((sum(norm_UVavrbynorm_Z)/Navr))
	
	! write(*,*)select_row(1:rank),select_column(1:rank)
	
	deallocate(row_R,column_R)
	deallocate(norm_row_R,norm_column_R)
	deallocate(norm_UVavrbynorm_Z)	

	n2 = OMP_get_wtime()	
	time_tmp = time_tmp + n2 - n1	
	
! ACA followed by SVD	
	
	allocate(QQ1(rankmax_r,rank))
	QQ1 = matU(1:rankmax_r,1:rank)
	! call copymatN(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
	allocate (tau_Q(rank))
	call geqrff90(QQ1,tau_Q,flop=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	
	
	allocate (RR1(rank,rank))
	RR1=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call un_or_gqrf90(QQ1,tau_Q,flop=flop)
	deallocate(tau_Q)
	stats%Flop_Fill = stats%Flop_Fill + flop

	allocate(QQ2(rankmax_c,rank))
	call copymatT(matV(1:rank,1:rankmax_c),QQ2,rank,rankmax_c)
	allocate (tau_Q(rank))
	call geqrff90(QQ2,tau_Q,flop=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	
	allocate (RR2(rank,rank))
	RR2=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR2(i,j)=QQ2(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call un_or_gqrf90(QQ2,tau_Q,flop=flop)
	deallocate(tau_Q)
	stats%Flop_Fill = stats%Flop_Fill + flop
	
	allocate(mattemp(rank,rank))
	mattemp=0
	call gemmf90(RR1,rank,RR2,rank,mattemp,rank,'N','T',rank,rank,rank,cone,czero,flop=flop)
	! call zgemm('N','T',rank,rank,rank, cone, RR1, rank,RR2,rank,czero,mattemp,rank)
	stats%Flop_Fill = stats%Flop_Fill + flop
	allocate(UUsml(rank,rank),VVsml(rank,rank),Singularsml(rank))
	call SVD_Truncate(mattemp,rank,rank,rank,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew,flop=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	! call zgemm('N','N',rankmax_r,ranknew,rank, cone, QQ1, rankmax_r,UUsml,rank,czero,matU,rankmax_r)
	call gemmf90(QQ1,rankmax_r,UUsml,rank,matU,rankmax_r,'N','N',rankmax_r,ranknew,rank,cone,czero,flop=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	! call zgemm('N','T',ranknew,rankmax_c,rank, cone, VVsml, rank,QQ2,rankmax_c,czero,matV,rmax) 
	call gemmf90(VVsml,rank,QQ2,rankmax_c,matV,rmax,'N','T',ranknew,rankmax_c,rank,cone,czero,flop=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	
	rank = ranknew
	Singular(1:ranknew) = Singularsml(1:ranknew)
	
	deallocate(mattemp,RR1,QQ1,UUsml,VVsml,Singularsml)
	deallocate(QQ2,RR2)

	deallocate(select_column)
	deallocate(select_row)
	deallocate(value_UV)
	
	
    return

end subroutine ACA_CompressionForward




subroutine BatchACA_CompressionForward(matU,matV,header_m,header_n,M,N,rmax,rank,tolerance,SVD_tolerance,bsize,msh,ker,stats,element_Zmn,ptree,error)
	
	
    use HODLR_DEFS
    implicit none

	integer rank, rankup, ranknew, row, column, rankmax,N,M,rmax
	DT,allocatable:: row_R(:,:),row_Rtmp(:,:),column_R(:,:),column_RT(:,:),fullmat(:,:),fullmat1(:,:)
	DT::matU(M,rmax),matV(rmax,N)
	DT, allocatable :: core(:,:),core_inv(:,:),tau(:),matUtmp(:,:),matVtmp(:,:)
	real(kind=8):: normA,normUV,flop
	integer itr,itrmax,r_est,Nqr,bsize
	integer,allocatable:: select_column(:), select_row(:), perms(:)
	integer,allocatable :: jpvt(:)
	
    integer i, j, ii, jj, indx, rank_1, rank_2
    real(kind=8) tolerance,SVD_tolerance
    integer edge_m, edge_n, header_m, header_n,mn
    real(kind=8) inner_UV,n1,n2,a,error

	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	type(Hstat)::stats
	
	n1 = OMP_get_wtime()
	
	! allocate(fullmat(M,N))
	! !$omp parallel do default(shared) private(i,j,edge_m,edge_n)
	! do i=1,M
	! do j=1,N
		! edge_m = header_m + i - 1
		! edge_n = header_n + j - 1
		! call element_Zmn(edge_m,edge_n,fullmat(i,j),msh,ker)
	! enddo
	! enddo
	! !$omp end parallel do		

	
	r_est=min(bsize,min(M,N))
	! r_est=min(M,N)
	itrmax = floor_safe(min(M,N)/dble(r_est))*2
	
	Nqr = max(max(M,N),r_est)
	allocate(jpvt(Nqr))
	allocate(tau(Nqr))
	
	
	allocate(select_column(r_est))
	allocate(select_row(r_est))
	select_column = 0
	select_row = 0	
	
	allocate(row_R(r_est,N))
	allocate(row_Rtmp(r_est,N))
	allocate(column_R(M,r_est))
	allocate(column_RT(r_est,M))
	row_R=0
	row_Rtmp=0
	column_R=0
	column_RT=0
	allocate(perms(max(M,N)))
	allocate(core(r_est,r_est))
	! allocate(core_inv(r_est,r_est))
	core=0
	! core_inv=0
	
	normA=0
	normUV=0
	itr=0
	rank=0
	
	do while (normUV>=tolerance*normA .and. itr<itrmax)
		
		!**** create random column index for the first iteration
		if(rank==0)then
			call rperm(N, perms)
			select_column = perms(1:r_est)
		endif
		
		!**** Compute columns column_R to find a new set of rows and columns 
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do j=1,r_est
		do i=1,M
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(j) - 1		
			call element_Zmn(edge_m,edge_n,column_R(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do	
		if(rank>0)then
			do j=1,r_est
<<<<<<< HEAD:SRC/Bplus_compress.f90
			

			call gemmf77('N','N',M,1,rank, -cone, matU, M,matV(1,select_column(j)),rmax,cone,column_R(1,j),M)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1,rank)
		
			! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-cone,cone)	
		
=======
			

			call gemmf77('N','N',M,1,rank, -cone, matU, M,matV(1,select_column(j)),rmax,cone,column_R(1,j),M)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1,rank)
		
			! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-cone,cone)	
		
			
			enddo
		endif
		
		
		!**** Find row pivots from the columns column_R
		call copymatT(column_R,column_RT,M,r_est)
		jpvt=0	
		! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(column_RT,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop		
		select_row(1:r_est) = jpvt(1:r_est)
		
	
	
		!**** Compute rows row_R in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do i=1,r_est
		do j=1,N
			edge_m = header_m + select_row(i) - 1
			edge_n = header_n + j - 1			
			call element_Zmn(edge_m,edge_n,row_R(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do		
		if(rank>0)then
			do i=1,r_est



			call gemmf77('N','N',1,N,rank, -cone, matU(select_row(i),1), M,matV,rmax,cone,row_R(i,1),r_est)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N,rank)
		
			! call gemmf90(matU(select_row(i),1), M,matV,rmax,row_R(i,1),r_est,'N','N',1,N,rank,-cone,cone)			
>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90
			
			enddo
		endif
		
<<<<<<< HEAD:SRC/Bplus_compress.f90
					
		!**** Find row pivots from the columns column_R
		call copymatT(column_R,column_RT,M,r_est)
		jpvt=0	
		! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(column_RT,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop		
		select_row(1:r_est) = jpvt(1:r_est)
		

	
		!**** Compute rows row_R in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do i=1,r_est
		do j=1,N
			edge_m = header_m + select_row(i) - 1
			edge_n = header_n + j - 1			
			call element_Zmn(edge_m,edge_n,row_R(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do		
		if(rank>0)then
			do i=1,r_est



			call gemmf77('N','N',1,N,rank, -cone, matU(select_row(i),1), M,matV,rmax,cone,row_R(i,1),r_est)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N,rank)
		
			! call gemmf90(matU(select_row(i),1), M,matV,rmax,row_R(i,1),r_est,'N','N',1,N,rank,-cone,cone)			
			
			enddo
		endif
		
		!**** Find column pivots from the rows row_R
		jpvt=0
		row_Rtmp = row_R		
		! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(row_Rtmp,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop	
		select_column(1:r_est) = jpvt(1:r_est)
		! write(*,*)sum(jpvt)
		
		!**** Compute columns column_R in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do j=1,r_est
		do i=1,M
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(j) - 1
			call element_Zmn(edge_m,edge_n,column_R(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do	
		if(rank>0)then
			do j=1,r_est
			

			call gemmf77('N','N',M,1,rank, -cone, matU, M,matV(1,select_column(j)),rmax,cone,column_R(1,j),M)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1,rank)
	
			! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-cone,cone)

			
			enddo
		endif		
		
		
		!**** Compute the skeleton matrix in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do i=1,r_est
		do j=1,r_est
			edge_m = header_m + select_row(i) - 1
			edge_n = header_n + select_column(j) - 1
			call element_Zmn(edge_m,edge_n,core(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do	
		if(rank>0)then
			allocate(matUtmp(r_est,rank))
			allocate(matVtmp(rank,r_est))
			do i=1,r_est
			matUtmp(i,1:rank) = matU(select_row(i),1:rank)
			enddo
			do j=1,r_est
			matVtmp(1:rank,j) = matV(1:rank,select_column(j))
			enddo		
			! call zgemm('N','N',r_est,r_est,rank, -cone, matUtmp, r_est,matVtmp,rank,cone,core,r_est)
			call gemmf90(matUtmp,r_est,matVtmp,rank,core,r_est,'N','N',r_est,r_est,rank, -cone,cone,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop
			deallocate(matUtmp)
			deallocate(matVtmp)
		endif	
			
			
		!**** generate column indices for the next iteration	
		row_Rtmp = row_R
		do j=1,r_est
		row_Rtmp(:,select_column(j))=0d0
		enddo
		
		
		jpvt=0	
		! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(row_Rtmp,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		select_column(1:r_est) = jpvt(1:r_est)	
=======
		!**** Find column pivots from the rows row_R
		jpvt=0
		row_Rtmp = row_R		
		! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(row_Rtmp,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop	
		select_column(1:r_est) = jpvt(1:r_est)
		! write(*,*)sum(jpvt)
		!**** Compute columns column_R in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do j=1,r_est
		do i=1,M
			edge_m = header_m + i - 1
			edge_n = header_n + select_column(j) - 1
			call element_Zmn(edge_m,edge_n,column_R(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do	
		if(rank>0)then
			do j=1,r_est
			

			call gemmf77('N','N',M,1,rank, -cone, matU, M,matV(1,select_column(j)),rmax,cone,column_R(1,j),M)
			stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1,rank)
	
			! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-cone,cone)

			
			enddo
		endif		
		
		
		!**** Compute the skeleton matrix in CUR
		!$omp parallel do default(shared) private(i,j,edge_m,edge_n)
		do i=1,r_est
		do j=1,r_est
			edge_m = header_m + select_row(i) - 1
			edge_n = header_n + select_column(j) - 1
			call element_Zmn(edge_m,edge_n,core(i,j),msh,ker)
		enddo
		enddo
		!$omp end parallel do	
		if(rank>0)then
			allocate(matUtmp(r_est,rank))
			allocate(matVtmp(rank,r_est))
			do i=1,r_est
			matUtmp(i,1:rank) = matU(select_row(i),1:rank)
			enddo
			do j=1,r_est
			matVtmp(1:rank,j) = matV(1:rank,select_column(j))
			enddo		
			! call zgemm('N','N',r_est,r_est,rank, -cone, matUtmp, r_est,matVtmp,rank,cone,core,r_est)
			call gemmf90(matUtmp,r_est,matVtmp,rank,core,r_est,'N','N',r_est,r_est,rank, -cone,cone,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop
			deallocate(matUtmp)
			deallocate(matVtmp)
		endif	
			
			
		!**** generate column indices for the next iteration	
		row_Rtmp = row_R
		do j=1,r_est
		row_Rtmp(:,select_column(j))=0d0
		enddo
		
		
		jpvt=0	
		! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
		call geqp3f90(row_Rtmp,jpvt,tau,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		select_column(1:r_est) = jpvt(1:r_est)	
		
		
		!**** form the LR update by CUR
		
		
		! call GeneralInverse(r_est,r_est,core,core_inv,tolerance)
		! call zgemm('N','N',r_est,N,r_est, cone, core_inv, r_est,row_R,r_est,czero,row_Rtmp,r_est)
		! ! call LeastSquare(r_est,r_est,N,core,row_R,row_Rtmp,tolerance)
>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90
		
		! call LR_ReCompression(column_R,row_Rtmp,M,N,r_est,rankup,SVD_tolerance)
		! ! rankup = r_est
		
<<<<<<< HEAD:SRC/Bplus_compress.f90
		!**** form the LR update by CUR
=======
		! if(rankup>0)then
>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90
		
		! call assert(rank+rankup<rmax,'try to increase rmax')
		
<<<<<<< HEAD:SRC/Bplus_compress.f90
		! call GeneralInverse(r_est,r_est,core,core_inv,tolerance)
		! call zgemm('N','N',r_est,N,r_est, cone, core_inv, r_est,row_R,r_est,czero,row_Rtmp,r_est)
		! ! call LeastSquare(r_est,r_est,N,core,row_R,row_Rtmp,tolerance)
		
		! call LR_ReCompression(column_R,row_Rtmp,M,N,r_est,rankup,SVD_tolerance)
		! ! rankup = r_est
		
		! if(rankup>0)then
		
		! call assert(rank+rankup<rmax,'try to increase rmax')
		
		! matU(:,rank+1:rank+rankup) = column_R(:,1:rankup)
		! matV(rank+1:rank+rankup,:) = row_Rtmp(1:rankup,:)
		! rank = rank + rankup
		
		! !**** update fnorm of UV and matUmatV 
		! call LR_Fnorm(column_R,row_Rtmp,M,N,rankup,normUV,tolerance*1e-2)		
		! call LR_Fnorm(matU,matV,M,N,rank,normA,tolerance*1e-2)
		! endif
		
		
		jpvt=0	
		call geqp3modf90(core,jpvt,tau,tolerance,SafeUnderflow,ranknew,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		rankup = ranknew
		
		if(rankup>0)then		
			row_Rtmp = row_R
			call un_or_mqrf90(core,tau,row_Rtmp,'L','C',r_est,N,rankup,flop=flop)	
			stats%Flop_Fill = stats%Flop_Fill + flop
			call trsmf90(core,row_Rtmp,'L','U','N','N',rankup,N,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop
			
			if(rank+rankup>min(M,N))rankup=min(M,N)-rank
			
			call assert(rank+rankup<=rmax,'try to increase rmax')
			do j=1,rankup
			matU(:,rank+j) = column_R(:,jpvt(j))
			enddo
			matV(rank+1:rank+rankup,:) = row_Rtmp(1:rankup,:)
			rank = rank + rankup
			
			if(rank==min(M,N))exit
			
			
			!**** update fnorm of UV and matUmatV 
			call LR_Fnorm(column_R,row_Rtmp,M,N,rankup,normUV,tolerance*1e-2,Flops=flop)	
			stats%Flop_Fill = stats%Flop_Fill + flop		
			call LR_Fnorm(matU,matV,M,N,rank,normA,tolerance*1e-2,Flops=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop
			
			if(normA>SafeUnderflow)then
				error = normUV/normA	
			else
				error = 0
				exit			
			endif			
		else 
			error = 0
			exit
		endif
		

		! write(*,*)itr,rank,rankup,error
=======
		! matU(:,rank+1:rank+rankup) = column_R(:,1:rankup)
		! matV(rank+1:rank+rankup,:) = row_Rtmp(1:rankup,:)
		! rank = rank + rankup
		
		! !**** update fnorm of UV and matUmatV 
		! call LR_Fnorm(column_R,row_Rtmp,M,N,rankup,normUV,tolerance*1e-2)		
		! call LR_Fnorm(matU,matV,M,N,rank,normA,tolerance*1e-2)
		! endif
		
		
		jpvt=0	
		call geqp3modf90(core,jpvt,tau,tolerance,SafeUnderflow,ranknew,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		rankup = ranknew
		
		if(rankup>0)then		
		row_Rtmp = row_R
		call un_or_mqrf90(core,tau,row_Rtmp,'L','C',r_est,N,rankup,flop=flop)	
		stats%Flop_Fill = stats%Flop_Fill + flop
		call trsmf90(core,row_Rtmp,'L','U','N','N',rankup,N,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		
		if(rank+rankup>min(M,N))rankup=min(M,N)-rank
		
		call assert(rank+rankup<=rmax,'try to increase rmax')
		do j=1,rankup
		matU(:,rank+j) = column_R(:,jpvt(j))
		enddo
		matV(rank+1:rank+rankup,:) = row_Rtmp(1:rankup,:)
		rank = rank + rankup
		
		if(rank==min(M,N))exit
		
		!**** update fnorm of UV and matUmatV 
		call LR_Fnorm(column_R,row_Rtmp,M,N,rankup,normUV,tolerance*1e-2,Flops=flop)	
		stats%Flop_Fill = stats%Flop_Fill + flop		
		call LR_Fnorm(matU,matV,M,N,rank,normA,tolerance*1e-2,Flops=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		
		endif
		
		
		write(*,*)itr,rank,rankup,normUV/normA
>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90
	
		itr = itr + 1
	enddo

	! write(*,*)normUV,normA
<<<<<<< HEAD:SRC/Bplus_compress.f90
	
	if(rank>0)then	
		call LR_ReCompression(matU,matV,M,N,rank,ranknew,SVD_tolerance,Flops=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		rank = ranknew
	else
		rank = 1
		matU(:,1)=0
		matV(1,:)=0
	endif

	deallocate(jpvt)
	deallocate(tau)
	deallocate(select_column)
	deallocate(select_row)	
	deallocate(row_R)
	deallocate(row_Rtmp)
	deallocate(column_R)
	deallocate(column_RT)
	deallocate(core)
	! deallocate(core_inv)
	deallocate(perms)
	

	n2 = OMP_get_wtime()	
	time_tmp = time_tmp + n2 - n1	
	
	
    return

end subroutine BatchACA_CompressionForward

=======
	error = normUV/normA
	
	call LR_ReCompression(matU,matV,M,N,rank,ranknew,SVD_tolerance,Flops=flop)
	stats%Flop_Fill = stats%Flop_Fill + flop
	rank = ranknew
	

	deallocate(jpvt)
	deallocate(tau)
	deallocate(select_column)
	deallocate(select_row)	
	deallocate(row_R)
	deallocate(row_Rtmp)
	deallocate(column_R)
	deallocate(column_RT)
	deallocate(core)
	! deallocate(core_inv)
	deallocate(perms)
	

	n2 = OMP_get_wtime()	
	time_tmp = time_tmp + n2 - n1	
	
	
    return

end subroutine BatchACA_CompressionForward

>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90





subroutine SeudoSkeleton_CompressionForward(blocks,header_m,header_n,M,N,rmaxc,rmaxr,rank,tolerance,SVD_tolerance,msh,ker,stats,element_Zmn,ptree,ctxt,pgno)
	
	
    use HODLR_DEFS
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2, rank
    integer index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real(kind=8) norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer ranknew, row, column, rankmax,N,M,rankmax_min,rmax,rmaxc,rmaxr,idxs_r,idxs_c,flag0,myAcols,myArows,npcol,nprow,rank_new,nproc
    DT value_Z,value_UV,maxvalue
    DT inner_U,inner_V,ctemp
    real(kind=8) inner_UV,flop
	DT,allocatable::matU(:,:),matV(:,:),UU(:,:),VV(:,:),MatrixSubselection(:,:)
	real(kind=8),allocatable::Singular(:)
    DT,allocatable:: row_R(:),column_R(:)
    real(kind=8),allocatable:: norm_row_R(:),norm_column_R(:)
	type(mesh)::msh
	type(kernelquant)::ker
	type(matrixblock)::blocks
	procedure(Z_elem)::element_Zmn
	type(Hstat)::stats
	
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:),matU2D(:,:),matV2D(:,:),matU1D(:,:),matV1D(:,:)
	real(kind=8), allocatable :: Singularsml(:)
	type(proctree)::ptree
	integer ctxt,ctxt1D,myrow,mycol,iproc,jproc,myi,myj,mnmax
	integer,optional::pgno
	integer LWORK,LRWORK,INFO,ierr
	DT:: TEMP(1)
	real(kind=8),allocatable::RWORK(:)
	integer::descsmatU(9),descsmatV(9),descsmatVQ(9),descsub(9),descUU(9),descVV(9),descButterU2D(9),descButterV2D(9),descButterU1D(9),descButterV1D(9),descQ(9),descR(9)
	integer :: numroc   ! blacs routine
	
	integer nb1Dr,nb1Dc
	integer,allocatable::select_col(:),select_row(:),N_p(:,:),M_p(:,:)
	DT,allocatable:: WORK(:)
	integer, allocatable :: ipiv(:),jpiv(:),JPERM(:)
	DT, allocatable :: tau(:)
	real(kind=8):: RTEMP(1)
	
	rank_new=0
	! ctxt = ptree%pgrp(pgno)%ctxt
	call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
	
	if(nprow*npcol==1)then

		rmax=min(rmaxc,rmaxr)
		allocate(select_col(rmaxc))
		allocate(select_row(rmaxr))
		call linspaceI(1,M,rmaxr,select_row)	
		call linspaceI(1,N,rmaxc,select_col)
		
		allocate(matV(N,rmaxr))
		matV=0
		do ii=1,N
			do jj=1,rmaxr
				edge_m = header_m + select_row(jj) - 1 
				edge_n = header_n + ii - 1 
				call element_Zmn(edge_m,edge_n,matV(ii,jj),msh,ker)					
			enddo
		enddo	

		allocate(MatrixSubselection(rmaxr,rmaxc))
		MatrixSubselection=0		
		
		do ii=1,rmaxr
			do jj=1,rmaxc
				edge_m = header_m + select_row(ii) - 1 
				edge_n = header_n + select_col(jj) - 1 
				call element_Zmn(edge_m,edge_n,MatrixSubselection(ii,jj),msh,ker)					
			enddo
		enddo	

		allocate(jpiv(rmaxc))
		jpiv=0
		allocate(tau(min(rmaxr,rmaxc)))
		tau=0
		call geqp3modf90(MatrixSubselection,jpiv,tau,tolerance,SafeUnderflow,rank_new,flop=flop)
		stats%Flop_Fill = stats%Flop_Fill + flop
		rank = rank_new
		if(rank>0)then
			matV = conjg(cmplx(matV,kind=8))
			call un_or_mqrf90(MatrixSubselection,tau,matV,'R','N',N, rmaxr, rank_new,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop				
			matV = conjg(cmplx(matV,kind=8))
			call trsmf90(MatrixSubselection,matV,'R', 'U', 'T', 'N', N, rank_new,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop	
		else 
			rank=1
			rank_new=1	
			matV=0
		endif
		allocate(blocks%ButterflyU%blocks(1)%matrix(M,rank))
		blocks%ButterflyU%blocks(1)%matrix=0
		do ii=1,M
			do jj=1,rank
				edge_m = header_m + ii - 1 
				edge_n = header_n + select_col(jpiv(jj)) - 1 
				call element_Zmn(edge_m,edge_n,blocks%ButterflyU%blocks(1)%matrix(ii,jj),msh,ker)					
			enddo
		enddo

		allocate(blocks%ButterflyV%blocks(1)%matrix(N,rank))
		blocks%ButterflyV%blocks(1)%matrix = matV(1:N,1:rank)
<<<<<<< HEAD:SRC/Bplus_compress.f90
		
		blocks%rankmax = rank
		blocks%rankmin = rank		
		
=======
		
		blocks%rankmax = rank
		blocks%rankmin = rank		
		
>>>>>>> 37a8bb5076fbc403962a70a6fb2317f74d01c3af:SRC/Bplus_compress.f90
		deallocate(matV)
		deallocate(jpiv)
		deallocate(tau)
		deallocate(select_col)
		deallocate(select_row)
		deallocate(MatrixSubselection)
	else 
		if(myrow/=-1 .and. mycol/=-1)then
			rmax=min(rmaxc,rmaxr)
			allocate(select_col(rmaxc))
			allocate(select_row(rmaxr))
			call linspaceI(1,M,rmaxr,select_row)	
			call linspaceI(1,N,rmaxc,select_col)

			
			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
			! nproc = npcol*nprow
			myArows = numroc(N, nbslpk, myrow, 0, nprow)
			myAcols = numroc(rmaxr, nbslpk, mycol, 0, npcol)	
			allocate(matV(myArows,myAcols))
			call descinit( descsmatV, N, rmaxr, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call assert(info==0,'descinit fail for descsmatV')
			matV=0	
			
			do myi=1,myArows
				call l2g(myi,myrow,N,nprow,nbslpk,ii)
				do myj=1,myAcols
					call l2g(myj,mycol,rmaxr,npcol,nbslpk,jj)
					edge_m = header_m + select_row(jj) - 1 
					edge_n = header_n + ii - 1 
					call element_Zmn(edge_m,edge_n,matV(myi,myj),msh,ker)					
				enddo
			enddo		

			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
			myArows = numroc(rmaxr, nbslpk, myrow, 0, nprow)
			myAcols = numroc(rmaxc, nbslpk, mycol, 0, npcol)
			
			allocate(MatrixSubselection(myArows,myAcols))
			call descinit( descsub, rmaxr, rmaxc, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call assert(info==0,'descinit fail for descsub')
			MatrixSubselection=0
			
			
			do myi=1,myArows
				call l2g(myi,myrow,rmaxr,nprow,nbslpk,ii)
				do myj=1,myAcols
					call l2g(myj,mycol,rmaxc,npcol,nbslpk,jj)
					edge_m = header_m + select_row(ii) - 1 
					edge_n = header_n + select_col(jj) - 1 
					call element_Zmn(edge_m,edge_n,MatrixSubselection(myi,myj),msh,ker)					
				enddo
			enddo	


			! Compute QR of MatrixSubselection*P
			allocate(ipiv(myAcols))
			ipiv=0
			allocate(tau(myAcols))
			tau=0
			allocate(jpiv(rmaxc))
			jpiv=0
			allocate(JPERM(rmaxc))
			JPERM=0
			call pgeqpfmodf90(rmaxr, rmaxc, MatrixSubselection, 1, 1, descsub, ipiv, tau, JPERM, jpiv, rank_new,tolerance, SafeUnderflow,flop=flop)
			stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

			
			rank = rank_new
			
			if(rank>0)then
			
				! Compute matV*conjg(Q)
				matV = conjg(cmplx(matV,kind=8))
				call pun_or_mqrf90('R', 'N', N, rmaxr, rank_new, MatrixSubselection, 1, 1, descsub, tau, matV, 1, 1, descsmatV)
				matV = conjg(cmplx(matV,kind=8))
				
				
				! Compute matV*conjg(Q)*(R^T)^-1
				call ptrsmf90('R', 'U', 'T', 'N', N, rank_new, cone, MatrixSubselection, 1, 1, descsub, matV, 1, 1, descsmatV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

				call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
				myArows = numroc(N, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank_new, nbslpk, mycol, 0, npcol)	
				allocate(matV2D(myArows,myAcols))
				call descinit( descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterV2D')
				matV2D(1:myArows,1:myAcols) = matV(1:myArows,1:myAcols)
				
			else 
				rank=1
				rank_new=1
				
				call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
				myArows = numroc(N, nbslpk, myrow, 0, nprow)
				myAcols = numroc(rank_new, nbslpk, mycol, 0, npcol)	
				allocate(matV2D(myArows,myAcols))
				call descinit( descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterV2D')
				matV2D(1:myArows,1:myAcols) = 0		
			endif
			
			
			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
			myArows = numroc(M, nbslpk, myrow, 0, nprow)
			myAcols = numroc(rank_new, nbslpk, mycol, 0, npcol)	
			allocate(matU2D(myArows,myAcols))
			call descinit( descButterU2D, M, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call assert(info==0,'descinit fail for descButterU2D')
			matU2D=0
			
			
			do myi=1,myArows
				call l2g(myi,myrow,M,nprow,nbslpk,ii)
				do myj=1,myAcols
					call l2g(myj,mycol,rank_new,npcol,nbslpk,jj)
					edge_m = header_m + ii - 1 
					edge_n = header_n + select_col(ipiv(myj)) - 1 
					call element_Zmn(edge_m,edge_n,matU2D(myi,myj),msh,ker)					
				enddo
			enddo		
		
			deallocate(select_col)
			deallocate(select_row)		
			deallocate(matV)
			deallocate(MatrixSubselection)	
			deallocate(ipiv)
			deallocate(jpiv)
			deallocate(JPERM)
			deallocate(tau)
		endif

		
		! call MPI_ALLREDUCE(rank_new,rank,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)


		if(.not. present(pgno))then
			myArows = numroc(M, nbslpk, myrow, 0, nprow)
			myAcols = numroc(rank, nbslpk, mycol, 0, npcol)	
			allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
			blocks%ButterflyU%blocks(1)%matrix(1:myArows,1:myAcols)=matU2D(1:myArows,1:myAcols)
			deallocate(matU2D)		

			myArows = numroc(N, nbslpk, myrow, 0, nprow)
			myAcols = numroc(rank, nbslpk, mycol, 0, npcol)	
			allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
			blocks%ButterflyV%blocks(1)%matrix(1:myArows,1:myAcols)=matV2D(1:myArows,1:myAcols)
			deallocate(matV2D)		

			blocks%rankmax = rank
			blocks%rankmin = rank
			
		else 
		
			! distribute UV factor into 1D grid

			nproc = ptree%pgrp(pgno)%nproc	
			ctxt1D = ptree%pgrp(pgno)%ctxt1D
			nb1Dc=rank
			nb1Dr=ceiling_safe(M/dble(nproc))
			allocate(M_p(nproc,2))
			do ii=1,nproc
				M_p(ii,1) = (ii-1)*nb1Dr+1
				M_p(ii,2) = ii*nb1Dr
			enddo	
			M_p(nproc,2) = M
			call blacs_gridinfo(ctxt1D, nprow, npcol, myrow, mycol)	
			myArows = numroc(M, nb1Dr, myrow, 0, nprow)
			myAcols = numroc(rank, nb1Dr, mycol, 0, npcol)	
			allocate(matU1D(myArows,myAcols))
			call descinit(descButterU1D, M, rank, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
			call assert(info==0,'descinit fail for descButterU1D')
			matU1D=0	

			nb1Dc=rank
			nb1Dr=ceiling_safe(N/dble(nproc))
			allocate(N_p(nproc,2))
			do ii=1,nproc
				N_p(ii,1) = (ii-1)*nb1Dr+1
				N_p(ii,2) = ii*nb1Dr
			enddo	
			N_p(nproc,2) = N	
			myArows = numroc(N, nb1Dr, myrow, 0, nprow)
			myAcols = numroc(rank, nb1Dr, mycol, 0, npcol)	
			allocate(matV1D(myArows,myAcols))
			call descinit(descButterV1D, N, rank, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
			call assert(info==0,'descinit fail for descButterV1D')
			matV1D=0	
				
			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)			

			if(myrow/=-1 .and. mycol/=-1)then
				call pgemr2df90(M, rank, matU2D, 1, 1, descButterU2D, matU1D, 1, 1, descButterU1D, ctxt1D)		
				call pgemr2df90(N, rank, matV2D, 1, 1, descButterV2D, matV1D, 1, 1, descButterV1D, ctxt1D)	
				deallocate(matU2D)
				deallocate(matV2D)
			else 
				! write(*,*)rank,'nima',ptree%MyID
				descButterU2D(2)=-1
				call pgemr2df90(M, rank, matU2D, 1, 1, descButterU2D, matU1D, 1, 1, descButterU1D, ctxt1D)
				descButterV2D(2)=-1
				call pgemr2df90(N, rank, matV2D, 1, 1, descButterV2D, matV1D, 1, 1, descButterV1D, ctxt1D)		
			endif	

			! redistribute from blacs 1D grid to 1D grid conformal to leaf sizes 
			allocate(blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc,rank))
			blocks%ButterflyU%blocks(1)%mdim=M;blocks%ButterflyU%blocks(1)%ndim=rank
			blocks%ButterflyU%blocks(1)%matrix=0

			allocate(blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc,rank))
			blocks%ButterflyV%blocks(1)%mdim=N;blocks%ButterflyV%blocks(1)%ndim=rank
			blocks%ButterflyV%blocks(1)%matrix=0	
			
			blocks%rankmax = max(blocks%rankmax,rank)
			blocks%rankmin = min(blocks%rankmin,rank)
				

			call Redistribute1Dto1D(matU1D,M_p,0,pgno,blocks%ButterflyU%blocks(1)%matrix,blocks%M_p,0,pgno,rank,ptree)

			call Redistribute1Dto1D(matV1D,N_p,0,pgno,blocks%ButterflyV%blocks(1)%matrix,blocks%N_p,0,pgno,rank,ptree)

			deallocate(N_p)
			deallocate(M_p)
			deallocate(matU1D)
			deallocate(matV1D)
		end if
	endif
	
	blocks%ButterflyV%blocks(1)%mdim=blocks%N;blocks%ButterflyV%blocks(1)%ndim=rank
	blocks%ButterflyU%blocks(1)%mdim=blocks%M;blocks%ButterflyV%blocks(1)%ndim=rank	
	
    return

end subroutine SeudoSkeleton_CompressionForward


subroutine LocalButterflySVD_Left(index_i_loc,index_j_loc,level_loc,level_butterflyL,level,index_i_m,blocks,option,ButterflyP_old,ButterflyP)
use misc
implicit none 
integer index_i_loc,index_j_loc,level_loc,level_butterflyL,index_i_m,index_i,index_j,level,group_m,mm,nn,nn1,nn2,j,i,mn,rank,mm1
type(butterfly_Kerl)ButterflyP_old,ButterflyP
DT, allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
DT, allocatable :: tau_Q(:)
real(kind=8),allocatable :: Singular(:)
type(matrixblock)::blocks
real(kind=8):: SVD_tolerance
type(Hoption)::option

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

	! write(*,*)'dddd',fnorm(QQ,mm,nn)

	mn=min(mm,nn)
	allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
	call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			
	! rank = min(rank,37)

	! rank = 7
	! write(*,*)'dddd', rank
	
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
		allocate(ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix(mm1,rank));ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%mdim=mm1;ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%ndim=rank
		allocate(ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix(mm-mm1,rank));ButterflyP%blocks(2*index_i_loc,index_j_loc)%mdim=mm-mm1;ButterflyP%blocks(2*index_i_loc,index_j_loc)%ndim=rank
		
		ButterflyP%blocks(2*index_i_loc-1,index_j_loc)%matrix=mat_tmp(1:mm1,1:rank)
		ButterflyP%blocks(2*index_i_loc,index_j_loc)%matrix=mat_tmp(1+mm1:mm,1:rank)						
	else 
		allocate (blocks%ButterflyU%blocks(index_i)%matrix(mm,rank));blocks%ButterflyU%blocks(index_i)%mdim=mm;blocks%ButterflyU%blocks(index_i)%ndim=rank
		! !$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, mm
				blocks%ButterflyU%blocks(index_i)%matrix(i,j)=mat_tmp(i,j)
				! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
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
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do i=1, rank
		do j=1, nn2
			blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=VV(i,j+nn1)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
		enddo
	enddo
	! !$omp end parallel do						

	deallocate(QQ,UU,VV,Singular,mat_tmp)

	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
	! if (level_loc==level_butterflyL) then
		! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_i)%matrix)/1024.0d3
	! endif


end subroutine  LocalButterflySVD_Left







subroutine LocalButterflySVD_Right(index_i_loc,index_j_loc,level_loc,level_butterflyR,level,level_butterfly,index_j_m,blocks,option,ButterflyP_old,ButterflyP)
use misc
implicit none 
integer index_i_loc,index_j_loc,level_loc,level_butterflyR,level_butterfly,index_j_m,index_i,index_j,level,group_n,mm,nn,nn1,mm1,mm2,j,i,mn,rank
type(butterfly_Kerl)ButterflyP_old,ButterflyP
DT, allocatable :: QQ(:,:), RR(:,:), UU(:,:), VV(:,:),mat_tmp(:,:),matU(:,:),matV(:,:)
DT, allocatable :: tau_Q(:)
real(kind=8),allocatable :: Singular(:)
type(matrixblock)::blocks
real(kind=8):: SVD_tolerance
type(Hoption)::option

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
	call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,rank)			
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
		allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix(rank,nn1));ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%mdim=rank;ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%ndim=nn1
		allocate(ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix(rank,nn-nn1));ButterflyP%blocks(index_i_loc,2*index_j_loc)%mdim=rank;ButterflyP%blocks(index_i_loc,2*index_j_loc)%ndim=nn-nn1
		
		ButterflyP%blocks(index_i_loc,2*index_j_loc-1)%matrix=mat_tmp(1:rank,1:nn1)
		ButterflyP%blocks(index_i_loc,2*index_j_loc)%matrix=mat_tmp(1:rank,1+nn1:nn)						
	else 
		allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn,rank));blocks%ButterflyV%blocks(index_j)%mdim=nn;blocks%ButterflyV%blocks(index_j)%ndim=rank
		! !$omp parallel do default(shared) private(i,j)
		do j=1, rank
			do i=1, nn
				blocks%ButterflyV%blocks(index_j)%matrix(i,j)=mat_tmp(j,i)
				! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
			enddo
		enddo
		! !$omp end parallel do

	end if		

	allocate (blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(mm1,rank)); blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%mdim=mm1; blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%ndim=rank
	allocate (blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(mm2,rank)); blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%mdim=mm2; blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%ndim=rank	
	
	! !$omp parallel do default(shared) private(i,j)
	do i=1, mm1
		do j=1, rank
			blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix(i,j)=UU(i,j)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
		enddo
	enddo
	! !$omp end parallel do						
	! !$omp parallel do default(shared) private(i,j)
	do i=1, mm2
		do j=1, rank
			blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix(i,j)=UU(i+mm1,j)
			! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
		enddo
	enddo
	! !$omp end parallel do						
	
	deallocate(QQ,UU,VV,Singular,mat_tmp)
  
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
	! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
	! if (level_loc==level_butterflyR) then
		! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
	! endif

end subroutine LocalButterflySVD_Right


!!!!!!! check error of LR compression of blocks by comparing the full block, assuming 2D block-cyclic distribution   
subroutine CheckLRError(blocks,option,msh,ker,stats,element_Zmn,ptree,pgno,gd)
use HODLR_DEFS
implicit none 

	type(matrixblock)::blocks
	type(mesh)::msh
	type(Hoption)::option
	type(kernelquant)::ker
	type(Hstat)::stats
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	integer pgno
	type(grid),pointer::gd
	integer Ntest
	integer nsproc1,nsproc2,nprow,npcol,nprow1D,npcol1D,myrow,mycol,nprow1,npcol1,myrow1,mycol1,nprow2,npcol2,myrow2,mycol2,myArows,myAcols,M1,N1,M2,N2,rank1,rank2,ierr,MyID	
	integer:: cridx,info
	DT,allocatable:: Vin(:,:),Vout1(:,:),Vout2(:,:),Vinter(:,:),Fullmat(:,:)
	integer::descVin(9),descVout(9),descVinter(9),descFull(9),descButterflyU(9),descButterflyV(9)	
	integer N,M,i,j,ii,jj,myi,myj,iproc,jproc,rmax
	integer :: numroc   ! blacs routine	
	integer edge_n,edge_m, rank
	real(kind=8):: fnorm1,fnorm0,rtemp1,rtemp0
	
	
	if(option%ErrFillFull==1)then
	
		Ntest=32
		! Ntest=blocks%N/2
		call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
		rank = blocks%rankmax

		
		myArows = numroc(blocks%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
		call descinit( descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
		call assert(info==0,'descinit fail for descButterflyV')		
		
		myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc(rank, nbslpk, mycol, 0, npcol)		
		call descinit( descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )	
		call assert(info==0,'descinit fail for descButterflyU')	

		myArows = numroc(blocks%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc(Ntest, nbslpk, mycol, 0, npcol)	
		call descinit( descVin, blocks%N, Ntest, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
		allocate(Vin(myArows,myAcols))
		do ii=1,myArows
		do jj=1,myAcols
			call random_dp_number(Vin(ii,jj))
		end do
		end do
		
		myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc(blocks%N, nbslpk, mycol, 0, npcol)			
		call descinit( descFull, blocks%M, blocks%N, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
		allocate(Fullmat(myArows,myAcols))
		do myi=1,myArows
			call l2g(myi,myrow,blocks%M,nprow,nbslpk,ii)
			do myj=1,myAcols
				call l2g(myj,mycol,blocks%N,npcol,nbslpk,jj)
				edge_m = blocks%headm + ii - 1 
				edge_n = blocks%headn + jj - 1 
				call element_Zmn(edge_m,edge_n,Fullmat(myi,myj),msh,ker)					
			enddo
		enddo			

		! compute the exact results
		myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc(Ntest, nbslpk, mycol, 0, npcol)			
		call descinit( descVout, blocks%M, Ntest, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )			
		allocate(Vout1(myArows,myAcols))
		allocate(Vout2(myArows,myAcols))
		Vout1=0
		Vout2=0
		call pgemmf90('N','N',blocks%M,Ntest,blocks%N,cone, Fullmat,1,1,descFull,Vin,1,1,descVin,czero,Vout1,1,1,descVout)
		
		! compute the approximate results
		myArows = numroc(rank, nbslpk, myrow, 0, nprow)
		myAcols = numroc(Ntest, nbslpk, mycol, 0, npcol)			
		call descinit( descVinter, rank, Ntest, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )			
		allocate(Vinter(myArows,myAcols))			
		Vinter=0
		call pgemmf90('T','N',rank,Ntest,blocks%N,cone, blocks%ButterflyV%blocks(1)%matrix,1,1,descButterflyV,Vin,1,1,descVin,czero,Vinter,1,1,descVinter)	
		call pgemmf90('N','N',blocks%M,Ntest,rank,cone, blocks%ButterflyU%blocks(1)%matrix,1,1,descButterflyU,Vinter,1,1,descVinter,czero,Vout2,1,1,descVout)				
		
		myArows = numroc(blocks%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc(Ntest, nbslpk, mycol, 0, npcol)			
		Vout2 = Vout2-Vout1
		rtemp1 = fnorm(Vout2,myArows,myAcols)**2d0
		rtemp0 = fnorm(Vout1,myArows,myAcols)**2d0
		
		 
		call MPI_ALLREDUCE(rtemp0,fnorm0,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%pgrp(pgno)%Comm,ierr)
		call MPI_ALLREDUCE(rtemp1,fnorm1,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%pgrp(pgno)%Comm,ierr)

		call MPI_Comm_rank(ptree%pgrp(pgno)%Comm,MyID,ierr)
		 
		deallocate(Vin,Vout1,Vout2,Vinter,Fullmat)
		
		if(MyID==0)then
			write(*,*)blocks%row_group,blocks%col_group,'BlockLR error:',sqrt(fnorm1/fnorm0)
		endif	
	endif
	
end subroutine CheckLRError




! integer function rank_approximate_func(group_m, group_n, flag,ker)

    ! use HODLR_DEFS
    ! implicit none

    ! integer i, j, k, mm, nn, edge_head, edge_tail, rank, group_m, group_n, flag
    ! real(kind=8) a, b, c, aa(2), bb(2), cc(2), angle, distance
	! type(kernelquant)::ker
	
	! if(group_m/=group_n)then
		
		! distance=group_dist(group_m,group_n)
		! distance=distance**2d0
		! ! distance=(basis_group(group_m)%center(1)-basis_group(group_n)%center(1))**2+(basis_group(group_m)%center(2)-basis_group(group_n)%center(2))**2+(basis_group(group_m)%center(3)-basis_group(group_n)%center(3))**2
		! ! ! distance=sqrt(distance)
		! angle=4*pi*(basis_group(group_m)%radius)**2/distance
		! rank=int(4*pi*(basis_group(group_n)%radius)**2*angle/ker%wavelength**2)+1
		! ! if(group_m==4 .and. group_n==24)write(*,*)int(rank*ker%rank_approximate_para1),rank,basis_group(group_n)%radius,basis_group(group_n)%radius,angle,distance
		! if (flag==1) then
			! rank_approximate_func=int(rank*ker%rank_approximate_para1**2)
		! elseif (flag==2) then
			! rank_approximate_func=int(rank*ker%rank_approximate_para2**2)
		! elseif (flag==3) then
			! rank_approximate_func=int(rank*ker%rank_approximate_para3**2)
		! endif
	! else 
		! rank_approximate_func = 100000
	! end if
    ! !if (rank==0) then
    ! !    pause
    ! !    continue
    ! !endif

    ! return

! end function rank_approximate_func





end module Bplus_compress
