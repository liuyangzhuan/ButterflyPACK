module matrices_fill
use Butterfly_rightmultiply
use Butterfly_exact
use Butterfly_compress_forward
contains 



subroutine matrices_filling(tolerance)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, kk, iii, jjj,ll
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, itemp
    real T0
	real*8:: tolerance, rtemp,rel_error,error,t1,t2,tim_tmp
    real*8 Memory_direct_forward,Memory_butterfly_forward
	integer mm,nn,header_m,header_n,edge_m,edge_n,group_m,group_n,group_m1,group_n1,group_m2,group_n2
	complex(kind=8)::ctemp,ctemp1,ctemp2
	type(matrixblock)::block_tmp,block_tmp1
	complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:),fullmat5(:,:),fullmat6(:,:)
	integer, allocatable :: ipiv(:)
	complex(kind=8), allocatable :: UU(:,:), VV(:,:),testin(:,:),testout(:,:),Vin(:,:),Vout1(:,:),Vout2(:,:)
	real*8,allocatable :: Singular(:)
	integer level_c,iter,level_cc
	
	
    Memory_direct_forward=0
    Memory_butterfly_forward=0
	tim_tmp = 0	
    !tolerance=0.001
    open (256,file='Info.txt',position='append')
    write (256,*) 'Forward ACA error threshold',tolerance
    write (256,*) 'Forward SVD error threshold',SVD_tolerance_forward
    write (256,*) ''
    close (256)
    write (*,*) ''
    write (*,*) 'ACA error threshold',tolerance
    write (*,*) 'SVD error threshold',SVD_tolerance_forward
    write (*,*) ''

    write(*,*) "Filling Leaf-blocks......"

    T0=secnds(0.0)
    level=0
    flag=0
	ForwardSymmetricFlag = 0
	allocate (rankmax_of_level(Maxlevel_for_blocks))
	rankmax_of_level = 0
	
	
	do level_c = 1,Maxlevel_for_blocks+1
		do ii =1,cascading_factors(level_c)%N_block_forward
            if (level_c/=Maxlevel_for_blocks+1) then
                if (cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
                    level=cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
                    write (*,*) 'level',level,'is filling...'
                endif
				if(level_c>=Maxlevel_for_blocks)t1=OMP_GET_WTIME()									  
				call Bplus_compress_N15(cascading_factors(level_c)%BP(ii),rtemp)				
                
				if(level_c>=Maxlevel_for_blocks)then
					t2=OMP_GET_WTIME()
					tim_tmp = tim_tmp + t2 - t1
				end if
				
				if(level==level_tmp)then
					call Bplus_randomized_Exact_test(cascading_factors(level_c)%BP(ii))
					stop
				end if
				
				
				Memory_butterfly_forward=Memory_butterfly_forward+rtemp
            else

                if (cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
                    level=cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
                    write (*,*) 'level',level,'is filling...'
                endif
                call full_filling(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1))
                Memory_direct_forward=Memory_direct_forward+SIZEOF(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
            endif
			! write(*,*)level_c,ii,cascading_factors(level_c)%N_block_forward
			
			! if(level==level_tmp)then
				! call Butterfly_compress_test(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1))
			! end if
			
		end do
		write(*,*)  'rankmax_of_level so far:',rankmax_of_level
	end do	
	
	
	if(preconditioner==1)then
		allocate(cascading_factors_copy(Maxlevel_for_blocks+1))
		do level_c = 1,Maxlevel_for_blocks+1
			cascading_factors_copy(level_c)%level = cascading_factors(level_c)%level
			cascading_factors_copy(level_c)%N_block_forward = cascading_factors(level_c)%N_block_forward
			allocate(cascading_factors_copy(level_c)%BP(cascading_factors_copy(level_c)%N_block_forward))
			do ii = 1, cascading_factors_copy(level_c)%N_block_forward
				call copy_Bplus(cascading_factors(level_c)%BP(ii),cascading_factors_copy(level_c)%BP(ii))
			end do
		end do		
	end if

	if(fullmatflag==-1 .and. explicitflag ==0)then
		allocate(cascading_factors_ref(Maxlevel_for_blocks+1))
		do level_c = 1,Maxlevel_for_blocks+1
			cascading_factors_ref(level_c)%level = cascading_factors(level_c)%level
			cascading_factors_ref(level_c)%N_block_forward = cascading_factors(level_c)%N_block_forward
			allocate(cascading_factors_ref(level_c)%BP(cascading_factors_ref(level_c)%N_block_forward))
			do ii = 1, cascading_factors_ref(level_c)%N_block_forward
				call copy_Bplus(cascading_factors(level_c)%BP(ii),cascading_factors_ref(level_c)%BP(ii))
				call delete_Bplus_onlyButterflydata(cascading_factors(level_c)%BP(ii))
			end do
		end do		
	end if
	
	write(*,*)  'rankmax_of_level:',rankmax_of_level
	write (*,*) ''
	write (*,*) 'Total filling time:',secnds(T0),'Seconds'

	write(*,*)''
	write(*,*)Memory_butterfly_forward,'MB costed for butterfly forward blocks'
	write(*,*)Memory_direct_forward,'MB costed for direct forward blocks'
	write(*,*)''

    return

end subroutine matrices_filling





! subroutine matrices_filling(tolerance)
	! ! use lapack95
	! ! use blas95
    ! use MODULE_FILE
    ! implicit none

    ! integer i, j, ii, jj, kk, iii, jjj,ll
    ! integer level, blocks, edge, patch, node, group
    ! integer rank, index_near, m, n, length, flag, itemp
    ! real T0
	! real*8:: tolerance, rtemp,rel_error
    ! real*8 Memory_direct_forward,Memory_butterfly_forward
	! integer mm,nn,header_m,header_n,edge_m,edge_n,group_m,group_n,group_m1,group_n1,group_m2,group_n2
	! complex(kind=8)::ctemp,ctemp1,ctemp2
	! type(matrixblock)::block_tmp,block_tmp1
	! complex(kind=8),allocatable::matrixtmp1(:,:),matrixtmp2(:,:),fullmat(:,:),fullmat_eye(:,:),fullmat1(:,:),fullmat2(:,:),fullmat3(:,:),fullmat4(:,:),fullmat5(:,:),fullmat6(:,:)
	! integer, allocatable :: ipiv(:)
	! complex(kind=8), allocatable :: UU(:,:), VV(:,:),testin(:,:),testout(:,:)
	! real*8,allocatable :: Singular(:)
	! integer level_c,iter,level_cc
	
	
    ! Memory_direct_forward=0
    ! Memory_butterfly_forward=0
    ! !tolerance=0.001
    ! open (256,file='Info.txt',position='append')
    ! write (256,*) 'Forward ACA error threshold',tolerance
    ! write (256,*) 'Forward SVD error threshold',SVD_tolerance_forward
    ! write (256,*) ''
    ! if (Forward_method==1) then
        ! write (256,*) 'Forward method: QRSVD'
    ! elseif (Forward_method==2) then
        ! write (256,*) 'Forward method: ACA'
    ! elseif (Forward_method==3) then
        ! write (256,*) 'Forward method: Regulation'
    ! elseif (Forward_method==4) then
        ! write (256,*) 'Forward method: FastSampling'
    ! endif
    ! close (256)
    ! write (*,*) ''
    ! write (*,*) 'ACA error threshold',tolerance
    ! write (*,*) 'SVD error threshold',SVD_tolerance_forward
    ! write (*,*) ''
    ! if (Forward_method==1) then
        ! write (*,*) 'Forward method: QRSVD'
    ! elseif (Forward_method==2) then
        ! write (*,*) 'Forward method: ACA'
    ! elseif (Forward_method==3) then
        ! write (*,*) 'Forward method: Regulation'
    ! elseif (Forward_method==4) then
        ! write (*,*) 'Forward method: FastSampling'
    ! endif

    ! write(*,*) "Filling Leaf-blocks......"

    ! T0=secnds(0.0)
    ! level=0
    ! flag=0
	! ForwardSymmetricFlag = 0
	! allocate (rankmax_of_level(Maxlevel_for_blocks))
	! rankmax_of_level = 0
	
	
	! do level_c = 1,Maxlevel_for_blocks+1
		! do ii =1,cascading_factors(level_c)%N_block_forward
            ! if (cascading_factors(level_c)%matrices_block(ii)%style==2) then
                ! if (cascading_factors(level_c)%matrices_block(ii)%level/=level) then
                    ! level=cascading_factors(level_c)%matrices_block(ii)%level
                    ! write (*,*) 'level',level,'is filling...'
                ! endif

				! ! call Butterfly_compress(cascading_factors(level_c)%matrices_block(ii),rtemp)
                 ! ! call Butterfly_compress_QRSVD(cascading_factors(level_c)%matrices_block(ii),rtemp)
				! ! call Butterfly_compress_QRSVD_Wrong(cascading_factors(level_c)%matrices_block(ii),rtemp)
                ! ! call Butterfly_compress_ID(cascading_factors(level_c)%matrices_block(ii),rtemp)
			! ! write(*,*)cascading_factors(level_c)%matrices_block(ii)%style,level_c,ii
				! call Butterfly_compress_N15(cascading_factors(level_c)%matrices_block(ii),rtemp)
				! call Butterfly_sym2asym(cascading_factors(level_c)%matrices_block(ii))
				
				! ! stop

				! ! ! call forward_butterfly_test(blocks)
				! ! ! pause

				
				! ! ! if(level==level_tmp)then
				! ! ! write(*,*)blocks,'gg'
				! ! ! ! if(level==4)then
				! ! ! ! if(blocks==354)then
				! ! ! ! write(*,*)blocks,'dfdf'
				! ! ! ! call Butterfly_randomized_Exact_test(cascading_factors(level_c)%matrices_block(ii))
				! ! ! call Butterfly_randomized_Exact_test_memfree(cascading_factors(level_c)%matrices_block(ii))
				
				! ! ! ! stop
				! ! ! end if	
				
                ! Memory_butterfly_forward=Memory_butterfly_forward+rtemp
				

				! ! if(level==level_tmp)then            
					! ! call Butterfly_compress_ForwardMatrix_acc(level_c,ii)
					! ! stop
				! ! end if					
	
				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              
					! ! call ACA_and_SeudoSkeleton_compress_ForwardMatrix_acc(level_c,ii)
					! ! stop
				! ! end if				
				
				! ! if(level==level_tmp)then
					! ! call Butterfly_compress_Sblock(level_c,ii)
					! ! ! stop
				! ! end if
							
				! ! if(level==level_tmp)then
					! ! call TFQMR_I_plusButter(level_c,ii)						
					! ! stop
				! ! end if								
									
				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              									
					! ! call Butterfly_compress_inv_of_IplusButter(level_c,ii)		
					! ! ! stop
				! ! end if				
				
				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              									
					! ! call Butterfly_compress_I_minus_BC_Schur(level_c,ii)
					! ! stop
				! ! end if																

				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              
					! ! call TFQMR_I_minus_BC(level_c,ii)		
					! ! stop
				! ! end if		
		
							
				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              
					! ! call TFQMR_recursive_A_and_I_minus_BC(level_c,ii)	
					! ! stop
				! ! end if									

								
				! ! if(level==level_tmp)then
					! ! ! write(*,*)rtemp              
					! ! call SVD_I_plus_Butter_and_OriginalMOMMatrix(level_c,ii)				
					! ! stop
				! ! end if							
							
							
				! ! ! ! if(cascading_factors(level_c)%matrices_block(ii)%level>4)then
					! ! ! ! call Butterfly_MVP_efficience_test(blocks)
				! ! ! ! end if
				
				! ! ! call forward_butterfly_test(blocks)
				! ! ! stop
                
				! !call Butterfly_get_vector_test(blocks)
            ! elseif (cascading_factors(level_c)%matrices_block(ii)%style==1) then
                ! if (cascading_factors(level_c)%matrices_block(ii)%level/=level) then
                    ! level=cascading_factors(level_c)%matrices_block(ii)%level
                    ! write (*,*) 'level',level,'is filling...'
                ! endif
                ! call full_filling(cascading_factors(level_c)%matrices_block(ii))
                ! Memory_direct_forward=Memory_direct_forward+SIZEOF(cascading_factors(level_c)%matrices_block(ii)%fullmat)/1024.0d3
            ! endif
		! end do
	! end do	
	
	
	! if(preconditioner==1)then
		! allocate(cascading_factors_copy(Maxlevel_for_blocks+1))
		! do level_c = 1,Maxlevel_for_blocks+1
			! cascading_factors_copy(level_c)%level = cascading_factors(level_c)%level
			! cascading_factors_copy(level_c)%N_block_forward = cascading_factors(level_c)%N_block_forward
			! allocate(cascading_factors_copy(level_c)%matrices_block(cascading_factors_copy(level_c)%N_block_forward))
			! do ii = 1, cascading_factors_copy(level_c)%N_block_forward
				! cascading_factors_copy(level_c)%matrices_block(ii)%level = cascading_factors(level_c)%matrices_block(ii)%level
				! cascading_factors_copy(level_c)%matrices_block(ii)%col_group = cascading_factors(level_c)%matrices_block(ii)%col_group
				! cascading_factors_copy(level_c)%matrices_block(ii)%row_group = cascading_factors(level_c)%matrices_block(ii)%row_group
				! cascading_factors_copy(level_c)%matrices_block(ii)%style = cascading_factors(level_c)%matrices_block(ii)%style						
				! call copy_butterfly(cascading_factors(level_c)%matrices_block(ii),cascading_factors_copy(level_c)%matrices_block(ii))
			! end do
		! end do		
	! end if

	
    ! write(*,*)  'rankmax_of_level:',rankmax_of_level
    ! write (*,*) ''
    ! write (*,*) 'Total filling time:',secnds(T0),'Seconds'

    ! write(*,*)''
    ! write(*,*)Memory_butterfly_forward,'MB costed for butterfly forward blocks'
    ! write(*,*)Memory_direct_forward,'MB costed for direct forward blocks'
    ! write(*,*)''
    
    ! return

! end subroutine matrices_filling


subroutine full_filling(blocks)

    use MODULE_FILE
    implicit none

    integer group_m, group_n, i, j
    integer mm, nn
    integer head_m, head_n, tail_m, tail_n
    complex(kind=8) value_Z
	type(matrixblock)::blocks

    group_m=blocks%row_group ! Note: row_group and col_group interchanged here   
    group_n=blocks%col_group
    head_m=basis_group(group_m)%head
    tail_m=basis_group(group_m)%tail
    head_n=basis_group(group_n)%head
    tail_n=basis_group(group_n)%tail
    mm=tail_m-head_m+1
    nn=tail_n-head_n+1

    allocate (blocks%fullmat(mm,nn))

    ! if (group_m/=group_n) then
        !$omp parallel do default(shared) private(j,i,value_Z)
        do j=head_n, tail_n
            do i=head_m, tail_m
                call element_Zmn(i,j,value_Z)
                blocks%fullmat(i-head_m+1,j-head_n+1)=value_Z
            enddo
        enddo
        !$omp end parallel do
    ! else
        ! !$omp parallel do default(shared) private(j,i,value_Z)
        ! do j=head_n, tail_n
            ! do i=j, tail_m
                ! call element_Zmn(i,j,value_Z)
                ! blocks%fullmat(i-head_m+1,j-head_n+1)=value_Z
                ! blocks%fullmat(j-head_n+1,i-head_m+1)=value_Z
            ! enddo
        ! enddo
        ! !$omp end parallel do
    ! endif

    return

end subroutine full_filling




subroutine Butterfly_compress_test(matrices_block1)

    use MODULE_FILE
	use Utilities	
    implicit none
    
    type(matrixblock) :: matrices_block1
    real*8 a, b, error
    integer i, j, k, ii, jj, iii,jjj,kk, group_m, group_n, mm, nn, mi, nj
    complex(kind=8) value1, value2, ctemp1, ctemp2
	complex(kind=8),allocatable:: Vin(:,:),Vout1(:,:),Vout2(:,:)
    
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	
	! write(*,*)'h1'

    group_m=matrices_block1%row_group
    group_n=matrices_block1%col_group
    mm=(basis_group(group_m)%tail-basis_group(group_m)%head+1)
    nn=(basis_group(group_n)%tail-basis_group(group_n)%head+1)
    

	allocate(Vin(nn,1))
	allocate(Vout1(mm,1))
	allocate(Vout2(mm,1))
	do ii=1,nn
		Vin(ii,1) = random_complex_number()
	end do
	
	
	! write(*,*)'h2'
	! write(*,*)matrices_block1%level,Maxlevel_for_blocks
	! write(*,*)'h22'
	
	if(matrices_block1%level==Maxlevel_for_blocks+1)then
		! write(*,*)'h3'
		call fullmat_block_MVP_randomized_dat(matrices_block1,'N',mm,1,Vin,Vout1,ctemp1,ctemp2)
		! write(*,*)'h4'
	else 
		call butterfly_block_MVP_randomized_dat(matrices_block1,'N',mm,nn,1,Vin,Vout1,ctemp1,ctemp2)
	end if	
	
	do ii=1,mm
		ctemp1 = 0d0
		do jj=1,nn
			ctemp1 = ctemp1 + matZ_glo(new2old(ii+basis_group(group_m)%head-1),new2old(jj+basis_group(group_n)%head-1))*Vin(jj,1)
		end do
		Vout2(ii,1) = ctemp1
	end do
	
	write(*,*)fnorm(Vout2,mm,1), fnorm(Vout2-Vout1,mm,1)/fnorm(Vout2,mm,1)
	
	
    ! do i=1,3
        ! call random_number(a)
        ! call random_number(b)
        ! mi=int(a*mm)+1
        ! nj=int(b*nn)+1
		
        ! ! iii=int((mi+1)/2)+basis_group(group_m)%head-1
        ! ! jjj=int((nj+1)/2)+basis_group(group_n)%head-1
        ! ! ii=2-mod(mi,2)
        ! ! jj=2-mod(nj,2)
		! ! call element_Zmn(iii,jjj,ii,jj,value1)

		! call element_Zmn(mi+basis_group(group_m)%head-1,nj+basis_group(group_n)%head-1,value1)
		
        ! call Butterfly_value(mi,nj,matrices_block1,value2)
        ! write (*,*) abs(value1), abs(value2), abs(value1-value2)/abs(value1)
    ! enddo
    
    return

end subroutine Butterfly_compress_test








end module matrices_fill
