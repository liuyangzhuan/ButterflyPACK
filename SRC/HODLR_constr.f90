#include "HODLR_config.fi"
module HODLR_constr

! use Butterfly_exact
use Bplus_compress
use Bplus_randomized


contains 


subroutine element_Zmn_user(edge_m,edge_n,value_e,msh,ker)
	
	use HODLR_DEFS
	implicit none
	
	integer edge_m, edge_n
	DT:: value_e
	type(mesh)::msh
	type(kernelquant)::ker	
	
	procedure(F_Z_elem), POINTER :: proc
	value_e=0
	proc => ker%FuncZmn
	call proc(msh%new2old(edge_m),msh%new2old(edge_n),value_e,ker%QuantZmn)

	return
	
end subroutine element_Zmn_user	
	


subroutine HODLR_construction(ho_bf1,option,stats,msh,ker,element_Zmn,ptree)
	
	
    use HODLR_DEFS
    implicit none
	real(kind=8) n1,n2
    integer i, j, ii, ii_inv, jj, kk, iii, jjj,ll
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, itemp,rank0_inner,rank0_outter,ierr
    real T0
	real(kind=8):: rtemp,rel_error,error,t1,t2,tim_tmp,rankrate_inner,rankrate_outter
	integer mm,nn,header_m,header_n,edge_m,edge_n,group_m,group_n,group_m1,group_n1,group_m2,group_n2
	type(matrixblock)::block_tmp,block_tmp1
	DT,allocatable::fullmat(:,:)
	integer level_c,iter,level_cc,level_butterfly,Bidxs,Bidxe
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Z_elem)::element_Zmn
	
	
    ! Memory_direct_forward=0
    ! Memory_butterfly_forward=0
	tim_tmp = 0	
    !tolerance=0.001
    if(ptree%MyID==Main_ID)write (*,*) ''
    ! write (*,*) 'ACA error threshold',tolerance
    if(ptree%MyID==Main_ID)write (*,*) 'SVD error threshold',option%tol_comp
    if(ptree%MyID==Main_ID)write (*,*) ''

    if(ptree%MyID==Main_ID)write(*,*) "constructing Leaf-blocks......"

    n1 = OMP_get_wtime()
    level=0
    flag=0
	! ForwardSymmetricFlag = 0
	allocate (stats%rankmax_of_level(ho_bf1%Maxlevel))
	stats%rankmax_of_level = 0
	allocate (stats%rankmax_of_level_global(ho_bf1%Maxlevel))
	stats%rankmax_of_level_global = 0	
	
	do level_c = 1,ho_bf1%Maxlevel+1
	! do level_c = 1,1
		if(level_c/=ho_bf1%Maxlevel+1)then
			Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
			Bidxe = ho_bf1%levels(level_c)%Bidxe*2
		else
			Bidxs = ho_bf1%levels(level_c)%Bidxs
			Bidxe = ho_bf1%levels(level_c)%Bidxe		
		endif
		
		do ii =Bidxs,Bidxe
		! do ii =Bidxs,Bidxs

			if(ptree%MyID >=ptree%pgrp(ho_bf1%levels(level_c)%BP(ii)%pgno)%head .and. ptree%MyID <=ptree%pgrp(ho_bf1%levels(level_c)%BP(ii)%pgno)%tail)then
				if (level_c/=ho_bf1%Maxlevel+1) then
					if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
						level=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
						if(ptree%MyID==Main_ID)write (*,*) 'constructing level',level
					endif
					if(level_c>=ho_bf1%Maxlevel)t1=OMP_GET_WTIME()	
					call Bplus_compress_N15(ho_bf1%levels(level_c)%BP(ii),option,rtemp,stats,msh,ker,element_Zmn,ptree)				
					
					if(level_c>=ho_bf1%Maxlevel)then
						t2=OMP_GET_WTIME()
						tim_tmp = tim_tmp + t2 - t1
					end if
					
					! if(level==option%level_check)then
						! ! call Bplus_randomized_Exact_test(ho_bf1%levels(level_c)%BP(ii))
						
						! rank0_inner=ho_bf1%levels(level_c)%BP(ii)%LL(2)%rankmax
						! rankrate_inner=1.2d0
						! rank0_outter=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%rankmax
						! rankrate_outter=1.2d0
						! level_butterfly=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level_butterfly
						! call Bplus_randomized_constr(level_butterfly,ho_bf1%levels(level_c)%BP(ii),ho_bf1%levels(level_c)%BP(ii),rank0_inner,rankrate_inner,Bplus_block_MVP_Exact_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Exact_dat,error,'Exact',option,stats,ptree)
						
						
						! stop
					! end if
					
					
					stats%Mem_Comp_for=stats%Mem_Comp_for+rtemp
				else

					if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
						level=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
						if(ptree%MyID==Main_ID)write (*,*) 'constructing level',level
					endif
					call full_construction(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn)
					stats%Mem_Direct_for=stats%Mem_Direct_for+SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
				endif
				! write(*,*)level_c,ii,ho_bf1%levels(level_c)%N_block_forward
				if(level==option%level_check)then
					call BF_compress_test(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn,ptree,stats)
					stop
				end if
			endif
		end do	
		if(ptree%MyID==Main_ID)write(*,*)  'rankmax_of_level so far:',stats%rankmax_of_level
	end do	
	n2 = OMP_get_wtime()
	stats%Time_Fill = stats%Time_Fill + n2-n1 	
	
	
	call MPI_ALLREDUCE(stats%rankmax_of_level,stats%rankmax_of_level_global,ho_bf1%Maxlevel,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write(*,*)  'rankmax_of_level:',stats%rankmax_of_level_global
	if(ptree%MyID==Main_ID)write (*,*) ''
	call MPI_ALLREDUCE(stats%Time_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) 'Total construction time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A26Es14.2)') 'Total construction flops:',rtemp

	if(ptree%MyID==Main_ID)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for butterfly forward blocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write(*,*)rtemp,'MB costed for direct forward blocks'
	if(ptree%MyID==Main_ID)write(*,*)''
	! stop
	
    return

end subroutine HODLR_construction


subroutine full_construction(blocks,msh,ker,element_Zmn)

    use HODLR_DEFS
    implicit none

    integer group_m, group_n, i, j
    integer mm, nn
    integer head_m, head_n, tail_m, tail_n
    DT value_Z
	type(matrixblock)::blocks
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Z_elem)::element_Zmn
	
    mm=blocks%M 
	head_m=blocks%headm
	tail_m=mm+head_m-1
    nn=blocks%N 
	head_n=blocks%headn
	tail_n=nn+head_n-1
	

    allocate (blocks%fullmat(mm,nn))

	!$omp parallel do default(shared) private(j,i,value_Z)
	do j=head_n, tail_n
		do i=head_m, tail_m
			call element_Zmn(i,j,value_Z,msh,ker)
			blocks%fullmat(i-head_m+1,j-head_n+1)=value_Z
		enddo
	enddo
	!$omp end parallel do
    return

end subroutine full_construction




subroutine BF_compress_test(blocks,msh,ker,element_Zmn,ptree,stats)

    use HODLR_DEFS
	use HODLR_Utilities	
    implicit none
    
    type(matrixblock) :: blocks
    real(kind=8) a, b, error,v1,v2
    integer i, j, k, ii, jj, iii,jjj,kk, group_m, group_n, mm, nn, mi, nj,head_m,head_n,Dimn,edge_m,edge_n
    DT value1, value2, ctemp1, ctemp2
	DT,allocatable:: Vin(:,:),Vout1(:,:),Vout2(:,:)
    type(kernelquant)::ker
    type(mesh)::msh
	procedure(Z_elem)::element_Zmn
	type(proctree)::ptree
	type(Hstat)::stats
	integer,allocatable::order_m(:),order_n(:)
	real(kind=8),allocatable::distance_m(:),distance_n(:),center(:)
	
	ctemp1=1.0d0 ; ctemp2=0.0d0
	
	! write(*,*)'h1'

	head_m = blocks%headm
	mm = blocks%M
	
	head_n = blocks%headn
	nn = blocks%N	
	


	! allocate(Vin(nn,1))
	! allocate(Vout1(mm,1))
	! allocate(Vout2(mm,1))
	! do ii=1,nn
		! Vin(ii,1) = random_complex_number()
	! end do
	
	
	! ! write(*,*)'h2'
	! ! write(*,*)blocks%level,Maxlevel_for_blocks
	! ! write(*,*)'h22'
	
	! if(allocated(blocks%fullmat))then
		! ! write(*,*)'h3'
		! call fullmat_block_MVP_dat(blocks,'N',mm,1,Vin,Vout1,ctemp1,ctemp2)
		! ! write(*,*)'h4'
	! else 
		! call BF_block_MVP_dat(blocks,'N',mm,nn,1,Vin,Vout1,ctemp1,ctemp2,ptree,stats)
	! end if	
	
	! do ii=1,mm
		! ctemp1 = 0d0
		! do jj=1,nn
			! ctemp1 = ctemp1 + ker%matZ_glo(msh%new2old(ii+head_m-1),msh%new2old(jj+head_n-1))*Vin(jj,1)
		! end do
		! Vout2(ii,1) = ctemp1
	! end do
	
	! write(*,*)fnorm(Vout2,mm,1), fnorm(Vout2-Vout1,mm,1)/fnorm(Vout2,mm,1)
	
	
	
	allocate(order_m(blocks%M))
	allocate(order_n(blocks%N))
	do ii=1,min(blocks%M,blocks%N)
		call random_number(a)
        call random_number(b)
        order_m(ii)=floor_safe(a*(mm-1))+1
        order_n(ii)=floor_safe(b*(nn-1))+1
		
		! order_m(ii)=ii
		! order_n(ii)=ii		
	enddo
	
	
	!!!!! The following picks the close points first, can be commented out if geometry info is not available
			allocate(distance_m(blocks%M))
			distance_m=Bigvalue
			allocate(distance_n(blocks%N))
			distance_n=Bigvalue

			Dimn = size(msh%xyz,1)
			allocate(center(Dimn))
			center = 0
			head_n = blocks%headn
			do j=1,blocks%N
			  edge_n = head_n-1+j
			  center = center + msh%xyz(1:Dimn,msh%new2old(edge_n))
			enddo
			center = center/blocks%N
			head_m = blocks%headm
			do i=1,blocks%M
				edge_m=head_m-1+i
				distance_m(i) = sum((msh%xyz(1:Dimn,msh%new2old(edge_m))-center(1:Dimn))**2d0)
			enddo		
			deallocate(center)
			call quick_sort(distance_m,order_m,blocks%M)     
			deallocate(distance_m)

			Dimn = size(msh%xyz,1)
			allocate(center(Dimn))
			center = 0
			head_m = blocks%headm
			do i=1,blocks%M
			  edge_m = head_m-1+i
			  center = center + msh%xyz(1:Dimn,msh%new2old(edge_m))
			enddo
			center = center/blocks%M
			head_n = blocks%headn
			do j=1,blocks%N
				edge_n=head_n-1+j
				distance_n(j) = sum((msh%xyz(1:Dimn,msh%new2old(edge_n))-center(1:Dimn))**2d0)
			enddo		
			deallocate(center)
			call quick_sort(distance_n,order_n,blocks%N)     
			deallocate(distance_n)	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	
	v1=0
	v2=0
    ! do i=1,min(mm,nn)
    do i=1,100
		do j=1,100
		
		mi = order_m(i)
		nj = order_n(j)
		
		
        ! iii=int((mi+1)/2)+basis_group(group_m)%head-1
        ! jjj=int((nj+1)/2)+basis_group(group_n)%head-1
        ! ii=2-mod(mi,2)
        ! jj=2-mod(nj,2)
		! call element_Zmn(iii,jjj,ii,jj,value1)

		call element_Zmn(mi+head_m-1,nj+head_n-1,value1,msh,ker)
		
        call Butterfly_value(mi,nj,blocks,value2)
        v1 =v1+abs(value1)**2d0
        v2 =v2+abs(value2)**2d0
		! if(abs(value1)>SafeUnderflow)write (*,*) abs(value1), abs(value2) !, abs(value1-value2)/abs(value1)
		enddo
    enddo
	
	deallocate(order_m)
	deallocate(order_n)
	
	write(*,*)'partial fnorm:',v1,v2,blocks%rankmax
    
    return

end subroutine BF_compress_test








end module HODLR_constr
