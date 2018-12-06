! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so. 

! Developers: Yang Liu, Xiaoye S. Li.
!             (Lawrence Berkeley National Lab, Computational Research Division).

#include "HODLR_config.fi"
module BPACK_constr

! use Butterfly_exact
use Bplus_compress
use Bplus_randomized


contains 


subroutine element_Zmn_user(edge_m,edge_n,value_e,msh,option,ker)
	
	use BPACK_DEFS
	implicit none
	
	integer edge_m, edge_n
	DT:: value_e
	type(mesh)::msh
	type(Hoption)::option
	type(kernelquant)::ker	
	
	procedure(F_Zelem), POINTER :: proc
	value_e=0
	proc => ker%FuncZmn
	call proc(msh%new2old(edge_m),msh%new2old(edge_n),value_e,ker%QuantApp)
	value_e =value_e*option%scale_factor
	return
	
end subroutine element_Zmn_user	
	

subroutine BPACK_construction(bmat,option,stats,msh,ker,element_Zmn,ptree)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	class(*)::bmat
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Zelem)::element_Zmn

	select TYPE(bmat)
    type is (hobf)
		call HODLR_construction(bmat,option,stats,msh,ker,element_Zmn,ptree)
    type is (Hmat)	
		call Hmat_construction(bmat,option,stats,msh,ker,element_Zmn,ptree)
	end select	
	
end subroutine BPACK_construction	
	
subroutine Hmat_construction(h_mat,option,stats,msh,ker,element_Zmn,ptree)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Zelem)::element_Zmn
	integer ierr
    integer i, j, ii, jj, iii, jjj, k, kk
    integer num_blocks
    real*8 T0, T1, T3, T4, rtemp1, rtemp2
    real*8 rtemp
    type(matrixblock), pointer :: blocks,blocks_copy
	integer Maxtmp
	DT:: ctemp
	real(kind=8):: scale_factor
	
    
    call MPI_barrier(ptree%Comm,ierr)

    T0=OMP_get_wtime()

	allocate(stats%rankmax_of_level(0:h_mat%Maxlevel))
	stats%rankmax_of_level = 0	
	allocate (stats%rankmax_of_level_global(0:h_mat%Maxlevel))
	stats%rankmax_of_level_global = 0	    
	
    num_blocks=2**msh%Dist_level
    if (ptree%MyID==ptree%nproc-1 .and. option%verbosity>=0) then
        write(*,*) "   "
        write (*,*) num_blocks*Rows_per_processor, 'total blocks'
    endif
    
	scale_factor = 0
	! compute the largest diagonal entry as the scaling factor
	do i=1, Rows_per_processor
		blocks=>h_mat%Local_blocks(ptree%MyID+1,i)
        do ii=blocks%headm,blocks%headm+blocks%M-1
			call element_Zmn(ii,ii,ctemp,msh,option,ker)
			scale_factor = max(scale_factor,abs(ctemp))
			! write(*,*)ii,abs(ctemp)
		enddo
	enddo
	option%scale_factor = 1d0/scale_factor
	call MPI_ALLREDUCE(MPI_IN_PLACE,option%scale_factor,1,MPI_DOUBLE_PRECISION,MPI_MIN,ptree%Comm,ierr)
	
	if (ptree%MyID==Main_ID .and. option%verbosity>=0) then 
		write(*,*)'element_Zmn is scaled by a factor of:', option%scale_factor
	endif
	
    do i=1, Rows_per_processor
        do j=1, num_blocks
            T3 = OMP_get_wtime()
            blocks=>h_mat%Local_blocks(j,i)
            rtemp1=0. ; rtemp2=0.
			call Hmat_block_construction(blocks,rtemp1,rtemp2,option,stats,msh,ker,element_Zmn,ptree)
			blocks_copy=>h_mat%Local_blocks_copy(j,i)
			
			call Hmat_block_copy('N',blocks_copy,blocks)
			stats%Mem_Comp_for=stats%Mem_Comp_for+rtemp1
			stats%Mem_Direct_for=stats%Mem_Direct_for+rtemp2
			call MPI_barrier(ptree%Comm,ierr)
			T4 = OMP_get_wtime()
            if (ptree%MyID==ptree%nproc-1 .and. option%verbosity>=0) then
				write (*,'(I6,A17,Es14.7,A8,Es14.7,A8,Es14.7)') (i-1)*num_blocks+j,'blocks finished',T4-T3,'secnds', rtemp1, 'Mbytes'
            endif
        enddo
    enddo    
	
	T1 = OMP_get_wtime()
	stats%Time_Fill = stats%Time_Fill + T1-T0 	

	call MPI_ALLREDUCE(stats%rankmax_of_level(0:h_mat%Maxlevel),stats%rankmax_of_level_global(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)
	
	stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
	
	
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)  'rankmax_of_level:',stats%rankmax_of_level_global
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
	call MPI_ALLREDUCE(stats%Time_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Total construction time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A26Es14.2)') 'Total construction flops:',rtemp

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for butterfly forward blocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for direct forward blocks'
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	! stop
    return

end subroutine Hmat_construction	
	

	

subroutine HODLR_construction(ho_bf1,option,stats,msh,ker,element_Zmn,ptree)
	
	
    use BPACK_DEFS
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
	procedure(Zelem)::element_Zmn
	
	
    ! Memory_direct_forward=0
    ! Memory_butterfly_forward=0
	tim_tmp = 0	
    !tolerance=0.001
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
    ! write (*,*) 'ACA error threshold',tolerance
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'SVD error threshold',option%tol_comp
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing Leaf-blocks......"

    n1 = OMP_get_wtime()
    level=0
    flag=0
	! ForwardSymmetricFlag = 0
	allocate (stats%rankmax_of_level(0:ho_bf1%Maxlevel))
	stats%rankmax_of_level = 0
	allocate (stats%rankmax_of_level_global(0:ho_bf1%Maxlevel))
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
			if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(ii)%pgno))then
				if (level_c/=ho_bf1%Maxlevel+1) then
					if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
						level=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
						if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'constructing level',level
					endif
					if(level_c>=ho_bf1%Maxlevel)t1=OMP_GET_WTIME()	
					
					! if(mod(ii,2)==1)then
						call Bplus_compress_N15(ho_bf1%levels(level_c)%BP(ii),option,rtemp,stats,msh,ker,element_Zmn,ptree)				
					! else
						! call BF_delete(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),1)
						! call BF_copy('T',ho_bf1%levels(level_c)%BP(ii-1)%LL(1)%matrices_block(1),ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
					! endif
					
					
					
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
						! call Bplus_randomized_constr(level_butterfly,ho_bf1%levels(level_c)%BP(ii),ho_bf1%levels(level_c)%BP(ii),rank0_inner,rankrate_inner,Bplus_block_MVP_Exact_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Exact_dat,error,'Exact',option,stats,ptree,msh)
						
						
						! stop
					! end if
					
					
					stats%Mem_Comp_for=stats%Mem_Comp_for+rtemp
				else

					if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
						level=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
						if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'constructing level',level
					endif
					call Full_construction(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,option,element_Zmn)
					stats%Mem_Direct_for=stats%Mem_Direct_for+SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
				endif
				! write(*,*)level_c,ii,ho_bf1%levels(level_c)%N_block_forward
				if(level==option%level_check)then
					call BF_compress_test(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn,ptree,option,stats)
					stop
				end if
			endif
		end do	
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)  'rankmax_of_level so far:',stats%rankmax_of_level
	end do	
	n2 = OMP_get_wtime()
	stats%Time_Fill = stats%Time_Fill + n2-n1 	
	
	
	call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel),stats%rankmax_of_level_global(0:ho_bf1%Maxlevel),ho_bf1%Maxlevel+1,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)  'rankmax_of_level:',stats%rankmax_of_level_global
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
	call MPI_ALLREDUCE(stats%Time_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Total construction time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A26Es14.2)') 'Total construction flops:',rtemp

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for butterfly forward blocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for direct forward blocks'
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	! stop
	
    return

end subroutine HODLR_construction


subroutine Full_construction(blocks,msh,ker,option,element_Zmn)

    use BPACK_DEFS
    implicit none

    integer group_m, group_n, i, j
    integer mm, nn
    integer head_m, head_n, tail_m, tail_n
    DT value_Z
	type(matrixblock)::blocks
	type(mesh)::msh
	type(Hoption)::option
	type(kernelquant)::ker
	procedure(Zelem)::element_Zmn
	
    mm=blocks%M 
	head_m=blocks%headm
	tail_m=mm+head_m-1
    nn=blocks%N 
	head_n=blocks%headn
	tail_n=nn+head_n-1
	

    allocate (blocks%fullmat(mm,nn))
	if(blocks%row_group==blocks%col_group)allocate(blocks%ipiv(mm))

	!$omp parallel do default(shared) private(j,i,value_Z)
	do j=head_n, tail_n
		do i=head_m, tail_m
			call element_Zmn(i,j,value_Z,msh,option,ker)
			blocks%fullmat(i-head_m+1,j-head_n+1)=value_Z
		enddo
	enddo
	!$omp end parallel do
    return

end subroutine Full_construction


recursive subroutine Hmat_block_construction(blocks,Memory_far,Memory_near,option,stats,msh,ker,element_Zmn,ptree)

    implicit none
    
    type(matrixblock), pointer :: blocks_son
    type(matrixblock) :: blocks
    real*8 Memory_far, Memory_near, rtemp, Memory_tmp
    integer i, j, k, flag, conv,m,n,ii,jj
	
    type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Zelem)::element_Zmn
	
    if (blocks%style==2) then
		if(option%forwardN15flag==1)then
			call BF_compress_N15(blocks,option,Memory_far,stats,msh,ker,element_Zmn,ptree)
			call BF_sym2asym(blocks)				
		else
			call BF_compress_NlogN(blocks,option,Memory_far,stats,msh,ker,element_Zmn,ptree)
		end if				
    elseif (blocks%style==1) then
		call Full_construction(blocks,msh,ker,option,element_Zmn)
		Memory_near=Memory_near+SIZEOF(blocks%fullmat)/1024.0d3
    elseif (blocks%style==4) then
		do ii=1,2
		do jj=1,2
			blocks_son=>blocks%sons(ii,jj)
			call Hmat_block_construction(blocks_son,Memory_far,Memory_near,option,stats,msh,ker,element_Zmn,ptree)		
		enddo
		enddo
    endif
    
    return                    

end subroutine Hmat_block_construction



end module BPACK_constr
