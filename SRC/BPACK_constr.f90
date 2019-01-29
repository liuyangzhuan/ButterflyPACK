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

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).

#include "ButterflyPACK_config.fi"
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


!**** Initialization of the construction phase
	! N is matrix dimension
	! P is the permutation vector returned
	! N_loc is the local number of rows/columns
	! bmat is the meta-data storing the compressed matrix
	! Coordinates(optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
	! clustertree(optional) is an array of leafsizes in a user-provided cluster tree. clustertree has length 2*nl with nl denoting level of the clustertree.
	! If clustertree is incomplete with 0 element, ButterflyPACK will adjust it to a complete tree and return a modified clustertree.
	! If the hierarchical matrix has more levels than clustertree, the code will generate more levels according to option%xyzsort, option%nogeo, and option%Nmin_leaf
subroutine BPACK_construction_Init(Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates,tree)
	implicit none
	integer Nunk,Ndim
	real(kind=8),optional:: Coordinates(:,:)

    real(kind=8) para
    real(kind=8) tolerance
    integer nn, mm,Maxlevel,give,need
    integer i,j,k,ii,edge,Dimn
	integer nlevel,level
	integer Permutation(Nunk)
	integer,optional:: tree(:)
	integer Nunk_loc
	integer groupm


	type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(Bmatrix)::bmat
	type(proctree)::ptree

	real(kind=8) t1,t2
	character(len=1024)  :: strings
	integer threads_num

	call assert(associated(ker%QuantApp),'ker%QuantApp is not assigned')
	call assert(associated(ker%FuncZmn) .or. associated(ker%FuncHMatVec),'neither ker%FuncZmn nor ker%FuncHMatVec is assigned')

	stats%Flop_Fill=0
	stats%Time_Fill=0


	!**** set thread number here
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',ptree%nproc
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)


	msh%Nunk = Nunk


	t1 = OMP_get_wtime()
	nlevel=0
	if(present(tree))then
		nlevel = ceiling_safe(log(dble(size(tree,1))) / log(2d0))
		Maxlevel=nlevel
		allocate(msh%pretree(2**Maxlevel))
		msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)

		!**** make 0-element node a 1-element node

		! write(*,*)'before adjustment:',msh%pretree
		need = 0
		do ii=1,2**Maxlevel
			if(msh%pretree(ii)==0)need=need+1
		enddo
		do while(need>0)
			give = ceiling_safe(need/dble(2**Maxlevel-need))
			do ii=1,2**Maxlevel
				nn = msh%pretree(ii)
				if(nn>1)then
					msh%pretree(ii) = msh%pretree(ii) - min(min(nn-1,give),need)
					need = need - min(min(nn-1,give),need)
				endif
			enddo
		enddo
		do ii=1,2**Maxlevel
			if(msh%pretree(ii)==0)msh%pretree(ii)=1
		enddo
		! write(*,*)'after adjustment:',msh%pretree
		tree(1:2**Maxlevel) = msh%pretree(1:2**Maxlevel)
	endif

	!**** copy geometry points if present
	if(option%nogeo==0)then
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "User-supplied kernel requiring reorder"
		call assert(present(Coordinates),'geometry points should be provided if option%nogeo==0')
		Ndim = size(Coordinates,1)
		Dimn = Ndim
		allocate (msh%xyz(Dimn,1:msh%Nunk))
		stats%Mem_Peak = stats%Mem_Peak + SIZEOF(msh%xyz)/1024.0d3
		msh%xyz=Coordinates
	endif


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()

	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format......"
    call Cluster_partition(bmat,option,msh,ker,stats,element_Zmn_user,ptree)
	call BPACK_structuring(bmat,option,msh,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()


	!**** return the permutation vector
	select case(option%format)
	case(HODLR)
		msh%idxs = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	case(HMAT)
		msh%idxs = bmat%h_mat%Local_blocks(1,1)%headm
		msh%idxe = bmat%h_mat%Local_blocks(1,1)%headm+bmat%h_mat%Local_blocks(1,1)%M-1
	end select

	Nunk_loc = msh%idxe-msh%idxs+1
	if(ptree%MyID==Main_ID)then
		do edge=1,Nunk
			Permutation(edge) = msh%new2old(edge)
		enddo
	endif


end subroutine BPACK_construction_Init





!**** Computation of the construction phase with matrix entry evaluation
subroutine BPACK_construction_Element(bmat,option,stats,msh,ker,element_Zmn,ptree)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Bmatrix)::bmat
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Zelem)::element_Zmn

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction......"

	select case(option%format)
    case(HODLR)
		call HODLR_construction(bmat%ho_bf,option,stats,msh,ker,element_Zmn,ptree)
    case(HMAT)
		call Hmat_construction(bmat%h_mat,option,stats,msh,ker,element_Zmn,ptree)
	end select

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction finished"

end subroutine BPACK_construction_Element

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
			stats%Mem_Peak = stats%Mem_Peak + rtemp1 + rtemp2 + rtemp1 + rtemp2
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




					! if(mod(ii,2)==1)then
						call Bplus_compress_N15(ho_bf1%levels(level_c)%BP(ii),option,rtemp,stats,msh,ker,element_Zmn,ptree)
					! else
						! call BF_delete(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),1)
						! call BF_copy('T',ho_bf1%levels(level_c)%BP(ii-1)%LL(1)%matrices_block(1),ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
					! endif


					if(level==option%level_check)then
						! call Bplus_randomized_Exact_test(ho_bf1%levels(level_c)%BP(ii))

						rank0_inner=ho_bf1%levels(level_c)%BP(ii)%LL(2)%rankmax
						rankrate_inner=1.2d0
						rank0_outter=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%rankmax
						rankrate_outter=1.2d0
						level_butterfly=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level_butterfly

						t1=OMP_GET_WTIME()
						call BF_randomized(level_butterfly,rank0_outter,rankrate_outter,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),ho_bf1%levels(level_c)%BP(ii),Bplus_block_MVP_Exact_dat,error,'Exact',option,stats,ptree,msh)

						t2=OMP_GET_WTIME()
						tim_tmp = tim_tmp + t2 - t1


						! call Bplus_randomized_constr(level_butterfly,ho_bf1%levels(level_c)%BP(ii),ho_bf1%levels(level_c)%BP(ii),rank0_inner,rankrate_inner,Bplus_block_MVP_Exact_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Exact_dat,error,'Exact',option,stats,ptree,msh)


						if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*)'time_tmp',time_tmp,'randomized_bf time,', tim_tmp,'stats%Time_random,',stats%Time_random



						stop
					end if


					stats%Mem_Comp_for=stats%Mem_Comp_for+rtemp
				else

					if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level/=level) then
						level=ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
						if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'constructing level',level
					endif
					call Full_construction(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,stats,option,element_Zmn)
					stats%Mem_Direct_for=stats%Mem_Direct_for+SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
				endif
				! write(*,*)level_c,ii,ho_bf1%levels(level_c)%N_block_forward
				! if(level>=option%level_check .and. level_c/=ho_bf1%Maxlevel+1)then
					! call BF_compress_test(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn,ptree,option,stats)
					! !stop
				! end if
			endif
		end do
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)  'rankmax_of_level so far:',stats%rankmax_of_level
	end do
	n2 = OMP_get_wtime()
	stats%Time_Fill = stats%Time_Fill + n2-n1

	stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
	stats%Mem_Peak = stats%Mem_Peak + stats%Mem_Fill

	call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel),stats%rankmax_of_level_global(0:ho_bf1%Maxlevel),ho_bf1%Maxlevel+1,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)  'rankmax_of_level:',stats%rankmax_of_level_global
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
	call MPI_ALLREDUCE(stats%Time_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Total construction time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Fill,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A26Es14.2)') 'Total construction flops:',rtemp
	call MPI_ALLREDUCE(time_tmp,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*)'time_tmp',time_tmp

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_Comp_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for butterfly forward blocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_for,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for direct forward blocks'
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	! stop

    return

end subroutine HODLR_construction


subroutine Full_construction(blocks,msh,ker,stats,option,element_Zmn)

    use BPACK_DEFS
    implicit none

    integer group_m, group_n, i, j
    integer mm, nn
    integer head_m, head_n, tail_m, tail_n
    DT value_Z
	type(matrixblock)::blocks
	type(mesh)::msh
	type(Hoption)::option
	type(Hstat)::stats
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
    real*8 Memory_far, Memory_near, rtemp, Memory_tmp,t1,t2
    integer i, j, k, flag, conv,m,n,ii,jj

    type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(Zelem)::element_Zmn
	! t1=OMP_GET_WTIME()
    if (blocks%style==2) then
		if(option%forwardN15flag==1)then
			call BF_compress_N15(blocks,option,Memory_tmp,stats,msh,ker,element_Zmn,ptree)
			call BF_sym2asym(blocks)
		else
			call BF_compress_NlogN(blocks,option,Memory_tmp,stats,msh,ker,element_Zmn,ptree)
		end if
		Memory_far=Memory_far+Memory_tmp
    elseif (blocks%style==1) then
		call Full_construction(blocks,msh,ker,stats,option,element_Zmn)
		Memory_near=Memory_near+SIZEOF(blocks%fullmat)/1024.0d3
    elseif (blocks%style==4) then
		do ii=1,2
		do jj=1,2
			blocks_son=>blocks%sons(ii,jj)
			call Hmat_block_construction(blocks_son,Memory_far,Memory_near,option,stats,msh,ker,element_Zmn,ptree)
		enddo
		enddo
    endif
! t2=OMP_GET_WTIME()
! if(blocks%level==1)write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,t2-t1
    return

end subroutine Hmat_block_construction



end module BPACK_constr
