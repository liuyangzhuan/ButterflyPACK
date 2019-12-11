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
module BPACK_factor
use Bplus_factor
use BPACK_DEFS
use MISC_Utilities
use omp_lib
use BPACK_block_sendrecv
use BPACK_Utilities
use Bplus_randomizedop
use BPACK_Solve_Mul
contains

subroutine BPACK_Factorization(bmat,option,stats,ptree,msh)
    implicit none
	type(Hoption)::option
	type(Hstat)::stats
	type(Bmatrix)::bmat
	type(proctree)::ptree
	type(mesh)::msh

	if(option%precon/=NOPRECON)then
	select case(option%format)
    case(HODLR)
		call HODLR_factorization(bmat%ho_bf,option,stats,ptree,msh)
    case(HMAT)
		call Hmat_Factorization(bmat%h_mat,option,stats,ptree,msh)
	end select
	endif

	if(option%ErrSol==1)then
		call BPACK_Test_Solve_error(bmat,msh%idxe-msh%idxs+1,option,ptree,stats)
	endif

end subroutine BPACK_Factorization


subroutine HODLR_factorization(ho_bf1,option,stats,ptree,msh)

    implicit none

    integer i, j, ii, jj, iii, jjj,index_ij,mm,nn
    integer level, blocks, edge, patch, node, group,level_c,groupm_diag
    integer rank, index_near, m, n, length, flag, itemp
    real T0
	real(kind=8)::rtemp=0
	real(kind=8) tmpfact
    real(kind=8) Memory, Memory_near
	integer,allocatable:: index_old(:),index_new(:)
	integer::block_num,block_num_new,level_butterfly
	integer, allocatable :: ipiv(:)
	integer rowblock,pgno1,pgno2,pgno,ierr,rowblock_inv
	type(matrixblock),pointer::block_o,block_off,block_off1,block_off2
	type(matrixblock)::block_tmp
	real(kind=8) n1,n2,nn1,nn2,flop
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree
	type(mesh)::msh

	nn1 = OMP_get_wtime()

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

	call MPI_barrier(ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Computing block inverse at level Maxlevel+1...'
	level_c = ho_bf1%Maxlevel+1
	do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe

		ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat
		nn = size(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,1)

#if 1
		allocate(ipiv(nn))
		call getrff90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		call getrif90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		deallocate(ipiv)
#else
		call GeneralInverse(nn,nn,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,SafeEps,Flops=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
#endif

		!!!!!!! the forward block BP can be deleted if not used in solution phase

		! write(*,*)fnorm(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,nn,nn)

		stats%Mem_Direct_inv=stats%Mem_Direct_inv+SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3


	end do


	call MPI_barrier(ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Computing block inverse at higher levels...'
	do level_c = ho_bf1%Maxlevel,1,-1

		!!!***** update the forward off-diagonal block by left multiplication of inverse of diagonal blocks in Z: Z_ij^l -> Z_ii^-1*Z_ij^l
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'update forward blocks at level:',level_c

		n1 = OMP_get_wtime()
		do rowblock_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		do rowblock=rowblock_inv*2-1,rowblock_inv*2

		if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(rowblock)%pgno))then
			call Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats,ptree,msh)

			call Bplus_ComputeMemory(ho_bf1%levels(level_c)%BP_inverse_update(rowblock),rtemp)
			stats%Mem_Sblock = stats%Mem_Sblock + rtemp


			! if(level_c==6)then
				! call BF_print_size_rank(ho_bf1%levels(level_c)%matrices_block(rowblock),option%tol_comp)
				! stop
			! end if

		end if
		end do
		end do
		n2 = OMP_get_wtime()
		stats%Time_Sblock=stats%Time_Sblock+n2-n1


		!!!***** compute the inverse of each block 2x2 submatrices whose two off-diagonal blocks are butterflies
		call MPI_barrier(ptree%Comm,ierr)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'compute block inverse at level:',level_c
		n1 = OMP_get_wtime()
		do rowblock = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe

			pgno =  ho_bf1%levels(level_c)%BP_inverse(rowblock)%pgno
			if(IOwnPgrp(ptree,pgno))then
				call Bplus_ReDistribute_Inplace(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2-1),stats,ptree,msh)
				call Bplus_ReDistribute_Inplace(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2),stats,ptree,msh)

				call Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,stats,ptree,msh)
				call Bplus_ComputeMemory(ho_bf1%levels(level_c)%BP_inverse_schur(rowblock),rtemp)
				stats%Mem_SMW=stats%Mem_SMW+rtemp

			endif
		end do
		n2 = OMP_get_wtime()
		stats%Time_Inv=stats%Time_Inv + n2-n1
	end do

	nn2 = OMP_get_wtime()
	stats%Time_Factor = nn2-nn1

	call MPI_ALLREDUCE(stats%Time_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'computing updated forward block time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Time_Inv,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'computing inverse block time:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Time_random(1),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Init:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(2),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_MVP:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(3),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Reconstruct:', rtemp
	call MPI_ALLREDUCE(stats%Time_random(4),rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Onesub:', rtemp
	call MPI_ALLREDUCE(stats%Time_SMW,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_SMW:', rtemp
	call MPI_ALLREDUCE(stats%Time_RedistB,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_RedistB:', rtemp
	call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_RedistV:', rtemp
	call MPI_ALLREDUCE(time_tmp,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*)'time_tmp',time_tmp
	call MPI_ALLREDUCE(stats%Flop_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A21Es14.2)') 'Factorization flops:',rtemp


	stats%Mem_Factor = stats%Mem_SMW + stats%Mem_Sblock + stats%Mem_Direct_inv
	stats%Mem_Peak = stats%Mem_Peak + stats%Mem_Factor

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''
	call MPI_ALLREDUCE(stats%Mem_SMW,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for butterfly inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_Sblock,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for butterfly Sblocks'
	call MPI_ALLREDUCE(stats%Mem_Direct_inv,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for direct inverse blocks'
	call MPI_ALLREDUCE(stats%Mem_int_vec,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)rtemp,'MB costed for storing intermidiate vectors'
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)''


    return

end subroutine HODLR_factorization




subroutine Hmat_Factorization(h_mat,option,stats,ptree,msh)
    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer blocks, level, flag,num_blocks,level_butterfly
    integer Primary_block, kk, i, j, k, intemp, ierr
    integer systime0(8), systime1(8)
    character (len=10) :: date1, time1, zone1
    real*8 rtemp1, rtemp2, rtemp
    real*8 nn1,nn2
    real*8 Memory, Memory_near
    type(matrixblock), pointer :: block
	integer mypgno


	nn1 = OMP_get_wtime()
    call Hmat_LU_TopLevel(h_mat%blocks_root,h_mat,option,stats,ptree,msh)
	nn2 = OMP_get_wtime()
	stats%Time_Factor = nn2-nn1


	call MPI_ALLREDUCE(stats%Time_Factor,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Factor:', rtemp
	call MPI_ALLREDUCE(stats%Time_Direct_LU,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Direct_LU:', rtemp
	call MPI_ALLREDUCE(stats%Time_Add_Multiply,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Add_Multiply:', rtemp
	call MPI_ALLREDUCE(stats%Time_Multiply,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Multiply:', rtemp
	call MPI_ALLREDUCE(stats%Time_XLUM,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_XLUM:', rtemp
	call MPI_ALLREDUCE(stats%Time_Split,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Split:', rtemp
	call MPI_ALLREDUCE(stats%Time_Comm,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Comm:', rtemp
	call MPI_ALLREDUCE(stats%Time_Idle,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_Idle:', rtemp

	call MPI_allreduce(MPI_IN_PLACE,stats%Add_random_CNT(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_INTEGER,MPI_sum,ptree%Comm,ierr)
	call MPI_allreduce(MPI_IN_PLACE,stats%Mul_random_CNT(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_INTEGER,MPI_sum,ptree%Comm,ierr)
	call MPI_allreduce(MPI_IN_PLACE,stats%XLUM_random_CNT(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_INTEGER,MPI_sum,ptree%Comm,ierr)
	call MPI_allreduce(MPI_IN_PLACE,stats%Add_random_Time(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_DOUBLE_PRECISION,MPI_max,ptree%Comm,ierr)
	call MPI_allreduce(MPI_IN_PLACE,stats%Mul_random_Time(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_DOUBLE_PRECISION,MPI_max,ptree%Comm,ierr)
	call MPI_allreduce(MPI_IN_PLACE,stats%XLUM_random_Time(0:h_mat%Maxlevel),h_mat%Maxlevel+1,MPI_DOUBLE_PRECISION,MPI_max,ptree%Comm,ierr)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0) then
        write (*,*) ''
		write(*,*) 'Randomized OPs: Count and Time'

        do level=0, h_mat%Maxlevel
			if(level>option%LRlevel)then
				level_butterfly=0 ! low rank below LRlevel
			else
				level_butterfly = h_mat%Maxlevel - level   ! butterfly
			endif

			if(stats%Add_random_CNT(level)+stats%Mul_random_CNT(level)+stats%XLUM_random_CNT(level)/=0)then
            write (*,'(A7,I5,A17,I5,A7,I8,Es10.2,A7,I8,Es10.2,A12,I8,Es10.2)') " level:",level,'level_butterfly:',level_butterfly,'add:',stats%Add_random_CNT(level),stats%Add_random_Time(level),'mul:',stats%Mul_random_CNT(level),stats%Mul_random_Time(level),'XLUM:',stats%XLUM_random_CNT(level),stats%XLUM_random_time(level)
			endif
		enddo
		! write(*,*)'max inverse butterfly rank:', butterflyrank_inverse
    endif


    if (ptree%MyID==Main_ID .and. option%verbosity>=0) then
        write (*,*) ''
        write (*,*) 'Unpacking all blocks...'
    endif


    num_blocks=2**msh%Dist_level
    do i=1, Rows_per_processor
        do j=1, num_blocks
            block=>h_mat%Local_blocks(j,i)
			if(option%ILU==0)then
			mypgno = msh%basis_group(block%row_group)%pgno
            call unpack_all_blocks_one_node(block,h_mat%Maxlevel,ptree,msh,mypgno)
			endif
			call Hmat_block_ComputeMemory(block,stats%Mem_Factor)
        enddo
    enddo

	stats%Mem_Peak = stats%Mem_Peak + stats%Mem_Factor

    call MPI_verbose_barrier('after printing',ptree)
    if (ptree%MyID==Main_ID .and. option%verbosity>=0) then
        write (*,*) 'Unpacking finished'
        write (*,*) ''
    endif

    return

end subroutine Hmat_Factorization



recursive subroutine Hmat_LU_TopLevel(blocks,h_mat,option,stats,ptree,msh)


    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level, id
    integer i, j, k, ii, jj, kk, flag
    real*8 rtemp1, rtemp2
    real*8 T0, T1,T2,T3
    character (len=10) :: date, time, zone
    integer systime0(8), systime1(8)
    type(global_matricesblock), pointer :: blocks, block_son1, block_son2, block_son3


    level=blocks%level
    if (associated(blocks,h_mat%First_block_eachlevel(level)%father).eqv. .true.) then
       T0 = OMP_get_wtime()
    endif

    if (level/=msh%Dist_level) then
        block_son1=>blocks%sons(1,1)
        call Hmat_LU_TopLevel(block_son1,h_mat,option,stats,ptree,msh)
		if(option%ILU==0)then
			block_son1=>blocks%sons(1,1)
			block_son2=>blocks%sons(1,2)
			block_son3=>blocks%sons(2,1)
			call Hmat_LXM_XUM_TopLevel(block_son1,block_son2,block_son3,h_mat,option,stats,ptree,msh)
			block_son1=>blocks%sons(2,1)
			block_son2=>blocks%sons(1,2)
			block_son3=>blocks%sons(2,2)
			call Hmat_add_multiply_TopLevel(block_son3,'-',block_son1,block_son2,h_mat,option,stats,ptree,msh)
        endif
		block_son1=>blocks%sons(2,2)
        call Hmat_LU_TopLevel(block_son1,h_mat,option,stats,ptree,msh)
    else
        call Hmat_LU_DistLevel(blocks,h_mat,option,stats,ptree,msh)
    endif

    if (associated(blocks,h_mat%First_block_eachlevel(level)%father).eqv. .true.) then

        T1 = OMP_get_wtime()

        if(ptree%MyID==Main_ID .and. option%verbosity>=0)then
            write (*,*) 'level:',level,'is completed',T1-T0,'secconds'
        endif
    endif

    return

end subroutine Hmat_LU_TopLevel

subroutine Hmat_LU_DistLevel(global_block,h_mat,option,stats,ptree,msh)


    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level, id, mm, nn
    integer i, j, k, ii, jj, kk, group_m, group_n, group_start
    type(global_matricesblock), pointer :: global_block
    type(matrixblock), pointer :: blocks
    integer mypgno


    group_start=2**msh%Dist_level-1
    group_m=global_block%row_group
    group_n=global_block%col_group
    mm=group_m-group_start
    nn=group_m-group_start
    id=int((mm-1)/Rows_per_processor)
    i=mm-id*Rows_per_processor
    j=nn
    if (ptree%MyID==id) then
        blocks=>h_mat%Local_blocks(j,i)
        if (associated(global_block,h_mat%First_block_eachlevel(msh%Dist_level)%father).neqv. .true.) then
			if(option%ILU==0)then
			mypgno = msh%basis_group(blocks%row_group)%pgno
			call unpack_all_blocks_one_node(blocks,h_mat%Maxlevel,ptree,msh,mypgno)
			endif
        endif
        if (blocks%row_group/=blocks%col_group) then
            write (*,*) 'Hmat_LU_DistLevel error1!'
        endif
        call Hmat_LU(blocks,h_mat,option,stats,ptree,msh)
		if(option%ILU==0)then
        call pack_all_blocks_one_node(blocks,msh)
		endif
    endif
    call MPI_verbose_barrier('before Hmat_LU_DistLevel returns',ptree)

    return

end subroutine Hmat_LU_DistLevel

subroutine Hmat_add_multiply_TopLevel(block3,chara,block1,block2,h_mat,option,stats,ptree,msh)



    implicit none


	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level_butterfly, num_rows, num_blocks, counts, ierr
    integer level, mm_group, nn_group, group_startm, group_startn
    integer style(3), mark(3), flag, id_33, id_12, send_ID, recv_ID, indices, numcols_perprocessor,color, Local_send_ID
    integer i, j, k, mm, nn, rank, level_blocks, ii, jj, kk, processor, num_processors
    integer iter, iter_total, iter_start, iter_end,groupm,groupn,recv,send,recv_max
    character chara

    real*8 T0, T1, T3, T4
    type(global_matricesblock), pointer :: block1, block2, block3
    type(matrixblock), pointer :: blocks11, blocks22, blocks33,blocks_s,blocks_r

    integer, allocatable :: array_statuses(:,:), cols_per_processor(:,:)
    integer, allocatable :: processor_blocks_map(:,:), rows_processor_map(:,:), cols_processor_map(:,:)
    integer Local_COMM,Local_Myid,Local_Np
	integer  send_count_ind, send_count_dat, S_request1,S_request2, success, recv_count_ind, recv_count_dat, send_count_tot,recv_count_tot
	integer,allocatable::Localmyid2myid(:)
	integer mypgno

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'in Hmat_add_multiply_TopLevel'

    level=block3%level

    if (level<=msh%Dist_level) then

        num_rows=2**(msh%Dist_level-level)
        i=block3%row_group
        group_startm=i*2**(msh%Dist_level-level)-2**msh%Dist_level
        group_startn=group_startm-num_rows

        num_blocks=num_rows*num_rows
        if (num_blocks<=ptree%nproc) then
            num_processors=num_blocks
            allocate(cols_per_processor(0:1,0:num_processors-1))
            numcols_perprocessor=1
            ii=0
            do i=1, num_rows
                do j=1, num_rows
                    ii=ii+1
                    cols_per_processor(0,ii-1)=i
                    cols_per_processor(1,ii-1)=j
                enddo
            enddo
        else
            if (mod(num_blocks,ptree%nproc)/=0) then
                write (*,*) 'Hmat_add_multiply_TopLevel error1!'
            endif
            num_processors=ptree%nproc
            numcols_perprocessor=num_blocks/ptree%nproc
            allocate(cols_per_processor(0:numcols_perprocessor,0:num_processors-1))
            ii=0
            do i=1, num_rows
                do j=1, num_rows
                    ii=ii+1
                    cols_per_processor(0,int((ii-1)/numcols_perprocessor))=i
                    jj=mod(ii,numcols_perprocessor)
                    if (jj==0) then
                        jj=numcols_perprocessor
                    endif
                    cols_per_processor(jj,int((ii-1)/numcols_perprocessor))=j
                enddo
            enddo
        endif

        jj=int(num_blocks/num_processors)

        if (group_startn==0) then
            do i=1, num_rows
                id_33=int((group_startm+i-1)/Rows_per_processor)
                ii=mod(i+group_startm,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (id_33==ptree%MyID) then
                    do j=1, num_rows
                        blocks33=>h_mat%Local_blocks(j+group_startm,ii)
                        call pack_all_blocks_one_node(blocks33,msh)
                    enddo
                endif
            enddo
        endif

        if (ptree%MyID<=num_processors-1) then
            allocate(h_mat%Computing_matricesblock_m(1,1))
            allocate(h_mat%Computing_matricesblock_l(1,num_rows))
            allocate(h_mat%Computing_matricesblock_u(num_rows,1))
        endif

		T3=OMP_get_wtime()
		call MPI_verbose_barrier('after allocate Computing_matricesblock',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3


        T0=OMP_get_wtime()
		do i=1, num_rows
			do j =1,num_rows
				send_ID=int((group_startm+i-1)/Rows_per_processor)
				ii=mod(i+group_startm,Rows_per_processor)
				if (ii==0) then
					ii=Rows_per_processor
				endif

				recv=0
				send=0
				if (ptree%MyID<=num_processors-1)then
				if(i==cols_per_processor(0,ptree%MyID))recv=1
				end if
				if(send_ID==ptree%MyID)send=1

				if(recv==1)blocks_r=>h_mat%Computing_matricesblock_l(1,j)
				if(send==1)blocks_s=>h_mat%Local_blocks(j+group_startn,ii)
				call blocks_partial_bcast(blocks_s,blocks_r,send,recv,send_ID,msh,ptree,option)

				T3=OMP_get_wtime()
				call MPI_verbose_barrier('bcast one blocks in L21',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
			end do
		end do


        ! do recv_ID=0, num_processors-1
            ! i=cols_per_processor(0,recv_ID)
            ! send_ID=int((group_startm+i-1)/Rows_per_processor)
            ! ii=mod(i+group_startm,Rows_per_processor)
            ! if (ii==0) then
                ! ii=Rows_per_processor
            ! endif
            ! do j=1, num_rows
                ! if (send_ID==recv_ID) then
                    ! if (send_ID==ptree%MyID) then
                        ! blocks11=>h_mat%Local_blocks(j+group_startn,ii)
                        ! blocks22=>h_mat%Computing_matricesblock_l(1,j)
                        ! call Hmat_block_copy_MPIdata(blocks22,blocks11,msh)
                    ! endif
                ! else
                    ! !indices=index_of_blocks_for_MPI(i,j)
                    ! if (send_ID==ptree%MyID) then
                        ! counts=0
                        ! indices=0
                        ! !request_global=MPI_request_null
                        ! blocks11=>h_mat%Local_blocks(j+group_startn,ii)
                        ! call blocks_send(blocks11,indices,recv_ID,counts,msh,ptree,option)
                    ! endif
                    ! if (recv_ID==ptree%MyID) then
                        ! counts=0
                        ! indices=0
                        ! blocks11=>h_mat%Computing_matricesblock_l(1,j)
                        ! call blocks_recv(blocks11,indices,send_ID,counts,msh,ptree,option)
                    ! endif
                ! endif
				! T3=OMP_get_wtime()
				! call MPI_verbose_barrier('1')
				! T4=OMP_get_wtime()
				! stats%Time_idle=stats%Time_idle+T4-T3
            ! enddo
        ! enddo



        if (ptree%MyID<=num_processors-1) then
            i=cols_per_processor(0,ptree%MyID)
            do j=1, num_rows
				blocks11=>h_mat%Local_blocks(1,1)
				mypgno = msh%basis_group(blocks11%row_group)%pgno
                blocks11=>h_mat%Computing_matricesblock_l(1,j)
			   call unpack_all_blocks_one_node(blocks11,h_mat%Maxlevel,ptree,msh,mypgno)
            enddo
        endif
        T1=OMP_get_wtime()
        stats%Time_Comm=stats%Time_Comm+T1-T0

        do iter = 1, numcols_perprocessor

            do recv_ID=0, num_processors-1
                i=cols_per_processor(0,recv_ID)
                j=cols_per_processor(iter,recv_ID)
                send_ID=int((group_startm+i-1)/Rows_per_processor)
                ii=mod(i+group_startm,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks11=>h_mat%Local_blocks(j+group_startm,ii)
                        blocks22=>h_mat%Computing_matricesblock_m(1,1)
                        call Hmat_block_copy_MPIdata(blocks22,blocks11,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks33=>h_mat%Local_blocks(j+group_startm,ii)
                        call blocks_send(blocks33,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks33=>h_mat%Computing_matricesblock_m(1,1)
                        call blocks_recv(blocks33,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send 1 blocks in Z_22',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo

            do recv_ID=0, num_processors-1
                i=cols_per_processor(0,recv_ID)
                j=cols_per_processor(iter,recv_ID)
                send_ID=int((group_startm+i-1)/Rows_per_processor)
                ii=mod(i+group_startm,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (send_ID==ptree%MyID) then
                    blocks33=>h_mat%Local_blocks(j+group_startm,ii)
                    call Hmat_block_delete(blocks33)
                endif
                if (recv_ID==ptree%MyID) then
					blocks33=>h_mat%Local_blocks(1,1)
					mypgno = msh%basis_group(blocks33%row_group)%pgno
                    blocks33=>h_mat%Computing_matricesblock_m(1,1)
					call unpack_all_blocks_one_node(blocks33,h_mat%Maxlevel,ptree,msh,mypgno)
                endif
            enddo


			do j =1,num_rows
				recv=0
				if(ptree%MyID<=num_processors-1)then
					if(j==cols_per_processor(iter,ptree%MyID))recv=1
				end if
				call MPI_allreduce(recv,recv_max,1,MPI_integer,MPI_max,ptree%Comm,ierr)
				if(recv_max==1)then  ! don't communicate if in this iteration no receiver exists

					do i=1, num_rows
						send_ID=int((group_startn+i-1)/Rows_per_processor)
						ii=mod(i+group_startm,Rows_per_processor)
						if (ii==0) then
							ii=Rows_per_processor
						endif
						send=0
						recv=0
						if(ptree%MyID<=num_processors-1)then
							if(j==cols_per_processor(iter,ptree%MyID))recv=1
						end if
						if(send_ID==ptree%MyID)send=1

						if(recv==1)blocks_r=>h_mat%Computing_matricesblock_u(i,1)
						if(send==1)blocks_s=>h_mat%Local_blocks(j+group_startm,ii)
						call blocks_partial_bcast(blocks_s,blocks_r,send,recv,send_ID,msh,ptree,option)
						T3=OMP_get_wtime()
						call MPI_verbose_barrier('bcast one blocks in U12',ptree)
						T4=OMP_get_wtime()
						stats%Time_idle=stats%Time_idle+T4-T3

					end do
				end if
			end do



            ! do recv_ID=0, num_processors-1
                ! j=cols_per_processor(iter,recv_ID)
                ! do i=1, num_rows
                    ! send_ID=int((group_startn+i-1)/Rows_per_processor)
                    ! ii=mod(i+group_startm,Rows_per_processor)
                    ! if (ii==0) then
                        ! ii=Rows_per_processor
                    ! endif
                    ! if (send_ID==recv_ID) then
                        ! if (send_ID==ptree%MyID) then
                            ! blocks11=>h_mat%Local_blocks(j+group_startm,ii)
                            ! blocks22=>h_mat%Computing_matricesblock_u(i,1)
                            ! call Hmat_block_copy_MPIdata(blocks22,blocks11,msh)
                        ! endif
                    ! else
                        ! !indices=index_of_blocks_for_MPI(i,j)
                        ! if (send_ID==ptree%MyID) then
                            ! counts=0
                            ! !request_global=MPI_request_null
                            ! indices=0
                            ! blocks22=>h_mat%Local_blocks(j+group_startm,ii)
                            ! call blocks_send(blocks22,indices,recv_ID,counts,msh,ptree,option)
                        ! endif
                        ! if (recv_ID==ptree%MyID) then
                            ! counts=0
                            ! indices=0
                            ! blocks22=>h_mat%Computing_matricesblock_u(i,1)
                            ! call blocks_recv(blocks22,indices,send_ID,counts,msh,ptree,option)
                        ! endif
                    ! endif
					! T3=OMP_get_wtime()
                    ! call MPI_verbose_barrier('3')
					! T4=OMP_get_wtime()
                    ! stats%Time_idle=stats%Time_idle+T4-T3
                ! enddo
            ! enddo



            if (ptree%MyID<=num_processors-1) then
                j=cols_per_processor(iter,ptree%MyID)
                do i=1, num_rows
					blocks22=>h_mat%Local_blocks(1,1)
					mypgno = msh%basis_group(blocks22%row_group)%pgno

					blocks22=>h_mat%Computing_matricesblock_u(i,1)

					call unpack_all_blocks_one_node(blocks22,h_mat%Maxlevel,ptree,msh,mypgno)
                enddo
            endif

            if (ptree%MyID<=num_processors-1) then
                i=cols_per_processor(0,ptree%MyID)
                j=cols_per_processor(iter,ptree%MyID)
                blocks33=>h_mat%Computing_matricesblock_m(1,1)
                do k=1, num_rows
                    blocks11=>h_mat%Computing_matricesblock_l(1,k)
                    blocks22=>h_mat%Computing_matricesblock_u(k,1)
                    call Hmat_add_multiply(blocks33,'-',blocks11,blocks22,h_mat,option,stats,ptree,msh)
                    call Hmat_block_delete(blocks22)
                enddo

            endif

            if (ptree%MyID<=num_processors-1) then
                j=cols_per_processor(iter,ptree%MyID)
                i=cols_per_processor(0,ptree%MyID)
                blocks33=>h_mat%Computing_matricesblock_m(1,1)
                call pack_all_blocks_one_node(blocks33,msh)
            endif

			T3=OMP_get_wtime()
			call MPI_verbose_barrier('before send back Z22',ptree)
			T4=OMP_get_wtime()
			stats%Time_idle=stats%Time_idle+T4-T3

            do send_ID=0, num_processors-1
                i=cols_per_processor(0,send_ID)
                j=cols_per_processor(iter,send_ID)
                recv_ID=int((group_startm+i-1)/Rows_per_processor)
                ii=mod(i+group_startm,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks11=>h_mat%Computing_matricesblock_m(1,1)
                        blocks22=>h_mat%Local_blocks(j+group_startm,ii)
                        call Hmat_block_copy_MPIdata(blocks22,blocks11,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks33=>h_mat%Computing_matricesblock_m(1,1)
                        call blocks_send(blocks33,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks33=>h_mat%Local_blocks(j+group_startm,ii)
                        call blocks_recv(blocks33,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send back 1 blocks in Z_22',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo

            if (ptree%MyID<=num_processors-1) then
                i=cols_per_processor(0,ptree%MyID)
                j=cols_per_processor(iter,ptree%MyID)
                blocks33=>h_mat%Computing_matricesblock_m(1,1)
                call Hmat_block_delete(blocks33)
            endif

        enddo

        if (ptree%MyID<=num_processors-1) then
            i=cols_per_processor(0,ptree%MyID)
            ! !$omp parallel do default(shared) private(j,blocks11)
            do j=1, num_rows
                blocks11=>h_mat%Computing_matricesblock_l(1,j)
                call Hmat_block_delete(blocks11)
            enddo
            ! !$omp end parallel do
        endif

        if (ptree%MyID<=num_processors-1) then
            deallocate(h_mat%Computing_matricesblock_l)
            deallocate(h_mat%Computing_matricesblock_u)
            deallocate(h_mat%Computing_matricesblock_m)
        endif

        deallocate (cols_per_processor)
		T3=OMP_get_wtime()
		call MPI_verbose_barrier('before Hmat_add_multiply_TopLevel finish',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3

    else

        !call Hmat_add_multiply(blocks3,chara,blocks1,blocks2)

    endif

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'out Hmat_add_multiply_TopLevel'
    return

end subroutine Hmat_add_multiply_TopLevel

subroutine Hmat_LXM_XUM_TopLevel(blocks_m,blocks_u,blocks_l,h_mat,option,stats,ptree,msh)



    implicit none


	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level, id_l, id_u, id_m, group_m, group_n, iteration_total, indices
    integer mm_group, nn_group, num_rows, num_groups, iter, iter_begin, iter_end,groupm,groupn,send,recv
    integer style(3), mark, style_m, group_startm, group_startn, mm, nn, num_processors
    integer i, j, k, blockson_m(2,2), send_ID, recv_ID,color, Local_send_ID
    integer rank, count1, count2, ii, jj, kk, counts, rows_per_processor_LUM

    real (kind=8) T0, T1,T3,T4
    type(global_matricesblock), pointer :: blocks_l, blocks_u, blocks_m
    type(matrixblock), pointer :: blocks_ll, blocks_uu, blocks_mm, blocks_kk, blocks_s,blocks_r
    INTEGER,allocatable :: array_statuses(:,:)
    integer Local_COMM,Local_Myid,Local_Np
	integer  send_count_ind, send_count_dat, S_request1,S_request2, success, recv_count_ind, recv_count_dat, send_count_tot,recv_count_tot
	integer,allocatable::Localmyid2myid(:)
	integer mypgno

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'in Hmat_LXM_XUM_TopLevel'

    level=blocks_m%level

    if (level<=msh%Dist_level) then

        num_rows=2**(msh%Dist_level-level)
        i=blocks_l%row_group
        group_startm=i*2**(msh%Dist_level-level)-2**msh%Dist_level
        j=blocks_l%col_group
        group_startn=j*2**(msh%Dist_level-level)-2**msh%Dist_level

        if (num_rows<=ptree%nproc/2) then
            num_processors=num_rows
            rows_per_processor_LUM=1
        else
            num_processors=ptree%nproc/2
            rows_per_processor_LUM=2*num_rows/ptree%nproc
        endif

        !********************************************************************************************************

        if (group_startn==0) then
            do k=1, num_rows
                id_l=int((group_startm+k-1)/Rows_per_processor)
                id_u=int((group_startn+k-1)/Rows_per_processor)
                ii=mod(k+group_startm,Rows_per_processor)
                jj=mod(k+group_startn,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (jj==0) then
                    jj=Rows_per_processor
                endif
                if (id_l==ptree%MyID) then
                    do i=1, num_rows
                        blocks_ll=>h_mat%Local_blocks(i+group_startn,ii)
						call pack_all_blocks_one_node(blocks_ll,msh)
                    enddo
                endif
                if (id_u==ptree%MyID) then
                    do j=1, num_rows
                        blocks_uu=>h_mat%Local_blocks(j+group_startm,jj)
						call pack_all_blocks_one_node(blocks_uu,msh)
                    enddo
                endif
            enddo
        endif

		T3=OMP_get_wtime()
		call MPI_verbose_barrier('before allocate Computing_matricesblock',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3

        if (ptree%MyID<=num_processors-1) then
            allocate(h_mat%Computing_matricesblock_m(num_rows,1))
            allocate(h_mat%Computing_matricesblock_l(1,num_rows))
        endif
        if (ptree%MyID>=ptree%nproc-num_processors) then
            allocate(h_mat%Computing_matricesblock_m(1,num_rows))
            allocate(h_mat%Computing_matricesblock_u(num_rows,1))
        endif

		T3=OMP_get_wtime()
		call MPI_verbose_barrier('after allocate Computing_matricesblock',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3

		T0=OMP_get_wtime()

        do i=1, num_rows
            send_ID=int((group_startm+i-1)/Rows_per_processor)
            ii=mod(i+group_startm,Rows_per_processor)
            if (ii==0) then
                ii=Rows_per_processor
            endif
            recv_ID=int((i-1)/rows_per_processor_LUM)
            do j=1, num_rows
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks_ll=>h_mat%Local_blocks(j+group_startn,ii)
                        blocks_mm=>h_mat%Computing_matricesblock_l(1,j)
                        call Hmat_block_copy_MPIdata(blocks_mm,blocks_ll,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_ll=>h_mat%Local_blocks(j+group_startn,ii)
                        call blocks_send(blocks_ll,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
                        call blocks_recv(blocks_ll,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send on 1 blocks in Z_21',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo
        enddo

        do j=1, num_rows
            recv_ID=int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors
            do i=1, num_rows
                send_ID=int((group_startn+i-1)/Rows_per_processor)
                ii=mod(i+group_startn,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                indices=0
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks_uu=>h_mat%Local_blocks(j+group_startm,ii)
                        blocks_mm=>h_mat%Computing_matricesblock_u(i,1)
                        call Hmat_block_copy_MPIdata(blocks_mm,blocks_uu,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        !request_global=MPI_request_null
                        indices=0
                        blocks_uu=>h_mat%Local_blocks(j+group_startm,ii)
                        call blocks_send(blocks_uu,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                        call blocks_recv(blocks_uu,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send on 1 blocks in Z_12',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo
        enddo
		T1=OMP_get_wtime()
        stats%Time_Comm=stats%Time_Comm+T1-T0

        do i=1, num_rows
            id_l=int((group_startm+i-1)/Rows_per_processor)
            id_u=int((group_startn+i-1)/Rows_per_processor)
            if (id_l==ptree%MyID) then
                ii=mod(i+group_startm,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                do j=1, num_rows
                    blocks_ll=>h_mat%Local_blocks(j+group_startn,ii)
                    call Hmat_block_delete(blocks_ll)
                enddo
            endif
            if (id_u==ptree%MyID) then
                ii=mod(i+group_startn,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                do j=1, num_rows
                    blocks_uu=>h_mat%Local_blocks(j+group_startm,ii)
                    call Hmat_block_delete(blocks_uu)
                enddo
            endif
        enddo

        do i=1, num_rows
            id_l=int((i-1)/rows_per_processor_LUM)
            if (id_l==ptree%MyID) then
                do j=1, num_rows
					blocks_ll=>h_mat%Local_blocks(1,1)
					mypgno = msh%basis_group(blocks_ll%row_group)%pgno
                    blocks_ll=>h_mat%Computing_matricesblock_l(1,j)

					call unpack_all_blocks_one_node(blocks_ll,h_mat%Maxlevel,ptree,msh,mypgno)
                enddo
            endif
        enddo
        do j=1, num_rows
            id_u=int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors
            if (id_u==ptree%MyID) then
                do i=1, num_rows
					blocks_uu=>h_mat%Local_blocks(1,1)
					mypgno = msh%basis_group(blocks_uu%row_group)%pgno
                    blocks_uu=>h_mat%Computing_matricesblock_u(i,1)

					call unpack_all_blocks_one_node(blocks_uu,h_mat%Maxlevel,ptree,msh,mypgno)
                enddo
            endif
        enddo


        iteration_total=int((num_rows-1)/Rows_per_processor)

        do iter=0, iteration_total

            iter_begin=iter*Rows_per_processor+1
            iter_end=min(iter_begin+Rows_per_processor-1,num_rows)

			T3=OMP_get_wtime()
			call MPI_verbose_barrier('before transfer one column/row in U11/L11',ptree)
			T4=OMP_get_wtime()
			stats%Time_idle=stats%Time_idle+T4-T3

			T0=OMP_get_wtime()

            do i=iter_begin, iter_end
                send_ID=int((group_startn+i-1)/Rows_per_processor)
                ii=mod(i+group_startn,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                do j=1, i

					recv=0
					send=0
					if (ptree%MyID>=ptree%nproc-num_processors .and. ptree%MyID<=ptree%nproc-1)recv=1
					if(send_ID==ptree%MyID)send=1

					if(recv==1)blocks_r=>h_mat%Computing_matricesblock_m(1,j)
					if(send==1)blocks_s=>h_mat%Local_blocks(j+group_startn,ii)
					call blocks_partial_bcast(blocks_s,blocks_r,send,recv,send_ID,msh,ptree,option)


					T3=OMP_get_wtime()
                    call MPI_verbose_barrier('bcast one blocks in L11',ptree)
					T4=OMP_get_wtime()
                    stats%Time_idle=stats%Time_idle+T4-T3
					! T0=OMP_get_wtime()

					! !indices=index_of_blocks_for_MPI(i,j)
                    ! do k=ptree%nproc-num_processors, ptree%nproc-1
                        ! if (send_ID==k) then
                            ! if (send_ID==ptree%MyID) then
                                ! blocks_kk=>h_mat%Local_blocks(j+group_startn,ii)
                                ! blocks_mm=>h_mat%Computing_matricesblock_m(1,j)
                                ! call Hmat_block_copy_MPIdata(blocks_mm,blocks_kk,msh)
                            ! endif
                        ! else
                            ! if (send_ID==ptree%MyID) then
                                ! counts=0
                                ! !request_global=MPI_request_null
                                ! indices=0
                                ! blocks_mm=>h_mat%Local_blocks(j+group_startn,ii)
                                ! call blocks_send(blocks_mm,indices,k,counts,msh,ptree,option)
                            ! endif
                            ! if (k==ptree%MyID) then
                                ! counts=0
                                ! indices=0
                                ! blocks_mm=>h_mat%Computing_matricesblock_m(1,j)
                                ! call blocks_recv(blocks_mm,indices,send_ID,counts,msh,ptree,option)
                            ! endif
                        ! endif
						! T3=OMP_get_wtime()
						! call MPI_verbose_barrier('4')
						! T4=OMP_get_wtime()
						! stats%Time_idle=stats%Time_idle+T4-T3
                    ! enddo

                enddo
            enddo
            do j=iter_begin, iter_end
                do i=1, j
                    send_ID=int((group_startn+i-1)/Rows_per_processor)
                    ii=mod(i+group_startn,Rows_per_processor)
                    if (ii==0) then
                        ii=Rows_per_processor
                    endif

					recv=0
					send=0
					if (ptree%MyID<=num_processors-1)recv=1
					if(send_ID==ptree%MyID)send=1

					if(recv==1)blocks_r=>h_mat%Computing_matricesblock_m(i,1)
					if(send==1)blocks_s=>h_mat%Local_blocks(j+group_startn,ii)
					call blocks_partial_bcast(blocks_s,blocks_r,send,recv,send_ID,msh,ptree,option)

					T3=OMP_get_wtime()
                    call MPI_verbose_barrier('bcast one blocks in U11',ptree)
					T4=OMP_get_wtime()
                    stats%Time_idle=stats%Time_idle+T4-T3
					! T0=OMP_get_wtime()

					! ! do k=0, num_processors-1
                        ! ! if (send_ID==k) then
                            ! ! if (send_ID==ptree%MyID) then
                                ! ! blocks_kk=>h_mat%Local_blocks(j+group_startn,ii)
                                ! ! blocks_mm=>h_mat%Computing_matricesblock_m(i,1)
                                ! ! call Hmat_block_copy_MPIdata(blocks_mm,blocks_kk,msh)
                            ! ! endif
                        ! ! else
                            ! ! if (send_ID==ptree%MyID) then
                                ! ! counts=0
                                ! ! !request_global=MPI_request_null
                                ! ! indices=0
                                ! ! blocks_mm=>h_mat%Local_blocks(j+group_startn,ii)
                                ! ! call blocks_send(blocks_mm,indices,k,counts,msh,ptree,option)
                            ! ! endif
                            ! ! if (k==ptree%MyID) then
                                ! ! counts=0
                                ! ! indices=0
                                ! ! blocks_mm=>h_mat%Computing_matricesblock_m(i,1)
                                ! ! call blocks_recv(blocks_mm,indices,send_ID,counts,msh,ptree,option)
                            ! ! endif
                        ! ! endif
						! ! T3=OMP_get_wtime()
						! ! call MPI_verbose_barrier('5')
						! ! T4=OMP_get_wtime()
						! ! stats%Time_idle=stats%Time_idle+T4-T3
                    ! ! enddo
                enddo
            enddo
			T1=OMP_get_wtime()
            stats%Time_Comm=stats%Time_Comm+T1-T0


            if (ptree%MyID<=num_processors-1) then
                do j=iter_begin, iter_end
                    do i=1, j
						blocks_mm=>h_mat%Local_blocks(1,1)
						mypgno = msh%basis_group(blocks_mm%row_group)%pgno
                        blocks_mm=>h_mat%Computing_matricesblock_m(i,1)

						call unpack_all_blocks_one_node(blocks_mm,h_mat%Maxlevel,ptree,msh,mypgno)
                    enddo
                enddo
            endif
            if (ptree%MyID>=ptree%nproc-num_processors) then
                do i=iter_begin, iter_end
                    do j=1, i
						blocks_mm=>h_mat%Local_blocks(1,1)
						mypgno = msh%basis_group(blocks_mm%row_group)%pgno
                        blocks_mm=>h_mat%Computing_matricesblock_m(1,j)

						call unpack_all_blocks_one_node(blocks_mm,h_mat%Maxlevel,ptree,msh,mypgno)
                    enddo
                enddo
            endif

            if (ptree%MyID<=num_processors-1) then
                do i=1, num_rows
                    if (int((i-1)/rows_per_processor_LUM)==ptree%MyID) then
                        do j=iter_begin, iter_end
                            blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
                            if (j>1) then
                                do k=1, j-1
                                    blocks_kk=>h_mat%Computing_matricesblock_l(1,k)
                                    blocks_mm=>h_mat%Computing_matricesblock_m(k,1)

                                    call Hmat_add_multiply(blocks_ll,'-',blocks_kk,blocks_mm,h_mat,option,stats,ptree,msh)
                                    call Hmat_block_delete(blocks_mm)

                                enddo
                            endif
                            blocks_mm=>h_mat%Computing_matricesblock_m(j,1)

                            call Hmat_XUM(blocks_mm,blocks_ll,h_mat,option,stats,ptree,msh)
                            call Hmat_block_delete(blocks_mm)

                        enddo
                    endif
                enddo
            endif
            if (ptree%MyID>=ptree%nproc-num_processors) then
                do j=1, num_rows
                    if (int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors==ptree%MyID) then
                        do i=iter_begin, iter_end
                            blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                            if (i>1) then
                                do k=1, i-1
                                    blocks_mm=>h_mat%Computing_matricesblock_m(1,k)
                                    blocks_kk=>h_mat%Computing_matricesblock_u(k,1)

                                    call Hmat_add_multiply(blocks_uu,'-',blocks_mm,blocks_kk,h_mat,option,stats,ptree,msh)
                                    call Hmat_block_delete(blocks_mm)

                                enddo
                            endif
                            blocks_mm=>h_mat%Computing_matricesblock_m(1,i)
                            call Hmat_LXM(blocks_mm,blocks_uu,h_mat,option,stats,ptree,msh)
                            call Hmat_block_delete(blocks_mm)
                        enddo
                    endif
                enddo
            endif
       enddo

        !**************************************write back******************************************************

        do i=1, num_rows
            if (int((i-1)/rows_per_processor_LUM)==ptree%MyID) then
                do j=1, num_rows
                    blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
					call pack_all_blocks_one_node(blocks_ll,msh)
                enddo
            endif
        enddo
        do j=1, num_rows
            if (int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors==ptree%MyID) then
                do i=1, num_rows
                    blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                    call pack_all_blocks_one_node(blocks_uu,msh)
                enddo
            endif
        enddo
		T3=OMP_get_wtime()
		call MPI_verbose_barrier('before send back Z_12/Z_21',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3

		T0=OMP_get_wtime()
        do i=1, num_rows
            send_ID=int((i-1)/rows_per_processor_LUM)
            recv_ID=int((group_startm+i-1)/Rows_per_processor)
            ii=mod(i+group_startm,Rows_per_processor)
            if (ii==0) then
                ii=Rows_per_processor
            endif
            do j=1, num_rows
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
                        blocks_mm=>h_mat%Local_blocks(j+group_startn,ii)
                        call Hmat_block_copy_MPIdata(blocks_mm,blocks_ll,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
                        call blocks_send(blocks_ll,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_ll=>h_mat%Local_blocks(j+group_startn,ii)
                        call blocks_recv(blocks_ll,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send back 1 blocks in Z_21',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo
        enddo
        do j=1, num_rows
            send_ID=int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors
            do i=1, num_rows
                recv_ID=int((group_startn+i-1)/Rows_per_processor)
                ii=mod(i+group_startn,Rows_per_processor)
                if (ii==0) then
                    ii=Rows_per_processor
                endif
                if (send_ID==recv_ID) then
                    if (send_ID==ptree%MyID) then
                        blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                        blocks_mm=>h_mat%Local_blocks(j+group_startm,ii)
                        call Hmat_block_copy_MPIdata(blocks_mm,blocks_uu,msh)
                    endif
                else
                    if (send_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                        call blocks_send(blocks_uu,indices,recv_ID,counts,msh,ptree,option)
                    endif
                    if (recv_ID==ptree%MyID) then
                        counts=0
                        indices=0
                        blocks_uu=>h_mat%Local_blocks(j+group_startm,ii)
                        call blocks_recv(blocks_uu,indices,send_ID,counts,msh,ptree,option)
                    endif
                endif
				T3=OMP_get_wtime()
				call MPI_verbose_barrier('p2p send back 1 blocks in Z_12',ptree)
				T4=OMP_get_wtime()
				stats%Time_idle=stats%Time_idle+T4-T3
            enddo
        enddo

		T1=OMP_get_wtime()
        stats%Time_Comm=stats%Time_Comm+T1-T0

        do i=1, num_rows
            id_l=int((i-1)/rows_per_processor_LUM)
            if (id_l==ptree%MyID) then
                do j=1, num_rows
                    blocks_ll=>h_mat%Computing_matricesblock_l(1,j)
                    call Hmat_block_delete(blocks_ll)
                enddo
            endif
        enddo
        do j=1, num_rows
            id_u=int((j-1)/rows_per_processor_LUM)+ptree%nproc-num_processors
            if (id_u==ptree%MyID) then
                do i=1, num_rows
                    blocks_uu=>h_mat%Computing_matricesblock_u(i,1)
                    call Hmat_block_delete(blocks_uu)
                enddo
            endif
        enddo

        if (ptree%MyID<=num_processors-1) then
            deallocate (h_mat%Computing_matricesblock_m)
            deallocate (h_mat%Computing_matricesblock_l)
        endif
        if (ptree%MyID>=ptree%nproc-num_processors) then
            deallocate (h_mat%Computing_matricesblock_m)
            deallocate (h_mat%Computing_matricesblock_u)
        endif

		T3=OMP_get_wtime()
		call MPI_verbose_barrier('before Hmat_LXM_XUM_TopLevel finish',ptree)
		T4=OMP_get_wtime()
		stats%Time_idle=stats%Time_idle+T4-T3

    else

        !call H_Solve_XLM(blocks_l,blocks_m)

    endif
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'out Hmat_LXM_XUM_TopLevel'

    return

end subroutine Hmat_LXM_XUM_TopLevel



recursive subroutine Hmat_LU(blocks,h_mat,option,stats,ptree,msh)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level, blocks_son(2,2), id
    integer i, j, k, ii, jj, kk, flag
    real*8 rtemp1, rtemp2
    real*8  T0, T1
    type(matrixblock):: blocks
    type(matrixblock), pointer :: block_son1, block_son2, block_son3



    if (blocks%style==4) then
        block_son1=>blocks%sons(1,1)
        call Hmat_LU(block_son1,h_mat,option,stats,ptree,msh)

		if(option%ILU==0)then
			block_son1=>blocks%sons(1,1)
			block_son2=>blocks%sons(2,1)
			call Hmat_XUM(block_son1,block_son2,h_mat,option,stats,ptree,msh)

			block_son1=>blocks%sons(1,1)
			block_son2=>blocks%sons(1,2)
			call Hmat_LXM(block_son1,block_son2,h_mat,option,stats,ptree,msh)

			block_son1=>blocks%sons(2,1)
			block_son2=>blocks%sons(1,2)
			block_son3=>blocks%sons(2,2)
			call Hmat_add_multiply(block_son3,'-',block_son1,block_son2,h_mat,option,stats,ptree,msh)
        endif

        block_son1=>blocks%sons(2,2)
        call Hmat_LU(block_son1,h_mat,option,stats,ptree,msh)
    else
        call Full_LU(blocks,option,stats)
    endif


    return

end subroutine Hmat_LU

recursive subroutine Hmat_add_multiply(block3,chara,block1,block2,h_mat,option,stats,ptree,msh)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

	integer rank0
    integer level_butterfly
    integer style(3), mark(3)
    integer ii,jj,i, j, k, level, mm, nn, rank, level_blocks,kk,found_flag,group_m,group_n,m,n
    character chara
    DT,allocatable:: matrixtemp(:,:)
    real*8 T0,T1,error
    type(matrixblock), pointer :: matrices_tmpblock
    type(matrixblock),target :: block2,block1, block3
    type(matrixblock), pointer :: block1_son, block2_son, block3_son

    style(3)=block3%style
    style(1)=block1%style
    style(2)=block2%style


    level_blocks=block3%level

	if (style(3)==4) then   !!! modified by Yang Liu, hybrid butterfly-LR treatment
		if((style(1)/=4 .or. style(2)/=4))then
			T0 = OMP_get_wtime()
			call Hmat_add_multiply_Hblock3(block3,chara,block1,block2,h_mat,option,stats,ptree,msh)
			T1 = OMP_get_wtime()
		else
			T0 = OMP_get_wtime()
			if (style(1)/=4) then
				allocate(block1%sons(2,2))
				call BF_split(block1,block1,ptree,stats,msh,option)
			endif
			if (style(2)/=4) then
				allocate(block2%sons(2,2))
				call BF_split(block2,block2,ptree,stats,msh,option)
			endif
			if (style(3)/=4) then
				allocate(block3%sons(2,2))
				call BF_split(block3,block3,ptree,stats,msh,option)
			endif
			T1 = OMP_get_wtime()
			stats%Time_split = stats%Time_split + T1-T0

            block1_son=>block1%sons(1,1)
            block2_son=>block2%sons(1,1)
            block3_son=>block3%sons(1,1)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(1,2)
            block2_son=>block2%sons(2,1)
            block3_son=>block3%sons(1,1)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(1,1)
            block2_son=>block2%sons(1,2)
            block3_son=>block3%sons(1,2)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(1,2)
            block2_son=>block2%sons(2,2)
            block3_son=>block3%sons(1,2)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(2,1)
            block2_son=>block2%sons(1,1)
            block3_son=>block3%sons(2,1)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(2,2)
            block2_son=>block2%sons(2,1)
            block3_son=>block3%sons(2,1)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(2,1)
            block2_son=>block2%sons(1,2)
            block3_son=>block3%sons(2,2)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
            block1_son=>block1%sons(2,2)
            block2_son=>block2%sons(2,2)
            block3_son=>block3%sons(2,2)
            call Hmat_add_multiply(block3_son,chara,block1_son,block2_son,h_mat,option,stats,ptree,msh)
		end if
    elseif(style(3)==1)then
            call Full_add_multiply(block3,chara,block1,block2,h_mat,option,stats,ptree,msh)
	elseif(style(3)==2)then
		T0 = OMP_get_wtime()

		h_mat%blocks_1 => block1
		h_mat%blocks_2 => block2
		rank0 = block3%rankmax
		call BF_randomized(block3%pgno,block3%level_butterfly,rank0,option%rankrate,block3,h_mat,BF_block_MVP_Add_Multiply_dat,error,'Add_Multiply',option,stats,ptree,msh,chara)
		T1 = OMP_get_wtime()
		stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
		stats%Time_Add_Multiply=stats%Time_Add_Multiply+T1-T0
		stats%Add_random_Time(level_blocks)=stats%Add_random_Time(level_blocks)+T1-T0
		stats%Add_random_CNT(level_blocks)=stats%Add_random_CNT(level_blocks)+1

    endif


#if 0
	if(dumping==0)then
	if(style(3)==2)then
	    group_m=block3%row_group
        group_n=block3%col_group
		write(*,*)group_m,group_n,m,n
		if(group_m==13 .and. group_n==11)then
			dumping=1
			m = msh%basis_group(group_m)%tail -msh%basis_group(group_m)%head + 1
			n = msh%basis_group(group_n)%tail -msh%basis_group(group_n)%head + 1
			allocate(matrixtemp(m,n))
			call butterfly_block_fullextract(block3,matrixtemp)
			write(*,*)'dumping:', group_m,group_n,m,n
			do ii=1,m
			do jj=1,n
				write(777,*)dble(matrixtemp(ii,jj)),aimag(matrixtemp(ii,jj))
			enddo
			enddo
			deallocate(matrixtemp)
			write(*,*)'dumping done'
			! stop
		endif
	endif
	endif
#endif

    return

end subroutine Hmat_add_multiply



recursive subroutine Hmat_LXM(blocks_l,blocks_m,h_mat,option,stats,ptree,msh)
    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer style(3), mark, style_m
    integer i, j, k, ii
    integer mm, nn, rank
    real(kind=8) T0,T1,error
    type(matrixblock), pointer :: blocks1, blocks2, blocks3
    type(matrixblock) :: blocks_l, blocks_m
    DT:: ctemp
	integer rank0

    if (blocks_m%style==4) then
		T0 = OMP_get_wtime()
        if (blocks_l%style/=4) then
			allocate(blocks_l%sons(2,2))
			call BF_split(blocks_l,blocks_l,ptree,stats,msh,option)
        endif
        T1 = OMP_get_wtime()
		stats%Time_Split=stats%Time_Split+T1-T0

        blocks1=>blocks_l%sons(1,1)
        blocks2=>blocks_m%sons(1,1)
        call Hmat_LXM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks2=>blocks_m%sons(1,2)
        call Hmat_LXM(blocks1,blocks2,h_mat,option,stats,ptree,msh)

        blocks2=>blocks_m%sons(1,1)
        blocks1=>blocks_l%sons(2,1)
        blocks3=>blocks_m%sons(2,1)
		! write(*,*)'1b',blocks1%style,blocks2%style,blocks3%style
        call Hmat_add_multiply(blocks3,'-',blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks2=>blocks_m%sons(1,2)
        blocks1=>blocks_l%sons(2,1)
        blocks3=>blocks_m%sons(2,2)
		! write(*,*)'2b',blocks1%style,blocks2%style,blocks3%style
        call Hmat_add_multiply(blocks3,'-',blocks1,blocks2,h_mat,option,stats,ptree,msh)

        blocks1=>blocks_l%sons(2,2)
        blocks2=>blocks_m%sons(2,1)
        call Hmat_LXM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks2=>blocks_m%sons(2,2)
        call Hmat_LXM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
    else if (blocks_m%style==2) then
		T0 = OMP_get_wtime()
		rank0 = blocks_m%rankmax
		call BF_randomized(blocks_m%pgno,blocks_m%level_butterfly,rank0,option%rankrate,blocks_m,blocks_l,BF_block_MVP_XLM_dat,error,'XLM',option,stats,ptree,msh)
		T1 = OMP_get_wtime()
		stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
		stats%Time_XLUM=stats%Time_XLUM+T1-T0
		stats%XLUM_random_Time(blocks_m%level)=stats%XLUM_random_Time(blocks_m%level)+T1-T0
		stats%XLUM_random_CNT(blocks_m%level)=stats%XLUM_random_CNT(blocks_m%level)+1
    else
		T0 = OMP_get_wtime()
        if (blocks_m%style==1) then
            mm=size(blocks_m%fullmat,1)
            nn=size(blocks_m%fullmat,2)
            do i=1, mm
                ii=blocks_l%ipiv(i)
                if (ii/=i) then
                    !$omp parallel do default(shared) private(j,ctemp)
                    do j=1, nn
                        ctemp=blocks_m%fullmat(i,j)
                        blocks_m%fullmat(i,j)=blocks_m%fullmat(ii,j)
                        blocks_m%fullmat(ii,j)=ctemp
                    enddo
                    !$omp end parallel do
                endif
            enddo
            call trsmf90(blocks_l%fullmat,blocks_m%fullmat,'L','L','N','U',mm,nn)
        else

			write(*,*)'should not come here H_Solve_XLM'

            mm=blocks_m%M
            rank=size(blocks_m%ButterflyU%blocks(1)%matrix,2)
            do i=1, mm
                ii=blocks_l%ipiv(i)
                if (ii/=i) then
                    !$omp parallel do default(shared) private(j,ctemp)
                    do j=1, rank
                        ctemp=blocks_m%ButterflyU%blocks(1)%matrix(i,j)
                        blocks_m%ButterflyU%blocks(1)%matrix(i,j)=blocks_m%ButterflyU%blocks(1)%matrix(ii,j)
                        blocks_m%ButterflyU%blocks(1)%matrix(ii,j)=ctemp
                    enddo
                    !$omp end parallel do
                endif
            enddo
            call trsmf90(blocks_l%fullmat,blocks_m%ButterflyU%blocks(1)%matrix,'L','L','N','U',mm,rank)
        endif
		T1 = OMP_get_wtime()
        stats%Time_XLUM=stats%Time_XLUM+T1-T0

    endif

	! write(*,*) 'out Hmat_LXM'

    return

end subroutine Hmat_LXM


recursive subroutine Hmat_XUM(blocks_u,blocks_m,h_mat,option,stats,ptree,msh)
    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer style(3), mark, style_m
    integer i, j, k, mm, nn, rank
    integer mm_1, mm_2, nn_1, nn_2, rank_1, rank_2, mm_3, nn_3, rank_3
    real (kind=8) T0,T1,error
    type(matrixblock) :: blocks_u,blocks_m
    type(matrixblock), pointer :: blocks1, blocks2, blocks3
	integer rank0

    if (blocks_m%style==4) then
		T0 = OMP_get_wtime()
        if (blocks_u%style/=4) then
			allocate(blocks_u%sons(2,2))
			call BF_split(blocks_u,blocks_u,ptree,stats,msh,option)
        endif
        T1 = OMP_get_wtime()
		stats%Time_Split=stats%Time_Split+T1-T0

        blocks1=>blocks_u%sons(1,1)
        blocks2=>blocks_m%sons(1,1)
        call Hmat_XUM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks2=>blocks_m%sons(2,1)
        call Hmat_XUM(blocks1,blocks2,h_mat,option,stats,ptree,msh)

        blocks1=>blocks_m%sons(1,1)
        blocks2=>blocks_u%sons(1,2)
        blocks3=>blocks_m%sons(1,2)
        call Hmat_add_multiply(blocks3,'-',blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks1=>blocks_m%sons(2,1)
        blocks2=>blocks_u%sons(1,2)
        blocks3=>blocks_m%sons(2,2)
        call Hmat_add_multiply(blocks3,'-',blocks1,blocks2,h_mat,option,stats,ptree,msh)

        blocks1=>blocks_u%sons(2,2)
        blocks2=>blocks_m%sons(1,2)
        call Hmat_XUM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
        blocks2=>blocks_m%sons(2,2)
        call Hmat_XUM(blocks1,blocks2,h_mat,option,stats,ptree,msh)
    else if (blocks_m%style==2) then
		T0 = OMP_get_wtime()
		rank0 = blocks_m%rankmax
		call BF_randomized(blocks_m%pgno,blocks_m%level_butterfly,rank0,option%rankrate,blocks_m,blocks_u,BF_block_MVP_XUM_dat,error,'XUM',option,stats,ptree,msh)
		T1 = OMP_get_wtime()
		stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
		stats%Time_XLUM=stats%Time_XLUM+T1-T0
		stats%XLUM_random_Time(blocks_m%level)=stats%XLUM_random_Time(blocks_m%level)+T1-T0
		stats%XLUM_random_CNT(blocks_m%level)=stats%XLUM_random_CNT(blocks_m%level)+1
    else
		T0 = OMP_get_wtime()
        if (blocks_m%style==1) then
			mm = size(blocks_m%fullmat,1)
			nn = size(blocks_m%fullmat,2)
            call trsmf90(blocks_u%fullmat,blocks_m%fullmat,'R','U','N','N',mm,nn)
        else

			write(*,*)'should not come here Hmat_XUM'
			mm = blocks_m%M
			rank = size(blocks_m%ButterflyV%blocks(1)%matrix,2)
            call trsmf90(blocks_u%fullmat,blocks_m%ButterflyV%blocks(1)%matrix,'L','U','T','N',mm,rank)
        endif
		T1 = OMP_get_wtime()
        stats%Time_XLUM=stats%Time_XLUM+T1-T0
    endif

    return

end subroutine Hmat_XUM


subroutine Hmat_add_multiply_Hblock3(blocks,chara,block1,block2,h_mat,option,stats,ptree,msh)

   implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer group_m, group_n, mm, nn
    integer level_blocks,level_butterfly
    character chara,charatmp
    real*8 T0,T1
    type(matrixblock),target :: block1, block2
    type(matrixblock) :: blocks
	real*8:: error_inout
	type(matrixblock),pointer::block_agent
	integer rank0


	if(block1%style==2)then
		level_butterfly=block1%level_butterfly
		rank0 = block1%rankmax
	endif
	if(block2%style==2)then
		level_butterfly=block2%level_butterfly
		rank0 = block2%rankmax
	endif

	allocate(block_agent)
	call BF_Init_randomized(level_butterfly,rank0,blocks%row_group,blocks%col_group,blocks,block_agent,msh,ptree,option,1)

	! block_agent%row_group = blocks%row_group
	! block_agent%col_group = blocks%col_group
	! block_agent%level = blocks%level
	! block_agent%style=2

	! block_agent%headm=blocks%headm
	! block_agent%M=blocks%M
	! block_agent%headn=blocks%headn
	! block_agent%N=blocks%N

	! block_agent%pgno=blocks%pgno
	! call ComputeParallelIndices(block_agent,block_agent%pgno,ptree,msh,0)
	! call ComputeParallelIndices(block_agent,block_agent%pgno,ptree,msh,1)	! is this needed?


	h_mat%blocks_1=>block1
	h_mat%blocks_2=>block2

	T0=OMP_get_wtime()
	call BF_randomized(block_agent%pgno,level_butterfly,rank0,option%rankrate,block_agent,h_mat,BF_block_MVP_Add_Multiply_dat,error_inout,'Multiply',option,stats,ptree,msh,'m')
	T1=OMP_get_wtime()
	stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
	stats%Time_Multiply=stats%Time_Multiply+T1-T0
	stats%Mul_random_Time(blocks%level)=stats%Mul_random_Time(blocks%level)+T1-T0
	stats%Mul_random_CNT(blocks%level)=stats%Mul_random_CNT(blocks%level)+1

	call Hmat_BF_add(blocks,chara,block_agent,h_mat,option,stats,ptree,msh)

	call Hmat_block_delete(block_agent)
	deallocate(block_agent)

end subroutine Hmat_add_multiply_Hblock3



recursive subroutine Hmat_BF_add(blocks_o,chara,blocks_1,h_mat,option,stats,ptree,msh)
implicit none

type(Hoption)::option
type(Hstat)::stats
type(Hmat)::h_mat
type(proctree)::ptree
type(mesh)::msh

type(matrixblock)::blocks_o
type(matrixblock),target::blocks_1
type(matrixblock),pointer::blocks_1_son,blocks_o_son
character chara
real(kind=8)error
real(kind=8)T0,T1
integer rank0

if(blocks_o%style==4)then
	T0=OMP_get_wtime()
	if (blocks_1%style/=4) then
		allocate(blocks_1%sons(2,2))
		call BF_split(blocks_1,blocks_1,ptree,stats,msh,option)
	endif
    T1 = OMP_get_wtime()
	stats%Time_Split=stats%Time_Split+T1-T0

	blocks_o_son => blocks_o%sons(1,1)
	blocks_1_son => blocks_1%sons(1,1)
	call Hmat_BF_add(blocks_o_son,chara,blocks_1_son,h_mat,option,stats,ptree,msh)
	blocks_o_son => blocks_o%sons(1,2)
	blocks_1_son => blocks_1%sons(1,2)
	call Hmat_BF_add(blocks_o_son,chara,blocks_1_son,h_mat,option,stats,ptree,msh)
	blocks_o_son => blocks_o%sons(2,1)
	blocks_1_son => blocks_1%sons(2,1)
	call Hmat_BF_add(blocks_o_son,chara,blocks_1_son,h_mat,option,stats,ptree,msh)
	blocks_o_son => blocks_o%sons(2,2)
	blocks_1_son => blocks_1%sons(2,2)
	call Hmat_BF_add(blocks_o_son,chara,blocks_1_son,h_mat,option,stats,ptree,msh)

	if (blocks_1%style/=4) then
		call BF_delete(blocks_1%sons(1,1),1)
		call BF_delete(blocks_1%sons(1,2),1)
		call BF_delete(blocks_1%sons(2,1),1)
		call BF_delete(blocks_1%sons(2,2),1)
		deallocate(blocks_1%sons)
	endif
else if(blocks_o%style==2)then
	T0=OMP_get_wtime()
	h_mat%blocks_1=>blocks_1
	rank0 = blocks_o%rankmax
	if(chara=='+')then
		call BF_randomized(blocks_o%pgno,blocks_o%level_butterfly,rank0,option%rankrate,blocks_o,h_mat,BF_block_MVP_Add_Multiply_dat,error,'Add',option,stats,ptree,msh,'a')
	elseif(chara=='-')then
		call BF_randomized(blocks_o%pgno,blocks_o%level_butterfly,rank0,option%rankrate,blocks_o,h_mat,BF_block_MVP_Add_Multiply_dat,error,'Add',option,stats,ptree,msh,'s')
	endif
    T1=OMP_get_wtime()
	stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
    stats%Time_Add_Multiply=stats%Time_Add_Multiply+T1-T0
	stats%Add_random_Time(blocks_o%level)=stats%Add_random_Time(blocks_o%level)+T1-T0
	stats%Add_random_CNT(blocks_o%level)=stats%Add_random_CNT(blocks_o%level)+1
else if(blocks_o%style==1)then
	call Full_add(blocks_o,chara,blocks_1,ptree,stats)
end if

end subroutine Hmat_BF_add


end module BPACK_factor
