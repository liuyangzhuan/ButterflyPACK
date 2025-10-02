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

!> @file BPACK_factor.f90
!> @brief Top-level subroutines for inverting/factorizing H, HOD-LR, HOD-BF and HSS-BF matrices


#include "ButterflyPACK_config.fi"
module BPACK_factor
    use Bplus_factor
    use BPACK_DEFS
    use MISC_Utilities
#ifdef HAVE_OPENMP
    use omp_lib
#endif
    use BPACK_block_sendrecv
    use BPACK_Utilities
    use Bplus_randomizedop
    use BPACK_Solve_Mul
contains


    subroutine BPACK_Factorization(bmat, option, stats, ptree, msh)
        implicit none
        type(Hoption)::option
        type(Hstat)::stats
        type(Bmatrix)::bmat
        type(proctree)::ptree
        type(mesh)::msh

        if (option%precon /= NOPRECON) then
            select case (option%format)
            case (HODLR)
                call HODLR_factorization(bmat%ho_bf, option, stats, ptree, msh)
            case (HMAT,BLR)
                call Hmat_Factorization(bmat%h_mat, option, stats, ptree, msh)
            case (HSS)
                call HSS_factorization(bmat%hss_bf, option, stats, ptree, msh)
            end select
        endif

        if (option%ErrSol == 1) then
            call BPACK_Test_Solve_error(bmat, msh%idxe - msh%idxs + 1, option, ptree, stats)
            if(ptree%MyID==Main_ID)write(*,*) 'RedistV time: ', stats%Time_RedistV
        endif

    end subroutine BPACK_Factorization

    subroutine HODLR_factorization(ho_bf1, option, stats, ptree, msh)

        implicit none

        integer i, j, ii, jj, iii, jjj, index_ij, mm, nn
        integer level, blocks, edge, patch, node, group, level_c, groupm_diag
        integer rank, index_near, m, n, length, flag, itemp
        real T0
        real(kind=8)::rtemp = 0
        real(kind=8) tmpfact, tol_used
        real(kind=8) Memory, Memory_near
        integer, allocatable:: index_old(:), index_new(:)
        integer::block_num, block_num_new, level_butterfly
        integer, allocatable :: ipiv(:)
        integer rowblock, pgno1, pgno2, pgno, ierr, rowblock_inv
        type(matrixblock), pointer::block_o, block_off, block_off1, block_off2
        type(matrixblock)::block_tmp
        real(kind=8) n1, n2, nn1, nn2, flop, norm
        type(Hoption)::option
        type(Hstat)::stats
        type(hobf)::ho_bf1
        type(proctree)::ptree
        type(mesh)::msh
        DT,allocatable::matrixtemp(:,:),UU(:,:),VV(:,:)
        DTR,allocatable::Singular(:)


        if (.not. allocated(stats%rankmax_of_level_global_factor)) allocate (stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel))
        stats%rankmax_of_level_global_factor = 0

        nn1 = MPI_Wtime()

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

        call MPI_barrier(ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Computing block inverse at level Maxlevel+1...'
        level_c = ho_bf1%Maxlevel + 1
        do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
#if HAVE_ZFP
            if(option%use_zfp==1)call ZFP_Decompress(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmatZFP, ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%M, ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%N, tol_used,1)
#endif
            nn = size(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat, 1)
            allocate(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat(nn,nn))
            ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat

            if (option%ErrSol == 0 .and. option%precon == 1) then
                call Bplus_delete(ho_bf1%levels(level_c)%BP(ii))
            else 
#if HAVE_ZFP
            if(option%use_zfp==1)call ZFP_Compress(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%FullmatZFP,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%M,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%N, option%tol_comp,1)
#endif
endif

            allocate(Singular(nn))
            allocate(UU(nn,nn))
            allocate(VV(nn,nn))
            call gesvd_robust(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat, Singular, UU, VV, nn, nn, nn)
            if(Singular(nn)/Singular(1)<option%jitter)then
                norm = Singular(1)
                ! write(*,*)norm*option%jitter,norm
                do iii=1,nn
                    ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat(iii,iii)=ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat(iii,iii)+norm*option%jitter
                enddo
            endif
            deallocate(UU)
            deallocate(VV)
            deallocate(Singular)


#if 1
            allocate (ipiv(nn))
            call getrff90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat, ipiv, flop=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop
            call getrif90(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat, ipiv, flop=flop)
            stats%Flop_Factor = stats%Flop_Factor + flop
            deallocate (ipiv)
#else
            allocate(matrixtemp(n,n))
            matrixtemp = ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat
            call GeneralInverse(nn, nn, matrixtemp, ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat, BPACK_SafeEps, Flops=flop)
            deallocate(matrixtemp)
            stats%Flop_Factor = stats%Flop_Factor + flop
#endif

            !!!!!!! the forward block BP can be deleted if not used in solution phase

            ! write(*,*)fnorm(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,nn,nn)
#if HAVE_ZFP
            if(option%use_zfp==1)then
                call ZFP_Compress(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat,ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%FullmatZFP,ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%M,ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%N, option%tol_comp,0)
                stats%Mem_Direct_inv = stats%Mem_Direct_inv + SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%FullmatZFP%buffer_r)/1024.0d3
#if DAT==0 || DAT==2
                stats%Mem_Direct_inv = stats%Mem_Direct_inv + SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%FullmatZFP%buffer_i)/1024.0d3
#endif
            else
                stats%Mem_Direct_inv = stats%Mem_Direct_inv + SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
            endif
#else
            stats%Mem_Direct_inv = stats%Mem_Direct_inv + SIZEOF(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
#endif
        end do

        call MPI_barrier(ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Computing block inverse at higher levels...'
        do level_c = ho_bf1%Maxlevel, 1, -1

            !!!>***** update the forward off-diagonal block by left multiplication of inverse of diagonal blocks in Z: Z_ij^l -> Z_ii^-1*Z_ij^l
            call MPI_barrier(ptree%Comm, ierr)
            if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'update forward blocks at level:', level_c

            n1 = MPI_Wtime()
            do rowblock_inv = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
            do rowblock = rowblock_inv*2 - 1, rowblock_inv*2

                if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(rowblock)%pgno)) then
                    call Bplus_Sblock_randomized_memfree(ho_bf1, level_c, rowblock, option, stats, ptree, msh)
                    if (option%ErrSol == 0 .and. option%precon == 1) then
                        call Bplus_delete(ho_bf1%levels(level_c)%BP(rowblock))
                    endif    
                    call Bplus_ComputeMemory(ho_bf1%levels(level_c)%BP_inverse_update(rowblock), rtemp, rank)
                    stats%Mem_Sblock = stats%Mem_Sblock + rtemp
                    stats%rankmax_of_level_global_factor(level_c) = max(stats%rankmax_of_level_global_factor(level_c),rank)
                    ! if(level_c==6)then
                    ! call BF_print_size_rank(ho_bf1%levels(level_c)%matrices_block(rowblock),option%tol_comp)
                    ! stop
                    ! end if

                end if
            end do
            end do
            n2 = MPI_Wtime()
            stats%Time_Sblock = stats%Time_Sblock + n2 - n1

            !!!>***** compute the inverse of each block 2x2 submatrices whose two off-diagonal blocks are butterflies
            call MPI_barrier(ptree%Comm, ierr)
            if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'compute block inverse at level:', level_c
            n1 = MPI_Wtime()
            do rowblock = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe

                pgno = ho_bf1%levels(level_c)%BP_inverse(rowblock)%pgno
                if (IOwnPgrp(ptree, pgno)) then
                    call Bplus_ReDistribute_Inplace(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2 - 1), stats, ptree, msh)
                    call Bplus_ReDistribute_Inplace(ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2), stats, ptree, msh)

                    call Bplus_inverse_schur_partitionedinverse(ho_bf1, level_c, rowblock, option, stats, ptree, msh)
                    call Bplus_ComputeMemory(ho_bf1%levels(level_c)%BP_inverse_schur(rowblock), rtemp,rank)
                    stats%Mem_SMW = stats%Mem_SMW + rtemp
                    stats%rankmax_of_level_global_factor(level_c) = max(stats%rankmax_of_level_global_factor(level_c),rank)

                endif
            end do
            n2 = MPI_Wtime()
            stats%Time_Inv = stats%Time_Inv + n2 - n1
        end do

        nn2 = MPI_Wtime()
        stats%Time_Factor = nn2 - nn1
        call MPI_ALLREDUCE(MPI_IN_PLACE, stats%rankmax_of_level_global_factor(0:ho_bf1%Maxlevel), ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
        call MPI_ALLREDUCE(stats%Time_Sblock, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
       if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'computing updated forward block time:', rtemp, 'Seconds'
        call MPI_ALLREDUCE(stats%Time_Inv, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'computing inverse block time:', rtemp, 'Seconds'
        call MPI_ALLREDUCE(stats%Time_random(1), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Init:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(2), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_MVP:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(3), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Reconstruct:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(4), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Onesub:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(5), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_OneBlock_XX:', rtemp
        call MPI_ALLREDUCE(stats%Time_SMW, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_SMW:', rtemp
        call MPI_ALLREDUCE(stats%Time_PartialUpdate, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_PartialUpdate:', rtemp
        call MPI_ALLREDUCE(stats%Time_RedistB, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_RedistB:', rtemp
        call MPI_ALLREDUCE(stats%Time_RedistV, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_RedistV:', rtemp
        call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', rtemp
        call MPI_ALLREDUCE(time_tmp1, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp1', rtemp
        call MPI_ALLREDUCE(time_tmp2, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp2', rtemp
        call MPI_ALLREDUCE(time_tmp3, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp3', rtemp
        call MPI_ALLREDUCE(time_tmp4, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp4', rtemp
        call MPI_ALLREDUCE(time_tmp5, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp5', rtemp
        call MPI_ALLREDUCE(stats%Flop_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A21,Es14.2)') 'Factorization flops:', rtemp

        stats%Mem_Factor = stats%Mem_SMW + stats%Mem_Sblock + stats%Mem_Direct_inv

        if (option%ErrSol == 0 .and. option%precon == 1) then        
            call LogMemory(stats, stats%Mem_Factor-stats%Mem_Fill)
        else 
        call LogMemory(stats, stats%Mem_Factor)
        endif        

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
        call MPI_ALLREDUCE(stats%Mem_SMW, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly inverse blocks'
        call MPI_ALLREDUCE(stats%Mem_Sblock, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly Sblocks'
        call MPI_ALLREDUCE(stats%Mem_Direct_inv, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for direct inverse blocks'
        call MPI_ALLREDUCE(stats%Mem_int_vec, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for storing intermidiate vectors'
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

        return

    end subroutine HODLR_factorization

    subroutine HSS_factorization(hss_bf1, option, stats, ptree, msh)

        implicit none

        integer i, j, ii, jj, iii, jjj, index_ij, mm, nn
        integer level, blocks, edge, patch, node, group, level_c, groupm_diag
        integer rank, index_near, m, n, length, flag, itemp
        real T0
        real(kind=8)::rtemp = 0
        real(kind=8) tmpfact
        real(kind=8) Memory, Memory_near
        integer, allocatable:: index_old(:), index_new(:)
        integer::block_num, block_num_new, level_butterfly
        integer, allocatable :: ipiv(:)
        integer rowblock, pgno1, pgno2, pgno, ierr, rowblock_inv
        type(matrixblock), pointer::block_o, block_off, block_off1, block_off2
        type(matrixblock)::block_tmp
        real(kind=8) n1, n2, nn1, nn2, flop
        type(Hoption)::option
        type(Hstat)::stats
        type(hssbf)::hss_bf1
        type(proctree)::ptree
        type(mesh)::msh

        if(.not. allocated(stats%rankmax_of_level_global_factor))allocate (stats%rankmax_of_level_global_factor(0:hss_bf1%Maxlevel))
        stats%rankmax_of_level_global_factor = 0

        nn1 = MPI_Wtime()

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

        call MPI_barrier(ptree%Comm, ierr)

        call Bplus_copy(hss_bf1%BP, hss_bf1%BP_inverse)
        call LogMemory(stats, stats%Mem_Fill)

        if (option%ErrSol == 0 .and. option%precon == 1) then
            call Bplus_delete(hss_bf1%BP)
            call LogMemory(stats, -stats%Mem_Fill)
        endif    

        call Bplus_inverse_schur_partitionedinverse_hss(hss_bf1%BP_inverse, option, stats, ptree, msh)
        call Bplus_ComputeMemory(hss_bf1%BP_inverse, rtemp,rank)
        stats%rankmax_of_level_global_factor(0)=rank
        stats%Mem_Factor = rtemp
        call LogMemory(stats, stats%Mem_Factor-stats%Mem_Fill)

        nn2 = MPI_Wtime()
        stats%Time_Inv = nn2 - nn1
        stats%Time_Factor = stats%Time_Inv
        call MPI_ALLREDUCE(MPI_IN_PLACE, stats%rankmax_of_level_global_factor(0:hss_bf1%Maxlevel), hss_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
        call MPI_ALLREDUCE(stats%Time_Inv, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'computing inverse block time:', rtemp, 'Seconds'
        call MPI_ALLREDUCE(stats%Time_random(1), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Init:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(2), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_MVP:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(3), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Reconstruct:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(4), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Onesub:', rtemp
        call MPI_ALLREDUCE(stats%Time_random(5), rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_OneBlock_XX:', rtemp
        call MPI_ALLREDUCE(stats%Time_SMW, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_SMW:', rtemp
        call MPI_ALLREDUCE(stats%Time_PartialUpdate, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_PartialUpdate:', rtemp
        call MPI_ALLREDUCE(stats%Time_RedistB, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_RedistB:', rtemp
        call MPI_ALLREDUCE(stats%Time_RedistV, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_RedistV:', rtemp
        call MPI_ALLREDUCE(stats%Time_Split, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Split:', rtemp
        call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', rtemp
        call MPI_ALLREDUCE(time_tmp1, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp1', rtemp
        call MPI_ALLREDUCE(time_tmp2, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp2', rtemp
        call MPI_ALLREDUCE(time_tmp3, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp3', rtemp
        call MPI_ALLREDUCE(time_tmp4, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp4', rtemp
        call MPI_ALLREDUCE(time_tmp5, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp5', rtemp
        call MPI_ALLREDUCE(stats%Flop_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A21,Es14.2)') 'Factorization flops:', rtemp

        

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
        call MPI_ALLREDUCE(stats%Mem_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly inverse blocks'
        call MPI_ALLREDUCE(stats%Mem_int_vec, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for storing intermidiate vectors'
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

        return

    end subroutine HSS_factorization

    subroutine Hmat_Factorization(h_mat, option, stats, ptree, msh)
        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer level, flag, num_blocks, level_butterfly
        integer Primary_block, kk, i, j, k, ij, intemp, ierr
        integer systime0(8), systime1(8)
        character(len=10) :: date1, time1, zone1
        real*8 rtemp1, rtemp2, rtemp
        real*8 nn1, nn2
        real*8 Memory, Memory_near
        type(matrixblock), pointer :: block
        integer mypgno
        type(nod), pointer::cur
        class(*), pointer::ptr
        type(matrixblock), pointer :: blocks, blocks_r, blocks_s, blocks_mm, blocks_uu, blocks_ll
        integer nprow, npcol, myrow, mycol, iproc, myi, jproc, myj, jproc1, myj1, iproc1, myi1, recv,send, send_ID
        real(kind=8) T0, T1, T3, T4
        integer             :: tid, nthreads
        double precision    :: t_start, t_stop
        double precision, allocatable :: t_total(:)

        nthreads=1
#ifdef HAVE_OPENMP
        nthreads = omp_get_max_threads()
#endif

        allocate(h_mat%blocks_1(nthreads))
        allocate(h_mat%blocks_2(nthreads))

        allocate(t_total(nthreads))
        t_total = 0.0d0

        if (.not. allocated(stats%rankmax_of_level_global_factor)) allocate (stats%rankmax_of_level_global_factor(0:h_mat%Maxlevel))
        stats%rankmax_of_level_global_factor = 0

#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
      !$omp parallel
      !$omp single
#endif
#endif
        nn1 = MPI_Wtime()

        if(option%ILU==1)then
            level=h_mat%Maxlevel
            cur => h_mat%lstblks(level)%head
            do i = 1, h_mat%lstblks(level)%num_nods
                select type (ptr=>cur%item)
                type is (block_ptr)
                    if(ptr%ptr%row_group==ptr%ptr%col_group)then
                        call Full_LU(ptr%ptr, option, stats)
                    endif
                end select
                cur => cur%next
            enddo
        else
            ! pack all blocks at Dist_level
            do j = 1, h_mat%myAcols
            do i = 1, h_mat%myArows
                blocks => h_mat%Local_blocks(j, i)
                ! call pack_all_blocks_one_node(blocks, msh,option)
                mypgno = blocks%pgno
            enddo
            enddo

            allocate (h_mat%Computing_matricesblock_m(1, 1))
            allocate (h_mat%Computing_matricesblock_l(max(1,h_mat%myAcols), max(1,h_mat%myArows)))
            allocate (h_mat%Computing_matricesblock_u(max(1,h_mat%myAcols), max(1,h_mat%myArows)))


            call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
            num_blocks = 2**h_mat%Dist_level
            do kk=1,num_blocks
                if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "starting panel ", kk

                T3 = MPI_Wtime()
                call MPI_verbose_barrier('   --diagonal factorization', ptree,option)
                T4 = MPI_Wtime()
                stats%Time_idle = stats%Time_idle + T4 - T3

                call g2l(kk, num_blocks, nprow, 1, iproc, myi)
                call g2l(kk, num_blocks, npcol, 1, jproc, myj)
                if(h_mat%myArows>0 .and. h_mat%myAcols>0)send_ID = blacs_pnum_wp(nprow,npcol, iproc, jproc)
                send=0
                recv=0
                if(iproc==myrow .and. jproc==mycol)then
                    blocks => h_mat%Local_blocks(myj, myi)
                    ! call unpack_all_blocks_one_node(blocks, h_mat%Maxlevel, ptree, msh, mypgno)
                    call Hmat_LU(blocks, h_mat, option, stats, ptree, msh)
                    call pack_all_blocks_one_node(blocks, msh, option)
                    send=1
                endif

                T3 = MPI_Wtime()
                call MPI_verbose_barrier('   --sending the digonal block', ptree, option)
                T4 = MPI_Wtime()
                stats%Time_idle = stats%Time_idle + T4 - T3

                T3 = MPI_Wtime()
                do j=kk+1,num_blocks
                    call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                    if(iproc==myrow .and. jproc1==mycol)then
                    recv=1
                    endif
                enddo
                do i=kk+1,num_blocks
                    call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                    if(iproc1==myrow .and. jproc==mycol)then
                    recv=1
                    endif
                enddo
                if (recv == 1) blocks_r => h_mat%Computing_matricesblock_m(1, 1)
                if (send == 1) blocks_s => h_mat%Local_blocks(myj, myi)
                call blocks_partial_bcast(blocks_s, blocks_r, send, recv, send_ID, msh, ptree, option)
                T4 = MPI_Wtime()
                stats%Time_Comm = stats%Time_Comm + T4 - T3

                T3 = MPI_Wtime()
                call MPI_verbose_barrier('   --L and U panel factorization', ptree, option)
                T4 = MPI_Wtime()
                stats%Time_idle = stats%Time_idle + T4 - T3

                if (recv == 1)then
                    call unpack_all_blocks_one_node(blocks_r, h_mat%Maxlevel, ptree, msh, mypgno, option)
                    Memory=0
                    call Hmat_block_ComputeMemory(blocks_r, Memory)
                    call LogMemory(stats, Memory)
                endif

#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
                !$omp taskloop default(shared) private(j,jproc1, myj1,blocks_uu,blocks_mm)
#endif
#endif
                do j=kk+1,num_blocks
                    call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                    if(iproc==myrow .and. jproc1==mycol)then
                        blocks_uu => h_mat%Local_blocks(myj1, myi)
                        blocks_mm => h_mat%Computing_matricesblock_m(1, 1)
                        call Hmat_LXM(blocks_mm, blocks_uu, h_mat, option, stats, ptree, msh)
                        call pack_all_blocks_one_node(blocks_uu, msh, option)
                    endif
                enddo
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
                !$omp end taskloop
#endif
#endif
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
                !$omp taskloop default(shared) private(i,iproc1, myi1,blocks_ll,blocks_mm)
#endif
#endif
                do i=kk+1,num_blocks
                    call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                    if(iproc1==myrow .and. jproc==mycol)then
                        blocks_ll => h_mat%Local_blocks(myj, myi1)
                        blocks_mm => h_mat%Computing_matricesblock_m(1, 1)
                        call Hmat_XUM(blocks_mm, blocks_ll, h_mat, option, stats, ptree, msh)
                        call pack_all_blocks_one_node(blocks_ll, msh, option)
                    endif
                enddo
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
                !$omp end taskloop
#endif
#endif
                if (recv == 1)then
                    Memory=0
                    call Hmat_block_ComputeMemory(blocks_r, Memory)
                    call LogMemory(stats, -Memory)
                    call Hmat_block_delete(blocks_r)
                endif

                T3 = MPI_Wtime()
                call MPI_verbose_barrier('   --sending L and U panels', ptree, option)
                T4 = MPI_Wtime()
                stats%Time_idle = stats%Time_idle + T4 - T3

                T3 = MPI_Wtime()
                do j=kk+1,num_blocks
                    send=0
                    recv=0
                    call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                    if(h_mat%myArows>0 .and. h_mat%myAcols>0)send_ID = blacs_pnum_wp(nprow,npcol, iproc, jproc1)
                    if(iproc==myrow .and. jproc1==mycol)then
                        send=1
                    endif
                    do i=kk+1,num_blocks
                        call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                        if(iproc1==myrow .and. jproc1==mycol)then
                            recv=1
                        endif
                    enddo
                    if (recv == 1) blocks_r => h_mat%Computing_matricesblock_m(1, 1)
                    if (send == 1) blocks_s => h_mat%Local_blocks(myj1, myi)
                    call blocks_partial_bcast(blocks_s, blocks_r, send, recv, send_ID, msh, ptree, option)
                    if (recv == 1)then
                        do i=kk+1,num_blocks
                            call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                            if(iproc1==myrow .and. jproc1==mycol)then
                                blocks_uu => h_mat%Computing_matricesblock_u(myj1, myi1)
                                call Hmat_block_copy_MPIdata(blocks_uu, blocks_r, msh)
                                call unpack_all_blocks_one_node(blocks_uu, h_mat%Maxlevel, ptree, msh, mypgno, option)
                                Memory=0
                                call Hmat_block_ComputeMemory(blocks_uu, Memory)
                                call LogMemory(stats, Memory)
                            endif
                        enddo
                    endif
                    if (recv == 1)call Hmat_block_delete(blocks_r)
                enddo

                do i=kk+1,num_blocks
                    send=0
                    recv=0
                    call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                    if(h_mat%myArows>0 .and. h_mat%myAcols>0)send_ID = blacs_pnum_wp(nprow,npcol, iproc1, jproc)
                    if(iproc1==myrow .and. jproc==mycol)then
                        send=1
                    endif
                    do j=kk+1,num_blocks
                        call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                        if(iproc1==myrow .and. jproc1==mycol)then
                            recv=1
                        endif
                    enddo
                    if (recv == 1) blocks_r => h_mat%Computing_matricesblock_m(1, 1)
                    if (send == 1) blocks_s => h_mat%Local_blocks(myj, myi1)
                    call blocks_partial_bcast(blocks_s, blocks_r, send, recv, send_ID, msh, ptree, option)
                    if (recv == 1)then
                        do j=kk+1,num_blocks
                            call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                            if(iproc1==myrow .and. jproc1==mycol)then
                                blocks_ll => h_mat%Computing_matricesblock_l(myj1, myi1)
                                call Hmat_block_copy_MPIdata(blocks_ll, blocks_r, msh)
                                call unpack_all_blocks_one_node(blocks_ll, h_mat%Maxlevel, ptree, msh, mypgno, option)
                                Memory=0
                                call Hmat_block_ComputeMemory(blocks_ll, Memory)
                                call LogMemory(stats, Memory)
                            endif
                        enddo
                    endif
                    if (recv == 1)call Hmat_block_delete(blocks_r)
                enddo
                T4 = MPI_Wtime()
                stats%Time_Comm = stats%Time_Comm + T4 - T3

                T3 = MPI_Wtime()
                call MPI_verbose_barrier('   --Schur updates', ptree, option)
                T4 = MPI_Wtime()
                stats%Time_idle = stats%Time_idle + T4 - T3


                if(kk+1<=num_blocks)then
                    T3 = MPI_Wtime()
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
             !$omp taskloop default(shared) private(i,j,iproc1, myi1,jproc1, myj1,blocks_ll,blocks_uu,blocks_mm,Memory,t_start, t_stop, tid)
#endif
#endif
                    do ij = 1, (num_blocks-kk)*(num_blocks-kk)
                        t_start = omp_get_wtime()
                        j = (ij - 1)/(num_blocks-kk) + 1
                        i = mod(ij - 1, (num_blocks-kk)) + 1
                        i = i + kk
                        j = j + kk
                        tid = omp_get_thread_num()
                        call g2l(i, num_blocks, nprow, 1, iproc1, myi1)
                        call g2l(j, num_blocks, npcol, 1, jproc1, myj1)
                        if(iproc1==myrow .and. jproc1==mycol)then
                            blocks_ll => h_mat%Computing_matricesblock_l(myj1, myi1)
                            blocks_uu => h_mat%Computing_matricesblock_u(myj1, myi1)
                            blocks_mm => h_mat%Local_blocks(myj1, myi1)
                            call Hmat_add_multiply(blocks_mm, '-', blocks_ll, blocks_uu, h_mat, option, stats, ptree, msh)
                            Memory=0
                            call Hmat_block_ComputeMemory(blocks_ll, Memory)
                            call LogMemory(stats, -Memory)
                            Memory=0
                            call Hmat_block_ComputeMemory(blocks_uu, -Memory)
                            call LogMemory(stats, -Memory)
                            call Hmat_block_delete(blocks_ll)
                            call Hmat_block_delete(blocks_uu)
                        endif
                        t_stop  = omp_get_wtime()
                        t_total(tid+1) = t_total(tid+1) + (t_stop - t_start)
                    enddo
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
                   !$omp end taskloop
#endif
#endif
                    T4 = MPI_Wtime()
                    ! time_tmp = time_tmp + T4 - T3
                    ! write(*,*)"total_time: ", T4 - T3
                    ! do tid = 1, nthreads
                    !     write(*,'(A,I3,2A,F8.4,A)') 'Thread ', tid-1, ': ', &
                    !         'spent ', t_total(tid), '  seconds in tasks'
                    ! end do

                endif

            enddo

            deallocate(h_mat%Computing_matricesblock_m)
            deallocate(h_mat%Computing_matricesblock_l)
            deallocate(h_mat%Computing_matricesblock_u)

            ! unpack all blocks at Dist_level
            do j = 1, h_mat%myAcols
                do i = 1, h_mat%myArows
                    blocks => h_mat%Local_blocks(j, i)
                    call unpack_all_blocks_one_node(blocks, h_mat%Maxlevel, ptree, msh, mypgno,option)
                    call Hmat_block_ComputeMemory(blocks, stats%Mem_Factor)
                enddo
            enddo
        endif

        nn2 = MPI_Wtime()
        stats%Time_Factor = nn2 - nn1
#ifdef HAVE_TOPLEVEL_OPENMP
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
#endif

        deallocate(t_total)
        deallocate(h_mat%blocks_1)
        deallocate(h_mat%blocks_2)

        call MPI_ALLREDUCE(MPI_IN_PLACE, stats%rankmax_of_level_global_factor(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)

        call MPI_ALLREDUCE(stats%Time_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Factor:', rtemp
        call MPI_ALLREDUCE(stats%Time_Direct_LU, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Direct_LU:', rtemp
        call MPI_ALLREDUCE(stats%Time_Add_Multiply, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Add_Multiply:', rtemp
        call MPI_ALLREDUCE(stats%Time_Multiply, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Multiply:', rtemp
        call MPI_ALLREDUCE(stats%Time_XLUM, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_XLUM:', rtemp
        call MPI_ALLREDUCE(stats%Time_Split, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Split:', rtemp
        call MPI_ALLREDUCE(stats%Time_Comm, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Comm:', rtemp
        call MPI_ALLREDUCE(stats%Time_Idle, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) '     Time_Idle:', rtemp
        call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', rtemp


        call MPI_allreduce(MPI_IN_PLACE, stats%Add_random_CNT(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_INTEGER, MPI_sum, ptree%Comm, ierr)
        call MPI_allreduce(MPI_IN_PLACE, stats%Mul_random_CNT(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_INTEGER, MPI_sum, ptree%Comm, ierr)
        call MPI_allreduce(MPI_IN_PLACE, stats%XLUM_random_CNT(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_INTEGER, MPI_sum, ptree%Comm, ierr)
        call MPI_allreduce(MPI_IN_PLACE, stats%Add_random_Time(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_DOUBLE_PRECISION, MPI_max, ptree%Comm, ierr)
        call MPI_allreduce(MPI_IN_PLACE, stats%Mul_random_Time(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_DOUBLE_PRECISION, MPI_max, ptree%Comm, ierr)
        call MPI_allreduce(MPI_IN_PLACE, stats%XLUM_random_Time(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_DOUBLE_PRECISION, MPI_max, ptree%Comm, ierr)

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            write (*, *) ''
            write (*, *) 'Randomized OPs: Count and Time'

            do level = 0, h_mat%Maxlevel
                if (level > option%LRlevel) then
                    level_butterfly = 0 ! low rank below LRlevel
                else
                    level_butterfly = h_mat%Maxlevel - level   ! butterfly
                endif

                if (stats%Add_random_CNT(level) + stats%Mul_random_CNT(level) + stats%XLUM_random_CNT(level) /= 0) then
                write (*, '(A7,I5,A17,I5,A7,I8,Es10.2,A7,I8,Es10.2,A12,I8,Es10.2)') " level:", level, 'level_butterfly:', level_butterfly, 'add:', stats%Add_random_CNT(level), stats%Add_random_Time(level), 'mul:', stats%Mul_random_CNT(level), stats%Mul_random_Time(level), 'XLUM:', stats%XLUM_random_CNT(level), stats%XLUM_random_time(level)
                endif
            enddo
            ! write(*,*)'max inverse butterfly rank:', butterflyrank_inverse
        endif

        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            write (*, *) ''
            write (*, *) 'Unpacking all blocks...'
        endif

        call LogMemory(stats, stats%Mem_Factor-stats%Mem_Fill)

        call MPI_verbose_barrier('after printing', ptree, option)
        if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            write (*, *) 'Unpacking finished'
            write (*, *) ''
        endif

        return

    end subroutine Hmat_Factorization


    recursive subroutine Hmat_LU(blocks, h_mat, option, stats, ptree, msh)

        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer level, blocks_son(2, 2), id
        integer i, j, k, ii, jj, kk, flag
        real*8 rtemp1, rtemp2
        real*8 T0, T1
        type(matrixblock):: blocks
        type(matrixblock), pointer :: block_son1, block_son2, block_son3

        if (blocks%style == 4) then
            block_son1 => blocks%sons(1, 1)
            call Hmat_LU(block_son1, h_mat, option, stats, ptree, msh)

            if (option%ILU == 0) then
                block_son1 => blocks%sons(1, 1)
                block_son2 => blocks%sons(2, 1)
                call Hmat_XUM(block_son1, block_son2, h_mat, option, stats, ptree, msh)

                block_son1 => blocks%sons(1, 1)
                block_son2 => blocks%sons(1, 2)
                call Hmat_LXM(block_son1, block_son2, h_mat, option, stats, ptree, msh)

                block_son1 => blocks%sons(2, 1)
                block_son2 => blocks%sons(1, 2)
                block_son3 => blocks%sons(2, 2)
                call Hmat_add_multiply(block_son3, '-', block_son1, block_son2, h_mat, option, stats, ptree, msh)
            endif

            block_son1 => blocks%sons(2, 2)
            call Hmat_LU(block_son1, h_mat, option, stats, ptree, msh)
        else
            call Full_LU(blocks, option, stats)
        endif

        return

    end subroutine Hmat_LU

    recursive subroutine Hmat_add_multiply(block3, chara, block1, block2, h_mat, option, stats, ptree, msh)

        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer rank0,tid
        integer level_butterfly
        integer style(3), mark(3)
        integer ii, jj, i, j, k, level, mm, nn, rank, level_blocks, kk, found_flag, group_m, group_n, m, n
        character chara
        DT, allocatable:: matrixtemp(:, :)
        real*8 T0, T1, error
        type(matrixblock), pointer :: matrices_tmpblock
        type(matrixblock), target :: block2, block1, block3
        type(matrixblock), pointer :: block1_son, block2_son, block3_son

        tid = 0
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num()
#endif

        style(3) = block3%style
        style(1) = block1%style
        style(2) = block2%style

        level_blocks = block3%level

        if ((style(1) ==2 .or. style(2) ==2)) then
            T0 = MPI_Wtime()
            call Hmat_add_multiply_Hblock3(block3, chara, block1, block2, h_mat, option, stats, ptree, msh)
            T1 = MPI_Wtime()
        else if (style(3) == 4) then   !!! modified by Yang Liu, hybrid butterfly-LR treatment
            T0 = MPI_Wtime()
            if (style(1) /= 4) then
                allocate (block1%sons(2, 2))
                call BF_split(block1, block1, ptree, stats, msh, option)
            endif
            if (style(2) /= 4) then
                allocate (block2%sons(2, 2))
                call BF_split(block2, block2, ptree, stats, msh, option)
            endif
            if (style(3) /= 4) then
                allocate (block3%sons(2, 2))
                call BF_split(block3, block3, ptree, stats, msh, option)
            endif
            T1 = MPI_Wtime()
            stats%Time_split = stats%Time_split + T1 - T0

            block1_son => block1%sons(1, 1)
            block2_son => block2%sons(1, 1)
            block3_son => block3%sons(1, 1)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(1, 2)
            block2_son => block2%sons(2, 1)
            block3_son => block3%sons(1, 1)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(1, 1)
            block2_son => block2%sons(1, 2)
            block3_son => block3%sons(1, 2)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(1, 2)
            block2_son => block2%sons(2, 2)
            block3_son => block3%sons(1, 2)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(2, 1)
            block2_son => block2%sons(1, 1)
            block3_son => block3%sons(2, 1)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(2, 2)
            block2_son => block2%sons(2, 1)
            block3_son => block3%sons(2, 1)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(2, 1)
            block2_son => block2%sons(1, 2)
            block3_son => block3%sons(2, 2)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
            block1_son => block1%sons(2, 2)
            block2_son => block2%sons(2, 2)
            block3_son => block3%sons(2, 2)
            call Hmat_add_multiply(block3_son, chara, block1_son, block2_son, h_mat, option, stats, ptree, msh)
        elseif (style(3) == 1 .or. (style(1) == 1 .and. style(2) == 1 .and. block3%level_butterfly==0)) then
        ! elseif (style(3) == 1) then
            call Full_add_multiply(block3, chara, block1, block2, h_mat, option, stats, ptree, msh)
        elseif (style(3) == 2) then
            T0 = MPI_Wtime()

            h_mat%blocks_1(tid+1)%ptr => block1
            h_mat%blocks_2(tid+1)%ptr => block2
            rank0 = block3%rankmax
            call BF_randomized(block3%pgno, block3%level_butterfly, rank0, option%rankrate, block3, h_mat, BF_block_MVP_Add_Multiply_dat, error, 'Add_Multiply', option, stats, ptree, msh, operand1=chara)
            T1 = MPI_Wtime()
            stats%rankmax_of_level_global_factor(block3%level)=max(stats%rankmax_of_level_global_factor(block3%level),block3%rankmax)
            stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
            stats%Time_Add_Multiply = stats%Time_Add_Multiply + T1 - T0
            ! time_tmp = time_tmp + T1 - T0
            stats%Add_random_Time(level_blocks) = stats%Add_random_Time(level_blocks) + T1 - T0
            stats%Add_random_CNT(level_blocks) = stats%Add_random_CNT(level_blocks) + 1

        endif

#if 0
        if (dumping == 0) then
        if (style(3) == 2) then
            group_m = block3%row_group
            group_n = block3%col_group
            write (*, *) group_m, group_n, m, n
            if (group_m == 13 .and. group_n == 11) then
                dumping = 1
                m = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                n = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                allocate (matrixtemp(m, n))
                call butterfly_block_fullextract(block3, matrixtemp)
                write (*, *) 'dumping:', group_m, group_n, m, n
                do ii = 1, m
                do jj = 1, n
                    write (777, *) dble(matrixtemp(ii, jj)), aimag(matrixtemp(ii, jj))
                enddo
                enddo
                deallocate (matrixtemp)
                write (*, *) 'dumping done'
                ! stop
            endif
        endif
        endif
#endif

        return

    end subroutine Hmat_add_multiply

    recursive subroutine Hmat_LXM(blocks_l, blocks_m, h_mat, option, stats, ptree, msh)
        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer style(3), mark, style_m
        integer i, j, k, ii
        integer mm, nn, rank
        real(kind=8) T0, T1, T2, T3, error, tol_used
        type(matrixblock), pointer :: blocks1, blocks2, blocks3
        type(matrixblock) :: blocks_l, blocks_m
        DT:: ctemp
        integer rank0

        if (blocks_m%style == 4) then
            T0 = MPI_Wtime()
            if (blocks_l%style /= 4) then
                allocate (blocks_l%sons(2, 2))
                call BF_split(blocks_l, blocks_l, ptree, stats, msh, option)
            endif
            T1 = MPI_Wtime()
            stats%Time_Split = stats%Time_Split + T1 - T0

            blocks1 => blocks_l%sons(1, 1)
            blocks2 => blocks_m%sons(1, 1)
            call Hmat_LXM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks2 => blocks_m%sons(1, 2)
            call Hmat_LXM(blocks1, blocks2, h_mat, option, stats, ptree, msh)

            blocks2 => blocks_m%sons(1, 1)
            blocks1 => blocks_l%sons(2, 1)
            blocks3 => blocks_m%sons(2, 1)
            ! write(*,*)'1b',blocks1%style,blocks2%style,blocks3%style
            call Hmat_add_multiply(blocks3, '-', blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks2 => blocks_m%sons(1, 2)
            blocks1 => blocks_l%sons(2, 1)
            blocks3 => blocks_m%sons(2, 2)
            ! write(*,*)'2b',blocks1%style,blocks2%style,blocks3%style
            call Hmat_add_multiply(blocks3, '-', blocks1, blocks2, h_mat, option, stats, ptree, msh)

            blocks1 => blocks_l%sons(2, 2)
            blocks2 => blocks_m%sons(2, 1)
            call Hmat_LXM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks2 => blocks_m%sons(2, 2)
            call Hmat_LXM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
        else if (blocks_m%style == 2) then
            T0 = MPI_Wtime()
            if(blocks_m%level_butterfly==0)then
                rank = size(blocks_m%butterflyU%blocks(1)%matrix,2)
                call Hmat_Lsolve(blocks_l, 'N', blocks_l%headm, rank, blocks_m%butterflyU%blocks(1)%matrix, blocks_l%M, ptree, stats)
            else
                T2 = MPI_Wtime()
                rank0 = blocks_m%rankmax
                call BF_randomized(blocks_m%pgno, blocks_m%level_butterfly, rank0, option%rankrate, blocks_m, blocks_l, BF_block_MVP_XLM_dat, error, 'XLM', option, stats, ptree, msh)
                T3 = MPI_Wtime()
                stats%XLUM_random_Time(blocks_m%level) = stats%XLUM_random_Time(blocks_m%level) + T3 - T2
                stats%XLUM_random_CNT(blocks_m%level) = stats%XLUM_random_CNT(blocks_m%level) + 1
            endif
            T1 = MPI_Wtime()
            stats%rankmax_of_level_global_factor(blocks_m%level)=max(stats%rankmax_of_level_global_factor(blocks_m%level),blocks_m%rankmax)
            stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
            stats%Time_XLUM = stats%Time_XLUM + T1 - T0
        else
            T0 = MPI_Wtime()
            if (blocks_m%style == 1) then
#if HAVE_ZFP
            if(option%use_zfp>=1)then
                 call ZFP_Decompress(blocks_m%fullmat,blocks_m%FullmatZFP,blocks_m%M,blocks_m%N,tol_used,0)
            endif
            if(option%use_zfp==1)then
                 call ZFP_Decompress(blocks_l%fullmat,blocks_l%FullmatZFP,blocks_l%M,blocks_l%N,tol_used,1)
            endif
#endif
                mm = size(blocks_m%fullmat, 1)
                nn = size(blocks_m%fullmat, 2)
                do i = 1, mm
                    ii = blocks_l%ipiv(i)
                    if (ii /= i) then
                        ! parallel do default(shared) private(j,ctemp)
                        do j = 1, nn
                            ctemp = blocks_m%fullmat(i, j)
                            blocks_m%fullmat(i, j) = blocks_m%fullmat(ii, j)
                            blocks_m%fullmat(ii, j) = ctemp
                        enddo
                        ! end parallel do
                    endif
                enddo
                call trsmf90(blocks_l%fullmat, blocks_m%fullmat, 'L', 'L', 'N', 'U', mm, nn)
#if HAVE_ZFP
                if(option%use_zfp>=1)then
                    call ZFP_Compress(blocks_m%fullmat,blocks_m%FullmatZFP,blocks_m%M,blocks_m%N,option%tol_comp,0)
                endif
                if(option%use_zfp==1)then
                    call ZFP_Compress(blocks_l%fullmat,blocks_l%FullmatZFP,blocks_l%M,blocks_l%N,option%tol_comp,1)
                endif

#endif

            else

                write (*, *) 'should not come here H_Solve_XLM'
                stop

                mm = blocks_m%M
                rank = size(blocks_m%ButterflyU%blocks(1)%matrix, 2)
                do i = 1, mm
                    ii = blocks_l%ipiv(i)
                    if (ii /= i) then
                        ! parallel do default(shared) private(j,ctemp)
                        do j = 1, rank
                            ctemp = blocks_m%ButterflyU%blocks(1)%matrix(i, j)
                            blocks_m%ButterflyU%blocks(1)%matrix(i, j) = blocks_m%ButterflyU%blocks(1)%matrix(ii, j)
                            blocks_m%ButterflyU%blocks(1)%matrix(ii, j) = ctemp
                        enddo
                        ! end parallel do
                    endif
                enddo
                call trsmf90(blocks_l%fullmat, blocks_m%ButterflyU%blocks(1)%matrix, 'L', 'L', 'N', 'U', mm, rank)
            endif
            T1 = MPI_Wtime()
            stats%Time_XLUM = stats%Time_XLUM + T1 - T0

        endif

        ! write(*,*) 'out Hmat_LXM'

        return

    end subroutine Hmat_LXM

    recursive subroutine Hmat_XUM(blocks_u, blocks_m, h_mat, option, stats, ptree, msh)
        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer style(3), mark, style_m
        integer i, j, k, mm, nn, rank
        integer mm_1, mm_2, nn_1, nn_2, rank_1, rank_2, mm_3, nn_3, rank_3
        real(kind=8) T0, T1, T2, T3, error, tol_used
        type(matrixblock) :: blocks_u, blocks_m
        type(matrixblock), pointer :: blocks1, blocks2, blocks3
        integer rank0

        if (blocks_m%style == 4) then
            T0 = MPI_Wtime()
            if (blocks_u%style /= 4) then
                allocate (blocks_u%sons(2, 2))
                call BF_split(blocks_u, blocks_u, ptree, stats, msh, option)
            endif
            T1 = MPI_Wtime()
            stats%Time_Split = stats%Time_Split + T1 - T0

            blocks1 => blocks_u%sons(1, 1)
            blocks2 => blocks_m%sons(1, 1)
            call Hmat_XUM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks2 => blocks_m%sons(2, 1)
            call Hmat_XUM(blocks1, blocks2, h_mat, option, stats, ptree, msh)

            blocks1 => blocks_m%sons(1, 1)
            blocks2 => blocks_u%sons(1, 2)
            blocks3 => blocks_m%sons(1, 2)
            call Hmat_add_multiply(blocks3, '-', blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks1 => blocks_m%sons(2, 1)
            blocks2 => blocks_u%sons(1, 2)
            blocks3 => blocks_m%sons(2, 2)
            call Hmat_add_multiply(blocks3, '-', blocks1, blocks2, h_mat, option, stats, ptree, msh)

            blocks1 => blocks_u%sons(2, 2)
            blocks2 => blocks_m%sons(1, 2)
            call Hmat_XUM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
            blocks2 => blocks_m%sons(2, 2)
            call Hmat_XUM(blocks1, blocks2, h_mat, option, stats, ptree, msh)
        else if (blocks_m%style == 2) then
            T0 = MPI_Wtime()
            if(blocks_m%level_butterfly==0)then
                rank = size(blocks_m%butterflyV%blocks(1)%matrix,2)
                call Hmat_Usolve(blocks_u, 'T', blocks_u%headm, rank, blocks_m%butterflyV%blocks(1)%matrix, blocks_m%N, ptree, stats)
            else
                T2 = MPI_Wtime()
                rank0 = blocks_m%rankmax
                call BF_randomized(blocks_m%pgno, blocks_m%level_butterfly, rank0, option%rankrate, blocks_m, blocks_u, BF_block_MVP_XUM_dat, error, 'XUM', option, stats, ptree, msh)
                T3 = MPI_Wtime()
                stats%XLUM_random_Time(blocks_m%level) = stats%XLUM_random_Time(blocks_m%level) + T3 - T2
                stats%XLUM_random_CNT(blocks_m%level) = stats%XLUM_random_CNT(blocks_m%level) + 1
            endif
            T1 = MPI_Wtime()
            stats%rankmax_of_level_global_factor(blocks_m%level)=max(stats%rankmax_of_level_global_factor(blocks_m%level),blocks_m%rankmax)
            stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
            stats%Time_XLUM = stats%Time_XLUM + T1 - T0
        else
            T0 = MPI_Wtime()
            if (blocks_m%style == 1) then
#if HAVE_ZFP
                if(option%use_zfp>=1)then
                    call ZFP_Decompress(blocks_m%fullmat,blocks_m%FullmatZFP,blocks_m%M,blocks_m%N,tol_used,0)
                endif
                if(option%use_zfp==1)then
                    call ZFP_Decompress(blocks_u%fullmat,blocks_u%FullmatZFP,blocks_u%M,blocks_u%N,tol_used,1)
                endif

#endif
                mm = size(blocks_m%fullmat, 1)
                nn = size(blocks_m%fullmat, 2)
                call trsmf90(blocks_u%fullmat, blocks_m%fullmat, 'R', 'U', 'N', 'N', mm, nn)
#if HAVE_ZFP
                if(option%use_zfp>=1)then
                    call ZFP_Compress(blocks_m%fullmat,blocks_m%FullmatZFP,blocks_m%M,blocks_m%N,option%tol_comp,0)
                endif
                if(option%use_zfp==1)then
                   call ZFP_Compress(blocks_u%fullmat,blocks_u%FullmatZFP,blocks_u%M,blocks_u%N,option%tol_comp,1)
                endif
#endif
            else

                write (*, *) 'should not come here Hmat_XUM'
                stop
                mm = blocks_m%M
                rank = size(blocks_m%ButterflyV%blocks(1)%matrix, 2)
                call trsmf90(blocks_u%fullmat, blocks_m%ButterflyV%blocks(1)%matrix, 'L', 'U', 'T', 'N', mm, rank)
            endif
            T1 = MPI_Wtime()
            stats%Time_XLUM = stats%Time_XLUM + T1 - T0
        endif

        return

    end subroutine Hmat_XUM

    subroutine Hmat_add_multiply_Hblock3(blocks, chara, block1, block2, h_mat, option, stats, ptree, msh)

        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        integer group_m, group_n, mm, nn,rank,tid
        integer level_blocks, level_butterfly
        character chara, charatmp
        real*8 T0, T1, T2, T3
        type(matrixblock), target :: block1, block2
        type(matrixblock) :: blocks
        real*8:: error_inout
        type(matrixblock), pointer::block_agent
        integer rank0

        tid = 0
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num()
#endif

        if (block1%style == 2) then
            level_butterfly = block1%level_butterfly
            rank0 = block1%rankmax
        endif
        if (block2%style == 2) then
            level_butterfly = block2%level_butterfly
            rank0 = block2%rankmax
        endif

        allocate (block_agent)
        call BF_Init_randomized(level_butterfly, rank0, blocks%row_group, blocks%col_group, blocks, block_agent, msh, ptree, option, 1)

        h_mat%blocks_1(tid+1)%ptr => block1
        h_mat%blocks_2(tid+1)%ptr => block2

        T0 = MPI_Wtime()
        if(level_butterfly==0)then ! use faster, deterministic schemes
            allocate(block_agent%butterflyU%blocks(1))
            allocate(block_agent%butterflyV%blocks(1))
            if (block1%style == 2) then
                rank = size(block1%butterflyU%blocks(1)%matrix,2)
                allocate(block_agent%butterflyU%blocks(1)%matrix(block_agent%M_loc,rank))
                block_agent%butterflyU%blocks(1)%matrix = block1%butterflyU%blocks(1)%matrix
                allocate(block_agent%butterflyV%blocks(1)%matrix(block_agent%N_loc,rank))
                block_agent%butterflyV%blocks(1)%matrix=0
                call Hmat_block_MVP_dat(block2, 'T', block2%headm, block2%headn, rank, block1%butterflyV%blocks(1)%matrix, block1%N_loc, block_agent%butterflyV%blocks(1)%matrix, block_agent%N_loc, BPACK_cone, ptree, stats)
            elseif (block2%style == 2) then
                rank = size(block2%butterflyV%blocks(1)%matrix,2)
                allocate(block_agent%butterflyV%blocks(1)%matrix(block_agent%N_loc,rank))
                block_agent%butterflyV%blocks(1)%matrix = block2%butterflyV%blocks(1)%matrix
                allocate(block_agent%butterflyU%blocks(1)%matrix(block_agent%M_loc,rank))
                block_agent%butterflyU%blocks(1)%matrix=0
                call Hmat_block_MVP_dat(block1, 'N', block1%headm, block1%headn, rank, block2%butterflyU%blocks(1)%matrix, block2%M_loc, block_agent%butterflyU%blocks(1)%matrix, block_agent%M_loc, BPACK_cone, ptree, stats)
            else
                write(*,*)'not supported style of block1 and block2 in Hmat_add_multiply_Hblock3'
                stop
            endif
        else
            T2 = MPI_Wtime()
            call BF_randomized(block_agent%pgno, level_butterfly, rank0, option%rankrate, block_agent, h_mat, BF_block_MVP_Add_Multiply_dat, error_inout, 'Multiply', option, stats, ptree, msh, operand1='m')
            T3 = MPI_Wtime()
            stats%Mul_random_Time(blocks%level) = stats%Mul_random_Time(blocks%level) + T3 - T2
            stats%Mul_random_CNT(blocks%level) = stats%Mul_random_CNT(blocks%level) + 1
        endif
        T1 = MPI_Wtime()

        stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
        stats%Time_Multiply = stats%Time_Multiply + T1 - T0

        call Hmat_BF_add(blocks, chara, block_agent, h_mat, option, stats, ptree, msh)

        call Hmat_block_delete(block_agent)
        deallocate (block_agent)

    end subroutine Hmat_add_multiply_Hblock3


    recursive subroutine Hmat_BF_add(blocks_o, chara, blocks_1, h_mat, option, stats, ptree, msh)
        implicit none

        type(Hoption)::option
        type(Hstat)::stats
        type(Hmat)::h_mat
        type(proctree)::ptree
        type(mesh)::msh

        type(matrixblock)::blocks_o
        type(matrixblock), target::blocks_1
        type(matrixblock), pointer::blocks_1_son, blocks_o_son
        character chara
        real(kind=8) error,flop
        real(kind=8) T0, T1, T2, T3
        integer rank0, tid
        integer rankmax,rankmax2,ranknew,i,j
        DT,allocatable::matUnew(:,:),matVnew(:,:)

        tid = 0
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num()
#endif

        if (blocks_o%style == 4) then
            T0 = MPI_Wtime()
            if (blocks_1%style /= 4) then
                allocate (blocks_1%sons(2, 2))
                call BF_split(blocks_1, blocks_1, ptree, stats, msh, option)
            endif
            T1 = MPI_Wtime()
            stats%Time_Split = stats%Time_Split + T1 - T0

            blocks_o_son => blocks_o%sons(1, 1)
            blocks_1_son => blocks_1%sons(1, 1)
            call Hmat_BF_add(blocks_o_son, chara, blocks_1_son, h_mat, option, stats, ptree, msh)
            blocks_o_son => blocks_o%sons(1, 2)
            blocks_1_son => blocks_1%sons(1, 2)
            call Hmat_BF_add(blocks_o_son, chara, blocks_1_son, h_mat, option, stats, ptree, msh)
            blocks_o_son => blocks_o%sons(2, 1)
            blocks_1_son => blocks_1%sons(2, 1)
            call Hmat_BF_add(blocks_o_son, chara, blocks_1_son, h_mat, option, stats, ptree, msh)
            blocks_o_son => blocks_o%sons(2, 2)
            blocks_1_son => blocks_1%sons(2, 2)
            call Hmat_BF_add(blocks_o_son, chara, blocks_1_son, h_mat, option, stats, ptree, msh)

            if (blocks_1%style /= 4) then
                call BF_delete(blocks_1%sons(1, 1), 1)
                call BF_delete(blocks_1%sons(1, 2), 1)
                call BF_delete(blocks_1%sons(2, 1), 1)
                call BF_delete(blocks_1%sons(2, 2), 1)
                deallocate (blocks_1%sons)
            endif
        else if (blocks_o%style == 2) then
            T0 = MPI_Wtime()
            h_mat%blocks_1(tid+1)%ptr => blocks_1
            rank0 = blocks_o%rankmax


            if(blocks_o%level_butterfly==0)then
                ! if(.false.)then
                    rankmax2 = size(blocks_1%butterflyU%blocks(1)%matrix,2)
                    rankmax = blocks_o%rankmax+rankmax2
                    allocate(matUnew(blocks_o%M,rankmax))
                    matUnew=0
                    allocate(matVnew(rankmax,blocks_o%N))
                    matVnew=0
                    call LR_Add(chara,blocks_o%butterflyU%blocks(1)%matrix,blocks_o%butterflyV%blocks(1)%matrix,blocks_1%butterflyU%blocks(1)%matrix,blocks_1%butterflyV%blocks(1)%matrix,blocks_o%rankmax,rankmax2,ranknew,matUnew,matVnew,blocks_o%M,blocks_o%N,option%tol_rand,flops=flop)
                    stats%Flop_Tmp = flop

                    deallocate(blocks_o%butterflyU%blocks(1)%matrix)
                    allocate(blocks_o%butterflyU%blocks(1)%matrix(blocks_o%M,ranknew))
                    !!$omp parallel do default(shared) private(i,j)
                    do j = 1, ranknew
                        do i = 1, blocks_o%M
                            blocks_o%ButterflyU%blocks(1)%matrix(i, j) = matUnew(i, j)
                        enddo
                     enddo
                     !!$omp end parallel do
                     deallocate(blocks_o%butterflyV%blocks(1)%matrix)
                     allocate(blocks_o%butterflyV%blocks(1)%matrix(blocks_o%N,ranknew))
                     !!$omp parallel do default(shared) private(i,j)
                     do j = 1, ranknew
                        do i = 1, blocks_o%N
                            blocks_o%butterflyV%blocks(1)%matrix(i, j) = matVnew(j, i)
                        enddo
                     enddo
                     !!$omp end parallel do
                     deallocate(matUnew)
                     deallocate(matVnew)
                     blocks_o%rankmax = ranknew

                else
                    T2 = MPI_Wtime()
                    if (chara == '+') then
                    call BF_randomized(blocks_o%pgno, blocks_o%level_butterfly, rank0, option%rankrate, blocks_o, h_mat, BF_block_MVP_Add_Multiply_dat, error, 'Add', option, stats, ptree, msh, operand1='a')
                    elseif (chara == '-') then
                    call BF_randomized(blocks_o%pgno, blocks_o%level_butterfly, rank0, option%rankrate, blocks_o, h_mat, BF_block_MVP_Add_Multiply_dat, error, 'Add', option, stats, ptree, msh, operand1='s')
                    endif
                    T3 = MPI_Wtime()
                    stats%Add_random_Time(blocks_o%level) = stats%Add_random_Time(blocks_o%level) + T3 - T2
                    stats%Add_random_CNT(blocks_o%level) = stats%Add_random_CNT(blocks_o%level) + 1
                endif
            T1 = MPI_Wtime()
            ! time_tmp = time_tmp + T1-T0
            stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
            stats%Time_Add_Multiply = stats%Time_Add_Multiply + T1 - T0
        else if (blocks_o%style == 1) then
            call Full_add(blocks_o, chara, blocks_1, ptree, stats, option)
        end if

    end subroutine Hmat_BF_add

end module BPACK_factor
