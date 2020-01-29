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
module Bplus_compress
   use Bplus_randomizedop
   use BPACK_structure
! use element_Z

contains

   subroutine BF_compress_NlogN(blocks, boundary_map, Nboundall, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)

      use BPACK_DEFS
      implicit none

      integer Nboundall, statflag
      integer boundary_map(*)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new, rank_new1
      integer group_m, group_n, mm, nn, index_i, index_i_loc, index_j_loc, index_j, ii, jj, ij
      integer level, length_1, length_2, level_blocks, index_ij
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e, flops, flops1
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, nnn1, ierr
      real(kind=8) Memory, flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, level_half, level_final, idx_r, inc_r, nr, idx_c, inc_c, nc
      integer passflag
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), select_row_pre(:), select_col_pre(:)
      integer::Nrow_pre, Ncol_pre
      integer::mrange_dummy(1), nrange_dummy(1)
      DT:: mat_dummy(1, 1)
      Memory = 0.

      blocks%rankmax = -100000
      blocks%rankmin = 100000

      group_m = blocks%row_group ! Note: row_group and col_group interchanged here
      group_n = blocks%col_group
      level_blocks = blocks%level

      level_butterfly = blocks%level_butterfly

      allocate (rankmax_for_butterfly(0:level_butterfly))
      rankmax_for_butterfly = -100000
      allocate (rankmin_for_butterfly(0:level_butterfly))
      rankmin_for_butterfly = 100000
      allocate (select_row_pre(blocks%M))
      select_row_pre = 0
      allocate (select_col_pre(blocks%N))
      select_col_pre = 0
      num_blocks = 2**level_butterfly

      blocks%level_half = BF_Switchlevel(level_butterfly, option)
      level_half = blocks%level_half

      if (level_butterfly == 0) then

         allocate (blocks%ButterflyU%blocks(num_blocks))
         allocate (blocks%ButterflyV%blocks(num_blocks))

         ! H-BACA
         leafsize = max(blocks%M, blocks%N)/option%LR_BLK_NUM

         call LR_HBACA(blocks, leafsize, rank_new, option, msh, ker, stats, ptree, blocks%pgno, 0)

         rankmax_for_butterfly(0) = max(blocks%rankmax, rankmax_for_butterfly(0))
         rankmin_for_butterfly(0) = min(blocks%rankmin, rankmin_for_butterfly(0))

      else
         if (level_butterfly /= 0) then
            allocate (blocks%ButterflyKerl(level_butterfly))
            allocate (blocks%ButterflySkel(0:level_butterfly + 1))
         endif

         do level = 0, level_half
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
            if (level == 0) then
               blocks%ButterflyV%idx = idx_c
               blocks%ButterflyV%inc = inc_c
               blocks%ButterflyV%nblk_loc = nc
               allocate (blocks%ButterflyV%blocks(blocks%ButterflyV%nblk_loc))
            elseif (level == level_butterfly + 1) then
               blocks%ButterflyU%idx = idx_r
               blocks%ButterflyU%inc = inc_r
               blocks%ButterflyU%nblk_loc = nr
               allocate (blocks%ButterflyU%blocks(blocks%ButterflyU%nblk_loc))
            else
               blocks%ButterflyKerl(level)%num_row = 2**level
               blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
               blocks%ButterflyKerl(level)%idx_r = idx_r
               blocks%ButterflyKerl(level)%inc_r = inc_r
               blocks%ButterflyKerl(level)%nr = nr
               blocks%ButterflyKerl(level)%idx_c = idx_c*2 - 1
               blocks%ButterflyKerl(level)%inc_c = inc_c
               blocks%ButterflyKerl(level)%nc = nc*2
               allocate (blocks%ButterflyKerl(level)%blocks(blocks%ButterflyKerl(level)%nr, blocks%ButterflyKerl(level)%nc))
            endif

            if (level < level_butterfly + 1) then
               blocks%ButterflySkel(level)%idx_r = idx_r
               blocks%ButterflySkel(level)%inc_r = inc_r
               blocks%ButterflySkel(level)%nr = nr
               blocks%ButterflySkel(level)%idx_c = idx_c
               blocks%ButterflySkel(level)%inc_c = inc_c
               blocks%ButterflySkel(level)%nc = nc
               if (level_half /= level) then ! the last level doesn't require doubling block columns
               if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                  blocks%ButterflySkel(level)%nc = 2
                  blocks%ButterflySkel(level)%idx_c = blocks%ButterflySkel(level)%idx_c - 1 + mod(blocks%ButterflySkel(level)%idx_c, 2)
               endif
               endif
               allocate (blocks%ButterflySkel(level)%inds(blocks%ButterflySkel(level)%nr, blocks%ButterflySkel(level)%nc))
            endif
            rank_new = 0
            flops = 0
            Nrow_pre = 0
            ! !$omp parallel do default(shared) private(index_ij,index_i,index_j,index_i_loc,index_j_loc,rank_new1,flops1) reduction(MAX:rank_new,flops)
            do index_ij = 1, nr*nc
               index_i_loc = (index_ij - 1)/nc + 1
               index_j_loc = mod(index_ij - 1, nc) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c

               ! if(option%sample_heuristic==1)then
               ! if(index_j_loc==1)then

               ! group_m=blocks%row_group    ! Note: row_group and col_group interchanged here
               ! if(level==level_butterfly+1)then
               ! group_m=group_m*2**level_butterfly-1+index_i
               ! else
               ! group_m=group_m*2**level-1+index_i
               ! endif
               ! header_m=msh%basis_group(group_m)%head

               ! ! Nrow_pre=msh%basis_group(group_m)%tail-msh%basis_group(group_m)%head+1
               ! Nrow_pre=0
               ! do ii=1,Nrow_pre
               ! select_row_pre(ii)=header_m+ii-1
               ! enddo
               ! endif
               ! endif

               call BF_compress_NlogN_oneblock_R(blocks, boundary_map, Nboundall, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, level, rank_new1, Nrow_pre, select_row_pre, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops = MAX(flops, flops1)
            enddo
            ! !$omp end parallel do

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
            enddo
            if (level /= level_butterfly + 1) then
            if (level_half == level) then
               call BF_all2all_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'R', 'C')
            else
               call BF_exchange_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'R', 'B')
            endif
            endif

            call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (level < level_butterfly + 1) then
            if (rank_new > rankmax_for_butterfly(level)) then
               rankmax_for_butterfly(level) = rank_new
               stats%Flop_Fill = stats%Flop_Fill + flops
            endif
            endif
         enddo

         level_final = level_half + 1
         do level = level_butterfly + 1, level_final, -1
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
            if (level == 0) then
               blocks%ButterflyV%idx = idx_c
               blocks%ButterflyV%inc = inc_c
               blocks%ButterflyV%nblk_loc = nc
               allocate (blocks%ButterflyV%blocks(blocks%ButterflyV%nblk_loc))
            elseif (level == level_butterfly + 1) then
               blocks%ButterflyU%idx = idx_r
               blocks%ButterflyU%inc = inc_r
               blocks%ButterflyU%nblk_loc = nr
               allocate (blocks%ButterflyU%blocks(blocks%ButterflyU%nblk_loc))
            else
               blocks%ButterflyKerl(level)%num_row = 2**level
               blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
               blocks%ButterflyKerl(level)%idx_r = idx_r*2 - 1
               blocks%ButterflyKerl(level)%inc_r = inc_r
               blocks%ButterflyKerl(level)%nr = nr*2
               blocks%ButterflyKerl(level)%idx_c = idx_c
               blocks%ButterflyKerl(level)%inc_c = inc_c
               blocks%ButterflyKerl(level)%nc = nc
               allocate (blocks%ButterflyKerl(level)%blocks(blocks%ButterflyKerl(level)%nr, blocks%ButterflyKerl(level)%nc))
            endif

            if (level > level_final) then
               blocks%ButterflySkel(level)%idx_r = idx_r
               blocks%ButterflySkel(level)%inc_r = inc_r
               blocks%ButterflySkel(level)%nr = nr
               blocks%ButterflySkel(level)%idx_c = idx_c
               blocks%ButterflySkel(level)%inc_c = inc_c
               blocks%ButterflySkel(level)%nc = nc

               if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                  blocks%ButterflySkel(level)%nr = 2
                  blocks%ButterflySkel(level)%idx_r = blocks%ButterflySkel(level)%idx_r - 1 + mod(blocks%ButterflySkel(level)%idx_r, 2)
               endif
               allocate (blocks%ButterflySkel(level)%inds(blocks%ButterflySkel(level)%nr, blocks%ButterflySkel(level)%nc))
            endif

            rank_new = 0
            flops = 0
            Ncol_pre = 0
            ! !$omp parallel do default(shared) private(index_ij,index_i,index_j,index_j_loc,index_i_loc,rank_new1,flops1) reduction(MAX:rank_new,flops)
            do index_ij = 1, nr*nc
               index_j_loc = (index_ij - 1)/nr + 1
               index_i_loc = mod(index_ij - 1, nr) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c

               ! if(option%sample_heuristic==1)then
               ! if(index_i_loc==1)then

               ! group_n=blocks%col_group
               ! if(level==0)then
               ! group_n=group_n*2**level_butterfly-1+index_j
               ! else
               ! group_n=group_n*2**(level_butterfly-level+1)-1+index_j
               ! endif
               ! header_n=msh%basis_group(group_n)%head

               ! ! Ncol_pre=msh%basis_group(group_n)%tail-msh%basis_group(group_n)%head+1
               ! Ncol_pre=0
               ! do ii=1,Ncol_pre
               ! select_col_pre(ii)=header_n+ii-1
               ! enddo
               ! endif
               ! endif

               call BF_compress_NlogN_oneblock_C(blocks, boundary_map, Nboundall, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, level, level_final, rank_new1, Ncol_pre, select_col_pre, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops = MAX(flops, flops1)
            enddo
            ! !$omp end parallel do

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
            enddo
            if (level > level_final) then
               call BF_exchange_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'C', 'B')
            endif

            call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (level > 0) then
            if (rank_new > rankmax_for_butterfly(level - 1)) then
               rankmax_for_butterfly(level - 1) = rank_new
               stats%Flop_Fill = stats%Flop_Fill + flops
            endif
            endif
         enddo

      endif

      if (statflag == 1) then
         if (allocated(stats%rankmax_of_level)) stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly), stats%rankmax_of_level(level_blocks))
      endif

      deallocate (rankmax_for_butterfly)
      deallocate (rankmin_for_butterfly)
      deallocate (select_row_pre)
      deallocate (select_col_pre)

      if (level_butterfly /= 0) then
         do level = 0, level_butterfly + 1
         if (allocated(blocks%ButterflySkel(level)%inds)) then
            nr = size(blocks%ButterflySkel(level)%inds, 1)
            nc = size(blocks%ButterflySkel(level)%inds, 2)
            do i = 1, nr
            do j = 1, nc
               if (allocated(blocks%ButterflySkel(level)%inds(i, j)%array)) deallocate (blocks%ButterflySkel(level)%inds(i, j)%array)
            enddo
            enddo
            deallocate (blocks%ButterflySkel(level)%inds)
         endif
         enddo
         deallocate (blocks%ButterflySkel)
      endif

      call BF_ComputeMemory(blocks, Memory)
      call BF_get_rank(blocks, ptree)
      return

   end subroutine BF_compress_NlogN

   subroutine BF_compress_NlogN_oneblock_R(blocks, boundary_map, Nboundall, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, level, rank_new, Nrow_pre, select_row_pre, flops)

      use BPACK_DEFS
      implicit none

      integer Nboundall
      integer boundary_map(*)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, jjj, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_r1, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, group_m_mid, group_n_mid, idxstart, idxend, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1, last
      real(kind=8) flop, flops
      DT ctemp
      type(matrixblock)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_row_tmp(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:), order(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, passflag, levelm, nrow, ncol
      integer Nrow_pre
      integer select_row_pre(:)

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:), mnmap(:, :)
      real(kind=8)::n2, n1

      flops = 0
      level_butterfly = blocks%level_butterfly
      group_m = blocks%row_group    ! Note: row_group and col_group interchanged here
      group_n = blocks%col_group
      if (level == level_butterfly + 1) then
         group_m = group_m*2**level_butterfly - 1 + index_i
         group_n = group_n - 1 + index_j
      else
         group_m = group_m*2**level - 1 + index_i
         group_n = group_n*2**(level_butterfly - level) - 1 + index_j
      endif

      mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
      if (level == 0) then
         nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
      elseif (level == level_butterfly + 1) then
         index_ii_loc = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
         index_jj_loc = (index_j - blocks%ButterflySkel(level - 1)%idx_c)/blocks%ButterflySkel(level - 1)%inc_c + 1
         nn = size(blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array, 1)
      else
         index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1
         index_ii_loc = (index_ii - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
         index_jj_loc = (index_jj - blocks%ButterflySkel(level - 1)%idx_c)/blocks%ButterflySkel(level - 1)%inc_c + 1
         nn1 = size(blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array, 1)
         nn2 = size(blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc + 1)%array, 1)
         nn = nn1 + nn2
      endif

      levelm = floor_safe(dble(level_butterfly)/2d0)

      ! select skeletons here, selection of at most (option%sample_para+option%knn)*nn columns, the first option%sample_para*nn are random, the next option%knn*nn are nearest points
      rankmax_r1 = min(ceiling_safe(option%sample_para*nn), mm)
      if (level == 0) rankmax_r1 = min(ceiling_safe(option%sample_para*nn), mm)
      rankmax_c = nn
      allocate (select_row(rankmax_r1 + nn*option%knn + Nrow_pre))
      allocate (select_column(rankmax_c))
      do i = 1, rankmax_c
         select_column(i) = i
      enddo
      call linspaceI(1, mm, rankmax_r1, select_row(1:rankmax_r1))
      header_m = msh%basis_group(group_m)%head
      header_n = msh%basis_group(group_n)%head
      if (2*group_n + 1 <= size(msh%basis_group, 1)) header_n2 = msh%basis_group(2*group_n + 1)%head
      if (level /= level_butterfly + 1) then
      do j = 1, nn
         if (level == 0) then
            edge_n = header_n + j - 1
         elseif (level < level_butterfly + 1) then
            if (j <= nn1) then
               edge_n = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array(j) + header_n - 1
            else
               edge_n = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc + 1)%array(j - nn1) + header_n2 - 1
            endif
         endif
         do jjj = 1, option%knn
         if (msh%nns(edge_n, jjj) >= msh%basis_group(group_m)%head .and. msh%nns(edge_n, jjj) <= msh%basis_group(group_m)%tail) then
            rankmax_r1 = rankmax_r1 + 1
            select_row(rankmax_r1) = msh%nns(edge_n, jjj) + 1 - header_m
         endif
         enddo
      enddo
      endif

      do ii = 1, Nrow_pre
         if (select_row_pre(ii) >= msh%basis_group(group_m)%head .and. select_row_pre(ii) <= msh%basis_group(group_m)%tail) then
         rankmax_r1 = rankmax_r1 + 1
         select_row(rankmax_r1) = select_row_pre(ii) + 1 - header_m
         endif
      enddo

      call remove_dup_int(select_row, rankmax_r1, rankmax_r)

      if (level == 0) then

         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head

         ! allocate (core(rankmax_r,rankmax_c))
         ! allocate(mrange(rankmax_r))
         ! allocate(nrange(rankmax_c))
         ! do i=1,rankmax_r
         ! edge_m=header_m+select_row(i)-1
         ! enddo
         ! do j=1,rankmax_c
         ! edge_n=header_n+select_column(j)-1
         ! enddo
         ! call element_Zmn_block_user(rankmax_r,rankmax_c,mrange,nrange,core,msh,option,ker,0,passflag,ptree,stats)
         ! deallocate(mrange)
         ! deallocate(nrange)

         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,k,edge_m,edge_n,ctemp)
         ! do ij=1,rankmax_c*rankmax_r
         ! j = (ij-1)/rankmax_r+1
         ! i = mod(ij-1,rankmax_r) + 1
         ! edge_m=header_m+select_row(i)-1
         ! edge_n=header_n+select_column(j)-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! core(i,j)=ctemp
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel

         allocate (matrix_V(rankmax_r, nn))
         allocate (mrange(rankmax_r))
         allocate (nrange(nn))
         do i = 1, rankmax_r
            mrange(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, nn
            nrange(j) = header_n + j - 1
         enddo

         if (Nboundall > 0) then

            allocate (mnmap(rankmax_r, nn))
            allocate (mmap(rankmax_r))
            allocate (nmap(nn))
            allocate (mrange1(rankmax_r))
            allocate (nrange1(nn))
            mnmap = 1
            do i = 1, rankmax_r
               group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
               group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
               if (group_n_mid /= -1) then
                  idxstart = msh%basis_group(group_n_mid)%head
                  idxend = msh%basis_group(group_n_mid)%tail
                  do j = 1, nn
                     if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                  enddo
               endif
            enddo

            nrow = 0
            do i = 1, rankmax_r
               if (sum(mnmap(i, :)) /= 0) then
                  nrow = nrow + 1
                  mmap(nrow) = i
                  mrange1(nrow) = mrange(i)
               endif
            enddo
            ncol = 0
            do j = 1, nn
               if (sum(mnmap(:, j)) /= 0) then
                  ncol = ncol + 1
                  nmap(ncol) = j
                  nrange1(ncol) = nrange(j)
               endif
            enddo
            allocate (matrix_tmp(nrow, ncol))
            call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

            matrix_V = 0
            do i = 1, nrow
            do j = 1, ncol
               matrix_V(mmap(i), nmap(j)) = matrix_tmp(i, j)
            enddo
            enddo
            deallocate (mnmap)
            deallocate (mmap)
            deallocate (nmap)
            deallocate (mrange1)
            deallocate (nrange1)
            deallocate (matrix_tmp)
         else
            call element_Zmn_block_user(rankmax_r, nn, mrange, nrange, matrix_V, msh, option, ker, 0, passflag, ptree, stats)
         endif

         deallocate (mrange)
         deallocate (nrange)

         allocate (core(rankmax_r, rankmax_c))
         allocate (core_tmp(rankmax_c, rankmax_r))
         do j = 1, rankmax_c
            core(:, j) = matrix_V(:, select_column(j))
         enddo
         call copymatT(core, core_tmp, rankmax_r, rankmax_c)

         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
         ! do ij=1,nn*rankmax_r
         ! j = (ij-1)/rankmax_r+1
         ! i = mod(ij-1,rankmax_r) + 1

         ! edge_m=header_m+select_row(i)-1
         ! edge_n=header_n+j-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! matrix_V(i,j)=ctemp
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel

         allocate (jpvt(max(rankmax_c, rankmax_r)))
         allocate (tau(max(rankmax_c, rankmax_r)))
         jpvt = 0
         call geqp3modf90(core, jpvt, tau, option%tol_comp, SafeUnderflow, rank_new, flop=flop)
         flops = flops + flop

         if (rank_new > 0) then
            call un_or_mqrf90(core, tau, matrix_V, 'L', 'C', rankmax_r, nn, rank_new, flop=flop)
            flops = flops + flop
            call trsmf90(core, matrix_V, 'L', 'U', 'N', 'N', rank_new, nn, flop=flop)
            flops = flops + flop
         else
            rank_new = 1
            matrix_V = 0
         endif

         index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1
         allocate (blocks%ButterflyV%blocks(index_j_loc_k)%matrix(nn, rank_new))
         call copymatT(matrix_V(1:rank_new, 1:nn), blocks%ButterflyV%blocks(index_j_loc_k)%matrix, rank_new, nn)

         index_j_loc_s = (index_j - blocks%ButterflySkel(0)%idx_c)/blocks%ButterflySkel(0)%inc_c + 1
         allocate (blocks%ButterflySkel(0)%inds(1, index_j_loc_s)%array(rank_new))
         do j = 1, rank_new
            blocks%ButterflySkel(0)%inds(1, index_j_loc_s)%array(j) = select_column(jpvt(j))
         enddo

         if (option%sample_heuristic == 1) then
            allocate (core_tmp1(rank_new, rankmax_r))
            do i = 1, rank_new
               core_tmp1(i, :) = core_tmp(jpvt(i), :)
            enddo
            jpvt = 0
            call geqp3modf90(core_tmp1, jpvt, tau, option%tol_comp, SafeUnderflow, Nrow_pre, flop=flop)
            flops = flops + flop
            Nrow_pre = rank_new
            do i = 1, Nrow_pre
               select_row_pre(i) = header_m + select_row(jpvt(i)) - 1
            enddo
            deallocate (core_tmp1)
         endif

      elseif (level == level_butterfly + 1) then
         index_i_loc_s = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
         rank_new = size(blocks%ButterflySkel(level - 1)%inds(index_i_loc_s, 1)%array)
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head
         index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1

         allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
         allocate (mrange(mm))
         allocate (nrange(rank_new))
         do i = 1, mm
            mrange(i) = i + header_m - 1
         enddo
         do j = 1, rank_new
            nrange(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s, 1)%array(j) + header_n - 1
         enddo

         if (Nboundall > 0) then

            allocate (mnmap(mm, rank_new))
            allocate (mmap(mm))
            allocate (nmap(rank_new))
            allocate (mrange1(mm))
            allocate (nrange1(rank_new))
            mnmap = 1

            do i = 1, mm
               group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
               group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
               if (group_n_mid /= -1) then
                  idxstart = msh%basis_group(group_n_mid)%head
                  idxend = msh%basis_group(group_n_mid)%tail
                  do j = 1, rank_new
                     if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                  enddo
               endif
            enddo

            nrow = 0
            do i = 1, mm
               if (sum(mnmap(i, :)) /= 0) then
                  nrow = nrow + 1
                  mmap(nrow) = i
                  mrange1(nrow) = mrange(i)
               endif
            enddo
            ncol = 0
            do j = 1, rank_new
               if (sum(mnmap(:, j)) /= 0) then
                  ncol = ncol + 1
                  nmap(ncol) = j
                  nrange1(ncol) = nrange(j)
               endif
            enddo
            allocate (matrix_tmp(nrow, ncol))
            call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

            blocks%ButterflyU%blocks(index_i_loc_k)%matrix = 0
            do i = 1, nrow
            do j = 1, ncol
               blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mmap(i), nmap(j)) = matrix_tmp(i, j)
            enddo
            enddo
            deallocate (mnmap)
            deallocate (mmap)
            deallocate (nmap)
            deallocate (mrange1)
            deallocate (nrange1)
            deallocate (matrix_tmp)
         else
            call element_Zmn_block_user(mm, rank_new, mrange, nrange, blocks%ButterflyU%blocks(index_i_loc_k)%matrix, msh, option, ker, 0, passflag, ptree, stats)
         endif

         deallocate (mrange)
         deallocate (nrange)

         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
         ! do ij=1,rank_new*mm
         ! j = (ij-1)/mm+1
         ! i = mod(ij-1,mm) + 1
         ! edge_m=i+header_m-1
         ! edge_n=blocks%ButterflySkel(level-1)%inds(index_i_loc_s,1)%array(j)+header_n-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! blocks%ButterflyU%blocks(index_i_loc_k)%matrix(i,j)=ctemp
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel
      else
         index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
         index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
         index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1
         index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1
         header_m = msh%basis_group(group_m)%head
         header_n1 = msh%basis_group(group_n)%head
         header_n2 = msh%basis_group(2*group_n + 1)%head
         nnn1 = msh%basis_group(2*group_n)%tail - msh%basis_group(2*group_n)%head + 1

         ! allocate (core(rankmax_r,rankmax_c))
         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,k,edge_m,edge_n,ctemp)
         ! do ij=1,rankmax_c*rankmax_r
         ! j = (ij-1)/rankmax_r+1
         ! i = mod(ij-1,rankmax_r) + 1
         ! if (select_column(j)<=nn1) then
         ! edge_m=select_row(i)+header_m-1
         ! edge_n=blocks%ButterflySkel(level-1)%inds(index_ii_loc,index_jj_loc)%array(select_column(j))+header_n1-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! core(i,j)=ctemp
         ! else
         ! edge_m=select_row(i)+header_m-1
         ! edge_n=blocks%ButterflySkel(level-1)%inds(index_ii_loc,index_jj_loc+1)%array(select_column(j)-nn1)+header_n2-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! core(i,j)=ctemp
         ! endif
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel

         allocate (matrix_V(rankmax_r, nn))
         allocate (mrange(rankmax_r))
         allocate (nrange(nn))
         do i = 1, rankmax_r
            mrange(i) = select_row(i) + header_m - 1
         enddo
         do j = 1, nn
            if (j <= nn1) then
               nrange(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array(j) + header_n1 - 1
            else
               nrange(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc + 1)%array(j - nn1) + header_n2 - 1
            endif
         enddo

         if (Nboundall > 0) then

            allocate (mnmap(rankmax_r, nn))
            allocate (mmap(rankmax_r))
            allocate (nmap(nn))
            allocate (mrange1(rankmax_r))
            allocate (nrange1(nn))
            mnmap = 1

            do i = 1, rankmax_r
               group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
               group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
               if (group_n_mid /= -1) then
                  idxstart = msh%basis_group(group_n_mid)%head
                  idxend = msh%basis_group(group_n_mid)%tail
                  do j = 1, nn
                     if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                  enddo
               endif
            enddo

            nrow = 0
            do i = 1, rankmax_r
               if (sum(mnmap(i, :)) /= 0) then
                  nrow = nrow + 1
                  mmap(nrow) = i
                  mrange1(nrow) = mrange(i)
               endif
            enddo
            ncol = 0
            do j = 1, nn
               if (sum(mnmap(:, j)) /= 0) then
                  ncol = ncol + 1
                  nmap(ncol) = j
                  nrange1(ncol) = nrange(j)
               endif
            enddo
            allocate (matrix_tmp(nrow, ncol))
            call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

            matrix_V = 0
            do i = 1, nrow
            do j = 1, ncol
               matrix_V(mmap(i), nmap(j)) = matrix_tmp(i, j)
            enddo
            enddo
            deallocate (mnmap)
            deallocate (mmap)
            deallocate (nmap)
            deallocate (mrange1)
            deallocate (nrange1)
            deallocate (matrix_tmp)
         else
            call element_Zmn_block_user(rankmax_r, nn, mrange, nrange, matrix_V, msh, option, ker, 0, passflag, ptree, stats)
         endif

         deallocate (mrange)
         deallocate (nrange)

         allocate (core(rankmax_r, rankmax_c))
         allocate (core_tmp(rankmax_c, rankmax_r))
         do j = 1, rankmax_c
            core(:, j) = matrix_V(:, select_column(j))
         enddo
         call copymatT(core, core_tmp, rankmax_r, rankmax_c)

         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
         ! do ij=1,nn*rankmax_r
         ! j = (ij-1)/rankmax_r+1
         ! i = mod(ij-1,rankmax_r) + 1
         ! if (j<=nn1) then
         ! edge_m=select_row(i)+header_m-1
         ! edge_n=blocks%ButterflySkel(level-1)%inds(index_ii_loc,index_jj_loc)%array(j)+header_n1-1
         ! else
         ! edge_m=select_row(i)+header_m-1
         ! edge_n=blocks%ButterflySkel(level-1)%inds(index_ii_loc,index_jj_loc+1)%array(j-nn1)+header_n2-1
         ! endif
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! matrix_V(i,j)=ctemp
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel

         allocate (jpvt(max(rankmax_c, rankmax_r)))
         allocate (tau(max(rankmax_c, rankmax_r)))
         jpvt = 0
         call geqp3modf90(core, jpvt, tau, option%tol_comp, SafeUnderflow, rank_new, flop=flop)
         flops = flops + flop

         if (rank_new > 0) then
            call un_or_mqrf90(core, tau, matrix_V, 'L', 'C', rankmax_r, nn, rank_new, flop=flop)
            flops = flops + flop
            call trsmf90(core, matrix_V, 'L', 'U', 'N', 'N', rank_new, nn, flop=flop)
            flops = flops + flop
         else
            rank_new = 1
            matrix_V = 0
         endif

         allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(rank_new, nn1))
         allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix(rank_new, nn2))

         blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrix_V(1:rank_new, 1:nn1)
         blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix = matrix_V(1:rank_new, 1 + nn1:nn)

         allocate (blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(rank_new))
         ! !$omp taskloop default(shared) private(j)
         do j = 1, rank_new
            if (select_column(jpvt(j)) <= nn1) then
               blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array(select_column(jpvt(j)))
            else
               blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc + 1)%array(select_column(jpvt(j)) - nn1) + nnn1
            endif
         enddo
         ! !$omp end taskloop

         if (option%sample_heuristic == 1) then
            allocate (core_tmp1(rank_new, rankmax_r))
            do i = 1, rank_new
               core_tmp1(i, :) = core_tmp(jpvt(i), :)
            enddo
            jpvt = 0
            call geqp3modf90(core_tmp1, jpvt, tau, option%tol_comp, SafeUnderflow, Nrow_pre, flop=flop)
            flops = flops + flop
            Nrow_pre = rank_new
            do i = 1, Nrow_pre
               select_row_pre(i) = header_m + select_row(jpvt(i)) - 1
            enddo
            deallocate (core_tmp1)
         endif

      endif

      if (allocated(core)) deallocate (core)
      if (allocated(core_tmp)) deallocate (core_tmp)
      if (allocated(tau)) deallocate (tau)
      if (allocated(jpvt)) deallocate (jpvt)
      if (allocated(matrix_V)) deallocate (matrix_V)
      if (allocated(select_row)) deallocate (select_row)
      if (allocated(select_column)) deallocate (select_column)

   end subroutine BF_compress_NlogN_oneblock_R

   subroutine BF_exchange_skel(blocks, skels, option, stats, msh, ptree, level, mode, collect)

      use BPACK_DEFS
      implicit none
      type(mesh)::msh
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(butterfly_skel)::skels

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, pgno_sub, pid0, tag, nproc, Ncol, Nskel, Nreqr, Nreqs, recvid, sendid

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive
      logical::sendflag
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)

      character mode, modetrans, collect

      real(kind=8)::n1, n2

      n1 = OMP_get_wtime()

      if (mode == 'R') modetrans = 'C'
      if (mode == 'C') modetrans = 'R'

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno

      ! allocation of communication quantities
      allocate (statuss(MPI_status_size, nproc))
      allocate (statusr(MPI_status_size, nproc))
      allocate (S_req(nproc))
      allocate (R_req(nproc))
      allocate (sendquant(nproc))
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      allocate (recvquant(nproc))
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo
      allocate (sendIDactive(nproc))
      allocate (recvIDactive(nproc))
      Nsendactive = 0
      Nrecvactive = 0

      ! calculate send buffer sizes in the first pass
      do ii = 1, skels%nr
      do jj = 1, skels%nc
         index_i = (ii - 1)*skels%inc_r + skels%idx_r
         index_j = (jj - 1)*skels%inc_c + skels%idx_c

         sendflag = .false.
         if (collect == 'R') then ! pair-wise reduction
            if (mode == 'R') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = 2*index_j - mod(index_i, 2)
            elseif (mode == 'C') then
               index_i0 = 2*index_i - mod(index_j, 2)
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, modetrans, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         elseif (collect == 'B') then ! pair-wise broadcast
            if (mode == 'R') then
               index_j0 = index_j + 2*mod(index_j, 2) - 1
               index_i0 = index_i
            elseif (mode == 'C') then
               index_i0 = index_i + 2*mod(index_i, 2) - 1
               index_j0 = index_j
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         endif
         sendflag = pid /= ptree%MyID

         if (sendflag) then
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif

            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 3 + size(skels%inds(ii, jj)%array)
         endif
      enddo
      enddo

      ! communicate receive buffer sizes
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat_i(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(tt), ierr)
      enddo

      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(tt), ierr)
      enddo
      if (Nsendactive > 0) then
         call MPI_waitall(Nsendactive, S_req, statuss, ierr)
      endif
      if (Nrecvactive > 0) then
         call MPI_waitall(Nrecvactive, R_req, statusr, ierr)
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat_i(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do ii = 1, skels%nr
      do jj = 1, skels%nc

         index_i = (ii - 1)*skels%inc_r + skels%idx_r
         index_j = (jj - 1)*skels%inc_c + skels%idx_c

         sendflag = .false.
         if (collect == 'R') then ! pair-wise reduction
            if (mode == 'R') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = 2*index_j - mod(index_i, 2)
            elseif (mode == 'C') then
               index_i0 = 2*index_i - mod(index_j, 2)
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, modetrans, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         elseif (collect == 'B') then ! pair-wise broadcast
            if (mode == 'R') then
               index_j0 = index_j + 2*mod(index_j, 2) - 1
               index_i0 = index_i
            elseif (mode == 'C') then
               index_i0 = index_i + 2*mod(index_i, 2) - 1
               index_j0 = index_j
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         endif
         sendflag = pid /= ptree%MyID

         if (sendflag) then
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            Nskel = size(skels%inds(ii, jj)%array)
            sendquant(pp)%dat_i(sendquant(pp)%size + 1, 1) = index_i
            sendquant(pp)%dat_i(sendquant(pp)%size + 2, 1) = index_j
            sendquant(pp)%dat_i(sendquant(pp)%size + 3, 1) = Nskel
            sendquant(pp)%size = sendquant(pp)%size + 3
            do i = 1, Nskel
               sendquant(pp)%dat_i(sendquant(pp)%size + i, 1) = skels%inds(ii, jj)%array(i)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + Nskel
         endif
      enddo
      enddo

      ! communicate the data buffer
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(tt), ierr)
      enddo

      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(tt), ierr)
      enddo

      if (Nsendactive > 0) then
         call MPI_waitall(Nsendactive, S_req, statuss, ierr)
      endif
      if (Nrecvactive > 0) then
         call MPI_waitall(Nrecvactive, R_req, statusr, ierr)
      endif

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            ii = (index_i - skels%idx_r)/skels%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            jj = (index_j - skels%idx_c)/skels%inc_c + 1
            i = i + 1
            Nskel = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            if (.not. allocated(skels%inds(ii, jj)%array)) then
               allocate (skels%inds(ii, jj)%array(Nskel))
               skels%inds(ii, jj)%array = 0
            endif
            skels%inds(ii, jj)%array = skels%inds(ii, jj)%array + NINT(dble(recvquant(pp)%dat_i(i + 1:i + Nskel, 1)))
            i = i + Nskel
         enddo
      enddo

      ! deallocation
      deallocate (S_req)
      deallocate (R_req)
      deallocate (statuss)
      deallocate (statusr)
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat_i)) deallocate (sendquant(pp)%dat_i)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat_i)) deallocate (recvquant(pp)%dat_i)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)
      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_exchange_skel

!*********** all to all communication of skeletons of one butterfly level from row-wise ordering to column-wise ordering or the reverse
   subroutine BF_all2all_skel(blocks, skels, option, stats, msh, ptree, level, mode, mode_new)

      use BPACK_DEFS
      implicit none
      type(mesh)::msh
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, pgno_sub, tag, nproc, Ncol, Nskel, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      character::mode, mode_new
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      logical all2all
      type(butterfly_skel)::skels
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      integer, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist

      n1 = OMP_get_wtime()

      call assert(mode /= mode_new, 'only row2col or col2row is supported')

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno

      ! mode_new and level_new determine the block range in the new mode
      if (mode_new == 'R') then
         level_new = max(level - 1, 0)
      elseif (mode_new == 'C') then
         level_new = min(level + 1, level_butterfly + 1)
      endif
      call GetLocalBlockRange(ptree, blocks%pgno, level_new, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, mode_new)

      ! allocation of communication quantities
      allocate (statuss(MPI_status_size, nproc))
      allocate (statusr(MPI_status_size, nproc))
      allocate (S_req(nproc))
      allocate (R_req(nproc))
      allocate (sendquant(nproc))
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      allocate (recvquant(nproc))
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo
      allocate (sendIDactive(nproc))
      allocate (recvIDactive(nproc))
      Nsendactive = 0
      Nrecvactive = 0

      ! calculate send buffer sizes in the first pass
      do ii = 1, nr
      do jj = 1, nc
         index_i = (ii - 1)*inc_r + idx_r
         index_j = (jj - 1)*inc_c + idx_c
         call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i, index_j, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         if (recvquant(pp)%active == 0) then
            recvquant(pp)%active = 1
            Nrecvactive = Nrecvactive + 1
            recvIDactive(Nrecvactive) = pp
         endif
      enddo
      enddo

      do ii = 1, skels%nr
      do jj = 1, skels%nc
         index_i = (ii - 1)*skels%inc_r + skels%idx_r
         index_j = (jj - 1)*skels%inc_c + skels%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         if (sendquant(pp)%active == 0) then
            sendquant(pp)%active = 1
            Nsendactive = Nsendactive + 1
            sendIDactive(Nsendactive) = pp
         endif
         sendquant(pp)%size = sendquant(pp)%size + 3 + size(skels%inds(ii, jj)%array)
      enddo
      enddo

      ! communicate receive buffer sizes
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat_i(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(tt), ierr)
      enddo

      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(tt), ierr)
      enddo
      if (Nsendactive > 0) then
         call MPI_waitall(Nsendactive, S_req, statuss, ierr)
      endif
      if (Nrecvactive > 0) then
         call MPI_waitall(Nrecvactive, R_req, statusr, ierr)
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat_i(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do ii = 1, skels%nr
      do jj = 1, skels%nc
         index_i = (ii - 1)*skels%inc_r + skels%idx_r
         index_j = (jj - 1)*skels%inc_c + skels%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head

         pp = pid - ptree%pgrp(blocks%pgno)%head + 1

         Nskel = size(skels%inds(ii, jj)%array)
         sendquant(pp)%dat_i(sendquant(pp)%size + 1, 1) = index_i
         sendquant(pp)%dat_i(sendquant(pp)%size + 2, 1) = index_j
         sendquant(pp)%dat_i(sendquant(pp)%size + 3, 1) = Nskel
         sendquant(pp)%size = sendquant(pp)%size + 3
         do i = 1, Nskel
            sendquant(pp)%dat_i(sendquant(pp)%size + i, 1) = skels%inds(ii, jj)%array(i)
         enddo
         deallocate (skels%inds(ii, jj)%array)
         sendquant(pp)%size = sendquant(pp)%size + Nskel
      enddo
      enddo

      deallocate (skels%inds)
      skels%idx_r = idx_r
      skels%idx_c = idx_c
      skels%inc_r = inc_r
      skels%inc_c = inc_c
      skels%nr = nr
      skels%nc = nc

      allocate (skels%inds(skels%nr, skels%nc))

      ! communicate the data buffer

      call MPI_ALLREDUCE(Nsendactive, Nsendactive_min, 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(blocks%pgno)%Comm, ierr)
      call MPI_ALLREDUCE(Nrecvactive, Nrecvactive_min, 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(blocks%pgno)%Comm, ierr)
#if 0
      all2all = (nproc == Nsendactive_min .and. Nsendactive_min == Nrecvactive_min)
#else
      all2all = .false.
#endif

      ! if(ptree%MyID==Main_ID)write(*,*)'Nactive: ',Nsendactive_min,' Nproc:', nproc

      if (all2all) then ! if truly all-to-all, use MPI_ALLTOALLV
         allocate (sdispls(nproc))
         allocate (sendcounts(nproc))
         dist = 0
         do pp = 1, nproc
            sendcounts(pp) = sendquant(pp)%size
            sdispls(pp) = dist
            dist = dist + sendquant(pp)%size
         enddo
         allocate (sendbufall2all(dist))
         do pp = 1, nproc
            if (sendquant(pp)%size > 0) sendbufall2all(sdispls(pp) + 1:sdispls(pp) + sendcounts(pp)) = sendquant(pp)%dat_i(:, 1)
         enddo

         allocate (rdispls(nproc))
         allocate (recvcounts(nproc))
         dist = 0
         do pp = 1, nproc
            recvcounts(pp) = recvquant(pp)%size
            rdispls(pp) = dist
            dist = dist + recvquant(pp)%size
         enddo
         allocate (recvbufall2all(dist))

         call MPI_ALLTOALLV(sendbufall2all, sendcounts, sdispls, MPI_INTEGER, recvbufall2all, recvcounts, rdispls, MPI_INTEGER, ptree%pgrp(blocks%pgno)%Comm, ierr)

         do pp = 1, nproc
            if (recvquant(pp)%size > 0) recvquant(pp)%dat_i(:, 1) = recvbufall2all(rdispls(pp) + 1:rdispls(pp) + recvcounts(pp))
         enddo

         deallocate (sdispls)
         deallocate (sendcounts)
         deallocate (sendbufall2all)
         deallocate (rdispls)
         deallocate (recvcounts)
         deallocate (recvbufall2all)

      else
         Nreqs = 0
         do tt = 1, Nsendactive
            pp = sendIDactive(tt)
            recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
            if (recvid /= ptree%MyID) then
               Nreqs = Nreqs + 1
               call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
            else
               if (sendquant(pp)%size > 0) recvquant(pp)%dat_i = sendquant(pp)%dat_i
            endif
         enddo

         Nreqr = 0
         do tt = 1, Nrecvactive
            pp = recvIDactive(tt)
            sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
            if (sendid /= ptree%MyID) then
               Nreqr = Nreqr + 1
               call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
            endif
         enddo

         if (Nreqs > 0) then
            call MPI_waitall(Nreqs, S_req, statuss, ierr)
         endif
         if (Nreqr > 0) then
            call MPI_waitall(Nreqr, R_req, statusr, ierr)
         endif
      endif

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            ii = (index_i - skels%idx_r)/skels%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            jj = (index_j - skels%idx_c)/skels%inc_c + 1
            i = i + 1
            Nskel = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            call assert(.not. allocated(skels%inds(ii, jj)%array), 'receiving dat_ia alreay exists locally')
            allocate (skels%inds(ii, jj)%array(Nskel))
            skels%inds(ii, jj)%array = NINT(dble(recvquant(pp)%dat_i(i + 1:i + Nskel, 1)))
            i = i + Nskel
         enddo
      enddo

      ! deallocation
      deallocate (S_req)
      deallocate (R_req)
      deallocate (statuss)
      deallocate (statusr)
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat_i)) deallocate (sendquant(pp)%dat_i)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat_i)) deallocate (recvquant(pp)%dat_i)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_skel

   subroutine BF_compress_NlogN_oneblock_C(blocks, boundary_map, Nboundall, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, level, level_final, rank_new, Ncol_pre, select_col_pre, flops)

      use BPACK_DEFS
      implicit none

      integer Nboundall
      integer boundary_map(*)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, iii, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_c1, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_j, index_i_loc_k, index_j_loc_k, index_i_loc_s, index_i_loc_s1, index_j_loc_s, index_j_loc_s1, ii, jj, ij
      integer level, level_final, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_m1, header_m2, mm1, mm2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, mmm1
      real(kind=8) flop, flops
      DT ctemp
      type(matrixblock)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, levelm, group_m_mid, group_n_mid, idxstart, idxend, nrow, ncol
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:), mnmap(:, :)
      integer::passflag = 0
      real(kind=8)::n2, n1
      integer::Ncol_pre
      integer::select_col_pre(:)

      flops = 0
      level_butterfly = blocks%level_butterfly
      group_m = blocks%row_group    ! Note: row_group and col_group interchanged here
      group_n = blocks%col_group
      if (level == 0) then
         group_m = group_m - 1 + index_i
         group_n = group_n*2**level_butterfly - 1 + index_j
      else
         group_m = group_m*2**(level - 1) - 1 + index_i
         group_n = group_n*2**(level_butterfly - level + 1) - 1 + index_j
      endif

      nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
      if (level == level_butterfly + 1) then
         mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
      elseif (level == 0) then
         index_ii_loc = (index_i - blocks%ButterflySkel(level + 1)%idx_r)/blocks%ButterflySkel(level + 1)%inc_r + 1
         index_jj_loc = (index_j - blocks%ButterflySkel(level + 1)%idx_c)/blocks%ButterflySkel(level + 1)%inc_c + 1
         mm = size(blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array, 1)
      else
         index_jj = int((index_j + 1)/2); index_ii = 2*index_i - 1
         index_ii_loc = (index_ii - blocks%ButterflySkel(level + 1)%idx_r)/blocks%ButterflySkel(level + 1)%inc_r + 1
         index_jj_loc = (index_jj - blocks%ButterflySkel(level + 1)%idx_c)/blocks%ButterflySkel(level + 1)%inc_c + 1
         mm1 = size(blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array, 1)
         mm2 = size(blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array, 1)
         mm = mm1 + mm2
      endif

      levelm = floor_safe(dble(level_butterfly)/2d0)

      ! select skeletons here, selection of at most (option%sample_para+option%knn)*mm rows, the first option%sample_para*mm are random, the next option%knn*mm are nearest points
      rankmax_r = mm
      rankmax_c1 = min(nn, ceiling_safe(option%sample_para*mm))
      if (level == level_butterfly + 1) rankmax_c1 = min(ceiling_safe(option%sample_para*mm), nn)
      allocate (select_row(rankmax_r))
      allocate (select_column(rankmax_c1 + option%knn*mm + Ncol_pre))
      do i = 1, rankmax_r
         select_row(i) = i
      enddo
      call linspaceI(1, nn, rankmax_c1, select_column(1:rankmax_c1))
      header_m = msh%basis_group(group_m)%head
      header_n = msh%basis_group(group_n)%head
      if (2*group_m + 1 <= size(msh%basis_group, 1)) header_m2 = msh%basis_group(2*group_m + 1)%head
      if (level /= 0) then
      do i = 1, mm
         if (level == level_butterfly + 1) then
            edge_m = header_m + i - 1
         elseif (level > 0) then
            if (i <= mm1) then
               edge_m = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(i) + header_m - 1
            else
               edge_m = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(i - mm1) + header_m2 - 1
            endif
         endif
         do iii = 1, option%knn
         if (msh%nns(edge_m, iii) >= msh%basis_group(group_n)%head .and. msh%nns(edge_m, iii) <= msh%basis_group(group_n)%tail) then
            rankmax_c1 = rankmax_c1 + 1
            select_column(rankmax_c1) = msh%nns(edge_m, iii) + 1 - header_n
         endif
         enddo
      enddo
      endif

      do ii = 1, Ncol_pre
         if (select_col_pre(ii) >= msh%basis_group(group_n)%head .and. select_col_pre(ii) <= msh%basis_group(group_n)%tail) then
         rankmax_c1 = rankmax_c1 + 1
         select_column(rankmax_c1) = select_col_pre(ii) + 1 - header_n
         endif
      enddo

      call remove_dup_int(select_column, rankmax_c1, rankmax_c)

      if (level == level_butterfly + 1) then
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head

         if (level == level_final) then
            index_i_loc_s1 = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
            index_j_loc_s1 = (index_j - blocks%ButterflySkel(level - 1)%idx_c)/blocks%ButterflySkel(level - 1)%inc_c + 1

            rank_new = size(blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array)
            index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
            allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))

            allocate (mrange(mm))
            allocate (nrange(rank_new))
            do i = 1, mm
               mrange(i) = header_m + i - 1
            enddo
            do j = 1, rank_new
               nrange(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (mnmap(mm, rank_new))
               allocate (mmap(mm))
               allocate (nmap(rank_new))
               allocate (mrange1(mm))
               allocate (nrange1(rank_new))
               mnmap = 1

               do i = 1, mm
                  group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, rank_new
                        if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                     enddo
                  endif
               enddo

               nrow = 0
               do i = 1, mm
                  if (sum(mnmap(i, :)) /= 0) then
                     nrow = nrow + 1
                     mmap(nrow) = i
                     mrange1(nrow) = mrange(i)
                  endif
               enddo
               ncol = 0
               do j = 1, rank_new
                  if (sum(mnmap(:, j)) /= 0) then
                     ncol = ncol + 1
                     nmap(ncol) = j
                     nrange1(ncol) = nrange(j)
                  endif
               enddo
               allocate (matrix_tmp(nrow, ncol))
               call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

               blocks%ButterflyU%blocks(index_i_loc_k)%matrix = 0
               do i = 1, nrow
               do j = 1, ncol
                  blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mmap(i), nmap(j)) = matrix_tmp(i, j)
               enddo
               enddo
               deallocate (mnmap)
               deallocate (mmap)
               deallocate (nmap)
               deallocate (mrange1)
               deallocate (nrange1)
               deallocate (matrix_tmp)
            else
               call element_Zmn_block_user(mm, rank_new, mrange, nrange, blocks%ButterflyU%blocks(index_i_loc_k)%matrix, msh, option, ker, 0, passflag, ptree, stats)
            endif

            deallocate (mrange)
            deallocate (nrange)

            ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
            ! do ij=1,mm*rank_new
            ! j = (ij-1)/mm+1
            ! i = mod(ij-1,mm) + 1
            ! edge_m=header_m+i-1
            ! edge_n=blocks%ButterflySkel(level-1)%inds(index_i_loc_s1,index_j_loc_s1)%array(j)+header_n-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! blocks%ButterflyU%blocks(index_i_loc_k)%matrix(i,j)=ctemp
            ! enddo
            ! !$omp end taskloop

         else
            ! allocate (core(rankmax_c,rankmax_r))
            ! ! !$omp parallel
            ! ! !$omp single
            ! !$omp taskloop default(shared) private(ij,i,j,k,edge_m,edge_n,ctemp)
            ! do ij=1,rankmax_c*rankmax_r
            ! j = (ij-1)/rankmax_r+1
            ! i = mod(ij-1,rankmax_r) + 1
            ! edge_m=header_m+select_row(i)-1
            ! edge_n=header_n+select_column(j)-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! core(j,i)=ctemp
            ! enddo
            ! !$omp end taskloop
            ! ! !$omp end single
            ! ! !$omp end parallel

            allocate (matrix_V(rankmax_c, mm))
            allocate (matrix_V_tmp(mm, rankmax_c))
            allocate (mrange(mm))
            allocate (nrange(rankmax_c))
            do i = 1, mm
               mrange(i) = header_m + i - 1
            enddo
            do j = 1, rankmax_c
               nrange(j) = header_n + select_column(j) - 1
            enddo

            if (Nboundall > 0) then

               allocate (mnmap(mm, rankmax_c))
               allocate (mmap(mm))
               allocate (nmap(rankmax_c))
               allocate (mrange1(mm))
               allocate (nrange1(rankmax_c))
               mnmap = 1

               do i = 1, mm
                  group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, rankmax_c
                        if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                     enddo
                  endif
               enddo

               nrow = 0
               do i = 1, mm
                  if (sum(mnmap(i, :)) /= 0) then
                     nrow = nrow + 1
                     mmap(nrow) = i
                     mrange1(nrow) = mrange(i)
                  endif
               enddo
               ncol = 0
               do j = 1, rankmax_c
                  if (sum(mnmap(:, j)) /= 0) then
                     ncol = ncol + 1
                     nmap(ncol) = j
                     nrange1(ncol) = nrange(j)
                  endif
               enddo
               allocate (matrix_tmp(nrow, ncol))
               call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

               matrix_V_tmp = 0
               do i = 1, nrow
               do j = 1, ncol
                  matrix_V_tmp(mmap(i), nmap(j)) = matrix_tmp(i, j)
               enddo
               enddo
               deallocate (mnmap)
               deallocate (mmap)
               deallocate (nmap)
               deallocate (mrange1)
               deallocate (nrange1)
               deallocate (matrix_tmp)
            else
               call element_Zmn_block_user(mm, rankmax_c, mrange, nrange, matrix_V_tmp, msh, option, ker, 0, passflag, ptree, stats)
            endif

            deallocate (mrange)
            deallocate (nrange)
            call copymatT(matrix_V_tmp, matrix_V, mm, rankmax_c)
            deallocate (matrix_V_tmp)

            allocate (core(rankmax_c, rankmax_r))
            allocate (core_tmp(rankmax_r, rankmax_c))
            do i = 1, rankmax_r
               core(:, i) = matrix_V(:, select_row(i))
            enddo
            call copymatT(core, core_tmp, rankmax_c, rankmax_r)

            ! ! !$omp parallel
            ! ! !$omp single
            ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
            ! do ij=1,mm*rankmax_c
            ! j = (ij-1)/mm+1
            ! i = mod(ij-1,mm) + 1

            ! edge_m=header_m+i-1
            ! edge_n=header_n+select_column(j)-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! matrix_V(j,i)=ctemp
            ! enddo
            ! !$omp end taskloop
            ! ! !$omp end single
            ! ! !$omp end parallel

            allocate (jpvt(max(rankmax_c, rankmax_r)))
            allocate (tau(max(rankmax_c, rankmax_r)))
            jpvt = 0
            call geqp3modf90(core, jpvt, tau, option%tol_comp, SafeUnderflow, rank_new, flop=flop)
            flops = flops + flop

            if (rank_new > 0) then
               call un_or_mqrf90(core, tau, matrix_V, 'L', 'C', rankmax_c, mm, rank_new, flop=flop)
               flops = flops + flop
               call trsmf90(core, matrix_V, 'L', 'U', 'N', 'N', rank_new, mm, flop=flop)
               flops = flops + flop
            else
               rank_new = 1
               matrix_V = 0
            endif

            index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
            allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
            call copymatT(matrix_V(1:rank_new, 1:mm), blocks%ButterflyU%blocks(index_i_loc_k)%matrix, rank_new, mm)

            index_i_loc_s = (index_i - blocks%ButterflySkel(level_butterfly + 1)%idx_r)/blocks%ButterflySkel(level_butterfly + 1)%inc_r + 1
            allocate (blocks%ButterflySkel(level_butterfly + 1)%inds(index_i_loc_s, 1)%array(rank_new))
            do j = 1, rank_new
               blocks%ButterflySkel(level_butterfly + 1)%inds(index_i_loc_s, 1)%array(j) = select_row(jpvt(j))
            enddo

            if (option%sample_heuristic == 1) then
               allocate (core_tmp1(rank_new, rankmax_c))
               do i = 1, rank_new
                  core_tmp1(i, :) = core_tmp(jpvt(i), :)
               enddo
               jpvt = 0
               call geqp3modf90(core_tmp1, jpvt, tau, option%tol_comp, SafeUnderflow, Ncol_pre, flop=flop)
               flops = flops + flop
               Ncol_pre = rank_new
               do i = 1, Ncol_pre
                  select_col_pre(i) = header_n + select_column(jpvt(i)) - 1
               enddo
               deallocate (core_tmp1)
            endif

         endif

      elseif (level == 0) then
         index_j_loc_s = (index_j - blocks%ButterflySkel(level + 1)%idx_c)/blocks%ButterflySkel(level + 1)%inc_c + 1
         rank_new = size(blocks%ButterflySkel(level + 1)%inds(1, index_j_loc_s)%array)
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head
         index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1

         allocate (blocks%ButterflyV%blocks(index_j_loc_k)%matrix(nn, rank_new))
         allocate (matrix_V_tmp(rank_new, nn))
         allocate (mrange(rank_new))
         allocate (nrange(nn))
         do i = 1, rank_new
            mrange(i) = blocks%ButterflySkel(level + 1)%inds(1, index_j_loc_s)%array(i) + header_m - 1
         enddo
         do j = 1, nn
            nrange(j) = j + header_n - 1
         enddo

         if (Nboundall > 0) then

            allocate (mnmap(rank_new, nn))
            allocate (mmap(rank_new))
            allocate (nmap(nn))
            allocate (mrange1(rank_new))
            allocate (nrange1(nn))
            mnmap = 1

            do i = 1, rank_new
               group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
               group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
               if (group_n_mid /= -1) then
                  idxstart = msh%basis_group(group_n_mid)%head
                  idxend = msh%basis_group(group_n_mid)%tail
                  do j = 1, nn
                     if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                  enddo
               endif
            enddo

            nrow = 0
            do i = 1, rank_new
               if (sum(mnmap(i, :)) /= 0) then
                  nrow = nrow + 1
                  mmap(nrow) = i
                  mrange1(nrow) = mrange(i)
               endif
            enddo
            ncol = 0
            do j = 1, nn
               if (sum(mnmap(:, j)) /= 0) then
                  ncol = ncol + 1
                  nmap(ncol) = j
                  nrange1(ncol) = nrange(j)
               endif
            enddo
            allocate (matrix_tmp(nrow, ncol))
            call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

            matrix_V_tmp = 0
            do i = 1, nrow
            do j = 1, ncol
               matrix_V_tmp(mmap(i), nmap(j)) = matrix_tmp(i, j)
            enddo
            enddo
            deallocate (mnmap)
            deallocate (mmap)
            deallocate (nmap)
            deallocate (mrange1)
            deallocate (nrange1)
            deallocate (matrix_tmp)
         else
            call element_Zmn_block_user(rank_new, nn, mrange, nrange, matrix_V_tmp, msh, option, ker, 0, passflag, ptree, stats)
         endif

         deallocate (mrange)
         deallocate (nrange)
         call copymatT(matrix_V_tmp, blocks%ButterflyV%blocks(index_j_loc_k)%matrix, rank_new, nn)
         deallocate (matrix_V_tmp)

         ! ! !$omp parallel
         ! ! !$omp single
         ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
         ! do ij=1,rank_new*nn
         ! j = (ij-1)/rank_new+1
         ! i = mod(ij-1,rank_new) + 1
         ! edge_m=blocks%ButterflySkel(level+1)%inds(1,index_j_loc_s)%array(i)+header_m-1
         ! edge_n=j+header_n-1
         ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
         ! blocks%ButterflyV%blocks(index_j_loc_k)%matrix(j,i)=ctemp
         ! enddo
         ! !$omp end taskloop
         ! ! !$omp end single
         ! ! !$omp end parallel

      else
         index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
         index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

         header_n = msh%basis_group(group_n)%head
         header_m1 = msh%basis_group(group_m)%head
         header_m2 = msh%basis_group(2*group_m + 1)%head
         mmm1 = msh%basis_group(2*group_m)%tail - msh%basis_group(2*group_m)%head + 1

         if (level == level_final) then

            index_i_loc_s1 = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
            index_j_loc_s1 = (index_j - blocks%ButterflySkel(level - 1)%idx_c)/blocks%ButterflySkel(level - 1)%inc_c + 1

            rank_new = size(blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array)
            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(mm1, rank_new))
            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix(mm2, rank_new))

            allocate (matrix_V_tmp(mm1 + mm2, rank_new))
            allocate (mrange(mm1 + mm2))
            allocate (nrange(rank_new))
            do i = 1, mm1 + mm2
               if (i <= mm1) then
                  mrange(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(i) + header_m1 - 1
               else
                  mrange(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(i - mm1) + header_m2 - 1
               endif
            enddo
            do j = 1, rank_new
               nrange(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (mnmap(mm1 + mm2, rank_new))
               allocate (mmap(mm1 + mm2))
               allocate (nmap(rank_new))
               allocate (mrange1(mm1 + mm2))
               allocate (nrange1(rank_new))
               mnmap = 1

               do i = 1, mm1 + mm2
                  group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, rank_new
                        if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                     enddo
                  endif
               enddo

               nrow = 0
               do i = 1, mm1 + mm2
                  if (sum(mnmap(i, :)) /= 0) then
                     nrow = nrow + 1
                     mmap(nrow) = i
                     mrange1(nrow) = mrange(i)
                  endif
               enddo
               ncol = 0
               do j = 1, rank_new
                  if (sum(mnmap(:, j)) /= 0) then
                     ncol = ncol + 1
                     nmap(ncol) = j
                     nrange1(ncol) = nrange(j)
                  endif
               enddo
               allocate (matrix_tmp(nrow, ncol))
               call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

               matrix_V_tmp = 0
               do i = 1, nrow
               do j = 1, ncol
                  matrix_V_tmp(mmap(i), nmap(j)) = matrix_tmp(i, j)
               enddo
               enddo
               deallocate (mnmap)
               deallocate (mmap)
               deallocate (nmap)
               deallocate (mrange1)
               deallocate (nrange1)
               deallocate (matrix_tmp)
            else
               call element_Zmn_block_user(mm1 + mm2, rank_new, mrange, nrange, matrix_V_tmp, msh, option, ker, 0, passflag, ptree, stats)
            endif

            deallocate (mrange)
            deallocate (nrange)
            if (mm1 > 0) blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrix_V_tmp(1:mm1, :)
            if (mm2 > 0) blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix = matrix_V_tmp(1 + mm1:mm1 + mm2, :)
            deallocate (matrix_V_tmp)

            ! ! !$omp parallel
            ! ! !$omp single
            ! !$omp taskloop default(shared) private(ij,i,j,k,edge_m,edge_n,ctemp)
            ! do ij=1,mm*rank_new
            ! j = (ij-1)/mm+1
            ! i = mod(ij-1,mm) + 1
            ! if (i<=mm1) then
            ! edge_n=blocks%ButterflySkel(level-1)%inds(index_i_loc_s1,index_j_loc_s1)%array(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc,index_jj_loc)%array(i)+header_m1-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! blocks%ButterflyKerl(level)%blocks(index_i_loc_k,index_j_loc_k)%matrix(i,j)=ctemp
            ! else
            ! edge_n=blocks%ButterflySkel(level-1)%inds(index_i_loc_s1,index_j_loc_s1)%array(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc+1,index_jj_loc)%array(i-mm1)+header_m2-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1,index_j_loc_k)%matrix(i-mm1,j)=ctemp
            ! endif
            ! enddo
            ! !$omp end taskloop
            ! ! !$omp end single
            ! ! !$omp end parallel

         else
            index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
            index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1

            ! allocate (core(rankmax_c,rankmax_r))
            ! ! !$omp parallel
            ! ! !$omp single
            ! !$omp taskloop default(shared) private(ij,i,j,k,edge_m,edge_n,ctemp)
            ! do ij=1,rankmax_c*rankmax_r
            ! j = (ij-1)/rankmax_r+1
            ! i = mod(ij-1,rankmax_r) + 1
            ! if (select_row(i)<=mm1) then
            ! edge_n=select_column(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc,index_jj_loc)%array(select_row(i))+header_m1-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! core(j,i)=ctemp
            ! else
            ! edge_n=select_column(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc+1,index_jj_loc)%array(select_row(i)-mm1)+header_m2-1
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! core(j,i)=ctemp
            ! endif
            ! enddo
            ! !$omp end taskloop
            ! ! !$omp end single
            ! ! !$omp end parallel

            ! allocate (matrix_V(rankmax_c,mm))
            ! ! !$omp parallel
            ! ! !$omp single
            ! !$omp taskloop default(shared) private(ij,i,j,edge_m,edge_n,ctemp)
            ! do ij=1,mm*rankmax_c
            ! j = (ij-1)/mm+1
            ! i = mod(ij-1,mm) + 1
            ! if (i<=mm1) then
            ! edge_n=select_column(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc,index_jj_loc)%array(i)+header_m1-1
            ! else
            ! edge_n=select_column(j)+header_n-1
            ! edge_m=blocks%ButterflySkel(level+1)%inds(index_ii_loc+1,index_jj_loc)%array(i-mm1)+header_m2-1
            ! endif
            ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
            ! matrix_V(j,i)=ctemp
            ! enddo
            ! !$omp end taskloop
            ! ! !$omp end single
            ! ! !$omp end parallel

            allocate (matrix_V(rankmax_c, mm))
            allocate (matrix_V_tmp(mm, rankmax_c))
            allocate (mrange(mm))
            allocate (nrange(rankmax_c))
            do i = 1, mm
               if (i <= mm1) then
                  mrange(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(i) + header_m1 - 1
               else
                  mrange(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(i - mm1) + header_m2 - 1
               endif
            enddo
            do j = 1, rankmax_c
               nrange(j) = select_column(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (mnmap(mm, rankmax_c))
               allocate (mmap(mm))
               allocate (nmap(rankmax_c))
               allocate (mrange1(mm))
               allocate (nrange1(rankmax_c))
               mnmap = 1

               do i = 1, mm
                  group_m_mid = findgroup(mrange(i), msh, levelm, blocks%row_group)
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, rankmax_c
                        if (nrange(j) >= idxstart .and. nrange(j) <= idxend) mnmap(i, j) = 0
                     enddo
                  endif
               enddo

               nrow = 0
               do i = 1, mm
                  if (sum(mnmap(i, :)) /= 0) then
                     nrow = nrow + 1
                     mmap(nrow) = i
                     mrange1(nrow) = mrange(i)
                  endif
               enddo
               ncol = 0
               do j = 1, rankmax_c
                  if (sum(mnmap(:, j)) /= 0) then
                     ncol = ncol + 1
                     nmap(ncol) = j
                     nrange1(ncol) = nrange(j)
                  endif
               enddo
               allocate (matrix_tmp(nrow, ncol))
               call element_Zmn_block_user(nrow, ncol, mrange1, nrange1, matrix_tmp, msh, option, ker, 0, passflag, ptree, stats)

               matrix_V_tmp = 0
               do i = 1, nrow
               do j = 1, ncol
                  matrix_V_tmp(mmap(i), nmap(j)) = matrix_tmp(i, j)
               enddo
               enddo
               deallocate (mnmap)
               deallocate (mmap)
               deallocate (nmap)
               deallocate (mrange1)
               deallocate (nrange1)
               deallocate (matrix_tmp)
            else
               call element_Zmn_block_user(mm, rankmax_c, mrange, nrange, matrix_V_tmp, msh, option, ker, 0, passflag, ptree, stats)
            endif

            deallocate (mrange)
            deallocate (nrange)
            call copymatT(matrix_V_tmp, matrix_V, mm, rankmax_c)
            deallocate (matrix_V_tmp)
            allocate (core(rankmax_c, rankmax_r))
            allocate (core_tmp(rankmax_r, rankmax_c))
            do i = 1, rankmax_r
               core(:, i) = matrix_V(:, select_row(i))
            enddo
            call copymatT(core, core_tmp, rankmax_c, rankmax_r)

            allocate (jpvt(max(rankmax_c, rankmax_r)))
            allocate (tau(max(rankmax_c, rankmax_r)))
            jpvt = 0
            call geqp3modf90(core, jpvt, tau, option%tol_comp, SafeUnderflow, rank_new, flop=flop)
            flops = flops + flop

            if (rank_new > 0) then
               call un_or_mqrf90(core, tau, matrix_V, 'L', 'C', rankmax_c, mm, rank_new, flop=flop)
               flops = flops + flop
               call trsmf90(core, matrix_V, 'L', 'U', 'N', 'N', rank_new, mm, flop=flop)
               flops = flops + flop
            else
               rank_new = 1
               matrix_V = 0
            endif

            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(mm1, rank_new))
            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix(mm2, rank_new))

            call copymatT(matrix_V(1:rank_new, 1:mm1), blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, rank_new, mm1)
            call copymatT(matrix_V(1:rank_new, 1 + mm1:mm), blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, rank_new, mm - mm1)

            allocate (blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(rank_new))
            ! !$omp taskloop default(shared) private(j)
            do j = 1, rank_new
               if (select_row(jpvt(j)) <= mm1) then
                  blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(j) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(select_row(jpvt(j)))
               else
                  blocks%ButterflySkel(level)%inds(index_i_loc_s, index_j_loc_s)%array(j) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(select_row(jpvt(j)) - mm1) + mmm1
               endif
            enddo
            ! !$omp end taskloop

            if (option%sample_heuristic == 1) then
               allocate (core_tmp1(rank_new, rankmax_c))
               do i = 1, rank_new
                  core_tmp1(i, :) = core_tmp(jpvt(i), :)
               enddo
               jpvt = 0
               call geqp3modf90(core_tmp1, jpvt, tau, option%tol_comp, SafeUnderflow, Ncol_pre, flop=flop)
               flops = flops + flop
               Ncol_pre = rank_new
               do i = 1, Ncol_pre
                  select_col_pre(i) = header_n + select_column(jpvt(i)) - 1
               enddo
               deallocate (core_tmp1)
            endif

         endif

      endif

      if (allocated(core)) deallocate (core)
      if (allocated(core_tmp)) deallocate (core_tmp)
      if (allocated(tau)) deallocate (tau)
      if (allocated(jpvt)) deallocate (jpvt)
      if (allocated(matrix_V)) deallocate (matrix_V)
      if (allocated(select_row)) deallocate (select_row)
      if (allocated(select_column)) deallocate (select_column)
   end subroutine BF_compress_NlogN_oneblock_C

   subroutine Bplus_compress_N15(bplus, option, Memory, stats, msh, ker, ptree)

      use BPACK_DEFS

      use MISC_Utilities
      implicit none

      type(blockplus)::bplus
      integer:: ii, ll, bb, ierr
      real(kind=8) Memory, rtemp
      integer:: level_butterfly, level_BP, levelm, groupm_start, Nboundall, statflag
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree

      Memory = 0
      do ll = 1, bplus%Lplus
         bplus%LL(ll)%rankmax = 0
         statflag = 0
         if (ll == 1 .or. option%bp_cnt_lr == 1) statflag = 1  !!! only record the rank of the top-layer butterfly in a bplus
         do bb = 1, bplus%LL(ll)%Nbound
            if (IOwnPgrp(ptree, bplus%LL(ll)%matrices_block(bb)%pgno)) then
               if (bplus%LL(ll)%matrices_block(bb)%style == 1) then
                  call Full_construction(bplus%LL(ll)%matrices_block(bb), msh, ker, stats, option, ptree)
                  Memory = Memory + SIZEOF(bplus%LL(ll)%matrices_block(bb)%fullmat)/1024.0d3
               else

                  level_butterfly = bplus%LL(ll)%matrices_block(bb)%level_butterfly
                  level_BP = bplus%level

                  bplus%LL(ll)%matrices_block(bb)%level_half = BF_Switchlevel(bplus%LL(ll)%matrices_block(bb)%level_butterfly, option)
                  levelm = bplus%LL(ll)%matrices_block(bb)%level_half
                  groupm_start = bplus%LL(ll)%matrices_block(1)%row_group*2**levelm
                  Nboundall = 0
                  if (allocated(bplus%LL(ll + 1)%boundary_map)) Nboundall = size(bplus%LL(ll + 1)%boundary_map, 1)
                  if (option%forwardN15flag == 1) then
                     call BF_compress_N15(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                     call BF_sym2asym(bplus%LL(ll)%matrices_block(bb))
                  else
                     call BF_compress_NlogN(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                  end if
                  Memory = Memory + rtemp
                  bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax, bplus%LL(ll)%matrices_block(bb)%rankmax)
               endif
            endif
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE, bplus%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(bplus%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      end do

      ! !!!!!!! check error
      if (option%ErrFillFull == 1) call Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
      ! !!!!!!! check error

      ! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'

      return

   end subroutine Bplus_compress_N15

   subroutine BF_compress_N15(blocks, boundary_map, Nboundall, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)

      use BPACK_DEFS

      use MISC_Utilities
      implicit none

      integer Nboundall, statflag
      integer boundary_map(*)
      integer groupm_start

      integer inc_c, inc_r, nc, nr, idx_c, idx_r
      integer blocks_idx, i, j, level_butterfly, level_butterflyL, level_butterflyR, num_blocks, k, attempt
      integer group_m, group_n, mm, nn, mn, index_i, index_j, index_ij_loc, index_i_m, index_j_m, index_i_loc, index_j_loc, ii, jj, nn1, nn2, mm1, mm2, idxs_m, idxs_n
      integer level, levelm, level_half, level_final, level_loc, length_1, length_2, level_blocks, index_ij, edge_m, edge_n
      integer rank, rankmax, butterflyB_inuse, rank1, rank2, rmax
      real(kind=8) rate, tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
      real(kind=8) Memory, n1, n2
      DT ctemp
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      real(kind=8), allocatable :: Singular(:)
      integer flag, tmpi, tmpj
      real(kind=8):: a, b, error
      real(kind=8):: maxvalue(1:20)
      real(kind=8):: minvalue(1:20)
      integer dimension_m, dimension_n, dimension_rank, num_col, num_row
      type(matrixblock)::blocks
      integer cnt_tmp, rankFar, rankNear, leafsize, frow, passflag
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:)
      integer emptyflag
      cnt_tmp = 0
      rankFar = 0
      rankNear = 0
      ! ForwardSymmetricFlag = 1
      maxvalue = 0
      minvalue = 10000
      Memory = 0

      ! write(*,*)blocks%row_group,blocks%col_group,'In BF_compress_N15'

      level_blocks = blocks%level
      level_butterfly = blocks%level_butterfly
      !level_butterfly=Maxlevel-level_blocks
      ! level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2

!     if (Maxlevel-level_blocks<8) then
!         level_butterfly=Maxlevel-level_blocks
!     endif
      blocks%rankmax = -100000
      blocks%rankmin = 100000
      ! blocks%level_half = BF_Switchlevel(blocks%level_butterfly,option)
      blocks%level_half = floor_safe(dble(level_butterfly)/2d0)
      levelm = blocks%level_half
      level_half = blocks%level_half

      ! blocks%level_butterfly=level_butterfly

      num_blocks = 2**level_butterfly

      allocate (blocks%ButterflyU%blocks(num_blocks))
      allocate (blocks%ButterflyV%blocks(num_blocks))
      if (level_butterfly /= 0) then
         allocate (blocks%ButterflyKerl(level_butterfly))
      endif

      memory_butterfly = 0.

      if (level_butterfly == 0) then
         allocate (rankmax_for_butterfly(0:level_butterfly))
         rankmax_for_butterfly = 0
         allocate (rankmin_for_butterfly(0:level_butterfly))
         rankmin_for_butterfly = 10000

         group_m = blocks%row_group  ! Note: row_group and col_group interchanged here
         group_n = blocks%col_group

         mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
         nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1

         ! H-BACA
         leafsize = max(blocks%M, blocks%N)/option%LR_BLK_NUM

         if (allocated(blocks%ButterflyU%blocks(1)%matrix)) deallocate (blocks%ButterflyU%blocks(1)%matrix)
         if (allocated(blocks%ButterflyV%blocks(1)%matrix)) deallocate (blocks%ButterflyV%blocks(1)%matrix)
         call LR_HBACA(blocks, leafsize, rank, option, msh, ker, stats, ptree, blocks%pgno, 0)

         rankmax_for_butterfly(0) = max(blocks%rankmax, rankmax_for_butterfly(0))
         rankmin_for_butterfly(0) = min(blocks%rankmin, rankmin_for_butterfly(0))

         memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyV%blocks(1)%matrix)/1024.0d3
         memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyU%blocks(1)%matrix)/1024.0d3

      else if (level_butterfly == 1) then
         allocate (rankmax_for_butterfly(0:level_butterfly))
         rankmax_for_butterfly = 0
         allocate (rankmin_for_butterfly(0:level_butterfly))
         rankmin_for_butterfly = 10000

         do level = 0, level_butterfly
            index_ij = 0
            if (level > 0) then
               blocks%ButterflyKerl(level)%num_row = 2**level
               blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
               allocate (blocks%ButterflyKerl(level)%blocks(2**level, 2**(level_butterfly - level + 1)))
            endif
            if (level /= level_butterfly) then
               allocate (ButterflyP%blocks(2**(level + 1), 2**(level_butterfly - level)))
            end if

            do index_i = 1, 2**level
               do index_j = 1, 2**(level_butterfly - level)

                  if (level == 0) then
                     group_m = blocks%row_group  ! Note: row_group and col_group interchanged here
                     group_n = blocks%col_group
                     group_n = group_n*2**level_butterfly - 1 + index_j

                     mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                     nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                     allocate (QQ(mm, nn))
                     allocate (mrange(mm))
                     do ii = 1, mm
                        mrange(ii) = msh%basis_group(group_m)%head + ii - 1
                     enddo
                     allocate (nrange(nn))
                     do jj = 1, nn
                        nrange(jj) = msh%basis_group(group_n)%head + jj - 1
                     enddo
                     call element_Zmn_block_user(mm, nn, mrange, nrange, QQ, msh, option, ker, 0, passflag, ptree, stats)

                     ! do ii=1,mm
                     ! do jj =1,nn
                     ! edge_m = msh%basis_group(group_m)%head + ii - 1
                     ! edge_n = msh%basis_group(group_n)%head + jj - 1
                     ! call element_Zmn(edge_m,edge_n,ctemp,msh,option,ker)
                     ! QQ(ii,jj) = ctemp
                     ! end do
                     ! end do

                     mn = min(mm, nn)
                     allocate (UU(mm, mn), VV(mn, nn), Singular(mn))
                     call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, rank)

                     ! rank = 7
                     ! write(*,*)rank
                     rankmax_for_butterfly(0) = max(rank, rankmax_for_butterfly(0))
                     rankmin_for_butterfly(0) = min(rank, rankmin_for_butterfly(0))
                     blocks%rankmax = max(blocks%rankmax, rank)
                     blocks%rankmin = min(blocks%rankmin, rank)

                     allocate (mat_tmp(mm, rank))
                     !$omp parallel do default(shared) private(i,j,k,ctemp)
                     do j = 1, rank
                        do i = 1, mm
                           mat_tmp(i, j) = UU(i, j)*Singular(j)
                        enddo
                     enddo
                     !$omp end parallel do

                     mm1 = msh%basis_group(group_m*2)%tail - msh%basis_group(group_m*2)%head + 1
                     allocate (ButterflyP%blocks(2*index_i - 1, index_j)%matrix(mm1, rank))
                     allocate (ButterflyP%blocks(2*index_i, index_j)%matrix(mm - mm1, rank))

                     ButterflyP%blocks(2*index_i - 1, index_j)%matrix = mat_tmp(1:mm1, 1:rank)
                     ButterflyP%blocks(2*index_i, index_j)%matrix = mat_tmp(1 + mm1:mm, 1:rank)

                     allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn, rank))

                     !$omp parallel do default(shared) private(i,j)
                     do j = 1, rank
                        do i = 1, nn
                           blocks%ButterflyV%blocks(index_j)%matrix(i, j) = VV(j, i)
                           ! blocks%ButterflyV%blocks(index_j)%matrix(i,j)=random_dp_number()
                        enddo
                     enddo
                     !$omp end parallel do
                     deallocate (QQ, UU, VV, Singular, mat_tmp, mrange, nrange)

                  else

                     group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
                     group_n = blocks%col_group
                     group_m = group_m*2**level - 1 + index_i
                     group_n = group_n*2**(level_butterfly - level) - 1 + index_j

                     mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1

                     ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
                     if (size(ButterflyP_old%blocks(index_i, 2*index_j - 1)%matrix, 1) /= mm) then
                        write (*, *) 'mm incorrect'
                        stop
                     end if
                     nn1 = size(ButterflyP_old%blocks(index_i, 2*index_j - 1)%matrix, 2)
                     nn2 = size(ButterflyP_old%blocks(index_i, 2*index_j)%matrix, 2)
                     nn = nn1 + nn2

                     allocate (QQ(mm, nn))
                     !$omp parallel do default(shared) private(i,j)
                     do j = 1, nn1
                        do i = 1, mm
                           QQ(i, j) = ButterflyP_old%blocks(index_i, 2*index_j - 1)%matrix(i, j)
                        enddo
                     enddo
                     !$omp end parallel do
                     !$omp parallel do default(shared) private(i,j)
                     do j = 1, nn2
                        do i = 1, mm
                           QQ(i, j + nn1) = ButterflyP_old%blocks(index_i, 2*index_j)%matrix(i, j)
                        enddo
                     enddo
                     !$omp end parallel do

                     mn = min(mm, nn)
                     allocate (UU(mm, mn), VV(mn, nn), Singular(mn))
                     call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, rank)

                     ! rank = 7
                     ! write(*,*)rank
                     rankmax_for_butterfly(level) = max(rank, rankmax_for_butterfly(level))
                     rankmin_for_butterfly(level) = min(rank, rankmin_for_butterfly(level))
                     blocks%rankmax = max(blocks%rankmax, rank)
                     blocks%rankmin = min(blocks%rankmin, rank)

                     allocate (mat_tmp(mm, rank))
                     !$omp parallel do default(shared) private(i,j,k,ctemp)
                     do j = 1, rank
                        do i = 1, mm
                           mat_tmp(i, j) = UU(i, j)*Singular(j)
                        enddo
                     enddo
                     !$omp end parallel do

                     if (level /= level_butterfly) then
                        mm1 = msh%basis_group(group_m*2)%tail - msh%basis_group(group_m*2)%head + 1
                        allocate (ButterflyP%blocks(2*index_i - 1, index_j)%matrix(mm1, rank))
                        allocate (ButterflyP%blocks(2*index_i, index_j)%matrix(mm - mm1, rank))

                        ButterflyP%blocks(2*index_i - 1, index_j)%matrix = mat_tmp(1:mm1, 1:rank)
                        ButterflyP%blocks(2*index_i, index_j)%matrix = mat_tmp(1 + mm1:mm, 1:rank)
                     else
                        allocate (blocks%ButterflyU%blocks(index_i)%matrix(mm, rank))
                        !$omp parallel do default(shared) private(i,j)
                        do j = 1, rank
                           do i = 1, mm
                              blocks%ButterflyU%blocks(index_i)%matrix(i, j) = mat_tmp(i, j)
                              ! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
                           enddo
                        enddo
                        !$omp end parallel do

                     end if

                     allocate (blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix(rank, nn1))
                     allocate (blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix(rank, nn2))

                     !$omp parallel do default(shared) private(i,j)
                     do i = 1, rank
                        do j = 1, nn1
                           blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix(i, j) = VV(i, j)
                           ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
                        enddo
                     enddo
                     !$omp end parallel do
                     !$omp parallel do default(shared) private(i,j)
                     do i = 1, rank
                        do j = 1, nn2
                           blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix(i, j) = VV(i, j + nn1)
                           ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
                        enddo
                     enddo
                     !$omp end parallel do

                     deallocate (QQ, UU, VV, Singular, mat_tmp)

                  end if

                  index_ij = index_ij + 1
                  if (level == 0) then
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyV%blocks(index_ij)%matrix)/1024.0d3
                  else
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix)/1024.0d3
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix)/1024.0d3
                  endif
                  if (level == level_butterfly) then
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyU%blocks(index_ij)%matrix)/1024.0d3
                  endif
               enddo
            enddo

            if (level /= level_butterfly) then
               if (allocated(ButterflyP_old%blocks)) then
                  do ii = 1, 2**(level)
                     do jj = 1, 2**(level_butterfly - level + 1)
                        deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                     end do
                  end do
                  deallocate (ButterflyP_old%blocks)
               end if
               allocate (ButterflyP_old%blocks(2**(level + 1), 2**(level_butterfly - level)))
               do ii = 1, 2**(level + 1)
                  do jj = 1, 2**(level_butterfly - level)
                     mm = size(ButterflyP%blocks(ii, jj)%matrix, 1)
                     nn = size(ButterflyP%blocks(ii, jj)%matrix, 2)
                     allocate (ButterflyP_old%blocks(ii, jj)%matrix(mm, nn))
                     ButterflyP_old%blocks(ii, jj)%matrix = ButterflyP%blocks(ii, jj)%matrix
                     deallocate (ButterflyP%blocks(ii, jj)%matrix)
                  end do
               end do
               deallocate (ButterflyP%blocks)
            else
               if (allocated(ButterflyP_old%blocks)) then
                  do ii = 1, 2**(level)
                     do jj = 1, 2**(level_butterfly - level + 1)
                        deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                     end do
                  end do
                  deallocate (ButterflyP_old%blocks)
               end if
            end if
         enddo
      else

         do level = 1, level_butterfly
            blocks%ButterflyKerl(level)%num_row = 2**level
            blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
            allocate (blocks%ButterflyKerl(level)%blocks(2**level, 2**(level_butterfly - level + 1)))
         end do

         level_butterflyL = level_butterfly - levelm
         level_butterflyR = level_butterfly - level_butterflyL

         allocate (rankmax_for_butterfly(-level_butterflyR:level_butterflyL))
         rankmax_for_butterfly = 0
         allocate (rankmin_for_butterfly(-level_butterflyR:level_butterflyL))
         rankmin_for_butterfly = 100000

         allocate (blocks%ButterflyMiddle(2**levelm, 2**(level_butterfly - levelm)))

         ! assign block indices, this only works in sequential

         ! construct the middle level and the left half
         do index_i_m = 1, 2**levelm

            level_loc = 0
            index_i_loc = 1
            allocate (ButterflyP_old%blocks(2**(level_loc + 1), 2**(level_butterflyL - level_loc)))

            do index_j_m = 1, 2**(level_butterfly - levelm)

               index_j_loc = index_j_m
               group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
               group_n = blocks%col_group
               group_m = group_m*2**levelm - 1 + index_i_m
               group_n = group_n*2**(level_butterfly - levelm) - 1 + index_j_m

               mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
               nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
               idxs_m = msh%basis_group(group_m)%head
               idxs_n = msh%basis_group(group_n)%head

               rmax = min(option%rmax, min(mm, nn))
               allocate (matU(mm, rmax))
               allocate (matV(rmax, nn))
               allocate (Singular(rmax))

               emptyflag = 0
               if (Nboundall > 0) then
                  if (boundary_map(group_m - groupm_start + 1) == group_n) emptyflag = 1
               endif

               if (emptyflag == 1) then
                  rank = 1
                  Singular(1:rank) = 0
                  rankmax_for_butterfly(level_loc) = max(rank, rankmax_for_butterfly(level_loc))
                  rankmin_for_butterfly(level_loc) = min(rank, rankmin_for_butterfly(level_loc))

                  allocate (blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(rank, rank))
                  blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix = 0
                  do ii = 1, rank
                     blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(ii, ii) = 0d0
                  end do
               else
                  frow = 1
                  call LR_ACA(matU, matV, Singular, idxs_m, idxs_n, mm, nn, frow, rmax, rank, option%tol_comp*0.1, option%tol_comp, msh, ker, stats, ptree, option, error)
                  rankmax_for_butterfly(level_loc) = max(rank, rankmax_for_butterfly(level_loc))
                  rankmin_for_butterfly(level_loc) = min(rank, rankmin_for_butterfly(level_loc))

                  allocate (blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(rank, rank))
                  blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix = 0
                  do ii = 1, rank
                     blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(ii, ii) = 1d0/Singular(ii)
                  end do
               end if

               blocks%rankmax = max(blocks%rankmax, rank)
               blocks%rankmin = min(blocks%rankmin, rank)

               allocate (mat_tmp(mm, rank))
               !$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, rank
                  do i = 1, mm
                     mat_tmp(i, j) = matU(i, j)*Singular(j)
                  enddo
               enddo
               !$omp end parallel do

               mm1 = msh%basis_group(group_m*2)%tail - msh%basis_group(group_m*2)%head + 1
               allocate (ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix(mm1, rank))
               allocate (ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix(mm - mm1, rank))

               ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix = mat_tmp(1:mm1, 1:rank)
               ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix = mat_tmp(1 + mm1:mm, 1:rank)

               deallocate (matU, matV, Singular, mat_tmp)

            end do

            n1 = OMP_get_wtime()

            do level_loc = 1, level_butterflyL
               level = level_loc + levelm

               if (level_loc /= level_butterflyL) then
                  allocate (ButterflyP%blocks(2**(level_loc + 1), 2**(level_butterflyL - level_loc)))
               end if

               ! do index_i_loc=1, 2**level_loc
               ! do index_j_loc=1, 2**(level_butterflyL-level_loc)

               ! write(*,*)'addaa111111'

               !$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
               do index_ij_loc = 1, 2**level_butterflyL
                  index_j_loc = mod(index_ij_loc - 1, 2**(level_butterflyL - level_loc)) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL - level_loc)))

                  call LocalButterflySVD_Left(index_i_loc, index_j_loc, level_loc, level_butterflyL, level, index_i_m, blocks, option, msh, ButterflyP_old, ButterflyP)
               enddo
               !$omp end parallel do

               ! write(*,*)'addaa1111112222'

               do index_ij_loc = 1, 2**level_butterflyL
                  index_j_loc = mod(index_ij_loc - 1, 2**(level_butterflyL - level_loc)) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL - level_loc)))
                  index_j = index_j_loc
                  index_i = (index_i_m - 1)*(2**level_loc) + index_i_loc

                  rank = size(blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 1)

                  rankmax_for_butterfly(level_loc) = max(rank, rankmax_for_butterfly(level_loc))
                  rankmin_for_butterfly(level_loc) = min(rank, rankmin_for_butterfly(level_loc))

                  blocks%rankmax = max(blocks%rankmax, rank)
                  blocks%rankmin = min(blocks%rankmin, rank)
                  memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix)/1024.0d3
                  memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix)/1024.0d3
                  if (level_loc == level_butterflyL) then
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyU%blocks(index_i)%matrix)/1024.0d3
                  endif
               enddo

               if (level_loc /= level_butterflyL) then
                  if (allocated(ButterflyP_old%blocks)) then
                     do ii = 1, 2**(level_loc)
                        do jj = 1, 2**(level_butterflyL - level_loc + 1)
                           deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                        end do
                     end do
                     deallocate (ButterflyP_old%blocks)
                  end if
                  allocate (ButterflyP_old%blocks(2**(level_loc + 1), 2**(level_butterflyL - level_loc)))
                  do ii = 1, 2**(level_loc + 1)
                     do jj = 1, 2**(level_butterflyL - level_loc)
                        mm = size(ButterflyP%blocks(ii, jj)%matrix, 1)
                        nn = size(ButterflyP%blocks(ii, jj)%matrix, 2)
                        allocate (ButterflyP_old%blocks(ii, jj)%matrix(mm, nn))
                        ButterflyP_old%blocks(ii, jj)%matrix = ButterflyP%blocks(ii, jj)%matrix
                        deallocate (ButterflyP%blocks(ii, jj)%matrix)
                     end do
                  end do
                  deallocate (ButterflyP%blocks)
               else
                  if (allocated(ButterflyP_old%blocks)) then
                     do ii = 1, 2**(level_loc)
                        do jj = 1, 2**(level_butterflyL - level_loc + 1)
                           deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                        end do
                     end do
                     deallocate (ButterflyP_old%blocks)
                  end if
               end if

            end do
            n2 = OMP_get_wtime()
            ! time_tmp = time_tmp + n2 - n1

         enddo

         ! write(*,*)'stats:',2**levelm,2**(level_butterfly-levelm), cnt_tmp,rankFar,rankNear

         ! construct the the right half
         do index_j_m = 1, 2**(level_butterfly - levelm)
            level_loc = 0
            index_j_loc = 1
            allocate (ButterflyP_old%blocks(2**(level_butterflyR - level_loc), 2**(level_loc + 1)))
            do index_i_m = 1, 2**levelm
               index_i_loc = index_i_m
               group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
               group_n = blocks%col_group
               group_m = group_m*2**levelm - 1 + index_i_m
               group_n = group_n*2**(level_butterfly - levelm) - 1 + index_j_m

               mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
               nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
               idxs_m = msh%basis_group(group_m)%head
               idxs_n = msh%basis_group(group_n)%head

               rmax = min(option%rmax, min(mm, nn))
               allocate (matU(mm, rmax))
               allocate (matV(rmax, nn))
               allocate (Singular(rmax))

               emptyflag = 0
               if (Nboundall > 0) then
                  if (boundary_map(group_m - groupm_start + 1) == group_n) emptyflag = 1
               endif

               if (emptyflag == 1) then
                  ! if(.not. near_or_far(group_m,group_n,2d0))then
                  rank = 1
                  Singular(1:rank) = 0
                  matV(1:rank, 1:nn) = 0
                  ! if(blocks==342)write(111,*)Singular(1:rank)
               else
                  frow = 1
                  call LR_ACA(matU, matV, Singular, idxs_m, idxs_n, mm, nn, frow, rmax, rank, option%tol_comp*0.1, option%tol_comp, msh, ker, stats, ptree, option, error)
               end if

               ! rank = min(rank,37)

               ! if(.not. near_or_far(group_m,group_n,2d0))then
               ! rank = 1
               ! Singular(1:rank) = 0
               ! ! if(blocks==342)write(111,*)Singular(1:rank)

               ! end if

               allocate (mat_tmp(rank, nn))
               !$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, nn
                  do i = 1, rank
                     mat_tmp(i, j) = matV(i, j)*Singular(i)
                  enddo
               enddo
               !$omp end parallel do

               nn1 = msh%basis_group(group_n*2)%tail - msh%basis_group(group_n*2)%head + 1
               allocate (ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix(rank, nn1))
               allocate (ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix(rank, nn - nn1))

               ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix = mat_tmp(1:rank, 1:nn1)
               ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix = mat_tmp(1:rank, 1 + nn1:nn)

               deallocate (matU, matV, Singular, mat_tmp)
            end do

            n1 = OMP_get_wtime()
            do level_loc = 1, level_butterflyR
               level = levelm + 1 - level_loc

               if (level_loc /= level_butterflyR) then
                  allocate (ButterflyP%blocks(2**(level_butterflyR - level_loc), 2**(level_loc + 1)))
               end if

               !$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
               do index_ij_loc = 1, 2**level_butterflyR
                  index_j_loc = mod(index_ij_loc - 1, 2**level_loc) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
                  call LocalButterflySVD_Right(index_i_loc, index_j_loc, level_loc, level_butterflyR, level, level_butterfly, index_j_m, blocks, option, msh, ButterflyP_old, ButterflyP)
               enddo
               !$omp end parallel do

               do index_ij_loc = 1, 2**level_butterflyR
                  index_j_loc = mod(index_ij_loc - 1, 2**level_loc) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))

                  index_i = index_i_loc
                  index_j = (index_j_m - 1)*(2**level_loc) + index_j_loc

                  rank = size(blocks%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, 2)
                  rankmax_for_butterfly(-level_loc) = max(rank, rankmax_for_butterfly(-level_loc))
                  rankmin_for_butterfly(-level_loc) = min(rank, rankmin_for_butterfly(-level_loc))
                  blocks%rankmax = max(blocks%rankmax, rank)
                  blocks%rankmin = min(blocks%rankmin, rank)

                  memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix)/1024.0d3
                  memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix)/1024.0d3
                  if (level_loc == level_butterflyR) then
                     memory_butterfly = memory_butterfly + SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
                  endif
               enddo

               ! do index_i_loc=1, 2**(level_butterflyR-level_loc)
               ! do index_j_loc=1,2**level_loc

               ! index_i = index_i_loc
               ! index_j = (index_j_m-1)*(2**level_loc)+index_j_loc
               ! group_n=blocks%col_group   ! Note: row_group and col_group interchanged here
               ! group_n=group_n*2**(level_butterfly-level+1)-1+index_j
               ! nn=msh%basis_group(group_n)%tail-msh%basis_group(group_n)%head+1

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
               ! nn1 = msh%basis_group(group_n*2)%tail-msh%basis_group(group_n*2)%head+1
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

               if (level_loc /= level_butterflyR) then
                  if (allocated(ButterflyP_old%blocks)) then
                     do ii = 1, 2**(level_butterflyR - level_loc + 1)
                        do jj = 1, 2**(level_loc)
                           deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                        end do
                     end do
                     deallocate (ButterflyP_old%blocks)
                  end if
                  allocate (ButterflyP_old%blocks(2**(level_butterflyR - level_loc), 2**(level_loc + 1)))
                  do ii = 1, 2**(level_butterflyR - level_loc)
                     do jj = 1, 2**(level_loc + 1)
                        mm = size(ButterflyP%blocks(ii, jj)%matrix, 1)
                        nn = size(ButterflyP%blocks(ii, jj)%matrix, 2)
                        allocate (ButterflyP_old%blocks(ii, jj)%matrix(mm, nn))
                        ButterflyP_old%blocks(ii, jj)%matrix = ButterflyP%blocks(ii, jj)%matrix
                        deallocate (ButterflyP%blocks(ii, jj)%matrix)
                     end do
                  end do
                  deallocate (ButterflyP%blocks)
               else
                  if (allocated(ButterflyP_old%blocks)) then
                     do ii = 1, 2**(level_butterflyR - level_loc + 1)
                        do jj = 1, 2**(level_loc)
                           deallocate (ButterflyP_old%blocks(ii, jj)%matrix)
                        end do
                     end do
                     deallocate (ButterflyP_old%blocks)
                  end if
               end if

            end do

            n2 = OMP_get_wtime()
            ! time_tmp = time_tmp + n2 - n1
         enddo

      end if

      do level = 0, level_half
         call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
         if (level == 0) then
            blocks%ButterflyV%idx = idx_c
            blocks%ButterflyV%inc = inc_c
            blocks%ButterflyV%nblk_loc = nc
         elseif (level == level_butterfly + 1) then
            blocks%ButterflyU%idx = idx_r
            blocks%ButterflyU%inc = inc_r
            blocks%ButterflyU%nblk_loc = nr
         else
            blocks%ButterflyKerl(level)%num_row = 2**level
            blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
            blocks%ButterflyKerl(level)%idx_r = idx_r
            blocks%ButterflyKerl(level)%inc_r = inc_r
            blocks%ButterflyKerl(level)%nr = nr
            blocks%ButterflyKerl(level)%idx_c = idx_c*2 - 1
            blocks%ButterflyKerl(level)%inc_c = inc_c
            blocks%ButterflyKerl(level)%nc = nc*2
         endif
      end do

      level_final = level_half + 1
      do level = level_butterfly + 1, level_final, -1
         call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
         if (level == 0) then
            blocks%ButterflyV%idx = idx_c
            blocks%ButterflyV%inc = inc_c
            blocks%ButterflyV%nblk_loc = nc
         elseif (level == level_butterfly + 1) then
            blocks%ButterflyU%idx = idx_r
            blocks%ButterflyU%inc = inc_r
            blocks%ButterflyU%nblk_loc = nr
         else
            blocks%ButterflyKerl(level)%num_row = 2**level
            blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
            blocks%ButterflyKerl(level)%idx_r = idx_r*2 - 1
            blocks%ButterflyKerl(level)%inc_r = inc_r
            blocks%ButterflyKerl(level)%nr = nr*2
            blocks%ButterflyKerl(level)%idx_c = idx_c
            blocks%ButterflyKerl(level)%inc_c = inc_c
            blocks%ButterflyKerl(level)%nc = nc
         endif
      end do

      ! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
      if (statflag == 1) then
         if (allocated(stats%rankmax_of_level)) stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly), stats%rankmax_of_level(level_blocks))
      endif
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

      Memory = memory_butterfly
      !write (*,*) memory_butterfly
      !pause

      return

   end subroutine BF_compress_N15

   subroutine BF_compress_test(blocks, msh, ker, element_Zmn, ptree, option, stats)

      use BPACK_DEFS
      use BPACK_Utilities
      implicit none

      type(matrixblock) :: blocks
      real(kind=8) a, b, error, v1, v2, v3
      integer i, j, k, ii, jj, iii, jjj, kk, group_m, group_n, mm, nn, mi, nj, head_m, head_n, Dimn, edge_m, edge_n
      DT value1, value2, ctemp1, ctemp2
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :)
      type(kernelquant)::ker
      type(mesh)::msh
      procedure(Zelem)::element_Zmn
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer, allocatable::order_m(:), order_n(:)
      real(kind=8), allocatable::distance_m(:), distance_n(:), center(:)

      ctemp1 = 1.0d0; ctemp2 = 0.0d0

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
      ! ! write(*,*)blocks%level,h_mat%Maxlevel
      ! ! write(*,*)'h22'

      ! if(allocated(blocks%fullmat))then
      ! ! write(*,*)'h3'
      ! call Full_block_MVP_dat(blocks,'N',mm,1,Vin,Vout1,ctemp1,ctemp2)
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

      allocate (order_m(blocks%M))
      allocate (order_n(blocks%N))
      do ii = 1, min(blocks%M, blocks%N)
         call random_number(a)
         call random_number(b)
         order_m(ii) = floor_safe(a*(mm - 1)) + 1
         order_n(ii) = floor_safe(b*(nn - 1)) + 1

         ! order_m(ii)=ii
         ! order_n(ii)=ii
      enddo

      ! !!!!! The following picks the close points first, can be commented out if geometry info is not available
      ! allocate(distance_m(blocks%M))
      ! distance_m=Bigvalue
      ! allocate(distance_n(blocks%N))
      ! distance_n=Bigvalue

      ! Dimn = 0
      ! if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
      ! allocate(center(Dimn))
      ! center = 0
      ! head_n = blocks%headn
      ! do j=1,blocks%N
      ! edge_n = head_n-1+j
      ! center = center + msh%xyz(1:Dimn,msh%new2old(edge_n))
      ! enddo
      ! center = center/blocks%N
      ! head_m = blocks%headm
      ! do i=1,blocks%M
      ! edge_m=head_m-1+i
      ! distance_m(i) = sum((msh%xyz(1:Dimn,msh%new2old(edge_m))-center(1:Dimn))**2d0)
      ! enddo
      ! deallocate(center)
      ! call quick_sort(distance_m,order_m,blocks%M)
      ! deallocate(distance_m)

      ! Dimn = 0
      ! if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
      ! allocate(center(Dimn))
      ! center = 0
      ! head_m = blocks%headm
      ! do i=1,blocks%M
      ! edge_m = head_m-1+i
      ! center = center + msh%xyz(1:Dimn,msh%new2old(edge_m))
      ! enddo
      ! center = center/blocks%M
      ! head_n = blocks%headn
      ! do j=1,blocks%N
      ! edge_n=head_n-1+j
      ! distance_n(j) = sum((msh%xyz(1:Dimn,msh%new2old(edge_n))-center(1:Dimn))**2d0)
      ! enddo
      ! deallocate(center)
      ! call quick_sort(distance_n,order_n,blocks%N)
      ! deallocate(distance_n)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      v1 = 0
      v2 = 0
      v3 = 0
      ! do i=1,min(mm,nn)
      do i = 1, 100
         do j = 1, 100

            mi = order_m(i)
            nj = order_n(j)

            ! iii=int((mi+1)/2)+msh%basis_group(group_m)%head-1
            ! jjj=int((nj+1)/2)+msh%basis_group(group_n)%head-1
            ! ii=2-mod(mi,2)
            ! jj=2-mod(nj,2)
            ! call element_Zmn(iii,jjj,ii,jj,value1)

            call element_Zmn(mi + head_m - 1, nj + head_n - 1, value1, msh, option, ker)

            call BF_value(mi, nj, blocks, value2)
            v1 = v1 + abs(value1)**2d0
            v2 = v2 + abs(value2)**2d0
            v3 = v3 + abs(value2 - value1)**2d0
            ! if(abs(value1)>SafeUnderflow)write (*,*) abs(value1), abs(value2) !, abs(value1-value2)/abs(value1)
         enddo
      enddo

      deallocate (order_m)
      deallocate (order_n)

      write (*, *) 'partial fnorm:', v1, v2, sqrt(v3/v1), ' rank: ', blocks%rankmax, ' level: ', blocks%level_butterfly

      return

   end subroutine BF_compress_test

   subroutine LR_HBACA(blocks, leafsize, rank, option, msh, ker, stats, ptree, pgno, cridx)
      use BPACK_DEFS
      implicit none
      integer rank, ranktmp, leafsize
      integer header_m, header_n
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, mn
      type(mesh)::msh
      type(Hoption)::option
      type(kernelquant)::ker
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno
      integer:: cridx, info
      integer::mrange_dummy(1), nrange_dummy(1)
      type(Hstat)::stats
      integer::passflag = 0
      integer::frow, rmax
      real(kind=8)::error
      DT:: mat_dummy(1, 1)

      call LR_HBACA_Leaflevel(blocks, leafsize, rank, option, msh, ker, stats, ptree, pgno, cridx)

      passflag = 0
      do while (passflag == 0)
         call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
      enddo

      call LR_HMerge(blocks, rank, option, msh, stats, ptree, pgno, cridx, 1)

   end subroutine LR_HBACA

   recursive subroutine LR_HMerge(blocks, rank, option, msh, stats, ptree, pgno, cridx, hbacaflag)
      use BPACK_DEFS
      implicit none
      integer rank, ranktmp
      integer header_m, header_n
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, rmax, mn
      type(mesh)::msh
      type(Hoption)::option
      type(matrixblock)::blocks, blockc(2)
      type(proctree)::ptree
      integer pgno, pgno1, pgno2
      integer:: cridx, info
      DT::tmp(1, 1)
      DT, allocatable:: UU(:, :), VV(:, :), matU(:, :), matV(:, :), matU1(:, :), matV1(:, :), matU2(:, :), matV2(:, :), matU1D(:, :), matV1D(:, :), Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Fullmat(:, :), QQ1(:, :), matU2D(:, :), matV2D(:, :)
      real(kind=8), allocatable::Singular(:)
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, M1, N1, M2, N2, rank1, rank2, ierr
      integer::descsmatU(9), descsmatV(9), descsmatU1(9), descsmatV1(9), descsmatU2(9), descsmatV2(9), descUU(9), descVV(9), descsmatU1c(9), descsmatU2c(9), descsmatV1c(9), descsmatV2c(9), descButterflyV(9), descButterflyU(9), descButterU1D(9), descButterV1D(9), descVin(9), descVout(9), descVinter(9), descFull(9)
      integer dims(6), dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
      DT:: TEMP(1)
      integer LWORK, mnmax, mnmin, rank_new
      DT, allocatable:: WORK(:)
      real(kind=8), allocatable::RWORK(:), center(:)
      real(kind=8):: rtemp, dist, error, rtemp1, rtemp0, fnorm1, fnorm0, flop
      integer :: nb1Dc, nb1Dr, frow, Dimn, edge_n, edge_m, MyID, Ntest, rmaxc, rmaxr
      integer, allocatable::M_p(:, :), N_p(:, :), mrange(:), nrange(:)
      type(Hstat)::stats
      integer::passflag = 0
      integer::Maxgrp
      integer::hbacaflag
      logical::built = .false.
      Maxgrp = 2**(ptree%nlevel) - 1

      rank = 0
      blocks%ButterflyU%idx = 1
      blocks%ButterflyU%inc = 1
      blocks%ButterflyU%nblk_loc = 1
      blocks%ButterflyV%idx = 1
      blocks%ButterflyV%inc = 1
      blocks%ButterflyV%nblk_loc = 1

      built = (hbacaflag == 1 .and. option%RecLR_leaf == ACANMERGE)
      if ((.not. (associated(blocks%sons))) .and. (.not. built)) then !  reach bottom level
         ! !!!!!!! check error
      else
         if (built) then  ! no need to do merge as LR is already built in parallel
            rank = blocks%rankmax
            goto 101
         endif

         if (pgno*2 > Maxgrp .or. hbacaflag /= 1) then
            pgno1 = pgno
            pgno2 = pgno
         else
            pgno1 = pgno*2
            pgno2 = pgno*2 + 1
         endif

         if (IOwnPgrp(ptree, pgno)) then

            call blacs_gridinfo(ptree%pgrp(pgno1)%ctxt, nprow1, npcol1, myrow1, mycol1)
            call blacs_gridinfo(ptree%pgrp(pgno2)%ctxt, nprow2, npcol2, myrow2, mycol2)

            dims_tmp(1:3) = 0
            ! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow1,npcol1,myrow1,mycol1,'dddd'
            if (IOwnPgrp(ptree, pgno1)) then
               call LR_HMerge(blocks%sons(1, 1), rank, option, msh, stats, ptree, pgno1, cridx + 1, hbacaflag)
               dims_tmp(1) = blocks%sons(1, 1)%M
               dims_tmp(2) = blocks%sons(1, 1)%N
               dims_tmp(3) = blocks%sons(1, 1)%rankmax
            endif

            dims_tmp(4:6) = 0

            ! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow2,npcol2,myrow2,mycol2,'dddd2'
            if (IOwnPgrp(ptree, pgno2)) then
               call LR_HMerge(blocks%sons(2, 1), rank, option, msh, stats, ptree, pgno2, cridx + 1, hbacaflag)
               dims_tmp(4) = blocks%sons(2, 1)%M
               dims_tmp(5) = blocks%sons(2, 1)%N
               dims_tmp(6) = blocks%sons(2, 1)%rankmax
            endif
            ! write(*,*)ptree%MyID,cridx+1
            call MPI_ALLREDUCE(dims_tmp, dims, 6, MPI_INTEGER, &
                               MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

            M1 = dims(1)
            N1 = dims(2)
            rank1 = dims(3)
            M2 = dims(4)
            N2 = dims(5)
            rank2 = dims(6)

            call blacs_gridinfo(ptree%pgrp(pgno)%ctxt, nprow, npcol, myrow, mycol)

            if (mod(cridx + 1, 2) == 0) then  ! merge along column dimension

               call assert(M1 == M2, 'M1/=M2 in column merge')

               if (nprow /= -1 .and. npcol /= -1) then
                  myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (matU(myArows, myAcols))
                  matU = 0
                  call descinit(descsmatU, M1, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU')

                  myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
                  allocate (matV1(myArows, myAcols))
                  call descinit(descsmatV1, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV1')

                  myArows = numroc_wp(N2, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
                  allocate (matV2(myArows, myAcols))
                  call descinit(descsmatV2, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV2')
               else
                  descsmatU(2) = -1
                  descsmatV1(2) = -1
                  descsmatV2(2) = -1
                  allocate (matU(1, 1))
                  matU = 0
                  allocate (matV1(1, 1))
                  allocate (matV2(1, 1))
               endif

               if (nprow1 /= -1 .and. npcol1 /= -1) then
                  ! redistribute U1
                  myArows = numroc_wp(M1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit(descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU1c')

                  call pgemr2df90(M1, rank1, blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V1
                  myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit(descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV1c')
                  call pgemr2df90(N1, rank1, blocks%sons(1, 1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyV%blocks(1)%matrix)
               else
                  descsmatU1c(2) = -1
                  call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  descsmatV1c(2) = -1
                  call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, ptree%pgrp(pgno)%ctxt1DCol)
               endif
               if (allocated(blocks%sons(1, 1)%ButterflyU%blocks)) deallocate (blocks%sons(1, 1)%ButterflyU%blocks)
               if (allocated(blocks%sons(1, 1)%ButterflyV%blocks)) deallocate (blocks%sons(1, 1)%ButterflyV%blocks)

               if (nprow2 /= -1 .and. npcol2 /= -1) then
                  ! redistribute U2
                  myArows = numroc_wp(M2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit(descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU2c')
                  call pgemr2df90(M2, rank2, blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU, 1, 1 + rank1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V2
                  myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit(descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV2c')
                  call pgemr2df90(N2, rank2, blocks%sons(2, 1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyV%blocks(1)%matrix)
               else
                  descsmatU2c(2) = -1
                  call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU, 1, 1 + rank1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  descsmatV2c(2) = -1
                  call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, ptree%pgrp(pgno)%ctxt1DCol)
               endif

               if (allocated(blocks%sons(2, 1)%ButterflyU%blocks)) deallocate (blocks%sons(2, 1)%ButterflyU%blocks)
               if (allocated(blocks%sons(2, 1)%ButterflyV%blocks)) deallocate (blocks%sons(2, 1)%ButterflyV%blocks)

               if (nprow /= -1 .and. npcol /= -1) then
                  ! compute truncated SVD on matU
                  mnmin = min(M1, rank1 + rank2)

                  myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
                  allocate (UU(myArows, myAcols))
                  call descinit(descUU, M1, mnmin, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descUU')
                  UU = 0
                  myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (VV(myArows, myAcols))
                  call descinit(descVV, mnmin, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descVV')
                  VV = 0

                  allocate (Singular(mnmin))
                  Singular = 0

                  call PSVD_Truncate(M1, rank1 + rank2, matU, descsmatU, UU, VV, descUU, descVV, Singular, option%tol_comp, rank, ptree%pgrp(pgno)%ctxt, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

                  do ii = 1, rank
                     call g2l(ii, rank, nprow, nbslpk, iproc, myi)
                     if (iproc == myrow) then
                        VV(myi, :) = VV(myi, :)*Singular(ii)
                     endif
                  enddo

                  ! compute butterfly U and V
                  blocks%rankmax = rank
                  blocks%rankmin = rank

                  myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyV%blocks(1)%matrix(myArows, myAcols))
                  blocks%ButterflyV%blocks(1)%matrix = 0

                  call descinit(descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descButterflyV')

                  myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyU%blocks(1)%matrix(myArows, myAcols))
                  call descinit(descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descButterflyU')
                  blocks%ButterflyU%blocks(1)%matrix = UU(1:myArows, 1:myAcols)

                  call pgemmf90('N', 'T', N1, rank, rank1, cone, matV1, 1, 1, descsmatV1, VV, 1, 1, descVV, czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descButterflyV, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
                  call pgemmf90('N', 'T', N2, rank, rank2, cone, matV2, 1, 1, descsmatV2, VV, 1, 1 + rank1, descVV, czero, blocks%ButterflyV%blocks(1)%matrix, 1 + N1, 1, descButterflyV, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
                  deallocate (UU, VV, Singular, matV1, matV2, matU)
               else
                  blocks%rankmax = 0
                  blocks%rankmin = 0
                  rank = 0
                  deallocate (matV1, matV2, matU)
               endif
               call BF_delete(blocks%sons(1, 1), 1)
               call BF_delete(blocks%sons(2, 1), 1)
               deallocate (blocks%sons)

               call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmin, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

            else
               if (nprow /= -1 .and. npcol /= -1) then
                  call assert(N1 == N2, 'N1/=N2 in row merge')
                  myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (matV(myArows, myAcols))
                  call descinit(descsmatV, N1, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV')

                  myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
                  allocate (matU1(myArows, myAcols))
                  call descinit(descsmatU1, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU1')

                  myArows = numroc_wp(M2, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
                  allocate (matU2(myArows, myAcols))
                  call descinit(descsmatU2, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU2')
               else
                  descsmatV(2) = -1
                  descsmatU1(2) = -1
                  descsmatU2(2) = -1
                  allocate (matV(1, 1))
                  allocate (matU1(1, 1))
                  allocate (matU2(1, 1))
               endif

               if (nprow1 /= -1 .and. npcol1 /= -1) then
                  ! redistribute U1
                  myArows = numroc_wp(M1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit(descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU1c')
                  ! write(*,*)shape(blocks%sons(1,1)%ButterflyU%blocks(1)%matrix),shape(matU1),rank1,M1,blocks%sons(1,1)%rankmax
                  call pgemr2df90(M1, rank1, blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V1
                  myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit(descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV1c')
                  call pgemr2df90(N1, rank1, blocks%sons(1, 1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyV%blocks(1)%matrix)
               else
                  descsmatU1c(2) = -1
                  call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, ptree%pgrp(pgno)%ctxt1DCol)
                  descsmatV1c(2) = -1
                  call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, ptree%pgrp(pgno)%ctxt1DCol)
               endif
               if (allocated(blocks%sons(1, 1)%ButterflyU%blocks)) deallocate (blocks%sons(1, 1)%ButterflyU%blocks)
               if (allocated(blocks%sons(1, 1)%ButterflyV%blocks)) deallocate (blocks%sons(1, 1)%ButterflyV%blocks)

               if (nprow2 /= -1 .and. npcol2 /= -1) then
                  ! redistribute U2
                  myArows = numroc_wp(M2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit(descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatU2c')
                  call pgemr2df90(M2, rank2, blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V2
                  myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit(descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descsmatV2c')
                  call pgemr2df90(N2, rank2, blocks%sons(2, 1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV, 1, 1 + rank1, descsmatV, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyV%blocks(1)%matrix)
               else
                  descsmatU2c(2) = -1
                  call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, ptree%pgrp(pgno)%ctxt1DCol)
                  descsmatV2c(2) = -1
                  call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV, 1, 1 + rank1, descsmatV, ptree%pgrp(pgno)%ctxt1DCol)
               endif

               if (allocated(blocks%sons(2, 1)%ButterflyU%blocks)) deallocate (blocks%sons(2, 1)%ButterflyU%blocks)
               if (allocated(blocks%sons(2, 1)%ButterflyV%blocks)) deallocate (blocks%sons(2, 1)%ButterflyV%blocks)
               if (nprow /= -1 .and. npcol /= -1) then
                  ! compute truncated SVD on matV
                  mnmax = max(N1, rank1 + rank2)
                  mnmin = min(N1, rank1 + rank2)

                  myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
                  allocate (UU(myArows, myAcols))
                  call descinit(descUU, N1, mnmin, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descUU')
                  UU = 0
                  myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (VV(myArows, myAcols))
                  call descinit(descVV, mnmin, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descVV')
                  VV = 0

                  allocate (Singular(mnmin))
                  Singular = 0

                  call PSVD_Truncate(N1, rank1 + rank2, matV, descsmatV, UU, VV, descUU, descVV, Singular, option%tol_comp, rank, ptree%pgrp(pgno)%ctxt, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
                  do ii = 1, rank
                     call g2l(ii, rank, nprow, nbslpk, iproc, myi)
                     if (iproc == myrow) then
                        VV(myi, :) = VV(myi, :)*Singular(ii)
                     endif
                  enddo

                  ! compute butterfly U and V
                  blocks%rankmax = rank
                  blocks%rankmin = rank

                  myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyV%blocks(1)%matrix(myArows, myAcols))
                  blocks%ButterflyV%blocks(1)%matrix = UU(1:myArows, 1:myAcols)

                  call descinit(descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descButterflyV')

                  myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyU%blocks(1)%matrix(myArows, myAcols))
                  call descinit(descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit fail for descButterflyU')
                  blocks%ButterflyU%blocks(1)%matrix = 0

                  call pgemmf90('N', 'T', M1, rank, rank1, cone, matU1, 1, 1, descsmatU1, VV, 1, 1, descVV, czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descButterflyU, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

                  call pgemmf90('N', 'T', M2, rank, rank2, cone, matU2, 1, 1, descsmatU2, VV, 1, 1 + rank1, descVV, czero, blocks%ButterflyU%blocks(1)%matrix, 1 + M1, 1, descButterflyU, flop=flop)
                  stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
                  deallocate (UU, VV, Singular, matU1, matU2, matV)
               else
                  rank = 0
                  blocks%rankmax = 0
                  blocks%rankmin = 0
                  deallocate (matU1, matU2, matV)
               endif
               call BF_delete(blocks%sons(1, 1), 1)
               call BF_delete(blocks%sons(2, 1), 1)
               deallocate (blocks%sons)
               call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmin, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
               call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
            endif
         endif

         ! write(*,*)ptree%MyID,cridx,rank,'hei',blocks%M,blocks%N

101      if (cridx == 0) then

            !! the following is needed when there is idle procs in the process grids
            ranktmp = rank
            call MPI_ALLREDUCE(ranktmp, rank, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
            blocks%rankmax = rank
            blocks%rankmin = rank

            ! distribute UV factor into 1D grid
            ! write(*,*)rank,'wocanide',ptree%MyID
            call blacs_gridinfo(ptree%pgrp(pgno)%ctxt, nprow, npcol, myrow, mycol)
            if (myrow /= -1 .and. mycol /= -1) then
               myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
               ! if(myArows>0 .and. myAcols>0)then
               allocate (matU2D(myArows, myAcols))
               matU2D = blocks%ButterflyU%blocks(1)%matrix
               ! endif
               myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
               ! if(myArows>0 .and. myAcols>0)then
               allocate (matV2D(myArows, myAcols))
               matV2D = blocks%ButterflyV%blocks(1)%matrix
               ! endif
            else
               allocate (matU2D(1, 1))  ! required for Redistribute2Dto1D
               matU2D = 0
               allocate (matV2D(1, 1))  ! required for Redistribute2Dto1D
               matV2D = 0
            endif

            if (allocated(blocks%ButterflyU%blocks(1)%matrix)) deallocate (blocks%ButterflyU%blocks(1)%matrix)
            allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc, rank))
            blocks%ButterflyU%blocks(1)%matrix = 0

            if (allocated(blocks%ButterflyV%blocks(1)%matrix)) deallocate (blocks%ButterflyV%blocks(1)%matrix)
            allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc, rank))
            blocks%ButterflyV%blocks(1)%matrix = 0
! write(*,*) blocks%row_group,blocks%col_group,isnanMat(matU2D,size(matU2D,1),size(matU2D,2)),isnanMat(matV2D,size(matV2D,1),size(matV2D,2)),'nima222bi',ptree%MyID
            ! write(*,*)size(matU2D,1),size(matU2D,2),ptree%MyID,'dddddddd',myrow,mycol,myArows,myAcols
            call Redistribute2Dto1D(matU2D, blocks%M, 0, pgno, blocks%ButterflyU%blocks(1)%matrix, blocks%M_p, 0, pgno, rank, ptree)
            call Redistribute2Dto1D(matV2D, blocks%N, 0, pgno, blocks%ButterflyV%blocks(1)%matrix, blocks%N_p, 0, pgno, rank, ptree)

            ! write(*,*) fnorm(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),fnorm(blocks%ButterflyV%blocks(1)%matrix,blocks%N_loc,rank),isnanMat(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),'nimabi'
            ! write(*,*) blocks%row_group,blocks%col_group,isnanMat(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),'nimabi',ptree%MyID

            if (allocated(matU2D)) deallocate (matU2D)
            if (allocated(matV2D)) deallocate (matV2D)

         endif

      endif

   end subroutine LR_HMerge

   recursive subroutine LR_HBACA_Leaflevel(blocks, leafsize, rank, option, msh, ker, stats, ptree, pgno, cridx)
      use BPACK_DEFS
      implicit none
      integer rank, ranktmp, leafsize
      integer header_m, header_n
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, rmax, mn
      type(mesh)::msh
      type(Hoption)::option
      type(kernelquant)::ker
      type(matrixblock)::blocks, blockc(2)
      type(proctree)::ptree
      integer pgno, pgno1, pgno2
      integer:: cridx, info
      DT, allocatable:: UU(:, :), VV(:, :), matU(:, :), matV(:, :), matU1(:, :), matV1(:, :), matU2(:, :), matV2(:, :), tmp(:, :), matU1D(:, :), matV1D(:, :), Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Fullmat(:, :), QQ1(:, :), matU2D(:, :), matV2D(:, :)
      real(kind=8), allocatable::Singular(:)
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, M1, N1, M2, N2, rank1, rank2, ierr
      integer::descsmatU(9), descsmatV(9), descsmatU1(9), descsmatV1(9), descsmatU2(9), descsmatV2(9), descUU(9), descVV(9), descsmatU1c(9), descsmatU2c(9), descsmatV1c(9), descsmatV2c(9), descButterflyV(9), descButterflyU(9), descButterU1D(9), descButterV1D(9), descVin(9), descVout(9), descVinter(9), descFull(9)
      integer dims(6), dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
      DT:: TEMP(1)
      integer LWORK, mnmax, mnmin, rank_new
      DT, allocatable:: WORK(:)
      real(kind=8), allocatable::RWORK(:), center(:)
      real(kind=8):: rtemp, dist, error, rtemp1, rtemp0, fnorm1, fnorm0, flop
      integer :: nb1Dc, nb1Dr, frow, Dimn, edge_n, edge_m, MyID, Ntest, rmaxc, rmaxr
      integer, allocatable::M_p(:, :), N_p(:, :), mrange(:), nrange(:)
      type(Hstat)::stats
      integer::passflag = 0
      integer::Maxgrp

      Maxgrp = 2**(ptree%nlevel) - 1

      if (option%RecLR_leaf == ACANMERGE) then
         frow = 1
         rmax = min(option%rmax, min(blocks%M, blocks%N))
         call LR_ACA_Parallel(blocks, blocks%headm, blocks%headn, blocks%M, blocks%N, frow, rmax, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, error, ptree%pgrp(pgno)%ctxt, pgno)
         ! call LR_SeudoSkeleton(blocks,blocks%headm,blocks%headn,blocks%M,blocks%N,min(blocks%N,1000),min(blocks%M,1000),rank,option%tol_comp,option%tol_comp,msh,ker,stats,ptree,option,gd%ctxt)

      else

         rank = 0
         blocks%ButterflyU%idx = 1
         blocks%ButterflyU%inc = 1
         blocks%ButterflyU%nblk_loc = 1
         blocks%ButterflyV%idx = 1
         blocks%ButterflyV%inc = 1
         blocks%ButterflyV%nblk_loc = 1

         if (mod(cridx, 2) == 0 .and. (min(blocks%M, blocks%N) <= leafsize .and. (pgno*2 > Maxgrp))) then ! reach bottom level, call sequential aca,rrqr etc.

            if (option%RecLR_leaf == PS) then
               ! !!!!! CUR

               ! rmaxc = min(blocks%M,blocks%N)
               ! rmaxr = min(blocks%M,blocks%N)
               rmaxc = blocks%N
               rmaxr = blocks%M

               call LR_SeudoSkeleton(blocks, blocks%headm, blocks%headn, blocks%M, blocks%N, rmaxc, rmaxr, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, ptree%pgrp(pgno)%ctxt)

            else if (option%RecLR_leaf == ACA) then
               !!!!! ACA-SVD
               rmax = min(option%rmax, min(blocks%M, blocks%N))
               allocate (UU(blocks%M, rmax))
               allocate (VV(rmax, blocks%N))
               allocate (Singular(rmax))
               frow = 1
               ! !!!!!!!!!!!! picking a good first row may requires some geometry information. The following can be commented out if geometry info is not available
               ! Dimn = 0
               ! if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
               ! dist = 1D300
               ! allocate(center(Dimn))
               ! center = 0
               ! header_n = blocks%headn
               ! do j=1,blocks%N
               ! edge_n = header_n-1+j
               ! center = center + msh%xyz(1:Dimn,msh%new2old(edge_n))
               ! enddo
               ! center = center/blocks%N

               ! header_m = blocks%headm
               ! do i=1,blocks%M
               ! edge_m=header_m-1+i
               ! rtemp = sum((msh%xyz(1:Dimn,msh%new2old(edge_m))-center(1:Dimn))**2d0)
               ! if(rtemp<dist)then
               ! dist = rtemp
               ! frow = i
               ! endif
               ! enddo
               ! deallocate(center)
               ! !!!!!!!!!!!!

               call LR_ACA(UU, VV, Singular, blocks%headm, blocks%headn, blocks%M, blocks%N, frow, rmax, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, error)

               ! if(error>option%tol_comp)then
               ! write(*,*)'niam',error
               ! deallocate (UU,VV,Singular)
               ! goto 100
               ! endif

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = VV(j, i)
                  enddo
               enddo
               !$omp end parallel do

               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = UU(i, j)*Singular(j)
                  enddo
               enddo
               !$omp end parallel do
               deallocate (UU, VV, Singular)

            else if (option%RecLR_leaf == BACA .or. option%RecLR_leaf == BACANOVER) then
               !!! blocked ACA

               rmax = min(option%rmax, min(blocks%M, blocks%N))
               allocate (UU(blocks%M, rmax))
               allocate (VV(rmax, blocks%N))

               if (option%RecLR_leaf == BACA) call LR_BACA(UU, VV, blocks%headm, blocks%headn, blocks%M, blocks%N, rmax, rank, option%tol_comp, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)
               if (option%RecLR_leaf == BACANOVER) call LR_BACA_noOverlap(UU, VV, blocks%headm, blocks%headn, blocks%M, blocks%N, rmax, rank, option%tol_comp, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = VV(j, i)
                  enddo
               enddo
               !$omp end parallel do

               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = UU(i, j)
                  enddo
               enddo
               !$omp end parallel do
               deallocate (UU, VV)

            else if (option%RecLR_leaf == SVD) then
               !!!!! SVD
               mn = min(blocks%M, blocks%N)
               allocate (UU(blocks%M, mn), VV(mn, blocks%N), Singular(mn))
               allocate (QQ1(blocks%M, blocks%N))
               allocate (mrange(blocks%M))
               allocate (nrange(blocks%N))
               do ii = 1, blocks%M
                  mrange(ii) = blocks%headm + ii - 1
               enddo
               do jj = 1, blocks%N
                  nrange(jj) = blocks%headn + jj - 1
               enddo
               call element_Zmn_block_user(blocks%M, blocks%N, mrange, nrange, QQ1, msh, option, ker, 0, passflag, ptree, stats)
               deallocate (mrange)
               deallocate (nrange)

               ! do ii=1,blocks%M
               ! do jj =1,blocks%N
               ! edge_m = blocks%headm +ii - 1
               ! edge_n = blocks%headn +jj - 1
               ! call element_Zmn(edge_m,edge_n,QQ1(ii,jj),msh,option,ker)
               ! end do
               ! end do

               call SVD_Truncate(QQ1, blocks%M, blocks%N, mn, UU, VV, Singular, option%tol_comp, rank, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop

               ! rank=blocks%N

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))

               ! blocks%ButterflyV%blocks(1)%matrix=0
               ! do j=1, rank
               ! blocks%ButterflyV%blocks(1)%matrix(j,j)=1d0
               ! enddo

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = VV(j, i)
                  enddo
               enddo
               !$omp end parallel do

               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))
               ! blocks%ButterflyU%blocks(1)%matrix = QQ1

               !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = UU(i, j)*Singular(j)
                  enddo
               enddo
               !$omp end parallel do

               deallocate (QQ1, UU, VV, Singular)
            else if (option%RecLR_leaf == RRQR) then
               !!!!! RRQR
               mn = min(blocks%M, blocks%N)
               allocate (UU(blocks%M, mn), VV(mn, blocks%N))

               allocate (QQ1(blocks%M, blocks%N))
               allocate (mrange(blocks%M))
               allocate (nrange(blocks%N))
               do ii = 1, blocks%M
                  mrange(ii) = blocks%headm + ii - 1
               enddo
               do jj = 1, blocks%N
                  nrange(jj) = blocks%headn + jj - 1
               enddo
               call element_Zmn_block_user(blocks%M, blocks%N, mrange, nrange, QQ1, msh, option, ker, 0, passflag, ptree, stats)

               ! write(*,*)'ni',blocks%row_group,blocks%col_group,fnorm(QQ1,blocks%M,blocks%N)

               deallocate (mrange)
               deallocate (nrange)

               ! allocate(QQ1(blocks%M,blocks%N))
               ! do ii=1,blocks%M
               ! do jj =1,blocks%N
               ! edge_m = blocks%headm +ii - 1
               ! edge_n = blocks%headn +jj - 1
               ! call element_Zmn(edge_m,edge_n,QQ1(ii,jj),msh,option,ker)
               ! end do
               ! end do

               call RRQR_LQ(QQ1, blocks%M, blocks%N, mn, UU, VV, option%tol_comp, rank, 'L', flops=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop

               ! rank=blocks%N

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))

               ! blocks%ButterflyV%blocks(1)%matrix=0
               ! do j=1, rank
               ! blocks%ButterflyV%blocks(1)%matrix(j,j)=1d0
               ! enddo

               ! !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = VV(j, i)
                  enddo
               enddo
               ! !$omp end parallel do

               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))
               ! blocks%ButterflyU%blocks(1)%matrix = QQ1

               ! !$omp parallel do default(shared) private(i,j)
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = UU(i, j)
                  enddo
               enddo
               ! !$omp end parallel do

               deallocate (QQ1, UU, VV)

            endif

         else

            if (pgno*2 > Maxgrp) then
               pgno1 = pgno
               pgno2 = pgno
            else
               pgno1 = pgno*2
               pgno2 = pgno*2 + 1
            endif

            if (IOwnPgrp(ptree, pgno)) then
               nsproc1 = ptree%pgrp(pgno1)%nproc
               nsproc2 = ptree%pgrp(pgno2)%nproc

               allocate (blocks%sons(2, 1))
               if (IOwnPgrp(ptree, pgno1)) then

                  ! proportional mapping along row or column dimensions
                  if (mod(cridx + 1, 2) == 0) then  ! split along column dimension
                     blocks%sons(1, 1)%headm = blocks%headm
                     blocks%sons(1, 1)%M = blocks%M
                     blocks%sons(1, 1)%headn = blocks%headn
                     blocks%sons(1, 1)%N = INT(blocks%N*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     call assert(blocks%sons(1, 1)%N > 0 .and. blocks%sons(1, 1)%N < blocks%N, 'column of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
                  else  ! split along row dimension
                     blocks%sons(1, 1)%headn = blocks%headn
                     blocks%sons(1, 1)%N = blocks%N
                     blocks%sons(1, 1)%headm = blocks%headm
                     blocks%sons(1, 1)%M = INT(blocks%M*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     call assert(blocks%sons(1, 1)%M > 0 .and. blocks%sons(1, 1)%M < blocks%M, 'row of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
                  endif
                  ! write(*,*)blocks%sons(1,1)%M,blocks%sons(1,1)%N,'ha1',nsproc1,nsproc2
                  allocate (blocks%sons(1, 1)%ButterflyU%blocks(1))
                  allocate (blocks%sons(1, 1)%ButterflyV%blocks(1))

                  call LR_HBACA_Leaflevel(blocks%sons(1, 1), leafsize, rank, option, msh, ker, stats, ptree, pgno1, cridx + 1)
               endif

               if (IOwnPgrp(ptree, pgno2)) then

                  ! proportional mapping along row or column dimensions
                  if (mod(cridx + 1, 2) == 0) then  ! split along column dimension
                     blocks%sons(2, 1)%headm = blocks%headm
                     blocks%sons(2, 1)%M = blocks%M
                     blocks%sons(2, 1)%headn = blocks%headn + INT(blocks%N*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     blocks%sons(2, 1)%N = blocks%N - INT(blocks%N*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     call assert(blocks%sons(2, 1)%N > 0 .and. blocks%sons(2, 1)%N < blocks%N, 'column of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
                  else  ! split along row dimension
                     blocks%sons(2, 1)%headn = blocks%headn
                     blocks%sons(2, 1)%N = blocks%N
                     blocks%sons(2, 1)%headm = blocks%headm + INT(blocks%M*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     blocks%sons(2, 1)%M = blocks%M - INT(blocks%M*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
                     call assert(blocks%sons(2, 1)%M > 0 .and. blocks%sons(2, 1)%M < blocks%M, 'row of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
                  endif
                  allocate (blocks%sons(2, 1)%ButterflyU%blocks(1))
                  allocate (blocks%sons(2, 1)%ButterflyV%blocks(1))
                  call LR_HBACA_Leaflevel(blocks%sons(2, 1), leafsize, rank, option, msh, ker, stats, ptree, pgno2, cridx + 1)
               endif
            endif
         endif
      endif

   end subroutine LR_HBACA_Leaflevel

   subroutine LR_ACA(matU, matV, Singular, header_m, header_n, rankmax_r, rankmax_c, frow, rmax, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, error)

      use BPACK_DEFS
      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance, dist
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n, Dimn, mn, Navr, itr
      integer rank, ranknew, row, column, rankmax, rankmax_c, rankmax_r, rankmax_min, rmax, idxs_r, idxs_c, frow
      DT value_Z, maxvalue
      DT inner_U, inner_V, ctemp, value_UVs
      real(kind=8) inner_UV, n1, n2, a, error, flop
      integer, allocatable:: select_column(:), select_row(:)
      DT::matU(rankmax_r, rmax), matV(rmax, rankmax_c)
      DT::matr(1, rankmax_c), matc(rankmax_r, 1)
      real(kind=8)::Singular(rmax)
      DT, allocatable:: row_R(:), column_R(:), value_UV(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:), norm_UVavrbynorm_Z(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      real(kind=8), allocatable :: Singularsml(:)
      integer::mrange(rankmax_r), nrange(rankmax_c)
      type(Hstat)::stats
      type(Hoption)::option
      integer:: passflag = 0

      Navr = 3 !5 !10
      itr = 1
      allocate (norm_UVavrbynorm_Z(Navr))
      norm_UVavrbynorm_Z = 0

      n1 = OMP_get_wtime()

      allocate (select_column(rankmax_c))
      allocate (select_row(rankmax_r))
      allocate (value_UV(max(rankmax_c, rankmax_r)))
      value_UV = 0

      rankmax_min = min(rankmax_r, rankmax_c)
      norm_Z = 0
      select_column = 0
      select_row = 0

      allocate (row_R(rankmax_c), column_R(rankmax_r))
      allocate (norm_row_R(rankmax_c), norm_column_R(rankmax_r))
      row_R = 0
      column_R = 0
      norm_row_R = 0
      norm_column_R = 0

      select_row(1) = frow

      mrange = select_row(1)
      do j = 1, rankmax_c
         nrange(j) = header_n + j - 1
      enddo
      call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)
      row_R = matr(1, :)
      norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))

      ! !$omp parallel do default(shared) private(j,value_Z,edge_m,edge_n)
      ! do j=1, rankmax_c
      ! ! value_Z=mat(select_row(1),j)
      ! edge_m = header_m + select_row(1) - 1
      ! edge_n = header_n + j - 1
      ! call element_Zmn(edge_m,edge_n,value_Z,msh,option,ker)
      ! row_R(j)=value_Z
      ! norm_row_R(j)=dble(value_Z*conjg(cmplx(value_Z,kind = 8)))
      ! enddo
      ! !$omp end parallel do

      select_column(1) = maxloc(norm_row_R, 1)
      maxvalue = row_R(select_column(1))

      if (abs(maxvalue) < SafeUnderflow) then

         do ii = 1, 100
            a = 0
            call random_number(a)
            select_row(1) = floor_safe(a*(rankmax_r - 1)) + 1

            mrange = select_row(1)
            do j = 1, rankmax_c
               nrange(j) = header_n + j - 1
            enddo
            call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)
            row_R = matr(1, :)
            norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))

            ! !$omp parallel do default(shared) private(j,value_Z,edge_m,edge_n)
            ! do j=1, rankmax_c
            ! ! value_Z=mat(select_row(1),j)
            ! edge_m = header_m + select_row(1) - 1
            ! edge_n = header_n + j - 1
            ! call element_Zmn(edge_m,edge_n,value_Z,msh,option,ker)
            ! row_R(j)=value_Z
            ! norm_row_R(j)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
            ! enddo
            ! !$omp end parallel do

            select_column(1) = maxloc(norm_row_R, 1)
            maxvalue = row_R(select_column(1))
            if (abs(maxvalue) > SafeUnderflow) exit
         end do
         if (abs(maxvalue) < SafeUnderflow) then
            rank = 1
            matU(:, 1) = 0
            matV(1, :) = 0
            Singular(1) = 0

            deallocate (select_column)
            deallocate (select_row)
            deallocate (value_UV)
            deallocate (row_R, column_R)
            deallocate (norm_row_R, norm_column_R)
            deallocate (norm_UVavrbynorm_Z)
            return
         endif
      end if

      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         row_R(j) = row_R(j)/maxvalue
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         matV(1, j) = row_R(j)
      enddo
      ! !$omp end parallel do

      nrange = header_n + select_column(1) - 1
      do i = 1, rankmax_r
         mrange(i) = header_m + i - 1
      enddo
      call element_Zmn_block_user(rankmax_r, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)
      column_R = matc(:, 1)
      norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))

      ! !$omp parallel do default(shared) private(i,value_Z,edge_m,edge_n)
      ! do i=1,rankmax_r
      ! edge_m = header_m + i - 1
      ! edge_n = header_n + select_column(1) - 1
      ! call element_Zmn(edge_m,edge_n,value_Z,msh,option,ker)
      ! ! value_Z=mat(i,select_column(1))
      ! column_R(i)=value_Z
      ! norm_column_R(i)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
      ! enddo
      ! !$omp end parallel do

      norm_column_R(select_row(1)) = 0

      ! !$omp parallel do default(shared) private(i)
      do i = 1, rankmax_r
         matU(i, 1) = column_R(i)
      enddo
      ! !$omp end parallel do

      norm_U = norm_vector(column_R, rankmax_r)
      norm_V = norm_vector(row_R, rankmax_c)
      norm_Z = norm_Z + norm_U*norm_V
      if (norm_Z > SafeUnderflow) then
         norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
      else
         norm_UVavrbynorm_Z(itr) = 0
      endif
      itr = mod(itr, Navr) + 1

      ! if(rankmax<2)write(*,*)'rankmax'
      select_row(2) = maxloc(norm_column_R, 1)

      rank = 1
      ! write(*,*)column_R,row_R
      ! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
      do while (tolerance**2 < (sum(norm_UVavrbynorm_Z)/min(Navr, rank)) .and. rank < rankmax_min)

         mrange(1) = header_m + select_row(rank + 1) - 1
         do j = 1, rankmax_c
            nrange(j) = header_n + j - 1
         enddo
         call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)
         row_R = matr(1, :)

         ! !$omp parallel do default(shared) private(j,i,value_Z,edge_m,edge_n)
         ! do j=1,rankmax_c
         ! edge_m = header_m + select_row(rank+1) - 1
         ! edge_n = header_n + j - 1
         ! call element_Zmn(edge_m,edge_n,row_R(j),msh,option,ker)
         ! enddo
         ! !$omp end parallel do

         call gemmf77('N', 'N', 1, rankmax_c, rank, cone, matU(select_row(rank + 1), 1), rankmax_r, matV, rmax, czero, value_UV, 1)

         ! call gemmf90(matU(select_row(rank+1),1), rankmax_r,matV,rmax,value_UV,1,'N','N',1,rankmax_c,rank,cone,czero)

         row_R = row_R - value_UV(1:rankmax_c)
         norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         if (abs(maxvalue) < SafeUnderflow) then
            ! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
            exit
            row_R = 0
         else
            do j = 1, rankmax_c
               row_R(j) = row_R(j)/maxvalue
            enddo
         endif
         ! !$omp parallel do default(shared) private(j)

         ! !$omp end parallel do
         ! !$omp parallel do default(shared) private(j)
         do j = 1, rankmax_c
            matV(rank + 1, j) = row_R(j)
         enddo
         ! !$omp end parallel do

         nrange(1) = header_n + select_column(rank + 1) - 1
         do i = 1, rankmax_r
            mrange(i) = header_m + i - 1
         enddo
         call element_Zmn_block_user(rankmax_r, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)
         column_R = matc(:, 1)

         ! !$omp parallel do default(shared) private(i,j,value_Z,value_UVs,edge_m,edge_n)
         ! do i=1,rankmax_r
         ! edge_m = header_m + i - 1
         ! edge_n = header_n + select_column(rank+1) - 1
         ! call element_Zmn(edge_m,edge_n,column_R(i),msh,option,ker)
         ! enddo
         ! !$omp end parallel do

         call gemmf77('N', 'N', rankmax_r, 1, rank, cone, matU, rankmax_r, matV(1, select_column(rank + 1)), rmax, czero, value_UV, rankmax_r)

         ! call gemmf90(matU, rankmax_r,matV(1,select_column(rank+1)),rmax,value_UV,rankmax_r,'N','N',rankmax_r,1,rank,cone,czero)

         column_R = column_R - value_UV(1:rankmax_r)
         norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_r

         do i = 1, rank + 1
            norm_column_R(select_row(i)) = 0
         enddo

         ! !$omp parallel do default(shared) private(i)
         do i = 1, rankmax_r
            matU(i, rank + 1) = column_R(i)
         enddo
         ! !$omp end parallel do

         norm_U = norm_vector(column_R, rankmax_r)
         norm_V = norm_vector(row_R, rankmax_c)

         inner_UV = 0
         !$omp parallel do default(shared) private(j,i,ctemp,inner_V,inner_U) reduction(+:inner_UV)
         do j = 1, rank
            inner_U = 0
            inner_V = 0
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i = 1, rankmax_r
               ctemp = matU(i, rank + 1)*conjg(cmplx(matU(i, j), kind=8))
               inner_U = inner_U + ctemp
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i = 1, rankmax_c
               ctemp = matV(rank + 1, i)*conjg(cmplx(matV(j, i), kind=8))
               inner_V = inner_V + ctemp
            enddo
            ! !$omp end parallel do
            inner_UV = inner_UV + 2*dble(inner_U*inner_V)
         enddo
         !$omp end parallel do

         norm_Z = norm_Z + inner_UV + norm_U*norm_V

         ! ! write(*,*)norm_Z,inner_UV,norm_U,norm_V,maxvalue,rank,'gan'
         ! if(isnan(sqrt(norm_Z)))then
         ! write(*,*)inner_UV,norm_U,norm_V,maxvalue
         ! stop
         ! endif

         if (norm_Z > SafeUnderflow) then
            norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
         else
            norm_UVavrbynorm_Z(itr) = 0
         endif
         itr = mod(itr, Navr) + 1

         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c + rank*rankmax_r

         rank = rank + 1
         if (rank > rmax) then
            ! write(*,*)'increase rmax',rank,rmax
            ! stop
            exit
         end if
         if (rank < rankmax_min) then
            select_row(rank + 1) = maxloc(norm_column_R, 1)
         endif

         if (norm_Z < 0) exit

         ! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)

      enddo
      ! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)

      error = sqrt((sum(norm_UVavrbynorm_Z)/Navr))

      ! write(*,*)select_row(1:rank),select_column(1:rank)

      deallocate (row_R, column_R)
      deallocate (norm_row_R, norm_column_R)
      deallocate (norm_UVavrbynorm_Z)

      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

! ACA followed by SVD

      allocate (QQ1(rankmax_r, rank))
      QQ1 = matU(1:rankmax_r, 1:rank)
      ! call copymatN(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
      allocate (tau_Q(rank))
      call geqrff90(QQ1, tau_Q, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop

      allocate (RR1(rank, rank))
      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, j
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ1, tau_Q, rankmax_r, rank, rank, flop=flop)
      deallocate (tau_Q)
      stats%Flop_Fill = stats%Flop_Fill + flop

      allocate (QQ2(rankmax_c, rank))
      call copymatT(matV(1:rank, 1:rankmax_c), QQ2, rank, rankmax_c)
      allocate (tau_Q(rank))
      call geqrff90(QQ2, tau_Q, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop

      allocate (RR2(rank, rank))
      RR2 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, j
            RR2(i, j) = QQ2(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ2, tau_Q, rankmax_c, rank, rank, flop=flop)
      deallocate (tau_Q)
      stats%Flop_Fill = stats%Flop_Fill + flop

      allocate (mattemp(rank, rank))
      mattemp = 0
      call gemmf90(RR1, rank, RR2, rank, mattemp, rank, 'N', 'T', rank, rank, rank, cone, czero, flop=flop)
      ! call zgemm('N','T',rank,rank,rank, cone, RR1, rank,RR2,rank,czero,mattemp,rank)
      stats%Flop_Fill = stats%Flop_Fill + flop
      allocate (UUsml(rank, rank), VVsml(rank, rank), Singularsml(rank))
      call SVD_Truncate(mattemp, rank, rank, rank, UUsml, VVsml, Singularsml, SVD_tolerance, ranknew, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop
      ! call zgemm('N','N',rankmax_r,ranknew,rank, cone, QQ1, rankmax_r,UUsml,rank,czero,matU,rankmax_r)
      call gemmf90(QQ1, rankmax_r, UUsml, rank, matU, rankmax_r, 'N', 'N', rankmax_r, ranknew, rank, cone, czero, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop
      ! call zgemm('N','T',ranknew,rankmax_c,rank, cone, VVsml, rank,QQ2,rankmax_c,czero,matV,rmax)
      call gemmf90(VVsml, rank, QQ2, rankmax_c, matV, rmax, 'N', 'T', ranknew, rankmax_c, rank, cone, czero, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop

      rank = ranknew
      Singular(1:ranknew) = Singularsml(1:ranknew)

      deallocate (mattemp, RR1, QQ1, UUsml, VVsml, Singularsml)
      deallocate (QQ2, RR2)

      deallocate (select_column)
      deallocate (select_row)
      deallocate (value_UV)

      return

   end subroutine LR_ACA

   subroutine LR_ACA_Parallel(blocks, header_m, header_n, M, N, frow, rmax, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, error, ctxt, pgno)

      use BPACK_DEFS
      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      type(matrixblock)::blocks
      integer index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance, dist
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n, Dimn, mn, Navr, itr
      integer rank, rank1, rank2, ranknew, row, column, rankmax, N, M, rankmax_min, rmax, idxs_r, idxs_c, frow, mn1, mn2
      DT value_Z, maxvalue
      DT inner_U, inner_V, ctemp, value_UVs
      real(kind=8) inner_UV, n1, n2, a, error, flop
      integer:: select_column(N), select_row(M)
      DT, allocatable::matU(:, :), matV(:, :), matU2D(:, :), matV2D(:, :)
      DT::tmpval(rmax)
      DT::matr(1, N), matc(M, 1)
      real(kind=8)::Singular(rmax)
      DT, allocatable:: value_UV(:)
      DT:: row_R(N), column_R(M)
      real(kind=8):: norm_row_R(N), norm_column_R(M)
      real(kind=8), allocatable:: norm_UVavrbynorm_Z(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), UUu(:, :), UUv(:, :), VVu(:, :), VVv(:, :)
      real(kind=8), allocatable :: Singularsml(:), Singularuv(:)
      integer::mrange(M), nrange(N)
      type(Hstat)::stats
      type(Hoption)::option
      integer:: passflag = 0
      integer ctxt, iproc, jproc, myi, myj, mnmax
      integer::pgno
      integer::headm_loc, headn_loc, pp
      integer::ierr
      integer myArows, myAcols, info, nprow, npcol, myrow, mycol, taun
      integer::descsMatU2D(9), descsMatV2D(9), descsMatSml(9), descsMatU2Dnew(9), descsMatV2Dnew(9), descsUUSml(9), descsVVSml(9), descsUU_u(9), descsVV_u(9), descsUU_v(9), descsVV_v(9)

      pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
      headm_loc = blocks%M_p(pp, 1)
      headn_loc = blocks%N_p(pp, 1)

      allocate (matU(blocks%M_loc, rmax))
      allocate (matV(rmax, blocks%N_loc))

      Navr = 3 !5 !10
      itr = 1
      allocate (norm_UVavrbynorm_Z(Navr))
      norm_UVavrbynorm_Z = 0

      n1 = OMP_get_wtime()

      allocate (value_UV(max(N, M)))
      value_UV = 0

      rankmax_min = min(M, N)
      norm_Z = 0
      select_column = 0
      select_row = 0

      row_R = 0
      column_R = 0
      norm_row_R = 0
      norm_column_R = 0

      select_row(1) = frow

      mrange = select_row(1)
      do j = 1, blocks%N_loc
         nrange(j) = header_n + j - 1 + headn_loc - 1
      enddo
      call element_Zmn_block_user(1, blocks%N_loc, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)
      row_R = 0
      row_R(headn_loc:headn_loc + blocks%N_loc - 1) = matr(1, 1:blocks%N_loc)

      call MPI_ALLREDUCE(MPI_IN_PLACE, row_R, N, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

      norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))

      select_column(1) = maxloc(norm_row_R, 1)
      maxvalue = row_R(select_column(1))

      ! write(*,*)'wori',select_column(1),sum(norm_row_R)

      ! !$omp parallel do default(shared) private(j)
      do j = 1, N
         row_R(j) = row_R(j)/maxvalue
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(j)
      do j = 1, blocks%N_loc
         matV(1, j) = row_R(j + headn_loc - 1)
      enddo
      ! !$omp end parallel do

      nrange = header_n + select_column(1) - 1
      do i = 1, blocks%M_loc
         mrange(i) = header_m + i - 1 + headm_loc - 1
      enddo
      call element_Zmn_block_user(blocks%M_loc, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)
      column_R = 0
      column_R(headm_loc:headm_loc + blocks%M_loc - 1) = matc(1:blocks%M_loc, 1)
      call MPI_ALLREDUCE(MPI_IN_PLACE, column_R, M, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
      norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))

      norm_column_R(select_row(1)) = 0

      ! !$omp parallel do default(shared) private(i)
      do i = 1, blocks%M_loc
         matU(i, 1) = column_R(i + headm_loc - 1)
      enddo
      ! !$omp end parallel do

      norm_U = norm_vector(column_R, M)
      norm_V = norm_vector(row_R, N)
      norm_Z = norm_Z + norm_U*norm_V
      if (norm_Z > SafeUnderflow) then
         norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
      else
         norm_UVavrbynorm_Z(itr) = 0
      endif
      itr = mod(itr, Navr) + 1

      ! if(rankmax<2)write(*,*)'rankmax'
      select_row(2) = maxloc(norm_column_R, 1)

      ! write(*,*)'wori',select_row(2),sum(norm_column_R)

      rank = 1
      ! write(*,*)column_R,row_R
      ! write(*,*)norm_Z,norm_U*norm_V,'hehe'
      do while (tolerance**2 < (sum(norm_UVavrbynorm_Z)/min(Navr, rank)) .and. rank < rankmax_min)

         mrange(1) = header_m + select_row(rank + 1) - 1
         do j = 1, blocks%N_loc
            nrange(j) = header_n + j - 1 + headn_loc - 1
         enddo
         call element_Zmn_block_user(1, blocks%N_loc, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)
         row_R = 0
         row_R(headn_loc:headn_loc + blocks%N_loc - 1) = matr(1, 1:blocks%N_loc)
         call MPI_ALLREDUCE(MPI_IN_PLACE, row_R, N, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         value_UV = 0
         if (select_row(rank + 1) >= headm_loc .and. select_row(rank + 1) <= headm_loc + blocks%M_loc - 1) then
            tmpval(1:rank) = matU(select_row(rank + 1) - headm_loc + 1, 1:rank)
         else
            tmpval(1:rank) = 0
         endif
         call MPI_ALLREDUCE(MPI_IN_PLACE, tmpval, rank, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         call gemmf77('N', 'N', 1, blocks%N_loc, rank, cone, tmpval(1), 1, matV, rmax, czero, value_UV(headn_loc), 1)
         call MPI_ALLREDUCE(MPI_IN_PLACE, value_UV, N, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         ! call gemmf90(matU(select_row(rank+1),1), M,matV,rmax,value_UV,1,'N','N',1,N,rank,cone,czero)

         row_R = row_R - value_UV(1:N)
         norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*blocks%N_loc

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         if (abs(maxvalue) < SafeUnderflow) then
            ! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
            exit
            row_R = 0
         else
            do j = 1, N
               row_R(j) = row_R(j)/maxvalue
            enddo
         endif
         ! !$omp parallel do default(shared) private(j)

         ! !$omp end parallel do
         ! !$omp parallel do default(shared) private(j)
         do j = 1, blocks%N_loc
            matV(rank + 1, j) = row_R(j + headn_loc - 1)
         enddo
         ! !$omp end parallel do

         nrange(1) = header_n + select_column(rank + 1) - 1
         do i = 1, blocks%M_loc
            mrange(i) = header_m + i - 1 + headm_loc - 1
         enddo
         call element_Zmn_block_user(blocks%M_loc, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)
         column_R = 0
         column_R(headm_loc:headm_loc + blocks%M_loc - 1) = matc(1:blocks%M_loc, 1)
         call MPI_ALLREDUCE(MPI_IN_PLACE, column_R, M, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         value_UV = 0
         if (select_column(rank + 1) >= headn_loc .and. select_column(rank + 1) <= headn_loc + blocks%N_loc - 1) then
            tmpval(1:rank) = matV(1:rank, select_column(rank + 1) - headn_loc + 1)
         else
            tmpval(1:rank) = 0
         endif
         ! write(*,*)'dd',ptree%MyID,sum(value_UV(1:M)),sum(tmpval(1:rank))
         call MPI_ALLREDUCE(MPI_IN_PLACE, tmpval, rank, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
         ! write(*,*)'dd3',ptree%MyID,sum(value_UV(1:M)),sum(tmpval(1:rank))
         call gemmf77('N', 'N', blocks%M_loc, 1, rank, cone, matU, blocks%M_loc, tmpval(1), rank, czero, value_UV(headm_loc), M)

         ! write(*,*)'dd4',ptree%MyID,sum(value_UV(1:M))

         call MPI_ALLREDUCE(MPI_IN_PLACE, value_UV, M, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         ! call gemmf90(matU, M,matV(1,select_column(rank+1)),rmax,value_UV,M,'N','N',M,1,rank,cone,czero)

         column_R = column_R - value_UV(1:M)
         norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*blocks%M_loc

         do i = 1, rank + 1
            norm_column_R(select_row(i)) = 0
         enddo

         ! !$omp parallel do default(shared) private(i)
         do i = 1, blocks%M_loc
            matU(i, rank + 1) = column_R(i + headm_loc - 1)
         enddo
         ! !$omp end parallel do

         norm_U = norm_vector(column_R, M)
         norm_V = norm_vector(row_R, N)

         inner_UV = 0
         ! !$omp parallel do default(shared) private(j,i,ctemp,inner_V,inner_U) reduction(+:inner_UV)
         do j = 1, rank
            inner_U = 0
            inner_V = 0
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i = 1, blocks%M_loc
               ctemp = matU(i, rank + 1)*conjg(cmplx(matU(i, j), kind=8))
               inner_U = inner_U + ctemp
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i = 1, blocks%N_loc
               ctemp = matV(rank + 1, i)*conjg(cmplx(matV(j, i), kind=8))
               inner_V = inner_V + ctemp
            enddo
            ! !$omp end parallel do
            call MPI_ALLREDUCE(MPI_IN_PLACE, inner_U, 1, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE, inner_V, 1, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
            inner_UV = inner_UV + 2*dble(inner_U*inner_V)

         enddo
         ! !$omp end parallel do

         norm_Z = norm_Z + inner_UV + norm_U*norm_V

         if (norm_Z > SafeUnderflow) then
            norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
         else
            norm_UVavrbynorm_Z(itr) = 0
         endif
         itr = mod(itr, Navr) + 1

         stats%Flop_Fill = stats%Flop_Fill + rank*blocks%N_loc + rank*blocks%M_loc

         ! if(ptree%MyID==Main_ID)write(*,*)rank,norm_U*norm_V,norm_Z,select_column(rank+1),select_row(rank+1)

         rank = rank + 1
         if (rank > rmax) then
            ! write(*,*)'increase rmax',rank,rmax
            ! stop
            exit
         end if
         if (rank < rankmax_min) then
            select_row(rank + 1) = maxloc(norm_column_R, 1)
         endif

         if (norm_Z < 0) exit

      enddo

      error = sqrt((sum(norm_UVavrbynorm_Z)/Navr))

      deallocate (norm_UVavrbynorm_Z)

      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

! ! ACA followed by SVD

      !!!!**** generate 2D grid blacs quantities for matU
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit(descsMatU2D, M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matU2D(myArows, myAcols))
         matU2D = 0

      else
         descsMatU2D(2) = -1
         allocate (matU2D(1, 1))
         matU2D = 0
      endif
      !!!!**** redistribution of input matrix
      call Redistribute1Dto2D(matU, blocks%M_p, 0, pgno, matU2D, M, 0, pgno, rank, ptree)

      !!!!**** generate 2D grid blacs quantities for matV transpose
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit(descsMatV2D, N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matV2D(myArows, myAcols))
         matV2D = 0
      else
         descsMatV2D(2) = -1
         allocate (matV2D(1, 1))
         matV2D = 0
      endif
      !!!!**** redistribution of input matrix
      allocate (matV1(blocks%N_loc, rmax))
      call copymatT(matV, matV1, rmax, blocks%N_loc)
      call Redistribute1Dto2D(matV1, blocks%N_p, 0, pgno, matV2D, N, 0, pgno, rank, ptree)
      deallocate (matV1)

      if (myrow /= -1 .and. mycol /= -1) then
         ! myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         ! myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         ! allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
         ! blocks%ButterflyU%blocks(1)%matrix=matU2D
         ! myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         ! myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         ! allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
         ! blocks%ButterflyV%blocks(1)%matrix=matV2D

! GCC 9 was segfaulting when using PXGEQRF when it calls PXGEQR2 ->PXLARFG->PXSCAL
#if __GNUC__ < 9

         mn = min(M, rank)
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn, nbslpk, mycol, 0, npcol)
         allocate (tau_Q(myAcols))
         call pgeqrff90(M, mn, matU2D, 1, 1, descsMatU2D, tau_Q, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (RR1(myArows, myAcols))
         RR1 = 0d0
         do myj = 1, myAcols
            call l2g(myj, mycol, rank, npcol, nbslpk, jj)
            do myi = 1, myArows
               call l2g(myi, myrow, rank, nprow, nbslpk, ii)
               if (ii <= jj) RR1(myi, myj) = matU2D(myi, myj)
            enddo
         enddo
         call pun_or_gqrf90(ctxt, matU2D, tau_Q, M, rank, rank, descsMatU2D, 1, 1, flop=flop)
         deallocate (tau_Q)
         stats%Flop_Fill = stats%Flop_Fill + flop

         mn = min(N, rank)
         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn, nbslpk, mycol, 0, npcol)
         allocate (tau_Q(myAcols))
         call pgeqrff90(N, mn, matV2D, 1, 1, descsMatV2D, tau_Q, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (RR2(myArows, myAcols))
         RR2 = 0d0
         do myj = 1, myAcols
            call l2g(myj, mycol, rank, npcol, nbslpk, jj)
            do myi = 1, myArows
               call l2g(myi, myrow, rank, nprow, nbslpk, ii)
               if (ii <= jj) RR2(myi, myj) = matV2D(myi, myj)
            enddo
         enddo
         call pun_or_gqrf90(ctxt, matV2D, tau_Q, N, rank, rank, descsMatV2D, 1, 1, flop=flop)
         deallocate (tau_Q)
         stats%Flop_Fill = stats%Flop_Fill + flop

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (mattemp(myArows, myAcols))
         mattemp = 0
         call descinit(descsMatSml, rank, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', rank, rank, rank, cone, RR1, 1, 1, descsMatSml, RR2, 1, 1, descsMatSml, czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         allocate (UUsml(myArows, myAcols), VVsml(myArows, myAcols), Singularsml(rank))
         call PSVD_Truncate(rank, rank, mattemp, descsMatSml, UUsml, VVsml, descsMatSml, descsMatSml, Singularsml, SVD_tolerance, ranknew, ctxt, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyU%blocks(1)%matrix(myArows, myAcols))
         blocks%ButterflyU%blocks(1)%matrix = 0
         call descinit(descsMatU2Dnew, M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', M, ranknew, rank, cone, matU2D, 1, 1, descsMatU2D, UUsml, 1, 1, descsMatSml, czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            blocks%ButterflyU%blocks(1)%matrix(:, myj) = blocks%ButterflyU%blocks(1)%matrix(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyV%blocks(1)%matrix(myArows, myAcols))
         blocks%ButterflyV%blocks(1)%matrix = 0
         call descinit(descsMatV2Dnew, N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', N, ranknew, rank, cone, matV2D, 1, 1, descsMatV2D, VVsml, 1, 1, descsMatSml, czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop

         rank = ranknew

         deallocate (mattemp, RR1, UUsml, VVsml, Singularsml)
         deallocate (RR2)

#else
         mn1 = min(M, rank)
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn1, nbslpk, mycol, 0, npcol)
         call descinit(descsUU_u, M, mn1, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUu(myArows, myAcols))
         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit(descsVV_u, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVu(myArows, myAcols))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(M, rank, matU2D, descsMatU2D, UUu, VVu, descsUU_u, descsVV_u, Singularuv, SVD_tolerance, rank1, ctxt, flop=flop)
         do ii = 1, rank1
            call g2l(ii, rank1, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VVu(myi, :) = VVu(myi, :)*Singularuv(ii)
            endif
         enddo
         deallocate (Singularuv)

         mn2 = min(N, rank)
         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
         call descinit(descsUU_v, N, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUv(myArows, myAcols))
         myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit(descsVV_v, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVv(myArows, myAcols))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(N, rank, matV2D, descsMatV2D, UUv, VVv, descsUU_v, descsVV_v, Singularuv, SVD_tolerance, rank2, ctxt, flop=flop)
         do ii = 1, rank2
            call g2l(ii, rank2, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VVv(myi, :) = VVv(myi, :)*Singularuv(ii)
            endif
         enddo
         deallocate (Singularuv)

         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         allocate (mattemp(myArows, myAcols))
         mattemp = 0
         call descinit(descsMatSml, rank1, rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', rank1, rank2, rank, cone, VVu, 1, 1, descsVV_u, VVv, 1, 1, descsVV_v, czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(min(rank1, rank2), nbslpk, mycol, 0, npcol)
         call descinit(descsUUSml, rank1, min(rank1, rank2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUsml(myArows, myAcols))
         myArows = numroc_wp(min(rank1, rank2), nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         call descinit(descsVVSml, min(rank1, rank2), rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVsml(myArows, myAcols))
         allocate (Singularsml(min(rank1, rank2)))
         call PSVD_Truncate(rank1, rank2, mattemp, descsMatSml, UUsml, VVsml, descsUUSml, descsVVSml, Singularsml, SVD_tolerance, ranknew, ctxt, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyU%blocks(1)%matrix(myArows, myAcols))
         blocks%ButterflyU%blocks(1)%matrix = 0
         call descinit(descsMatU2Dnew, M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', M, ranknew, rank1, cone, UUu, 1, 1, descsUU_u, UUsml, 1, 1, descsUUSml, czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            blocks%ButterflyU%blocks(1)%matrix(:, myj) = blocks%ButterflyU%blocks(1)%matrix(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyV%blocks(1)%matrix(myArows, myAcols))
         blocks%ButterflyV%blocks(1)%matrix = 0
         call descinit(descsMatV2Dnew, N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', N, ranknew, rank2, cone, UUv, 1, 1, descsUU_v, VVsml, 1, 1, descsVVSml, czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop
         rank = ranknew

         deallocate (mattemp, UUsml, VVsml, Singularsml)
         deallocate (UUu, VVu, UUv, VVv)

#endif

      endif
      call MPI_Bcast(rank, 1, MPI_INTEGER, Main_ID, ptree%pgrp(pgno)%Comm, ierr)

      deallocate (matU2D)
      deallocate (matV2D)

      blocks%rankmax = max(blocks%rankmax, rank)
      blocks%rankmin = min(blocks%rankmin, rank)

      deallocate (value_UV)
      deallocate (matU)
      deallocate (matV)

      return

   end subroutine LR_ACA_Parallel

   subroutine LR_BACA(matU, matV, header_m, header_n, M, N, rmax, rank, tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option, error)

      use BPACK_DEFS
      implicit none

      integer rank, rankup, ranknew, row, column, rankmax, N, M, rmax
      DT, allocatable:: row_R(:, :), row_Rtmp(:, :), column_R(:, :), column_RT(:, :), fullmat(:, :), fullmat1(:, :)
      DT::matU(M, rmax), matV(rmax, N)
      DT, allocatable :: core(:, :), core_inv(:, :), tau(:), matUtmp(:, :), matVtmp(:, :)
      real(kind=8):: normA, normUV, flop
      integer itr, itrmax, r_est, Nqr, bsize
      integer, allocatable:: select_column(:), select_row(:), perms(:)
      integer, allocatable :: jpvt(:)

      integer i, j, ii, jj, indx, rank_1, rank_2
      real(kind=8) tolerance, SVD_tolerance
      integer edge_m, edge_n, header_m, header_n, mn
      real(kind=8) inner_UV, n1, n2, a, error

      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer::mrange(M), nrange(N)
      integer::passflag = 0

      n1 = OMP_get_wtime()

      r_est = min(bsize, min(M, N))
      ! r_est=min(M,N)
      itrmax = floor_safe(min(M, N)/dble(r_est))*2

      Nqr = max(max(M, N), r_est)
      allocate (jpvt(Nqr))
      allocate (tau(Nqr))

      allocate (select_column(r_est))
      allocate (select_row(r_est))
      select_column = 0
      select_row = 0

      allocate (row_R(r_est, N))
      allocate (row_Rtmp(r_est, N))
      allocate (column_R(M, r_est))
      allocate (column_RT(r_est, M))
      row_R = 0
      row_Rtmp = 0
      column_R = 0
      column_RT = 0
      allocate (perms(max(M, N)))
      allocate (core(r_est, r_est))
      ! allocate(core_inv(r_est,r_est))
      core = 0
      ! core_inv=0

      normA = 0
      normUV = 0
      itr = 0
      rank = 0

      do while (normUV >= tolerance*normA .and. itr < itrmax)

         !**** create random column index for the first iteration
         if (rank == 0) then
            call rperm(N, perms)
            select_column = perms(1:r_est)
         endif

         !**** Compute columns column_R to find a new set of rows and columns
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo
         call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         ! !$omp parallel do default(shared) private(i,j,edge_m,edge_n)
         ! do i=1,M
         ! do j=1,r_est
         ! edge_m = header_m + i - 1
         ! edge_n = header_n + select_column(j) - 1
         ! call element_Zmn(edge_m,edge_n,column_R(i,j),msh,option,ker)
         ! enddo
         ! enddo
         ! !$omp end parallel do

         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -cone, matU, M, matV(1, select_column(j)), rmax, cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
            enddo
         endif

         !**** Find row pivots from the columns column_R
         call copymatT(column_R, column_RT, M, r_est)
         jpvt = 0
         ! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
         call geqp3f90(column_RT, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_row(1:r_est) = jpvt(1:r_est)

         !**** Compute rows row_R in CUR
         do i = 1, r_est
            mrange(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, N
            nrange(j) = header_n + j - 1
         enddo
         call element_Zmn_block_user(r_est, N, mrange, nrange, row_R, msh, option, ker, 0, passflag, ptree, stats)

         ! !$omp parallel do default(shared) private(i,j,edge_m,edge_n)
         ! do j=1,N
         ! do i=1,r_est
         ! edge_m = header_m + select_row(i) - 1
         ! edge_n = header_n + j - 1
         ! call element_Zmn(edge_m,edge_n,row_R(i,j),msh,option,ker)
         ! enddo
         ! enddo
         ! !$omp end parallel do

         if (rank > 0) then
            do i = 1, r_est
               call gemmf77('N', 'N', 1, N, rank, -cone, matU(select_row(i), 1), M, matV, rmax, cone, row_R(i, 1), r_est)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N, rank)
            enddo
         endif

         !**** Find column pivots from the rows row_R
         jpvt = 0
         row_Rtmp = row_R
         ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
         call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_column(1:r_est) = jpvt(1:r_est)

         !**** Compute columns column_R in CUR
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo
         call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         ! !$omp parallel do default(shared) private(i,j,edge_m,edge_n)
         ! do i=1,M
         ! do j=1,r_est
         ! edge_m = header_m + i - 1
         ! edge_n = header_n + select_column(j) - 1
         ! call element_Zmn(edge_m,edge_n,column_R(i,j),msh,option,ker)
         ! enddo
         ! enddo
         ! !$omp end parallel do
         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -cone, matU, M, matV(1, select_column(j)), rmax, cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
            enddo
         endif

         !**** Compute the skeleton matrix in CUR
         do i = 1, r_est
            core(i, :) = column_R(select_row(i), :)
         enddo

#if 1
         !**** generate column indices for the next iteration
         row_Rtmp = row_R
         do j = 1, r_est
            row_Rtmp(:, select_column(j)) = 0d0
         enddo

         jpvt = 0
         ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
         call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_column(1:r_est) = jpvt(1:r_est)
#endif

         !**** form the LR update by CUR

         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tolerance, SafeUnderflow, ranknew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rankup = ranknew

         if (rankup > 0) then
            row_Rtmp = row_R
            call un_or_mqrf90(core, tau, row_Rtmp, 'L', 'C', r_est, N, rankup, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call trsmf90(core, row_Rtmp, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop

            if (rank + rankup > rmax) rankup = rmax - rank

            ! call assert(rank+rankup<=rmax,'try to increase rmax')
            do j = 1, rankup
               matU(:, rank + j) = column_R(:, jpvt(j))
            enddo
            matV(rank + 1:rank + rankup, :) = row_Rtmp(1:rankup, :)
            rank = rank + rankup

            if (rank == rmax) exit

            !**** update fnorm of UV and matUmatV
            call LR_Fnorm(column_R, row_Rtmp, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            ! if(rankup<8)then ! update fnorm seems more efficienct than recompute fnorm when block size is small
            call LR_FnormUp(matU, matV, M, N, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)
            ! else
            ! call LR_Fnorm(matU,matV,M,N,rank,normA,tolerance*1e-2,Flops=flop)
            ! endif
            stats%Flop_Fill = stats%Flop_Fill + flop

            if (normA > SafeUnderflow) then
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

         itr = itr + 1
      enddo

      ! write(*,*)normUV,normA

      if (rank > 0) then
         call LR_ReCompression(matU, matV, M, N, rank, ranknew, SVD_tolerance, Flops=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rank = ranknew
      else
         rank = 1
         matU(:, 1) = 0
         matV(1, :) = 0
      endif

      deallocate (jpvt)
      deallocate (tau)
      deallocate (select_column)
      deallocate (select_row)
      deallocate (row_R)
      deallocate (row_Rtmp)
      deallocate (column_R)
      deallocate (column_RT)
      deallocate (core)
      ! deallocate(core_inv)
      deallocate (perms)

      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

      return

   end subroutine LR_BACA

   subroutine LR_BACA_noOverlap(matU, matV, header_m, header_n, M, N, rmax, rank, tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option, error)

      use BPACK_DEFS
      implicit none

      integer rank, rankup, ranknew, row, column, rankmax, N, M, rmax
      DT, allocatable:: row_R(:, :), row_R_knn(:, :), row_Rtmp(:, :), row_Rtmp_knn(:, :), column_R(:, :), column_R_knn(:, :), column_RT(:, :), fullmat(:, :), fullmat1(:, :)
      DT::matU(M, rmax), matV(rmax, N)
      DT, allocatable :: core(:, :), core_knn(:, :), core_inv(:, :), tau(:), matUtmp(:, :), matVtmp(:, :)
      real(kind=8):: normA, normUV, flop, maxvalue
      integer itr, itrmax, r_est, r_est_knn_r, r_est_knn, r_est_knn_c, r_est_tmp, Nqr, bsize
      integer, allocatable:: select_column(:), select_column_knn(:), select_row_knn(:), select_column1(:), select_row(:), perms(:), rows(:), columns(:)
      integer, allocatable :: jpvt(:)

      integer i, j, ii, jj, iii, jjj, indx, rank_1, rank_2
      real(kind=8) tolerance, SVD_tolerance
      integer edge_m, edge_n, header_m, header_n, mn
      real(kind=8) inner_UV, n1, n2, a, error

      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer::mrange(M), nrange(N)
      integer::passflag = 0

      n1 = OMP_get_wtime()

      r_est = min(bsize, min(M, N))
      ! r_est=min(M,N)
      itrmax = floor_safe(min(M, N)/dble(r_est))*2

      Nqr = max(max(M, N), r_est)
      allocate (jpvt(Nqr))
      allocate (tau(Nqr))

      allocate (columns(N))
      columns = 0
      allocate (rows(M))
      rows = 0

      allocate (select_column(r_est))
      allocate (select_column1(r_est))
      allocate (select_row(r_est))
      select_column = 0
      select_row = 0

      allocate (row_R(r_est, N))
      allocate (row_Rtmp(r_est, N))
      allocate (column_R(M, r_est))
      allocate (column_RT(r_est, M))
      row_R = 0
      row_Rtmp = 0
      column_R = 0
      column_RT = 0
      allocate (perms(max(M, N)))
      allocate (core(r_est, r_est))
      ! allocate(core_inv(r_est,r_est))
      core = 0
      ! core_inv=0

      normA = 0
      normUV = 0
      itr = 0
      rank = 0

      !**** if nearest neighbour is available, select them first
      if (option%knn > 0) then
         allocate (select_column_knn(M*option%knn))
         allocate (select_row_knn(N*option%knn))
         r_est_knn_r = 0
         r_est_knn_c = 0

         do i = 1, M
            edge_m = header_m + i - 1
            do iii = 1, option%knn
               if (msh%nns(edge_m, iii) >= header_n .and. msh%nns(edge_m, iii) <= header_n + N - 1) then
                  r_est_knn_c = r_est_knn_c + 1
                  select_column_knn(r_est_knn_c) = msh%nns(edge_m, iii) + 1 - header_n
               endif
            enddo
         enddo
         r_est_tmp = r_est_knn_c
         if (r_est_knn_c > 0) call remove_dup_int(select_column_knn, r_est_tmp, r_est_knn_c)
         allocate (column_R_knn(M, r_est_knn_c))

         do j = 1, N
            edge_n = header_n + j - 1
            do jjj = 1, option%knn
               if (msh%nns(edge_n, jjj) >= header_m .and. msh%nns(edge_n, jjj) <= header_m + M - 1) then
                  r_est_knn_r = r_est_knn_r + 1
                  select_row_knn(r_est_knn_r) = msh%nns(edge_n, jjj) + 1 - header_m
               endif
            enddo
         enddo
         r_est_tmp = r_est_knn_r
         if (r_est_knn_r > 0) call remove_dup_int(select_row_knn, r_est_tmp, r_est_knn_r)
         allocate (row_R_knn(r_est_knn_r, N))
         allocate (row_Rtmp_knn(r_est_knn_r, N))
         allocate (core_knn(r_est_knn_r, r_est_knn_c))
         if (r_est_knn_r > 0 .and. r_est_knn_c > 0) then
            do i = 1, M
               mrange(i) = header_m + i - 1
            enddo
            do j = 1, r_est_knn_c
               nrange(j) = header_n + select_column_knn(j) - 1
            enddo
            call element_Zmn_block_user(M, r_est_knn_c, mrange, nrange, column_R_knn, msh, option, ker, 0, passflag, ptree, stats)

            do i = 1, r_est_knn_r
               mrange(i) = header_m + select_row_knn(i) - 1
            enddo
            do j = 1, N
               nrange(j) = header_n + j - 1
            enddo
            call element_Zmn_block_user(r_est_knn_r, N, mrange, nrange, row_R_knn, msh, option, ker, 0, passflag, ptree, stats)

            !**** Compute the skeleton matrix in CUR
            do i = 1, r_est_knn_r
               core_knn(i, :) = column_R_knn(select_row_knn(i), :)
            enddo

            jpvt = 0
            call geqp3modf90(core_knn, jpvt, tau, tolerance, SafeUnderflow, ranknew, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            rankup = ranknew
            if (rankup > 0) then

               if (rank + rankup > rmax) rankup = rmax - rank

               row_Rtmp_knn = row_R_knn
               call un_or_mqrf90(core_knn, tau, row_Rtmp_knn, 'L', 'C', r_est_knn_r, N, rankup, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop
               call trsmf90(core_knn, row_Rtmp_knn, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop

               columns(rank + 1:rankup + rank) = select_column_knn(jpvt(1:rankup))
               rows(rank + 1:rankup + rank) = select_row_knn(1:rankup)

               ! call assert(rank+rankup<=rmax,'try to increase rmax')
               do j = 1, rankup
                  matU(:, rank + j) = column_R_knn(:, jpvt(j))
               enddo
               matV(rank + 1:rank + rankup, :) = row_Rtmp_knn(1:rankup, :)

               rank = rank + rankup

               if (rank == rmax) goto 10 !*** skip ACA iteration

               !**** update fnorm of UV and matUmatV
               call LR_Fnorm(column_R_knn, row_Rtmp_knn, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop
               call LR_FnormUp(matU, matV, M, N, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)

               stats%Flop_Fill = stats%Flop_Fill + flop

               if (normA > SafeUnderflow) then
                  error = normUV/normA
               else
                  error = 0
               endif

               !**** Find column pivots for the next iteration
               jpvt = 0
               row_Rtmp_knn = row_R_knn
               if (rank > 0) row_Rtmp_knn(:, columns(1:rank)) = 0
               ! call geqp3modf90(row_Rtmp_knn,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
               call geqp3f90(row_Rtmp_knn, jpvt, tau, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop
               select_column(1:r_est) = jpvt(1:r_est)
            else
               error = 0
               goto 20   !*** no effective rank found using KNN, go to ACA iteration
            endif
         endif
      endif

20    do while (normUV >= tolerance*normA .and. itr < itrmax)

         !**** create random column index for the first iteration
         if (rank == 0) then
            call rperm(N, perms)
            select_column = perms(1:r_est)
         endif

         select_column1 = select_column

         !**** Compute columns column_R to find a new set of rows and columns
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo
         call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -cone, matU, M, matV(1, select_column(j)), rmax, cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
               ! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-cone,cone)
            enddo
         endif

         !**** Find row pivots from the columns column_R
         call copymatT(column_R, column_RT, M, r_est)
         if (rank > 0) column_RT(:, rows(1:rank)) = 0
         jpvt = 0
         ! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
         call geqp3f90(column_RT, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_row(1:r_est) = jpvt(1:r_est)

         !**** Compute rows row_R in CUR
         ! !$omp end parallel do
         do i = 1, r_est
            mrange(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, N
            nrange(j) = header_n + j - 1
         enddo
         call element_Zmn_block_user(r_est, N, mrange, nrange, row_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do i = 1, r_est
               call gemmf77('N', 'N', 1, N, rank, -cone, matU(select_row(i), 1), M, matV, rmax, cone, row_R(i, 1), r_est)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N, rank)
               ! call gemmf90(matU(select_row(i),1), M,matV,rmax,row_R(i,1),r_est,'N','N',1,N,rank,-cone,cone)
            enddo
         endif

         !**** Compute the skeleton matrix in CUR
         do i = 1, r_est
            core(i, :) = column_R(select_row(i), :)
         enddo
         maxvalue = abs(core(1, 1))

         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tolerance, SafeUnderflow, ranknew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rankup = ranknew

         if (rankup > 0) then

            row_Rtmp = row_R
            call un_or_mqrf90(core, tau, row_Rtmp, 'L', 'C', r_est, N, rankup, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call trsmf90(core, row_Rtmp, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop

            if (rank + rankup > rmax) rankup = rmax - rank

            columns(rank + 1:rankup + rank) = select_column1(jpvt(1:rankup))
            rows(rank + 1:rankup + rank) = select_row(1:rankup)

            ! call assert(rank+rankup<=rmax,'try to increase rmax')
            do j = 1, rankup
               matU(:, rank + j) = column_R(:, jpvt(j))
            enddo
            matV(rank + 1:rank + rankup, :) = row_Rtmp(1:rankup, :)

            rank = rank + rankup

            if (rank == rmax) exit

            !**** update fnorm of UV and matUmatV
            call LR_Fnorm(column_R, row_Rtmp, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call LR_FnormUp(matU, matV, M, N, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)

            stats%Flop_Fill = stats%Flop_Fill + flop

            if (normA > SafeUnderflow) then
               error = normUV/normA
            else
               error = 0
               exit
            endif

            !**** Find column pivots for the next iteration
            jpvt = 0
            row_Rtmp = row_R
            if (rank > 0) row_Rtmp(:, columns(1:rank)) = 0
            ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,SafeUnderflow,ranknew)
            call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            select_column(1:r_est) = jpvt(1:r_est)
         else
            error = 0
            exit
         endif

         ! write(*,*)itr,rank,rankup,error

         itr = itr + 1
      enddo

      ! write(*,*)normUV,normA

10    if (allocated(row_R_knn)) deallocate (row_R_knn)
      if (allocated(row_Rtmp_knn)) deallocate (row_Rtmp_knn)
      if (allocated(column_R_knn)) deallocate (column_R_knn)
      if (allocated(core_knn)) deallocate (core_knn)
      if (allocated(select_column_knn)) deallocate (select_column_knn)
      if (allocated(select_row_knn)) deallocate (select_row_knn)

      if (rank > 0) then
         call LR_ReCompression(matU, matV, M, N, rank, ranknew, SVD_tolerance, Flops=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rank = ranknew
      else
         rank = 1
         matU(:, 1) = 0
         matV(1, :) = 0
      endif

      deallocate (jpvt)
      deallocate (tau)
      deallocate (select_column)
      deallocate (select_column1)
      deallocate (select_row)
      deallocate (row_R)
      deallocate (row_Rtmp)
      deallocate (column_R)
      deallocate (column_RT)
      deallocate (core)
      ! deallocate(core_inv)
      deallocate (perms)
      deallocate (columns)
      deallocate (rows)

      n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

      return

   end subroutine LR_BACA_noOverlap

   subroutine LR_SeudoSkeleton(blocks, header_m, header_n, M, N, rmaxc, rmaxr, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, ctxt, pgno)

      use BPACK_DEFS
      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2, rank
      integer index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n
      integer ranknew, row, column, rankmax, N, M, rankmax_min, rmax, rmaxc, rmaxr, idxs_r, idxs_c, flag0, myAcols, myArows, npcol, nprow, rank_new, nproc
      DT value_Z, value_UV, maxvalue
      DT inner_U, inner_V, ctemp
      real(kind=8) inner_UV, flop
      DT, allocatable::matU(:, :), matV(:, :), matV_tmp(:, :), UU(:, :), VV(:, :), MatrixSubselection(:, :)
      real(kind=8), allocatable::Singular(:)
      DT, allocatable:: row_R(:), column_R(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(matrixblock)::blocks
      type(Hstat)::stats

      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), matU2D(:, :), matV2D(:, :), matU1D(:, :), matV1D(:, :)
      real(kind=8), allocatable :: Singularsml(:)
      type(proctree)::ptree
      type(Hoption)::option
      integer ctxt, myrow, mycol, iproc, jproc, myi, myj, mnmax
      integer, optional::pgno
      integer LWORK, LRWORK, INFO, ierr
      DT:: TEMP(1)
      real(kind=8), allocatable::RWORK(:)
      integer::descsmatU(9), descsmatV(9), descsmatVQ(9), descsub(9), descUU(9), descVV(9), descButterU2D(9), descButterV2D(9), descButterU1D(9), descButterV1D(9), descQ(9), descR(9)

      integer nb1Dr, nb1Dc
      integer, allocatable::select_col(:), select_row(:), N_p(:, :), M_p(:, :)
      DT, allocatable:: WORK(:)
      integer, allocatable :: ipiv(:), jpiv(:), JPERM(:)
      integer:: mrange(M), nrange(N)
      DT, allocatable :: tau(:)
      real(kind=8):: RTEMP(1)
      integer:: passflag = 0
      rank_new = 0
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)

      if (nprow*npcol == 1) then

         rmax = min(rmaxc, rmaxr)
         allocate (select_col(rmaxc))
         allocate (select_row(rmaxr))
         call linspaceI(1, M, rmaxr, select_row)
         call linspaceI(1, N, rmaxc, select_col)

         allocate (matV(N, rmaxr))
         allocate (matV_tmp(rmaxr, N))
         do ii = 1, rmaxr
            mrange(ii) = header_m + select_row(ii) - 1
         enddo
         do jj = 1, N
            nrange(jj) = header_n + jj - 1
         enddo
         call element_Zmn_block_user(rmaxr, N, mrange, nrange, matV_tmp, msh, option, ker, 0, passflag, ptree, stats)
         call copymatT(matV_tmp, matV, rmaxr, N)
         deallocate (matV_tmp)

         ! matV=0
         ! do ii=1,N
         ! do jj=1,rmaxr
         ! edge_m = header_m + select_row(jj) - 1
         ! edge_n = header_n + ii - 1
         ! call element_Zmn(edge_m,edge_n,matV(ii,jj),msh,option,ker)
         ! enddo
         ! enddo

         allocate (MatrixSubselection(rmaxr, rmaxc))
         MatrixSubselection = 0
         do ii = 1, rmaxr
            mrange(ii) = header_m + select_row(ii) - 1
         enddo
         do jj = 1, rmaxc
            nrange(jj) = header_n + select_col(jj) - 1
         enddo
         call element_Zmn_block_user(rmaxr, rmaxc, mrange, nrange, MatrixSubselection, msh, option, ker, 0, passflag, ptree, stats)

         ! do ii=1,rmaxr
         ! do jj=1,rmaxc
         ! edge_m = header_m + select_row(ii) - 1
         ! edge_n = header_n + select_col(jj) - 1
         ! call element_Zmn(edge_m,edge_n,MatrixSubselection(ii,jj),msh,option,ker)
         ! enddo
         ! enddo

         allocate (jpiv(rmaxc))
         jpiv = 0
         allocate (tau(min(rmaxr, rmaxc)))
         tau = 0
         call geqp3modf90(MatrixSubselection, jpiv, tau, tolerance, SafeUnderflow, rank_new, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rank = rank_new
         if (rank > 0) then
            matV = conjg(cmplx(matV, kind=8))
            call un_or_mqrf90(MatrixSubselection, tau, matV, 'R', 'N', N, rmaxr, rank_new, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            matV = conjg(cmplx(matV, kind=8))
            call trsmf90(MatrixSubselection, matV, 'R', 'U', 'T', 'N', N, rank_new, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
         else
            rank = 1
            rank_new = 1
            matV = 0
         endif
         allocate (blocks%ButterflyU%blocks(1)%matrix(M, rank))
         blocks%ButterflyU%blocks(1)%matrix = 0

         do ii = 1, M
            mrange(ii) = header_m + ii - 1
         enddo
         do jj = 1, rank
            nrange(jj) = header_n + select_col(jpiv(jj)) - 1
         enddo
         call element_Zmn_block_user(M, rank, mrange, nrange, blocks%ButterflyU%blocks(1)%matrix, msh, option, ker, 0, passflag, ptree, stats)

         ! do ii=1,M
         ! do jj=1,rank
         ! edge_m = header_m + ii - 1
         ! edge_n = header_n + select_col(jpiv(jj)) - 1
         ! call element_Zmn(edge_m,edge_n,blocks%ButterflyU%blocks(1)%matrix(ii,jj),msh,option,ker)
         ! enddo
         ! enddo

         allocate (blocks%ButterflyV%blocks(1)%matrix(N, rank))
         blocks%ButterflyV%blocks(1)%matrix = matV(1:N, 1:rank)

         blocks%rankmax = rank
         blocks%rankmin = rank

         deallocate (matV)
         deallocate (jpiv)
         deallocate (tau)
         deallocate (select_col)
         deallocate (select_row)
         deallocate (MatrixSubselection)
      else
         if (myrow /= -1 .and. mycol /= -1) then
            rmax = min(rmaxc, rmaxr)
            allocate (select_col(rmaxc))
            allocate (select_row(rmaxr))
            call linspaceI(1, M, rmaxr, select_row)
            call linspaceI(1, N, rmaxc, select_col)

            call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
            ! nproc = npcol*nprow
            myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rmaxr, nbslpk, mycol, 0, npcol)
            allocate (matV(myArows, myAcols))
            allocate (matV_tmp(myAcols, myArows))
            call descinit(descsmatV, N, rmaxr, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit fail for descsmatV')
            matV = 0
            do myj = 1, myAcols
               call l2g(myj, mycol, rmaxr, npcol, nbslpk, jj)
               mrange(myj) = header_m + select_row(jj) - 1
            enddo
            do myi = 1, myArows
               call l2g(myi, myrow, N, nprow, nbslpk, ii)
               nrange(myi) = header_n + ii - 1
            enddo
            call element_Zmn_block_user(myAcols, myArows, mrange, nrange, matV_tmp, msh, option, ker, 0, passflag, ptree, stats)
            call copymatT(matV_tmp, matV, myAcols, myArows)
            deallocate (matV_tmp)

            ! do myi=1,myArows
            ! call l2g(myi,myrow,N,nprow,nbslpk,ii)
            ! do myj=1,myAcols
            ! call l2g(myj,mycol,rmaxr,npcol,nbslpk,jj)
            ! edge_m = header_m + select_row(jj) - 1
            ! edge_n = header_n + ii - 1
            ! call element_Zmn(edge_m,edge_n,matV(myi,myj),msh,option,ker)
            ! enddo
            ! enddo

            call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
            myArows = numroc_wp(rmaxr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rmaxc, nbslpk, mycol, 0, npcol)

            allocate (MatrixSubselection(myArows, myAcols))
            call descinit(descsub, rmaxr, rmaxc, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit fail for descsub')
            MatrixSubselection = 0
            do myi = 1, myArows
               call l2g(myi, myrow, rmaxr, nprow, nbslpk, ii)
               mrange(myi) = header_m + select_row(ii) - 1
            enddo
            do myj = 1, myAcols
               call l2g(myj, mycol, rmaxc, npcol, nbslpk, jj)
               nrange(myj) = header_n + select_col(jj) - 1
            enddo
            call element_Zmn_block_user(myArows, myAcols, mrange, nrange, MatrixSubselection, msh, option, ker, 0, passflag, ptree, stats)

            ! do myi=1,myArows
            ! call l2g(myi,myrow,rmaxr,nprow,nbslpk,ii)
            ! do myj=1,myAcols
            ! call l2g(myj,mycol,rmaxc,npcol,nbslpk,jj)
            ! edge_m = header_m + select_row(ii) - 1
            ! edge_n = header_n + select_col(jj) - 1
            ! call element_Zmn(edge_m,edge_n,MatrixSubselection(myi,myj),msh,option,ker)
            ! enddo
            ! enddo

            ! Compute QR of MatrixSubselection*P
            allocate (ipiv(myAcols))
            ipiv = 0
            allocate (tau(myAcols))
            tau = 0
            allocate (jpiv(rmaxc))
            jpiv = 0
            allocate (JPERM(rmaxc))
            JPERM = 0
            call pgeqpfmodf90(rmaxr, rmaxc, MatrixSubselection, 1, 1, descsub, ipiv, tau, JPERM, jpiv, rank_new, tolerance, SafeUnderflow, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

            rank = rank_new

            if (rank > 0) then

               ! Compute matV*conjg(Q)
               matV = conjg(cmplx(matV, kind=8))
               call pun_or_mqrf90('R', 'N', N, rmaxr, rank_new, MatrixSubselection, 1, 1, descsub, tau, matV, 1, 1, descsmatV)
               matV = conjg(cmplx(matV, kind=8))

               ! Compute matV*conjg(Q)*(R^T)^-1
               call ptrsmf90('R', 'U', 'T', 'N', N, rank_new, cone, MatrixSubselection, 1, 1, descsub, matV, 1, 1, descsmatV, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

               call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
               myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
               allocate (matV2D(myArows, myAcols))
               call descinit(descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
               call assert(info == 0, 'descinit fail for descButterV2D')
               matV2D(1:myArows, 1:myAcols) = matV(1:myArows, 1:myAcols)

            else
               rank = 1
               rank_new = 1

               call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
               myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
               allocate (matV2D(myArows, myAcols))
               call descinit(descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
               call assert(info == 0, 'descinit fail for descButterV2D')
               matV2D(1:myArows, 1:myAcols) = 0
            endif

            call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
            myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
            allocate (matU2D(myArows, myAcols))
            call descinit(descButterU2D, M, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit fail for descButterU2D')
            matU2D = 0

            do myi = 1, myArows
               call l2g(myi, myrow, M, nprow, nbslpk, ii)
               mrange(myi) = header_m + ii - 1
            enddo
            do myj = 1, myAcols
               call l2g(myj, mycol, rank_new, npcol, nbslpk, jj)
               nrange(myj) = header_n + select_col(ipiv(myj)) - 1
            enddo
            call element_Zmn_block_user(myArows, myAcols, mrange, nrange, matU2D, msh, option, ker, 0, passflag, ptree, stats)

            ! do myi=1,myArows
            ! call l2g(myi,myrow,M,nprow,nbslpk,ii)
            ! do myj=1,myAcols
            ! call l2g(myj,mycol,rank_new,npcol,nbslpk,jj)
            ! edge_m = header_m + ii - 1
            ! edge_n = header_n + select_col(ipiv(myj)) - 1
            ! call element_Zmn(edge_m,edge_n,matU2D(myi,myj),msh,option,ker)
            ! enddo
            ! enddo

            deallocate (select_col)
            deallocate (select_row)
            deallocate (matV)
            deallocate (MatrixSubselection)
            deallocate (ipiv)
            deallocate (jpiv)
            deallocate (JPERM)
            deallocate (tau)
         endif

         ! call MPI_ALLREDUCE(rank_new,rank,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)

         if (.not. present(pgno)) then
            myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (blocks%ButterflyU%blocks(1)%matrix(myArows, myAcols))
            blocks%ButterflyU%blocks(1)%matrix(1:myArows, 1:myAcols) = matU2D(1:myArows, 1:myAcols)
            deallocate (matU2D)

            myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (blocks%ButterflyV%blocks(1)%matrix(myArows, myAcols))
            blocks%ButterflyV%blocks(1)%matrix(1:myArows, 1:myAcols) = matV2D(1:myArows, 1:myAcols)
            deallocate (matV2D)

            blocks%rankmax = rank
            blocks%rankmin = rank

         else

            ! distribute UV factor from 2D grid into 1D grid conformal to leaf sizes

            allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc, rank))
            blocks%ButterflyU%blocks(1)%matrix = 0

            allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc, rank))
            blocks%ButterflyV%blocks(1)%matrix = 0

            blocks%rankmax = max(blocks%rankmax, rank)
            blocks%rankmin = min(blocks%rankmin, rank)

            if (.not. allocated(matU2D)) allocate (matU2D(1, 1))  ! Required for Redistribute2Dto1D
            if (.not. allocated(matV2D)) allocate (matV2D(1, 1))

            call Redistribute2Dto1D(matU2D, M, 0, pgno, blocks%ButterflyU%blocks(1)%matrix, blocks%M_p, 0, pgno, rank, ptree)
            call Redistribute2Dto1D(matV2D, N, 0, pgno, blocks%ButterflyV%blocks(1)%matrix, blocks%N_p, 0, pgno, rank, ptree)

            if (allocated(matU2D)) deallocate (matU2D)
            if (allocated(matV2D)) deallocate (matV2D)

         end if
      endif

      return

   end subroutine LR_SeudoSkeleton

!!!!!!! check error of Bplus compression of blocks by comparing the full block, assuming final distribution
   subroutine Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
      use BPACK_DEFS
      implicit none

      type(blockplus)::bplus
      type(matrixblock), pointer::blocks
      type(mesh)::msh
      type(Hoption)::option
      type(kernelquant)::ker
      type(Hstat)::stats
      type(proctree)::ptree
      integer pgno, pp
      integer Ntest
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, M1, N1, M2, N2, rank1, rank2, ierr, MyID
      integer:: cridx, info
      DT, allocatable:: Vin_glo(:, :), Vin(:, :), Vout1(:, :), Vout2(:, :), Vout_glo(:, :), Vinter(:, :), Fullmat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank, head_n
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      integer, allocatable::mrange(:), nrange(:)
      integer::mrange_dummy(1), nrange_dummy(1)
      DT:: mat_dummy(1, 1)
      integer:: passflag = 0

      blocks => bplus%LL(1)%matrices_block(1)
      pgno = blocks%pgno
      pp = ptree%MyID - ptree%pgrp(pgno)%head + 1

      Ntest = 32
      ! Ntest=blocks%N/2
      allocate (Vin_glo(blocks%N, Ntest))
      if (pp == 1) then
      do ii = 1, blocks%N
      do jj = 1, Ntest
         call random_dp_number(Vin_glo(ii, jj))
      end do
      end do
      endif
      call MPI_Bcast(Vin_glo, blocks%N*Ntest, MPI_DT, Main_ID, ptree%pgrp(pgno)%Comm, ierr)

      allocate (Vin(blocks%N_loc, Ntest))
      head_n = blocks%N_p(pp, 1) - 1

      Vin(:, :) = Vin_glo(head_n + 1:head_n + blocks%N_loc, :)

      allocate (Vout1(blocks%M_loc, Ntest))
      allocate (Vout2(blocks%M_loc, Ntest))
      Vout1 = 0
      Vout2 = 0

      call Bplus_block_MVP_dat(bplus, 'N', blocks%M_loc, blocks%N_loc, Ntest, Vin, Vout2, cone, czero, ptree, stats)

      allocate (Fullmat(blocks%M_loc, blocks%N))
      allocate (mrange(blocks%M_loc))
      allocate (nrange(blocks%N))
      do myi = 1, blocks%M_loc
         mrange(myi) = blocks%M_p(pp, 1) + blocks%headm - 1 + myi - 1
      enddo
      do myj = 1, blocks%N
         nrange(myj) = blocks%headn + myj - 1
      enddo
      call element_Zmn_block_user(blocks%M_loc, blocks%N, mrange, nrange, Fullmat, msh, option, ker, 0, passflag, ptree, stats)
      deallocate (mrange)
      deallocate (nrange)

      passflag = 0
      do while (passflag == 0)
         call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
      enddo

      call gemmf90(Fullmat, blocks%M_loc, Vin_glo, blocks%N, Vout1, blocks%M_loc, 'N', 'N', blocks%M_loc, Ntest, blocks%N, cone, czero)

      Vout2 = Vout2 - Vout1

      rtemp1 = fnorm(Vout2, blocks%M_loc, Ntest)**2d0
      rtemp0 = fnorm(Vout1, blocks%M_loc, Ntest)**2d0

      deallocate (Vout2, Vout1, Vin_glo)

      call MPI_ALLREDUCE(rtemp0, fnorm0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(rtemp1, fnorm1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

      call MPI_Comm_rank(ptree%pgrp(pgno)%Comm, MyID, ierr)

      if (MyID == Main_ID) then
         write (*, *) blocks%row_group, blocks%col_group, 'CheckError_Full error:', sqrt(fnorm1/fnorm0)
      endif

   end subroutine Bplus_CheckError_Full

   subroutine LocalButterflySVD_Left(index_i_loc, index_j_loc, level_loc, level_butterflyL, level, index_i_m, blocks, option, msh, ButterflyP_old, ButterflyP)
      use MISC_Utilities
      implicit none
      integer index_i_loc, index_j_loc, level_loc, level_butterflyL, index_i_m, index_i, index_j, level, group_m, mm, nn, nn1, nn2, j, i, mn, rank, mm1
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      real(kind=8), allocatable :: Singular(:)
      type(matrixblock)::blocks
      real(kind=8):: SVD_tolerance
      type(Hoption)::option
      type(mesh)::msh

      index_j = index_j_loc
      index_i = (index_i_m - 1)*(2**level_loc) + index_i_loc
      group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
      group_m = group_m*2**level - 1 + index_i
      mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1

      ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
      if (size(ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix, 1) /= mm) then
         write (*, *) 'mm incorrect'
         stop
      end if
      nn1 = size(ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix, 2)
      nn2 = size(ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix, 2)
      nn = nn1 + nn2

      allocate (QQ(mm, nn))
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn1
         do i = 1, mm
            QQ(i, j) = ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn2
         do i = 1, mm
            QQ(i, j + nn1) = ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      ! write(*,*)'dddd',fnorm(QQ,mm,nn)

      mn = min(mm, nn)
      allocate (UU(mm, mn), VV(mn, nn), Singular(mn))
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, rank)
      ! rank = min(rank,37)

      ! rank = 7
      ! write(*,*)'dddd', rank

      ! rankmax_for_butterfly(level_loc)=max(rank,rankmax_for_butterfly(level_loc))
      ! rankmin_for_butterfly(level_loc)=min(rank,rankmin_for_butterfly(level_loc))

      ! blocks%rankmax = max(blocks%rankmax,rank)
      ! blocks%rankmin = min(blocks%rankmin,rank)

      allocate (mat_tmp(mm, rank))
      ! !$omp parallel do default(shared) private(i,j,k,ctemp)
      do j = 1, rank
         do i = 1, mm
            mat_tmp(i, j) = UU(i, j)*Singular(j)
         enddo
      enddo
      ! !$omp end parallel do

      if (level_loc /= level_butterflyL) then
         mm1 = msh%basis_group(group_m*2)%tail - msh%basis_group(group_m*2)%head + 1
         allocate (ButterflyP%blocks(2*index_i_loc - 1, index_j_loc)%matrix(mm1, rank))
         allocate (ButterflyP%blocks(2*index_i_loc, index_j_loc)%matrix(mm - mm1, rank))

         ButterflyP%blocks(2*index_i_loc - 1, index_j_loc)%matrix = mat_tmp(1:mm1, 1:rank)
         ButterflyP%blocks(2*index_i_loc, index_j_loc)%matrix = mat_tmp(1 + mm1:mm, 1:rank)
      else
         allocate (blocks%ButterflyU%blocks(index_i)%matrix(mm, rank))
         ! !$omp parallel do default(shared) private(i,j)
         do j = 1, rank
            do i = 1, mm
               blocks%ButterflyU%blocks(index_i)%matrix(i, j) = mat_tmp(i, j)
               ! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
            enddo
         enddo
         ! !$omp end parallel do

      end if

      allocate (blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix(rank, nn1))
      allocate (blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix(rank, nn2))

      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, rank
         do j = 1, nn1
            blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix(i, j) = VV(i, j)
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, rank
         do j = 1, nn2
            blocks%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix(i, j) = VV(i, j + nn1)
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
         enddo
      enddo
      ! !$omp end parallel do

      deallocate (QQ, UU, VV, Singular, mat_tmp)

      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix)/1024.0d3
      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix)/1024.0d3
      ! if (level_loc==level_butterflyL) then
      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyU%blocks(index_i)%matrix)/1024.0d3
      ! endif

   end subroutine LocalButterflySVD_Left

   subroutine LocalButterflySVD_Right(index_i_loc, index_j_loc, level_loc, level_butterflyR, level, level_butterfly, index_j_m, blocks, option, msh, ButterflyP_old, ButterflyP)
      use MISC_Utilities
      implicit none
      integer index_i_loc, index_j_loc, level_loc, level_butterflyR, level_butterfly, index_j_m, index_i, index_j, level, group_n, mm, nn, nn1, mm1, mm2, j, i, mn, rank
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      real(kind=8), allocatable :: Singular(:)
      type(matrixblock)::blocks
      real(kind=8):: SVD_tolerance
      type(Hoption)::option
      type(mesh)::msh

      index_i = index_i_loc
      index_j = (index_j_m - 1)*(2**level_loc) + index_j_loc
      group_n = blocks%col_group   ! Note: row_group and col_group interchanged here
      group_n = group_n*2**(level_butterfly - level + 1) - 1 + index_j
      nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1

      ! call assert(size(ButterflyP_old%blocks(index_i,2*index_j-1)%matrix,1)==mm,'mm incorrect')
      if (size(ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix, 2) /= nn) then
         write (*, *) 'nn incorrect'
         stop
      end if
      mm1 = size(ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix, 1)
      mm2 = size(ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix, 1)
      mm = mm1 + mm2

      !!!!!!!!!!!!!!!!!!

      allocate (QQ(mm, nn))
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, mm1
            QQ(i, j) = ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, mm2
            QQ(i + mm1, j) = ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      mn = min(mm, nn)
      allocate (UU(mm, mn), VV(mn, nn), Singular(mn))
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, rank)
      ! rank = min(rank,37)

      ! rank = 7
      ! write(*,*)rank

      ! rankmax_for_butterfly(-level_loc)=max(rank,rankmax_for_butterfly(-level_loc))
      ! rankmin_for_butterfly(-level_loc)=min(rank,rankmin_for_butterfly(-level_loc))
      ! blocks%rankmax = max(blocks%rankmax,rank)
      ! blocks%rankmin = min(blocks%rankmin,rank)

      allocate (mat_tmp(rank, nn))
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, rank
            mat_tmp(i, j) = VV(i, j)*Singular(i)
         enddo
      enddo
      ! !$omp end parallel do

      if (level_loc /= level_butterflyR) then
         nn1 = msh%basis_group(group_n*2)%tail - msh%basis_group(group_n*2)%head + 1
         allocate (ButterflyP%blocks(index_i_loc, 2*index_j_loc - 1)%matrix(rank, nn1))
         allocate (ButterflyP%blocks(index_i_loc, 2*index_j_loc)%matrix(rank, nn - nn1))

         ButterflyP%blocks(index_i_loc, 2*index_j_loc - 1)%matrix = mat_tmp(1:rank, 1:nn1)
         ButterflyP%blocks(index_i_loc, 2*index_j_loc)%matrix = mat_tmp(1:rank, 1 + nn1:nn)
      else
         allocate (blocks%ButterflyV%blocks(index_j)%matrix(nn, rank))
         ! !$omp parallel do default(shared) private(i,j)
         do j = 1, rank
            do i = 1, nn
               blocks%ButterflyV%blocks(index_j)%matrix(i, j) = mat_tmp(j, i)
               ! blocks%ButterflyU%blocks(index_i)%matrix(i,j)=random_dp_number()
            enddo
         enddo
         ! !$omp end parallel do

      end if

      allocate (blocks%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix(mm1, rank))
      allocate (blocks%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix(mm2, rank))

      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, mm1
         do j = 1, rank
            blocks%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix(i, j) = UU(i, j)
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix(i,j)=random_dp_number()
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, mm2
         do j = 1, rank
            blocks%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix(i, j) = UU(i + mm1, j)
            ! blocks%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix(i,j)=random_dp_number()
         enddo
      enddo
      ! !$omp end parallel do

      deallocate (QQ, UU, VV, Singular, mat_tmp)

      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix)/1024.0d3
      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix)/1024.0d3
      ! if (level_loc==level_butterflyR) then
      ! memory_butterfly=memory_butterfly+SIZEOF(blocks%ButterflyV%blocks(index_j)%matrix)/1024.0d3
      ! endif

   end subroutine LocalButterflySVD_Right

   subroutine Full_construction(blocks, msh, ker, stats, option, ptree)

      use BPACK_DEFS
      implicit none

      integer group_m, group_n, i, j
      integer mm, nn
      integer head_m, head_n, tail_m, tail_n
      DT value_Z
      type(matrixblock)::blocks
      type(mesh)::msh
      type(Hoption)::option
      type(proctree)::ptree
      type(Hstat)::stats
      type(kernelquant)::ker
      integer, allocatable::mrange(:), nrange(:)
      integer passflag

      mm = blocks%M
      head_m = blocks%headm
      tail_m = mm + head_m - 1
      nn = blocks%N
      head_n = blocks%headn
      tail_n = nn + head_n - 1
      allocate (mrange(mm))
      do i = head_m, tail_m
         mrange(i - head_m + 1) = i
      enddo
      allocate (nrange(nn))
      do j = head_n, tail_n
         nrange(j - head_n + 1) = j
      enddo

      allocate (blocks%fullmat(mm, nn))
      if (blocks%row_group == blocks%col_group) allocate (blocks%ipiv(mm))

      call element_Zmn_block_user(mm, nn, mrange, nrange, blocks%fullmat, msh, option, ker, 0, passflag, ptree, stats)

      deallocate (mrange)
      deallocate (nrange)

      return

   end subroutine Full_construction

end module Bplus_compress
