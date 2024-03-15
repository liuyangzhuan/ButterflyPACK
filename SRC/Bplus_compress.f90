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

!> @file Bplus_compress.f90
!> @brief Block-level subroutines for constructing a LR/Butterfly/Dense block from entry evaluations


#include "ButterflyPACK_config.fi"
module Bplus_compress
   use BPACK_DEFS
   use MISC_Utilities
   use Bplus_randomizedop
   use BPACK_structure
   use BPACK_Utilities
! use element_Z

contains

   subroutine BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)


      implicit none

      integer Nboundall, Ninadmissible, statflag
      integer boundary_map(:,:)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new, rank_new1
      integer group_m, group_n, mm, nn, index_i, index_i_loc, index_j_loc, index_j, ii, jj, ij
      integer level, length_1, length_2, level_blocks, index_ij
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e, flops, flops1
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, nnn1, ierr
      real(kind=8) Memory, flop,n2,n1
      DT ctemp
      DT,target,allocatable:: alldat_loc_in(:)
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
      integer*8 ::nnz_loc
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), select_row_pre(:), select_col_pre(:)
      integer::mrange_dummy(1), nrange_dummy(1)
      type(intersect), allocatable :: submats(:)
      type(intersect) :: submats_dummy(1)


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
      blocks%level_half = BF_Switchlevel(level_butterfly, option%pat_comp)
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

            allocate(submats(nr*nc))
            do index_ij = 1, nr*nc
               submats(index_ij)%nr=0
               submats(index_ij)%nc=0
            enddo
            n1 = MPI_Wtime()

            nnz_loc=0
            do index_ij = 1, nr*nc
               index_i_loc = (index_ij - 1)/nc + 1
               index_j_loc = mod(index_ij - 1, nc) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               call BF_compress_NlogN_oneblock_R_sample(submats,blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij,level, nnz_loc, flops1)
            enddo
            n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2 - n1
            if(Nboundall==0)then ! Nboundall>0 means there are intersections with masks, which cannot use contiguous buffers yet.
               allocate(alldat_loc_in(nnz_loc))
               call LogMemory(stats, SIZEOF(alldat_loc_in)/1024.0d3)
               call element_Zmn_blocklist_user(submats, nr*nc, msh, option, ker, 0, passflag, ptree, stats, alldat_loc_in)
            else
               call element_Zmn_blocklist_user(submats, nr*nc, msh, option, ker, 0, passflag, ptree, stats)
            endif

            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp/max(1,blocks%level_butterfly/2)
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,index_i,index_j,index_i_loc,index_j_loc,rank_new1,flops1) reduction(MAX:rank_new) reduction(+:flops)
#endif
            do index_ij = 1, nr*nc
               index_i_loc = (index_ij - 1)/nc + 1
               index_j_loc = mod(index_ij - 1, nc) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               call BF_compress_NlogN_oneblock_R_rankreveal(submats,blocks, option, stats, msh, ker, ptree, index_i, index_j, index_ij,level, rank_new1, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp*max(1,blocks%level_butterfly/2)

            do index_ij = 1, nr*nc
               if(Nboundall>0)then
                  call LogMemory(stats, -SIZEOF(submats(index_ij)%dat)/1024.0d3)
                  if(associated(submats(index_ij)%dat))deallocate(submats(index_ij)%dat)
               endif
               if(allocated(submats(index_ij)%rows))deallocate(submats(index_ij)%rows)
               if(allocated(submats(index_ij)%cols))deallocate(submats(index_ij)%cols)
               if(allocated(submats(index_ij)%masks))then
                  call LogMemory(stats, -SIZEOF(submats(index_ij)%masks)/1024.0d3)
                  deallocate(submats(index_ij)%masks)
               endif
            enddo
            deallocate(submats)
            if(allocated(alldat_loc_in))then
               call LogMemory(stats, -SIZEOF(alldat_loc_in)/1024.0d3)
               deallocate(alldat_loc_in)
            endif

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
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
            endif
            endif
            stats%Flop_Fill = stats%Flop_Fill + flops
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

            n1 = MPI_Wtime()
            allocate(submats(nr*nc))
            do index_ij = 1, nr*nc
               submats(index_ij)%nr=0
               submats(index_ij)%nc=0
            enddo
            nnz_loc=0
            do index_ij = 1, nr*nc
               index_j_loc = (index_ij - 1)/nr + 1
               index_i_loc = mod(index_ij - 1, nr) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               call BF_compress_NlogN_oneblock_C_sample(submats,blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, level_final, nnz_loc)
            enddo
            ! !$omp end parallel do
            n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2 - n1
            if(Nboundall==0)then ! Nboundall>0 means there are intersections with masks, which cannot use contiguous buffers yet.
               allocate(alldat_loc_in(nnz_loc))
               call LogMemory(stats, SIZEOF(alldat_loc_in)/1024.0d3)
               call element_Zmn_blocklist_user(submats, nr*nc, msh, option, ker, 0, passflag, ptree, stats, alldat_loc_in)
            else
               call element_Zmn_blocklist_user(submats, nr*nc, msh, option, ker, 0, passflag, ptree, stats)
            endif

            rank_new = 0
            flops = 0
            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp/max(1,blocks%level_butterfly/2)
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,index_i,index_j,index_j_loc,index_i_loc,rank_new1,flops1) reduction(MAX:rank_new) reduction(+:flops)
#endif
            do index_ij = 1, nr*nc
               index_j_loc = (index_ij - 1)/nr + 1
               index_i_loc = mod(index_ij - 1, nr) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c

               call BF_compress_NlogN_oneblock_C_rankreveal(submats, blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, level_final, rank_new1, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops = flops+flops1

            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp*max(1,blocks%level_butterfly/2)
            do index_ij = 1, nr*nc
               if(Nboundall>0)then
                  if(associated(submats(index_ij)%dat))then
                     call LogMemory(stats, -SIZEOF(submats(index_ij)%dat)/1024.0d3)
                     deallocate(submats(index_ij)%dat)
                  endif
               endif
               if(allocated(submats(index_ij)%rows))deallocate(submats(index_ij)%rows)
               if(allocated(submats(index_ij)%cols))deallocate(submats(index_ij)%cols)
               if(allocated(submats(index_ij)%masks))deallocate(submats(index_ij)%masks)
            enddo
            deallocate(submats)
            if(allocated(alldat_loc_in))then
               call LogMemory(stats, -SIZEOF(alldat_loc_in)/1024.0d3)
               deallocate(alldat_loc_in)
            endif

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
            enddo
            if (level > level_final) then
               call BF_exchange_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'C', 'B')
            endif

            call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (level > 0) then
            if (rank_new > rankmax_for_butterfly(level - 1)) then
               rankmax_for_butterfly(level - 1) = rank_new
            endif
            endif
            stats%Flop_Fill = stats%Flop_Fill + flops
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


   subroutine BF_compress_NlogN_oneblock_R_sample(submats, blocks, boundary_map, Nboundall,Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, nnz_loc, flops)


      implicit none

      integer Nboundall,Ninadmissible
      integer*8 nnz_loc
      integer boundary_map(:,:)
      integer groupm_start
      type(intersect) :: submats(:)
      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, jjj, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_r1, rankmax_c, rankmax_min
      integer group_m, group_n, group_m_mid, group_n_mid, idxstart, idxend, mm, nn, index_i, index_ij, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2, inter
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
      integer Nlayer, passflag, levelm, nrow, ncol,rank_new

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      real(kind=8)::n2, n1,overrate

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

      overrate=1
      if(level<=levelm)then
      group_n_mid = group_n
      do i=1,levelm-level
         group_n_mid = floor_safe(group_n_mid/2d0)
      enddo
      idxstart = msh%basis_group(group_n_mid)%head
      idxend = msh%basis_group(group_n_mid)%tail
      inter = min(msh%basis_group(group_n_mid)%tail,msh%basis_group(group_m)%tail)-max(msh%basis_group(group_n_mid)%head,msh%basis_group(group_m)%head)+1
      if(inter>0)then
         if(inter==mm)then
         overrate = 1
         else
         overrate = dble(mm)/dble(mm-inter)
         endif
      endif
      endif
      ! select skeletons here, selection of at most (option%sample_para+option%knn)*nn columns, the first option%sample_para*nn are random, the next option%knn*nn are nearest points
      rankmax_r1 = min(ceiling_safe(option%sample_para*nn*overrate), mm)
      if (level == 0) rankmax_r1 = min(ceiling_safe(option%sample_para_outer*nn*overrate), mm)
      rankmax_c = nn
      allocate (select_row(rankmax_r1 + nn*option%knn))

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

      call remove_dup_int(select_row, rankmax_r1, rankmax_r)

      if (level == 0) then

         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head

         if (Nboundall > 0) then
            allocate (submats(index_ij)%dat(rankmax_r, nn))
            submats(index_ij)%dat=0
            call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
         endif
         nnz_loc = nnz_loc + rankmax_r*nn

         allocate (submats(index_ij)%rows(rankmax_r))
         allocate (submats(index_ij)%cols(nn))
         submats(index_ij)%nr = rankmax_r
         submats(index_ij)%nc = nn
         do i = 1, rankmax_r
            submats(index_ij)%rows(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, nn
            submats(index_ij)%cols(j) = header_n + j - 1
         enddo

         if (Nboundall > 0) then
            allocate (submats(index_ij)%masks(rankmax_r, nn))
            call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
            submats(index_ij)%masks = 1
            do i = 1, rankmax_r
               group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
               do jj=1,Ninadmissible
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, nn
                        if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                     enddo
                  endif
               enddo
            enddo
         endif


      elseif (level == level_butterfly + 1) then
         index_i_loc_s = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
         rank_new = size(blocks%ButterflySkel(level - 1)%inds(index_i_loc_s, 1)%array)
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head
         index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1

         ! allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
         if (Nboundall > 0) then
            allocate (submats(index_ij)%dat(mm, rank_new))
            submats(index_ij)%dat=0
            call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
         endif
         nnz_loc = nnz_loc + mm*rank_new

         allocate (submats(index_ij)%rows(mm))
         allocate (submats(index_ij)%cols(rank_new))
         submats(index_ij)%nr = mm
         submats(index_ij)%nc = rank_new
         do i = 1, mm
            submats(index_ij)%rows(i) = i + header_m - 1
         enddo
         do j = 1, rank_new
            submats(index_ij)%cols(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s, 1)%array(j) + header_n - 1
         enddo

         if (Nboundall > 0) then
            allocate (submats(index_ij)%masks(mm, rank_new))
            call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
            submats(index_ij)%masks = 1
            do i = 1, mm
               group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
               do jj=1,Ninadmissible
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, rank_new
                        if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                     enddo
                  endif
               enddo
            enddo
         endif

      else
         index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
         index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
         index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1
         index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1
         header_m = msh%basis_group(group_m)%head
         header_n1 = msh%basis_group(group_n)%head
         header_n2 = msh%basis_group(2*group_n + 1)%head
         nnn1 = msh%basis_group(2*group_n)%tail - msh%basis_group(2*group_n)%head + 1

         if (Nboundall > 0) then
            allocate (submats(index_ij)%dat(rankmax_r, nn))
            submats(index_ij)%dat=0
            call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
         endif
         nnz_loc = nnz_loc + rankmax_r*nn

         allocate (submats(index_ij)%rows(rankmax_r))
         allocate (submats(index_ij)%cols(nn))
         submats(index_ij)%nr = rankmax_r
         submats(index_ij)%nc = nn
         do i = 1, rankmax_r
            submats(index_ij)%rows(i) = select_row(i) + header_m - 1
         enddo
         do j = 1, nn
            if (j <= nn1) then
               submats(index_ij)%cols(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc)%array(j) + header_n1 - 1
            else
               submats(index_ij)%cols(j) = blocks%ButterflySkel(level - 1)%inds(index_ii_loc, index_jj_loc + 1)%array(j - nn1) + header_n2 - 1
            endif
         enddo

         if (Nboundall > 0) then

            allocate (submats(index_ij)%masks(rankmax_r, nn))
            call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
            submats(index_ij)%masks = 1

            do i = 1, rankmax_r
               group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
               do jj=1,Ninadmissible
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, nn
                        if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                     enddo
                  endif
               enddo
            enddo
         endif

      endif

      if (allocated(core)) deallocate (core)
      if (allocated(core_tmp)) deallocate (core_tmp)
      if (allocated(tau)) deallocate (tau)
      if (allocated(jpvt)) deallocate (jpvt)
      if (allocated(matrix_V)) deallocate (matrix_V)
      if (allocated(select_row)) deallocate (select_row)
      if (allocated(select_column)) deallocate (select_column)

   end subroutine BF_compress_NlogN_oneblock_R_sample


   subroutine BF_MD_compress_N_oneblock_R_sample(Ndim, dim_MD, subtensors, blocks, boundary_map, Nboundall,Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_ij, bb_m, level, flops)


      implicit none

      integer Nboundall,Ninadmissible, Ndim, bb_m, dim,dim_i,flag,index_ii_scalar,dim_ii
      integer dim_MD(Ndim+2),idx_MD(Ndim+2), dim_MD1(Ndim*2-1),idx_MD1(Ndim*2-1), idx_c_m(Ndim), dims_m(Ndim), dims_row(Ndim),dims_col(Ndim),idx_m(Ndim),idx_n(Ndim), dims_bm(Ndim), group_scalar, dim_product, Nsample
      integer boundary_map(:,:,:)
      integer groupm_start(Ndim)
      type(intersect_MD) :: subtensors(:)
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      integer i, j,j1, level_butterfly, header_m(Ndim), header_n(Ndim), sampleidx(Ndim*2), sampleidx1(Ndim*2)
      integer group_m(Ndim), group_n(Ndim), group_m_mid(Ndim), group_n_mid(Ndim), idxstart(Ndim), idxend(Ndim), mm(Ndim), nn(Ndim), nn_scalar, mmnn(Ndim*2), index_i(Ndim), index_ij, index_ij1, index_j, jj
      integer level
      integer header_n1, header_n2, nn1, nn2, index_ii(Ndim), index_jj, edge_n,jjj
      real(kind=8) flops
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row_tmp(:), select_column(:), column_pivot(:), row_pivot(:)
      type(iarray),allocatable::select_idx(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:), order(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:),p(:)
      integer Nlayer, passflag, levelm, nrow, ncol,rank_new, Nnear_dim(Ndim), Nmoresample

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:), tmpidx(:)
      real(kind=8)::n2, n1,overrate

      flops = 0
      level_butterfly = blocks%level_butterfly
      levelm = floor_safe(dble(level_butterfly)/2d0)

      call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
      index_i = idx_MD(1:Ndim)
      dim = idx_MD(Ndim+2)
      index_j = idx_MD(Ndim+1)
      ! dims_m = 2**(level_butterfly-levelm)
      call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb_m, idx_c_m)


      group_m = blocks%row_group
      group_n = blocks%col_group
      group_m = group_m*2**level - 1 + index_i
      group_n = group_n*2**(level_butterfly-levelm) + blocks%idx_c_m - 1 + idx_c_m -1
      group_n(dim) = group_n(dim)*2**(levelm - level) + index_j - 1
      do dim_i=1,Ndim
         mm(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
      enddo
      do dim_i=1,Ndim
         nn(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
      enddo

      if (level == 0) then
         nn_scalar = nn(dim)
      elseif (level == level_butterfly + 1) then
         write(*,*)"should not come to level == level_butterfly + 1 in BF_MD_compress_N_oneblock_R_sample"
      else
         do dim_i=1,Ndim
         index_ii(dim_i) = int((index_i(dim_i) + 1)/2)
         enddo
         index_jj = 2*index_j - 1
         dims_row = 2**(level-1)
         call MultiIndexToSingleIndex(Ndim, dims_row, index_ii_scalar, index_ii)

         nn1 = size(blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj, dim)%array, 1)
         nn2 = size(blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj + 1, dim)%array, 1)
         nn_scalar = nn1 + nn2
      endif
      mmnn(1:Ndim)=mm
      mmnn(1+Ndim:Ndim*2)=nn


      do dim_i =1,Ndim
         header_m(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%head
         header_n(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%head
         if (level > 0 .and. dim==dim_i) then
            header_n1 = header_n(dim_i)
            header_n2 = msh(dim_i)%basis_group(2*group_n(dim_i)+1)%head
         endif
      enddo

      overrate=1
      ! select skeletons here, selection of at most (option%sample_para+option%knn)*nn columns, the first option%sample_para*nn are random, the next option%knn*nn are nearest points
      allocate (select_idx(Ndim*2))
      do dim_i=1,Ndim*2
         if(option%fastsample_tensor==2 .and. dim_i>Ndim .and. dim_i/=Ndim+dim)then
            sampleidx(dim_i) =2
            allocate(select_idx(dim_i)%dat(sampleidx(dim_i)))
            select_idx(dim_i)%dat(1)=1
            select_idx(dim_i)%dat(2)=mmnn(dim_i)
         else
            sampleidx(dim_i) = min(ceiling_safe(option%sample_para*nn_scalar*overrate), mmnn(dim_i))
            if (level == 0) sampleidx(dim_i) = min(ceiling_safe(option%sample_para_outer*nn_scalar*overrate), mmnn(dim_i))
            allocate(select_idx(dim_i)%dat(sampleidx(dim_i)))
            call linspaceI(1, mmnn(dim_i), sampleidx(dim_i), select_idx(dim_i)%dat(1:sampleidx(dim_i)))
         endif
      enddo

      if(option%knn>0)then
         Nnear_dim=0
         do dim_i=1,Ndim ! first loop the Ndim-1 dimensions on the column dimension
            do j = 1, nn(dim_i)
               edge_n = header_n(dim_i) + j - 1
               do jjj = 1, option%knn
               if (msh(dim_i)%nns(edge_n, jjj) >= msh(dim_i)%basis_group(group_m(dim_i))%head .and. msh(dim_i)%nns(edge_n, jjj) <= msh(dim_i)%basis_group(group_m(dim_i))%tail) then
                  Nnear_dim(dim_i) = 1
               endif
               enddo
            enddo
         enddo
         if(sum(Nnear_dim)==Ndim)then ! group_m and group_n are closeby
            do dim_i=1,Ndim
               if(dim_i/=dim)then
                  Nmoresample=0
                  do j = 1, nn(dim_i)
                     edge_n = header_n(dim_i) + j - 1
                     do jjj = 1, option%knn
                     if (msh(dim_i)%nns(edge_n, jjj) >= msh(dim_i)%basis_group(group_m(dim_i))%head .and. msh(dim_i)%nns(edge_n, jjj) <= msh(dim_i)%basis_group(group_m(dim_i))%tail) then
                        Nmoresample = Nmoresample + 1
                        exit
                     endif
                     enddo
                  enddo
                  if(Nmoresample>0)then
                     sampleidx1(dim_i+Ndim) = sampleidx(dim_i+Ndim) + Nmoresample
                     allocate(tmpidx(sampleidx(dim_i+Ndim)))
                     tmpidx=select_idx(dim_i+Ndim)%dat
                     deallocate(select_idx(dim_i+Ndim)%dat)
                     allocate(select_idx(dim_i+Ndim)%dat(sampleidx1(dim_i+Ndim)))
                     select_idx(dim_i+Ndim)%dat(1:sampleidx(dim_i+Ndim))=tmpidx
                     deallocate(tmpidx)
                     Nmoresample=0

                     do j = 1, nn(dim_i)
                        edge_n = header_n(dim_i) + j - 1
                        do jjj = 1, option%knn
                        if (msh(dim_i)%nns(edge_n, jjj) >= msh(dim_i)%basis_group(group_m(dim_i))%head .and. msh(dim_i)%nns(edge_n, jjj) <= msh(dim_i)%basis_group(group_m(dim_i))%tail) then
                           Nmoresample = Nmoresample + 1
                           select_idx(dim_i+Ndim)%dat(sampleidx(dim_i+Ndim)+Nmoresample)=j
                           exit
                        endif
                        enddo
                     enddo
                     call remove_dup_int(select_idx(dim_i+Ndim)%dat, sampleidx1(dim_i+Ndim), sampleidx(dim_i+Ndim))
                  endif
               endif
            enddo

            do dim_i=1,Ndim ! then loop the Ndim dimensions on the row dimension
               Nmoresample=0
               do j = 1, nn(dim_i)
                  edge_n = header_n(dim_i) + j - 1
                  do jjj = 1, option%knn
                  if (msh(dim_i)%nns(edge_n, jjj) >= msh(dim_i)%basis_group(group_m(dim_i))%head .and. msh(dim_i)%nns(edge_n, jjj) <= msh(dim_i)%basis_group(group_m(dim_i))%tail) then
                     Nmoresample = Nmoresample + 1
                  endif
                  enddo
               enddo

               if(Nmoresample>0)then
                  sampleidx1(dim_i) = sampleidx(dim_i) + Nmoresample
                  allocate(tmpidx(sampleidx(dim_i)))
                  tmpidx=select_idx(dim_i)%dat
                  deallocate(select_idx(dim_i)%dat)
                  allocate(select_idx(dim_i)%dat(sampleidx1(dim_i)))
                  select_idx(dim_i)%dat(1:sampleidx(dim_i))=tmpidx
                  deallocate(tmpidx)
                  Nmoresample=0
                  do j = 1, nn(dim_i)
                     edge_n = header_n(dim_i) + j - 1
                     do jjj = 1, option%knn
                     if (msh(dim_i)%nns(edge_n, jjj) >= msh(dim_i)%basis_group(group_m(dim_i))%head .and. msh(dim_i)%nns(edge_n, jjj) <= msh(dim_i)%basis_group(group_m(dim_i))%tail) then
                        Nmoresample = Nmoresample + 1
                        select_idx(dim_i)%dat(sampleidx(dim_i)+Nmoresample)=msh(dim_i)%nns(edge_n, jjj) + 1 - header_m(dim_i)
                     endif
                     enddo
                  enddo
                  call remove_dup_int(select_idx(dim_i)%dat, sampleidx1(dim_i), sampleidx(dim_i))
               endif
            enddo
         endif
      endif

      dims_row = sampleidx(1:Ndim)
      dims_col = sampleidx(1+Ndim:2*Ndim)
      dims_col(dim) = nn_scalar




      allocate (subtensors(index_ij)%dat(product(dims_row),product(dims_col)))
      subtensors(index_ij)%dat=0
      call LogMemory(stats, SIZEOF(subtensors(index_ij)%dat)/1024.0d3)
      allocate (subtensors(index_ij)%masks(product(dims_row),product(dims_col)))
      call LogMemory(stats, SIZEOF(subtensors(index_ij)%masks)/1024.0d3)
      subtensors(index_ij)%masks=1


      if(option%fastsample_tensor==1)then ! uniform sampling the row dimension of the unfolded tensor
         subtensors(index_ij)%masks=0
         dim_product = product(dims_row)
         do dim_i =1,Ndim
            if (dim/=dim_i) then
               dim_product = dim_product*dims_col(dim_i)
            endif
         enddo
         Nsample = min(dim_product,ceiling_safe(nn_scalar*option%sample_para*Ndim))
         dim_MD1(1:Ndim) = dims_row
         dim_ii=0
         do dim_i=1,Ndim
            if(dim_i/=dim)then
               dim_ii=dim_ii+1
               dim_MD1(Ndim+dim_ii) = dims_col(dim_i)
            endif
         enddo
         allocate(p(dim_product))
         call rperm(dim_product, p)
         do index_ij1=1,Nsample
            call SingleIndexToMultiIndex(Ndim*2-1,dim_MD1, p(index_ij1), idx_MD1)
            idx_m = idx_MD1(1:Ndim)
            dim_ii=0
            do dim_i=1,Ndim
               if(dim_i/=dim)then
                  dim_ii=dim_ii+1
                  idx_n(dim_i) = idx_MD1(Ndim+dim_ii)
               endif
            enddo
            do j1=1,dims_col(dim)
               idx_n(dim)=j1
               call MultiIndexToSingleIndex(Ndim,dims_row, i, idx_m)
               call MultiIndexToSingleIndex(Ndim,dims_col, j, idx_n)
               subtensors(index_ij)%masks(i,j)=1
            enddo
         enddo
         deallocate(p)
      endif

      subtensors(index_ij)%nr = dims_row
      subtensors(index_ij)%nc = dims_col
      allocate (subtensors(index_ij)%rows(Ndim))
      do dim_i=1,Ndim
         allocate (subtensors(index_ij)%rows(dim_i)%dat(dims_row(dim_i)))
         do i = 1, dims_row(dim_i)
            subtensors(index_ij)%rows(dim_i)%dat(i) = header_m(dim_i) + select_idx(dim_i)%dat(i) - 1
         enddo
      enddo

      allocate (subtensors(index_ij)%cols(Ndim))
      do dim_i=1,Ndim
         allocate (subtensors(index_ij)%cols(dim_i)%dat(dims_col(dim_i)))
         do j = 1, dims_col(dim_i)
            if(dim==dim_i)then
               if (level == 0) then
                  subtensors(index_ij)%cols(dim_i)%dat(j) = header_n(dim_i) + j - 1
               elseif (level == level_butterfly + 1) then
                  write(*,*)"should not come here in BF_MD_compress_N_oneblock_R_sample"
               else
                  if (j <= nn1) then
                     subtensors(index_ij)%cols(dim_i)%dat(j) = blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj, dim)%array(j) + header_n1 - 1
                  else
                     subtensors(index_ij)%cols(dim_i)%dat(j) = blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj + 1, dim)%array(j - nn1) + header_n2 - 1
                  endif
               endif
            else
               subtensors(index_ij)%cols(dim_i)%dat(j) = header_n(dim_i) + select_idx(dim_i+Ndim)%dat(j) - 1
            endif
         enddo
      enddo


      if (Nboundall > 0) then
         allocate (subtensors(index_ij)%masks(product(dims_row),product(dims_col)))
         call LogMemory(stats, SIZEOF(subtensors(index_ij)%masks)/1024.0d3)
         subtensors(index_ij)%masks = 1
         do i = 1, product(dims_row)
            call SingleIndexToMultiIndex(Ndim,dims_row, i, idx_m)
            do dim_i = 1,Ndim
               group_m_mid(dim_i) = findgroup(subtensors(index_ij)%rows(dim_i)%dat(idx_m(dim_i)), msh(dim_i), levelm, blocks%row_group(dim_i))
            enddo

            dims_bm= Nboundall
            call MultiIndexToSingleIndex(Ndim,dims_bm, group_scalar, group_m_mid - groupm_start + 1)
            do jj=1,Ninadmissible
               group_n_mid = boundary_map(group_scalar,jj,:)
               if (ALL(group_n_mid /= -1)) then
                  do dim_i=1,Ndim
                     idxstart(dim_i) = msh(dim_i)%basis_group(group_n_mid(dim_i))%head
                     idxend(dim_i) = msh(dim_i)%basis_group(group_n_mid(dim_i))%tail
                  enddo
                  do j = 1, product(dims_col)
                     call SingleIndexToMultiIndex(Ndim,dims_col, j, idx_n)
                     flag=Ndim
                     do dim_i=1,Ndim
                        if (subtensors(index_ij)%cols(dim_i)%dat(idx_n(dim_i)) >= idxstart(dim_i) .and. subtensors(index_ij)%cols(dim_i)%dat(idx_n(dim_i)) <= idxend(dim_i))flag = flag-1
                     enddo
                     if(flag==0)subtensors(index_ij)%masks(i, j) = 0
                  enddo
               endif
            enddo
         enddo
      endif


      do dim_i=1,Ndim*2
         deallocate(select_idx(dim_i)%dat)
      enddo
      deallocate(select_idx)

   end subroutine BF_MD_compress_N_oneblock_R_sample




   subroutine BF_compress_NlogN_oneblock_R_rankreveal(submats, blocks, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, rank_new, flops)


      implicit none

      type(intersect) :: submats(:)
      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, jjj, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_r1, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, group_m_mid, group_n_mid, idxstart, idxend, mm, nn, index_i, index_ij, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2, inter
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

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      real(kind=8)::n2, n1,overrate

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
      rankmax_c = nn
      allocate (select_column(rankmax_c))
      do i = 1, rankmax_c
         select_column(i) = i
      enddo

      rankmax_r = size(submats(index_ij)%dat,1)
      if (level == 0) then

         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head

         allocate (core(rankmax_r, rankmax_c))
         do j = 1, rankmax_c
            core(:, j) = submats(index_ij)%dat(:, select_column(j))
         enddo

         allocate (jpvt(max(rankmax_c, rankmax_r)))
         allocate (tau(max(rankmax_c, rankmax_r)))
         jpvt = 0
         call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
         flops = flops + flop

         if (rank_new > 0) then
            call un_or_mqrf90(core, tau, submats(index_ij)%dat, 'L', 'C', rankmax_r, nn, rank_new, flop=flop)
            flops = flops + flop
            call trsmf90(core, submats(index_ij)%dat, 'L', 'U', 'N', 'N', rank_new, nn, flop=flop)
            flops = flops + flop
         else
            rank_new = 1
            jpvt(1) = 1
            submats(index_ij)%dat = 0
         endif

         index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1
         allocate (blocks%ButterflyV%blocks(index_j_loc_k)%matrix(nn, rank_new))
         call copymatT(submats(index_ij)%dat(1:rank_new, 1:nn), blocks%ButterflyV%blocks(index_j_loc_k)%matrix, rank_new, nn)

         index_j_loc_s = (index_j - blocks%ButterflySkel(0)%idx_c)/blocks%ButterflySkel(0)%inc_c + 1
         allocate (blocks%ButterflySkel(0)%inds(1, index_j_loc_s)%array(rank_new))
         do j = 1, rank_new
            blocks%ButterflySkel(0)%inds(1, index_j_loc_s)%array(j) = select_column(jpvt(j))
         enddo

      elseif (level == level_butterfly + 1) then
         index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
         mm = size(submats(index_ij)%dat,1)
         rank_new = size(submats(index_ij)%dat,2)
         allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
         blocks%ButterflyU%blocks(index_i_loc_k)%matrix = submats(index_ij)%dat
      else
         index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
         index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
         index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1
         index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1
         header_m = msh%basis_group(group_m)%head
         header_n1 = msh%basis_group(group_n)%head
         header_n2 = msh%basis_group(2*group_n + 1)%head
         nnn1 = msh%basis_group(2*group_n)%tail - msh%basis_group(2*group_n)%head + 1

         allocate (core(rankmax_r, rankmax_c))
         do j = 1, rankmax_c
            core(:, j) = submats(index_ij)%dat(:, select_column(j))
         enddo

         allocate (jpvt(max(rankmax_c, rankmax_r)))
         allocate (tau(max(rankmax_c, rankmax_r)))
         jpvt = 0
         call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
         flops = flops + flop

         if (rank_new > 0) then
            call un_or_mqrf90(core, tau, submats(index_ij)%dat, 'L', 'C', rankmax_r, nn, rank_new, flop=flop)
            flops = flops + flop
            call trsmf90(core, submats(index_ij)%dat, 'L', 'U', 'N', 'N', rank_new, nn, flop=flop)
            flops = flops + flop
         else
            rank_new = 1
            jpvt(1) = 1
            submats(index_ij)%dat = 0
         endif

         allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(rank_new, nn1))
         allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix(rank_new, nn2))

         blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = submats(index_ij)%dat(1:rank_new, 1:nn1)
         blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix = submats(index_ij)%dat(1:rank_new, 1 + nn1:nn)

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

      endif

      if (allocated(core)) deallocate (core)
      if (allocated(core_tmp)) deallocate (core_tmp)
      if (allocated(tau)) deallocate (tau)
      if (allocated(jpvt)) deallocate (jpvt)
      if (allocated(matrix_V)) deallocate (matrix_V)
      if (allocated(select_row)) deallocate (select_row)
      if (allocated(select_column)) deallocate (select_column)

   end subroutine BF_compress_NlogN_oneblock_R_rankreveal




   subroutine BF_MD_compress_N_oneblock_R_rankreveal(Ndim, dim_MD, subtensors, blocks, option, stats, msh, ker, ptree, index_ij, bb_m, level, rank_new, flops)


      implicit none

      integer Ndim, bb_m, dim,dim_i,dim_ii,flag
      integer dim_MD(Ndim+2),idx_MD(Ndim+2), idx_c_m(Ndim), dims_m(Ndim), dims_row(Ndim),dims_col(Ndim),idx_m(Ndim),idx_n(Ndim),idx_n1(Ndim-1),dims_col1(Ndim-1)
      type(intersect_MD) :: subtensors(:)
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      integer i, j, j1, j_dim, level_butterfly, header_m(Ndim), header_n(Ndim), rankmax_r, rankmax_c
      integer group_m(Ndim), group_n(Ndim), group_m_mid(Ndim), group_n_mid(Ndim), idxstart(Ndim), idxend(Ndim), mm(Ndim), nn(Ndim), nn_scalar, mmnn(Ndim*2), index_i(Ndim), index_i_scalar, index_ij, index_j, index_j_s, index_j_k, jj,nnn1
      integer level, cnt
      integer header_n1, header_n2, nn1, nn2, index_ii(Ndim), index_ii_scalar, index_jj
      real(kind=8) flops,flop
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row_tmp(:), select_column(:), column_pivot(:), row_pivot(:)
      type(iarray),allocatable::select_idx(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:), order(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, passflag, levelm, nrow, ncol,rank_new

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      real(kind=8)::n2, n1,overrate

      flops = 0
      level_butterfly = blocks%level_butterfly
      levelm = floor_safe(dble(level_butterfly)/2d0)

      call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
      index_i = idx_MD(1:Ndim)
      dim = idx_MD(Ndim+2)
      index_j = idx_MD(Ndim+1)
      ! dims_m = 2**(level_butterfly-levelm)
      call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb_m, idx_c_m)


      group_m = blocks%row_group
      group_n = blocks%col_group
      group_m = group_m*2**level - 1 + index_i
      group_n = group_n*2**(level_butterfly-levelm) + blocks%idx_c_m - 1 + idx_c_m -1
      group_n(dim) = group_n(dim)*2**(levelm - level) + index_j - 1
      do dim_i=1,Ndim
         mm(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
      enddo
      do dim_i=1,Ndim
         nn(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
      enddo

      if (level == 0) then
         nn_scalar = nn(dim)
      elseif (level == level_butterfly + 1) then
         write(*,*)"should not come to level == level_butterfly + 1 in BF_MD_compress_N_oneblock_R_rankreveal"
      else
         do dim_i=1,Ndim
         index_ii(dim_i) = int((index_i(dim_i) + 1)/2)
         enddo
         index_jj = 2*index_j - 1
         dims_row = 2**level
         call MultiIndexToSingleIndex(Ndim, dims_row, index_i_scalar, index_i)
         dims_row = 2**(level-1)
         call MultiIndexToSingleIndex(Ndim, dims_row, index_ii_scalar, index_ii)

         nn1 = size(blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj, dim)%array, 1)
         nn2 = size(blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj + 1, dim)%array, 1)
         nn_scalar = nn1 + nn2
      endif
      mmnn(1:Ndim)=mm
      mmnn(1+Ndim:Ndim*2)=nn

      dims_row = subtensors(index_ij)%nr
      dims_col = subtensors(index_ij)%nc

      do dim_i =1,Ndim
         header_m(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%head
         header_n(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%head
         if (level > 0 .and. dim==dim_i) then
            header_n1 = header_n(dim_i)
            header_n2 = msh(dim_i)%basis_group(2*group_n(dim_i)+1)%head
         endif
      enddo


      rankmax_c = dims_col(dim)
      rankmax_r=product(dims_row)
      do dim_i=1,Ndim
         if(dim_i/=dim)rankmax_r = rankmax_r*dims_col(dim_i)
      enddo
      allocate (core(rankmax_r, rankmax_c))
core=0
      allocate (core_tmp(rankmax_r, rankmax_c))

      ! reshaping
do j_dim = 1,dims_col(dim)
         cnt=0
      do j = 1, product(dims_col)
         call SingleIndexToMultiIndex(Ndim, dims_col, j, idx_n)
         if(idx_n(dim)==j_dim)then ! this makes sure that all the other 2*Ndim-1 indices are traversed before moving to the next idx_n(dim) value
            do i=1,product(dims_row)
               if(subtensors(index_ij)%masks(i,j)==1)then
                  cnt = cnt+1
                  core(cnt,idx_n(dim)) = subtensors(index_ij)%dat(i,j)
         endif
         enddo
         endif
         enddo
      enddo
      core_tmp=core

      allocate (jpvt(max(rankmax_c, rankmax_r)))
      allocate (tau(max(rankmax_c, rankmax_r)))
      jpvt = 0
      call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
      flops = flops + flop

      if (rank_new > 0) then
         call un_or_mqrf90(core, tau, core_tmp, 'L', 'C', rankmax_r, rankmax_c, rank_new, flop=flop)
         flops = flops + flop
         call trsmf90(core, core_tmp, 'L', 'U', 'N', 'N', rank_new, rankmax_c, flop=flop)
         flops = flops + flop
      else
         rank_new = 1
         jpvt(1) = 1
         core_tmp = 0
      endif

      if (level == 0) then
         index_j_k=index_j
         index_j_s=index_j
         allocate (blocks%ButterflyV(bb_m)%blocks(index_j_k,dim)%matrix(rankmax_c, rank_new))
         call copymatT(core_tmp(1:rank_new, 1:rankmax_c), blocks%ButterflyV(bb_m)%blocks(index_j_k,dim)%matrix, rank_new, rankmax_c)
         allocate (blocks%ButterflySkel_R(bb_m,level)%inds(1, index_j_s,dim)%array(rank_new))
         do j = 1, rank_new
            blocks%ButterflySkel_R(bb_m,level)%inds(1, index_j_s,dim)%array(j) = jpvt(j)
         enddo

      elseif (level == level_butterfly + 1) then
         write(*,*)"should not arrive at level == level_butterfly + 1 in BF_MD_compress_N_oneblock_R_rankreveal"
      else
         index_j_k=2*index_j-1
         index_j_s=index_j
         nnn1 = msh(dim)%basis_group(2*group_n(dim))%tail - msh(dim)%basis_group(2*group_n(dim))%head + 1
         allocate (blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_scalar, index_j_k,dim)%matrix(rank_new, nn1))
         allocate (blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_scalar, index_j_k+1,dim)%matrix(rank_new, nn2))
         blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_scalar, index_j_k,dim)%matrix = core_tmp(1:rank_new, 1:nn1)
         blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_scalar, index_j_k+1,dim)%matrix = core_tmp(1:rank_new, 1 + nn1:rankmax_c)

         allocate (blocks%ButterflySkel_R(bb_m,level)%inds(index_i_scalar, index_j_s,dim)%array(rank_new))
         do j = 1, rank_new
            if (jpvt(j) <= nn1) then
               blocks%ButterflySkel_R(bb_m,level)%inds(index_i_scalar, index_j_s,dim)%array(j) = blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj,dim)%array(jpvt(j))
            else
               blocks%ButterflySkel_R(bb_m,level)%inds(index_i_scalar, index_j_s,dim)%array(j) = blocks%ButterflySkel_R(bb_m,level - 1)%inds(index_ii_scalar, index_jj + 1,dim)%array(jpvt(j) - nn1) + nnn1
            endif
         enddo

      endif

      deallocate(core)
      deallocate(core_tmp)
      deallocate(jpvt)
      deallocate(tau)

   end subroutine BF_MD_compress_N_oneblock_R_rankreveal


   subroutine BF_MD_compress_N_oneblock_C_sample(Ndim, dim_MD, subtensors, blocks, boundary_map, Nboundall,Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_ij, bb_m, level, flops)


      implicit none

      integer Nboundall,Ninadmissible, Ndim, bb_m, dim,dim_i,flag,index_jj_scalar, dim_ii
      integer dim_MD(Ndim+2),idx_MD(Ndim+2), dim_MD1(Ndim*2-1),idx_MD1(Ndim*2-1), idx_r_m(Ndim), dims_m(Ndim), dims_row(Ndim),dims_col(Ndim),idx_m(Ndim),idx_n(Ndim), dims_bm(Ndim), group_scalar, dim_product, Nsample
      integer boundary_map(:,:,:)
      integer groupm_start(Ndim)
      type(intersect_MD) :: subtensors(:)
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      integer i, j, i1, level_butterfly, header_m(Ndim), header_n(Ndim), sampleidx(Ndim*2), sampleidx1(Ndim*2)
      integer group_m(Ndim), group_n(Ndim), group_m_mid(Ndim), group_n_mid(Ndim), idxstart(Ndim), idxend(Ndim), mm(Ndim), nn(Ndim), mm_scalar, mmnn(Ndim*2), index_j(Ndim), index_ij, index_ij1, index_i, jj
      integer level
      integer header_m1, header_m2, mm1, mm2, index_jj(Ndim), index_ii
      real(kind=8) flops
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row_tmp(:), select_column(:), column_pivot(:), row_pivot(:)
      type(iarray),allocatable::select_idx(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:), order(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:),p(:)
      integer Nlayer, passflag, levelm, nrow, ncol,rank_new, Nnear_dim(Ndim), edge_m, iii, Nmoresample

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:), tmpidx(:)
      real(kind=8)::n2, n1,overrate

      flops = 0
      level_butterfly = blocks%level_butterfly
      levelm = floor_safe(dble(level_butterfly)/2d0)

      call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
      index_i = idx_MD(1)
      dim = idx_MD(Ndim+2)
      index_j = idx_MD(2:Ndim+1)
      ! dims_m = 2**(levelm)
      call SingleIndexToMultiIndex(Ndim, blocks%nr_m, bb_m, idx_r_m)


      group_m = blocks%row_group
      group_n = blocks%col_group
      group_n = group_n*2**(level_butterfly - level + 1) - 1 + index_j

      group_m = group_m*2**(levelm) + blocks%idx_r_m - 1 + idx_r_m -1
      group_m(dim) = group_m(dim)*2**(level - levelm -1) + index_i - 1

      do dim_i=1,Ndim
         mm(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
      enddo
      do dim_i=1,Ndim
         nn(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
      enddo

      if (level == 0) then
         write(*,*)"should not come to level == 0 in BF_MD_compress_N_oneblock_C_sample"
      elseif (level == level_butterfly + 1) then
         mm_scalar = mm(dim)
      else
         do dim_i=1,Ndim
         index_jj(dim_i) = int((index_j(dim_i) + 1)/2)
         enddo
         index_ii = 2*index_i - 1
         dims_col = 2**(level_butterfly - level)
         call MultiIndexToSingleIndex(Ndim, dims_col, index_jj_scalar, index_jj)

         mm1 = size(blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii, index_jj_scalar, dim)%array, 1)
         mm2 = size(blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii+1, index_jj_scalar, dim)%array, 1)
         mm_scalar = mm1 + mm2
      endif
      mmnn(1:Ndim)=mm
      mmnn(1+Ndim:Ndim*2)=nn

      do dim_i =1,Ndim
         header_m(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%head
         header_n(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%head
         if (level < level_butterfly+1 .and. dim==dim_i) then
            header_m1 = header_m(dim_i)
            header_m2 = msh(dim_i)%basis_group(2*group_m(dim_i)+1)%head
         endif
      enddo

      overrate=1
      ! select skeletons here, selection of at most (option%sample_para+option%knn)*nn columns, the first option%sample_para*nn are random, the next option%knn*nn are nearest points
      allocate (select_idx(Ndim*2))
      do dim_i=1,Ndim*2

         if(option%fastsample_tensor==2 .and. dim_i<=Ndim .and. dim_i/=dim)then
            sampleidx(dim_i) =2
            allocate(select_idx(dim_i)%dat(sampleidx(dim_i)))
            select_idx(dim_i)%dat(1)=1
            select_idx(dim_i)%dat(2)=mmnn(dim_i)
         else
            sampleidx(dim_i) = min(ceiling_safe(option%sample_para*mm_scalar*overrate), mmnn(dim_i))
            if (level == level_butterfly+1) sampleidx(dim_i) = min(ceiling_safe(option%sample_para_outer*mm_scalar*overrate), mmnn(dim_i))
            allocate(select_idx(dim_i)%dat(sampleidx(dim_i)))
            call linspaceI(1, mmnn(dim_i), sampleidx(dim_i), select_idx(dim_i)%dat(1:sampleidx(dim_i)))
         endif
      enddo


      if(option%knn>0)then
         Nnear_dim=0
         do dim_i=1,Ndim ! first loop the Ndim-1 dimensions on the row dimension
            do i = 1, mm(dim_i)
               edge_m = header_m(dim_i) + i - 1
               do iii = 1, option%knn
               if (msh(dim_i)%nns(edge_m, iii) >= msh(dim_i)%basis_group(group_n(dim_i))%head .and. msh(dim_i)%nns(edge_m, iii) <= msh(dim_i)%basis_group(group_n(dim_i))%tail) then
                  Nnear_dim(dim_i) = 1
               endif
               enddo
            enddo
         enddo
         if(sum(Nnear_dim)==Ndim)then ! group_m and group_n are closeby
            do dim_i=1,Ndim
               if(dim_i/=dim)then
                  Nmoresample=0
                  do i = 1, mm(dim_i)
                     edge_m = header_m(dim_i) + i - 1
                     do iii = 1, option%knn
                     if (msh(dim_i)%nns(edge_m, iii) >= msh(dim_i)%basis_group(group_n(dim_i))%head .and. msh(dim_i)%nns(edge_m, iii) <= msh(dim_i)%basis_group(group_n(dim_i))%tail) then
                        Nmoresample = Nmoresample + 1
                        exit
                     endif
                     enddo
                  enddo
                  if(Nmoresample>0)then
                     sampleidx1(dim_i) = sampleidx(dim_i) + Nmoresample
                     allocate(tmpidx(sampleidx(dim_i)))
                     tmpidx=select_idx(dim_i)%dat
                     deallocate(select_idx(dim_i)%dat)
                     allocate(select_idx(dim_i)%dat(sampleidx1(dim_i)))
                     select_idx(dim_i)%dat(1:sampleidx(dim_i))=tmpidx
                     deallocate(tmpidx)
                     Nmoresample=0

                     do i = 1, mm(dim_i)
                        edge_m = header_m(dim_i) + i - 1
                        do iii = 1, option%knn
                        if (msh(dim_i)%nns(edge_m, iii) >= msh(dim_i)%basis_group(group_n(dim_i))%head .and. msh(dim_i)%nns(edge_m, iii) <= msh(dim_i)%basis_group(group_n(dim_i))%tail) then
                           Nmoresample = Nmoresample + 1
                           select_idx(dim_i)%dat(sampleidx(dim_i)+Nmoresample)=i
                           exit
                        endif
                        enddo
                     enddo
                     call remove_dup_int(select_idx(dim_i)%dat, sampleidx1(dim_i), sampleidx(dim_i))
                  endif
               endif
            enddo

            do dim_i=1,Ndim ! then loop the Ndim dimensions on the column dimension
               Nmoresample=0
               do i = 1, mm(dim_i)
                  edge_m = header_m(dim_i) + i - 1
                  do iii = 1, option%knn
                  if (msh(dim_i)%nns(edge_m, iii) >= msh(dim_i)%basis_group(group_n(dim_i))%head .and. msh(dim_i)%nns(edge_m, iii) <= msh(dim_i)%basis_group(group_n(dim_i))%tail) then
                     Nmoresample = Nmoresample + 1
                  endif
                  enddo
               enddo

               if(Nmoresample>0)then
                  sampleidx1(dim_i+Ndim) = sampleidx(dim_i+Ndim) + Nmoresample
                  allocate(tmpidx(sampleidx(dim_i+Ndim)))
                  tmpidx=select_idx(dim_i+Ndim)%dat
                  deallocate(select_idx(dim_i+Ndim)%dat)
                  allocate(select_idx(dim_i+Ndim)%dat(sampleidx1(dim_i+Ndim)))
                  select_idx(dim_i+Ndim)%dat(1:sampleidx(dim_i+Ndim))=tmpidx
                  deallocate(tmpidx)
                  Nmoresample=0
                  do i = 1, mm(dim_i)
                     edge_m = header_m(dim_i) + i - 1
                     do iii = 1, option%knn
                     if (msh(dim_i)%nns(edge_m, iii) >= msh(dim_i)%basis_group(group_n(dim_i))%head .and. msh(dim_i)%nns(edge_m, iii) <= msh(dim_i)%basis_group(group_n(dim_i))%tail) then
                        Nmoresample = Nmoresample + 1
                        select_idx(dim_i+Ndim)%dat(sampleidx(dim_i+Ndim)+Nmoresample)=msh(dim_i)%nns(edge_m, iii) + 1 - header_n(dim_i)
                     endif
                     enddo
                  enddo
                  call remove_dup_int(select_idx(dim_i+Ndim)%dat, sampleidx1(dim_i+Ndim), sampleidx(dim_i+Ndim))
               endif
            enddo
         endif
      endif


      dims_row = sampleidx(1:Ndim)
      dims_col = sampleidx(1+Ndim:2*Ndim)
      dims_row(dim) = mm_scalar




      allocate (subtensors(index_ij)%dat(product(dims_row),product(dims_col)))
      subtensors(index_ij)%dat=0
      call LogMemory(stats, SIZEOF(subtensors(index_ij)%dat)/1024.0d3)
      allocate (subtensors(index_ij)%masks(product(dims_row),product(dims_col)))
      call LogMemory(stats, SIZEOF(subtensors(index_ij)%masks)/1024.0d3)
      subtensors(index_ij)%masks=1

      if(option%fastsample_tensor==1)then ! uniform sampling the row dimension of the unfolded tensor
         subtensors(index_ij)%masks=0
         dim_product = product(dims_col)
         do dim_i =1,Ndim
            if (dim/=dim_i) then
               dim_product = dim_product*dims_row(dim_i)
            endif
         enddo
         Nsample = min(dim_product,ceiling_safe(mm_scalar*option%sample_para*Ndim))

         dim_ii=0
         do dim_i=1,Ndim
            if(dim_i/=dim)then
               dim_ii=dim_ii+1
               dim_MD1(dim_ii) = dims_row(dim_i)
            endif
         enddo
         dim_MD1(Ndim:2*Ndim-1) = dims_col

         allocate(p(dim_product))
         call rperm(dim_product, p)
         do index_ij1=1,Nsample
            call SingleIndexToMultiIndex(Ndim*2-1,dim_MD1, p(index_ij1), idx_MD1)
            dim_ii=0
            do dim_i=1,Ndim
               if(dim_i/=dim)then
                  dim_ii=dim_ii+1
                  idx_m(dim_i) = idx_MD1(dim_ii)
               endif
            enddo
            idx_n = idx_MD1(Ndim:Ndim*2-1)

            do i1=1,dims_row(dim)
               idx_m(dim)=i1
               call MultiIndexToSingleIndex(Ndim,dims_row, i, idx_m)
               call MultiIndexToSingleIndex(Ndim,dims_col, j, idx_n)
               subtensors(index_ij)%masks(i,j)=1
            enddo
         enddo
         deallocate(p)
      endif


      subtensors(index_ij)%nr = dims_row
      subtensors(index_ij)%nc = dims_col
      allocate (subtensors(index_ij)%rows(Ndim))
      do dim_i=1,Ndim
         allocate (subtensors(index_ij)%rows(dim_i)%dat(dims_row(dim_i)))
         do i = 1, dims_row(dim_i)
            if(dim==dim_i)then
               if (level == 0) then
                    write(*,*)"should not come here in BF_MD_compress_N_oneblock_C_sample"
               elseif (level == level_butterfly + 1) then
                  subtensors(index_ij)%rows(dim_i)%dat(i) = header_m(dim_i) + i - 1
               else
                  if (i <= mm1) then
                     subtensors(index_ij)%rows(dim_i)%dat(i) = blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii, index_jj_scalar, dim)%array(i) + header_m1 - 1
                  else
                     subtensors(index_ij)%rows(dim_i)%dat(i) = blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii+1, index_jj_scalar, dim)%array(i - mm1) + header_m2 - 1
                  endif
               endif
            else
                subtensors(index_ij)%rows(dim_i)%dat(i) = header_m(dim_i) + select_idx(dim_i)%dat(i) - 1
            endif
         enddo
      enddo

      allocate (subtensors(index_ij)%cols(Ndim))
      do dim_i=1,Ndim
         allocate (subtensors(index_ij)%cols(dim_i)%dat(dims_col(dim_i)))
         do j = 1, dims_col(dim_i)
               subtensors(index_ij)%cols(dim_i)%dat(j) = header_n(dim_i) + select_idx(dim_i+Ndim)%dat(j) - 1
         enddo
      enddo


      if (Nboundall > 0) then
         allocate (subtensors(index_ij)%masks(product(dims_row),product(dims_col)))
         call LogMemory(stats, SIZEOF(subtensors(index_ij)%masks)/1024.0d3)
         subtensors(index_ij)%masks = 1
         do i = 1, product(dims_row)
            call SingleIndexToMultiIndex(Ndim,dims_row, i, idx_m)
            do dim_i = 1,Ndim
               group_m_mid(dim_i) = findgroup(subtensors(index_ij)%rows(dim_i)%dat(idx_m(dim_i)), msh(dim_i), levelm, blocks%row_group(dim_i))
            enddo

            dims_bm= Nboundall
            call MultiIndexToSingleIndex(Ndim,dims_bm, group_scalar, group_m_mid - groupm_start + 1)
            do jj=1,Ninadmissible
               group_n_mid = boundary_map(group_scalar,jj,:)
               if (ALL(group_n_mid /= -1)) then
                  do dim_i=1,Ndim
                     idxstart(dim_i) = msh(dim_i)%basis_group(group_n_mid(dim_i))%head
                     idxend(dim_i) = msh(dim_i)%basis_group(group_n_mid(dim_i))%tail
                  enddo
                  do j = 1, product(dims_col)
                     call SingleIndexToMultiIndex(Ndim,dims_col, j, idx_n)
                     flag=Ndim
                     do dim_i=1,Ndim
                        if (subtensors(index_ij)%cols(dim_i)%dat(idx_n(dim_i)) >= idxstart(dim_i) .and. subtensors(index_ij)%cols(dim_i)%dat(idx_n(dim_i)) <= idxend(dim_i))flag = flag-1
                     enddo
                     if(flag==0)subtensors(index_ij)%masks(i, j) = 0
                  enddo
               endif
            enddo
         enddo
      endif


      do dim_i=1,Ndim*2
         deallocate(select_idx(dim_i)%dat)
      enddo
      deallocate(select_idx)

   end subroutine BF_MD_compress_N_oneblock_C_sample



   subroutine BF_MD_compress_N_oneblock_C_rankreveal(Ndim, dim_MD, subtensors, blocks, option, stats, msh, ker, ptree, index_ij, bb_m, level, rank_new, flops)


      implicit none

      integer bb_m, dim,dim_i,dim_ii,flag,index_jj_scalar,Ndim, cnt
      integer dim_MD(Ndim+2),idx_MD(Ndim+2), idx_r_m(Ndim), dims_m(Ndim), dims_row(Ndim), dims_row1(Ndim), dims_col(Ndim),idx_m(Ndim),idx_m1(Ndim),idx_n(Ndim), dims_bm(Ndim), group_scalar
      type(intersect_MD) :: subtensors(:)
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      integer i, j, i1, i_dim, level_butterfly, header_m(Ndim), header_n(Ndim), rankmax_r, rankmax_c
      integer group_m(Ndim), group_n(Ndim), group_m_mid(Ndim), group_n_mid(Ndim), idxstart(Ndim), idxend(Ndim), mm(Ndim), nn(Ndim), mm_scalar, mmnn(Ndim*2), index_j(Ndim), index_j_scalar, index_ij, index_i, jj, index_i_s, index_i_k, mmm1
      integer level
      integer header_m1, header_m2, mm1, mm2, index_jj(Ndim), index_ii
      real(kind=8) flops,flop
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row_tmp(:), select_column(:), column_pivot(:), row_pivot(:)
      type(iarray),allocatable::select_idx(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:), order(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_tmp(:, :), matrix_little_cc(:, :), core(:, :), core_tmp(:, :), core_tmp1(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, passflag, levelm, nrow, ncol,rank_new

      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      real(kind=8)::n2, n1,overrate

      flops = 0
      level_butterfly = blocks%level_butterfly
      levelm = floor_safe(dble(level_butterfly)/2d0)

      call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
      index_i = idx_MD(1)
      dim = idx_MD(Ndim+2)
      index_j = idx_MD(2:Ndim+1)
      ! dims_m = 2**(levelm)
      call SingleIndexToMultiIndex(Ndim, blocks%nr_m, bb_m, idx_r_m)


      group_m = blocks%row_group
      group_n = blocks%col_group
      group_n = group_n*2**(level_butterfly - level + 1) - 1 + index_j

      group_m = group_m*2**(levelm) + blocks%idx_r_m - 1 + idx_r_m -1
      group_m(dim) = group_m(dim)*2**(level - levelm -1) + index_i - 1

      do dim_i=1,Ndim
         mm(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
      enddo
      do dim_i=1,Ndim
         nn(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
      enddo

      if (level == 0) then
         write(*,*)"should not come to level == 0 in BF_MD_compress_N_oneblock_C_sample"
      elseif (level == level_butterfly + 1) then
         mm_scalar = mm(dim)
      else
         do dim_i=1,Ndim
         index_jj(dim_i) = int((index_j(dim_i) + 1)/2)
         enddo
         index_ii = 2*index_i - 1
         dims_col = 2**(level_butterfly - level+1)
         call MultiIndexToSingleIndex(Ndim, dims_col, index_j_scalar, index_j)
         dims_col = 2**(level_butterfly - level)
         call MultiIndexToSingleIndex(Ndim, dims_col, index_jj_scalar, index_jj)

         mm1 = size(blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii, index_jj_scalar, dim)%array, 1)
         mm2 = size(blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii+1, index_jj_scalar, dim)%array, 1)
         mm_scalar = mm1 + mm2
      endif
      mmnn(1:Ndim)=mm
      mmnn(1+Ndim:Ndim*2)=nn

      dims_row = subtensors(index_ij)%nr
      dims_col = subtensors(index_ij)%nc

      do dim_i =1,Ndim
         header_m(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%head
         header_n(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%head
         if (level < level_butterfly+1 .and. dim==dim_i) then
            header_m1 = header_m(dim_i)
            header_m2 = msh(dim_i)%basis_group(2*group_m(dim_i)+1)%head
         endif
      enddo


      rankmax_r = dims_row(dim)
      rankmax_c=product(dims_col)
      do dim_i=1,Ndim
         if(dim_i/=dim)rankmax_c = rankmax_c*dims_row(dim_i)
      enddo
      allocate (core(rankmax_c, rankmax_r))
core=0
      allocate (core_tmp(rankmax_c, rankmax_r))
      allocate (matrix_tmp(product(dims_col), product(dims_row)))
      call copymatT(subtensors(index_ij)%dat, matrix_tmp, product(dims_row), product(dims_col))

      ! reshaping
do i_dim = 1,dims_row(dim)
         cnt=0
      do i = 1, product(dims_row)
         call SingleIndexToMultiIndex(Ndim, dims_row, i, idx_m)
         if(idx_m(dim)==i_dim)then ! this makes sure that all the other 2*Ndim-1 indices are traversed before moving to the next idx_m(dim) value
            do j=1,product(dims_col)
               if(subtensors(index_ij)%masks(i,j)==1)then
                  cnt = cnt+1
                  core(cnt,idx_m(dim)) = matrix_tmp(j,i)
         endif
         enddo
         endif
         enddo
      enddo
      core_tmp=core


      allocate (jpvt(max(rankmax_c, rankmax_r)))
      allocate (tau(max(rankmax_c, rankmax_r)))
      jpvt = 0
      call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
      flops = flops + flop

      if (rank_new > 0) then
         call un_or_mqrf90(core, tau, core_tmp, 'L', 'C', rankmax_c, rankmax_r, rank_new, flop=flop)
         flops = flops + flop
         call trsmf90(core, core_tmp, 'L', 'U', 'N', 'N', rank_new, rankmax_r, flop=flop)
         flops = flops + flop
      else
         rank_new = 1
         jpvt(1) = 1
         core_tmp = 0
      endif

      if (level == 0) then
        write(*,*)"should not arrive at level == 0 in BF_MD_compress_N_oneblock_C_rankreveal"
      elseif (level == level_butterfly + 1) then
         index_i_k=index_i
         index_i_s=index_i
         allocate (blocks%ButterflyU(bb_m)%blocks(index_i_k,dim)%matrix(rankmax_r, rank_new))
         call copymatT(core_tmp(1:rank_new, 1:rankmax_r), blocks%ButterflyU(bb_m)%blocks(index_i_k,dim)%matrix, rank_new, rankmax_r)
         allocate (blocks%ButterflySkel_L(bb_m,level)%inds(index_i_s,1,dim)%array(rank_new))
         do j = 1, rank_new
            blocks%ButterflySkel_L(bb_m,level)%inds(index_i_s,1,dim)%array(j) = jpvt(j)
         enddo
      else
         index_i_k=2*index_i-1
         index_i_s=index_i
         mmm1 = msh(dim)%basis_group(2*group_m(dim))%tail - msh(dim)%basis_group(2*group_m(dim))%head + 1
         allocate (blocks%ButterflyKerl_L(bb_m,level)%blocks(index_i_k,index_j_scalar,dim)%matrix(mm1,rank_new))
         allocate (blocks%ButterflyKerl_L(bb_m,level)%blocks(index_i_k+1,index_j_scalar,dim)%matrix(mm2,rank_new))
         call copymatT(core_tmp(1:rank_new, 1:mm1), blocks%ButterflyKerl_L(bb_m,level)%blocks(index_i_k,index_j_scalar,dim)%matrix, rank_new, mm1)
         call copymatT(core_tmp(1:rank_new, 1 + mm1:rankmax_r), blocks%ButterflyKerl_L(bb_m,level)%blocks(index_i_k+1,index_j_scalar,dim)%matrix, rank_new, rankmax_r - mm1)

         allocate (blocks%ButterflySkel_L(bb_m,level)%inds(index_i_s,index_j_scalar,dim)%array(rank_new))
         do j = 1, rank_new
            if (jpvt(j) <= mm1) then
               blocks%ButterflySkel_L(bb_m,level)%inds(index_i_s,index_j_scalar,dim)%array(j) = blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii, index_jj_scalar,dim)%array(jpvt(j))
            else
               blocks%ButterflySkel_L(bb_m,level)%inds(index_i_s,index_j_scalar,dim)%array(j) = blocks%ButterflySkel_L(bb_m,level + 1)%inds(index_ii+1, index_jj_scalar,dim)%array(jpvt(j) - mm1) + mmm1
            endif
         enddo

      endif

      deallocate(core)
      deallocate(core_tmp)
      deallocate(jpvt)
      deallocate(tau)
      deallocate(matrix_tmp)


   end subroutine BF_MD_compress_N_oneblock_C_rankreveal

   subroutine BF_exchange_skel(blocks, skels, option, stats, msh, ptree, level, mode, collect)


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

      n1 = MPI_Wtime()

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
         sendquant(ii)%size_i = 0
         sendquant(ii)%active = 0
      enddo
      allocate (recvquant(nproc))
      do ii = 1, nproc
         recvquant(ii)%size_i = 0
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
            sendquant(pp)%size_i = sendquant(pp)%size_i + 3 + size(skels%inds(ii, jj)%array)
         endif
      enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat_i(sendquant(pp)%size_i, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         else
            recvquant(pp)%size_i = sendquant(pp)%size_i
         endif
      enddo
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size_i = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat_i(recvquant(pp)%size_i, 1))
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
            sendquant(pp)%dat_i(sendquant(pp)%size_i + 1, 1) = index_i
            sendquant(pp)%dat_i(sendquant(pp)%size_i + 2, 1) = index_j
            sendquant(pp)%dat_i(sendquant(pp)%size_i + 3, 1) = Nskel
            sendquant(pp)%size_i = sendquant(pp)%size_i + 3
            do i = 1, Nskel
               sendquant(pp)%dat_i(sendquant(pp)%size_i + i, 1) = skels%inds(ii, jj)%array(i)
            enddo
            sendquant(pp)%size_i = sendquant(pp)%size_i + Nskel
         endif
      enddo
      enddo

      ! communicate the data buffer
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr =Nreqr+1
            call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         else
            if(recvquant(pp)%size_i>0)recvquant(pp)%dat_i = sendquant(pp)%dat_i
         endif
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         i = 0
         do while (i < recvquant(pp)%size_i)
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
      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_exchange_skel

!>*********** all to all communication of skeletons of one butterfly level from row-wise ordering to column-wise ordering or the reverse
   subroutine BF_all2all_skel(blocks, skels, option, stats, msh, ptree, level, mode, mode_new)


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

      n1 = MPI_Wtime()

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
         sendquant(ii)%size_i = 0
         sendquant(ii)%active = 0
      enddo
      allocate (recvquant(nproc))
      do ii = 1, nproc
         recvquant(ii)%size_i = 0
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
         sendquant(pp)%size_i = sendquant(pp)%size_i + 3 + size(skels%inds(ii, jj)%array)
      enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat_i(sendquant(pp)%size_i, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid/=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size_i = sendquant(pp)%size_i
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid/=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size_i = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat_i(recvquant(pp)%size_i, 1))
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
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 1, 1) = index_i
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 2, 1) = index_j
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 3, 1) = Nskel
         sendquant(pp)%size_i = sendquant(pp)%size_i + 3
         do i = 1, Nskel
            sendquant(pp)%dat_i(sendquant(pp)%size_i + i, 1) = skels%inds(ii, jj)%array(i)
         enddo
         deallocate (skels%inds(ii, jj)%array)
         sendquant(pp)%size_i = sendquant(pp)%size_i + Nskel
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
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size_i > 0) recvquant(pp)%dat_i = sendquant(pp)%dat_i
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(blocks%pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size_i)
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
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

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

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_skel



!>*********** all to all communication of skeletons at the middle butterfly level from row-wise ordering to column-wise ordering
   subroutine BF_MD_all2all_skel(Ndim, blocks, ButterflySkel_R_transposed, option, stats, msh, ptree)


      implicit none
      integer Ndim
      type(mesh)::msh(Ndim)
      integer i, j, levelm, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, pgno_sub, tag, nproc, Ncol, Nskel(Ndim), Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new, dim_i, bb

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      character::mode, mode_new
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      integer, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer,allocatable::proc_of_groupc(:),proc_of_groupr(:)
      integer:: dims_r(Ndim),dims_c(Ndim),idx_r_m(Ndim),idx_c_m(Ndim),idx_r_scalar,idx_c_scalar
      type(butterfly_skel_MD):: ButterflySkel_R_transposed(:)

      n1 = MPI_Wtime()

      levelm = blocks%level_half
      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno
      dims_r = 2**levelm
      dims_c = 2**(level_butterfly-levelm)

      ! computation of the process ID of each middle level group
      allocate(proc_of_groupr(product(dims_r)))
      proc_of_groupr=0
      allocate(proc_of_groupc(product(dims_c)))
      proc_of_groupc=0
      pp = ptree%MyID - ptree%pgrp(blocks%pgno)%head + 1
      do bb=1, product(blocks%nr_m)
         call SingleIndexToMultiIndex(Ndim, blocks%nr_m, bb, idx_r_m)
         idx_r_m = idx_r_m + blocks%idx_r_m - 1
         call MultiIndexToSingleIndex(Ndim, dims_r, idx_r_scalar, idx_r_m)
         proc_of_groupr(idx_r_scalar)=pp
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, proc_of_groupr, product(dims_r), MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)

      do bb=1, product(blocks%nc_m)
         call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb, idx_c_m)
         idx_c_m = idx_c_m + blocks%idx_c_m - 1
         call MultiIndexToSingleIndex(Ndim, dims_c, idx_c_scalar, idx_c_m)
         proc_of_groupc(idx_c_scalar)=pp
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, proc_of_groupc, product(dims_c), MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)



      do bb=1, product(blocks%nr_m)
         allocate(ButterflySkel_R_transposed(bb)%nc(Ndim))
         ButterflySkel_R_transposed(bb)%nc = blocks%ButterflySkel_L(bb,levelm+1)%nc
         allocate(ButterflySkel_R_transposed(bb)%nr(Ndim))
         ButterflySkel_R_transposed(bb)%nr = blocks%ButterflySkel_L(bb,levelm+1)%nr
         allocate(ButterflySkel_R_transposed(bb)%idx_r(Ndim))
         ButterflySkel_R_transposed(bb)%idx_r=blocks%ButterflySkel_L(bb,levelm+1)%idx_r
         allocate(ButterflySkel_R_transposed(bb)%inc_r(Ndim))
         ButterflySkel_R_transposed(bb)%inc_r=blocks%ButterflySkel_L(bb,levelm+1)%inc_r
         allocate(ButterflySkel_R_transposed(bb)%idx_c(Ndim))
         ButterflySkel_R_transposed(bb)%idx_c=blocks%ButterflySkel_L(bb,levelm+1)%idx_c
         allocate(ButterflySkel_R_transposed(bb)%inc_c(Ndim))
         ButterflySkel_R_transposed(bb)%inc_c=blocks%ButterflySkel_L(bb,levelm+1)%inc_c
         allocate (ButterflySkel_R_transposed(bb)%inds(ButterflySkel_R_transposed(bb)%nr(1), product(ButterflySkel_R_transposed(bb)%nc),Ndim))
      enddo



      ! allocation of communication quantities
      allocate (statuss(MPI_status_size, nproc))
      allocate (statusr(MPI_status_size, nproc))
      allocate (S_req(nproc))
      allocate (R_req(nproc))
      allocate (sendquant(nproc))
      do ii = 1, nproc
         sendquant(ii)%size_i = 0
         sendquant(ii)%active = 0
      enddo
      allocate (recvquant(nproc))
      do ii = 1, nproc
         recvquant(ii)%size_i = 0
         recvquant(ii)%active = 0
      enddo
      allocate (sendIDactive(nproc))
      allocate (recvIDactive(nproc))
      Nsendactive = 0
      Nrecvactive = 0

      ! calculate send buffer sizes in the first pass
      do bb=1, product(blocks%nr_m)
         call assert(ButterflySkel_R_transposed(bb)%nr(1)==1,'ButterflySkel_R_transposed(bb)%nr(1) should be 1')
         do ii = 1, ButterflySkel_R_transposed(bb)%nr(1)
            do jj=1,product(ButterflySkel_R_transposed(bb)%nc)
               pp = proc_of_groupc(jj)
               if (recvquant(pp)%active == 0) then
                  recvquant(pp)%active = 1
                  Nrecvactive = Nrecvactive + 1
                  recvIDactive(Nrecvactive) = pp
               endif
            enddo
         enddo
      enddo

      do bb=1, product(blocks%nc_m)
         do ii = 1, product(blocks%ButterflySkel_R(bb,levelm)%nr)
            call assert(blocks%ButterflySkel_R(bb,levelm)%nc(1)==1,'blocks%ButterflySkel_R(bb,levelm)%nc(1) should be 1')
            do jj=1, blocks%ButterflySkel_R(bb,levelm)%nc(1)
               pp = proc_of_groupr(ii)
               if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
               endif
               sendquant(pp)%size_i = sendquant(pp)%size_i + 2 + Ndim
               do dim_i=1,Ndim
                  sendquant(pp)%size_i = sendquant(pp)%size_i + size(blocks%ButterflySkel_R(bb,levelm)%inds(ii, jj, dim_i)%array)
               enddo
            enddo
         enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat_i(sendquant(pp)%size_i, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid/=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size_i = sendquant(pp)%size_i
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid/=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size_i, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size_i = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat_i(recvquant(pp)%size_i, 1))
      enddo

      ! pack the send buffer in the second pass
      do bb=1, product(blocks%nc_m)
         do ii = 1, product(blocks%ButterflySkel_R(bb,levelm)%nr)
            do jj=1, blocks%ButterflySkel_R(bb,levelm)%nc(1)
               pp = proc_of_groupr(ii)
               call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb, idx_c_m)
               idx_c_m = idx_c_m + blocks%idx_c_m - 1
               call MultiIndexToSingleIndex(Ndim, dims_c, idx_c_scalar, idx_c_m)
               sendquant(pp)%dat_i(sendquant(pp)%size_i + 1, 1) = ii ! global row index at the middle level
               sendquant(pp)%dat_i(sendquant(pp)%size_i + 2, 1) = idx_c_scalar ! global column index at the middle level
               do dim_i=1,Ndim
                  sendquant(pp)%dat_i(sendquant(pp)%size_i + 2 + dim_i, 1) = size(blocks%ButterflySkel_R(bb,levelm)%inds(ii, jj, dim_i)%array)
               enddo
               sendquant(pp)%size_i = sendquant(pp)%size_i + 2 + Ndim

               do dim_i=1,Ndim
                  do i = 1, size(blocks%ButterflySkel_R(bb,levelm)%inds(ii, jj, dim_i)%array)
                     sendquant(pp)%dat_i(sendquant(pp)%size_i + i, 1) = blocks%ButterflySkel_R(bb,levelm)%inds(ii, jj, dim_i)%array(i)
                  enddo
                  sendquant(pp)%size_i = sendquant(pp)%size_i + size(blocks%ButterflySkel_R(bb,levelm)%inds(ii, jj, dim_i)%array)
               enddo
            enddo
         enddo
      enddo

      ! communicate the data buffer
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size_i > 0) recvquant(pp)%dat_i = sendquant(pp)%dat_i
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size_i, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(blocks%pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size_i)
            i = i + 1
            idx_r_scalar = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            call SingleIndexToMultiIndex(Ndim, dims_r, idx_r_scalar, idx_r_m)
            idx_r_m = idx_r_m - blocks%idx_r_m + 1
            call MultiIndexToSingleIndex(Ndim, blocks%nr_m, bb, idx_r_m)

            i = i + 1
            jj = NINT(dble(recvquant(pp)%dat_i(i, 1)))

            do dim_i=1,Ndim
               i = i + 1
               Nskel(dim_i) = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            enddo
            do dim_i=1,Ndim
               call assert(.not. allocated(ButterflySkel_R_transposed(bb)%inds(1,jj,dim_i)%array), 'receiving dat_ia alreay exists locally')
               allocate (ButterflySkel_R_transposed(bb)%inds(1,jj,dim_i)%array(Nskel(dim_i)))
               ButterflySkel_R_transposed(bb)%inds(1,jj,dim_i)%array = NINT(dble(recvquant(pp)%dat_i(i + 1:i + Nskel(dim_i), 1)))
               i = i + Nskel(dim_i)
            enddo
         enddo
      enddo
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

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

      deallocate (proc_of_groupr)
      deallocate (proc_of_groupc)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_MD_all2all_skel




   subroutine BF_compress_NlogN_oneblock_C_sample(submats, blocks, boundary_map, Nboundall, Ninadmissible,groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, level_final, nnz_loc)


      implicit none

      type(intersect) :: submats(:)
      integer Nboundall,Ninadmissible
      integer*8 nnz_loc
      integer boundary_map(:,:)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, iii, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_c1, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_j, index_ij, index_i_loc_k, index_j_loc_k, index_i_loc_s, index_i_loc_s1, index_j_loc_s, index_j_loc_s1, ii, jj, ij
      integer level, level_final, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2, inter
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
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      integer::passflag = 0
      real(kind=8)::n2, n1,overrate


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
      overrate=1
      if(level>levelm)then
      group_m_mid = group_m
      do i=1,level-levelm
         group_m_mid = floor_safe(group_m_mid/2d0)
      enddo
      idxstart = msh%basis_group(group_m_mid)%head
      idxend = msh%basis_group(group_m_mid)%tail
      inter = min(msh%basis_group(group_m_mid)%tail,msh%basis_group(group_n)%tail)-max(msh%basis_group(group_m_mid)%head,msh%basis_group(group_n)%head)+1
      if(inter>0)then
         if(inter==nn)then
         overrate = 1
         else
         overrate = dble(nn)/dble(nn-inter)
         endif
      endif
      endif



      ! select skeletons here, selection of at most (option%sample_para+option%knn)*mm rows, the first option%sample_para*mm are random, the next option%knn*mm are nearest points
      rankmax_r = mm
      rankmax_c1 = min(nn, ceiling_safe(option%sample_para*mm*overrate))
      if (level == level_butterfly + 1) rankmax_c1 = min(ceiling_safe(option%sample_para_outer*mm*overrate), nn)
      allocate (select_column(rankmax_c1 + option%knn*mm))
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


      call remove_dup_int(select_column, rankmax_c1, rankmax_c)

      if (level == level_butterfly + 1) then
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head

         if (level == level_final) then
            index_i_loc_s1 = (index_i - blocks%ButterflySkel(level - 1)%idx_r)/blocks%ButterflySkel(level - 1)%inc_r + 1
            index_j_loc_s1 = (index_j - blocks%ButterflySkel(level - 1)%idx_c)/blocks%ButterflySkel(level - 1)%inc_c + 1

            rank_new = size(blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array)
            index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
            ! allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
            if (Nboundall > 0) then
               allocate (submats(index_ij)%dat(mm, rank_new))
               submats(index_ij)%dat=0
               call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
            endif
            nnz_loc = nnz_loc + mm*rank_new

            allocate (submats(index_ij)%rows(mm))
            allocate (submats(index_ij)%cols(rank_new))
            submats(index_ij)%nr = mm
            submats(index_ij)%nc = rank_new
            do i = 1, mm
               submats(index_ij)%rows(i) = header_m + i - 1
            enddo
            do j = 1, rank_new
               submats(index_ij)%cols(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (submats(index_ij)%masks(mm, rank_new))
               call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
               submats(index_ij)%masks = 1

               do i = 1, mm
                  group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
                  do jj=1,Ninadmissible
                     group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                     if (group_n_mid /= -1) then
                        idxstart = msh%basis_group(group_n_mid)%head
                        idxend = msh%basis_group(group_n_mid)%tail
                        do j = 1, rank_new
                           if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                        enddo
                     endif
                  enddo
               enddo
            endif

         else
            if (Nboundall > 0) then
               allocate (submats(index_ij)%dat(mm,rankmax_c))
               submats(index_ij)%dat=0
               call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
            endif
            nnz_loc = nnz_loc + mm*rankmax_c

            allocate (submats(index_ij)%rows(mm))
            allocate (submats(index_ij)%cols(rankmax_c))
            submats(index_ij)%nr = mm
            submats(index_ij)%nc = rankmax_c
            do i = 1, mm
               submats(index_ij)%rows(i) = header_m + i - 1
            enddo
            do j = 1, rankmax_c
               submats(index_ij)%cols(j) = header_n + select_column(j) - 1
            enddo

            if (Nboundall > 0) then

               allocate (submats(index_ij)%masks(mm, rankmax_c))
               call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
               submats(index_ij)%masks = 1

               do i = 1, mm
                  group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
                  do jj=1,Ninadmissible
                     group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                     if (group_n_mid /= -1) then
                        idxstart = msh%basis_group(group_n_mid)%head
                        idxend = msh%basis_group(group_n_mid)%tail
                        do j = 1, rankmax_c
                           if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                        enddo
                     endif
                  enddo
               enddo
            endif
         endif

      elseif (level == 0) then
         index_j_loc_s = (index_j - blocks%ButterflySkel(level + 1)%idx_c)/blocks%ButterflySkel(level + 1)%inc_c + 1
         rank_new = size(blocks%ButterflySkel(level + 1)%inds(1, index_j_loc_s)%array)
         header_m = msh%basis_group(group_m)%head
         header_n = msh%basis_group(group_n)%head
         index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1

         ! allocate (blocks%ButterflyV%blocks(index_j_loc_k)%matrix(nn, rank_new))
         ! allocate (matrix_V_tmp(rank_new, nn))
         if (Nboundall > 0) then
            allocate (submats(index_ij)%dat(rank_new, nn))
            call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
            submats(index_ij)%dat=0
         endif
         nnz_loc = nnz_loc + rank_new*nn

         allocate (submats(index_ij)%rows(rank_new))
         allocate (submats(index_ij)%cols(nn))
         submats(index_ij)%nr = rank_new
         submats(index_ij)%nc = nn

         do i = 1, rank_new
            submats(index_ij)%rows(i) = blocks%ButterflySkel(level + 1)%inds(1, index_j_loc_s)%array(i) + header_m - 1
         enddo
         do j = 1, nn
            submats(index_ij)%cols(j) = j + header_n - 1
         enddo

         if (Nboundall > 0) then

            allocate (submats(index_ij)%masks(rank_new, nn))
            call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)
            submats(index_ij)%masks = 1

            do i = 1, rank_new
               group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
               do jj=1,Ninadmissible
                  group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                  if (group_n_mid /= -1) then
                     idxstart = msh%basis_group(group_n_mid)%head
                     idxend = msh%basis_group(group_n_mid)%tail
                     do j = 1, nn
                        if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                     enddo
                  endif
               enddo
            enddo
         endif

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
            if (Nboundall > 0) then
               allocate (submats(index_ij)%dat(mm1 + mm2, rank_new))
               submats(index_ij)%dat=0
               call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
            endif
            nnz_loc = nnz_loc + (mm1 + mm2)*rank_new

            allocate (submats(index_ij)%rows(mm1 + mm2))
            allocate (submats(index_ij)%cols(rank_new))
            submats(index_ij)%nr = mm1+mm2
            submats(index_ij)%nc = rank_new
            do i = 1, mm1 + mm2
               if (i <= mm1) then
                  submats(index_ij)%rows(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(i) + header_m1 - 1
               else
                  submats(index_ij)%rows(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(i - mm1) + header_m2 - 1
               endif
            enddo
            do j = 1, rank_new
               submats(index_ij)%cols(j) = blocks%ButterflySkel(level - 1)%inds(index_i_loc_s1, index_j_loc_s1)%array(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (submats(index_ij)%masks(mm1 + mm2, rank_new))
               submats(index_ij)%masks = 1
               call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)

               do i = 1, mm1 + mm2
                  group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
                  do jj=1,Ninadmissible
                     group_n_mid = boundary_map(group_m_mid - groupm_start + 1, jj)
                     if (group_n_mid /= -1) then
                        idxstart = msh%basis_group(group_n_mid)%head
                        idxend = msh%basis_group(group_n_mid)%tail
                        do j = 1, rank_new
                           if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                        enddo
                     endif
                  enddo
               enddo
            endif

            ! deallocate (matrix_V_tmp)
         else
            index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
            index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1
            if (Nboundall > 0) then
               allocate (submats(index_ij)%dat(mm,rankmax_c))
               submats(index_ij)%dat=0
               call LogMemory(stats, SIZEOF(submats(index_ij)%dat)/1024.0d3)
            endif
            nnz_loc = nnz_loc + mm*rankmax_c
            allocate (submats(index_ij)%rows(mm))
            allocate (submats(index_ij)%cols(rankmax_c))
            submats(index_ij)%nr = mm
            submats(index_ij)%nc = rankmax_c
            do i = 1, mm
               if (i <= mm1) then
                  submats(index_ij)%rows(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc, index_jj_loc)%array(i) + header_m1 - 1
               else
                  submats(index_ij)%rows(i) = blocks%ButterflySkel(level + 1)%inds(index_ii_loc + 1, index_jj_loc)%array(i - mm1) + header_m2 - 1
               endif
            enddo
            do j = 1, rankmax_c
               submats(index_ij)%cols(j) = select_column(j) + header_n - 1
            enddo

            if (Nboundall > 0) then

               allocate (submats(index_ij)%masks(mm, rankmax_c))
               submats(index_ij)%masks = 1
               call LogMemory(stats, SIZEOF(submats(index_ij)%masks)/1024.0d3)

               do i = 1, mm
                  group_m_mid = findgroup(submats(index_ij)%rows(i), msh, levelm, blocks%row_group)
                  do jj=1,Ninadmissible
                     group_n_mid = boundary_map(group_m_mid - groupm_start + 1,jj)
                     if (group_n_mid /= -1) then
                        idxstart = msh%basis_group(group_n_mid)%head
                        idxend = msh%basis_group(group_n_mid)%tail
                        do j = 1, rankmax_c
                           if (submats(index_ij)%cols(j) >= idxstart .and. submats(index_ij)%cols(j) <= idxend) submats(index_ij)%masks(i, j) = 0
                        enddo
                     endif
                  enddo
               enddo
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
   end subroutine BF_compress_NlogN_oneblock_C_sample

   subroutine BF_compress_NlogN_oneblock_C_rankreveal(submats, blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_i, index_j, index_ij, level, level_final, rank_new,flops)


      implicit none

      type(intersect) :: submats(:)
      integer Nboundall,Ninadmissible
      integer boundary_map(:,:)
      integer groupm_start

      type(mesh)::msh
      type(kernelquant)::ker
      integer i, j, iii, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_c1, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_j, index_ij, index_i_loc_k, index_j_loc_k, index_i_loc_s, index_i_loc_s1, index_j_loc_s, index_j_loc_s1, ii, jj, ij
      integer level, level_final, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2, inter
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
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), mrange(:), nrange(:), mrange1(:), nrange1(:), mmap(:), nmap(:)
      integer::passflag = 0
      real(kind=8)::n2, n1,overrate


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

      rankmax_r = mm
      rankmax_c = size(submats(index_ij)%dat,2)
      allocate (select_row(rankmax_r))
      do i = 1, rankmax_r
         select_row(i) = i
      enddo

      if (level == level_butterfly + 1) then
         if (level == level_final) then
            index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
            mm = size(submats(index_ij)%dat,1)
            rank_new = size(submats(index_ij)%dat,2)
            allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
            blocks%ButterflyU%blocks(index_i_loc_k)%matrix = submats(index_ij)%dat

         else
            allocate (core(rankmax_c, rankmax_r))
            allocate (core_tmp(rankmax_c, rankmax_r))
            call copymatT(submats(index_ij)%dat, core, rankmax_r, rankmax_c)
            call copymatT(submats(index_ij)%dat, core_tmp, rankmax_r, rankmax_c)
            allocate (jpvt(max(rankmax_c, rankmax_r)))
            allocate (tau(max(rankmax_c, rankmax_r)))
            jpvt = 0
            call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
            flops = flops + flop

            if (rank_new > 0) then
               call un_or_mqrf90(core, tau, core_tmp, 'L', 'C', rankmax_c, mm, rank_new, flop=flop)
               flops = flops + flop
               call trsmf90(core, core_tmp, 'L', 'U', 'N', 'N', rank_new, mm, flop=flop)
               flops = flops + flop
            else
               rank_new = 1
               jpvt(1) = 1
               core_tmp = 0
            endif

            index_i_loc_k = (index_i - blocks%ButterflyU%idx)/blocks%ButterflyU%inc + 1
            allocate (blocks%ButterflyU%blocks(index_i_loc_k)%matrix(mm, rank_new))
            call copymatT(core_tmp(1:rank_new, 1:mm), blocks%ButterflyU%blocks(index_i_loc_k)%matrix, rank_new, mm)

            index_i_loc_s = (index_i - blocks%ButterflySkel(level_butterfly + 1)%idx_r)/blocks%ButterflySkel(level_butterfly + 1)%inc_r + 1
            allocate (blocks%ButterflySkel(level_butterfly + 1)%inds(index_i_loc_s, 1)%array(rank_new))
            do j = 1, rank_new
               blocks%ButterflySkel(level_butterfly + 1)%inds(index_i_loc_s, 1)%array(j) = select_row(jpvt(j))
            enddo
         endif
      elseif (level == 0) then
         index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1
         rank_new = size(submats(index_ij)%dat,1)
         nn = size(submats(index_ij)%dat,2)
         allocate (blocks%ButterflyV%blocks(index_j_loc_k)%matrix(nn, rank_new))
         call copymatT(submats(index_ij)%dat, blocks%ButterflyV%blocks(index_j_loc_k)%matrix, rank_new, nn)
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
            if (mm1 > 0) blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = submats(index_ij)%dat(1:mm1, :)
            if (mm2 > 0) blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix = submats(index_ij)%dat(1 + mm1:mm1 + mm2, :)
         else
            index_i_loc_s = (index_i - blocks%ButterflySkel(level)%idx_r)/blocks%ButterflySkel(level)%inc_r + 1
            index_j_loc_s = (index_j - blocks%ButterflySkel(level)%idx_c)/blocks%ButterflySkel(level)%inc_c + 1
            allocate (core(rankmax_c, rankmax_r))
            allocate (core_tmp(rankmax_c, rankmax_r))
            call copymatT(submats(index_ij)%dat, core, rankmax_r, rankmax_c)
            call copymatT(submats(index_ij)%dat, core_tmp, rankmax_r, rankmax_c)
            allocate (jpvt(max(rankmax_c, rankmax_r)))
            allocate (tau(max(rankmax_c, rankmax_r)))
            jpvt = 0
            call geqp3modf90(core, jpvt, tau, option%tol_comp, BPACK_SafeUnderflow, rank_new, flop=flop)
            flops = flops + flop

            if (rank_new > 0) then
               call un_or_mqrf90(core, tau, core_tmp, 'L', 'C', rankmax_c, mm, rank_new, flop=flop)
               flops = flops + flop
               call trsmf90(core, core_tmp, 'L', 'U', 'N', 'N', rank_new, mm, flop=flop)
               flops = flops + flop
            else
               rank_new = 1
               jpvt(1) = 1
               core_tmp = 0
            endif

            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(mm1, rank_new))
            allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix(mm2, rank_new))

            call copymatT(core_tmp(1:rank_new, 1:mm1), blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, rank_new, mm1)
            call copymatT(core_tmp(1:rank_new, 1 + mm1:mm), blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, rank_new, mm - mm1)

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

         endif

      endif

      if (allocated(core)) deallocate (core)
      if (allocated(core_tmp)) deallocate (core_tmp)
      if (allocated(tau)) deallocate (tau)
      if (allocated(jpvt)) deallocate (jpvt)
      if (allocated(matrix_V)) deallocate (matrix_V)
      if (allocated(select_row)) deallocate (select_row)
      if (allocated(select_column)) deallocate (select_column)
   end subroutine BF_compress_NlogN_oneblock_C_rankreveal


  subroutine BF_MD_delete_subtensors(Ndim, dims, subtensors, stats)
  implicit none
  integer Ndim,index_ij,dim_i
  integer dims(:)
  type(Hstat)::stats
  type(intersect_MD) :: subtensors(:)

   do index_ij = 1, product(dims)
      call LogMemory(stats, -SIZEOF(subtensors(index_ij)%dat)/1024.0d3)

      if(associated(subtensors(index_ij)%dat))deallocate(subtensors(index_ij)%dat)
      if(allocated(subtensors(index_ij)%rows))then
         do dim_i=1,Ndim
            call iarray_finalizer(subtensors(index_ij)%rows(dim_i))
         enddo
         deallocate(subtensors(index_ij)%rows)
      endif
      if(allocated(subtensors(index_ij)%nr))deallocate(subtensors(index_ij)%nr)
      if(allocated(subtensors(index_ij)%cols))then
         do dim_i=1,Ndim
            call iarray_finalizer(subtensors(index_ij)%cols(dim_i))
         enddo
         deallocate(subtensors(index_ij)%cols)
      endif
      if(allocated(subtensors(index_ij)%nc))deallocate(subtensors(index_ij)%nc)
      if(allocated(subtensors(index_ij)%masks))then
         call LogMemory(stats, -SIZEOF(subtensors(index_ij)%masks)/1024.0d3)
         deallocate(subtensors(index_ij)%masks)
      endif
   enddo

  end subroutine BF_MD_delete_subtensors



   subroutine BF_MD_compress_N(Ndim,blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)


      implicit none

      integer Nboundall, Ninadmissible, statflag, Ndim, dim_i
      integer boundary_map(:,:,:)
      integer groupm_start(Ndim)

      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      integer i, j, bb, bb_g, level_butterfly, level_butterflyL, level_butterflyR, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new, rank_new1
      integer group_m(Ndim), group_n(Ndim), mm, nn, index_i, index_i_loc, index_j_loc, index_j, ii, jj, ij
      integer level, length_1, length_2, level_blocks, index_ij
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e, flops, flops1
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, nnn1, ierr
      real(kind=8) Memory, flop,n2,n1, Memory_dense, Memory_comp
      DT ctemp
      DT,target,allocatable:: alldat_loc_in(:)
      type(matrixblock_MD)::blocks
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer Nlayer, level_half, level_final, dim_MD(Ndim+2), idx_MD(Ndim+2), idx_r(Ndim), inc_r(Ndim), nr(Ndim), idx_c(Ndim), inc_c(Ndim), nc(Ndim), idx_c_scalar, idx_r_scalar, dim_subtensor(Ndim*2),idx_subtensor(Ndim*2)
      integer passflag,use_zfp
      integer*8 nnz_loc
      integer, allocatable :: rankmax_for_butterfly(:), rankmin_for_butterfly(:), select_row_pre(:), select_col_pre(:)
      integer::mrange_dummy(1), nrange_dummy(1)
      type(intersect_MD), allocatable :: subtensors(:)
      type(intersect_MD) :: subtensors_dummy(1)
      type(butterfly_skel_MD), allocatable :: ButterflySkel_R_transposed(:)

      ! DT:: mat_dummy(1, 1)
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

      call assert(option%pat_comp==3, 'option%pat_comp can only be set to 3 for BF_MD_compress_N')
      num_blocks = 2**level_butterfly
      blocks%level_half = BF_Switchlevel(level_butterfly, option%pat_comp)
      level_half = blocks%level_half




      if (level_butterfly == 0 .and. ptree%pgrp(blocks%pgno)%nproc>1) then
            write(*,*)'level_butterfly = 0 has not been implemented in parallel for BF_MD_compress_N'
      else
         level_butterflyL = level_butterfly-blocks%level_half
         level_butterflyR = blocks%level_half
         allocate(blocks%nr_m(Ndim))
         allocate(blocks%nc_m(Ndim))
         allocate(blocks%idx_r_m(Ndim))
         allocate(blocks%idx_c_m(Ndim))
         call GetLocalBlockRange_MD(ptree, blocks%pgno, blocks%level_half+1, level_butterfly, Ndim, 1, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
         call assert(product(inc_r)==1,'inc_r has to be 1 for matrixblock_MD type')
         blocks%nr_m = nr
         blocks%idx_r_m = idx_r
         call GetLocalBlockRange_MD(ptree, blocks%pgno, blocks%level_half, level_butterfly, Ndim, 1, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
         call assert(product(inc_c)==1,'inc_c has to be 1 for matrixblock_MD type')
         blocks%nc_m = nc
         blocks%idx_c_m = idx_c

         allocate(blocks%ButterflyU(product(blocks%nr_m)))
         allocate(blocks%ButterflyV(product(blocks%nc_m)))
         if (level_butterfly /= 0) then
            allocate (blocks%ButterflyKerl_L(product(blocks%nr_m),blocks%level_half+1:level_butterfly))
            allocate (blocks%ButterflyKerl_R(product(blocks%nc_m),1:blocks%level_half))
         endif
         allocate (blocks%ButterflySkel_L(product(blocks%nr_m),blocks%level_half+1:level_butterfly+1))
         allocate (blocks%ButterflySkel_R(product(blocks%nc_m),0:blocks%level_half))
         allocate (ButterflySkel_R_transposed(product(blocks%nr_m)))

         do bb=1, product(blocks%nc_m)
         do level = 0, level_half
            nr=2**(level)
            nc=2**(blocks%level_half-level)
            idx_r=1
            idx_c=1
            inc_r=1
            inc_c=1

            if (level == 0) then
               blocks%ButterflyV(bb)%num_blk = nc(1)
               blocks%ButterflyV(bb)%nblk_loc = nc(1)
               blocks%ButterflyV(bb)%idx = idx_c(1)
               blocks%ButterflyV(bb)%inc = inc_c(1)
               allocate (blocks%ButterflyV(bb)%blocks(blocks%ButterflyV(bb)%nblk_loc,Ndim))
            elseif (level == level_butterfly + 1) then
               write(*,*)"should not reach level == level_butterfly + 1 for the right half of matrixblock_MD"
            else
               allocate(blocks%ButterflyKerl_R(bb,level)%num_col(Ndim))
               blocks%ButterflyKerl_R(bb,level)%num_col = nc*2
               allocate(blocks%ButterflyKerl_R(bb,level)%nc(Ndim))
               blocks%ButterflyKerl_R(bb,level)%nc = nc*2
               allocate(blocks%ButterflyKerl_R(bb,level)%num_row(Ndim))
               blocks%ButterflyKerl_R(bb,level)%num_row = nr
               allocate(blocks%ButterflyKerl_R(bb,level)%nr(Ndim))
               blocks%ButterflyKerl_R(bb,level)%nr = nr
               allocate(blocks%ButterflyKerl_R(bb,level)%idx_r(Ndim))
               blocks%ButterflyKerl_R(bb,level)%idx_r=idx_r
               allocate(blocks%ButterflyKerl_R(bb,level)%inc_r(Ndim))
               blocks%ButterflyKerl_R(bb,level)%inc_r=inc_r
               allocate(blocks%ButterflyKerl_R(bb,level)%idx_c(Ndim))
               blocks%ButterflyKerl_R(bb,level)%idx_c=idx_c
               allocate(blocks%ButterflyKerl_R(bb,level)%inc_c(Ndim))
               blocks%ButterflyKerl_R(bb,level)%inc_c=inc_c
               allocate (blocks%ButterflyKerl_R(bb,level)%blocks(product(blocks%ButterflyKerl_R(bb,level)%nr), blocks%ButterflyKerl_R(bb,level)%nc(1),Ndim))
            endif

            allocate(blocks%ButterflySkel_R(bb,level)%nc(Ndim))
            blocks%ButterflySkel_R(bb,level)%nc = nc
            allocate(blocks%ButterflySkel_R(bb,level)%nr(Ndim))
            blocks%ButterflySkel_R(bb,level)%nr = nr
            allocate(blocks%ButterflySkel_R(bb,level)%idx_r(Ndim))
            blocks%ButterflySkel_R(bb,level)%idx_r=idx_r
            allocate(blocks%ButterflySkel_R(bb,level)%inc_r(Ndim))
            blocks%ButterflySkel_R(bb,level)%inc_r=inc_r
            allocate(blocks%ButterflySkel_R(bb,level)%idx_c(Ndim))
            blocks%ButterflySkel_R(bb,level)%idx_c=idx_c
            allocate(blocks%ButterflySkel_R(bb,level)%inc_c(Ndim))
            blocks%ButterflySkel_R(bb,level)%inc_c=inc_c
            allocate (blocks%ButterflySkel_R(bb,level)%inds(product(blocks%ButterflySkel_R(bb,level)%nr), blocks%ButterflySkel_R(bb,level)%nc(1),Ndim))

            rank_new = 0
            flops = 0

            dim_MD(1:Ndim)=nr
            dim_MD(1+Ndim)=nc(1)
            dim_MD(2+Ndim)=Ndim

            allocate(subtensors(product(dim_MD)))
            do index_ij = 1, product(dim_MD)
               allocate(subtensors(index_ij)%nr(Ndim))
               allocate(subtensors(index_ij)%nc(Ndim))
               subtensors(index_ij)%nr=0
               subtensors(index_ij)%nc=0
            enddo
            n1 = MPI_Wtime()

            do index_ij = 1, product(dim_MD)
               call BF_MD_compress_N_oneblock_R_sample(Ndim, dim_MD, subtensors, blocks, boundary_map, Nboundall,Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_ij, bb, level, flops1)
            enddo
            n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2 - n1

            call element_Zmn_tensorlist_user(Ndim, subtensors, product(dim_MD), msh, option, ker, 0, passflag, ptree, stats)


            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp/max(1,blocks%level_butterfly/2)
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,rank_new1,flops1) reduction(MAX:rank_new) reduction(+:flops)
#endif
            do index_ij = 1, product(dim_MD)
               call BF_MD_compress_N_oneblock_R_rankreveal(Ndim, dim_MD, subtensors, blocks, option, stats, msh, ker, ptree, index_ij, bb, level, rank_new1, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp*max(1,blocks%level_butterfly/2)

            call BF_MD_delete_subtensors(Ndim, dim_MD, subtensors, stats)
            deallocate(subtensors)

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_tensorlist_user(Ndim, subtensors_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
            enddo
            ! call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (level < level_butterfly + 1) then
            if (rank_new > rankmax_for_butterfly(level)) then
               rankmax_for_butterfly(level) = rank_new
            endif
            endif
            if (level /= level_butterfly + 1) then
            if (level_half == level) then
               nc=2**(level_butterfly-level_half)
               call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb, idx_c)
               idx_c=idx_c + blocks%idx_c_m -1
               call MultiIndexToSingleIndex(Ndim, nc, bb_g, idx_c)
               if(option%verbosity>=1)write(*,*)"R: good until here!!",bb_g,product(nc),'rank_new',rank_new
               ! call BF_all2all_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'R', 'C')
            endif
            endif
            stats%Flop_Fill = stats%Flop_Fill + flops
         enddo
         enddo

         level_final = level_half + 1
         do bb=1, product(blocks%nr_m)
         do level = level_butterfly + 1, level_final, -1
            nr=2**(level - blocks%level_half - 1 )
            nc=2**(level_butterfly + 1-level)
            idx_r=1
            idx_c=1
            inc_r=1
            inc_c=1
            if (level == 0) then
               write(*,*)"should not reach level == 0 for the left half of matrixblock_MD"
            elseif (level == level_butterfly + 1) then
               blocks%ButterflyU(bb)%num_blk = nr(1)
               blocks%ButterflyU(bb)%nblk_loc = nr(1)
               blocks%ButterflyU(bb)%idx = idx_r(1)
               blocks%ButterflyU(bb)%inc = inc_r(1)
               allocate (blocks%ButterflyU(bb)%blocks(blocks%ButterflyU(bb)%nblk_loc,Ndim))
            else
               allocate(blocks%ButterflyKerl_L(bb,level)%num_col(Ndim))
               blocks%ButterflyKerl_L(bb,level)%num_col = nc
               allocate(blocks%ButterflyKerl_L(bb,level)%nc(Ndim))
               blocks%ButterflyKerl_L(bb,level)%nc = nc
               allocate(blocks%ButterflyKerl_L(bb,level)%num_row(Ndim))
               blocks%ButterflyKerl_L(bb,level)%num_row = nr*2
               allocate(blocks%ButterflyKerl_L(bb,level)%nr(Ndim))
               blocks%ButterflyKerl_L(bb,level)%nr = nr*2
               allocate(blocks%ButterflyKerl_L(bb,level)%idx_r(Ndim))
               blocks%ButterflyKerl_L(bb,level)%idx_r=idx_r
               allocate(blocks%ButterflyKerl_L(bb,level)%inc_r(Ndim))
               blocks%ButterflyKerl_L(bb,level)%inc_r=inc_r
               allocate(blocks%ButterflyKerl_L(bb,level)%idx_c(Ndim))
               blocks%ButterflyKerl_L(bb,level)%idx_c=idx_c
               allocate(blocks%ButterflyKerl_L(bb,level)%inc_c(Ndim))
               blocks%ButterflyKerl_L(bb,level)%inc_c=inc_c
               allocate (blocks%ButterflyKerl_L(bb,level)%blocks(blocks%ButterflyKerl_L(bb,level)%nr(1), product(blocks%ButterflyKerl_L(bb,level)%nc),Ndim))
            endif

            allocate(blocks%ButterflySkel_L(bb,level)%nc(Ndim))
            blocks%ButterflySkel_L(bb,level)%nc = nc
            allocate(blocks%ButterflySkel_L(bb,level)%nr(Ndim))
            blocks%ButterflySkel_L(bb,level)%nr = nr
            allocate(blocks%ButterflySkel_L(bb,level)%idx_r(Ndim))
            blocks%ButterflySkel_L(bb,level)%idx_r=idx_r
            allocate(blocks%ButterflySkel_L(bb,level)%inc_r(Ndim))
            blocks%ButterflySkel_L(bb,level)%inc_r=inc_r
            allocate(blocks%ButterflySkel_L(bb,level)%idx_c(Ndim))
            blocks%ButterflySkel_L(bb,level)%idx_c=idx_c
            allocate(blocks%ButterflySkel_L(bb,level)%inc_c(Ndim))
            blocks%ButterflySkel_L(bb,level)%inc_c=inc_c
            allocate (blocks%ButterflySkel_L(bb,level)%inds(blocks%ButterflySkel_L(bb,level)%nr(1), product(blocks%ButterflySkel_L(bb,level)%nc),Ndim))

            rank_new = 0
            flops = 0

            dim_MD(1)=nr(1)
            dim_MD(2:1+Ndim)=nc
            dim_MD(2+Ndim)=Ndim

            allocate(subtensors(product(dim_MD)))
            do index_ij = 1, product(dim_MD)
               allocate(subtensors(index_ij)%nr(Ndim))
               allocate(subtensors(index_ij)%nc(Ndim))
               subtensors(index_ij)%nr=0
               subtensors(index_ij)%nc=0
            enddo
            n1 = MPI_Wtime()

            do index_ij = 1, product(dim_MD)
               call BF_MD_compress_N_oneblock_C_sample(Ndim, dim_MD, subtensors, blocks, boundary_map, Nboundall,Ninadmissible, groupm_start, option, stats, msh, ker, ptree, index_ij, bb, level, flops1)
            enddo
            n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2 - n1

            call element_Zmn_tensorlist_user(Ndim, subtensors, product(dim_MD), msh, option, ker, 0, passflag, ptree, stats)


            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp/max(1,blocks%level_butterfly/2)
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,rank_new1,flops1) reduction(MAX:rank_new) reduction(+:flops)
#endif
            do index_ij = 1, product(dim_MD)
               call BF_MD_compress_N_oneblock_C_rankreveal(Ndim, dim_MD, subtensors, blocks, option, stats, msh, ker, ptree, index_ij, bb, level, rank_new1, flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            if(option%format==3 .and. option%near_para<=0.1d0)option%tol_comp = option%tol_comp*max(1,blocks%level_butterfly/2)


            call BF_MD_delete_subtensors(Ndim, dim_MD, subtensors, stats)
            deallocate(subtensors)

            passflag = 0
            do while (passflag == 0)
               call element_Zmn_tensorlist_user(Ndim, subtensors_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
            enddo

            ! call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (level > 0) then
            if (rank_new > rankmax_for_butterfly(level - 1)) then
               rankmax_for_butterfly(level - 1) = rank_new
            endif
            endif
            if (level_final == level) then

               nr=2**(level_half)
               call SingleIndexToMultiIndex(Ndim, blocks%nr_m, bb, idx_r)
               idx_r=idx_r + blocks%idx_r_m -1
               call MultiIndexToSingleIndex(Ndim, nr, bb_g, idx_r)
               if(option%verbosity>=1)write(*,*)"L: good until here!!",bb_g,product(nr),'rank_new',rank_new
               ! call BF_all2all_skel(blocks, blocks%ButterflySkel(level), option, stats, msh, ptree, level, 'R', 'C')
            endif

            stats%Flop_Fill = stats%Flop_Fill + flops
         enddo
         enddo

         call BF_MD_all2all_skel(Ndim, blocks, ButterflySkel_R_transposed, option, stats, msh, ptree)

         nc=2**(level_butterfly -level_half)
         dim_subtensor(1:Ndim)=blocks%nr_m
         dim_subtensor(1+Ndim:Ndim*2)=nc
         allocate(subtensors(product(dim_subtensor)))
         do index_ij = 1, product(dim_subtensor)
            call SingleIndexToMultiIndex(Ndim*2, dim_subtensor, index_ij, idx_subtensor)
            call MultiIndexToSingleIndex(Ndim, blocks%nr_m, idx_r_scalar, idx_subtensor(1:Ndim))
            call MultiIndexToSingleIndex(Ndim, nc, idx_c_scalar, idx_subtensor(1+Ndim:2*Ndim))

            group_m = blocks%row_group
            group_n = blocks%col_group
            group_m = group_m*2**(level_half) + blocks%idx_r_m - 1 + idx_subtensor(1:Ndim) -1
            group_n = group_n*2**(level_butterfly - level_half) + idx_subtensor(1+Ndim:2*Ndim) -1

            allocate(subtensors(index_ij)%nr(Ndim))
            allocate(subtensors(index_ij)%nc(Ndim))
            allocate(subtensors(index_ij)%rows(Ndim))
            allocate(subtensors(index_ij)%cols(Ndim))
            do dim_i=1,Ndim
               subtensors(index_ij)%nr(dim_i)=size(blocks%ButterflySkel_L(idx_r_scalar,level_final)%inds(1, idx_c_scalar, dim_i)%array)
               allocate (subtensors(index_ij)%rows(dim_i)%dat(subtensors(index_ij)%nr(dim_i)))
               subtensors(index_ij)%rows(dim_i)%dat = blocks%ButterflySkel_L(idx_r_scalar,level_final)%inds(1, idx_c_scalar, dim_i)%array + msh(dim_i)%basis_group(group_m(dim_i))%head - 1
            enddo

            do dim_i=1,Ndim
               subtensors(index_ij)%nc(dim_i)=size(ButterflySkel_R_transposed(idx_r_scalar)%inds(1, idx_c_scalar, dim_i)%array)
               allocate (subtensors(index_ij)%cols(dim_i)%dat(subtensors(index_ij)%nc(dim_i)))
               subtensors(index_ij)%cols(dim_i)%dat = ButterflySkel_R_transposed(idx_r_scalar)%inds(1, idx_c_scalar, dim_i)%array + msh(dim_i)%basis_group(group_n(dim_i))%head - 1
            enddo
            use_zfp=0
#if HAVE_ZFP
            if(option%use_zfp==1)use_zfp=1
#endif
            if(use_zfp==0)then
               allocate(subtensors(index_ij)%dat(product(subtensors(index_ij)%nr),product(subtensors(index_ij)%nc)))
               call LogMemory(stats, SIZEOF(subtensors(index_ij)%dat)/1024.0d3)
            endif
         enddo

         if(use_zfp==0)then
            call element_Zmn_tensorlist_user(Ndim, subtensors, product(dim_subtensor), msh, option, ker, 0, passflag, ptree, stats)
         else
            allocate(blocks%MiddleZFP(product(dim_subtensor)))
            call element_Zmn_tensorlist_user(Ndim, subtensors, product(dim_subtensor), msh, option, ker, 0, passflag, ptree, stats,blocks%MiddleZFP)
         endif

         allocate(blocks%ButterflyMiddle(product(dim_subtensor)))
         do index_ij = 1, product(dim_subtensor)
            blocks%ButterflyMiddle(index_ij)%matrix => subtensors(index_ij)%dat
            subtensors(index_ij)%dat=>null()
            allocate(blocks%ButterflyMiddle(index_ij)%dims_m(Ndim))
            blocks%ButterflyMiddle(index_ij)%dims_m = subtensors(index_ij)%nr
            allocate(blocks%ButterflyMiddle(index_ij)%dims_n(Ndim))
            blocks%ButterflyMiddle(index_ij)%dims_n = subtensors(index_ij)%nc
         enddo

         call BF_MD_delete_subtensors(Ndim, dim_subtensor, subtensors, stats)
         deallocate(subtensors)

         passflag = 0
         do while (passflag == 0)
            call element_Zmn_tensorlist_user(Ndim, subtensors_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
         enddo


      endif

      call MPI_ALLREDUCE(MPI_IN_PLACE, rankmax_for_butterfly(0:level_butterfly), level_butterfly + 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)

      if (statflag == 1) then
         if (allocated(stats%rankmax_of_level)) stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly), stats%rankmax_of_level(level_blocks))
      endif

      deallocate (rankmax_for_butterfly)
      deallocate (rankmin_for_butterfly)



      do bb=1, product(blocks%nc_m)
      do level = 0, blocks%level_half
      if (allocated(blocks%ButterflySkel_R(bb,level)%inds)) then
         do i = 1, size(blocks%ButterflySkel_R(bb,level)%inds, 1)
         do j = 1, size(blocks%ButterflySkel_R(bb,level)%inds, 2)
         do k = 1, Ndim
            if (allocated(blocks%ButterflySkel_R(bb,level)%inds(i, j,k)%array)) deallocate (blocks%ButterflySkel_R(bb,level)%inds(i, j,k)%array)
         enddo
         enddo
         enddo
         deallocate (blocks%ButterflySkel_R(bb,level)%inds)
      endif
      enddo
      enddo
      deallocate (blocks%ButterflySkel_R)


      do bb=1, product(blocks%nr_m)
      do level = blocks%level_half+1,level_butterfly+1
      if (allocated(blocks%ButterflySkel_L(bb,level)%inds)) then
         do i = 1, size(blocks%ButterflySkel_L(bb,level)%inds, 1)
         do j = 1, size(blocks%ButterflySkel_L(bb,level)%inds, 2)
         do k = 1, Ndim
            if (allocated(blocks%ButterflySkel_L(bb,level)%inds(i, j,k)%array)) deallocate (blocks%ButterflySkel_L(bb,level)%inds(i, j,k)%array)
         enddo
         enddo
         enddo
         deallocate (blocks%ButterflySkel_L(bb,level)%inds)
      endif
      enddo
      enddo
      deallocate (blocks%ButterflySkel_L)


      do bb=1, product(blocks%nr_m)
      if (allocated(ButterflySkel_R_transposed(bb)%inds)) then
         do i = 1, size(ButterflySkel_R_transposed(bb)%inds, 1)
         do j = 1, size(ButterflySkel_R_transposed(bb)%inds, 2)
         do k = 1, Ndim
            if (allocated(ButterflySkel_R_transposed(bb)%inds(i, j,k)%array)) deallocate (ButterflySkel_R_transposed(bb)%inds(i, j,k)%array)
         enddo
         enddo
         enddo
         deallocate (ButterflySkel_R_transposed(bb)%inds)
      endif
      enddo
      deallocate(ButterflySkel_R_transposed)


      call BF_MD_ComputeMemory(Ndim, blocks, memory_dense,memory_comp)
      Memory = memory_dense + memory_comp
      call BF_MD_get_rank(Ndim, blocks, ptree)
      if(option%verbosity>=0 .and. ptree%MyID==Main_ID)write(*,*)"After BF_MD_compress_N: myID",ptree%MyID, "rankmin",blocks%rankmin,"rankmax",blocks%rankmax,"memory_subtensor",memory_dense,"memory_factormat",memory_comp
      return

   end subroutine BF_MD_compress_N







   subroutine BF_compress_N15_seq(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)

      implicit none

      integer Nboundall, Ninadmissible, statflag
      integer boundary_map(:,:)
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
      DTR, allocatable :: Singular(:)
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
      type(intersect)::submats(1)
      type(SVD_quant)::SVD_Q

      level_butterfly = blocks%level_butterfly
      if(level_butterfly<=1)then
         call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)
      else

         cnt_tmp = 0
         rankFar = 0
         rankNear = 0
         maxvalue = 0
         minvalue = 10000
         Memory = 0

         ! write(*,*)blocks%row_group,blocks%col_group,'In BF_compress_N15'

         level_blocks = blocks%level
         level_butterfly = blocks%level_butterfly
         !level_butterfly=Maxlevel-level_blocks
         ! level_butterfly=int((Maxlevel_for_blocks-level_blocks)/2)*2

         blocks%rankmax = -100000
         blocks%rankmin = 100000
         call assert(option%pat_comp==3,'forwardN15flag=1 requires pat_comp=3')
         blocks%level_half = BF_Switchlevel(blocks%level_butterfly,option%pat_comp)
         levelm = blocks%level_half
         level_half = blocks%level_half

         num_blocks = 2**level_butterfly

         allocate (blocks%ButterflyU%blocks(num_blocks))
         allocate (blocks%ButterflyV%blocks(num_blocks))
         if (level_butterfly /= 0) then
            allocate (blocks%ButterflyKerl(level_butterfly))
         endif

         memory_butterfly = 0.

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

               rmax = min(option%BACA_Batch, min(mm, nn))
               allocate (SVD_Q%matU(mm, rmax))
               allocate (SVD_Q%matV(rmax, nn))
               allocate (SVD_Q%Singular(rmax))
               SVD_Q%matU=0
               SVD_Q%matV=0
               SVD_Q%Singular=0

               emptyflag = 0
               if (Nboundall > 0) then
                  do jj=1,Ninadmissible
                     if (boundary_map(group_m - groupm_start + 1,jj) == group_n) emptyflag = 1
                  enddo
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
                  call LR_ACA(SVD_Q, idxs_m, idxs_n, mm, nn, frow, rank, option%tol_comp*0.1, option%tol_comp, msh, ker, stats, ptree, option, error)
                  rankmax_for_butterfly(level_loc) = max(rank, rankmax_for_butterfly(level_loc))
                  rankmin_for_butterfly(level_loc) = min(rank, rankmin_for_butterfly(level_loc))

                  allocate (blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(rank, rank))
                  blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix = 0
                  do ii = 1, rank
                     blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix(ii, ii) = 1d0/SVD_Q%Singular(ii)
                  end do
               end if

               blocks%rankmax = max(blocks%rankmax, rank)
               blocks%rankmin = min(blocks%rankmin, rank)

               allocate (mat_tmp(mm, rank))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j,k,ctemp)
#endif
               do j = 1, rank
                  do i = 1, mm
                     mat_tmp(i, j) = SVD_Q%matU(i, j)*SVD_Q%Singular(j)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               mm1 = msh%basis_group(group_m*2)%tail - msh%basis_group(group_m*2)%head + 1
               allocate (ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix(mm1, rank))
               allocate (ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix(mm - mm1, rank))

               ButterflyP_old%blocks(2*index_i_loc - 1, index_j_loc)%matrix = mat_tmp(1:mm1, 1:rank)
               ButterflyP_old%blocks(2*index_i_loc, index_j_loc)%matrix = mat_tmp(1 + mm1:mm, 1:rank)

               deallocate (SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, mat_tmp)

            end do

            n1 = MPI_Wtime()

            do level_loc = 1, level_butterflyL
               level = level_loc + levelm

               if (level_loc /= level_butterflyL) then
                  allocate (ButterflyP%blocks(2**(level_loc + 1), 2**(level_butterflyL - level_loc)))
               end if

               ! do index_i_loc=1, 2**level_loc
               ! do index_j_loc=1, 2**(level_butterflyL-level_loc)

               ! write(*,*)'addaa111111'
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
#endif
               do index_ij_loc = 1, 2**level_butterflyL
                  index_j_loc = mod(index_ij_loc - 1, 2**(level_butterflyL - level_loc)) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**(level_butterflyL - level_loc)))

                  call LocalButterflySVD_Left(index_i_loc, index_j_loc, level_loc, level_butterflyL, level, index_i_m, blocks, option, msh, ButterflyP_old, ButterflyP)
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
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
            n2 = MPI_Wtime()
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

               rmax = min(option%BACA_Batch, min(mm, nn))
               allocate (SVD_Q%matU(mm, rmax))
               allocate (SVD_Q%matV(rmax, nn))
               allocate (SVD_Q%Singular(rmax))
               SVD_Q%matU=0
               SVD_Q%matV=0
               SVD_Q%Singular=0

               emptyflag = 0
               if (Nboundall > 0) then
                  do jj=1,Ninadmissible
                     if (boundary_map(group_m - groupm_start + 1,jj) == group_n) emptyflag = 1
                  enddo
               endif

               if (emptyflag == 1) then
                  ! if(.not. near_or_far(group_m,group_n,2d0))then
                  rank = 1
                  SVD_Q%Singular(1:rank) = 0
                  SVD_Q%matV(1:rank, 1:nn) = 0
                  ! if(blocks==342)write(111,*)Singular(1:rank)
               else
                  frow = 1
                  call LR_ACA(SVD_Q, idxs_m, idxs_n, mm, nn, frow, rank, option%tol_comp*0.1, option%tol_comp, msh, ker, stats, ptree, option, error)
               end if

               ! rank = min(rank,37)

               ! if(.not. near_or_far(group_m,group_n,2d0))then
               ! rank = 1
               ! Singular(1:rank) = 0
               ! ! if(blocks==342)write(111,*)Singular(1:rank)

               ! end if

               allocate (mat_tmp(rank, nn))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j,k,ctemp)
#endif
               do j = 1, nn
                  do i = 1, rank
                     mat_tmp(i, j) = SVD_Q%matV(i, j)*SVD_Q%Singular(i)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               nn1 = msh%basis_group(group_n*2)%tail - msh%basis_group(group_n*2)%head + 1
               allocate (ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix(rank, nn1))
               allocate (ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix(rank, nn - nn1))

               ButterflyP_old%blocks(index_i_loc, 2*index_j_loc - 1)%matrix = mat_tmp(1:rank, 1:nn1)
               ButterflyP_old%blocks(index_i_loc, 2*index_j_loc)%matrix = mat_tmp(1:rank, 1 + nn1:nn)

               deallocate (SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, mat_tmp)
            end do

            n1 = MPI_Wtime()
            do level_loc = 1, level_butterflyR
               level = levelm + 1 - level_loc

               if (level_loc /= level_butterflyR) then
                  allocate (ButterflyP%blocks(2**(level_butterflyR - level_loc), 2**(level_loc + 1)))
               end if
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(index_ij_loc,index_i_loc,index_j_loc)
#endif
               do index_ij_loc = 1, 2**level_butterflyR
                  index_j_loc = mod(index_ij_loc - 1, 2**level_loc) + 1
                  index_i_loc = ceiling_safe(dble(index_ij_loc)/dble(2**level_loc))
                  call LocalButterflySVD_Right(index_i_loc, index_j_loc, level_loc, level_butterflyR, level, level_butterfly, index_j_m, blocks, option, msh, ButterflyP_old, ButterflyP)
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
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
               ! call SVD_Truncate(QQ,mm,nn,mn,UU,VV,Singular,option%tol_comp,BPACK_SafeUnderflow,rank)
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

            n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2 - n1
         enddo


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

         deallocate (rankmax_for_butterfly)
         deallocate (rankmin_for_butterfly)

         Memory = memory_butterfly

         return



      endif

   end subroutine BF_compress_N15_seq


   subroutine BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)

      implicit none

      integer Nboundall, Ninadmissible, statflag
      integer boundary_map(:,:)
      integer groupm_start
      logical finish
      integer ranknew, rankup, Nqr, bsize, r_est, r_est_knn_r, r_est_knn_c, r_est_tmp, inc_c, inc_r, nc, nr, nrc, idx_c, idx_r
      integer blocks_idx, i, j, level_butterfly, num_blocks, k, attempt
      integer group_m, group_n, M,N, mm, nn, mn, index_i, index_j, index_i_m, index_j_m, index_ij_loc, index_i_loc, index_i_loc_s, index_j_loc, index_j_loc_s, ii, jj, iii, jjj, nn1, nn2, mm1, mm2, idxs_m, idxs_n, header_m, header_n
      integer level, levelm, level_half, level_final, level_loc, length_1, length_2, level_blocks, index_ij, edge_m, edge_n
      integer rank, rank0, rank_new1, rank_new, rankmax, butterflyB_inuse, rank1, rank2, rmax
      real(kind=8) rate, tolerance, SVD_tolerance, memory_butterfly, rtemp, norm_1, norm_2, norm_e
      real(kind=8) Memory, n1, n2
      DT ctemp
      type(butterfly_kerl) ButterflyMiddle, ButterflyP_old, ButterflyP_old1, ButterflyP_old2, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:), tau(:)
      integer, allocatable :: jpvt(:)
      DTR, allocatable :: Singular(:)
      integer flag, tmpi, tmpj
      real(kind=8):: a, b, error, flops, flops1, flop
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
      integer, allocatable :: rankmax_for_butterfly(:), mrange(:), nrange(:)
      integer emptyflag
      type(intersect),allocatable::submats(:),submats_full(:)
      type(acaquant),allocatable::acaquants(:)
      integer ierr
      DT, allocatable:: row_R(:, :), row_R_knn(:, :), row_Rtmp(:, :), row_Rtmp_knn(:, :), column_R(:, :), column_R_knn(:, :), column_RT(:, :), core_knn(:,:)
      integer, allocatable:: select_column(:), select_column_knn(:), select_row_knn(:), select_column1(:), select_row(:), perms(:)
      type(intersect) :: submats_dummy(1)
      type(SVD_quant)::SVD_Q
      integer :: fullmatflag

      ! write(*,*)'In: ',ptree%MyID,blocks%row_group,blocks%col_group
      flops=0

      level_butterfly = blocks%level_butterfly
      if(level_butterfly==0)then ! this calls LR_HBACA in BF_compress_NlogN
         call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, statflag)
      else
         Memory = 0
         level_blocks = blocks%level
         level_butterfly = blocks%level_butterfly

         call assert(option%pat_comp==3,'forwardN15flag=1 requires pat_comp=3')
         blocks%level_half = BF_Switchlevel(blocks%level_butterfly,option%pat_comp)
         level_half = blocks%level_half

         if (level_butterfly /= 0) then
            allocate (blocks%ButterflyKerl(level_butterfly))
         endif

         allocate (rankmax_for_butterfly(0:level_butterfly))
         rankmax_for_butterfly = 0

         call GetLocalBlockRange(ptree, blocks%pgno, level_half, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
         ButterflyMiddle%idx_r = idx_r
         ButterflyMiddle%inc_r = inc_r
         ButterflyMiddle%nr = nr
         ButterflyMiddle%idx_c = idx_c
         ButterflyMiddle%inc_c = inc_c
         ButterflyMiddle%nc = nc

         ButterflyP_old%idx_r = idx_r
         ButterflyP_old%inc_r = inc_r
         ButterflyP_old%nr = nr
         ButterflyP_old%idx_c = idx_c
         ButterflyP_old%inc_c = inc_c
         ButterflyP_old%nc = nc

         ButterflyP_old1%idx_r = idx_r
         ButterflyP_old1%inc_r = inc_r
         ButterflyP_old1%nr = nr
         ButterflyP_old1%idx_c = idx_c
         ButterflyP_old1%inc_c = inc_c
         ButterflyP_old1%nc = nc
         if (nc == 1 .and. 2**(level_butterfly - level_half) > 1) then ! double the number of local block columns used for MPI communication
            ButterflyMiddle%nc = 2
            ButterflyMiddle%idx_c = ButterflyMiddle%idx_c - 1 + mod(ButterflyMiddle%idx_c, 2)
            ButterflyP_old%nc = 2
            ButterflyP_old%idx_c = ButterflyP_old%idx_c - 1 + mod(ButterflyP_old%idx_c, 2)
         endif
         allocate (ButterflyMiddle%blocks(ButterflyMiddle%nr, ButterflyMiddle%nc))
         allocate (ButterflyP_old%blocks(ButterflyP_old%nr, ButterflyP_old%nc))
         allocate (ButterflyP_old1%blocks(ButterflyP_old1%nr, ButterflyP_old1%nc))




if(option%elem_extract>=1)then ! advancing multiple acas for entry extraction

         tolerance = option%tol_comp*0.1
         ! tolerance = option%tol_comp
         SVD_tolerance = option%tol_comp
         bsize = option%BACA_Batch

         if(allocated(stats%rankmax_of_level))then
            if(level_blocks>=1)then
            if(stats%rankmax_of_level(level_blocks-1)/=0)then
               bsize = max(64,ceiling_safe(stats%rankmax_of_level(level_blocks-1)/3d0))
            endif
            endif
         endif


         fullmatflag=0

         if(fullmatflag==1)then
            nrc=nr*nc
            allocate(submats_full(max(1,nrc)))
            do index_i_loc = 1, nr
               do index_j_loc = 1, nc
                  index_ij_loc = (index_j_loc-1)*nr+index_i_loc
                  index_i = (index_i_loc - 1)*inc_r + idx_r
                  index_j = (index_j_loc - 1)*inc_c + idx_c
                  group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
                  group_n = blocks%col_group
                  group_m = group_m*2**level_half - 1 + index_i
                  group_n = group_n*2**(level_butterfly - level_half) - 1 + index_j

                  M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                  N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                  header_m = msh%basis_group(group_m)%head
                  header_n = msh%basis_group(group_n)%head

                  submats_full(index_ij_loc)%nr=0
                  submats_full(index_ij_loc)%nc=0

                  emptyflag = 0
                  if (Nboundall > 0) then
                     do jj=1,Ninadmissible
                     if (boundary_map(group_m - groupm_start + 1,jj) == group_n) emptyflag = 1
                     enddo
                  endif

                  if (emptyflag == 1) then
                  else
                     submats_full(index_ij_loc)%nr=M
                     submats_full(index_ij_loc)%nc=N
                     allocate(submats_full(index_ij_loc)%rows(submats_full(index_ij_loc)%nr))
                     allocate(submats_full(index_ij_loc)%cols(submats_full(index_ij_loc)%nc))
                     do i=1,M
                        submats_full(index_ij_loc)%rows(i)=header_m + i - 1
                     enddo
                     do j = 1, N
                        submats_full(index_ij_loc)%cols(j) = header_n + j - 1
                     enddo
                     allocate(submats_full(index_ij_loc)%dat(submats_full(index_ij_loc)%nr,submats_full(index_ij_loc)%nc))
                     call LogMemory(stats, SIZEOF(submats_full(index_ij_loc)%dat)/1024.0d3)
                  endif
               enddo
            enddo
            call element_Zmn_blocklist_user(submats_full, nrc, msh, option, ker, 0, passflag, ptree, stats)
         endif


         n1 = MPI_Wtime()

         nrc=nr*nc
         allocate(submats(max(1,nrc*2)))  ! odd for columns, even for rows
         allocate(acaquants(nrc))
         do index_i_loc = 1, nr
            do index_j_loc = 1, nc
               index_ij_loc = (index_j_loc-1)*nr+index_i_loc
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
               group_n = blocks%col_group
               group_m = group_m*2**level_half - 1 + index_i
               group_n = group_n*2**(level_butterfly - level_half) - 1 + index_j

               M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
               N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
               header_m = msh%basis_group(group_m)%head
               header_n = msh%basis_group(group_n)%head
               acaquants(index_ij_loc)%M=M
               acaquants(index_ij_loc)%N=N
               acaquants(index_ij_loc)%header_m=header_m
               acaquants(index_ij_loc)%header_n=header_n
               submats(index_ij_loc*2-1)%nr=0
               submats(index_ij_loc*2-1)%nc=0
               submats(index_ij_loc*2)%nr=0
               submats(index_ij_loc*2)%nc=0

               emptyflag = 0
               if (Nboundall > 0) then
                  do jj=1,Ninadmissible
                     if (boundary_map(group_m - groupm_start + 1,jj) == group_n) emptyflag = 1
                  enddo
               endif

               if (emptyflag == 1) then
                  rank = 1
                  acaquants(index_ij_loc)%rank=rank
                  allocate(acaquants(index_ij_loc)%Singular(rank))
                  acaquants(index_ij_loc)%Singular = 1
                  allocate(acaquants(index_ij_loc)%matU(M,rank))
                  acaquants(index_ij_loc)%matU = 0
                  allocate(acaquants(index_ij_loc)%matV(rank,N))
                  acaquants(index_ij_loc)%matV = 0
                  acaquants(index_ij_loc)%finish = .true.
               else
                  r_est = min(bsize, min(M, N))
                  acaquants(index_ij_loc)%itrmax = floor_safe(min(M, N)/dble(r_est))*2
                  allocate(acaquants(index_ij_loc)%select_column(r_est))
                  acaquants(index_ij_loc)%select_column=0
                  allocate(acaquants(index_ij_loc)%select_row(r_est))
                  acaquants(index_ij_loc)%select_row=0
                  allocate(acaquants(index_ij_loc)%rows(M))
                  allocate(acaquants(index_ij_loc)%columns(N))

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
                     if (r_est_knn_r > 0 .and. r_est_knn_c > 0) then
                        submats(index_ij_loc*2-1)%nr=M
                        submats(index_ij_loc*2-1)%nc=r_est_knn_c
                        allocate(submats(index_ij_loc*2-1)%rows(submats(index_ij_loc*2-1)%nr))
                        allocate(submats(index_ij_loc*2-1)%cols(submats(index_ij_loc*2-1)%nc))
                        do i=1,M
                           submats(index_ij_loc*2-1)%rows(i)=header_m + i - 1
                        enddo
                        do j = 1, r_est_knn_c
                           submats(index_ij_loc*2-1)%cols(j) = header_n + select_column_knn(j) - 1
                        enddo
                        allocate(submats(index_ij_loc*2-1)%dat(submats(index_ij_loc*2-1)%nr,submats(index_ij_loc*2-1)%nc))
                        call LogMemory(stats, SIZEOF(submats(index_ij_loc*2-1)%dat)/1024.0d3)
                        if(fullmatflag==1)then
                           do i=1,M
                              do j = 1, r_est_knn_c
                                 submats(index_ij_loc*2-1)%dat(i,j) = submats_full(index_ij_loc)%dat(i,select_column_knn(j))
                              enddo
                           enddo
                        endif

                        submats(index_ij_loc*2)%nr=r_est_knn_r
                        submats(index_ij_loc*2)%nc=N
                        allocate(submats(index_ij_loc*2)%rows(submats(index_ij_loc*2)%nr))
                        allocate(submats(index_ij_loc*2)%cols(submats(index_ij_loc*2)%nc))
                        do i = 1, r_est_knn_r
                           submats(index_ij_loc*2)%rows(i) = header_m + select_row_knn(i) - 1
                        enddo
                        do j = 1, N
                           submats(index_ij_loc*2)%cols(j) = header_n + j - 1
                        enddo
                        allocate(submats(index_ij_loc*2)%dat(submats(index_ij_loc*2)%nr,submats(index_ij_loc*2)%nc))
                        call LogMemory(stats, SIZEOF(submats(index_ij_loc*2)%dat)/1024.0d3)
                        if(fullmatflag==1)then
                           do i=1,r_est_knn_r
                              do j = 1, N
                                 submats(index_ij_loc*2)%dat(i,j) = submats_full(index_ij_loc)%dat(select_row_knn(i),j)
                              enddo
                           enddo
                        endif
                     endif
                     deallocate(select_column_knn)
                     deallocate(select_row_knn)
                  endif
               endif
            enddo
         enddo

         ! write(*,*)'myid ',ptree%MyID,'before element_Zmn_blocklist_user'

         ! the rows and columns from KNN
         if(fullmatflag==0)call element_Zmn_blocklist_user(submats, nrc*2, msh, option, ker, 0, passflag, ptree, stats)

         ! write(*,*)'myid ',ptree%MyID,'after element_Zmn_blocklist_user'

         do index_i_loc = 1, nr
            do index_j_loc = 1, nc
               index_ij_loc = (index_j_loc-1)*nr+index_i_loc
               M = acaquants(index_ij_loc)%M
               N = acaquants(index_ij_loc)%N
               header_m = acaquants(index_ij_loc)%header_m
               header_n = acaquants(index_ij_loc)%header_n

               if(acaquants(index_ij_loc)%finish .eqv. .false.)then
                  r_est = min(bsize, min(M, N))
                  if (option%knn > 0) then

                     r_est_knn_r = submats(index_ij_loc*2)%nr
                     r_est_knn_c = submats(index_ij_loc*2-1)%nc

                     if (r_est_knn_r > 0 .and. r_est_knn_c > 0) then
                        allocate (row_R_knn(r_est_knn_r, N))
                        allocate (row_Rtmp_knn(r_est_knn_r, N))
                        allocate (core_knn(r_est_knn_r, r_est_knn_c))
                        allocate (column_R_knn(M, r_est_knn_c))

                        column_R_knn = submats(index_ij_loc*2-1)%dat
                        row_R_knn = submats(index_ij_loc*2)%dat
                        do i = 1, r_est_knn_r
                           core_knn(i, :) = column_R_knn(submats(index_ij_loc*2)%rows(i) - header_m + 1, :)
                        enddo

                        r_est = min(bsize, min(M, N))
                        Nqr = max(max(M, N), r_est)
                        allocate (jpvt(Nqr))
                        allocate (tau(Nqr))
                        jpvt = 0
                        call geqp3modf90(core_knn, jpvt, tau, tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
                        stats%Flop_Fill = stats%Flop_Fill + flop
                        rankup = ranknew
                        if (rankup > 0) then
                           rank = acaquants(index_ij_loc)%rank

                           row_Rtmp_knn = row_R_knn
                           call un_or_mqrf90(core_knn, tau, row_Rtmp_knn, 'L', 'C', r_est_knn_r, N, rankup, flop=flop)
                           stats%Flop_Fill = stats%Flop_Fill + flop
                           call trsmf90(core_knn, row_Rtmp_knn, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
                           stats%Flop_Fill = stats%Flop_Fill + flop

                           acaquants(index_ij_loc)%columns(rank + 1:rankup + rank) = submats(index_ij_loc*2-1)%cols(jpvt(1:rankup)) - header_n + 1
                           acaquants(index_ij_loc)%rows(rank + 1:rankup + rank) = submats(index_ij_loc*2)%rows(1:rankup) - header_m + 1

                           allocate(acaquants(index_ij_loc)%matU(M,rankup))
                           allocate(acaquants(index_ij_loc)%matV(rankup,N))
                           call LogMemory(stats, SIZEOF(acaquants(index_ij_loc)%matU)/1024.0d3 + SIZEOF(acaquants(index_ij_loc)%matV)/1024.0d3)

                           do j = 1, rankup
                              acaquants(index_ij_loc)%matU(:, rank + j) = column_R_knn(:, jpvt(j))
                           enddo
                           acaquants(index_ij_loc)%matV(rank + 1:rank + rankup, :) = row_Rtmp_knn(1:rankup, :)

                           rank = rank + rankup
                           rank0 = rank
                           acaquants(index_ij_loc)%rank = rank
                           acaquants(index_ij_loc)%rank0 = rank0

                           !>**** Find column pivots for the next iteration
                           jpvt = 0
                           row_Rtmp_knn = row_R_knn
                           if (rank > 0) row_Rtmp_knn(:, acaquants(index_ij_loc)%columns(1:rank)) = 0
                           ! call geqp3modf90(row_Rtmp_knn,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
                           call geqp3f90(row_Rtmp_knn, jpvt, tau, flop=flop)
                           stats%Flop_Fill = stats%Flop_Fill + flop
                           acaquants(index_ij_loc)%select_column(1:r_est) = jpvt(1:r_est)
                        endif
                        deallocate(row_R_knn)
                        deallocate(row_Rtmp_knn)
                        deallocate(core_knn)
                        deallocate(column_R_knn)
                        deallocate(jpvt)
                        deallocate(tau)
                     endif
                  endif
                  if(acaquants(index_ij_loc)%rank==0)then
                     allocate(perms(N))
                     call rperm(N, perms)
                     acaquants(index_ij_loc)%select_column = perms(1:r_est)
                     ! do i=1,r_est
                     !    acaquants(index_ij_loc)%select_column(i)=i
                     ! enddo
                     deallocate(perms)
                  endif
               endif
            enddo
         enddo

         if(fullmatflag==0)then
         passflag = 0
         do while (passflag == 0)
            call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
         enddo
         endif

         finish=.false.
         do while(.not. finish)
            do index_i_loc = 1, nr
               do index_j_loc = 1, nc
                  index_ij_loc = (index_j_loc-1)*nr+index_i_loc
                  M = acaquants(index_ij_loc)%M
                  N = acaquants(index_ij_loc)%N
                  header_m = acaquants(index_ij_loc)%header_m
                  header_n = acaquants(index_ij_loc)%header_n
                  r_est = min(bsize, min(M, N))

                  ! need to reset submats
                  submats(index_ij_loc*2-1)%nr=0
                  submats(index_ij_loc*2-1)%nc=0
                  if(associated(submats(index_ij_loc*2-1)%dat))then
                     call LogMemory(stats, -SIZEOF(submats(index_ij_loc*2-1)%dat)/1024.0d3)
                     deallocate(submats(index_ij_loc*2-1)%dat)
                     deallocate(submats(index_ij_loc*2-1)%rows)
                     deallocate(submats(index_ij_loc*2-1)%cols)
                  endif
                  submats(index_ij_loc*2)%nr=0
                  submats(index_ij_loc*2)%nc=0
                  if(associated(submats(index_ij_loc*2)%dat))then
                     call LogMemory(stats, -SIZEOF(submats(index_ij_loc*2)%dat)/1024.0d3)
                     deallocate(submats(index_ij_loc*2)%dat)
                     deallocate(submats(index_ij_loc*2)%rows)
                     deallocate(submats(index_ij_loc*2)%cols)
                  endif

                  if(acaquants(index_ij_loc)%finish .eqv. .false.)then
                     submats(index_ij_loc*2-1)%nr = M
                     submats(index_ij_loc*2-1)%nc = r_est
                     allocate(submats(index_ij_loc*2-1)%rows(submats(index_ij_loc*2-1)%nr))
                     do i=1,M
                        submats(index_ij_loc*2-1)%rows(i) = header_m + i - 1
                     enddo
                     allocate(submats(index_ij_loc*2-1)%cols(submats(index_ij_loc*2-1)%nc))
                     submats(index_ij_loc*2-1)%cols = acaquants(index_ij_loc)%select_column + header_n -1
                     allocate(submats(index_ij_loc*2-1)%dat(submats(index_ij_loc*2-1)%nr,submats(index_ij_loc*2-1)%nc))
                     call LogMemory(stats, SIZEOF(submats(index_ij_loc*2-1)%dat)/1024.0d3)
                     if(fullmatflag==1)then
                        do i=1,M
                           do j = 1, r_est
                              submats(index_ij_loc*2-1)%dat(i,j) = submats_full(index_ij_loc)%dat(i,acaquants(index_ij_loc)%select_column(j))
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo

            ! the columns
            if(fullmatflag==0)call element_Zmn_blocklist_user(submats, nrc*2, msh, option, ker, 0, passflag, ptree, stats)

#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij_loc,M,N,header_m,header_n,flops1) reduction(+:flops)
#endif
            do index_ij_loc = 1, nr*nc
               M = acaquants(index_ij_loc)%M
               N = acaquants(index_ij_loc)%N
               header_m = acaquants(index_ij_loc)%header_m
               header_n = acaquants(index_ij_loc)%header_n
               flops1 = 0
               if(acaquants(index_ij_loc)%finish .eqv. .false.)then
                  call LR_BACA_noOverlap_Oneiteration(header_m, header_n, M,N,acaquants(index_ij_loc),submats(index_ij_loc*2-1),submats(index_ij_loc*2), tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option,0, flops1)
                  ! write(*,*)acaquants(index_ij_loc)%itr, index_ij_loc,'col', acaquants(index_ij_loc)%finish, acaquants(index_ij_loc)%rank, acaquants(index_ij_loc)%normUV, acaquants(index_ij_loc)%normA
               endif
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            do index_i_loc = 1, nr
               do index_j_loc = 1, nc
                  index_ij_loc = (index_j_loc-1)*nr+index_i_loc
                  M = acaquants(index_ij_loc)%M
                  N = acaquants(index_ij_loc)%N
                  header_m = acaquants(index_ij_loc)%header_m
                  header_n = acaquants(index_ij_loc)%header_n
                  r_est = min(bsize, min(M, N))

                  ! need to reset submats, the column data is not deleted since it's still needed, just set nr=nc=0
                  submats(index_ij_loc*2-1)%nr=0
                  submats(index_ij_loc*2-1)%nc=0
                  ! if(associated(submats(index_ij_loc*2-1)%dat))then
                  !    deallocate(submats(index_ij_loc*2-1)%dat)
                  !    deallocate(submats(index_ij_loc*2-1)%rows)
                  !    deallocate(submats(index_ij_loc*2-1)%cols)
                  ! endif
                  submats(index_ij_loc*2)%nr=0
                  submats(index_ij_loc*2)%nc=0
                  if(associated(submats(index_ij_loc*2)%dat))then
                     call LogMemory(stats, -SIZEOF(submats(index_ij_loc*2)%dat)/1024.0d3)
                     deallocate(submats(index_ij_loc*2)%dat)
                     deallocate(submats(index_ij_loc*2)%rows)
                     deallocate(submats(index_ij_loc*2)%cols)
                  endif

                  if(acaquants(index_ij_loc)%finish .eqv. .false.)then
                     submats(index_ij_loc*2)%nr = r_est
                     submats(index_ij_loc*2)%nc = N
                     allocate(submats(index_ij_loc*2)%rows(submats(index_ij_loc*2)%nr))
                     submats(index_ij_loc*2)%rows = acaquants(index_ij_loc)%select_row + header_m -1
                     allocate(submats(index_ij_loc*2)%cols(submats(index_ij_loc*2)%nc))
                     do j=1,N
                        submats(index_ij_loc*2)%cols(j) = header_n + j - 1
                     enddo
                     allocate(submats(index_ij_loc*2)%dat(submats(index_ij_loc*2)%nr,submats(index_ij_loc*2)%nc))
                     call LogMemory(stats, SIZEOF(submats(index_ij_loc*2)%dat)/1024.0d3)
                     if(fullmatflag==1)then
                        do i=1,r_est
                           do j = 1, N
                              submats(index_ij_loc*2)%dat(i,j) = submats_full(index_ij_loc)%dat(acaquants(index_ij_loc)%select_row(i),j)
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
            ! the rows
            if(fullmatflag==0)call element_Zmn_blocklist_user(submats, nrc*2, msh, option, ker, 0, passflag, ptree, stats)
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij_loc,M,N,header_m,header_n,flops1) reduction(+:flops)
#endif
            do index_ij_loc = 1, nr*nc
               M = acaquants(index_ij_loc)%M
               N = acaquants(index_ij_loc)%N
               header_m = acaquants(index_ij_loc)%header_m
               header_n = acaquants(index_ij_loc)%header_n
               flops1 = 0
               if(acaquants(index_ij_loc)%finish .eqv. .false.)then
                  call LR_BACA_noOverlap_Oneiteration(header_m, header_n, M,N,acaquants(index_ij_loc),submats(index_ij_loc*2-1),submats(index_ij_loc*2), tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option,1,flops1)
                  ! write(*,*)acaquants(index_ij_loc)%itr, index_ij_loc,'row', acaquants(index_ij_loc)%finish, acaquants(index_ij_loc)%rank, acaquants(index_ij_loc)%normUV, acaquants(index_ij_loc)%normA
               endif
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            finish=.true.
            do index_i_loc = 1, nr
               do index_j_loc = 1, nc
                  index_ij_loc = (index_j_loc-1)*nr+index_i_loc
                  finish = finish .and. acaquants(index_ij_loc)%finish
               enddo
            enddo
         enddo

         do index_i_loc = 1, nr
            do index_j_loc = 1, nc
               index_ij_loc = (index_j_loc-1)*nr+index_i_loc
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               index_i_loc_s = (index_i - ButterflyP_old%idx_r)/ButterflyP_old%inc_r + 1
               index_j_loc_s = (index_j - ButterflyP_old%idx_c)/ButterflyP_old%inc_c + 1
               M = acaquants(index_ij_loc)%M
               N = acaquants(index_ij_loc)%N

               ! delete submats
               submats(index_ij_loc*2-1)%nr=0
               submats(index_ij_loc*2-1)%nc=0
               if(associated(submats(index_ij_loc*2-1)%dat))then
                  call LogMemory(stats, -SIZEOF(submats(index_ij_loc*2-1)%dat)/1024.0d3)
                  deallocate(submats(index_ij_loc*2-1)%dat)
                  deallocate(submats(index_ij_loc*2-1)%rows)
                  deallocate(submats(index_ij_loc*2-1)%cols)
               endif
               submats(index_ij_loc*2)%nr=0
               submats(index_ij_loc*2)%nc=0
               if(associated(submats(index_ij_loc*2)%dat))then
                  call LogMemory(stats, -SIZEOF(submats(index_ij_loc*2)%dat)/1024.0d3)
                  deallocate(submats(index_ij_loc*2)%dat)
                  deallocate(submats(index_ij_loc*2)%rows)
                  deallocate(submats(index_ij_loc*2)%cols)
               endif

               if(fullmatflag==1)then
                  submats_full(index_ij_loc)%nr=0
                  submats_full(index_ij_loc)%nc=0
                  if(associated(submats_full(index_ij_loc)%dat))then
                     call LogMemory(stats, -SIZEOF(submats_full(index_ij_loc)%dat)/1024.0d3)
                     deallocate(submats_full(index_ij_loc)%dat)
                     deallocate(submats_full(index_ij_loc)%rows)
                     deallocate(submats_full(index_ij_loc)%cols)
                  endif
               endif

               rank = acaquants(index_ij_loc)%rank
               rankmax_for_butterfly(level_half) = max(rank, rankmax_for_butterfly(level_half))
               allocate (ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank, rank))
               ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
               do ii = 1, rank
                  ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(ii, ii) = 1d0/acaquants(index_ij_loc)%Singular(ii)
               end do

               allocate (mat_tmp(M, rank))
               !!$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, rank
                  do i = 1, M
                     mat_tmp(i, j) = acaquants(index_ij_loc)%matU(i, j)*acaquants(index_ij_loc)%Singular(j)
                  enddo
               enddo
               !!$omp end parallel do
               allocate (ButterflyP_old%blocks(index_i_loc_s, index_j_loc_s)%matrix(M, rank))
               ButterflyP_old%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat_tmp(1:M, 1:rank)
               deallocate(mat_tmp)

               allocate (mat_tmp(rank, N))
               !!$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, N
                  do i = 1, rank
                     mat_tmp(i, j) = acaquants(index_ij_loc)%matV(i, j)*acaquants(index_ij_loc)%Singular(i)
                  enddo
               enddo
               !!$omp end parallel do
               allocate (ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix(rank, N))
               ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix = mat_tmp(1:rank, 1:N)
               deallocate (mat_tmp)

               if(allocated(acaquants(index_ij_loc)%Singular))deallocate(acaquants(index_ij_loc)%Singular)
               if(allocated(acaquants(index_ij_loc)%matU))then
                  call LogMemory(stats, -SIZEOF(acaquants(index_ij_loc)%matU)/1024.0d3)
                  deallocate(acaquants(index_ij_loc)%matU)
               endif
               if(allocated(acaquants(index_ij_loc)%matV))then
                  call LogMemory(stats, -SIZEOF(acaquants(index_ij_loc)%matV)/1024.0d3)
                  deallocate(acaquants(index_ij_loc)%matV)
               endif
               if(allocated(acaquants(index_ij_loc)%select_column))deallocate(acaquants(index_ij_loc)%select_column)
               if(allocated(acaquants(index_ij_loc)%select_row))deallocate(acaquants(index_ij_loc)%select_row)
               if(allocated(acaquants(index_ij_loc)%rows))deallocate(acaquants(index_ij_loc)%rows)
               if(allocated(acaquants(index_ij_loc)%columns))deallocate(acaquants(index_ij_loc)%columns)
            enddo
         enddo

         passflag = 0
         do while (passflag == 0)
            call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
         enddo
         deallocate(submats)
         deallocate(acaquants)

else
         allocate(submats(1))
         ! construct the middle level and the left half
         do index_i_loc = 1, nr
            do index_j_loc = 1, nc
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               index_i_loc_s = (index_i - ButterflyP_old%idx_r)/ButterflyP_old%inc_r + 1
               index_j_loc_s = (index_j - ButterflyP_old%idx_c)/ButterflyP_old%inc_c + 1

               group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
               group_n = blocks%col_group
               group_m = group_m*2**level_half - 1 + index_i
               group_n = group_n*2**(level_butterfly - level_half) - 1 + index_j

               mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
               nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
               idxs_m = msh%basis_group(group_m)%head
               idxs_n = msh%basis_group(group_n)%head

               rmax = min(option%BACA_Batch, min(mm, nn))
               allocate (SVD_Q%matU(mm, rmax))
               allocate (SVD_Q%matV(rmax, nn))
               allocate (SVD_Q%Singular(rmax))
               SVD_Q%matU=0
               SVD_Q%matV=0
               SVD_Q%Singular=0

               emptyflag = 0
               if (Nboundall > 0) then
                  do jj=1,Ninadmissible
                     if (boundary_map(group_m - groupm_start + 1,jj) == group_n) emptyflag = 1
                  enddo
               endif

               if (emptyflag == 1) then
                  rank = 1
                  SVD_Q%Singular(1:rank) = 0
                  SVD_Q%matU(1:mm, 1:rank) = 0
                  SVD_Q%matV(1:rank, 1:nn) = 0
                  rankmax_for_butterfly(level_half) = max(rank, rankmax_for_butterfly(level_half))
                  allocate (ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank, rank))
                  ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                  do ii = 1, rank
                     ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(ii, ii) = 0d0
                  end do
               else
                  if(option%RecLR_leaf == ACA)then
                     frow = 1
                     call LR_ACA(SVD_Q, idxs_m, idxs_n, mm, nn, frow, rank, option%tol_comp*0.1, option%tol_comp, msh, ker, stats, ptree, option, error)
                  elseif(option%RecLR_leaf == BACANOVER)then
                     call LR_BACA_noOverlap(SVD_Q, idxs_m, idxs_n, mm, nn, rank, option%tol_comp*0.1, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)
                  elseif(option%RecLR_leaf == BACA)then
                     call LR_BACA(SVD_Q, idxs_m, idxs_n, mm, nn, rank, option%tol_comp*0.1, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)
                  elseif (option%RecLR_leaf == RRQR) then
                        !!!!! RRQR
                        mn = min(mm, nn)
                        allocate (QQ(mm, nn))
                        allocate (mrange(mm))
                        allocate (nrange(nn))
                        do ii = 1, mm
                           mrange(ii) = idxs_m + ii - 1
                        enddo
                        do jj = 1, nn
                           nrange(jj) = idxs_n + jj - 1
                        enddo
                        submats(1)%nr = mm
                        submats(1)%nc = nn
                        allocate(submats(1)%rows(submats(1)%nr))
                        submats(1)%rows = mrange
                        allocate(submats(1)%cols(submats(1)%nc))
                        submats(1)%cols = nrange
                        allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
                        call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
                        QQ = submats(1)%dat
                        deallocate(submats(1)%rows)
                        deallocate(submats(1)%cols)
                        deallocate(submats(1)%dat)
                        deallocate (mrange)
                        deallocate (nrange)

                        call RRQR_SVD(QQ, mm, nn, mn, rmax, SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, option%tol_comp, rank, flops1)
                        stats%Flop_Fill = stats%Flop_Fill + flops1

                        deallocate(QQ)

                  else
                     write(*,*)'unsupported LR compression in BF_compress_N15, try RRQR, ACA, BACA, or BACANOVER instead'
                     stop
                  endif

                  rankmax_for_butterfly(level_half) = max(rank, rankmax_for_butterfly(level_half))
                  allocate (ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank, rank))
                  ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                  do ii = 1, rank
                     ButterflyMiddle%blocks(index_i_loc_s, index_j_loc_s)%matrix(ii, ii) = 1d0/SVD_Q%Singular(ii)
                  end do
               end if

               allocate (mat_tmp(mm, rank))
               !!$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, rank
                  do i = 1, mm
                     mat_tmp(i, j) = SVD_Q%matU(i, j)*SVD_Q%Singular(j)
                  enddo
               enddo
               !!$omp end parallel do
               allocate (ButterflyP_old%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, rank))
               ButterflyP_old%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat_tmp(1:mm, 1:rank)
               deallocate(mat_tmp)

               allocate (mat_tmp(rank, nn))
               !!$omp parallel do default(shared) private(i,j,k,ctemp)
               do j = 1, nn
                  do i = 1, rank
                     mat_tmp(i, j) = SVD_Q%matV(i, j)*SVD_Q%Singular(i)
                  enddo
               enddo
               !!$omp end parallel do
               allocate (ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix(rank, nn))
               ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix = mat_tmp(1:rank, 1:nn)
               deallocate (SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, mat_tmp)
            end do
         enddo
         passflag = 0
         do while (passflag == 0)
            call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 1, passflag, ptree, stats)
         enddo
         deallocate(submats)
endif



n2 = MPI_Wtime()
time_tmp = time_tmp + n2 - n1



         call BF_exchange_matvec(blocks, ButterflyP_old, stats, ptree, level_half, 'R', 'B')
         call BF_exchange_matvec(blocks, ButterflyMiddle, stats, ptree, level_half, 'R', 'B')
         call BF_all2all_vec_n_ker(blocks, ButterflyP_old1, stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'R', 'C', 0)


         !>*** needs to copy ButterflyP_old1 to ButterflyP_old2 with rows possibly doubled, before calling BF_exchange_matvec
         call GetLocalBlockRange(ptree, blocks%pgno, level_half+1, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
         ButterflyP_old2%idx_r = idx_r
         ButterflyP_old2%inc_r = inc_r
         ButterflyP_old2%nr = nr
         ButterflyP_old2%idx_c = idx_c
         ButterflyP_old2%inc_c = inc_c
         ButterflyP_old2%nc = nc

         if (nr == 1 .and. 2**(level_half) > 1) then ! double the number of local block rows used for MPI communication
            ButterflyP_old2%nr = 2
            ButterflyP_old2%idx_r = ButterflyP_old2%idx_r - 1 + mod(ButterflyP_old2%idx_r, 2)
         endif
         allocate (ButterflyP_old2%blocks(ButterflyP_old2%nr, ButterflyP_old2%nc))

         do index_i_loc = 1, nr
            do index_j_loc = 1, nc
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               index_i_loc_s = (index_i - ButterflyP_old2%idx_r)/ButterflyP_old2%inc_r + 1
               index_j_loc_s = (index_j - ButterflyP_old2%idx_c)/ButterflyP_old2%inc_c + 1

               mm = size(ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix, 1)
               nn = size(ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix, 2)
               allocate (ButterflyP_old2%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, nn))
               ButterflyP_old2%blocks(index_i_loc_s, index_j_loc_s)%matrix = ButterflyP_old1%blocks(index_i_loc, index_j_loc)%matrix
            enddo
         enddo
         if(level_half/=0)call BF_exchange_matvec(blocks, ButterflyP_old2, stats, ptree, level_half+1, 'C', 'B')
         call BF_delete_ker_onelevel(ButterflyP_old1)


         n1 = MPI_Wtime()
         ! construct the the left half
         do level = level_half+1, level_butterfly
            !>*** Caution!! the left half is row-wise, the right half is column-wise
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
            blocks%ButterflyKerl(level)%num_row = 2**level
            blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
            blocks%ButterflyKerl(level)%idx_r = idx_r
            blocks%ButterflyKerl(level)%inc_r = inc_r
            blocks%ButterflyKerl(level)%nr = nr
            blocks%ButterflyKerl(level)%idx_c = idx_c*2 - 1
            blocks%ButterflyKerl(level)%inc_c = inc_c
            blocks%ButterflyKerl(level)%nc = nc*2
            allocate (blocks%ButterflyKerl(level)%blocks(blocks%ButterflyKerl(level)%nr, blocks%ButterflyKerl(level)%nc))

            ButterflyP%idx_r = idx_r
            ButterflyP%inc_r = inc_r
            ButterflyP%nr = nr
            ButterflyP%idx_c = idx_c
            ButterflyP%inc_c = inc_c
            ButterflyP%nc = nc
            if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
               ButterflyP%nc = 2
               ButterflyP%idx_c = ButterflyP%idx_c - 1 + mod(ButterflyP%idx_c, 2)
            endif
            allocate (ButterflyP%blocks(ButterflyP%nr, ButterflyP%nc))

            rank_new=0
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,index_i, index_i_loc,index_j,index_j_loc,rank_new1,flops1) reduction(MAX:rank_new)  reduction(+:flops)
#endif
            do index_ij = 1, nr*nc
               index_i_loc = (index_ij - 1)/nc + 1
               index_j_loc = mod(index_ij - 1, nc) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               call ButterflySVD_Left(index_i, index_j, level, level_butterfly, level_half, blocks, option, msh, ButterflyP_old, ButterflyP, ButterflyMiddle, rank_new1,flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif

            call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (rank_new > rankmax_for_butterfly(level)) then
               rankmax_for_butterfly(level) = rank_new
            endif

            call BF_delete_ker_onelevel(ButterflyP_old)
            if(level/=level_butterfly)call BF_exchange_matvec(blocks, ButterflyP, stats, ptree, level, 'R', 'B')
            call BF_copy_ker_onelevel(ButterflyP, ButterflyP_old)
            call BF_delete_ker_onelevel(ButterflyP)
            call BF_all2all_vec_n_ker(blocks, blocks%ButterflyKerl(level), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level, 'R', 'C', 1)

         enddo
         n2 = MPI_Wtime()
         ! time_tmp = time_tmp + n2 - n1
         call BF_delete_ker_onelevel(ButterflyMiddle)

         !>**** need converting ButterflyP_old to column-wise layout before copy into ButterflyU
         call BF_all2all_vec_n_ker(blocks, ButterflyP_old, stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_butterfly+1, 'R', 'C', 0)
         blocks%ButterflyU%idx = ButterflyP_old%idx_r
         blocks%ButterflyU%inc = ButterflyP_old%inc_r
         blocks%ButterflyU%nblk_loc = ButterflyP_old%nr
         allocate (blocks%ButterflyU%blocks(blocks%ButterflyU%nblk_loc))
         do ii = 1, ButterflyP_old%nr
            mm = size(ButterflyP_old%blocks(ii, 1)%matrix, 1)
            nn = size(ButterflyP_old%blocks(ii, 1)%matrix, 2)
            allocate (blocks%ButterflyU%blocks(ii)%matrix(mm, nn))
            blocks%ButterflyU%blocks(ii)%matrix = ButterflyP_old%blocks(ii, 1)%matrix
            deallocate (ButterflyP_old%blocks(ii, 1)%matrix)
         end do
         deallocate (ButterflyP_old%blocks)

         ! construct the the right half
         n1 = MPI_Wtime()
         do level = level_half, 1, -1
            !>*** Caution!! the left half is row-wise, the right half is column-wise
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
            blocks%ButterflyKerl(level)%num_row = 2**level
            blocks%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
            blocks%ButterflyKerl(level)%idx_r = idx_r*2 - 1
            blocks%ButterflyKerl(level)%inc_r = inc_r
            blocks%ButterflyKerl(level)%nr = nr*2
            blocks%ButterflyKerl(level)%idx_c = idx_c
            blocks%ButterflyKerl(level)%inc_c = inc_c
            blocks%ButterflyKerl(level)%nc = nc
            allocate (blocks%ButterflyKerl(level)%blocks(blocks%ButterflyKerl(level)%nr, blocks%ButterflyKerl(level)%nc))

            ButterflyP%idx_r = idx_r
            ButterflyP%inc_r = inc_r
            ButterflyP%nr = nr
            ButterflyP%idx_c = idx_c
            ButterflyP%inc_c = inc_c
            ButterflyP%nc = nc
            if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
               ButterflyP%nr = 2
               ButterflyP%idx_r = ButterflyP%idx_r - 1 + mod(ButterflyP%idx_r, 2)
            endif
            allocate (ButterflyP%blocks(ButterflyP%nr, ButterflyP%nc))

            rank_new=0
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_ij,index_i,index_i_loc,index_j,index_j_loc,rank_new1,flops1) reduction(MAX:rank_new) reduction(+:flops)
#endif
            do index_ij = 1, nr*nc
               index_i_loc = (index_ij - 1)/nc + 1
               index_j_loc = mod(index_ij - 1, nc) + 1
               index_i = (index_i_loc - 1)*inc_r + idx_r
               index_j = (index_j_loc - 1)*inc_c + idx_c
               call ButterflySVD_Right(index_i, index_j, level, level_butterfly, blocks, option, msh, ButterflyP_old2, ButterflyP,rank_new1,flops1)
               rank_new = MAX(rank_new, rank_new1)
               flops =flops+flops1
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            call MPI_ALLREDUCE(MPI_IN_PLACE, rank_new, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if (rank_new > rankmax_for_butterfly(level)) then
               rankmax_for_butterfly(level) = rank_new
            endif

            call BF_delete_ker_onelevel(ButterflyP_old2)
            if(level/=1)call BF_exchange_matvec(blocks, ButterflyP, stats, ptree, level, 'C', 'B')
            call BF_copy_ker_onelevel(ButterflyP, ButterflyP_old2)
            call BF_delete_ker_onelevel(ButterflyP)
            call BF_all2all_vec_n_ker(blocks, blocks%ButterflyKerl(level), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level, 'C', 'R', 1)
         end do

         n2 = MPI_Wtime()
         ! time_tmp = time_tmp + n2 - n1


         !>**** need converting ButterflyP_old2 to row-wise layout before copy into ButterflyV
         call BF_all2all_vec_n_ker(blocks, ButterflyP_old2, stats, ptree, ptree%pgrp(blocks%pgno)%nproc, 0, 'C', 'R', 0)
         blocks%ButterflyV%idx = ButterflyP_old2%idx_c
         blocks%ButterflyV%inc = ButterflyP_old2%inc_c
         blocks%ButterflyV%nblk_loc = ButterflyP_old2%nc
         allocate (blocks%ButterflyV%blocks(blocks%ButterflyV%nblk_loc))
         do jj = 1, ButterflyP_old2%nc
            mm = size(ButterflyP_old2%blocks(1, jj)%matrix, 1)
            nn = size(ButterflyP_old2%blocks(1, jj)%matrix, 2)
            allocate (blocks%ButterflyV%blocks(jj)%matrix(nn, mm))
            do i=1,mm
            do j=1,nn
               blocks%ButterflyV%blocks(jj)%matrix(j,i)=ButterflyP_old2%blocks(1, jj)%matrix(i,j)
            enddo
            enddo
         end do


         call BF_delete_ker_onelevel(ButterflyP_old2)


         ! write(*,*)rankmax_for_butterfly,level_butterfly,blocks%level,Maxlevel_for_blocks
         if (statflag == 1) then
            if (allocated(stats%rankmax_of_level)) stats%rankmax_of_level(level_blocks) = max(maxval(rankmax_for_butterfly), stats%rankmax_of_level(level_blocks))
         endif
         ! write(*,*)stats%rankmax_of_level,'nitaima',rankmax_for_butterfly
         ! stop

         deallocate (rankmax_for_butterfly)


         call BF_ComputeMemory(blocks, Memory)
         call BF_get_rank(blocks, ptree)
         stats%Flop_Fill = stats%Flop_Fill + flops

         ! stop
         ! return
      endif

      ! write(*,*)'Out: ',ptree%MyID,blocks%row_group,blocks%col_group

   end subroutine BF_compress_N15



   subroutine BF_compress_test(blocks, msh, ker, element_Zmn, ptree, option, stats)



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
      ! distance_m=BPACK_Bigvalue
      ! allocate(distance_n(blocks%N))
      ! distance_n=BPACK_Bigvalue

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
            ! if(abs(value1)>BPACK_SafeUnderflow)write (*,*) abs(value1), abs(value2) !, abs(value1-value2)/abs(value1)
         enddo
      enddo

      deallocate (order_m)
      deallocate (order_n)

      write (*, *) 'partial fnorm:', v1, v2, sqrt(v3/v1), ' rank: ', blocks%rankmax, ' level: ', blocks%level_butterfly

      return

   end subroutine BF_compress_test

   subroutine LR_HBACA(blocks, leafsize, rank, option, msh, ker, stats, ptree, pgno, cridx)

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
      type(intersect)::submats(1)

      call LR_HBACA_Leaflevel(blocks, leafsize, rank, option, msh, ker, stats, ptree, pgno, cridx)

      passflag = 0
      do while (passflag == 0)

         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 1, passflag, ptree, stats)
         ! write(*,*)ptree%MyID,'after leaf',blocks%row_group,blocks%col_group,'passflag',passflag
         ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
      enddo
      stats%Flop_tmp=0
      call LR_HMerge(blocks, rank, option, msh, stats, ptree, pgno, cridx, 1)
      stats%Flop_Fill=stats%Flop_Fill+stats%Flop_tmp

   end subroutine LR_HBACA

   recursive subroutine LR_HMerge(blocks, rank, option, msh, stats, ptree, pgno, cridx, hbacaflag)

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
      DTR, allocatable::Singular(:)
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

            call blacs_gridinfo_wrp(ptree%pgrp(pgno1)%ctxt, nprow1, npcol1, myrow1, mycol1)
            call blacs_gridinfo_wrp(ptree%pgrp(pgno2)%ctxt, nprow2, npcol2, myrow2, mycol2)

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

            call blacs_gridinfo_wrp(ptree%pgrp(pgno)%ctxt, nprow, npcol, myrow, mycol)

            if (mod(cridx + 1, 2) == 0) then  ! merge along column dimension

               call assert(M1 == M2, 'M1/=M2 in column merge')

               if (nprow /= -1 .and. npcol /= -1) then
                  myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (matU(max(1,myArows), max(1,myAcols)))
                  matU = 0
                  call descinit_wp(descsmatU, M1, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU')

                  myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
                  allocate (matV1(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descsmatV1, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV1')

                  myArows = numroc_wp(N2, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
                  allocate (matV2(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descsmatV2, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV2')
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
                  call descinit_wp(descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU1c')

                  call pgemr2df90(M1, rank1, blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V1
                  myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit_wp(descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV1c')
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
                  call descinit_wp(descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU2c')
                  call pgemr2df90(M2, rank2, blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU, 1, 1 + rank1, descsmatU, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V2
                  myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit_wp(descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV2c')
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
                  allocate (UU(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descUU, M1, mnmin, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descUU')
                  UU = 0
                  myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (VV(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descVV, mnmin, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descVV')
                  VV = 0

                  allocate (Singular(mnmin))
                  Singular = 0

                  call PSVD_Truncate(M1, rank1 + rank2, matU, descsmatU, UU, VV, descUU, descVV, Singular, option%tol_comp, rank, ptree%pgrp(pgno)%ctxt, BPACK_SafeUnderflow, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)

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
                  allocate (blocks%ButterflyV%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
                  blocks%ButterflyV%blocks(1)%matrix = 0

                  call descinit_wp(descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descButterflyV')

                  myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyU%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descButterflyU')
                  blocks%ButterflyU%blocks(1)%matrix = UU(1:max(1,myArows), 1:max(1,myAcols))

                  call pgemmf90('N', 'T', N1, rank, rank1, BPACK_cone, matV1, 1, 1, descsmatV1, VV, 1, 1, descVV, BPACK_czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descButterflyV, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)
                  call pgemmf90('N', 'T', N2, rank, rank2, BPACK_cone, matV2, 1, 1, descsmatV2, VV, 1, 1 + rank1, descVV, BPACK_czero, blocks%ButterflyV%blocks(1)%matrix, 1 + N1, 1, descButterflyV, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)
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
                  allocate (matV(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descsmatV, N1, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV')

                  myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
                  allocate (matU1(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descsmatU1, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU1')

                  myArows = numroc_wp(M2, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
                  allocate (matU2(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descsmatU2, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU2')
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
                  call descinit_wp(descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU1c')
                  ! write(*,*)shape(blocks%sons(1,1)%ButterflyU%blocks(1)%matrix),shape(matU1),rank1,M1,blocks%sons(1,1)%rankmax
                  call pgemr2df90(M1, rank1, blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(1, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V1
                  myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
                  myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
                  call descinit_wp(descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno1)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV1c')
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
                  call descinit_wp(descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatU2c')
                  call pgemr2df90(M2, rank2, blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, ptree%pgrp(pgno)%ctxt1DCol)
                  deallocate (blocks%sons(2, 1)%ButterflyU%blocks(1)%matrix)

                  ! redistribute V2
                  myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
                  myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
                  call descinit_wp(descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno2)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descsmatV2c')
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
                  allocate (UU(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descUU, N1, mnmin, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descUU')
                  UU = 0
                  myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank1 + rank2, nbslpk, mycol, 0, npcol)
                  allocate (VV(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descVV, mnmin, rank1 + rank2, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descVV')
                  VV = 0

                  allocate (Singular(mnmin))
                  Singular = 0

                  call PSVD_Truncate(N1, rank1 + rank2, matV, descsmatV, UU, VV, descUU, descVV, Singular, option%tol_comp, rank, ptree%pgrp(pgno)%ctxt, BPACK_SafeUnderflow, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)
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
                  allocate (blocks%ButterflyV%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
                  blocks%ButterflyV%blocks(1)%matrix = UU(1:max(1,myArows), 1:max(1,myAcols))

                  call descinit_wp(descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descButterflyV')

                  myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
                  allocate (blocks%ButterflyU%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
                  call descinit_wp(descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ptree%pgrp(pgno)%ctxt, max(myArows, 1), info)
                  call assert(info == 0, 'descinit_wp fail for descButterflyU')
                  blocks%ButterflyU%blocks(1)%matrix = 0

                  call pgemmf90('N', 'T', M1, rank, rank1, BPACK_cone, matU1, 1, 1, descsmatU1, VV, 1, 1, descVV, BPACK_czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descButterflyU, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)

                  call pgemmf90('N', 'T', M2, rank, rank2, BPACK_cone, matU2, 1, 1, descsmatU2, VV, 1, 1 + rank1, descVV, BPACK_czero, blocks%ButterflyU%blocks(1)%matrix, 1 + M1, 1, descButterflyU, flop=flop)
                  stats%Flop_tmp = stats%Flop_tmp + flop/dble(nprow*npcol)
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
            call blacs_gridinfo_wrp(ptree%pgrp(pgno)%ctxt, nprow, npcol, myrow, mycol)
            if (myrow /= -1 .and. mycol /= -1) then
               myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
               ! if(myArows>0 .and. myAcols>0)then
               allocate (matU2D(max(1,myArows), max(1,myAcols)))
               matU2D = blocks%ButterflyU%blocks(1)%matrix
               ! endif
               myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
               ! if(myArows>0 .and. myAcols>0)then
               allocate (matV2D(max(1,myArows), max(1,myAcols)))
               matV2D = blocks%ButterflyV%blocks(1)%matrix
               ! endif
            else
               allocate (matU2D(1, 1))  ! required for Redistribute2Dto1D
               matU2D = 0
               allocate (matV2D(1, 1))  ! required for Redistribute2Dto1D
               matV2D = 0
            endif

            if (associated(blocks%ButterflyU%blocks(1)%matrix)) deallocate (blocks%ButterflyU%blocks(1)%matrix)
            allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc, rank))
            blocks%ButterflyU%blocks(1)%matrix = 0

            if (associated(blocks%ButterflyV%blocks(1)%matrix)) deallocate (blocks%ButterflyV%blocks(1)%matrix)
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
      DT, allocatable:: UU(:,:),VV(:,:),matU(:, :), matV(:, :), matU1(:, :), matV1(:, :), matU2(:, :), matV2(:, :), tmp(:, :), matU1D(:, :), matV1D(:, :), Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Fullmat(:, :), QQ1(:, :), matU2D(:, :), matV2D(:, :)
      DTR, allocatable::Singular(:)
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
      type(intersect)::submats(1)
      type(SVD_quant)::SVD_Q

      Maxgrp = 2**(ptree%nlevel) - 1

      if (option%RecLR_leaf == ACANMERGE) then
         frow = 1
         call LR_ACA_Parallel(blocks, blocks%headm, blocks%headn, blocks%M, blocks%N, frow, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, error, ptree%pgrp(pgno)%ctxt, pgno)
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

               rmaxc = min(option%rmax,blocks%N)
               rmaxr = min(blocks%M,option%rmax)
               ! rmaxc = blocks%N
               ! rmaxr = blocks%M


               call LR_SeudoSkeleton(blocks, blocks%headm, blocks%headn, blocks%M, blocks%N, rmaxc, rmaxr, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, ptree%pgrp(pgno)%ctxt)

            else if (option%RecLR_leaf == ACA) then
               !!!!! ACA-SVD
               rmax = min(option%BACA_Batch, min(blocks%M, blocks%N))
               allocate (SVD_Q%MatU(blocks%M, rmax))
               allocate (SVD_Q%MatV(rmax, blocks%N))
               allocate (SVD_Q%Singular(rmax))
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

               call LR_ACA(SVD_Q, blocks%headm, blocks%headn, blocks%M, blocks%N, frow, rank, option%tol_comp, option%tol_comp, msh, ker, stats, ptree, option, error)

               ! if(error>option%tol_comp)then
               ! write(*,*)'niam',error
               ! deallocate (UU,VV,Singular)
               ! goto 100
               ! endif

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = SVD_Q%MatV(j, i)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = SVD_Q%MatU(i, j)*SVD_Q%Singular(j)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               deallocate (SVD_Q%MatU, SVD_Q%MatV, SVD_Q%Singular)

            else if (option%RecLR_leaf == BACA .or. option%RecLR_leaf == BACANOVER) then
               !!! blocked ACA

               rmax = min(option%BACA_Batch, min(blocks%M, blocks%N))
               allocate (SVD_Q%Singular(rmax))
               allocate (SVD_Q%MatU(blocks%M, rmax))
               allocate (SVD_Q%MatV(rmax, blocks%N))

               if (option%RecLR_leaf == BACA) call LR_BACA(SVD_Q, blocks%headm, blocks%headn, blocks%M, blocks%N, rank, option%tol_comp, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)
               if (option%RecLR_leaf == BACANOVER) call LR_BACA_noOverlap(SVD_Q, blocks%headm, blocks%headn, blocks%M, blocks%N, rank, option%tol_comp, option%tol_comp, option%BACA_Batch, msh, ker, stats, ptree, option, error)

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = SVD_Q%MatV(j, i)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = SVD_Q%MatU(i, j)*SVD_Q%Singular(j)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               deallocate (SVD_Q%MatU, SVD_Q%MatV, SVD_Q%Singular)

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

               ! call element_Zmn_block_user(blocks%M, blocks%N, mrange, nrange, QQ1, msh, option, ker, 0, passflag, ptree, stats)

               submats(1)%nr = blocks%M
               submats(1)%nc = blocks%N
               allocate(submats(1)%rows(submats(1)%nr))
               submats(1)%rows = mrange
               allocate(submats(1)%cols(submats(1)%nc))
               submats(1)%cols = nrange
               allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
               call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
               QQ1 = submats(1)%dat
               deallocate(submats(1)%rows)
               deallocate(submats(1)%cols)
               deallocate(submats(1)%dat)



               deallocate (mrange)
               deallocate (nrange)

               ! do ii=1,blocks%M
               ! do jj =1,blocks%N
               ! edge_m = blocks%headm +ii - 1
               ! edge_n = blocks%headn +jj - 1
               ! call element_Zmn(edge_m,edge_n,QQ1(ii,jj),msh,option,ker)
               ! end do
               ! end do

               call SVD_Truncate(QQ1, blocks%M, blocks%N, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop

               ! rank=blocks%N

               blocks%rankmax = rank
               blocks%rankmin = rank

               allocate (blocks%ButterflyV%blocks(1)%matrix(blocks%N, rank))

               ! blocks%ButterflyV%blocks(1)%matrix=0
               ! do j=1, rank
               ! blocks%ButterflyV%blocks(1)%matrix(j,j)=1d0
               ! enddo
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%N
                     blocks%ButterflyV%blocks(1)%matrix(i, j) = VV(j, i)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               allocate (blocks%ButterflyU%blocks(1)%matrix(blocks%M, rank))
               ! blocks%ButterflyU%blocks(1)%matrix = QQ1
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(i,j)
#endif
               do j = 1, rank
                  do i = 1, blocks%M
                     blocks%ButterflyU%blocks(1)%matrix(i, j) = UU(i, j)*Singular(j)
                  enddo
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
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

               submats(1)%nr = blocks%M
               submats(1)%nc = blocks%N
               allocate(submats(1)%rows(submats(1)%nr))
               submats(1)%rows = mrange
               allocate(submats(1)%cols(submats(1)%nc))
               submats(1)%cols = nrange
               allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
               call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
               QQ1 = submats(1)%dat
               deallocate(submats(1)%rows)
               deallocate(submats(1)%cols)
               deallocate(submats(1)%dat)

               ! call element_Zmn_block_user(blocks%M, blocks%N, mrange, nrange, QQ1, msh, option, ker, 0, passflag, ptree, stats)

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

   subroutine LR_ACA(SVD_Q,header_m, header_n, rankmax_r, rankmax_c, frow, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, error)


      implicit none

      type(SVD_quant)::SVD_Q
      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance, dist
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n, Dimn, mn, Navr, itr
      integer rank, ranknew, row, column, rankmax, rankmax_c, rankmax_r, rankmax_min, rmax, rmax0,idxs_r, idxs_c, frow
      DT value_Z, maxvalue
      DT inner_U, inner_V, ctemp, value_UVs
      real(kind=8) inner_UV, n1, n2, a, error, flop
      integer, allocatable:: select_column(:), select_row(:)
      DT,allocatable::matU(:, :), matV(:, :)
      DT::matr(1, rankmax_c), matc(rankmax_r, 1)
      DTR,allocatable::Singular(:)

      DT, allocatable:: row_R(:), column_R(:), value_UV(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:), norm_UVavrbynorm_Z(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:)
      integer::mrange(rankmax_r), nrange(rankmax_c)
      type(Hstat)::stats
      type(Hoption)::option
      integer:: passflag = 0
      type(intersect)::submats(1)

      Navr = 3 !5 !10
      itr = 1
      allocate (norm_UVavrbynorm_Z(Navr))
      norm_UVavrbynorm_Z = 0

      n1 = MPI_Wtime()

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

      submats(1)%nr = 1
      submats(1)%nc = rankmax_c
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      matr = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)

      ! call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)

      rmax=size(SVD_Q%MatU,2)

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

      if (abs(maxvalue) < BPACK_SafeUnderflow) then

         do ii = 1, 100
            a = 0
            call random_number(a)
            select_row(1) = floor_safe(a*(rankmax_r - 1)) + 1

            mrange = select_row(1)
            do j = 1, rankmax_c
               nrange(j) = header_n + j - 1
            enddo

            submats(1)%nr = 1
            submats(1)%nc = rankmax_c
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            matr = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)


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
            if (abs(maxvalue) > BPACK_SafeUnderflow) exit
         end do
         if (abs(maxvalue) < BPACK_SafeUnderflow) then
            rank = 1
            SVD_Q%matU(:, 1) = 0
            SVD_Q%matV(1, :) = 0
            SVD_Q%Singular(1) = 0

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
         SVD_Q%matV(1, j) = row_R(j)
      enddo
      ! !$omp end parallel do

      nrange = header_n + select_column(1) - 1
      do i = 1, rankmax_r
         mrange(i) = header_m + i - 1
      enddo

      submats(1)%nr = rankmax_r
      submats(1)%nc = 1
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      matc = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)


      ! call element_Zmn_block_user(rankmax_r, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)



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
         SVD_Q%matU(i, 1) = column_R(i)
      enddo
      ! !$omp end parallel do

      norm_U = norm_vector(column_R, rankmax_r)
      norm_V = norm_vector(row_R, rankmax_c)
      norm_Z = norm_Z + norm_U*norm_V
      if (norm_Z > BPACK_SafeUnderflow) then
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

         submats(1)%nr = 1
         submats(1)%nc = rankmax_c
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         matr = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(1, rankmax_c, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)

         row_R = matr(1, :)

         ! !$omp parallel do default(shared) private(j,i,value_Z,edge_m,edge_n)
         ! do j=1,rankmax_c
         ! edge_m = header_m + select_row(rank+1) - 1
         ! edge_n = header_n + j - 1
         ! call element_Zmn(edge_m,edge_n,row_R(j),msh,option,ker)
         ! enddo
         ! !$omp end parallel do

         call gemmf77('N', 'N', 1, rankmax_c, rank, BPACK_cone, SVD_Q%matU(select_row(rank + 1), 1), rankmax_r, SVD_Q%matV, rmax, BPACK_czero, value_UV, 1)

         ! call gemmf90(SVD_Q%matU(select_row(rank+1),1), rankmax_r,SVD_Q%matV,rmax,value_UV,1,'N','N',1,rankmax_c,rank,BPACK_cone,BPACK_czero)

         row_R = row_R - value_UV(1:rankmax_c)
         norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         if (abs(maxvalue) < BPACK_SafeUnderflow) then
            ! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
            exit
            row_R = 0
         else
            do j = 1, rankmax_c
               row_R(j) = row_R(j)/maxvalue
            enddo
         endif
         ! !$omp parallel do default(shared) private(j)

         if (rank+1 > rmax) then

            allocate(matU(rankmax_r,rmax))
            allocate(matV(rmax,rankmax_c))
            allocate(Singular(rmax))
            matU=SVD_Q%matU
            matV=SVD_Q%matV
            Singular=SVD_Q%Singular

            rmax0=rmax
            rmax = min(2*rmax,rankmax_min)
            deallocate(SVD_Q%matU,SVD_Q%matV,SVD_Q%Singular)
            allocate(SVD_Q%matU(rankmax_r,rmax))
            allocate(SVD_Q%matV(rmax,rankmax_c))
            allocate(SVD_Q%Singular(rmax))
            SVD_Q%matU=0
            SVD_Q%matV=0
            SVD_Q%Singular=0
            SVD_Q%matU(:,1:rmax0) = matU(:,1:rmax0)
            SVD_Q%matV(1:rmax0,:) = matV(1:rmax0,:)
            SVD_Q%Singular(1:rmax0) = Singular(1:rmax0)
            deallocate(matU)
            deallocate(matV)
            deallocate(Singular)
         end if


         ! !$omp end parallel do
         ! !$omp parallel do default(shared) private(j)
         do j = 1, rankmax_c
            SVD_Q%matV(rank + 1, j) = row_R(j)
         enddo
         ! !$omp end parallel do

         nrange(1) = header_n + select_column(rank + 1) - 1
         do i = 1, rankmax_r
            mrange(i) = header_m + i - 1
         enddo


         submats(1)%nr = rankmax_r
         submats(1)%nc = 1
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         matc = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(rankmax_r, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)


         column_R = matc(:, 1)

         ! !$omp parallel do default(shared) private(i,j,value_Z,value_UVs,edge_m,edge_n)
         ! do i=1,rankmax_r
         ! edge_m = header_m + i - 1
         ! edge_n = header_n + select_column(rank+1) - 1
         ! call element_Zmn(edge_m,edge_n,column_R(i),msh,option,ker)
         ! enddo
         ! !$omp end parallel do

         call gemmf77('N', 'N', rankmax_r, 1, rank, BPACK_cone, SVD_Q%matU, rankmax_r, SVD_Q%matV(1, select_column(rank + 1)), rmax, BPACK_czero, value_UV, rankmax_r)

         ! call gemmf90(matU, rankmax_r,SVD_Q%matV(1,select_column(rank+1)),rmax,value_UV,rankmax_r,'N','N',rankmax_r,1,rank,BPACK_cone,BPACK_czero)

         column_R = column_R - value_UV(1:rankmax_r)
         norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_r

         do i = 1, rank + 1
            norm_column_R(select_row(i)) = 0
         enddo

         ! !$omp parallel do default(shared) private(i)
         do i = 1, rankmax_r
            SVD_Q%matU(i, rank + 1) = column_R(i)
         enddo
         ! !$omp end parallel do

         norm_U = norm_vector(column_R, rankmax_r)
         norm_V = norm_vector(row_R, rankmax_c)

         inner_UV = 0
#ifdef HAVE_OPENMP
         !$omp parallel do default(shared) private(j,i,ctemp,inner_V,inner_U) reduction(+:inner_UV)
#endif
         do j = 1, rank
            inner_U = 0
            inner_V = 0
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i = 1, rankmax_r
               ctemp = SVD_Q%matU(i, rank + 1)*conjg(cmplx(SVD_Q%matU(i, j), kind=8))
               inner_U = inner_U + ctemp
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i = 1, rankmax_c
               ctemp = SVD_Q%matV(rank + 1, i)*conjg(cmplx(SVD_Q%matV(j, i), kind=8))
               inner_V = inner_V + ctemp
            enddo
            ! !$omp end parallel do
            inner_UV = inner_UV + 2*dble(inner_U*inner_V)
         enddo
#ifdef HAVE_OPENMP
         !$omp end parallel do
#endif
         norm_Z = norm_Z + inner_UV + norm_U*norm_V

         ! ! write(*,*)norm_Z,inner_UV,norm_U,norm_V,maxvalue,rank,'gan'
         ! if(myisnan(sqrt(norm_Z)))then
         ! write(*,*)inner_UV,norm_U,norm_V,maxvalue
         ! stop
         ! endif

         if (norm_Z > BPACK_SafeUnderflow) then
            norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
         else
            norm_UVavrbynorm_Z(itr) = 0
         endif
         itr = mod(itr, Navr) + 1

         stats%Flop_Fill = stats%Flop_Fill + rank*rankmax_c + rank*rankmax_r

         rank = rank + 1

         if(rankmax_min<=rmax)then
            exit
         endif
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

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

! ACA followed by SVD

      allocate (QQ1(rankmax_r, rank))
      QQ1 = SVD_Q%matU(1:rankmax_r, 1:rank)
      ! call copymatN(SVD_Q%matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
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
      call copymatT(SVD_Q%matV(1:rank, 1:rankmax_c), QQ2, rank, rankmax_c)
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
      call gemmf90(RR1, rank, RR2, rank, mattemp, rank, 'N', 'T', rank, rank, rank, BPACK_cone, BPACK_czero, flop=flop)
      ! call zgemm('N','T',rank,rank,rank, BPACK_cone, RR1, rank,RR2,rank,BPACK_czero,mattemp,rank)
      stats%Flop_Fill = stats%Flop_Fill + flop
      allocate (UUsml(rank, rank), VVsml(rank, rank), Singularsml(rank))
      call SVD_Truncate(mattemp, rank, rank, rank, UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop
      ! call zgemm('N','N',rankmax_r,ranknew,rank, BPACK_cone, QQ1, rankmax_r,UUsml,rank,BPACK_czero,SVD_Q%matU,rankmax_r)
      call gemmf90(QQ1, rankmax_r, UUsml, rank, SVD_Q%matU, rankmax_r, 'N', 'N', rankmax_r, ranknew, rank, BPACK_cone, BPACK_czero, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop
      ! call zgemm('N','T',ranknew,rankmax_c,rank, BPACK_cone, VVsml, rank,QQ2,rankmax_c,BPACK_czero,SVD_Q%matV,rmax)
      call gemmf90(VVsml, rank, QQ2, rankmax_c, SVD_Q%matV, rmax, 'N', 'T', ranknew, rankmax_c, rank, BPACK_cone, BPACK_czero, flop=flop)
      stats%Flop_Fill = stats%Flop_Fill + flop

      rank = ranknew
      SVD_Q%Singular(1:ranknew) = Singularsml(1:ranknew)

      deallocate (mattemp, RR1, QQ1, UUsml, VVsml, Singularsml)
      deallocate (QQ2, RR2)

      deallocate (select_column)
      deallocate (select_row)
      deallocate (value_UV)

      return

   end subroutine LR_ACA

   subroutine LR_ACA_Parallel(blocks, header_m, header_n, M, N, frow, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, error, ctxt, pgno)


      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      type(matrixblock)::blocks
      integer index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance, dist
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n, Dimn, mn, Navr, itr
      integer rank, rank1, rank2, ranknew, row, column, rankmax, N, M, rankmax_min, rmax0, rmax, idxs_r, idxs_c, frow, mn1, mn2
      DT value_Z, maxvalue
      DT inner_U, inner_V, ctemp, value_UVs
      real(kind=8) inner_UV, n1, n2, a, error, flop
      integer:: select_column(N), select_row(M)
      DT, allocatable::matU(:, :), matV(:, :),matU0(:, :), matV0(:, :), matU2D(:, :), matV2D(:, :)
      DT,allocatable::tmpval(:)
      DT::matr(1, blocks%N_loc), matc(blocks%M_loc, 1)
      DT, allocatable:: value_UV(:)
      DT:: row_R(N), column_R(M)
      real(kind=8):: norm_row_R(N), norm_column_R(M)
      real(kind=8), allocatable:: norm_UVavrbynorm_Z(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), UUu(:, :), UUv(:, :), VVu(:, :), VVv(:, :)
      DTR, allocatable :: Singularsml(:), Singularuv(:)
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
      type(intersect)::submats(1)

      pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
      headm_loc = blocks%M_p(pp, 1)
      headn_loc = blocks%N_p(pp, 1)

      rmax = option%BACA_Batch

      allocate (matU(blocks%M_loc, rmax))
      allocate (matV(rmax, blocks%N_loc))
      allocate (tmpval(rmax))

      Navr = 3 !5 !10
      itr = 1
      allocate (norm_UVavrbynorm_Z(Navr))
      norm_UVavrbynorm_Z = 0

      n1 = MPI_Wtime()

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

      submats(1)%nr = 1
      submats(1)%nc = blocks%N_loc
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      matr = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)

      ! call element_Zmn_block_user(1, blocks%N_loc, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)


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

      submats(1)%nr = blocks%M_loc
      submats(1)%nc = 1
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      matc = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)


      ! call element_Zmn_block_user(blocks%M_loc, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)

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
      if (norm_Z > BPACK_SafeUnderflow) then
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

         submats(1)%nr = 1
         submats(1)%nc = blocks%N_loc
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         matr = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(1, blocks%N_loc, mrange, nrange, matr, msh, option, ker, 0, passflag, ptree, stats)

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

         call gemmf77('N', 'N', 1, blocks%N_loc, rank, BPACK_cone, tmpval(1), 1, matV, rmax, BPACK_czero, value_UV(headn_loc), 1)
         call MPI_ALLREDUCE(MPI_IN_PLACE, value_UV, N, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         ! call gemmf90(matU(select_row(rank+1),1), M,matV,rmax,value_UV,1,'N','N',1,N,rank,BPACK_cone,BPACK_czero)

         row_R = row_R - value_UV(1:N)
         norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))
         stats%Flop_Fill = stats%Flop_Fill + rank*blocks%N_loc

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         if (abs(maxvalue) < BPACK_SafeUnderflow) then
            ! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
            exit
            row_R = 0
         else
            do j = 1, N
               row_R(j) = row_R(j)/maxvalue
            enddo
         endif
         ! !$omp parallel do default(shared) private(j)


         if (rank+1 > rmax) then
            rmax0=rmax
            rmax = rmax*2
            deallocate(tmpval)
            allocate(tmpval(rmax))

            allocate(matU0(blocks%M_loc,rmax0))
            allocate(matV0(rmax0,blocks%N_loc))
            matU0=matU
            matV0=matV
            deallocate(matU)
            deallocate(matV)
            allocate(matU(blocks%M_loc,rmax))
            allocate(matV(rmax,blocks%N_loc))
            matU(:,1:rmax0)=matU0(:,1:rmax0)
            matV(1:rmax0,:)=matV0(1:rmax0,:)
            deallocate(matU0)
            deallocate(matV0)

            ! write(*,*)'increase rmax',rank,rmax
            ! stop
            ! exit
         end if

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

         submats(1)%nr = blocks%M_loc
         submats(1)%nc = 1
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         matc = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(blocks%M_loc, 1, mrange, nrange, matc, msh, option, ker, 0, passflag, ptree, stats)
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
         call gemmf77('N', 'N', blocks%M_loc, 1, rank, BPACK_cone, matU, blocks%M_loc, tmpval(1), rank, BPACK_czero, value_UV(headm_loc), M)

         ! write(*,*)'dd4',ptree%MyID,sum(value_UV(1:M))

         call MPI_ALLREDUCE(MPI_IN_PLACE, value_UV, M, MPI_DT, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)

         ! call gemmf90(matU, M,matV(1,select_column(rank+1)),rmax,value_UV,M,'N','N',M,1,rank,BPACK_cone,BPACK_czero)

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

         if (norm_Z > BPACK_SafeUnderflow) then
            norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
         else
            norm_UVavrbynorm_Z(itr) = 0
         endif
         itr = mod(itr, Navr) + 1

         stats%Flop_Fill = stats%Flop_Fill + rank*blocks%N_loc + rank*blocks%M_loc

         ! if(ptree%MyID==Main_ID)write(*,*)rank,norm_U*norm_V,norm_Z,select_column(rank+1),select_row(rank+1)

         rank = rank + 1
         if (rank < rankmax_min) then
            select_row(rank + 1) = maxloc(norm_column_R, 1)
         endif

         if (norm_Z < 0) exit

      enddo

      error = sqrt((sum(norm_UVavrbynorm_Z)/Navr))

      deallocate (norm_UVavrbynorm_Z)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

! ! ACA followed by SVD

      !!!!>**** generate 2D grid blacs quantities for matU
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatU2D, M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matU2D(max(1,myArows), max(1,myAcols)))
         matU2D = 0

      else
         descsMatU2D(2) = -1
         allocate (matU2D(1, 1))
         matU2D = 0
      endif
      !!!!>**** redistribution of input matrix
      call Redistribute1Dto2D(matU, blocks%M_p, 0, pgno, matU2D, M, 0, pgno, rank, ptree)

      !!!!>**** generate 2D grid blacs quantities for matV transpose
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatV2D, N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (matV2D(max(1,myArows), max(1,myAcols)))
         matV2D = 0
      else
         descsMatV2D(2) = -1
         allocate (matV2D(1, 1))
         matV2D = 0
      endif
      !!!!>**** redistribution of input matrix
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
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (RR1(max(1,myArows), max(1,myAcols)))
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
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         mn = min(N, rank)
         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn, nbslpk, mycol, 0, npcol)
         allocate (tau_Q(myAcols))
         call pgeqrff90(N, mn, matV2D, 1, 1, descsMatV2D, tau_Q, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (RR2(max(1,myArows), max(1,myAcols)))
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
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         allocate (mattemp(max(1,myArows), max(1,myAcols)))
         mattemp = 0
         call descinit_wp(descsMatSml, rank, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', rank, rank, rank, BPACK_cone, RR1, 1, 1, descsMatSml, RR2, 1, 1, descsMatSml, BPACK_czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         allocate (UUsml(max(1,myArows), max(1,myAcols)), VVsml(max(1,myArows), max(1,myAcols)), Singularsml(rank))
         call PSVD_Truncate(rank, rank, mattemp, descsMatSml, UUsml, VVsml, descsMatSml, descsMatSml, Singularsml, SVD_tolerance, ranknew, ctxt, BPACK_SafeUnderflow, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyU%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
         blocks%ButterflyU%blocks(1)%matrix = 0
         call descinit_wp(descsMatU2Dnew, M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', M, ranknew, rank, BPACK_cone, matU2D, 1, 1, descsMatU2D, UUsml, 1, 1, descsMatSml, BPACK_czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            blocks%ButterflyU%blocks(1)%matrix(:, myj) = blocks%ButterflyU%blocks(1)%matrix(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyV%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
         blocks%ButterflyV%blocks(1)%matrix = 0
         call descinit_wp(descsMatV2Dnew, N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', N, ranknew, rank, BPACK_cone, matV2D, 1, 1, descsMatV2D, VVsml, 1, 1, descsMatSml, BPACK_czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

         rank = ranknew

         deallocate (mattemp, RR1, UUsml, VVsml, Singularsml)
         deallocate (RR2)

#else
         mn1 = min(M, rank)
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn1, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUU_u, M, mn1, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUu(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVV_u, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVu(max(1,myArows), max(1,myAcols)))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(M, rank, matU2D, descsMatU2D, UUu, VVu, descsUU_u, descsVV_u, Singularuv, SVD_tolerance, rank1, ctxt, BPACK_SafeUnderflow, flop=flop)
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
         call descinit_wp(descsUU_v, N, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUv(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVV_v, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVv(max(1,myArows), max(1,myAcols)))
         allocate (Singularuv(mn1))
         call PSVD_Truncate(N, rank, matV2D, descsMatV2D, UUv, VVv, descsUU_v, descsVV_v, Singularuv, SVD_tolerance, rank2, ctxt, BPACK_SafeUnderflow, flop=flop)
         do ii = 1, rank2
            call g2l(ii, rank2, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VVv(myi, :) = VVv(myi, :)*Singularuv(ii)
            endif
         enddo
         deallocate (Singularuv)

         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         allocate (mattemp(max(1,myArows), max(1,myAcols)))
         mattemp = 0
         call descinit_wp(descsMatSml, rank1, rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', rank1, rank2, rank, BPACK_cone, VVu, 1, 1, descsVV_u, VVv, 1, 1, descsVV_v, BPACK_czero, mattemp, 1, 1, descsMatSml, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(min(rank1, rank2), nbslpk, mycol, 0, npcol)
         call descinit_wp(descsUUSml, rank1, min(rank1, rank2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UUsml(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(min(rank1, rank2), nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsVVSml, min(rank1, rank2), rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VVsml(max(1,myArows), max(1,myAcols)))
         allocate (Singularsml(min(rank1, rank2)))
         call PSVD_Truncate(rank1, rank2, mattemp, descsMatSml, UUsml, VVsml, descsUUSml, descsVVSml, Singularsml, SVD_tolerance, ranknew, ctxt, BPACK_SafeUnderflow, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyU%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
         blocks%ButterflyU%blocks(1)%matrix = 0
         call descinit_wp(descsMatU2Dnew, M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'N', M, ranknew, rank1, BPACK_cone, UUu, 1, 1, descsUU_u, UUsml, 1, 1, descsUUSml, BPACK_czero, blocks%ButterflyU%blocks(1)%matrix, 1, 1, descsMatU2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
         do myj = 1, myAcols
            call l2g(myj, mycol, ranknew, npcol, nbslpk, jj)
            blocks%ButterflyU%blocks(1)%matrix(:, myj) = blocks%ButterflyU%blocks(1)%matrix(:, myj)*Singularsml(jj)
         enddo

         myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
         allocate (blocks%ButterflyV%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
         blocks%ButterflyV%blocks(1)%matrix = 0
         call descinit_wp(descsMatV2Dnew, N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         call pgemmf90('N', 'T', N, ranknew, rank2, BPACK_cone, UUv, 1, 1, descsUU_v, VVsml, 1, 1, descsVVSml, BPACK_czero, blocks%ButterflyV%blocks(1)%matrix, 1, 1, descsMatV2Dnew, flop=flop)
         stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)
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
      deallocate (tmpval)

      return

   end subroutine LR_ACA_Parallel

   subroutine LR_BACA(SVD_Q, header_m, header_n, M, N, rank, tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option, error)


      implicit none

      integer rank, rankup, ranknew, row, column, rankmax, N, M, rmax, rmax0
      DT, allocatable:: row_R(:, :), row_Rtmp(:, :), column_R(:, :), column_RT(:, :), fullmat(:, :), fullmat1(:, :)
      DT,allocatable::matU(:, :), matV(:, :)
      DTR,allocatable::Singular(:)
      DT, allocatable :: core(:, :), core_inv(:, :), tau(:), matUtmp(:, :), matVtmp(:, :)
      real(kind=8):: normA, normUV, flop
      integer itr, itrmax, r_est, Nqr, bsize
      integer, allocatable:: select_column(:), select_row(:), perms(:)
      integer, allocatable :: jpvt(:)

      integer i, j, ii, jj, indx, rank_1, rank_2, rankmax_min
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
      type(intersect)::submats(1)
      type(SVD_quant)::SVD_Q

      n1 = MPI_Wtime()

      rmax = size(SVD_Q%MatU,2)
      rankmax_min = min(M, N)
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

         !>**** create random column index for the first iteration
         if (rank == 0) then
            call rperm(N, perms)
            select_column = perms(1:r_est)
         endif

         !>**** Compute columns column_R to find a new set of rows and columns
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo

         submats(1)%nr = M
         submats(1)%nc = r_est
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         column_R = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -BPACK_cone, SVD_Q%matU, M, SVD_Q%matV(1, select_column(j)), rmax, BPACK_cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
            enddo
         endif

         !>**** Find row pivots from the columns column_R
         call copymatT(column_R, column_RT, M, r_est)
         jpvt = 0
         ! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
         call geqp3f90(column_RT, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_row(1:r_est) = jpvt(1:r_est)

         !>**** Compute rows row_R in CUR
         do i = 1, r_est
            mrange(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, N
            nrange(j) = header_n + j - 1
         enddo

         submats(1)%nr = r_est
         submats(1)%nc = N
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         row_R = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)


         ! call element_Zmn_block_user(r_est, N, mrange, nrange, row_R, msh, option, ker, 0, passflag, ptree, stats)

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
               call gemmf77('N', 'N', 1, N, rank, -BPACK_cone, SVD_Q%matU(select_row(i), 1), M, SVD_Q%matV, rmax, BPACK_cone, row_R(i, 1), r_est)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N, rank)
            enddo
         endif

         !>**** Find column pivots from the rows row_R
         jpvt = 0
         row_Rtmp = row_R
         ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
         call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_column(1:r_est) = jpvt(1:r_est)

         !>**** Compute columns column_R in CUR
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo

         submats(1)%nr = M
         submats(1)%nc = r_est
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         column_R = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)


         ! call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -BPACK_cone, SVD_Q%matU, M, SVD_Q%matV(1, select_column(j)), rmax, BPACK_cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
            enddo
         endif

         !>**** Compute the skeleton matrix in CUR
         do i = 1, r_est
            core(i, :) = column_R(select_row(i), :)
         enddo

#if 1
         !>**** generate column indices for the next iteration
         row_Rtmp = row_R
         do j = 1, r_est
            row_Rtmp(:, select_column(j)) = 0d0
         enddo

         jpvt = 0
         ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
         call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_column(1:r_est) = jpvt(1:r_est)
#endif

         !>**** form the LR update by CUR

         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rankup = ranknew

         if (rankup > 0) then
            row_Rtmp = row_R
            call un_or_mqrf90(core, tau, row_Rtmp, 'L', 'C', r_est, N, rankup, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call trsmf90(core, row_Rtmp, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            if(rank + rankup > rankmax_min)then
               rankup=rankmax_min-rank
               if(rankup==0)exit
            endif
            if (rank + rankup > rmax)then
               allocate(matU(M,rmax))
               allocate(matV(rmax,N))
               allocate(Singular(rmax))
               matU=SVD_Q%matU
               matV=SVD_Q%matV
               Singular=SVD_Q%Singular
               rmax0=rmax
               rmax = min(max(rank + rankup,2*rmax),rankmax_min)
               deallocate(SVD_Q%matU,SVD_Q%matV,SVD_Q%Singular)
               allocate(SVD_Q%matU(M,rmax))
               allocate(SVD_Q%matV(rmax,N))
               allocate(SVD_Q%Singular(rmax))
               SVD_Q%matU=0
               SVD_Q%matV=0
               SVD_Q%Singular=0
               SVD_Q%matU(:,1:rmax0) = matU(:,1:rmax0)
               SVD_Q%matV(1:rmax0,:) = matV(1:rmax0,:)
               SVD_Q%Singular(1:rmax0) = Singular(1:rmax0)
               deallocate(matU)
               deallocate(matV)
               deallocate(Singular)
            endif

            ! call assert(rank+rankup<=rmax,'try to increase rmax')
            do j = 1, rankup
               SVD_Q%matU(:, rank + j) = column_R(:, jpvt(j))
            enddo
            SVD_Q%matV(rank + 1:rank + rankup, :) = row_Rtmp(1:rankup, :)
            rank = rank + rankup

            !>**** update fnorm of UV and matUmatV
            call LR_Fnorm(column_R, row_Rtmp, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            ! if(rankup<8)then ! update fnorm seems more efficienct than recompute fnorm when block size is small
            call LR_FnormUp(SVD_Q%matU, SVD_Q%matV, M, N, 0, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)
            ! else
            ! call LR_Fnorm(SVD_Q%matU,SVD_Q%matV,M,N,rank,normA,tolerance*1e-2,Flops=flop)
            ! endif
            stats%Flop_Fill = stats%Flop_Fill + flop

            if (normA > BPACK_SafeUnderflow) then
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
         call LR_ReCompression(SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, M, N, rank, ranknew, SVD_tolerance, Flops=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rank = ranknew
      else
         rank = 1
         SVD_Q%matU(:, 1) = 0
         SVD_Q%matV(1, :) = 0
         SVD_Q%Singular(1) = 0
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

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

      return

   end subroutine LR_BACA

   subroutine LR_BACA_noOverlap(SVD_Q, header_m, header_n, M, N, rank, tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option, error)


      implicit none

      integer rank, rank0, rankup, ranknew, row, column, rankmax, N, M, rmax, rmax0
      DT, allocatable:: row_R(:, :), row_R_knn(:, :), row_Rtmp(:, :), row_Rtmp_knn(:, :), column_R(:, :), column_R_knn(:, :), column_RT(:, :), fullmat(:, :), fullmat1(:, :)
      DT,allocatable::matU(:, :), matV(:, :)
      DTR,allocatable::Singular(:)
      DT, allocatable :: core(:, :), core_knn(:, :), core_inv(:, :), tau(:), matUtmp(:, :), matVtmp(:, :)
      real(kind=8):: normA, normUV, flop, maxvalue
      integer itr, itrmax, r_est, r_est_knn_r, r_est_knn, r_est_knn_c, r_est_tmp, Nqr, bsize
      integer, allocatable:: select_column(:), select_column_knn(:), select_row_knn(:), select_column1(:), select_row(:), perms(:), rows(:), columns(:)
      integer, allocatable :: jpvt(:)

      integer i, j, ii, jj, iii, jjj, indx, rank_1, rank_2, rankmax_min
      real(kind=8) tolerance, SVD_tolerance
      integer edge_m, edge_n, header_m, header_n, mn
      real(kind=8) inner_UV, n1, n2, a, error
      type(intersect)::submats(1)

      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer::mrange(M), nrange(N)
      integer::passflag = 0
      type(SVD_quant)::SVD_Q

      n1 = MPI_Wtime()

      rmax = size(SVD_Q%MatU,2)
      rankmax_min = min(M, N)
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
      rank0 = 0
      error=1.0

      !>**** if nearest neighbour is available, select them first
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

            submats(1)%nr = M
            submats(1)%nc = r_est_knn_c
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            column_R_knn = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)


            ! call element_Zmn_block_user(M, r_est_knn_c, mrange, nrange, column_R_knn, msh, option, ker, 0, passflag, ptree, stats)

            do i = 1, r_est_knn_r
               mrange(i) = header_m + select_row_knn(i) - 1
            enddo
            do j = 1, N
               nrange(j) = header_n + j - 1
            enddo

            submats(1)%nr = r_est_knn_r
            submats(1)%nc = N
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            row_R_knn = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)


            ! call element_Zmn_block_user(r_est_knn_r, N, mrange, nrange, row_R_knn, msh, option, ker, 0, passflag, ptree, stats)

            !>**** Compute the skeleton matrix in CUR
            do i = 1, r_est_knn_r
               core_knn(i, :) = column_R_knn(select_row_knn(i), :)
            enddo

            jpvt = 0
            call geqp3modf90(core_knn, jpvt, tau, tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            rankup = ranknew
            if (rankup > 0) then
               if(rank + rankup > rankmax_min)then
                  rankup=rankmax_min-rank
                  call assert(rankup>0,"rankup should not be zero here")
               endif
               if (rank + rankup > rmax)then
                  allocate(matU(M,rmax))
                  allocate(matV(rmax,N))
                  allocate(Singular(rmax))
                  matU=SVD_Q%matU
                  matV=SVD_Q%matV
                  Singular=SVD_Q%Singular
                  rmax0=rmax
                  rmax = min(max(rank + rankup,2*rmax),rankmax_min)
                  deallocate(SVD_Q%matU,SVD_Q%matV,SVD_Q%Singular)
                  allocate(SVD_Q%matU(M,rmax))
                  allocate(SVD_Q%matV(rmax,N))
                  allocate(SVD_Q%Singular(rmax))
                  SVD_Q%matU=0
                  SVD_Q%matV=0
                  SVD_Q%Singular=0
                  SVD_Q%matU(:,1:rmax0) = matU(:,1:rmax0)
                  SVD_Q%matV(1:rmax0,:) = matV(1:rmax0,:)
                  SVD_Q%Singular(1:rmax0) = Singular(1:rmax0)
                  deallocate(matU)
                  deallocate(matV)
                  deallocate(Singular)
               endif

               row_Rtmp_knn = row_R_knn
               call un_or_mqrf90(core_knn, tau, row_Rtmp_knn, 'L', 'C', r_est_knn_r, N, rankup, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop
               call trsmf90(core_knn, row_Rtmp_knn, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop

               columns(rank + 1:rankup + rank) = select_column_knn(jpvt(1:rankup))
               rows(rank + 1:rankup + rank) = select_row_knn(1:rankup)

               ! call assert(rank+rankup<=rmax,'try to increase rmax')
               do j = 1, rankup
                  SVD_Q%matU(:, rank + j) = column_R_knn(:, jpvt(j))
               enddo
               SVD_Q%matV(rank + 1:rank + rankup, :) = row_Rtmp_knn(1:rankup, :)

               rank = rank + rankup
               rank0 = rank

               ! if (rank == rmax)then
               !    goto 10 !>*** skip ACA iteration
               ! endif

               ! !>**** update fnorm of UV and matUmatV (this is commented out to leave normUV=normA=0 for aca iteration)
               ! call LR_Fnorm(column_R_knn, row_Rtmp_knn, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
               ! stats%Flop_Fill = stats%Flop_Fill + flop
               ! call LR_FnormUp(SVD_Q%matU, SVD_Q%matV, M, N, 0, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)

               ! stats%Flop_Fill = stats%Flop_Fill + flop

               ! if (normA > BPACK_SafeUnderflow) then
               !    error = normUV/normA
               ! else
               !    error = 0
               ! endif

               !>**** Find column pivots for the next iteration
               jpvt = 0
               row_Rtmp_knn = row_R_knn
               if (rank > 0) row_Rtmp_knn(:, columns(1:rank)) = 0
               ! call geqp3modf90(row_Rtmp_knn,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
               call geqp3f90(row_Rtmp_knn, jpvt, tau, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop
               select_column(1:r_est) = jpvt(1:r_est)
            ! else
            !    error = 0
               ! goto 20   !>*** no effective rank found using KNN, go to ACA iteration
            endif
         endif
      endif

      do while (normUV >= tolerance*normA .and. itr < itrmax)

         !>**** create random column index for the first iteration
         if (rank == 0) then
            call rperm(N, perms)
            select_column = perms(1:r_est)
            ! do i=1,r_est
            !    select_column(i)=i
            ! enddo
         endif

         select_column1 = select_column

         !>**** Compute columns column_R to find a new set of rows and columns
         do i = 1, M
            mrange(i) = header_m + i - 1
         enddo
         do j = 1, r_est
            nrange(j) = header_n + select_column(j) - 1
         enddo

         submats(1)%nr = M
         submats(1)%nc = r_est
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         column_R = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(M, r_est, mrange, nrange, column_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -BPACK_cone, SVD_Q%matU, M, SVD_Q%matV(1, select_column(j)), rmax, BPACK_cone, column_R(1, j), M)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(M, 1, rank)
               ! call gemmf90(SVD_Q%matU, M,SVD_Q%matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-BPACK_cone,BPACK_cone)
            enddo
         endif

         !>**** Find row pivots from the columns column_R
         call copymatT(column_R, column_RT, M, r_est)
         if (rank > 0) column_RT(:, rows(1:rank)) = 0
         jpvt = 0
         ! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
         call geqp3f90(column_RT, jpvt, tau, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         select_row(1:r_est) = jpvt(1:r_est)

         ! write(*,*)itr, 'BACA ref col', select_column(1:r_est), 'row', select_row(1:r_est), normUV,normA, header_m, header_n, fnorm(column_R,M,r_est)


         !>**** Compute rows row_R in CUR
         ! !$omp end parallel do
         do i = 1, r_est
            mrange(i) = header_m + select_row(i) - 1
         enddo
         do j = 1, N
            nrange(j) = header_n + j - 1
         enddo

         submats(1)%nr = r_est
         submats(1)%nc = N
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         row_R = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(r_est, N, mrange, nrange, row_R, msh, option, ker, 0, passflag, ptree, stats)

         if (rank > 0) then
            do i = 1, r_est
               call gemmf77('N', 'N', 1, N, rank, -BPACK_cone, SVD_Q%matU(select_row(i), 1), M, SVD_Q%matV, rmax, BPACK_cone, row_R(i, 1), r_est)
               stats%Flop_Fill = stats%Flop_Fill + flops_gemm(1, N, rank)
               ! call gemmf90(SVD_Q%matU(select_row(i),1), M,SVD_Q%matV,rmax,row_R(i,1),r_est,'N','N',1,N,rank,-BPACK_cone,BPACK_cone)
            enddo
         endif

         !>**** Compute the skeleton matrix in CUR
         do i = 1, r_est
            core(i, :) = column_R(select_row(i), :)
         enddo
         maxvalue = abs(core(1, 1))

         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rankup = ranknew

         if (rankup > 0) then

            row_Rtmp = row_R
            call un_or_mqrf90(core, tau, row_Rtmp, 'L', 'C', r_est, N, rankup, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call trsmf90(core, row_Rtmp, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            if(rank + rankup > rankmax_min)then
               rankup=rankmax_min-rank
               if(rankup==0)exit
            endif
            if (rank + rankup > rmax)then
               allocate(matU(M,rmax))
               allocate(matV(rmax,N))
               allocate(Singular(rmax))
               matU=SVD_Q%matU
               matV=SVD_Q%matV
               Singular=SVD_Q%Singular
               rmax0=rmax
               rmax = min(max(rank + rankup,2*rmax),rankmax_min)
               deallocate(SVD_Q%matU,SVD_Q%matV,SVD_Q%Singular)
               allocate(SVD_Q%matU(M,rmax))
               allocate(SVD_Q%matV(rmax,N))
               allocate(SVD_Q%Singular(rmax))
               SVD_Q%matU=0
               SVD_Q%matV=0
               SVD_Q%Singular=0
               SVD_Q%matU(:,1:rmax0) = matU(:,1:rmax0)
               SVD_Q%matV(1:rmax0,:) = matV(1:rmax0,:)
               SVD_Q%Singular(1:rmax0) = Singular(1:rmax0)
               deallocate(matU)
               deallocate(matV)
               deallocate(Singular)
            endif

            columns(rank + 1:rankup + rank) = select_column1(jpvt(1:rankup))
            rows(rank + 1:rankup + rank) = select_row(1:rankup)

            ! call assert(rank+rankup<=rmax,'try to increase rmax')
            do j = 1, rankup
               SVD_Q%matU(:, rank + j) = column_R(:, jpvt(j))
            enddo
            SVD_Q%matV(rank + 1:rank + rankup, :) = row_Rtmp(1:rankup, :)

            rank = rank + rankup

            !>**** update fnorm of UV and matUmatV
            call LR_Fnorm(column_R, row_Rtmp, M, N, rankup, normUV, tolerance*1e-2, Flops=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop
            call LR_FnormUp(SVD_Q%matU, SVD_Q%matV, M, N, rank0, rank - rankup, rankup, rmax, normA, normUV, tolerance*1e-2, Flops=flop)

            stats%Flop_Fill = stats%Flop_Fill + flop

            if (normA > BPACK_SafeUnderflow) then
               error = normUV/normA
            else
               error = 0
               exit
            endif

            !>**** Find column pivots for the next iteration
            jpvt = 0
            row_Rtmp = row_R
            if (rank > 0) row_Rtmp(:, columns(1:rank)) = 0
            ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
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

      if (allocated(row_R_knn)) deallocate (row_R_knn)
      if (allocated(row_Rtmp_knn)) deallocate (row_Rtmp_knn)
      if (allocated(column_R_knn)) deallocate (column_R_knn)
      if (allocated(core_knn)) deallocate (core_knn)
      if (allocated(select_column_knn)) deallocate (select_column_knn)
      if (allocated(select_row_knn)) deallocate (select_row_knn)

      if (rank > 0) then
         call LR_ReCompression(SVD_Q%matU, SVD_Q%matV, SVD_Q%Singular, M, N, rank, ranknew, SVD_tolerance, Flops=flop)
         stats%Flop_Fill = stats%Flop_Fill + flop
         rank = ranknew
      else
         rank = 1
         SVD_Q%matU(:, 1) = 0
         SVD_Q%matV(1, :) = 0
         SVD_Q%Singular(1) = 0
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

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

      return

   end subroutine LR_BACA_noOverlap




   subroutine LR_BACA_noOverlap_Oneiteration(header_m, header_n, M,N,acaquants,submatc,submatr, tolerance, SVD_tolerance, bsize, msh, ker, stats, ptree, option,rcflag,flops)


      implicit none
      type(acaquant)::acaquants
      integer rank, rank0, rankup, ranknew, row, column, rankmax, N, M
      DT, allocatable:: row_R(:, :), row_R_knn(:, :), row_Rtmp(:, :), row_Rtmp_knn(:, :), column_R(:, :), column_R_knn(:, :), column_RT(:, :), fullmat(:, :), fullmat1(:, :)
      DT, allocatable:: matU(:,:), matV(:,:)
      DTR, allocatable::Singular(:)
      DT, allocatable :: core(:, :), core_knn(:, :), core_inv(:, :), tau(:), matUtmp(:, :), matVtmp(:, :)
      real(kind=8):: normA, normUV, flop, maxvalue,flops
      integer itr, itrmax, r_est, r_est_knn_r, r_est_knn, r_est_knn_c, r_est_tmp, Nqr, bsize
      integer, allocatable:: select_column(:), select_column_knn(:), select_row_knn(:), select_column1(:), select_row(:), perms(:)
      integer, allocatable :: jpvt(:)

      integer i, j, ii, jj, iii, jjj, indx, rank_1, rank_2
      real(kind=8) tolerance, SVD_tolerance
      integer edge_m, edge_n, header_m, header_n, mn
      real(kind=8) inner_UV, n1, n2, a, error
      type(intersect)::submatc,submatr

      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      integer rcflag ! rows or columns passed in

      n1 = MPI_Wtime()


      flops=0
      r_est = min(bsize, min(M, N))
      ! r_est=min(M,N)

      Nqr = max(max(M, N), r_est)
      allocate (jpvt(Nqr))
      allocate (tau(Nqr))
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

      normA = acaquants%normA
      normUV = acaquants%normUV
      itr = acaquants%itr
      rank = acaquants%rank
      rank0 = acaquants%rank0
      select_column=acaquants%select_column
      select_row=acaquants%select_row
      select_column1 = select_column

      column_R = submatc%dat


      if(rcflag==0)then


         if (rank > 0) then
            do j = 1, r_est
               call gemmf77('N', 'N', M, 1, rank, -BPACK_cone, acaquants%matU, M, acaquants%matV(1, select_column(j)), rank, BPACK_cone, column_R(1, j), M)
               flops = flops + flops_gemm(M, 1, rank)
               ! call gemmf90(matU, M,matV(1,select_column(j)),rmax,column_R(1,j),M,'N','N',M,1,rank,-BPACK_cone,BPACK_cone)
            enddo
            submatc%dat = column_R
         endif

         !>**** Find row pivots from the columns column_R
         call copymatT(column_R, column_RT, M, r_est)
         if (rank > 0) column_RT(:, acaquants%rows(1:rank)) = 0
         jpvt = 0
         ! call geqp3modf90(column_RT,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
         call geqp3f90(column_RT, jpvt, tau, flop=flop)
         flops = flops + flop
         acaquants%select_row(1:r_est) = jpvt(1:r_est)

         ! write(*,*)itr, 'BACA col', acaquants%select_column(1:r_est), 'row', acaquants%select_row(1:r_est), normUV,normA, header_m, header_n, fnorm(column_R,M,r_est)


      else
         row_R = submatr%dat
         if (rank > 0) then
            do i = 1, r_est
               call gemmf77('N', 'N', 1, N, rank, -BPACK_cone, acaquants%matU(select_row(i), 1), M, acaquants%matV, rank, BPACK_cone, row_R(i, 1), r_est)
               flops = flops + flops_gemm(1, N, rank)
               ! call gemmf90(matU(select_row(i),1), M,matV,rmax,row_R(i,1),r_est,'N','N',1,N,rank,-BPACK_cone,BPACK_cone)
            enddo
         endif

         !>**** Compute the skeleton matrix in CUR
         ! write(*,*)select_column(1:r_est),'jiba',header_m, header_n
         do j = 1, r_est
            core(:, j) = row_R(:, select_column(j))
         enddo
         maxvalue = abs(core(1, 1))


         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
         flops = flops + flop
         rankup = ranknew
         if(rank+rankup>min(M,N))rankup = min(M,N)-rank
         if (rankup > 0) then
            row_Rtmp = row_R
            call un_or_mqrf90(core, tau, row_Rtmp, 'L', 'C', r_est, N, rankup, flop=flop)
            flops = flops + flop
            call trsmf90(core, row_Rtmp, 'L', 'U', 'N', 'N', rankup, N, flop=flop)
            flops = flops + flop

            acaquants%columns(rank + 1:rankup + rank) = select_column1(jpvt(1:rankup))
            acaquants%rows(rank + 1:rankup + rank) = select_row(1:rankup)

            ! call assert(rank+rankup<=rmax,'try to increase rmax')
            allocate(matU(M,rank+rankup))
            allocate(matV(rank+rankup,N))
            if(rank>0)then
               matU(:,1:rank) = acaquants%matU(:,1:rank)
               deallocate(acaquants%matU)
            endif
            do j = 1, rankup
               matU(:, rank + j) = column_R(:, jpvt(j))
            enddo
            if(rank>0)then
               matV(1:rank,:) = acaquants%matV(1:rank,:)
               deallocate(acaquants%matV)
            endif
            matV(rank + 1:rank + rankup, :) = row_Rtmp(1:rankup, :)
            allocate(acaquants%matU(M,rank+rankup))
            acaquants%matU=matU
            deallocate(matU)
            allocate(acaquants%matV(rank+rankup,N))
            acaquants%matV=matV
            deallocate(matV)

            rank = rank + rankup
            acaquants%rank=rank


            !>**** update fnorm of UV and matUmatV
            call LR_Fnorm(column_R, row_Rtmp, M, N, rankup, acaquants%normUV, tolerance*1e-2, Flops=flop)
            flops = flops + flop
            call LR_FnormUp(acaquants%matU, acaquants%matV, M, N, acaquants%rank0, rank - rankup, rankup, rank, acaquants%normA, acaquants%normUV, tolerance*1e-2, Flops=flop)

            ! write(*,*)acaquants%normUV,acaquants%normA,'gana',rank,rankup

            flops = flops + flop

            if (acaquants%normA < BPACK_SafeUnderflow) then
               acaquants%finish = .true.
            else
               !>**** Find column pivots for the next iteration
               jpvt = 0
               row_Rtmp = row_R
               if (rank > 0) row_Rtmp(:, acaquants%columns(1:rank)) = 0
               ! call geqp3modf90(row_Rtmp,jpvt,tau,tolerance*1e-2,BPACK_SafeUnderflow,ranknew)
               call geqp3f90(row_Rtmp, jpvt, tau, flop=flop)
               flops = flops + flop
               acaquants%select_column(1:r_est) = jpvt(1:r_est)
            endif
         else
            acaquants%finish = .true.
         endif
         acaquants%itr = acaquants%itr + 1
         if(acaquants%finish .eqv. .false.)acaquants%finish = .not. (acaquants%normUV >= tolerance*acaquants%normA .and. acaquants%itr < acaquants%itrmax)
         if(acaquants%finish .eqv. .true.)then
            if (rank > 0) then
               if(allocated(acaquants%Singular))deallocate(acaquants%Singular)
               allocate(acaquants%Singular(rank))
               call LR_ReCompression(acaquants%matU, acaquants%matV, acaquants%Singular, M, N, rank, ranknew, SVD_tolerance, Flops=flop)
               flops = flops + flop
               rank = ranknew
               acaquants%rank = rank
            else
               rank = 1
               acaquants%rank = rank
               allocate(acaquants%matU(M,rank))
               acaquants%matU=0
               allocate(acaquants%matV(rank,N))
               acaquants%matV=0
               allocate(acaquants%Singular(rank))
               acaquants%Singular =1
            endif
         endif
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
      deallocate (perms)


      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

      return

   end subroutine LR_BACA_noOverlap_Oneiteration






   subroutine LR_SeudoSkeleton(blocks, header_m, header_n, M, N, rmaxc, rmaxr, rank, tolerance, SVD_tolerance, msh, ker, stats, ptree, option, ctxt, pgno)


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
      DTR, allocatable::Singular(:)
      DT, allocatable:: row_R(:), column_R(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:)
      type(mesh)::msh
      type(kernelquant)::ker
      type(matrixblock)::blocks
      type(Hstat)::stats

      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), matU2D(:, :), matV2D(:, :), matU1D(:, :), matV1D(:, :)
      DTR, allocatable :: Singularsml(:)
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
      type(intersect)::submats(1)

      rank_new = 0
      ! ctxt = ptree%pgrp(pgno)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)

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

         submats(1)%nr = rmaxr
         submats(1)%nc = N
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         matV_tmp = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         allocate (MatrixSubselection(rmaxr, rmaxc))
         MatrixSubselection = 0
         do jj = 1, rmaxc
            MatrixSubselection(:,jj)=matV_tmp(:,select_col(jj))
         enddo


         ! call element_Zmn_block_user(rmaxr, N, mrange, nrange, matV_tmp, msh, option, ker, 0, passflag, ptree, stats)
         call copymatT(matV_tmp, matV, rmaxr, N)
         deallocate (matV_tmp)


         allocate (jpiv(rmaxc))
         jpiv = 0
         allocate (tau(min(rmaxr, rmaxc)))
         tau = 0
         call geqp3modf90(MatrixSubselection, jpiv, tau, tolerance, BPACK_SafeUnderflow, rank_new, flop=flop)
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

         submats(1)%nr = M
         submats(1)%nc = rank
         allocate(submats(1)%rows(submats(1)%nr))
         submats(1)%rows = mrange(1:submats(1)%nr)
         allocate(submats(1)%cols(submats(1)%nc))
         submats(1)%cols = nrange(1:submats(1)%nc)
         allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
         call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
         blocks%ButterflyU%blocks(1)%matrix = submats(1)%dat
         deallocate(submats(1)%rows)
         deallocate(submats(1)%cols)
         deallocate(submats(1)%dat)

         ! call element_Zmn_block_user(M, rank, mrange, nrange, blocks%ButterflyU%blocks(1)%matrix, msh, option, ker, 0, passflag, ptree, stats)


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

            call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
            ! nproc = npcol*nprow
            myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rmaxr, nbslpk, mycol, 0, npcol)
            allocate (matV(max(1,myArows), max(1,myAcols)))
            allocate (matV_tmp(max(1,myAcols), max(1,myArows)))
            call descinit_wp(descsmatV, N, rmaxr, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descsmatV')
            matV = 0
            do myj = 1, myAcols
               call l2g(myj, mycol, rmaxr, npcol, nbslpk, jj)
               mrange(myj) = header_m + select_row(jj) - 1
            enddo
            do myi = 1, myArows
               call l2g(myi, myrow, N, nprow, nbslpk, ii)
               nrange(myi) = header_n + ii - 1
            enddo

            submats(1)%nr = myAcols
            submats(1)%nc = myArows
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            matV_tmp = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(myAcols, myArows, mrange, nrange, matV_tmp, msh, option, ker, 0, passflag, ptree, stats)
            call copymatT(matV_tmp, matV, myAcols, myArows)
            deallocate (matV_tmp)

            call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
            myArows = numroc_wp(rmaxr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rmaxc, nbslpk, mycol, 0, npcol)

            allocate (MatrixSubselection(max(1,myArows), max(1,myAcols)))
            call descinit_wp(descsub, rmaxr, rmaxc, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descsub')
            MatrixSubselection = 0
            do myi = 1, myArows
               call l2g(myi, myrow, rmaxr, nprow, nbslpk, ii)
               mrange(myi) = header_m + select_row(ii) - 1
            enddo
            do myj = 1, myAcols
               call l2g(myj, mycol, rmaxc, npcol, nbslpk, jj)
               nrange(myj) = header_n + select_col(jj) - 1
            enddo

            submats(1)%nr = myArows
            submats(1)%nc = myAcols
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            MatrixSubselection = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(myArows, myAcols, mrange, nrange, MatrixSubselection, msh, option, ker, 0, passflag, ptree, stats)

            ! Compute QR of MatrixSubselection*P
            allocate (ipiv(myAcols))
            ipiv = 0
            allocate (tau(myAcols))
            tau = 0
            allocate (jpiv(rmaxc))
            jpiv = 0
            allocate (JPERM(rmaxc))
            JPERM = 0
            call pgeqpfmodf90(rmaxr, rmaxc, MatrixSubselection, 1, 1, descsub, ipiv, tau, JPERM, jpiv, rank_new, tolerance, BPACK_SafeUnderflow, flop=flop)
            stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

            rank = rank_new

            if (rank > 0) then

               ! Compute matV*conjg(Q)
               matV = conjg(cmplx(matV, kind=8))
               call pun_or_mqrf90('R', 'N', N, rmaxr, rank_new, MatrixSubselection, 1, 1, descsub, tau, matV, 1, 1, descsmatV)
               matV = conjg(cmplx(matV, kind=8))

               ! Compute matV*conjg(Q)*(R^T)^-1
               call ptrsmf90('R', 'U', 'T', 'N', N, rank_new, BPACK_cone, MatrixSubselection, 1, 1, descsub, matV, 1, 1, descsmatV, flop=flop)
               stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

               call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
               myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
               allocate (matV2D(max(1,myArows), max(1,myAcols)))
               call descinit_wp(descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
               call assert(info == 0, 'descinit_wp fail for descButterV2D')
               matV2D(1:max(1,myArows), 1:max(1,myAcols)) = matV(1:max(1,myArows), 1:max(1,myAcols))

            else
               rank = 1
               rank_new = 1

               call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
               myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
               myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
               allocate (matV2D(max(1,myArows), max(1,myAcols)))
               call descinit_wp(descButterV2D, N, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
               call assert(info == 0, 'descinit_wp fail for descButterV2D')
               matV2D(1:max(1,myArows), 1:max(1,myAcols)) = 0
            endif

            call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
            myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank_new, nbslpk, mycol, 0, npcol)
            allocate (matU2D(max(1,myArows), max(1,myAcols)))
            call descinit_wp(descButterU2D, M, rank_new, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
            call assert(info == 0, 'descinit_wp fail for descButterU2D')
            matU2D = 0

            do myi = 1, myArows
               call l2g(myi, myrow, M, nprow, nbslpk, ii)
               mrange(myi) = header_m + ii - 1
            enddo
            do myj = 1, myAcols
               call l2g(myj, mycol, rank_new, npcol, nbslpk, jj)
               nrange(myj) = header_n + select_col(ipiv(myj)) - 1
            enddo

            submats(1)%nr = myArows
            submats(1)%nc = myAcols
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            matU2D = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(myArows, myAcols, mrange, nrange, matU2D, msh, option, ker, 0, passflag, ptree, stats)

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
            allocate (blocks%ButterflyU%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
            blocks%ButterflyU%blocks(1)%matrix(1:max(1,myArows), 1:max(1,myAcols)) = matU2D(1:max(1,myArows), 1:max(1,myAcols))
            deallocate (matU2D)

            myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
            allocate (blocks%ButterflyV%blocks(1)%matrix(max(1,myArows), max(1,myAcols)))
            blocks%ButterflyV%blocks(1)%matrix(1:max(1,myArows), 1:max(1,myAcols)) = matV2D(1:max(1,myArows), 1:max(1,myAcols))
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
      type(intersect)::submats(1)

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
      call Bplus_block_MVP_dat(bplus, 'N', blocks%M_loc, blocks%N_loc, Ntest, Vin, blocks%N_loc, Vout2, blocks%M_loc, BPACK_cone, BPACK_czero, ptree, stats)

      allocate (Fullmat(blocks%M_loc, blocks%N))
      allocate (mrange(blocks%M_loc))
      allocate (nrange(blocks%N))
      do myi = 1, blocks%M_loc
         mrange(myi) = blocks%M_p(pp, 1) + blocks%headm - 1 + myi - 1
      enddo
      do myj = 1, blocks%N
         nrange(myj) = blocks%headn + myj - 1
      enddo

      submats(1)%nr = blocks%M_loc
      submats(1)%nc = blocks%N
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      Fullmat = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)

      ! call element_Zmn_block_user(blocks%M_loc, blocks%N, mrange, nrange, Fullmat, msh, option, ker, 0, passflag, ptree, stats)
      deallocate (mrange)
      deallocate (nrange)

      passflag = 0
      do while (passflag == 0)
         ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 1, passflag, ptree, stats)
         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 1, passflag, ptree, stats)
      enddo

      call gemmf90(Fullmat, blocks%M_loc, Vin_glo, blocks%N, Vout1, blocks%M_loc, 'N', 'N', blocks%M_loc, Ntest, blocks%N, BPACK_cone, BPACK_czero)

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


   subroutine ButterflySVD_Left(index_i, index_j, level, level_butterfly, levelm, blocks, option, msh, ButterflyP_old, ButterflyP, ButterflyMiddle, rank,flops)

      implicit none
      integer  offsetm, level_butterfly, index_i, index_j, index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc_k,index_j_loc_k,index_i_loc_s,index_j_loc_s,level, group_m, group_m0, mm, mm0, nn,nn1,nn2, j, i, mn, rank, mm1, levelm
      type(butterfly_kerl) ButterflyP_old, ButterflyP, ButterflyMiddle
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      DTR, allocatable :: Singular(:)
      type(matrixblock)::blocks
      real(kind=8):: SVD_tolerance, flops, flop
      type(Hoption)::option
      type(mesh)::msh

      flops = 0

      index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1
      index_ii_loc = (index_ii - ButterflyP_old%idx_r)/ButterflyP_old%inc_r + 1
      index_jj_loc = (index_jj - ButterflyP_old%idx_c)/ButterflyP_old%inc_c + 1
      index_i_loc_s = (index_i - ButterflyP%idx_r)/ButterflyP%inc_r + 1
      index_j_loc_s = (index_j - ButterflyP%idx_c)/ButterflyP%inc_c + 1

      group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
      group_m = group_m*2**level - 1 + index_i
      group_m0 = int(group_m/2)

      mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
      mm0 = msh%basis_group(group_m0)%tail - msh%basis_group(group_m0)%head + 1

      if (size(ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix, 1) /= mm0) then
         write (*, *) 'mm0 incorrect'
         stop
      end if
      nn1 = size(ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
      nn2 = size(ButterflyP_old%blocks(index_ii_loc, index_jj_loc+1)%matrix, 2)
      nn= nn1+nn2

      allocate (QQ(mm, nn))
      if(mod(index_i,2)==1)offsetm=0
      if(mod(index_i,2)==0)offsetm=mm0-mm

      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn1
         do i = 1, mm
            QQ(i, j) = ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix(i+offsetm, j)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn2
         do i = 1, mm
            QQ(i, j+nn1) = ButterflyP_old%blocks(index_ii_loc, index_jj_loc+1)%matrix(i+offsetm, j)
         enddo
      enddo
      ! !$omp end parallel do


      ! write(*,*)'dddd',fnorm(QQ,mm,nn)

      mn = min(mm, nn)
      allocate (UU(mm, mn), VV(mn, nn), Singular(mn))

#if 1
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank, flop=flop)
      flops = flops + flop

      if(level==level_butterfly)then
         do i=1,rank
            VV(i, :) =  Singular(i)*VV(i, :)
         enddo
      else
         do j=1,rank
            UU(:, j) =  UU(:, j)*Singular(j)
         enddo
      endif
#else
      if(level==level_butterfly)then
         call RRQR_LQ(QQ, mm, nn, mn, UU, VV, option%tol_comp, rank, 'L', flops=flop)
      else
         call RRQR_LQ(QQ, mm, nn, mn, UU, VV, option%tol_comp, rank, 'R', flops=flop)
      endif
      flops = flops + flop
#endif

      if(levelm+1==level)then ! merge ButterflyMiddle into ButterflyKerl
         do j = 1, nn1
               VV(:, j) = VV(:, j)* ButterflyMiddle%blocks(index_ii_loc, index_jj_loc)%matrix(j, j)
         enddo
         do j = 1, nn2
            VV(:, j+nn1) = VV(:, j+nn1)* ButterflyMiddle%blocks(index_ii_loc, index_jj_loc+1)%matrix(j, j)
         enddo
      endif

      allocate (mat_tmp(mm, rank))
      ! !$omp parallel do default(shared) private(i,j,k,ctemp)
      do j = 1, rank
         do i = 1, mm
            mat_tmp(i, j) = UU(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      allocate (ButterflyP%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, rank))
      ButterflyP%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat_tmp(1:mm, 1:rank)

      index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
      index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1
      allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(rank, nn1))
      allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(rank, nn2))

      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, rank
         do j = 1, nn1
            blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(i, j) = VV(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, rank
         do j = 1, nn2
            blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(i, j) = VV(i, j + nn1)
         enddo
      enddo
      ! !$omp end parallel do

      deallocate (QQ, UU, VV, Singular, mat_tmp)

   end subroutine ButterflySVD_Left


   subroutine LocalButterflySVD_Left(index_i_loc, index_j_loc, level_loc, level_butterflyL, level, index_i_m, blocks, option, msh, ButterflyP_old, ButterflyP)

      implicit none
      integer index_i_loc, index_j_loc, level_loc, level_butterflyL, index_i_m, index_i, index_j, level, group_m, mm, nn, nn1, nn2, j, i, mn, rank, mm1
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      DTR, allocatable :: Singular(:)
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
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
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



   subroutine ButterflySVD_Right(index_i, index_j, level, level_butterfly, blocks, option, msh, ButterflyP_old, ButterflyP,rank,flops)

      implicit none
      integer level_butterfly, index_j_m, index_i, index_j, index_ii, index_ii_loc,index_jj,index_jj_loc, index_i_loc_k, index_j_loc_k,index_i_loc_s,index_j_loc_s,level, group_n, group_n0, mm, nn, nn0, offsetn, nn1, mm1, mm2, j, i, mn, rank
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      DTR, allocatable :: Singular(:)
      type(matrixblock)::blocks
      real(kind=8):: SVD_tolerance,flops,flop
      type(Hoption)::option
      type(mesh)::msh

      flops=0

      index_jj = int((index_j + 1)/2); index_ii = 2*index_i - 1
      index_ii_loc = (index_ii - ButterflyP_old%idx_r)/ButterflyP_old%inc_r + 1
      index_jj_loc = (index_jj - ButterflyP_old%idx_c)/ButterflyP_old%inc_c + 1
      index_i_loc_s = (index_i - ButterflyP%idx_r)/ButterflyP%inc_r + 1
      index_j_loc_s = (index_j - ButterflyP%idx_c)/ButterflyP%inc_c + 1

      group_n = blocks%col_group   ! Note: row_group and col_group interchanged here
      group_n = group_n*2**(level_butterfly - level + 1) - 1 + index_j
      group_n0 = int(group_n/2)
      nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
      nn0 = msh%basis_group(group_n0)%tail - msh%basis_group(group_n0)%head + 1

      if (size(ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix, 2) /= nn0) then
         write (*, *) 'nn0 incorrect'
         stop
      end if
      mm1 = size(ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix, 1)
      mm2 = size(ButterflyP_old%blocks(index_ii_loc+1, index_jj_loc)%matrix, 1)
      mm = mm1 + mm2

      !!!!!!!!!!!!!!!!!!

      allocate (QQ(mm, nn))
      if(mod(index_j,2)==1)offsetn=0
      if(mod(index_j,2)==0)offsetn=nn0-nn

      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, mm1
            QQ(i, j) = ButterflyP_old%blocks(index_ii_loc, index_jj_loc)%matrix(i, j+offsetn)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, mm2
            QQ(i + mm1, j) = ButterflyP_old%blocks(index_ii_loc+1, index_jj_loc)%matrix(i, j+offsetn)
         enddo
      enddo
      ! !$omp end parallel do

      mn = min(mm, nn)
      allocate (UU(mm, mn), VV(mn, nn), Singular(mn))

#if 1
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank, flop=flop)
      flops = flops + flop

      if(level==1)then
         do j=1,rank
            UU(:, j) =  UU(:, j)*Singular(j)
         enddo
      else
         do i=1,rank
            VV(i, :) =  Singular(i)*VV(i, :)
         enddo
      endif
#else
      if(level==1)then
         call RRQR_LQ(QQ, mm, nn, mn, UU, VV, option%tol_comp, rank, 'L', flops=flop)
      else
         call RRQR_LQ(QQ, mm, nn, mn, UU, VV, option%tol_comp, rank, 'R', flops=flop)
      endif
      flops = flops + flop
#endif

      allocate (mat_tmp(rank, nn))
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, nn
         do i = 1, rank
            mat_tmp(i, j) = VV(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      allocate (ButterflyP%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank, nn))
      ButterflyP%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat_tmp(1:rank, 1:nn)

      index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1
      index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

      allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(mm1, rank))
      allocate (blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(mm2, rank))

      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, mm1
         do j = 1, rank
            blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(i, j) = UU(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(i,j)
      do i = 1, mm2
         do j = 1, rank
            blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(i, j) = UU(i + mm1, j)
         enddo
      enddo
      ! !$omp end parallel do

      deallocate (QQ, UU, VV, Singular, mat_tmp)

   end subroutine ButterflySVD_Right


   subroutine LocalButterflySVD_Right(index_i_loc, index_j_loc, level_loc, level_butterflyR, level, level_butterfly, index_j_m, blocks, option, msh, ButterflyP_old, ButterflyP)

      implicit none
      integer index_i_loc, index_j_loc, level_loc, level_butterflyR, level_butterfly, index_j_m, index_i, index_j, level, group_n, mm, nn, nn1, mm1, mm2, j, i, mn, rank
      type(butterfly_kerl) ButterflyP_old, ButterflyP
      DT, allocatable :: QQ(:, :), RR(:, :), UU(:, :), VV(:, :), mat_tmp(:, :), matU(:, :), matV(:, :)
      DT, allocatable :: tau_Q(:)
      DTR, allocatable :: Singular(:)
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
      call SVD_Truncate(QQ, mm, nn, mn, UU, VV, Singular, option%tol_comp, BPACK_SafeUnderflow, rank)
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


   subroutine Full_construction(blocks, msh, ker, stats, option, ptree, memory)


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
      type(intersect)::submats(1)
      real(kind=8)::memory

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


      submats(1)%nr = mm
      submats(1)%nc = nn
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      blocks%fullmat = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)

#if HAVE_ZFP
      if(option%use_zfp==1)then
         call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N, option%tol_comp,0)
         memory = SIZEOF(blocks%FullmatZFP%buffer_r)/1024.0d3
#if DAT==0 || DAT==2
         memory = memory + SIZEOF(blocks%FullmatZFP%buffer_i)/1024.0d3
#endif
      else
         memory = SIZEOF(blocks%fullmat)/1024.0d3
      endif
#else
      memory = SIZEOF(blocks%fullmat)/1024.0d3
#endif

      deallocate (mrange)
      deallocate (nrange)

      return

   end subroutine Full_construction

end module Bplus_compress
