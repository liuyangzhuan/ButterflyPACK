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

!> @file Bplus_utilities.f90
!> @brief Block-level utility subroutines not included in other files

#include "ButterflyPACK_config.fi"

module Bplus_Utilities
   use BPACK_DEFS
   use MISC_Utilities
   use magma_utilities
contains

   subroutine Bplus_delete(bplus)


      implicit none
      type(matrixblock), pointer::block
      type(blockplus)::bplus

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall
      real(kind=8)::rtemp

      if (associated(bplus%LL)) then
      do ll = 1, LplusMax
         if (bplus%LL(ll)%Nbound > 0) then
            if (associated(bplus%LL(ll)%matrices_block)) then
            do bb = 1, bplus%LL(ll)%Nbound
               ! write(*,*)ll,bplus%Lplus,bb,bplus%LL(ll)%Nbound,'fff'
               call BF_delete(bplus%LL(ll)%matrices_block(bb), 1)
            end do
            deallocate (bplus%LL(ll)%matrices_block)
            endif
            if (allocated(bplus%LL(ll)%boundary_map)) deallocate (bplus%LL(ll)%boundary_map)
         end if
      end do
      deallocate (bplus%LL)
      endif

   end subroutine Bplus_delete

   subroutine Bplus_copy(bplus_i, bplus_o, memory)


      implicit none
      type(matrixblock), pointer::block_i, block_o
      type(blockplus)::bplus_i, bplus_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall, Ninadmissible
      real(kind=8), optional::memory
      real(kind=8)::rtemp

      call Bplus_delete(bplus_o)

      if (present(memory)) memory = 0

      allocate (bplus_o%LL(LplusMax))
      bplus_o%Lplus = bplus_i%Lplus
      bplus_o%boundary = bplus_i%boundary
      bplus_o%level = bplus_i%level
      bplus_o%col_group = bplus_i%col_group
      bplus_o%row_group = bplus_i%row_group
      bplus_o%pgno = bplus_i%pgno

      do ll = 1, LplusMax
         bplus_o%LL(ll)%Nbound = bplus_i%LL(ll)%Nbound
         bplus_o%LL(ll)%rankmax = bplus_i%LL(ll)%rankmax

         if (bplus_i%LL(ll)%Nbound > 0) then
            allocate (bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
            do bb = 1, bplus_i%LL(ll)%Nbound
               call BF_copy('N', bplus_i%LL(ll)%matrices_block(bb), bplus_o%LL(ll)%matrices_block(bb), rtemp)
               if (present(memory)) memory = memory + rtemp
            end do
            if (allocated(bplus_i%LL(ll)%boundary_map)) then
               Nboundall = size(bplus_i%LL(ll)%boundary_map,1)
               Ninadmissible = size(bplus_i%LL(ll)%boundary_map,2)
               allocate (bplus_o%LL(ll)%boundary_map(Nboundall,Ninadmissible))
               if (present(memory)) memory = memory + SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
               bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
            endif
         end if
      end do

   end subroutine Bplus_copy

   subroutine Bplus_copy_delete(bplus_i, bplus_o, memory)


      implicit none
      type(matrixblock), pointer::block_i, block_o
      type(blockplus)::bplus_i, bplus_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall, Ninadmissible
      real(kind=8), optional::memory
      real(kind=8)::rtemp

      if (present(memory)) memory = 0

      allocate (bplus_o%LL(LplusMax))
      bplus_o%Lplus = bplus_i%Lplus
      bplus_o%boundary = bplus_i%boundary
      bplus_o%level = bplus_i%level
      bplus_o%col_group = bplus_i%col_group
      bplus_o%row_group = bplus_i%row_group

      do ll = 1, LplusMax
         bplus_o%LL(ll)%Nbound = bplus_i%LL(ll)%Nbound
         bplus_o%LL(ll)%rankmax = bplus_i%LL(ll)%rankmax
         if (bplus_i%LL(ll)%Nbound > 0) then
            allocate (bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
            do bb = 1, bplus_i%LL(ll)%Nbound
               call BF_copy_delete(bplus_i%LL(ll)%matrices_block(bb), bplus_o%LL(ll)%matrices_block(bb), rtemp)
               if (present(memory)) memory = memory + rtemp
            end do
            deallocate (bplus_i%LL(ll)%matrices_block)
            Nboundall = size(bplus_i%LL(ll)%boundary_map,1)
            Ninadmissible = size(bplus_i%LL(ll)%boundary_map,2)
            allocate (bplus_o%LL(ll)%boundary_map(Nboundall,Ninadmissible))
            if (present(memory)) memory = memory + SIZEOF(bplus_o%LL(ll)%boundary_map)/1024.0d3
            bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
            deallocate (bplus_i%LL(ll)%boundary_map)
         end if
      end do

      deallocate (bplus_i%LL)
   end subroutine Bplus_copy_delete

   subroutine Bplus_extract_partial(bplus_i, ll_s, row_group, agent_bplus, msh)


      implicit none
      type(matrixblock), pointer::block_i, block_o
      type(blockplus)::bplus_i, agent_bplus

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb, bb_o
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall
      real(kind=8)::rtemp
      integer row_group, ll_s, idx_s, idx_e
      type(mesh)::msh

      call assert(bplus_i%row_group == bplus_i%col_group, 'only works for square matrix')

      idx_s = msh%basis_group(row_group)%head
      idx_e = msh%basis_group(row_group)%tail

! allocate(agent_bplus)
      allocate (agent_bplus%LL(LplusMax))
      do ll = 1, LplusMax
         agent_bplus%LL(ll)%Nbound = 0
      end do

      agent_bplus%Lplus = bplus_i%Lplus - ll_s + 1
      agent_bplus%row_group = row_group
      agent_bplus%col_group = row_group
      agent_bplus%level = GetTreelevel(row_group) - 1

      do ll = 1, agent_bplus%Lplus
         agent_bplus%LL(ll)%Nbound = 0
         agent_bplus%LL(ll)%rankmax = bplus_i%LL(ll + ll_s - 1)%rankmax
         do bb = 1, bplus_i%LL(ll + ll_s - 1)%Nbound
            if (msh%basis_group(bplus_i%LL(ll + ll_s - 1)%matrices_block(bb)%row_group)%head >= idx_s .and. msh%basis_group(bplus_i%LL(ll + ll_s - 1)%matrices_block(bb)%row_group)%tail <= idx_e) then
               agent_bplus%LL(ll)%Nbound = agent_bplus%LL(ll)%Nbound + 1
            end if
         end do
         if (agent_bplus%LL(ll)%Nbound > 0) then
            allocate (agent_bplus%LL(ll)%matrices_block(agent_bplus%LL(ll)%Nbound))
         end if
      end do

      do ll = 1, agent_bplus%Lplus
         bb_o = 0
         do bb = 1, bplus_i%LL(ll + ll_s - 1)%Nbound
            if (msh%basis_group(bplus_i%LL(ll + ll_s - 1)%matrices_block(bb)%row_group)%head >= idx_s .and. msh%basis_group(bplus_i%LL(ll + ll_s - 1)%matrices_block(bb)%row_group)%tail <= idx_e) then
               bb_o = bb_o + 1
               call BF_copy('N', bplus_i%LL(ll + ll_s - 1)%matrices_block(bb), agent_bplus%LL(ll)%matrices_block(bb_o))
            end if
         end do
      end do

   end subroutine Bplus_extract_partial

   subroutine Bplus_ComputeMemory(bplus_i, memory, rank)


      implicit none
      type(matrixblock), pointer::block_i, block_o
      type(blockplus)::bplus_i, bplus_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall
      real(kind=8)::memory
      real(kind=8)::rtemp

      memory = 0
      rank = 0

      do ll = 1, LplusMax
         if (bplus_i%LL(ll)%Nbound > 0) then
            do bb = 1, bplus_i%LL(ll)%Nbound
               call BF_ComputeMemory(bplus_i%LL(ll)%matrices_block(bb), rtemp)
               memory = memory + rtemp
               rank = max(rank,bplus_i%LL(ll)%matrices_block(bb)%rankmax)
            end do
         end if
      end do

   end subroutine Bplus_ComputeMemory

   logical function Bplus_checkNAN(bplus_i)


      implicit none
      type(matrixblock), pointer::block_i, block_o
      type(blockplus)::bplus_i, bplus_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, ll, bb
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, Nboundall
      real(kind=8)::rtemp

      Bplus_checkNAN = .false.

#ifdef NDEBUG
      Bplus_checkNAN = .false.
#else
      do ll = 1, LplusMax
         if (bplus_i%LL(ll)%Nbound > 0) then
            do bb = 1, bplus_i%LL(ll)%Nbound
               if (BF_checkNAN(bplus_i%LL(ll)%matrices_block(bb))) then
                  Bplus_checkNAN = .true.
                  return
               end if
            end do
         end if
      end do
#endif
   end function Bplus_checkNAN




   subroutine Hmat_parallelblock_MVP_dat(blocks_1, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)



      implicit none

      integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b, ctemp1, ctemp2
      character chara
      type(matrixblock), pointer::blocks
      type(matrixblock)::blocks_1
      integer ll, bb, Maxlevel
      type(proctree)::ptree
      type(Hstat)::stats
      integer:: ldi, ldo
      DT :: random1(ldi, *), random2(ldo, *)
      DT, allocatable :: Vout(:, :), Vin_loc(:, :), Vout_loc(:, :)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)
      real(kind=8)::n2,n1
      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)
      type(vectorsblock),allocatable:: Vin_locs(:), Vout_locs(:)

      integer idx_start_m, idx_start_n, idx_start_n_loc, idx_start_m_loc, idx_end_n_loc, idx_end_m_loc, idx_start_i_loc, idx_start_o_loc, idx_end_i_loc, idx_end_o_loc

      type(nod), pointer::cur
      class(*), pointer::ptr

      if (chara == 'N') allocate (Vout(M, Nrnd))
      if (chara == 'T') allocate (Vout(N, Nrnd))
      Vout = 0

      ctemp1 = 1.0d0; ctemp2 = 0.0d0
      Maxlevel = size(blocks_1%lstblks)-1

      do level = 0, Maxlevel

         allocate(Vin_locs(blocks_1%lstblks(level)%num_nods))
         allocate(Vout_locs(blocks_1%lstblks(level)%num_nods))
         ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
         cur => blocks_1%lstblks(level)%head
         do bb = 1, blocks_1%lstblks(level)%num_nods
            select type (ptr=>cur%item)
            type is (block_ptr)
               blocks => ptr%ptr
               n1 = MPI_Wtime()
               if (chara == 'N') then
                  if (blocks%M_loc > 0) then
                     allocate (Vout_locs(bb)%vector(blocks%M_loc, Nrnd))
                     Vout_locs(bb)%vector=0
                  endif
                  if (blocks%N_loc > 0) allocate (Vin_locs(bb)%vector(blocks%N_loc, Nrnd))
                  call Redistribute1Dto1D(random1, ldi, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Vin_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Nrnd, ptree)
               else
                  if (blocks%N_loc > 0)then
                     allocate (Vout_locs(bb)%vector(blocks%N_loc, Nrnd))
                     Vout_locs(bb)%vector=0
                  endif
                  if (blocks%M_loc > 0) allocate (Vin_locs(bb)%vector(blocks%M_loc, Nrnd))
                  call Redistribute1Dto1D(random1, ldi, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Vin_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Nrnd, ptree)
               endif
               n2 = MPI_Wtime()
               stats%Time_RedistV = stats%Time_RedistV + n2-n1

            end select
            cur => cur%next
         enddo


         cur => blocks_1%lstblks(level)%head
         do bb = 1, blocks_1%lstblks(level)%num_nods
            select type (ptr=>cur%item)
            type is (block_ptr)
               blocks => ptr%ptr
               if (blocks%N_loc > 0 .or. blocks%M_loc > 0) then
                  call BF_block_MVP_dat(blocks, chara, blocks%M_loc, blocks%N_loc, Nrnd,&
                  &Vin_locs(bb)%vector, size(Vin_locs(bb)%vector,1), Vout_locs(bb)%vector, size(Vout_locs(bb)%vector,1), ctemp1, ctemp2, ptree, stats)
               endif
            end select
            cur => cur%next
         enddo


         ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
         cur => blocks_1%lstblks(level)%head
         do bb = 1, blocks_1%lstblks(level)%num_nods
            select type (ptr=>cur%item)
            type is (block_ptr)
               blocks => ptr%ptr
               n1 = MPI_Wtime()
               if (chara == 'N') then
                  call Redistribute1Dto1D(Vout_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Vout, M, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Nrnd, ptree, 1)
                  if (blocks%M_loc > 0) deallocate (Vout_locs(bb)%vector)
                  if (blocks%N_loc > 0) deallocate (Vin_locs(bb)%vector)
               else
                  call Redistribute1Dto1D(Vout_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Vout, N, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Nrnd, ptree, 1)
                  if (blocks%N_loc > 0) deallocate (Vout_locs(bb)%vector)
                  if (blocks%M_loc > 0) deallocate (Vin_locs(bb)%vector)
               endif
               n2 = MPI_Wtime()
               stats%Time_RedistV = stats%Time_RedistV + n2-n1
            end select
            cur => cur%next
         enddo

         deallocate(Vin_locs)
         deallocate(Vout_locs)
      end do

      random2(1:size(Vout,1),1:Nrnd) = random2(1:size(Vout,1),1:Nrnd)*b + Vout*a
      deallocate (Vout)

   end subroutine Hmat_parallelblock_MVP_dat

   subroutine BP_Mult(BP, chara, xin, xout, Ninloc, Noutloc, Ncol, ptree, stats)
      implicit none
      type(blockplus)::BP
      integer Ninloc, Noutloc, Ncol
      character chara
      type(proctree)::ptree
      type(Hstat)::stats
      DT::xin(Ninloc, Ncol), xout(Noutloc, Ncol)
      if (chara == 'N') then
         call Bplus_block_MVP_dat(BP, chara, Noutloc, Ninloc, Ncol, xin, Ninloc, xout, Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      else
         call Bplus_block_MVP_dat(BP, chara, Ninloc, Noutloc, Ncol, xin, Ninloc, xout, Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      endif
   end subroutine BP_Mult


   subroutine Bplus_block_MVP_dat(bplus, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats, level_start, level_end)



      implicit none

      integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b, ctemp1, ctemp2
      character chara
      type(matrixblock), pointer::blocks, blocks_1
      type(blockplus)::bplus
      integer ll, bb
      integer, optional:: level_start, level_end
      integer:: level_s, level_e
      type(proctree)::ptree
      type(Hstat)::stats
      integer:: ldi, ldo
      DT :: random1(ldi, *), random2(ldo, *)
      DT, allocatable :: Vout(:, :), Vin_loc(:, :), Vout_loc(:, :)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)
      real(kind=8)::n2,n1
      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)
      type(vectorsblock_oneL),allocatable::Vin_locs(:),Vout_locs(:)

      integer idx_start_m, idx_start_n, idx_start_n_loc, idx_start_m_loc, idx_end_n_loc, idx_end_m_loc, idx_start_i_loc, idx_start_o_loc, idx_end_i_loc, idx_end_o_loc

      blocks_1 => bplus%LL(1)%matrices_block(1)
      if(blocks_1%style==4)then
         call Hmat_parallelblock_MVP_dat(blocks_1, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
      else
         level_s = 1
         level_e = bplus%Lplus
         if (present(level_start)) level_s = level_start
         if (present(level_end)) level_e = level_end

         if (chara == 'N') allocate (Vout(M, Nrnd))
         if (chara == 'T') allocate (Vout(N, Nrnd))
         Vout = 0

         ctemp1 = 1.0d0; ctemp2 = 0.0d0


         allocate(Vin_locs(level_s:level_e))
         allocate(Vout_locs(level_s:level_e))

         do ll = level_s, level_e
            allocate(Vin_locs(ll)%vs(bplus%LL(ll)%Nbound))
            allocate(Vout_locs(ll)%vs(bplus%LL(ll)%Nbound))
            ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)
               if (chara == 'N') then
                  if (blocks%M_loc > 0) then
                     allocate (Vout_locs(ll)%vs(bb)%vector(blocks%M_loc, Nrnd))
                     Vout_locs(ll)%vs(bb)%vector=0
                  endif
                  if (blocks%N_loc > 0) allocate (Vin_locs(ll)%vs(bb)%vector(blocks%N_loc, Nrnd))
               else
                  if (blocks%N_loc > 0)then
                     allocate (Vout_locs(ll)%vs(bb)%vector(blocks%N_loc, Nrnd))
                     Vout_locs(ll)%vs(bb)%vector=0
                  endif
                  if (blocks%M_loc > 0) allocate (Vin_locs(ll)%vs(bb)%vector(blocks%M_loc, Nrnd))
               endif
            enddo
         enddo

         n1 = MPI_Wtime()
         if (chara == 'N') then
            call Bplus_vec_1Dto1D(bplus, 0, 1, level_s, level_e, ldi, random1, Nrnd, Vin_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         else
            call Bplus_vec_1Dto1D(bplus, 1, 1, level_s, level_e, ldi, random1, Nrnd, Vin_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         endif
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1


         do ll = level_s, level_e
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)

               if (blocks%N_loc > 0 .or. blocks%M_loc > 0) then
                  call BF_block_MVP_dat(blocks, chara, blocks%M_loc, blocks%N_loc, Nrnd,&
                  &Vin_locs(ll)%vs(bb)%vector, size(Vin_locs(ll)%vs(bb)%vector,1), Vout_locs(ll)%vs(bb)%vector, size(Vout_locs(ll)%vs(bb)%vector,1), ctemp1, ctemp2, ptree, stats)
               endif
            enddo
         enddo


         n1 = MPI_Wtime()
         if (chara == 'N') then
            call Bplus_vec_1Dto1D(bplus, 1, 0, level_s, level_e, M, Vout, Nrnd, Vout_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         else
            call Bplus_vec_1Dto1D(bplus, 0, 0, level_s, level_e, N, Vout, Nrnd, Vout_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         endif
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1


         do ll = level_s, level_e
            ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)
               if (blocks%M_loc > 0) then
                  deallocate (Vin_locs(ll)%vs(bb)%vector)
                  deallocate (Vout_locs(ll)%vs(bb)%vector)
               endif
            enddo
            deallocate(Vin_locs(ll)%vs)
            deallocate(Vout_locs(ll)%vs)
         enddo
         deallocate(Vin_locs)
         deallocate(Vout_locs)

         random2(1:size(Vout,1),1:Nrnd) = random2(1:size(Vout,1),1:Nrnd)*b + Vout*a
         deallocate (Vout)
      endif
   end subroutine Bplus_block_MVP_dat



   !!!!! this version is no more used as the communication is too slow
   ! subroutine Bplus_block_MVP_dat(bplus, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats, level_start, level_end)



   !    implicit none

   !    integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
   !    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
   !    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
   !    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
   !    integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
   !    DT ctemp, a, b, ctemp1, ctemp2
   !    character chara
   !    type(matrixblock), pointer::blocks, blocks_1
   !    type(blockplus)::bplus
   !    integer ll, bb
   !    integer, optional:: level_start, level_end
   !    integer:: level_s, level_e
   !    type(proctree)::ptree
   !    type(Hstat)::stats
   !    integer:: ldi, ldo
   !    DT :: random1(ldi, *), random2(ldo, *)
   !    DT, allocatable :: Vout(:, :), Vin_loc(:, :), Vout_loc(:, :)
   !    DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)
   !    real(kind=8)::n2,n1
   !    integer, allocatable:: arr_acc_m(:), arr_acc_n(:)
   !    type(vectorsblock),allocatable:: Vin_locs(:), Vout_locs(:)

   !    integer idx_start_m, idx_start_n, idx_start_n_loc, idx_start_m_loc, idx_end_n_loc, idx_end_m_loc, idx_start_i_loc, idx_start_o_loc, idx_end_i_loc, idx_end_o_loc

   !    blocks_1 => bplus%LL(1)%matrices_block(1)
   !    if(blocks_1%style==4)then
   !       call Hmat_parallelblock_MVP_dat(blocks_1, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
   !    else
   !       level_s = 1
   !       level_e = bplus%Lplus
   !       if (present(level_start)) level_s = level_start
   !       if (present(level_end)) level_e = level_end

   !       if (chara == 'N') allocate (Vout(M, Nrnd))
   !       if (chara == 'T') allocate (Vout(N, Nrnd))
   !       Vout = 0

   !       ctemp1 = 1.0d0; ctemp2 = 0.0d0

   !       do ll = level_s, level_e

   !          allocate(Vin_locs(bplus%LL(ll)%Nbound))
   !          allocate(Vout_locs(bplus%LL(ll)%Nbound))
   !          ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
   !          do bb = 1, bplus%LL(ll)%Nbound
   !             blocks => bplus%LL(ll)%matrices_block(bb)

   !             n1 = MPI_Wtime()
   !             if (chara == 'N') then
   !                if (blocks%M_loc > 0) then
   !                   allocate (Vout_locs(bb)%vector(blocks%M_loc, Nrnd))
   !                   Vout_locs(bb)%vector=0
   !                endif
   !                if (blocks%N_loc > 0) allocate (Vin_locs(bb)%vector(blocks%N_loc, Nrnd))
   !                call Redistribute1Dto1D(random1, ldi, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Vin_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Nrnd, ptree)
   !             else
   !                if (blocks%N_loc > 0)then
   !                   allocate (Vout_locs(bb)%vector(blocks%N_loc, Nrnd))
   !                   Vout_locs(bb)%vector=0
   !                endif
   !                if (blocks%M_loc > 0) allocate (Vin_locs(bb)%vector(blocks%M_loc, Nrnd))
   !                call Redistribute1Dto1D(random1, ldi, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Vin_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Nrnd, ptree)
   !             endif
   !             n2 = MPI_Wtime()
   !             stats%Time_RedistV = stats%Time_RedistV + n2-n1
   !          enddo

   !          do bb = 1, bplus%LL(ll)%Nbound
   !             blocks => bplus%LL(ll)%matrices_block(bb)

   !             if (blocks%N_loc > 0 .or. blocks%M_loc > 0) then
   !                call BF_block_MVP_dat(blocks, chara, blocks%M_loc, blocks%N_loc, Nrnd,&
   !                &Vin_locs(bb)%vector, size(Vin_locs(bb)%vector,1), Vout_locs(bb)%vector, size(Vout_locs(bb)%vector,1), ctemp1, ctemp2, ptree, stats)
   !             endif
   !          enddo

   !          ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
   !          do bb = 1, bplus%LL(ll)%Nbound
   !             blocks => bplus%LL(ll)%matrices_block(bb)
   !             n1 = MPI_Wtime()
   !             if (chara == 'N') then
   !                call Redistribute1Dto1D(Vout_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Vout, M, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Nrnd, ptree, 1)
   !                if (blocks%M_loc > 0) deallocate (Vout_locs(bb)%vector)
   !                if (blocks%N_loc > 0) deallocate (Vin_locs(bb)%vector)
   !             else
   !                call Redistribute1Dto1D(Vout_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Vout, N, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Nrnd, ptree, 1)
   !                if (blocks%N_loc > 0) deallocate (Vout_locs(bb)%vector)
   !                if (blocks%M_loc > 0) deallocate (Vin_locs(bb)%vector)
   !             endif
   !             n2 = MPI_Wtime()
   !             stats%Time_RedistV = stats%Time_RedistV + n2-n1

   !          end do
   !          deallocate(Vin_locs)
   !          deallocate(Vout_locs)
   !       end do

   !       random2(1:size(Vout,1),1:Nrnd) = random2(1:size(Vout,1),1:Nrnd)*b + Vout*a
   !       deallocate (Vout)
   !    endif
   ! end subroutine Bplus_block_MVP_dat



  subroutine Bplus_MD_block_MVP_dat(Ndim, bplus, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats, msh, option,level_start, level_end)



      implicit none

      integer Ndim
      integer M(Ndim), N(Ndim), Nrnd, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b, ctemp1, ctemp2
      character chara
      type(matrixblock_MD), pointer::blocks, blocks_1
      type(blockplus_MD)::bplus
      integer ll, bb
      integer, optional:: level_start, level_end
      integer:: level_s, level_e
      type(proctree)::ptree
      type(Hstat)::stats
      type(Hoption)::option
      type(mesh)::msh(Ndim)
      integer:: ldi(Ndim), ldo(Ndim)
      integer:: dim_in(Ndim), dim_out(Ndim)
      DT :: random1(product(ldi), *), random2(product(ldo), *)
      DT, allocatable :: Vout(:, :), Vin_loc(:, :), Vout_loc(:, :)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)
      real(kind=8)::n2,n1
      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)
      type(vectorsblock_oneL),allocatable::Vin_locs(:),Vout_locs(:)

      integer idx_start_m, idx_start_n, idx_start_n_loc, idx_start_m_loc, idx_end_n_loc, idx_end_m_loc, idx_start_i_loc, idx_start_o_loc, idx_end_i_loc, idx_end_o_loc

      blocks_1 => bplus%LL(1)%matrices_block(1)
      if(blocks_1%style==4)then
         write(*,*)"not supported style==4 in Bplus_MD_block_MVP_dat"
         stop
      else
         level_s = 1
         level_e = bplus%Lplus
         if (present(level_start)) level_s = level_start
         if (present(level_end)) level_e = level_end

         if (chara == 'N') allocate (Vout(product(M), Nrnd))
         if (chara == 'T') allocate (Vout(product(N), Nrnd))
         Vout = 0

         ctemp1 = 1.0d0; ctemp2 = 0.0d0


         allocate(Vin_locs(level_s:level_e))
         allocate(Vout_locs(level_s:level_e))

         do ll = level_s, level_e
            allocate(Vin_locs(ll)%vs(bplus%LL(ll)%Nbound))
            allocate(Vout_locs(ll)%vs(bplus%LL(ll)%Nbound))
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)
               if (chara == 'N') then
                  if (ALL(blocks%M_loc > 0)) then
                     allocate (Vout_locs(ll)%vs(bb)%vector(product(blocks%M_loc), Nrnd))
                     Vout_locs(ll)%vs(bb)%vector=0
                  endif
                  if (ALL(blocks%N_loc > 0)) allocate (Vin_locs(ll)%vs(bb)%vector(product(blocks%N_loc), Nrnd))
               else
                  if (ALL(blocks%N_loc > 0))then
                     allocate (Vout_locs(ll)%vs(bb)%vector(product(blocks%N_loc), Nrnd))
                     Vout_locs(ll)%vs(bb)%vector=0
                  endif
                  if (ALL(blocks%M_loc > 0)) allocate (Vin_locs(ll)%vs(bb)%vector(product(blocks%M_loc), Nrnd))
               endif
            enddo
         enddo

         n1 = MPI_Wtime()
         if (chara == 'N') then
            call Bplus_MD_vec_1Dto1D(Ndim, bplus, 0, 1, level_s, level_e, ldi, random1, Nrnd, Vin_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         else
            call Bplus_MD_vec_1Dto1D(Ndim, bplus, 1, 1, level_s, level_e, ldi, random1, Nrnd, Vin_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         endif
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1

         do ll = level_s, level_e
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)
               if (chara == 'N') then
                  dim_in = blocks%N_loc
                  dim_out = blocks%M_loc
               else
                  dim_in = blocks%M_loc
                  dim_out = blocks%N_loc
               endif

               if (ALL(blocks%N_loc > 0) .or. ALL(blocks%M_loc > 0)) then
                  if(blocks%style==1)then
                     call Full_block_MD_MVP_dat(blocks, chara, product(blocks%M_loc), Nrnd, Vin_locs(ll)%vs(bb)%vector, product(dim_in), Vout_locs(ll)%vs(bb)%vector, product(dim_out), ctemp1, ctemp2)
                  else
                     call BF_MD_block_mvp(chara, Vin_locs(ll)%vs(bb)%vector, dim_in, Vout_locs(ll)%vs(bb)%vector, dim_out, Nrnd, blocks, Ndim, ptree, stats,msh,option)
                  endif
               endif
            enddo
         enddo


         n1 = MPI_Wtime()
         if (chara == 'N') then
            call Bplus_MD_vec_1Dto1D(Ndim, bplus, 1, 0, level_s, level_e, M, Vout, Nrnd, Vout_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         else
            call Bplus_MD_vec_1Dto1D(Ndim, bplus, 0, 0, level_s, level_e, N, Vout, Nrnd, Vout_locs, ptree,ptree%pgrp(blocks_1%pgno)%nproc)
         endif
         n2 = MPI_Wtime()
         stats%Time_RedistV = stats%Time_RedistV + n2-n1


         do ll = level_s, level_e
            do bb = 1, bplus%LL(ll)%Nbound
               blocks => bplus%LL(ll)%matrices_block(bb)
               if (ALL(blocks%M_loc > 0))then
                  deallocate (Vout_locs(ll)%vs(bb)%vector)
                  deallocate (Vin_locs(ll)%vs(bb)%vector)
               endif
            end do
            deallocate(Vin_locs(ll)%vs)
            deallocate(Vout_locs(ll)%vs)
         end do
         deallocate(Vin_locs)
         deallocate(Vout_locs)

         random2(1:size(Vout,1),1:Nrnd) = random2(1:size(Vout,1),1:Nrnd)*b + Vout*a
         deallocate (Vout)
      endif
   end subroutine Bplus_MD_block_MVP_dat



   !!!!! this version is no more used as the communication is too slow
!   subroutine Bplus_MD_block_MVP_dat(Ndim, bplus, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats, msh, option,level_start, level_end)



!       implicit none

!       integer Ndim
!       integer M(Ndim), N(Ndim), Nrnd, index_i, index_j, na, nb, index_start, num_vectors
!       integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
!       integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
!       integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
!       integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
!       DT ctemp, a, b, ctemp1, ctemp2
!       character chara
!       type(matrixblock_MD), pointer::blocks, blocks_1
!       type(blockplus_MD)::bplus
!       integer ll, bb
!       integer, optional:: level_start, level_end
!       integer:: level_s, level_e
!       type(proctree)::ptree
!       type(Hstat)::stats
!       type(Hoption)::option
!       type(mesh)::msh(Ndim)
!       integer:: ldi(Ndim), ldo(Ndim)
!       integer:: dim_in(Ndim), dim_out(Ndim)
!       DT :: random1(product(ldi), *), random2(product(ldo), *)
!       DT, allocatable :: Vout(:, :), Vin_loc(:, :), Vout_loc(:, :)
!       DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)
!       real(kind=8)::n2,n1
!       integer, allocatable:: arr_acc_m(:), arr_acc_n(:)
!       type(vectorsblock),allocatable:: Vin_locs(:), Vout_locs(:)

!       integer idx_start_m, idx_start_n, idx_start_n_loc, idx_start_m_loc, idx_end_n_loc, idx_end_m_loc, idx_start_i_loc, idx_start_o_loc, idx_end_i_loc, idx_end_o_loc

!       blocks_1 => bplus%LL(1)%matrices_block(1)
!       if(blocks_1%style==4)then
!          write(*,*)"not supported style==4 in Bplus_MD_block_MVP_dat"
!          stop
!       else
!          level_s = 1
!          level_e = bplus%Lplus
!          if (present(level_start)) level_s = level_start
!          if (present(level_end)) level_e = level_end

!          if (chara == 'N') allocate (Vout(product(M), Nrnd))
!          if (chara == 'T') allocate (Vout(product(N), Nrnd))
!          Vout = 0

!          ctemp1 = 1.0d0; ctemp2 = 0.0d0

!          do ll = level_s, level_e

!             allocate(Vin_locs(bplus%LL(ll)%Nbound))
!             allocate(Vout_locs(bplus%LL(ll)%Nbound))
!             ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
!             do bb = 1, bplus%LL(ll)%Nbound
!                blocks => bplus%LL(ll)%matrices_block(bb)

!                n1 = MPI_Wtime()
!                if (chara == 'N') then
!                   if (ALL(blocks%M_loc > 0)) then
!                      allocate (Vout_locs(bb)%vector(product(blocks%M_loc), Nrnd))
!                      Vout_locs(bb)%vector=0
!                   endif
!                   if (ALL(blocks%N_loc > 0)) allocate (Vin_locs(bb)%vector(product(blocks%N_loc), Nrnd))
!                   call Redistribute1Dto1D_MD(Ndim, random1, ldi, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Vin_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Nrnd, ptree)
!                else
!                   if (ALL(blocks%N_loc > 0))then
!                      allocate (Vout_locs(bb)%vector(product(blocks%N_loc), Nrnd))
!                      Vout_locs(bb)%vector=0
!                   endif
!                   if (ALL(blocks%M_loc > 0)) allocate (Vin_locs(bb)%vector(product(blocks%M_loc), Nrnd))
!                   call Redistribute1Dto1D_MD(Ndim, random1, ldi, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Vin_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Nrnd, ptree)
!                endif
!                n2 = MPI_Wtime()
!                stats%Time_RedistV = stats%Time_RedistV + n2-n1
!             enddo

!             do bb = 1, bplus%LL(ll)%Nbound
!                blocks => bplus%LL(ll)%matrices_block(bb)
!                if (chara == 'N') then
!                   dim_in = blocks%N_loc
!                   dim_out = blocks%M_loc
!                else
!                   dim_in = blocks%M_loc
!                   dim_out = blocks%N_loc
!                endif

!                if (ALL(blocks%N_loc > 0) .or. ALL(blocks%M_loc > 0)) then
!                   if(blocks%style==1)then
!                      call Full_block_MD_MVP_dat(blocks, chara, product(blocks%M_loc), Nrnd, Vin_locs(bb)%vector, product(dim_in), Vout_locs(bb)%vector, product(dim_out), ctemp1, ctemp2)
!                   else
!                      call BF_MD_block_mvp(chara, Vin_locs(bb)%vector, dim_in, Vout_locs(bb)%vector, dim_out, Nrnd, blocks, Ndim, ptree, stats,msh,option)
!                   endif
!                endif
!             enddo

!             ! The followin loop for input redistrbution should be rewritten in a more efficient way rather than using Redistribute1Dto1D
!             do bb = 1, bplus%LL(ll)%Nbound
!                blocks => bplus%LL(ll)%matrices_block(bb)
!                n1 = MPI_Wtime()
!                if (chara == 'N') then
!                   call Redistribute1Dto1D_MD(Ndim,Vout_locs(bb)%vector, blocks%M_loc, blocks%M_p, blocks%headm, blocks%pgno, Vout, M, blocks_1%M_p, blocks_1%headm, blocks_1%pgno, Nrnd, ptree, 1)
!                   if (ALL(blocks%M_loc > 0)) deallocate (Vout_locs(bb)%vector)
!                   if (ALL(blocks%N_loc > 0)) deallocate (Vin_locs(bb)%vector)
!                else
!                   call Redistribute1Dto1D_MD(Ndim,Vout_locs(bb)%vector, blocks%N_loc, blocks%N_p, blocks%headn, blocks%pgno, Vout, N, blocks_1%N_p, blocks_1%headn, blocks_1%pgno, Nrnd, ptree, 1)
!                   if (ALL(blocks%N_loc > 0)) deallocate (Vout_locs(bb)%vector)
!                   if (ALL(blocks%M_loc > 0)) deallocate (Vin_locs(bb)%vector)
!                endif
!                n2 = MPI_Wtime()
!                stats%Time_RedistV = stats%Time_RedistV + n2-n1

!             end do
!             deallocate(Vin_locs)
!             deallocate(Vout_locs)
!          end do

!          random2(1:size(Vout,1),1:Nrnd) = random2(1:size(Vout,1),1:Nrnd)*b + Vout*a
!          deallocate (Vout)
!       endif
!    end subroutine Bplus_MD_block_MVP_dat





! redistribute Bplus
   subroutine Bplus_ReDistribute_Inplace(bplus_o, stats, ptree, msh)
      implicit none

      integer nproc_i, nproc_o, idxs_i, idxs_o, idxe_i, idxe_o, ii, jj, iii, jjj, level
      type(proctree)::ptree
      type(mesh)::msh
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o, rank, rankmax
      type(blockplus)::bplus_o
      type(matrixblock), pointer::blocks
      DT, pointer::dat_new(:, :), dat_old(:, :)
      real(kind=8)::n1, n2
      type(Hstat)::stats
      integer pgno, ll, bb

      if (associated(bplus_o%LL)) then
      do ll = 1, LplusMax
         if (bplus_o%LL(ll)%Nbound > 0) then
            if (associated(bplus_o%LL(ll)%matrices_block)) then
            do bb = 1, bplus_o%LL(ll)%Nbound
               pgno = bplus_o%LL(ll)%matrices_block(bb)%pgno_db
               if (IOwnPgrp(ptree, pgno)) call BF_ReDistribute_Inplace(bplus_o%LL(ll)%matrices_block(bb), pgno, stats, ptree, msh)
            end do
            endif
            call MPI_ALLREDUCE(MPI_IN_PLACE, bplus_o%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(bplus_o%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
         end if
      end do
      endif

   end subroutine Bplus_ReDistribute_Inplace

! redistribute Butterfly
   subroutine BF_ReDistribute_Inplace(blocks, pgno_new, stats, ptree, msh)
      implicit none

      integer nproc_i, nproc_o, idxs_i, idxs_o, idxe_i, idxe_o, ii, jj, iii, jjj, level
      type(proctree)::ptree
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o, rank, rankmax, pgno_new
      type(matrixblock)::blocks, blocks_dummy
      DT, pointer::dat_new(:, :), dat_old(:, :)
      real(kind=8)::n1, n2
      type(Hstat)::stats
      type(mesh)::msh

      dat_new => null()
      dat_old => null()

      ! call MPI_barrier(ptree%pgrp(pgno_new)%Comm,ierr)
      n1 = MPI_Wtime()

      if (blocks%level_butterfly == 0) then

         if (blocks%pgno /= pgno_new) then
            ! communicate block sizes first
            if (blocks%M_loc > 0) then
               rank = blocks%rankmax
            else
               rank = 0
            endif
            call assert(MPI_COMM_NULL /= ptree%pgrp(pgno_new)%Comm, 'communicator should not be null 4')
            rankmax = 0
            call MPI_ALLREDUCE(rank, rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno_new)%Comm, ierr)
            rank = rankmax

            blocks_dummy%level = blocks%level
            blocks_dummy%M = blocks%M
            blocks_dummy%N = blocks%N
            blocks_dummy%row_group = blocks%row_group
            blocks_dummy%col_group = blocks%col_group
            call ComputeParallelIndices(blocks_dummy, pgno_new, ptree, msh)

            ! redistribute U
            allocate (dat_old(max(1, blocks%M_loc), rank))
            dat_old = 0
            if (blocks%M_loc > 0) dat_old = blocks%ButterflyU%blocks(1)%matrix

            allocate (dat_new(max(1, blocks_dummy%M_loc), rank))
            dat_new = 0

            call Redistribute1Dto1D(dat_old, max(1, blocks%M_loc), blocks%M_p, 0, blocks%pgno, dat_new, max(1, blocks_dummy%M_loc), blocks_dummy%M_p, 0, pgno_new, rank, ptree)
            if (blocks%M_loc > 0) then
               deallocate (blocks%ButterflyU%blocks(1)%matrix)

            endif
            deallocate (dat_old)
            if (blocks_dummy%M_loc > 0) then
               if (.not. allocated(blocks%ButterflyU%blocks)) allocate (blocks%ButterflyU%blocks(1))
               if (.not. associated(blocks%ButterflyU%blocks(1)%matrix)) allocate (blocks%ButterflyU%blocks(1)%matrix(blocks_dummy%M_loc, rank))
               blocks%ButterflyU%blocks(1)%matrix = dat_new

            endif
            deallocate (dat_new)

            ! redistribute V
            allocate (dat_old(max(1, blocks%N_loc), rank))
            dat_old = 0
            if (blocks%N_loc > 0) dat_old = blocks%ButterflyV%blocks(1)%matrix

            allocate (dat_new(max(1, blocks_dummy%N_loc), rank))
            dat_new = 0

            call Redistribute1Dto1D(dat_old,max(1, blocks%N_loc), blocks%N_p, 0, blocks%pgno, dat_new, max(1, blocks_dummy%N_loc),blocks_dummy%N_p, 0, pgno_new, rank, ptree)
            if (blocks%N_loc > 0) then
               deallocate (blocks%ButterflyV%blocks(1)%matrix)

            endif
            deallocate (dat_old)

            if (blocks_dummy%N_loc > 0) then
               if (.not. allocated(blocks%ButterflyV%blocks)) allocate (blocks%ButterflyV%blocks(1))
               if (.not. associated(blocks%ButterflyV%blocks(1)%matrix)) allocate (blocks%ButterflyV%blocks(1)%matrix(blocks_dummy%N_loc, rank))
               blocks%ButterflyV%blocks(1)%matrix = dat_new

               blocks%rankmax = rank
            endif
            deallocate (dat_new)

            blocks%pgno = pgno_new
            if (IOwnPgrp(ptree, blocks%pgno)) call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)

            deallocate (blocks_dummy%M_p)
            deallocate (blocks_dummy%N_p)

         endif
      else
         if (blocks%pgno /= pgno_new) then

            !>*** make sure every process has blocks%ButterflyKerl allocated
            if (.not. allocated(blocks%ButterflyKerl)) then
               allocate (blocks%ButterflyKerl(blocks%level_butterfly))
            endif

            do level = 0, blocks%level_butterfly + 1
               if (level == 0) then
                  call BF_all2all_UV(blocks, blocks%pgno, blocks%ButterflyV, level, 0, blocks, pgno_new, blocks%ButterflyV, level, stats, ptree)
               elseif (level == blocks%level_butterfly + 1) then
                  call BF_all2all_UV(blocks, blocks%pgno, blocks%ButterflyU, level, 0, blocks, pgno_new, blocks%ButterflyU, level, stats, ptree)
               else
                  call BF_all2all_ker(blocks, blocks%pgno, blocks%ButterflyKerl(level), level, 0, 0, blocks, pgno_new, blocks%ButterflyKerl(level), level, stats, ptree)
               endif
            enddo

            !>*** delete dummy blocks%ButterflyKerl if I don't share the output butterfly
            if (.not. IOwnPgrp(ptree, pgno_new)) then
               deallocate (blocks%ButterflyKerl)
            endif

            blocks%pgno = pgno_new
            if (IOwnPgrp(ptree, blocks%pgno)) call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)
            call BF_get_rank(blocks, ptree)

         endif
      endif

      n2 = MPI_Wtime()
      stats%Time_RedistB = stats%Time_RedistB + n2 - n1

   end subroutine BF_ReDistribute_Inplace



! change the pattern of a butterfly
   subroutine BF_ChangePattern(blocks, pat_i, pat_o, stats, ptree)
      implicit none

      integer pat_i, pat_o,level
      type(proctree)::ptree
      integer ierr
      type(matrixblock)::blocks
      real(kind=8)::n1, n2
      type(Hstat)::stats

      ! call MPI_barrier(ptree%pgrp(pgno_new)%Comm,ierr)
      n1 = MPI_Wtime()

      if (blocks%level_butterfly > 0) then
         if (pat_i /= pat_o) then

            !>*** make sure every process has blocks%ButterflyKerl allocated
            if (.not. allocated(blocks%ButterflyKerl)) then
               allocate (blocks%ButterflyKerl(blocks%level_butterfly))
            endif

            do level = 1, blocks%level_butterfly
               call BF_all2all_ker_pattern(blocks, blocks%ButterflyKerl(level), pat_i, blocks,blocks%ButterflyKerl(level), pat_o, level, blocks%pgno,stats, ptree)
            enddo
            blocks%level_half=BF_Switchlevel(blocks%level_butterfly, pat_o)


            !>*** delete dummy blocks%ButterflyKerl if I don't share the output butterfly
            if (.not. IOwnPgrp(ptree, blocks%pgno)) then
               deallocate (blocks%ButterflyKerl)
            endif

         endif
      endif

      n2 = MPI_Wtime()
      stats%Time_RedistB = stats%Time_RedistB + n2 - n1

   end subroutine BF_ChangePattern



   subroutine BF_delete(blocks, allflag)


      implicit none

      integer butterflyB_inuse, num_col, num_row
      integer i, j, mm, nn, rank, num_blocks, level, level_butterfly, index_i_m, index_j_m, levelm
      real(kind=8) memory_butterfly, rtemp
      type(matrixblock)::blocks
      integer allflag

      level_butterfly = blocks%level_butterfly

      if (allocated(blocks%ButterflyU%blocks)) then
         ! !$omp parallel do default(shared) private(i)
         do i = 1, blocks%ButterflyU%nblk_loc
            if (associated(blocks%ButterflyU%blocks(i)%matrix)) deallocate (blocks%ButterflyU%blocks(i)%matrix)
         enddo
         ! !$omp end parallel do
         deallocate (blocks%ButterflyU%blocks)
      end if

      if (allocated(blocks%ButterflyV%blocks)) then
         ! !$omp parallel do default(shared) private(i)
         do i = 1, blocks%ButterflyV%nblk_loc
            if (associated(blocks%ButterflyV%blocks(i)%matrix)) deallocate (blocks%ButterflyV%blocks(i)%matrix)
         enddo
         ! !$omp end parallel do
         deallocate (blocks%ButterflyV%blocks)
      end if

      if (allocated(blocks%ButterflyKerl)) then
         if (level_butterfly /= 0) then
            ! !$omp parallel do default(shared) private(level,i,j,num_col,num_row)
            do level = 1, level_butterfly
               if (allocated(blocks%ButterflyKerl(level)%blocks)) then
               do j = 1, blocks%ButterflyKerl(level)%nc
                  do i = 1, blocks%ButterflyKerl(level)%nr
                     if (associated(blocks%ButterflyKerl(level)%blocks(i, j)%matrix))then
                        deallocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix)
                     endif
                     enddo
               enddo
               deallocate (blocks%ButterflyKerl(level)%blocks)
               endif
            enddo
            ! !$omp end parallel do
         endif
         deallocate (blocks%ButterflyKerl)
      end if

      if (allocated(blocks%ButterflyMiddle)) then
         write (*, *) 'warning: this part has not been parallelized in BF_delete'
         levelm = ceiling_safe(dble(level_butterfly)/2d0)
         do index_i_m = 1, 2**levelm
            do index_j_m = 1, 2**(level_butterfly - levelm)
               if (associated(blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix)) deallocate (blocks%ButterflyMiddle(index_i_m, index_j_m)%matrix)
            end do
         end do
         deallocate (blocks%ButterflyMiddle)
      end if

      ! blocks%level_butterfly=0
      blocks%rankmax = -1000
      blocks%rankmin = 1000

      if (associated(blocks%fullmat)) deallocate (blocks%fullmat)
#if HAVE_ZFP
      if (allocated(blocks%FullmatZFP%buffer_r)) then
         deallocate (blocks%FullmatZFP%buffer_r)
         call zFORp_stream_close(blocks%FullmatZFP%stream_r)
      endif
      if (allocated(blocks%FullmatZFP%buffer_i))then
         deallocate (blocks%FullmatZFP%buffer_i)
         call zFORp_stream_close(blocks%FullmatZFP%stream_i)
      endif
#endif
      if (allocated(blocks%fullmat_MPI)) deallocate (blocks%fullmat_MPI)
      if (allocated(blocks%ipiv)) deallocate (blocks%ipiv)
      if (allocated(blocks%Butterfly_data_MPI)) deallocate (blocks%Butterfly_data_MPI)
      if (allocated(blocks%Butterfly_index_MPI)) deallocate (blocks%Butterfly_index_MPI)
      if (associated(blocks%ms)) deallocate (blocks%ms)
      if (associated(blocks%ns)) deallocate (blocks%ns)

      if (allocated(blocks%lstblks)) then
         do level = 0, size(blocks%lstblks)-1
            call list_finalizer(blocks%lstblks(level))
         enddo
         deallocate (blocks%lstblks)
      endif

      if (allflag == 1) then
         if (associated(blocks%N_p)) deallocate (blocks%N_p)
         if (associated(blocks%M_p)) deallocate (blocks%M_p)
      endif
      return

   end subroutine BF_delete

   subroutine BF_copy(trans, block_i, block_o, memory)


      implicit none
      type(matrixblock)::block_i, block_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      character::trans
      real(kind=8), optional::memory
      real(kind=8)::tol_used

      if (present(memory)) memory = 0

      if (trans == 'N') then

         block_o%level = block_i%level
         block_o%col_group = block_i%col_group
         block_o%row_group = block_i%row_group
         block_o%style = block_i%style
         block_o%level_butterfly = block_i%level_butterfly
         block_o%level_half = block_i%level_half
         block_o%rankmax = block_i%rankmax
         block_o%rankmin = block_i%rankmin
         block_o%dimension_rank = block_i%dimension_rank
         block_o%M = block_i%M
         block_o%N = block_i%N
         block_o%headm = block_i%headm
         block_o%headn = block_i%headn

         block_o%M_loc = block_i%M_loc
         block_o%N_loc = block_i%N_loc
         block_o%pgno = block_i%pgno
         block_o%pgno_db = block_i%pgno_db

         if (associated(block_i%N_p)) then
            if (associated(block_o%N_p)) deallocate (block_o%N_p)
            allocate (block_o%N_p(size(block_i%N_p, 1), 2))
            if (present(memory)) memory = memory + SIZEOF(block_o%N_p)/1024.0d3
            block_o%N_p = block_i%N_p
         endif
         if (associated(block_i%M_p)) then
            if (associated(block_o%M_p)) deallocate (block_o%M_p)
            allocate (block_o%M_p(size(block_i%M_p, 1), 2))
            if (present(memory)) memory = memory + SIZEOF(block_o%M_p)/1024.0d3
            block_o%M_p = block_i%M_p
         endif
         if (associated(block_i%ns)) then
            if (associated(block_o%ns)) deallocate (block_o%ns)
            allocate (block_o%ns(size(block_i%ns, 1)))
            if (present(memory)) memory = memory + SIZEOF(block_o%ns)/1024.0d3
            block_o%ns = block_i%ns
         endif
         if (associated(block_i%ms)) then
            if (associated(block_o%ms)) deallocate (block_o%ms)
            allocate (block_o%ms(size(block_i%ms, 1)))
            if (present(memory)) memory = memory + SIZEOF(block_o%ms)/1024.0d3
            block_o%ms = block_i%ms
         endif


         level_butterfly = block_i%level_butterfly
         num_blocks = 2**level_butterfly

         if (block_i%style == 2) then
            if (allocated(block_i%ButterflyU%blocks)) then
               if (level_butterfly /= 0) then
                  allocate (block_o%ButterflyKerl(level_butterfly))
               end if

               do level = 0, level_butterfly + 1
                  if (level == 0) then
                     block_o%ButterflyV%num_blk = block_i%ButterflyV%num_blk
                     block_o%ButterflyV%nblk_loc = block_i%ButterflyV%nblk_loc
                     block_o%ButterflyV%idx = block_i%ButterflyV%idx
                     block_o%ButterflyV%inc = block_i%ButterflyV%inc
                     allocate (block_o%ButterflyV%blocks(block_o%ButterflyV%nblk_loc))

                     do jj = 1, block_o%ButterflyV%nblk_loc
                        nn = size(block_i%ButterflyV%blocks(jj)%matrix, 1)
                        rank = size(block_i%ButterflyV%blocks(jj)%matrix, 2)
                        allocate (block_o%ButterflyV%blocks(jj)%matrix(nn, rank))
                        block_o%ButterflyV%blocks(jj)%matrix = block_i%ButterflyV%blocks(jj)%matrix
                        if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyV%blocks(jj)%matrix)/1024.0d3
                     enddo

                  else if (level == level_butterfly + 1) then
                     block_o%ButterflyU%num_blk = block_i%ButterflyU%num_blk
                     block_o%ButterflyU%nblk_loc = block_i%ButterflyU%nblk_loc
                     block_o%ButterflyU%idx = block_i%ButterflyU%idx
                     block_o%ButterflyU%inc = block_i%ButterflyU%inc
                     allocate (block_o%ButterflyU%blocks(block_o%ButterflyU%nblk_loc))

                     do ii = 1, block_o%ButterflyU%nblk_loc
                        nn = size(block_i%ButterflyU%blocks(ii)%matrix, 1)
                        rank = size(block_i%ButterflyU%blocks(ii)%matrix, 2)
                        allocate (block_o%ButterflyU%blocks(ii)%matrix(nn, rank))
                        block_o%ButterflyU%blocks(ii)%matrix = block_i%ButterflyU%blocks(ii)%matrix
                        if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyU%blocks(ii)%matrix)/1024.0d3
                     enddo
                  else
                     block_o%ButterflyKerl(level)%num_row = block_i%ButterflyKerl(level)%num_row
                     block_o%ButterflyKerl(level)%num_col = block_i%ButterflyKerl(level)%num_col
                     block_o%ButterflyKerl(level)%nc = block_i%ButterflyKerl(level)%nc
                     block_o%ButterflyKerl(level)%nr = block_i%ButterflyKerl(level)%nr
                     block_o%ButterflyKerl(level)%idx_c = block_i%ButterflyKerl(level)%idx_c
                     block_o%ButterflyKerl(level)%idx_r = block_i%ButterflyKerl(level)%idx_r
                     block_o%ButterflyKerl(level)%inc_c = block_i%ButterflyKerl(level)%inc_c
                     block_o%ButterflyKerl(level)%inc_r = block_i%ButterflyKerl(level)%inc_r
                     if(block_o%ButterflyKerl(level)%nr>0 .and. block_o%ButterflyKerl(level)%nc>0)then
                        allocate (block_o%ButterflyKerl(level)%blocks(block_o%ButterflyKerl(level)%nr, block_o%ButterflyKerl(level)%nc))
                        do ii = 1, block_o%ButterflyKerl(level)%nr
                           do jj = 1, block_o%ButterflyKerl(level)%nc
                              nn = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 2)
                              rank = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 1)
                              allocate (block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix(rank, nn))
                              block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix = block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix
                              if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix)/1024.0d3
                           enddo
                        enddo
                     endif
                  endif
               enddo
            endif

            if (allocated(block_i%ButterflyMiddle)) then
               write (*, *) 'ButterflyMiddle is not yet handled when the butterfly is distributed'
               levelm = ceiling_safe(dble(level_butterfly)/2d0)
               allocate (block_o%ButterflyMiddle(2**levelm, 2**(level_butterfly - levelm)))

               do index_i_m = 1, 2**levelm
                  do index_j_m = 1, 2**(level_butterfly - levelm)
                     rank = size(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, 1)
                     allocate (block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix(rank, rank))
                     block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix
                  end do
               end do
            end if

         else if (block_i%style == 1) then
            if (associated(block_i%fullmat)) then
               mm = size(block_i%fullmat, 1)
               nn = size(block_i%fullmat, 2)
               allocate (block_o%fullmat(mm, nn))
               block_o%fullmat = block_i%fullmat
               if(allocated(block_i%ipiv))then
                  allocate(block_o%ipiv(size(block_i%ipiv,1)))
                  block_o%ipiv=block_i%ipiv
               endif
               if (present(memory)) memory = memory + SIZEOF(block_o%fullmat)/1024.0d3
            endif
#if HAVE_ZFP
            if(allocated(block_i%FullmatZFP%buffer_r))then
               call ZFP_Decompress(block_i%fullmat,block_i%FullmatZFP,block_i%M,block_i%N,tol_used,1)
               mm = size(block_i%fullmat, 1)
               nn = size(block_i%fullmat, 2)
               allocate (block_o%fullmat(mm, nn))
               block_o%fullmat = block_i%fullmat
               call ZFP_Compress(block_i%fullmat,block_i%FullmatZFP,block_i%M,block_i%N,tol_used,1)
               call ZFP_Compress(block_o%fullmat,block_o%FullmatZFP,block_o%M,block_o%N,tol_used,0)
               if (present(memory)) memory = memory + SIZEOF(block_o%FullmatZFP%buffer_r)/1024.0d3
               if (present(memory) .and. allocated(block_o%FullmatZFP%buffer_i)) memory = memory + SIZEOF(block_o%FullmatZFP%buffer_i)/1024.0d3

               ! mm = size(block_i%FullmatZFP%buffer_r,1)
               ! allocate(block_o%FullmatZFP%buffer_r(mm))
               ! block_o%FullmatZFP%buffer_r = block_i%FullmatZFP%buffer_r
               ! block_o%FullmatZFP%stream_r = block_i%FullmatZFP%stream_r
               ! if (present(memory)) memory = memory + SIZEOF(block_o%FullmatZFP%buffer_r)/1024.0d3
            endif
            ! if(allocated(block_i%FullmatZFP%buffer_i))then
            !    mm = size(block_i%FullmatZFP%buffer_i,1)
            !    allocate(block_o%FullmatZFP%buffer_i(mm))
            !    block_o%FullmatZFP%buffer_i = block_i%FullmatZFP%buffer_i
            !    block_o%FullmatZFP%FullmatZFP%stream_i = block_i%FullmatZFP%stream_i
            !    if (present(memory)) memory = memory + SIZEOF(block_o%FullmatZFP%buffer_i)/1024.0d3
            ! endif
#endif
         else
            ! write(*,*)'block style not implemented'
            ! stop
         end if
      else if (trans == 'T') then
         write (*, *) 'transposed copy is not well tested if the butterfly is distributed'
         stop

         block_o%level = block_i%level
         block_o%col_group = block_i%row_group
         block_o%row_group = block_i%col_group
         block_o%style = block_i%style
         block_o%level_butterfly = block_i%level_butterfly
         block_o%rankmax = block_i%rankmax
         block_o%rankmin = block_i%rankmin
         block_o%dimension_rank = block_i%dimension_rank
         block_o%M = block_i%N
         block_o%N = block_i%M
         block_o%headm = block_i%headn
         block_o%headn = block_i%headm

         block_o%M_loc = block_i%N_loc
         block_o%N_loc = block_i%M_loc
         block_o%pgno = block_i%pgno
         block_o%pgno_db = block_i%pgno_db

         if (associated(block_i%N_p)) then
            if (associated(block_o%M_p)) deallocate (block_o%M_p)
            allocate (block_o%M_p(size(block_i%N_p, 1), 2))
            if (present(memory)) memory = memory + SIZEOF(block_o%M_p)/1024.0d3
            block_o%M_p = block_i%N_p
         endif
         if (associated(block_i%M_p)) then
            if (associated(block_o%N_p)) deallocate (block_o%N_p)
            allocate (block_o%N_p(size(block_i%M_p, 1), 2))
            if (present(memory)) memory = memory + SIZEOF(block_o%N_p)/1024.0d3
            block_o%N_p = block_i%M_p
         endif
         if (associated(block_i%ns)) then
            if (associated(block_o%ms)) deallocate (block_o%ms)
            allocate (block_o%ms(size(block_i%ns, 1)))
            if (present(memory)) memory = memory + SIZEOF(block_o%ms)/1024.0d3
            block_o%ms = block_i%ns
         endif
         if (associated(block_i%ms)) then
            if (associated(block_o%ns)) deallocate (block_o%ns)
            allocate (block_o%ns(size(block_i%ms, 1)))
            if (present(memory)) memory = memory + SIZEOF(block_o%ns)/1024.0d3
            block_o%ns = block_i%ms
         endif



         level_butterfly = block_i%level_butterfly
         num_blocks = 2**level_butterfly

         if (block_i%style == 2) then
            if (allocated(block_i%ButterflyU%blocks)) then
               allocate (block_o%ButterflyU%blocks(num_blocks))
               allocate (block_o%ButterflyV%blocks(num_blocks))
               if (level_butterfly /= 0) then
                  allocate (block_o%ButterflyKerl(level_butterfly))
               end if

               do level = 0, level_butterfly
                  index_ij = 0
                  if (level > 0) then
                     block_o%ButterflyKerl(level)%num_col = 2**level
                     block_o%ButterflyKerl(level)%num_row = 2**(level_butterfly - level + 1)
                     allocate (block_o%ButterflyKerl(level)%blocks(2**(level_butterfly - level + 1), 2**level))
                  endif
                  do index_i = 1, 2**level
                     do index_j = 1, 2**(level_butterfly - level)
                        index_ij = index_ij + 1
                        if (level == 0) then
                           nn = size(block_i%ButterflyV%blocks(index_ij)%matrix, 1)
                           rank = size(block_i%ButterflyV%blocks(index_ij)%matrix, 2)
                           allocate (block_o%ButterflyU%blocks(index_ij)%matrix(nn, rank))
                           block_o%ButterflyU%blocks(index_ij)%matrix = block_i%ButterflyV%blocks(index_ij)%matrix
                           if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyU%blocks(index_ij)%matrix)/1024.0d3
                        else
                           nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 2)
                           rank = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 1)
                           allocate (block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j - 1, index_i)%matrix(nn, rank))
                           call copymatT(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j - 1, index_i)%matrix, rank, nn)
                           if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j - 1, index_i)%matrix)/1024.0d3
                           nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, 2)
                           allocate (block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j, index_i)%matrix(nn, rank))
                           call copymatT(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j, index_i)%matrix, rank, nn)
                           if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly - level + 1)%blocks(2*index_j, index_i)%matrix)/1024.0d3
                        endif
                        if (level == level_butterfly) then
                           mm = size(block_i%ButterflyU%blocks(index_ij)%matrix, 1)
                           rank = size(block_i%ButterflyU%blocks(index_ij)%matrix, 2)
                           allocate (block_o%ButterflyV%blocks(index_ij)%matrix(mm, rank))
                           block_o%ButterflyV%blocks(index_ij)%matrix = block_i%ButterflyU%blocks(index_ij)%matrix
                           if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyV%blocks(index_ij)%matrix)/1024.0d3
                        endif
                     enddo
                  enddo
               enddo
            endif

            if (allocated(block_i%ButterflyMiddle)) then
               levelm = ceiling_safe(dble(level_butterfly)/2d0)
               allocate (block_o%ButterflyMiddle(2**(level_butterfly - levelm), 2**levelm))

               do index_i_m = 1, 2**levelm
                  do index_j_m = 1, 2**(level_butterfly - levelm)
                     rank = size(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, 1)
                     allocate (block_o%ButterflyMiddle(index_j_m, index_i_m)%matrix(rank, rank))
                     call copymatT(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, block_o%ButterflyMiddle(index_j_m, index_i_m)%matrix, rank, rank)
                  end do
               end do
            end if
         else if (block_i%style == 1) then
            if (associated(block_i%fullmat)) then
               mm = size(block_i%fullmat, 1)
               nn = size(block_i%fullmat, 2)
               allocate (block_o%fullmat(nn, mm))
               call copymatT(block_i%fullmat, block_o%fullmat, mm, nn)
               if(allocated(block_i%ipiv))then
                  allocate(block_o%ipiv(size(block_i%ipiv,1)))
                  block_o%ipiv=block_i%ipiv
               endif
               if (present(memory)) memory = memory + SIZEOF(block_o%fullmat)/1024.0d3
            endif
#if HAVE_ZFP
         write(*,*)'not yet impleted for ZFP data. Can be implemented: copy ZFP data, decompress, transpose, recompress'
         stop
#endif

         else
            ! write(*,*)'block style not implemented'
            ! stop
         end if
      endif

   end subroutine BF_copy


   subroutine BF_delete_ker_onelevel(ker_o)

      implicit none
      type(butterfly_kerl)::ker_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      character::trans

      if(ker_o%nr>0 .and. ker_o%nc>0 .and. allocated(ker_o%blocks))then
         do ii = 1, ker_o%nr
            do jj = 1, ker_o%nc
               if(associated(ker_o%blocks(ii, jj)%matrix))then
                  deallocate(ker_o%blocks(ii, jj)%matrix)
                  ker_o%blocks(ii, jj)%matrix => null()
               endif
            enddo
         enddo
         if(allocated(ker_o%blocks))deallocate(ker_o%blocks)
      endif
   end subroutine BF_delete_ker_onelevel


   subroutine BF_copy_ker_onelevel(ker_i, ker_o)

      implicit none
      type(butterfly_kerl)::ker_i, ker_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      character::trans

      if(ker_o%nr>0 .and. ker_o%nc>0 .and. allocated(ker_o%blocks))then
         do ii = 1, ker_o%nr
            do jj = 1, ker_o%nc
               if(associated(ker_o%blocks(ii, jj)%matrix))then
                  deallocate(ker_o%blocks(ii, jj)%matrix)
                  ker_o%blocks(ii, jj)%matrix => null()
               endif
            enddo
         enddo
         if(allocated(ker_o%blocks))deallocate(ker_o%blocks)
      endif

      ker_o%num_row = ker_i%num_row
      ker_o%num_col = ker_i%num_col
      ker_o%nc = ker_i%nc
      ker_o%nr = ker_i%nr
      ker_o%idx_c = ker_i%idx_c
      ker_o%idx_r = ker_i%idx_r
      ker_o%inc_c = ker_i%inc_c
      ker_o%inc_r = ker_i%inc_r
      if(ker_o%nr>0 .and. ker_o%nc>0)then
         allocate (ker_o%blocks(ker_o%nr, ker_o%nc))
         do ii = 1, ker_o%nr
            do jj = 1, ker_o%nc
               nn = size(ker_i%blocks(ii, jj)%matrix, 2)
               rank = size(ker_i%blocks(ii, jj)%matrix, 1)
               allocate (ker_o%blocks(ii, jj)%matrix(rank, nn))
               ker_o%blocks(ii, jj)%matrix = ker_i%blocks(ii, jj)%matrix
            enddo
         enddo
      endif

   end subroutine BF_copy_ker_onelevel


   subroutine BF_copy_delete(block_i, block_o, memory)


      implicit none
      type(matrixblock)::block_i, block_o

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      real(kind=8), optional::memory
      if (present(memory)) memory = 0
      block_o%level = block_i%level
      block_o%col_group = block_i%col_group
      block_o%row_group = block_i%row_group
      block_o%style = block_i%style
      block_o%level_butterfly = block_i%level_butterfly
      block_o%level_half = block_i%level_half
      block_o%rankmax = block_i%rankmax
      block_o%rankmin = block_i%rankmin
      block_o%M = block_i%M
      block_o%N = block_i%N
      block_o%headm = block_i%headm
      block_o%headn = block_i%headn

      block_o%M_loc = block_i%M_loc
      block_o%N_loc = block_i%N_loc
      block_o%pgno = block_i%pgno
      block_o%pgno_db = block_i%pgno_db

      if (associated(block_i%N_p)) then
         allocate (block_o%N_p(size(block_i%N_p, 1), 2))
         block_o%N_p = block_i%N_p
         deallocate (block_i%N_p)
      endif
      if (associated(block_i%M_p)) then
         allocate (block_o%M_p(size(block_i%M_p, 1), 2))
         block_o%M_p = block_i%M_p
         deallocate (block_i%M_p)
      endif

      if (associated(block_i%ms)) then
         allocate (block_o%ms(size(block_i%ms, 1)))
         block_o%ms = block_i%ms
         deallocate (block_i%ms)
      endif
      if (associated(block_i%ns)) then
         allocate (block_o%ns(size(block_i%ns, 1)))
         block_o%ns = block_i%ns
         deallocate (block_i%ns)
      endif


      level_butterfly = block_i%level_butterfly
      num_blocks = 2**level_butterfly

      if (block_i%style == 2) then

         if (level_butterfly /= 0) then
            if (allocated(block_i%ButterflyKerl)) allocate (block_o%ButterflyKerl(level_butterfly))
         end if

         do level = 0, level_butterfly + 1
            if (level == 0) then
               block_o%ButterflyV%num_blk = block_i%ButterflyV%num_blk
               block_o%ButterflyV%nblk_loc = block_i%ButterflyV%nblk_loc
               block_o%ButterflyV%idx = block_i%ButterflyV%idx
               block_o%ButterflyV%inc = block_i%ButterflyV%inc
               if (allocated(block_i%ButterflyV%blocks)) allocate (block_o%ButterflyV%blocks(block_o%ButterflyV%nblk_loc))

               do jj = 1, block_o%ButterflyV%nblk_loc
               if (associated(block_i%ButterflyV%blocks(jj)%matrix)) then
                  nn = size(block_i%ButterflyV%blocks(jj)%matrix, 1)
                  rank = size(block_i%ButterflyV%blocks(jj)%matrix, 2)
                  allocate (block_o%ButterflyV%blocks(jj)%matrix(nn, rank))
                  block_o%ButterflyV%blocks(jj)%matrix = block_i%ButterflyV%blocks(jj)%matrix
                  if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyV%blocks(jj)%matrix)/1024.0d3
                  deallocate (block_i%ButterflyV%blocks(jj)%matrix)
               endif
               enddo

            else if (level == level_butterfly + 1) then
               block_o%ButterflyU%num_blk = block_i%ButterflyU%num_blk
               block_o%ButterflyU%nblk_loc = block_i%ButterflyU%nblk_loc
               block_o%ButterflyU%idx = block_i%ButterflyU%idx
               block_o%ButterflyU%inc = block_i%ButterflyU%inc
               if (allocated(block_i%ButterflyU%blocks)) allocate (block_o%ButterflyU%blocks(block_o%ButterflyU%nblk_loc))

               do ii = 1, block_o%ButterflyU%nblk_loc
                  if (associated(block_i%ButterflyU%blocks(ii)%matrix)) then
                     nn = size(block_i%ButterflyU%blocks(ii)%matrix, 1)
                     rank = size(block_i%ButterflyU%blocks(ii)%matrix, 2)
                     allocate (block_o%ButterflyU%blocks(ii)%matrix(nn, rank))
                     block_o%ButterflyU%blocks(ii)%matrix = block_i%ButterflyU%blocks(ii)%matrix
                     if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyU%blocks(ii)%matrix)/1024.0d3
                     deallocate (block_i%ButterflyU%blocks(ii)%matrix)
                  endif
               enddo
            else
               block_o%ButterflyKerl(level)%num_row = block_i%ButterflyKerl(level)%num_row
               block_o%ButterflyKerl(level)%num_col = block_i%ButterflyKerl(level)%num_col
               block_o%ButterflyKerl(level)%nc = block_i%ButterflyKerl(level)%nc
               block_o%ButterflyKerl(level)%nr = block_i%ButterflyKerl(level)%nr
               block_o%ButterflyKerl(level)%idx_c = block_i%ButterflyKerl(level)%idx_c
               block_o%ButterflyKerl(level)%idx_r = block_i%ButterflyKerl(level)%idx_r
               block_o%ButterflyKerl(level)%inc_c = block_i%ButterflyKerl(level)%inc_c
               block_o%ButterflyKerl(level)%inc_r = block_i%ButterflyKerl(level)%inc_r

               if (allocated(block_i%ButterflyKerl(level)%blocks)) then
                  allocate (block_o%ButterflyKerl(level)%blocks(block_o%ButterflyKerl(level)%nr, block_o%ButterflyKerl(level)%nc))
                  do ii = 1, block_o%ButterflyKerl(level)%nr
                     do jj = 1, block_o%ButterflyKerl(level)%nc
                        if (associated(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix)) then
                           nn = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 2)
                           rank = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 1)
                           allocate (block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix(rank, nn))
                           block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix = block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix
                           if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(ii, jj)%matrix)/1024.0d3

                           ! write(*,*)'jdi',shape(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix),ii,jj,level,block_i%row_group,block_i%col_group
                           deallocate (block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix)
                        endif
                     enddo
                  enddo
                  deallocate (block_i%ButterflyKerl(level)%blocks)
               endif
            endif
         enddo
         if (allocated(block_i%ButterflyU%blocks)) deallocate (block_i%ButterflyU%blocks)
         if (allocated(block_i%ButterflyV%blocks)) deallocate (block_i%ButterflyV%blocks)
         if (allocated(block_i%ButterflyKerl)) then
            ! if(level_butterfly/=0)deallocate(block_i%ButterflyKerl)
            deallocate (block_i%ButterflyKerl)
         endif
         if (allocated(block_i%ButterflyMiddle)) then
            write (*, *) 'ButterflyMiddle is not yet handled when the butterfly is distributed'
            levelm = ceiling_safe(dble(level_butterfly)/2d0)
            allocate (block_o%ButterflyMiddle(2**levelm, 2**(level_butterfly - levelm)))

            do index_i_m = 1, 2**levelm
               do index_j_m = 1, 2**(level_butterfly - levelm)
                  rank = size(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, 1)
                  allocate (block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix(rank, rank))
                  block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix = block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix
                  deallocate (block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix)
               end do
            end do
            deallocate (block_i%ButterflyMiddle)
         end if

      else if (block_i%style == 1) then
         mm = size(block_i%fullmat, 1)
         nn = size(block_i%fullmat, 2)
         if(associated(block_i%fullmat))then
            allocate (block_o%fullmat(mm, nn))
            block_o%fullmat = block_i%fullmat
            deallocate (block_i%fullmat)
         endif
#if HAVE_ZFP
         if(allocated(block_i%FullmatZFP%buffer_r))then
            mm = size(block_i%FullmatZFP%buffer_r,1)
            allocate(block_o%FullmatZFP%buffer_r(mm))
            block_o%FullmatZFP%buffer_r = block_i%FullmatZFP%buffer_r
            block_o%FullmatZFP%stream_r = block_i%FullmatZFP%stream_r
            deallocate(block_i%FullmatZFP%buffer_r)
            call zFORp_stream_close(block_i%FullmatZFP%stream_r)
         endif
         if(allocated(block_i%FullmatZFP%buffer_i))then
            mm = size(block_i%FullmatZFP%buffer_i,1)
            allocate(block_o%FullmatZFP%buffer_i(mm))
            block_o%FullmatZFP%buffer_i = block_i%FullmatZFP%buffer_i
            block_o%FullmatZFP%stream_i = block_i%FullmatZFP%stream_i
            deallocate(block_i%FullmatZFP%buffer_i)
            call zFORp_stream_close(block_i%FullmatZFP%stream_i)
         endif
#endif
      else
         write (*, *) 'block style not implemented'
         stop
      end if

   end subroutine BF_copy_delete

   subroutine BF_ComputeMemory(block_i, memory)


      implicit none
      type(matrixblock)::block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      real(kind=8)::memory
      memory = 0

      level_butterfly = block_i%level_butterfly
      num_blocks = 2**level_butterfly

      if (block_i%style == 2) then
         if (block_i%M_loc > 0) then
            do level = 0, level_butterfly + 1
               if (level == 0) then
                  do jj = 1, block_i%ButterflyV%nblk_loc
                     memory = memory + SIZEOF(block_i%ButterflyV%blocks(jj)%matrix)/1024.0d3
                  enddo
               elseif (level == level_butterfly + 1) then
                  do jj = 1, block_i%ButterflyU%nblk_loc
                     memory = memory + SIZEOF(block_i%ButterflyU%blocks(jj)%matrix)/1024.0d3
                  enddo
               else
                  do ii = 1, block_i%ButterflyKerl(level)%nr
                     do jj = 1, block_i%ButterflyKerl(level)%nc
                        memory = memory + SIZEOF(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix)/1024.0d3
                     enddo
                  enddo
               endif
            enddo
         endif

      else if (block_i%style == 1) then
         if(associated(block_i%fullmat))memory = memory + SIZEOF(block_i%fullmat)/1024.0d3
#if HAVE_ZFP
         if(allocated(block_i%FullmatZFP%buffer_r))memory = memory + SIZEOF(block_i%FullmatZFP%buffer_r)/1024.0d3
         if(allocated(block_i%FullmatZFP%buffer_i))memory = memory + SIZEOF(block_i%FullmatZFP%buffer_i)/1024.0d3
#endif
      else
         write (*, *) 'block style not implemented'
         stop
      end if

   end subroutine BF_ComputeMemory




   subroutine BF_MD_ComputeMemory(Ndim, blocks, memory_dense,memory_comp)


      implicit none
      integer Ndim,dim_i,dim, dim_MD(Ndim+2),dim_subtensor(Ndim*2),nc(Ndim),idx_MD(Ndim+2),index_r_vector(Ndim),index_r_scalar,index_c_vector(Ndim),index_c_scalar
      type(matrixblock_MD)::blocks

      integer bb, i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, edge, patch, node, group, level_c
      integer::level_butterfly, level_half, level_final

      real(kind=8)::memory_dense,memory_comp
      memory_dense = 0
      memory_comp = 0

      level_butterfly = blocks%level_butterfly
      level_half = blocks%level_half

      if (blocks%style == 2) then
         if (allocated(blocks%M_loc) .and. allocated(blocks%N_loc)) then
         if (product(blocks%M_loc)/=0 .and. product(blocks%N_loc)/=0) then

            do bb=1, product(blocks%nc_m)
            do level = 0, level_half
               if (level == 0) then
                  do jj = 1, blocks%ButterflyV(bb)%nblk_loc
                     do dim_i=1,Ndim
                        memory_comp = memory_comp + SIZEOF(blocks%ButterflyV(bb)%blocks(jj,dim_i)%matrix)/1024.0d3
                     enddo
                  enddo
               elseif (level == level_butterfly + 1) then
                  write(*,*)"should not reach level == level_butterfly + 1 for the right half of matrixblock_MD"
               else
                  dim_MD(1:Ndim)=blocks%ButterflyKerl_R(bb,level)%nr
                  dim_MD(1+Ndim)=blocks%ButterflyKerl_R(bb,level)%nc(1)
                  dim_MD(2+Ndim)=Ndim
                  do index_ij=1,product(dim_MD)
                     call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
                     index_r_vector = idx_MD(1:Ndim)
                     dim = idx_MD(Ndim+2)
                     index_j = idx_MD(Ndim+1)
                     call MultiIndexToSingleIndex(Ndim,dim_MD(1:Ndim), index_r_scalar, index_r_vector)
                     memory_comp = memory_comp + SIZEOF(blocks%ButterflyKerl_R(bb,level)%blocks(index_r_scalar,index_j,dim)%matrix)/1024.0d3
                  enddo
               endif
            enddo
            enddo

            level_final = level_half + 1
            do bb=1, product(blocks%nr_m)
            do level = level_butterfly + 1, level_final, -1
               if (level == 0) then
                  write(*,*)"should not reach level == 0 for the left half of matrixblock_MD"
               elseif (level == level_butterfly + 1) then
                  do ii = 1, blocks%ButterflyU(bb)%nblk_loc
                     do dim_i=1,Ndim
                        memory_comp = memory_comp + SIZEOF(blocks%ButterflyU(bb)%blocks(ii,dim_i)%matrix)/1024.0d3
                     enddo
                  enddo
               else
                  dim_MD(1)=blocks%ButterflyKerl_L(bb,level)%nr(1)
                  dim_MD(2:1+Ndim)=blocks%ButterflyKerl_L(bb,level)%nc
                  dim_MD(2+Ndim)=Ndim
                  do index_ij=1,product(dim_MD)
                     call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
                     index_i = idx_MD(1)
                     dim = idx_MD(Ndim+2)
                     index_c_vector = idx_MD(2:Ndim+1)
                     call MultiIndexToSingleIndex(Ndim,dim_MD(2:1+Ndim), index_c_scalar, index_c_vector)
                     memory_comp = memory_comp + SIZEOF(blocks%ButterflyKerl_L(bb,level)%blocks(index_i,index_c_scalar,dim)%matrix)/1024.0d3
                  enddo
               endif
            enddo
            enddo

            nc=2**(level_butterfly -level_half)
            dim_subtensor(1:Ndim)=blocks%nr_m
            dim_subtensor(1+Ndim:Ndim*2)=nc
            do index_ij=1,product(dim_subtensor)
               if(associated(blocks%ButterflyMiddle(index_ij)%matrix))memory_dense = memory_dense + SIZEOF(blocks%ButterflyMiddle(index_ij)%matrix)/1024.0d3
#if HAVE_ZFP
               if(allocated(blocks%MiddleZFP))then
                  if(allocated(blocks%MiddleZFP(index_ij)%buffer_r))memory_dense = memory_dense + SIZEOF(blocks%MiddleZFP(index_ij)%buffer_r)/1024.0d3
                  if(allocated(blocks%MiddleZFP(index_ij)%buffer_i))memory_dense = memory_dense + SIZEOF(blocks%MiddleZFP(index_ij)%buffer_i)/1024.0d3
               endif
#endif
               if(allocated(blocks%MiddleQTT))then
                  if(allocated(blocks%MiddleQTT(index_ij)%core))memory_dense = memory_dense + SIZEOF(blocks%MiddleQTT(index_ij)%core)/1024.0d3
                  if(allocated(blocks%MiddleQTT(index_ij)%coreZFP%buffer_r))memory_dense = memory_dense + SIZEOF(blocks%MiddleQTT(index_ij)%coreZFP%buffer_r)/1024.0d3
                  if(allocated(blocks%MiddleQTT(index_ij)%coreZFP%buffer_i))memory_dense = memory_dense + SIZEOF(blocks%MiddleQTT(index_ij)%coreZFP%buffer_i)/1024.0d3
               endif
            enddo
         endif
         endif

      else if (blocks%style == 1) then
         if(associated(blocks%fullmat))memory_dense = memory_dense + SIZEOF(blocks%fullmat)/1024.0d3
#if HAVE_ZFP
         if(allocated(blocks%FullmatZFP%buffer_r))memory_dense = memory_dense + SIZEOF(blocks%FullmatZFP%buffer_r)/1024.0d3
         if(allocated(blocks%FullmatZFP%buffer_i))memory_dense = memory_dense + SIZEOF(blocks%FullmatZFP%buffer_i)/1024.0d3
#endif
         if(allocated(blocks%FullmatQTT%core))memory_dense = memory_dense + SIZEOF(blocks%FullmatQTT%core)/1024.0d3
         if(allocated(blocks%FullmatQTT%coreZFP%buffer_r))memory_dense = memory_dense + SIZEOF(blocks%FullmatQTT%coreZFP%buffer_r)/1024.0d3
         if(allocated(blocks%FullmatQTT%coreZFP%buffer_i))memory_dense = memory_dense + SIZEOF(blocks%FullmatQTT%coreZFP%buffer_i)/1024.0d3
      else
         write (*, *) 'block style not implemented'
         stop
      end if

   end subroutine BF_MD_ComputeMemory





   subroutine BF_MD_delete(Ndim, blocks, allflag)


      implicit none

      integer Ndim
      integer butterflyB_inuse, num_col, num_row, bb, index_i_k, index_j_k, dim_i, nc(Ndim), dim_subtensor(Ndim*2)
      integer index_ij, i, j, mm, nn, rank, num_blocks, level, level_butterfly, index_i_m, index_j_m, levelm
      real(kind=8) memory_butterfly, rtemp
      type(matrixblock_MD)::blocks
      integer allflag

      level_butterfly = blocks%level_butterfly


      if(allocated(blocks%nr_m))then

         if(allocated(blocks%ButterflyU))then
            do bb=1, product(blocks%nr_m)
               if(allocated(blocks%ButterflyU(bb)%blocks))then
                  do index_i_k=1,blocks%ButterflyU(bb)%num_blk
                     do dim_i=1,Ndim
                        if(associated(blocks%ButterflyU(bb)%blocks(index_i_k,dim_i)%matrix))then
                           deallocate(blocks%ButterflyU(bb)%blocks(index_i_k,dim_i)%matrix)
                        endif
                     enddo
                  enddo
                  deallocate(blocks%ButterflyU(bb)%blocks)
               endif
            enddo
            deallocate(blocks%ButterflyU)
         endif

         if(allocated(blocks%ButterflyKerl_L))then
            do bb=1, product(blocks%nr_m)
               do level=blocks%level_half+1,level_butterfly
                  if(allocated(blocks%ButterflyKerl_L(bb,level)%blocks))then
                     do index_i_k=1,blocks%ButterflyKerl_L(bb,level)%nr(1)
                        do index_j_k=1,product(blocks%ButterflyKerl_L(bb,level)%nc)
                           do dim_i=1,Ndim
                              if(associated(blocks%ButterflyKerl_L(bb,level)%blocks(index_i_k,index_j_k,dim_i)%matrix))then
                                 deallocate(blocks%ButterflyKerl_L(bb,level)%blocks(index_i_k,index_j_k,dim_i)%matrix)
                              endif
                           enddo
                        enddo
                     enddo
                     deallocate(blocks%ButterflyKerl_L(bb,level)%num_row)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%num_col)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%nr)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%nc)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%idx_r)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%idx_c)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%inc_r)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%inc_c)
                     deallocate(blocks%ButterflyKerl_L(bb,level)%blocks)
                  endif
               enddo
            enddo
            deallocate(blocks%ButterflyKerl_L)
         endif

         nc=2**(level_butterfly -blocks%level_half)
         dim_subtensor(1:Ndim)=blocks%nr_m
         dim_subtensor(1+Ndim:Ndim*2)=nc
         if(allocated(blocks%ButterflyMiddle))then
            do index_ij = 1, product(dim_subtensor)
               if(associated(blocks%ButterflyMiddle(index_ij)%matrix))then
                  deallocate(blocks%ButterflyMiddle(index_ij)%matrix)
               endif
               if(allocated(blocks%ButterflyMiddle(index_ij)%dims_m))then
                  deallocate(blocks%ButterflyMiddle(index_ij)%dims_m)
               endif
               if(allocated(blocks%ButterflyMiddle(index_ij)%dims_n))then
                  deallocate(blocks%ButterflyMiddle(index_ij)%dims_n)
               endif
            enddo
            deallocate(blocks%ButterflyMiddle)
         endif
      endif

      if(allocated(blocks%nc_m))then

         if(allocated(blocks%ButterflyV))then
            do bb=1, product(blocks%nc_m)
               if(allocated(blocks%ButterflyV(bb)%blocks))then
                  do index_i_k=1,blocks%ButterflyV(bb)%num_blk
                     do dim_i=1,Ndim
                        if(associated(blocks%ButterflyV(bb)%blocks(index_i_k,dim_i)%matrix))then
                           deallocate(blocks%ButterflyV(bb)%blocks(index_i_k,dim_i)%matrix)
                        endif
                     enddo
                  enddo
                  deallocate(blocks%ButterflyV(bb)%blocks)
               endif
            enddo
            deallocate(blocks%ButterflyV)
         endif

         if(allocated(blocks%ButterflyKerl_R))then
            do bb=1, product(blocks%nc_m)
               do level=1,blocks%level_half
                  if(allocated(blocks%ButterflyKerl_R(bb,level)%blocks))then
                     do index_j_k=1,blocks%ButterflyKerl_R(bb,level)%nc(1)
                        do index_i_k=1,product(blocks%ButterflyKerl_R(bb,level)%nr)
                           do dim_i=1,Ndim
                              if(associated(blocks%ButterflyKerl_R(bb,level)%blocks(index_i_k,index_j_k,dim_i)%matrix))then
                                 deallocate(blocks%ButterflyKerl_R(bb,level)%blocks(index_i_k,index_j_k,dim_i)%matrix)
                              endif
                           enddo
                        enddo
                     enddo
                     deallocate(blocks%ButterflyKerl_R(bb,level)%num_row)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%num_col)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%nr)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%nc)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%idx_r)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%idx_c)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%inc_r)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%inc_c)
                     deallocate(blocks%ButterflyKerl_R(bb,level)%blocks)
                  endif
               enddo
            enddo
            deallocate(blocks%ButterflyKerl_R)
         endif
      endif

      blocks%rankmax = -1000
      blocks%rankmin = 1000

      if (associated(blocks%fullmat)) deallocate (blocks%fullmat)
#if HAVE_ZFP
      if (allocated(blocks%FullmatZFP%buffer_r)) then
         deallocate (blocks%FullmatZFP%buffer_r)
         call zFORp_stream_close(blocks%FullmatZFP%stream_r)
      endif
      if (allocated(blocks%FullmatZFP%buffer_i))then
         deallocate (blocks%FullmatZFP%buffer_i)
         call zFORp_stream_close(blocks%FullmatZFP%stream_i)
      endif
      if(allocated(blocks%MiddleZFP))then
         do index_ij = 1, size(blocks%MiddleZFP,1)
            if (allocated(blocks%MiddleZFP(index_ij)%buffer_r)) then
               deallocate (blocks%MiddleZFP(index_ij)%buffer_r)
               call zFORp_stream_close(blocks%MiddleZFP(index_ij)%stream_r)
            endif
            if (allocated(blocks%MiddleZFP(index_ij)%buffer_i))then
               deallocate (blocks%MiddleZFP(index_ij)%buffer_i)
               call zFORp_stream_close(blocks%MiddleZFP(index_ij)%stream_i)
            endif
         enddo
         deallocate(blocks%MiddleZFP)
      endif
#endif
      if(allocated(blocks%MiddleQTT))then
         call TT_Delete(blocks%MiddleQTT(1))
         deallocate(blocks%MiddleQTT)
      endif
      if (allocated(blocks%FullmatQTT%core)) then
         call TT_Delete(blocks%FullmatQTT)
      endif
      if (allocated(blocks%fullmat_MPI)) deallocate (blocks%fullmat_MPI)
      if (allocated(blocks%ipiv)) deallocate (blocks%ipiv)
      if (allocated(blocks%Butterfly_data_MPI)) deallocate (blocks%Butterfly_data_MPI)
      if (allocated(blocks%Butterfly_index_MPI)) deallocate (blocks%Butterfly_index_MPI)
      if (associated(blocks%ms)) deallocate (blocks%ms)
      if (associated(blocks%ns)) deallocate (blocks%ns)


      if (allflag == 1) then
         if (associated(blocks%N_p)) deallocate (blocks%N_p)
         if (associated(blocks%M_p)) deallocate (blocks%M_p)
         if (allocated(blocks%row_group)) deallocate(blocks%row_group)
         if (allocated(blocks%col_group)) deallocate(blocks%col_group)
         if (allocated(blocks%M)) deallocate(blocks%M)
         if (allocated(blocks%N)) deallocate(blocks%N)
         if (allocated(blocks%M_loc)) deallocate(blocks%M_loc)
         if (allocated(blocks%N_loc)) deallocate(blocks%N_loc)
         if (allocated(blocks%headm)) deallocate(blocks%headm)
         if (allocated(blocks%headn)) deallocate(blocks%headn)
         if (allocated(blocks%nr_m)) deallocate(blocks%nr_m)
         if (allocated(blocks%nc_m)) deallocate(blocks%nc_m)
         if (allocated(blocks%idx_r_m)) deallocate(blocks%idx_r_m)
         if (allocated(blocks%idx_c_m)) deallocate(blocks%idx_c_m)
      endif
      return

   end subroutine BF_MD_delete





   integer function BF_Switchlevel(level_butterfly, pat_comp)


      implicit none
      integer::pat_comp
      integer::level_butterfly

      if (pat_comp == 1) BF_Switchlevel = level_butterfly ! from right to left until the second last level
      if (pat_comp == 2) BF_Switchlevel = 0 ! from left to right until the second last level
      if (pat_comp == 3) BF_Switchlevel = floor_safe(dble(level_butterfly)/2d0)  ! from outer to inner

   end function BF_Switchlevel

   logical function BF_checkNAN(block_i)


      implicit none
      type(matrixblock)::block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      DTR:: temp

#ifdef NDEBUG
      BF_checkNAN=.false.
#else
      level_butterfly = block_i%level_butterfly
      num_blocks = 2**level_butterfly
      temp = 0

      if (block_i%style == 2) then
         do level = 0, level_butterfly + 1
            if (level == 0) then
               if (allocated(block_i%ButterflyV%blocks)) then
               do index_j = 1, block_i%ButterflyV%nblk_loc
                  if (associated(block_i%ButterflyV%blocks(index_j)%matrix)) then
                     mm = size(block_i%ButterflyV%blocks(index_j)%matrix, 1)
                     nn = size(block_i%ButterflyV%blocks(index_j)%matrix, 2)
                     temp = temp + fnorm(block_i%ButterflyV%blocks(index_j)%matrix, mm, nn)
                     if (myisnan(temp)) write (*, *) 'V', level_butterfly, index_j, fnorm(block_i%ButterflyV%blocks(index_j)%matrix, mm, nn), mm, nn
                  endif
               enddo
               endif
            elseif (level == level_butterfly + 1) then
               if (allocated(block_i%ButterflyU%blocks)) then
               do index_i = 1, block_i%ButterflyU%nblk_loc
                  if (associated(block_i%ButterflyU%blocks(index_i)%matrix)) then
                     mm = size(block_i%ButterflyU%blocks(index_i)%matrix, 1)
                     nn = size(block_i%ButterflyU%blocks(index_i)%matrix, 2)
                     temp = temp + fnorm(block_i%ButterflyU%blocks(index_i)%matrix, mm, nn)
                     if (myisnan(temp)) write (*, *) 'U', level_butterfly, index_i, fnorm(block_i%ButterflyU%blocks(index_i)%matrix, mm, nn), mm, nn
                  endif
               enddo
               endif
            else
               if (allocated(block_i%ButterflyKerl)) then
               if (allocated(block_i%ButterflyKerl(level)%blocks)) then
               do index_i = 1, block_i%ButterflyKerl(level)%nr
                  do index_j = 1, block_i%ButterflyKerl(level)%nc
                     if (associated(block_i%ButterflyKerl(level)%blocks(index_i, index_j)%matrix)) then
                        mm = size(block_i%ButterflyKerl(level)%blocks(index_i, index_j)%matrix, 1)
                        nn = size(block_i%ButterflyKerl(level)%blocks(index_i, index_j)%matrix, 2)
                        temp = temp + fnorm(block_i%ButterflyKerl(level)%blocks(index_i, index_j)%matrix, mm, nn)
                        if (myisnan(temp)) write (*, *) 'Ker', level_butterfly, level, index_i, index_j, fnorm(block_i%ButterflyKerl(level)%blocks(index_i, index_j)%matrix, mm, nn), mm, nn
                     endif
                  enddo
               enddo
               endif
               endif
            endif
         enddo
      else if (block_i%style == 1) then
         if (associated(block_i%fullmat)) then
            mm = size(block_i%fullmat, 1)
            nn = size(block_i%fullmat, 2)
            temp = temp + fnorm(block_i%fullmat, mm, nn)
         endif
      else
         write (*, *) 'block style not implemented'
         stop
      end if

      BF_checkNAN = myisnan(temp)
#endif
   end function BF_checkNAN




   subroutine BF_print_size(block_i)


      implicit none
      type(matrixblock)::block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, index_i_loc, index_j_loc, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c, inc_r, inc_c, idx_r, idx_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      DTR:: temp

      level_butterfly = block_i%level_butterfly
      num_blocks = 2**level_butterfly

      if (block_i%style == 2) then
         do level = 0, level_butterfly + 1
            if (level == 0) then
               if (allocated(block_i%ButterflyV%blocks)) then
               do index_j_loc = 1, block_i%ButterflyV%nblk_loc
                  if (associated(block_i%ButterflyV%blocks(index_j_loc)%matrix)) then
                     idx_c = block_i%ButterflyV%idx
                     inc_c = block_i%ButterflyV%inc
                     index_j = (index_j_loc - 1)*inc_c + idx_c
                     mm = size(block_i%ButterflyV%blocks(index_j_loc)%matrix, 1)
                     nn = size(block_i%ButterflyV%blocks(index_j_loc)%matrix, 2)
                     write (*, *) "BF",block_i%row_group,block_i%col_group,'L:',level_butterfly, 'V idx:', index_j, 'size:', mm, nn
                  endif
               enddo
               endif
            elseif (level == level_butterfly + 1) then
               if (allocated(block_i%ButterflyU%blocks)) then
               do index_i_loc = 1, block_i%ButterflyU%nblk_loc
                  if (associated(block_i%ButterflyU%blocks(index_i_loc)%matrix)) then
                     idx_r = block_i%ButterflyU%idx
                     inc_r = block_i%ButterflyU%inc
                     index_i = (index_i_loc - 1)*inc_r + idx_r

                     mm = size(block_i%ButterflyU%blocks(index_i_loc)%matrix, 1)
                     nn = size(block_i%ButterflyU%blocks(index_i_loc)%matrix, 2)
                     write (*, *) "BF",block_i%row_group,block_i%col_group,'L:',level_butterfly, 'U idx:', index_i, 'size:', mm, nn
                  endif
               enddo
               endif
            else
               if (allocated(block_i%ButterflyKerl)) then
               if (allocated(block_i%ButterflyKerl(level)%blocks)) then
               do index_i_loc = 1, block_i%ButterflyKerl(level)%nr
                  do index_j_loc = 1, block_i%ButterflyKerl(level)%nc

                     inc_r = block_i%ButterflyKerl(level)%inc_r
                     idx_r = block_i%ButterflyKerl(level)%idx_r
                     inc_c = block_i%ButterflyKerl(level)%inc_c
                     idx_c = block_i%ButterflyKerl(level)%idx_c
                     index_i = (index_i_loc - 1)*inc_r + idx_r
                     index_j = (index_j_loc - 1)*inc_c + idx_c

                     if (associated(block_i%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix)) then
                        mm = size(block_i%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix, 1)
                        nn = size(block_i%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix, 2)
                        write (*, *) "BF",block_i%row_group,block_i%col_group,'L:', level_butterfly, "ker level:", level,' idx:', index_i, index_j, 'size:', mm, nn
                     endif
                  enddo
               enddo
               endif
               endif
            endif
         enddo
      else if (block_i%style == 1) then
         if (associated(block_i%fullmat)) then
            mm = size(block_i%fullmat, 1)
            nn = size(block_i%fullmat, 2)
            write (*, *) "Dense",block_i%row_group,block_i%col_group, 'size:', mm, nn
         endif
      else
         write (*, *) 'block style not implemented'
         stop
      end if

   end subroutine BF_print_size




   subroutine BF_print_size_rank(block_i, tolerance)


      implicit none
      type(matrixblock)::block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, truerank, index_i, index_j, levelm, index_i_m, index_j_m, mm1, mm2, nn1, nn2
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      DT, allocatable::matrixtemp(:, :), mat11(:, :), mat12(:, :), mat21(:, :), mat22(:, :)
      real(kind=8)::tolerance

      level_butterfly = block_i%level_butterfly
      num_blocks = 2**level_butterfly

      do level = 0, level_butterfly + 1
         ! write(*,*)level
         if (level == 0) then
            do index_ij = 1, 2**level_butterfly
               nn = size(block_i%ButterflyV%blocks(index_ij)%matrix, 1)
               rank = size(block_i%ButterflyV%blocks(index_ij)%matrix, 2)
               allocate (matrixtemp(nn, rank))
               matrixtemp = block_i%ButterflyV%blocks(index_ij)%matrix
               call GetRank(nn, rank, matrixtemp, truerank, tolerance)
               write (*, *) level, index_ij, nn, rank, truerank
               deallocate (matrixtemp)
            end do
         else if (level == level_butterfly + 1) then
            do index_ij = 1, 2**level_butterfly
               mm = size(block_i%ButterflyU%blocks(index_ij)%matrix, 1)
               rank = size(block_i%ButterflyU%blocks(index_ij)%matrix, 2)
               allocate (matrixtemp(mm, rank))
               matrixtemp = block_i%ButterflyU%blocks(index_ij)%matrix
               call GetRank(mm, rank, matrixtemp, truerank, tolerance)
               write (*, *) level, index_ij, mm, rank, truerank
               deallocate (matrixtemp)
            end do
         else
            do index_i = 1, 2**(level - 1)
               do index_j = 1, 2**(level_butterfly - level)

                  mm1 = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, 2*index_j - 1)%matrix, 1)
                  nn1 = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, 2*index_j - 1)%matrix, 2)

                  mm2 = size(block_i%ButterflyKerl(level)%blocks(2*index_i, 2*index_j)%matrix, 1)
                  nn2 = size(block_i%ButterflyKerl(level)%blocks(2*index_i, 2*index_j)%matrix, 2)

                  allocate (mat11(mm1, nn1))
                  mat11 = block_i%ButterflyKerl(level)%blocks(2*index_i - 1, 2*index_j - 1)%matrix
                  allocate (mat12(mm1, nn2))
                  mat12 = block_i%ButterflyKerl(level)%blocks(2*index_i - 1, 2*index_j)%matrix
                  allocate (mat21(mm2, nn1))
                  mat21 = block_i%ButterflyKerl(level)%blocks(2*index_i, 2*index_j - 1)%matrix
                  allocate (mat22(mm2, nn2))
                  mat22 = block_i%ButterflyKerl(level)%blocks(2*index_i, 2*index_j)%matrix
                  allocate (matrixtemp(mm1 + mm2, nn1 + nn2))
                  matrixtemp(1:mm1, 1:nn1) = mat11
                  matrixtemp(1:mm1, 1 + nn1:nn2 + nn1) = mat12
                  matrixtemp(1 + mm1:mm2 + mm1, 1:nn1) = mat21
                  matrixtemp(1 + mm1:mm2 + mm1, 1 + nn1:nn2 + nn1) = mat22
                  call GetRank(mm1 + mm2, nn1 + nn2, matrixtemp, truerank, tolerance)
                  write (*, *) level, index_i, index_j, (mm1 + mm2), (nn1 + nn2), truerank
                  deallocate (mat11)
                  deallocate (mat12)
                  deallocate (mat21)
                  deallocate (mat22)
                  deallocate (matrixtemp)

               enddo
            enddo

         end if

      enddo

   end subroutine BF_print_size_rank

   subroutine BF_extract_partial(block_o, level_butterfly_loc, ij_loc, head, group,LR, agent_block,pgno,ptree)


      implicit none

      type(proctree)::ptree
      type(matrixblock)::block_o, agent_block
      integer level_butterfly, level_butterfly_loc, ij_loc, index_i, index_i_start, index_j_start, index_j, level, i, ii,jj, nn, mm, num_blocks, rank, pgno, index_i_loc, index_j_loc
      character LR
      integer idx_r, inc_r, nr, idx_c, inc_c, nc, nblk_loc, idx, inc, pp, ierr,  head, group
      type(butterfly_skel)::sizes
      real(kind=8)::n1,n2
      integer reqm, reqn
      integer statusm(MPI_status_size), statusn(MPI_status_size)
      ! call MPI_Barrier(ptree%pgrp(pgno)%Comm, ierr)
      ! allocate(agent_block)
      n1=MPI_Wtime()
      call assert(level_butterfly_loc >= 1, 'level_butterfly_loc cannot be zero')

      agent_block%style = block_o%style
      agent_block%level_butterfly = level_butterfly_loc
      agent_block%pgno = pgno
      agent_block%pgno_db = pgno
      agent_block%level = -1
      agent_block%M_loc = 0
      agent_block%N_loc = 0

      level_butterfly = block_o%level_butterfly
      num_blocks = 2**level_butterfly_loc

      allocate (agent_block%ButterflyKerl(level_butterfly_loc))

      if (LR == 'L') then
         agent_block%headm = head
         agent_block%headn = 0
         agent_block%row_group = group
         agent_block%col_group = -1
         agent_block%level_half=BF_Switchlevel(agent_block%level_butterfly, 2)
         do level = 0, level_butterfly_loc+1
            if(level==0)then
               index_i_start = (ij_loc - 1)*2

               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
               sizes%idx_r = idx_r
               sizes%idx_c = idx_c
               sizes%inc_r = inc_r
               sizes%inc_c = inc_c
               sizes%nr = nr
               sizes%nc = nc
               idx = idx_c
               inc = inc_c
               nblk_loc = nc
               allocate (sizes%inds(sizes%nr, sizes%nc))

               do jj = 1, nblk_loc
                  index_i = 1
                  index_j = (jj - 1)*inc + idx
                  index_i_loc = (index_i+index_i_start - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + 1)%idx_r)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + 1)%inc_r + 1 !index_i_loc is local index in block_o
                  index_j_loc = (index_j - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + 1)%idx_c)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + 1)%inc_c + 1
                  nn = size(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + 1)%blocks(index_i_loc,index_j_loc)%matrix, 2)
                  sizes%inds(1,jj)%size = nn
               enddo
               call BF_all2all_sizes(agent_block, sizes, ptree, level, 'C', 'R')


               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
               idx = idx_c
               inc = inc_c
               nblk_loc = nc
               agent_block%ButterflyV%idx = idx
               agent_block%ButterflyV%inc = inc
               agent_block%ButterflyV%nblk_loc = nblk_loc
               agent_block%ButterflyV%num_blk = 2**level_butterfly_loc
               allocate (agent_block%ButterflyV%blocks(nblk_loc))
               do jj = 1, nblk_loc
                  nn = sizes%inds(1,jj)%size
                  allocate (agent_block%ButterflyV%blocks(jj)%matrix(nn, nn))
                  agent_block%ButterflyV%blocks(jj)%matrix = 0
                  do ii = 1, nn
                     agent_block%ButterflyV%blocks(jj)%matrix(ii, ii) = 1
                  end do
                  agent_block%N_loc =agent_block%N_loc + nn
               enddo
               deallocate(sizes%inds)
            elseif(level==level_butterfly_loc+1)then
               index_i_start = (ij_loc - 1)*2**level_butterfly_loc
               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
               idx = idx_r
               inc = inc_r
               nblk_loc = nr
               agent_block%ButterflyU%idx = idx
               agent_block%ButterflyU%inc = inc
               agent_block%ButterflyU%nblk_loc = nblk_loc
               agent_block%ButterflyU%num_blk = 2**level_butterfly_loc

               allocate (agent_block%ButterflyU%blocks(nblk_loc))

               do ii = 1, nblk_loc
                  index_i = (ii - 1)*inc + idx
                  index_i_loc = (index_i+index_i_start - block_o%ButterflyU%idx)/block_o%ButterflyU%inc + 1 !index_i_loc is local index in block_o
                  mm = size(block_o%ButterflyU%blocks(index_i_loc)%matrix, 1)
                  nn = size(block_o%ButterflyU%blocks(index_i_loc)%matrix, 2)
                  allocate (agent_block%ButterflyU%blocks(ii)%matrix(mm, nn))
                  agent_block%ButterflyU%blocks(ii)%matrix = block_o%ButterflyU%blocks(index_i_loc)%matrix
                  agent_block%M_loc =agent_block%M_loc + mm
               enddo
            else
               index_i_start = (ij_loc - 1)*2**level
               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
               idx_r = idx_r*2 - 1
               inc_r = inc_r
               nr = nr*2
               if(nr>0 .and. nc>0)then
                  allocate (agent_block%ButterflyKerl(level)%blocks(nr, nc))

                  agent_block%ButterflyKerl(level)%num_row = 2**level
                  agent_block%ButterflyKerl(level)%num_col = 2**(level_butterfly_loc - level + 1)
                  agent_block%ButterflyKerl(level)%idx_r = idx_r
                  agent_block%ButterflyKerl(level)%idx_c = idx_c
                  agent_block%ButterflyKerl(level)%inc_r = inc_r
                  agent_block%ButterflyKerl(level)%inc_c = inc_c
                  agent_block%ButterflyKerl(level)%nr = nr
                  agent_block%ButterflyKerl(level)%nc = nc

                  do ii = 1, nr
                     do jj = 1, nc
                        index_i = (ii - 1)*inc_r + idx_r
                        index_j = (jj - 1)*inc_c + idx_c

                        index_i_loc = (index_i+index_i_start - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%idx_r)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%inc_r + 1 !index_i_loc is local index in block_o
                        index_j_loc = (index_j - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%idx_c)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%inc_c + 1
                        nn = size(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix, 2)
                        mm = size(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix, 1)
                        allocate (agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix(mm, nn))
                        agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix = block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix
                     enddo
                  enddo
               endif
            endif
         enddo
      else if (LR == 'R') then
         agent_block%headm = 0
         agent_block%headn = head
         agent_block%row_group = -1
         agent_block%col_group = group
         agent_block%level_half=BF_Switchlevel(agent_block%level_butterfly, 1)
         do level = 0, level_butterfly_loc+1
            if(level==0)then
               index_j_start = (ij_loc - 1)*2**(level_butterfly_loc)

               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
               idx = idx_c
               inc = inc_c
               nblk_loc = nc
               agent_block%ButterflyV%idx = idx
               agent_block%ButterflyV%inc = inc
               agent_block%ButterflyV%nblk_loc = nblk_loc
               agent_block%ButterflyV%num_blk = 2**level_butterfly_loc
               allocate (agent_block%ButterflyV%blocks(nblk_loc))

               do jj = 1, nblk_loc
                  index_j = (jj - 1)*inc + idx
                  index_j_loc = (index_j+index_j_start - block_o%ButterflyV%idx)/block_o%ButterflyV%inc + 1 !index_i_loc is local index in block_o
                  mm = size(block_o%ButterflyV%blocks(index_j_loc)%matrix, 1)
                  nn = size(block_o%ButterflyV%blocks(index_j_loc)%matrix, 2)
                  allocate (agent_block%ButterflyV%blocks(jj)%matrix(mm, nn))
                  agent_block%ButterflyV%blocks(jj)%matrix = block_o%ButterflyV%blocks(index_j_loc)%matrix
                  agent_block%N_loc =agent_block%N_loc + mm
               enddo

            elseif(level==level_butterfly_loc+1)then
               index_j_start = (ij_loc - 1)*2

               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
               sizes%idx_r = idx_r
               sizes%idx_c = idx_c
               sizes%inc_r = inc_r
               sizes%inc_c = inc_c
               sizes%nr = nr
               sizes%nc = nc
               idx = idx_r
               inc = inc_r
               nblk_loc = nr
               allocate (sizes%inds(sizes%nr, sizes%nc))
               do ii = 1, nblk_loc
                  index_i = (ii - 1)*inc + idx
                  index_j = 1
                  index_i_loc = (index_i - block_o%ButterflyKerl(level-1)%idx_r)/block_o%ButterflyKerl(level-1)%inc_r + 1 !index_i_loc is local index in block_o
                  index_j_loc = (index_j + index_j_start - block_o%ButterflyKerl(level-1)%idx_c)/block_o%ButterflyKerl(level-1)%inc_c + 1
                  mm = size(block_o%ButterflyKerl(level-1)%blocks(index_i_loc,index_j_loc)%matrix, 1)
                  sizes%inds(ii,1)%size = mm
               enddo
               call BF_all2all_sizes(agent_block, sizes, ptree, level, 'R', 'C')


               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
               idx = idx_r
               inc = inc_r
               nblk_loc = nr
               agent_block%ButterflyU%idx = idx
               agent_block%ButterflyU%inc = inc
               agent_block%ButterflyU%nblk_loc = nblk_loc
               agent_block%ButterflyU%num_blk = 2**level_butterfly_loc
               allocate (agent_block%ButterflyU%blocks(nblk_loc))
               do ii = 1, nblk_loc
                  mm = sizes%inds(ii,1)%size
                  allocate (agent_block%ButterflyU%blocks(ii)%matrix(mm, mm))
                  agent_block%ButterflyU%blocks(ii)%matrix = 0
                  do i = 1, mm
                     agent_block%ButterflyU%blocks(ii)%matrix(i, i) = 1
                  end do
                  agent_block%M_loc =agent_block%M_loc + mm
               enddo
               deallocate(sizes%inds)
            else
               index_j_start = (ij_loc - 1)*2**(level_butterfly_loc - level + 1)

               call GetLocalBlockRange(ptree, pgno, level, level_butterfly_loc, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
               idx_c = idx_c*2 - 1
               inc_c = inc_c
               nc = nc*2
               if(nr>0 .and. nc>0)then
                  allocate (agent_block%ButterflyKerl(level)%blocks(nr, nc))

                  agent_block%ButterflyKerl(level)%num_row = 2**level
                  agent_block%ButterflyKerl(level)%num_col = 2**(level_butterfly_loc - level + 1)
                  agent_block%ButterflyKerl(level)%idx_r = idx_r
                  agent_block%ButterflyKerl(level)%idx_c = idx_c
                  agent_block%ButterflyKerl(level)%inc_r = inc_r
                  agent_block%ButterflyKerl(level)%inc_c = inc_c
                  agent_block%ButterflyKerl(level)%nr = nr
                  agent_block%ButterflyKerl(level)%nc = nc

                  do ii = 1, nr
                     do jj = 1, nc
                        index_i = (ii - 1)*inc_r + idx_r
                        index_j = (jj - 1)*inc_c + idx_c

                        index_i_loc = (index_i - block_o%ButterflyKerl(level)%idx_r)/block_o%ButterflyKerl(level)%inc_r + 1 !index_i_loc is local index in block_o
                        index_j_loc = (index_j+index_j_start - block_o%ButterflyKerl(level)%idx_c)/block_o%ButterflyKerl(level)%inc_c + 1
                        nn = size(block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix, 2)
                        mm = size(block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix, 1)
                        allocate (agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix(mm, nn))
                        agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix = block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix
                     enddo
                  enddo
               endif
            endif
         enddo
      end if

      ! n2=MPI_Wtime()
      ! time_tmp = time_tmp+n2-n1


      ! n1=MPI_Wtime()


      !>*********** compute M_p, N_p, ms, ns, M and N
      allocate(agent_block%M_p(ptree%pgrp(agent_block%pgno)%nproc,2))
      allocate(agent_block%N_p(ptree%pgrp(agent_block%pgno)%nproc,2))

#ifdef HAVE_MPI3
      call MPI_IALLGATHER(agent_block%M_loc, 1, MPI_INTEGER, agent_block%M_p, 1, MPI_INTEGER, ptree%pgrp(agent_block%pgno)%Comm, reqm, ierr)
      call MPI_IALLGATHER(agent_block%N_loc, 1, MPI_INTEGER, agent_block%N_p, 1, MPI_INTEGER, ptree%pgrp(agent_block%pgno)%Comm, reqn, ierr)
      call MPI_Wait(reqm, statusm, ierr)
      call MPI_Wait(reqn, statusn, ierr)
#else
      call MPI_ALLGATHER(agent_block%M_loc, 1, MPI_INTEGER, agent_block%M_p, 1, MPI_INTEGER, ptree%pgrp(agent_block%pgno)%Comm, ierr)
      call MPI_ALLGATHER(agent_block%N_loc, 1, MPI_INTEGER, agent_block%N_p, 1, MPI_INTEGER, ptree%pgrp(agent_block%pgno)%Comm, ierr)
#endif

      agent_block%M = 0
      agent_block%N = 0
      do pp=1,ptree%pgrp(agent_block%pgno)%nproc
         agent_block%M = agent_block%M + agent_block%M_p(pp,1)
         agent_block%N = agent_block%N + agent_block%N_p(pp,1)
         if(pp==1)then
            agent_block%M_p(pp,2) = agent_block%M_p(pp,1)
            agent_block%M_p(pp,1) = 1
            agent_block%N_p(pp,2) = agent_block%N_p(pp,1)
            agent_block%N_p(pp,1) = 1
         else
            agent_block%M_p(pp,2) = agent_block%M_p(pp-1,2) + agent_block%M_p(pp,1)
            agent_block%M_p(pp,1) = agent_block%M_p(pp-1,2) + 1
            agent_block%N_p(pp,2) = agent_block%N_p(pp-1,2) + agent_block%N_p(pp,1)
            agent_block%N_p(pp,1) = agent_block%N_p(pp-1,2) + 1
         endif
      enddo

      allocate(agent_block%ms(agent_block%ButterflyU%nblk_loc))
      do ii=1,agent_block%ButterflyU%nblk_loc
         if(ii==1)then
            agent_block%ms(ii) = size(agent_block%ButterflyU%blocks(ii)%matrix,1)
         else
            agent_block%ms(ii) = agent_block%ms(ii-1) + size(agent_block%ButterflyU%blocks(ii)%matrix,1)
         endif
      enddo
      allocate(agent_block%ns(agent_block%ButterflyV%nblk_loc))
      do jj=1,agent_block%ButterflyV%nblk_loc
         if(jj==1)then
            agent_block%ns(jj) = size(agent_block%ButterflyV%blocks(jj)%matrix,1)
         else
            agent_block%ns(jj) = agent_block%ns(jj-1) + size(agent_block%ButterflyV%blocks(jj)%matrix,1)
         endif
      enddo

      call BF_get_rank(agent_block, ptree)
      n2=MPI_Wtime()
      ! time_tmp = time_tmp+n2-n1
   end subroutine BF_extract_partial


   subroutine BF_copyback_partial(block_o, level_butterfly_loc, ij_loc, LR, agent_block,pgno,ptree)


      implicit none

      type(proctree)::ptree
      type(matrixblock)::block_o, agent_block
      integer level_butterfly, level_butterfly_loc, ij_loc, index_i, index_i_start, index_j_start, index_j, level, i, ii,jj, nn, mm, num_blocks, rank, pgno, index_i_loc, index_j_loc
      character LR
      integer idx_r, inc_r, nr, idx_c, inc_c, nc, nblk_loc, idx, inc, pp, ierr


      call assert(level_butterfly_loc >= 1, 'level_butterfly_loc cannot be zero')
      level_butterfly = block_o%level_butterfly

      if (LR == 'L') then
         do level = 0, level_butterfly_loc+1
            if(level==0)then

            elseif(level==level_butterfly_loc+1)then
               index_i_start = (ij_loc - 1)*2**level_butterfly_loc
               idx = agent_block%ButterflyU%idx
               inc = agent_block%ButterflyU%inc
               nblk_loc = agent_block%ButterflyU%nblk_loc
               do ii = 1, nblk_loc
                  index_i = (ii - 1)*inc + idx
                  index_i_loc = (index_i+index_i_start - block_o%ButterflyU%idx)/block_o%ButterflyU%inc + 1 !index_i_loc is local index in block_o
                  mm = size(agent_block%ButterflyU%blocks(ii)%matrix, 1)
                  nn = size(agent_block%ButterflyU%blocks(ii)%matrix, 2)
                  deallocate(block_o%ButterflyU%blocks(index_i_loc)%matrix)
                  allocate (block_o%ButterflyU%blocks(index_i_loc)%matrix(mm, nn))
                  block_o%ButterflyU%blocks(index_i_loc)%matrix = agent_block%ButterflyU%blocks(ii)%matrix
               enddo
            else
               index_i_start = (ij_loc - 1)*2**level

               idx_r = agent_block%ButterflyKerl(level)%idx_r
               idx_c = agent_block%ButterflyKerl(level)%idx_c
               inc_r = agent_block%ButterflyKerl(level)%inc_r
               inc_c = agent_block%ButterflyKerl(level)%inc_c
               nr = agent_block%ButterflyKerl(level)%nr
               nc = agent_block%ButterflyKerl(level)%nc

               do ii = 1, nr
                  do jj = 1, nc
                     index_i = (ii - 1)*inc_r + idx_r
                     index_j = (jj - 1)*inc_c + idx_c

                     index_i_loc = (index_i+index_i_start - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%idx_r)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%inc_r + 1 !index_i_loc is local index in block_o
                     index_j_loc = (index_j - block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%idx_c)/block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%inc_c + 1
                     nn = size(agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix, 2)
                     mm = size(agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix, 1)
                     deallocate(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix)
                     allocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix(mm, nn))
                     block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i_loc, index_j_loc)%matrix = agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix
                  enddo
               enddo
            endif
         enddo
      else if (LR == 'R') then
         do level = 0, level_butterfly_loc+1
            if(level==0)then
               index_j_start = (ij_loc - 1)*2**(level_butterfly_loc)
               idx = agent_block%ButterflyV%idx
               inc = agent_block%ButterflyV%inc
               nblk_loc = agent_block%ButterflyV%nblk_loc

               do jj = 1, nblk_loc
                  index_j = (jj - 1)*inc + idx
                  index_j_loc = (index_j+index_j_start - block_o%ButterflyV%idx)/block_o%ButterflyV%inc + 1 !index_i_loc is local index in block_o
                  mm = size(agent_block%ButterflyV%blocks(jj)%matrix, 1)
                  nn = size(agent_block%ButterflyV%blocks(jj)%matrix, 2)
                  deallocate (block_o%ButterflyV%blocks(index_j_loc)%matrix)
                  allocate (block_o%ButterflyV%blocks(index_j_loc)%matrix(mm, nn))
                  block_o%ButterflyV%blocks(index_j_loc)%matrix = agent_block%ButterflyV%blocks(jj)%matrix
               enddo

            elseif(level==level_butterfly_loc+1)then
            else
               index_j_start = (ij_loc - 1)*2**(level_butterfly_loc - level + 1)

               idx_r = agent_block%ButterflyKerl(level)%idx_r
               idx_c = agent_block%ButterflyKerl(level)%idx_c
               inc_r = agent_block%ButterflyKerl(level)%inc_r
               inc_c = agent_block%ButterflyKerl(level)%inc_c
               nr = agent_block%ButterflyKerl(level)%nr
               nc = agent_block%ButterflyKerl(level)%nc

               do ii = 1, nr
                  do jj = 1, nc
                     index_i = (ii - 1)*inc_r + idx_r
                     index_j = (jj - 1)*inc_c + idx_c

                     index_i_loc = (index_i - block_o%ButterflyKerl(level)%idx_r)/block_o%ButterflyKerl(level)%inc_r + 1 !index_i_loc is local index in block_o
                     index_j_loc = (index_j+index_j_start - block_o%ButterflyKerl(level)%idx_c)/block_o%ButterflyKerl(level)%inc_c + 1
                     nn = size(agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix, 2)
                     mm = size(agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix, 1)
                     deallocate(block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix)
                     allocate (block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix(mm, nn))
                     block_o%ButterflyKerl(level)%blocks(index_i_loc, index_j_loc)%matrix = agent_block%ButterflyKerl(level)%blocks(ii,jj)%matrix
                  enddo
               enddo
            endif
         enddo
      end if

      block_o%rankmax = max(block_o%rankmax,agent_block%rankmax)
   end subroutine BF_copyback_partial



   subroutine BF_copy_partial(block_i, block_o, level_butterfly_loc, ij_loc, LR, memory)


      implicit none
      type(matrixblock)::block_o, block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m, dimension_m, dimension_n, index_i_start, index_j_start
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly, level_butterfly_loc, ij_loc
      character LR
      real(kind=8), optional::memory
      if (present(memory)) memory = 0

!!!!! be careful here, may need changes later
      block_o%rankmax = max(block_o%rankmax, block_i%rankmax)
      block_o%rankmin = max(block_o%rankmin, block_i%rankmin)

      call assert(level_butterfly_loc >= 1, 'level_butterfly_loc cannot be zero')
      call assert(level_butterfly_loc == block_i%level_butterfly, 'level_butterfly_loc/=block_i%level_butterfly')

      level_butterfly = block_o%level_butterfly
      num_blocks = 2**level_butterfly_loc

      if (LR == 'L') then

         do level = 1, level_butterfly_loc
            do index_i = 1, 2**level
               do index_j = 1, 2**(level_butterfly_loc - level)
                  index_i_start = (ij_loc - 1)*2**level

                  if (level == 1) then
                     dimension_n = size(block_i%ButterflyV%blocks(2*index_j - 1)%matrix, 1)
                     nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 2)
                     rank = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 1)
                     deallocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix)
                     allocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix(rank, dimension_n))
                     ! call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j-1)%matrix, block_i%ButterflyV%blocks(2*index_j-1)%matrix, &
                     ! &block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j-1)%matrix,rank,dimension_n,nn)
                     call gemmf90(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, rank, block_i%ButterflyV%blocks(2*index_j - 1)%matrix, dimension_n, block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix, rank, 'N', 'T', rank, dimension_n, nn, BPACK_cone, BPACK_czero)

#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix, rank, dimension_n))) then
                        write (*, *) 'NAN in L 1'
                     end if
#endif
                     dimension_n = size(block_i%ButterflyV%blocks(2*index_j)%matrix, 1)
                     nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, 2)
                     rank = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, 1)
                     deallocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix)
                     allocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix(rank, dimension_n))
                     ! call gemmNT_omp(block_i%ButterflyKerl(level)%blocks(index_i,2*index_j)%matrix, block_i%ButterflyV%blocks(2*index_j)%matrix, &
                     ! &block_o%ButterflyKerl(level_butterfly-level_butterfly_loc+level)%blocks(index_i+index_i_start,2*index_j)%matrix,rank,dimension_n,nn)
                     call gemmf90(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, rank, block_i%ButterflyV%blocks(2*index_j)%matrix, dimension_n, block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix, rank, 'N', 'T', rank, dimension_n, nn, BPACK_cone, BPACK_czero)

#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix, rank, dimension_n))) then
                        write (*, *) 'NAN in L 2'
                     end if
#endif
                  else
                     nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 2)
                     rank = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix, 1)
                     deallocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix)
                     allocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix(rank, nn))
                     block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix = block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j - 1)%matrix
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix, rank, nn))) then
                        write (*, *) 'NAN in L 3'
                     end if
#endif
                     nn = size(block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix, 2)
                     deallocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix)
                     allocate (block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix(rank, nn))
                     block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix = block_i%ButterflyKerl(level)%blocks(index_i, 2*index_j)%matrix
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix, rank, nn))) then
                        write (*, *) 'NAN in L 4'
                     end if
#endif
                  end if

                  if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j - 1)%matrix)/1024.0d3
                  if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level_butterfly - level_butterfly_loc + level)%blocks(index_i + index_i_start, 2*index_j)%matrix)/1024.0d3

                  if (level == level_butterfly_loc) then
                     index_i_start = (ij_loc - 1)*2**level
                     mm = size(block_i%ButterflyU%blocks(index_i)%matrix, 1)
                     rank = size(block_i%ButterflyU%blocks(index_i)%matrix, 2)
                     deallocate (block_o%ButterflyU%blocks(index_i + index_i_start)%matrix)
                     allocate (block_o%ButterflyU%blocks(index_i + index_i_start)%matrix(mm, rank))
                     block_o%ButterflyU%blocks(index_i + index_i_start)%matrix = block_i%ButterflyU%blocks(index_i)%matrix
                     if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyU%blocks(index_i + index_i_start)%matrix)/1024.0d3
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyU%blocks(index_i + index_i_start)%matrix, mm, rank))) then
                        write (*, *) 'NAN in L 5'
                     end if
#endif
                  endif
               enddo
            enddo
         enddo

      else if (LR == 'R') then

         do level = 1, level_butterfly_loc
            do index_i = 1, 2**(level - 1)
               do index_j = 1, 2**(level_butterfly_loc - level + 1)
                  ! write(*,*)level,index_i,index_j
                  index_j_start = (ij_loc - 1)*2**(level_butterfly_loc - level + 1)
                  if (level == level_butterfly_loc) then
                     ! write(*,*)'good 1'
                     dimension_m = size(block_i%ButterflyU%blocks(2*index_i - 1)%matrix, 1)
                     mm = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, 1)
                     rank = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, 2)
                     ! write(*,*)dimension_m,mm,rank,'d'
                     deallocate (block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix)
                     allocate (block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix(dimension_m, rank))
                     ! call gemm_omp(block_i%ButterflyU%blocks(2*index_i-1)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i-1,index_j)%matrix,&
                     ! &block_o%ButterflyKerl(level)%blocks(2*index_i-1,index_j+index_j_start)%matrix,dimension_m,rank,mm)

                     call gemmf90(block_i%ButterflyU%blocks(2*index_i - 1)%matrix, dimension_m, block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, mm, block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix, dimension_m, 'N', 'N', dimension_m, rank, mm, BPACK_cone, BPACK_czero)

! write(*,*)'good 1.1'
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix, dimension_m, rank))) then
                        write (*, *) 'NAN in R 1'
                     end if
#endif
                     dimension_m = size(block_i%ButterflyU%blocks(2*index_i)%matrix, 1)
                     mm = size(block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix, 1)
                     rank = size(block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix, 2)
                     deallocate (block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix)
                     allocate (block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix(dimension_m, rank))
                     ! call gemm_omp(block_i%ButterflyU%blocks(2*index_i)%matrix, block_i%ButterflyKerl(level)%blocks(2*index_i,index_j)%matrix,&
                     ! &block_o%ButterflyKerl(level)%blocks(2*index_i,index_j+index_j_start)%matrix,dimension_m,rank,mm)

                     call gemmf90(block_i%ButterflyU%blocks(2*index_i)%matrix, dimension_m, block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix, mm, block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix, dimension_m, 'N', 'N', dimension_m, rank, mm, BPACK_cone, BPACK_czero)

! write(*,*)'good 2'
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix, dimension_m, rank))) then
                        write (*, *) 'NAN in R 2'
                     end if
#endif
                  else
                     ! write(*,*)'good 3'
                     mm = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, 1)
                     rank = size(block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix, 2)
                     deallocate (block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix)
                     allocate (block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix(mm, rank))
                     block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i - 1, index_j)%matrix

#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix, mm, rank))) then
                        write (*, *) 'NAN in R 3'
                     end if
#endif
                     mm = size(block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix, 1)
                     rank = size(block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix, 2)
                     deallocate (block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix)
                     allocate (block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix(mm, rank))
                     block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix = block_i%ButterflyKerl(level)%blocks(2*index_i, index_j)%matrix
                     ! write(*,*)'good 4'
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix, mm, rank))) then
                        write (*, *) 'NAN in R 4'
                     end if
#endif
                  end if

                  if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i - 1, index_j + index_j_start)%matrix)/1024.0d3
                  if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyKerl(level)%blocks(2*index_i, index_j + index_j_start)%matrix)/1024.0d3

                  if (level == 1) then
                     index_j_start = (ij_loc - 1)*2**(level_butterfly_loc - level + 1)
                     nn = size(block_i%ButterflyV%blocks(index_j)%matrix, 1)
                     rank = size(block_i%ButterflyV%blocks(index_j)%matrix, 2)
                     deallocate (block_o%ButterflyV%blocks(index_j + index_j_start)%matrix)
                     allocate (block_o%ButterflyV%blocks(index_j + index_j_start)%matrix(nn, rank))
                     block_o%ButterflyV%blocks(index_j + index_j_start)%matrix = block_i%ButterflyV%blocks(index_j)%matrix
                     if (present(memory)) memory = memory + SIZEOF(block_o%ButterflyV%blocks(index_j + index_j_start)%matrix)/1024.0d3
#ifndef NDEBUG
                     if (myisnan(fnorm(block_o%ButterflyV%blocks(index_j + index_j_start)%matrix, nn, rank))) then
                        write (*, *) 'NAN in R 5'
                     end if
#endif
                  endif
               end do
            end do
         end do

      end if

   end subroutine BF_copy_partial

   subroutine BF_Partial_MVP_Half(block_rand, chara, level_start, level_end, random, num_vect_sub, nth_s, nth_e, Ng)


      implicit none

      integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
      integer i, j, ii, jj, ij, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, level_end, level_start
      DT ctemp, a, b
      character chara
      integer num_vect_sub, num_vect_subsub, nth_s, nth_e, Ng, nth, dimension_rank, level_butterfly

      type(RandomBlock) :: random

      type(matrixblock)::block_rand
      ! write(*,*)'in '

      level_butterfly = block_rand%level_butterfly
      dimension_rank = block_rand%dimension_rank
      num_vect_subsub = num_vect_sub/(nth_e - nth_s + 1)

      if (chara == 'N') then

         num_blocks = 2**level_butterfly

         do level = level_start, level_end
            if (level == 0) then
               num_groupn = num_blocks
               do nth = nth_s, nth_e
                  ! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
                  do j = (nth - 1)*Ng + 1, nth*Ng
                     rank = size(block_rand%ButterflyV%blocks(j)%matrix, 2)
                     nn = size(block_rand%ButterflyV%blocks(j)%matrix, 1)
                     if (.not. associated(random%RandomVectorRR(1)%blocks(1, j)%matrix)) allocate (random%RandomVectorRR(1)%blocks(1, j)%matrix(rank, num_vect_subsub))
                     random%RandomVectorRR(1)%blocks(1, j)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                     !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                     do jj = 1, num_vect_subsub
                        do ii = 1, rank
                           ctemp = 0d0
                           do kk = 1, nn
                              ctemp = ctemp + block_rand%ButterflyV%blocks(j)%matrix(kk, ii)*random%RandomVectorRR(0)%blocks(1, j)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                           enddo
                           random%RandomVectorRR(1)%blocks(1, j)%matrix(ii, jj + (nth - nth_s)*num_vect_subsub) = ctemp
                        enddo
                     enddo
#ifdef HAVE_OPENMP
                     !$omp end parallel do
#endif
                  enddo
                  ! !$omp end parallel do
               enddo
            elseif (level == level_butterfly + 1) then

            else
               num_groupm = block_rand%ButterflyKerl(level)%num_row
               num_groupn = block_rand%ButterflyKerl(level)%num_col

               ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm,nth)
               do ij = 1, (num_groupm/2)*(num_groupn/2)
                  index_i = (ij - 1)/(num_groupn/2) + 1
                  index_j = mod(ij - 1, (num_groupn/2)) + 1
                  i = index_i*2 - 1
                  j = index_j*2 - 1

                  do nth = nth_s, nth_e

                     if ((j >= (nth - 1)*Ng/2**(level - 1) + 1 .and. j <= nth*Ng/2**(level - 1)) .or. &
                     & (j + 1 >= (nth - 1)*Ng/2**(level - 1) + 1 .and. j + 1 <= nth*Ng/2**(level - 1))) then
                        nn1 = size(block_rand%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
                        nn2 = size(block_rand%ButterflyKerl(level)%blocks(i, j + 1)%matrix, 2)
                        mm = size(block_rand%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                        ! write(*,*)ij,i,j,level,'ha',index_i
                        if (.not. associated(random%RandomVectorRR(level + 1)%blocks(i, index_j)%matrix)) allocate (random%RandomVectorRR(level + 1)%blocks(i, index_j)%matrix(mm, num_vect_subsub))
                        random%RandomVectorRR(level + 1)%blocks(i, index_j)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                        do jj = 1, num_vect_subsub
                           do ii = 1, mm
                              ctemp = 0d0
                              do kk = 1, nn1
                                 ctemp = ctemp + block_rand%ButterflyKerl(level)%blocks(i, j)%matrix(ii, kk)*random%RandomVectorRR(level)%blocks(index_i, j)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                              enddo
                              do kk = 1, nn2
                                 ctemp = ctemp + block_rand%ButterflyKerl(level)%blocks(i, j + 1)%matrix(ii, kk)*random%RandomVectorRR(level)%blocks(index_i, j + 1)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                              enddo
                              random%RandomVectorRR(level + 1)%blocks(i, index_j)%matrix(ii, jj + (nth - nth_s)*num_vect_subsub) = ctemp
                           enddo
                        enddo
#ifdef HAVE_OPENMP
                        !$omp end parallel do
#endif
                        nn1 = size(block_rand%ButterflyKerl(level)%blocks(i + 1, j)%matrix, 2)
                        nn2 = size(block_rand%ButterflyKerl(level)%blocks(i + 1, j + 1)%matrix, 2)
                        mm = size(block_rand%ButterflyKerl(level)%blocks(i + 1, j)%matrix, 1)
                        ! write(*,*)ij,i,j,level,'ha',index_i
                        if (.not. associated(random%RandomVectorRR(level + 1)%blocks(i + 1, index_j)%matrix)) allocate (random%RandomVectorRR(level + 1)%blocks(i + 1, index_j)%matrix(mm, num_vect_subsub))
                        random%RandomVectorRR(level + 1)%blocks(i + 1, index_j)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                        do jj = 1, num_vect_subsub
                           do ii = 1, mm
                              ctemp = 0d0
                              do kk = 1, nn1
                                 ctemp = ctemp + block_rand%ButterflyKerl(level)%blocks(i + 1, j)%matrix(ii, kk)*random%RandomVectorRR(level)%blocks(index_i, j)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                              enddo
                              do kk = 1, nn2
                                 ctemp = ctemp + block_rand%ButterflyKerl(level)%blocks(i + 1, j + 1)%matrix(ii, kk)*random%RandomVectorRR(level)%blocks(index_i, j + 1)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                              enddo
                              random%RandomVectorRR(level + 1)%blocks(i + 1, index_j)%matrix(ii, jj + (nth - nth_s)*num_vect_subsub) = ctemp
                           enddo
                        enddo
#ifdef HAVE_OPENMP
                        !$omp end parallel do
#endif
                        ! write(*,*)ij,i,j,level,'ha done0',index_i
                        deallocate (random%RandomVectorRR(level)%blocks(index_i, j)%matrix)
                        deallocate (random%RandomVectorRR(level)%blocks(index_i, j + 1)%matrix)
                        ! write(*,*)ij,i,j,level,'ha done',index_i
                     end if
                  end do
               enddo
               ! !$omp end parallel do
            endif
         enddo

      elseif (chara == 'T') then

         num_blocks = 2**level_butterfly

         do level = level_start, level_end
            if (level == 0) then
               num_groupm = num_blocks
               do nth = nth_s, nth_e
                  ! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
                  do i = (nth - 1)*Ng + 1, nth*Ng
                     rank = size(block_rand%ButterflyU%blocks(i)%matrix, 2)
                     mm = size(block_rand%ButterflyU%blocks(i)%matrix, 1)
                     if (.not. associated(random%RandomVectorLL(1)%blocks(i, 1)%matrix)) allocate (random%RandomVectorLL(1)%blocks(i, 1)%matrix(rank, num_vect_subsub))
                     random%RandomVectorLL(1)%blocks(i, 1)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                     !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                     do jj = 1, num_vect_subsub
                        do ii = 1, rank
                           ctemp = 0d0
                           do kk = 1, mm
                              ctemp = ctemp + block_rand%ButterflyU%blocks(i)%matrix(kk, ii)*random%RandomVectorLL(0)%blocks(i, 1)%matrix(kk, jj + (nth - nth_s)*num_vect_subsub)
                           enddo
                           random%RandomVectorLL(1)%blocks(i, 1)%matrix(ii, jj + (nth - nth_s)*num_vect_subsub) = ctemp
                        enddo
                     enddo
#ifdef HAVE_OPENMP
                     !$omp end parallel do
#endif
                  end do
                  ! !$omp end parallel do
               enddo
            elseif (level == level_butterfly + 1) then
            else
               num_groupm = block_rand%ButterflyKerl(level_butterfly - level + 1)%num_row
               num_groupn = block_rand%ButterflyKerl(level_butterfly - level + 1)%num_col

               ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn,nth)
               do ij = 1, (num_groupn/2)*(num_groupm/2)
                  index_j = (ij - 1)/(num_groupm/2) + 1
                  index_i = mod(ij - 1, (num_groupm/2)) + 1
                  j = 2*index_j - 1
                  i = 2*index_i - 1

                  do nth = nth_s, nth_e

                     if ((i >= (nth - 1)*Ng/2**(level - 1) + 1 .and. i <= nth*Ng/2**(level - 1)) .or. &
                     & (i + 1 >= (nth - 1)*Ng/2**(level - 1) + 1 .and. i + 1 <= nth*Ng/2**(level - 1))) then
                        mm1 = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j)%matrix, 1)
                        mm2 = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i + 1, j)%matrix, 1)
                        nn = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j)%matrix, 2)
                        if (.not. associated(random%RandomVectorLL(level + 1)%blocks(index_i, j)%matrix)) allocate (random%RandomVectorLL(level + 1)%blocks(index_i, j)%matrix(nn, num_vect_subsub))
                        random%RandomVectorLL(level + 1)%blocks(index_i, j)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                        do ii = 1, num_vect_subsub
                           do jj = 1, nn
                              ctemp = 0d0
                              do kk = 1, mm1
                                 ctemp = ctemp + random%RandomVectorLL(level)%blocks(i, index_j)%matrix(kk, ii + (nth - nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j)%matrix(kk, jj)
                              enddo
                              do kk = 1, mm2
                                 ctemp = ctemp + random%RandomVectorLL(level)%blocks(i + 1, index_j)%matrix(kk, ii + (nth - nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i + 1, j)%matrix(kk, jj)
                              enddo
                              random%RandomVectorLL(level + 1)%blocks(index_i, j)%matrix(jj, ii + (nth - nth_s)*num_vect_subsub) = ctemp
                           enddo
                        enddo
#ifdef HAVE_OPENMP
                        !$omp end parallel do
#endif
                        mm1 = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j + 1)%matrix, 1)
                        mm2 = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i + 1, j + 1)%matrix, 1)
                        nn = size(block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j + 1)%matrix, 2)
                        if (.not. associated(random%RandomVectorLL(level + 1)%blocks(index_i, j + 1)%matrix)) allocate (random%RandomVectorLL(level + 1)%blocks(index_i, j + 1)%matrix(nn, num_vect_subsub))
                        random%RandomVectorLL(level + 1)%blocks(index_i, j + 1)%matrix(:, (nth - nth_s)*num_vect_subsub + 1:(nth - nth_s + 1)*num_vect_subsub) = 0
#ifdef HAVE_OPENMP
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
#endif
                        do ii = 1, num_vect_subsub
                           do jj = 1, nn
                              ctemp = 0d0
                              do kk = 1, mm1
                                 ctemp = ctemp + random%RandomVectorLL(level)%blocks(i, index_j)%matrix(kk, ii + (nth - nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i, j + 1)%matrix(kk, jj)
                              enddo
                              do kk = 1, mm2
                                 ctemp = ctemp + random%RandomVectorLL(level)%blocks(i + 1, index_j)%matrix(kk, ii + (nth - nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly - level + 1)%blocks(i + 1, j + 1)%matrix(kk, jj)
                              enddo
                              random%RandomVectorLL(level + 1)%blocks(index_i, j + 1)%matrix(jj, ii + (nth - nth_s)*num_vect_subsub) = ctemp
                           enddo
                        enddo
#ifdef HAVE_OPENMP
                        !$omp end parallel do
#endif
                        deallocate (random%RandomVectorLL(level)%blocks(i, index_j)%matrix)
                        deallocate (random%RandomVectorLL(level)%blocks(i + 1, index_j)%matrix)
                     end if
                  end do
               enddo
               ! !$omp end parallel do

            endif
         enddo

      endif
      ! write(*,*)'out '
      return

   end subroutine BF_Partial_MVP_Half

   subroutine BF_exchange_extraction(blocks, kerls, stats, ptree, level, collect)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i1, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j1, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hstat)::stats
      type(proctree)::ptree
      type(butterfly_kerl)::kerls

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_subs, pgno_subr, pgno_sub, pids, pidr, pid, pid0, tag, nproc, Ncol, Ncol1, Nrow, Nreqr, Nreqs, recvid, sendid

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive
      logical::sendflag, recvflag
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer rr, cc
      character mode, modetrans, collect
      DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)

      real(kind=8)::n1, n2

      n1 = MPI_Wtime()

      mode = 'R'
      modetrans = 'C'

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno+level*10

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

      do nn = 1, size(kerls%index, 1)
         ii = kerls%index(nn, 1)
         jj = kerls%index(nn, 2)
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c

         sendflag = .false.
         recvflag = .false.
         if (collect == 'R') then ! pair-wise reduction
            if (mode == 'R') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = 2*index_j - mod(index_i, 2)
               index_i1 = floor_safe((index_i - 1)/2d0) + 1
               index_j1 = 2*index_j - mod(index_i - 1, 2)
            elseif (mode == 'C') then
               write (*, *) 'mode=C not needed in BF_exchange_extraction'
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, modetrans, pgno_subs)
            pids = ptree%pgrp(pgno_subs)%head
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i1, index_j1, modetrans, pgno_subr)
            pidr = ptree%pgrp(pgno_subr)%head
         elseif (collect == 'B') then ! pair-wise broadcast
            if (mode == 'R') then
               index_j0 = index_j + 2*mod(index_j, 2) - 1
               index_i0 = index_i
               index_j1 = index_j
               index_i1 = index_i
            elseif (mode == 'C') then
               write (*, *) 'mode=C not needed in BF_exchange_extraction'
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, mode, pgno_subs)
            pids = ptree%pgrp(pgno_subs)%head
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i1, index_j1, mode, pgno_subr)
            pidr = ptree%pgrp(pgno_subr)%head
         endif
         sendflag = pids /= ptree%MyID
         recvflag = pidr /= ptree%MyID

         if (recvflag) then
            pp = pidr - ptree%pgrp(blocks%pgno)%head + 1
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         endif

         if (sendflag) then
            pp = pids - ptree%pgrp(blocks%pgno)%head + 1
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            if (associated(kerls%blocks(ii, jj)%matrix)) then
               sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls%blocks(ii, jj)%matrix, 1)*size(kerls%blocks(ii, jj)%matrix, 2)
            endif
         endif
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid /=ptree%MyID)then
            Nreqs = Nreqs+1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /=ptree%MyID)then
            Nreqr = Nreqr+1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do nn = 1, size(kerls%index, 1)
         ii = kerls%index(nn, 1)
         jj = kerls%index(nn, 2)
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c

         sendflag = .false.
         if (collect == 'R') then ! pair-wise reduction
            if (mode == 'R') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = 2*index_j - mod(index_i, 2)
            elseif (mode == 'C') then
               write (*, *) 'mode=C not needed in BF_exchange_extraction'
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, modetrans, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         elseif (collect == 'B') then ! pair-wise broadcast

            if (mode == 'R') then
               index_j0 = index_j + 2*mod(index_j, 2) - 1
               index_i0 = index_i
            elseif (mode == 'C') then
               write (*, *) 'mode=C not needed in BF_exchange_extraction'
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
         endif
         sendflag = pid /= ptree%MyID

         if (sendflag) then
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            if (associated(kerls%blocks(ii, jj)%matrix)) then
               Nrow = size(kerls%blocks(ii, jj)%matrix, 1)
               Ncol = size(kerls%blocks(ii, jj)%matrix, 2)
               sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
               sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
               sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
               sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
               sendquant(pp)%size = sendquant(pp)%size + 4
               do i = 1, Nrow*Ncol
                  rr = mod(i - 1, Nrow) + 1
                  cc = (i - 1)/Nrow + 1
                  sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls%blocks(ii, jj)%matrix(rr, cc)
               enddo
               sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
            endif
         endif
      enddo

      ! communicate the data buffer
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid /=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if(sendquant(pp)%size>0)recvquant(pp)%dat= sendquant(pp)%dat
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - kerls%idx_r)/kerls%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jj = (index_j - kerls%idx_c)/kerls%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))

            if (associated(kerls%blocks(ii, jj)%matrix)) then
               Ncol1 = size(kerls%blocks(ii, jj)%matrix, 2)
               allocate (mat(Nrow, Ncol1 + Ncol))
               mat(1:Nrow, 1:Ncol1) = kerls%blocks(ii, jj)%matrix
               deallocate (kerls%blocks(ii, jj)%matrix)
            else
               Ncol1 = 0
               allocate (mat(Nrow, Ncol1 + Ncol))
            endif
            allocate (kerls%blocks(ii, jj)%matrix(Nrow, Ncol1 + Ncol))

            do j = 1, Nrow*Ncol
               rr = mod(j - 1, Nrow) + 1
               cc = (j - 1)/Nrow + 1
               mat(rr, cc + Ncol1) = recvquant(pp)%dat(i + j, 1)
            enddo
            kerls%blocks(ii, jj)%matrix = mat
            deallocate (mat)
            i = i + Nrow*Ncol
         enddo
      enddo

      ! deallocation
      deallocate (S_req)
      deallocate (R_req)
      deallocate (statuss)
      deallocate (statusr)
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)
      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_exchange_extraction

   subroutine BF_exchange_matvec(blocks, kerls, stats, ptree, level, mode, collect)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hstat)::stats
      type(proctree)::ptree
      type(butterfly_kerl)::kerls

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_sub, pgno_sub_mine, pid, pid0, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive
      logical::sendflag
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer rr, cc
      character mode, modetrans, collect

      real(kind=8)::n1, n2

      n1 = MPI_Wtime()

      if (mode == 'R') modetrans = 'C'
      if (mode == 'C') modetrans = 'R'

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno+level*10

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

      call GetPgno_Sub(ptree, blocks%pgno, level_butterfly, pgno_sub_mine)

      ! calculate send buffer sizes in the first pass
      do ii = 1, kerls%nr
      do jj = 1, kerls%nc
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c

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
         sendflag = pid /= ptree%MyID .and. (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID)

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
            sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls%blocks(ii, jj)%matrix, 1)*size(kerls%blocks(ii, jj)%matrix, 2)
         endif
      enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid /= ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do ii = 1, kerls%nr
      do jj = 1, kerls%nc

         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c

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
         sendflag = pid /= ptree%MyID .and. (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID)

         if (sendflag) then
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            Nrow = size(kerls%blocks(ii, jj)%matrix, 1)
            Ncol = size(kerls%blocks(ii, jj)%matrix, 2)
            sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
            sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
            sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
            sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
            sendquant(pp)%size = sendquant(pp)%size + 4
            do i = 1, Nrow*Ncol
               rr = mod(i - 1, Nrow) + 1
               cc = (i - 1)/Nrow + 1
               sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls%blocks(ii, jj)%matrix(rr, cc)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
         endif
      enddo
      enddo

      ! communicate the data buffer
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid /= ptree%MyID)then
            Nreqs=Nreqs+1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if(sendquant(pp)%size>0) recvquant(pp)%dat=sendquant(pp)%dat
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr=Nreqr+1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - kerls%idx_r)/kerls%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jj = (index_j - kerls%idx_c)/kerls%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            if (.not. associated(kerls%blocks(ii, jj)%matrix)) then
               allocate (kerls%blocks(ii, jj)%matrix(Nrow, Ncol))
               kerls%blocks(ii, jj)%matrix = 0
            endif
            do j = 1, Nrow*Ncol
               rr = mod(j - 1, Nrow) + 1
               cc = (j - 1)/Nrow + 1
               kerls%blocks(ii, jj)%matrix(rr, cc) = kerls%blocks(ii, jj)%matrix(rr, cc) + recvquant(pp)%dat(i + j, 1)
            enddo
            i = i + Nrow*Ncol
         enddo
      enddo

      ! deallocation
      deallocate (S_req)
      deallocate (R_req)
      deallocate (statuss)
      deallocate (statusr)
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)
      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_exchange_matvec



!>*********** all to all communication of sizes of one butterfly level from row-wise ordering to column-wise ordering or the reverse
   subroutine BF_all2all_sizes(blocks, sizes, ptree, level, mode, mode_new)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
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
      type(butterfly_skel)::sizes
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      integer, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist

      n1 = MPI_Wtime()

      call assert(mode /= mode_new, 'only row2col or col2row is supported')

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno

      if(nproc>1)then

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

      do ii = 1, sizes%nr
      do jj = 1, sizes%nc
         index_i = (ii - 1)*sizes%inc_r + sizes%idx_r
         index_j = (jj - 1)*sizes%inc_c + sizes%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         if (sendquant(pp)%active == 0) then
            sendquant(pp)%active = 1
            Nsendactive = Nsendactive + 1
            sendIDactive(Nsendactive) = pp
         endif
         sendquant(pp)%size_i = sendquant(pp)%size_i + 3
      enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs = 0
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
      do ii = 1, sizes%nr
      do jj = 1, sizes%nc
         index_i = (ii - 1)*sizes%inc_r + sizes%idx_r
         index_j = (jj - 1)*sizes%inc_c + sizes%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head

         pp = pid - ptree%pgrp(blocks%pgno)%head + 1

         Nskel = sizes%inds(ii, jj)%size
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 1, 1) = index_i
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 2, 1) = index_j
         sendquant(pp)%dat_i(sendquant(pp)%size_i + 3, 1) = Nskel
         sendquant(pp)%size_i = sendquant(pp)%size_i + 3
      enddo
      enddo

      deallocate (sizes%inds)
      sizes%idx_r = idx_r
      sizes%idx_c = idx_c
      sizes%inc_r = inc_r
      sizes%inc_c = inc_c
      sizes%nr = nr
      sizes%nc = nc

      allocate (sizes%inds(sizes%nr, sizes%nc))

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
            ii = (index_i - sizes%idx_r)/sizes%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            jj = (index_j - sizes%idx_c)/sizes%inc_c + 1
            i = i + 1
            Nskel = NINT(dble(recvquant(pp)%dat_i(i, 1)))
            sizes%inds(ii, jj)%size = Nskel
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
      endif


      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_sizes



! !>*********** all to all communication of sizes of one butterfly level from row-wise ordering to column-wise ordering or the reverse
!    subroutine BF_all2all_sizes(blocks, sizes, ptree, nproc, level, mode, mode_new)


!       implicit none
!       integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
!       integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
!       integer level, length_1, length_2, level_blocks
!       integer rank, rankmax, butterflyB_inuse, rank1, rank2
!       real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
!       integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
!       real(kind=8) flop
!       DT ctemp
!       type(matrixblock)::blocks
!       type(proctree)::ptree

!       integer ierr, nsendrecv, pid, pgno_sub, tag, nproc, Ncol, Nskel, Nreqr, Nreqs, recvid, sendid, tmpi
!       integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

!       type(commquant1D)::sendquant(nproc), recvquant(nproc)
!       integer ::S_req(nproc), R_req(nproc)
!       integer :: statuss(MPI_status_size, nproc), statusr(MPI_status_size, nproc)
!       character::mode, mode_new
!       real(kind=8)::n1, n2
!       integer::sendIDactive(nproc), recvIDactive(nproc)
!       integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
!       logical all2all
!       type(butterfly_skel)::sizes
!       integer::dist


!       call assert(mode /= mode_new, 'only row2col or col2row is supported')

!       level_butterfly = blocks%level_butterfly
!       ! nproc = ptree%pgrp(blocks%pgno)%nproc
!       tag = blocks%pgno

!       if(nproc>1)then

!       ! mode_new and level_new determine the block range in the new mode
!       if (mode_new == 'R') then
!          level_new = max(level - 1, 0)
!       elseif (mode_new == 'C') then
!          level_new = min(level + 1, level_butterfly + 1)
!       endif
!       call GetLocalBlockRange(ptree, blocks%pgno, level_new, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, mode_new)

!       ! allocation of communication quantities
!       do ii = 1, nproc
!          sendquant(ii)%size = 0
!          sendquant(ii)%active = 0
!       enddo
!       do ii = 1, nproc
!          recvquant(ii)%size = 0
!          recvquant(ii)%active = 0
!       enddo
!       Nsendactive = 0
!       Nrecvactive = 0

!       ! calculate send buffer sizes in the first pass
!       do ii = 1, nr
!       do jj = 1, nc
!          index_i = (ii - 1)*inc_r + idx_r
!          index_j = (jj - 1)*inc_c + idx_c
!          call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i, index_j, mode, pgno_sub)
!          pid = ptree%pgrp(pgno_sub)%head
!          pp = pid - ptree%pgrp(blocks%pgno)%head + 1
!          if (recvquant(pp)%active == 0) then
!             recvquant(pp)%active = 1
!             Nrecvactive = Nrecvactive + 1
!             recvIDactive(Nrecvactive) = pp
!          endif
!       enddo
!       enddo

!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(tt), ierr)
!       enddo

!       do ii = 1, sizes%nr
!       do jj = 1, sizes%nc
!          index_i = (ii - 1)*sizes%inc_r + sizes%idx_r
!          index_j = (jj - 1)*sizes%inc_c + sizes%idx_c
!          call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
!          pid = ptree%pgrp(pgno_sub)%head
!          pp = pid - ptree%pgrp(blocks%pgno)%head + 1
!          if (sendquant(pp)%active == 0) then
!             sendquant(pp)%active = 1
!             Nsendactive = Nsendactive + 1
!             sendIDactive(Nsendactive) = pp
!          endif
!          sendquant(pp)%size = sendquant(pp)%size + 3
!       enddo
!       enddo


!       n1 = MPI_Wtime()

!       ! communicate receive buffer sizes
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          allocate (sendquant(pp)%dat_i(sendquant(pp)%size, 1))
!          recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(tt), ierr)
!       enddo



!       if (Nrecvactive > 0) then
!          call MPI_waitall(Nrecvactive, R_req, statusr, ierr)
!       endif
!       n2 = MPI_Wtime()
!       ! time_tmp = time_tmp + n2 - n1

!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          sendquant(pp)%size = 0
!       enddo
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          allocate (recvquant(pp)%dat_i(recvquant(pp)%size, 1))
!       enddo

!       n1 = MPI_Wtime()

!       ! pack the send buffer in the second pass
!       do ii = 1, sizes%nr
!       do jj = 1, sizes%nc
!          index_i = (ii - 1)*sizes%inc_r + sizes%idx_r
!          index_j = (jj - 1)*sizes%inc_c + sizes%idx_c
!          call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
!          pid = ptree%pgrp(pgno_sub)%head

!          pp = pid - ptree%pgrp(blocks%pgno)%head + 1

!          Nskel = sizes%inds(ii, jj)%size
!          sendquant(pp)%dat_i(sendquant(pp)%size + 1, 1) = index_i
!          sendquant(pp)%dat_i(sendquant(pp)%size + 2, 1) = index_j
!          sendquant(pp)%dat_i(sendquant(pp)%size + 3, 1) = Nskel
!          sendquant(pp)%size = sendquant(pp)%size + 3
!       enddo
!       enddo




!       deallocate (sizes%inds)
!       sizes%idx_r = idx_r
!       sizes%idx_c = idx_c
!       sizes%inc_r = inc_r
!       sizes%inc_c = inc_c
!       sizes%nr = nr
!       sizes%nc = nc

!       allocate (sizes%inds(sizes%nr, sizes%nc))
!       if (Nsendactive > 0) then
!          call MPI_waitall(Nsendactive, S_req, statuss, ierr)
!       endif
!       ! communicate the data buffer
!       Nreqs = 0
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          if (recvid /= ptree%MyID) then
!             Nreqs = Nreqs + 1
!             call MPI_Isend(sendquant(pp)%dat_i, sendquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
!          else
!             if (sendquant(pp)%size > 0) recvquant(pp)%dat_i = sendquant(pp)%dat_i
!          endif
!       enddo

!       Nreqr = 0
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          if (sendid /= ptree%MyID) then
!             Nreqr = Nreqr + 1
!             call MPI_Irecv(recvquant(pp)%dat_i, recvquant(pp)%size, MPI_INTEGER, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
!          endif
!       enddo

!       ! copy data from buffer to target
!       do tt = 1, Nrecvactive
!          if(tt==1 .and. Nreqr+1== Nrecvactive)then
!             pp = ptree%MyID + 1 - ptree%pgrp(blocks%pgno)%head
!          else
!             call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
!             pp = statusr(MPI_SOURCE, 1) + 1
!          endif
!          ! n1 = MPI_Wtime()
!          i = 0
!          do while (i < recvquant(pp)%size)
!             i = i + 1
!             index_i = NINT(dble(recvquant(pp)%dat_i(i, 1)))
!             ii = (index_i - sizes%idx_r)/sizes%inc_r + 1
!             i = i + 1
!             index_j = NINT(dble(recvquant(pp)%dat_i(i, 1)))
!             jj = (index_j - sizes%idx_c)/sizes%inc_c + 1
!             i = i + 1
!             Nskel = NINT(dble(recvquant(pp)%dat_i(i, 1)))
!             sizes%inds(ii, jj)%size = Nskel
!          enddo
!          ! n2 = MPI_Wtime()
!          ! time_tmp = time_tmp + n2 - n1
!       enddo
!       if (Nreqs > 0) then
!          call MPI_waitall(Nreqs, S_req, statuss, ierr)
!       endif

!       ! deallocation
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          if (allocated(sendquant(pp)%dat_i)) deallocate (sendquant(pp)%dat_i)
!       enddo
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          if (allocated(recvquant(pp)%dat_i)) deallocate (recvquant(pp)%dat_i)
!       enddo
!       endif

!       n2 = MPI_Wtime()
!       ! time_tmp = time_tmp + n2 - n1


!    end subroutine BF_all2all_sizes




!>*********** all to all communication of extraction results of one butterfly level from row-wise ordering to column-wise ordering or the reverse
   subroutine BF_all2all_extraction(blocks, kerls, kerls1, stats, ptree, level, mode, mode_new)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, pgno_sub, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      character::mode, mode_new
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls, kerls1
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist

      n1 = MPI_Wtime()

      call assert(mode /= mode_new, 'only row2col or col2row is supported')

      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno+level*10

      ! mode_new and level_new determine the block range in the new mode
      if (mode_new == 'R') then
         level_new = max(level - 1, 0)
      elseif (mode_new == 'C') then
         level_new = min(level + 1, level_butterfly + 1)
      endif

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

      do nn = 1, size(kerls1%index, 1)
         ii = kerls1%index(nn, 1)
         jj = kerls1%index(nn, 2)
         index_i = (ii - 1)*kerls1%inc_r + kerls1%idx_r
         index_j = (jj - 1)*kerls1%inc_c + kerls1%idx_c
         call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i, index_j, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         if (recvquant(pp)%active == 0) then
            recvquant(pp)%active = 1
            Nrecvactive = Nrecvactive + 1
            recvIDactive(Nrecvactive) = pp
         endif
      enddo

      do nn = 1, size(kerls%index, 1)
         ii = kerls%index(nn, 1)
         jj = kerls%index(nn, 2)
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         if (sendquant(pp)%active == 0) then
            sendquant(pp)%active = 1
            Nsendactive = Nsendactive + 1
            sendIDactive(Nsendactive) = pp
         endif
         sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls%blocks(ii, jj)%matrix, 1)*size(kerls%blocks(ii, jj)%matrix, 2)
      enddo

      ! communicate receive buffer sizes
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid /=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid /=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do nn = 1, size(kerls%index, 1)
         ii = kerls%index(nn, 1)
         jj = kerls%index(nn, 2)
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head

         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         Nrow = size(kerls%blocks(ii, jj)%matrix, 1)
         Ncol = size(kerls%blocks(ii, jj)%matrix, 2)

         sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
         sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
         sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
         sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
         sendquant(pp)%size = sendquant(pp)%size + 4
         do i = 1, Nrow*Ncol
            rr = mod(i - 1, Nrow) + 1
            cc = (i - 1)/Nrow + 1
            sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls%blocks(ii, jj)%matrix(rr, cc)
         enddo
         sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
         deallocate (kerls%blocks(ii, jj)%matrix)
         if (allocated(kerls%blocks(ii, jj)%index)) deallocate (kerls%blocks(ii, jj)%index)
      enddo

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - kerls1%idx_r)/kerls1%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jj = (index_j - kerls1%idx_c)/kerls1%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(kerls1%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
            allocate (kerls1%blocks(ii, jj)%matrix(Nrow, Ncol))
            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  kerls1%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
               enddo
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)


      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_extraction


   !>***** switching the matvecs/temporary buffer/kernels from row/col distributions to col/row distributions, kerflag=1: kernels, kerflag=0: matvecs and buffer
   subroutine BF_all2all_vec_n_ker(blocks, kerls, stats, ptree, nproc, level, mode, mode_new, kerflag)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i,index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_sub, pgno_sub_mine, pid, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      integer::sendIDactive(nproc), recvIDactive(nproc)
      character::mode, mode_new
      real(kind=8)::n1, n2
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist, kerflag

      n1 = MPI_Wtime()

      call assert(mode /= mode_new, 'only row2col or col2row is supported')

      level_butterfly = blocks%level_butterfly
      ! nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno+level*10

      ! mode_new and level_new determine the block range in the new mode
      if(kerflag==0)then
         if (mode_new == 'R') then
            level_new = max(level - 1, 0)
         elseif (mode_new == 'C') then
            level_new = min(level + 1, level_butterfly + 1)
         endif
      else
         level_new=level
      endif

      call GetLocalBlockRange(ptree, blocks%pgno, level_new, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, mode_new)

      if(kerflag==1)then
         if (mode_new == 'R') then
            idx_c = idx_c*2 - 1
            inc_c = inc_c
            nc = nc*2
         elseif (mode_new == 'C') then
            idx_r = idx_r*2 - 1
            inc_r = inc_r
            nr = nr*2
         endif
      endif


      ! allocation of communication quantities
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo

      Nsendactive = 0
      Nrecvactive = 0

      call GetPgno_Sub(ptree, blocks%pgno, level_butterfly, pgno_sub_mine)

      if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
         ! calculate send buffer sizes in the first pass
         do ii = 1, nr
         do jj = 1, nc
            index_i = (ii - 1)*inc_r + idx_r
            index_j = (jj - 1)*inc_c + idx_c
            index_i0 = index_i
            index_j0 = index_j
            if(kerflag==1)then
            if (mode == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
            endif
            if (mode == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
            endif
            endif
            call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         enddo
         enddo

         do ii = 1, kerls%nr
         do jj = 1, kerls%nc
            index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
            index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
            index_i0 = index_i
            index_j0 = index_j
            if(kerflag==1)then
            if (mode_new == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
            endif
            if (mode_new == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
            endif
            endif
            call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i0, index_j0, mode_new, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            pp = pid - ptree%pgrp(blocks%pgno)%head + 1
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls%blocks(ii, jj)%matrix, 1)*size(kerls%blocks(ii, jj)%matrix, 2)
         enddo
         enddo
      endif

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
      do ii = 1, kerls%nr
      do jj = 1, kerls%nc
         index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
         index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
         index_i0 = index_i
         index_j0 = index_j
         if(kerflag==1)then
            if (mode_new == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
            endif
            if (mode_new == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
            endif
         endif
         call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i0, index_j0, mode_new, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head

         pp = pid - ptree%pgrp(blocks%pgno)%head + 1
         Nrow = size(kerls%blocks(ii, jj)%matrix, 1)
         Ncol = size(kerls%blocks(ii, jj)%matrix, 2)

         sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
         sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
         sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
         sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
         sendquant(pp)%size = sendquant(pp)%size + 4
         do i = 1, Nrow*Ncol
            rr = mod(i - 1, Nrow) + 1
            cc = (i - 1)/Nrow + 1
            sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls%blocks(ii, jj)%matrix(rr, cc)
         enddo
         sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
         deallocate (kerls%blocks(ii, jj)%matrix)
      enddo
      enddo
      deallocate (kerls%blocks)

      kerls%idx_r = idx_r
      kerls%idx_c = idx_c
      kerls%inc_r = inc_r
      kerls%inc_c = inc_c
      kerls%nr = nr
      kerls%nc = nc

      allocate (kerls%blocks(kerls%nr, kerls%nc))
      endif


      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - kerls%idx_r)/kerls%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jj = (index_j - kerls%idx_c)/kerls%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(kerls%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
            allocate (kerls%blocks(ii, jj)%matrix(Nrow, Ncol))
            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  kerls%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
               enddo
            enddo
         enddo
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_vec_n_ker


! !>*********** all to all communication of matvec results of one butterfly level from row-wise ordering to column-wise ordering or the reverse
!    subroutine BF_all2all_vec_n_ker(blocks, kerls, stats, ptree, nproc, level, mode, mode_new)


!       implicit none
!       integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
!       integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
!       integer level, length_1, length_2, level_blocks
!       integer rank, rankmax, butterflyB_inuse, rank1, rank2
!       real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
!       integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
!       real(kind=8) flop
!       DT ctemp
!       type(matrixblock)::blocks
!       type(Hstat)::stats
!       type(proctree)::ptree

!       integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
!       integer, allocatable:: select_row_rr(:), select_column_rr(:)
!       DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

!       integer, allocatable::jpvt(:)
!       integer ierr, nsendrecv, pgno_sub, pgno_sub_mine, pid, tag, nproc, Ncol, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
!       integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new


!       type(commquant1D)::sendquant(nproc), recvquant(nproc)
!       integer::sendactive(nproc), recvactive(nproc)
!       integer::S_req(nproc),R_req(nproc)
!       integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
!       integer::sendIDactive(nproc), recvIDactive(nproc)
!       character::mode, mode_new
!       real(kind=8)::n1, n2
!       integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
!       type(butterfly_kerl)::kerls
!       integer rr, cc, cnt
!       logical all2all
!       integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
!       DT, allocatable::sendbufall2all(:), recvbufall2all(:)
!       integer::dist

!       n1 = MPI_Wtime()

!       call assert(mode /= mode_new, 'only row2col or col2row is supported')

!       level_butterfly = blocks%level_butterfly
!       ! nproc = ptree%pgrp(blocks%pgno)%nproc
!       tag = blocks%pgno+level*10

!       ! mode_new and level_new determine the block range in the new mode
!       if (mode_new == 'R') then
!          level_new = max(level - 1, 0)
!       elseif (mode_new == 'C') then
!          level_new = min(level + 1, level_butterfly + 1)
!       endif
!       call GetLocalBlockRange(ptree, blocks%pgno, level_new, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, mode_new)

!       ! allocation of communication quantities
!       do ii = 1, nproc
!          sendquant(ii)%size = 0
!          sendquant(ii)%active = 0
!       enddo
!       do ii = 1, nproc
!          recvquant(ii)%size = 0
!          recvquant(ii)%active = 0
!       enddo
!       Nsendactive = 0
!       Nrecvactive = 0

!       call GetPgno_Sub(ptree, blocks%pgno, level_butterfly, pgno_sub_mine)

!       if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
!          ! calculate send buffer sizes in the first pass
!          do ii = 1, nr
!          do jj = 1, nc
!             index_i = (ii - 1)*inc_r + idx_r
!             index_j = (jj - 1)*inc_c + idx_c
!             call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i, index_j, mode, pgno_sub)
!             pid = ptree%pgrp(pgno_sub)%head
!             pp = pid - ptree%pgrp(blocks%pgno)%head + 1
!             if (recvquant(pp)%active == 0) then
!                recvquant(pp)%active = 1
!                Nrecvactive = Nrecvactive + 1
!                recvIDactive(Nrecvactive) = pp
!             endif
!          enddo
!          enddo

!          do ii = 1, kerls%nr
!          do jj = 1, kerls%nc
!             index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
!             index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
!             call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
!             pid = ptree%pgrp(pgno_sub)%head
!             pp = pid - ptree%pgrp(blocks%pgno)%head + 1
!             if (sendquant(pp)%active == 0) then
!                sendquant(pp)%active = 1
!                Nsendactive = Nsendactive + 1
!                sendIDactive(Nsendactive) = pp
!             endif
!             sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls%blocks(ii, jj)%matrix, 1)*size(kerls%blocks(ii, jj)%matrix, 2)
!          enddo
!          enddo
!       endif
!       n2 = MPI_Wtime()
!       time_tmp = time_tmp + n2 - n1

!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
!       enddo
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          sendquant(pp)%size = 0
!       enddo

!       Nreqr = 0
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          if (sendid /= ptree%MyID) then
!             Nreqr = Nreqr + 1
!          endif
!       enddo

!       n1 = MPI_Wtime()
!       ! pack the send buffer in the second pass
!       if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
!       do ii = 1, kerls%nr
!       do jj = 1, kerls%nc
!          index_i = (ii - 1)*kerls%inc_r + kerls%idx_r
!          index_j = (jj - 1)*kerls%inc_c + kerls%idx_c
!          call GetBlockPID(ptree, blocks%pgno, level_new, level_butterfly, index_i, index_j, mode_new, pgno_sub)
!          pid = ptree%pgrp(pgno_sub)%head

!          pp = pid - ptree%pgrp(blocks%pgno)%head + 1
!          Nrow = size(kerls%blocks(ii, jj)%matrix, 1)
!          Ncol = size(kerls%blocks(ii, jj)%matrix, 2)

!          sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
!          sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
!          sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
!          sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
!          sendquant(pp)%size = sendquant(pp)%size + 4
!          do i = 1, Nrow*Ncol
!             rr = mod(i - 1, Nrow) + 1
!             cc = (i - 1)/Nrow + 1
!             sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls%blocks(ii, jj)%matrix(rr, cc)
!          enddo
!          sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
!          deallocate (kerls%blocks(ii, jj)%matrix)
!       enddo
!       enddo
!       deallocate (kerls%blocks)

!       kerls%idx_r = idx_r
!       kerls%idx_c = idx_c
!       kerls%inc_r = inc_r
!       kerls%inc_c = inc_c
!       kerls%nr = nr
!       kerls%nc = nc

!       allocate (kerls%blocks(kerls%nr, kerls%nc))
!       endif
!       n2 = MPI_Wtime()
!       time_tmp = time_tmp + n2 - n1

!       Nreqs = 0
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
!          if (recvid /= ptree%MyID) then
!             Nreqs = Nreqs + 1
!             call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
!          endif
!       enddo

!       cnt = 0
!       do tt = 1, Nrecvactive
!          if(tt==1 .and. Nreqr+1== Nrecvactive)then
!             pp = ptree%MyID + 1 - ptree%pgrp(blocks%pgno)%head
!             recvquant(pp)%size=sendquant(pp)%size
!             if(recvquant(pp)%size>0)then
!             allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
!             recvquant(pp)%dat = sendquant(pp)%dat
!             endif
!          else
!             call MPI_Probe(MPI_ANY_SOURCE, tag+1, ptree%pgrp(blocks%pgno)%Comm, statusr(:,1),ierr)
!             pp = statusr(MPI_SOURCE, 1) + 1
!             call MPI_Get_count(statusr(:,1), MPI_DT, recvquant(pp)%size,ierr)
!             allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
!             cnt = cnt + 1
!             call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(cnt), ierr)
!          endif
!       enddo

!       if (Nreqr > 0) then
!          call MPI_waitall(Nreqr, R_req, statusr, ierr)
!       endif

!       ! copy data from buffer to target
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          n1 = MPI_Wtime()
!          i = 0
!          do while (i < recvquant(pp)%size)
!             i = i + 1
!             index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
!             ii = (index_i - kerls%idx_r)/kerls%inc_r + 1
!             i = i + 1
!             index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
!             jj = (index_j - kerls%idx_c)/kerls%inc_c + 1
!             i = i + 1
!             Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
!             i = i + 1
!             Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
!             call assert(.not. associated(kerls%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
!             allocate (kerls%blocks(ii, jj)%matrix(Nrow, Ncol))
!             do cc = 1, Ncol
!                do rr = 1, Nrow
!                   i = i + 1
!                   kerls%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
!                enddo
!             enddo
!          enddo
!          n2 = MPI_Wtime()
!          time_tmp = time_tmp + n2 - n1
!       enddo

!       if (Nreqs > 0) then
!          call MPI_waitall(Nreqs, S_req, statuss, ierr)
!       endif

!       ! deallocation
!       do tt = 1, Nsendactive
!          pp = sendIDactive(tt)
!          if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
!       enddo
!       do tt = 1, Nrecvactive
!          pp = recvIDactive(tt)
!          if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
!       enddo

!       ! n2 = MPI_Wtime()
!       ! time_tmp = time_tmp + n2 - n1

!    end subroutine BF_all2all_vec_n_ker

!>*********** all to all communication of one level of a butterfly from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
   subroutine BF_all2all_ker(block_i, pgno_i, kerls_i, level_i, offset_r, offset_c, block_o, pgno_o, kerls_o, level_o, stats, ptree)


      implicit none
      integer pgno_i, pgno_o, pgno, level_i, level_o
      integer i, j, level_butterfly_i, level_butterfly_o, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pgno_sub, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, num_row, num_col, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls_i, kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist
      character::mode
      integer offset_r, offset_c

      n1 = MPI_Wtime()

      nproc = max(ptree%pgrp(pgno_i)%nproc, ptree%pgrp(pgno_o)%nproc)
      pgno = min(pgno_i, pgno_o)
      tag = pgno+level_i*10

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_i = block_i%level_butterfly
      else
         level_butterfly_i = -1
         block_i%level_half = -1
      endif
      if (IOwnPgrp(ptree, pgno_o)) then
         level_butterfly_o = block_o%level_butterfly
      else
         level_butterfly_o = -1
         block_o%level_half = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_i, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_o%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      if (level_i <= block_i%level_half) then
         call assert(level_o <= block_o%level_half, 'row-wise ordering is only redistributed to row-wise ordering')
         mode = 'R'
      endif

      if (level_i > block_i%level_half) then
         call assert(level_o > block_o%level_half, 'column-wise ordering is only redistributed to column-wise ordering')
         mode = 'C'
      endif

      call assert((ptree%pgrp(pgno_i)%head <= ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail >= ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head <= ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail >= ptree%pgrp(pgno_i)%tail), 'pgno_i or pgno_o should be contained in the other')

      call GetLocalBlockRange(ptree, pgno_o, level_o, level_butterfly_o, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

      num_row = 2**level_o
      num_col = 2**(level_butterfly_o - level_o + 1)

      if (mode == 'R') then
         idx_c = idx_c*2 - 1
         inc_c = inc_c
         nc = nc*2
      elseif (mode == 'C') then
         idx_r = idx_r*2 - 1
         inc_r = inc_r
         nr = nr*2
      endif

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
         index_i = (ii - 1)*inc_r + idx_r - offset_r
         index_j = (jj - 1)*inc_c + idx_c - offset_c
         if (mode == 'R') then
            index_j0 = floor_safe((index_j - 1)/2d0) + 1
            index_i0 = index_i
         endif
         if (mode == 'C') then
            index_i0 = floor_safe((index_i - 1)/2d0) + 1
            index_j0 = index_j
         endif
         call GetBlockPID(ptree, pgno_i, level_i, level_butterfly_i, index_i0, index_j0, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         endif
      enddo
      enddo

      do ii = 1, kerls_i%nr
      do jj = 1, kerls_i%nc
         index_i = (ii - 1)*kerls_i%inc_r + kerls_i%idx_r + offset_r
         index_j = (jj - 1)*kerls_i%inc_c + kerls_i%idx_c + offset_c
         if (mode == 'R') then
            index_j0 = floor_safe((index_j - 1)/2d0) + 1
            index_i0 = index_i
         endif
         if (mode == 'C') then
            index_i0 = floor_safe((index_i - 1)/2d0) + 1
            index_j0 = index_j
         endif
         call GetBlockPID(ptree, pgno_o, level_o, level_butterfly_o, index_i0, index_j0, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls_i%blocks(ii, jj)%matrix, 1)*size(kerls_i%blocks(ii, jj)%matrix, 2)
         endif
      enddo
      enddo

      ! communicate receive buffer sizes
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if(recvid /= ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do ii = 1, kerls_i%nr
      do jj = 1, kerls_i%nc
         index_i = (ii - 1)*kerls_i%inc_r + kerls_i%idx_r + offset_r
         index_j = (jj - 1)*kerls_i%inc_c + kerls_i%idx_c + offset_c

         if (mode == 'R') then
            index_j0 = floor_safe((index_j - 1)/2d0) + 1
            index_i0 = index_i
         endif
         if (mode == 'C') then
            index_i0 = floor_safe((index_i - 1)/2d0) + 1
            index_j0 = index_j
         endif

         call GetBlockPID(ptree, pgno_o, level_o, level_butterfly_o, index_i0, index_j0, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            Nrow = size(kerls_i%blocks(ii, jj)%matrix, 1)
            Ncol = size(kerls_i%blocks(ii, jj)%matrix, 2)

            sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
            sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
            sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
            sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
            sendquant(pp)%size = sendquant(pp)%size + 4
            do i = 1, Nrow*Ncol
               rr = mod(i - 1, Nrow) + 1
               cc = (i - 1)/Nrow + 1
               sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls_i%blocks(ii, jj)%matrix(rr, cc)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
            deallocate (kerls_i%blocks(ii, jj)%matrix)
         endif
      enddo
      enddo
      if (allocated(kerls_i%blocks)) deallocate (kerls_i%blocks)

      if (nr > 0 .and. nc > 0) then
         kerls_o%idx_r = idx_r
         kerls_o%idx_c = idx_c
         kerls_o%inc_r = inc_r
         kerls_o%inc_c = inc_c
         kerls_o%nr = nr
         kerls_o%nc = nc
         kerls_o%num_row = num_row
         kerls_o%num_col = num_col
         if (.not. allocated(kerls_o%blocks)) allocate (kerls_o%blocks(kerls_o%nr, kerls_o%nc))
      endif


      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo


      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - kerls_o%idx_r)/kerls_o%inc_r + 1
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jj = (index_j - kerls_o%idx_c)/kerls_o%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(kerls_o%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
            allocate (kerls_o%blocks(ii, jj)%matrix(Nrow, Ncol))
            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  kerls_o%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
               enddo
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_ker









!>*********** all to all communication of one level of a butterfly from an old pattern pat_i to an new pattern pat_o
   subroutine BF_all2all_ker_pattern(block_i, kerls_i, pat_i, block_o,kerls_o, pat_o, level, pgno,stats, ptree)


      implicit none
      integer pgno
      integer level_butterfly
      integer index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt, i, j
      integer level
      integer index_ii, index_jj, index_ii_loc, index_jj_loc
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree
      integer ierr, nsendrecv, pgno_sub, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, num_row, num_col, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls_i, kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist,level_half_i,level_half_o,pat_i,pat_o
      character::mode_i,mode_o

      level_butterfly = block_i%level_butterfly

      level_half_i = BF_Switchlevel(level_butterfly, pat_i)
      if(level<=level_half_i)then
         mode_i='R'
      else
         mode_i='C'
      endif
      level_half_o = BF_Switchlevel(level_butterfly, pat_o)
      if(level<=level_half_o)then
         mode_o='R'
      else
         mode_o='C'
      endif

      if(mode_i/=mode_o)then

         n1 = MPI_Wtime()

         nproc = ptree%pgrp(pgno)%nproc
         tag = pgno+level*10


         call GetLocalBlockRange(ptree, pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, mode_o)

         num_row = 2**level
         num_col = 2**(level_butterfly - level + 1)

         if (mode_o == 'R') then
            idx_c = idx_c*2 - 1
            inc_c = inc_c
            nc = nc*2
         elseif (mode_o == 'C') then
            idx_r = idx_r*2 - 1
            inc_r = inc_r
            nr = nr*2
         endif

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
            if (mode_i == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
               index_i0 = index_i
            endif
            if (mode_i == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = index_j
            endif
            call GetBlockPID(ptree, pgno, level, level_butterfly, index_i0, index_j0, mode_i, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            if (pid /= -1) then
               pp = pid - ptree%pgrp(pgno)%head + 1
               if (recvquant(pp)%active == 0) then
                  recvquant(pp)%active = 1
                  Nrecvactive = Nrecvactive + 1
                  recvIDactive(Nrecvactive) = pp
               endif
            endif
         enddo
         enddo

         do ii = 1, kerls_i%nr
         do jj = 1, kerls_i%nc
            index_i = (ii - 1)*kerls_i%inc_r + kerls_i%idx_r
            index_j = (jj - 1)*kerls_i%inc_c + kerls_i%idx_c
            if (mode_o == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
               index_i0 = index_i
            endif
            if (mode_o == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = index_j
            endif
            call GetBlockPID(ptree, pgno, level, level_butterfly, index_i0, index_j0, mode_o, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            if (pid /= -1) then
               pp = pid - ptree%pgrp(pgno)%head + 1
               if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
               endif
               sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls_i%blocks(ii, jj)%matrix, 1)*size(kerls_i%blocks(ii, jj)%matrix, 2)
            endif
         enddo
         enddo

         ! communicate receive buffer sizes
         Nreqs = 0
         do tt = 1, Nsendactive
            pp = sendIDactive(tt)
            allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
            recvid = pp - 1 + ptree%pgrp(pgno)%head
            if(recvid /=ptree%MyID)then
               Nreqs = Nreqs + 1
               call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
            else
               recvquant(pp)%size = sendquant(pp)%size
            endif
         enddo

         Nreqr = 0
         do tt = 1, Nrecvactive
            pp = recvIDactive(tt)
            sendid = pp - 1 + ptree%pgrp(pgno)%head
            if(sendid /=ptree%MyID)then
               Nreqr = Nreqr + 1
               call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
            sendquant(pp)%size = 0
         enddo
         do tt = 1, Nrecvactive
            pp = recvIDactive(tt)
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
         enddo

         ! pack the send buffer in the second pass
         do ii = 1, kerls_i%nr
         do jj = 1, kerls_i%nc
            index_i = (ii - 1)*kerls_i%inc_r + kerls_i%idx_r
            index_j = (jj - 1)*kerls_i%inc_c + kerls_i%idx_c

            if (mode_o == 'R') then
               index_j0 = floor_safe((index_j - 1)/2d0) + 1
               index_i0 = index_i
            endif
            if (mode_o == 'C') then
               index_i0 = floor_safe((index_i - 1)/2d0) + 1
               index_j0 = index_j
            endif

            call GetBlockPID(ptree, pgno, level, level_butterfly, index_i0, index_j0, mode_o, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            if (pid /= -1) then
               pp = pid - ptree%pgrp(pgno)%head + 1
               Nrow = size(kerls_i%blocks(ii, jj)%matrix, 1)
               Ncol = size(kerls_i%blocks(ii, jj)%matrix, 2)

               sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
               sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
               sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
               sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
               sendquant(pp)%size = sendquant(pp)%size + 4
               do i = 1, Nrow*Ncol
                  rr = mod(i - 1, Nrow) + 1
                  cc = (i - 1)/Nrow + 1
                  sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls_i%blocks(ii, jj)%matrix(rr, cc)
               enddo
               sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
               deallocate (kerls_i%blocks(ii, jj)%matrix)
            endif
         enddo
         enddo
         if (allocated(kerls_i%blocks)) deallocate (kerls_i%blocks)

         if (nr > 0 .and. nc > 0) then
            kerls_o%idx_r = idx_r
            kerls_o%idx_c = idx_c
            kerls_o%inc_r = inc_r
            kerls_o%inc_c = inc_c
            kerls_o%nr = nr
            kerls_o%nc = nc
            kerls_o%num_row = num_row
            kerls_o%num_col = num_col
            if (.not. allocated(kerls_o%blocks)) allocate (kerls_o%blocks(kerls_o%nr, kerls_o%nc))
         endif

         Nreqs = 0
         do tt = 1, Nsendactive
            pp = sendIDactive(tt)
            recvid = pp - 1 + ptree%pgrp(pgno)%head
            if (recvid /= ptree%MyID) then
               Nreqs = Nreqs + 1
               call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
            else
               if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
            endif
         enddo

         Nreqr = 0
         do tt = 1, Nrecvactive
            pp = recvIDactive(tt)
            sendid = pp - 1 + ptree%pgrp(pgno)%head
            if (sendid /= ptree%MyID) then
               Nreqr = Nreqr + 1
               call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
            endif
         enddo

         ! copy data from buffer to target
         do tt = 1, Nrecvactive
            if(tt==1 .and. Nreqr+1== Nrecvactive)then
               pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
            else
               call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
               pp = statusr(MPI_SOURCE, 1) + 1
            endif
            i = 0
            do while (i < recvquant(pp)%size)
               i = i + 1
               index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
               ii = (index_i - kerls_o%idx_r)/kerls_o%inc_r + 1
               i = i + 1
               index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
               jj = (index_j - kerls_o%idx_c)/kerls_o%inc_c + 1
               i = i + 1
               Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
               i = i + 1
               Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
               call assert(.not. associated(kerls_o%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
               allocate (kerls_o%blocks(ii, jj)%matrix(Nrow, Ncol))
               do cc = 1, Ncol
                  do rr = 1, Nrow
                     i = i + 1
                     kerls_o%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
                  enddo
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
            if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
         enddo
         deallocate (sendquant)
         do tt = 1, Nrecvactive
            pp = recvIDactive(tt)
            if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
         enddo
         deallocate (recvquant)
         deallocate (sendIDactive)
         deallocate (recvIDactive)

         n2 = MPI_Wtime()
         ! time_tmp = time_tmp + n2 - n1
      endif
   end subroutine BF_all2all_ker_pattern




!>*********** convert blocks in block_i%sons to block_o%sons, this is a local function without MPI communication, it is assumed block_i%sons has L levels, and block_o%sons will have max(L-2,0) levels
   subroutine BF_convert_to_smallBF(block_i, block_o, stats, ptree)


      implicit none
      integer i, j, level_butterfly_i, level_butterfly_o, level_butterfly_c_o, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_ic, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_jc, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt, iii, jjj
      integer mm1, nn1, mm2, nn2, M1, N1, kk
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(matrixblock), pointer::block_c_i, block_c_o
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi, pgno_sub_mine
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, num_row, num_col, level_new, idx, level_butterfly_o_true

      real(kind=8)::t1, t2
      DT, allocatable::matrixtemp1(:, :), matrixtemp2(:, :)

      t1 = MPI_Wtime()

      do iii = 1, 2
         do jjj = 1, 2
         if (IOwnPgrp(ptree, block_o%sons(iii, jjj)%pgno)) then
            level_butterfly_o_true = max(block_i%level_butterfly - 2, 0)
            call GetPgno_Sub(ptree, block_o%sons(iii, jjj)%pgno, level_butterfly_o_true, pgno_sub_mine)
            if (ptree%pgrp(pgno_sub_mine)%head == ptree%MyID) then
               block_c_o => block_o%sons(iii, jjj)
               block_c_i => block_i%sons(iii, jjj)
               if (block_i%level_butterfly == 1) then  ! l-level butterfly becomes 0-level butterfly
                  block_c_o%level_butterfly = 0
                  block_c_o%level_half = 0
                  allocate (block_c_o%ButterflyU%blocks(1))
                  allocate (block_c_o%ButterflyV%blocks(1))
                  block_c_o%ButterflyU%nblk_loc = 1
                  block_c_o%ButterflyU%inc = 1
                  block_c_o%ButterflyU%idx = 1
                  block_c_o%ButterflyV%nblk_loc = 1
                  block_c_o%ButterflyV%inc = 1
                  block_c_o%ButterflyV%idx = 1

                  call assert(block_c_i%ButterflyU%nblk_loc == 1, 'parent block has more than one ButterflyU block')
                  call assert(block_c_i%ButterflyV%nblk_loc == 1, 'parent block has more than one ButterflyV block')
                  call assert(block_c_i%ButterflyKerl(1)%nr == 1 .and. block_c_i%ButterflyKerl(1)%nc == 1, 'parent block has more than one ButterflyKerl block')

                  mm1 = size(block_c_i%ButterflyKerl(1)%blocks(1, 1)%matrix, 1)
                  nn1 = size(block_c_i%ButterflyKerl(1)%blocks(1, 1)%matrix, 2)
                  M1 = size(block_c_i%ButterflyU%blocks(1)%matrix, 1)
                  N1 = size(block_c_i%ButterflyV%blocks(1)%matrix, 1)

                  allocate (block_c_o%ButterflyU%blocks(1)%matrix(M1, nn1))
                  allocate (block_c_o%ButterflyV%blocks(1)%matrix(N1, nn1))
                  call gemmf90(block_c_i%ButterflyU%blocks(1)%matrix, M1, block_c_i%ButterflyKerl(1)%blocks(1, 1)%matrix, mm1, block_c_o%ButterflyU%blocks(1)%matrix, M1, 'N', 'N', M1, nn1, mm1, BPACK_cone, BPACK_czero)
                  block_c_o%ButterflyV%blocks(1)%matrix = block_c_i%ButterflyV%blocks(1)%matrix
               else ! L-level butterfly becomes (L-2)-level butterfly
                  block_c_o%level_butterfly = block_i%level_butterfly - 2
                  block_c_o%level_half = floor_safe(dble(block_c_o%level_butterfly)/2d0) ! from outer to inner
                  call assert(block_c_o%level_butterfly >= 0, 'negative level_butterfly!')
                  if (block_c_o%level_butterfly > 0) then
                     allocate (block_c_o%ButterflyKerl(block_c_o%level_butterfly))
                  endif
                  do level = 0, block_c_o%level_butterfly + 1
                     if (level == 0) then
                        block_c_o%ButterflyV%num_blk = 2**block_c_o%level_butterfly
                        call assert(mod(block_c_i%ButterflyV%nblk_loc, 2) == 0, 'parent block should have even number of ButterflyV blocks')
                        block_c_o%ButterflyV%nblk_loc = block_c_i%ButterflyV%nblk_loc/2
                        block_c_o%ButterflyV%inc = block_c_i%ButterflyV%inc
                        idx = block_c_i%ButterflyV%idx
                        if (idx > block_c_o%ButterflyV%num_blk*2) idx = idx - block_c_o%ButterflyV%num_blk*2
                        block_c_o%ButterflyV%idx = ceiling_safe(idx/2d0)
                        allocate (block_c_o%ButterflyV%blocks(block_c_o%ButterflyV%nblk_loc))
                        do ii = 1, block_c_o%ButterflyV%nblk_loc
                           mm1 = size(block_c_i%ButterflyV%blocks(2*ii - 1)%matrix, 1)
                           nn1 = size(block_c_i%ButterflyV%blocks(2*ii - 1)%matrix, 2)
                           mm2 = size(block_c_i%ButterflyV%blocks(2*ii)%matrix, 1)
                           nn2 = size(block_c_i%ButterflyV%blocks(2*ii)%matrix, 2)
                           kk = size(block_c_i%ButterflyKerl(1)%blocks(1, 2*ii - 1)%matrix, 1)
                           allocate (block_c_o%ButterflyV%blocks(ii)%matrix(mm1 + mm2, kk))
                           N1 = N1 + mm1 + mm2
                           allocate (matrixtemp1(mm1, kk))
                           matrixtemp1 = 0
                           allocate (matrixtemp2(mm2, kk))
                           matrixtemp2 = 0
                           call gemmf90(block_c_i%ButterflyV%blocks(2*ii - 1)%matrix, mm1, block_c_i%ButterflyKerl(1)%blocks(1, 2*ii - 1)%matrix, kk, matrixtemp1, mm1, 'N', 'T', mm1, kk, nn1, BPACK_cone, BPACK_czero)
                           call gemmf90(block_c_i%ButterflyV%blocks(2*ii)%matrix, mm2, block_c_i%ButterflyKerl(1)%blocks(1, 2*ii)%matrix, kk, matrixtemp2, mm2, 'N', 'T', mm2, kk, nn2, BPACK_cone, BPACK_czero)
                           block_c_o%ButterflyV%blocks(ii)%matrix(1:mm1, 1:kk) = matrixtemp1
                           block_c_o%ButterflyV%blocks(ii)%matrix(1 + mm1:mm1 + mm2, 1:kk) = matrixtemp2
                           deallocate (matrixtemp1)
                           deallocate (matrixtemp2)
                        enddo
                     elseif (level == block_c_o%level_butterfly + 1) then
                        block_c_o%ButterflyU%num_blk = 2**block_c_o%level_butterfly
                        call assert(mod(block_c_i%ButterflyU%nblk_loc, 2) == 0, 'parent block should have even number of ButterflyU blocks')
                        block_c_o%ButterflyU%nblk_loc = block_c_i%ButterflyU%nblk_loc/2
                        block_c_o%ButterflyU%inc = block_c_i%ButterflyU%inc
                        idx = block_c_i%ButterflyU%idx
                        if (idx > block_c_o%ButterflyU%num_blk*2) idx = idx - block_c_o%ButterflyU%num_blk*2
                        block_c_o%ButterflyU%idx = ceiling_safe(idx/2d0)
                        allocate (block_c_o%ButterflyU%blocks(block_c_o%ButterflyU%nblk_loc))
                        do ii = 1, block_c_o%ButterflyU%nblk_loc
                           mm1 = size(block_c_i%ButterflyU%blocks(2*ii - 1)%matrix, 1)
                           nn1 = size(block_c_i%ButterflyU%blocks(2*ii - 1)%matrix, 2)
                           mm2 = size(block_c_i%ButterflyU%blocks(2*ii)%matrix, 1)
                           nn2 = size(block_c_i%ButterflyU%blocks(2*ii)%matrix, 2)
                           kk = size(block_c_i%ButterflyKerl(block_c_o%level_butterfly + 2)%blocks(2*ii - 1, 1)%matrix, 2)
                           allocate (block_c_o%ButterflyU%blocks(ii)%matrix(mm1 + mm2, kk))
                           M1 = M1 + mm1 + mm2
                           allocate (matrixtemp1(mm1, kk))
                           matrixtemp1 = 0
                           allocate (matrixtemp2(mm2, kk))
                           matrixtemp2 = 0
                           call gemmf90(block_c_i%ButterflyU%blocks(2*ii - 1)%matrix, mm1, block_c_i%ButterflyKerl(block_c_i%level_butterfly)%blocks(2*ii - 1, 1)%matrix, nn1, matrixtemp1, mm1, 'N', 'N', mm1, kk, nn1, BPACK_cone, BPACK_czero)
                           call gemmf90(block_c_i%ButterflyU%blocks(2*ii)%matrix, mm2, block_c_i%ButterflyKerl(block_c_i%level_butterfly)%blocks(2*ii, 1)%matrix, nn2, matrixtemp2, mm2, 'N', 'N', mm2, kk, nn2, BPACK_cone, BPACK_czero)
                           block_c_o%ButterflyU%blocks(ii)%matrix(1:mm1, 1:kk) = matrixtemp1
                           block_c_o%ButterflyU%blocks(ii)%matrix(1 + mm1:mm1 + mm2, 1:kk) = matrixtemp2
                           deallocate (matrixtemp1)
                           deallocate (matrixtemp2)
                        end do
                     else
                        num_col = block_c_i%ButterflyKerl(level + 1)%num_col
                        num_row = block_c_i%ButterflyKerl(level + 1)%num_row
                        block_c_o%ButterflyKerl(level)%num_row = num_row/2
                        block_c_o%ButterflyKerl(level)%num_col = num_col/2
                        block_c_o%ButterflyKerl(level)%nr = block_c_i%ButterflyKerl(level + 1)%nr
                        block_c_o%ButterflyKerl(level)%inc_r = block_c_i%ButterflyKerl(level + 1)%inc_r
                        block_c_o%ButterflyKerl(level)%idx_r = block_c_i%ButterflyKerl(level + 1)%idx_r
                        if (block_c_o%ButterflyKerl(level)%idx_r > block_c_o%ButterflyKerl(level)%num_row) block_c_o%ButterflyKerl(level)%idx_r = block_c_o%ButterflyKerl(level)%idx_r - block_c_o%ButterflyKerl(level)%num_row
                        block_c_o%ButterflyKerl(level)%nc = block_c_i%ButterflyKerl(level + 1)%nc
                        block_c_o%ButterflyKerl(level)%inc_c = block_c_i%ButterflyKerl(level + 1)%inc_c
                        block_c_o%ButterflyKerl(level)%idx_c = block_c_i%ButterflyKerl(level + 1)%idx_c
                        if (block_c_o%ButterflyKerl(level)%idx_c > block_c_o%ButterflyKerl(level)%num_col) block_c_o%ButterflyKerl(level)%idx_c = block_c_o%ButterflyKerl(level)%idx_c - block_c_o%ButterflyKerl(level)%num_col

                        if(block_c_o%ButterflyKerl(level)%nr>0 .and. block_c_o%ButterflyKerl(level)%nc>0)then
                           allocate (block_c_o%ButterflyKerl(level)%blocks(block_c_o%ButterflyKerl(level)%nr, block_c_o%ButterflyKerl(level)%nc))

                           do ii = 1, block_c_o%ButterflyKerl(level)%nr
                           do jj = 1, block_c_o%ButterflyKerl(level)%nc
                              mm = size(block_c_i%ButterflyKerl(level + 1)%blocks(ii, jj)%matrix, 1)
                              nn = size(block_c_i%ButterflyKerl(level + 1)%blocks(ii, jj)%matrix, 2)
                              allocate (block_c_o%ButterflyKerl(level)%blocks(ii, jj)%matrix(mm, nn))
                              block_c_o%ButterflyKerl(level)%blocks(ii, jj)%matrix = block_c_i%ButterflyKerl(level + 1)%blocks(ii, jj)%matrix
                           enddo
                           enddo
                        endif
                     endif
                  enddo
               endif
            endif
         endif
         enddo
      enddo
      t2 = MPI_Wtime()
      ! time_tmp = time_tmp + t2 - t1
   end subroutine BF_convert_to_smallBF

!>*********** all to all communication of one level of a butterfly into four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
   subroutine BF_all2all_ker_split(block_i, pgno_i, level_i, block_o, pgnos_o, level_o, stats, ptree)


      implicit none
      integer pgno_sub, pgno_i, pgno_o, pgno, level_i, level_o
      integer i, j, level_butterfly_i, level_butterfly_o, level_butterfly_c_o, level_butterfly_o_true, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_ic, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_jc, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt, iii, jjj
      integer pgnos_o(2,2)
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1, pgno_sub_i_mine, pgno_sub_o_mine
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, num_row, num_col, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      ! type(butterfly_kerl)::kerls_i,kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist
      character::mode

      n1 = MPI_Wtime()

      nproc = ptree%pgrp(pgno_i)%nproc
      pgno = pgno_i
      do iii=1,2
      do jjj=1,2
      nproc = max(nproc, ptree%pgrp(pgnos_o(iii,jjj))%nproc)
      pgno = min(pgno_i, pgnos_o(iii,jjj))
      enddo
      enddo
      tag = pgno+level_i*10

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_i = block_i%level_butterfly
      else
         level_butterfly_i = -1
         block_i%level_half = -1
      endif
      if (IOwnPgrp(ptree, pgnos_o(1,1))) then
         level_butterfly_o = block_o%level_butterfly
      else
         level_butterfly_o = -1
         block_o%level_half = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_i, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_o%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      level_butterfly_o_true = max(level_butterfly_i - 2, 0)

      if (level_i <= block_i%level_half) then
         call assert(level_o <= block_o%level_half, 'row-wise ordering is only redistributed to row-wise ordering')
         mode = 'R'
      endif

      if (level_i > block_i%level_half) then
         call assert(level_o > block_o%level_half, 'column-wise ordering is only redistributed to column-wise ordering')
         mode = 'C'
      endif

      do iii=1,2
         do jjj=1,2
             call assert((ptree%pgrp(pgno_i)%head <= ptree%pgrp(pgnos_o(iii,jjj))%head .and. ptree%pgrp(pgno_i)%tail >= ptree%pgrp(pgnos_o(iii,jjj))%tail) .or. (ptree%pgrp(pgnos_o(iii,jjj))%head <= ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgnos_o(iii,jjj))%tail >= ptree%pgrp(pgno_i)%tail), 'pgno_i or pgno_o should be contained in the other')
         enddo
      enddo
      call GetPgno_Sub(ptree, pgno_i, level_butterfly_i, pgno_sub_i_mine)

      num_row = 2**level_o
      num_col = 2**(level_butterfly_o - level_o + 1)
      level_butterfly_c_o = max(level_butterfly_o - 2, 0)

      do iii = 1, 2
         do jjj = 1, 2
             pgno_o=pgnos_o(iii,jjj)
             if (IOwnPgrp(ptree, pgno_o)) then
               call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
               call GetLocalBlockRange(ptree, pgno_o, level_o - 1, level_butterfly_c_o, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)
               if (mode == 'R') then
                  idx_c = idx_c*2 - 1
                  inc_c = inc_c
                  if (level_butterfly_o > 1) nc = nc*2
               elseif (mode == 'C') then
                  idx_r = idx_r*2 - 1
                  inc_r = inc_r
                  if (level_butterfly_o > 1) nr = nr*2
               endif
               if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
                  if (nr > 0 .and. nc > 0) then
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_r = idx_r + (iii - 1)*num_row/2
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_c = idx_c + (jjj - 1)*num_col/2
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_r = inc_r
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_c = inc_c
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nr = nr
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nc = nc
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%num_row = num_row
                     block_o%sons(iii, jjj)%ButterflyKerl(level_o)%num_col = num_col
                     allocate (block_o%sons(iii, jjj)%ButterflyKerl(level_o)%blocks(block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nr, block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nc))
                  endif
               endif
            endif
         enddo
      enddo

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
      do iii = 1, 2
         do jjj = 1, 2
             pgno_o=pgnos_o(iii,jjj)
             if (IOwnPgrp(ptree, pgno_o)) then
               call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
               if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
                  do ii = 1, block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nr
                  do jj = 1, block_o%sons(iii, jjj)%ButterflyKerl(level_o)%nc
                     index_i = (ii - 1)*block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_r + block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_r
                     index_j = (jj - 1)*block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_c + block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_c
                     if (mode == 'R') then
                        index_j0 = floor_safe((index_j - 1)/2d0) + 1
                        index_i0 = index_i
                     endif
                     if (mode == 'C') then
                        index_i0 = floor_safe((index_i - 1)/2d0) + 1
                        index_j0 = index_j
                     endif
                     call GetBlockPID(ptree, pgno_i, level_i, level_butterfly_i, index_i0, index_j0, mode, pgno_sub)
                     pid = ptree%pgrp(pgno_sub)%head
                     pp = pid - ptree%pgrp(pgno)%head + 1
                     if (recvquant(pp)%active == 0) then
                        recvquant(pp)%active = 1
                        Nrecvactive = Nrecvactive + 1
                        recvIDactive(Nrecvactive) = pp
                     endif
                  enddo
                  enddo
               endif
            endif
         enddo
      enddo


      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then

         do ii = 1, block_i%ButterflyKerl(level_i)%nr
         do jj = 1, block_i%ButterflyKerl(level_i)%nc
            index_i = (ii - 1)*block_i%ButterflyKerl(level_i)%inc_r + block_i%ButterflyKerl(level_i)%idx_r
            index_j = (jj - 1)*block_i%ButterflyKerl(level_i)%inc_c + block_i%ButterflyKerl(level_i)%idx_c
            index_ic = index_i
            index_jc = index_j

            iii = 1
            jjj = 1
            if (index_ic > num_row/2)then
                index_ic = index_ic - num_row/2
                iii=2
            endif
            if (index_jc > num_col/2)then
               index_jc = index_jc - num_col/2
               jjj=2
            endif
            pgno_o=pgnos_o(iii,jjj)

            if (mode == 'R') then
               index_j0 = floor_safe((index_jc - 1)/2d0) + 1
               index_i0 = index_ic
            endif
            if (mode == 'C') then
               index_i0 = floor_safe((index_ic - 1)/2d0) + 1
               index_j0 = index_jc
            endif
            call GetBlockPID(ptree, pgno_o, level_o - 1, level_butterfly_c_o, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head
            pp = pid - ptree%pgrp(pgno)%head + 1
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 4 + size(block_i%ButterflyKerl(level_i)%blocks(ii, jj)%matrix, 1)*size(block_i%ButterflyKerl(level_i)%blocks(ii, jj)%matrix, 2)
         enddo
         enddo
      endif
      endif

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if(recvid /=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr=0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if(sendid /=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then
         ! pack the send buffer in the second pass
         do ii = 1, block_i%ButterflyKerl(level_i)%nr
         do jj = 1, block_i%ButterflyKerl(level_i)%nc
            index_i = (ii - 1)*block_i%ButterflyKerl(level_i)%inc_r + block_i%ButterflyKerl(level_i)%idx_r
            index_j = (jj - 1)*block_i%ButterflyKerl(level_i)%inc_c + block_i%ButterflyKerl(level_i)%idx_c

            index_ic = index_i
            index_jc = index_j
            iii = 1
            jjj = 1
            if (index_ic > num_row/2)then
                index_ic = index_ic - num_row/2
                iii=2
            endif
            if (index_jc > num_col/2)then
               index_jc = index_jc - num_col/2
               jjj=2
            endif
            pgno_o=pgnos_o(iii,jjj)

            if (mode == 'R') then
               index_j0 = floor_safe((index_jc - 1)/2d0) + 1
               index_i0 = index_ic
            endif
            if (mode == 'C') then
               index_i0 = floor_safe((index_ic - 1)/2d0) + 1
               index_j0 = index_jc
            endif
            call GetBlockPID(ptree, pgno_o, level_o - 1, level_butterfly_c_o, index_i0, index_j0, mode, pgno_sub)
            pid = ptree%pgrp(pgno_sub)%head

            pp = pid - ptree%pgrp(pgno)%head + 1
            Nrow = size(block_i%ButterflyKerl(level_i)%blocks(ii, jj)%matrix, 1)
            Ncol = size(block_i%ButterflyKerl(level_i)%blocks(ii, jj)%matrix, 2)

            sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
            sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
            sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
            sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
            sendquant(pp)%size = sendquant(pp)%size + 4
            do i = 1, Nrow*Ncol
               rr = mod(i - 1, Nrow) + 1
               cc = (i - 1)/Nrow + 1
               sendquant(pp)%dat(sendquant(pp)%size + i, 1) = block_i%ButterflyKerl(level_i)%blocks(ii, jj)%matrix(rr, cc)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
            ! deallocate(block_i%ButterflyKerl(level_i)%blocks(ii,jj)%matrix)
         enddo
         enddo
      endif
      endif
      ! if(allocated(block_i%ButterflyKerl(level_i)%blocks))deallocate(block_i%ButterflyKerl(level_i)%blocks)

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            iii = 1
            if (index_i > num_row/2) iii = 2
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            jjj = 1
            if (index_j > num_col/2) jjj = 2

            ii = (index_i - block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_r)/block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_r + 1
            jj = (index_j - block_o%sons(iii, jjj)%ButterflyKerl(level_o)%idx_c)/block_o%sons(iii, jjj)%ButterflyKerl(level_o)%inc_c + 1
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(block_o%sons(iii, jjj)%ButterflyKerl(level_o)%blocks(ii, jj)%matrix), 'receiving dat alreay exists locally')
            allocate (block_o%sons(iii, jjj)%ButterflyKerl(level_o)%blocks(ii, jj)%matrix(Nrow, Ncol))
            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  block_o%sons(iii, jjj)%ButterflyKerl(level_o)%blocks(ii, jj)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
               enddo
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_ker_split

!>*********** all to all communication of one level of a butterfly from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
   subroutine BF_all2all_UV(block_i, pgno_i, kerls_i, level_i, offset, block_o, pgno_o, kerls_o, level_o, stats, ptree)


      implicit none
      integer pgno_sub, pgno_i, pgno_o, pgno, level_i, level_o
      integer i, j, level_butterfly_i, level_butterfly_o, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree
      integer offset

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, idx, inc_r, inc_c, inc, nr, nc, nblk_loc, num_blk, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_UV)::kerls_i, kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist
      character::mode

      n1 = MPI_Wtime()

      nproc = max(ptree%pgrp(pgno_i)%nproc, ptree%pgrp(pgno_o)%nproc)
      pgno = min(pgno_i, pgno_o)
      tag = pgno+level_i*10

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_i = block_i%level_butterfly
      else
         level_butterfly_i = -1
         block_i%level_half = -1
      endif
      if (IOwnPgrp(ptree, pgno_o)) then
         level_butterfly_o = block_o%level_butterfly
      else
         level_butterfly_o = -1
         block_o%level_half = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_i, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_o%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      if (level_i <= block_i%level_half) then
         call assert(level_o <= block_o%level_half, 'row-wise ordering is only redistributed to row-wise ordering')
         mode = 'R'
      endif

      if (level_i > block_i%level_half) then
         call assert(level_o > block_o%level_half, 'column-wise ordering is only redistributed to column-wise ordering')
         mode = 'C'
      endif

      call assert((ptree%pgrp(pgno_i)%head <= ptree%pgrp(pgno_o)%head .and. ptree%pgrp(pgno_i)%tail >= ptree%pgrp(pgno_o)%tail) .or. (ptree%pgrp(pgno_o)%head <= ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgno_o)%tail >= ptree%pgrp(pgno_i)%tail), 'pgno_i or pgno_o should be contained in the other')

      call GetLocalBlockRange(ptree, pgno_o, level_o, level_butterfly_o, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

      if (mode == 'R') then
         num_blk = 2**level_butterfly_o
         idx = idx_c
         inc = inc_c
         nblk_loc = nc
      elseif (mode == 'C') then
         num_blk = 2**level_butterfly_o
         idx = idx_r
         inc = inc_r
         nblk_loc = nr
      endif

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
      do ii = 1, nblk_loc
         ! convert indices from output to input
         if (mode == 'R') then
            index_i = 1
            index_j = (ii - 1)*inc + idx - offset
         elseif (mode == 'C') then
            index_i = (ii - 1)*inc + idx - offset
            index_j = 1
         endif
         call GetBlockPID(ptree, pgno_i, level_i, level_butterfly_i, index_i, index_j, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         endif
      enddo

      do ii = 1, kerls_i%nblk_loc
         ! convert indices from input to output
         if (mode == 'R') then
            index_i = 1
            index_j = (ii - 1)*kerls_i%inc + kerls_i%idx + offset
         elseif (mode == 'C') then
            index_i = (ii - 1)*kerls_i%inc + kerls_i%idx + offset
            index_j = 1
         endif
         call GetBlockPID(ptree, pgno_o, level_o, level_butterfly_o, index_i, index_j, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 4 + size(kerls_i%blocks(ii)%matrix, 1)*size(kerls_i%blocks(ii)%matrix, 2)
         endif
      enddo

      ! communicate receive buffer sizes
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if(recvid /= ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do ii = 1, kerls_i%nblk_loc
         ! convert indices from input to output
         if (mode == 'R') then
            index_i = 1
            index_j = (ii - 1)*kerls_i%inc + kerls_i%idx + offset
         elseif (mode == 'C') then
            index_i = (ii - 1)*kerls_i%inc + kerls_i%idx + offset
            index_j = 1
         endif

         call GetBlockPID(ptree, pgno_o, level_o, level_butterfly_o, index_i, index_j, mode, pgno_sub)
         pid = ptree%pgrp(pgno_sub)%head
         if (pid /= -1) then
            pp = pid - ptree%pgrp(pgno)%head + 1
            Nrow = size(kerls_i%blocks(ii)%matrix, 1)
            Ncol = size(kerls_i%blocks(ii)%matrix, 2)

            sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
            sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
            sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
            sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
            sendquant(pp)%size = sendquant(pp)%size + 4
            do i = 1, Nrow*Ncol
               rr = mod(i - 1, Nrow) + 1
               cc = (i - 1)/Nrow + 1
               sendquant(pp)%dat(sendquant(pp)%size + i, 1) = kerls_i%blocks(ii)%matrix(rr, cc)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
            deallocate (kerls_i%blocks(ii)%matrix)
         endif
      enddo
      if (allocated(kerls_i%blocks)) deallocate (kerls_i%blocks)

      if (nblk_loc > 0) then
         kerls_o%idx = idx
         kerls_o%inc = inc
         kerls_o%nblk_loc = nblk_loc
         kerls_o%num_blk = num_blk
         if (.not. allocated(kerls_o%blocks)) allocate (kerls_o%blocks(kerls_o%nblk_loc))
      endif

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo


      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))

            if (mode == 'R') then
               ii = (index_j - kerls_o%idx)/kerls_o%inc + 1
            elseif (mode == 'C') then
               ii = (index_i - kerls_o%idx)/kerls_o%inc + 1
            endif
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(kerls_o%blocks(ii)%matrix), 'receiving dat alreay exists locally')
            allocate (kerls_o%blocks(ii)%matrix(Nrow, Ncol))
            do cc = 1, Ncol
               do rr = 1, Nrow
                  i = i + 1
                  kerls_o%blocks(ii)%matrix(rr, cc) = recvquant(pp)%dat(i, 1)
               enddo
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_UV

!>*********** all to all communication of one level of a butterfly to four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
   subroutine BF_all2all_U_split(block_i, pgno_i, level_i, block_o, pgnos_o, level_o, stats, ptree)


      implicit none
      integer pgno_sub, pgno_i, pgno_o, pgno, level_i, level_o, level_c_o
      integer i, j, level_butterfly_i, level_butterfly_o, level_butterfly_o_true, level_butterfly_c_o, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      integer pgnos_o(2,2)
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1, iii, jjj, si, sj, ni, nj, pgno_sub_i_mine, pgno_sub_o_mine
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, idx, inc_r, inc_c, inc, nr, nc, nblk_loc, num_blk, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      ! type(butterfly_UV)::kerls_i,kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist
      character::mode

      n1 = MPI_Wtime()
      nproc = ptree%pgrp(pgno_i)%nproc
      pgno = pgno_i
      do iii=1,2
      do jjj=1,2
      nproc = max(nproc, ptree%pgrp(pgnos_o(iii,jjj))%nproc)
      pgno = min(pgno_i, pgnos_o(iii,jjj))
      enddo
      enddo

      tag = pgno+level_i*10

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_i = block_i%level_butterfly
      else
         level_butterfly_i = -1
         block_i%level_half = -1
      endif
      if (IOwnPgrp(ptree, pgnos_o(1,1))) then
         level_butterfly_o = block_i%level_butterfly
      else
         level_butterfly_o = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_i, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      level_butterfly_c_o = max(level_butterfly_o - 1, 0)
      level_butterfly_o_true = max(level_butterfly_i - 2, 0)
      num_blk = 2**level_butterfly_o

      ! if(level_i<=block_i%level_half)then
      ! mode='R'
      ! level_c_o=0
      ! endif

      ! if(level_i>block_i%level_half)then
      mode = 'C'
      level_c_o = level_butterfly_c_o + 1
      ! endif
      do iii=1,2
      do jjj=1,2
         call assert((ptree%pgrp(pgno_i)%head <= ptree%pgrp(pgnos_o(iii,jjj))%head .and. ptree%pgrp(pgno_i)%tail >= ptree%pgrp(pgnos_o(iii,jjj))%tail) .or. (ptree%pgrp(pgnos_o(iii,jjj))%head <= ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgnos_o(iii,jjj))%tail >= ptree%pgrp(pgno_i)%tail), 'pgno_i or pgno_o should be contained in the other')
      enddo
      enddo

      call GetPgno_Sub(ptree, pgno_i, level_butterfly_i, pgno_sub_i_mine)


      do iii = 1, 2
      do jjj = 1, 2
         pgno_o=pgnos_o(iii,jjj)
         if (IOwnPgrp(ptree, pgno_o)) then
            call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
            call GetLocalBlockRange(ptree, pgno_o, level_butterfly_o_true + 1, level_butterfly_o_true, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

            ! if(mode=='R')then
            ! num_blk=2**level_butterfly_o
            ! idx=idx_c
            ! inc=inc_c
            ! nblk_loc=nc
            ! elseif(mode=='C')then
            if (level_butterfly_c_o == 0) then
               idx = idx_r
               inc = inc_r
               nblk_loc = nr
            else
               idx = idx_r*2 - 1
               inc = inc_r
               nblk_loc = nr*2
            endif
            ! endif
            if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
               if (nblk_loc > 0) then
                  block_o%sons(iii, jjj)%ButterflyU%idx = idx + (iii - 1)*num_blk/2
                  block_o%sons(iii, jjj)%ButterflyU%inc = inc
                  block_o%sons(iii, jjj)%ButterflyU%nblk_loc = nblk_loc
                  block_o%sons(iii, jjj)%ButterflyU%num_blk = num_blk
                  allocate (block_o%sons(iii, jjj)%ButterflyU%blocks(block_o%sons(iii, jjj)%ButterflyU%nblk_loc))
               endif
            endif
         endif
      enddo
      enddo



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
      do iii = 1, 2
      do jjj = 1, 2
         pgno_o=pgnos_o(iii,jjj)
         if (IOwnPgrp(ptree, pgno_o)) then
            call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
            call GetLocalBlockRange(ptree, pgno_o, level_butterfly_o_true + 1, level_butterfly_o_true, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

            ! if(mode=='R')then
            ! num_blk=2**level_butterfly_o
            ! idx=idx_c
            ! inc=inc_c
            ! nblk_loc=nc
            ! elseif(mode=='C')then
            if (level_butterfly_c_o == 0) then
               idx = idx_r
               inc = inc_r
               nblk_loc = nr
            else
               idx = idx_r*2 - 1
               inc = inc_r
               nblk_loc = nr*2
            endif
            if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
               do ii = 1, nblk_loc
                  ! convert indices from output to input
                  ! if(mode=='R')then
                  ! index_i = 1
                  ! index_j = (ii-1)*block_o%sons(iii,jjj)%ButterflyU%inc+block_o%sons(iii,jjj)%ButterflyU%idx
                  ! elseif(mode=='C')then
                  index_i = (ii - 1)*block_o%sons(iii, jjj)%ButterflyU%inc + block_o%sons(iii, jjj)%ButterflyU%idx
                  index_j = 1
                  ! endif
                  call GetBlockPID(ptree, pgno_i, level_i, level_butterfly_i, index_i, index_j, mode, pgno_sub)
                  pid = ptree%pgrp(pgno_sub)%head
                  pp = pid - ptree%pgrp(pgno)%head + 1
                  if (recvquant(pp)%active == 0) then
                     recvquant(pp)%active = 1
                     Nrecvactive = Nrecvactive + 1
                     recvIDactive(Nrecvactive) = pp
                  endif
               enddo
            endif
         endif
      enddo
      enddo


      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then
      do ii = 1, block_i%ButterflyU%nblk_loc
         ! convert indices from input to output
         ! if(mode=='R')then
         ! index_i = 1
         ! index_j = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
         ! index_i0 = 1
         ! index_j0 = index_j
         ! if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
         ! elseif(mode=='C')then
         index_i = (ii - 1)*block_i%ButterflyU%inc + block_i%ButterflyU%idx
         index_j = 1
         index_j0 = 1
         index_i0 = index_i
         si=1
         if (index_i0 > num_blk/2)then
            index_i0 = index_i0 - num_blk/2
            si=2
         endif
         ! endif
         do iii = si, si
            do jjj = 1, 2
               pgno_o=pgnos_o(iii,jjj)
               call GetBlockPID(ptree, pgno_o, level_c_o, level_butterfly_o_true, ceiling_safe(index_i0/2d0), index_j0, mode, pgno_sub)
               pid = ptree%pgrp(pgno_sub)%head

               pp = pid - ptree%pgrp(pgno)%head + 1
               if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
               endif
               sendquant(pp)%size = sendquant(pp)%size + 6 + size(block_i%ButterflyU%blocks(ii)%matrix, 1)*size(block_i%ButterflyU%blocks(ii)%matrix, 2)
            enddo
         enddo
      enddo
      endif
      endif

      ! communicate receive buffer sizes
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if(recvid /= ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if(sendid /= ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then
         ! pack the send buffer in the second pass
         do ii = 1, block_i%ButterflyU%nblk_loc
            ! convert indices from input to output
            ! if(mode=='R')then
            ! index_i = 1
            ! index_j = (ii-1)*block_i%ButterflyU%inc+block_i%ButterflyU%idx
            ! index_i0 = 1
            ! index_j0 = index_j
            ! if(index_j0>num_blk/2)index_j0=index_j0-num_blk/2
            ! elseif(mode=='C')then
            index_i = (ii - 1)*block_i%ButterflyU%inc + block_i%ButterflyU%idx
            index_j = 1
            index_j0 = 1
            index_i0 = index_i
            si=1
            if (index_i0 > num_blk/2)then
               index_i0 = index_i0 - num_blk/2
               si=2
            endif
            ! endif
            do iii = si, si
            do jjj = 1, 2
               pgno_o=pgnos_o(iii,jjj)
               call GetBlockPID(ptree, pgno_o, level_c_o, level_butterfly_o_true, ceiling_safe(index_i0/2d0), index_j0, mode, pgno_sub)
               pid = ptree%pgrp(pgno_sub)%head

               pp = pid - ptree%pgrp(pgno)%head + 1
               Nrow = size(block_i%ButterflyU%blocks(ii)%matrix, 1)
               Ncol = size(block_i%ButterflyU%blocks(ii)%matrix, 2)

               sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
               sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
               sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
               sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
               sendquant(pp)%dat(sendquant(pp)%size + 5, 1) = iii
               sendquant(pp)%dat(sendquant(pp)%size + 6, 1) = jjj
               sendquant(pp)%size = sendquant(pp)%size + 6
               do i = 1, Nrow*Ncol
                  rr = mod(i - 1, Nrow) + 1
                  cc = (i - 1)/Nrow + 1
                  sendquant(pp)%dat(sendquant(pp)%size + i, 1) = block_i%ButterflyU%blocks(ii)%matrix(rr, cc)
               enddo
               sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol
               ! deallocate(block_i%ButterflyU%blocks(ii)%matrix)
            enddo
            enddo
         enddo
         ! if(allocated(block_i%ButterflyU%blocks))deallocate(block_i%ButterflyU%blocks)
      endif
      endif

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            iii = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            jjj = NINT(dble(recvquant(pp)%dat(i, 1)))
            ii = (index_i - block_o%sons(iii, jjj)%ButterflyU%idx)/block_o%sons(iii, jjj)%ButterflyU%inc + 1
            call assert(.not. associated(block_o%sons(iii, jjj)%ButterflyU%blocks(ii)%matrix), 'receiving dat alreay exists locally')
            allocate (block_o%sons(iii, jjj)%ButterflyU%blocks(ii)%matrix(Nrow, Ncol))
            j=0
            do cc = 1, Ncol
               do rr = 1, Nrow
                  j = j + 1
                  block_o%sons(iii, jjj)%ButterflyU%blocks(ii)%matrix(rr, cc) = recvquant(pp)%dat(i+j, 1)
               enddo
            enddo
            i = i + Nrow*Ncol
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_U_split

!>*********** all to all communication of one level of a butterfly to four children butterflies from an old process pgno_i to an new process group pgno_o
!**  it is also assummed row-wise ordering mapped to row-wise ordering, column-wise ordering mapped to column-wise ordering
   subroutine BF_all2all_V_split(block_i, pgno_i, level_i, block_o, pgnos_o, level_o, stats, ptree)


      implicit none
      integer pgno_sub, pgno_i, pgno_o, pgno, level_i, level_o, level_c_o
      integer i, j, level_butterfly_i, level_butterfly_o, level_butterfly_c_o, level_butterfly_o_true, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_loc_k, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer pgnos_o(2,2)
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1, iii, jjj, si, sj, ni, nj, pgno_sub_i_mine, pgno_sub_o_mine
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::block_i, block_o
      type(Hstat)::stats
      type(proctree)::ptree

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Ncol, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, idx, inc_r, inc_c, inc, nr, nc, nblk_loc, num_blk, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      real(kind=8)::n1, n2
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      ! type(butterfly_UV)::kerls_i,kerls_o
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist
      character::mode

      n1 = MPI_Wtime()

      nproc = ptree%pgrp(pgno_i)%nproc
      pgno = pgno_i
      do iii=1,2
      do jjj=1,2
      nproc = max(nproc, ptree%pgrp(pgnos_o(iii,jjj))%nproc)
      pgno = min(pgno_i, pgnos_o(iii,jjj))
      enddo
      enddo
      tag = pgno+level_i*10

      if (IOwnPgrp(ptree, pgno_i)) then
         level_butterfly_i = block_i%level_butterfly
      else
         level_butterfly_i = -1
         block_i%level_half = -1
      endif
      if (IOwnPgrp(ptree, pgnos_o(1,1))) then
         level_butterfly_o = block_i%level_butterfly
      else
         level_butterfly_o = -1
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_i, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, level_butterfly_o, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%level_half, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

      level_butterfly_c_o = max(level_butterfly_o - 1, 0)
      level_butterfly_o_true = max(level_butterfly_i - 2, 0)
      num_blk = 2**level_butterfly_o

      ! if(level_i<=block_i%level_half)then
      mode = 'R'
      level_c_o = 0
      ! endif

      ! if(level_i>block_i%level_half)then
      ! mode='C'
      ! level_c_o=level_butterfly_c_o+1
      ! endif

      do iii=1,2
      do jjj=1,2
          call assert((ptree%pgrp(pgno_i)%head <= ptree%pgrp(pgnos_o(iii,jjj))%head .and. ptree%pgrp(pgno_i)%tail >= ptree%pgrp(pgnos_o(iii,jjj))%tail) .or. (ptree%pgrp(pgnos_o(iii,jjj))%head <= ptree%pgrp(pgno_i)%head .and. ptree%pgrp(pgnos_o(iii,jjj))%tail >= ptree%pgrp(pgno_i)%tail), 'pgno_i or pgno_o should be contained in the other')
      enddo
      enddo
      call GetPgno_Sub(ptree, pgno_i, level_butterfly_i, pgno_sub_i_mine)


  do iii = 1, 2
  do jjj = 1, 2
      pgno_o=pgnos_o(iii,jjj)
      if (IOwnPgrp(ptree, pgno_o)) then
          call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
          call GetLocalBlockRange(ptree, pgno_o, level_c_o, level_butterfly_o_true, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

          if (level_butterfly_c_o == 0) then
              idx = idx_c
              inc = inc_c
              nblk_loc = nc
          else
              idx = idx_c*2 - 1
              inc = inc_c
              nblk_loc = nc*2
          endif

          if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
              if (nblk_loc > 0) then
                  block_o%sons(iii, jjj)%ButterflyV%idx = idx + (jjj - 1)*num_blk/2
                  block_o%sons(iii, jjj)%ButterflyV%inc = inc
                  block_o%sons(iii, jjj)%ButterflyV%nblk_loc = nblk_loc
                  block_o%sons(iii, jjj)%ButterflyV%num_blk = num_blk
                  allocate (block_o%sons(iii, jjj)%ButterflyV%blocks(block_o%sons(iii, jjj)%ButterflyV%nblk_loc))
              endif
          endif
      endif
  enddo
  enddo



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
      do iii = 1, 2
      do jjj = 1, 2
          pgno_o=pgnos_o(iii,jjj)
          if (IOwnPgrp(ptree, pgno_o)) then
              call GetPgno_Sub(ptree, pgno_o, level_butterfly_o_true, pgno_sub_o_mine)
              call GetLocalBlockRange(ptree, pgno_o, level_c_o, level_butterfly_o_true, idx_r, inc_r, nr, idx_c, inc_c, nc, mode)

              if (level_butterfly_c_o == 0) then
                  idx = idx_c
                  inc = inc_c
                  nblk_loc = nc
              else
                  idx = idx_c*2 - 1
                  inc = inc_c
                  nblk_loc = nc*2
              endif
              if (ptree%pgrp(pgno_sub_o_mine)%head == ptree%MyID) then
                  do ii = 1, nblk_loc
                      ! convert indices from output to input
                      ! if(mode=='R')then
                      index_i = 1
                      index_j = (ii - 1)*block_o%sons(iii, jjj)%ButterflyV%inc + block_o%sons(iii, jjj)%ButterflyV%idx
                      ! elseif(mode=='C')then
                      ! index_i = (ii-1)*block_o%sons(iii,jjj)%ButterflyV%inc+block_o%sons(iii,jjj)%ButterflyV%idx
                      ! index_j = 1
                      ! endif
                      call GetBlockPID(ptree, pgno_i, level_i, level_butterfly_i, index_i, index_j, mode, pgno_sub)
                      pid = ptree%pgrp(pgno_sub)%head
                      pp = pid - ptree%pgrp(pgno)%head + 1
                      if (recvquant(pp)%active == 0) then
                          recvquant(pp)%active = 1
                          Nrecvactive = Nrecvactive + 1
                          recvIDactive(Nrecvactive) = pp
                      endif
                  enddo
              endif
          endif
      enddo
      enddo

      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then
      do ii = 1, block_i%ButterflyV%nblk_loc
         ! convert indices from input to output
         ! if(mode=='R')then
         index_i = 1
         index_j = (ii - 1)*block_i%ButterflyV%inc + block_i%ButterflyV%idx
         index_i0 = 1
         index_j0 = index_j
         sj=1
         if (index_j0 > num_blk/2)then
            index_j0 = index_j0 - num_blk/2
            sj=2
         endif
         ! elseif(mode=='C')then
         ! index_i = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
         ! index_j = 1
         ! index_j0 = 1
         ! index_i0 = index_i
         ! if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
         ! endif
         do iii = 1, 2
          do jjj = sj, sj
             pgno_o=pgnos_o(iii,jjj)
              call GetBlockPID(ptree, pgno_o, level_c_o, level_butterfly_o_true, index_i0, ceiling_safe(index_j0/2d0), mode, pgno_sub)
              pid = ptree%pgrp(pgno_sub)%head

              pp = pid - ptree%pgrp(pgno)%head + 1
              if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
              endif
              sendquant(pp)%size = sendquant(pp)%size + 6 + size(block_i%ButterflyV%blocks(ii)%matrix, 1)*size(block_i%ButterflyV%blocks(ii)%matrix, 2)
              enddo
          enddo
      enddo
      endif
      endif

      ! communicate receive buffer sizes
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if(recvid /=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if(sendid /=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      if (IOwnPgrp(ptree, pgno_i)) then
      if (ptree%pgrp(pgno_sub_i_mine)%head == ptree%MyID) then
         ! pack the send buffer in the second pass
         do ii = 1, block_i%ButterflyV%nblk_loc
            ! convert indices from input to output
            ! if(mode=='R')then
            index_i = 1
            index_j = (ii - 1)*block_i%ButterflyV%inc + block_i%ButterflyV%idx
            index_i0 = 1
            index_j0 = index_j
            sj=1
            if (index_j0 > num_blk/2)then
               index_j0 = index_j0 - num_blk/2
               sj=2
            endif
            ! elseif(mode=='C')then
            ! index_i = (ii-1)*block_i%ButterflyV%inc+block_i%ButterflyV%idx
            ! index_j = 1
            ! index_j0 = 1
            ! index_i0 = index_i
            ! if(index_i0>num_blk/2)index_i0=index_i0-num_blk/2
            ! endif
            do iii = 1, 2
              do jjj = sj, sj
                 pgno_o=pgnos_o(iii,jjj)
                  call GetBlockPID(ptree, pgno_o, level_c_o, level_butterfly_o_true, index_i0, ceiling_safe(index_j0/2d0), mode, pgno_sub)
                  pid = ptree%pgrp(pgno_sub)%head

                  pp = pid - ptree%pgrp(pgno)%head + 1
                  Nrow = size(block_i%ButterflyV%blocks(ii)%matrix, 1)
                  Ncol = size(block_i%ButterflyV%blocks(ii)%matrix, 2)

                  sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = index_i
                  sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = index_j
                  sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = Nrow
                  sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = Ncol
                  sendquant(pp)%dat(sendquant(pp)%size + 5, 1) = iii
                  sendquant(pp)%dat(sendquant(pp)%size + 6, 1) = jjj
                  sendquant(pp)%size = sendquant(pp)%size + 6
                  do i = 1, Nrow*Ncol
                      rr = mod(i - 1, Nrow) + 1
                      cc = (i - 1)/Nrow + 1
                      sendquant(pp)%dat(sendquant(pp)%size + i, 1) = block_i%ButterflyV%blocks(ii)%matrix(rr, cc)
                  enddo
                  sendquant(pp)%size = sendquant(pp)%size + Nrow*Ncol

                  ! write(*,*)pid,iii,jjj,'send', index_i, index_j,Nrow,Ncol

                  ! deallocate(block_i%ButterflyV%blocks(ii)%matrix)
              enddo
              enddo
          enddo
      endif
      endif
      ! if(allocated(block_i%ButterflyV%blocks))deallocate(block_i%ButterflyV%blocks)

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            index_i = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            index_j = NINT(dble(recvquant(pp)%dat(i, 1)))

            i = i + 1
            Nrow = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            Ncol = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            iii = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            jjj = NINT(dble(recvquant(pp)%dat(i, 1)))
            ! write(*,*)ptree%MyID,iii,jjj,'recv', index_i, index_j,Nrow,Ncol,block_o%sons(iii, jjj)%ButterflyV%idx, block_o%sons(iii, jjj)%ButterflyV%inc,block_i%level_butterfly,num_blk
          ii = (index_j - block_o%sons(iii, jjj)%ButterflyV%idx)/block_o%sons(iii, jjj)%ButterflyV%inc + 1
              call assert(.not. associated(block_o%sons(iii, jjj)%ButterflyV%blocks(ii)%matrix), 'receiving dat alreay exists locally V')
              allocate (block_o%sons(iii, jjj)%ButterflyV%blocks(ii)%matrix(Nrow, Ncol))
              j=0
              do cc = 1, Ncol
              do rr = 1, Nrow
                  j = j + 1
                  block_o%sons(iii, jjj)%ButterflyV%blocks(ii)%matrix(rr, cc) = recvquant(pp)%dat(i+j, 1)
              enddo
              enddo
            i = i + Nrow*Ncol
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_all2all_V_split







   !>***** redistribute the input and output vectors in Bplus_block_MVP_dat between layout of bplus%LL(1)%matrices_block(1) to the layout of bplus%LL(xx)%matrices_block(yy)
   subroutine Bplus_vec_1Dto1D(BP, rowcol, one2all, level_s, level_e, ld1, dat_1, Nrnd, vecs, ptree,nproc)


      implicit none
      type(blockplus)::BP
      integer::rowcol ! whether the distribution is along row or column dimension
      integer::one2all ! whether the redistribution is from BP to the individual BFs or the other way
      integer:: level_s, level_e
      integer ld1,Nrnd
      DT::dat_1(ld1, *)
      type(vectorsblock_oneL)::vecs(level_s:level_e)
      type(proctree)::ptree
      integer ll,bb,offset_s,offset_r,sizen, head_i, head_o
      type(matrixblock),pointer::blocks,blocks_1
      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      integer::sendIDactive(nproc), recvIDactive(nproc)
      real(kind=8)::n1, n2
      integer Nsendactive, Nrecvactive
      integer rr, cc
      integer i, j, ii, jj, ij, pp, tt
      integer ierr, tag, nproc, nproc_i, nproc_o, Nreqr, Nreqs, recvid, sendid, idxs_i, idxe_i, idxs_o, idxe_o


      n1 = MPI_Wtime()


      blocks_1 => BP%LL(1)%matrices_block(1)
      tag = blocks_1%pgno+1000

      ! allocation of communication quantities
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo

      Nsendactive = 0
      Nrecvactive = 0

      if(one2all==1)then
         nproc_i = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_i==nproc,'nproc should equal nproc_i')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_i = blocks_1%headm
            idxs_i = blocks_1%M_p(ii, 1) + head_i
            idxe_i = blocks_1%M_p(ii, 2) + head_i
         else
            head_i = blocks_1%headn
            idxs_i = blocks_1%N_p(ii, 1) + head_i
            idxe_i = blocks_1%N_p(ii, 2) + head_i
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_o
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1) + head_o
                     idxe_o = blocks%M_p(jj, 2) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1) + head_o
                     idxe_o = blocks%N_p(jj, 2) + head_o
                  endif
                  if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     if (sendquant(pp)%active == 0) then
                        sendquant(pp)%active = 1
                        Nsendactive = Nsendactive + 1
                        sendIDactive(Nsendactive) = pp
                     endif
                     offset_s = max(idxs_i, idxs_o) - idxs_i
                     sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                     sendquant(pp)%size = sendquant(pp)%size + 4 + sizen*Nrnd
                  endif
               enddo
            enddo
         enddo

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1) + head_o
                     idxe_o = blocks%M_p(jj, 2) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1) + head_o
                     idxe_o = blocks%N_p(jj, 2) + head_o
                  endif
                  do ii = 1, nproc_i
                     if(rowcol==1)then
                        head_i = blocks_1%headm
                        idxs_i = blocks_1%M_p(ii, 1) + head_i
                        idxe_i = blocks_1%M_p(ii, 2) + head_i
                     else
                        head_i = blocks_1%headn
                        idxs_i = blocks_1%N_p(ii, 1) + head_i
                        idxe_i = blocks_1%N_p(ii, 2) + head_i
                     endif
                     if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                        pp = ii
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        recvquant(pp)%size = recvquant(pp)%size + 4 + sizen*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      else

         nproc_o = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_o==nproc,'nproc should equal nproc_o')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_o = blocks_1%headm
            idxs_o = blocks_1%M_p(ii, 1) + head_o
            idxe_o = blocks_1%M_p(ii, 2) + head_o
         else
            head_o = blocks_1%headn
            idxs_o = blocks_1%N_p(ii, 1) + head_o
            idxe_o = blocks_1%N_p(ii, 2) + head_o
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_i
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1) + head_i
                     idxe_i = blocks%M_p(jj, 2) + head_i
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1) + head_i
                     idxe_i = blocks%N_p(jj, 2) + head_i
                  endif
                  if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     if (recvquant(pp)%active == 0) then
                        recvquant(pp)%active = 1
                        Nrecvactive = Nrecvactive + 1
                        recvIDactive(Nrecvactive) = pp
                     endif
                     sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                     recvquant(pp)%size = recvquant(pp)%size + 4 + sizen*Nrnd
                  endif
               enddo
            enddo
         enddo

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1) + head_i
                     idxe_i = blocks%M_p(jj, 2) + head_i
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1) + head_i
                     idxe_i = blocks%N_p(jj, 2) + head_i
                  endif
                  do ii = 1, nproc_o
                     if(rowcol==1)then
                        head_o = blocks_1%headm
                        idxs_o = blocks_1%M_p(ii, 1) + head_o
                        idxe_o = blocks_1%M_p(ii, 2) + head_o
                     else
                        head_o = blocks_1%headn
                        idxs_o = blocks_1%N_p(ii, 1) + head_o
                        idxe_o = blocks_1%N_p(ii, 2) + head_o
                     endif
                     if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                        pp = ii
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        offset_s = max(idxs_i, idxs_o) - idxs_i
                        sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                        sendquant(pp)%size = sendquant(pp)%size + 4 + sizen*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass


      if(one2all==1)then
         nproc_i = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_i==nproc,'nproc should equal nproc_i')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_i = blocks_1%headm
            idxs_i = blocks_1%M_p(ii, 1) + head_i
            idxe_i = blocks_1%M_p(ii, 2) + head_i
         else
            head_i = blocks_1%headn
            idxs_i = blocks_1%N_p(ii, 1) + head_i
            idxe_i = blocks_1%N_p(ii, 2) + head_i
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_o
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1) + head_o
                     idxe_o = blocks%M_p(jj, 2) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1) + head_o
                     idxe_o = blocks%N_p(jj, 2) + head_o
                  endif
                  if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     offset_s = max(idxs_i, idxs_o) - idxs_i
                     offset_r = max(idxs_i, idxs_o) - idxs_o
                     sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1

                     sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = ll
                     sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = bb
                     sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = offset_r
                     sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = sizen
                     sendquant(pp)%size = sendquant(pp)%size + 4
                     do i = 1, sizen*Nrnd
                        rr = mod(i - 1, sizen) + 1
                        cc = (i - 1)/sizen + 1
                        sendquant(pp)%dat(sendquant(pp)%size + i, 1) = dat_1(rr+offset_s, cc)
                     enddo
                     sendquant(pp)%size = sendquant(pp)%size + sizen*Nrnd
                  endif
               enddo
            enddo
         enddo
      else
         nproc_o = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_o==nproc,'nproc should equal nproc_o')

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1) + head_i
                     idxe_i = blocks%M_p(jj, 2) + head_i
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1) + head_i
                     idxe_i = blocks%N_p(jj, 2) + head_i
                  endif
                  do ii = 1, nproc_o
                     if(rowcol==1)then
                        head_o = blocks_1%headm
                        idxs_o = blocks_1%M_p(ii, 1) + head_o
                        idxe_o = blocks_1%M_p(ii, 2) + head_o
                     else
                        head_o = blocks_1%headn
                        idxs_o = blocks_1%N_p(ii, 1) + head_o
                        idxe_o = blocks_1%N_p(ii, 2) + head_o
                     endif
                     if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                        pp = ii
                        offset_s = max(idxs_i, idxs_o) - idxs_i
                        offset_r = max(idxs_i, idxs_o) - idxs_o
                        sizen = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1

                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = ll
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = bb
                        sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = offset_r
                        sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = sizen
                        sendquant(pp)%size = sendquant(pp)%size + 4
                        do i = 1, sizen*Nrnd
                           rr = mod(i - 1, sizen) + 1
                           cc = (i - 1)/sizen + 1
                           sendquant(pp)%dat(sendquant(pp)%size + i, 1) = vecs(ll)%vs(bb)%vector(rr+offset_s, cc)
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + sizen*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks_1%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks_1%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks_1%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks_1%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(blocks_1%pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            ll = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            bb = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            offset_r = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            sizen = NINT(dble(recvquant(pp)%dat(i, 1)))
            if(one2all==1)then
               do cc = 1, Nrnd
                  do rr = 1, sizen
                     i = i + 1
                     vecs(ll)%vs(bb)%vector(rr+offset_r, cc)= recvquant(pp)%dat(i, 1)
                  enddo
               enddo
            else
               do cc = 1, Nrnd
                  do rr = 1, sizen
                     i = i + 1
                     dat_1(rr+offset_r, cc)= dat_1(rr+offset_r, cc) + recvquant(pp)%dat(i, 1)
                  enddo
               enddo
            endif
         enddo
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine Bplus_vec_1Dto1D




   !>***** redistribute the input and output vectors in Bplus_MD_block_MVP_dat between layout of bplus%LL(1)%matrices_block(1) to the layout of bplus%LL(xx)%matrices_block(yy)
   subroutine Bplus_MD_vec_1Dto1D(Ndim, BP, rowcol, one2all, level_s, level_e, ld1, dat_1, Nrnd, vecs, ptree,nproc)


      implicit none
      integer Ndim
      type(blockplus_MD)::BP
      integer::rowcol ! whether the distribution is along row or column dimension
      integer::one2all ! whether the redistribution is from BP to the individual BFs or the other way
      integer:: level_s, level_e, dim_i
      integer ld1(Ndim),Nrnd
      DT::dat_1(product(ld1), *)
      type(vectorsblock_oneL)::vecs(level_s:level_e)
      type(proctree)::ptree
      integer ll,bb,offset_s(Ndim),offset_r(Ndim),sizen(Ndim), head_i(Ndim), head_o(Ndim)
      type(matrixblock_MD),pointer::blocks,blocks_1
      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      integer::sendIDactive(nproc), recvIDactive(nproc)
      real(kind=8)::n1, n2
      integer Nsendactive, Nrecvactive
      integer rr, cc
      integer i, j, ii, jj, ij, pp, tt, iii
      integer ierr, tag, nproc, nproc_i, nproc_o, Nreqr, Nreqs, recvid, sendid, idxs_i(Ndim), idxe_i(Ndim), idxs_o(Ndim), idxe_o(Ndim), idx_MD(Ndim), idx_s_MD(Ndim),idx_r_MD(Ndim), dims_md(Ndim),dims_md_r(Ndim),idx_r_scalar,idx_s_scalar


      n1 = MPI_Wtime()


      blocks_1 => BP%LL(1)%matrices_block(1)
      tag = blocks_1%pgno+1000

      ! allocation of communication quantities
      do ii = 1, nproc
         sendquant(ii)%size = 0
         sendquant(ii)%active = 0
      enddo
      do ii = 1, nproc
         recvquant(ii)%size = 0
         recvquant(ii)%active = 0
      enddo

      Nsendactive = 0
      Nrecvactive = 0

      if(one2all==1)then
         nproc_i = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_i==nproc,'nproc should equal nproc_i')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_i = blocks_1%headm
            idxs_i = blocks_1%M_p(ii, 1, :) + head_i
            idxe_i = blocks_1%M_p(ii, 2, :) + head_i
         else
            head_i = blocks_1%headn
            idxs_i = blocks_1%N_p(ii, 1, :) + head_i
            idxe_i = blocks_1%N_p(ii, 2, :) + head_i
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_o
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1, :) + head_o
                     idxe_o = blocks%M_p(jj, 2, :) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1, :) + head_o
                     idxe_o = blocks%N_p(jj, 2, :) + head_o
                  endif
                  if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     if (sendquant(pp)%active == 0) then
                        sendquant(pp)%active = 1
                        Nsendactive = Nsendactive + 1
                        sendIDactive(Nsendactive) = pp
                     endif
                     do dim_i=1,Ndim
                        offset_s(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_i(dim_i)
                        sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                     enddo
                     sendquant(pp)%size = sendquant(pp)%size + 2 + 2*Ndim + product(sizen)*Nrnd
                  endif
               enddo
            enddo
         enddo

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1, :) + head_o
                     idxe_o = blocks%M_p(jj, 2, :) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1, :) + head_o
                     idxe_o = blocks%N_p(jj, 2, :) + head_o
                  endif
                  do ii = 1, nproc_i
                     if(rowcol==1)then
                        head_i = blocks_1%headm
                        idxs_i = blocks_1%M_p(ii, 1, :) + head_i
                        idxe_i = blocks_1%M_p(ii, 2, :) + head_i
                     else
                        head_i = blocks_1%headn
                        idxs_i = blocks_1%N_p(ii, 1, :) + head_i
                        idxe_i = blocks_1%N_p(ii, 2, :) + head_i
                     endif
                     if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                        pp = ii
                        if (recvquant(pp)%active == 0) then
                           recvquant(pp)%active = 1
                           Nrecvactive = Nrecvactive + 1
                           recvIDactive(Nrecvactive) = pp
                        endif
                        do dim_i=1,Ndim
                           sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                        enddo
                        recvquant(pp)%size = recvquant(pp)%size + 2 + 2*Ndim + product(sizen)*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      else

         nproc_o = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_o==nproc,'nproc should equal nproc_o')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_o = blocks_1%headm
            idxs_o = blocks_1%M_p(ii, 1, :) + head_o
            idxe_o = blocks_1%M_p(ii, 2, :) + head_o
         else
            head_o = blocks_1%headn
            idxs_o = blocks_1%N_p(ii, 1, :) + head_o
            idxe_o = blocks_1%N_p(ii, 2, :) + head_o
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_i
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1, :) + head_i
                     idxe_i = blocks%M_p(jj, 2, :) + head_i
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1, :) + head_i
                     idxe_i = blocks%N_p(jj, 2, :) + head_i
                  endif
                  if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     if (recvquant(pp)%active == 0) then
                        recvquant(pp)%active = 1
                        Nrecvactive = Nrecvactive + 1
                        recvIDactive(Nrecvactive) = pp
                     endif
                     do dim_i=1,Ndim
                        sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                     enddo
                     recvquant(pp)%size = recvquant(pp)%size + 2 + 2*Ndim + product(sizen)*Nrnd
                  endif
               enddo
            enddo
         enddo

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1, :) + head_i
                     idxe_i = blocks%M_p(jj, 2, :) + head_i
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1, :) + head_i
                     idxe_i = blocks%N_p(jj, 2, :) + head_i
                  endif
                  do ii = 1, nproc_o
                     if(rowcol==1)then
                        head_o = blocks_1%headm
                        idxs_o = blocks_1%M_p(ii, 1, :) + head_o
                        idxe_o = blocks_1%M_p(ii, 2, :) + head_o
                     else
                        head_o = blocks_1%headn
                        idxs_o = blocks_1%N_p(ii, 1, :) + head_o
                        idxe_o = blocks_1%N_p(ii, 2, :) + head_o
                     endif
                     if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                        pp = ii
                        if (sendquant(pp)%active == 0) then
                           sendquant(pp)%active = 1
                           Nsendactive = Nsendactive + 1
                           sendIDactive(Nsendactive) = pp
                        endif
                        do dim_i=1,Ndim
                           sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + 2 + 2*Ndim + product(sizen)*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass


      if(one2all==1)then
         nproc_i = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_i==nproc,'nproc should equal nproc_i')

         ii = ptree%myid - ptree%pgrp(blocks_1%pgno)%head + 1
         if(rowcol==1)then
            head_i = blocks_1%headm
            idxs_i = blocks_1%M_p(ii, 1, :) + head_i
            idxe_i = blocks_1%M_p(ii, 2, :) + head_i
         else
            head_i = blocks_1%headn
            idxs_i = blocks_1%N_p(ii, 1, :) + head_i
            idxe_i = blocks_1%N_p(ii, 2, :) + head_i
         endif
         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_o = ptree%pgrp(blocks%pgno)%nproc
               do jj=1,nproc_o
                  if(rowcol==1)then
                     head_o = blocks%headm
                     idxs_o = blocks%M_p(jj, 1, :) + head_o
                     idxe_o = blocks%M_p(jj, 2, :) + head_o
                  else
                     head_o = blocks%headn
                     idxs_o = blocks%N_p(jj, 1, :) + head_o
                     idxe_o = blocks%N_p(jj, 2, :) + head_o
                  endif
                  if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                     pp = ptree%pgrp(blocks%pgno)%head + jj - ptree%pgrp(blocks_1%pgno)%head
                     do dim_i=1,Ndim
                        offset_s(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_i(dim_i)
                        offset_r(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_o(dim_i)
                        sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                     enddo

                     sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = ll
                     sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = bb
                     sendquant(pp)%dat(sendquant(pp)%size + 3:sendquant(pp)%size + 2 + Ndim, 1) = offset_r
                     sendquant(pp)%dat(sendquant(pp)%size + 3 + Ndim: sendquant(pp)%size + 2 + 2*Ndim, 1) = sizen
                     sendquant(pp)%size = sendquant(pp)%size + 2 + 2*Ndim

                     do iii=1,product(sizen)
                        call SingleIndexToMultiIndex(Ndim, sizen, iii, idx_MD)
                        idx_s_MD = idx_MD + offset_s
                        call MultiIndexToSingleIndex(Ndim, ld1, idx_s_scalar, idx_s_MD)
                        do cc=1,Nrnd
                           sendquant(pp)%dat(sendquant(pp)%size+(iii-1)*Nrnd+cc,1) = dat_1(idx_s_scalar, cc)
                        enddo
                     enddo
                     sendquant(pp)%size = sendquant(pp)%size + product(sizen)*Nrnd
                  endif
               enddo
            enddo
         enddo
      else
         nproc_o = ptree%pgrp(blocks_1%pgno)%nproc
         call assert(nproc_o==nproc,'nproc should equal nproc_o')

         do ll=level_s,level_e
            do bb = 1, BP%LL(ll)%Nbound
               blocks => BP%LL(ll)%matrices_block(bb)
               nproc_i = ptree%pgrp(blocks%pgno)%nproc
               if (IOwnPgrp(ptree, blocks%pgno)) then
                  jj = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
                  if(rowcol==1)then
                     head_i = blocks%headm
                     idxs_i = blocks%M_p(jj, 1, :) + head_i
                     idxe_i = blocks%M_p(jj, 2, :) + head_i
                     dims_md = blocks%M_loc
                  else
                     head_i = blocks%headn
                     idxs_i = blocks%N_p(jj, 1, :) + head_i
                     idxe_i = blocks%N_p(jj, 2, :) + head_i
                     dims_md = blocks%N_loc
                  endif
                  do ii = 1, nproc_o
                     if(rowcol==1)then
                        head_o = blocks_1%headm
                        idxs_o = blocks_1%M_p(ii, 1, :) + head_o
                        idxe_o = blocks_1%M_p(ii, 2, :) + head_o
                     else
                        head_o = blocks_1%headn
                        idxs_o = blocks_1%N_p(ii, 1, :) + head_o
                        idxe_o = blocks_1%N_p(ii, 2, :) + head_o
                     endif
                     if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                        pp = ii
                        do dim_i=1,Ndim
                           offset_s(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_i(dim_i)
                           offset_r(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_o(dim_i)
                           sizen(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                        enddo

                        sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = ll
                        sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = bb
                        sendquant(pp)%dat(sendquant(pp)%size + 3:sendquant(pp)%size + 2 + Ndim, 1) = offset_r
                        sendquant(pp)%dat(sendquant(pp)%size + 3 + Ndim: sendquant(pp)%size + 2 + 2*Ndim, 1) = sizen
                        sendquant(pp)%size = sendquant(pp)%size + 2 + 2*Ndim

                        do iii=1,product(sizen)
                           call SingleIndexToMultiIndex(Ndim, sizen, iii, idx_MD)
                           idx_s_MD = idx_MD + offset_s
                           call MultiIndexToSingleIndex(Ndim, dims_md, idx_s_scalar, idx_s_MD)
                           do cc=1,Nrnd
                              sendquant(pp)%dat(sendquant(pp)%size+(iii-1)*Nrnd+cc,1) = vecs(ll)%vs(bb)%vector(idx_s_scalar, cc)
                           enddo
                        enddo
                        sendquant(pp)%size = sendquant(pp)%size + product(sizen)*Nrnd
                     endif
                  enddo
               endif
            enddo
         enddo
      endif

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks_1%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks_1%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks_1%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks_1%pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      ! copy data from buffer to target
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(blocks_1%pgno)%head
         else
            call MPI_waitany(Nreqr, R_req, sendid, statusr(:,1), ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
         endif
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            ll = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            bb = NINT(dble(recvquant(pp)%dat(i, 1)))
            do dim_i=1,Ndim
               i = i + 1
               offset_r(dim_i) = NINT(dble(recvquant(pp)%dat(i, 1)))
            enddo
            do dim_i=1,Ndim
               i = i + 1
               sizen(dim_i) = NINT(dble(recvquant(pp)%dat(i, 1)))
            enddo

            if(one2all==1)then
               blocks => BP%LL(ll)%matrices_block(bb)
               if(rowcol==1)then
                  dims_md_r = blocks%M_loc
               else
                  dims_md_r = blocks%N_loc
               endif
               do iii=1,product(sizen)
                  call SingleIndexToMultiIndex(Ndim, sizen, iii, idx_MD)
                  idx_r_MD = idx_MD + offset_r
                  call MultiIndexToSingleIndex(Ndim, dims_md_r, idx_r_scalar, idx_r_MD)
                  do cc=1,Nrnd
                     i = i + 1
                     vecs(ll)%vs(bb)%vector(idx_r_scalar, cc) = recvquant(pp)%dat(i,1)
                  enddo
               enddo
            else
               do iii=1,product(sizen)
                  call SingleIndexToMultiIndex(Ndim, sizen, iii, idx_MD)
                  idx_r_MD = idx_MD + offset_r
                  call MultiIndexToSingleIndex(Ndim, ld1, idx_r_scalar, idx_r_MD)
                  do cc=1,Nrnd
                     i = i + 1
                     dat_1(idx_r_scalar, cc) = dat_1(idx_r_scalar, cc) + recvquant(pp)%dat(i,1)
                  enddo
               enddo
            endif
         enddo
      enddo

      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif

      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine Bplus_MD_vec_1Dto1D





   subroutine BF_Mult(chara, xin, xout, Ninloc, Noutloc, Ncol, blocks, option, stats, ptree)
      implicit none
      real(kind=8) t1, t2
      character chara
      integer Ninloc, Noutloc, Ncol
      DT::xin(Ninloc, Ncol), xout(Noutloc, Ncol)

      type(Hstat)::stats
      type(Hoption)::option
      type(matrixblock)::blocks

      type(proctree)::ptree

      t1 = MPI_Wtime()

      stats%Flop_Tmp = 0
      stats%Flop_C_Mult = 0
      stats%Time_C_Mult = 0
      xout = 0
      if (chara == 'N') then
         call BF_block_MVP_dat(blocks, chara, Noutloc, Ninloc, Ncol, xin, Ninloc, xout,Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      else
         call BF_block_MVP_dat(blocks, chara, Ninloc, Noutloc, Ncol, xin, Ninloc, xout, Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      endif

      t2 = MPI_Wtime()

      xout = xout/option%scale_factor

      stats%Time_C_Mult = stats%Time_C_Mult + t2 - t1
      stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

   end subroutine BF_Mult




   subroutine BF_block_MVP_dat(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
      implicit none
      type(matrixblock)::blocks
      character chara
      integer M, N, Nrnd
      integer ldi,ldo
      DT :: random1(ldi, *), random2(ldo, *)
      DT :: a, b
      type(proctree)::ptree
      type(Hstat)::stats

#ifdef HAVE_MAGMA
   call BF_block_MVP_dat_batch_magma(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
#else

! #ifdef HAVE_MKL
#if 0
      call BF_block_MVP_dat_batch_mkl(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
#else
      call BF_block_MVP_dat_nonbatch(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)
#endif
#endif
   end subroutine BF_block_MVP_dat



   subroutine BF_block_MVP_dat_nonbatch(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)

      implicit none

      integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, level_half, level_final, pgno_sub
      integer idx_r, inc_r, nr, idx_c, inc_c, nc
      integer idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno, comm, ierr
      type(Hstat)::stats
      real(kind=8)::flop, flops,n1,n2,n3,n4,n5,n6
      integer index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc, index_i_loc_s, index_i_loc_k, index_j_loc, index_j_loc0, index_i_loc0, index_j_loc_s, index_j_loc_k

      type(butterfly_vec) :: BFvec
      integer ldi,ldo
      DT :: random1(ldi, *), random2(ldo, *)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), Vout_tmp(:, :)

      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)

      n1 = MPI_Wtime()

      level_butterfly = blocks%level_butterfly
      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm
      if (comm == MPI_COMM_NULL) then
         write (*, *) 'ninin', pgno, comm == MPI_COMM_NULL, ptree%MyID
      endif

      call assert(IOwnPgrp(ptree, pgno), 'I do not share this block!')

      if (blocks%style == 1) then
         call Full_block_MVP_dat(blocks, chara, M, Nrnd, random1, ldi, random2, ldo, a, b)
         return
      endif

      if (level_butterfly == 0) then
         rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
         call assert(rank > 0, 'rank incorrect in blocks%ButterflyU')
         allocate (matrixtemp(rank, Nrnd))
         matrixtemp = 0
         allocate (matrixtemp1(rank, Nrnd))
         matrixtemp1 = 0
         ! for implementation simplicity, MPI_ALLREDUCE is used even when nproc==1
         if (chara == 'N') then !Vout=U*V^T*Vin
            allocate (Vout_tmp(M,Nrnd))
            Vout_tmp = 0
            call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyV%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 2')
            call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
            call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, M, 'N', 'N', size(blocks%ButterflyU%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            random2(1:M,1:Nrnd) = b*random2(1:M,1:Nrnd) + a*Vout_tmp
         else if (chara == 'T') then !Vout=V*U^T*Vin
            allocate (Vout_tmp(N,Nrnd))
            Vout_tmp = 0
            call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyU%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 3')
            call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
            call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, N, 'N', 'N', size(blocks%ButterflyV%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            random2(1:N,1:Nrnd) = b*random2(1:N,1:Nrnd) + a*Vout_tmp
         endif

         deallocate (matrixtemp)
         deallocate (matrixtemp1)
         deallocate (Vout_tmp)

      else
#ifdef HAVE_TASKLOOP
         !$omp parallel
         !$omp single
#endif
         allocate (arr_acc_n(blocks%ButterflyV%nblk_loc))
         allocate (arr_acc_m(blocks%ButterflyU%nblk_loc))
         k1 = 0
         do i = 1, blocks%ButterflyV%nblk_loc
            arr_acc_n(i) = k1
            nn = size(blocks%ButterflyV%blocks(i)%matrix, 1)
            k1 = k1 + nn
         enddo

         k2 = 0
         do i = 1, blocks%ButterflyU%nblk_loc
            arr_acc_m(i) = k2
            mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
            k2 = k2 + mm
         enddo

         num_vectors = Nrnd

         if (BF_checkNAN(blocks)) then
            write (*, *) 'NAN in 0 BF_block_MVP_dat'
            stop
         end if

         if (chara == 'N') then
            n5 = MPI_Wtime()
            if (isnanMat(random1(1:N,1:1),N,1)) then
               write (*, *) 'NAN in 1 BF_block_MVP_dat'
               stop
            end if

            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(0:level_butterfly + 2))

            allocate (BFvec%vec(0)%blocks(1, blocks%ButterflyV%nblk_loc))
            BFvec%vec(0)%num_row = 1
            BFvec%vec(0)%num_col = num_blocks
            BFvec%vec(0)%idx_r = 1
            BFvec%vec(0)%inc_r = 1
            BFvec%vec(0)%nr = 1
            BFvec%vec(0)%idx_c = blocks%ButterflyV%idx
            BFvec%vec(0)%inc_c = blocks%ButterflyV%inc
            BFvec%vec(0)%nc = blocks%ButterflyV%nblk_loc

            do level = 0, level_half
               ! n1 = MPI_Wtime()
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

               BFvec%vec(level + 1)%idx_r = idx_r
               BFvec%vec(level + 1)%inc_r = inc_r
               BFvec%vec(level + 1)%nr = nr
               BFvec%vec(level + 1)%idx_c = idx_c
               BFvec%vec(level + 1)%inc_c = inc_c
               BFvec%vec(level + 1)%nc = nc

               if (nr > 0 .and. nc > 0) then
                  if (level /= level_butterfly + 1) then
                     BFvec%vec(level + 1)%num_row = 2**level
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
                  else
                     BFvec%vec(level + 1)%num_row = 2**level_butterfly
                     BFvec%vec(level + 1)%num_col = 1
                  endif
                  if (level_half /= level) then ! the last level doesn't require doubling block columns
                  if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                     BFvec%vec(level + 1)%nc = 2
                     BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

                  if (level == 0) then
                     flops = 0
                     n3 = MPI_Wtime()
#ifdef HAVE_TASKLOOP
                     !$omp taskloop default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s)
#endif
                     do j = 1, blocks%ButterflyV%nblk_loc
                        index_j = (j - 1)*inc_c + idx_c
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                        nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                        allocate (BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(rank, num_vectors))
                        BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix = 0

                        call gemmf90(blocks%ButterflyV%blocks(j)%matrix, nn, random1(1 + arr_acc_n(j), 1), ldi, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank, 'T', 'N', rank, num_vectors, nn, BPACK_cone, BPACK_czero, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif
                     enddo
#ifdef HAVE_TASKLOOP
                     !$omp end taskloop
#endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                        index_j = idx_c
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                     endif

                  elseif (level == level_butterfly + 1) then
                     write(*,*)'should not arrive here'

                  else
                     n3 = MPI_Wtime()
                     flops = 0
#ifdef HAVE_TASKLOOP
                     !$omp taskloop default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                        index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        nn2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                        allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                        BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0

                        call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn1, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn1, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif
                        call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix, nn2, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn2, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif
                     enddo
#ifdef HAVE_TASKLOOP
                     !$omp end taskloop
#endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif
               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo

               ! n2 = MPI_Wtime()
               ! time_tmp = time_tmp + n2-n1

               if (level_half /= level) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'B')
               endif
            enddo
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

            if (level_half + 1 /= 0) then
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'R', 'C', 0)
            else
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'R', 'C', 0)
            endif

            n5 = MPI_Wtime()
            do level = level_half + 1, level_butterfly + 1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C')


               ! convert the local column-wise kernel block ranges to local row-wise output vector ranges
               if (level /= 0 .and. level /= level_butterfly + 1) then
                  idx_r = idx_r0*2 - 1
                  nr = nr0*2
                  inc_r = inc_r0
                  idx_c = ceiling_safe(idx_c0/2d0)
                  if (inc_c0 > 1) then
                     nc = nc0
                  else
                     nc = ceiling_safe(nc0/2d0)
                  endif
                  inc_c = ceiling_safe(inc_c0/2d0)
               else
                  idx_r = idx_r0
                  nr = nr0
                  inc_r = inc_r0
                  idx_c = idx_c0
                  nc = nc0
                  inc_c = inc_c0
               endif

               BFvec%vec(level + 1)%idx_r = idx_r
               BFvec%vec(level + 1)%inc_r = inc_r
               BFvec%vec(level + 1)%nr = nr
               BFvec%vec(level + 1)%idx_c = idx_c
               BFvec%vec(level + 1)%inc_c = inc_c
               BFvec%vec(level + 1)%nc = nc

               if (nr > 0 .and. nc > 0) then

                  if (level /= level_butterfly + 1) then
                     BFvec%vec(level + 1)%num_row = 2**level
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
                  else
                     BFvec%vec(level + 1)%num_row = 2**level_butterfly
                     BFvec%vec(level + 1)%num_col = 1
                  endif

                  allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

                  if (level == 0) then
                     write (*, *) 'should not arrive here'
                  elseif (level == level_butterfly + 1) then
                     flops = 0
                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        mm = size(blocks%ButterflyU%blocks(1)%matrix, 1)
                        rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                        allocate (matrixtemp(rank, num_vectors))
                        matrixtemp = 0
                        if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                           index_i = idx_r0
                           index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1
                           matrixtemp = BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix
                        endif
                        call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                        allocate (BFvec%vec(level + 1)%blocks(1, 1)%matrix(mm, num_vectors))
                        BFvec%vec(level + 1)%blocks(1, 1)%matrix = 0
                        call gemmf90(blocks%ButterflyU%blocks(1)%matrix, mm, matrixtemp, rank, random2(1 + arr_acc_m(1), 1), ldo, 'N', 'N', mm, num_vectors, rank, a, b, flop=flop)
                        flops = flops + flop
                        deallocate (matrixtemp)
                     else
                     n3 = MPI_Wtime()
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(i,index_i,index_i_loc_s,rank,mm,flop)
#endif
                        do i = 1, nr0
                           index_i = (i - 1)*inc_r0 + idx_r0
                           index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1

                           rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                           mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                           allocate (BFvec%vec(level + 1)%blocks(i, 1)%matrix(mm, num_vectors))
                           BFvec%vec(level + 1)%blocks(i, 1)%matrix = 0

                           call gemmf90(blocks%ButterflyU%blocks(i)%matrix, mm, BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix, rank, random2(1 + arr_acc_m(i), 1), ldo, 'N', 'N', mm, num_vectors, rank, a, b, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif
                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif
                        n4 = MPI_Wtime()
                        ! time_tmp = time_tmp + n4-n3

                     endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                  else

                     n3 = MPI_Wtime()
                     flops = 0

                     if (nc0 > 1 .and. inc_c0 == 1) then  ! this special treatment makes sure two threads do not write to the same address simultaneously
                        ! n1 = MPI_Wtime()
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,index_j_loc0,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                        do index_ij = 1, nr0*nc0/2
                           index_j_loc0 = (index_ij - 1)/nr0 + 1
                           do jj = 1, 2
                              index_j_loc = 2*(index_j_loc0 - 1) + jj       !index_i_loc is local index of column-wise ordering at current level
                              index_i_loc = mod(index_ij - 1, nr0) + 1
                              index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                              index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                              index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                              index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                              index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                              index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                              index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                              index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                              index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                              if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                                 BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                              endif
                              call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                              !$omp atomic
#endif
                              flops = flops + flop
#ifdef HAVE_TASKLOOP
                              !$omp end atomic
#endif

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 2)
                              ! !$omp critical
                              if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix(mm, num_vectors))
                                 BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix = 0
                              endif
                              ! !$omp end critical
                              call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                              !$omp atomic
#endif
                              flops = flops + flop
#ifdef HAVE_TASKLOOP
                              !$omp end atomic
#endif
                           enddo
                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif

   ! n2 = MPI_Wtime()
   ! time_tmp = time_tmp + n2-n1

                     else
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                        do index_ij = 1, nr0*nc0
                           index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                           index_i_loc = mod(index_ij - 1, nr0) + 1
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                           index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 2)
                           ! !$omp critical
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix = 0
                           endif
                           ! !$omp end critical
                           call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, BFvec%vec(level + 1)%blocks(index_i_loc_s + 1, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif

                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif

                     endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif

               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level /= level_butterfly + 1) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'R')
               endif
            enddo

            if (isnanMat(random2(1:M,1:1),M,1)) then
               write (*, *) ptree%MyID, 'NAN in 2 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if
            !deallocate (BFvec%vec)
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

         elseif (chara == 'T') then
            if (isnanMat(random1(1:M,1:1),M,1)) then
               write (*, *) 'NAN in 3 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if

            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(0:level_butterfly + 2))
            allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc, 1))
            BFvec%vec(0)%num_row = num_blocks
            BFvec%vec(0)%num_col = 1
            BFvec%vec(0)%idx_r = blocks%ButterflyU%idx
            BFvec%vec(0)%inc_r = blocks%ButterflyU%inc
            BFvec%vec(0)%nr = blocks%ButterflyU%nblk_loc
            BFvec%vec(0)%idx_c = 1
            BFvec%vec(0)%inc_c = 1
            BFvec%vec(0)%nc = 1

            n5 = MPI_Wtime()
            do level = level_butterfly + 1, level_half + 1, -1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

               BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
               BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
               BFvec%vec(level_butterfly - level + 2)%nr = nr
               BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
               BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
               BFvec%vec(level_butterfly - level + 2)%nc = nc

               if (nr > 0 .and. nc > 0) then

                  if (level /= 0) then
                     BFvec%vec(level_butterfly - level + 2)%num_row = 2**(level - 1)
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**(level_butterfly - level + 1)
                  else
                     BFvec%vec(level_butterfly - level + 2)%num_row = 1
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**level_butterfly
                  endif
                  if (level_half + 1 /= level) then ! the last level doesn't require doubling block rows
                  if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                     BFvec%vec(level_butterfly - level + 2)%nr = 2
                     BFvec%vec(level_butterfly - level + 2)%idx_r = BFvec%vec(level_butterfly - level + 2)%idx_r - 1 + mod(BFvec%vec(level_butterfly - level + 2)%idx_r, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

                  if (level == level_butterfly + 1) then
                     flops = 0
                     n3 = MPI_Wtime()
#ifdef HAVE_TASKLOOP
                     !$omp taskloop default(shared) private(i,rank,mm,flop,index_i,index_i_loc_s)
#endif
                     do i = 1, blocks%ButterflyU%nblk_loc
                        index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                        index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                        rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                        mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                        allocate (BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(rank, num_vectors))
                        BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix = 0

                        call gemmf90(blocks%ButterflyU%blocks(i)%matrix, mm, random1(1 + arr_acc_m(i), 1), ldi, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank, 'T', 'N', rank, num_vectors, mm, BPACK_cone, BPACK_czero, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif

                     enddo
#ifdef HAVE_TASKLOOP
                     !$omp end taskloop
#endif
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                     stats%Flop_Tmp = stats%Flop_Tmp + flops

                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                        index_i = blocks%ButterflyU%idx
                        index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1
                        call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                     endif

                  elseif (level == 0) then
                     write(*,*)'should not arrive here'
                  else
                     n3=MPI_Wtime()
                     flops = 0
#ifdef HAVE_TASKLOOP
                     !$omp taskloop default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                        index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        mm2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                        BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0

                        call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm1, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm1, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm1, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif
                        call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, mm2, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc + 1, index_jj_loc)%matrix, mm2, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm2, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                        !$omp atomic
#endif
                        flops = flops + flop
#ifdef HAVE_TASKLOOP
                        !$omp end atomic
#endif

                     enddo
#ifdef HAVE_TASKLOOP
                     !$omp end taskloop
#endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif
               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level_half + 1 /= level) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'B')
               endif
            enddo
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

            if (level_half /= level_butterfly + 1) then
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'C', 'R', 0)
            else
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'C', 'R', 0)
            endif

            n5 = MPI_Wtime()
            do level = level_half, 0, -1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'R')

               ! convert the local row-wise kernel block ranges to local column-wise output vector ranges
               if (level /= 0 .and. level /= level_butterfly + 1) then
                  idx_r = ceiling_safe(idx_r0/2d0)
                  if (inc_r0 > 1) then
                     nr = nr0
                  else
                     nr = ceiling_safe(nr0/2d0)
                  endif
                  inc_r = ceiling_safe(inc_r0/2d0)
                  idx_c = idx_c0*2 - 1
                  nc = nc0*2
                  inc_c = inc_c0
               else
                  idx_r = idx_r0
                  nr = nr0
                  inc_r = inc_r0
                  idx_c = idx_c0
                  nc = nc0
                  inc_c = inc_c0
               endif

               BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
               BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
               BFvec%vec(level_butterfly - level + 2)%nr = nr
               BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
               BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
               BFvec%vec(level_butterfly - level + 2)%nc = nc

               if (nr > 0 .and. nc > 0) then

                  if (level /= 0) then
                     BFvec%vec(level + 1)%num_row = 2**(level - 1)
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level + 1)
                  else
                     BFvec%vec(level + 1)%num_row = 1
                     BFvec%vec(level + 1)%num_col = 2**level_butterfly
                  endif

                  allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

                  if (level == level_butterfly + 1) then
                     write(*,*)'should not arrive here'


                  elseif (level == 0) then
                     flops = 0

                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                     n3=MPI_Wtime()
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        nn = size(blocks%ButterflyV%blocks(1)%matrix, 1)
                        rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                        allocate (matrixtemp(rank, num_vectors))
                        matrixtemp = 0
                        if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                           index_j = blocks%ButterflyV%idx
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1

                           matrixtemp = BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix
                        endif
                        call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                        allocate (BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix(nn, num_vectors))
                        BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix = 0
                        call gemmf90(blocks%ButterflyV%blocks(1)%matrix, nn, matrixtemp, rank, random2(1 + arr_acc_n(1), 1), ldo, 'N', 'N', nn, num_vectors, rank, a, b, flop=flop)
                        flops = flops + flop
                        deallocate (matrixtemp)
                     else
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s)
#endif
                        do j = 1, blocks%ButterflyV%nblk_loc
                           index_j = (j - 1)*blocks%ButterflyV%inc + blocks%ButterflyV%idx
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1

                           nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                           rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                           allocate (BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix(nn, num_vectors))
                           BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix = 0
                           call gemmf90(blocks%ButterflyV%blocks(j)%matrix, nn, BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix, rank, random2(1 + arr_acc_n(j), 1), ldo, 'N', 'N', nn, num_vectors, rank, a, b, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif
                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif
                     endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  else

                     flops = 0
                     n3 = MPI_Wtime()
                     if (nr0 > 1 .and. inc_r0 == 1) then ! this special treatment makes sure two threads do not write to the same address simultaneously
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,index_i_loc0,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                        do index_ij = 1, nr0*nc0/2
                           index_i_loc0 = (index_ij - 1)/nc0 + 1
                           do ii = 1, 2
                              index_i_loc = (index_i_loc0 - 1)*2 + ii
                              index_j_loc = mod(index_ij - 1, nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
                              index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                              index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                              index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                              index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                              index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                              index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                              index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                              index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                              index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                              if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                                 BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                              endif
                              ! write(*,*)index_ii_loc,index_jj_loc,shape(BFvec%vec(level_butterfly-level+1)%blocks),index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level_butterfly-level+2)%blocks),'lv:',level,shape(blocks%ButterflyKerl(level)%blocks)
                              call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                              !$omp atomic
#endif
                              flops = flops + flop
#ifdef HAVE_TASKLOOP
                              !$omp end atomic
#endif

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                              ! !$omp critical
                              if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                                 allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                                 BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                              endif
                              ! !$omp end critical
                              call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, mm, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix, nn, 'T', 'N', nn, num_vectors, mm, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                              !$omp atomic
#endif
                              flops = flops + flop
#ifdef HAVE_TASKLOOP
                              !$omp end atomic
#endif
                           enddo
                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif
                     else
#ifdef HAVE_TASKLOOP
                        !$omp taskloop default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop)
#endif
                        do index_ij = 1, nr0*nc0
                           index_j_loc = (index_ij - 1)/nr0 + 1
                           index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                           index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           ! write(*,*)index_ii_loc,index_jj_loc,shape(BFvec%vec(level_butterfly-level+1)%blocks),index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level_butterfly-level+2)%blocks),'lv:',level,shape(blocks%ButterflyKerl(level)%blocks)
                           call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                           endif
                           call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, mm, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix, nn, 'T', 'N', nn, num_vectors, mm, BPACK_cone, BPACK_cone, flop=flop)
#ifdef HAVE_TASKLOOP
                           !$omp atomic
#endif
                           flops = flops + flop
#ifdef HAVE_TASKLOOP
                           !$omp end atomic
#endif
                        enddo
#ifdef HAVE_TASKLOOP
                        !$omp end taskloop
#endif

                     endif

                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif

               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level /= 0) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'R')
               endif
            enddo

            if (isnanMat(random2(1:N,1:1),N,1)) then
               write (*, *) 'NAN in 4 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5
         endif

         do level = 0, level_butterfly + 2
            do j = 1, BFvec%vec(level)%nc
               do i = 1, BFvec%vec(level)%nr
                  if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
               enddo
            enddo

            if (allocated(BFvec%vec(level)%blocks))deallocate (BFvec%vec(level)%blocks)
         enddo

         deallocate (BFvec%vec)
         deallocate (arr_acc_m, arr_acc_n)
#ifdef HAVE_TASKLOOP
         !$omp end single
         !$omp end parallel
#endif

      endif

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2-n1
      time_tmp2 = time_tmp2 + n2-n1
      return

   end subroutine BF_block_MVP_dat_nonbatch



#ifdef HAVE_MKL

   subroutine BF_block_MVP_dat_batch_mkl(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)



      implicit none

      integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, level_half, level_final, pgno_sub
      integer idx_r, inc_r, nr, idx_c, inc_c, nc
      integer idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno, comm, ierr
      type(Hstat)::stats
      real(kind=8)::flop, flops,n1,n2,n3,n4,n5,n6
      integer index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc, index_i_loc_s, index_i_loc_k, index_j_loc, index_j_loc0, index_i_loc0, index_j_loc_s, index_j_loc_k
      integer*8:: cnta,cntb,cntc
      character*1,allocatable::transa_array(:),transb_array(:),transc_array(:)
      DT,allocatable::alpha_array(:),beta_array(:)
      integer,allocatable::group_size(:),m_array(:),n_array(:),k_array(:),lda_array(:),ldb_array(:),ldc_array(:)
      integer(c_intptr_t),allocatable::a_array(:),b_array(:),c_array(:)
      integer group_count,cnt

      type(butterfly_vec) :: BFvec
      integer ldi,ldo
      DT,target :: random1(ldi, *), random2(ldo, *)
      DT, pointer :: random1_p(:, :),random2_p(:, :)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), Vout_tmp(:, :)

      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)

      n1 = MPI_Wtime()

      level_butterfly = blocks%level_butterfly
      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm
      if (comm == MPI_COMM_NULL) then
         write (*, *) 'ninin', pgno, comm == MPI_COMM_NULL, ptree%MyID
      endif

      call assert(IOwnPgrp(ptree, pgno), 'I do not share this block!')

      if (blocks%style == 1) then
         call Full_block_MVP_dat(blocks, chara, M, Nrnd, random1, ldi, random2, ldo, a, b)
         return
      endif

      if (level_butterfly == 0) then
         rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
         call assert(rank > 0, 'rank incorrect in blocks%ButterflyU')
         allocate (matrixtemp(rank, Nrnd))
         matrixtemp = 0
         allocate (matrixtemp1(rank, Nrnd))
         matrixtemp1 = 0
         ! for implementation simplicity, MPI_ALLREDUCE is used even when nproc==1
         if (chara == 'N') then !Vout=U*V^T*Vin
            allocate (Vout_tmp(M,Nrnd))
            Vout_tmp = 0
            call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyV%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 2')
            call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
            call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, M, 'N', 'N', size(blocks%ButterflyU%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            random2(1:M,1:Nrnd) = b*random2(1:M,1:Nrnd) + a*Vout_tmp
         else if (chara == 'T') then !Vout=V*U^T*Vin
            allocate (Vout_tmp(N,Nrnd))
            Vout_tmp = 0
            call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyU%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 3')
            call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
            call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, N, 'N', 'N', size(blocks%ButterflyV%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
            stats%Flop_Tmp = stats%Flop_Tmp + flop
            random2(1:N,1:Nrnd) = b*random2(1:N,1:Nrnd) + a*Vout_tmp
         endif

         deallocate (matrixtemp)
         deallocate (matrixtemp1)
         deallocate (Vout_tmp)

      else



         allocate (arr_acc_n(blocks%ButterflyV%nblk_loc))
         allocate (arr_acc_m(blocks%ButterflyU%nblk_loc))
         k1 = 0
         do i = 1, blocks%ButterflyV%nblk_loc
            arr_acc_n(i) = k1
            nn = size(blocks%ButterflyV%blocks(i)%matrix, 1)
            k1 = k1 + nn
         enddo

         k2 = 0
         do i = 1, blocks%ButterflyU%nblk_loc
            arr_acc_m(i) = k2
            mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
            k2 = k2 + mm
         enddo

         num_vectors = Nrnd

         if (BF_checkNAN(blocks)) then
            write (*, *) 'NAN in 0 BF_block_MVP_dat'
            stop
         end if

         if (chara == 'N') then

            if (isnanMat(random1(1:N,1:1),N,1)) then
               write (*, *) 'NAN in 1 BF_block_MVP_dat'
               stop
            end if

            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(0:level_butterfly + 2))

            allocate (BFvec%vec(0)%blocks(1, blocks%ButterflyV%nblk_loc))
            BFvec%vec(0)%num_row = 1
            BFvec%vec(0)%num_col = num_blocks
            BFvec%vec(0)%idx_r = 1
            BFvec%vec(0)%inc_r = 1
            BFvec%vec(0)%nr = 1
            BFvec%vec(0)%idx_c = blocks%ButterflyV%idx
            BFvec%vec(0)%inc_c = blocks%ButterflyV%inc
            BFvec%vec(0)%nc = blocks%ButterflyV%nblk_loc

            n5 = MPI_Wtime()
            do level = 0, level_half
               ! n1 = MPI_Wtime()
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

               BFvec%vec(level + 1)%idx_r = idx_r
               BFvec%vec(level + 1)%inc_r = inc_r
               BFvec%vec(level + 1)%nr = nr
               BFvec%vec(level + 1)%idx_c = idx_c
               BFvec%vec(level + 1)%inc_c = inc_c
               BFvec%vec(level + 1)%nc = nc

               if (nr > 0 .and. nc > 0) then
                  if (level /= level_butterfly + 1) then
                     BFvec%vec(level + 1)%num_row = 2**level
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
                  else
                     BFvec%vec(level + 1)%num_row = 2**level_butterfly
                     BFvec%vec(level + 1)%num_col = 1
                  endif
                  if (level_half /= level) then ! the last level doesn't require doubling block columns
                  if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                     BFvec%vec(level + 1)%nc = 2
                     BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

                  if (level == 0) then
                     flops = 0

                     n3 = MPI_Wtime()
                     group_count=blocks%ButterflyV%nblk_loc
                     allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_array='T'
                     transb_array='N'
                     alpha_array=BPACK_cone
                     beta_array=BPACK_cone
                     group_size=1

                     cnt=0
                     do j = 1, blocks%ButterflyV%nblk_loc
                        index_j = (j - 1)*inc_c + idx_c
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                        nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                        allocate (BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(rank, num_vectors))
                        BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix = 0
                        cnt=cnt+1
                        m_array(cnt)=rank
                        n_array(cnt)=num_vectors
                        k_array(cnt)=nn
                        lda_array(cnt)=nn
                        ldb_array(cnt)=ldi
                        ldc_array(cnt)=rank
                        a_array(cnt)=LOC(blocks%ButterflyV%blocks(j)%matrix(1,1))
                        b_array(cnt)=LOC(random1(1 + arr_acc_n(j), 1))
                        c_array(cnt)=LOC(BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(1,1))
                     enddo
                     call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     stats%Flop_Tmp = stats%Flop_Tmp + flops

                     deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3

                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                        index_j = idx_c
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                     endif

                  elseif (level == level_butterfly + 1) then
                     write(*,*)'should not arrive here'
                  else

                     n3 = MPI_Wtime()
                     group_count=nr*nc
                     allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_array='N'
                     transb_array='N'
                     alpha_array=BPACK_cone
                     beta_array=BPACK_cone
                     group_size=1

                     do j=1,2
                        cnt=0
                        do index_ij = 1, nr*nc
                           index_j_loc = (index_ij - 1)/nr + 1
                           index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c + idx_c

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                           index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix, 2)
                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                           if(.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix))then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           cnt=cnt+1
                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=nn1
                           lda_array(cnt)=mm
                           ldb_array(cnt)=nn1
                           ldc_array(cnt)=mm
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc+j-1)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        stats%Flop_Tmp = stats%Flop_Tmp + flops
                     enddo
                     deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif
               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo

               ! n2 = MPI_Wtime()
               ! time_tmp = time_tmp + n2-n1

               if (level_half /= level) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'B')
               endif
            enddo

            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

            if (level_half + 1 /= 0) then
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'R', 'C', 0)
            else
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'R', 'C', 0)
            endif

            n5 = MPI_Wtime()
            do level = level_half + 1, level_butterfly + 1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C')


               ! convert the local column-wise kernel block ranges to local row-wise output vector ranges
               if (level /= 0 .and. level /= level_butterfly + 1) then
                  idx_r = idx_r0*2 - 1
                  nr = nr0*2
                  inc_r = inc_r0
                  idx_c = ceiling_safe(idx_c0/2d0)
                  if (inc_c0 > 1) then
                     nc = nc0
                  else
                     nc = ceiling_safe(nc0/2d0)
                  endif
                  inc_c = ceiling_safe(inc_c0/2d0)
               else
                  idx_r = idx_r0
                  nr = nr0
                  inc_r = inc_r0
                  idx_c = idx_c0
                  nc = nc0
                  inc_c = inc_c0
               endif

               BFvec%vec(level + 1)%idx_r = idx_r
               BFvec%vec(level + 1)%inc_r = inc_r
               BFvec%vec(level + 1)%nr = nr
               BFvec%vec(level + 1)%idx_c = idx_c
               BFvec%vec(level + 1)%inc_c = inc_c
               BFvec%vec(level + 1)%nc = nc
               if (level /= level_butterfly + 1) then
                  BFvec%vec(level + 1)%num_row = 2**level
                  BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
               else
                  BFvec%vec(level + 1)%num_row = 2**level_butterfly
                  BFvec%vec(level + 1)%num_col = 1
               endif
               if (nr > 0 .and. nc > 0) then
                  allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

                  if (level == 0) then
                     write (*, *) 'should not arrive here'
                  elseif (level == level_butterfly + 1) then
                     flops = 0
                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        mm = size(blocks%ButterflyU%blocks(1)%matrix, 1)
                        rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                        allocate (matrixtemp(rank, num_vectors))
                        matrixtemp = 0
                        if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                           index_i = idx_r0
                           index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1
                           matrixtemp = BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix
                        endif
                        call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                        allocate (BFvec%vec(level + 1)%blocks(1, 1)%matrix(mm, num_vectors))
                        BFvec%vec(level + 1)%blocks(1, 1)%matrix = 0
                        call gemmf90(blocks%ButterflyU%blocks(1)%matrix, mm, matrixtemp, rank, random2(arr_acc_m(1)+1, 1), ldo, 'N', 'N', mm, num_vectors, rank, a, b, flop=flop)

                        flops = flops + flop
                        deallocate (matrixtemp)
                     else


                        n3 = MPI_Wtime()
                        group_count=nr0
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='N'
                        transb_array='N'
                        alpha_array=a
                        beta_array=b
                        group_size=1
                        cnt=0
                        do i = 1, nr0

                           index_i = (i - 1)*inc_r0 + idx_r0
                           index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1

                           rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                           mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                           allocate (BFvec%vec(level + 1)%blocks(i, 1)%matrix(mm, num_vectors))
                           BFvec%vec(level + 1)%blocks(i, 1)%matrix = 0

                           cnt=cnt+1

                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=rank
                           lda_array(cnt)=mm
                           ldb_array(cnt)=rank
                           ldc_array(cnt)=ldo
                           a_array(cnt)=LOC(blocks%ButterflyU%blocks(i)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix(1,1))
                           c_array(cnt)=LOC(random2(arr_acc_m(i)+1, 1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)

                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                        n4 = MPI_Wtime()
                        ! time_tmp = time_tmp + n4-n3

                     endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                  else
                     n3=MPI_Wtime()
                     if (nc0 > 1 .and. inc_c0 == 1) then
                        group_count=nr0*nc0
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='N'
                        transb_array='N'
                        alpha_array=BPACK_cone
                        beta_array=BPACK_cone
                        group_size=1


                        do jj=1,2
                           cnt=0
                           do index_ij = 1, nr0*nc0/2
                              index_j_loc0 = (index_ij - 1)/nr0 + 1

                              index_j_loc = 2*(index_j_loc0 - 1) + jj       !index_i_loc is local index of column-wise ordering at current level
                              index_i_loc = mod(index_ij - 1, nr0) + 1
                              index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                              index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                              index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                              index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                              index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                              index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                              index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                              index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                              index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                              if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                                 BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                              endif
                              cnt = cnt +1
                              m_array(cnt)=mm
                              n_array(cnt)=num_vectors
                              k_array(cnt)=nn
                              lda_array(cnt)=mm
                              ldb_array(cnt)=nn
                              ldc_array(cnt)=mm
                              a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                              b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                              c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))


                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                              if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(mm, num_vectors))
                                 BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix = 0
                              endif
                              cnt = cnt +1
                              m_array(cnt)=mm
                              n_array(cnt)=num_vectors
                              k_array(cnt)=nn
                              lda_array(cnt)=mm
                              ldb_array(cnt)=nn
                              ldc_array(cnt)=mm
                              a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                              b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                              c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                           enddo
                           call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                           stats%Flop_Tmp = stats%Flop_Tmp + flops
                        enddo
                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)

                     else

                        group_count=nr0*nc0*2
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='N'
                        transb_array='N'
                        alpha_array=BPACK_cone
                        beta_array=BPACK_cone
                        group_size=1

                        cnt=0
                        do index_ij = 1, nr0*nc0


                           index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                           index_i_loc = mod(index_ij - 1, nr0) + 1
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                           index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           ! !$omp critical
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           cnt = cnt +1
                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=nn
                           lda_array(cnt)=mm
                           ldb_array(cnt)=nn
                           ldc_array(cnt)=mm
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))


                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                           ! !$omp critical
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix = 0
                           endif
                           cnt = cnt +1
                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=nn
                           lda_array(cnt)=mm
                           ldb_array(cnt)=nn
                           ldc_array(cnt)=mm
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        stats%Flop_Tmp = stats%Flop_Tmp + flops
                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)

                     endif
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif

               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level /= level_butterfly + 1) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'R')
               endif
            enddo

            if (isnanMat(random2(1:M,1:1),M,1)) then
               write (*, *) ptree%MyID, 'NAN in 2 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if
            !deallocate (BFvec%vec)

            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5



         elseif (chara == 'T') then
            if (isnanMat(random1(1:M,1:1),M,1)) then
               write (*, *) 'NAN in 3 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if

            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(0:level_butterfly + 2))
            allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc, 1))
            BFvec%vec(0)%num_row = num_blocks
            BFvec%vec(0)%num_col = 1
            BFvec%vec(0)%idx_r = blocks%ButterflyU%idx
            BFvec%vec(0)%inc_r = blocks%ButterflyU%inc
            BFvec%vec(0)%nr = blocks%ButterflyU%nblk_loc
            BFvec%vec(0)%idx_c = 1
            BFvec%vec(0)%inc_c = 1
            BFvec%vec(0)%nc = 1

            n5 = MPI_Wtime()
            do level = level_butterfly + 1, level_half + 1, -1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

               BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
               BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
               BFvec%vec(level_butterfly - level + 2)%nr = nr
               BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
               BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
               BFvec%vec(level_butterfly - level + 2)%nc = nc

               if (nr > 0 .and. nc > 0) then

                  if (level /= 0) then
                     BFvec%vec(level_butterfly - level + 2)%num_row = 2**(level - 1)
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**(level_butterfly - level + 1)
                  else
                     BFvec%vec(level_butterfly - level + 2)%num_row = 1
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**level_butterfly
                  endif
                  if (level_half + 1 /= level) then ! the last level doesn't require doubling block rows
                  if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                     BFvec%vec(level_butterfly - level + 2)%nr = 2
                     BFvec%vec(level_butterfly - level + 2)%idx_r = BFvec%vec(level_butterfly - level + 2)%idx_r - 1 + mod(BFvec%vec(level_butterfly - level + 2)%idx_r, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

                  if (level == level_butterfly + 1) then
                     n3=MPI_Wtime()
                     group_count=blocks%ButterflyU%nblk_loc
                     allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_array='T'
                     transb_array='N'
                     alpha_array=BPACK_cone
                     beta_array=BPACK_cone
                     group_size=1
                     cnt=0
                     do i = 1, blocks%ButterflyU%nblk_loc

                        index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                        index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                        rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                        mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                        allocate (BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(rank, num_vectors))
                        BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix = 0
                        cnt=cnt+1
                        m_array(cnt)=rank
                        n_array(cnt)=num_vectors
                        k_array(cnt)=mm
                        lda_array(cnt)=mm
                        ldb_array(cnt)=ldi
                        ldc_array(cnt)=rank
                        a_array(cnt)=LOC(blocks%ButterflyU%blocks(i)%matrix(1,1))
                        b_array(cnt)=LOC(random1(1 + arr_acc_m(i), 1))
                        c_array(cnt)=LOC(BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(1,1))
                     enddo
                     call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                     deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                        index_i = blocks%ButterflyU%idx
                        index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1
                        call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                     endif

                  elseif (level == 0) then
                     write(*,*)'should not arrive here'
                  else
                     n3=MPI_Wtime()
                     group_count=nr*nc
                     allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_array='T'
                     transb_array='N'
                     alpha_array=BPACK_cone
                     beta_array=BPACK_cone
                     group_size=1

                     do i=1,2
                        cnt=0
                        do index_ij = 1, nr*nc

                           index_j_loc = (index_ij - 1)/nr + 1
                           index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c + idx_c

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                           index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k +i-1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if(.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix))then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           cnt=cnt+1
                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=mm1
                           lda_array(cnt)=mm1
                           ldb_array(cnt)=mm1
                           ldc_array(cnt)=nn
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+i-1, index_j_loc_k)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+i-1, index_jj_loc)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        stats%Flop_Tmp = stats%Flop_Tmp + flops
                     enddo
                     deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif
               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level_half + 1 /= level) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'B')
               endif
            enddo
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

            if (level_half /= level_butterfly + 1) then
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'C', 'R', 0)
            else
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'C', 'R', 0)
            endif

            n5 = MPI_Wtime()
            do level = level_half, 0, -1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'R')

               ! convert the local row-wise kernel block ranges to local column-wise output vector ranges
               if (level /= 0 .and. level /= level_butterfly + 1) then
                  idx_r = ceiling_safe(idx_r0/2d0)
                  if (inc_r0 > 1) then
                     nr = nr0
                  else
                     nr = ceiling_safe(nr0/2d0)
                  endif
                  inc_r = ceiling_safe(inc_r0/2d0)
                  idx_c = idx_c0*2 - 1
                  nc = nc0*2
                  inc_c = inc_c0
               else
                  idx_r = idx_r0
                  nr = nr0
                  inc_r = inc_r0
                  idx_c = idx_c0
                  nc = nc0
                  inc_c = inc_c0
               endif

               BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
               BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
               BFvec%vec(level_butterfly - level + 2)%nr = nr
               BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
               BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
               BFvec%vec(level_butterfly - level + 2)%nc = nc

               if (nr > 0 .and. nc > 0) then
                  if (level /= 0) then
                     BFvec%vec(level + 1)%num_row = 2**(level - 1)
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level + 1)
                  else
                     BFvec%vec(level + 1)%num_row = 1
                     BFvec%vec(level + 1)%num_col = 2**level_butterfly
                  endif

                  allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

                  if (level == level_butterfly + 1) then
                     write(*,*)'should not arrive here'

                  elseif (level == 0) then
                     flops = 0

                     call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                     if (ptree%pgrp(pgno_sub)%nproc > 1) then
                        nn = size(blocks%ButterflyV%blocks(1)%matrix, 1)
                        rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                        allocate (matrixtemp(rank, num_vectors))
                        matrixtemp = 0
                        if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                           index_j = blocks%ButterflyV%idx
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1

                           matrixtemp = BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix
                        endif
                        call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                        allocate (BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix(nn, num_vectors))
                        BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix = 0
                        call gemmf90(blocks%ButterflyV%blocks(1)%matrix, nn, matrixtemp, rank, random2(arr_acc_n(1)+1, 1), ldo, 'N', 'N', nn, num_vectors, rank, a, b, flop=flop)
                        flops = flops + flop
                        deallocate (matrixtemp)
                     else
                        n3=MPI_Wtime()
                        group_count=blocks%ButterflyV%nblk_loc
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='N'
                        transb_array='N'
                        alpha_array=a
                        beta_array=b
                        group_size=1
                        cnt=0
                        do j = 1, blocks%ButterflyV%nblk_loc
                           index_j = (j - 1)*blocks%ButterflyV%inc + blocks%ButterflyV%idx
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1
                           nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                           rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                           allocate (BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix(nn, num_vectors))
                           BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix = 0
                           cnt=cnt+1

                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=rank
                           lda_array(cnt)=nn
                           ldb_array(cnt)=rank
                           ldc_array(cnt)=ldo
                           a_array(cnt)=LOC(blocks%ButterflyV%blocks(j)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix(1,1))
                           c_array(cnt)=LOC(random2(arr_acc_n(j)+1, 1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                        n4 = MPI_Wtime()
                        ! time_tmp = time_tmp + n4-n3
                     endif
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                  else
                     n3=MPI_Wtime()
                     if (nr0 > 1 .and. inc_r0 == 1) then ! this special treatment makes sure two threads do not write to the same address simultaneously
                        group_count=nr0*nc0
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='T'
                        transb_array='N'
                        alpha_array=BPACK_cone
                        beta_array=BPACK_cone
                        group_size=1

                        do ii = 1, 2
                           cnt=0
                           do index_ij = 1, nr0*nc0/2
                              index_i_loc0 = (index_ij - 1)/nc0 + 1
                              index_i_loc = (index_i_loc0 - 1)*2 + ii
                              index_j_loc = mod(index_ij - 1, nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
                              index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                              index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                              index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                              index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                              index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                              index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                              index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                              index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                              index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                              if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                                 allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                                 BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                              endif

                              cnt = cnt +1
                              m_array(cnt)=nn
                              n_array(cnt)=num_vectors
                              k_array(cnt)=mm
                              lda_array(cnt)=mm
                              ldb_array(cnt)=mm
                              ldc_array(cnt)=nn
                              a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                              b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                              c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))

                              mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                              nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                              if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                                 allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                                 BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                              endif

                              cnt = cnt +1
                              m_array(cnt)=nn
                              n_array(cnt)=num_vectors
                              k_array(cnt)=mm
                              lda_array(cnt)=mm
                              ldb_array(cnt)=mm
                              ldc_array(cnt)=nn
                              a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                              b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                              c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))
                           enddo
                           call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                           stats%Flop_Tmp = stats%Flop_Tmp + flops
                        enddo
                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)

                     else
                        group_count=nr0*nc0*2
                        allocate(transa_array(group_count),transb_array(group_count),alpha_array(group_count),beta_array(group_count),group_size(group_count),m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_array='T'
                        transb_array='N'
                        alpha_array=BPACK_cone
                        beta_array=BPACK_cone
                        group_size=1

                        cnt=0
                        do index_ij = 1, nr0*nc0
                           index_j_loc = (index_ij - 1)/nr0 + 1
                           index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                           index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           cnt = cnt + 1
                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=mm
                           lda_array(cnt)=mm
                           ldb_array(cnt)=mm
                           ldc_array(cnt)=nn
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))


                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                           endif
                           cnt = cnt + 1
                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=mm
                           lda_array(cnt)=mm
                           ldb_array(cnt)=mm
                           ldc_array(cnt)=nn
                           a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                           b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))
                        enddo
                        call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        stats%Flop_Tmp = stats%Flop_Tmp + flops
                        deallocate(transa_array,transb_array,alpha_array,beta_array,group_size,m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array)
                     endif
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
               endif

               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo

               if (level /= 0) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'R')
               endif
            enddo
            n6 = MPI_Wtime()
            time_tmp1 = time_tmp1 + n6-n5

            if (isnanMat(random2(1:N,1:1),N,1)) then
               write (*, *) 'NAN in 4 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
               stop
            end if

         endif

         do level = 0, level_butterfly + 2
            do j = 1, BFvec%vec(level)%nc
               do i = 1, BFvec%vec(level)%nr
                  if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
               enddo
            enddo
            if (allocated(BFvec%vec(level)%blocks)) deallocate (BFvec%vec(level)%blocks)
         enddo

         deallocate (BFvec%vec)
         deallocate (arr_acc_m, arr_acc_n)


      endif

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2-n1
      time_tmp2 = time_tmp2 + n2-n1
      return

   end subroutine BF_block_MVP_dat_batch_mkl

#endif


#ifdef HAVE_MAGMA
subroutine BF_block_MVP_dat_batch_magma(blocks, chara, M, N, Nrnd, random1, ldi, random2, ldo, a, b, ptree, stats)


   implicit none

   integer M, N, Nrnd, index_i, index_j, na, nb, index_start, num_vectors
   integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
   integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
   integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
   integer vector_a, vector_b, nn1, nn2, mm1, mm2, level_half, level_final, pgno_sub
   integer idx_r, inc_r, nr, idx_c, inc_c, nc
   integer idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0
   DT ctemp, a, b
   character chara
   type(matrixblock)::blocks
   type(proctree)::ptree
   integer pgno, comm, ierr
   type(Hstat)::stats
   real(kind=8)::flop, flops,n1,n2,n3,n4,n5,n6,n7,n8
   integer index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc, index_i_loc_s, index_i_loc_k, index_j_loc, index_j_loc0, index_i_loc0, index_j_loc_s, index_j_loc_k
   integer*8:: cnta,cntb,cntc
   integer,target,allocatable::group_size(:),m_array(:),n_array(:),k_array(:),lda_array(:),ldb_array(:),ldc_array(:)
   type(c_ptr),target,allocatable::a_array(:),b_array(:),c_array(:)
   integer group_count,cnt
   type(c_ptr) :: dm_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dn_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dk_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dlda_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dldb_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dldc_array  !! on GPU, int*     (array of batchcount integers)
   type(c_ptr) :: dA_array,dB_array,dC_array  !! on GPU, double** (array of batchcount DT* pointers)
   DT,target,allocatable::AA(:),BB(:),CC(:)
   integer :: asize,bsize,csize

   type(butterfly_vec) :: BFvec
   integer ldi,ldo
   DT,target :: random1(ldi, *), random2(ldo, *)
   DT, pointer :: random1_p(:, :),random2_p(:, :)
   DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), Vout_tmp(:, :)

   integer, allocatable:: arr_acc_m(:), arr_acc_n(:)

   type(c_ptr) :: dA,dB,dC           !! on GPU, double*  (array of n * n * batchcount doubles)
   integer info
   integer :: dev
   type(c_ptr) :: queue  !! magma_queue_t
   integer(c_int):: transa_magma,transb_magma
   DT:: alpha_magma,beta_magma

   n1 = MPI_Wtime()

   call MAGMA_init()
   dev = 0
   call MAGMA_queue_create( dev, queue )


   level_butterfly = blocks%level_butterfly
   pgno = blocks%pgno
   comm = ptree%pgrp(pgno)%comm
   if (comm == MPI_COMM_NULL) then
      write (*, *) 'ninin', pgno, comm == MPI_COMM_NULL, ptree%MyID
   endif

   call assert(IOwnPgrp(ptree, pgno), 'I do not share this block!')

   if (blocks%style == 1) then
      call Full_block_MVP_dat(blocks, chara, M, Nrnd, random1, ldi, random2, ldo, a, b)
      return
   endif

   if (level_butterfly == 0) then
      rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
      call assert(rank > 0, 'rank incorrect in blocks%ButterflyU')
      allocate (matrixtemp(rank, Nrnd))
      matrixtemp = 0
      allocate (matrixtemp1(rank, Nrnd))
      matrixtemp1 = 0
      ! for implementation simplicity, MPI_ALLREDUCE is used even when nproc==1
      if (chara == 'N') then !Vout=U*V^T*Vin
         allocate (Vout_tmp(M,Nrnd))
         Vout_tmp = 0
         call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyV%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 2')
         call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
         call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, M, 'N', 'N', size(blocks%ButterflyU%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         random2(1:M,1:Nrnd) = b*random2(1:M,1:Nrnd) + a*Vout_tmp
      else if (chara == 'T') then !Vout=V*U^T*Vin
         allocate (Vout_tmp(N,Nrnd))
         Vout_tmp = 0
         call gemmf90(blocks%ButterflyU%blocks(1)%matrix, size(blocks%ButterflyU%blocks(1)%matrix, 1), random1, ldi, matrixtemp, rank, 'T', 'N', rank, Nrnd, size(blocks%ButterflyU%blocks(1)%matrix, 1), BPACK_cone, BPACK_czero, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         call assert(MPI_COMM_NULL /= comm, 'communicator should not be null 3')
         call MPI_ALLREDUCE(matrixtemp, matrixtemp1, rank*Nrnd, MPI_DT, MPI_SUM, comm, ierr)
         call gemmf90(blocks%ButterflyV%blocks(1)%matrix, size(blocks%ButterflyV%blocks(1)%matrix, 1), matrixtemp1, rank, Vout_tmp, N, 'N', 'N', size(blocks%ButterflyV%blocks(1)%matrix, 1), Nrnd, rank, BPACK_cone, BPACK_czero, flop)
         stats%Flop_Tmp = stats%Flop_Tmp + flop
         random2(1:N,1:Nrnd) = b*random2(1:N,1:Nrnd) + a*Vout_tmp
      endif

      deallocate (matrixtemp)
      deallocate (matrixtemp1)
      deallocate (Vout_tmp)

   else

      allocate (arr_acc_n(blocks%ButterflyV%nblk_loc))
      allocate (arr_acc_m(blocks%ButterflyU%nblk_loc))
      k1 = 0
      do i = 1, blocks%ButterflyV%nblk_loc
         arr_acc_n(i) = k1
         nn = size(blocks%ButterflyV%blocks(i)%matrix, 1)
         k1 = k1 + nn
      enddo

      k2 = 0
      do i = 1, blocks%ButterflyU%nblk_loc
         arr_acc_m(i) = k2
         mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
         k2 = k2 + mm
      enddo

      num_vectors = Nrnd

      if (BF_checkNAN(blocks)) then
         write (*, *) 'NAN in 0 BF_block_MVP_dat'
         stop
      end if

      if (chara == 'N') then

         if (isnanMat(random1(1:N,1:1),N,1)) then
            write (*, *) 'NAN in 1 BF_block_MVP_dat'
            stop
         end if

         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         level_half = blocks%level_half

         allocate (BFvec%vec(0:level_butterfly + 2))

         allocate (BFvec%vec(0)%blocks(1, blocks%ButterflyV%nblk_loc))
         BFvec%vec(0)%num_row = 1
         BFvec%vec(0)%num_col = num_blocks
         BFvec%vec(0)%idx_r = 1
         BFvec%vec(0)%inc_r = 1
         BFvec%vec(0)%nr = 1
         BFvec%vec(0)%idx_c = blocks%ButterflyV%idx
         BFvec%vec(0)%inc_c = blocks%ButterflyV%inc
         BFvec%vec(0)%nc = blocks%ButterflyV%nblk_loc

         n5 = MPI_Wtime()
         do level = 0, level_half
            ! n1 = MPI_Wtime()
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

            BFvec%vec(level + 1)%idx_r = idx_r
            BFvec%vec(level + 1)%inc_r = inc_r
            BFvec%vec(level + 1)%nr = nr
            BFvec%vec(level + 1)%idx_c = idx_c
            BFvec%vec(level + 1)%inc_c = inc_c
            BFvec%vec(level + 1)%nc = nc

            if (nr > 0 .and. nc > 0) then
               if (level /= level_butterfly + 1) then
                  BFvec%vec(level + 1)%num_row = 2**level
                  BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
               else
                  BFvec%vec(level + 1)%num_row = 2**level_butterfly
                  BFvec%vec(level + 1)%num_col = 1
               endif
               if (level_half /= level) then ! the last level doesn't require doubling block columns
               if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                  BFvec%vec(level + 1)%nc = 2
                  BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
               endif
               endif
               allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

               if (level == 0) then
                  flops = 0
                  call magmaf_wtime(n7)

                  n3 = MPI_Wtime()
                  group_count=blocks%ButterflyV%nblk_loc
                  allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                  transa_magma=MagmaTrans
                  transb_magma=MagmaNoTrans
                  alpha_magma=BPACK_cone
                  beta_magma=BPACK_cone
                  cnt=0
                  asize=0
                  bsize=0
                  csize=0
                  do j = 1, blocks%ButterflyV%nblk_loc
                     index_j = (j - 1)*inc_c + idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                     nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                     allocate (BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(rank, num_vectors))
                     BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix = 0
                     cnt=cnt+1
                     m_array(cnt)=rank
                     n_array(cnt)=num_vectors
                     k_array(cnt)=nn
                     lda_array(cnt)=nn
                     ldb_array(cnt)=nn
                     ldc_array(cnt)=rank
                     ! a_array(cnt)=LOC(blocks%ButterflyV%blocks(j)%matrix(1,1))
                     ! b_array(cnt)=LOC(random1(1 + arr_acc_n(j), 1))
                     ! c_array(cnt)=LOC(BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(1,1))
                     asize=asize+lda_array(cnt)*m_array(cnt)
                     bsize=bsize+ldb_array(cnt)*n_array(cnt)
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo
                  allocate(AA(asize))
                  allocate(BB(bsize))
                  allocate(CC(csize))
                  info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                  info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                  info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                  info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                  info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                  info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                  info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                  cnt=0
                  asize=0
                  bsize=0
                  csize=0
                  do j = 1, blocks%ButterflyV%nblk_loc
                     index_j = (j - 1)*inc_c + idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                     nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                     BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix = 0
                     cnt=cnt+1

                     call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix,ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                     call copymatf77(' ',ldb_array(cnt),n_array(cnt),random1(1 + arr_acc_n(j), 1),ldi,BB(bsize+1),ldb_array(cnt))
                     call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyV%blocks(j)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                     a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                     b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                     c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                     asize=asize+lda_array(cnt)*m_array(cnt)
                     bsize=bsize+ldb_array(cnt)*n_array(cnt)
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo

                  call magmaf_wtime(n8)
                  time_tmp5 = time_tmp5 + n8-n7

                  call magmaf_wtime(n7)
                  call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                  call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                  call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                  call MAGMA_gemm_vbatched( &
                  transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                  dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)


                  ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                  ! stats%Flop_Tmp = stats%Flop_Tmp + flops

                  call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )

                  call magmaf_wtime(n8)
                  time_tmp5 = time_tmp5 + n8-n7

                  cnt=0
                  csize=0
                  do j = 1, blocks%ButterflyV%nblk_loc
                     index_j = (j - 1)*inc_c + idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     cnt=cnt+1
                     call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix,ldc_array(cnt))
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo

                  deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                  info = MAGMA_free( dA )
                  info = MAGMA_free( dB )
                  info = MAGMA_free( dC )
                  info = MAGMA_free( dA_array )
                  info = MAGMA_free( dB_array )
                  info = MAGMA_free( dC_array )
                  info = MAGMA_free( dm_array )
                  info = MAGMA_free( dn_array )
                  info = MAGMA_free( dk_array )
                  info = MAGMA_free( dlda_array )
                  info = MAGMA_free( dldb_array )
                  info = MAGMA_free( dldc_array )

                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3

                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                     index_j = idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                  endif

               elseif (level == level_butterfly + 1) then
                  write(*,*)'should not arrive here'
               else

                  n3 = MPI_Wtime()
                  do j=1,2
                     group_count=nr*nc

                     call magmaf_wtime(n7)

                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaNoTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=BPACK_cone
                     beta_magma=BPACK_cone

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                        index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix, 2)
                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                        if(.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix))then
                           allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                           BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                        endif
                        cnt=cnt+1
                        m_array(cnt)=mm
                        n_array(cnt)=num_vectors
                        k_array(cnt)=nn1
                        lda_array(cnt)=mm
                        ldb_array(cnt)=nn1
                        ldc_array(cnt)=mm

                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc+j-1)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                     enddo

                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                        index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix, 2)
                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc+j-1)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7

                     call magmaf_wtime(n7)

                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)


                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     ! stats%Flop_Tmp = stats%Flop_Tmp + flops

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7

                     cnt=0
                     csize=0
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                        index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+j-1)%matrix, 2)
                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))

                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )

                  enddo
                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3
               endif
            endif
            do j = 1, BFvec%vec(level)%nc
               do i = 1, BFvec%vec(level)%nr
                  if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
               enddo
            enddo

            ! n2 = MPI_Wtime()
            ! time_tmp = time_tmp + n2-n1

            if (level_half /= level) then
               call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'B')
            endif
         enddo

         n6 = MPI_Wtime()
         time_tmp1 = time_tmp1 + n6-n5

         if (level_half + 1 /= 0) then
            call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'R', 'C', 0)
         else
            call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'R', 'C', 0)
         endif

         n5 = MPI_Wtime()
         do level = level_half + 1, level_butterfly + 1
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C')


            ! convert the local column-wise kernel block ranges to local row-wise output vector ranges
            if (level /= 0 .and. level /= level_butterfly + 1) then
               idx_r = idx_r0*2 - 1
               nr = nr0*2
               inc_r = inc_r0
               idx_c = ceiling_safe(idx_c0/2d0)
               if (inc_c0 > 1) then
                  nc = nc0
               else
                  nc = ceiling_safe(nc0/2d0)
               endif
               inc_c = ceiling_safe(inc_c0/2d0)
            else
               idx_r = idx_r0
               nr = nr0
               inc_r = inc_r0
               idx_c = idx_c0
               nc = nc0
               inc_c = inc_c0
            endif

            BFvec%vec(level + 1)%idx_r = idx_r
            BFvec%vec(level + 1)%inc_r = inc_r
            BFvec%vec(level + 1)%nr = nr
            BFvec%vec(level + 1)%idx_c = idx_c
            BFvec%vec(level + 1)%inc_c = inc_c
            BFvec%vec(level + 1)%nc = nc

            if (nr > 0 .and. nc > 0) then

               if (level /= level_butterfly + 1) then
                  BFvec%vec(level + 1)%num_row = 2**level
                  BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
               else
                  BFvec%vec(level + 1)%num_row = 2**level_butterfly
                  BFvec%vec(level + 1)%num_col = 1
               endif

               allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

               if (level == 0) then
                  write (*, *) 'should not arrive here'
               elseif (level == level_butterfly + 1) then
                  flops = 0
                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     mm = size(blocks%ButterflyU%blocks(1)%matrix, 1)
                     rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                     allocate (matrixtemp(rank, num_vectors))
                     matrixtemp = 0
                     if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                        index_i = idx_r0
                        index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1
                        matrixtemp = BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix
                     endif
                     call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                     allocate (BFvec%vec(level + 1)%blocks(1, 1)%matrix(mm, num_vectors))
                     BFvec%vec(level + 1)%blocks(1, 1)%matrix = 0
                     call gemmf90(blocks%ButterflyU%blocks(1)%matrix, mm, matrixtemp, rank, random2(arr_acc_m(1)+1, 1), ldo, 'N', 'N', mm, num_vectors, rank, a, b, flop=flop)

                     flops = flops + flop
                     deallocate (matrixtemp)
                  else

                     call magmaf_wtime(n7)
                     n3 = MPI_Wtime()
                     group_count=nr0

                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaNoTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=a
                     beta_magma=b

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     cnt=0
                     do i = 1, nr0

                        index_i = (i - 1)*inc_r0 + idx_r0
                        index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1

                        rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                        mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                        allocate (BFvec%vec(level + 1)%blocks(i, 1)%matrix(mm, num_vectors))
                        BFvec%vec(level + 1)%blocks(i, 1)%matrix = 0

                        cnt=cnt+1

                        m_array(cnt)=mm
                        n_array(cnt)=num_vectors
                        k_array(cnt)=rank
                        lda_array(cnt)=mm
                        ldb_array(cnt)=rank
                        ldc_array(cnt)=mm

                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     cnt=0
                     do i = 1, nr0

                        index_i = (i - 1)*inc_r0 + idx_r0
                        index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1

                        rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                        mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)


                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),random2(arr_acc_m(i)+1, 1),ldo,CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyU%blocks(i)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                        ! a_array(cnt)=LOC(blocks%ButterflyU%blocks(i)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_i_loc_s, 1)%matrix(1,1))
                        ! c_array(cnt)=LOC(random2(arr_acc_m(i)+1, 1))
                     enddo
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7


                     call magmaf_wtime(n7)
                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7
                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)

                     csize=0
                     cnt=0
                     do i = 1, nr0

                        index_i = (i - 1)*inc_r0 + idx_r0
                        index_i_loc_s = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),random2(arr_acc_m(i)+1, 1),ldo)

                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )

                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3

                  endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops
               else
                  n3=MPI_Wtime()
                  if (nc0 > 1 .and. inc_c0 == 1) then
                     group_count=nr0*nc0

                     do jj=1,2
                        call magmaf_wtime(n7)
                        allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_magma=MagmaNoTrans
                        transb_magma=MagmaNoTrans
                        alpha_magma=BPACK_cone
                        beta_magma=BPACK_cone

                        cnt=0
                        asize=0
                        bsize=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_j_loc0 = (index_ij - 1)/nr0 + 1

                           index_j_loc = 2*(index_j_loc0 - 1) + jj       !index_i_loc is local index of column-wise ordering at current level
                           index_i_loc = mod(index_ij - 1, nr0) + 1
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                           index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif
                           cnt = cnt +1
                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=nn
                           lda_array(cnt)=mm
                           ldb_array(cnt)=nn
                           ldc_array(cnt)=mm


                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(mm, num_vectors))
                              BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix = 0
                           endif
                           cnt = cnt +1
                           m_array(cnt)=mm
                           n_array(cnt)=num_vectors
                           k_array(cnt)=nn
                           lda_array(cnt)=mm
                           ldb_array(cnt)=nn
                           ldc_array(cnt)=mm

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)


                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                        enddo

                        allocate(AA(asize))
                        allocate(BB(bsize))
                        allocate(CC(csize))
                        info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                        info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                        info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                        info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                        info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                        info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                        info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                        cnt=0
                        asize=0
                        bsize=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_j_loc0 = (index_ij - 1)/nr0 + 1

                           index_j_loc = 2*(index_j_loc0 - 1) + jj       !index_i_loc is local index of column-wise ordering at current level
                           index_i_loc = mod(index_ij - 1, nr0) + 1
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                           index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           cnt = cnt +1

                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                           call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                           call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                           a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                           b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                           c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                           cnt = cnt +1

                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                           call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                           call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                           a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                           b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                           c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)


                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                        enddo
                        call magmaf_wtime(n8)
                        time_tmp5 = time_tmp5 + n8-n7

                        ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        ! stats%Flop_Tmp = stats%Flop_Tmp + flops
                        call magmaf_wtime(n7)
                        call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                        call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                        call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                        call MAGMA_gemm_vbatched( &
                        transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                        dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                        call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                        call magmaf_wtime(n8)
                        time_tmp5 = time_tmp5 + n8-n7
                        cnt=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_j_loc0 = (index_ij - 1)/nr0 + 1

                           index_j_loc = 2*(index_j_loc0 - 1) + jj       !index_i_loc is local index of column-wise ordering at current level
                           index_i_loc = mod(index_ij - 1, nr0) + 1
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                           index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                           index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                           index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           cnt = cnt +1

                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))
                           c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                           cnt = cnt +1

                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1),ldc_array(cnt))

                           asize=asize+lda_array(cnt)*k_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)
                        enddo

                        deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                        info = MAGMA_free( dA )
                        info = MAGMA_free( dB )
                        info = MAGMA_free( dC )
                        info = MAGMA_free( dA_array )
                        info = MAGMA_free( dB_array )
                        info = MAGMA_free( dC_array )
                        info = MAGMA_free( dm_array )
                        info = MAGMA_free( dn_array )
                        info = MAGMA_free( dk_array )
                        info = MAGMA_free( dlda_array )
                        info = MAGMA_free( dldb_array )
                        info = MAGMA_free( dldc_array )

                     enddo


                  else

                     group_count=nr0*nc0*2
                     call magmaf_wtime(n7)

                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaNoTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=BPACK_cone
                     beta_magma=BPACK_cone

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                        index_i_loc = mod(index_ij - 1, nr0) + 1
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                        index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        ! !$omp critical
                        if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                           allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                           BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                        endif
                        cnt = cnt +1
                        m_array(cnt)=mm
                        n_array(cnt)=num_vectors
                        k_array(cnt)=nn
                        lda_array(cnt)=mm
                        ldb_array(cnt)=nn
                        ldc_array(cnt)=mm

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                        ! !$omp critical
                        if (.not. associated(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix)) then
                           allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(mm, num_vectors))
                           BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix = 0
                        endif
                        cnt = cnt +1
                        m_array(cnt)=mm
                        n_array(cnt)=num_vectors
                        k_array(cnt)=nn
                        lda_array(cnt)=mm
                        ldb_array(cnt)=nn
                        ldc_array(cnt)=mm

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                        index_i_loc = mod(index_ij - 1, nr0) + 1
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                        index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)

                        cnt = cnt +1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                        cnt = cnt +1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7


                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     ! stats%Flop_Tmp = stats%Flop_Tmp + flops

                     call magmaf_wtime(n7)
                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7

                     cnt=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                        index_i_loc = mod(index_ij - 1, nr0) + 1
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                        index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2)  !index_ii is global index in BFvec%vec(level+1)

                        index_i_loc_s = (index_ii - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)

                        cnt = cnt +1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))

                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                        cnt = cnt +1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level + 1)%blocks(index_i_loc_s+1, index_j_loc_s)%matrix(1,1),ldc_array(cnt))

                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )

                  endif
                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3
               endif
            endif

            do j = 1, BFvec%vec(level)%nc
               do i = 1, BFvec%vec(level)%nr
                  if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
               enddo
            enddo

            if (level /= level_butterfly + 1) then
               call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'R')
            endif
         enddo

         if (isnanMat(random2(1:M,1:1),M,1)) then
            write (*, *) ptree%MyID, 'NAN in 2 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
            stop
         end if
         !deallocate (BFvec%vec)

         n6 = MPI_Wtime()
         time_tmp1 = time_tmp1 + n6-n5



      elseif (chara == 'T') then

         if (isnanMat(random1(1:M,1:1),M,1)) then
            write (*, *) 'NAN in 3 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
            stop
         end if

         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         level_half = blocks%level_half

         allocate (BFvec%vec(0:level_butterfly + 2))
         allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc, 1))
         BFvec%vec(0)%num_row = num_blocks
         BFvec%vec(0)%num_col = 1
         BFvec%vec(0)%idx_r = blocks%ButterflyU%idx
         BFvec%vec(0)%inc_r = blocks%ButterflyU%inc
         BFvec%vec(0)%nr = blocks%ButterflyU%nblk_loc
         BFvec%vec(0)%idx_c = 1
         BFvec%vec(0)%inc_c = 1
         BFvec%vec(0)%nc = 1

         n5 = MPI_Wtime()
         do level = level_butterfly + 1, level_half + 1, -1
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

            BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
            BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
            BFvec%vec(level_butterfly - level + 2)%nr = nr
            BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
            BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
            BFvec%vec(level_butterfly - level + 2)%nc = nc

            if (nr > 0 .and. nc > 0) then

               if (level /= 0) then
                  BFvec%vec(level_butterfly - level + 2)%num_row = 2**(level - 1)
                  BFvec%vec(level_butterfly - level + 2)%num_col = 2**(level_butterfly - level + 1)
               else
                  BFvec%vec(level_butterfly - level + 2)%num_row = 1
                  BFvec%vec(level_butterfly - level + 2)%num_col = 2**level_butterfly
               endif
               if (level_half + 1 /= level) then ! the last level doesn't require doubling block rows
               if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                  BFvec%vec(level_butterfly - level + 2)%nr = 2
                  BFvec%vec(level_butterfly - level + 2)%idx_r = BFvec%vec(level_butterfly - level + 2)%idx_r - 1 + mod(BFvec%vec(level_butterfly - level + 2)%idx_r, 2)
               endif
               endif
               allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

               if (level == level_butterfly + 1) then
                  n3=MPI_Wtime()
                  group_count=blocks%ButterflyU%nblk_loc

                  call magmaf_wtime(n7)

                  allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                  transa_magma=MagmaTrans
                  transb_magma=MagmaNoTrans
                  alpha_magma=BPACK_cone
                  beta_magma=BPACK_cone

                  cnt=0
                  asize=0
                  bsize=0
                  csize=0
                  do i = 1, blocks%ButterflyU%nblk_loc

                     index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                     rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                     mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                     allocate (BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(rank, num_vectors))
                     BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix = 0
                     cnt=cnt+1
                     m_array(cnt)=rank
                     n_array(cnt)=num_vectors
                     k_array(cnt)=mm
                     lda_array(cnt)=mm
                     ldb_array(cnt)=mm
                     ldc_array(cnt)=rank
                     ! a_array(cnt)=LOC(blocks%ButterflyU%blocks(i)%matrix(1,1))
                     ! b_array(cnt)=LOC(random1(1 + arr_acc_m(i), 1))
                     ! c_array(cnt)=LOC(BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(1,1))
                     asize=asize+lda_array(cnt)*m_array(cnt)
                     bsize=bsize+ldb_array(cnt)*n_array(cnt)
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo

                  allocate(AA(asize))
                  allocate(BB(bsize))
                  allocate(CC(csize))
                  info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                  info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                  info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                  info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                  info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                  info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                  info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                  info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                  cnt=0
                  asize=0
                  bsize=0
                  csize=0
                  do i = 1, blocks%ButterflyU%nblk_loc

                     index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                     rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                     mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                     cnt=cnt+1

                     call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                     call copymatf77(' ',ldb_array(cnt),n_array(cnt),random1(1 + arr_acc_m(i), 1),ldi,BB(bsize+1),ldb_array(cnt))
                     call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyU%blocks(i)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                     a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                     b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                     c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                     ! a_array(cnt)=LOC(blocks%ButterflyU%blocks(i)%matrix(1,1))
                     ! b_array(cnt)=LOC(random1(1 + arr_acc_m(i), 1))
                     ! c_array(cnt)=LOC(BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(1,1))
                     asize=asize+lda_array(cnt)*m_array(cnt)
                     bsize=bsize+ldb_array(cnt)*n_array(cnt)
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo
                  call magmaf_wtime(n8)
                  time_tmp5 = time_tmp5 + n8-n7

                  ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                  ! stats%Flop_Tmp = stats%Flop_Tmp + flops

                  call magmaf_wtime(n7)

                  call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                  call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                  call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                  call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                  call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                  call MAGMA_gemm_vbatched( &
                  transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                  dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                  call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                  call magmaf_wtime(n8)
                  time_tmp5 = time_tmp5 + n8-n7
                  cnt=0
                  csize=0
                  do i = 1, blocks%ButterflyU%nblk_loc

                     index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                     rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                     mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                     cnt=cnt+1

                     call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(1,1),ldc_array(cnt))
                     csize=csize+ldc_array(cnt)*n_array(cnt)
                  enddo
                  deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                  info = MAGMA_free( dA )
                  info = MAGMA_free( dB )
                  info = MAGMA_free( dC )
                  info = MAGMA_free( dA_array )
                  info = MAGMA_free( dB_array )
                  info = MAGMA_free( dC_array )
                  info = MAGMA_free( dm_array )
                  info = MAGMA_free( dn_array )
                  info = MAGMA_free( dk_array )
                  info = MAGMA_free( dlda_array )
                  info = MAGMA_free( dldb_array )
                  info = MAGMA_free( dldc_array )


                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3
                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                     index_i = blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1
                     call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                  endif

               elseif (level == 0) then
                  write(*,*)'should not arrive here'
               else
                  n3=MPI_Wtime()
                  group_count=nr*nc
                  do i=1,2
                     call magmaf_wtime(n7)

                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=BPACK_cone
                     beta_magma=BPACK_cone

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr*nc

                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                        index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k +i-1, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        if(.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix))then
                           allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                           BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                        endif
                        cnt=cnt+1
                        m_array(cnt)=nn
                        n_array(cnt)=num_vectors
                        k_array(cnt)=mm1
                        lda_array(cnt)=mm1
                        ldb_array(cnt)=mm1
                        ldc_array(cnt)=nn
                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+i-1, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+i-1, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))

                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr*nc

                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                        index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+i-1, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k+i-1, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+i-1, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+i-1, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))

                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7


                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     ! stats%Flop_Tmp = stats%Flop_Tmp + flops
                     call magmaf_wtime(n7)
                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7
                     cnt=0
                     csize=0
                     do index_ij = 1, nr*nc

                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                        index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )
                  enddo
                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3
               endif
            endif
            do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
               do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                  if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
               enddo
            enddo

            if (level_half + 1 /= level) then
               call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'B')
            endif
         enddo
         n6 = MPI_Wtime()
         time_tmp1 = time_tmp1 + n6-n5

         if (level_half /= level_butterfly + 1) then
            call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'C', 'R', 0)
         else
            call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'C', 'R', 0)
         endif

         n5 = MPI_Wtime()
         do level = level_half, 0, -1
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'R')

            ! convert the local row-wise kernel block ranges to local column-wise output vector ranges
            if (level /= 0 .and. level /= level_butterfly + 1) then
               idx_r = ceiling_safe(idx_r0/2d0)
               if (inc_r0 > 1) then
                  nr = nr0
               else
                  nr = ceiling_safe(nr0/2d0)
               endif
               inc_r = ceiling_safe(inc_r0/2d0)
               idx_c = idx_c0*2 - 1
               nc = nc0*2
               inc_c = inc_c0
            else
               idx_r = idx_r0
               nr = nr0
               inc_r = inc_r0
               idx_c = idx_c0
               nc = nc0
               inc_c = inc_c0
            endif

            BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
            BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
            BFvec%vec(level_butterfly - level + 2)%nr = nr
            BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
            BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
            BFvec%vec(level_butterfly - level + 2)%nc = nc

            if (nr > 0 .and. nc > 0) then

               if (level /= 0) then
                  BFvec%vec(level + 1)%num_row = 2**(level - 1)
                  BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level + 1)
               else
                  BFvec%vec(level + 1)%num_row = 1
                  BFvec%vec(level + 1)%num_col = 2**level_butterfly
               endif

               allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

               if (level == level_butterfly + 1) then
                  write(*,*)'should not arrive here'

               elseif (level == 0) then
                  flops = 0

                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     nn = size(blocks%ButterflyV%blocks(1)%matrix, 1)
                     rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                     allocate (matrixtemp(rank, num_vectors))
                     matrixtemp = 0
                     if (ptree%MyID == ptree%pgrp(pgno_sub)%head) then
                        index_j = blocks%ButterflyV%idx
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1

                        matrixtemp = BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix
                     endif
                     call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                     allocate (BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix(nn, num_vectors))
                     BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix = 0
                     call gemmf90(blocks%ButterflyV%blocks(1)%matrix, nn, matrixtemp, rank, random2(arr_acc_n(1)+1, 1), ldo, 'N', 'N', nn, num_vectors, rank, a, b, flop=flop)
                     flops = flops + flop
                     deallocate (matrixtemp)
                  else
                     n3=MPI_Wtime()
                     group_count=blocks%ButterflyV%nblk_loc
                     call magmaf_wtime(n7)
                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaNoTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=a
                     beta_magma=b

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do j = 1, blocks%ButterflyV%nblk_loc
                        index_j = (j - 1)*blocks%ButterflyV%inc + blocks%ButterflyV%idx
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1
                        nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                        rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                        allocate (BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix(nn, num_vectors))
                        BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix = 0
                        cnt=cnt+1

                        m_array(cnt)=nn
                        n_array(cnt)=num_vectors
                        k_array(cnt)=rank
                        lda_array(cnt)=nn
                        ldb_array(cnt)=rank
                        ldc_array(cnt)=nn
                        ! a_array(cnt)=LOC(blocks%ButterflyV%blocks(j)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix(1,1))
                        ! c_array(cnt)=LOC(random2(arr_acc_n(j)+1, 1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do j = 1, blocks%ButterflyV%nblk_loc
                        index_j = (j - 1)*blocks%ButterflyV%inc + blocks%ButterflyV%idx
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1

                        cnt=cnt+1

                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),random2(arr_acc_n(j)+1, 1),ldo,CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),k_array(cnt),blocks%ButterflyV%blocks(j)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyV%blocks(j)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly + 1)%blocks(1, index_j_loc_s)%matrix(1,1))
                        ! c_array(cnt)=LOC(random2(arr_acc_n(j)+1, 1))
                        asize=asize+lda_array(cnt)*k_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7

                     call magmaf_wtime(n7)
                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)


                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7
                     cnt=0
                     csize=0
                     do j = 1, blocks%ButterflyV%nblk_loc
                        index_j = (j - 1)*blocks%ButterflyV%inc + blocks%ButterflyV%idx
                        index_j_loc_s = (index_j - BFvec%vec(level_butterfly + 1)%idx_c)/BFvec%vec(level_butterfly + 1)%inc_c + 1
                        cnt=cnt+1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),random2(arr_acc_n(j)+1, 1),ldo)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )
                     n4 = MPI_Wtime()
                     ! time_tmp = time_tmp + n4-n3
                  endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops
               else
                  n3=MPI_Wtime()
                  if (nr0 > 1 .and. inc_r0 == 1) then ! this special treatment makes sure two threads do not write to the same address simultaneously
                     group_count=nr0*nc0
                     do ii = 1, 2

                        call magmaf_wtime(n7)

                        allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                        transa_magma=MagmaTrans
                        transb_magma=MagmaNoTrans
                        alpha_magma=BPACK_cone
                        beta_magma=BPACK_cone

                        cnt=0
                        asize=0
                        bsize=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_i_loc0 = (index_ij - 1)/nc0 + 1
                           index_i_loc = (index_i_loc0 - 1)*2 + ii
                           index_j_loc = mod(index_ij - 1, nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                           index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                           endif

                           cnt = cnt +1
                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=mm
                           lda_array(cnt)=mm
                           ldb_array(cnt)=mm
                           ldc_array(cnt)=nn
                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                           asize=asize+lda_array(cnt)*m_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)


                           mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                           nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                           if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                              allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                              BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                           endif

                           cnt = cnt +1
                           m_array(cnt)=nn
                           n_array(cnt)=num_vectors
                           k_array(cnt)=mm
                           lda_array(cnt)=mm
                           ldb_array(cnt)=mm
                           ldc_array(cnt)=nn
                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))

                           asize=asize+lda_array(cnt)*m_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                        enddo

                        allocate(AA(asize))
                        allocate(BB(bsize))
                        allocate(CC(csize))
                        info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                        info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                        info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                        info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                        info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                        info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                        info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                        info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                        cnt=0
                        asize=0
                        bsize=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_i_loc0 = (index_ij - 1)/nc0 + 1
                           index_i_loc = (index_i_loc0 - 1)*2 + ii
                           index_j_loc = mod(index_ij - 1, nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                           index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           cnt = cnt +1

                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                           call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                           call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))

                           a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                           b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                           c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )
                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                           asize=asize+lda_array(cnt)*m_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                           cnt = cnt +1
                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                           call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                           call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))
                           a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                           b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                           c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )
                           ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                           ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                           ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))

                           asize=asize+lda_array(cnt)*m_array(cnt)
                           bsize=bsize+ldb_array(cnt)*n_array(cnt)
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                        enddo
                        call magmaf_wtime(n8)
                        time_tmp5 = time_tmp5 + n8-n7

                        call magmaf_wtime(n7)
                        ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                        ! stats%Flop_Tmp = stats%Flop_Tmp + flops

                        call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                        call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                        call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                        call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                        call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                        call MAGMA_gemm_vbatched( &
                        transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                        dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                        call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                        call magmaf_wtime(n8)
                        time_tmp5 = time_tmp5 + n8-n7
                        cnt=0
                        csize=0
                        do index_ij = 1, nr0*nc0/2
                           index_i_loc0 = (index_ij - 1)/nc0 + 1
                           index_i_loc = (index_i_loc0 - 1)*2 + ii
                           index_j_loc = mod(index_ij - 1, nc0) + 1  !index_i_loc is local index of row-wise ordering at current level
                           index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                           index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                           index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                           index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                           index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                           index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                           cnt = cnt +1
                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                           cnt = cnt +1
                           call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1),ldc_array(cnt))
                           csize=csize+ldc_array(cnt)*n_array(cnt)

                        enddo
                        deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                        info = MAGMA_free( dA )
                        info = MAGMA_free( dB )
                        info = MAGMA_free( dC )
                        info = MAGMA_free( dA_array )
                        info = MAGMA_free( dB_array )
                        info = MAGMA_free( dC_array )
                        info = MAGMA_free( dm_array )
                        info = MAGMA_free( dn_array )
                        info = MAGMA_free( dk_array )
                        info = MAGMA_free( dlda_array )
                        info = MAGMA_free( dldb_array )
                        info = MAGMA_free( dldc_array )
                     enddo
                  else

                     call magmaf_wtime(n7)

                     group_count=nr0*nc0*2
                     allocate(m_array(group_count),n_array(group_count),k_array(group_count),lda_array(group_count), ldb_array(group_count), ldc_array(group_count),a_array(group_count),b_array(group_count), c_array(group_count))
                     transa_magma=MagmaTrans
                     transb_magma=MagmaNoTrans
                     alpha_magma=BPACK_cone
                     beta_magma=BPACK_cone

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1
                        index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                        index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
                           allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                           BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0
                        endif
                        cnt = cnt + 1
                        m_array(cnt)=nn
                        n_array(cnt)=num_vectors
                        k_array(cnt)=mm
                        lda_array(cnt)=mm
                        ldb_array(cnt)=mm
                        ldc_array(cnt)=nn
                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 1)
                        nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                        if (.not. associated(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix)) then
                           allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix(nn, num_vectors))
                           BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s + 1)%matrix = 0
                        endif
                        cnt = cnt + 1
                        m_array(cnt)=nn
                        n_array(cnt)=num_vectors
                        k_array(cnt)=mm
                        lda_array(cnt)=mm
                        ldb_array(cnt)=mm
                        ldc_array(cnt)=nn
                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))
                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     allocate(AA(asize))
                     allocate(BB(bsize))
                     allocate(CC(csize))
                     info = MAGMA_malloc( dA,          asize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dB,          bsize * C_SIZEOF_DT )
                     info = MAGMA_malloc( dC,          csize * C_SIZEOF_DT )

                     info = MAGMA_malloc(    dA_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dB_array, group_count * sizeof_ptr )
                     info = MAGMA_malloc(    dC_array, group_count * sizeof_ptr )

                     info = MAGMA_malloc( dm_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dn_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dk_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dlda_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldb_array, (group_count+1) * sizeof_int )
                     info = MAGMA_malloc( dldc_array, (group_count+1) * sizeof_int )

                     cnt=0
                     asize=0
                     bsize=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1
                        index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                        index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        cnt = cnt + 1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))
                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1))
                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)


                        cnt = cnt + 1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1),ldc_array(cnt),CC(csize+1),ldc_array(cnt))
                        call copymatf77(' ',ldb_array(cnt),n_array(cnt),BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1),ldb_array(cnt),BB(bsize+1),ldb_array(cnt))
                        call copymatf77(' ',lda_array(cnt),m_array(cnt),blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1),lda_array(cnt),AA(asize+1),lda_array(cnt))
                        a_array(cnt)    = MAGMA_offset_1d( dA, 1, asize+1 )
                        b_array(cnt)    = MAGMA_offset_1d( dB, 1, bsize+1 )
                        c_array(cnt)    = MAGMA_offset_1d( dC, 1, csize+1 )

                        ! a_array(cnt)=LOC(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(1,1))
                        ! b_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix(1,1))
                        ! c_array(cnt)=LOC(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1))
                        asize=asize+lda_array(cnt)*m_array(cnt)
                        bsize=bsize+ldb_array(cnt)*n_array(cnt)
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7
                     call magmaf_wtime(n7)
                     ! call gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size,flop=flops)
                     ! stats%Flop_Tmp = stats%Flop_Tmp + flops
                     call MAGMA_setvector(asize, int(C_SIZEOF_DT), c_loc(AA), 1, dA, 1, queue )
                     call MAGMA_setvector(bsize, int(C_SIZEOF_DT), c_loc(BB), 1, dB, 1, queue )
                     call MAGMA_setvector(csize, int(C_SIZEOF_DT), c_loc(CC), 1, dC, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int), c_loc(m_array), 1, dm_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(n_array), 1, dn_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(k_array), 1, dk_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(lda_array), 1, dlda_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldb_array), 1, dldb_array, 1, queue )
                     call MAGMA_setvector(group_count, int(sizeof_int),c_loc(ldc_array), 1, dldc_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(a_array), 1,    dA_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(b_array), 1,    dB_array, 1, queue )
                     call MAGMA_setvector( group_count, int(sizeof_ptr),   c_loc(c_array), 1,    dC_array, 1, queue )


                     call MAGMA_gemm_vbatched( &
                     transa_magma, transb_magma, dm_array, dn_array, dk_array, alpha_magma, dA_array, &
                     dlda_array, dB_array, dldb_array,beta_magma,dC_array, dldc_array, group_count, queue)

                     call MAGMA_getvector(csize, int(C_SIZEOF_DT), dC, 1, c_loc(CC), 1, queue )
                     call magmaf_wtime(n8)
                     time_tmp5 = time_tmp5 + n8-n7
                     cnt=0
                     csize=0
                     do index_ij = 1, nr0*nc0
                        index_j_loc = (index_ij - 1)/nr0 + 1
                        index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                        index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                        index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                        index_i_loc_s = (index_ii - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_s = (index_jj - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        cnt = cnt + 1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1,1),ldc_array(cnt))
                        csize=csize+ldc_array(cnt)*n_array(cnt)

                        cnt = cnt + 1
                        call copymatf77(' ',ldc_array(cnt),n_array(cnt),CC(csize+1),ldc_array(cnt),BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s+1)%matrix(1,1),ldc_array(cnt))
                        csize=csize+ldc_array(cnt)*n_array(cnt)
                     enddo

                     deallocate(m_array,n_array,k_array,lda_array,ldb_array,ldc_array,a_array,b_array,c_array,AA,BB,CC)
                     info = MAGMA_free( dA )
                     info = MAGMA_free( dB )
                     info = MAGMA_free( dC )
                     info = MAGMA_free( dA_array )
                     info = MAGMA_free( dB_array )
                     info = MAGMA_free( dC_array )
                     info = MAGMA_free( dm_array )
                     info = MAGMA_free( dn_array )
                     info = MAGMA_free( dk_array )
                     info = MAGMA_free( dlda_array )
                     info = MAGMA_free( dldb_array )
                     info = MAGMA_free( dldc_array )
                  endif
                  n4 = MPI_Wtime()
                  ! time_tmp = time_tmp + n4-n3
               endif
            endif

            do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
               do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                  if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
               enddo
            enddo

            if (level /= 0) then
               call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'R')
            endif
         enddo
         n6 = MPI_Wtime()
         time_tmp1 = time_tmp1 + n6-n5

         if (isnanMat(random2(1:N,1:1),N,1)) then
            write (*, *) 'NAN in 4 BF_block_MVP_dat', blocks%row_group, blocks%col_group, blocks%level, blocks%level_butterfly
            stop
         end if

      endif

      do level = 0, level_butterfly + 2
         do j = 1, BFvec%vec(level)%nc
            do i = 1, BFvec%vec(level)%nr
               if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
            enddo
         enddo
         if (allocated(BFvec%vec(level)%blocks)) deallocate (BFvec%vec(level)%blocks)
      enddo

      deallocate (BFvec%vec)
      deallocate (arr_acc_m, arr_acc_n)


   endif

   n2 = MPI_Wtime()
   ! time_tmp = time_tmp + n2-n1
   time_tmp2 = time_tmp2 + n2-n1


   call MAGMA_queue_destroy( queue )
   call MAGMA_finalize()

   return

end subroutine BF_block_MVP_dat_batch_magma
#endif



!>**** Matvec of partial levels of BF with vectors
! if chara=='N', out=BF(level_end:0)*vec, if chara=='T', out=vec*BF(level_butterfly+1:level_end).
   !blocks: working BF
   !chara: 'N' or 'T'
   !num_vectors: number of vectors
   !VectIn: dimension (mnloc,num_vectors) the local input vectors
   !BFvec: storing the result of partial matvec
   !level_end: the last level of the operator
   !ptree: process tree
   !stats: statistics
   subroutine BF_block_MVP_partial(blocks, chara, num_vectors, VectIn, BFvec, level_end, ptree, stats)



      implicit none

      integer index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, level_half, level_final
      integer idx_r, inc_r, nr, idx_c, inc_c, nc
      integer idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, pgno_sub
      DT ctemp
      character chara
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno, comm, ierr
      type(Hstat)::stats
      real(kind=8)::flop, flops
      integer index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc, index_i_loc_s, index_i_loc_k, index_j_loc, index_j_loc_s, index_j_loc_k, level_end

      type(butterfly_vec) :: BFvec
      DT :: VectIn(:, :)
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), Vout_tmp(:, :)

      integer, allocatable:: arr_acc_m(:), arr_acc_n(:)

      integer idxs, groupn_start, groupm_start

      level_butterfly = blocks%level_butterfly
      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm
      if (comm == MPI_COMM_NULL) then
         write (*, *) 'ninin', pgno, comm == MPI_COMM_NULL, ptree%MyID
      endif

      call assert(IOwnPgrp(ptree, pgno), 'I do not share this block!')

      if (BF_checkNAN(blocks)) then
         write (*, *) 'NAN in 0 BF_block_MVP_partial'
         stop
      end if

      if (chara == 'N') then
         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         level_half = blocks%level_half
         call assert(level_half >= level_end, 'partial matvec with chara=N requires row-wise ordering')

         if (.not. allocated(BFvec%vec)) allocate (BFvec%vec(0:level_butterfly + 2))

         allocate (BFvec%vec(0)%blocks(1, blocks%ButterflyV%nblk_loc))
         BFvec%vec(0)%num_row = 1
         BFvec%vec(0)%num_col = num_blocks
         BFvec%vec(0)%idx_r = 1
         BFvec%vec(0)%inc_r = 1
         BFvec%vec(0)%nr = 1
         BFvec%vec(0)%idx_c = blocks%ButterflyV%idx
         BFvec%vec(0)%inc_c = blocks%ButterflyV%inc
         BFvec%vec(0)%nc = blocks%ButterflyV%nblk_loc

         groupn_start = blocks%col_group*2**level_butterfly
         call GetLocalBlockRange(ptree, blocks%pgno, 0, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
         if (BFvec%vec(0)%nc > 1) then
            call assert(associated(blocks%ns),'blocks%ns not computed in BF_block_MVP_partial')
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(i,nn,ii,jj,idxs)
#endif
            do i = 1, BFvec%vec(0)%nc
               if(i==1)then
                  idxs = 0
                  nn = blocks%ns(i)
               else
                  idxs = blocks%ns(i-1)
                  nn = blocks%ns(i) - blocks%ns(i-1)
               endif
               allocate (BFvec%vec(0)%blocks(1, i)%matrix(nn, num_vectors))
               do ii = 1, nn
                  do jj = 1, num_vectors
                     BFvec%vec(0)%blocks(1, i)%matrix(ii, jj) = VectIn(ii + idxs, jj)
                  enddo
               enddo
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
         else
            allocate (BFvec%vec(0)%blocks(1, 1)%matrix(blocks%N_loc, num_vectors))
            BFvec%vec(0)%blocks(1, 1)%matrix = VectIn
         endif

         do level = 0, level_end
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

            BFvec%vec(level + 1)%idx_r = idx_r
            BFvec%vec(level + 1)%inc_r = inc_r
            BFvec%vec(level + 1)%nr = nr
            BFvec%vec(level + 1)%idx_c = idx_c
            BFvec%vec(level + 1)%inc_c = inc_c
            BFvec%vec(level + 1)%nc = nc

            if (nr > 0 .and. nc > 0) then

               if (level /= level_butterfly + 1) then
                  BFvec%vec(level + 1)%num_row = 2**level
                  BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
               else
                  BFvec%vec(level + 1)%num_row = 2**level_butterfly
                  BFvec%vec(level + 1)%num_col = 1
               endif
               if (level_half /= level) then ! the last level doesn't require doubling block columns
               if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                  BFvec%vec(level + 1)%nc = 2
                  BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
               endif
               endif
               allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))

               if (level == 0) then
                  flops = 0
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(j,rank,nn,flop,index_j,index_j_loc_s,idxs) reduction(+:flops)
#endif
                  do j = 1, blocks%ButterflyV%nblk_loc
                     index_j = (j - 1)*inc_c + idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                     nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                     allocate (BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix(rank, num_vectors))
                     BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix = 0
                     call gemmf90(blocks%ButterflyV%blocks(j)%matrix, nn, BFvec%vec(0)%blocks(1, j)%matrix, nn, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank, 'T', 'N', rank, num_vectors, nn, BPACK_cone, BPACK_czero, flop=flop)
                     flops = flops + flop
                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops

                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                     index_j = idx_c
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(1, index_j_loc_s)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                  endif

               elseif (level == level_butterfly + 1) then
                  flops = 0
                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     mm = size(blocks%ButterflyU%blocks(1)%matrix, 1)
                     rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                     allocate (matrixtemp(rank, num_vectors))
                     matrixtemp = 0
                     if (ptree%MyID == ptree%pgrp(pgno_sub)%head) matrixtemp = BFvec%vec(level)%blocks(1, 1)%matrix
                     call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                     allocate (BFvec%vec(level + 1)%blocks(1, 1)%matrix(mm, num_vectors))
                     BFvec%vec(level + 1)%blocks(1, 1)%matrix = 0
                     call gemmf90(blocks%ButterflyU%blocks(1)%matrix, mm, matrixtemp, rank, BFvec%vec(level + 1)%blocks(1, 1)%matrix, mm, 'N', 'N', mm, num_vectors, rank, BPACK_cone, BPACK_czero, flop=flop)
                     flops = flops + flop
                     deallocate (matrixtemp)
                  else
#ifdef HAVE_OPENMP
                     !$omp parallel do default(shared) private(i,rank,mm,flop) reduction(+:flops)
#endif
                     do i = 1, blocks%ButterflyU%nblk_loc
                        rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                        mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                        allocate (BFvec%vec(level + 1)%blocks(i, 1)%matrix(mm, num_vectors))
                        BFvec%vec(level + 1)%blocks(i, 1)%matrix = 0

                        call gemmf90(blocks%ButterflyU%blocks(i)%matrix, mm, BFvec%vec(level)%blocks(i, 1)%matrix, rank, BFvec%vec(level + 1)%blocks(i, 1)%matrix, mm, 'N', 'N', mm, num_vectors, rank, BPACK_cone, BPACK_czero, flop=flop)
                        flops = flops + flop
                     enddo
#ifdef HAVE_OPENMP
                     !$omp end parallel do
#endif
                  endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops

               else
                  flops = 0
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
#endif
                  do index_ij = 1, nr*nc
                     index_j_loc = (index_ij - 1)/nr + 1
                     index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                     index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                     index_j = (index_j_loc - 1)*inc_c + idx_c

                     index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

                     index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                     index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                     index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                     index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                     index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1
                     index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                     nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                     nn2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                     mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                     allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, num_vectors))
                     BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0

                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn1, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn1, BPACK_cone, BPACK_cone, flop=flop)
                     flops = flops + flop

                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix, nn2, BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, mm, 'N', 'N', mm, num_vectors, nn2, BPACK_cone, BPACK_cone, flop=flop)
                     flops = flops + flop
                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops
               endif
            endif

            do j = 1, BFvec%vec(level)%nc
               do i = 1, BFvec%vec(level)%nr
                  if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
               enddo
            enddo
            if (level_half /= level) then
               call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'B')
            endif
         enddo
      elseif (chara == 'T') then
         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         level_half = blocks%level_half
         call assert(level_half + 1 <= level_end, 'partial matvec with chara=T requires column-wise ordering')

         if (.not. allocated(BFvec%vec)) allocate (BFvec%vec(0:level_butterfly + 2))
         allocate (BFvec%vec(0)%blocks(blocks%ButterflyU%nblk_loc, 1))
         BFvec%vec(0)%num_row = num_blocks
         BFvec%vec(0)%num_col = 1
         BFvec%vec(0)%idx_r = blocks%ButterflyU%idx
         BFvec%vec(0)%inc_r = blocks%ButterflyU%inc
         BFvec%vec(0)%nr = blocks%ButterflyU%nblk_loc
         BFvec%vec(0)%idx_c = 1
         BFvec%vec(0)%inc_c = 1
         BFvec%vec(0)%nc = 1

         groupm_start = blocks%row_group*2**level_butterfly
         call GetLocalBlockRange(ptree, blocks%pgno, level_butterfly + 1, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')
         if (BFvec%vec(0)%nr > 1) then
            call assert(associated(blocks%ms),'blocks%ms not computed in BF_block_MVP_partial')
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(i,mm,ii,jj,idxs)
#endif
            do i = 1, BFvec%vec(0)%nr
               if(i==1)then
                  idxs = 0
                  mm = blocks%ms(i)
               else
                  idxs = blocks%ms(i-1)
                  mm = blocks%ms(i) - blocks%ms(i-1)
               endif
               allocate (BFvec%vec(0)%blocks(i, 1)%matrix(mm, num_vectors))
               do ii = 1, mm
                  do jj = 1, num_vectors
                     BFvec%vec(0)%blocks(i, 1)%matrix(ii, jj) = VectIn(ii + idxs, jj)
                  enddo
               enddo
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
         else
            allocate (BFvec%vec(0)%blocks(1, 1)%matrix(blocks%M_loc, num_vectors))
            BFvec%vec(0)%blocks(1, 1)%matrix = VectIn
         endif

         do level = level_butterfly + 1, level_end, -1
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

            BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
            BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
            BFvec%vec(level_butterfly - level + 2)%nr = nr
            BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
            BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
            BFvec%vec(level_butterfly - level + 2)%nc = nc

            if (nr > 0 .and. nc > 0) then

               if (level /= 0) then
                  BFvec%vec(level_butterfly - level + 2)%num_row = 2**(level - 1)
                  BFvec%vec(level_butterfly - level + 2)%num_col = 2**(level_butterfly - level + 1)
               else
                  BFvec%vec(level_butterfly - level + 2)%num_row = 1
                  BFvec%vec(level_butterfly - level + 2)%num_col = 2**level_butterfly
               endif
               if (level_half + 1 /= level) then ! the last level doesn't require doubling block rows
               if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                  BFvec%vec(level_butterfly - level + 2)%nr = 2
                  BFvec%vec(level_butterfly - level + 2)%idx_r = BFvec%vec(level_butterfly - level + 2)%idx_r - 1 + mod(BFvec%vec(level_butterfly - level + 2)%idx_r, 2)
               endif
               endif
               allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))

               if (level == level_butterfly + 1) then
                  flops = 0
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(i,rank,mm,flop,index_i,index_i_loc_s,idxs) reduction(+:flops)
#endif
                  do i = 1, blocks%ButterflyU%nblk_loc
                     index_i = (i - 1)*blocks%ButterflyU%inc + blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1

                     rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
                     mm = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                     allocate (BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix(rank, num_vectors))
                     BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix = 0
                     call gemmf90(blocks%ButterflyU%blocks(i)%matrix, mm, BFvec%vec(0)%blocks(i, 1)%matrix, mm, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank, 'T', 'N', rank, num_vectors, mm, BPACK_cone, BPACK_czero, flop=flop)
                     flops = flops + flop

                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops

                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, idx_r, 1, 'C', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
                     index_i = blocks%ButterflyU%idx
                     index_i_loc_s = (index_i - BFvec%vec(1)%idx_r)/BFvec%vec(1)%inc_r + 1
                     call MPI_ALLREDUCE(MPI_IN_PLACE, BFvec%vec(1)%blocks(index_i_loc_s, 1)%matrix, rank*num_vectors, MPI_DT, MPI_SUM, ptree%pgrp(pgno_sub)%Comm, ierr)
                  endif

               elseif (level == 0) then
                  flops = 0

                  call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, 1, idx_c, 'R', pgno_sub)
                  if (ptree%pgrp(pgno_sub)%nproc > 1) then
                     nn = size(blocks%ButterflyV%blocks(1)%matrix, 1)
                     rank = size(blocks%ButterflyV%blocks(1)%matrix, 2)
                     allocate (matrixtemp(rank, num_vectors))
                     matrixtemp = 0
                     if (ptree%MyID == ptree%pgrp(pgno_sub)%head) matrixtemp = BFvec%vec(level_butterfly + 1)%blocks(1, 1)%matrix
                     call MPI_Bcast(matrixtemp, rank*num_vectors, MPI_DT, Main_ID, ptree%pgrp(pgno_sub)%Comm, ierr)
                     allocate (BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix(nn, num_vectors))
                     BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix = 0
                     call gemmf90(blocks%ButterflyV%blocks(1)%matrix, nn, matrixtemp, rank, BFvec%vec(level_butterfly + 2)%blocks(1, 1)%matrix, nn, 'N', 'N', nn, num_vectors, rank, BPACK_cone, BPACK_czero, flop=flop)
                     flops = flops + flop
                     deallocate (matrixtemp)
                  else
#ifdef HAVE_OPENMP
                     !$omp parallel do default(shared) private(j,rank,nn,flop) reduction(+:flops)
#endif
                     do j = 1, blocks%ButterflyV%nblk_loc
                        nn = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                        rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                        allocate (BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix(nn, num_vectors))
                        BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix = 0
                        call gemmf90(blocks%ButterflyV%blocks(j)%matrix, nn, BFvec%vec(level_butterfly + 1)%blocks(1, j)%matrix, rank, BFvec%vec(level_butterfly + 2)%blocks(1, j)%matrix, nn, 'N', 'N', nn, num_vectors, rank, BPACK_cone, BPACK_czero, flop=flop)
                        flops = flops + flop
                     enddo
#ifdef HAVE_OPENMP
                     !$omp end parallel do
#endif
                  endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops
               else
                  flops = 0
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop) reduction(+:flops)
#endif
                  do index_ij = 1, nr*nc
                     index_j_loc = (index_ij - 1)/nr + 1
                     index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                     index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                     index_j = (index_j_loc - 1)*inc_c + idx_c

                     index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                     index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                     index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                     index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                     index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                     index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1
                     index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                     mm1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                     mm2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 1)
                     nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                     allocate (BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(nn, num_vectors))
                     BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix = 0

                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm1, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm1, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm1, BPACK_cone, BPACK_cone, flop=flop)
                     flops = flops + flop
                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, mm2, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc + 1, index_jj_loc)%matrix, mm2, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, nn, 'T', 'N', nn, num_vectors, mm2, BPACK_cone, BPACK_cone, flop=flop)
                     flops = flops + flop

                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
                  stats%Flop_Tmp = stats%Flop_Tmp + flops
               endif
            endif

            do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
               do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                  if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
               enddo
            enddo

            if (level_half + 1 /= level) then
               call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'B')
            endif
         enddo
      endif

      return

   end subroutine BF_block_MVP_partial

   subroutine Full_block_extraction(blocks, inters, ptree, msh, stats,option)

      implicit none
      type(matrixblock)::blocks
      type(Hoption)::option
      type(proctree)::ptree
      type(mesh)::msh
      type(Hstat)::stats
      integer nn, ri, ci, nng, headm, headn, pp, row_group, col_group
      type(intersect)::inters(:)
      integer ii, jj
      real(kind=8)::tol_used
#if HAVE_ZFP
      if(option%use_zfp==1)call ZFP_Decompress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,1)
#endif
      headm = msh%basis_group(blocks%row_group)%head
      headn = msh%basis_group(blocks%col_group)%head
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         do ii = 1, blocks%inters(nn)%nr_loc
            ri = inters(nng)%rows(blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii))) - headm + 1
            do jj = 1, blocks%inters(nn)%nc
               ci = inters(nng)%cols(blocks%inters(nn)%cols(jj)) - headn + 1
               blocks%inters(nn)%dat_loc(ii, jj) = blocks%fullmat(ri, ci)
            enddo
         enddo
      enddo
#if HAVE_ZFP
      if(option%use_zfp==1)call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,option%tol_comp,1)
#endif
   end subroutine Full_block_extraction

   subroutine LR_block_extraction(blocks, inters, ptree, msh, stats)

      implicit none
      type(matrixblock)::blocks
      type(proctree)::ptree
      type(mesh)::msh
      type(Hstat)::stats
      integer nn, ri, ci, nng, headm, headn, head, tail, pp, row_group, col_group
      type(intersect)::inters(:)
      integer ii, jj, rank, ncol, nrow, iidx, pgno, ierr, nr_loc
      DT, allocatable::Vpartial(:, :), matU(:, :)
      real(kind=8)::t1, t2, t3, t4

      headm = blocks%headm
      headn = blocks%headn
      pgno = blocks%pgno

      pp = ptree%myid - ptree%pgrp(pgno)%head + 1
      head = blocks%N_p(pp, 1)
      tail = blocks%N_p(pp, 2)

      t1 = MPI_Wtime()
      rank = size(blocks%ButterflyU%blocks(1)%matrix, 2)
      ncol = 0
      do nn = 1, size(blocks%inters, 1)
         ncol = ncol + blocks%inters(nn)%nc
      enddo
      allocate (Vpartial(rank, ncol))
      t2 = MPI_Wtime()
      call LR_all2all_extraction(blocks, inters, Vpartial, rank, ncol, stats, ptree, msh)
      t3 = MPI_Wtime()

      nr_loc = 0
      do nn = 1, size(blocks%inters, 1)
         nr_loc = max(nr_loc, blocks%inters(nn)%nr_loc)
      enddo
      allocate (matU(nr_loc, rank))

      iidx = 0
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         if (blocks%inters(nn)%nr_loc > 0) then
            do ii = 1, blocks%inters(nn)%nr_loc
               ri = inters(nng)%rows(blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii))) - headm + 1 - blocks%M_p(pp, 1) + 1
               matU(ii, :) = blocks%ButterflyU%blocks(1)%matrix(ri, :)
            enddo
            if (blocks%inters(nn)%nc > 0) then
               call gemmf77('N', 'N', blocks%inters(nn)%nr_loc, blocks%inters(nn)%nc, rank, BPACK_cone, matU, nr_loc, Vpartial(1, iidx + 1), rank, BPACK_czero, blocks%inters(nn)%dat_loc(1, 1), blocks%inters(nn)%nr_loc)
            endif
            stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(blocks%inters(nn)%nr_loc, blocks%inters(nn)%nc, rank)
         endif
         iidx = iidx + blocks%inters(nn)%nc
      enddo

      deallocate (Vpartial)
      deallocate (matU)

      t4 = MPI_Wtime()
      ! time_tmp = time_tmp + t3 - t2

   end subroutine LR_block_extraction

!>*********** all to all communication of columns in the V factor from the 1D block column layout to that needed by the 1D block row layout
   subroutine LR_all2all_extraction(blocks, inters, Vpartial, rank, ncol, stats, ptree, msh)


      implicit none
      integer i, j, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt
      integer level, length_1, length_2, level_blocks
      integer level_p, ncol, rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      real(kind=8) flop
      DT ctemp
      type(matrixblock)::blocks
      DT::Vpartial(rank, ncol)
      type(Hstat)::stats
      type(proctree)::ptree
      type(mesh)::msh
      type(intersect)::inters(:)
      integer ri, ci, nng, head, tail, row_group, col_group
      integer, allocatable:: select_row(:), select_column(:), column_pivot(:), row_pivot(:)
      integer, allocatable:: select_row_rr(:), select_column_rr(:)
      DT, allocatable:: UU(:, :), VV(:, :), matrix_little(:, :), matrix_little_inv(:, :), matrix_U(:, :), matrix_V(:, :), matrix_V_tmp(:, :), matrix_little_cc(:, :), core(:, :), tau(:)

      integer, allocatable::jpvt(:)
      integer ierr, nsendrecv, pid, tag, nproc, Nrow, Nreqr, Nreqs, recvid, sendid, tmpi
      integer idx_r, idx_c, inc_r, inc_c, nr, nc, level_new

      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      character::mode, mode_new
      real(kind=8)::n1, n2, n3, n4
      integer, allocatable::sendIDactive(:), recvIDactive(:)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      type(butterfly_kerl)::kerls, kerls1
      integer rr, cc
      logical all2all
      integer, allocatable::sdispls(:), sendcounts(:), rdispls(:), recvcounts(:), col_idx_loc(:), activeproc(:, :), row_idx(:, :)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist, pgno, ncol_loc, iidx
      integer::headm, headn, nc_loc, ncmax, nrmax, head1, tail1

      n1 = MPI_Wtime()

      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno
      level_p = ptree%nlevel - GetTreelevel(blocks%pgno)
      headm = blocks%headm
      headn = blocks%headn
      pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
      head = blocks%N_p(pp, 1)
      tail = blocks%N_p(pp, 2)

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

      ncmax = 0
      nrmax = 0
      do nn = 1, size(blocks%inters, 1)
         ncmax = max(blocks%inters(nn)%nc, ncmax)
         nrmax = max(blocks%inters(nn)%nr, nrmax)
      enddo
      allocate (col_idx_loc(ncmax))
      allocate (row_idx(nrmax, 1))
      allocate (activeproc(nproc, size(blocks%inters, 1)))
      activeproc = 0

      ! iidx=0
      do nn = 1, size(blocks%inters, 1)

         nng = blocks%inters(nn)%idx
         nc_loc = 0
         do jj = 1, blocks%inters(nn)%nc
            ci = inters(nng)%cols(blocks%inters(nn)%cols(jj)) - headn + 1
            if (ci >= head .and. ci <= tail) then
               nc_loc = nc_loc + 1
               col_idx_loc(nc_loc) = jj
            endif
         enddo

         row_idx(1:blocks%inters(nn)%nr, 1) = 0
         do ii = 1, blocks%inters(nn)%nr
            ri = inters(nng)%rows(blocks%inters(nn)%rows(ii)) - headm + 1
            row_idx(ii, 1) = ri
         enddo
         call PIKSRT_INT_Multi(blocks%inters(nn)%nr, 1, row_idx)

         ii = 1
         do pp = 1, nproc
            head1 = blocks%M_p(pp, 1)
            tail1 = blocks%M_p(pp, 2)
            do while (row_idx(ii, 1) < head1)
               ii = ii + 1
               if (ii > blocks%inters(nn)%nr) exit
            enddo
            if (ii <= blocks%inters(nn)%nr) then
               if (tail1 >= row_idx(ii, 1)) activeproc(pp, nn) = 1
            else
               exit
            endif
         enddo

         do pp = 1, nproc
            if (activeproc(pp, nn) == 1) then
            do jj = 1, nc_loc
               ci = inters(nng)%cols(blocks%inters(nn)%cols(col_idx_loc(jj))) - headn + 1
               if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
               endif
               sendquant(pp)%size = sendquant(pp)%size + (1 + rank)
            enddo
            endif
         enddo

      enddo
      deallocate (row_idx)

      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         if (blocks%inters(nn)%nr_loc > 0) then
            do jj = 1, blocks%inters(nn)%nc
               ci = inters(nng)%cols(blocks%inters(nn)%cols(jj))
               pgno = findpggroup(ci, msh, ptree, blocks%col_group, blocks%pgno)
               pp = ptree%pgrp(pgno)%head - ptree%pgrp(blocks%pgno)%head + 1
               if (recvquant(pp)%active == 0) then
                  recvquant(pp)%active = 1
                  Nrecvactive = Nrecvactive + 1
                  recvIDactive(Nrecvactive) = pp
               endif
               recvquant(pp)%size = recvquant(pp)%size + (1 + rank)
            enddo
         endif
      enddo

      n2 = MPI_Wtime()

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      iidx = 0
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         nc_loc = 0
         do jj = 1, blocks%inters(nn)%nc
            ci = inters(nng)%cols(blocks%inters(nn)%cols(jj)) - headn + 1
            if (ci >= head .and. ci <= tail) then
               nc_loc = nc_loc + 1
               col_idx_loc(nc_loc) = jj
            endif
         enddo

         do pp = 1, nproc
            if (activeproc(pp, nn) == 1) then
            do jj = 1, nc_loc
               ci = inters(nng)%cols(blocks%inters(nn)%cols(col_idx_loc(jj))) - headn + 1
               sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = iidx + col_idx_loc(jj)
               sendquant(pp)%size = sendquant(pp)%size + 1
               sendquant(pp)%dat(sendquant(pp)%size + 1:sendquant(pp)%size + rank, 1) = blocks%ButterflyV%blocks(1)%matrix(ci - head + 1, :)
               sendquant(pp)%size = sendquant(pp)%size + rank
            enddo
            endif
         enddo
         iidx = iidx + blocks%inters(nn)%nc
      enddo
      deallocate (col_idx_loc)
      deallocate (activeproc)

      n3 = MPI_Wtime()

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            iidx = NINT(dble(recvquant(pp)%dat(i, 1)))
            Vpartial(:, iidx) = recvquant(pp)%dat(i + 1:i + rank, 1)
            i = i + rank
         enddo
      enddo

      n4 = MPI_Wtime()
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)

      ! time_tmp = time_tmp + n4 - n1

   end subroutine LR_all2all_extraction



   subroutine BF_MD_block_extraction(blocks, Ndim, Ninter, inters, ptree, msh, stats, option)

      implicit none
      integer Ndim,dim_i,Ninter
      type(matrixblock_MD)::blocks
      type(proctree)::ptree
      type(mesh)::msh(Ndim)
      type(Hstat)::stats
      type(Hoption)::option
      integer ri, ci, nng, headm, headn, head, tail, pp, row_group, col_group
      type(intersect_MD)::inters(Ninter)
      integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
      integer level_butterfly, num_blocks, level_half, levelm
      integer mm, nn, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1(Ndim), nvec2(Ndim)
      integer level
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      integer nsendrecv, pid, pgno_sub, pid0, tag, nproc, Nreqr, Nreqs, recvid, sendid
      integer idx_r_m(Ndim),idx_c_m(Ndim),idx_r_m_loc(Ndim),idx_c_m_loc(Ndim)
      type(butterfly_vec),allocatable :: BFvec(:,:), BFvec1(:,:), BFvec_transposed(:,:)
      type(butterfly_kerl), allocatable::g_idx_m(:), g_idx_n(:)
      integer, allocatable::group_ms(:), group_ns(:), group_ms1(:), group_ns1(:)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(ipair):: p
      DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :), data_buffer(:)
      DT, allocatable::Vpartial(:, :)
      DT::val
      integer idxr, idxc, sender, receiver
      integer, allocatable::num_nods_i(:), num_nods_j(:)
      real(kind=8)::n1, n2, n3, n4, n5
      integer dims_c_m(Ndim),dims_r_m(Ndim),bbm_r,bbm_c,group_n(Ndim),group_m(Ndim),dims_MD_old(Ndim*2),dims_MD_new(Ndim*2),dim_subtensor(Ndim*2),idx_subtensor(Ndim*2), index_ij, index_scalar,index_vector(Ndim),index_scalar_old,index_vector_old(Ndim)
      real(kind=8)::tol_used

      n1 = MPI_Wtime()

      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm

      pp = ptree%myid - ptree%pgrp(pgno)%head + 1


      level_butterfly = blocks%level_butterfly
      level_half = blocks%level_half
      levelm = level_half
      dims_c_m = blocks%nc_m
      dims_r_m = blocks%nr_m

      ! write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,'inin'


      !>**** Step 1: multiply the factor matrices out from outtermost to middle level along each dimension
      !>******* preallocate BFvec and BFvec1, number of blocks are exactly those in BF_block_MVP_dat, BFvec for the first 0:level_half+1 levels and BFvec for the next level_half+1 to level_butterfly+2 levels,note that level level_half+1 is duplicated for all2all communication
      allocate(BFvec(Ninter,Ndim))
      allocate(BFvec1(Ninter,Ndim))
      allocate(BFvec_transposed(Ninter,Ndim))
      ! write(*,*)ptree%MyID,blocks%idx_r_m,blocks%idx_c_m,'wogan'
      do nn=1,Ninter
         idx_r_m = inters(nn)%idx_r_m
         idx_r_m_loc = idx_r_m - blocks%idx_r_m + 1
         idx_c_m = inters(nn)%idx_c_m
         idx_c_m_loc = idx_c_m - blocks%idx_c_m + 1

         receiver = inters(nn)%receiver
         sender = inters(nn)%sender

         if(sender==pp)then
            call MultiIndexToSingleIndex(Ndim, dims_c_m, bbm_c, idx_c_m_loc)
            call assert(0<bbm_c .and. bbm_c<=product(dims_c_m),"bbm_c incorrect")
            ! write(*,*)"inter #",nn,"sender c",idx_c_m,"myID",ptree%MyID
         endif

         if(receiver==pp)then
            call MultiIndexToSingleIndex(Ndim, dims_r_m, bbm_r, idx_r_m_loc)
            call assert(0<bbm_r .and. bbm_r<=product(dims_r_m),"bbm_c incorrect")
            ! write(*,*)"inter #",nn,"receiver r",idx_r_m,"myID",ptree%MyID
         endif

         group_n = blocks%col_group
         group_m = blocks%row_group
         do dim_i=1,Ndim
            group_n(dim_i) = group_n(dim_i)*2**(level_butterfly-levelm) - 1 + idx_c_m(dim_i)
            group_m(dim_i) = group_m(dim_i)*2**levelm - 1 + idx_r_m(dim_i)
         enddo

         if(sender==pp)then
            num_blocks = 2**levelm
            do dim_i=1,Ndim
               allocate (BFvec(nn,dim_i)%vec(1:level_half + 1))
               allocate (BFvec(nn,dim_i)%vec(1)%blocks(1, num_blocks))
               BFvec(nn,dim_i)%vec(1)%num_row = 1
               BFvec(nn,dim_i)%vec(1)%num_col = num_blocks
               BFvec(nn,dim_i)%vec(1)%idx_r = 1
               BFvec(nn,dim_i)%vec(1)%inc_r = 1
               BFvec(nn,dim_i)%vec(1)%nr = 1
               BFvec(nn,dim_i)%vec(1)%idx_c = 1
               BFvec(nn,dim_i)%vec(1)%inc_c = 1
               BFvec(nn,dim_i)%vec(1)%nc = num_blocks

               do level = 1, level_half
                  BFvec(nn,dim_i)%vec(level + 1)%idx_r = 1
                  BFvec(nn,dim_i)%vec(level + 1)%inc_r = 1
                  BFvec(nn,dim_i)%vec(level + 1)%nr = 1
                  BFvec(nn,dim_i)%vec(level + 1)%idx_c = 1
                  BFvec(nn,dim_i)%vec(level + 1)%inc_c = 1
                  BFvec(nn,dim_i)%vec(level + 1)%nc = BFvec(nn,dim_i)%vec(level)%nc/2
                  BFvec(nn,dim_i)%vec(level + 1)%num_row = 1
                  BFvec(nn,dim_i)%vec(level + 1)%num_col = BFvec(nn,dim_i)%vec(level)%nc/2
                  allocate (BFvec(nn,dim_i)%vec(level + 1)%blocks(BFvec(nn,dim_i)%vec(level + 1)%nr, BFvec(nn,dim_i)%vec(level + 1)%nc))
               enddo

               !>**** multiply BF with BFvec
               do level = 0, level_half
                  call BF_MD_block_extraction_multiply_oneblock_right(blocks, bbm_c, Ndim, BFvec(nn,dim_i), idx_r_m, level , dim_i, ptree,stats)
               enddo

            enddo
         endif

         if(receiver==pp)then
            num_blocks = 2**(level_butterfly-levelm)
            do dim_i=1,Ndim
               allocate (BFvec1(nn,dim_i)%vec(1:level_butterfly-level_half + 1))
               allocate (BFvec1(nn,dim_i)%vec(1)%blocks(num_blocks,1))
               BFvec1(nn,dim_i)%vec(1)%num_row = num_blocks
               BFvec1(nn,dim_i)%vec(1)%num_col = 1
               BFvec1(nn,dim_i)%vec(1)%idx_r = 1
               BFvec1(nn,dim_i)%vec(1)%inc_r = 1
               BFvec1(nn,dim_i)%vec(1)%nr = num_blocks
               BFvec1(nn,dim_i)%vec(1)%idx_c = 1
               BFvec1(nn,dim_i)%vec(1)%inc_c = 1
               BFvec1(nn,dim_i)%vec(1)%nc = 1

               do level = 1, level_butterfly-level_half
                  BFvec1(nn,dim_i)%vec(level + 1)%idx_r = 1
                  BFvec1(nn,dim_i)%vec(level + 1)%inc_r = 1
                  BFvec1(nn,dim_i)%vec(level + 1)%nr = BFvec1(nn,dim_i)%vec(level )%nr/2
                  BFvec1(nn,dim_i)%vec(level + 1)%idx_c = 1
                  BFvec1(nn,dim_i)%vec(level + 1)%inc_c = 1
                  BFvec1(nn,dim_i)%vec(level + 1)%nc = 1
                  BFvec1(nn,dim_i)%vec(level + 1)%num_row = BFvec1(nn,dim_i)%vec(level )%nr/2
                  BFvec1(nn,dim_i)%vec(level + 1)%num_col = 1
                  allocate (BFvec1(nn,dim_i)%vec(level + 1)%blocks(BFvec1(nn,dim_i)%vec(level + 1)%nr, BFvec1(nn,dim_i)%vec(level + 1)%nc))
               enddo

               allocate (BFvec_transposed(nn,dim_i)%vec(1))
               BFvec_transposed(nn,dim_i)%vec(1)%idx_r=1
               BFvec_transposed(nn,dim_i)%vec(1)%inc_r=1
               BFvec_transposed(nn,dim_i)%vec(1)%nr = BFvec1(nn,dim_i)%vec(level_butterfly-level_half + 1)%nr
               call assert(BFvec_transposed(nn,dim_i)%vec(1)%nr==1,"BFvec_transposed(nn,dim_i)%vec(1)%nr should be 1")
               BFvec_transposed(nn,dim_i)%vec(1)%num_row = BFvec1(nn,dim_i)%vec(level_butterfly-level_half + 1)%num_row
               BFvec_transposed(nn,dim_i)%vec(1)%idx_c=1
               BFvec_transposed(nn,dim_i)%vec(1)%inc_c=1
               BFvec_transposed(nn,dim_i)%vec(1)%nc = 1
               BFvec_transposed(nn,dim_i)%vec(1)%num_col = 1
               allocate (BFvec_transposed(nn,dim_i)%vec(1)%blocks(BFvec_transposed(nn,dim_i)%vec(1)%nr,BFvec_transposed(nn,dim_i)%vec(1)%nc))

               !>**** multiply BF with BFvec1
               do level = 0, level_butterfly-level_half
                  call BF_MD_block_extraction_multiply_oneblock_left(blocks, bbm_r, Ndim, BFvec1(nn,dim_i), idx_c_m, level , dim_i, ptree,stats)
               enddo

            enddo
         endif
      enddo

      !>**** Step 2: all to all communication of the right multiplication results
      n2 = MPI_Wtime()
      call BF_MD_all2all_extraction(Ndim, blocks, Ninter, inters, BFvec, BFvec_transposed, stats, msh, ptree)
      n3 = MPI_Wtime()

      !>**** Step 3: tensor matrix contraction
      nc=2**(level_butterfly -level_half)
      dim_subtensor(1:Ndim)=blocks%nr_m
      dim_subtensor(1+Ndim:Ndim*2)=nc

      do nn=1,Ninter
         receiver = inters(nn)%receiver
         if(receiver==pp)then

            idx_r_m = inters(nn)%idx_r_m
            idx_r_m_loc = idx_r_m - blocks%idx_r_m + 1
            idx_c_m = inters(nn)%idx_c_m
            idx_subtensor(1:Ndim) = idx_r_m_loc
            idx_subtensor(1+Ndim:Ndim*2) = idx_c_m
            call MultiIndexToSingleIndex(Ndim*2, dim_subtensor, index_ij, idx_subtensor)

            do dim_i=1,Ndim
               dims_MD_old(dim_i) = size(BFvec1(nn,dim_i)%vec(level_butterfly-level_half + 1)%blocks(1,1)%matrix,2)
            enddo
            do dim_i=1,Ndim
               dims_MD_old(dim_i+Ndim) = size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1)
            enddo


            allocate(mat1(product(dims_MD_old(1:Ndim)),product(dims_MD_old(1+Ndim:2*Ndim))))

#if HAVE_ZFP
            if(option%use_zfp==1 .and. option%use_qtt==0)call ZFP_Decompress(blocks%ButterflyMiddle(index_ij)%matrix,blocks%MiddleZFP(index_ij),product(blocks%ButterflyMiddle(index_ij)%dims_m),product(blocks%ButterflyMiddle(index_ij)%dims_n),tol_used,1)
#endif
            if(option%use_qtt==1)then
               call QTT_Decompress(blocks%MiddleQTT(index_ij), data_buffer)
               mat1 = reshape(data_buffer,[product(dims_MD_old(1:Ndim)),product(dims_MD_old(1+Ndim:2*Ndim))])
               deallocate(data_buffer)
            else
               mat1 = blocks%ButterflyMiddle(index_ij)%matrix
            endif

#if HAVE_ZFP
            if(option%use_zfp==1 .and. option%use_qtt==0)deallocate(blocks%ButterflyMiddle(index_ij)%matrix)
#endif



            call assert(size(blocks%ButterflyMiddle(index_ij)%matrix,1)==product(dims_MD_old(1:Ndim)) .and. size(blocks%ButterflyMiddle(index_ij)%matrix,2)==product(dims_MD_old(1+Ndim:2*Ndim)),"dims_MD_old seems wrong")
            ! write(*,*)dims_MD_old,'dims',shape(mat1),shape(blocks%ButterflyMiddle(index_ij)%matrix),'wori'

            ! looping all the Ndim*2 dimensions and performing contraction along each dimension
            do dim_i=1,Ndim*2
               dims_MD_new = dims_MD_old
               if(dim_i<=Ndim)then
                  dims_MD_new(dim_i) = size(BFvec1(nn,dim_i)%vec(level_butterfly-level_half + 1)%blocks(1,1)%matrix,1)
               else
                  dims_MD_new(dim_i)=size(BFvec_transposed(nn,dim_i-Ndim)%vec(1)%blocks(1,1)%matrix,2)
               endif
               ! write(*,*)"contracting dim",dim_i,"from", dims_MD_old,'to', dims_MD_new

               allocate(mat2(product(dims_MD_new(1:Ndim)),product(dims_MD_new(1+Ndim:2*Ndim))))
               mat2=0


               if(dim_i<=Ndim)then
                  do index_scalar=1,product(dims_MD_new(1:Ndim))
                     call SingleIndexToMultiIndex(Ndim, dims_MD_new(1:Ndim), index_scalar, index_vector)
                     index_vector_old = index_vector
                     do ii=1,dims_MD_old(dim_i)
                        index_vector_old(dim_i) = ii
                        call MultiIndexToSingleIndex(Ndim, dims_MD_old(1:Ndim), index_scalar_old, index_vector_old)
                        mat2(index_scalar,:) = mat2(index_scalar,:) + BFvec1(nn,dim_i)%vec(level_butterfly-level_half + 1)%blocks(1,1)%matrix(index_vector(dim_i),ii)*mat1(index_scalar_old,:)
                     enddo
                  enddo
               else
                  do index_scalar=1,product(dims_MD_new(1+Ndim:2*Ndim))
                     call SingleIndexToMultiIndex(Ndim, dims_MD_new(1+Ndim:2*Ndim), index_scalar, index_vector)
                     index_vector_old = index_vector
                     do ii=1,dims_MD_old(dim_i)
                        index_vector_old(dim_i-Ndim) = ii
                        call MultiIndexToSingleIndex(Ndim, dims_MD_old(1+Ndim:2*Ndim), index_scalar_old, index_vector_old)
                        mat2(:,index_scalar) = mat2(:,index_scalar) + mat1(:,index_scalar_old)*BFvec_transposed(nn,dim_i-Ndim)%vec(1)%blocks(1,1)%matrix(ii,index_vector(dim_i-Ndim))
                     enddo
                  enddo
               endif

               deallocate(mat1)
               allocate(mat1(product(dims_MD_new(1:Ndim)),product(dims_MD_new(1+Ndim:2*Ndim))))
               mat1 = mat2
               deallocate(mat2)

               dims_MD_old = dims_MD_new
            enddo

            inters(nn)%dat = mat1
            deallocate(mat1)

         endif
      enddo


      n4 = MPI_Wtime()

      do nn=1,Ninter
         do dim_i=1,Ndim
            if(allocated(BFvec(nn,dim_i)%vec))then
               do level=1,level_half + 1
                  if(allocated(BFvec(nn,dim_i)%vec(level)%blocks))then
                  do ii=1,BFvec(nn,dim_i)%vec(level)%nr
                     do jj=1,BFvec(nn,dim_i)%vec(level)%nc
                        if(associated(BFvec(nn,dim_i)%vec(level)%blocks(ii,jj)%matrix))deallocate(BFvec(nn,dim_i)%vec(level)%blocks(ii,jj)%matrix)
                     enddo
                  enddo
                  deallocate(BFvec(nn,dim_i)%vec(level)%blocks)
                  endif
               enddo
               deallocate(BFvec(nn,dim_i)%vec)
            endif

            if(allocated(BFvec1(nn,dim_i)%vec))then
               do level=1,level_butterfly-level_half + 1
                  if(allocated(BFvec1(nn,dim_i)%vec(level)%blocks))then
                  do ii=1,BFvec1(nn,dim_i)%vec(level)%nr
                     do jj=1,BFvec1(nn,dim_i)%vec(level)%nc
                        if(associated(BFvec1(nn,dim_i)%vec(level)%blocks(ii,jj)%matrix))deallocate(BFvec1(nn,dim_i)%vec(level)%blocks(ii,jj)%matrix)
                     enddo
                  enddo
                  deallocate(BFvec1(nn,dim_i)%vec(level)%blocks)
                  endif
               enddo
               deallocate(BFvec1(nn,dim_i)%vec)
            endif

            if(allocated(BFvec_transposed(nn,dim_i)%vec))then
               if(allocated(BFvec_transposed(nn,dim_i)%vec(1)%blocks))then
               do ii=1,BFvec_transposed(nn,dim_i)%vec(1)%nr
                  do jj=1,BFvec_transposed(nn,dim_i)%vec(1)%nc
                     if(associated(BFvec_transposed(nn,dim_i)%vec(1)%blocks(ii,jj)%matrix))deallocate(BFvec_transposed(nn,dim_i)%vec(1)%blocks(ii,jj)%matrix)
                  enddo
               enddo
               deallocate(BFvec_transposed(nn,dim_i)%vec(1)%blocks)
               endif
               deallocate(BFvec_transposed(nn,dim_i)%vec)
            endif

         enddo
      enddo
      deallocate(BFvec)
      deallocate(BFvec1)
      deallocate(BFvec_transposed)

      n5 = MPI_Wtime()
      ! time_tmp = time_tmp + n5 - n1
      ! stats%Time_Entry_BF = stats%Time_Entry_BF + n5-n2


   end subroutine BF_MD_block_extraction






!>*********** all to all communication of matvec results at the middle butterfly level from row-wise ordering to column-wise ordering
   subroutine BF_MD_all2all_extraction(Ndim, blocks, Ninter, inters, BFvec, BFvec_transposed, stats, msh, ptree)


      implicit none
      integer Ndim,Ninter
      type(intersect_MD)::inters(Ninter)
      type(mesh)::msh(Ndim)
      type(butterfly_vec) :: BFvec(:,:), BFvec_transposed(:,:)
      integer iijj, i, j, levelm, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, ij, pp, tt,mypp
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn
      real(kind=8) flop
      DT ctemp
      type(matrixblock_MD)::blocks
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
      integer:: dims_r(Ndim),dims_c(Ndim),idx_r_m(Ndim),idx_c_m(Ndim),idx_r_scalar,idx_c_scalar,sender, receiver


      n1 = MPI_Wtime()

      levelm = blocks%level_half
      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno
      dims_r = 2**levelm
      dims_c = 2**(level_butterfly-levelm)
      mypp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1


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

      do nn=1,Ninter
         pp=inters(nn)%sender
         receiver=inters(nn)%receiver
         ! write(*,*)nn,'sender',sender,'receiver',receiver
         if(mypp==receiver)then
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         endif
      enddo


      do nn=1,Ninter
         sender=inters(nn)%sender
         pp=inters(nn)%receiver
         if(mypp==sender)then
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 1 + Ndim*2
            do dim_i=1,Ndim
               sendquant(pp)%size = sendquant(pp)%size + size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)*size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,2)
            enddo
         endif
      enddo

      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid/=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid/=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do nn=1,Ninter
         sender=inters(nn)%sender
         pp=inters(nn)%receiver
         if(mypp==sender)then
            sendquant(pp)%dat(sendquant(pp)%size + 1, 1)=nn
            do dim_i=1,Ndim
               sendquant(pp)%dat(sendquant(pp)%size + 1 + dim_i*2-1, 1) = size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)
               sendquant(pp)%dat(sendquant(pp)%size + 1 + dim_i*2, 1) = size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,2)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + 1 + Ndim*2
            do dim_i=1,Ndim
               do iijj=1,size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)*size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,2)
                  ii = mod(iijj - 1, size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)) + 1
                  jj = ceiling_safe(dble(iijj)/dble(size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)))
                  sendquant(pp)%dat(sendquant(pp)%size + iijj, 1) = BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix(ii,jj)
               enddo
               sendquant(pp)%size = sendquant(pp)%size + size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,1)*size(BFvec(nn,dim_i)%vec(levelm+1)%blocks(1,1)%matrix,2)
            enddo
         endif
      enddo


      ! communicate the data buffer
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         ! write(*,*)ptree%MyID,recvid,'sending?',sendquant(pp)%size
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         ! write(*,*)ptree%MyID,sendid,'receiving?',recvquant(pp)%size
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            nn = NINT(dble(recvquant(pp)%dat(i, 1)))
            do dim_i=1,Ndim
               i = i + 1
               mmm = NINT(dble(recvquant(pp)%dat(i, 1)))
               i = i + 1
               nnn = NINT(dble(recvquant(pp)%dat(i, 1)))
               call assert(.not. associated(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix), 'receiving data alreay exists locally')
               allocate(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix(mmm,nnn))
            enddo

            do dim_i=1,Ndim
               do iijj=1,size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1)*size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,2)
                  ii = mod(iijj - 1, size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1)) + 1
                  jj = ceiling_safe(dble(iijj)/dble(size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1)))
                  BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix(ii,jj) = recvquant(pp)%dat(i+iijj, 1)
               enddo
               i = i + size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1)*size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,2)
            enddo
         enddo
      enddo
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      ! do nn=1,Ninter
      ! do dim_i=1,Ndim
      ! write(*,*)'myID',ptree%MyID,'zala',nn,dim_i,size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,1),size(BFvec_transposed(nn,dim_i)%vec(1)%blocks(1,1)%matrix,2)
      ! enddo
      ! enddo


      ! deallocation
      deallocate (S_req)
      deallocate (R_req)
      deallocate (statuss)
      deallocate (statusr)
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)


      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_MD_all2all_extraction



!>*********** all to all communication of matvec results at the middle butterfly level from row-wise ordering to column-wise ordering
   subroutine BF_MD_all2all_mvp(Ndim, blocks, BFvec, BFvec_transposed, stats, msh, ptree)

      implicit none
      integer Ndim
      type(mesh)::msh(Ndim)
      type(butterfly_vec) :: BFvec(:), BFvec_transposed(:)
      integer iijj, i, j, levelm, level_butterfly, num_blocks, k, attempt, edge_m, edge_n, header_m, header_n, leafsize, nn_start, rankmax_r, rankmax_c, rankmax_min, rank_new
      integer group_m, group_n, mm, nn, index_i, index_i_loc_k, index_i_loc_s, index_j, index_j_loc_k, index_j_loc_s, ii, jj, rr, ij, pp, tt,mypp
      integer level, length_1, length_2, level_blocks
      integer rank, rankmax, butterflyB_inuse, rank1, rank2
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn
      real(kind=8) flop
      DT ctemp
      type(matrixblock_MD)::blocks
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
      integer:: dims_r(Ndim),dims_c(Ndim),idx_r_m(Ndim),idx_c_m(Ndim),idx_r_scalar,idx_c_scalar,sender, receiver


      n1 = MPI_Wtime()

      levelm = blocks%level_half
      level_butterfly = blocks%level_butterfly
      nproc = ptree%pgrp(blocks%pgno)%nproc
      tag = blocks%pgno
      dims_r = 2**levelm
      dims_c = 2**(level_butterfly-levelm)
      mypp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1

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
      do bb=1, product(blocks%nr_m)
         do jj = 1, BFvec_transposed(bb)%vec(1)%nc
            pp = proc_of_groupc(jj)
            if (recvquant(pp)%active == 0) then
               recvquant(pp)%active = 1
               Nrecvactive = Nrecvactive + 1
               recvIDactive(Nrecvactive) = pp
            endif
         enddo
      enddo

      do bb=1, product(blocks%nc_m)
         do ii = 1, BFvec(bb)%vec(levelm+1)%nr
            pp = proc_of_groupr(ii)
            if (sendquant(pp)%active == 0) then
               sendquant(pp)%active = 1
               Nsendactive = Nsendactive + 1
               sendIDactive(Nsendactive) = pp
            endif
            sendquant(pp)%size = sendquant(pp)%size + 4
            sendquant(pp)%size = sendquant(pp)%size + size(BFvec(bb)%vec(levelm+1)%blocks(ii,1)%matrix,1)*size(BFvec(bb)%vec(levelm+1)%blocks(ii,1)%matrix,2)
         enddo
      enddo


      ! communicate receive buffer sizes
      Nreqs=0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(recvid/=ptree%MyID)then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            recvquant(pp)%size = sendquant(pp)%size
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         if(sendid/=ptree%MyID)then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%size, 1, MPI_INTEGER, pp - 1, tag, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
      enddo

      ! pack the send buffer in the second pass
      do bb=1, product(blocks%nc_m)
         do rr = 1, BFvec(bb)%vec(levelm+1)%nr
            pp = proc_of_groupr(rr)
            call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb, idx_c_m)
            idx_c_m = idx_c_m + blocks%idx_c_m - 1
            call MultiIndexToSingleIndex(Ndim, dims_c, idx_c_scalar, idx_c_m)
            sendquant(pp)%dat(sendquant(pp)%size + 1, 1) = rr ! global row index at the middle level
            sendquant(pp)%dat(sendquant(pp)%size + 2, 1) = idx_c_scalar ! global column index at the middle level
            sendquant(pp)%dat(sendquant(pp)%size + 3, 1) = size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,1)
            sendquant(pp)%dat(sendquant(pp)%size + 4, 1) = size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,2)
            sendquant(pp)%size = sendquant(pp)%size + 4

            do iijj=1,size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,1)*size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,2)
               ii = mod(iijj - 1, size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,1)) + 1
               jj = ceiling_safe(dble(iijj)/dble(size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,1)))
               sendquant(pp)%dat(sendquant(pp)%size + iijj, 1) = BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix(ii,jj)
            enddo
            sendquant(pp)%size = sendquant(pp)%size + size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,1)*size(BFvec(bb)%vec(levelm+1)%blocks(rr,1)%matrix,2)
         enddo
      enddo

      ! communicate the data buffer
      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         ! write(*,*)ptree%MyID,recvid,'sending?',sendquant(pp)%size
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, S_req(Nreqs), ierr)
         else
            if (sendquant(pp)%size > 0) recvquant(pp)%dat = sendquant(pp)%dat
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(blocks%pgno)%head
         ! write(*,*)ptree%MyID,sendid,'receiving?',recvquant(pp)%size
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(blocks%pgno)%Comm, R_req(Nreqr), ierr)
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
         do while (i < recvquant(pp)%size)
            i = i + 1
            idx_r_scalar = NINT(dble(recvquant(pp)%dat(i, 1)))
            call SingleIndexToMultiIndex(Ndim, dims_r, idx_r_scalar, idx_r_m)
            idx_r_m = idx_r_m - blocks%idx_r_m + 1
            call MultiIndexToSingleIndex(Ndim, blocks%nr_m, bb, idx_r_m)
            i = i + 1
            rr = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            mmm = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            nnn = NINT(dble(recvquant(pp)%dat(i, 1)))
            call assert(.not. associated(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix), 'receiving data alreay exists locally')
            allocate(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix(mmm,nnn))
            do iijj=1,size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,1)*size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,2)
               ii = mod(iijj - 1, size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,1)) + 1
               jj = ceiling_safe(dble(iijj)/dble(size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,1)))
               BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix(ii,jj) = recvquant(pp)%dat(i+iijj, 1)
            enddo
            i = i + size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,1)*size(BFvec_transposed(bb)%vec(1)%blocks(1,rr)%matrix,2)
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
         if (allocated(sendquant(pp)%dat)) deallocate (sendquant(pp)%dat)
      enddo
      deallocate (sendquant)
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat)) deallocate (recvquant(pp)%dat)
      enddo
      deallocate (recvquant)
      deallocate (sendIDactive)
      deallocate (recvIDactive)
      deallocate (proc_of_groupr)
      deallocate (proc_of_groupc)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

   end subroutine BF_MD_all2all_mvp


   subroutine BF_MD_block_mvp(chara, xin, Ninloc, xout, Noutloc, Nvec, blocks, Ndim, ptree, stats,msh,option)

      implicit none
      integer Ndim,dim_i,bb_m,num_threads
      character chara
      integer Ninloc(Ndim), Noutloc(Ndim), Nvec
      DT::xin(:, :), xout(:, :)
      DT,allocatable::xout1(:, :),data_buffer(:)
      type(matrixblock_MD)::blocks
      type(proctree)::ptree
      type(Hstat)::stats
      type(mesh)::msh(Ndim)
      type(Hoption)::option
      integer ri, ci, nng, headm, headn, head, tail, pp, row_group, col_group
      integer k1, i, j, ii, iii,iii1, iii2, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc0, nr0, inc_c, inc_c0, inc_r, inc_r0
      integer level_butterfly, num_blocks, level_half, levelm, nr(Ndim), nc(Ndim)
      integer mm, nn, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1(Ndim), nvec2(Ndim)
      integer level
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      integer nsendrecv, pid, pgno_sub, pid0, tag, nproc, Nreqr, Nreqs, recvid, sendid
      integer idx_r_m(Ndim),idx_c_m(Ndim),idx_r_m_loc(Ndim),idx_c_m_loc(Ndim)
      type(butterfly_vec),allocatable :: BFvec(:), BFvec1(:), BFvec_transposed(:)
      type(butterfly_kerl), allocatable::g_idx_m(:), g_idx_n(:)
      integer, allocatable::group_ms(:), group_ns(:), group_ms1(:), group_ns1(:)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(ipair):: p
      DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)
      DT, allocatable::Vpartial(:, :)
      DT::val
      type(vectorsblock),allocatable::mats(:),mats1(:)
      type(vectorset),allocatable::mats1_1D(:)
      integer sender, receiver
      integer, allocatable::num_nods_i(:), num_nods_j(:)
      real(kind=8)::n1, n2, n3, n4, n5, n6, n7, t1, t2
      integer dims_c_m(Ndim),dims_r_m(Ndim),bbm_r,bbm_c,group_n(Ndim),group_n1(Ndim),group_m(Ndim),group_m1(Ndim), dims_MD_old(Ndim*2),dims_MD_new(Ndim*2),dims_ref(Ndim+1),dims_ref_old(Ndim+1),dims_ref_new(Ndim+1),offsets_ref(Ndim+1), dims_in(Ndim),idx_in(Ndim), dims_MD3(Ndim),idx_MD3(Ndim), idxr(Ndim), dim_subtensor(Ndim*2),idx_subtensor(Ndim*2), index_ij, index_scalar,index_vector(Ndim),index_scalar_old,index_vector_old(Ndim)
      real(kind=8)::flop,flops,tol_used
      integer, save:: my_tid = 0
#ifdef HAVE_OPENMP
      !$omp threadprivate(my_tid)
#endif

#ifdef HAVE_OPENMP
!$omp parallel default(shared)
!$omp master
      num_threads = omp_get_num_threads()
!$omp end master
      my_tid = omp_get_thread_num()
!$omp end parallel
#else
      num_threads = 1
      my_tid = 0
#endif


      stats%Flop_Tmp = 0
      stats%Flop_C_Mult = 0
      stats%Time_C_Mult = 0

      n1 = MPI_Wtime()
      t1 = MPI_Wtime()

      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm

      pp = ptree%myid - ptree%pgrp(pgno)%head + 1


      level_butterfly = blocks%level_butterfly
      level_half = blocks%level_half
      levelm = level_half
      dims_c_m = 2**(level_butterfly-levelm)
      dims_r_m = 2**levelm


      ! write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,'inin'
      if(chara=='N')then

         !>******* preallocate BFvec and BFvec1, number of blocks are exactly those in BF_block_MVP_dat, BFvec for the first 0:level_half+1 levels and BFvec for the next level_half+1 to level_butterfly+2 levels,note that level level_half+1 is duplicated for all2all communication
         allocate(BFvec(product(blocks%nc_m)))
         allocate(BFvec1(product(blocks%nr_m)))

         allocate(BFvec_transposed(product(blocks%nr_m)))

         group_n1 = blocks%col_group
         group_n1 = group_n1*2**(level_butterfly-levelm) + blocks%idx_c_m - 1

         ! write(*,*)ptree%MyID,blocks%idx_r_m,blocks%idx_c_m,'wogan'
         do nn=1,product(blocks%nc_m)
            call SingleIndexToMultiIndex(Ndim, blocks%nc_m, nn, idx_c_m)
            group_n = blocks%col_group
            group_n = group_n*2**(level_butterfly-levelm) - 1 + idx_c_m + blocks%idx_c_m - 1

            num_blocks = 2**levelm

            allocate (BFvec(nn)%vec(0:level_half + 1))

            allocate (BFvec(nn)%vec(0)%blocks(1, 1))
            allocate (BFvec(nn)%vec(0)%index_MD(Ndim,1,num_blocks+1))
            BFvec(nn)%vec(0)%num_row = 1
            BFvec(nn)%vec(0)%num_col = num_blocks
            BFvec(nn)%vec(0)%idx_r = 1
            BFvec(nn)%vec(0)%inc_r = 1
            BFvec(nn)%vec(0)%nr = 1
            BFvec(nn)%vec(0)%idx_c = 1
            BFvec(nn)%vec(0)%inc_c = 1
            BFvec(nn)%vec(0)%nc = num_blocks
            do dim_i=1,Ndim
               k1 = 0
               do i = 1, num_blocks
                  BFvec(nn)%vec(0)%index_MD(dim_i,1,i) = k1
                  nnn = size(blocks%ButterflyV(nn)%blocks(i,dim_i)%matrix, 1)
                  k1 = k1 + nnn
               enddo
               BFvec(nn)%vec(0)%index_MD(dim_i,1,num_blocks+1)=k1
            enddo


            dims_ref_old(1:Ndim) = Ninloc
            dims_ref_old(Ndim+1) = Nvec
            do dim_i=1,Ndim
               dims_ref_new(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
               offsets_ref(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%head - msh(dim_i)%basis_group(group_n1(dim_i))%head
            enddo
            dims_ref_new(Ndim+1) = Nvec
            offsets_ref(Ndim+1) = 0
            call TensorUnfoldingReshape(Ndim+1,dims_ref_old,dims_ref_new,offsets_ref,Ndim+1,1,xin, 'T', mat,'N',1)
            allocate(BFvec(nn)%vec(0)%blocks(1,1)%matrix(size(mat,1),size(mat,2)))
            BFvec(nn)%vec(0)%blocks(1,1)%matrix = mat
            deallocate(mat)

            allocate (BFvec(nn)%vec(1)%blocks(1, 1))
            allocate (BFvec(nn)%vec(1)%index_MD(Ndim,1,num_blocks+1))
            BFvec(nn)%vec(1)%num_row = 1
            BFvec(nn)%vec(1)%num_col = num_blocks
            BFvec(nn)%vec(1)%idx_r = 1
            BFvec(nn)%vec(1)%inc_r = 1
            BFvec(nn)%vec(1)%nr = 1
            BFvec(nn)%vec(1)%idx_c = 1
            BFvec(nn)%vec(1)%inc_c = 1
            BFvec(nn)%vec(1)%nc = num_blocks
            do dim_i=1,Ndim
               k1 = 0
               do i = 1, num_blocks
                  BFvec(nn)%vec(1)%index_MD(dim_i,1,i) = k1
                  nnn = size(blocks%ButterflyV(nn)%blocks(i,dim_i)%matrix, 2)
                  k1 = k1 + nnn
               enddo
               BFvec(nn)%vec(1)%index_MD(dim_i,1,num_blocks+1) = k1
            enddo

            do level = 1, level_half
               nr=2**(level)
               BFvec(nn)%vec(level + 1)%idx_r = 1
               BFvec(nn)%vec(level + 1)%inc_r = 1
               BFvec(nn)%vec(level + 1)%nr = product(nr)
               BFvec(nn)%vec(level + 1)%idx_c = 1
               BFvec(nn)%vec(level + 1)%inc_c = 1
               BFvec(nn)%vec(level + 1)%nc = BFvec(nn)%vec(level)%nc/2
               BFvec(nn)%vec(level + 1)%num_row = product(nr)
               BFvec(nn)%vec(level + 1)%num_col = BFvec(nn)%vec(level)%nc/2
               allocate (BFvec(nn)%vec(level + 1)%blocks(BFvec(nn)%vec(level + 1)%nr,1))
               allocate (BFvec(nn)%vec(level + 1)%index_MD(Ndim, product(nr), BFvec(nn)%vec(level + 1)%nc+1))
               do dim_i=1,Ndim
                  do i=1,BFvec(nn)%vec(level + 1)%nr
                     k1 = 0
                     do j = 1, BFvec(nn)%vec(level + 1)%nc
                        BFvec(nn)%vec(level + 1)%index_MD(dim_i,i,j) = k1
                        nnn = size(blocks%ButterflyKerl_R(nn,level)%blocks(i, 2*j-1, dim_i)%matrix, 1)
                        k1 = k1 + nnn
                     enddo
                     BFvec(nn)%vec(level + 1)%index_MD(dim_i,i, BFvec(nn)%vec(level + 1)%nc+1) = k1
                  enddo
               enddo
            enddo
         enddo

         group_m1 = blocks%row_group
         group_m1 = group_m1*2**levelm + blocks%idx_r_m - 1
         do mm=1,product(blocks%nr_m)
            call SingleIndexToMultiIndex(Ndim, blocks%nr_m, mm, idx_r_m)
            group_m = blocks%row_group
            group_m = group_m*2**levelm - 1 + idx_r_m + blocks%idx_r_m - 1

            num_blocks = 2**(level_butterfly-levelm)

            allocate (BFvec1(mm)%vec(0:level_butterfly-level_half + 1))
            allocate (BFvec1(mm)%vec(0)%blocks(1, 1))
            allocate (BFvec1(mm)%vec(0)%index_MD(Ndim,num_blocks+1,1))

            BFvec1(mm)%vec(0)%num_row = num_blocks
            BFvec1(mm)%vec(0)%num_col = 1
            BFvec1(mm)%vec(0)%idx_r = 1
            BFvec1(mm)%vec(0)%inc_r = 1
            BFvec1(mm)%vec(0)%nr = num_blocks
            BFvec1(mm)%vec(0)%idx_c = 1
            BFvec1(mm)%vec(0)%inc_c = 1
            BFvec1(mm)%vec(0)%nc = 1
            do dim_i=1,Ndim
               k1 = 0
               do i = 1, num_blocks
                  BFvec1(mm)%vec(0)%index_MD(dim_i,i,1) = k1
                  mmm = size(blocks%ButterflyU(mm)%blocks(i,dim_i)%matrix, 1)
                  k1 = k1 + mmm
               enddo
               BFvec1(mm)%vec(0)%index_MD(dim_i,num_blocks+1,1)=k1
            enddo

            allocate (BFvec1(mm)%vec(1)%blocks(1, 1))
            allocate (BFvec1(mm)%vec(1)%index_MD(Ndim,num_blocks+1,1))
            BFvec1(mm)%vec(1)%num_row = num_blocks
            BFvec1(mm)%vec(1)%num_col = 1
            BFvec1(mm)%vec(1)%idx_r = 1
            BFvec1(mm)%vec(1)%inc_r = 1
            BFvec1(mm)%vec(1)%nr = num_blocks
            BFvec1(mm)%vec(1)%idx_c = 1
            BFvec1(mm)%vec(1)%inc_c = 1
            BFvec1(mm)%vec(1)%nc = 1
            do dim_i=1,Ndim
               k1 = 0
               do i = 1, num_blocks
                  BFvec1(mm)%vec(1)%index_MD(dim_i,i,1) = k1
                  mmm = size(blocks%ButterflyU(mm)%blocks(i,dim_i)%matrix, 2)
                  k1 = k1 + mmm
               enddo
               BFvec1(mm)%vec(1)%index_MD(dim_i,num_blocks+1,1) = k1
            enddo

            do level = 1, level_butterfly-level_half
               nc=2**(level)
               BFvec1(mm)%vec(level + 1)%idx_r = 1
               BFvec1(mm)%vec(level + 1)%inc_r = 1
               BFvec1(mm)%vec(level + 1)%nc = product(nc)
               BFvec1(mm)%vec(level + 1)%idx_c = 1
               BFvec1(mm)%vec(level + 1)%inc_c = 1
               BFvec1(mm)%vec(level + 1)%nr = BFvec1(mm)%vec(level)%nr/2
               BFvec1(mm)%vec(level + 1)%num_col = product(nc)
               BFvec1(mm)%vec(level + 1)%num_row = BFvec1(mm)%vec(level)%nr/2
               allocate (BFvec1(mm)%vec(level + 1)%blocks(1,BFvec1(mm)%vec(level + 1)%nc))
               allocate (BFvec1(mm)%vec(level + 1)%index_MD(Ndim, BFvec1(mm)%vec(level + 1)%nr+1, product(nc)))
               do dim_i=1,Ndim
                  do j=1,BFvec1(mm)%vec(level + 1)%nc
                     k1 = 0
                     do i = 1, BFvec1(mm)%vec(level + 1)%nr
                        BFvec1(mm)%vec(level + 1)%index_MD(dim_i,i,j) = k1
                        nnn = size(blocks%ButterflyKerl_L(mm,level_butterfly-level+1)%blocks(2*i-1, j, dim_i)%matrix, 2)
                        k1 = k1 + nnn
                     enddo
                     BFvec1(mm)%vec(level + 1)%index_MD(dim_i,BFvec1(mm)%vec(level + 1)%nr+1,j) = k1
                  enddo
               enddo
            enddo
            allocate (BFvec_transposed(mm)%vec(1))
            level = level_butterfly-level_half
            nc=2**(level)
            BFvec_transposed(mm)%vec(1)%idx_r = 1
            BFvec_transposed(mm)%vec(1)%inc_r = 1
            BFvec_transposed(mm)%vec(1)%nc = product(nc)
            BFvec_transposed(mm)%vec(1)%idx_c = 1
            BFvec_transposed(mm)%vec(1)%inc_c = 1
            BFvec_transposed(mm)%vec(1)%nr = BFvec1(mm)%vec(level_butterfly-level_half + 1)%num_row
            BFvec_transposed(mm)%vec(1)%num_col = product(nc)
            BFvec_transposed(mm)%vec(1)%num_row = BFvec1(mm)%vec(level_butterfly-level_half + 1)%nr
            allocate (BFvec_transposed(mm)%vec(1)%blocks(1,BFvec_transposed(mm)%vec(1)%nc))
            allocate (BFvec_transposed(mm)%vec(1)%index_MD(Ndim, BFvec_transposed(mm)%vec(1)%nr+1, product(nc)))
            do dim_i=1,Ndim
               do j=1,BFvec_transposed(mm)%vec(1)%nc
                  k1 = 0
                  do i = 1, BFvec_transposed(mm)%vec(1)%nr
                     BFvec_transposed(mm)%vec(1)%index_MD(dim_i,i,j) = k1
                     if(level_butterfly>0)then
                        nnn = size(blocks%ButterflyKerl_L(mm,level_butterfly-level+1)%blocks(2*i-1, j, dim_i)%matrix, 2)
                     else
                        nnn = size(blocks%ButterflyU(mm)%blocks(i,dim_i)%matrix, 2)
                     endif
                     k1 = k1 + nnn
                  enddo
                  BFvec_transposed(mm)%vec(1)%index_MD(dim_i,BFvec_transposed(mm)%vec(1)%nr+1,j) = k1
               enddo
            enddo
         enddo

         n2 = MPI_Wtime()
         !>**** Step 1: multiply the factor matrices out from outtermost to middle level along each dimension
#ifdef HAVE_TASKLOOP
         !$omp parallel
         !$omp single
#endif
         do nn=1,product(blocks%nc_m)
#ifdef HAVE_TASKLOOP
!$omp task default(shared) private(level) firstprivate(nn)
#endif
            do level = 0, level_half
               call BF_MD_block_mvp_multiply_right(blocks, nn, Ndim, BFvec(nn), Nvec, level, ptree,stats)
            enddo
#ifdef HAVE_TASKLOOP
!$omp end task
#endif
            ! write(*,*)nn,product(blocks%nc_m),ptree%MyID,'done'
         enddo
#ifdef HAVE_TASKLOOP
         !$omp end single
         !$omp end parallel
#endif
         !>**** Step 2: all to all communication of the right multiplication results
         n3 = MPI_Wtime()
         call BF_MD_all2all_mvp(Ndim, blocks, BFvec, BFvec_transposed, stats, msh, ptree)
         n4 = MPI_Wtime()

         !>**** Step 3: contraction with the middle level core tensors
         dim_subtensor(1:Ndim) = blocks%nr_m
         dim_subtensor(1+Ndim:2*Ndim) = 2**(level_butterfly-level_half)
         allocate(mats(num_threads))
         allocate(mats1(num_threads))
         allocate(mats1_1D(num_threads))

         flops = 0
#ifdef HAVE_TASKLOOP
         !$omp parallel
         !$omp single
         !$omp taskloop default(shared) private(index_ij,idx_subtensor,mm,j,dims_ref,offsets_ref,tol_used) reduction(+:flops)
#else
#ifdef HAVE_OPENMP
         !$omp parallel do default(shared) private(index_ij,idx_subtensor,mm,j,dims_ref,offsets_ref,tol_used) reduction(+:flops)
#endif
#endif
         do index_ij = 1,product(dim_subtensor)
            call SingleIndexToMultiIndex(Ndim*2, dim_subtensor, index_ij, idx_subtensor)
            call MultiIndexToSingleIndex(Ndim,dim_subtensor(1:Ndim),mm,idx_subtensor(1:Ndim))
            call MultiIndexToSingleIndex(Ndim,dim_subtensor(1+Ndim:2*Ndim),j,idx_subtensor(1+Ndim:Ndim*2))

            dims_ref(1:Ndim) = blocks%ButterflyMiddle(index_ij)%dims_n
            dims_ref(1+Ndim) = Nvec
            offsets_ref = 0
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,Ndim,Ndim+1,BFvec_transposed(mm)%vec(1)%blocks(1,j)%matrix, 'N', mats(my_tid+1)%vector, 'T',1)

            allocate(mats1(my_tid+1)%vector(product(blocks%ButterflyMiddle(index_ij)%dims_m),Nvec))
            mats1(my_tid+1)%vector=0

#if HAVE_ZFP
            if(allocated(blocks%MiddleZFP))then
            if(allocated(blocks%MiddleZFP(index_ij)%buffer_r))then
               call ZFP_Decompress(blocks%ButterflyMiddle(index_ij)%matrix,blocks%MiddleZFP(index_ij),product(blocks%ButterflyMiddle(index_ij)%dims_m),product(blocks%ButterflyMiddle(index_ij)%dims_n),tol_used,1)
            endif
            endif
#endif

            if(allocated(blocks%MiddleQTT))then

               allocate(mats1_1D(my_tid+1)%vector(product(blocks%ButterflyMiddle(index_ij)%dims_m)*Nvec))
               mats1_1D(my_tid+1)%vector=0
               call QTT_Apply_Fullvec(blocks%MiddleQTT(index_ij),reshape(mats(my_tid+1)%vector,[product(blocks%ButterflyMiddle(index_ij)%dims_n)*Nvec]),mats1_1D(my_tid+1)%vector)
               mats1(my_tid+1)%vector = reshape(mats1_1D(my_tid+1)%vector,[product(blocks%ButterflyMiddle(index_ij)%dims_m),Nvec])
               deallocate(mats1_1D(my_tid+1)%vector)
               ! write(*,*)'fnorm of output (QTT_Apply):', sqrt(sum(abs(mats1(my_tid+1)%vector)**2d0))


               ! mats1(my_tid+1)%vector=0
               ! call QTT_Decompress(blocks%MiddleQTT(index_ij), mats1_1D(my_tid+1)%vector)
               ! allocate(blocks%ButterflyMiddle(index_ij)%matrix(product(blocks%ButterflyMiddle(index_ij)%dims_m),product(blocks%ButterflyMiddle(index_ij)%dims_n)))
               ! blocks%ButterflyMiddle(index_ij)%matrix = reshape(mats1_1D(my_tid+1)%vector,[product(blocks%ButterflyMiddle(index_ij)%dims_m),product(blocks%ButterflyMiddle(index_ij)%dims_n)])
               ! deallocate(mats1_1D(my_tid+1)%vector)
               ! call gemmf90(blocks%ButterflyMiddle(index_ij)%matrix, product(blocks%ButterflyMiddle(index_ij)%dims_m), mats(my_tid+1)%vector, product(blocks%ButterflyMiddle(index_ij)%dims_n), mats1(my_tid+1)%vector, product(blocks%ButterflyMiddle(index_ij)%dims_m), 'N', 'N', product(blocks%ButterflyMiddle(index_ij)%dims_m), Nvec, product(blocks%ButterflyMiddle(index_ij)%dims_n), BPACK_cone, BPACK_czero, flop)
               ! deallocate(blocks%ButterflyMiddle(index_ij)%matrix)
               ! write(*,*)'fnorm of output (gemm):', sqrt(sum(abs(mats1(my_tid+1)%vector)**2d0))




            else
               call gemmf90(blocks%ButterflyMiddle(index_ij)%matrix, product(blocks%ButterflyMiddle(index_ij)%dims_m), mats(my_tid+1)%vector, product(blocks%ButterflyMiddle(index_ij)%dims_n), mats1(my_tid+1)%vector, product(blocks%ButterflyMiddle(index_ij)%dims_m), 'N', 'N', product(blocks%ButterflyMiddle(index_ij)%dims_m), Nvec, product(blocks%ButterflyMiddle(index_ij)%dims_n), BPACK_cone, BPACK_czero, flop)
            endif
            flops = flops + flop
#if HAVE_ZFP
            if(allocated(blocks%MiddleZFP))then
            if(allocated(blocks%MiddleZFP(index_ij)%buffer_r))then
               call ZFP_Compress(blocks%ButterflyMiddle(index_ij)%matrix,blocks%MiddleZFP(index_ij),product(blocks%ButterflyMiddle(index_ij)%dims_m),product(blocks%ButterflyMiddle(index_ij)%dims_n),tol_used,1)
            endif
            endif
#endif
            dims_ref(1:Ndim) = blocks%ButterflyMiddle(index_ij)%dims_m
            dims_ref(1+Ndim) = Nvec
            offsets_ref = 0
            deallocate(mats(my_tid+1)%vector)
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,Ndim+1,Ndim,mats1(my_tid+1)%vector, 'T', mats(my_tid+1)%vector, 'N',1)
            allocate(BFvec1(mm)%vec(level_butterfly-level_half + 1)%blocks(1,j)%matrix(size(mats(my_tid+1)%vector,1),size(mats(my_tid+1)%vector,2)))
            BFvec1(mm)%vec(level_butterfly-level_half + 1)%blocks(1,j)%matrix = mats(my_tid+1)%vector
            deallocate(mats1(my_tid+1)%vector)
            deallocate(mats(my_tid+1)%vector)
         enddo
#ifdef HAVE_TASKLOOP
         !$omp end taskloop
         !$omp end single
         !$omp end parallel
#else
#ifdef HAVE_OPENMP
         !$omp end parallel do
#endif
#endif
         stats%Flop_Tmp = stats%Flop_Tmp + flops
         deallocate(mats)
         deallocate(mats1)

         n5 = MPI_Wtime()

         !>**** Step 4: multiply the factor matrices out from middle level to the outtermost level along each dimension
         allocate(xout1(size(xout,1),size(xout,2)))
         xout1 = 0
#ifdef HAVE_TASKLOOP
         !$omp parallel
         !$omp single
#endif
         do mm=1,product(blocks%nr_m)
#ifdef HAVE_TASKLOOP
!$omp task default(shared) private(level,group_m,idx_r_m,dims_ref_new,dim_i,dims_ref_old,offsets_ref) firstprivate(mm)
#endif
            call SingleIndexToMultiIndex(Ndim, blocks%nr_m, mm, idx_r_m)
            group_m = blocks%row_group
            group_m = group_m*2**levelm - 1 + idx_r_m + blocks%idx_r_m - 1

            do level = level_butterfly-level_half, 0, -1
               call BF_MD_block_mvp_multiply_left(blocks, mm, Ndim, BFvec1(mm), Nvec, level, ptree,stats)
            enddo

            dims_ref_new(1:Ndim) = Noutloc
            dims_ref_new(Ndim+1) = Nvec

            do dim_i=1,Ndim
               dims_ref_old(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
               offsets_ref(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%head - msh(dim_i)%basis_group(group_m1(dim_i))%head
            enddo
            dims_ref_old(Ndim+1) = Nvec
            offsets_ref(Ndim+1) = 0

            call TensorUnfoldingReshape(Ndim+1,dims_ref_old,dims_ref_new,offsets_ref,Ndim,Ndim+1,BFvec1(mm)%vec(0)%blocks(1, 1)%matrix, 'N', xout1,'T',0)
#ifdef HAVE_TASKLOOP
!$omp end task
#endif
            ! write(*,*)mm,product(blocks%nr_m),ptree%MyID,'done L'
         enddo
#ifdef HAVE_TASKLOOP
         !$omp end single
         !$omp end parallel
#endif
         xout=xout1
         deallocate(xout1)

         n6 = MPI_Wtime()

         do nn=1,product(blocks%nc_m)
            if(allocated(BFvec(nn)%vec))then
               do level=0,level_half + 1
                  if(allocated(BFvec(nn)%vec(level)%blocks))then
                  do ii=1,BFvec(nn)%vec(level)%nr
                     if(associated(BFvec(nn)%vec(level)%blocks(ii,1)%matrix))deallocate(BFvec(nn)%vec(level)%blocks(ii,1)%matrix)
                  enddo
                  deallocate(BFvec(nn)%vec(level)%blocks)
                  deallocate(BFvec(nn)%vec(level)%index_MD)
                  endif
               enddo
               deallocate(BFvec(nn)%vec)
            endif
         enddo

         do nn=1,product(blocks%nr_m)
            if(allocated(BFvec1(nn)%vec))then
               do level=0,level_butterfly-level_half + 1
                  if(allocated(BFvec1(nn)%vec(level)%blocks))then
                  do jj=1,BFvec1(nn)%vec(level)%nc
                     if(associated(BFvec1(nn)%vec(level)%blocks(1,jj)%matrix))deallocate(BFvec1(nn)%vec(level)%blocks(1,jj)%matrix)
                  enddo
                  deallocate(BFvec1(nn)%vec(level)%blocks)
                  deallocate(BFvec1(nn)%vec(level)%index_MD)
                  endif
               enddo
               deallocate(BFvec1(nn)%vec)
            endif

            if(allocated(BFvec_transposed(nn)%vec))then
               if(allocated(BFvec_transposed(nn)%vec(1)%blocks))then
               do jj=1,BFvec_transposed(nn)%vec(1)%nc
                  if(associated(BFvec_transposed(nn)%vec(1)%blocks(1,jj)%matrix))deallocate(BFvec_transposed(nn)%vec(1)%blocks(1,jj)%matrix)
               enddo
               deallocate(BFvec_transposed(nn)%vec(1)%blocks)
               deallocate(BFvec_transposed(nn)%vec(1)%index_MD)
               endif
               deallocate(BFvec_transposed(nn)%vec)
            endif
         enddo

         deallocate(BFvec)
         deallocate(BFvec1)
         deallocate(BFvec_transposed)
   elseif(chara=='T')then
      write(*,*)'chara==T not yet implemented in BF_MD_block_mvp'
   endif
      n7 = MPI_Wtime()
      ! time_tmp = time_tmp + n5 - n1
      if(ptree%MyID==Main_ID .and. option%verbosity>=2)write(*,*)'BF_MD_block_mvp: MyID',ptree%MyID,'allocate',n2-n1,'right_mult',n3-n2,'all2all',n4-n3,'tensor_contract',n5-n4,'left_mult',n6-n5,'deallocate',n7-n6

      t2 = MPI_Wtime()
      stats%Time_C_Mult = stats%Time_C_Mult + t2 - t1
      stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

   end subroutine BF_MD_block_mvp








   subroutine BF_block_extraction(blocks, inters, ptree, msh, stats)

      implicit none
      type(matrixblock)::blocks
      type(proctree)::ptree
      type(mesh)::msh
      type(Hstat)::stats
      integer ri, ci, nng, headm, headn, head, tail, pp, row_group, col_group
      type(intersect)::inters(:)
      integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
      integer level_butterfly, num_blocks, level_half
      integer group_m, group_n, mm, nn, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
      integer level
      real(kind=8) rate, tolerance, rtemp, norm_1, norm_2, norm_e
      integer header_n1, header_n2, nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
      integer nsendrecv, pid, pgno_sub, pid0, tag, nproc, Nreqr, Nreqs, recvid, sendid
      type(butterfly_vec) :: BFvec, BFvec1
      type(butterfly_kerl), allocatable::g_idx_m(:), g_idx_n(:)
      integer, allocatable::group_ms(:), group_ns(:), group_ms1(:), group_ns1(:)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(ipair):: p
      DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)
      DT, allocatable::Vpartial(:, :)
      DT::val
      integer idxr, idxc
      integer, allocatable::num_nods_i(:), num_nods_j(:)
      real(kind=8)::n1, n2, n3, n4, n5


      n1 = MPI_Wtime()

      headm = blocks%headm
      headn = blocks%headn
      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm

      pp = ptree%myid - ptree%pgrp(pgno)%head + 1
      head = blocks%N_p(pp, 1)
      tail = blocks%N_p(pp, 2)

      level_butterfly = blocks%level_butterfly
      num_blocks = 2**level_butterfly
      level_half = blocks%level_half

      ! write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,'inin'

      !>******* preallocate BFvec and BFvec1, number of blocks are exactly those in BF_block_MVP_dat, BFvec for the first 0:level_half+1 levels and BFvec for the next level_half+1 to level_butterfly+2 levels,note that level level_half+1 is duplicated for all2all communication
      allocate (BFvec%vec(0:level_half + 1))
      allocate (BFvec%vec(0)%blocks(1, blocks%ButterflyV%nblk_loc))
      BFvec%vec(0)%num_row = 1
      BFvec%vec(0)%num_col = num_blocks
      BFvec%vec(0)%idx_r = 1
      BFvec%vec(0)%inc_r = 1
      BFvec%vec(0)%nr = 1
      BFvec%vec(0)%idx_c = blocks%ButterflyV%idx
      BFvec%vec(0)%inc_c = blocks%ButterflyV%inc
      BFvec%vec(0)%nc = blocks%ButterflyV%nblk_loc

      do level = 0, level_half
         call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

         BFvec%vec(level + 1)%idx_r = idx_r
         BFvec%vec(level + 1)%inc_r = inc_r
         BFvec%vec(level + 1)%nr = nr
         BFvec%vec(level + 1)%idx_c = idx_c
         BFvec%vec(level + 1)%inc_c = inc_c
         BFvec%vec(level + 1)%nc = nc
         if (level /= level_butterfly + 1) then
            BFvec%vec(level + 1)%num_row = 2**level
            BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
         else
            BFvec%vec(level + 1)%num_row = 2**level_butterfly
            BFvec%vec(level + 1)%num_col = 1
         endif
         if (level_half /= level) then ! the last level doesn't require doubling block columns
         if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
            BFvec%vec(level + 1)%nc = 2
            BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
         endif
         endif
         allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))
      enddo

      allocate (BFvec1%vec(level_half + 1:level_butterfly + 2))
      do level = level_half, level_butterfly + 1
         if (level == level_half) then
            call GetLocalBlockRange(ptree, blocks%pgno, level_half + 1, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C') ! this is the same as the target block range used in BF_all2all_vec_n_ker
            idx_r = idx_r0
            nr = nr0
            inc_r = inc_r0
            idx_c = idx_c0
            nc = nc0
            inc_c = inc_c0
         else
            call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C')
            ! convert the local column-wise kernel block ranges to local row-wise output vector ranges
            if (level /= 0 .and. level /= level_butterfly + 1) then
               idx_r = idx_r0*2 - 1
               nr = nr0*2
               inc_r = inc_r0
               idx_c = ceiling_safe(idx_c0/2d0)
               if (inc_c0 > 1) then
                  nc = nc0
               else
                  nc = ceiling_safe(nc0/2d0)
               endif
               inc_c = ceiling_safe(inc_c0/2d0)
            else
               idx_r = idx_r0
               nr = nr0
               inc_r = inc_r0
               idx_c = idx_c0
               nc = nc0
               inc_c = inc_c0
            endif
         endif

         BFvec1%vec(level + 1)%idx_r = idx_r
         BFvec1%vec(level + 1)%inc_r = inc_r
         BFvec1%vec(level + 1)%nr = nr
         BFvec1%vec(level + 1)%idx_c = idx_c
         BFvec1%vec(level + 1)%inc_c = inc_c
         BFvec1%vec(level + 1)%nc = nc
         if (level /= level_butterfly + 1) then
            BFvec1%vec(level + 1)%num_row = 2**level
            BFvec1%vec(level + 1)%num_col = 2**(level_butterfly - level)
         else
            BFvec1%vec(level + 1)%num_row = 2**level_butterfly
            BFvec1%vec(level + 1)%num_col = 1
         endif
         allocate (BFvec1%vec(level + 1)%blocks(BFvec1%vec(level + 1)%nr, BFvec1%vec(level + 1)%nc))
      enddo

      !>**** compute group_ms and group_ns which stores the leaf block number (from 1 to 2^L) of each index. group_ms1 and group_ns1 are used to track the block number in each butterfly level.
      nrow = 0
      ncol = 0
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         nrow = nrow + blocks%inters(nn)%nr
         ncol = ncol + blocks%inters(nn)%nc
      enddo
      allocate (group_ms(nrow))
      allocate (group_ns(ncol))
      iidx = 0
      jidx = 0
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         do ii = 1, blocks%inters(nn)%nr
            iidx = iidx + 1
            group_ms(iidx) = findgroup(inters(nng)%rows(blocks%inters(nn)%rows(ii)), msh, level_butterfly, blocks%row_group) - blocks%row_group*2**level_butterfly + 1
         enddo
         do jj = 1, blocks%inters(nn)%nc
            jidx = jidx + 1
            group_ns(jidx) = findgroup(inters(nng)%cols(blocks%inters(nn)%cols(jj)), msh, level_butterfly, blocks%col_group) - blocks%col_group*2**level_butterfly + 1
         enddo
      enddo
      allocate (group_ms1(nrow))
      group_ms1 = group_ms
      allocate (group_ns1(ncol))
      group_ns1 = group_ns

      !>**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
      allocate (g_idx_m(0:level_half + 1))
      allocate (g_idx_n(0:level_half + 1))
      do level = 0, level_half + 1
         g_idx_m(level)%idx_r = BFvec%vec(level)%idx_r
         g_idx_m(level)%inc_r = BFvec%vec(level)%inc_r
         g_idx_m(level)%nr = BFvec%vec(level)%nr
         g_idx_m(level)%idx_c = 1
         g_idx_m(level)%inc_c = 1
         g_idx_m(level)%nc = 1
         allocate (g_idx_m(level)%blocks(g_idx_m(level)%nr, g_idx_m(level)%nc))
         allocate (g_idx_m(level)%index(iidx, 1))

         g_idx_n(level)%idx_c = BFvec%vec(level)%idx_c
         g_idx_n(level)%inc_c = BFvec%vec(level)%inc_c
         g_idx_n(level)%nc = BFvec%vec(level)%nc
         g_idx_n(level)%idx_r = 1
         g_idx_n(level)%inc_r = 1
         g_idx_n(level)%nr = 1
         allocate (g_idx_n(level)%blocks(g_idx_n(level)%nr, g_idx_n(level)%nc))
         allocate (g_idx_n(level)%index(jidx, 1))
      enddo

      iidx = 0
      jidx = 0
      allocate (num_nods_i(0:level_half + 1))
      allocate (num_nods_j(0:level_half + 1))
      do nn = 1, size(blocks%inters, 1)
         !>*** for each index in each intersection, traverse the row and column tree
         num_nods_i = 0
         num_nods_j = 0
         do level = level_butterfly + 2, 0, -1
            do ii = 1, blocks%inters(nn)%nr

               if (level <= level_half + 1) then
                  index_i_s = group_ms1(iidx + ii)
                  index_i_loc_s = (index_i_s - g_idx_m(level)%idx_r)/g_idx_m(level)%inc_r + 1
                  if (index_i_s >= g_idx_m(level)%idx_r .and. mod(index_i_s - g_idx_m(level)%idx_r, g_idx_m(level)%inc_r) == 0 .and. index_i_loc_s <= g_idx_m(level)%nr) then
                     if (g_idx_m(level)%blocks(index_i_loc_s, 1)%lst%idx /= nn) then ! not seen this index_i_loc_s before
                        num_nods_i(level) = num_nods_i(level) + 1
                        g_idx_m(level)%index(num_nods_i(level), 1) = index_i_loc_s
                        g_idx_m(level)%blocks(index_i_loc_s, 1)%lst%idx = nn
                     endif
                  endif
               endif

            enddo
            if (level > 0 .and. level < level_butterfly + 2) group_ms1(iidx + 1:iidx + blocks%inters(nn)%nr) = floor((group_ms1(iidx + 1:iidx + blocks%inters(nn)%nr) + 1)/2d0)
         enddo
         iidx = iidx + blocks%inters(nn)%nr

         do level = 0, level_butterfly + 2
            do jj = 1, blocks%inters(nn)%nc
               if (level <= level_half + 1) then
                  index_j_s = group_ns1(jidx + jj)
                  index_j_loc_s = (index_j_s - g_idx_n(level)%idx_c)/g_idx_n(level)%inc_c + 1
                  if (index_j_s >= g_idx_n(level)%idx_c .and. mod(index_j_s - g_idx_n(level)%idx_c, g_idx_n(level)%inc_c) == 0 .and. index_j_loc_s <= g_idx_n(level)%nc) then
                     if (g_idx_n(level)%blocks(1, index_j_loc_s)%lst%idx /= nn) then ! not seen this index_j_loc_s before
                        num_nods_j(level) = num_nods_j(level) + 1
                        g_idx_n(level)%index(num_nods_j(level), 1) = index_j_loc_s
                        g_idx_n(level)%blocks(1, index_j_loc_s)%lst%idx = nn
                     endif
                     if (level == 0) BFvec%vec(level)%blocks(1, index_j_loc_s)%ndim = BFvec%vec(level)%blocks(1, index_j_loc_s)%ndim + 1
                  endif
               endif

            enddo
            if (level > 0 .and. level < level_butterfly + 2) group_ns1(jidx + 1:jidx + blocks%inters(nn)%nc) = floor((group_ns1(jidx + 1:jidx + blocks%inters(nn)%nc) + 1)/2d0)

         enddo
         jidx = jidx + blocks%inters(nn)%nc

         !>**** construct BFvec%vec(level)%lst for active BF blocks and BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst for list of intersection#s
         do level = 0, level_half + 1
            do ii = 1, num_nods_i(level)
               index_i_loc_s = g_idx_m(level)%index(ii, 1)
               do jj = 1, num_nods_j(level)
                  index_j_loc_s = g_idx_n(level)%index(jj, 1)
                  if (BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods == 0) then ! first see this block
                     p%i = index_i_loc_s
                     p%j = index_j_loc_s
                     ! write(*,*)p%i,p%j,BFvec%vec(level)%lst%num_nods
                     call append(BFvec%vec(level)%lst, p)
                  endif
                  call append(BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst, nn)
               enddo
            enddo
         enddo

      enddo

      do level = 0, level_half + 1
         deallocate (g_idx_m(level)%index)
         deallocate (g_idx_n(level)%index)
      enddo
      deallocate (g_idx_m)
      deallocate (g_idx_n)
      deallocate (group_ms1)
      deallocate (group_ns1)
      deallocate (num_nods_i)
      deallocate (num_nods_j)

      !>**** copy *%lst to *%index
      do level = 0, level_half + 1
         allocate (BFvec%vec(level)%index(BFvec%vec(level)%lst%num_nods, 2))
         cur => BFvec%vec(level)%lst%head
         do nn = 1, BFvec%vec(level)%lst%num_nods
            select type (ptr=>cur%item)
            type is (ipair)
               BFvec%vec(level)%index(nn, 1) = ptr%i
               BFvec%vec(level)%index(nn, 2) = ptr%j
            end select
            cur => cur%next
         enddo
         call list_finalizer(BFvec%vec(level)%lst)

         do nn = 1, size(BFvec%vec(level)%index, 1)
            index_i_loc_s = BFvec%vec(level)%index(nn, 1)
            index_j_loc_s = BFvec%vec(level)%index(nn, 2)
            ! if(ptree%MyID>=2)write(*,*)ptree%MyID,'my',blocks%row_group,blocks%col_group,'BFvec range',level,index_i_loc_s,index_j_loc_s,(index_i_loc_s-1)*BFvec%vec(level)%inc_r+BFvec%vec(level)%idx_r,(index_j_loc_s-1)*BFvec%vec(level)%inc_c+BFvec%vec(level)%idx_c
            allocate (BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods, 1))
            cur => BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%head
            do ii = 1, BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods
               select type (ptr=>cur%item)
               type is (integer)
                  BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 1) = ptr
               end select
               cur => cur%next
            enddo
            call list_finalizer(BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst)
         enddo
      enddo

      allocate (group_ms1(nrow))
      group_ms1 = group_ms
      allocate (group_ns1(ncol))
      group_ns1 = group_ns

      !>**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
      allocate (g_idx_m(level_half + 1:level_butterfly + 2))
      allocate (g_idx_n(level_half + 1:level_butterfly + 2))
      do level = level_half + 1, level_butterfly + 2
         g_idx_m(level)%idx_r = BFvec1%vec(level)%idx_r
         g_idx_m(level)%inc_r = BFvec1%vec(level)%inc_r
         g_idx_m(level)%nr = BFvec1%vec(level)%nr
         g_idx_m(level)%idx_c = 1
         g_idx_m(level)%inc_c = 1
         g_idx_m(level)%nc = 1
         allocate (g_idx_m(level)%blocks(g_idx_m(level)%nr, g_idx_m(level)%nc))
         allocate (g_idx_m(level)%index(iidx, 1))

         g_idx_n(level)%idx_c = BFvec1%vec(level)%idx_c
         g_idx_n(level)%inc_c = BFvec1%vec(level)%inc_c
         g_idx_n(level)%nc = BFvec1%vec(level)%nc
         g_idx_n(level)%idx_r = 1
         g_idx_n(level)%inc_r = 1
         g_idx_n(level)%nr = 1
         allocate (g_idx_n(level)%blocks(g_idx_n(level)%nr, g_idx_n(level)%nc))
         allocate (g_idx_n(level)%index(jidx, 1))
      enddo

      !>**** compute g_idx_m and g_idx_n which stores the local (with halo) blocks for the block rows and columns
      iidx = 0
      jidx = 0
      allocate (num_nods_i(level_half + 1:level_butterfly + 2))
      allocate (num_nods_j(level_half + 1:level_butterfly + 2))
      do nn = 1, size(blocks%inters, 1)
         num_nods_i = 0
         num_nods_j = 0
         !>*** for each index in each intersection, traverse the row and column tree
         do level = level_butterfly + 2, 0, -1
            do ii = 1, blocks%inters(nn)%nr

               if (level >= level_half + 1) then
                  index_i_s = group_ms1(iidx + ii)
                  index_i_loc_s = (index_i_s - g_idx_m(level)%idx_r)/g_idx_m(level)%inc_r + 1
                  if (index_i_s >= g_idx_m(level)%idx_r .and. mod(index_i_s - g_idx_m(level)%idx_r, g_idx_m(level)%inc_r) == 0 .and. index_i_loc_s <= g_idx_m(level)%nr) then
                     if (g_idx_m(level)%blocks(index_i_loc_s, 1)%lst%idx /= nn) then ! not seen this index_i_loc_s before
                        num_nods_i(level) = num_nods_i(level) + 1
                        g_idx_m(level)%index(num_nods_i(level), 1) = index_i_loc_s
                        g_idx_m(level)%blocks(index_i_loc_s, 1)%lst%idx = nn
                     endif
                  endif
               endif

            enddo
            if (level > 0 .and. level < level_butterfly + 2) group_ms1(iidx + 1:iidx + blocks%inters(nn)%nr) = floor((group_ms1(iidx + 1:iidx + blocks%inters(nn)%nr) + 1)/2d0)
         enddo
         iidx = iidx + blocks%inters(nn)%nr

         do level = 0, level_butterfly + 2
            do jj = 1, blocks%inters(nn)%nc
               if (level >= level_half + 1) then
                  index_j_s = group_ns1(jidx + jj)
                  index_j_loc_s = (index_j_s - g_idx_n(level)%idx_c)/g_idx_n(level)%inc_c + 1
                  if (index_j_s >= g_idx_n(level)%idx_c .and. mod(index_j_s - g_idx_n(level)%idx_c, g_idx_n(level)%inc_c) == 0 .and. index_j_loc_s <= g_idx_n(level)%nc) then
                     if (g_idx_n(level)%blocks(1, index_j_loc_s)%lst%idx /= nn) then ! not seen this index_j_loc_s before
                        num_nods_j(level) = num_nods_j(level) + 1
                        g_idx_n(level)%index(num_nods_j(level), 1) = index_j_loc_s
                        g_idx_n(level)%blocks(1, index_j_loc_s)%lst%idx = nn

                     endif
                     if (level == 0) BFvec1%vec(level)%blocks(1, index_j_loc_s)%ndim = BFvec1%vec(level)%blocks(1, index_j_loc_s)%ndim + 1
                  endif
               endif

            enddo
            if (level > 0 .and. level < level_butterfly + 2) group_ns1(jidx + 1:jidx + blocks%inters(nn)%nc) = floor((group_ns1(jidx + 1:jidx + blocks%inters(nn)%nc) + 1)/2d0)
         enddo
         jidx = jidx + blocks%inters(nn)%nc

         !>**** construct BFvec1%vec(level)%lst for active BF blocks and BFvec1%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%lst for list of intersection#s
         do level = level_half + 1, level_butterfly + 2
            do ii = 1, num_nods_i(level)
               index_i_loc_s = g_idx_m(level)%index(ii, 1)
               do jj = 1, num_nods_j(level)
                  index_j_loc_s = g_idx_n(level)%index(jj, 1)
                  if (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods == 0) then ! first see this block
                     p%i = index_i_loc_s
                     p%j = index_j_loc_s
                     call append(BFvec1%vec(level)%lst, p)
                  endif
                  call append(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst, nn)
               enddo
            enddo
         enddo
      enddo
      do level = level_half + 1, level_butterfly + 2
         deallocate (g_idx_m(level)%index)
         deallocate (g_idx_n(level)%index)
      enddo
      deallocate (g_idx_m)
      deallocate (g_idx_n)
      deallocate (group_ms1)
      deallocate (group_ns1)
      deallocate (num_nods_i)
      deallocate (num_nods_j)

      !>**** copy *%lst to *%index
      do level = level_half + 1, level_butterfly + 2
         allocate (BFvec1%vec(level)%index(BFvec1%vec(level)%lst%num_nods, 2))
         cur => BFvec1%vec(level)%lst%head
         do nn = 1, BFvec1%vec(level)%lst%num_nods
            select type (ptr=>cur%item)
            type is (ipair)
               BFvec1%vec(level)%index(nn, 1) = ptr%i
               BFvec1%vec(level)%index(nn, 2) = ptr%j
            end select
            cur => cur%next
         enddo
         call list_finalizer(BFvec1%vec(level)%lst)

         do nn = 1, size(BFvec1%vec(level)%index, 1)
            index_i_loc_s = BFvec1%vec(level)%index(nn, 1)
            index_j_loc_s = BFvec1%vec(level)%index(nn, 2)

! if(ptree%MyID>=2)write(*,*)ptree%MyID,'my',blocks%row_group,blocks%col_group,'BFvec1 range',level,index_i_loc_s,index_j_loc_s,(index_i_loc_s-1)*BFvec1%vec(level)%inc_r+BFvec1%vec(level)%idx_r,(index_j_loc_s-1)*BFvec1%vec(level)%inc_c+BFvec1%vec(level)%idx_c

            if (level /= level_butterfly + 2) then ! last level doesn't need a list of intersection#
               allocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods, 1))
               cur => BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%head
               do ii = 1, BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods
                  select type (ptr=>cur%item)
                  type is (integer)
                     BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 1) = ptr
                  end select
                  cur => cur%next
               enddo
            endif
            call list_finalizer(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst)
         enddo
      enddo

      !>**** create a list of row_loc indices for each block of BFvec1%vec(level_butterfly+2),note that BFvec1%vec(level)%blocks(index_i_loc_s,1)%lst and BFvec1%vec(level)%blocks(index_i_loc_s,1)%index are used at this level compared to other levels
      iidx = 0
      do nn = 1, size(blocks%inters, 1)
         !>*** for each index in each intersection, traverse the row and column tree
         do ii = 1, blocks%inters(nn)%nr
            iidx = iidx + 1
            level = level_butterfly + 2
            index_i_s = group_ms(iidx)
            index_i_loc_s = (index_i_s - BFvec1%vec(level)%idx_r)/BFvec1%vec(level)%inc_r + 1
            if (index_i_s >= BFvec1%vec(level)%idx_r .and. mod(index_i_s - BFvec1%vec(level)%idx_r, BFvec1%vec(level)%inc_r) == 0 .and. index_i_loc_s <= BFvec1%vec(level)%nr) then
               p%i = nn
               p%j = blocks%inters(nn)%glo2loc(ii)
               call append(BFvec1%vec(level)%blocks(index_i_loc_s, 1)%lst, p)
               ! write(*,*)ptree%MyID,'idtarg',index_i_s,index_i_loc_s
            endif
         enddo
      enddo

      !>**** copy *%lst to *%index
      level = level_butterfly + 2
      do nn = 1, size(BFvec1%vec(level)%index, 1)
         index_i_loc_s = BFvec1%vec(level)%index(nn, 1)
         index_j_loc_s = BFvec1%vec(level)%index(nn, 2)

         allocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods, 2))
         cur => BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%head
         do ii = 1, BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst%num_nods
            select type (ptr=>cur%item)
            type is (ipair)
               BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 1) = ptr%i
               BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 2) = ptr%j
            end select
            cur => cur%next
         enddo
         call list_finalizer(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%lst)
      enddo

      n2 = MPI_Wtime()

      !>**** generate data in BFvec%vec(0)
      do nn = 1, size(BFvec%vec(0)%index, 1)
         index_j_loc_s = BFvec%vec(0)%index(nn, 2)
         allocate (BFvec%vec(0)%blocks(1, index_j_loc_s)%matrix(3, BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim))
         BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim = 0
      enddo
      jidx = 0
      do nn = 1, size(blocks%inters, 1)
         nng = blocks%inters(nn)%idx
         do jj = 1, blocks%inters(nn)%nc
            jidx = jidx + 1
            index_j_s = group_ns(jidx)
            index_j_loc_s = (index_j_s - BFvec%vec(0)%idx_c)/BFvec%vec(0)%inc_c + 1
            if (index_j_s >= BFvec%vec(0)%idx_c .and. mod(index_j_s - BFvec%vec(0)%idx_c, BFvec%vec(0)%inc_c) == 0 .and. index_j_loc_s <= BFvec%vec(0)%nc) then
               BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim = BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim + 1
               BFvec%vec(0)%blocks(1, index_j_loc_s)%matrix(1, BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim) = nn
               BFvec%vec(0)%blocks(1, index_j_loc_s)%matrix(2, BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim) = jj
               BFvec%vec(0)%blocks(1, index_j_loc_s)%matrix(3, BFvec%vec(0)%blocks(1, index_j_loc_s)%ndim) = inters(nng)%cols(blocks%inters(nn)%cols(jj)) - msh%basis_group(blocks%col_group*2**level_butterfly + index_j_s - 1)%head + 1
            endif
         enddo
      enddo

      !>**** multiply BF with BFvec
      do level = 0, level_half
#ifdef HAVE_TASKLOOP
         !$omp taskloop default(shared) private(nn)
#endif
         do nn = 1, size(BFvec%vec(level + 1)%index, 1)
            call BF_block_extraction_multiply_oneblock_right(blocks, BFvec, level,nn,ptree,stats)
         enddo
#ifdef HAVE_TASKLOOP
         !$omp end taskloop
#endif

         do nn = 1, size(BFvec%vec(level)%index, 1)
            index_i_loc_s = BFvec%vec(level)%index(nn, 1)
            index_j_loc_s = BFvec%vec(level)%index(nn, 2)
            if (allocated(BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)) deallocate (BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)
            if (associated(BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) deallocate (BFvec%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)
         enddo
         if (level_half /= level) then
            call BF_exchange_extraction(blocks, BFvec%vec(level + 1), stats, ptree, level, 'B')
         endif
      enddo

      n3 = MPI_Wtime()

      ! write(*,*)blocks%row_group,blocks%col_group,'myid',ptree%MyID,'level',level,'before all2all'
      ! all2all communication from BFvec%vec(level_half+1) to BFvec1%vec(level_half+1)
      call BF_all2all_extraction(blocks, BFvec%vec(level_half + 1), BFvec1%vec(level_half + 1), stats, ptree, level_half, 'R', 'C')
      ! write(*,*)blocks%row_group,blocks%col_group,'myid',ptree%MyID,'level',level,'after all2all'

      n4 = MPI_Wtime()

      !>**** multiply BF with BFvec1
      do level = level_half + 1, level_butterfly + 1
         if (level == 0) then
            write (*, *) 'should not come here as level_half>=0'
            stop
         elseif (level == level_butterfly + 1) then
#ifdef HAVE_TASKLOOP
            !$omp taskloop default(shared) private(nn)
#endif
            do nn = 1, size(BFvec1%vec(level + 1)%index, 1)
               call BF_block_extraction_multiply_oneblock_last(blocks, BFvec1, inters, level,nn,ptree, msh,stats)
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
#endif

            do nn = 1, size(BFvec1%vec(level)%index, 1)
               index_i_loc_s = BFvec1%vec(level)%index(nn, 1)
               index_j_loc_s = BFvec1%vec(level)%index(nn, 2)
               if (allocated(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)) deallocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)
               if (associated(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) deallocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)
            enddo

            do nn = 1, size(BFvec1%vec(level + 1)%index, 1)
               index_i_loc_s = BFvec1%vec(level + 1)%index(nn, 1)
               index_j_loc_s = BFvec1%vec(level + 1)%index(nn, 2)
               if (allocated(BFvec1%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index)) deallocate (BFvec1%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index)
               if (associated(BFvec1%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) deallocate (BFvec1%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)
            enddo

         else
#ifdef HAVE_TASKLOOP
            !$omp taskloop default(shared) private(nn)
#endif
            do nn = 1, size(BFvec1%vec(level + 1)%index, 1)
               call BF_block_extraction_multiply_oneblock_left(blocks, BFvec1, level,nn,ptree,stats)
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
#endif

            do nn = 1, size(BFvec1%vec(level)%index, 1)
               index_i_loc_s = BFvec1%vec(level)%index(nn, 1)
               index_j_loc_s = BFvec1%vec(level)%index(nn, 2)
               if (allocated(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)) deallocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%index)
               if (associated(BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) deallocate (BFvec1%vec(level)%blocks(index_i_loc_s, index_j_loc_s)%matrix)
            enddo

            call BF_exchange_extraction(blocks, BFvec1%vec(level + 1), stats, ptree, level, 'R')

            !!!!!! sort the intersection# here
#ifdef HAVE_TASKLOOP
            !$omp taskloop default(shared) private(nn)
#endif
            do nn = 1, size(BFvec1%vec(level + 1)%index, 1)
               call BF_block_extraction_sort_oneblock(blocks, BFvec1, level,nn,ptree)
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
#endif
         endif
      enddo

      deallocate (group_ms)
      deallocate (group_ns)

      do level = level_half + 1, level_butterfly + 2
         if (allocated(BFvec1%vec(level)%index)) deallocate (BFvec1%vec(level)%index)
      enddo
      deallocate (BFvec1%vec)

      do level = 0, level_half + 1
         if (allocated(BFvec%vec(level)%index)) deallocate (BFvec%vec(level)%index)
      enddo
      deallocate (BFvec%vec)

      n5 = MPI_Wtime()
      ! time_tmp = time_tmp + n5 - n1
      ! stats%Time_Entry_BF = stats%Time_Entry_BF + n5-n2


   end subroutine BF_block_extraction



subroutine BF_block_extraction_multiply_oneblock_right(blocks, BFvec, level,nn,ptree,stats)

    implicit none
    type(matrixblock)::blocks
    type(proctree)::ptree
    type(Hstat):: stats
    type(butterfly_vec) :: BFvec
    integer level,nn

    integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
    integer level_butterfly, num_blocks, level_half
    integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
    integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
    integer pid, pgno_sub
    DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)


    level_butterfly = blocks%level_butterfly
    index_i_loc_s = BFvec%vec(level + 1)%index(nn, 1)
    index_i_s = (index_i_loc_s - 1)*BFvec%vec(level + 1)%inc_r + BFvec%vec(level + 1)%idx_r
    index_j_loc_s = BFvec%vec(level + 1)%index(nn, 2)
    index_j_s = (index_j_loc_s - 1)*BFvec%vec(level + 1)%inc_c + BFvec%vec(level + 1)%idx_c
    call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i_s, index_j_s, 'R', pgno_sub)
    pid = ptree%pgrp(pgno_sub)%head
    call assert(ptree%pgrp(pgno_sub)%nproc==1,'pgno_sub can only has one proc here')
    if (pid == ptree%MyID) then
        if (level == 0) then
            index_j = index_j_s
            index_j_loc_k = (index_j - blocks%ButterflyV%idx)/blocks%ButterflyV%inc + 1
            index_jj_loc = index_j_loc_k
            rank = size(blocks%ButterflyV%blocks(index_j_loc_k)%matrix, 2)
            allocate (BFvec%vec(level + 1)%blocks(1, index_j_loc_s)%matrix(rank + 2, BFvec%vec(level)%blocks(1, index_jj_loc)%ndim))
            BFvec%vec(level + 1)%blocks(1, index_j_loc_s)%matrix(1:2, :) = BFvec%vec(level)%blocks(1, index_jj_loc)%matrix(1:2, :)
            do nnn = 1, BFvec%vec(level)%blocks(1, index_jj_loc)%ndim
                jjj = NINT(dble(BFvec%vec(level)%blocks(1, index_jj_loc)%matrix(3, nnn)))
                BFvec%vec(level + 1)%blocks(1, index_j_loc_s)%matrix(3:rank + 2, nnn) = blocks%ButterflyV%blocks(index_j_loc_k)%matrix(jjj, :)
            enddo
        else if (level == level_butterfly + 1) then
            write (*, *) 'should not come here as level_half<=level_butterfly'
            stop
        else
            index_i = index_i_s
            index_j = index_j_s
            index_ii = floor_safe((index_i + 1)/2d0); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)

            index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
            index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

            index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
            index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

            nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
            nn2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
            mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
            if (associated(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix)) then
                nvec1 = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
            else
                nvec1 = 0
            endif

            if (associated(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix)) then
                nvec2 = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix, 2)
            else
                nvec2 = 0
            endif

            allocate (mat1(mm + 2, nvec1))
            if (nvec1 > 0) then
                mat1 = 0
                mat1(1:2, :) = BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1:2, :)
                call gemmf77('N', 'N', mm, nvec1, nn1, BPACK_cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(3, 1), nn1 + 2, BPACK_czero, mat1(3, 1), mm + 2)
#ifdef HAVE_TASKLOOP
                !$omp atomic
#endif
                stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm, nvec1, nn1)
#ifdef HAVE_TASKLOOP
                !$omp end atomic
#endif
            endif

            allocate (mat2(mm + 2, nvec2))
            if (nvec2 > 0) then
                mat2 = 0
                mat2(1:2, :) = BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix(1:2, :)
                call gemmf77('N', 'N', mm, nvec2, nn2, BPACK_cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc + 1)%matrix(3, 1), nn2 + 2, BPACK_czero, mat2(3, 1), mm + 2)
#ifdef HAVE_TASKLOOP
                !$omp atomic
#endif
                stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm, nvec2, nn2)
#ifdef HAVE_TASKLOOP
                !$omp end atomic
#endif
            endif

            !>**** filter out the columns of mat1 and mat2 that are not specified by BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index
            allocate (mat(mm + 2, nvec1 + nvec2))
            idx1 = 0
            idx2 = 0
            idx = 0
            gg1 = 0
            gg2 = 0
            do ii = 1, size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index, 1)
                gg = NINT(dble(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 1)))
                do while (gg1 <= gg)
                if (gg1 == gg) then
                    idx = idx + 1
                    mat(:, idx) = mat1(:, idx1)
                endif
                idx1 = idx1 + 1
                if (idx1 > nvec1) exit
                gg1 = NINT(dble(mat1(1, idx1)))
                enddo
                do while (gg2 <= gg)
                if (gg2 == gg) then
                    idx = idx + 1
                    mat(:, idx) = mat2(:, idx2)
                endif
                idx2 = idx2 + 1
                if (idx2 > nvec2) exit
                gg2 = NINT(dble(mat2(1, idx2)))
                enddo
            enddo
            allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm + 2, idx))
            BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat(:, 1:idx)
            deallocate (mat1)
            deallocate (mat2)
            deallocate (mat)

        endif
    endif

end subroutine BF_block_extraction_multiply_oneblock_right




subroutine BF_MD_block_extraction_multiply_oneblock_right(blocks, bb_m, Ndim, BFvec, idx_r_m, level ,dim_i, ptree,stats)

    implicit none
    type(matrixblock_MD)::blocks
    type(proctree)::ptree
    type(Hstat):: stats
    type(butterfly_vec) :: BFvec
    integer bb_m, level,nn, Ndim, idx_r_m(Ndim),idx_r_m_1(Ndim), dims(Ndim)

    integer ii, iii, jj, jjj, ll, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
    integer level_butterfly, num_blocks, level_half
    integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
    integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
    integer pid, pgno_sub, index_r_scalar, dim_i
    DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)


    level_butterfly = blocks%level_butterfly
    level_half = blocks%level_half

    idx_r_m_1=idx_r_m
    dims = 2**level_half
    do ll= level_half-1,level,-1
      idx_r_m_1=int((idx_r_m_1 + 1)/2)
      dims=dims/2
    enddo
    call MultiIndexToSingleIndex(Ndim,dims, index_r_scalar, idx_r_m_1)

    do index_j_s=1,BFvec%vec(level+1)%nc
      if (level == 0) then
         index_j_k = index_j_s
         rank = size(blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_i)%matrix, 2)
         nn = size(blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_i)%matrix, 1)
         allocate (BFvec%vec(level + 1)%blocks(1, index_j_s)%matrix(rank, nn))
         do nnn = 1, nn
            BFvec%vec(level + 1)%blocks(1, index_j_s)%matrix(:, nnn) = blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_i)%matrix(nnn, :)
         enddo
      else
         index_j_k = 2*index_j_s - 1 !index_ii is global index in BFvec%vec(level)

            nn1 = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_r_scalar, index_j_k,dim_i)%matrix, 2)
            nn2 = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_r_scalar, index_j_k+1,dim_i)%matrix, 2)
            mm = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_r_scalar, index_j_k,dim_i)%matrix, 1)
            nvec1 = size(BFvec%vec(level)%blocks(1, index_j_k)%matrix, 2)
            nvec2 = size(BFvec%vec(level)%blocks(1, index_j_k + 1)%matrix, 2)


            allocate (mat1(mm, nvec1))
            if (nvec1 > 0) then
                mat1 = 0
                call gemmf77('N', 'N', mm, nvec1, nn1, BPACK_cone, blocks%ButterflyKerl_R(bb_m,level)%blocks(index_r_scalar, index_j_k,dim_i)%matrix, mm, BFvec%vec(level)%blocks(1, index_j_k)%matrix, nn1, BPACK_czero, mat1, mm)
#ifdef HAVE_TASKLOOP
                !$omp atomic
#endif
                stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm, nvec1, nn1)
#ifdef HAVE_TASKLOOP
                !$omp end atomic
#endif
            endif

            allocate (mat2(mm, nvec2))
            if (nvec2 > 0) then
                mat2 = 0
                call gemmf77('N', 'N', mm, nvec2, nn2, BPACK_cone, blocks%ButterflyKerl_R(bb_m,level)%blocks(index_r_scalar, index_j_k+1,dim_i)%matrix, mm, BFvec%vec(level)%blocks(1, index_j_k + 1)%matrix, nn2, BPACK_czero, mat2, mm)
#ifdef HAVE_TASKLOOP
                !$omp atomic
#endif
                stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm, nvec2, nn2)
#ifdef HAVE_TASKLOOP
                !$omp end atomic
#endif
            endif
            allocate (BFvec%vec(level + 1)%blocks(1, index_j_s)%matrix(mm, nvec1+nvec2))
            BFvec%vec(level + 1)%blocks(1, index_j_s)%matrix(:,1:nvec1) = mat1
            BFvec%vec(level + 1)%blocks(1, index_j_s)%matrix(:,1+nvec1:nvec1+nvec2) = mat2
            deallocate (mat1)
            deallocate (mat2)
      endif
    enddo

end subroutine BF_MD_block_extraction_multiply_oneblock_right




subroutine BF_block_extraction_multiply_oneblock_left(blocks, BFvec, level,nn,ptree,stats)

   implicit none
   type(matrixblock)::blocks
   type(proctree)::ptree
   type(Hstat):: stats
   type(butterfly_vec) :: BFvec
   integer level,nn

   integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
   integer level_butterfly, num_blocks, level_half
   integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
   integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
   integer pid, pgno_sub
   DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)

   level_butterfly = blocks%level_butterfly
   index_i_loc_s = BFvec%vec(level + 1)%index(nn, 1)
   index_i_s = (index_i_loc_s - 1)*BFvec%vec(level + 1)%inc_r + BFvec%vec(level + 1)%idx_r
   index_j_loc_s = BFvec%vec(level + 1)%index(nn, 2)
   index_j_s = (index_j_loc_s - 1)*BFvec%vec(level + 1)%inc_c + BFvec%vec(level + 1)%idx_c
   do jjj = 0, 1
      index_i = floor_safe((index_i_s - 1)/2d0) + 1
      index_j = 2*index_j_s - jjj
      call GetBlockPID(ptree, blocks%pgno, level, level_butterfly, index_i, index_j, 'C', pgno_sub)
      call assert(ptree%pgrp(pgno_sub)%nproc==1,'pgno_sub can only has one proc here')
      pid = ptree%pgrp(pgno_sub)%head
      if (pid == ptree%MyID) then
         index_i_k = 2*index_i - mod(index_i_s, 2)
         index_j_k = index_j

         index_i_loc_k = (index_i_k - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
         index_j_loc_k = (index_j_k - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

         index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
         index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

         if (associated(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix)) then
            mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
            nn2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
            nvec2 = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)

            allocate (mat2(mm + 2, nvec2))
            mat2 = 0
            mat2(1:2, :) = BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1:2, :)
            call gemmf77('N', 'N', mm, nvec2, nn2, BPACK_cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(3, 1), nn2 + 2, BPACK_czero, mat2(3, 1), mm + 2)

#ifdef HAVE_TASKLOOP
            !$omp atomic
#endif
            stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(mm, nvec2, nn2)
#ifdef HAVE_TASKLOOP
            !$omp end atomic
#endif

            if (associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then
               nvec1 = size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, 2)
               allocate (mat1(mm + 2, nvec1))
               mat1 = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix
               deallocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)
            else
               nvec1 = 0
               allocate (mat1(mm + 2, nvec1))
            endif

            allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm + 2, nvec1 + nvec2))
            if (nvec1 > 0) BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(:, 1:nvec1) = mat1
            if (nvec2 > 0) BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(:, 1 + nvec1:nvec1 + nvec2) = mat2

            deallocate (mat1)
            deallocate (mat2)

         endif

      endif
   enddo

end subroutine BF_block_extraction_multiply_oneblock_left







subroutine BF_MD_block_mvp_multiply_right(blocks, bb_m, Ndim, BFvec, Nvec, level, ptree,stats)

    implicit none
    type(matrixblock_MD)::blocks
    type(proctree)::ptree
    type(Hstat):: stats
    type(butterfly_vec) :: BFvec
    integer bb_m, level,nn, Ndim, idx_c_m(Ndim),idx_c_m_1(Ndim), dims(Ndim),dims_MD_old(Ndim*2),dims_MD_new(Ndim*2),dims_ref(Ndim+1),idx_ref(Ndim+1),offsets_ref(Ndim+1),idx_ref_scalar, dims_in(Ndim),idx_in(Ndim), dims_MD3(Ndim),idx_MD3(Ndim), idx_MD(Ndim), dims_ref_scalar,dims_MD3_scalar

    integer i, j, ii, iii, jj, jjj, ll, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr(Ndim), idxr(Ndim), nr0(Ndim), idxr0(Ndim), inc_c, inc_c0, inc_r, inc_r0
    integer level_butterfly, num_blocks, level_half, levelm
    integer group_m(Ndim), group_n(Ndim), mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_ii_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, Nvec
    integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
    integer pid, pgno_sub, index_c_scalar, dim_i, dim_ii
    DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)
    real(kind=8)::flop,flops


    level_butterfly = blocks%level_butterfly
    level_half = blocks%level_half
    levelm = level_half
    offsets_ref=0

    call SingleIndexToMultiIndex(Ndim, blocks%nc_m, bb_m, idx_c_m)
    group_n = blocks%col_group
    group_n = group_n*2**(level_butterfly-levelm) - 1 + idx_c_m + blocks%idx_c_m - 1

   if (level == 0) then
      do dim_ii=1,Ndim
         if(dim_ii==1)then
            do dim_i=1,Ndim
               dims_ref(dim_i) = BFvec%vec(0)%index_MD(dim_i,1,BFvec%vec(0)%num_col+1)
            enddo
            dims_ref(Ndim+1)=Nvec
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,dim_ii,dim_ii,BFvec%vec(level)%blocks(1,1)%matrix, 'N', mat, 'N',1)
         else
            do dim_i=1,Ndim
               if(dim_i<dim_ii)then
                  dims_ref(dim_i) = BFvec%vec(1)%index_MD(dim_i,1,BFvec%vec(1)%num_col+1)
               else
                  dims_ref(dim_i) = BFvec%vec(0)%index_MD(dim_i,1,BFvec%vec(0)%num_col+1)
               endif
            enddo
            dims_ref(Ndim+1)=Nvec
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,dim_ii-1,dim_ii,BFvec%vec(level+1)%blocks(1,1)%matrix, 'N', mat, 'N',1)
         endif

         if(associated(BFvec%vec(level+1)%blocks(1,1)%matrix))deallocate(BFvec%vec(level+1)%blocks(1,1)%matrix)
         allocate(BFvec%vec(level+1)%blocks(1,1)%matrix(BFvec%vec(level+1)%index_MD(dim_ii,1,BFvec%vec(level+1)%num_col+1),size(mat,2)))
         BFvec%vec(level+1)%blocks(1,1)%matrix=0

         flops=0
#ifdef HAVE_TASKLOOP
         ! !$omp taskloop default(shared) private(index_j_s,index_j_k,mm,nn1,nvec1) reduction(+:flops)
         !$omp taskgroup task_reduction(+:flops)
#else
#ifdef HAVE_OPENMP
         !$omp parallel do default(shared) private(index_j_s,index_j_k,mm,nn1,nvec1) reduction(+:flops)
#endif
#endif
         do index_j_s=1,BFvec%vec(level+1)%nc
            index_j_k = index_j_s
            mm = size(blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_ii)%matrix,2)
            nn1 = size(blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_ii)%matrix,1)
            nvec1 = size(mat,2)
#ifdef HAVE_TASKLOOP
            !$omp task default(shared) firstprivate(mm,nn1,nvec1,index_j_k)
#endif
            call gemmf77('T', 'N', mm, nvec1, nn1, BPACK_cone, blocks%ButterflyV(bb_m)%blocks(index_j_k,dim_ii)%matrix, nn1, mat(BFvec%vec(level)%index_MD(dim_ii,1,index_j_s)+1,1), size(mat,1), BPACK_czero, BFvec%vec(level+1)%blocks(1,1)%matrix(BFvec%vec(level+1)%index_MD(dim_ii,1,index_j_s)+1,1), BFvec%vec(level+1)%index_MD(dim_ii,1,BFvec%vec(level+1)%num_col+1))
            flops = flops + flops_gemm(mm, nvec1, nn1)
#ifdef HAVE_TASKLOOP
            !$omp end task
#endif
         enddo
#ifdef HAVE_TASKLOOP
         !   !$omp end taskloop
         !$omp end taskgroup
#else
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
#endif
#ifdef HAVE_TASKLOOP
         !$omp atomic
#endif
         stats%Flop_Tmp = stats%Flop_Tmp + flops
#ifdef HAVE_TASKLOOP
         !$omp end atomic
#endif
         deallocate(mat)
      enddo
   else
      nr=2**(level)
      nr0=2**(level-1)

      do index_i_s=1,product(nr)
         call SingleIndexToMultiIndex(Ndim, nr, index_i_s, idxr)
         idxr0=int((idxr + 1)/2)
         call MultiIndexToSingleIndex(Ndim, nr0, index_ii_s, idxr0)
         do dim_ii=1,Ndim
            if(dim_ii==1)then
               do dim_i=1,Ndim
                  dims_ref(dim_i) = BFvec%vec(level)%index_MD(dim_i,index_ii_s,BFvec%vec(level)%num_col+1)
               enddo
               dims_ref(Ndim+1)=Nvec
               call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,Ndim,dim_ii,BFvec%vec(level)%blocks(index_ii_s,1)%matrix, 'N', mat, 'N',1)
            else
               do dim_i=1,Ndim
                  if(dim_i<dim_ii)then
                     dims_ref(dim_i) = BFvec%vec(level+1)%index_MD(dim_i,index_i_s,BFvec%vec(level+1)%num_col+1)
                  else
                     dims_ref(dim_i) = BFvec%vec(level)%index_MD(dim_i,index_ii_s,BFvec%vec(level)%num_col+1)
                  endif
               enddo
               dims_ref(Ndim+1)=Nvec
               call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,dim_ii-1,dim_ii,BFvec%vec(level+1)%blocks(index_i_s,1)%matrix, 'N', mat, 'N',1)
            endif

            if(associated(BFvec%vec(level+1)%blocks(index_i_s,1)%matrix))deallocate(BFvec%vec(level+1)%blocks(index_i_s,1)%matrix)
            allocate(BFvec%vec(level+1)%blocks(index_i_s,1)%matrix(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,BFvec%vec(level+1)%num_col+1),size(mat,2)))
            BFvec%vec(level+1)%blocks(index_i_s,1)%matrix=0

            flops=0
#ifdef HAVE_TASKLOOP
            !   !$omp taskloop default(shared) private(index_j_s,index_j_k,mm,nn1,nn2,nvec1) reduction(+:flops)
            !$omp taskgroup task_reduction(+:flops)
#else
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(index_j_s,index_j_k,mm,nn1,nn2,nvec1) reduction(+:flops)
#endif
#endif
            do index_j_s=1,BFvec%vec(level+1)%nc
               index_j_k = 2*index_j_s - 1
               nn1 = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_s, index_j_k,dim_ii)%matrix, 2)
               nn2 = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_s, index_j_k+1,dim_ii)%matrix, 2)
               mm = size(blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_s, index_j_k,dim_ii)%matrix, 1)
               nvec1 = size(mat,2)
#ifdef HAVE_TASKLOOP
               !$omp task default(shared) firstprivate(mm,nn1,nn2,nvec1,index_j_k)
#endif
               call gemmf77('N', 'N', mm, nvec1, nn1, BPACK_cone, blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_s, index_j_k,dim_ii)%matrix, mm, mat(BFvec%vec(level)%index_MD(dim_ii,index_ii_s,index_j_k)+1,1), size(mat,1), BPACK_cone, BFvec%vec(level+1)%blocks(index_i_s,1)%matrix(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,index_j_s)+1,1), BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,BFvec%vec(level+1)%num_col+1))
               flops = flops + flops_gemm(mm, nvec1, nn1)
               call gemmf77('N', 'N', mm, nvec1, nn2, BPACK_cone, blocks%ButterflyKerl_R(bb_m,level)%blocks(index_i_s, index_j_k+1,dim_ii)%matrix, mm, mat(BFvec%vec(level)%index_MD(dim_ii,index_ii_s,index_j_k+1)+1,1), size(mat,1), BPACK_cone, BFvec%vec(level+1)%blocks(index_i_s,1)%matrix(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,index_j_s)+1,1), BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,BFvec%vec(level+1)%num_col+1))
               flops = flops + flops_gemm(mm, nvec1, nn2)
#ifdef HAVE_TASKLOOP
               !$omp end task
#endif
            enddo

#ifdef HAVE_TASKLOOP
            !  !$omp end taskloop
            !$omp end taskgroup
#else
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
#endif
#ifdef HAVE_TASKLOOP
            !$omp atomic
#endif
            stats%Flop_Tmp = stats%Flop_Tmp + flops
#ifdef HAVE_TASKLOOP
            !$omp end atomic
#endif
            deallocate(mat)
         enddo
      enddo
   endif

end subroutine BF_MD_block_mvp_multiply_right



subroutine BF_MD_block_mvp_multiply_left(blocks, bb_m, Ndim, BFvec, Nvec, level, ptree,stats)

    implicit none
    type(matrixblock_MD)::blocks
    type(proctree)::ptree
    type(Hstat):: stats
    type(butterfly_vec) :: BFvec
    integer bb_m, level,nn, Ndim, idx_r_m(Ndim),idx_r_m_1(Ndim), dims(Ndim),dims_MD_old(Ndim*2),dims_MD_new(Ndim*2),dims_ref(Ndim+1),idx_ref(Ndim+1),offsets_ref(Ndim+1),idx_ref_scalar, dims_in(Ndim),idx_in(Ndim), dims_MD3(Ndim),idx_MD3(Ndim), idx_MD(Ndim), dims_ref_scalar,dims_MD3_scalar

    integer i, j, ii, iii, jj, jjj, ll, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_r0, nc(Ndim), idxc(Ndim), nc0(Ndim), idxc0(Ndim), idx_aux(Ndim), dims_aux(Ndim), inc_c, inc_c0, inc_r, inc_r0
    integer level_butterfly, num_blocks, level_half, levelm
    integer group_m(Ndim), group_n(Ndim), mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_jj_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, Nvec
    integer nn1, nn2, mm1, mm2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
    integer pid, pgno_sub, index_c_scalar, dim_i, dim_ii
    DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)
    real(kind=8)::flop,flops

    level_butterfly = blocks%level_butterfly
    level_half = blocks%level_half
    levelm = level_half
    offsets_ref=0

    call SingleIndexToMultiIndex(Ndim, blocks%nr_m, bb_m, idx_r_m)
    group_m = blocks%row_group
    group_m = group_m*2**levelm - 1 + idx_r_m + blocks%idx_r_m - 1

   if (level == 0) then
      do dim_ii=1,Ndim
         if(dim_ii==1)then
            do dim_i=1,Ndim
               dims_ref(dim_i) = BFvec%vec(1)%index_MD(dim_i,BFvec%vec(1)%num_row+1,1)
            enddo
            dims_ref(Ndim+1)=Nvec
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,Ndim,dim_ii,BFvec%vec(1)%blocks(1,1)%matrix, 'N', mat, 'N',1)
         else
            do dim_i=1,Ndim
               if(dim_i<dim_ii)then
                  dims_ref(dim_i) = BFvec%vec(0)%index_MD(dim_i,BFvec%vec(0)%num_row+1,1)
               else
                  dims_ref(dim_i) = BFvec%vec(1)%index_MD(dim_i,BFvec%vec(1)%num_row+1,1)
               endif
            enddo
            dims_ref(Ndim+1)=Nvec
            call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,dim_ii-1,dim_ii,BFvec%vec(0)%blocks(1,1)%matrix, 'N', mat, 'N',1)
         endif

         if(associated(BFvec%vec(level)%blocks(1,1)%matrix))deallocate(BFvec%vec(level)%blocks(1,1)%matrix)
         allocate(BFvec%vec(level)%blocks(1,1)%matrix(BFvec%vec(level)%index_MD(dim_ii,BFvec%vec(level)%num_row+1,1),size(mat,2)))
         BFvec%vec(level)%blocks(1,1)%matrix=0

         flops=0
#ifdef HAVE_TASKLOOP
       !  !$omp taskloop default(shared) private(index_i_s,index_i_k,mm,nn1,nvec1) reduction(+:flops)
         !$omp taskgroup task_reduction(+:flops)
#else
#ifdef HAVE_OPENMP
         !$omp parallel do default(shared) private(index_i_s,index_i_k,mm,nn1,nvec1) reduction(+:flops)
#endif
#endif
         do index_i_s=1,BFvec%vec(level)%num_row
            index_i_k = index_i_s
            mm = size(blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_ii)%matrix,1)
            nn1 = size(blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_ii)%matrix,2)
            nvec1 = size(mat,2)
#ifdef HAVE_TASKLOOP
            !$omp task default(shared) firstprivate(mm,nn1,nvec1,index_i_k)
#endif
            call gemmf77('N', 'N', mm, nvec1, nn1, BPACK_cone, blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_ii)%matrix, mm, mat(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,1)+1,1), size(mat,1), BPACK_czero, BFvec%vec(level)%blocks(1,1)%matrix(BFvec%vec(level)%index_MD(dim_ii,index_i_s,1)+1,1), BFvec%vec(level)%index_MD(dim_ii,BFvec%vec(level+1)%num_row+1,1))
            flops = flops + flops_gemm(mm, nvec1, nn1)
#ifdef HAVE_TASKLOOP
            !$omp end task
#endif
         enddo
#ifdef HAVE_TASKLOOP
         !   !$omp end taskloop
         !$omp end taskgroup
#else
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
#endif
#ifdef HAVE_TASKLOOP
            !$omp atomic
#endif
            stats%Flop_Tmp = stats%Flop_Tmp + flops
#ifdef HAVE_TASKLOOP
            !$omp end atomic
#endif
         deallocate(mat)
      enddo
   else
      nc=2**(level-1)
      nc0=2**(level)

      do index_j_s=1,product(nc)
         call SingleIndexToMultiIndex(Ndim, nc, index_j_s, idxc)
         dims_aux=2
         do j=1,product(dims_aux)
            call SingleIndexToMultiIndex(Ndim, dims_aux, j, idx_aux)
            idxc0=idxc*2-1 + idx_aux -1
            call MultiIndexToSingleIndex(Ndim, nc0, index_jj_s, idxc0)

            do dim_ii=1,Ndim
               if(dim_ii==1)then
                  do dim_i=1,Ndim
                     dims_ref(dim_i) = BFvec%vec(level+1)%index_MD(dim_i,BFvec%vec(level+1)%num_row+1,index_jj_s)
                  enddo
                  dims_ref(Ndim+1)=Nvec
                  call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,Ndim,dim_ii,BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix, 'N', mat, 'N',1)
               else
                  do dim_i=1,Ndim
                     if(dim_i<dim_ii)then
                        dims_ref(dim_i) = BFvec%vec(level)%index_MD(dim_i,BFvec%vec(level)%num_row+1,index_j_s)
                     else
                        dims_ref(dim_i) = BFvec%vec(level+1)%index_MD(dim_i,BFvec%vec(level+1)%num_row+1,index_jj_s)
                     endif
                  enddo
                  dims_ref(Ndim+1)=Nvec
                  call TensorUnfoldingReshape(Ndim+1,dims_ref,dims_ref,offsets_ref,dim_ii-1,dim_ii,BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix, 'N', mat, 'N',1)
               endif
               deallocate(BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix)
               allocate(BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix(BFvec%vec(level)%index_MD(dim_ii,BFvec%vec(level)%num_row+1,index_j_s),size(mat,2)))
               BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix=0

               flops=0
#ifdef HAVE_TASKLOOP
            !   !$omp taskloop default(shared) private(index_i_s,index_i_k,mm1,mm2,nn,nvec1) reduction(+:flops)
               !$omp taskgroup task_reduction(+:flops)
#else
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(index_i_s,index_i_k,mm1,mm2,nn,nvec1) reduction(+:flops)
#endif
#endif
               do index_i_s=1,BFvec%vec(level+1)%num_row
                  index_i_k = 2*index_i_s - 1
                  mm1 = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k, index_jj_s,dim_ii)%matrix, 1)
                  mm2 = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k+1, index_jj_s,dim_ii)%matrix, 1)
                  nn = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k, index_jj_s,dim_ii)%matrix, 2)
                  nvec1 = size(mat,2)
#ifdef HAVE_TASKLOOP
                  !$omp task default(shared) firstprivate(mm1,mm2,nn,nvec1,index_i_k)
#endif
                  call gemmf77('N', 'N', mm1, nvec1, nn, BPACK_cone, blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k, index_jj_s,dim_ii)%matrix, mm1, mat(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,index_jj_s)+1,1), size(mat,1), BPACK_czero, BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix(BFvec%vec(level)%index_MD(dim_ii,index_i_k,index_j_s)+1,1), BFvec%vec(level)%index_MD(dim_ii,BFvec%vec(level)%num_row+1,index_j_s))
                  flops = flops + flops_gemm(mm1, nvec1, nn)
                  call gemmf77('N', 'N', mm2, nvec1, nn, BPACK_cone, blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k+1, index_jj_s,dim_ii)%matrix, mm2, mat(BFvec%vec(level+1)%index_MD(dim_ii,index_i_s,index_jj_s)+1,1), size(mat,1), BPACK_czero, BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix(BFvec%vec(level)%index_MD(dim_ii,index_i_k+1,index_j_s)+1,1), BFvec%vec(level)%index_MD(dim_ii,BFvec%vec(level)%num_row+1,index_j_s))
                  flops = flops + flops_gemm(mm2, nvec1, nn)
#ifdef HAVE_TASKLOOP
                  !$omp end task
#endif
               enddo
#ifdef HAVE_TASKLOOP
          !  !$omp end taskloop
          !$omp end taskgroup
#else
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
#endif
#ifdef HAVE_TASKLOOP
               !$omp atomic
#endif
               stats%Flop_Tmp = stats%Flop_Tmp + flops
#ifdef HAVE_TASKLOOP
               !$omp end atomic
#endif
               deallocate(mat)
            enddo

            if(.not. associated(BFvec%vec(level)%blocks(1,index_j_s)%matrix))then
               allocate(BFvec%vec(level)%blocks(1,index_j_s)%matrix(size(BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix,1),size(BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix,2)))
               BFvec%vec(level)%blocks(1,index_j_s)%matrix=0
            endif
            BFvec%vec(level)%blocks(1,index_j_s)%matrix = BFvec%vec(level)%blocks(1,index_j_s)%matrix + BFvec%vec(level+1)%blocks(1,index_jj_s)%matrix
         enddo
      enddo
   endif

end subroutine BF_MD_block_mvp_multiply_left



subroutine BF_MD_block_extraction_multiply_oneblock_left(blocks, bb_m, Ndim, BFvec, idx_c_m, level ,dim_i, ptree,stats)

    implicit none
    type(matrixblock_MD)::blocks
    type(proctree)::ptree
    type(Hstat):: stats
    type(butterfly_vec) :: BFvec
    integer bb_m, level,nn, Ndim, idx_c_m(Ndim),idx_c_m_1(Ndim), dims(Ndim)

    integer ii, iii, jj, jjj, ll, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
    integer level_butterfly, num_blocks, level_half
    integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
    integer mm1, mm2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
    integer pid, pgno_sub, index_c_scalar, dim_i
    DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)


    level_butterfly = blocks%level_butterfly
    level_half = blocks%level_half

    idx_c_m_1=idx_c_m
    dims = 2**(level_butterfly-level_half)
    do ll= level_butterfly-level_half-1,level,-1
      idx_c_m_1=int((idx_c_m_1 + 1)/2)
      dims=dims/2
    enddo
    call MultiIndexToSingleIndex(Ndim,dims, index_c_scalar, idx_c_m_1)

    do index_i_s=1,BFvec%vec(level+1)%nr
      if (level == 0) then
         index_i_k = index_i_s
         rank = size(blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_i)%matrix, 2)
         mm = size(blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_i)%matrix, 1)
         allocate (BFvec%vec(level + 1)%blocks(index_i_s,1)%matrix(mm,rank))
         BFvec%vec(level + 1)%blocks(index_i_s,1)%matrix = blocks%ButterflyU(bb_m)%blocks(index_i_k,dim_i)%matrix
      else
         index_i_k = 2*index_i_s - 1 !index_ii is global index in BFvec%vec(level)

         mm1 = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k,index_c_scalar,dim_i)%matrix, 1)
         mm2 = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k+1,index_c_scalar,dim_i)%matrix, 1)
         nn = size(blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k,index_c_scalar,dim_i)%matrix, 2)
         nvec1 = size(BFvec%vec(level)%blocks(index_i_k,1)%matrix, 1)
         nvec2 = size(BFvec%vec(level)%blocks(index_i_k+1,1)%matrix, 1)


         allocate (mat1(nvec1,nn))
         if (nvec1 > 0) then
               mat1 = 0
               call gemmf77('N', 'N', nvec1, nn, mm1, BPACK_cone, BFvec%vec(level)%blocks(index_i_k,1)%matrix, nvec1, blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k,index_c_scalar,dim_i)%matrix, mm1, BPACK_czero, mat1, nvec1)
#ifdef HAVE_TASKLOOP
               !$omp atomic
#endif
               stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(nvec1, nn, mm1)
#ifdef HAVE_TASKLOOP
               !$omp end atomic
#endif
         endif

         allocate (mat2(nvec2,nn))
         if (nvec2 > 0) then
               mat2 = 0
               call gemmf77('N', 'N', nvec2, nn, mm2, BPACK_cone, BFvec%vec(level)%blocks(index_i_k+1,1)%matrix, nvec2, blocks%ButterflyKerl_L(bb_m,level_butterfly-level+1)%blocks(index_i_k+1,index_c_scalar,dim_i)%matrix, mm2, BPACK_czero, mat2, nvec2)
#ifdef HAVE_TASKLOOP
               !$omp atomic
#endif
               stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(nvec2, nn, mm2)
#ifdef HAVE_TASKLOOP
               !$omp end atomic
#endif
         endif
         allocate (BFvec%vec(level + 1)%blocks(index_i_s,1)%matrix(nvec1+nvec2, nn))
         BFvec%vec(level + 1)%blocks(index_i_s, 1)%matrix(1:nvec1,:) = mat1
         BFvec%vec(level + 1)%blocks(index_i_s, 1)%matrix(1+nvec1:nvec1+nvec2,:) = mat2
         deallocate (mat1)
         deallocate (mat2)
      endif
    enddo

end subroutine BF_MD_block_extraction_multiply_oneblock_left






subroutine BF_block_extraction_sort_oneblock(blocks, BFvec, level,nn,ptree)

   implicit none
   type(matrixblock)::blocks
   type(proctree)::ptree
   type(butterfly_vec) :: BFvec
   integer level,nn

   integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0
   integer level_butterfly, num_blocks, level_half
   integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
   integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
   integer pid, pgno_sub
   DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)

   level_butterfly = blocks%level_butterfly


   index_i_loc_s = BFvec%vec(level + 1)%index(nn, 1)
   index_j_loc_s = BFvec%vec(level + 1)%index(nn, 2)

   if (associated(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)) then

! if(ptree%MyID>=2)write(*,*)ptree%MyID,'fanibef',level+1,index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)

      gg1 = 0
      nvec1 = 0
      do jj = 1, size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, 2)
         gg = NINT(dble(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(1, jj)))

         if (gg < gg1) then
            nvec1 = jj - 1
            exit
         endif
         gg1 = gg
      enddo
      nvec2 = size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, 2) - nvec1

      mm = size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix, 1)
      allocate (mat1(mm, nvec1))
      if (nvec1 > 0) mat1 = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(:, 1:nvec1)
      allocate (mat2(mm, nvec2))
      if (nvec2 > 0) mat2 = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(:, 1 + nvec1:nvec1 + nvec2)
      deallocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix)

      !>**** filter out the columns of mat1 and mat2 that are not specified by BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%index
      allocate (mat(mm, nvec1 + nvec2))
      idx1 = 0
      idx2 = 0
      idx = 0
      gg1 = 0
      gg2 = 0
      do ii = 1, size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index, 1)
         gg = NINT(dble(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(ii, 1)))
         do while (gg1 <= gg)
            if (gg1 == gg) then
               idx = idx + 1
               mat(:, idx) = mat1(:, idx1)
            endif
            idx1 = idx1 + 1
            if (idx1 > nvec1) exit
            gg1 = NINT(dble(mat1(1, idx1)))
         enddo
         do while (gg2 <= gg)
            if (gg2 == gg) then
               idx = idx + 1
               mat(:, idx) = mat2(:, idx2)
            endif
            idx2 = idx2 + 1
            if (idx2 > nvec2) exit
            gg2 = NINT(dble(mat2(1, idx2)))
         enddo
      enddo

      allocate (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(mm, idx))
      BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = mat(:, 1:idx)

      deallocate (mat)
      deallocate (mat1)
      deallocate (mat2)
   endif

   ! if(ptree%MyID>=2)write(*,*)ptree%MyID,'fani',level+1,index_i_loc_s,index_j_loc_s,shape(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%matrix)

end subroutine BF_block_extraction_sort_oneblock



subroutine BF_block_extraction_multiply_oneblock_last(blocks, BFvec, inters, level,nn,ptree,msh,stats)

   implicit none
   type(matrixblock)::blocks
   type(proctree)::ptree
   type(Hstat):: stats
   type(butterfly_vec) :: BFvec
   integer level,nn
   type(intersect)::inters(:)
   type(mesh)::msh

   integer ii, iii, jj, jjj, nnn, gg, gg1, gg2, rank, ncol, nrow, iidx, jidx, pgno, comm, ierr, idx1, idx2, idx, idx_r, idx_c, idx_r0, idx_c0, nc, nc0, nr, nr0, inc_c, inc_c0, inc_r, inc_r0, idxc, idxr, nng, ri
   integer level_butterfly, num_blocks, level_half
   integer group_m, group_n, mm, index_i, index_i0, index_i_k, index_j_k, index_i_loc_k, index_i_s, index_j_s, index_i_loc_s, index_j, index_j0, index_j_loc_k, index_j_loc_s, ij, tt, nvec1, nvec2
   integer nn1, nn2, mmm, index_ii, index_jj, index_ii_loc, index_jj_loc, nnn1
   integer pid, pgno_sub
   DT, allocatable::mat1(:, :), mat2(:, :), mat(:, :)
   DT, allocatable::Vpartial(:, :)

   level_butterfly = blocks%level_butterfly
   index_i_loc_s = BFvec%vec(level + 1)%index(nn, 1)
   index_j_loc_s = BFvec%vec(level + 1)%index(nn, 2)

   index_i_s = (index_i_loc_s - 1)*BFvec%vec(level + 1)%inc_r + BFvec%vec(level + 1)%idx_r
   index_j_s = (index_j_loc_s - 1)*BFvec%vec(level + 1)%inc_c + BFvec%vec(level + 1)%idx_c

   index_ii_loc = (index_i_s - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
   index_jj_loc = (index_j_s - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

   idxr = 0
   idxc = 0

   do while (idxr < size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index, 1))
      idx = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(idxr + 1, 1)
      nr = 0
      do while (BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(idxr + 1 + nr, 1) == idx)
         nr = nr + 1
         if (idxr + 1 + nr > size(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index, 1)) exit
      enddo
      if (blocks%inters(idx)%nc == 0) then  ! this takes care of the case where the intersection idx has 0 columns (resulting from the 2D block-cyclic distribution before element_Zmn_block_user in BPACK_CheckError and BF_CheckError)
         nc = 0
         goto 111
      endif

      ! if(ptree%MyID>=2)write(*,*)'id',ptree%MyID,index_i_loc_s,index_j_loc_s,idx,'shape',shape(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix),allocated(BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix),idxc
      if(.not. associated(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix))then
         write(*,*)ptree%MyID,level,index_ii_loc, index_jj_loc,'nanina',BFvec%vec(level)%index
         endif
         idx1 = NINT(dble(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1, idxc + 1)))

      if (idx /= idx1) then
         if (blocks%inters(idx1)%nr == 0) then  ! this takes care of the case where the intersection idx has 0 rows (resulting from the 2D block-cyclic distribution before element_Zmn_block_user in BPACK_CheckError and BF_CheckError)
            nr = 0
            goto 111
         endif
      endif
      if(idx/=idx1)then
         write(*,*)ptree%MyID,idx,idx1,blocks%inters(idx1)%nr,blocks%inters(idx)%nc,'nani',BFvec%vec(level)%index
      endif
      call assert(idx == idx1, 'row and column intersection# not match')
      nc = 0
      do while (NINT(dble(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(1, idxc + 1 + nc))) == idx)
         nc = nc + 1
         if (idxc + 1 + nc > size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)) exit
      enddo

      rank = size(blocks%ButterflyU%blocks(index_i_loc_s)%matrix, 2)
      allocate (Vpartial(nr, nc))
      Vpartial = 0
      allocate (mat(nr, rank))

      nng = blocks%inters(idx)%idx
      do ii = 1, nr
         iii = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(idxr + ii, 2)
         mmm = inters(nng)%rows(blocks%inters(idx)%rows(blocks%inters(idx)%rows_loc(iii)))
         ri = mmm - msh%basis_group(findgroup(mmm, msh, level_butterfly, blocks%row_group))%head + 1
         mat(ii, :) = blocks%ButterflyU%blocks(index_i_loc_s)%matrix(ri, :)
      enddo

      call gemmf77('N', 'N', nr, nc, rank, BPACK_cone, mat, nr, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(3, idxc + 1), rank + 2, BPACK_czero, Vpartial, nr)
#ifdef HAVE_TASKLOOP
      !$omp atomic
#endif
      stats%Flop_Tmp = stats%Flop_Tmp + flops_gemm(nr, nc, rank)
#ifdef HAVE_TASKLOOP
      !$omp end atomic
#endif

      do ii = 1, nr
         iii = BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%index(idxr + ii, 2)
         do jj = 1, nc
            jjj = NINT(dble(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix(2, idxc + jj)))
            blocks%inters(idx)%dat_loc(iii, jjj) = Vpartial(ii, jj)
         enddo
      enddo

      ! do ii=1,nr
      ! iii = BFvec%vec(level+1)%blocks(index_i_loc_s,index_j_loc_s)%index(idxr+ii,2)
      ! iii=inters(nng)%rows(blocks%inters(idx)%rows(blocks%inters(idx)%rows_loc(iii)))-msh%basis_group(blocks%row_group)%head+1
      ! do jj=1,nc
      ! jjj = BFvec%vec(level)%blocks(index_i_loc_s,index_j_loc_s)%matrix(2,idxc+jj)
      ! jjj=inters(nng)%cols(blocks%inters(idx)%cols(jjj))-msh%basis_group(blocks%col_group)%head+1
      ! call BF_value(iii,jjj,blocks,val)
      ! write(*,*)'seq:',abs(val)
      ! enddo
      ! enddo

      deallocate (Vpartial)
      deallocate (mat)
111               idxr = idxr + nr
      idxc = idxc + nc
   enddo

end subroutine BF_block_extraction_multiply_oneblock_last




!>*** Find the group index of point idx at the (group%level+level) level
   integer function findgroup(idx, msh, level, group)

      implicit none
      integer idx, level, ll, group, group1
      type(mesh)::msh
      group1 = group
      if (idx < msh%basis_group(group1)%head .or. idx > msh%basis_group(group1)%tail) then
         findgroup = -1
      else
         do ll = 1, level
            if (idx <= msh%basis_group(2*group1)%tail) then
               group1 = group1*2
            else
               group1 = group1*2 + 1
            endif
         enddo
         findgroup = group1
      endif

   end function findgroup

!>*** Find the process group index of point idx in a group
   integer function findpggroup(idx, msh, ptree, group, pgno)

      implicit none
      integer idx, ll, group, group1
      type(mesh)::msh
      type(proctree)::ptree
      integer level_p, pgno, pgno1

      level_p = ptree%nlevel - GetTreelevel(pgno)

      group1 = group
      pgno1 = pgno
      if (idx < msh%basis_group(group1)%head .or. idx > msh%basis_group(group1)%tail) then
         findpggroup = -1
      else
         do ll = 1, level_p
            if (idx <= msh%basis_group(2*group1)%tail) then
               group1 = group1*2
               pgno1 = pgno1*2
            else
               group1 = group1*2 + 1
               pgno1 = pgno1*2 + 1
            endif
         enddo
         findpggroup = pgno1
      endif

   end function findpggroup

   subroutine BF_value(mi, nj, blocks, value)


      implicit none

      integer mm, nn, mi, nj, groupm_start, groupn_start, level_butterfly, flag
      integer i, j, ii, jj, rank, group_m, group_n, header_mm, header_nn, k, kk
      integer group, level, mii, njj, rank1, rank2, index_ij, level_blocks, flag1
      DT ctemp, value

      type(matrixblock) :: blocks
      type(vectorset), allocatable:: vectors_set(:)
      integer, allocatable :: group_index_mm(:), group_index_nn(:)

      level_butterfly = blocks%level_butterfly

      allocate (group_index_mm(0:level_butterfly), group_index_nn(0:level_butterfly))

      flag = 0; i = 0; k = 0
      do while (flag == 0)
         i = i + 1
         if (size(blocks%ButterflyU%blocks(i)%matrix, 1) + k >= mi) then
            flag = 1
         endif
         k = k + size(blocks%ButterflyU%blocks(i)%matrix, 1)
      enddo
      group_index_mm(0) = i
      mii = mi - k + size(blocks%ButterflyU%blocks(i)%matrix, 1)

      flag = 0; j = 0; k = 0
      do while (flag == 0)
         j = j + 1
         if (size(blocks%ButterflyV%blocks(j)%matrix, 1) + k >= nj) then
            flag = 1
         endif
         k = k + size(blocks%ButterflyV%blocks(j)%matrix, 1)
      enddo
      group_index_nn(0) = j
      njj = nj - k + size(blocks%ButterflyV%blocks(j)%matrix, 1)

      if (level_butterfly > 0) then
         group_index_mm(1) = group_index_mm(0)
         group_index_nn(1) = group_index_nn(0)
         do level = 1, level_butterfly - 1
            group_index_mm(level + 1) = int((group_index_mm(level) + 1)/2)
            group_index_nn(level + 1) = int((group_index_nn(level) + 1)/2)
         enddo
      endif

!     if (group_index_mm(0)/=group_m .or. group_index_nn(0)/=group_n) then
!         write (*,*) 'BF_value_func error1!'
!         pause
!         continue
!     endif

!     do level=0, level_butterfly
!         group_index_mm(level)=group_index_mm(level)-group_m*2**level+1
!         group_index_nn(level)=group_index_nn(level)-group_n*2**level+1
!     enddo

      allocate (vectors_set(0:level_butterfly))
      do level = 0, level_butterfly
         if (level == 0) then
            rank = size(blocks%ButterflyV%blocks(group_index_nn(0))%matrix, 2)
            allocate (vectors_set(level)%vector(rank))
            !!$omp parallel do default(shared) private(i)
            do i = 1, rank
               vectors_set(level)%vector(i) = blocks%ButterflyV%blocks(group_index_nn(0))%matrix(njj, i)
            enddo
            !!$omp end parallel do
            ! write(*,*)'seq: ',level, abs(sum(vectors_set(level)%vector))
         else
            rank1 = size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly - level + 1), group_index_nn(level))%matrix, 2)
            rank2 = size(blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly - level + 1), group_index_nn(level))%matrix, 1)
            allocate (vectors_set(level)%vector(rank2))
            !!$omp parallel do default(shared) private(i,j,ctemp)
            do i = 1, rank2
               ctemp = 0
               do j = 1, rank1
                  ctemp = ctemp + blocks%ButterflyKerl(level)%blocks(group_index_mm(level_butterfly - level + 1), group_index_nn(level))%matrix(i, j)*vectors_set(level - 1)%vector(j)
               enddo
               vectors_set(level)%vector(i) = ctemp
            enddo
            !!$omp end parallel do
            deallocate (vectors_set(level - 1)%vector)
            ! write(*,*)'seq: ',level, abs(sum(vectors_set(level)%vector))
         endif
         if (level == level_butterfly) then
            rank = size(vectors_set(level)%vector, 1)
            ctemp = 0
            !!$omp parallel do default(shared) private(i) reduction(+:ctemp)
            ! write(*,*)'seq: ', level, abs(sum(blocks%ButterflyU%blocks(group_index_mm(0))%matrix(mii,1:rank))),mii
            do i = 1, rank
               ctemp = ctemp + blocks%ButterflyU%blocks(group_index_mm(0))%matrix(mii, i)*vectors_set(level)%vector(i)
            enddo
            !!$omp end parallel do
            value = ctemp
            deallocate (vectors_set(level)%vector)
         endif
      enddo
      deallocate (vectors_set)

      return

   end subroutine BF_value

   subroutine BF_get_rank(block_i, ptree, level_o)


      implicit none
      type(matrixblock)::block_i

      integer i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, blocks, edge, patch, node, group, level_c
      integer::block_num, block_num_new, num_blocks, level_butterfly
      integer::ierr, levels, levelo
      integer, optional:: level_o
      type(proctree)::ptree
      integer lve, lvs

      block_i%rankmin = 100000
      block_i%rankmax = -100000

      if (IOwnPgrp(ptree, block_i%pgno)) then

         level_butterfly = block_i%level_butterfly
         num_blocks = 2**level_butterfly

         lvs = 0
         lve = level_butterfly + 1
         if (present(level_o)) then
            lvs = level_o
            lve = level_o
         endif

         do level = lvs, lve
            if (level == 0) then
               do jj = 1, block_i%ButterflyV%nblk_loc
                  nn = size(block_i%ButterflyV%blocks(jj)%matrix, 1)
                  rank = size(block_i%ButterflyV%blocks(jj)%matrix, 2)
                  block_i%rankmin = min(block_i%rankmin, rank)
                  block_i%rankmax = max(block_i%rankmax, rank)
               enddo
            elseif (level == level_butterfly + 1) then
               do jj = 1, block_i%ButterflyU%nblk_loc
                  mm = size(block_i%ButterflyU%blocks(jj)%matrix, 1)
                  rank = size(block_i%ButterflyU%blocks(jj)%matrix, 2)
                  block_i%rankmin = min(block_i%rankmin, rank)
                  block_i%rankmax = max(block_i%rankmax, rank)
               enddo
            else
               do ii = 1, block_i%ButterflyKerl(level)%nr
                  do jj = 1, block_i%ButterflyKerl(level)%nc
                     nn = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 2)
                     rank = size(block_i%ButterflyKerl(level)%blocks(ii, jj)%matrix, 1)
                     block_i%rankmin = min(block_i%rankmin, rank)
                     block_i%rankmax = max(block_i%rankmax, rank)
                  enddo
               enddo
            endif
         enddo

         call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(block_i%pgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, block_i%rankmin, 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(block_i%pgno)%Comm, ierr)

      endif

   end subroutine BF_get_rank





   subroutine BF_MD_get_rank(Ndim, blocks, ptree)


      implicit none
      type(proctree)::ptree
      integer Ndim,dim_i,dim, dim_MD(Ndim+2),dim_subtensor(Ndim*2),nc(Ndim),idx_MD(Ndim+2),index_r_vector(Ndim),index_r_scalar,index_c_vector(Ndim),index_c_scalar
      type(matrixblock_MD)::blocks

      integer bb, i, j, ii, jj, iii, jjj, index_ij, mm, nn, rank, index_i, index_j, levelm, index_i_m, index_j_m
      integer level, edge, patch, node, group, level_c
      integer::level_butterfly, level_half, level_final,ierr

      real(kind=8)::memory
      memory = 0

      blocks%rankmin = 100000
      blocks%rankmax = -100000

      if (IOwnPgrp(ptree, blocks%pgno)) then

         level_butterfly = blocks%level_butterfly
         level_half = blocks%level_half

         if (blocks%style == 2) then
            if (allocated(blocks%M_loc) .and. allocated(blocks%N_loc)) then
            if (product(blocks%M_loc)/=0 .and. product(blocks%N_loc)/=0) then

               do bb=1, product(blocks%nc_m)
               do level = 0, level_half
                  if (level == 0) then
                     do jj = 1, blocks%ButterflyV(bb)%nblk_loc
                        do dim_i=1,Ndim
                           rank = size(blocks%ButterflyV(bb)%blocks(jj,dim_i)%matrix, 2)
                           blocks%rankmin = min(blocks%rankmin, rank)
                           blocks%rankmax = max(blocks%rankmax, rank)
                        enddo
                     enddo
                  elseif (level == level_butterfly + 1) then
                     write(*,*)"should not reach level == level_butterfly + 1 for the right half of matrixblock_MD"
                  else
                     dim_MD(1:Ndim)=blocks%ButterflyKerl_R(bb,level)%nr
                     dim_MD(1+Ndim)=blocks%ButterflyKerl_R(bb,level)%nc(1)
                     dim_MD(2+Ndim)=Ndim
                     do index_ij=1,product(dim_MD)
                        call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
                        index_r_vector = idx_MD(1:Ndim)
                        dim = idx_MD(Ndim+2)
                        index_j = idx_MD(Ndim+1)
                        call MultiIndexToSingleIndex(Ndim,dim_MD(1:Ndim), index_r_scalar, index_r_vector)
                        rank = size(blocks%ButterflyKerl_R(bb,level)%blocks(index_r_scalar,index_j,dim)%matrix, 1)
                        blocks%rankmin = min(blocks%rankmin, rank)
                        blocks%rankmax = max(blocks%rankmax, rank)
                     enddo
                  endif
               enddo
               enddo

               level_final = level_half + 1
               do bb=1, product(blocks%nr_m)
               do level = level_butterfly + 1, level_final, -1
                  if (level == 0) then
                     write(*,*)"should not reach level == 0 for the left half of matrixblock_MD"
                  elseif (level == level_butterfly + 1) then
                     do ii = 1, blocks%ButterflyU(bb)%nblk_loc
                        do dim_i=1,Ndim
                           rank = size(blocks%ButterflyU(bb)%blocks(ii,dim_i)%matrix, 2)
                           blocks%rankmin = min(blocks%rankmin, rank)
                           blocks%rankmax = max(blocks%rankmax, rank)
                        enddo
                     enddo
                  else
                     dim_MD(1)=blocks%ButterflyKerl_L(bb,level)%nr(1)
                     dim_MD(2:1+Ndim)=blocks%ButterflyKerl_L(bb,level)%nc
                     dim_MD(2+Ndim)=Ndim
                     do index_ij=1,product(dim_MD)
                        call SingleIndexToMultiIndex(Ndim+2,dim_MD, index_ij, idx_MD)
                        index_i = idx_MD(1)
                        dim = idx_MD(Ndim+2)
                        index_c_vector = idx_MD(2:Ndim+1)
                        call MultiIndexToSingleIndex(Ndim,dim_MD(2:1+Ndim), index_c_scalar, index_c_vector)
                        rank = size(blocks%ButterflyKerl_L(bb,level)%blocks(index_i,index_c_scalar,dim)%matrix, 2)
                        blocks%rankmin = min(blocks%rankmin, rank)
                        blocks%rankmax = max(blocks%rankmax, rank)
                     enddo
                  endif
               enddo
               enddo
            endif
            endif
         endif
         call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, blocks%rankmin, 1, MPI_INTEGER, MPI_MIN, ptree%pgrp(blocks%pgno)%Comm, ierr)


      endif
   end subroutine BF_MD_get_rank








   subroutine BF_sym2asym(blocks)




      implicit none

      integer M, N, Nrnd, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      integer:: dimension_n, num_row, num_col, mn_min

      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)

      DTR, allocatable :: Singular(:)
      DT, allocatable :: UU(:, :), VV(:, :)

      if (allocated(blocks%ButterflyMiddle)) then

         group_m = blocks%row_group ! Note: row_group and col_group interchanged here
         group_n = blocks%col_group
         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         levelm = floor_safe(dble(level_butterfly)/2d0)

         call assert(level_butterfly >= 2, 'level_butterfly not correct')

         level = levelm
         num_groupm = blocks%ButterflyKerl(level)%num_row
         num_groupn = blocks%ButterflyKerl(level)%num_col

         ! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
         do ij = 1, num_groupm*(num_groupn/2)
            i = (ij - 1)/(num_groupn/2) + 1
            j = (mod(ij - 1, (num_groupn/2)) + 1)*2 - 1
            index_i = int((i + 1)/2)
            index_j = int((j + 1)/2)

            nn1 = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
            nn2 = size(blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix, 2)
            mm = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 1)

            allocate (matrixtemp(mm, nn1))
            matrixtemp = blocks%ButterflyKerl(level)%blocks(i, j)%matrix
            ! call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm,nn1,mm)
            call gemmf90(blocks%ButterflyMiddle(i, index_j)%matrix, mm, matrixtemp, mm, blocks%ButterflyKerl(level)%blocks(i, j)%matrix, mm, 'N', 'N', mm, nn1, mm, BPACK_cone, BPACK_czero)
            deallocate (matrixtemp)

            allocate (matrixtemp(mm, nn2))
            matrixtemp = blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix
            ! call gemm_omp(blocks%ButterflyMiddle(i,index_j)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,nn2,mm)
            call gemmf90(blocks%ButterflyMiddle(i, index_j)%matrix, mm, matrixtemp, mm, blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix, mm, 'N', 'N', mm, nn2, mm, BPACK_cone, BPACK_czero)
            deallocate (matrixtemp)

            deallocate (blocks%ButterflyMiddle(i, index_j)%matrix)
         enddo
         ! !$omp end parallel do

         deallocate (blocks%ButterflyMiddle)

         do level = 0, levelm - 1
            if (level == 0) then

               iijj = 0
               do j = 1, num_blocks
                  iijj = iijj + 1
                  dimension_n = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                  rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
                  mn_min = min(dimension_n, rank)

                  allocate (matrixtemp(rank, dimension_n))
                  allocate (UU(rank, mn_min))
                  allocate (VV(mn_min, dimension_n))
                  allocate (Singular(mn_min))

                  call copymatT(blocks%ButterflyV%blocks(j)%matrix, matrixtemp, dimension_n, rank)

                  call gesvd_robust(matrixtemp, Singular, UU, VV, rank, dimension_n, mn_min)
                  do ii = 1, mn_min
                     UU(:, ii) = UU(:, ii)*Singular(ii)
                  end do

                  deallocate (blocks%ButterflyV%blocks(j)%matrix)
                  allocate (blocks%ButterflyV%blocks(j)%matrix(dimension_n, mn_min))
                  call copymatT(VV, blocks%ButterflyV%blocks(j)%matrix, mn_min, dimension_n)

                  index_j = mod(iijj - 1, blocks%ButterflyKerl(level + 1)%num_col) + 1
                  index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level + 1)%num_col))
                  mm1 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, 1)
                  allocate (matrixtemp1(mm1, mn_min))
                  matrixtemp1 = 0
                  ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
                  call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, mm1, UU, rank, matrixtemp1, mm1, 'N', 'N', mm1, mn_min, rank, BPACK_cone, BPACK_czero)

                  deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix)
                  allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix(mm1, mn_min))
                  blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix = matrixtemp1
                  deallocate (matrixtemp1)

                  mm2 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, 1)
                  allocate (matrixtemp1(mm2, mn_min))
                  ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
                  call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, mm2, UU, rank, matrixtemp1, mm2, 'N', 'N', mm2, mn_min, rank, BPACK_cone, BPACK_czero)
                  deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix)
                  allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix(mm2, mn_min))
                  blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix = matrixtemp1
                  deallocate (matrixtemp1)

                  deallocate (matrixtemp)
                  deallocate (UU)
                  deallocate (VV)
                  deallocate (Singular)

               enddo
            else
               num_row = blocks%ButterflyKerl(level)%num_row
               num_col = blocks%ButterflyKerl(level)%num_col

               iijj = 0
               do i = 1, num_row
                  do j = 1, num_col, 2
                     iijj = iijj + 1
                     rank = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                     nn1 = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
                     nn2 = size(blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix, 2)
                     mn_min = min(nn1 + nn2, rank)

                     allocate (matrixtemp(rank, nn1 + nn2))
                     allocate (UU(rank, mn_min))
                     allocate (VV(mn_min, nn1 + nn2))
                     allocate (Singular(mn_min))

                     ! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
                     matrixtemp(1:rank, 1:nn1) = blocks%ButterflyKerl(level)%blocks(i, j)%matrix
                     ! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
                     matrixtemp(1:rank, 1 + nn1:nn2 + nn1) = blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix

                     call gesvd_robust(matrixtemp, Singular, UU, VV, rank, nn1 + nn2, mn_min)
                     do ii = 1, mn_min
                        UU(:, ii) = UU(:, ii)*Singular(ii)
                     end do

                     deallocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix)
                     allocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix(mn_min, nn1))
                     ! call copymatN(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
                     blocks%ButterflyKerl(level)%blocks(i, j)%matrix = VV(1:mn_min, 1:nn1)
                     deallocate (blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix)
                     allocate (blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix(mn_min, nn2))
                     ! call copymatN(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
                     blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix = VV(1:mn_min, 1 + nn1:nn2 + nn1)

                     index_j = mod(iijj - 1, blocks%ButterflyKerl(level + 1)%num_col) + 1
                     index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level + 1)%num_col))

                     mm1 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, 1)
                     allocate (matrixtemp1(mm1, mn_min))
                     ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
                     call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, mm1, UU, rank, matrixtemp1, mm1, 'N', 'N', mm1, mn_min, rank, BPACK_cone, BPACK_czero)
                     deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix)
                     allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix(mm1, mn_min))
                     blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix = matrixtemp1
                     deallocate (matrixtemp1)

                     mm2 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, 1)
                     allocate (matrixtemp1(mm2, mn_min))
                     ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
                     call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, mm2, UU, rank, matrixtemp1, mm2, 'N', 'N', mm2, mn_min, rank, BPACK_cone, BPACK_czero)
                     deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix)
                     allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix(mm2, mn_min))
                     blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix = matrixtemp1
                     deallocate (matrixtemp1)

                     deallocate (matrixtemp)
                     deallocate (UU)
                     deallocate (VV)
                     deallocate (Singular)

                  end do
               end do
            end if
         end do

      end if

   end subroutine BF_sym2asym




   subroutine BF_MoveSingular_Ker(blocks, chara, level_start, level_end, ptree, stats, tolerance)

      implicit none

      integer level_start, level_end, index_i, index_j, na, nb, index_start, num_vectors, num_vectors1, num_vectors2
      integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
      integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, ranknew, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, level_half, level_final, pgno_sub
      integer idx_r, inc_r, nr, idx_c, inc_c, nc
      integer idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0
      DT ctemp
      character chara
      type(matrixblock)::blocks
      type(proctree)::ptree
      integer pgno, comm, ierr, mn_min
      type(Hstat)::stats
      real(kind=8)::flop, flops,n1,n2, tolerance
      integer index_ii, index_jj, index_ii_loc, index_jj_loc, index_i_loc, index_i_loc_s, index_i_loc_k, index_j_loc, index_j_loc0, index_i_loc0, index_j_loc_s, index_j_loc_k
      DTR, allocatable :: Singular(:)
      DT, allocatable :: UU(:, :), VV(:, :)

      type(butterfly_vec) :: BFvec
      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :), Vout_tmp(:, :)

      n1 = MPI_Wtime()

      level_butterfly = blocks%level_butterfly
      pgno = blocks%pgno
      comm = ptree%pgrp(pgno)%comm
      if (comm == MPI_COMM_NULL) then
         write (*, *) 'ninin', pgno, comm == MPI_COMM_NULL, ptree%MyID
      endif
      if(level_butterfly<=1)return

      call assert(IOwnPgrp(ptree, pgno), 'I do not share this block!')
      call assert(blocks%style == 2 .and.level_butterfly>=1, 'BF should have at least one level in BF_MoveSingular_Ker')
      call assert(1<=level_start .and. level_start<=level_butterfly, 'it should be 1<=level_start<=level_butterfly in BF_MoveSingular_Ker')
      call assert(1<=level_end .and. level_end<=level_butterfly, 'it should be 1<=level_end<=level_butterfly in BF_MoveSingular_Ker')

#ifdef HAVE_TASKLOOP
      !$omp parallel
      !$omp single
#endif
      if (BF_checkNAN(blocks)) then
         write (*, *) 'NAN in 0 BF_block_MVP_dat'
         stop
      end if

      if (chara == 'N') then
         if(level_start<level_end)then
            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(level_start+1:level_end))

            do level = level_start, min(level_half,level_end)
               ! n1 = MPI_Wtime()
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')

               if(level/=level_end)then
               BFvec%vec(level + 1)%idx_r = idx_r
               BFvec%vec(level + 1)%inc_r = inc_r
               BFvec%vec(level + 1)%nr = nr
               BFvec%vec(level + 1)%idx_c = idx_c
               BFvec%vec(level + 1)%inc_c = inc_c
               BFvec%vec(level + 1)%nc = nc
               endif

               if (nr > 0 .and. nc > 0) then
                  if(level/=level_end)then
                  if (level /= level_butterfly + 1) then
                     BFvec%vec(level + 1)%num_row = 2**level
                     BFvec%vec(level + 1)%num_col = 2**(level_butterfly - level)
                  else
                     BFvec%vec(level + 1)%num_row = 2**level_butterfly
                     BFvec%vec(level + 1)%num_col = 1
                  endif
                  if (level_half /= level) then ! the last level doesn't require doubling block columns
                  if (nc == 1 .and. 2**(level_butterfly - level) > 1) then ! double the number of local block columns used for MPI communication
                     BFvec%vec(level + 1)%nc = 2
                     BFvec%vec(level + 1)%idx_c = BFvec%vec(level + 1)%idx_c - 1 + mod(BFvec%vec(level + 1)%idx_c, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level + 1)%blocks(BFvec%vec(level + 1)%nr, BFvec%vec(level + 1)%nc))
                  endif


                  if (level == 0) then
                  elseif (level == level_butterfly + 1) then
                  else
                     flops = 0
                     !!$omp taskloop default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop)
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of row-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of row-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level)


                        index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        nn1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                        nn2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix, 2)
                        rank = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)

                        if(level==level_start)then
                           allocate (matrixtemp(rank, nn1 + nn2))
                           matrixtemp(1:rank, 1:nn1) = blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix
                           matrixtemp(1:rank, 1 + nn1:nn2 + nn1) = blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k + 1)%matrix
                           num_vectors1 = nn1
                           num_vectors2 = nn2
                        else
                           index_ii_loc = (index_ii - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level)
                           index_jj_loc = (index_jj - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1
                           num_vectors1 = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                           num_vectors2 = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc+1)%matrix, 2)

                           allocate (matrixtemp(rank, num_vectors1 + num_vectors2))
                           matrixtemp=0
                           call gemmf77('N', 'N', rank, num_vectors1, nn1, BPACK_cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, rank, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn1, BPACK_czero, matrixtemp(1, 1), rank)
                           call gemmf77('N', 'N', rank, num_vectors2, nn2, BPACK_cone, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix, rank, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc+1)%matrix, nn2, BPACK_czero, matrixtemp(1, 1+num_vectors1), rank)
                        endif

                        if(level==level_end)then
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(rank,num_vectors1))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrixtemp(1:rank, 1:num_vectors1)
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(rank,num_vectors2))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix = matrixtemp(1:rank, 1+num_vectors1:num_vectors1+num_vectors2)
                        else
                           index_i_loc_s = (index_i - BFvec%vec(level + 1)%idx_r)/BFvec%vec(level + 1)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level+1)
                           index_j_loc_s = (index_j - BFvec%vec(level + 1)%idx_c)/BFvec%vec(level + 1)%inc_c + 1

                           mn_min = min(num_vectors1 + num_vectors2, rank)
                           allocate (UU(rank, mn_min))
                           allocate (VV(mn_min, num_vectors1 + num_vectors2))
                           allocate (Singular(mn_min))
                           call assert(.not. myisnan(fnorm(matrixtemp, rank, num_vectors1 + num_vectors2)), 'matrixtemp NAN at 4')

                           ! call SVD_Truncate(matrixtemp, rank, num_vectors1 + num_vectors2, mn_min, UU, VV, Singular, tolerance, BPACK_SafeUnderflow, ranknew)
                           call gesvd_robust(matrixtemp, Singular, UU, VV, rank, num_vectors1 + num_vectors2, mn_min)
                           ranknew = mn_min
                           call assert(.not. myisnan(sum(Singular)), 'Singular NAN at 4')

                           do ii = 1, ranknew
                              UU(:, ii) = UU(:, ii)*Singular(ii)
                           end do

                           allocate(BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank, ranknew))
                           BFvec%vec(level + 1)%blocks(index_i_loc_s, index_j_loc_s)%matrix = UU(1:rank,1:ranknew)
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(ranknew,num_vectors1))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = VV(1:ranknew, 1:num_vectors1)
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(ranknew,num_vectors2))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix = VV(1:ranknew, 1+num_vectors1:num_vectors1+num_vectors2)

                           deallocate(UU,VV,Singular)
                        endif
                        deallocate(matrixtemp)
                     enddo
                     !!$omp end taskloop
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                  endif
               endif
               if(level/=level_start) then
               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo
               endif

               ! n2 = MPI_Wtime()
               ! time_tmp = time_tmp + n2-n1

               if (level_half /= level .and. level/=level_end) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level + 1), stats, ptree, level, 'R', 'B')
               endif
            enddo

            call assert(level_end<=level_half + 1,'level_end can at most one more level over level_half')
            if(level_half + 1==level_end)then
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half, 'R', 'C', 0)
               level = level_end
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'C')
               if (level == 0) then
                  write (*, *) 'should not arrive here'
               elseif (level == level_butterfly + 1) then
               else
                  flops = 0

                  !!$omp taskloop default(shared) private(index_ij,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc,index_i_loc_s,index_i_loc_k, index_j_loc,index_j_loc_s,index_j_loc_k,ij,ii,jj,kk,i,j,index_i,index_j,mm,mm1,mm2,nn,nn1,nn2,flop)
                  do index_ij = 1, nr0*nc0
                     index_j_loc = (index_ij - 1)/nr0 + 1       !index_i_loc is local index of column-wise ordering at current level
                     index_i_loc = mod(index_ij - 1, nr0) + 1
                     index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of column-wise ordering at current level
                     index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                     index_ii_loc = (index_i - BFvec%vec(level)%idx_r)/BFvec%vec(level)%inc_r + 1  !index_ii_loc is local index in BFvec%vec(level)
                     index_jj_loc = (index_j - BFvec%vec(level)%idx_c)/BFvec%vec(level)%inc_c + 1

                     index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                     index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                     mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                     nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                     num_vectors = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                     allocate (matrixtemp(mm, num_vectors))
                     matrixtemp = 0
                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, matrixtemp, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_czero, flop=flop)
                     deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                     allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(mm, num_vectors))
                     blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrixtemp
                     deallocate(matrixtemp)

                     mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 1)
                     nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, 2)
                     num_vectors = size(BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                     allocate (matrixtemp(mm, num_vectors))
                     matrixtemp = 0
                     call gemmf90(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, mm, BFvec%vec(level)%blocks(index_ii_loc, index_jj_loc)%matrix, nn, matrixtemp, mm, 'N', 'N', mm, num_vectors, nn, BPACK_cone, BPACK_czero, flop=flop)
                     deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix)
                     allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(mm, num_vectors))
                     blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix = matrixtemp
                     deallocate(matrixtemp)

                  enddo
                  !!$omp end taskloop

               endif

               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo
            endif

            do level = level_start+1, level_end
               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo
               if (allocated(BFvec%vec(level)%blocks)) deallocate (BFvec%vec(level)%blocks)
            enddo
         endif

      elseif (chara == 'T') then
         level_butterfly = blocks%level_butterfly
         num_blocks = 2**level_butterfly
         level_half = blocks%level_half

         if(level_start>level_end)then
            level_butterfly = blocks%level_butterfly
            num_blocks = 2**level_butterfly
            level_half = blocks%level_half

            allocate (BFvec%vec(level_butterfly - level_start + 2:level_butterfly - level_end+1))

            do level = level_start, max(level_half+1,level_end), -1
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

               if(level/=level_end)then
               BFvec%vec(level_butterfly - level + 2)%idx_r = idx_r
               BFvec%vec(level_butterfly - level + 2)%inc_r = inc_r
               BFvec%vec(level_butterfly - level + 2)%nr = nr
               BFvec%vec(level_butterfly - level + 2)%idx_c = idx_c
               BFvec%vec(level_butterfly - level + 2)%inc_c = inc_c
               BFvec%vec(level_butterfly - level + 2)%nc = nc
               endif

               if (nr > 0 .and. nc > 0) then
                  if(level/=level_end)then
                  if (level /= 0) then
                     BFvec%vec(level_butterfly - level + 2)%num_row = 2**(level - 1)
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**(level_butterfly - level + 1)
                  else
                     BFvec%vec(level_butterfly - level + 2)%num_row = 1
                     BFvec%vec(level_butterfly - level + 2)%num_col = 2**level_butterfly
                  endif
                  if (level_half + 1 /= level) then ! the last level doesn't require doubling block rows
                  if (nr == 1 .and. 2**(level - 1) > 1) then ! double the number of local block rows used for MPI communication
                     BFvec%vec(level_butterfly - level + 2)%nr = 2
                     BFvec%vec(level_butterfly - level + 2)%idx_r = BFvec%vec(level_butterfly - level + 2)%idx_r - 1 + mod(BFvec%vec(level_butterfly - level + 2)%idx_r, 2)
                  endif
                  endif
                  allocate (BFvec%vec(level_butterfly - level + 2)%blocks(BFvec%vec(level_butterfly - level + 2)%nr, BFvec%vec(level_butterfly - level + 2)%nc))
                  endif

                  if (level == level_butterfly + 1) then
                  elseif (level == 0) then
                  else
                     flops = 0
                     !!$omp taskloop default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop)
                     do index_ij = 1, nr*nc
                        index_j_loc = (index_ij - 1)/nr + 1
                        index_i_loc = mod(index_ij - 1, nr) + 1  !index_i_loc is local index of column-wise ordering at current level
                        index_i = (index_i_loc - 1)*inc_r + idx_r  !index_i is global index of column-wise ordering at current level
                        index_j = (index_j_loc - 1)*inc_c + idx_c

                        index_ii = 2*index_i - 1; index_jj = int((index_j + 1)/2) !index_ii is global index in BFvec%vec(level_butterfly-level+1)

                        index_i_loc_k = (2*index_i - 1 - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                        index_j_loc_k = (index_j - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                        mm1 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                        mm2 = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k + 1, index_j_loc_k)%matrix, 1)
                        rank = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)

                        if(level==level_start)then
                           allocate (matrixtemp(mm1 + mm2,rank))
                           matrixtemp(1:mm1,1:rank) = blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix
                           matrixtemp(1 + mm1:mm2 + mm1,1:rank) = blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix
                           num_vectors1 = mm1
                           num_vectors2 = mm2
                        else
                           index_ii_loc = (index_ii - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                           index_jj_loc = (index_jj - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                           num_vectors1 = size(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                           num_vectors2 = size(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+1, index_jj_loc)%matrix, 2)


                           allocate (matrixtemp(num_vectors1 + num_vectors2, rank))
                           matrixtemp=0
                           call gemmf77('T', 'N', num_vectors1,rank, mm1, BPACK_cone, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm1, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm1, BPACK_czero, matrixtemp(1, 1), num_vectors1 + num_vectors2)

                           call gemmf77('T', 'N', num_vectors2,rank, mm2, BPACK_cone, BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc+1, index_jj_loc)%matrix, mm2, blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix, mm2, BPACK_czero, matrixtemp(1+num_vectors1, 1), num_vectors1 + num_vectors2)
                        endif

                        if(level==level_end)then
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(num_vectors1,rank))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrixtemp(1:num_vectors1,1:rank)

                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(num_vectors2,rank))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix = matrixtemp(1+num_vectors1:num_vectors1+num_vectors2,1:rank)
                        else
                           index_i_loc_s = (index_i - BFvec%vec(level_butterfly - level + 2)%idx_r)/BFvec%vec(level_butterfly - level + 2)%inc_r + 1 !index_i_loc_s is local index in BFvec%vec(level_butterfly-level+2)
                           index_j_loc_s = (index_j - BFvec%vec(level_butterfly - level + 2)%idx_c)/BFvec%vec(level_butterfly - level + 2)%inc_c + 1



                           mn_min = min(num_vectors1 + num_vectors2, rank)
                           allocate (UU(num_vectors1 + num_vectors2, mn_min))
                           allocate (VV(mn_min,rank))
                           allocate (Singular(mn_min))

                           call assert(.not. myisnan(fnorm(matrixtemp, num_vectors1 + num_vectors2,rank)), 'matrixtemp NAN at 4')

                           ! call SVD_Truncate(matrixtemp, num_vectors1 + num_vectors2, rank, mn_min, UU, VV, Singular, tolerance, BPACK_SafeUnderflow, ranknew)
                           call gesvd_robust(matrixtemp, Singular, UU, VV, num_vectors1 + num_vectors2, rank, mn_min)
                           ranknew=mn_min


                           call assert(.not. myisnan(sum(Singular)), 'Singular NAN at 4')

                           do ii = 1, ranknew
                              VV(ii, :) = VV(ii, :)*Singular(ii)
                           end do

                           allocate(BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix(rank,ranknew))
                           call copymatT(VV, BFvec%vec(level_butterfly - level + 2)%blocks(index_i_loc_s, index_j_loc_s)%matrix, ranknew,rank)
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(num_vectors1,ranknew))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = UU(1:num_vectors1,1:ranknew)
                           deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix)
                           allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix(num_vectors2,ranknew))
                           blocks%ButterflyKerl(level)%blocks(index_i_loc_k+1, index_j_loc_k)%matrix = UU(1+num_vectors1:num_vectors1+num_vectors2,1:ranknew)

                           deallocate(UU,VV,Singular)
                        endif
                        deallocate(matrixtemp)
                     enddo
                     !!$omp end taskloop
                     stats%Flop_Tmp = stats%Flop_Tmp + flops
                  endif
               endif

               if(level/=level_start) then
               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo
               endif

               if (level_half + 1 /= level .and. level/=level_end) then
                  call BF_exchange_matvec(blocks, BFvec%vec(level_butterfly - level + 2), stats, ptree, level, 'C', 'B')
               endif
            enddo

            call assert(level_end>=level_half,'level_end can be at most level_half')
            if(level_half ==level_end)then
               level = level_half
               call BF_all2all_vec_n_ker(blocks, BFvec%vec(level_butterfly - level_half + 1), stats, ptree, ptree%pgrp(blocks%pgno)%nproc, level_half + 1, 'C', 'R', 0)
               call GetLocalBlockRange(ptree, blocks%pgno, level, level_butterfly, idx_r0, inc_r0, nr0, idx_c0, inc_c0, nc0, 'R')
               if (level == level_butterfly + 1) then
               elseif (level == 0) then
               else

                  flops = 0

                  !!$omp taskloop default(shared) private(index_ij,ii,jj,kk,ctemp,i,j,index_i,index_j,index_i_loc,index_j_loc,index_ii,index_jj,index_ii_loc,index_jj_loc,index_i_loc_s,index_j_loc_s,index_i_loc_k,index_j_loc_k,mm,mm1,mm2,nn,nn1,nn2,flop)
                  do index_ij = 1, nr0*nc0
                     index_j_loc = (index_ij - 1)/nr0 + 1
                     index_i_loc = mod(index_ij - 1, nr0) + 1  !index_i_loc is local index of row-wise ordering at current level
                     index_i = (index_i_loc - 1)*inc_r0 + idx_r0  !index_i is global index of row-wise ordering at current level
                     index_j = (index_j_loc - 1)*inc_c0 + idx_c0

                     index_ii = int((index_i + 1)/2); index_jj = 2*index_j - 1 !index_ii is global index in BFvec%vec(level_butterfly-level+2)

                     index_ii_loc = (index_i - BFvec%vec(level_butterfly - level + 1)%idx_r)/BFvec%vec(level_butterfly - level + 1)%inc_r + 1 !index_ii_loc is local index in BFvec%vec(level_butterfly-level+1)
                     index_jj_loc = (index_j - BFvec%vec(level_butterfly - level + 1)%idx_c)/BFvec%vec(level_butterfly - level + 1)%inc_c + 1

                     index_i_loc_k = (index_i - blocks%ButterflyKerl(level)%idx_r)/blocks%ButterflyKerl(level)%inc_r + 1 !index_i_loc_k is local index of kernels at current level
                     index_j_loc_k = (2*index_j - 1 - blocks%ButterflyKerl(level)%idx_c)/blocks%ButterflyKerl(level)%inc_c + 1

                     mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 1)
                     nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, 2)
                     num_vectors = size(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                     allocate (matrixtemp(num_vectors,nn))
                     matrixtemp = 0
                     call gemmf90(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix, mm, matrixtemp, num_vectors, 'T', 'N', num_vectors,nn, mm, BPACK_cone, BPACK_cone, flop=flop)
                     deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix)
                     allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix(num_vectors,nn))
                     blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k)%matrix = matrixtemp
                     deallocate(matrixtemp)

                     mm = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix, 1)
                     nn = size(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix, 2)
                     num_vectors = size(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, 2)
                     allocate (matrixtemp(num_vectors,nn))
                     matrixtemp = 0
                     call gemmf90(BFvec%vec(level_butterfly - level + 1)%blocks(index_ii_loc, index_jj_loc)%matrix, mm, blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix, mm, matrixtemp, num_vectors, 'T', 'N', num_vectors,nn, mm, BPACK_cone, BPACK_cone, flop=flop)
                     deallocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix)
                     allocate(blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix(num_vectors,nn))
                     blocks%ButterflyKerl(level)%blocks(index_i_loc_k, index_j_loc_k+1)%matrix = matrixtemp
                     deallocate(matrixtemp)
                  enddo
                  !!$omp end taskloop
               endif

               do j = 1, BFvec%vec(level_butterfly - level + 1)%nc
                  do i = 1, BFvec%vec(level_butterfly - level + 1)%nr
                     if (associated(BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level_butterfly - level + 1)%blocks(i, j)%matrix)
                  enddo
               enddo
            endif

            do level = level_butterfly - level_start + 2, level_butterfly - level_end+1
               do j = 1, BFvec%vec(level)%nc
                  do i = 1, BFvec%vec(level)%nr
                     if (associated(BFvec%vec(level)%blocks(i, j)%matrix)) deallocate (BFvec%vec(level)%blocks(i, j)%matrix)
                  enddo
               enddo
               if (allocated(BFvec%vec(level)%blocks)) deallocate (BFvec%vec(level)%blocks)
            enddo
         endif
      endif

      if(allocated(BFvec%vec))deallocate (BFvec%vec)
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2-n1

      return

   end subroutine BF_MoveSingular_Ker





   subroutine BF_MoveSingulartoLeft(blocks)




      implicit none

      integer M, N, Nrnd, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      integer:: dimension_n, dimension_m, num_row, num_col, mn_min

      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)

      DTR, allocatable :: Singular(:)
      DT, allocatable :: UU(:, :), VV(:, :)

      group_m = blocks%row_group ! Note: row_group and col_group interchanged here
      group_n = blocks%col_group
      level_butterfly = blocks%level_butterfly
      num_blocks = 2**level_butterfly

      do level = 0, level_butterfly
         if (level == 0) then
            iijj = 0
            do j = 1, num_blocks
               iijj = iijj + 1
               dimension_n = size(blocks%ButterflyV%blocks(j)%matrix, 1)
               rank = size(blocks%ButterflyV%blocks(j)%matrix, 2)
               mn_min = min(dimension_n, rank)

               allocate (matrixtemp(rank, dimension_n))
               allocate (UU(rank, mn_min))
               allocate (VV(mn_min, dimension_n))
               allocate (Singular(mn_min))

               call copymatT(blocks%ButterflyV%blocks(j)%matrix, matrixtemp, dimension_n, rank)
               call assert(.not. myisnan(fnorm(matrixtemp, rank, dimension_n)), 'matrixtemp NAN at 3')

               call gesvd_robust(matrixtemp, Singular, UU, VV, rank, dimension_n, mn_min)
               call assert(.not. myisnan(sum(Singular)), 'Singular NAN at 3')

               do ii = 1, mn_min
                  UU(:, ii) = UU(:, ii)*Singular(ii)
               end do

               deallocate (blocks%ButterflyV%blocks(j)%matrix)
               allocate (blocks%ButterflyV%blocks(j)%matrix(dimension_n, mn_min))
               call copymatT(VV, blocks%ButterflyV%blocks(j)%matrix, mn_min, dimension_n)

               index_j = mod(iijj - 1, blocks%ButterflyKerl(level + 1)%num_col) + 1
               index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level + 1)%num_col))

               mm1 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, 1)
               allocate (matrixtemp1(mm1, mn_min))
               ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
               call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, mm1, UU, rank, matrixtemp1, mm1, 'N', 'N', mm1, mn_min, rank, BPACK_cone, BPACK_czero)

               deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix)
               allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix(mm1, mn_min))
               blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix = matrixtemp1
               deallocate (matrixtemp1)

               mm2 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, 1)
               allocate (matrixtemp1(mm2, mn_min))
               ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
               call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, mm2, UU, rank, matrixtemp1, mm2, 'N', 'N', mm2, mn_min, rank, BPACK_cone, BPACK_czero)
               deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix)
               allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix(mm2, mn_min))
               blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix = matrixtemp1
               deallocate (matrixtemp1)

               deallocate (matrixtemp)
               deallocate (UU)
               deallocate (VV)
               deallocate (Singular)

            enddo
         else
            num_row = blocks%ButterflyKerl(level)%num_row
            num_col = blocks%ButterflyKerl(level)%num_col

            iijj = 0
            do i = 1, num_row
               do j = 1, num_col, 2
                  iijj = iijj + 1
                  rank = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                  nn1 = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
                  nn2 = size(blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix, 2)
                  mn_min = min(nn1 + nn2, rank)

                  allocate (matrixtemp(rank, nn1 + nn2))
                  allocate (UU(rank, mn_min))
                  allocate (VV(mn_min, nn1 + nn2))
                  allocate (Singular(mn_min))

                  ! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:rank,1:nn1),rank,nn1)
                  matrixtemp(1:rank, 1:nn1) = blocks%ButterflyKerl(level)%blocks(i, j)%matrix
                  ! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,matrixtemp(1:rank,1+nn1:nn2+nn1),rank,nn2)
                  matrixtemp(1:rank, 1 + nn1:nn2 + nn1) = blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix
                  call assert(.not. myisnan(fnorm(matrixtemp, rank, nn1 + nn2)), 'matrixtemp NAN at 4')
                  call gesvd_robust(matrixtemp, Singular, UU, VV, rank, nn1 + nn2, mn_min)
                  call assert(.not. myisnan(sum(Singular)), 'Singular NAN at 4')

                  do ii = 1, mn_min
                     UU(:, ii) = UU(:, ii)*Singular(ii)
                  end do

                  deallocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix)
                  allocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix(mn_min, nn1))
                  ! call copymatN(VV(1:mn_min,1:nn1),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mn_min,nn1)
                  blocks%ButterflyKerl(level)%blocks(i, j)%matrix = VV(1:mn_min, 1:nn1)
                  deallocate (blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix)
                  allocate (blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix(mn_min, nn2))
                  ! call copymatN(VV(1:mn_min,1+nn1:nn2+nn1),blocks%ButterflyKerl(level)%blocks(i,j+1)%matrix,mn_min,nn2)
                  blocks%ButterflyKerl(level)%blocks(i, j + 1)%matrix = VV(1:mn_min, 1 + nn1:nn2 + nn1)

                  if (level /= level_butterfly) then
                     index_j = mod(iijj - 1, blocks%ButterflyKerl(level + 1)%num_col) + 1
                     index_i = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level + 1)%num_col))

                     mm1 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, 1)
                     allocate (matrixtemp1(mm1, mn_min))
                     ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2-1,index_j)%matrix,UU,matrixtemp1,mm1,mn_min,rank)

                     call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix, mm1, UU, rank, matrixtemp1, mm1, 'N', 'N', mm1, mn_min, rank, BPACK_cone, BPACK_czero)

                     deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix)
                     allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix(mm1, mn_min))
                     blocks%ButterflyKerl(level + 1)%blocks(index_i*2 - 1, index_j)%matrix = matrixtemp1
                     deallocate (matrixtemp1)

                     mm2 = size(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, 1)
                     allocate (matrixtemp1(mm2, mn_min))
                     ! call gemm_omp(blocks%ButterflyKerl(level+1)%blocks(index_i*2,index_j)%matrix,UU,matrixtemp1,mm2,mn_min,rank)
                     call gemmf90(blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix, mm2, UU, rank, matrixtemp1, mm2, 'N', 'N', mm2, mn_min, rank, BPACK_cone, BPACK_czero)

                     deallocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix)
                     allocate (blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix(mm2, mn_min))
                     blocks%ButterflyKerl(level + 1)%blocks(index_i*2, index_j)%matrix = matrixtemp1
                     deallocate (matrixtemp1)
                  else
                     mm1 = size(blocks%ButterflyU%blocks(i)%matrix, 1)
                     allocate (matrixtemp1(mm1, mn_min))
                     ! call gemm_omp(blocks%ButterflyU%blocks(i)%matrix,UU,matrixtemp1,mm1,mn_min,rank)
                     call gemmf90(blocks%ButterflyU%blocks(i)%matrix, mm1, UU, rank, matrixtemp1, mm1, 'N', 'N', mm1, mn_min, rank, BPACK_cone, BPACK_czero)
                     deallocate (blocks%ButterflyU%blocks(i)%matrix)
                     allocate (blocks%ButterflyU%blocks(i)%matrix(mm1, mn_min))
                     blocks%ButterflyU%blocks(i)%matrix = matrixtemp1
                     deallocate (matrixtemp1)
                  end if

                  deallocate (matrixtemp)
                  deallocate (UU)
                  deallocate (VV)
                  deallocate (Singular)

               end do
            end do
         end if
      end do

   end subroutine BF_MoveSingulartoLeft

   subroutine BF_MoveSingulartoRight(blocks)




      implicit none

      integer M, N, Nrnd, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, ij, iijj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, mm1, mm2, levelm
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      integer:: dimension_n, dimension_m, num_row, num_col, mn_min

      DT, allocatable::matrixtemp(:, :), matrixtemp1(:, :)

      DTR, allocatable :: Singular(:)
      DT, allocatable :: UU(:, :), VV(:, :)

      group_m = blocks%row_group ! Note: row_group and col_group interchanged here
      group_n = blocks%col_group
      level_butterfly = blocks%level_butterfly
      num_blocks = 2**level_butterfly

      do level = level_butterfly + 1, 1, -1
         if (level == level_butterfly + 1) then
            iijj = 0
            do i = 1, num_blocks
               iijj = iijj + 1
               dimension_m = size(blocks%ButterflyU%blocks(i)%matrix, 1)
               rank = size(blocks%ButterflyU%blocks(i)%matrix, 2)
               mn_min = min(dimension_m, rank)

               allocate (matrixtemp(dimension_m, rank))
               allocate (UU(dimension_m, mn_min))
               allocate (VV(mn_min, rank))
               allocate (Singular(mn_min))

               ! call copymatN(blocks%ButterflyU%blocks(i)%matrix,matrixtemp,dimension_m,rank)
               matrixtemp = blocks%ButterflyU%blocks(i)%matrix
               call assert(.not. myisnan(fnorm(matrixtemp, dimension_m, rank)), 'matrixtemp NAN at 1')

               call gesvd_robust(matrixtemp, Singular, UU, VV, dimension_m, rank, mn_min)
               call assert(.not. myisnan(sum(Singular)), 'Singular NAN at 1')

               do ii = 1, mn_min
                  VV(ii, :) = VV(ii, :)*Singular(ii)
               end do

               deallocate (blocks%ButterflyU%blocks(i)%matrix)
               allocate (blocks%ButterflyU%blocks(i)%matrix(dimension_m, mn_min))
               ! call copymatN(UU,blocks%ButterflyU%blocks(i)%matrix,dimension_m,mn_min)
               blocks%ButterflyU%blocks(i)%matrix = UU

               index_i = mod(iijj - 1, blocks%ButterflyKerl(level - 1)%num_row) + 1
               index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level - 1)%num_row))

               nn1 = size(blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix, 2)
               allocate (matrixtemp1(mn_min, nn1))
               ! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,nn1,rank)
               call gemmf90(VV, mn_min, blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix, rank, matrixtemp1, mn_min, 'N', 'N', mn_min, nn1, rank, BPACK_cone, BPACK_czero)

               deallocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix)
               allocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix(mn_min, nn1))
               blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix = matrixtemp1
               deallocate (matrixtemp1)

               nn2 = size(blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix, 2)
               allocate (matrixtemp1(mn_min, nn2))
               ! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,nn2,rank)
               call gemmf90(VV, mn_min, blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix, rank, matrixtemp1, mn_min, 'N', 'N', mn_min, nn2, rank, BPACK_cone, BPACK_czero)
               deallocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix)
               allocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix(mn_min, nn2))
               blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix = matrixtemp1
               deallocate (matrixtemp1)

               deallocate (matrixtemp)
               deallocate (UU)
               deallocate (VV)
               deallocate (Singular)

            enddo
         else
            num_row = blocks%ButterflyKerl(level)%num_row
            num_col = blocks%ButterflyKerl(level)%num_col

            iijj = 0
            do j = 1, num_col
               do i = 1, num_row, 2
                  iijj = iijj + 1
                  rank = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 2)

                  mm1 = size(blocks%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                  mm2 = size(blocks%ButterflyKerl(level)%blocks(i + 1, j)%matrix, 1)
                  mn_min = min(mm1 + mm2, rank)

                  allocate (matrixtemp(mm1 + mm2, rank))
                  allocate (UU(mm1 + mm2, mn_min))
                  allocate (VV(mn_min, rank))
                  allocate (Singular(mn_min))

                  ! call copymatN(blocks%ButterflyKerl(level)%blocks(i,j)%matrix,matrixtemp(1:mm1,1:rank),mm1,rank)
                  matrixtemp(1:mm1, 1:rank) = blocks%ButterflyKerl(level)%blocks(i, j)%matrix
                  ! call copymatN(blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,matrixtemp(1+mm1:mm2+mm1,1:rank),mm2,rank)
                  matrixtemp(1 + mm1:mm2 + mm1, 1:rank) = blocks%ButterflyKerl(level)%blocks(i + 1, j)%matrix
                  call assert(.not. myisnan(fnorm(matrixtemp, mm1 + mm2, rank)), 'matrixtemp NAN at 2')

                  call gesvd_robust(matrixtemp, Singular, UU, VV, mm1 + mm2, rank, mn_min)
                  ! if(myisnan(sum(Singular)).and. mm1+mm2<rank)then
                  ! write(*,*)mm1+mm2,rank,mm1+mm2>=rank,'rank too large?'
                  ! end if

                  ! call assert(.not. myisnan(sum(Singular)),'Singular NAN at 2')
                  if (myisnan(sum(Singular))) then
                     write (*, *) 'Singular NAN at 2', mm1 + mm2, rank
                     do ii = 1, mm1 + mm2
                        do jj = 1, rank
                           write (777, *) dble(matrixtemp(ii, jj)), aimag(cmplx(matrixtemp(ii, jj), kind=8)), abs(matrixtemp(ii, jj))
                        end do
                     end do
                     stop
                  end if

                  do ii = 1, mn_min
                     VV(ii, :) = VV(ii, :)*Singular(ii)
                  end do

                  deallocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix)
                  allocate (blocks%ButterflyKerl(level)%blocks(i, j)%matrix(mm1, mn_min))
                  ! call copymatN(UU(1:mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i,j)%matrix,mm1,mn_min)
                  blocks%ButterflyKerl(level)%blocks(i, j)%matrix = UU(1:mm1, 1:mn_min)
                  deallocate (blocks%ButterflyKerl(level)%blocks(i + 1, j)%matrix)
                  allocate (blocks%ButterflyKerl(level)%blocks(i + 1, j)%matrix(mm2, mn_min))
                  ! call copymatN(UU(1+mm1:mm2+mm1,1:mn_min),blocks%ButterflyKerl(level)%blocks(i+1,j)%matrix,mm2,mn_min)
                  blocks%ButterflyKerl(level)%blocks(i + 1, j)%matrix = UU(1 + mm1:mm2 + mm1, 1:mn_min)

                  if (level /= 1) then
                     index_i = mod(iijj - 1, blocks%ButterflyKerl(level - 1)%num_row) + 1
                     index_j = ceiling_safe(dble(iijj)/dble(blocks%ButterflyKerl(level - 1)%num_row))
                     nn1 = size(blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix, 2)

                     allocate (matrixtemp1(mn_min, nn1))
                     ! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2-1)%matrix,matrixtemp1,mn_min,nn1,rank)
                     call gemmf90(VV, mn_min, blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix, rank, matrixtemp1, mn_min, 'N', 'N', mn_min, nn1, rank, BPACK_cone, BPACK_czero)

                     deallocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix)
                     allocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix(mn_min, nn1))
                     blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2 - 1)%matrix = matrixtemp1
                     deallocate (matrixtemp1)

                     nn2 = size(blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix, 2)
                     allocate (matrixtemp1(mn_min, nn2))
                     ! call gemm_omp(VV,blocks%ButterflyKerl(level-1)%blocks(index_i,index_j*2)%matrix,matrixtemp1,mn_min,nn2,rank)
                     call gemmf90(VV, mn_min, blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix, rank, matrixtemp1, mn_min, 'N', 'N', mn_min, nn2, rank, BPACK_cone, BPACK_czero)

                     deallocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix)
                     allocate (blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix(mn_min, nn2))
                     blocks%ButterflyKerl(level - 1)%blocks(index_i, index_j*2)%matrix = matrixtemp1
                     deallocate (matrixtemp1)
                  else
                     nn1 = size(blocks%ButterflyV%blocks(j)%matrix, 1)
                     allocate (matrixtemp1(nn1, mn_min))
                     ! call gemmNT_omp(blocks%ButterflyV%blocks(j)%matrix,VV,matrixtemp1,nn1,mn_min,rank)
                     call gemmf90(blocks%ButterflyV%blocks(j)%matrix, nn1, VV, mn_min, matrixtemp1, nn1, 'N', 'T', nn1, mn_min, rank, BPACK_cone, BPACK_czero)
                     deallocate (blocks%ButterflyV%blocks(j)%matrix)
                     allocate (blocks%ButterflyV%blocks(j)%matrix(nn1, mn_min))
                     blocks%ButterflyV%blocks(j)%matrix = matrixtemp1
                     deallocate (matrixtemp1)
                  end if

                  deallocate (matrixtemp)
                  deallocate (UU)
                  deallocate (VV)
                  deallocate (Singular)

               end do
            end do
         end if
      end do

   end subroutine BF_MoveSingulartoRight

   subroutine BF_Init_blocks(level_butterfly, groupm, groupn, pgno, block_rand, msh, ptree)


      implicit none

      integer level_c, rowblock, kover
      integer i, j, k, level, num_blocks, blocks3, num_row, num_col, ii, jj, kk, level_butterfly, mm, nn
      integer dimension_max, dimension_m, dimension_n, blocks, groupm, groupm_start, groupn_start, groupn, index_j, index_i
      real(kind=8) a, b, c, d
      DT ctemp
      type(matrixblock)::block, block_rand
      DT, allocatable::matrixtemp1(:, :), UU(:, :), VV(:, :)
      DTR, allocatable:: Singular(:)
      type(mesh)::msh
      ! type(Hoption)::option
      type(proctree)::ptree
      integer level_final, level_half
      integer idx_r, inc_r, nr, idx_c, inc_c, nc
      integer pgno

      block_rand%level_butterfly = level_butterfly
      num_blocks = 2**level_butterfly

      ! level_half = BF_Switchlevel(level_butterfly,option%pat_comp)
      level_half = floor_safe(dble(level_butterfly)/2d0) ! from outer to inner
      block_rand%level_half = level_half

      block_rand%style = 2
      block_rand%row_group = groupm
      block_rand%col_group = groupn

      block_rand%M = msh%basis_group(block_rand%row_group)%tail - msh%basis_group(block_rand%row_group)%head + 1
      block_rand%N = msh%basis_group(block_rand%col_group)%tail - msh%basis_group(block_rand%col_group)%head + 1
      block_rand%headm = msh%basis_group(block_rand%row_group)%head
      block_rand%headn = msh%basis_group(block_rand%col_group)%head

      block_rand%pgno = pgno

      groupm_start = groupm*2**level_butterfly
      groupn_start = groupn*2**level_butterfly
      if (level_butterfly /= 0) then
         allocate (block_rand%ButterflyKerl(level_butterfly))
      endif

      !>****** row-wise ordering from right side
      do level = 0, level_half
         if (level_butterfly == 0) then
            if (level == 0) then
               block_rand%ButterflyV%idx = 1
               block_rand%ButterflyV%inc = 1
               block_rand%ButterflyV%nblk_loc = 1
               block_rand%ButterflyV%num_blk = num_blocks
            elseif (level == level_butterfly + 1) then
               block_rand%ButterflyU%idx = 1
               block_rand%ButterflyU%inc = 1
               block_rand%ButterflyU%nblk_loc = 1
               block_rand%ButterflyU%num_blk = num_blocks
            endif
         else
            call GetLocalBlockRange(ptree, block_rand%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'R')
            if (level == 0) then
               block_rand%ButterflyV%idx = idx_c
               block_rand%ButterflyV%inc = inc_c
               block_rand%ButterflyV%nblk_loc = nc
               block_rand%ButterflyV%num_blk = num_blocks
            elseif (level == level_butterfly + 1) then
               block_rand%ButterflyU%idx = idx_r
               block_rand%ButterflyU%inc = inc_r
               block_rand%ButterflyU%nblk_loc = nr
               block_rand%ButterflyU%num_blk = num_blocks
            else
               block_rand%ButterflyKerl(level)%num_row = 2**level
               block_rand%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
               block_rand%ButterflyKerl(level)%idx_r = idx_r
               block_rand%ButterflyKerl(level)%inc_r = inc_r
               block_rand%ButterflyKerl(level)%nr = nr
               block_rand%ButterflyKerl(level)%idx_c = idx_c*2 - 1
               block_rand%ButterflyKerl(level)%inc_c = inc_c
               block_rand%ButterflyKerl(level)%nc = nc*2
            endif
         endif
      enddo

      !>****** column-wise ordering from left side
      level_final = level_half + 1
      do level = level_butterfly + 1, level_final, -1
         if (level_butterfly == 0) then
            if (level == 0) then
               block_rand%ButterflyV%idx = 1
               block_rand%ButterflyV%inc = 1
               block_rand%ButterflyV%nblk_loc = 1
               block_rand%ButterflyV%num_blk = num_blocks
            elseif (level == level_butterfly + 1) then
               block_rand%ButterflyU%idx = 1
               block_rand%ButterflyU%inc = 1
               block_rand%ButterflyU%nblk_loc = 1
               block_rand%ButterflyU%num_blk = num_blocks
            endif
         else
            call GetLocalBlockRange(ptree, block_rand%pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, 'C')

            if (level == 0) then
               block_rand%ButterflyV%idx = idx_c
               block_rand%ButterflyV%inc = inc_c
               block_rand%ButterflyV%nblk_loc = nc
               block_rand%ButterflyV%num_blk = num_blocks
            elseif (level == level_butterfly + 1) then
               block_rand%ButterflyU%idx = idx_r
               block_rand%ButterflyU%inc = inc_r
               block_rand%ButterflyU%nblk_loc = nr
               block_rand%ButterflyU%num_blk = num_blocks
            else
               block_rand%ButterflyKerl(level)%num_row = 2**level
               block_rand%ButterflyKerl(level)%num_col = 2**(level_butterfly - level + 1)
               block_rand%ButterflyKerl(level)%idx_r = idx_r*2 - 1
               block_rand%ButterflyKerl(level)%inc_r = inc_r
               block_rand%ButterflyKerl(level)%nr = nr*2
               block_rand%ButterflyKerl(level)%idx_c = idx_c
               block_rand%ButterflyKerl(level)%inc_c = inc_c
               block_rand%ButterflyKerl(level)%nc = nc
            endif
         endif
      enddo

      return

   end subroutine BF_Init_blocks

   recursive subroutine Hmat_block_copy(trans, block2, block1, memory)


      implicit none

      integer blocks, flag_recv, count1, count2, recv_count, mm, nn, length
      integer i, ii, j, jj, style, send_ID, group_m, group_n, indices, requests
      character chara

      type(matrixblock), pointer :: block1, block2, blocks_son1, blocks_son2
      character::trans
      real(kind=8), optional::memory
      real(kind=8)::memory_tmp

      block2%style = block1%style

      block2%level = block1%level
      block2%row_group = block1%row_group
      block2%col_group = block1%col_group
      block2%level_butterfly = 0
      group_m = block2%row_group
      group_n = block2%col_group
      block2%pgno = block1%pgno
      block2%M = block1%M
      block2%N = block1%N
      block2%headm = block1%headm
      block2%headn = block1%headn
      block2%M_loc = block1%M_loc
      block2%N_loc = block1%N_loc

      if (associated(block1%N_p)) then
         if (associated(block2%N_p)) deallocate (block2%N_p)
         allocate (block2%N_p(size(block1%N_p, 1), 2))
         block2%N_p = block1%N_p
      endif
      if (associated(block1%M_p)) then
         if (associated(block2%M_p)) deallocate (block2%M_p)
         allocate (block2%M_p(size(block1%M_p, 1), 2))
         block2%M_p = block1%M_p
      endif

      style = block2%style
      if (style == 4) then
         allocate (block2%sons(2, 2))
         do j = 1, 2
            do i = 1, 2
               block2%sons(i, j)%father => block2
            enddo
         enddo

         blocks_son1 => block1%sons(1, 1)
         blocks_son2 => block2%sons(1, 1)
         call Hmat_block_copy(trans, blocks_son2, blocks_son1, memory)
         blocks_son1 => block1%sons(2, 1)
         blocks_son2 => block2%sons(2, 1)
         call Hmat_block_copy(trans, blocks_son2, blocks_son1, memory)
         blocks_son1 => block1%sons(1, 2)
         blocks_son2 => block2%sons(1, 2)
         call Hmat_block_copy(trans, blocks_son2, blocks_son1, memory)
         blocks_son1 => block1%sons(2, 2)
         blocks_son2 => block2%sons(2, 2)
         call Hmat_block_copy(trans, blocks_son2, blocks_son1, memory)

      else
         call BF_copy(trans, block1, block2, memory_tmp)
         if (present(memory)) memory = memory + memory_tmp
      endif

      return

   end subroutine Hmat_block_copy

   recursive subroutine Hmat_block_delete(blocks)

      implicit none

      integer level_actual, num_col, num_row
      integer i, j, mm, nn, rank, num_blocks, level, level_butterfly
      real*8 memory_butterfly, rtemp
      type(matrixblock) :: blocks
      type(matrixblock), pointer :: blocks_son

      if (blocks%style == 4) then

         blocks_son => blocks%sons(1, 1)
         call Hmat_block_delete(blocks_son)
         blocks_son => blocks%sons(2, 1)
         call Hmat_block_delete(blocks_son)
         blocks_son => blocks%sons(1, 2)
         call Hmat_block_delete(blocks_son)
         blocks_son => blocks%sons(2, 2)
         call Hmat_block_delete(blocks_son)

         deallocate (blocks%sons)
         call BF_delete(blocks, 1)
      else
         call BF_delete(blocks, 1)
      endif

      return

   end subroutine Hmat_block_delete

   recursive subroutine Hmat_block_ComputeMemory(blocks, memory)

      implicit none

      integer level_actual, num_col, num_row
      integer i, j, mm, nn, rank, num_blocks, level, level_butterfly
      real*8 memory_butterfly, rtemp, memory
      type(matrixblock) :: blocks
      type(matrixblock), pointer :: blocks_son

      if (blocks%style == 4) then

         blocks_son => blocks%sons(1, 1)
         call Hmat_block_ComputeMemory(blocks_son, memory)
         blocks_son => blocks%sons(2, 1)
         call Hmat_block_ComputeMemory(blocks_son, memory)
         blocks_son => blocks%sons(1, 2)
         call Hmat_block_ComputeMemory(blocks_son, memory)
         blocks_son => blocks%sons(2, 2)
         call Hmat_block_ComputeMemory(blocks_son, memory)
      else
         call BF_ComputeMemory(blocks, rtemp)
         memory = memory + rtemp
      endif

      return

   end subroutine Hmat_block_ComputeMemory

   recursive subroutine Hmat_Lsolve(blocks_l, trans, idx_start, nvec, Vinout, ld, ptree, stats)
      implicit none

      ! integer vectors_y
      integer style(3)
      integer i, j, k, ii, zfpflag
      integer mm, nn, nvec, idxs_m, idx_start ! idx_start means the global indice of the first element of Vinout
      integer head, tail
      DT ctemp
      integer ld
      DT:: Vinout(ld, *)
      type(matrixblock) :: blocks_l !!!! modified by Yang Liu. passing pointer is dangerous, blocks_u row/row_group becomes different once in this subroutine
      character trans ! 'N' means multiple L^-1 from left, 'T' means multiple L^-1 from right
      type(proctree)::ptree
      type(Hstat)::stats
      real(kind=8)::tol_used

      if (blocks_l%style == 4) then
         if (trans == 'N') then
            call Hmat_Lsolve(blocks_l%sons(1, 1), trans, idx_start, nvec, Vinout, ld, ptree, stats)
            call Hmat_block_MVP_dat(blocks_l%sons(2, 1), trans, idx_start, idx_start, nvec, Vinout, ld, Vinout, ld, -BPACK_cone, ptree, stats)
            call Hmat_Lsolve(blocks_l%sons(2, 2), trans, idx_start, nvec, Vinout, ld, ptree, stats)
         else
            call Hmat_Lsolve(blocks_l%sons(2, 2), trans, idx_start, nvec, Vinout, ld, ptree, stats)
            call Hmat_block_MVP_dat(blocks_l%sons(2, 1), trans, idx_start, idx_start, nvec, Vinout, ld, Vinout, ld, -BPACK_cone, ptree, stats)
            call Hmat_Lsolve(blocks_l%sons(1, 1), trans, idx_start, nvec, Vinout, ld, ptree, stats)
         end if
      else
         mm = blocks_l%M
         idxs_m = blocks_l%headm - idx_start + 1

         if (trans == 'N') then
            do i = 1, mm
               ii = blocks_l%ipiv(i)
               if (ii /= i) then
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(j,ctemp)
#endif
                  do j = 1, nvec
                     ctemp = Vinout(idxs_m + i - 1, j)
                     Vinout(idxs_m + i - 1, j) = Vinout(idxs_m + ii - 1, j)
                     Vinout(idxs_m + ii - 1, j) = ctemp
                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
               endif
            enddo
         endif
         ! write(*,*)blocks_l%level,blocks_l%pgno,ptree%MyID,blocks_l%headm,mm,idx_start,'daha'
#if HAVE_ZFP
         zfpflag=0
         if(allocated(blocks_l%FullmatZFP%buffer_r))zfpflag=1
         if(zfpflag==1)call ZFP_Decompress(blocks_l%fullmat,blocks_l%FullmatZFP,blocks_l%M,blocks_l%N,tol_used,1)
#endif
         call trsmf90(blocks_l%fullmat, Vinout(idxs_m:idxs_m + mm - 1, 1:nvec), 'L', 'L', trans, 'U', mm, nvec)
#if HAVE_ZFP
         if(zfpflag==1)call ZFP_Compress(blocks_l%fullmat,blocks_l%FullmatZFP,blocks_l%M,blocks_l%N,tol_used,1)
#endif
         if (trans /= 'N') then
            do i = mm, 1, -1
               ii = blocks_l%ipiv(i)
               if (ii /= i) then
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(j,ctemp)
#endif
                  do j = 1, nvec
                     ctemp = Vinout(idxs_m + i - 1, j)
                     Vinout(idxs_m + i - 1, j) = Vinout(idxs_m + ii - 1, j)
                     Vinout(idxs_m + ii - 1, j) = ctemp
                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
               endif
            enddo
         end if
      endif

      return

   end subroutine Hmat_Lsolve

   recursive subroutine Hmat_Usolve(blocks_u, trans, idx_start, nvec, Vinout, ld, ptree, stats)
      implicit none

      type(proctree)::ptree
      type(Hstat)::stats

      integer vectors_x, vectors_y
      integer style(3), mark
      integer i, j, k, ii, zfpflag
      integer mm, nn, nvec
      integer head, tail
      integer ld
      DT Vinout(ld, *)
      type(matrixblock) :: blocks_u, blocks !!!! modified by Yang Liu. passing pointer is dangerous, blocks_u row/row_group becomes different once in this subroutine
      character trans
      integer idx_start, idxs_m
      real(kind=8)::tol_used

      mark = 0
      if (blocks_u%style == 4) then
         if (trans == 'N') then
            call Hmat_Usolve(blocks_u%sons(2, 2), trans, idx_start, nvec, Vinout, ld, ptree, stats)
            call Hmat_block_MVP_dat(blocks_u%sons(1, 2), trans, idx_start, idx_start, nvec, Vinout, ld, Vinout, ld, -BPACK_cone, ptree, stats)
            call Hmat_Usolve(blocks_u%sons(1, 1), trans, idx_start, nvec, Vinout, ld, ptree, stats)
         else
            call Hmat_Usolve(blocks_u%sons(1, 1), trans, idx_start, nvec, Vinout, ld, ptree, stats)
            call Hmat_block_MVP_dat(blocks_u%sons(1, 2), trans, idx_start, idx_start, nvec, Vinout, ld, Vinout, ld, -BPACK_cone, ptree, stats)
            call Hmat_Usolve(blocks_u%sons(2, 2), trans, idx_start, nvec, Vinout, ld, ptree, stats)
         end if

      else
         mm = blocks_u%M
         idxs_m = blocks_u%headm - idx_start + 1
#if HAVE_ZFP
         zfpflag=0
         if(allocated(blocks_u%FullmatZFP%buffer_r))zfpflag=1
         if(zfpflag==1)call ZFP_Decompress(blocks_u%fullmat,blocks_u%FullmatZFP,blocks_u%M,blocks_u%N,tol_used,1)
#endif
         call trsmf90(blocks_u%fullmat, Vinout(idxs_m:idxs_m + mm - 1, 1:nvec), 'L', 'U', trans, 'N', mm, nvec)
#if HAVE_ZFP
         if(zfpflag==1)call ZFP_Compress(blocks_u%fullmat,blocks_u%FullmatZFP,blocks_u%M,blocks_u%N,tol_used,1)
#endif
      endif

      return

   end subroutine Hmat_Usolve

   recursive subroutine Hmat_block_MVP_dat(blocks, trans, idx_start_m, idx_start_n, Nrnd, Vin, ldi, Vout, ldo, a, ptree, stats,level_start,level_end)

      implicit none
      integer,optional:: level_start, level_end
      integer idx_start_m, idx_start_n
      integer Nrnd
      integer mm, nn, idxs_m, idxs_n, zfpflag
      DT a
      character trans
      type(matrixblock)::blocks
      type(matrixblock), pointer::blocks_son
      integer:: style
      DT, allocatable::Vintmp(:, :), Vouttmp(:, :)
      integer ldi, ldo, flag
      DT::Vin(ldi, *), Vout(ldo, *)
      real(kind=8)::tol_used
      type(proctree)::ptree
      type(Hstat)::stats

      style = blocks%style
      mm = blocks%M
      idxs_m = blocks%headm - idx_start_m + 1
      nn = blocks%N
      idxs_n = blocks%headn - idx_start_n + 1

      if (style == 4) then
         blocks_son => blocks%sons(1, 1)
         call Hmat_block_MVP_dat(blocks_son, trans, idx_start_m, idx_start_n, Nrnd, Vin, ldi, Vout, ldo, a, ptree, stats,level_start,level_end)
         blocks_son => blocks%sons(1, 2)
         call Hmat_block_MVP_dat(blocks_son, trans, idx_start_m, idx_start_n, Nrnd, Vin, ldi, Vout, ldo, a, ptree, stats,level_start,level_end)
         blocks_son => blocks%sons(2, 1)
         call Hmat_block_MVP_dat(blocks_son, trans, idx_start_m, idx_start_n, Nrnd, Vin, ldi, Vout, ldo, a, ptree, stats,level_start,level_end)
         blocks_son => blocks%sons(2, 2)
         call Hmat_block_MVP_dat(blocks_son, trans, idx_start_m, idx_start_n, Nrnd, Vin, ldi, Vout, ldo, a, ptree, stats,level_start,level_end)
      else
         flag=1
         if(present(level_start) .and. present(level_end))then
         if(blocks%level>=level_start .and. blocks%level<=level_end)then
            if(style==1 .and. blocks%level==level_end)then
               flag=0    ! only multiply with dense blocks when level_end = h_mat%Maxlevel+1
            else
               flag=1
            endif
         else
            flag=0
         endif
         endif
         if(flag==1)then
         if (style == 1) then
#if HAVE_ZFP
            zfpflag=0
            if(allocated(blocks%FullmatZFP%buffer_r))zfpflag=1
            if(zfpflag==1)call ZFP_Decompress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,1)
#endif
            if (trans == 'N') then
               allocate (Vintmp(nn, Nrnd))
               Vintmp = Vin(idxs_n:idxs_n + nn - 1, 1:Nrnd)
               allocate (Vouttmp(mm, Nrnd))
               Vouttmp = 0
               call gemmf90(blocks%fullmat, mm, Vintmp, nn, Vouttmp, mm, trans, 'N', mm, Nrnd, nn, a, BPACK_czero)
               Vout(idxs_m:idxs_m + mm - 1, 1:Nrnd) = Vout(idxs_m:idxs_m + mm - 1, 1:Nrnd) + Vouttmp
               deallocate (Vintmp)
               deallocate (Vouttmp)
            else
               allocate (Vintmp(mm, Nrnd))
               Vintmp = Vin(idxs_m:idxs_m + mm - 1, 1:Nrnd)
               allocate (Vouttmp(nn, Nrnd))
               Vouttmp = 0
               call gemmf90(blocks%fullmat, mm, Vintmp, mm, Vouttmp, nn, trans, 'N', nn, Nrnd, mm, a, BPACK_czero)
               Vout(idxs_n:idxs_n + nn - 1, 1:Nrnd) = Vout(idxs_n:idxs_n + nn - 1, 1:Nrnd) + Vouttmp
               deallocate (Vintmp)
               deallocate (Vouttmp)
            endif
#if HAVE_ZFP
            if(zfpflag==1)call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,1)
#endif
         else
            if (trans == 'N') then
               call BF_block_MVP_dat(blocks, trans, mm, nn, Nrnd, Vin(idxs_n, 1), ldi, Vout(idxs_m, 1), ldo, a, BPACK_cone, ptree, stats)
            else
               call BF_block_MVP_dat(blocks, trans, mm, nn, Nrnd, Vin(idxs_m, 1), ldi, Vout(idxs_n, 1), ldo, a, BPACK_cone, ptree, stats)
            endif
         endif
         endif
      endif

   end subroutine Hmat_block_MVP_dat



   !>**** Multiply with dense blocks.
   subroutine Full_block_MVP_dat(blocks, chara, M, num_vectors, random1, ldi, random2, ldo, a, b)



      implicit none

      integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2
      DT ctemp, a, b
      character chara
      type(matrixblock)::blocks
      integer M, M1, N1, zfpflag
      integer ldi,ldo
      real(kind=8)::tol_used
      DT :: random1(ldi, *), random2(ldo, *)
      DT:: al, be
      DT, allocatable :: random2tmp(:, :)
#if HAVE_ZFP
      zfpflag=0
      if(allocated(blocks%FullmatZFP%buffer_r))zfpflag=1
      if(zfpflag==1)call ZFP_Decompress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,1)
#endif
      M1=size(blocks%fullmat, 1)
      N1=size(blocks%fullmat, 2)


      al = 1d0
      be = 0d0

      call assert(M1 == M, 'M not equal fullmat dim')

      if (chara == 'N') then
         allocate (random2tmp(M1, num_vectors))
         random2tmp = random2(1:M1, 1:num_vectors)
         call gemmf90(blocks%fullmat, M, random1, ldi, random2tmp, M1, 'N', 'N', M1, num_vectors, N1, BPACK_cone, BPACK_czero)
         random2(1:M1, 1:num_vectors) = a*random2tmp + b*random2(1:M1, 1:num_vectors)
      elseif (chara == 'T') then
         allocate (random2tmp(N1, num_vectors))
         random2tmp = random2(1:N1, 1:num_vectors)
         call gemmf90(blocks%fullmat, M, random1, ldi, random2tmp, N1, 'T', 'N', N1, num_vectors, M1, al, be)
         random2(1:N1, 1:num_vectors) = a*random2tmp + b*random2(1:N1, 1:num_vectors)
      end if
#if HAVE_ZFP
      if(zfpflag==1)call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,blocks%M,blocks%N,tol_used,1)
#endif
      ! write(*,*)'wo cao ni ma'
      deallocate (random2tmp)
   end subroutine Full_block_MVP_dat


   !>**** Multiply with dense blocks (as tensor). This is the same as Full_block_MVP_dat, except that blocks needs to be type(matrixblock_MD)
   subroutine Full_block_MD_MVP_dat(blocks, chara, M, num_vectors, random1, ldi, random2, ldo, a, b)



      implicit none

      integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start, num_vectors
      integer i, j, ii, jj, level, level_butterfly, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
      integer header_m, header_n, tailer_m, tailer_n, vector_inuse, mm, nn, num_blocks, level_define, col_vector
      integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
      integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2
      DT ctemp, a, b
      character chara
      type(matrixblock_MD)::blocks
      integer M, M1, N1, zfpflag,qttflag
      integer ldi,ldo
      real(kind=8)::tol_used
      DT :: random1(ldi, *), random2(ldo, *)
      DT:: al, be
      DT, allocatable :: random2tmp(:, :),random1tmp(:, :),random2tmp_1D(:)

      if(allocated(blocks%FullmatQTT%core))qttflag=1
      zfpflag=0
#if HAVE_ZFP
      if(allocated(blocks%FullmatZFP%buffer_r))zfpflag=1
      if(zfpflag==1)call ZFP_Decompress(blocks%fullmat,blocks%FullmatZFP,product(blocks%M),product(blocks%N),tol_used,1)
#endif
      M1=size(blocks%fullmat, 1)
      N1=size(blocks%fullmat, 2)


      al = 1d0
      be = 0d0

      call assert(M1 == M, 'M not equal fullmat dim')

      if (chara == 'N') then
         allocate (random2tmp(M1, num_vectors))
         random2tmp = random2(1:M1, 1:num_vectors)
         if(qttflag==1)then
            allocate(random1tmp(N1, num_vectors))
            random1tmp = random1(1:N1, 1:num_vectors)
            allocate(random2tmp_1D(M1*num_vectors))
            random2tmp_1D = 0
            call QTT_Apply_Fullvec(blocks%FullmatQTT,reshape(random1tmp,[N1*num_vectors]),random2tmp_1D)
            random2tmp = reshape(random2tmp_1D,[M1,num_vectors])
            deallocate(random1tmp)
            deallocate(random2tmp_1D)
         else
            call gemmf90(blocks%fullmat, M, random1, ldi, random2tmp, M1, 'N', 'N', M1, num_vectors, N1, BPACK_cone, BPACK_czero)
         endif
         random2(1:M1, 1:num_vectors) = a*random2tmp + b*random2(1:M1, 1:num_vectors)

      elseif (chara == 'T') then
         allocate (random2tmp(N1, num_vectors))
         random2tmp = random2(1:N1, 1:num_vectors)
         if(qttflag==1)then
            write(*,*)"QTT_Apply hasn't been implmented for transposed multiply"
            stop
         else
            call gemmf90(blocks%fullmat, M, random1, ldi, random2tmp, N1, 'T', 'N', N1, num_vectors, M1, al, be)
         endif
         random2(1:N1, 1:num_vectors) = a*random2tmp + b*random2(1:N1, 1:num_vectors)
      end if
      if(zfpflag==1)then
#if HAVE_ZFP
      call ZFP_Compress(blocks%fullmat,blocks%FullmatZFP,product(blocks%M),product(blocks%N),tol_used,1)
#endif
      endif

      ! write(*,*)'wo cao ni ma'
      deallocate (random2tmp)
   end subroutine Full_block_MD_MVP_dat



! compute arrays M_p(1:P+1,1:2,1:Ndim) and N_p(1:P+1,1:2,1:Ndim) the holds the start and end column/row of each process sharing this block
   subroutine ComputeParallelIndices_MD(block, pgno, Ndim, ptree, msh)
      implicit none
      type(matrixblock_MD)::block
      integer found, Ndim, pp, MyID, pgno, pgno1, level, level_p(Ndim), ith(Ndim), nleaf(Ndim), level_butterfly, nproc, num_blocks, gg, ii, ii_new, Maxlevel, dim_i, dim_s
      type(proctree)::ptree
      integer, pointer::M_p(:, :, :), N_p(:, :, :)
      type(mesh)::msh(Ndim)
      dim_s=1
      level_butterfly = block%level_butterfly

      if(.not. allocated(block%M_loc))allocate(block%M_loc(Ndim))
      block%M_loc = 0
      if(.not. allocated(block%N_loc))allocate(block%N_loc(Ndim))
      block%N_loc = 0

      Maxlevel = GetTreelevel(msh(1)%Maxgroup) - 1
      ! write(*,*)msh%Maxgroup,GetTreelevel(msh%Maxgroup),Maxlevel-block%level,block%level,ptree%nlevel-GetTreelevel(pgno),pgno,Maxlevel-block%level>=ptree%nlevel-GetTreelevel(pgno),pgno
      call assert((Maxlevel - block%level)*Ndim >= ptree%nlevel - GetTreelevel(pgno), 'too many process sharing this group')

      nproc = ptree%pgrp(pgno)%nproc

      if (associated(block%M_p)) deallocate (block%M_p)
      if (associated(block%N_p)) deallocate (block%N_p)
      allocate (block%M_p(nproc, 2, Ndim))
      allocate (block%N_p(nproc, 2, Ndim))
      M_p => block%M_p
      N_p => block%N_p

      do dim_i=1,Ndim
         M_p(:, 1, dim_i) = block%M(dim_i) + 1
         N_p(:, 1, dim_i) = block%N(dim_i) + 1
         M_p(:, 2, dim_i) = -block%M(dim_i) - 1
         N_p(:, 2, dim_i) = -block%N(dim_i) - 1
      enddo

      do pp=1,nproc
         MyID=ptree%pgrp(pgno)%head+pp-1
         found = 0
         level_p = 0
         pgno1 = pgno
         ith = 1
         dim_i = dim_s

         do while (found == 0)
            if (ptree%pgrp(pgno1)%nproc == 1) then
               found = 1
               exit
            endif
            if (MyID >= ptree%pgrp(pgno1*2)%head .and. MyID <= ptree%pgrp(pgno1*2)%tail) then
               pgno1 = pgno1*2
               ith(dim_i) = ith(dim_i)*2
            elseif (MyID >= ptree%pgrp(pgno1*2+1)%head .and. MyID <= ptree%pgrp(pgno1*2+1)%tail) then
               pgno1 = pgno1*2 + 1
               ith(dim_i) = ith(dim_i)*2 + 1
            endif
            level_p(dim_i) = level_p(dim_i) + 1
            dim_i = dim_i + 1
            dim_i = mod(dim_i-1,ndim)+1 ! reset dim to 1 if dim=ndim+1
         enddo
         ith = ith - 2**level_p + 1
         nleaf = 2**(level_butterfly - level_p)

         do dim_i=1,ndim
            gg = block%row_group(dim_i)*2**level_butterfly + (ith(dim_i) - 1)*nleaf(dim_i)
            M_p(pp, 1, dim_i) = min(M_p(pp, 1, dim_i), msh(dim_i)%basis_group(gg)%head - msh(dim_i)%basis_group(block%row_group(dim_i))%head + 1)

            gg = block%row_group(dim_i)*2**level_butterfly + ith(dim_i)*nleaf(dim_i) -1
            M_p(pp, 2, dim_i) = max(M_p(pp, 2, dim_i), msh(dim_i)%basis_group(gg)%tail - msh(dim_i)%basis_group(block%row_group(dim_i))%head + 1)

            gg = block%col_group(dim_i)*2**level_butterfly + (ith(dim_i) - 1)*nleaf(dim_i)
            N_p(pp, 1, dim_i) = min(N_p(pp, 1, dim_i), msh(dim_i)%basis_group(gg)%head - msh(dim_i)%basis_group(block%col_group(dim_i))%head + 1)

            gg = block%col_group(dim_i)*2**level_butterfly + ith(dim_i)*nleaf(dim_i) -1
            N_p(pp, 2, dim_i) = max(N_p(pp, 2, dim_i), msh(dim_i)%basis_group(gg)%tail - msh(dim_i)%basis_group(block%col_group(dim_i))%head + 1)
         enddo
      enddo

      if (IOwnPgrp(ptree, pgno)) then
         ii = ptree%myid - ptree%pgrp(pgno)%head + 1
         block%M_loc = M_p(ii, 2, :) - M_p(ii, 1, :) + 1
         block%N_loc = N_p(ii, 2, :) - N_p(ii, 1, :) + 1
      endif
      ! write(*,*)level_butterfly,level_p,block%M_loc,block%N_loc,'nima',M_p,N_p,block%M,block%N,block%row_group,block%col_group,nleaf
      ! endif
   end subroutine ComputeParallelIndices_MD


! compute arrays M_p(1:P+1) and N_p(1:P+1) the holds the start and end column/row of each process sharing this block
   subroutine ComputeParallelIndices(block, pgno, ptree, msh)
      implicit none
      type(matrixblock)::block
      integer pgno, level, level_p, level_butterfly, nproc, num_blocks, proc, gg, ii, ii_new, Maxlevel
      type(proctree)::ptree
      integer, pointer::M_p(:, :), N_p(:, :)
      type(mesh)::msh

      block%M_loc = 0
      block%N_loc = 0

      Maxlevel = GetTreelevel(msh%Maxgroup) - 1
      ! write(*,*)msh%Maxgroup,GetTreelevel(msh%Maxgroup),Maxlevel-block%level,block%level,ptree%nlevel-GetTreelevel(pgno),pgno,Maxlevel-block%level>=ptree%nlevel-GetTreelevel(pgno),pgno
      call assert(Maxlevel - block%level >= ptree%nlevel - GetTreelevel(pgno), 'too many process sharing this group')

      ! if(IOwnPgrp(ptree,pgno))then

      ! level_butterfly = block%level_butterfly
      level_p = ptree%nlevel - GetTreelevel(pgno)
      nproc = ptree%pgrp(pgno)%nproc
      num_blocks = 2**level_p

      if (associated(block%M_p)) deallocate (block%M_p)
      if (associated(block%N_p)) deallocate (block%N_p)
      allocate (block%M_p(nproc, 2))
      allocate (block%N_p(nproc, 2))
      M_p => block%M_p
      N_p => block%N_p

      M_p(:, 1) = block%M + 1
      N_p(:, 1) = block%N + 1
      M_p(:, 2) = -block%M - 1
      N_p(:, 2) = -block%N - 1

      do ii = 1, num_blocks

         ! if(flag==1)then  ! compute optimal renumbering of data pieces among the twice many processes
         ! if(mod(ii,2)==1)then
         ! ii_new=ceiling_safe(ii/2d0)
         ! else
         ! ii_new=ii/2+num_blocks/2
         ! endif
         ! else
         ii_new = ii
         ! endif

         gg = block%row_group*2**level_p + ii_new - 1
         proc = ptree%pgrp(pgno*2**level_p + ii - 1)%head - ptree%pgrp(pgno)%head
         M_p(proc + 1, 1) = min(M_p(proc + 1, 1), msh%basis_group(gg)%head - msh%basis_group(block%row_group)%head + 1)
         M_p(proc + 1, 2) = max(M_p(proc + 1, 2), msh%basis_group(gg)%tail - msh%basis_group(block%row_group)%head + 1)
         gg = block%col_group*2**level_p + ii_new - 1
         N_p(proc + 1, 1) = min(N_p(proc + 1, 1), msh%basis_group(gg)%head - msh%basis_group(block%col_group)%head + 1)
         N_p(proc + 1, 2) = max(N_p(proc + 1, 2), msh%basis_group(gg)%tail - msh%basis_group(block%col_group)%head + 1)
      enddo

      if (IOwnPgrp(ptree, pgno)) then
         ii = ptree%myid - ptree%pgrp(pgno)%head + 1
         block%M_loc = M_p(ii, 2) - M_p(ii, 1) + 1
         block%N_loc = N_p(ii, 2) - N_p(ii, 1) + 1
      endif
      ! write(*,*)level_butterfly,level_p,block%M_loc,block%N_loc,'nima',M_p,N_p,block%M,block%N,block%row_group,block%col_group
      ! endif
   end subroutine ComputeParallelIndices

! compute arrays M_p(1:P+1) or N_p(1:P+1) the holds the start and end column/row of each process sharing this block
   subroutine ComputeParallelIndicesSub(base_group, pgno, ptree, msh, MN_p)
      implicit none

      integer pgno, level, level_p, level_butterfly, nproc, num_blocks, proc, gg, ii, ii_new, Maxlevel, Baselevel, base_group
      type(proctree)::ptree
      integer::MN_p(:, :)
      type(mesh)::msh

      Maxlevel = GetTreelevel(msh%Maxgroup) - 1
      Baselevel = GetTreelevel(base_group) - 1

      call assert(Maxlevel - Baselevel >= ptree%nlevel - GetTreelevel(pgno), 'too many process sharing this group')

      level_p = ptree%nlevel - GetTreelevel(pgno)
      nproc = ptree%pgrp(pgno)%nproc
      num_blocks = 2**level_p

      MN_p(:, 1) = BPACK_BigINT
      MN_p(:, 2) = -BPACK_BigINT

      do ii = 1, num_blocks
         ii_new = ii
         gg = base_group*2**level_p + ii_new - 1
         proc = ptree%pgrp(pgno*2**level_p + ii - 1)%head - ptree%pgrp(pgno)%head
         MN_p(proc + 1, 1) = min(MN_p(proc + 1, 1), msh%basis_group(gg)%head - msh%basis_group(base_group)%head + 1)
         MN_p(proc + 1, 2) = max(MN_p(proc + 1, 2), msh%basis_group(gg)%tail - msh%basis_group(base_group)%head + 1)
      enddo

   end subroutine ComputeParallelIndicesSub

   function node_score_block_ptr_row(this) result(score)
      implicit none
      type(nod)::this
      real(kind=8)::score
      class(*), pointer::ptr

      select TYPE (ptr=>this%item)
      type is (block_ptr)
         score = dble(ptr%ptr%row_group)
      class default
         write (*, *) 'unexpected item type in node_score_dble'
         stop
      end select
   end function node_score_block_ptr_row

   function nod_score_ipair(this) result(score)
      implicit none
      type(nod)::this
      real(kind=8)::score
      class(*), pointer::ptr

      select TYPE (ptr=>this%item)
      type is (ipair)
         score = dble(ptr%i)
      class default
         write (*, *) 'unexpected item type in nod_score_ipair'
         stop
      end select
   end function nod_score_ipair


   subroutine element_Zmn_block_user(nrow, ncol, mrange, nrange, values, msh, option, ker, myflag, passflag, ptree, stats)


      implicit none

      integer ii, jj, nn, pp, ij, i, j, nrow, ncol, passflag, myflag, Ninter, idx, nc, nr, pgno, ctxt, nprow, npcol, myrow, mycol
      integer mrange(nrow)
      integer nrange(ncol)
      DT:: value_e, values(nrow, ncol)
      type(mesh)::msh
      type(proctree)::ptree
      type(Hoption)::option
      type(Hstat)::stats
      type(kernelquant)::ker
      integer ierr, idx_row, idx_col
      integer*8 idx_dat
      integer, allocatable:: flags(:), dests(:), colidx1(:), rowidx1(:), colidx(:), rowidx(:), allrows(:), allcols(:), disps(:), pgidx(:), pmaps(:, :)
      procedure(F_Zelem_block), POINTER :: proc
      procedure(C_Zelem_block), POINTER :: proc_c
      procedure(F_Zelem), POINTER :: proc1
      procedure(C_Zelem), POINTER :: proc1_c
      DT, allocatable::alldat_loc(:)
      type(intersect), allocatable::inters(:)
      integer myArows, myAcols, Npmap
      real(kind=8)::t1, t2, t3, t4
      integer reqm, reqn
      integer statusm(MPI_status_size), statusn(MPI_status_size)

      if (option%elem_extract == 0) then

         t1 = MPI_Wtime()

         if (option%cpp == 1) then
            call c_f_procpointer(ker%C_FuncZmn, proc1_C)
#ifdef HAVE_TASKLOOP
            !$omp parallel
            !$omp single
            !$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
#endif
            do ij = 1, ncol*nrow
               jj = (ij - 1)/nrow + 1
               ii = mod(ij - 1, nrow) + 1
               value_e = 0
               call proc1_C(msh%new2old(mrange(ii)), msh%new2old(nrange(jj)), value_e, ker%C_QuantApp)
               value_e = value_e*option%scale_factor
               values(ii, jj) = value_e
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
            !$omp end single
            !$omp end parallel
#else
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
#endif
         else
            proc1 => ker%FuncZmn
            if (nrow*ncol > 0) then
#ifdef HAVE_TASKLOOP
               !$omp parallel
               !$omp single
               !$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
#endif
               do ij = 1, ncol*nrow
                  jj = (ij - 1)/nrow + 1
                  ii = mod(ij - 1, nrow) + 1
                  value_e = 0
                  call proc1(msh%new2old(mrange(ii)), msh%new2old(nrange(jj)), value_e, ker%QuantApp)
                  value_e = value_e*option%scale_factor
                  values(ii, jj) = value_e
               enddo
#ifdef HAVE_TASKLOOP
               !$omp end taskloop
               !$omp end single
               !$omp end parallel
#else
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
#endif

            endif
         endif

         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1

         passflag = 2
      else if (option%elem_extract >= 1) then

         allocate (flags(ptree%nproc))

#ifdef HAVE_MPI3
         call MPI_IALLGATHER(myflag, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
         call MPI_Wait(reqm, statusm, ierr)
#else
         call MPI_ALLGATHER(myflag, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif

         passflag = minval(flags)

         t1 = MPI_Wtime()

         if (passflag == 0) then

            allocate (colidx1(ptree%nproc))
            allocate (rowidx1(ptree%nproc))
            allocate (disps(ptree%nproc))

#ifdef HAVE_MPI3
            call MPI_IALLGATHER(nrow, 1, MPI_INTEGER, rowidx1, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
            call MPI_IALLGATHER(ncol, 1, MPI_INTEGER, colidx1, 1, MPI_INTEGER, ptree%Comm, reqn, ierr)
#else
            call MPI_ALLGATHER(nrow, 1, MPI_INTEGER, rowidx1, 1, MPI_INTEGER, ptree%Comm, ierr)
            call MPI_ALLGATHER(ncol, 1, MPI_INTEGER, colidx1, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif

            Npmap = ptree%nproc
            allocate (pmaps(Npmap, 3))
            do pp = 1, Npmap
               pmaps(pp, 1) = 1
               pmaps(pp, 2) = 1
               pmaps(pp, 3) = pp - 1
            enddo

            Ninter = 0
            do pp = 1, ptree%nproc
            if (flags(pp) == 0) then
               Ninter = Ninter + 1
            endif
            enddo

            allocate (colidx(Ninter))
            allocate (rowidx(Ninter))
            allocate (pgidx(Ninter))

#ifdef HAVE_MPI3
            call MPI_Wait(reqm, statusm, ierr)
            call MPI_Wait(reqn, statusn, ierr)
#endif
            !>***** Count number of active intersections Ninter
            Ninter = 0
            do pp = 1, ptree%nproc
            if (flags(pp) == 0) then
               Ninter = Ninter + 1
               pgidx(Ninter) = pp
               rowidx(Ninter) = rowidx1(pp)
               colidx(Ninter) = colidx1(pp)
            endif
            enddo

            !>***** count number of local data
            idx_dat = 0
            do nn = 1, Ninter
               nr = rowidx(nn)
               nc = colidx(nn)
               ! datidx(nn)=ntot_loc
               nprow = pmaps(pgidx(nn), 1)
               npcol = pmaps(pgidx(nn), 2)
               call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
               if (myrow /= -1 .and. mycol /= -1) then
                  myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
                  idx_dat = idx_dat + myArows*myAcols
               endif
            enddo
            allocate (alldat_loc(idx_dat))
            if (idx_dat > 0) alldat_loc = 0

            !>***** Broadcast mrange and nrange for each intersection
            idx_row = sum(rowidx)
            allocate (allrows(idx_row))
            idx = 0
            do pp = 1, ptree%nproc
               disps(pp) = idx
               idx = idx + rowidx1(pp)
            enddo

#ifdef HAVE_MPI3
            call MPI_IALLGATHERV(mrange, nrow, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, reqm, ierr)
#else
            call MPI_ALLGATHERV(mrange, nrow, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif

            idx_col = sum(colidx)
            allocate (allcols(idx_col))
            idx = 0
            do pp = 1, ptree%nproc
               disps(pp) = idx
               idx = idx + colidx1(pp)
            enddo
#ifdef HAVE_MPI3
            call MPI_IALLGATHERV(nrange, ncol, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, reqn, ierr)
            call MPI_Wait(reqm, statusm, ierr)
            call MPI_Wait(reqn, statusn, ierr)
#else
            call MPI_ALLGATHERV(nrange, ncol, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif
            if (option%cpp == 1) then
               call c_f_procpointer(ker%C_FuncZmnBlock, proc_C)
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               pgidx = pgidx - 1
               call proc_C(Ninter, idx_row, idx_col, idx_dat, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%C_QuantApp)
            else
               proc => ker%FuncZmnBlock
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               call proc(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%QuantApp)
            endif

            idx = 0
            do jj = 1, ncol ! note that alldat_loc has column major
            do ii = 1, nrow
               idx = idx + 1
               values(ii, jj) = alldat_loc(idx)*option%scale_factor
            enddo
            enddo

            deallocate (allrows)
            deallocate (allcols)
            deallocate (colidx1)
            deallocate (colidx)
            deallocate (rowidx1)
            deallocate (rowidx)
            deallocate (disps)
            deallocate (alldat_loc)
            deallocate (pgidx)
            deallocate (pmaps)

         endif
         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1
         deallocate (flags)
      endif

      return

   end subroutine element_Zmn_block_user





   subroutine element_Zmn_blocklist_user(submats, Nsub, msh, option, ker, myflag, passflag, ptree, stats,alldat_loc_in)


      implicit none


      integer ii, jj, nn, pp, ij, i, j, nrow, ncol, passflag, myflag,myflag1, Ninter, Ninter_loc,Nsub,Nzero, idx, nc, nr, pgno, ctxt, nprow, npcol, myrow, mycol
      type(mesh)::msh
      type(proctree)::ptree
      type(Hoption)::option
      type(Hstat)::stats
      type(kernelquant)::ker
      integer ierr, idx_row, idx_col
      integer*8 idx_dat
      integer, allocatable:: flags(:), dests(:), colidx(:), rowidx(:), colidx1(:), rowidx1(:), allrows(:), allcols(:), disps(:), pgidx(:), pmaps(:, :),nsubs(:),mrange(:),nrange(:),nrows_loc(:),ncols_loc(:)
      procedure(F_Zelem_block), POINTER :: proc
      procedure(C_Zelem_block), POINTER :: proc_c
      procedure(F_Zelem), POINTER :: proc1
      procedure(C_Zelem), POINTER :: proc1_c
      DT,pointer::alldat_loc(:)
      DT,target,optional::alldat_loc_in(:)
      type(intersect), allocatable::inters(:)
      integer myArows, myAcols, Npmap
      real(kind=8)::t1, t2, t3, t4
      integer reqm, reqn,empty
      integer statusm(MPI_status_size), statusn(MPI_status_size)
      type(intersect) :: submats(*)
      type(intersect),allocatable :: submats1(:)
      integer,allocatable:: new2old_sub(:)
      DT:: value_e

      allocate(submats1(Nsub))
      allocate(new2old_sub(Nsub))

      t1 = MPI_Wtime()

      !!! copy submats into submats1 as submats can have zero elements in the masks
      Ninter_loc=0
      new2old_sub=-1
      nr=0
      nc=0
      idx_dat=0
      do nn=1,Nsub

         nrow = 0
         do i = 1, submats(nn)%nr
            empty = 0
            if(allocated(submats(nn)%masks))then
               if (sum(submats(nn)%masks(i, :)) == 0) then
                  empty=1
               endif
            endif
            if (empty == 0) then
               nrow = nrow + 1
            endif
         enddo
         ncol = 0
         do j = 1, submats(nn)%nc
            empty = 0
            if(allocated(submats(nn)%masks))then
               if (sum(submats(nn)%masks(:, j)) == 0) then
                  empty=1
               endif
            endif
            if (empty == 0) then
               ncol = ncol + 1
            endif
         enddo
         if(ncol*nrow>0)then
            Ninter_loc = Ninter_loc + 1
            new2old_sub(Ninter_loc) = nn

            allocate(submats1(Ninter_loc)%rows(nrow))
            allocate(submats1(Ninter_loc)%mmap(nrow))
            nrow = 0
            do i = 1, submats(nn)%nr
               empty = 0
               if(allocated(submats(nn)%masks))then
                  if (sum(submats(nn)%masks(i, :)) == 0) then
                     empty=1
                  endif
               endif
               if (empty == 0) then
                  nrow = nrow + 1
                  nr =nr +1
                  submats1(Ninter_loc)%mmap(nrow) = i
                  submats1(Ninter_loc)%rows(nrow) = submats(nn)%rows(i)
               endif
            enddo
            submats1(Ninter_loc)%nr = nrow

            allocate(submats1(Ninter_loc)%cols(ncol))
            allocate(submats1(Ninter_loc)%nmap(ncol))
            ncol = 0
            do j = 1, submats(nn)%nc
               empty = 0
               if(allocated(submats(nn)%masks))then
                  if (sum(submats(nn)%masks(:, j)) == 0) then
                     empty=1
                  endif
               endif
               if (empty == 0) then
                  ncol = ncol + 1
                  nc =nc+1
                  submats1(Ninter_loc)%nmap(ncol) = j
                  submats1(Ninter_loc)%cols(ncol) = submats(nn)%cols(j)
               endif
            enddo
            submats1(Ninter_loc)%nc = ncol


            if(present(alldat_loc_in))then
               call Array1DtoPointer2D(alldat_loc_in(idx_dat+1:idx_dat + nrow*ncol), submats(nn)%dat, nrow, ncol)
               idx_dat = idx_dat + nrow*ncol
            endif

         endif
      enddo

      t2 = MPI_Wtime()
      stats%Time_Entry = stats%Time_Entry + t2 - t1

      myflag1 = myflag
      if(Ninter_loc==0)myflag1=max(myflag,1)

      if (option%elem_extract == 0) then

         t1 = MPI_Wtime()

         do nn=1,Ninter_loc
            if (option%cpp == 1) then
               call c_f_procpointer(ker%C_FuncZmn, proc1_C)
#ifdef HAVE_TASKLOOP
               !$omp parallel
               !$omp single
               !$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
#endif
               do ij = 1, submats1(nn)%nc*submats1(nn)%nr
                  jj = (ij - 1)/submats1(nn)%nr + 1
                  ii = mod(ij - 1, submats1(nn)%nr) + 1
                  value_e = 0
                  call proc1_C(msh%new2old(submats1(nn)%rows(ii)), msh%new2old(submats1(nn)%cols(jj)), value_e, ker%C_QuantApp)
                  value_e = value_e*option%scale_factor
                  submats(new2old_sub(nn))%dat(submats1(nn)%mmap(ii), submats1(nn)%nmap(jj)) = value_e
               enddo
#ifdef HAVE_TASKLOOP
               !$omp end taskloop
               !$omp end single
               !$omp end parallel
#else
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
#endif
            else
               proc1 => ker%FuncZmn
               if (submats1(nn)%nr*submats1(nn)%nc > 0) then
#ifdef HAVE_TASKLOOP
                  !$omp parallel
                  !$omp single
                  !$omp taskloop default(shared) private(ij,ii,jj,value_e)
#else
#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(ij,ii,jj,value_e)
#endif
#endif
                  do ij = 1, submats1(nn)%nc*submats1(nn)%nr
                     jj = (ij - 1)/submats1(nn)%nr + 1
                     ii = mod(ij - 1, submats1(nn)%nr) + 1
                     value_e = 0
                     call proc1(msh%new2old(submats1(nn)%rows(ii)), msh%new2old(submats1(nn)%cols(jj)), value_e, ker%QuantApp)
                     value_e = value_e*option%scale_factor
                     submats(new2old_sub(nn))%dat(submats1(nn)%mmap(ii), submats1(nn)%nmap(jj)) = value_e
                  enddo
#ifdef HAVE_TASKLOOP
                  !$omp end taskloop
                  !$omp end single
                  !$omp end parallel
#else
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif
#endif

               endif
            endif
         enddo
         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1

         passflag = 2
      else if (option%elem_extract == 1) then


         ! write(*,*)ptree%MyID,'in',Nsub,nr,nc
         allocate (flags(ptree%nproc))

#ifdef HAVE_MPI3
         call MPI_IALLGATHER(myflag1, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
         call MPI_Wait(reqm, statusm, ierr)
#else
         call MPI_ALLGATHER(myflag1, 1, MPI_INTEGER, flags, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif

         passflag = minval(flags)

         t1 = MPI_Wtime()

         if (passflag == 0) then

            allocate (disps(ptree%nproc))
            allocate (nsubs(ptree%nproc))
#ifdef HAVE_MPI3
            call MPI_IALLGATHER(Ninter_loc, 1, MPI_INTEGER, nsubs, 1, MPI_INTEGER, ptree%Comm, reqm, ierr)
            call MPI_Wait(reqm, statusm, ierr)
#else
            call MPI_ALLGATHER(Ninter_loc, 1, MPI_INTEGER, nsubs, 1, MPI_INTEGER, ptree%Comm, ierr)
#endif


            allocate(nrows_loc(max(1,Ninter_loc)))
            allocate(ncols_loc(max(1,Ninter_loc)))
            allocate(mrange(max(1,nr)))
            allocate(nrange(max(1,nc)))
            nr=0
            nc=0
            do nn=1,Ninter_loc
               do i = 1, submats1(nn)%nr
                  nr =nr +1
                  mrange(nr)=submats1(nn)%rows(i)
               enddo
               nrows_loc(nn) = submats1(nn)%nr
               do i = 1, submats1(nn)%nc
                  nc =nc +1
                  nrange(nc)=submats1(nn)%cols(i)
               enddo
               ncols_loc(nn) = submats1(nn)%nc
            enddo

            idx = 0
            do pp = 1, ptree%nproc
               disps(pp) = idx
               idx = idx + nsubs(pp)
            enddo

            Ninter = sum(nsubs)

            allocate (colidx(Ninter))
            allocate (rowidx(Ninter))


#ifdef HAVE_MPI3
            call MPI_IALLGATHERV(nrows_loc, Ninter_loc, MPI_INTEGER, rowidx, nsubs, disps, MPI_INTEGER, ptree%comm, reqm, ierr)
            call MPI_IALLGATHERV(ncols_loc, Ninter_loc, MPI_INTEGER, colidx, nsubs, disps, MPI_INTEGER, ptree%comm, reqn, ierr)
            call MPI_Wait(reqm, statusm, ierr)
            call MPI_Wait(reqn, statusn, ierr)
#else
            call MPI_ALLGATHERV(nrows_loc, Ninter_loc, MPI_INTEGER, rowidx, nsubs, disps, MPI_INTEGER, ptree%comm,ierr)
            call MPI_ALLGATHERV(ncols_loc, Ninter_loc, MPI_INTEGER, colidx, nsubs, disps, MPI_INTEGER, ptree%comm,ierr)
#endif

            !>***** Generate pmaps and pgidx for all intersections
            Npmap = ptree%nproc
            allocate (pmaps(Npmap, 3))
            do pp = 1, Npmap
               pmaps(pp, 1) = 1
               pmaps(pp, 2) = 1
               pmaps(pp, 3) = pp - 1
            enddo

            allocate (pgidx(Ninter))
            Ninter = 0
            do pp = 1, ptree%nproc
               do nn = 1, nsubs(pp)
                  Ninter = Ninter + 1
                  pgidx(Ninter) = pp
               enddo
            enddo

            !>***** count number of local data
            idx_dat = 0
            do nn = 1, Ninter
               nrow = rowidx(nn)
               ncol = colidx(nn)
               ! datidx(nn)=ntot_loc
               nprow = pmaps(pgidx(nn), 1)
               npcol = pmaps(pgidx(nn), 2)
               call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
               if (myrow /= -1 .and. mycol /= -1) then
                  myArows = numroc_wp(nrow, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(ncol, nbslpk, mycol, 0, npcol)
                  idx_dat = idx_dat + myArows*myAcols
               endif
            enddo
            if(present(alldat_loc_in))then
               alldat_loc => alldat_loc_in
            else
               allocate (alldat_loc(idx_dat))
               call LogMemory(stats, SIZEOF(alldat_loc)/1024.0d3)
            endif
            if (idx_dat > 0) alldat_loc(1:idx_dat) = 0

            allocate (rowidx1(ptree%nproc))
            allocate (colidx1(ptree%nproc))

            Ninter = 0
            do pp = 1, ptree%nproc
               rowidx1(pp) = 0
               colidx1(pp) = 0
               do nn = 1, nsubs(pp)
                  Ninter = Ninter + 1
                  rowidx1(pp) = rowidx1(pp) + rowidx(Ninter)
                  colidx1(pp) = colidx1(pp) + colidx(Ninter)
               enddo
            enddo


            !>***** Broadcast mrange and nrange for each intersection
            idx_row = sum(rowidx1)
            allocate (allrows(idx_row))
            idx = 0
            do pp = 1, ptree%nproc
               disps(pp) = idx
               idx = idx + rowidx1(pp)
            enddo

#ifdef HAVE_MPI3
            call MPI_IALLGATHERV(mrange, nr, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, reqm, ierr)
#else
            call MPI_ALLGATHERV(mrange, nr, MPI_INTEGER, allrows, rowidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif

            idx_col = sum(colidx1)
            allocate (allcols(idx_col))
            idx = 0
            do pp = 1, ptree%nproc
               disps(pp) = idx
               idx = idx + colidx1(pp)
            enddo
#ifdef HAVE_MPI3
            call MPI_IALLGATHERV(nrange, nc, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, reqn, ierr)
            call MPI_Wait(reqm, statusm, ierr)
            call MPI_Wait(reqn, statusn, ierr)
#else
            call MPI_ALLGATHERV(nrange, nc, MPI_INTEGER, allcols, colidx1, disps, MPI_INTEGER, ptree%comm, ierr)
#endif
            if (option%cpp == 1) then
               call c_f_procpointer(ker%C_FuncZmnBlock, proc_C)
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               pgidx = pgidx - 1
               call proc_C(Ninter, idx_row, idx_col, idx_dat, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%C_QuantApp)
            else
               proc => ker%FuncZmnBlock
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               call proc(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%QuantApp)
            endif

            if(.not. present(alldat_loc_in))then
               idx = 0
               do nn=1,Ninter_loc
               do jj = 1, submats1(nn)%nc ! note that alldat_loc has column major
               do ii = 1, submats1(nn)%nr
                  idx = idx + 1
                  submats(new2old_sub(nn))%dat(submats1(nn)%mmap(ii), submats1(nn)%nmap(jj)) = alldat_loc(idx)*option%scale_factor
               enddo
               enddo
               enddo
               deallocate (alldat_loc)
            endif

            deallocate (allrows)
            deallocate (allcols)
            deallocate (colidx)
            deallocate (rowidx)
            deallocate (colidx1)
            deallocate (rowidx1)
            deallocate (disps)
            deallocate (pgidx)
            deallocate (pmaps)
            deallocate (nsubs)
            deallocate (nrows_loc)
            deallocate (ncols_loc)
            deallocate (mrange)
            deallocate (nrange)

         endif
         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1
         deallocate (flags)
         ! write(*,*)ptree%MyID,'out'


      else if (option%elem_extract == 2) then
         if(Ninter_loc>0)then
            allocate(rowidx(Ninter_loc))  ! rowidx
            allocate(colidx(Ninter_loc))  ! colidx
            allocate(allrows(nr))             ! allrows
            allocate(allcols(nc))             ! allcols
            idx_row=nr
            idx_col=nc

            nr=0
            nc=0
            do nn=1,Ninter_loc
               do i = 1, submats1(nn)%nr
                  nr =nr +1
                  allrows(nr)=submats1(nn)%rows(i)
               enddo
               rowidx(nn) = submats1(nn)%nr
               do i = 1, submats1(nn)%nc
                  nc =nc +1
                  allcols(nc)=submats1(nn)%cols(i)
               enddo
               colidx(nn) = submats1(nn)%nc
            enddo

            !>***** Generate pmaps and pgidx for all intersections
            Npmap = 1
            Ninter = Ninter_loc
            allocate (pmaps(Npmap, 3))
            do pp = 1, Npmap
               pmaps(pp, 1) = 1
               pmaps(pp, 2) = 1
               pmaps(pp, 3) = ptree%MyID
            enddo

            allocate (pgidx(Ninter))
            pgidx=1

            !>***** count number of local data
            idx_dat = 0
            do nn = 1, Ninter
               nrow = rowidx(nn)
               ncol = colidx(nn)
               ! datidx(nn)=ntot_loc
               nprow = pmaps(pgidx(nn), 1)
               npcol = pmaps(pgidx(nn), 2)
               call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
               if (myrow /= -1 .and. mycol /= -1) then
                  myArows = numroc_wp(nrow, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(ncol, nbslpk, mycol, 0, npcol)
                  idx_dat = idx_dat + myArows*myAcols
               endif
            enddo
            if(present(alldat_loc_in))then
               alldat_loc => alldat_loc_in
            else
               allocate (alldat_loc(idx_dat))
               call LogMemory(stats, SIZEOF(alldat_loc)/1024.0d3)
            endif
            if (idx_dat > 0) alldat_loc(1:idx_dat) = 0

            if (option%cpp == 1) then
               call c_f_procpointer(ker%C_FuncZmnBlock, proc_C)
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               pgidx = pgidx - 1
               call proc_C(Ninter, idx_row, idx_col, idx_dat, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%C_QuantApp)
            else
               proc => ker%FuncZmnBlock
               ! !>***** parallel extraction of the data
               do ii = 1, idx_row
                  allrows(ii) = abs(msh%new2old(allrows(ii)))
               enddo
               do jj = 1, idx_col
                  allcols(jj) = abs(msh%new2old(allcols(jj)))
               enddo
               call proc(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, ker%QuantApp)
            endif

            if(.not. present(alldat_loc_in))then
               idx = 0
               do nn=1,Ninter_loc
               do jj = 1, submats1(nn)%nc ! note that alldat_loc has column major
               do ii = 1, submats1(nn)%nr
                  idx = idx + 1
                  submats(new2old_sub(nn))%dat(submats1(nn)%mmap(ii), submats1(nn)%nmap(jj)) = alldat_loc(idx)*option%scale_factor
               enddo
               enddo
               enddo
               deallocate (alldat_loc)
            endif

            deallocate (allrows)
            deallocate (allcols)
            deallocate (colidx)
            deallocate (rowidx)
            deallocate (pgidx)
            deallocate (pmaps)
         endif
         passflag=2
         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1
      endif


      do nn = 1, Ninter_loc
         if(allocated(submats1(nn)%mmap))deallocate(submats1(nn)%mmap)
         if(allocated(submats1(nn)%nmap))deallocate(submats1(nn)%nmap)
         if(allocated(submats1(nn)%rows))deallocate(submats1(nn)%rows)
         if(allocated(submats1(nn)%cols))deallocate(submats1(nn)%cols)
      enddo
      deallocate(submats1)
      deallocate(new2old_sub)

      return

   end subroutine element_Zmn_blocklist_user







   subroutine element_Zmn_tensorlist_user(Ndim, subtensors, Nsub, msh, option, ker, myflag, passflag, ptree, stats, zfpquants, qttquants)


      implicit none


      integer Ndim, ii, jj, ij, nn, pp, i, j, nrow, ncol, passflag, myflag,myflag1, Ninter, Ninter_loc,Nsub,Nzero, idx, nc, nr, pgno, ctxt, nprow, npcol, myrow, mycol, dim_i
      type(mesh)::msh(Ndim)
      integer:: dims(2*Ndim),num_threads
      integer*8:: dims8(2*Ndim)
      integer,allocatable::idxs(:,:)
      type(proctree)::ptree
      type(Hoption)::option
      type(Hstat)::stats
      type(kernelquant)::ker
      integer ierr, idx_row, idx_col
      integer*8 idx_dat
      integer, allocatable:: flags(:), dests(:), colidx(:), rowidx(:), colidx1(:), rowidx1(:), allrows(:), allcols(:), disps(:), pgidx(:), pmaps(:, :),nsubs(:),mrange(:),nrange(:),nrows_loc(:),ncols_loc(:)
      procedure(F_Zelem_block), POINTER :: proc
      procedure(C_Zelem_block), POINTER :: proc_c
      procedure(F_Zelem_MD), POINTER :: proc1
      procedure(C_Zelem_MD), POINTER :: proc1_c
      DT,pointer::alldat_loc(:)
      type(intersect), allocatable::inters(:)
      integer myArows, myAcols, Npmap
      real(kind=8)::t1, t2, t3, t4, tmpmem
      integer reqm, reqn,empty
      integer statusm(MPI_status_size), statusn(MPI_status_size)
      type(intersect_MD) :: subtensors(*)
      type(zfpquant),optional:: zfpquants(*)
      type(TTtype),optional:: qttquants(*)
      DT:: value_e
integer, save:: my_tid = 0
#ifdef HAVE_OPENMP
!$omp threadprivate(my_tid)
#endif

      myflag1 = myflag
      if(Nsub==0)myflag1=max(myflag,1)

#ifdef HAVE_OPENMP
!$omp parallel default(shared)
!$omp master
      num_threads = omp_get_num_threads()
!$omp end master
      my_tid = omp_get_thread_num()
!$omp end parallel
#else
      num_threads = 1
      my_tid = 0
#endif

      if (option%elem_extract == 0) then

         t1 = MPI_Wtime()
         if (option%cpp == 1) then
            call c_f_procpointer(ker%C_FuncZmn_MD, proc1_C)
         else
            proc1 => ker%FuncZmn_MD
         endif
         do nn=1,Nsub
            dims(1:Ndim) = subtensors(nn)%nr
            dims(1+Ndim:2*Ndim) = subtensors(nn)%nc
            dims8 = dims
            if (product(dims8)> 0) then
               if(present(zfpquants) .or. present(qttquants))then
                  allocate(subtensors(nn)%dat(product(subtensors(nn)%nr),product(subtensors(nn)%nc)))

#if 0
! not sure why the following is causing compiling error for certain intel compilers
                  call LogMemory(stats, SIZEOF(subtensors(nn)%dat)/1024.0d3)
#else
#if DAT==0
                  call LogMemory(stats, 16*SIZE(subtensors(nn)%dat)/1024.0d3)
#elif DAT==1
                  call LogMemory(stats, 8*SIZE(subtensors(nn)%dat)/1024.0d3)
#elif DAT==2
                  call LogMemory(stats, 8*SIZE(subtensors(nn)%dat)/1024.0d3)
#elif DAT==3
                  call LogMemory(stats, 4*SIZE(subtensors(nn)%dat)/1024.0d3)
#endif
#endif
               endif
               allocate(idxs(2*Ndim,num_threads))
#ifdef HAVE_TASKLOOP
               !$omp parallel
               !$omp single
               !$omp taskloop default(shared) private(ij,i,j,value_e,dim_i,empty)
#else
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(ij,i,j,value_e,dim_i,empty)
#endif
#endif
               do ij = 1, product(dims8)
                  call SingleIndexToMultiIndex(2*Ndim,dims, ij, idxs(:,my_tid+1))

                  call MultiIndexToSingleIndex(Ndim,dims(1:Ndim), i, idxs(1:Ndim,my_tid+1))
                  call MultiIndexToSingleIndex(Ndim,dims(1+Ndim:2*Ndim), j, idxs(1+Ndim:2*Ndim,my_tid+1))

                  empty = 0
                  if(allocated(subtensors(nn)%masks))then
                     if (subtensors(nn)%masks(i,j) == 0) then
                        empty=1
                     endif
                  endif
                  if(empty==0)then
                     do dim_i=1,Ndim
                        idxs(dim_i,my_tid+1) = msh(dim_i)%new2old(subtensors(nn)%rows(dim_i)%dat(idxs(dim_i,my_tid+1)))
                     enddo
                     do dim_i=1,Ndim
                        idxs(dim_i+Ndim,my_tid+1) = msh(dim_i)%new2old(subtensors(nn)%cols(dim_i)%dat(idxs(dim_i+Ndim,my_tid+1)))
                     enddo

                     value_e = 0
                     if (option%cpp == 1) then
                        call proc1_C(Ndim, idxs(1:Ndim,my_tid+1), idxs(1+Ndim:2*Ndim,my_tid+1),value_e, ker%C_QuantApp)
                     else
                        call proc1(Ndim, idxs(1:Ndim,my_tid+1), idxs(1+Ndim:2*Ndim,my_tid+1),value_e, ker%QuantApp)
                     endif
                     value_e = value_e*option%scale_factor
                     subtensors(nn)%dat(i, j) = value_e
                  else
                     subtensors(nn)%dat(i, j) = 0
                  endif
               enddo
#ifdef HAVE_TASKLOOP
               !$omp end taskloop
               !$omp end single
               !$omp end parallel
#else
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
#endif
#if HAVE_ZFP
               if(present(zfpquants))then
                  tmpmem = SIZEOF(subtensors(nn)%dat)/1024.0d3
                  call ZFP_Compress(subtensors(nn)%dat,zfpquants(nn),product(subtensors(nn)%nr),product(subtensors(nn)%nc),option%tol_comp,0)
                  if(allocated(zfpquants(nn)%buffer_r))call LogMemory(stats, SIZEOF(zfpquants(nn)%buffer_r)/1024.0d3)
                  if(allocated(zfpquants(nn)%buffer_i))call LogMemory(stats, SIZEOF(zfpquants(nn)%buffer_i)/1024.0d3)
                  call LogMemory(stats, -tmpmem)
               endif
#endif
               if(present(qttquants))then
                  tmpmem = SIZEOF(subtensors(nn)%dat)/1024.0d3

                  qttquants(nn)%d_org = Ndim
                  qttquants(nn)%mpo = 1
                  allocate(qttquants(nn)%m_n_org(qttquants(nn)%d_org,2))
                  qttquants(nn)%m_n_org(:,1)=dims(1:Ndim)
                  qttquants(nn)%m_n_org(:,2)=dims(1+Ndim:2*Ndim)
                  call QTT_Compress_SVD(reshape(subtensors(nn)%dat,[product(dims)]),option%tol_comp,qttquants(nn),option%use_zfp)
                  deallocate(subtensors(nn)%dat)
                  if(allocated(qttquants(nn)%core))call LogMemory(stats, SIZEOF(qttquants(nn)%core)/1024.0d3)

                  if(allocated(qttquants(nn)%coreZFP%buffer_r))call LogMemory(stats, SIZEOF(qttquants(nn)%coreZFP%buffer_r)/1024.0d3)
                  if(allocated(qttquants(nn)%coreZFP%buffer_i))call LogMemory(stats, SIZEOF(qttquants(nn)%coreZFP%buffer_i)/1024.0d3)

                  call LogMemory(stats, -tmpmem)
               endif
               deallocate(idxs)
            endif
         enddo
         t2 = MPI_Wtime()
         stats%Time_Entry = stats%Time_Entry + t2 - t1

         passflag = 2
      else if (option%elem_extract == 1) then
         write(*,*)"elem_extract == 1 not implemented in element_Zmn_tensorlist_user"
      else if (option%elem_extract == 2) then
         write(*,*)"elem_extract == 2 not implemented in element_Zmn_tensorlist_user"
      endif

      return

   end subroutine element_Zmn_tensorlist_user







end module Bplus_Utilities
